import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable
import numba as nb
from mendeleev import element
from functools import lru_cache
import multiprocessing as mp

try:
    import cupy as cp
    from numba import cuda
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False

@dataclass
class ElectrodeConfig:
    start: float
    width: float
    ap_start: float
    ap_width: float
    outer_diameter: float
    voltage: float

@nb.njit(parallel=True)
def solve_field_cpu(potential, mask, maxiter, thresh, omega=1.9):
    for _ in range(maxiter):
        vold = potential.copy()
        
        for i in nb.prange(1, potential.shape[0]-1):
            for j in range(1, potential.shape[1]-1):
                if not mask[i, j]:
                    v_new = 0.25 * (
                        potential[i+1, j] + potential[i-1, j] +
                        potential[i, j+1] + potential[i, j-1]
                    )
                    potential[i, j] = vold[i, j] + omega * (v_new - vold[i, j])
        
        maxdiff = np.max(np.abs(potential - vold))
        
        if maxdiff < thresh:
            break
            
    return potential

if GPU_AVAILABLE:
    @cuda.jit
    def _solve_field_kernel(potential, new_potential, mask, omega):
        i, j = cuda.grid(2)
        
        if i > 0 and i < potential.shape[0]-1 and j > 0 and j < potential.shape[1]-1:
            if not mask[i, j]:
                v_new = 0.25 * (
                    potential[i+1, j] + potential[i-1, j] +
                    potential[i, j+1] + potential[i, j-1]
                )
                new_potential[i, j] = potential[i, j] + omega * (v_new - potential[i, j])
            else:
                new_potential[i, j] = potential[i, j]
                
    def solve_field_gpu(potential, mask, maxiter, thresh, omega=1.9):
        d_potential = cp.asarray(potential)
        d_mask = cp.asarray(mask)
        d_new_potential = cp.zeros_like(d_potential)
        
        threads_per_block = (16, 16)
        blocks_per_grid_x = (potential.shape[0] + threads_per_block[0] - 1) // threads_per_block[0]
        blocks_per_grid_y = (potential.shape[1] + threads_per_block[1] - 1) // threads_per_block[1]
        blocks_per_grid = (blocks_per_grid_x, blocks_per_grid_y)
        
        for _ in range(maxiter):
            d_new_potential[0, :] = d_potential[0, :]
            d_new_potential[-1, :] = d_potential[-1, :]
            d_new_potential[:, 0] = d_potential[:, 0]
            d_new_potential[:, -1] = d_potential[:, -1]
            
            _solve_field_kernel[blocks_per_grid, threads_per_block](
                d_potential, d_new_potential, d_mask, omega)
            
            maxdiff = cp.max(cp.abs(d_new_potential - d_potential))
            
            d_potential, d_new_potential = d_new_potential, d_potential
            
            if maxdiff < thresh:
                break
                
        return cp.asnumpy(d_potential)

def solve_field(potential, mask, maxiter, thresh, omega=1.9):
    if GPU_AVAILABLE:
        try:
            return solve_field_gpu(potential, mask, maxiter, thresh, omega)
        except Exception:
            return solve_field_cpu(potential, mask, maxiter, thresh, omega)
    else:
        return solve_field_cpu(potential, mask, maxiter, thresh, omega)

@nb.njit
def get_field(z, r, Ez, Er, size, dz, dr, nz, nr):
    if 0 <= z < size and 0 <= r < size:
        i = int(min(max(1, z / dz), nz - 2))
        j = int(min(max(1, r / dr), nr - 2))
        return Ez[i, j], Er[i, j]
    else:
        return 0.0, 0.0

@nb.njit
def calc_dynamics(z, r, vz, vr, Ez, Er, qm, mass, c):
    vsq = vz**2 + vr**2
    csq = c**2
    
    gamma = 1.0 / np.sqrt(1.0 - (vsq/csq))
    
    Fz = qm * mass * Ez
    Fr = qm * mass * Er
    
    factor = gamma/(gamma + 1.0) * (1.0/csq)
    vdotF = vz*Fz + vr*Fr
    
    az = Fz/(gamma * mass) - factor * vz * vdotF/mass
    ar = Fr/(gamma * mass) - factor * vr * vdotF/mass
    
    return np.array([vz, vr, az, ar])

class PotentialField:    
    def __init__(self, nz: float, nr: float, physical_size: float):
        self.nz = int(nz)
        self.nr = int(nr)
        self.size = physical_size
        self.dz = physical_size / self.nz
        self.dr = physical_size / self.nr
        self.potential = np.zeros((self.nz, self.nr))
        self.electrode_mask = np.zeros((self.nz, self.nr), dtype=bool)
        self.Ez = None
        self.Er = None
        self._field_cache = {}
    
    def add_electrode(self, config: ElectrodeConfig):
        start, width = int(config.start), int(config.width)
        ap_start, ap_width = int(config.ap_start), int(config.ap_width)
        outer_diameter = config.outer_diameter
        voltage = config.voltage
        
        ap_center = ap_start + ap_width / 2
        
        r_min = max(0, int(ap_center - outer_diameter / 2))
        r_max = min(int(ap_center + outer_diameter / 2), self.nr)
        
        self.potential[start:start+width, r_min:r_max] = voltage
        self.electrode_mask[start:start+width, r_min:r_max] = True
        
        if ap_width > 0:
            self.potential[start:start+width, ap_start:ap_start+ap_width] = 0
            self.electrode_mask[start:start+width, ap_start:ap_start+ap_width] = False
        
        self._field_cache.clear()
    
    def solve_potential(self, max_iterations: float = 2000, convergence_threshold: float = 1e-6):
        self.potential = solve_field(self.potential, self.electrode_mask, 
                                    int(max_iterations), convergence_threshold)
        
        if GPU_AVAILABLE:
            try:
                d_potential = cp.asarray(self.potential)
                d_Ez = cp.zeros_like(d_potential)
                d_Er = cp.zeros_like(d_potential)
                
                d_Ez[1:-1, :] = -(d_potential[2:, :] - d_potential[:-2, :]) / (2 * self.dz)
                d_Er[:, 1:-1] = -(d_potential[:, 2:] - d_potential[:, :-2]) / (2 * self.dr)
                
                d_Ez[0, :] = -(d_potential[1, :] - d_potential[0, :]) / self.dz
                d_Ez[-1, :] = -(d_potential[-1, :] - d_potential[-2, :]) / self.dz
                d_Er[:, 0] = -(d_potential[:, 1] - d_potential[:, 0]) / self.dr
                d_Er[:, -1] = -(d_potential[:, -1] - d_potential[:, -2]) / self.dr
                
                self.Ez = cp.asnumpy(d_Ez)
                self.Er = cp.asnumpy(d_Er)
            except Exception:
                self.Ez, self.Er = np.gradient(-self.potential, self.dz, self.dr)
        else:
            dz, dr = self.dz, self.dr
            Ez = np.zeros_like(self.potential)
            Er = np.zeros_like(self.potential)
            
            Ez[1:-1, :] = -(self.potential[2:, :] - self.potential[:-2, :]) / (2 * dz)
            Er[:, 1:-1] = -(self.potential[:, 2:] - self.potential[:, :-2]) / (2 * dr)
            
            Ez[0, :] = -(self.potential[1, :] - self.potential[0, :]) / dz
            Ez[-1, :] = -(self.potential[-1, :] - self.potential[-2, :]) / dz
            Er[:, 0] = -(self.potential[:, 1] - self.potential[:, 0]) / dr
            Er[:, -1] = -(self.potential[:, -1] - self.potential[:, -2]) / dr
            
            self.Ez, self.Er = Ez, Er
            
        self._field_cache.clear()
        return self.potential
    
    def get_field_at_position(self, z: float, r: float) -> Tuple[float, float]:
        cache_key = (z, r)
        if cache_key in self._field_cache:
            return self._field_cache[cache_key]
        
        result = get_field(z, r, self.Ez, self.Er, self.size, 
                          self.dz, self.dr, self.nz, self.nr)
        
        if len(self._field_cache) < 100000:
            self._field_cache[cache_key] = result
            
        return result

class ParticleTracer:
    ELECTRON_CHARGE = -1.602e-19 
    ELECTRON_MASS = 9.11e-31
    SPEED_OF_LIGHT = 299792458.0

    def __init__(self, potential_field: PotentialField):
        self.field = potential_field
        self.current_ion = {
            'symbol': 'e-',
            'atomic_number': 0,
            'mass': self.ELECTRON_MASS,
            'charge': self.ELECTRON_CHARGE,
            'charge_mass_ratio': self.ELECTRON_CHARGE / self.ELECTRON_MASS
        }
        self.q_m = self.current_ion['charge_mass_ratio']

    def set_ion(self, symbol: str = 'e-', charge_state: float = 1):
        if symbol == 'e-':
            self.current_ion = {
                'symbol': 'e-',
                'atomic_number': 0,
                'mass': self.ELECTRON_MASS,
                'charge': self.ELECTRON_CHARGE,
                'charge_mass_ratio': self.ELECTRON_CHARGE / self.ELECTRON_MASS
            }
        else:
            elem = element(symbol)
            
            isotope_mass = elem.mass
            electron_charge = 1.602e-19
            
            ion_charge = charge_state * electron_charge
            
            self.current_ion = {
                'symbol': f"{symbol}{'+' if charge_state > 0 else '-'}{abs(charge_state)}",
                'atomic_number': elem.atomic_number,
                'mass': isotope_mass * 1.66053906660e-27,
                'charge': ion_charge,
                'charge_mass_ratio': ion_charge / (isotope_mass * 1.66053906660e-27)
            }
        
        self.q_m = self.current_ion['charge_mass_ratio']
        return self

    def get_velocity_from_energy(self, energy_eV: float) -> float:
        energy_joules = energy_eV * 1.602e-19
        mass = self.current_ion['mass']
        rest_energy = mass * self.SPEED_OF_LIGHT**2
        total_energy = rest_energy + energy_joules
        return self.SPEED_OF_LIGHT * np.sqrt(1 - (rest_energy/total_energy)**2)

    def particle_dynamics(self, t: float, state: List[float]) -> List[float]:
        z, r, vz, vr = state
        Ez, Er = self.field.get_field_at_position(z, r)
        return calc_dynamics(
            z, r, vz, vr, 
            Ez, Er, 
            self.q_m, 
            self.current_ion['mass'], 
            self.SPEED_OF_LIGHT
        )

    def trace_trajectory(self, 
                   initial_position: Tuple[float, float],
                   initial_velocity: Tuple[float, float],
                   simulation_time: float,
                   method: str = 'BDF',
                   rtol: float = 1e-8,
                   atol: float = 1e-10,
                   max_step: float = None) -> dict:
        initial_state = [
            initial_position[0], 
            initial_position[1],
            initial_velocity[0], 
            initial_velocity[1]
        ]
    
        solution = solve_ivp(
            self.particle_dynamics,
            (0, simulation_time),
            initial_state,
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=max_step)
    
        return solution

class EinzelLens:
    def __init__(self, 
                position: float, 
                width: float, 
                aperture_center: float,
                aperture_width: float,
                outer_diameter: float,
                focus_voltage: float,
                electrode_thickness: float = 2):
        center_thickness = width - 2 * electrode_thickness
        
        self.electrode1 = ElectrodeConfig(
            start=position - width/2,
            width=electrode_thickness,
            ap_start=aperture_center - aperture_width/2,
            ap_width=aperture_width,
            outer_diameter=outer_diameter,
            voltage=0
        )
        
        self.electrode2 = ElectrodeConfig(
            start=position - width/2 + electrode_thickness,
            width=center_thickness,
            ap_start=aperture_center - aperture_width/2,
            ap_width=aperture_width,
            outer_diameter=outer_diameter,
            voltage=focus_voltage
        )
        
        self.electrode3 = ElectrodeConfig(
            start=position + width/2 - electrode_thickness,
            width=electrode_thickness,
            ap_start=aperture_center - aperture_width/2,
            ap_width=aperture_width,
            outer_diameter=outer_diameter,
            voltage=0 
        )
    
    def add_to_field(self, field: PotentialField):
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)

class IonOpticsSystem:
    
    def __init__(self, nr: float, nz: float, physical_size: float = 0.1):
        self.field = PotentialField(nz, nr, physical_size)
        self.tracer = ParticleTracer(self.field)
        self.elements = []
        
    def add_electrode(self, config: ElectrodeConfig):
        self.field.add_electrode(config)
        
    def add_einzel_lens(self, 
                       position: float, 
                       width: float, 
                       aperture_center: float,
                       aperture_width: float,
                       outer_diameter: float,
                       focus_voltage: float,
                       electrode_thickness: float = 2):
        lens = EinzelLens(
            position, width, aperture_center, aperture_width, 
            outer_diameter, focus_voltage, electrode_thickness
        )
        lens.add_to_field(self.field)
        self.elements.append(lens)
        
    def solve_fields(self):
        return self.field.solve_potential()
    
    def _simulate_particle(self, params):
        start_z, r_pos, angle, velocity_magnitude, simulation_time = params
        
        vz = velocity_magnitude * np.cos(angle)
        vr = velocity_magnitude * np.sin(angle)
        
        sol = self.tracer.trace_trajectory(
            initial_position=(start_z, r_pos),
            initial_velocity=(vz, vr),
            simulation_time=simulation_time
        )
        return sol

    def simulate_beam(self, energy_eV: float, start_z: float,
                         r_range: Tuple[float, float],
                         angle_range: tuple,
                         num_particles: float,
                         simulation_time: float):
        velocity_magnitude = self.tracer.get_velocity_from_energy(energy_eV)
        min_angle_rad = np.radians(angle_range[0])
        max_angle_rad = np.radians(angle_range[1])
        angles = np.linspace(min_angle_rad, max_angle_rad, int(num_particles))
        r_positions = np.linspace(r_range[0], r_range[1], int(num_particles))
        
        if mp.cpu_count() > 1:
            params = [
                (start_z, r_pos, angle, velocity_magnitude, simulation_time) 
                for r_pos, angle in zip(r_positions, angles)
            ]
            
            with mp.Pool() as pool:
                trajectories = pool.map(self._simulate_particle, params)
        else:
            trajectories = []
            for r_pos, angle in zip(r_positions, angles):
                vz = velocity_magnitude * np.cos(angle)
                vr = velocity_magnitude * np.sin(angle)
                
                sol = self.tracer.trace_trajectory(
                    initial_position=(start_z, r_pos),
                    initial_velocity=(vz, vr),
                    simulation_time=simulation_time
                )
                trajectories.append(sol)
    
        return trajectories
        
    def visualize_system(self, 
                       trajectories=None, 
                       r_limits=None,
                       figsize=(15, 6),
                       title="Electron Trajectories"):
        plt.figure(figsize=figsize)
        plt.title(title)
        
        if trajectories:
            colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
            for i, sol in enumerate(trajectories):
                z_traj = sol.y[0]
                r_traj = sol.y[1]
                plt.plot(z_traj, r_traj, lw=1.5, color=colors[i])
        
        plt.xlabel('z position (meters)')
        plt.ylabel('r position (meters)')
        plt.grid(True, alpha=0.3)
        
        if r_limits:
            plt.ylim(r_limits)
            
        plt.tight_layout()
        return plt.gcf()