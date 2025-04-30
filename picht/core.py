import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable


@dataclass
class ElectrodeConfig:
    start: int       
    width: int          
    ap_start: int        
    ap_width: int        
    voltage: float


class PotentialField:    
    def __init__(self, nr: int, nz: int, physical_size: float):
        self.nr = nr
        self.nz = nz
        self.size = physical_size
        self.dr = physical_size / nr
        self.dz = physical_size / nz
        self.potential = np.zeros((nr, nz))
        self.electrode_mask = np.zeros((nr, nz), dtype=bool)
        self.Er = None
        self.Ez = None
    
    def add_electrode(self, config: ElectrodeConfig):
        start, width = config.start, config.width
        ap_start, ap_width = config.ap_start, config.ap_width
        voltage = config.voltage
        
        self.potential[start:start+width, ap_start:ap_start+ap_width] = voltage
        self.electrode_mask[start:start+width, ap_start:ap_start+ap_width] = True
    
    def solve_potential(self, max_iterations: int = 2000, convergence_threshold: float = 1e-6):
        for _ in range(max_iterations):
            V_old = self.potential.copy()
            self.potential[1:-1, 1:-1] = 0.25 * (
                V_old[2:, 1:-1] + V_old[:-2, 1:-1] +
                V_old[1:-1, 2:] + V_old[1:-1, :-2]
            )
            
            self.potential[self.electrode_mask] = V_old[self.electrode_mask]
            
            if np.max(np.abs(self.potential - V_old)) < convergence_threshold:
                break
                
        self.Er, self.Ez = np.gradient(-self.potential, self.dr, self.dz)
        
        return self.potential
    
    def get_field_at_position(self, r: float, z: float) -> Tuple[float, float]:
        if 0 <= r < self.size and 0 <= z < self.size:
            i = int(min(max(1, r / self.dr), self.nr - 2))
            j = int(min(max(1, z / self.dz), self.nz - 2))
            return self.Er[i, j], self.Ez[i, j]
        else:
            return 0, 0


class ParticleTracer:
    ELECTRON_CHARGE = -1.602e-19 
    ELECTRON_MASS = 9.11e-31
    SPEED_OF_LIGHT = 299792458.0
    
    def __init__(self, potential_field: PotentialField):
        self.field = potential_field
        self.q_m = self.ELECTRON_CHARGE / self.ELECTRON_MASS
    
    def set_charge_mass_ratio(self, q: float, m: float):
        self.q_m = q / m
    
    def get_velocity_from_energy(self, energy_eV: float, mass: float = ELECTRON_MASS) -> float:
        energy_joules = energy_eV * 1.602e-19
        rest_energy = mass * self.SPEED_OF_LIGHT**2
        total_energy = rest_energy + energy_joules
        return self.SPEED_OF_LIGHT * np.sqrt(1 - (rest_energy/total_energy)**2)
    
    def particle_dynamics(self, t: float, state: List[float]) -> List[float]:
        r, z, vr, vz = state
        Er, Ez = self.field.get_field_at_position(r, z)
        
        v = np.sqrt(vr**2 + vz**2)
        
        gamma = 1.0 / np.sqrt(1.0 - (v/self.SPEED_OF_LIGHT)**2)
        
        ar = self.q_m * Er / gamma
        az = self.q_m * Ez / gamma
        
        return [vr, vz, ar, az]
    
    def trace_trajectory(self, 
                   initial_position: Tuple[float, float],
                   initial_velocity: Tuple[float, float],
                   simulation_time: float,
                   method: str = 'BDF',
                   rtol: float = 1e-8,
                   atol: float = 1e-10) -> dict:
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
            atol=atol)
    
        return solution


class EinzelLens:
    def __init__(self, 
                position: float, 
                width: float, 
                aperture_center: float,
                aperture_width: float,
                focus_voltage: float,
                electrode_thickness: float = 2):
        center_thickness = width - 2 * electrode_thickness
        
        self.electrode1 = ElectrodeConfig(
            start=int(position - width/2),
            width=electrode_thickness,
            ap_start=int(aperture_center - aperture_width/2),
            ap_width=int(aperture_width),
            voltage=0
        )
        
        self.electrode2 = ElectrodeConfig(
            start=int(position - width/2 + electrode_thickness),
            width=int(center_thickness),
            ap_start=int(aperture_center - aperture_width/2),
            ap_width=int(aperture_width),
            voltage=focus_voltage
        )
        
        self.electrode3 = ElectrodeConfig(
            start=int(position + width/2 - electrode_thickness),
            width=electrode_thickness,
            ap_start=int(aperture_center - aperture_width/2),
            ap_width=int(aperture_width),
            voltage=0 
        )
    
    def add_to_field(self, field: PotentialField):
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)


class IonOpticsSystem:
    
    def __init__(self, nr: int, nz: int, physical_size: float = 0.1):
        self.field = PotentialField(nr, nz, physical_size)
        self.tracer = ParticleTracer(self.field)
        self.elements = []
        
    def add_electrode(self, config: ElectrodeConfig):
        self.field.add_electrode(config)
        
    def add_einzel_lens(self, 
                       position: float, 
                       width: float, 
                       aperture_center: float,
                       aperture_width: float,
                       focus_voltage: float,
                       electrode_thickness: float = 2):
        lens = EinzelLens(
            position, width, aperture_center, aperture_width, 
            focus_voltage, electrode_thickness
        )
        lens.add_to_field(self.field)
        self.elements.append(lens)
        
    def solve_fields(self):
        return self.field.solve_potential()
    def simulate_beam(self, energy_eV: float, start_r: float,
                                z_range: Tuple[float, float],
                                  angle_range: tuple,
                                  num_particles: int,
                                  simulation_time: float):
        velocity_magnitude = self.tracer.get_velocity_from_energy(energy_eV)
        min_angle_rad = np.radians(angle_range[0])
        max_angle_rad = np.radians(angle_range[1])
        angles = np.linspace(min_angle_rad, max_angle_rad, num_particles)
        z_position = (z_range[0] + z_range[1]) / 2
        trajectories = []
        for angle in angles:
            vr = velocity_magnitude * np.cos(angle)
            vz = velocity_magnitude * np.sin(angle)
       
            sol = self.tracer.trace_trajectory(
                initial_position=(start_r, z_position),
                initial_velocity=(vr, vz),
                simulation_time=simulation_time
            )
            trajectories.append(sol)
       
        return trajectories
        
    def visualize_system(self, 
                       trajectories=None, 
                       z_limits=None,
                       figsize=(15, 6),
                       title="Electron Trajectories"):
        plt.figure(figsize=figsize)
        plt.title(title)
        
        if trajectories:
            colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
            for i, sol in enumerate(trajectories):
                r_traj = sol.z[0] / self.field.dr
                z_traj = sol.z[1] / self.field.dz
                r_plot = r_traj/dr
                z_plot = z_traj/dz
                plt.plot(r_traj, z_traj, lw=1.5, color=colors[i])
        
        plt.xlabel('radial position (meters)')
        plt.zlabel('z position (meters)')
        plt.grid(True, alpha=0.3)
        
        if z_limits:
            plt.zlim(z_limits)
            
        plt.tight_layout()
        return plt.gcf()