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
    def __init__(self, nx: int, ny: int, physical_size: float):
        self.nx = nx
        self.ny = ny
        self.size = physical_size
        self.dx = physical_size / nx
        self.dy = physical_size / ny
        self.potential = np.zeros((nx, ny))
        self.electrode_mask = np.zeros((nx, ny), dtype=bool)
        self.Ex = None
        self.Ey = None
    
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
                
        self.Ex, self.Ey = np.gradient(-self.potential, self.dx, self.dy)
        
        return self.potential
    
    def get_field_at_position(self, x: float, y: float) -> Tuple[float, float]:
        if 0 <= x < self.size and 0 <= y < self.size:
            i = int(min(max(1, x / self.dx), self.nx - 2))
            j = int(min(max(1, y / self.dy), self.ny - 2))
            return self.Ex[i, j], self.Ey[i, j]
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
        x, y, vx, vy = state
        Ex, Ey = self.field.get_field_at_position(x, y)
        
        v = np.sqrt(vx**2 + vy**2)
        gamma = 1.0 / np.sqrt(1.0 - (v/self.SPEED_OF_LIGHT)**2)
        
        ax = self.q_m * Ex / gamma
        ay = self.q_m * Ey / gamma
        
        return [vx, vy, ax, ay]
    
    def trace_trajectory(self, 
                       initial_position: Tuple[float, float],
                       initial_velocity: Tuple[float, float],
                       simulation_time: float,
                       method: str = 'RK45',
                       rtol: float = 1e-8) -> dict:
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
            rtol=rtol
        )
        
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
    
    def __init__(self, nx: int, ny: int, physical_size: float = 0.1):
        self.field = PotentialField(nx, ny, physical_size)
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
    def simulate_beam(self, energy_eV: float, start_x: float,
                                y_range: Tuple[float, float],
                                  angle_range: tuple,
                                  num_particles: int,
                                  simulation_time: float):
        velocity_magnitude = self.tracer.get_velocity_from_energy(energy_eV)
        min_angle_rad = np.radians(angle_range[0])
        max_angle_rad = np.radians(angle_range[1])
        angles = np.linspace(min_angle_rad, max_angle_rad, num_particles)
        y_position = (y_range[0] + y_range[1]) / 2
        trajectories = []
        for angle in angles:
            vx = velocity_magnitude * np.cos(angle)
            vy = velocity_magnitude * np.sin(angle)
       
            sol = self.tracer.trace_trajectory(
                initial_position=(start_x, y_position),
                initial_velocity=(vx, vy),
                simulation_time=simulation_time
            )
            trajectories.append(sol)
       
        return trajectories
        
    def visualize_system(self, 
                       trajectories=None, 
                       y_limits=None,
                       figsize=(15, 6),
                       title="Electron Trajectories"):
        plt.figure(figsize=figsize)
        plt.title(title)
        
        if trajectories:
            colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
            for i, sol in enumerate(trajectories):
                x_traj = sol.y[0] / self.field.dx
                y_traj = sol.y[1] / self.field.dy
                plt.plot(x_traj, y_traj, lw=1.5, color=colors[i])
        
        plt.xlabel('x position (grid units)')
        plt.ylabel('y position (grid units)')
        plt.grid(True, alpha=0.3)
        
        if y_limits:
            plt.ylim(y_limits)
            
        plt.tight_layout()
        return plt.gcf()