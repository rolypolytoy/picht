import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable


@dataclass
class ElectrodeConfig:
    """Configuration for a single electrode in an einzel lens or other electrostatic element"""
    start: int           # Starting position (grid units)
    width: int           # Width of electrode (grid units)
    ap_start: int        # Aperture starting position (grid units)
    ap_width: int        # Aperture width (grid units)
    voltage: float       # Electrode voltage (V)


class PotentialField:
    """Calculates and stores the electric potential field and resulting forces"""
    
    def __init__(self, nx: int, ny: int, physical_size: float):
        """
        Initialize the potential field calculator
        
        Args:
            nx: Number of grid points in x direction
            ny: Number of grid points in y direction
            physical_size: Physical size of the system in meters
        """
        self.nx = nx
        self.ny = ny
        self.size = physical_size
        self.dx = physical_size / nx
        self.dy = physical_size / ny
        
        # Initialize potential and mask arrays
        self.potential = np.zeros((nx, ny))
        self.electrode_mask = np.zeros((nx, ny), dtype=bool)
        
        # Fields will be calculated after solving
        self.Ex = None
        self.Ey = None
    
    def add_electrode(self, config: ElectrodeConfig):
        """Add an electrode to the system"""
        start, width = config.start, config.width
        ap_start, ap_width = config.ap_start, config.ap_width
        voltage = config.voltage
        
        self.potential[start:start+width, ap_start:ap_start+ap_width] = voltage
        self.electrode_mask[start:start+width, ap_start:ap_start+ap_width] = True
    
    def solve_potential(self, max_iterations: int = 2000, convergence_threshold: float = 1e-6):
        """Solve the potential field using finite difference method"""
        for _ in range(max_iterations):
            V_old = self.potential.copy()
            
            # Apply finite difference (Laplace equation solver)
            self.potential[1:-1, 1:-1] = 0.25 * (
                V_old[2:, 1:-1] + V_old[:-2, 1:-1] +
                V_old[1:-1, 2:] + V_old[1:-1, :-2]
            )
            
            # Maintain fixed values at electrodes
            self.potential[self.electrode_mask] = V_old[self.electrode_mask]
            
            # Check for convergence
            if np.max(np.abs(self.potential - V_old)) < convergence_threshold:
                break
                
        # Calculate electric field components
        self.Ex, self.Ey = np.gradient(-self.potential, self.dx, self.dy)
        
        return self.potential
    
    def get_field_at_position(self, x: float, y: float) -> Tuple[float, float]:
        """Get electric field at a specific position"""
        if 0 <= x < self.size and 0 <= y < self.size:
            i = int(min(max(1, x / self.dx), self.nx - 2))
            j = int(min(max(1, y / self.dy), self.ny - 2))
            return self.Ex[i, j], self.Ey[i, j]
        else:
            return 0, 0


class ParticleTracer:
    """Simulates charged particle trajectories through electric fields"""
    
    # Physical constants
    ELECTRON_CHARGE = -1.602e-19  # C
    ELECTRON_MASS = 9.11e-31      # kg
    
    def __init__(self, potential_field: PotentialField):
        """
        Initialize the particle tracer
        
        Args:
            potential_field: The electric potential field to trace particles through
        """
        self.field = potential_field
        self.q_m = self.ELECTRON_CHARGE / self.ELECTRON_MASS  # Default to electron
    
    def set_charge_mass_ratio(self, q: float, m: float):
        """Set custom charge and mass for the particle"""
        self.q_m = q / m
    
    @staticmethod
    def get_velocity_from_energy(energy_eV: float, mass: float = ELECTRON_MASS) -> float:
        """Convert energy in eV to velocity"""
        energy_joules = energy_eV * 1.602e-19
        return np.sqrt(2 * energy_joules / mass)
    
    def particle_dynamics(self, t: float, state: List[float]) -> List[float]:
        """State evolution function for the particle dynamics"""
        x, y, vx, vy = state
        Ex, Ey = self.field.get_field_at_position(x, y)
        return [vx, vy, self.q_m * Ex, self.q_m * Ey]
    
    def trace_trajectory(self, 
                       initial_position: Tuple[float, float],
                       initial_velocity: Tuple[float, float],
                       simulation_time: float,
                       method: str = 'RK45',
                       rtol: float = 1e-8) -> dict:
        """
        Trace a particle trajectory through the field
        
        Args:
            initial_position: Starting (x, y) position
            initial_velocity: Starting (vx, vy) velocity
            simulation_time: Total time to simulate
            method: Integration method for solve_ivp
            rtol: Relative tolerance for the solver
            
        Returns:
            Solution object from solve_ivp
        """
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
    """Representation of an einzel lens with three electrodes"""
    
    def __init__(self, 
                position: float, 
                width: float, 
                aperture_center: float,
                aperture_width: float,
                focus_voltage: float,
                electrode_thickness: float = 2):
        """
        Initialize an einzel lens
        
        Args:
            position: Center position of the lens (grid units)
            width: Total width of the lens assembly (grid units)
            aperture_center: Y-position of the aperture center (grid units)
            aperture_width: Width of the aperture (grid units)
            focus_voltage: Voltage for the center electrode (V)
            electrode_thickness: Thickness of the outer electrodes (grid units)
        """
        center_thickness = width - 2 * electrode_thickness
        
        # Create the three electrodes
        self.electrode1 = ElectrodeConfig(
            start=int(position - width/2),
            width=electrode_thickness,
            ap_start=int(aperture_center - aperture_width/2),
            ap_width=int(aperture_width),
            voltage=0  # Outer electrodes are grounded
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
            voltage=0  # Outer electrodes are grounded
        )
    
    def add_to_field(self, field: PotentialField):
        """Add this einzel lens to a potential field"""
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)


class IonOpticsSystem:
    """Complete ion optics system with multiple lens elements"""
    
    def __init__(self, nx: int, ny: int, physical_size: float = 0.1):
        """
        Initialize an ion optics system
        
        Args:
            nx: Number of grid points in x direction
            ny: Number of grid points in y direction
            physical_size: Physical size of the system in meters
        """
        self.field = PotentialField(nx, ny, physical_size)
        self.tracer = ParticleTracer(self.field)
        self.elements = []
        
    def add_electrode(self, config: ElectrodeConfig):
        """Add a custom electrode to the system"""
        self.field.add_electrode(config)
        
    def add_einzel_lens(self, 
                       position: float, 
                       width: float, 
                       aperture_center: float,
                       aperture_width: float,
                       focus_voltage: float,
                       electrode_thickness: float = 2):
        """Add an einzel lens to the system"""
        lens = EinzelLens(
            position, width, aperture_center, aperture_width, 
            focus_voltage, electrode_thickness
        )
        lens.add_to_field(self.field)
        self.elements.append(lens)
        
    def solve_fields(self):
        """Solve the potential field for the entire system"""
        return self.field.solve_potential()
    
    def simulate_beam(self, 
                    energy_eV: float,
                    start_x: float,
                    y_range: Tuple[float, float],
                    num_particles: int,
                    simulation_time: float):
        """
        Simulate a beam of particles through the system
        
        Args:
            energy_eV: Particle energy in eV
            start_x: Starting x position for all particles
            y_range: (min_y, max_y) for particle starting positions
            num_particles: Number of particles to simulate
            simulation_time: Duration of simulation in seconds
            
        Returns:
            List of trajectory solutions
        """
        velocity = self.tracer.get_velocity_from_energy(energy_eV)
        y_positions = np.linspace(y_range[0], y_range[1], num_particles)
        
        trajectories = []
        for y_pos in y_positions:
            sol = self.tracer.trace_trajectory(
                initial_position=(start_x, y_pos),
                initial_velocity=(velocity, 0),
                simulation_time=simulation_time
            )
            trajectories.append(sol)
            
        return trajectories
        
    def visualize_system(self, 
                       trajectories=None, 
                       y_limits=None,
                       figsize=(15, 6),
                       title="Ion Optics System Simulation"):
        """
        Visualize the system and particle trajectories
        
        Args:
            trajectories: List of trajectory solutions (optional)
            y_limits: (min_y, max_y) for plot y-axis (grid units)
            figsize: Figure size for matplotlib
            title: Plot title
        """
        plt.figure(figsize=figsize)
        plt.title(title)
        
        # Plot potential field
        plt.contourf(self.field.potential.T, levels=20, cmap='Blues', alpha=0.2)
        
        # Plot particle trajectories if provided
        if trajectories:
            colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
            for i, sol in enumerate(trajectories):
                x_traj = sol.y[0] / self.field.dx
                y_traj = sol.y[1] / self.field.dy
                plt.plot(x_traj, y_traj, lw=1.5, color=colors[i])
        
        # Plot electrodes
        for i in range(self.field.nx):
            for j in range(self.field.ny):
                if self.field.electrode_mask[i, j]:
                    plt.axvspan(i, i+1, ymin=j/self.field.ny, ymax=(j+1)/self.field.ny,
                              color='gray', alpha=0.4)
        
        plt.xlabel('x position (grid units)')
        plt.ylabel('y position (grid units)')
        plt.grid(True, alpha=0.3)
        
        if y_limits:
            plt.ylim(y_limits)
            
        plt.tight_layout()
        return plt.gcf()
