import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable
import numba as nb
from mendeleev import element

@dataclass
class ElectrodeConfig:
    start: float
    width: float
    ap_start: float
    ap_width: float
    outer_diameter: float
    voltage: float

@nb.njit
def solve_field(potential, mask, maxiter, thresh, omega=1.9):
    """
    This solves the electrostatic potential using successive over relaxation. This
    codebase previously used the Jacobi method but we switched to SOR due to greater
    performance without that much overhead.
    
    We implement a numerical solution to Laplace's equation for electrostatic
    potential in cylindrical coordinates using the SOR technique with Dirichlet boundary
    conditions.

    I highly recommend not adjusting anything here, including the value for omega, or
    removing the @nb.njit tag, because this is a performance-critical function- perhaps
    the most performance-critical one.
    """
    for _ in range(maxiter):
        vold = potential.copy()
        
        for i in range(1, potential.shape[0]-1):
            for j in range(1, potential.shape[1]-1):
                if not mask[i, j]:
                    v_new = 0.25 * (
                        potential[i+1, j] + potential[i-1, j] +
                        potential[i, j+1] + potential[i, j-1]
                    )
                    potential[i, j] = vold[i, j] + omega * (v_new - vold[i, j])
        
        maxdiff = 0.0
        for i in range(potential.shape[0]):
            for j in range(potential.shape[1]):
                diff = abs(potential[i, j] - vold[i, j])
                if diff > maxdiff:
                    maxdiff = diff
        
        if maxdiff < thresh:
            break
            
    return potential

#retrieve electric field components at a given position using interpolation to obtain electric field values at fractional grid positions
@nb.njit
def get_field(z, r, Ez, Er, axial_size, radial_size, dz, dr, nz, nr):
    if 0 <= z < axial_size and 0 <= r < radial_size:
        i = int(min(max(1, z / dz), nz - 2))
        j = int(min(max(1, r / dr), nr - 2))
        return Ez[i, j], Er[i, j]
    else:
        return 0.0, 0.0

@nb.njit
def calc_dynamics(z, r, vz, vr, Ez, Er, qm, mass, c):
    """    
    Computes relativistic dynamics of charged particles using
    the Lorentz factor (special relativity) and the Lorentz force 
    for electrostatics. You can swap the Lorentz force handling
    with the paraxial ray equation for easy speedups- I've not 
    done so intentionally to maintain more physical realism.
    
    Arguments:
        z: Z-axis position in meters
        r: R-axis position in meters
        vz: Z-axis velocity in meters per second
        vr: R-axis velocity in meters per second
        Ez: Z-axis electric field in volts per meter
        Er: R-axis electric field in volts per meter
        qm: charge/mass ratio of the particle
        mass: particle mass- self explanatory
        c: speed of light in a vacuum
        
    Returns:
        ndarray: Array containing [vz, vr, az, ar], for complete kinematic and kinetic information about the particle.

    Feel free to modify this. What the code does here is intuitive, and if you want to simulate non-relativistic physics simply force gamma = 1.0.
    """
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
    """
    Attributes:
        nz: Grid points in the z-axis
        nr: Grid points in the r-axis
        axial_size: Size of the z-axis in meters. 
        radial_size: Size of the r-axis in meters. 

        Note: In a lot of places in Picht we use "grid units", where one unit in the z-axis is found by dz = axial_size/nz,
        and one unit in the r-axis is found by dr = radial_size/nr. Remember this when initializing electrodes or einzel lenses.

        dz: Grid spacing in z-axis
        dr: Grid spacing in r-axis
        potential: 2D array with electric potential values
        electrode_mask: Boolean mask for electrode positions
        Ez: Z-component of electric field
        Er: R-component of electric field
    """
    
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        self.nz = nz
        self.nr = nr
        self.axial_size = axial_size
        self.radial_size = radial_size
        self.dz = axial_size / nz
        self.dr = radial_size / nr
        self.potential = np.zeros((int(nz), int(nr)))
        self.electrode_mask = np.zeros((int(nz), int(nr)), dtype=bool)
        self.Ez = None
        self.Er = None
    
    #single electrode initialization
    def add_electrode(self, config: ElectrodeConfig):
        start, width = config.start, config.width
        ap_start, ap_width = config.ap_start, config.ap_width
        outer_diameter = config.outer_diameter
        voltage = config.voltage
        
        ap_center = ap_start + ap_width / 2
        
        r_min = max(0, ap_center - outer_diameter / 2)
        r_max = min(ap_center + outer_diameter / 2, self.nr)
        
        self.potential[int(start):int(start+width), int(r_min):int(r_max)] = voltage
        self.electrode_mask[int(start):int(start+width), int(r_min):int(r_max)] = True
        
        if ap_width > 0:
            self.potential[int(start):int(start+width), int(ap_start):int(ap_start+ap_width)] = 0
            self.electrode_mask[int(start):int(start+width), int(ap_start):int(ap_start+ap_width)] = False

    #electrostatic potential solver using successive over relaxation with parameterizable convergence thresholds and repetition. 
    #tweak as needed to modify performance/accuracy tradeoffs    
    def solve_potential(self, max_iterations: float = 2000, convergence_threshold: float = 1e-6):
        self.potential = solve_field(self.potential, self.electrode_mask, int(max_iterations), convergence_threshold)
        self.Ez, self.Er = np.gradient(-self.potential, self.dz, self.dr)
        return self.potential
    
    def get_field_at_position(self, z: float, r: float) -> Tuple[float, float]:
        return get_field(z, r, self.Ez, self.Er, self.axial_size, self.radial_size, 
                         self.dz, self.dr, int(self.nz), int(self.nr))

#uses the particle dynamics functions previously initialized to trace particle dynamics. fairly boilerplate/convergence-like
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

    #code for calling up information about ions based on what you specify. defaults to electrons
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

    #relativistic velocity calculator to support keV, MeV, GeV and above energies
    def get_velocity_from_energy(self, energy_eV: float) -> float:
        kinetic_energy = energy_eV * 1.602e-19
        mass = self.current_ion['mass']
        rest_energy = mass * self.SPEED_OF_LIGHT**2
        total_energy = rest_energy + kinetic_energy
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

    #numerical integration solver (using scipy) to solve ODEs for particle trajectories
    #uses BDF instead of RK45 due to better support for stiff problems, which electrodynamics is
    #if you want to modify speed, increase rtol and atol for better speeds or reduce them for better accuracy
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

#einzel/unipotential lens class. probably the magnum opus of this library due to its generality and usability
class EinzelLens:
    def __init__(self, 
                position: float, 
                width: float, 
                aperture_center: float,
                aperture_width: float,
                outer_diameter: float,
                focus_voltage: float,
                gap_size: int):
        """
        Arguments:
            position: Position where the lens begins on the z-axis in grid units
            width: Width of the full lens assembly on the z-axis in grid units
            aperture_center: Center of the aperture of all the einzel lenses on the r-axis in grid units
            aperture_width: Size of the aperture on the r-axis in grid units
            outer_diameter: Top-to-bottom diameter of the full electrode in the r-axis, in grid units
            focus_voltage: Voltage applied to the center electrode in volts. Make it the same polarity as the charge of the particle to focus, and the opposite to defocus
            gap_size: Size of gaps between electrodes on the z-axis in grid units. Important for fringing fields but defaults to a reasonable value so it's optional.
        """
        electrode_thickness = (width - 3 * gap_size)/3.0 
        
        self.electrode1 = ElectrodeConfig(
            start=position,
            width=electrode_thickness,
            ap_start=aperture_center - aperture_width/2,
            ap_width=aperture_width,
            outer_diameter=outer_diameter,
            voltage=0
        )
        
        self.electrode2 = ElectrodeConfig(
            start=position + electrode_thickness + gap_size,
            width=electrode_thickness,
            ap_start=aperture_center - aperture_width/2,
            ap_width=aperture_width,
            outer_diameter=outer_diameter,
            voltage=focus_voltage
        )
        
        self.electrode3 = ElectrodeConfig(
            start=position + 2 * electrode_thickness + 2 * gap_size,
            width=electrode_thickness,
            ap_start=aperture_center - aperture_width/2,
            ap_width=aperture_width,
            outer_diameter=outer_diameter,
            voltage=0 
        )
    
    #computes the electric fields of the three electrodes in an einzel lens
    def add_to_field(self, field: PotentialField):
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)


#this class can be thought of as the main, or top-level class. This is where you can actually initialize the beams, what runs the simulations, and what enables its visualization. Most of the code here is boilerplate or connective in nature.
class IonOpticsSystem:    
    def __init__(self, nr: float, nz: float, axial_size: float = 0.1, radial_size: float = 0.1):
        self.field = PotentialField(nz, nr, axial_size, radial_size)
        self.tracer = ParticleTracer(self.field)
        self.elements = []

    #boilerplate code to add a lone electrode to the system using the function
    def add_electrode(self, config: ElectrodeConfig):
        self.field.add_electrode(config)
        
    def add_einzel_lens(self, 
                       position: float, 
                       width: float, 
                       aperture_center: float,
                       aperture_width: float,
                       outer_diameter: float,
                       focus_voltage: float,
                       gap_size: int = 1):

        #boilerplate code to initialize an Einzel/Unipotential lens in the system using the pre-existing class
        lens = EinzelLens(
            position, width, aperture_center, aperture_width, 
            outer_diameter, focus_voltage, gap_size
        )
        lens.add_to_field(self.field)
        self.elements.append(lens)
        
    def solve_fields(self):
        return self.field.solve_potential()

    def simulate_beam(self, energy_eV: float, start_z: float,
                         r_range: Tuple[float, float],
                         angle_range: tuple,
                         num_particles: float,
                         simulation_time: float):
        """        
        Arguments:
            energy_eV: Kinetic energy of each particle in the beam in electronvolts, the standard unit of particle physics
            start_z: Starting position for all particles on the z-axis. Importantly, this is in meters, not in grid units
            r_range: The range of initial beam values radially- this determines how large the diameter of your beam is
            angle_range: Range of initial angles from the horizontal in radians, not degrees. It accepts negative 
            num_particles: Number of electrons/ions you want to simulate. I'd recommend 3-6 for prototyping and 10-100 for final visualization. Massively affects how long it takes to compute.
            simulation_time: Total simulation time in seconds. Usually should be within 1e-10 to 1e-7. 
            
        Returns:
            list: Solved trajectories for all individual charged particles, based on the field initialized using successive over relaxation
        """
        velocity_magnitude = self.tracer.get_velocity_from_energy(energy_eV)
        min_angle_rad = np.radians(angle_range[0])
        max_angle_rad = np.radians(angle_range[1])
        angles = np.linspace(min_angle_rad, max_angle_rad, int(num_particles))
        r_positions = np.linspace(r_range[0], r_range[1], int(num_particles))
    
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
        
    #uses matplotlib to visualize electron dynamics with (z, r) where z is the axial direction and r is the radial
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