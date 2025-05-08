import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable
import numba as nb
from mendeleev import element

@dataclass
class ElectrodeConfig:
    """
    Configures the simplest kind of electrostatic lens: the cylindrical electrode.
    
    Attributes:
        start (float): Position where the electrode begins on the z-axis in grid units.
        width (float): Width of the electrode on the z-axis in grid units.
        ap_start (float): Starting position of the aperture on the r-axis in grid units.
        ap_width (float): Width of the aperture on the r-axis in grid units.
        outer_diameter (float): Full diameter of the electrode on the r-axis in grid units.
        voltage (float): Voltage of the electrode in volts.
    """
    start: float
    width: float
    ap_start: float
    ap_width: float
    outer_diameter: float
    voltage: float

@nb.njit
def solve_field(potential, mask, maxiter, thresh, omega=1.9):
    """
    Implements a numerical solution to Laplace's equation for electrostatic
    potential in cylindrical coordinates using the SOR technique with Dirichlet boundary
    conditions. It iteratively refines the potential field until convergence is achieved, as defined
    by the thresh argument.
    
    Parameters:
        potential (ndarray): 2D numpy array containing initial potential values.
        mask (ndarray): Boolean mask indicating electrode positions, True where electrodes are.
        maxiter (int): Maximum number of iterations to perform for the successive over relaxation solver.
        thresh (float): Convergence threshold for the solution.
        omega (float, optional): Relaxation parameter for SOR- omega defaults to 1.9 and must be kept between 1 and 2 otherwise SOR is invalid.
    
    Returns:
        ndarray: 2D numpy array containing the solved potential field.
        
    Note:
        This is a performance-critical function optimized with Numba's JIT compilation.
        Modifying the omega parameter or removing the @nb.njit tag may significantly
        impact performance.
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

@nb.njit
def get_field(z, r, Ez, Er, axial_size, radial_size, dz, dr, nz, nr):
    """
    Provides electric field values at fractional grid positions by
    interpolating between grid points. It handles boundary conditions by returning
    zero field values for positions outside the simulation domain, which means
    outside the defined domain, no forces can be imposed on the charged particles.
    
    Parameters:
        z (float): Z-axis position in meters.
        r (float): R-axis position in meters.
        Ez (ndarray): Z-component of the electric field array.
        Er (ndarray): R-component of the electric field array.
        axial_size (float): Total size of z-axis in meters.
        radial_size (float): Total size of r-axis in meters.
        dz (float): Grid spacing in z-axis, found by dividing nz/z.
        dr (float): Grid spacing in r-axis, found by dividing nr/r.
        nz (int): Number of grid points in z-axis.
        nr (int): Number of grid points in r-axis.
        
    Returns:
        tuple: (Ez, Er) Electric field components at the specified position in V/m.
    """
    if 0 <= z < axial_size and 0 <= r < radial_size:
        i = int(min(max(1, z / dz), nz - 2))
        j = int(min(max(1, r / dr), nr - 2))
        return Ez[i, j], Er[i, j]
    else:
        return 0.0, 0.0

@nb.njit
def calc_dynamics(z, r, vz, vr, Ez, Er, qm, mass, c):
    """    
    Calculates the acceleration of charged particles by applying
    the Lorentz force with special-relativistic corrections. It calculates the Lorentz factor
    through the variable gamma, and provides correct velocities and accelerations with
    minimal overhead.
    
    Parameters:
        z (float): Z-axis position in meters.
        r (float): R-axis position in meters.
        vz (float): Z-axis velocity in meters per second.
        vr (float): R-axis velocity in meters per second.
        Ez (float): Z-axis electric field in volts per meter.
        Er (float): R-axis electric field in volts per meter.
        qm (float): Charge-to-mass ratio of the particle.
        mass (float): Particle mass in kg.
        c (float): Speed of light in a vacuum in meters per second.
        
    Returns:
        ndarray: Array containing [vz, vr, az, ar], representing velocity and acceleration
                components in the z and r directions for complete kinematic information.
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
    This class handles initializing the electric potential field, adding electrodes and einzel lenses, 
    and solving for the electric field. Calculating the electric field is necessary at the start
    of the particle trajectory calculations, as this is the bulk of the information particles use
    for trajectory calculations.
    
    Attributes:
        nz (int): Number of grid points in the z-axis.
        nr (int): Number of grid points in the r-axis.
        axial_size (float): Size of the z-axis in meters.
        radial_size (float): Size of the r-axis in meters.
        dz (float): Grid spacing in z-axis per z-unit, in meters, found by axial_size/nz.
        dr (float): Grid spacing in r-axis per r-unit, in meters, found by radial_size/nr.
        potential (ndarray): 2D array with electric potential values in Volts.
        electrode_mask (ndarray): Boolean mask for electrode positions.
        Ez (ndarray): Z-component of electric field in Volts per meter.
        Er (ndarray): R-component of electric field in Volts per meter.
    """
    
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        """
        Initializes the potential field with specified dimensions. This parameterization is extremely
        important, as finer meshes increase accuracy at the cost of performance, and scales O(n^2)
        with increases in nz and nr.
        
        Parameters:
            nz (float): Number of grid points in the z-axis.
            nr (float): Number of grid points in the r-axis.
            axial_size (float): Physical length of the z-axis in meters.
            radial_size (float): Physical length of the r-axis in meters.
        """
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
    
    def add_electrode(self, config: ElectrodeConfig):
        """
        Adds a single electrode to the potential field.
        
        This method places an electrode with the specified configuration into the simulation
        domain, and handling the potential and field calculations cursorily.
        
        Parameters:
            config (ElectrodeConfig): Configuration parameters for the electrode.
        """
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
    
    def solve_potential(self, max_iterations: float = 2000, convergence_threshold: float = 1e-6):
        """
        Solves the electrostatic potential field using successive over-relaxation.
        
        This method computes the potential field throughout the full domain, then calculates the
        electric field using E = -∇V, using numpy's gradient-finding tool. Tweaking here isn't necessary
        for perf.
        
        Parameters:
            max_iterations (float, optional): Maximum number of iterations for the solver. Defaults to 2000.
            convergence_threshold (float, optional): Convergence criterion for the solution. Defaults to 1e-6.
            
        Returns:
            ndarray: The solved potential field.
        """
        self.potential = solve_field(self.potential, self.electrode_mask, int(max_iterations), convergence_threshold)
        self.Ez, self.Er = np.gradient(-self.potential, self.dz, self.dr)
        return self.potential
    
    def get_field_at_position(self, z: float, r: float) -> Tuple[float, float]:
        """
        Returns the electric field components at a specific position.
        
        Parameters:
            z (float): Position along the z-axis in meters.
            r (float): Position along the r-axis in meters.
            
        Returns:
            Tuple[float, float]: The electric field components (Ez, Er) at the specified position.
            
        Raises:
            ValueError: If called prior to solve_potential() being executed. You can't ask for field values
            before solving for the electric field. Reference the tutorial examples for proper syntax.
        """
        if self.Ez is None or self.Er is None:
            raise ValueError("Electric field components not calculated. Call solve_potential() first.")
            
        return get_field(z, r, self.Ez, self.Er, self.axial_size, self.radial_size, 
                         self.dz, self.dr, int(self.nz), int(self.nr))

class ParticleTracer:
    """
    Handles the computation of particle dynamics, including relativistic effects,
    for various ion species in the calculated electromagnetic fields. It supports tracking
    particles through the simulation domain and calculating their trajectories.
    
    Attributes:
        field (PotentialField): The potential field in which particles are traced.
        current_ion (dict): Information about the currently selected ion species.
        q_m (float): Charge-to-mass ratio of the current ion species.
        
    Constants:
        ELECTRON_CHARGE (float): Elementary charge in Coulombs.
        ELECTRON_MASS (float): Electron rest mass in kilograms.
        SPEED_OF_LIGHT (float): Speed of light in vacuum in meters per second.
    """
    
    ELECTRON_CHARGE = -1.602e-19 
    ELECTRON_MASS = 9.11e-31
    SPEED_OF_LIGHT = 299792458.0

    def __init__(self, potential_field: PotentialField):
        """
        Initializes the particle tracer with a potential field.
        
        Parameters:
            potential_field (PotentialField): The electric potential field in which particles will be traced.
        """
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
        """        
        Configures the particle tracer to simulate a specific element and charge.
        It integrates with the mendeleev library to get atomic species data and automatically calculates
        trajectories based on this.
        
        Parameters:
            symbol (str, optional): Chemical symbol of the element, or 'e-' for electrons. Defaults to 'e-' for backwards compatibility on versions prior to ion optics support.
            charge_state (float, optional): Charge state of the ion (positive for cations, negative for anions). Defaults to 1.
            
        Returns:
            ParticleTracer: Uses self-reference to allow for for method chaining and greater conciseness.
        """
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
        """
        Converts particle energy in electronvolts to velocity in meters per second,
        accounting for relativistic effects. It's accurate for all energy scales from single-digit eV to GeV.
        
        Parameters:
            energy_eV (float): Kinetic energy of the particle in electronvolts.
            
        Returns:
            float: Particle velocity in meters per second.
            
        Note:
            This method accounts for relativistic effects, making it accurate for high-energy particles
            where classical calculations would break down.
        """
        kinetic_energy = energy_eV * 1.602e-19
        mass = self.current_ion['mass']
        rest_energy = mass * self.SPEED_OF_LIGHT**2
        total_energy = rest_energy + kinetic_energy
        return self.SPEED_OF_LIGHT * np.sqrt(1 - (rest_energy/total_energy)**2)

    def particle_dynamics(self, t: float, state: List[float]) -> List[float]:
        """        
        Differential equation function used by the ODE solver to calculate particle trajectories.
        It calculates acceleration given position and velocity, and is relativistically compliant. A four-vector
        approach isn't used, so proper time doesn't need to be tracked.
        
        Parameters:
            t (float): Current time in the simulation.
            state (List[float]): Current state [z, r, vz, vr] with position and velocity components.
            
        Returns:
            List[float]: Derivatives of the state vector [vz, vr, az, ar] representing velocities and accelerations.
        """
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
                   atol: float = 1e-10) -> dict:
        """        
        Solves the equations of motion for a charged particle
        in the electric field, by using an ODE solver from scipy.
        
        Parameters:
            initial_position (Tuple[float, float]): Initial (z, r) position in meters.
            initial_velocity (Tuple[float, float]): Initial (vz, vr) velocity in meters per second.
            simulation_time (float): Total simulation time in seconds- should be between 1e-7 and 1e-10 as typical values.
            method (str, optional): Integration method for solve_ivp. Defaults to 'BDF', due to stiffness-related errors with 'RK45'.
            rtol (float, optional): Relative tolerance for the ODE solver. Defaults to 1e-8.
            atol (float, optional): Absolute tolerance for the ODE solver. Defaults to 1e-10.
            
        Returns:
            dict: Solution object from scipy.integrate.solve_ivp containing the trajectory information.
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
            rtol=rtol,
            atol=atol)
    
        return solution

class EinzelLens:
    """
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using the pre-existing
    ElectrodeConfig class in the archetypal arrangement of the unipotential lens.
    
    Attributes:
        electrode1 (ElectrodeConfig): Configures the first electrode at 0V.
        electrode2 (ElectrodeConfig): Configures for the second electrode held at (focus_voltage) V.
        electrode3 (ElectrodeConfig): Configuration for the third electrode at 0V.
    """
    
    def __init__(self, 
                position: float, 
                width: float, 
                aperture_center: float,
                aperture_width: float,
                outer_diameter: float,
                focus_voltage: float,
                gap_size: int = 1):
        """
        Creates a parameterizable unipotential/einzel lens with custom geometries, voltages, gap sizes, etc etc.
        
        Parameters:
            position (float): Position where the lens begins on the z-axis in grid units (dz).
            width (float): Width of the full lens assembly on the z-axis in grid units (dz).
            aperture_center (float): Center of the aperture on the r-axis in grid units (dr).
            aperture_width (float): Size of the aperture on the r-axis in grid units (dr).
            outer_diameter (float): Full diameter of the electrodes in grid units on the r-axis (dr).
            focus_voltage (float): Voltage applied to the center electrode in volts.
            gap_size (int, optional): Size of gaps between electrodes in grid units- very important parameter for fringing behavior. Defaults to 1.

        Note:
            If you want to add a focusing lens, make the polarity identical to the charge, and if you want a defocusing lens, make it the opposite. 
            E.g. for electrons a -10000V lens will be converging and a 10000V lens diverging, whereas for positive ions the inverse is true.
            You want unipotential lenses' voltage to be at, above, or near the eV value you give your particles in terms of speed. Note that, since 
            cylindrical electrodes aren't unipotential, you can have accelerating or decelerating behavior from them, but einzel lenses will
            keep your electrons at the same energy.
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
    
    def add_to_field(self, field: PotentialField):
        """        
        Parameters:
            field (PotentialField): An empty potential field initialization before adding electrodes to the field calculations.
        """
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)

class IonOpticsSystem:    
    """
    The top-level class which synthesizes all prior functionality.
    Initializes the potential field, particle tracing, and visualization
    components to provide a complete environment for designing and analyzing
    ion/electron optics systems. 
    
    Attributes:
        field (PotentialField): The potential field initialized in the simulation.
        tracer (ParticleTracer): The trajectory calcuation tracer.
        elements (list): List of all electrodes and lenses inside the system.
    """
    
    def __init__(self, nr: float, nz: float, axial_size: float = 0.1, radial_size: float = 0.1):
        """
        Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in meters, and with
        a grid of (nz, nr) units.
        
        Parameters:
            nr (float): Number of grid points in the r-axis direction.
            nz (float): Number of grid points in the z-axis direction.
            axial_size (float, optional): Length of the system in the z-direction in meters.
            radial_size (float, optional): Length of the system in the r-direction in meters.
        """
        self.field = PotentialField(nz, nr, axial_size, radial_size)
        self.tracer = ParticleTracer(self.field)
        self.elements = []

    def add_electrode(self, config: ElectrodeConfig):
        """
        Adds a single electrode to the ion optics system.
        
        Parameters:
            config (ElectrodeConfig): Configuration parameters for the electrode.
        """
        self.field.add_electrode(config)
        
    def add_einzel_lens(self, 
                       position: float, 
                       width: float, 
                       aperture_center: float,
                       aperture_width: float,
                       outer_diameter: float,
                       focus_voltage: float,
                       gap_size: int = 1):
        """        
        Inserts the einzel/unipotential lenses to the system, using the pre-existing EinzelLens class.
        
        Parameters:
            position (float): Position where the lens begins on the z-axis in grid units (dz).
            width (float): Width of the full lens assembly on the z-axis in grid units (dz).
            aperture_center (float): Center of the aperture on the r-axis in grid units (dr).
            aperture_width (float): Size of the aperture on the r-axis in grid units (dr).
            outer_diameter (float): Full diameter of the electrodes in grid units on the r-axis (dr).
            focus_voltage (float): Voltage applied to the center electrode in volts.
            gap_size (int, optional): Size of gaps between electrodes in grid units- very important parameter for fringing behavior. Defaults to 1.
        """
        lens = EinzelLens(
            position, width, aperture_center, aperture_width, 
            outer_diameter, focus_voltage, gap_size
        )
        lens.add_to_field(self.field)
        self.elements.append(lens)
        
    def solve_fields(self):
        """
        Solves the potential field, calculates all electric field components and properly
        initializes things like electrode masks and boundary conditions.
        Call this after all electrodes and lenses have been added
        to the system, but before simulating particle trajectories.
        
        Returns:
            ndarray: The solved potential field.
        """
        return self.field.solve_potential()

    def simulate_beam(self, energy_eV: float, start_z: float,
                         r_range: Tuple[float, float],
                         angle_range: tuple,
                         num_particles: float,
                         simulation_time: float):
        """                
        Initializes a parameterized beam, where you can define each particle's
        kinetic energy in eV, its angular spread, and how many individual particles you
        simulate. Note that we don't implement Coulomb repulsion between particles, so the
        number of particles is only for better, more granular understanding.
        
        Parameters:
            energy_eV (float): Kinetic energy of each particle in electronvolts.
            start_z (float): Starting position for all particles on the z-axis in meters, not in grid units.
            r_range (Tuple[float, float]): The range of initial radial positions in meters- effectively the beam width.
            angle_range (tuple): Range of initial angles from the horizontal in radians.
            num_particles (float): Number of particles to simulate in the beam- keep at 3-6 for prototyping and 10-100 for full simulation.
            simulation_time (float): Total simulation time in seconds, typically between 1e-10 to 1e-7.
            
        Returns:
            list: List of trajectory solutions for all simulated particles.            
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
        
    def visualize_system(self, 
                       trajectories=None, 
                       r_limits=None,
                       figsize=(15, 6),
                       title="Electron Trajectories"):
        """
        Visualizes the charged particle trajectories, by creating a plot in
        cylindrical coordinates.
        
        Parameters:
            trajectories (list, optional): List of trajectory solutions to visualize. None by default.
            r_limits (tuple, optional): Y-axis limits (min, max) for the plot in meters. None by default.
            figsize (tuple, optional): Figure size as (width, height) in inches. (15, 6) by default.
            title (str, optional): Plot title, "Electron Trajectories" by default.
            
        Returns:
            matplotlib.figure.Figure: The created figure object.

        Example Usage:
        figure = system.visualize_system(
        trajectories=trajectories,
        r_limits = (0.049, 0.051))

        This forcefully limits the view of the trajectories to be within 0.049 meters and 0.051 meters
        on the radial axis. If you don't specify r_limits, it's completely fine, because it auto-sizes
        it based on how divergent the beams are. Specify r_limits to 'zoom in' to better understand 
        focusing, or ultra-fine beam behavior.
        """
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