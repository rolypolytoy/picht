import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable
import numba as nb
from mendeleev import element
<<<<<<< Updated upstream
=======
import pyamg
from joblib import Parallel, delayed
import multiprocessing
import os
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

@dataclass
class ElectrodeConfig:
    """
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    Configures the simplest kind of electrostatic lens: the cylindrical electrode.
=======
    Configures a simple cylindrically symmetric electrostatic lens.
>>>>>>> Stashed changes
=======
    Configures a simple cylindrically symmetric electrostatic lens.
>>>>>>> Stashed changes
=======
    Configures a simple cylindrically symmetric electrostatic lens.
>>>>>>> Stashed changes
=======
    Configures a simple cylindrically symmetric electrostatic lens.
>>>>>>> Stashed changes
    
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    interpolating between grid points. It handles boundary conditions by returning
    zero field values for positions outside the simulation domain, which means
    outside the defined domain, no forces can be imposed on the charged particles.
=======
    picking the value at the nearest neighbor. Necessary because 
    FDM is discretizing by nature.
>>>>>>> Stashed changes
=======
    picking the value at the nearest neighbor. Necessary because 
    FDM is discretizing by nature.
>>>>>>> Stashed changes
=======
    picking the value at the nearest neighbor. Necessary because 
    FDM is discretizing by nature.
>>>>>>> Stashed changes
=======
    picking the value at the nearest neighbor. Necessary because 
    FDM is discretizing by nature.
>>>>>>> Stashed changes
    
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
<<<<<<< Updated upstream
def calc_dynamics(z, r, vz, vr, Ez, Er, qm, mass, c):
    """    
    Calculates the acceleration of charged particles by applying
    the Lorentz force with special-relativistic corrections. It calculates the Lorentz factor
    through the variable gamma, and provides correct velocities and accelerations with
    minimal overhead.
=======
def calc_dynamics(z, r, pz, pr, Ez, Er, q, m, c):
    """    
    Calculates the acceleration of charged particles by applying
    the Lorentz force with special-relativistic corrections. Uses energy
    momentum formalism for full eV to TeV support.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    
    Parameters:
        z (float): Z-axis position in meters.
        r (float): R-axis position in meters.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        vz (float): Z-axis velocity in meters per second.
        vr (float): R-axis velocity in meters per second.
        Ez (float): Z-axis electric field in volts per meter.
        Er (float): R-axis electric field in volts per meter.
        qm (float): Charge-to-mass ratio of the particle.
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        pz (float): Z-axis momentum in kg * m/s.
        pr (float): R-axis momentum in kg * m/s.
        Ez (float): Z-axis electric field in volts per meter.
        Er (float): R-axis electric field in volts per meter.
        q (float): Charge of the particle.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        mass (float): Particle mass in kg.
        c (float): Speed of light in a vacuum in meters per second.
        
    Returns:
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        ndarray: Array containing [vz, vr, dpz_dt, dpr_dt], representing velocity and force
                components in the z and r directions for complete kinematic information (force is dp/dt).
    """
    p_sq = pz**2 + pr**2
    E = np.sqrt((p_sq * c**2) + (m * c**2)**2)
    vz = pz * c**2 / E
    vr = pr * c**2 / E
    dpz_dt = q * Ez
    dpr_dt = q * Er
    return np.array([vz, vr, dpz_dt, dpr_dt])

class PotentialField:
    """
    This class handles everything related to the electric potential field, including
    adding electrodes and einzel lenses, and solving for the electric field.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        """
        Initializes the potential field with specified dimensions. This parameterization is extremely
        important, as finer meshes increase accuracy at the cost of performance, and scales O(n^2)
        with increases in nz and nr.
=======
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        """
        Initializes the potential field with the specific dimensions. This decision is 
        extremely important, as finer meshes increase accuracy at the cost of performance.
>>>>>>> Stashed changes
=======
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        """
        Initializes the potential field with the specific dimensions. This decision is 
        extremely important, as finer meshes increase accuracy at the cost of performance.
>>>>>>> Stashed changes
=======
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        """
        Initializes the potential field with the specific dimensions. This decision is 
        extremely important, as finer meshes increase accuracy at the cost of performance.
>>>>>>> Stashed changes
=======
    def __init__(self, nz: float, nr: float, axial_size: float, radial_size: float):
        """
        Initializes the potential field with the specific dimensions. This decision is 
        extremely important, as finer meshes increase accuracy at the cost of performance.
>>>>>>> Stashed changes
        
        Parameters:
            nz (float): Number of grid points in the z-axis.
            nr (float): Number of grid points in the r-axis.
            axial_size (float): Physical length of the z-axis in meters.
            radial_size (float): Physical length of the r-axis in meters.
        """
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        self.nz = nz
        self.nr = nr
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        self.nz = int(nz)
        self.nr = int(nr)
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        Adds a single electrode to the potential field.
        
        This method places an electrode with the specified configuration into the simulation
        domain, and handling the potential and field calculations cursorily.
        
=======
        Adds a single electrode to the electric field and handles all
        necessary calculations.
                
>>>>>>> Stashed changes
=======
        Adds a single electrode to the electric field and handles all
        necessary calculations.
                
>>>>>>> Stashed changes
=======
        Adds a single electrode to the electric field and handles all
        necessary calculations.
                
>>>>>>> Stashed changes
=======
        Adds a single electrode to the electric field and handles all
        necessary calculations.
                
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    
    def solve_potential(self, max_iterations: float = 2000, convergence_threshold: float = 1e-6):
        """
        Solves the electrostatic potential field using successive over-relaxation.
        
        This method computes the potential field throughout the full domain, then calculates the
        electric field using E = -∇V, using numpy's gradient-finding tool. Tweaking here isn't necessary
        for perf.
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

    def build_laplacian_matrix(self, mask, dirichlet_values=None):
        """
        Builds a sparse matrix for the Laplacian ∇²V = 0, 
        and implements Dirichlet (0 Volts) boundary conditions at all boundaries,
        to simulate a grounded metal boundary. First uses List of Lists matrices for
        the linear system, and converts to Compressed Sparse Row for rapid solving
        with PyAMG. 

        Parameters:
            mask (ndarray): Boolean array, true where electrodes exist, false elsewhere.
            dirichlet_values (ndarray): Electric potential values where mask is True.

        Returns:
            A (scipy.sparse.csr_matrix): A sparse matrix representation of the Laplacian operator
            on voltage. Its shape is (nz * nr, nz * nr). 
            b (ndarray): Contains the potential values at electrode positions (where mask is True)
            and 0 elsewhere. 
        """
        n = self.nz * self.nr
        A = lil_matrix((n, n))
        b = np.zeros(n)

        def idx(i, j):
            return i * self.nr + j

        for i in range(self.nz):
            for j in range(self.nr):
                k = idx(i, j)
                if mask[i, j]:
                    A[k, k] = 1.0
                    if dirichlet_values is not None:
                        b[k] = dirichlet_values[i, j]
                    
                elif i == 0 or i == self.nz - 1 or j == 0 or j == self.nr - 1:
                    A[k, k] = 1.0
                    b[k] = 0.0
                else:
                    A[k, k] = -4.0
                    A[k, idx(i-1, j)] = 1.0
                    A[k, idx(i+1, j)] = 1.0
                    A[k, idx(i, j-1)] = 1.0
                    A[k, idx(i, j+1)] = 1.0
        return A.tocsr(), b

    def solve_potential(self, max_iterations: float = 500, convergence_threshold: float = 1e-6):
        """
        Solves the electrostatic potential field using Multigrid methods with PyAMG. O(N) complexity
        with doubling in nz AND nr values- so a 1000x1000 grid costs only 10 times more time to solve
        than a 100x100 one. 

        This method first creates a CSR matrix for the Laplacian values for the potential field, uses PyAMG
        to actually solve for Laplace's equation ∇²V = 0. Then, it finds the gradient E = -∇V and thus
        finds the electric field.
        
        Parameters:
            max_iterations (float, optional): Maximum number of iterations for the solver. Defaults to 2000.
            convergence_threshold (float, optional): Convergence criterion for the solution. Defaults to 1e-6.
            
        Returns:
            ndarray: The solved potential field.
        """
        A, b = self.build_laplacian_matrix(self.electrode_mask, self.potential)
        
        scale_factor = 1.0 / max(np.max(np.abs(A.data)), 1e-10)
        A_scaled = A * scale_factor
        b_scaled = b * scale_factor
        
        ml = pyamg.smoothed_aggregation_solver(
            A_scaled, 
            max_coarse=10,
            strength='symmetric',
            smooth='jacobi',
            improve_candidates=None
        )
        
        x0 = np.zeros_like(b_scaled)
        x = ml.solve(b_scaled, x0=x0, tol=convergence_threshold, maxiter=int(max_iterations))
        self.potential = x.reshape((self.nz, self.nr))
>>>>>>> Stashed changes
        
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            
        Raises:
            ValueError: If called prior to solve_potential() being executed. You can't ask for field values
            before solving for the electric field. Reference the tutorial examples for proper syntax.
        """
        if self.Ez is None or self.Er is None:
            raise ValueError("Electric field components not calculated. Call solve_potential() first.")
            
=======
        """
>>>>>>> Stashed changes
=======
        """
>>>>>>> Stashed changes
=======
        """
>>>>>>> Stashed changes
=======
        """
>>>>>>> Stashed changes
        return get_field(z, r, self.Ez, self.Er, self.axial_size, self.radial_size, 
                         self.dz, self.dr, int(self.nz), int(self.nr))

class ParticleTracer:
    """
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    Handles trajectory dynamics calculations and visualizations.
    
    Attributes:
        field (PotentialField): The potential field class, which has most of the information required
        to calculate particle trajectories.
        current_ion (dict): A dictionary with selected information about the charged particle.
        q_m (float): Charge-to-mass ratio of the charged particle.
        
    Constants:
        ELECTRON_CHARGE (float): Electron's charge in Coulombs.
        ELECTRON_MASS (float): Electron mass in kilograms.
        SPEED_OF_LIGHT (float): Speed of light in vacuum in meters per second.
    """
    ELECTRON_CHARGE = -1.60217663e-19 
    ELECTRON_MASS = 9.1093837e-31
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    SPEED_OF_LIGHT = 299792458.0

    def __init__(self, potential_field: PotentialField):
        """
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        Initializes the particle tracer with a potential field.
        
        Parameters:
            potential_field (PotentialField): The electric potential field in which particles will be traced.
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        Initializes a potential field with the class.
        
        Parameters:
            potential_field (PotentialField): The electric potential field required to calculate particle dynamics.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        Configures the particle tracer to simulate a specific element and charge.
        It integrates with the mendeleev library to get atomic species data and automatically calculates
        trajectories based on this.
        
        Parameters:
            symbol (str, optional): Chemical symbol of the element, or 'e-' for electrons. Defaults to 'e-' for backwards compatibility on versions prior to ion optics support.
            charge_state (float, optional): Charge state of the ion (positive for cations, negative for anions). Defaults to 1.
            
        Returns:
            ParticleTracer: Uses self-reference to allow for for method chaining and greater conciseness.
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        Configures the particle tracer to pick a specific charged particle. 
        Integrates with the Mendeleev library to automatically retrieve
        information about any and all atoms, and natively supports electrons.

        Parameters:
            symbol (str, optional): Chemical symbol of the element, or 'e-' for electrons. Defaults to 'e-' for backwards compatibility.
            charge_state (float, optional): Charge of the ion, defaults to 1.
            
        Returns:
            ParticleTracer: Uses self-reference for method chaining to support Picht's general style of power but conciseness.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            electron_charge = 1.602e-19
            
=======
            electron_charge = 1.60217663e-19
>>>>>>> Stashed changes
=======
            electron_charge = 1.60217663e-19
>>>>>>> Stashed changes
=======
            electron_charge = 1.60217663e-19
>>>>>>> Stashed changes
=======
            electron_charge = 1.60217663e-19
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            energy_eV (float): Kinetic energy of the particle in electronvolts.
            
        Returns:
            float: Particle velocity in meters per second.
            
        Note:
            This method accounts for relativistic effects, making it accurate for high-energy particles
            where classical calculations would break down.
        """
        kinetic_energy = energy_eV * 1.602e-19
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
            energy_eV (float): Kinetic energy of the particle in electronvolts, the standard unit of energy in
            particle physics. 1eV is approximately 1.6022e-19 Joules, and is the kinetic energy an electron
            has after being accelerated through an electric field with potential difference of 1 Volt.
            
        Returns:
            float: Particle velocity in meters per second.
        """
        kinetic_energy = energy_eV * 1.60217663e-19
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            state (List[float]): Current state [z, r, vz, vr] with position and velocity components.
            
        Returns:
            List[float]: Derivatives of the state vector [vz, vr, az, ar] representing velocities and accelerations.
        """
        z, r, vz, vr = state
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
            state (List[float]): Current state [z, r, pz, pr] with position and momentum components.
            
        Returns:
            List[float]: Derivatives of the state vector [vz, vr, dpz_dt, dpr_dt] representing velocities and force components.
        """
        z, r, pz, pr = state
>>>>>>> Stashed changes
        Ez, Er = self.field.get_field_at_position(z, r)
        return calc_dynamics(
            z, r, vz, vr, 
            Ez, Er, 
            self.q_m, 
            self.current_ion['mass'], 
            self.SPEED_OF_LIGHT
        )

    def trace_trajectory(self, 
<<<<<<< Updated upstream
                   initial_position: Tuple[float, float],
                   initial_velocity: Tuple[float, float],
                   simulation_time: float,
                   method: str = 'BDF',
                   rtol: float = 1e-8,
                   atol: float = 1e-10) -> dict:
        """        
        Solves the equations of motion for a charged particle
        in the electric field, by using an ODE solver from scipy.
=======
                    initial_position: Tuple[float, float],
                    initial_velocity: Tuple[float, float],
                    simulation_time: float,
                    method: str = 'BDF',
                    rtol: float = 1e-9,
                    atol: float = 1e-12) -> dict:
        """        
        Solves the equations of motion for a charged particle
        in the electric field, by using an ODE solver from scipy.
        
        Parameters:
            initial_position (Tuple[float, float]): Initial (z, r) position in meters.
            initial_velocity (Tuple[float, float]): Initial (vz, vr) velocity in meters per second.
            simulation_time (float): Total simulation time in seconds- should be between 1e-7 and 1e-10 as typical values.
            method (str, optional): Integration method for solve_ivp. Defaults to 'BDF', due to stiffness-related errors with 'RK45'.
            rtol (float, optional): Relative tolerance for the ODE solver. Defaults to 1e-9.
            atol (float, optional): Absolute tolerance for the ODE solver. Defaults to 1e-12.
            
        Returns:
            dict: Trajectory information produced by scipy.integrate.solve_ivp's BDF solver.
        """
        vsq = initial_velocity[0]**2 + initial_velocity[1]**2
        gamma = 1.0 / np.sqrt(1.0 - vsq/self.SPEED_OF_LIGHT**2)
        pz = gamma * self.m * initial_velocity[0]
        pr = gamma * self.m * initial_velocity[1]
>>>>>>> Stashed changes
        
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using the pre-existing
    ElectrodeConfig class in the archetypal arrangement of the unipotential lens.
=======
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using 
    the pre-existing ElectrodeConfig class in the geometry of the unipotential lens.
>>>>>>> Stashed changes
=======
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using 
    the pre-existing ElectrodeConfig class in the geometry of the unipotential lens.
>>>>>>> Stashed changes
=======
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using 
    the pre-existing ElectrodeConfig class in the geometry of the unipotential lens.
>>>>>>> Stashed changes
=======
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using 
    the pre-existing ElectrodeConfig class in the geometry of the unipotential lens.
>>>>>>> Stashed changes
    
    Attributes:
        electrode1 (ElectrodeConfig): Configures the first electrode at 0V.
        electrode2 (ElectrodeConfig): Configures for the second electrode held at (focus_voltage) V.
        electrode3 (ElectrodeConfig): Configuration for the third electrode at 0V.
    """
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    def __init__(self, 
                position: float, 
                width: float, 
                aperture_center: float,
                aperture_width: float,
                outer_diameter: float,
                focus_voltage: float,
                gap_size: int = 1):
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        """
        Creates a parameterizable unipotential/einzel lens with custom geometries, voltages, gap sizes, etc etc.
=======

        """
        Creates a parameterizable einzel lens.
>>>>>>> Stashed changes
=======

        """
        Creates a parameterizable einzel lens.
>>>>>>> Stashed changes
=======

        """
        Creates a parameterizable einzel lens.
>>>>>>> Stashed changes
        
        Parameters:
            position (float): Position where the lens begins on the z-axis in grid units (dz).
            width (float): Width of the full lens assembly on the z-axis in grid units (dz).
            aperture_center (float): Center of the aperture on the r-axis in grid units (dr).
            aperture_width (float): Size of the aperture on the r-axis in grid units (dr).
            outer_diameter (float): Full diameter of the electrodes in grid units on the r-axis (dr).
            focus_voltage (float): Voltage applied to the center electrode in volts.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            gap_size (int, optional): Size of gaps between electrodes in grid units- very important parameter for fringing behavior. Defaults to 1.

        Note:
            If you want to add a focusing lens, make the polarity identical to the charge, and if you want a defocusing lens, make it the opposite. 
            E.g. for electrons a -10000V lens will be converging and a 10000V lens diverging, whereas for positive ions the inverse is true.
            You want unipotential lenses' voltage to be at, above, or near the eV value you give your particles in terms of speed. Note that, since 
            cylindrical electrodes aren't unipotential, you can have accelerating or decelerating behavior from them, but einzel lenses will
            keep your electrons at the same energy.
        """
        electrode_thickness = (width - 3 * gap_size)/3.0 
=======

        """
        Creates a parameterizable einzel lens.
>>>>>>> Stashed changes
        
        Parameters:
            position (float): Position where the lens begins on the z-axis in grid units (dz).
            width (float): Width of the full lens assembly on the z-axis in grid units (dz).
            aperture_center (float): Center of the aperture on the r-axis in grid units (dr).
            aperture_width (float): Size of the aperture on the r-axis in grid units (dr).
            outer_diameter (float): Full diameter of the electrodes in grid units on the r-axis (dr).
            focus_voltage (float): Voltage applied to the center electrode in volts.
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
            gap_size (int, optional): Size of gaps between electrodes in grid units. 
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            field (PotentialField): An empty potential field initialization before adding electrodes to the field calculations.
=======
            field (PotentialField): Initializes a PotentialField so it can add three electrodes to it as defined by __init__().
>>>>>>> Stashed changes
=======
            field (PotentialField): Initializes a PotentialField so it can add three electrodes to it as defined by __init__().
>>>>>>> Stashed changes
=======
            field (PotentialField): Initializes a PotentialField so it can add three electrodes to it as defined by __init__().
>>>>>>> Stashed changes
=======
            field (PotentialField): Initializes a PotentialField so it can add three electrodes to it as defined by __init__().
>>>>>>> Stashed changes
        """
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)

<<<<<<< Updated upstream
class IonOpticsSystem:    
    """
    The top-level class which synthesizes all prior functionality.
    Initializes the potential field, particle tracing, and visualization
    components to provide a complete environment for designing and analyzing
    ion/electron optics systems. 
=======
class IonOpticsSystem:
    """
    Initializes the potential field, particle tracing, and visualization
    components to provide a complete environment for designing and analyzing
    ion/electron optics systems. Contains a lot of referential code for this
    reason.
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    
    Attributes:
        field (PotentialField): The potential field initialized in the simulation.
        tracer (ParticleTracer): The trajectory calcuation tracer.
        elements (list): List of all electrodes and lenses inside the system.
    """
    
    def __init__(self, nr: float, nz: float, axial_size: float = 0.1, radial_size: float = 0.1):
        """
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in meters, and with
        a grid of (nz, nr) units.
=======
        Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in 
        meters, and with a grid of (nz, nr) units.
>>>>>>> Stashed changes
=======
        Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in 
        meters, and with a grid of (nz, nr) units.
>>>>>>> Stashed changes
=======
        Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in 
        meters, and with a grid of (nz, nr) units.
>>>>>>> Stashed changes
=======
        Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in 
        meters, and with a grid of (nz, nr) units.
>>>>>>> Stashed changes
        
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        Adds a single electrode to the ion optics system.
=======
        Adds a single electrode to the ion optics system when ElectrodeConfig is used
        by calling on the pre-existing add_electrode() function.
>>>>>>> Stashed changes
=======
        Adds a single electrode to the ion optics system when ElectrodeConfig is used
        by calling on the pre-existing add_electrode() function.
>>>>>>> Stashed changes
=======
        Adds a single electrode to the ion optics system when ElectrodeConfig is used
        by calling on the pre-existing add_electrode() function.
>>>>>>> Stashed changes
=======
        Adds a single electrode to the ion optics system when ElectrodeConfig is used
        by calling on the pre-existing add_electrode() function.
>>>>>>> Stashed changes
        
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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
=======
        Solves the potential field using the pre-existing
        PyAMG-based multigrid solver.
        
=======
        Solves the potential field using the pre-existing
        PyAMG-based multigrid solver.
        
>>>>>>> Stashed changes
=======
        Solves the potential field using the pre-existing
        PyAMG-based multigrid solver.
        
>>>>>>> Stashed changes
=======
        Solves the potential field using the pre-existing
        PyAMG-based multigrid solver.
        
>>>>>>> Stashed changes
        Returns:
            ndarray: The solved potential field.
        """
        return self.field.solve_potential()
        
    def simulate_beam(self, energy_eV: float, start_z: float,
                              r_range: Tuple[float, float],
                              angle_range: tuple,
                              num_particles: float,
                              simulation_time: float,
                              n_jobs: int = -1):
        """                
        Initializes a parameterized beam, where you can define each particle's
        kinetic energy, angular spread, number of particles. Uses joblib for 
        embarassingly parallel simultaneous trajectory calculations.
        
        Parameters:
            energy_eV (float): Kinetic energy of each particle in electronvolts.
            start_z (float): Starting position for all particles on the z-axis in meters, not in grid units.
            r_range (Tuple[float, float]): The range of initial radial positions in meters- effectively the beam width.
            angle_range (tuple): Range of initial angles from the horizontal in radians.
            num_particles (float): Number of particles to simulate in the beam- keep at 3-6 for prototyping and 10-100 for full simulation.
            simulation_time (float): Total simulation time in seconds, typically between 1e-10 to 1e-7.
            n_jobs (float, optional): Number of parallel jobs (cores) used by joblib. Defaults to all available.

        Returns:
            list: List of trajectory solutions for all simulated particles.            
        """
        velocity_magnitude = self.tracer.get_velocity_from_energy(energy_eV)
        min_angle_rad = np.radians(angle_range[0])
        max_angle_rad = np.radians(angle_range[1])
        angles = np.linspace(min_angle_rad, max_angle_rad, int(num_particles))
        r_positions = np.linspace(r_range[0], r_range[1], int(num_particles))
        particle_params = []
        for r_pos, angle in zip(r_positions, angles):
            vz = velocity_magnitude * np.cos(angle)
            vr = velocity_magnitude * np.sin(angle)
            particle_params.append((start_z, r_pos, vz, vr))
        
        def trace_particle(params):
            """
            Calculates a single particle trajectory in a way amenable to parallel processing.
            Called by simulate_beam for parallel processing of beam trajectories.

            Parameters:
                params(tuple): (z0, r0, vz0, vr0) are the initial positions and velocities.

            Returns:
                trajectories (scipy.integrate.OdeResult): Contains trajectory data at time points and 
                state vectors at each timepoint [z, r, pz, pr]
            """
            z0, r0, vz0, vr0 = params
            return self.tracer.trace_trajectory(
                initial_position=(z0, r0),
                initial_velocity=(vz0, vr0),
                simulation_time=simulation_time
            )
        trajectories = Parallel(n_jobs=n_jobs)(
            delayed(trace_particle)(params) for params in particle_params
        )
        return trajectories

    def visualize_system(self, trajectories=None, r_limits=None, figsize=(16, 10), title="Picht"):
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
        """
        fig = plt.figure(figsize=figsize)
        ax_main = fig.add_axes([0.1, 0.1, 0.7, 0.8])
        ax_main.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax_main.set_xlabel('Axial position (m)', fontsize=12)
        ax_main.set_ylabel('Radial position (m)', fontsize=12)
        ax_main.grid(True, alpha=0.3)
        ax_checkbox = fig.add_axes([0.82, 0.4, 0.15, 0.2])
        ax_checkbox.set_facecolor('white')
        ax_checkbox.set_title('Display Options', fontsize=12, fontweight='bold', pad=10)
        checkbox = CheckButtons(ax_checkbox, 
                               ['Lenses', 'Electric Field', 'Animate'],
                               [True, False, False])
        for text in checkbox.labels:
            text.set_fontsize(11)
        self.animation = None
        self.trajectory_lines = []
        self.lens_patches = []
        def draw_lenses():
            self.lens_patches.clear()
            if checkbox.get_status()[0]:  
                for element in self.elements:
                    if isinstance(element, ElectrodeConfig):
                        z_start = element.start * self.field.dz
                        z_end = (element.start + element.width) * self.field.dz
                        ap_center = element.ap_start + element.ap_width / 2
                        r_inner = element.ap_start * self.field.dr
                        r_outer = (element.ap_start + element.ap_width) * self.field.dr
                        inner_electrode = max(0, (ap_center - element.outer_diameter / 2) * self.field.dr)
                        outer_electrode = min((ap_center + element.outer_diameter / 2) * self.field.dr, self.field.radial_size)
                        patch1 = ax_main.fill_between([z_start, z_end], [inner_electrode, inner_electrode], 
                                                     [r_inner, r_inner], color='#000000', alpha=1, linewidth=0)
                        patch2 = ax_main.fill_between([z_start, z_end], [r_outer, r_outer], 
                                                     [outer_electrode, outer_electrode], color='#000000', alpha=1, linewidth=0)
                        self.lens_patches.extend([patch1, patch2])
                    elif isinstance(element, EinzelLens):
                        for electrode in [element.electrode1, element.electrode2, element.electrode3]:
                            z_start = electrode.start * self.field.dz
                            z_end = (electrode.start + electrode.width) * self.field.dz                           
                            ap_center = electrode.ap_start + electrode.ap_width / 2
                            r_inner = electrode.ap_start * self.field.dr
                            r_outer = (electrode.ap_start + electrode.ap_width) * self.field.dr                           
                            inner_electrode = max(0, (ap_center - electrode.outer_diameter / 2) * self.field.dr)
                            outer_electrode = min((ap_center + electrode.outer_diameter / 2) * self.field.dr, self.field.radial_size)                           
                            patch1 = ax_main.fill_between([z_start, z_end], [inner_electrode, inner_electrode], 
                                                         [r_inner, r_inner], color='#000000', alpha=1, linewidth=0)
                            patch2 = ax_main.fill_between([z_start, z_end], [r_outer, r_outer], 
                                                         [outer_electrode, outer_electrode], color='#000000', alpha=1, linewidth=0)
                            self.lens_patches.extend([patch1, patch2])
        
        def draw_trajectories():
            """
            Draws/redraws static trajectories for particles when the animated tickbox is left unchecked.
            """
            self.trajectory_lines.clear()           
            if trajectories and not checkbox.get_status()[2]:  
                colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
                for i, sol in enumerate(trajectories):
                    z_traj = sol.y[0]
                    r_traj = sol.y[1]
                    line, = ax_main.plot(z_traj, r_traj, lw=2, color=colors[i])
                    self.trajectory_lines.append(line)
        
        def setup_animation():
            """
            Handles particle animation logic- animates trajectories when the animate tickbox is checked
            and reverts to static when it's unchecked.
            """
            if self.animation is not None:
                self.animation.event_source.stop()
                self.animation = None
            self.trajectory_lines.clear()
            if checkbox.get_status()[2] and trajectories:  
                colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
                for i in range(len(trajectories)):
                    line, = ax_main.plot([], [], lw=2, color=colors[i])
                    self.trajectory_lines.append(line)
                
                def init():
                    """
                    Initializes line objects for matplotlib's FuncAnimation to work correctly.

                    Returns:
                        list: Empty trajectory line objects
                    """
                    for line in self.trajectory_lines:
                        line.set_data([], [])
                    return self.trajectory_lines
                
                def animate(frame):
                    """
                    Animates particle positions based on the frame it's currently at by 
                    checking what frame the particle is at, and using trajectory data to
                    extend the path a bit more.
    
    
                    Parameters:
                        frame (int): Current frame number, ranges from 0 to 99, loops after 99 to 0
        
                    Returns:
                        list: Updated line objects for visualization
                    """
                    for i, (line, sol) in enumerate(zip(self.trajectory_lines, trajectories)):
                        max_index = len(sol.t)
                        current_index = int((frame / 100.0) * max_index) % max_index 
                        z_data = sol.y[0][:current_index]
                        r_data = sol.y[1][:current_index]
                        line.set_data(z_data, r_data)
                    
                    return self.trajectory_lines
                
                self.animation = FuncAnimation(fig, animate, init_func=init, frames=100, 
                                             interval=int(20/0.3), blit=True, repeat=True)
            else:
                draw_trajectories()

        def on_check(label):
            """
            Handles the checkbox logic from UI interactions. Identifies whether its
            state is checked or unchecked and handles the visualization changes based
            on state.

            Parameters:
                label (str): Name of the checkbox that was clicked
            """
            if self.animation is not None:
                self.animation.event_source.stop()
                self.animation = None
            ax_main.clear()
            ax_main.set_title(title, fontsize=16, fontweight='bold', pad=20)
            ax_main.set_xlabel('z position (meters)', fontsize=12)
            ax_main.set_ylabel('r position (meters)', fontsize=12)
            ax_main.grid(True, alpha=0.3)
            if r_limits:
                ax_main.set_ylim(r_limits)
            else:
                ax_main.set_ylim(0, self.field.radial_size)
            ax_main.set_xlim(0, self.field.axial_size)
            if checkbox.get_status()[1]:
                z_coords = np.linspace(0, self.field.axial_size, 200)
                r_coords = np.linspace(0, self.field.radial_size, 100)
                Z, R = np.meshgrid(z_coords, r_coords)
                E_magnitude = np.zeros_like(Z)
                for i in range(Z.shape[0]):
                    for j in range(Z.shape[1]):
                        Ez, Er = self.field.get_field_at_position(Z[i,j], R[i,j])
                        E_magnitude[i,j] = np.sqrt(Ez**2 + Er**2)
                ax_main.contourf(Z, R, E_magnitude, levels=50, cmap='turbo', alpha=1)
            draw_lenses()    
            if checkbox.get_status()[2]:
                setup_animation()
            else:
                draw_trajectories()
            fig.canvas.draw_idle()
        checkbox.on_clicked(on_check)
        draw_lenses()
        draw_trajectories()
        if r_limits:
            ax_main.set_ylim(r_limits)
        else:
            ax_main.set_ylim(0, self.field.radial_size)
        ax_main.set_xlim(0, self.field.axial_size)
        
        fig.patch.set_facecolor('white')
        
        return fig
class Export:
    """
    Exports trajectory data to the HDF5 format for data analysis, and electrode
    data to STEP for integration with Parametric CAD software and manufacturing.
    Entirely separate from other classes, ensures absolute backwards compatibility
    and minimal interference.
    
    Attributes:
        system: Reference to the ion optics system
        field: Shortcut to the system's potential field
    """
    def __init__(self, ion_optics_system):
        """
        Creates the export handler by initializing/pulling from the IonOpticsSystem.
        
        Parameters:
            ion_optics_system: IonOpticsSystem instance to export data from
        """
        self.system = ion_optics_system
        self.field = ion_optics_system.field
   
    def export_traj(self, trajectories):
        """
        Exports simulation data to HDF5 using the H5PY library, including electrode
        configuration data, particle trajectory data, and electric field data.
        
        Parameters:
            trajectories: List of scipy.integrate.OdeResult objects from simulate_beam
            
        Output:
            simulation_results.h5, which contains:
            - Electric field data (electric potential, Ez, Er)
            - Particle trajectories (pz and pr, momentum vs time)
            - System configuration (grid information, electrodes, ion/electron properties)
        """
        import h5py
        
        with h5py.File("simulation_results.h5", 'w') as f:
            field_group = f.create_group('fields')
            field_group.create_dataset('potential', data=self.field.potential)
            field_group.create_dataset('Ez', data=self.field.Ez)
            field_group.create_dataset('Er', data=self.field.Er)
            grid_group = f.create_group('grid')
            grid_group.attrs['dz'] = self.field.dz
            grid_group.attrs['dr'] = self.field.dr
            grid_group.attrs['nz'] = self.field.nz
            grid_group.attrs['nr'] = self.field.nr
            grid_group.attrs['axial_size'] = self.field.axial_size
            grid_group.attrs['radial_size'] = self.field.radial_size
            traj_group = f.create_group('trajectories')
            for i, sol in enumerate(trajectories):
                particle = traj_group.create_group(f'particle_{i}')
                particle.create_dataset('time', data=sol.t)
                particle.create_dataset('z', data=sol.y[0])
                particle.create_dataset('r', data=sol.y[1])
                particle.create_dataset('pz', data=sol.y[2])
                particle.create_dataset('pr', data=sol.y[3])
            ion_group = f.create_group('ion')
            ion_group.attrs['mass'] = self.system.tracer.m
            ion_group.attrs['charge'] = self.system.tracer.q
            ion_group.attrs['symbol'] = self.system.tracer.current_ion['symbol']
            electrode_group = f.create_group('electrodes')
            for i, element in enumerate(self.system.elements):
                if isinstance(element, ElectrodeConfig):
                    elec = electrode_group.create_group(f'electrode_{i}')
                    elec.attrs['start'] = element.start
                    elec.attrs['width'] = element.width
                    elec.attrs['ap_start'] = element.ap_start
                    elec.attrs['ap_width'] = element.ap_width
                    elec.attrs['outer_diameter'] = element.outer_diameter
                    elec.attrs['voltage'] = element.voltage
                elif isinstance(element, EinzelLens):
                    lens = electrode_group.create_group(f'einzel_lens_{i}')
                    for j, electrode in enumerate([element.electrode1, element.electrode2, element.electrode3]):
                        sub_elec = lens.create_group(f'electrode_{j}')
                        sub_elec.attrs['start'] = electrode.start
                        sub_elec.attrs['width'] = electrode.width
                        sub_elec.attrs['ap_start'] = electrode.ap_start
                        sub_elec.attrs['ap_width'] = electrode.ap_width
                        sub_elec.attrs['outer_diameter'] = electrode.outer_diameter
                        sub_elec.attrs['voltage'] = electrode.voltage
   
    def cad_export(self):
        """
        Exports electrode geometry to the .step file format, by converting the
        full IonOpticsSystem geometry into actual CAD files with cadquery.
        
        Output:
            save.step file, saved in the same directory as the script run.
        """
        import cadquery as cq
        import os
        import sys
        
        shapes = []
        for element in self.system.elements:
            if isinstance(element, ElectrodeConfig):
                electrode_shape = self._create_electrode_shape(element)
                shapes.append(electrode_shape)
                
            elif isinstance(element, EinzelLens):
                for electrode in [element.electrode1, element.electrode2, element.electrode3]:
                    electrode_shape = self._create_electrode_shape(electrode)
                    shapes.append(electrode_shape)
       
        if shapes:
            combined = shapes[0]
            for shape in shapes[1:]:
                combined = combined.union(shape)
           
            output_path = os.path.join(os.getcwd(), "save.step")
            cq.exporters.export(combined, output_path)
            print(f"Exported to: {output_path}")
           
            sys.exit(0)
   
    def _create_electrode_shape(self, electrode_config):
        """
        Converts 2D electrode config data to 3D CAD geometric data.
        
        Parameters:
            electrode_config: ElectrodeConfig data in 2D axisymmetric data
            
        Returns:
            3D CAD solid representation of the electrode
        """
        import cadquery as cq
        
        z_start_mm = electrode_config.start * self.field.dz * 1000
        z_length_mm = electrode_config.width * self.field.dz * 1000
       
        outer_radius_mm = (electrode_config.outer_diameter * self.field.dr * 1000) / 2
        inner_radius_mm = (electrode_config.ap_width * self.field.dr * 1000) / 2
       
        ap_center_mm = (electrode_config.ap_start + electrode_config.ap_width / 2) * self.field.dr * 1000
       
        if inner_radius_mm > 0 and inner_radius_mm < outer_radius_mm:
            electrode = (cq.Workplane("XY")
                        .workplane(offset=z_start_mm)
                        .center(0, ap_center_mm)
                        .circle(outer_radius_mm)
                        .circle(inner_radius_mm)
                        .extrude(z_length_mm))
        else:
            electrode = (cq.Workplane("XY")
                        .workplane(offset=z_start_mm)
                        .center(0, ap_center_mm)
                        .circle(outer_radius_mm)
                        .extrude(z_length_mm))
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        return electrode
>>>>>>> Stashed changes
=======
        return electrode
>>>>>>> Stashed changes
=======
        return electrode
>>>>>>> Stashed changes
=======
        return electrode
>>>>>>> Stashed changes
