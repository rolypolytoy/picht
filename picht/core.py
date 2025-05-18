import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import CheckButtons, Slider
from scipy.integrate import solve_ivp
from scipy.sparse import diags, csr_matrix, lil_matrix
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable
import numba as nb
from mendeleev import element
import pyamg
from joblib import Parallel, delayed
import multiprocessing
import os

@nb.njit
def get_field(z, r, Ez, Er, axial_size, radial_size, dz, dr, nz, nr):
    """
    Provides electric field values at fractional grid positions by
    picking the value at the nearest neighbor. Necessary because 
    FDM is discretizing by nature.
    
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
        i = int(min(max(0, z / dz), nz - 1))
        j = int(min(max(0, r / dr), nr - 1))
        return Ez[i, j], Er[i, j]
    else:
        return 0.0, 0.0


@nb.njit
def calc_dynamics(z, r, pz, pr, Ez, Er, Bz, Br, q, m, c, r_axis=0.0):
    """    
    Calculates the acceleration of charged particles by applying
    the Lorentz force with special-relativistic corrections. Uses energy
    momentum formalism for full eV to TeV support.

    Note: F = q^2 Bz^2 /(4m) * r is the paraxial ray equation for small-angle magnetostatics.
    This is near-necessary if you're simulating things axisymmetrically because you can't simulate 
    azimuthal velocity, so the paraxial ray equation is necessary for keeping both axisymmetry and
    accuracy. Since Fr is equivalent to dpr/dt we basically solve for this.
    
    Parameters:
        z (float): Z-axis position in meters.
        r (float): R-axis position in meters.
        pz (float): Z-axis momentum in kg * m/s.
        pr (float): R-axis momentum in kg * m/s.
        Ez (float): Z-axis electric field in volts per meter.
        Er (float): R-axis electric field in volts per meter.
        Bz (float): Z-axis magnetic field in Tesla.
        Br (float): R-axis electric field in Tesla.
        q (float): Charge of the particle.
        mass (float): Particle mass in kg.
        c (float): Speed of light in a vacuum in meters per second.
        
    Returns:
        ndarray: Array containing [vz, vr, dpz_dt, dpr_dt], representing velocity and force
                components in the z and r directions for complete kinematic information (force is dp/dt).
    """
    p_sq = pz**2 + pr**2
    E = np.sqrt((p_sq * c**2) + (m * c**2)**2)
    vz = pz * c**2 / E
    vr = pr * c**2 / E    
    dpz_dt = q * Ez
    dpr_dt = q * Er    
    r_from_axis = r - r_axis
    dpr_dt += -((q**2) * Bz**2 / (4 * m)) * r_from_axis
    return np.array([vz, vr, dpz_dt, dpr_dt])

@dataclass
class MagneticLensConfig:
    """
    Parameterizes a simple magnetic lens.
    
    Attributes:
        start (float): Position where the magnetic lens begins on the z-axis in grid units.
        length (float): Length of the magnetic lens on the z-axis in grid units.
        ap_start (float): Starting position of the aperture on the r-axis in grid units.
        ap_width (float): Width of the aperture in the lens on the r-axis in grid units.
        outer_diameter (float): Diameter of the magnetic lens on the r-axis in grid units.
        mu_r (float): Relative magnetic permeability of the lens material in dimensionless units.
        mmf (float): Magnetomotive force of the lens in ampere-turns.
    """
    start: float
    length: float
    ap_start: float
    ap_width: float
    outer_diameter: float
    mu_r: float
    mmf: float

class MagneticField:
   """
   This class handles everything related to the magnetic vector potential field, including
   adding magnetic lenses, and solving for the magnetic field. Analagous to ElectricField
   but solves the magnetic field instead of the electric field.
   
   Attributes:
       nz (int): Number of grid points in the z-axis.
       nr (int): Number of grid points in the r-axis.
       axial_size (float): Size of the z-axis in meters.
       radial_size (float): Size of the r-axis in meters.
       dz (float): Grid spacing in z-axis per z-unit, in meters, found by axial_size/nz.
       dr (float): Grid spacing in r-axis per r-unit, in meters, found by radial_size/nr.
       vector_potential (ndarray): 2D array with magnetic vector potential values in Tesla-meters.
       magnetic_mask (ndarray): Boolean mask for magnetic material positions.
       mu_r (ndarray): 2D array of relative permeability values.
       current_density (ndarray): 2D array of current density values in amperes per square meter.
       Bz (ndarray): Z-component of magnetic field in Tesla.
       Br (ndarray): R-component of magnetic field in Tesla.
       lens_config (MagneticLensConfig): Configuration of the magnetic lens.
   """
  
   def __init__(self, ion_optics_system):
       """
       Initializes the magnetic field with the specific dimensions. Generally when
       calculating magnetic lenses, you need to have finer meshes than with solely
       electrostatic lenses, because there are more discretization artefacts
       because of the paraxial ray equation. Since we use multigrid methods 
       this isn't exceptionally costly.
       
       Parameters:
           ion_optics_system (IonOpticsSystem): The parent ion optics system containing
               field dimensions and grid parameters.
       """
       self.nz = ion_optics_system.field.nz
       self.nr = ion_optics_system.field.nr
       self.axial_size = ion_optics_system.field.axial_size
       self.radial_size = ion_optics_system.field.radial_size
       self.dz = ion_optics_system.field.dz
       self.dr = ion_optics_system.field.dr
       self.vector_potential = np.zeros((self.nz, self.nr))
       self.magnetic_mask = np.zeros((self.nz, self.nr), dtype=bool)
       self.mu_r = np.ones((self.nz, self.nr))
       self.current_density = np.zeros((self.nz, self.nr))
       self.Bz = np.zeros((self.nz, self.nr))
       self.Br = np.zeros((self.nz, self.nr))
       self.lens_config = None
      
   def add_magnetic_lens(self, config):
       """
       Adds a magnetic lens to the field and handles all necessary calculations.
       Analagous to add_electrode() in ElectricField.

       Parameters:
           config (MagneticLensConfig): Configuration parameters for the magnetic lens.
       """
       self.lens_config = config
       start = int(config.start)
       end = int(config.start + config.length)
       ap_center = config.ap_start + config.ap_width / 2
       bore_radius = int(config.ap_width / 2)
       outer_radius = int(config.outer_diameter / 2)
       
       area = 0
       for i in range(start, end):
           for j in range(self.nr):
               r_from_axis = abs(j - ap_center)
               if bore_radius <= r_from_axis <= outer_radius:
                   area += 1
       
       area_physical = area * self.dz * self.dr
       current_density = config.mmf / area_physical if area_physical > 0 else 0
       
       for i in range(start, end):
           for j in range(self.nr):
               r_from_axis = abs(j - ap_center)
               if bore_radius <= r_from_axis <= outer_radius:
                   self.mu_r[i, j] = config.mu_r
                   self.current_density[i, j] = current_density
                   self.magnetic_mask[i, j] = True
  
   def build_laplacian_matrix(self, mask, dirichlet_values=None):
       """
       Builds a sparse matrix for the Laplacian ∇²A = -μ₀μᵣJ, 
       and implements Neumann boundary conditions at all boundaries, to simulate
       how magnetic fields behave at metal boundaries. Identical implementation
       details as build_laplacian_matrix() in ElectricField.

       Parameters:
           mask (ndarray): Boolean array, true where magnetic materials exist, false elsewhere.
           dirichlet_values (ndarray, optional): Vector potential values where mask is True.

       Returns:
           A (scipy.sparse.csr_matrix): A sparse matrix representation of the Laplacian operator
           on vector potential. Its shape is (nz * nr, nz * nr).
           b (ndarray): Contains the source terms from current density and boundary conditions.
       """
       n = self.nz * self.nr
       A = lil_matrix((n, n))
       b = np.zeros(n)
       mu_0 = 4*np.pi*1e-7

       def idx(i, j):
           return i * self.nr + j

       def harmonic_mean(a, b):
           if a == 0 or b == 0:
               return 0
           return 2 * a * b / (a + b)

       ap_center = (self.lens_config.ap_start + self.lens_config.ap_width / 2) if self.lens_config else 0
       epsilon = 1e-10

       for i in range(self.nz):
           for j in range(self.nr):
               k = idx(i, j)
               
               if i < 2 or i >= self.nz - 2:
                   A[k, k] = -1.0
                   if i == 0:
                       A[k, idx(min(i+1, self.nz-1), j)] = 1.0
                   elif i == 1:
                       A[k, idx(min(i+1, self.nz-1), j)] = 1.0
                   elif i == self.nz - 2:
                       A[k, idx(max(i-1, 0), j)] = 1.0
                   else:
                       A[k, idx(max(i-1, 0), j)] = 1.0
                   b[k] = 0.0
               elif j < 2 or j >= self.nr - 2:
                   A[k, k] = -1.0
                   if j == 0:
                       A[k, idx(i, min(j+1, self.nr-1))] = 1.0
                   elif j == 1:
                       A[k, idx(i, min(j+1, self.nr-1))] = 1.0
                   elif j == self.nr - 2:
                       A[k, idx(i, max(j-1, 0))] = 1.0
                   else:
                       A[k, idx(i, max(j-1, 0))] = 1.0
                   b[k] = 0.0
               else:
                   mu_center = self.mu_r[i, j]
                   mu_left2 = self.mu_r[i-2, j]
                   mu_left1 = self.mu_r[i-1, j]
                   mu_right1 = self.mu_r[i+1, j]
                   mu_right2 = self.mu_r[i+2, j]
                   mu_down2 = self.mu_r[i, j-2]
                   mu_down1 = self.mu_r[i, j-1]
                   mu_up1 = self.mu_r[i, j+1]
                   mu_up2 = self.mu_r[i, j+2]
                   
                   mu_interface_left2 = harmonic_mean(mu_center, mu_left2)
                   mu_interface_left1 = harmonic_mean(mu_center, mu_left1)
                   mu_interface_right1 = harmonic_mean(mu_center, mu_right1)
                   mu_interface_right2 = harmonic_mean(mu_center, mu_right2)
                   mu_interface_down2 = harmonic_mean(mu_center, mu_down2)
                   mu_interface_down1 = harmonic_mean(mu_center, mu_down1)
                   mu_interface_up1 = harmonic_mean(mu_center, mu_up1)
                   mu_interface_up2 = harmonic_mean(mu_center, mu_up2)
                   
                   r_from_axis = (j - ap_center) * self.dr
                   r_abs = max(abs(r_from_axis), epsilon)
                   
                   if r_from_axis > 0:
                       r_deriv_coeff_down2 = -1/(12*r_abs*self.dr)
                       r_deriv_coeff_down1 = 8/(12*r_abs*self.dr)
                       r_deriv_coeff_up1 = -8/(12*r_abs*self.dr)
                       r_deriv_coeff_up2 = 1/(12*r_abs*self.dr)
                   else:
                       r_deriv_coeff_down2 = 1/(12*r_abs*self.dr)
                       r_deriv_coeff_down1 = -8/(12*r_abs*self.dr)
                       r_deriv_coeff_up1 = 8/(12*r_abs*self.dr)
                       r_deriv_coeff_up2 = -1/(12*r_abs*self.dr)
                   
                   coeff_left2 = -mu_interface_left2 / (12*self.dz**2)
                   coeff_left1 = 16*mu_interface_left1 / (12*self.dz**2)
                   coeff_right1 = 16*mu_interface_right1 / (12*self.dz**2)
                   coeff_right2 = -mu_interface_right2 / (12*self.dz**2)
                   
                   coeff_down2 = -mu_interface_down2 / (12*self.dr**2) + mu_interface_down2 * r_deriv_coeff_down2
                   coeff_down1 = 16*mu_interface_down1 / (12*self.dr**2) + mu_interface_down1 * r_deriv_coeff_down1
                   coeff_up1 = 16*mu_interface_up1 / (12*self.dr**2) + mu_interface_up1 * r_deriv_coeff_up1
                   coeff_up2 = -mu_interface_up2 / (12*self.dr**2) + mu_interface_up2 * r_deriv_coeff_up2
                   
                   coeff_center = -(coeff_left2 + coeff_left1 + coeff_right1 + coeff_right2 + 
                                  coeff_down2 + coeff_down1 + coeff_up1 + coeff_up2 + 
                                  mu_center/r_abs**2)
                   
                   A[k, k] = coeff_center
                   A[k, idx(i-2, j)] = coeff_left2
                   A[k, idx(i-1, j)] = coeff_left1
                   A[k, idx(i+1, j)] = coeff_right1
                   A[k, idx(i+2, j)] = coeff_right2
                   A[k, idx(i, j-2)] = coeff_down2
                   A[k, idx(i, j-1)] = coeff_down1
                   A[k, idx(i, j+1)] = coeff_up1
                   A[k, idx(i, j+2)] = coeff_up2
                   
                   if self.current_density[i, j] != 0:
                       b[k] = -mu_0 * mu_center * self.current_density[i, j]
                   else:
                       b[k] = 0.0
                      
       return A.tocsr(), b

   def solve_vector_potential(self, max_iterations: float = 500, convergence_threshold: float = 1e-6):
       """
       Solves the magnetic vector potential field using Multigrid methods with PyAMG. It first 
       creates a CSR matrix for the Laplacian operator for the vector potential field, and uses
       PyAMG to solve for ∇²A = -μ₀μᵣJ. Then, it calculates B = ∇ × A to find the magnetic field.
       This is functionally identical to how solve_potential() in ElectricField calculates the 
       electric field.
       
       Parameters:
           max_iterations (float, optional): Maximum number of iterations for the solver. Defaults to 500.
           convergence_threshold (float, optional): Convergence criterion for the solution. Defaults to 1e-6.
           
       Returns:
           ndarray: The solved vector potential field.
       """
       A, b = self.build_laplacian_matrix(self.magnetic_mask, self.vector_potential)
       
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
       self.vector_potential = x.reshape((self.nz, self.nr))
       
       self.calculate_b_from_a()
       
       return self.vector_potential
  
   def calculate_b_from_a(self):
       """
       Calculates the magnetic field components from the vector potential using B = ∇ × A, with special
       handling of differentiation at boundaries and at r = 0.        
       """
       ap_center = (self.lens_config.ap_start + self.lens_config.ap_width / 2) if self.lens_config else 0
       epsilon = 1e-10
       
       for i in range(self.nz):
           for j in range(self.nr):
               if 2 <= i <= self.nz - 3:
                   self.Br[i, j] = -((-self.vector_potential[i+2, j] + 8*self.vector_potential[i+1, j] - 
                                    8*self.vector_potential[i-1, j] + self.vector_potential[i-2, j]) / (12*self.dz))
               elif i == 0:
                   self.Br[i, j] = -(self.vector_potential[i+1, j] - self.vector_potential[i, j]) / self.dz
               elif i == 1:
                   self.Br[i, j] = -(self.vector_potential[i+1, j] - self.vector_potential[i-1, j]) / (2*self.dz)
               elif i == self.nz - 2:
                   self.Br[i, j] = -(self.vector_potential[i+1, j] - self.vector_potential[i-1, j]) / (2*self.dz)
               else:
                   self.Br[i, j] = -(self.vector_potential[i, j] - self.vector_potential[i-1, j]) / self.dz
               
               r_dist = max(abs((j - ap_center) * self.dr), epsilon)
               r_signed = (j - ap_center) * self.dr
               
               if 2 <= j <= self.nr - 3:
                   dA_dr = (-self.vector_potential[i, j+2] + 8*self.vector_potential[i, j+1] - 
                           8*self.vector_potential[i, j-1] + self.vector_potential[i, j-2]) / (12*self.dr)
               elif j == 0:
                   if j+2 < self.nr:
                       dA_dr = (-3*self.vector_potential[i, j] + 4*self.vector_potential[i, j+1] - 
                               self.vector_potential[i, j+2]) / (2*self.dr)
                   else:
                       dA_dr = (self.vector_potential[i, j+1] - self.vector_potential[i, j]) / self.dr
               elif j == 1:
                   dA_dr = (self.vector_potential[i, j+1] - self.vector_potential[i, j-1]) / (2*self.dr)
               elif j == self.nr - 2:
                   dA_dr = (self.vector_potential[i, j+1] - self.vector_potential[i, j-1]) / (2*self.dr)
               else:
                   if j-2 >= 0:
                       dA_dr = (3*self.vector_potential[i, j] - 4*self.vector_potential[i, j-1] + 
                               self.vector_potential[i, j-2]) / (2*self.dr)
                   else:
                       dA_dr = (self.vector_potential[i, j] - self.vector_potential[i, j-1]) / self.dr
               
               if r_signed < 0:
                   dA_dr = -dA_dr
               
               self.Bz[i, j] = dA_dr + self.vector_potential[i, j] / r_dist
  
   def get_field_at_position(self, z: float, r: float):
       """
       Returns the magnetic field components at a specific position.
       
       Parameters:
           z (float): Position along the z-axis in meters.
           r (float): Position along the r-axis in meters.
           
       Returns:
           Tuple[float, float]: The magnetic field components (Bz, Br) at the specified position.
       """
       if 0 <= z < self.axial_size and 0 <= r < self.radial_size:
           i = int(min(max(0, z / self.dz), self.nz - 1))
           j = int(min(max(0, r / self.dr), self.nr - 1))
           return self.Bz[i, j], self.Br[i, j]
       else:
           return 0.0, 0.0
                                          
@dataclass
class ElectrodeConfig:
    """
    Parameterizes a simple electrostatic lens.
    
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

class ElectricField:
    """
    This class handles everything related to the electric potential field, including
    adding electrodes and einzel lenses, and solving for the electric field.
    
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
        Initializes the potential field with the specific dimensions. This decision is 
        extremely important, as finer meshes increase accuracy at the cost of performance.
        
        Parameters:
            nz (float): Number of grid points in the z-axis.
            nr (float): Number of grid points in the r-axis.
            axial_size (float): Physical length of the z-axis in meters.
            radial_size (float): Physical length of the r-axis in meters.
        """
        self.nz = int(nz)
        self.nr = int(nr)
        self.axial_size = axial_size
        self.radial_size = radial_size
        self.dz = axial_size / nz
        self.dr = radial_size / nr
        self.potential = np.zeros((self.nz, self.nr))
        self.electrode_mask = np.zeros((self.nz, self.nr), dtype=bool)
        self.Ez = np.zeros((self.nz, self.nr))
        self.Er = np.zeros((self.nz, self.nr))
    
    def add_electrode(self, config: ElectrodeConfig):
        """
        Adds a single electrode to the electric field and handles all
        necessary calculations.
                
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
        than a 100x100 one. It first creates a CSR matrix for the Laplacian values for the potential field, 
        uses PyAMG to actually solve for Laplace's equation ∇²V = 0. Then, it finds the gradient E = -∇V and 
        thus finds the electric field.
        
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
        """
        return get_field(z, r, self.Ez, self.Er, self.axial_size, self.radial_size, 
                         self.dz, self.dr, self.nz, self.nr)

class ParticleTracer:
    """
    Handles trajectory dynamics calculations and visualizations.
    
    Attributes:
        field (ElectricField): The electric field class, which has most of the information required
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
    SPEED_OF_LIGHT = 299792458.0

    def __init__(self, electric_field: ElectricField):
        """
        Initializes an electric field with the class.
        
        Parameters:
            electric_field (ElectricField): The electric field required to calculate particle dynamics.
        """
        self.field = electric_field
        self.current_ion = {
            'symbol': 'e-',
            'atomic_number': 0,
            'mass': self.ELECTRON_MASS,
            'charge': self.ELECTRON_CHARGE,
            'charge_mass_ratio': self.ELECTRON_CHARGE / self.ELECTRON_MASS
        }
        self.q = self.current_ion['charge']
        self.m = self.current_ion['mass']

    def set_ion(self, symbol: str = 'e-', charge_state: float = 1):
        """        
        Configures the particle tracer to pick a specific charged particle. 
        Integrates with the Mendeleev library to automatically retrieve
        information about any and all atoms, and natively supports electrons.

        Parameters:
            symbol (str, optional): Chemical symbol of the element, or 'e-' for electrons. Defaults to 'e-' for backwards compatibility.
            charge_state (float, optional): Charge of the ion, defaults to 1.
            
        Returns:
            ParticleTracer: Uses self-reference for method chaining to support Picht's general style of power but conciseness.
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
            electron_charge = 1.60217663e-19
            ion_charge = charge_state * electron_charge
            
            self.current_ion = {
                'symbol': f"{symbol}{'+' if charge_state > 0 else '-'}{abs(charge_state)}",
                'atomic_number': elem.atomic_number,
                'mass': isotope_mass * 1.66053906660e-27,
                'charge': ion_charge,
                'charge_mass_ratio': ion_charge / (isotope_mass * 1.66053906660e-27)
            }
        
        self.q = self.current_ion['charge']
        self.m = self.current_ion['mass']
        return self

    def get_velocity_from_energy(self, energy_eV: float) -> float:
        """
        Converts particle energy in electronvolts to velocity in meters per second,
        accounting for relativistic effects. It's accurate for all energy scales from single-digit eV to GeV.
        
        Parameters:
            energy_eV (float): Kinetic energy of the particle in electronvolts, the standard unit of energy in
            particle physics. 1eV is approximately 1.6022e-19 Joules, and is the kinetic energy an electron
            has after being accelerated through an electric field with potential difference of 1 Volt.
            
        Returns:
            float: Particle velocity in meters per second.
        """
        kinetic_energy = energy_eV * 1.60217663e-19
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
            state (List[float]): Current state [z, r, pz, pr] with position and momentum components.
            
        Returns:
            List[float]: Derivatives of the state vector [vz, vr, dpz_dt, dpr_dt] representing velocities and force components.
        """
        z, r, pz, pr = state
        Ez, Er = self.field.get_field_at_position(z, r)
        Bz, Br = 0.0, 0.0
        r_axis = 0.0
        if hasattr(self, 'magnetic_lenses') and self.magnetic_lenses is not None:
            Bz, Br = self.magnetic_lenses.get_field_at_position(z, r)
            r_axis = (self.magnetic_lenses.lens_config.ap_start + 
                self.magnetic_lenses.lens_config.ap_width / 2) * self.field.dr

    
        return calc_dynamics(
            z, r, pz, pr, Ez, Er, Bz, Br,
            self.q, self.m, self.SPEED_OF_LIGHT, r_axis)
    
    def trace_trajectory(self, 
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
        
        initial_state = [
            initial_position[0], 
            initial_position[1],
            pz,
            pr
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
    Implements an Einzel (unipotential) lens for charged particle focusing. Implements three electrodes using 
    the pre-existing ElectrodeConfig class in the geometry of the unipotential lens.
    
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
        Creates a parameterizable einzel lens.
        
        Parameters:
            position (float): Position where the lens begins on the z-axis in grid units (dz).
            width (float): Width of the full lens assembly on the z-axis in grid units (dz).
            aperture_center (float): Center of the aperture on the r-axis in grid units (dr).
            aperture_width (float): Size of the aperture on the r-axis in grid units (dr).
            outer_diameter (float): Full diameter of the electrodes in grid units on the r-axis (dr).
            focus_voltage (float): Voltage applied to the center electrode in volts.
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
    
    def add_to_field(self, field: ElectricField):
        """        
        Parameters:
            field (ElectricField): Initializes an ElectricField so it can add three electrodes to it as defined by __init__().
        """
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)

class IonOpticsSystem:
  """
  Initializes the potential field, particle tracing, and visualization
  components to provide a complete environment for designing and analyzing
  ion/electron optics systems. Contains a lot of referential code for this
  reason.
  
  Attributes:
      field (ElectricField): The electric field initialized in the simulation.
      tracer (ParticleTracer): The trajectory calcuation tracer.
      elements (list): List of all electrodes and lenses inside the system.
  """
  
  def __init__(self, nr: float, nz: float, axial_size: float = 0.1, radial_size: float = 0.1):
      """
      Initializes an ion optics system with dimensions with the domain (axial_size, radial_size) in 
      meters, and with a grid of (nz, nr) units.
      
      Parameters:
          nr (float): Number of grid points in the r-axis direction.
          nz (float): Number of grid points in the z-axis direction.
          axial_size (float, optional): Length of the system in the z-direction in meters.
          radial_size (float, optional): Length of the system in the r-direction in meters.
      """
      self.field = ElectricField(nz, nr, axial_size, radial_size)
      self.tracer = ParticleTracer(self.field)
      self.elements = []
      self.magnetic_lenses = None

  def add_magnetic_lens(self, config: MagneticLensConfig):
      if self.magnetic_lenses is None:
          self.magnetic_lenses = MagneticField(self)
      self.magnetic_lenses.add_magnetic_lens(config)
      self.elements.append(config)

  def add_electrode(self, config: ElectrodeConfig):
      """
      Adds a single electrode to the ion optics system when ElectrodeConfig is used
      by calling on the pre-existing add_electrode() function.
      
      Parameters:
          config (ElectrodeConfig): Configuration parameters for the electrode.
      """
      self.field.add_electrode(config)
      self.elements.append(config)
      
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
      Solves the potential field using the pre-existing
      PyAMG-based multigrid solver.
      
      Returns:
          ndarray: The solved potential field.
      """
      result = {}
      
      has_electrostatic = any(isinstance(element, ElectrodeConfig) or 
                             isinstance(element, EinzelLens) 
                             for element in self.elements)
      
      has_magnetic = any(isinstance(element, MagneticLensConfig) 
                         for element in self.elements)
      
      if has_electrostatic:
          potential = self.field.solve_potential()
          result['potential'] = potential
      
      if has_magnetic and self.magnetic_lenses is not None:
          vector_potential = self.magnetic_lenses.solve_vector_potential()
          result['vector_potential'] = vector_potential
      
      if has_magnetic:
          self.tracer.magnetic_lenses = self.magnetic_lenses
      else:
          self.tracer.magnetic_lenses = None
      
      if has_electrostatic:
          return potential
      elif has_magnetic:
          return vector_potential
      else:
          return result
      
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

  def visualize_system(self, trajectories=None, r_limits=None, figsize=(16, 10), title="Picht", display_options=None):
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
      ax_checkbox = fig.add_axes([0.82, 0.4, 0.15, 0.25])
      ax_checkbox.set_facecolor('white')
      ax_checkbox.set_title('Display Options', fontsize=12, fontweight='bold', pad=10)
      
      if display_options is None:
          display_options = [True, True, True, False]
      
      checkbox = CheckButtons(ax_checkbox, 
                             ['Lenses', 'Electric Field', 'Magnetic Field', 'Animate'],
                             display_options)
      for text in checkbox.labels:
          text.set_fontsize(11)
      self.animation = None
      self.trajectory_lines = []
      self.lens_patches = []
      self.magnetic_patches = []
      
      def draw_lenses():
          self.lens_patches.clear()
          self.magnetic_patches.clear()
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
                  elif isinstance(element, MagneticLensConfig):
                      z_start = element.start * self.field.dz
                      z_end = (element.start + element.length) * self.field.dz
                      ap_center = element.ap_start + element.ap_width / 2
                      r_axis = ap_center * self.field.dr
                      bore_radius = (element.ap_width / 2) * self.field.dr
                      outer_radius = (element.outer_diameter / 2) * self.field.dr
                      
                      r_inner_bottom = max(0, r_axis - bore_radius)
                      r_outer_bottom = max(0, r_axis - outer_radius)
                      r_inner_top = min(self.field.radial_size, r_axis + bore_radius)
                      r_outer_top = min(self.field.radial_size, r_axis + outer_radius)
                      
                      patch1 = ax_main.fill_between([z_start, z_end], [r_outer_bottom, r_outer_bottom], 
                                                   [r_inner_bottom, r_inner_bottom], color='#FF0000', alpha=1.0, linewidth=0)
                      patch2 = ax_main.fill_between([z_start, z_end], [r_inner_top, r_inner_top], 
                                                   [r_outer_top, r_outer_top], color='#FF0000', alpha=1.0, linewidth=0)
                      self.magnetic_patches.extend([patch1, patch2])
      
      def draw_trajectories():
          """
          Draws/redraws static trajectories for particles when the animated tickbox is left unchecked.
          """
          self.trajectory_lines.clear()           
          if trajectories and not checkbox.get_status()[3]:  
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
          if checkbox.get_status()[3] and trajectories:  
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
          ax_main.set_title(title, fontsize=20, fontweight='bold', pad=20)
          ax_main.set_xlabel('z position (m)', fontsize=12)
          ax_main.set_ylabel('r position (m)', fontsize=12)
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
          
          if checkbox.get_status()[2] and self.magnetic_lenses is not None:  
              z_coords = np.linspace(0, self.field.axial_size, 200)
              r_coords = np.linspace(0, self.field.radial_size, 100)
              Z, R = np.meshgrid(z_coords, r_coords)
              B_z = np.zeros_like(Z)
              for i in range(Z.shape[0]):
                  for j in range(Z.shape[1]):
                      Bz, Br = self.magnetic_lenses.get_field_at_position(Z[i,j], R[i,j])
                      B_z[i,j] = Bz
              ax_main.contourf(Z, R, B_z, levels=50, cmap='RdBu', alpha=1)
              
          draw_lenses()    
          if checkbox.get_status()[3]:
              setup_animation()
          else:
              draw_trajectories()
          fig.canvas.draw_idle()
          
      checkbox.on_clicked(on_check)      
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
      
      if checkbox.get_status()[2] and self.magnetic_lenses is not None:  
          z_coords = np.linspace(0, self.field.axial_size, 200)
          r_coords = np.linspace(0, self.field.radial_size, 100)
          Z, R = np.meshgrid(z_coords, r_coords)
          B_z = np.zeros_like(Z)
          for i in range(Z.shape[0]):
              for j in range(Z.shape[1]):
                  Bz, Br = self.magnetic_lenses.get_field_at_position(Z[i,j], R[i,j])
                  B_z[i,j] = Bz
          ax_main.contourf(Z, R, B_z, levels=50, cmap='RdBu', alpha=1)
          
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
        field: Electric field
        magnetic_lenses: Magnetic field
    """
    def __init__(self, ion_optics_system):
        """
        Creates the export handler by initializing/pulling from the IonOpticsSystem.
        
        Parameters:
            ion_optics_system: IonOpticsSystem instance to export data from
        """
        self.system = ion_optics_system
        self.field = ion_optics_system.field
        self.magnetic_lenses = ion_optics_system.magnetic_lenses
   
    def export_traj(self, trajectories):
        """
        Exports simulation data to HDF5 using the H5PY library, including lens
        configuration data, particle trajectory data, and E/B-field data.
        
        Parameters:
            trajectories: List of scipy.integrate.OdeResult objects from simulate_beam
            
        Output:
            simulation_results.h5, which contains
            - Electric field data (electric potential, Ez, Er)
            - Magnetic field data (vector potential, Bz, Br)
            - Particle trajectories (pz and pr, momentum vs time)
            - System configuration (grid information, lenses, particle properties)
        """
        import h5py
        
        with h5py.File("simulation_results.h5", 'w') as f:
            field_group = f.create_group('fields')
            field_group.create_dataset('potential', data=self.field.potential)
            field_group.create_dataset('Ez', data=self.field.Ez)
            field_group.create_dataset('Er', data=self.field.Er)
            
            if self.magnetic_lenses is not None:
                mag_field_group = f.create_group('magnetic_fields')
                mag_field_group.create_dataset('vector_potential', data=self.magnetic_lenses.vector_potential)
                mag_field_group.create_dataset('Bz', data=self.magnetic_lenses.Bz)
                mag_field_group.create_dataset('Br', data=self.magnetic_lenses.Br)
                mag_field_group.create_dataset('mu_r', data=self.magnetic_lenses.mu_r)
                mag_field_group.create_dataset('current_density', data=self.magnetic_lenses.current_density)
            
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
            mag_lens_group = f.create_group('magnetic_lenses')
            
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
                        
                elif isinstance(element, MagneticLensConfig):
                    mag_lens = mag_lens_group.create_group(f'magnetic_lens_{i}')
                    mag_lens.attrs['start'] = element.start
                    mag_lens.attrs['length'] = element.length
                    mag_lens.attrs['ap_start'] = element.ap_start
                    mag_lens.attrs['ap_width'] = element.ap_width
                    mag_lens.attrs['outer_diameter'] = element.outer_diameter
                    mag_lens.attrs['mu_r'] = element.mu_r
                    mag_lens.attrs['mmf'] = element.mmf
                    
                    if hasattr(element, 'r_axis'):
                        mag_lens.attrs['r_axis'] = element.r_axis
                    
                    if hasattr(element, 'bore_radius'):
                        mag_lens.attrs['bore_radius'] = element.bore_radius
                    if hasattr(element, 'outer_radius'):
                        mag_lens.attrs['outer_radius'] = element.outer_radius
   
    def cad_export(self):
        """
        Exports lens geometry to the .step file format, by converting the
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
                    
            elif isinstance(element, MagneticLensConfig):
                magnetic_shape = self._create_magnetic_lens_shape(element)
                shapes.append(magnetic_shape)
       
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
        return electrode
        
    def _create_magnetic_lens_shape(self, magnetic_config):
        """
        Converts 2D magnetic lens config data to 3D CAD geometric data. 
        Analagous to create_electrode_shape().
        
        Parameters:
            magnetic_config: MagneticLensConfig data in 2D axisymmetric format
            
        Returns:
            3D CAD solid representation of the magnetic lens
        """
        import cadquery as cq
        
        z_start_mm = magnetic_config.start * self.field.dz * 1000
        z_length_mm = magnetic_config.length * self.field.dz * 1000
        
        outer_radius_mm = (magnetic_config.outer_diameter * self.field.dr * 1000) / 2
        
        ap_center_mm = (magnetic_config.ap_start + magnetic_config.ap_width / 2) * self.field.dr * 1000
        inner_radius_mm = (magnetic_config.ap_width * self.field.dr * 1000) / 2
        
        r_axis = 0
        if hasattr(magnetic_config, 'r_axis'):
            r_axis = magnetic_config.r_axis * self.field.dr * 1000
        
        if inner_radius_mm > 0 and inner_radius_mm < outer_radius_mm:
            mag_lens = (cq.Workplane("XY")
                       .workplane(offset=z_start_mm)
                       .center(0, r_axis)
                       .circle(outer_radius_mm)
                       .circle(inner_radius_mm)
                       .extrude(z_length_mm))
        else:
            mag_lens = (cq.Workplane("XY")
                       .workplane(offset=z_start_mm)
                       .center(0, r_axis)
                       .circle(outer_radius_mm)
                       .extrude(z_length_mm))
                       
        return mag_lens