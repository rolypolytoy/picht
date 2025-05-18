API Documentation
=================

Classes
------------

ElectrodeConfig
~~~~~~~~~~~~~~~

.. class:: ElectrodeConfig

   Dataclass for configuring a cylindrically symmetric electrostatic lens element.
   
   **Parameters**
   
   :param float start: Z-axis starting position in grid units
   :param float width: Z-axis width in grid units
   :param float ap_start: R-axis aperture starting position in grid units
   :param float ap_width: R-axis aperture width in grid units
   :param float outer_diameter: Full electrode diameter in grid units
   :param float voltage: Applied voltage in volts

MagneticLensConfig
~~~~~~~~~~~~~~~~

.. class:: MagneticLensConfig

   Dataclass for configuring a cylindrically symmetric magnetic lens element.
   
   **Parameters**
   
   :param float start: Position where the magnetic lens begins on the z-axis in grid units
   :param float length: Length of the magnetic lens on the z-axis in grid units
   :param float ap_start: Starting position of the aperture on the r-axis in grid units
   :param float ap_width: Width of the aperture in the lens on the r-axis in grid units
   :param float outer_diameter: Diameter of the magnetic lens on the r-axis in grid units
   :param float mu_r: Relative magnetic permeability of the lens material in dimensionless units
   :param float mmf: Magnetomotive force of the lens in ampere-turns

ElectricField
~~~~~~~~~~~~~~

.. class:: ElectricField(nz, nr, axial_size, radial_size)

   Manages electric potential fields, electrode configurations, and field solving operations.
   
   **Parameters**
   
   :param float nz: Number of z-axis grid points
   :param float nr: Number of r-axis grid points
   :param float axial_size: Physical z-axis length in meters
   :param float radial_size: Physical r-axis length in meters
   
   **Attributes**
   
   .. attribute:: nz
      :type: int
      
      Number of z-axis grid points
   
   .. attribute:: nr
      :type: int
      
      Number of r-axis grid points
   
   .. attribute:: axial_size
      :type: float
      
      Z-axis size in meters
   
   .. attribute:: radial_size
      :type: float
      
      R-axis size in meters
   
   .. attribute:: dz
      :type: float
      
      Z-axis grid spacing (meters), calculated as axial_size/nz
   
   .. attribute:: dr
      :type: float
      
      R-axis grid spacing (meters), calculated as radial_size/nr
   
   .. attribute:: potential
      :type: numpy.ndarray
      
      2D array of electric potential values in volts
   
   .. attribute:: electrode_mask
      :type: numpy.ndarray
      
      Boolean mask indicating electrode positions
   
   .. attribute:: Ez
      :type: numpy.ndarray
      
      Z-component of electric field in V/m
   
   .. attribute:: Er
      :type: numpy.ndarray
      
      R-component of electric field in V/m
   
   **Methods**
   
   .. method:: add_electrode(config)
   
      Adds an electrode to the field configuration.
      
      :param ElectrodeConfig config: Electrode configuration parameters
   
   .. method:: build_laplacian_matrix(mask, dirichlet_values=None)
   
      Constructs sparse matrix representation of the Laplacian operator for solving Laplace's equation.
      Implements Dirichlet boundary conditions (0V) at all boundaries to simulate grounded metal boundaries.
      
      :param numpy.ndarray mask: Boolean array marking electrode positions
      :param numpy.ndarray dirichlet_values: Optional potential values at masked positions
      :returns: Sparse matrix A and right-hand side vector b
      :rtype: tuple(scipy.sparse.csr_matrix, numpy.ndarray)
   
   .. method:: solve_potential(max_iterations=500, convergence_threshold=1e-6)
   
      Solves the electrostatic potential field using PyAMG multigrid methods.
      Complexity is O(N) with grid size.
      
      :param float max_iterations: Maximum solver iterations (default: 500)
      :param float convergence_threshold: Convergence criterion (default: 1e-6)
      :returns: Solved potential field array
      :rtype: numpy.ndarray
   
   .. method:: get_field_at_position(z, r)
   
      Returns electric field components at a specific position.
      
      :param float z: Z-axis position in meters
      :param float r: R-axis position in meters
      :returns: Electric field components (Ez, Er)
      :rtype: tuple(float, float)

MagneticField
~~~~~~~~~~~~~

.. class:: MagneticField(electron_optics)

   Manages magnetic vector potential fields, magnetic lens configurations, and field solving operations.
   
   **Parameters**
   
   :param ElectronOptics electron_optics: The parent electron optics system containing field dimensions and grid parameters
   
   **Attributes**
   
   .. attribute:: nz
      :type: int
      
      Number of z-axis grid points
   
   .. attribute:: nr
      :type: int
      
      Number of r-axis grid points
   
   .. attribute:: axial_size
      :type: float
      
      Z-axis size in meters
   
   .. attribute:: radial_size
      :type: float
      
      R-axis size in meters
   
   .. attribute:: dz
      :type: float
      
      Z-axis grid spacing (meters)
   
   .. attribute:: dr
      :type: float
      
      R-axis grid spacing (meters)
   
   .. attribute:: vector_potential
      :type: numpy.ndarray
      
      2D array of magnetic vector potential values in Tesla-meters
   
   .. attribute:: magnetic_mask
      :type: numpy.ndarray
      
      Boolean mask indicating magnetic material positions
   
   .. attribute:: mu_r
      :type: numpy.ndarray
      
      2D array of relative permeability values
   
   .. attribute:: current_density
      :type: numpy.ndarray
      
      2D array of current density values in amperes per square meter
   
   .. attribute:: Bz
      :type: numpy.ndarray
      
      Z-component of magnetic field in Tesla
   
   .. attribute:: Br
      :type: numpy.ndarray
      
      R-component of magnetic field in Tesla
   
   .. attribute:: lens_config
      :type: MagneticLensConfig
      
      Configuration of the magnetic lens
   
   **Methods**
   
   .. method:: add_magnetic_lens(config)
   
      Adds a magnetic lens to the field and handles all necessary calculations.
      
      :param MagneticLensConfig config: Configuration parameters for the magnetic lens
   
   .. method:: build_laplacian_matrix(mask, dirichlet_values=None)
   
      Builds a sparse matrix for the Laplacian ∇²A = -μ₀μᵣJ, and implements Neumann 
      boundary conditions at all boundaries.
      
      :param numpy.ndarray mask: Boolean array, true where magnetic materials exist
      :param numpy.ndarray dirichlet_values: Vector potential values where mask is True
      :returns: Sparse matrix A and right-hand side vector b
      :rtype: tuple(scipy.sparse.csr_matrix, numpy.ndarray)
   
   .. method:: solve_vector_potential(max_iterations=500, convergence_threshold=1e-6)
   
      Solves the magnetic vector potential field using Multigrid methods with PyAMG.
      
      :param float max_iterations: Maximum number of iterations for the solver (default: 500)
      :param float convergence_threshold: Convergence criterion for the solution (default: 1e-6)
      :returns: The solved vector potential field
      :rtype: numpy.ndarray
   
   .. method:: calculate_b_from_a()
   
      Calculates the magnetic field components from the vector potential using B = ∇ × A,
      with special handling of differentiation at boundaries and at r = 0.
   
   .. method:: get_field_at_position(z, r)
   
      Returns the magnetic field components at a specific position.
      
      :param float z: Position along the z-axis in meters
      :param float r: Position along the r-axis in meters
      :returns: The magnetic field components (Bz, Br) at the specified position
      :rtype: tuple(float, float)

ParticleTracer
~~~~~~~~~~~~~~

.. class:: ParticleTracer(electric_field)

   Handles charged particle trajectory calculations and dynamics simulations.
   
   **Parameters**
   
   :param ElectricField electric_field: Electric field for particle simulation
   
   **Class Constants**
   
   .. attribute:: ELECTRON_CHARGE
      :value: -1.60217663e-19
      
      Elementary charge in Coulombs
   
   .. attribute:: ELECTRON_MASS
      :value: 9.1093837e-31
      
      Electron rest mass in kilograms
   
   .. attribute:: SPEED_OF_LIGHT
      :value: 299792458.0
      
      Speed of light in vacuum (m/s)
   
   **Attributes**
   
   .. attribute:: field
      :type: ElectricField
      
      Associated potential field instance
   
   .. attribute:: current_ion
      :type: dict
      
      Current particle properties including symbol, atomic_number, mass, charge, and charge_mass_ratio
   
   .. attribute:: q
      :type: float
      
      Particle charge in Coulombs
   
   .. attribute:: m
      :type: float
      
      Particle mass in kilograms
   
   **Methods**
   
   .. method:: set_ion(symbol='e-', charge_state=1)
   
      Configures the tracer for a specific ion or electron.
      Integrates with Mendeleev library for atomic data.
      
      :param str symbol: Chemical symbol or 'e-' for electrons
      :param float charge_state: Ion charge state
      :returns: Self for method chaining
      :rtype: ParticleTracer
   
   .. method:: get_velocity_from_energy(energy_eV)
   
      Converts kinetic energy to velocity with relativistic corrections.
      Accurate for all energy scales from single-digit eV to GeV.
      
      :param float energy_eV: Kinetic energy in electronvolts
      :returns: Particle velocity in m/s
      :rtype: float
   
   .. method:: particle_dynamics(t, state)
   
      Calculates particle dynamics for ODE solver.
      Uses relativistic equations of motion.
      
      :param float t: Current time
      :param list state: State vector [z, r, pz, pr]
      :returns: State derivatives [vz, vr, dpz_dt, dpr_dt]
      :rtype: list
   
   .. method:: trace_trajectory(initial_position, initial_velocity, simulation_time, method='BDF', rtol=1e-9, atol=1e-12)
   
      Solves particle equations of motion in the electric field.
      
      :param tuple initial_position: Initial (z, r) position in meters
      :param tuple initial_velocity: Initial (vz, vr) velocity in m/s
      :param float simulation_time: Total simulation time in seconds
      :param str method: ODE solver method (default: 'BDF')
      :param float rtol: Relative tolerance (default: 1e-9)
      :param float atol: Absolute tolerance (default: 1e-12)
      :returns: ODE solution object
      :rtype: scipy.integrate.OdeResult

EinzelLens
~~~~~~~~~~

.. class:: EinzelLens(position, width, aperture_center, aperture_width, outer_diameter, focus_voltage, gap_size=1)

   Three-electrode einzel (unipotential) lens system for particle focusing.
   
   **Parameters**
   
   :param float position: Starting z-axis position in grid units
   :param float width: Total lens width in grid units
   :param float aperture_center: Aperture center on r-axis in grid units
   :param float aperture_width: Aperture size in grid units
   :param float outer_diameter: Electrode diameter in grid units
   :param float focus_voltage: Central electrode voltage in volts
   :param int gap_size: Inter-electrode gap size in grid units (default: 1)
   
   **Attributes**
   
   .. attribute:: electrode1
      :type: ElectrodeConfig
      
      First electrode configuration (0V)
   
   .. attribute:: electrode2
      :type: ElectrodeConfig
      
      Central electrode configuration (focus_voltage)
   
   .. attribute:: electrode3
      :type: ElectrodeConfig
      
      Third electrode configuration (0V)
   
   **Methods**
   
   .. method:: add_to_field(field)
   
      Adds all three electrodes to the electric field.
      
      :param ElectricField field: Target electric field

ElectronOptics
~~~~~~~~~~~~~~~

.. class:: ElectronOptics(nr, nz, axial_size=0.1, radial_size=0.1)

   Complete electron optics simulation environment integrating field calculations,
   particle tracing, and visualization capabilities.
   
   **Parameters**
   
   :param float nr: Number of r-axis grid points
   :param float nz: Number of z-axis grid points
   :param float axial_size: Z-axis size in meters (default: 0.1)
   :param float radial_size: R-axis size in meters (default: 0.1)
   
   **Attributes**
   
   .. attribute:: field
      :type: ElectricField
      
      Electric field instance
   
   .. attribute:: tracer
      :type: ParticleTracer
      
      Particle trajectory calculator
   
   .. attribute:: elements
      :type: list
      
      List of all system electrodes and lenses
   
   .. attribute:: magnetic_lenses
      :type: MagneticField
      
      Magnetic field instance
   
   **Methods**
   
   .. method:: add_magnetic_lens(config)
   
      Adds a magnetic lens to the system.
      
      :param MagneticLensConfig config: Magnetic lens configuration
   
   .. method:: add_electrode(config)
   
      Adds an electrode to the system.
      
      :param ElectrodeConfig config: Electrode configuration
   
   .. method:: add_einzel_lens(position, width, aperture_center, aperture_width, outer_diameter, focus_voltage, gap_size=1)
   
      Adds an einzel lens to the system.
      
      :param float position: Starting z-axis position in grid units
      :param float width: Total lens width in grid units
      :param float aperture_center: Aperture center on r-axis in grid units
      :param float aperture_width: Aperture size in grid units
      :param float outer_diameter: Electrode diameter in grid units
      :param float focus_voltage: Central electrode voltage in volts
      :param int gap_size: Inter-electrode gap size (default: 1)
   
   .. method:: solve_fields()
   
      Solves the potential and/or vector potential fields using PyAMG multigrid solver.
      
      :returns: Solved potential field, vector potential field, or a dictionary containing both
      :rtype: numpy.ndarray or dict
   
   .. method:: simulate_beam(energy_eV, start_z, r_range, angle_range, num_particles, simulation_time, n_jobs=-1)
   
      Simulates a particle beam with specified parameters.
      Uses parallel processing for multiple trajectories.
      
      :param float energy_eV: Particle kinetic energy in electronvolts
      :param float start_z: Starting z-position in meters
      :param tuple r_range: Range of initial r-positions in meters
      :param tuple angle_range: Range of initial angles in degrees
      :param float num_particles: Number of particles to simulate
      :param float simulation_time: Total simulation time in seconds
      :param int n_jobs: Number of parallel jobs (default: -1 for all cores)
      :returns: List of trajectory solutions
      :rtype: list
   
   .. method:: visualize_system(trajectories=None, r_limits=None, figsize=(16, 10), title="Picht")
   
      Creates interactive visualization of the system and particle trajectories.
      
      :param list trajectories: Optional trajectory solutions to display
      :param tuple r_limits: Optional y-axis limits in meters
      :param tuple figsize: Figure dimensions in inches
      :param str title: Plot title
      :returns: Matplotlib figure object
      :rtype: matplotlib.figure.Figure

Export
~~~~~~

.. class:: Export(electron_optics)

   Handles data export to HDF5 and CAD formats.
   
   **Parameters**
   
   :param ElectronOptics electron_optics: System instance to export
   
   **Attributes**
   
   .. attribute:: system
      
      Reference to the electron optics system
   
   .. attribute:: field
      
      Shortcut to system's potential field
   
   .. attribute:: magnetic_lenses
      
      Shortcut to system's magnetic field
   
   **Methods**
   
   .. method:: export_traj(trajectories)
   
      Exports simulation data to HDF5 format including fields, trajectories,
      system configuration, and particle properties.
      
      :param list trajectories: List of trajectory solutions
      :output: Creates 'simulation_results.h5' file
   
   .. method:: cad_export()
   
      Exports electrode geometry to STEP format for CAD integration.
      
      :output: Creates 'save.step' file
   
   .. method:: _create_electrode_shape(electrode_config)
   
      Converts 2D electrode config data to 3D CAD geometric data.
      
      :param ElectrodeConfig electrode_config: ElectrodeConfig data in 2D axisymmetric data
      :returns: 3D CAD solid representation of the electrode
   
   .. method:: _create_magnetic_lens_shape(magnetic_config)
   
      Converts 2D magnetic lens config data to 3D CAD geometric data.
      
      :param MagneticLensConfig magnetic_config: MagneticLensConfig data in 2D axisymmetric format
      :returns: 3D CAD solid representation of the magnetic lens

Misc. Functions
-----------------

.. function:: get_field(z, r, Ez, Er, axial_size, radial_size, dz, dr, nz, nr)

   Retrieves electric field values at fractional grid positions using
   nearest-neighbor interpolation.
   
   :param float z: Z-axis position in meters
   :param float r: R-axis position in meters
   :param numpy.ndarray Ez: Z-component field array
   :param numpy.ndarray Er: R-component field array
   :param float axial_size: Total z-axis size in meters
   :param float radial_size: Total r-axis size in meters
   :param float dz: Z-axis grid spacing
   :param float dr: R-axis grid spacing
   :param int nz: Number of z-axis grid points
   :param int nr: Number of r-axis grid points
   :returns: Electric field components (Ez, Er) in V/m
   :rtype: tuple(float, float)
   
.. function:: calc_dynamics(z, r, pz, pr, Ez, Er, Bz, Br, q, m, c, r_axis=0.0)

   Calculates charged particle acceleration by applying the Lorentz force
   with special-relativistic corrections. Uses energy momentum formalism for
   full eV to TeV support.
   
   :param float z: Z-axis position in meters
   :param float r: R-axis position in meters
   :param float pz: Z-axis momentum in kg⋅m/s
   :param float pr: R-axis momentum in kg⋅m/s
   :param float Ez: Z-axis electric field in V/m
   :param float Er: R-axis electric field in V/m
   :param float Bz: Z-axis magnetic field in Tesla
   :param float Br: R-axis magnetic field in Tesla
   :param float q: Particle charge in Coulombs
   :param float m: Particle mass in kilograms
   :param float c: Speed of light in m/s
   :param float r_axis: Reference r-axis position (default: 0.0)
   :returns: Array [vz, vr, dpz_dt, dpr_dt] with velocities and forces
   :rtype: numpy.ndarray
