===============
Getting Started
===============

Picht is a Python library for simulating charged particle optics. It currently supports electrons and ions from 0 eV to several TeV in energy, and enables the rapid creation of electrostatic and magnetostatic lenses, and the visualization of particle trajectories and fieldlines based on these parameterizations.
Uses the Lorentz force equation and the paraxial ray equation for speed and accuracy, has a PyAMG multigrid solver to enable rapid solutions using the finite difference method, and has extensive CPU parallelization that works on any operating system. 

It's particularly useful as an educational tool for electrodynamics, or to design electron optics systems, like electron microscopes (SEM/TEM/STEM), focused ion beam systems (Ga+/He+/Ne+), particle accelerators (LINACs, particle injectors) and mass spectrometers. Currently axisymmetric, with export options to .hdf5 and .step formats.

Installation
============

Install using PyPi:

.. code-block:: bash

   pip install picht

Example
===================

Here's a simple example that creates an Einzel lens system and simulates electron trajectories:

.. code-block:: python

    import numpy as np
    from picht import IonOpticsSystem, ElectrodeConfig
    import matplotlib.pyplot as plt

    system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)

    system.add_einzel_lens(
        position=20.0,
        width=60.0,
        aperture_center=50.0,
        aperture_width=48.0,
        outer_diameter=50.0,
        focus_voltage=-5000
    )

    potential = system.solve_fields()

    trajectories = system.simulate_beam(
        energy_eV= 100,  
        start_z=0,
        r_range=(0.0499925, 0.0500075),
        angle_range=(0, 0),
        num_particles=6,
        simulation_time=2e-9
    )

    figure = system.visualize_system(
        trajectories=trajectories)

    plt.show()


The API generally follows this principle, because systems are intended to be created once, and parameters adjusted one way or another until it finally works. The tweaking, 'turning of a knob' experience is the archetypal experience of using a computational physics tool and Picht aims to make this as painless as possible.

Basic Concepts
==============

Grid System
-----------

Picht uses a cylindrically symmetric grid system:

- **z-axis**: The optical axis, which is the direction the particle beam propagates
- **r-axis**: The radial distance from the optical axis
- Grid dimensions are specified in grid units, with physical size determined by ``axial_size`` and ``radial_size`` parameters. They're the system's reference frame and the smallest discrete units where unique potential field values are calculated.

The axisymmetric view, rather than a full 3D treatment, is done to increase how fast Picht runs while maintaining physical accuracy for rotationally symmetric electron optics systems. All relevant physics use the axisymmetric variants, which ensures full adherence to realistic behavior.

Electrode Configuration
-----------------------

Electrodes are defined using the ``ElectrodeConfig`` dataclass:

.. code-block:: python

   from picht import ElectrodeConfig

   electrode = ElectrodeConfig(
       start=10,              # Starting position on z-axis (grid units)
       width=20,              # Width on z-axis (grid units)
       ap_start=30,           # Aperture start on r-axis (grid units)
       ap_width=40,           # Aperture width on r-axis (grid units)
       outer_diameter=100,    # Full diameter (grid units)
       voltage=500            # Applied voltage (V)
   )

   system.add_electrode(electrode)

This creates a simple electrode with an inner diameter equal to ap_width * (radial_size/nr) in meters, and an outer diameter equal to outer_diameter * (radial_size/nr). Using grid units rather than meters feels clunky now, but as you use the system more you'll find it more intuitive.
This creates a rudimentary focusing/defocusing effect based on the geometry and magnitude of the voltage.

Particle Types
--------------

Picht can simulate electrons, as well as any ion species. The syntax for electrons, as well as for any elements is shown below.

.. code-block:: python

   system.tracer.set_ion('e-')
   system.tracer.set_ion('H', charge_state=1)
   system.tracer.set_ion('He', charge_state=2)
   system.tracer.set_ion('F', charge_state=-1)

You can simulate any combination of any element and charge state- simply use its chemical symbol (with the correct capitalization) and the charge you want to use. 
It uses a Mendeleev backend to make sure even physically unrealistic combinations can be used, so feel free to use H2000+ or the like. 

Advanced Features
=================

Custom Lens Configurations
--------------------------

Create complex electrode arrangements by stacking multiple electrodes:

.. code-block:: python

    #Wehnelt Cylinders- important components in electron guns
    wehnelt1 = ElectrodeConfig(
        start=0,
        width=30,
        ap_start=30,
        ap_width=40,
        outer_diameter = 50,
        voltage=-5100 #biased at -100V in relation to the cathode
    )
    wehnelt2 = ElectrodeConfig(
        start=30,
        width=5,
        ap_start=45,
        ap_width=10,
        outer_diameter = 50,
        voltage=-5100 #biased at -100V in relation to the cathode
    )
    system.add_electrode(wehnelt1)
    system.add_electrode(wehnelt2)


This creates a sequential electrode geometry of a hollow cylinder with a base and aperture. If you're unsure how the system works, just copy-paste this code, run it, tweak some values, and see how the UI changes in response.

Einzel Lenses
----------------------
You can also add one or multiple einzel lenses. Einzel lenses are three cylindrical electrodes in sequence- the first and third at 0V, and the middle at a fixed voltage you select. 
Einzel lenses are superior in focusing properties than single cylindrical electrodes, because they keep the beam's energy identical, and provide smoother focusing/defocusing. 
You can, of course, chain them in sequence, to make analogous 'electron optics' systems with full concave and convex lens analogues with defocusing and focusing Einzel lenses:

.. code-block:: python

    #Condenser Lens
    system.add_einzel_lens(
        position= 70.0,
        width=70.0,
        aperture_center=50.0,
        aperture_width=48.0,
        outer_diameter=50.0,
        focus_voltage=-7200
    )

    #Objective Lens
    system.add_einzel_lens(
        position= 142.0,
        width=63.0,
        aperture_center=50.0,
        aperture_width=48.0,
        outer_diameter=50.0,
        focus_voltage=-10000
    )

Interactive Visualization
-------------------------

The visualization uses Matplotlib and provides interactive checkboxes for greater control of your own user experience:

- **Lenses**: Toggles electrode visibility
- **Electric Field**: Shows electric field magnitude as a color map (turbo)
- **Magnetic Field**: Shows magnetic field magnitude as a color map (RdBu)
- **Animate**: Animates particle trajectories

You can also zoom into or out of any points by clicking the four-arrow button, right clicking the trajectories, and dragging outwards. You can use this for fine-grained aberration and focal length analysis. As you get more familiar with Matplotlib's visualization system you'll find analysis to be a lot easier.

Exporting Results
=================

Trajectory Data (HDF5)
----------------------

The following code saves your trajectory information as a .hdf5 file in the same directory as the script it was run in, for data analysis of the field and particle information. Make sure to put this after the trajectories are computed.

.. code-block:: python

   from picht import Export

   exporter = Export(system)
   exporter.export_traj(trajectories)

CAD Export (STEP)
-----------------

The following code saves your electrodes as parametric CAD files in the .step format, to directly manufacture the prototyped designs you have. Make sure to define the electrodes prior to running this. This often creates error messages but successfully exports it, so if you see a save.step file, ignore any prints in the command line.

.. code-block:: python

   exporter = Export(system) 
   exporter.cad_export()

Examples
========

Here's a complex example where we initialize an electron beam at only 0.5 eV with a high divergence. I use empirical values for thermionic emission, create a Wehnelt Cylinder, and a cathode and anode with voltages characteristic of real systems. 
Wehnelt1 and Wehnelt2 accelerate and "push" the electrons out from where they are, and the fields between the cathode at -5000V and the anode at 0V creates an incredible accelerant that effectively accelerates the electrons to 5keV of speed. 
This effectively demonstrates how Picht can be used outside its 'intended' parameterizations- i.e. the beam's initial energy and parameterization is not at all its sole state, and it calculates all physical properties realistically.

The beauty of Picht is because it runs fast and is physically realistic, a lot of the process is tweaking the same number a few hundred volts, or a few grid units here or there, and seeing the often substantial differences in results. You're now ready to play around with Picht.
Note that Picht will always be backwards compatible, but the UI might change as improvements in the interface happen. Good luck!

Electron Extraction
--------------------------

.. code-block:: python
    
    #Wehnelt Cylinders- responsible for the first crossover
    wehnelt1 = ElectrodeConfig(
        start=0,
        width=30,
        ap_start=30,
        ap_width=40,
        outer_diameter = 50,
        voltage=-5100 #biased at -100V in relation to the cathode
        )
    wehnelt2 = ElectrodeConfig(
        start=30,
        width=5,
        ap_start=45,
        ap_width=10,
        outer_diameter = 50,
        voltage=-5100 #biased at -100V in relation to the cathode
    )
    system.add_electrode(wehnelt1)
    system.add_electrode(wehnelt2)

    #Anode- +5000V in relation to the cathode, to provide acceleration
    anode = ElectrodeConfig(
        start=50,
        width = 2,
        ap_start=49,
        ap_width=2,
        outer_diameter = 50,
        voltage=0
    )
    #Cathode- represents the thermionic tungsten filament electrons boil off from
    cathode = ElectrodeConfig(
        start=24,
        width = 1,
        ap_start=50,
        ap_width=0,
        outer_diameter = 2,
        voltage=-5000
    )
    system.add_electrode(anode)
    system.add_electrode(cathode)


API Reference
=============

For detailed API documentation, see the API documentation. For a better understanding of the physics, refer to Computational Physics. Or, just check out the Gallery with several tutorial examples with pre-generated visualizations.