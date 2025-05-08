Picht
=====

An electron optics library that uses the finite difference method (FDM) to simulate electron and ion trajectories through electrostatic lenses. Currently supports Dirichlet boundary conditions, relativistic energies, non-paraxial beam configurations, and provides pre-packaged support for cylindrical and unipotential (einzel) lenses, as well as scripting for custom lens geometries.

It exists to provide a tool that's free and open-source, easily modifiable, and just as powerful as commercial tools, but with architectural decisions that enable even greater power and accuracy, through intelligent architectural decisions and a focus constrained to electron optics — a branch of computational physics with a relatively small open-source community. Capable of simulating electrostatic lens arrays, electron microscopes (SEM, TEM, STEM), focused ion beams (Ga/He/Ne), and mass spectrometers, and any and all applications, from the classical to relativistic regime that benefit from axisymmetric electron optics simulations without magnetic field handling.

Installation
------------

To install Picht via pip:

.. code-block:: bash

   pip install picht

.. image:: https://img.shields.io/pypi/v/picht.svg 
   :target: https://pypi.org/project/picht/ 
.. image:: https://github.com/rolypolytoy/picht/actions/workflows/tests.yml/badge.svg 
   :target: https://github.com/rolypolytoy/picht/actions/workflows/tests.yml 
.. image:: https://img.shields.io/badge/License-MIT-yellow.svg 
   :target: https://opensource.org/licenses/MIT 

Basic Usage
-----------

You can initialize an einzel lens as follows:

.. code-block:: python

   import numpy as np
   from picht import IonOpticsSystem, ElectrodeConfig
   import matplotlib.pyplot as plt

   system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1) # all grid units are in mm.

   system.add_einzel_lens(
       position=20.0,
       width=60.0,
       aperture_center=50.0,
       aperture_width=48.0,
       outer_diameter=50.0,
       focus_voltage=-7000,
       gap_size = 4
   )
   potential = system.solve_fields()

   trajectories = system.simulate_beam(
       energy_eV= 10000,  
       start_z=0,
       r_range=(0.0499925, 0.0500075),
       angle_range=(0, 0),
       num_particles=6,
       simulation_time=2e-9
   )

   figure = system.visualize_system(
       trajectories=trajectories
   )

   plt.show()

This produces, in less than 30 seconds, a physically realistic map of the particle trajectories.

.. image:: https://github.com/user-attachments/assets/d5f92b58-d0d4-4d68-8d23-6b07bb790105 

In this we can observe several realistic behaviors, including how the fringing fields in einzel lenses first mildly defocus and then focus the beam, the beam crossover, and spherical aberration in the lens.

Advanced Use Case: SEM Simulation
--------------------------------

Here's a full simulation of an electrostatic lens-only scanning electron microscope (SEM):

.. code-block:: python

   import numpy as np
   from picht import IonOpticsSystem, ElectrodeConfig
   import matplotlib.pyplot as plt

   system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1) # all grid units are in mm.

   wehnelt1 = ElectrodeConfig(
       start=0,
       width=30,
       ap_start=30,
       ap_width=40,
       outer_diameter = 50,
       voltage=-5100
   )
   wehnelt2 = ElectrodeConfig(
       start=30,
       width=5,
       ap_start=40,
       ap_width=20,
       outer_diameter = 50,
       voltage=-5100
   )
   system.add_electrode(wehnelt1)
   system.add_electrode(wehnelt2)
   anode = ElectrodeConfig(
       start=40,
       width = 2,
       ap_start=48,
       ap_width=4,
       outer_diameter = 50,
       voltage=0
   )
   cathode = ElectrodeConfig(
       start=22,
       width = 2,
       ap_start=50,
       ap_width=0,
       outer_diameter = 2,
       voltage=-5000
   )

   system.add_electrode(anode)
   system.add_einzel_lens(
       position=80.0,
       width=60.0,
       aperture_center=50.0,
       aperture_width=48.0,
       outer_diameter=50.0,
       focus_voltage=-7000
   )

   system.add_einzel_lens(
       position=160.0,
       width=60.0,
       aperture_center=50.0,
       aperture_width=48.0,
       outer_diameter=50.0,
       focus_voltage=-6500
   )
   potential = system.solve_fields()

   trajectories = system.simulate_beam(
       energy_eV= 10,  
       start_z=0.025,
       r_range=(0.0499925, 0.0500075),
       angle_range=(-2, 2),
       num_particles=6,
       simulation_time=1e-8
   )

   figure = system.visualize_system(
       trajectories=trajectories
   )

   plt.show()

.. image:: https://github.com/user-attachments/assets/8e4bc3db-832a-4892-869d-d16839526ebe 

You can also zoom in on the focal point using:

.. code-block:: python

   figure = system.visualize_system(
       trajectories=trajectories,
       r_limits = (0.049, 0.051)
   )

.. image:: https://github.com/user-attachments/assets/5d8518e4-04b8-4677-aba3-23a68ba41b8d 

Ion Support
-----------

You can simulate ions by specifying them before computing trajectories:

.. code-block:: python

   system.tracer.set_ion('Na', charge_state=1)

For example:

.. code-block:: python

   system.tracer.set_ion('H', charge_state=1)
   system.tracer.set_ion('Ga', charge_state=1)

Supports helium, neon, hydrogen, sodium, gallium, and more via integration with the Mendeleev library. You can specify any combination of atoms and charges to produce all physically realistic, and even physically unrealistic ions.

Internals
---------

The code has several architectural decisions that make it powerful, Pythonic, and performant:

- Uses **Numba** for compiled language-level speeds.
- Solves the full Laplacian for voltage (∇²V = 0).
- Calculates electric field via gradient: E = -∇V.
- Uses **BDF solver** for stiff differential equations instead of RK45.
- Avoids paraxial ray approximation; uses Lorentz force directly.
- Includes relativistic corrections using Lorentz factor.
- Uses **Finite Difference Method (FDM)** instead of Boundary Element Method (BEM) for grounded boundaries.
- Vectorized operations and JIT compilation ensure speed.
- Fully unit-tested with 15+ tests covering edge cases and physical realism.