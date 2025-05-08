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

Electrode Creation
-----------
Before creating any electrodes, we need to familiarize ourselves with a few concepts. First- whenever you initialize a system, you need to pick 
nr, nz, axial_size and radial_size. nr is the grid resolution in the r-dimension, and nz is the same for the z-dimension. Axial_size is the total 
length of the domain in meters in the z-dimension, and radial_size is the same for the r-dimension.

This means the grid resolution in the z-dimension is axial_size/nz, and in the r-dimension is radial_size/nr. You should, if using nonstandard values for these,
ensure you're aware of these so you can properly convert from grid units to meters.

You can initialize a cylindrical electrode as as follows:

.. code-block:: python

   import numpy as np
   from picht import IonOpticsSystem, ElectrodeConfig
   import matplotlib.pyplot as plt

   system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)

   electrode = ElectrodeConfig(
        start=30,
        width=5,
        ap_start=40,
        ap_width=20,
        outer_diameter = 50,
        voltage=-5000
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

.. image:: https://github.com/user-attachments/assets/493d152d-b57e-425d-bcfd-545ca7eefac4 


Note the syntax for configuring an electrode, and how we parameterize the beam's energy in eV, its start_z, and its radial range. Keep angle_range = (0, 0) for parallel beams, or (-θ, θ), where θ is your intended half-angle of divergence.
As we can see- this produces some degree of focusing behavior, but after the electrode ends, the beam divergences once more, showcasing why we don't often use cylindrical-only lens elements, despite their simplicity.

Now, we'll build an Einzel Lens at -5000V with a beam energy of 10keV electrons. 

Unipotential Lenses
--------------------------------
Einzel lenses are three cylindrical lenses arranged in a specific pattern: the first electrode is grounded (0V), the second is at its focus_voltage, and the third is also grounded at (0V). 
They're called unipotential lenses because they can provide beam focusing and defocusing without affecting the beam's net energy, which reduces (sometimes!) the behavior we saw in cylindrical electrodes.
Always make sure the aperture_center (in grid units) aligns with the center of the r_range() when you parameterize your beam.

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

.. image:: https://github.com/user-attachments/assets/c14983e2-4fa1-4191-93e1-01c6a47082c2 

Here, we can see valid focusing behavior, which happens when the focus voltage is negative, and the particle is negatively charged, or vice versa. 
We can thus make an einzel lens that causes divergent behavior by re-using the same code but flipping the sign on focus_voltage from -5000 to 5000:

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
       focus_voltage=5000
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

.. image:: https://github.com/user-attachments/assets/2f5890ed-85ba-47f0-988e-8695df49adb8 

We can see clearly visible defocusing, demonstrating how einzel lenses are cleanly usable for both focusing and defocusing with single-character changes in the code.
Some problems that may occur when you set your focus voltage drastically above your energy in eV can be demonstrated by this example:

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
       trajectories=trajectories
   )

   plt.show()

.. image:: https://github.com/user-attachments/assets/7d54258c-4a83-4e33-b315-efc9c0d04d37 

Here the beam energy is only 100 eV but we revert to a -5kV focus voltage. Here, we see beam reflection, due to the extremely strong fields coming from the unipotential lens.
Whenever chaining several einzel lenses, this problem becomes especially pertinent, so carefully tuning electron energies and einzel lens focus voltages are important. Minimize the difference
between aperture_width and outer_diameter as well, for cleaner field configurations.

Advanced Use Case: SEM Simulation
--------------------------------

Here's a full simulation of an electrostatic lens-only scanning electron microscope (SEM):

.. code-block:: python

   import numpy as np
   from picht import IonOpticsSystem, ElectrodeConfig
   import matplotlib.pyplot as plt

   system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)

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
       num_particles=100,
       simulation_time=1e-8
   )

   figure = system.visualize_system(
       trajectories=trajectories
   )

   plt.show()

.. image:: https://github.com/user-attachments/assets/8e4bc3db-832a-4892-869d-d16839526ebe 

This synthesizes electrodes, einzel lens chaining, massively larger particle amoutns and advanced beam parameterization to provide a fully functional description of a production-level SEM system.
You can also zoom in radially by system.visualize_system() to be:

.. code-block:: python

   figure = system.visualize_system(
       trajectories=trajectories,
       r_limits = (0.049, 0.051)
   )

And you get the following as a result (num_particles to 6 electrons rather than the previous 100 for speed reasons)

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
