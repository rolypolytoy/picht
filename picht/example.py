<<<<<<< Updated upstream
=======
"""
Unipotential Lenses: Focusing
--------------------------------
Einzel lenses are three cylindrical lenses arranged in a specific pattern: the first electrode is grounded (0V), the second is at its focus_voltage, and the third is also grounded at (0V). 
They're called unipotential lenses because they can provide beam focusing and defocusing without affecting the beam's net energy, which reduces (sometimes!) the behavior we saw in cylindrical electrodes.
Always make sure the aperture_center (in grid units) aligns with the center of the r_range() when you parameterize your beam.

We'll first examine an example of an einzel lens used for focusing, which happens when the polarity of the focus_voltage is the same as the polarity of the charge-
either negative to negative, or positive to positive. We'll thus set it to -5000V, with initial electron energies at 10keV.
"""

<<<<<<< Updated upstream
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
=======
>>>>>>> Stashed changes
import numpy as np
from core import IonOpticsSystem, ElectrodeConfig
import matplotlib.pyplot as plt

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1) #all grid units are in mm.


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
=======
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)


system.add_einzel_lens(
    position=20.0,
>>>>>>> Stashed changes
=======
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)


system.add_einzel_lens(
    position=20.0,
>>>>>>> Stashed changes
=======
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)


system.add_einzel_lens(
    position=20.0,
>>>>>>> Stashed changes
=======
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)


system.add_einzel_lens(
    position=20.0,
>>>>>>> Stashed changes
=======
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)


system.add_einzel_lens(
    position=20.0,
>>>>>>> Stashed changes
    width=60.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    focus_voltage=-7000
)

system.add_einzel_lens(
    position=160.0,
    width=60.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-6500
=======
    focus_voltage=-5000
>>>>>>> Stashed changes
=======
    focus_voltage=-5000
>>>>>>> Stashed changes
=======
    focus_voltage=-5000
>>>>>>> Stashed changes
=======
    focus_voltage=-5000
>>>>>>> Stashed changes
=======
    focus_voltage=-5000
>>>>>>> Stashed changes
)
potential = system.solve_fields()

trajectories = system.simulate_beam(
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    energy_eV= 10,  
    start_z=0.025,
    r_range=(0.0499925, 0.0500075),
    angle_range=(-2, 2),
    num_particles=6,
    simulation_time=1e-8
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    energy_eV= 10000,  
    start_z=0,
    r_range=(0.0499925, 0.0500075),
    angle_range=(0, 0),
    num_particles=6,
    simulation_time=2e-9
<<<<<<< Updated upstream
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
=======
>>>>>>> Stashed changes
)

figure = system.visualize_system(
    trajectories=trajectories,
    r_limits = (0.049, 0.051))

plt.show()
