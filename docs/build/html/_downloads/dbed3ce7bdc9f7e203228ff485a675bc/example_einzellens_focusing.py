"""
Unipotential Lenses: Focusing
--------------------------------
Einzel lenses are three cylindrical lenses arranged in a specific pattern: the first electrode is grounded (0V), the second is at its focus_voltage, and the third is also grounded at (0V). 
They're called unipotential lenses because they can provide beam focusing and defocusing without affecting the beam's net energy, which reduces (sometimes!) the behavior we saw in cylindrical electrodes.
Always make sure the aperture_center (in grid units) aligns with the center of the r_range() when you parameterize your beam.

We'll first examine an example of an einzel lens used for focusing, which happens when the polarity of the focus_voltage is the same as the polarity of the charge-
either negative to negative, or positive to positive. We'll thus set it to -5000V, with initial electron energies at 10keV.
"""

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
    trajectories=trajectories)

plt.show()