"""
Advanced Use Case: SEM Simulation
----------------------------------

Here's a full simulation of an electrostatic lens-only scanning electron microscope (SEM), where we combine
electrostatic lenses, einzel lenses, and complex acceleration, focusing and defocusing behaviors in one instance:
"""

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
    num_particles=6,
    simulation_time=1e-8
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()