import numpy as np
import matplotlib.pyplot as plt
from core import IonOpticsSystem

system = IonOpticsSystem(nr=100, nz = 500, physical_size=0.4)

system.add_einzel_lens(
    position=20, 
    width=6, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-2000
)

system.add_einzel_lens(
    position=70, 
    width=6, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-4000
)

system.add_einzel_lens(
    position=110, 
    width=6, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-4750
)


system.solve_fields()
trajectories = system.simulate_beam(
    energy_eV=10000,
    start_z=0,
    r_range=(0.1999925, 0.2000075),
    angle_range=(0, 0),
    num_particles=10,
    simulation_time=3e-9
)

system.visualize_system(
    trajectories=trajectories)

plt.show()