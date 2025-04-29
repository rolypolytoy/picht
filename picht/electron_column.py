import numpy as np
import matplotlib.pyplot as plt
from picht import IonOpticsSystem
#100 to 500 nm spot size.

system = IonOpticsSystem(nx=500, ny=100, physical_size=0.4)

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
    start_x=0,
    y_range=(0.1999925, 0.2000075),
    num_particles=100,
    simulation_time=6e-9
)

system.visualize_system(
    trajectories=trajectories,
    y_limits=(49.998125, 50.001875)
)

plt.show()