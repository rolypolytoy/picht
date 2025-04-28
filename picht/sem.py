import numpy as np
import matplotlib.pyplot as plt
from picht import IonOpticsSystem
#100 to 500 nm spot size.

system = IonOpticsSystem(nx=500, ny=100, physical_size=0.1)

system.add_einzel_lens(
    position=20, 
    width=8, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-5000
)

system.add_einzel_lens(
    position=70, 
    width=10, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-5000
)

system.add_einzel_lens(
    position=190, 
    width=10, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-3000
)


system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV=10000,
    start_x=0,
    y_range=(0.0499925, 0.0500075),
    num_particles=50,
    simulation_time=2e-9
)

system.visualize_system(
    trajectories=trajectories,
    y_limits=(49.9925, 50.0075)
)

plt.show()