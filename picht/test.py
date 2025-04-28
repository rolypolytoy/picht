import numpy as np
import matplotlib.pyplot as plt
from core import IonOpticsSystem

system = IonOpticsSystem(nx=200, ny=100, physical_size=0.1)

#Parameters of Unipotential Lens in mm.

system.add_einzel_lens(
    position=100, 
    width=10, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-500
)

system.solve_fields()

#parameters of the Electron Beam
trajectories = system.simulate_beam(
    energy_eV=10000,
    start_x=0,
    y_range=(0.0499925, 0.0500075),
    num_particles=20,
    simulation_time=2e-9
)

system.visualize_system(
    trajectories=trajectories,
    y_limits=(49.9925, 50.0075)
)

plt.show()