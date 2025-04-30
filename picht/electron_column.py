import numpy as np
import matplotlib.pyplot as plt
from core import IonOpticsSystem

nr = 500
nz = 100
physical_size = 0.4

system = IonOpticsSystem(nr=nr, nz=nz, physical_size=physical_size)

system.add_einzel_lens(
    position=2,
    width=6,
    aperture_center=0.1,
    aperture_width=10,
    focus_voltage=-2000
)

system.add_einzel_lens(
    position=70 * dz,
    width=6 * dz,
    aperture_center=0.1,
    aperture_width=10 * dr,
    focus_voltage=-4000
)

system.add_einzel_lens(
    position=110 * dz,
    width=6 * dz,
    aperture_center=0.1,
    aperture_width=10 * dr,
    focus_voltage=-4750
)

system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV=10000,
    start_z=0.0,
    r_range=(0, 0.2),
    angle_range=(0, 0),
    num_particles=5000,
    simulation_time=6e-9
)

system.visualize_system(trajectories=trajectories)

plt.show()
