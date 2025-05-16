import numpy as np
from core import IonOpticsSystem, ElectrodeConfig, MagneticLensConfig
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=200, nz=400, axial_size=0.1, radial_size=0.1)

mag_config = MagneticLensConfig(
    start=100,
    length=50,  
    ap_start=80,
    ap_width=40,
    outer_diameter = 100,
    mu_r=1000,
    mmf=200
)
system.add_magnetic_lens(mag_config)

system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV=10000,
    start_z=0,
    r_range=(0.042, 0.058),
    angle_range=(0, 0),
    num_particles=20,
    simulation_time=1e-9
)

fig = system.visualize_system(trajectories=trajectories)
plt.show()