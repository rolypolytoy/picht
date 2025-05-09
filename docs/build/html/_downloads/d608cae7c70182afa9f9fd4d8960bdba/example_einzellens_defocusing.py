"""
Unipotential Lenses: Defocusing
--------------------------------
Here, we can see valid focusing behavior, which happens when the focus voltage is negative, and the particle is negatively charged, or vice versa. 
We can thus make an einzel lens that causes divergent behavior by re-using the same code but flipping the sign on focus_voltage from -5000 to 5000.
We can see clearly visible defocusing, demonstrating how einzel lenses are cleanly usable for both focusing and defocusing with single-character changes in the code:
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
    trajectories=trajectories)

plt.show()