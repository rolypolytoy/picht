"""
Unipotential Lenses: Deflection
--------------------------------
Some problems that may occur when you set your focus voltage drastically above your energy in eV can be demonstrated by this example.
Here the beam energy is only 100 eV but we revert to a -5kV focus voltage. Here, we see beam reflection, due to the extremely strong fields coming from the unipotential lens.
Whenever chaining several einzel lenses, this problem becomes especially pertinent, so carefully tuning electron energies and einzel lens focus voltages are important. Minimize the difference
between aperture_width and outer_diameter as well, for cleaner field configurations.
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
    energy_eV= 100,  
    start_z=0,
    r_range=(0.0499925, 0.0500075),
    angle_range=(0, 0),
    num_particles=6,
    simulation_time=2e-9
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()