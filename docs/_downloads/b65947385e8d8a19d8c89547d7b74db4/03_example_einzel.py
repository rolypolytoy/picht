"""
Unipotential Lenses: Deflection
--------------------------------
Some problems that may occur when you set your focus voltage drastically above your energy in eV can be demonstrated by this example.
Here the beam energy is only 1000 eV but we set the einzel lens to a -50kV focus voltage. Here, we see beam reflection, due to the extremely strong fields coming from the unipotential lens.
Whenever chaining several einzel lenses, this problem becomes especially pertinent, so carefully tuning electron energies and einzel lens focus voltages are important. Minimize the difference
between aperture_width and outer_diameter as well, for cleaner field configurations.
"""
import numpy as np
from picht import ElectronOptics, ElectrodeConfig
import matplotlib.pyplot as plt

system = ElectronOptics(nr=100, nz=400, axial_size=0.4, radial_size = 0.1)


system.add_einzel_lens(
    position=100.0,
    width=10.0,
    aperture_center=50.0,
    aperture_width=40.0,
    outer_diameter=80.0,
    focus_voltage=-50000
)

potential = system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV= 1000,  
    start_z=0,
    r_range=(0.045, 0.055),
    angle_range=(0, 0),
    num_particles=10,
    simulation_time=1e-8
)

figure = system.visualize_system(
    trajectories=trajectories,  
    display_options=[True, False, False, False])

plt.show()