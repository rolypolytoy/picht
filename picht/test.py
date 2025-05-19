import numpy as np
from core import ElectronOptics, ElectrodeConfig, MagneticLensConfig, Export
import matplotlib.pyplot as plt

#Condenser Lens
system = ElectronOptics(nr=100, nz=400, axial_size=0.4, radial_size = 0.1)
system.add_einzel_lens(
    position= 60.0,
    width=60.0,
    aperture_center=50.0,
    aperture_width=80.0,
    outer_diameter=100.0,
    focus_voltage=-8000,
    gap_size = 5
)
system.solve_fields()

#Electron Gun's Beam.
trajectories = system.simulate_beam(
    energy_eV= 9800,  
    start_z=0.0456210,
    r_range=(0.00496, 0.00504),
    angle_range=(-0.15, 0.15),
    num_particles=10, 
    simulation_time=1e-8
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()

#exporter = Export(system)
#exporter.cad_export()
