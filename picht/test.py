import numpy as np
from picht import ElectronOptics, ElectrodeConfig, MagneticLensConfig, Export
import matplotlib.pyplot as plt

system = ElectronOptics(nr=100, nz=400, axial_size=0.4, radial_size = 0.1)

mag_config = MagneticLensConfig(
    start=100,
    length=100,  
    ap_start=44,
    ap_width=10,
    outer_diameter = 100,
    mu_r=1000,
    mmf=200
)
system.add_magnetic_lens(mag_config)

system.solve_fields()

#Electron Gun's Beam.
trajectories = system.simulate_beam(
    energy_eV= 9800,  
    start_z=0.0456210,
    r_range=(0.04996, 0.05004),
    angle_range=(-0.15, 0.15),
    num_particles=10, 
    simulation_time=2e-9
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()

#exporter = Export(system)
#exporter.cad_export()
