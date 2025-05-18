import numpy as np
from core import IonOpticsSystem, ElectrodeConfig, MagneticLensConfig, Export
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=1000, nz=4000, axial_size=0.4, radial_size = 0.01)

#Electron Gun- Finished.
wehnelt1 = ElectrodeConfig(
    start=0,
    width=300,
    ap_start=200,
    ap_width=600,
    outer_diameter = 1000,
    voltage=-10000
)
wehnelt2 = ElectrodeConfig(
    start=300,
    width=50,
    ap_start=450,
    ap_width=100,
    outer_diameter = 1000,
    voltage=-10000
)
system.add_electrode(wehnelt1)
system.add_electrode(wehnelt2)
anode = ElectrodeConfig(
    start=450,
    width = 10,
    ap_start=490,
    ap_width=20,
    outer_diameter = 1000,
    voltage=0
)
cathode = ElectrodeConfig(
    start=248,
    width = 2,
    ap_start=500,
    ap_width=0,
    outer_diameter = 6,
    voltage=-9800
)
system.add_electrode(anode)
system.add_electrode(cathode)
#Electron Gun- Finished.


system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV= 0.5,  
    start_z=0.025,
    r_range=(0.0049875, 0.0050125),
    angle_range=(-2, 2),
    num_particles=4, 
    simulation_time=1e-7
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()

#exporter = Export(system)
#exporter.cad_export()
