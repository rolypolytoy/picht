import numpy as np
from core import IonOpticsSystem, ElectrodeConfig, MagneticLensConfig, Export
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=100, nz=400, axial_size=0.4, radial_size = 0.1)

wehnelt1 = ElectrodeConfig(
    start=30,
    width=30,
    ap_start=30,
    ap_width=40,
    outer_diameter = 50,
    voltage=-10000
)
wehnelt2 = ElectrodeConfig(
    start=60,
    width=5,
    ap_start=45,
    ap_width=10,
    outer_diameter = 50,
    voltage=-10000
)
system.add_electrode(wehnelt1)
system.add_electrode(wehnelt2)

anode = ElectrodeConfig(
    start=90,
    width = 1,
    ap_start=48,
    ap_width=4,
    outer_diameter = 50,
    voltage=0
)
cathode = ElectrodeConfig(
    start=54,
    width = 1,
    ap_start=50,
    ap_width=0,
    outer_diameter = 2,
    voltage=-9800
)
system.add_electrode(anode)
system.add_electrode(cathode)
condenser_core = MagneticLensConfig(
    start=130,
    length=50,  
    ap_start=45,
    ap_width=10,
    outer_diameter = 20,
    mu_r=3183,
    mmf=0
)
system.add_magnetic_lens(condenser_core)
condenser_coil = MagneticLensConfig(
    start=130,
    length=50,  
    ap_start=40,
    ap_width=20,
    outer_diameter = 26,
    mu_r=1,
    mmf=2000
)
system.add_magnetic_lens(condenser_coil)

objective_core = MagneticLensConfig(
    start=200,
    length=50,  
    ap_start=45,
    ap_width=10,
    outer_diameter = 20,
    mu_r=3183,
    mmf=0
)
system.add_magnetic_lens(objective_core)
objective_coil = MagneticLensConfig(
    start=200,
    length=50,  
    ap_start=40,
    ap_width=20,
    outer_diameter = 26,
    mu_r=1,
    mmf=2000
)
system.add_magnetic_lens(objective_coil)


system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV= 0.5,  
    start_z=0.055,
    r_range=(0.0499875, 0.0500125),
    angle_range=(-2, 2),
    num_particles=100, 
    simulation_time=2e-8
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()

#exporter = Export(system)
#exporter.cad_export()
