import numpy as np
from picht import IonOpticsSystem, ElectrodeConfig
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=400, nz=400, physical_size=0.4) #all grid units are in mm.


wehnelt1 = ElectrodeConfig(
    start=0,
    width=30,
    ap_start=180,
    ap_width=40,
    outer_diameter = 50,
    voltage=-10200
)
wehnelt2 = ElectrodeConfig(
    start=30,
    width=5,
    ap_start=190,
    ap_width=20,
    outer_diameter = 50,
    voltage=-10200
)
system.add_electrode(wehnelt1)
system.add_electrode(wehnelt2)
anode = ElectrodeConfig(
    start=40,
    width = 2,
    ap_start=198,
    ap_width=4,
    outer_diameter = 50,
    voltage=0
)
cathode = ElectrodeConfig(
    start=22,
    width = 2,
    ap_start=200,
    ap_width=0,
    outer_diameter = 2,
    voltage=-10000
)

system.add_electrode(anode)
system.add_einzel_lens(
    position=80.0,
    width=60.0,
    aperture_center=200.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7500
)
system.add_einzel_lens(
    position=150.0,
    width=60.0,
    aperture_center=200.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7500
)
system.add_einzel_lens(
    position=220.0,
    width=50.0,
    aperture_center=200.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7000
)
potential = system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV= 10,  
    start_z=0.025,
    r_range=(0.1999925, 0.2000075),
    angle_range=(0, 0),
    num_particles=100,
    simulation_time=2e-8
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()