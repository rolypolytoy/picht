import numpy as np
from picht import IonOpticsSystem, ElectrodeConfig
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=200, nz=200, physical_size=1)

#Anode
anode = ElectrodeConfig(
    start=40,
    width=5,
    ap_start=95,
    ap_width=10,
    outer_diameter = 60,
    voltage=10000
)
system.add_electrode(anode)

#Wehnelt Cap
wehnelt1 = ElectrodeConfig(
    start=20,
    width=20,
    ap_start=30,
    ap_width=40,
    outer_diameter = 0,
    voltage=0
)
wehnelt2 = ElectrodeConfig(
    start=40,
    width=5,
    ap_start=45,
    ap_width=10,
    outer_diameter = 0,
    voltage=0
)
system.add_electrode(wehnelt1)
system.add_electrode(wehnelt2)


#Condenser Lens 1
system.add_einzel_lens(
    position=20.0,
    width=6.0,
    aperture_center=50.0,
    aperture_width=10.0,
    outer_diameter=30.0,
    focus_voltage=0.0
)

#Condenser Lens 2
system.add_einzel_lens(
    position=70.0,
    width=6.0,
    aperture_center=50.0,
    aperture_width=10.0,
    outer_diameter=30.0,
    focus_voltage=0.0
)

#Objective Lens
system.add_einzel_lens(
    position=110.0,
    width=6.0,
    aperture_center=50.0,
    aperture_width=10.0,
    outer_diameter=30.0,
    focus_voltage=0.0
)

potential = system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV=100,  
    start_z=0,
    r_range=(0.4999925, 0.5000075),
    angle_range=(0, 0),
    num_particles=4,
    simulation_time=1e-9
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()