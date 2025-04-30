import numpy as np
import matplotlib.pyplot as plt
import core
from core import IonOpticsSystem
system = IonOpticsSystem(nr=500, nz=100, physical_size=0.4)

#Wehnelt Cap
electrode1 = core.ElectrodeConfig(
    start=20,
    width=10,
    ap_start=40,
    ap_width=20,
    voltage=1000
)
electrode2 = core.ElectrodeConfig(
    start=40,
    width=10,
    ap_start=40,
    ap_width=20,
    voltage=1000
)

system.add_electrode(electrode1)
system.add_electrode(electrode2)

#Electron Lenses
system.add_einzel_lens(
    position=20, 
    width=6, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-2000
)

system.add_einzel_lens(
    position=70, 
    width=6, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-4000
)

system.add_einzel_lens(
    position=110, 
    width=6, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-4750
)


system.solve_fields()

#Thermionic Filament Beam Spread
trajectories = system.simulate_beam(
    energy_eV=10000,
    start_r=0,
    z_range=(0.1999925, 0.2000075),
    angle_range=(-2, 2),
    num_particles=100,
    simulation_time=6e-9
)
system.visualize_system(
    trajectories=trajectories,
)

plt.show()