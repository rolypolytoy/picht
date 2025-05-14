"""
Electrode Creation
-------------------
Before creating any electrodes, we need to familiarize ourselves with a few concepts. First- whenever you initialize a system, you need to pick 
nr, nz, axial_size and radial_size. nr is the grid resolution in the r-dimension, and nz is the same for the z-dimension. Axial_size is the total 
length of the domain in meters in the z-dimension, and radial_size is the same for the r-dimension.

This means the grid resolution in the z-dimension is axial_size/nz, and in the r-dimension is radial_size/nr. You should, if using nonstandard values for these,
ensure you're aware of these so you can properly convert from grid units to meters.

You can initialize a cylindrical electrode with -5000V focus voltage and 10keV energy per electron as follows:
"""
import numpy as np
from picht import IonOpticsSystem, ElectrodeConfig
import matplotlib.pyplot as plt

#Specifies the grid size, resolution, and domain
system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1)

#Configures electrode parameterization
electrode = ElectrodeConfig(
    start=30,
    width=5,
    ap_start=40,
    ap_width=20,
    outer_diameter = 50,
    voltage=-5000
)

system.add_electrode(electrode)

#FDM field solving stage
potential = system.solve_fields()

#Beam parameterization
trajectories = system.simulate_beam(
    energy_eV= 10000,  
    start_z=0,
    r_range=(0.0499925, 0.0500075),
    angle_range=(0, 0),
    num_particles=6,
    simulation_time=2e-9
)

#Visualization control
figure = system.visualize_system(
    trajectories=trajectories)

plt.show()