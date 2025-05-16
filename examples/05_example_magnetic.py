"""
Magnetic Lens Creation
-------------------
Initialize a magnetic lens with the same parameters as electrodes, but with a few new concepts.

Mu_r is the relative permeability- a dimensionless constant specific to what material the magnet is. Iron has ~1000, air has 1, higher values increase the flux density (magnetic field strength).
MMF is the magnetomotive force- sometimes referred to as ampere-turns. This is the magnetic analogue to voltage. By specifying MMF and Mu_r you control the properties of the magnet, especially since these are applicable
both to DC electromagnets and permanent magnets.

Since the paraxial ray equation can sometimes create harmless artefacts at the middle particle's trajectories, feel free to increase nr and nz more than you would for electrostatic lenses. The multigrid handler can take it.
"""

import numpy as np
from picht import IonOpticsSystem, ElectrodeConfig, MagneticLensConfig
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=200, nz=400, axial_size=0.1, radial_size=0.1)

mag_config = MagneticLensConfig(
    start=100,
    length=50,  
    ap_start=80,
    ap_width=40,
    outer_diameter = 100,
    mu_r=1000,
    mmf=200
)
system.add_magnetic_lens(mag_config)

system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV=10000,
    start_z=0,
    r_range=(0.042, 0.058),
    angle_range=(0, 0),
    num_particles=20,
    simulation_time=1e-9
)

fig = system.visualize_system(trajectories=trajectories)
plt.show()