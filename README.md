# Picht

A Python library for simulating the electrodynamics of electrons and ions through electrostatic lenses. Supports cylinder lenses, unipotential lenses, and variable beam diameters, as well as custom amounts of initial divergence/convergence. Parameterize the electrodes with the voltages you want, control the mesh size, and obtain accurate trajectories in the axisymmetric view.

You can use it to simulate electron optics, to help design electron microscopes, focused ion beam systems, or simply better understand how electrostatic lenses work. Most implementations use the paraxial ray equation for simplification, but this approximation deviates significantly from reality with angles disparate from the z-axis. We, instead, compute the electric fields in full, and then compute particle trajectories using an ODE solver. This, combined with relativistic corrections at every step, allows for a significant degree of accuracy for very high particle energies (keV, MeV scale and above). It's capable of simulating multi-lens systems with 100+ charged particles in seconds to minutes, with easy Pythonic syntax, and installation is extremely easy, with PyPi handling the extremely small amount of dependencies in its entirety.

## Installation
```bash
pip install picht
```

## Documentation

You can initialize an einzel lens in the following manner:

```python
import numpy as np
from picht import IonOpticsSystem, ElectrodeConfig
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1) #all grid units are in mm.

system.add_einzel_lens(
    position=20.0,
    width=60.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7000,
    gap_size = 4
)
potential = system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV= 10000,  
    start_z=0,
    r_range=(0.0499925, 0.0500075),
    angle_range=(0, 0),
    num_particles=6,
    simulation_time=2e-9
)

figure = system.visualize_system(
    trajectories=trajectories)

plt.show()
```
It produces, in less than 30 seconds, a physically realistic map of the particle trajectories:
![Einzel_Lens](https://github.com/user-attachments/assets/d5f92b58-d0d4-4d68-8d23-6b07bb790105)

In this we can observe several realistic behaviors, including how the fringing fields in einzel lenses first mildly defocus and then focus the beam, the beam crossover, and spherical aberration in the lens. By default, we assume Dirichlet boundary conditions, to better simulate real electrostatic lens systems with metal boundaries.

You can also specify ions by, prior to computing the trajectories, where the below syntax is for an Na+ ion:

```python
system.tracer.set_ion('Na', charge_state=1)
```

You can use the chemical symbol of any element, and use negative charge_state values to represent ions with a negative charge. If you don't include this, by default, it assumes it's an electron.
