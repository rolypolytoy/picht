# Picht
An open-source electrodynamics and electron optics tool that's intuitive, performant, and physically accurate. Great for simulating the dynamics of electrons and ions through electrostatic and magnetostatic lenses of variable geometries. Calculates and visualizes electric and magnetic fields, and electrodynamics through those fields, using multigrid methods for FDM.

## Installation
```bash
pip install picht
```
[![PyPI version](https://img.shields.io/pypi/v/picht.svg)](https://pypi.org/project/picht/) ![tests](https://github.com/rolypolytoy/picht/actions/workflows/tests.yml/badge.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15376240.svg)](https://doi.org/10.5281/zenodo.15376240)

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

In this we can observe several realistic behaviors, including how the fringing fields in einzel lenses first mildly defocus and then focus the beam, the beam crossover, and spherical aberration in the lens. By default, we assume Dirichlet boundary conditions, to better simulate real electrostatic lens systems with metal boundaries. Neumann boundary conditions might provide more idealized behavior, and are the defaults in most commercial electron optics solvers (including those used in COMSOL and ANSYS Maxwell), however since Dirichlet boundary conditions effectively simulate grounded boundaries rather than infinitely extending ones, for real-life systems, this is significantly more accurate, and reduces the insidious simulation-experimental gap.

How to get images like this, and what this represents is all in the [documentation website](https://rolypolytoy.github.io/picht/), which I recommend you check out for a longer intro to Picht.

You can also specify ions by, prior to computing the trajectories, where the below syntax is for an Na+ ion:

```python
system.tracer.set_ion('Na', charge_state=1)
```

This extends the utility of Picht from electrostatic electron optics and scanning electron microscope prototyping and simulation, to focused ion beam simulations, as well as certain varieties of LINACs, and proton injectors. For example, you can use:

```python
system.tracer.set_ion('H', charge_state=1)
```

For protons, or:

```python
system.tracer.set_ion('Ga', charge_state=1)
```
For gallium (Ga+) ions like those used in gallium FIB. It also supports helium ions, neon ions, and any combination of elements and ionic charge that exists, due to integration with the Mendeleev library, and automatic parsing of charges and atomic weights. 

Additional information, API documentation, and tutorials can be found at the [official website](https://rolypolytoy.github.io/picht/).
