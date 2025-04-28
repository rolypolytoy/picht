# Picht

A Python library for simulating electron beam trajectories through unipotential lens systems for applications in electron microscopy. Currently an extremely basic prototype, and it supports custom electrode geometries and unipotential lenses.

## Installation

```bash
pip install picht
```

## How Do I Make Unipotential Lenses?

Unipotential (or einzel) lenses are amongst the simplest kind of lenses to compute the electrodynamics of. To make unipotential lenses using Picht, use (or repurpose) the following example code. Adjust the parameters of system.add_einzel_lens() to adjust the dimensions of your unipotential lens, and adjust system.simulate_beam() to adjust the parameters of the electron beam. Note that, for unipotential lenses, only the middle electrode is adjustable- the first and third electrodes are at ground. 

```python
import numpy as np
import matplotlib.pyplot as plt
from picht import IonOpticsSystem

system = IonOpticsSystem(nx=400, ny=100, physical_size=0.1)

#Parameters of Unipotential Lens in mm.

system.add_einzel_lens(
    position=100, 
    width=10, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-1000
)

system.solve_fields()

#parameters of the Electron Beam
trajectories = system.simulate_beam(
    energy_eV=10000,
    start_x=0,
    y_range=(0.0499925, 0.0500075),
    num_particles=10,
    simulation_time=1e-9
)

system.visualize_system(
    trajectories=trajectories,
    y_limits=(49.9925, 50.0075)
)

plt.show()
```

## Why Did You Make This?

Mathematical modelling is fun, and electron optics are really hard to solve analytically, and really easy to solve numerically. This is also small enough to be coded in a few hours and useful enough to package in this manner, so that it's easily re-usable.
