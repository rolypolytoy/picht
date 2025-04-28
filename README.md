# Picht

A Python library for simulating electron beam trajectories through unipotential lens systems for applications in electron microscopy. Currently an extremely basic prototype, and it supports custom electrode geometries and unipotential lenses.

## Features

- Field calculation using finite difference method
- Predefined optical elements (unipotential lenses)
- Customizable electrode configurations
- Visualizations of particle trajectories

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

#Parameters of Unipotential Lens in mm. Adjust to make your custom system, or add more elements.
system.add_einzel_lens(
    position=100, 
    width=10, 
    aperture_center=50, 
    aperture_width=10, 
    focus_voltage=-1000
)

system.solve_fields()

#Parameters of the Electron Beam
trajectories = system.simulate_beam(
    energy_eV=10000,
    start_x=0,
    y_range=(0.0499, 0.0501),
    num_particles=10,
    simulation_time=1e-9
)

system.visualize_system(trajectories=trajectories)
plt.show()
```

## Core Classes

- `PotentialField`: This class calculates the electron potential fields. 
- `ParticleTracer`: This class simulates the electron trajectories inside the electron potential fields.
- `ElectrodeConfig`: This is the class which configures a single electrode.
- `EinzelLens`: This class specifically exists to represent einzel (unipotential) lenses, to speed up their simulation.
- `IonOpticsSystem`: Example class of a multi-element ion/electron optics system.

More classes will be added in the future based on utility.
