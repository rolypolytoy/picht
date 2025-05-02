# Picht

A Python library for simulating electron and ion beam trajectories through electrostatic lenses. Supports unipotential (einzel) lenses, custom electrode geometries, and initial emission sizes beam divergence/convergence, for electrodynamic simulations. Includes a a solver for Laplace's equation for electrostatics ∇²V = 0 using Successive Over Relaxation. It then calculates electric fields (E = -∇V) numerically, solves for the non-magnetic Lorentz force equation, and uses BDF to solve for trajectories. Supports electrons and every ion of every element.

## Installation
```bash
pip install picht
```

## Documentation

I've made a reference example that includes several electrodes and einzel lenses in the configuration of an actual scanning electron microscope's internals. 

```python
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
    position=100.0,
    width=40.0,
    aperture_center=200.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7500
)
system.add_einzel_lens(
    position=220.0,
    width=40.0,
    aperture_center=200.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-8000
)
system.add_einzel_lens(
    position=300.0,
    width=50.0,
    aperture_center=200.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7100
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
```

The visualization is beautiful, and shows cross-over based demagnification of spherical aberrations, like in real-world electrostatic lenses. In other words, the 'thickness' of the focal point reduces, until the last focal point, which is the thinnest. However, rogue electrons do increase in every crossover, which is also indicative of real-life systems. Apertures are used to mitigate this. 

![FullElectronOptics](https://github.com/user-attachments/assets/14148cd5-3fac-41ae-9801-5afd87e06fbd)

I very highly recommend first simulating a system with 5-10 electrons, and when you're satisfied with how the beam optics looks like, increase the number of electrons to get a more accurate image. The former might take less than 30 seconds, and the latter might take a few minutes.

You can also specify ions by, prior to computing the trajectories, where the below syntax is for an Na+ ion:

```python
system.tracer.set_ion('Na', charge_state=1)
```

You can use the chemical symbol of any element, and use negative charge_state values to represent ions with a negative charge. If you don't include this, by default, it assumes it's an electron.

You can also define individual electrodes by, instead of system.add_einzel_lens(), you use:

```python
electrode1 = ElectrodeConfig(
    start=10,           
    width=5,       
    ap_start=50,        
    ap_width=20,
    outer_diameter = 30,        
    voltage=1000        
)
system.add_electrode(electrode1)
```

## Why Did You Make This?

There are very few ways to simulate electron optics that are open-source, easy to use out of the box, and powerful.
