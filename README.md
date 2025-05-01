# Picht

A Python library for simulating electron and ion beam trajectories through electrostatic lenses. Supports unipotential (einzel) lenses, custom electrode geometries, and initial emission sizes beam divergence/convergence, for electrodynamic simulations. Includes a a solver for Laplace's equation for electrostatics ∇²V = 0 using Successive Over Relaxation. It then calculates electric fields (E = -∇V) numerically, solves for the non-magnetic Lorentz force equation, and uses BDF to solve for trajectories. Supports electrons and every ion of every element.

## Installation
```bash
pip install picht
```

## How Do I Make Unipotential Lenses?

Unipotential (or einzel) lenses are amongst the simplest kind of lenses to compute the electrodynamics of. To make unipotential lenses using Picht, use (or repurpose) the following example code. Adjust the parameters of system.add_einzel_lens() to adjust the dimensions of your unipotential lens, and adjust system.simulate_beam() to adjust the parameters of the electron beam. Note that, for unipotential lenses, only the middle electrode is adjustable- the first and third electrodes are at ground. 

For a complex example, we have a system of three unipotential lenses in an array, with 100 electrons being tracked. We have a 15-micron thin beam initially- note the r-axis has 1e-5 meters as its multiplier in this case- and the z-axis is in meters:

```python
import numpy as np
import matplotlib.pyplot as plt
from picht import IonOpticsSystem

system = IonOpticsSystem(nr=100, nz = 500, physical_size=0.4)

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

trajectories = system.simulate_beam(
    energy_eV=10000,
    start_z=0,
    r_range=(0.1999925, 0.2000075),
    angle_range=(0, 0),
    num_particles=100,
    simulation_time=3e-9
)

system.visualize_system(
    trajectories=trajectories)

plt.show()
```

The visualization is beautiful, and shows cross-over based demagnification of spherical aberrations, like in real-world electrostatic lenses. In other words, the 'thickness' of the focal point reduces, until the last focal point, which is the thinnest:
![electron_dynamics](https://github.com/user-attachments/assets/c767f92c-fd64-4da2-a13d-962f2af2c863)

Over 80 out of the 100 electrons focus at the focal point around 140 units, with a spot size of a few hundred nanometers, coming from a 15 micrometer diameter spread initially. This is comparable to the demagnification actual scanning electron microscopes do internally, and so this system is a good reference for actual SEM geometries. 

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
    voltage=1000        
)
system.add_electrode(electrode1)
```

## Why Did You Make This?

There are very few ways to simulate electron optics that are open-source, easy to use out of the box, and powerful.