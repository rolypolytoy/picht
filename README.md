# Picht

A Python library for simulating electron beam trajectories through electrostatic lenses in the axisymmetric view, currently in its alpha. Supports unipotential (einzel) lenses, custom electrode geometries, and initial emission sizes beam divergence/convergence, for electrodynamic simulations. Currently supports electrons, is relativistic, and has a solver for Laplace's equation for electrostatics ∇²V = 0 using the Jacobi method. It then calculates electric fields (E = -∇V) numerically, solves for the non-magnetic Lorentz force equation, and uses BDF to solve for trajectories.

## Installation
```bash
pip install picht
```
That's it.

## How Do I Make Unipotential Lenses?

Unipotential (or einzel) lenses are amongst the simplest kind of lenses to compute the electrodynamics of. To make unipotential lenses using Picht, use (or repurpose) the following example code. Adjust the parameters of system.add_einzel_lens() to adjust the dimensions of your unipotential lens, and adjust system.simulate_beam() to adjust the parameters of the electron beam. Note that, for unipotential lenses, only the middle electrode is adjustable- the first and third electrodes are at ground. 

The system below shows a unipotential lens with 500V of focusing power, and how it impacts electrons moving with an energy of 10 keV. nx represents the amount of grid-cells in the x-dimension, ny represents the amount of grid-cells in the y-dimension and physical_size represents the dimensions of the system in meters (0.1 represents a 10cm x 10cm system). 

```python
import numpy as np
import matplotlib.pyplot as plt
from picht import IonOpticsSystem

system = IonOpticsSystem(nx=200, ny=100, physical_size=0.1)
#each x-unit is 0.1m/200 = 0.5 mm. each y-unit is 0.1m/100 = 1 mm

#Parameters of Unipotential Lens.

system.add_einzel_lens(
    position=100, #number of grid-units from x = 0
    width=10, #x-axis thickness of the full einzel lens
    aperture_center=50, #y-midpoint of the aperture of the einzel lens
    aperture_width=10, #y-thickness of the aperture of the einzel lens
    focus_voltage=-500 #voltage of the middle lens in volts.
)

system.solve_fields()

#parameters of the Electron Beam
trajectories = system.simulate_beam(
    energy_eV=10000, #energy of the electron beam
    start_x=0, #where it begins in x-coordinates
    y_range=(0.0499925, 0.0500075), #lower and upper limit of the electrons, in the y-axis, in meters
    angle_range=(0, 0), #lower and upper limit of angles of the electrons, in radians
    num_particles=20, #number of discrete electrons
    simulation_time=2e-9 #timescale of the simulation in seconds
)

system.visualize_system(
    trajectories=trajectories,
    y_limits=(49.9925, 50.0075) #y-limits of the trajectory in grid units
)

plt.show()
```

This will then display the electron trajectories in a matplotlib-style image:
![Figure_1](https://github.com/user-attachments/assets/4cf887fa-c9cb-4e6a-9aec-a8a68c11b858)
You can see the spherical aberrations in this image- how the focal length of inner and outer electrons differ. This is a major reason as to why crossovers are used in real electron optics- to reduce the size of spherical aberrations.

This is an example of a system of lenses with three unipotential lenses in an array, with 100 electrons being tracked:
```python
import numpy as np
import matplotlib.pyplot as plt
from picht import IonOpticsSystem
#100 to 500 nm spot size.

system = IonOpticsSystem(nx=500, ny=100, physical_size=0.4)

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
    start_x=0,
    y_range=(0.1999925, 0.2000075),
    angle_range=(0, 0),
    num_particles=100,
    simulation_time=6e-9
)

system.visualize_system(
    trajectories=trajectories,
    y_limits=(49.998125, 50.001875)
)

plt.show()
```

The visualization is beautiful, and shows cross-over based demagnification of spherical aberrations, like in real-world electrostatic lenses:
![Electron_Trajectory](https://github.com/user-attachments/assets/c7624809-dc87-4094-83ca-65bb778f4e36)
Over 80 out of the 100 electrons focus at the focal point around 140 units, with a spot size of a few hundred nanometers, coming from a 15 micrometer diameter spread initially. This is comparable to the demagnification actual scanning electron microscopes do internally, and so this system is a good reference for actual SEM geometries. 

I very highly recommend first simulating a system with 5-10 electrons, and when you're satisfied with how the beam optics looks like, increase the number of electrons to get a more accurate image.

## Why Did You Make This?

The paraxial ray equation is a second-order ODE which means it's relatively difficult to solve it analytically, and relatively easy to solve numerically (numerical integration solvers like BDF make it quite easy). There are very few open-source options for simulating electron optics, that are easy to set up, are on Python, are Pythonic in their syntax, and are easily customizable. You can vary the discretization of the grid by varying nx and ny, vary the dimensions by varying physical_size, modify the parameters of the einzel lenses, and vary the timescales you observe the particle trajectories in. This is an open-source alternative to commercial multiphysics systems, if identifying charged particle trajectories in electrostatic lenses is your problem. 
