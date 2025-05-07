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

Here, we demonstrate a complete simulation of an electrostatic lens-only scanning electron microscope with full accounting of divergence post-acceleration, Wehnelt cylinders, and one condenser lens and objective lens, with the final focal length ~8.7 mm after the final lens- an entirely physically plausible number, with tight convergence. I've increased the amount of particles from 6 to 100, and increased the initial beam divergence from 0 radians to +-2 radians to more accurately model the physical 'boiloff' process of thermionic sources. Regardless, the initial beam is quite straight due to acceleration between the cathode and anode, and we get this in just a few minutes:


We can see why we need two lenses- between the first and second lens we can place a beam-limiting aperture to thin the electron beam's width, and the second lens reduces the beam spot size considerably, and also has a focal point after its own final lens, which is necessary, since you need the focus to be outside the electron column to be able to get clear samples.

If you want to more accurate identify the focal length- you can modify the system.visualize_system to be:
```python
figure = system.visualize_system(
    trajectories=trajectories,
    r_limits = (0.049, 0.051))
```
Or in general modify r_limits = (), to provide hardcoded limits on the r-axis, to better see extremely thin beams. For example- here- you can see how the final lens has a diamond-esque shape at its focus- typical of present, but low spherical aberrations. Since you can in matplotlib, hover over a point to see its coordinates, I can use this to find that the focal spot- even without beam limiting apertures- is from r = 0.04997 to 0.05003, which is approximately 60 microns, compared to, in the first crossover, a focal spot from r = 0.048 to 0.052, which is a spot size of 4 mm, or 4000 microns. This means the first to second crossover has a demagnification of 4000/60 = 66.67, which is impressive but not abnormal for electrostatic lens systems, which often have less components and require greater voltages than electromagnetic lens systems. This means, for example, if we have a beam limiting aperture at the first crossover of 300 micrometers diameter, the final beam spot size will be 300/66.67 = 4.5 microns.  

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
system.tracer.set_ion('Ga', charge_state=3)
```
For gallium ions. It also supports helium ions, neon ions, and any combination of elements and ionic charge that exists, due to integration with the Mendeleev library, and automatic parsing of charges and atomic weights. 