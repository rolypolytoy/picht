# Picht
An electron optics library that uses the finite difference method (FDM) to simulate electron and ion trajectories through electrostatic lenses. Currently supports Dirichlet boundary conditions, relativistic energies, non-paraxial beam configurations, and provides pre-packaged support for cylindrical and unipotential (einzel) lenses, as well as scripting for custom lens geometries.

It exists to provide a tool that's free and open-source, easily modifiable, and just as powerful as commercial tools, but with architectural decisions that enable even greater power and accuracy, through intelligent architectural decisions and a focus constrained to electron optics- a branch of computational physics with a relatively small open-source community.

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

In this we can observe several realistic behaviors, including how the fringing fields in einzel lenses first mildly defocus and then focus the beam, the beam crossover, and spherical aberration in the lens. By default, we assume Dirichlet boundary conditions, to better simulate real electrostatic lens systems with metal boundaries. Neumann boundary conditions might provide more idealized behavior, and are the defaults in most commercial electron optics solvers (including those used in COMSOL and ANSYS Maxwell), however since Dirichlet boundary conditions effectively simulate grounded boundaries rather than infinitely extending ones, for real-life systems, this is significantly more accurate, and reduces the insidious simulation-experimental gap.

Here, we demonstrate a complete simulation of an electrostatic lens-only scanning electron microscope with full accounting of divergence post-acceleration, Wehnelt cylinders, and one condenser lens and objective lens, with the final focal length ~8.7 mm after the final lens- an entirely physically plausible number, with tight convergence. I've increased the amount of particles from 6 to 100, and increased the initial beam divergence from 0 radians to +-2 radians to more accurately model the physical 'boiloff' process of thermionic sources. Regardless, the initial beam is quite straight due to acceleration between the cathode and anode, and we get this in just a few minutes:

```python
import numpy as np
from picht import IonOpticsSystem, ElectrodeConfig
import matplotlib.pyplot as plt

system = IonOpticsSystem(nr=100, nz=600, axial_size=0.6, radial_size = 0.1) #all grid units are in mm.


wehnelt1 = ElectrodeConfig(
    start=0,
    width=30,
    ap_start=30,
    ap_width=40,
    outer_diameter = 50,
    voltage=-5100
)
wehnelt2 = ElectrodeConfig(
    start=30,
    width=5,
    ap_start=40,
    ap_width=20,
    outer_diameter = 50,
    voltage=-5100
)
system.add_electrode(wehnelt1)
system.add_electrode(wehnelt2)
anode = ElectrodeConfig(
    start=40,
    width = 2,
    ap_start=48,
    ap_width=4,
    outer_diameter = 50,
    voltage=0
)
cathode = ElectrodeConfig(
    start=22,
    width = 2,
    ap_start=50,
    ap_width=0,
    outer_diameter = 2,
    voltage=-5000
)

system.add_electrode(anode)
system.add_einzel_lens(
    position=80.0,
    width=60.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-7000
)

system.add_einzel_lens(
    position=160.0,
    width=60.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-6500
)
potential = system.solve_fields()

trajectories = system.simulate_beam(
    energy_eV= 10,  
    start_z=0.025,
    r_range=(0.0499925, 0.0500075),
    angle_range=(-2, 2),
    num_particles=6,
    simulation_time=1e-8
)

figure = system.visualize_system(
    trajectories=trajectories,
    r_limits = (0.049, 0.051))

plt.show()
```

![SEM](https://github.com/user-attachments/assets/8e4bc3db-832a-4892-869d-d16839526ebe)

We can see why we need two lenses- between the first and second lens we can place a beam-limiting aperture to thin the electron beam's width, and the second lens reduces the beam spot size considerably, and also has a focal point after its own final lens, which is necessary, since you need the focus to be outside the electron column to be able to get clear samples.

If you want to more accurate identify the focal length- you can modify the system.visualize_system to be:
```python
figure = system.visualize_system(
    trajectories=trajectories,
    r_limits = (0.049, 0.051))
```
![focus](https://github.com/user-attachments/assets/5d8518e4-04b8-4677-aba3-23a68ba41b8d)

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

## Internals

The code has several architectural decisions that make it powerful, Pythonic, and performant. The library uses Numba for compiled language-level speeds, calculates particle dynamics using the BDF solver because RK45 isn't good for stiff problems, and does not use the paraxial ray equation, but instead the Lorentz force for electrostatics, and indeed to do this, solves the full Laplacian for voltage (∇²V = 0), followed by the electric field by solving for the negative of the gradient field of voltage (E = -∇V). I also calculate relativistic corrections not using a full four-vector treatment, but by using the Lorentz factor for velocity and acceleration corrections, to get both rapid computation and accurate results in the high-keV or MeV energy regime.

In addition- we use the finite difference method (FDM) instead of the boundary element method (BEM) to allow support for non-infinite problems (ie problems with grounded boundaries), and we've got it to be computationally quick too, using vectorization instead of recursive methods, and Numba's JIT for all the computations with the greatest overhead. 


