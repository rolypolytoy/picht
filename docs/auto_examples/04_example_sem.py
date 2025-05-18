"""
Full-Scale SEM Simulation
----------------------------------

Here's a full simulation of an electrostatic lens-only scanning electron microscope (SEM), where we combine
electrostatic lenses, einzel lenses, and complex acceleration, focusing and defocusing behaviors in one instance. We chain together several
electrodes, a condenser (einzel) lens, and an objective (einzel) lens and observe two full crossovers- that's where the beams make an X.
Note how because of the clean-ness of the design, it looks like two clean lines- this isn't because we've parameterized the beam this way, but because
of the electron optics at play. Tinker with the parameters here, see how things change. Note: You'll have to zoom in for this, because the default system is WAY too small to see.
Learn to use matplotlib's visualization tools to better understand the system.

Some design decisions we've made for full physical realism include: 0.1 eV beam initialization to mimic thermionic emission from tungsten
5kV accelerating voltage from a hairpin cathode, with -100V biased Wehnelt Cylinders.
A -7200V condenser lens and a -10,000V objective lens.
Three total crossover points of increasing tightness.
"""
import numpy as np
from picht import ElectronOptics, ElectrodeConfig
import matplotlib.pyplot as plt

system = ElectronOptics(nr=100, nz=400, axial_size=0.4, radial_size = 0.1)


#Wehnelt Cylinders- responsible for the first crossover
wehnelt1 = ElectrodeConfig(
    start=0,
    width=30,
    ap_start=30,
    ap_width=40,
    outer_diameter = 50,
    voltage=-5100 #biased at -100V in relation to the cathode
)
wehnelt2 = ElectrodeConfig(
    start=30,
    width=5,
    ap_start=45,
    ap_width=10,
    outer_diameter = 50,
    voltage=-5100 #biased at -100V in relation to the cathode
)
system.add_electrode(wehnelt1)
system.add_electrode(wehnelt2)

#Anode- +5000V in relation to the cathode, to provide acceleration
anode = ElectrodeConfig(
    start=50,
    width = 2,
    ap_start=49,
    ap_width=2,
    outer_diameter = 50,
    voltage=0
)
#Cathode- represents the thermionic tungsten filament electrons boil off from
cathode = ElectrodeConfig(
    start=24,
    width = 1,
    ap_start=50,
    ap_width=0,
    outer_diameter = 2,
    voltage=-5000
)
system.add_electrode(anode)
system.add_electrode(cathode)

#Condenser Lens- In between the first and second crossover point, provides initial focusing
system.add_einzel_lens(
    position= 70.0,
    width=70.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-8000
)

#Objective Lens- Provides final focusing mere millimeters after its end

system.add_einzel_lens(
    position= 210.0,
    width=57.0,
    aperture_center=50.0,
    aperture_width=48.0,
    outer_diameter=50.0,
    focus_voltage=-10000
)

potential = system.solve_fields()

#Notice how we initialize it at only 0.1 eV- the acceleration happens from the field lines between the cathode and anode
trajectories = system.simulate_beam(
    energy_eV= 0.1,  
    start_z=0.025, #We begin at z = 0.025, or 25 grid units in the z-direction so that there's a bit of Wehnelt Cylinder behind this
    r_range=(0.0499925, 0.0500075), #15 micron thick beam, which is a realistic amount
    angle_range=(-2, 2), #very high initial angular divergence to mimic thermionic emission
    num_particles=100, #increasing this won't improve visualization, because the beams are artificially forced into an axisymmetric path because of the electrode configurations
    simulation_time=2e-8 #empirically found value for when the full simulation completes
)

figure = system.visualize_system(
    trajectories=trajectories,  
    display_options=[True, False, False, False])

plt.show()