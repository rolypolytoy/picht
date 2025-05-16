# Contributing

If you want to contribute to the codebase, here's how you go about it. First, clone the repository.

```bash
git clone https://github.com/rolypolytoy/picht.git
```

Then, navigate to /picht/core.py, because this is where the active code is. I'll go over the codebase systemically to provide transparency as to what everything does.

## Codebase

I'll go over, quickly, some common changes you might want to make and how you would go about it.

You might see a function beginning like:

```python
@nb.njit
def solve_field(potential, mask, maxiter, thresh, omega=1.9):
    for _ in range(maxiter):
        vold = potential.copy()
	...
```

This is a finite-difference-method solver using the successive over relaxation (SOR) method. You can upgrade the speed of the system by integrating multigrid methods, if computational speedups are what you're looking for. Ensure you return potential, though, as calls to this function expect that.

Additionally- if you want to implement the paraxial ray equation instead of full Lorentz handling, or want to extend the electric-field-only Lorentz force handling to the full electromagnetic treatment, look for:

```python
@nb.njit
def calc_dynamics(z, r, vz, vr, Ez, Er, qm, mass, c):
    vsq = vz**2 + vr**2
    csq = c**2
    
    gamma = 1.0 / np.sqrt(1.0 - (vsq/csq))
    
    Fz = qm * mass * Ez
    Fr = qm * mass * Er
    
    factor = gamma/(gamma + 1.0) * (1.0/csq)
    vdotF = vz*Fz + vr*Fr
    
    az = Fz/(gamma * mass) - factor * vz * vdotF/mass
    ar = Fr/(gamma * mass) - factor * vr * vdotF/mass
    
    return np.array([vz, vr, az, ar])
```

Note how we calculate gamma at every timestep to properly account for the Lorentz force. To remove relativistic dynamics, simply replace gamma with 1 to force classical dynamics, or modify the equations for Fz and Fr to add magnetic field components. Note you'll have to calculate magnetic fields entirely separately, add new classes for magnetic lenses, etc etc, but integrating it with the existing solver is relatively straightforward due to the modular structure of the codebase.

If you want to modify what numerical integration method it uses, it's pretty simple. Look inside the ParticleTracer class, at around line 173, for the trace_trajectory() class:

```python
    def trace_trajectory(self, 
                   initial_position: Tuple[float, float],
                   initial_velocity: Tuple[float, float],
                   simulation_time: float,
                   method: str = 'BDF',
                   rtol: float = 1e-8,
                   atol: float = 1e-10) -> dict:
        initial_state = [
            initial_position[0], 
            initial_position[1],
            initial_velocity[0], 
            initial_velocity[1]
        ]
```
As you can currently see, it uses the BDF numerical integrator- this is good for stiff problems, but you can replace it with RK45 for improved performance, or relax rtol and atol to 1e-6 and 1e-8 respectively for greater speed, or 1e-10 and 1e-12 for more accuracy respectively. 

Finally, if you want to make modifications to the UI, in the IonOpticsSystems class at around line 293, look for visualize_system():

```python
    def visualize_system(self, 
                       trajectories=None, 
                       r_limits=None,
                       figsize=(15, 6),
                       title="Electron Trajectories"):
        plt.figure(figsize=figsize)
        plt.title(title)
        
        if trajectories:
            colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
            for i, sol in enumerate(trajectories):
                z_traj = sol.y[0]
                r_traj = sol.y[1]
                plt.plot(z_traj, r_traj, lw=1.5, color=colors[i])
        
        plt.xlabel('z position (meters)')
        plt.ylabel('r position (meters)')
        plt.grid(True, alpha=0.3)
        
        if r_limits:
            plt.ylim(r_limits)
            
        plt.tight_layout()
        return plt.gcf()
```

Since the visualization relies quite heavily on matplotlib.pyplot, modify the display however you wish, for example by changing the x and y-labels, on which axis r and x are plotted, the scaling of the image, and other details. 

These are the first places to start, and the places where the greatest improvements could be made with minimal overhead. Begin at these places, modify something, run the unit tests, and see when things break, when they change, to get a better feel for the code. After that, feel free to modify larger pieces of code in your own forks and make pull requests.

Note: some code snippets are from old versions, and may not reflect how the current codebase looks like.
