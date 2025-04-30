import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Tuple
import numba as nb


@dataclass
class ElectrodeConfig:
    start: float
    width: float
    ap_start: float
    ap_width: float
    voltage: float


@nb.njit(parallel=True)
def laplace_iteration(potential, electrode_mask, electrode_values):
    V_old = potential.copy()
    nr, nz = potential.shape
    omega = 1.95

    for i in range(1, nr - 1):
        for j in range(1, nz - 1):
            if not electrode_mask[i, j]:
                residual = 0.25 * (
                    potential[i + 1, j] + potential[i - 1, j] +
                    potential[i, j + 1] + potential[i, j - 1]
                ) - potential[i, j]
                potential[i, j] += omega * residual

    for i in range(nr):
        for j in range(nz):
            if electrode_mask[i, j]:
                potential[i, j] = electrode_values[i, j]

    max_diff = np.max(np.abs(potential - V_old))
    return potential, max_diff


@nb.njit
def get_field_values(Er, Ez, r, z, dr, dz, nr, nz):
    i_f = r / dr
    j_f = z / dz
    i = int(i_f)
    j = int(j_f)

    if i < 0 or i >= nr - 1 or j < 0 or j >= nz - 1:
        return 0.0, 0.0

    di = i_f - i
    dj = j_f - j

    Er_interp = (
        Er[i, j] * (1 - di) * (1 - dj) +
        Er[i + 1, j] * di * (1 - dj) +
        Er[i, j + 1] * (1 - di) * dj +
        Er[i + 1, j + 1] * di * dj
    )

    Ez_interp = (
        Ez[i, j] * (1 - di) * (1 - dj) +
        Ez[i + 1, j] * di * (1 - dj) +
        Ez[i, j + 1] * (1 - di) * dj +
        Ez[i + 1, j + 1] * di * dj
    )

    return Er_interp, Ez_interp


@nb.njit
def relativity(t, state, q_m, c, Er, Ez):
    r, z, vr, vz = state
    v2 = vr * vr + vz * vz
    gamma = 1.0 / np.sqrt(1.0 - v2 / (c * c))
    ar = q_m * Er / gamma
    az = q_m * Ez / gamma
    return np.array([vr, vz, ar, az])


class PotentialField:
    def __init__(self, nr: int, nz: int, physical_size: float):
        self.nr = nr
        self.nz = nz
        self.size = physical_size
        self.dr = physical_size / nr
        self.dz = physical_size / nz
        self.potential = np.zeros((nr, nz))
        self.electrode_mask = np.zeros((nr, nz), dtype=bool)
        self.electrode_values = np.zeros((nr, nz))
        self.Er = None
        self.Ez = None

    def add_electrode(self, config: ElectrodeConfig):
        i_start = int(config.ap_start / self.dr)
        i_end = int((config.ap_start + config.ap_width) / self.dr)
        j_start = int(config.start / self.dz)
        j_end = int((config.start + config.width) / self.dz)

        self.potential[i_start:i_end, j_start:j_end] = config.voltage
        self.electrode_mask[i_start:i_end, j_start:j_end] = True
        self.electrode_values[i_start:i_end, j_start:j_end] = config.voltage

    def solve_potential(self, max_iterations: int = 2000, convergence_threshold: float = 1e-8):
        for _ in range(max_iterations):
            self.potential, max_diff = laplace_iteration(
                self.potential, self.electrode_mask, self.electrode_values
            )
            if max_diff < convergence_threshold:
                break

        self.Er, self.Ez = np.gradient(-self.potential, self.dr, self.dz)
        return self.potential

    def get_field_at_position(self, r: float, z: float) -> Tuple[float, float]:
        return get_field_values(self.Er, self.Ez, r, z, self.dr, self.dz, self.nr, self.nz)


class ParticleTracer:
    ELECTRON_CHARGE = -1.602e-19
    ELECTRON_MASS = 9.11e-31
    SPEED_OF_LIGHT = 299792458.0

    def __init__(self, potential_field: PotentialField):
        self.field = potential_field
        self.q_m = self.ELECTRON_CHARGE / self.ELECTRON_MASS
        self.rest_mass = self.ELECTRON_MASS
        self.charge = self.ELECTRON_CHARGE
        self.c = self.SPEED_OF_LIGHT

    def set_charge_mass_ratio(self, q: float, m: float):
        self.q_m = q / m
        self.rest_mass = m
        self.charge = q

    def get_velocity_from_energy(self, energy_eV: float) -> float:
        energy_j = energy_eV * 1.602e-19
        rest_energy = self.rest_mass * self.c**2
        total_energy = rest_energy + energy_j
        return self.c * np.sqrt(1 - (rest_energy / total_energy)**2)

    def particle_dynamics(self, t, state):
        r, z, vr, vz = state
        Er, Ez = self.field.get_field_at_position(r, z)
        return relativity(t, state, self.q_m, self.c, Er, Ez)

    def trace_trajectory(self, initial_position, initial_velocity, simulation_time):
        r0, z0 = initial_position
        vr0, vz0 = initial_velocity
        state0 = [r0, z0, vr0, vz0]

        return solve_ivp(
            self.particle_dynamics,
            t_span=(0, simulation_time),
            y0=state0,
            method='BDF',
            rtol=1e-8,
            atol=1e-10
        )


class EinzelLens:
    def __init__(self, position, width, aperture_center, aperture_width, focus_voltage, electrode_thickness=0.002):
        center_thickness = width - 2 * electrode_thickness

        self.electrode1 = ElectrodeConfig(
            start=position,
            width=electrode_thickness,
            ap_start=aperture_center,
            ap_width=aperture_width,
            voltage=0
        )

        self.electrode2 = ElectrodeConfig(
            start=position + electrode_thickness,
            width=center_thickness,
            ap_start=aperture_center,
            ap_width=aperture_width,
            voltage=focus_voltage
        )

        self.electrode3 = ElectrodeConfig(
            start=position + width - electrode_thickness,
            width=electrode_thickness,
            ap_start=aperture_center,
            ap_width=aperture_width,
            voltage=0
        )

    def add_to_field(self, field: PotentialField):
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)


class IonOpticsSystem:
    def __init__(self, nr: int, nz: int, physical_size: float):
        self.field = PotentialField(nr, nz, physical_size)
        self.tracer = ParticleTracer(self.field)

    def add_einzel_lens(self, position, width, aperture_center, aperture_width, focus_voltage):
        lens = EinzelLens(position, width, aperture_center, aperture_width, focus_voltage)
        lens.add_to_field(self.field)

    def solve_fields(self):
        self.field.solve_potential()

    def simulate_beam(self, energy_eV, start_z, r_range, angle_range, num_particles, simulation_time):
        v = self.tracer.get_velocity_from_energy(energy_eV)
        r_values = np.linspace(r_range[0], r_range[1], num_particles)
        angles = np.radians(np.linspace(angle_range[0], angle_range[1], num_particles))
        results = []

        for r, angle in zip(r_values, angles):
            vr = v * np.sin(angle)
            vz = v * np.cos(angle)
            sol = self.tracer.trace_trajectory((r, start_z), (vr, vz), simulation_time)
            results.append(sol)

        return results

    def visualize_system(self, trajectories):
        for sol in trajectories:
            r, z = sol.y[0], sol.y[1]
            plt.plot(z, r, lw=0.8)
        plt.xlabel('z (m)')
        plt.ylabel('r (m)')
        plt.title('Particle Trajectories')