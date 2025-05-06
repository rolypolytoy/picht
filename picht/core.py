import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Union, Callable
import numba as nb
from mendeleev import element
from scipy.special import ellipk, ellipe
import math

@dataclass
class ElectrodeConfig:
    start: float
    width: float
    ap_start: float
    ap_width: float
    outer_diameter: float
    voltage: float

@nb.njit
def solve_field(potential, mask, maxiter, thresh, omega=1.9):
    for _ in range(maxiter):
        vold = potential.copy()
        for i in range(1, potential.shape[0]-1):
            for j in range(1, potential.shape[1]-1):
                if not mask[i, j]:
                    v_new = 0.25 * (
                        potential[i+1, j] + potential[i-1, j] +
                        potential[i, j+1] + potential[i, j-1]
                    )
                    potential[i, j] = vold[i, j] + omega * (v_new - vold[i, j])
        maxdiff = 0.0
        for i in range(potential.shape[0]):
            for j in range(potential.shape[1]):
                diff = abs(potential[i, j] - vold[i, j])
                if diff > maxdiff:
                    maxdiff = diff
        if maxdiff < thresh:
            break
    return potential

@nb.njit
def get_field(z, r, Ez, Er, size, dz, dr, nz, nr):
    if 0 <= z < size and 0 <= r < size:
        i = int(min(max(1, z / dz), nz - 2))
        j = int(min(max(1, r / dr), nr - 2))
        return Ez[i, j], Er[i, j]
    else:
        return 0.0, 0.0

@nb.njit
def get_bfield(z, r, Bz, size, dz, dr, nz, nr):
    if 0 <= z < size and 0 <= r < size:
        i = int(min(max(1, z / dz), nz - 2))
        j = int(min(max(1, r / dr), nr - 2))
        return Bz[i, j]
    else:
        return 0.0

@nb.njit
def calc_dynamics(z, r, vz, vr, Ez, Er, Bz, qm, mass, c):
    vsq = vz**2 + vr**2
    csq = c**2
    gamma = 1.0 / np.sqrt(1.0 - (vsq/csq))
    Fz_e = qm * mass * Ez
    Fr_e = qm * mass * Er
    Fz_m = qm * mass * vr * Bz
    Fr_m = qm * mass * (-vz * Bz)
    Fz = Fz_e + Fz_m
    Fr = Fr_e + Fr_m
    factor = gamma/(gamma + 1.0) * (1.0/csq)
    vdotF = vz*Fz + vr*Fr
    az = Fz/(gamma * mass) - factor * vz * vdotF/mass
    ar = Fr/(gamma * mass) - factor * vr * vdotF/mass
    return np.array([vz, vr, az, ar])

class PotentialField:
    def __init__(self, nz: float, nr: float, physical_size: float):
        self.nz = nz
        self.nr = nr
        self.size = physical_size
        self.dz = physical_size / nz
        self.dr = physical_size / nr
        self.potential = np.zeros((int(nz), int(nr)))
        self.electrode_mask = np.zeros((int(nz), int(nr)), dtype=bool)
        self.Ez = None
        self.Er = None
        self.Bz = None

    def add_electrode(self, config: ElectrodeConfig):
        start, width = config.start, config.width
        ap_start, ap_width = config.ap_start, config.ap_width
        outer_diameter = config.outer_diameter
        voltage = config.voltage
        ap_center = ap_start + ap_width / 2
        r_min = max(0, ap_center - outer_diameter / 2)
        r_max = min(ap_center + outer_diameter / 2, self.nr)
        self.potential[int(start):int(start+width), int(r_min):int(r_max)] = voltage
        self.electrode_mask[int(start):int(start+width), int(r_min):int(r_max)] = True
        if ap_width > 0:
            self.potential[int(start):int(start+width), int(ap_start):int(ap_start+ap_width)] = 0
            self.electrode_mask[int(start):int(start+width), int(ap_start):int(ap_start+ap_width)] = False

    def solve_potential(self, max_iterations: float = 2000, convergence_threshold: float = 1e-6):
        self.potential = solve_field(self.potential, self.electrode_mask, int(max_iterations), convergence_threshold)
        self.Ez, self.Er = np.gradient(-self.potential, self.dz, self.dr)
        return self.potential

    def get_field_at_position(self, z: float, r: float) -> Tuple[float, float]:
        return get_field(z, r, self.Ez, self.Er, self.size, self.dz, self.dr, int(self.nz), int(self.nr))

class ParticleTracer:
    ELECTRON_CHARGE = -1.602e-19
    ELECTRON_MASS = 9.11e-31
    SPEED_OF_LIGHT = 299792458.0

    def __init__(self, potential_field: PotentialField):
        self.field = potential_field
        self.current_ion = {
            'symbol': 'e-',
            'atomic_number': 0,
            'mass': self.ELECTRON_MASS,
            'charge': self.ELECTRON_CHARGE,
            'charge_mass_ratio': self.ELECTRON_CHARGE / self.ELECTRON_MASS
        }
        self.q_m = self.current_ion['charge_mass_ratio']

    def set_ion(self, symbol: str = 'e-', charge_state: float = 1):
        if symbol == 'e-':
            self.current_ion = {
                'symbol': 'e-',
                'atomic_number': 0,
                'mass': self.ELECTRON_MASS,
                'charge': self.ELECTRON_CHARGE,
                'charge_mass_ratio': self.ELECTRON_CHARGE / self.ELECTRON_MASS
            }
        else:
            elem = element(symbol)
            isotope_mass = elem.mass
            electron_charge = 1.602e-19
            ion_charge = charge_state * electron_charge
            self.current_ion = {
                'symbol': f"{symbol}{'+' if charge_state > 0 else '-'}{abs(charge_state)}",
                'atomic_number': elem.atomic_number,
                'mass': isotope_mass * 1.66053906660e-27,
                'charge': ion_charge,
                'charge_mass_ratio': ion_charge / (isotope_mass * 1.66053906660e-27)
            }
        self.q_m = self.current_ion['charge_mass_ratio']
        return self

    def get_velocity_from_energy(self, energy_eV: float) -> float:
        energy_joules = energy_eV * 1.602e-19
        mass = self.current_ion['mass']
        rest_energy = mass * self.SPEED_OF_LIGHT**2
        total_energy = rest_energy + energy_joules
        return self.SPEED_OF_LIGHT * np.sqrt(1 - (rest_energy/total_energy)**2)

    def particle_dynamics(self, t: float, state: List[float]) -> List[float]:
        z, r, vz, vr = state
        Ez, Er = self.field.get_field_at_position(z, r)
        Bz = get_bfield(z, r, self.field.Bz, self.field.size, self.field.dz, self.field.dr, int(self.field.nz), int(self.field.nr))
        return calc_dynamics(
            z, r, vz, vr,
            Ez, Er, Bz,
            self.q_m,
            self.current_ion['mass'],
            self.SPEED_OF_LIGHT
        )

    def trace_trajectory(self, initial_position: Tuple[float, float], initial_velocity: Tuple[float, float], simulation_time: float, method: str = 'BDF', rtol: float = 1e-8, atol: float = 1e-10):
        initial_state = [initial_position[0], initial_position[1], initial_velocity[0], initial_velocity[1]]
        solution = solve_ivp(self.particle_dynamics, (0, simulation_time), initial_state, method=method, rtol=rtol, atol=atol)
        return solution

class EinzelLens:
    def __init__(self, position: float, width: float, aperture_center: float, aperture_width: float, outer_diameter: float, focus_voltage: float, gap_size: int):
        electrode_thickness = (width - 3 * gap_size)/3.0
        self.electrode1 = ElectrodeConfig(position, electrode_thickness, aperture_center - aperture_width/2, aperture_width, outer_diameter, 0)
        self.electrode2 = ElectrodeConfig(position + electrode_thickness + gap_size, electrode_thickness, aperture_center - aperture_width/2, aperture_width, outer_diameter, focus_voltage)
        self.electrode3 = ElectrodeConfig(position + 2 * electrode_thickness + 2 * gap_size, electrode_thickness, aperture_center - aperture_width/2, aperture_width, outer_diameter, 0)

    def add_to_field(self, field: PotentialField):
        field.add_electrode(self.electrode1)
        field.add_electrode(self.electrode2)
        field.add_electrode(self.electrode3)

@nb.njit
def solve_vector_potential(A_phi, J_phi, mu_r_map, mask, dz, dr, maxiter=1000, tol=1e-6, omega=1.9):
    mu0 = 4 * np.pi * 1e-7
    nz, nr = A_phi.shape
    for _ in range(maxiter):
        maxdiff = 0.0
        for i in range(1, nz - 1):
            for j in range(1, nr - 1):
                if not mask[i, j]:
                    mu_eff = mu0 * mu_r_map[i, j]
                    dA_dz = A_phi[i+1, j] + A_phi[i-1, j] - 2 * A_phi[i, j]
                    dA_dr = A_phi[i, j+1] + A_phi[i, j-1] - 2 * A_phi[i, j]
                    dA_dr_linear = (A_phi[i, j+1] - A_phi[i, j-1]) / (2 * dr)
                    laplacian = dA_dz / dz**2 + dA_dr / dr**2 + dA_dr_linear / (dr + 1e-12)
                    source = -mu_eff * J_phi[i, j]
                    residual = source - laplacian
                    delta = omega * residual * dz**2 * dr**2 / (2 * (dz**2 + dr**2))
                    A_phi[i, j] += delta
                    maxdiff = max(maxdiff, abs(delta))
        if maxdiff < tol:
            break
    return A_phi

@nb.njit
def compute_Bz_from_A_phi(A_phi, dz, dr):
    nz, nr = A_phi.shape
    Bz = np.zeros_like(A_phi)
    for i in range(nz):
        for j in range(nr):
            r = j * dr
            if j + 1 >= nr or j - 1 < 0:
                dA_dr = 0
            else:
                dA_dr = (A_phi[i, j+1] - A_phi[i, j-1]) / (2 * dr)
            Bz[i, j] = A_phi[i, j] / (r + 1e-12) + dA_dr
    return Bz

class AirCoreSolenoidLens:
    def __init__(self, position_z: int, length: int, radius: int, bore_radius: int, current_density: float, turns: int = 100):
        self.position_z = position_z
        self.length = length
        self.radius = radius
        self.bore_radius = bore_radius
        self.current_density = current_density
        self.turns = turns
        self.spacing = length / turns

    def apply_field(self, field: PotentialField):
        nz = int(field.nz)
        nr = int(field.nr)
        dz = field.dz
        dr = field.dr
        J_phi = np.zeros((nz, nr))
        mu_r_map = np.ones((nz, nr))
        mask = np.zeros((nz, nr), dtype=np.bool_)
        for i in range(self.turns):
            z_index = self.position_z + int(i * self.spacing)
            for di in range(1):
                for dj in range(self.radius - self.bore_radius, self.radius):
                    zi = z_index + di
                    if zi < nz and dj < nr:
                        J_phi[zi, dj] = self.current_density
                        mask[zi, dj] = True
        A_phi = np.zeros((nz, nr))
        A_phi = solve_vector_potential(A_phi, J_phi, mu_r_map, mask, dz, dr)
        Bz = compute_Bz_from_A_phi(A_phi, dz, dr)
        if field.Bz is None:
            field.Bz = np.zeros_like(Bz)
        field.Bz += Bz

class IronCoreSolenoidLens:
    def __init__(self, position_z: int, length: int, radius: int, bore_radius: int, current_density: float, turns: int = 100, mu_r: float = 5000):
        self.position_z = position_z
        self.length = length
        self.radius = radius
        self.bore_radius = bore_radius
        self.current_density = current_density
        self.turns = turns
        self.mu_r = mu_r
        self.spacing = length / turns

    def apply_field(self, field: PotentialField):
        nz = int(field.nz)
        nr = int(field.nr)
        dz = field.dz
        dr = field.dr
        J_phi = np.zeros((nz, nr))
        mu_r_map = np.ones((nz, nr))
        mask = np.zeros((nz, nr), dtype=np.bool_)
        for i in range(self.turns):
            z_index = self.position_z + int(i * self.spacing)
            for di in range(1):
                for dj in range(self.radius - self.bore_radius, self.radius):
                    zi = z_index + di
                    if zi < nz and dj < nr:
                        J_phi[zi, dj] = self.current_density
                        mu_r_map[zi, dj] = self.mu_r
                        mask[zi, dj] = True
        A_phi = np.zeros((nz, nr))
        A_phi = solve_vector_potential(A_phi, J_phi, mu_r_map, mask, dz, dr)
        Bz = compute_Bz_from_A_phi(A_phi, dz, dr)
        if field.Bz is None:
            field.Bz = np.zeros_like(Bz)
        field.Bz += Bz

class IonOpticsSystem:
    def __init__(self, nr: float, nz: float, physical_size: float = 0.1):
        self.field = PotentialField(nz, nr, physical_size)
        self.tracer = ParticleTracer(self.field)
        self.elements = []

    def add_electrode(self, config: ElectrodeConfig):
        self.field.add_electrode(config)

    def add_einzel_lens(self, position: float, width: float, aperture_center: float, aperture_width: float, outer_diameter: float, focus_voltage: float, gap_size: int = 1):
        lens = EinzelLens(position, width, aperture_center, aperture_width, outer_diameter, focus_voltage, gap_size)
        lens.add_to_field(self.field)
        self.elements.append(lens)

    def add_air_core_solenoid_lens(self, position_z: int, length: int, radius: int, bore_radius: int, current_density: float, turns: int = 100):
        lens = AirCoreSolenoidLens(position_z, length, radius, bore_radius, current_density, turns)
        lens.apply_field(self.field)
        self.elements.append(lens)

    def add_iron_core_solenoid_lens(self, position_z: int, length: int, radius: int, bore_radius: int, current_density: float, turns: int = 100, mu_r: float = 5000):
        lens = IronCoreSolenoidLens(position_z, length, radius, bore_radius, current_density, turns, mu_r)
        lens.apply_field(self.field)
        self.elements.append(lens)

    def solve_fields(self):
        return self.field.solve_potential()

    def simulate_beam(self, energy_eV: float, start_z: float, r_range: Tuple[float, float], angle_range: tuple, num_particles: float, simulation_time: float):
        velocity_magnitude = self.tracer.get_velocity_from_energy(energy_eV)
        min_angle_rad = np.radians(angle_range[0])
        max_angle_rad = np.radians(angle_range[1])
        angles = np.linspace(min_angle_rad, max_angle_rad, int(num_particles))
        r_positions = np.linspace(r_range[0], r_range[1], int(num_particles))
        trajectories = []
        for r_pos, angle in zip(r_positions, angles):
            vz = velocity_magnitude * np.cos(angle)
            vr = velocity_magnitude * np.sin(angle)
            sol = self.tracer.trace_trajectory((start_z, r_pos), (vz, vr), simulation_time)
            trajectories.append(sol)
        return trajectories

    def visualize_system(self, trajectories=None, r_limits=None, figsize=(15, 6), title="Electron Trajectories"):
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