import pytest
import numpy as np
from unittest.mock import patch, MagicMock
import sys
import os

if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main(["-v", __file__]))

#we always pull the local version of picht, not the pypi one
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
picht_dir = os.path.join(project_root, 'picht')

if picht_dir not in sys.path:
    sys.path.insert(0, picht_dir)

from core import (
    ElectricField, 
    ElectrodeConfig, 
    EinzelLens, 
    ElectronOptics, 
    ParticleTracer,
    MagneticField,
    MagneticLensConfig,
    Export
)

@pytest.fixture
def electric_field():
    return ElectricField(100, 50, 0.1, 0.05)
@pytest.fixture
def electron_optics():
    return ElectronOptics(50, 100, 0.1, 0.05)
@pytest.fixture
def particle_tracer(electric_field):
    return ParticleTracer(electric_field)

#tests if the ElectricField class is properly using the discretization values from electron_optics
def test_electric_field_initialization(electric_field):
    assert electric_field.nz == 100
    assert electric_field.nr == 50
    assert electric_field.axial_size == 0.1
    assert electric_field.radial_size == 0.05
    assert electric_field.dz == 0.001
    assert electric_field.dr == 0.001
    assert electric_field.potential.shape == (100, 50)
    assert electric_field.Ez.shape == (100, 50)
    assert electric_field.Er.shape == (100, 50)

#tests if the dirichlet masking of electrodes is correctly done
def test_electrode_addition(electric_field):
    config = ElectrodeConfig(
        start=10,
        width=5,
        ap_start=10,
        ap_width=5,
        outer_diameter=30,
        voltage=100
    )
    electric_field.add_electrode(config)
    assert np.any(electric_field.potential != 0)
    assert np.any(electric_field.electrode_mask)

#tests if the grounded, voltage, grounded electrode configuration of einzel lenses is correctly done
def test_einzel_lens_creation():
    lens = EinzelLens(
        position=20,
        width=15,
        aperture_center=25,
        aperture_width=5,
        outer_diameter=30,
        focus_voltage=1000,
        gap_size=1
    )
    assert lens.electrode1.voltage == 0
    assert lens.electrode2.voltage == 1000
    assert lens.electrode3.voltage == 0

#tests if electrodes with nonzero voltage cause nonzero values in electric field
def test_potential_solving(electric_field):
    config = ElectrodeConfig(
        start=10,
        width=5,
        ap_start=10,
        ap_width=5,
        outer_diameter=30,
        voltage=100
    )
    electric_field.add_electrode(config)
    result = electric_field.solve_potential(max_iterations=10, convergence_threshold=1e-3)
    assert np.any(electric_field.Ez != 0)
    assert np.any(electric_field.Er != 0)
    assert result.shape == (100, 50)

#does a basic sanity test on velocity calculations- makes sure an electron at 1000eV is faster than 100eV, and 100eV is faster than 0, as well as makes sure heavier objects move slower than lighter ones at fixed energies
def test_particle_tracer_velocity_calculation(particle_tracer):
    particle_tracer.set_ion('e-')
    v1 = particle_tracer.get_velocity_from_energy(100)
    v2 = particle_tracer.get_velocity_from_energy(1000)
    assert v2 > v1
    assert v1 > 0
    
    particle_tracer.set_ion('H', 1)
    v3 = particle_tracer.get_velocity_from_energy(1000)
    assert v3 < v2

#checks if scipy's solve_ivp is being used to solve the trajectory calculations
@patch('core.solve_ivp')
def test_trace_trajectory(mock_solve_ivp, particle_tracer):
    mock_solution = MagicMock()
    mock_solve_ivp.return_value = mock_solution
    
    result = particle_tracer.trace_trajectory(
        initial_position=(0.01, 0.01),
        initial_velocity=(1e6, 0),
        simulation_time=1e-9
    )
    
    assert mock_solve_ivp.called
    assert result == mock_solution

#tests if the previously initialized electronoptics has the correct discretization and size values and is empty
def test_electron_optics_initialization(electron_optics):
    assert electron_optics.field.nz == 100
    assert electron_optics.field.nr == 50
    assert electron_optics.field.axial_size == 0.1
    assert electron_optics.field.radial_size == 0.05
    assert isinstance(electron_optics.tracer, ParticleTracer)
    assert len(electron_optics.elements) == 0

#tests if magneticfield is using the electronoptics's discretization initialization for the shape of the vector potential A's grid
def test_magnetic_field_initialization(electron_optics):
    magnetic_field = MagneticField(electron_optics)
    assert magnetic_field.nz == 100
    assert magnetic_field.nr == 50
    assert magnetic_field.axial_size == 0.1
    assert magnetic_field.radial_size == 0.05
    assert magnetic_field.vector_potential.shape == (100, 50)
    assert magnetic_field.Bz.shape == (100, 50)
    assert magnetic_field.Br.shape == (100, 50)

#tests if magnetic lens's dirichlet masking affects permeability and magnetomotive force properly
def test_magnetic_lens_addition(electron_optics):
    config = MagneticLensConfig(
        start=10,
        length=10,
        ap_start=10,
        ap_width=5,
        outer_diameter=30,
        mu_r=1000,
        mmf=1000
    )
    electron_optics.add_magnetic_lens(config)
    assert electron_optics.magnetic_lenses is not None
    assert np.any(electron_optics.magnetic_lenses.magnetic_mask)
    assert np.any(electron_optics.magnetic_lenses.mu_r != 1)
    assert np.any(electron_optics.magnetic_lenses.current_density != 0)

#tests if beam parameterization happens correctly
@patch('core.Parallel')
def test_simulate_beam(mock_parallel, electron_optics):
    mock_results = [MagicMock() for _ in range(5)]
    
    def side_effect(*args, **kwargs):
        mock_callable = MagicMock()
        mock_callable.return_value = mock_results
        return mock_callable
    
    mock_parallel.side_effect = side_effect
    
    trajectories = electron_optics.simulate_beam(
        energy_eV=1000,
        start_z=0.01,
        r_range=(0.01, 0.02),
        angle_range=(-1, 1),
        num_particles=5,
        simulation_time=1e-9
    )
    
    assert mock_parallel.called
    assert trajectories == mock_results
    assert len(trajectories) == 5