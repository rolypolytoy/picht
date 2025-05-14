import pytest
import numpy as np
from scipy.constants import elementary_charge, electron_mass
import matplotlib.pyplot as plt
from picht import (
    PotentialField, 
    ParticleTracer, 
    ElectrodeConfig, 
    EinzelLens, 
    IonOpticsSystem
)

class TestPotentialField:
    def test_init(self):
        field = PotentialField(100, 50, 0.1, 0.05)
        assert field.nz == 100
        assert field.nr == 50
        assert field.axial_size == 0.1
        assert field.radial_size == 0.05
        assert field.dz == 0.001
        assert field.dr == 0.001
        assert field.potential.shape == (100, 50)
        assert not np.any(field.electrode_mask)
        assert field.Ez is None
        assert field.Er is None

    def test_add_electrode(self):
        field = PotentialField(100, 50, 0.1, 0.05)
        config = ElectrodeConfig(
            start=10, 
            width=5, 
            ap_start=20, 
            ap_width=5, 
            outer_diameter=30, 
            voltage=1000
        )
        field.add_electrode(config)
        assert np.any(field.potential == 1000)
        assert np.any(field.electrode_mask)
        aperture_region = field.potential[10:15, 20:25]
        assert np.all(aperture_region == 0)
        aperture_mask = field.electrode_mask[10:15, 20:25]
        assert not np.any(aperture_mask)

    def test_solve_potential(self):
        field = PotentialField(50, 30, 0.1, 0.05)
        config = ElectrodeConfig(
            start=10, 
            width=5, 
            ap_start=10, 
            ap_width=5, 
            outer_diameter=20, 
            voltage=1000
        )
        field.add_electrode(config)
        
        potential = field.solve_potential(max_iterations=500)
        assert field.Ez is not None
        assert field.Er is not None
        electrode_region = field.potential[10:15, 2:10]
        assert np.all(electrode_region >= 950)
        electrode_region2 = field.potential[10:15, 15:19]
        assert np.all(electrode_region2 >= 950)
        center_field = field.potential[20:30, 10:20]
        assert np.any(center_field > 0)

    def test_get_field_at_position(self):
        field = PotentialField(50, 30, 0.1, 0.06)
        config = ElectrodeConfig(
            start=10, 
            width=5, 
            ap_start=10, 
            ap_width=5, 
            outer_diameter=20, 
            voltage=1000
        )
        field.add_electrode(config)
        field.solve_potential(max_iterations=200)
        Ez, Er = field.get_field_at_position(0.05, 0.03)
        assert isinstance(Ez, float)
        assert isinstance(Er, float)
        Ez, Er = field.get_field_at_position(0.2, 0.1)
        assert Ez == 0.0
        assert Er == 0.0


class TestParticleTracer:
    @pytest.fixture
    def field(self):
        field = PotentialField(50, 30, 0.1, 0.05) 
        config = ElectrodeConfig(
            start=20, 
            width=5, 
            ap_start=10, 
            ap_width=10, 
            outer_diameter=30, 
            voltage=1000
        )
        field.add_electrode(config)
        field.solve_potential(max_iterations=100)
        return field
    
    def test_init(self, field):
        tracer = ParticleTracer(field)
        assert tracer.field == field
        assert tracer.current_ion['symbol'] == 'e-'
        
    def test_set_ion(self, field):
        tracer = ParticleTracer(field)
        tracer.set_ion('e-')
        assert tracer.current_ion['symbol'] == 'e-'
        tracer.set_ion('H', 1)
        assert tracer.current_ion['symbol'] == 'H+1'
        assert tracer.current_ion['atomic_number'] == 1
        tracer.set_ion('C', -2)
        assert tracer.current_ion['symbol'] == 'C-2'
        assert tracer.current_ion['atomic_number'] == 6
        assert tracer.current_ion['charge'] < 0
        
    def test_get_velocity_from_energy(self, field):
        tracer = ParticleTracer(field)
        v_1eV = tracer.get_velocity_from_energy(1)
        v_1keV = tracer.get_velocity_from_energy(1000)
        v_1MeV = tracer.get_velocity_from_energy(1e6)
        assert v_1eV < v_1keV < v_1MeV
        assert v_1MeV < tracer.SPEED_OF_LIGHT
        tracer.set_ion('Ar', 1)
        v_ar_1keV = tracer.get_velocity_from_energy(1000)
        assert v_ar_1keV < v_1keV
        

class TestEinzelLens:
    def test_init(self):
        lens = EinzelLens(
            position=20, 
            width=15, 
            aperture_center=15,
            aperture_width=10,
            outer_diameter=30,
            focus_voltage=-2000,
            gap_size=1
        )
        
        assert lens.electrode1.start == 20
        assert lens.electrode1.voltage == 0
        assert lens.electrode2.start > lens.electrode1.start
        assert lens.electrode2.voltage == -2000
        assert lens.electrode3.start > lens.electrode2.start
        assert lens.electrode3.voltage == 0
        assert lens.electrode1.ap_width == 10
        assert lens.electrode2.ap_width == 10
        assert lens.electrode3.ap_width == 10
    
    def test_add_to_field(self):
        field = PotentialField(100, 50, 0.1, 0.05)
        lens = EinzelLens(
            position=20, 
            width=15, 
            aperture_center=25,
            aperture_width=10,
            outer_diameter=40,
            focus_voltage=-3000,
            gap_size=1
        )
        
        lens.add_to_field(field)
        
        unique_voltages = np.unique(field.potential)
        assert 0 in unique_voltages
        assert -3000 in unique_voltages
        assert np.any(field.electrode_mask)


class TestIonOpticsSystem:
    def test_init(self):
        system = IonOpticsSystem(nr=50, nz=100)
        assert isinstance(system.field, PotentialField)
        assert isinstance(system.tracer, ParticleTracer)
        assert len(system.elements) == 0
    
    def test_add_electrode(self):
        system = IonOpticsSystem(nr=50, nz=100)
        config = ElectrodeConfig(
            start=10, 
            width=5, 
            ap_start=20, 
            ap_width=5, 
            outer_diameter=30, 
            voltage=1000
        )
        system.add_electrode(config)
        assert np.any(system.field.potential == 1000)
    
    def test_add_einzel_lens(self):
        system = IonOpticsSystem(nr=50, nz=100)
        system.add_einzel_lens(
            position=20, 
            width=15, 
            aperture_center=25,
            aperture_width=10,
            outer_diameter=40,
            focus_voltage=-2000
        )
        
        assert len(system.elements) == 1
        unique_voltages = np.unique(system.field.potential)
        assert 0 in unique_voltages
        assert -2000 in unique_voltages
    
    def test_solve_fields(self):
        system = IonOpticsSystem(nr=30, nz=50)
        system.add_einzel_lens(
            position=10, 
            width=9, 
            aperture_center=15,
            aperture_width=5,
            outer_diameter=25,
            focus_voltage=-1000
        )
        
        potential = system.solve_fields()
        assert system.field.Ez is not None
        assert system.field.Er is not None
    
    def test_simulate_beam(self):
        system = IonOpticsSystem(nr=30, nz=80, axial_size=0.08, radial_size=0.04)
        system.add_einzel_lens(
            position=20, 
            width=12, 
            aperture_center=15,
            aperture_width=6,
            outer_diameter=25,
            focus_voltage=-200
        )
        system.solve_fields()
        
        trajectories = system.simulate_beam(
            energy_eV=1000,
            start_z=0.005,
            r_range=(0.005, 0.015),
            angle_range=(-5, 5),
            num_particles=3,
            simulation_time=1e-8
        )
        
        assert len(trajectories) == 3
        for traj in trajectories:
            assert hasattr(traj, 'y')
            assert hasattr(traj, 't')
            assert traj.y.shape[0] == 4
            assert traj.y[0, -1] > traj.y[0, 0]
    
    def test_visualize_system(self, monkeypatch):
        mock_called = False
        def mock_figure(*args, **kwargs):
            nonlocal mock_called
            mock_called = True
            return plt.Figure()
        
        monkeypatch.setattr(plt, 'figure', mock_figure)
        system = IonOpticsSystem(nr=20, nz=40)
        system.add_electrode(ElectrodeConfig(
            start=10, width=5, ap_start=8, ap_width=4, 
            outer_diameter=15, voltage=100
        ))
        fig = system.visualize_system()
        assert mock_called
        class MockSol:
            def __init__(self):
                self.y = [
                    np.linspace(0, 0.05, 10),
                    np.linspace(0.01, 0.02, 10),
                    np.zeros(10),
                    np.zeros(10)
                ]
                self.t = np.linspace(0, 1e-8, 10)
        fig = system.visualize_system(trajectories=[MockSol(), MockSol()])
        assert mock_called
