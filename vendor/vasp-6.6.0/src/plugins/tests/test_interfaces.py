import functools
from dataclasses import dataclass, field
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import vasp
import vasp._machine_learning
from vasp.typing import IndexArray

from .util import mock_vasp_plugin


@dataclass(frozen=True)
class MockConstants:
    ion_types: IndexArray = field(default_factory=lambda: np.arange(1, 4))


class MockCalculator:
    def __init__(self, shape):
        self._shape = shape

    def get_potential_energy(self, *args, **kwargs):
        return 1.0

    def get_forces(self, *args, **kwargs):
        return np.ones(self._shape)

    def get_stress(self, *args, **kwargs):
        return np.ones((3, 3))


@pytest.fixture
def mock_interface(function_name):
    interface = getattr(vasp, f"interface_{function_name}")
    return functools.partial(interface, MockConstants(), MagicMock())


def test_interface(mock_interface):
    with pytest.warns(UserWarning):
        assert mock_interface() == vasp.SUCCESS


def test_vasp_plugin(mock_interface, function_name):
    with mock_vasp_plugin(function_name) as function:
        assert mock_interface() == vasp.SUCCESS
        function.assert_called_once_with(*mock_interface.args)


def test_vasp_plugin_without_function(mock_interface, function_name):
    with pytest.warns(UserWarning), mock_vasp_plugin(function_name, function=None):
        assert mock_interface() == vasp.SUCCESS


def test_entry_points(mock_entry_points, mock_interface):
    plugins = mock_entry_points.return_value
    functions = [plugin.load.return_value for plugin in plugins]
    assert mock_interface() == 0
    for function in functions:
        function.assert_called_once_with(*mock_interface.args)
        function.reset_mock()


def test_interface_force_and_stress():
    constants = vasp.ConstantsForceAndStress(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        forces=np.random.random((3, 5)),
        stress=np.random.random((3, 3)),
    )
    with mock_vasp_plugin("force_and_stress") as force_and_stress:
        assert vasp.interface_force_and_stress(constants, None) == vasp.SUCCESS
    force_and_stress.assert_called_once_with(constants, None)
    assert all(constants.ion_types == [0, 1, 2, 2, 2])  # output: Python style indexing
    array_names = """shape_grid,ion_types,atomic_numbers,lattice_vectors,positions,\
ZVAL,POMASS,forces,stress,charge_density"""
    check_arrays_not_writeable(constants, array_names)


def test_interface_with_neighbor_list():
    constants = vasp.ConstantsForceAndStress(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        forces=np.random.random((3, 5)),
        stress=np.random.random((3, 3)),
        neighbor_list=mock_SrTiO3_neighbor_list(),
    )
    with mock_vasp_plugin("force_and_stress") as force_and_stress:
        assert vasp.interface_force_and_stress(constants, None) == vasp.SUCCESS
    force_and_stress.assert_called_once_with(constants, None)
    original_list = mock_SrTiO3_neighbor_list()
    assert len(constants.neighbor_list) == len(original_list)
    array_names = "neighbors,distances,directions"
    for original, adjusted in zip(original_list, constants.neighbor_list):
        assert all(
            original.neighbors == adjusted.neighbors + 1
        )  # output: Python style indexing
        assert np.allclose(original.distances, adjusted.distances)
        assert np.allclose(original.directions, adjusted.directions)
        check_arrays_not_writeable(adjusted, array_names)


def test_interface_machine_learning():
    vasp._machine_learning.get_calculators.cache_clear()
    constants = mock_force_and_stress_constants()
    additions = vasp.AdditionsForceAndStress(
        total_energy=2.0,
        forces=np.full((3, 5), 3.0),
        stress=np.full((3, 3), 4.0),
    )
    with mock_vasp_plugin("machine_learning") as machine_learning:
        machine_learning.return_value = MockCalculator(constants.forces.shape)
        assert vasp.interface_machine_learning(constants, additions) == vasp.SUCCESS
        assert all(
            constants.ion_types == [0, 1, 2, 2, 2]
        )  # output: Python style indexing
        array_names = """shape_grid,ion_types,atomic_numbers,lattice_vectors,positions,\
ZVAL,POMASS,forces,stress,charge_density"""
        check_arrays_not_writeable(constants, array_names)
        assert np.allclose(additions.total_energy, 3.0)
        assert np.allclose(additions.forces, 4.0)
        assert np.allclose(additions.stress, 5.0)
        # check calling the interface for a second time does not create a new calculator
        constants = mock_force_and_stress_constants()
        assert vasp.interface_machine_learning(constants, additions) == vasp.SUCCESS
        machine_learning.assert_called_once()
        check_arrays_not_writeable(constants, array_names)
        assert np.allclose(additions.total_energy, 4.0)
        assert np.allclose(additions.forces, 5.0)
        assert np.allclose(additions.stress, 6.0)


def test_entry_points_machine_learning(mock_entry_points):
    vasp._machine_learning.get_calculators.cache_clear()
    constants = mock_force_and_stress_constants()
    additions = vasp.AdditionsForceAndStress(
        total_energy=2.0,
        forces=np.full((3, 5), 3.0),
        stress=np.full((3, 3), 4.0),
    )
    shape = constants.forces.shape
    first_plugin, second_plugin = mock_entry_points.return_value
    first_interface = first_plugin.load.return_value
    first_interface.return_value = MockCalculator(shape)
    second_interface = second_plugin.load.return_value
    second_interface.return_value = MockCalculator(shape)
    assert vasp.interface_machine_learning(constants, additions) == vasp.SUCCESS
    first_interface.assert_called_once()
    second_interface.assert_called_once()
    assert np.allclose(additions.total_energy, 4.0)
    assert np.allclose(additions.forces, 5.0)
    assert np.allclose(additions.stress, 6.0)


def test_interface_local_potential():
    constants = vasp.ConstantsLocalPotential(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
    )
    with mock_vasp_plugin("local_potential") as local_potential:
        assert vasp.interface_local_potential(constants, None) == vasp.SUCCESS
    local_potential.assert_called_once_with(constants, None)
    array_names = """shape_grid,ion_types,atomic_numbers,lattice_vectors,positions,\
ZVAL,charge_density"""
    check_arrays_not_writeable(constants, array_names)


def test_interface_structure():
    shape_grid = [60, 100, 80]
    constants = vasp.ConstantsStructure(
        **mock_SrTiO3_structure(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        total_energy=np.random.random(),
        forces=np.random.random((3, 5)),
        stress=np.random.random((3, 3)),
        shape_grid=shape_grid,
        charge_density=np.random.random(shape_grid),
    )
    with mock_vasp_plugin("structure") as structure:
        assert vasp.interface_structure(constants, None) == vasp.SUCCESS
    structure.assert_called_once_with(constants, None)
    array_names = """ion_types,atomic_numbers,lattice_vectors,positions,POMASS,\
forces,stress"""
    check_arrays_not_writeable(constants, array_names)


def test_interface_occupancies():
    constants = vasp.ConstantsOccupancies(
        NELECT=27.0, EFERMI=2.5, NUPDOWN=0.5, ISMEAR=1, SIGMA=0.1, EMIN=-5.2, EMAX=6.3
    )
    with mock_vasp_plugin("occupancies") as occupancies:
        assert vasp.interface_occupancies(constants, None) == vasp.SUCCESS
    occupancies.assert_called_once_with(constants, None)


def test_error_in_plugin(mock_interface, function_name):
    def broken_interface(*args):
        raise NotImplementedError("There is no implementation here.")

    with mock_vasp_plugin(function_name, function=broken_interface):
        assert mock_interface() == vasp.ERROR_IN_PLUGIN


def mock_force_and_stress_constants():
    return vasp.ConstantsForceAndStress(
        **mock_SrTiO3_structure(),
        **mock_electron_data(),
        POMASS=np.array([87.62, 47.88, 16.0]),
        forces=np.zeros((3, 5)),
        stress=np.zeros((3, 3)),
    )


def mock_SrTiO3_structure():
    return {
        "number_ions": 5,
        "number_ion_types": 3,
        "ion_types": np.array([1, 2, 3, 3, 3]),  # input: Fortran style indexing
        "atomic_numbers": np.array([38, 22, 8, 8, 8]),
        "lattice_vectors": 4.0 * np.eye(3),
        "positions": np.array(
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ]
        ),
    }


def mock_electron_data():
    shape = (10, 12, 14)
    return {
        "NELECT": 40.0,
        "ENCUT": 250.0,
        "ZVAL": np.array([10.0, 12.0, 6.0]),
        "shape_grid": np.array(shape),
        "charge_density": np.random.random(shape),
    }


def mock_SrTiO3_neighbor_list():
    neighbors_Sr = vasp.Neighbors(
        neighbors=np.zeros(0, np.int32),
        distances=np.zeros(0),
        directions=np.zeros((0, 3)),
    )
    neighbors_Ti = vasp.Neighbors(
        neighbors=np.array([3, 4, 5, 5, 4, 3]),
        distances=2 * np.ones(6),
        directions=np.concatenate((np.eye(3), -np.eye(3))),
    )
    neighbors_O1 = vasp.Neighbors(
        neighbors=np.array([2, 2]),
        distances=np.array([2.0, 2.0]),
        directions=np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]),
    )
    neighbors_O2 = vasp.Neighbors(
        neighbors=np.array([2, 2]),
        distances=np.array([2.0, 2.0]),
        directions=np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]),
    )
    neighbors_O3 = vasp.Neighbors(
        neighbors=np.array([2, 2]),
        distances=np.array([2.0, 2.0]),
        directions=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]),
    )
    return [neighbors_Sr, neighbors_Ti, neighbors_O1, neighbors_O2, neighbors_O3]


def check_arrays_not_writeable(constants, array_names):
    for array_name in array_names.split(","):
        array = getattr(constants, array_name)
        assert not array.flags.writeable
