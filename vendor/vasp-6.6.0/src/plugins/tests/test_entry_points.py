import pytest
import vasp._entry_points as entry_points

from .util import mock_vasp_plugin


def test_nonexisting_plugin(function_name):
    with pytest.warns(UserWarning):
        assert entry_points.load(function_name) == []


def test_vasp_plugin_with_function(function_name):
    with mock_vasp_plugin(function_name):
        functions = entry_points.load(function_name)
        assert len(functions) == 1
        assert functions[0].name == f"mock_{function_name}"


def test_vasp_plugin_without_function(function_name):
    with pytest.warns(UserWarning), mock_vasp_plugin(function_name, function=None):
        assert entry_points.load(function_name) == []


def test_entry_points(mock_entry_points, function_name):
    functions = mock_entry_points.return_value
    expected = [function.load.return_value for function in functions]
    assert entry_points.load(function_name) == expected
    for function in functions:
        function.load.assert_called_once_with()
        function.reset_mock()
