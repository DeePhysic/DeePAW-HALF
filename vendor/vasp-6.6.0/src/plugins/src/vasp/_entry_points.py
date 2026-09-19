import warnings
from importlib import import_module
from importlib.metadata import entry_points
from importlib.util import find_spec


def load(function_name):
    functions = _load_entry_points(function_name) + _load_vasp_plugin(function_name)
    _raise_warning_if_no_functions(function_name, functions)
    return functions


def _load_entry_points(function_name):
    return [
        entry_point.load()
        for entry_point in entry_points(name=function_name, group="vasp")
    ]


def _load_vasp_plugin(function_name):
    if not find_spec("vasp_plugin"):
        return []
    vasp_plugin = import_module("vasp_plugin")
    function = getattr(vasp_plugin, function_name, None)
    if function:
        return [function]
    else:
        return []


def _raise_warning_if_no_functions(function_name, functions):
    if functions:
        return
    message = f"""\
No plugin for {function_name} found. Please make sure that you have either installed \
Python packages that provide entry points for VASP, installed a package called \
vasp_plugin, or have a file called vasp_plugin.py in the directory in which you have \
launched your VASP calculation."""
    warnings.warn(message)
