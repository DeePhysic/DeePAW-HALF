import sys
import types
from contextlib import contextmanager
from unittest.mock import MagicMock


@contextmanager
def mock_vasp_plugin(function_name, function="mock"):
    vasp_plugin = types.ModuleType("vasp_plugin")
    vasp_plugin.__spec__ = "mock module spec"
    if function == "mock":
        function = MagicMock()
        function.name = f"mock_{function_name}"
        setattr(vasp_plugin, function_name, function)
    elif function is not None:
        setattr(vasp_plugin, function_name, function)
    sys.modules["vasp_plugin"] = vasp_plugin
    yield function
    del sys.modules["vasp_plugin"]
