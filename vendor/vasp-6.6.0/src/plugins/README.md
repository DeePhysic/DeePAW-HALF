`vasp (python)`
--------------

This Python package links VASP to Python through a C++ interface. Through a selected number of interfaces, modify interior working of VASP through Python scripting. This README provides minimal documentation to get started with this interface.

## Installation
1. Create a new [conda environment](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html). Enter that environment and navigate to the directory `src/plugins`. Here you will need to do `pip install .` to install the `vasp` python package.
2. [Compile VASP](https://www.vasp.at/wiki/index.php/Installing_VASP.6.X.X) as you would normally do. You will have to make just one adjustment to the `makefile.include` file: add the following lines at the bottom:
```bash
CPP_OPTIONS+= -DPLUGINS
LLIBS      += $(shell python3-config --ldflags --embed) -lstdc++
CXX_FLAGS   = -Wall -Wextra  $(shell python3 -m pybind11 --includes) -std=c++11
```
NOTE: Make sure to be within the conda environment when you compile VASP.
3. When running VASP with the python interface you will need to add the `lib` directory of your python to `LD_LIBRARY_PATH`. You can do this somewhat easily with conda by running `python3-config --ldflags`, followed by `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path-from-earlier-command>`

## INCAR tags
You are now ready to work with the python interface of vasp. The available options that can be set in the INCAR file are:

```bash
PLUGINS/LOCAL_POTENTIAL = T         ! Modified the local potential every SCF step
PLUGINS/FORCE_AND_STRESS = T        ! Modifies the force and stress at the end of the SCF loop
PLUGINS/STRUCTURE = T               ! Modifies the structure at the end of the SCF loop
PLUGINS/OCCUPANCIES = T             ! Modifies NELECT, EFERMI, SIGMA, ISMEAR, EMIN, EMAX, NUPDOWN at the end of the SCF loop
```

## Python scripting (getting started)
1. Create a file called `vasp_plugin.py` in the folder in which you are running VASP. An alternative route through entry-points is also possible (see below for further instructions).
2. The relevant lines of python code to access data from vasp through the any of the interfaces listed above are:
```python
def local_potential(constants, additions):
	"""Defines the PLUGINS/LOCAL_POTENTIAL interface"""

def force_and_stress(constants, additions):
	"""Defines the PLUGINS/FORCE_AND_STRESS interface"""

def structure(constants, additions):
	"""Defines the PLUGINS/STRUCTURE interface"""

def occupancies(constants, additions):
	"""Defines the PLUGINS/OCCUPANCIES interface"""
```
As the names of the input variables suggest,
- `constants` provides a dataclass with quantities computed in VASP that must stay constant through the interface. Typically these include "metadata" such as the lattice vectors, number of grid points, etc. Some interfaces also pass in quantities such as the charge density as `constants`. Please check `src/plugins/src` for for details on what quantities are available as constants. Alternatively, you may print the contents of `constants`.
- `additions` stores quantities that are allowed to change in the interface. NOTE: You must use the syntax `additions.<quantity-name> += ...` in order to change additions. That is, quantities must be "added" (or subtracted) to `additions` instead of them being assigned a new value.


## Entry points (Advanced)
Consider using [entry points](https://packaging.python.org/en/latest/specifications/entry-points/) in cases where you do not want to move `vasp_plugin.py` to each new VASP calculation directory. Within the group "vasp", add the following lines to your `pyproject.toml` file for each interface you would like to introduce new functions for,
```bash
force_and_stress = "/path/to/python_function:force_and_stress"
```
