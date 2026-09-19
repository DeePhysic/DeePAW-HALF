import functools
import time
import traceback
from functools import wraps

import numpy as np
import vasp
import vasp._entry_points as entry_points
from ase import Atoms
from ase.calculators.calculator import Calculator, PropertyNotImplementedError
from vasp._adjust_dataclass import make_necessary_adjustments
from vasp._force_and_stress import *

__all__ = ["interface_machine_learning"]


def monitor_time(function):
    @wraps(function)
    def time_function(*args, **kwargs):
        t1 = time.perf_counter(), time.process_time()
        value = function(*args, **kwargs)
        t2 = time.perf_counter(), time.process_time()
        print(f"Real time: {t2[0] - t1[0]:.2f} seconds")
        print(f"CPU time: {t2[1] - t1[1]:.2f} seconds")
        return value

    return time_function


def _make_additions_through_calc(
    calculator: Calculator,
    constants: ConstantsForceAndStress,
    additions: AdditionsForceAndStress,
):
    atoms = Atoms(cell=constants.lattice_vectors, scaled_positions=constants.positions)
    atoms.set_atomic_numbers(constants.atomic_numbers[constants.ion_types])
    atoms.set_pbc(True)
    atoms.calc = calculator
    # In case the machine learning calculation runs on float32, then we need to
    # explicitly upgrade the arrays here in order to make them work the plugins.cpp
    additions.total_energy += np.array(atoms.get_potential_energy(), dtype=np.float64)
    additions.forces += np.array(atoms.get_forces(), dtype=np.float64)
    try:
        stress = atoms.get_stress(voigt=False)
    except PropertyNotImplementedError:
        stress = np.zeros((3, 3))
    additions.stress += np.array(stress, dtype=np.float64)


def interface_machine_learning(
    constants: ConstantsForceAndStress, additions: AdditionsForceAndStress
) -> int:
    make_necessary_adjustments(constants)
    try:
        for calculator in get_calculators():
            _make_additions_through_calc(calculator, constants, additions)
        return vasp.SUCCESS
    except:
        traceback.print_exc()
        return vasp.ERROR_IN_PLUGIN


@functools.cache
def get_calculators():
    return [
        machine_learning() for machine_learning in entry_points.load("machine_learning")
    ]
