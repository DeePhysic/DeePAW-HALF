from dataclasses import dataclass, field
from typing import List, Optional

from vasp._apply_interface import apply_interface
from vasp._neighbors import Neighbors
from vasp.typing import DoubleArray, IndexArray, IntArray

__all__ = (
    "ConstantsForceAndStress",
    "AdditionsForceAndStress",
    "interface_force_and_stress",
)


@dataclass(frozen=True)
class ConstantsForceAndStress:
    ENCUT: float
    NELECT: float
    shape_grid: IntArray
    number_ions: int
    number_ion_types: int
    ion_types: IndexArray
    atomic_numbers: IntArray
    lattice_vectors: DoubleArray
    positions: DoubleArray
    ZVAL: DoubleArray
    POMASS: DoubleArray
    forces: DoubleArray
    stress: DoubleArray
    charge_density: DoubleArray
    neighbor_list: List[Neighbors] = field(default_factory=list)


@dataclass
class AdditionsForceAndStress:
    total_energy: float
    forces: DoubleArray
    stress: DoubleArray


def interface_force_and_stress(
    constants: ConstantsForceAndStress, additions: AdditionsForceAndStress
) -> int:
    return apply_interface(constants, additions, "force_and_stress")
