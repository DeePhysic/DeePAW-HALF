from dataclasses import dataclass, field
from typing import List, Optional

from vasp._apply_interface import apply_interface
from vasp._neighbors import Neighbors
from vasp.typing import DoubleArray, IndexArray, IntArray

__all__ = ("ConstantsStructure", "AdditionsStructure", "interface_structure")


@dataclass(frozen=True)
class ConstantsStructure:
    number_ions: int
    number_ion_types: int
    ion_types: IndexArray
    atomic_numbers: IntArray
    lattice_vectors: DoubleArray
    positions: DoubleArray
    POMASS: DoubleArray
    total_energy: float
    forces: DoubleArray
    stress: DoubleArray
    shape_grid: IntArray
    charge_density: DoubleArray
    neighbor_list: List[Neighbors] = field(default_factory=list)


@dataclass
class AdditionsStructure:
    lattice_vectors: DoubleArray
    positions: DoubleArray


def interface_structure(
    constants: ConstantsStructure, additions: AdditionsStructure
) -> int:
    return apply_interface(constants, additions, "structure")
