from dataclasses import dataclass, field
from typing import List, Optional

from vasp._apply_interface import apply_interface
from vasp._neighbors import Neighbors
from vasp.typing import DoubleArray, IndexArray, IntArray

__all__ = (
    "ConstantsLocalPotential",
    "AdditionsLocalPotential",
    "interface_local_potential",
)


@dataclass(frozen=True)
class ConstantsLocalPotential:
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
    charge_density: Optional[DoubleArray] = None
    hartree_potential: Optional[DoubleArray] = None
    ion_potential: Optional[DoubleArray] = None
    dipole_moment: Optional[DoubleArray] = None
    neighbor_list: List[Neighbors] = field(default_factory=list)


@dataclass
class AdditionsLocalPotential:
    total_energy: float
    total_potential: DoubleArray


def interface_local_potential(
    constants: ConstantsLocalPotential, additions: AdditionsLocalPotential
) -> int:
    return apply_interface(constants, additions, "local_potential")
