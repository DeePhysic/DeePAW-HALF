from dataclasses import dataclass
from typing import Optional

from vasp._apply_interface import apply_interface

__all__ = (
    "ConstantsOccupancies",
    "AdditionsOccupancies",
    "interface_occupancies",
)


@dataclass(frozen=True)
class ConstantsOccupancies:
    NELECT: float
    EFERMI: float
    NUPDOWN: float
    ISMEAR: int
    SIGMA: float
    EMIN: float
    EMAX: float


@dataclass
class AdditionsOccupancies:
    NELECT: float
    EFERMI: float
    NUPDOWN: float
    ISMEAR: int
    SIGMA: float
    EMIN: float
    EMAX: float


def interface_occupancies(
    constants: ConstantsOccupancies, additions: AdditionsOccupancies
) -> int:
    return apply_interface(constants, additions, "occupancies")
