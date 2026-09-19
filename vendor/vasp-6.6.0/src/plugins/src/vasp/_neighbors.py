from dataclasses import dataclass

from vasp.typing import DoubleArray, IndexArray

__all__ = ("Neighbors",)


@dataclass(frozen=True)
class Neighbors:
    neighbors: IndexArray
    distances: DoubleArray
    directions: DoubleArray
