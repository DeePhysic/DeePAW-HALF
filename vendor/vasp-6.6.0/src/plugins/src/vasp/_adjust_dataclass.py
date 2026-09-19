from dataclasses import fields

import numpy as np
from vasp.typing import Dataclass, IndexArray


def make_necessary_adjustments(constants: Dataclass):
    adjust_indexing(constants)
    freeze_arrays(constants)
    adjust_neighbor_list(constants)


def adjust_indexing(constants: Dataclass):
    """Fortran is 1-indexed while Python is 0-indexed. Adjust the index such that it
    it follows the right indexing for Python."""
    for field in fields(constants):
        if field.type is IndexArray:
            quantity = getattr(constants, field.name)
            quantity -= 1


def freeze_arrays(constants: Dataclass):
    """We make all arrays in the dataclass immutable. This is not strictly true because
    the user can always undo this but it should at least prevent accidental modification.
    """
    for field in fields(constants):
        possible_array = getattr(constants, field.name)
        if isinstance(possible_array, np.ndarray):
            possible_array.flags.writeable = False


def adjust_neighbor_list(constants: Dataclass):
    """The neighbor list is an optional argument present on classes with structural data
    if it is present, we need to adjust all dataclasses in the list."""
    neighbor_list = getattr(constants, "neighbor_list", None)
    if neighbor_list is None:
        return
    for neighbor in neighbor_list:
        make_necessary_adjustments(neighbor)
