import traceback

import vasp
import vasp._entry_points as entry_points
from vasp._adjust_dataclass import make_necessary_adjustments


def apply_interface(constants, additions, interface_name: str) -> int:
    make_necessary_adjustments(constants)
    try:
        for interface in entry_points.load(interface_name):
            interface(constants, additions)
        return vasp.SUCCESS
    except:
        traceback.print_exc()
        return vasp.ERROR_IN_PLUGIN
