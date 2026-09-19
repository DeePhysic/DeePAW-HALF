from typing import ClassVar, Dict, Protocol

import numpy as np
import numpy.typing as np_typing

IntArray = np_typing.NDArray[np.int32]
IndexArray = np_typing.NDArray[np.int32]
DoubleArray = np_typing.NDArray[np.float64]


class Dataclass(Protocol):
    __dataclass_fields__: ClassVar[Dict]
