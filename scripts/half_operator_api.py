"""CPU matrix-free HALF H/S adapter for the stochastic prototype.

Requires a libhalf built from large-scalebranch.  The C ABI is unchanged; its
CPU half_apply_hs implementation now caches PAW projectors and uses FFT matvec.
"""

from __future__ import annotations

import ctypes as ct
from pathlib import Path

import numpy as np

from stochastic_half import GeneralizedOperator, estimate_spectral_bounds


class HalfOperatorAPI:
    def __init__(
        self, library: Path, charge: Path, potential: Path, encut: float,
        *, xc: int = 2, use_uspp: bool = True, backend: int = 1,
    ) -> None:
        self.lib = ct.CDLL(str(library))
        self.error = ct.create_string_buffer(2048)
        self.handle = ct.c_int64()
        self.kpoint = np.zeros(3, dtype=np.float64)
        ptr = ct.c_void_p
        self.lib.half_create_from_files_backend.argtypes = [
            ct.c_char_p, ct.c_char_p, ct.c_double, ct.c_int, ct.c_int,
            ct.c_int, ct.c_int, ct.POINTER(ct.c_int64), ptr, ct.c_int,
        ]
        self.lib.half_get_basis_size.argtypes = [
            ct.c_int64, ptr, ct.POINTER(ct.c_int64), ptr, ct.c_int,
        ]
        self.lib.half_apply_hs.argtypes = [
            ct.c_int64, ptr, ct.c_int, ct.c_int64, ptr, ptr, ptr, ptr, ct.c_int,
        ]
        self.lib.half_destroy.argtypes = [ct.c_int64, ptr, ct.c_int]
        self._check(self.lib.half_create_from_files_backend(
            str(charge).encode(), str(potential).encode(), encut,
            xc, int(use_uspp), backend, 1,
            ct.byref(self.handle), self.error, len(self.error),
        ))
        size = ct.c_int64()
        try:
            self._check(self.lib.half_get_basis_size(
                self.handle, self.kpoint.ctypes.data,
                ct.byref(size), self.error, len(self.error),
            ))
            self.size = size.value
        except BaseException:
            self.close()
            raise

    def _check(self, status: int) -> None:
        if status:
            raise RuntimeError(f"HALF status {status}: {self.error.value.decode()}")

    def apply(self, vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(vector, dtype=np.complex128)
        if values.shape != (self.size,):
            raise ValueError(f"expected one complex vector of length {self.size}")
        states = np.asfortranarray(values.reshape(self.size, 1))
        hstates = np.empty_like(states, order="F")
        sstates = np.empty_like(states, order="F")
        self._check(self.lib.half_apply_hs(
            self.handle, self.kpoint.ctypes.data, 1, self.size,
            states.ctypes.data, hstates.ctypes.data, sstates.ctypes.data,
            self.error, len(self.error),
        ))
        return hstates[:, 0].copy(), sstates[:, 0].copy()

    def generalized(self, lower_ev: float, upper_ev: float) -> GeneralizedOperator:
        return GeneralizedOperator(
            self.size, lambda v: self.apply(v)[0], lambda v: self.apply(v)[1],
            lower_ev, upper_ev,
        )

    def estimated_generalized(self, padding_ev: float = 2.0) -> GeneralizedOperator:
        overlap_solver = self.generalized(-1, 1)
        bounds = estimate_spectral_bounds(
            self.size, lambda v: self.apply(v)[0],
            lambda v: self.apply(v)[1],
            overlap_solver.solve_s,
            padding_ev=padding_ev,
        )
        return self.generalized(*bounds)

    def close(self) -> None:
        if self.handle.value:
            self._check(self.lib.half_destroy(self.handle, self.error, len(self.error)))
            self.handle.value = 0

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
