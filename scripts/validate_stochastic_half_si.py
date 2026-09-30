#!/usr/bin/env python3
"""Validate stochastic HALF against a small actual Si HALF H/S oracle.

The existing v1 API assembles dense H/S, so this is a correctness check only.
"""

from __future__ import annotations

import argparse
import ctypes as ct
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from stochastic_half import (GeneralizedOperator, StochasticHalf, exact_reference,
                             polynomial_errors)


def read_half_matrices(library: Path, charge: Path, potential: Path, encut: float):
    lib = ct.CDLL(str(library))
    c_i = ct.c_int
    c_i64 = ct.c_int64
    c_d = ct.c_double
    ptr = ct.c_void_p
    lib.half_create_from_files_backend.argtypes = [
        ct.c_char_p, ct.c_char_p, c_d, c_i, c_i, c_i, c_i,
        ct.POINTER(c_i64), ptr, c_i,
    ]
    lib.half_get_basis_size.argtypes = [c_i64, ptr, ct.POINTER(c_i64), ptr, c_i]
    lib.half_assemble_hs.argtypes = [c_i64, ptr, c_i64, ptr, ptr, ptr, c_i]
    lib.half_apply_hs.argtypes = [c_i64, ptr, c_i, c_i64, ptr, ptr, ptr, ptr, c_i]
    lib.half_get_basis.argtypes = [c_i64, ptr, c_i64, ptr, ptr, ptr, c_i]
    lib.half_get_system_info.argtypes = [
        c_i64, ct.POINTER(c_i), ct.POINTER(c_i), ptr, ptr, ct.POINTER(c_d), ptr, c_i,
    ]
    lib.half_destroy.argtypes = [c_i64, ptr, c_i]
    error = ct.create_string_buffer(2048)
    handle = c_i64()
    def check(code: int):
        if code:
            raise RuntimeError(f"HALF API status {code}: {error.value.decode()}")
    check(lib.half_create_from_files_backend(
        str(charge).encode(), str(potential).encode(), encut,
        2, 1, 1, 1, ct.byref(handle), error, len(error),
    ))
    try:
        kpoint = np.zeros(3, dtype=np.float64)
        n = c_i64()
        check(lib.half_get_basis_size(handle, kpoint.ctypes.data, ct.byref(n), error, len(error)))
        count = n.value
        h = np.empty((count, count), dtype=np.complex128, order="F")
        s = np.empty_like(h, order="F")
        check(lib.half_assemble_hs(
            handle, kpoint.ctypes.data, count,
            h.ctypes.data, s.ctypes.data, error, len(error),
        ))
        gvec = np.empty((count, 3), dtype=np.int32)
        kinetic = np.empty(count, dtype=np.float64)
        check(lib.half_get_basis(
            handle, kpoint.ctypes.data, count,
            gvec.ctypes.data, kinetic.ctypes.data, error, len(error),
        ))
        ni, nt = c_i(), c_i()
        grid = np.empty(3, dtype=np.int32)
        lattice = np.empty(9, dtype=np.float64)
        volume = c_d()
        check(lib.half_get_system_info(
            handle, ct.byref(ni), ct.byref(nt), grid.ctypes.data,
            lattice.ctypes.data, ct.byref(volume), error, len(error),
        ))
        rng = np.random.default_rng(2)
        psi = np.asfortranarray(rng.normal(size=(count, 2)) + 1j*rng.normal(size=(count, 2)))
        hpsi = np.empty_like(psi, order="F")
        spsi = np.empty_like(psi, order="F")
        check(lib.half_apply_hs(
            handle, kpoint.ctypes.data, 2, count, psi.ctypes.data,
            hpsi.ctypes.data, spsi.ctypes.data, error, len(error),
        ))
        apply_error = {
            "h_max_abs": float(np.max(np.abs(hpsi - h @ psi))),
            "s_max_abs": float(np.max(np.abs(spsi - s @ psi))),
        }
        return h, s, gvec, grid, ni.value, apply_error
    finally:
        check(lib.half_destroy(handle, error, len(error)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--charge", type=Path, required=True)
    parser.add_argument("--potcar", type=Path, required=True)
    parser.add_argument("--encut", type=float, default=50.0)
    parser.add_argument("--kbt", type=float, default=0.5)
    parser.add_argument("--degree", type=int, default=256)
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    h, s, gvec, grid, nions, apply_error = read_half_matrices(
        args.library, args.charge, args.potcar, args.encut,
    )
    from scipy.linalg import eigh
    eigenvalues = eigh(h, s, eigvals_only=True)
    target_per_spin = 2 * nions  # Si POTCAR ZVAL=4; two identical spin channels.
    mu = brentq(
        lambda x: np.sum(1 / (1 + np.exp(np.clip((eigenvalues-x)/args.kbt, -700, 700)))) - target_per_spin,
        eigenvalues[0] - 20*args.kbt, eigenvalues[-1] + 20*args.kbt,
    )
    oracle = exact_reference(h, s, mu, args.kbt)
    pad = max(1.0, 0.05 * (eigenvalues[-1] - eigenvalues[0]))
    operator = GeneralizedOperator(
        len(h), lambda v: h @ v, lambda v: s @ v,
        eigenvalues[0] - pad, eigenvalues[-1] + pad,
        inverse_overlap=lambda v: np.linalg.solve(s, v),
    )
    exhaustive = StochasticHalf(
        operator, args.degree, len(h), probes=np.sqrt(len(h))*np.eye(len(h), dtype=complex),
    )
    stochastic = StochasticHalf(operator, args.degree, args.samples, seed=20261001)

    def pw_to_grid(vector):
        coefficients = np.zeros(tuple(grid), dtype=np.complex128)
        indices = gvec % grid[None, :]
        coefficients[indices[:, 0], indices[:, 1], indices[:, 2]] = vector
        return np.fft.ifftn(coefficients, norm="ortho").ravel()

    observed = exhaustive.observables(mu, args.kbt, real_space_map=pw_to_grid)
    random_result = stochastic.trace(mu, args.kbt)
    # The exact real-space diagonal is formed by the small-system oracle only.
    exact_density = np.zeros(int(np.prod(grid)), dtype=np.float64)
    for column in range(len(h)):
        basis_vector = np.zeros(len(h), dtype=complex)
        basis_vector[column] = 1
        mapped_kernel = pw_to_grid(oracle["density_kernel"][:, column])
        exact_density += np.real(mapped_kernel * np.conj(pw_to_grid(basis_vector)))
    density_error = observed["smooth_density"] - exact_density
    report = {
        "system": "Si primitive cell; fixed VASP charge and Si POTCAR",
        "encut_ev": args.encut,
        "kbt_ev": args.kbt,
        "atoms": nions,
        "basis_size": len(h),
        "exact_spectral_bounds_ev": [float(eigenvalues[0]), float(eigenvalues[-1])],
        "fft_grid": grid.tolist(),
        "target_electrons_per_spin": target_per_spin,
        "chemical_potential_ev": mu,
        "chebyshev_degree": args.degree,
        "random_probes": args.samples,
        "exact": {k: oracle[k] for k in ("electrons", "band_energy_ev")},
        "exhaustive_polynomial": {
            "electrons": observed["trace"]["electrons"],
            "band_energy_ev": observed["trace"]["band_energy_ev"],
        },
        "random_estimate": random_result,
        "polynomial_max_error": polynomial_errors(
            operator.lower_ev, operator.upper_ev, mu, args.kbt, args.degree,
        ),
        "matrix_free_apply_error": apply_error,
        "density_max_abs_error": float(np.max(np.abs(density_error))),
        "density_relative_l1_error": float(np.sum(np.abs(density_error)) / np.sum(np.abs(exact_density))),
        "note": "The existing HALF v1 API supplies dense H/S here. This validates the generalized Fermi polynomial against actual Si PAW operators, not large-system scaling.",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
