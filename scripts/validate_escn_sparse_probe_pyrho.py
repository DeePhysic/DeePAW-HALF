#!/usr/bin/env python3
"""Compare shifted-probe interlacing with pyRho Fourier interpolation."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import requests
from ase.io import read
from pyrho.utils import interpolate_fourier

from validate_escn_sparse_probe import comparison, predict


def spectral_metrics(reference: np.ndarray, candidate: np.ndarray, sparse_shape: tuple[int, ...]) -> dict[str, float]:
    ref_fft = np.fft.fftn(reference.astype(np.float64))
    cand_fft = np.fft.fftn(candidate.astype(np.float64))
    total_ref = float(np.sum(np.abs(ref_fft) ** 2))
    delta = cand_fft - ref_fft

    high_mask = np.zeros(reference.shape, dtype=bool)
    for axis, (dense_n, sparse_n) in enumerate(zip(reference.shape, sparse_shape)):
        mode = np.fft.fftfreq(dense_n) * dense_n
        # Modes strictly outside the sparse-grid Nyquist range cannot be
        # recovered by zero-padding that sparse FFT.  Keep the even-grid
        # Nyquist plane itself out of this mask because pyRho retains it.
        axis_mask = np.abs(mode) > sparse_n / 2
        reshape = [1] * reference.ndim
        reshape[axis] = dense_n
        high_mask |= axis_mask.reshape(reshape)

    high_ref = float(np.sum(np.abs(ref_fft[high_mask]) ** 2))
    high_cand = float(np.sum(np.abs(cand_fft[high_mask]) ** 2))
    return {
        "relative_spectral_l2": float(np.sqrt(np.sum(np.abs(delta) ** 2) / total_ref)),
        "reference_high_frequency_power_fraction": high_ref / total_ref,
        "candidate_high_frequency_power_fraction": high_cand / total_ref,
        "retained_reference_high_frequency_power_fraction": high_cand / high_ref if high_ref else 1.0,
    }


def density_summary(array: np.ndarray) -> dict[str, float | int | bool]:
    return {
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "mean": float(np.mean(array, dtype=np.float64)),
        "negative_points": int(np.count_nonzero(array < 0)),
        "finite": bool(np.isfinite(array).all()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("structure", type=Path)
    parser.add_argument("--url", default="http://127.0.0.1:8265")
    parser.add_argument("--grid", nargs=3, type=int, required=True, metavar=("NX", "NY", "NZ"))
    parser.add_argument("--stride", nargs=3, type=int, required=True, metavar=("SX", "SY", "SZ"))
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    shape = tuple(args.grid)
    stride = tuple(args.stride)
    if min(shape) < 1 or min(stride) < 1:
        raise SystemExit("grid and stride entries must be positive")
    if any(n % s for n, s in zip(shape, stride)):
        raise SystemExit("each stride must divide its global grid dimension exactly")
    sparse_shape = tuple(n // s for n, s in zip(shape, stride))

    atoms = read(args.structure)
    positions = np.asarray(atoms.positions, dtype=float)
    cell = np.asarray(atoms.cell, dtype=float)
    numbers = atoms.numbers.tolist()
    session = requests.Session()

    reference_arrays, reference_seconds = predict(
        session, args.url, numbers, positions, cell, shape, False, args.timeout
    )
    reference = reference_arrays["density"]
    direct = np.empty(shape, dtype=np.float32)
    fft_candidates: list[np.ndarray] = []
    residue_records = []

    for residue in itertools.product(*(range(s) for s in stride)):
        delta = np.array([r / n for r, n in zip(residue, shape)]) @ cell
        arrays, elapsed = predict(
            session,
            args.url,
            numbers,
            positions - delta,
            cell,
            sparse_shape,
            False,
            args.timeout,
        )
        sparse = arrays["density"]
        target = tuple(slice(r, None, s) for r, s in zip(residue, stride))
        direct[target] = sparse

        # interpolate_fourier treats the first sparse sample as fractional origin
        # zero.  The physical origin of this residue is r/N, so roll the dense
        # interpolation by +r after upsampling.
        interpolated = interpolate_fourier(sparse, list(shape))
        aligned = np.roll(interpolated, shift=residue, axis=tuple(range(len(shape))))
        sample_error = float(np.max(np.abs(aligned[target] - sparse)))
        fft_candidates.append(aligned)
        residue_records.append(
            {
                "residue": list(residue),
                "elapsed_seconds": elapsed,
                "known_sample_max_abs_error": sample_error,
                "comparison": comparison(reference, aligned),
                "spectral": spectral_metrics(reference, aligned, sparse_shape),
            }
        )

    fft_consensus = np.mean(np.stack(fft_candidates, axis=0), axis=0)
    oracle_fft_candidates = []
    for residue in itertools.product(*(range(s) for s in stride)):
        target = tuple(slice(r, None, s) for r, s in zip(residue, stride))
        oracle_sparse = reference[target]
        oracle_interpolated = interpolate_fourier(oracle_sparse, list(shape))
        oracle_fft_candidates.append(
            np.roll(oracle_interpolated, shift=residue, axis=tuple(range(len(shape))))
        )
    oracle_fft_single = oracle_fft_candidates[0]
    oracle_fft_consensus = np.mean(np.stack(oracle_fft_candidates, axis=0), axis=0)

    report = {
        "structure": str(args.structure.resolve()),
        "global_grid": list(shape),
        "stride": list(stride),
        "sparse_grid": list(sparse_shape),
        "requests": len(residue_records),
        "reference_seconds": reference_seconds,
        "sparse_total_seconds": sum(record["elapsed_seconds"] for record in residue_records),
        "methods": {
            "direct_interlacing": {
                "comparison": comparison(reference, direct),
                "spectral": spectral_metrics(reference, direct, sparse_shape),
                "density": density_summary(direct),
            },
            "pyrho_fft_consensus": {
                "definition": "mean of phase-aligned pyRho Fourier interpolants from every shifted sparse residue",
                "comparison": comparison(reference, fft_consensus),
                "spectral": spectral_metrics(reference, fft_consensus, sparse_shape),
                "density": density_summary(fft_consensus),
            },
            "oracle_downsample_pyrho_fft_single": {
                "definition": "pyRho interpolation of reference[::stride], isolating interpolation error from eSCN batching differences",
                "comparison": comparison(reference, oracle_fft_single),
                "spectral": spectral_metrics(reference, oracle_fft_single, sparse_shape),
                "density": density_summary(oracle_fft_single),
            },
            "oracle_downsample_pyrho_fft_consensus": {
                "definition": "mean of phase-aligned pyRho interpolants from every exact reference residue",
                "comparison": comparison(reference, oracle_fft_consensus),
                "spectral": spectral_metrics(reference, oracle_fft_consensus, sparse_shape),
                "density": density_summary(oracle_fft_consensus),
            },
        },
        "reference_density": density_summary(reference),
        "per_residue_pyrho_fft": residue_records,
    }

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    args.output_prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    np.savez_compressed(
        args.output_prefix.with_suffix(".npz"),
        reference=reference,
        direct_interlacing=direct,
        pyrho_fft_consensus=fft_consensus,
        oracle_downsample_pyrho_fft_single=oracle_fft_single,
        oracle_downsample_pyrho_fft_consensus=oracle_fft_consensus,
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
