#!/usr/bin/env python3
"""Reconstruct a dense eSCN grid from shifted sparse probe sublattices."""

from __future__ import annotations

import argparse
import base64
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests
from ase.io import read


def decode(encoded: str, shape: tuple[int, int, int]) -> np.ndarray:
    values = np.frombuffer(base64.b64decode(encoded, validate=True), dtype="<f4")
    if values.size != int(np.prod(shape)):
        raise ValueError(f"response has {values.size} values, expected {shape}")
    return values.reshape(shape)


def predict(
    session: requests.Session,
    url: str,
    numbers: list[int],
    positions: np.ndarray,
    cell: np.ndarray,
    shape: tuple[int, int, int],
    include_uncertainty: bool,
    timeout: float,
) -> tuple[dict[str, np.ndarray], float]:
    payload = {
        "atoms": {
            "numbers": numbers,
            "positions": positions.tolist(),
            "cell": cell.tolist(),
            "pbc": [True, True, True],
        },
        "grid_shape": list(shape),
        "include_uncertainty": include_uncertainty,
    }
    started = time.perf_counter()
    response = session.post(f"{url.rstrip('/')}/v1/predict", json=payload, timeout=timeout)
    if not response.ok:
        detail = response.text.strip()
        raise requests.HTTPError(
            f"{response.status_code} response from eSCN service: {detail}",
            response=response,
        )
    result = response.json()
    arrays = {"density": decode(result["density_b64"], shape)}
    if include_uncertainty:
        for name in ("nu", "alpha", "beta", "risk"):
            arrays[name] = decode(result["uncertainty"][name], shape)
    return arrays, time.perf_counter() - started


def comparison(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
    delta = candidate.astype(np.float64) - reference.astype(np.float64)
    denominator = float(np.sum(np.abs(reference), dtype=np.float64))
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "mae": float(np.mean(np.abs(delta))),
        "rmse": float(np.sqrt(np.mean(delta * delta))),
        "relative_l1": float(np.sum(np.abs(delta)) / denominator),
        "reference_mean": float(np.mean(reference, dtype=np.float64)),
        "reconstructed_mean": float(np.mean(candidate, dtype=np.float64)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("structure", type=Path)
    parser.add_argument("--url", default="http://127.0.0.1:8265")
    parser.add_argument("--grid", nargs=3, type=int, required=True, metavar=("NX", "NY", "NZ"))
    parser.add_argument("--stride", nargs=3, type=int, required=True, metavar=("SX", "SY", "SZ"))
    parser.add_argument("--include-uncertainty", action="store_true")
    parser.add_argument("--no-oracle", action="store_true")
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()

    shape = tuple(args.grid)
    stride = tuple(args.stride)
    if min(shape) < 1 or min(stride) < 1:
        raise SystemExit("grid and stride entries must be positive")
    if any(n % s for n, s in zip(shape, stride)):
        raise SystemExit("each stride must divide its global grid dimension exactly")
    sparse_shape = tuple(n // s for n, s in zip(shape, stride))

    atoms = read(args.structure)
    cell = np.asarray(atoms.cell, dtype=float)
    positions = np.asarray(atoms.positions, dtype=float)
    numbers = atoms.numbers.tolist()
    fields = ("density", "nu", "alpha", "beta", "risk") if args.include_uncertainty else ("density",)
    assembled = {name: np.empty(shape, dtype=np.float32) for name in fields}
    session = requests.Session()
    tile_records = []

    for residue in itertools.product(*(range(s) for s in stride)):
        # rho_{R-delta}(x) = rho_R(x+delta): a negative rigid atom shift
        # moves the fixed coarse probes onto this dense-grid residue class.
        fractional_delta = np.array([r / n for r, n in zip(residue, shape)])
        delta = fractional_delta @ cell
        arrays, elapsed = predict(
            session, args.url, numbers, positions - delta, cell, sparse_shape,
            args.include_uncertainty, args.timeout,
        )
        target = tuple(slice(r, None, s) for r, s in zip(residue, stride))
        for name in fields:
            assembled[name][target] = arrays[name]
        tile_record = {"residue": list(residue), "elapsed_seconds": elapsed}
        tile_records.append(tile_record)
        print(json.dumps(tile_record), file=sys.stderr, flush=True)

    report = {
        "structure": str(args.structure.resolve()),
        "global_grid": list(shape),
        "stride": list(stride),
        "sparse_grid": list(sparse_shape),
        "requests": int(np.prod(stride)),
        "sparse_points_per_request": int(np.prod(sparse_shape)),
        "sparse_total_seconds": sum(x["elapsed_seconds"] for x in tile_records),
        "tiles": tile_records,
    }

    if not args.no_oracle:
        reference, elapsed = predict(
            session, args.url, numbers, positions, cell, shape,
            args.include_uncertainty, args.timeout,
        )
        report["oracle_seconds"] = elapsed
        report["comparison"] = {
            name: comparison(reference[name], assembled[name]) for name in fields
        }

    if args.output_prefix:
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.output_prefix.with_suffix(".npz"), **assembled)
        args.output_prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
