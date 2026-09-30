#!/usr/bin/env python3
"""Benchmark eSCN latency as a function of cubic probe-grid size."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import requests
from ase.io import read

from validate_escn_sparse_probe import predict


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("structure", type=Path)
    parser.add_argument("--url", default="http://127.0.0.1:8265")
    parser.add_argument("--sizes", nargs="+", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=3000)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if min(args.sizes) < 1:
        raise SystemExit("probe sizes must be positive")

    atoms = read(args.structure)
    positions = np.asarray(atoms.positions, dtype=float)
    cell = np.asarray(atoms.cell, dtype=float)
    numbers = atoms.numbers.tolist()
    session = requests.Session()
    measurements = []

    for n in args.sizes:
        arrays, elapsed = predict(
            session,
            args.url,
            numbers,
            positions,
            cell,
            (n, n, n),
            False,
            args.timeout,
        )
        density = arrays["density"]
        record = {
            "n": n,
            "points": int(density.size),
            "batches": math.ceil(density.size / args.batch_size),
            "elapsed_seconds": elapsed,
            "points_per_second": float(density.size / elapsed),
            "density_mean": float(np.mean(density, dtype=np.float64)),
            "density_finite": bool(np.isfinite(density).all()),
        }
        measurements.append(record)
        print(json.dumps(record), flush=True)

    points = np.asarray([row["points"] for row in measurements], dtype=float)
    seconds = np.asarray([row["elapsed_seconds"] for row in measurements], dtype=float)
    design = np.column_stack((np.ones_like(points), points))
    intercept, seconds_per_probe = np.linalg.lstsq(design, seconds, rcond=None)[0]
    fitted = intercept + seconds_per_probe * points
    r_squared = 1.0 - np.sum((seconds - fitted) ** 2) / np.sum(
        (seconds - np.mean(seconds)) ** 2
    )

    report = {
        "structure": str(args.structure.resolve()),
        "atoms": len(atoms),
        "service_batch_size": args.batch_size,
        "measurements": measurements,
        "linear_fit": {
            "intercept_seconds": float(intercept),
            "seconds_per_probe": float(seconds_per_probe),
            "asymptotic_probes_per_second": float(1.0 / seconds_per_probe),
            "r_squared": float(r_squared),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
