#!/usr/bin/env python3
"""Time the HAPPY slice already implemented by HALF."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("chgcar", type=Path)
    parser.add_argument("--happy-src", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=5)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("--iterations must be positive")

    sys.path.insert(0, str(args.happy_src.resolve()))
    from happy.basis import PlaneWaveBasis  # noqa: PLC0415
    from happy.chgcar import read_chgcar  # noqa: PLC0415

    def run_once() -> tuple[int, float, float]:
        crystal, grid, shape = read_chgcar(str(args.chgcar))
        basis = PlaneWaveBasis.from_crystal(crystal, 400.0, *shape)
        return basis.npw, float(np.mean(grid)), float(np.sum(basis.kinetic))

    run_once()
    samples = []
    result = None
    for _ in range(args.iterations):
        start = time.perf_counter()
        result = run_once()
        samples.append(time.perf_counter() - start)
    assert result is not None
    print(json.dumps({
        "implementation": "Python HAPPY",
        "iterations": args.iterations,
        "samples_seconds": samples,
        "median_seconds": float(np.median(samples)),
        "gamma_plane_waves": result[0],
        "smooth_electrons": result[1],
        "kinetic_checksum_eV": result[2],
    }, indent=2))


if __name__ == "__main__":
    main()
