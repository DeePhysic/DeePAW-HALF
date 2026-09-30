#!/usr/bin/env python3
"""Compare CPU/CUDA matrix-free HALF H/S against a small-cell dense oracle."""

import argparse
import json
from pathlib import Path

import numpy as np

from half_operator_api import HalfOperatorAPI
from validate_stochastic_half_si import read_half_matrices


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--charge", type=Path, required=True)
    parser.add_argument("--potcar", type=Path, required=True)
    parser.add_argument("--encut", type=float, required=True)
    parser.add_argument("--backend", choices=("cpu", "cuda"), default="cuda")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    h, s, _, _, atoms, _ = read_half_matrices(
        args.library, args.charge, args.potcar, args.encut,
    )
    rng = np.random.default_rng(3)
    vector = rng.normal(size=len(h)) + 1j*rng.normal(size=len(h))
    with HalfOperatorAPI(
        args.library, args.charge, args.potcar, args.encut,
        backend=1 if args.backend == "cpu" else 2,
    ) as api:
        hvector, svector = api.apply(vector)
    report = {
        "atoms": atoms,
        "plane_waves": len(h),
        "backend": args.backend,
        "h_psi_max_abs_error": float(np.max(np.abs(hvector - h @ vector))),
        "s_psi_max_abs_error": float(np.max(np.abs(svector - s @ vector))),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
