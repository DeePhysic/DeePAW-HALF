#!/usr/bin/env python3
"""Compare a HALF basis CSV with HAPPY's basis for the same CHGCAR."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("half_csv", type=Path)
    parser.add_argument("chgcar", type=Path)
    parser.add_argument("--encut", type=float, default=400.0)
    parser.add_argument("--happy-src", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(args.happy_src.resolve()))
    from happy.basis import PlaneWaveBasis  # noqa: PLC0415
    from happy.chgcar import read_chgcar  # noqa: PLC0415

    crystal, _, shape = read_chgcar(str(args.chgcar))
    reference = PlaneWaveBasis.from_crystal(crystal, args.encut, *shape)
    actual = np.genfromtxt(args.half_csv, delimiter=",", names=True)

    half_g = np.column_stack([actual["n1"], actual["n2"], actual["n3"]]).astype(int)
    half_q = np.column_stack([actual["qx"], actual["qy"], actual["qz"]])
    half_fft = actual["fft_index_1based"].astype(np.int64) - 1
    half_kinetic = actual["kinetic_eV"]

    report = {
        "half_npw": int(len(actual)),
        "happy_npw": int(reference.npw),
        "gvectors_exact": bool(np.array_equal(half_g, reference.gvectors)),
        "fft_indices_exact": bool(np.array_equal(half_fft, reference.fft_index)),
        "qvector_max_abs_error_1_per_A": float(
            np.max(np.abs(half_q - reference.gvectors_cart))
        ),
        "kinetic_max_abs_error_eV": float(
            np.max(np.abs(half_kinetic - reference.kinetic))
        ),
    }
    report["passed"] = bool(
        report["half_npw"] == report["happy_npw"]
        and report["gvectors_exact"]
        and report["fft_indices_exact"]
        and report["qvector_max_abs_error_1_per_A"] < 1e-13
        and report["kinetic_max_abs_error_eV"] < 1e-11
    )
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
