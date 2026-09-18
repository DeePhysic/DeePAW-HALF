#!/usr/bin/env python3
"""Run once the Python HAPPY slice implemented by half-inspect."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("chgcar", type=Path)
    parser.add_argument("--happy-src", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.happy_src.resolve()))

    from happy.basis import PlaneWaveBasis  # noqa: PLC0415
    from happy.chgcar import read_chgcar  # noqa: PLC0415

    crystal, grid, shape = read_chgcar(str(args.chgcar))
    basis = PlaneWaveBasis.from_crystal(crystal, 400.0, *shape)
    print(f"gamma_plane_waves: {basis.npw}")
    print(f"smooth_electrons: {float(np.mean(grid)):.16e}")
    print(f"kinetic_checksum_eV: {float(np.sum(basis.kinetic)):.16e}")


if __name__ == "__main__":
    main()
