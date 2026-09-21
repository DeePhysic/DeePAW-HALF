#!/usr/bin/env python3
"""Displace one CHGCAR header coordinate while leaving its density untouched."""

from __future__ import annotations

import argparse
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--atom", type=int, required=True, help="one-based atom index")
    parser.add_argument("--axis", type=int, choices=(1, 2, 3), required=True)
    parser.add_argument("--delta", type=float, required=True, help="Cartesian displacement in Angstrom")
    args = parser.parse_args()

    lines = args.input.read_text().splitlines(keepends=True)
    scale = float(lines[1].split()[0])
    lattice = [[float(x) * scale for x in lines[i].split()[:3]] for i in range(2, 5)]
    cursor = 5
    try:
        counts = [int(x) for x in lines[cursor].split()]
    except ValueError:
        cursor += 1
        counts = [int(x) for x in lines[cursor].split()]
    cursor += 1
    if lines[cursor].lstrip().lower().startswith("s"):
        cursor += 1
    direct = lines[cursor].lstrip().lower().startswith("d")
    cursor += 1
    if not 1 <= args.atom <= sum(counts):
        raise SystemExit("atom index is outside the CHGCAR structure")
    row = cursor + args.atom - 1
    fields = lines[row].split()
    xyz = [float(x) for x in fields[:3]]
    if direct:
        inverse = inverse3(lattice)
        fractional_delta = [args.delta * inverse[args.axis - 1][j] for j in range(3)]
        xyz = [(xyz[j] + fractional_delta[j]) % 1.0 for j in range(3)]
    else:
        xyz[args.axis - 1] += args.delta / scale
    newline = "  ".join(f"{x:20.16f}" for x in xyz)
    if len(fields) > 3:
        newline += "  " + "  ".join(fields[3:])
    lines[row] = newline + "\n"
    args.output.write_text("".join(lines))


def inverse3(matrix: list[list[float]]) -> list[list[float]]:
    a, b, c = matrix
    determinant = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    if abs(determinant) < 1.0e-15:
        raise SystemExit("singular CHGCAR lattice")
    return [
        [(b[1] * c[2] - b[2] * c[1]) / determinant,
         (a[2] * c[1] - a[1] * c[2]) / determinant,
         (a[1] * b[2] - a[2] * b[1]) / determinant],
        [(b[2] * c[0] - b[0] * c[2]) / determinant,
         (a[0] * c[2] - a[2] * c[0]) / determinant,
         (a[2] * b[0] - a[0] * b[2]) / determinant],
        [(b[0] * c[1] - b[1] * c[0]) / determinant,
         (a[1] * c[0] - a[0] * c[1]) / determinant,
         (a[0] * b[1] - a[1] * b[0]) / determinant],
    ]


if __name__ == "__main__":
    main()
