"""Finite-difference audit of the experimental fixed-charge onsite force.

The density grid is copied byte-for-byte; only one ion's CHGCAR position is
displaced.  This tests the derivative of the *same* frozen-density functional.
"""

import argparse
import json
import os
import subprocess
from pathlib import Path

def displaced_charge(source: Path, target: Path, shift_angstrom: float) -> None:
    with source.open("r") as stream:
        header = [next(stream) for _ in range(9)]
        remaining = stream.read()
    lattice = [[float(x) for x in line.split()] for line in header[2:5]]
    if not header[7].strip().lower().startswith("direct"):
        raise ValueError("This audit requires Direct CHGCAR coordinates")
    fractional = [float(x) for x in header[8].split()[:3]]
    a, b, c = lattice
    determinant = (a[0] * (b[1] * c[2] - b[2] * c[1])
                   - a[1] * (b[0] * c[2] - b[2] * c[0])
                   + a[2] * (b[0] * c[1] - b[1] * c[0]))
    if abs(determinant) < 1e-12:
        raise ValueError("Singular CHGCAR lattice")
    # The rows of lattice are the cell vectors.  For an x-Cartesian shift,
    # delta_f is the first row of the inverse lattice.
    inverse_first_row = [(b[1] * c[2] - b[2] * c[1]) / determinant,
                         (a[2] * c[1] - a[1] * c[2]) / determinant,
                         (a[1] * b[2] - a[2] * b[1]) / determinant]
    fractional = [value + shift_angstrom * inverse_first_row[i] for i, value in enumerate(fractional)]
    header[8] = "  " + "  ".join(f"{x:.16f}" for x in fractional) + "\n"
    with target.open("w") as stream:
        stream.writelines(header)
        stream.write(remaining)


def run_half(binary: Path, charge: Path, potcar: Path, output: Path, encut: float, kspacing: float,
             bands: int, onsite_lmax: int, backend: str = "cpu") -> dict:
    command = [
        str(binary), "energy", str(charge), str(potcar), "--backend", backend,
        "--onsite-lmax", str(onsite_lmax), "--encut", str(encut), "--kspacing", str(kspacing),
        "--bands", str(bands), "--no-kpoint-symmetry", "--output", str(output),
    ]
    if "center" in output.stem:
        command.append("--forces")
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    subprocess.run(command, check=True, env=environment, stdout=subprocess.DEVNULL)
    return json.loads(output.read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("charge", type=Path)
    parser.add_argument("potcar", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--center-shift", type=float, default=0.02)
    parser.add_argument("--step", type=float, default=0.002)
    parser.add_argument("--encut", type=float, default=100.0)
    parser.add_argument("--kspacing", type=float, default=2.0)
    parser.add_argument("--bands", type=int, default=8)
    parser.add_argument("--onsite-lmax", type=int, default=0)
    parser.add_argument("--backend", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, displacement in (
        ("center", args.center_shift),
        ("plus", args.center_shift + args.step),
        ("minus", args.center_shift - args.step),
    ):
        charge = args.workdir / f"{name}.CHGCAR"
        displaced_charge(args.charge, charge, displacement)
        results[name] = run_half(args.binary, charge, args.potcar, args.workdir / f"{name}.json",
                                 args.encut, args.kspacing, args.bands, args.onsite_lmax, args.backend)
    force = results["center"]["forces_eV_per_Angstrom"][0][0]
    difference = -(results["plus"]["free_energy_eV"] - results["minus"]["free_energy_eV"]) / (2 * args.step)
    report = {"analytic_force_x_eV_per_A": force, "finite_difference_force_x_eV_per_A": difference,
              "absolute_error_eV_per_A": abs(force - difference), "step_A": args.step}
    (args.workdir / "force_fd.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
