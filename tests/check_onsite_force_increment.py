"""Compare the onsite force increment with a fixed-density energy derivative.

The existing Harris force has other contributions, so subtract the -1 reference
on the same displaced structures to isolate the new one-centre term.
"""

import argparse
import json
from pathlib import Path

from check_onsite_force_fd import displaced_charge, run_half


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("charge", type=Path)
    parser.add_argument("potcar", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--onsite-lmax", type=int, required=True)
    parser.add_argument("--reference-lmax", type=int, default=-1)
    parser.add_argument("--center-shift", type=float, default=0.02)
    parser.add_argument("--step", type=float, default=0.002)
    parser.add_argument("--encut", type=float, default=200.0)
    parser.add_argument("--kspacing", type=float, default=1.5)
    parser.add_argument("--bands", type=int, default=8)
    parser.add_argument("--backend", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    if args.onsite_lmax < 0:
        parser.error("--onsite-lmax must be non-negative")
    if args.reference_lmax < -1 or args.reference_lmax >= args.onsite_lmax:
        parser.error("--reference-lmax must be -1 or a lower non-negative channel")
    args.workdir.mkdir(parents=True, exist_ok=True)
    results = {}
    for lmax in (args.reference_lmax, args.onsite_lmax):
        case = {}
        for name, displacement in (
            ("center", args.center_shift),
            ("plus", args.center_shift + args.step),
            ("minus", args.center_shift - args.step),
        ):
            charge = args.workdir / f"{name}.CHGCAR"
            displaced_charge(args.charge, charge, displacement)
            output = args.workdir / f"lmax_{lmax}_{name}.json"
            case[name] = run_half(args.binary, charge, args.potcar, output,
                                  args.encut, args.kspacing, args.bands, lmax, args.backend)
        results[lmax] = case
    target, reference = results[args.onsite_lmax], results[args.reference_lmax]
    analytic = (target["center"]["forces_eV_per_Angstrom"][0][0]
                - reference["center"]["forces_eV_per_Angstrom"][0][0])
    delta_plus = target["plus"]["free_energy_eV"] - reference["plus"]["free_energy_eV"]
    delta_minus = target["minus"]["free_energy_eV"] - reference["minus"]["free_energy_eV"]
    finite_difference = -(delta_plus - delta_minus) / (2.0 * args.step)
    report = {
        "onsite_lmax": args.onsite_lmax,
        "reference_lmax": args.reference_lmax,
        "backend": args.backend,
        "analytic_increment_eV_per_A": analytic,
        "finite_difference_increment_eV_per_A": finite_difference,
        "absolute_increment_error_eV_per_A": abs(analytic - finite_difference),
        "step_A": args.step,
    }
    (args.workdir / "onsite_force_increment.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
