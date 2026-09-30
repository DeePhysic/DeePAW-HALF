#!/usr/bin/env python3
"""Run stochastic finite-T HALF traces through libhalf's matrix-free CPU API."""

import argparse
import json
from pathlib import Path

from half_operator_api import HalfOperatorAPI
from stochastic_half import StochasticHalf, polynomial_errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--charge", type=Path, required=True)
    parser.add_argument("--potcar", type=Path, required=True)
    parser.add_argument("--encut", type=float, required=True)
    parser.add_argument("--backend", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--electrons-per-spin", type=float, required=True)
    parser.add_argument("--kbt", type=float, default=0.5)
    parser.add_argument("--degree", type=int, default=256)
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--lower-ev", type=float)
    parser.add_argument("--upper-ev", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.lower_ev is None) != (args.upper_ev is None):
        parser.error("supply both --lower-ev and --upper-ev, or neither")
    with HalfOperatorAPI(args.library, args.charge, args.potcar, args.encut,
                         backend=1 if args.backend == "cpu" else 2) as api:
        operator = (api.estimated_generalized() if args.lower_ev is None else
                    api.generalized(args.lower_ev, args.upper_ev))
        solver = StochasticHalf(operator, args.degree, args.samples, args.seed)
        mu = solver.find_mu(args.electrons_per_spin, args.kbt)
        report = {
            "basis_size": api.size,
            "encut_ev": args.encut,
            "backend": args.backend,
            "spectral_bounds_ev": [operator.lower_ev, operator.upper_ev],
            "kbt_ev": args.kbt,
            "degree": args.degree,
            "samples": args.samples,
            "seed": args.seed,
            "target_electrons_per_spin": args.electrons_per_spin,
            "mu_ev": mu,
            "trace": solver.trace(mu, args.kbt),
            "polynomial_max_error": polynomial_errors(
                operator.lower_ev, operator.upper_ev, mu, args.kbt, args.degree,
            ),
            "note": "Finite-temperature fixed-density band trace, not a PAW total energy",
        }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
