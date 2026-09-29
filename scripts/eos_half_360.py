#!/usr/bin/env python3
"""Run a resumable DeepAW-charge -> HALF Harris E/F and band benchmark.

The source tree is read-only.  VASP OUTCAR is used only as the reference;
HALF always receives CHGCAR.deepaw, never the self-consistent CHGCAR.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np
from ase.io import read
from ase.io.vasp import read_vasp


def save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def parse_incar(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text().splitlines():
        line = line.split("!", 1)[0].split("#", 1)[0]
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip().upper()] = value.strip()
    return result


def charge_structure(path: Path, natoms: int):
    with path.open() as handle:
        lines = [next(handle) for _ in range(8 + natoms)]
    return read_vasp(io.StringIO("".join(lines)))


def material_and_kind(label: str) -> tuple[str, str, str]:
    match = re.fullmatch(r"nc-(.+?)-(mp\d+)-(.*)", label)
    if not match:
        raise ValueError(f"unexpected structure label: {label}")
    material, mpid, suffix = match.groups()
    return material, mpid, "probe" if suffix.startswith("probe") else "eos"


def onsite_lmax_for_potcar(potcar: Path, setting: str | None) -> int | None:
    """Use the highest PAW multipole supported by the actual POTCAR channels."""
    if setting is None:
        return None  # Preserve the original linear benchmark command line.
    if setting != "auto":
        return int(setting)
    channels = [int(value) for value in re.findall(
        r"(?im)^\s*non local part\s*\n\s*(\d+)\s+", potcar.read_text(errors="replace"))]
    if not channels:
        raise ValueError(f"Could not read projector angular momenta from {potcar}")
    return 2 * max(channels)


def reference(case: Path) -> dict:
    calc = case / "scf" / "calc"
    required = ["CHGCAR.deepaw", "POTCAR", "POSCAR", "OUTCAR", "INCAR", "IBZKPT"]
    for name in required:
        if not (calc / name).is_file():
            raise ValueError(f"missing {calc / name}")
    incar = parse_incar(calc / "INCAR")
    expected = {"ENCUT": 520, "KSPACING": 0.35, "SIGMA": 0.02, "ISMEAR": 0, "ISPIN": 1}
    for key, value in expected.items():
        if key not in incar or not math.isclose(float(incar[key]), value, rel_tol=0, abs_tol=1e-7):
            raise ValueError(f"{case.name}: {key} differs from the benchmark protocol: {incar.get(key)}")
    atoms = read(calc / "POSCAR", format="vasp")
    charge_atoms = charge_structure(calc / "CHGCAR.deepaw", len(atoms))
    fractional_delta = atoms.get_scaled_positions(wrap=True) - charge_atoms.get_scaled_positions(wrap=True)
    fractional_delta -= np.rint(fractional_delta)
    if (not np.array_equal(atoms.numbers, charge_atoms.numbers)
            or not np.allclose(atoms.cell.array, charge_atoms.cell.array, atol=2e-5, rtol=0)
            or np.max(np.abs(fractional_delta)) > 2e-5):
        raise ValueError(f"{case.name}: DeepAW CHGCAR geometry does not match POSCAR")
    text = (calc / "OUTCAR").read_text(errors="replace")
    matches = re.findall(r"number of bands\s+NBANDS=\s*(\d+)", text)
    if not matches:
        raise ValueError(f"{case.name}: missing VASP NBANDS")
    vasp = read(calc / "OUTCAR", index=-1)
    kpoint_lines = (calc / "IBZKPT").read_text().splitlines()
    nk = int(kpoint_lines[1])
    weights = np.array([float(line.split()[3]) for line in kpoint_lines[3:3 + nk]], dtype=float)
    weights /= weights.sum()
    material, mpid, kind = material_and_kind(case.name)
    return {
        "label": case.name, "material": material, "mpid": mpid, "kind": kind,
        "natoms": len(atoms), "volume_A3": float(atoms.get_volume()),
        "vasp_energy_eV": float(vasp.get_potential_energy()),
        "vasp_forces_eV_per_A": np.asarray(vasp.get_forces(), dtype=float).tolist(),
        "vasp_nbands": int(matches[-1]), "vasp_ir_kpoints": len(weights),
        "vasp_weight_spectrum": np.sort(weights).tolist(),
        "source": str(calc),
    }


def run_command(command: list[str], out: Path, env: dict, timeout: int) -> tuple[int, float]:
    start = time.monotonic()
    with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
        try:
            proc = subprocess.run(command, stdout=stdout, stderr=stderr, env=env, timeout=timeout, check=False)
            code = proc.returncode
        except subprocess.TimeoutExpired:
            code = 124
            stderr.write(f"\nTIMEOUT after {timeout} seconds\n")
    return code, time.monotonic() - start


def valid_energy(path: Path, natoms: int) -> bool:
    if not path.is_file():
        return False
    try:
        data = json.loads(path.read_text())
        force = np.asarray(data["forces_eV_per_Angstrom"], dtype=float)
        return (math.isfinite(float(data["internal_energy_eV"]))
                and force.shape == (natoms, 3) and np.isfinite(force).all())
    except (ValueError, KeyError, TypeError):
        return False


def reconcile_energy_convention(record: dict, result_dir: Path) -> dict:
    """Match VASP and HALF Gaussian sigma->0 energies for EOS comparison."""
    half = json.loads((result_dir / "harris.json").read_text())
    vasp_text = (Path(record["source"]) / "OUTCAR").read_text(errors="replace")
    free_matches = re.findall(r"free\s+energy\s+TOTEN\s*=\s*([-+\d.Ee]+)", vasp_text)
    internal_matches = re.findall(r"energy\s+without entropy=\s*([-+\d.Ee]+)", vasp_text)
    if not free_matches or not internal_matches:
        raise ValueError(f"{record['label']}: cannot recover VASP smearing energies")
    vasp_free = float(free_matches[-1])
    vasp_internal = float(internal_matches[-1])
    vasp_sigma0 = float(record["vasp_energy_eV"])
    if abs((vasp_internal + vasp_free) / 2 - vasp_sigma0) > 2e-5:
        raise ValueError(f"{record['label']}: VASP sigma->0 identity failed")
    half_internal = float(half["internal_energy_eV"])
    half_free = float(half["free_energy_eV"])
    record.update(vasp_internal_energy_eV=vasp_internal, vasp_free_energy_eV=vasp_free,
                  half_internal_energy_eV=half_internal, half_free_energy_eV=half_free,
                  half_energy_eV=(half_internal + half_free) / 2,
                  energy_convention="Gaussian sigma->0 = (internal + free) / 2")
    record["energy_difference_eV"] = record["half_energy_eV"] - vasp_sigma0
    return record


def run_energy(case: Path, out: Path, half: Path, env: dict,
               onsite_setting: str | None = None, onsite_tol: float | None = None,
               case_timeout: int = 3600, onsite_mix: float | None = None,
               onsite_max_iter: int | None = None) -> dict:
    ref = reference(case)
    out.mkdir(parents=True, exist_ok=True)
    if valid_energy(out / "harris.json", ref["natoms"]) and (out / "result.json").is_file():
        return json.loads((out / "result.json").read_text())
    calc = case / "scf" / "calc"
    onsite_lmax = onsite_lmax_for_potcar(calc / "POTCAR", onsite_setting)
    kcheck = subprocess.run([str(half), "kpoints", str(calc / "CHGCAR.deepaw"), "--kspacing", "0.35"],
                            capture_output=True, text=True, env=env, timeout=120, check=True)
    kgrid = json.loads(kcheck.stdout)
    half_weights = np.sort(np.array([row[3] for row in kgrid["kpoints"]], dtype=float))
    same_mesh = (len(half_weights) == ref["vasp_ir_kpoints"]
                 and np.allclose(half_weights, ref["vasp_weight_spectrum"], atol=1e-8, rtol=0))
    mode = "irreducible" if same_mesh else "full_mesh_mismatch"
    bands = ref["vasp_nbands"]
    attempts = []
    for _ in range(4):
        command = [str(half), "energy", str(calc / "CHGCAR.deepaw"), str(calc / "POTCAR"),
                   "--encut", "520", "--kspacing", "0.35", "--bands", str(bands),
                   "--ismear", "0", "--sigma", "0.02", "--xc", "pbe", "--backend", "cuda",
                   "--solver", "evd", "--forces", "--output-prefix", str(out / "harris")]
        if onsite_lmax is not None:
            command.extend(["--onsite-lmax", str(onsite_lmax)])
        if onsite_tol is not None:
            command.extend(["--onsite-tol", str(onsite_tol)])
        if onsite_mix is not None:
            command.extend(["--onsite-mix", str(onsite_mix)])
        if onsite_max_iter is not None:
            command.extend(["--onsite-max-iter", str(onsite_max_iter)])
        if mode != "irreducible":
            command.append("--no-kpoint-symmetry")
        code, seconds = run_command(command, out, env, case_timeout)
        stderr = (out / "stderr.log").read_text(errors="replace")
        attempts.append({"mode": mode, "bands": bands, "exit_code": code, "seconds": seconds,
                         "error_tail": stderr[-500:] if code else ""})
        if code == 0 and valid_energy(out / "harris.json", ref["natoms"]):
            break
        if "FFT grid is incompatible with a crystal symmetry translation" in stderr and mode == "irreducible":
            mode = "full_grid_incompatible"
        elif "highest computed band is occupied" in stderr:
            bands += 8
        else:
            raise RuntimeError(f"{case.name}: HALF energy failed; {stderr[-700:]}")
    else:
        raise RuntimeError(f"{case.name}: HALF energy failed after {len(attempts)} attempts")
    result = json.loads((out / "harris.json").read_text())
    forces = np.asarray(result["forces_eV_per_Angstrom"], dtype=float)
    vasp_forces = np.asarray(ref.pop("vasp_forces_eV_per_A"), dtype=float)
    ref.pop("vasp_weight_spectrum")
    record = {**ref, "half_energy_eV": float(result["internal_energy_eV"]),
              "half_free_energy_eV": float(result["free_energy_eV"]),
              "energy_difference_eV": float(result["internal_energy_eV"] - ref["vasp_energy_eV"]),
              "half_force_component_mae_eV_per_A": float(np.mean(np.abs(forces - vasp_forces))),
              "half_force_max_atom_error_eV_per_A": float(np.max(np.linalg.norm(forces - vasp_forces, axis=1))),
              "half_force_max_norm_eV_per_A": float(np.max(np.linalg.norm(forces, axis=1))),
              "vasp_force_max_norm_eV_per_A": float(np.max(np.linalg.norm(vasp_forces, axis=1))),
              "half_ir_kpoints": int(result["irreducible_kpoint_count"]),
              "half_full_kpoints": int(result["full_kpoint_count"]),
              "half_nbands": bands, "kpoint_mode": mode, "attempts": attempts,
              "half_binary": str(half), "onsite_lmax": onsite_lmax,
              "onsite_tol_eV": onsite_tol, "onsite_mix": onsite_mix,
              "onsite_max_iter": onsite_max_iter, "case_timeout_seconds": case_timeout}
    reconcile_energy_convention(record, out)
    save_json(out / "result.json", record)
    return record


def run_bands(case: Path, out: Path, half: Path, env: dict, npoints: int, solver: str,
              onsite_setting: str | None = None) -> dict:
    ref = reference(case)
    out.mkdir(parents=True, exist_ok=True)
    if all((out / name).is_file() for name in ("bands.json", "band_result.json", "bands_fermi.png")):
        return json.loads((out / "band_result.json").read_text())
    calc = case / "scf" / "calc"
    atoms = read(calc / "POSCAR", format="vasp")
    path = atoms.cell.bandpath(npoints=npoints)
    if len(path.kpts) < 2 or not np.isfinite(path.kpts).all():
        raise ValueError(f"{case.name}: invalid ASE band path")
    lines = [f"ASE standard path {path.path}", str(len(path.kpts)), "Reciprocal lattice"]
    lines += ["  %.12f %.12f %.12f 1" % tuple(point) for point in path.kpts]
    (out / "KPOINTS.band").write_text("\n".join(lines) + "\n")
    save_json(out / "band_path.json", {"path": path.path,
              "special_points": {key: np.asarray(val).tolist() for key, val in path.special_points.items()},
              "kpoint_count": len(path.kpts), "reference_cell_A": atoms.cell.array.tolist()})
    attempts = []
    for method in dict.fromkeys([solver, "evd"]):
        command = [str(half), "bands", str(calc / "CHGCAR.deepaw"), str(calc / "POTCAR"),
                   str(out / "KPOINTS.band"), "--encut", "520", "--bands", str(ref["vasp_nbands"]),
                   "--xc", "pbe", "--backend", "cuda", "--solver", method,
                   "--output-prefix", str(out / "bands")]
        onsite_lmax = onsite_lmax_for_potcar(calc / "POTCAR", onsite_setting)
        if onsite_lmax is not None:
            command.extend(["--onsite-lmax", str(onsite_lmax)])
        code, seconds = run_command(command, out, env, 3600)
        attempts.append({"solver": method, "exit_code": code, "seconds": seconds,
                         "error_tail": (out / "stderr.log").read_text(errors="replace")[-500:] if code else ""})
        if code == 0 and (out / "bands.json").is_file() and (out / "bands.npz").is_file():
            data = json.loads((out / "bands.json").read_text())
            plot_bands_fermi(path, out)
            result = {"label": case.name, "material": ref["material"], "mpid": ref["mpid"],
                      "path": path.path, "kpoint_count": len(path.kpts), "bands": ref["vasp_nbands"],
                      "solver": method, "attempts": attempts,
                      "sampled_band_gap_eV": data.get("sampled_band_gap_eV")}
            save_json(out / "band_result.json", result)
            return result
    raise RuntimeError(f"{case.name}: HALF bands failed: {attempts[-1]['error_tail']}")


def plot_bands_fermi(path, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bands = np.load(out / "bands.npz")
    eigenvalues = np.asarray(bands["eigenvalues_eV"], dtype=float)
    energy = json.loads((out / "harris.json").read_text())
    fermi = float(energy["fermi_level_eV"])
    x, ticks, labels = path.get_linear_kpoint_axis()
    cuts = [0] + [i for i, delta in enumerate(np.diff(x), 1) if abs(delta) < 1e-10] + [len(x)]
    fig, ax = plt.subplots(figsize=(4.7, 3.3), layout="constrained")
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        for band in (eigenvalues[lo:hi] - fermi).T:
            ax.plot(x[lo:hi], band, color="#202020", lw=0.62)
    ax.axhline(0, color="#B04A3F", lw=0.8, ls="--")
    for tick in ticks:
        ax.axvline(tick, color="#C8C8C8", lw=0.55)
    ax.set_xticks(ticks, ["Γ" if label == "G" else label for label in labels], fontsize=8)
    ax.set_xlim(float(x[0]), float(x[-1]))
    ax.set_ylim(-7, 7)
    ax.set_ylabel(r"Energy − $E_\mathrm{F}$ (eV)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "bands_fermi.png", dpi=210)
    plt.close(fig)


def summarize(root: Path, output: Path) -> None:
    records = []
    errors = []
    for case in sorted(root.glob("nc-*-mp*")):
        result = output / case.name / "result.json"
        if result.is_file():
            record = json.loads(result.read_text())
            if record.get("energy_convention") != "Gaussian sigma->0 = (internal + free) / 2":
                reconcile_energy_convention(record, output / case.name)
                save_json(result, record)
            records.append(record)
        elif (output / case.name / "error.json").is_file():
            errors.append(json.loads((output / case.name / "error.json").read_text()))
    if records:
        flat = [{k: v for k, v in row.items() if k != "attempts"} for row in records]
        fieldnames = list(dict.fromkeys(key for row in flat for key in row))
        with (output / "harris_360.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(flat)
    summary = {"source_count": len(list(root.glob("nc-*-mp*"))), "energy_force_completed": len(records),
               "energy_force_failed": len(errors),
               "bands_completed": len(list(output.glob("nc-*-mp*/band_result.json"))),
               "mean_force_component_mae_eV_per_A": float(np.mean([r["half_force_component_mae_eV_per_A"] for r in records])) if records else None,
               "probe_force_component_mae_eV_per_A": float(np.mean([r["half_force_component_mae_eV_per_A"] for r in records if r["kind"] == "probe"])) if any(r["kind"] == "probe" for r in records) else None,
               "max_force_atom_error_eV_per_A": float(np.max([r["half_force_max_atom_error_eV_per_A"] for r in records])) if records else None}
    save_json(output / "summary.json", summary)
    if len(records) < 30:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    grouped = {}
    for row in records:
        if row["kind"] == "eos":
            grouped.setdefault((row["material"], row["mpid"]), []).append(row)
    curve_rows = []
    all_relative_errors = []
    for key, rows in sorted(grouped.items()):
        if len(rows) < 3:
            continue
        rows.sort(key=lambda r: r["volume_A3"])
        volumes = np.array([r["volume_A3"] for r in rows])
        ev = np.array([r["vasp_energy_eV"] for r in rows])
        eh = np.array([r["half_energy_eV"] for r in rows])
        jv = int(np.argmin(ev))
        jh = int(np.argmin(eh))
        delta_mev_atom = ((eh - eh[jv]) - (ev - ev[jv])) * 1000 / rows[0]["natoms"]
        all_relative_errors.extend(np.abs(delta_mev_atom).tolist())
        curve_rows.append({"material": key[0], "mpid": key[1], "eos_points": len(rows),
                           "relative_energy_mae_meV_per_atom": float(np.mean(np.abs(delta_mev_atom))),
                           "relative_energy_max_error_meV_per_atom": float(np.max(np.abs(delta_mev_atom))),
                           "discrete_minimum_same_point": jv == jh,
                           "vasp_min_volume_A3": float(volumes[jv]),
                           "half_min_volume_A3": float(volumes[jh]),
                           "minimum_volume_difference_percent": float(100 * (volumes[jh] / volumes[jv] - 1))})
    if curve_rows:
        with (output / "eos_curve_errors.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(curve_rows[0]))
            writer.writeheader()
            writer.writerows(curve_rows)
    summary["eos_relative_energy_mae_meV_per_atom"] = float(np.mean(all_relative_errors)) if all_relative_errors else None
    summary["eos_materials_with_curve"] = len(curve_rows)
    save_json(output / "summary.json", summary)
    fig, axes = plt.subplots(6, 5, figsize=(16, 17), constrained_layout=True)
    for ax, (key, rows) in zip(axes.flat, sorted(grouped.items())):
        if len(rows) < 3:
            ax.set_visible(False)
            continue
        rows.sort(key=lambda r: r["volume_A3"])
        volume = np.array([r["volume_A3"] for r in rows])
        ev = np.array([r["vasp_energy_eV"] for r in rows])
        eh = np.array([r["half_energy_eV"] for r in rows])
        natoms = rows[0]["natoms"]
        reference_index = int(np.argmin(ev))
        v0 = volume[reference_index]
        x = volume / v0
        ax.plot(x, (ev - ev[reference_index]) * 1000 / natoms, "o-", color="#242424", lw=1.35, ms=3, label="VASP SCF")
        ax.plot(x, (eh - eh[reference_index]) * 1000 / natoms, "s-", color="#0072B2", lw=1.35, ms=3, label="DeePAW-HALF")
        ax.set_title(f"{key[0]}  {key[1]}", fontsize=10)
        ax.tick_params(labelsize=8, length=3)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.flat[-5:]:
        ax.set_xlabel(r"$V/V_{0,\mathrm{VASP}}$", fontsize=9)
    for row in axes:
        row[0].set_ylabel("Relative energy (meV/atom)", fontsize=9)
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(output / "eos_relative_energy_30.png", dpi=220)
    fig.savefig(output / "eos_relative_energy_30.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--half", type=Path, required=True)
    parser.add_argument("--gpu", default="0")
    parser.add_argument("--onsite-lmax", help="HALF onsite maximum: integer, or auto=2*max(POTCAR projector l); omitted preserves linear")
    parser.add_argument("--onsite-tol", type=float, help="Explicit onsite D_ij convergence tolerance in eV")
    parser.add_argument("--onsite-mix", type=float, help="Fixed onsite damping in (0,1] for recovery cases")
    parser.add_argument("--onsite-max-iter", type=int, help="Maximum onsite iterations for recovery cases")
    parser.add_argument("--case-timeout", type=int, default=3600, help="Per-energy-case wall-time limit in seconds")
    parser.add_argument("--phase", choices=["energy", "bands", "both", "summarize"], default="energy")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--only")
    parser.add_argument("--band-points", type=int, default=60)
    parser.add_argument("--band-solver", choices=["evd", "evx", "acc"], default="evx")
    args = parser.parse_args()
    if args.case_timeout <= 0:
        parser.error("--case-timeout must be positive")
    if args.onsite_mix is not None and not (0 < args.onsite_mix <= 1):
        parser.error("--onsite-mix must be in (0,1]")
    if args.onsite_max_iter is not None and args.onsite_max_iter < 2:
        parser.error("--onsite-max-iter must be at least 2")
    args.output.mkdir(parents=True, exist_ok=True)
    cases = sorted(p for p in args.source.glob("nc-*-mp*") if p.is_dir())
    if args.only:
        cases = [p for p in cases if args.only in p.name]
    if args.limit is not None:
        cases = cases[:args.limit]
    if args.phase == "summarize":
        summarize(args.source, args.output)
        return
    env = os.environ.copy()
    env.update(CUDA_VISIBLE_DEVICES=args.gpu, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for index, case in enumerate(cases, 1):
        out = args.output / case.name
        for phase in (["energy", "bands"] if args.phase == "both" else [args.phase]):
            try:
                result = run_energy(case, out, args.half, env, args.onsite_lmax,
                                    args.onsite_tol, args.case_timeout, args.onsite_mix,
                                    args.onsite_max_iter) if phase == "energy" else run_bands(
                    case, out, args.half, env, args.band_points, args.band_solver, args.onsite_lmax)
                if (out / "error.json").exists():
                    (out / "error.json").rename(out / f"error.resolved.{int(time.time())}.json")
                print(f"{index}/{len(cases)} {case.name} {phase} OK", flush=True)
                if phase == "energy":
                    print(f"  E={result['half_energy_eV']:.8f} eV; F-MAE={result['half_force_component_mae_eV_per_A']:.6f} eV/A; k={result['kpoint_mode']}", flush=True)
            except Exception as exc:
                out.mkdir(parents=True, exist_ok=True)
                save_json(out / "error.json", {"label": case.name, "phase": phase, "error": str(exc)})
                print(f"{index}/{len(cases)} {case.name} {phase} FAILED: {exc}", flush=True)
        if index % 10 == 0 or index == len(cases):
            summarize(args.source, args.output)


if __name__ == "__main__":
    main()
