# HfO2 VASP/HALF `NWAVE` boundary mismatch

## Summary

The MP-85 `mp-352 / HfO2` case failed at k point 34 with `ENCUT=520 eV`:

```text
HALF: host and HALF plane-wave bases have different sizes
DeePAW-HALF k-point solve failed on the root rank
```

Both codes use `HSQDTM |G+k|^2 < ENCUT`, matching constants, and the same
strict less-than test. The mismatch is caused by lattice precision, not a
different cutoff formula.

## Reproduction and evidence

- Case: `/data/limusen/deepawchgs/054_mp-352_HfO2`
- `KSPACING=0.35`, `KGAMMA=.TRUE.`, `EDIFF=1E-4`
- `ISPIN=1`, `LMAXMIX=-1`, `ENCUT=520 eV`
- k point 34: `(0.0, 0.5, 0.5)`
- VASP basis: 3672 plane waves
- HALF basis: 3676 plane waves

VASP uses the roughly 16-digit POSCAR lattice, whereas the DeepAW CHGCAR
header stores only six decimal places. Re-enumeration reproduces both counts.
The rounded CHGCAR lattice adds `(-1,9,-1)`, `(-1,-10,-1)`, `(1,9,0)`, and
`(1,-10,0)`, each evaluated at `519.9999852228437 eV`.

The adapter sends VASP's full-precision lattice through
`half_set_request_geometry()`, but `context_make_basis()` still builds from
`self%crystal`, parsed from CHGCAR. The retained request geometry therefore
does not yet control the file-backed solve.

## Workaround validation

The isolated retry is in
`/data/limusen/deepawchgs_retry_encut521/054_mp-352_HfO2`. At 521 eV both
lattices produce 3680 plane waves at k point 34. HALF initialized all 36 k
points, accepted the DeepAW CHGCAR, converged to `EDIFF=1E-4` in six SCF loops,
and terminated normally in 132 seconds.

## Required permanent fix

Changing ENCUT is diagnostic only. The mapped solve must treat the host
request geometry and G-vector set as authoritative. Basis generation, PAW
positions, projectors, and local-potential assembly must consistently use the
request geometry. Diagnostics should report both basis sizes and the first
differing G vector. A regression test must obtain VASP's 3672 vectors at
520 eV without changing ENCUT.

