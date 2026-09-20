# HfO2 VASP-like NBANDS validation

Date: 2026-09-20

The `deepaw-half-vasp` branch adds `HALF_MODE=VASP_LIKE` (`VASP` was used as
its accepted short alias in this validation). VASP passes
`WDES%NB_TOT` through the existing C ABI as `nbands`; HALF then uses the CUDA
generalized-Hermitian index-range solver (`cusolverDnZhegvdx`) to compute only
the lowest requested eigenpairs. `HALF_MODE=TRADITIONAL` is the default and
computes the full spectrum (`DENSE` is its accepted short alias).

The comparison ran on an NVIDIA RTX PRO 6000 Blackwell with NVHPC 26.5 and its
CUDA 13.2 toolkit. Both runs used the same HfO2 input:

- 36 k-points, 52 bands, approximately 3,700 plane waves per k-point
- `ENCUT=521`, `KSPACING=0.35`, `EDIFF=1E-4`, `ALGO=All`
- `ISPIN=1`, `LMAXMIX=-1`, `ICHARG=1`
- DeepAW CHGCAR input; `LHALF_INIT=.TRUE.`, `LHALF_API=.FALSE.`
- one MPI rank and one CPU thread; the same GPU and executable

| HALF mode | VASP elapsed time | SCF LOOP count | Final EIGENVAL values |
|---|---:|---:|---:|
| `VASP` | 105.990 s | 6 | reference |
| `DENSE` | 135.934 s | 6 | max absolute difference 0.0 eV |

For this end-to-end run, VASP-like mode is **1.283x faster** and reduces wall
time by **22.0%**. The SCF trajectory and LOOP count are unchanged. The output
explicitly records `requested NBANDS=52 (VASP-like partial spectrum)` for all
36 k-points.

Run artifacts are stored outside the repository at
`/data/limusen/half-vasp-like-validation/{vasp,dense}` on the Pro 6000 host.
