# Validation contract

HAPPY is the numerical oracle for HALF. Comparisons use identical CHGCAR,
POTCAR, XC, ENCUT, k points and band counts.

| Quantity | Acceptance threshold |
|---|---:|
| smooth electron count | `5e-8 electron` |
| reciprocal lattice | `1e-12 1/Angstrom` |
| G-vector integer list | exact |
| kinetic energies | `1e-11 eV` max |
| local potential grid | `1e-9 eV` max |
| overlap eigenvalues | `1e-10` max |
| reconstructed eigenvalues | `1e-7 eV` max vs HAPPY |
| PAW metric residual | `1e-10` |
| fixed-density energy terms | `1e-7 eV/cell` |

The existing Si case is the first gate. HfO2 is required before declaring
multi-species support. DeepAW-predicted densities are tested only after the
reconstructor matches HAPPY on converged VASP densities.

Performance reports must state GPU, NVHPC/CUDA versions, precision, FFT grid,
plane-wave count, bands, k points, and whether deterministic reductions were
enabled. Iteration counts and wall time are reported separately.

## Completed results

The first Si input/basis gate passed on 2026-09-18 using NVHPC 24.5.1 and an
RTX 3080. HALF and HAPPY produced identical ordered G-vector and FFT-index
arrays for all 725 Gamma-point plane waves. The maximum reciprocal-coordinate
error was `1.78e-15 1/Angstrom`; the maximum kinetic-energy error was
`1.14e-13 eV`. The machine-readable record is
[`validation/si_gamma_basis.json`](validation/si_gamma_basis.json).

The same strict basis comparison also passed with the bare-metal CUDA 13.0.88
and NVHPC 25.9 build. Its machine-readable record is
[`validation/si_gamma_basis_cuda13.json`](validation/si_gamma_basis_cuda13.json).

The later full-operator result below supersedes this slice-only milestone for
the Si Gamma case.

## Full Si Gamma PBE operator

The native Fortran path now reads the matching POTCAR, constructs NLCC,
Hartree, ionic and PBE potentials on the complete `56^3` grid, builds the PAW
projectors and DION/QPAW matrices, assembles the dense 725-by-725 H/S pair, and
solves the generalized Hermitian eigenproblem. The first eight PBE
eigenvalues agree with Python HAPPY to `1.21e-11 eV` on CPU and `1.62e-11 eV`
with the CUDA 13 cuSOLVER path. The record is
[`validation/si_gamma_full_operator.json`](validation/si_gamma_full_operator.json).

The later CUDA MIMIC_US validation below supersedes this Gamma/DION gate.
Arbitrary k points, occupations, total energy, and wavefunction output remain
separate delivery gates.

On bare-metal `cmu001`, five fresh-process PBE runs gave medians of `5.21 s`
for Python HAPPY, `0.53 s` for CPU Fortran, `0.63 s` for CUDA 13, and `0.66 s`
for CUDA 12. Thus CPU Fortran is `9.83x` and CUDA 13 is `8.27x` faster than
Python for this complete small Gamma workload. The GPU loses to CPU Fortran at
725 plane waves because startup, transfer, and small dense-solver overhead
dominate.

Port 6666 is `cmu006` with an RTX 4090 24 GB. After its driver was upgraded to
595.84, the CUDA 13.0/NVHPC 25.9 `cc89` build passed the same full Gamma
comparison. HALF 0.4's device-resident potential, projector, H/S and cuSOLVER
pipeline has a `2.45e-11 eV` maximum error. Five fresh-process runs gave a
`0.36 s` median, versus `7.06 s` Python HAPPY and `0.56 s` CPU Fortran:
CUDA 13 is `19.61x` faster than Python and `1.56x` faster than CPU. The older
hybrid CUDA path, retained in the raw historical record, had a `0.69 s` median.
The isolated CUDA 13 density kernel completes 5000 iterations in about `0.01240 s`,
`41.86x` faster than one CPU thread. Full samples and environment details are
in [`validation/full_pbe_benchmark.json`](validation/full_pbe_benchmark.json).

Port 2222 at `hpc.icqms.group` is `BH000`, a Rocky Linux 10.1 workstation with
two RTX PRO 6000 Blackwell 96 GB GPUs. The native CUDA 13.1/NVHPC 26.3 `cc120`
build passed on both devices. The full GPU PBE path has a `2.46e-11 eV`
maximum error. Five fresh-process medians are `6.74 s` for Python HAPPY and
`0.55 s` for both CPU Fortran and CUDA 13.1: both Fortran implementations are
`12.25x` faster than Python. The older hybrid CUDA path took `0.76 s`. The
isolated CUDA density kernel has a `0.011865 s` median for 5000
iterations. The complete raw samples are included in the same benchmark JSON.
An 800 eV scaling probe increases the basis from 725 to 2085 plane waves. The
full GPU internal medians are `0.564 s` on the RTX 4090 and `0.693 s` on
Blackwell. GPU H/S assembly itself is only `0.0148 s` and `0.0097 s`;
cuSOLVER takes `0.398 s` and `0.503 s`, respectively, and explains the
remaining ordering. The older hybrid path favored Blackwell because its host
CPU assembled H/S faster, illustrating why device-resident comparisons are
needed for meaningful GPU conclusions.

## HfO2 large-matrix, single-core comparison

The HfO2 fixture uses its smooth `CHGCAR`, the matching concatenated PBE
`Hf_pv`+`O` POTCAR, PBE, 520 eV, 60 bands, and Gamma point. This creates a
3407-by-3407 complex128 generalized Hermitian problem. On `hpc.icqms.group`,
both CPU executables were pinned to one logical CPU with `taskset -c 0`; all
OpenMP and BLAS thread controls were set to one. HALF CPU Fortran took
`43.94 s`; HAPPY Python with `--uspp-dij` took `64.23 s`, so HALF CPU was
`1.46x` faster.

After the CUDA QDEP implementation and the CHGCAR Fortran-to-CUDA grid-layout
fix, HALF with `--uspp-dij` reproduces all 60 recorded HAPPY HfO2 eigenvalues
to `7.9e-12 eV` maximum and `2.4e-12 eV` RMS error; the overlap extrema also
agree to roundoff. Five fresh
processes on the RTX PRO 6000 took `1.818336`, `1.695874`, `1.693474`,
`1.693574`, and `1.712614 s`, for a `1.695874 s` median. Median component
times were `0.257585 s` for the effective potential, `0.349113 s` for GPU H/S
assembly including QDEP, and `1.061438 s` for cuSOLVER. Against the matched
`64.23 s` one-core HAPPY `--uspp-dij` run, this is a `37.87x` speedup.
The machine-readable record is
[`validation/hfo2_cuda_uspp_parity.json`](validation/hfo2_cuda_uspp_parity.json).

The corrected one-core CPU Fortran MIMIC_US run took `43.790434 s`, reproducing
the 60 HAPPY bands with `7.78e-12 eV` maximum and `3.04e-12 eV` RMS error. It
is `1.47x` faster than HAPPY; CUDA is `25.82x` faster than CPU Fortran. The
earlier RTX 4090 timing used HALF's DION-only operator and therefore is not
cited as a MIMIC_US speedup.

## HALF versus Python HAPPY

For the implemented Si slice (CHGCAR read, smooth electron count, and 400 eV
Gamma basis construction), single-thread CPU Fortran HALF is `24.35x` faster
than Python HAPPY in warm in-process timing (`0.02494 s` versus `0.60715 s`).
Fresh command-line processes give `23.60x` (`0.02900 s` versus `0.68440 s`).
Both produced 725 plane waves and matching electron and kinetic checksums; the
strict vector-level comparison is recorded separately.

This earlier number remains a basis-construction microbenchmark, not the
full-operator speedup reported above. Raw samples and methodology are in
[`validation/half_vs_happy_si_slice.json`](validation/half_vs_happy_si_slice.json).

## Backend microbenchmark

The CPU and CUDA implementations run the same repeated FP64 density update and
produce identical checksums. On the Si `56^3` grid, both CUDA 12.4 and CUDA
13.0 on an RTX 3080 are about 17 times faster than one Xeon Platinum 8581C
thread. With 16 replicated grids, both are about 14 times faster than one CPU
thread. Sixty OpenMP cores are competitive on the small elementwise case and
faster on the replicated case.

CUDA 13 was compiled and timed bare metal on `cmu001`, with the mamba CUDA
13.0.3 environment and NVHPC 25.9 from NVIDIA's Rocky/RHEL RPM. Its
one-thread-comparison GPU medians differ from
CUDA 12 by only `0.11%` (small) and `0.09%` (replicated), so the versions are
effectively tied for this kernel.

These are kernel-only timings and exclude allocation and PCIe transfers. They
must not be cited as the speedup of the future PAW operator or eigensolver. The
machine-readable runs and medians are in
[`validation/backend_benchmark_cuda12.json`](validation/backend_benchmark_cuda12.json)
and
[`validation/backend_benchmark_cuda13.json`](validation/backend_benchmark_cuda13.json).
