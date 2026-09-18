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

This closes the Gamma/DION validation gate. Potential-dependent MIMIC_US D,
arbitrary k points, occupations, total energy, and wavefunction output remain
separate delivery gates; the result is not yet a complete replacement for all
HAPPY commands.

On bare-metal `cmu001`, five fresh-process PBE runs gave medians of `5.21 s`
for Python HAPPY, `0.53 s` for CPU Fortran, `0.63 s` for CUDA 13, and `0.66 s`
for CUDA 12. Thus CPU Fortran is `9.83x` and CUDA 13 is `8.27x` faster than
Python for this complete small Gamma workload. The GPU loses to CPU Fortran at
725 plane waves because startup, transfer, and small dense-solver overhead
dominate.

Port 6666 is `cmu006` with an RTX 4090 24 GB. After its driver was upgraded to
595.84, the CUDA 13.0/NVHPC 25.9 `cc89` build passed the same full Gamma
comparison with a `1.62e-11 eV` maximum error. Five fresh-process runs gave a
`0.69 s` median, versus the existing `7.06 s` Python HAPPY and `0.56 s` CPU
Fortran medians: CUDA 13 is therefore `10.23x` faster than Python but remains
`1.23x` slower than CPU Fortran for this small 725-plane-wave workload. The
isolated CUDA 13 density kernel completes 5000 iterations in about `0.01240 s`,
`41.86x` faster than one CPU thread. Full samples and environment details are
in [`validation/full_pbe_benchmark.json`](validation/full_pbe_benchmark.json).

Port 2222 at `hpc.icqms.group` is `BH000`, a Rocky Linux 10.1 workstation with
two RTX PRO 6000 Blackwell 96 GB GPUs. The native CUDA 13.1/NVHPC 26.3 `cc120`
build passed on both devices, and GPU0/GPU1 produced identical first-eight
eigenvalues. Against HAPPY run on the same node and inputs, the maximum
eigenvalue difference is `8.31e-12 eV`. Five fresh-process medians are `6.74 s`
for Python HAPPY, `0.55 s` for CPU Fortran, and `0.76 s` for CUDA 13.1: the
Fortran implementations are respectively `12.25x` and `8.87x` faster than
Python. The isolated CUDA density kernel has a `0.011865 s` median for 5000
iterations. The complete raw samples are included in the same benchmark JSON.
An 800 eV scaling probe increases the basis from 725 to 2085 plane waves. At
that size the Blackwell node completes the internal workflow in `2.011 s`
versus `2.275 s` on the RTX 4090, reversing the small-case ordering. The
Blackwell cuSOLVER stage alone remains slower (`0.517 s` versus `0.436 s`),
but its deficit narrows from about 62% at 725 plane waves to about 19% at 2085.
This confirms that fixed overhead and small-matrix utilization dominate the
original comparison; the different host CPUs also affect H/S assembly.

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
