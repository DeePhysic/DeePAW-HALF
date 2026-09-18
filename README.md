# DeePAW–HALF

**HALF** is the CPU/CUDA Fortran implementation of the fixed-density electronic-
structure reconstruction performed by DeePAW–HAPPY. The project name expands
to **H**arris **A**ssociative **L**inearized Augmented Plane Wave **F**ortran.

The scientific contract is intentionally the same as HAPPY:

```text
structure -> DeepAW -> smooth CHGCAR -> HALF -> H(k), S(k), bands, waves, energy
```

HALF is not a new learned Hamiltonian. It reconstructs the local potential and
linearized PAW/USPP-like operator from the predicted smooth density and the
matching VASP PAW dataset. HAPPY remains the executable numerical oracle while
the port is developed.

## Current milestone

Version 0.4 provides a numerically closed Si Gamma fixed-density path:

- VASP CHGCAR structure and smooth-grid reader that stops after exactly
  `NGX*NGY*NGZ` values;
- lattice, reciprocal-lattice and plane-wave basis data models;
- Gamma-point plane-wave selection using HAPPY's cutoff convention;
- smooth-electron-count validation;
- multi-dataset text POTCAR parsing;
- complete-grid Hartree, ionic, NLCC, LDA and PBE potentials on CPU and GPU;
- PAW reciprocal projectors and DION/QPAW overlap matrices on CPU and GPU;
- device-resident cuFFT potential construction, H/S assembly, Hermitian
  cleanup and cuSOLVER eigensolution without full-matrix host transfers;
- dense generalized H/S solvers using MKL on CPU and cuSOLVER on GPU;
- independent CPU Fortran and CUDA Fortran backends;
- CUDA 12.4 and CUDA 13.0 build presets;
- a correctness-checked backend benchmark;
- a unified `half` CLI, HAPPY-compatible command aliases, and parity-oriented
  tests.

The Si Gamma/DION path matches HAPPY to about `1.6e-11 eV`. Arbitrary k
points, potential-dependent MIMIC_US D, energy/forces and `vaspwave.h5` output
are tracked in [`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md); HALF is not
yet a complete replacement for every HAPPY workflow.

## Numerical model, derived step by step

HAPPY's target model, and HALF's long-term porting target, is a
**fixed-density** reconstruction: it reads the smooth valence density
`rho~(r)` from `CHGCAR` and does not run an SCF cycle. At a chosen cutoff, the
wavefunction is expanded in plane waves `|k+G>` and the following generalized
Hermitian problem is solved:

```text
H(k) c_n = epsilon_n S(k) c_n
```

For the full HAPPY `MIMIC_US` operator, the dense matrices have the form

```text
H_GG' = |k+G|^2/2 * delta_GG' + Veff(G-G')
        + sum_(I,i,j) beta_i^I(G) D_ij^I beta_j^I*(G')

S_GG' = delta_GG' + sum_(I,i,j) beta_i^I(G) Q_ij^I beta_j^I*(G')
```

The first term is kinetic energy. `Veff(G-G')` is the Fourier coefficient of
the fixed-density local effective potential. The last terms are the PAW/USPP
augmentation corrections: `beta` are reciprocal PAW projectors, `Q` augments
the overlap, and `D` is the onsite Hamiltonian matrix.

The potential is derived from the fixed smooth density as

```text
Veff(r) = Vion_local(r) + VH[rho~](r) + Vxc[rho~ + rho_core](r)
VH(G)   = 4*pi*e^2*rho~_G/(Omega*|G|^2),  G != 0
D_ij^I = DION_ij^I + integral Veff(r) QDEP_ij^I(r) dr
```

Here `rho~_G` is HAPPY's electron-number-normalized density coefficient,
`Omega` is the cell volume, and `e^2` is the electrostatic conversion factor
in the eV/angstrom convention. `rho_core` is the POTCAR partial-core density
used by NLCC only; it is not added to the Hartree density. `DION` is the frozen
onsite term from POTCAR, whereas `QDEP` supplies the potential-dependent
MIMIC_US contribution. `VH(G=0)` is a potential gauge and is set to zero.

### What HALF implements today

HALF 0.4 implements the **DION subset** of the equation above:

```text
D_ij^I = DION_ij^I
DeltaD_ij^I = integral Veff(r) QDEP_ij^I(r) dr   # not implemented
```

Consequently, the completed numerical claim is Si Gamma/DION parity with
HAPPY, not full MIMIC_US parity. The CLI rejects `--uspp-dij`; the
potential-dependent `QDEP`/`DeltaD` construction, matrix-free operator,
arbitrary-k bands, energy and forces remain planned. The authoritative status
is [`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md).

For the implemented DION problem, `S` is positive definite and the generalized
problem can be reduced through `S = L L^H` to `L^-1 H L^-H y = epsilon y`, then
`c = L^-H y`. MKL and cuSOLVER perform this dense generalized-Hermitian solve.

The execution path follows the derivation directly:

```text
CHGCAR + POTCAR
  -> parse structure, smooth density, PAW datasets
  -> select the cutoff plane-wave basis
  -> FFT construction of Veff (Hartree + local ionic + XC/NLCC)
  -> reciprocal PAW projectors and DION/QPAW matrices
  -> dense H and S assembly
  -> generalized Hermitian EVD (MKL on CPU, cuSOLVER on GPU)
  -> eigenvalues and a JSON validation report
```

The CUDA path keeps FFT potential construction, projector work, dense matrix
assembly, and eigensolution resident on the device. It avoids transfers of the
full complex H/S matrices; the cubic dense eigensolve consequently dominates
as the plane-wave count grows. This device-resident path does not change the
current DION-only scientific scope.

## Large-matrix performance

The HfO2 benchmark uses PBE, 520 eV, 60 bands and a 3407-by-3407 complex128
generalized eigensystem. CPU commands were pinned to one logical CPU and all
BLAS/OpenMP thread controls were set to one.

| implementation | resource | wall time |
| --- | --- | ---: |
| HALF CUDA | RTX PRO 6000 | 1.60 s |
| HALF CUDA | RTX 4090 | 1.795 s |
| HALF CPU Fortran | one logical CPU | 43.94 s |
| HAPPY Python (`--uspp-dij`) | one logical CPU | 64.23 s |

Thus HALF CPU is 1.46x faster than HAPPY Python on one core. The RTX 4090 and
RTX PRO 6000 are respectively 24.47x and 27.41x faster than the one-core HALF
CPU path. The input, cutoff and basis size are matched, but this is not a
numerically equivalent HfO2 comparison: HALF uses DION while the timed HAPPY
command used `--uspp-dij`. HfO2 numerical parity remains an independent
validation gate. Raw measurements and exact environment controls are in
[`docs/validation/hfo2_single_core_benchmark.json`](docs/validation/hfo2_single_core_benchmark.json).

## Requirements

- Linux x86-64 with an NVIDIA GPU;
- NVIDIA HPC SDK (`nvfortran`);
- CMake >= 3.24 and Make;
- GPU compute capability selected with `HALF_GPU_CC` (default `86`, RTX 3080).

## Build

The recommended entry point is the repository's CMake driver. It configures an
architecture-specific build directory, builds, and runs CTest:

```bash
tools/half-cmake cpu all
tools/half-cmake cuda12 all       # GPU architecture auto-detected
tools/half-cmake cuda13 all
HALF_GPU_CC=89 tools/half-cmake cuda12 all
```

Set `MKLROOT` (or pass `-DHALF_MKL_ROOT=...`) to enable the CPU complete-grid
FFT and dense Gamma solver. The CUDA Gamma pipeline uses cuFFT and cuSOLVER and
does not require oneMKL; without MKL, `half gamma --backend cuda` remains fully
available while the CPU Gamma backend is disabled.

The equivalent direct CMake commands are:

```bash
# Pure CPU Fortran: no .cuf source or CUDA runtime linkage.
cmake --preset cpu-release
cmake --build --preset cpu-release

# CUDA 12.4 with NVHPC 24.5 or newer.
cmake --preset cuda12-release
cmake --build --preset cuda12-release

# CUDA 13.0 requires NVHPC 25.9 or newer.
cmake --preset cuda13-release
cmake --build --preset cuda13-release
```

## Command-line interface

The unified CLI follows the argument style of HAPPY while allowing explicit
CPU/CUDA backend selection:

```bash
# Inspect a DeepAW or VASP smooth-density input.
./build/cuda12-cc89-release/half inspect CHGCAR.smooth --encut 400

# Reconstruct and validate the Gamma eigenspectrum; a JSON report is written.
./build/cuda12-cc89-release/half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 8 --xc pbe --backend cuda \
  --reference-eigenval EIGENVAL --output gamma_validation.json

# Inspect parsed POTCAR or PAW data.
./build/cpu-release/half potcar POTCAR
./build/cpu-release/half paw CHGCAR.smooth POTCAR --encut 400
```

`half-validate-gamma` is a compatibility alias for `half gamma`. The build also
creates `half-bands` and `half-energy` aliases so scripts can adopt the HAPPY-
style command names now; in version 0.4 those two commands exit with a clear
unsupported-feature error because arbitrary-k bands and total energy are not
yet numerically complete. Run `half --help` or `half gamma --help` for the
complete option list.

The older focused executables remain available during the transition. For
example:

```bash
./build/cuda12-release/half-inspect CHGCAR.smooth 400
```

## Accuracy policy

Every production module is accepted only after comparison with HAPPY on the
same input. Default tolerances and the required observables are documented in
[`docs/VALIDATION.md`](docs/VALIDATION.md). Performance work begins only after
the corresponding numerical row passes.
