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

Version 0.4 provides a numerically closed Gamma fixed-density path on CUDA:

- VASP CHGCAR structure and smooth-grid reader that stops after exactly
  `NGX*NGY*NGZ` values;
- lattice, reciprocal-lattice and plane-wave basis data models;
- Gamma-point plane-wave selection using HAPPY's cutoff convention;
- smooth-electron-count validation;
- multi-dataset text POTCAR parsing;
- complete-grid Hartree, ionic, NLCC, LDA and PBE potentials on CPU and GPU;
- PAW reciprocal projectors and DION/QPAW overlap matrices on CPU and GPU;
- atom-dependent MIMIC_US QDEP matrices on CPU and GPU, including periodic cubic
  B-spline sampling, two-Bessel compensation functions and Gaunt transforms;
- device-resident cuFFT potential construction, H/S assembly, Hermitian
  cleanup and cuSOLVER eigensolution without full-matrix host transfers;
- dense generalized H/S solvers using MKL on CPU and cuSOLVER on GPU;
- CUDA Gamma CLI solver selection: `--solver evd` (default divide-and-conquer)
  or `--solver evj` (Jacobi, through a device-pointer CUDA C++ bridge);
- independent CPU Fortran and CUDA Fortran backends;
- CUDA 12.4 and CUDA 13.0 build presets;
- a correctness-checked backend benchmark;
- a unified `half` CLI, HAPPY-compatible command aliases, and parity-oriented
  tests.

The CUDA MIMIC_US path reproduces HAPPY for both Si (725 plane waves) and HfO2
(3407 plane waves): the tested eigenvalues agree to about `1e-11 eV` or better.
Explicit multi-k bands, Gamma-centered full/irreducible spglib meshes,
occupations, Ewald and the fixed-density Harris total energy are also native
Fortran features. The Si total energy and every reported component agree with
HAPPY to better than `4e-12 eV`. Central finite-difference forces reproduce
HAPPY within `2.3e-9 eV/Angstrom` in the validated Si case. Native `vaspwave.h5` output and HDF5
charge/structure/embedded-POTCAR input are supported. Forces and the matrix-free
solver remain tracked in [`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md).

## Numerical model, derived step by step

HAPPY's target model, and HALF's long-term porting target, is a
**fixed-density** reconstruction: it reads the smooth valence density
`rho~(r)` from `CHGCAR` and does not run an SCF cycle. At a chosen cutoff, the
wavefunction is expanded in plane waves `|k+G>` and the following generalized
Hermitian problem is solved:

$$
H(\mathbf{k})\,\mathbf{c}_n
=\varepsilon_n S(\mathbf{k})\,\mathbf{c}_n .
$$

For the full HAPPY `MIMIC_US` operator, the dense matrices have the form

$$
H_{\mathbf G\mathbf G'}(\mathbf k)=
\frac{\hbar^2|\mathbf k+\mathbf G|^2}{2m_e}\delta_{\mathbf G\mathbf G'}
+V_{\mathrm{eff}}(\mathbf G-\mathbf G')
+\sum_{Iij}\beta_i^I(\mathbf G)D_{ij}^I\beta_j^{I*}(\mathbf G'),
$$

$$
S_{\mathbf G\mathbf G'}(\mathbf k)=
\delta_{\mathbf G\mathbf G'}
+\sum_{Iij}\beta_i^I(\mathbf G)Q_{ij}^I\beta_j^{I*}(\mathbf G').
$$

The first term is kinetic energy. `Veff(G-G')` is the Fourier coefficient of
the fixed-density local effective potential. The last terms are the PAW/USPP
augmentation corrections: `beta` are reciprocal PAW projectors, `Q` augments
the overlap, and `D` is the onsite Hamiltonian matrix.

The potential is derived from the fixed smooth density as

$$
V_{\mathrm{eff}}(\mathbf r)=V_{\mathrm{ion}}^{\mathrm{local}}(\mathbf r)
+V_H[\widetilde\rho](\mathbf r)
+V_{\mathrm{xc}}[\widetilde\rho+\rho_{\mathrm{core}}](\mathbf r),
$$

$$
V_H(\mathbf G)=\frac{4\pi e^2\widetilde\rho_{\mathbf G}}
{\Omega|\mathbf G|^2}\quad(\mathbf G\ne0),
\qquad
D_{ij}^I=D_{ij}^{I,\mathrm{ION}}
+\int V_{\mathrm{eff}}(\mathbf r)Q_{ij}^{I,\mathrm{DEP}}(\mathbf r-\mathbf R_I)\,d^3r.
$$

Here `rho~_G` is HAPPY's electron-number-normalized density coefficient,
`Omega` is the cell volume, and `e^2` is the electrostatic conversion factor
in the eV/angstrom convention. `rho_core` is the POTCAR partial-core density
used by NLCC only; it is not added to the Hartree density. `DION` is the frozen
onsite term from POTCAR, whereas `QDEP` supplies the potential-dependent
MIMIC_US contribution. `VH(G=0)` is a potential gauge and is set to zero.

### What HALF implements today

HALF now evaluates the full arbitrary-k MIMIC_US expression above on CPU and CUDA.
`--uspp-dij` enables the potential-dependent term. The effective potential is
prefiltered for periodic cubic B-spline interpolation on the GPU, sampled on
concentric angular grids around every atom, projected onto real spherical
harmonics, radially integrated against VASP's two-Bessel compensation
functions, and contracted with the AE-minus-PS multipole moments. Explicit
multi-point bands and symmetry-reduced fixed-density energies are implemented;
matrix-free application, forces and wavefunction export remain planned. The authoritative status is
[`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md).

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
  -> optional GPU QDEP construction and atom-dependent D = DION + DeltaD
  -> dense H and S assembly
  -> generalized Hermitian EVD (MKL on CPU, cuSOLVER on GPU)
  -> eigenvalues and a JSON validation report
```

The CUDA path keeps FFT potential construction, projector work, dense matrix
assembly, and eigensolution resident on the device. It avoids transfers of the
full complex H/S matrices; the cubic dense eigensolve consequently dominates
as the plane-wave count grows. This device-resident path does not change the
Gamma-point scientific scope.

## Large-matrix performance

The HfO2 benchmark uses PBE, 520 eV, 60 bands and a 3407-by-3407 complex128
generalized eigensystem. CPU commands were pinned to one logical CPU and all
BLAS/OpenMP thread controls were set to one.

| implementation | resource | wall time |
| --- | --- | ---: |
| HALF CUDA (`--uspp-dij`) | RTX PRO 6000 | 1.696 s median |
| HALF CPU Fortran (`--uspp-dij`) | one logical CPU | 43.790 s |
| HAPPY Python (`--uspp-dij`) | one logical CPU | 64.23 s |

The numerically equivalent CPU Fortran result is `1.47x` faster than HAPPY;
CUDA is `37.87x` faster than HAPPY and `25.82x` faster than CPU Fortran. Its
five fresh-process samples have a `1.696 s` median, split into
`0.258 s` potential construction, `0.349 s` H/S assembly including QDEP, and
`1.061 s` cuSOLVER time. The earlier 4090 DION-only measurement remains
historical data and is not presented as a MIMIC_US speedup. Raw measurements
and exact environment controls are in
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
  --encut 400 --bands 8 --xc pbe --backend cuda --uspp-dij \
  --reference-eigenval EIGENVAL --output gamma_validation.json

# Solve an arbitrary fractional reciprocal k point with the same operator.
./build/cuda12-cc89-release/half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 8 --backend cuda --uspp-dij \
  --kpoint 0.125 0.25 0.375 --output kpoint_validation.json

# Reconstruct several k points from a VASP explicit reciprocal KPOINTS file.
./build/cuda12-cc89-release/half bands CHGCAR.smooth POTCAR KPOINTS \
  --encut 400 --bands 12 --backend cuda --uspp-dij --output bands.json

# Or generate the ASE-compatible automatic cubic high-symmetry path.
./build/cuda12-cc89-release/half bands CHGCAR.smooth POTCAR \
  --encut 400 --bands 12 --npoints 60 --backend cuda --output-prefix bands

# Evaluate the symmetry-reduced fixed-density Harris energy.
./build/cuda12-cc89-release/half energy CHGCAR.smooth POTCAR \
  --encut 400 --kspacing 0.5 --bands 12 --backend cuda \
  --vaspwave-h5 vaspwave.h5 --forces --force-step 0.001 --output energy.json

# Inspect parsed POTCAR or PAW data.
./build/cpu-release/half potcar POTCAR
./build/cpu-release/half paw CHGCAR.smooth POTCAR --encut 400
```

`half-validate-gamma`, `half-bands`, and `half-energy` are compatibility aliases
for their corresponding subcommands. Run `half COMMAND --help` for the complete
option list.

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
