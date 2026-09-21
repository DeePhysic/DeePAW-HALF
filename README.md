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

Version 0.6 provides a numerically closed fixed-density path plus a reusable
library interface:

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
- CUDA solver selection: dense `evd`, Jacobi `evj`, index-range `evx`, or the
  matrix-free all-band `acc` path;
- independent CPU Fortran and CUDA Fortran backends;
- optional CPU MPI distribution over independent k points;
- CUDA 12.4 and CUDA 13.0 build presets;
- a correctness-checked backend benchmark;
- a unified `half` CLI, HAPPY-compatible command aliases, and parity-oriented
  tests;
- installable `libhalf.so` with a versioned C ABI, portable Fortran bindings,
  CMake/pkg-config metadata, and runtime CPU/CUDA backend selection.

The CUDA MIMIC_US path reproduces HAPPY for both Si (725 plane waves) and HfO2
(3407 plane waves): the tested eigenvalues agree to about `1e-11 eV` or better.
Explicit multi-k bands, Gamma-centered full/irreducible spglib meshes,
occupations, Ewald and the fixed-density Harris total energy are also native
Fortran features. The Si total energy and every reported component agree with
HAPPY to better than `4e-12 eV`. The old central finite-difference force path
is retained only as a numerical oracle. A native analytic force implementation
now contains the Ewald, reciprocal-space local, and generalized PAW
`D-epsilon Q` Hellmann-Feynman derivatives; augmentation, NLCC, and fixed-density
Harris corrections are still being completed and must not yet be presented as
a production total force. Native `vaspwave.h5` output and HDF5
charge/structure/embedded-POTCAR input are supported. Dense H/S assembly,
matrix application, and k-point solution are exposed through `libhalf`.

The derivation of the matrix-free PAW operator, constrained all-band
minimization, residual preconditioning, S-orthogonalization, restarted
Rayleigh-Ritz update, complexity, and direct VASP handoff is given in
[Matrix-free Harris all-band acceleration](docs/HARRIS_ACC_THEORY.md).

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
matrix application and wavefunction export are implemented. Finite-difference
forces are validation-only while the analytic PAW force is completed. The authoritative status is
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
- an MPI Fortran implementation for the optional multi-process CPU build;
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

# CPU MPI build (opt-in, so serial builds have no MPI startup overhead).
cmake -S . -B build/cpu-mpi -DHALF_ENABLE_CUDA=OFF \
  -DHALF_ENABLE_MPI=ON -DHALF_MKL_ROOT="$MKLROOT"
cmake --build build/cpu-mpi -j

# Intel oneAPI 2023: use the MPI wrapper as the compiler.  mpiifort is faster
# than that release's mpiifx for the validated HfO2 case.
module load intel/oneapi2023
cmake -S . -B build/cpu-oneapi-mpi -DCMAKE_Fortran_COMPILER=mpiifort \
  -DHALF_ENABLE_CUDA=OFF -DHALF_ENABLE_MPI=ON \
  -DHALF_MKL_ROOT="$MKLROOT"
cmake --build build/cpu-oneapi-mpi -j

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

# Matrix-free Harris initialization: block H*Psi/S*Psi, all-band Ritz solve,
# and early stopping at an initial-wavefunction residual tolerance.
./build/cuda13-cc89-release/half gamma CHGCAR.smooth POTCAR \
  --encut 400 --bands 64 --backend cuda --solver acc --uspp-dij \
  --acc-tol 1e-3 --acc-max-iter 40 --acc-block-size 16 \
  --output gamma_acc.json

# Reconstruct several k points from a VASP explicit reciprocal KPOINTS file.
./build/cuda12-cc89-release/half bands CHGCAR.smooth POTCAR KPOINTS \
  --encut 400 --bands 12 --backend cuda --uspp-dij --output bands.json

# Or generate the ASE-compatible automatic cubic high-symmetry path.
./build/cuda12-cc89-release/half bands CHGCAR.smooth POTCAR \
  --encut 400 --bands 12 --npoints 60 --backend cuda --output-prefix bands

# Evaluate the symmetry-reduced fixed-density Harris energy.
./build/cuda12-cc89-release/half energy CHGCAR.smooth POTCAR \
  --encut 400 --kspacing 0.5 --bands 12 --backend cuda \
  --vaspwave-h5 vaspwave.h5 --forces --force-step 0.001 --output-prefix energy

# Alternatively consume a matching VASP EIGENVAL on an explicit mesh.
./build/cuda12-cc89-release/half energy CHGCAR.smooth POTCAR \
  --kpoints-file KPOINTS --reference-eigenval EIGENVAL --bands 60 --output energy.json

# Distribute k points over four single-thread CPU ranks.
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun --bind-to core -np 4 build/cpu-mpi/half bands \
  CHGCAR.smooth POTCAR KPOINTS --backend cpu --bands 60 --output bands.json

# Inspect parsed POTCAR or PAW data.
./build/cpu-release/half potcar POTCAR
./build/cpu-release/half paw CHGCAR.smooth POTCAR --encut 400
```

`--solver acc` is the production large-basis path. It never forms the
`NPL x NPL` dense H/S matrices: local-potential FFTs, kinetic terms, PAW
projector contractions, residuals, preconditioning, S-orthogonalization and
block rotations remain on the GPU. cuSOLVER handles only the `NBANDS` and
`2*NBANDS` Rayleigh-Ritz problems. The JSON report records the iteration count
and final maximum residual. Use `--solver evx` as the dense partial-spectrum
reference. For VASP, `HALF_MODE=ACC` selects the same path and uses VASP's
`NBANDS`; its default residual tolerance is `1e-4 eV`.

On the 20-atom CsPbBr3 validation case (`NPL=20640--20780`, 12 k points,
105 bands), integrated VASP `HALF_MODE=ACC` retained the same five SCF loops
and final energy as dense `VASP_LIKE`, while reducing total elapsed time from
1861.0 to 275.5 seconds (6.75x). See the
[RTX PRO 6000 record](docs/validation/cspbbr3_vasp_acc_pro6000.json).

MPI is implemented for the CPU `bands` and `energy` paths. Rank `r` solves
k points `r+1, r+1+nranks, ...`; collective reductions restore the ordered
eigenvalue and plane-wave arrays, and only rank 0 writes JSON. Finite-difference
forces use the same distribution for every displaced structure. Use one
BLAS/OpenMP thread per rank unless deliberately testing hybrid MPI+OpenMP.
Multi-rank CUDA and multi-rank `vaspwave.h5` output are rejected explicitly.

On the 36-k-point HfO2 CPU case, 1/2/4 ranks took 116.59/58.44/29.87 seconds,
or 1.995x and 3.903x speedups. The 2- and 4-rank eigenvalues were bit-for-bit
identical to serial. See
[`docs/validation/hfo2_cpu_mpi.json`](docs/validation/hfo2_cpu_mpi.json).

`half-validate-gamma`, `half-bands`, and `half-energy` are compatibility aliases
for their corresponding subcommands. Run `half COMMAND --help` for the complete
option list.

`bands --output-prefix NAME` writes `NAME.json`, `NAME.csv`, `NAME.npz`, and a
dependency-free `NAME.png` plot.  The report includes VBM, CBM, sampled gap,
Gamma direct gap, path distance, basis size, and overlap diagnostics.
`energy --output-prefix NAME` writes JSON plus a HAPPY-compatible NPZ containing
k points, weights, eigenvalues, occupations, and forces.

## Library API and VASP integration

The same implementation is available as installable `libhalf.so`.  Its stable
C ABI uses opaque handles and runtime integer selectors in the libxc style;
`half_api.f90` provides compiler-neutral `bind(C)` interfaces for Fortran hosts.
It exposes basis queries, dense H/S assembly, block H/S application, and direct
k-point solution.  CPU contexts support all three compute levels; CUDA contexts
provide the fully device-resident direct solve in ABI v1.  Installed CMake
clients link `HALF::half`; pkg-config clients use `half.pc`.  CPU/CUDA and
EVD/EVJ are selected when the context is created.

The VASP adapter also transfers its lattice, fractional positions, species
indices, and dense charge-grid dimensions directly into HALF memory. HALF
deep-copies this request metadata. It can now use the same descriptor to call a
DeePAW-eSCN endpoint, translate its C-order float32 grid, apply explicit POTCAR
valence normalization, and construct the HALF context without changing the
VASP-to-HALF structure ABI. Select this path explicitly in INCAR with
`LHALF_API=.TRUE.` and `HALF_ESCN_URL=http://127.0.0.1:8265`.

See [the API guide](docs/API.md), the
[tested VASP 6.6.0 integration](docs/VASP_INTEGRATION.md), the standalone
[C example](examples/api/half_c_example.c), and the original
[VASP-side adapter module](examples/vasp/half_vasp_init.F).
The combined article on direct DeepAW-to-band-structure reconstruction,
energy/force evaluation, and HALF-initialized VASP SCF acceleration includes
pointwise band validation, comparisons with reported machine-learning
potentials, the MP-85 mean iteration speedup, and large-basis ACC scaling in the
[validation report](docs/validation/HFO2_VASP_HALF_VS_SAD_KSPACING035.md).

The private `deepaw-half-acc` branch vendors the complete VASP 6.6.0
source under `vendor/vasp-6.6.0`. On the Pro 6000 host, one command configures
and builds HALF, the oneMKL FFTW wrapper, and `vasp_std`, then runs the HALF
tests:

```bash
tools/half-cmake vasp all
```

The tool detects the newest NVHPC under `/opt/nvidia/hpc_sdk/Linux_x86_64`, its
bundled CUDA version, GPU compute capability, HPC-X MPI, and oneMKL under
`/opt/intel/oneapi`. The executable is written to
`build/vasp-cuda<CUDA>-cc<CC>-release/vasp-6.6.0/bin/vasp_std`.

To additionally build VASP's own OpenACC GPU port, use:

```bash
tools/half-cmake vasp all -- -DHALF_VASP_ENABLE_OPENACC=ON
```

The default is the configuration validated on the Pro 6000: the VASP driver
uses its NVHPC/HPC-X CPU path, while HALF performs potential construction, H/S
assembly, and the generalized eigensolve on the CUDA GPU.

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
