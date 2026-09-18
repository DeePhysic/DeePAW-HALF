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

Version 0.3 provides a numerically closed Si Gamma fixed-density path:

- VASP CHGCAR structure and smooth-grid reader that stops after exactly
  `NGX*NGY*NGZ` values;
- lattice, reciprocal-lattice and plane-wave basis data models;
- Gamma-point plane-wave selection using HAPPY's cutoff convention;
- smooth-electron-count validation;
- multi-dataset text POTCAR parsing;
- complete-grid Hartree, ionic, NLCC, LDA and PBE potentials;
- PAW reciprocal projectors and DION/QPAW overlap matrices;
- dense generalized H/S solvers using MKL on CPU and cuSOLVER on GPU;
- independent CPU Fortran and CUDA Fortran backends;
- CUDA 12.4 and CUDA 13.0 build presets;
- a correctness-checked backend benchmark;
- an `half-inspect` executable and parity-oriented tests.

The Si Gamma/DION path matches HAPPY to about `1.6e-11 eV`. Arbitrary k
points, potential-dependent MIMIC_US D, energy/forces and `vaspwave.h5` output
are tracked in [`docs/PORTING_MATRIX.md`](docs/PORTING_MATRIX.md); HALF is not
yet a complete replacement for every HAPPY workflow.

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

Set `MKLROOT` (or pass `-DHALF_MKL_ROOT=...`) to enable the complete-grid FFT
and dense Gamma solver. Without oneMKL, the input, basis, PAW inspection, and
backend benchmark targets remain available.

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

Inspect a DeepAW or VASP smooth-density input:

```bash
./build/cuda12-release/half-inspect CHGCAR.smooth 400
```

For the configured project environment, use:

```bash
scripts/remote_build.sh hz.icqms.group 8122
scripts/remote_benchmark.sh hz.icqms.group 8122 5000
scripts/remote_cuda13_build.sh hz.icqms.group 8122 5000
```

The optional fourth argument overrides the GPU compute capability; by default
it is detected with `nvidia-smi` (for example, `86` for RTX 3080 and `89` for
RTX 4090). Architecture-specific build directories prevent collisions when
nodes share the same user directory.

The CUDA 12 helpers synchronize HALF and the licensed-data-free Si CHGCAR,
then build with the cluster's NVHPC 24.5 installation and run on its RTX 3080.
The CUDA 13 helper is container-free: it uses the existing mamba installation
to create the persistent `half-cuda13` environment under the shared user
installation, verifies NVIDIA's native RHEL/Rocky NVHPC 25.9 RPM, and extracts
NVHPC under `/share/home/limusen/app`. No mamba environment is installed in
`/tmp`.
POTCAR files are never copied by these scripts.

## Accuracy policy

Every production module is accepted only after comparison with HAPPY on the
same input. Default tolerances and the required observables are documented in
[`docs/VALIDATION.md`](docs/VALIDATION.md). Performance work begins only after
the corresponding numerical row passes.
