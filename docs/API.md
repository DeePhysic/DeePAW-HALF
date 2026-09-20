# libhalf integration API

`libhalf` is the stable, compiler-neutral integration layer for VASP and other
electronic-structure programs.  It follows the same model as libxc: callers
select an implementation by integer identifiers, create an opaque context,
query its capabilities and basis, and evaluate through a C ABI.  The public ABI
is versioned independently through `HALF_ABI_VERSION`.

## Build and consume

```bash
cmake -S . -B build -DHALF_ENABLE_CUDA=OFF -DHALF_BUILD_SHARED_LIBRARY=ON
cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

Downstream CMake projects use the exported target:

```cmake
find_package(HALF 0.6 CONFIG REQUIRED)
target_link_libraries(mycode PRIVATE HALF::half)
```

The installation also provides `half.pc`, `half.h`, and the portable Fortran
interface source `half_api.f90`.  See `examples/api/half_c_example.c` and
`examples/cmake/CMakeLists.txt` for a complete external build.

CUDA clients must perform their final link with `nvfortran -cuda` (CMake:
`LINKER_LANGUAGE Fortran` plus the installed CUDA version/CC flags) so NVHPC
registers CUDA Fortran device code.  The exported example does this
automatically.  CPU-only `libhalf` can be linked by an ordinary C compiler.

## Lifecycle and operations

1. Check `half_get_abi_version()` and `half_get_capabilities()`.
2. Create a context with `half_create_from_files_backend()` or obtain its
   density from DeePAW-eSCN with `half_create_from_escn()`. Select LDA/PBE,
   CPU/CUDA/AUTO, EVD/EVJ/VASP/ACC, ENCUT, and MIMIC_US QDEP at runtime. The VASP
   selector is the CUDA index-range solve for the `NBANDS` requested by the
   solve call. ACC is the matrix-free block H*Psi path and likewise consumes
   the requested `NBANDS`. The shorter
   `half_create_from_files()` is ABI-compatible shorthand for AUTO + EVD.
3. A host may attach its already parsed structure with
   `half_set_request_geometry()`: lattice vectors, fractional positions,
   one-based species indices, and the three-dimensional real-space grid.
4. For each fractional reciprocal k point, call `half_get_basis_size()` and
   optionally `half_get_basis()` to obtain HALF's exact G-vector ordering.
5. Choose one of three integration levels:
   `half_assemble_hs()` returns dense H/S; `half_apply_hs()` applies them to a
   state block; `half_solve_kpoint()` returns eigenpairs.
6. Release the opaque handle with `half_destroy()`.

All matrices and state blocks are complex128, column-major BLAS arrays.  Energy
units are eV, lengths are Angstrom, and k points are fractional reciprocal
coordinates.  Passing a null eigenvector pointer requests eigenvalues only.
Each routine returns a status code and fills a caller-owned error buffer.

`half_set_request_geometry()` takes a `half_request_geometry_v1` descriptor and
deep-copies every array before returning; the caller may immediately release
its buffers. `lattice` is row-major as `[lattice_vector][Cartesian_component]`,
and `positions_fractional` is atom-major. `half_get_request_geometry()` can
query or copy the retained values. The file-backed path retains this metadata
for host-side auditing. The eSCN constructor consumes the same structure
directly.

## DeePAW-eSCN remote density

Version 0.6 adds two levels of remote inference support:

- `half_escn_health()` and `half_escn_predict()` expose the HTTP service
  directly. The caller supplies atomic numbers, Cartesian positions in
  Angstrom, row-major cell vectors, and `[nx,ny,nz]`. The result contains the
  server's C-order float32 density and optionally `nu`, `alpha`, `beta`, and
  `risk`.
- `half_create_from_escn()` accepts the existing
  `half_request_geometry_v1`, derives atomic numbers and valence charge from
  the matching POTCAR, calls `/v1/predict`, converts C-order data to HALF's
  CHGCAR-compatible internal order, and creates a normal reusable context.

Use `HALF_DENSITY_NORMALIZE_VALENCE` for the numerical HALF path. It rescales
the returned grid so its mean equals the sum of POTCAR ZVAL values. Use
`HALF_DENSITY_RAW` only when the model output is already in the exact
normalization expected downstream. Neither mode clips negative predictions.

The API capability is reported as `HALF_CAP_ESCN_API`. CMake uses libcurl when
available; otherwise Unix builds use a dependency-free HTTP transport. The
built-in transport supports `http://`, which is the expected scheme behind the
SSH tunnel. HTTPS requires libcurl. Build-time disablement is available through
`-DHALF_ENABLE_ESCN_API=OFF`.

For the documented server, establish the tunnel first and then run:

```bash
ssh -N -L 8265:127.0.0.1:8265 cmu-pro6000
./build/cpu-release/half-escn-example http://127.0.0.1:8265
```

`half_escn_predict()` enforces the server's 2,000,000-point and 2 MiB request
limits. The high-level constructor deliberately requests density only because
uncertainty is not consumed by the Hamiltonian. Applications that need NIG
fields should call the lower-level function and set
`HALF_ESCN_INCLUDE_UNCERTAINTY`.

In a CUDA build, AUTO selects CUDA.  `half_solve_kpoint()` runs potential
construction, MIMIC_US/QDEP, dense H/S assembly, and the generalized
eigensolution on the GPU, copying back only the requested eigenpairs.  ABI v1
host H/S export and H/S application are CPU interfaces; CUDA contexts return
`HALF_ERROR_UNAVAILABLE` for those two calls.  A future device-pointer extension
can add them without changing existing symbols.  In a CPU build, AUTO selects
oneMKL. EVJ, VASP partial-spectrum mode, and ACC are CUDA-only, while EVD is
available on both. ACC keeps wavefunction blocks, H*Psi/S*Psi, residuals, and
projected algebra on the device. cuSOLVER sees only the `NBANDS`/`2*NBANDS`
Rayleigh-Ritz matrices, not dense `NPL x NPL` H/S.

## VASP adapter pattern

A minimal VASP patch creates one context after the fixed smooth density and
matching POTCAR are available, then transfers `LATT_CUR%A`, `T_INFO%POSION`,
`T_INFO%ITYP`, and `GRIDC%NGPTAR` directly from memory into the context. At
each VASP k point it queries `npw` and the
G-vector map, translates VASP's coefficient ordering once, and then either:

- uses a CPU context and calls `half_apply_hs()` inside an iterative VASP
  eigensolver; or
- calls `half_solve_kpoint()` to use HALF's CPU or fully GPU-resident dense EVD.
- Host codes with their own plane-wave ordering can call
  `half_solve_kpoint_mapped()` with `(Gx,Gy,Gz)` triples. HALF validates exact
  basis equality and returns coefficients in the requested order. See
  [VASP_INTEGRATION.md](VASP_INTEGRATION.md).

The context is reusable across all k points.  Destroy it when the density,
structure, POTCAR, ENCUT, XC, or backend changes.  MPI processes should own
separate contexts; k-point distribution remains the caller's responsibility.
Do not invoke two mutating operations concurrently on the same context.

The file constructor consumes CHGCAR/vaspwave.h5 plus POTCAR (or the POTCAR
embedded in vaspwave.h5). The remote constructor consumes the existing VASP
geometry descriptor plus POTCAR; no `vasp2half` structure change is required.
The VASP adapter selects it explicitly through INCAR:

```text
LHALF_INIT = .TRUE.
LHALF_API = .TRUE.
HALF_MODE = VASP_LIKE
# or: HALF_MODE = ACC
HALF_ESCN_URL = http://127.0.0.1:8265
```

`HALF_MODE=TRADITIONAL` is the default. `HALF_MODE=VASP_LIKE` selects the
CUDA partial-spectrum solver and forwards VASP's `NBANDS` through the existing
mapped k-point ABI, so HALF returns only the requested lowest eigenpairs.
`DENSE` and `VASP` are accepted as shorter aliases.

`HALF_MODE=ACC` (aliases `MATRIX_FREE` and `MATRIX-FREE`) instead selects the
matrix-free block solver. It uses VASP's `NBANDS` and the initialization-oriented
default maximum residual of `1e-4 eV`, avoiding dense H/S allocation entirely.

With `LHALF_API=.FALSE.` (the default), HALF uses the established CHGCAR path
even if a `HALF_ESCN_URL` environment variable exists. When `LHALF_API=.TRUE.`,
the INCAR URL takes precedence; the environment variable is accepted only as a
backward-compatible fallback when the INCAR URL is empty. A missing endpoint is
an error, and a missing CHGCAR never silently enables network access. In the
current VASP integration this selects the density used inside HALF to
reconstruct the initial waves; VASP still reads CHGCAR for its own `ICHARG=1`
host density. Replacing that second copy requires a separate VASP charge-grid
injection and is not hidden in this constructor.
