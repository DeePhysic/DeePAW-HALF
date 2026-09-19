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
find_package(HALF 0.5 CONFIG REQUIRED)
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
2. Create a context with `half_create_from_files_backend()`.  Select LDA/PBE,
   CPU/CUDA/AUTO, EVD/EVJ, ENCUT, and MIMIC_US QDEP at runtime.  The shorter
   `half_create_from_files()` is ABI-compatible shorthand for AUTO + EVD.
3. For each fractional reciprocal k point, call `half_get_basis_size()` and
   optionally `half_get_basis()` to obtain HALF's exact G-vector ordering.
4. Choose one of three integration levels:
   `half_assemble_hs()` returns dense H/S; `half_apply_hs()` applies them to a
   state block; `half_solve_kpoint()` returns eigenpairs.
5. Release the opaque handle with `half_destroy()`.

All matrices and state blocks are complex128, column-major BLAS arrays.  Energy
units are eV, lengths are Angstrom, and k points are fractional reciprocal
coordinates.  Passing a null eigenvector pointer requests eigenvalues only.
Each routine returns a status code and fills a caller-owned error buffer.

In a CUDA build, AUTO selects CUDA.  `half_solve_kpoint()` runs potential
construction, MIMIC_US/QDEP, dense H/S assembly, and the generalized
eigensolution on the GPU, copying back only the requested eigenpairs.  ABI v1
host H/S export and H/S application are CPU interfaces; CUDA contexts return
`HALF_ERROR_UNAVAILABLE` for those two calls.  A future device-pointer extension
can add them without changing existing symbols.  In a CPU build, AUTO selects
oneMKL.  EVJ is CUDA-only, while EVD is available on both.

## VASP adapter pattern

A minimal VASP patch creates one context after the fixed smooth density and
matching POTCAR are available.  At each VASP k point it queries `npw` and the
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

The current file constructor consumes CHGCAR/vaspwave.h5 plus POTCAR (or the
POTCAR embedded in vaspwave.h5).  This keeps the first VASP adapter small and
auditable.  New constructors can be added without changing existing symbols;
the reserved opaque context and ABI-version query are the extension boundary.
