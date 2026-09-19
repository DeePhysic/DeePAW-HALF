# HALF architecture

HALF separates the frozen-density physics from the GPU execution policy.

```text
CHGCAR/HDF5   POTCAR/vaspout.h5
     |               |
     v               v
 charge_io         paw_data
     |               |
     +----> local_potential <---- XC/NLCC
                    |
             FFT-grid Veff(G)
                    |
 k mesh -> basis -> H/S operator -> block eigensolver
                    |                    |
                    +------ observables-+
                              |
                 bands / energy / vaspwave.h5
```

## Execution rules

1. Large FFT-grid, wavefunction and projector arrays remain device-resident.
2. cuFFT performs real/reciprocal transforms.
3. CUDA Fortran kernels construct PAW projectors and dense DION/QPAW terms;
   cuBLAS contraction is reserved for larger multi-species workloads.
4. cuSOLVER is used only by the small-system dense reference path.
5. The production path applies H and S matrix-free and uses a block iterative
   eigensolver.
6. Optional CPU MPI distributes complete k-point solves cyclically across ranks;
   band and FFT-domain decomposition remain future scaling extensions.
7. All reductions that affect regression results have a deterministic mode.

## Build variants

| Preset | Compiler/toolkit | CUDA source linked | Intended use |
|---|---|---:|---|
| `cpu-release` | NVFORTRAN or GNU Fortran + OpenMP | no | portable CPU reference and fallback |
| `cuda12-release` | NVHPC 24.5+, CUDA 12.4 | yes | current `cmu001` production baseline |
| `cuda13-release` | NVHPC 25.9+, CUDA 13.0 | yes | CUDA 13 compatibility and performance |

The CPU preset excludes `half_cuda.cuf` at CMake configure time. It is not a
CUDA binary with GPU calls disabled. CUDA 12 and CUDA 13 compile the same
Fortran sources so version comparisons do not mix in algorithm changes.

## Source layout

| HALF module | HAPPY oracle |
|---|---|
| `half_chgcar` | `happy/chgcar.py` |
| `half_basis` | `happy/basis.py` |
| `half_potcar` | `happy/potcar.py` |
| `half_potential` | `happy/potential.py`, `happy/xc.py` |
| `half_cuda_potential` | GPU Hartree, ionic, NLCC, LDA and PBE pipeline |
| `half_paw` | `happy/paw.py`, `happy/uspp_dij.py` |
| `half_cuda_assembly` | GPU projector and device-resident H/S construction |
| `half_operator` | `happy/hamiltonian.py` |
| `half_solver` | SciPy generalized `eigh` call sites |
| `half_energy` | `happy/total_energy.py`, `happy/ewald.py` |
| `half_vaspwave` | `happy/vaspwave.py` |
| `half_parallel` | MPI lifecycle, cyclic k-point ownership and reductions |

The public command names planned for compatibility are `half-bands`,
`half-energy`, and `half-validate-gamma`. The current bootstrap command is
`half-inspect` so incomplete physics cannot be mistaken for a finished result.
