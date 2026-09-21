# HAPPY to HALF porting matrix

Statuses are `done`, `active`, `planned`, and `blocked`.

| Capability | HAPPY source | HALF target | Status |
|---|---|---|---|
| CHGCAR header and smooth grid | `chgcar.py` | `half_chgcar.F90` | done |
| Crystal and reciprocal lattice | `crystal.py` | `half_types.F90` | done |
| Plane-wave basis and FFT indices | `basis.py` | `half_basis.F90` | done; exact Si parity |
| CUDA device bootstrap/kernel | n/a | `half_cuda.cuf` | done |
| HDF5 charge/structure and embedded POTCAR input | `chgcar.py`, `hdf5io.py`, `potcar.py` | `half_vaspwave.F90`, `half_hdf5_bridge.c` | done; native HDF5 C ABI avoids compiler-specific Fortran module files |
| Multi-dataset POTCAR parser and structure validation | `potcar.py` | `half_potcar.F90` | done; dataset count/order validation, Si and HfO2 parsed |
| LDA/PBE and NLCC | `xc.py`, `potential.py` | `half_xc.F90`, `half_cuda_potential.cuf` | done on CPU/GPU; Si and HfO2 CUDA parity |
| Local ionic and Hartree potential | `potential.py` | `half_potential.F90`, `half_cuda_potential.cuf` | done on CPU/GPU; explicit CHGCAR Fortran-to-CUDA grid reorder |
| PAW projectors and overlap | `paw.py` | `half_paw.F90`, `half_cuda_assembly.cuf` | done on CPU/GPU; Si and HfO2 parity |
| Potential-dependent MIMIC_US D | `uspp_dij.py` | `half_uspp.F90`, `half_cuda_uspp.cuf` | done on CPU/CUDA: periodic cubic B-spline sampling, QDEP multipoles and atom-dependent D |
| Dense H/S parity solver | `reconstruct_bandstructure.py` | `half_dense_solver.F90`, `half_cuda_assembly.cuf`, `half_cuda_solver.cuf` | done at arbitrary k points; MKL CPU and device-resident CUDA/cuSOLVER, DION or MIMIC_US |
| Matrix-free H/S application | `hamiltonian.py` | `half_cuda_assembly.cuf` | done; batched GPU H*Psi/S*Psi with cuFFT local potential, kinetic kernel and cuBLAS PAW contractions; no NPL-by-NPL matrices |
| Block iterative eigensolver | n/a | `half_cuda_iterative_solver.cuf` | done; restarted all-band optimization, GPU residual/preconditioner/S-orthogonalization, NBANDS/2NBANDS cuSOLVER Rayleigh-Ritz and configurable early stop |
| Arbitrary single k point | `basis.py`, `reconstruct_bandstructure.py` | `half_basis.F90`, `half_cli.F90` | done on CPU/CUDA; `--kpoint KX KY KZ` |
| Band paths, k meshes and symmetry reduction | `kpoints.py`, ASE | `half_kpoints.F90`, `half_cli.F90` | done; explicit reciprocal KPOINTS, ASE-identical automatic SC/FCC/BCC paths, full Gamma meshes and spglib irreducible meshes |
| Occupations, EIGENVAL, Ewald and energy | `occupations.py`, `ewald.py`, `total_energy.py` | `half_energy.F90`, `half_cli.F90` | done on CPU/CUDA; zero/finite-T occupations, reconstructed and reference-EIGENVAL paths; Si/HfO2 component parity |
| Analytic forces | not present in HAPPY | `half_energy.F90`, `half_local_forces.F90`, `half_forces.F90` | in progress: analytic Ewald, reciprocal local and generalized PAW `(D-epsilon Q)` terms implemented and derivative-tested; augmentation, NLCC and fixed-density Harris corrections remain before production use |
| Finite-difference force oracle | `total_energy.py` | `half_cli.F90` | validation only; Si parity within 2.3e-9 eV/Angstrom, not the production force algorithm |
| `vaspwave.h5` output | `vaspwave.py` | `half_vaspwave.F90`, `half_hdf5_bridge.c` | done for EVD bands/energy; VASP FFT coefficient order, float32 complex packing and charge round trip validated |
| MPI k-point distribution | n/a | `half_parallel.F90`, `half_cli.F90` | done for CPU bands, energy and finite-difference forces; cyclic distribution, rank-0 output and collective reconstruction; 3.903x on 4 ranks for the 36-k-point HfO2 case |
| HAPPY-style CLI artifacts | `happy-bands`, `happy-energy` | `half_cli.F90`, `half_artifacts.F90` | done; bands JSON/CSV/NPZ/PNG and energy JSON/NPZ, including gap and overlap diagnostics |
| Reusable application API | n/a | `half_library.F90`, `half_c_api.F90`, `include/half.h` | done; ABI v1 opaque contexts, C/Fortran bindings, CPU H/S/apply/solve and fully device-resident CUDA solve, CMake/pkg-config exports; CUDA device-pointer apply is an additive future extension |

## Ordered delivery gates

1. Input and basis parity.
2. Local potential parity on the complete FFT grid.
3. Gamma-point H/S matrix and eigenvalue parity for Si.
4. Arbitrary-k and high-symmetry band path parity.
5. Uniform k-mesh energy parity and VASP wavefunction output.
6. Matrix-free block solver optimization for larger systems (completed as the
   CUDA `acc` solver).
7. Optional multi-GPU, band and FFT-domain scaling beyond the completed CPU
   MPI k-point layer.
