# HAPPY to HALF porting matrix

Statuses are `done`, `active`, `planned`, and `blocked`.

| Capability | HAPPY source | HALF target | Status |
|---|---|---|---|
| CHGCAR header and smooth grid | `chgcar.py` | `half_chgcar.F90` | done |
| Crystal and reciprocal lattice | `crystal.py` | `half_types.F90` | done |
| Plane-wave basis and FFT indices | `basis.py` | `half_basis.F90` | done; exact Si parity |
| CUDA device bootstrap/kernel | n/a | `half_cuda.cuf` | done |
| HDF5 charge/structure input | `chgcar.py`, `hdf5io.py` | `half_hdf5.F90` | planned |
| Multi-dataset POTCAR parser | `potcar.py` | `half_potcar.F90` | done; Si and HfO2 parsed |
| LDA/PBE and NLCC | `xc.py`, `potential.py` | `half_xc.F90`, `half_cuda_potential.cuf` | done on CPU/GPU; Si and HfO2 CUDA parity |
| Local ionic and Hartree potential | `potential.py` | `half_potential.F90`, `half_cuda_potential.cuf` | done on CPU/GPU; explicit CHGCAR Fortran-to-CUDA grid reorder |
| PAW projectors and overlap | `paw.py` | `half_paw.F90`, `half_cuda_assembly.cuf` | done on CPU/GPU; Si and HfO2 parity |
| Potential-dependent MIMIC_US D | `uspp_dij.py` | `half_uspp.F90`, `half_cuda_uspp.cuf` | done on CPU/CUDA: periodic cubic B-spline sampling, QDEP multipoles and atom-dependent D |
| Dense H/S parity solver | `reconstruct_bandstructure.py` | `half_dense_solver.F90`, `half_cuda_assembly.cuf`, `half_cuda_solver.cuf` | done at Gamma; MKL CPU and device-resident CUDA/cuSOLVER, DION or MIMIC_US |
| Matrix-free H/S application | `hamiltonian.py` | `half_operator.cuf` | planned |
| Block iterative eigensolver | n/a | `half_lobpcg.cuf` | planned |
| k meshes and symmetry reduction | `kpoints.py` | `half_kpoints.F90` | planned |
| Occupations, Ewald and energy | `occupations.py`, `ewald.py`, `total_energy.py` | `half_energy.*` | planned |
| Analytic forces | theory only | `half_forces.cuf` | planned |
| Finite-difference force oracle | `total_energy.py` | test driver | planned |
| `vaspwave.h5` output | `vaspwave.py` | `half_vaspwave.F90` | planned |
| MPI k-point distribution | n/a | `half_mpi.F90` | blocked on single-GPU parity |

## Ordered delivery gates

1. Input and basis parity.
2. Local potential parity on the complete FFT grid.
3. Gamma-point H/S matrix and eigenvalue parity for Si.
4. Arbitrary-k and high-symmetry band path parity.
5. Uniform k-mesh energy parity and VASP wavefunction output.
6. Matrix-free block solver parity.
7. Multi-GPU and MPI scaling.
