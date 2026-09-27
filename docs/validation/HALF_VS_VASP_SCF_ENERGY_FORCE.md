# HALF energy and force versus fully self-consistent VASP

This is the project-level energy/force comparison. The reference is a fully
self-consistent VASP calculation with `LMAXMIX=-1`; a frozen-density
`ICHARG=11` calculation is retained only as an implementation regression and
is not used in the accuracy table below.

Both calculations use PBE, `KSPACING=0.35`, `ISPIN=1`, `ISMEAR=0`, and
`SIGMA=0.02 eV`. VASP uses `EDIFF=1E-4`, `ALGO=All`, and converges the
electronic loop. Si uses 520 eV and HfO2 uses 521 eV.

## Total energy

| System | Atoms | HALF Harris (eV/cell) | self-consistent VASP (eV/cell) | HALF - VASP (eV/cell) | difference (meV/atom) |
|---|---:|---:|---:|---:|---:|
| Si | 8 | `-43.3066227505` | `-43.2957423300` | `-0.0108804205` | `-1.360` |
| HfO2 | 12 | `-121.8844407617` | `-122.0024819000` | `+0.1180411383` | `+9.837` |

The HfO2 absolute-energy difference is dominated by the position-independent
POTCAR atomic reference (`EATOM`) reconstruction. It shifts the total energy
but does not enter the Hamiltonian eigenproblem or atomic forces.

## Atomic force

| System | component MAE (meV/Angstrom) | component RMSE (meV/Angstrom) | maximum component error (meV/Angstrom) | HALF max force norm (meV/Angstrom) | VASP max force norm (meV/Angstrom) |
|---|---:|---:|---:|---:|---:|
| Si | `0.414` | `0.520` | `0.913` | `0.0158` | `1.096` |
| HfO2 | `10.988` | `14.566` | `28.154` | `95.31` | `122.19` |

This comparison includes both the remaining HALF functional error and the
difference between the DeepAW input density and the final self-consistent VASP
density. It is therefore the relevant end-to-end accuracy result.

## Consequence for VASP initialization

The workflow distinguishes three accuracy levels:

| Mode | Electronic update | Intended result |
|---|---|---|
| HALF direct | none | SCF-free Harris energy, bands, and analytic forces |
| one-step approximate SCF | VASP `ICHARG=1`, `NELM=1` | near-self-consistent energy after one electronic refinement |
| ACC-SCF | continue VASP to `EDIFF` | fully self-consistent VASP accuracy with fewer loops |

For the one-step approximate-SCF mode, the energy differences from fully
self-consistent VASP are `0.303 meV/atom` for Si and `2.492 meV/atom` for
HfO2. The corresponding force-component MAEs are `0.949` and
`27.986 meV/Angstrom`; one step is therefore an energy approximation, while
ACC-SCF is the higher-accuracy route for converged density and forces.

`EATOM`, Ewald bookkeeping, atomic PAW double counting, and the analytic force
response are evaluated after the eigenproblem. Improving them changes reported
energies or forces but does not change $H$, $S$, eigenvectors, the wavefunctions
passed to VASP, or the number of VASP SCF iterations.

Only changes on the wavefunction-generation path can reduce SCF loops:

1. a DeepAW density closer to the self-consistent density;
2. closer VASP parity in the local/NLCC potential and PAW SETDIJ operator;
3. a smaller occupied-subspace residual and better all-band convergence;
4. consistent plane-wave ordering, k points, occupancies, and basis cutoff at
   the direct memory handoff.

The new POTCAR interpolation is on this path because it enters
$V_{\mathrm{eff}}$; it can slightly improve the initial subspace. The completed
Harris-force response and future `EATOM` correction are not on this path and
cannot by themselves reduce SCF loops. Any loop-count claim therefore requires
a controlled MP-85 A/B rerun with identical VASP inputs.

## Data

- HALF: `/data/limusen/deepaw_half_force_compare/final_vasp_interp_gpu`
- self-consistent Si:
  `/data/limusen/deepawchgs/020_mp-149_Si`
- self-consistent HfO2:
  `/data/limusen/deepawchgs_retry_encut521/054_mp-352_HfO2`
- machine-readable summary:
  [`half_vs_vasp_scf_energy_force.json`](half_vs_vasp_scf_energy_force.json)
