# DeePAW-HALF: band structures, energies, analytic-force development, and accelerated VASP SCF

## Abstract

DeePAW-HALF converts a smooth electron density predicted by DeepAW into a
plane-wave electronic structure. It provides three related
capabilities:

1. **DeepAW → HALF → band structure**: HALF reconstructs the fixed-density PAW
   Hamiltonian from the smooth density and matching POTCAR, and computes band
   structures, gaps, eigenvalues, energies, and wavefunctions without starting
   VASP.
2. **DeepAW → HALF → VASP SCF**: HALF writes the fixed-density eigenstates
   directly into VASP's wavefunction memory, replacing SAD/random initial
   orbitals. VASP then follows its normal `ALGO=All` SCF path with fewer
   electronic iterations.
3. **DeepAW → HALF → energy and analytic-force development**: HALF evaluates
   the Harris total energy directly from the learned density. Its native
   Hellmann--Feynman implementation currently contains the Ewald, reciprocal
   local, and generalized PAW `D-epsilon Q` derivatives. The former
   finite-difference result is retained only as a derivative oracle; it is not
   a production total-force implementation.

For HfO2, the standalone HALF bands have a **0.840 meV** mean absolute
difference from VASP and the sampled gap differs by **0.557 meV**. For Si,
the CUDA and Python implementations agree within
$3.425\times10^{-12}$ eV/cell. The
previous $2.238\times10^{-9}$ eV/Angstrom finite-difference comparison and
**45.31x** timing are validation-oracle results, not claims for the unfinished
analytic total force. In VASP SCF, HALF initialization reduces the electronic iterations from 16 to 4
and gives an end-to-end speedup of **1.38x** in the primary comparison. Across an 85-material
benchmark, mean SCF iterations fall from 25.365 to 11.435, a **2.218x mean-loop
speedup**. For a larger 20-atom CsPbBr3 case, the matrix-free ACC solver
preserves five SCF loops and the same final energy while reducing total time
from 1861.0 to 275.5 seconds, a **6.75x speedup**.

## 1. Role of DeePAW-HALF

DeepAW primarily predicts a smooth valence density, whereas conventional band
structures and VASP SCF calculations require wavefunctions. HALF is the
deterministic physics layer that maps density to wavefunctions:

```text
                         +-> HALF fixed-density solve -> bands / gap / energy
structure -> DeepAW rho -+
                         +-> HALF Harris energy -> analytic-force development
                         +-> HALF initial waves -> VASP SCF -> converged properties
```

HALF is not another machine-learned band model. It reads a DeepAW smooth
`CHGCAR` or eSCN API density and a VASP PAW `POTCAR` with exactly matching
elements and datasets, then reconstructs the DeePAW-HALF fixed-density
PAW/MIMIC_US operator.

For a DeepAW density $\widetilde\rho_0$, HALF constructs the Hamiltonian once,

$$
\widehat H_0=\widehat H[\widetilde\rho_0],
$$

and solves the generalized PAW eigenproblem

$$
H(\mathbf k)C(\mathbf k)
=S(\mathbf k)C(\mathbf k)\varepsilon(\mathbf k).
$$

In plane-wave form,

$$
H=T+V_{\mathrm{eff}}+BDB^{\dagger},
$$

$$
S=I+BQB^{\dagger}.
$$

$B$ is the PAW projector matrix, $D$ includes DION and the
potential-dependent QDEP correction, and $Q$ is the PAW overlap augmentation.
The complete derivation, matrix-free $H\Psi/S\Psi$, S-metric
orthogonalization, and all-band Rayleigh-Ritz algorithm are documented in
[Matrix-free Harris all-band acceleration](../HARRIS_ACC_THEORY.md).

## 2. Capability one: DeepAW directly to band structure

### 2.1 Workflow

The direct band path does not enter VASP:

```text
POSCAR / structure
  -> DeepAW smooth density
  -> CHGCAR or in-memory density + POTCAR
  -> HALF Veff, PAW projectors, D and Q
  -> solve H(k)C=S(k)C epsilon along a high-symmetry path
  -> JSON / CSV / NPZ / PNG / vaspwave.h5
```

HALF supports explicit `KPOINTS`, automatic high-symmetry paths, arbitrary
single k points, and uniform meshes. The CPU build can distribute independent
k points with MPI. CUDA dense solvers serve small and medium bases, while the
matrix-free ACC solver targets $N_{\mathrm{PL}}\gg N_{\mathrm{BANDS}}$.

### 2.2 Standalone HfO2 band structure

The validation system contains 12 HfO2 atoms. The standalone HALF run used
PBE, 500 eV, MIMIC_US, 60 bands, and 100 high-symmetry-path samples. Energies
are aligned to the valence-band maximum (VBM).

![HfO2 band structure computed directly by DeePAW-HALF](assets/hfo2_half_direct_bs.png)

HALF gives:

- sampled band gap: **4.576341 eV**;
- direct Gamma-point gap: **4.629091 eV**.

### 2.3 Point-by-point comparison with VASP

VASP 6.6.0 and HALF use the same 100 k points and 60 bands. Each spectrum is
aligned to its own VBM.

![HfO2 band-structure comparison between HALF and VASP](assets/hfo2_half_vasp_band_comparison.png)

| Metric | DeePAW-HALF | VASP 6.6.0 | Difference |
|---|---:|---:|---:|
| Sampled gap | 4.576341 eV | 4.575784 eV | 0.557 meV |
| Direct Gamma gap | 4.629091 eV | 4.631752 eV | 2.661 meV |
| Band MAE | -- | -- | 0.840 meV |
| Band RMSE | -- | -- | 1.135 meV |
| Maximum pointwise band difference | -- | -- | 3.250 meV |

HALF therefore reproduces the HfO2 dispersion and gap without entering the
VASP executable. This validates electronic-structure reconstruction for a
given smooth density. Ultimate physical accuracy also depends on the DeepAW
density, POTCAR compatibility, XC functional, and basis settings.

### 2.4 Numerical parity and implementation speed

The HfO2 Gamma/MIMIC_US benchmark uses 520 eV, 3407 plane waves, and 60 bands.
All paths use the same CHGCAR, POTCAR, PBE functional, and complex128 precision.

| DeePAW-HALF implementation | Hardware/parallelism | Time | Relative time |
|---|---|---:|---:|
| Python implementation | one logical CPU | 64.23 s | 1.00x |
| CPU Fortran implementation | one logical CPU | 43.790 s | 1.47x |
| CUDA Fortran implementation | RTX PRO 6000 | 1.696 s | **37.87x** |

The 60 CUDA eigenvalues differ from the Python implementation by at most
$7.80\times10^{-12}$ eV, with an RMS difference of
$2.31\times10^{-12}$ eV. CUDA acceleration therefore preserves the
DeePAW-HALF numerical model. This is a single-k-point operator benchmark, not an
end-to-end timing for the 100-point path.

### 2.5 CLI entry point

```bash
half bands CHGCAR.deepaw POTCAR KPOINTS \
  --encut 500 --bands 60 --backend cuda --solver acc --uspp-dij \
  --output-prefix hfo2_bands
```

This writes `hfo2_bands.json`, `.csv`, `.npz`, and `.png`. For a small system,
`--solver evx` can replace `--solver acc` to use the dense partial-spectrum
reference solver.

## 3. Capability two: DeePAW-HALF acceleration of VASP SCF

### 3.1 Why the initial wavefunctions matter

Native VASP SAD begins from a superposition of isolated-atom densities and
default initial orbitals. In systems with bonding, charge transfer, or complex
PAW channels, that initial subspace can be far from the converged low-energy
subspace. A DeepAW density already contains structural-environment information.
The Harris eigenstates reconstructed by HALF are consequently a more informed
starting point for VASP.

The returned HALF orbitals do not need final-SCF eigenstate accuracy. VASP
continues with its own `ORTHCH`, `PROALL`, `ALGO=All`, density update, and
mixing. ACC can therefore stop at a loose residual or a fixed iteration cap.

### 3.2 Direct in-memory handoff

The adapted VASP does not use WAVECAR as an intermediate:

1. VASP passes lattice vectors, fractional positions, species, FFT grid,
   ENCUT, and NBANDS to libhalf.
2. At every k point, VASP passes its own ordered $(G_x,G_y,G_z)$ list.
3. HALF verifies exact basis equality and returns complex128 coefficients
   already permuted into VASP order.
4. VASP distributes global bands through `DIS_PW_BAND`, fills `W%CELTOT`, and
   skips random `WFINIT`.
5. VASP starts its unmodified SCF path from those orbitals.

The INCAR selection is

```text
LHALF_INIT = .TRUE.
LHALF_API  = .FALSE.  # use an existing DeepAW CHGCAR
HALF_MODE  = ACC
```

With `LHALF_API=.TRUE.`, HALF instead obtains density from the DeePAW-eSCN
service configured by `HALF_ESCN_URL`; the wavefunction handoff is unchanged.

### 3.3 HfO2: HALF initialization versus SAD

The VASP 6.6.0 full-complex `vasp_std` benchmark uses `KSPACING=0.35`, 36
irreducible k points, `PREC=High` (actual ENCUT 500 eV), `ALGO=All`,
`EDIFF=1E-4`, and ISPIN=1.

| Metric | DeePAW-HALF | VASP SAD | HALF effect |
|---|---:|---:|---:|
| SCF loops | **4** | 16 | 75% fewer; 4.0x iteration reduction |
| End-to-end wall time | **193.44 s** | 266.96 s | **1.38x**; 27.54% less |
| Final TOTEN | -121.101882671 eV | -121.101878858 eV | 3.813 micro-eV/cell |

The micro-eV final-energy difference shows that HALF changes the route to
convergence rather than the final VASP solution.

An independent Si eSCN validation shows the same trend: HALF-generated waves
required four SCF loops, versus 12 for the original SAD baseline. A
same-machine SAD rerun took 11 loops, and `LOOP+` time fell from 2.4036 to
0.9251 seconds, a 2.60x acceleration of the SCF section.

### 3.4 MP-85 benchmark: average SCF acceleration

The broader benchmark covers 85 Materials Project structures and uses
precomputed DeepAW CHGCAR densities, `LHALF_API=.FALSE.`, `EDIFF=1E-4`,
`KSPACING=0.35`, `ALGO=All`, and ISPIN=1. Both HALF and SAD
sets reached EDIFF and terminated normally for all 85 structures. Every HALF
`vasp.out` contains the initialization banner and initialized-k-point records.

| Aggregate over 85 materials | DeePAW-HALF | VASP SAD | HALF effect |
|---|---:|---:|---:|
| Total SCF loops | 972 | 2156 | 1184 fewer |
| Mean SCF loops | **11.435** | 25.365 | **54.92% fewer** |
| Mean-loop speedup | -- | -- | **2.218x** |

Here the reported average speedup is defined explicitly as

$$
S_{\mathrm{mean}}=
\frac{\langle N_{\mathrm{SCF}}^{\mathrm{SAD}}\rangle}
{\langle N_{\mathrm{SCF}}^{\mathrm{HALF}}\rangle}
=\frac{25.3647}{11.4353}=2.2181.
$$

This statistic measures iteration acceleration, not wall time. It avoids
mixing material-dependent system sizes and makes the initialization effect
comparable across the full set. HfO2 and XeF2 use the separately validated
521 eV retry because their VASP/HALF basis boundary differed at exactly
520 eV.

### 3.5 Why large bases also require ACC

Better initial orbitals reduce SCF iterations, but dense
$N_{\mathrm{PL}}\times N_{\mathrm{PL}}$ H/S construction can make HALF
initialization itself the bottleneck. ACC uses block $H\Psi/S\Psi$ and a
globally coupled all-band Rayleigh-Ritz step. cuSOLVER sees only
$N_{\mathrm{BANDS}}$ or $2N_{\mathrm{BANDS}}$ projected matrices.

The 20-atom CsPbBr3 test contains 12 k points, 105 bands, and
$N_{\mathrm{PL}}=20640$--$20780$. At $N_{\mathrm{PL}}=20767$, two dense
complex128 H/S matrices alone require about 13.8 GB.

| HALF initialization | Total VASP time | SCF loops | Final E0 |
|---|---:|---:|---:|
| Dense baseline | 1861.047 s | 5 | -63.791680 eV |
| Matrix-free `ACC` | **275.522 s** | 5 | -63.791680 eV |

ACC gives a **6.75x total speedup**. After subtracting VASP's `LOOP+`, HALF
initialization plus other peripheral phases falls from about 1573.7 to 53.8
seconds, an estimated **29.3x speedup**. This subtraction is not an
independently instrumented HALF-only timer, but it clearly exposes removal of
the dense-initialization bottleneck.

After 40 ACC iterations, the per-k-point maximum residuals remain
0.022--0.410 eV. Nevertheless, VASP still uses five SCF loops and reaches the
identical final E0. This directly validates the design choice that an initial
subspace need not be converged to final-SCF eigenstate accuracy.

## 4. Capability three: direct energy and force evaluation

### 4.1 Harris total energy

After solving the occupied fixed-density eigenstates, HALF evaluates

$$
E_{\mathrm{HALF}}=
\sum_{n\mathbf{k}}w_{\mathbf{k}}f_{n\mathbf{k}}\varepsilon_{n\mathbf{k}}
-E_{\mathrm H}
-\int \widetilde\rho(\mathbf r)v_{\mathrm{xc}}(\mathbf r)\,d\mathbf r
+E_{\mathrm{xc}}
+E_{\mathrm{Ewald}}
+E_{G=0}
+E_{\mathrm{atom}}
+E_{\mathrm{PAW}}.
$$

The implementation includes occupations and entropy, Hartree and XC
double-counting corrections, Ewald ion-ion energy, local-potential $G=0$ and
atomic reference terms. The spherical atomic PAW correction is reconstructed
from POTCAR atomic occupancies, AE/PS partial waves, core density, `DEXC`, and
compensation charge. It is not read from a CHGCAR augmentation tail or copied
from VASP output. An explicit k-point set or a symmetry-reduced mesh can be
used, with MPI distributing independent k points in the CPU build.

### 4.2 Atomic forces

HALF evaluates production forces analytically. In compact form,

$$
F_{I\alpha}=-\frac{\partial E_{\mathrm{HALF}}}{\partial R_{I\alpha}}
=F^{\mathrm{Ewald}}+F^{\mathrm{local}}+F^{D-\varepsilon Q}
+F^{\mathrm{aug}}+F^{\mathrm{NLCC}}+F^{\mathrm{Harris\ response}}.
$$

The generalized nonlocal term is evaluated as

$$
F^{D-\varepsilon Q}_{I\alpha}
=-\sum_{n\mathbf k}w_{\mathbf k}f_{n\mathbf k}
\left\langle\psi_{n\mathbf k}\left|
\partial_{I\alpha}H-\varepsilon_{n\mathbf k}\partial_{I\alpha}S
\right|\psi_{n\mathbf k}\right\rangle .
$$

The DeepAW smooth input density is frozen. The Harris convergence force uses
the full response

$$
g_{\mathrm{Hxc}}=v_H[n_{\mathrm{out}}-n_{\mathrm{in}}]
+f_{\mathrm{xc}}[n_{\mathrm{in}}+n_{\mathrm{core}}]
(n_{\mathrm{out}}-n_{\mathrm{in}}),
$$

contracted with the translated POTCAR `PSPRHO`; NLCC separately differentiates
`PSPCOR`. Central finite differences are only a developer oracle and never the
production CLI. A single call produces energy components and analytic forces:

```bash
half energy CHGCAR.deepaw POTCAR \
  --encut 520 --kspacing 0.35 --bands 24 --backend cuda \
  --forces --output-prefix si_energy_force
```

Against the strict frozen-density VASP `ICHARG=11` oracle, the current
force-component MAE is `6.99e-6 eV/Angstrom` for Si and
`1.007e-3 eV/Angstrom` for HfO2. This is a like-for-like comparison in which
VASP retains its Harris convergence correction; a one-step `ICHARG=1` run does
not. Full formulas, $L$-channel diagnostics, and raw-data locations are in
[`HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.md`](HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.md).

### 4.3 Comparison with reported machine-learning potentials

DeepAW-HALF and universal machine-learning potentials share the objective of
using learned information to replace or substantially reduce conventional
SCF-DFT work. The following table places HALF's current energy/force validation
beside representative published results.

| Method | Reported energy result | Reported force result |
|---|---:|---:|
| DeePAW-HALF, Si/HfO2 current validation | 2.59/7.17 meV/atom absolute energy difference vs initial VASP MIMIC_US state | 0.00699/1.007 meV/Angstrom component MAE vs fixed-density VASP `ICHARG=11` |
| [M3GNet](https://doi.org/10.1038/s43588-022-00349-3) | 35 meV/atom MAE | 72 meV/Angstrom MAE |
| [CHGNet](https://doi.org/10.1038/s42256-023-00716-3) | 30 meV/atom MAE | 77 meV/Angstrom MAE |
| [MACE-MP-0 medium](https://doi.org/10.1063/5.0297006) | 20 meV/atom MAE | 45 meV/Angstrom MAE |

These results position DeePAW-HALF as an electronic-structure route from a
learned density to bands, energy, forces, and wavefunctions, while retaining a
direct path into VASP when a fully self-consistent result is required.

## 5. Relationship between the three capabilities

| Item | Direct band structure | Direct energy and forces | VASP SCF acceleration |
|---|---|---|---|
| Input | DeepAW density + POTCAR + k path | DeepAW density + POTCAR + k mesh | DeepAW density + POTCAR + VASP memory metadata |
| HALF output | Fixed-density eigenvalues and vectors | Harris energy and atomic forces | Initial VASP eigensubspace |
| Enters VASP | No | No | Yes |
| Primary value | Screening, bands, and gaps | Energy/force evaluation without SCF | Fewer SCF loops and lower wall time |
| Large-system path | ACC or k-point MPI | ACC plus k-point MPI | ACC followed by VASP `ALGO=All` |

All three paths share the same CHGCAR/POTCAR parser, effective potential, PAW
projectors, MIMIC_US terms, basis construction, and $H\Psi/S\Psi$
implementation. Direct-band validation therefore also provides the numerical
foundation for VASP initial wavefunctions.

## 6. Scope and limitations

- HALF computes fixed-density electronic structure for a supplied DeepAW
  density. Physical accuracy depends jointly on the density model, functional,
  POTCAR, and basis.
- CHGCAR and POTCAR must have consistent element order, PAW dataset versions,
  and valence definitions.
- The current VASP adapter supports ISPIN=1, no SOC/noncollinear spin,
  full-complex storage, and KPAR=1. These are adapter restrictions, not
  limitations of the PAW/ACC equations.
- HfO2 and CsPbBr3 wall times are single bare-metal observations for each
  corresponding path. Iteration counts, final energies, and eigenvalue errors
  are more stable than wall times measured in different periods.
- ACC tolerance and iteration limits should follow the task: standalone bands
  need tighter convergence, while VASP initialization can stop earlier.
- Central finite differences are retained only as a validation oracle. Native
  analytic Ewald, local-potential, generalized nonlocal PAW, augmentation,
  NLCC, and full frozen-density Hartree-plus-XC Harris response are
  implemented. Against fixed-density VASP `ICHARG=11`, the Si/HfO2
  force-component MAEs are `6.99e-6` and `1.007e-3 eV/Angstrom`.

## 7. Conclusion

DeePAW-HALF now provides three usable routes from a learned density to
plane-wave electronic structure:

- **Direct calculation:** HALF produces high-symmetry band structures from a
  DeepAW density without VASP. HfO2 bands have a 0.840 meV MAE relative to
  VASP, while the independent Gamma/MIMIC_US implementations agree at about
  $10^{-11}$ eV.
- **SCF acceleration:** HALF injects environment-aware initial wavefunctions
  directly into VASP. HfO2 SCF loops fall from 16 to 4; across MP-85, mean
  loops fall from 25.365 to 11.435 (2.218x); and ACC reduces the large-basis
  CsPbBr3 job time by 6.75x relative to dense HALF initialization.
- **Energy and analytic forces:** HALF evaluates the Harris total energy
  without a preceding SCF run. All production force components are analytic
  and use POTCAR-derived augmentation. Against fixed-density VASP `ICHARG=11`,
  the Si/HfO2 force-component MAEs are `6.99e-6` and
  `1.007e-3 eV/Angstrom`, respectively.

DeePAW-HALF is a reusable density-to-electronic-structure layer: DeepAW
supplies density, while HALF produces bands, energies, forces, and wavefunctions
directly, or gives
VASP a substantially better SCF starting point.

## 8. Data and reproducibility records

- HfO2 standalone bands:
  [`assets/hfo2_half_direct_bs.json`](assets/hfo2_half_direct_bs.json)
- HALF/VASP pointwise comparison:
  [`assets/hfo2_half_vasp_band_comparison.json`](assets/hfo2_half_vasp_band_comparison.json)
- HfO2 HALF/SAD SCF record:
  [`hfo2_vasp_half_vs_sad_kspacing035.json`](hfo2_vasp_half_vs_sad_kspacing035.json)
- HfO2 cross-implementation parity and performance:
  [`hfo2_cuda_uspp_parity.json`](hfo2_cuda_uspp_parity.json)
- CsPbBr3 ACC/VASP record:
  [`cspbbr3_vasp_acc_pro6000.json`](cspbbr3_vasp_acc_pro6000.json)
- MP-85 aggregate HALF/SAD record:
  [`mp85_deepaw_half_vs_sad_ediff1e4.json`](mp85_deepaw_half_vs_sad_ediff1e4.json)
- Si total-energy parity:
  [`si_total_energy_parity.json`](si_total_energy_parity.json)
- Si finite-difference-force parity and performance:
  [`si_force_parity.json`](si_force_parity.json)

中文版：
[`HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md`](HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md).
