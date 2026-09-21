# DeePAW-HALF: direct band structures from DeepAW and accelerated VASP SCF

## Abstract

DeePAW-HALF converts a smooth electron density predicted by DeepAW into a
plane-wave electronic structure. It provides two related but distinct
capabilities:

1. **DeepAW → HALF → band structure**: HALF reconstructs the fixed-density PAW
   Hamiltonian from the smooth density and matching POTCAR, and computes band
   structures, gaps, eigenvalues, energies, and wavefunctions without starting
   VASP.
2. **DeepAW → HALF → VASP SCF**: HALF writes the fixed-density eigenstates
   directly into VASP's wavefunction memory, replacing SAD/random initial
   orbitals. VASP then follows its normal `ALGO=All` SCF path with fewer
   electronic iterations.

For HfO2, the standalone HALF bands have a **0.840 meV** mean absolute
difference from VASP and the sampled gap differs by **0.557 meV**. In VASP
SCF, HALF initialization reduces the electronic iterations from 16 to 4--5
and gives an end-to-end speedup of **1.38--1.80x**. Across an 85-material
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
                         +-> HALF initial waves -> VASP SCF -> energy / forces
```

HALF is not another machine-learned band model. It reads a DeepAW smooth
`CHGCAR` or eSCN API density and a VASP PAW `POTCAR` with exactly matching
elements and datasets, then reconstructs the same fixed-density PAW/MIMIC_US
operator used by HAPPY.

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

### 2.4 Numerical parity and speed versus Python HAPPY

The HfO2 Gamma/MIMIC_US benchmark uses 520 eV, 3407 plane waves, and 60 bands.
All paths use the same CHGCAR, POTCAR, PBE functional, and complex128 precision.

| Implementation | Hardware/parallelism | Time | Relative to HAPPY |
|---|---|---:|---:|
| HAPPY Python | one logical CPU | 64.23 s | 1.00x |
| HALF CPU Fortran | one logical CPU | 43.790 s | 1.47x |
| HALF CUDA Fortran | RTX PRO 6000 | 1.696 s | **37.87x** |

The 60 HALF CUDA eigenvalues differ from HAPPY by at most
$7.80\times10^{-12}$ eV, with an RMS difference of
$2.31\times10^{-12}$ eV. CUDA acceleration therefore preserves the HAPPY
numerical model. This is a single-k-point operator benchmark, not an
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

#### `LMAXPAW=-1`

| Metric | DeePAW-HALF | VASP SAD | HALF effect |
|---|---:|---:|---:|
| SCF loops | **4** | 16 | 75% fewer; 4.0x iteration reduction |
| End-to-end wall time | **193.44 s** | 266.96 s | **1.38x**; 27.54% less |
| Final TOTEN | -121.101882671 eV | -121.101878858 eV | 3.813 micro-eV/cell |

#### Default VASP LMAX

| Metric | DeePAW-HALF | VASP SAD | HALF effect |
|---|---:|---:|---:|
| SCF loops | **5** | 16 | 68.75% fewer; 3.2x iteration reduction |
| End-to-end wall time | **145.76 s** | 262.70 s | **1.80x**; 44.51% less |
| Final TOTEN | -121.040359303 eV | -121.040376497 eV | 17.194 micro-eV/cell |

These results show that HALF changes the route to convergence rather than the
final VASP solution. Single-run wall times across the two different LMAX
settings are affected by node load and cache state, so the cross-setting
ablation should emphasize iteration counts. HALF and SAD wall times within
the same setting are the corresponding end-to-end comparisons.

An independent Si eSCN validation shows the same trend: HALF-generated waves
required four SCF loops, versus 12 for the original SAD baseline. A
same-machine SAD rerun took 11 loops, and `LOOP+` time fell from 2.4036 to
0.9251 seconds, a 2.60x acceleration of the SCF section.

### 3.4 MP-85 benchmark: average SCF acceleration

The broader benchmark covers 85 Materials Project structures and uses
precomputed DeepAW CHGCAR densities, `LHALF_API=.FALSE.`, `EDIFF=1E-4`,
`KSPACING=0.35`, `ALGO=All`, `LMAXMIX=-1`, and ISPIN=1. Both HALF and SAD
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
| Dense `VASP_LIKE` | 1861.047 s | 5 | -63.791680 eV |
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

## 4. Relationship between the two capabilities

| Item | Direct band structure | VASP SCF acceleration |
|---|---|---|
| Input | DeepAW density + POTCAR + k path | DeepAW density + POTCAR + VASP memory metadata |
| HALF output | Fixed-density eigenvalues and vectors | Initial VASP eigensubspace |
| Enters VASP | No | Yes |
| Self-consistent | No; one-shot Harris density | Yes; completed by VASP |
| Primary value | Screening, bands, and gaps | Fewer SCF loops and lower wall time |
| Large-system solver | ACC or k-point MPI | ACC followed by VASP `ALGO=All` |

Both paths share the same CHGCAR/POTCAR parser, effective potential, PAW
projectors, MIMIC_US terms, basis construction, and $H\Psi/S\Psi$
implementation. Direct-band validation therefore also provides the numerical
foundation for VASP initial wavefunctions.

## 5. Scope and limitations

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

## 6. Conclusion

DeePAW-HALF now provides two usable routes from a learned density to
plane-wave electronic structure:

- **Direct calculation:** HALF produces high-symmetry band structures from a
  DeepAW density without VASP. HfO2 bands have a 0.840 meV MAE relative to
  VASP, while Gamma/MIMIC_US eigenvalues agree with HAPPY at about
  $10^{-11}$ eV.
- **SCF acceleration:** HALF injects environment-aware initial wavefunctions
  directly into VASP. HfO2 SCF loops fall from 16 to 4--5; across MP-85, mean
  loops fall from 25.365 to 11.435 (2.218x); and ACC reduces the large-basis
  CsPbBr3 job time by 6.75x relative to dense HALF initialization.

HALF is therefore more than a Fortran/CUDA reproduction of HAPPY. It is a
reusable density-to-wavefunction layer: DeepAW supplies density, while HALF
either produces bands directly or gives VASP a substantially better SCF
starting point.

## 7. Data and reproducibility records

- HfO2 standalone bands:
  [`assets/hfo2_half_direct_bs.json`](assets/hfo2_half_direct_bs.json)
- HALF/VASP pointwise comparison:
  [`assets/hfo2_half_vasp_band_comparison.json`](assets/hfo2_half_vasp_band_comparison.json)
- HfO2 HALF/SAD SCF record:
  [`hfo2_vasp_half_vs_sad_kspacing035.json`](hfo2_vasp_half_vs_sad_kspacing035.json)
- HfO2 CUDA/HAPPY parity and performance:
  [`hfo2_cuda_uspp_parity.json`](hfo2_cuda_uspp_parity.json)
- CsPbBr3 ACC/VASP record:
  [`cspbbr3_vasp_acc_pro6000.json`](cspbbr3_vasp_acc_pro6000.json)
- Si eSCN/HALF/SAD record:
  [`si_escn_half_vs_sad_ediff1e4_lmaxmixm1.json`](si_escn_half_vs_sad_ediff1e4_lmaxmixm1.json)
- MP-85 aggregate HALF/SAD record:
  [`mp85_deepaw_half_vs_sad_ediff1e4_lmaxmixm1.json`](mp85_deepaw_half_vs_sad_ediff1e4_lmaxmixm1.json)

中文版：
[`HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md`](HFO2_VASP_HALF_VS_SAD_KSPACING035.zh-CN.md).
