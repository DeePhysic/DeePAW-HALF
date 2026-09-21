# Matrix-free Harris all-band acceleration

This document describes the theory and implementation of the CUDA `acc`
solver in DeePAW-HALF and its direct VASP 6.6 integration. It is intended to
make the numerical contract auditable: every large operation below has a
corresponding implementation in `half_cuda_assembly.cuf` or
`half_cuda_iterative_solver.cuf`.

## 1. Fixed-density Harris problem

DeepAW supplies a smooth valence density, denoted by
$\widetilde{\rho}_0(\mathbf r)$. HALF does not update that density. It builds a
single fixed Hamiltonian

$$
\widehat H_0 = \widehat H[\widetilde{\rho}_0]
$$

and obtains initial orbitals from the generalized PAW eigenproblem

$$
\widehat H_0 |\widetilde\psi_n\rangle
= \varepsilon_n \widehat S |\widetilde\psi_n\rangle.
$$

This differs from a VASP SCF calculation. HALF solves the one-shot Harris
problem; VASP subsequently uses the returned orbitals as its initial subspace
and continues to update density and orbitals normally.

The fixed effective local potential is

$$
V_{\mathrm{eff}}(\mathbf r)
= V_{\mathrm{ion}}^{\mathrm{loc}}(\mathbf r)
+ V_H[\widetilde\rho_0](\mathbf r)
+ V_{\mathrm{xc}}[\widetilde\rho_0+\rho_{\mathrm{core}}](\mathbf r).
$$

The partial-core density enters only the nonlinear core correction in the XC
term. It is not added to the Hartree density.

## 2. PAW generalized eigenproblem

The PAW transformation from a smooth pseudo wavefunction to its all-electron
counterpart is

$$
\widehat{\mathcal T}
= 1 + \sum_{Ii}
\left(|\phi_i^I\rangle-|\widetilde\phi_i^I\rangle\right)
\langle\widetilde p_i^I|.
$$

Consequently, the pseudo-space overlap operator is not the identity:

$$
\widehat S
= \widehat{\mathcal T}^{\dagger}\widehat{\mathcal T}
= 1 + \sum_{Iij}
|\widetilde p_i^I\rangle Q_{ij}^I
\langle\widetilde p_j^I|.
$$

For a plane-wave basis $|\mathbf k+\mathbf G\rangle$, define the projector
matrix

$$
B_{\mathbf G,Ii}
= \langle\mathbf k+\mathbf G|\widetilde p_i^I\rangle.
$$

Let $T$ be the diagonal plane-wave kinetic matrix and let $V$ denote the local
potential convolution. With block-diagonal PAW matrices $D$ and $Q$, the dense
operators would be

$$
H = T + V + B D B^{\dagger},
$$

$$
S = I + B Q B^{\dagger}.
$$

For the MIMIC_US path, the onsite Hamiltonian is atom dependent:

$$
D_{ij}^{I}
= D_{ij}^{I,\mathrm{ION}}
+ \int V_{\mathrm{eff}}(\mathbf r)
Q_{ij}^{I,\mathrm{DEP}}(\mathbf r-\mathbf R_I)\,d^3r.
$$

The dense `evd`, `evj`, and `evx` solvers explicitly build $H$ and $S$. The
`acc` solver never constructs either $N_{\mathrm{PL}}\times N_{\mathrm{PL}}$
matrix.

## 3. Why all-band Rayleigh-Ritz is the correct formulation

Collect the requested $m=N_{\mathrm{BANDS}}$ orbitals as columns of $X$. The
fixed-density occupied-subspace problem can be written as a constrained trace
minimization:

$$
\min_X \operatorname{Tr}(X^{\dagger} H X)
\quad\text{subject to}\quad
X^{\dagger} S X = I_m.
$$

Introduce a Hermitian matrix of Lagrange multipliers $\Lambda$:

$$
\mathcal L(X,\Lambda)
= \operatorname{Tr}(X^{\dagger} H X)
- \operatorname{Tr}\left[\Lambda(X^{\dagger}S X-I_m)\right].
$$

Stationarity with respect to $X^{\dagger}$ gives

$$
H X = S X \Lambda.
$$

Therefore the bands are coupled through one common variational subspace. A
band block is an execution unit for FFTs and contractions, not an independent
eigenproblem. This is the same essential distinction as VASP `ALGO=All`:
execution is blocked, while optimization and orthogonality remain global over
all requested bands.

## 4. Matrix-free application of H and S

For an input block $X_b$ containing a few bands, HALF evaluates

$$
Y_b = B^{\dagger}X_b,
$$

$$
H X_b
= T X_b
+ \mathcal F\!\left[
V_{\mathrm{eff}}(\mathbf r)
\mathcal F^{-1}[X_b]
\right]
+ B D Y_b,
$$

$$
S X_b = X_b + B Q Y_b.
$$

Here $\mathcal F$ and $\mathcal F^{-1}$ include the FFT normalization used by
the plane-wave representation. The implementation performs the following work
without returning wavefunction-sized arrays to the CPU:

1. scatter plane-wave coefficients into the FFT grid;
2. inverse cuFFT to real space;
3. multiply by $V_{\mathrm{eff}}(\mathbf r)$ with a CUDA kernel;
4. forward cuFFT and gather back to plane-wave order;
5. add the diagonal kinetic term;
6. evaluate $B^{\dagger}X_b$, $BDY_b$, and $BQY_b$ with cuBLAS and CUDA
   projector kernels.

`--acc-block-size` controls only the number of bands in one H/S application
batch. It does not remove coupling between different batches.

## 5. Restarted block iteration

### 5.1 Initial subspace

HALF starts from unit vectors associated with the $m$ lowest kinetic-energy
plane waves. At arbitrary VASP k points these indices are obtained by an
$O(N_{\mathrm{PL}}\log N_{\mathrm{PL}})$ heap sort. The previous insertion
sort was unacceptable for large shifted-k bases because its worst-case cost
was $O(N_{\mathrm{PL}}^2)$.

After applying H and S, an initial Rayleigh-Ritz solve produces an
S-orthonormal Ritz basis.

### 5.2 Residual

For Ritz vectors $X$ and diagonal Ritz values $\Lambda$, the generalized
residual is

$$
R = H X - S X\Lambda.
$$

The reported convergence number is

$$
r_{\max} = \max_n \|R_n\|_2.
$$

Because each Ritz vector is normalized in the S metric, $r_{\max}$ has energy
units and is a direct measure of eigen-equation error. `--acc-tol` is the
target value; `--acc-max-iter` is a hard limit and may stop first.

### 5.3 Kinetic preconditioner

The plane-wave kinetic diagonal supplies an inexpensive approximation to the
inverse of $H-\varepsilon_n S$. HALF uses

$$
P_{\mathbf G n}
= -\frac{R_{\mathbf G n}}
{\left|T_{\mathbf G}-\varepsilon_n\right|+\delta},
$$

with $\delta=1\ \mathrm{eV}$ in the current implementation. This damps
high-kinetic-energy residual components without forming or factoring H.

### 5.4 S-orthogonalization

The correction block must be orthogonal to the current Ritz space in the PAW
metric. First compute

$$
C = X^{\dagger} S P,
$$

then update the correction and its already-computed operator images:

$$
P \leftarrow P-XC,
$$

$$
HP \leftarrow HP-(HX)C,
$$

$$
SP \leftarrow SP-(SX)C.
$$

Next form the correction metric

$$
M=P^{\dagger}SP.
$$

If


$$
M=U\,\operatorname{diag}(\mu)\,U^{\dagger},
$$

then the normalized correction is

$$
P \leftarrow P U\,\operatorname{diag}(\mu^{-1/2}).
$$

The same right transformation is applied to $HP$ and $SP$. Only this
$m\times m$ metric eigensolve uses cuSOLVER.

### 5.5 Global Rayleigh-Ritz update

Construct the restarted trial space and its operator images:

$$
Z=[X\;P],\qquad HZ=[HX\;HP],\qquad SZ=[SX\;SP].
$$

Project the operators:

$$
H_Z=Z^{\dagger}HZ,
$$

$$
S_Z=Z^{\dagger}SZ.
$$

Solve the small generalized problem

$$
H_ZY=S_ZY\Theta
$$

with cuSOLVER, retain its lowest $m$ columns $Y_m$, and rotate all bands
together:

$$
X\leftarrow ZY_m,\qquad
HX\leftarrow HZY_m,\qquad
SX\leftarrow SZY_m.
$$

The subspace dimension never exceeds $2m$. Thus cuSOLVER handles only
$m\times m$ and $2m\times2m$ matrices; it never receives a dense
$N_{\mathrm{PL}}\times N_{\mathrm{PL}}$ Hamiltonian in `acc` mode.

## 6. Early stopping for VASP initialization

A standalone band-structure calculation needs tightly converged eigenpairs.
VASP initialization does not: it needs a subspace with useful overlap with the
low-energy eigenspace. VASP subsequently performs its own orthogonalization,
all-band optimization, density update, and SCF convergence.

For that reason HALF permits either a loose residual tolerance or termination
at the iteration limit. This is an explicit accuracy/performance tradeoff, not
a claim that the returned states are final SCF orbitals. In the validated
20-atom CsPbBr3 case, 40 ACC iterations ended with per-k-point maximum
residuals from 0.022 to 0.410 eV. Nevertheless, ACC and dense VASP_LIKE both
required five VASP SCF loops and reached the same reported final energy.

## 7. Complexity and memory

Let $N=N_{\mathrm{PL}}$, $m=N_{\mathrm{BANDS}}$, $N_r$ be the real-space FFT
grid size, $b$ the H/S execution block size, and $N_p$ the total number of PAW
projectors.

The dense path stores at least H and S:

$$
M_{\mathrm{dense}} \ge 2\times16N^2\ \mathrm{bytes}
=32N^2\ \mathrm{bytes}
$$

for complex128 data, before solver workspaces and copies. Dense assembly is
$O(N^2)$ and the full generalized eigensolve is $O(N^3)$; index-range EVX
reduces eigenvector work but not dense H/S storage and assembly.

The ACC state storage is linear in the plane-wave count:

$$
M_{\mathrm{ACC}}=O(Nm+N_rb+NN_p+m^2).
$$

One matrix-free H/S application costs approximately

$$
O\!\left(mN_r\log N_r + NmN_p + Nm\right),
$$

while projected cuBLAS products cost $O(Nm^2)$ and the small Rayleigh-Ritz
solve costs $O(m^3)$. This replaces scaling in the full basis dimension with
scaling in the requested band count.

For the CsPbBr3 validation, $N$ was 20640--20780 and $m=105$. Two dense
complex128 matrices alone require about 13.8 GB at $N=20767$, whereas ACC does
not allocate them.

## 8. Direct VASP handoff

The adapted VASP source is retained under `vendor/vasp-6.6.0`. The integration
sequence is:

1. VASP reads INCAR, POSCAR, POTCAR, and the DeepAW smooth CHGCAR.
2. After `ALLOCW`, VASP passes lattice vectors, fractional positions, species,
   FFT-grid dimensions, ENCUT, and NBANDS directly to libhalf.
3. For every VASP k point, the adapter builds the exact serial
   $(G_x,G_y,G_z)$ list from `wavedes1`.
4. `half_solve_kpoint_mapped()` verifies exact basis equality, runs ACC, and
   returns coefficients already permuted into VASP's ordering.
5. VASP distributes every returned global band with `DIS_PW_BAND`, stores the
   Ritz values in `W%CELTOT`, skips random `WFINIT`, and continues through its
   normal `PROALL`, `ORTHCH`, and SCF path.

The relevant INCAR selection is

```text
LHALF_INIT = .TRUE.
LHALF_API  = .FALSE.
HALF_MODE  = ACC
```

`LHALF_API=.FALSE.` uses the existing DeepAW CHGCAR. Setting it to `.TRUE.`
selects the DeePAW-eSCN density API configured by `HALF_ESCN_URL`; the
wavefunction handoff is otherwise unchanged.

The current adapter requires ISPIN=1, no SOC/noncollinear spin, the
full-complex VASP build, and KPAR=1. These are adapter limitations, not
limitations of the matrix-free equations.

## 9. Implementation map

| Operation | Source |
|---|---|
| GPU fixed-density potential | `src/half_cuda_potential.cuf` |
| PAW/MIMIC_US data and matrix-free H/S | `src/half_cuda_assembly.cuf` |
| Residual, preconditioner, S-orthogonalization, restarted Ritz solve | `src/half_cuda_iterative_solver.cuf` |
| Dense/ACC dispatch | `src/half_cuda_solver.cuf` |
| C and Fortran solver selector | `include/half.h`, `include/half_api.f90` |
| CLI controls and JSON diagnostics | `app/half_cli.F90` |
| VASP memory adapter | `vendor/vasp-6.6.0/src/half_vasp_init.F` |
| VASP INCAR validation | `vendor/vasp-6.6.0/src/reader.F` |

Measured Pro 6000 results and the exact calculation contract are recorded in
[`validation/cspbbr3_vasp_acc_pro6000.json`](validation/cspbbr3_vasp_acc_pro6000.json).
