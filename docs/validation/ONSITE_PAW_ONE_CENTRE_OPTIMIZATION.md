# Onsite PAW / one-centre alignment and optimization

This note records the `force-accuracy-experiments` implementation that aligns
HALF's MIMIC_US onsite construction with the VASP 6.6.0 algorithm. It describes
the mathematics and validation without reproducing proprietary source.

## Direct-grid SETDIJ

For a frozen-density MIMIC_US calculation,

$$
H=H_{\mathrm{smooth}}+\sum_{Iij}|\widetilde p_i^I\rangle
D_{ij}^I\langle\widetilde p_j^I|,
\qquad
S=1+\sum_{Iij}|\widetilde p_i^I\rangle q_{ij}^I
\langle\widetilde p_j^I|,
$$

with

$$
D_{ij}^I=D_{ij}^{\mathrm{ion}}+
\sum_{LM}C_{ij}^{LM}Q_{ij}^{L}D_{LM}^{I}.
$$

The old HALF implementation interpolated the three-dimensional potential onto
radial shells and then performed radial and angular quadrature. The new default
uses the same discrete target as VASP SETDIJ:

$$
D_{LM}^{I}\simeq\frac{\Omega}{N_{\mathrm{grid}}}
\sum_{\mathbf r_g}V_{\mathrm{eff}}(\mathbf r_g)
g_L(r_{gI})Y_{LM}(\widehat{\mathbf r_{gI}}).
$$

Here $g_L=A_1j_L(q_1r)+A_2j_L(q_2r)$ is zero with zero first derivative at
the augmentation radius and is normalized by

$$
\int_0^{r_c}g_L(r)r^{L+2}\,dr=1.
$$

CPU and CUDA Fortran implement this direct FFT-grid contraction. On CUDA, one
thread block reduces each `(ion,LM)` pair, so the effective potential and
onsite contraction stay on the GPU. The legacy shell implementation remains
available with `-DHALF_QDEP_DIRECT_GRID=OFF` for controlled A/B tests.

## Shared one-centre radial moments

After POTCAR parsing, HALF now reconstructs every angular channel once:

$$
Q_{ij}^{L}=\int
\left(\phi_i\phi_j-\widetilde\phi_i\widetilde\phi_j\right)r^L\,dr.
$$

The weights follow Simpson integration in the logarithmic coordinate
$x=\ln r$, including $dr=r\,dx$. The resulting `qpaw_l(i,j,L)` is shared by
the overlap, SETDIJ, augmentation density, and atomic PAW double-counting.
The active $Q^0$ reconstruction differs from the tabulated Si and HfO2 values
by at most `2.48e-7` and `3.71e-7`, respectively.

## Force consistency and comparison

### First-stage result (historical)

The position derivative acts on the same compact grid kernel used for the
energy. `half-force-check` uses a displaced-ion central difference only as an
independent development oracle. The maximum errors in `dD/dR` are
`2.75e-9 eV/Angstrom` for Si and `2.43e-9 eV/Angstrom` for HfO2, four to five
orders of magnitude below the previous shell-quadrature errors.

| system | HALF free energy (eV) | HALF - VASP (eV/cell) | force-component MAE (eV/Angstrom) |
|---|---:|---:|---:|
| Si | -43.306625404 | -0.013307174 | 0.000739611 |
| HfO2 | -121.884427193 | +0.088155727 | 0.068379790 |

These values record the first stage after the direct-grid implementation. The
reference used an ordinary VASP `ICHARG=1`, `NELM=1` run. VASP regarded the
single exact diagonalization as electronically converged and therefore did not
add the frozen-density Harris convergence correction to its final force. It is
not the strict reference for the current HALF Harris force.

The analytic-force pass also caches the k-independent `D` and `dD/dR` arrays.
In the contended Pro 6000 run, HfO2 wall time decreased from `121.87` to
`97.78 s`; this is evidence that repeated construction was removed, not an
exclusive-GPU performance claim.

Machine-readable results are in
[`onsite_paw_direct_grid.json`](onsite_paw_direct_grid.json).

## Subsequent closure

The follow-up implementation includes the full Hartree-plus-XC density
response, differentiates POTCAR `PSPRHO` on the fixed-density Harris path, and
matches VASP interpolation for the local potential, `PSPCOR`, and `PSPRHO`.
The strict oracle is now VASP `ICHARG=11`, which retains the convergence-force
correction. The resulting force-component MAEs are `6.99e-6 eV/Angstrom` for
Si and `1.007e-3 eV/Angstrom` for HfO2. See
[`HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.md`](HARRIS_FORCE_L_RESPONSE_OPTIMIZATION.md)
for the derivation, $L$-channel diagnostic, and raw-data locations.
