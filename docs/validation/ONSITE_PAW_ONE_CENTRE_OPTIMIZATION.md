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

The position derivative acts on the same compact grid kernel used for the
energy. `half-force-check` uses a displaced-ion central difference only as an
independent development oracle. The maximum errors in `dD/dR` are
`2.75e-9 eV/Angstrom` for Si and `2.43e-9 eV/Angstrom` for HfO2, four to five
orders of magnitude below the previous shell-quadrature errors.

| system | HALF free energy (eV) | HALF - VASP (eV/cell) | force-component MAE (eV/Angstrom) |
|---|---:|---:|---:|
| Si | -43.306625404 | -0.013307174 | 0.000739611 |
| HfO2 | -121.884427193 | +0.088155727 | 0.068379790 |

The Si energy difference decreases from about `20.7` to `13.3 meV/cell`.
HfO2 `dD/dR` is now numerically converged but its VASP force difference is
essentially unchanged. This isolates the remaining HfO2 discrepancy away from
QDEP angular integration and radial moments; the next comparison must target
the atomic AE/PS double-counting and the output-augmentation contribution to
the local force.

The analytic-force pass also caches the k-independent `D` and `dD/dR` arrays.
In the contended Pro 6000 run, HfO2 wall time decreased from `121.87` to
`97.78 s`; this is evidence that repeated construction was removed, not an
exclusive-GPU performance claim.

Machine-readable results are in
[`onsite_paw_direct_grid.json`](onsite_paw_direct_grid.json).
