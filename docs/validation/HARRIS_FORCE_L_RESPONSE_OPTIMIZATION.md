# Harris-force L-channel and density-response optimization

This note records the PAW augmentation angular-momentum diagnostics and the
completed frozen-density Harris-force response on `force-accuracy-experiments`.
Tests used bare-metal RTX Pro 6000, CUDA 13.2, NVHPC 26.5, PBE, and VASP
`IBZKPT`: 520 eV/24 bands for Si and 521 eV/80 bands for HfO2.

## 1. Consistent L-channel truncation

The augmentation kernel is

$$
Q_{ij}^{I}(\mathbf r)=\sum_{LM}C_{ij}^{LM}Q_{ij}^{L}
g_L(r_I)Y_{LM}(\widehat{\mathbf r_I}),
\qquad \mathbf r_I=\mathbf r-\mathbf R_I.
$$

The same expansion enters the onsite strength, augmentation density, and its
analytic force:

$$
D_{ij}^{I}=D_{ij}^{\mathrm{ion}}+
\int V_{\mathrm{eff}}(\mathbf r)Q_{ij}^{I}(\mathbf r)\,d\mathbf r,
$$

$$
\rho_{\mathrm{aug}}=\sum_{Iij}P_{ij}^{I}Q_{ij}^{I},
\qquad
\mathbf F_I^{\mathrm{aug}}=-\sum_{ij}P_{ij}^{I}
\frac{\partial D_{ij}^{I}}{\partial\mathbf R_I}.
$$

The experimental `HALF_QDEP_LMAX` therefore truncates SETDIJ, augmentation
density, and analytic `dD/dR` together. Its production default is `-1`, meaning
all POTCAR channels. HfO2 produced bitwise-identical energies and forces for
`Lmax=0,2,4`: symmetry removes the nonspherical onsite multipoles after full
k-point weighting. Higher-L quadrature was therefore not the source of the
force discrepancy.

## 2. Missing Harris response

With $\Delta n=n_{\mathrm{out}}-n_{\mathrm{in}}$, the required response field is

$$
g_{H\!xc}=v_H[\Delta n]+f_{xc}[n_{\mathrm{in}}+n_{\mathrm{core}}]\Delta n,
$$

and the ionic correction is

$$
\mathbf F_I^{\mathrm{Harris}}=-\int g_{H\!xc}(\mathbf r)
\frac{\partial n_{\mathrm{atom}}^I(\mathbf r-\mathbf R_I)}
{\partial\mathbf R_I}\,d\mathbf r.
$$

Here $n_{\mathrm{atom}}$ is POTCAR `PSPRHO`; `PSPCOR` belongs only to the
separate NLCC/`FORCOR` term. The old implementation kept only the XC kernel and
contracted it with `PSPCOR`, omitting the dominant Hartree response. The new
implementation evaluates

$$
g_{H\!xc}\approx
\frac{V_{\mathrm{eff}}[n_{\mathrm{in}}+t\Delta n]
-V_{\mathrm{eff}}[n_{\mathrm{in}}-t\Delta n]}{2t},
\qquad t=10^{-4}.
$$

No atom is displaced; this is a density-direction derivative inside the
analytic force, not a finite difference of total energies.

## 3. POTCAR interpolation

CPU and CUDA now share the POTCAR conventions: the local-potential table uses
a cubic spline with zero first derivative at $G=0$ and a natural far boundary,
while `PSPCOR` and `PSPRHO` use VASP's four-point cubic formula. CUDA evaluates
the four-point formula on device and no longer allocates unused core-spline
second derivatives. Unit tests cover cubic exactness and spline boundaries.

## 4. Fixed-density validation

`ICHARG=1, NELM=1, ALGO=All` reaches zero eigensolver residual after one direct
diagonalization, so VASP marks the electronic step converged and omits its
convergence correction. A like-for-like Harris-force oracle must explicitly
freeze the same CHGCAR with `ICHARG=11`.

| System | Old HALF MAE | New HALF MAE | New maximum error | MAE improvement |
|---|---:|---:|---:|---:|
| Si | `2.01117e-4` | `6.99076e-6` | `1.51799e-5` | `28.8x` |
| HfO2 | `3.64625e-2` | `1.00706e-3` | `2.57082e-3` | `36.2x` |

All force values are in eV/Angstrom. For Hf atom 1, VASP's convergence
correction is approximately `(0.192, 0.0448, 0.0307)` eV/Angstrom and HALF
gives `(0.1906, 0.0437, 0.0299)` eV/Angstrom.

| System | HALF free energy (eV) | HALF–VASP Harris-force MAE | Maximum error |
|---|---:|---:|---:|
| Si | `-43.3066227505` | `6.99076e-6` | `1.51799e-5` |
| HfO2 | `-121.8844407617` | `1.00706e-3` | `2.57082e-3` |

CUDA direct-grid, CPU direct-grid, and legacy-shell builds all pass `9/9`
CTest. Final data are under
`/data/limusen/deepaw_half_force_compare/final_vasp_interp_gpu`; the VASP
oracles are under
`/data/limusen/mimic_us_energy_force_compare/{Si,HfO2}_harris_icharg11`.
The same values are available in
[`harris_force_l_response.json`](harris_force_l_response.json).
