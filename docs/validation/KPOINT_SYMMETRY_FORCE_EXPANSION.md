# Irreducible-k space-group expansion for analytic forces

This note records the implementation on branch `force-kpoint-symmetry` that
allows standalone HALF analytic forces to use irreducible k points without
breaking crystal symmetry.

## Method

For every irreducible representative $\mathbf k$, HALF diagonalizes the PAW
Hamiltonian only once. Its multiplicity already includes its full k-star and,
for a nonmagnetic scalar calculation, the time-reversal partner. The scalar
density contribution is expanded with the spatial operations
$g=(R_g,\tau_g)$:

$$
\overline n_{\mathbf k}(\mathbf r)=
\frac{1}{|G|}\sum_{g\in G}
n_{\mathbf k}\!\left(R_g^{-1}(\mathbf r-\tau_g)\right).
$$

HALF first reconstructs the smooth wavefunction density and POTCAR PAW
augmentation density from the representative onsite occupations, then applies
this expansion to their sum. This is equivalent to rotating every wavefunction
and PAW occupation matrix through the star, but does not store full-mesh
wavefunctions or explicit Wigner matrices.

For an atom $I$ mapped to $gI$, a Cartesian force transforms as

$$
\overline{\mathbf F}_{gI}\mathrel{+}=
\frac{1}{|G|}\,C_g\mathbf F_I,
\qquad
C_g=L^{T}R_g(L^{T})^{-1}.
$$

The same atom permutation and Cartesian rotation are applied separately to
the Ewald, local, projector, PAW augmentation, NLCC, and Harris-response force
components. Time reversal needs no additional transformation because charge
density and nonmagnetic atomic force are time-reversal even.

The space-group operations and fractional translations come from spglib.
When an irreducible mesh is selected, the input DeepAW density is symmetrized
with the same operations before $V_{\mathrm{eff}}$ is built. This makes the
Hamiltonian and k-point reduction use one consistent symmetry group.

## Si primitive-cell validation

The CLI test used a two-atom diamond-Si primitive cell on bare-metal RTX PRO
6000, PBE, 520 eV, `KSPACING=0.35`, 12 bands, CUDA EVD, and one GPU.

| quantity | full mesh | irreducible + expansion |
|---|---:|---:|
| k points solved | 216 | 16 |
| spatial operations | 1 | 48 |
| free energy (eV/cell) | `-10.8288461353` | `-10.8288461451` |
| energy difference | reference | `-9.78e-9 eV` |
| maximum force norm | `1.65e-7 eV/Angstrom` | `2.45e-18 eV/Angstrom` |
| maximum force-component difference | reference | `1.06e-7 eV/Angstrom` |
| wall time | `35.41 s` | `4.42 s` |

The force difference is at the numerical-noise level of the explicit full
mesh. The irreducible implementation is **8.01x faster** and uses 13.5 times
fewer diagonalizations.

The JSON output records `space_group_operation_count` and sets
`force_kpoint_expansion` to `observable-space-group` when this path is active.

## Regression

CUDA direct-grid, CPU direct-grid, and legacy-shell builds each pass `9/9`
CTest. Machine-readable results are in
[`kpoint_symmetry_force_expansion.json`](kpoint_symmetry_force_expansion.json).
Raw runs are under
`/data/limusen/deepaw_half_si_primitive_cli_timing_20260921`.
