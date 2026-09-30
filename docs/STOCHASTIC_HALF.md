# Stochastic HALF on `large-scalebranch`

## Numerical contract

The input is the fixed DeePAW/VASP smooth density and its PAW potential. HALF
constructs a generalized Hermitian problem `H psi = epsilon S psi`, with
positive-definite `S`. For one spin channel, define `A=S^-1 H`,
`f(E)=1/(1+exp((E-mu)/kBT))`, and `P=f(A) S^-1`. The implementation estimates:

- electron number `Tr f(A)`;
- frozen-Hamiltonian band energy `Tr[A f(A)]`;
- smooth real-space density `diag(F P F*)`;
- onsite PAW matrix `B* P B`.

Random complex Z4 vectors have covariance `E[zz*]=I`. A Chebyshev recurrence
applies the finite-temperature Fermi operator without diagonalizing `H`.
`S^-1` is applied by conjugate gradients, with a caller-supplied exact inverse
allowed for small-system validation. The polynomial interval must contain the
complete generalized spectrum. A generalized Lanczos estimate with explicit
padding is supplied; the maximum finite-temperature polynomial error is
reported separately from Monte Carlo standard error.

`fixed_spin_channels` holds `N_up+N_down` and `N_up-N_down=NUPDOWN` fixed,
using common probes in both channels. `paired_band_energy_difference` also
uses common probes in two structures and reports the paired and independent
standard errors. No random sample is normalized individually. Statistical
density values can be locally negative and must not be clipped.

## Implementation and verified scope

`scripts/stochastic_half.py` is a reusable external solver. The CPU
`half_apply_hs` C ABI applies local potential through FFTs and nonlocal PAW
terms through projector contractions on CPU or CUDA. It no longer allocates
dense `N_PW x N_PW` matrices per application. CPU caches the PAW projectors;
CUDA's stochastic interface builds them in tiles of eight atoms, so its
projector scratch scales as `O(N_PW x 8 x channels_per_atom)` rather than
`O(N_PW x N_projectors)`. The ABI signature is unchanged. `scripts/half_operator_api.py`
and `scripts/run_stochastic_half.py` connect this operator to the stochastic
solver. The current adapter still copies each random vector to/from GPU;
resident GPU batches are required for production throughput.

For a small-cell smoke run (with a library built from this branch):

```bash
python scripts/run_stochastic_half.py \
  --library /path/to/libhalf.so --charge /path/to/CHGCAR \
  --potcar /path/to/POTCAR --encut 50 --backend cuda \
  --electrons-per-spin 4 --kbt 1.0 --degree 128 --samples 64
```

The script estimates the spectral interval by generalized Lanczos when bounds
are omitted. A production calculation must verify adequate padding and report
the Chebyshev polynomial error and stochastic standard error together.
On the 20,341-dimensional Si load test, automatic extremal generalized
Lanczos was interrupted after more than six minutes. Explicit broad bounds
were used for the load test; they have not been certified as production
bounds. A GPU-resident spectral estimator is still needed at this scale.

The real Si primitive-cell gate used VASP CHGCAR, Si POTCAR, a 60^3 charge
grid, `ENCUT=50 eV`, `N_PW=27`, `kBT=0.5 eV`, and a degree-256 polynomial.
Matrix-free versus the dense HALF operator differed by at most
`5.71e-13 eV` for `H psi` and `3.56e-15` for `S psi`. Exhaustive basis probes
against the exact generalized eigensolution differed by `2.12e-7` electrons
per spin, `1.17e-6 eV` in band energy, and `5.63e-8` relative L1 in the
smooth density. These measure polynomial/operator bias, not stochastic noise.
CUDA matrix-free `H psi` and `S psi` differed from the dense CPU oracle by
`4.86e-12 eV` and `3.58e-15`, respectively. A direct GPU stochastic run
completed with degree 96 and four probes.
For an 18-atom `3x3x1` Si supercell (293 plane waves, 180x180x60 grid), the
eight-atom tile boundary and final two-atom tile matched the dense HALF
operator within `5.26e-12 eV` for `H psi` and `8.59e-15` for `S psi`.
`scripts/check_matrix_free_hs.py` reproduces this operator check for any
small case with a CHGCAR and POTCAR.
At 64 Z4 probes, the random electron estimate was `3.828 +/- 0.154` versus
exact `4.000`, and band energy was `13.872 +/- 0.953 eV` versus exact
`14.768 eV`. A direct four-probe run through the new CPU matrix-free C API
also completed, with expected large statistical error. See
`docs/validation/stochastic_half_si_20261001.json`.

For the same Si geometry with VASP and DeePAW charge inputs, 16 paired CUDA
probes gave a band-energy difference of `-28.995 +/- 1.800 meV`; the exact
generalized eigensolution gave `-29.535 meV`. Treating those two stochastic
traces as independent would give a `1.479 eV` standard error. This is a
fixed-Hamiltonian band-energy test, not an NRR reaction energy.

An 8000-atom, 60 A load test used a uniformly strained Si supercell and one
coarse 110^3 eSCN density grid. At `ENCUT=12 eV`, the basis had 20,341 plane
waves and could hold the 16,000 electrons per spin. The tiled CUDA operator
completed its first H/S application in `7.23 s` and a cached application in
`0.45 s`; both outputs were finite. A complete one-probe Fermi trace also
finished. The broad test bounds and degree-32 polynomial gave a maximum
occupation approximation error of `0.004`, and one probe has no statistical
error estimate. These values establish execution and memory feasibility only;
see `docs/validation/stochastic_half_si8000_scale_only_20261001.json`.
With two probes and degree 64, the polynomial occupation error fell to
`2.61e-5`, but the per-spin band trace had a `166 eV` standard error.
This confirms that direct absolute energy estimates need many probes or
control variates even when the operator fits in GPU memory.

## Production gates still required

The present path does not yet constitute a large-system PAW total-energy or
force engine. The CPU operator retains `O(N_PW N_projectors)` projectors, and
CPU FFT plans are recreated per call. The CUDA tile path removes the full
projector allocation but recomputes projectors for every H/S application.
A 10,648-atom Si run still requires batched resident GPU random vectors,
operator timing and peak-memory measurements at its physical cutoff.
For NRR energetics, the onsite PAW double counting, augmentation, Hartree/XC
and ionic terms must be connected to the stochastic occupancy and density,
then forces and independent-seed uncertainty validated. A band trace alone
cannot be reported as a reaction energy or barrier.

Small-system convergence gates should independently vary Chebyshev degree,
temperature, spectral padding, overlap CG tolerance, probe count and seed.
For reaction differences, paired probes and deterministic/local control
variates are needed until the confidence interval is well below the chemical
energy scale of interest.
