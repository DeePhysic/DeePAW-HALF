"""Numerical checks against a generalized PAW-like eigenproblem."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from stochastic_half import (GeneralizedOperator, StochasticHalf, exact_reference,
                             estimate_spectral_bounds, fixed_spin_channels,
                             paired_band_energy_difference)


def test_generalized_fermi_operator_electron_energy_density_and_onsite():
    rng = np.random.default_rng(17)
    n = 9
    trial = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = 0.15 * (trial + trial.conj().T) + np.diag(np.linspace(-2, 3, n))
    projectors = (rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))) / 5
    q = np.diag([0.4, -0.15, 0.2])
    s = np.eye(n) + projectors @ q @ projectors.conj().T
    mu, kbt = 0.15, 0.65
    oracle = exact_reference(h, s, mu, kbt)
    lower, upper = oracle["eigenvalue_bounds"]
    operator = GeneralizedOperator(
        n, lambda v: h @ v, lambda v: s @ v,
        lower - 0.3, upper + 0.3,
    )
    probes = np.sqrt(n) * np.eye(n, dtype=complex)
    solver = StochasticHalf(operator, degree=96, samples=n, probes=probes)
    observed = solver.observables(
        mu, kbt, real_space_map=lambda v: np.fft.fft(v, norm="ortho"),
        paw_projectors=projectors,
    )
    f = np.fft.fft(np.eye(n), axis=0, norm="ortho")
    density_ref = np.diag(f @ oracle["density_kernel"] @ f.conj().T).real
    onsite_ref = projectors.conj().T @ oracle["density_kernel"] @ projectors
    assert abs(observed["trace"]["electrons"] - oracle["electrons"]) < 1e-9
    assert abs(observed["trace"]["band_energy_ev"] - oracle["band_energy_ev"]) < 1e-8
    np.testing.assert_allclose(observed["smooth_density"], density_ref, atol=1e-9)
    np.testing.assert_allclose(observed["paw_onsite_occupancy"], onsite_ref, atol=1e-9)
    assert abs(solver.find_mu(4.25, kbt) - mu) > 1e-3  # independent electron target
    assert abs(solver.trace(solver.find_mu(4.25, kbt), kbt)["electrons"] - 4.25) < 1e-7


def test_complex_z4_probe_and_overlap_failure():
    n = 6
    h = np.diag(np.linspace(-1, 1, n))
    s = np.eye(n)
    operator = GeneralizedOperator(n, lambda v: h @ v, lambda v: s @ v, -1.1, 1.1)
    solver = StochasticHalf(operator, degree=50, samples=12, seed=9)
    assert np.all(np.isin(solver.probes, [1, -1, 1j, -1j]))
    assert np.isfinite(solver.trace(0, 0.4)["electrons_stderr"])
    broken = GeneralizedOperator(n, lambda v: h @ v, lambda v: -s @ v, -1.1, 1.1)
    try:
        broken.solve_s(np.ones(n, dtype=complex))
    except ValueError:
        pass
    else:
        raise AssertionError("indefinite PAW overlap must fail")


def test_fixed_nupdown_uses_shared_probes_and_channel_targets():
    n = 8
    h = np.diag(np.linspace(-2, 3, n))
    s = np.eye(n)
    up = GeneralizedOperator(n, lambda v: h @ v, lambda v: s @ v, -2.2, 3.2)
    result = fixed_spin_channels(up, up, total_electrons=7, nupdown=1,
                                 temperature_ev=0.4, degree=90, samples=16)
    assert abs(result["up"]["trace"]["electrons"] - 4) < 1e-7
    assert abs(result["down"]["trace"]["electrons"] - 3) < 1e-7
    np.testing.assert_array_equal(result["up"]["solver"].probes,
                                  result["down"]["solver"].probes)


def test_generalized_lanczos_bounds_and_correlated_energy_difference():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(8, 8))
    h = np.diag(np.linspace(-2, 3, 8)) + 0.15 * (x + x.T)
    s = np.eye(8) + 0.05 * x.T @ x
    bounds = estimate_spectral_bounds(8, lambda v: h @ v, lambda v: s @ v,
                                      lambda v: np.linalg.solve(s, v), padding_ev=0.2)
    exact = np.linalg.eigvals(np.linalg.solve(s, h)).real
    assert bounds[0] < min(exact) and bounds[1] > max(exact)
    first = GeneralizedOperator(8, lambda v: h @ v, lambda v: s @ v, *bounds,
                                inverse_overlap=lambda v: np.linalg.solve(s, v))
    shifted_h = h + 0.03 * s
    second = GeneralizedOperator(8, lambda v: shifted_h @ v, lambda v: s @ v,
                                 bounds[0], bounds[1]+0.1,
                                 inverse_overlap=lambda v: np.linalg.solve(s, v))
    result = paired_band_energy_difference(first, second, 4, 0.4, 100, 64, 3)
    assert abs(result["band_energy_difference_ev"] - 0.12) < 0.03
    assert result["paired_stderr_ev"] < result["unpaired_stderr_ev"] * 0.1
