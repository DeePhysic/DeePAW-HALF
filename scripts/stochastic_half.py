#!/usr/bin/env python3
"""Matrix-free stochastic Fermi operator for a fixed generalized PAW H/S.

This is an external numerical prototype.  The caller supplies applications of
H and S; no eigenvectors or dense Hamiltonian are formed by this module.
Energy is the frozen-Hamiltonian band trace, not a PAW total energy.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.fft import dct
from scipy.special import expit
from scipy.sparse.linalg import LinearOperator, cg


Vector = np.ndarray
Apply = Callable[[Vector], Vector]


@dataclass(frozen=True)
class GeneralizedOperator:
    size: int
    apply_h: Apply
    apply_s: Apply
    lower_ev: float
    upper_ev: float
    overlap_tolerance: float = 1e-11
    overlap_maxiter: int = 200
    inverse_overlap: Apply | None = None

    def __post_init__(self) -> None:
        if self.size < 1 or not self.lower_ev < self.upper_ev:
            raise ValueError("positive basis size and ordered spectral bounds required")

    def solve_s(self, rhs: Vector) -> Vector:
        norm_s = np.vdot(rhs, self.apply_s(rhs))
        if norm_s.real <= 0 or abs(norm_s.imag) > 1e-9 * abs(norm_s.real):
            raise ValueError("PAW overlap must be Hermitian positive definite")
        if self.inverse_overlap is not None:
            return np.asarray(self.inverse_overlap(rhs), dtype=np.complex128)
        metric = LinearOperator(
            (self.size, self.size), matvec=self.apply_s, dtype=np.complex128
        )
        result, info = cg(
            metric, rhs, rtol=self.overlap_tolerance, atol=0.0,
            maxiter=self.overlap_maxiter,
        )
        if info != 0:
            raise RuntimeError(f"PAW overlap CG failed (info={info})")
        return result

    def apply_a(self, vector: Vector) -> Vector:
        # S^-1 H is self-adjoint in the S metric and similar to S^-1/2 H S^-1/2.
        return self.solve_s(self.apply_h(vector))

    def apply_scaled(self, vector: Vector) -> Vector:
        width = self.upper_ev - self.lower_ev
        return (2.0 * self.apply_a(vector) -
                (self.upper_ev + self.lower_ev) * vector) / width


def estimate_spectral_bounds(
    size: int, apply_h: Apply, apply_s: Apply, inverse_overlap: Apply,
    *, padding_ev: float = 1.0, tolerance: float = 1e-6,
) -> tuple[float, float]:
    """Generalized Lanczos estimates with a conservative explicit padding."""
    from scipy.sparse.linalg import eigsh
    if size < 3 or padding_ev <= 0:
        raise ValueError("basis size >= 3 and positive spectral padding required")
    shape = (size, size)
    h = LinearOperator(shape, matvec=apply_h, dtype=np.complex128)
    s = LinearOperator(shape, matvec=apply_s, dtype=np.complex128)
    si = LinearOperator(shape, matvec=inverse_overlap, dtype=np.complex128)
    minimum = float(eigsh(h, k=1, M=s, Minv=si, which="SA", tol=tolerance,
                          return_eigenvectors=False)[0])
    maximum = float(eigsh(h, k=1, M=s, Minv=si, which="LA", tol=tolerance,
                          return_eigenvectors=False)[0])
    return minimum - padding_ev, maximum + padding_ev


def fermi_coefficients(
    lower_ev: float, upper_ev: float, mu_ev: float, temperature_ev: float,
    degree: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Chebyshev coefficients of f and epsilon*f (finite T only)."""
    if temperature_ev <= 0 or degree < 1:
        raise ValueError("finite positive kBT and polynomial degree required")
    nodes = max(4 * (degree + 1), 128)
    theta = np.pi * (np.arange(nodes) + 0.5) / nodes
    energies = 0.5 * (upper_ev + lower_ev) + 0.5 * (upper_ev - lower_ev) * np.cos(theta)
    occupations = expit((mu_ev - energies) / temperature_ev)
    def coefficients(values: np.ndarray) -> np.ndarray:
        result = dct(values, type=2)[:degree + 1] / nodes
        result[0] *= 0.5
        return result
    return coefficients(occupations), coefficients(energies * occupations)


def polynomial_errors(
    lower_ev: float, upper_ev: float, mu_ev: float, temperature_ev: float,
    degree: int,
) -> dict[str, float]:
    """Resolve finite-T truncation separately from Monte Carlo noise."""
    from numpy.polynomial.chebyshev import chebval
    coefficients_n, coefficients_e = fermi_coefficients(
        lower_ev, upper_ev, mu_ev, temperature_ev, degree,
    )
    x = np.linspace(-1, 1, max(4 * degree, 1001))
    energies = 0.5 * (upper_ev + lower_ev) + 0.5 * (upper_ev-lower_ev) * x
    exact = expit((mu_ev - energies) / temperature_ev)
    return {
        "occupation_max_abs": float(np.max(np.abs(chebval(x, coefficients_n)-exact))),
        "energy_weighted_max_abs_ev": float(np.max(np.abs(
            chebval(x, coefficients_e)-energies*exact))),
    }


def probe_vectors(size: int, samples: int, seed: int) -> np.ndarray:
    """Independent complex Z4 vectors satisfying E[zz*]=I."""
    if samples < 1:
        raise ValueError("at least one stochastic probe is required")
    rng = np.random.default_rng(seed)
    phases = np.array([1, 1j, -1, -1j], dtype=np.complex128)
    return phases[rng.integers(0, 4, (samples, size))]


class StochasticHalf:
    def __init__(
        self, operator: GeneralizedOperator, degree: int, samples: int,
        seed: int = 0, probes: np.ndarray | None = None,
    ) -> None:
        if degree < 1:
            raise ValueError("degree must be positive")
        self.operator = operator
        self.degree = degree
        self.probes = (probe_vectors(operator.size, samples, seed) if probes is None
                       else np.asarray(probes, dtype=np.complex128))
        if self.probes.ndim != 2 or self.probes.shape[1] != operator.size:
            raise ValueError("probe array must have shape (samples, basis size)")
        self.moments = np.empty((len(self.probes), degree + 1), dtype=np.float64)
        for i, probe in enumerate(self.probes):
            for k, vector in enumerate(self._vectors(probe)):
                self.moments[i, k] = np.vdot(probe, vector).real

    def _vectors(self, vector: Vector):
        previous = np.asarray(vector, dtype=np.complex128)
        yield previous
        current = self.operator.apply_scaled(previous)
        yield current
        for _ in range(2, self.degree + 1):
            following = 2.0 * self.operator.apply_scaled(current) - previous
            yield following
            previous, current = current, following

    def _coefficient_pair(self, mu_ev: float, temperature_ev: float):
        return fermi_coefficients(
            self.operator.lower_ev, self.operator.upper_ev,
            mu_ev, temperature_ev, self.degree,
        )

    def trace(self, mu_ev: float, temperature_ev: float) -> dict[str, float | None]:
        counts, energies = self.trace_samples(mu_ev, temperature_ev)
        def stderr(values: np.ndarray) -> float | None:
            return float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else None
        return {
            "electrons": float(counts.mean()),
            "electrons_stderr": stderr(counts),
            "band_energy_ev": float(energies.mean()),
            "band_energy_stderr_ev": stderr(energies),
        }

    def trace_samples(self, mu_ev: float, temperature_ev: float) -> tuple[np.ndarray, np.ndarray]:
        occupation, energy = self._coefficient_pair(mu_ev, temperature_ev)
        counts = self.moments @ occupation
        energies = self.moments @ energy
        return counts, energies

    def find_mu(
        self, target_electrons: float, temperature_ev: float,
        tolerance: float = 1e-7,
    ) -> float:
        from scipy.optimize import brentq
        if not 0 < target_electrons < self.operator.size:
            raise ValueError("one spin channel needs 0 < electrons < basis size")
        width = self.operator.upper_ev - self.operator.lower_ev
        lo = self.operator.lower_ev - width - 20 * temperature_ev
        hi = self.operator.upper_ev + width + 20 * temperature_ev
        root = lambda mu: self.trace(mu, temperature_ev)["electrons"] - target_electrons
        if root(lo) >= 0 or root(hi) <= 0:
            raise RuntimeError("sampled electron count does not bracket the target; increase probes/order")
        return float(brentq(root, lo, hi, xtol=tolerance))

    def observables(
        self, mu_ev: float, temperature_ev: float,
        *,
        real_space_map: Apply | None = None,
        paw_projectors: np.ndarray | None = None,
    ) -> dict[str, np.ndarray | dict[str, float]]:
        """Estimate P=f(S^-1 H) S^-1, real-space diagonal, and B* P B.

        The map sends plane-wave coefficients to real-space grid values.  The
        projector matrix has basis vectors as rows and PAW channels as columns.
        Outputs are raw spin-channel occupations; no per-sample normalization.
        """
        coefficients, _ = self._coefficient_pair(mu_ev, temperature_ev)
        density_sum = None
        onsite_sum = None
        projectors = None if paw_projectors is None else np.asarray(paw_projectors)
        if projectors is not None and projectors.shape[0] != self.operator.size:
            raise ValueError("PAW projector first dimension must equal basis size")
        for probe in self.probes:
            # Generalized density kernel: P=f(S^-1 H) S^-1.
            source = self.operator.solve_s(probe)
            result = np.zeros_like(source)
            for coefficient, vector in zip(coefficients, self._vectors(source)):
                result += coefficient * vector
            if real_space_map is not None:
                density_sample = np.real(
                    real_space_map(result) * np.conj(real_space_map(probe))
                )
                if density_sum is None:
                    density_sum = np.zeros_like(density_sample)
                density_sum += density_sample
            if projectors is not None:
                left = projectors.conj().T @ result
                right = projectors.conj().T @ probe
                if onsite_sum is None:
                    onsite_sum = np.zeros((len(left), len(left)), dtype=np.complex128)
                onsite_sum += np.outer(left, right.conj())
        output: dict[str, np.ndarray | dict[str, float]] = {"trace": self.trace(mu_ev, temperature_ev)}
        if density_sum is not None:
            output["smooth_density"] = density_sum / len(self.probes)
        if onsite_sum is not None:
            matrix = onsite_sum / len(self.probes)
            output["paw_onsite_occupancy"] = 0.5 * (matrix + matrix.conj().T)
        return output


def fixed_spin_channels(
    up: GeneralizedOperator, down: GeneralizedOperator,
    total_electrons: float, nupdown: float, temperature_ev: float,
    degree: int, samples: int, seed: int = 0,
) -> dict[str, object]:
    """Constrain N_up-N_down with shared Z4 probes in both spin channels."""
    if up.size != down.size:
        raise ValueError("spin channels must use the same plane-wave basis")
    targets = ((total_electrons + nupdown) / 2, (total_electrons - nupdown) / 2)
    shared_probes = probe_vectors(up.size, samples, seed)
    channels = []
    for operator, target in zip((up, down), targets):
        solver = StochasticHalf(operator, degree, samples, probes=shared_probes)
        mu = solver.find_mu(target, temperature_ev)
        channels.append({"target_electrons": target, "mu_ev": mu,
                         "trace": solver.trace(mu, temperature_ev), "solver": solver})
    return {"nupdown": nupdown, "up": channels[0], "down": channels[1]}


def paired_band_energy_difference(
    first: GeneralizedOperator, second: GeneralizedOperator,
    electrons_per_spin: float, temperature_ev: float,
    degree: int, samples: int, seed: int = 0,
) -> dict[str, float | None]:
    """Correlated structural energy difference using identical random probes."""
    if first.size != second.size:
        raise ValueError("paired states require equal plane-wave basis size")
    probes = probe_vectors(first.size, samples, seed)
    a = StochasticHalf(first, degree, samples, probes=probes)
    b = StochasticHalf(second, degree, samples, probes=probes)
    mu_a = a.find_mu(electrons_per_spin, temperature_ev)
    mu_b = b.find_mu(electrons_per_spin, temperature_ev)
    _, energy_a = a.trace_samples(mu_a, temperature_ev)
    _, energy_b = b.trace_samples(mu_b, temperature_ev)
    paired = energy_b - energy_a
    return {
        "first_mu_ev": mu_a,
        "second_mu_ev": mu_b,
        "band_energy_difference_ev": float(paired.mean()),
        "paired_stderr_ev": float(paired.std(ddof=1) / np.sqrt(samples)) if samples > 1 else None,
        "unpaired_stderr_ev": float(np.sqrt((energy_a.var(ddof=1) + energy_b.var(ddof=1)) / samples))
        if samples > 1 else None,
    }


def exact_reference(h: np.ndarray, s: np.ndarray, mu_ev: float, temperature_ev: float):
    """Small-system generalized diagonalization oracle for external validation."""
    from scipy.linalg import eigh
    eigenvalues, eigenvectors = eigh(h, s)
    occupation = expit((mu_ev - eigenvalues) / temperature_ev)
    density_kernel = (eigenvectors * occupation) @ eigenvectors.conj().T
    return {
        "electrons": float(occupation.sum()),
        "band_energy_ev": float(occupation @ eigenvalues),
        "density_kernel": density_kernel,
        "eigenvalue_bounds": (float(eigenvalues[0]), float(eigenvalues[-1])),
    }
