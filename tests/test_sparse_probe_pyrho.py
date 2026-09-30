"""Contract tests for Fourier reconstruction of shifted sparse probes."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

pyrho = pytest.importorskip("pyrho")
from pyrho.utils import interpolate_fourier  # noqa: E402


def _aligned_fft(sparse: np.ndarray, residue: tuple[int, ...], shape: tuple[int, ...]) -> np.ndarray:
    interpolated = interpolate_fourier(sparse, list(shape))
    return np.roll(interpolated, shift=residue, axis=tuple(range(len(shape))))


def test_shifted_pyrho_fft_recovers_band_limited_periodic_field() -> None:
    shape = (12, 12, 12)
    stride = (2, 2, 2)
    x, y, z = np.meshgrid(
        *(np.arange(n) / n for n in shape), indexing="ij"
    )
    dense = (
        0.3
        + 0.2 * np.sin(2 * np.pi * (2 * x + y))
        - 0.1 * np.cos(2 * np.pi * (y - 2 * z))
    )

    for residue in itertools.product(*(range(s) for s in stride)):
        target = tuple(slice(r, None, s) for r, s in zip(residue, stride))
        reconstructed = _aligned_fft(dense[target], residue, shape)
        assert np.max(np.abs(reconstructed - dense)) < 1.0e-12


def test_direct_interlacing_preserves_high_frequency_that_fft_loses() -> None:
    shape = (12, 12, 12)
    stride = (2, 2, 2)
    x, y, z = np.meshgrid(
        *(np.arange(n) / n for n in shape), indexing="ij"
    )
    dense = (
        0.3
        + 0.2 * np.sin(2 * np.pi * (4 * x + y))
        - 0.1 * np.cos(2 * np.pi * (y - 4 * z))
    )
    direct = np.empty(shape)
    fft_candidates = []

    for residue in itertools.product(*(range(s) for s in stride)):
        target = tuple(slice(r, None, s) for r, s in zip(residue, stride))
        sparse = dense[target]
        direct[target] = sparse
        fft_candidates.append(_aligned_fft(sparse, residue, shape))

    fft_consensus = np.mean(fft_candidates, axis=0)
    assert np.max(np.abs(direct - dense)) == 0.0
    assert np.sqrt(np.mean((fft_consensus - dense) ** 2)) > 0.1
