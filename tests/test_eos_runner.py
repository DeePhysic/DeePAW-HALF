"""Small, data-free checks for the resumable EOS benchmark protocol."""

import json
from pathlib import Path

import pytest

from scripts.eos_half_360 import onsite_lmax_for_potcar, valid_energy


@pytest.mark.parametrize("projector_l,expected", [(1, 2), (2, 4), (3, 6)])
def test_auto_onsite_channel_from_potcar(tmp_path: Path, projector_l: int, expected: int):
    potcar = tmp_path / "POTCAR"
    potcar.write_text(f"non local part\n  0  1.0\nnon local part\n  {projector_l}  1.0\n")
    assert onsite_lmax_for_potcar(potcar, "auto") == expected
    assert onsite_lmax_for_potcar(potcar, None) is None
    assert onsite_lmax_for_potcar(potcar, "0") == 0


def test_auto_onsite_channel_rejects_missing_projectors(tmp_path: Path):
    potcar = tmp_path / "POTCAR"
    potcar.write_text("no projector information\n")
    with pytest.raises(ValueError, match="Could not read projector"):
        onsite_lmax_for_potcar(potcar, "auto")


def test_valid_energy_requires_finite_analytic_forces(tmp_path: Path):
    result = tmp_path / "harris.json"
    result.write_text(json.dumps({"internal_energy_eV": 1.0,
                                  "forces_eV_per_Angstrom": [[0.0, 0.1, 0.2]]}))
    assert valid_energy(result, 1)
    assert not valid_energy(result, 2)
    result.write_text('{"internal_energy_eV": 1.0, "forces_eV_per_Angstrom": [[NaN, 0, 0]]}')
    assert not valid_energy(result, 1)
