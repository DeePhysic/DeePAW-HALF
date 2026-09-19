from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture(
    params=["force_and_stress", "local_potential", "structure", "occupancies"]
)
def function_name(request):
    return request.param


@pytest.fixture
def mock_entry_points():
    return_value = (MagicMock(), MagicMock())
    with patch("vasp._entry_points.entry_points", return_value=return_value) as mock:
        yield mock
