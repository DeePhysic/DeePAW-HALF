from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HAPPY = ROOT.parent / "happy"


def test_expected_project_surface_exists():
    expected = [
        "CMakeLists.txt",
        "CMakePresets.json",
        "src/half_chgcar.F90",
        "src/half_basis.F90",
        "src/half_cuda.cuf",
        "src/half_cpu.F90",
        "app/half_inspect.F90",
        "app/half_inspect_cpu.F90",
        "app/half_benchmark.F90",
        "docs/PORTING_MATRIX.md",
        "scripts/remote_build.sh",
        "scripts/remote_benchmark.sh",
        "scripts/remote_cuda13_build.sh",
        "tests/compare_happy_basis.py",
    ]
    assert all((ROOT / item).is_file() for item in expected)


def test_build_rejects_non_nvhpc_compilers():
    source = (ROOT / "CMakeLists.txt").read_text()
    assert 'MATCHES "NVHPC|PGI"' in source
    assert "HALF CUDA backend requires NVIDIA HPC SDK" in source
    assert "HALF_ENABLE_CUDA" in source


def test_cuda_fortran_kernel_is_real_cuda_code():
    source = (ROOT / "src/half_cuda.cuf").read_text().lower()
    assert "use cudafor" in source
    assert "attributes(global)" in source
    assert "<<<blocks,threads>>>" in source
    assert "cudadevicesynchronize" in source


def test_physical_constants_follow_happy_definition():
    source = (ROOT / "src/half_constants.F90").read_text().lower()
    assert "autoa = 0.529177249_dp" in source
    assert "rytoev = 13.605826_dp" in source
    assert "hsqdtm = rytoev*autoa*autoa" in source


def test_si_reference_input_contract():
    # Independent lightweight read: confirms the bootstrap target and expected
    # oracle values without importing HAPPY from the source tree.
    path = HAPPY / "examples/si/CHGCAR.smooth"
    lines = path.read_text().splitlines()
    counts = [int(v) for v in lines[6].split()]
    coordinate_line = 7
    if lines[coordinate_line].lower().startswith("s"):
        coordinate_line += 1
    cursor = coordinate_line + 1 + sum(counts)
    while not lines[cursor].strip():
        cursor += 1
    shape = tuple(int(v) for v in lines[cursor].split())
    values = np.fromstring(" ".join(lines[cursor + 1 :]), sep=" ", count=np.prod(shape))
    assert shape == (56, 56, 56)
    assert values.size == np.prod(shape)
    assert abs(values.mean() - 8.000000050054025) < 5e-12

    scale = float(lines[1])
    lattice = scale * np.array([[float(v) for v in lines[i].split()] for i in range(2, 5)])
    reciprocal = 2.0 * np.pi * np.linalg.inv(lattice).T
    hsqdtm = 13.605826 * 0.529177249**2
    npw = 0
    for n1 in range(-shape[0] // 2, (shape[0] - 1) // 2 + 1):
        for n2 in range(-shape[1] // 2, (shape[1] - 1) // 2 + 1):
            for n3 in range(-shape[2] // 2, (shape[2] - 1) // 2 + 1):
                g = np.array([n1, n2, n3]) @ reciprocal
                npw += hsqdtm * np.dot(g, g) < 400.0
    np.testing.assert_allclose(abs(np.linalg.det(lattice)), 40.32952684742029, rtol=1e-12)
    assert npw == 725
