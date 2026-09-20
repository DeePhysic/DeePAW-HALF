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
        "src/half_kpoints.F90",
        "src/half_energy.F90",
        "src/half_parallel.F90",
        "src/half_library.F90",
        "src/half_c_api.F90",
        "src/half_escn_client.c",
        "src/half_escn_socket.c",
        "examples/api/half_escn_example.c",
        "src/half_artifacts.c",
        "include/half.h",
        "include/half_api.f90",
        "docs/API.md",
        "docs/API.zh-CN.md",
        "src/half_vaspwave.F90",
        "src/half_hdf5_bridge.c",
        "src/half_cuda.cuf",
        "src/half_cuda_potential.cuf",
        "src/half_cuda_uspp.cuf",
        "src/half_uspp.F90",
        "src/half_cuda_assembly.cuf",
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


def test_stable_library_api_and_cli_artifacts_are_exposed():
    header = (ROOT / "include/half.h").read_text()
    cmake = (ROOT / "CMakeLists.txt").read_text()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    library = (ROOT / "src/half_library.F90").read_text().lower()
    cuda_solver = (ROOT / "src/half_cuda_solver.cuf").read_text().lower()
    assert "HALF_ABI_VERSION 1" in header
    assert "half_create_from_files_backend" in header
    assert "half_assemble_hs" in header
    assert "half_apply_hs" in header
    assert "half_solve_kpoint" in header
    assert "half_solve_kpoint_mapped" in header
    assert "half_request_geometry_v1" in header
    assert "half_set_request_geometry" in header
    assert "half_get_request_geometry" in header
    assert "half_escn_predict" in header
    assert "half_create_from_escn" in header
    assert "HALF_CAP_ESCN_API" in header
    assert "EXPORT_NAME half" in cmake
    assert "HALFConfig.cmake" in cmake
    assert "write_bands_artifacts" in cli
    assert "write_energy_npz" in cli
    assert "solve_dense_gamma_cuda_full" in library
    assert "assemble_dense_gamma_cuda_full" in cuda_solver


def test_private_vasp_one_click_build_is_wired():
    cmake = (ROOT / "CMakeLists.txt").read_text()
    integration = (ROOT / "cmake" / "HALFVaspIntegration.cmake").read_text()
    tool = (ROOT / "tools" / "half-cmake").read_text()
    adapter = ROOT / "vendor" / "vasp-6.6.0" / "src" / "half_vasp_init.F"
    assert 'option(HALF_BUILD_VASP' in cmake
    assert "add_custom_target(vasp-half ALL" in integration
    assert "PrepareMklFftw.cmake" in integration
    assert 'vasp) preset="base"' in tool
    assert adapter.is_file()
    adapter_source = adapter.read_text()
    vasp_main = (ROOT / "vendor" / "vasp-6.6.0" / "src" / "main.F").read_text()
    vasp_reader = (ROOT / "vendor" / "vasp-6.6.0" / "src" / "reader.F").read_text()
    vasp_base = (ROOT / "vendor" / "vasp-6.6.0" / "src" / "base.F").read_text()
    assert "half_set_request_geometry" in adapter_source
    assert "T_INFO%POSION" in adapter_source
    assert "HALF_ESCN_URL" in adapter_source
    assert "USE_ESCN_API" in adapter_source
    assert "half_create_from_escn" in adapter_source
    assert "INFO%LHALF_API" in vasp_main
    assert "INFO%HALF_ESCN_URL" in vasp_main
    assert "'LHALF_API'" in vasp_reader
    assert "'HALF_ESCN_URL'" in vasp_reader
    assert "LOGICAL LHALF_API" in vasp_base
    assert "CHARACTER(255) HALF_ESCN_URL" in vasp_base
    assert "GRIDC%NGPTAR" in vasp_main


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


def test_full_cuda_gamma_pipeline_is_device_resident():
    potential = (ROOT / "src/half_cuda_potential.cuf").read_text().lower()
    assembly = (ROOT / "src/half_cuda_assembly.cuf").read_text().lower()
    solver = (ROOT / "src/half_cuda_solver.cuf").read_text().lower()
    assert "use cufft" in potential
    assert "evaluate_pbe_kernel" in potential
    assert "projector_kernel" in assembly
    assert "assemble_dense_gamma_cuda_full" in assembly
    assert "cusolverdnzhegvd" in solver


def test_cuda_mimic_us_pipeline_is_device_resident_and_cli_enabled():
    uspp = (ROOT / "src/half_cuda_uspp.cuf").read_text().lower()
    potential = (ROOT / "src/half_cuda_potential.cuf").read_text().lower()
    assembly = (ROOT / "src/half_cuda_assembly.cuf").read_text().lower()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "sample_spheres_kernel" in uspp
    assert "project_vlm_kernel" in uspp
    assert "radial_v_kernel" in uspp
    assert "finish_dij_kernel" in uspp
    assert "reorder_chgcar_f_to_c_kernel" in potential
    assert "build_uspp_dij_cuda_species" in assembly
    assert "case('--uspp-dij')" in cli
    assert "--uspp-dij is not implemented" not in cli


def test_cpu_mimic_us_pipeline_is_available():
    source = (ROOT / "src/half_uspp.F90").read_text().lower()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "build_uspp_dij_cpu" in source
    assert "cubic_sample" in source
    assert "call build_uspp_dij_cpu" in cli
    assert "--uspp-dij requires a cuda build" not in cli


def test_arbitrary_kpoint_cli_is_enabled():
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "case('--kpoint')" in cli
    assert "build_plane_wave_basis(crystal,rho%shape,encut,kpoint,basis)" in cli


def test_kmesh_bands_and_total_energy_are_native_fortran_features():
    kpoints = (ROOT / "src/half_kpoints.F90").read_text().lower()
    energy = (ROOT / "src/half_energy.F90").read_text().lower()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "spg_get_ir_reciprocal_mesh" in kpoints
    assert "read_explicit_kpoints" in kpoints
    assert "compute_occupations" in energy
    assert "ewald_energy" in energy
    assert "call compute_occupations" in cli
    assert "internal_energy=band_energy-eh-exv+exc+ewald" in cli
    assert "reserved for happy-compatible" not in cli


def test_cpu_mpi_distributes_kpoints_and_has_validated_scaling():
    import json

    cmake = (ROOT / "CMakeLists.txt").read_text()
    parallel = (ROOT / "src/half_parallel.F90").read_text().lower()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    record = json.loads((ROOT / "docs/validation/hfo2_cpu_mpi.json").read_text())
    assert "HALF_ENABLE_MPI" in cmake
    assert "mpi_allreduce" in parallel
    assert "parallel_owns" in parallel
    assert "if(.not.parallel_owns(ik))cycle" in cli
    assert record["max_abs_eigenvalue_error_eV"]["2_ranks"] == 0.0
    assert record["max_abs_eigenvalue_error_eV"]["4_ranks"] == 0.0
    assert record["speedup"]["4_ranks"] > 3.8


def test_si_total_energy_parity_record():
    import json

    record = json.loads((ROOT / "docs/validation/si_total_energy_parity.json").read_text())
    assert record["max_abs_component_error_eV"] < 1e-10
    assert record["total_energy_abs_error_eV"] < 1e-10


def test_native_hdf5_input_and_vaspwave_output_are_enabled():
    bridge = (ROOT / "src/half_hdf5_bridge.c").read_text().lower()
    writer = (ROOT / "src/half_vaspwave.F90").read_text().lower()
    chgcar = (ROOT / "src/half_chgcar.F90").read_text().lower()
    potcar = (ROOT / "src/half_potcar.F90").read_text().lower()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "half_write_vaspwave_h5_c" in bridge
    assert "half_probe_vaspwave_c" in bridge
    assert "input/potcar/content" in bridge
    assert "vasp_permutation" in writer
    assert "read_vaspwave_h5" in chgcar
    assert "extract_hdf5_potcar" in potcar
    assert "case('--vaspwave-h5')" in cli


def test_si_finite_difference_force_parity_record():
    import json

    record = json.loads((ROOT / "docs/validation/si_force_parity.json").read_text())
    assert record["max_abs_force_error_eV_per_Angstrom"] < 1e-8
    assert record["speedup"] > 40.0
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "case('--forces')" in cli
    assert "evaluate_free_energy" in cli


def test_si_automatic_band_path_parity_record():
    import json

    record = json.loads((ROOT / "docs/validation/si_band_path_parity.json").read_text())
    assert record["path"] == "GXWKGLUWLK,UX"
    assert record["kpoint_max_abs_error"] == 0.0
    assert record["max_abs_eigenvalue_error_eV"] < 1e-10
    source = (ROOT / "src/half_kpoints.F90").read_text().lower()
    assert "generate_cubic_band_path" in source


def test_hfo2_reference_eigenval_energy_parity_record():
    import json

    record = json.loads(
        (ROOT / "docs/validation/hfo2_reference_eigenval_energy.json").read_text()
    )
    assert record["kpoints"] == 36
    assert record["internal_energy_abs_error_eV"] < 1e-9
    assert record["max_abs_component_error_eV"] < 1e-9
    energy = (ROOT / "src/half_energy.F90").read_text().lower()
    assert "read_vasp_eigenval" in energy


def test_potcar_structure_order_is_validated():
    source = (ROOT / "src/half_potcar.F90").read_text().lower()
    cli = (ROOT / "app/half_cli.F90").read_text().lower()
    assert "validate_potcar_structure" in source
    assert "dataset order does not match" in source
    assert "call validate_potcar_structure" in cli


def test_rtx4090_full_workflow_record():
    import json

    record = json.loads((ROOT / "docs/validation/rtx4090_full_workflow.json").read_text())
    assert record["gpu"] == "NVIDIA GeForce RTX 4090"
    assert record["ctest"]["failed"] == 0
    assert record["si_energy"]["abs_error_eV"] < 1e-10
    assert record["si_automatic_bands"]["path"] == "GXWKGLUWLK,UX"
    assert record["vaspwave_h5"]["wave_shape"] == [8, 259, 2]


def test_si_arbitrary_kpoint_parity_record():
    import json

    record = json.loads(
        (ROOT / "docs/validation/si_arbitrary_kpoint.json").read_text()
    )
    assert record["plane_waves"] == 733
    assert record["max_abs_eigenvalue_error_eV"] < 1e-10


def test_hfo2_cuda_mimic_us_parity_record():
    import json

    record = json.loads(
        (ROOT / "docs/validation/hfo2_cuda_uspp_parity.json").read_text()
    )
    assert record["plane_waves"] == 3407
    assert record["bands_compared"] == 60
    assert record["max_abs_eigenvalue_error_eV"] < 1e-10
    assert record["speedup_vs_happy_single_core"] > 30.0


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
