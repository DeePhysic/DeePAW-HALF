cmake_minimum_required(VERSION 3.24)

project(HALF VERSION 0.6.0 LANGUAGES Fortran C)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(HALF_ENABLE_CUDA "Build the CUDA Fortran backend" ON)
option(HALF_ENABLE_EVJ "Build the optional CUDA C++ cuSOLVER Jacobi bridge" OFF)
option(HALF_ENABLE_GPU_DEBUG "Enable CUDA bounds checks and synchronous launches" OFF)
option(HALF_ENABLE_OPENMP "Parallelize the CPU backend with OpenMP" ON)
option(HALF_ENABLE_MPI "Distribute CPU k points across MPI ranks" OFF)
option(HALF_BUILD_SHARED_LIBRARY "Build the reusable libhalf shared library" ON)
option(HALF_BUILD_EXAMPLES "Build libhalf API examples" ON)
option(HALF_BUILD_VASP "Build the bundled private VASP 6.6.0 integration" OFF)
option(HALF_ENABLE_ESCN_API "Enable the DeePAW-eSCN HTTP inference client" ON)
set(HALF_GPU_CC "86" CACHE STRING "NVIDIA GPU compute capability, without the decimal point")
set(HALF_CUDA_VERSION "12.4" CACHE STRING "CUDA toolkit selected by NVHPC")

if(HALF_ENABLE_CUDA AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  message(FATAL_ERROR
    "The HALF CUDA backend requires NVIDIA HPC SDK (nvfortran).")
endif()

if(HALF_ENABLE_CUDA AND HALF_ENABLE_EVJ)
  # The EVJ bridge is CUDA C++ because NVHPC does not expose Zhegvj in its
  # CUDA Fortran module.  Prefer the compiler from the selected HPC SDK.
  if(NOT CMAKE_CUDA_COMPILER)
    find_program(HALF_NVCC NAMES nvcc HINTS "$ENV{NVHPC}/cuda/${HALF_CUDA_VERSION}/bin")
    if(HALF_NVCC)
      set(CMAKE_CUDA_COMPILER "${HALF_NVCC}" CACHE FILEPATH "CUDA compiler" FORCE)
    endif()
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED YES)
  set(CMAKE_CUDA_ARCHITECTURES "${HALF_GPU_CC}")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_Fortran_PREPROCESS ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(HALFDependencies)
include(HALFTargetOptions)

set(HALF_CORE_SOURCES
  src/half_kinds.F90
  src/half_parallel.F90
  src/half_constants.F90
  src/half_math.F90
  src/half_types.F90
  src/half_chgcar.F90
  src/half_potcar.F90
  src/half_basis.F90
  src/half_kpoints.F90
  src/half_paw.F90
  src/half_energy.F90
  src/half_artifacts.F90
  src/half_artifacts.c
  src/half_library.F90
  src/half_c_api.F90
  src/half_cpu.F90
)
if(HALF_HAVE_ESCN_SOCKET)
  list(APPEND HALF_CORE_SOURCES src/half_escn_socket.c)
else()
  # libcurl implementation when available; otherwise ABI-compatible stubs.
  list(APPEND HALF_CORE_SOURCES src/half_escn_client.c)
endif()
if(HALF_HAVE_MKL)
  list(APPEND HALF_CORE_SOURCES
    src/half_fft.F90 src/half_uspp.F90 src/half_xc.F90
    src/half_potential.F90 src/half_dense_solver.F90)
endif()
if(HALF_HAVE_HDF5)
  list(APPEND HALF_CORE_SOURCES src/half_vaspwave.F90 src/half_hdf5_bridge.c)
endif()
if(HALF_ENABLE_CUDA)
  list(APPEND HALF_CORE_SOURCES src/half_cuda.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_potential.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_uspp.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_assembly.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_solver.cuf)
  set_source_files_properties(src/half_cuda.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_potential.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_uspp.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_assembly.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_solver.cuf PROPERTIES LANGUAGE Fortran)
endif()

if(HALF_BUILD_SHARED_LIBRARY)
  add_library(half_core SHARED ${HALF_CORE_SOURCES})
else()
  add_library(half_core STATIC ${HALF_CORE_SOURCES})
endif()
set_target_properties(half_core PROPERTIES
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
  OUTPUT_NAME half
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  POSITION_INDEPENDENT_CODE YES
  EXPORT_NAME half
)
target_include_directories(half_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(HALF_HAVE_CURL)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_ESCN_API=1)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_CURL=1)
  target_link_libraries(half_core PUBLIC CURL::libcurl)
elseif(HALF_HAVE_ESCN_SOCKET)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_ESCN_API=1 HALF_HAVE_ESCN_SOCKET=1)
endif()
if(HALF_HAVE_MPI)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_MPI=1)
  # Conda HDF5 prefixes can also contain an mpi_f08.mod built by another
  # Fortran compiler.  Put the selected MPI wrapper's module directory first.
  target_include_directories(half_core BEFORE PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
  target_link_libraries(half_core PUBLIC MPI::MPI_Fortran)
endif()
if(HALF_HAVE_MKL)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_MKL=1)
  half_link_mkl(half_core)
endif()
if(HALF_HAVE_SPGLIB)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_SPGLIB=1)
  target_link_libraries(half_core PUBLIC "${HALF_SPGLIB_LIBRARY}")
endif()
if(HALF_HAVE_HDF5)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_HDF5=1)
  # HALF calls HDF5 through its C ABI; keeping this include directory off
  # Fortran command lines avoids picking up incompatible mpi*.mod files from
  # mixed compiler prefixes.
  target_include_directories(half_core PRIVATE $<$<COMPILE_LANGUAGE:C>:${HALF_HDF5_INCLUDE_DIR}>)
  target_link_libraries(half_core PUBLIC "${HALF_HDF5_LIBRARY}")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  target_compile_options(half_core PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3;-Mbackslash;-Mextend>)
  if(HALF_ENABLE_OPENMP)
    target_compile_options(half_core PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-mp>)
    target_link_options(half_core PRIVATE -mp)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(half_core PRIVATE -O3 -ffree-line-length-none)
  if(HALF_ENABLE_OPENMP)
    target_compile_options(half_core PRIVATE -fopenmp)
    target_link_options(half_core PUBLIC -fopenmp)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
  target_compile_options(half_core PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)
  if(HALF_ENABLE_OPENMP)
    target_compile_options(half_core PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-qopenmp>)
    target_link_options(half_core PUBLIC $<$<LINK_LANGUAGE:Fortran>:-qopenmp>)
  endif()
endif()

if(HALF_ENABLE_CUDA)
  target_compile_definitions(half_core PRIVATE HALF_HAVE_CUDA=1)
  if(HALF_ENABLE_EVJ)
    add_library(half_evj_bridge OBJECT src/half_cuda_evj_bridge.cu)
    set_target_properties(half_evj_bridge PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES)
    get_filename_component(HALF_NVFORTRAN_BIN "${CMAKE_Fortran_COMPILER}" DIRECTORY)
    get_filename_component(HALF_NVHPC_ROOT "${HALF_NVFORTRAN_BIN}/../.." ABSOLUTE)
    target_include_directories(half_evj_bridge PRIVATE
      "${HALF_NVHPC_ROOT}/math_libs/${HALF_CUDA_VERSION}/targets/x86_64-linux/include")
    target_sources(half_core PRIVATE $<TARGET_OBJECTS:half_evj_bridge>)
    target_compile_definitions(half_core PRIVATE HALF_HAVE_EVJ=1)
  endif()
  target_compile_options(half_core PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-cuda;-gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC},lineinfo>)
  target_link_options(half_core PRIVATE
    -cuda -cudalib=cusolver,cufft -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
endif()

if(HALF_ENABLE_GPU_DEBUG)
  target_compile_options(half_core PRIVATE -Mbounds -gpu=debug)
  target_compile_definitions(half_core PRIVATE HALF_GPU_DEBUG=1)
endif()

# Unified command-line interface.  It remains buildable without oneMKL so that
# input inspection and help/version are always available; Gamma reconstruction
# is enabled automatically when MKL is discovered.
add_executable(half app/half_cli.F90)
half_set_fortran_options(half)
target_link_libraries(half PRIVATE half_core)
if(HALF_HAVE_MKL)
  target_compile_definitions(half PRIVATE HALF_CLI_HAVE_MKL=1)
  half_link_mkl(half)
  set_target_properties(half PROPERTIES INSTALL_RPATH "${HALF_MKL_LIBRARY_DIR}")
endif()
if(HALF_HAVE_HDF5)
  target_compile_definitions(half PRIVATE HALF_CLI_HAVE_HDF5=1)
endif()
if(HALF_ENABLE_CUDA)
  target_compile_definitions(half PRIVATE HALF_CLI_CUDA=1)
  target_compile_options(half PRIVATE -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC},lineinfo)
  target_link_options(half PRIVATE -cuda -cudalib=cusolver,cufft -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
endif()
if(UNIX)
  add_custom_command(TARGET half POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_NAME:half> $<TARGET_FILE_DIR:half>/half-validate-gamma
    COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_NAME:half> $<TARGET_FILE_DIR:half>/half-bands
    COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_NAME:half> $<TARGET_FILE_DIR:half>/half-energy
    VERBATIM)
endif()

if(HALF_ENABLE_CUDA)
  add_executable(half-inspect app/half_inspect.F90)
else()
  add_executable(half-inspect app/half_inspect_cpu.F90)
endif()
set_target_properties(half-inspect PROPERTIES
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
)
target_link_libraries(half-inspect PRIVATE half_core)
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  target_compile_options(half-inspect PRIVATE -O3 -Mextend)
endif()

add_executable(half-slice-benchmark app/half_slice_benchmark.F90)
set_target_properties(half-slice-benchmark PROPERTIES
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
)
target_link_libraries(half-slice-benchmark PRIVATE half_core)
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  target_compile_options(half-slice-benchmark PRIVATE -O3 -Mextend)
endif()

add_executable(half-potcar-inspect app/half_potcar_inspect.F90)
set_target_properties(half-potcar-inspect PROPERTIES
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
)
target_link_libraries(half-potcar-inspect PRIVATE half_core)

add_executable(half-paw-inspect app/half_paw_inspect.F90)
set_target_properties(half-paw-inspect PROPERTIES Fortran_STANDARD 2008 Fortran_STANDARD_REQUIRED YES)
target_link_libraries(half-paw-inspect PRIVATE half_core)

add_executable(half-energy-check app/half_energy_check.F90)
set_target_properties(half-energy-check PROPERTIES Fortran_STANDARD 2008 Fortran_STANDARD_REQUIRED YES)
target_link_libraries(half-energy-check PRIVATE half_core)

if(HALF_HAVE_MKL)
  add_executable(half-uspp-inspect app/half_uspp_inspect.F90)
  half_set_fortran_options(half-uspp-inspect)
  target_link_libraries(half-uspp-inspect PRIVATE half_core)
  half_link_mkl(half-uspp-inspect)
endif()

if(HALF_HAVE_MKL)
  add_executable(half-veff-inspect app/half_veff_inspect.F90)
  half_set_fortran_options(half-veff-inspect)
  target_link_libraries(half-veff-inspect PRIVATE half_core)
  half_link_mkl(half-veff-inspect)
  add_executable(half-gamma app/half_gamma_lda.F90)
  half_set_fortran_options(half-gamma)
  target_link_libraries(half-gamma PRIVATE half_core)
  half_link_mkl(half-gamma)
endif()
if(HALF_ENABLE_CUDA)
  target_compile_options(half-inspect PRIVATE
    -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC},lineinfo)
  target_link_options(half-inspect PRIVATE
    -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})

  add_executable(half-benchmark app/half_benchmark.F90)
  target_link_libraries(half-benchmark PRIVATE half_core)
  target_compile_options(half-benchmark PRIVATE
    -O3 -Mextend -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC},lineinfo)
  target_link_options(half-benchmark PRIVATE
    -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
  add_executable(half-gamma-cuda app/half_gamma_cuda.cuf)
  set_source_files_properties(app/half_gamma_cuda.cuf PROPERTIES LANGUAGE Fortran)
  target_link_libraries(half-gamma-cuda PRIVATE half_core)
  target_compile_options(half-gamma-cuda PRIVATE -O3 -Mextend -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
  target_link_options(half-gamma-cuda PRIVATE -cuda -cudalib=cusolver,cufft -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
endif()

include(CTest)
if(BUILD_TESTING)
  add_test(NAME half_inspect_help COMMAND half-inspect --help)
  add_test(NAME half_cli_help COMMAND half --help)
  add_test(NAME half_cli_version COMMAND half --version)
  add_test(NAME half_cli_gamma_help COMMAND half gamma --help)
  add_test(NAME half_cli_bands_help COMMAND half bands --help)
  add_test(NAME half_cli_energy_help COMMAND half energy --help)
  add_test(NAME half_energy_numerics COMMAND half-energy-check)
  add_executable(half-c-api-check tests/test_c_api.c)
  target_link_libraries(half-c-api-check PRIVATE half_core)
  if(HALF_ENABLE_CUDA)
    set_target_properties(half-c-api-check PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_options(half-c-api-check PRIVATE -Mnomain -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
  endif()
  add_test(NAME half_c_api COMMAND half-c-api-check)
  if(HALF_HAVE_ESCN_API AND UNIX)
    add_executable(half-escn-api-check tests/test_escn_client.c)
    target_link_libraries(half-escn-api-check PRIVATE half_core)
    if(HALF_ENABLE_CUDA)
      set_target_properties(half-escn-api-check PROPERTIES LINKER_LANGUAGE Fortran)
      target_link_options(half-escn-api-check PRIVATE -Mnomain -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
    endif()
    add_test(NAME half_escn_api COMMAND half-escn-api-check)
  endif()
endif()
if(HALF_BUILD_EXAMPLES)
  add_executable(half-c-example examples/api/half_c_example.c)
  target_link_libraries(half-c-example PRIVATE half_core)
  add_executable(half-mapped-example examples/api/half_mapped_example.c)
  target_link_libraries(half-mapped-example PRIVATE half_core)
  if(HALF_HAVE_ESCN_API)
    add_executable(half-escn-example examples/api/half_escn_example.c)
    target_link_libraries(half-escn-example PRIVATE half_core)
  endif()
  if(HALF_ENABLE_CUDA)
    set_target_properties(half-c-example PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_options(half-c-example PRIVATE -Mnomain -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
    set_target_properties(half-mapped-example PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_options(half-mapped-example PRIVATE -Mnomain -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
    if(TARGET half-escn-example)
      set_target_properties(half-escn-example PROPERTIES LINKER_LANGUAGE Fortran)
      target_link_options(half-escn-example PRIVATE -Mnomain -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
    endif()
  endif()
endif()

if(HALF_BUILD_VASP)
  include(HALFVaspIntegration)
endif()

configure_package_config_file(cmake/HALFConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/HALFConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/HALF")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/HALFConfigVersion.cmake"
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
configure_file(cmake/half.pc.in "${CMAKE_CURRENT_BINARY_DIR}/half.pc" @ONLY)

install(TARGETS half_core EXPORT HALFTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS half half-inspect
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES include/half.h include/half_api.f90 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/half"
  FILES_MATCHING PATTERN "*.mod")
install(EXPORT HALFTargets FILE HALFTargets.cmake NAMESPACE HALF::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/HALF")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/HALFConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/HALFConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/HALF")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/half.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
if(UNIX)
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink half \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/half-validate-gamma\")")
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink half \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/half-bands\")")
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink half \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/half-energy\")")
endif()
