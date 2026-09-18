cmake_minimum_required(VERSION 3.24)

project(HALF VERSION 0.4.0 LANGUAGES Fortran)
include(GNUInstallDirs)

option(HALF_ENABLE_CUDA "Build the CUDA Fortran backend" ON)
option(HALF_ENABLE_GPU_DEBUG "Enable CUDA bounds checks and synchronous launches" OFF)
option(HALF_ENABLE_OPENMP "Parallelize the CPU backend with OpenMP" ON)
set(HALF_GPU_CC "86" CACHE STRING "NVIDIA GPU compute capability, without the decimal point")
set(HALF_CUDA_VERSION "12.4" CACHE STRING "CUDA toolkit selected by NVHPC")

if(HALF_ENABLE_CUDA AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  message(FATAL_ERROR
    "The HALF CUDA backend requires NVIDIA HPC SDK (nvfortran).")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_Fortran_PREPROCESS ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(HALFDependencies)
include(HALFTargetOptions)

set(HALF_CORE_SOURCES
  src/half_kinds.F90
  src/half_constants.F90
  src/half_math.F90
  src/half_types.F90
  src/half_chgcar.F90
  src/half_potcar.F90
  src/half_basis.F90
  src/half_paw.F90
  src/half_fft.F90
  src/half_xc.F90
  src/half_potential.F90
  src/half_dense_solver.F90
  src/half_cpu.F90
)
if(HALF_ENABLE_CUDA)
  list(APPEND HALF_CORE_SOURCES src/half_cuda.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_potential.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_assembly.cuf)
  list(APPEND HALF_CORE_SOURCES src/half_cuda_solver.cuf)
  set_source_files_properties(src/half_cuda.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_potential.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_assembly.cuf PROPERTIES LANGUAGE Fortran)
  set_source_files_properties(src/half_cuda_solver.cuf PROPERTIES LANGUAGE Fortran)
endif()

add_library(half_core ${HALF_CORE_SOURCES})
set_target_properties(half_core PROPERTIES
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  target_compile_options(half_core PRIVATE -O3 -Mbackslash -Mextend)
  if(HALF_ENABLE_OPENMP)
    target_compile_options(half_core PRIVATE -mp)
    target_link_options(half_core PUBLIC -mp)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(half_core PRIVATE -O3 -ffree-line-length-none)
  if(HALF_ENABLE_OPENMP)
    target_compile_options(half_core PRIVATE -fopenmp)
    target_link_options(half_core PUBLIC -fopenmp)
  endif()
endif()

if(HALF_ENABLE_CUDA)
  target_compile_options(half_core PRIVATE
    -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC},lineinfo)
  target_link_options(half_core PUBLIC
    -cuda -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
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
endif()

install(TARGETS half_core half half-inspect)
if(UNIX)
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink half \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/half-validate-gamma\")")
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink half \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/half-bands\")")
  install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink half \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/half-energy\")")
endif()
