include_guard(GLOBAL)

if(NOT HALF_ENABLE_CUDA)
  message(FATAL_ERROR "HALF_BUILD_VASP requires HALF_ENABLE_CUDA=ON")
endif()
if(NOT HALF_BUILD_SHARED_LIBRARY)
  message(FATAL_ERROR "HALF_BUILD_VASP requires HALF_BUILD_SHARED_LIBRARY=ON")
endif()
if(NOT HALF_HAVE_MKL OR NOT HALF_MKL_ROOT)
  message(FATAL_ERROR "HALF_BUILD_VASP requires HALF_MKL_ROOT or MKLROOT")
endif()
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  message(FATAL_ERROR "HALF_BUILD_VASP requires NVHPC nvfortran")
endif()

set(HALF_VASP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/vendor/vasp-6.6.0"
  CACHE PATH "Private VASP 6.6.0 source tree")
set(HALF_VASP_JOBS "8" CACHE STRING "Parallel jobs passed to the VASP make build")
option(HALF_VASP_ENABLE_OPENACC "Build VASP itself with its OpenACC GPU port" OFF)

if(NOT EXISTS "${HALF_VASP_SOURCE_DIR}/src/main.F")
  message(FATAL_ERROR "HALF_VASP_SOURCE_DIR does not contain src/main.F: ${HALF_VASP_SOURCE_DIR}")
endif()

find_program(HALF_VASP_MAKE NAMES gmake make REQUIRED)
find_program(HALF_VASP_MPIFC NAMES mpif90 mpifort REQUIRED)
find_program(HALF_VASP_MPICC NAMES mpicc REQUIRED)
get_filename_component(_half_nvhpc_bin "${CMAKE_Fortran_COMPILER}" DIRECTORY)
find_program(HALF_VASP_NVC NAMES nvc HINTS "${_half_nvhpc_bin}" REQUIRED NO_DEFAULT_PATH)
find_program(HALF_VASP_NVCXX NAMES nvc++ HINTS "${_half_nvhpc_bin}" REQUIRED NO_DEFAULT_PATH)

if(HALF_VASP_ENABLE_OPENACC)
  set(HALF_VASP_ARCH_INCLUDE "${HALF_VASP_SOURCE_DIR}/arch/makefile.include.nvhpc_omp_acc")
  set(HALF_VASP_COMPILE_FLAGS "-acc -mp $(GPU) -gpu=tripcount:host")
else()
  set(HALF_VASP_ARCH_INCLUDE "${HALF_VASP_SOURCE_DIR}/arch/makefile.include.nvhpc")
  set(HALF_VASP_COMPILE_FLAGS "")
endif()
if(NOT EXISTS "${HALF_VASP_ARCH_INCLUDE}")
  message(FATAL_ERROR "Missing VASP NVHPC template: ${HALF_VASP_ARCH_INCLUDE}")
endif()

set(HALF_VASP_WORK_DIR "${CMAKE_BINARY_DIR}/vasp-6.6.0")
set(HALF_VASP_FFTW_ROOT "${CMAKE_BINARY_DIR}/mkl-fftw")
file(MAKE_DIRECTORY "${HALF_VASP_WORK_DIR}")
file(COPY "${HALF_VASP_SOURCE_DIR}/" DESTINATION "${HALF_VASP_WORK_DIR}")
file(MAKE_DIRECTORY "${HALF_VASP_WORK_DIR}/bin")
configure_file("${CMAKE_CURRENT_LIST_DIR}/vasp.makefile.include.in"
  "${HALF_VASP_WORK_DIR}/makefile.include" @ONLY)

file(GLOB_RECURSE _half_vasp_sources CONFIGURE_DEPENDS
  "${HALF_VASP_SOURCE_DIR}/src/*"
  "${HALF_VASP_SOURCE_DIR}/arch/*")
set(HALF_VASP_EXECUTABLE "${HALF_VASP_WORK_DIR}/bin/vasp_std")

add_custom_command(
  OUTPUT "${HALF_VASP_EXECUTABLE}"
  COMMAND "${CMAKE_COMMAND}"
    -DHALF_MKL_ROOT=${HALF_MKL_ROOT}
    -DHALF_VASP_FFTW_ROOT=${HALF_VASP_FFTW_ROOT}
    -DHALF_VASP_MAKE=${HALF_VASP_MAKE}
    -P "${CMAKE_CURRENT_LIST_DIR}/PrepareMklFftw.cmake"
  COMMAND "${HALF_VASP_MAKE}" DEPS=1 std -j${HALF_VASP_JOBS}
  DEPENDS half_core ${_half_vasp_sources}
    "${CMAKE_CURRENT_LIST_DIR}/vasp.makefile.include.in"
    "${CMAKE_CURRENT_LIST_DIR}/PrepareMklFftw.cmake"
  WORKING_DIRECTORY "${HALF_VASP_WORK_DIR}"
  COMMENT "Building private VASP 6.6.0 with DeePAW-HALF"
  USES_TERMINAL
  VERBATIM)

add_custom_target(vasp-half ALL DEPENDS "${HALF_VASP_EXECUTABLE}")
add_dependencies(vasp-half half_core)
install(PROGRAMS "${HALF_VASP_EXECUTABLE}" TYPE BIN OPTIONAL)

message(STATUS "VASP-HALF source: ${HALF_VASP_SOURCE_DIR}")
message(STATUS "VASP-HALF build:  ${HALF_VASP_WORK_DIR}")
message(STATUS "VASP-HALF MPI:    ${HALF_VASP_MPIFC}")
message(STATUS "VASP-HALF OpenACC: ${HALF_VASP_ENABLE_OPENACC}")
