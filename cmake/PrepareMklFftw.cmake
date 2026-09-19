if(NOT DEFINED HALF_MKL_ROOT OR NOT DEFINED HALF_VASP_FFTW_ROOT OR NOT DEFINED HALF_VASP_MAKE)
  message(FATAL_ERROR "HALF_MKL_ROOT, HALF_VASP_FFTW_ROOT, and HALF_VASP_MAKE are required")
endif()

set(_fftw_interface "${HALF_MKL_ROOT}/share/mkl/interfaces/fftw3xf")
set(_fftw_archive "${HALF_VASP_FFTW_ROOT}/lib/libfftw3.a")
if(NOT EXISTS "${_fftw_interface}/makefile")
  message(FATAL_ERROR "oneMKL FFTW3 wrapper sources not found under ${_fftw_interface}")
endif()

file(MAKE_DIRECTORY "${HALF_VASP_FFTW_ROOT}/lib" "${HALF_VASP_FFTW_ROOT}/include")
if(NOT EXISTS "${_fftw_archive}")
  execute_process(
    COMMAND "${HALF_VASP_MAKE}" -C "${_fftw_interface}" libintel64
      compiler=gnu MKLROOT=${HALF_MKL_ROOT}
      INSTALL_DIR=${HALF_VASP_FFTW_ROOT}/lib INSTALL_LIBNAME=libfftw3.a
    RESULT_VARIABLE _fftw_result)
  if(NOT _fftw_result EQUAL 0)
    message(FATAL_ERROR "Building the oneMKL FFTW3 wrapper failed: ${_fftw_result}")
  endif()
endif()

if(NOT EXISTS "${HALF_VASP_FFTW_ROOT}/lib/libfftw3_omp.a")
  file(CREATE_LINK "libfftw3.a" "${HALF_VASP_FFTW_ROOT}/lib/libfftw3_omp.a" SYMBOLIC COPY_ON_ERROR)
endif()
file(COPY "${HALF_MKL_ROOT}/include/fftw/fftw3.f" DESTINATION "${HALF_VASP_FFTW_ROOT}/include")
