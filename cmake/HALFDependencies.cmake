include_guard(GLOBAL)

set(HALF_SPGLIB_ROOT "" CACHE PATH "spglib installation prefix")
if(NOT HALF_SPGLIB_ROOT AND DEFINED ENV{CONDA_PREFIX})
  set(HALF_SPGLIB_ROOT "$ENV{CONDA_PREFIX}" CACHE PATH "spglib installation prefix" FORCE)
endif()

set(HALF_HDF5_ROOT "" CACHE PATH "HDF5 Fortran installation prefix")
if(NOT HALF_HDF5_ROOT AND DEFINED ENV{CONDA_PREFIX})
  set(HALF_HDF5_ROOT "$ENV{CONDA_PREFIX}" CACHE PATH "HDF5 Fortran installation prefix" FORCE)
endif()
find_path(HALF_HDF5_INCLUDE_DIR hdf5.h HINTS "${HALF_HDF5_ROOT}/include")
find_library(HALF_HDF5_LIBRARY NAMES hdf5 HINTS "${HALF_HDF5_ROOT}/lib")
if(HALF_HDF5_INCLUDE_DIR AND HALF_HDF5_LIBRARY)
  set(HALF_HAVE_HDF5 TRUE)
  message(STATUS "HALF HDF5: ${HALF_HDF5_LIBRARY}")
else()
  set(HALF_HAVE_HDF5 FALSE)
  message(STATUS "HALF HDF5 disabled: set HALF_HDF5_ROOT for vaspwave.h5 support")
endif()
find_path(HALF_SPGLIB_INCLUDE_DIR spglib.h
  HINTS "${HALF_SPGLIB_ROOT}/include")
find_library(HALF_SPGLIB_LIBRARY NAMES symspg
  HINTS "${HALF_SPGLIB_ROOT}/lib")
if(HALF_SPGLIB_INCLUDE_DIR AND HALF_SPGLIB_LIBRARY)
  set(HALF_HAVE_SPGLIB TRUE)
  message(STATUS "HALF spglib: ${HALF_SPGLIB_LIBRARY}")
else()
  set(HALF_HAVE_SPGLIB FALSE)
  message(STATUS "HALF spglib disabled: set HALF_SPGLIB_ROOT for symmetry reduction")
endif()

set(HALF_MKL_ROOT "" CACHE PATH "Intel oneMKL installation prefix")
if(NOT HALF_MKL_ROOT AND DEFINED ENV{MKLROOT})
  set(HALF_MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "Intel oneMKL installation prefix" FORCE)
endif()

set(_half_mkl_hints)
if(HALF_MKL_ROOT)
  list(APPEND _half_mkl_hints
    "${HALF_MKL_ROOT}/lib"
    "${HALF_MKL_ROOT}/lib/intel64")
endif()

find_library(HALF_MKL_RT
  NAMES mkl_rt
  HINTS ${_half_mkl_hints}
)

if(HALF_MKL_RT)
  get_filename_component(HALF_MKL_LIBRARY_DIR "${HALF_MKL_RT}" DIRECTORY)
  set(HALF_HAVE_MKL TRUE)
  message(STATUS "HALF dense solver: ${HALF_MKL_RT}")
else()
  set(HALF_HAVE_MKL FALSE)
  message(STATUS "HALF dense solver disabled: set MKLROOT or HALF_MKL_ROOT")
endif()
