include_guard(GLOBAL)

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
