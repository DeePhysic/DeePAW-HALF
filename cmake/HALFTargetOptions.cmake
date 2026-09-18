include_guard(GLOBAL)

function(half_set_fortran_options target)
  set_target_properties(${target} PROPERTIES
    Fortran_STANDARD 2008
    Fortran_STANDARD_REQUIRED YES)

  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    target_compile_options(${target} PRIVATE -O3 -Mextend)
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE -O3 -ffree-line-length-none)
  endif()
endfunction()

function(half_link_mkl target)
  if(NOT HALF_HAVE_MKL)
    message(FATAL_ERROR "half_link_mkl called without a discovered oneMKL runtime")
  endif()
  target_link_libraries(${target} PRIVATE "${HALF_MKL_RT}" pthread m dl)
  set_target_properties(${target} PROPERTIES BUILD_RPATH "${HALF_MKL_LIBRARY_DIR}")
endfunction()
