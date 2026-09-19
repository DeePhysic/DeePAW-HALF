cmake_minimum_required(VERSION 3.24)
project(HALFDownstreamExample LANGUAGES C Fortran)

find_package(HALF 0.5 CONFIG REQUIRED)
add_executable(half-downstream ../api/half_c_example.c)
target_link_libraries(half-downstream PRIVATE HALF::half)
if(HALF_WITH_CUDA)
  # NVHPC's Fortran driver supplies CUDA Fortran device-registration startup.
  set_target_properties(half-downstream PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_options(half-downstream PRIVATE -Mnomain -cuda
    -gpu=cuda${HALF_CUDA_VERSION},cc${HALF_GPU_CC})
endif()
