#include "half.h"
#include <stdio.h>
#include <string.h>

int main(void) {
  char error[128];
  half_handle handle = 0;
  const double gamma[3] = {0.0, 0.0, 0.0};
  const int32_t gzero[3] = {0, 0, 0};
  double eigenvalue = 0.0, smin = 0.0, smax = 0.0;
  double _Complex eigenvector = 0.0;
  if (half_get_abi_version() != HALF_ABI_VERSION) return 1;
  if (strcmp(half_get_version_string(), "0.6.0") != 0) return 2;
  if (half_get_capabilities() & ~(HALF_CAP_CPU | HALF_CAP_CUDA | HALF_CAP_MPI |
                                  HALF_CAP_HDF5 | HALF_CAP_SPGLIB |
                                  HALF_CAP_ESCN_API)) return 3;
  if (half_destroy(12345, error, sizeof error) != HALF_ERROR_INVALID_HANDLE) return 4;
  if (strstr(error, "invalid") == NULL) return 5;
  if (half_solve_kpoint_mapped(12345, gamma, 1, 1, gzero, &eigenvalue,
                               &eigenvector, 1, &smin, &smax, error,
                               sizeof error) != HALF_ERROR_INVALID_HANDLE)
    return 6;
  if (half_create_from_files_backend("", "", 400.0, HALF_XC_PBE, 1,
                                     99, HALF_SOLVER_EVD, &handle, error,
                                     sizeof error) != HALF_ERROR_INVALID_ARGUMENT)
    return 7;
  if (half_set_request_geometry(12345, NULL, error, sizeof error) !=
      HALF_ERROR_INVALID_HANDLE)
    return 8;
  if (half_get_request_geometry(12345, NULL, NULL, NULL, NULL, 0, NULL, 0,
                                NULL, error, sizeof error) !=
      HALF_ERROR_INVALID_HANDLE)
    return 9;
  puts("HALF C API 0.6.0: ok");
  return 0;
}
