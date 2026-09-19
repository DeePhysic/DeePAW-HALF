#include "half.h"
#include <stdio.h>
#include <string.h>

int main(void) {
  char error[128];
  half_handle handle = 0;
  if (half_get_abi_version() != HALF_ABI_VERSION) return 1;
  if (strcmp(half_get_version_string(), "0.5.0") != 0) return 2;
  if (!(half_get_capabilities() & (HALF_CAP_CPU | HALF_CAP_CUDA))) return 3;
  if (half_destroy(12345, error, sizeof error) != HALF_ERROR_INVALID_HANDLE) return 4;
  if (strstr(error, "invalid") == NULL) return 5;
  if (half_create_from_files_backend("", "", 400.0, HALF_XC_PBE, 1,
                                     99, HALF_SOLVER_EVD, &handle, error,
                                     sizeof error) != HALF_ERROR_INVALID_ARGUMENT)
    return 6;
  puts("HALF C API 0.5.0: ok");
  return 0;
}
