#include "half.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  char error[512];
  half_handle context;
  int64_t npw;
  const double gamma[3] = {0.0, 0.0, 0.0};
  double eigenvalues[8], smin, smax;
  int status, backend = HALF_BACKEND_AUTO;
  if (argc < 3 || argc > 4) {
    fprintf(stderr, "usage: %s CHGCAR POTCAR [auto|cpu|cuda]\n", argv[0]);
    return 2;
  }
  if (argc == 4) {
    if (strcmp(argv[3], "cpu") == 0) backend = HALF_BACKEND_CPU;
    else if (strcmp(argv[3], "cuda") == 0) backend = HALF_BACKEND_CUDA;
    else if (strcmp(argv[3], "auto") != 0) return 2;
  }
  if (backend == HALF_BACKEND_AUTO)
    status = half_create_from_files(argv[1], argv[2], 200.0, HALF_XC_PBE, 1,
                                    &context, error, sizeof error);
  else
    status = half_create_from_files_backend(
        argv[1], argv[2], 200.0, HALF_XC_PBE, 1, backend, HALF_SOLVER_EVD,
        &context, error, sizeof error);
  if (status != HALF_SUCCESS) {
    fprintf(stderr, "HALF create failed (%d): %s\n", status, error);
    return status;
  }
  status = half_get_basis_size(context, gamma, &npw, error, sizeof error);
  if (status == HALF_SUCCESS)
    status = half_solve_kpoint(context, gamma, 8, eigenvalues, NULL, npw,
                               &smin, &smax, error, sizeof error);
  if (status != HALF_SUCCESS) {
    fprintf(stderr, "HALF solve failed (%d): %s\n", status, error);
    half_destroy(context, error, sizeof error);
    return status;
  }
  printf("npw=%lld overlap=[%.12g, %.12g]\n", (long long)npw, smin, smax);
  for (int band = 0; band < 8; ++band)
    printf("band %d: %.15g eV\n", band + 1, eigenvalues[band]);
  return half_destroy(context, error, sizeof error);
}
