#include "half.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  char error[512];
  half_handle context = 0;
  const double gamma[3] = {0.0, 0.0, 0.0};
  int backend = HALF_BACKEND_AUTO, status;
  int64_t npw;
  int32_t *gvec = NULL, *reversed = NULL;
  double *kinetic = NULL, *eigenvalues = NULL;
  double _Complex *vectors = NULL;
  double smin, smax;
  const int nbands = 8;

  if (argc < 3 || argc > 4) {
    fprintf(stderr, "usage: %s CHGCAR POTCAR [auto|cpu|cuda]\n", argv[0]);
    return 2;
  }
  if (argc == 4) {
    if (strcmp(argv[3], "cpu") == 0) backend = HALF_BACKEND_CPU;
    else if (strcmp(argv[3], "cuda") == 0) backend = HALF_BACKEND_CUDA;
    else if (strcmp(argv[3], "auto") != 0) return 2;
  }
  status = half_create_from_files_backend(
      argv[1], argv[2], 200.0, HALF_XC_PBE, 1, backend, HALF_SOLVER_EVD,
      &context, error, sizeof error);
  if (status != HALF_SUCCESS) goto fail;
  status = half_get_basis_size(context, gamma, &npw, error, sizeof error);
  if (status != HALF_SUCCESS) goto fail;

  gvec = malloc((size_t)(3 * npw) * sizeof *gvec);
  reversed = malloc((size_t)(3 * npw) * sizeof *reversed);
  kinetic = malloc((size_t)npw * sizeof *kinetic);
  eigenvalues = malloc((size_t)nbands * sizeof *eigenvalues);
  vectors = malloc((size_t)(npw * nbands) * sizeof *vectors);
  if (!gvec || !reversed || !kinetic || !eigenvalues || !vectors) {
    status = HALF_ERROR_INTERNAL;
    snprintf(error, sizeof error, "allocation failed");
    goto fail;
  }
  status = half_get_basis(context, gamma, npw, gvec, kinetic, error,
                          sizeof error);
  if (status != HALF_SUCCESS) goto fail;
  for (int64_t i = 0; i < npw; ++i)
    for (int axis = 0; axis < 3; ++axis)
      reversed[3 * i + axis] = gvec[3 * (npw - 1 - i) + axis];

  status = half_solve_kpoint_mapped(
      context, gamma, nbands, npw, reversed, eigenvalues, vectors, npw,
      &smin, &smax, error, sizeof error);
  if (status != HALF_SUCCESS) goto fail;
  for (int band = 0; band < nbands; ++band)
    if (!isfinite(eigenvalues[band])) {
      status = HALF_ERROR_INTERNAL;
      snprintf(error, sizeof error, "non-finite mapped eigenvalue");
      goto fail;
    }
  printf("mapped reverse-order solve: npw=%lld overlap=[%.12g, %.12g]\n",
         (long long)npw, smin, smax);
  for (int band = 0; band < nbands; ++band)
    printf("band %d: %.15g eV\n", band + 1, eigenvalues[band]);

  free(vectors); free(eigenvalues); free(kinetic); free(reversed); free(gvec);
  return half_destroy(context, error, sizeof error);

fail:
  fprintf(stderr, "HALF mapped example failed (%d): %s\n", status, error);
  free(vectors); free(eigenvalues); free(kinetic); free(reversed); free(gvec);
  if (context) half_destroy(context, error, sizeof error);
  return status ? status : 1;
}
