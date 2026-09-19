#ifndef DEEPAW_HALF_H
#define DEEPAW_HALF_H

#include <complex.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HALF_ABI_VERSION 1
#define HALF_VERSION_MAJOR 0
#define HALF_VERSION_MINOR 5
#define HALF_VERSION_PATCH 0

enum half_status {
  HALF_SUCCESS = 0,
  HALF_ERROR_INVALID_ARGUMENT = 1,
  HALF_ERROR_IO = 2,
  HALF_ERROR_INVALID_HANDLE = 3,
  HALF_ERROR_CAPACITY = 4,
  HALF_ERROR_UNAVAILABLE = 5,
  HALF_ERROR_INTERNAL = 6
};

enum half_xc { HALF_XC_LDA = 1, HALF_XC_PBE = 2 };
enum half_backend { HALF_BACKEND_AUTO = 0, HALF_BACKEND_CPU = 1, HALF_BACKEND_CUDA = 2 };
enum half_solver { HALF_SOLVER_EVD = 1, HALF_SOLVER_EVJ = 2 };

enum half_capability {
  HALF_CAP_CPU = 1 << 0,
  HALF_CAP_CUDA = 1 << 1,
  HALF_CAP_MPI = 1 << 2,
  HALF_CAP_HDF5 = 1 << 3,
  HALF_CAP_SPGLIB = 1 << 4
};

typedef int64_t half_handle;

int half_get_abi_version(void);
const char *half_get_version_string(void);
int half_get_capabilities(void);

int half_create_from_files(const char *charge_path, const char *potential_path,
                           double encut_eV, int xc, int use_uspp_dij,
                           half_handle *handle, char *error, int error_capacity);
/* Extensible constructor with explicit execution policy.  The original
 * constructor is equivalent to AUTO + EVD and remains ABI-stable. */
int half_create_from_files_backend(const char *charge_path,
                                   const char *potential_path,
                                   double encut_eV, int xc, int use_uspp_dij,
                                   int backend, int solver, half_handle *handle,
                                   char *error, int error_capacity);
int half_destroy(half_handle handle, char *error, int error_capacity);

int half_get_system_info(half_handle handle, int *nions, int *ntypes,
                         int grid[3], double lattice[9], double *volume_A3,
                         char *error, int error_capacity);

/* Query/fill the k-point plane-wave basis. G vectors are packed as
 * gvec[3*i + axis]; kinetic_eV[i] has npw entries. */
int half_get_basis_size(half_handle handle, const double kpoint[3], int64_t *npw,
                        char *error, int error_capacity);
int half_get_basis(half_handle handle, const double kpoint[3], int64_t capacity,
                   int32_t *gvec, double *kinetic_eV,
                   char *error, int error_capacity);

/* Matrices and state blocks use Fortran/BLAS column-major storage. */
int half_assemble_hs(half_handle handle, const double kpoint[3], int64_t ld,
                     double _Complex *h, double _Complex *s,
                     char *error, int error_capacity);
int half_apply_hs(half_handle handle, const double kpoint[3], int nstates,
                  int64_t ld, const double _Complex *psi,
                  double _Complex *hpsi, double _Complex *spsi,
                  char *error, int error_capacity);
int half_solve_kpoint(half_handle handle, const double kpoint[3], int nbands,
                      double *eigenvalues_eV, double _Complex *eigenvectors,
                      int64_t ld, double *overlap_min, double *overlap_max,
                      char *error, int error_capacity);

/* Solve and return eigenvectors in a host application's G-vector order.
 * host_gvec is packed as host_gvec[3*i + axis].  The host and HALF bases
 * must contain exactly the same G vectors; ordering may differ. */
int half_solve_kpoint_mapped(half_handle handle, const double kpoint[3],
                            int nbands, int64_t host_npw,
                            const int32_t *host_gvec,
                            double *eigenvalues_eV,
                            double _Complex *eigenvectors, int64_t ld,
                            double *overlap_min, double *overlap_max,
                            char *error, int error_capacity);

#ifdef __cplusplus
}
#endif
#endif
