#ifndef DEEPAW_HALF_H
#define DEEPAW_HALF_H

#include <complex.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HALF_ABI_VERSION 1
#define HALF_VERSION_MAJOR 0
#define HALF_VERSION_MINOR 6
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
  HALF_CAP_SPGLIB = 1 << 4,
  HALF_CAP_ESCN_API = 1 << 5
};

enum half_escn_flag { HALF_ESCN_INCLUDE_UNCERTAINTY = 1 << 0 };
enum half_density_normalization {
  HALF_DENSITY_RAW = 0,
  HALF_DENSITY_NORMALIZE_VALENCE = 1
};

typedef int64_t half_handle;

/* VASP-to-HALF request metadata. ABI v1 stores a deep copy in the context and
 * is also the input geometry for half_create_from_escn(). lattice is
 * row-major, positions_fractional is atom-major [nions][3], and
 * species contains one-based type indices. Set struct_size to sizeof this
 * structure, set flags to zero, and zero every reserved field. */
typedef struct half_request_geometry_v1 {
  uint32_t struct_size;
  uint32_t flags;
  int32_t nions;
  int32_t ntypes;
  int32_t grid[3];
  int32_t reserved_i32;
  double lattice[9];
  const int32_t *species;
  const double *positions_fractional;
  uint64_t reserved[8];
} half_request_geometry_v1;

/* Direct DeePAW-eSCN wire request. Cell vectors and Cartesian positions are
 * in angstrom and stored row-major/atom-major. The API is three-dimensionally
 * periodic, so pbc is intentionally not configurable here. */
typedef struct half_escn_request_v1 {
  uint32_t struct_size;
  uint32_t flags;
  int32_t nions;
  int32_t grid[3];
  int32_t reserved_i32;
  double cell[9];
  const int32_t *atomic_numbers;
  const double *positions_cartesian;
  uint64_t reserved[8];
} half_escn_request_v1;

/* Caller-owned output buffers. capacity is measured in float elements and
 * must be at least grid[0]*grid[1]*grid[2]. Uncertainty pointers may be null
 * unless HALF_ESCN_INCLUDE_UNCERTAINTY was requested. Returned arrays retain
 * the server's C-order [nx,ny,nz] layout. */
typedef struct half_escn_result_v1 {
  uint32_t struct_size;
  uint32_t flags;
  int32_t grid[3];
  int32_t reserved_i32;
  int64_t capacity;
  float *density;
  float *nu;
  float *alpha;
  float *beta;
  float *risk;
  double elapsed_seconds;
  uint64_t reserved[8];
} half_escn_result_v1;

int half_get_abi_version(void);
const char *half_get_version_string(void);
int half_get_capabilities(void);

int half_escn_health(const char *base_url, int timeout_seconds,
                     char *error, int error_capacity);
int half_escn_predict(const char *base_url,
                      const half_escn_request_v1 *request,
                      half_escn_result_v1 *result, int timeout_seconds,
                      char *error, int error_capacity);

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
/* Obtain the smooth density from DeePAW-eSCN and create a normal HALF context.
 * geometry uses fractional positions and one-based POTCAR type indices. RAW
 * preserves the model output; NORMALIZE_VALENCE rescales its mean to the sum
 * of POTCAR ZVAL values. Negative values are never clipped. */
int half_create_from_escn(const char *base_url, const char *potential_path,
                          const half_request_geometry_v1 *geometry,
                          double encut_eV, int xc, int use_uspp_dij,
                          int backend, int solver, int normalization,
                          int timeout_seconds, half_handle *handle,
                          char *error, int error_capacity);
int half_destroy(half_handle handle, char *error, int error_capacity);

int half_set_request_geometry(half_handle handle,
                              const half_request_geometry_v1 *geometry,
                              char *error, int error_capacity);
/* Query dimensions first with null species/positions, then provide capacities
 * measured in int32_t and double elements, respectively. lattice and grid are
 * always returned when their pointers are non-null. */
int half_get_request_geometry(half_handle handle, int32_t *nions,
                              int32_t *ntypes, int32_t grid[3],
                              double lattice[9], int64_t species_capacity,
                              int32_t *species, int64_t positions_capacity,
                              double *positions_fractional,
                              char *error, int error_capacity);

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
