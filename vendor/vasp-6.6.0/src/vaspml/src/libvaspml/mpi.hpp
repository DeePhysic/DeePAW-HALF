#ifndef MPI_HPP
#define MPI_HPP

// OpenMPI: avoid loading deprecated C++ bindings.
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

// MPICH: avoid loading deprecated C++ bindings.
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#include <mpi.h> // IWYU pragma: export

#endif
