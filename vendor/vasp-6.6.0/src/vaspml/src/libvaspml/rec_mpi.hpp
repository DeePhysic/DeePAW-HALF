#ifndef REC_MPI_HPP
#define REC_MPI_HPP

#include "types.hpp"

namespace vaspml
{
class MlMPI;
struct Record;
}

namespace vaspml::rec::mpi
{

void      bcast(Record& record, MlMPI& mlmpi, const Int root = 0);
Vec1ShRec allGatherVVector(Record& send, MlMPI& mlmpi);
ShRec     allGatherV(Record& send, MlMPI& mlmpi);
ShRec     gatherV(Record& data, MlMPI& mlmpi, const Int root = 0);

} // namespace vaspml::rec::mpi

#endif // REC_MPI_HPP
