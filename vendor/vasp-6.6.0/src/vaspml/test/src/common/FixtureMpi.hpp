#ifndef FIXTUREMPI_HPP
#define FIXTUREMPI_HPP

#include "mpi.hpp"

#include "Fixture.hpp"

#include "types.hpp"

namespace vaspml
{

struct FixtureMpi : public Fixture
{
    FixtureMpi();
    ~FixtureMpi();

    static Int rank;
    static Int numberRanks;
    static MPI_Comm communicator;
};

}

#endif
