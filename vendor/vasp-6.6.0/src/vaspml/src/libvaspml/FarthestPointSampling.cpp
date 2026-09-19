#include "FarthestPointSampling.hpp"

#include "Random.hpp"

#include <random>

using namespace vaspml;

Int fps::chooseFirstIndexOMP(const Size nPoints, const Size seedIn)
{
    Size                      seedValue = (seedIn == 0) ? std::random_device{}() : seedIn;
    UniformRandomIntGenerator rand(seedValue);
    Int                       firstIndex = rand.uniformInt(0, nPoints - 1);

    return firstIndex;
}

Int fps::chooseFirstIndexMPI(const Size                          nPoints,
                             const Size                          seedIn,
                             const std::shared_ptr<const MlMPI>& mlmpi)
{
    Size seedValue;
    if (seedIn == 0 and mlmpi and mlmpi->get_rank() == 0)
    {
        std::random_device rd;
        seedValue = rd();
    }
    else seedValue = seedIn;
    UniformRandomIntGenerator rand(seedValue);
    Int                       firstIndex = rand.uniformInt(0, nPoints - 1);
    if (mlmpi) mlmpi->bcast(firstIndex, 0);

    return firstIndex;
}
