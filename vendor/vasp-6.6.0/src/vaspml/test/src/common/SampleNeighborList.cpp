#include "SampleNeighborList.hpp"

using namespace vaspml;

SampleNeighborList::SampleNeighborList(Real cutoff, bool typeSort, bool distSort, String sample) :
    NearestNeighborNSquare(cutoff, typeSort, distSort),
    structure(sample)
{
    if (structure.isDirect()) computeNearestNeighborsDirectCoordinates(structure);
    else computeNearestNeighborsCartesianCoordinates(structure);
}
