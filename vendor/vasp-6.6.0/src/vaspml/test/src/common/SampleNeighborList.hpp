#ifndef SAMPLENEIGHBORLIST_HPP
#define SAMPLENEIGHBORLIST_HPP

#include "SampleStructure.hpp"

#include "nearest_neighbor.hpp"
#include "types.hpp"

namespace vaspml
{

/// Sample neighbor list for given sample structure.
struct SampleNeighborList : public NearestNeighborNSquare
{
    /** Set up neighbor list for the given sample structure.
     *
     *  @param sample String representation of sample structure, see #SampleStructure.
     *
     *  The neighbor list is automatically computed in the constructor.
     */
    SampleNeighborList(Real cutoff, bool typeSort, bool distSort, String sample);

    SampleStructure structure;
};

} //namespace vaspml

#endif
