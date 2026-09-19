#ifndef REFERENCEBATCH_HPP
#define REFERENCEBATCH_HPP

#include "Record.hpp"
#include "types.hpp"

#include <algorithm> // for transform, for_each
#include <iterator>  // for distance

namespace vaspml
{
/*******************************************************************************************
 * Class to collect entries from an ML_AB Record to form batches over structures
 * to be used for a batch learning.
 * A std::vector<ReferenceBatch> can be set up by the following code example.
 * @code
 * std::vector<ReferenceBatch> batches( nFrames );
 * Vec1ShCRec structures = dataMLAB->vcget<Vec1ShRec>( "structures" );
 * Size nStrucs = structures.size();
 * for ( Size batchNum = 0; batchNum < nFrames; batchNum++ )
 * {
 *     Size start  =  batchNum*batchSize;
 *     Size end    =  batchSize*(batchNum+1);
 *     if ( nStrucs < end )
 *     {
 *         end = nStrucs;
 *     }
 *     batches[batchNum].initData( structures.begin()+start,
 *                                 structures.begin()+end );
 * }
 * @endcode
 *******************************************************************************************/
class ReferenceBatch
{
  public:
    /*******************************************************************************************
     * Create ReferenceBatch class
     *
     * @param dataIn user can supply ShRec for adding needed memory or using already
     * prepared memory.
     *
     * If no dataIn is supplied the class will create it's own memory.
     *******************************************************************************************/
    ReferenceBatch(const ShRec& dataIn = nullptr);
    /*******************************************************************************************
     * Used for filling data arrays with the data needed from the corresponding batch.
     *
     * @param startIter is the starting iterator which has to contain a ShCRec
     * @param startIter is the end iterator which has to contain a ShCRec
     *
     * Both iterators can be obtained from a Vec1ShCRec with the cbegin() and cend() command.
     *******************************************************************************************/
    template<typename Iterator>
    void initData(const Iterator startIter, const Iterator endIter);
    /*******************************************************************************************
     * Get the number of force entries which will 3 * total number in the batch.
     *******************************************************************************************/
    Size getNForceEntries();
    /*******************************************************************************************
     * Get number of energy entries.
     *
     * This will equal the number of structures from which batch
     * was created.
     *******************************************************************************************/
    Size getNEnergyEntries();
    /*******************************************************************************************
     * Get force vector.
     *******************************************************************************************/
    const Vec1Real& get_forces();
    /*******************************************************************************************
     * Get energy vector.
     *******************************************************************************************/
    const Vec1Real& get_energies();

  private:
    /*******************************************************************************************
     * Record where the data arrays are created.
     *******************************************************************************************/
    ShRec data;
    /*******************************************************************************************
     * Storage for force values over the whole batch. Forces are ordered according to the batches.
     *******************************************************************************************/
    Vec1Real& forces;
    /*******************************************************************************************
     * Storage array for energy entries.
     *******************************************************************************************/
    Vec1Real& energies;
};

template<typename Iterator>
void ReferenceBatch::initData(const Iterator startIter, const Iterator endIter)
{
    energies.resize(std::distance(startIter, endIter));
    std::transform(startIter,
                   endIter,
                   energies.begin(),
                   [](const ShCRec& struc) { return struc->cget<Real>("energy"); });
    std::for_each(startIter,
                  endIter,
                  [&](const ShCRec& struc)
                  {
                      const auto& values = struc->cget<Vec1Real>("forces");
                      forces.insert(forces.end(), values.cbegin(), values.cend());
                  });

    return;
}

ShRec makeDataMapReferenceBatch(const ShRec& data);

} // namespace vaspml

#endif
