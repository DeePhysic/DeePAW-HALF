#ifndef BATCHTOOLS_HPP
#define BATCHTOOLS_HPP

#include "nearest_neighbor.hpp"
#include "types.hpp"

#include <iterator>

namespace vaspml
{

class AtomBatchMap;
class BatchManager;
class Descriptor;
class DescriptorCollector;
class DescriptorSHS2;
class DescriptorSHS3;
class Frame;

namespace batch_tools
{

/*******************************************************************************************
 * Extract neighbor lists of given key from frame class iterators
 *
 * @param startFrame Iterator to a std::vector<Frame> where to start
 * @param endFrame Iterator to a std::vector<Frame> where to end
 * @param key of neighbor list which should be extracted. Either 2-body or 3-body.
 *******************************************************************************************/
template<typename Iterator>
std::vector<std::shared_ptr<NearestNeighborNSquare>> makeNeighborListVectorFromFrame(
    const Iterator startFrame,
    const Iterator endFrame,
    const String   key)
{
    Size                                                 n = std::distance(startFrame, endFrame);
    std::vector<std::shared_ptr<NearestNeighborNSquare>> nLists(n);
    Iterator                                             frame = startFrame;
    for (Size i = 0; i < n; i++)
    {
        nLists[i] = std::make_shared<NearestNeighborNSquare>();
        nLists[i] = frame->get_neighborLists().at(key);
        frame++;
    }
    return nLists;
}
/*******************************************************************************************
 * Extract neighbor lists of given key from frame class vector
 *
 * @param key of neighbor list which should be extracted. Either 2-body or 3-body.
 *******************************************************************************************/
std::vector<std::shared_ptr<NearestNeighborNSquare>> makeNeighborListVectorFromFrame(
    const std::vector<Frame>& frames,
    const String              key);
/*******************************************************************************************
 * Create a NearestNeighborNSquare object from vector of frame classes.
 *
 * @param frames Vector of frame classes containing set up neighbor lists and SHS2 and SHS3
 * @param strucIdx Original structure indices-1 as in ML_AB file. Format [type][structire index]
 * @param atomIdx Original atom indices-1 as in ML_AB file. Format [type][atom index]
 * @param atomBatchMap by supplying the original structure and atom index of ML_AB
 * class will give back the atom index in the corresponding type ordered frame.
 * @param strucToBatch mapping of the original ML_AB structure index to batch index
 * (entry in frame vector)
 * @param key of neighbor list in frame class "2-body" or "3-body" that should be collected
 *
 * Function will construct a functional NearestNeighborNSquare.
 * The neighbor list corresponding to keys have to be fully set up in the frames
 * vector. The frames can be set up in a batched form. For this the maps atomBatchMap
 * and structToBatch have to be supplied in a functional way. The 2D arrays
 * strucIndex and atomIdx are used to select the needed central atoms.
 *******************************************************************************************/
std::shared_ptr<NearestNeighborNSquare> generateNeighborListFromBatchListFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const String&                    key);
/*******************************************************************************************
 * Create a NeighborListMap object from vector of frame classes.
 *
 * @param frames Vector of frame classes containing set up neighbor lists and SHS2 and SHS3
 * @param strucIdx Original structure indices-1 as in ML_AB file. Format [type][structire index]
 * @param atomIdx Original atom indices-1 as in ML_AB file. Format [type][atom index]
 * @param atomBatchMap by supplying the original structure and atom index of ML_AB
 * class will give back the atom index in the corresponding type ordered frame.
 * @param strucToBatch mapping of the original ML_AB structure index to batch index
 * (entry in frame vector)
 * @param key is a vector of keys which can contain {"2-body", "3-body"}
 *
 * Class will super neighbor list from a set of frame classes
 * which have to be set up and contain ready to use neighbor lists in type sorted format,
 * By the supplied strucIdx and atomIdx in combination with the set up classes
 * atomBatchMap strucToBatch the user can chose the central atoms which are selected from
 * supplied (maybe batched) frame vector.
 *******************************************************************************************/
std::map<String, std::shared_ptr<NearestNeighborNSquare>> generateNeighborListFromBatchListFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const Vec1String&                keys = {});
/*******************************************************************************************
 * Create a DescriptorSHS2 object from vector of frame classes.
 *
 * @param frames Vector of frame classes containing set up neighbor lists and SHS2 and SHS3
 * @param strucIdx Original structure indices-1 as in ML_AB file. Format [type][structire index]
 * @param atomIdx Original atom indices-1 as in ML_AB file. Format [type][atom index]
 * @param atomBatchMap by supplying the original structure and atom index of ML_AB
 * class will give back the atom index in the corresponding type ordered frame.
 * @param strucToBatch mapping of the original ML_AB structure index to batch index
 * (entry in frame vector)
 * @param key of descriptor key "SHS3-3-body" that should be collected
 *
 * Function will construct a functional DescriptorSHS2.
 * The descriptors corresponding to keys have to be fully set up in the frames
 * vector. The frames can be set up in a batched form. For this the maps atomBatchMap
 * and structToBatch have to be supplied in a functional way. The 2D arrays
 * strucIndex and atomIdx are used to select the needed central atoms.
 *******************************************************************************************/
std::shared_ptr<DescriptorSHS2> generateDescriptorSHS2FromBatchFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const String&                    key = "SHS2-2-body");
/*******************************************************************************************
 * Create a DescriptorSHS3 object from vector of frame classes.
 *
 * @param frames Vector of frame classes containing set up neighbor lists and SHS2 and SHS3
 * @param strucIdx Original structure indices-1 as in ML_AB file. Format [type][structire index]
 * @param atomIdx Original atom indices-1 as in ML_AB file. Format [type][atom index]
 * @param atomBatchMap by supplying the original structure and atom index of ML_AB
 * class will give back the atom index in the corresponding type ordered frame.
 * @param strucToBatch mapping of the original ML_AB structure index to batch index
 * (entry in frame vector)
 * @param key of descriptor key "SHS3-3-body" that should be collected
 *
 * Function will construct a functional DescriptorSHS3.
 * The descriptors corresponding to keys have to be fully set up in the frames
 * vector. The frames can be set up in a batched form. For this the maps atomBatchMap
 * and structToBatch have to be supplied in a functional way. The 2D arrays
 * strucIndex and atomIdx are used to select the needed central atoms.
 *******************************************************************************************/
std::shared_ptr<DescriptorSHS3> generateDescriptorSHS3FromBatchFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const String&                    key = "SHS3-3-body");
/*******************************************************************************************
 * Create a DescriptorMap object from vector of frame classes.
 *
 * @param frames Vector of frame classes containing set up neighbor lists and SHS2 and SHS3
 * @param strucIdx Original structure indices-1 as in ML_AB file. Format [type][structire index]
 * @param atomIdx Original atom indices-1 as in ML_AB file. Format [type][atom index]
 * @param atomBatchMap by supplying the original structure and atom index of ML_AB
 * class will give back the atom index in the corresponding type ordered frame.
 * @param strucToBatch mapping of the original ML_AB structure index to batch index
 * (entry in frame vector)
 * @param key of descriptor from "SHS2-2-body", "SHS3-3-body" set which should be collected
 *
 * Function will construct a functional DescriptorMap for the descriptors supplied with
 * keys. The descriptors corresponding to keys have to be fully set up in the frames
 * vector. The frames can be set up in a batched form. For this the maps atomBatchMap
 * and structToBatch have to be supplied in a functional way. The 2D arrays
 * strucIndex and atomIdx are used to select the needed central atoms.
 *******************************************************************************************/
std::map<String, std::shared_ptr<Descriptor>> generateDescriptorFromBatchListFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const Vec1String&                keys = {});
/*******************************************************************************************
 * Create a DescriptorCollector object from vector of frame classes.
 *
 * @param frames Vector of frame classes containing set up neighbor lists and SHS2 and SHS3
 * @param strucIdx Original structure indices-1 as in ML_AB file. Format [type][structire index]
 * @param atomIdx Original atom indices-1 as in ML_AB file. Format [type][atom index]
 * @param atomBatchMap by supplying the original structure and atom index of ML_AB
 * class will give back the atom index in the corresponding type ordered frame.
 * @param strucToBatch mapping of the original ML_AB structure index to batch index
 * (entry in frame vector)
 *
 * Class will generate local refernce configurations from a set of frame classes
 * which have to be set up contain ready to use neighbor lists in type sorted format,
 * DescriptorSHS2 and DescriptorSHS3 classes corresponding the neighbor lists.
 * By the supplied strucIdx and atomIdx in combination with the set up classes
 * atomBatchMap strucToBatch the user can chose local reference configurations from the
 * supplied (maybe batched) frame vector.
 *******************************************************************************************/
std::shared_ptr<DescriptorCollector> createDescriptorCollectorReferenceConfigs(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch);

} //namespace batch_tools

} //namespace vaspml

#endif
