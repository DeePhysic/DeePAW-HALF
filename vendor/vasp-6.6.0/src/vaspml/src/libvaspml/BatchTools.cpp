#include "BatchTools.hpp"

#include "BatchMap.hpp"
#include "Descriptor.hpp"
#include "DescriptorCollector.hpp"
#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"
#include "Frame.hpp"
#include "constants.hpp"

using namespace vaspml;

std::vector<std::shared_ptr<NearestNeighborNSquare>>
vaspml::batch_tools::makeNeighborListVectorFromFrame(const std::vector<Frame>& frames,
                                                     const String              key)
{
    return makeNeighborListVectorFromFrame(frames.begin(), frames.end(), key);
}

std::shared_ptr<NearestNeighborNSquare>
vaspml::batch_tools::generateNeighborListFromBatchListFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const String&                    key)
{
    std::shared_ptr<NearestNeighborNSquare> nList = std::make_shared<NearestNeighborNSquare>();
    for (Size type = 0; type < strucIdx.size(); type++)
    {
        for (Size locRef = 0; locRef < atomIdx[type].size(); locRef++)
        {
            const Size& batch = strucToBatch.getBatchIndexForStructure(strucIdx[type][locRef]);
            const Size& strucInBatch =
                strucToBatch.getStructureIndexWithinBatch(strucIdx[type][locRef]);
            const Size& atomInBatchTO =
                atomBatchMap[batch].get_mapStrucAtom_BatchTO(strucInBatch, atomIdx[type][locRef]);
            neighbor_list::addElementData(*frames[batch].get_neighborLists().at(key),
                                          *nList,
                                          atomInBatchTO);
        }
    }
    nList->computeNumberAtomsPerType(atomBatchMap[0].get_typeOrder());
    nList->compute_centralAtomIndexPerType();

    return nList;
}

std::map<String, std::shared_ptr<NearestNeighborNSquare>>
vaspml::batch_tools::generateNeighborListFromBatchListFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const Vec1String&                keys)
{
    Vec1String tempKeys;
    if (keys.empty()) { tempKeys = constants::bodyKeyList; }
    else { tempKeys = keys; }

    std::map<String, std::shared_ptr<NearestNeighborNSquare>> nList;
    for (const auto& key : tempKeys)
    {
        nList[key] = generateNeighborListFromBatchListFrames(frames,
                                                             strucIdx,
                                                             atomIdx,
                                                             atomBatchMap,
                                                             strucToBatch,
                                                             key);
        neighbor_list::copyParameters(*frames[0].get_neighborLists().at(key), *nList[key]);
    }

    return nList;
}

std::shared_ptr<DescriptorSHS2> vaspml::batch_tools::generateDescriptorSHS2FromBatchFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const String&                    key)
{
    const DescriptorSHS2* desc =
        static_cast<const DescriptorSHS2*>(frames[0].get_descriptors().at(key).get());
    const Real                      weight = desc->get_weight();
    std::shared_ptr<DescriptorSHS2> descriptorNew = std::make_shared<DescriptorSHS2>();
    // if descriptor not used, return empty descriptor
    if (weight <= 0) return descriptorNew;
    for (Size type = 0; type < strucIdx.size(); type++)
    {
        for (Size locRef = 0; locRef < atomIdx[type].size(); locRef++)
        {
            const Size& batch = strucToBatch.getBatchIndexForStructure(strucIdx[type][locRef]);
            const Size& strucInBatch =
                strucToBatch.getStructureIndexWithinBatch(strucIdx[type][locRef]);
            const Size& atomInBatchTO =
                atomBatchMap[batch].get_mapStrucAtom_BatchTO(strucInBatch, atomIdx[type][locRef]);
            const DescriptorSHS2* desc =
                static_cast<const DescriptorSHS2*>(frames[batch].get_descriptors().at(key).get());
            addElementData(*desc, *descriptorNew, atomInBatchTO);
        }
    }

    return descriptorNew;
}

std::shared_ptr<DescriptorSHS3> vaspml::batch_tools::generateDescriptorSHS3FromBatchFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const String&                    key)
{
    const DescriptorSHS3* desc =
        static_cast<const DescriptorSHS3*>(frames[0].get_descriptors().at(key).get());
    const Real                      weight = desc->get_weight();
    std::shared_ptr<DescriptorSHS3> descriptorNew = std::make_shared<DescriptorSHS3>();
    // if descriptor not used, return empty descriptor
    if (weight <= 0) return descriptorNew;
    for (Size type = 0; type < strucIdx.size(); type++)
    {
        for (Size locRef = 0; locRef < atomIdx[type].size(); locRef++)
        {
            const Size& batch = strucToBatch.getBatchIndexForStructure(strucIdx[type][locRef]);
            const Size& strucInBatch =
                strucToBatch.getStructureIndexWithinBatch(strucIdx[type][locRef]);
            const Size& atomInBatchTO =
                atomBatchMap[batch].get_mapStrucAtom_BatchTO(strucInBatch, atomIdx[type][locRef]);
            const DescriptorSHS3* desc =
                static_cast<const DescriptorSHS3*>(frames[batch].get_descriptors().at(key).get());
            addElementData(*desc, *descriptorNew, atomInBatchTO);
        }
    }

    return descriptorNew;
}

std::map<String, std::shared_ptr<Descriptor>>
vaspml::batch_tools::generateDescriptorFromBatchListFrames(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch,
    const Vec1String&                keys)
{
    Vec1String tempKeys;
    if (keys.empty()) { tempKeys = constants::descriptorKeyList; }
    else { tempKeys = keys; }

    std::map<String, std::shared_ptr<Descriptor>> descriptorMap;
    for (const auto& key : tempKeys)
    {
        if (key == "SHS2-2-body")
        {
            descriptorMap[key] = generateDescriptorSHS2FromBatchFrames(frames,
                                                                       strucIdx,
                                                                       atomIdx,
                                                                       atomBatchMap,
                                                                       strucToBatch,
                                                                       key);
            const DescriptorSHS2* desc1 =
                static_cast<const DescriptorSHS2*>(frames[0].get_descriptors().at(key).get());
            DescriptorSHS2* desc2 = static_cast<DescriptorSHS2*>(descriptorMap[key].get());
            copyParameters(*desc1, *desc2);
        }
        else if (key == "SHS3-3-body")
        {
            descriptorMap[key] = generateDescriptorSHS3FromBatchFrames(frames,
                                                                       strucIdx,
                                                                       atomIdx,
                                                                       atomBatchMap,
                                                                       strucToBatch,
                                                                       key);
            const DescriptorSHS3* desc1 =
                static_cast<const DescriptorSHS3*>(frames[0].get_descriptors().at(key).get());
            DescriptorSHS3* desc2 = static_cast<DescriptorSHS3*>(descriptorMap[key].get());
            copyParameters(*desc1, *desc2);
        }
    }

    return descriptorMap;
}

std::shared_ptr<DescriptorCollector> vaspml::batch_tools::createDescriptorCollectorReferenceConfigs(
    const std::vector<Frame>&        frames,
    const Vec2Int&                   strucIdx,
    const Vec2Int&                   atomIdx,
    const std::vector<AtomBatchMap>& atomBatchMap,
    const BatchManager&              strucToBatch)
{
    auto descriptors = generateDescriptorFromBatchListFrames(frames,
                                                             strucIdx,
                                                             atomIdx,
                                                             atomBatchMap,
                                                             strucToBatch);
    // filter for needed neighbor lists.
    Vec1String neighborKeys;
    for (const auto& [key, item] : descriptors)
    {
        if (item->get_weight() > 0) neighborKeys.push_back(key.substr(5));
    }
    auto neighborLists = generateNeighborListFromBatchListFrames(frames,
                                                                 strucIdx,
                                                                 atomIdx,
                                                                 atomBatchMap,
                                                                 strucToBatch,
                                                                 neighborKeys);
    for (auto& [key, list] : neighborLists)
    {
        const String& descKey = constants::bodyDescMap.at(key);
        if (descriptors[descKey]->get_weight() > 0)
        {
            descriptors[descKey]->set_neighborList(list);
        }
    }

    std::shared_ptr<DescriptorCollector> collector =
        std::make_shared<DescriptorCollector>(0, std::map<String, ShVec2Real>(), nullptr);
    collector->setDescriptorMap(descriptors);
    collector->updateCollector();

    return collector;
}
