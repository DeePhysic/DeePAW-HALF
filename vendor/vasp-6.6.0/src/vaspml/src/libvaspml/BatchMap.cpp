#include "BatchMap.hpp"

#include "Structure.hpp"
#include "Tutor.hpp"

#include <cmath>          // for ceil
#include <numeric>        // for iota
#include <unordered_set>  // for unordered_set

using namespace vaspml;

Int AtomBatchMap::get_mapStrucAtom_Batch(const Int& strucIdx, const Int& atomIdx) const
{
    return mapStrucAtom_Batch[strucIdx][atomIdx];
}

Int AtomBatchMap::get_mapStrucAtom_BatchTO(const Int& strucIdx, const Int& atomIdx) const
{
    return mapStrucAtom_BatchTO[strucIdx][atomIdx];
}

Int AtomBatchMap::get_strucIndex(const Int& flatBatchIdx) const
{
    return strucIndex[flatBatchIdx];
}

Int AtomBatchMap::get_strucIndexTO(const Int& flatBatchIdx) const
{
    return strucIndexTO[flatBatchIdx];
}

Int AtomBatchMap::get_atomIndex(const Int& flatBatchIdx) const
{
    return atomIndex[flatBatchIdx];
}

const Vec1Int& AtomBatchMap::get_atomIndex() const
{
    return atomIndex;
}

Int AtomBatchMap::get_atomIndexTO(const Int& flatBatchIdx) const
{
    return atomIndexTO[flatBatchIdx];
}

const Vec1Int& AtomBatchMap::get_atomIndexTO() const
{
    return atomIndexTO;
}

Int AtomBatchMap::get_mapTO_OrigOrder(const Int& flatBatchIdx) const
{
    return mapTO_OrigOrder[flatBatchIdx];
}

Int AtomBatchMap::get_mapOrigOrder_TO(const Int& flatBatchIdx) const
{
    return mapOrigOrder_TO[flatBatchIdx];
}

Int AtomBatchMap::get_mapTO(const Int& inputIdx) const
{
    return mapTO[inputIdx];
}

Int AtomBatchMap::get_typesIdx(const Int& atomIndx) const
{
    return typesIdx[atomIndx];
}

String AtomBatchMap::get_types(const Int& atomIndx) const
{
    return types[atomIndx];
}

Int AtomBatchMap::get_typesIdxTO(const Int& atomIndx) const
{
    return typesIdxTO[atomIndx];
}

String AtomBatchMap::get_typesTO(const Int& atomIndx) const
{
    return typesTO[atomIndx];
}

const Vec1Int& AtomBatchMap::get_typesIdx() const
{
    return typesIdx;
}

const Vec1String& AtomBatchMap::get_types() const
{
    return types;
}

const Vec1Int& AtomBatchMap::get_typesIdxTO() const
{
    return typesIdxTO;
}

const Vec1String& AtomBatchMap::get_typesTO() const
{
    return typesTO;
}

const Vec1String& AtomBatchMap::get_typeOrder() const
{
    return typeOrder;
}

Int AtomBatchMap::get_numberAtomsPerStruc(const Int nStruc) const
{
    return numberAtomsPerStruc[nStruc];
}

const Vec1Int& AtomBatchMap::get_numberAtomsPerStruc() const
{
    return numberAtomsPerStruc;
}

void AtomBatchMap::makeMap(const std::vector<Structure>& structures, const Vec1String& typeOrder)
{
    nStruc = structures.size();
    set_numberAtomsPerStruc(structures);
    makeFlatMaps(structures);
    makeUniqueTypeOrder(typeOrder);
    makeIntTypes();
    makeTypeOrderMaps();
}

void AtomBatchMap::set_numberAtomsPerStruc(const std::vector<Structure>& structures)
{
    numberAtomsPerStruc.resize(structures.size());
    for (Size i = 0; i < structures.size(); i++)
    {
        numberAtomsPerStruc[i] = structures[i].get_numAtoms();
    }
}

void AtomBatchMap::makeFlatMaps(const std::vector<Structure>& structures)
{
    Size totalAtoms = 0;
    std::for_each(structures.begin(),
                  structures.end(),
                  [&totalAtoms](const auto& struc) { totalAtoms += struc.get_numAtoms(); });
    types.resize(totalAtoms);
    strucIndex.resize(totalAtoms);
    atomIndex.resize(totalAtoms);
    atomIndexTO.resize(totalAtoms);
    mapStrucAtom_Batch.resize(structures.size());
    mapStrucAtom_BatchTO.resize(structures.size());
    nAtomsInStruc.resize(structures.size());
    Size counter = 0;
    for (Size nStruc = 0; nStruc < structures.size(); nStruc++)
    {
        const auto& keys = structures[nStruc].get_atomTypeNames();
        mapStrucAtom_Batch[nStruc].resize(keys.size());
        mapStrucAtom_BatchTO[nStruc].resize(keys.size());
        nAtomsInStruc[nStruc] = keys.size();
        for (Size atom = 0; atom < keys.size(); atom++)
        {
            types[counter] = keys[atom];
            strucIndex[counter] = nStruc;
            mapStrucAtom_Batch[nStruc][atom] = counter;
            atomIndex[counter] = atom;
            counter++;
        }
    }
}

void AtomBatchMap::makeTypeOrderMaps()
{
    mapTO.resize(types.size());
    std::iota(mapTO.begin(), mapTO.end(), 0);
    std::stable_sort(mapTO.begin(),
                     mapTO.end(),
                     [&](const Size& idx1, const Size& idx2)
                     {
                         const String& el1 = types[idx1];
                         const String& el2 = types[idx2];
                         Int           r1 =
                             std::distance(typeOrder.cbegin(),
                                           std::find(typeOrder.cbegin(), typeOrder.cend(), el1));
                         Int r2 =
                             std::distance(typeOrder.cbegin(),
                                           std::find(typeOrder.cbegin(), typeOrder.cend(), el2));
                         return r1 < r2;
                     });
    // Build a map from string to its position in custom_order
    strucIndexTO.resize(types.size());
    mapTO_OrigOrder.resize(types.size());
    mapOrigOrder_TO.resize(types.size());
    typesTO.resize(types.size());
    for (Size i = 0; i < types.size(); i++)
    {
        strucIndexTO[i] = strucIndex[mapTO[i]];
        mapTO_OrigOrder[i] = mapTO[i];
        mapOrigOrder_TO[mapTO[i]] = i;
        typesTO[i] = types[mapTO[i]];
        typesIdxTO[i] = typesIdx[mapTO[i]];
        atomIndexTO[i] = atomIndex[mapTO[i]];
        mapStrucAtom_BatchTO[strucIndex[mapTO[i]]][atomIndexTO[i]] = i;
    }
}

void AtomBatchMap::makeUniqueTypeOrder(const Vec1String& typeOrder)
{
    if (typeOrder.empty())
    {
        std::unordered_set<String> seen;
        for (String value : types)
        {
            if (seen.find(value) == seen.end())
            {
                seen.insert(value);
                this->typeOrder.push_back(value);
            }
        }
    }
    else this->typeOrder = typeOrder;
}

void AtomBatchMap::makeIntTypes()
{
    typesIdx.resize(types.size());
    typesIdxTO.resize(types.size());
    Int currentIdx = 0;
    std::transform(types.cbegin(),
                   types.cend(),
                   typesIdx.begin(),
                   [&](const String& s)
                   {
                       auto it = strTypeToIntType.find(s);
                       if (it == strTypeToIntType.end())
                       {
                           int id = currentIdx++;
                           strTypeToIntType[s] = id;
                           return id;
                       }
                       else { return it->second; }
                   });
}

BatchManager::BatchManager(const Size& batchSize) : batchSize(batchSize)
{
    if (batchSize == 0)
    {
        global::tutor.bug("Error: BatchManager::BatchManager( const Size& batchSize )\n."
                          "batchSize can not be zero!");
    }
}

Size BatchManager::getBatchIndexForStructure(const Size& structureIndex) const
{
    if (structureIndex >= numStructures)
    {
        global::tutor.bug(
            "Error: BatchManager::getBatchIndexForStructure( const Size& structureIndex )\n."
            "batchSize can not be zero!");
    }
    return structureIndex / batchSize;
}

std::vector<Size> BatchManager::getStructureIndicesInBatch(const Size& batchIndex) const
{
    Size              startIndex = batchIndex * batchSize;
    Size              endIndex = std::min(startIndex + batchSize, numStructures);
    std::vector<Size> structureIndices;
    for (Size i = startIndex; i < endIndex; ++i) { structureIndices.push_back(i); }
    return structureIndices;
}

Size BatchManager::getStructureIndexWithinBatch(const Size& structureIndex) const
{
    if (structureIndex >= numStructures)
    {
        global::tutor.bug(
            "Error: BatchManager::getStructureIndexWithinBatch( const Size& structureIndex )\n."
            "structureIndex is out of bounds.");
    }
    return structureIndex % batchSize;
}

Size BatchManager::getNumberOfBatches() const
{
    return static_cast<Size>(std::ceil(static_cast<Real>(numStructures) / batchSize));
}

void BatchManager::setNumberOfStructures(const Size& numStructures)
{
    this->numStructures = numStructures;
}
