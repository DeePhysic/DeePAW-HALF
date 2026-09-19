#include "MlMPI.hpp"

#include "dataset.hpp"

#include "Item.hpp"
#include "Record.hpp"
#include "Structure.hpp"
#include "Tutor.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <algorithm> // for for_each, transform, find, max, all_of, rem...
#include <iostream>  // for basic_ostream, cout
#include <iterator>  // for distance

using namespace vaspml;

void dataset::getPerStructureLrc(const Record& set, Vec1Int& lrcStructure, Vec1Int& lrcAtom)
{
    lrcStructure.clear();
    lrcAtom.clear();

    Vec1ShCRec structures = set.vcget<Vec1ShRec>("structures");
    for (Size i = 0; i < structures.size(); ++i)
    {
        const Record& s = *structures[i];
        if (s.contains("lrc"))
        {
            const Vec1Int& lrc = s.cget<Vec1Int>("lrc");
            lrcStructure.resize(lrc.size() + lrcStructure.size(), i);
            lrcAtom.insert(lrcAtom.end(), lrc.begin(), lrc.end());
        }
    }

    return;
}

void dataset::setPerStructureLrc(const Vec2Int& lrcStructure, const Vec2Int& lrcAtom, Record& set)
{
    Vec1ShRec& structures = set.get<Vec1ShRec>("structures");

    // First, clear all per-structure lrc vectors.
    for (Size i = 0; i < structures.size(); ++i)
    {
        Record& s = *structures[i];
        if (s.contains("lrc")) s.get<Vec1Int>("lrc").clear();
        else s.add("lrc", ItemIndex::VEC1INT);
    }

    // Add individual lrcs for each structure,atom-pair.
    for (Size t = 0; t < lrcStructure.size(); ++t)
    {
        for (Size i = 0; i < lrcStructure[t].size(); ++i)
        {
            structures[lrcStructure[t][i]]->get<Vec1Int>("lrc").push_back(lrcAtom[t][i]);
        }
    }

    return;
}

void dataset::concat(Record& set1, const Record& set2)
{
    set1.get<Int>("numStructures") += set2.cget<Int>("numStructures");
    set1.get<Int>("maxTypes") = std::max(set1.cget<Int>("maxTypes"), set2.cget<Int>("maxTypes"));

    // Determine which types must be added from set2.
    Vec1String&       types1 = set1.get<Vec1String>("typeNames");
    const Vec1String& types2 = set2.cget<Vec1String>("typeNames");
    for (Size i = 0; i < types2.size(); ++i)
    {
        Size indexInType1 =
            std::distance(types1.begin(), std::find(types1.begin(), types1.end(), types2[i]));
        if (indexInType1 >= types1.size()) types1.push_back(types2[i]);
    }

    set1.get<Int>("maxAtoms") = std::max(set1.cget<Int>("maxAtoms"), set2.cget<Int>("maxAtoms"));

    for (Size i = 0; i < types1.size(); ++i)
    {
        std::cout << str("%zu %s\n", i, types1.at(i).c_str());
    }

    return;
}

std::vector<Int> dataset::makeTypeMap(const Vec1String& subTypes, const Vec1String& types)
{
    VASPML_DEBUG_L1(
        if (std::all_of(subTypes.cbegin(),
                        subTypes.cend(),
                        [&types](const String& t)
                        { return std::find(types.cbegin(), types.cend(), t) == types.end(); }))
        {
            global::tutor.error(
                "Structure types not found in ML_AB types. Your ML_AB files is wrong.\n");
        }
    );
    Vec1String typesInStruc = types;
    typesInStruc.erase(std::remove_if(typesInStruc.begin(),
                                      typesInStruc.end(),
                                      [&](const String& t)
                                      {
                                          return std::find(subTypes.cbegin(), subTypes.cend(), t)
                                              == subTypes.cend();
                                      }),
                       typesInStruc.end());
    Vec1Int indexMap(subTypes.size());
    std::transform(subTypes.cbegin(),
                   subTypes.cend(),
                   indexMap.begin(),
                   [&](const String& t)
                   {
                       return std::distance(
                           typesInStruc.cbegin(),
                           std::find(typesInStruc.cbegin(), typesInStruc.cend(), t));
                   });

    return indexMap;
}

Vec2Int dataset::makeTypeMapSet(Record& set)
{
    const Vec1String& mlabTypes = set.cget<Vec1String>("typeNames");
    Vec1ShRec&        MLABstrucs = set.get<Vec1ShRec>("structures");
    Vec2Int           typeMaps;
    // generate type maps
    std::for_each(MLABstrucs.cbegin(),
                  MLABstrucs.cend(),
                  [&mlabTypes, &typeMaps](const ShRec& struc)
                  {
                      const Vec1String& types = struc->cget<Vec1String>("typeNames");
                      typeMaps.push_back(makeTypeMap(types, mlabTypes));
                  });

    return typeMaps;
}

Vec2Int dataset::createReorderMap(const Record& set, const Vec2Int& typeMaps)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");
    Vec2Int    atomOrder(MLABstrucs.size());
    for (Size nStruc = 0; nStruc < MLABstrucs.size(); nStruc++)
    {
        const Int& nAtoms = MLABstrucs[nStruc]->cget<Int>("numAtoms");
        atomOrder[nStruc].resize(nAtoms);
        const Size&    nTypes = typeMaps[nStruc].size();
        const Vec1Int& numTypesPerAtom = MLABstrucs[nStruc]->cget<Vec1Int>("numAtomsPerType");
        Vec1Int        numTypesPerAtomReord(nTypes);
        // temporarily reorder numAtomsPerType for computing type Offset in new ordering
        std::transform(typeMaps[nStruc].cbegin(),
                       typeMaps[nStruc].cend(),
                       numTypesPerAtomReord.begin(),
                       [&numTypesPerAtom](const Int& indx) { return numTypesPerAtom[indx]; });
        Size atom = 0;
        for (Size itype = 0; itype < nTypes; itype++)
        {
            const Int& newType = typeMaps[nStruc][itype];
            //Int typeOffset = std::reduce( VASPML_POLICY_SEQ
            //                              numTypesPerAtomReord.cbegin(),
            //                              numTypesPerAtomReord.cbegin()+newType,0);
            Int typeOffset = 0;
            for (Size i = 0; i < (Size)newType; ++i) typeOffset += numTypesPerAtomReord[i];
            for (Int atomInType = 0; atomInType < numTypesPerAtom[itype]; atomInType++)
            {
                atomOrder[nStruc][atom] = typeOffset + atomInType;
                atom++;
            }
        }
    }

    return atomOrder;
}

void dataset::reorderSingleStructure(Record&        structure,
                                     const Vec1Int& atomMap,
                                     const Vec1Int& typeMap)
{
    // reorder type names
    Vec1String&      typeNames = structure.get<Vec1String>("typeNames");
    const Vec1String typeNamesCopy = structure.get<Vec1String>("typeNames");
    std::transform(typeMap.cbegin(),
                   typeMap.cend(),
                   typeNames.begin(),
                   [&typeNamesCopy](const Int& indx) { return typeNamesCopy[indx]; });
    // reorder number atoms per type
    Vec1Int&      numAtomsPerType = structure.get<Vec1Int>("numAtomsPerType");
    const Vec1Int numAtomsPerTypeCopy = numAtomsPerType;
    std::transform(typeMap.cbegin(),
                   typeMap.cend(),
                   numAtomsPerType.begin(),
                   [&numAtomsPerTypeCopy](const Int& indx) { return numAtomsPerTypeCopy[indx]; });
    // reorder positions
    Vec1Real&      positions = structure.get<Vec1Real>("positions");
    const Vec1Real positionsCopy = positions;
    std::for_each(atomMap.cbegin(),
                  atomMap.cend(),
                  [&positions, &positionsCopy](const Int& indx)
                  {
                      positions[3 * indx] = positionsCopy[3 * indx];
                      positions[3 * indx + 1] = positionsCopy[3 * indx + 1];
                      positions[3 * indx + 2] = positionsCopy[3 * indx + 2];
                  });

    Vec1Real&      forces = structure.get<Vec1Real>("forces");
    const Vec1Real forcesCopy = forces;
    std::for_each(atomMap.cbegin(),
                  atomMap.cend(),
                  [&forces, &forcesCopy](const Int& indx)
                  {
                      forces[3 * indx] = forcesCopy[3 * indx];
                      forces[3 * indx + 1] = forcesCopy[3 * indx + 1];
                      forces[3 * indx + 2] = forcesCopy[3 * indx + 2];
                  });
}

void dataset::reorderStructures(Vec1ShRec&     structures,
                                const Vec2Int& atomMap,
                                const Vec2Int& typeMap)
{
    for (Size i = 0; i < structures.size(); i++)
    {
        reorderSingleStructure(*structures[i], atomMap[i], typeMap[i]);
    }
}

void dataset::reorderLocalReferenceConfigurations(Record& set, const Vec2Int& atomMap)
{
    Vec2Int& lrcStructure = set.get<Vec2Int>("lrcStructure");
    Vec2Int& lrcAtom = set.get<Vec2Int>("lrcAtom");
    for (Size type = 0; type < lrcStructure.size(); type++)
    {
        for (Size nlrc = 0; nlrc < lrcAtom[type].size(); nlrc++)
        {
            // in c++ arrays start at zero !!
            // reordering part.
            const Int& nstruc = lrcStructure[type][nlrc];
            const Int& newAtom = atomMap[nstruc][lrcAtom[type][nlrc]];
            lrcAtom[type][nlrc] = newAtom;
        }
    }
}

void dataset::reorderRecord(Record& set, const Vec2Int& atomMap, const Vec2Int& typeMap)
{
    Vec1ShRec structures = set.get<Vec1ShRec>("structures");
    reorderStructures(structures, atomMap, typeMap);
    reorderLocalReferenceConfigurations(set, atomMap);
}

void dataset::reorderSet(Record& set)
{
    Vec2Int typeMaps = makeTypeMapSet(set);
    Vec2Int atomMaps = createReorderMap(set, typeMaps);
    reorderRecord(set, atomMaps, typeMaps);
}

std::vector<Structure> dataset::prepareStructures(Record& mlab)
{
    Vec1ShRec MLABstrucs = mlab.get<Vec1ShRec>("structures");
    //std::vector<Structure> strucs( MLABstrucs.size() );
    std::vector<Structure> strucs;
    const Vec1String&      mlabTypes = mlab.cget<Vec1String>("typeNames");
    reorderSet(mlab);
    // prepare structure classes
    //std::transform( MLABstrucs.begin(), MLABstrucs.end(), strucs.begin(),
    //        [&mlabTypes]( const ShRec& inputStruc ){
    //           Structure struc = Structure( inputStruc, false );
    //           struc.set_typeMap( mlabTypes );
    //           return struc;
    //        });
    for (Size i = 0; i < MLABstrucs.size(); i++)
    {
        strucs.emplace_back(Structure(MLABstrucs[i], false));
        strucs[i].set_typeMap(mlabTypes);
    }

    return strucs;
}

std::vector<Structure> dataset::prepareStructures(Record& mlab, const std::shared_ptr<MlMPI>& mlmpi)
{
    std::vector<Structure> strucs = dataset::prepareStructures(mlab);
    for (Structure& struc : strucs)
    {
        struc.set_atomIndx(mlmpi::getRoundRobinIndexes(struc.get_numAtoms(),
                                                       mlmpi->get_numberRanks(),
                                                       mlmpi->get_rank()));
    }

    return strucs;
}

Vec1Real dataset::extractEnergies(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractEnergies(MLABstrucs);
}

Vec1Real dataset::extractEnergies(const Vec1ShCRec& MLABstrucs)
{
    Vec1Real energies;
    energies.reserve(MLABstrucs.size());
    for (const auto& struc : MLABstrucs) { energies.push_back(struc->cget<Real>("energy")); }

    return energies;
}

Vec1Int dataset::extractNumberAtoms(const Record& set)
{
    const Vec1ShCRec& MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractNumberAtoms(MLABstrucs);
}

Vec1Int dataset::extractNumberAtoms(const Vec1ShCRec& MLABstrucs)
{
    Vec1Int numAtoms(MLABstrucs.size());
    std::transform(MLABstrucs.cbegin(),
                   MLABstrucs.cend(),
                   numAtoms.begin(),
                   [](const std::shared_ptr<const Record>& struc)
                   { return struc->cget<Int>("numAtoms"); });

    return numAtoms;
}

Vec2Int dataset::extractNumberAtomsPerType(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractNumberAtomsPerType(MLABstrucs);
}

Vec2Int dataset::extractNumberAtomsPerType(const Vec1ShCRec& MLABstrucs)
{
    Vec2Int numAtoms(MLABstrucs.size());
    std::transform(MLABstrucs.cbegin(),
                   MLABstrucs.cend(),
                   numAtoms.begin(),
                   [](const std::shared_ptr<const Record>& struc)
                   { return struc->cget<Vec1Int>("numAtomsPerType"); });

    return numAtoms;
}

Vec1Int dataset::countTotalAtomsPerType(const Record& set)
{
    const Vec1String& superTypes = set.cget<Vec1String>("typeNames");
    Vec1Int           totalAtomsPerType(superTypes.size());
    const Vec1ShCRec& MLABstrucs = set.vcget<Vec1ShRec>("structures");
    std::for_each(MLABstrucs.cbegin(),
                  MLABstrucs.cend(),
                  [&](const ShCRec& struc)
                  {
                      const Vec1String& subTypes = struc->cget<Vec1String>("typeNames");
                      const Vec1Int&    nAtomsStruc = struc->cget<Vec1Int>("numAtomsPerType");

                      TypeMap typeMap(superTypes, subTypes);

                      for (Size ffType = 0; ffType < superTypes.size(); ffType++)
                      {
                          Int type = typeMap.toSubType(ffType);
                          if (type >= 0) totalAtomsPerType[ffType] += nAtomsStruc[type];
                      }
                  });

    return totalAtomsPerType;
}

Vec2Real dataset::extractForceComponent(const Record& set, const Int xyzIndex)
{

    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractForceComponent(MLABstrucs, xyzIndex);
}

Vec2Real dataset::extractForceComponent(const Vec1ShCRec& MLABstrucs, const Int xyzIndex)
{
    VASPML_DEBUG_L1(
        if (xyzIndex < 0 || xyzIndex > 2)
        {
            global::tutor.error("ERROR: Vec2Real dataset::extractForceComponent( const Record& set,"
                                "const Int xyzIndex )\n xyzIndex must be 0, 1 or 2.\n");
        }
    );
    Vec2Real forces(MLABstrucs.size());
    std::transform(MLABstrucs.cbegin(),
                   MLABstrucs.cend(),
                   forces.begin(),
                   [&xyzIndex](const std::shared_ptr<const Record>& struc)
                   {
                       const Vec1Real& forceStruc = struc->cget<Vec1Real>("forces");
                       Vec1Real        forceX(struc->cget<Int>("numAtoms"));
                       for (Int i = 0; i < struc->cget<Int>("numAtoms"); i++)
                       {
                           forceX[i] = forceStruc[3 * i + xyzIndex];
                       }
                       return forceX;
                   });

    return forces;
}

Vec2Real dataset::extractForcesX(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractForcesX(MLABstrucs);
}

Vec2Real dataset::extractForcesX(const Vec1ShCRec& MLABstrucs)
{
    return extractForceComponent(MLABstrucs, 0);
}

Vec2Real dataset::extractForcesY(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractForcesY(MLABstrucs);
}

Vec2Real dataset::extractForcesY(const Vec1ShCRec& MLABstrucs)
{
    return extractForceComponent(MLABstrucs, 1);
}

Vec2Real dataset::extractForcesZ(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractForcesZ(MLABstrucs);
}

Vec2Real dataset::extractForcesZ(const Vec1ShCRec& MLABstrucs)
{
    return extractForceComponent(MLABstrucs, 2);
}

Vec2Real dataset::extractForces(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractForces(MLABstrucs);
}

Vec2Real dataset::extractForces(const Vec1ShCRec& MLABstrucs)
{
    Vec2Real forces(MLABstrucs.size());
    std::transform(MLABstrucs.cbegin(),
                   MLABstrucs.cend(),
                   forces.begin(),
                   [](const std::shared_ptr<const Record>& struc)
                   { return struc->cget<Vec1Real>("forces"); });

    return forces;
}

Vec2Real dataset::extractStress(const Record& set)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractStress(MLABstrucs);
}

Vec2Real dataset::extractStress(const Vec1ShCRec& MLABstrucs)
{
    Vec2Real stress(MLABstrucs.size());
    std::transform(MLABstrucs.cbegin(),
                   MLABstrucs.cend(),
                   stress.begin(),
                   [](const std::shared_ptr<const Record>& struc)
                   { return struc->cget<Vec1Real>("stress"); });

    return stress;
}

Vec1Real dataset::extractStressComponent(const Record& set, const Int xyzIndex)
{
    Vec1ShCRec MLABstrucs = set.vcget<Vec1ShRec>("structures");

    return extractStressComponent(MLABstrucs, xyzIndex);
}

Vec1Real dataset::extractStressComponent(const Vec1ShCRec& MLABstrucs, const Int xyzIndex)
{
    VASPML_DEBUG_L1(
        if (xyzIndex < 0 || xyzIndex > 5)
        {
            global::tutor.error(
                "ERROR: Vec1Real dataset::extractStressComponent( const Record& set,"
                "const Int xyzIndex ).\n xyzIndex must be 0, 1, 2, 3, 4 or 5.\n");
        }
    );
    Vec1Real stress(MLABstrucs.size());
    std::transform(MLABstrucs.cbegin(),
                   MLABstrucs.cend(),
                   stress.begin(),
                   [&xyzIndex](const std::shared_ptr<const Record>& struc)
                   { return struc->cget<Vec1Real>("stress")[xyzIndex]; });

    return stress;
}

ShRec dataset::createMlffFromSetup(const Record& setup)
{

    ShRec   mlffPtr = std::make_shared<Record>();
    Record& mlffRecord = *mlffPtr;
    createMlffFromSetup(setup, mlffRecord);

    return mlffPtr;
}

void dataset::createMlffFromSetup(const Record& setup, Record& mlffRecord)
{
    // create header
    mlffRecord.add("header", "ShRec");
    Record& header = mlffRecord.dget<ShRec>("header");
    std::for_each(constants::incarToMlffHeader.cbegin(),
                  constants::incarToMlffHeader.cend(),
                  [&header, &setup](const auto& key)
                  {
                      if (setup.contains(key)) { header.copyFrom(setup, key); }
                      else
                      {
                          global::tutor.error("Key \"" + key
                                              + "\" not found in setup (INCAR) Record.");
                      }
                  });

    std::for_each(constants::incarToMlff.cbegin(),
                  constants::incarToMlff.cend(),
                  [&mlffRecord, &setup](const auto& key)
                  {
                      if (setup.contains(key)) { mlffRecord.copyFrom(setup, key); }
                      else
                      {
                          global::tutor.error("Key \"" + key
                                              + "\" not found in setup (INCAR) Record.");
                      }
                  });

    ShRec            unusedTags = createUnusedMlffTags();
    const Vec1String unusedKeys = unusedTags->keys();
    std::for_each(unusedKeys.cbegin(),
                  unusedKeys.cend(),
                  [&mlffRecord, &unusedTags](const auto& key)
                  { mlffRecord.copyFrom(*unusedTags, key); });

    // adding ML_NHYP and ML_NHYP2. Only supervector exists in c++ code
    mlffRecord.put<Int>("ML_NHYP", setup.cget<Int>("ML_NHYP"));
    mlffRecord.put<Int>("ML_NHYP2", setup.cget<Int>("ML_NHYP"));
    mlffRecord.put<bool>("WRITECMAT", false);
    header.put<String>("_version", "0.2.1");
}

ShRec dataset::createUnusedMlffTags()
{
    ShRec record = std::make_shared<Record>();
    record->put<Real>("RATIO_DESC_LOW_RES", 0.0);

    record->put<bool>("ML_LWINDOW1", false);
    record->put<Int>("ML_IWINDOW1", 0);
    record->put<bool>("ML_LWINDOW2", false);
    record->put<Int>("ML_IWINDOW2", 0);

    record->put<bool>("ML_LMETRIC1", false);
    record->put<Int>("ML_NMETRIC1", 0);
    record->put<Real>("ML_RMETRIC1", 0.0);
    record->put<bool>("ML_LMETRIC2", false);
    record->put<Int>("ML_NMETRIC2", 0);
    record->put<Real>("ML_RMETRIC2", 0.0);

    record->put<bool>("ML_LVARTRAN1", false);
    record->put<Int>("ML_NVARTRAN1", 0);
    record->put<bool>("ML_LVARTRAN2", false);
    record->put<Int>("ML_NVARTRAN2", 0);
    record->put<Int>("ML_LMAX1", 0);
    record->put<Int>("ML_ISCALE_TOTEN", 2);
    return record;
}

void dataset::addMlabToMlff(const Record& mlab, Record& mlff)
{
    mlff.put<Int>("numTypes", mlab.cget<Int>("maxTypes"));
    mlff.put<Vec1String>("types", mlab.cget<Vec1String>("typeNames"));
    mlff.put<Int>("training_structures", mlab.cget<Int>("numStructures"));
    mlff.put<Vec1Int>("local_reference_cfgs", mlab.cget<Vec1Int>("numLrc"));
    mlff.put<Vec1Real>("atomicMass", mlab.cget<Vec1Real>("atomicMass"));
    mlff.put<Vec1Int>("numLrc", mlab.cget<Vec1Int>("numLrc"));
}

void dataset::mlabStatistics(Record& mlabStat, const Record& mlab)
{
    const Vec1Real energies = dataset::extractEnergies(mlab);
    Real           meanE = math::stats::calculateMean(energies);
    mlabStat.put<Real>("average-energy-per-atom", meanE);
    mlabStat.put<Real>("variance-training-energies",
                       math::stats::calculateVariance(energies, meanE));

    const Real meanX = math::stats::calculateMean(vector_tools::Vec2ToVec1(extractForcesX(mlab)));
    const Real varX =
        math::stats::calculateVariance(vector_tools::Vec2ToVec1(extractForcesX(mlab)), meanX);
    const Real meanY = math::stats::calculateMean(vector_tools::Vec2ToVec1(extractForcesX(mlab)));
    const Real varY =
        math::stats::calculateVariance(vector_tools::Vec2ToVec1(extractForcesX(mlab)), meanY);
    const Real meanZ = math::stats::calculateMean(vector_tools::Vec2ToVec1(extractForcesX(mlab)));
    const Real varZ =
        math::stats::calculateVariance(vector_tools::Vec2ToVec1(extractForcesX(mlab)), meanZ);

    mlabStat.put<Vec1Real>("variance-training-force", {varX, varY, varZ});
    Vec1Real stressVar(6, 0.0);
    for (Int i = 0; i < 6; i++)
    {
        const Real meanStress =
            math::stats::calculateMean(dataset::extractStressComponent(mlab, i));
        stressVar[i] =
            math::stats::calculateVariance(dataset::extractStressComponent(mlab, i), meanStress);
    }
    mlabStat.put<Vec1Real>("variance-training-stress", stressVar);
}

void dataset::rescaleEnergy(Record& mlab, const Real& energyFactor)
{
    Vec1ShRec& MLABstrucs = mlab.get<Vec1ShRec>("structures");
    std::for_each(MLABstrucs.begin(),
                  MLABstrucs.end(),
                  [&energyFactor](ShRec& struc) { struc->get<Real>("energy") *= energyFactor; });
}

void dataset::rescaleForces(Record& mlab, const Real& forceFactor)
{
    Vec1ShRec& MLABstrucs = mlab.get<Vec1ShRec>("structures");
    std::for_each(MLABstrucs.begin(),
                  MLABstrucs.end(),
                  [&forceFactor](ShRec& struc)
                  {
                      Vec1Real& forces = struc->get<Vec1Real>("forces");
                      std::for_each(forces.begin(),
                                    forces.end(),
                                    [&forceFactor](Real& f) { f *= forceFactor; });
                  });
}

void dataset::rescaleStress(Record& mlab, const Real& stressFactor)
{
    Vec1ShRec& MLABstrucs = mlab.get<Vec1ShRec>("structures");
    std::for_each(MLABstrucs.begin(),
                  MLABstrucs.end(),
                  [&stressFactor](ShRec& struc)
                  {
                      Vec1Real& stress = struc->get<Vec1Real>("stress");
                      std::for_each(stress.begin(),
                                    stress.end(),
                                    [&stressFactor](Real& s) { s *= stressFactor; });
                  });
}

Real dataset::meanEnergy(const Vec1ShCRec& data)
{
    const Vec1Real energies = extractEnergies(data);

    return math::stats::calculateMean(energies);
}

Real dataset::meanForce(const Vec1ShCRec& data)
{
    const Vec1Real forces = vector_tools::Vec2ToVec1(extractForces(data));

    return math::stats::calculateMean(forces);
}

Real dataset::meanStress(const Vec1ShCRec& data)
{
    const Vec1Real stress = vector_tools::Vec2ToVec1(extractStress(data));

    return math::stats::calculateMean(stress);
}

std::tuple<Real, Real, Real> dataset::meanDatsets(const Vec1ShCRec& data)
{
    return {meanEnergy(data), meanForce(data), meanStress(data)};
}

Real dataset::varianceEnergy(const Vec1ShCRec& data)
{
    const Vec1Real energies = extractEnergies(data);
    const Real     mean = math::stats::calculateMean(energies);

    return math::stats::calculateVariance(energies, mean);
}

Real dataset::varianceForce(const Vec1ShCRec& data)
{
    const Vec1Real forces = vector_tools::Vec2ToVec1(extractForces(data));
    const Real     mean = math::stats::calculateMean(forces);

    return math::stats::calculateVariance(forces, mean);
}

Real dataset::varianceStress(const Vec1ShCRec& data)
{
    const Vec1Real stress = vector_tools::Vec2ToVec1(extractStress(data));
    const Real     mean = math::stats::calculateMean(stress);

    return math::stats::calculateVariance(stress, mean);
}

std::tuple<Real, Real, Real> dataset::varianceDatsets(const Vec1ShCRec& data)
{
    return {varianceEnergy(data), varianceForce(data), varianceStress(data)};
}

Real dataset::rmseEnergy(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predEnergies = extractEnergies(predicted);
    const Vec1Real refEnergies = extractEnergies(reference);

    return math::stats::calculateRMSE(predEnergies, refEnergies);
}

Real dataset::rmseForce(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predForce = vector_tools::Vec2ToVec1(extractForces(predicted));
    const Vec1Real refForce = vector_tools::Vec2ToVec1(extractForces(reference));

    return math::stats::calculateRMSE(predForce, refForce);
}

Real dataset::rmseStress(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predStress = vector_tools::Vec2ToVec1(extractStress(predicted));
    const Vec1Real refStress = vector_tools::Vec2ToVec1(extractStress(reference));

    return math::stats::calculateRMSE(predStress, refStress);
}

std::tuple<Real, Real, Real> dataset::rmseDatsets(const Vec1ShCRec& predicted,
                                                  const Vec1ShCRec& reference)
{
    Real rmseE = rmseEnergy(predicted, reference);
    Real rmseF = rmseForce(predicted, reference);
    Real rmseS = rmseStress(predicted, reference);

    return {rmseE, rmseF, rmseS};
}

Real dataset::maeEnergy(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predEnergies = extractEnergies(predicted);
    const Vec1Real refEnergies = extractEnergies(reference);

    return math::stats::calculateMAE(predEnergies, refEnergies);
}

Real dataset::maeForce(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predForce = vector_tools::Vec2ToVec1(extractForces(predicted));
    const Vec1Real refForce = vector_tools::Vec2ToVec1(extractForces(reference));

    return math::stats::calculateRMSE(predForce, refForce);
}

Real dataset::maeStress(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predStress = vector_tools::Vec2ToVec1(extractStress(predicted));
    const Vec1Real refStress = vector_tools::Vec2ToVec1(extractStress(reference));

    return math::stats::calculateMAE(predStress, refStress);
}

std::tuple<Real, Real, Real> dataset::maeDatsets(const Vec1ShCRec& predicted,
                                                 const Vec1ShCRec& reference)
{
    Real maeE = maeEnergy(predicted, reference);
    Real maeF = maeForce(predicted, reference);
    Real maeS = maeStress(predicted, reference);

    return {maeE, maeF, maeS};
}

Real dataset::pearsonCorrEnergy(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predEnergies = extractEnergies(predicted);
    const Vec1Real refEnergies = extractEnergies(reference);

    return math::stats::pearsonCorr(predEnergies, refEnergies);
}

Real dataset::pearsonCorrForce(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predForce = vector_tools::Vec2ToVec1(extractForces(predicted));
    const Vec1Real refForce = vector_tools::Vec2ToVec1(extractForces(reference));

    return math::stats::pearsonCorr(predForce, refForce);
}

Real dataset::pearsonCorrStress(const Vec1ShCRec& predicted, const Vec1ShCRec& reference)
{
    const Vec1Real predStress = vector_tools::Vec2ToVec1(extractStress(predicted));
    const Vec1Real refStress = vector_tools::Vec2ToVec1(extractStress(reference));

    return math::stats::pearsonCorr(predStress, refStress);
}

std::tuple<Real, Real, Real> dataset::pearsonCorrDatsets(const Vec1ShCRec& predicted,
                                                         const Vec1ShCRec& reference)
{
    Real pearsonEnergy = pearsonCorrEnergy(predicted, reference);
    Real pearsonForce = pearsonCorrForce(predicted, reference);
    Real pearsonStress = pearsonCorrStress(predicted, reference);

    return {pearsonEnergy, pearsonForce, pearsonStress};
}

Real dataset::computeEnergyShiftIsolatedAtom(const Record& mlab, const Vec1Real& atomRefEnergy)
{
    Real energy = 0.0;

    Vec1ShCRec     MLABstrucs = mlab.vcget<Vec1ShRec>("structures");
    const Vec1Real ab_initioEnergies = extractEnergies(mlab);
    for (Size config = 0; config < MLABstrucs.size(); config++)
    {
        ShCRec&        struc = MLABstrucs[config];
        const Vec1Int& numAtomsPerType = struc->cget<Vec1Int>("numAtomsPerType");
        Real           refEnergy = 0;
        for (Size i = 0; i < numAtomsPerType.size(); i++)
        {
            refEnergy += (atomRefEnergy[i] * numAtomsPerType[i]);
        }
        energy += (ab_initioEnergies[config] - refEnergy) / (Real)(struc->cget<Int>("numAtoms"));
    }

    return energy / (Real)(MLABstrucs.size());
}

Real dataset::computeEnergyShiftAverageTrainEnergy(const Record&   mlab,
                                                   const Vec1Real& atomRefEnergy)
{
    Real energy = 0.0;

    Vec1ShCRec     MLABstrucs = mlab.vcget<Vec1ShRec>("structures");
    const Vec1Real ab_initioEnergies = extractEnergies(mlab);
    for (Size config = 0; config < MLABstrucs.size(); config++)
    {
        ShCRec&        struc = MLABstrucs[config];
        const Vec1Int& numAtomsPerType = struc->cget<Vec1Int>("numAtomsPerType");
        Real           refEnergy = 0;
        for (Size i = 0; i < numAtomsPerType.size(); i++)
        {
            refEnergy += (atomRefEnergy[i] * numAtomsPerType[i]);
        }
        energy += (ab_initioEnergies[config] - refEnergy);
    }

    return energy / (Real)(MLABstrucs.size());
}

Vec1Real dataset::prepareEnergyRefitIsolatedAtom(const Record& mlab, const Vec1Real& atomRefEnergy)
{
    Vec1ShCRec MLABstrucs = mlab.vcget<Vec1ShRec>("structures");
    Vec1Real   energy = extractEnergies(mlab);
    Real       energyShift = computeEnergyShiftIsolatedAtom(mlab, atomRefEnergy);
    for (Size config = 0; config < MLABstrucs.size(); config++)
    {
        ShCRec&        struc = MLABstrucs[config];
        const Vec1Int& numAtomsPerType = struc->cget<Vec1Int>("numAtomsPerType");
        Real           refEnergy = 0;
        for (Size i = 0; i < numAtomsPerType.size(); i++)
        {
            refEnergy += (atomRefEnergy[i] * numAtomsPerType[i]);
        }
        energy[config] =
            (energy[config] - refEnergy) / (Real)(struc->cget<Int>("numAtoms")) - energyShift;
    }

    return energy;
}

std::vector<TypeMap> dataset::prepareTypeMaps(const Record& mlab)
{
    Vec1ShCRec           MLABstrucs = mlab.vcget<Vec1ShRec>("structures");
    std::vector<TypeMap> typeMaps(MLABstrucs.size());
    const Vec1String&    allTypes = mlab.cget<Vec1String>("typeNames");
    for (Size i = 0; i < MLABstrucs.size(); i++)
    {
        typeMaps[i].update(allTypes, MLABstrucs[i]->cget<Vec1String>("typeNames"));
    }

    return typeMaps;
}
