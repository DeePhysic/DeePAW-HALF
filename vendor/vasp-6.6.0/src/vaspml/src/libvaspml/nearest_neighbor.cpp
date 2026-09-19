#include "nearest_neighbor.hpp"

#include "Lattice.hpp"
#include "ParallelEnvironment.hpp"
#include "Record.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "utils.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_map>

using namespace vaspml;

MapString vaspml::neighbor_list::dataMapNearestNeighborNSquare()
{
    return MapString{
        {"globalIndex",                "Vec2Int"   },
        {"typeIndex",                  "Vec2Int"   },
        {"typeIndexCentral",           "Vec1Int"   },
        {"typeKeyCentral",             "Vec1String"},
        {"distances",                  "Vec2Real"  },
        {"connectionVector",           "Vec2Real"  },
        {"connectionVectorNormalized", "Vec2Real"  },
        {"numberNeighbors",            "Vec1Int"   },
        {"numberNeighborsType",        "Vec2Int"   },
        {"nAtomsType",                 "Vec1Int"   },
        {"centralAtomIndexPerType",    "Vec1Int"   }
    };
}

NearestNeighborNSquare::NearestNeighborNSquare(ShRec record) :
    data(assignOrMakeRecord(record, neighbor_list::dataMapNearestNeighborNSquare())),
    globalIndex(data->get<Vec2Int>("globalIndex")),
    typeIndex(data->get<Vec2Int>("typeIndex")),
    typeIndexCentral(data->get<Vec1Int>("typeIndexCentral")),
    typeKeyCentral(data->get<Vec1String>("typeKeyCentral")),
    distances(data->get<Vec2Real>("distances")),
    connectionVector(data->get<Vec2Real>("connectionVector")),
    connectionVectorNormalized(data->get<Vec2Real>("connectionVectorNormalized")),
    numberNeighbors(data->get<Vec1Int>("numberNeighbors")),
    numberNeighborsType(data->get<Vec2Int>("numberNeighborsType")),
    nAtomsType(data->get<Vec1Int>("nAtomsType")),
    centralAtomIndexPerType(data->get<Vec1Int>("centralAtomIndexPerType"))
{}

NearestNeighborNSquare::NearestNeighborNSquare(Real  Input_CutOff,
                                               bool  isTypeSorted_in,
                                               bool  isDistSorted_in,
                                               ShRec record) :
    data(assignOrMakeRecord(record, neighbor_list::dataMapNearestNeighborNSquare())),
    globalIndex(data->get<Vec2Int>("globalIndex")),
    typeIndex(data->get<Vec2Int>("typeIndex")),
    typeIndexCentral(data->get<Vec1Int>("typeIndexCentral")),
    typeKeyCentral(data->get<Vec1String>("typeKeyCentral")),
    distances(data->get<Vec2Real>("distances")),
    connectionVector(data->get<Vec2Real>("connectionVector")),
    connectionVectorNormalized(data->get<Vec2Real>("connectionVectorNormalized")),
    numberNeighbors(data->get<Vec1Int>("numberNeighbors")),
    numberNeighborsType(data->get<Vec2Int>("numberNeighborsType")),
    nAtomsType(data->get<Vec1Int>("nAtomsType")),
    centralAtomIndexPerType(data->get<Vec1Int>("centralAtomIndexPerType"))
{
    cutOff = Input_CutOff;
    cutOffSquared = Input_CutOff * Input_CutOff;
    isTypeSorted = isTypeSorted_in;
    isDistSorted = isDistSorted_in;
    eps = (Real)1e-6;
    if (isTypeSorted && isDistSorted)
    {
        global::tutor.bug("ERROR: Type- and distance-sorted neighbor list not yet implemented.");
    }
    if (record != nullptr)
    {
        compute_unique_types(typeIndexCentral);
        compute_centralAtomIndexPerType();
        nAtoms = centralAtomIndexPerType.size();
    }
}

void NearestNeighborNSquare::init(const Real cutOff,
                                  const bool isTypeSorted_in,
                                  const bool isDistSorted_in)
{
    this->cutOff = cutOff;
    cutOffSquared = cutOff * cutOff;
    isTypeSorted = isTypeSorted_in;
    isDistSorted = isDistSorted_in;
    eps = (Real)1e-6;

    return;
}

void NearestNeighborNSquare::computeNearestNeighborsDirectCoordinates(const Structure& positions,
                                                                      const Vec1Int    centralAtoms)
{
    VASPML_DEBUG_L1(
        if (!positions.isDirect())
            global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                              + "Positions are not in direct coordinates!");
    );
    nAtoms = centralAtoms.size();
    allocateStorageArrays(nAtoms);
    typeIndexCentral = positions.get_types();
    typeKeyCentral = positions.get_atomTypeNames();
    compute_unique_types(typeIndexCentral);
    for (Size indx = 0; indx < centralAtoms.size(); indx++)
    {
        computeNearestNeighborsSingleAtomDirect(centralAtoms[indx], indx, positions);
    }
    typeIndexCentral = vector_tools::getElements(positions.get_types(), centralAtoms);
    typeKeyCentral = vector_tools::getElements(positions.get_atomTypeNames(), centralAtoms);
    compute_unique_types(typeIndexCentral);
    compute_centralAtomIndexPerType();
    totalNumberTypes = positions.get_typeNames().size();

    return;
}

void NearestNeighborNSquare::computeNearestNeighborsDirectCoordinates(const Structure& positions)
{
    nAtoms = positions.get_numAtoms();
    allocateStorageArrays(nAtoms);
    typeIndexCentral = positions.get_types();
    typeKeyCentral = positions.get_atomTypeNames();
    compute_unique_types(typeIndexCentral);
    for (Size i = 0; i < nAtoms; i++) { computeNearestNeighborsSingleAtomDirect(i, i, positions); }
    compute_centralAtomIndexPerType();

    return;
}

void NearestNeighborNSquare::computeNearestNeighborsCartesianCoordinates(const Structure& positions,
                                                                         const Vec1Int centralAtoms)
{
    VASPML_DEBUG_L1(
        if (positions.isDirect())
            global::tutor.bug("ERROR:" + flf(VASPML_FLF)
                              + "Positions are not in Cartesian coordinates!");
    );
    nAtoms = centralAtoms.size();
    allocateStorageArrays(nAtoms);

    typeIndexCentral = positions.get_types();
    typeKeyCentral = positions.get_atomTypeNames();
    compute_unique_types(typeIndexCentral);

    for (Size indx = 0; indx < centralAtoms.size(); indx++)
    {
        computeNearestNeighborsSingleAtomCartesian(centralAtoms[indx], indx, positions);
    }

    typeIndexCentral = vector_tools::getElements(positions.get_types(), centralAtoms);
    typeKeyCentral = vector_tools::getElements(positions.get_atomTypeNames(), centralAtoms);
    compute_centralAtomIndexPerType();
    totalNumberTypes = positions.get_typeNames().size();

    return;
}

void NearestNeighborNSquare::computeNearestNeighborsCartesianCoordinates(const Structure& positions)
{
    nAtoms = positions.get_numAtoms();
    allocateStorageArrays(nAtoms);
    typeIndexCentral = positions.get_types();
    typeKeyCentral = positions.get_atomTypeNames();
    compute_unique_types(typeIndexCentral);
    for (Size i = 0; i < nAtoms; i++)
    {
        computeNearestNeighborsSingleAtomCartesian(i, i, positions);
    }
    compute_centralAtomIndexPerType();
    totalNumberTypes = positions.get_typeNames().size();

    return;
}

void NearestNeighborNSquare::compute_centralAtomIndexPerType()
{
    //Size     counter = 0;
    //for (Size type = 0; type < (Size)nTypes; type++)
    //{
    //    for (Size atomPerType = 0; atomPerType < (Size)nAtomsType[type]; atomPerType++)
    //    {
    //        centralAtomIndexPerType[counter] = atomPerType;
    //        counter++;
    //    }
    //}

    centralAtomIndexPerType.resize(typeIndexCentral.size());
    std::unordered_map<Int, Int> occurrenceMap;
    std::transform(typeIndexCentral.cbegin(),
                   typeIndexCentral.cend(),
                   centralAtomIndexPerType.begin(),
                   [&occurrenceMap](const Int val) mutable
                   {
                       Int count = occurrenceMap[val];
                       occurrenceMap[val]++;
                       return count;
                   });

    return;
}

void NearestNeighborNSquare::allocateStorageArrays(Size size)
{
    globalIndex.resize(size);
    typeIndex.resize(size);
    typeKey.resize(size);
    distances.resize(size);
    connectionVector.resize(size);
    connectionVectorNormalized.resize(size);
    numberNeighbors.resize(size);

    typeIndexCentral.resize(size);
    typeKeyCentral.resize(size);
    typeEnd = std::make_shared<Vec2Int>();
    typeEnd->resize(size);
    typeStart = std::make_shared<Vec2Int>();
    typeStart->resize(size);
    centralAtomIndexPerType.resize(size);

    return;
}

void NearestNeighborNSquare::allocateWorkArrays(Vec1Int&    nn_list,
                                                Vec1Int&    types_index,
                                                Vec1String& types_key,
                                                Vec1Real&   dist,
                                                Vec1Real&   dist_vecs,
                                                Vec1Real&   normed_dist,
                                                Size        N)
{
    nn_list.resize(N);
    types_index.resize(N);
    types_key.resize(N);
    dist.resize(N);
    dist_vecs.resize(3 * N);
    normed_dist.resize(3 * N);

    return;
}

void NearestNeighborNSquare::computeNearestNeighborsSingleAtomDirect(const Size       indx,
                                                                     const Size       entryIndex,
                                                                     const Structure& positions)
{
    Vec1Int xyz_min(3);
    Vec1Int xyz_max(3);

    auto& [pos_x, pos_y, pos_z] = positions.getAtom(indx);
    Vec1Real pos_central(3);
    pos_central[0] = pos_x;
    pos_central[1] = pos_y;
    pos_central[2] = pos_z;

    const Lattice& lattice = positions.get_lattice();
    const Lattice& inverse_lattice = positions.get_inverseLattice();
    latticeVolume = lattice.get_volume();

    lattice.timesVectorInPlace(pos_central);

    set_periodic_size_images(pos_central, inverse_lattice, xyz_min, xyz_max);

    // position central atom
    // positions neighbor atom
    Real neighbor_x;
    Real neighbor_y;
    Real neighbor_z;
    // distance vector between central atom and neighbor
    // nearest image convention
    Real delta_x;
    Real delta_y;
    Real delta_z;
    // compute distance vector in cartesian coordinates
    Real distance;

    // define temporrary arrays to store nearest neighbor data
    // for certain atom
    Vec1Int    nn_list;
    Vec1Int    types_index;
    Vec1String types_key;
    Vec1Real   dist;
    Vec1Real   dist_vecs;
    Vec1Real   normed_dist;

    Size nAtomsLoc = positions.get_numAtoms();

    Size multi = (xyz_max[0] - xyz_min[0] + 1) * (xyz_max[1] - xyz_min[1] + 1)
               * (xyz_max[2] - xyz_min[2] + 1);
    // allocate work arrays
    allocateWorkArrays(nn_list,
                       types_index,
                       types_key,
                       dist,
                       dist_vecs,
                       normed_dist,
                       multi * nAtomsLoc);

    //std::cout <<  "COMPUTE NEIGHBOR LIST COMPUTE NEIGHBOR LIST " << std::endl;
    //std::cout <<  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    //std::cout <<  xyz_min[0] << "  " << xyz_max[0] << std::endl;
    //std::cout <<  xyz_min[1] << "  " << xyz_max[1] << std::endl;
    //std::cout <<  xyz_min[2] << "  " << xyz_max[2] << std::endl;
    //std::cout <<  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    Size offset = positions.get_indexOffset();
    Size counter = 0;
    for (Int ix = xyz_min[0]; ix <= xyz_max[0]; ix++)
    {
        for (Int iy = xyz_min[1]; iy <= xyz_max[1]; iy++)
        {
            for (Int iz = xyz_min[2]; iz <= xyz_max[2]; iz++)
            {
                for (Size n_indx = 0; n_indx < nAtomsLoc; n_indx++)
                {
                    //NOTE following workaround is needed to get
                    //compiled code with nc++ (NCC) 5.0.1 (Build 16:03:10 May  8 2023)
                    //auto [ x, y, z ] = positions.get_Atom( n_indx );
                    std::tuple<const Real&, const Real&, const Real&> tmp =
                        positions.getAtom(n_indx);
                    auto [x, y, z] = tmp;

                    neighbor_x = x + (Real)ix;
                    neighbor_y = y + (Real)iy;
                    neighbor_z = z + (Real)iz;

                    delta_x = neighbor_x - pos_x;
                    delta_y = neighbor_y - pos_y;
                    delta_z = neighbor_z - pos_z;
                    lattice.timesVectorInPlace(delta_x, delta_y, delta_z);

                    distance = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

                    // collecting data
                    if (distance <= cutOffSquared && distance > eps)
                    {
                        //std::cout << n_indx << "   " << std::sqrt( distance ) << std::endl;
                        nn_list[counter] = n_indx + offset;
                        types_index[counter] = typeIndexCentral[n_indx + offset];
                        types_key[counter] = typeKeyCentral[n_indx + offset];
                        dist[counter] = std::sqrt(distance);
                        dist_vecs[3 * counter] = delta_x;
                        dist_vecs[3 * counter + 1] = delta_y;
                        dist_vecs[3 * counter + 2] = delta_z;
                        normed_dist[3 * counter] = delta_x / dist[counter];
                        normed_dist[3 * counter + 1] = delta_y / dist[counter];
                        normed_dist[3 * counter + 2] = delta_z / dist[counter];
                        counter++;
                    }
                }
            }
        }
    }
    fillNeighborArrays(indx + offset,
                       entryIndex + offset,
                       counter,
                       nn_list,
                       types_index,
                       types_key,
                       dist,
                       dist_vecs,
                       normed_dist);

    return;
}

void NearestNeighborNSquare::computeNearestNeighborsSingleAtomCartesian(const Size       indx,
                                                                        const Size       entryIndex,
                                                                        const Structure& positions)
{
    Vec1Int xyz_min(3);
    Vec1Int xyz_max(3);
    auto& [pos_x, pos_y, pos_z] = positions.getAtom(indx);

    Vec1Real pos_central(3);
    pos_central[0] = pos_x;
    pos_central[1] = pos_y;
    pos_central[2] = pos_z;
    const Lattice& lattice = positions.get_lattice();
    const Lattice& inverse_lattice = positions.get_inverseLattice();
    latticeVolume = lattice.get_volume();

    set_periodic_size_images(pos_central, inverse_lattice, xyz_min, xyz_max);

    Real neighbor_x;
    Real neighbor_y;
    Real neighbor_z;
    // distance vector between central atom and neighbor
    // nearest image convention
    Real delta_x;
    Real delta_y;
    Real delta_z;
    // compute distance vector in cartesian coordinates
    Real distance;

    // define temporrary arrays to store nearest neighbor data
    // for certain atom
    Vec1Int    nn_list;
    Vec1Int    types_index;
    Vec1String types_key;
    Vec1Real   dist;
    Vec1Real   dist_vecs;
    Vec1Real   normed_dist;

    Real shift_x;
    Real shift_y;
    Real shift_z;

    Size nAtomsLoc = positions.get_numAtoms();
    Size multi = (xyz_max[0] - xyz_min[0] + 1) * (xyz_max[1] - xyz_min[1] + 1)
               * (xyz_max[2] - xyz_min[2] + 1);

    // allocate work arrays
    allocateWorkArrays(nn_list,
                       types_index,
                       types_key,
                       dist,
                       dist_vecs,
                       normed_dist,
                       multi * nAtomsLoc);
    Size offset = positions.get_indexOffset();
    Size counter = 0;
    for (Int ix = xyz_min[0]; ix <= xyz_max[0]; ix++)
    {
        for (Int iy = xyz_min[1]; iy <= xyz_max[1]; iy++)
        {
            for (Int iz = xyz_min[2]; iz <= xyz_max[2]; iz++)
            {
                for (Size n_indx = 0; n_indx < nAtomsLoc; n_indx++)
                {
                    //NOTE following workaround is needed to get
                    //compiled code with nc++ (NCC) 5.0.1 (Build 16:03:10 May  8 2023)
                    //auto [ x, y, z ] = positions.get_Atom( n_indx );
                    std::tuple<const Real&, const Real&, const Real&> tmp =
                        positions.getAtom(n_indx);
                    auto [x, y, z] = tmp;

                    shift_x = (Real)ix;
                    shift_y = (Real)iy;
                    shift_z = (Real)iz;

                    lattice.timesVectorInPlace(shift_x, shift_y, shift_z);

                    neighbor_x = x + shift_x;
                    neighbor_y = y + shift_y;
                    neighbor_z = z + shift_z;

                    delta_x = neighbor_x - pos_x;
                    delta_y = neighbor_y - pos_y;
                    delta_z = neighbor_z - pos_z;
                    distance = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
                    if (distance <= cutOffSquared && distance > eps)
                    {
                        nn_list[counter] = n_indx + offset;
                        types_index[counter] = typeIndexCentral[n_indx + offset];
                        types_key[counter] = typeKeyCentral[n_indx + offset];
                        dist[counter] = std::sqrt(distance);
                        dist_vecs[3 * counter] = delta_x;
                        dist_vecs[3 * counter + 1] = delta_y;
                        dist_vecs[3 * counter + 2] = delta_z;
                        normed_dist[3 * counter] = delta_x / dist[counter];
                        normed_dist[3 * counter + 1] = delta_y / dist[counter];
                        normed_dist[3 * counter + 2] = delta_z / dist[counter];
                        counter++;
                    }
                }
            }
        }
    }
    fillNeighborArrays(indx + offset,
                       entryIndex + offset,
                       counter,
                       nn_list,
                       types_index,
                       types_key,
                       dist,
                       dist_vecs,
                       normed_dist);

    return;
}

void NearestNeighborNSquare::set_periodic_size_images(const Vec1Real& pos,
                                                      const Lattice&  inverse_lattice,
                                                      Vec1Int&        imin,
                                                      Vec1Int&        imax)
{
    Vec1Real rr(3);
    Vec1Real xyz_max(9);
    Vec1Real xyz_min(9);

    for (Size ixyz = 0; ixyz < 3; ixyz++)
    {
        rr[ixyz] =
            std::sqrt(inverse_lattice.component(0, ixyz) * inverse_lattice.component(0, ixyz)
                      + inverse_lattice.component(1, ixyz) * inverse_lattice.component(1, ixyz)
                      + inverse_lattice.component(2, ixyz) * inverse_lattice.component(2, ixyz));
    }
    Size counter = 0;
    for (Size ixyz = 0; ixyz < 3; ixyz++)
    {
        for (Size jxyz = 0; jxyz < 3; jxyz++)
        {
            xyz_max[counter] =
                pos[jxyz] + cutOff * inverse_lattice.component(jxyz, ixyz) / rr[jxyz];
            xyz_min[counter] =
                pos[jxyz] - cutOff * inverse_lattice.component(jxyz, ixyz) / rr[jxyz];
            counter++;
        }
    }

    for (Size ixyz = 0; ixyz < 3; ixyz++)
    {
        imax[ixyz] = ceil(inverse_lattice.component(0, ixyz) * xyz_max[3 * ixyz]
                          + inverse_lattice.component(1, ixyz) * xyz_max[3 * ixyz + 1]
                          + inverse_lattice.component(2, ixyz) * xyz_max[3 * ixyz + 2])
                   - 1;
        imin[ixyz] = ceil(inverse_lattice.component(0, ixyz) * xyz_min[3 * ixyz]
                          + inverse_lattice.component(1, ixyz) * xyz_min[3 * ixyz + 1]
                          + inverse_lattice.component(2, ixyz) * xyz_min[3 * ixyz + 2])
                   - 1;
    }

    return;
}

void NearestNeighborNSquare::fillNeighborArrays(const Size        index,
                                                const Size        entryIndex,
                                                const Size        nelements,
                                                const Vec1Int&    nn_list,
                                                const Vec1Int&    types_index,
                                                const Vec1String& types_key,
                                                const Vec1Real&   dist,
                                                const Vec1Real&   dist_vecs,
                                                const Vec1Real&   normed_dist)
{

    Vec1Size indx_array;
    if (isTypeSorted) { indx_array = vector_tools::argSortSlice(types_index, 0, nelements); }
    else if (isDistSorted) { indx_array = vector_tools::argSortSlice(dist, 0, nelements); }
    else { indx_array = vector_tools::generateIntSequence((Size)0, (Size)nelements); }

    globalIndex[entryIndex].resize(nelements);
    typeIndex[entryIndex].resize(nelements);
    typeKey[entryIndex].resize(nelements);
    distances[entryIndex].resize(nelements);
    connectionVector[entryIndex].resize(3 * nelements);
    connectionVectorNormalized[entryIndex].resize(3 * nelements);

    for (Size i = 0; i < (Size)indx_array.size(); i++)
    {
        Size myElement = indx_array[i];
        globalIndex[entryIndex][i] = nn_list[myElement];
        typeIndex[entryIndex][i] = types_index[myElement];
        typeKey[entryIndex][i] = types_key[myElement];
        distances[entryIndex][i] = dist[myElement];
        connectionVector[entryIndex][i * 3] = dist_vecs[myElement * 3];
        connectionVector[entryIndex][i * 3 + 1] = dist_vecs[myElement * 3 + 1];
        connectionVector[entryIndex][i * 3 + 2] = dist_vecs[myElement * 3 + 2];
        connectionVectorNormalized[entryIndex][i * 3] = normed_dist[myElement * 3];
        connectionVectorNormalized[entryIndex][i * 3 + 1] = normed_dist[myElement * 3 + 1];
        connectionVectorNormalized[entryIndex][i * 3 + 2] = normed_dist[myElement * 3 + 2];
        numberNeighborsType[index][unique_types[typeIndex[entryIndex][i]]]++;
    }

    numberNeighbors[entryIndex] = nelements;
    for (Size i = 1; i < (*typeStart)[entryIndex].size(); i++)
    {
        (*typeStart)[entryIndex][i] = (*typeEnd)[entryIndex][i - 1];
    }

    return;
}

void NearestNeighborNSquare::compute_unique_types(const Vec1Int& types)
{
    Vec1Int copy_types = vector_tools::get_unique(types);
    numberNeighborsType.resize(types.size());
    unique_types.clear();
    for (Size i = 0; i < copy_types.size(); i++) { unique_types.insert({copy_types[i], (Int)i}); }

    nTypes = copy_types.size();
    totalNumberTypes = nTypes;
    for (Size i = 0; i < numberNeighborsType.size(); i++) { numberNeighborsType[i].resize(nTypes); }
    nAtomsType.resize(nTypes);
    // count number of atoms per type
    copy_types.resize(types.size());
    for (Size i = 0; i < (Size)types.size(); i++) { copy_types[i] = types[i]; }

    for (Size i = 0; i < (Size)nTypes; i++)
    {
        nAtomsType[i] = vector_tools::count_element(copy_types, unique_types[i]);
    }

    return;
}

void NearestNeighborNSquare::computeNumberAtomsPerType(const Vec1String& unique)
{
    nTypes = unique.size();
    nAtomsType.resize(nTypes);
    for (Int type = 0; type < nTypes; type++)
    {
        nAtomsType[type] = std::count(typeKeyCentral.cbegin(), typeKeyCentral.cend(), unique[type]);
    }
    numberNeighborsType.resize(typeIndex.size());
    for (Size i = 0; i < numberNeighborsType.size(); i++)
        numberNeighborsType[i].resize(unique.size());

    return;
}

void NearestNeighborNSquare::writeListToScreen() const
{
    for (Size i = 0; i < globalIndex.size(); i++)
    {
        // Vec1Real d = distances[i];
        // Vec1Int t = typeIndex[i];
        // std::sort(d.begin(), d.end());
        std::cout << "Central atom " << std::setw(7) << i << " " << globalIndex[i].size()
                  << " |   ";
        for (Size j = 0; j < globalIndex[i].size(); j++)
        {
            std::cout << std::setw(7) << globalIndex[i][j];
            //std::cout << str(" %5.3f ", (*distances)[i][j]);
            //std::cout << "Central atom " << std::setw(7) << i
            //          << str(" - %3zu %3d %24.16E\n", j, t[j], d[j]);
            //std::cout << "Central atom " << std::setw(7) << i
            //          << str(" - %3zu %24.16E\n", j, d[j]);
            //std::cout << str(" %24.16E ", d[j]);
            //std::cout << str(" %5.3f ", d[j]);
            //std::cout << str(" %3d ", d[j]);
        }
        std::cout << std::endl;
    }

    return;
}

const Int& NearestNeighborNSquare::get_globalIndex(const Size atomIndx, const Size indxNeigh) const
{
    return globalIndex[atomIndx][indxNeigh];
}

const Vec1Int& NearestNeighborNSquare::get_globalIndex(const Size atomIndx) const
{
    return globalIndex[atomIndx];
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec2Int& NearestNeighborNSquare::get_globalIndex() const
{
    return globalIndex;
}

const Real& NearestNeighborNSquare::get_distances(const Size atomIndx, const Size indxNeigh) const
{
    return distances[atomIndx][indxNeigh];
}

const Vec2Real& NearestNeighborNSquare::get_distances() const
{
    return distances;
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& NearestNeighborNSquare::get_distances(Size atomIndx) const
{
    return distances[atomIndx];
}

const Real& NearestNeighborNSquare::get_connectionVector_x(const Size atomIndx,
                                                           const Size indxNeigh) const
{
    return connectionVector[atomIndx][3 * indxNeigh];
}

const Real& NearestNeighborNSquare::get_connectionVector_y(const Size atomIndx,
                                                           const Size indxNeigh) const
{
    return connectionVector[atomIndx][3 * indxNeigh + 1];
}

const Real& NearestNeighborNSquare::get_connectionVector_z(const Size atomIndx,
                                                           const Size indxNeigh) const
{
    return connectionVector[atomIndx][3 * indxNeigh + 2];
}

const std::tuple<const Real&, const Real&, const Real&>
NearestNeighborNSquare::get_connectionVector(const Size atomIndx, const Size indxNeigh) const
{
    return std::tie(get_connectionVector_x(atomIndx, indxNeigh),
                    get_connectionVector_y(atomIndx, indxNeigh),
                    get_connectionVector_z(atomIndx, indxNeigh));
}

const Vec2Real& NearestNeighborNSquare::get_connectionVector() const
{
    return connectionVector;
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& NearestNeighborNSquare::get_connectionVector(Size atomIndx) const
{
    return connectionVector[atomIndx];
}

const Real& NearestNeighborNSquare::get_connectionVectorNormalized_x(const Size atomIndx,
                                                                     Size       indxNeigh) const
{
    return connectionVectorNormalized[atomIndx][3 * indxNeigh];
}

const Real& NearestNeighborNSquare::get_connectionVectorNormalized_y(const Size atomIndx,
                                                                     Size       indxNeigh) const
{
    return connectionVectorNormalized[atomIndx][3 * indxNeigh + 1];
}

const Real& NearestNeighborNSquare::get_connectionVectorNormalized_z(const Size atomIndx,
                                                                     Size       indxNeigh) const
{
    return connectionVectorNormalized[atomIndx][3 * indxNeigh + 2];
}

const std::tuple<const Real&, const Real&, const Real&>
NearestNeighborNSquare::get_connectionVectorNormalized(const Size atomIndx,
                                                       const Size indxNeigh) const
{
    return std::tie(get_connectionVectorNormalized_x(atomIndx, indxNeigh),
                    get_connectionVectorNormalized_y(atomIndx, indxNeigh),
                    get_connectionVectorNormalized_z(atomIndx, indxNeigh));
}

const Vec2Real& NearestNeighborNSquare::get_connectionVectorNormalized() const
{
    return connectionVectorNormalized;
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Real& NearestNeighborNSquare::get_connectionVectorNormalized(Size atomIndx) const
{
    return connectionVectorNormalized[atomIndx];
}

VASPML_EXEC_SPACE_SPECIFIER
Size NearestNeighborNSquare::get_size(const Size atomIndx) const
{
    return numberNeighbors[atomIndx];
}

const Vec1Int& NearestNeighborNSquare::get_size() const
{
    return numberNeighbors;
}

VASPML_EXEC_SPACE_SPECIFIER
const Int& NearestNeighborNSquare::get_typeIndexCentral(const Size atomIndex) const
{
    return typeIndexCentral[atomIndex];
}

const Vec1Int& NearestNeighborNSquare::get_typeIndexCentral() const
{
    return typeIndexCentral;
}

const String& NearestNeighborNSquare::get_typeKeyCentral(const Size atomIndex) const
{
    return typeKeyCentral[atomIndex];
}

const Vec1String& NearestNeighborNSquare::get_typeKeyCentral() const
{
    return typeKeyCentral;
}

VASPML_EXEC_SPACE_SPECIFIER
const Int& NearestNeighborNSquare::get_typeIndex(const Size atomIndex, const Size indxNeigh) const
{
    return typeIndex[atomIndex][indxNeigh];
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Int& NearestNeighborNSquare::get_typeIndex(const Size atomIndex) const
{
    return typeIndex[atomIndex];
}

const Vec2Int& NearestNeighborNSquare::get_typeIndex() const
{
    return typeIndex;
}

const String& NearestNeighborNSquare::get_typeKey(const Size atomIndex, const Size indxNeigh) const
{
    return typeKey[atomIndex][indxNeigh];
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1String& NearestNeighborNSquare::get_typeKey(const Size atomIndex) const
{
    return typeKey[atomIndex];
}

const Vec2String& NearestNeighborNSquare::get_typeKey() const
{
    return typeKey;
}

const std::tuple<const Int&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&>
NearestNeighborNSquare::get_neighborData(const Size atomIndex, const Size indxNeigh) const
{

    return std::tie(typeIndex[atomIndex][indxNeigh],
                    distances[atomIndex][indxNeigh],
                    get_connectionVector_x(atomIndex, indxNeigh),
                    get_connectionVector_y(atomIndex, indxNeigh),
                    get_connectionVector_z(atomIndex, indxNeigh),
                    get_connectionVectorNormalized_x(atomIndex, indxNeigh),
                    get_connectionVectorNormalized_y(atomIndex, indxNeigh),
                    get_connectionVectorNormalized_z(atomIndex, indxNeigh));
}

Size NearestNeighborNSquare::get_numberAtoms() const
{
    return typeIndex.size();
}

bool NearestNeighborNSquare::is_typeSorted() const
{
    return isTypeSorted;
}

Size NearestNeighborNSquare::get_numberTypes() const
{
    return totalNumberTypes;
}

const Vec1Int& NearestNeighborNSquare::get_nAtomsType() const
{
    return nAtomsType;
}

const Int& NearestNeighborNSquare::get_nAtomsType(const Size type) const
{
    return nAtomsType[type];
}

VASPML_EXEC_SPACE_SPECIFIER
const Size& NearestNeighborNSquare::get_nAtoms() const
{
    return nAtoms;
}

const Vec1Int& NearestNeighborNSquare::get_numberNeighbors() const
{
    return numberNeighbors;
}

const Int& NearestNeighborNSquare::get_numberNeighbors(const Size atomIdx) const
{
    return numberNeighbors[atomIdx];
}
/// number of nearest neighbors per type for central atom
const Vec2Int& NearestNeighborNSquare::get_numberNeighborsType() const
{
    return numberNeighborsType;
}

const Vec1Int& NearestNeighborNSquare::get_numberNeighborsType(const Size atomIdx) const
{
    return numberNeighborsType[atomIdx];
}

const Real& NearestNeighborNSquare::get_latticeVolume() const
{
    return latticeVolume;
}

VASPML_EXEC_SPACE_SPECIFIER
const Vec1Int& NearestNeighborNSquare::get_centralAtomIndexPerType() const
{
    return centralAtomIndexPerType;
}

const Int& NearestNeighborNSquare::get_centralAtomIndexPerType(const Size atomIdx) const
{
    return centralAtomIndexPerType[atomIdx];
}

Real NearestNeighborNSquare::get_cutOff() const
{
    return cutOff;
}
bool NearestNeighborNSquare::get_isTypeSorted() const
{
    return isTypeSorted;
}
bool NearestNeighborNSquare::get_isDistSorted() const
{
    return isDistSorted;
}

void NearestNeighborNSquare::set_nAtoms(const Int nAtoms)
{
    this->nAtoms = nAtoms;
}

void NearestNeighborNSquare::set_nTypes(const Int nTypes)
{
    this->nTypes = nTypes;
}

void NearestNeighborNSquare::set_totalNumberTypes(const Int totalNumberTypes)
{
    this->totalNumberTypes = totalNumberTypes;
}

Vec1String NearestNeighborNSquare::getUniqueCentralTypes() const
{
    return vector_tools::get_unique(typeKeyCentral);
}

VASPML_EXEC_SPACE_SPECIFIER
void NearestNeighborNSquare::detachNeighborArrays()
{
    //this->globalIndex.reset();
    //this->typeIndex.reset();
    //this->typeIndexCentral.reset();
    //this->distances.reset();
    //this->connectionVector.reset();
    //this->connectionVectorNormalized.reset();
    //this->numberNeighbors.reset();
    //this->numberNeighborsType.reset();
    //this->centralAtomIndexPerType.reset();

    return;
}

VASPML_EXEC_SPACE_SPECIFIER
void NearestNeighborNSquare::resetNeighborArrays(const ShVec2Int& /* globalIndex */,
                                                 const ShVec2Int& /* typeIndex */,
                                                 const ShVec1Int& /* typeIndexCentral */,
                                                 const ShVec2Real& /* distances */,
                                                 const ShVec2Real& /* connectionVector */,
                                                 const ShVec2Real& /* connectionVectorNormalized */,
                                                 const ShVec1Int& /* numberNeighbors */,
                                                 const ShVec2Int& /* numberNeighborsType */,
                                                 const ShVec1Int& /* centralAtomIndexPerType */)
{
    //this->globalIndex = globalIndex;
    //this->typeIndex = typeIndex;
    //this->typeIndexCentral = typeIndexCentral;
    //this->distances = distances;
    //this->connectionVector = connectionVector;
    //this->connectionVectorNormalized = connectionVectorNormalized;
    //this->numberNeighbors = numberNeighbors;
    //this->numberNeighborsType = numberNeighborsType;
    //this->centralAtomIndexPerType = centralAtomIndexPerType;

    return;
}

void NearestNeighborNSquare::allocateArraysRefConfWrapper(const Vec1Size& sizes)
{
    Size nAtoms = 0;
    for (auto s : sizes) nAtoms += s;
    globalIndex.resize(nAtoms);
    typeIndex.resize(nAtoms);
    typeIndexCentral.resize(nAtoms);
    distances.resize(nAtoms);
    connectionVector.resize(nAtoms);
    connectionVectorNormalized.resize(nAtoms);
    numberNeighbors.resize(nAtoms);
    numberNeighborsType.resize(nAtoms);
    centralAtomIndexPerType.resize(nAtoms);

    return;
}

void NearestNeighborNSquare::reorderNeighborList(const Vec1Int& newOrder)
{
    vector_tools::inplaceReorderVector(globalIndex, newOrder);
    vector_tools::inplaceReorderVector(typeIndex, newOrder);
    vector_tools::inplaceReorderVector(typeIndexCentral, newOrder);
    vector_tools::inplaceReorderVector(distances, newOrder);
    vector_tools::inplaceReorderVector(connectionVector, newOrder);
    vector_tools::inplaceReorderVector(connectionVectorNormalized, newOrder);
    vector_tools::inplaceReorderVector(numberNeighbors, newOrder);
    vector_tools::inplaceReorderVector(numberNeighborsType, newOrder);
    compute_centralAtomIndexPerType();

    return;
}

void NearestNeighborNSquare::setParameters(Real cutOff, bool isTypeSorted_in, bool isDistSorted_in)
{
    this->cutOff = cutOff;
    this->isTypeSorted = isTypeSorted_in;
    this->isDistSorted = isDistSorted_in;

    return;
}

void NearestNeighborNSquare::setData(const Vec2Int&    globalIndex,
                                     const Vec2Int&    typeIndex,
                                     const Vec2String& typeKey,
                                     const Vec1Int&    typeIndexCentral,
                                     const Vec1String& typeKeyCentral,
                                     const Vec2Real&   distances,
                                     const Vec2Real&   connectionVector,
                                     const Vec2Real&   connectionVectorNormalized,
                                     const Vec1Int&    numberNeighbors,
                                     const Vec2Int&    numberNeighborsType,
                                     const Vec1Int&    centralAtomIndexPerType)
{
    this->globalIndex = globalIndex;
    this->typeIndex = typeIndex;
    this->typeKey = typeKey;
    this->typeIndexCentral = typeIndexCentral;
    this->typeKeyCentral = typeKeyCentral;
    this->distances = distances;
    this->connectionVector = connectionVector;
    this->connectionVectorNormalized = connectionVectorNormalized;
    this->numberNeighbors = numberNeighbors;
    this->numberNeighborsType = numberNeighborsType;
    this->centralAtomIndexPerType = centralAtomIndexPerType;
    nAtoms = globalIndex.size();

    return;
}

void NearestNeighborNSquare::extendData(const Vec2Int&    globalIndex,
                                        const Vec2Int&    typeIndex,
                                        const Vec2String& typeKey,
                                        const Vec1Int&    typeIndexCentral,
                                        const Vec1String& typeKeyCentral,
                                        const Vec2Real&   distances,
                                        const Vec2Real&   connectionVector,
                                        const Vec2Real&   connectionVectorNormalized,
                                        const Vec1Int&    numberNeighbors,
                                        const Vec2Int&    numberNeighborsType,
                                        const Vec1Int&    centralAtomIndexPerType)
{
    this->globalIndex.insert(this->globalIndex.end(), globalIndex.cbegin(), globalIndex.cend());
    this->typeIndex.insert(this->typeIndex.end(), typeIndex.cbegin(), typeIndex.cend());
    this->typeKey.insert(this->typeKey.end(), typeKey.cbegin(), typeKey.cend());
    this->typeIndexCentral.insert(this->typeIndexCentral.end(),
                                  typeIndexCentral.cbegin(),
                                  typeIndexCentral.cend());
    this->typeKeyCentral.insert(this->typeKeyCentral.end(),
                                typeKeyCentral.cbegin(),
                                typeKeyCentral.cend());
    this->distances.insert(this->distances.end(), distances.cbegin(), distances.cend());
    this->connectionVector.insert(this->connectionVector.end(),
                                  connectionVector.cbegin(),
                                  connectionVector.cend());
    this->connectionVectorNormalized.insert(this->connectionVectorNormalized.end(),
                                            connectionVectorNormalized.cbegin(),
                                            connectionVectorNormalized.cend());
    this->numberNeighbors.insert(this->numberNeighbors.end(),
                                 numberNeighbors.cbegin(),
                                 numberNeighbors.cend());
    this->numberNeighborsType.insert(this->numberNeighborsType.end(),
                                     numberNeighborsType.cbegin(),
                                     numberNeighborsType.cend());
    this->centralAtomIndexPerType.insert(this->centralAtomIndexPerType.end(),
                                         centralAtomIndexPerType.cbegin(),
                                         centralAtomIndexPerType.cend());
    nAtoms += globalIndex.size();

    return;
}

void NearestNeighborNSquare::addElement(const Vec1Int&    globalIndex,
                                        const Vec1Int&    typeIndex,
                                        const Vec1String& typeKey,
                                        const Int&        typeIndexCentral,
                                        const String&     typeKeyCentral,
                                        const Vec1Real&   distances,
                                        const Vec1Real&   connectionVector,
                                        const Vec1Real&   connectionVectorNormalized,
                                        const Int&        numberNeighbors,
                                        const Vec1Int&    numberNeighborsType,
                                        const Int&        centralAtomIndexPerType)
{
    this->globalIndex.push_back(globalIndex);
    this->typeIndex.push_back(typeIndex);
    this->typeKey.push_back(typeKey);
    this->typeIndexCentral.push_back(typeIndexCentral);
    this->typeKeyCentral.push_back(typeKeyCentral);
    this->distances.push_back(distances);
    this->connectionVector.push_back(connectionVector);
    this->connectionVectorNormalized.push_back(connectionVectorNormalized);
    this->numberNeighbors.push_back(numberNeighbors);
    this->numberNeighborsType.push_back(numberNeighborsType);
    this->centralAtomIndexPerType.push_back(centralAtomIndexPerType);
    nAtoms++;

    return;
}

bool NearestNeighborNSquare::operator<(const NearestNeighborNSquare& other) const
{
    return this->cutOff < other.cutOff;
}

bool NearestNeighborNSquare::operator>(const NearestNeighborNSquare& other) const
{
    return this->cutOff > other.cutOff;
}

bool NearestNeighborNSquare::operator<=(const NearestNeighborNSquare& other) const
{
    return this->cutOff <= other.cutOff;
}

bool NearestNeighborNSquare::operator==(const NearestNeighborNSquare& other) const
{
    return this->cutOff == other.cutOff;
}

bool NearestNeighborNSquare::operator!=(const NearestNeighborNSquare& other) const
{
    return this->cutOff != other.cutOff;
}

Vec2Int NearestNeighborNSquare::computeNumberNeighborsPerType() const
{
    Vec2Int numberNeighborsPerType(nAtoms);
    std::for_each(numberNeighborsPerType.begin(),
                  numberNeighborsPerType.end(),
                  [&](Vec1Int& slice)
                  {
                      slice.resize(nAtomsType.size());
                      std::fill(slice.begin(), slice.end(), 0);
                  });
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size neighbor = 0; neighbor < typeIndex[atom].size(); neighbor++)
        {
            numberNeighborsPerType[atom][typeIndex[atom][neighbor]]++;
        }
    }
    return numberNeighborsPerType;
}

Vec1Int NearestNeighborNSquare::computeMaxNeighborsPerType() const
{

    Vec1Int maxNeighborsPerType(nAtomsType.size(), -1000);
    Vec2Int neighborsPerType = computeNumberNeighborsPerType();
    for (Size atom = 0; atom < nAtoms; atom++)
    {
        for (Size type = 0; type < neighborsPerType.size(); type++)
        {
            if (maxNeighborsPerType[type] < neighborsPerType[atom][type])
                maxNeighborsPerType[type] = neighborsPerType[atom][type];
        }
    }
    return maxNeighborsPerType;
}

const ShRec& NearestNeighborNSquare::getNeighborRecord() const
{
    return data;
}

ShRec& NearestNeighborNSquare::getNeighborRecord()
{
    return data;
}

namespace vaspml
{

void neighbor_list::copyParameters(const NearestNeighborNSquare& listA,
                                   NearestNeighborNSquare&       listB)
{
    listB.setParameters(listA.get_cutOff(), listA.get_isTypeSorted(), listA.get_isDistSorted());

    return;
}

NearestNeighborNSquare* neighbor_list::copyParameters(const NearestNeighborNSquare& listA)
{
    return new NearestNeighborNSquare(listA.get_cutOff(),
                                      listA.get_isTypeSorted(),
                                      listA.get_isDistSorted());
}

void neighbor_list::copyData(const NearestNeighborNSquare& listA, NearestNeighborNSquare& listB)
{
    listB.setData(listA.get_globalIndex(),
                  listA.get_typeIndex(),
                  listA.get_typeKey(),
                  listA.get_typeIndexCentral(),
                  listA.get_typeKeyCentral(),
                  listA.get_distances(),
                  listA.get_connectionVector(),
                  listA.get_connectionVectorNormalized(),
                  listA.get_numberNeighbors(),
                  listA.get_numberNeighborsType(),
                  listA.get_centralAtomIndexPerType());

    return;
}

void neighbor_list::extendData(const NearestNeighborNSquare& listA, NearestNeighborNSquare& listB)
{
    listB.extendData(listA.get_globalIndex(),
                     listA.get_typeIndex(),
                     listA.get_typeKey(),
                     listA.get_typeIndexCentral(),
                     listA.get_typeKeyCentral(),
                     listA.get_distances(),
                     listA.get_connectionVector(),
                     listA.get_connectionVectorNormalized(),
                     listA.get_numberNeighbors(),
                     listA.get_numberNeighborsType(),
                     listA.get_centralAtomIndexPerType());

    return;
}

void neighbor_list::addElementData(const NearestNeighborNSquare& listA,
                                   NearestNeighborNSquare&       listB,
                                   const Size                    atomIndx)
{
    listB.addElement(listA.get_globalIndex(atomIndx),
                     listA.get_typeIndex(atomIndx),
                     listA.get_typeKey(atomIndx),
                     listA.get_typeIndexCentral(atomIndx),
                     listA.get_typeKeyCentral(atomIndx),
                     listA.get_distances(atomIndx),
                     listA.get_connectionVector(atomIndx),
                     listA.get_connectionVectorNormalized(atomIndx),
                     listA.get_numberNeighbors(atomIndx),
                     listA.get_numberNeighborsType(atomIndx),
                     listA.get_centralAtomIndexPerType(atomIndx));

    return;
}

std::shared_ptr<const NearestNeighborNSquare> neighbor_list::getMaxNeighborListPtr(
    const std::map<String, std::shared_ptr<const NearestNeighborNSquare>>& neighborList)
{

    return std::max_element(neighborList.cbegin(),
                            neighborList.cend(),
                            [](const auto& a, const auto& b) { return (*a.second) < (*b.second); })
        ->second;
}

ShRec neighbor_list::extractNeighborListFromBatchList(const NearestNeighborNSquare& neighborList,
                                                      const AtomBatchMap&           batchMap,
                                                      const Int                     nStruc)
{
    const Int nAtoms = batchMap.get_numberAtomsPerStruc(nStruc);
    ShRec     result;

    // extract
    result = assignOrMakeRecord(result, neighbor_list::dataMapNearestNeighborNSquare());
    Vec2Int& globalIndex = result->get<Vec2Int>("globalIndex");
    globalIndex.resize(nAtoms);
    Vec2Int& typeIndex = result->get<Vec2Int>("typeIndex");
    typeIndex.resize(nAtoms);
    Vec1Int& typeIndexCentral = result->get<Vec1Int>("typeIndexCentral");
    typeIndexCentral.resize(nAtoms);
    Vec1String& typeKeyCentral = result->get<Vec1String>("typeKeyCentral");
    typeKeyCentral.resize(nAtoms);
    Vec2Real& distances = result->get<Vec2Real>("distances");
    distances.resize(nAtoms);
    Vec2Real& connectionVector = result->get<Vec2Real>("connectionVector");
    connectionVector.resize(nAtoms);
    Vec2Real& connectionVectorNormalized = result->get<Vec2Real>("connectionVectorNormalized");
    connectionVectorNormalized.resize(nAtoms);
    Vec1Int& numberNeighbors = result->get<Vec1Int>("numberNeighbors");
    numberNeighbors.resize(nAtoms);
    Vec2Int& numberNeighborsType = result->get<Vec2Int>("numberNeighborsType");
    numberNeighborsType.resize(nAtoms);
    Vec1Int& centralAtomIndexPerType = result->get<Vec1Int>("centralAtomIndexPerType");
    centralAtomIndexPerType.resize(nAtoms);

    for (Int i = 0; i < nAtoms; i++)
    {
        const Size entry = batchMap.get_mapStrucAtom_Batch(nStruc, i);
        globalIndex[i] = neighborList.get_globalIndex(entry);
        typeIndex[i] = neighborList.get_typeIndex(entry);
        typeIndexCentral[i] = neighborList.get_typeIndexCentral(entry);
        distances[i] = neighborList.get_distances(entry);
        connectionVector[i] = neighborList.get_connectionVector(entry);
        connectionVectorNormalized[i] = neighborList.get_connectionVectorNormalized(entry);
        numberNeighbors[i] = neighborList.get_numberNeighbors(entry);
        numberNeighborsType[i] = neighborList.get_numberNeighborsType(entry);
    }
    const Vec1Int& tic = typeIndexCentral;
    Vec1Int        unique = vector_tools::get_unique(tic);
    Vec1Int&       nAtomsType = result->get<Vec1Int>("nAtomsType");
    nAtomsType.resize(unique.size());
    std::transform(unique.cbegin(),
                   unique.cend(),
                   nAtomsType.begin(),
                   [&tic](const Int val) { return std::count(tic.begin(), tic.end(), val); });
    // computing the ordinal vector.
    std::unordered_map<Int, Int> occurrenceMap;
    std::transform(typeIndexCentral.cbegin(),
                   typeIndexCentral.cend(),
                   centralAtomIndexPerType.begin(),
                   [&occurrenceMap](const Int val) mutable
                   {
                       Int count = occurrenceMap[val];
                       occurrenceMap[val]++;
                       return count;
                   });

    return result;
}

NearestNeighborNSquare neighbor_list::getNeighborElements(
    const NearestNeighborNSquare& neighborList,
    const Vec1Int&                atoms)
{
    ShRec result = nullptr;
    result = assignOrMakeRecord(result, neighbor_list::dataMapNearestNeighborNSquare());
    const ShRec& nRecord = neighborList.getNeighborRecord();

    result->put("globalIndex",
                vector_tools::filterVector(nRecord->get<Vec2Int>("globalIndex"), atoms));
    result->put("typeIndex", vector_tools::filterVector(nRecord->get<Vec2Int>("typeIndex"), atoms));
    result->put("typeIndexCentral",
                vector_tools::filterVector(nRecord->get<Vec1Int>("typeIndexCentral"), atoms));
    result->put("typeKeyCentral",
                vector_tools::filterVector(nRecord->get<Vec1String>("typeKeyCentral"), atoms));
    result->put("distances",
                vector_tools::filterVector(nRecord->get<Vec2Real>("distances"), atoms));
    result->put("connectionVector",
                vector_tools::filterVector(nRecord->get<Vec2Real>("connectionVector"), atoms));
    result->put(
        "connectionVectorNormalized",
        vector_tools::filterVector(nRecord->get<Vec2Real>("connectionVectorNormalized"), atoms));
    result->put("numberNeighbors",
                vector_tools::filterVector(nRecord->get<Vec1Int>("numberNeighbors"), atoms));
    result->put("numberNeighborsType",
                vector_tools::filterVector(nRecord->get<Vec2Int>("numberNeighborsType"), atoms));

    NearestNeighborNSquare nListResult(neighborList.get_cutOff(),
                                       neighborList.get_isTypeSorted(),
                                       neighborList.get_isDistSorted(),
                                       result);
    return nListResult;
}

} //namespace vaspml
