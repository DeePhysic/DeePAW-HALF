#include "MlMPI.hpp"

#include "InterfaceLAMMPS.hpp"

#include "BasisFunctions.hpp"
#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "Predictor.hpp"
#include "Record.hpp"
#include "TypeMap.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "nearest_neighbor.hpp"
#include "settings.hpp"
#include "setup.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm> // for find, max, sort
#include <cmath>     // for sqrt
#include <iterator>  // for distance
#include <sstream>   // for basic_stringstream, basic_ostream
#include <stdexcept> // for runtime_error

using namespace vaspml;

struct InterfaceLAMMPS::Impl : public Frame
{
    /***********************************************************************************************
     * Map from LAMMPS types to VASP type strings.
     *
     * LAMMPS types which are not supposed to be handled by VASP have an empty string entry.
     **********************************************************************************************/
    Vec1String typesToNames;
    /***********************************************************************************************
     * Map from LAMMPS types to VASP types.
     *
     * LAMMPS types which are not supposed to be handled by VASP have a -1 entry.
     **********************************************************************************************/
    Vec1Int typesLtoV;
    /***********************************************************************************************
     * List of all LAMMPS types used.
     *
     * LAMMPS types which are not supposed to be handled by VASP are not in this list.
     **********************************************************************************************/
    Vec1Int typesUsed;
    /***********************************************************************************************
     * Type-ordered sequence of LAMMPS neighbor lists.
     *
     * In LAMMPS on each processor there are `list->inum` neighbor lists for a (sub)set of local
     * atoms. These neighbor lists, as well as the local atoms, are not sorted according to types.
     * However, VASPML requires its neighbor lists to be type sorted for both, the central atoms and
     * also within each neighbor list for the neighbor atoms. This array gives a sequence of
     * neighbor lists for type-sorted local atoms, e.g. in
     *
     *     for (int ii = 0; ii < inum; ++ii)
     *     {
     *         const int i = ilist[order[ii]];
     *     }
     *
     * the local atom indices `i` are type-sorted according to ascending VASP force field type
     * numbers.
     **********************************************************************************************/
    Vec1Size order;
    //**********************************************************************************************
    // Copies of pointers to LAMMPS data.
    //**********************************************************************************************
    /// Atom positions.
    const double* const* x = nullptr;
    /// Atom forces.
    double* const* f = nullptr;
    /// Atom types.
    const int* type = nullptr;
    /// Atom global tags.
    const int* tag = nullptr;
    /// Atom global tags (8-byte integer version, LAMMPS_BIGBIG).
    const int64_t* tag64 = nullptr;
    /// Number of central atoms in neighbor list.
    int inum = 0;
    /// Central atom indices in x, type and tag list.
    const int* ilist = nullptr;
    /// Number of neighbors for each central atom in neighbor list.
    const int* numneigh = nullptr;
    /// Indices of neighbors in x, type and tag arrays.
    const int* const* firstneigh = nullptr;
    /// Integer mask to remove special neighbor bits.
    int neighmask = 0;
    /***********************************************************************************************
     * Force field used for predictions.
     **********************************************************************************************/
    Record forceFieldData;
    /***********************************************************************************************
     * Basis functions required for predictions.
     **********************************************************************************************/
    BasisFunctionMap basisFunctions;
    /***********************************************************************************************
     * Kernel instance required for predictions.
     **********************************************************************************************/
    std::unique_ptr<KernelPolynomial> kernel;
    /***********************************************************************************************
     * Predictor class instance required for predictions.
     **********************************************************************************************/
    std::unique_ptr<Predictor> predictor;
    /*******************************************************************************************
     * shared pointer to vaspml main communicator. This coommunicator is a copy or lammps world
     * communicator
     *******************************************************************************************/
    std::shared_ptr<MlMPI> mpiMain;

    /***********************************************************************************************
     * Set up internal type map.
     **********************************************************************************************/
    void setTypeMap(const Vec1String& types, const Vec1String& subtypes);
    /***********************************************************************************************
     * Compute neighbor lists required for all descriptors.
     *
     * @param x 2-dimensional array containing the positions of all local + ghost atoms.
     * @param type 1-dimensional array containing the LAMMPS types of all local + ghost atoms.
     * @param inum Number of atoms for which neighbor list entries are present (this may be smaller
     *             than `nlocal` in case there are atoms which are not mapped to types in the force
     *             field, i.e., have `NULL` entry in `pair_coeff` line).
     * @param ilist 1-dimensional array containing the indices of central atoms (in the x and type
     *              arrays) for which neighbor lists are present.
     * @param numneigh 1-dimensional array containing the number of neighbors for each central atom.
     * @param firstneigh 2-dimensional array containing the indices of neighbor atoms (in the x and
     *                   type arrays) for each central atom (index 1: central atom index (result of
     *                   ilist), index 2: neighbor index).
     * @param neighmask Integer mask to remove bits which mark special neighbors in LAMMPS.
     **********************************************************************************************/
    void computeNeighborLists(const double* const* const& x,
                              const int* const&           type,
                              const int&                  inum,
                              const int* const&           ilist,
                              const int* const&           numneigh,
                              const int* const* const&    firstneigh,
                              const int                   neighmask);
    /***********************************************************************************************
     * Sum up forces from pair and central forces contributions.
     **********************************************************************************************/
    void getForces(double* const* const& f) const;
};

InterfaceLAMMPS::InterfaceLAMMPS() : pImpl{std::make_unique<Impl>()}
{}

InterfaceLAMMPS::~InterfaceLAMMPS() = default;

// Deleted in header, we do not allow copying.
//InterfaceLAMMPS::InterfaceLAMMPS(InterfaceLAMMPS const& other) :
//    pImpl{std::make_unique<Impl>(*other.pImpl)}
//{
//}

// Deleted in header, we do not allow copying.
//InterfaceLAMMPS& InterfaceLAMMPS::operator=(InterfaceLAMMPS const& other)
//{
//    *pImpl = *other.pImpl;
//    return *this;
//}

InterfaceLAMMPS::InterfaceLAMMPS(InterfaceLAMMPS&& other) :
    pImpl{std::make_unique<Impl>(std::move(*other.pImpl))}
{}

InterfaceLAMMPS& InterfaceLAMMPS::operator=(InterfaceLAMMPS&& other)
{
    *pImpl = std::move(*other.pImpl);
    return *this;
}

void InterfaceLAMMPS::Impl::setTypeMap(const Vec1String& types, const Vec1String& subtypes)
{
    if (typeMap != nullptr)
    {
        throw std::runtime_error("ERROR: TypeMap already set, reassignment is not supported.");
    }
    typeMap = std::make_shared<TypeMap>(types, subtypes);

    return;
}

double InterfaceLAMMPS::setupForceField(std::string forceFieldFileName, const MPI_Comm world)
{
    pImpl->mpiMain = std::make_shared<MlMPI>(world, false);

    Record& ffd = pImpl->forceFieldData;
    io::readMlffAndConvertUnits(ffd, forceFieldFileName);
    // NOTE: The feature support check is acting directly on the force field record here.
    settings::checkFeatureSupport(ffd);

    pImpl->basisFunctions = setup::makeBasisFunctions(ffd);

    pImpl->init(ffd, pImpl->basisFunctions);
    pImpl->kernel = std::make_unique<KernelPolynomial>(ffd, pImpl->mpiMain);
    pImpl->predictor = std::make_unique<Predictor>(ffd, pImpl->mpiMain);

    return std::max(ffd.cget<Real>("ML_RCUT1"), ffd.cget<Real>("ML_RCUT2")) * constants::AUTOA;
}

std::string InterfaceLAMMPS::setupTypes(std::vector<std::string> typesToNames)
{
    Vec1String namesFFD = pImpl->forceFieldData.cget<Vec1String>("types");
    string_tools::trimVector(namesFFD);
    Vec1String        namesUsed;
    Vec1Int&          typesUsed = pImpl->typesUsed;
    typesUsed.clear();
    Vec1Int& typesLtoV = pImpl->typesLtoV;
    typesLtoV.clear();

    // Loop over all LAMMPS type names, check if they are present and store map to vaspml types.
    for (Size i = 0; i < typesToNames.size(); ++i)
    {
        auto const& name = string_tools::trim(typesToNames[i]);
        // Ignore unmapped LAMMPS types, log -1 in mapping.
        if (name == "")
        {
            typesLtoV.push_back(-1);
            continue;
        }
        // Search for name of LAMMPS type in force field type names.
        auto const it = std::find(namesFFD.begin(), namesFFD.end(), name);
        // If not found, throw an error.
        if (it == namesFFD.end())
        {
            throw std::runtime_error("ERROR: Type with name \"" + name
                                     + "\" is not present in force field file.");
        }
        // Compute index in force field type list and add value to mapping.
        typesLtoV.push_back(std::distance(namesFFD.begin(), it));
        // Store LAMMPS type index in list of all used indices.
        typesUsed.push_back(i);
        // If type is not yet in list of all used types, add it.
        if (std::find(namesUsed.begin(), namesUsed.end(), name) == namesUsed.end())
        {
            namesUsed.push_back(name);
        }
    }

    // At this point we have two pieces of information (see header for example details).
    // * A list of unique type names in LAMMPS type order, e.g. {"Cs", "Br", "Pb"}.
    // * An "uncompressed" map to all vaspml types defined in the force field, in this example
    //   {-1, 4, -1, 3, 1, 3}. The order of indices corresponds to the one in the force field file.
    // Now, first we must bring the first list into the order of types in the force field file, this
    // will be our subtype list. Second, we need to update the map so it refers to the subtypes
    // instead of the force field types. Finally, for the example the result would look like this:
    // * Subtype list = {"Pb", "Br", "Cs"}
    // * Map from LAMMPS types to VASP subtypes = {-1, 2, -1, 1, 0, 1}

    // First, reorder subtypes.
    std::sort(namesUsed.begin(),
              namesUsed.end(),
              [&ffd = namesFFD](String a, String b)
              {
                  return std::distance(ffd.begin(), std::find(ffd.begin(), ffd.end(), a))
                       < std::distance(ffd.begin(), std::find(ffd.begin(), ffd.end(), b));
              });

    // Second, compress mapping indices to refer to subtype list.
    // Loop over all subtype indices.
    for (Int i = 0; i < static_cast<Int>(namesUsed.size()); ++i)
    {
        // If index i is not found (i.e. there is a gap), decrease all indices greater than i and
        // reset loop index to check again. Entries containing -1 are ignored.
        if (std::find(typesLtoV.begin(), typesLtoV.end(), i) == typesLtoV.end())
        {
            for (auto& m : typesLtoV)
            {
                if (m > i) m--;
            }
            i--;
        }
    }

    // Finally, set up the vaspml TypeMap.
    pImpl->setTypeMap(namesFFD, namesUsed);
    const TypeMap& typeMap = *pImpl->get_typeMap();

    std::stringstream report;
    report << "   LAMMPS       pair_coeff      VASP      |             VASP force field\n";
    report << "    types       names           subtypes  |     types       names        subtypes\n";
    report << "----------------------------------------- | -------------------------------------\n";
    for (Size i = 0; i < std::max(typesToNames.size() - 1, namesFFD.size()); ++i)
    {
        if (i + 1 < typesToNames.size())
        {
            report << str("%9zu <---> ", i + 1);
            if (typesToNames[i + 1].empty()) report << "unmapped! <---> unmapped!";
            else
            {
                report << str("%-9s <---> %-9d", typesToNames[i + 1].c_str(), typesLtoV[i + 1] + 1);
            }
        }
        else report << String(41, ' ');

        if (i < namesFFD.size())
        {
            report << " | " << str("%9zu <---> %-6s <---> ", i + 1, namesFFD[i].c_str());
            if (typeMap.toSubType(i) < 0) report << "unused!\n";
            else report << str("%-9d\n", typeMap.toSubType(i) + 1);
        }
        else report << " |\n";
    }

    // Move LAMMPS types list to internal memory.
    pImpl->typesToNames = std::move(typesToNames);

    return report.str();
}

void InterfaceLAMMPS::setTags(const int* const& tag)
{
    pImpl->tag = tag;

    return;
}

void InterfaceLAMMPS::setTags(const int64_t* const& tag)
{
    pImpl->tag64 = tag;

    return;
}

void InterfaceLAMMPS::computeNeighborLists(const double* const* const& x,
                                           const int* const&           type,
                                           const int&                  inum,
                                           const int* const&           ilist,
                                           const int* const&           numneigh,
                                           const int* const* const&    firstneigh,
                                           const int                   neighmask)
{
    // Store LAMMPS data and pointers for later use.
    pImpl->x = x;
    pImpl->type = type;
    pImpl->inum = inum;
    pImpl->ilist = ilist;
    pImpl->numneigh = numneigh;
    pImpl->firstneigh = firstneigh;
    pImpl->neighmask = neighmask;
    // Forward LAMMPS data and compute neighbor lists.
    pImpl->computeNeighborLists(x, type, inum, ilist, numneigh, firstneigh, neighmask);

    return;
}

void InterfaceLAMMPS::Impl::computeNeighborLists(const double* const* const& x,
                                                 const int* const&           type,
                                                 const int&                  inum,
                                                 const int* const&           ilist,
                                                 const int* const&           numneigh,
                                                 const int* const* const&    firstneigh,
                                                 const int                   neighmask)
{
    // Local atoms need to be passed to VASPML in a type-ordered way. Hence, we need to create an
    // index array ordering the LAMMPS neighbor lists.
    Vec1Int typesNeighborLists;
    for (int ii = 0; ii < inum; ++ii) typesNeighborLists.push_back(typesLtoV[type[ilist[ii]]]);
    // For debugging purposes it may be comfortable to sort according to tag.
    // for (int ii = 0; ii < inum; ++ii) typesNeighborLists.push_back(tag[ilist[ii]]);
    order = vector_tools::argSort(typesNeighborLists);

    // Here we bypass the neighbor list class functionality and instead enter the neighbor data
    // directly in the Frame memory. Looping over neighbor lists but actually only use their names
    // to access memory.
    for (auto const& nl : neighborLists)
    {
        // Name of current neighbor list.
        Record& dataNList = data->dget<ShRec>(nl.first + "-neighbor-list");
        // Create aliases for the current neighbor list's data.
        Vec2Int&  globalIndex = dataNList.get<Vec2Int>("globalIndex");
        Vec2Int&  typeIndex = dataNList.get<Vec2Int>("typeIndex");
        Vec1Int&  typeIndexCentral = dataNList.get<Vec1Int>("typeIndexCentral");
        Vec2Real& distances = dataNList.get<Vec2Real>("distances");
        Vec2Real& connectionVector = dataNList.get<Vec2Real>("connectionVector");
        Vec2Real& connectionVectorNormalized =
            dataNList.get<Vec2Real>("connectionVectorNormalized");
        Vec1Int& numberNeighbors = dataNList.get<Vec1Int>("numberNeighbors");
        Vec2Int& numberNeighborsType = dataNList.get<Vec2Int>("numberNeighborsType");
        Vec1Int& nAtomsType = dataNList.get<Vec1Int>("nAtomsType");
        Vec1Int& n_centralAtomIndexPerType = dataNList.get<Vec1Int>("centralAtomIndexPerType");

        // Set internal variables of neighbor lists.
        nl.second->set_nAtoms(inum);
        nl.second->set_nTypes(typeMap->countStructureTypes());
        nl.second->set_totalNumberTypes(typeMap->countStructureTypes());
        // Update all first dimensions with the number of available LAMMPS neighbor lists.
        globalIndex.resize(inum);
        n_centralAtomIndexPerType.resize(inum);
        typeIndex.resize(inum);
        typeIndexCentral.resize(inum);
        distances.resize(inum);
        connectionVector.resize(inum);
        connectionVectorNormalized.resize(inum);
        numberNeighbors.resize(inum);
        numberNeighborsType.resize(inum);
        // The number of atoms per type needs to be reinitialized to zero.
        nAtomsType.clear();
        nAtomsType.resize(typeMap->countStructureTypes(), 0);

        // Use cutoff in LAMMPS units.
        Real rc2 = nl.second->get_cutOff() * constants::AUTOA;
        rc2 *= rc2;
        for (int ii = 0; ii < inum; ++ii)
        {
            // Index of central atom in LAMMPS x and type arrays.
            const int i = ilist[order[ii]];
            // Type index converted from LAMMPS to vasmpl.
            const int ti = typesLtoV[type[i]];
            // Update number of atoms per type array.
            nAtomsType[ti]++;
            // Store central atom type.
            typeIndexCentral[ii] = ti;
            // Reset total number of neighbors.
            numberNeighbors[ii] = 0;
            // Reset number of neighbors per type.
            numberNeighborsType[ii].clear();
            numberNeighborsType[ii].resize(typeMap->countStructureTypes(), 0);
            // Reserve space in per-neighbor quantities according to LAMMPS numneigh list. Actually
            // not all LAMMPS neighbors will be within the cutoff but we use the number as an
            // estimate (upper boundary) for the memory required. Apply clear() so we can use
            // push_back() below when filling in the actual neighbor's data.
            globalIndex[ii].reserve(numneigh[i]);
            globalIndex[ii].clear();
            typeIndex[ii].reserve(numneigh[i]);
            typeIndex[ii].clear();
            distances[ii].reserve(numneigh[i]);
            distances[ii].clear();
            connectionVector[ii].reserve(3 * numneigh[i]);
            connectionVector[ii].clear();
            connectionVectorNormalized[ii].reserve(3 * numneigh[i]);
            connectionVectorNormalized[ii].clear();
            // In VASPML neighbor atoms need to be type-sorted. Hence, we need to create an
            // index array ordering the neighbors within one LAMMPS neighbor list.
            Vec1Int typesNeighbors;
            for (int jj = 0; jj < numneigh[i]; ++jj)
            {
                int j = firstneigh[i][jj];
                j &= neighmask;
                typesNeighbors.push_back(type[j]);
            }
            Vec1Size orderNeighbors = vector_tools::argSort(typesNeighbors);
            // Loop over neighbors in type-sorted order.
            for (int jj = 0; jj < numneigh[i]; ++jj)
            {
                int j = firstneigh[i][orderNeighbors[jj]];
                j &= neighmask;
                const Real dx = x[j][0] - x[i][0];
                const Real dy = x[j][1] - x[i][1];
                const Real dz = x[j][2] - x[i][2];
                const Real d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= rc2)
                {
                    const Real d = std::sqrt(d2);
                    const int  tj = typesLtoV[type[j]];
                    globalIndex[ii].push_back(j);
                    typeIndex[ii].push_back(tj);
                    connectionVector[ii].push_back(dx / constants::AUTOA);
                    connectionVector[ii].push_back(dy / constants::AUTOA);
                    connectionVector[ii].push_back(dz / constants::AUTOA);
                    connectionVectorNormalized[ii].push_back(dx / d);
                    connectionVectorNormalized[ii].push_back(dy / d);
                    connectionVectorNormalized[ii].push_back(dz / d);
                    distances[ii].push_back(d / constants::AUTOA);
                    numberNeighbors[ii]++;
                    numberNeighborsType[ii][tj]++;
                }
            }
        }
        //nl.second->writeListToScreen();
    }

    return;
}

void InterfaceLAMMPS::predict()
{
    if (pImpl->inum <= 0) return;
    pImpl->update();
    pImpl->kernel->updatePolynomialKernel(*pImpl);
    pImpl->predictor->update(*pImpl->kernel);

    return;
}

double InterfaceLAMMPS::getLocalEnergy() const
{
    if (pImpl->inum <= 0) return 0.0;
    return pImpl->predictor->get_totalEnergy() * constants::EUNIT;
}

void InterfaceLAMMPS::getForces(double* const* const& f) const
{
    pImpl->f = f;
    pImpl->getForces(f);

    return;
}

void InterfaceLAMMPS::Impl::getForces(double* const* const& f) const
{
    for (const String& key : constants::descriptorKeyList)
    {
        // Extract neighbor list key out of descriptor key.
        // This is a quick fix and only works for "SHS2-" and "SHS3-" prefixes.
        const String nlkey = key.substr(5);

        Vec2Int const& globalIndex =
            data->dcget<ShRec>(nlkey + "-neighbor-list").cget<Vec2Int>("globalIndex");
        for (Size i = 0; i < globalIndex.size(); ++i)
        {
            // Convert VASPML central atom index i to LAMMPS neighbor list index (order[i]).
            // Then, look up actual index in force array (ilist[order[i]]).
            const int gi = ilist[order[i]];
            f[gi][0] += predictor->get_centralForcesX(key, i) * constants::FUNIT;
            f[gi][1] += predictor->get_centralForcesY(key, i) * constants::FUNIT;
            f[gi][2] += predictor->get_centralForcesZ(key, i) * constants::FUNIT;
            for (Size j = 0; j < globalIndex[i].size(); ++j)
            {
                // Neighbor index in force array was stored in globalIndex array during neighbor
                // list construction.
                const int gj = globalIndex[i][j];
                f[gj][0] += predictor->get_pairForcesX(key, i, j) * constants::FUNIT;
                f[gj][1] += predictor->get_pairForcesY(key, i, j) * constants::FUNIT;
                f[gj][2] += predictor->get_pairForcesZ(key, i, j) * constants::FUNIT;
            }
        }
    }

    return;
}
