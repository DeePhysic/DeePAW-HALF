#include "NeighborListGrace.hpp"

#include "Structure.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <cmath>
#include <iostream>

using namespace vaspml;

/*
 * signature_def['serving_default']:
  The given SavedModel SignatureDef contains the following input(s):
    +inputs['batch_tot_nat'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat:0
    +inputs['batch_tot_nat_real'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_batch_tot_nat_real:0
    +inputs['map_atoms_to_structure'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_map_atoms_to_structure:0
    +inputs['n_struct_total'] tensor_info:
        dtype: DT_INT32
        shape: ()
        name: serving_default_n_struct_total:0


  The given SavedModel SignatureDef contains the following output(s):
    outputs['atomic_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 1)
        name: StatefulPartitionedCall:0
    outputs['total_energy'] tensor_info:
        dtype: DT_DOUBLE
        shape: (1, 1)
        name: StatefulPartitionedCall:1
    outputs['total_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:2
    outputs['virial'] tensor_info:
        dtype: DT_DOUBLE
        shape: (6)
        name: StatefulPartitionedCall:3
    outputs['z_pair_f'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: StatefulPartitionedCall:4
        */

int32_t NeighborListGrace::pickTotalNeighbors(int32_t tot_neighbors_real)
{
    int32_t candidate = static_cast<int32_t>(std::round(tot_neighbors_real * padding_fraction));
    if (auto search = neighbor_list_sizes.upper_bound(tot_neighbors_real);
        search != neighbor_list_sizes.end())
    {
        return *search;
    }
    else
    {
        neighbor_list_sizes.insert(candidate);
        return candidate;
    }
}

void NeighborListGrace::getNeighborSize(Vec1Int numNeighbors)
{
    tot_neighbors_real = 0;
    for (Int n : numNeighbors) tot_neighbors_real += n;
    // TODO: additional logic to reduce the number of recompilations
    tot_neighbors = pickTotalNeighbors(tot_neighbors_real);
    //std::cout << "Picked " << tot_neighbors << " out of" << std::endl;
    //for (auto n : neighbor_list_sizes) {
    //    std::cout << " " << n;
    //}
    //std::cout << std::endl;
    tot_nat_real = numNeighbors.size();
    tot_nat = tot_nat_real + 1;

    atomic_mu_i.resize(tot_nat);
    structure_index.resize(tot_nat);

    ind_i.resize(tot_neighbors);
    ind_j.resize(tot_neighbors);
    mu_i.resize(tot_neighbors);
    mu_j.resize(tot_neighbors);
    bond_vector.resize(3 * tot_neighbors);
    padding.resize(tot_neighbors);

    return;
}

void NeighborListGrace::printTable()
{
    std::cout << "P    i    I    J   MI   MJ B" << std::endl;
    for (size_t i = 0; i < ind_i.size(); ++i)
    {
        std::cout << str("%s %4i %4i %4i %4i %4i %2.4f %2.4f %2.4f",
                         (padding[i] ? "v" : "x"),
                         i,
                         ind_i[i],
                         ind_j[i],
                         mu_i[i],
                         mu_j[i],
                         bond_vector[3 * i + 0],
                         bond_vector[3 * i + 1],
                         bond_vector[3 * i + 2])
                  << std::endl;
    }

    return;
}

bool NeighborListGrace::checkCutoff(std::string type1, std::string type2, double distance)
{
    return distance <= type_cutoff.at(type1) + type_cutoff.at(type2);
}

void NeighborListGrace::copyFromVaspml(const Structure& structure, const NearestNeighborNSquare& nl)
{
    copyFromVaspml(structure.get_typeNames(), nl);

    return;
}

void NeighborListGrace::copyFromVaspml(const Vec1String&             typeNames,
                                       const NearestNeighborNSquare& nl)
{
    /*============================================================================================+
     | Here you can use all kind of getters to retrieve the original structure and neighbor list
     | data (all references, no copying involved).
     +============================================================================================*/
    // From neighbor list:
    const Vec1Int&  numNeighbors = nl.get_size();
    const Vec2Int&  globalIndex = nl.get_globalIndex();
    const Vec2Int&  typeIndex = nl.get_typeIndex();
    const Vec1Int&  typeIndexCentral = nl.get_typeIndexCentral();
    const Vec2Real& distances = nl.get_distances();
    const Vec2Real& connectionVector = nl.get_connectionVector();

    getNeighborSize(numNeighbors);

    double bond_padding = 2 * nl.get_cutOff();

    auto tot_neighbors_real = 0;
    for (Size i = 0; i < globalIndex.size(); ++i)
    {
        atomic_mu_i[i] = types2indices.at(typeNames[typeIndexCentral[i]]);
        for (Size j = 0; j < globalIndex[i].size(); ++j)
        {
            if (!checkCutoff(typeNames[typeIndexCentral[i]],
                             typeNames[typeIndex[i][j]],
                             distances[i][j] * vaspml::constants::AUTOA))
            {
                continue;
            };

            ind_i[tot_neighbors_real] = i;
            ind_j[tot_neighbors_real] = globalIndex[i][j];
            mu_i[tot_neighbors_real] =
                types2indices.at(typeNames[typeIndexCentral[ind_i[tot_neighbors_real]]]);
            mu_j[tot_neighbors_real] =
                types2indices.at(typeNames[typeIndexCentral[ind_j[tot_neighbors_real]]]);
            for (auto n = 0; n < 3; ++n)
            {
                bond_vector[3 * tot_neighbors_real + n] =
                    connectionVector[i][3 * j + n] * vaspml::constants::AUTOA;
            }

            tot_neighbors_real += 1;
        }
    }

    // pad the atom arrays
    for (auto p = tot_nat_real; p < tot_nat; ++p)
    {
        atomic_mu_i[p] = 0;
        structure_index[p] = 1;
    }
    // pad the bond arrays
    for (auto p = tot_neighbors_real; p < tot_neighbors; ++p)
    {
        mu_i[p] = 0;
        mu_j[p] = 0;
        ind_i[p] = numNeighbors.size();
        ind_j[p] = numNeighbors.size();
        bond_vector[3 * p + 0] = bond_padding;
        bond_vector[3 * p + 1] = bond_padding;
        bond_vector[3 * p + 2] = bond_padding;
        padding[p] = true;
    }

    return;
}

std::vector<std::tuple<std::string, cppflow::tensor>> NeighborListGrace::exportToList(
    const std::string& prefix)
{
    std::vector<std::tuple<std::string, cppflow::tensor>> inputs;

    inputs.emplace_back(prefix + "batch_tot_nat:0", cppflow::tensor(tot_nat));
    inputs.emplace_back(prefix + "batch_tot_nat_real:0", cppflow::tensor(tot_nat_real));
    inputs.emplace_back(prefix + "atomic_mu_i" + ":0", cppflow::tensor(atomic_mu_i, {tot_nat}));
    // inputs.emplace_back(prefix + "map_atoms_to_structure" + ":0",
    // cppflow::tensor(structure_index, {tot_nat}));
    inputs.emplace_back(prefix + "bond_vector" + ":0",
                        cppflow::tensor(bond_vector, {tot_neighbors, 3}));
    inputs.emplace_back(prefix + "ind_i" + ":0", cppflow::tensor(ind_i, {tot_neighbors}));
    inputs.emplace_back(prefix + "ind_j" + ":0", cppflow::tensor(ind_j, {tot_neighbors}));
    inputs.emplace_back(prefix + "mu_i" + ":0", cppflow::tensor(mu_i, {tot_neighbors}));
    inputs.emplace_back(prefix + "mu_j" + ":0", cppflow::tensor(mu_j, {tot_neighbors}));

    return inputs;
}
