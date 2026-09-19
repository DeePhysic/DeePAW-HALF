#ifndef NEIGHBORLISTGRACE_HPP
#define NEIGHBORLISTGRACE_HPP

#include "types.hpp"

#include <set>   // for set
#include <tuple> // for tuple

#ifdef VASPML_ENABLE_GRACE
#include <cppflow/cppflow.h>
#else
#include "cppflow_stub.hpp"
#endif

namespace vaspml
{

class NearestNeighborNSquare;
class Structure;

/* maps element string to numerical type index */
typedef std::map<std::string, int32_t> IndexMap;
/* maps element string to half the cutoff radius per element, i.e.
 * a bond is included iff bond_vector[i][j] <= type_cutoff[type[i]] + type_cutoff[type[j]]
 */
typedef std::map<std::string, double> CutoffMap;

class NeighborListGrace
{
  public:
    NeighborListGrace(double    padding_fraction,
                      IndexMap  types2indices = {},
                      CutoffMap type_cutoff = {}) :
        padding_fraction(padding_fraction),
        types2indices(types2indices),
        type_cutoff(type_cutoff) {};

    int32_t pickTotalNeighbors(int32_t tot_neighbors_real);
    void    copyFromVaspml(const Vec1String& typeNames, const NearestNeighborNSquare& nl);
    void    copyFromVaspml(const Structure& structure, const NearestNeighborNSquare& nl);
    std::vector<std::tuple<std::string, cppflow::tensor>> exportToList(const std::string& prefix);
    void                                                  printTable();

  private:
    double            padding_fraction;
    std::set<int32_t> neighbor_list_sizes;

    IndexMap  types2indices;
    CutoffMap type_cutoff;

    // scalars
    int32_t tot_nat;
    int32_t tot_nat_real;
    int32_t tot_neighbors;
    int32_t tot_neighbors_real;

    // per atom arrays
    std::vector<int32_t> structure_index;
    /* +inputs['atomic_mu_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_atomic_mu_i:0 */
    std::vector<int32_t> atomic_mu_i;

    // per bond arrays
    /* +inputs['ind_i'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_i:0 */
    std::vector<int32_t> ind_i;
    /* +inputs['ind_j'] tensor_info:
        dtype: DT_INT32
        shape: (-1)
        name: serving_default_ind_j:0 */
    std::vector<int32_t> ind_j;
    std::vector<int32_t> mu_i;
    std::vector<int32_t> mu_j;
    /* +inputs['bond_vector'] tensor_info:
        dtype: DT_DOUBLE
        shape: (-1, 3)
        name: serving_default_bond_vector:0 */
    std::vector<double> bond_vector;

    // diagnostic only
    std::vector<bool> padding;

    void getNeighborSize(Vec1Int numNeighbors);
    bool checkCutoff(std::string type1, std::string type2, double distance);
};

} // namespace vaspml

#endif
