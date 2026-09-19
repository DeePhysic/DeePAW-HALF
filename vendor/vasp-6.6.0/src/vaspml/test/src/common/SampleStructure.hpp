#ifndef SAMPLESTRUCTURE_HPP
#define SAMPLESTRUCTURE_HPP

#include "Structure.hpp"
#include "types.hpp"

namespace vaspml
{

/** Collection of sample structures
 *
 * The following structures can be selected by passing the string to the setup
 * function or the constructor:
 *   - `CsPbBr3_40` ... Cesium Lead Bromide structure CsPbBr3 with 40 atoms (bulk).
 *   - `CaO_16` ... Calcium oxide CaO structure with 16 atoms (bulk).
 *   - `HCN_24` ... Azobenzene C12H10N2 with 24 atoms (molecule).
 */
struct SampleStructure : public Structure
{
    explicit SampleStructure(String sample);

    void setup(String sample);
};

} //namespace vaspml

#endif
