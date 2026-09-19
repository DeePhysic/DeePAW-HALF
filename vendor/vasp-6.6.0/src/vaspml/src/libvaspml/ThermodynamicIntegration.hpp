#ifndef THERMODYNAMICINTEGRATION_HPP
#define THERMODYNAMICINTEGRATION_HPP

#include "types.hpp"

namespace vaspml
{

class ThermodynIntegration
{
  public:
    ThermodynIntegration(const bool&       thermodynIntegrationOnOff_in = false,
                         const ShVec1Real& couplingConstants_in = nullptr,
                         const ShVec1Int&  atomList_in = nullptr,
                         const ShVec2Real& clnmCoupling = nullptr);

    bool            get_StatusOnOff() const;
    void            set_clnmCoupling(const Size& centralAtom, const Size& lnm, const Real& value);
    const Vec1Int&  get_atomList() const;
    const Vec1Real& get_couplingConstants() const;

    void allocate_clnmCoupling(const Size nBasis);

  private:
    /// check if thermodynamic integration is switched on or off
    bool thermodynIntegrationOnOff;

    /// coupling parameters
    ShVec2Real clnmCoupling;

    /// check if atom included in thermodynamic coupling
    ShVec1Int atomList;
    // coupling constants for every atom in structure
    ShVec1Real couplingConstants;
};

} //namespace vaspml

#endif

/*
 *
 * notes to myself (jL)
 * I don't need Ferencian mask, because I will put the thermodynamic
 * integration into a seperate function, so I am going to only loop
 * over the needed atoms. These atoms have to be sorted according to the
 * atom number because this makes the memory access nicer
 *
 *
 */
