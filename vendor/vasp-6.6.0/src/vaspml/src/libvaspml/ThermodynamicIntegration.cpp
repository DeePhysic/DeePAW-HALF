#include "ThermodynamicIntegration.hpp"

using namespace vaspml;

ThermodynIntegration::ThermodynIntegration(const bool&       thermodynIntegrationOnOff_in,
                                           const ShVec1Real& couplingConstants_in,
                                           const ShVec1Int&  atomList_in,
                                           const ShVec2Real& clnmCoupling_in)
{
    thermodynIntegrationOnOff = thermodynIntegrationOnOff_in;
    if (!thermodynIntegrationOnOff) return;

    if (couplingConstants_in == nullptr) { couplingConstants = std::make_shared<Vec1Real>(); }
    else { couplingConstants = couplingConstants_in; }

    if (atomList_in == nullptr) { atomList = std::make_shared<Vec1Int>(); }
    else { atomList = atomList_in; }

    if (clnmCoupling == nullptr) { clnmCoupling = std::make_shared<Vec2Real>(); }
    else { clnmCoupling = clnmCoupling_in; }
}

bool ThermodynIntegration::get_StatusOnOff() const
{
    return thermodynIntegrationOnOff;
}

void ThermodynIntegration::set_clnmCoupling(const Size& centralAtom,
                                            const Size& lnm,
                                            const Real& value)
{
    (*clnmCoupling)[centralAtom][lnm] = value;

    return;
}

const Vec1Int& ThermodynIntegration::get_atomList() const
{
    return *atomList;
}

const Vec1Real& ThermodynIntegration::get_couplingConstants() const
{
    return *couplingConstants;
}

void ThermodynIntegration::allocate_clnmCoupling(const Size nsize)
{
    clnmCoupling->resize(atomList->size());
    for (Size i = 0; i < atomList->size(); i++) { (*clnmCoupling)[i].resize(nsize); }

    return;
}
