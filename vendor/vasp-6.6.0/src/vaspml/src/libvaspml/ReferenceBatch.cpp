#include "ReferenceBatch.hpp"

using namespace vaspml;

ShRec vaspml::makeDataMapReferenceBatch(const ShRec& dataIn)
{
    ShRec data;
    if (dataIn == nullptr) data = std::make_shared<Record>();
    else data = dataIn;
    data = assignOrMakeRecord(data,
                              MapString{
                                  {"forces",   "Vec1Real"},
                                  {"energies", "Vec1Real"}
    });

    return data;
}

ReferenceBatch::ReferenceBatch(const ShRec& dataIn) :
    data(makeDataMapReferenceBatch(dataIn)),
    forces(data->get<Vec1Real>("forces")),
    energies(data->get<Vec1Real>("energies"))
{}

Size ReferenceBatch::getNForceEntries()
{
    return forces.size();
}

Size ReferenceBatch::getNEnergyEntries()
{
    return energies.size();
}

const Vec1Real& ReferenceBatch::get_forces()
{
    return forces;
}

const Vec1Real& ReferenceBatch::get_energies()
{
    return energies;
}
