#include "settings.hpp"

#include "Item.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "setup.hpp"

#include <algorithm> // for max

using namespace vaspml;

void settings::updateFromMlff(Record& mlff, Record& settings)
{
    setup::rescaleMlffUnits(mlff,
                            constants::EUNIT,
                            constants::AUTOA,
                            constants::FUNIT,
                            constants::SUNIT);

    Vec1String keys = settings.keys();

    for (const auto& k : keys)
    {
        if (mlff.contains(k))
        {
            ItemIndex i = settings.itemIndexOf(k);
            if (i == ItemIndex::REAL) updateSingleKey<Real>(k, mlff, k, settings);
            else if (i == ItemIndex::INT) updateSingleKey<Int>(k, mlff, k, settings);
            else if (i == ItemIndex::STRING) updateSingleKey<String>(k, mlff, k, settings);
            else if (i == ItemIndex::BOOL) updateSingleKey<bool>(k, mlff, k, settings);
            else if (i == ItemIndex::VEC1REAL) updateSingleKey<Vec1Real>(k, mlff, k, settings);
            else if (i == ItemIndex::VEC1INT) updateSingleKey<Vec1Int>(k, mlff, k, settings);
            else if (i == ItemIndex::VEC1STRING) updateSingleKey<Vec1String>(k, mlff, k, settings);
            else if (i == ItemIndex::VEC2REAL) updateSingleKey<Vec2Real>(k, mlff, k, settings);
            else if (i == ItemIndex::VEC2INT) updateSingleKey<Vec2Int>(k, mlff, k, settings);
            else if (i == ItemIndex::VEC2STRING) updateSingleKey<Vec1String>(k, mlff, k, settings);
            else { global::tutor.warning("Cannot check/update value for key \"" + k + "\"."); }
        }
    }

    Vec1Int numLrc = mlff.cget<Vec1Int>("numLrc");
    Int     maxLrc = 0;
    for (const Int& i : numLrc) maxLrc = std::max(maxLrc, i);
    updateSingleKeyFromValue(maxLrc, "ML_MB", settings);

    updateSingleKey<Int>("numTypes", mlff, "NTYP", settings);
    updateSingleKey<Vec1String>("types", mlff, "TYPE", settings);

    setup::rescaleMlffUnits(mlff,
                            1.0 / constants::EUNIT,
                            1.0 / constants::AUTOA,
                            1.0 / constants::FUNIT,
                            1.0 / constants::SUNIT);

    return;
}

void settings::checkFeatureSupport(const Record& settings)
{
    /*============================================================================================+
     | NOTE: Usually this function takes in complete setttings from the INCAR reader, modified by
     | force field data, etc. However, it should be designed such that also only the force field
     | data can be processed on its own. This way we can perform checks even if there is no log
     | file necessary (e.g. LAMMPS interface).
     +============================================================================================*/
    if (!settings.cget<bool>("ML_LFAST"))
    {
        global::tutor.error("This force field file does not support the fast prediction mode. "
                            "Please perform a refit with ML_MODE = refit.");
    }

    return;
}
