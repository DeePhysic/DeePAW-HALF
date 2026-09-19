#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include "Record.hpp"
#include "io.hpp"
#include "types.hpp"

namespace vaspml::settings
{

// TODO: ML_FF record argument should be const, but currently needs to be unit-converted internally.
void updateFromMlff(Record& mlff, Record& settings);
void checkFeatureSupport(const Record& settings);

template<typename T>
void updateSingleKeyFromValue(T value, String key, Record& settings)
{
    T& fromSettings = settings.get<T>(key);
    if (value != fromSettings)
    {
        fromSettings = value;

        Record& states = settings.dget<ShRec>("states");
        if (!states.contains(key)) return;

        Int&      state = states.get<Int>(key);
        const Int stateIncar = static_cast<Int>(io::TagState::INCAR);
        const Int stateIncarAlt = static_cast<Int>(io::TagState::INCAR_ALT);
        if (state == stateIncar || state == stateIncarAlt)
        {
            state = static_cast<Int>(io::TagState::OVERRIDE);
        }
    }

    return;
}

template<typename T>
void updateSingleKey(String key1, const Record& mlff, String key2, Record& settings)
{
    updateSingleKeyFromValue<T>(mlff.cget<T>(key1), key2, settings);

    return;
}

} // namespace vaspml::settings

#endif
