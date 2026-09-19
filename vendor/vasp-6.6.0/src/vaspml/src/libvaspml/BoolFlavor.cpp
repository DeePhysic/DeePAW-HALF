#include "BoolFlavor.hpp"

#include "debug.hpp"

#include <algorithm>
#include <stdexcept>

using namespace vaspml;

text::BoolFlavor::BoolFlavor(Vec1String trues, Vec1String falses, bool fortranStyle) :
    trues(trues),
    falses(falses),
    fortranStyle(fortranStyle)
{}

bool text::BoolFlavor::isBool(const String& input) const
{
    VASPML_DEBUG_L1(
        if (input.empty())
        {
            throw std::runtime_error("Cannot check string for boolean, input string is empty.");
        }
    );
    if (std::find(trues.begin(), trues.end(), input) != trues.end()) return true;
    if (std::find(falses.begin(), falses.end(), input) != falses.end()) return true;
    if (fortranStyle)
    {
        String sub = input.substr(0, 2);
        if (sub[0] == '.') sub.erase(sub.begin());
        if (sub.empty()) return false;
        if (sub[0] == 't' || sub[0] == 'T' || sub[0] == 'f' || sub[0] == 'F') return true;
    }

    return false;
}

bool text::BoolFlavor::convert(const String& input) const
{
    VASPML_DEBUG_L1(
        if (input.empty())
        {
            throw std::runtime_error("Cannot convert boolean, input string is empty.");
        }
    );
    if (std::find(trues.begin(), trues.end(), input) != trues.end()) return true;
    if (std::find(falses.begin(), falses.end(), input) != falses.end()) return false;
    if (fortranStyle)
    {
        String sub = input.substr(0, 2);
        if (sub[0] == '.') sub.erase(sub.begin());
        if (!sub.empty())
        {
            if (sub[0] == 't' || sub[0] == 'T') return true;
            if (sub[0] == 'f' || sub[0] == 'F') return false;
        }
    }

    throw std::runtime_error("String \"" + input + "\" cannot be parsed as a boolean ("
                             + flf(VASPML_FLF) + ").");
}

text::BoolFlavorJson::BoolFlavorJson() : BoolFlavor({"true"}, {"false"}, false)
{}

text::BoolFlavorExtXyz::BoolFlavorExtXyz() :
    BoolFlavor({"T", "TRUE", "True", "true"}, {"F", "FALSE", "False", "false"}, false)
{}

text::BoolFlavorIncar::BoolFlavorIncar() : BoolFlavor({}, {}, true)
{}
