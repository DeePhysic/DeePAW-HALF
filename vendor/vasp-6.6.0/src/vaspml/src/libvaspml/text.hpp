#ifndef TEXT_HPP
#define TEXT_HPP

#include "Item.hpp"
#include "debug.hpp"
#include "meta.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <exception>    // for exception
#include <stdexcept>    // for invalid_argument, runtime_error
#include <type_traits>  // for enable_if_t, is_base_of, is_base_of_v, is_same_v

namespace vaspml::text
{

class BoolFlavor;

/**************************************************************************************************
 * Determine type (Int, Real, bool, String) from input string (ignores leading/trailing
 * whitespaces, requires BoolFlavor choice).
 *
 * @param input Input string to check for type.
 * @return Type index of detected type.
 **************************************************************************************************/
template<typename T>
ItemIndex detectType(const String& input)
{
    static_assert(std::is_base_of<BoolFlavor, T>::value, "T must inherit from BoolFlavor");

    if (input.empty()) return ItemIndex::STRING;

    String trimmed = string_tools::trim(input);
    if (trimmed.empty()) return ItemIndex::STRING;

    Size pos;

    try
    {
        std::stoi(trimmed, &pos);
        if (pos != trimmed.size()) throw std::invalid_argument("trailing characters");
        return ItemIndex::INT;
    }
    catch (...)
    {}

    try
    {
        std::stod(trimmed, &pos);
        if (pos != trimmed.size()) throw std::invalid_argument("trailing characters");
        return ItemIndex::REAL;
    }
    catch (...)
    {}

    if (T().isBool(trimmed)) return ItemIndex::BOOL;

    return ItemIndex::STRING;
}
/**************************************************************************************************
 * Convert string to floating-point or integer number.
 *
 * @param input Input string (outer whitespaces will be ignored).
 * @param forbidRealToIntCast If true, throw if there are extra characters after a successfully
 * decoded int value which could be a decimal point and digits from a floating-point number.
 **************************************************************************************************/
template<typename T,
         typename std::enable_if_t<std::negation_v<std::is_base_of<BoolFlavor, T>>>* = nullptr>
T convert(const String& input, bool forbidRealToIntCast = true)
{
    T      result;
    String trimmed = string_tools::trim(input);

    if constexpr (std::is_same_v<T, Real>)
    {
        try
        {
            if (trimmed.empty()) throw std::invalid_argument("empty input string");
            Size pos;
            result = stod(trimmed, &pos);
            if (pos != trimmed.size()) throw std::invalid_argument("stod, trailing characters");
        }
        catch (std::exception& e)
        {
            throw std::runtime_error("Cannot convert string \"" + trimmed
                                     + "\" to floating-point number (" + flf(VASPML_FLF) + ", "
                                     + e.what() + ").");
        }
    }
    else if constexpr (std::is_same_v<T, Int>)
    {
        try
        {
            if (trimmed.empty()) throw std::invalid_argument("empty input string");
            Size pos;
            result = stoi(trimmed, &pos);
            if (forbidRealToIntCast && pos != trimmed.size())
            {
                throw std::invalid_argument("stoi, trailing characters");
            }
        }
        catch (std::exception& e)
        {
            throw std::runtime_error("Cannot convert string \"" + trimmed + "\" to integer number ("
                                     + flf(VASPML_FLF) + ", " + e.what() + ").");
        }
    }
    else static_assert(alwaysFalse<T>, "Conversion from string not implemented.");

    return result;
}
/**************************************************************************************************
 * Convert string to bool (requires BoolFlavor choice).
 *
 * @param input Input string (outer whitespaces will be ignored).
 **************************************************************************************************/
template<typename T, typename std::enable_if_t<std::is_base_of_v<BoolFlavor, T>>* = nullptr>
bool convert(const String& input)
{
    String trimmed = string_tools::trim(input);

    try
    {
        if (trimmed.empty()) throw std::invalid_argument("empty input string");
        T boolFlavor;
        return boolFlavor.convert(trimmed);
    }
    catch (std::exception& e)
    {
        throw std::runtime_error("Cannot convert string \"" + trimmed + "\" to bool ("
                                 + flf(VASPML_FLF) + ", " + e.what() + ").");
    }
}

} //namespace vaspml::text

#endif
