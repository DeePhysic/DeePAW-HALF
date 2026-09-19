#ifndef SMARTENUM_HPP
#define SMARTENUM_HPP

#include "debug.hpp"
#include "types.hpp"

#include <sstream>     // for basic_ostringstream, ostream, ostringstream
#include <stdexcept>   // for runtime_error
#include <type_traits> // for is_same_v, enable_if_t, false_type, true_type
#include <typeinfo>    // for type_info

namespace vaspml
{

/*================================================================================================+
 | NOTE: Steps to add new SmartEnums here:
 | 1. Add forward declaration right below.
 | 2. Add enum to metafunction isSmart below.
 | 3. Add declaration of static "mapEnumsNames" member function at the end of this file.
 +================================================================================================*/
// Forward declare all enums that should be treated as smart enum.
// Add aditional enums here:
namespace math
{
enum class CutoffType;
}
enum class ItemIndex : int8_t;
enum class DescriptorStorage;
enum class TotalEnergyType;
namespace io
{
enum class TagState : Int;
namespace detail::extxyz
{
enum class CommentLineMode;
enum class ValueType;
} //namespace detail::extxyz
} //namespace io
// Define metafunction to check whether a given enum is a smart enum.
// Only if listed here the smart enum features work.
// First, the general case that defines the negative outcome.
// Add additional enums here with || (OR) combination:
// clang-format off
template<typename T, bool smart =
    std::is_same_v<T, math::CutoffType> ||
    std::is_same_v<T, ItemIndex> ||
    std::is_same_v<T, DescriptorStorage> ||
    std::is_same_v<T, TotalEnergyType> ||
    std::is_same_v<T, io::TagState> ||
    std::is_same_v<T, io::detail::extxyz::CommentLineMode> ||
    std::is_same_v<T, io::detail::extxyz::ValueType>
>
// clang-format on
struct isSmartEnum : std::false_type
{};
// Now, the partial specialization if the above expression is true.
template<typename T>
struct isSmartEnum<T, true> : std::true_type
{};

/** Provides static functions to convert enums to strings and vice-versa.
 *
 * Strings corresponding to enums need to be defined first in a map, see implementation for
 * math::CutoffType. In addition, the metafunction isSmartEnum needs to be updated (requires a
 * forward declaration in the header of SmartEnum. Also, an explicit specialization of the "names"
 * static member needs to be added to the header.
 */
template<typename T>
class SmartEnum
{
    static_assert(isSmartEnum<T>::value, "Template argument is not a smart enum.");

  public:
    typedef std::map<T, String> EnumMap;

    static T toEnum(String stringToConvert)
    {
        EnumMap m = mapEnumsNames();
        for (const auto& n : m)
        {
            if (n.second == stringToConvert) return n.first;
        }
        throw std::runtime_error("ERROR: Provided string \"" + stringToConvert
                                 + "\" does not correspond to an enum.");
    };

    static String toString(T enumToConvert)
    {
        EnumMap m = mapEnumsNames();
        VASPML_DEBUG_L1(
            if (m.find(enumToConvert) == m.end())
            {
                throw std::runtime_error("ERROR: Invalid value ("
                                         + std::to_string(static_cast<Int64>(enumToConvert))
                                         + ") encountered during name lookup for SmartEnum \""
                                         + typeid(T).name() + "\".");
            }
        );
        return m.at(enumToConvert);
    };

    static std::vector<T> listEnums()
    {
        EnumMap        m = mapEnumsNames();
        std::vector<T> keys;
        for (const auto& n : m) keys.push_back(n.first);
        return keys;
    };

    static Vec1String listNames()
    {
        EnumMap    m = mapEnumsNames();
        Vec1String values;
        for (const auto& n : m) values.push_back(n.second);
        return values;
    };

    static EnumMap mapEnumsNames();
};

template<typename T, typename std::enable_if_t<isSmartEnum<T>::value>* = nullptr>
std::ostream& operator<<(std::ostream& os, const T& enumToStream)
{
    os << SmartEnum<T>::toString(enumToStream);
    return os;
}

template<typename T, typename std::enable_if_t<isSmartEnum<T>::value>* = nullptr>
String toString(const T& enumToConvert)
{
    std::ostringstream ss;
    ss << SmartEnum<T>::toString(enumToConvert);
    return ss.str();
}

// Explicit specialization of static member function "mapEnumsNames".
template<>
SmartEnum<ItemIndex>::EnumMap SmartEnum<ItemIndex>::mapEnumsNames();
template<>
SmartEnum<math::CutoffType>::EnumMap SmartEnum<math::CutoffType>::mapEnumsNames();
template<>
SmartEnum<DescriptorStorage>::EnumMap SmartEnum<DescriptorStorage>::mapEnumsNames();
template<>
SmartEnum<TotalEnergyType>::EnumMap SmartEnum<TotalEnergyType>::mapEnumsNames();
template<>
SmartEnum<io::TagState>::EnumMap SmartEnum<io::TagState>::mapEnumsNames();
template<>
SmartEnum<io::detail::extxyz::CommentLineMode>::EnumMap
SmartEnum<io::detail::extxyz::CommentLineMode>::mapEnumsNames();
template<>
SmartEnum<io::detail::extxyz::ValueType>::EnumMap
SmartEnum<io::detail::extxyz::ValueType>::mapEnumsNames();

} //namespace vaspml

#endif
