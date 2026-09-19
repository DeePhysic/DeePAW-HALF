#ifndef ITEM_HPP
#define ITEM_HPP

#include "meta.hpp"
#include "types.hpp"

#include <variant>
#include <type_traits>

namespace vaspml
{

/**************************************************************************************************
 * Variant specialization which can represent all data types typically used in VASPml.
 *
 * If this list gets modified all RecordFunctor and all derived classes must also be updated and
 * checked for compatibility. For more details, refer to Record.
 **************************************************************************************************/
using Item = std::variant<Real,
                          Int,
                          String,
                          bool,
                          ShRec,
                          Vec1Real,
                          Vec1Int,
                          Vec1String,
                          Vec1ShRec,
                          Vec2Real,
                          Vec2Int,
                          Vec2String>;

enum class ItemIndex : int8_t
{
    // clang-format off
    // Actual type indices.
    UNDEFINED    =  0,
    REAL         =  1,
    INT          =  2,
    STRING       =  3,
    BOOL         =  4,
    SHREC        =  5,
    VEC1REAL     =  6,
    VEC1INT      =  7,
    VEC1STRING   =  8,
    VEC1SHREC    =  9,
    VEC2REAL     =  10,
    VEC2INT      =  11,
    VEC2STRING   =  12,
    // Control markers for binary serialization.
    RECORD_BEGIN = -1,
    RECORD_ITEM  = -2,
    RECORD_END   = -3
    // clang-format on
};

/// Metafunction to determine if a type is a scalar type in the Item variant, general case (false).
template<typename T,
         bool scalar =
             // clang-format off
    std::is_same_v<T, Real>   ||
    std::is_same_v<T, Int>    ||
    std::is_same_v<T, String> ||
    std::is_same_v<T, bool>
         // clang-format on
         >
struct isScalar : std::false_type
{};
/// Specialization if type is actually a scalar type in the Item variant (true).
template<typename T>
struct isScalar<T, true> : std::true_type
{};

template<typename T>
constexpr ItemIndex itemIndex()
{
    // clang-format off
    if      constexpr (std::is_same_v<T, Real      >) return ItemIndex::REAL;
    else if constexpr (std::is_same_v<T, Int       >) return ItemIndex::INT;
    else if constexpr (std::is_same_v<T, String    >) return ItemIndex::STRING;
    else if constexpr (std::is_same_v<T, bool      >) return ItemIndex::BOOL;
    else if constexpr (std::is_same_v<T, ShRec     >) return ItemIndex::SHREC;
    else if constexpr (std::is_same_v<T, Vec1Real  >) return ItemIndex::VEC1REAL;
    else if constexpr (std::is_same_v<T, Vec1Int   >) return ItemIndex::VEC1INT;
    else if constexpr (std::is_same_v<T, Vec1String>) return ItemIndex::VEC1STRING;
    else if constexpr (std::is_same_v<T, Vec1ShRec >) return ItemIndex::VEC1SHREC;
    else if constexpr (std::is_same_v<T, Vec2Real  >) return ItemIndex::VEC2REAL;
    else if constexpr (std::is_same_v<T, Vec2Int   >) return ItemIndex::VEC2INT;
    else if constexpr (std::is_same_v<T, Vec2String>) return ItemIndex::VEC2STRING;
    else static_assert(alwaysFalse<T>, "Provided type is not an Index type.");
    // clang-format on

    // Dummy return to silence compiler warnings about missing return statement.
    return ItemIndex::UNDEFINED;
}

} // namespace vaspml

#endif
