#ifndef META_HPP
#define META_HPP

#include <memory>
#include <type_traits>
#include <vector>

namespace vaspml
{

/*================================================================================================+
 | Metafunction isSharedPtr()
 +================================================================================================*/
/// Metafunction to determine if a type is a shared_ptr, general case (false).
template<typename T>
struct isSharedPtr : std::false_type
{};
/// Specialization if type is actually a shared_ptr (true).
template<typename T>
struct isSharedPtr<std::shared_ptr<T>> : std::true_type
{};

/*================================================================================================+
 | Metafunction isVector()
 +================================================================================================*/
/// Metafunction to determine if a type is a std::vector, general case (false).
template<typename T>
struct isVector : std::false_type
{};
/// Specialization if type is actually a std::vector (true).
template<typename T>
struct isVector<std::vector<T>> : std::true_type
{};

/*================================================================================================+
 | Metafunction isVectorSharedPtr()
 +================================================================================================*/
/// Metafunction to determine if a type is a std::vector of std::shared_ptr, general case (false).
template<typename T>
struct isVectorSharedPtr : std::false_type
{};
/// Specialization if type is actually a std::vector of std::shared_ptr (true).
template<typename T>
struct isVectorSharedPtr<std::vector<std::shared_ptr<T>>> : std::true_type
{};

/*================================================================================================+
 | Metafunction isVector1D()
 +================================================================================================*/
/// Metafunction to determine if a type is strictly a 1D-std::vector, general case (false).
template<typename T>
struct isVector1D : std::false_type
{};
/// Specialization if type is actually a 1D-std::vector (true).
template<typename U, typename Alloc>
struct isVector1D<std::vector<U, Alloc>> : std::integral_constant<bool, !isVector1D<U>::value>
{};

/*================================================================================================+
 | Metafunction isVector2D()
 +================================================================================================*/
/// Metafunction to determine if a type is strictly a 2D-std::vector, general case (false).
template<typename T>
struct isVector2D : std::false_type
{};
/// Metafunction to determine if a type is strictly a 2D-std::vector, 2 arguments case (false).
template<typename T, typename AllocOuter>
struct isVector2D<std::vector<T, AllocOuter>> : std::false_type
{};
/// Specialization if type is actually a 2D-std::vector (true).
template<typename T, typename AllocInner, typename AllocOuter>
struct isVector2D<std::vector<std::vector<T, AllocInner>, AllocOuter>> :
    std::integral_constant<bool, !isVector2D<T>::value && !isVector1D<T>::value>
{};

/**************************************************************************************************
 * Always-false template expression to be used in static_asserts.
 *
 * Use in `if constexpr ...` blocks with `else static_assert(alwaysFalse<T>, "..")` at the end to
 * disallow compilation if the type `T` is not covered.
 **************************************************************************************************/
template<typename T>
constexpr std::false_type alwaysFalse{};

} //namespace vaspml

#endif
