#ifndef IO_DETAIL_EXTXYZ_HPP
#define IO_DETAIL_EXTXYZ_HPP

#include "meta.hpp"
#include "text.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <type_traits> // for is_same_v

namespace vaspml
{
struct Record;
}

namespace vaspml::io::detail::extxyz
{

enum class CommentLineMode
{
    Unknown,
    BareString,
    NewArray,
    LegacyBraces,
    LegacyQuotes,
};

enum class ValueType
{
    Unknown,
    Integer,
    Float,
    Bool,
    BareString,
    QuotedString
};

/**************************************************************************************************
 * Intermediate structure to hold data retrieved from parsing a comment line value.
 **************************************************************************************************/
struct CommentLineValue
{
    CommentLineMode mode = CommentLineMode::Unknown;
    Int             dim = -1;
    Vec2String      parts;
    String          full;
    ValueType       type = ValueType::Unknown;
};

ValueType     valueType(const String& input);
void          checkBareString(const String& input);
String        parseBareString(const String&           input,
                              String::const_iterator& pos,
                              String                  extraSearchChars = "");
inline String illegalCharsInBareString()
{
    return String("=\",[]{}");
}
Size             countIllegalCharsInBareString(const String& input);
String           removeMeaninglessEscapes(const String& input);
String           removeAllEscapes(const String& input);
void             parseCommentLine(Record& record, String input);
CommentLineValue parseCommentLineValue(const String& input, String::const_iterator& pos);
CommentLineValue parseCommentLineNewArray(const String&           input,
                                          String::const_iterator& pos,
                                          String                  enclosures,
                                          String                  delimiters);
CommentLineValue parseCommentLineLegacyArray(const String& input, String::const_iterator& pos);
void             interpretCommentLineValue(CommentLineValue& value);
void convertCommentLineValue(Record& record, const String& key, const CommentLineValue& value);
void parseProperties(Record& record);

template<typename T>
T readScalar(String& line)
{
    // clang-format off
    if      constexpr (std::is_same_v<T, Real  >) return stod(line);
    else if constexpr (std::is_same_v<T, Int   >) return stoi(line);
    else if constexpr (std::is_same_v<T, String>) return string_tools::trim(line);
    // clang-format on
    else
    {
        static_assert(alwaysFalse<T>,
                      "Reading an extended xyz scalar value of this type is not implemented.");
        return T{};
    }
}

template<typename T>
String valueToString(const T& value, String fmt, bool trim = false)
{
    String result;
    if constexpr (std::is_same_v<T, Real> || std::is_same_v<T, Int>)
    {
        result = str(fmt.c_str(), value);
    }
    else if constexpr (std::is_same_v<T, String>) result = str(fmt.c_str(), value.c_str());
    else if constexpr (isVector<T>::value)
    {
        using T1 = typename T::value_type;
        for (Size i = 0; i < value.size(); ++i)
        {
            if constexpr (std::is_same_v<T1, Real> || std::is_same_v<T1, Int>)
            {
                result += str(fmt.c_str(), value[i]);
            }
            else if constexpr (std::is_same_v<T1, String>)
            {
                result += str(fmt.c_str(), value[i].c_str());
            }
            else static_assert(alwaysFalse<T>, "Writing provided vector type is not implemented.");
            if (i < value.size() - 1) result += " ";
        }
    }
    else static_assert(alwaysFalse<T>, "Writing provided type is not implemented.");

    if (trim) result = string_tools::trim(result);

    return result;
}

template<typename T>
T stringsToValues(const Vec1String& parts)
{
    T result(parts.size());
    using T1 = typename T::value_type;
    for (Size i = 0; i < parts.size(); ++i) result[i] = text::convert<T1>(parts[i]);

    return result;
}

template<typename T>
T stringsToValues(const Vec2String& parts)
{
    T result(parts.size());
    using T1 = typename T::value_type;
    for (Size i = 0; i < parts.size(); ++i) result[i] = stringsToValues<T1>(parts[i]);

    return result;
}

template<typename T>
T removeAllEscapesAndQuotes(const T& input)
{
    T result = input;

    if constexpr (isVector<T>())
    {
        for (auto& i : result) i = removeAllEscapesAndQuotes(i);
    }
    else
    {
        if (result.empty()) return result;
        if (result.back() == '"') result.pop_back();
        if (result.front() == '"') result.erase(result.begin());
        result = removeAllEscapes(result);
    }

    return result;
}

} // namespace vaspml::io::detail::extxyz

#endif
