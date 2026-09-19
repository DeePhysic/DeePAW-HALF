#ifndef JSONER_HPP
#define JSONER_HPP

#include "Item.hpp"
#include "RecordFunctor.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "meta.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <ostream>     // for basic_ostream, operator<<, ostream
#include <regex>       // for regex, regex_replace
#include <type_traits> // for is_same_v, enable_if_t

namespace vaspml
{
struct Record;
}

namespace vaspml::rec::detail
{

class Jsoner : public RecordFunctor
{
  public:
    /**********************************************************************************************
     * Constructor.
     *
     * @param output Existing and opened output stream reference.
     * @param indentSpaces Specifies indentation level given in number of spaces. Special meaning of
     * negative numbers (see below).
     * @param hideMeta If true, metadata, i.e. all entries with keys starting with "_" will not be
     * written.
     *
     * If `indentSpaces` is negative the formatting will be special:
     *
     * * `indentSpaces == -1`: Removes all newlines but keeps some spaces for readability.
     * * `indentSpaces == -2`: Removes all newlines and also most spaces (real numbers with
     * scientific notation may have some leading spaces remaining).
     **********************************************************************************************/
    Jsoner(std::ostream& output, Int indentSpaces, bool hideMeta);

    template<typename T>
    void writeJsonScalar(T item)
    {
        // clang-format off
        if      constexpr (std::is_same_v<T, Real  >) output << str("%24.16E", item);
        else if constexpr (std::is_same_v<T, Int   >) output << str("%d", item);
        else if constexpr (std::is_same_v<T, String>)
        {
            item = std::regex_replace(item, std::regex("\\\\"), "\\\\");
            item = std::regex_replace(item, std::regex("\""), "\\\"");
            item = std::regex_replace(item, std::regex("[\b]"), "\\b");
            item = std::regex_replace(item, std::regex("\\f"), "\\f");
            item = std::regex_replace(item, std::regex("\\n"), "\\n");
            item = std::regex_replace(item, std::regex("\\r"), "\\r");
            item = std::regex_replace(item, std::regex("\\t"), "\\t");
            output << '"' + item + '"';
        }
        else if constexpr (std::is_same_v<T, bool  >) output << (item ? "true" : "false");
        else output << String{"null"};
        // clang-format on

        return;
    }

    template<typename T>
    void writeJsonVector(T vec)
    {
        // Early exit if vector empty.
        if (vec.empty())
        {
            output << '[' + sp + ']';
            return;
        }

        output << '[' + nl;
        using T1 = typename T::value_type;
        auto last = vec.cend() - 1;
        for (auto item = vec.cbegin(); item != vec.cend(); ++item)
        {
            output << indentBase << indentSingle;
            if constexpr (isScalar<T1>::value) writeJsonScalar<T1>(*item);
            else if constexpr (isVector<T1>::value)
            {
                // Remember original indentation.
                String indentBaseOuter = indentBase;
                // Increase indent.
                indentBase += indentSingle;
                writeJsonVector(*item);
                // Restore original indentation.
                indentBase = indentBaseOuter;
            }
            if (item != last) output << ',' + nl;
        }
        output << nl << indentBase << ']';

        return;
    }

    template<typename T>
    void add(const String& key, const T& item)
    {
        // Do not write entries with metadata.
        if (hideMeta && key[0] == '_') return;

        //output << indentBase + '"' + key + '"' + sp + ':' + sp;
        output << indentBase;
        writeJsonScalar(key);
        output << sp + ':' + sp;

        if constexpr (isScalar<T>::value) writeJsonScalar(item);
        else if constexpr (isVector<T>::value) writeJsonVector(item);
        else static_assert(alwaysFalse<T>, "Type is not supported by Jsoner.");

        return;
    }

    // clang-format off
    void operator()(const Record& record, bool entry) override;
    void operator()(const String& key, const Int&        item) override;
    void operator()(const String& key, const Real&       item) override;
    void operator()(const String& key, const String&     item) override;
    void operator()(const String& key, const bool&       item) override;
    void operator()(const String& key, const Vec1Real&   item) override;
    void operator()(const String& key, const Vec1Int&    item) override;
    void operator()(const String& key, const Vec1String& item) override;
    void operator()(const String& key, const Vec2Real&   item) override;
    void operator()(const String& key, const Vec2Int&    item) override;
    void operator()(const String& key, const Vec2String& item) override;
    // clang-format on

    void down(const String& key, Int arrayIndex = -1, Size arraySize = 1) override;
    void up(const String& key, Int arrayIndex = -1, Size arraySize = 1) override;
    void next(const String& key, bool last) override;
    bool shouldLoop(const String& key, Size arraySize, bool skipRecordArrays) override;
    void postLoop() override;

    std::ostream& output;

  private:
    Int    indentSpaces;
    String indentSingle;
    String indentBase;
    String sp;
    String nl;
    bool   hideMeta;
};

void dejson(const String& input, Record& record, String::const_iterator& pos);
/// Determine type of next content.
ItemIndex dejsonGetType(const String& input, String::const_iterator pos, Int& level);
/**************************************************************************************************
 * Read one scalar value from JSON string.
 *
 * @param delimiters Set of expected characters which terminate the value.
 **************************************************************************************************/
template<typename T, typename std::enable_if_t<isScalar<T>::value>* = nullptr>
T dejsonRead(const String& input, String::const_iterator& pos, String delimiters = ",}")
{
    T result = T();

    if constexpr (std::is_same_v<T, String>)
    {
        pos = scan(input, pos, flf(VASPML_FLF), "\"");
        result = string_tools::parseQuotedString(input, pos);
        // Process escaped characters.
        result = std::regex_replace(result, std::regex("\\\\\\\\"), "\\");
        result = std::regex_replace(result, std::regex("\\\\\""), "\"");
        // Process escaped forward slash as slash (CAN be escaped according to JSON standard).
        result = std::regex_replace(result, std::regex("\\\\/"), "/");
        result = std::regex_replace(result, std::regex("\\\\b"), "\b");
        result = std::regex_replace(result, std::regex("\\\\f"), "\f");
        result = std::regex_replace(result, std::regex("\\\\n"), "\n");
        result = std::regex_replace(result, std::regex("\\\\r"), "\r");
        result = std::regex_replace(result, std::regex("\\\\t"), "\t");
        pos = scan(input, pos, flf(VASPML_FLF), delimiters);
    }
    else if constexpr (std::is_same_v<T, bool>)
    {
        // Note: accepting also first letter uppercase variants "True" and "False" although this is
        // not according to JSON specs because VASP 6.4.X (ML_FF versions 0.2.0 - 0.2.3) incorrectly
        // used this format.
        pos = scan(input, pos, flf(VASPML_FLF), "tTfF");
        String::const_iterator pos2 = find(input, pos + 1, delimiters, flf(VASPML_FLF));
        String                 valueString(pos, pos2);
        String                 sub;
        if (sub = valueString.substr(0, 4); sub == "true" || sub == "True") result = true;
        else if (sub = valueString.substr(0, 5); sub == "false" || sub == "False") result = false;
        else
        {
            global::tutor.error(
                "Invalid JSON: Unknown boolean value, expected \"true\" or \"false\" but got \""
                + valueString + "\".");
        }
        pos = pos2;
    }
    else
    {
        String::const_iterator pos2 = find(input, pos + 1, delimiters, flf(VASPML_FLF));
        String                 valueString(pos, pos2);
        if constexpr (std::is_same_v<T, Real>) result = std::stod(valueString);
        if constexpr (std::is_same_v<T, Int>) result = std::stoi(valueString);
        pos = pos2;
    }

    return result;
}

} // namespace vaspml::rec::detail

#endif
