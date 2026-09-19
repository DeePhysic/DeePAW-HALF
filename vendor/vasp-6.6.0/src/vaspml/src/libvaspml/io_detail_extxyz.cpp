#include "io_detail_extxyz.hpp"

#include "BoolFlavor.hpp"
#include "Item.hpp"
#include "Record.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "debug.hpp"

#include <algorithm> // for find, count_if, max
#include <regex>     // for regex_replace, regex
#include <stdexcept> // for runtime_error

using namespace vaspml;
using namespace vaspml::io::detail::extxyz;

template<>
SmartEnum<CommentLineMode>::EnumMap SmartEnum<CommentLineMode>::mapEnumsNames()
{
    return SmartEnum<CommentLineMode>::EnumMap{
        {CommentLineMode::Unknown,      "Unknown"     },
        {CommentLineMode::BareString,   "BareString"  },
        {CommentLineMode::NewArray,     "NewArray"    },
        {CommentLineMode::LegacyBraces, "LegacyBraces"},
        {CommentLineMode::LegacyQuotes, "LegacyQuotes"}
    };
}

template<>
SmartEnum<ValueType>::EnumMap SmartEnum<ValueType>::mapEnumsNames()
{
    return SmartEnum<ValueType>::EnumMap{
        {ValueType::Unknown,      "Unknown"     },
        {ValueType::Integer,      "Integer"     },
        {ValueType::Float,        "Float"       },
        {ValueType::Bool,         "Bool"        },
        {ValueType::BareString,   "BareString"  },
        {ValueType::QuotedString, "QuotedString"}
    };
}

enum class LastStop
{
    Open,
    Close,
    Delimiter,
    QuotedString
};

ValueType io::detail::extxyz::valueType(const String& input)
{
    ItemIndex i = text::detectType<text::BoolFlavorExtXyz>(input);

    if (i == ItemIndex::INT) return ValueType::Integer;
    else if (i == ItemIndex::REAL) return ValueType::Float;
    else if (i == ItemIndex::BOOL) return ValueType::Bool;
    else if (i == ItemIndex::STRING)
    {
        if (!input.empty() && input.front() == '"' && input.back() == '"')
        {
            return ValueType::QuotedString;
        }
        else return ValueType::BareString;
    }
    else
    {
        throw std::runtime_error(
            "Unexpected type \"" + toString(i)
            + "\" encountered when converting ItemIndex to extended xyz value type.");
    }
}

void io::detail::extxyz::checkBareString(const String& input)
{
    Size count = countIllegalCharsInBareString(input);
    if (count > 0)
    {
        throw std::runtime_error(
            "Bare string '" + input + "' in extended xyz comment line contains "
            + std::to_string(count)
            + " character(s) which require double quotes around string (any whitespace + any of "
            + illegalCharsInBareString() + ").");
    }

    return;
}

String io::detail::extxyz::parseBareString(const String&           input,
                                           String::const_iterator& pos,
                                           String                  extraSearchChars)
{
    String::const_iterator pos2 = pos;
    if (pos2 != input.end()) ++pos2;

    pos2 = find(input, pos2, string_tools::WHITESPACE + extraSearchChars, flf(VASPML_FLF), false);

    String result(pos, pos2);

    checkBareString(result);

    pos = pos2;

    return result;
}

Size io::detail::extxyz::countIllegalCharsInBareString(const String& input)
{
    String illegalChars = illegalCharsInBareString() + string_tools::WHITESPACE;

    return std::count_if(
        input.begin(),
        input.end(),
        [&illegalChars](char c)
        { return std::find(illegalChars.begin(), illegalChars.end(), c) != illegalChars.end(); });
}

String io::detail::extxyz::removeMeaninglessEscapes(const String& input)
{
    String result = input;

    // Search for "standalone" backslashes, i.e., which are not part of \\ (escaped backslash)
    // or \n (newline representation) or \" (escaped double quote), and remove them.
    String::iterator it = result.begin();
    Int              n = 0;
    while (it != result.end())
    {
        // If backslash found, increment counter, go to next character.
        if (*it == '\\')
        {
            ++n;
            ++it;
        }
        else
        {
            // Something other than backslash found after a sequence of backslashes.
            if (n != 0)
            {
                // If there is an odd number of backslashes and current character is 'n' or '"'
                // do nothing, otherwise, remove the last backslash.
                if (n % 2 == 1 && *it != '"' && *it != 'n') it = result.erase(it - 1);
                else ++it;
                // Reset backslash counter.
                n = 0;
            }
            // Just go to next character if there are no preceeding backslashes.
            else ++it;
        }
    }
    // If there is an odd number of backslashes at the end, remove the last one.
    if (n % 2 == 1) result.pop_back();

    return result;
}

String io::detail::extxyz::removeAllEscapes(const String& input)
{
    String result = input;

    // First, remove escapes which do not actually escape anything.
    result = removeMeaninglessEscapes(result);
    // At this point all remaining slashes are part of valid escape sequences: \", \n or \\.
    // Hence, it is safe to replace now \" and \n with consecutive regex_replace() calls.
    result = std::regex_replace(result, std::regex("\\\\\""), "\"");
    result = std::regex_replace(result, std::regex("\\\\n"), "\n");
    // Finally, replace escaped backslashes with just single backslash.
    // Note: do not put this in front of the \" and \n replacements!
    result = std::regex_replace(result, std::regex("\\\\\\\\"), "\\");

    return result;
}

void io::detail::extxyz::parseCommentLine(Record& record, String input)
{
    using namespace io::detail;

    input = string_tools::trim(input);

    String::const_iterator pos = input.begin();

    while (pos != input.end())
    {
        String key;
        // Expect key, can be
        // * quoted string, search for " quote
        // * bare string, search for letter
        // both could be with or without actual value (assume boolean with "true" value).
        pos = scan(input, pos, flf(VASPML_FLF), "", false);
        // If at end of string, exit loop.
        if (pos == input.end()) break;
        // Parse quoted or simple string.
        if (*pos == '"')
        {
            key = string_tools::parseQuotedString(input, pos);
            key = removeAllEscapes(key);
        }
        else key = parseBareString(input, pos, "=");
        if (key.empty())
        {
            global::tutor.warning("Empty key encountered while parsing extended xyz comment line.");
        }
        // Now at end of key string, search for equal separator.
        pos = scan(input, pos, flf(VASPML_FLF), "", false);
        // Instead of "= value" there could be a key without a value (implicit boolean).
        if (pos == input.end() || *pos != '=')
        {
            record.put(key, true);
            continue;
        }
        // There is an equal sign, hence we expect a value now.
        // Skip equal sign and search for next non-whitespace character.
        pos++;
        pos = scan(input, pos, flf(VASPML_FLF), "", false);
        if (pos == input.end())
        {
            global::tutor.error("Unexpected end of string while parsing value for key \"" + key
                                + "\" in extended xyz comment line.");
        }
        CommentLineValue value = parseCommentLineValue(input, pos);
        interpretCommentLineValue(value);
        convertCommentLineValue(record, key, value);
    }

    return;
}

CommentLineValue io::detail::extxyz::parseCommentLineValue(const String&           input,
                                                           String::const_iterator& pos)
{
    CommentLineValue result;

    // List of possible value types:
    //
    // Primitive types
    // ---------------
    // * integer
    // * float
    // * boolean
    // * bare string
    // * quoted string
    //
    // Array types
    // -----------
    // * new-style 1D or 2D array:
    //   - opens/closes with brackets [], [[],[],...,[]]
    //   - all primitive types allowed
    //   - comma separator
    //   - optional whitespaces
    // * legacy 1D array:
    //   - opens/closes with quotes "" or braces {}
    //   - quotes: all primitive types allowed, except strings (also bare strings)
    //   - braces: all primitive types allowed
    //   - whitespace separator
    //   - single element is converted to scalar

    if (*pos == '"')
    {
        // Store entire string, need to decide later whether this is a quoted string or an array of
        // non-string data.
        result.mode = CommentLineMode::LegacyQuotes;
        result.dim = -1;
        // Surrounding quotes are not included here yet.
        result.full = string_tools::parseQuotedString(input, pos);
        result.parts.resize(1);
        Vec1String& value = result.parts.back();
        value = string_tools::split(result.full, string_tools::WHITESPACE, true, true);
        if (value.empty()) value.push_back(String{});
        result.full = '"' + result.full + '"';
    }
    else if (*pos == '{') { result = parseCommentLineLegacyArray(input, pos); }
    else if (*pos == '[') { result = parseCommentLineNewArray(input, pos, "[]", ","); }
    // Anything else must be a simple string.
    else
    {
        result.mode = CommentLineMode::BareString;
        result.dim = 0;
        result.full = parseBareString(input, pos);
        result.parts.resize(1);
        result.parts.back().push_back(result.full);
    }

    return result;
}

CommentLineValue io::detail::extxyz::parseCommentLineNewArray(const String&           input,
                                                              String::const_iterator& pos,
                                                              String                  enclosures,
                                                              String                  delimiters)
{
    CommentLineValue result;
    result.mode = CommentLineMode::NewArray;
    Int&                         dim = result.dim;
    Vec2String&                  parts = result.parts;
    const String::const_iterator startPos = pos;

    LastStop last = LastStop::Open;
    dim = -1;
    parts.push_back(Vec1String());
    pos++;
    String::const_iterator pos2 = pos;
    String::const_iterator pos3;
    do {
        pos3 = find(input, pos2, enclosures + delimiters + "\"", flf(VASPML_FLF));
        // Opening (,[,{.
        if (*pos3 == enclosures[0])
        {
            if (dim == 1)
            {
                throw std::runtime_error(
                    "Mixed dimensionality found while parsing array in comment line of extended "
                    "xyz file (switch from 1D to 2D).");
            }
            if ((dim == -1 || last == LastStop::Delimiter)
                && !string_tools::trim(String(pos2, pos3)).empty())
            {
                throw std::runtime_error(
                    "Non-whitespace characters found before opening enclosure while parsing array "
                    "in comment line of extended xyz file.");
            }
            if (last == LastStop::Close)
            {
                throw std::runtime_error("Malformatted 2D array found while parsing extended xyz "
                                         "file: missing comma before opening enclosure.");
            }
            dim = 2;
            pos = pos3 + 1;
            pos2 = pos;
            last = LastStop::Open;
            do {
                pos3 = find(input, pos2, enclosures[1] + delimiters + "\"", flf(VASPML_FLF));
                if (*pos3 == '"')
                {
                    if (last == LastStop::QuotedString)
                    {
                        throw std::runtime_error(
                            "Malformatted 2D array content while parsing extended xyz file: "
                            "back-to-back quoted strings without delimiter encountered.");
                    }
                    if (!string_tools::trim(String(pos2, pos3)).empty())
                    {
                        throw std::runtime_error(
                            "Non-whitespace characters found before start of quoted string while "
                            "parsing array in comment line of extended xyz file.");
                    }
                    string_tools::parseQuotedString(input, pos3);
                    pos2 = pos3;
                    last = LastStop::QuotedString;
                }
                else
                {
                    bool   doBreak = *pos3 == enclosures[1];
                    String newItem = string_tools::trim(String(pos, pos3));
                    if (last == LastStop::QuotedString && newItem.back() != '"')
                    {
                        throw std::runtime_error(
                            "Malformatted 2D array content while parsing extended xyz file: "
                            "additional content after closing double quote detected.");
                    }
                    parts.back().push_back(newItem);
                    pos2 = pos3 + 1;
                    pos = pos2;
                    if (doBreak)
                    {
                        last = LastStop::Close;
                        break;
                    }
                    else last = LastStop::Delimiter;
                }
            }
            while (true);
        }
        // Any of the delimiters.
        else if (std::find(delimiters.begin(), delimiters.end(), *pos3) != delimiters.end())
        {
            if (dim == -1) dim = 1;
            if (dim == 2)
            {
                if (last == LastStop::Delimiter)
                {
                    throw std::runtime_error(
                        "Malformatted 2D array while parsing extended xyz file: expected opening "
                        "enclosure after delimiter, got another delimiter.");
                }
                parts.push_back(Vec1String());
            }
            else
            {
                String newItem = string_tools::trim(String(pos, pos3));
                if (last == LastStop::QuotedString && newItem.back() != '"')
                {
                    throw std::runtime_error(
                        "Malformatted 1D array content while parsing extended xyz file: additional "
                        "content after closing double quote detected (before delimiter).");
                }
                parts.back().push_back(newItem);
            }
            pos2 = pos3 + 1;
            pos = pos2;
            last = LastStop::Delimiter;
        }
        // Closing ),],}.
        else if (*pos3 == enclosures[1])
        {
            if (dim == -1) dim = 1;
            if (dim == 1)
            {
                String newItem = string_tools::trim(String(pos, pos3));
                if (last == LastStop::QuotedString && newItem.back() != '"')
                {
                    throw std::runtime_error(
                        "Malformatted 1D array content while parsing extended xyz file: additional "
                        "content after closing double quote detected (before closing enclosure).");
                }
                parts.back().push_back(newItem);
            }
            else if (dim == 2 && last == LastStop::Delimiter)
            {
                throw std::runtime_error(
                    "Malformatted 2D array while parsing extended xyz file: expected opening "
                    "enclosure after delimiter, got closing enclosure.");
            }
            last = LastStop::Close;
            break;
        }
        // Double quote \" => skip quoted string.
        else
        {
            if (dim == -1) dim = 1;
            if (dim == 2)
            {
                throw std::runtime_error(
                    "Mixed dimensionality found while parsing array in comment line of extended "
                    "xyz file (switch from 2D to 1D, got quoted string).");
            }
            if (last == LastStop::QuotedString)
            {
                throw std::runtime_error(
                    "Malformatted 1D array content while parsing extended xyz file: "
                    "back-to-back quoted strings without delimiter encountered.");
            }
            if (!string_tools::trim(String(pos2, pos3)).empty())
            {
                throw std::runtime_error(
                    "Non-whitespace characters found before start of quoted string while parsing "
                    "array in comment line of extended xyz file.");
            }
            string_tools::parseQuotedString(input, pos3);
            pos2 = pos3;
            last = LastStop::QuotedString;
        }
    }
    while (true);
    pos = pos3 + 1;
    result.full = String(startPos, pos);

    return result;
}

CommentLineValue io::detail::extxyz::parseCommentLineLegacyArray(const String&           input,
                                                                 String::const_iterator& pos)
{
    CommentLineValue result;
    result.mode = CommentLineMode::LegacyBraces;
    result.parts.push_back(Vec1String());
    Vec1String&                  vec = result.parts.back();
    const String::const_iterator startPos = pos;

    LastStop last = LastStop::Open;
    pos++;
    String::const_iterator pos2 = pos;
    String::const_iterator pos3;
    do {
        pos3 = find(input, pos2, "}\"" + string_tools::WHITESPACE, flf(VASPML_FLF));
        // Double quote \" => skip quoted string.
        if (*pos3 == '"')
        {
            if (last == LastStop::QuotedString)
            {
                throw std::runtime_error(
                    "Malformatted legacy array content while parsing extended xyz file: "
                    "back-to-back quoted strings without delimiter encountered.");
            }
            if (!String(pos2, pos3).empty())
            {
                throw std::runtime_error(
                    "Non-whitespace characters found before start of quoted string while "
                    "parsing legacy array in comment line of extended xyz file.");
            }
            string_tools::parseQuotedString(input, pos3);
            pos2 = pos3;
            last = LastStop::QuotedString;
        }
        // Closing brace or any of the delimiters.
        else
        {
            String newItem = string_tools::trim(String(pos, pos3));
            if (last == LastStop::QuotedString && newItem.back() != '"')
            {
                throw std::runtime_error(
                    "Malformatted legacy array content while parsing extended xyz file: additional "
                    "content after closing double quote detected.");
            }
            vec.push_back(newItem);
            if (*pos3 == '}')
            {
                last = LastStop::Close;
                break;
            }
            else last = LastStop::Delimiter;
            pos2 = pos3 + 1;
            pos = pos2;
        }
    }
    while (true);
    pos = pos3 + 1;
    result.full = String(startPos, pos);

    // Remove empty entries but leave at least a single empty string.
    vector_tools::removeEmptyElements(vec);
    if (vec.empty()) vec.push_back(String{});

    if (vec.size() == 1) result.dim = 0;
    else result.dim = 1;

    return result;
}

void io::detail::extxyz::interpretCommentLineValue(CommentLineValue& value)
{
    using VT = ValueType;
    using CLM = CommentLineMode;
    const Vec2String& parts = value.parts;

    VT&  type = value.type;
    bool hasIntOrFloat = false;
    for (const auto& p1 : parts)
    {
        for (const auto& p2 : p1)
        {
            VT currentType = valueType(p2);
            if (currentType == VT::Integer || currentType == VT::Float) hasIntOrFloat = true;
            type = std::max(type, currentType);
            if (value.mode != CLM::LegacyQuotes && type == VT::BareString) checkBareString(p2);
            if (type == VT::QuotedString) break;
        }
        if (type == VT::QuotedString) break;
    }

    // Mixing bools with integers or floats makes little sense.
    if (hasIntOrFloat && type == VT::Bool)
    {
        throw std::runtime_error(
            "Detected unsupported array containing a mixture of booleans and integers and/or "
            "floating-point numbers in extended xyz comment line.");
    }

    // Only legacy arrays surrounded by quotes still require determination of actual dimensionality.
    if (value.mode == CLM::LegacyQuotes)
    {
        // Only if all tokens inside legacy quote array are numbers it is actually a one-dimensional
        // array. Otherwise, it is just a quoted string.
        if (type == VT::BareString) value.dim = 0;
        // Single element legacy arrays are interpreted as a scalar value.
        else if (parts.size() == 1 && parts[0].size() == 1) value.dim = 0;
        else value.dim = 1;
        // Collapse parts into single quoted string.
        if (value.dim == 0)
        {
            value.type = VT::QuotedString;
            value.parts.resize(1);
            value.parts.front() = {value.full};
        }
    }

    if (type == VT::Unknown)
    {
        throw std::runtime_error("Could not determine value type while parsing extended xyz file.");
    }

    return;
}

void io::detail::extxyz::convertCommentLineValue(Record&                 record,
                                                 const String&           key,
                                                 CommentLineValue const& value)
{
    using VT = ValueType;
    const VT& type = value.type;

    VASPML_DEBUG_L3(
        if (type == VT::Bool && value.dim > 0)
        {
            global::tutor.warning(
                "Converting boolean vectors is not yet implemented, storing data from extended xyz "
                "comment line as string vector.");
        }
    );

    if (value.dim == 0)
    {
        const String& s = value.parts.front().front();
        if (type == VT::Integer) record.put(key, text::convert<Int>(s));
        else if (type == VT::Float) record.put(key, text::convert<Real>(s));
        else if (type == VT::Bool) record.put(key, text::convert<text::BoolFlavorExtXyz>(s));
        else if (type == VT::BareString) record.put(key, s);
        else if (type == VT::QuotedString) record.put(key, removeAllEscapesAndQuotes(s));
    }
    else if (value.dim == 1)
    {
        const Vec1String& v1 = value.parts.front();
        if (type == VT::Integer) record.put(key, stringsToValues<Vec1Int>(v1));
        else if (type == VT::Float) record.put(key, stringsToValues<Vec1Real>(v1));
        else if (type == VT::Bool) record.put(key, v1);
        else if (type == VT::BareString) record.put(key, v1);
        else if (type == VT::QuotedString) record.put(key, removeAllEscapesAndQuotes(v1));
    }
    else if (value.dim == 2)
    {
        const Vec2String& v2 = value.parts;
        if (type == VT::Integer) record.put(key, stringsToValues<Vec2Int>(v2));
        else if (type == VT::Float) record.put(key, stringsToValues<Vec2Real>(v2));
        else if (type == VT::Bool) record.put(key, v2);
        else if (type == VT::BareString) record.put(key, v2);
        else if (type == VT::QuotedString) record.put(key, removeAllEscapesAndQuotes(v2));
    }

    return;
}

void io::detail::extxyz::parseProperties(Record& record)
{
    if (!record.contains("Properties"))
    {
        record.add("Properties", "ShRec");
        record.get<ShRec>("Properties")->put("raw", String("species:S:1:pos:R:3"));
    }

    return;
}
