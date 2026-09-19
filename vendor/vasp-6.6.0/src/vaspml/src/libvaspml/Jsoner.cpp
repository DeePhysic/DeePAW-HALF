#include "Jsoner.hpp"

#include "Record.hpp"
#include "SmartEnum.hpp"

#include <algorithm> // for find, find_if_not, count_if
#include <cctype>    // for isdigit

using namespace vaspml;
using namespace vaspml::rec::detail;

Jsoner::Jsoner(std::ostream& output, Int indentSpaces, bool hideMeta) :
    output(output),
    indentSpaces(indentSpaces),
    indentSingle(""),
    hideMeta(hideMeta)
{
    if (indentSpaces >= 0)
    {
        indentSingle = String(indentSpaces, ' ');
        sp = " ";
        nl = "\n";
    }
    else if (indentSpaces == -1)
    {
        sp = " ";
        nl = "";
        indentBase = " ";
    }
    else if (indentSpaces == -2)
    {
        sp = "";
        nl = "";
        indentBase = "";
    }
}

void Jsoner::operator()(const Record& /*record*/, bool entry)
{
    // Top-level opening brace.
    if (level == 0 && entry)
    {
        output << '{' + nl;
        level++;
    }
    // Top-level closing brace.
    else if (level == 1 && !entry)
    {
        level--;
        if (indentSpaces == -1) output << sp;
        output << '}' + nl;
    }
    if (indentSpaces >= 0) indentBase = String(indentSpaces * level, ' ');

    return;
}

// clang-format off
void Jsoner::operator()(const String& key, const Int&        item) { add(key, item); }
void Jsoner::operator()(const String& key, const Real&       item) { add(key, item); }
void Jsoner::operator()(const String& key, const String&     item) { add(key, item); }
void Jsoner::operator()(const String& key, const bool&       item) { add(key, item); }
void Jsoner::operator()(const String& key, const Vec1Real&   item) { add(key, item); }
void Jsoner::operator()(const String& key, const Vec1Int&    item) { add(key, item); }
void Jsoner::operator()(const String& key, const Vec1String& item) { add(key, item); }
void Jsoner::operator()(const String& key, const Vec2Real&   item) { add(key, item); }
void Jsoner::operator()(const String& key, const Vec2Int&    item) { add(key, item); }
void Jsoner::operator()(const String& key, const Vec2String& item) { add(key, item); }
// clang-format on

void Jsoner::down(const String& key, Int arrayIndex, Size /*arraySize*/)
{
    this->key = key;

    // If single entry, write key, colon and open brace.
    if (arrayIndex < 0) output << indentBase + '"' + key + '"' + sp + ':' + sp + '{' + nl;

    // Increase indentation.
    level++;
    if (indentSpaces >= 0) indentBase = String(indentSpaces * level, ' ');

    // For any Record vector entry, write opening brace.
    if (arrayIndex >= 0)
    {
        output << indentBase << '{' + nl;
        // Increase indentation again.
        level++;
        indentBase += indentSingle;
    }

    return;
}

void Jsoner::up(const String& /*key*/, Int arrayIndex, Size arraySize)
{
    // Check if this is the last entry in a Record vector.
    bool lastInArray = arrayIndex >= 0 && (arrayIndex + 1 == static_cast<Int>(arraySize));

    // Decrease indentation.
    level--;
    if (indentSpaces >= 0) indentBase = String(indentSpaces * level, ' ');

    // For single or vector entry, write closing brace.
    output << indentBase << '}';
    // Single entry is completed, return immediately.
    if (arrayIndex < 0) return;

    // If there is another vector entry, write comma and continue in next line.
    if (arrayIndex >= 0 && !lastInArray) output << ',' + nl;
    // If last vector entry, no comma, just newline.
    else if (lastInArray) output << nl;

    // Decrease indentation again for vector entries.
    level--;
    if (indentSpaces >= 0) indentBase = String(indentSpaces * level, ' ');

    return;
}

void Jsoner::next(const String& key, bool last)
{
    if (hideMeta && key[0] == '_') return;

    if (!last) output << ',';
    output << nl;

    return;
}

bool Jsoner::shouldLoop(const String& key, Size /* arraySize */, bool /* skipRecordArrays */)
{
    // Before descending into loop, write key and open brackets.
    output << indentBase + '"' + key + '"' + sp + ':' + sp << '[' + nl;

    return true;
}

void Jsoner::postLoop()
{
    // After finishing the loop, write closing bracket.
    output << indentBase << ']';

    return;
}

void rec::detail::dejson(const String& input, Record& record, String::const_iterator& pos)
{
    pos = scan(input, pos, flf(VASPML_FLF), "{");
    pos++;

    // A JSON object can now start with the first entry or it is empty, hence search for quote or
    // closing brace.
    pos = scan(input, pos, flf(VASPML_FLF), "\"}");
    if (*pos == '}')
    {
        pos++;
        return;
    }

    while (true)
    {
        pos = scan(input, pos, flf(VASPML_FLF), "\"");
        pos++;
        String::const_iterator pos2 = find(input, pos, "\"", flf(VASPML_FLF));
        String                 key(pos, pos2);
        key = string_tools::trim(key);
        if (key == "") global::tutor.error("Invalid JSON: key cannot be empty.");
        pos = pos2 + 1;

        // Search for colon separating key and value.
        pos = scan(input, pos, flf(VASPML_FLF), ":");
        pos++;

        // Find value type and add empty entry in Record.
        Int       level = 0;
        ItemIndex type = dejsonGetType(input, pos, level);
        record.add(key, type);

        if (type == ItemIndex::REAL) record.get<Real>(key) = dejsonRead<Real>(input, pos);
        else if (type == ItemIndex::INT) record.get<Int>(key) = dejsonRead<Int>(input, pos);
        else if (type == ItemIndex::STRING)
        {
            record.get<String>(key) = dejsonRead<String>(input, pos);
        }
        else if (type == ItemIndex::BOOL) record.get<bool>(key) = dejsonRead<bool>(input, pos);
        // For JSON sub-object descend into next level of dejson().
        else if (type == ItemIndex::SHREC) dejson(input, record.dget<ShRec>(key), pos);
        else if (type == ItemIndex::VEC1REAL || type == ItemIndex::VEC1INT
                 || type == ItemIndex::VEC1STRING || type == ItemIndex::VEC1SHREC)
        {
            pos = scan(input, pos, flf(VASPML_FLF), "[");
            do {
                pos++;
                if (type == ItemIndex::VEC1REAL)
                {
                    record.get<Vec1Real>(key).push_back(dejsonRead<Real>(input, pos, ",]"));
                }
                else if (type == ItemIndex::VEC1INT)
                {
                    record.get<Vec1Int>(key).push_back(dejsonRead<Int>(input, pos, ",]"));
                }
                else if (type == ItemIndex::VEC1STRING)
                {
                    record.get<Vec1String>(key).push_back(dejsonRead<String>(input, pos, ",]"));
                }
                // For vector of JSON sub-objects descend into each of them.
                else if (type == ItemIndex::VEC1SHREC)
                {
                    pos = scan(input, pos, flf(VASPML_FLF), "{]");
                    // The vector could be empty, break loop immediately.
                    if (*pos == ']') break;
                    record.get<Vec1ShRec>(key).push_back(std::make_shared<Record>());
                    dejson(input, *(record.get<Vec1ShRec>(key).back()), pos);
                    // Search for comma signalling next item or closing bracket for end of array.
                    pos = scan(input, pos, flf(VASPML_FLF), ",]");
                }
            }
            while (*pos == ',');
            pos++;
        }
        else if (type == ItemIndex::VEC2REAL || type == ItemIndex::VEC2INT
                 || type == ItemIndex::VEC2STRING)
        {
            pos = scan(input, pos, flf(VASPML_FLF), "[");
            do {
                pos++;
                if (type == ItemIndex::VEC2REAL) record.get<Vec2Real>(key).push_back(Vec1Real());
                else if (type == ItemIndex::VEC2INT) record.get<Vec2Int>(key).push_back(Vec1Int());
                else if (type == ItemIndex::VEC2STRING)
                {
                    record.get<Vec2String>(key).push_back(Vec1String());
                }
                pos = scan(input, pos, flf(VASPML_FLF), "[");
                // Sanity check: is there actually any content in this array?
                pos2 = find(input, pos + 1, "]", flf(VASPML_FLF));
                String arrayContent(pos + 1, pos2);
                arrayContent = string_tools::trim(arrayContent);
                // Without content, skip reading of values (would throw).
                if (arrayContent.empty()) pos = pos2;
                else
                {
                    do {
                        pos++;
                        if (type == ItemIndex::VEC2REAL)
                        {
                            record.get<Vec2Real>(key).back().push_back(
                                dejsonRead<Real>(input, pos, ",]"));
                        }
                        else if (type == ItemIndex::VEC2INT)
                        {
                            record.get<Vec2Int>(key).back().push_back(
                                dejsonRead<Int>(input, pos, ",]"));
                        }
                        else if (type == ItemIndex::VEC2STRING)
                        {
                            record.get<Vec2String>(key).back().push_back(
                                dejsonRead<String>(input, pos, ",]"));
                        }
                    }
                    while (*pos == ',');
                }
                pos++;
                pos = scan(input, pos, flf(VASPML_FLF), ",]");
            }
            while (*pos == ',');
            pos++;
        }

        // Search for comma signalling next item or closing brace for end of current Record.
        pos = scan(input, pos, flf(VASPML_FLF), ",}");

        if (*pos == ',') pos++;
        else
        {
            pos++;
            break;
        }
    }

    return;
}

ItemIndex rec::detail::dejsonGetType(const String& input, String::const_iterator pos, Int& level)
{
    // Search for any content following whitespaces.
    pos = std::find_if_not(pos,
                           input.end(),
                           [](const char& c)
                           {
                               const String& s = string_tools::WHITESPACE;
                               return std::find(s.begin(), s.end(), c) != s.end();
                           });
    if (pos == input.end())
    {
        global::tutor.error("Invalid JSON: No content found after colon marker.");
    }

    // If value is array, descend down until value type can be determined.
    if (*pos == '[')
    {
        bool firstCall = level == 0;
        pos++;
        level++;
        ItemIndex valueType = dejsonGetType(input, pos, level);
        // If this is not the root call immediately return the value type found below.
        if (!firstCall) return valueType;
        // Depending on how many levels down the recursion ended, return now the corresponding
        // vector dimension and value type.
        if (level == 1)
        {
            if (valueType == ItemIndex::REAL) return ItemIndex::VEC1REAL;
            if (valueType == ItemIndex::INT) return ItemIndex::VEC1INT;
            if (valueType == ItemIndex::STRING) return ItemIndex::VEC1STRING;
            if (valueType == ItemIndex::SHREC) return ItemIndex::VEC1SHREC;
        }
        else if (level == 2)
        {
            if (valueType == ItemIndex::REAL) return ItemIndex::VEC2REAL;
            if (valueType == ItemIndex::INT) return ItemIndex::VEC2INT;
            if (valueType == ItemIndex::STRING) return ItemIndex::VEC2STRING;
        }
        // If the type was not covered yet, it is not implemented.
        global::tutor.error(str("Cannot read JSON: %d-dimensional vector of \"", level)
                            + toString(valueType) + "\" is not implemented.");
        return ItemIndex::UNDEFINED;
    }
    if (*pos == '{') return ItemIndex::SHREC;
    // Note: accepting also first letter uppercase variants "True" and "False" although this is not
    // according to JSON specs because VASP 6.4.X (ML_FF versions 0.2.0 - 0.2.2) incorrectly used
    // this format.
    if (*pos == 't' || *pos == 'T' || *pos == 'f' || *pos == 'F') return ItemIndex::BOOL;
    if (*pos == '"') return ItemIndex::STRING;
    if (*pos == 'n' && String(pos, pos + 4) == "null")
    {
        global::tutor.error("JSON \"null\" value is not implemented.");
    }

    // At this point value must be integer or floating-point variable.
    // A number may start with a plus or minus sign.
    if (*pos == '+' || *pos == '-') pos++;
    // Search for any non-digit character.
    String::const_iterator old = pos;
    pos = std::find_if_not(pos, input.end(), [](unsigned char c) { return std::isdigit(c); });
    if (pos == input.end())
    {
        global::tutor.error("Invalid JSON: Unexpected end while determining numeric type.");
    }
    // A real number may have a decimal point, could be written in scientific notation (e, E), be
    // infinite (Inf, inf) or not-a-number (NaN, nan).
    if (*pos == '.' || *pos == 'e' || *pos == 'E' || *pos == 'i' || *pos == 'I' || *pos == 'n'
        || *pos == 'N')
    {
        return ItemIndex::REAL;
    }
    // Since only integer is left as value option, check if there is at least one digit.
    if (std::count_if(old, pos, [](unsigned char c) { return std::isdigit(c); }) <= 0)
    {
        // If there are no digits, it could still be an empty array. Then, assume an empty
        // Vec1ShRec.
        if (*pos == ']') return ItemIndex::SHREC;
        else
        {
            global::tutor.error("Invalid JSON: Unknown value encountered (expected digits), "
                                "cannot deduce value type.");
        }
    }

    // Now there is only integer left as option, but perform a final sanity check before returning.
    // A value can only be followed by a comma, closing brace or bracket.
    String sane = string_tools::WHITESPACE + ",}]";
    if (std::find(sane.begin(), sane.end(), *pos) != sane.end()) return ItemIndex::INT;
    else
    {
        global::tutor.error("Invalid JSON: Unexpected character \"" + String(1, *pos)
                            + "\" encountered.");
    }

    return ItemIndex::UNDEFINED;
}
