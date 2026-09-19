#ifndef IO_DETAIL_HPP
#define IO_DETAIL_HPP

#include "BoolFlavor.hpp"
#include "Item.hpp"
#include "Record.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "buffer.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "io.hpp"
#include "meta.hpp"
#include "text.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm>   // for find
#include <iosfwd>      // for istream, ostream
#include <regex>       // for regex_replace, regex
#include <type_traits> // for is_same_v

namespace vaspml::io::detail
{

Size bytesLeftInStream(std::istream& strm);
void checkLine(const String& line, const String& search, const Int n);
Int  readLine(std::istream& input, String& buffer, Int skip = 0, bool expectEOF = false);

template<typename T>
T readMlabScalar(String& line)
{
    // clang-format off
    if      constexpr (std::is_same_v<T, Real  >) return stod(line);
    else if constexpr (std::is_same_v<T, Int   >) return stoi(line);
    else if constexpr (std::is_same_v<T, String>) return string_tools::trim(line);
    // clang-format on
    else
    {
        static_assert(alwaysFalse<T>,
                      "Reading an ML_AB scalar value of this type is not implemented.");
        return T{};
    }
}

template<typename T>
T readMlabVector(std::istream& input, Int& n, bool expectEOF = false)
{
    String line;
    T      result;
    using scalarT = typename T::value_type;
    using namespace string_tools;
    while (readLine(input, line, 0, expectEOF) > 0)
    {
        // Increase line counter.
        n++;
        // Skip empty lines (actually this should not happen).
        if (line.size() < 1) continue;
        // If first character is one of the separator characters we reached the end of the vector.
        if (checkForSubstring("*-=", line.substr(0, 1))) break;
        // Split current line and transform elements into desired type.
        for (auto element : split(trim(line), " "))
        {
            result.push_back(readMlabScalar<scalarT>(element));
        }
    }

    return result;
}

template<typename T>
String mlabItemString(const String& label, const T& value, char sep1, char sep2, Int group = 3)
{
    String result;
    result += String(50, sep1) + '\n';
    result += label + '\n';
    result += String(50, sep2) + '\n';
    if constexpr (std::is_same_v<T, Real>) result += str("%24.16E\n", value);
    else if constexpr (std::is_same_v<T, Int>) result += str("%d\n", value);
    else if constexpr (std::is_same_v<T, String>) result += value + "\n";
    else if constexpr (isVector<T>::value)
    {
        using T1 = typename T::value_type;
        for (Size i = 0; i < value.size(); ++i)
        {
            if constexpr (std::is_same_v<T1, Real>) result += str("%24.16E", value[i]);
            else if constexpr (std::is_same_v<T1, Int>) result += str("%d", value[i]);
            else if constexpr (std::is_same_v<T1, String>) result += value[i];
            else static_assert(alwaysFalse<T>, "Writing provided vector type is not implemented.");
            if ((i + 1) % group == 0 || i + 1 == value.size()) result += '\n';
            else result += " ";
        }
    }
    else static_assert(alwaysFalse<T>, "Writing provided type is not implemented.");

    return result;
}

void processMlffHeader(Buffer& buffer, Size& position, bool doWrite, Record& rec);

template<typename T>
void processMlffEntry(Buffer& buffer,
                      Size&   position,
                      bool    doWrite,
                      Record& record,
                      String  key,
                      Int     nchars = -1,
                      bool    assumeEntryPresent = false,
                      Int     createVecWithSize = -1)
{
    if (!doWrite)
    {
        if constexpr (isScalar<T>::value) record.add(key, itemIndex<T>(), assumeEntryPresent);
        else if constexpr (isVector<T>::value)
        {
            if (createVecWithSize >= 0) record.put<T>(key, T(createVecWithSize));
            else record.add(key, itemIndex<T>(), assumeEntryPresent);
        }
    }

    // Now read or write this entry.
    process(buffer, record.get<T>(key), position, doWrite, nchars);

    return;
}

void setupKnownTags(Record& record);
void readIncarToStrings(Record& record, const String& input);
void processIncarTags(Record& record);
void addVaspInterfaceDataToSetup(Record& record,
                                 Int     ntyp,
                                 Int     nsw = 1,
                                 Int     ibrion = 0,
                                 bool    mlabExists = false,
                                 bool    mlffExists = false);

template<typename T>
T convertIncar(const String& strValue, const String& tag)
{
    if constexpr (std::is_same_v<T, Real>) return text::convert<Real>(strValue);
    else if constexpr (std::is_same_v<T, Int>)
    {
        ItemIndex detectedType = text::detectType<text::BoolFlavorIncar>(strValue);
        if (detectedType == ItemIndex::REAL)
        {
            global::tutor.warning("Expected an integer value for tag \"" + tag
                                  + "\", but a floating-point number was provided (\"" + strValue
                                  + "\"). Please check your input!");
        }
        return text::convert<Int>(strValue, false);
    }
    else if constexpr (std::is_same_v<T, String>) return strValue;
    else if constexpr (std::is_same_v<T, bool>)
    {
        return text::convert<text::BoolFlavorIncar>(strValue);
    }
    else if constexpr (isVector<T>::value)
    {
        using T1 = typename T::value_type;
        String tmp = strValue;
        while (true)
        {
            // Check if Fortran repeat specifier is present (e.g. 3*0.0).
            String::const_iterator pos = find(tmp, tmp.begin(), "*", flf(VASPML_FLF), false);
            // If not found, quit loop.
            if (pos == tmp.end()) break;
            // Search for previous separator (Fortran allows spaces and commas).
            String::const_iterator posStart = rfind(tmp, pos, ", ", flf(VASPML_FLF), false);
            // If not found, this must be the first item.
            if (posStart == tmp.end()) posStart = tmp.begin();
            // Else increment iterator to point to number after separator.
            else posStart++;
            // Search for next separator.
            String::const_iterator posEnd = find(tmp, pos, ", ", flf(VASPML_FLF), false);
            // Split into repeat specifier and actual value.
            Vec1String parts = string_tools::split(String{posStart, posEnd}, "*");
            if (parts.size() != 2)
            {
                global::tutor.error("Cannot process vector with repeat specifier, got \""
                                    + String{posStart, posEnd} + "\".");
            }
            Int repeat = std::stoi(parts[0]);
            if (repeat <= 0)
            {
                global::tutor.error("Negative repeat specifier encountered ("
                                    + std::to_string(repeat) + ").");
            }
            String r = parts[1];
            // Before parsing replace "a*b" with a repetitions of "b".
            for (Int i = 2; i <= repeat; ++i) r += " " + parts[1];
            tmp.replace(posStart, posEnd, r);
            // Replace all commas with spaces for splitting below.
            tmp = std::regex_replace(tmp, std::regex(","), " ");
        }
        T result;
        for (auto t : string_tools::split(tmp, " "))
        {
            if constexpr (std::is_same_v<T1, Real>) result.push_back(convertIncar<Real>(t, tag));
            else if constexpr (std::is_same_v<T1, Int>) result.push_back(convertIncar<Int>(t, tag));
            else
            {
                static_assert(alwaysFalse<T>,
                              "INCAR value conversion not implemented for this vector type.");
            }
        }
        return result;
    }
    else
    {
        static_assert(alwaysFalse<T>, "INCAR value conversion not implemented for this type.");
        return T{};
    }
}

template<typename T>
void setTagValue(Record& setup, const Record& incar, String name, T def)
{
    const Record&     types = incar.dcget<ShRec>("types");
    const Record&     tags = incar.dcget<ShRec>("tags");
    const Vec1String& alt = incar.cget<Vec1String>("tag_state_alt");

    if (itemIndex<T>() != SmartEnum<ItemIndex>::toEnum(types.cget<String>(name)))
    {
        global::tutor.bug("Inconsistent types for tag \"" + name + "\" encountered.");
    }
    setup.add(name, itemIndex<T>());
    if (tags.contains(name))
    {
        setup.get<T>(name) = tags.cget<T>(name);
        Int& state = setup.dget<ShRec>("states").get<Int>(name);
        if (std::find(alt.begin(), alt.end(), name) == alt.end())
        {
            state = static_cast<Int>(TagState::INCAR);
        }
        else { state = static_cast<Int>(TagState::INCAR_ALT); }
    }
    else
    {
        setup.get<T>(name) = def;
        setup.dget<ShRec>("states").get<Int>(name) = static_cast<Int>(TagState::DEFAULT);
    }

    if constexpr (std::is_same_v<T, String>)
    {
        const Vec1String& list = constants::caseSensitiveTags;
        if (std::find(list.begin(), list.end(), name) == list.end())
        {
            setup.get<T>(name) = string_tools::makeLowerCase(setup.cget<T>(name));
        }
    }

    return;
}

void   writeMllogfileKernel(const Record& setup, std::ostream& output);
void   writeMllogfileGrace(const Record& setup, std::ostream& output);
String mllogfileSection(char fill = '-', String title = "");
String mllogfileTag(String        tag,
                    const Record& setup,
                    String        vectorCommon = "",
                    Vec1String    vectorIndex = Vec1String{});

template<typename T>
String mllogfileTagManual(String   description,
                          T        value,
                          TagState state = TagState::DEFAULT,
                          String   tag = "")
{
    String result = description + " : ";
    if constexpr (std::is_same_v<T, Real>) result += str("%13.5E", value);
    else if constexpr (std::is_same_v<T, Int>) result += str("%13d", value);
    else
    {
        static_assert(
            alwaysFalse<T>,
            "Writing manual INCAR tag parameter section entry is not imlemented for this type.");
    }
    if (state == TagState::DEFAULT && tag.empty()) return result + "\n";
    tag.resize(24, ' ');
    result += " " + SmartEnum<TagState>::toString(state) + " " + tag + "\n";

    return result;
}

} // namespace vaspml::io::detail

#endif
