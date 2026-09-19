#include "utils.hpp"

#include <cctype>    // for tolower, toupper
#include <chrono>    // for duration, system_clock, time_point_cast, millis...
#include <cstdarg>   // for va_end, va_list, va_start
#include <cstdio>    // for vsnprintf
#include <ctime>     // for tm, localtime, time_t
#include <iostream>  // for cout, cerr
#include <stdexcept> // for runtime_error, logic_error

namespace vaspml
{

String str(const char* format, ...)
{
    char buffer[STR_MAX];

    va_list args;
    va_start(args, format);
    int n = vsnprintf(buffer, STR_MAX, format, args);
    if (n < 0)
    {
        std::cerr << "WARNING: There was an encoding error, maybe some output was lost "
                     "(format string: \""
                  << format << "\").\n";
    }
    // n is the number of characters which should have been written (without null terminator). If n
    // is STR_MAX, already one character is lost.
    else if (static_cast<Size>(n) >= STR_MAX)
    {
        std::cerr << "WARNING: Write buffer was too small, " << std::to_string(n - STR_MAX + 1)
                  << " characters were lost.\n";
    }
    va_end(args);

    // Construct resulting string from buffer contents.
    String result(buffer);

    return result;
}

void now(Vec1Int& values)
{
    using namespace std::chrono;

    const auto  cNow = system_clock::now();
    std::time_t tNow = system_clock::to_time_t(cNow);
    std::tm*    tmNow = std::localtime(&tNow);
    Int         year = tmNow->tm_year + 1900;
    Int         month = tmNow->tm_mon + 1;
    Int         day = tmNow->tm_mday;
    Int         hour = tmNow->tm_hour;
    Int         min = tmNow->tm_min;
    Int         sec = tmNow->tm_sec;
    Int         msec = time_point_cast<milliseconds>(cNow).time_since_epoch().count() % 1000;

    values = Vec1Int{year, month, day, hour, min, sec, msec};

    return;
}

String now()
{
    Vec1Int values;

    now(values);
    Int year = values[0];
    Int month = values[1];
    Int day = values[2];
    Int hour = values[3];
    Int min = values[4];
    Int sec = values[5];
    Int msec = values[6];

    return str("%4d-%02d-%02dT%02d:%02d:%02d.%03d", year, month, day, hour, min, sec, msec);
}

String::const_iterator scan(const String&          input,
                            String::const_iterator pos,
                            String                 sFLF,
                            String                 expect,
                            bool                   errorOnEnd)
{
    // Search for any content following whitespaces.
    pos = std::find_if_not(pos,
                           input.end(),
                           [](const char& c)
                           {
                               const String& s = string_tools::WHITESPACE;
                               return std::find(s.begin(), s.end(), c) != s.end();
                           });
    // Check if unexpected end of string was reached.
    if (pos == input.end())
    {
        if (errorOnEnd)
        {
            throw std::runtime_error(
                "Invalid string: Content ended but expected non-whitespace character (" + sFLF
                + ").");
        }
        else return pos;
    }
    if ((!expect.empty()) && std::find(expect.begin(), expect.end(), *pos) == expect.end())
    {
        throw std::runtime_error("Invalid string: Unexpected character \"" + String(1, *pos)
                                 + "\" encountered, expected one of \"" + expect + "\" (" + sFLF
                                 + ").");
    }

    return pos;
}

String::const_iterator find(const String&          input,
                            String::const_iterator pos,
                            String                 value,
                            String                 sFLF,
                            bool                   errorOnEnd)
{
    pos = std::find_if(pos,
                       input.end(),
                       [&value](const char& c)
                       { return std::find(value.begin(), value.end(), c) != value.end(); });
    if (errorOnEnd && pos == input.end())
    {
        throw std::runtime_error(
            "Invalid string: Unexpected end of content encountered, expected one character of \""
            + value + "\" (" + sFLF + ").");
    }

    return pos;
}

String::const_iterator rfind(const String&          input,
                             String::const_iterator pos,
                             String                 value,
                             String                 sFLF,
                             bool                   errorOnEnd)
{
    // Create a reverse iterator from the given forward iterator. The result rpos will not
    // correspond to pos in the sense that rpos will dereference to the character "left" of *pos.
    // To ensure that the search starts with the character at the input iterator we need to
    // decrement rpos once. However, input.rbegin() cannot be decremented, i.e. pos = input.end()
    // and pos = input.end() - 1 both start the search at the last element.
    String::const_reverse_iterator rpos(pos);
    if (rpos != input.rbegin()) rpos--;

    rpos = std::find_if(rpos,
                        input.rend(),
                        [&value](const char& c)
                        { return std::find(value.begin(), value.end(), c) != value.end(); });
    if (errorOnEnd && rpos == input.rend())
    {
        throw std::runtime_error(
            "Invalid string: Unexpected end of content encountered, expected one character of \""
            + value + "\" (" + sFLF + ").");
    }

    // If no match was found, signal this by returning input.end().
    if (rpos == input.rend()) return input.end();

    // The base() of a reverse iterators will give a forward iterator actually pointing to the next
    // item. Hence, we need to decrement the forward iterator before returning.
    return --(rpos.base());
}

} // namespace vaspml

namespace vaspml::string_tools
{

String ltrim(const String& s)
{
    Size start = s.find_first_not_of(WHITESPACE);

    return (start == String::npos) ? "" : s.substr(start);
}

String rtrim(const String& s)
{
    Size end = s.find_last_not_of(WHITESPACE);

    return (end == String::npos) ? "" : s.substr(0, end + 1);
}

String trim(const String& s)
{
    return rtrim(ltrim(s));
}

void trimVector(Vec1String& stringVector)
{
    std::transform(stringVector.begin(),
                   stringVector.end(),
                   stringVector.begin(),
                   [](String& str) { return trim(str); });
}

String makeLowerCase(const String& in_str)
{
    String out_str = in_str;
    std::transform(out_str.begin(),
                   out_str.end(),
                   out_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    return out_str;
}

String makeUpperCase(const String& in_str)
{
    String out_str = in_str;
    std::transform(out_str.begin(),
                   out_str.end(),
                   out_str.begin(),
                   [](unsigned char c) { return std::toupper(c); });

    return out_str;
}

bool checkForSubstring(const String& in_string, const String& sub_string)
{

    if (in_string.find(sub_string) != String::npos) return true;
    else return false;
}

Vec1String split(const String& inString,
                 const String& delimiter,
                 bool          removeEmpty,
                 bool          anyDelimiter)
{
    Vec1String parts;
    Size       startPos = 0;
    Size       endPos = 0;

    if (anyDelimiter)
    {
        auto findFirstDelimiter = [](const String& input, const String& delimiters, Size start)
        {
            Size firstPos = String::npos;
            for (const char c : delimiters) { firstPos = std::min(firstPos, input.find(c, start)); }
            return firstPos;
        };
        while ((endPos = findFirstDelimiter(inString, delimiter, startPos)) != String::npos)
        {
            parts.push_back(inString.substr(startPos, endPos - startPos));
            startPos = endPos + 1;
        }
    }
    else
    {
        while ((endPos = inString.find(delimiter, startPos)) != String::npos)
        {
            parts.push_back(inString.substr(startPos, endPos - startPos));
            startPos = endPos + delimiter.length();
        }
    }
    // Add the last token.
    parts.push_back(inString.substr(startPos));

    if (removeEmpty) vector_tools::removeEmptyElements(parts);

    return parts;
}

Size maxLength(const Vec1String& input)
{
    Size result = 0;
    for (const auto& i : input) result = std::max(result, i.length());

    return result;
}

String parseQuotedString(const String& input, String::const_iterator& pos)
{
    if (pos == input.end())
    {
        throw std::logic_error("Expected quoted string but iterator is already at end of input.");
    }

    // Start searching after opening quotes.
    String::const_iterator it = pos + 1;

    Int n = 0;
    while (it != input.end())
    {
        // If backslash found, increment counter, go to next character.
        if (*it == '\\')
        {
            ++n;
            ++it;
        }
        else if (*it == '"')
        {
            // No preceeding backslashes => found end of quoted string.
            if (n == 0) break;
            // Even number of preceeding backslashes => all escaped backslashes => found end.
            if (n % 2 == 0) break;
            // Odd number of preceeding backslashes => quote is escaped => continue searching.
            ++it;
            // Reset backslash counter.
            n = 0;
        }
        else
        {
            // Something other than backslash found after a sequence of backslashes.
            if (n != 0) n = 0;
            // Continue searching.
            ++it;
        }
    }

    if (it == input.end())
    {
        throw std::runtime_error("Reached end of string while parsing quoted string: unable to "
                                 "find closing double quotes.");
    }

    // Prepare result without surrounding quotes.
    String result(pos + 1, it);

    // Go to character after closing quote.
    pos = it + 1;

    return result;
}

} // namespace vaspml::string_tools

namespace vaspml::vector_tools
{

Vec1Size generateIntSequence(Size start, Size end)
{
    // assert( ( end <= start ) &&
    //         "ERROR in function VectorTools::generateIntSequence. end <= start " );
    Vec1Size indexArray(end - start);
    std::iota(indexArray.begin(), indexArray.end(), start);

    return indexArray;
}

Vec1Int invertIndex(const Vec1Int& index)
{
    Vec1Int inverted(index.size());
    for (Size i = 0; i < index.size(); ++i) { inverted[index[i]] = i; }

    return inverted;
}

} // namespace vaspml::vector_tools

namespace vaspml::file_io
{

/// second argument declares whether it is a binary file or not
std::ifstream openFileI(const String fname, const std::ifstream::openmode file_mode)
{
    std::ifstream input_stream;

    input_stream.open(fname, file_mode);
    if (!input_stream)
    {
        std::cout << "Opening infile stream failed " << std::endl;
        std::cout << "Check your file name " << fname << std::endl;
        throw;
    }

    return input_stream;
}

std::ofstream openFileO(const String fname, const std::ofstream::openmode file_mode)
{
    std::ofstream output_stream;

    output_stream.open(fname, file_mode);
    if (!output_stream)
    {
        std::cout << "Opening outfile stream failed " << std::endl;
        std::cout << "Check your file name " << fname << std::endl;
        throw;
    }

    return output_stream;
}

void writeRealVec2D(const Vec2Real& data, const String& fname)
{
    auto file = file_io::openFileO(fname);
    for (const auto& x : data)
        for (const auto& y : x) file << str("%24.16E ", y) << std::endl;
    file.close();

    return;
}

} // namespace vaspml::file_io
