#ifndef BUFFER_HPP
#define BUFFER_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "meta.hpp"
#include "types.hpp"

#include <algorithm>    // for copy
#include <cstring>      // for memcpy
#include <iterator>     // for copy
#include <type_traits>  // for is_arithmetic_v, enable_if_t, is_same_v
#include <typeinfo>     // for type_info

namespace vaspml
{
enum class ItemIndex : int8_t;
}

namespace vaspml::io::detail
{

/**************************************************************************************************
 * Write given input to buffer.
 *
 * @param nchars Important for vectors and string input, determines output method (see below).
 *
 * If `nchars < 0` the sizes of vectors and strings are written before the actual contents. This
 * approach is used for the serialization of Records. If `nchars >= 0`, the sizes will not be
 * written to the output buffer. For strings, `nchars` then determines how many characters are
 * written out, too long strings will be cut, too short strings will be padded with spaces. This
 * method is used for the ML_FF reader/writer.
 **************************************************************************************************/
template<typename T>
void write(const T& input, Buffer& buffer, Int nchars = -1)
{
    // Store current position in buffer.
    Size position = buffer.size();

    // Logicals in Fortran are usually 4 bytes long.
    if constexpr (std::is_same_v<T, bool>)
    {
        // Cast bool to integer.
        int32_t logical = input;
        // Add space for single item.
        buffer.resize(buffer.size() + sizeof(int32_t));
        // Copy item to previous end of buffer.
        memcpy(&buffer[position], &logical, sizeof(int32_t));
    }
    else if constexpr (std::is_arithmetic_v<T>)
    {
        // Add space for single item.
        buffer.resize(buffer.size() + sizeof(T));
        // Copy item to previous end of buffer.
        memcpy(&buffer[position], &input, sizeof(T));
    }
    else if constexpr (std::is_same_v<T, String>)
    {
        Size newSize = buffer.size();
        // Either use input string length or get length from argument.
        if (nchars < 0)
        {
            // Compute new size from input and add one char for null terminator.
            newSize += input.size() + 1;
            buffer.resize(newSize);
            // Copy characters to previous end of buffer.
            input.copy(reinterpret_cast<char*>(&buffer[position]), input.size());
            // Add null terminator at the end of buffer.
            buffer[position + input.size()] = static_cast<Byte>('\0');
        }
        else
        {
            String tmp(input);
            // If input string is too short, pad it with spaces to the right.
            if (nchars > static_cast<Int>(tmp.size())) tmp.resize(nchars, ' ');
            // Set new size according to argument.
            newSize += nchars;
            buffer.resize(newSize);
            // Copy characters to previous end of buffer.
            tmp.copy(reinterpret_cast<char*>(&buffer[position]), nchars);
        }
    }
    else if constexpr (isVector<T>::value)
    {
        using T1 = typename T::value_type;
        Size n = input.size();
        if (nchars < 0)
        {
            // First, copy n to previous end of buffer.
            write(n, buffer);
            // Update position in buffer.
            position = buffer.size();
        }

        if constexpr (std::is_arithmetic_v<T1>)
        {
            // Add space for vector data.
            buffer.resize(buffer.size() + n * sizeof(T1));
            // Then, copy the data from the vector to buffer.
            memcpy(&buffer[position], input.data(), n * sizeof(T1));
        }
        else if constexpr (std::is_same_v<T1, String>)
        {
            // Recursively call function for each string.
            for (String const& s : input) write(s, buffer, nchars);
        }
        else if constexpr (isVector<T1>::value)
        {
            // Recursively call function for vector items.
            for (auto const& sub : input) write(sub, buffer, nchars);
        }
        else static_assert(alwaysFalse<T>, "Function is not implemented for provided vector type.");
    }
    else static_assert(alwaysFalse<T>, "Function is not implemented for provided type.");

    return;
}
/**************************************************************************************************
 * Return an error if out of buffer.
 **************************************************************************************************/
template<typename T>
void errorOutOfBuffer(const char* file, Int line, const char* func)
{
    global::tutor.bug("ERROR at " + String(file) + ":" + std::to_string(line) + " in function \""
                      + String(func)
                      + "\": Out of buffer upon attempt to read in variable of type \""
                      + typeid(T).name() + "\".");
    return;
}
/**************************************************************************************************
 * Read arithmetic value from buffer.
 **************************************************************************************************/
template<typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
void read(Buffer const& buffer, T& output, Size& position, [[maybe_unused]] Int unused = -1)
{
    if (buffer.size() - position < sizeof(T)) errorOutOfBuffer<T>(VASPML_FLF);
    memcpy(&output, &buffer[position], sizeof(T));
    position += sizeof(T);

    return;
}
/**************************************************************************************************
 * Read bool from buffer.
 **************************************************************************************************/
void read(Buffer const& buffer, bool& output, Size& position, Int unused = -1);
/**************************************************************************************************
 * Read string from buffer.
 *
 * These options are supported:
 * * nchars < 0: Determine length of string from buffer, strings are null-terminated.
 * * nchars > 0: String has nchars characters, no null-terminator required.
 * * nchars = 0: Empty string is returned.
 **************************************************************************************************/
void read(Buffer const& buffer, String& output, Size& position, Int nchars = -1);
/**************************************************************************************************
 * Read 1D-vector of arithmetic values from buffer.
 **************************************************************************************************/
template<typename T, typename std::enable_if_t<isVector1D<T>::value>* = nullptr>
void read(Buffer const& buffer, T& output, Size& position, Int nelements = -1)
{
    using T1 = typename T::value_type;
    if constexpr (!std::is_arithmetic_v<T1>)
    {
        static_assert(alwaysFalse<T>, "Function is not implemented for provided vector type.");
    }

    Size n;
    if (nelements < 0) read(buffer, n, position);
    else n = output.size();

    if (n * sizeof(T1) > buffer.size() - position) errorOutOfBuffer<T>(VASPML_FLF);
    output.resize(n);
    std::copy(reinterpret_cast<const T1*>(&buffer[position]),
              reinterpret_cast<const T1*>(&buffer[position]) + n,
              output.begin());
    position += n * sizeof(T1);

    return;
}
/**************************************************************************************************
 * Read vector of strings from buffer.
 *
 * These options are supported:
 * * nchars < 0: Determine number of elements from buffer, strings are null-terminated.
 * * nchars > 0: Determine number of elements from size of output vector, each strings has nchars
 *               characters.
 * * nchars = 0: Empty vector of strings is returned.
 **************************************************************************************************/
void read(Buffer const& buffer, Vec1String& output, Size& position, Int nchars = -1);
/**************************************************************************************************
 * Read 2D-vector of arithmetic values from buffer.
 **************************************************************************************************/
template<typename T, typename std::enable_if_t<isVector2D<T>::value>* = nullptr>
void read(Buffer const& buffer, T& output, Size& position, Int nelements = -1)
{
    using T2 = typename T::value_type::value_type;
    if constexpr (!(std::is_arithmetic_v<T2> || std::is_same_v<T2, String>))
    {
        static_assert(alwaysFalse<T>, "Function is not implemented for provided vector type.");
    }

    Size n;
    if (nelements < 0) read(buffer, n, position);
    else n = output.size();

    output.resize(n);
    for (auto& sub : output) read(buffer, sub, position, nelements);

    return;
}
/**************************************************************************************************
 * Combined read/write function for buffer operations.
 **************************************************************************************************/
template<typename T>
void process(Buffer& buffer, T& field, Size& position, bool doWrite, Int nchars = -1)
{
    if (doWrite) write(field, buffer, nchars);
    else read(buffer, field, position, nchars);

    return;
}
/**************************************************************************************************
 * Write an ItemIndex to a binary buffer.
 **************************************************************************************************/
void write(const vaspml::ItemIndex input, Buffer& buffer);
/**************************************************************************************************
 * Read an ItemIndex from a given buffer and update the position marker.
 **************************************************************************************************/
vaspml::ItemIndex readItemIndex(Buffer const& buffer, Size& position);

} // namespace vaspml::io::detail

#endif
