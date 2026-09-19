#ifndef UTILS_HPP
#define UTILS_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <algorithm>   // for sort, transform, equal, move, unique, copy
#include <fstream>     // for basic_istream, basic_ifstream, ifstream, basi...
#include <iterator>    // for distance
#include <numeric>     // for iota, accumulate
#include <sstream>     // for basic_stringstream
#include <type_traits> // for is_same, is_integral, is_arithmetic

namespace vaspml
{

/// Buffer length for vsnprintf() call in str() function.
const Size STR_MAX = 1024;
/***************************************************************************************************
 * Alternative to `printf()`, use format specifiers as usual but return string.
 *
 * Under the hood this uses `vsnprintf()` to avoid buffer overflows.
 *
 * @param format Format specifier forwarded to `vsnprintf()`.
 * @return String constructed from output of `vsnprintf()`.
 *
 * @warning As this function is calling `vsnprintf()` it inherits its problems: if the format
 * specifier does not match with actual arguments, undefined behavior may occur. Also an
 * intermediate fixed length buffer (size vaspml::STR_MAX) is used to construct the string result.
 * Hence, if the output exceeds the buffer, some characters may be lost.
 *
 * @note Originally it was intended to pass the format as a std::string. However, the NEC compiler
 * (5.0.1) does not give correct results (contains garbage strings) when `format` is passed as
 * std::string and `format.c_str()` is used as argument for `vsnprintf()`. Even when the contents of
 * `format` are copied first to a char array which is passed to `vsnprintf()` this does not help.
 *
 * @todo Test newer versions of the NEC compiler and eventually change the argument type back to
 * `std::string`.
 **************************************************************************************************/
String str(const char* format, ...);
/**************************************************************************************************
 * Stores current date and time as individual integer numbers in vector.
 *
 * @param values Vector of integers which is used as storage for individual values for year, month,
 * day, hour, minute, second and millisecond (imitate Fortran's `date_and_time()`
 * function.
 **************************************************************************************************/
void now(Vec1Int& values);
/*******************************************************************************************
 * @brief Checks if a type is a numeric type (excluding `char`).
 *
 * This function determines whether the given type `T` is a numeric type.
 * Numeric types include all arithmetic types (e.g., integers, floating-point types)
 * except for `char`, which is explicitly excluded.
 *
 * @tparam T The type to check.
 * @return true if `T` is a numeric type (excluding `char`), false otherwise.
 *
 * @note This function uses `std::is_arithmetic` to check for arithmetic types
 *       and `std::is_same` to exclude `char`.
 *******************************************************************************************/
template<typename T>
bool isNumberType()
{
    if constexpr (std::is_arithmetic<T>::value && !std::is_same<T, char>::value) return true;
    else return false;
}
/**************************************************************************************************
 * Returns current date and time as string (e.g. "2023-01-10T22:17:13.021").
 **************************************************************************************************/
String now();
/**************************************************************************************************
 * Scan for next non-whitespace character in given string.
 *
 * @param input Input string to search.
 * @param pos Starting position from where to search.
 * @param sFLF Debugging string containing file, line and function (see flf()).
 * @param expect If not empty, gives a list of characters which are expected.
 * @param errorOnEnd If true, gives an error when the end of the string is reached.
 *
 * If `expect` argument is used, and the next non-whitespace character is not part of the list an
 * error will be thrown.
 **************************************************************************************************/
String::const_iterator scan(const String&          input,
                            String::const_iterator pos,
                            String                 sFLF,
                            String                 expect = "",
                            bool                   errorOnEnd = true);
/**************************************************************************************************
 * Search for any of the given characters in given string.
 *
 * @param input Input string to search.
 * @param pos Starting position from where to search.
 * @param value String containing all characters which should be searched for.
 * @param sFLF Debugging string containing file, line and function (see flf()).
 * @param errorOnEnd If true, gives an error when the end of the string is reached.
 * @return String iterator corresponding to position where one of the desired characters was found.
 **************************************************************************************************/
String::const_iterator find(const String&          input,
                            String::const_iterator pos,
                            String                 value,
                            String                 sFLF,
                            bool                   errorOnEnd = true);
/**************************************************************************************************
 * Search for any of the given characters in given string in reverse direction.
 *
 * @param input Input string to search.
 * @param pos Starting position from where to search in reverse direction.
 * @param value String containing all characters which should be searched for.
 * @param sFLF Debugging string containing file, line and function (see flf()).
 * @param errorOnEnd If true, gives an error when the end of the string is reached.
 * @return String iterator corresponding to position where one of the desired characters was found.
 *
 * @important If `errorOnEnd` is `false` and the search was unsuccessful `input.end()` is returned.
 * Although this is in the "wrong" direction it is the only way to signal that no match was found
 * (input.begin() would be ambiguous).
 **************************************************************************************************/
String::const_iterator rfind(const String&          input,
                             String::const_iterator pos,
                             String                 value,
                             String                 sFLF,
                             bool                   errorOnEnd = true);
} // namespace vaspml

namespace vaspml::string_tools
{

/// define a symbol for white-space
inline const String WHITESPACE = " \n\r\t\f\v";

/** left trim an input string
 * @param s -> input string to trim will not change
 */
String ltrim(const String& s);

/** right trim an input string
 * @param s -> input string to trim will not change
 */
String rtrim(const String& s);

/** apply left and right trim to remove trailing
 *  and ending white spaces
 * @param s -> input string to trim will not change
 */
String trim(const String& s);

/*******************************************************************************************
 * trimming the elements of string vector element wise.
 * @param stringVector the leading and trailing spaces of this vector will be removed
 *******************************************************************************************/
void trimVector(Vec1String& stringVector);

/** extract single value from string
 * @param data -> string containing single value to extract
 */
template<typename T>
T extractValue(const String& data)
{
    T value;
    if constexpr (std::is_same<T, int>::value) { value = std::stoi(trim(data).c_str()); }
    else if constexpr (std::is_same<T, float>::value) { value = std::stof(trim(data).c_str()); }
    else if constexpr (std::is_same<T, double>::value) { value = std::stod(trim(data).c_str()); }
    else if constexpr (std::is_same<T, String>::value) { value = trim(data).c_str(); }

    return value;
}

/** extract Vector from a string separated by spaces
 * @param data -> input string containing numbers separated by space
 */
template<typename T>
std::vector<T> extractData(const String& input)
{
    String            temp;
    std::vector<T>    data;
    std::stringstream strStream(input.c_str());
    while (strStream >> temp)
    {
        if constexpr (std::is_same<T, int>::value) { data.push_back(std::stoi(temp)); }
        else if constexpr (std::is_same<T, float>::value) { data.push_back(std::stof(temp)); }
        else if constexpr (std::is_same<T, double>::value) { data.push_back(std::stod(temp)); }
        else if constexpr (std::is_same<T, String>::value) { data.push_back(temp); }
    }

    return std::vector<T>(data.begin(), data.end());
}
/*
 * convert a string to lower case letters ( copy of string will be made )
 * @param in_str input string to be converted to lower case
 *
 */
String makeLowerCase(const String& in_str);
/*
 * convert a string to upper case letters ( copy of string will be made )
 * @param in_str input string to be converted to upper case case
 *
 */
String makeUpperCase(const String& in_str);
/*
 *
 * check if a string contains a substring ( function if case sensitive )
 * @param in_string string which will be searched to contain
 * @param sub_string the subsrting for which will be searched to be conatined
 * in in_string
 *
 */
bool checkForSubstring(const String& in_string, const String& sub_string);

/***************************************************************************************************
 * Split string at given delimeter string and return vector of string parts.
 *
 * @param inString Input string which will be split. Will not be destroyed.
 * @param delimiter String at which the input string will be split (can be also multiple
 *                  characters).
 * @param removeEmpty Removes empty string parts from the return vector at the end of the function.
 * @param anyDelimiter If true, any of the characters in the delimiter string may be used to split.
 *                     If false, only the delimiter as a whole splits the string.
 *
 * @note If `removeEmpty` is set to `false`, multiple consecutive delimiters may result in unwanted
 *       empty return vector entries. For example, splitting `a__b` (two underscores between
 *       `a` and `b`) with underscore as delimiter, results in a vector containing `"a"`, `""`
 *       and `"b"`!
 ***************************************************************************************************/
Vec1String split(const String& inString,
                 const String& delimiter,
                 bool          removeEmpty = true,
                 bool          anyDelimiter = false);
/**************************************************************************************************
 * Get maximum length of vector of strings.
 *
 * @param input Vector of strings.
 * @return Length of the longest string in input vector.
 **************************************************************************************************/
Size maxLength(const Vec1String& input);
/**************************************************************************************************
 * Search for end of quoted string, return it (without quotes) and advance to character after
 * closing quote.
 *
 * This function is useful if the string itself is allowed to contain an escaped double quote (\")
 * and also backslashes can appear if escaped (\\). In this case it is not trivial to find the
 * actual end of the string. Used for example when deserializing a JSON file or parsing the extended
 * XYZ comment line.
 *
 * @param input Input string to search.
 * @param pos Position in input string with start of string (first quote). On exit, points to
 *            character right of closing quote. Must not be same as input.end().
 * @return String between quotes (without the surrounding quotes).
 **************************************************************************************************/
String parseQuotedString(const String& input, String::const_iterator& pos);
} // namespace vaspml::string_tools

namespace vaspml::vector_tools
{
/** sort vector input_vector and return the indexes in sorted order
 *  input_vector will not be destroyed
 * @param input_vector  -->  vector containing numerical T to sort
 * @param sorted_indx -> contains the indices of the vector in sorted order
 */
template<typename T>
Vec1Size argSort(const std::vector<T>& input_vector)
{
    Vec1Size sorted_indx(input_vector.size());
    std::iota(sorted_indx.begin(), sorted_indx.end(), 0);

    // sort the sorted_indx array based on lambda funtcion
    // based on input_vector
    std::sort(sorted_indx.begin(),
              sorted_indx.end(),
              // def lambda function that defines the order
              [&input_vector](Size indx1, Size indx2) -> bool
              { return input_vector[indx1] < input_vector[indx2]; });

    return sorted_indx;
}
/**
 *check if two vectors are equal
 *
 *if the size of the vectors does not match they are said to be not equal
 *
 *@param vec1 first vector for comparison
 *@param vec2 second vector for comprison
 *
 *@note the tmeplate type T of the functions has to implement elementwise
 *comparison operator
 */
template<typename T>
bool checkEqual(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    if (vec1.size() != vec2.size()) { return false; }

    return std::equal(vec1.begin(), vec1.end(), vec2.begin(), vec2.end());
}

template<typename T>
Vec1Size argSortSlice(const std::vector<T>& input_vector, const Size start, const Size end)
{
    Vec1Size sorted_indx(end - start);
    std::iota(sorted_indx.begin(), sorted_indx.end(), start);
    // sort the sorted_indx array based on lambda funtcion
    // based on input_vector
    std::sort(sorted_indx.begin(),
              sorted_indx.end(),
              // def lambda function that defines the order
              [&input_vector](Size indx1, Size indx2) -> bool
              { return input_vector[indx1] < input_vector[indx2]; });

    return sorted_indx;
}

/** return slice of a vector, note that the copy constructor will be called, can be expensive
 * @param array -> array from which a slice will be extracted
 * @param indx1 -> start indx from which slice starts
 * @param indx2 -> exclusive indx at which slice ends
 */
template<typename T>
std::vector<T> sliceVector(const std::vector<T>& array, Size indx1, Size indx2)
{
    return std::vector<T>(array.begin() + indx1, array.begin() + indx2);
}
/**
 * allocate a 1d array, all elements are set to zero
 * @param array -> input std::vector<T> to resize
 * @param Int n -> allocation size
 */
template<typename T>
void allocate_vector(std::vector<T>& array, const Size n)
{
    array.resize(n);

    return;
}
/**
 * allocate a 2d array, all elements are set to zero
 * @param array  ->  input std::vector<std::vector<T>> to resize
 * @param Int n0 ->  allocation size for first dimension
 * @param Int n1 ->  allocation size for second dimension
 */
template<typename T>
void allocate_vector(std::vector<std::vector<T>>& array, const Size n0, const Size n1)
{
    array.resize(n0);
    for (Size i = 0; i < n0; i++) { array[i].resize(n1); }

    return;
}
/**
 * allocate a 3d array, all elements are set to zero
 * @param array  ->  input std::vector<std::vector<T>> to resize
 * @param Int n0 ->  allocation size for first dimension
 * @param Int n1 ->  allocation size for second dimension
 * @param Int n2 ->  allocation size for second dimension
 */
template<typename T>
void allocate_vector(std::vector<std::vector<std::vector<T>>>& array,
                     const Size                                n0,
                     const Size                                n1,
                     const Size                                n2)
{
    array.resize(n0);
    for (Size i = 0; i < n0; i++)
    {
        array[i].resize(n1);
        for (Size j = 0; j < n1; j++) { array[i][j].resize(n2); }
    }

    return;
}
/*******************************************************************************************
 * allocate a std::vector<std::vector<T>> with varaible lengths for every entry
 *
 * @param array will be resized to the desired lengths on output
 * @param sizes contains the array lengths for every entry in the first dimension of the array
 *
 * array[ 0 ].resize( sizes[0] ), array[1].resize( sizes[1] ) ....
 *******************************************************************************************/
template<typename T>
void allocate_vector(std::vector<std::vector<T>>& array, const Vec1Int& sizes)
{
    array.resize(sizes.size());
    for (Size i = 0; i < sizes.size(); i++) { array[i].resize(sizes[i]); }

    return;
}
/**
 * generate a list of integer starting at start and ending at end
 * @param start start integer
 * @param end integer
 */
Vec1Size generateIntSequence(Size start, Size end);

/**
 * generate a list of integer starting at start and ending at end
 * @param start start integer
 * @param end integer
 */
Vec1Size generateIntSequence(Size start, Size end);

/*
 * get unique elements of vector without destruction
 *@param input_vector -> vector from which unique elements are taken
 *@return unique elements of vector
 */
template<typename T>
std::vector<T> get_unique(const std::vector<T>& input_vector)
{
    std::vector<T>                    output_vector = input_vector;
    typename std::vector<T>::iterator ip;
    std::sort(output_vector.begin(), output_vector.end());
    ip = std::unique(output_vector.begin(), output_vector.end());
    output_vector.resize(std::distance(output_vector.begin(), ip));

    return output_vector;
}

/*
 * get unique elements of vector with destruction of input vector
 *@param input_vector -> vector from which unique elements are taken
 *                       on output vector contains unique elements
 */
template<typename T>
void get_unique_in_place(std::vector<T>& vector)
{
    typename std::vector<T>::iterator ip;
    std::sort(vector.begin(), vector.end());
    ip = std::unique(vector.begin(), vector.end());
    vector.resize(std::distance(vector.begin(), ip));

    return;
}

/**
 *
 * count the number of occurences of
 *@param element in the std::vector
 *@param vec_in
 */
template<typename T>
T count_element(const std::vector<T>& vec_in, const T element)
{
    return count(vec_in.begin(), vec_in.end(), element);
}

/**
 * function return the maximum value in your std::vector
 *@param vec input vector from which to extract max value
 */
template<typename T>
T get_max(const std::vector<T>& vec)
{
    return *max_element(vec.begin(), vec.end());
}

/**
 * function return the minimum value in your std::vector
 *@param vec input vector from which to extract min value
 */
template<typename T>
T get_min(const std::vector<T>& vec)
{
    return *min_element(vec.begin(), vec.end());
}

/*******************************************************************************************
 * summing all elements of e vector
 * @param input vector which has to  be summed over
 *******************************************************************************************/
template<typename T>
T sum(const std::vector<T>& input)
{
    return std::accumulate(input.begin(), input.end(), (T)0);
}

/*******************************************************************************************
 * making a sum over a vector slice
 *
 * @param input input vector over which the sum is taken
 * @param start index of vector can at least be zero and maximal end -1
 * @param end can at least be start + 1 and maximimal input.size()-1
 *******************************************************************************************/
template<typename T>
T sum(const std::vector<T>& input, const Size start, const Size end)
{
    VASPML_DEBUG_L1(
        if (end < start)
        {
            global::tutor.bug("ERROR: sum( const std::vector<T>& input, const Size start, const"
                              " Size end )\n"
                              "       end index has to be larger than start index");
        }
    );

    return std::accumulate(input.begin() + start, input.end() + end, T::value_type(0));
}

/*******************************************************************************************
 * Makes the Hadamard product of two input vectors
 * @param vec1 input vector one which will be multiplied
 * @param vec2 input vector two which will be multiplied with vec1
 *
 * @return result \f[=vec1[i]*vec2[i]\f]
 *******************************************************************************************/
template<typename T>
std::vector<T> elementWiseProduct(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    VASPML_DEBUG_L1(
        if (vec1.size() != vec2.size())
        {
            global::tutor.bug(
                "std::vector<T> elementWiseProduct( const std::vector<T>& vec1, const"
                " std::vector<T>& vec2 )\n"
                "       size of vec1 and vec2 do not agree. vec1.size() != vec2.size()");
        }
    );

    std::vector<T> result(vec1.size());
    std::transform(vec1.begin(),
                   vec1.end(),
                   vec2.begin(),
                   result.begin(),
                   [](const T& a, const T& b) { return a * b; });

    return result;
}

/**
 * @brief Converts a 1D vector of real values into a 2D vector based on specified sizes.
 *
 * This function takes a 1D vector (`dataIn`) and splits it into a 2D vector (`dataOut`)
 * where each inner vector's size is determined by the corresponding value in the `sizes` vector.
 *
 * @tparam T The type of the elements in the input and output vectors.
 * @param dataIn The input 1D vector containing real values to be split.
 * @param sizes A vector of integers specifying the sizes of the inner vectors in the output.
 * @return A 2D vector where each inner vector has the size specified by the corresponding value in
 * `sizes`.
 *
 * @throws std::runtime_error If the sum of `sizes` does not match the size of `dataIn`.
 *
 * @note The function performs a debug-level check to ensure the sum of `sizes` matches the size of
 * `dataIn`. If the check fails, a debug message is logged, and an error is raised.
 */
template<typename T>
std::vector<std::vector<T>> Vec1ToVec2(const std::vector<T>& dataIn, const Vec1Int& sizes)
{
    VASPML_DEBUG_L1(
        if (std::accumulate(sizes.begin(), sizes.end(), 0) != (Int)dataIn.size())
        {
            global::tutor.bug("vaspml::vector_tools::Vec1RealToVec2Real: "
                              "The sum of sizes does not match the size of dataIn.");
        }
    );
    std::vector<std::vector<T>> dataOut(sizes.size());
    Size                        start = 0;
    for (Size i = 0; i < sizes.size(); i++)
    {
        dataOut[i].resize(sizes[i]);
        std::copy(dataIn.begin() + start, dataIn.begin() + start + sizes[i], dataOut[i].begin());
        start += sizes[i];
    }

    return dataOut;
}

template<typename T>
std::vector<T> Vec2ToVec1(const std::vector<std::vector<T>>& dataIn)
{
    Int size = 0;
    for (const auto& sub : dataIn) size += sub.size();
    std::vector<T> dataOut(size);
    Size counter = 0;
    std::for_each(dataIn.cbegin(), dataIn.cend(), [&](const std::vector<T>& row) 
            {
               std::copy(row.begin(), row.end(), dataOut.begin() + counter);
               counter += row.size();
            });
    return dataOut;
}

/*******************************************************************************************
 * @brief Computes the sizes of inner vectors in a 2D vector.
 *
 * @tparam T The type of elements stored in the inner vectors.
 * @param input A 2D vector whose inner vector sizes are to be computed.
 * @return A vector containing the sizes of each inner vector in the input.
 *
 * This function takes a 2D vector as input and returns a 1D vector where each
 * element represents the size of the corresponding inner vector in the input.
 *******************************************************************************************/
template<typename T>
std::vector<Int> getSizes(const std::vector<std::vector<T>>& input)
{
    std::vector<Int> sizes(input.size());
    std::transform(input.cbegin(),
                   input.cend(),
                   sizes.begin(),
                   [](const std::vector<T>& vec) { return vec.size(); });

    return sizes;
}
/*******************************************************************************************
 * @brief Filters elements from an input vector based on specified indices.
 *
 * This function creates a new vector containing elements from the input vector `vecIn`
 * at the positions specified by the indices in `indx`. If any index in `indx` is out
 * of bounds for `vecIn`, an error is logged.
 *
 * @tparam T The type of elements in the input vector.
 * @param vecIn The input vector from which elements are to be filtered.
 * @param indx A vector of indices specifying the positions of elements to extract from `vecIn`.
 * @return A vector containing the elements from `vecIn` at the specified indices.
 *
 * @note If `indx` contains indices that are out of bounds for `vecIn`, an error is logged
 *       using the `VASPML_DEBUG_L1` macro.
 *
 * @warning Ensure that all indices in `indx` are within the valid range of `vecIn` to avoid
 *          undefined behavior.
 *******************************************************************************************/
template<typename T, typename IndexType>
std::vector<T> filterVector(const std::vector<T>& vecIn, const std::vector<IndexType>& indx)
{
    VASPML_DEBUG_L1(
        static_assert(std::is_integral<IndexType>::value, "IndexType must be an integral type.");
        if (!indx.empty() && (Size)*std::max_element(indx.begin(), indx.end()) >= vecIn.size())
        {
            global::tutor.bug(
                "Error: " + flf(VASPML_FLF)
                + ": One or more indices in indx are out of bounds for vector vecIn.");
        }
    );
    std::vector<T> filtered(indx.size());
    std::transform(indx.cbegin(),
                   indx.cend(),
                   filtered.begin(),
                   [&](const IndexType index) { return vecIn[index]; });

    return filtered;
}

// TODO DELETE
/*******************************************************************************************
 * @brief Extracts elements from a vector based on specified indexes.
 *
 * @tparam T The type of elements in the input vector.
 * @param data The input vector containing elements to extract.
 * @param indexes A vector of integer indexes specifying which elements to extract.
 * @return A vector containing the elements from the input vector at the specified indexes.
 *
 * @note Assumes that all indexes in the `indexes` vector are valid and within bounds of the `data`
 *vector.
 *******************************************************************************************/
template<typename T>
std::vector<T> getElements(const std::vector<T>& data, const Vec1Int& indexes)
{
    std::vector<T> result(indexes.size());

    std::transform(indexes.cbegin(),
                   indexes.cend(),
                   result.begin(),
                   [&data](const Int index) { return data[index]; });

    return result;
}

/**************************************************************************************************
 * Remove empty elements from vector.
 **************************************************************************************************/
template<typename T>
void removeEmptyElements(std::vector<T>& vec)
{
    auto it = vec.begin();
    while (it != vec.end())
    {
        if (it->empty()) it = vec.erase(it);
        else ++it;
    }

    return;
}
/*******************************************************************************************
 * @brief Reorders the elements of a vector based on a given index mapping.
 *
 * This function rearranges the elements of the input vector `data` according to
 * the order specified by the `index` vector. Each element in `index` specifies
 * the new position of the corresponding element in `data`.
 *
 * @tparam T The type of elements in the vector.
 * @param data Reference to the vector to be reordered.
 * @param index A vector of integers specifying the new order of elements.
 *
 * @note The size of `index` must match the size of `data`, and `index` must
 *       contain a valid permutation of indices (i.e., no duplicates and all
 *       indices in the range [0, data.size()-1]).
 *******************************************************************************************/
template<typename T>
void reorderVector(std::vector<T>& data, const Vec1Int& index)
{
    std::vector<T> temp(data.size());
    for (Size i = 0; i < data.size(); i++) { temp[i] = data[index[i]]; }
    data = std::move(temp);

    return;
}

/*******************************************************************************************
 * @brief Reorders the elements of a vector in-place based on the given index mapping.
 *
 * This function rearranges the elements of the input vector `data` according to the
 * mapping specified in the `index` vector. The reordering is performed in-place
 * without allocating additional memory for the reordered vector, except for a
 * temporary buffer to track visited indices.
 *
 * @tparam T The type of elements in the vector.
 * @param[in,out] data The vector to be reordered. The elements are rearranged in-place.
 * @param[in] index A vector of integers specifying the new order of elements. Each
 *                  element at position `i` in `index` indicates the target position
 *                  of the element at position `i` in `data`. The size of `index`
 *                  must match the size of `data`.
 *
 * @note The function assumes that the `index` vector represents a valid permutation
 *       of indices (i.e., it contains all integers from `0` to `data.size() - 1`).
 *       Undefined behavior may occur if this condition is not met.
 *
 * @warning The function modifies the input vector `data` directly. Ensure that
 *          the input vector is not const.
 *
 *******************************************************************************************/
template<typename T>
void inplaceReorderVector(std::vector<T>& data, const Vec1Int& index)
{
    std::vector<bool> done(data.size(), false);
    for (Size i = 0; i < data.size(); ++i)
    {
        if (done[i] || (Size)index[i] == i) continue;
        Size j = i;
        T    temp = std::move(data[i]);
        while (!done[j])
        {
            done[j] = true;
            Size next = index[j];
            if (next == i) data[j] = std::move(temp);
            else data[j] = std::move(data[next]);
            j = next;
        }
    }

    return;
}
/*******************************************************************************************
 * @brief Inverts the given index vector.
 *
 * This function takes a vector of integers, where each element represents an index mapping,
 * and returns a new vector where the positions and values are inverted. For example, if the
 * input vector represents a permutation, the output vector will represent its inverse.
 *
 * @param index The input vector of integers representing the index mapping.
 * @return A vector of integers where the positions and values are inverted.
 *
 * @note The input vector must contain a valid permutation of indices (0 to n-1) for the function
 * to work correctly.
 *******************************************************************************************/
Vec1Int invertIndex(const Vec1Int& index);
} // namespace vaspml::vector_tools

namespace vaspml::file_io
{

/** open file for input reading
 * @param fname -> string containing filename
 */
std::ifstream openFileI(const String                  fname,
                        const std::ifstream::openmode file_mode = std::ifstream::in);
/** open file for input writing
 * @param fname -> string containing filename
 */
std::ofstream openFileO(const String                  fname,
                        const std::ofstream::openmode file_mode = std::ifstream::out);

/*******************************************************************************************
 * dumb a 2d Real vector to a file in a single column format. File format second index fast
 *
 * @param[in] data array which will be written to the file.
 * @param[in] fname file name to which data will be written
 *******************************************************************************************/
void writeRealVec2D(const Vec2Real& data, const String& fname);

} // namespace vaspml::file_io

namespace vaspml::read_values
{

/// scalar reader
template<typename T, typename V>
void io_scalar(V& io_stream, T& scalar)
{
    // Do reading if input stream
    if constexpr (std::is_same<V, std::ifstream>::value)
    {

        if constexpr (std::is_integral<T>::value && std::is_same<T, bool>::value)
        {
            Int int_value;
            io_stream.read((char*)&int_value, sizeof(int_value));
            scalar = T(int_value);
        }
        else { io_stream.read((char*)&scalar, sizeof(scalar)); }
    }
    // Do writing
    else
    {
        if constexpr (std::is_integral<T>::value && std::is_same<T, bool>::value)
        {
            Int int_value = (Int)scalar;
            io_stream.write((char*)&int_value, sizeof(int_value));
        }
        else { io_stream.write((char*)&scalar, sizeof(scalar)); }
    }

    return;
}

/// array reader
template<typename T, typename V>
void io_array_binary(V& io_stream, std::vector<T>& array, const Int& number_entries)
{
    // Do reading if input stream
    if constexpr (std::is_same<V, std::ifstream>::value)
    {
        if constexpr (std::is_same<T, String>::value) // value here is member variable of
                                                      // is_same!!!
        {
            // Character reading
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                char value[3];

                io_stream.read((char*)&value, sizeof(value) - 1);
                value[2] = '\0';

                array.push_back(value);
            }
        }
        else
        {
            // Non-character reading
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                T value;

                io_stream.read((char*)&value, sizeof(value));
                array.push_back(value);
            }
        }
    }
    // Otherwise do writing
    else
    {
        if constexpr (std::is_same<T, String>::value) // value here is member variable of
                                                      // is_same!!!
        {
            // Character writing
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                io_stream.write(array[i_entry].c_str(), 2);
            }
        }
        else
        {
            for (Int i_entry = 0; i_entry < number_entries; i_entry++)
            {
                io_stream.write((char*)&array[i_entry], sizeof(array[i_entry]));
            }
        }
    }

    return;
}

} //namespace vaspml::read_values

#endif
