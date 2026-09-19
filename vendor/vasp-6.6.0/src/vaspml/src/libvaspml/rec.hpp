#ifndef REC_HPP
#define REC_HPP

#include "types.hpp"

#include <limits>
#include <ostream>

namespace vaspml
{
struct Record;
}

/// Functions dealing with general operations on Records.
namespace vaspml::rec
{

/***************************************************************************************************
 * Check if a given Record is of the correct "type" (contains "_type" metadata entry).
 *
 * @param record Record to check (const, will not be modified).
 * @param desiredType The desired "type" of the Record.
 * @param sFLF Debugging string containing file, line and function (see flf()).
 * @return True if correct type specifier is found.
 **************************************************************************************************/
bool checkType(const Record& record, String desiredType, String sFLF);
/**************************************************************************************************
 * Compare two Records and return a "diff" Record.
 *
 * @param lhs First input Record.
 * @param rhs Second input Record.
 * @param threshold Up to this threshold floating-point numbers are considered identical (compared
 * to absolute value of difference).
 * @return Record containing entries for "only-in-left", "only-in-right", "diff".
 **************************************************************************************************/
Record diff(const Record& lhs,
            const Record& rhs,
            bool          logVectorValues = true,
            Real          threshold = 10 * std::numeric_limits<Real>::epsilon());
/**************************************************************************************************
 * Convert given Record to JSON (String variant), see ostream variant for argument details.
 *
 * Generates a string output, internally calls `std::ostream` version of toJson().
 **************************************************************************************************/
String toJson(const Record& record, Int indentSpaces = 2, bool hideMeta = true);
/**************************************************************************************************
 * Convert given Record to JSON (std::ostream variant).
 *
 * @param indentSpaces Specifies indentation level given in number of spaces. Special meaning of
 * negative numbers (see below).
 * @param hideMeta If true, metadata, i.e. all entries with keys starting with "_" will not be
 *
 * If `indentSpaces` is negative the formatting will be special:
 *
 * * `indentSpaces == -1`: Removes all newlines but keeps some spaces for readability.
 * * `indentSpaces == -2`: Removes all newlines and also most spaces (real numbers with
 * scientific notation may have some leading spaces remaining).
 *
 * Writes JSON representation to output stream. Internally uses ::Jsoner functor and traverse()
 * function to walk through Record.
 **************************************************************************************************/
void toJson(const Record& record, std::ostream& strm, Int indentSpaces = 2, bool hideMeta = true);
/**************************************************************************************************
 * Convert input string to Record.
 **************************************************************************************************/
void fromJson(const String& input, Record& record);
/**************************************************************************************************
 * Convert given Record to binary representation.
 **************************************************************************************************/
void toBuffer(const Record& record, Buffer& buffer);
/**************************************************************************************************
 * Convert buffer contents to Record.
 **************************************************************************************************/
void fromBuffer(const Buffer& buffer, Record& record);
/**************************************************************************************************
 * Compute memory footprint of given Record and return "mem_info" Record.
 **************************************************************************************************/
Record computeMemory(const Record& record);
/**************************************************************************************************
 * Create string report from Record containing memory information of other Record.
 **************************************************************************************************/
String printMemoryUsage(const Record& memInfo);
/*******************************************************************************************
 * @brief Merges the contents of one Record object into another.
 *
 * This function merges the data from the `other` Record into the `first` Record.
 * If a key exists in both records, the merging behavior depends on the type of the value:
 * - **REAL, INT, STRING**: Converts the value to a vector (if not already) and appends the new
 *value.
 * - **VEC1REAL, VEC1INT, VEC1STRING**: Appends the elements of the vector from `other` to the
 *vector in `first`.
 * - **BOOL**: Logs a warning about potential creation of `std::vector<bool>`.
 * - **SHREC**: Recursively merges shared Record objects.
 * - **VEC1SHREC**: Appends shared Record pointers, converting them to non-const pointers.
 * - **VEC2REAL, VEC2INT, VEC2STRING**: Dynamically appends 2D vector elements using `std::visit`.
 * - Unsupported types throw a `std::runtime_error`.
 *
 * If a key exists in `other` but not in `first`, the key-value pair is directly added to `first`.
 *
 * @param first The Record object to merge into. This object is modified.
 * @param other The Record object to merge from. This object remains unchanged.
 * @throws std::runtime_error If an unsupported type is encountered during merging.
 *******************************************************************************************/
void merge(Record& first, const Record& other);

} // namespace vaspml::rec

#endif
