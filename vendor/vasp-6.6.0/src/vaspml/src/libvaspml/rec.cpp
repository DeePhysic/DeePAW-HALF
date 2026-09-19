#include "rec.hpp"

#include "Item.hpp"
#include "Jsoner.hpp"
#include "MemoryCounter.hpp"
#include "Record.hpp"
#include "RecordDiff.hpp"
#include "RecordFunctor.hpp"
#include "Serializer.hpp"
#include "Tutor.hpp"
#include "debug.hpp"
#include "utils.hpp"

#include <algorithm>   // for max_element
#include <sstream>     // for basic_stringstream, basic_iostream, bas...
#include <stdexcept>   // for runtime_error
#include <type_traits> // for is_same_v, decay_t
#include <variant>     // for get, holds_alternative, variant, visit

using namespace vaspml;
using namespace vaspml::rec::detail;

bool rec::checkType(const Record& record, String desiredType, String sFLF)
{
    if (!record.contains("_type"))
    {
        throw std::runtime_error("ERROR: Record used in " + sFLF
                                 + " is not typed (i.e. does not have a \"_type\" field).");
    }
    const String& actualType = record.cget<String>("_type");
    if (actualType != desiredType)
    {
        throw std::runtime_error("ERROR: Record used in " + sFLF + " does not have desired type \""
                                 + desiredType + "\", got \"" + actualType + "\" instead.");
    }

    return true;
}

Record rec::diff(const Record& lhs, const Record& rhs, bool logVectorValues, Real threshold)
{
    if (&lhs == &rhs) return rec::detail::makeRecordDiffResult();

    RecordDiff recordDiff(lhs, rhs, logVectorValues, threshold);
    traverse(lhs, recordDiff);
    traverse(rhs, recordDiff);

    return recordDiff.result;
}

String rec::toJson(const Record& record, Int indentSpaces, bool hideMeta)
{
    std::stringstream strm;
    toJson(record, strm, indentSpaces, hideMeta);

    return strm.str();
}

void rec::toJson(const Record& record, std::ostream& strm, Int indentSpaces, bool hideMeta)
{
    Jsoner f(strm, indentSpaces, hideMeta);
    traverse(record, f);

    return;
}

void rec::fromJson(const String& input, Record& record)
{
    if (input.empty()) return;

    String::const_iterator pos = input.cbegin();
    dejson(input, record, pos);

    return;
}

void rec::toBuffer(const Record& record, Buffer& buffer)
{
    Serializer f(buffer);
    traverse(record, f);

    return;
}

void rec::fromBuffer(const Buffer& buffer, Record& record)
{
    deserialize(buffer, record);

    return;
}

Record rec::computeMemory(const Record& record)
{
    Record        memInfo;
    MemoryCounter f(memInfo);
    traverse(record, f);

    return memInfo;
}

String rec::printMemoryUsage(const Record& memInfo)
{
    String result;

    rec::checkType(memInfo, "memory_info", flf(VASPML_FLF));

    Vec1String keys = {"Memory group", "overhead", "data", "total"};
    Real       total = 0.0;
    // Compute format specifier for key column.
    const String keyFormat = str("%%-%zus",
                                 std::max_element(keys.begin(),
                                                  keys.end(),
                                                  [](const String& a, const String& b)
                                                  { return a.length() < b.length(); })
                                     ->length());

    // Prepare table separator line.
    const String separator = str(("|-" + keyFormat + "-|-%16s-|-%12s-|-%8s-|\n").c_str(),
                                 String(keys.at(0).length(), '-').c_str(),
                                 String(16, '-').c_str(),
                                 String(12, '-').c_str(),
                                 String(8, '-').c_str());

    // Table header.
    result += "\n Memory usage estimation\n";
    result += " -----------------------\n";
    result += str(("| " + keyFormat + " | %16s | %12s | %8s |\n").c_str(),
                  keys.at(0).c_str(),
                  "kB",
                  "MB",
                  "GB");
    keys.erase(keys.begin());
    result += separator;

    // Fill table with all key-value pairs.
    for (auto k : keys)
    {
        Real mem;
        if (k == "total")
        {
            mem = total;
            result += separator;
        }
        else
        {
            mem = memInfo.cget<Real>(k);
            total += mem;
        }
        result += str(("| " + keyFormat + " | %16.3lf | %12.2lf | %8.1lf |\n").c_str(),
                      k.c_str(),
                      mem / 1.0E3,
                      mem / 1.0E6,
                      mem / 1.0E9);
    }
    result += "\n";

    return result;
}

void rec::merge(Record& first, const Record& other)
{
    RecordMap& dict = first.dict;
    for (const String& key : other.keys())
    {
        if (first.contains(key))
        {
            // If the key exists, handle merging based on the type
            const auto& type = first.itemIndexOf(key);

            if (type == ItemIndex::REAL)
            {
                // Convert to Vec1Real if not already a vector
                if (!std::holds_alternative<Vec1Real>(dict[key]))
                {
                    Real currentValue = std::get<Real>(dict[key]);
                    dict[key] = Vec1Real{currentValue};
                }
                std::get<Vec1Real>(dict[key]).push_back(other.cget<Real>(key));
            }
            else if (type == ItemIndex::INT)
            {
                // Convert to Vec1Int if not already a vector
                if (!std::holds_alternative<Vec1Int>(dict[key]))
                {
                    Int currentValue = std::get<Int>(dict[key]);
                    dict[key] = Vec1Int{currentValue};
                }
                std::get<Vec1Int>(dict[key]).push_back(other.cget<Int>(key));
            }
            else if (type == ItemIndex::STRING)
            {
                // Convert to Vec1String if not already a vector
                if (!std::holds_alternative<Vec1String>(dict[key]))
                {
                    String currentValue = std::get<String>(dict[key]);
                    dict[key] = Vec1String{currentValue};
                }
                std::get<Vec1String>(dict[key]).push_back(other.cget<String>(key));
            }
            else if (type == ItemIndex::BOOL)
            {
                global::tutor.warning(flf(VASPML_FLF)
                                      + " Cannot create a std::vector<bool> during Record merge.");
            }
            else if (type == ItemIndex::VEC1REAL)
            {
                // Append to existing Vec1Real
                auto&       currentVec = std::get<Vec1Real>(dict[key]);
                const auto& otherVec = other.cget<Vec1Real>(key);
                currentVec.insert(currentVec.end(), otherVec.begin(), otherVec.end());
            }
            else if (type == ItemIndex::VEC1INT)
            {
                // Append to existing Vec1Int
                auto&       currentVec = std::get<Vec1Int>(dict[key]);
                const auto& otherVec = other.cget<Vec1Int>(key);
                currentVec.insert(currentVec.end(), otherVec.begin(), otherVec.end());
            }
            else if (type == ItemIndex::VEC1STRING)
            {
                // Append to existing Vec1String
                auto&       currentVec = std::get<Vec1String>(dict[key]);
                const auto& otherVec = other.cget<Vec1String>(key);
                currentVec.insert(currentVec.end(), otherVec.begin(), otherVec.end());
            }
            else if (type == ItemIndex::VEC1SHREC)
            {
                // Append to existing Vec1ShRec
                auto&       currentVec = std::get<Vec1ShRec>(dict[key]);
                const auto& otherVec = other.vcget<Vec1ShRec>(key);
                //currentVec.insert(currentVec.end(), otherVec.begin(), otherVec.end());
                for (const auto& elem : otherVec)
                {
                    currentVec.push_back(
                        std::const_pointer_cast<Record>(elem)); // Explicit conversion
                }
            }
            else if (type == ItemIndex::SHREC)
            {
                // Recursively merge shared records
                auto&       currentRec = *std::get<ShRec>(dict[key]);
                const auto& otherRec = *other.cget<ShRec>(key);
                rec::merge(currentRec, otherRec);
            }
            else if (type == ItemIndex::VEC2REAL || type == ItemIndex::VEC2INT
                     || type == ItemIndex::VEC2STRING)
            {
                // Use std::visit to handle all 2D vector types dynamically
                std::visit(
                    [&](auto& currentVec)
                    {
                        using T = std::decay_t<decltype(currentVec)>;
                        if constexpr (std::is_same_v<T, Vec2Real> || std::is_same_v<T, Vec2Int>
                                      || std::is_same_v<T, Vec2String>)
                        {
                            const auto& otherVec = other.cget<T>(key);
                            currentVec.insert(currentVec.end(), otherVec.begin(), otherVec.end());
                        }
                        else
                        {
                            throw std::runtime_error("Unsupported type for merging in key: " + key);
                        }
                    },
                    dict[key]);
            }
            else { throw std::runtime_error("Unsupported type for merging in key: " + key); }
        }
        else { first.copyFrom(other, key); }
    }

    return;
}
