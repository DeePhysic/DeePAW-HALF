#include "MemoryCounter.hpp"

#include "Item.hpp"
#include "utils.hpp"

using namespace vaspml;
using namespace vaspml::rec::detail;

MemoryCounter::MemoryCounter(Record& record) : memInfo(record)
{
    memInfo["_type"] = String{"memory_info"};
    memInfo["overhead"] = Real(0.0);
    memInfo["data"] = Real(0.0);
}

void MemoryCounter::operator()(const Record& record, bool entry)
{
    if (entry) memInfo.get<Real>("overhead") += sizeof(record);

    return;
}

// clang-format off
void MemoryCounter::operator()(const String& key, const Int&        item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Real&       item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const String&     item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const bool&       item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec1Real&   item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec1Int&    item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec1String& item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec2Real&   item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec2Int&    item) { process(key, item); }
void MemoryCounter::operator()(const String& key, const Vec2String& item) { process(key, item); }
// clang-format on
//
Size vaspml::rec::detail::memoryUsageString(const String& input)
{
    // Check if input string is smaller than SSO buffer. If yes, only the the stack memory is
    // required. If no, add stack and heap memory together.
    if (input.length() <= String{}.capacity()) return sizeof(String);
    else return sizeof(String) + input.length() * sizeof(String::value_type);
}

Size vaspml::rec::detail::memoryUsageMapEntry(const String& input)
{
    constexpr Size mem = sizeof(Item)       // Stack size of std::variant specialization.
                       + 2 * sizeof(Item*); // Pointers for double-linked list (assumption).

    return mem + memoryUsageString(input);
}

String vaspml::rec::detail::printRecordByteSizes()
{
    String result;

    result += str("sizeof(Record)     = %zu\n", sizeof(Record));
    result += str("sizeof(RecordMap)  = %zu\n\n", sizeof(RecordMap));

    result += str("sizeof(Item)       = %zu\n", sizeof(Item));
    result += str("sizeof(Real)       = %zu\n", sizeof(Real));
    result += str("sizeof(Int)        = %zu\n", sizeof(Int));
    result += str("sizeof(String)     = %zu\n", sizeof(String));
    result += str("sizeof(bool)       = %zu\n", sizeof(bool));
    result += str("sizeof(ShRec)      = %zu\n", sizeof(ShRec));
    result += str("sizeof(Vec1Real)   = %zu\n", sizeof(Vec1Real));
    result += str("sizeof(Vec1Int)    = %zu\n", sizeof(Vec1Int));
    result += str("sizeof(Vec1String) = %zu\n", sizeof(Vec1String));
    result += str("sizeof(Vec1ShRec)  = %zu\n", sizeof(Vec1ShRec));
    result += str("sizeof(Vec2Real)   = %zu\n", sizeof(Vec2Real));
    result += str("sizeof(Vec2Int)    = %zu\n", sizeof(Vec2Int));

    return result;
}
