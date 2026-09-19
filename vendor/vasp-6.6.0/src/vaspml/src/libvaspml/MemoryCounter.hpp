#ifndef MEMORYCOUNTER_HPP
#define MEMORYCOUNTER_HPP

#include "Record.hpp"
#include "RecordFunctor.hpp"
#include "meta.hpp"
#include "types.hpp"

namespace vaspml::rec::detail
{

/**************************************************************************************************
 * Compute memory usage of map entry.
 **************************************************************************************************/
Size memoryUsageMapEntry(const String& input);

struct MemoryCounter : public RecordFunctor
{
    MemoryCounter(Record& record);

    template<typename T>
    void process(const String& key, const T& item)
    {
        const Real overhead = memoryUsageMapEntry(key);
        Real       data = sizeof(T);

        if constexpr (isVector<T>::value)
        {
            using T1 = typename T::value_type;
            data += item.size() * sizeof(T1);
            if constexpr (isVector<T1>::value)
            {
                using T2 = typename T1::value_type;
                for (const auto& i : item) data += i.size() * sizeof(T2);
            }
        }

        memInfo.get<Real>("overhead") += overhead;
        memInfo.get<Real>("data") += data;

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

    Record& memInfo;
};

/**************************************************************************************************
 * Compute memory usage of a given string.
 **************************************************************************************************/
Size memoryUsageString(const String& input);
/**************************************************************************************************
 * Create string with overview of byte sizes for value types in a Record.
 **************************************************************************************************/
String printRecordByteSizes();

} // namespace vaspml::rec::detail

#endif
