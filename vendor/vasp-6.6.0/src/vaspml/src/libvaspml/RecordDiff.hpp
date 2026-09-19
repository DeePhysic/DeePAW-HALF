#ifndef RECORDDIFF_HPP
#define RECORDDIFF_HPP

#include "Item.hpp"
#include "Record.hpp"
#include "RecordFunctor.hpp"
#include "meta.hpp"
#include "types.hpp"

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace vaspml::rec::detail
{

struct RecordDiff : public RecordFunctor
{
    RecordDiff(const Record& lhs, const Record& rhs, bool logVectorValues, Real threshold);

    template<typename T>
    void diffEntry(String name, const T& leftValue, const T& rightValue)
    {
        currLocalDiff->get<Vec1String>("KEYS").push_back(name);
        currLocalDiff->put(name + "-1", leftValue);
        currLocalDiff->put(name + "-2", rightValue);

        return;
    }

    template<typename T>
    void diffEntryPush(String name, const T& leftValue, const T& rightValue)
    {
        currLocalDiff->get<std::vector<T>>(name + "-1").push_back(leftValue);
        currLocalDiff->get<std::vector<T>>(name + "-2").push_back(rightValue);

        return;
    }

    template<typename T>
    bool hasSameValue(const T& leftValue, const T& rightValue)
    {
        if constexpr (std::is_same_v<T, Real>)
        {
            if (std::abs(leftValue - rightValue) < threshold) return true;
        }
        else
        {
            if (leftValue == rightValue) return true;
        }

        return false;
    }

    template<typename T>
    void cmp(const String& key, const T& item)
    {
        ItemIndex type = itemIndex<T>();
        if (!isPresent(key, type)) return;

        // When traversing the right-hand-side there is no need to check everything again, diff was
        // already created. At this point, only "only-in-2" is left to fill and we can return here.
        if (traverseRun == 2) return;

        if (!hasSameType(key, type)) return;

        const auto& itemOther = currOther->cget<T>(key);
        if constexpr (isVector<T>::value)
        {
            if (item.size() != itemOther.size())
            {
                diffEntry("SIZE1:" + key,
                          static_cast<Int>(item.size()),
                          static_cast<Int>(itemOther.size()));
                // If size is different it makes little sense to continue checking all elements.
                return;
            }

            Vec1String& keysList = currLocalDiff->get<Vec1String>("KEYS");
            using T1 = typename T::value_type;
            if constexpr (isVector<T1>::value)
            {
                bool stopSearching = false;
                using T2 = typename T1::value_type;
                for (Size i = 0; i < item.size(); ++i)
                {
                    if (item[i].size() != itemOther[i].size())
                    {
                        diffEntry("SIZE2:" + key + ";" + std::to_string(i),
                                  static_cast<Int>(item[i].size()),
                                  static_cast<Int>(itemOther[i].size()));
                        // If size is different it makes little sense to continue checking all
                        // elements.
                        continue;
                    }
                    for (Size j = 0; j < item[i].size(); ++j)
                    {
                        if (!hasSameValue(item[i][j], itemOther[i][j]))
                        {
                            if (std::find(keysList.begin(), keysList.end(), key) == keysList.end())
                            {
                                keysList.push_back(key);
                                if (!logVectorValues)
                                {
                                    stopSearching = true;
                                    break;
                                }
                            }
                            String indexKey1 = "INDEX1:" + key;
                            String indexKey2 = "INDEX2:" + key;
                            if (!currLocalDiff->contains(indexKey1))
                            {
                                currLocalDiff->put<Vec1Int>(indexKey1);
                                currLocalDiff->put<Vec1Int>(indexKey2);
                                currLocalDiff->put<std::vector<T2>>(key + "-1");
                                currLocalDiff->put<std::vector<T2>>(key + "-2");
                            }
                            currLocalDiff->get<Vec1Int>(indexKey1).push_back(i);
                            currLocalDiff->get<Vec1Int>(indexKey2).push_back(j);
                            diffEntryPush(key, item[i][j], itemOther[i][j]);
                        }
                    }
                    if (stopSearching) break;
                }
            }
            else
            {
                for (Size i = 0; i < item.size(); ++i)
                {
                    if (!hasSameValue(item[i], itemOther[i]))
                    {
                        if (std::find(keysList.begin(), keysList.end(), key) == keysList.end())
                        {
                            keysList.push_back(key);
                            if (!logVectorValues) break;
                        }
                        String indexKey = "INDEX1:" + key;
                        if (!currLocalDiff->contains(indexKey))
                        {
                            currLocalDiff->put<Vec1Int>(indexKey);
                            currLocalDiff->put<std::vector<T1>>(key + "-1");
                            currLocalDiff->put<std::vector<T1>>(key + "-2");
                        }
                        currLocalDiff->get<Vec1Int>(indexKey).push_back(i);
                        diffEntryPush(key, item[i], itemOther[i]);
                    }
                }
            }
        }
        else
        {
            if (!hasSameValue(item, itemOther)) diffEntry(key, item, itemOther);
        }

        return;
    }

    bool isPresent(const String& key, ItemIndex shRecType);
    bool hasSameType(const String& key, ItemIndex shRecType);
    void updateCurrPointers();

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

    void down(const String& key, Int arrayIndex = -1, Size arraySize = 1) override;
    void up(const String& key, Int arrayIndex = -1, Size arraySize = 1) override;
    bool shouldLoop(const String& key, Size arraySize, bool skipRecordArrays) override;
    bool shouldDescend(const Record& record,
                       const String& key,
                       Int           arrayIndex = -1,
                       Size          arraySize = 1) override;

    /// Output record.
    Record result;

    const Record& lhs;
    const Record& rhs;
    /// Pointer to current result record.
    Record* currResult;
    /// Pointer to current local-diff sub-record of currResult.
    Record* currLocalDiff;
    /// Pointer to current only-in sub-record of currResult.
    Record* currOnlyIn;
    /// Pointer to current common-sub sub-record of currResult.
    Record* currCommonSub;
    /// Pointer to current "other" record.
    const Record* currOther;

    /// History of record levels of result records.
    std::vector<Record*> histResult;
    /// History of record levels of "other" record.
    std::vector<const Record*> histOther;

    /// Number of current record being traversed (1 = lhs, 2 = rhs).
    Int traverseRun;
    /// Up to this threshold floating-point numbers are considered identical.
    Real threshold;
    /// Log individual values of differing vectors (if false, show only keys).
    bool logVectorValues;
};

Record makeRecordDiffResult();
bool   emptyRecordDiffLocal(const Record& diff);
bool   emptyRecordDiffRecursive(const Record& diff);

} // namespace vaspml::rec::detail

#endif
