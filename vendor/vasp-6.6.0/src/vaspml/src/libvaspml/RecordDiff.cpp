#include "RecordDiff.hpp"

#include "SmartEnum.hpp"
#include "debug.hpp"
#include "rec.hpp"

#include <stdexcept>

using namespace vaspml::rec::detail;

RecordDiff::RecordDiff(const Record& lhs, const Record& rhs, bool logVectorValues, Real threshold) :
    result(makeRecordDiffResult()),
    lhs(lhs),
    rhs(rhs),
    currResult(nullptr),
    currLocalDiff(nullptr),
    currOnlyIn(nullptr),
    currCommonSub(nullptr),
    currOther(nullptr),
    traverseRun(0),
    threshold(threshold),
    logVectorValues(logVectorValues)
{}

bool RecordDiff::isPresent(const String& key, ItemIndex type)
{
    if (!currOther->contains(key))
    {
        currOnlyIn->put(key, SmartEnum<ItemIndex>::toString(type));
        return false;
    }
    else return true;
}

bool RecordDiff::hasSameType(const String& key, ItemIndex type)
{
    if (currOther->itemIndexOf(key) != type)
    {
        // First pass has already created diff entry.
        if (traverseRun < 2)
        {
            diffEntry("TYPE:" + key,
                      SmartEnum<ItemIndex>::toString(type),
                      SmartEnum<ItemIndex>::toString(currOther->itemIndexOf(key)));
        }
        return false;
    }
    else return true;
}

void RecordDiff::updateCurrPointers()
{
    currLocalDiff = &currResult->dget<ShRec>("local-diff");
    currCommonSub = &currResult->dget<ShRec>("common-sub");
    if (traverseRun == 1) currOnlyIn = &currResult->dget<ShRec>("only-in-1");
    else currOnlyIn = &currResult->dget<ShRec>("only-in-2");

    return;
}

void RecordDiff::operator()(const Record& /*record*/, bool entry)
{
    if (entry && level == 0)
    {
        traverseRun++;
        if (traverseRun > 2)
        {
            throw std::runtime_error("RecordDiff instance can only be used once per side of "
                                     "comparison");
        }
        else
        {
            currResult = &result;
            updateCurrPointers();
            if (traverseRun == 1) currOther = &rhs;
            else currOther = &lhs;
        }
    }

    return;
}

// clang-format off
void RecordDiff::operator()(const String& key, const Int&        item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Real&       item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const String&     item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const bool&       item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Vec1Real&   item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Vec1Int&    item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Vec1String& item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Vec2Real&   item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Vec2Int&    item) { cmp(key, item); }
void RecordDiff::operator()(const String& key, const Vec2String& item) { cmp(key, item); }
// clang-format on

void RecordDiff::down(const String& key, Int arrayIndex, Size /* arraySize */)
{
    level++;

    // Save currResult pointer to history.
    histResult.push_back(currResult);
    // Point currResult to corresponding sub-record in common-sub.
    if (arrayIndex < 0) currResult = &currCommonSub->dget<ShRec>(key);
    else currResult = &*currCommonSub->get<Vec1ShRec>(key).back();

    // Update all pointers below currResult.
    updateCurrPointers();

    // Save currOther pointer to history.
    histOther.push_back(currOther);
    // Point currOther to corresponding sub-record for given key.
    if (arrayIndex < 0) currOther = &currOther->dcget<ShRec>(key);
    else currOther = &*currOther->vcget<Vec1ShRec>(key)[arrayIndex];

    return;
}

void RecordDiff::up(const String& /*key*/, Int /*arrayIndex*/, Size /*arraySize*/)
{
    level--;

    // Retrieve last currResult pointer from history.
    currResult = histResult.back();
    histResult.pop_back();

    // Update all pointers below currResult.
    updateCurrPointers();

    // Retrieve last currOther pointer from history.
    currOther = histOther.back();
    histOther.pop_back();

    return;
}

bool RecordDiff::shouldDescend(const Record& /* record */,
                               const String& key,
                               Int           arrayIndex,
                               Size          arraySize)
{
    // Only single ShRec entries are handled here, Vec1ShRec in shouldLoop().
    if (arrayIndex < 0)
    {
        // Only Vec1ShRec is handled here, single ShRec entries in shouldDescend().
        if (!isPresent(key, ItemIndex::SHREC)) return false;
        if (!hasSameType(key, ItemIndex::SHREC)) return false;
    }

    if (traverseRun < 2)
    {
        if (arrayIndex < 0)
        {
            currCommonSub->add(key, ItemIndex::SHREC);
            currCommonSub->get<ShRec>(key) = std::make_shared<Record>(makeRecordDiffResult());
        }
        else if (arraySize > 0)
        {
            currCommonSub->get<Vec1ShRec>(key).push_back(
                std::make_shared<Record>(makeRecordDiffResult()));
        }
    }

    return true;
}

bool RecordDiff::shouldLoop(const String& key, Size arraySize, bool skipRecordArrays)
{
    // Only Vec1ShRec is handled here, single ShRec entries in shouldDescend().
    if (!isPresent(key, ItemIndex::VEC1SHREC)) return false;
    if (!hasSameType(key, ItemIndex::VEC1SHREC)) return false;

    const auto& itemOther = currOther->vcget<Vec1ShRec>(key);
    Size        itemOtherSize = itemOther.size();
    if (!itemOther.empty() && skipRecordArrays) itemOtherSize = 1;
    if (arraySize != itemOtherSize)
    {
        // First pass has already created diff entry.
        if (traverseRun < 2)
        {
            diffEntry("SIZE1:" + key, static_cast<Int>(arraySize), static_cast<Int>(itemOtherSize));
        }
        // If size is different it makes little sense to continue checking all elements.
        return false;
    }

    if (arraySize > 0 && traverseRun < 2) { currCommonSub->add(key, ItemIndex::VEC1SHREC); }

    return true;
}

vaspml::Record vaspml::rec::detail::makeRecordDiffResult()
{
    Record res(vaspml::makeRecord({
        {"_type",      "String"},
        {"only-in-1",  "ShRec" },
        {"only-in-2",  "ShRec" },
        {"common-sub", "ShRec" },
        {"local-diff", "ShRec" }
    }));

    res.get<String>("_type") = "RecordDiffResult";
    res.get<ShRec>("local-diff")->put<Vec1String>("KEYS");

    return res;
}

bool vaspml::rec::detail::emptyRecordDiffLocal(const Record& diff)
{
    rec::checkType(diff, "RecordDiffResult", flf(VASPML_FLF));

    return diff.cget<ShRec>("only-in-1")->empty() && diff.cget<ShRec>("only-in-2")->empty()
        && diff.cget<ShRec>("local-diff")->size() == 1
        && diff.cget<ShRec>("local-diff")->cget<Vec1String>("KEYS").empty();
}

bool vaspml::rec::detail::emptyRecordDiffRecursive(const Record& diff)
{
    if (!emptyRecordDiffLocal(diff)) return false;

    const Record& common = diff.dcget<ShRec>("common-sub");

    for (const auto& k : common.keys())
    {
        if (common.itemIndexOf(k) == ItemIndex::SHREC)
        {
            if (!emptyRecordDiffRecursive(common.dcget<ShRec>(k))) return false;
        }
        else if (common.itemIndexOf(k) == ItemIndex::VEC1SHREC)
        {
            for (const auto& subDiff : common.vcget<Vec1ShRec>(k))
            {
                if (!emptyRecordDiffRecursive(*subDiff)) return false;
            }
        }
    }

    return true;
}
