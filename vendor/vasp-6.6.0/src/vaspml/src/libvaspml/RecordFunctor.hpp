#ifndef RECORDFUNCTOR_HPP
#define RECORDFUNCTOR_HPP

#include "types.hpp"

namespace vaspml
{
struct Record;
}

namespace vaspml::rec::detail
{

/// Base class for all functors used to traverse() a Record.
struct RecordFunctor
{
    Int    level = 0;
    String key;

    // clang-format off
    virtual void operator()(const Record& record, bool entry);
    virtual void operator()(const String& key, const Real&       item);
    virtual void operator()(const String& key, const Int&        item);
    virtual void operator()(const String& key, const String&     item);
    virtual void operator()(const String& key, const bool&       item);
    virtual void operator()(const String& key, const Vec1Real&   item);
    virtual void operator()(const String& key, const Vec1Int&    item);
    virtual void operator()(const String& key, const Vec1String& item);
    virtual void operator()(const String& key, const Vec2Real&   item);
    virtual void operator()(const String& key, const Vec2Int&    item);
    virtual void operator()(const String& key, const Vec2String& item);
    // clang-format on

    virtual void down(const String& key, Int arrayIndex = -1, Size arraySize = 1);
    virtual void up(const String& key, Int arrayIndex = -1, Size arraySize = 1);
    virtual void next(const String& key, bool last);
    virtual bool shouldLoop(const String& key, Size arraySize, bool skipRecordArrays);
    virtual bool shouldDescend(const Record& record,
                               const String& key,
                               Int           arrayIndex = -1,
                               Size          arraySize = 1);
    virtual void postLoop();
};

/**************************************************************************************************
 * Traverse a Record and apply functor to all Items.
 *
 * @param record Record to traverse (const, will not be modified).
 * @param f A functor which applies its operator() to all Items of the traversed Record.
 * @param skipRecordArrays If true, only traverse the first entry of each Vec1ShRec Item.
 *
 * This function recursively traverses all key-Item pairs of a given Record and applies the
 * functor's operator() to all Items.
 **************************************************************************************************/
void traverse(const Record& record, RecordFunctor& f, bool skipRecordArrays = false);

} // namespace vaspml::rec::detail

#endif
