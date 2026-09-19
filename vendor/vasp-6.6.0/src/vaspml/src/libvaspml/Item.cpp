#include "Item.hpp"

#include "SmartEnum.hpp"

using namespace vaspml;

template<>
SmartEnum<ItemIndex>::EnumMap SmartEnum<ItemIndex>::mapEnumsNames()
{
    return SmartEnum<ItemIndex>::EnumMap{
        {ItemIndex::UNDEFINED,    "UNDEFINED"   },
        {ItemIndex::REAL,         "Real"        },
        {ItemIndex::INT,          "Int"         },
        {ItemIndex::STRING,       "String"      },
        {ItemIndex::BOOL,         "bool"        },
        {ItemIndex::SHREC,        "ShRec"       },
        {ItemIndex::VEC1REAL,     "Vec1Real"    },
        {ItemIndex::VEC1INT,      "Vec1Int"     },
        {ItemIndex::VEC1STRING,   "Vec1String"  },
        {ItemIndex::VEC1SHREC,    "Vec1ShRec"   },
        {ItemIndex::VEC2REAL,     "Vec2Real"    },
        {ItemIndex::VEC2INT,      "Vec2Int"     },
        {ItemIndex::VEC2STRING,   "Vec2String"  },
        {ItemIndex::RECORD_BEGIN, "RECORD_BEGIN"},
        {ItemIndex::RECORD_ITEM,  "RECORD_ITEM" },
        {ItemIndex::RECORD_END,   "RECORD_END"  }
    };
}
