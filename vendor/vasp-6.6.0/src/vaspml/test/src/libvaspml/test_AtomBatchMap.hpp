#ifndef TEST_ATOMBATCHMAP_HPP
#define TEST_ATOMBATCHMAP_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "BatchMap.hpp"
#include "Structure.hpp"
#include "types.hpp"

namespace vaspml
{

struct TestCase_AtomBatchMap : public TestCase
{
    Vec1String order;
    std::vector<Structure> structures;
    TestCase_AtomBatchMap(String name) : TestCase(name) {}
    AtomBatchMap map;
    Vec1String types;
    Vec1Int nTypes;

    Vec1String initOrder;
    Vec1Int typeOrderInt;
    Vec1String typeOrder;
    Vec2Int strucList;
    Vec2Int atomList;
    Vec1String typeBatchTest;
};

template<>
inline void TestCaseContainer<TestCase_AtomBatchMap>::setup()
{
    TestCase_AtomBatchMap* tc = nullptr;
    testCases.push_back(TestCase_AtomBatchMap("TestCase_AtomBatchMap"));
    tc = &(testCases.back());
    tc->order = { "Pb", "Br", "Cs" };
    tc->types = {"Pb","Br","Cs"};
    tc->nTypes = {2,4,2};
    tc->structures.resize(2);
    tc->structures[0].set_types( tc->types, tc->nTypes );
    tc->structures[1].set_types( tc->types, tc->nTypes );
    tc->map.makeMap( tc->structures, tc->types );

    tc->initOrder = {"Pb","Pb","Br","Br","Br","Br","Cs","Cs",
                     "Pb","Pb","Br","Br","Br","Br","Cs","Cs"};
    tc->typeOrderInt = {0,1,8,9,2,3,4,5,10,11,12,13,6,7,14,15};
    tc->typeOrder    = {"Pb","Pb","Pb","Pb","Br","Br","Br","Br",
                        "Br","Br","Br","Br","Cs","Cs","Cs","Cs"};
    tc->strucList = {{0,1,1},{0,1,1},{0,0,1}};
    tc->atomList  = {{1,2,3},{2,4,6},{3,5,7}};
    tc->typeBatchTest = {"Pb","Br","Br","Br","Br","Cs","Br","Br","Cs"};
}

} // namespace vaspml

#endif
