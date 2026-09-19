#ifndef TEST_DESCRIPTORSHS3REDUCEDLINELEM_HPP
#define TEST_DESCRIPTORSHS3REDUCEDLINELEM_HPP

#include "SampleNeighborList.hpp"
#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "BasisFunctions.hpp"
#include "BasisFunctionsAngular.hpp"
#include "TypeMap.hpp"
#include "cutoff.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <limits>

using CT = vaspml::math::CutoffType;

namespace vaspml
{

// this crazy number belongs to a check 1E-6
Real const defaultTolerance = 100000000000.0 * std::numeric_limits<Real>::epsilon();

struct TestCase_DescriptorSHS3ReducedLinElem_CsPbBr3 : public TestCase
{
    Real                                   tolerance;
    Real                                   cutoff;
    Vec2Int                                descriptorListMLFF;
    Vec1Real                               targetSHS3;
    TypeMap                                typeMap;
    std::shared_ptr<SampleNeighborList>    nn_list;
    std::shared_ptr<BasisFunctionsAngular> descriptor;

    TestCase_DescriptorSHS3ReducedLinElem_CsPbBr3(String name) :
        TestCase(name),
        tolerance(defaultTolerance)
    {}
};

template<>
void TestCaseContainer<TestCase_DescriptorSHS3ReducedLinElem_CsPbBr3>::setup()
{
    String                                         structureName = "";
    TestCase_DescriptorSHS3ReducedLinElem_CsPbBr3* tc = nullptr;

    //////////////////////////
    // Structure CsPbBr3_40 //
    //////////////////////////
    structureName = "CsPbBr3_40";
    testCases.push_back(TestCase_DescriptorSHS3ReducedLinElem_CsPbBr3(
        "TestCase_descriptorSHS3ReducedLinElem_" + structureName));

    tc = &(testCases.back());

    tc->nn_list = std::make_shared<SampleNeighborList>(5.0, true, false, structureName);

    tc->descriptor = std::make_shared<BasisFunctionsAngular>(CT::BP,
                                                             5.0,
                                                             0.5,
                                                             1000,
                                                             2,
                                                             4,
                                                             BasisFunctionType::bodyOrder3);

    Vec1String ff_types = {"Pb", "Br", "Cs"};
    Vec1String struc_types = {"Pb", "Br", "Cs"};

    tc->typeMap.update(ff_types, struc_types);

    vector_tools::allocate_vector(tc->descriptorListMLFF, 3, 66);

    // here is the list of active descriptors
    tc->descriptorListMLFF[0][0] = 1;
    tc->descriptorListMLFF[0][1] = 2;
    tc->descriptorListMLFF[0][2] = 3;
    tc->descriptorListMLFF[0][3] = 4;
    tc->descriptorListMLFF[0][4] = 5;
    tc->descriptorListMLFF[0][5] = 6;
    tc->descriptorListMLFF[0][6] = 7;
    tc->descriptorListMLFF[0][7] = 8;
    tc->descriptorListMLFF[0][8] = 9;
    tc->descriptorListMLFF[0][9] = 10;
    tc->descriptorListMLFF[0][10] = 11;
    tc->descriptorListMLFF[0][11] = 12;
    tc->descriptorListMLFF[0][12] = 13;
    tc->descriptorListMLFF[0][13] = 14;
    tc->descriptorListMLFF[0][14] = 15;
    tc->descriptorListMLFF[0][15] = 16;
    tc->descriptorListMLFF[0][16] = 17;
    tc->descriptorListMLFF[0][17] = 18;
    tc->descriptorListMLFF[0][18] = 19;
    tc->descriptorListMLFF[0][19] = 20;
    tc->descriptorListMLFF[0][20] = 21;
    tc->descriptorListMLFF[0][21] = 22;
    tc->descriptorListMLFF[0][22] = 23;
    tc->descriptorListMLFF[0][23] = 24;
    tc->descriptorListMLFF[0][24] = 25;
    tc->descriptorListMLFF[0][25] = 26;
    tc->descriptorListMLFF[0][26] = 27;
    tc->descriptorListMLFF[0][27] = 28;
    tc->descriptorListMLFF[0][28] = 29;
    tc->descriptorListMLFF[0][29] = 30;
    tc->descriptorListMLFF[0][30] = 31;
    tc->descriptorListMLFF[0][31] = 32;
    tc->descriptorListMLFF[0][32] = 33;
    tc->descriptorListMLFF[0][33] = 34;
    tc->descriptorListMLFF[0][34] = 35;
    tc->descriptorListMLFF[0][35] = 36;
    tc->descriptorListMLFF[0][36] = 37;
    tc->descriptorListMLFF[0][37] = 38;
    tc->descriptorListMLFF[0][38] = 39;
    tc->descriptorListMLFF[0][39] = 40;
    tc->descriptorListMLFF[0][40] = 41;
    tc->descriptorListMLFF[0][41] = 42;
    tc->descriptorListMLFF[0][42] = 43;
    tc->descriptorListMLFF[0][43] = 44;
    tc->descriptorListMLFF[0][44] = 45;
    tc->descriptorListMLFF[0][45] = 46;
    tc->descriptorListMLFF[0][46] = 47;
    tc->descriptorListMLFF[0][47] = 48;
    tc->descriptorListMLFF[0][48] = 49;
    tc->descriptorListMLFF[0][49] = 50;
    tc->descriptorListMLFF[0][50] = 51;
    tc->descriptorListMLFF[0][51] = 52;
    tc->descriptorListMLFF[0][52] = 53;
    tc->descriptorListMLFF[0][53] = 54;
    tc->descriptorListMLFF[0][54] = 55;
    tc->descriptorListMLFF[0][55] = 56;
    tc->descriptorListMLFF[0][56] = 57;
    tc->descriptorListMLFF[0][57] = 58;
    tc->descriptorListMLFF[0][58] = 59;
    tc->descriptorListMLFF[0][59] = 60;
    tc->descriptorListMLFF[0][60] = 61;
    tc->descriptorListMLFF[0][61] = 62;
    tc->descriptorListMLFF[0][62] = 63;
    tc->descriptorListMLFF[0][63] = 64;
    tc->descriptorListMLFF[0][64] = 65;
    tc->descriptorListMLFF[0][65] = 66;
    tc->descriptorListMLFF[1][0] = 1;
    tc->descriptorListMLFF[1][1] = 2;
    tc->descriptorListMLFF[1][2] = 3;
    tc->descriptorListMLFF[1][3] = 4;
    tc->descriptorListMLFF[1][4] = 5;
    tc->descriptorListMLFF[1][5] = 6;
    tc->descriptorListMLFF[1][6] = 7;
    tc->descriptorListMLFF[1][7] = 8;
    tc->descriptorListMLFF[1][8] = 9;
    tc->descriptorListMLFF[1][9] = 10;
    tc->descriptorListMLFF[1][10] = 11;
    tc->descriptorListMLFF[1][11] = 12;
    tc->descriptorListMLFF[1][12] = 13;
    tc->descriptorListMLFF[1][13] = 14;
    tc->descriptorListMLFF[1][14] = 15;
    tc->descriptorListMLFF[1][15] = 16;
    tc->descriptorListMLFF[1][16] = 17;
    tc->descriptorListMLFF[1][17] = 18;
    tc->descriptorListMLFF[1][18] = 19;
    tc->descriptorListMLFF[1][19] = 20;
    tc->descriptorListMLFF[1][20] = 21;
    tc->descriptorListMLFF[1][21] = 22;
    tc->descriptorListMLFF[1][22] = 23;
    tc->descriptorListMLFF[1][23] = 24;
    tc->descriptorListMLFF[1][24] = 25;
    tc->descriptorListMLFF[1][25] = 26;
    tc->descriptorListMLFF[1][26] = 27;
    tc->descriptorListMLFF[1][27] = 28;
    tc->descriptorListMLFF[1][28] = 29;
    tc->descriptorListMLFF[1][29] = 30;
    tc->descriptorListMLFF[1][30] = 31;
    tc->descriptorListMLFF[1][31] = 32;
    tc->descriptorListMLFF[1][32] = 33;
    tc->descriptorListMLFF[1][33] = 34;
    tc->descriptorListMLFF[1][34] = 35;
    tc->descriptorListMLFF[1][35] = 36;
    tc->descriptorListMLFF[1][36] = 37;
    tc->descriptorListMLFF[1][37] = 38;
    tc->descriptorListMLFF[1][38] = 39;
    tc->descriptorListMLFF[1][39] = 40;
    tc->descriptorListMLFF[1][40] = 41;
    tc->descriptorListMLFF[1][41] = 42;
    tc->descriptorListMLFF[1][42] = 43;
    tc->descriptorListMLFF[1][43] = 44;
    tc->descriptorListMLFF[1][44] = 45;
    tc->descriptorListMLFF[1][45] = 46;
    tc->descriptorListMLFF[1][46] = 47;
    tc->descriptorListMLFF[1][47] = 48;
    tc->descriptorListMLFF[1][48] = 49;
    tc->descriptorListMLFF[1][49] = 50;
    tc->descriptorListMLFF[1][50] = 51;
    tc->descriptorListMLFF[1][51] = 52;
    tc->descriptorListMLFF[1][52] = 53;
    tc->descriptorListMLFF[1][53] = 54;
    tc->descriptorListMLFF[1][54] = 55;
    tc->descriptorListMLFF[1][55] = 56;
    tc->descriptorListMLFF[1][56] = 57;
    tc->descriptorListMLFF[1][57] = 58;
    tc->descriptorListMLFF[1][58] = 59;
    tc->descriptorListMLFF[1][59] = 60;
    tc->descriptorListMLFF[1][60] = 61;
    tc->descriptorListMLFF[1][61] = 62;
    tc->descriptorListMLFF[1][62] = 63;
    tc->descriptorListMLFF[1][63] = 64;
    tc->descriptorListMLFF[1][64] = 65;
    tc->descriptorListMLFF[1][65] = 66;
    tc->descriptorListMLFF[2][0] = 1;
    tc->descriptorListMLFF[2][1] = 2;
    tc->descriptorListMLFF[2][2] = 3;
    tc->descriptorListMLFF[2][3] = 4;
    tc->descriptorListMLFF[2][4] = 5;
    tc->descriptorListMLFF[2][5] = 6;
    tc->descriptorListMLFF[2][6] = 7;
    tc->descriptorListMLFF[2][7] = 8;
    tc->descriptorListMLFF[2][8] = 9;
    tc->descriptorListMLFF[2][9] = 10;
    tc->descriptorListMLFF[2][10] = 11;
    tc->descriptorListMLFF[2][11] = 12;
    tc->descriptorListMLFF[2][12] = 13;
    tc->descriptorListMLFF[2][13] = 14;
    tc->descriptorListMLFF[2][14] = 15;
    tc->descriptorListMLFF[2][15] = 16;
    tc->descriptorListMLFF[2][16] = 17;
    tc->descriptorListMLFF[2][17] = 18;
    tc->descriptorListMLFF[2][18] = 19;
    tc->descriptorListMLFF[2][19] = 20;
    tc->descriptorListMLFF[2][20] = 21;
    tc->descriptorListMLFF[2][21] = 22;
    tc->descriptorListMLFF[2][22] = 23;
    tc->descriptorListMLFF[2][23] = 24;
    tc->descriptorListMLFF[2][24] = 25;
    tc->descriptorListMLFF[2][25] = 26;
    tc->descriptorListMLFF[2][26] = 27;
    tc->descriptorListMLFF[2][27] = 28;
    tc->descriptorListMLFF[2][28] = 29;
    tc->descriptorListMLFF[2][29] = 30;
    tc->descriptorListMLFF[2][30] = 31;
    tc->descriptorListMLFF[2][31] = 32;
    tc->descriptorListMLFF[2][32] = 33;
    tc->descriptorListMLFF[2][33] = 34;
    tc->descriptorListMLFF[2][34] = 35;
    tc->descriptorListMLFF[2][35] = 36;
    tc->descriptorListMLFF[2][36] = 37;
    tc->descriptorListMLFF[2][37] = 38;
    tc->descriptorListMLFF[2][38] = 39;
    tc->descriptorListMLFF[2][39] = 40;
    tc->descriptorListMLFF[2][40] = 41;
    tc->descriptorListMLFF[2][41] = 42;
    tc->descriptorListMLFF[2][42] = 43;
    tc->descriptorListMLFF[2][43] = 44;
    tc->descriptorListMLFF[2][44] = 45;
    tc->descriptorListMLFF[2][45] = 46;
    tc->descriptorListMLFF[2][46] = 47;
    tc->descriptorListMLFF[2][47] = 48;
    tc->descriptorListMLFF[2][48] = 49;
    tc->descriptorListMLFF[2][49] = 50;
    tc->descriptorListMLFF[2][50] = 51;
    tc->descriptorListMLFF[2][51] = 52;
    tc->descriptorListMLFF[2][52] = 53;
    tc->descriptorListMLFF[2][53] = 54;
    tc->descriptorListMLFF[2][54] = 55;
    tc->descriptorListMLFF[2][55] = 56;
    tc->descriptorListMLFF[2][56] = 57;
    tc->descriptorListMLFF[2][57] = 58;
    tc->descriptorListMLFF[2][58] = 59;
    tc->descriptorListMLFF[2][59] = 60;
    tc->descriptorListMLFF[2][60] = 61;
    tc->descriptorListMLFF[2][61] = 62;
    tc->descriptorListMLFF[2][62] = 63;
    tc->descriptorListMLFF[2][63] = 64;
    tc->descriptorListMLFF[2][64] = 65;
    tc->descriptorListMLFF[2][65] = 66;

    // descriptors for first atom
    tc->targetSHS3.resize(66);
    tc->targetSHS3[0] = 0.0;
    tc->targetSHS3[1] = 0.0;
    tc->targetSHS3[2] = 0.0;
    tc->targetSHS3[3] = 0.0;
    tc->targetSHS3[4] = 0.0;
    tc->targetSHS3[5] = 0.0;
    tc->targetSHS3[6] = 0.0;
    tc->targetSHS3[7] = 0.0;
    tc->targetSHS3[8] = 0.0;
    tc->targetSHS3[9] = 0.0;
    tc->targetSHS3[10] = 0.0;
    tc->targetSHS3[11] = 0.0;
    tc->targetSHS3[12] = 0.0;
    tc->targetSHS3[13] = 0.0;
    tc->targetSHS3[14] = 0.0;
    tc->targetSHS3[15] = 0.0;
    tc->targetSHS3[16] = 0.0;
    tc->targetSHS3[17] = 0.0;
    tc->targetSHS3[18] = 0.0;
    tc->targetSHS3[19] = 0.0;
    tc->targetSHS3[20] = 0.0;
    tc->targetSHS3[21] = 0.0;
    tc->targetSHS3[22] = 0.10526038087255735;
    tc->targetSHS3[23] = -8.1678037880686075E-002;
    tc->targetSHS3[24] = -5.9054584224806965E-002;
    tc->targetSHS3[25] = 7.1077124761860844E-002;
    tc->targetSHS3[26] = 3.1637357569606353E-002;
    tc->targetSHS3[27] = 3.2349198191625843E-002;
    tc->targetSHS3[28] = -3.8934962052387036E-002;
    tc->targetSHS3[29] = 1.6644776424954939E-002;
    tc->targetSHS3[30] = -2.8331475877888843E-002;
    tc->targetSHS3[31] = 2.4095789071622839E-002;
    tc->targetSHS3[32] = 1.0306955241880938E-004;
    tc->targetSHS3[33] = 9.2453320870626460E-005;
    tc->targetSHS3[34] = -1.2245300688135713E-004;
    tc->targetSHS3[35] = 2.0669939711720611E-004;
    tc->targetSHS3[36] = -1.6855856344839414E-004;
    tc->targetSHS3[37] = 9.7495183646042046E-005;
    tc->targetSHS3[38] = 3.0487567762586338E-004;
    tc->targetSHS3[39] = 3.1638482211004314E-004;
    tc->targetSHS3[40] = -2.0607344501230704E-004;
    tc->targetSHS3[41] = 3.0710986131121312E-004;
    tc->targetSHS3[42] = -1.6682011222998136E-004;
    tc->targetSHS3[43] = 7.5486697950884047E-005;
    tc->targetSHS3[44] = 7.8713243955954928E-005;
    tc->targetSHS3[45] = -6.1078472909291236E-005;
    tc->targetSHS3[46] = -4.4160755036908725E-005;
    tc->targetSHS3[47] = 5.3151157298603376E-005;
    tc->targetSHS3[48] = 7.5859628163443292E-005;
    tc->targetSHS3[49] = 7.7566469981039027E-005;
    tc->targetSHS3[50] = -9.3357725510215755E-005;
    tc->targetSHS3[51] = -6.6592952399812432E-005;
    tc->targetSHS3[52] = 1.1334947231397213E-004;
    tc->targetSHS3[53] = -8.0411899594322455E-005;
    tc->targetSHS3[54] = -1.3354024449713259E-007;
    tc->targetSHS3[55] = -9.8411308277933071E-007;
    tc->targetSHS3[56] = 4.6928687196384564E-007;
    tc->targetSHS3[57] = 9.7835721430687626E-007;
    tc->targetSHS3[58] = -6.5883617382773794E-007;
    tc->targetSHS3[59] = 4.8792254707784505E-007;
    tc->targetSHS3[60] = 6.9979915467594864E-006;
    tc->targetSHS3[61] = 3.2170192681685378E-006;
    tc->targetSHS3[62] = -4.3914424247641056E-006;
    tc->targetSHS3[63] = -2.8925469840164312E-006;
    tc->targetSHS3[64] = 5.5841575326411814E-006;
    tc->targetSHS3[65] = -4.0018619936387622E-006;
}

} //namespace vaspml
#endif
