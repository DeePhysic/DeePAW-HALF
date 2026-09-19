#include "SampleStructure.hpp"

#include <stdexcept>

using namespace vaspml;

SampleStructure::SampleStructure(String sample)
{
    setup(sample);
}

void SampleStructure::setup(String sample)
{
    bool       directLattice;
    Vec1Real   lattice;
    Vec1String types;
    Vec1Int    number_atoms;
    Vec1Real   positions;

    // Cesium Lead Bromide structure CsPbBr3 with 40 atoms.
    // Data taken originally from /fsc/home/jona/ForumTests/NNList/mlff_bug/CsPbI3/NNTest
    if (sample == "CsPbBr3_40")
    {
        directLattice = true;

        lattice.resize(9);
        lattice[0] = 11.720599;
        lattice[1] = 0.352943;
        lattice[2] = -0.392053;
        lattice[3] = 0.000003;
        lattice[4] = 11.756821;
        lattice[5] = 0.172043;
        lattice[6] = 0.000032;
        lattice[7] = 0.000023;
        lattice[8] = 11.878922;

        types.resize(3);
        types[0] = "Pb";
        types[1] = "Br";
        types[2] = "Cs";

        number_atoms.resize(3);
        number_atoms[0] = 8;
        number_atoms[1] = 24;
        number_atoms[2] = 8;

        positions.resize(3 * 40);
        positions[0] = 0.56551602;
        positions[1] = 0.87156268;
        positions[2] = 0.21391052;
        positions[3] = 0.57021909;
        positions[4] = 0.88082219;
        positions[5] = 0.71694288;
        positions[6] = 0.06881313;
        positions[7] = 0.37553489;
        positions[8] = 0.20907002;
        positions[9] = 0.06458391;
        positions[10] = 0.38752891;
        positions[11] = 0.70536099;
        positions[12] = 0.06752362;
        positions[13] = 0.88197489;
        positions[14] = 0.21262166;
        positions[15] = 0.06689049;
        positions[16] = 0.88915338;
        positions[17] = 0.70549603;
        positions[18] = 0.56757310;
        positions[19] = 0.37342848;
        positions[20] = 0.21303976;
        positions[21] = 0.56426728;
        positions[22] = 0.38436588;
        positions[23] = 0.71898604;
        positions[24] = 0.57977371;
        positions[25] = 0.83254086;
        positions[26] = 0.45552409;
        positions[27] = 0.03564442;
        positions[28] = 0.41185722;
        positions[29] = 0.96396452;
        positions[30] = 0.10877164;
        positions[31] = 0.33948173;
        positions[32] = 0.46046055;
        positions[33] = 0.52963781;
        positions[34] = 0.92296378;
        positions[35] = 0.96299987;
        positions[36] = 0.01964462;
        positions[37] = 0.92660674;
        positions[38] = 0.46011073;
        positions[39] = 0.61376908;
        positions[40] = 0.34438977;
        positions[41] = 0.96064813;
        positions[42] = 0.53704924;
        positions[43] = 0.40773114;
        positions[44] = 0.46345572;
        positions[45] = 0.09097253;
        positions[46] = 0.85453977;
        positions[47] = 0.95725728;
        positions[48] = 0.31179584;
        positions[49] = 0.92267176;
        positions[50] = 0.23889681;
        positions[51] = 0.31559319;
        positions[52] = 0.92425459;
        positions[53] = 0.67975043;
        positions[54] = 0.81784281;
        positions[55] = 0.41647750;
        positions[56] = 0.24737584;
        positions[57] = 0.81605860;
        positions[58] = 0.43569019;
        positions[59] = 0.68109170;
        positions[60] = 0.32069238;
        positions[61] = 0.32941746;
        positions[62] = 0.17646188;
        positions[63] = 0.31935606;
        positions[64] = 0.33038226;
        positions[65] = 0.73773697;
        positions[66] = 0.81695813;
        positions[67] = 0.82978129;
        positions[68] = 0.17012218;
        positions[69] = 0.81680048;
        positions[70] = 0.83950003;
        positions[71] = 0.74181736;
        positions[72] = 0.03004165;
        positions[73] = 0.12410075;
        positions[74] = 0.16184966;
        positions[75] = 0.01369403;
        positions[76] = 0.13570767;
        positions[77] = 0.74181787;
        positions[78] = 0.51736305;
        positions[79] = 0.62070171;
        positions[80] = 0.17174387;
        positions[81] = 0.52728509;
        positions[82] = 0.63106244;
        positions[83] = 0.76381618;
        positions[84] = 0.61453648;
        positions[85] = 0.12070975;
        positions[86] = 0.24158969;
        positions[87] = 0.60363380;
        positions[88] = 0.12560660;
        positions[89] = 0.68056718;
        positions[90] = 0.10357399;
        positions[91] = 0.62786407;
        positions[92] = 0.24182939;
        positions[93] = 0.11474066;
        positions[94] = 0.63539900;
        positions[95] = 0.69119129;
        positions[96] = 0.79062716;
        positions[97] = 0.08312205;
        positions[98] = 0.96560138;
        positions[99] = 0.82918114;
        positions[100] = 0.18111056;
        positions[101] = 0.45622823;
        positions[102] = 0.28171489;
        positions[103] = 0.58655477;
        positions[104] = 0.97009963;
        positions[105] = 0.33581972;
        positions[106] = 0.65670823;
        positions[107] = 0.46128160;
        positions[108] = 0.26595342;
        positions[109] = 0.10685946;
        positions[110] = 0.94409577;
        positions[111] = 0.36668353;
        positions[112] = 0.13721076;
        positions[113] = 0.46565333;
        positions[114] = 0.77217197;
        positions[115] = 0.61411500;
        positions[116] = 0.95719585;
        positions[117] = 0.85184039;
        positions[118] = 0.66571424;
        positions[119] = 0.46245759;
    }
    // Calcium oxide CaO structure with 16 atoms (bulk).
    else if (sample == "CaO_16")
    {
        directLattice = true;

        lattice.resize(9);
        lattice[0] = 4.84190736;
        lattice[1] = 4.84190736;
        lattice[2] = 0.00000000;
        lattice[3] = 0.00000000;
        lattice[4] = 4.84190736;
        lattice[5] = 4.84190736;
        lattice[6] = 4.84190736;
        lattice[7] = 0.00000000;
        lattice[8] = 4.84190736;

        types.resize(2);
        types[0] = "Ca";
        types[1] = "O";

        number_atoms.resize(2);
        number_atoms[0] = 8;
        number_atoms[1] = 8;

        positions.resize(3 * 16);
        positions[0] = 0.00000000;
        positions[1] = 0.00000000;
        positions[2] = 0.00000000;
        positions[3] = 0.50000000;
        positions[4] = 0.00000000;
        positions[5] = 0.00000000;
        positions[6] = 0.00000000;
        positions[7] = 0.50000000;
        positions[8] = 0.00000000;
        positions[9] = 0.50000000;
        positions[10] = 0.50000000;
        positions[11] = 0.00000000;
        positions[12] = 0.00000000;
        positions[13] = 0.00000000;
        positions[14] = 0.50000000;
        positions[15] = 0.50000000;
        positions[16] = 0.00000000;
        positions[17] = 0.50000000;
        positions[18] = 0.00000000;
        positions[19] = 0.50000000;
        positions[20] = 0.50000000;
        positions[21] = 0.50000000;
        positions[22] = 0.50000000;
        positions[23] = 0.50000000;
        positions[24] = 0.25000000;
        positions[25] = 0.25000000;
        positions[26] = 0.25000000;
        positions[27] = 0.75000000;
        positions[28] = 0.25000000;
        positions[29] = 0.25000000;
        positions[30] = 0.25000000;
        positions[31] = 0.75000000;
        positions[32] = 0.25000000;
        positions[33] = 0.75000000;
        positions[34] = 0.75000000;
        positions[35] = 0.25000000;
        positions[36] = 0.25000000;
        positions[37] = 0.25000000;
        positions[38] = 0.75000000;
        positions[39] = 0.75000000;
        positions[40] = 0.25000000;
        positions[41] = 0.75000000;
        positions[42] = 0.25000000;
        positions[43] = 0.75000000;
        positions[44] = 0.75000000;
        positions[45] = 0.75000000;
        positions[46] = 0.75000000;
        positions[47] = 0.75000000;
    }
    // Azobenzene C12H10N2 with 24 atoms (molecule).
    else if (sample == "HCN_24")
    {

        directLattice = false;

        lattice.resize(9);
        lattice[0] = 20.0000000;
        lattice[1] = 0.00000000;
        lattice[2] = 0.00000000;
        lattice[3] = 0.00000000;
        lattice[4] = 20.0000000;
        lattice[5] = 0.00000000;
        lattice[6] = 0.00000000;
        lattice[7] = 0.00000000;
        lattice[8] = 20.0000000;

        types.resize(3);
        types[0] = "H";
        types[1] = "C";
        types[2] = "N";

        number_atoms.resize(3);
        number_atoms[0] = 10;
        number_atoms[1] = 12;
        number_atoms[2] = 2;

        positions.resize(3 * 24);
        positions[0] = 2.313030180000;
        positions[1] = 0.715871250000;
        positions[2] = 2.271780990000;
        positions[3] = 5.582469200000;
        positions[4] = -0.706829660000;
        positions[5] = -0.271662220000;
        positions[6] = 3.982078020000;
        positions[7] = -0.611727320000;
        positions[8] = -2.359525930000;
        positions[9] = 1.584390280000;
        positions[10] = -0.253462150000;
        positions[11] = -2.061452000000;
        positions[12] = -2.099752190000;
        positions[13] = -1.925564420000;
        positions[14] = -1.074219680000;
        positions[15] = -4.716584560000;
        positions[16] = -1.702503910000;
        positions[17] = -0.468009660000;
        positions[18] = -5.577836410000;
        positions[19] = 0.378325580000;
        positions[20] = 0.230899760000;
        positions[21] = -4.009345050000;
        positions[22] = 2.365335710000;
        positions[23] = 0.746591170000;
        positions[24] = -1.546996870000;
        positions[25] = 1.894044660000;
        positions[26] = 0.810782410000;
        positions[27] = 4.705817060000;
        positions[28] = 0.236224650000;
        positions[29] = 1.811314790000;
        positions[30] = 2.652579080000;
        positions[31] = 0.308396760000;
        positions[32] = 1.254638730000;
        positions[33] = 3.981985820000;
        positions[34] = 0.061995520000;
        positions[35] = 1.066697640000;
        positions[36] = 4.483371900000;
        positions[37] = -0.408320760000;
        positions[38] = -0.178261560000;
        positions[39] = 2.269070140000;
        positions[40] = -0.308899040000;
        positions[41] = -1.175242670000;
        positions[42] = 1.778310620000;
        positions[43] = 0.041513910000;
        positions[44] = 0.189297230000;
        positions[45] = -2.650534190000;
        positions[46] = -1.005804020000;
        positions[47] = -0.603692940000;
        positions[48] = -4.016280040000;
        positions[49] = -0.913193310000;
        positions[50] = -0.144731030000;
        positions[51] = -4.502558330000;
        positions[52] = 0.281120000000;
        positions[53] = 0.287105040000;
        positions[54] = -3.651175970000;
        positions[55] = 1.370910890000;
        positions[56] = 0.386281140000;
        positions[57] = -2.226685480000;
        positions[58] = 1.122910680000;
        positions[59] = 0.460748770000;
        positions[60] = -1.786987870000;
        positions[61] = -0.006568090000;
        positions[62] = -0.229997490000;
        positions[63] = 3.610819030000;
        positions[64] = -0.503277730000;
        positions[65] = -1.315586850000;
        positions[66] = -0.351746430000;
        positions[67] = -0.182240950000;
        positions[68] = -0.453871580000;
        positions[69] = 0.385919690000;
        positions[70] = 0.119222680000;
        positions[71] = 0.482413680000;
    }
    else { throw std::invalid_argument("ERROR: Unknown sample structure: \"" + sample + "\"."); }

    // Finally set up the structure.
    importPoscarData(lattice, types, number_atoms, positions, directLattice);

    return;
}
