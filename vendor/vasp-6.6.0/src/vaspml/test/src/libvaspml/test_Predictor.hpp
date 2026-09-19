#ifndef TEST_PREDICTOR_HPP
#define TEST_PREDICTOR_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "Record.hpp"
#include "io.hpp"
#include "types.hpp"

namespace vaspml
{

struct TestCase_Predictor : public TestCase
{
    Real     toleranceEnergy;
    Real     toleranceStress;
    Real     toleranceForces;
    String   structure;
    Record   ff;
    Real     energy;
    Vec1Real stress;
    Vec1Real forces;

    TestCase_Predictor(String name) : TestCase(name) {}
};

template<>
inline void TestCaseContainer<TestCase_Predictor>::setup()
{
    String              structureName = "";
    String              version = "";
    TestCase_Predictor* tc = nullptr;

    structureName = "MAPbI3.cubic";
    version = "6.5.1";
    testCases.push_back(TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff, "../../data/ff/" + version + "/refit/ML_FF.MAPbI3");
    // Last modifications: AMD AOCC 5.0
    // GNU, PALGO = off    : 7.0E-12
    // GNU, PALGO = serial : 4.0E-11
    // Intel Classic 2022  : 5.0E-11
    tc->toleranceEnergy = 9.0E-11;
    // GNU, PALGO = off    : 5.0E-10
    // GNU, PALGO = serial : 8.0E-10
    tc->toleranceStress = 2.0E-09;
    // Intel Classic 2022  : 5.0E-11
    // Intel oneAPI 2025.3 : 6.0E-11
    tc->toleranceForces = 7.0E-11;

    // clang-format off
    tc->energy = -4.9588457435512410E+01;

    tc->stress.push_back( 1.8549046846740708E+01); // xx
    tc->stress.push_back( 1.1276432332078421E-04); // xy
    tc->stress.push_back(-5.2292265908173807E-01); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 1.9337960571939657E+01); // yy
    tc->stress.push_back( 2.9044881615321820E-04); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back( 1.2448032898765842E+01); // zz

    tc->forces.insert(tc->forces.end(), { 3.8487326963378887E-01, -4.9245506792619359E-05,  6.5970680566089610E-02});
    tc->forces.insert(tc->forces.end(), {-1.7468865170210496E-01,  2.4651115423396056E-05, -9.2520098462522668E-02});
    tc->forces.insert(tc->forces.end(), {-1.6034706906042953E-01,  1.3121768102383977E-05,  2.1950178374564924E-01});
    tc->forces.insert(tc->forces.end(), { 6.4750411586353473E-02,  4.8652469945159689E-06, -3.2914215384931689E-01});
    tc->forces.insert(tc->forces.end(), { 2.0094017500315081E-01, -3.2135925220159706E-05,  2.3709409575287105E-01});
    tc->forces.insert(tc->forces.end(), {-5.0960333926956269E-01,  6.6993794873849711E-05, -1.5979888661051483E-01});
    tc->forces.insert(tc->forces.end(), { 7.6905318374988249E-02,  6.4834495742199364E-06, -5.1830636069934244E-01});
    tc->forces.insert(tc->forces.end(), {-2.9637987963817103E-01,  3.2268017271165050E-01,  1.5497707681424694E-01});
    tc->forces.insert(tc->forces.end(), {-2.9645820769621589E-01, -3.2261751535026667E-01,  1.5495701518021457E-01});
    tc->forces.insert(tc->forces.end(), { 3.8725131088039788E-01,  6.2519469157954677E-01, -1.0584372427671189E-01});
    tc->forces.insert(tc->forces.end(), { 3.8709882478301322E-01, -6.2528352435367485E-01, -1.0588256831510966E-01});
    tc->forces.insert(tc->forces.end(), {-6.4342162895208416E-02, -8.5585302113069409E-06,  4.7899314015444705E-01});
    // clang-format on


    structureName = "MAPbI3.orthorhombic";
    version = "6.5.1";
    testCases.push_back(TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff, "../../data/ff/" + version + "/refit/ML_FF.MAPbI3");
    // Last modifications: GNU 11.2 + MKL 2023.2
    // Intel oneAPI 2025.3: 3.0E-10
    tc->toleranceEnergy = 4.0E-10;
    // GNU                : 2.0E-10
    // Intel Classic 2022 : 4.0E-10
    // AMD AOCC 5.0       : 5.0E-10
    tc->toleranceStress = 2.0E-09;
    // GNU                : 8.0E-11
    // Intel Classic 2022 : 1.0E-10
    tc->toleranceForces = 2.0E-10;

    // clang-format off
    tc->energy = -1.9688725929424024E+02;

    tc->stress.push_back( 4.6532990804609369E+01); // xx
    tc->stress.push_back( 2.3773582703824548E-06); // xy
    tc->stress.push_back( 1.2253889544379973E-04); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 3.8660259828483952E+01); // yy
    tc->stress.push_back( 2.2004050007720416E-11); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back( 2.7838815619349678E+01); // zz

    tc->forces.insert(tc->forces.end(), {-7.4815280881438428E-08, -2.1294818316445688E-07, -3.7830851297479162E-07});
    tc->forces.insert(tc->forces.end(), {-1.7289421437310816E-07,  2.7884185572547216E-07,  7.8935045897048237E-07});
    tc->forces.insert(tc->forces.end(), {-4.9562488423366616E-08,  2.0781230255982659E-07, -3.4306548727435207E-07});
    tc->forces.insert(tc->forces.end(), {-1.4764034961059388E-07, -2.8397792401726694E-07,  7.5410745487263923E-07});
    tc->forces.insert(tc->forces.end(), { 1.4078462048594467E-01, -4.5774154048402567E-09,  5.2836428656138867E-01});
    tc->forces.insert(tc->forces.end(), {-1.4078642774547948E-01,  3.0973206233588428E-08, -5.2836416717196899E-01});
    tc->forces.insert(tc->forces.end(), {-1.4078671474378360E-01,  3.0971485978701167E-08,  5.2836405264380271E-01});
    tc->forces.insert(tc->forces.end(), { 1.4078679017828397E-01, -4.5755926906230309E-09, -5.2836233407327216E-01});
    tc->forces.insert(tc->forces.end(), {-1.2483890444367436E-01,  6.3739960488781403E-02, -2.8563724678576619E-01});
    tc->forces.insert(tc->forces.end(), { 1.2483863630420759E-01, -6.3743700774258408E-02,  2.8563732875617470E-01});
    tc->forces.insert(tc->forces.end(), { 1.2483570385866578E-01, -6.3740962132276577E-02, -2.8564045473756683E-01});
    tc->forces.insert(tc->forces.end(), {-1.2483659111867250E-01,  6.3740751166354900E-02,  2.8563863848721904E-01});
    tc->forces.insert(tc->forces.end(), { 1.2483887713293330E-01,  6.3744098570625934E-02,  2.8563721803453729E-01});
    tc->forces.insert(tc->forces.end(), {-1.2483895812240908E-01, -6.3740033528514387E-02, -2.8563721384759222E-01});
    tc->forces.insert(tc->forces.end(), {-1.2483664479840763E-01, -6.3740824198951787E-02,  2.8563860555118131E-01});
    tc->forces.insert(tc->forces.end(), { 1.2483594468734231E-01,  6.3741359929410601E-02, -2.8564034401577926E-01});
    tc->forces.insert(tc->forces.end(), { 3.2166236779104401E-01,  3.6816913742416705E-06,  4.5057213809615138E-01});
    tc->forces.insert(tc->forces.end(), {-3.2166309504384433E-01, -5.3594085221747403E-08, -4.5058145125624294E-01});
    tc->forces.insert(tc->forces.end(), {-3.2166348707874048E-01, -5.3578961931138168E-08,  4.5058162102836535E-01});
    tc->forces.insert(tc->forces.end(), { 3.2162295846074929E-01,  3.6816229749338116E-06, -4.5040689149695112E-01});
    tc->forces.insert(tc->forces.end(), { 4.9860696545711913E-01,  3.3805550172952846E-06, -6.4454965655060059E-01});
    tc->forces.insert(tc->forces.end(), {-4.9860021410214023E-01,  5.5859603280344137E-08,  6.4456640364530249E-01});
    tc->forces.insert(tc->forces.end(), {-4.9859942637945109E-01,  5.5883537012743684E-08, -6.4456620328916125E-01});
    tc->forces.insert(tc->forces.end(), { 4.9858079519964454E-01,  3.3805581959362492E-06,  6.4456299272838358E-01});
    tc->forces.insert(tc->forces.end(), {-1.1272756172228304E+00, -4.5303448515707756E-01, -7.6849005926838032E-01});
    tc->forces.insert(tc->forces.end(), { 1.1272700187700255E+00,  4.5302376628908747E-01,  7.6848499679964855E-01});
    tc->forces.insert(tc->forces.end(), { 1.1272695185475987E+00,  4.5302408309414505E-01, -7.6848431663970462E-01});
    tc->forces.insert(tc->forces.end(), {-1.1272754728566154E+00, -4.5302388627602236E-01,  7.6849099285172517E-01});
    tc->forces.insert(tc->forces.end(), { 1.1272700036877050E+00, -4.5302372058790091E-01,  7.6848503091189757E-01});
    tc->forces.insert(tc->forces.end(), {-1.1272743195442785E+00,  4.5303070225436509E-01, -7.6848776943673691E-01});
    tc->forces.insert(tc->forces.end(), {-1.1272741751708943E+00,  4.5302010336741128E-01,  7.6848870303016370E-01});
    tc->forces.insert(tc->forces.end(), { 1.1272695034773306E+00, -4.5302403741268271E-01, -7.6848435076226962E-01});
    tc->forces.insert(tc->forces.end(), { 1.8318888406499685E+00, -1.7208878396705230E-08, -6.7232123643183028E-01});
    tc->forces.insert(tc->forces.end(), {-1.8318845151918217E+00, -1.1891263627228624E-08,  6.7231740327306944E-01});
    tc->forces.insert(tc->forces.end(), {-1.8318849650582809E+00, -1.1887914508073987E-08, -6.7231724057143949E-01});
    tc->forces.insert(tc->forces.end(), { 1.8319025840150689E+00, -1.7209541413180732E-08,  6.7232337280961452E-01});
    tc->forces.insert(tc->forces.end(), {-1.8724194378490425E-01, -1.4150653174992383E-03, -9.7492172121341703E-02});
    tc->forces.insert(tc->forces.end(), { 1.8724199020372631E-01,  1.4071405182263523E-03,  9.7492843254014180E-02});
    tc->forces.insert(tc->forces.end(), { 1.8724262792453786E-01,  1.4078170361121389E-03, -9.7494264218411716E-02});
    tc->forces.insert(tc->forces.end(), {-1.8722002921893108E-01, -1.5052041010490567E-03,  9.7406114507174929E-02});
    tc->forces.insert(tc->forces.end(), { 1.8724199351778087E-01, -1.4070432897454149E-03,  9.7492840855752036E-02});
    tc->forces.insert(tc->forces.end(), {-1.8724342147002679E-01,  1.4113018645029820E-03, -9.7494391056537663E-02});
    tc->forces.insert(tc->forces.end(), {-1.8722150687274652E-01,  1.5014407228555147E-03,  9.7408333407915326E-02});
    tc->forces.insert(tc->forces.end(), { 1.8724263122792467E-01, -1.4077198328067257E-03, -9.7494261795749951E-02});
    tc->forces.insert(tc->forces.end(), { 3.7388341392462404E-01, -1.3715928523707947E-09, -3.8772275949617063E-01});
    tc->forces.insert(tc->forces.end(), {-3.7387627622364711E-01,  2.3364380021413511E-08,  3.8771897753133538E-01});
    tc->forces.insert(tc->forces.end(), {-3.7387487512154755E-01,  2.3368437040326350E-08, -3.8771970398333933E-01});
    tc->forces.insert(tc->forces.end(), { 3.7388124072323409E-01, -1.3741103948322268E-09,  3.8771477215808314E-01});
    // clang-format on


    structureName = "MAPbI3.tetragonal";
    version = "6.5.1";
    testCases.push_back(TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff, "../../data/ff/" + version + "/refit/ML_FF.MAPbI3");
    // Last modifications: NEC 5.0.1
    // GNU                        : 5.0E-11
    // Intel 2022, PALGO = serial : 5.0E-10
    tc->toleranceEnergy = 6.0E-10;
    // GNU                        : 3.0E-10
    tc->toleranceStress = 2.0E-09;
    tc->toleranceForces = 2.0E-10;

    // clang-format off
    tc->energy = -1.9702722728307296E+02;

    tc->stress.push_back( 3.8449920205712743E+01); // xx
    tc->stress.push_back(-6.9529306446392869E-01); // xy
    tc->stress.push_back(-2.0979852267565686E+00); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 3.8455817011513155E+01); // yy
    tc->stress.push_back(-2.0811420852182261E+00); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back( 2.5943786295592449E+01); // zz

    tc->forces.insert(tc->forces.end(), {-3.1973167025491134E-01, -3.7697357314257923E-01, -1.2572378280232935E-01});
    tc->forces.insert(tc->forces.end(), {-3.7736666116353285E-01, -3.1929378749036852E-01, -1.2565734822444971E-01});
    tc->forces.insert(tc->forces.end(), { 5.6655839356550830E-02,  7.7316580846429325E-02, -1.1942709479417896E-01});
    tc->forces.insert(tc->forces.end(), { 7.7222974764913979E-02,  5.6733485924968179E-02, -1.1945082811573654E-01});
    tc->forces.insert(tc->forces.end(), { 3.1018038801811826E-01, -1.1133417151527734E-02,  3.1658388029406508E-02});
    tc->forces.insert(tc->forces.end(), { 1.4402483361584156E-02, -9.9699335170838925E-04, -3.2028010964489152E-04});
    tc->forces.insert(tc->forces.end(), { 1.2096260944901494E-01, -7.6692781575023364E-02, -1.6750049231020975E-01});
    tc->forces.insert(tc->forces.end(), {-2.9973706949751899E-01, -4.0242405741888471E-01,  4.2839776507913846E-02});
    tc->forces.insert(tc->forces.end(), {-4.0260734109162782E-01, -2.9947374234193913E-01,  4.2958412364263468E-02});
    tc->forces.insert(tc->forces.end(), {-7.6791640400301925E-02,  1.2121512042498138E-01, -1.6727248398770225E-01});
    tc->forces.insert(tc->forces.end(), {-1.0855759042858643E-02,  3.1015222727659503E-01,  3.2028540133411491E-02});
    tc->forces.insert(tc->forces.end(), {-9.8618798712482060E-04,  1.4403638143189687E-02, -3.0256494611848267E-04});
    tc->forces.insert(tc->forces.end(), {-3.9410880124582720E-01, -1.7679978365366070E-01,  1.4934615047311409E-01});
    tc->forces.insert(tc->forces.end(), {-1.7693396737262848E-01, -3.9414343026303988E-01,  1.4909573462353737E-01});
    tc->forces.insert(tc->forces.end(), { 2.6539682047849167E-01, -3.1701354367917367E-01, -1.1350598588662575E-02});
    tc->forces.insert(tc->forces.end(), {-3.1682000327724519E-01,  2.6565577510371174E-01, -1.0679289255305818E-02});
    tc->forces.insert(tc->forces.end(), { 2.4621118966132122E-01,  2.3483078626636245E-01, -1.0571293209935412E-01});
    tc->forces.insert(tc->forces.end(), { 1.8910140853606430E-01, -3.5578111968692883E-01, -2.0157750386387893E-01});
    tc->forces.insert(tc->forces.end(), {-3.5586602328123085E-01,  1.8960981845114697E-01, -2.0094911045519356E-01});
    tc->forces.insert(tc->forces.end(), { 2.3490021565226782E-01,  2.4615054193211566E-01, -1.0569990970929184E-01});
    tc->forces.insert(tc->forces.end(), { 4.8117827444218542E-01, -1.9338337183316359E-01,  9.0444446017040525E-01});
    tc->forces.insert(tc->forces.end(), { 5.7903451997024835E-01,  2.6832545042365347E-01,  8.5532526056481017E-01});
    tc->forces.insert(tc->forces.end(), { 2.6976089968702949E-01,  5.7783949816718472E-01,  8.5568187204479751E-01});
    tc->forces.insert(tc->forces.end(), {-1.9196695771507979E-01,  4.8028503749728219E-01,  9.0522058310178266E-01});
    tc->forces.insert(tc->forces.end(), { 1.3921616422231980E-01, -4.5140331837367521E-03, -2.6772249286623540E-01});
    tc->forces.insert(tc->forces.end(), {-1.0204535087239912E-01, -7.1637500184254757E-02,  1.2033003988199834E-01});
    tc->forces.insert(tc->forces.end(), {-2.6152630542131250E-01,  3.9494731174851000E-01, -1.8123190744614026E-01});
    tc->forces.insert(tc->forces.end(), {-4.1129005925115703E-01,  3.3162191570406996E-01, -1.5625372923188446E-01});
    tc->forces.insert(tc->forces.end(), { 2.5766863994428766E-02,  1.2217899594491858E-02, -8.7655547369335868E-02});
    tc->forces.insert(tc->forces.end(), { 4.0227049564351997E-01,  9.0430753478638926E-01, -1.0305620546952841E-01});
    tc->forces.insert(tc->forces.end(), { 9.6974195913309469E-02,  5.8000162855238423E-02, -1.6192701602139326E-01});
    tc->forces.insert(tc->forces.end(), {-3.5127044185805062E-01, -3.3989136031250972E-01, -1.7565037937025110E-01});
    tc->forces.insert(tc->forces.end(), { 3.8115336239328490E-01, -9.0472649099074665E-01, -1.7355155529944916E-01});
    tc->forces.insert(tc->forces.end(), { 2.6247856239764628E-01, -5.0862246239894521E-02, -2.1861520900534265E-01});
    tc->forces.insert(tc->forces.end(), {-1.2273982205840821E-01, -4.5863688477383996E-01, -1.3552048414677402E-02});
    tc->forces.insert(tc->forces.end(), { 2.0988262891159026E-02,  1.2934429274484521E-01,  1.6535179391757074E-01});
    tc->forces.insert(tc->forces.end(), { 9.0450065844755445E-01,  4.0168684077361805E-01, -1.0363527871272706E-01});
    tc->forces.insert(tc->forces.end(), { 1.2136950787704364E-02,  2.5858366778338093E-02, -8.7639847039548274E-02});
    tc->forces.insert(tc->forces.end(), { 3.3112170497661858E-01, -4.1136686387493898E-01, -1.5710975757388254E-01});
    tc->forces.insert(tc->forces.end(), { 3.9453471701531395E-01, -2.6162358757981963E-01, -1.8198840504430422E-01});
    tc->forces.insert(tc->forces.end(), {-7.1578147199769529E-02, -1.0212840542272460E-01,  1.2029492490186092E-01});
    tc->forces.insert(tc->forces.end(), {-4.7141869778940179E-03,  1.3952819568577893E-01, -2.6755649649088159E-01});
    tc->forces.insert(tc->forces.end(), {-3.4036646699733458E-01, -3.5080357147212365E-01, -1.7566292353568527E-01});
    tc->forces.insert(tc->forces.end(), { 5.7888906545069167E-02,  9.7115776248004659E-02, -1.6188197610869973E-01});
    tc->forces.insert(tc->forces.end(), {-9.0462920414890113E-01,  3.8205538402295453E-01, -1.7206902242576100E-01});
    tc->forces.insert(tc->forces.end(), {-5.0909916784802285E-02,  2.6276986600086050E-01, -2.1825387609406013E-01});
    tc->forces.insert(tc->forces.end(), { 1.2955100617393428E-01,  2.0697134876465531E-02,  1.6522662024518070E-01});
    tc->forces.insert(tc->forces.end(), {-4.5874749087973687E-01, -1.2236809665459956E-01, -1.3164664178208642E-02});
    // clang-format on


    structureName = "CsPbBr3.triclinic";
    version = "6.5.1";
    testCases.push_back(TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff, "../../data/ff/" + version + "/refit/ML_FF.CsPbBr3");
    // Last modifications: Intel oneAPI 2025.3
    tc->toleranceEnergy = 7.0E-10;
    // GNU                : 9.0E-10
    // Intel Classic 2022 : 1.0E-09
    tc->toleranceStress = 2.0E-09;
    // GNU                : 3.0E-11
    tc->toleranceForces = 4.0E-11;

    // clang-format off
    tc->energy = -1.2753979769733697E+02;

    tc->stress.push_back( 2.4557303326638058E+00); // xx
    tc->stress.push_back( 5.8196146343187083E-01); // xy
    tc->stress.push_back( 3.0828852252729879E+00); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 1.2318373685395829E+00); // yy
    tc->stress.push_back(-1.6886600871636328E+00); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back(-1.3627755249290183E+00); // zz

    tc->forces.insert(tc->forces.end(), {-8.8615545751082811E-03,  3.1747493195434504E-02, -2.8887001327838679E-01});
    tc->forces.insert(tc->forces.end(), {-1.8146434830475175E-01, -2.5239321433377433E-02,  3.3358743503910243E-02});
    tc->forces.insert(tc->forces.end(), {-4.4766651076594581E-03, -2.5393056662038260E-03,  9.0055183762551810E-02});
    tc->forces.insert(tc->forces.end(), { 5.4031250428496810E-02, -1.3017725766698968E-01,  7.1119726291029084E-02});
    tc->forces.insert(tc->forces.end(), {-1.0329467074922374E-01, -9.0685042613822903E-02, -1.2187166955482953E-01});
    tc->forces.insert(tc->forces.end(), {-2.1540985651260353E-02,  3.2729578715023316E-03,  1.3755256052970630E-01});
    tc->forces.insert(tc->forces.end(), { 1.0675669994577692E-02,  4.4018551097865122E-02, -3.3523059028234305E-02});
    tc->forces.insert(tc->forces.end(), { 5.3952528974676277E-02, -9.0452083875018074E-02, -8.5718197419609954E-02});
    tc->forces.insert(tc->forces.end(), { 1.0084877342672841E-01,  2.2124519612099356E-02,  3.6286633374670957E-01});
    tc->forces.insert(tc->forces.end(), {-3.8275446278822928E-02,  5.8862742192133371E-03, -2.2468651520665048E-01});
    tc->forces.insert(tc->forces.end(), {-5.2303510052106386E-02,  9.5196817887476595E-03,  1.1508060570373627E-01});
    tc->forces.insert(tc->forces.end(), { 1.9938895139836425E-02, -4.8070530690676405E-02, -1.8082815287520362E-01});
    tc->forces.insert(tc->forces.end(), {-2.1065199176613847E-03, -5.0942735509237148E-02, -1.5295093485136732E-01});
    tc->forces.insert(tc->forces.end(), {-8.5629815394674229E-03,  4.9887303636290932E-03,  1.7043408885843014E-01});
    tc->forces.insert(tc->forces.end(), {-1.0667025156061879E-02, -2.3357108151631140E-02, -1.1096424898348660E-01});
    tc->forces.insert(tc->forces.end(), { 2.8226922598800615E-02, -3.8835111076183752E-02,  1.2823968748764331E-01});
    tc->forces.insert(tc->forces.end(), { 1.9313595473956205E-01,  2.5223085708290419E-02,  6.6280404531378873E-02});
    tc->forces.insert(tc->forces.end(), {-5.1002426208592890E-02, -1.0894572214969090E-03, -5.4678398663315600E-02});
    tc->forces.insert(tc->forces.end(), { 5.5138968481022274E-02,  1.1392912658666862E-01,  4.8841861302638997E-02});
    tc->forces.insert(tc->forces.end(), {-4.8291967868604666E-02,  4.1675199567271719E-02, -9.9110776708807330E-02});
    tc->forces.insert(tc->forces.end(), {-9.4655382126920928E-02, -6.8091236892569851E-02, -6.3261347313065736E-02});
    tc->forces.insert(tc->forces.end(), {-3.0980678744233336E-02,  4.6698946973768224E-02,  5.5390207898560265E-02});
    tc->forces.insert(tc->forces.end(), {-1.4254953669271689E-01, -1.5394076332111161E-02,  7.5380430882907315E-02});
    tc->forces.insert(tc->forces.end(), { 8.7762840399329786E-02, -8.0016001932886105E-02,  1.3658593680752907E-02});
    tc->forces.insert(tc->forces.end(), {-7.4026083132358589E-02,  2.5486700798021322E-01,  2.1878329339071440E-02});
    tc->forces.insert(tc->forces.end(), {-5.2793051539107703E-02,  7.2521925723868480E-02,  1.2221152657549536E-02});
    tc->forces.insert(tc->forces.end(), { 1.6533191953030803E-02,  1.3617666525759198E-01,  1.6041751540299999E-02});
    tc->forces.insert(tc->forces.end(), {-1.7860785316657810E-01, -4.0453888355105973E-02, -3.6567792106050635E-02});
    tc->forces.insert(tc->forces.end(), { 4.4388797804383665E-03, -6.4965371279070816E-02,  6.8406503079958811E-02});
    tc->forces.insert(tc->forces.end(), { 1.3001683130701649E-01,  1.3477362550799146E-01, -6.2268049156493289E-02});
    tc->forces.insert(tc->forces.end(), { 6.5083227792370052E-02, -1.2986643307562809E-01,  7.6090496199498175E-02});
    tc->forces.insert(tc->forces.end(), { 2.8497813802910908E-02, -4.6371776114255955E-02, -1.6084192703236991E-01});
    tc->forces.insert(tc->forces.end(), { 7.5394913230365637E-02,  4.5576829108029293E-02, -6.2564062803558279E-02});
    tc->forces.insert(tc->forces.end(), { 4.3841934561966237E-02, -1.2406176114084635E-02,  1.3737773142221443E-01});
    tc->forces.insert(tc->forces.end(), { 7.0762249062000349E-02,  1.0973840352855671E-01, -6.8517431664357489E-02});
    tc->forces.insert(tc->forces.end(), {-8.5557119547857524E-02, -7.4881839186296734E-02,  1.8482022203979265E-02});
    tc->forces.insert(tc->forces.end(), { 9.1938955733091909E-02,  4.5000306757345794E-02,  1.4169887071283027E-01});
    tc->forces.insert(tc->forces.end(), {-8.6169530970348185E-02,  7.5885754036490199E-02, -8.8635594055472749E-02});
    tc->forces.insert(tc->forces.end(), { 1.4794796714433253E-01, -7.3407586318204682E-02,  8.4406977798832092E-02});
    tc->forces.insert(tc->forces.end(), {-1.9804312211112764E-03, -1.1638274537972623E-01, -4.9004092432929539E-02});
    // clang-format on

    structureName = "CsPbBr3.triclinic";
    version = "6.6.0";
    testCases.push_back(TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff, "../../data/ff/" + version + "/refit/ML_FF.CsPbBr3");
    // TODO: Current tolerances / 100 work for GNU toolchain, actually determine tolerances for
    // other compiler manually!
    tc->toleranceEnergy = 7.0E-10 * 100;
    tc->toleranceStress = 2.0E-09 * 100;
    tc->toleranceForces = 5.0E-11 * 100;

    // clang-format off
    tc->energy = -1.2753731818831979E+02;

    tc->stress.push_back( 2.4639778410476767E+00); // xx
    tc->stress.push_back( 5.7113760006126468E-01); // xy
    tc->stress.push_back( 3.0842324854801753E+00); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 1.2278547094965879E+00); // yy
    tc->stress.push_back(-1.6717369120985357E+00); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back(-1.3851950275378919E+00); // zz

    tc->forces.insert(tc->forces.end(), {-8.8971295520471685E-03,  3.2361114956681675E-02, -2.8416020153817245E-01});
    tc->forces.insert(tc->forces.end(), {-1.8202911111249867E-01, -2.4235257207313919E-02,  3.2679096004149588E-02});
    tc->forces.insert(tc->forces.end(), {-3.5665539243691781E-03, -4.2513529558547069E-03,  8.7225405051414018E-02});
    tc->forces.insert(tc->forces.end(), { 5.1764717827602925E-02, -1.2823209387057580E-01,  7.0154846875924728E-02});
    tc->forces.insert(tc->forces.end(), {-1.0381117172497871E-01, -8.7537521639803675E-02, -1.1703239009088094E-01});
    tc->forces.insert(tc->forces.end(), {-1.6510203694656042E-02,  6.9581798269557080E-03,  1.4204274830186689E-01});
    tc->forces.insert(tc->forces.end(), { 1.1790982103907060E-02,  4.4104736701089553E-02, -3.5737924966434173E-02});
    tc->forces.insert(tc->forces.end(), { 5.6396405797894317E-02, -9.4290619206904056E-02, -9.1063108220760400E-02});
    tc->forces.insert(tc->forces.end(), { 1.0479339686336428E-01,  2.0701603541694251E-02,  3.6318882132306402E-01});
    tc->forces.insert(tc->forces.end(), {-3.7353862348298870E-02,  7.5406338899794067E-03, -2.2814724255959171E-01});
    tc->forces.insert(tc->forces.end(), {-4.9648206838759083E-02,  8.0292232356809819E-03,  1.1448532259051343E-01});
    tc->forces.insert(tc->forces.end(), { 2.0490987119342725E-02, -4.3224771280710512E-02, -1.8431826523643816E-01});
    tc->forces.insert(tc->forces.end(), {-8.1165481087788641E-04, -5.2430322209128744E-02, -1.5039938299499331E-01});
    tc->forces.insert(tc->forces.end(), {-7.6603955490370159E-03,  5.8276949408768821E-03,  1.6737961119268407E-01});
    tc->forces.insert(tc->forces.end(), {-8.8766635098657034E-03, -2.6846177610290137E-02, -1.1365052327259070E-01});
    tc->forces.insert(tc->forces.end(), { 2.7157700761624628E-02, -3.7662743179548513E-02,  1.2729745806556478E-01});
    tc->forces.insert(tc->forces.end(), { 1.9398781756717762E-01,  2.2484959006690126E-02,  6.5107892311750279E-02});
    tc->forces.insert(tc->forces.end(), {-5.1406356225708126E-02, -4.2585825397065146E-04, -5.4036416763742892E-02});
    tc->forces.insert(tc->forces.end(), { 5.1780081699652029E-02,  1.0887664014993557E-01,  4.9031135253155521E-02});
    tc->forces.insert(tc->forces.end(), {-5.0325050573869608E-02,  3.7979190696133550E-02, -9.7794095607553033E-02});
    tc->forces.insert(tc->forces.end(), {-9.6024828848581836E-02, -6.6918971804820190E-02, -6.3347770419485125E-02});
    tc->forces.insert(tc->forces.end(), {-3.0214929012424837E-02,  4.4534846399067553E-02,  5.7267571397070241E-02});
    tc->forces.insert(tc->forces.end(), {-1.4209381521414982E-01, -1.4647969479188555E-02,  7.5021572626737060E-02});
    tc->forces.insert(tc->forces.end(), { 8.9735914181038490E-02, -8.1607155377391868E-02,  1.2100491942247570E-02});
    tc->forces.insert(tc->forces.end(), {-7.6172080168503437E-02,  2.5362143346782146E-01,  2.4048953375501208E-02});
    tc->forces.insert(tc->forces.end(), {-5.0186597323548118E-02,  7.2783380367977588E-02,  1.5025437327196824E-02});
    tc->forces.insert(tc->forces.end(), { 9.1478826952577242E-03,  1.3523783044120688E-01,  1.7385281255852945E-02});
    tc->forces.insert(tc->forces.end(), {-1.7763471139233139E-01, -4.4504042004846728E-02, -3.6549685898062123E-02});
    tc->forces.insert(tc->forces.end(), { 3.0860987593983400E-03, -6.6907588082984415E-02,  6.7977092983520115E-02});
    tc->forces.insert(tc->forces.end(), { 1.2368607952163151E-01,  1.3828043400289666E-01, -5.6584998120152268E-02});
    tc->forces.insert(tc->forces.end(), { 7.5113094311376707E-02, -1.3153181780703496E-01,  7.6887613114504233E-02});
    tc->forces.insert(tc->forces.end(), { 2.9013757662070528E-02, -4.7061252595620752E-02, -1.6332220962838673E-01});
    tc->forces.insert(tc->forces.end(), { 7.1681993847584169E-02,  4.6376945884959719E-02, -6.0599415467090377E-02});
    tc->forces.insert(tc->forces.end(), { 4.3155286522486992E-02, -1.1932612734499793E-02,  1.3866493929614965E-01});
    tc->forces.insert(tc->forces.end(), { 7.4342768945675153E-02,  1.1340444638694412E-01, -6.8022876243610600E-02});
    tc->forces.insert(tc->forces.end(), {-8.3580371829230754E-02, -7.2343360915636964E-02,  1.6619583006105925E-02});
    tc->forces.insert(tc->forces.end(), { 8.8804636877218318E-02,  4.6092624838767783E-02,  1.3961730840983522E-01});
    tc->forces.insert(tc->forces.end(), {-8.4541763904730366E-02,  8.1958107255875670E-02, -9.2498268394879415E-02});
    tc->forces.insert(tc->forces.end(), { 1.4348635568577875E-01, -7.5039940350177692E-02,  8.5459939788559258E-02});
    tc->forces.insert(tc->forces.end(), {-8.0705011916156755E-03, -1.1552259742493266E-01, -4.7403346070543265E-02});
    // clang-format on

    /*============================================================================================+
     | ML_DESC_TYPE = 1
     +============================================================================================*/
    structureName = "MAPbI3.tetragonal";
    version = "6.5.1";
    testCases.push_back(
        TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version + ".DESC_TYPE_1"));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff,
                                "../../data/ff/" + version + "/refit/ML_FF.MAPbI3.DESC_TYPE_1");
    tc->toleranceEnergy = 7.0E-10;
    tc->toleranceStress = 7.0E-09;
    tc->toleranceForces = 6.0E-10;

    // clang-format off
    tc->energy = -1.9900393500057123E+02;

    tc->stress.push_back( 2.1103707788776699E+01); // xx
    tc->stress.push_back(-3.5071143544592553E+00); // xy
    tc->stress.push_back( 1.9114546259798972E+00); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 2.1104794852484623E+01); // yy
    tc->stress.push_back( 1.8736849883107771E+00); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back( 5.6104510824483704E+01); // zz

    tc->forces.insert(tc->forces.end(), {-4.7406964483064223E-01, -9.4677143142453846E-01, -4.5123859325589649E-01});
    tc->forces.insert(tc->forces.end(), {-9.4765896309352249E-01, -4.7281396044584334E-01, -4.5069251609150118E-01});
    tc->forces.insert(tc->forces.end(), { 2.7901389639194135E-01,  4.7038682237699841E-01, -3.9923788226176232E-01});
    tc->forces.insert(tc->forces.end(), { 4.7014348205374795E-01,  2.7910879809008093E-01, -3.9945830376850960E-01});
    tc->forces.insert(tc->forces.end(), { 6.4881690178588014E-01, -1.6561937375681175E+00,  1.9479335807300124E+00});
    tc->forces.insert(tc->forces.end(), {-1.9062319276560817E+00,  2.0794897312499288E+00,  1.2027357710739779E+00});
    tc->forces.insert(tc->forces.end(), { 1.4406725760228156E+00,  2.3192595850033624E+00,  4.1615450818579278E+00});
    tc->forces.insert(tc->forces.end(), {-1.1363172118768803E+00, -9.0377042409673181E-01, -1.5172029876027482E-01});
    tc->forces.insert(tc->forces.end(), {-9.0482797971920559E-01, -1.1354397287127491E+00, -1.5198738287866398E-01});
    tc->forces.insert(tc->forces.end(), { 2.3251722270311146E+00,  1.4340726025337951E+00,  4.1605253153219852E+00});
    tc->forces.insert(tc->forces.end(), {-1.6534432598697621E+00,  6.4785589196004700E-01,  1.9505881178274216E+00});
    tc->forces.insert(tc->forces.end(), { 2.0793911585903375E+00, -1.9092311186285922E+00,  1.1981403908284238E+00});
    tc->forces.insert(tc->forces.end(), { 3.7381848779539928E-01,  2.2280312970166078E+00, -1.0540493236996846E+00});
    tc->forces.insert(tc->forces.end(), { 2.2271054346091463E+00,  3.7330406190759496E-01, -1.0561862285621284E+00});
    tc->forces.insert(tc->forces.end(), {-2.0401298278133448E+00,  7.8058084906389902E-01, -1.0335151972643852E+00});
    tc->forces.insert(tc->forces.end(), { 7.7780304487496044E-01, -2.0395420811320988E+00, -1.0367643101369719E+00});
    tc->forces.insert(tc->forces.end(), {-7.8306308319403539E-01,  2.7732717048319261E-01, -8.4703742632068979E-01});
    tc->forces.insert(tc->forces.end(), {-9.8373062425377134E-01, -3.1035151392320232E-01, -8.5137001766093823E-01});
    tc->forces.insert(tc->forces.end(), {-3.1209693781171871E-01, -9.8250695121106368E-01, -8.5214438667140924E-01});
    tc->forces.insert(tc->forces.end(), { 2.7574199896284346E-01, -7.8230094401173611E-01, -8.4825817915964608E-01});
    tc->forces.insert(tc->forces.end(), {-8.9115557000269674E-01,  2.3708844027438691E+00, -3.5561015268530322E-01});
    tc->forces.insert(tc->forces.end(), {-5.1473562897159719E-01, -2.3286094406339481E+00, -7.1524788735904221E-01});
    tc->forces.insert(tc->forces.end(), {-2.3298314524330164E+00, -5.1210296860223858E-01, -7.1315484365229509E-01});
    tc->forces.insert(tc->forces.end(), { 2.3697793982333062E+00, -8.9258522968798737E-01, -3.5937011032557864E-01});
    tc->forces.insert(tc->forces.end(), { 1.5634642990283018E-01, -6.0308681811737197E-01,  5.9645656991753816E-02});
    tc->forces.insert(tc->forces.end(), { 1.2883605543913974E+00,  1.5576213655558293E+00,  9.3880173126418887E-01});
    tc->forces.insert(tc->forces.end(), { 1.1248905164419565E+00,  1.7127306588840513E-01,  9.1072434187371776E-01});
    tc->forces.insert(tc->forces.end(), { 9.5441407966591074E-01, -1.5288120864518679E+00, -9.1150880279224988E-01});
    tc->forces.insert(tc->forces.end(), {-2.7369975530107410E-01, -7.9359758084997900E-01, -1.7077712783089585E+00});
    tc->forces.insert(tc->forces.end(), {-2.0830659980477870E+00, -1.4962412843026116E+00, -6.1324291614063282E-02});
    tc->forces.insert(tc->forces.end(), {-2.4294664592301446E-01,  8.9530186293331449E-01, -1.6386350449797211E+00});
    tc->forces.insert(tc->forces.end(), { 7.8371077690951874E-01,  1.2339367628407056E+00,  6.9724472453870628E-03});
    tc->forces.insert(tc->forces.end(), {-1.7924537196278736E+00,  1.5524592682425589E+00,  8.7098980412440782E-02});
    tc->forces.insert(tc->forces.end(), {-9.1722361610715647E-02,  6.6949716049677543E-01, -4.8017304206252459E-01});
    tc->forces.insert(tc->forces.end(), { 3.9783640019157557E-01, -2.1262723128474972E-01,  3.6034483908107307E-01});
    tc->forces.insert(tc->forces.end(), { 3.0251832597094730E-01, -3.5881915386762864E-01,  9.8893556045292119E-01});
    tc->forces.insert(tc->forces.end(), {-1.4979300240440607E+00, -2.0818316964431665E+00, -6.1998749221121946E-02});
    tc->forces.insert(tc->forces.end(), {-7.9577720435433530E-01, -2.7111425638319542E-01, -1.7071694217958826E+00});
    tc->forces.insert(tc->forces.end(), {-1.5291187257372110E+00,  9.5665083015923402E-01, -9.0864535479177921E-01});
    tc->forces.insert(tc->forces.end(), { 1.7319679455897055E-01,  1.1237061633943650E+00,  9.1182160136077051E-01});
    tc->forces.insert(tc->forces.end(), { 1.5597033710348744E+00,  1.2860677139651702E+00,  9.3848888591837842E-01});
    tc->forces.insert(tc->forces.end(), {-6.0289607696087644E-01,  1.5674572746649651E-01,  6.0521087711047146E-02});
    tc->forces.insert(tc->forces.end(), { 1.2345530995929761E+00,  7.8274414008575321E-01,  6.4525565685639849E-03});
    tc->forces.insert(tc->forces.end(), { 8.9322410655506979E-01, -2.4175205789369716E-01, -1.6399451197396078E+00});
    tc->forces.insert(tc->forces.end(), { 1.5511645170505177E+00, -1.7937575740424041E+00,  8.3243738423551650E-02});
    tc->forces.insert(tc->forces.end(), { 6.6887181390431605E-01, -9.1688063513984272E-02, -4.8105005889685332E-01});
    tc->forces.insert(tc->forces.end(), {-3.5744417687341673E-01,  3.0165627934725531E-01,  9.8969660315636188E-01});
    tc->forces.insert(tc->forces.end(), {-2.1190259251571206E-01,  3.9758578042426496E-01,  3.6104791661754015E-01});
    // clang-format on

    structureName = "MAPbI3.tetragonal";
    version = "6.6.0";
    testCases.push_back(
        TestCase_Predictor("TestCase_Predictor_" + structureName + "_" + version + ".DESC_TYPE_1"));
    tc = &(testCases.back());
    tc->structure = "../../data/structure/POSCAR." + structureName;
    io::readMlffAndConvertUnits(tc->ff,
                                "../../data/ff/" + version + "/refit/ML_FF.MAPbI3.DESC_TYPE_1");
    // TODO: Current tolerances / 100 work for GNU toolchain, actually determine tolerances for
    // other compiler manually!
    tc->toleranceEnergy = 7.0E-10 * 100;
    tc->toleranceStress = 7.0E-09 * 100;
    tc->toleranceForces = 6.0E-10 * 100;

    // clang-format off
    tc->energy = -2.0170313339543520E+02;

    tc->stress.push_back( 1.4363801678801160E+01); // xx
    tc->stress.push_back(-7.6550932461830268E-01); // xy
    tc->stress.push_back(-1.8214583214998981E-01); // xz
    tc->stress.push_back(          tc->stress[1]); // yx
    tc->stress.push_back( 1.4365430225450632E+01); // yy
    tc->stress.push_back(-1.9928255107763940E-01); // yz
    tc->stress.push_back(          tc->stress[2]); // zx
    tc->stress.push_back(          tc->stress[5]); // zy
    tc->stress.push_back( 3.0127056791169522E+01); // zz

    tc->forces.insert(tc->forces.end(), {-4.1552946751492925E-01, -5.2011451909617412E-01, -1.4593454873766759E-01});
    tc->forces.insert(tc->forces.end(), {-5.2060525340073116E-01, -4.1495704679233830E-01, -1.4581340775172924E-01});
    tc->forces.insert(tc->forces.end(), { 2.5386093052229852E-01,  5.8668966229199516E-02, -1.4527890162366258E-01});
    tc->forces.insert(tc->forces.end(), { 5.8698922344467774E-02,  2.5398264996297076E-01, -1.4505389750781766E-01});
    tc->forces.insert(tc->forces.end(), {-1.0157162980865115E+00,  1.2849240648251241E+00,  4.1173573552774428E-01});
    tc->forces.insert(tc->forces.end(), { 1.0680281865401882E+00, -8.5416649444167758E-01,  7.5878952028378721E-01});
    tc->forces.insert(tc->forces.end(), {-9.2084282104872051E-01, -9.9741512179354630E-01,  2.3072733606141413E-01});
    tc->forces.insert(tc->forces.end(), { 1.0388227043503571E+00,  7.3147556325017649E-01,  5.1054835913506891E-01});
    tc->forces.insert(tc->forces.end(), { 7.3287082003157999E-01,  1.0376651876215566E+00,  5.1090130872477657E-01});
    tc->forces.insert(tc->forces.end(), {-9.9786440334925075E-01, -9.2033377361271129E-01,  2.3081618525163902E-01});
    tc->forces.insert(tc->forces.end(), { 1.2846074071943512E+00, -1.0171874912167542E+00,  4.0908345157687342E-01});
    tc->forces.insert(tc->forces.end(), {-8.5246089101379263E-01,  1.0678152470013917E+00,  7.6100373976119517E-01});
    tc->forces.insert(tc->forces.end(), {-5.4002204901877549E-01, -6.4325932703630428E-01,  8.7890478406779479E-02});
    tc->forces.insert(tc->forces.end(), {-6.4357731104624127E-01, -5.3962350207144005E-01,  8.8009841894953866E-02});
    tc->forces.insert(tc->forces.end(), { 7.2445269348767005E-01, -6.8949332098915439E-01,  1.7688408972796851E-01});
    tc->forces.insert(tc->forces.end(), {-6.8872547813299789E-01,  7.2478334712224801E-01,  1.7851355157834728E-01});
    tc->forces.insert(tc->forces.end(), {-9.3526900240402666E-02,  1.9697852821232179E-01,  3.8033187949327035E-01});
    tc->forces.insert(tc->forces.end(), {-3.9950401431418892E-01, -2.7686419148015384E-01,  4.3117039948242186E-01});
    tc->forces.insert(tc->forces.end(), {-2.7667786761751684E-01, -3.9978573131833672E-01,  4.3102867919285670E-01});
    tc->forces.insert(tc->forces.end(), { 1.9734380111758412E-01, -9.4118317170292065E-02,  3.7999651959777281E-01});
    tc->forces.insert(tc->forces.end(), { 3.3223931941401141E-01, -2.0043210684463800E-01,  1.3380109431041162E-01});
    tc->forces.insert(tc->forces.end(), { 6.1799489641406979E-01,  2.1632164381326519E-01,  2.9026761599465745E-01});
    tc->forces.insert(tc->forces.end(), { 2.1713634042408250E-01,  6.1749185129663087E-01,  2.9072987551309015E-01});
    tc->forces.insert(tc->forces.end(), {-2.0001934259329165E-01,  3.3224015786711164E-01,  1.3441482027612237E-01});
    tc->forces.insert(tc->forces.end(), { 1.0375395515800388E-01, -3.7495016205316356E-01, -2.9922172950659687E-01});
    tc->forces.insert(tc->forces.end(), {-1.6994700579415459E-01, -5.7505810932630942E-02,  2.6518032530958291E-01});
    tc->forces.insert(tc->forces.end(), {-3.6859387183137179E-02,  5.1175560836817435E-01, -3.0221773695057902E-01});
    tc->forces.insert(tc->forces.end(), { 3.4675113125114904E-01,  4.8881499732965117E-01, -3.4792666286320834E-01});
    tc->forces.insert(tc->forces.end(), {-2.0824003051576409E-01, -2.2907174285976487E-01, -6.6946953811192511E-01});
    tc->forces.insert(tc->forces.end(), {-2.2917380639835364E-01, -3.2982904856029771E-01, -3.9810692822957977E-01});
    tc->forces.insert(tc->forces.end(), { 2.8540955130882034E-01,  4.2644589285966927E-01, -5.1970901506196676E-01});
    tc->forces.insert(tc->forces.end(), { 4.6169527759357820E-01, -4.6609743952100730E-01, -6.7440349679599088E-01});
    tc->forces.insert(tc->forces.end(), {-2.1931957157158305E-01,  4.0442431163656845E-01, -2.2236584485789790E-01});
    tc->forces.insert(tc->forces.end(), { 3.8906939702111054E-01,  3.2053036788457501E-01, -3.7751977808264653E-01});
    tc->forces.insert(tc->forces.end(), { 2.4713093725900337E-01, -4.9695361321679976E-01,  4.8461416262441472E-02});
    tc->forces.insert(tc->forces.end(), {-5.5094064633821216E-01,  4.2539461670447298E-01,  3.7513261517754226E-01});
    tc->forces.insert(tc->forces.end(), {-3.3046552786621208E-01, -2.2845859956572834E-01, -3.9799003182426174E-01});
    tc->forces.insert(tc->forces.end(), {-2.3000478221288975E-01, -2.0729004872521772E-01, -6.6944443952466448E-01});
    tc->forces.insert(tc->forces.end(), { 4.8868338269591960E-01,  3.4677240973018769E-01, -3.4809023749677470E-01});
    tc->forces.insert(tc->forces.end(), { 5.1137832005457484E-01, -3.6908204319265572E-02, -3.0284991009928580E-01});
    tc->forces.insert(tc->forces.end(), {-5.7332419574389171E-02, -1.7020790524944779E-01,  2.6505039399596664E-01});
    tc->forces.insert(tc->forces.end(), {-3.7521393610483128E-01,  1.0438993326995351E-01, -2.9866954604971468E-01});
    tc->forces.insert(tc->forces.end(), {-4.6651500262629153E-01,  4.6283403740316936E-01, -6.7333318085511917E-01});
    tc->forces.insert(tc->forces.end(), { 4.2606859371287925E-01,  2.8567761714660778E-01, -5.1987122947816466E-01});
    tc->forces.insert(tc->forces.end(), { 4.0399718723258765E-01, -2.1937679278093722E-01, -2.2308454446272130E-01});
    tc->forces.insert(tc->forces.end(), { 3.2039759950911034E-01,  3.8925554922424599E-01, -3.7744057745432313E-01});
    tc->forces.insert(tc->forces.end(), { 4.2539803521936576E-01, -5.5170285620884985E-01,  3.7400688524380438E-01});
    tc->forces.insert(tc->forces.end(), {-4.9670517629359345E-01,  2.4746061909735981E-01,  4.9319065545809650E-02});
    // clang-format on

    return;
}

} //namespace vaspml

#endif
