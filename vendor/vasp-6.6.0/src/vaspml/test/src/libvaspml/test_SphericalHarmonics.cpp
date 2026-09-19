#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE SphericalHarmonics

#include "test_SphericalHarmonics.hpp"

#include "boost_helpers.hpp"

#include "boost_unit_test_case.hpp"

using namespace vaspml;
namespace bdata = boost::unit_test::data;

const Real tol = 10.0 * std::numeric_limits<Real>::epsilon();

TestCaseContainer<TestCaseIntegrate3SphericalHarmonics> containerIntegrate3SphericalHarmonics;
TestCaseContainer<TestCaseSphericalHarmonics>           containerSphericalHarmonics;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_DATA_TEST_CASE(Integrate3Ylms_CorrectData,
                     bdata::make(containerIntegrate3SphericalHarmonics.testCases),
                     testCase)
{
    math::SphericalHarmonics::Integrate3SphericalHarmonicsData data;
    math::SphericalHarmonics::integrate3SphericalHarmonics(3, data);

    REQUIRE_CLOSE_COLLECTIONS(data.ylm3, testCase.data.ylm3, testCase.tolerance, "ylm3");
    REQUIRE_CLOSE_COLLECTIONS(data.ylm3i, testCase.data.ylm3i, testCase.tolerance, "ylm3i");
    REQUIRE_EQUAL_COLLECTIONS(data.jl, testCase.data.jl, "jl");
    REQUIRE_EQUAL_COLLECTIONS(data.js, testCase.data.js, "js");
    REQUIRE_EQUAL_COLLECTIONS(data.indcg, testCase.data.indcg, "indcg");
    REQUIRE_EQUAL_COLLECTIONS_2D(data.lookup, testCase.data.lookup, "lookup");
}

BOOST_DATA_TEST_CASE(ComputeSphericalHarmonicsAndDerivatives_CorrectResults,
                     bdata::make(containerSphericalHarmonics.testCases),
                     testCase)
{
    math::SphericalHarmonics calculator(testCase.jmax);

    Vec1Real ylm(testCase.ylm.size(), 9999);
    Vec1Real ylmd(testCase.ylmd.size(), 9999);
    calculator.computeSphericalHarmonicsAndGradients(testCase.rhat, testCase.rnorm, ylm, ylmd);

    REQUIRE_CLOSE_COLLECTIONS(ylm, testCase.ylm, testCase.tolerance, "ylm");
    REQUIRE_CLOSE_COLLECTIONS(ylmd, testCase.ylmd, testCase.tolerance, "ylmd");
}

BOOST_DATA_TEST_CASE(ComputeSphericalHarmonicAndDerivative_CorrectResults,
                     bdata::make(containerSphericalHarmonics.testCases),
                     testCase)
{
    math::SphericalHarmonics calculator(testCase.jmax);

    if (testCase.jmax > calculator.get_lhc())
    {
        BOOST_REQUIRE(true);
        return;
    }

    Vec1Real ylm(testCase.ylm.size(), 9999);
    Vec1Real ylmd(testCase.ylmd.size(), 9999);
    Int      iylm = 0;
    for (Size i = 0; i < testCase.rnorm.size(); ++i)
    {
        for (Int l = 0; l <= testCase.jmax; ++l)
        {
            for (Int m = -l; m <= l; ++m)
            {
                calculator.computeSphericalHarmonicAndGradient(testCase.rhat[3 * i],
                                                               testCase.rhat[3 * i + 1],
                                                               testCase.rhat[3 * i + 2],
                                                               testCase.rnorm[i],
                                                               l,
                                                               m,
                                                               ylm[iylm],
                                                               ylmd[3 * iylm],
                                                               ylmd[3 * iylm + 1],
                                                               ylmd[3 * iylm + 2]);
                iylm++;
            }
        }
    }

    REQUIRE_CLOSE_COLLECTIONS(ylm, testCase.ylm, testCase.tolerance, "ylm");
    REQUIRE_CLOSE_COLLECTIONS(ylmd, testCase.ylmd, testCase.tolerance, "ylmd");
}

BOOST_AUTO_TEST_CASE(ComputeYlms_CorrectResults)
{
    math::SphericalHarmonics calculator(3);

    // Set up test vectors and their norms.
    Vec1Real rn; // Must be normalized vector!
    Vec1Real norm;

    rn.insert(rn.end(), {-2.1229201078614458E-002, -0.71305220599719088, -0.70078946377931761});
    rn.insert(rn.end(), {-0.77081456423851702, -0.60817539586965652, -0.18965124680995826});
    rn.insert(rn.end(), {-2.0766665167792766E-002, -0.69751642345205456, 0.71626781627577185});
    rn.insert(rn.end(), {-2.1124923847718552E-002, 0.71642209068291818, -0.69734720589814136});
    rn.insert(rn.end(), {-0.67758773667324901, 0.13071474614968065, -0.72373234987048107});
    rn.insert(rn.end(), {-0.54699846732837210, 0.82624448650790994, -0.13458352520160935});
    rn.insert(rn.end(), {-2.0669025445139558E-002, 0.70096093735194798, 0.71290010218390831});
    rn.insert(rn.end(), {-0.77393646145972272, 0.14930156291254454, 0.61541156711185940});
    rn.insert(rn.end(), {0.84380404793814245, -0.51231971116995889, -0.15975995189638365});
    rn.insert(rn.end(), {0.62247915380642793, 0.13910291081008214, -0.77017535878564414});
    rn.insert(rn.end(), {0.72578469769038045, 0.16218818486175768, 0.66852940495518320});

    norm.push_back(5.3105783689561186);
    norm.push_back(4.2934472421041727);
    norm.push_back(5.4288608752241272);
    norm.push_back(5.3367925418810662);
    norm.push_back(6.0276700416518230);
    norm.push_back(6.0502028116595667);
    norm.push_back(5.4545066160737230);
    norm.push_back(5.2772746915588522);
    norm.push_back(5.0967568086521489);
    norm.push_back(5.6641903090315910);
    norm.push_back(4.8579701415373142);

    // Set up result vectors for spherical harmonics and their derivatives.
    Vec1Real ylm(norm.size() * calculator.get_ldim());
    Vec1Real ylmd(norm.size() * 3 * calculator.get_ldim());

    // Run spherical harmonics calculator.
    calculator.computeSphericalHarmonicsAndGradients(rn, norm, ylm, ylmd);

    // Set up reference data.
    Vec1Real ylm_ref;
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(-0.34839909896814575);
    ylm_ref.push_back(-0.34240749231767487);
    ylm_ref.push_back(-1.0372640972703200E-002);
    ylm_ref.push_back(1.6538483181251055E-002);
    ylm_ref.push_back(0.54594587508936154);
    ylm_ref.push_back(0.14928038428664819);
    ylm_ref.push_back(1.6254061992703309E-002);
    ylm_ref.push_back(-0.27750335153414135);
    ylm_ref.push_back(0.21334952465016629);
    ylm_ref.push_back(-3.0664243832327701E-002);
    ylm_ref.push_back(-0.47435340312052698);
    ylm_ref.push_back(0.14238889499378954);
    ylm_ref.push_back(-1.4122589752159609E-002);
    ylm_ref.push_back(0.51452302744291700);
    ylm_ref.push_back(1.9100877557411737E-002);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(-0.29715602609946684);
    ylm_ref.push_back(-9.2664075576866220E-002);
    ylm_ref.push_back(-0.37662193229829399);
    ylm_ref.push_back(0.51217627342622130);
    ylm_ref.push_back(0.12601587119428639);
    ylm_ref.push_back(-0.28135993660211422);
    ylm_ref.push_back(0.15971522278184827);
    ylm_ref.push_back(0.12251716844226712);
    ylm_ref.push_back(-0.50690718266738943);
    ylm_ref.push_back(-0.25699470658817597);
    ylm_ref.push_back(0.22797552475605307);
    ylm_ref.push_back(0.19959233295674358);
    ylm_ref.push_back(0.28894107845419303);
    ylm_ref.push_back(-6.1475443884204374E-002);
    ylm_ref.push_back(0.23444618462571604);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(-0.34080827659221458);
    ylm_ref.push_back(0.34997025422756123);
    ylm_ref.push_back(-1.0146644764830416E-002);
    ylm_ref.push_back(1.5825662362725937E-002);
    ylm_ref.push_back(-0.54584655408455218);
    ylm_ref.push_back(0.17003350764868555);
    ylm_ref.push_back(-1.6251104978385569E-002);
    ylm_ref.push_back(-0.26554275248668585);
    ylm_ref.push_back(0.19970595237909050);
    ylm_ref.push_back(2.9990682805234564E-002);
    ylm_ref.push_back(-0.49898032623587069);
    ylm_ref.push_back(-0.11622007190839864);
    ylm_ref.push_back(-1.4855789787677442E-002);
    ylm_ref.push_back(-0.50322117826892387);
    ylm_ref.push_back(1.7879388061136391E-002);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(0.35004563309041525);
    ylm_ref.push_back(-0.34072559647031453);
    ylm_ref.push_back(-1.0321690855753180E-002);
    ylm_ref.push_back(-1.6535023569653862E-002);
    ylm_ref.push_back(-0.54583167110461883);
    ylm_ref.push_back(0.14472668488640011);
    ylm_ref.push_back(1.6094775183113936E-002);
    ylm_ref.push_back(-0.28013725642624954);
    ylm_ref.push_back(-0.21639972804928415);
    ylm_ref.push_back(3.0507238931624485E-002);
    ylm_ref.push_back(0.46871582311121129);
    ylm_ref.push_back(0.14795167408782064);
    ylm_ref.push_back(-1.3820883244968796E-002);
    ylm_ref.push_back(0.51685527870246961);
    ylm_ref.push_back(1.9187242066293637E-002);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(6.3867553311486483E-002);
    ylm_ref.push_back(-0.35361744409211987);
    ylm_ref.push_back(-0.33107107017316367);
    ylm_ref.push_back(-9.6767789107145508E-002);
    ylm_ref.push_back(-0.10335780240966466);
    ylm_ref.push_back(0.18020387285837774);
    ylm_ref.push_back(0.53577719014264913);
    ylm_ref.push_back(0.24147439643114851);
    ylm_ref.push_back(0.10491555720057578);
    ylm_ref.push_back(0.18529249282266747);
    ylm_ref.push_back(9.6719879992103117E-002);
    ylm_ref.push_back(0.10291590808835496);
    ylm_ref.push_back(-0.50136810501940077);
    ylm_ref.push_back(-0.46237899284889844);
    ylm_ref.push_back(-0.16306743911781985);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(0.40370513155370297);
    ylm_ref.push_back(-6.5757848474256253E-002);
    ylm_ref.push_back(-0.26726482514368988);
    ylm_ref.push_back(-0.49378214444845780);
    ylm_ref.push_back(-0.12149017895070215);
    ylm_ref.push_back(-0.29825376294437511);
    ylm_ref.push_back(8.0430118163151612E-002);
    ylm_ref.push_back(-0.20948120907119366);
    ylm_ref.push_back(0.10478858928335301);
    ylm_ref.push_back(0.17582324908050642);
    ylm_ref.push_back(-0.34343188734160096);
    ylm_ref.push_back(0.14612174709245029);
    ylm_ref.push_back(0.22736214168461805);
    ylm_ref.push_back(7.4590924792046043E-002);
    ylm_ref.push_back(0.56444083986836002);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(0.34249127473598701);
    ylm_ref.push_back(0.34832478066290584);
    ylm_ref.push_back(-1.0098937751080554E-002);
    ylm_ref.push_back(-1.5829037720426518E-002);
    ylm_ref.push_back(0.54596297432197627);
    ylm_ref.push_back(0.16547954145699928);
    ylm_ref.push_back(-1.6098646881800658E-002);
    ylm_ref.push_back(-0.26817640629227879);
    ylm_ref.push_back(-0.20268950133292463);
    ylm_ref.push_back(-2.9856040485819708E-002);
    ylm_ref.push_back(0.49373463774469006);
    ylm_ref.push_back(-0.12207586350221035);
    ylm_ref.push_back(-1.4558605546899243E-002);
    ylm_ref.push_back(-0.50582263969664532);
    ylm_ref.push_back(1.7971614666195217E-002);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(7.2949118670101079E-002);
    ylm_ref.push_back(0.30069163754496686);
    ylm_ref.push_back(-0.37814729912247785);
    ylm_ref.push_back(-0.12624388734655262);
    ylm_ref.push_back(0.10038543526389761);
    ylm_ref.push_back(4.2954499016359927E-002);
    ylm_ref.push_back(-0.52036929175178237);
    ylm_ref.push_back(0.31502906888134802);
    ylm_ref.push_back(0.15633607418583598);
    ylm_ref.push_back(-0.20555357473594685);
    ylm_ref.push_back(6.0981054490881169E-002);
    ylm_ref.push_back(-0.25407945288432698);
    ylm_ref.push_back(-0.31610828854082729);
    ylm_ref.push_back(0.51293850827436682);
    ylm_ref.push_back(-0.24298897074415826);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(-0.25032069777502031);
    ylm_ref.push_back(-7.8059113798062699E-002);
    ylm_ref.push_back(0.41228477737642821);
    ylm_ref.push_back(-0.47230589631142467);
    ylm_ref.push_back(8.9423092315649863E-002);
    ylm_ref.push_back(-0.29124211930086491);
    ylm_ref.push_back(-0.14728218655256786);
    ylm_ref.push_back(0.24556871557538074);
    ylm_ref.push_back(-0.56635551630106784);
    ylm_ref.push_back(0.19963666604519584);
    ylm_ref.push_back(0.20427178029071127);
    ylm_ref.push_back(0.17124759478022880);
    ylm_ref.push_back(-0.33644099832741331);
    ylm_ref.push_back(-0.10379823763653513);
    ylm_ref.push_back(-3.7544814744463348E-002);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(6.7966031634813961E-002);
    ylm_ref.push_back(-0.37630961490839826);
    ylm_ref.push_back(0.30414487815702468);
    ylm_ref.push_back(9.4602307007958339E-002);
    ylm_ref.push_back(-0.11704868395393900);
    ylm_ref.push_back(0.24585095782798419);
    ylm_ref.push_back(-0.52378749889195753);
    ylm_ref.push_back(0.20110029677434982);
    ylm_ref.push_back(9.3821251754928128E-002);
    ylm_ref.push_back(-0.19277040818602337);
    ylm_ref.push_back(0.12498169453461047);
    ylm_ref.push_back(9.8149733421487653E-003);
    ylm_ref.push_back(0.55928735784268702);
    ylm_ref.push_back(-0.40978055949799064);
    ylm_ref.push_back(0.12099681819360295);
    ylm_ref.push_back(0.28209479177387814);
    ylm_ref.push_back(7.9245554524429926E-002);
    ylm_ref.push_back(0.32664514654206683);
    ylm_ref.push_back(0.35462022639222124);
    ylm_ref.push_back(0.12860792116465375);
    ylm_ref.push_back(0.11846237221910590);
    ylm_ref.push_back(0.10748377255994163);
    ylm_ref.push_back(0.53011368911991485);
    ylm_ref.push_back(0.27338753072447108);
    ylm_ref.push_back(0.14871343643213311);
    ylm_ref.push_back(0.22747687454378837);
    ylm_ref.push_back(9.1522009892582523E-002);
    ylm_ref.push_back(-0.19093706918663142);
    ylm_ref.push_back(0.40955680180108123);
    ylm_ref.push_back(0.48355762588547724);
    ylm_ref.push_back(0.19178866508759387);

    Vec1Real ylmd_ref;
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.3927361604978473E-003, 4.5225914981358313E-002, -4.5975108695190350E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.3687845278262797E-003, -4.5975108695190350E-002, 4.6821067630053892E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {9.1964052705994365E-002, -1.3927361604978473E-003, -1.3687845278262797E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.14656442612209611, 7.3753456610808649E-005, 4.3648710008926551E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.3648710008926551E-003, 2.4346486294402302E-003, -2.6094774461198375E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.7150809448643321E-003, 0.12478315379763374, -0.12707920979122408});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.14404387887362036, 4.3648710008926551E-003, -7.7690133973171597E-005});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.5861525413199954E-003, 7.2175813630206725E-002, -7.3239263751530281E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.2649946150971578E-002, -8.3385566269730371E-002, 8.4461477784321753E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.27162466729152263, -4.2540402820101650E-003, -3.8999186112419495E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-8.2942979179552751E-003, -0.15332351177923570, 0.15625769975647033});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-4.5649466900585080E-003, -0.15332867664458702, 0.15614998077326786});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.12502069253867532, -8.2942979179552768E-003, 4.6521504320609901E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.4268303561388503E-002, -6.4737143798517277E-002, 6.5437711153396638E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.16909592150150859, -2.3972985615599300E-003, 7.5617152094513669E-003});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.3349250579460682E-002, 7.1709167653340983E-002, -1.3126051787522817E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.6636240067396665E-002, -1.3126051787522817E-002, 0.10970873237152644});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.6185927087008749E-002, -5.3349250579460682E-002, -1.6636240067396665E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {2.9143199087092671E-002, -5.1046679832223357E-002, 4.5247962005540418E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.5247962005540418E-002, -1.2559527187000905E-002, -0.14362886168835526});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.2219586513718739E-002, 9.6412966362228684E-003, -8.0582832484579980E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {9.0879112842021520E-003, 4.5247962005540418E-002, -0.18203830537418075});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.15215657023120172, 0.18947127640022884, 1.0823719235601529E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.11353266155095115, -0.12294652753582899, -6.7173700144976270E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.0762160840254667E-002, -1.0790120283918557E-002, 0.28156231558516231});
    ylmd_ref.insert(ylmd_ref.end(),
                    {2.2979937620566413E-002, -6.9176589187661730E-002, 0.12843695636391969});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.1263276105949561E-002, 2.4666834546259026E-002, -0.20616764229596030});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.8182573906876624E-002, 2.2979937620566410E-002, 0.16278375814629958});
    ylmd_ref.insert(ylmd_ref.end(),
                    {6.5310749566286921E-002, -0.10377912864965752, 6.7352236789061798E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.21873905695122592, -0.28692273063491308, 3.1068050019159334E-002});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.3036715305604211E-003, 4.6212851553522685E-002, 4.4965234080222558E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.3387182422373333E-003, 4.4965234080222558E-002, 4.3826888847293410E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.9962113812384101E-002, -1.3036715305604211E-003, 1.3387182422373333E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.14025284471129815, -1.1259611025879705E-004, -4.1759819904021410E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-4.1759819904021410E-003, 3.8835782907163116E-003, 3.6608351795406023E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.7137293383207410E-003, 0.12473775566776636, 0.12157987975684928});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.14402327376128332, -4.1759819904021410E-003, 1.0899146722847305E-004});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.2107803220434032E-003, 7.2138450547845290E-002, 7.0069847735283505E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.1737763865476512E-002, -8.1520290690620320E-002, -7.9045834676287785E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.26567374811310451, 3.6399162868488580E-003, -4.1580206637722135E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-8.1650961761396308E-003, -0.14248035047201732, -0.13898704371419196});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.8010706413285954E-003, 0.16125967243281369, 0.15717720331789473});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.13152803707419167, -8.1650961761396326E-003, -4.1379633545967649E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.3694780756948613E-002, 7.2051786863802841E-002, 6.9768470137158034E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.15829158189148945, -2.5544044973102322E-003, -7.0768604338392986E-003});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.3856051709394415E-003, 4.4562737963584668E-002, 4.5739710183002096E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.3487131496902551E-003, 4.5739710183002096E-002, 4.7031708138667624E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {9.1512732251996229E-002, 1.3856051709394415E-003, -1.3487131496902551E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.14653506248418227, 1.1469621678560771E-004, -4.3211919501348305E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-4.3211919501348305E-003, 3.7861952487544281E-003, 4.0206639433124549E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.6426235117260879E-003, -0.12353445487678737, -0.12702390119195836});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.14263353645620705, -4.3211919501348305E-003, -1.1855611478811225E-004});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.5424613267282715E-003, -7.1453553791050811E-002, -7.3209865852480738E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.2609424649266903E-002, -8.2942669998951696E-002, -8.4829461438515663E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.27023749917765810, -4.3069703322134324E-003, 3.7615951979528489E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.1582619382684681E-003, -0.15408449261365695, -0.15854638273035207});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-4.4236439861773495E-003, 0.15002166615414148, 0.15425928998693783});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.12235093028922898, 8.1582619382684681E-003, 4.6750097534224268E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.4116785908283021E-002, 6.2448603498193020E-002, 6.3729147654473253E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.16986466407197756, 2.3124578620822003E-003, 7.5214669891056074E-003});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {7.1795354749260428E-003, 7.9674912124585842E-002, 7.6684712532694276E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-3.9751154581264513E-002, 7.6684712532694276E-003, 3.8601669720987809E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.3843278241865104E-002, 7.1795354749260428E-003, -3.9751154581264513E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.9368771469624494E-003, -0.11861954277073303, -2.3237496053481205E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.3237496053481205E-002, -0.12669768914578697, -1.1272315210958930E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.11142261898699594, -2.1494750521497689E-002, -0.10820060908647604});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.0723967190177576E-002, -2.3237496053481205E-002, 5.8432447568806797E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.8526849670251910E-002, -3.4165891331193285E-002, 5.7986860977804679E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.6639097537143158E-002, 0.12298689840102453, 3.7791111105644822E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.7120499389571994E-002, 0.22311636345004668, 2.4268582689802608E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.6049361164340545E-002, 0.11387223737242531, -2.2546597243970054E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.14745545458752876, 2.8445913749010652E-002, 0.14319148252435074});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.11595096714218524, 4.6049361164340545E-002, 0.11687509057877099});
    ylmd_ref.insert(ylmd_ref.end(),
                    {7.9238942113411134E-002, 7.5448439632084419E-002, -6.0559835161375745E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {7.4819786116391451E-002, 6.2629412144473262E-002, -5.8737711297030419E-002});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.6498956330341906E-002, 2.5626144710873593E-002, 8.9802047035136574E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.9451630713786112E-003, 8.9802047035136574E-003, 7.9295290386139194E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {5.6594642004176646E-002, 3.6498956330341906E-002, -5.9451630713786112E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {5.9917985543577343E-002, 3.6089241725336979E-002, -2.1967839343638523E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.1967839343638523E-002, 8.8793985053949240E-003, 0.14379864460441322});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.0988553236153404E-003, -4.6808396705905267E-003, -4.1331857478171953E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-9.7597912592624642E-003, -2.1967839343638523E-002, -9.5198987087895454E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.13665560367903948, -9.1988111192382635E-002, -9.3195948823286950E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.23603827360803861, -0.15512532558642955, 6.9928983468847518E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.4390998056593103E-003, -3.6861719651611703E-002, -0.20419745527422192});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.4865673206262675E-002, -3.1141100314011919E-002, -9.0120074408507264E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.2388414943822038E-002, -1.8712775547432028E-002, -0.16523398073361195});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.2239033263512515E-002, -2.4865673206262675E-002, 5.9662174310315834E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {5.5403341941996455E-002, 2.2568132578678943E-002, -8.6628344708341340E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.0899145372456688E-002, 3.3211882376322832E-002, 3.7667053665133572E-002});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.2978187342179805E-003, 4.5564067370774598E-002, -4.4763363937788299E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.3199239200623874E-003, -4.4763363937788299E-002, 4.4051966032467124E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.9539496617776729E-002, 1.2978187342179805E-003, 1.3199239200623874E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.14028389445683462, -7.1653974169579897E-005, 4.1376877516697769E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.1376877516697769E-003, 2.4714336746529614E-003, -2.3100803138050392E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.6443762341848212E-003, -0.12359389599465007, 0.12162969064933077});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.14267328942870955, 4.1376877516697769E-003, 6.8116647079834129E-005});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.1724741782552290E-003, -7.1476932597477838E-002, 7.0100927877028291E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.1707773472865859E-002, -8.1172989223926670E-002, 7.9474117302946748E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.26448422760591545, 3.7016617085556805E-003, 4.0284900513240162E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.0407940769584400E-003, -0.14355712476701940, 0.14138605368247525});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.6608873729228834E-003, -0.15806744205175674, 0.15555536884538707});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.12889800451886266, 8.0407940769584400E-003, -4.1690082648409765E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.3558989718299981E-002, -6.9813192455600304E-002, 6.8250894579279728E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.15911182607526214, 2.4749708542814866E-003, -7.0466314373064621E-003});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.0698321779696164E-002, 9.0522309789273100E-002, -8.5069915939764672E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.4097803421966494E-002, -8.5069915939764672E-003, 5.7520864037925244E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {3.7129113180886790E-002, 1.0698321779696164E-002, 4.4097803421966494E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.1187466163976087E-003, -0.15308398637274290, 2.9443966096555637E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {2.9443966096555637E-002, 0.12172789452985747, 7.4967619643532520E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.10510617739943827, -2.0276233694819887E-002, 0.13709975714408348});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.5221085235144935E-002, 2.9443966096555637E-002, -3.8861062897888127E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.7826315970320017E-002, -4.8734975526318926E-002, -7.3474489882501201E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-8.7345055655464899E-003, 0.18016632615591957, -5.4693587524284153E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-4.0108073622151533E-002, -0.24343991947850094, 8.6200197751335889E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.6844178704363346E-002, 6.8359620376338498E-002, 4.2326503329131744E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {9.0295781877149692E-002, -1.7419131970138172E-002, 0.11778118158979414});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.16543036042862186, 4.6844178704363353E-002, -0.21940844806627960});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-3.5211791053132317E-002, -9.3863308014596442E-002, -2.1510402290169137E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.6528901384728754E-002, 9.8140210412443862E-002, 8.5008777456146489E-002});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.1442357560929149E-002, 7.0703448849515863E-002, -7.8464058884108734E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.2923237006218474E-002, -7.8464058884108734E-003, 9.3418581563317990E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {2.7608723180697781E-002, 4.1442357560929149E-002, 1.2923237006218474E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.6565525265448750E-002, 8.5927887866716540E-002, -2.9609247648229853E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.9609247648229853E-002, -1.6268984802441763E-002, -0.10421561367120531});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-7.9962223094018105E-003, 4.8549450716830781E-003, -5.7802526228516317E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.4520827354952163E-002, -2.9609247648229853E-002, 0.17164585854661354});
    ylmd_ref.insert(ylmd_ref.end(),
                    {9.9567833988396706E-002, 0.15919015046876531, 1.5394906078696962E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.8985954031913756E-002, -1.4662598541725205E-002, -5.3257944279533745E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.2733812330930302E-002, -1.6253206502940425E-002, -0.22640321192845794});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.3924351524754962E-002, -6.3704226176330356E-002, 7.7926075629732389E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.5832030663258832E-002, 1.5684042427471817E-002, -0.18673275606572753});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-3.8825975361379168E-002, -2.3924351524754962E-002, -0.12834629748315193});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-2.4901309533310535E-002, -7.7720975903505601E-002, 0.11771512147024465});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.17477290325325964, 0.28895653285324352, -3.5305732779039424E-003});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-7.4692825543238549E-003, 8.4592538902215533E-002, 9.2415261394225394E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.1355406135259509E-002, 9.2415261394225394E-003, 3.5093827781422442E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {5.2836972132799649E-002, -7.4692825543238549E-003, 4.1355406135259509E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {6.0380631655176315E-003, 0.11542156501640595, 2.5726665866303115E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {2.5726665866303119E-002, -0.14280774657998363, -4.9997185613139209E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.12335806241901053, -2.7566329650474109E-002, -0.10468054850733846});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-3.3431129785789264E-002, 2.5726665866303119E-002, -2.2373511533242532E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {7.5867355083330770E-002, -3.6708501805982029E-002, 5.4688308397097180E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {2.3187990378968796E-002, 0.10813319783098228, 3.8271392179489182E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.8812292535326767E-003, -0.23045953896376892, -3.4445678274022146E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-5.5179195967109469E-002, 0.14629458979337479, -1.8174946494485165E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {0.18627836215351920, 4.1626875756464225E-002, 0.15807415212995396});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-8.8299111686277951E-002, -5.5179195967109469E-002, -8.1332052999312399E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.10956056563905194, 8.4864153332158382E-002, -7.3222697654898597E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {7.5153803924575574E-002, -6.3034523902567099E-002, 4.9356799178651899E-002});
    ylmd_ref.insert(ylmd_ref.end(), {0.0000000000000000, 0.0000000000000000, 0.0000000000000000});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.1839350419642377E-002, 9.7931812134540441E-002, -1.0905374439949865E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-4.8801051062047379E-002, -1.0905374439949865E-002, 5.5626242768165364E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {4.7596953328077854E-002, -1.1839350419642377E-002, -4.8801051062047379E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-1.9524359472205685E-003, 0.15464021802047609, -3.5396749878549530E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-3.5396749878549530E-002, 0.14244104796415979, 3.8714293964151876E-003});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.12635584010313025, -2.8236230962466049E-002, 0.14402764866305173});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-8.0478170658425247E-003, -3.5396749878549530E-002, 1.7324469205335893E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {8.1538986076530526E-002, -5.4730451984876423E-002, -7.5244432515001139E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {1.9130517931736594E-002, 0.16746051081857183, -6.1395584330631352E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-3.7438626001697145E-002, 0.26592728366433283, -2.3870126683837770E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-6.3169785041682386E-002, 0.10204233594333853, 4.3823954304454342E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.13805604931730611, -3.0850829618222667E-002, 0.15736421957744806});
    ylmd_ref.insert(ylmd_ref.end(),
                    {-0.16652326083403271, -6.3169785041682386E-002, 0.19611018801132876});
    ylmd_ref.insert(ylmd_ref.end(),
                    {7.1979201571242851E-002, -0.11294922806836390, -5.0741781167739190E-002});
    ylmd_ref.insert(ylmd_ref.end(),
                    {9.6395202075149414E-002, -0.10499341205422515, -7.9178972953248980E-002});

    // Finally, compare computation results to reference data.
    REQUIRE_CLOSE_COLLECTIONS(ylm, ylm_ref, tol, "ylm");
    REQUIRE_CLOSE_COLLECTIONS(ylmd, ylmd_ref, tol, "ylmd");
}

BOOST_AUTO_TEST_SUITE_END()
