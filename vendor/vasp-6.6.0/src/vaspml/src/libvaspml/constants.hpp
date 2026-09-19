#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "types.hpp"

namespace vaspml::constants
{

// Pi variants
constexpr Real PI = 3.14159265358979323846264338327950288419716939937510;
constexpr Real PI2 = 2 * PI;
constexpr Real PI3 = 3 * PI;
constexpr Real PI4 = 4 * PI;
constexpr Real SQRT_PI = 1.77245385090551602729816748334114518279754945612239;

// Euler's number
constexpr Real E = 2.71828182845904523536028747135266249775724709369995;

// Other constants
constexpr Real SQRT2 = 1.41421356237309504880168872420969807856967187537695;
constexpr Real EPS_TOL = 1.0e-10;

/// conversion factor from angstroem to bohr
constexpr Real AUTOA = 0.529177249;
constexpr Real RYTOEV = 13.605826;
constexpr Real EVTOJ = 1.60217733e-19;
constexpr Real AMTOKG = 1.660540e-27;

constexpr Real EUNIT = 2.0 * RYTOEV;
constexpr Real MUNIT =
    1.054571726e-034 * 1.054571726e-034 / 2.0 / (RYTOEV * EVTOJ * AUTOA * AUTOA * 1.0e-20) / AMTOKG;
constexpr Real FUNIT = 2.0 * RYTOEV / AUTOA;
constexpr Real SUNIT =
    2.0 * RYTOEV / ((AUTOA * 1.0e-10) * (AUTOA * 1.0e-10) * (AUTOA * 1.0e-10)) * EVTOJ / 1.0e8;

// key list which is needed in Descriptor collector and Kernel routines
// for extracting needed descriptors from descriptor map
const Vec1String descriptorKeyList = {"SHS2-2-body", "SHS3-3-body"};
const Vec1String bodyKeyList = {"2-body", "3-body"};
const MapString  bodyDescMap = {
    {"2-body", "SHS2-2-body"},
    {"3-body", "SHS3-3-body"}
};

constexpr Size mlffHeaderSize = 4096;
const Int      mlffVersionMajor = 0;
const Int      mlffVersionMinor = 2;
const Int      mlffVersionPatch = 5;

/**************************************************************************************************
 * INCAR tag information
 *
 * The first entry is current INCAR tag name. The next-to-last entry is type of INCAR tag, the last
 * is a short description. Entries in between the first and the next-to-ast are alternative
 * (deprecated) names of INCAR tag.
 **************************************************************************************************/
const Vec2String tagInfo{
    // clang-format off
    /* General, mode of operation: */
    {"ML_LMLFF", "ML_FF_LMLFF", "bool",
     "Main control tag to enable machine learning"},
    {"ML_MODE", "ML_FF_MODE", "String",
     "Machine learning operation mode in strings (supertag)"},
    {"ML_ISTART", "ML_FF_ISTART", "Int",
     "Machine learning operation mode"},
    {"ML_TYPE", "String",
     "Machine learning force field type"},
    {"ML_LFAST", "ML_FF_LFAST", "bool",
     "Precontraction of weights for fast execution, but no Bayesian error estimation"},
    {"ML_LFORCESLOW", "ML_FF_LFORCESLOW", "bool",
     "Forcing slow algorithm althoug force field was trained with ML_LFAST=.TRUE. No error estimation"},
    {"ML_LIB", "bool",
     "Enable experimental VASPml C++ library, requires build with -Dlibvaspml"},
    {"ML_PALGO", "String",
     "Inside VASPml, use C++ Parallel Algorithms code path"},
    {"ML_NCSHMEM", "ML_FF_NCSHMEM", "Int",
     "Number of MPI tasks sharing memory segments when using MPI shared memory (-Duse_shmem)"},

    /* Descriptors: */
    {"ML_RCUT1", "ML_FF_RCUT1_MB", "Real",
     "Cutoff radius of radial descriptors"},
    {"ML_SION1", "ML_FF_SION1_MB", "Real",
     "Gaussian width for broadening the atomic distribution for radial descriptors"},
    {"ML_MRB1", "ML_FF_MRB1_MB", "Int",
     "Number of radial basis functions for atomic distribution for radial descriptors"},
    {"ML_RCUT2", "ML_FF_RCUT2_MB", "Real",
     "Cutoff radius of angular descriptors"},
    {"ML_SION2", "ML_FF_SION2_MB", "Real",
     "Gaussian width for broadening the atomic distribution for angular descriptors"},
    {"ML_MRB2", "ML_FF_MRB2_MB", "Int",
     "Number of radial basis functions for atomic distribution for angular descriptors"},
    {"ML_MRB2_MAX", "ML_FF_MRB2_MAX_MB", "Int",
     "Distance for summation of summed up descriptor"},
    {"ML_MRB2_RED", "ML_FF_MRB2_MB_RED", "Int",
     "Number of radial basis functions for reduced (color)atomic distribution for angular descriptors"},
    {"ML_MRB2_MAX_RED", "ML_FF_MRB2_MAX_MB_RED", "Int",
     "Distance for summation of summed up reduced descriptor"},
    {"ML_LMAX2", "ML_FF_LMAX2_MB", "Int",
     "Maximum angular momentum quantum number of spherical harmonics used to expand atomic distributions"},
    {"ML_LMAX2_RED", "ML_FF_LMAX2_RED_MB", "Int",
     "Maximum angular momentum for reduced (type color) descriptor"},
    {"ML_LAFILT2", "ML_FF_LAFILT2_MB", "bool",
     "Angular filtering enabled"},
    {"ML_AFILT2", "ML_FF_AFILT2_MB", "Real",
     "Angular filtering parameter a_FILT"},
    {"ML_IAFILT2", "ML_FF_IAFILT2_MB", "Int",
     "Angular filtering function type"},
    {"ML_LSPARSDES", "ML_FF_LSPARSDES_MB", "bool",
     "Enable sparsification of angular descriptors"},
    {"ML_NRANK_SPARSDES", "ML_FF_NRANK_SPARSDES_MB", "Int",
     "Number of highest eigenvalues relevant in the sparsification algorithm of the angular descriptors"},
    {"ML_RDES_SPARSDES", "ML_FF_RDES_SPARSDES_MB", "Real",
     "Desired ratio of selected to all descriptors resulting from the angular descriptor sparsification"},
    {"ML_DESC_TYPE", "ML_FF_DESC_TYPE", "Int",
     "Descriptor type (standard, linear-scaling with element types, ...)"},
    {"ML_DESC_RATIO_DUAL", "ML_FF_DESC_RATIO_DUAL", "Real",
     "Ratio between high-resolution type-independent and low-resolution type-dependent descriptor"},
    {"ML_DESC_FACTORE_TESTE", "ML_FF_DESC_FACTORE_TESTE", "Real",
     "Factori di high-resolution type-independent and low-resolution type-dependent descriptor"},
    {"ML_BASIS_TYPE2", "ML_FF_BASIS_TYPE2", "Int",
     "Basis set type for angular descriptor"},
    {"ML_BASIS_TYPE2_RED", "ML_FF_BASIS_TYPE2_RED", "Int",
     "Basis set type for reduced angular descriptor"},

    /* Kernel: */
    {"ML_W1", "ML_FF_W1_MB", "Real",
     "Weight of radial descriptors in the kernel (the angular counterpart is chosen so that the sum is 1.0)"},
    {"ML_NHYP", "ML_FF_NHYP1_MB", "Int",
     "Power of the polynomial kernel"},

    /* Controlling Bayesian criterion and Spilling Factor: */
    {"ML_ICRITERIA", "ML_FF_ICRITERIA", "Int",
     "Enable automatic updating of the Bayesian error estimation threshold during on-the-fly training"},
    {"ML_IUPDATE_CRITERIA", "ML_FF_IUPDATE_CRITERIA", "Int",
     "Decides whether update of threshold is done in the same MD step or the next MD step"},
    {"ML_CTIFOR", "ML_FF_CTIFOR", "Real",
     "Bayesian error estimation threshold (initial or static value depending on other settings)"},
    {"ML_SCLC_CTIFOR", "ML_FF_SCLC_CTIFOR", "Real",
     "Scaling factor for ML_CTIFOR. The interval 0<ML_SCLC_CTIFOR<1 increases num. of local configurations"},
    {"ML_CSIG", "ML_FF_CSIG", "Real",
     "Standard error parameter required for the automatic update of the Bayesian error estimation threshold"},
    {"ML_CSLOPE", "ML_FF_CSLOPE", "Real",
     "Slope parameter required for the automatic update of the Bayesian error estimation threshold"},
    {"ML_CX", "ML_FF_XMIX", "Real",
     "Additional parameter controlling the update of the Bayesian error estimation threshold"},
    {"ML_CDOUB", "ML_FF_CDOUB", "Real",
     "Factor controlling the occurence of enforced first principle calculations"},
    {"ML_CALGO", "ML_FF_CALGO", "Int",
     "Decides what to use for learning criterion. 0: Bayesian error estimate; 1: Spilling factor"},
    {"ML_ESTBLOCK", "ML_IERR", "Int",
     "Time steps between error estimations (prediction modes only)"},
    {"ML_SF_LIMIT", "Real",
     "Maximum allowed spilling factor, higher values trigger graceful exit (fast prediction mode only)"},

    /* Sparsification and regression: */
    {"ML_EPS_LOW", "ML_FF_EPS_LOW", "Real",
     "Threshold for the CUR algorithm used in the sparsification of local reference configurations"},
    {"ML_EPS_REG", "ML_FF_EPS_REG", "Real",
     "Convergence criterion for the optimization of parameters within the Bayesian linear regression"},
    {"ML_LBASIS_DISCARD", "ML_FF_LBASIS_DISCARD", "bool",
     "If ML_MB is exceeded do not stop but instead discard older local reference configurations"},
    {"ML_IALGO_LINREG", "ML_FF_IALGO_LINREG", "Int",
     "Linear regression algorithm"},
    {"ML_ISVD", "ML_FF_ISVD", "Int",
     "Leverage scoring calculation mode (for sparsification of local configurations)"},
    {"ML_IREG", "ML_FF_IREG_MB", "Int",
     "Time dependency of regularization parameters"},
    {"ML_SIGV0", "ML_FF_SIGV0_MB", "Real",
     "Initial regularization parameter (noise, reversed and squared)"},
    {"ML_SIGW0", "ML_FF_SIGW0_MB", "Real",
     "Initial regularization parameter (precision, reversed and squared)"},
    {"ML_SAVECMAT", "ML_FF_SAVECMAT_MB", "bool",
     "Decides whether CMAT is saved so it does not have to be calculated at each evidence approx. step"},

    /* Weighting and sampling: */
    {"ML_WTOTEN", "ML_FF_WTOTEN", "Real",
     "Scaling weight for total energies in the training data"},
    {"ML_WTIFOR", "ML_FF_WTIFOR", "Real",
     "Scaling weight for forces in the training data"},
    {"ML_WTSIF", "ML_FF_WTSIF", "Real",
     "Scaling weight for stresses in the training data"},
    {"ML_MHIS", "ML_FF_MHIS", "Int",
     "Number of estimated errors stored to determine the threshold for the Bayesian error"},
    {"ML_NMDINT", "ML_FF_NMDINT", "Int",
     "Minimum number of MD steps between potential training samples"},
    {"ML_IWEIGHT", "ML_FF_IWEIGHT", "Int",
     "Weighting method for energies, forces and stresses"},
    {"ML_LUSE_NAMES", "ML_FF_LUSE_NAMES", "bool",
     "Use also structure names for distinction into groups for weighting"},

    /* Sizes of static arrays: */
    {"ML_MB", "ML_FF_MB_MB", "Int",
     "Maximum number of local reference configurations of the ML force field"},
    {"ML_MB_MIN", "Int",
     "Minimum number of local reference configurations required (otherwise training is blocked)"},
    {"ML_MCONF", "ML_FF_MCONF", "Int",
     "Maximum number of structures stored for training"},
    {"ML_MCONF_NEW", "ML_FF_MCONF_NEW", "Int",
     "Maximum number of configurations stored temporarily as training candidates"},

    /* Special features: */
    {"ML_LHEAT", "ML_FF_LHEAT_MB", "bool",
     "Enable heat flux calculation (output written to ML_HEAT)"},
    {"ML_LCOUPLE", "ML_FF_LCOUPLE_MB", "bool",
     "Enable thermodynamic integration (t.i.)"},
    {"ML_NATOM_COUPLED", "ML_FF_NATOM_COUPLED_MB", "Int",
     "Number of atoms whose interaction is controlled by the t.i. coupling parameter"},
    {"ML_ICOUPLE", "ML_FF_ICOUPLE_MB", "Vec1Int",
     "List of atoms whose interaction is controlled by the t.i. coupling parameter"},
    {"ML_RCOUPLE", "ML_FF_RCOUPLE_MB", "Real",
     "Thermodynamic integration coupling parameter"},
    {"ML_LEMPPOT", "bool",
     "Model functions enabled for thermodynamic integration"},
    {"ML_EMPPOT_RCUT", "Real",
     "Cutoff radius of the model potential"},
    {"ML_SRPOT_B0", "Real",
     "Height of the soft-repulsive potential used in the model potential"},
    {"ML_SRPOT_N0", "Int",
     "Power parameter used in the model potential"},
    {"ML_SRPOT_S0", "Real",
     "Radius of the soft-repulsive potential used in the model potential"},
    {"ML_MOPOT_NM", "Int",
     "Number of pairs which are coupled by Morse model potential"},
    {"ML_MOPOT_DM", "Vec1Real",
     "Binding strength of Morse model potential"},
    {"ML_MOPOT_RM", "Vec1Real",
     "Equilibrium length of Morse model potential"},
    {"ML_MOPOT_RKM", "Vec1Real",
     "Force constant of Morse model potential"},
    {"ML_MOPOT_IJM", "Vec1Int",
     "List of pairs which are coupled by Morse model potential"},

    /* Reference energies: */
    {"ML_ISCALE_TOTEN", "ML_FF_ISCALE_TOTEN_MB", "Int",
     "Scaling mode of energies"},
    {"ML_EATOM_REF", "ML_FF_EATOM", "Vec1Real",
     "List of reference energies of isolated atoms for each species in the system"},

    /* Output options: */
    {"ML_OUTPUT_MODE", "ML_FF_OUTPUT_MODE", "Int",
     "Controls the verbosity of the output at each MD step when machine learning is used"},
    {"ML_OUTBLOCK", "ML_FF_OUTBLOCK", "Int",
     "Sets the output frequency of file output (only relevant for ML_MODE = run)"},
    {"ML_LEATOM", "ML_FF_LEATOM_MB", "bool",
     "Enable output of kin. and pot. energy for each atom at each MD time step (written to ML_EATOM)"},
    {"ML_GRACE_MODEL", "String",
     "GRACE model name or path"},

    /* Not important, description not used and not yet reviewed: */
    {"ML_IBROAD1", "ML_FF_IBROAD1_MB", "Int",
     "Specifies whether broadening of the atomic distribution is used for the radial descriptor"},
    {"ML_IBROAD2", "ML_FF_IBROAD2_MB", "Int",
     "Specifies whether broadening of the atomic distribution is used for the angular descriptor"},
    {"ML_TOTNUM_INSTANCES", "ML_FF_TOTNUM_INSTANCES", "Int",
     "Specifies the number of machine learning instances"},
    {"ML_FILESTREAM_START", "ML_FF_FILESTREAM_START", "Int",
     "Specifies the starting value for the filehandle numbers within Fortran"},
    {"ML_IMAT_SPARS", "ML_FF_IMAT_SPARS", "Int",
     "Specifies the type of matrix for sparsification"},
    {"ML_LCONF_DISCARD", "ML_FF_LCONF_DISCARD", "bool",
     "Discard training structures if the maximum number in storage is exceeded"},
    {"ML_LMLMB", "ML_FF_LMLMB", "bool",
     "Turns many body interaction on. Code won't work without this tag"},
    {"ML_NDIM_SCALAPACK", "ML_FF_NDIM_SCALAPACK", "Int",
     "Specifies the dimension of scalapack grid. Should never be touched"},
    {"ML_FF_NWRITE", "Int",
     "Deprecated, in past it specified contents of the ML_LOGFILE"},
    {"ML_LNMDINT_RANDOM", "ML_FF_LNMDINT_RANDOM", "bool",
     "Decides wether ML_NMDINT times a random number between 0 and 1 is used instead of ML_NMDINT"},
    {"ML_LTEST", "ML_FF_LTEST", "bool",
     "Whether ab initio calculations are executed at every MD step regardless of the estimated error"},
    {"ML_LTRJ", "ML_FF_LTRJ", "bool",
     "DESCRIPTION OF ML_LTRJ"},
    {"ML_LTOTEN_SYSTEM", "ML_FF_LTOTEN_SYSTEM", "bool",
     "DESCRIPTION OF ML_LTOTEN_SYSTEM"},
    {"ML_NTEST", "ML_FF_NTEST", "Int",
     "This parameter determines how often the test ab initio calculation will be executed"},
    {"ML_INVERSE_SOAP", "ML_FF_INVERSE_SOAP_MB", "Int",
     "Decides how the inverse matrix is calculated for ridge regression"},
    {"ML_LSUPERVEC", "ML_FF_LSUPERVEC_MB", "bool",
     "Specifies whether super-vector is used for kernel or not"},
    {"ML_LSIC", "ML_FF_LSIC", "bool",
     "Specifies whether self-interactions are substracted or not"},
    {"ML_LMETRIC1", "ML_FF_LMETRIC1_MB", "bool",
     "Specifies whether metric function is used for radial descriptor or not"},
    {"ML_LMETRIC2", "ML_FF_LMETRIC2_MB", "bool",
     "Specifies whether metric function is used for angular descriptor or not"},
    {"ML_LVARTRAN1", "ML_FF_LVARTRAN1_MB", "bool",
     "Specifies whether variable transform is done for radial coordinate of radial descriptor or not"},
    {"ML_LVARTRAN2", "ML_FF_LVARTRAN2_MB", "bool",
     "Specifies whether variable transform is done for radial coordinate of angular descriptor or not"},
    {"ML_NMETRIC1", "ML_FF_NMETRIC1_MB", "Int",
     "Specifies the polynomial parameter for the metric function for the radial descriptor"},
    {"ML_NMETRIC2", "ML_FF_NMETRIC2_MB", "Int",
     "Specifies the polynomial parameter for the metric function for the angular descriptor"},
    {"ML_NVARTRAN1", "ML_FF_NVARTRAN1_MB", "Int",
     "Polynomial parameter for the variable transform of the radial coordinate for the radial descriptor"},
    {"ML_NVARTRAN2", "ML_FF_NVARTRAN2_MB", "Int",
     "Polynomial parameter for the variable transform of the radial coordinate for the angular descriptor"},
    {"ML_LWINDOW1", "ML_FF_LWINDOW1_MB", "bool",
     "Specifies wheter a window function is used for the radial descriptor or not"},
    {"ML_LWINDOW2", "ML_FF_LWINDOW2_MB", "bool",
     "Specifies wheter a window function is used for the angular descriptor or not"},
    {"ML_IWINDOW1", "ML_FF_IWINDOW1_MB", "Int",
     "Specifies the type of window function for radial descriptor if ML_LWINDOW1=.TRUE."},
    {"ML_IWINDOW2", "ML_FF_IWINDOW2_MB", "Int",
     "Specifies the type of window function for angular descriptor if ML_LWINDOW2=.TRUE."},
    {"ML_MSPL1", "ML_FF_MSPL1_MB", "Int",
     "Maximum number of radial grid points used in spline-interpolation for radial descriptor"},
    {"ML_MSPL2", "ML_FF_MSPL2_MB", "Int",
     "Maximum number of radial grid points used in spline-interpolation for angular descriptor"},
    {"ML_LNORM1", "ML_FF_LNORM1_MB", "bool",
     "Specifies whether normalization is taken into account for radial descriptor"},
    {"ML_LNORM2", "ML_FF_LNORM2_MB", "bool",
     "Specifies whether normalization is taken into account for angular descriptor"},
    {"ML_NHYP2", "ML_FF_NHYP2_MB", "Int",
     "Sets the hyper parameter (power) of SOAP"},
    {"ML_NR1", "ML_FF_NR1_MB", "Int",
     "Number of radial mesh points for calculating radial descriptor"},
    {"ML_NR2", "ML_FF_NR2_MB", "Int",
     "Number of radial mesh points for calculating angular descriptor"},
    {"ML_RMETRIC1", "ML_FF_RMETRIC1_MB", "Real",
     "Parameter used for metric function for radial descriptor"},
    {"ML_RMETRIC2", "ML_FF_RMETRIC2_MB", "Real",
     "Parameter used for metric function for angular descriptor"},
    {"ML_RANDOM_SEED", "ML_FF_RANDOM_SEED", "Vec1Int",
     "Sets the random seed for ML_LNMDINT_RANDOM"},
    {"ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS", "ML_FF_LDISCARD_STRUCTURES_NOT_GIVING_BASIS", "bool",
     "Specifies whether structures not giving new local configurations are discarded for ML_ISTART=3"},
    {"ML_ICUT1", "ML_FF_ICUT1_MB", "Int",
     "Specifies the type of cut off function for the radial descriptor"},
    {"ML_ICUT2", "ML_FF_ICUT2_MB", "Int",
     "Specifies the type of cut off function for the angular descriptor"},
    {"ML_NBLOCK", "ML_FF_NBLOCK_FFN", "Int",
     "Minimum blocksize for printing of ML_FFN file during learning"},

    /* Temporary, to be renamed. */
    {"ML_RALGO", "String",
     "Can be used to switch between different refitting algorithms."},
    {"ML_SALGO", "String",
     "Can be used to switch between different selecting algorithms."},
    {"ML_SMETRIC", "String",
     "Can be used to choose the norm function for farthest points sampling algorithms."},
    {"ML_SNCONF", "Vec1Int",
     "Number of LRCs per type selected by farthest point sampling or random selection"},
    {"ML_STRESH", "Vec1Real",
     "Number of LRCs per type selected by farthest point sampling with random selection"},
    {"ML_SPAR", "String",
     "Determines if the selection algorithms are executed with openmp or mpi parallelization."}
    // clang-format on
};

/**************************************************************************************************
 * Data from INCAR reader which should be copied to ML_FF header record.
 *
 * @TODO Large overlap with headerSort variable in io::detail::processMlffHeader()!
 **************************************************************************************************/
const Vec1String incarToMlffHeader{"ML_LFAST",
                                   "ML_DESC_TYPE",
                                   "ML_IALGO_LINREG",
                                   "ML_RCUT1",
                                   "ML_RCUT2",
                                   "ML_W1",
                                   "ML_SION1",
                                   "ML_SION2",
                                   "ML_LMAX2",
                                   "ML_MRB1",
                                   "ML_MRB2",
                                   "ML_IWEIGHT",
                                   "ML_WTOTEN",
                                   "ML_WTIFOR",
                                   "ML_WTSIF"};

const Vec1String incarToMlff{"ML_DESC_TYPE",
                             "ML_IWEIGHT",
                             "ML_WTOTEN",
                             "ML_WTIFOR",
                             "ML_WTSIF",
                             "ML_LSIC",
                             "ML_LSUPERVEC",
                             "ML_W1",
                             "ML_W2",
                             "ML_ICUT1",
                             "ML_ICUT2",
                             "ML_RCUT1",
                             "ML_RCUT2",
                             "ML_IBROAD1",
                             "ML_IBROAD2",
                             "ML_SION1",
                             "ML_SION2",
                             "ML_MRB1",
                             "ML_MRB2",
                             "ML_BASIS_TYPE2",
                             "ML_MRB2_MAX",
                             "ML_MRB2_RED",
                             "ML_BASIS_TYPE2_RED",
                             "ML_NR1",
                             "ML_NR2",
                             "ML_MSPL1",
                             "ML_MSPL2",
                             "ML_LMAX2",
                             "ML_LMAX2_RED",
                             "ML_DESC_RATIO_DUAL",
                             "ML_DESC_FACTORE_TESTE",
                             "ML_NHYP",
                             "ML_LNORM1",
                             "ML_LNORM2",
                             "ML_LAFILT2",
                             "ML_IAFILT2",
                             "ML_AFILT2",
                             "ML_EATOM_REF",
                             "ML_LMLMB"};
/**************************************************************************************************
 * List of case-sensitive (string-type) INCAR tags.
 *
 * Normally, all string-type INCAR tags are automatically converted to lowercase. However, for some
 * tags this is not desirable (e.g. file paths). This list contains all INCAR tags which will not be
 * automatically converted to lowercase.
 **************************************************************************************************/
const Vec1String caseSensitiveTags{"ML_GRACE_MODEL"};

} // namespace vaspml::constants

#endif
