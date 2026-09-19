#include "symbol.inc"
      MODULE setexm_struct_def
      USE prec
#ifdef USELIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE
      !
      ! the table has two section one for low densities
      !   (1...NEXCHF(1),0...RHOEXC(1))
      ! and one for high densities
      !   (...NEXCHF(2),...RHOEXC(2))
      !
      INTEGER, PARAMETER :: NEXCH=4000

      TYPE exctable
        REAL(q) :: EXCTAB(NEXCH,5,6) ! table including spline coeff.
        REAL(q) :: RHOEXC(2)         ! maximal densities
        INTEGER :: NEXCHF(2)         ! number of points
      END TYPE
      !
      ! NSMA is the size of the table EXCTAB
      !
      INTEGER, PARAMETER :: NSMA=2000
      !
      ! SXCTAG is the tag (GGA, METAGGA or XC) used in INCAR to specify the functional
      ! SZXCNAM is the list of the functional components specified in INCAR with the GGA, METAGGA or XC tag
      ! SZXC_C is the string of the coefficients (XC_C) multiplying the functional components
      !
      CHARACTER(LEN=7), SAVE :: SXCTAG
      CHARACTER(LEN=400), SAVE :: SZXCNAM
      CHARACTER(LEN=400), SAVE :: SZXC_C
      !
      ! stores the exchange correlation interpolation table
      !
      TYPE (exctable), SAVE :: EXCTAB
      !
      ! parameters that determine the spherical integration mesh
      ! for potentials in the GGA/LDA case
      !
      INTEGER, SAVE :: PHPTS_FACT=3, THPTS_FACT=3
      !
      ! for meta-GGA functionals
      !
      INTEGER, SAVE :: LMAXTAU=6
      LOGICAL, SAVE :: LMIXTAU=.FALSE.
      !
      ! The maximum number of parameter of a functional
      !
      INTEGER, SAVE :: NPARAMMAX=200
      !
      ! xc related variables
      !
      TYPE xc_info
         INTEGER :: NXC                            !< Number of components of the functional
         INTEGER :: IVDW_NL                        !< ID number of the nonlocal vdW functional
         INTEGER :: LFCI                           !< Interpolation of correlation from paramagnetic to 
                                                   !< ferromagnetic case according to 
                                                   !< Vosko, Wilk and Nusair, CAN. J. PHYS. 58, 1200 (1980).
         INTEGER :: ITAUNCL                        !< Definition of the two-component tau in non-collinear calculations
         INTEGER :: ILAPNCL                        !< Definition of the two-component density Laplacian in non-collinear calculations
         INTEGER, ALLOCATABLE :: ID(:)             !< ID number of the components of the functional
         INTEGER, ALLOCATABLE :: FAMILY(:)         !< Family (LDA, GGA or MGGA) of the components of the functional
         INTEGER, ALLOCATABLE :: KIND(:)           !< Kind (exchange, correlation or exchange-correlation) of the components of the functional
         INTEGER, ALLOCATABLE :: ID_LDA_TABLE(:)   !< ID number of the LDA-table functional for each component of the functional
         INTEGER, ALLOCATABLE :: NPARAM(:)         !< Number of parameters for each component of the functional
         REAL(q), ALLOCATABLE :: COEFF(:)          !< Multiplying coefficient for each component of the functional
         REAL(q), ALLOCATABLE :: PARAM(:,:)        !< Parameters for each component of the functional
         REAL(q) :: PARAM1                         !< mu in the exchange GGA optB86b/B86R (GGA=MK) or
                                                   !< beta in the exchange GGA optB88 (GGA=BO)
         REAL(q) :: PARAM2                         !< kappa in the exchange GGA optB86b/B86R (GGA=MK) or
                                                   !< mu in the exchange GGA optB88 (GGA=BO).
         REAL(q) :: Zab_VDW                        !< Z_ab in the kernel of the nonlocal vdW-DF functional of Dion et al.
         REAL(q) :: GAMMA_VDW                      !< gamma in the kernel of the nonlocal vdW-DF3 functional
         REAL(q) :: ALPHA_VDW                      !< alpha in the kernel of the nonlocal vdW-DF3 functional
         REAL(q) :: BPARAM                         !< b in the kernel of the nonlocal rVV10 functional
         REAL(q) :: CPARAM                         !< C in the kernel of the nonlocal rVV10 functional
         REAL(q) :: ALDAX                          !< Fraction of LDA exchange in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: BLDAX                          !< Fraction of LDA to the short-range exchange in a Hartree-Fock/DFT hybrid functional type calculation with LMODELHF
         REAL(q) :: ALDAC                          !< Fraction of LDA correlation in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AGGAX                          !< Fraction of gradient corrections to the exchange in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: BGGAX                          !< Fraction of gradient corrections to the short-range exchange in a Hartree-Fock/DFT hybrid functional type calculation with LMODELHF
         REAL(q) :: AGGAC                          !< Fraction of gradient corrections to the correlation in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AMGGAX                         !< Fraction of meta-GGA exchange in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AMGGAC                         !< Fraction of meta-GGA correlation in a Hartree-Fock/DFT hybrid functional type calculation
         REAL(q) :: AEXX                           !< Fraction of Hartree-Fock exchange (WARINING: to be used only for the PUSH_XC_TYPE and POP_XC_TYPE subroutines)
         REAL(q) :: BEXX                           !< Fraction of short-range Hartree-Fock exchange (WARINING: to be used only for the PUSH_XC_TYPE and POP_XC_TYPE subroutines)
         REAL(q) :: LDASCREEN                      !< Screened LDA exchange parameter
         REAL(q) :: LDASCREENC                     !< Screened LDA correlation parameter
         REAL(q) :: NABMAG_REG                     !< Regularization paramater for the norm of the magnetization when calculating m nabla m / |m|
         LOGICAL :: LSFBXC                         !< .TRUE. if the sources and drains are removed from the XC B field
         LOGICAL :: LNOXC                          !< .TRUE. if no functional is chosen
         LOGICAL :: LDOLDA                         !< .TRUE. if at least one component of the functional is a LDA
         LOGICAL :: LDOGGA                         !< .TRUE. if at least one component of the functional is a GGA
         LOGICAL :: LDOMETAGGA                     !< .TRUE. if at least one component of the functional is a MGGA
         LOGICAL :: LUSE_LIBXC                     !< .TRUE. if at least one component of the functional is from Libxc
         LOGICAL :: LTAU                           !< .TRUE. if at least one component of the functional depends on tau 
         LOGICAL :: LLAP                           !< .TRUE. if the Laplacian has to be computed (e.g., for MBJ and LL-MGGA)
         LOGICAL :: LMU                            !< .TRUE. if the nonlocal component de_xc/dtau of a MGGA has to be computed
         LOGICAL :: LDLAP                          !< .TRUE. if the component LAP(de_xc/dLAP(rho)) of a LL-MGGA has to be computed
         LOGICAL :: LUSE_VDW                       !< .TRUE. is a nonlocal vdW functional is used
         LOGICAL :: LSPIN_VDW                      !< .TRUE. is the spin-polarized formulation of the vdW functional of Dion et al. is used
         LOGICAL :: LUSE_LONGRANGE_HF              !< Short range LDA exchange interaction only
                                                   !< long range contribution is done in HF
                                                   !< the default is that HF treats short range, and LDA long range
                                                   !< but using this flag the behavior can inverted
                                                   !< the flag should be identical to LRHFCALC in fock.F
         LOGICAL :: LRANGE_SEPARATED_CORR          !< .TRUE. for range-separated LDA correlation
                                                   !< .FALSE. for complete LDA correlation (default)
         LOGICAL :: LUSE_THOMAS_FERMI              !< Thomas-Fermi screening in local exchange.
                                                   !< Should be identical to L_THOMAS_FERMI in fock.F.
         LOGICAL :: LUSE_MODEL_HF                  !< No local exchange in the short range
                                                   !< a fraction LDAX of the local exchange in the long range limit
                                                   !< should be identical to L_MODEL_HF
         LOGICAL :: LNO_AUG_XC                     !< Compute XC from soft charge density
         CHARACTER(LEN=40), ALLOCATABLE :: NAME(:) !< Name of the components of the functionals
      END TYPE xc_info
      TYPE (xc_info), SAVE :: XC_DATA, &
     &                        XC_DATA_PP, &
     &                        XC_DATA_NOXC, &
     &                        XC_DATA_STACK_2(5), &
     &                        XC_DATA_TABLE, &
     &                        XC_DATA_WPBE_TABLE
      !
      ! libxc related variables
      !
      CHARACTER(LEN=35) :: LIBXC_VERSION='     F    Libxc'
#ifdef USELIBXC
      INTEGER, PARAMETER :: &
     &   LIBXC_FAMILY_UNKNOWN = -1, &
     &   LIBXC_FAMILY_NONE = 0, &
     &   LIBXC_FAMILY_LDA = 1, &
     &   LIBXC_FAMILY_GGA = 2, &
     &   LIBXC_FAMILY_MGGA = 4, &
     &   LIBXC_FAMILY_LCA = 8, &
     &   LIBXC_FAMILY_OEP = 16, &
     &   LIBXC_FAMILY_HYB_GGA = 32, &
     &   LIBXC_FAMILY_HYB_MGGA = 64, &
     &   LIBXC_FAMILY_HYB_LDA = 128
      INTEGER, ALLOCATABLE :: LIBXC_ID(:)
      INTEGER, ALLOCATABLE :: LIBXC_FAMILY(:)
      INTEGER, ALLOCATABLE :: LIBXC_KIND(:)
      REAL(q), ALLOCATABLE :: LIBXC_AMGGAXC(:)
      CHARACTER(LEN=128), ALLOCATABLE :: LIBXC_PARAM_DESCR(:,:)
      LOGICAL, SAVE :: LTBOUNDLIBXC
      LOGICAL, ALLOCATABLE :: LIBXC_HYB(:)
      LOGICAL, ALLOCATABLE :: LIBXC_NOEXC(:)
      TYPE (xc_f03_func_t), ALLOCATABLE :: LIBXC_FUNC_NSP(:)
      TYPE (xc_f03_func_t), ALLOCATABLE :: LIBXC_FUNC_SPN(:)
      TYPE (xc_f03_func_info_t), ALLOCATABLE :: LIBXC_INFO_NSP(:)
      TYPE (xc_f03_func_info_t), ALLOCATABLE :: LIBXC_INFO_SPN(:)
#endif

GPACC !$ACC DECLARE CREATE(EXCTAB, &
GPACC !$ACC&               LMAXTAU, &
GPACC !$ACC&               LMIXTAU, &
GPACC !$ACC&               EXCTAB%RHOEXC, &
GPACC !$ACC&               EXCTAB%NEXCHF, &
GPACC !$ACC&               EXCTAB%EXCTAB)

GPOMP !$OMP DECLARE TARGET(EXCTAB, &
GPOMP !$OMP&               LMAXTAU, &
GPOMP !$OMP&               LMIXTAU)
! In OpenMP offloading, as long as the derived type does not contain
! pointers, all components of the derived type are automatically
! declared as target data. In fact, it is not possible to declare
! components of a derived type as target data.

#ifdef USELIBXC
GPACC !$ACC DECLARE CREATE(LIBXC_ID, &
GPACC !$ACC&               LIBXC_FAMILY, &
GPACC !$ACC&               LIBXC_HYB, &
GPACC !$ACC&               LIBXC_NOEXC, &
GPACC !$ACC&               LIBXC_KIND, &
GPACC !$ACC&               LIBXC_FUNC_NSP, &
GPACC !$ACC&               LIBXC_FUNC_SPN, &
GPACC !$ACC&               LTBOUNDLIBXC, &
GPACC !$ACC&               LIBXC_AMGGAXC)

GPOMP !$OMP DECLARE TARGET(LIBXC_ID, &
GPOMP !$OMP&               LIBXC_FAMILY, &
GPOMP !$OMP&               LIBXC_HYB, &
GPOMP !$OMP&               LIBXC_NOEXC, &
GPOMP !$OMP&               LIBXC_KIND, &
GPOMP !$OMP&               LIBXC_FUNC_NSP, &
GPOMP !$OMP&               LIBXC_FUNC_SPN, &
GPOMP !$OMP&               LTBOUNDLIBXC, &
GPOMP !$OMP&               LIBXC_AMGGAXC)
#endif

    CONTAINS

      SUBROUTINE ALLOCATE_XC_DATA(XC_DATA_TEMP,NXC)
      IMPLICIT NONE
      INTEGER :: NXC
      TYPE (xc_info) :: XC_DATA_TEMP

      ALLOCATE(XC_DATA_TEMP%ID(NXC), &
     &         XC_DATA_TEMP%FAMILY(NXC), &
     &         XC_DATA_TEMP%KIND(NXC), &
     &         XC_DATA_TEMP%ID_LDA_TABLE(NXC), &
     &         XC_DATA_TEMP%NPARAM(NXC), &
     &         XC_DATA_TEMP%COEFF(NXC), &
     &         XC_DATA_TEMP%PARAM(NXC,NPARAMMAX), &
     &         XC_DATA_TEMP%NAME(NXC))

      END SUBROUTINE ALLOCATE_XC_DATA

      END MODULE setexm_struct_def

      MODULE xc_name
      ENUM, BIND(C)
         ENUMERATOR :: ID_XC_NOXC,ID_XC_COHSM,ID_XC_SLATER,ID_XC_HL,ID_XC_PZ,ID_XC_VWN5, &
        &   ID_XC_WI,ID_XC_PW91,ID_XC_AM05,ID_XC_OPTB88PBE,ID_XC_OPTB86BPBE,ID_XC_OPTPBE, &
        &   ID_XC_PBE,ID_XC_PBE_X,ID_XC_PBE_C,ID_XC_PZ_C, &
        &   ID_XC_PBESOL,ID_XC_REVPBE,ID_XC_RPBE,ID_XC_B3,ID_XC_B5,ID_XC_ACFDT_RA, &
        &   ID_XC_ACFDT_PL,ID_XC_ACFDT_03,ID_XC_ACFDT_05,ID_XC_ACFDT_10,ID_XC_ACFDT_20, &
        &   ID_XC_BEEF,ID_XC_PW86RPBE,ID_XC_CX,ID_XC_M06L,ID_XC_MBJ,ID_XC_MS0, &
        &   ID_XC_MS1,ID_XC_MS2,ID_XC_R2SCAN,ID_XC_RSCAN,ID_XC_RTPSS,ID_XC_SCAN, &
        &   ID_XC_TPSS,ID_XC_LMBJ,ID_XC_RSCANZT,ID_XC_RSCANZT_X,ID_XC_RSCANZT_C, &
        &   ID_XC_RTPSS_X,ID_XC_RTPSS_C,ID_XC_TPSS_X, &
        &   ID_XC_TPSS_C,ID_XC_M06L_X,ID_XC_M06L_C,ID_XC_MS0_X,ID_XC_MS1_X,ID_XC_MS2_X, &
        &   ID_XC_MS0_C,ID_XC_MS1_C,ID_XC_MS2_C,ID_XC_SCAN_X,ID_XC_SCAN_C,ID_XC_RSCAN_X, &
        &   ID_XC_RSCAN_C,ID_XC_R2SCAN_X,ID_XC_R2SCAN_C,ID_XC_LIBXC,ID_XC_PBETEST, &
        &   ID_XC_SCANL,ID_XC_RSCANL,ID_XC_R2SCANL,ID_XC_OFR2,ID_XC_SREGTM1,ID_XC_SREGTM2, &
        &   ID_XC_SREGTM3,ID_XC_SREGTM2L,ID_XC_TASK_X,ID_XC_CC_C, &
        &   ID_XC_LAK,ID_XC_LAK_X,ID_XC_LAK_C,ID_XC_PW92,ID_XC_PW92_C, &
        &   ID_XC_MSPBEL,ID_XC_MSRPBEL,ID_XC_MSB86BL,ID_XC_RMSPBEL,ID_XC_RMSRPBEL,ID_XC_RMSB86BL
      END ENUM
      END MODULE xc_name

      MODULE xc_lda_table_name
      ENUM, BIND(C)
         ENUMERATOR :: ID_XC_LDA_TABLE_NOXC,ID_XC_LDA_TABLE_SLATER,ID_XC_LDA_TABLE_HL, &
        &   ID_XC_LDA_TABLE_PZ,ID_XC_LDA_TABLE_VWN5,ID_XC_LDA_TABLE_WI,ID_XC_LDA_TABLE_PZ_C, &
        &   ID_XC_LDA_TABLE_PW92,ID_XC_LDA_TABLE_B3,ID_XC_LDA_TABLE_B5,ID_XC_LDA_TABLE_ACFDT_RA, &
        &   ID_XC_LDA_TABLE_ACFDT_03,ID_XC_LDA_TABLE_ACFDT_05,ID_XC_LDA_TABLE_ACFDT_10, &
        &   ID_XC_LDA_TABLE_ACFDT_20,ID_XC_LDA_TABLE_ACFDT_PL,ID_XC_LDA_TABLE_X,ID_XC_LDA_TABLE_PW92_C, &
        &   ID_XC_LDA_TABLE_PW92_RELA
      END ENUM
      END MODULE xc_lda_table_name

      MODULE xc_family_name
      ENUM, BIND(C)
      !XC_FAM_XXX is here used because XC_FAMILY_XXX already exist (imported from Libxc).
         ENUMERATOR :: XC_FAM_UNKNOWN,XC_FAM_NONE,XC_FAM_LDA,XC_FAM_GGA,XC_FAM_MGGA
      END ENUM
      END MODULE xc_family_name
