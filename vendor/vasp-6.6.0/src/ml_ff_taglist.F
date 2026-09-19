#if defined(MPI) || defined(MPI_CHAIN)
!==============================================================================
!> Contains list of MLFF related tags, descriptions and taglist type.
!==============================================================================
module ml_ff_taglist

   implicit none
   private
   public :: taglist,                       &
             taglist_get_tag_width,         &
             taglist_description,           &
             taglist_get_all_tags,          &
             taglist_get_alt_tags,          &
             taglist_get_description_width, &
             taglist_state_symbol

   !===========================================================================
   !> Describes state of a ML tag (default value, user-defined,...)
   !===========================================================================
   enum, bind(c)
      !> Dummy state for usage as kind function argument.
      enumerator :: TL_STATE     = 0
      !> Tag is not set yet.
      enumerator :: TL_UNSET     = 1
      !> Tag is set to default value.
      enumerator :: TL_DEFAULT   = 2
      !> Tag is set to user-defined value from INCAR file.
      enumerator :: TL_INCAR     = 3
      !> Tag is set to user-defined value from INCAR file (alternative tag
      !> name).
      enumerator :: TL_INCAR_ALT = 4
      !> Tag was set by user in INCAR but overridden by program, e.g., taken
      !> from ML_FF file.
      enumerator :: TL_OVERRIDE  = 5
   end enum
   ! Note: There is no way to make enums public without listing them.
   public :: TL_STATE,     &
             TL_UNSET,     &
             TL_DEFAULT,   &
             TL_INCAR,     &
             TL_INCAR_ALT, &
             TL_OVERRIDE

   !===========================================================================
   !> Stores state of all ML tags (see #TL_STATE enum).
   !===========================================================================
   type taglist

      private
      !> State for each tag.
      integer(kind(TL_STATE)), allocatable :: state(:)
      !> Actual tag name version used.
      integer                , allocatable :: alt_index(:)

      contains

      procedure, public, pass(this) :: init => taglist_init
      procedure, public, pass(this) :: get_state => taglist_get_state
      procedure, public, pass(this) :: set_state => taglist_set_state
      procedure, public, pass(this) :: get_deprecated => taglist_get_deprecated
      procedure, public, pass(this) :: symbol => taglist_state_symbol_from_tag

   end type taglist

   !> Number of tags listed here.
   integer, parameter                  :: ntags = 128
   !> Maximum number of alternative tags listed here.
   integer, parameter                  :: nalt = 2
   !> Static list of all known MLFF tags.
   !>
   !> The first entry is considered to be the "preferred" tag, all others are
   !> deprecated.
   !> @warning If this list needs to be extended with more alternatives, make
   !> sure to fill it left to right, i.e. empty strings should only appear at
   !> the end of each line (e.g. "A", "B", " " and NOT "A", " ", "B").
   character(len=*), parameter         :: tags(nalt, ntags) = reshape([                            &
      ! General, mode of operation:
      !========================================================================
      "ML_LMLFF                                  ", "ML_FF_LMLFF                               ",  &
      "ML_MODE                                   ", "ML_FF_MODE                                ",  &
      "ML_ISTART                                 ", "ML_FF_ISTART                              ",  &
      "ML_TYPE                                   ", "                                          ",  &
      "ML_LFAST                                  ", "ML_FF_LFAST                               ",  &
      "ML_LFORCESLOW                             ", "ML_FF_LFORCESLOW                          ",  &
      "ML_LIB                                    ", "                                          ",  &
      "ML_PALGO                                  ", "                                          ",  &
      "ML_NCSHMEM                                ", "ML_FF_NCSHMEM                             ",  &

      ! Descriptors:
      !========================================================================
      "ML_RCUT1                                  ", "ML_FF_RCUT1_MB                            ",  &
      "ML_SION1                                  ", "ML_FF_SION1_MB                            ",  &
      "ML_MRB1                                   ", "ML_FF_MRB1_MB                             ",  &
      "ML_RCUT2                                  ", "ML_FF_RCUT2_MB                            ",  &
      "ML_SION2                                  ", "ML_FF_SION2_MB                            ",  &
      "ML_MRB2                                   ", "ML_FF_MRB2_MB                             ",  &
      "ML_MRB2_MAX                               ", "ML_FF_MRB2_MAX_MB                         ",  &
      "ML_MRB2_RED                               ", "ML_FF_MRB2_MB_RED                         ",  &
      "ML_MRB2_MAX_RED                           ", "ML_FF_MRB2_MAX_MB_RED                     ",  &
      "ML_LMAX2                                  ", "ML_FF_LMAX2_MB                            ",  &
      "ML_LMAX2_RED                              ", "ML_FF_LMAX2_RED_MB                        ",  &
      "ML_LAFILT2                                ", "ML_FF_LAFILT2_MB                          ",  &
      "ML_AFILT2                                 ", "ML_FF_AFILT2_MB                           ",  &
      "ML_IAFILT2                                ", "ML_FF_IAFILT2_MB                          ",  &
      "ML_LSPARSDES                              ", "ML_FF_LSPARSDES_MB                        ",  &
      "ML_NRANK_SPARSDES                         ", "ML_FF_NRANK_SPARSDES_MB                   ",  &
      "ML_RDES_SPARSDES                          ", "ML_FF_RDES_SPARSDES_MB                    ",  &
      "ML_DESC_TYPE                              ", "ML_FF_DESC_TYPE                           ",  &
      "ML_DESC_RATIO_DUAL                        ", "ML_FF_DESC_RATIO_DUAL                     ",  &
      "ML_DESC_FACTORE_TESTE                     ", "ML_FF_DESC_FACTORE_TESTE                  ",  &
      "ML_BASIS_TYPE2                            ", "ML_FF_BASIS_TYPE2                         ",  &
      "ML_BASIS_TYPE2_RED                        ", "ML_FF_BASIS_TYPE2_RED                     ",  &

      ! Kernel:
      !========================================================================
      "ML_W1                                     ", "ML_FF_W1_MB                               ",  &
      "ML_NHYP                                   ", "ML_FF_NHYP1_MB                            ",  &

      ! Controlling Bayesian criterion and Spilling Factor:
      !========================================================================
      "ML_ICRITERIA                              ", "ML_FF_ICRITERIA                           ",  &
      "ML_IUPDATE_CRITERIA                       ", "ML_FF_IUPDATE_CRITERIA                    ",  &
      "ML_CTIFOR                                 ", "ML_FF_CTIFOR                              ",  &
      "ML_SCLC_CTIFOR                            ", "ML_FF_SCLC_CTIFOR                         ",  &
      "ML_CSIG                                   ", "ML_FF_CSIG                                ",  &
      "ML_CSLOPE                                 ", "ML_FF_CSLOPE                              ",  &
      "ML_CX                                     ", "ML_FF_XMIX                                ",  &
      "ML_CDOUB                                  ", "ML_FF_CDOUB                               ",  &
      "ML_CALGO                                  ", "ML_FF_CALGO                               ",  &
      "ML_ESTBLOCK                               ", "ML_IERR                                   ",  &
      "ML_SF_LIMIT                               ", "                                          ",  &

      ! Sparsification and regression:
      !========================================================================
      "ML_EPS_LOW                                ", "ML_FF_EPS_LOW                             ",  &
      "ML_EPS_REG                                ", "ML_FF_EPS_REG                             ",  &
      "ML_LBASIS_DISCARD                         ", "ML_FF_LBASIS_DISCARD                      ",  &
      "ML_IALGO_LINREG                           ", "ML_FF_IALGO_LINREG                        ",  &
      "ML_ISVD                                   ", "ML_FF_ISVD                                ",  &
      "ML_IREG                                   ", "ML_FF_IREG_MB                             ",  &
      "ML_SIGV0                                  ", "ML_FF_SIGV0_MB                            ",  &
      "ML_SIGW0                                  ", "ML_FF_SIGW0_MB                            ",  &
      "ML_SAVECMAT                               ", "ML_FF_SAVECMAT_MB                         ",  &

      ! Weighting and sampling:
      !========================================================================
      "ML_WTOTEN                                 ", "ML_FF_WTOTEN                              ",  &
      "ML_WTIFOR                                 ", "ML_FF_WTIFOR                              ",  &
      "ML_WTSIF                                  ", "ML_FF_WTSIF                               ",  &
      "ML_MHIS                                   ", "ML_FF_MHIS                                ",  &
      "ML_NMDINT                                 ", "ML_FF_NMDINT                              ",  &
      "ML_IWEIGHT                                ", "ML_FF_IWEIGHT                             ",  &
      "ML_LUSE_NAMES                             ", "ML_FF_LUSE_NAMES                          ",  &

      ! Sizes of static arrays:
      !========================================================================
      "ML_MB                                     ", "ML_FF_MB_MB                               ",  &
      "ML_MB_MIN                                 ", "                                          ",  &
      "ML_MCONF                                  ", "ML_FF_MCONF                               ",  &
      "ML_MCONF_NEW                              ", "ML_FF_MCONF_NEW                           ",  &

      ! Special features:
      !========================================================================
      "ML_LHEAT                                  ", "ML_FF_LHEAT_MB                            ",  &
      "ML_LCOUPLE                                ", "ML_FF_LCOUPLE_MB                          ",  &
      "ML_NATOM_COUPLED                          ", "ML_FF_NATOM_COUPLED_MB                    ",  &
      "ML_ICOUPLE                                ", "ML_FF_ICOUPLE_MB                          ",  &
      "ML_RCOUPLE                                ", "ML_FF_RCOUPLE_MB                          ",  &
      "ML_LEMPPOT                                ", "                                          ",  &
      "ML_EMPPOT_RCUT                            ", "                                          ",  &
      "ML_SRPOT_B0                               ", "                                          ",  &
      "ML_SRPOT_N0                               ", "                                          ",  &
      "ML_SRPOT_S0                               ", "                                          ",  &
      "ML_MOPOT_NM                               ", "                                          ",  &
      "ML_MOPOT_DM                               ", "                                          ",  &
      "ML_MOPOT_RM                               ", "                                          ",  &
      "ML_MOPOT_RKM                              ", "                                          ",  &
      "ML_MOPOT_IJM                              ", "                                          ",  &

      ! Reference energies:
      !========================================================================
      "ML_ISCALE_TOTEN                           ", "ML_FF_ISCALE_TOTEN_MB                     ",  &
      "ML_EATOM_REF                              ", "ML_FF_EATOM                               ",  &

      ! Input/output options:
      !========================================================================
      "ML_OUTPUT_MODE                            ", "ML_FF_OUTPUT_MODE                         ",  &
      "ML_OUTBLOCK                               ", "ML_FF_OUTBLOCK                            ",  &
      "ML_LEATOM                                 ", "ML_FF_LEATOM_MB                           ",  &
      "ML_GRACE_MODEL                            ", "                                          ",  &

      ! Not important, description not used and not yet reviewed:
      !==========================================================
      "ML_IBROAD1                                ", "ML_FF_IBROAD1_MB                          ",  &
      "ML_IBROAD2                                ", "ML_FF_IBROAD2_MB                          ",  &
      "ML_TOTNUM_INSTANCES                       ", "ML_FF_TOTNUM_INSTANCES                    ",  &
      "ML_FILESTREAM_START                       ", "ML_FF_FILESTREAM_START                    ",  &
      "ML_IMAT_SPARS                             ", "ML_FF_IMAT_SPARS                          ",  &
      "ML_LCONF_DISCARD                          ", "ML_FF_LCONF_DISCARD                       ",  &
      "ML_LMLMB                                  ", "ML_FF_LMLMB                               ",  &
      "ML_NDIM_SCALAPACK                         ", "ML_FF_NDIM_SCALAPACK                      ",  &
      "ML_FF_NWRITE                              ", "                                          ",  & !  MAYBE DELETE, SET DEFAULT!?
      "ML_LNMDINT_RANDOM                         ", "ML_FF_LNMDINT_RANDOM                      ",  &
      "ML_LTEST                                  ", "ML_FF_LTEST                               ",  &
      "ML_LTRJ                                   ", "ML_FF_LTRJ                                ",  &
      "ML_LTOTEN_SYSTEM                          ", "ML_FF_LTOTEN_SYSTEM                       ",  &
      "ML_NTEST                                  ", "ML_FF_NTEST                               ",  &
      "ML_INVERSE_SOAP                           ", "ML_FF_INVERSE_SOAP_MB                     ",  &
      "ML_LSUPERVEC                              ", "ML_FF_LSUPERVEC_MB                        ",  &
      "ML_LSIC                                   ", "ML_FF_LSIC                                ",  &
      "ML_LMETRIC1                               ", "ML_FF_LMETRIC1_MB                         ",  &
      "ML_LMETRIC2                               ", "ML_FF_LMETRIC2_MB                         ",  &
      "ML_LVARTRAN1                              ", "ML_FF_LVARTRAN1_MB                        ",  &
      "ML_LVARTRAN2                              ", "ML_FF_LVARTRAN2_MB                        ",  &
      "ML_NMETRIC1                               ", "ML_FF_NMETRIC1_MB                         ",  &
      "ML_NMETRIC2                               ", "ML_FF_NMETRIC2_MB                         ",  &
      "ML_NVARTRAN1                              ", "ML_FF_NVARTRAN1_MB                        ",  &
      "ML_NVARTRAN2                              ", "ML_FF_NVARTRAN2_MB                        ",  &
      "ML_LWINDOW1                               ", "ML_FF_LWINDOW1_MB                         ",  &
      "ML_LWINDOW2                               ", "ML_FF_LWINDOW2_MB                         ",  &
      "ML_IWINDOW1                               ", "ML_FF_IWINDOW1_MB                         ",  &
      "ML_IWINDOW2                               ", "ML_FF_IWINDOW2_MB                         ",  &
      "ML_MSPL1                                  ", "ML_FF_MSPL1_MB                            ",  &
      "ML_MSPL2                                  ", "ML_FF_MSPL2_MB                            ",  &
      "ML_LNORM1                                 ", "ML_FF_LNORM1_MB                           ",  &
      "ML_LNORM2                                 ", "ML_FF_LNORM2_MB                           ",  &
      "ML_NHYP2                                  ", "ML_FF_NHYP2_MB                            ",  &
      "ML_NR1                                    ", "ML_FF_NR1_MB                              ",  &
      "ML_NR2                                    ", "ML_FF_NR2_MB                              ",  &
      "ML_RMETRIC1                               ", "ML_FF_RMETRIC1_MB                         ",  &
      "ML_RMETRIC2                               ", "ML_FF_RMETRIC2_MB                         ",  &
      "ML_RANDOM_SEED                            ", "ML_FF_RANDOM_SEED                         ",  &
      "ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS   ", "ML_FF_LDISCARD_STRUCTURES_NOT_GIVING_BASIS",  &
      "ML_ICUT1                                  ", "ML_FF_ICUT1_MB                            ",  &
      "ML_ICUT2                                  ", "ML_FF_ICUT2_MB                            ",  &
      "ML_NBLOCK                                 ", "ML_FF_NBLOCK_FFN                          "], &
      [nalt, ntags])
   !> Short description of each tag.
   character(len=*), parameter :: descriptions(ntags) = [                                                      &
      ! General, mode of operation:
      !========================================================================
      "Main control tag to enable machine learning                                                          ", & ! ML_LMLFF
      "Machine learning operation mode in strings (supertag)                                                ", & ! ML_MODE
      "Machine learning operation mode                                                                      ", & ! ML_ISTART
      "Machine learning force field type                                                                    ", & ! ML_TYPE
      "Precontraction of weights for fast execution, but no Bayesian error estimation                       ", & ! ML_LFAST
      "Forcing slow algorithm althoug force field was trained with ML_LFAST=.TRUE. No error estimation      ", & ! ML_LFORCESLOW
      "Enable experimental VASPml C++ library, requires build with -Dlibvaspml                              ", & ! ML_LIB
      "Inside VASPml, use C++ Parallel Algorithms code path                                                 ", & ! ML_PALGO
      "Number of MPI tasks sharing memory segments when using MPI shared memory (-Duse_shmem)               ", & ! ML_NCSHMEM

      ! Descriptors:
      !========================================================================
      "Cutoff radius of radial descriptors                                                                  ", & ! ML_RCUT1
      "Gaussian width for broadening the atomic distribution for radial descriptors                         ", & ! ML_SION1
      "Number of radial basis functions for atomic distribution for radial descriptors                      ", & ! ML_MRB1                                   
      "Cutoff radius of angular descriptors                                                                 ", & ! ML_RCUT2
      "Gaussian width for broadening the atomic distribution for angular descriptors                        ", & ! ML_SION2
      "Number of radial basis functions for atomic distribution for angular descriptors                     ", & ! ML_MRB2
      "Distance for summation of summed up descriptor                                                       ", & ! ML_MRB2_MAX
      "Number of radial basis functions for reduced (color)atomic distribution for angular descriptors      ", & ! ML_MRB2_RED
      "Distance for summation of summed up reduced descriptor                                               ", & ! ML_MRB2_MAX_RED
      "Maximum angular momentum quantum number of spherical harmonics used to expand atomic distributions   ", & ! ML_LMAX2
      "Maximum angular momentum for reduced (type color) descriptor                                         ", & ! ML_LMAX2_RED
      "Angular filtering enabled                                                                            ", & ! ML_LAFILT2
      "Angular filtering parameter a_FILT                                                                   ", & ! ML_AFILT2
      "Angular filtering function type                                                                      ", & ! ML_IAFILT2
      "Enable sparsification of angular descriptors                                                         ", & ! ML_LSPARSDES
      "Number of highest eigenvalues relevant in the sparsification algorithm of the angular descriptors    ", & ! ML_NRANK_SPARSDES
      "Desired ratio of selected to all descriptors resulting from the angular descriptor sparsification    ", & ! ML_RDES_SPARSDES
      "Descriptor type (standard, linear-scaling with element types, ...)                                   ", & ! ML_DESC_TYPE
      "Ratio between high-resolution type-independent and low-resolution type-dependent descriptor          ", & ! ML_DESC_RATIO_DUAL
      "Factori di high-resolution type-independent and low-resolution type-dependent descriptor             ", & ! ML_DESC_FACTORE_TESTE
      "Basis set type for angular descriptor                                                                ", & ! ML_BASIS_TYPE2
      "Basis set type for reduced angular descriptor                                                        ", & ! ML_BASIS_TYPE2_RED

      ! Kernel:
      !========================================================================
      "Weight of radial descriptors in the kernel (the angular counterpart is chosen so that the sum is 1.0)", & ! ML_W1
      "Power of the polynomial kernel                                                                       ", & ! ML_NHYP

      ! Controlling Bayesian criterion:
      !========================================================================
      "Enable automatic updating of the Bayesian error estimation threshold during on-the-fly training      ", & ! ML_ICRITERIA
      "Decides whether update of threshold is done in the same MD step or the next MD step                  ", & ! ML_IUPDATE_CRITERIA
      "Bayesian error estimation threshold (initial or static value depending on other settings)            ", & ! ML_CTIFOR
      "Scaling factor for ML_CTIFOR. The interval 0<ML_SCLC_CTIFOR<1 increases num. of local configurations ", & ! ML_SCLC_CTIFOR
      "Standard error parameter required for the automatic update of the Bayesian error estimation threshold", & ! ML_CSIG
      "Slope parameter required for the automatic update of the Bayesian error estimation threshold         ", & ! ML_CSLOPE
      "Additional parameter controlling the update of the Bayesian error estimation threshold               ", & ! ML_CX
      "Factor controlling the occurence of enforced first principle calculations                            ", & ! ML_CDOUB
      "Decides what to use for learning criterion. 0: Bayesian error estimate; 1: Spilling factor           ", & ! ML_CALGO
      "Time steps between error estimations (prediction modes only)                                         ", & ! ML_ESTBLOCK
      "Maximum allowed spilling factor, higher values trigger graceful exit (fast prediction mode only)     ", & ! ML_SF_LIMIT

      ! Sparsification and regression:
      !========================================================================
      "Threshold for the CUR algorithm used in the sparsification of local reference configurations         ", & ! ML_EPS_LOW
      "Convergence criterion for the optimization of parameters within the Bayesian linear regression       ", & ! ML_EPS_REG
      !Specifies whether the calculation is continued if number of local configurations (NB) is larger than ML_MB and basis sets are thrown away to satisfy NB<=ML_NB", & ! ML_LBASIS_DISCARD
      "If ML_MB is exceeded do not stop but instead discard older local reference configurations            ", & ! ML_LBASIS_DISCARD
      "Linear regression algorithm                                                                          ", & ! ML_IALGO_LINREG
      ! ML_FF_IALGO_LINREG=1: Solve L2-norm (ridge regression). ML_IALGO_LINREG=2: Solve L1-norm with QR factorization. ML_INPUT_HANDLE%ML_FF_IALGO_LINREG = 3 : Solve L1-norm with truncated SVD. ML_INPUT_HANDLE%ML_FF_IALGO_LINREG = 4 : SVD + Tihonov regularization
      "Leverage scoring calculation mode (for sparsification of local configurations)                       ", & ! ML_ISVD
      ! ML_ISVD=1: The values are given from the eigenvalues provided by the diagonalization. ML_ISVD=2: These values are provided from the eigenvalues given from SVD
      "Time dependency of regularization parameters                                                         ", & ! ML_IREG
      ! ML_IREG_MB=1: The (initial) precision (ML_FF_SIGV0_MB) and noise (ML_FF_SIGW0_MB) parameters are kept constant. ML_IREG_MB=2: The parameters are optimized
      "Initial regularization parameter (noise, reversed and squared)                                       ", & ! ML_SIGV0
      "Initial regularization parameter (precision, reversed and squared)                                   ", & ! ML_SIGW0
      "Decides whether CMAT is saved so it does not have to be calculated at each evidence approx. step     ", & ! ML_SAVECMAT

      ! Weighting and sampling:
      !========================================================================
      "Scaling weight for total energies in the training data                                               ", & ! ML_WTOTEN
      "Scaling weight for forces in the training data                                                       ", & ! ML_WTIFOR
      "Scaling weight for stresses in the training data                                                     ", & ! ML_WTISIF
      "Number of estimated errors stored to determine the threshold for the Bayesian error                  ", & ! ML_MHIS
      "Minimum number of MD steps between potential training samples                                        ", & ! ML_NMDINT
      "Weighting method for energies, forces and stresses                                                   ", & ! ML_IWEIGHT
      ! ML_IWEIGHT=1: The unnormalized energies, forces and stress tensor training data are divided by the weights determined by the flags ML_WTOTEN (eV/atom), ML_WTIFOR (eV/Angstrom) and ML_WTSIF (kBar), respectively.. ML_IWEIGHT=2: The training data are normalized by using their standard deviations. Then, the normalized energy, forces and stress tensor are multiplied by unitless ML_WTOTEN, ML_WTIFOR and ML_WTSIF, respectively. ML_IWEIGHT=3: Same as ML_FF_IWEIGHT=2 but the training data is divided into individual subsets.
      "Use also structure names for distinction into groups for weighting                                   ", & ! ML_USE_NAMES

      ! Sizes of static arrays:
      !========================================================================
      "Maximum number of local reference configurations of the ML force field                               ", & ! ML_MB
      "Minimum number of local reference configurations required (otherwise training is blocked)            ", & ! ML_MB_MIN
      ! This number is needed for allocation of static arrays, especially the number of columns in the design matrix.
      "Maximum number of structures stored for training                                                     ", & ! ML_MCONF
      ! This number is needed for allocation of static arrays, especially the number of rows in the design matrix.
      "Maximum number of configurations stored temporarily as training candidates                           ", & ! ML_MCONF_NEW

      ! Special features:
      !========================================================================
      "Enable heat flux calculation (output written to ML_HEAT)                                             ", & ! ML_LHEAT
      "Enable thermodynamic integration (t.i.)                                                              ", & ! ML_LCOUPLE
      "Number of atoms whose interaction is controlled by the t.i. coupling parameter                       ", & ! ML_NATOM_COUPLED
      "List of atoms whose interaction is controlled by the t.i. coupling parameter                         ", & ! ML_ICOUPLE
      "Thermodynamic integration coupling parameter                                                         ", & ! ML_RCOUPLE
      "Model functions enabled for thermodynamic integration                                                ", & ! ML_LEMPPOT
      "Cutoff radius of the model potential                                                                 ", & ! ML_EMPPOT_RCUT
      "Height of the soft-repulsive potential used in the model potential                                   ", & ! ML_SRPOT_B0
      "Power parameter used in the model potential                                                          ", & ! ML_SRPOT_N0
      "Radius of the soft-repulsive potential used in the model potential                                   ", & ! ML_SRPOT_S0
      "Number of pairs which are coupled by Morse model potential                                           ", & ! ML_MOPOT_NM
      "Binding strength of Morse model potential                                                            ", & ! ML_MOPOT_DM
      "Equilibrium length of Morse model potential                                                          ", & ! ML_MOPOT_RM
      "Force constant of Morse model potential                                                              ", & ! ML_MOPOT_RKM
      "List of pairs which are coupled by Morse model potential                                             ", & ! ML_MOPOT_IJM

      ! Reference energies:
      !========================================================================
      "Scaling mode of energies                                                                             ", & ! ML_ISCALE_TOTEN
      ! Specifies how to scale the energy data. ML_ISCALE_TOTEN=1: The total energy is scaled to the total energy of the isolated atoms given by ML_EATOM. ML_ISCALE_TOTEN=2: The total energy is scaled to the average of the training data.
      "List of reference energies of isolated atoms for each species in the system                          ", & ! ML_EATOM_REF

      ! Input/output options:
      !========================================================================
      "Controls the verbosity of the output at each MD step when machine learning is used                   ", & ! ML_OUTPUT_MODE
      "Sets the output frequency of file output (only relevant for ML_MODE = run)                           ", & ! ML_OUTBLOCK
      "Enable output of kin. and pot. energy for each atom at each MD time step (written to ML_EATOM)       ", & ! ML_LEATOM
      "GRACE model name or path                                                                             ", & ! ML_GRACE_MODEL

      ! Not important, description not used and not yet reviewed:
      !========================================================================
      "Specifies whether broadening of the atomic distribution is used for the radial descriptor            ", & ! ML_IBROAD1
      ! ML_IBROAD1=1: No broadening. ML_IBROAD1=2: Element-independent broadening.
      "Specifies whether broadening of the atomic distribution is used for the angular descriptor           ", & ! ML_IBROAD2
      ! ML_IBROAD2=1: No broadening. ML_IBROAD2=2: Element-independent broadening.
      "Specifies the number of machine learning instances                                                   ", & ! ML_TOTNUM_INSTANCES
      "Specifies the starting value for the filehandle numbers within Fortran                               ", & ! ML_FILESTREAM_START
      "Specifies the type of matrix for sparsification                                                      ", & ! ML_IMAT_SPARS
      ! ML_IMAT_SPARS=1: Covariance matrix. ML_IMAT_SPARS=2: SOAP matrix.
      !"ML_LCONF_DISCARD                          ", Specifies whether the calculation is continued if number of training structures (NCONF) is larger than ML_MCONF and training structures are thrown away to satisfy NCONF<=ML_MCONF, or not.
      "Discard training structures if the maximum number in storage is exceeded                             ", & ! ML_LCONF_DISCARD
      ! Specifies whether the calculation is continued if number of training structures (NCONF) is larger than ML_MCONF and training structures are thrown away to satisfy NCONF<=ML_MCONF, or not.
      "Turns many body interaction on. Code won't work without this tag                                     ", & ! ML_LMLMB
      ! Turning off can possibly accelerate the MD when only force fields are used.
      "Specifies the dimension of scalapack grid. Should never be touched                                   ", & ! ML_NDIM_SCALAPACK
      "Deprecated, in past it specified contents of the ML_LOGFILE                                          ", & ! ML_NWRITE
      "Decides wether ML_NMDINT times a random number between 0 and 1 is used instead of ML_NMDINT          ", & ! ML_LNMDINT_RANDOM
      "Whether ab initio calculations are executed at every MD step regardless of the estimated error       ", & ! ML_LTEST
      "DESCRIPTION OF ML_LTRJ                                                                               ", & ! ML_LTRJ
      ! ML_LTRJ=.TRUE.: MD trajectory is always obtained using ab initio forces and stress tensors. ML_LTRJ=.FALSE.: MD trajectory is obtained using generated force field when the estimated error is small. Otherwise, ab initio forces and stress tensors are used.
      "DESCRIPTION OF ML_LTOTEN_SYSTEM                                                                      ", & !ML_LTOTEN_SYSTEM
      ! ML_LTOTEN_SYSTEM=.TRUE.:  Total energy of system is used in regression. ML_LTOTEN_SYSTEM=.FALSE.: Energy per atom is used in regression.
      "This parameter determines how often the test ab initio calculation will be executed                  ", & ! ML_NTEST
      ! The test data will be obtained at every ML_NTEST MD step.
      "Decides how the inverse matrix is calculated for ridge regression                                    ", & ! ML_INVERSE_SOAP
      ! ML_INVERSE_SOAP=1: Diagonalization. ML_INVERSE_SOAP=2: LU factorization.
      "Specifies whether super-vector is used for kernel or not                                             ", & ! ML_LSUPERVEC
      "Specifies whether self-interactions are substracted or not                                           ", & ! ML_LSIC
      "Specifies whether metric function is used for radial descriptor or not                               ", & ! ML_LMETRIC1
      "Specifies whether metric function is used for angular descriptor or not                              ", & ! ML_LMETRIC2
      "Specifies whether variable transform is done for radial coordinate of radial descriptor or not       ", & ! ML_LVARTRAN1
      "Specifies whether variable transform is done for radial coordinate of angular descriptor or not      ", & ! ML_LVARTRAN2
      "Specifies the polynomial parameter for the metric function for the radial descriptor                 ", & ! ML_NMETRIC1
      "Specifies the polynomial parameter for the metric function for the angular descriptor                ", & ! ML_NMETRIC2
      "Polynomial parameter for the variable transform of the radial coordinate for the radial descriptor   ", & ! ML_NVARTRAN1
      "Polynomial parameter for the variable transform of the radial coordinate for the angular descriptor  ", & ! ML_NVARTRAN2
      "Specifies wheter a window function is used for the radial descriptor or not                          ", & ! ML_LWINDOW1
      "Specifies wheter a window function is used for the angular descriptor or not                         ", & ! ML_LWINDOW2
      "Specifies the type of window function for radial descriptor if ML_LWINDOW1=.TRUE.                    ", & ! ML_IWINDOW1
      ! ML_IWINDOW1=1: Hann window. ML_IWINDOW1=2: Hamming window. ML_IWINDOW1=3: Blackman window. ML_IWINDOW1=4: Nuttall window. ML_IWINDOW1=5: Blackman-Nutall. ML_IWINDOW1=6: Blackman-Harris.
      "Specifies the type of window function for angular descriptor if ML_LWINDOW2=.TRUE.                   ", & ! ML_IWINDOW2
      ! ML_IWINDOW2=1: Hann window. ML_IWINDOW2=2: Hamming window. ML_IWINDOW2=3: Blackman window. ML_IWINDOW2=4: Nuttall window. ML_IWINDOW2=5: Blackman-Nutall. ML_IWINDOW2=6: Blackman-Harris.
      "Maximum number of radial grid points used in spline-interpolation for radial descriptor              ", & ! ML_MSPL1
      "Maximum number of radial grid points used in spline-interpolation for angular descriptor             ", & ! ML_MSPL2
      "Specifies whether normalization is taken into account for radial descriptor                          ", & ! ML_LNORM1
      "Specifies whether normalization is taken into account for angular descriptor                         ", & ! ML_LNORM2
      "Sets the hyper parameter (power) of SOAP                                                             ", & ! ML_NHYP2
      "Number of radial mesh points for calculating radial descriptor                                       ", & ! ML_NR1
      "Number of radial mesh points for calculating angular descriptor                                      ", & ! ML_NR2
      "Parameter used for metric function for radial descriptor                                             ", & ! ML_RMETRIC1
      "Parameter used for metric function for angular descriptor                                            ", & ! ML_RMETRIC2
      "Sets the random seed for ML_LNMDINT_RANDOM                                                           ", & ! ML_RANDOM_SEED
      "Specifies whether structures not giving new local configurations are discarded for ML_ISTART=3       ", & ! ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS
      ! Specifies whether structures not giving new local configurations are thrown away as trainings structures for ML_ISTART=3
      "Specifies the type of cut off function for the radial descriptor                                     ", & ! ML_ICUT1
      ! ML_ICUT1=1: Behler-Parinello type cutoff function.. ML_ICUT1=2: Miwa-Ohno type cutoff function.
      "Specifies the type of cut off function for the angular descriptor                                    ", & ! ML_ICUT2
      ! ML_ICUT2=1: Behler-Parinello type cutoff function.. ML_ICUT2=2: Miwa-Ohno type cutoff function.
      "Minimum blocksize for printing of ML_FFN file during learning                                        "]   ! ML_NBLOCK
   !> Maximum number of characters for tags.
   integer, parameter                  :: tag_width = len(tags)
   !> Maximum number of characters for tag description.
   integer, parameter                  :: description_width = len(descriptions)
   !> Tag state symbol width.
   integer, parameter                  :: state_symbol_width = 3

   contains

   !===========================================================================
   !> Initialize taglist type.
   !>
   !> All states are set to #TL_UNSET.
   !===========================================================================
   subroutine taglist_init(this)

      class(taglist), intent(inout) :: this

      ! Initialize state array to TL_UNSET.
      allocate(this%state(ntags))
      this%state = TL_UNSET
      ! Initialize alternative index to 1 (preferred tag).
      allocate(this%alt_index(ntags))
      this%alt_index = 1

   end subroutine taglist_init

   !===========================================================================
   !> Get current tag state.
   !===========================================================================
   function taglist_get_state(this, tag) result(state)

      class(taglist),   intent(in) :: this
      !> Tag for which the state should be returned
      character(len=*), intent(in) :: tag
      !> State of tag.
      integer(kind(TL_STATE))      :: state

      state = this%state(taglist_find_index(tag))

   end function taglist_get_state

   !===========================================================================
   !> Set tag state.
   !===========================================================================
   subroutine taglist_set_state(this, tag, in_state, alt_index)

      class(taglist),          intent(inout)          :: this
      !> Tag for which state should be set.
      character(len=*),        intent(in)             :: tag
      !> State which should be set.
      integer(kind(TL_STATE)), intent(in)             :: in_state
      !> Optional alternative tag name index.
      integer,                 intent(in),   optional :: alt_index

      integer :: t

      t = taglist_find_index(tag)
      this%state(t) = in_state
      if (present(alt_index)) this%alt_index(t) = alt_index

   end subroutine taglist_set_state

   !===========================================================================
   !> Get list of currently used deprecated tags and their preferred version.
   !===========================================================================
   subroutine taglist_get_deprecated(this, list)

      class(taglist),           intent(in)                 :: this
      !> Output list of deprecated tags.
      !>
      !> - `list(:,1)` contains deprecated tags.
      !> - `list(:,2)` contains preferred tag names.
      character(len=tag_width), intent(inout), allocatable :: list(:,:)

      integer :: i
      integer :: c

      ! Count how many deprecated tags were found.
      c = 0
      do i = 1, ntags
         if (this%alt_index(i) > 1) c = c + 1
      end do

      if (allocated(list)) deallocate(list)
      allocate(character(len=tag_width) :: list(c, 2))

      ! Fill list with (deprecated, preferred) tags.
      c = 1
      do i = 1, ntags
         if (this%alt_index(i) > 1) then
            list(c, 1) = tags(this%alt_index(i), i)
            list(c, 2) = tags(1, i)
            c = c + 1
         end if
      end do

   end subroutine taglist_get_deprecated

   !===========================================================================
   !> Get tag state symbol (for display e.g. in ML_LOGFILE).
   !===========================================================================
   function taglist_state_symbol_from_tag(this, tag) result(symbol)

      class(taglist),   intent(in)      :: this
      !> Tag for which symbol should be returned.
      character(len=*), intent(in)      :: tag
      !> Returned state symbol.
      character(len=state_symbol_width) :: symbol

      symbol = taglist_state_symbol(this%get_state(tag))

   end function taglist_state_symbol_from_tag

   !=======================!
   ! Non-member procedures !
   !=======================!

   !===========================================================================
   !> Get maximum length of tag name.
   !===========================================================================
   pure function taglist_get_tag_width() result(width)

      integer :: width

      width = tag_width

   end function taglist_get_tag_width

   !===========================================================================
   !> Get description for given tag.
   !===========================================================================
   function taglist_description(tag) result(description)

      !> Tag for which description should be returned.
      character(len=*), intent(in)  :: tag
      !> Returned tag description.
      character(len=:), allocatable :: description

      description = descriptions(taglist_find_index(tag))

   end function taglist_description

   !===========================================================================
   !> Find index of given tag in list of all ML tags.
   !>
   !> Stops if input is not in list of all tags.
   !===========================================================================
   function taglist_find_index(tag) result(tag_index)

      use iso_fortran_env, only: error_unit
      !> Tag for which index should be returned.
      character(len=*), intent(in) :: tag
      !> Returned tag index.
      integer(kind(TL_STATE))      :: tag_index

      integer :: t ! Tag loop index.
      integer :: a ! Alternative tag loop index.

      tag_loop: &
      do t = 1, ntags
         do a = 1, nalt
            ! Avoid comparing with potential empty alternative string, as this
            ! is asking for trouble...
            if (a > 1 .and. adjustl(trim(tags(a, t))) == "") cycle
            if (adjustl(trim(tag)) == adjustl(trim(tags(a, t)))) exit tag_loop
         end do
      end do tag_loop
      if (t >= ntags + 1) then
         write(error_unit, '(A,1X,A)') "ERROR: Unknown tag: ", tag
         stop 1
      end if
      tag_index = t

   end function taglist_find_index

   !===========================================================================
   !> Get list of all preferred tags.
   !===========================================================================
   function taglist_get_all_tags() result(list)

      !> Returned list of tags.
      character(len=:), allocatable :: list(:)

      list = tags(1, :)

   end function taglist_get_all_tags

   !===========================================================================
   !> Get list of alternative tags for given tag.
   !>
   !> Input tag can be preferred or alternative tag.
   !===========================================================================
   function taglist_get_alt_tags(tag) result(list)

      !> Input tag of interest.
      character(len=*)              :: tag
      !> Returned list of all versions of input tag.
      character(len=:), allocatable :: list(:)

      integer :: t ! Tag loop index.
      integer :: a ! Alternative tag loop index.
      integer :: n ! Number of non-empty alternative strings.
      integer :: c ! Counter.

      t = taglist_find_index(tag)
      n = 0
      do a = 1, nalt
         if (.not. (adjustl(trim(tags(a, t))) == "")) n = n + 1
      end do
      allocate(character(len=tag_width) :: list(n))
      c = 1
      do a = 1, nalt
         if (.not. (adjustl(trim(tags(a, t))) == "")) then
            list(c) = tags(a, t)
            c = c + 1
         end if
      end do

   end function taglist_get_alt_tags

   !===========================================================================
   !> Get maximum length of tag description.
   !===========================================================================
   pure function taglist_get_description_width() result(width)

      !> Returned description width.
      integer :: width

      width = description_width

   end function taglist_get_description_width

   !===========================================================================
   !> Get state symbol for given state.
   !===========================================================================
   pure function taglist_state_symbol(state) result(symbol)

      !> Input state of interest.
      integer(kind(TL_STATE)), intent(in) :: state
      !> Returned state symbol.
      character(len=state_symbol_width)   :: symbol

      select case(state)
         case(TL_UNSET)
            symbol = "(?)"
         case(TL_DEFAULT)
            symbol = "   "
         case(TL_INCAR)
            symbol = "(I)"
         case(TL_INCAR_ALT)
            symbol = "(i)"
         case(TL_OVERRIDE)
            symbol = "(!)"
      end select

   end function taglist_state_symbol

end module ml_ff_taglist
#endif
