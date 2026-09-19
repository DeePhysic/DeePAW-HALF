#include "symbol.inc"
!************************************************************************
!
!>  This module contains all the types used in the machine-learning force field part.
!
!***********************************************************************
      MODULE ML_FF_STRUCT
      USE ML_FF_PREC
#if defined(MPI) || defined(MPI_CHAIN)
      USE ML_FF_TAGLIST
      USE VERLET_CELL_LIST_MOD, only : VERLET_CELL_LIST_TYPE
      USE MPI_HELP, ONLY : ML_MPI_PAR
#endif
#ifdef use_shmem
#ifndef sysv
      USE mpi_f08, ONLY : MPI_Win
#endif
      USE iso_c_binding
#endif
      IMPLICIT NONE

!-------------------------------------------------------------------------
!
!> This structure contains all the input variables read in and also some
!> related variables.
!
!-------------------------------------------------------------------------
      TYPE ML_FF_INPUT
         !> real space primitive lattice vectors
         REAL(q)                       :: A(1:3,1:3)        = 0.0_q
         !> reciprocal space lattice vectors
         REAL(q)                       :: B(1:3,1:3)        = 0.0_q
         !> communicator passed from VASP. Usually COMM_WORLD.
         INTEGER                       :: COMM_VASP         = 0
         !> location of directory passed from VASP
         CHARACTER(LEN=10)             :: DIR_APP           = ""
         !> length of dir_app
         INTEGER                       :: DIR_LEN           = 0
         !> this arrays specify the atomic energies in isolated states.
         REAL(q),          ALLOCATABLE :: EATOM_VASP(:)
         LOGICAL                       :: LABORT            = .FALSE.
         !> maximum number of atoms per species
         INTEGER                       :: MITYP             = 0
         !> the number of MD steps
         INTEGER                       :: NSTEP             = 0
         !> the number of element types
         INTEGER                       :: NTYP              = 0
         !> the number of atoms for each species
         INTEGER,          ALLOCATABLE :: NITYP(:)
         !> the number of atoms
         INTEGER                       :: NIONS             = 0
         !> this variable specifies the dimension of scaLAPACK grid
         INTEGER                       :: NDIM_SCALAPACK_FF = 0
         !> this variable specifies the contents of ML_LOGFILE.
         INTEGER                       :: NWRITE_FF         = 0
         !> cell volume
         REAL(q)                       :: OMEGA             = 0.0_q
         !> the mass of atoms
         REAL(q),          ALLOCATABLE :: POMASS(:)
         !> atomic positions
         REAL(q),          ALLOCATABLE :: POSION(:,:)
         !> the length of each MD step
         REAL(q)                       :: POTIM             = 0.0_q
         !> internal name of current model system
         CHARACTER(LEN=40)             :: SZNAM2            = ""
         !> name of current model system
         CHARACTER(LEN=40)             :: SZNAM2_INPUT      = ""
         !> total energy, force and stress tensor provided by ab initio calculations.
         REAL(q)                       :: TOTEN_ML          = 0.0_q
         REAL(q),          ALLOCATABLE :: TIFOR_ML(:,:)
         REAL(q),          ALLOCATABLE :: TSIF_ML(:,:)
         !> type information for each species
         CHARACTER(LEN=2), ALLOCATABLE :: TYPE(:)
         !> velocity.
         REAL(q),          ALLOCATABLE :: VEL_ML(:,:)
         !> Weights used for scaling the training data.
         REAL(q)                       :: WTIFOR_FF         = 0.0_q
         REAL(q)                       :: WTOTEN_FF         = 0.0_q
         REAL(q)                       :: WTSIF_FF          = 0.0_q
      END TYPE

!-------------------------------------------------------------------------------------
!
!> This structure contains all the variables for ab-initio stuff (mainly
! from the ML_AB file.
!
!-------------------------------------------------------------------------------------

      TYPE ABINITIO
         REAL(q),           ALLOCATABLE :: A(:,:,:)               !< Lattice vectors. 1d:,2d:,3d:
         CHARACTER(LEN=9)               :: ABCAR          = ""    !< Ab initio data file name.
         CHARACTER(LEN=10)              :: ABNCAR         = ""    !< New ab initio data file name.
         REAL(q),           ALLOCATABLE :: B(:,:,:)               !< Reciprocal lattice vectors
         REAL(q),           ALLOCATABLE :: CTIFOR_ALLCONF_OLD(:)  !< old CTIFOR (Bayesian error criterion) for each training configuration
         REAL(q),           ALLOCATABLE :: CTIFOR_ALLCONF_NEW(:)  !< new CTIFOR (Bayesian error criterion) for each training configuration
         LOGICAL,           ALLOCATABLE :: CTIFOR_PRESENT(:)      !< Checks if the
         REAL(q),           ALLOCATABLE :: EATOM(:)               !< Atomic reference energies read from ML_AB file
         CHARACTER(LEN=10)              :: HISCAR         = ""    !< Histogram file
         INTEGER,           ALLOCATABLE :: ISYS(:)                !< Model system number.
         INTEGER                        :: IU1            = 0     !< Code number of ab initio data file.
         INTEGER                        :: IU2            = 0     !< Code number of new ab initio data file.
         INTEGER                        :: IU3            = 0     !< Code number of output file summarizing the regression results
         INTEGER                        :: IU4            = 0     !< Code number of histogram file
         INTEGER,           ALLOCATABLE :: LADD(:,:)              !< Address of ion type in FF
         INTEGER,           ALLOCATABLE :: LB(:,:)                !< Label of basis set
         INTEGER,           ALLOCATABLE :: LNCONF(:,:)            !< Conf number of each basis set - needed only for helping arrays
         INTEGER,           ALLOCATABLE :: LNITYP(:,:)            !< List of atom number
         INTEGER,           ALLOCATABLE :: LNIONS(:,:,:)          !< List of total atom number
         INTEGER,           ALLOCATABLE :: LNIONS_FF(:,:)         !< Ion number of each basis set - needed only for helping arrays
         INTEGER,           ALLOCATABLE :: LNTYP(:,:)             !< List of element number
         INTEGER,           ALLOCATABLE :: LSYS(:)                !< Variable to identify whether the stored configuration is the same as the system handled by VASP currently
         INTEGER                        :: MCONF          = 0     !< Maximum number of configuration
         INTEGER                        :: MCONF_ALLOCATE = 0
         INTEGER                        :: MITYP          = 0     !< Maximum number of atoms per specie
         INTEGER                        :: MITYP_ALLOCATE = 0     !< Size of MITYP used for allocation
         INTEGER,           ALLOCATABLE :: MITYP_CONF(:)
         INTEGER                        :: MIONS          = 0     !< Maximum number of atoms per system
         INTEGER                        :: MIONS_ALLOCATE = 0     !< Size of MIONS used for allocation
         INTEGER,           ALLOCATABLE :: MIONS_CONF(:)
         INTEGER                        :: MTYP           = 0     !< Maximum number of species
         INTEGER                        :: MTYP_ALLOCATE  = 0     !< Size of MTYP used for allocation
         INTEGER,           ALLOCATABLE :: MTYP_CONF(:)
         INTEGER,           ALLOCATABLE :: NBASIS(:)              !< The number of basis sets per atom type, needed only for helping arrays
         INTEGER                        :: NB_MAX         = 0     !< Maximum number of possible basis sets on the ML_AB file
         INTEGER                        :: NCONF          = 0     !< Number of configuration
         INTEGER,           ALLOCATABLE :: NCONF_SYS(:)           !< The number of configuration per model system
         INTEGER,           ALLOCATABLE :: NITYP(:,:)             !< Number of atoms for each specie
         INTEGER,           ALLOCATABLE :: NIONS(:)               !< Total number of atom
         INTEGER                        :: NSYS           = 0
         INTEGER                        :: NSYS_ALLOCATE  = 0
         INTEGER,           ALLOCATABLE :: NTYP(:)                !< Number of species
         REAL(q),           ALLOCATABLE :: POMASS(:)              !< Atomic mass
         REAL(q),           ALLOCATABLE :: POSION(:,:,:)          !< Atom positions
         CHARACTER(LEN=10)              :: REGCAR         = ""    !< Output file summarizing the regression results.
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2(:)              !< Name of model system.
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_ADD(:)          !< Name of model system.
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_INPUT(:)        !< Name of model system.
         REAL(q),           ALLOCATABLE :: TIFOR(:,:,:)           !< Ab initio atom forces
         REAL(q),           ALLOCATABLE :: TIFOR_MB(:,:,:)        !< Machine-learning many-body interaction force
         REAL(q),           ALLOCATABLE :: TOTEN(:)               !< Ab initio total energies
         REAL(q),           ALLOCATABLE :: TOTEN_MB(:)            !< Machine-learning many-body interaction energy
         REAL(q),           ALLOCATABLE :: TSIF(:,:,:)            !< Ab initio stress
         REAL(q),           ALLOCATABLE :: TSIF_MB(:,:,:)         !< Machine-learning many-body interaction stress
         CHARACTER(LEN=2),  ALLOCATABLE :: TYPE(:)                !< Element symbols
         CHARACTER(LEN=2),  ALLOCATABLE :: TYPE_CONF(:,:)         !< Element symbols for each type within each local reference configuration
         REAL(q),           ALLOCATABLE :: OMEGA(:)               !< Volme of cell
         REAL(q)                        :: OMEGA_MAX      = 0.0_q !< Maximum volume amongst training structures
         REAL(q)                        :: VERSION        = 0.0_q !< Version number of ML_AB file
      END TYPE

      TYPE ABINITIONEW
         REAL(q),           ALLOCATABLE :: A(:,:,:)
         REAL(q),           ALLOCATABLE :: B(:,:,:)
         REAL(q),           ALLOCATABLE :: BE(:,:)
         REAL(q),           ALLOCATABLE :: BEAV(:,:)          !< The average error in force estimated by Bayesian theorem.
         REAL(q),           ALLOCATABLE :: BEMAX(:,:)         !< The maximum error in force estimated by Bayesian theorem.
         REAL(q),           ALLOCATABLE :: BEAV_TYPE(:,:,:)   !< Av. Bayesian force error. Careful: ordering of type is according to ff.
         REAL(q),           ALLOCATABLE :: BEMAX_TYPE(:,:,:)  !< Max Bayesian force error. Careful: ordering of type is according to ff.
         REAL(q),           ALLOCATABLE :: CTIFOR_ALLCONF(:)
         REAL(q),           ALLOCATABLE :: EATOM_KIN(:)       !< Atomic kinetic energies. If FFM%LEATOM=.TRUE. or FFM%LHEAT=.TRUE., they are calculated.
         REAL(q),           ALLOCATABLE :: EATOM_POT(:)       !< Atomic potential energies. If FFM%LEATOM=.TRUE. or FFM%LHEAT=.TRUE., they are calculated.
         INTEGER                        :: ISYS           = 0 !< Model system number.
         INTEGER,           ALLOCATABLE :: LADD(:)
         INTEGER,           ALLOCATABLE :: LNITYP(:)
         INTEGER,           ALLOCATABLE :: LNIONS(:,:)
         INTEGER,           ALLOCATABLE :: LNTYP(:)
         INTEGER,           ALLOCATABLE :: LSYS(:)
         INTEGER                        :: MCONF          = 0 !< If the number of new configuration reaches this number, execute force field generation
         INTEGER                        :: MITYP          = 0
         INTEGER                        :: MITYP_ALLOCATE = 0
         INTEGER                        :: MIONS          = 0
         INTEGER                        :: MIONS_ALLOCATE = 0
         INTEGER                        :: MTYP           = 0
         INTEGER                        :: MTYP_ALLOCATE  = 0
         INTEGER                        :: NCONF              !< The number number of new configuration
         INTEGER,           ALLOCATABLE :: NITYP(:)
         INTEGER                        :: NIONS          = 0
         INTEGER,           ALLOCATABLE :: NSTEP(:)
         INTEGER                        :: NSYS           = 0
         INTEGER                        :: NSYS_ALLOCATE  = 0
         INTEGER                        :: NTYP           = 0
         REAL(q),           ALLOCATABLE :: POMASS(:)
         REAL(q),           ALLOCATABLE :: POSION(:,:,:)
         REAL(q)                        :: POTIM          = 0
         REAL(q),           ALLOCATABLE :: QHEAT(:,:)
         REAL(q)          , ALLOCATABLE :: SF(:,:)
         REAL(q)          , ALLOCATABLE :: SFMAX(:)
         REAL(q)          , ALLOCATABLE :: SFMAX_TYPE(:,:)      !< MAx SF per type. Careful: ordering of type is according to ff. 
         REAL(q)          , ALLOCATABLE :: SFMEAN(:)
         REAL(q)          , ALLOCATABLE :: SFMEAN_TYPE(:,:)     !< Mean SF per typ. Careful: ordering of type is according to ff.
         REAL(q)          , ALLOCATABLE :: SFMIN(:)
         REAL(q)          , ALLOCATABLE :: SFVAR(:)
         CHARACTER(LEN=40)              :: SZNAM2         = ""
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_ADD(:)
         CHARACTER(LEN=40)              :: SZNAM2_INPUT   = ""
         REAL(q),           ALLOCATABLE :: TIFOR(:,:,:)
         REAL(q),           ALLOCATABLE :: TIFOR_MB(:,:,:)
         REAL(q),           ALLOCATABLE :: TOTEN(:)
         REAL(q),           ALLOCATABLE :: TOTEN_MB(:)
         REAL(q),           ALLOCATABLE :: TSIF(:,:,:)
         REAL(q),           ALLOCATABLE :: TSIF_MB(:,:,:)
         CHARACTER(LEN=2),  ALLOCATABLE :: TYPE(:)
         INTEGER,           ALLOCATABLE :: TYPE_OF_EACH_ATOM_IN_ABN(:)
         REAL(q),           ALLOCATABLE :: OMEGA(:)
         REAL(q),           ALLOCATABLE :: VEL(:,:,:)
      END TYPE

!--------------------------------------------------------------------------------------
!
!> This is the handle for force field variables.
!
!--------------------------------------------------------------------------------------

      TYPE FF_PAR
         REAL(q),           ALLOCATABLE :: BEHIS(:)                                       !< History of Bayesian error.
         !> If the maximum Bayesian error is larger than CDOUB*CTIFOR sampling is stopped, and force field is calculated.
         REAL(q)                        :: CDOUB                                = 0.0_q
         REAL(q)                        :: CSIG                                 = 0.0_q   !< A criteia for variance to judge the convergence.
         REAL(q)                        :: CSIG_LOG                             = 0.0_q   !< Current criteria sigma, for logging only.
         REAL(q)                        :: CSLOPE                               = 0.0_q   !< A criteria for slope to judge the convergence.
         REAL(q)                        :: CSLOPE_LOG                           = 0.0_q   !< Current criteria slope, for logging only.
         REAL(q)                        :: CTIFOR                               = 0.0_q   !< A criteria for force. 
         REAL(q)                        :: CTIFOR_PREV                          = 0.0_q   !< The threshold before the last update,
                                                                                          !< used for logging.
         REAL(q)                        :: SCALE_CTIFOR                         = 0.0_q   !< Scaling factor for the threshold.
         REAL(q)                        :: CX                                   = 0.0_q   !< Parameter x for criteria update.
         LOGICAL                        :: DELTADFT                             = .FALSE. !< Decides whether delta is run or not
         CHARACTER(LEN=10)              :: DIR_APP                              = ""      !< location of directory passed from VASP
         INTEGER                        :: DIR_LEN                              = 0       !< location of directory passed from VASP
         REAL(q),           ALLOCATABLE :: EATOM_REF(:)                                   !< Reference energy for atoms.
         REAL(q),           ALLOCATABLE :: ERRHIS(:)
         CHARACTER(LEN=9)               :: FFCAR                                = ""      !< Input file summarizing old FF parameters. If ISTART=2, this file is read.
         CHARACTER(LEN=15)              :: FFNCAR                               = ""      !< Output file summarizing new FF parameters
         LOGICAL                        :: FORCING_UPDATE                       = .FALSE. !< Forces learning, creation of new force field and eventually ab initio calculation
         CHARACTER(LEN=:), ALLOCATABLE  :: GRACE_MODEL                                    !< GRACE model name or path.
         !> Algorithm for linear regression
         !>
         !> IALGO_LINREG=1, L2-norm with ridge regression
         !> IALGO_LINREG=2, L1-norm with QR factorization
         !> IALGO_LINREG=3, L1-norm with truncated SVD
         !> IALGO_LINREG=4, SVD + Tikhonov regularization
         INTEGER                        :: IALGO_LINREG                         = 0
         INTEGER                        :: ICRITERIA                            = 0       !< Parameter controlling update of threshold for Bayesian learning
         INTEGER                        :: IERR                                 = 0       !< Frequency for error estimation of BEE
         INTEGER                        :: IERR_COUNTER                         = 0       !< Counter since last calculation of BEE
         INTEGER                        :: CALGO                                = 0       !< Chooses error criterion type
         !> CALGO=0 : BEE used
         !> CALGO=1 : Spilling factor used
         LOGICAL                        :: IFF                                  = .FALSE. !< Flag for judging the presence of FF.
         LOGICAL                        :: IFF_OLD                              = .FALSE. !< Old IFF data.
         !> Type of matrix used the sparsification.
         !>
         !> IMAT_SPARS=1: Covariance matrix.
         !> IMAT_SPARS=2: SOAP matrix.
         INTEGER                        :: IMAT_SPARS                           = 0
         !> Controls the verbosity of the output
         !>
         !> OUTPUT_MODE = 0: No XML, XDATCAR
         INTEGER                        :: OUTPUT_MODE                          = 0
         INTEGER                        :: OUTBLOCK                             = 0       !< Sets the output frequency for ML_ISTART=2
         LOGICAL                        :: ISPARS                               = .FALSE. !< ISPARS=.TRUE. means the successful execution of sparsification.
         !> ISTART determines the mode of operation:
         !>
         !> ISTART = 0 : Start the calculation without any previous data.
         !> ISTART = 1 : Restart the calculation after reading the previous database file.
         !> ISTART = 2 : Only force field is used without learning.
         !> ISTART = 3 : New force-field created only from ML_AB file.
         INTEGER                        :: ISTART                               = 0
         INTEGER                        :: IU1                                  = 0       !< File code number for the old FF file
         INTEGER                        :: IU2                                  = 0       !< File code number for the new FF file
         INTEGER                        :: IU3                                  = 0       !< File code number for the local atomic energies
         INTEGER                        :: IU4                                  = 0       !< File code number for the local heat flux.
         !> Type of the SVD for the sparsification.
         !>
         !> ISVD=1: Diagonalization.
         !> ISVD=2: Singular value decomposition.
         INTEGER                        :: ISVD                                 = 0
         INTEGER                        :: IUPDATE_CRITERIA                     = 0       !< Decides whether Bayesian criterion is updated in same step or next
         !> This variable specifies the way to scale the training data.
         !>
         !> If IWEIGHT = 1, scale the training data by using constants, WTIFOR, WTOTEN and WTSIF that are
         !> provided from INCAR file.
         !> If IWEIGHT = 2, scale the training data by using variances in all training data.
         !> If IWEIGHT = 3, scale the training data by using averages of variances in system-dependent training data.
         !> If IWEIGHT = 4, scale the training data by using variances in all training data.
         !> In this scaling, the number of structure datasets are also taken into account.
         !> If IWEIGHT = 5, scale the training data by using averages of variances in
         !> system-dependent training data. In this scaling, the number of structure datasets are also taken into account.
         LOGICAL                        :: LIB                                  = .FALSE. !< Enable VASPml if available (-Dlibvaspml).
         INTEGER                        :: IWEIGHT                              = 0
         LOGICAL                        :: LABINITIO                            = .FALSE. !< Shows whether ab initio calculations 
                                                                                          !< are provided this turn or not
         INTEGER,           ALLOCATABLE :: LADD(:)                                        !< Address of element
         LOGICAL                        :: LBASIS_DISCARD                       = .FALSE. !< If LBASIS_DISCARD=.TRUE., 
                                                                                          !< basis sets are discarded when its number 
                                                                                          !< exceeds MB_FFM.
         LOGICAL                        :: LCONF_DISCARD                        = .FALSE. !< If LCONF_DISCARD=.TRUE., 
                                                                                          !< configurations that do not provide any 
                                                                                          !< local reference configurations are discarded.
         LOGICAL                        :: LCTIFOR_PRESENT_IN_MLAB              = .FALSE. !< True if all structures in ML_AB file 
                                                                                          !< contain a CTIFOR value.
         LOGICAL                        :: LDISCARD_STRUCTURES_NOT_GIVING_BASIS = .FALSE.
         LOGICAL                        :: LTHRESHOLD                           = .FALSE. !< If .TRUE., marks that the CTIFOR 
                                                                                          !< threshold was exceeded (used for logging).
         LOGICAL                        :: LDOUB                                = .FALSE. !< If LDOUB=.TRUE., the sampling is stopped, 
                                                                                          !< and a new force field is generated.
         LOGICAL                        :: LFAST                                = .FALSE. !< Fast execution mode for ISTART=2 on or not
         LOGICAL                        :: LFAST_FROM_FILE                      = .FALSE. !< Checks if LFAST was read from ML_LOGFILE
         LOGICAL                        :: LFAST_READ                           = .FALSE. !< Value of LFAST read from the ML_FF file
         LOGICAL                        :: LFMAT                                = .FALSE. !< If LFMAT=.TRUE., design matrix elements 
                                                                                          !< are calculated
         LOGICAL                        :: LFORCESLOW ! To force slow algorithm, if LFAST=.TRUE. was used for training
         LOGICAL                        :: LGENFF                               = .FALSE. !< LGENFF=.TRUE. means that force field 
                                                                                          !< generation is executed.
         LOGICAL                        :: LLIST                                = .FALSE. !< If LLIST=.TRUE., the machine generates lists 
                                                                                          !<storing the locations of design-matrix elements.
         !> If LLIST_NEW=.TRUE., the machine generates lists for storing the locations of design-matrix elements related to the candidates of new training data.
         LOGICAL                        :: LLIST_NEW                            = .FALSE.
         LOGICAL                        :: LLOG_THRESHOLD_HISTORY               = .FALSE. !< If .TRUE. write threshold history to log 
                                                                                          !< file on next occasion.
         LOGICAL                        :: LLOG_THRESHOLD_UPDATE                = .FALSE. !< If .TRUE. write threshold update info to 
                                                                                          !< log file on next occasion.
         LOGICAL                        :: LMLABEXIST ! Variable neccessary for checking of existance of ML_AB file
         LOGICAL                        :: LMLFFEXIST ! Variable neccessary for checking of existance of ML_FF file
         LOGICAL                        :: LMLONLY                              = .FALSE. !< If LMLONLY=.TRUE., only the predictions 
                                                                                          !< by the machine-learning force field 
                                                                                          !< will be executed.
         !> If LMLMB=.TRUE., the program generates force fields including many-body terms.
         !>
         !> Actually, current program needs LMLMB=.TRUE.
         LOGICAL                        :: LMLMB                                = .FALSE.
         LOGICAL                        :: LNMDINT_RANDOM                       = .FALSE. !< Decids wether NMDINT times a random number between 0 and 1 is used instead of NMDINT
         LOGICAL                        :: LOVER                                = .FALSE. !< LOVER=.TRUE. means that sparsification is executed to prevent the overcompleteness.
         LOGICAL                        :: LPS                                  = .FALSE. !< LPS=.TRUE. means that angular descriptor is calculated.
         LOGICAL                        :: LRUN_VASP_TO_ML                      = .FALSE. !< Decides whether ML_TO_VASP is run or not.
         LOGICAL                        :: LREAD_SINV_FROM_FILE                 = .FALSE. !< Decides whether to read SINV from file
         LOGICAL                        :: LSAMPLE                              = .FALSE. !< If LSAMPLE=.TRUE., sampling is executed.
         LOGICAL                        :: LSECOND_TIME_ML_TO_VASP              = .FALSE. !< Helping variable to prevent incrementation of control variables
         LOGICAL                        :: LSIG                                 = .FALSE. !< If LSIG=.TRUE., the machine calculates the variance of training data.
         LOGICAL                        :: LSOAP                                = .FALSE. !< If LSOAP=.TRUE., SOAP matrix elements are calculated.
         LOGICAL                        :: LTEST                                = .FALSE. !< If LTEST=.TRUE., ab initio calculation is executed at every MD step regardless the estimated error.
         LOGICAL                        :: LTOTEN_SYSTEM                        = .FALSE. !< IF LTOTEN_SYSTEM=.TRUE. total energy of system is learned instead of energy/atom
         LOGICAL                        :: LTRJ                                 = .FALSE. !< If LTRJ=.TRUE., MD is executed using ab initio forces.
         LOGICAL                        :: LUSE_NAMES                           = .FALSE. !< If LUSE_NAMES=.TRUE., names of structures are also considered within grouping process
         INTEGER                        :: MCONF                                = 0       !< The number of maximum training configurations
         INTEGER                        :: MCONF_INIT                           = 0       !< MCONF value provided by INCAR file
         INTEGER                        :: MCONF_NEW                            = 0       !< The number of maximum new ab initio configurations
         INTEGER                        :: MHIS                                 = 0       !< The maximum number of history of Bayesian errors stored in the memory.
         INTEGER                        :: MITYP                                = 0       !< Maximum number of atoms for each species
         INTEGER                        :: MITYP_ALLOCATE                       = 0       !< Size of MITYP used for allocation
         INTEGER                        :: MIONS                                = 0       !< Maximum number of atoms for each systeem
         INTEGER                        :: MIONS_ALLOCATE                       = 0       !< Size for allocation of MIONS
         CHARACTER(LEN=12)              :: ML_EATOM                             = ""      !< Output file containing the local atomic energies.
         REAL(q)                        :: ML_EPS_LOW                           = 0.0_q   !< Error tolerance EPS_LOW in SPARS_ICUR
         REAL(q)                        :: ML_EPS_REG                           = 0.0_q   !< Parameter for regularization in Bayesian linear regression in BLEA_MB
         CHARACTER(LEN=11)              :: ML_HEAT                              = ""      !< Output file containing the local heat flux.
         LOGICAL                        :: ML_FF_BINARY                         = .FALSE. !< Keeps track whether ML_FF file is in binary format or not
         LOGICAL                        :: ML_FF_HAS_HEADER                     = .FALSE. !< If the ML_FF file contains an ASCII header.
         INTEGER                        :: ML_FF_VERSION_READ(1:3)              = 0       !< Version number of ML_FF file read from ML_FF file, x.y.z format (x: major, y: minor, z:patch)
         CHARACTER(LEN=23)              :: ML_FF_CREATION_DATETIME              = ""      !< Creation date and time string of last written ML_FFN file.
         CHARACTER(LEN=:), ALLOCATABLE  :: ML_TYPE                                        !< ML force field type (kernel-based, grace, ...).
         CHARACTER(LEN=:), ALLOCATABLE  :: MODE                                           !< ML operation mode (it is a super variable)
         INTEGER                        :: MTYP                                 = 0       !< The number of atom type in force field file
         INTEGER                        :: MTYP_ALLOCATE                        = 0       !< Size for MTYP used for allocation
         INTEGER                        :: NCSHMEM                              = 0       !< number of compute core sharing memory for shared
                                                                                          !< memory mpi
         INTEGER                        :: NHIS                                 = 0       !< The number of history of Bayesian errors stored in the memory.
         INTEGER                        :: NBLOCK_FFN                           = 0       !< Block size for which ML_FFN file is printed out
         INTEGER                        :: NCOUNT_GENFF                         = 0       !< Counter how often force field was updated since last print of ML_FFN
         INTEGER                        :: NMDINT                               = 0       !< Interval of no-sampling.
         INTEGER                        :: NMDINT_RANDOM                        = 0       !< NMDINT multiplied with a random number
         INTEGER                        :: NSTEP                                = 0       !< step number of the MD
         INTEGER                        :: NSTEP_FROM_GENFF                     = 0       !< The number of MD step from the previous FF generation.
         INTEGER                        :: NSTEP_FROM_SAMPLING                  = 0       !< The number of MD step from the previous sampling.
         INTEGER                        :: NSYS                                 = 0       !< The number of different systems in the traiing dataset.
         INTEGER                        :: NSYS_ALLOCATE                        = 0
         INTEGER                        :: NSYS_AB_ONLY                         = 0       !< Number of different systems only from AB file
         INTEGER                        :: NSW                                  = 0       !< Number of force field steps
         INTEGER                        :: NTEST                                = 0       !< This parameter determines how often the test ab initio calculation will be executed.
         INTEGER                        :: NTIFOR                               = 0       !< The number of ab initio force data.
         INTEGER                        :: NTOTEN                               = 0       !< The number of ab initio energy data.
         INTEGER                        :: NTSIF                                = 0       !< The number of ab initio stress tensor data.
         CHARACTER(LEN=:), ALLOCATABLE  :: PALGO                                          !< Use C++ Parallel Algorithms code path in VASPml.
         REAL(q),           ALLOCATABLE :: POMASS(:)                                      !< Mass of atom
         REAL(q)                        :: RELTOTEN                             = 0.0_q   !< Rel error energy
         REAL(q)                        :: RELTOTEN_LAST                        = 0.0_q   !< Rel error energy last configuration
         REAL(q)                        :: RELTIFOR                             = 0.0_q   !< Rel error forces
         REAL(q)                        :: RELTIFOR_LAST                        = 0.0_q   !< Rel error forces last configuration
         REAL(q),           ALLOCATABLE :: RELTIFOR_TYPE(:)                               !< Rel error forces for each type
         REAL(q)                        :: RELTSIF                              = 0.0_q   !< Rel error stress
         REAL(q)                        :: RELTSIF_LAST                         = 0.0_q   !< Rel error stress last configuration
         INTEGER                        :: SEED(1:3)                            = 0       !< Random seed used for random numbers
         REAL(q)                        :: SF_LIMIT                             = 0.0_q   !< Maximum allowed spilling factor, higher values trigger graceful exit.
         LOGICAL                        :: SF_LIMIT_EXCEEDED                    = .FALSE. !< Signals that the maximum spilling factor was exceeded.
         REAL(q)                        :: SIG(1:10)                            = 0.0_q   !< Variance in the ab initio data.
         REAL(q),           ALLOCATABLE :: SIG_SYS(:,:)                                   !< Variance in the ab initio data for individual system.
         REAL(q)                        :: STIFOR                               = 0.0_q   !< Root mean square error in predicted force.
         REAL(q)                        :: STIFOR_LAST                          = 0.0_q   !< Root mean square error in predicted 
                                                                                          !< force for last added configuration
         REAL(q),           ALLOCATABLE :: STIFOR_TYPE(:)                                 !< Root mean square error in predicted                                                                                                          !< force for each type
         REAL(q)                        :: STOTEN                               = 0.0_q   !< Root mean square error in predicted total energy.
         REAL(q)                        :: STOTEN_LAST                          = 0.0_q   !< Root mean square error in predicted total energy for last added configuration
         REAL(q)                        :: STSIF                                = 0.0_q   !< Root mean square error in predicted stress.
         REAL(q)                        :: STSIF_LAST                           = 0.0_q   !< Root mean square error in predicted stress for last added configuration
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2(:)                                      !< Name of the model system in the training dataset.
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_INPUT(:)                                !< Name of the model system in the training dataset.
         REAL(q)                        :: TIFORAV(1:3)                         = 0.0_q   !< Average total forces.
         REAL(q)                        :: TOTENAV                              = 0.0_q   !< Average total energies.
         REAL(q)                        :: TSIFAV(1:3,1:3)                      = 0.0_q   !< Average total stress tensors.
         REAL(q)                        :: WTIFOR                               = 0.0_q   !< Weight for force.
         REAL(q)                        :: WTOTEN                               = 0.0_q   !< Weight for energy.
         REAL(q)                        :: WTSIF                                = 0.0_q   !< Weight for stress.
         CHARACTER(LEN=2),  ALLOCATABLE :: TYPE(:)                                        !< Type of ions
      END TYPE


      TYPE REDUCED_VARS
         ! Variables for spline grid and basis functions
         INTEGER                      :: BASIS_TYPE2_RED     = 0
         REAL(q), ALLOCATABLE         :: FNL2_RED(:,:,:,:,:)
         REAL(q), ALLOCATABLE         :: QGRID2_RED(:,:)
         REAL(q), ALLOCATABLE         :: RB2_RED(:,:,:)
         REAL(q), ALLOCATABLE         :: RG2_RED(:)
         REAL(q), ALLOCATABLE         :: RNORM2_RED(:,:)               !< Norm of Bessel function.
 
         ! Variables for CLM and DCLM
         REAL(q), ALLOCATABLE         :: DCLM2_RED(:,:,:,:)
         
         ! This entity contains variables needed for reduced descriptor.
         REAL(q), ALLOCATABLE         :: DPS_UGANG_RED(:,:,:,:) 
         INTEGER                      :: IVAR_OFFSET2_RED
         INTEGER, ALLOCATABLE         :: OFFSET_MAP2(:,:,:)
         INTEGER, ALLOCATABLE         :: L_LNRB2_RED(:)                !< List of angular momentum number used for reduced descriptor.
         LOGICAL, ALLOCATABLE         :: LFLAG_VAR_RED(:,:,:,:,:,:)
         INTEGER                      :: LMAX2_RED       = 0       !< Maximum angular momentum for color descriptor.
         INTEGER, ALLOCATABLE         :: LVAR_RED(:,:,:,:,:,:)         !< Label of angular descriptor.
         INTEGER                      :: MLNRB2_RED          = 0       !< (LMAX2+1)*MRB2
         INTEGER                      :: MMVAR2_RED      = 0       !< Maximum total number of variables for 2-point descriptor.
         INTEGER                      :: MMVAR2_ALLOCATE_RED = 0
         !INTEGER                      :: NR2             = 0       !< The number of radial mesh used.
         INTEGER                      :: MRB2_RED        = 0       !< Maximum number of radial basis sets.
         INTEGER                      :: MRB2_MAX_RED   = 0       !< Maximum distance to sum up for summed descriptor. 
         INTEGER                      :: NLNRB2_RED          = 0       !< SUM_L[NRB2(L)]
         INTEGER, ALLOCATABLE         :: NRB2_RED(:)                   !< The number of radial basis set.
         INTEGER, ALLOCATABLE         :: NRB_LNRB2_RED(:)              !< List of basis set number used for reduced descriptor.
         INTEGER, ALLOCATABLE         :: NNVAR2_RED(:)                 !< Total number of variables.
         INTEGER, ALLOCATABLE         :: OFFSET_MAP2_RED(:,:,:)
         REAL(q), ALLOCATABLE         :: PS_RED(:,:,:)                 !< Angular descriptor for basis set.
         REAL(q), ALLOCATABLE         :: PS_UGANG_RED(:,:,:)           !< Array for angular descriptor for basis set.
         REAL(q), POINTER, CONTIGUOUS :: PS_ALL_RED(:,:)     => NULL() !< Non-distributed angular descriptor.
         REAL(q), ALLOCATABLE         :: PS_NEW_RED(:,:,:,:)           !< Angular descriptor for basis set.
         REAL(q)                      :: RATIO_DESC          = 0.0_q
         REAL(q)                      :: RATIO_DESC_LOW_RES  = 0.0_q
         REAL(q)                      :: FACTOR_DESC         = 0.0_q
         REAL(q), ALLOCATABLE         :: RNORM_PS_RED(:,:)
         REAL(q), ALLOCATABLE         :: WVAR_RED(:,:,:,:,:,:)         !< Weight for angular descriptor.
      END TYPE 

!--------------------------------------------------------------------------------------
!
!> Handle for many-body variables of force field.
!
!--------------------------------------------------------------------------------------

      TYPE FFM_PAR ! Parameters used for many-body interaction fitting.
         TYPE (REDUCED_VARS)          :: RV
         REAL(q)                      :: AFILT2          = 0.0_q   !< Parameter for angular filtering
         INTEGER                      :: BASIS_TYPE2     = 0
         REAL(q), ALLOCATABLE         :: C00(:,:,:)                !< Expansion coefficients with l=m=0 used as pairwise descriptors for basis set.
         REAL(q), POINTER, CONTIGUOUS :: C00_ALL(:,:)    => NULL() !< Non-distributed C00 array.
         REAL(q), POINTER, CONTIGUOUS :: C00_W_ALL(:,:)  => NULL()
         REAL(q), ALLOCATABLE         :: C00_NEW(:,:,:,:)          !< New expansion coefficients with l=m=0 used as pairwise descriptors.
         REAL(q), ALLOCATABLE         :: CMAT(:,:)                 !< Covariance matrix that can be used for Bayesian error estimation.
         REAL(q), POINTER, CONTIGUOUS :: CMAT_ALL(:,:)   => NULL() !< Covariance matrix that can be used for Bayesian error estimation.
         REAL(q), ALLOCATABLE         :: CLM1(:,:,:)
         REAL(q), ALLOCATABLE         :: CLM2(:,:,:)
         REAL(q), ALLOCATABLE         :: DCLM1(:,:,:,:,:)
         REAL(q), ALLOCATABLE         :: DCLM2(:,:,:,:,:)
         REAL(q), ALLOCATABLE         :: DCLM1_HEAD(:,:,:,:,:)
         REAL(q), ALLOCATABLE         :: DCLM2_HEAD(:,:,:,:,:)
         INTEGER                      :: DESC_TYPE       = 0       !< Descriptor type (standard, linear-scaling with element types, etc ...)
         REAL(q)                      :: DHDRCOUPLE      = 0.0_q   !< Derivative of Hamiltonian with respect to the coupling parameter.
         REAL(q)                      :: DR1             = 0.0_q   !< Increment for radial mesh used for calculations of 1-point correlation (radial) descriptor.
         REAL(q)                      :: DR2             = 0.0_q   !< Increment for radial mesh used for calculations of 2-points correlation (angular) descriptor.
         REAL(q)                      :: DRSPL1          = 0.0_q   !< Increment for radial mesh used in spline-interpolation for 1-point correlation (radial) descriptor.
         REAL(q)                      :: DRSPL2          = 0.0_q   !< Increment for radial mesh used in spline-interpolation for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE         :: DMAT_FMAT(:)              !< Scaling vector for FMAT.
         REAL(q), ALLOCATABLE         :: DMAT_YMAT(:)              !< Scaling vector for FMAT.
         REAL(q), ALLOCATABLE         :: FMAT(:,:,:)               !< Design matrix.
         REAL(q), ALLOCATABLE         :: FMAT_NEW(:,:,:)           !< Design matrix for new configurations.
         REAL(q), ALLOCATABLE         :: FNL1(:,:,:,:,:)           !< Spline-interpolated radial functions used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE         :: FNL2(:,:,:,:,:)           !< Spline-interpolated radial functions used for calculating 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE         :: FWIN1(:,:)                !< Window function used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE         :: FWIN2(:,:)                !< Window function used for calculating 2-points correlation (angular) descriptor.
         INTEGER                      :: IAFILT2         = 0       !< Type of angular filtering.
         INTEGER                      :: IBROAD1         = 0       !< Type of broadening method for 1-point correlation (radial) descriptor.
         INTEGER                      :: IBROAD2         = 0       !< Type of broadening method for 2-points correlation (angular) descriptor.
         INTEGER, ALLOCATABLE         :: ICOUPLE(:)                !< The list of atoms where the coupling parameter is introduced for the calculatoin of the chemical potential.
         !> Cutoff function type (radial).
         !>
         !> If ICUT1=1, use Behler-Parinello cutoff function.
         !> If ICUT1=2, use Miwa's cutoff function.
         !> If ICUT1=3, used Jinnouchi's cutoff function.
         INTEGER                      :: ICUT1           = 0
         !> Cutoff function type (angular).
         !>
         !> If ICUT2=1, use Behler-Parinello cutoff function.
         !> If ICUT2=2, use Miwa's cutoff function.
         !> If ICUT2=3, use Jinnouchi's cutoff function.
         INTEGER                      :: ICUT2           = 0
         INTEGER, ALLOCATABLE         :: INDEX_MAP_TO_PS(:,:,:)    !< Mapping of the small number of indices to large number for PS
         INTEGER, ALLOCATABLE         :: INDEX_NUM2(:,:)           !< Actual number of descriptor indices for each C_nlm^iJ used in the calculation of dP/dC
         INTEGER, ALLOCATABLE         :: INDEX_PAIR_HELP1(:,:,:)   !< Index swapping help for derivatives of C00 with respect to coefficients
         INTEGER, ALLOCATABLE         :: INDEX_PAIR_HELP2(:,:,:)   !< Index swapping help for derivatives of PS with respect to coefficients
         INTEGER, ALLOCATABLE         :: INDEX_NUM2_B(:,:)           !< Actual number of descriptor indices for each C_nlm^iJ used in the calculation of dP/dC for linear-scaling descriptor
         INTEGER, ALLOCATABLE         :: INDEX_MAP_TO_PS_B(:,:,:)    !< Mapping of the small number of indices to large number for PS for linear-scaling descriptor
         INTEGER, ALLOCATABLE         :: INDEX_PAIR_HELP2_B(:,:,:,:) !< Index swapping help for derivatives of PS w.r.t. coeff. for linear-scaling descriptor
         !> If INVERSE_SOAP=1, diagonalization method is used.
         !> If INVERSE_SOAP=2, LU factorization is used.
         INTEGER                      :: INVERSE_SOAP    = 0
         INTEGER                      :: IREG            = 0       !< Type of method to optimize the regularization parameters
         !> Type of the way to scale the total energy for the many-body term.
         !>
         !> ISCALE_TOTEN=1: Scale to the total energies of isolated atoms.
         !> ISCALE_TOTEN=2: Scale to the average of the training energy data.
         INTEGER                      :: ISCALE_TOTEN    = 0
         INTEGER                      :: IVAR_OFFSET1    = 0
         INTEGER                      :: IVAR_OFFSET2    = 0
         !> Type of the similarity measure.
         !>
         !> ML_FF_LSUPERVEC_MB = .FALSE. : The similarity measure written in "Jinnouchi et al., Phys. Rev. B,
         !> "On-the-fly machine learning force field generation: Application to melting points" (2019).
         !> ML_FF_LSUPERVEC_MB = .TRUE. : A new similarity measure using super vector.
         LOGICAL                      :: LSUPERVEC       = .FALSE.
         !> Type of window function used for calculating 1-point correlation (radial) descriptor.
         !>
         !> IWINDOW1=1 : Hann window.
         !> IWINDOW1=2 : Hamming window.
         !> IWINDOW1=3 : Blackman window.
         !> IWINDOW1=4 : Nuttall window
         !> IWINDOW1=5 : Blackman-Nutall window
         !> IWINDOW1=6 : Blackman-Harris window
         INTEGER                      :: IWINDOW1        = 0
         !> Type of window function used for calculating 2-points correlation (angular) descriptor.
         !> IWINDOW2=1 : Hann window.
         !> IWINDOW2=2 : Hamming window.
         !> IWINDOW2=3 : Blackman window.
         !> IWINDOW2=4 : Nuttall window
         !> IWINDOW2=5 : Blackman-Nutall window
         !> IWINDOW2=6 : Blackman-Harris window
         INTEGER                      :: IWINDOW2        = 0
         LOGICAL                      :: LAFILT2         = .FALSE. !< If LAFILT2=.TRUE., angular filtering is executed.
         INTEGER, ALLOCATABLE         :: LB(:,:)                   !< Label of basis set.
         INTEGER, ALLOCATABLE         :: LBHEAD(:)                 !< List of begining of basis sets
         LOGICAL                      :: LCOUPLE         = .FALSE. !< If LCOUPLE=.TRUE., the coupling parameter is introduced for the calculation of the chemical potential.
         LOGICAL, ALLOCATABLE         :: LIONS_COUPLE(:)           !< For coupled atoms, LIONS_COUPLE is set to .TRUE.
         INTEGER, ALLOCATABLE         :: L_LNRB1(:)                !< List of angular momentum number used for calculating 1-point correlation (radial) descriptor.
         INTEGER, ALLOCATABLE         :: L_LNRB2(:)                !< List of angular momentum number used for calculating 2-points correlation (angular) descriptor.
         LOGICAL                      :: LEATOM          = .FALSE. !< If LEATOM=.TRUE., atomic energies are calculated.
         LOGICAL, ALLOCATABLE         :: LFLAG_VAR(:,:,:,:,:,:)
         LOGICAL, ALLOCATABLE         :: LFLAG_VAR_SIC(:,:,:,:,:)
         LOGICAL                      :: LHEAT           = .FALSE. !< It LHEAT=.TRUE., calculate heat flux
         INTEGER                      :: LMAX1           = 0       !< Maximum angular momentum number for 1-point correlation (radial) descriptor. This is set as zero.
         INTEGER                      :: LMAX2           = 0       !< Maximum angular momentum number for 2-points correlation (angular) descriptor.
         !> A logical variable that switches on or off the metric function used for calculating 1-point correlation (radial) descriptor.
         LOGICAL                      :: LMETRIC1        = .FALSE.
         !> A logical variable that switches on or off the metric function used for calculating 2-points correlation (angular) descriptor.
         LOGICAL                      :: LMETRIC2        = .FALSE.
         !> Maps the local reference configuration index (NB(TYPE)) and a given atom type, to which training structure it belongs to.
         !>
         !> First index: local reference configuration index (NB), depends also on a type
         !> Second index: Atom type
         !> Result: training structure number (usually ICONF)
         INTEGER, ALLOCATABLE         :: LNCONF(:,:)
         !> Maps the local reference configuration index (NB(TYPE)) and a given atom type, to which atom it belongs to.
         !>
         !> This is linked to LNCONF, so the atom should be one of the atoms
         !> of the corresponding LNCONF training structure.
         !> First index: local reference configuration index (NB), depends also on a type
         !> Second index: Atom type
         !> Result: atom number (usually INIONS) belonging to ICONF (from LNCONF(NB,TYPE))
         INTEGER, ALLOCATABLE         :: LNIONS(:,:)
         LOGICAL                      :: LNORM1          = .FALSE. !< If LNORM1=.TRUE., normalization is executed for 1-point correlation (radial) descriptor.
         LOGICAL                      :: LNORM2          = .FALSE. !< If LNORM2=.TRUE., normalization is executed for 2-points correlation (angular) descriptor.
         LOGICAL                      :: LSIC            = .FALSE. !< If LSIC=.TRUE., self-interaction is removed from the angular descriptor.
         LOGICAL                      :: LSPARSDES       = .FALSE. !< If LSPARSDES=.TRUE., descriptors are sparsified. This can be .TRUE., only when ISTART=1.
         !> List of training force data provided by the stored ab initio data.
         !>
         !> For a specific training set element it gives the map to it's position in the design matrix for example.
         !> First index: cartesian direction of force.
         !> Second index: respective atom for force.
         !> Third index: training configuration.
         INTEGER, ALLOCATABLE         :: LTIFOR(:,:,:)
         INTEGER, ALLOCATABLE         :: LTIFOR_NEW(:,:,:)         !< List of training force data provided by the new configurations.
         !> List of training energy data provided by the stored ab initio data.
         !> Similar as LTIFOR, read description of that.
         !> First dimension: training configuration.
         INTEGER, ALLOCATABLE         :: LTOTEN(:)
         INTEGER, ALLOCATABLE         :: LTOTEN_NEW(:)             !< List of training energy data provided by the new configurations.
         !> List of training stress tensor data provided by the stored ab initio data.
         !>
         !> Similar as LTIFOR, read description of that.
         !> First dimension: First tensor direction.
         !> Second dimension: Second tensor direction.
         !> Third index: training configuration.
         INTEGER, ALLOCATABLE         :: LTSIF(:,:,:)
         INTEGER, ALLOCATABLE         :: LTSIF_NEW(:,:,:)          !< List of training stress tensor data provided by the new configurations.
         INTEGER, ALLOCATABLE         :: LVAR(:,:,:,:,:,:)         !< Label of angular descriptor.
         INTEGER, ALLOCATABLE         :: LVAR_SIC(:,:,:,:,:)       !< Label of angular descriptor.
         LOGICAL                      :: LVARTRAN1       = .FALSE. !< If LVARTRAN1=.TRUE., the radial coordinate is transformed.
         LOGICAL                      :: LVARTRAN2       = .FALSE. !< If LVARTRAN1=.TRUE., the radial coordinate is transformed.
         !> If LWINDOW1=.TRUE., window function is used to smoothen the expanded atomic distribution in calculations of 1-point correlation (radial) descriptor.
         LOGICAL                      :: LWINDOW1        = .FALSE.
         !> If LWINDOW2=.TRUE., window function is used to smoothen the expanded atomic distribution in calculations of 2-points correlation (angular) descriptor.
         LOGICAL                      :: LWINDOW2        = .FALSE.
         INTEGER                      :: MB              = 0       !< Maximum number of basis sets provided by the stored ab initio data.
         INTEGER                      :: MB_INPUT        = 0
         INTEGER                      :: MB_ALLOCATE     = 0
         INTEGER                      :: MB_TARGET       = 0
         INTEGER                      :: MB_NEW          = 0       !< Maximum number of basis sets provided by the new configurations.
         INTEGER                      :: MB_MIN          = 0       !< Minimum number of local reference configurations required (after sparsification) to trigger training.
         INTEGER                      :: MBTOT           = 0       !< Maximum total number of basis sets provided by the stored ab initio data.
         INTEGER                      :: MBTOT_ALLOCATE  = 0
         INTEGER                      :: MBTOT_NEW       = 0       !< Maximum total number of basis sets provided by the new configurations.
         INTEGER                      :: MDATA           = 0       !< Maximum number of training data provided by the stored ab initio data.
         INTEGER                      :: MDATA_ALLOCATE
         INTEGER                      :: MDATA_NEW       = 0       !< Maximum number of training data provided by the new configurations.
         INTEGER                      :: MDATA_NEW_ALLOCATE
         INTEGER                      :: MLNRB1          = 0       !< (LMAX1+1)*MRB1
         INTEGER                      :: MLNRB2          = 0       !< (LMAX2+1)*MRB2
         INTEGER                      :: MRB1            = 0       !< Maximum number of radial basis sets used in calculations of 1-point correlation (radial) descriptor.
         INTEGER                      :: MRB2            = 0       !< Maximum number of radial basis sets used in calculations of 2-points correlation (angular) descriptor.
         INTEGER                      :: MRB2_MAX        = 0       !< Maximum distance to sum up for summed descriptor. 
         INTEGER                      :: MSPL1           = 0       !< Maximum number of radial grids used in spline-interpolation for calculations of 1-point correlation (radial) descriptor.
         INTEGER                      :: MSPL2           = 0       !< Maximum number of radial grids used in spline-interpolation for calculations of 2-points correlation (angular) descriptor.
         INTEGER                      :: MMVAR1          = 0       !< Maximum total number of variables for 1-point descriptor.
         INTEGER                      :: MMVAR1_ALLOCATE
         INTEGER                      :: MMVAR2          = 0       !< Maximum total number of variables for 2-point descriptor.
         INTEGER                      :: MMVAR2_ALLOCATE
         INTEGER                      :: MMVAR2_SIC      = 0       !< Maximum total number of variables for self-interaction terms in 2-point (angular) descriptor.
         INTEGER                      :: MMVAR2_SIC_ALLOCATE
         INTEGER                      :: MCOL_DCLM1_MAX  = 0       !< Maximum number of neighbors for derivatives of radial descriptor with respect to r
         INTEGER                      :: MCOL_DCLM2_MAX  = 0       !< Maximum number of neighbors for derivatives of angular descriptor with respect to r
         INTEGER                      :: MCOL_DCLM2_MAX_RED  = 0       !< Maximum number of neighbors for derivatives of angular descriptor with respect to r
         INTEGER                      :: NATOM_COUPLED   = 0       !< The number of atom where the coupling parameter is introduced.
         INTEGER, ALLOCATABLE         :: NB(:)                     !< Number of basis set for each species provided by the stored ab initio data.
         INTEGER, ALLOCATABLE         :: NB_NEW(:)                 !< The number of basis set for each species provided by the new configurations.
         INTEGER                      :: NBTOT           = 0       !< The total number of basis sets provided by the stored ab initio data.
         INTEGER                      :: NBTOT_NEW       = 0       !< The total number of basis sets provided by the new configurations.
         INTEGER                      :: NDATA           = 0       !< The number of training data provided by the stored ab initio data.
         INTEGER                      :: NDATA_NEW       = 0       !< The number of training data provided by the new configurations.
         INTEGER                      :: NHYP1           = 0       !< Hyper-parameter for SOAP 1-point kernel.
         INTEGER                      :: NHYP2           = 0       !< Hyper-parameter for SOAP 2-point kernel.
         INTEGER                      :: NLNRB1          = 0       !< SUM_L[NRB1(L)]
         INTEGER                      :: NLNRB2          = 0       !< SUM_L[NRB2(L)]
         INTEGER                      :: NMETRIC1        = 0       !< Polynomical parameter for metric function used for calculating 1-point correlation (radial) descriptor.
         INTEGER                      :: NMETRIC2        = 0       !< Polynomical parameter for metric function used for calculating 2-points correlation (angular) descriptor.
         INTEGER                      :: NR1             = 0       !< The number of radial mesh used for calculating 1-point correlation (radial) descriptor.
         INTEGER                      :: NR2             = 0       !< The number of radial mesh used for calculating 2-points correlation (angular) descriptor.
         INTEGER                      :: NRANK_SPARSDES  = 0       !< Rank of the matrix used in the sparsification of the descriptors.
         INTEGER, ALLOCATABLE         :: NRB1(:)                   !< The number of radial basis set used for calculating 1-point correlation (radial) descriptor.
         INTEGER, ALLOCATABLE         :: NRB2(:)                   !< The number of radial basis set used for calculating 2-points correlation (angular) descriptor.
         INTEGER, ALLOCATABLE         :: NRB_LNRB1(:)              !< List of basis set number used for calculating 1-point correlation (radial) descriptor.
         INTEGER, ALLOCATABLE         :: NRB_LNRB2(:)              !< List of basis set number used for calculating 2-points correlation (angular) descriptor.
         INTEGER                      :: NSPL1           = 0       !< The number of radial grids used for spline interpolation used for calculating 1-point correlation (radial) descriptor.
         INTEGER                      :: NSPL2           = 0       !< The number of radial grids used for spline interpolation used for calculating 2-points correlation (angular) descriptor.
         INTEGER                      :: NTIFOR          = 0       !< The total number of training force data.
         INTEGER                      :: NTOTEN          = 0       !< The total number of training energy data.
         INTEGER                      :: NTSIF           = 0       !< The total number of training stress-tensor data.
         INTEGER                      :: NNVAR1          = 0       !< Total number of variables for 1-point descriptor.
         INTEGER                      :: NNVAR1_ALLOCATE
         INTEGER, ALLOCATABLE         :: NNVAR2(:)                 !< Total number of variables for 2-point descriptor.
         INTEGER, ALLOCATABLE         :: NNVAR2_SIC(:)             !< Total number of variables for self-interaction 2-point (angular) descriptor.
         INTEGER                      :: NVARTRAN1       = 0       !< Polynomial parameter used for radial variable transform.
         INTEGER                      :: NVARTRAN2       = 0       !< Polynomial parameter used for radial variable transform.
         INTEGER, ALLOCATABLE         :: OFFSET_MAP1(:,:,:)
         INTEGER, ALLOCATABLE         :: OFFSET_MAP2(:,:,:)
         REAL(q), ALLOCATABLE         :: PREFACTOR_HELP1(:,:,:)    !< Prefactor for C_nlm related variables for radial descriptor
         REAL(q), ALLOCATABLE         :: PREFACTOR_HELP2(:,:,:)    !< Prefactor for C_nlm related variables for angular descriptor
         REAL(q), ALLOCATABLE         :: PREFACTOR_HELP2_B(:,:,:)  !< Prefactor for C_nlm related variables for angular descriptor
         REAL(q), ALLOCATABLE         :: PS(:,:,:)                 !< Angular descriptor for basis set.
         REAL(q), POINTER, CONTIGUOUS :: PS_ALL(:,:)     => NULL() !< Non-distributed angular descriptor.
         REAL(q), POINTER, CONTIGUOUS :: PS_W_ALL(:,:)   => NULL()
         REAL(q), ALLOCATABLE         :: PS_NEW(:,:,:,:)           !< Angular descriptor for basis set.
         REAL(q), ALLOCATABLE         :: QGRID1(:,:)               !< Q-grid for spherical Bessel functions used as radial basis sets for 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE         :: QGRID2(:,:)               !< Q-grid for spherical Bessel functions used as radial basis sets for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE         :: RB1(:,:,:)                !< Radial basis set for 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE         :: RB2(:,:,:)                !< Radial basis set for 2-points correlation (angular) descriptor.
         REAL(q)                      :: RCOUPLE         = 0.0_q   !< Coupliing parameter for the calculation of the chemical potential.
         REAL(q), ALLOCATABLE         :: RCOUPLE_IONS(:)           !< Coupling parameters for each ions.
         !> Radial grid for 1-point correlation (radial) descriptor.
         !>
         !> When metric is introduced, this grid corresponds to the equispaced grid over R^N.
         !> Otherwise, this grid corresponds to the equispaced over R.
         REAL(q), ALLOCATABLE         :: RG1(:)
         !> Radial grid for 2-points correlation (angular) descriptor.
         !>
         !> When metric is introduced, this grid corresponds to the equispaced grid over R^N.
         !> Otherwise, this grid corresponds to the equispaced over R.
         REAL(q), ALLOCATABLE         :: RG2(:)
         REAL(q)                      :: RCUT1           = 0.0_q   !< Cutoff radius for 1-point correlation (radial) descriptor.
         REAL(q)                      :: RCUT2           = 0.0_q   !< Cutoff radius 2-points correlation (angular) descriptor.
         REAL(q)                      :: RDES_SPARSDES   = 0.0_q   !< Ratio of the number of selecting descriptors to all.
         REAL(q)                      :: RMETRIC1        = 0.0_q   !< A parameter for metric function for 1-point correlation (radial) descriptor.
         REAL(q)                      :: RMETRIC2        = 0.0_q   !< A parameter for metric function for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE         :: RNORM1(:,:)               !< Norm of Bessel function used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE         :: RNORM2(:,:)               !< Norm of Bessel function used for calculating 2-points correlation (angular) descriptor.
         LOGICAL                      :: SAVECMAT        = .FALSE. !< Decides whether CMAT is recalculated for each evidence step or not
         REAL(q)                      :: SIG(1:10)       = 0.0_q   !< Variations in first principles data.
         REAL(q), ALLOCATABLE         :: SIG_SYS(:,:)              !< Variations in first principles data for individual system.
         REAL(q)                      :: SIGV ,SIGW      = 0.0_q   !< Pricision parameters determined by the evidence approximation.
         REAL(q)                      :: SIGV0,SIGW0     = 0.0_q   !< (Initial) pricision parameters determined by the evidence approximation.
         REAL(q)                      :: SION1           = 0.0_q   !< A parameter that specifies the broadening width for 1-point correlation (radial) descriptor.
         REAL(q)                      :: SION2           = 0.0_q   !< A parameter that specifies the broadening width for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE         :: SOAP(:,:,:)               !< SOAP matrix used for sparsification.
         REAL(q), ALLOCATABLE         :: SINV(:,:,:)               ! Inverse matrix of SOAP.
         REAL(q), POINTER, CONTIGUOUS :: SINV_ALL(:,:,:)
         REAL(q)                      :: TIFORAV(1:3)    = 0.0_q   !< Average total force given from many-body interactions.
         REAL(q)                      :: TOTENAV         = 0.0_q   !< Average total energy given from many-body interactions.
         REAL(q)                      :: TSIFAV(1:3,1:3) = 0.0_q   !< Average total stress given from many-body interactions.
         REAL(q)                      :: W1              = 0.0_q   !< Weight for 1-point descriptor in SOAP descriptors.
         REAL(q)                      :: W2              = 0.0_q   !< Weight for 2-point descriptor in SOAP descriptors.
         REAL(q), ALLOCATABLE         :: WION1(:)                  !< Broadening width used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE         :: WION2(:)                  !< Broadening width used for calculating 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE         :: WMAT(:,:,:)               !< Regression coefficients.
         REAL(q), ALLOCATABLE         :: WMAT_ALL(:)               !< Regression coefficients owned by all processes.
         LOGICAL                      :: WRITECMAT       = .FALSE. !< Decides what is written at end of file (CMAT or SINV)
         REAL(q), ALLOCATABLE         :: WVAR(:,:,:,:,:,:)         !< Weight for angular descriptor.
         REAL(q), ALLOCATABLE         :: YMAT(:,:)                 !< Data vector
         ! Model potentials:
         LOGICAL                      :: LEMPPOT         = .FALSE. !< Enables model potentials for thermodynamic integration
      END TYPE


!--------------------------------------------------------------------------------------
!
!> Important variables used for the logfile.
!
!--------------------------------------------------------------------------------------

      TYPE LOG_PAR
         INTEGER           :: IU             = 0
         CHARACTER(LEN=14) :: ML_LOGFILE     = ""
         INTEGER           :: NWRITE         = 0
         LOGICAL           :: SKIP_LOG_LINES = .FALSE.
      END TYPE


!--------------------------------------------------------------------------------------
!
!> Variables needed for gathering information about memory usage.
!
!--------------------------------------------------------------------------------------

      TYPE MEMORY
         REAL(q) :: RMEM(1:11)    = 0.0_q
         REAL(q) :: RMEM_CMAT     = 0.0_q !< Memory used for storing covariance matrix of many-body interactions.
         REAL(q) :: RMEM_CMAT_ALL = 0.0_q !< Memory used for storing covariance matrix of many-body interactions.
         REAL(q) :: RMEM_FMAT     = 0.0_q !< Memory used for storing design matrix of many-body interactions.
         REAL(q) :: RMEM_PS       = 0.0_q !< Memory used for storing angular descriptor.
         REAL(q) :: RMEM_PS_ALL   = 0.0_q !< Memory used for storing non-distributed angular descriptor.
         REAL(q) :: RMEM_RST      = 0.0_q !< Memory necessary for restructuring of the arrays in sparsification process
         REAL(q) :: RMEM_SINV     = 0.0_q !< Memory used for storing inverse SOAP matrix necessary for error estimation
         REAL(q) :: RMEM_SINV_ALL = 0.0_q !< Memory used for storing inverse SOAP matrix necessary for error estimation
         REAL(q) :: RMEM_SOAP     = 0.0_q !< Memory used for storing inverse SOAP matrix necessary for error estimation
         REAL(q) :: RMEM_SOAP_ALL = 0.0_q !< Memory used for storing inverse SOAP matrix necessary for error estimation
         REAL(q) :: RMEM_SVD      = 0.0_q !< Memory used for SVD decomposition
         REAL(q) :: RMEM_WMAT_ALL = 0.0_q !< Memory used for storing non-distributed regression coefficients.
      END TYPE MEMORY


!--------------------------------------------------------------------------------------
!
! A lot of parallel stuff.
!
!--------------------------------------------------------------------------------------

      TYPE SCALAPACK_GRID
         INTEGER :: ICTXT  = 0 !< context number
         INTEGER :: NDIM   = 0 !< dimension of scaLAPACK grid
         INTEGER :: NP_ROW = 0 !< number of rows in processor grid
         INTEGER :: NP_COL = 0 !< number of columns in processor grid
         INTEGER :: ME_ROW = 0 !< row coordinate in processor grid
         INTEGER :: ME_COL = 0 !< column coordinate in processor grid
         INTEGER :: NPROCS = 0 !< number of processors in context
      END TYPE

      TYPE PARALLEL_DESCRIPTOR
         INTEGER :: DESC(1:9) = 0 !< descriptor number
         INTEGER :: M         = 0 !< number of rows in global array
         INTEGER :: N         = 0 !< number of columns in global array
         INTEGER :: MB        = 0 !< blocking factor for rows
         INTEGER :: NB        = 0 !< blocking factor for columns
         INTEGER :: IRSRC     = 0 !< process row of first row of matrix distribution
         INTEGER :: ICSRC     = 0 !< process column of first column of matrix distribution
         INTEGER :: LLD       = 0 !< Leading dimension of local array(MAX(1, LOCr(M))
         INTEGER :: INFO      = 0 !< info
         INTEGER :: LOC_M     = 0 !< local size of M specific to one process
         INTEGER :: LOC_N     = 0 !< local size of N specific to one process
      END TYPE

      TYPE SCALA_WORK
         INTEGER              :: INFO     = 0
         INTEGER, ALLOCATABLE :: IPIV(:)
         INTEGER, ALLOCATABLE :: IWORK(:)
         INTEGER              :: LOC_M    = 0
         INTEGER              :: LWORK    = 0
         INTEGER              :: LIWORK   = 0
         REAL(q), ALLOCATABLE :: WORK(:)
      END TYPE

#ifdef use_shmem
      TYPE M_shmem_desc
#ifndef sysv
         TYPE (MPI_Win)           :: shmem_win
#else
         INTEGER(KIND=c_int)      :: shmem_win = 0
#endif
         TYPE (c_ptr)             :: shmem_ptr
#if defined(MPI) || defined(MPI_CHAIN)
         TYPE (ML_MPI_PAR), POINTER  :: comm => NULL()
#endif
         INTEGER                  :: dim1      = 0
         INTEGER                  :: dim2      = 0
         INTEGER                  :: dim3      = 0
      END TYPE M_shmem_desc

      TYPE M_shmem
         TYPE (M_shmem_desc) :: desc
         COMPLEX(q), POINTER :: fptr_c(:) => NULL()
         REAL(q), POINTER    :: fptr_r(:) => NULL()
      END TYPE M_shmem

#endif


      TYPE PARALLEL_SUPER
#if defined(MPI) || defined(MPI_CHAIN)
         TYPE(ML_MPI_PAR)          :: COMM_WORLD
#ifdef use_shmem
         TYPE(ML_MPI_PAR)          :: COMM_ZERO
         TYPE(ML_MPI_PAR)          :: COMM_SHMEM
#endif
         TYPE(ML_MPI_PAR)          :: COMM_NB
         TYPE(ML_MPI_PAR)          :: COMM_NIONS
         TYPE(ML_MPI_PAR)          :: COMM_NEIB
#endif
         TYPE(SCALAPACK_GRID)      :: CONTEXT_WORLD
         TYPE(SCALAPACK_GRID)      :: CONTEXT_SINGLE
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_A_SPARSDES_CUR
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_BE_ENERGY
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_BE_FORCE
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_BE_STRESS
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_C00
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_C00_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_C00_NEW_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_CMAT
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_CMAT0
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_CMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DC00
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DCLM1
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DCLM2
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DCLM2_RED
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DPS
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DSOAP1
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_DSOAP2
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_EATOM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_ENERGY_HELP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_FORCE_HELP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_FMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_FMAT_NEW_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_FMAT_FFM_TRANSPOSED
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_GMAT
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_GENERAL_SINGLE
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_HEAT1
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_HEAT2
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_PS
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_PS_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_PS_FFM_RED
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_PS_NEW_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_RBS
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SOAP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SOAP_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SOAP_NEW
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SOAP_NEW_TRANSPOSED
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SOAP_HELP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SOAP_TRANSPOSED
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_STRESS_HELP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_V_SPARSDES_CUR
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_ENERGY
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_FORCE1
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_FORCE2
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_HEAT1
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_HEAT2
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_STRESS1
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_VEC_LOC_STRESS2
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_WMAT
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_WMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_YMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_ZMAT
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_ZMAT_HEAT
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_FMAT_HELP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_FMAT_TRANS_HELP
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SING_VT
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SING_U
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SING_SIGMA
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SING_VSIG
         TYPE(PARALLEL_DESCRIPTOR) :: DESC_SING_B1_SVD
         TYPE(SCALA_WORK)          :: SCAWORK
         REAL(q),ALLOCATABLE       :: GENERAL_SINGLE(:,:) !< Array mainly to reduce distributed
      END TYPE PARALLEL_SUPER

!--------------------------------------------------------------------------------------
!
!> All important shared memory variables used in ml_ff_ff.F.
!
!--------------------------------------------------------------------------------------

      TYPE ML_SHMEM
         REAL(q), POINTER, CONTIGUOUS :: CLM1(:,:)                                  => NULL()
         REAL(q), POINTER, CONTIGUOUS :: CLM2(:,:)                                  => NULL()
         REAL(q), POINTER, CONTIGUOUS :: CLM2_RED(:,:)                              => NULL()
         REAL(q), POINTER, CONTIGUOUS :: CLM_COUPLE_HELP1_MANYBODY(:,:)             => NULL()
         REAL(q), POINTER, CONTIGUOUS :: CLM_COUPLE_HELP2_MANYBODY(:,:)             => NULL()
         REAL(q), POINTER, CONTIGUOUS :: PS_SIC(:,:)                                => NULL()
         REAL(q), POINTER, CONTIGUOUS :: PS_SIC_COUPLE(:,:)                         => NULL()
         REAL(q), POINTER, CONTIGUOUS :: SOAP_ALL(:,:)                              => NULL()
         REAL(q), POINTER, CONTIGUOUS :: SOAP_C00(:,:)                              => NULL()
         REAL(q), POINTER, CONTIGUOUS :: SOAP_PS(:,:)                               => NULL()
         REAL(q), POINTER, CONTIGUOUS :: SOAP_PS_RED(:,:)                           => NULL()
         INTEGER                      :: MCOL_CLM1_GLOBAL                           = 0
         INTEGER                      :: MCOL_CLM2_GLOBAL                           = 0
         INTEGER                      :: MCOL_CLM2_GLOBAL_RED                           = 0
         INTEGER                      :: MCOL_CLM_COUPLE_HELP1_MANYBODY             = 0
         INTEGER                      :: MCOL_CLM_COUPLE_HELP2_MANYBODY             = 0
         INTEGER                      :: MCOL_PS_SIC_GLOBAL                         = 0
         INTEGER                      :: MCOL_PS_SIC_COUPLE_GLOBAL                  = 0
         INTEGER                      :: MCOL_SOAP_ALL_GLOBAL                       = 0
         INTEGER                      :: MCOL_SOAP_C00_GLOBAL                       = 0
         INTEGER                      :: MCOL_SOAP_PS_GLOBAL                        = 0
         INTEGER                      :: MROW_CLM1_GLOBAL                           = 0
         INTEGER                      :: MROW_CLM2_GLOBAL                           = 0
         INTEGER                      :: MROW_CLM2_GLOBAL_RED                           = 0
         INTEGER                      :: MROW_CLM_COUPLE_HELP1_MANYBODY             = 0
         INTEGER                      :: MROW_CLM_COUPLE_HELP2_MANYBODY             = 0
         INTEGER                      :: MROW_PS_SIC_GLOBAL                         = 0
         INTEGER                      :: MROW_PS_SIC_COUPLE_GLOBAL                  = 0
         INTEGER                      :: MROW_SOAP_ALL_GLOBAL                       = 0
         INTEGER                      :: MROW_SOAP_C00_GLOBAL                       = 0
         INTEGER                      :: MROW_SOAP_PS_GLOBAL                        = 0
#ifdef use_shmem
         TYPE (M_shmem)               :: MLFF_SHM_2D_C00_ALL
         TYPE (M_shmem)               :: MLFF_SHM_2D_C00_W_ALL
         TYPE (M_shmem)               :: MLFF_SHM_2D_PS_ALL
         TYPE (M_shmem)               :: MLFF_SHM_2D_PS_ALL_RED
         TYPE (M_shmem)               :: MLFF_SHM_2D_PS_W_ALL
         TYPE (M_shmem)               :: MLFF_SHM_2D_CLM1
         TYPE (M_shmem)               :: MLFF_SHM_2D_CLM2
         TYPE (M_shmem)               :: MLFF_SHM_2D_CLM2_RED
         TYPE (M_shmem)               :: MLFF_SHM_2D_CLM_COUPLE_HELP1_MANYBODY
         TYPE (M_shmem)               :: MLFF_SHM_2D_CLM_COUPLE_HELP2_MANYBODY
         TYPE (M_shmem)               :: MLFF_SHM_2D_PS_SIC
         TYPE (M_shmem)               :: MLFF_SHM_2D_PS_SIC_COUPLE
         TYPE (M_shmem)               :: MLFF_SHM_2D_SOAP_ALL
         TYPE (M_shmem)               :: MLFF_SHM_3D_SINV_ALL
         TYPE (M_shmem)               :: MLFF_SHM_2D_SOAP_C00
         TYPE (M_shmem)               :: MLFF_SHM_2D_SOAP_PS
         TYPE (M_shmem)               :: MLFF_SHM_2D_SOAP_PS_RED
         TYPE (M_shmem)               :: MLFF_SHM_2D_CMAT_ALL
         LOGICAL                      :: MLFF_SHMEM_ALLOC_C00_ALL                   = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_C00_W_ALL                 = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_PS_ALL                    = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_PS_ALL_RED                = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_PS_W_ALL                  = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_CLM1                      = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_CLM2                      = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_CLM2_RED                      = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_CLM_COUPLE_HELP1_MANYBODY = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_CLM_COUPLE_HELP2_MANYBODY = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_PS_SIC                    = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_PS_SIC_COUPLE             = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_SOAP_ALL                  = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_SINV_ALL                  = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_SOAP_C00                  = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_SOAP_PS                   = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_SOAP_PS_RED               = .FALSE.
         LOGICAL                      :: MLFF_SHMEM_ALLOC_CMAT_ALL                  = .FALSE.
#else
         REAL(q),POINTER              :: C00_ALL_NOSHMEM(:,:)                       => NULL()
         REAL(q),POINTER              :: C00_W_ALL_NOSHMEM(:,:)                     => NULL()
         REAL(q),POINTER              :: PS_ALL_NOSHMEM(:,:)                        => NULL()
         REAL(q),POINTER              :: PS_ALL_NOSHMEM_RED(:,:)                    => NULL()
         REAL(q),POINTER              :: PS_W_ALL_NOSHMEM(:,:)                      => NULL()
         REAL(q),POINTER              :: CLM1_NOSHMEM(:,:)                          => NULL()
         REAL(q),POINTER              :: CLM2_NOSHMEM(:,:)                          => NULL()
         REAL(q),POINTER              :: CLM2_NOSHMEM_RED(:,:)                          => NULL()
         REAL(q),POINTER              :: CLM_COUPLE_HELP1_MANYBODY_NOSHMEM(:,:)     => NULL()
         REAL(q),POINTER              :: CLM_COUPLE_HELP2_MANYBODY_NOSHMEM(:,:)     => NULL()
         REAL(q),POINTER              :: PS_SIC_NOSHMEM(:,:)                        => NULL()
         REAL(q),POINTER              :: PS_SIC_COUPLE_NOSHMEM(:,:)                 => NULL()
         REAL(q),POINTER              :: SOAP_ALL_NOSHMEM(:,:)                      => NULL()
         REAL(q),POINTER              :: SINV_ALL_NOSHMEM(:,:,:)                    => NULL()
         REAL(q),POINTER              :: SOAP_C00_NOSHMEM(:,:)                      => NULL()
         REAL(q),POINTER              :: SOAP_PS_NOSHMEM(:,:)                       => NULL()
         REAL(q),POINTER              :: SOAP_PS_NOSHMEM_RED(:,:)                   => NULL()
         REAL(q),POINTER              :: CMAT_ALL_NOSHMEM(:,:)                      => NULL()
#endif
         INTEGER                      :: ISTATUS_SHMEM                              = 0
         INTEGER                      :: MROW_CLM1_GLOBAL_ALLOCATE                  = 0
         INTEGER                      :: MCOL_CLM1_GLOBAL_ALLOCATE                  = 0
         INTEGER                      :: MROW_CLM_COUPLE_HELP1_MANYBODY_ALLOCATE    = 0
         INTEGER                      :: MCOL_CLM_COUPLE_HELP1_MANYBODY_ALLOCATE    = 0
         INTEGER                      :: MROW_CLM2_GLOBAL_ALLOCATE                  = 0
         INTEGER                      :: MROW_CLM2_GLOBAL_ALLOCATE_RED                  = 0
         INTEGER                      :: MCOL_CLM2_GLOBAL_ALLOCATE                  = 0
         INTEGER                      :: MCOL_CLM2_GLOBAL_ALLOCATE_RED                  = 0
         INTEGER                      :: MROW_CLM_COUPLE_HELP2_MANYBODY_ALLOCATE    = 0
         INTEGER                      :: MCOL_CLM_COUPLE_HELP2_MANYBODY_ALLOCATE    = 0
         INTEGER                      :: MROW_PS_SIC_GLOBAL_ALLOCATE                = 0
         INTEGER                      :: MCOL_PS_SIC_GLOBAL_ALLOCATE                = 0
         INTEGER                      :: MROW_PS_SIC_COUPLE_GLOBAL_ALLOCATE         = 0
         INTEGER                      :: MCOL_PS_SIC_COUPLE_GLOBAL_ALLOCATE         = 0
         INTEGER                      :: MROW_SOAP_ALL_GLOBAL_ALLOCATE              = 0
         INTEGER                      :: MCOL_SOAP_ALL_GLOBAL_ALLOCATE              = 0
         INTEGER                      :: MROW_SOAP_C00_GLOBAL_ALLOCATE              = 0
         INTEGER                      :: MCOL_SOAP_C00_GLOBAL_ALLOCATE              = 0
         INTEGER                      :: MROW_SOAP_PS_GLOBAL_ALLOCATE               = 0
         INTEGER                      :: MCOL_SOAP_PS_GLOBAL_ALLOCATE               = 0
      END TYPE ML_SHMEM

!--------------------------------------------------------------------------------------
!
!> Helping arrays neccessary for LAPACK operations.
!
!--------------------------------------------------------------------------------------
      TYPE LAPACK_WORK
         INTEGER              :: INFO    = 0
         INTEGER, ALLOCATABLE :: IPIV(:)
         INTEGER              :: LWORK   = 0
         REAL(q), ALLOCATABLE :: WORK(:)
      END TYPE

!--------------------------------------------------------------------------------------
!
!> Variables storing timing information.
!
!--------------------------------------------------------------------------------------
      TYPE TIME
         REAL(q) :: CPUTIM(1:12) = 0.0_q
         REAL(q) :: VPUTIM(1:12) = 0.0_q
      END TYPE TIME

!--------------------------------------------------------------------------------------
!
!> Variables storing timing information.
!
!--------------------------------------------------------------------------------------
      
      TYPE ASA_STUFF
         INTEGER :: NCG=13000 ! LMAX= 6 (f electrons)
!   (LMAX=8: NCG=45000, LMAX=10: NCG=123000, LMAX=14: NCG=588000, LMAX=20: NCG=3191000)

         REAL(q),ALLOCATABLE :: YLM3(:)  ! table which contains the intregral of three
                                 ! real spherical harmonics
         REAL(q),ALLOCATABLE :: YLM3I(:) ! inverse table
         INTEGER,ALLOCATABLE :: JL(:)    ! index L for each element in the array YLM3
         INTEGER,ALLOCATABLE :: JS(:)    ! compound index L,M for each element in the array YLM3
                                 ! JS =  L*(L+1)+M+1 (where M=-L,...,L)
         INTEGER,ALLOCATABLE :: INDCG(:) ! index into array YLM3 which gives the starting
                                 ! position for one l,m ,lp,mp  quadruplet
         INTEGER :: LMAXCG=-1    ! maximum l
         INTEGER,ALLOCATABLE  :: YLM3LOOKUP_TABLE(:,:) !
      END TYPE ASA_STUFF

#if defined(MPI) || defined(MPI_CHAIN)
!--------------------------------------------------------------------------------------
!
!> Variables storing timing information.
!
!--------------------------------------------------------------------------------------
      TYPE EMP_POT_TYPE
         TYPE (VERLET_CELL_LIST_TYPE) :: NEAREST_NEIGHBOR
         INTEGER, ALLOCATABLE :: LNIONS_EST_MODEL(:, :, :)
         INTEGER, ALLOCATABLE :: NNEIB_EST_MODEL(:, :)
         REAL(q), ALLOCATABLE :: XYZ_EST_MODEL(:, :, :, :)
         REAL(q), ALLOCATABLE :: XYZ_EST_MODEL_NORM(:, :, :, :)
         REAL(q), ALLOCATABLE :: NORMALIZATION_XYZ_MODEL(:, :, :)
         INTEGER              :: MNEIB_EST_MODEL
         REAL(q)              :: EMPPOT_RCUT = 0.0_q
         ! Parameters for repulsive exponential function
         REAL(q)              :: SRPOT_S0    = 2.0_q !< Parameter for repulsive exponential function
         INTEGER              :: SRPOT_N0    = 10    !< Parameter for repulsive exponential function
         REAL(q)              :: SRPOT_B0    = 1.0_q !< Parameter for repulsive exponential function
         ! Parameters for Morse potential
         INTEGER              :: MOPOT_NM    = 0     !< Parameter for Morse potential
         INTEGER, ALLOCATABLE :: IM_COUPLE(:)        !< Parameter for Morse potential
         INTEGER, ALLOCATABLE :: JM_COUPLE(:)        !< Parameter for Morse potential
         REAL(q), ALLOCATABLE :: MOPOT_DM(:)         !< Parameter for Morse potential
         REAL(q), ALLOCATABLE :: MOPOT_RM(:)         !< Parameter for Morse potential
         REAL(q), ALLOCATABLE :: MOPOT_RKM(:)        !< Parameter for Morse potential
         INTEGER, ALLOCATABLE :: MOPOT_IJM(:)        !< Parameter for Morse potential
      END TYPE EMP_POT_TYPE
#endif


!--------------------------------------------------------------------------------------
!
!> This is the handle super type used from main.
!
!--------------------------------------------------------------------------------------
      TYPE ML_SUPER_TYPE
         TYPE(ABINITIO)              :: AB
         TYPE(ABINITIO)              :: AB_HELP
         TYPE(ABINITIONEW)           :: ABN
         TYPE(ML_FF_INPUT)           :: ML_INPUT_HANDLE
         TYPE(FF_PAR)                :: FF
         TYPE(FFM_PAR)               :: FFM
         TYPE(LOG_PAR)               :: LOGF
         TYPE(MEMORY)                :: MEM0
         TYPE(MEMORY)                :: MEM1
         TYPE(ML_SHMEM)              :: MLSHM
         TYPE(PARALLEL_SUPER)        :: PAR_SUP_HANDLE
         TYPE(LAPACK_WORK)           :: LAPACKWORK
         TYPE(TIME)                  :: TIM
         INTEGER                     :: INSTANCE                 = 0
         !> Copy of IO%IU0 from main.F, for consistent screen output from inside ML routines.
         INTEGER                     :: IO_IU0                   = -1
         !> Copy of TIU0 from main.F, for consistent screen output from inside ML routines.
         INTEGER                     :: TIU0                     = -1
         LOGICAL                     :: SKIP_ML_TO_VASP          = .FALSE.
         LOGICAL                     :: SKIP_VASP_TO_ML          = .FALSE.
         TYPE(ASA_STUFF)             :: ASA
#if defined(MPI) || defined(MPI_CHAIN)
         TYPE(EMP_POT_TYPE)          :: EMP
#endif
         !> class object for nearest neighbor list computation
#if defined(MPI) || defined(MPI_CHAIN)
         TYPE(TAGLIST)               :: TAG_LIST
         TYPE(VERLET_CELL_LIST_TYPE) :: NEAREST_NEIGHBOR_RADIAL
         TYPE(VERLET_CELL_LIST_TYPE) :: NEAREST_NEIGHBOR_ANGULAR
#endif
      END TYPE

! Global variables
      TYPE (ML_SUPER_TYPE), ALLOCATABLE :: ML_SUPER_HANDLE_MAIN(:)
#if defined(MPI) || defined(MPI_CHAIN)
      TYPE (ML_MPI_PAR)                 :: ML_COMM_WORLD_GLOBAL
#endif
      LOGICAL                           :: ML_IO_WRITE             = .FALSE.
      INTEGER                           :: ML_FILESTREAM_GLOBAL    = 0
      INTEGER                           :: ML_TOTNUM_INSTANCES     = 0
      INTEGER                           :: ML_MB_DEFAULT           = 1500
      INTEGER                           :: ML_MCONF_DEFAULT        = 1500

      END MODULE
