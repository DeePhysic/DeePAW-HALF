#include "symbol.inc"
#define usgrid
!****************** PROGRAM VASP  Version 5.0 (f90)********************
! RC:  $Id: main.F,v 1.18 2003/06/27 13:22:18 kresse Exp kresse $
!> Vienna Ab initio total energy and Molecular-dynamics Program
!            written  by   Kresse Georg
!                     and  Juergen Furthmueller
! Georg Kresse                       email: Georg.Kresse@univie.ac.at
! Juergen Furthmueller               email: furth@ifto.physik.uni-jena.de
! Institut fuer Materialphysik         voice: +43-1-4277-51402
! Uni Wien, Sensengasse 8/12           fax:   +43-1-4277-9514 (or 9513)
! A-1090 Wien, AUSTRIA                 http://cms.mpi.univie.ac.at/kresse
!>
!> This program comes without any waranty.
!> No part of this program must be distributed, modified, or supplied
!> to any other person for any reason whatsoever
!> without prior written permission of the Institut of Materials Science
!> University Vienna
!>
!> This program performs total energy calculations using
!> a selfconsistency cylce (i.e. mixing + iterative matrix diagonal.)
!> or a direct optimisation of the one electron wavefunctions
!> most of the algorithms implemented are described in
!> * G. Kresse and J. Furthmueller,
!>   Efficiency of ab-initio total energy calculations for
!>   metals and semiconductors using a plane-wave basis set,
!>   *Comput. Mat. Sci.* **6**, 15-50 (1996)
!> * G. Kresse and J. Furthmueller,
!>   Efficient iterative schemes for ab-initio total energy
!>   calculations using a plane-wave basis set,
!>   *Phys. Rev. B* **54**, 11169 (1996)
!>
!> The iterative matrix diagonalization is based
!> * on the conjugated gradient eigenvalue minimisation proposed by
!>   D.M. Bylander, L. Kleinmann, S. Lee, *Phys Rev. B* **42**, 1394 (1990)
!>   and is a variant of an algorithm proposed by
!>   M.P. Teter, M.C. Payne and D.C. Allan, *Phys. Rev. B* **40**,12255 (1989),
!>   T.A. Arias, M.C. Payne, J.D. Joannopoulos, *Phys Rev. B* **45**,1538(1992).
!> * or the residual vector minimization method (RMM-DIIS) proposed by
!>   P. Pulay, *Chem. Phys. Lett.* **73**, 393 (1980).
!>   D. M. Wood and A. Zunger, *J. Phys. A*, 1343 (1985)
!>
!> For the mixing a Broyden/Pulay like method is used (see for instance):
!> D. Johnson, *Phys. Rev. B* **38**, 12807 (1988)
!>
!> The program can use normconserving PP,
!> generalised ultrasoft-PP (Vanderbilt-PP Vanderbilt *Phys Rev B* **40**,
!> 12255 (1989)) and PAW (P.E. Bloechl, *Phys. Rev. B* **50**, 17953 (1994))
!> datasets. Partial core corrections can be handled.
!> Spin and GGA and exact exchange functionals are implemented
!>
!> The units used in the programs are electron-volts and angstroms.
!> The unit cell is arbitrary, and arbitrary species of ions are handled.
!> A full featured symmetry-code is included, and calculation of
!> Monkhorst-Pack special-points is possible (tetrahedron method can be
!> used as well).
! This part was written by J. Furthmueller.
!
! The original version was written by  M.C. Payne
! at Professor J. Joannopoulos research  group at the MIT
! (3000 lines, excluding FFT, July 1989)
! The program was completely rewritten and vasply extended by
! Kresse Georg (gK) and Juergen Furthmueller. Currently the
! code has about 200000 source lines
!
!** The following parts have been taken from other programs
! - Tetrahedron method (original author unknown)
!
! please refer to the README file to learn about new features
! notes on singe-precision:
! USAGE NOT RECOMMENDED DUE TO FINITE DIFFERENCES IN FEW SPECIFIC
! FORCE-SUBROUTINE
! (except for native 64-bit-REAL machines like CRAY style machines)
!**********************************************************************

      PROGRAM VAMP
      USE prec

      USE charge
      USE pseudo
      USE lattice
      USE steep
      USE us
      USE pawm
      USE pot
      USE force
      USE fileio
      USE nonl_high
      USE rmm_diis
      USE ini
      USE ebs
      USE wave_high
      USE choleski
      USE mwavpre
      USE mwavpre_noio
      USE msphpro
      USE broyden
      USE msymmetry
      USE subrot
      USE melf
      USE base
      USE mpimy
      USE mgrid
      USE mkpoints
      USE constant
      USE setexm
      USE poscar
      USE wave
      USE hamil
      USE main_mpi
      USE command_line
      USE chain
      USE pardens
      USE finite_differences
      USE checkpointing
      USE LDAPLUSU_MODULE
      USE cl
      USE Constrained_M_modular
      USE writer
      USE sym_prec
      USE elpol
      USE mdipol
      USE wannier
      USE vaspxml
      USE full_kpoints
      USE kpoints_change
      USE fock
      USE compat_gga
      USE mlr_main
      USE mlrf_main
      USE mlr_optic
      USE pwkli
      USE gridq
      USE twoelectron4o
      USE dfast
      USE aedens
      USE chi_glb
      USE afqmc, ONLY: afqmc_settings, afqmc_reader, afqmc_propagation
      USE xi
      USE subrotscf
      USE pead
      USE egrad
      USE hamil_struct_def
      USE morbitalmag
      USE mconverse
      USE relativistic
      USE rhfatm
      USE tau_mu
      USE mbj, ONLY : CREATE_CMBJ_AUX,CREATE_CMBJ_AUG,UPDATE_CMBJ,PRINT_CMBJ_AUX
      USE mkproj
      USE classicfields
      USE rpa_force
      USE rpa_high
      USE wap, ONLY: getwave
      USE density_of_states, ONLY: DENINI
      USE half_vasp_init, ONLY: HALF_INITIALIZE_WAVEFUNCTIONS
      USE elphon_base, ONLY: elphon_pw_settings => elph_settings, ELECTRON_PHONON_READER
      USE elphon_driver, ONLY: ELECTRON_PHONON_DRIVER
#ifdef vasp6
      USE chi_super
      USE chi_GG
      USE time_propagation
#endif
! Thomas Bucko's code contributions
      USE random_seeded
#ifdef tbdyn
      USE dynconstr
#endif
      USE vdwforcefield
      USE vdw_sd3_d4, ONLY: VDW_WRITER
      USE plugins, ONLY: PLUGINS_READER, WRITE_PLUGINS, STOP_PLUGINS, &
          PLUGINS_MACHINE_LEARNING, PLUGINS_OCCUPANCIES, PLUGINS_STRUCTURE
      USE version, ONLY: VASP, RELEASE_DATE
      USE license_check, ONLY: CHECK_LICENSE_AND_STOP_IF_INVALID
      USE dimer_heyden
      USE dvvtrajectory

      USE mlwf
      USE umco
      USE dmft
      USE crpa
      USE chgfit
      USE stockholder
      USE mlr_main_nmr
      USE hyperfine
      USE wannier_interpolation
      USE wave_interpolate
      USE auger
      USE dmatrix

      USE lcao
! solvation__
      USE solvation
! solvation__
      USE elphon_util
      USE elphon
      USE phonon
      USE locproj
      USE embed
      USE elphon_derivative
      USE ml_interface
      USE ml_interface_writer
#ifdef ML_AVAILABLE
      USE ml_ff_struct, ONLY : ML_SUPER_TYPE, ML_SUPER_HANDLE_MAIN, ML_TOTNUM_INSTANCES, ML_COMM_WORLD_GLOBAL
      USE ml_reader
#ifdef libvaspml
      USE mvaspml, ONLY : vaspml_wrapper
#endif
#endif
#ifdef PROFILING
      USE profiling
#endif
#ifdef BUILD_INFO
      USE build_info
#endif
! bexternal__
      USE bexternal
! bexternal__
      USE mopenmp, ONLY : INIT_OMP
#ifdef VASP_HDF5
      USE vhdf5
#endif
      USE reader_tags
      USE incar_reader, ONLY: INCAR_FROM_FILE, UNUSED_INCAR_TAGS, str
      USE core_con_mat

#ifdef CUSOLVERMP
      USE mcusolvermp, ONLY : CUSOLVERMP_FREE
#endif
#ifdef CUBLASMP
      USE mcublasmp, ONLY : CUBLASMP_FREE
#endif

__GPU USE moffload
__GPU USE wpbe

#ifdef MinimaxTest
      USE minimax_struct, ONLY : imag_grid_handle, MINIMAX_INTEGRATION_TEST
      USE minimax, ONLY : SET_IMAG_GRID_HANDLE, DEALLOCATE_IMAG_GRID_HANDLE
#endif
      USE tutor, ONLY: vtutor, isError, isAlert, isWarning, isAdvice, newArgument, &
          ISYM_MD, VASP_4_4, DoubleCounting, FFTGridIsNotSufficient, &
          HighestBandsOcccupied, PartialDOS, LVTOT, NBANDSchanged, NumberOfElectrons, &
          NPAREfficiency, LongLattice, IALGO8, LOPTICSAndGW, state_order,&
          NCOREANDGWR, NonCollinearAndGWR, RPAFORCES_COMPAT, &
          CheckpointContinue, CheckpointSingle
      USE scpc, ONLY: SCPC_READER, SCPC_SET
! embedding__
      USE mextpot
! embedding__
      USE esf, ONLY: ESF_READER
      USE coulomb_cutoff, only: change_grid_for_kernel_truncation, &
                                create_grids_for_kernel_truncation, &
                                setup_padding_and_truncation_factor, &
                                validate_kernel_dimensionality, &
                                validate_surface_normal_direction, &
                                write_coulomb_cutoff
      USE msdgw
#ifdef VASP_DEBUG
      USE debug
#endif
      IMPLICIT NONE

!=======================================================================
!  a small set of parameters might be set here
!  but this is really rarely necessary :->
!=======================================================================
!----I/O-related things (adapt on installation or for special purposes)
!     IU6    overall output ('console protocol'/'OUTCAR' I/O-unit)
!     IU0    very important output ('standard [error] output I/O-unit')
!     IU5    input-parameters ('standard input'/INCAR I/O-unit)
!     ICMPLX size of complex items (in bytes/complex item)
!     MRECL  maximum record length for direct access files
!            (if no restictions set 0 or very large number)
      INTEGER,PARAMETER :: ICMPLX=16,MRECL=10000000

!=======================================================================
!  structures
!=======================================================================
      TYPE (potcar),ALLOCATABLE :: P(:)
      TYPE (wavedes)     WDES
      TYPE (nonlr_struct), TARGET :: NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W          ! wavefunction
      TYPE (wavespin)    W_F        ! wavefunction for all bands simultaneous
      TYPE (wavespin)    W_G        ! same as above
      TYPE (wavefun)     WUP
      TYPE (wavefun)     WDW
      TYPE (wavefun)     WTMP       ! temporary
      TYPE (latt)        LATT_CUR
      TYPE (latt)        LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      TYPE (mixing)      MIX
      TYPE (kpoints_struct) KPOINTS
      TYPE (symmetry)    SYMM
      TYPE (grid_3d)     GRID        ! grid for wavefunctions
      TYPE (grid_3d)     GRID_SOFT   ! grid for soft chargedensity
      TYPE (grid_3d)     GRIDC       ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS      ! very find grid temporarily used in us.F
      TYPE (grid_3d)     GRIDB       ! Broyden grid
      TYPE (grid_3d)     GRID_COARSE ! Coarse grid for coulomb cutoff
      TYPE (grid_3d)     GRID_FINE   ! Fine grid for coulomb cutoff
      TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
      TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRID_SOFT and GRIDC
      TYPE (transit)     TABLE_CUTOFF ! index table for coulomb cutoff
      TYPE (prediction)  PRED
      TYPE (smear_struct) SMEAR_LOOP
      TYPE (paco_struct) PACO
      TYPE (energy)      E
      TYPE (ham_handle)  HAMILTONIAN
      TYPE (tau_handle)  KINEDEN
      TYPE (afqmc_settings) AFQMC_SET
      TYPE (prim_cell_t) PRIM_CELL
      TYPE (polar_data_t) POLAR_DATA
      TYPE (phonon_settings) PHON_SETTINGS
      TYPE (elphon_settings) ELPH_SETTINGS
      TYPE (elphon_interpolator) ELPH_INTERPOL
      TYPE (elphon_pw_settings) ELPH_PW_SETTINGS
      TYPE (wave_interpolator_settings) WAVE_INTERP_SET
      TYPE (elphon_deriv_settings) ELPH_DERIV_SETTINGS
      TYPE (checkpoint_settings_struct) CHECKPOINT_FD_SETTINGS
      TYPE (checkpoint_struct) FINDIFF_CHECKPOINT
      TYPE (nh_chains) NHC     ! Nose-Hoover chains
#ifdef libvaspml
      TYPE (vaspml_wrapper) VASPML
#endif
#ifdef VASP_DEBUG
      TYPE (debug_settings) DBG
#endif

#ifdef tbdyn
       TYPE (gadget_io)   g_io
       INTEGER :: SEED(3),SEED_INIT(3)
       INTEGER, PARAMETER :: K_SEED=3
#else
       INTEGER :: MDALGO_GL 
       LOGICAL :: LATTICE_CONSTRAINTS_GL(3)=(/.TRUE.,.TRUE.,.TRUE./)
#endif

      INTEGER :: NGX,NGY,NGZ,NGXC,NGYC,NGZC
      INTEGER :: NRPLWV,LDIM,LMDIM,LDIM2,LMYDIM
      INTEGER :: IRMAX,IRDMAX,ISPIND
      INTEGER :: NPLWV,MPLWV,NPLWVC,MPLWVC,NTYPD,NIOND,NIONPD,NTYPPD
      INTEGER :: NEDOS
      LOGICAL :: LNBANDS=.TRUE.
      INTEGER :: TIU6, TIU0
      LOGICAL :: ALL_LDO_AB_INITIO_AMONG_IMAGES
      LOGICAL :: ANY_LDO_AB_INITIO_AMONG_IMAGES
      LOGICAL :: ANY_MLFF_AMONG_IMAGES
      !INTEGER :: MDALGO=0           ! dublicates MDALGO in tbdyn
!=======================================================================
!  begin array dimensions ...
!=======================================================================
      COMPLEX(q),ALLOCATABLE:: CHTOT(:,:)    ! charge-density in real / reciprocal space
      RGRID  ,ALLOCATABLE:: SV(:,:)          ! soft part of local potential
      COMPLEX(q),ALLOCATABLE:: CHTOTL(:,:)   ! old charge-density
      RGRID     ,ALLOCATABLE:: DENCOR(:)     ! partial core
      COMPLEX(q),ALLOCATABLE:: CVTOT(:,:)    ! local potential
      COMPLEX(q),ALLOCATABLE:: CVTOT_DERIV(:,:,:,:)   ! derivative of local potential
      COMPLEX(q),ALLOCATABLE:: MU_DERIV(:,:,:,:)   ! derivative of tau (kinetic energy density)
      COMPLEX(q),ALLOCATABLE:: CSTRF(:,:)    ! structure-factor
!-----non-local pseudopotential parameters
      OVERLAP,ALLOCATABLE:: CDIJ(:,:,:,:)    ! strength of PP
      OVERLAP,ALLOCATABLE:: CDIJ_DERIV(:,:,:,:,:,:)   ! derivative of PP strength
      OVERLAP,ALLOCATABLE:: CQIJ(:,:,:,:)    ! overlap of PP
      OVERLAP,ALLOCATABLE:: CRHODE(:,:,:,:)  ! augmentation occupancies
!-----elements required for mixing in PAW method
      REAL(q)   ,ALLOCATABLE::   RHOLM(:,:),RHOLM_LAST(:,:)
!-----charge-density and potential on small grid
      COMPLEX(q),ALLOCATABLE:: CHDEN(:,:)    ! pseudo charge density
!-----description how to go from one grid to the second grid
!-----density of states
      REAL(q)   ,ALLOCATABLE::  DOS(:,:),DOSI(:,:)
      REAL(q)   ,ALLOCATABLE::  DDOS(:,:),DDOSI(:,:)
!-----local l-projected wavefunction characters
      REAL(q)   ,ALLOCATABLE::   PAR(:,:,:,:,:),DOSPAR(:,:,:,:)
!  all-band-simultaneous-update arrays
      GDEF   ,POINTER::   CHF(:,:,:,:),CHAM(:,:,:,:)
!  optics stuff
      GDEF   ,ALLOCATABLE::   NABIJ(:,:)
!
      LOGICAL :: LVCADER

!-----remaining mainly work arrays
      COMPLEX(q), ALLOCATABLE,TARGET :: CWORK1(:),CWORK2(:),CWORK(:,:)
      TYPE (wavefun1)    W1            ! current wavefunction
      TYPE (wavedes1)    WDES1         ! descriptor for one k-point

      GDEF, ALLOCATABLE  ::  CPROTM(:),CMAT(:,:)

      REAL(q), ALLOCATABLE :: FORCE_CONST(:,:)
      REAL(q), ALLOCATABLE :: ORBITAL_MOMENT(:,:,:)
!=======================================================================
!  a few fixed size (or small) arrays
!=======================================================================
!-----Forces and stresses
      REAL(q)   VTMP(3), XCSIF(3,3), EWSIF(3,3), TSIF(3,3), D2SIF(3,3)
!-----forces on ions
      REAL(q)   ,ALLOCATABLE::  EWIFOR(:,:), TIFOR(:,:)
!-----data for STM simulation (Bardeen)
      REAL(q)  STM(7)
!-----Temporary data for tutorial messages ...
      INTEGER,PARAMETER :: NTUTOR=1000
      REAL(q)     RTUT(NTUTOR),RDUM
      INTEGER  ITUT(NTUTOR),IDUM
      CHARACTER(20) :: CTUT(1)
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
!=======================================================================
!  end array dimensions ...
!=======================================================================
      INTEGER NTYP_PP      ! number of types on POTCAR file

      INTEGER ITYP
      INTEGER I,J,N,NT,K
      INTEGER NPAR, NODE_ME, NCPU, IONODE, JOBPAR, JOBPAR_
      INTEGER MPLMAX, LMAX_TABLE, MALLOC, LMAX_CALC, LPAR, LPRO, POUT
      INTEGER NBANDS, NPPSTR, NDEGREES_OF_FREEDOM, NCDIJ, NELM, NN, NK
      INTEGER NMAG, NI, NKP, NTYP, NIOND_LOC, N_MIX_PAW, NSTEP, NCOUNT
      INTEGER NPL, NPRO, N1, N2, NRPLWV_RED, NPROD_RED
      INTEGER IRDMAA, ITMP, IOCCUP, ISP, INDMAX, IGPAR, IOCCVS, IUXML_SET
      INTEGER ISCALE, IFLAG, IOCC, IU0, IUT, IERROR
      INTEGER OUTPUT_MODE, LBLOCK_HELP, NSW_ACT
      REAL(q) ENMAXI, EKIN, EKIN_LAT, ECONV, EFERMI, EENTROPY, ES, E1TEST
      REAL(q) EVALUE, EACC, ETOTAL, EN, EXHF, TEIN, TOTEN, TOTENG
      REAL(q) XAU, YAU, ZAU, AOMEGA, OMEGA_OLD, FACT, FACTSI, DISMAX
      REAL(q) QE, QF, SQQ, SQQAU, WOSZI, ELEKTR, RHOMAG, TMP
      REAL(q) TMEAN, TMEAN0, SMEAN, SMEAN0, SCALEE, SCALEQ
      REAL(q) SUMNEL, EPS, PCFAK, PRESS, DELTAE
!---- used for creation of param.inc
      REAL(q)    WFACT,PSRMX,PSDMX
      REAL(q)    XCUTOF,YCUTOF,ZCUTOF
#if defined ML_AVAILABLE || !defined ML_LEARNING_NOT_ACTIVE
      REAL(q)             :: TOTEN_ML, TSIF_ML(3,3)
      REAL(q),ALLOCATABLE :: TIFOR_ML(:,:)
#endif

!---- timing information
      INTEGER IERR

      INTEGER NORDER   !   order of smearing
!---- a few logical and string variables
      LOGICAL    LTMP,LSTOP2
      LOGICAL    LPAW           ! paw is used
      LOGICAL    LREALLOCATE    ! reallocation of proj operators required
      LOGICAL    L_NO_US        ! no ultrasoft PP
      LOGICAL    LADDGRID       ! additional support grid
      LOGICAL    LCHAIN         ! initialize chains
      LOGICAL    LBERRY         ! calculate electronic polarisation
      LOGICAL    LSKIP_ML_ISTART3
      LOGICAL    :: LDELTADFT = .FALSE.

#ifdef libbeef
      LOGICAL    LBEEFENS       ! calculate BEEF-vdW ensemble XC energies
      LOGICAL    LBEEFBAS       ! only print basis energies
#endif

      CHARACTER (LEN=40)  SZ
      CHARACTER (LEN=1)   CHARAC
      CHARACTER (LEN=5)   IDENTIFY
      CHARACTER (LEN=:), ALLOCATABLE :: UNUSED_TAGS

!-----parameters for sphpro.f
      INTEGER :: LDIMP,LMDIMP,LTRUNC=3
#ifdef MinimaxTest
      TYPE(imag_grid_handle) :: IMAG_GRIDS
#endif
!=======================================================================
! All COMMON blocks
!=======================================================================
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX
      COMMON /WAVCUT/ IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX

      REAL(q)  RHOTOT(4)
      INTEGER(qi8) IL,I1,I2_0,I3,I4

      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      LOGICAL, EXTERNAL :: USE_OEP_IN_GW
      REAL(q), EXTERNAL :: ENCUTGW_IN_CHI, RHO0

      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      TYPE (xc_info) :: XC_DATA_TEMP

#ifdef libbeef
      LOGICAL  LBEEFCALCBASIS
      COMMON /BEEFENS/ LBEEFCALCBASIS,LBEEFBAS
#endif
      LOGICAL :: LEXIST_POT
#ifdef VASP_HDF5
      INTEGER :: IH5ERR
      CHARACTER(4) :: POTENTIALTYP
#endif
#ifdef PROFILING
!=======================================================================
!  initialise profiling
!=======================================================================
      CALL INIT_PROFILING
#endif
!=======================================================================
!  initialise / set constants and parameters ...
!=======================================================================

! Variables which require zero-initialization:
      TOTEN = 0.0_q
      TSIF = 0.0_q
#if defined ML_AVAILABLE || !defined ML_LEARNING_NOT_ACTIVE
      TOTEN_ML = 0.0_q
      TSIF_ML = 0.0_q
#endif

#ifdef libbeef
      LBEEFCALCBASIS = .FALSE.
#endif

      IO%LOPEN =.TRUE.  ! open all files with file names
      IO%IU0   = 6
      IO%IU6   = 8
      IO%IU5   = 7

      IO%ICMPLX=ICMPLX
      IO%MRECL =MRECL
      PRED%ICMPLX=ICMPLX

#ifdef tbdyn
      g_io%REPORT=66
      g_io%REFCOORD=67
      g_io%CONSTRAINTS=69
      g_io%STRUCTINPUT=533
      g_io%PENALTY=534
      g_io%IRCPOINTS=535
      g_io%HESSEMAT=536
#endif

! Initiliaze output verbosity flag (mainly for molecular dynamics)
! And output frequency helping variable
      OUTPUT_MODE = 1
      LBLOCK_HELP = 1

! Initialize variable for skipping symmetry routines
      LSKIP_ML_ISTART3 = .FALSE. 
!
! with all the dynamic libraries (blas, lapack, scalapack etc.
! VASP has a size of at least 30 Mbyte
!
      CALL INIT_FINAL_TIMING
      CALL REGISTER_ALLOCATE(30000000._q, "base")

! switch off kill
!     CALL sigtrp()

      NPAR=1
      IUXML_SET=20

#if defined(MPI) || defined(MPI_CHAIN)
      CALL INIT_MPI(NPAR,IO)

      NODE_ME= COMM%NODE_ME
      IONODE = COMM%IONODE

      IF (NODE_ME/=IONODE.OR.KIMAGES>0) THEN
         IUXML_SET=-1
      ENDIF
#else
      CALL PARSE_COMMAND_LINE(IO, SHOULD_WRITE=.TRUE.)
      INCAR_F = INCAR_FROM_FILE(INCAR)
#endif
#ifdef _OPENMP
!-----------------------------------------------------------------------
!  initialise OpenMP: read thread related INCAR tags
!-----------------------------------------------------------------------
      CALL INIT_OMP(IO)
#endif
#ifdef OFFLOADING
!-----------------------------------------------------------------------
!  initialise OpenACC runtime
!-----------------------------------------------------------------------
      CALL INIT_OFFLOAD(IO)
#endif
#ifdef GDBDEBUG
!-----------------------------------------------------------------------
!  command line argument -g or --debug might be passed in order to 
!  wait for debugger to attach to vasp
!-----------------------------------------------------------------------
      IF ( IO%HALT ) THEN
         CALL WAIT_ON_DEBUGGER()
      ENDIF
#endif
#ifdef VASP_HDF5
      ! variable to check if one is on the IO node
      LWRITEH5=(IUXML_SET/=-1)
#if !defined(MPI) && !defined(MPI_CHAIN)
      CALL START_H5()
#endif
      CALL PROCESS_INCAR(INCAR_F, "LSYNCH5", IO%LSYNCH5)
      CALL OPEN_VASPOUT(IO%LSYNCH5, SUBDIR=DIR_APP(1:DIR_LEN))
      CALL PROCESS_INCAR(INCAR_F, "LH5", IO%LH5)
      CALL PROCESS_INCAR(INCAR_F, "LCHARGH5", IO%WRT_DENSITY%LCHARGH5)
      CALL OPEN_VASPWAVE(IO, COMM, SUBDIR=DIR_APP(1:DIR_LEN))
      INCAR_F%TO_HDF5 = .TRUE.
#endif
      TIU6 = IO%IU6
      TIU0 = IO%IU0
      vtutor%unitOut = IO%IU6
      vtutor%unitErr = IO%IU0
      vtutor%discardOutputToNegativeUnit = .TRUE.

      CALL START_XML( IUXML_SET, DIR_APP(1:DIR_LEN)//"vasprun.xml" )

!-----------------------------------------------------------------------
!  open Files
!-----------------------------------------------------------------------
      IF (IO%IU0>=0) WRITE(TIU0,*) VASP()
#if defined(makeparam)
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTPAR',STATUS='UNKNOWN')
#elif defined(makekpoints)
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTKPT',STATUS='UNKNOWN')
#else
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTCAR',STATUS='UNKNOWN')
#endif
      ! before we do anything else check license
      CALL CHECK_LICENSE_AND_STOP_IF_INVALID(IO, COMM_WORLD, RELEASE_DATE)

#ifdef logflow
      OPEN(UNIT=19,FILE=DIR_APP(1:DIR_LEN)//'JOBFLOW',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE(19,'(A)') 'SECTION = "----------------------------------------------------------------------------------------------------"'
#endif

#ifndef makekpoints
      CALL OPENWAV(IO, COMM)

#ifndef makeparam
      io_begin
      IF (KIMAGES==0) THEN
      OPEN(UNIT=22,FILE=DIR_APP(1:DIR_LEN)//'EIGENVAL',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE=DIR_APP(1:DIR_LEN)//'CONTCAR',STATUS='UNKNOWN')
      OPEN(UNIT=16,FILE=DIR_APP(1:DIR_LEN)//'DOSCAR',STATUS='UNKNOWN')
      OPEN(UNIT=17,FILE=DIR_APP(1:DIR_LEN)//'OSZICAR',STATUS='UNKNOWN')
      OPEN(UNIT=60,FILE=DIR_APP(1:DIR_LEN)//'PCDAT',STATUS='UNKNOWN')
      OPEN(UNIT=61,FILE=DIR_APP(1:DIR_LEN)//'XDATCAR',STATUS='UNKNOWN')
      OPEN(UNIT=70,FILE=DIR_APP(1:DIR_LEN)//'CHG',STATUS='UNKNOWN')
#ifdef tbdyn
      OPEN(UNIT=g_io%REPORT,FILE=DIR_APP(1:DIR_LEN)//'REPORT',STATUS='UNKNOWN')
#endif
      ENDIF
      io_end
#endif
      IF (IO%IU6>=0) WRITE(IO%IU6,*) VASP()
      CALL XML_GENERATOR
#if defined(MPI)
      CALL PARSE_GENERATOR_XML(VASP()//" parallel")
#else
      CALL PARSE_GENERATOR_XML(VASP()//" serial")
#endif
      CALL MY_DATE_AND_TIME(IO%IU6)
      CALL XML_CLOSE_TAG

      CALL WRT_DISTR(IO%IU6)
#endif

! unit for extrapolation of wavefunction
      PRED%IUDIR =21
! unit for broyden mixing
      MIX%IUBROY=23
! unit for total potential
      IO%IUVTOT=62

#ifdef BUILD_INFO
      CALL dump_build_info(IO)
#endif

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

 140  FORMAT (5X, //, &
     &'--------------------------------------- Ionic step ', &
     &I8,'  -------------------------------------------'//)
!-----------------------------------------------------------------------
! read header of POSCAR file to get NTYPD, NTYPDD, NIOND and NIONPD
!-----------------------------------------------------------------------
      CALL RD_POSCAR_HEAD(LATT_CUR, T_INFO, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

      ALLOCATE(T_INFO%ATOMOM(3*NIOND),T_INFO%RWIGS(NTYPPD),T_INFO%ROPT(NTYPD),T_INFO%POMASS(NTYPD), &
               T_INFO%DARWIN_V(NTYPD),T_INFO%DARWIN_R(NTYPD),T_INFO%VCA(NTYPD),DYN%EFOR(3,NIOND))
      IF (IO%IU6>=0) THEN
         WRITE(TIU6,130)
         IF (INCAR_F%ERROR /= "") &
            CALL vtutor%alert("Found issues when parsing the INCAR file\n" // INCAR_F%ERROR // &
                "\nPlease check if the interpretation in the OUTCAR file is correct.")
         WRITE(TIU6,*)'INCAR:'
         WRITE(TIU6,'(a)') str(INCAR_F, INDENTATION=3)
      ENDIF
!  first scan of POTCAR to get LDIM, LMDIM, LDIM2 ...
      LDIM =16
      LDIM2=(LDIM*(LDIM+1))/2
      LMDIM=64

      ALLOCATE(P(NTYPD))
      T_INFO%POMASS=0
      T_INFO%RWIGS=0
      T_INFO%VCA=0

      INFO%NLSPLINE=.FALSE.
!-----------------------------------------------------------------------
! read pseudopotentials
!-----------------------------------------------------------------------
      CALL ALLOCATE_XC_DATA(XC_DATA_PP,1)
  ! check whether generation of potcar.h5 is requested
#ifdef VASP_HDF5
      CALL POTCAR_H5_READER(IO%IU5, IO%IU0, POTENTIALTYP)
#endif
      CALL RD_PSEUDO(INFO,P, &
     &           NTYP_PP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           T_INFO%POMASS,T_INFO%RWIGS,0,T_INFO%TYPE,T_INFO%VCA, &
     &           IO%IU0,IO%IU6,-1,LPAW &
#ifdef VASP_HDF5
     &           , POTENTIALTYP, T_INFO%TYPEF, T_INFO%SHA256 &
#endif
     &)
!-----------------------------------------------------------------------
! read INCAR
!-----------------------------------------------------------------------
      CALL XML_TAG("incar")
      WRITEXMLINCAR = .TRUE.
      INCAR_F%TO_XML = .TRUE.
      CALL READER( &
          IO%IU5,IO%IU0,IO%INTERACTIVE,INFO%SZNAM1,INFO%ISTART,INFO%LHALF_INIT,INFO%LHALF_API,INFO%HALF_ESCN_URL,INFO%IALGO,INFO%IALGO_COMPAT,INFO%LROBUST,INFO%ISEARCH,MIX%IMIX,MIX%MAXMIX,MIX%MREMOVE, &
          MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
          MIX%WC,MIX%INIMIX,MIX%MIXPRE,MIX%MIXFIRST,IO%LFOUND,INFO%LDIAG,INFO%LSUBROT,INFO%LREAL,IO%LREALD,IO%LPDENS, &
          DYN%IBRION,INFO%ICHARG,INFO%INIWAV,INFO%NELM,INFO%NELMALL,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,DYN%EDIFFG, &
          DYN%NSW,DYN%ISIF,PRED%IWAVPR,PRED%MAXDIS,SYMM%ISYM,DYN%NBLOCK,DYN%KBLOCK,INFO%ENMAX,DYN%POTIM,DYN%TEBEG, &
          DYN%TEEND,DYN%NFREE,DYN%EFOR, &
          PACO%NPACO,PACO%APACO,T_INFO%NTYP,NTYPD,DYN%SMASS,SCALEE,T_INFO%POMASS, &
          T_INFO%DARWIN_V,T_INFO%DARWIN_R,T_INFO%VCA,LVCADER, &
          T_INFO%RWIGS,INFO%NELECT,INFO%NUP_DOWN,INFO%TIME, &
          KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%EFERMI,KPOINTS%FERMI_METHOD,KPOINTS%ISMEAR, &
          KPOINTS%BANDGAP,KPOINTS%SPACING,KPOINTS%LGAMMA,KPOINTS%LKBLOWUP, &
          DYN%PSTRESS,INFO%NDAV,INFO%NKRYLOV,INFO%DAVITER,INFO%NSTRICT,INFO%FBREAK, &
          KPOINTS%SIGMA,KPOINTS%EFERMI_NEDOS,KPOINTS%LTET,INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,IO%NWRITE,INFO%LCORR, &
          IO%IDIOT,T_INFO%NIONS,T_INFO%NTYPP,INFO%RADIAL_FDORDER,IO%LMUSIC,IO%LOPTICS,STM, &
          INFO%ISPIN,T_INFO%ATOMOM,NIOND,IO%LWAVE,IO%LDOWNSAMPLE,INFO%SZPREC, &
          INFO%ENAUG,IO%LORBIT,IO%LELF,T_INFO%ROPT,INFO%ENINI, &
          NGX,NGY,NGZ,NGXC,NGYC,NGZC,NBANDS,WDES%NBANDSLOW,WDES%NBANDSHIGH,NEDOS,NBLK,LATT_CUR, &
          LPLANE_WISE,LCOMPAT,LMAX_CALC,SET_LMAX_MIX_TO,WDES%NSIM,LPAW,LADDGRID, &
          WDES%LNONCOLLINEAR,WDES%LSORBIT,WDES%SAXIS, &
          WDES%LSPIRAL,WDES%LZEROZ,WDES%QSPIRAL,WDES%LORBITALREAL, &
          INFO%LASPH,INFO%TURBO,INFO%IRESTART,INFO%NREBOOT,INFO%NMIN,INFO%EREF, &
          INFO%NLSPLINE,FFTW_PLAN_EFFORT, &
          IO%LH5,IO%WRT_POTENTIAL,IO%WRT_DENSITY,IO%LWAP,IO%VELOCITY, &
          HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE, &
          EWALD_PARAMS%MAX_ARGUMENT &
#ifdef libbeef
         ,LBEEFENS,LBEEFBAS &
#endif
         )

#ifdef VASP_DEBUG
      CALL DEBUG_READER(DBG, IO%LOPEN, IO%IU5, IO%IU0)
#endif

      CALL OMP_READER(IO%IU5,IO%IU0)
      CALL MPI_READER(IO%IU5,IO%IU0)

      !c read in the MD related input parameters (INCAR) and data (ICONST, PENALTYPOT)
      !IF (DYN%IBRION==0) THEN
      !  CALL DYNCONSTR_INIT(DYN,T_INFO,INFO,LATT_CUR, IO,g_io)
      !ENDIF

      CALL RANDOM_READER(INCAR_F, COMM, SEED, WDES%RANDOM_NUMBER_GENERATOR, WDES%PCG_SEED)
#ifdef tbdyn
      ! in case molecular dynamics is selected with symmetry operations on
      IF (DYN%IBRION==0 .AND. SYMM%ISYM>0 ) THEN
         CALL vtutor%write(isAlert, ISYM_MD)
      ENDIF
      CALL RANE_ION(RDUM,PUT=SEED(1:K_SEED))
      SEED_INIT=SEED
#endif

      CALL KPOINTS_SET_NKP_FROM_KSPACING(KPOINTS,LATT_CUR%BNORM,KPOINTS%SPACING)

      CALL GGA_COMPAT_MODE(IO%IU5, IO%IU0, LCOMPAT)

      IF (LCOMPAT) THEN
         CALL vtutor%write(isAlert, VASP_4_4)
      ENDIF
! The meaning of LVTOT has changed w.r.t. previous VASP version,
! therefore we write a warning
      IF (IO%WRT_POTENTIAL%LVTOT) THEN
         CALL vtutor%write(isWarning, LVTOT)
      ENDIF

#if defined(MPI) || defined(MPI_CHAIN)
#ifndef _OPENMP
      IF ( REAL(COMM%NCPU,q)/COMM_INB%NCPU/COMM_INB%NCPU>4 .OR. &
           REAL(COMM%NCPU,q)/COMM_INB%NCPU/COMM_INB%NCPU<0.25_q) THEN
         CALL vtutor%write(isAlert, NPAREfficiency)
      ENDIF
#endif
#endif
!-----------------------------------------------------------------------
! core level shift related items (parses INCAR)
!-----------------------------------------------------------------------
      CALL INIT_CL_SHIFT(IO%IU5,IO%IU0, T_INFO%NIONS, T_INFO%NTYP )
! Berry phase read INCAR
      CALL READER_ADD_ON(IO%IU5,IO%IU0,LBERRY,IGPAR,NPPSTR, &
            INFO%ICHARG,KPOINTS%ISMEAR,KPOINTS%SIGMA)
! Do we want to write the AE-densities?
      CALL INIT_AEDENS(IO%IU0,IO%IU5)

      ISPIND=INFO%ISPIN

      DYN%TEMP =DYN%TEBEG
      INFO%RSPIN=3-INFO%ISPIN

      CALL RESPONSE_READER(INFO, IO)
#ifdef vasp6
      ! euclidean, aka low-scaling GW requires NCORE=1 and Collinear setup 
      IF( ICHIREAL > 1 ) THEN
#ifdef MPI
         IF( NCORE > 1 ) THEN
             CALL VTUTOR%WRITE(ISERROR, NCOREANDGWR, NEWARGUMENT( IVAL=[NCORE]) )
         ENDIF
#endif 
         IF( WDES%LNONCOLLINEAR ) THEN
             CALL vtutor%write(isError, NonCollinearAndGWR )
         ENDIF
      ELSE
         ! euclidean QFT code must be selected 
         IF( ICHIREAL < 1 .AND. LRPAFORCE) THEN
            CALL VTUTOR%WRITE(ISERROR, RPAFORCES_COMPAT)
         ENDIF 
      ENDIF
#endif
      CALL AFQMC_READER(INCAR_F, MAX(IMAGES, 1), AFQMC_SET)
      CALL DMATRIX_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL XC_FOCK_READER(XC_DATA, IO%IU5, IO%IU0, IO%IU6, INFO%SZPREC, DYN%ISIF, SYMM%ISYM, INFO%IALGO, INFO%ICHARG, &
         LATT_CUR%OMEGA, T_INFO%NTYP, T_INFO%NIONS, MIX%IMIX, MIX%AMIX, MIX%AMIX_MAG, MIX%BMIX, &
         MIX%BMIX_MAG, IO%WRT_POTENTIAL%LVTOT, AFQMC_SET, INFO%LASPH, WDES%LNONCOLLINEAR, WDES%SAXIS)
      CALL MGGA_READER(INCAR_F,XC_DATA,IO)
      CALL EGRAD_READER(IO%IU5, IO%IU6, IO%IU0, T_INFO%NTYP)
      CALL CLASSICFIELDS_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL HYPERFINE_READER(IO%IU5, IO%IU6, IO%IU0, T_INFO%NTYP)
      CALL CORE_CON_MAT_READER(IO%IU5,IO%IU0)
      CALL MONTE_CARLO_READER(IO%IU5,IO%IU0)
! First read if machine learning is exectued or not and the number of instances
      CALL MACHINE_LEARNING_GENERAL_AND_INSTANCE_READER(INCAR_F)
#ifdef ML_AVAILABLE
      IF (ML_LMLFF) THEN
! Allocate all instances of the machine learning
         CALL MACHINE_LEARNING_HANDLE_ALLOCATE(ML_TOTNUM_INSTANCES)
! Loop over all instances
         DO I = 1, ML_TOTNUM_INSTANCES
! Set the machine learning instance variable
            ML_SUPER_HANDLE_MAIN(I)%INSTANCE = I
! Pass units for consistent screen output from inside ML routines
            ML_SUPER_HANDLE_MAIN(I)%IO_IU0 = IO%IU0
            ML_SUPER_HANDLE_MAIN(I)%TIU0 = TIU0
#ifndef ML_LEARNING_NOT_ACTIVE
            CALL INQUIRE_ML_AB_FILE(DIR_APP, DIR_LEN, ML_SUPER_HANDLE_MAIN(I))
#endif
            CALL INQUIRE_ML_FF_FILE(DIR_APP, DIR_LEN, ML_SUPER_HANDLE_MAIN(I))
! Read incar tags for all instances
            CALL MACHINE_LEARNING_READER(               &
               ML_LMLFF,                                &
               INCAR_F,                                 &
               ML_SUPER_HANDLE_MAIN(I)%TAG_LIST,        &
               ML_SUPER_HANDLE_MAIN(I)%FFM,             &
               ML_SUPER_HANDLE_MAIN(I)%FF,              &
               ML_SUPER_HANDLE_MAIN(I)%EMP,             &
               ML_SUPER_HANDLE_MAIN(I)%ML_INPUT_HANDLE, &
               T_INFO%NTYP,                             &
               PRED%IWAVPR,                             &
               ML_SUPER_HANDLE_MAIN(I)%INSTANCE,        &
               ML_TOTNUM_INSTANCES,                     &
               DYN%NSW,                                 &
               DYN%IBRION)
! After calling the Fortran ML reader it is safe to assume that FF%LIB is only .TRUE. if
! the code was compiled with the VASPml library.
! Here set skipping of a symmetry routine if needed
            IF (ML_SUPER_HANDLE_MAIN(I)%FF%ISTART.EQ.3) THEN
               LSKIP_ML_ISTART3 = .TRUE.
               SYMM%ISYM = -1
               IF (T_INFO%NIONS.GT.1000) THEN
                  CALL vtutor%alert("You are using ML_MODE=select with a POSCAR file that &
                       contains " // str(T_INFO%NIONS) // " atoms. Since the current POSCAR file &
                       is not relevant for ML_MODE=select, but at the same time for other &
                       parts in VASP, this can lead to unnecessary waiting times. So please &
                       consider using a small POSCAR file.")
               ENDIF
            ENDIF
         ENDDO
#ifdef libvaspml
! If VASPml library is requested call also its reader here.
! Note: The Fortran reader is executed first, because it determines whether the library is compiled,
! is requested and the desired features are available in the library.
         IF (ML_SUPER_HANDLE_MAIN(1)%FF%LIB) THEN
            CALL VASPML%reader(INCAR_F,                               &
                               T_INFO%NTYP,                           &
                               DYN%NSW,                               &
                               DYN%IBRION,                            &
                               ML_SUPER_HANDLE_MAIN(1)%FF%LMLABEXIST, &
                               ML_SUPER_HANDLE_MAIN(1)%FF%LMLFFEXIST)
         ENDIF
#endif
! Here we overwrite the output verbosity flag by the one from the first instance
         OUTPUT_MODE=ML_SUPER_HANDLE_MAIN(1)%FF%OUTPUT_MODE
         IF (DYN%IBRION==0) THEN
            LBLOCK_HELP=ML_SUPER_HANDLE_MAIN(1)%FF%OUTBLOCK
         ELSE
            IF (ML_SUPER_HANDLE_MAIN(1)%FF%OUTBLOCK>1) &
               CALL vtutor%bug("ML_OUTBLOCK>0 is not allowed for IBRION/=0.", __FILE__, __LINE__)
         ENDIF
! If ISTART_FF=2 (this automatically means ML_FF_LMLONLY=.FALSE.), the program does not execute any ab initio calculations.
! In this case, LDO_AB_INITIO is set as .FALSE., here, too.
! For the moment only one instance of machine learning is run, so we set
! ML_FF_LMLONLY to LMLONLY of that instance
         ML_FF_LMLONLY = ML_SUPER_HANDLE_MAIN(1)%FF%LMLONLY
         LDELTADFT = ML_SUPER_HANDLE_MAIN(1)%FF%DELTADFT
! Here we need to reset VASP parameters ML_ISTART=2
         IF (ML_FF_LMLONLY .AND. (.NOT.LDELTADFT)) THEN
            INFO%ENMAX=1 
            NGX=1
            NGY=1
            NGZ=1
            NGXC=1
            NGYC=1
            NGZC=1 
         ENDIF

         ! Check here the type information on POTCAR and POSCAR
         ! They must be the same
         DO ITYP = 1, NTYPD
            IF (T_INFO%TYPE(ITYP) .NE. P(ITYP)%ELEMENT) THEN
               IF (LDELTADFT) THEN
                  CALL vtutor%error("Type information on POSCAR: " &
                       //T_INFO%TYPE(ITYP)// " and POTCAR: " &
                       //P(ITYP)%ELEMENT// " do not agree. &
                       They must be the same for &
                       ML_MODE = delta. &
                       Please fix either the POSCAR file &
                       or the entry VRHFIN on the POTCAR file.") 
               ENDIF
            ENDIF
         ENDDO
      ENDIF
#endif

      CALL PLUGINS_READER(INCAR_F, DYN%IBRION, INFO%PLUGIN)

#ifndef NO_MACHINE_LEARNING
      IF (INFO%PLUGIN%FORCE_FIELD) THEN
         IF (INFO%PLUGIN%ML_MODE .EQ. "delta") THEN
            LDELTADFT = .TRUE. 
         ELSE
            INFO%ENMAX=1 
            NGX=1
            NGY=1
            NGZ=1
            NGXC=1
            NGYC=1
            NGZC=1 
         ENDIF
      ENDIF
#endif

!-----------------------------------------------------------------------
! Parameters and validators for the kernel truncation method
!-----------------------------------------------------------------------
      IF (HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%ACTIVE) THEN
          CALL SETUP_PADDING_AND_TRUNCATION_FACTOR(LATT_CUR,HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE)
          CALL VALIDATE_KERNEL_DIMENSIONALITY(HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE)
          CALL VALIDATE_SURFACE_NORMAL_DIRECTION(HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE)
     END IF

!-----------------------------------------------------------------------
! loop over different smearing parameters
!-----------------------------------------------------------------------
#ifndef makekpoints
      SMEAR_LOOP%ISMCNT=0
      IF (KPOINTS%ISMEAR==-3) THEN
        IF(IO%IU6>=0)   WRITE(TIU6,7219)
 7219   FORMAT('   Loop over smearing-parameters in INCAR')
        ALLOCATE(SMEAR_LOOP%SMEARS(2*DYN%NSW))
        CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'SMEARINGS', SMEAR_LOOP%SMEARS, 2*DYN%NSW, IERR, WRITEXMLINCAR)
        IF ((IERR==0).AND. (MOD(DYN%NSW,2)/=0)) THEN
           CALL vtutor%error("Error reading item 'SMEARINGS' from file INCAR.")
        ENDIF
        SMEAR_LOOP%ISMCNT=DYN%NSW/2+1
        DYN%NSW   =SMEAR_LOOP%ISMCNT+1
        DYN%KBLOCK=DYN%NSW
        KPOINTS%LTET  =.TRUE.
        DYN%IBRION=-1
        KPOINTS%ISMEAR=-5
      ENDIF
!=======================================================================
!  now read in Pseudopotential
!  modify the potential if required (core level shifts)
!=======================================================================
      LMDIM=0
      LDIM=0
      DO NT=1,NTYP_PP
        LMDIM=MAX(LMDIM,P(NT)%LMMAX)
        LDIM =MAX(LDIM ,P(NT)%LMAX)
      END DO
      CALL DEALLOC_PP(P,NTYP_PP)

      LDIM2=(LDIM*(LDIM+1))/2
      LMYDIM=9
! second scan with correct setting
      CALL RD_PSEUDO(INFO,P, &
     &           NTYP_PP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           T_INFO%POMASS,T_INFO%RWIGS,INFO%RADIAL_FDORDER,T_INFO%TYPE,T_INFO%VCA, &
     &           IO%IU0,IO%IU6,IO%NWRITE,LPAW &
#ifdef VASP_HDF5
     &           , POTENTIALTYP, T_INFO%TYPEF, T_INFO%SHA256 &
#endif
      )

      CALL CL_MODIFY_PP( NTYP_PP, P, INFO%ENAUG )
! now check everything
      CALL POST_PSEUDO(NTYPD,NTYP_PP,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,T_INFO%VCA,P,XC_DATA,INFO, &
     &        IO%LREALD,T_INFO%ROPT, IO%IDIOT,IO%IU6,IO%IU0,LMAX_CALC,L_NO_US,WDES%LSPIRAL)
      CALL LDIM_PSEUDO(IO%LORBIT, NTYPD, P, LDIMP, LMDIMP)
! check the difference between ENINI and ENMAX for spin spiral calculations
      IF (WDES%LSPIRAL) CALL CHECK_SPIRAL_ENCUT(WDES,INFO,LATT_CUR,IO)
!-----------------------------------------------------------------------
! LDA+U initialisation (parses INCAR)
!-----------------------------------------------------------------------
      CALL LDAU_READER(T_INFO%NTYP,IO%IU5,IO%IU0)
      IF ( (LGW0 .OR. LGW .OR. LscQPGW) .AND. USELDApU()) THEN
         CALL vtutor%write(isAlert, DoubleCounting)
      ENDIF

      IF (USELDApU().OR.LCALC_ORBITAL_MOMENT()) &
     &   CALL INITIALIZE_LDAU(T_INFO%NIONS,T_INFO%NTYP,P,WDES%LNONCOLLINEAR,IO%IU0,IO%IDIOT)

      CALL SET_PAW_AUG(T_INFO%NTYP, P, IO%IU6, LMAX_CALC, LCOMPAT)
!-----------------------------------------------------------------------
! optics initialisation (parses INCAR)
!-----------------------------------------------------------------------
      IF (IO%LOPTICS) CALL SET_NABIJ_AUG(P,T_INFO%NTYP)
#endif
!-----------------------------------------------------------------------
!  Read UNIT=15: POSCAR Startjob and Continuation-job
!-----------------------------------------------------------------------
      CALL RD_POSCAR(LATT_CUR, T_INFO, DYN, NHC, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, &
     &           IO%IU0,IO%IU6)

!-----------------------------------------------------------------------
! diverse INCAR readers
!-----------------------------------------------------------------------
      CALL RESPONSE_SET_ENCUT(INFO%ENMAX)
      CALL CONSTRAINED_M_READER(T_INFO,WDES,IO%IU0,IO%IU5)
      CALL WRITER_READER(IO%IU0,IO%IU5)
      CALL WANNIER_READER(IO%IU0,IO%IU5,IO%IU6,IO%IDIOT)
      CALL FIELD_READER(T_INFO,P,LATT_CUR,INFO%NELECT,IO%IU0,IO%IU5,IO%IU6)
      CALL LR_READER(INFO%EDIFF,IO%IU0,IO%IU5,IO%IU6)
#ifdef gammareal
      IF (LEPSILON) &
        CALL vtutor%error("The Gamma-only version (vasp_gam) does not &
        &support the use of LEPSILON = .TRUE.. That is because some linear &
        &response routines require a complex shift to obtain stable &
        &convergence. Instead, use the standard version (vasp_std) to &
        &run this calculation.")
#endif
      CALL ORBITALMAG_READER(IO%IU0, IO%IU5, IO%IU6, T_INFO%NIONS, XC_DATA, WDES%SAXIS)
      CALL RHFATM_READER(IO)
      CALL CHGFIT_READER(IO,T_INFO%NTYP)
      CALL STOCKHOLDER_READER(IO)
      CALL CRPA_READER( INFO, COMM, IO )
      CALL MLWF_READER(INFO,IO%IU5,IO%IU6,IO%IU0)
      CALL WAVE_INTERPOLATOR_READER(INCAR_F,WAVE_INTERP_SET,IO)
      CALL LPRJ_READER(T_INFO,P,WDES%LNONCOLLINEAR,IO,SYMM%ISYM)
      CALL PEAD_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL UMCO_READER(INFO, IO)
! embedding__
      CALL EXTPT_READER(NGXC,NGYC,NGZC,IO%IU0,IO%IU5,IO%IU6)
! embedding__
      CALL SCPC_READER(INCAR_F, SCPC_SET)
      !c read in the IVDW related input parameters (INCAR)
      !c this includes also the D3 and D4 methods 
      CALL IVDW_INIT(IO,XC_DATA,DYN,T_INFO,LATT_CUR,P%ELEMENT,KPOINTS)
      CALL EMBED_READER(IO)
      CALL AUGER_READER(IO)
      CALL WANNIER_INTERPOL_READER(IO)

      CALL FINDIFF_READER(INCAR_F, CHECKPOINT_FD_SETTINGS)
      CALL PHON_READER(INCAR_F, PHON_SETTINGS, POLAR_DATA)
      CALL ELPH_READER(INCAR_F, ELPH_SETTINGS)
! electronic structure factor
      CALL ESF_READER(INCAR_F)
      CALL ELPHON_DERIV_READER(INCAR_F, ELPH_DERIV_SETTINGS)
      CALL ELPH_CHECK_INPUT_CONSISTENCY(INCAR_F, PHON_SETTINGS, ELPH_SETTINGS, POLAR_DATA)
      CALL ELECTRON_PHONON_READER(INCAR_F,ELPH_PW_SETTINGS,XC_DATA,WAVE_INTERP_SET,IO)
! solvation__
      CALL SOL_READER(T_INFO%NIONS,INFO%EDIFF,IO)
! solvation__
! bexternal__
      CALL BEXT_READER(IO%IU0,IO%IU5,INFO%ISPIN,WDES%LNONCOLLINEAR)
! bexternal__
      CALL CLASSICFIELDS_WRITE(IO%IU6)
#ifdef tbdyn
      !c read in the MD related input parameters (INCAR) and data (ICONST, PENALTYPOT)
      !c MDALGO is read in here.
      !c Read BBK structure optimization parameters, which is done in STEP_tb
      IF (DYN%IBRION==0 .OR. DYN%IBRION==4) THEN
        CALL DYNCONSTR_INIT(DYN,T_INFO,INFO,LATT_CUR,NHC,IO,g_io)
      ELSE IF ( DYN%IBRION==1 .OR. DYN%IBRION==2 ) THEN
        !! reading LATTICE_CONSTRAINTS for modes IBRION 1 and 2 to do selective cell shape optimization
        CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'LATTICE_CONSTRAINTS', LATTICE_CONSTRAINTS_GL, 3, IERR, WRITEXMLINCAR)
        IF (IERR==3) LATTICE_CONSTRAINTS_GL = .TRUE.
 
        !> 
        !> if all three directions are restricted but cell is allowed to change
        !> 
        IF( (DYN%ISIF==3 .OR. DYN%ISIF==4 .OR. DYN%ISIF==5) .AND. .NOT. ANY(LATTICE_CONSTRAINTS_GL) ) THEN
           CALL VTUTOR%ERROR( "ISIF=3, 4 or 5 and LATTICE_CONSTRAINTS= F F F &
              &contradict each other. Please review your INCAR tags ISIF and LATTICE_CONSTRAINTS.")
        ENDIF

        IF ( DYN%ISIF>5 .AND. .NOT. ALL(LATTICE_CONSTRAINTS_GL) ) THEN
           CALL VTUTOR%ERROR( "ISIF>5 and constraining a lattice vector via LATTICE_CONSTRAINTS is not supported.&
                              & Please remove the LATTICE_CONSTRAINTS tag from INCAR file." )
        END IF

      ENDIF
#endif

      IF (ML_LMLFF .AND. INFO%PLUGIN%FORCE_FIELD) THEN
         CALL vtutor%error("ML_LMLFF and PLUGINS/MACHINE_LEARNING are both set to .TRUE.. &
              There can only be one enabled at a time.")
      ENDIF

#ifndef NO_MACHINE_LEARNING
! Set here LBLOCK_HELP if ML plugin is true
      IF (INFO%PLUGIN%FORCE_FIELD) THEN
         IF (INFO%PLUGIN%OUTBLOCK > MAX(DYN%NSW,1)) THEN
            CALL vtutor%alert("PLUGINS/ML_OUTBLOCK > max(NSW,1), and will be reset to &
                              &PLUGINS/ML_OUTBLOCK = "//str(MAX(DYN%NSW,1))//".") 
            INFO%PLUGIN%OUTBLOCK = MAX(DYN%NSW,1)
         ENDIF
         IF (INFO%PLUGIN%OUTBLOCK>1 .AND. DYN%IBRION/=0) THEN
            CALL vtutor%alert("INFO%PLUGIN%OUTBLOCK > 1 is only supported for IBRION = 0, &
                              &defaulting back to INFO%PLUGIN%OUTBLOCK = 1 ...")
            INFO%PLUGIN%OUTBLOCK = 1
         ENDIF
         LBLOCK_HELP = INFO%PLUGIN%OUTBLOCK
         OUTPUT_MODE = INFO%PLUGIN%OUTPUT_MODE
      ENDIF
#endif
!-----------------------------------------------------------------------
! exchange correlation table
!-----------------------------------------------------------------------
      CALL PUSH_XC_TYPE_FOR_GW(XC_DATA) ! switch now to AEXX=1.0 ; ALDAC = 0.0
 ! switch now to AEXX=1.0 ; ALDAC = 0.0
      IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
         CALL SETUP_LDA_XC(XC_DATA,2,IO%IU6,IO%IU0,IO%IDIOT)
      ELSE
         CALL SETUP_LDA_XC(XC_DATA,1,IO%IU6,IO%IU0,IO%IDIOT)
      ENDIF
!-----------------------------------------------------------------------
! init all chains (INCAR reader)
!-----------------------------------------------------------------------
      LCHAIN = IMAGES > 0 .AND. .NOT.AFQMC_SET % ACTIVE
      IF (LCHAIN) CALL chain_init( T_INFO, IO)
!-----------------------------------------------------------------------
!xml finish copying parameters from INCAR to xml file
! no INCAR reading from here
      CALL XML_CLOSE_TAG("incar")
      writexmlincar = .FALSE.
      INCAR_F%TO_XML = .FALSE.
      UNUSED_TAGS = UNUSED_INCAR_TAGS(INCAR_F)
      IF (UNUSED_TAGS /= "" .AND. IO%NWRITE>=3) &
         CALL vtutor%alert("The following INCAR tags have not been used:\n" // UNUSED_TAGS &
            // "\nPlease check whether the tags are spelled correctly, have the correct &
            &format, and do not appear multiple times in the file.")

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Check for variables between images
!-----------------------------------------------------------------------
#ifdef MPI
      IF (LVCAIMAGES) THEN
         CALL SYNC_VCAIMAGES(DYN, LBLOCK_HELP, COMM_WORLD)
#ifdef ML_AVAILABLE
         IF (ML_LMLFF) THEN
            IF (.NOT. ML_SUPER_HANDLE_MAIN(1)%FF%LIB) THEN
               DO I = 1, ML_TOTNUM_INSTANCES
                  CALL MACHINE_LEARNING_POST_SYNC_VCAIMAGES(ML_SUPER_HANDLE_MAIN(I), LBLOCK_HELP)
               ENDDO
            ENDIF
         ENDIF
#endif
      ENDIF
#endif
!-----------------------------------------------------------------------

__GPU CALL BLAS_OFFLOAD_INIT()
__GPU CALL LAPACK_OFFLOAD_INIT()

      CALL COUNT_DEGREES_OF_FREEDOM( T_INFO, NDEGREES_OF_FREEDOM, &
          IO%IU6, IO%IU0, DYN%IBRION)

!-----for static calculations or relaxation jobs DYN%VEL is meaningless
      IF (DYN%INIT == -1) THEN
        CALL INITIO(T_INFO%NIONS,T_INFO%LSDYN,NDEGREES_OF_FREEDOM, &
               T_INFO%NTYP,T_INFO%ITYP,DYN%TEMP, &
               T_INFO%POMASS,DYN%POTIM, &
               DYN%POSION,DYN%VEL,T_INFO%LSFOR,LATT_CUR%A,LATT_CUR%B,DYN%INIT, COMM, IO%IU6)
#ifdef tbdyn
        ! TB dynamics will initialize the velocities for MDALGO_GL>0
        IF (MDALGO_GL==0) DYN%INIT=0
#else
        DYN%INIT=0
#endif
      ENDIF

      IF (DYN%IBRION/=0 .AND. DYN%IBRION/=40 .AND. DYN%IBRION/=44) THEN
          DYN%VEL=0._q
      ENDIF
      IF (IO%IU6>=0) THEN
         WRITE(TIU6,*)
         WRITE(TIU6,130)
      ENDIF

#ifdef tbdyn
      !c some MD methods (e.q. Langevin dynamics) do not conserve total momentum!
      IF ( T_INFO%LSDYN ) THEN
         CALL SET_SELECTED_VEL_ZERO(T_INFO, DYN%VEL,LATT_CUR)
      ELSE IF (MDALGO_GL/=3 .AND. DYN%IBRION/=44 .AND. DYN%IBRION/=40) THEN
         CALL SYMVEL_WARNING( T_INFO%NIONS, T_INFO%NTYP, T_INFO%ITYP, &
         T_INFO%POMASS, DYN%VEL, IO%IU6, IO%IU0 )
      ENDIF
#else
      IF ( T_INFO%LSDYN ) THEN
         CALL SET_SELECTED_VEL_ZERO(T_INFO, DYN%VEL,LATT_CUR)
      ELSE
         !tb beg
         !c do not remove CoM motion if iDM or DVV is to be used
         IF (DYN%IBRION/=44 .AND. DYN%IBRION/=40) THEN
           CALL SYMVEL_WARNING( T_INFO%NIONS, T_INFO%NTYP, T_INFO%ITYP, &
             T_INFO%POMASS, DYN%VEL, IO%IU6, IO%IU0 )
         ENDIF
         !tb end
      ENDIF
#endif

      IF (.NOT.LSKIP_ML_ISTART3) THEN
         CALL NEAREST_NEIGHBOUR(IO%IU6, IO%IU0, T_INFO, LATT_CUR, P%RWIGS)
      ENDIF
!-----------------------------------------------------------------------
!  initialize the symmetry stuff
!-----------------------------------------------------------------------
      ALLOCATE(SYMM%ROTMAP(NIOND,1,1), &
               SYMM%TAU(NIOND,3), &
     &         SYMM%TAUROT(NIOND,3),SYMM%WRKROT(3*(NIOND+2)), &
     &         SYMM%PTRANS(NIOND+2,3),SYMM%INDROT(NIOND+2))
      IF (INFO%ISPIN==2) THEN
         ALLOCATE(SYMM%MAGROT(48,NIOND))
      ELSE
         ALLOCATE(SYMM%MAGROT(1,1))
      ENDIF
      ! break symmetry parallel to IGPAR
      IF (LBERRY) THEN
         LATT_CUR%A(:,IGPAR)=LATT_CUR%A(:,IGPAR)*(1+TINY*10)
         CALL LATTIC(LATT_CUR)
      ENDIF
! Rotate the initial magnetic moments to counter the clockwise
! rotation brought on by the spin spiral
      IF (WDES%LSPIRAL) CALL COUNTER_SPIRAL(WDES%QSPIRAL,T_INFO%NIONS,T_INFO%POSION,T_INFO%ATOMOM)

      NCDIJ=INFO%ISPIN
      IF (WDES%LNONCOLLINEAR) NCDIJ=4
! In case we would like to exploit symmetry (to some degree)
      IF (SYMM%ISYM>0) THEN
#ifdef debug
! First we check the sensitivity of the symmetry detection to the tolerance parameter SYMPREC
         CALL CHKTOL(LATT_CUR%A,DYN%POSION,DYN%VEL,DYN%EFOR,T_INFO%LSFOR, &
                     T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,NIOND, &
                     SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP, &
                     SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                     SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6,IO%NWRITE)
! then we do the final analysis
#endif
          CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,DYN%EFOR,T_INFO%LSFOR, &
                      T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,NIOND, &
                      SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP, &
                      SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                      SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
         IF (IO%NWRITE>=3) CALL WRTSYM(T_INFO%NIONS,NIOND,SYMM%PTRANS,SYMM%ROTMAP,SYMM%MAGROT,NCDIJ,IO%IU6)
      ELSE
! In case we do not want to exploit symmetry at all
         CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,NIOND,SYMM%PTRANS, &
              SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)
      END IF
      CALL SYMM_SET(SYMM)
      CALL PRIM_CELL%INIT_FROM_SYMM(LATT_CUR, DYN%POSION, SYMM)
      !CALL PRIM_CELL%WRITE_POSCAR(P(:)%ELEMENT,T_INFO,IO)
      CALL PRIM_CELL%WRITE_OUTCAR(IO)
      CALL PRIM_CELL%WRITE_XML(IO)
      IF (CHECKPOINT_FD_SETTINGS%DOES_FINISH()) THEN
         CALL POLAR_DATA%CHECK(PRIM_CELL)
         CALL ELPH_DERIV_SETTINGS%CHECK(DYN%IBRION, DYN%NFREE, PRIM_CELL)
      ENDIF

      ! Access checkpoint file
      IF (CHECKPOINT_FD_SETTINGS%DO_OVERWRITE_CHECKPOINT()) THEN
         CALL FINDIFF_CHECKPOINT%INIT('vaspcheckfd.h5', 'w', COMM)
      ELSE IF (CHECKPOINT_FD_SETTINGS%DO_WRITE_CHECKPOINT()) THEN
         CALL FINDIFF_CHECKPOINT%INIT('vaspcheckfd.h5', 'a', COMM)
      ELSE
         CALL FINDIFF_CHECKPOINT%INIT_DISABLED()
      ENDIF
      IF (FINDIFF_CHECKPOINT%SHOULD_READ()) CALL vtutor%write(isAlert, CheckpointContinue)
      IF (CHECKPOINT_FD_SETTINGS%DO_SINGLE()) CALL vtutor%write(isAlert, CheckpointSingle)

      IF (DYN%IBRION==5 .OR. DYN%IBRION==6 ) THEN
         DYN%NSW=4*(3*T_INFO%NIONS+9)+1
         DYN%KBLOCK=DYN%NSW
      ENDIF

      ANY_MLFF_AMONG_IMAGES = ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD
      CALLMPI( M_or(COMM_WORLD, ANY_MLFF_AMONG_IMAGES, 1) )

!-----be carefull about division by 0
      IF (ANY_MLFF_AMONG_IMAGES .AND. (LBLOCK_HELP /= 1)) THEN
         IF (DYN%NBLOCK/=1) CALL vtutor%alert("NBLOCK is set to ML_OUTBLOCK ( = "//str(LBLOCK_HELP)//" ).")
         DYN%NBLOCK=LBLOCK_HELP
      ELSE
         DYN%NBLOCK=MAX(1,DYN%NBLOCK)
      ENDIF
      DYN%KBLOCK=MAX(1,DYN%KBLOCK)
      ! Large numbers of NBLOCK, KBLOCK, NSW and ML_OUTBLOCK can lead to an integer overflow.
      ! Converting one factor + automatic type coercion in mixed expressions avoids this.
      IF (DYN%NSW < REAL(DYN%KBLOCK,q) * DYN%NBLOCK) DYN%KBLOCK = 1
      IF (DYN%NSW < REAL(DYN%KBLOCK,q) * DYN%NBLOCK) DYN%NBLOCK = MAX(DYN%NSW,1)
 
      DYN%NSW=INT(DYN%NSW/DYN%NBLOCK/DYN%KBLOCK)*DYN%NBLOCK*DYN%KBLOCK

      ! First initialization of LDO_AB_INITIO
      LDO_AB_INITIO = .TRUE.
      LDO_AB_INITIO_NEW = .TRUE.

#ifndef NO_MACHINE_LEARNING
      IF (INFO%PLUGIN%FORCE_FIELD) THEN
         IF (INFO%PLUGIN%ML_MODE .NE. "delta") THEN
            LDO_AB_INITIO = .FALSE.
            LDO_AB_INITIO_NEW = .FALSE.
         ENDIF
      END IF
#endif

#ifdef ML_AVAILABLE
! For safety NSW is set to 1 here if any of the ML_MODEs is set that
! requires NSW=1
      IF (ML_LMLFF) THEN
         DO I = 1, ML_TOTNUM_INSTANCES
            IF (ML_SUPER_HANDLE_MAIN(I)%FF%MODE.EQ."refit" .OR. &
                ML_SUPER_HANDLE_MAIN(I)%FF%MODE.EQ."REFIT" .OR. &
                ML_SUPER_HANDLE_MAIN(I)%FF%MODE.EQ."refitbayesian" .OR. &
                ML_SUPER_HANDLE_MAIN(I)%FF%MODE.EQ."REFITBAYESIAN" .OR. &
                ML_SUPER_HANDLE_MAIN(I)%FF%MODE.EQ."select" .OR. &
                ML_SUPER_HANDLE_MAIN(I)%FF%MODE.EQ."SELECT") THEN
               DYN%NSW=1
            ENDIF
         ENDDO
      ENDIF

! If ML_LMLFF=.TRUE., initialize the machine-learning variables.
      IF(ML_LMLFF) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*) 'Machine learning feature enabled'
         IF (IO%IU0>=0) WRITE(TIU0,*) 'Passing communicator to machine learning'
! set here the global communicator for all machine learning instances
         CALL MACHINE_LEARNING_SET_COMM_WORLD_GLOBAL(COMM)
#ifdef use_shmem
! Set here first the communicators for intra and inter band
! communication needed for shared memory
#endif
         IF (ML_SUPER_HANDLE_MAIN(1)%FF%LIB) THEN
#ifdef libvaspml
            CALL VASPML%init(ML_COMM_WORLD_GLOBAL, &
                             T_INFO%POSION, &
                             LATT_CUR%A, &
                             LATT_CUR%B, &
                             T_INFO%NIONS, &
                             T_INFO%NTYP, &
                             T_INFO%ITYP, &
                             T_INFO%TYPE)
            IF (LVCAIMAGES) CALL VASPML%post_sync_vcaimages(LBLOCK_HELP)
#endif
         ELSE
            DO I = 1, ML_TOTNUM_INSTANCES
! Initialize necessary variables and arrays for machine learning
               CALL MACHINE_LEARNING_INIT(ML_SUPER_HANDLE_MAIN(I), DIR_APP, DIR_LEN, &
                                          DYN, INFO, LATT_CUR, T_INFO, IO, COMM%MPI_COMM)
            ENDDO
         ENDIF
! If ISTART_FF=2 (this automatically means ML_FF_LMLONLY=.FALSE.), the program does not execute any ab initio calculations.
! In this case, LDO_AB_INITIO is set as .FALSE., here, too.
! For the moment only one instance of machine learning is run, so we set
! ML_FF_LMLONLY to LMLONLY of that instance
         IF(ML_FF_LMLONLY .AND. (.NOT.LDELTADFT)) THEN
            LDO_AB_INITIO=.FALSE.
            LDO_AB_INITIO_NEW=.FALSE.
         ELSE
            LDO_AB_INITIO=.TRUE.
            LDO_AB_INITIO_NEW=.TRUE.
         ENDIF
! Also set LDO_AB_INITIO to .FALSE. if we refit
         DO I = 1, ML_TOTNUM_INSTANCES
            IF (ML_SUPER_HANDLE_MAIN(I)%FF%ISTART.EQ.4) THEN
               LDO_AB_INITIO=.FALSE.
               LDO_AB_INITIO_NEW=.FALSE.
               EXIT
            ENDIF
         ENDDO
! Check if any of the instances has an option to create new ML_ABN file from old one (ML_FF_ISTART=3)
! If yes then no ab initio calculations have to be done
         DO I = 1, ML_TOTNUM_INSTANCES
            IF (ML_SUPER_HANDLE_MAIN(I)%FF%ISTART.EQ.3) THEN
               LDO_AB_INITIO=.FALSE.
               LDO_AB_INITIO_NEW=.FALSE.
               EXIT
            ENDIF
         ENDDO
      ENDIF
#endif

      NSW_ACT = DYN%NSW

!=======================================================================
!  Mini driver to compute the phonon dispersion and quit
!=======================================================================
      IF (PHON_SETTINGS%READ_FORCE_CONSTANTS) THEN
        IF(IO%IU0>=0)  WRITE(IO%IU0,*) 'Start phonon computation'
        ALLOCATE(FORCE_CONST(3*T_INFO%NIONS,3*T_INFO%NIONS))
#ifdef VASP_HDF5
        IF (.NOT.HDF5_FOUND) THEN
            CALL vtutor%error('Could not find ' // VASPIN // ' file from which to read the force constants')
        ENDIF
        CALL VH5_READ_FORCE_CONSTANTS(IH5INFILEID,FORCE_CONST)
#endif
        IF (PHON_SETTINGS%DO_INIT()) THEN
          CALL PHON_INIT(PRIM_CELL, -FORCE_CONST, IO)
        ENDIF

        IF (PHON_SETTINGS%DO_POLAR) THEN
          CALL PHON_INIT_POLAR(PHON_SETTINGS, PRIM_CELL, POLAR_DATA)
        ENDIF

        CALL PHONON_DRIVER(PRIM_CELL, T_INFO, PHON_SETTINGS, IO)

        IF(IO%IU0>=0)  WRITE(IO%IU0,*) 'Finished phonon computation'
        GOTO 5100
      ENDIF

!=======================================================================
!  Read UNIT=14: KPOINTS
!  number of k-points and k-points in reciprocal lattice
!=======================================================================
      IF(IO%IU6>=0)  WRITE(TIU6,*)
      ! use full k-point grid if finite differences are used or
      ! linear response is applied
      IF (DYN%IBRION==5 .OR. DYN%IBRION==6 .OR. DYN%IBRION==7.OR. DYN%IBRION==8 &
          .OR. LEPSILON .OR. LVEL .OR. LMAGBLOCH  &
          .OR. LCHIMAG .OR. LTIME_EVOLUTION) CALL USE_FULL_KPOINTS
      ! apply preferentially time inversion symmetry to generate orbitals at -k
      IF (WDES%LORBITALREAL) CALL USE_TIME_INVERSION

      IF (LBERRY) THEN
         CALL RD_KPOINTS_BERRY(KPOINTS,NPPSTR,IGPAR, &
        &   LATT_CUR, &
        &   SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
        &   IO%IU6,IO%IU0)
          IF (LBERRY) THEN
            LATT_CUR%A(:,IGPAR)=LATT_CUR%A(:,IGPAR)/(1+TINY*10)
            CALL LATTIC(LATT_CUR)
         ENDIF
      ELSE
#ifdef oldsym
         CALL SETUP_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND. &
            .NOT.WDES%LSORBIT.AND. &
            .NOT.WDES%LSPIRAL, &
            SYMM%ISYM<0,IO%IU6,IO%IU0)

         CALL SETUP_FULL_KPOINTS(KPOINTS,LATT_CUR,T_INFO%NIOND, &
            SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM, &
            SYMM%ISYM>=0.AND. &
            .NOT.WDES%LSORBIT.AND. &
            .NOT.WDES%LSPIRAL, &
            IO%IU6,IO%IU0,LSYMGRAD)
#else
         CALL SETUP_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            SYMM%ISYM<0,IO%IU6,IO%IU0)

         CALL SETUP_FULL_KPOINTS(KPOINTS,LATT_CUR,T_INFO%NIOND, &
            SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            IO%IU6,IO%IU0, LSYMGRAD)
#endif
      ENDIF
      CALL SETUP_ORIG_KPOINTS

!=======================================================================
!  at this point we have enough information to
!  create a param.inc file
!=======================================================================
      XCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
!
!  setup NGX, NGY, NGZ if required
!
! high precision do not allow for wrap around
      IF (INFO%SZPREC(1:7)=='singlen') THEN
        WFACT=3
      ELSE IF (INFO%SZPREC(1:1)=='h' .OR. INFO%SZPREC(1:1)=='a' .OR. INFO%SZPREC(1:1)=='s' ) THEN
        WFACT=4
      ELSE
! medium-low precision allow for small wrap around
        WFACT=3
      ENDIF
      GRID%NGPTAR(1)=XCUTOF*WFACT+0.5_q
      GRID%NGPTAR(2)=YCUTOF*WFACT+0.5_q
      GRID%NGPTAR(3)=ZCUTOF*WFACT+0.5_q
      IF (NGX /= -1)   GRID%NGPTAR(1)=  NGX
      IF (NGY /= -1)   GRID%NGPTAR(2)=  NGY
      IF (NGZ /= -1)   GRID%NGPTAR(3)=  NGZ
      CALL FFTCHK(GRID%NGPTAR)
      IF ( HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%ACTIVE .and. &
           HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%COARSEN_BEFORE_PAD ) THEN
           !
           ! Coarsened grids must be divisible by the padding factor
           !
          CALL CHANGE_GRID_FOR_KERNEL_TRUNCATION(GRID, HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE)
      ENDIF
!
!  setup NGXC, NGYC, NGZC if required
!
      IF (INFO%LOVERL) THEN
        IF (INFO%ENAUG==0) INFO%ENAUG=INFO%ENMAX*1.5_q
        IF (INFO%SZPREC(1:1)=='h') THEN
          WFACT=16._q/3._q
        ELSE IF (INFO%SZPREC(1:1)=='l') THEN
          WFACT=3
        ELSE
          WFACT=4
        ENDIF
        XCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
        YCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
        ZCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
        GRIDC%NGPTAR(1)=XCUTOF*WFACT
        GRIDC%NGPTAR(2)=YCUTOF*WFACT
        GRIDC%NGPTAR(3)=ZCUTOF*WFACT
        ! prec Single no double grid technique
        IF (INFO%SZPREC(1:1)=='s') THEN
           GRIDC%NGPTAR(1)=GRID%NGPTAR(1)
           GRIDC%NGPTAR(2)=GRID%NGPTAR(2)
           GRIDC%NGPTAR(3)=GRID%NGPTAR(3)
        ELSE IF (INFO%SZPREC(1:1)=='a' .OR. INFO%SZPREC(1:1)=='n') THEN
           GRIDC%NGPTAR(1)=GRID%NGPTAR(1)*2
           GRIDC%NGPTAR(2)=GRID%NGPTAR(2)*2
           GRIDC%NGPTAR(3)=GRID%NGPTAR(3)*2
        ENDIF
        IF (NGXC /= -1)  GRIDC%NGPTAR(1)=NGXC
        IF (NGYC /= -1)  GRIDC%NGPTAR(2)=NGYC
        IF (NGZC /= -1)  GRIDC%NGPTAR(3)=NGZC
        CALL FFTCHK(GRIDC%NGPTAR)
      ELSE
        GRIDC%NGPTAR(1)= 1
        GRIDC%NGPTAR(2)= 1
        GRIDC%NGPTAR(3)= 1
      ENDIF

      GRIDC%NGPTAR(1)=MAX(GRIDC%NGPTAR(1),GRID%NGPTAR(1))
      GRIDC%NGPTAR(2)=MAX(GRIDC%NGPTAR(2),GRID%NGPTAR(2))
      GRIDC%NGPTAR(3)=MAX(GRIDC%NGPTAR(3),GRID%NGPTAR(3))
      GRIDUS%NGPTAR=GRIDC%NGPTAR
      IF (LADDGRID) GRIDUS%NGPTAR=GRIDC%NGPTAR*2

      NGX = GRID %NGPTAR(1); NGY = GRID %NGPTAR(2); NGZ = GRID %NGPTAR(3)
      NGXC= GRIDC%NGPTAR(1); NGYC= GRIDC%NGPTAR(2); NGZC= GRIDC%NGPTAR(3)

      ! NBANDS was not read from INCAR, 
      ! so a default value must be selected for NBANDS 
      IF (NBANDS == -1 ) THEN
         LNBANDS=.FALSE.
         NBANDS=GET_NBANDS_SERIAL(WDES%LNONCOLLINEAR,INFO%ISPIN,INFO%NELECT,T_INFO)
         NBANDS=((NBANDS+NPAR-1)/NPAR)*NPAR
      ELSE IF (NBANDS == -2 ) THEN
         CALL LPRJ_FROM_PSEUDO(T_INFO,P,WDES%LNONCOLLINEAR,LPRJ_functions,NBANDS)
         NBANDS=MAX(INT(NBANDS*1.2),NBANDS+4)
      ENDIF

      IF (NBANDS/=((NBANDS+NPAR-1)/NPAR)*NPAR) THEN
         ITUT(1)=NBANDS
         ITUT(2)=((NBANDS+NPAR-1)/NPAR)*NPAR
         CALL vtutor%write(isAlert, NBANDSchanged, newArgument(ival = ITUT))
      ENDIF

      NBANDS=((NBANDS+NPAR-1)/NPAR)*NPAR

      IF (INFO%EBREAK == -1) INFO%EBREAK=0.25_q*MIN(INFO%EDIFF,ABS(DYN%EDIFFG)/10)/NBANDS

      IF(INFO%TURBO==0)THEN
         IF (((.NOT. WDES%LNONCOLLINEAR).AND.  INFO%NELECT>REAL(NBANDS*2,KIND=q)).OR. &
                 ((WDES%LNONCOLLINEAR).AND.((INFO%NELECT*2)>REAL(NBANDS*2,KIND=q)))) THEN
            ITUT(1)=INFO%NELECT ; ITUT(2)=NBANDS ; CTUT(1) = "NBANDS"
            CALL vtutor%write(isError, NumberOfElectrons, newArgument(ival = ITUT, cval = CTUT))
         ENDIF
      ELSE
         IF (( (.NOT. WDES%LNONCOLLINEAR).AND.  INFO%NELECT>REAL(NBANDS*2,KIND=q)).OR. &
                    ((WDES%LNONCOLLINEAR).AND.((INFO%NELECT*2)>REAL(NBANDS*2,KIND=q)))) THEN
            IF(KPOINTS%EFERMI==HUGE(1.0_q)) THEN
               ITUT(1)=INFO%NELECT ; ITUT(2)=NBANDS ; CTUT(1) = "NBANDS"
               CALL vtutor%write(isError, NumberOfElectrons, newArgument(ival = ITUT, cval = CTUT))
            ENDIF
         ENDIF
      ENDIF

      NRPLWV=4*PI*SQRT(INFO%ENMAX /RYTOEV)**3/3* &
     &     LATT_CUR%OMEGA/AUTOA**3/(2*PI)**3*1.1_q+50
#ifdef gammareal
      NRPLWV=NRPLWV/2
#endif
      WDES%NRPLWV=NRPLWV
      PSRMX=0
      PSDMX=0
      DO NT=1,T_INFO%NTYP
        PSRMX=MAX(PSRMX,P(NT)%PSRMAX)
        PSDMX=MAX(PSDMX,P(NT)%PSDMAX)
      ENDDO
      IF (INFO%LREAL) THEN
       IRMAX=4*PI*PSRMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRID%NGPTAR(1)*GRID%NGPTAR(2)*GRID%NGPTAR(3)))+50
      ELSE
       IRMAX=1
      ENDIF
      IRDMAX=1
      IF (INFO%LOVERL) THEN
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDC%NGPTAR(1)*GRIDC%NGPTAR(2)*GRIDC%NGPTAR(3)))+200
      ENDIF
#ifdef usgrid
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDUS%NGPTAR(1)*GRIDUS%NGPTAR(2)*GRIDUS%NGPTAR(3)))+200
#endif

      NPLWV =NGX *NGY *NGZ;
      MPLWV =NGX *NGY *NGZ
      NPLWVC=NGXC*NGYC*NGZC;
      MPLWVC=xm(NGXC)*NGYC*zm(NGZC)
!=======================================================================
!  set the basic quantities in WDES
!  and set the grids
!=======================================================================

      WDES%ENMAX =INFO%ENMAX

      WDES%NB_PAR=NPAR
      WDES%NB_TOT=NBANDS
      WDES%NBANDS=NBANDS/NPAR
#ifdef MPI
      WDES%NB_LOW=COMM_INTER%NODE_ME
#else
      WDES%NB_LOW=1
#endif
      CALL INIT_KPOINT_WDES(WDES, KPOINTS )
      WDES%ISPIN =INFO%ISPIN
      WDES%COMM  =>COMM
      WDES%COMM_INB    =>COMM_INB
      WDES%COMM_INTER  =>COMM_INTER
      WDES%COMM_KIN    =>COMM_KIN
      WDES%COMM_KINTER =>COMM_KINTER
      NULLIFY( WDES%COMM_SHMEM )
#ifdef use_shmem
      WDES%COMM_SHMEM  =>COMM_SHMEM
      WDES%COMM_intra_node =>COMM_intra_node
      WDES%COMM_inter_node =>COMM_inter_node
#endif

      CALL  SET_FULL_KPOINTS(WDES%NKPTS_FOR_GEN_LAYOUT,WDES%VKPT)
      CALL  SET_FOCK_KPOINTS(WDES%NKDIM)

      IF (WDES%LNONCOLLINEAR) then
         WDES%NRSPINORS = 2
         INFO%RSPIN = 1
      ELSE
         WDES%NRSPINORS = 1
      ENDIF
      WDES%RSPIN = INFO%RSPIN

      CALL WDES_SET_NPRO(WDES,T_INFO,P,INFO%LOVERL)
!
! set up the descriptor for the initial wavefunctions
! (read from file)
      LATT_INI=LATT_CUR
! get header from WAVECAR file (LATT_INI is important)
! also set INFO%ISTART
      IF (INFO%ISTART > 0) THEN
#ifdef VASP_HDF5
        IF (WAVECAR_FOUND) THEN
#endif
        CALL INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI,INFO%ISTART, IO%IU0)
#ifdef VASP_HDF5
        ELSE IF (VASPWAVE_FOUND) THEN
           IH5ERR = VH5_READ_WAVEFUNCTION_HEADER(IH5WAVEFILEID, WDES, LATT_INI, LATT_CUR, ENMAXI, INFO%ISTART, IO%IU0)
           ! if wave dataset is not present, initialize wavefunction later on with random numbers
           IF( IH5ERR /=0 ) INFO%ISTART=0 
        ELSE
          ! When neither a WAVECAR nor a vaspwave.h5 were found set INFO%ISTART=0
          ! (when running without HDF5 support this is done by INWAV_HEAD).
          INFO%ISTART=0
        ENDIF
#endif
        IF (INFO%ISTART == 0 .AND. INFO%ICHARG == 0) INFO%ICHARG=2
      ENDIF

      CALL INIT_SCALAAWARE( WDES%NB_TOT, NRPLWV, WDES%COMM_KIN )

!=======================================================================
!  Write all important information
!=======================================================================
      IF (DYN%IBRION==5 .OR. DYN%IBRION==6 ) THEN
         IF (DYN%NFREE /= 1 .AND. DYN%NFREE /= 2 .AND. DYN%NFREE /= 4)  DYN%NFREE =2
      ENDIF

      IF (IO%IU6>=0) THEN

      WRITE(TIU6,130)
      WRITE(TIU6,7205) KPOINTS%NKPTS,WDES%NKDIM,WDES%NB_TOT,NEDOS, &
     &              T_INFO%NIONS,LDIM,LMDIM, &
     &              NPLWV,IRMAX,IRDMAX, &
     &              NGX,NGY,NGZ, &
     &              NGXC,NGYC,NGZC,GRIDUS%NGPTAR,T_INFO%NITYP
      IF (WDES%NBANDSLOW/=-1 .OR. WDES%NBANDSHIGH/=-1) THEN
        WRITE(TIU6,"('   lowest and highest band optimized NBANDSLOW= ',I7,' NBANDSHIGH= ',I7)" ) WDES%NBANDSLOW, WDES%NBANDSHIGH
      ENDIF

      XAU= (NGX*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YAU= (NGY*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZAU= (NGZ*PI/(LATT_CUR%ANORM(3)/AUTOA))
      WRITE(TIU6,7211) XAU,YAU,ZAU
      XAU= (NGXC*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YAU= (NGYC*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZAU= (NGZC*PI/(LATT_CUR%ANORM(3)/AUTOA))
      WRITE(TIU6,7212) XAU,YAU,ZAU

      ENDIF

 7211 FORMAT('   NGX,Y,Z   is equivalent  to a cutoff of ', &
     &           F6.2,',',F6.2,',',F6.2,' a.u.')
 7212 FORMAT('   NGXF,Y,Z  is equivalent  to a cutoff of ', &
     &           F6.2,',',F6.2,',',F6.2,' a.u.'/)

      XCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
! high precision do not allow for wrap around
      IF (INFO%SZPREC(1:7)=='singlen') THEN
        WFACT=3
      ELSEIF (INFO%SZPREC(1:1)=='s'.OR.INFO%SZPREC(1:1)=='h'.OR.INFO%SZPREC(1:1)=='a') THEN
        WFACT=4
      ELSE
! medium-low precision allow for small wrap around
        WFACT=3
      ENDIF
      ITUT(1)=XCUTOF*WFACT+0.5_q
      ITUT(2)=YCUTOF*WFACT+0.5_q
      ITUT(3)=ZCUTOF*WFACT+0.5_q
      CALL FFTCHK(ITUT(1:3))
      IF (IO%IU6>=0 .AND. (ITUT(1)/=NGX .OR. ITUT(2)/=NGY .OR. ITUT(3)/=NGZ)) &
        WRITE(TIU6,72111) ITUT(1),ITUT(2),ITUT(3)

72111 FORMAT(' Based on PREC, I would recommend the following setting:'/ &
     &       '   dimension x,y,z NGX = ',I5,' NGY =',I5,' NGZ =',I5 /)

      IF (.NOT. ML_FF_LMLONLY .AND. (NGX<ITUT(1) .OR. NGY<ITUT(2) .OR. NGZ<ITUT(3))) THEN
         CALL vtutor%write(isAlert, FFTGridIsNotSufficient, newArgument(ival = ITUT))
      ENDIF


      AOMEGA=LATT_CUR%OMEGA/T_INFO%NIONS
      QF=(3._q*PI*PI*INFO%NELECT/(LATT_CUR%OMEGA))**(1._q/3._q)*AUTOA

! chose the mass so that the typical Nose frequency is 40 timesteps
!-----just believe this (or look out  for all this factors in  STEP)
      IF (DYN%SMASS==0)  DYN%SMASS= &
         ((40._q*DYN%POTIM*1E-15_q/2._q/PI/LATT_CUR%ANORM(1))**2)* &
         2.E20_q*BOLKEV*EVTOJ/AMTOKG*NDEGREES_OF_FREEDOM*MAX(DYN%TEBEG,DYN%TEEND)
!      IF (DYN%SMASS<0)  DYN%SMASS= &
!         ((ABS(DYN%SMASS)*DYN%POTIM*1E-15_q/2._q/PI/LATT_CUR%ANORM(1))**2)* &
!         2.E20_q*BOLKEV*EVTOJ/AMTOKG*NDEGREES_OF_FREEDOM*MAX(DYN%TEBEG,DYN%TEEND)

      SQQ=  DYN%SMASS*(AMTOKG/EVTOJ)*(1E-10_q*LATT_CUR%ANORM(1))**2
      SQQAU=SQQ/RYTOEV
      IF (DYN%SMASS>0) THEN
!        WOSZI= SQRT(2**DYN%TEMP*NDEGREES_OF_FREEDOM/SQQ) ! gK: who put this in?
        WOSZI= SQRT(2*BOLKEV*DYN%TEMP*NDEGREES_OF_FREEDOM/SQQ)
      ELSE
        WOSZI=1E-30_q
      ENDIF
!-----initial temperature
      CALL EKINC(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,DYN%POTIM,LATT_CUR%A,DYN%VEL)
      TEIN = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM
      IF (IO%IU6>=0) THEN

      WRITE(TIU6,7210) INFO%SZNAM1,T_INFO%SZNAM2

      WRITE(TIU6,7206) IO%NWRITE,INFO%SZPREC,INFO%ISTART,INFO%ICHARG,WDES%ISPIN,WDES%LNONCOLLINEAR, &
     &      WDES%LSORBIT,INFO%INIWAV,INFO%LASPH, &
     &      INFO%ENMAX,INFO%ENMAX/RYTOEV,SQRT(INFO%ENMAX/RYTOEV), &
     &      XCUTOF,YCUTOF,ZCUTOF,INFO%ENINI, &
     &      INFO%ENAUG, &
     &      INFO%NELM,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,INFO%LREAL,INFO%NLSPLINE,LCOMPAT,GGA_COMPAT, &
     &      LMAX_CALC,SET_LMAX_MIX_TO,XC_DATA%LFCI, &
     &      T_INFO%ROPT

      WRITE(TIU6,7226) EWALD_PARAMS%MAX_ARGUMENT

      WRITE(TIU6,7204) &
     &      DYN%EDIFFG,DYN%NSW,DYN%NBLOCK,DYN%KBLOCK, &
     &      DYN%IBRION,DYN%NFREE,DYN%ISIF,PRED%IWAVPR,PRED%MAXDIS,SYMM%ISYM,INFO%LCORR

      TMP=0
      IF (DYN%POTIM>0) TMP=1/(WOSZI*(DYN%POTIM*1E-15_q)/2._q/PI)

      WRITE(TIU6,7207) &
     &      DYN%POTIM,TEIN,DYN%TEBEG,DYN%TEEND, &
     &      DYN%SMASS,WOSZI,TMP,SQQAU,SCALEE, &
     &      PACO%NPACO,PACO%APACO,DYN%PSTRESS

      WRITE(TIU6,7215) (T_INFO%POMASS(NI),NI=1,T_INFO%NTYP)
      RTUT(1:T_INFO%NTYP)=P(1:T_INFO%NTYP)%ZVALF ! work around IBM bug
      WRITE(TIU6,7216) (RTUT(NI),NI=1,T_INFO%NTYP)
      WRITE(TIU6,7203) (T_INFO%RWIGS(NI),NI=1,T_INFO%NTYP)
      WRITE(TIU6,72031) (T_INFO%VCA(NI),NI=1,T_INFO%NTYP)

      IF (KPOINTS%EFERMI == HUGE(1._q)) THEN
          ! In this case the EFERMI tag is not set and so the Fermi energy is written out as zero
          WRITE(TIU6,7208) &
         &      INFO%NELECT,INFO%NUP_DOWN, &
         &      KPOINTS%EMIN,KPOINTS%EMAX, &
         &      0._q,KPOINTS%FERMI_METHOD, &
         &      KPOINTS%ISMEAR,KPOINTS%SIGMA,KPOINTS%EFERMI_NEDOS
      ELSE
          WRITE(TIU6,7208) &
         &      INFO%NELECT,INFO%NUP_DOWN, &
         &      KPOINTS%EMIN,KPOINTS%EMAX, &
         &      KPOINTS%EFERMI,KPOINTS%FERMI_METHOD, &
         &      KPOINTS%ISMEAR,KPOINTS%SIGMA,KPOINTS%EFERMI_NEDOS
      END IF

      WRITE(TIU6,7209) &
     &      INFO%IALGO,INFO%LDIAG,INFO%LSUBROT,INFO%LROBUST,INFO%ISEARCH, &
     &      INFO%TURBO,INFO%IRESTART,INFO%NREBOOT,INFO%NMIN,INFO%EREF,&
     &      MIX%IMIX,MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
     &      MIX%WC,MIX%INIMIX,MIX%MIXPRE,MIX%MAXMIX, &
     &      INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,INFO%TIME, &
     &      AOMEGA,AOMEGA/(AUTOA)**3, &
     &      QF,QF/AUTOA,QF**2*RYTOEV,QF**2,SQRT(QF*4/PI)/AUTOA
      WRITE(TIU6,*)
      WRITE(TIU6,7224) IO%LH5, IO%LWAVE,IO%LDOWNSAMPLE,IO%WRT_DENSITY%LCHARG,&
            IO%WRT_DENSITY%LCHARGH5,IO%WRT_DENSITY%LTAU,IO%WRT_POTENTIAL%LVTOT,&
            IO%WRT_POTENTIAL%LVHAR,STR(IO%WRT_POTENTIAL),STR(IO%WRT_DENSITY),&
            IO%LELF,IO%LORBIT

      CALL WRITE_PLUGINS(TIU6, INFO%PLUGIN)
      WRITE(TIU6,*)

      CALL WRITE_EFIELD(TIU6)
      WRITE(TIU6,*)

      CALL WRITE_COULOMB_CUTOFF(TIU6, HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE)
      WRITE(TIU6,*)

      ENDIF

 7210 FORMAT( &
     &       ' SYSTEM =  ',A40/ &
     &       ' POSCAR =  ',A40/)

 7205 FORMAT(//' Dimension of arrays:'/ &
     &       '   k-points           NKPTS = ',I6, &
     &       '   k-points in BZ     NKDIM = ',I6, &
     &       '   number of bands    NBANDS= ',I6/ &
     &       '   number of dos      NEDOS = ',I6, &
     &       '   number of ions     NIONS = ',I6/ &
     &       '   non local maximal  LDIM  = ',I6, &
     &       '   non local SUM 2l+1 LMDIM = ',I6/ &
     &       '   total plane-waves  NPLWV = ',I6/ &
     &       '   max r-space proj   IRMAX = ',I6, &
     &       '   max aug-charges    IRDMAX= ',I6/ &
     &       '   dimension x,y,z NGX = ',I5,' NGY =',I5,' NGZ =',I5/ &
     &       '   dimension x,y,z NGXF= ',I5,' NGYF=',I5,' NGZF=',I5/ &
     &       '   support grid    NGXF= ',I5,' NGYF=',I5,' NGZF=',I5/ &
     &       '   ions per type =            ',999I4/)

 7206 FORMAT(' Startparameter for this run:'/ &
     &       '   NWRITE = ',I6,  '    write-flag & timer' / &
     &       '   PREC   = ',A6,  '    normal or accurate (medium, high low for compatibility)'/ &
     &       '   ISTART = ',I6,  '    job   : 0-new  1-cont  2-samecut'/ &
     &       '   ICHARG = ',I6,  '    charge: 1-file 2-atom 10-const'/ &
     &       '   ISPIN  = ',I6,  '    spin polarized calculation?'/ &
     &       '   LNONCOLLINEAR = ',L6, ' non collinear calculations'/ &
     &       '   LSORBIT = ',L6, '    spin-orbit coupling'/ &
     &       '   INIWAV = ',I6,  '    electr: 0-lowe 1-rand  2-diag'/ &
     &       '   LASPH  = ',L6,  '    aspherical Exc in radial PAW'/ &
     &       ' Electronic Relaxation 1'/ &
     &       '   ENCUT  = ', &
     &              F6.1,' eV ',F6.2,' Ry  ',F6.2,' a.u. ', &
     &              3F6.2,'*2*pi/ulx,y,z'/ &
     &       '   ENINI  = ',F6.1,'     initial cutoff'/ &
     &       '   ENAUG  = ',F6.1,' eV  augmentation charge cutoff'/ &
     &       '   NELM   = ',I6,  ';   NELMIN=',I3,'; NELMDL=',I3, &
     &         '     # of ELM steps '    / &
     &       '   EDIFF  = ',E7.1,'   stopping-criterion for ELM'/ &
     &       '   LREAL  = ',L6,  '    real-space projection'     / &
     &       '   NLSPLINE    = ',L1,'    spline interpolate recip. space projectors'/ &
     &       '   LCOMPAT= ',L6,  '    compatible to vasp.4.4'/&
     &       '   GGA_COMPAT  = ',L1,'    GGA compatible to vasp.4.4-vasp.4.6'/&
     &       '   LMAXPAW     = ',I4,' max onsite density'/&
     &       '   LMAXMIX     = ',I4,' max onsite mixed and CHGCAR'/&
     &       '   VOSKOWN= ',I6,  '    Vosko Wilk Nusair interpolation'/&
     &      ('   ROPT   = ',4F10.5))

 7226 FORMAT( &
     &      /' Ewald summation'/ &
     &       '   EWALD_CUTOFF = ',F6.2/)

 7204 FORMAT( &
     &       ' Ionic relaxation'/ &
     &       '   EDIFFG = ',E7.1,'   stopping-criterion for IOM'/ &
     &       '   NSW    = ',I6,  '    number of steps for IOM' / &
     &       '   NBLOCK = ',I6,  ';   KBLOCK = ',I6, &
     &         '    inner block; outer block '/ &
     &       '   IBRION = ',I6, &
     &         '    ionic relax: 0-MD 1-quasi-New 2-CG'/ &
     &       '   NFREE  = ',I6,  &
     &         '    steps in history (QN), initial steepest desc. (CG)'/ &
     &       '   ISIF   = ',I6,  '    stress and relaxation' / &
     &       '   IWAVPR = ',I6, &
     &         '    prediction:  0-non 1-charg 2-wave 3-comb' / &
     &       '   MAXDIS = ',F6.2, &
     &         '    maximum displacement; if above reset density' / &
     &       '   ISYM   = ',I6, &
     &         '    0-nonsym 1-usesym 2-fastsym' / &
     &       '   LCORR  = ',L6, &
     &         '    Harris-Foulkes like correction to forces' /)

 7207 FORMAT( &
     &       '   POTIM  =' ,F7.4,'    time-step for ionic-motion'/ &
     &       '   TEIN   = ',F6.1,'    initial temperature'       / &
     &       '   TEBEG  = ',F6.1,';   TEEND  =',F6.1, &
     &               ' temperature during run'/ &
     &       '   SMASS  = ',F6.2,'    Nose mass-parameter (am)'/ &
     &       '   estimated Nose-frequenzy (Omega)   = ',E9.2, &
     &           ' period in steps =',E9.2,' mass=',E12.3,'a.u.'/ &
     &       '   SCALEE = ',F6.4,'    scale energy and forces'       / &
     &       '   NPACO  = ',I6,  ';   APACO  = ',F4.1, &
     &       '  distance and # of slots for P.C.'  / &
     &       '   PSTRESS= ',F6.1,' pullay stress'/)

!    &       '   damping for Cell-Motion     SIDAMP = ',F6.2/
!    &       '   mass for Cell-Motion        SIMASS = ',F6.2/

 7215 FORMAT('  Mass of Ions in am'/ &
     &       ('   POMASS = ',8F6.2))
 7216 FORMAT('  Ionic Valenz'/ &
     &       ('   ZVAL   = ',8F6.2))
 7203 FORMAT('  Atomic Wigner-Seitz radii'/ &
     &       ('   RWIGS  = ',8F6.2))
72031 FORMAT('  virtual crystal weights '/ &
     &       ('   VCA    = ',8F6.2))

 7208 FORMAT( &
     &       '   NELECT = ',F12.4,  '    total number of electrons'/ &
     &       '   NUPDOWN= ',F12.4,  '    fix difference up-down'// &
     &       ' DOS related values:'/ &
     &       '   EMIN   = ',F6.2,';   EMAX   =',F6.2, &
     &       '  energy-range for DOS'/ &
     &       '   EFERMI = ',F6.2,';   METHOD = ',A,/ &
     &       '   ISMEAR =',I6,';   SIGMA  = ',F6.2, &
     &       '  broadening in eV -4-tet -1-fermi 0-gaus'/ &
     &       '   EFERMI_NEDOS =',I5,'; number of points for ismear=-14 and -15'/)

 7209 FORMAT( &
     &       ' Electronic relaxation 2 (details)'/  &
     &       '   IALGO  = ',I6,  '    algorithm'            / &
     &       '   LDIAG  = ',L6,  '    sub-space diagonalisation (order eigenvalues)' / &
     &       '   LSUBROT= ',L6,  '    optimize rotation matrix (better conditioning)' / &
     &       '   LROBUST= ',L6,  '    robust break criterion during scf (extrapolated energy error)' / &
     &       '   ISEARCH= ',I6,  '    line search for conjugate gradient (0, 1=use gradient)' / &
     &       '   TURBO    = ',I6,  '    0=normal 1=particle mesh'/ &
     &       '   IRESTART = ',I6,  '    0=no restart 2=restart with 2 vectors'/ &
     &       '   NREBOOT  = ',I6,  '    no. of reboots'/ &
     &       '   NMIN     = ',I6,  '    reboot dimension'/ &
     &       '   EREF     = ',F6.2,  '    reference energy to select bands'/ &
     &       '   IMIX   = ',I6,  '    mixing-type and parameters'/ &
     &       '     AMIX     = ',F6.2,';   BMIX     =',F6.2/ &
     &       '     AMIX_MAG = ',F6.2,';   BMIX_MAG =',F6.2/ &
     &       '     AMIN     = ',F6.2/ &
     &       '     WC   = ',F6.0,';   INIMIX=',I4,';  MIXPRE=',I4,';  MAXMIX=',I4// &
     &       ' Intra band minimization:'/ &
     &       '   WEIMIN = ',F6.4,'     energy-eigenvalue tresh-hold'/ &
     &       '   EBREAK = ',E9.2,'  absolut break condition' / &
     &       '   DEPER  = ',F6.2,'     relativ break condition  ' // &
     &       '   TIME   = ',F6.2,'     timestep for ELM'          // &
     &       '  volume/ion in A,a.u.               = ',F10.2,3X,F10.2/ &
     &       '  Fermi-wavevector in a.u.,A,eV,Ry     = ',4F10.6/ &
     &       '  Thomas-Fermi vector in A             = ',2F10.6/)


 7224 FORMAT( &
     &       ' Write flags'/  &
     &       '   LH5          = ',L6,  '    write vaspwave.h5 (instead of legacy files WAVECAR, CHGCAR, etc.)' / &
     &       '   LWAVE        = ',L6,  '    write Kohn-Sham orbitals (WAVECAR/vaspwave.h5)' / &
     &       '   LDOWNSAMPLE  = ',L6,  '    k-point downsampling of WAVECAR' / &
     &       '   LCHARG       = ',L6,  '    write restart densities (CHGCAR/vaspwave.h5)' / &
     &       '   LCHARGH5     = ',L6,  '    write charge density to vaspwave.h5' / &
     &       '   LTAU         = ',L6,  '    evaluate kinetic energy density (TAUCAR/vaspwave.h5)' / &
     &       '   LVTOT        = ',L6,  '    write LOCPOT, total local potential' / &
     &       '   LVHAR        = ',L6,  '    write LOCPOT, Hartree potential only' / &
     &       '   WRT_POTENTIAL= ',A,   '    ! write potential to vaspout.h5' / &
     &       '   WRT_DENSITY  = ',A,   '    ! write densities to vaspout.h5' / &
     &       '   LELF         = ',L6,  '    write electronic localiz. function (ELF)'/&
     &       '   LORBIT       = ',I6,  '    0 simple, 1 ext, 2 COOP (PROOUT), +10 PAW based schemes'//)

       CALL WRITE_CL_SHIFT(IO%IU6)
       IF (USELDApU()) CALL WRITE_LDApU(IO%IU6)

       CALL WRITE_FOCK(XC_DATA,IO%IU6)
       CALL VDW_WRITER(VDW_SET_GL)

       CALL WRITE_BERRY_PARA(IO%IU6,LBERRY,IGPAR,NPPSTR)
       CALL LR_WRITER(IO%IU6)
       CALL WAVE_INTERPOLATOR_WRITE_OUTCAR(WAVE_INTERP_SET,IO%IU6)
#ifdef gammareal
       !> the gamma-only version does not support KPOINTS_OPT 
       IF ( WAVE_INTERP_SET%LKPOINTS_OPT ) THEN
           CALL VTUTOR%ERROR( 'The gamma-only version does not support band structure calculations. &
             &Please use vasp_std or vasp_ncl in combination with KPOINTS_OPT. If you do not intent to &
             &run a band structure calculation, please delete KPOINTS_OPT from the working directory.')
       ENDIF
#endif
       CALL WRITE_ORBITALMAG(IO%IU6)
       CALL WRITE_RESPONSE(IO%IU6)
       CALL PEAD_WRITER(IO%IU6)
       CALL CHGFIT_WRITER(IO)
       CALL WRITE_RANDOM(IO%IU6,WDES%RANDOM_NUMBER_GENERATOR,WDES%PCG_SEED)
! solvation__
       CALL SOL_WRITER(IO)
! solvation__

#ifdef logflow
       CLOSE(19)
#endif
       CALL XML_TAG("parameters")
       CALL XML_WRITER( &
          NPAR, &
          INFO%SZNAM1,INFO%ISTART,INFO%IALGO,MIX%IMIX,MIX%MAXMIX,MIX%MREMOVE, &
          MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
          MIX%WC,MIX%INIMIX,MIX%MIXPRE,MIX%MIXFIRST,IO%LFOUND,INFO%LDIAG,INFO%LSUBROT,INFO%LROBUST,INFO%ISEARCH,INFO%LREAL,IO%LREALD,IO%LPDENS, &
          DYN%IBRION,INFO%ICHARG,INFO%INIWAV,INFO%NELM,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,DYN%EDIFFG, &
          DYN%NSW,DYN%ISIF,PRED%IWAVPR,PRED%MAXDIS,SYMM%ISYM,DYN%NBLOCK,DYN%KBLOCK,INFO%ENMAX,DYN%POTIM,DYN%TEBEG, &
          DYN%TEEND,DYN%NFREE, &
          PACO%NPACO,PACO%APACO,T_INFO%NTYP,NTYPD,DYN%SMASS,SCALEE,T_INFO%POMASS,T_INFO%DARWIN_V,T_INFO%DARWIN_R,  &
          T_INFO%RWIGS,INFO%NELECT,INFO%NUP_DOWN,INFO%TIME,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%EFERMI, &
          KPOINTS%ISMEAR,KPOINTS%SPACING,KPOINTS%LGAMMA,KPOINTS%LKBLOWUP,DYN%PSTRESS,INFO%NDAV, &
          KPOINTS%SIGMA,KPOINTS%LTET,INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,IO%NWRITE,INFO%LCORR, &
          IO%IDIOT,T_INFO%NIONS,T_INFO%NTYPP,IO%LMUSIC,IO%LOPTICS,STM, &
          INFO%ISPIN,T_INFO%ATOMOM,NIOND,IO%LWAVE,IO%LDOWNSAMPLE,IO%WRT_DENSITY%LCHARG,IO%WRT_POTENTIAL%LVTOT,IO%WRT_POTENTIAL%LVHAR,INFO%SZPREC, &
          INFO%ENAUG,IO%LORBIT,IO%LELF,T_INFO%ROPT,INFO%ENINI, &
          NGX,NGY,NGZ,NGXC,NGYC,NGZC,NBANDS,WDES%NBANDSLOW,WDES%NBANDSHIGH,NEDOS,NBLK,LATT_CUR, &
          LPLANE_WISE,LCOMPAT,LMAX_CALC,SET_LMAX_MIX_TO,WDES%NSIM,IO%WRT_DENSITY%LPARD,LPAW,LADDGRID, &
          WDES%LNONCOLLINEAR,WDES%LSORBIT,WDES%SAXIS, &
          WDES%LSPIRAL,WDES%LZEROZ,WDES%QSPIRAL, &
          INFO%LASPH,WDES%LORBITALREAL, &
          INFO%TURBO,INFO%IRESTART,INFO%NREBOOT,INFO%NMIN,INFO%EREF, &
          INFO%NLSPLINE,MDALGO_GL,PHON_NSTRUCT)

       CALL  XML_WRITE_GGA_COMPAT_MODE
       CALL  XML_WRITE_BERRY(LBERRY, IGPAR, NPPSTR)
       CALL  XML_WRITE_CL_SHIFT
       CALL  XML_WRITE_LDAU
       CALL  XML_WRITE_CONSTRAINED_M(T_INFO%NIONS)
       CALL  XML_WRITE_XC_FOCK(XC_DATA)
       CALL  XML_WRITE_LR
       CALL  XML_WRITE_ORBITALMAG
       CALL  XML_WRITE_RESPONSE
       CALL  XML_WRITE_CLASSICFIELDS
       CALL  XML_WRITE_WAVE_INTERPOLATOR(WAVE_INTERP_SET)
! solvation__
       CALL XML_WRITE_SOL
! solvation__

       CALL XML_CLOSE_TAG("parameters")
!=======================================================================
!  set some important flags and write out text information
!  DYN%IBRION        selects dynamic
!  INFO%LCORR =.TRUE. calculate Harris corrections to forces
!=======================================================================
       IF (MIX%AMIN>=0.1_q .AND. MAXVAL(LATT_CUR%ANORM(:))>50) THEN
          CALL vtutor%write(isAlert, LongLattice)
       ENDIF
!---- relaxation related information
      IF (DYN%IBRION==10) THEN
         INFO%NELMDL=ABS(INFO%NELM)
         INFO%LCORR=.TRUE.
         IF (DYN%POTIM <= 0.0001_q ) DYN%POTIM=1E-20_q
      ENDIF

      IF (IO%IU6>=0) THEN

      WRITE(TIU6,130)

      IF (DYN%IBRION == -1) THEN
        WRITE(TIU6,*)'Static calculation'
      ELSE IF (DYN%IBRION==0) THEN
        WRITE(TIU6,*)'molecular dynamics for ions'
        IF (DYN%SMASS>0) THEN
          WRITE(TIU6,*)'  using nose mass (canonical ensemble)'
        ELSE IF (DYN%SMASS==-3) THEN
          WRITE(TIU6,*)'  using a microcanonical ensemble'
        ELSE IF (DYN%SMASS==-1) THEN
          WRITE(TIU6,*)'  scaling velocities every NBLOCK steps'
        ELSE IF (DYN%SMASS==-2) THEN
          WRITE(TIU6,*)'  keeping initial velocities unchanged'
        ENDIF
      ELSE IF (DYN%IBRION==4) THEN
         WRITE(TIU6,*)'Beeman-Brooks-Karplus method for relaxation of ions'
      ELSE IF (DYN%IBRION==1) THEN
           WRITE(TIU6,*)'quasi-Newton-method for relaxation of ions'
      ELSE IF (DYN%IBRION==2) THEN
           WRITE(TIU6,*)'conjugate gradient relaxation of ions'
      ELSE IF (DYN%IBRION==3) THEN
              WRITE(TIU6,*)'quickmin algorithm: (dynamic with friction)'
      ELSE IF (DYN%IBRION==5) THEN
              WRITE(TIU6,*)'finite differences'
      ELSE IF (DYN%IBRION==6) THEN
              WRITE(TIU6,*)'finite differences with symmetry'
      ELSE IF (DYN%IBRION==7) THEN
              WRITE(TIU6,*)'linear response'
      ELSE IF (DYN%IBRION==8) THEN
              WRITE(TIU6,*)'linear response using symmetry'
      ELSE IF (DYN%IBRION==10) THEN
           WRITE(TIU6,*)'relaxation of ions and charge simultaneously'
      ELSE IF (DYN%IBRION==11) THEN
           WRITE(TIU6,*)'interactive mode, write forces and read positions'
      ELSE IF (DYN%IBRION==12) THEN
           WRITE(TIU6,*)'interactive mode, call external program to update positions'
      !tb beg
      ELSE IF (DYN%IBRION==44) THEN
          WRITE(TIU6,*)'improved dimer method for transition state relaxation'
      !tb end
      ENDIF

      IF (DYN%IBRION/=-1 .AND. T_INFO%LSDYN) THEN
        WRITE(TIU6,*)'  using selective dynamics as specified on POSCAR'
        IF (.NOT.T_INFO%LDIRCO) THEN
          WRITE(TIU6,*)'  WARNING: If single coordinates had been '// &
     &                'selected the selection of coordinates'
          WRITE(TIU6,*)'           is made according to the '// &
     &                'corresponding   d i r e c t   coordinates!'
          WRITE(TIU6,*)'           Don''t support selection of '// &
     &                'single cartesian coordinates -- sorry ... !'
        ENDIF
      ENDIF
      ENDIF

      IF (INFO%ICHARG>=10) THEN
        INFO%LCHCON=.TRUE.
        IF(IO%IU6>=0)  WRITE(TIU6,*)'charge density and potential remain constant during run'
        MIX%IMIX=0
      ELSE
        INFO%LCHCON=.FALSE.
        IF(IO%IU6>=0)  WRITE(TIU6,*)'charge density and potential will be updated during run'
      ENDIF

      IF ((WDES%ISPIN==2).AND.INFO%LCHCON.AND.(DYN%IBRION/=-1)) THEN
         IF (IO%IU0>=0)  &
         WRITE(TIU0,*) &
          'Spin polarized Harris functional dynamics is a good joke ...'
         CALL vtutor%error("Spin polarized Harris functional dynamics is a good joke ...")
      ENDIF
      IF (IO%IU6>=0) THEN
        IF (WDES%ISPIN==1 .AND. .NOT. WDES%LNONCOLLINEAR ) THEN
          WRITE(TIU6,*)'non-spin polarized calculation'
        ELSE IF ( WDES%LNONCOLLINEAR ) THEN
          WRITE(TIU6,*)'non collinear spin polarized calculation'
        ELSE
          WRITE(TIU6,*)'spin polarized calculation'
        ENDIF
      ENDIF

! paritial dos
      JOBPAR=1
!     IF (DYN%IBRION>=0) JOBPAR=0
      DO NT=1,T_INFO%NTYP
         IF (T_INFO%RWIGS(NT)<=0._q) JOBPAR=0
      ENDDO

!  INFO%LCDIAG  call EDDIAG after  eigenvalue optimization
!  INFO%LPDIAG  call EDDIAG before eigenvalue optimization
!  INFO%LDIAG   perform sub space rotation (when calling EDDIAG)
!  INFO%LORTHO  orthogonalization of wavefunctions within optimization
!                     no Gramm-Schmidt required
!  INFO%LRMM    use RMM-DIIS minimization
!  INFO%LDAVID  use blocked Davidson
!  INFO%LCHCON  charge constant during run
!  INFO%LCHCOS  charge constant during band minimisation
!  INFO%LONESW  all band simultaneous

      INFO%LCHCOS=.TRUE.
      INFO%LONESW=.FALSE.
      INFO%LONESW_AUTO=.FALSE.
      INFO%LDAVID=.FALSE.
      INFO%LRMM  =.FALSE.
      INFO%LORTHO=.TRUE.
      INFO%LCDIAG=.FALSE.
      INFO%LPDIAG=.TRUE.
      INFO%LPRECONDH=.FALSE.
      INFO%IHARMONIC=0
      INFO%LEXACT_DIAG=.FALSE.
      INFO%LDAVID_FULL=.FALSE.

!  Davidson (non-blocked)
      IF (INFO%IALGO>=110) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Davidson (non-blocked)'
        INFO%IALGO=INFO%IALGO-110
        INFO%LCDIAG=.FALSE.  ! routine does the diagonalisation itself
        INFO%LPDIAG=.FALSE.  ! hence LPDIAG and LCDIAG are set to .FALSE.
        INFO%LDAVID_FULL=.TRUE.
!  all bands conjugate gradient (CG) or damped orbital optimization
!  with fall back to DIIS when convergence is save
      ELSE IF (INFO%IALGO>=100) THEN
        IF (INFO%LDIAG) THEN
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Freysoldt, et al. PRB 79, 241103 (2009))'
        ELSE
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Kresse, et al. variant)'
        ENDIF
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Fall back to RMM-DIIS when convergence is save'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LCHCOS=INFO%LCHCON
        INFO%LONESW=.TRUE.
        INFO%LONESW_AUTO=.TRUE.
        INFO%LRMM  =.TRUE.    ! this is tricky, set some usefull defaults for calling
        INFO%LORTHO=.FALSE.   !  routines in electron.F  (fall back is DIIS)
!  exact diagonalization
      ELSE IF (INFO%IALGO>=90) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Exact diagonalization'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LEXACT_DIAG=.TRUE.
        INFO%LORTHO=.TRUE.
        INFO%LCDIAG=.FALSE.
        INFO%LPDIAG=.FALSE.
!  routines implemented in david_inner (Harmonic Jacobi Davidson, and Davidson)
      ELSE IF (INFO%IALGO>=80 .OR. (INFO%IALGO>=70 .AND. INFO%EREF==0)) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Davidson algorithm suitable for deep iteration (NRMM>8)'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%IHARMONIC=2
        INFO%LORTHO=.FALSE.
        INFO%LCDIAG=.TRUE.
        INFO%LPDIAG=.FALSE.
      ELSE IF (INFO%IALGO>=70) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Jacobi Davidson Harmonic (JDH) for inner eigenvalue problems'
        INFO%IALGO=INFO%IALGO-70
        INFO%IHARMONIC=1
        INFO%LORTHO=.FALSE.
        INFO%LCDIAG=.TRUE.
        INFO%LPDIAG=.FALSE.
!  RMM-DIIS + Davidson
      ELSE IF (INFO%IALGO>=60) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'RMM-DIIS sequential band-by-band and'
        IF(IO%IU6>=0)  WRITE(TIU6,*) ' variant of blocked Davidson during initial phase'
        INFO%IALGO=INFO%IALGO-60
        INFO%LRMM   =.TRUE.
        INFO%LDAVID =.TRUE.
        INFO%LORTHO =.FALSE.
        INFO%LDIAG  =.TRUE.        ! subspace rotation is always selected
!  all bands conjugate gradient (CG) or damped orbital optimization
      ELSE IF (INFO%IALGO>=50) THEN
        IF (INFO%LDIAG) THEN
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Freysoldt, et al. PRB 79, 241103 (2009))'
        ELSE
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Kresse, et al. variant)'
        ENDIF
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LCHCOS=INFO%LCHCON
        INFO%LONESW=.TRUE.
!  RMM-DIIS
      ELSE IF (INFO%IALGO>=40) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'RMM-DIIS sequential band-by-band'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LRMM  =.TRUE.
        INFO%LORTHO=.FALSE.
!  blocked Davidson (Liu)
      ELSE IF (INFO%IALGO>=30) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Variant of blocked Davidson'
        INFO%IALGO=INFO%IALGO-30
        IF (INFO%LDIAG) THEN    ! if LDIAG is set
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Davidson routine will perform the subspace rotation'
           INFO%LCDIAG=.FALSE.  ! routine does the diagonalisation itself
           INFO%LPDIAG=.FALSE.  ! hence LPDIAG and LCDIAG are set to .FALSE.
        ENDIF
        INFO%LDAVID=.TRUE.
      ELSE IF (INFO%IALGO>=20) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient sequential band-by-band (Teter, Alan, Payne)'
        INFO%IALGO  =INFO%IALGO-20
        INFO%LORTHO=.FALSE.
        CALL vtutor%write(isError, IALGO8)
      ELSE IF (INFO%IALGO>=10) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Compatibility mode'
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient sequential band-by-band (Teter, Alan, Payne)'
        INFO%IALGO=INFO%IALGO-10
        INFO%LCDIAG=.TRUE.
        INFO%LPDIAG=.FALSE.
        CALL vtutor%write(isError, IALGO8)
      ELSE IF (INFO%IALGO>=5 .OR. INFO%IALGO==0) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient sequential band-by-band (Teter, Alan, Payne)'
        CALL vtutor%write(isError, IALGO8)
      ELSE IF (INFO%IALGO<0) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Performance tests'
      ELSE IF (INFO%IALGO <1) THEN
        CALL vtutor%error("Algorithms no longer implemented")
      ELSE IF (INFO%IALGO==2) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'None: do nothing, only one-electron occupancies are recalculated'
        INFO%LDIAG =.FALSE.
      ELSE IF (INFO%IALGO==3) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Eigenval: update one-electron energies, occupancies fixed'
        INFO%LDIAG =.FALSE.
      ELSE IF (INFO%IALGO==4) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Subrot: only subspace diagonalization (rotation)'
      ENDIF

      SZ=' '
      IF (INFO%LCHCOS) THEN
         SZ='   charged. constant during bandupdate'
      ELSE
         INFO%LCORR=.FALSE.
!        IMIX used to be set to 0, removed
!         MIX%IMIX=0
      ENDIF

      IF (IO%IU6>=0) THEN

      IF (.NOT. INFO%LRMM .AND. .NOT. INFO%LDAVID) THEN
        IF (INFO%IALGO==5) THEN
          WRITE(TIU6,*)'steepest descent',SZ
        ELSEIF (INFO%IALGO==6) THEN
          WRITE(TIU6,*)'conjugated gradient',SZ
        ELSEIF (INFO%IALGO==7) THEN
          WRITE(TIU6,*)'preconditioned steepest descent',SZ
        ELSEIF (INFO%IALGO==8) THEN
          WRITE(TIU6,*)'preconditioned conjugated gradient',SZ
        ELSEIF (INFO%IALGO==0) THEN
          WRITE(TIU6,*)'preconditioned conjugated gradient (Jacobi prec)',SZ
        ENDIF
        IF (.NOT.INFO%LONESW) THEN
          WRITE(TIU6,*)'   band-by band algorithm'
        ENDIF
      ENDIF

      IF (INFO%LDIAG) THEN
        WRITE(TIU6,*)'perform sub-space diagonalisation'
      ELSE
        WRITE(TIU6,*)'perform Loewdin sub-space diagonalisation'
        WRITE(TIU6,*)'   ordering is kept fixed'
      ENDIF

      IF (INFO%LPDIAG) THEN
        WRITE(TIU6,*)'   before iterative eigenvector-optimisation'
      ELSE
        WRITE(TIU6,*)'   after iterative eigenvector-optimisation'
      ENDIF

      IF (MIX%IMIX==1 .OR. MIX%IMIX==2 .OR. MIX%IMIX==3) THEN
        WRITE(TIU6,*)'Kerker-like  mixing scheme'
      ELSE IF (MIX%IMIX==4) THEN
       WRITE(TIU6,'(A,F10.1)')' modified Broyden-mixing scheme, WC = ',MIX%WC
       IF (MIX%INIMIX==1) THEN
         WRITE(TIU6,'(A,F8.4,A,F12.4)') &
     &     ' initial mixing is a Kerker type mixing with AMIX =',MIX%AMIX, &
     &     ' and BMIX =',MIX%BMIX
       ELSE IF (MIX%INIMIX==2) THEN
         WRITE(TIU6,*)'initial mixing equals unity matrix (no mixing!)'
       ELSE
         WRITE(TIU6,'(A,F8.4)') &
     &     ' initial mixing is a simple linear mixing with ALPHA =',MIX%AMIX
       ENDIF
       IF (MIX%MIXPRE==1) THEN
         WRITE(TIU6,*)'Hartree-type preconditioning will be used'
       ELSE IF (MIX%MIXPRE==2) THEN
         WRITE(TIU6,'(A,A,F12.4)') &
     &     ' (inverse) Kerker-type preconditioning will be used', &
     &     ' corresponding to BMIX =',MIX%BMIX
       ELSE
         WRITE(TIU6,*)'no preconditioning will be used'
       ENDIF
      ELSE
        WRITE(TIU6,*)'no mixing'
      ENDIF
      IF (WDES%NB_TOT*2==NINT(INFO%NELECT)) THEN
        WRITE(TIU6,*)'2*number of bands equal to number of electrons'
        IF (MIX%IMIX/=0 .AND..NOT.INFO%LCHCOS) THEN
          WRITE(TIU6,*) &
     &      'WARNING: mixing without additional bands will not converge'
        ELSE IF (MIX%IMIX/=0) THEN
          WRITE(TIU6,*) 'WARNING: mixing has no effect'
        ENDIF

      ELSE
        WRITE(TIU6,*)'using additional bands ',INT(WDES%NB_TOT-INFO%NELECT/2)
        IF (KPOINTS%SIGMA<=0) THEN
          WRITE(TIU6,*) &
     &  'WARNING: no broadening specified (might cause bad convergence)'
        ENDIF
      ENDIF

      IF (INFO%LREAL) THEN
        WRITE(TIU6,*)'real space projection scheme for non local part'
      ELSE
        WRITE(TIU6,*)'reciprocal scheme for non local part'
      ENDIF

      IF (INFO%LCORE) THEN
        WRITE(TIU6,*)'use partial core corrections'
      ENDIF

      IF (INFO%LCORR) THEN
        WRITE(TIU6,*)'calculate Harris-corrections to forces ', &
     &              '  (improved forces if not selfconsistent)'
      ELSE
        WRITE(TIU6,*)'no Harris-corrections to forces '
      ENDIF

      IF (XC_DATA%LDOGGA.OR.XC_DATA%LDOMETAGGA) THEN
        WRITE(TIU6,*)'use gradient corrections '
        IF (INFO%LCHCON) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*)'WARNING: stress and forces are not correct'
           WRITE(TIU6,*)'WARNING: stress and forces are not correct'
           WRITE(TIU6,*)' (second derivative of E(xc) not defined)'
        ENDIF
      ENDIF

      IF (INFO%LOVERL) THEN
         WRITE(TIU6,*)'use of overlap-Matrix (Vanderbilt PP)'
      ENDIF
      IF (KPOINTS%ISMEAR==-1) THEN
        WRITE(TIU6,7213) KPOINTS%SIGMA
 7213 FORMAT(' Fermi-smearing in eV        SIGMA  = ',F6.2)

      ELSE IF (KPOINTS%ISMEAR==-2) THEN
        WRITE(TIU6,7214)
 7214 FORMAT(' partial occupancies read from INCAR or WAVECAR (fixed during run)')
      ELSE IF (KPOINTS%ISMEAR==-4) THEN
        WRITE(TIU6,7222)
 7222 FORMAT(' Fermi weights with tetrahedron method witout', &
     &       ' Bloechl corrections')
      ELSE IF (KPOINTS%ISMEAR==-5) THEN
        WRITE(TIU6,7223)
 7223 FORMAT(' Fermi weights with tetrahedron method with', &
     &       ' Bloechl corrections')

      ELSE IF (KPOINTS%ISMEAR>0) THEN
        WRITE(TIU6,7217) KPOINTS%ISMEAR,KPOINTS%SIGMA
 7217 FORMAT(' Methfessel and Paxton  Order N=',I2, &
     &       ' SIGMA  = ',F6.2)
      ELSE
        WRITE(TIU6,7218) KPOINTS%SIGMA
 7218 FORMAT(' Gauss-broadening in eV      SIGMA  = ',F6.2)
      ENDIF

      WRITE(TIU6,130)
!=======================================================================
!  write out the lattice parameters
!=======================================================================
      WRITE(TIU6,7220) INFO%ENMAX,LATT_CUR%OMEGA, &
     &    ((LATT_CUR%A(I,J),I=1,3),(LATT_CUR%B(I,J),I=1,3),J=1,3), &
     &    (LATT_CUR%ANORM(I),I=1,3),(LATT_CUR%BNORM(I),I=1,3)

      WRITE(TIU6,*)

      IF (INFO%ISTART==1 .OR.INFO%ISTART==2) THEN

      WRITE(TIU6,*)'old parameters found on file WAVECAR:'
      WRITE(TIU6,7220) ENMAXI,LATT_INI%OMEGA, &
     &    ((LATT_INI%A(I,J),I=1,3),(LATT_INI%B(I,J),I=1,3),J=1,3)


      WRITE(TIU6,*)
 7220 FORMAT('  energy-cutoff  :  ',F10.2/ &
     &       '  volume of cell :  ',F10.2/ &
     &       '      direct lattice vectors',17X,'reciprocal lattice vectors'/ &
     &       3(2(3X,3F13.9)/) / &
     &       '  length of vectors'/ &
     &        (2(3X,3F13.9)/) /)

      ENDIF
!=======================================================================
!  write out k-points,weights,size, positions & external forces
!=======================================================================

 7104 FORMAT(' k-points in units of 2pi/SCALE and weight: ',A40)
 7105 FORMAT(' k-points in reciprocal lattice and weights: ',A40)
 7016 FORMAT(' position of ions in fractional coordinates (direct lattice) ')
 7017 FORMAT(' position of ions in cartesian coordinates  (Angst):')
 7018 FORMAT(' external forces on ions in cartesian coordinates  (eV/Angst):')
 7009 FORMAT(1X,3F12.8,F12.3)
 7007 FORMAT(1X,3F12.8)

      WRITE(TIU6,7104) KPOINTS%SZNAMK

      DO NKP=1,KPOINTS%NKPTS
        VTMP(1)=WDES%VKPT(1,NKP)
        VTMP(2)=WDES%VKPT(2,NKP)
        VTMP(3)=WDES%VKPT(3,NKP)
        CALL DIRKAR(1,VTMP,LATT_CUR%B)
        WRITE(TIU6,7009) VTMP(1)*LATT_CUR%SCALE,VTMP(2)*LATT_CUR%SCALE, &
                  VTMP(3)*LATT_CUR%SCALE,KPOINTS%WTKPT(NKP)
      ENDDO

      WRITE(TIU6,*)
      WRITE(TIU6,7105) KPOINTS%SZNAMK
      DO NKP=1,KPOINTS%NKPTS
        WRITE(TIU6,7009) WDES%VKPT(1,NKP),WDES%VKPT(2,NKP),WDES%VKPT(3,NKP),KPOINTS%WTKPT(NKP)
      ENDDO
      WRITE(TIU6,*)

      WRITE(TIU6,7016)
      WRITE(TIU6,7007) ((DYN%POSION(I,J),I=1,3),J=1,T_INFO%NIONS)
      WRITE(TIU6,*)
      WRITE(TIU6,7017)

      DO J=1,T_INFO%NIONS
        VTMP(1)=DYN%POSION(1,J)
        VTMP(2)=DYN%POSION(2,J)
        VTMP(3)=DYN%POSION(3,J)
        CALL  DIRKAR(1,VTMP,LATT_CUR%A)
        WRITE(TIU6,7007) (VTMP(I),I=1,3)
      ENDDO

      ! write out external forces if any are present
      IF (ANY(ABS(DYN%EFOR) > 1.E-5_q)) THEN
         WRITE(TIU6,*)
         WRITE(TIU6,7018)
         DO J=1,T_INFO%NIONS
            VTMP(1)=DYN%EFOR(1,J)
            VTMP(2)=DYN%EFOR(2,J)
            VTMP(3)=DYN%EFOR(3,J)
            WRITE(TIU6,7007) (VTMP(I),I=1,3)
         ENDDO
      ENDIF

      WRITE(TIU6,*)
      CALL  WRITE_EULER(IO%IU6, WDES%LNONCOLLINEAR, WDES%SAXIS)

      WRITE(TIU6,130)
      ENDIF

      io_begin
      IF (KIMAGES==0) THEN
!=======================================================================
!  write out initial header for PCDAT, XDATCAR
!=======================================================================
      IF (OUTPUT_MODE.GT.0) THEN
         CALL PCDAT_HEAD(60,T_INFO, LATT_CUR, DYN, PACO, INFO%SZNAM1)
      ENDIF
      ! for DYN%ISIF>=3 header is written at every ionic loop step
      IF ( DYN%ISIF<3 ) THEN
        CALL XDAT_HEAD(61, T_INFO, LATT_CUR%A, INFO%SZNAM1)
      ENDIF
#ifdef VASP_HDF5
      ! write the initial data to the history for methods that do not relax ions
      CALL VH5_WRITE_LATTICE_ION_HIST_START(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,T_INFO%NIONS, DYN)
      CALL VH5_WRITE_LATTICE_ION_HIST_STEP(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY, 1, LATT_CUR%A, DYN%POSIOC)
      CALL VH5_WRITE_LATTICE_ION(IH5OUTFILEID, INFO%SZNAM1, T_INFO, LATT_CUR%SCALE, LATT_CUR%A, T_INFO%LSDYN, DYN%POSIOC)
      CALL VH5_WRITE_ENERGIES_HIST_START(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,DYN%IBRION)
      CALL VH5_WRITE_STRESS_FORCES_HIST_START(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,T_INFO%NIONS, DYN)
      CALL VH5_WRITE_OSZICAR_START(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY)
      IF (IO%LORBIT >= 10) &
         CALL VH5_WRITE_MOMENTS_HIST_START(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,[LDIMP,T_INFO%NIONS,NCDIJ],"spin", DYN)
      IF (LCALC_ORBITAL_MOMENT().AND.WDES%LNONCOLLINEAR) &
         CALL VH5_WRITE_MOMENTS_HIST_START(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,[LDIMP-1,T_INFO%NIONS,3],"orbital", DYN)
      IF (IO%VELOCITY) &
         CALL VH5_WRITE_VELOCITY_HIST_START(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY, T_INFO%NIONS, DYN)
      IF (DYN%IBRION == 0) &
         CALL VH5_WRITE_PAIR_CORRELATION_START(IH5INTERMEDIATEGROUP_ID, GRP_PAIR_CORRELATION, T_INFO, DYN, PACO)
#endif
!=======================================================================
!  write out initial header for DOS
!=======================================================================
      JOBPAR_=JOBPAR
      IF (IO%LORBIT>=10 ) JOBPAR_=1

      WRITE(16,'(4I4)') T_INFO%NIONP,T_INFO%NIONS,JOBPAR_,WDES%NCDIJ
      WRITE(16,'(5E15.7)')AOMEGA,((LATT_CUR%ANORM(I)*1E-10),I=1,3),DYN%POTIM*1E-15
      WRITE(16,*) DYN%TEMP
      WRITE(16,*) ' CAR '
      WRITE(16,*) INFO%SZNAM1

!=======================================================================
!  write out initial header for EIGENVALUES
!=======================================================================
      WRITE(22,'(4I5)') T_INFO%NIONS,T_INFO%NIONS,DYN%NBLOCK*DYN%KBLOCK,WDES%ISPIN
      WRITE(22,'(5E15.7)') &
     &         AOMEGA,((LATT_CUR%ANORM(I)*1E-10_q),I=1,3),DYN%POTIM*1E-15_q
      WRITE(22,*) DYN%TEMP
      WRITE(22,*) ' CAR '
      WRITE(22,*) INFO%SZNAM1
      WRITE(22,'(3I7)') NINT(INFO%NELECT),KPOINTS%NKPTS,WDES%NB_TOT
      ENDIF
!=======================================================================
!  write hdf5 header info
!=======================================================================

#ifdef VASP_HDF5
      CALL VH5_WRITE_DOSHEADER(IH5OUTFILEID, T_INFO, JOBPAR_, DYN, WDES, LATT_CUR, INFO, KPOINTS)
      CALL VH5_WRITE_EIGENVALHEADER(IH5OUTFILEID, T_INFO, DYN, WDES, LATT_CUR, INFO, KPOINTS)
#endif
      io_end

#ifdef VASP_HDF5
      IF (INCAR_FOUND .AND. POSCAR_FOUND .AND. KPOINTS_FOUND) THEN
#endif
      IF (IO%IU0>=0) &
      WRITE(TIU0,*)'POSCAR, INCAR and KPOINTS ok, starting setup'
#ifdef VASP_HDF5
      ELSE
        IF (INCAR_FOUND) THEN
          IF (IO%IU0>=0) &
            WRITE(TIU0,*)'INCAR ok, starting setup'
        ENDIF
        IF (POSCAR_FOUND) THEN
          IF (IO%IU0>=0) &
            WRITE(TIU0,*)'POSCAR ok, starting setup'
        ENDIF
        IF (KPOINTS_FOUND) THEN
          IF (IO%IU0>=0) &
            WRITE(TIU0,*)'KPOINTS ok, starting setup'
        ENDIF
        IF (HDF5_FOUND) THEN
          IF (IO%IU0>=0) &
            WRITE(TIU0,*) VASPIN // ' ok, starting setup'
        ENDIF
      ENDIF
#endif

!=======================================================================
! initialize the required grid structures
!=======================================================================
      CALL INILGRD(NGX,NGY,NGZ,GRID)
      CALL INILGRD(NGX,NGY,NGZ,GRID_SOFT)
      CALL INILGRD(NGXC,NGYC,NGZC,GRIDC)
      CALL INILGRD(GRIDUS%NGPTAR(1),GRIDUS%NGPTAR(2),GRIDUS%NGPTAR(3),GRIDUS)
#ifdef MPI
! only wavefunction grid uses local communication
      GRID%COMM     =>COMM_INB
! all other grids use world wide communication at the moment set their
! communication boards to the world wide communicator

!PK Charge density grids are replicated per set of distributed k-points
      GRID_SOFT%COMM=>COMM_KIN
      GRIDC%COMM    =>COMM_KIN
      GRIDUS%COMM   =>COMM_KIN
      GRIDB%COMM    =>COMM_KIN

      GRIDC%COMM_KIN    =>COMM_KIN
      GRIDC%COMM_KINTER =>COMM_KINTER
#endif
      CALL GEN_RC_GRID(GRIDUS)
#ifdef usgrid
      CALL GEN_RC_SUB_GRID(GRIDC,GRIDUS, C_TO_US, .TRUE.,.TRUE.)
#else
      CALL GEN_RC_GRID(GRIDC)
#endif
      CALL GEN_RC_SUB_GRID(GRID_SOFT, GRIDC, SOFT_TO_C, .TRUE.,.TRUE.)
!=======================================================================
!  allocate work arrays
!=======================================================================
!
! GEN_LAYOUT determines the data layout (distribution) of the columns on parallel
! computers and allocates all required arrays of the WDES descriptor
      IF (INFO%ISTART==1) THEN
         CALL GEN_LAYOUT(GRID,WDES, LATT_CUR%B,LATT_CUR%B,IO%IU6,.TRUE.)
         CALL GEN_INDEX(GRID,WDES, LATT_CUR%B,LATT_CUR%B,IO%IU6,IO%IU0,.TRUE.)
! all other cases use LATT_INI for setup of GENSP
      ELSE
         DWRITE0 'call to genlay'
         CALL GEN_LAYOUT(GRID,WDES, LATT_CUR%B,LATT_INI%B,IO%IU6,.TRUE.)
         DWRITE0 'call to genind'
         CALL GEN_INDEX(GRID,WDES, LATT_CUR%B,LATT_INI%B,IO%IU6,IO%IU0,.TRUE.)
      ENDIF
      CALL SET_NBLK_NSTRIP( WDES)
!
! non local projection operators
!

      IF(LDO_AB_INITIO) THEN
         CALL NONL_ALLOC(NONL_S,T_INFO,P,WDES,INFO%LREAL)
#ifndef shmem_rproj
         CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WDES%LSPIRAL)
#else
         CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WDES%LSPIRAL,WDES%COMM_intra_node)
#endif
!  optimize grid for real space representation and calculate IRMAX, IRALLOC
         NONLR_S%IRMAX=0 ; NONLR_S%IRALLOC=0
         CALL REAL_OPTLAY(GRID,LATT_CUR,NONLR_S,LPLANE_WISE,LREALLOCATE, IO%IU6, IO%IU0)
! allign GRID_SOFT with GRID in real space
         CALL SET_RL_GRID(GRID_SOFT,GRID)
! allocate real space projectors
         CALL NONLR_ALLOC(NONLR_S)
!  init FFT
         CALL FFTINI(WDES%NINDPW(1,1),WDES%NGVECTOR(1),KPOINTS%NKPTS,WDES%NGDIM,GRID)
#ifdef MPI
         CALL MAPSET(GRID)   ! generate the communication maps (patterns) for FFT
         IF (IO%IU6 >=0) THEN
            IF (GRID%RL%NFAST==1) THEN
               WRITE(TIU6,"(/' serial   3D FFT for wavefunctions')")
            ELSE IF (GRID%IN_RL%LOCAL  ) THEN
               WRITE(TIU6,"(/' parallel 3D FFT for wavefunctions:'/'    minimum data exchange during FFTs selected (reduces bandwidth)')")
            ENDIF
         ENDIF
         DWRITE0  'mapset aug done'
         CALL MAPSET(GRID_SOFT)
         DWRITE0  'mapset soft done'
         CALL MAPSET(GRIDC)
         IF (GRIDC%IN_RL%LOCAL .AND. IO%IU6 >=0) THEN
            WRITE(TIU6,"(' parallel 3D FFT for charge:'/'    minimum data exchange during FFTs selected (reduces bandwidth)'/)")
         ENDIF
         DWRITE0  'mapset wave done'
#ifdef usgrid
         CALL MAPSET(GRIDUS)
#endif
#endif
      ENDIF

      IF (IO%DRY_RUN) THEN
         IF (IO%IU6 >= 0) THEN
            WRITE(TIU0,*) 'Dry run completed, exiting VASP'
            WRITE(TIU6,*) 'Dry run completed, exiting VASP'
         ENDIF
         CALL VTUTOR%STOPCODE()
      ENDIF

      IF(LDO_AB_INITIO) THEN
#ifdef OFFLOADING
         CALL FFT_OFFLOAD_MAKEPLANS(GRID,WDES%NSIM*2)
         CALL FFT_OFFLOAD_MAKEPLANS(GRIDC)
         CALL FFT_OFFLOAD_MAKEPLANS(GRID_SOFT)
#endif
      ENDIF

! allocate all other arrays
!
      IF(LDO_AB_INITIO) THEN
         ISP     =WDES%ISPIN
         NCDIJ   =WDES%NCDIJ
         LMDIM   =P(1)%LMDIM
      ENDIF
      NTYP    =T_INFO%NTYP
      NIOND   =T_INFO%NIOND
      NIOND_LOC=WDES%NIONS

      IF(LDO_AB_INITIO) THEN
         ALLOCATE(CHTOT(GRIDC%MPLWV,NCDIJ),SV(DIMREAL(GRID%MPLWV),NCDIJ))
         ALLOCATE(CHTOTL(GRIDC%MPLWV,NCDIJ),DENCOR(GRIDC%RL%NP), &
                  CVTOT(GRIDC%MPLWV,NCDIJ),CSTRF(GRIDC%MPLWV,NTYP), &
! small grid quantities
                  CHDEN(GRID_SOFT%MPLWV,NCDIJ), &
! non local things
                  CDIJ(LMDIM,LMDIM,NIOND_LOC,NCDIJ), &
                  CQIJ(LMDIM,LMDIM,NIOND_LOC,NCDIJ), &
                  CRHODE(LMDIM,LMDIM,NIOND_LOC,NCDIJ), &
! forces (depend on NIOND)
                  EWIFOR(3,NIOND), &
! dos
                  DOS(NEDOS,NCDIJ),DOSI(NEDOS,NCDIJ), &
                  DDOS(NEDOS,NCDIJ),DDOSI(NEDOS,NCDIJ), &
                  PAR(1,1,1,1,NCDIJ),DOSPAR(1,1,1,NCDIJ))
DOESI    DOS = 0.0_q; DOSI = 0.0_q ! Breaks bulk_CoO_wurzite_SOC_G0W0_sym (intel).
DOESI    SV = 0 ! Breaks SiC_OEP (intel), array is used as a whole in OEP_GW.

         IF (IO%IU0>=0) WRITE(TIU0,*)'FFT: planning ... GRIDC'
         CALL FFTGRIDPLAN(GRIDC)

         IF (IO%IU0>=0) WRITE(TIU0,*)'FFT: planning ... GRID_SOFT'
         CALL FFTGRIDPLAN(GRID_SOFT)
         ! If the coulomb kernel truncation is used, setup the associated grids here
         IF (HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%ACTIVE) THEN
            !
            ! Create the fine and coarse grids
            ! The coarse and fine grids are passed to main.F to ensure that they stay in scope
            !
            CALL CREATE_GRIDS_FOR_KERNEL_TRUNCATION(GRIDC, HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE, COMM_KIN, &
                                                    GRID_COARSE, GRID_FINE, TABLE_CUTOFF)
            IF (IO%IU0>=0) WRITE(TIU0,*)'Generated coarse and fine grid for coulomb truncation...'
            IF (HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%COARSEN_BEFORE_PAD) THEN
__GPU           CALL FFT_OFFLOAD_MAKEPLANS(HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%COARSE_GRID)
            ELSE
__GPU           CALL FFT_OFFLOAD_MAKEPLANS(HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%FINE_GRID)
            END IF
            !
         ENDIF
         !
#ifdef OFFLOADING
         CALL OFFLOAD_ON_PUSH(.TRUE.)
GPACC !$ACC ENTER DATA COPYIN(HAMILTONIAN) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:HAMILTONIAN) IF(OFFLOAD_ON)
         CALL COPYIN_TYPED_VAR(HAMILTONIAN%COULOMB_POT)
GPOMP IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:HAMILTONIAN%COULOMB_POT)
         CALL OFFLOAD_ON_POP
#endif
         !
      ELSE
          ALLOCATE(EWIFOR(3,NIOND))
          EWIFOR=0.0_q
      ENDIF

! TIFOR and SIPACO must be allocated here because these two arrays are used even
! when LDO_AB_INITIO=.FALSE.
      ALLOCATE (TIFOR(3,NIOND))
      TIFOR = 0.0_q
#if defined ML_AVAILABLE || !defined ML_LEARNING_NOT_ACTIVE
      ALLOCATE (TIFOR_ML(3,NIOND))
      TIFOR_ML = 0.0_q
#endif
      IF (OUTPUT_MODE.GT.0) THEN
         CALL PCDAT_ALLOCATE_PACO( T_INFO, PACO)
      ENDIF


      IF(LDO_AB_INITIO) THEN
         CALL REGISTER_ALLOCATE(16._q*(SIZE(CHTOT)+SIZE(CHTOTL)+SIZE(DENCOR)+SIZE(CVTOT)+SIZE(CSTRF)+SIZE(CHDEN)+SIZE(SV)), "grid")
         IF (WDES%LGAMMA) THEN
            CALL REGISTER_ALLOCATE(8._q*(SIZE(CDIJ)+SIZE(CQIJ)+SIZE(CRHODE)), "one-center")
         ELSE
            CALL REGISTER_ALLOCATE(16._q*(SIZE(CDIJ)+SIZE(CQIJ)+SIZE(CRHODE)), "one-center")
         ENDIF

         CALL ALLOCATE_AVEC(HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT, GRID, GRIDC)

         CALL ALLOCATE_MU(XC_DATA%LMU, HAMILTONIAN%MU, HAMILTONIAN%MUTOT, GRID, GRIDC, WDES)
! test_
!     CALL GENERATE_TAU_HANDLE(XC_DATA%LTAU, KINEDEN, GRIDC, WDES%NCDIJ)
         CALL ALLOCATE_TAU(XC_DATA%LTAU, KINEDEN%TAU,KINEDEN%TAUC,KINEDEN%TAUL,GRIDC,WDES%NCDIJ)
! test_
         CALL CREATE_CMBJ_AUX(GRIDC,T_INFO,LATT_CUR,XC_DATA)
         CALL CREATE_CMBJ_AUG(T_INFO,XC_DATA)

         DWRITE0 'allocation done'
!
         IF (IO%IU0>=0) WRITE(TIU0,*)'FFT: planning ... GRID'
         CALL FFTGRIDPLAN(GRID)

         MPLMAX=MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV,GRID%MPLWV)
#ifdef MPI
! give T3D opportunity to allocate all required shmem workspace
         CALL SHM_MAX(WDES, MPLMAX, MALLOC)
         CALL SHM_ALLOC(MALLOC)
#endif
         MIX%NEIG=0
! calculate required numbers elements which must be mixed in PAW
! set table for Clebsch-Gordan coefficients, maximum L is 2*3 (f states)
         LMAX_TABLE=6;  CALL YLM3ST_(LMAX_TABLE)

         N_MIX_PAW=0
         CALL SET_RHO_PAW_ELEMENTS(WDES, P , T_INFO, INFO%LOVERL, N_MIX_PAW )
         ALLOCATE( RHOLM(N_MIX_PAW,WDES%NCDIJ), RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ))
! solve pseudo atomic problem
         IF (LPRJ_LCAO()) CALL LCAO_INIT(P,XC_DATA,INFO,IO%IU0,IO%IU6)
! setup fock (requires YLM)
         CALL SETUP_FOCK(T_INFO, P, XC_DATA, WDES, GRID, LATT_CUR, LMDIM,  INFO%SZPREC, IO%IU6, IO%IU0 )
! setup atomic PAW terms
         IF (LRHFATM()) THEN
            CALL SET_PAW_ATOM_POT_RHF(P,XC_DATA,T_INFO,INFO%LOVERL,INFO%EALLAT,IO)
         ELSE
           CALL RHFATM_CROP_PSEUDO(P,T_INFO,INFO%LOVERL,XC_DATA,IO)
           CALL SET_PAW_ATOM_POT(P, XC_DATA, T_INFO, INFO%LOVERL,  &
                LMDIM, INFO%EALLAT, IO%IU6  )
         ENDIF
! possibly setup of the scf determination of the subspace rotation
         CALL SETUP_SUBROT_SCF(INFO,WDES,LATT_CUR,GRID,GRIDC,GRID_SOFT,SOFT_TO_C,IO%IU0,IO%IU5,IO%IU6)
! setup pead
         CALL PEAD_SETUP(WDES,GRID,GRIDC,GRIDUS,C_TO_US,KPOINTS,LATT_CUR,LATT_INI,T_INFO,P,LMDIM,INFO%LOVERL,IRDMAX,IO)
      ENDIF
!=======================================================================
! allocate wavefunctions
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         CALL ALLOCW(WDES,W,WUP,WDW)
         IF (INFO%LONESW) THEN
           CALL ALLOCW(WDES,W_F,WTMP,WTMP)
           CALL ALLOCW(WDES,W_G,WTMP,WTMP)
           ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
                    CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
         ELSE
           ALLOCATE(CHAM(1,1,1,1),CHF(1,1,1,1))
         ENDIF
         CALL DUMP_ALLOCATE(IO%IU6)
      ENDIF
!=======================================================================
! now read in wavefunctions
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         W%CELTOT=0
         CALL START_TIMING("G")
         IF (INFO%LHALF_INIT) THEN
           IF (INFO%ICHARG/=1) CALL vtutor%error("LHALF_INIT requires ICHARG=1")
           IF (IO%IU0>=0) WRITE(TIU0,*)'initializing wavefunctions with DeePAW-HALF'
           CALL HALF_INITIALIZE_WAVEFUNCTIONS(WDES,W,INFO%ENMAX, &
                DIR_APP(1:DIR_LEN)//'CHGCAR',DIR_APP(1:DIR_LEN)//'POTCAR', &
                LATT_CUR,T_INFO,GRIDC%NGPTAR,IO%IU0,INFO%LHALF_API,INFO%HALF_ESCN_URL)
         ELSE IF (INFO%ISTART>0) THEN
#ifdef VASP_HDF5
           IF (WAVECAR_FOUND) THEN
#endif
           IF (IO%IU0>=0) WRITE(TIU0,*)'reading WAVECAR'
#ifdef VASP_HDF5
           ELSE IF (VASPWAVE_FOUND) THEN
             IF (IO%IU0>=0) WRITE(TIU0,*)'reading vaspwave.h5'
           ENDIF
#endif

#ifdef VASP_HDF5
           IF (WAVECAR_FOUND) THEN
#endif
           IF (IO%LDOWNSAMPLE ) THEN
             CALL INWAV_DOWNFOLD(IO, WDES, W, SYMM, KPOINTS, GRID, LATT_CUR, LATT_INI, INFO%ISTART,  EFERMI )
           ELSE
             CALL INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, INFO%ISTART,  EFERMI )
           ENDIF
#ifdef VASP_HDF5
           ELSE IF (VASPWAVE_FOUND) THEN
               VH5_CHECK( VH5_READ_WAVEFUNCTIONS(IH5WAVEFILEID, IO, WDES, W, SYMM, KPOINTS, GRID, LATT_CUR, LATT_INI, INFO%ISTART, EFERMI) )
           ENDIF
#endif
         ELSE
            IF (IO%IU0>=0) WRITE(TIU0,*)'WAVECAR not read'
         ENDIF
         CALL CLOSEWAV

         IF (WDES%LSPIRAL.AND.(INFO%ISTART>0)) CALL CLEANWAV(WDES,W,INFO%ENINI)
         CALL STOP_TIMING("G",IO%IU6,'INWAV')


         IF (INFO%ISTART/=2) LATT_INI=LATT_CUR

      ENDIF
!=======================================================================
! At this very point everything has been read in
! and we are ready to write all important information
! to the xml file
!=======================================================================
      CALL XML_ATOMTYPES(T_INFO%NIONS, T_INFO%NTYP, T_INFO%NITYP, T_INFO%ITYP, P%ELEMENT, P%POMASS, P%ZVALF, P%SZNAMP )

      CALL XML_TAG("structure","initialpos")
      CALL XML_CRYSTAL(LATT_CUR%A, LATT_CUR%B, LATT_CUR%OMEGA)
      CALL XML_POSITIONS(T_INFO%NIONS, DYN%POSION)
      IF (T_INFO%LSDYN) CALL XML_LSDYN(T_INFO%NIONS,T_INFO%LSFOR(1,1))
      IF (DYN%IBRION<=0 .AND. DYN%NSW>0 ) CALL XML_VEL(T_INFO%NIONS, DYN%VEL)
      IF (T_INFO%LSDYN) CALL XML_NOSE(DYN%SMASS)
      CALL XML_CLOSE_TAG("structure")
!=======================================================================
! initialize index tables for broyden mixing
!=======================================================================
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (((MIX%IMIX==4).AND.(.NOT.INFO%LCHCON)).OR.DYN%IBRION==10) THEN
! Use a reduced mesh but only if using preconditioning ... :
            IF (EXXOEP>=1) THEN
               CALL BRGRID(GRIDC,GRIDB,MAX(INFO%ENMAX*2,ENCUTGW_IN_CHI()),IO%IU6,LATT_CUR%B)
            ELSE
               CALL BRGRID(GRIDC,GRIDB,INFO%ENMAX,IO%IU6,LATT_CUR%B)
            ENDIF
            CALL INILGRD(GRIDB%NGX,GRIDB%NGY,GRIDB%NGZ,GRIDB)
            CALL GEN_RC_SUB_GRID(GRIDB,GRIDC, B_TO_C, .FALSE.,.TRUE.)
         ENDIF
! calculate the structure-factor, initialize some arrays
         IF (INFO%TURBO==0) CALL STUFAK(GRIDC,T_INFO,CSTRF)
         CHTOT=0 ; CHDEN=0; CVTOT=0
      ENDIF
!=======================================================================
! construct initial  charge density:  a bit of heuristic is used
!  to get sensible defaults if the user specifies stupid values in the
!  INCAR files
! for the initial charge density there are several possibilties
! if INFO%ICHARG= 1 read in charge-density from a File
! if INFO%ICHARG= 2-3 construct atomic charge-densities of overlapping atoms
! if INFO%ICHARG= 4 read potential from file
! if INFO%ICHARG= 5 read GAMMA file and update orbitals accordingly
! if INFO%ICHARG >=10 keep chargedensity constant
!
!=======================================================================
      ! subtract 10 from ICHARG (10 means fixed charge density)
      IF(LDO_AB_INITIO) THEN
         IF (INFO%ICHARG>10) THEN
           INFO%INICHG= INFO%ICHARG-10
         ELSE
           INFO%INICHG= INFO%ICHARG
         ENDIF

         IF (INFO%INICHG==5) THEN
            ! no convergence correction for forces (reread GAMMA in last electronic step)
            IF (.NOT. MIX%MIXFIRST) INFO%LCORR=.FALSE.
            CALL CREATE_VASP_LOCK(COMM)
            ! and reset INICHG to 1 so that CHGCAR is read as well
            INFO%INICHG=1
         ENDIF

      ! then initialize CRHODE and than RHOLM (PAW related occupancies)
         CALL DEPATO(WDES, LMDIM, CRHODE, INFO%LOVERL, P, T_INFO)
         CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
              CRHODE, RHOLM)

! MM: Stuff for Janos Angyan: write the AE-charge density for the
         IF (LWRT_AECHG()) THEN
      ! overlapping atomic charges to AECCAR1
         ! add overlapping atomic charges on dense regular grid
            CALL RHOATO_WORK(.TRUE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
         ! add (n_ae - n_ps - n_comp) as defined on radial grid
            CALL AUGCHG(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        &                LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
        &                 LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX,.FALSE.,.FALSE.,IO%IU0)

#ifdef MPI
            IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
         io_begin
         ! write AECCAR1
            OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'AECCAR1',STATUS='UNKNOWN')
         ! write header
            CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
           &             .FALSE., DYN%POSION)
         io_end
         ! write AE charge density
            CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
         io_begin
            CLOSE(99)
         io_end
         ! reset charge densities to zero again
            CHTOT=0; CHDEN=0
      ! AE core density to AECCAR0
         ! add partial core density on dense regular grid
            IF (INFO%LCORE) THEN
               CALL RHOATO_WORK(.TRUE.,.TRUE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
            ENDIF
         ! add core densities (nc_ae-nc_ps) as defined on radial grid
            CALL AUGCHG(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        &                LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
        &                 LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX,.FALSE.,.TRUE.,IO%IU0)
         io_begin
         ! write AECCAR0
            OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'AECCAR0',STATUS='UNKNOWN')
         ! write header
            CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
        &                .FALSE., DYN%POSION)
         io_end
         ! write AE core charge density
            CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
         io_begin
            CLOSE(99)
         io_end
         ! and reset charge densities to zero again

#ifdef MPI
            END IF
#endif
            CHTOT=0; CHDEN=0
         ENDIF
      ENDIF

     ! initial set of wavefunctions from diagonalization of Hamiltonian
     ! set INICHG to 2

      IF(LDO_AB_INITIO) THEN
         IF (INFO%INICHG==0 .AND. INFO%INIWAV==2) THEN
           INFO%INICHG=2
           IF (IO%IU6>=0) &
           WRITE(TIU6,*)'WARNING: no initial charge-density supplied,', &
                          ' atomic charge-density will be used'
         ENDIF

         IF (INFO%INICHG==1 .OR.INFO%INICHG==2 .OR.INFO%INICHG==3) THEN
            IF (IO%IU6>=0) WRITE(TIU6,*)'initial charge density was supplied:'
         ENDIF

         IF (INFO%INICHG==1) THEN

!PK Only first k-point group reads charge density
!PK Broadcast results outside of READCH due to awkward error logic in READCH
#ifdef MPI
            IF (COMM_KINTER%NODE_ME.EQ.1 ) THEN
#endif
              CALL READCH(GRIDC, INFO%LOVERL, T_INFO, CHTOT, RHOLM, INFO%INICHG, WDES%NCDIJ, &
                LATT_CUR, P, CSTRF(1,1), 18, IO%IU0, INFO%LHALF_INIT)

              CALL READTAU(XC_DATA%LTAU,GRIDC,T_INFO,LATT_CUR,WDES,IO,KINEDEN)
#ifdef MPI
            END IF
#endif
            CALLMPI( M_bcast_i( COMM_KINTER, INFO%INICHG, 1))
            IF (INFO%INICHG.NE.0) THEN
               CALLMPI( M_bcast_z( COMM_KINTER, CHTOT, GRIDC%MPLWV*WDES%NCDIJ))
               CALLMPI( M_bcast_d( COMM_KINTER, RHOLM, SIZE(RHOLM)) )
            END IF
            IF (XC_DATA%LTAU) THEN
              CALLMPI( M_bcast_i(COMM_KINTER,IO%INITAU,1) )
              IF (IO%INITAU.NE.0) THEN
                CALLMPI( M_bcast_z( COMM_KINTER, KINEDEN%TAU, GRIDC%MPLWV*WDES%NCDIJ))
              ENDIF
            ENDIF

            IF (INFO%ICHARG>10 .AND. INFO%INICHG==0) THEN
               CALL vtutor%error("ERROR: charge density could not be read from file CHGCAR for ICHARG>10")
            ENDIF
        ! error on reading CHGCAR, set INFO%INICHG to 2
           IF (INFO%INICHG==0)  INFO%INICHG=2
        ! no magnetization density set it according to MAGMOM
           IF (INFO%INICHG==-1) THEN
              IF (WDES%NCDIJ>1) &
                   CALL MRHOATO(.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CHTOT(1,2),WDES%NCDIJ-1)
              INFO%INICHG=1
              IF( TIU0 >=0) WRITE(TIU0,*)'magnetization density of overlapping atoms calculated'
           ENDIF

         ENDIF

         IF (INFO%INICHG==1) THEN
         ELSE IF (INFO%INICHG==2 .OR.INFO%INICHG==3) THEN
            IF(INFO%TURBO==0)THEN
               CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
            ELSE
               CALL RHOATO_PARTICLE_MESH(.FALSE.,.FALSE.,GRIDC,LATT_CUR,T_INFO,INFO,P,CHTOT,IO%IU6)
            ENDIF
            IF (WDES%NCDIJ>1) &
                 CALL MRHOATO(.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CHTOT(1,2),WDES%NCDIJ-1)

            IF (IO%IU6>=0) WRITE(TIU6,*)'charge density of overlapping atoms calculated'
         ELSE IF (INFO%INICHG==4) THEN
            IF (IO%IU6>=0) WRITE(TIU6,*)'potential read from file POT'
            INFO%INICHG=4
         ELSE
            INFO%INICHG=0
         ENDIF

         IF (INFO%INICHG==1 .OR.INFO%INICHG==2 .OR.INFO%INICHG==3) THEN
            DO I=1,WDES%NCDIJ
               RHOTOT(I) =RHO0(GRIDC, CHTOT(1,I))
            ENDDO
            IF(IO%IU6>=0)  WRITE(TIU6,200) RHOTOT(1:WDES%NCDIJ)
 200        FORMAT(' number of electron ',F15.7,' magnetization ',3F15.7)
         ENDIF

      ! set the partial core density
         DENCOR=0
         IF (INFO%LCORE) CALL RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENCOR,IO%IU6)
         IF (XC_DATA%LTAU) CALL TAUPAR(GRIDC,T_INFO,LATT_CUR%B,LATT_CUR%OMEGA,P,CSTRF,KINEDEN%TAUC)

         IF (INFO%INIWAV==2) THEN
            CALL vtutor%error("ERROR: this version does not support INIWAV=2")
         ENDIF

         IF (IO%IU6>=0) THEN
           IF (INFO%INICHG==0 .OR. INFO%INICHG==4 .OR. (.NOT.INFO%LCHCOS.AND. INFO%NELMDL==0) ) THEN
             WRITE(TIU6,*)'charge density for first step will be calculated', &
              ' from the start-wavefunctions'
           ELSE
             WRITE(TIU6,*)'keeping initial charge density in first step'
           ENDIF
           WRITE(TIU6,130)
         ENDIF

      ! add Gaussian "charge-transfer" charges, if required
         CALL RHOADD_GAUSSIANS(T_INFO,LATT_CUR,P,GRIDC,NCDIJ,CHTOT,CSTRF)
         CALL RHOADD_GAUSSIANS_LIST(LATT_CUR,GRIDC,NCDIJ,CHTOT)

         DWRITE0 'atomic charge done'
      ENDIF
!========================subroutine SPHER ==============================
! RSPHER calculates the real space projection operators
!    (VOLUME)^.5 Y(L,M)  VNLR(L) EXP(-i r k)
! subroutine SPHER calculates the nonlocal pseudopotential
! multiplied by the spherical harmonics and (1/VOLUME)^.5:
!    1/(VOLUME)^.5 Y(L,M)  VNL(L)
! (routine must be called if the size of the unit cell is changed)
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (INFO%LREAL) THEN
            CALL RSPHER(GRID,NONLR_S,LATT_CUR)

            INDMAX=0
#ifndef _OPENMP
            DO NI=1,T_INFO%NIONS
               INDMAX=MAX(INDMAX,NONLR_S%NLIMAX(NI))
            ENDDO
#else
            INDMAX=MAXVAL(NONLR_S%NLIMAX)
#endif
            IF (IO%IU6>=0) &
            WRITE(TIU6,*)'Maximum index for non-local projection operator ',INDMAX
         ELSE
            CALL SPHER(GRID, NONL_S, P, WDES, LATT_CUR,  1)
            CALL PHASE(WDES,NONL_S,0)
         ENDIF

         DWRITE0 'non local setup done'
      ENDIF
!=======================================================================
! set the coefficients of the plane wave basis states to zero before
! initialising the wavefunctions by calling WFINIT (usually random)
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (INFO%ISTART<=0) THEN
           IF (.NOT.INFO%LHALF_INIT) THEN
             W%CPTWFP=0
             CALL WFINIT(WDES, W, INFO%ENINI)
           ENDIF
         IF (.NOT.INFO%LHALF_INIT .AND. INFO%INIWAV==1 .AND. INFO%NELMDL==0 .AND. INFO%INICHG/=0) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*) 'WARNING: random wavefunctions but no delay for ', &
                           'mixing, default for NELMDL'
           INFO%NELMDL=-5
           IF (INFO%LRMM .AND. .NOT. INFO%LDAVID) INFO%NELMDL=-12
         ENDIF

      ! initialize the occupancies, in the spin polarized case
      ! we try to set up and down occupancies individually
         ELEKTR=INFO%NELECT
         CALL DENINI(W%FERTOT(:,:,1),WDES%NB_TOT,KPOINTS%NKPTS,ELEKTR,WDES%LNONCOLLINEAR)

         IF (WDES%ISPIN==2) THEN
           RHOMAG=RHO0(GRIDC,CHTOT(1,2))
           ELEKTR=INFO%NELECT+RHOMAG
           CALL DENINI(W%FERTOT(:,:,1),WDES%NB_TOT,KPOINTS%NKPTS,ELEKTR,WDES%LNONCOLLINEAR)
           ELEKTR=INFO%NELECT-RHOMAG
           CALL DENINI(W%FERTOT(:,:,2),WDES%NB_TOT,KPOINTS%NKPTS,ELEKTR,WDES%LNONCOLLINEAR)
         ENDIF

         ENDIF
         DWRITE0 'wavefunctions initialized'
      ENDIF
!=======================================================================
! INFO%LONESW initialize W_F%CELEN fermi-weights and augmentation charge
!=======================================================================
      E%EENTROPY=0

      IF(KPOINTS%ISMEAR == -2) THEN
        IF(INFO%LDIAG .OR. .NOT. INFO%LONESW) THEN
         CALL vtutor%write(isAlert, state_order)
        ENDIF
      ENDIF

      IF(LDO_AB_INITIO) THEN
         IF (INFO%LONESW .OR. (KPOINTS%ISMEAR/=-2 .AND. INFO%IALGO==3) ) THEN
            IF (INFO%LONESW) W_F%CELTOT = W%CELTOT
            CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
                  INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
                  NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
         ENDIF
      ENDIF
!=======================================================================
!  read Fermi-weigths from INCAR if supplied
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         CALL OPEN_INCAR_IF_FOUND(IO%IU5)

         CALL PROCESS_INCAR(.FALSE., IO%IU0, IO%IU5, 'FERWE', W%FERTOT, KPOINTS%NKPTS*WDES%NB_TOT, IERR, WRITEXMLINCAR)

         IF (WDES%ISPIN==2) THEN
            CALL PROCESS_INCAR(.FALSE., IO%IU0, IO%IU5, 'FERDO', W%FERTOT(:,:,INFO%ISPIN), KPOINTS%NKPTS*WDES%NB_TOT, IERR, WRITEXMLINCAR)
         ENDIF

         CALL CLOSE_INCAR_IF_FOUND(IO%IU5)
      ENDIF

!=======================================================================
! write out STM file if wavefunctions exist
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (STM(1) < STM(2) .AND. STM(3) /= 0 .AND. STM(4) < 0 .AND. INFO%ISTART == 1) THEN
            IF (STM(7) == 0) THEN
               CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
                    INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
                    NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
            ELSEIF (STM(7) /= 0) THEN
               EFERMI = STM(7)
               do_io WRITE(*,*)'-------------  VASP_IETS  --------------'
               do_io write(*,*) 'Efermi FROM INPUT= ',EFERMI
               do_io WRITE(*,*)'STM(7)=',STM(7)
               do_io WRITE(*,*)'-------------  VASP_IETS  --------------'
            ELSE
               CALL vtutor%error("------------- VASP_IETS -------------- \n ERROR WITH STM(7) STOP NOW \n &
                    &STM(7)= " // str(STM(7)))
            ENDIF
            do_io WRITE(*,*) 'Writing STM file'
            CALL  WRT_STM_FILE(GRID, WDES, WUP, WDW, EFERMI, LATT_CUR, &
                  STM, T_INFO)
            CALL vtutor%error("STM File written, ...")
         ENDIF
      ENDIF

!=======================================================================
! calculate the projections of the wavefunctions onto the projection
! operators using real-space projection scheme or reciprocal scheme
! then perform an orthogonalisation of the wavefunctions
!=======================================================================
      IF(LDO_AB_INITIO) THEN
      ! first call SETDIJ to set the array CQIJ
         CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                     LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

         CALL CALC_JPAW_HAMIL(T_INFO, P)
      ! set vector potential
         CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, &
                    LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)
      ! first call to SETDIJ_AVEC only sets phase twisted projectors
         CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                    LMDIM, CDIJ, HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX )

         DWRITE0 'setdij done'

         IF (IO%IU6>=0) &
         WRITE(TIU6,*)'Maximum index for augmentation-charges ', &
                        IRDMAA,'(set IRDMAX)'
         IF (WDES%LORBITALREAL) CALL WVREAL_PRECISE(W)

         CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
         DWRITE0 'proall done'

         CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)

         CALL REDIS_PW_OVER_BANDS(WDES, W)
         DWRITE0 'orthch done'

         IF (IO%IU6>=0) WRITE(TIU6,130)

      ! strictly required for HF, since SET_CHARGE fails for ISYM=3
         CALL KPAR_SYNC_ALL(WDES,W)
      ENDIF
!=======================================================================
! partial band decomposed chargedensities PARDENS
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (IO%WRT_DENSITY%LPARD) THEN
              CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
                   INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
                   NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
              CALL PARCHG(W,WUP,WDW,WDES,CHDEN,CHTOT,CRHODE,INFO,GRID, &
                   GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                   LATT_CUR,P,T_INFO,SOFT_TO_C,SYMM,IO, &
                   DYN,EFERMI,LMDIM,IRDMAX,NIOND)
         ENDIF
      ENDIF
!=======================================================================
! calculate initial chargedensity if
! ) we do not have any chargedensity
! ) we use a selfconsistent minimization scheme and do not have a delay
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (INFO%INICHG==0 .OR. INFO%INICHG==4 .OR. (.NOT.INFO%LCHCOS .AND. INFO%NELMDL==0)  ) THEN
            IF (IO%IU6>=0) &
                 WRITE(TIU6,*)'initial charge from wavefunction'

            IF (IO%IU0>=0) &
                 WRITE(TIU0,*)'initial charge from wavefunction'

            CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
                 GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
                 LATT_CUR, P, SYMM, T_INFO, &
                 CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

            CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
                 XC_DATA%LTAU,T_INFO%NIONS,W,WDES,KINEDEN)
         ENDIF
      ENDIF

!=======================================================================
! initialise the predictor for the wavefunction with
! the first available ionic positions
! PRED%INIPRE=2 continue with existing file
! PRED%INIPRE=1 initialise
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (INFO%ISTART==3) THEN
           PRED%INIPRE=2
         ELSE
           PRED%INIPRE=1
         ENDIF

         IF (DYN%IBRION/=-1 .AND. PRED%IWAVPR >= 12 ) THEN
           CALL WAVPRE_NOIO(GRIDC,P,XC_DATA,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
              CHTOT,RHOLM,N_MIX_PAW, CSTRF, KINEDEN, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
         ELSE IF (DYN%IBRION/=-1 .AND. PRED%IWAVPR >= 2 .AND. PRED%IWAVPR <10 )  THEN
           CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
              CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
         ENDIF
         DWRITE0 "wavpre is ok"
      ENDIF

!=======================================================================
! if INFO%IALGO < 0 make some performance tests
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (INFO%IALGO<0) THEN
            CALL PERFORMANCE_TEST
            GOTO 5100
         ENDIF
      ENDIF

!======================== SUBROUTINE FEWALD ============================
! calculate ewald energy forces, and stress
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         CALL START_TIMING("G")
         IF(INFO%TURBO==0) THEN
            CALL FEWALD_DATA_WRITER(GRIDC,LATT_CUR,T_INFO,HAMILTONIAN%COULOMB_POT,IO%IU6)
            CALL FEWALD(GRIDC,LATT_CUR,P,T_INFO,HAMILTONIAN%COULOMB_POT,E%TEWEN,EWIFOR,EWSIF)
         ELSE
            ALLOCATE(CWORK1(GRIDC%MPLWV))
            CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK1,E%PSCENC,E%TEWEN,EWIFOR,HAMILTONIAN%COULOMB_POT)
            DEALLOCATE(CWORK1)
         ENDIF
         CALL STOP_TIMING("G",TIU6,'FEWALD')
      ENDIF
!=======================================================================
! Set INFO%LPOTOK to false: this requires  a recalculation of the
! total local potential
!=======================================================================
      INFO%LPOTOK=.FALSE.
      TOTENG=0
      INFO%LSTOP=.FALSE.
      INFO%LSOFT=.FALSE.
      INFO%LHARD=.FALSE.

      IF(LDO_AB_INITIO) THEN
         IF (INFO%INICHG==4) THEN
            INQUIRE( FILE=DIR_APP(1:DIR_LEN)//'POT', EXIST=LEXIST_POT )

            ! prioritize legacy file 
            IF( LEXIST_POT ) THEN
               ! read legacy file POT
               IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE=DIR_APP(1:DIR_LEN)//'POT',STATUS='UNKNOWN')
               REWIND IO%IUVTOT
               CALL READPOT(GRIDC, T_INFO, CVTOT, INFO%INICHG, WDES, &
                    LATT_CUR, P, LMDIM, CDIJ, N_MIX_PAW, IO%IUVTOT, IO%IU0)
               IF (IO%LOPEN) CLOSE(IO%IUVTOT)
            ELSE
#ifndef VASP_HDF5
               CALL vtutor%error("POT file required to continue!")
#else
               ! vaspwave.h5 required to continue 
               IF( .NOT. VASPWAVE_FOUND ) THEN
                   CALL vtutor%error("vaspwave.h5 with a pot dataset required to continue!")
               ENDIF
           
               ! for legacy reasons the return code should be 4, so catch this here 
               ! without using the VH5 macro: 
               IH5ERR = VH5_READ_POT( IH5WAVEFILEID, GRIDC, T_INFO, CVTOT, WDES, &
                    LATT_CUR, P, LMDIM, CDIJ, N_MIX_PAW, IO)
               IF( IH5ERR == -1 ) THEN
                   INFO%INICHG = -1
               ELSEIF( IH5ERR /= 0 ) THEN
                   INFO%INICHG = 0
               ENDIF
#endif 
            ENDIF

            IF (INFO%INICHG==4) THEN
               DO ISP=1,WDES%NCDIJ
                  CALL FFT_RC_SCALE(CVTOT(1,ISP),CVTOT(1,ISP),GRIDC)
               ENDDO
               CALL SET_SV( GRID, GRIDC, GRID_SOFT, WDES%COMM_INTER, SOFT_TO_C, WDES%NCDIJ, SV, CVTOT)
               INFO%LPOTOK=.TRUE.
            ELSE
               CALL vtutor%error("reading POT failed can not continue")
            ENDIF
         ENDIF
      ENDIF
!=======================================================================
! some special calculation routines such as GW, MP2, MP2NO, ...
!=======================================================================
      IF(LDO_AB_INITIO) THEN
         IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. LBRACKETST) THEN
            CALL PSI_CORR_MAIN(P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI,  &
             & T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,MIX%AMIX,MIX%BMIX)
         ENDIF

         IF(FOURORBIT==1) THEN
            CALL TWOELECTRON4O_MAIN( &
                 P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                 T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID, LMDIM)
            CALL STOP_TIMING("G",IO%IU6,'TWOE4O')
         ENDIF
         IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. &
            & LBRACKETST .OR. FOURORBIT==1 ) THEN
            CALL vtutor%stopCode()
         ENDIF

         IF (L2E4W) THEN
            CALL DMFT_MATRIX_ELEMENTS( &
        &      WDES,W,KPOINTS,GRID,T_INFO,INFO,P, &
        &      NONL_S,NONLR_S,SYMM,LATT_CUR,CQIJ,LMDIM,IO)
         ! stop execution of VASP (wave functions have been messed up so let's quit)
            GOTO 5100
         ENDIF
#ifdef vasp6
         IF (LTIME_EVOLUTION) THEN
            CALL XML_TAG("calculation")
            CALL CALCULATE_TIME_PROPAGATION( &
                 KINEDEN, HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
                 GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,&
                 CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
                 T_INFO, DYN, INFO, XC_DATA, IO, KPOINTS, SYMM,  &
                 LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, CHDEN, SV, &
                 EFERMI, NEDOS, DOS, DOSI , IEPSILON, &
                 NBANDSO, NBANDSV, OMEGAMAX, SHIFT, INFO%NELM, &
                 LHARTREE, LADDER, LTRIPLET, LFXC, LGWLF, ANTIRES)
            CALL XML_CLOSE_TAG("calculation")
         ! stop execution of VASP (wave functions have been messed up so let's quit)
            GOTO 5100
         ENDIF
#endif
         IF (AFQMC_SET % ACTIVE) THEN
            CALL AFQMC_PROPAGATION(AFQMC_SET, KINEDEN, HAMILTONIAN, P, W, WDES, NONLR_S, NONL_S, &
               LATT_CUR, T_INFO, INFO, IO, GRID, GRID_SOFT, GRIDC, GRIDUS, C_TO_US, &
               SOFT_TO_C, SYMM, CHTOT, CVTOT, CSTRF, CHDEN, DENCOR, SV, CDIJ, CQIJ, &
               CRHODE, RHOLM, LMDIM, IRDMAX, N_MIX_PAW, KPOINTS)
            GOTO 5100
         END IF

         IF (LCHI) THEN
            IF (.NOT. LRPAFORCE) THEN
#ifdef MinimaxTest
               CALL MM_INT_TEST
               GOTO 5100
#endif
               CALL DO_RPA
               INFO%LSTOP=.TRUE.   ! stop ionic loop immediately
            ENDIF
         ENDIF
      ENDIF


!***********************************************************************
!***********************************************************************
!
! ++++++++++++ do 1000 n=1,required number of timesteps ++++++++++++++++
!
! this is the main loop of the program during which one complete step of
! the electron dynamics and ionic movements is performed.
! NSTEP           loop counter for ionic movement
! N               loop counter for self-consistent loop
! INFO%NELM            number of electronic movement loops
!***********************************************************************
!***********************************************************************

      IF (IO%IU0>=0) WRITE(TIU0,*)'entering main loop'

      NSTEP = 0
      IF (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD) THEN
         IF(LBLOCK_HELP>=2) THEN
            CALL START_TIMING("BLOCK+")
         END IF
      END IF
      CALL START_TIMING("LOOP+")

      CALL CONSTRAINED_M_INIT(T_INFO,GRIDC,LATT_CUR)
      CALL INIT_WRITER(P,T_INFO,WDES)

!!      PROFILING_START('ion_loop_main')
!=======================================================================
      ion: DO
      IF (INFO%LSTOP) EXIT ion
!=======================================================================

!  reset broyden mixing
      MIX%LRESET=.TRUE.
io_begin
      IF (((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) .OR. (MOD(NSTEP,MAX(DYN%NBLOCK,100))==0)) THEN
        IF (IO%LOPEN) CALL WFORCE(IO%IU6)            ! write of OUTCAR
        IF (IO%LOPEN .AND. NSTEP/=0) CALL WFORCE(17) ! write of OSZICAR
      ENDIF
io_end
!=======================================================================
! initialize pair-correlation funtion to zero
! also set TMEAN und SMEAN to zero
! TMEAN/SMEAN  is the mean temperature
!=======================================================================
      IF (MOD(NSTEP,DYN%NBLOCK*DYN%KBLOCK)==0) THEN
        IF (OUTPUT_MODE.GT.0) THEN
           CALL  PCDAT_CLEAR_PACO( PACO)
        ENDIF
        TMEAN=0
        TMEAN0=0
        SMEAN=0
        SMEAN0=0
        IF(LDO_AB_INITIO) THEN
           DDOSI=0._q
           DDOS =0._q
        ENDIF
      ENDIF

! embedding__
      IF(EXTPT_LEXTPOT().AND.LPAW) CALL EXTPT_CALC_VLM(LATT_CUR, GRIDC, T_INFO, P, LMDIM, WDES, IO%IU0)
! embedding__

      NSTEP = NSTEP + 1

#ifdef VASP_DEBUG
      !=============================================================================================
      ! Hooks for debugging ionic loop:
      !
      ! * Writing STOPCAR file:
      IF (NSTEP == DBG%WRITE_STOPCAR_AT_STEP) CALL WRITE_STOPCAR(COMM, IO%IU0)
      ! * Trigger internal soft stop:
      IF (NSTEP == DBG%REQUEST_EXIT_AT_STEP) THEN
         IF (DBG%REQUEST_EXIT_ON_RANK < 0 .OR. COMM%NODE_ME == DBG%REQUEST_EXIT_ON_RANK) THEN
            CALL TRIGGER_SOFT_STOP(INFO%LSOFT_TRIGGER, .TRUE.)
         ENDIF
      ENDIF
      !=============================================================================================
#endif

      IF (IO%IU6>=0) THEN
         IF (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD) THEN
            IF (NSTEP.EQ.1 .AND. MOD(NSTEP,LBLOCK_HELP)==0) WRITE(TIU6,140) NSTEP
         ELSE
            WRITE(TIU6,140) NSTEP
         ENDIF
      ENDIF
!    IF ((INFO%NELM==0.OR. ML_LMLFF) .AND. IO%IU6>=0) WRITE(TIU6,140) NSTEP

!***********************************************************************
!***********************************************************************
! this part performs the total energy minimisation
! INFO%NELM loops are made, if the difference between two steps
! is less then INFO%EDIFF the loop is aborted
!***********************************************************************
!***********************************************************************

! If ML_LMLFF=.TRUE., get energy, force and stress tensor from machine-learning library.
    IF (NSTEP==1) THEN
       IF (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD) THEN
          CALL START_TIMING("G")
       ENDIF
       IF (ML_LMLFF) THEN
#ifdef ML_AVAILABLE
          DO I = 1, ML_TOTNUM_INSTANCES
             IF (ML_SUPER_HANDLE_MAIN(1)%FF%LIB .AND. (ML_SUPER_HANDLE_MAIN(1)%FF%ISTART == 2)) THEN
#ifdef libvaspml
                 CALL VASPML%compute(T_INFO%POSION, &
                                     LATT_CUR%A, &
                                     LATT_CUR%B, &
                                     T_INFO%NIONS, &
                                     T_INFO%NTYP, &
                                     T_INFO%ITYP, &
                                     T_INFO%TYPE, &
                                     LATT_CUR%OMEGA,&
                                     TOTEN,&
                                     TIFOR,&
                                     TSIF)
#endif
             ELSE
                CALL ML_TO_VASP_MACHINE_LEARNING(ML_SUPER_HANDLE_MAIN(I), DYN, LATT_CUR, &
                     NIOND, NSTEP, T_INFO, TOTEN, TIFOR, TSIF, COMM%MPI_COMM)
                CALL TRIGGER_SOFT_STOP(INFO%LSOFT_TRIGGER, ML_REQUESTS_EXIT(ML_SUPER_HANDLE_MAIN(I)))
             ENDIF
             IF (SYMM%ISYM>0) THEN
                CALL FORSYM(TIFOR,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
                IF (DYN%ISIF/=0) CALL TSYM(TSIF,ISYMOP,NROTK,LATT_CUR%A)
             ENDIF
             PRESS=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))/3._q -DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
          ENDDO
          IF (.NOT. LDO_AB_INITIO_NEW .AND. INFO%NELMDL==0) INFO%NELMDL=-5
          LDO_AB_INITIO=LDO_AB_INITIO_NEW
          !CALL STOP_TIMING("G",TIU6,'MLFF',WRITE_FILES=MOD(NSTEP,LBLOCK_HELP)==0)
#endif
       ELSE IF (INFO%PLUGIN%FORCE_FIELD) THEN
#ifndef NO_MACHINE_LEARNING
          CALL PLUGINS_MACHINE_LEARNING(INFO, T_INFO, LATT_CUR, P, GRIDC, DYN, &
             SYMM, ISYMOP, NROTK, LDO_AB_INITIO, TOTEN, TIFOR, TSIF, PRESS)
#endif
       ENDIF
       IF (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD) THEN
          CALL STOP_TIMING("G",TIU6,'MLFF',WRITE_FILES=MOD(NSTEP,LBLOCK_HELP)==0)
       ENDIF
    ENDIF

    ! ML force and stress writing is only done here if not running ML_MODE = delta
    IF (LDELTADFT) THEN
       TOTEN_ML = TOTEN
       TIFOR_ML = TIFOR
       TSIF_ML  = TSIF
    ELSE
       IF (ML_LMLFF) THEN
#ifdef ML_AVAILABLE
          DO I = 1, ML_TOTNUM_INSTANCES
             IF (.NOT.ML_SUPER_HANDLE_MAIN(I)%SKIP_ML_TO_VASP) THEN 
                CALL ML_TO_VASP_WRITING( DYN, GRIDC, IO, LATT_CUR, .TRUE., &
                     (NSTEP<=NSW_ACT.OR.IO%NWRITE>=1), NSTEP, T_INFO, TOTEN, TIFOR, TSIF, .FALSE. )
             ENDIF
          ENDDO
#endif
       ELSE IF (INFO%PLUGIN%FORCE_FIELD) THEN
#ifndef NO_MACHINE_LEARNING
          CALL ML_TO_VASP_WRITING( DYN, GRIDC, IO, LATT_CUR, .TRUE., &
               (NSTEP<=NSW_ACT.OR.IO%NWRITE>=1), NSTEP, T_INFO, TOTEN, TIFOR, TSIF, .FALSE. )
#endif
       ENDIF
    ENDIF


! Here we need to check if any of the cores has LDO_AB_INITIO set. If yes, then
! at every step a check needs to be made for the STOPCAR file.
    ALL_LDO_AB_INITIO_AMONG_IMAGES = LDO_AB_INITIO
    CALLMPI(M_and(COMM_WORLD,ALL_LDO_AB_INITIO_AMONG_IMAGES,1))
! Here we check the reverse
    ANY_LDO_AB_INITIO_AMONG_IMAGES = LDO_AB_INITIO
    CALLMPI(M_or(COMM_WORLD,ANY_LDO_AB_INITIO_AMONG_IMAGES,1))
! If both are true then some cores have LDO_AB_INITIO at true and some not
! (DFT vs MLFF). In that case a warning is necessary, since the DFT should
! be put on image 1, because only image 1 can write to standard out.
    IF (.NOT.ALL_LDO_AB_INITIO_AMONG_IMAGES .AND. ANY_LDO_AB_INITIO_AMONG_IMAGES .AND. NSTEP.EQ.1) THEN
       CALL vtutor%alert("Please note that you are combining DFT and non-DFT&
               & calculations within your images. Only the first image will&
               & output to standard out. To ensure you see the electronic&
               & steps in the standard output, please assign the DFT calculation&
               & to the first image. Regardless, the output files (OUTCAR, OSZICAR,&
               & etc.) will still be fully generated as usual.")
    ENDIF

    IF ((.NOT. LJ_ONLY .AND. .NOT. LCHI) .AND. (LDO_AB_INITIO)) THEN

      CALL XML_TAG("calculation")

      CALL ELECTRONIC_OPTIMIZATION

      ! possibly orbitals have been updated
      ! force them to be real again be calling WVREAL_PRECISE
      ! costs some extra work, but it is ultimately safer to do it here
      IF (W%WDES%LORBITALREAL .AND. (INFO%LDIAG .OR. INFO%LEXACT_DIAG)) THEN
        CALL WVREAL_PRECISE(W)
        CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
        CALL STOP_TIMING("G",IO%IU6,"WVREAL",XMLTAG="wvreal")
        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
        CALL STOP_TIMING("G",IO%IU6,"ORTHCH",XMLTAG="orth")
      ENDIF

      IF (NSTEP==1) THEN
         IF (.NOT.INFO%LONESW.AND.(LPEAD_CALC_EPS().OR.LPEAD_NONZERO_EFIELD())) THEN
            CALL ALLOCW(WDES,W_F,WTMP,WTMP)
            CALL ALLOCW(WDES,W_G,WTMP,WTMP)
            DEALLOCATE(CHAM, CHF)
            ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
              CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
!           ! setup for scf subspace rotation
!           INFO%LONESW=.TRUE.
!           CALL SETUP_SUBROT_SCF(INFO,WDES,LATT_CUR,GRID,GRIDC,GRID_SOFT,SOFT_TO_C,IO%IU0,IO%IU5,IO%IU6)
!           INFO%LONESW=.FALSE.
         ENDIF
      ENDIF


      IF (NSTEP==1.OR.LPEAD_RELAX_STRUCTURE()) THEN
         CALL PEAD_ELMIN( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET_GL,DOS,DOSI,CHF,CHAM,ECONV, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)
      ENDIF

      IF (NSTEP==1.AND..NOT.LPEAD_RELAX_STRUCTURE()) THEN
         IF (.NOT.INFO%LONESW.AND.(LPEAD_CALC_EPS().OR.LPEAD_NONZERO_EFIELD())) THEN
            CALL DEALLOCW(W_F)
            CALL DEALLOCW(W_G)
            DEALLOCATE(CHAM, CHF)
            ALLOCATE(CHAM(1,1,1,1),CHF (1,1,1,1))
         ENDIF
      ENDIF
      

      IF ((.NOT.L2E4W) .AND. MLWFRUN()) THEN
         IF (IO%LORBIT==14) THEN
            CALL SPHPRO_OPTPROJ(P,T_INFO,W,WDES,LDIMP,LMDIMP,-1,KPOINTS%EMIN,KPOINTS%EMAX)
         ENDIF

         CALL MLWF_FREE(MLWF_GLOBAL)
         IF (NSTEP == 1) CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
         IF (LWANNIERINTERPOL.OR.LCRPAPLOT) THEN
            BLOCK
            TYPE (wannier_interpolator) :: WANN_INTERP
            CALL WANNIER_INTERPOL_SETUP_KPOINTS(WANN_INTERP,KPOINTS,LATT_CUR,SYMM,W)
            CALL WANNIER_INTERPOL_FROM_MLWF(WANN_INTERP,MLWF_GLOBAL,IO)
            CALL WANNIER_INTERPOL_SETUP_HAM(WANN_INTERP,W,IO)
            CALL WANNIER_INTERPOL_SETUP_CDER(WANN_INTERP,HAMILTONIAN,W,KPOINTS,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX,LMDIM,T_INFO,&
                                             NONLR_S,NONL_S,P,SV,CQIJ,CDIJ,LATT_INI,LATT_CUR,SYMM,INFO,IO)
            CALL WANNIER_INTERPOLATE(WANN_INTERP,W,KPOINTS,SYMM,LATT_INI,LATT_CUR,NEDOS,T_INFO,INFO,IO)
            CALL WANNIER_INTERPOL_FREE(WANN_INTERP)
            END BLOCK
         ENDIF
      ENDIF

      ! Prepare finite-difference operations
      IF (NSTEP == 1) THEN
         IF (ELPH_SETTINGS%DO_INIT()) THEN
            CALL ELPH_INTERPOL%INIT(ELPH_SETTINGS, MLWF_GLOBAL, W, CQIJ, CDIJ, P, PRIM_CELL, T_INFO, INFO, &
                 GRID, GRIDC, GRIDUS, C_TO_US, NONL_S, NONLR_S, IRDMAX, &
                 HAMILTONIAN, SV, SYMM, EXHF, IO)
         ENDIF
         IF (ELPH_DERIV_SETTINGS%GENERATE) THEN
            ALLOCATE(CVTOT_DERIV(GRIDC%MPLWV, WDES%NCDIJ, 3, PRIM_CELL%NUM_PRIM_ATOMS))
            ALLOCATE(CDIJ_DERIV(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ, 3, PRIM_CELL%NUM_PRIM_ATOMS))
            IF (XC_DATA%LTAU) ALLOCATE(MU_DERIV(GRIDC%MPLWV, WDES%NCDIJ, 3, PRIM_CELL%NUM_PRIM_ATOMS))
         ENDIF
      ENDIF
    ELSE IF (LCHI) THEN
#ifdef MinimaxTest
       CALL MM_INT_TEST
       GOTO 5100
#endif
       CALL DO_RPA
    END IF ! LJ_ONLY AND LDO_AB_INITIO

    IF(LDO_AB_INITIO) THEN
       ! embedding, symmetry has to be switched off,
       ! in case some states have been declared as active in embedding scheme,
       ! rotate them
       CALL ROTATE_ACTIVE_OCCUPIED_ORBITALS(WDES,W,NTARGET_STATES,&
           SYMM,LATT_CUR,LATT_INI,T_INFO,KPOINTS,GRID,NONL_S, P, DYN, IO)
    ENDIF

    !===============================================================================================
    ! This block is responsible for assessing the need for a 'soft stop', i.e., an early graceful
    ! exit from the ionic loop. This may be necessary because
    !
    ! * of a user request:
    !   via STOPCAR file containing LSTOP or LABORT set to .TRUE.)
    !
    ! * an internal routine requests it:
    !   via TRIGGER_SOFT_STOP (which updates INFO%LSOFT_TRIGGER),
    !   e.g. because continuing would produce unreliable results.
    !
    ! If a 'soft stop' was requested once (INFO%LSOFT = .TRUE.) it will not be checked again.
    IF (.NOT. INFO%LSOFT .AND. ( &
           ! If this is an ab-initio iteration
           ((.NOT. LJ_ONLY) .AND. ANY_LDO_AB_INITIO_AMONG_IMAGES) &
           ! or NSTEP is a is a multiple of MAX(NBLOCK, 100)
           .OR. (MOD(NSTEP, MAX(DYN%NBLOCK, 100)) == 0) &
           ! or ML_OUTBLOCK > 1 and the next iteration (NSTEP + 1) is a multiple of ML_OUTBLOCK
           .OR. (LBLOCK_HELP > 1 .AND. MOD(NSTEP + 1, LBLOCK_HELP) == 0))) THEN
       ! Check whether a STOPCAR file exists and contains LSTOP = .TRUE.
       ! No file check needed if stop was requested internally via LSOFT_TRIGGER.
       LTMP = INFO%LSOFT_TRIGGER
       IF (.NOT. LTMP) THEN
          CALL RDATAB(IO%LOPEN,'STOPCAR',99,'LSTOP','=','#',';','L', &
                      IDUM,RDUM,CDUM,LTMP,CHARAC,NCOUNT,1,IERR)
       ENDIF
       IF (.NOT. LTMP) THEN
          CALL RDATAB(IO%LOPEN,'STOPCAR',99,'LABORT','=','#',';','L', &
                      IDUM,RDUM,CDUM,LTMP,CHARAC,NCOUNT,1,IERR)
       ENDIF
       ! Check if one or more ranks request a soft stop.
       ITMP = 0
       IF (LTMP) ITMP = 1
       CALLMPI ( M_sum_i(COMM_WORLD, ITMP, 1))
       ! If so ...
       IF (ITMP > 0) THEN
         ! Signal that a 'soft stop' is scheduled, prevents checking again.
         INFO%LSOFT = .TRUE.
         IF (DYN%IBRION == 0) THEN
            ! In case of an MD run we reset NSW_ACT so the ionic loop stops after the current block
            ! defined by ML_OUTBLOCK.
            NSW_ACT = MIN(NSW_ACT, (NSTEP / LBLOCK_HELP + 1) * LBLOCK_HELP)
         ELSE
            ! All non-MD runs stop at the end of this loop iteration (actually beginning of next
            ! iteration)
            NSW_ACT = NSTEP
         ENDIF
         ! Report ...
         IF (LBLOCK_HELP == 1) THEN
            ! that we are quitting after this ionic step
            IF (IO%IU0>=0) WRITE(TIU0,*) 'soft stop encountered!  aborting job ...'
            IF (IO%IU6>=0) WRITE(TIU6,*) 'soft stop encountered!  aborting job ...'
         ELSE
            ! that we are quitting after step NSW_ACT
            IF (IO%IU0>=0) WRITE(TIU0,*) 'soft stop encountered!  aborting job ... after iteration '//str(NSW_ACT)
            IF (IO%IU6>=0) WRITE(TIU6,*) 'soft stop encountered!  aborting job ... after iteration '//str(NSW_ACT)
         ENDIF
       ENDIF
    ENDIF
    !===============================================================================================
     

    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
!=======================================================================
! do some check for occupation numbers (do we have enough bands]:
!=======================================================================
      IOCCUP=0
      IOCCVS=0
      SUMNEL=0._q
      DO ISP=1,WDES%ISPIN
      DO NN=1,KPOINTS%NKPTS
         IF (ABS(W%FERTOT(WDES%NB_TOT,NN,ISP))>1.E-2_q) THEN
            IOCCUP=IOCCUP+1
! total occupancy ('number of electrons in this band'):
            SUMNEL=SUMNEL+WDES%RSPIN*W%FERTOT(WDES%NB_TOT,NN,ISP)*KPOINTS%WTKPT(NN)
            IF (IOCCUP<NTUTOR) THEN
               ITUT(IOCCUP)=NN
               RTUT(IOCCUP)=W%FERTOT(WDES%NB_TOT,NN,ISP)
            ENDIF
         ENDIF
! count seperately 'seriously large occupations' ...
         IF (ABS(W%FERTOT(WDES%NB_TOT,NN,ISP))>2.E-1_q) IOCCVS=IOCCVS+1
      ENDDO
      ENDDO

      IF (((IOCCUP/=0).AND.(SUMNEL>1E-3_q)).OR.(IOCCVS/=0)) THEN
         IOCC=MIN(IOCCUP+1,NTUTOR)
         ITUT(IOCC)=WDES%NB_TOT
         RTUT(IOCC)=SUMNEL
         CALL vtutor%write(isAlert, HighestBandsOcccupied, newArgument(rval = RTUT))
      ELSE IF (IOCCUP/=0) THEN
! for less serious occupancies give just some 'good advice' ...
         IOCC=MIN(IOCCUP+1,NTUTOR)
         ITUT(IOCC)=WDES%NB_TOT
         RTUT(IOCC)=SUMNEL
         CALL vtutor%write(isAdvice, HighestBandsOcccupied, newArgument(rval = RTUT))
      ENDIF

! Okay when we arrive here first then we got it for the first ionic
! configuration and if more steps follow we might not need any further
! 'delay' switching on the selfconsistency? (if NELMDL<0 switch off!).
      IF (INFO%NELMDL<0) INFO%NELMDL=0

! calculation partial DOS
      IF (IO%LORBIT>=10) THEN
         CALL SPHPRO_FAST( &
          GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
          INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,NSTEP, IO%LORBIT,PAR, &
          EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, SYMM)
      ENDIF
!***********************************************************************
!***********************************************************************
! Now perform the ion movements:
!***********************************************************************
!***********************************************************************

!====================== FORCES+STRESS ==================================
!
!=======================================================================
      NORDER=0
      IF (KPOINTS%ISMEAR>=0) NORDER=KPOINTS%ISMEAR

 7261 FORMAT(/ &
     &        A/ &
     &        '  ---------------------------------------------------'/ &
     &        '  free  energy   TOTEN  = ',F18.8,' eV'// &
     &        '  energy  without entropy=',F18.8, &
     &        '  energy(sigma->0) =',F18.8)

      ! no forces for OEP methods (EXXOEP/=0)
      ! no forces for IALGO=1-4
      ! no forces if potential was read in INFO%INICHG==4
      ! no forces if inner eigenvalue solver are used
      ! no forces for RPA (has already been done)
      IF ( EXXOEP==0 .AND.  (INFO%LONESW .OR. INFO%LDIAG .OR. INFO%IALGO>4) &
           .AND. INFO%INICHG/=4 .AND. INFO%IHARMONIC==0 .AND. (.NOT. LCHI)) THEN
        CALL FORCE_AND_STRESS( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
          T_INFO,T_INFO,DYN,INFO,XC_DATA,IO,MIX,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,&
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET_GL, &
          LMDIM, IRDMAX, (NSTEP==1 .OR.NSTEP==DYN%NSW).OR.IO%NWRITE>=1, &
          DYN%ISIF/=0, (DYN%ISIF/=0 .AND. HAMILTONIAN%COULOMB_POT%KERNEL_TRUNCATE%ACTIVE .eqv. .FALSE.), &
          .TRUE., XCSIF, EWSIF, TSIF, EWIFOR, TIFOR, PRESS, TOTEN, KPOINTS)
      ENDIF

!     trigger calculation of ensemble energies in next call to POTLOK
#ifdef libbeef
      IF(LBEEFENS) LBEEFCALCBASIS = .TRUE.
#endif

      IF (IO%IU6>=0) THEN
         WRITE(TIU6,130)
         WRITE(TIU6,7261) '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)', &
            TOTEN,TOTEN-E%EENTROPY,TOTEN-E%EENTROPY/(2+NORDER)
      ENDIF

      ! Deltawriting is possibly done here
      IF (LDELTADFT) THEN
         IF (ML_LMLFF) THEN
#ifdef ML_AVAILABLE
            TOTEN = TOTEN + TOTEN_ML
            TIFOR = TIFOR + TIFOR_ML
            TSIF  = TSIF + TSIF_ML
            DO I = 1, ML_TOTNUM_INSTANCES
               IF (.NOT.ML_SUPER_HANDLE_MAIN(I)%SKIP_ML_TO_VASP) THEN 
                  CALL ML_TO_VASP_WRITING( DYN, GRIDC, IO, LATT_CUR, .TRUE., &
                       (NSTEP<=NSW_ACT.OR.IO%NWRITE>=1), NSTEP, T_INFO, TOTEN, TIFOR, TSIF, LDELTADFT )
               ENDIF
            ENDDO
#endif
         ELSE IF (INFO%PLUGIN%FORCE_FIELD) THEN
#ifndef NO_MACHINE_LEARNING
            TOTEN = TOTEN + TOTEN_ML
            TIFOR = TIFOR + TIFOR_ML
            TSIF  = TSIF + TSIF_ML
            CALL ML_TO_VASP_WRITING( DYN, GRIDC, IO, LATT_CUR, .TRUE., &
                 (NSTEP<=NSW_ACT.OR.IO%NWRITE>=1), NSTEP, T_INFO, TOTEN, TIFOR, TSIF, LDELTADFT )
#endif
         ENDIF
      ENDIF

    ELSE
      NORDER = 1
    ENDIF ! LJ_ONLY and LDO_AB_INITIO

    IF (OUTPUT_MODE .GT.0 .AND. MOD(NSTEP,LBLOCK_HELP)==0) THEN
       CALL XML_TAG("structure")
       CALL XML_CRYSTAL(LATT_CUR%A, LATT_CUR%B, LATT_CUR%OMEGA)
       CALL XML_POSITIONS(T_INFO%NIONS, DYN%POSION)
       CALL XML_CLOSE_TAG("structure")
       CALL XML_FORCES(T_INFO%NIONS, TIFOR)
       IF (DYN%ISIF/=0) THEN
          CALL XML_STRESS(TSIF*EVTOJ*1E22_q/LATT_CUR%OMEGA)
       ENDIF
    ENDIF
#ifdef VASP_HDF5
    IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
       CALL VH5_WRITE_STRESS_FORCES_HIST_STEP(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY, NSTEP/LBLOCK_HELP, &
                                              TSIF*EVTOJ*1E22_q/LATT_CUR%OMEGA, TIFOR)
       CALL VH5_WRITE_ENERGIES_HIST_STEP(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY, NSTEP/LBLOCK_HELP, &
                                         [TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER)])
    ENDIF

    IF (LCALC_ORBITAL_MOMENT().AND.WDES%LNONCOLLINEAR) THEN
      ALLOCATE(ORBITAL_MOMENT(LDIMP-1,T_INFO%NIONS,3))
      CALL GET_ORBITAL_MOMENT(WDES,ORBITAL_MOMENT)
      CALL VH5_WRITE_MOMENTS_HIST_STEP(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,NSTEP,ORBITAL_MOMENT,"orbital")
      DEALLOCATE(ORBITAL_MOMENT)
    END IF
#endif
! check the consistency of forces and total energy
    CALL CHECK(T_INFO%NIONS,DYN%POSION,TIFOR,EWIFOR,TOTEN,E%TEWEN,LATT_CUR%A,IO%IU6,MOD(NSTEP,LBLOCK_HELP)==0)

    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
      IF (DYN%PSTRESS/=0) THEN
         TOTEN=TOTEN+DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA

         IF (IO%IU6>=0) &
              WRITE(TIU6,7264) TOTEN,DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
 7264    FORMAT ('  enthalpy is  TOTEN    = ',F18.8,' eV   P V=',F18.8/)

      ELSE
         IF (IO%IU6>=0) WRITE(TIU6,*)
      ENDIF
      IF (IO%IU6>=0) WRITE(IO%IU6,130)

!-----------------------------------------------------------------------
! we have maybe to update CVTOT here (which might have been destroyed
!  by the force routines)
!-----------------------------------------------------------------------
      IF (LVCADER) THEN
         ! derivative with respect to VCA parameters
         CALL VCA_DER( &
              HAMILTONIAN,KINEDEN, &
              P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
              T_INFO,INFO,XC_DATA,IO,KPOINTS,GRID,GRID_SOFT, &
              GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,SYMM, &
              CHTOT,DENCOR,CVTOT,CSTRF, &
              CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
              CHDEN,SV,DOS,DOSI,CHAM, &
              LMDIM,IRDMAX,NEDOS)
          CALL STOP_TIMING("G",TIU6,'VCADER')
      ENDIF

      IF (.NOT.INFO%LPOTOK) THEN

      CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,  &
           INFO,XC_DATA,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF,&
           HAMILTONIAN%COULOMB_POT)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)

      CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ)

      CALL SET_DD_PAW(WDES, P, XC_DATA, T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ, RHOLM, CRHODE, &
         E, LASPH =INFO%LASPH, LCOREL=.FALSE.)

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,XC_DATA,IO%IU6)

      CALL STOP_TIMING("G",TIU6,'POTLOK')

      ENDIF

      CALL EGRAD_EFG_PW_HAR_ONLY(T_INFO,LATT_CUR,WDES,GRIDC,P,CSTRF,CHTOT(:,1))
      CALL EGRAD_WRITE_EFG(T_INFO,WDES,IO)

      IF (WDES%NCDIJ>1) THEN
         CALL HYPERFINE_PW(T_INFO,LATT_CUR,GRIDC,CHTOT(:,2))
         CALL HYPERFINE_WRITE_TENSORS(T_INFO,WDES,IO)
      ENDIF

      IF (WDES%ISPIN==2) THEN
         CALL DMATRIX_CALC(WDES,W,T_INFO,LATT_CUR,P,NONL_S,NONLR_S,GRID,GRIDC,CHTOT,LMDIM,CRHODE,IO%IU6,IO%IU0)
      ENDIF
      

!=======================================================================
!
!  electronic polarisation
!
!=======================================================================
      IF (LBERRY) THEN
      io_begin
         IF (IO%IU6>=0) THEN
            WRITE(TIU6,7230) IGPAR,NPPSTR
         ENDIF
      io_end
         CALL BERRY(WDES,W,GRID,GRIDC,GRIDUS,C_TO_US,LATT_CUR,KPOINTS, &
           IGPAR,NPPSTR,P,T_INFO,LMDIM,CQIJ,IRDMAX,INFO%LOVERL,IO%IU6, &
           DIP%POSCEN)
      ENDIF

 7230 FORMAT(/' Berry-Phase calculation of electronic polarization'// &
     &        '   IGPAR = ',I1,'   NPPSTR = ',I4/)
      IF (IO%IU6>=0) WRITE(TIU6,130)

      IF (LWANNIER()) THEN
         CALL LOCALIZE(W,WDES,GRID,GRIDC,GRIDUS,C_TO_US,LATT_CUR,T_INFO,P, &
        &     LMDIM,INFO%LOVERL,IRDMAX,IO%IU0,IO%IU6)
      ENDIF
      IF ( (DYN%IBRION /= -1 .AND. DYN%IBRION /=8 .AND. DYN%IBRION /=7) .AND. LEPSILON ) THEN
         ! if LEPSILON is set LR_SKELETON is called every NBLOCK*KBLOCK steps
         ! to calculate autocorrelation functions of the polarization
         IF (MOD(NSTEP-1,DYN%NBLOCK*DYN%KBLOCK)==0) THEN
            CALL LR_SKELETON( &
              KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
              T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
              GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
              CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
              CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
              CHDEN,SV,VDW_SET_GL,DOS,DOSI,CHAM, &
              DYN%IBRION,LMDIM,IRDMAX,NEDOS, &
              TOTEN,EFERMI, TIFOR, &
              .FALSE., LFAST=.TRUE. )
         ENDIF

         CALL STORE_DIPOLE_AND_AUTOCORRELATE( T_INFO, LATT_CUR, DYN%POTIM, DYN%NSW, NSTEP, DYN%TEBEG, IO)
      ENDIF
    ENDIF ! LJ_ONLY AND LDO_AB_INITIO



#if defined(ML_AVAILABLE)
#ifndef ML_LEARNING_NOT_ACTIVE
! Provide ab initio data to machine-learning library
    IF(ML_LMLFF) THEN
       !< We do not want to learn external forces, so they should be removed here if added in ab-initio step.
       IF (LDO_AB_INITIO) TIFOR = TIFOR - DYN%EFOR
       DO I = 1, ML_TOTNUM_INSTANCES
          IF (ML_SUPER_HANDLE_MAIN(I)%FF%LRUN_VASP_TO_ML) THEN
             CALL VASP_TO_ML_MACHINE_LEARNING( ML_SUPER_HANDLE_MAIN(I), DIR_APP,DIR_LEN, INFO, LATT_CUR, NSTEP, T_INFO, DYN, IO, GRIDC, NIOND, TOTEN, TIFOR, TSIF )
             IF (ML_SUPER_HANDLE_MAIN(I)%FF%IUPDATE_CRITERIA.EQ.1) THEN
                IF (ML_SUPER_HANDLE_MAIN(I)%FF%LGENFF) THEN
                   CALL ML_TO_VASP_MACHINE_LEARNING( ML_SUPER_HANDLE_MAIN(I), DYN, LATT_CUR, &
                        NIOND, NSTEP, T_INFO, TOTEN, TIFOR, TSIF, COMM%MPI_COMM )
                   CALL TRIGGER_SOFT_STOP(INFO%LSOFT_TRIGGER, ML_REQUESTS_EXIT(ML_SUPER_HANDLE_MAIN(I)))
                ENDIF
             ENDIF
             CALL STOP_TIMING("G",TIU6,'MLFFLR')
          ENDIF
       ENDDO
       !< Add external forces for the MD step
       TIFOR = TIFOR + DYN%EFOR
    ENDIF
#endif
#endif


    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
! add chain forces and constrain forces
!=======================================================================
      CALL VCA_FORCE(T_INFO, TIFOR)

      IF (DYN%IBRION /=5 .AND. DYN%IBRION /=6  .AND. DYN%IBRION /=7 ) &
      CALL SET_SELECTED_FORCES_ZERO(T_INFO,DYN%VEL,TIFOR,LATT_CUR)

      ! scale forces and total energy by scaling argument
      ! this allows  integration from ideal gas to liquid state
      ! alternatively if the DYNMATFULL file exists add forces from
      ! second order Hessian matrix (allows thermodynamic integration from
      ! harmonic model)
      ! maybe this should go over to the force.F file
      ! of course it is unique since it  changes radically energy
      CALL DYNMATFULL_ENERGY_FORCE(SCALEE, T_INFO%NIONS, DYN%POSION, TOTEN, TIFOR, LATT_CUR%A , IO%IU0)
      IF (SCALEE/=1) THEN
         EENTROPY=E%EENTROPY*SCALEE
         IF (IO%IU6>=0) WRITE(TIU6,7261) '  SCALED FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)', &
            TOTEN,TOTEN-EENTROPY,TOTEN-EENTROPY/(2+NORDER)
      ENDIF

      IF ( DYN%ISIF >= 3 ) THEN
         CALL CHAIN_STRESS( TSIF )
      END IF
      CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
           LATT_CUR%A,LATT_CUR%B,IO%IU6)

      CALL PARALLEL_TEMPERING(NSTEP,T_INFO%NIONS,DYN%POSION,DYN%VEL,TOTEN,TIFOR,DYN%TEBEG,DYN%TEEND, &
           LATT_CUR%A,LATT_CUR%B,IO%IU6)

      IF (DYN%IBRION /=5 .AND. DYN%IBRION /=6 .AND. DYN%IBRION /=7 ) &
         CALL SET_SELECTED_FORCES_ZERO(T_INFO,DYN%VEL,TIFOR,LATT_CUR)

      EKIN=0
      EKIN_LAT=0
      TEIN=0
      ES =0
      EPS=0

!     CYCLE ion ! if uncommented VASP iterates forever without ionic upd

      IF(LDO_AB_INITIO) THEN
         DO I=1,WDES%NCDIJ
            RHOTOT(I)=RHO0(GRIDC, CHTOT(1,I))
         END DO
      ENDIF ! LDO_AB_INITIO
    END IF ! LJ_ONLY

      IF (LJ_IS_ACTIVE()) THEN
        !FIXME: TOTEN should be declared explicitly
        TOTEN = 0.0_q
        CALL FORCE_AND_STRESS_LJ(GRIDC,XC_DATA,IO,WDES,P,LATT_CUR,DYN,T_INFO,TSIF,TIFOR,TOTEN)

        IF (DYN%IBRION /=5 .AND. DYN%IBRION /=6 .AND. DYN%IBRION /=7 ) &
          CALL SET_SELECTED_FORCES_ZERO(T_INFO,DYN%VEL,TIFOR,LATT_CUR)
      END IF

      DYN%AC=  LATT_CUR%A

    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
!
! write statement for static calculations with possibly extensive
! post processing (GW, dielectric response, ...)
!
!=======================================================================
      IF (DYN%IBRION/=0) THEN
         CALL XML_TAG("energy")
         IF(LDO_AB_INITIO) THEN
            EENTROPY=E%EENTROPY*SCALEE
         ELSE
            EENTROPY=0.0
         ENDIF
         CALL XML_ENERGY(TOTEN, TOTEN-EENTROPY, TOTEN-EENTROPY/(2+NORDER))
         CALL XML_CLOSE_TAG
      ENDIF

      IF (DYN%IBRION==-1 .OR. DYN%IBRION==5 .OR. DYN%IBRION==6 .OR. &
           DYN%IBRION==7.OR.DYN%IBRION==8) THEN
         IF(LDO_AB_INITIO) THEN
            EENTROPY=E%EENTROPY*SCALEE
         ELSE
            EENTROPY=0.0D0
         ENDIF
         io_begin
         IF(LDO_AB_INITIO) THEN
            WRITE(17 ,7281,ADVANCE='NO') NSTEP,TOTEN, &
                 TOTEN-EENTROPY/(2+NORDER),EENTROPY
            IF (IO%IU0>=0) &
                 WRITE(TIU0,7281,ADVANCE='NO') NSTEP,TOTEN, &
                 TOTEN-EENTROPY/(2+NORDER),EENTROPY

            IF ( WDES%NCDIJ>=2 ) THEN
              WRITE(17,77281) RHOTOT(2:WDES%NCDIJ)
              IF (IO%IU0>=0) WRITE(TIU0,77281) RHOTOT(2:WDES%NCDIJ)
            ELSE
              WRITE(17,*)
              IF (IO%IU0>=0) WRITE(TIU0,*)
            ENDIF
         ELSE
            WRITE(17 ,7281,ADVANCE='NO') NSTEP,TOTEN, &
                 TOTEN,0.0D0
            IF (IO%IU0>=0) &
                 WRITE(TIU0,7281,ADVANCE='NO') NSTEP,TOTEN, &
                 TOTEN,0.0D0

              WRITE(17,*)
              IF (IO%IU0>=0) WRITE(TIU0,*)
         ENDIF

         IF (M_CONSTRAINED()) CALL WRITE_CONSTRAINED_M(17,.TRUE.)

         io_end
      ENDIF
    END IF ! LJ_ONLY
!=======================================================================
! IBRION = 0  molecular dynamics
! ------------------------------
! ) calculate the accelerations in reduced units
!   the scaling is brain-damaging, so here it is in a more native form
!    convert dE/dr from EV/Angst to J/m   EVTOJ/1E-10
!    divide by mass                       1/T_INFO%POMASS*AMTOKG
!    multiply by timestep*2               (DYN%POTIM*1E-15)**2
! ) transform to direct mesh KARDIR
! ) integrate the equations of motion for the ions using the algorithm
!   of choice
! after this, the coordinates and lattice constants at which the calculations
! have been performed are stored in DYN%AC and  DYN%POSIOC
! the updated one go to T_INFO%POSION = DYN%POSION and LATT_CUR%A
!=======================================================================
ibrion: IF (DYN%IBRION==0 .OR. DYN%IBRION==4) THEN

        FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
        NI=1
        DO NT=1,T_INFO%NTYP
        DO NI=NI,T_INFO%NITYP(NT)+NI-1
          DYN%D2C(1,NI)=TIFOR(1,NI)*FACT/2/T_INFO%POMASS(NT)
          DYN%D2C(2,NI)=TIFOR(2,NI)*FACT/2/T_INFO%POMASS(NT)
          DYN%D2C(3,NI)=TIFOR(3,NI)*FACT/2/T_INFO%POMASS(NT)
        ENDDO
        ENDDO

        CALL KARDIR(T_INFO%NIONS,DYN%D2C,LATT_CUR%B)

        ISCALE=0
        IF (DYN%SMASS==-1 .AND. MOD(NSTEP-1,DYN%NBLOCK)==0 ) THEN
          ISCALE=1
        ENDIF

        DYN%TEMP=DYN%TEBEG+(DYN%TEEND-DYN%TEBEG)*NSTEP/ABS(DYN%NSW)

#ifdef tbdyn
        IF (IO%IU6>0) THEN
           IF (NSTEP == 1 .AND. MOD(NSTEP,LBLOCK_HELP)==0) THEN
             WRITE(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '        RANDOM_SEED = '
             DO i=1,K_SEED
               WRITE(IO%IU6,ADVANCE='NO',FMT='(X,I16)') SEED_INIT(i)
             ENDDO
             WRITE(IO%IU6,ADVANCE='NO',FMT='(/)')
           ENDIF

           IF (NSTEP == 1) THEN
             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '        RANDOM_SEED = '
             DO i=1,K_SEED
               WRITE(g_io%REPORT,ADVANCE='NO',FMT='(X,I16)') SEED_INIT(i)
             ENDDO
             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(/)')
           ENDIF
        ENDIF

        CALL STEP_tb(DYN,T_INFO,INFO,LATT_CUR,NHC,EKIN,EKIN_LAT,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,&
                          IO,TOTEN,g_io,WDES,TSIF,ISCALE,TEIN,TIFOR)


        CALL RANE_ION(RDUM,GET=SEED(1:K_SEED))

        IF (IO%IU6>0) THEN
           IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
              WRITE(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '        RANDOM_SEED = '
              DO i=1,K_SEED
                 WRITE(IO%IU6,ADVANCE='NO',FMT='(X,I16)') SEED(i)
              ENDDO
              WRITE(IO%IU6,ADVANCE='NO',FMT='(/)')
           ENDIF

           IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(/,3X,A22)')  '        RANDOM_SEED = '
             DO i=1,K_SEED
               WRITE(g_io%REPORT,ADVANCE='NO',FMT='(X,I16)') SEED(i)
             ENDDO
             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(/)')
           ENDIF
        ENDIF
#else
        IF ( .NOT. DYN&IBRION==4 ) THEN
           IF (.NOT. T_INFO%LSDYN) THEN
              CALL SYMVEL (T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS, &
              DYN%POSION,DYN%D2C,LATT_CUR%A,LATT_CUR%B)
           ENDIF

           CALL STEP(DYN%INIT,ISCALE,T_INFO%NIONS,LATT_CUR%A,LATT_CUR%ANORM,DYN%D2C,DYN%SMASS,DYN%POSION,DYN%POSIOC, &
                DYN%POTIM,T_INFO%POMASS,T_INFO%NTYP,T_INFO%ITYP,DYN%TEMP,DYN%VEL,DYN%D2,DYN%D3,DYN%SNOSE, &
                EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM, IO%IU6)

           TEIN = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM
        ENDIF
#endif
        CALL STOP_TIMING("G",TIU6,'IONSTEP',WRITE_FILES=MOD(NSTEP,LBLOCK_HELP)==0)

! sum energy of images along chain
        CALLMPI_C( sum_chain( TOTEN ))
        CALLMPI_C( sum_chain( EKIN ))
        CALLMPI_C( sum_chain( EKIN_LAT ))
        CALLMPI_C( sum_chain( ES ))
        CALLMPI_C( sum_chain( EPS ))

       ETOTAL=TOTEN+EKIN+ES+EPS


!  report  energy  of electrons + kinetic energy + nose-energy
        io_begin

        IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
           WRITE(TIU6,7260) TOTEN,EKIN-EKIN_LAT,EKIN_LAT,TEIN,ES,EPS,ETOTAL
#ifdef VASP_HDF5
           CALL VH5_WRITE_ENERGIES_HIST_STEP(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY, NSTEP/LBLOCK_HELP, &
                                             [TOTEN, EKIN-EKIN_LAT, EKIN_LAT, TEIN, ES, EPS, ETOTAL])
#endif
           CALL OFIELD_WRITE(TIU6)
           IF (LJ_IS_ACTIVE()) THEN
              WRITE(TIU6, 7262) TEIN * BOLKEV / LJ_EPSILON
           END IF
        ENDIF
        IF(LDO_AB_INITIO) THEN
           EENTROPY = E%EENTROPY*SCALEE
        ELSE
           EENTROPY = 0.0
        END IF
        IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
           CALL XML_TAG("energy")
           CALL XML_ENERGY(TOTEN, TOTEN-EENTROPY, TOTEN-EENTROPY/(2+NORDER))
           CALL XML_TAG_REAL("kinetic",EKIN-EKIN_LAT)
           CALL XML_TAG_REAL("lattice kinetic",EKIN_LAT)
           CALL XML_TAG_REAL("nosepot",ES)
           CALL XML_TAG_REAL("nosekinetic",EPS)
           CALL XML_TAG_REAL("total",ETOTAL)
           CALL XML_CLOSE_TAG
        ENDIF

! If LDO_AB_INITIO=.TRUE., you can use NORDER to compute the entropic
! contributions from the electrons.
        IF(LDO_AB_INITIO) THEN
           WRITE(17,7280,ADVANCE='NO') NSTEP,TEIN,ETOTAL,TOTEN, &
                TOTEN-EENTROPY/(2+NORDER),EKIN,ES,EPS
           IF (IO%IU0>=0) &
                WRITE(TIU0,7280,ADVANCE='NO')NSTEP,TEIN,ETOTAL,TOTEN, &
                TOTEN-EENTROPY/(2+NORDER),EKIN,ES,EPS
           IF (WDES%NCDIJ>=2) THEN
              WRITE(17,77280) RHOTOT(2:WDES%NCDIJ)
              IF (IO%IU0>=0) WRITE(TIU0,77280) RHOTOT(2:WDES%NCDIJ)
           ELSE
              WRITE(17,*)
              IF (IO%IU0>=0) WRITE(TIU0,*)
           ENDIF
! Otherwise, you do not need to take into account the entropic
! contributions.
        ELSE IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
           WRITE(17,7280,ADVANCE='NO') NSTEP,TEIN,ETOTAL,TOTEN, &
                TOTEN,EKIN,ES,EPS
           IF (IO%IU0>=0) &
                WRITE(TIU0,7280,ADVANCE='NO')NSTEP,TEIN,ETOTAL,TOTEN, &
                TOTEN,EKIN,ES,EPS
           WRITE(17,*)
           IF (IO%IU0>=0) WRITE(TIU0,*)
        ENDIF

        IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
            WRITE(TIU6,7270) DISMAX
        END IF         ! LJ_ONLY AND LDO_AB_INITIO
        io_end

7260  FORMAT(/ &
     &        '  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '% ion-electron   TOTEN  = ',F16.6,'  see above'/ &
     &        '  kinetic energy EKIN   = ',F16.6/ &
     &        '  kin. lattice  EKIN_LAT= ',F16.6, &
     &        '  (temperature',F8.2,' K)'/ &
     &        '  nose potential ES     = ',F16.6/ &
     &        '  nose kinetic   EPS    = ',F16.6/ &
     &        '  ---------------------------------------------------'/ &
     &        '  total energy   ETOTAL = ',F16.6,' eV'/)

7262  FORMAT('  reduced temperature T*= ',F16.6, ' eps/DOF'/)

 7270 FORMAT( '  maximum distance moved by ions :',E14.2/)

 7280 FORMAT(I7,' T= ',F6.0,' E= ',E14.8, &
     &   ' F= ',E14.8,' E0= ',E14.8,1X,' EK= ',E11.5, &
     &   ' SP= ',E8.2,' SK= ',E8.2)
77280 FORMAT(' mag=',3F11.3)

!=======================================================================
! DYN%IBRION =
! 1  quasi-Newton algorithm
! 2  conjugate gradient
! 3  quickmin
! 4  Beeman-Brooks-Karplus structure optimization done in STEP_tb
! 5  finite differences
! 6  finite differences with symmetry
! 7  linear response
!  IBRION ==5 finite differences
!=======================================================================
      ELSE IF (DYN%IBRION==5) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
        DYN%POSIOC=DYN%POSION
        CALL FINITE_DIFF( INFO%LSTOP, DYN%POTIM, T_INFO%NIONS, T_INFO%NTYP, &
             T_INFO%NITYP, T_INFO%POMASS, DYN%POSION, TIFOR, DYN%NFREE, &
             T_INFO%LSDYN,T_INFO%LSFOR, LATT_CUR%A, LATT_CUR%B,  &
             IO%IU6, IO%IU0 , IO%NWRITE, &
             FORCE_CONST, PHON_SETTINGS%DO_INIT())
        CALL LATTIC(LATT_CUR)
        ! we need to reinitialise the symmetry code at this point
        ! the number of k-points is changed on the fly
        IF (SYMM%ISYM>0) THEN
          CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,DYN%EFOR,T_INFO%LSFOR, &
             T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,NIOND, &
             SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
             SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#endif
          IF(LDO_AB_INITIO) CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          IF(LDO_AB_INITIO) CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          IF (INFO%LONESW.AND.LDO_AB_INITIO) THEN
             CALL DEALLOCW(W_F)
             CALL DEALLOCW(W_G)
             DEALLOCATE(CHAM, CHF)
             CALL ALLOCW(WDES,W_F,WTMP,WTMP)
             CALL ALLOCW(WDES,W_G,WTMP,WTMP)
             ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
               CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
          ENDIF
        ENDIF
       END IF !LJ_ONLY
      ELSE IF (DYN%IBRION==6) THEN ibrion
        IF (.NOT. LJ_ONLY) THEN
           BLOCK
             LOGICAL :: DO_ALL_DISP, RETURN_FC

             DO_ALL_DISP = ELPH_SETTINGS%DO_PERTURB() .OR. ELPH_DERIV_SETTINGS%GENERATE
             RETURN_FC = DO_ALL_DISP .OR. PHON_SETTINGS%DO_INIT()

             DYN%POSIOC=DYN%POSION

             CALL FINITE_DIFF_ID( INFO%LSTOP, DYN%POTIM, T_INFO,  &
               DYN%POSION, TOTEN, TIFOR, TSIF, DYN%NFREE, &
               DYN%ISIF>=3, LATT_CUR%A, LATT_CUR%B, WDES%SAXIS, SYMM, PRIM_CELL, NCDIJ, &
               DYN%TEBEG, DO_ALL_DISP, IO%IU6, IO%IU0 , IO%NWRITE, FORCE_CONST, RETURN_FC, &
               CHECKPOINT_FD_SETTINGS, FINDIFF_CHECKPOINT, &
               GRIDC, WDES, CVTOT, CVTOT_DERIV, CDIJ, CDIJ_DERIV, &
               HAMILTONIAN%MUTOT, MU_DERIV, LDO_AB_INITIO, ELPH_DERIV_SETTINGS%GENERATE, XC_DATA%LTAU, &
               ELPH_INTERPOL, ELPH_SETTINGS, GRID, P, HAMILTONIAN, &
               CQIJ, SV, NONL_S, NONLR_S, E%EXHF, IO)
           END BLOCK
           CALL LATTIC(LATT_CUR)
        ! we need to reinitialise the symmetry code at this point
        ! the number of k-points is changed on the fly
           IF (SYMM%ISYM>0 .AND. .NOT. CHECKPOINT_FD_SETTINGS%DO_SINGLE()) THEN
             CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,DYN%EFOR,T_INFO%LSFOR, &
                T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,NIOND, &
                SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
#ifdef oldsym
             CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
                  SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
                  T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#else
             CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
                  SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
                  T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#endif
             IF(LDO_AB_INITIO) CALL KPAR_SYNC_ALL(WDES,W)
             CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
             IF(LDO_AB_INITIO) CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
             IF (INFO%LONESW.AND.LDO_AB_INITIO) THEN
                CALL DEALLOCW(W_F)
                CALL DEALLOCW(W_G)
                DEALLOCATE(CHAM, CHF)
                CALL ALLOCW(WDES,W_F,WTMP,WTMP)
                CALL ALLOCW(WDES,W_G,WTMP,WTMP)
                ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
                CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
             ENDIF
           ENDIF
        END IF !LJ_ONLY

      ELSE IF (DYN%IBRION==7.OR.DYN%IBRION==8) THEN ibrion

       IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
         CALL LR_SKELETON( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET_GL,DOS,DOSI,CHAM, &
          DYN%IBRION,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, TIFOR, &
          PHON_SETTINGS%DO_INIT(), FORCE_CONST)

       END IF ! LJ_ONLY AND LDO_AB_INITIO
!=======================================================================
! meaning of DYN%ISIF :
!  DYN%ISIF  calculate                           relax
!        force     stress                    ions      lattice
!   0     X                                   X
!   1     X        uniform                    X
!   2     X          X                        X
!   3     X          X                        X          X
!   4     X          X                        X          X **
!   5     X          X                                   X **
!   6     X          X                                   X
!   7     X          X                                  uniform
!   8     X          X                        X         uniform
!
!   **  for DYN%ISIF=4 & DYN%ISIF=5 isotropic pressure will be subtracted
!       (-> cell volume constant, optimize only the cell shape)
!
!=======================================================================
      ELSE IF (DYN%IBRION>0) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
! sum energy of images along chain
        EENTROPY=E%EENTROPY*SCALEE

        IF (IMAGES==0) THEN
           CALLMPI_C( sum_chain( TOTEN ))
           CALLMPI_C( sum_chain( EENTROPY ))
        ENDIF

        io_begin
! If LDO_AB_INITIO=.TRUE., write the information on the electronic
! entropy, too.
        IF(LDO_AB_INITIO) THEN
           WRITE(17,7281,ADVANCE='NO') NSTEP,TOTEN, &
                TOTEN-EENTROPY/(2+NORDER),TOTEN-TOTENG
           IF (IO%IU0>=0) &
                WRITE(TIU0,7281,ADVANCE='NO')NSTEP,TOTEN, &
                TOTEN-EENTROPY/(2+NORDER),TOTEN-TOTENG

           IF ( WDES%NCDIJ>=2 ) THEN
              WRITE(17,77281) RHOTOT(2:WDES%NCDIJ)
              IF (IO%IU0>=0) WRITE(TIU0,77281) RHOTOT(2:WDES%NCDIJ)
           ELSE
              WRITE(17,*)
              IF (IO%IU0>=0) WRITE(TIU0,*)
           ENDIF
        ELSE
           WRITE(17,7281,ADVANCE='NO') NSTEP,TOTEN, &
                TOTEN,TOTEN-TOTENG
           IF (IO%IU0>=0) &
                WRITE(TIU0,7281,ADVANCE='NO')NSTEP,TOTEN, &
                TOTEN,TOTEN-TOTENG

                WRITE(17,*)
                IF (IO%IU0>=0) WRITE(TIU0,*)
        ENDIF

        io_end
 7281 FORMAT(I4,' F= ',E14.8,' E0= ',E14.8,1X,' d E =',E12.6)
77281 FORMAT('  mag=',3F11.4)
!-----------------------------------------------------------------------
!  set DYN%D2C to forces in cartesian coordinates multiplied by FACT
!  FACT is determined from timestep in a way, that a stable timestep
!   gives a good trial step
!-----------------------------------------------------------------------
        FACT=0
        IF ( DYN%ISIF<5 .OR. DYN%ISIF .EQ. 8 ) FACT=10*DYN%POTIM*EVTOJ/AMTOKG *1E-10_q
        LSTOP2=.TRUE.

        NI=1
        DO NT=1,T_INFO%NTYP
        DO NI=NI,T_INFO%NITYP(NT)+NI-1
           DYN%D2C(1,NI)=TIFOR(1,NI)*FACT
           DYN%D2C(2,NI)=TIFOR(2,NI)*FACT
           DYN%D2C(3,NI)=TIFOR(3,NI)*FACT
           IF (SQRT(TIFOR(1,NI)**2+TIFOR(2,NI)**2+TIFOR(3,NI)**2) &
                &       >ABS(DYN%EDIFFG)) LSTOP2=.FALSE.
        ENDDO
        ENDDO
! for all DYN%ISIF greater or equal 3 cell shape optimisations will be done
        FACTSI = 0
        IF (DYN%ISIF>=3) FACTSI=10*DYN%POTIM*EVTOJ/AMTOKG/T_INFO%NIONS *1E-10_q

        DO I=1,3
        DO K=1,3
           D2SIF(I,K)=TSIF(I,K)*FACTSI
        ENDDO
        D2SIF(I,I)=D2SIF(I,I)-DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA*FACTSI
        ENDDO
! For DYN%ISIF =4 or =5 we take only pure shear stresses: subtract pressure
        IF ((DYN%ISIF==4).OR.(DYN%ISIF==5)) THEN
           DO I=1,3
              D2SIF(I,I)=D2SIF(I,I)-PRESS*FACTSI
           ENDDO
        ENDIF
! For DYN%ISIF =7 take only pressure (volume relaxation)
        IF ( DYN%ISIF==7 .OR. DYN%ISIF .EQ. 8 ) THEN
           DO I=1,3
           DO J=1,3
              D2SIF(J,I)=0
           ENDDO
           D2SIF(I,I)=PRESS*FACTSI
           ENDDO
        ENDIF

        CALL CONSTR_CELL_RELAX(D2SIF)
#ifdef tbdyn
        IF (DYN%IBRION==1.OR. DYN%IBRION==2) &
           CALL APPLY_STRESS_CONSTRAINT( D2SIF, TSIF, LATTICE_CONSTRAINTS_GL )
#endif

        IF (FACTSI/=0) THEN
           DO I=1,3
           DO J=1,3
              IF (FACTSI/=0) THEN
                 IF (ABS(D2SIF(J,I))/FACTSI/T_INFO%NIONS>ABS(DYN%EDIFFG)) LSTOP2=.FALSE.
              ENDIF
           ENDDO
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!  do relaxations using diverse algorithms
!-----------------------------------------------------------------------
        ! change of the energy between two ionic step used as stopping criterion
        INFO%LSTOP=(ABS(TOTEN-TOTENG)<DYN%EDIFFG)

        CALLMPI_C ( and_chain( LSTOP2 ))
        CALLMPI_C ( and_chain( INFO%LSTOP ))

        ! IFLAG=0 means no reinit of wavefunction prediction
        IFLAG=0
        IF (DYN%IBRION==1) THEN
           CALL BRIONS(T_INFO%NIONS,DYN%POSION,DYN%POSIOC,DYN%D2C,LATT_CUR%A,LATT_CUR%B,D2SIF, &
                MAX(DYN%NSW+1,DYN%NFREE+1),DYN%NFREE,IO%IU6,IO%IU0,FACT,FACTSI,E1TEST)
! Sometimes there is the danger that the optimisation scheme (especially
! the Broyden scheme) fools itself by performing a 'too small' step - to
! avoid this use a second break condition ('trial step energy change'):
           INFO%LSTOP=INFO%LSTOP .AND. (ABS(E1TEST) < DYN%EDIFFG)
! if we have very small forces (small trial energy change) we can stop
           INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==3 .AND. DYN%SMASS<=0) THEN
           CALL ION_VEL_QUENCH(T_INFO%NIONS,LATT_CUR%A,LATT_CUR%B,IO%IU6,IO%IU0, &
                T_INFO%LSDYN, &
                DYN%POSION,DYN%POSIOC,FACT,DYN%D2C,FACTSI,D2SIF,DYN%D2,E1TEST)
           IF (IFLAG==1) INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==3) THEN
           CALL IONDAMPED(T_INFO%NIONS,LATT_CUR%A,LATT_CUR%B,IO%IU6,IO%IU0, &
                T_INFO%LSDYN, &
                DYN%POSION,DYN%POSIOC,FACT,DYN%D2C,FACTSI,D2SIF,DYN%D2,E1TEST,DYN%SMASS)
           IF (IFLAG==1) INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==2) THEN
           IFLAG=1
           IF (NSTEP==1) IFLAG=0
!  set accuracy of energy (determines whether cubic interpolation is used)
           EACC=MAX(ABS(INFO%EDIFF),ABS(ECONV))

           CALLMPI_C( sum_chain( EACC))

           IF ( LHYPER_NUDGE() ) EACC=1E10    ! energy not very accurate, use only force
           CALL IONCGR(IFLAG,T_INFO%NIONS,TOTEN,LATT_CUR%A,LATT_CUR%B,DYN%NFREE,DYN%POSION,DYN%POSIOC, &
                FACT,DYN%D2C,FACTSI,D2SIF,DYN%D2,DYN%D3,DISMAX,IO%IU6,IO%IU0, &
                EACC,DYN%EDIFFG,E1TEST,LSTOP2)

!    if IFLAG=1 new trial step -> reinit of waveprediction
           INFO%LSTOP=.FALSE.
           IF (IFLAG==1) THEN
              INFO%LSTOP=(ABS(TOTEN-TOTENG)<DYN%EDIFFG)
              TOTENG=TOTEN
           ENDIF
           IF (IFLAG==2) INFO%LSTOP=.TRUE.
! if we have very small forces (small trial energy change) we can stop
           IF (IFLAG==1) INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
!-----------------------------------------------------------------------
! dimer method
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==44) THEN
           CALL dimer(DYN,T_INFO,INFO,LATT_CUR,IO,TOTEN,FACT,LSTOP2)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
! damped velocity Verler to calculate IRC (not optimization!!!)
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==40) THEN
           FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
           NI=1
           DO NT=1,T_INFO%NTYP
             DO NI=NI,T_INFO%NITYP(NT)+NI-1
               DYN%D2C(1,NI)=TIFOR(1,NI)*FACT/2/T_INFO%POMASS(NT)
               DYN%D2C(2,NI)=TIFOR(2,NI)*FACT/2/T_INFO%POMASS(NT)
               DYN%D2C(3,NI)=TIFOR(3,NI)*FACT/2/T_INFO%POMASS(NT)
             ENDDO
           ENDDO
           CALL KARDIR(T_INFO%NIONS,DYN%D2C,LATT_CUR%B)

           CALL dvv(DYN,T_INFO,INFO,LATT_CUR,IO,TOTEN,FACT)
           FACT=10*DYN%POTIM*EVTOJ/AMTOKG *1E-10_q
           NI=1
           DO NT=1,T_INFO%NTYP
             DO NI=NI,T_INFO%NITYP(NT)+NI-1
               DYN%D2C(1,NI)=TIFOR(1,NI)*FACT
               DYN%D2C(2,NI)=TIFOR(2,NI)*FACT
               DYN%D2C(3,NI)=TIFOR(3,NI)*FACT
             ENDDO
           ENDDO
           TOTENG=TOTEN
!-----------------------------------------------------------------------
! interactive mode
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==11) THEN
io_begin
           IF (IO%LOPEN) CALL WFORCE(IO%IU6)            ! write of OUTCAR
           IF (IO%LOPEN) CALL WFORCE(17)                ! write of OSZICAR
io_end
           DYN%POSIOC(:,:) = DYN%POSION(:,:)
           CALL INPOS(LATT_CUR, T_INFO, DYN, IO%IU6, IO%IU0, INFO%LSTOP, WDES%COMM)
           IF ( INFO%NELMDL==0) INFO%NELMDL=-5 ! gK if the positions are widely different we need a delay before mixing
        ELSE IF (DYN%IBRION==12) THEN
           DYN%POSIOC(:,:) = DYN%POSION(:,:)
           CALL PLUGINS_STRUCTURE(GRIDC, INFO, T_INFO, LATT_CUR, P, COMM, TOTEN, TIFOR, TSIF, CHTOT)
           TOTENG=TOTEN
        ENDIF

        CALLMPI_C ( and_chain( LSTOP2 ))
        CALLMPI_C ( and_chain( INFO%LSTOP ))

! restrict volume for constant volume relaxation
        IF (DYN%ISIF==4 .OR. DYN%ISIF==5) THEN
           OMEGA_OLD=LATT_CUR%OMEGA
           CALL LATTIC(LATT_CUR)
           SCALEQ=(ABS(OMEGA_OLD) / ABS(LATT_CUR%OMEGA))**(1._q/REAL(COUNT(LATTICE_CONSTRAINTS_GL),KIND=q))
           DO I=1,3
              IF ( LATTICE_CONSTRAINTS_GL( I ) ) THEN
                 LATT_CUR%A(1,I)=SCALEQ*LATT_CUR%A(1,I)
                 LATT_CUR%A(2,I)=SCALEQ*LATT_CUR%A(2,I)
                 LATT_CUR%A(3,I)=SCALEQ*LATT_CUR%A(3,I)
              END IF
           ENDDO
        ENDIF
        CALL LATTIC(LATT_CUR)
!-----------------------------------------------------------------------
!  reinitialize the prediction algorithm for the wavefunction if needed
!-----------------------------------------------------------------------
        PRED%INIPRE=3
        IF ( PRED%IWAVPR >=12 .AND. &
             &     (ABS(TOTEN-TOTENG)/T_INFO%NIONS>1.0_q .OR. IFLAG==1)) THEN
           CALL WAVPRE_NOIO(GRIDC,P,XC_DATA,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
                CHTOT,RHOLM,N_MIX_PAW, CSTRF, KINEDEN, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)

        ELSE IF ( PRED%IWAVPR >=2 .AND. PRED%IWAVPR <10   .AND. &
             &     (ABS(TOTEN-TOTENG)/T_INFO%NIONS>1.0_q .OR. IFLAG==1)) THEN
           CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
                CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
        ENDIF

        ! use forces as stopping criterion if EDIFFG<0
        IF (DYN%EDIFFG<0) INFO%LSTOP=LSTOP2
        io_begin
        WRITE(TIU6,130)

        IF (INFO%LSTOP) THEN
           IF (IO%IU0>=0) &
                WRITE(TIU0,*) 'reached required accuracy - stopping ', &
                'structural energy minimisation'
           WRITE(TIU6,*) ' '
           WRITE(TIU6,*) 'reached required accuracy - stopping ', &
                'structural energy minimisation'
        ENDIF
        io_end
       END IF !LJ_ONLY
      ELSE IF ( DYN%IBRION == -1 .AND. (LEPSILON  .OR. WAVE_INTERP_SET%LKPOINTS_OPT .OR. LMAGBLOCH)) THEN ibrion
       IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
         IF (WAVE_INTERP_SET%LKPOINTS_OPT) THEN
           CALL INTERPOLATE_OPTICS_DENSE( WAVE_INTERP_SET, HAMILTONIAN, E, KPOINTS, GRID, LATT_CUR, LATT_INI, &
                T_INFO, NONLR_S, NONL_S, W, NEDOS, &
                LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, GRIDC, GRIDUS, C_TO_US, IRDMAX)
         ELSE
           CALL LR_SKELETON( &
                KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
                GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
                CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
                CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
                CHDEN,SV,VDW_SET_GL,DOS,DOSI,CHAM, &
                DYN%IBRION,LMDIM,IRDMAX,NEDOS, &
                TOTEN,EFERMI, TIFOR, &
                .FALSE.)
         ENDIF
       END IF !LJ_ONLY AND LDO_AB_INITIO
      ELSE IF ( DYN%IBRION == -1 .AND. LAUGER) THEN ibrion
       IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
         CALL CALCULATE_AUGER(W,WDES,LATT_CUR,SYMM,T_INFO,GRID,P,NONL_S,KPOINTS,IO)
       END IF !LJ_ONLY AND LDO_AB_INITIO
      ELSE IF ( LCHIMAG ) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
         CALL MLR_SKELETON( &
              HAMILTONIAN,W,WDES,GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,KPOINTS,LATT_CUR,LATT_INI, &
              T_INFO,DYN,SYMM,P,XC_DATA,NONL_S,NONLR_S,LMDIM,CDIJ,CQIJ,SV,E,INFO,IRDMAX,IO)
       END IF !LJ_ONLY
      ENDIF ibrion
!=======================================================================
!  update of ionic positions performed
!  in any case POSION should now hold the new positions and
!  POSIOC the old one
!=======================================================================

!   reached required number of time steps
      IF (NSTEP>=NSW_ACT) INFO%LSTOP=.TRUE.
#ifdef tbdyn
      IF (INFO%LSTOP) THEN
        CLOSE(g_io%REPORT)
      END IF
#endif

!   if we need to pull the brake, then POSION is reset to POSIOC
!   except for molecular dynamics, where the next electronic
!   step should indeed correspond to POSIOC
      IF ( INFO%LSTOP .AND. DYN%IBRION>0 ) THEN
         LATT_CUR%A=DYN%AC
         DYN%POSION=DYN%POSIOC
      ENDIF

!      CALLMPI ( M_bcast_d(WDES%COMM, DYN%POSION , T_INFO%NIONS*3))
!      CALLMPI ( M_bcast_d(WDES%COMM, LATT_CUR%A(1,1) , 9))
      CALL LATTIC(LATT_CUR)

!------------------------------------------------------------------------
! For VCAIMAGES the ionic positions and Bravais matrix have to be in
! sync on every image/process throughout the whole MD.
! So here we broadcast them and hence sync them.
!------------------------------------------------------------------------
      IF (LVCAIMAGES) THEN
         CALLMPI ( M_bcast_d(COMM_WORLD, DYN%POSION , T_INFO%NIONS*3))
         CALLMPI ( M_bcast_d(COMM_WORLD, LATT_CUR%A(1,1) , 9))
         CALLMPI ( M_bcast_d(COMM_WORLD, LATT_CUR%ANORM(1) , 3))
         CALLMPI ( M_bcast_d(COMM_WORLD, LATT_CUR%B(1,1) , 9))
         CALLMPI ( M_bcast_d(COMM_WORLD, LATT_CUR%BNORM(1) , 3))
         CALLMPI ( M_bcast_d(COMM_WORLD, LATT_CUR%OMEGA , 1))
      ENDIF

!=======================================================================
!  update mean temperature mean energy
!=======================================================================
      SMEAN =SMEAN +1._q/DYN%SNOSE(1)
      SMEAN0=SMEAN0+DYN%SNOSE(1)
      TMEAN =TMEAN +TEIN/DYN%SNOSE(1)
      TMEAN0=TMEAN0+TEIN

    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
!=======================================================================
!  SMEAR_LOOP%ISMCNT != 0 Loop over several KPOINTS%SIGMA-values
!  set new smearing parameters and continue main loop
!=======================================================================
      IF (SMEAR_LOOP%ISMCNT/=0) THEN
         KPOINTS%ISMEAR=NINT(SMEAR_LOOP%SMEARS(2*SMEAR_LOOP%ISMCNT-1))
         KPOINTS%SIGMA=SMEAR_LOOP%SMEARS(2*SMEAR_LOOP%ISMCNT)

         SMEAR_LOOP%ISMCNT=SMEAR_LOOP%ISMCNT-1

         io_begin
         IF (IO%IU0>=0) &
              WRITE(TIU0,7283) KPOINTS%ISMEAR,KPOINTS%SIGMA
         WRITE(17,7283) KPOINTS%ISMEAR,KPOINTS%SIGMA
         io_end

7283     FORMAT('ISMEAR = ',I4,' SIGMA = ',F10.3)

         KPOINTS%LTET=((KPOINTS%ISMEAR<=-4).OR.(KPOINTS%ISMEAR>=30))
         IF (KPOINTS%ISMEAR==-6) KPOINTS%ISMEAR=-1
         IF (KPOINTS%ISMEAR>=0)  KPOINTS%ISMEAR=MOD(KPOINTS%ISMEAR,30)

         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      ENDIF

      IF (INFO%PLUGIN%OCCUPANCIES) THEN
         CALL PLUGINS_OCCUPANCIES(INFO, COMM, INFO%NELECT, EFERMI, KPOINTS%EFERMI, KPOINTS%SIGMA, &
              KPOINTS%ISMEAR, KPOINTS%EMIN, KPOINTS%EMAX, INFO%NUP_DOWN)
4488     FORMAT('PLUGINS/OCCUPANCIES has triggered the following new values'/, &
                & 'NELECT  = ', F10.3, / &
                & 'EFERMI  = ', F10.3, / &
                & 'SIGMA   = ', F10.3, / &
                & 'ISMEAR  = ', I4,    / &
                & 'EMIN    = ', F10.3, / &
                & 'EMAX    = ', F10.3, / &
                & 'NUPDOWN = ', F10.3, / &
                & "Recomputing occupancies with these values."/)
         io_begin
         IF (IO%IU0>=0) THEN
             IF (KPOINTS%EFERMI == HUGE(1.0_q)) THEN
                 WRITE(TIU6, 4488) INFO%NELECT, EFERMI, KPOINTS%SIGMA, KPOINTS%ISMEAR, &
                                   KPOINTS%EMIN, KPOINTS%EMAX, INFO%NUP_DOWN
             ELSE
                 WRITE(TIU6, 4488) INFO%NELECT, KPOINTS%EFERMI, KPOINTS%SIGMA, KPOINTS%ISMEAR, &
                                   KPOINTS%EMIN, KPOINTS%EMAX, INFO%NUP_DOWN
             END IF
         END IF
         io_end
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      END IF
!=======================================================================
! tasks which are done all DYN%NBLOCK   steps:
!=======================================================================
!-----------------------------------------------------------------------
!  write out the position of the IONS to XDATCAR
!-----------------------------------------------------------------------
      IF (INFO%INICHG==3) THEN
        IF (MOD(NSTEP,DYN%NBLOCK)==1) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*) 'non selfconsistent'
           INFO%LCHCON=.FALSE.
        ELSE
           IF (IO%IU0>=0) &
           WRITE(TIU0,*) 'selfconsistent'
           INFO%LCHCON=.TRUE.
        ENDIF
      ENDIF
    END IF !LJ_ONLY AND LDO_AB_INITIO

   nblock: IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
        io_begin
        ! tb: write lattice vectors to XDATCAR
        IF (DYN%ISIF>=3) THEN
           CALL XDAT_HEAD(61, T_INFO, DYN%AC, INFO%SZNAM1)
        ENDIF
        WRITE(61,'(A,I12)') 'Direct configuration=',NSTEP
        WRITE(61,7007) ((DYN%POSIOC(I,J),I=1,3),J=1,T_INFO%NIONS)
#ifdef VASP_HDF5
        CALL VH5_WRITE_LATTICE_ION_HIST_STEP(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,NSTEP/DYN%NBLOCK,DYN%AC,DYN%POSIOC)
        IF (IO%VELOCITY) &
           CALL VH5_WRITE_VELOCITY_HIST_STEP(IH5INTERMEDIATEGROUP_ID,GRP_HISTORY,NSTEP/DYN%NBLOCK,LATT_CUR,DYN)
#endif
        IF (IO%LOPEN) CALL WFORCE(61)
        io_end

    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
!-----------------------------------------------------------------------
! acummulate dos
!-----------------------------------------------------------------------
        DO ISP=1,WDES%NCDIJ
        DO I=1,NEDOS
          DDOSI(I,ISP)=DDOSI(I,ISP)+DOSI(I,ISP)
          DDOS (I,ISP)=DDOS (I,ISP)+DOS (I,ISP)
        ENDDO
        ENDDO
    END IF !LJ_ONLY AND LDO_AB_INITIO

!-----------------------------------------------------------------------
! evaluate the pair-correlation function using the positions
! considered in DFT
!-----------------------------------------------------------------------
       IF (OUTPUT_MODE.GT.0) THEN
          CALL PCDAT_PACO(1._q, T_INFO, PACO, DYN%AC, DYN%POSIOC, COMM)
       ENDIF
    ENDIF nblock

!=======================================================================
! tasks which are done all DYN%NBLOCK*DYN%KBLOCK   steps
!=======================================================================
!-----------------------------------------------------------------------
! write  pair-correlation and density of states
! quantities are initialized to 0 at the beginning of the main-loop
!-----------------------------------------------------------------------
    wrtpair: IF (MOD(NSTEP,DYN%NBLOCK*DYN%KBLOCK)==0 .AND. DYN%IBRION ==0) THEN
       IF (OUTPUT_MODE.GT.0) THEN
          CALL PCDAT_FINALIZE(T_INFO, PACO, COMM)
       ENDIF

 io_begin
       IF (OUTPUT_MODE.GT.0) THEN
          CALL PCDAT_WRITE(60, PACO, TMEAN0/(DYN%NBLOCK*DYN%KBLOCK), TMEAN/SMEAN)
#ifdef VASP_HDF5
          CALL VH5_WRITE_PAIR_CORRELATION_STEP(IH5INTERMEDIATEGROUP_ID, GRP_PAIR_CORRELATION, NSTEP, DYN, PACO)
#endif
       ENDIF
       WRITE(TIU6,8022) SMEAN0/(DYN%NBLOCK*DYN%KBLOCK),TMEAN0/(DYN%NBLOCK*DYN%KBLOCK), &
            &                TMEAN/SMEAN
8022   FORMAT(/' mean value of Nose-termostat <S>:',F10.3, &
            &        ' mean value of <T> :',F10.3/ &
            &        ' mean temperature <T/S>/<1/S>  :',F10.3/)


    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
       DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(NEDOS-1)

       WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
       DO I=1,NEDOS
          EN=KPOINTS%EMIN+DELTAE*(I-1)
          WRITE(16,7062) EN,(DDOS(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN),(DDOSI(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN)
       ENDDO
       IF (IO%LOPEN) CALL WFORCE(16)
7062   FORMAT(3X,F8.3,8E12.4)

#ifdef VASP_HDF5
       CALL VH5_WRITE_DOS(IH5OUTFILEID, WDES, KPOINTS, DDOS, DDOSI, DOSPAR, EFERMI, T_INFO%NIONP, -1, SCALE=1.d0/DYN%KBLOCK)
#endif
       CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .FALSE., &
            DDOS, DDOSI, DOSPAR, NEDOS, LPAR, T_INFO%NIONP, WDES%NCDIJ)
    END IF !LJ_ONLY AND LDO_AB_INITIO
 io_end
    ENDIF wrtpair
!-----------------------------------------------------------------------
!  update file CONTCAR
!-----------------------------------------------------------------------

!-----write out positions (only done on IONODE)
#ifdef MPI
    IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
      io_begin

      IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
         CALL OUTPOS(13,.TRUE.,T_INFO%SZNAM2,T_INFO,LATT_CUR%SCALE,DYN%AC,T_INFO%LSDYN,DYN%POSIOC,LEGACY=.FALSE.)
         CALL OUTPOS_TRAIL(13,IO%LOPEN, LATT_CUR, T_INFO, DYN, NHC)
      ENDIF

#ifdef VASP_HDF5
      IF (MOD(NSTEP,LBLOCK_HELP)==0) THEN
         CALL VH5_WRITE_LATTICE_ION(IH5OUTFILEID, INFO%SZNAM1, T_INFO, LATT_CUR%SCALE, DYN%AC, T_INFO%LSDYN, DYN%POSIOC)
         CALL VH5_WRITE_POSITIONS_TRAILER(IH5OUTFILEID, LATT_CUR, T_INFO, DYN, NHC)
      ENDIF
#endif

      io_end

    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
!=======================================================================
!  append new chargedensity to file CHG
!=======================================================================
      IF (.NOT.IO%LH5 .AND. IO%WRT_DENSITY%LCHARG .AND. MOD(NSTEP,10)==1) THEN

      io_begin
      CALL OUTPOS(70,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSIOC)
      io_end

      DO ISP=1,WDES%NCDIJ
         CALL OUTCHG(GRIDC,70,.FALSE.,CHTOT(:,ISP))
      ENDDO
      io_begin
      IF (IO%LOPEN) CALL WFORCE(70)
      io_end
      ENDIF
    END IF !LJ_ONLY AND LDO_AB_INITIO

#ifdef MPI
    ENDIF
#endif

    ! For ML and related (ML-Plugins) we stop the LOOP+ timer here because the next "Ionic step" in
    ! OUTCAR starts just below (see write(...,140)). In constrast to pure AIMD, timings of all
    ! following preparation work ("WAVPRE",...) are accounted to the next ionic step.
    IF (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD) THEN
       CALL STOP_TIMING("LOOP+",IO%IU6,XMLTAG='totalsc',WRITE_FILES=MOD(NSTEP,LBLOCK_HELP)==0)
       IF(MOD(NSTEP,LBLOCK_HELP)==0 .AND. LBLOCK_HELP>=2 ) THEN
          CALL STOP_TIMING("BLOCK+",IO%IU6,NAME="BLOCK+",XMLTAG='totalblock')
       ENDIF
    ENDIF

#ifdef ML_AVAILABLE
! Machine-learning prediction on the new configuration is executed here.
! Must call ML_TO_VASP here to decide if an ab-initio calculation 
! needs to be done in the next step. If so then the wave functions will
! be predicted below. 
    IF (.NOT.INFO%LSTOP) THEN
       IF((NSTEP+1<=NSW_ACT)) THEN
          IF (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD) THEN
             ! IO for next ionic step here if ML is .TRUE., otherwise IO at beginning of ion loop
             IF (IO%IU6>=0 .AND. MOD(NSTEP+1,LBLOCK_HELP)==0) WRITE(TIU6,140) NSTEP+1
          ENDIF
          IF(ML_LMLFF) THEN
             DO I = 1, ML_TOTNUM_INSTANCES
                IF (ML_SUPER_HANDLE_MAIN(1)%FF%LIB .AND. (ML_SUPER_HANDLE_MAIN(1)%FF%ISTART == 2)) THEN
#ifdef libvaspml
                   CALL VASPML%compute(T_INFO%POSION, &
                                       LATT_CUR%A, &
                                       LATT_CUR%B, &
                                       T_INFO%NIONS, &
                                       T_INFO%NTYP, &
                                       T_INFO%ITYP, &
                                       T_INFO%TYPE, &
                                       LATT_CUR%OMEGA, &
                                       TOTEN, &
                                       TIFOR, &
                                       TSIF)
#endif
                ELSE
                   CALL ML_TO_VASP_MACHINE_LEARNING(ML_SUPER_HANDLE_MAIN(I), DYN, LATT_CUR, &
                        NIOND, NSTEP+1, T_INFO, TOTEN, TIFOR, TSIF, COMM%MPI_COMM)
                   CALL TRIGGER_SOFT_STOP(INFO%LSOFT_TRIGGER, ML_REQUESTS_EXIT(ML_SUPER_HANDLE_MAIN(I)))
                ENDIF
                IF (SYMM%ISYM>0) THEN
                   CALL FORSYM(TIFOR,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
                   IF (DYN%ISIF/=0) CALL TSYM(TSIF,ISYMOP,NROTK,LATT_CUR%A)
                ENDIF
                PRESS=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))/3._q -DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
             ENDDO
             IF (.NOT. LDO_AB_INITIO_NEW .AND. INFO%NELMDL==0) INFO%NELMDL=-5
          ELSE IF (INFO%PLUGIN%FORCE_FIELD) THEN
               CALL PLUGINS_MACHINE_LEARNING(INFO, T_INFO, LATT_CUR, P, GRIDC, DYN, &
                     SYMM, ISYMOP, NROTK, LDO_AB_INITIO, TOTEN, TIFOR, TSIF, PRESS)
          ENDIF
          CALL STOP_TIMING("G",TIU6,'MLFF',WRITE_FILES=MOD(NSTEP+1,LBLOCK_HELP)==0)
       ENDIF
    ENDIF
#endif

    IF (.NOT. LJ_ONLY.AND.LDO_AB_INITIO_NEW) THEN
!=======================================================================
! if ions were moved recalculate some quantities
!=======================================================================
!=======================================================================
! WAVPRE prediction of the new wavefunctions and charge-density
! if charge-density constant during ELM recalculate the charge-density
! according to overlapping atoms
! for relaxation jobs do not predict in the last step
!=======================================================================
      INFO%LPOTOK=.FALSE.
  prepare_next_step: &
    & IF ( .NOT. INFO%LSTOP .OR. DYN%IBRION==0 .OR. DYN%IBRION==4) THEN

      ! extrapolate charge by subtracting  atomic charges at old positions
      ! and adding the one at the new positions
      IF (PRED%IWAVPR==1 .OR.PRED%IWAVPR==11) PRED%INIPRE=5
      ! quadratic extrapolation of wavefunctions and charge
      IF (PRED%IWAVPR==2 .OR.PRED%IWAVPR==12) PRED%INIPRE=0
      ! mixed mode
      IF (PRED%IWAVPR==3 .OR.PRED%IWAVPR==13) PRED%INIPRE=4
      PRED%IPRE=0

      IF (PRED%IWAVPR >=11) THEN
        CALL WAVPRE_NOIO(GRIDC,P,XC_DATA,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, KINEDEN, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
      ELSE IF (PRED%IWAVPR >=1 .AND. PRED%IWAVPR<10 ) THEN
        CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
      ENDIF

      IF (PRED%IPRE<0) THEN
        IF (IO%IU0>=0) WRITE(TIU0,*)'bond charge predicted'
        PRED%IPRE=ABS(PRED%IPRE)
      ELSE
      ! PRED%IPRE < 0 then WAVPRE calculated new structure factor
      !     in all other cases we have to recalculate s.f.
         IF (INFO%TURBO==0) CALL STUFAK(GRIDC,T_INFO,CSTRF)
      ENDIF
      IF (PRED%IPRE>1) THEN
        io_begin
        IF (IO%IU0>=0) WRITE(TIU0,*)'prediction of wavefunctions'
        IF (TIU6>=0) WRITE(TIU6,2450) PRED%ALPHA,PRED%BETA
        io_end
 2450   FORMAT(' Prediction of Wavefunctions ALPHA=',F6.3,' BETA=',F6.3)
      ENDIF

      IF (INFO%LCHCON.AND.INFO%INICHG==2) THEN
        IF (IO%IU0>=0)  WRITE(TIU0,*)'charge from overlapping atoms'
        ! then initialize CRHODE and than RHOLM (PAW related occupancies)
__GPU   CALL OFFLOAD_ON_PUSH(.FALSE.)
        CALL DEPATO(WDES, LMDIM, CRHODE, INFO%LOVERL, P, T_INFO)
        CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)

        CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
        ! set magnetization to 0
        DO ISP=2,WDES%NCDIJ
           CALL RC_ADD(CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),GRIDC)
        ENDDO
        ! add Gaussian "charge-transfer" charges, if required
        CALL RHOADD_GAUSSIANS(T_INFO,LATT_CUR,P,GRIDC,NCDIJ,CHTOT,CSTRF)
        CALL RHOADD_GAUSSIANS_LIST(LATT_CUR,GRIDC,NCDIJ,CHTOT)
__GPU   CALL OFFLOAD_ON_POP
      ENDIF

      IF (INFO%LCORE) CALL RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENCOR,IO%IU6)
      IF (XC_DATA%LTAU) CALL TAUPAR(GRIDC,T_INFO,LATT_CUR%B,LATT_CUR%OMEGA,P,CSTRF,KINEDEN%TAUC)

      CALL STOP_TIMING("G",TIU6,'WAVPRE')
!-----------------------------------------------------------------------
! call the ewald program to get the energy of the new ionic
! configuration
!-----------------------------------------------------------------------
      IF(INFO%TURBO==0) THEN
         CALL FEWALD(GRIDC,LATT_CUR,P,T_INFO,HAMILTONIAN%COULOMB_POT,E%TEWEN,EWIFOR,EWSIF)
      ELSE
         ALLOCATE(CWORK1(GRIDC%MPLWV))
         CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK1,E%PSCENC,E%TEWEN,EWIFOR,HAMILTONIAN%COULOMB_POT)
         DEALLOCATE(CWORK1)
      ENDIF

      CALL STOP_TIMING("G",TIU6,'FEWALD')

! volume might have changed restet IRDMAX
      IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &     (GRIDC%NGPTAR(1)*GRIDC%NGPTAR(2)*GRIDC%NGPTAR(3)))+200
#ifdef usgrid
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDUS%NGPTAR(1)*GRIDUS%NGPTAR(2)*GRIDUS%NGPTAR(3)))+200
#endif
!-----------------------------------------------------------------------
! if basis cell changed recalculate kinetic-energy array and tables
!-----------------------------------------------------------------------
      IF (DYN%ISIF>=3) THEN
        CALL GEN_INDEX(GRID,WDES, LATT_CUR%B,LATT_INI%B,-1,-1,.TRUE.)
        CALL STOP_TIMING("G",TIU6,'GENKIN')
        IF (WDES%LSPIRAL) CALL CLEANWAV(WDES,W,INFO%ENINI)
      ENDIF
!-----------------------------------------------------------------------
!  recalculate the real-space projection operators
!  if volume changed also recalculate reciprocal projection operators
!  and reset the cache for the phase-factor
!-----------------------------------------------------------------------
      IF (INFO%LREAL) THEN
! reset IRMAX, IRALLOC if required (no redistribution of GRIDS allowed)

        CALL REAL_OPTLAY(GRID,LATT_CUR,NONLR_S,.TRUE.,LREALLOCATE, IO%IU6,IO%IU0)
        IF (LREALLOCATE) THEN
           ! reallocate real space projectors
           CALL NONLR_DEALLOC(NONLR_S)
           CALL NONLR_ALLOC(NONLR_S)
        END IF
        CALL RSPHER(GRID,NONLR_S,LATT_CUR)

      ELSE
        IF (DYN%ISIF>=3) THEN
          CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR, 1)
        ENDIF
        CALL PHASE(WDES,NONL_S,0)
      ENDIF

      CALL RESETUP_FOCK( WDES, LATT_CUR)
!-----------------------------------------------------------------------
! recalculate projections and perform Gramm-Schmidt orthogonalization
!-----------------------------------------------------------------------
      CALL WVREAL(WDES,GRID,W) ! only for gamma some action
      CALL START_TIMING("G")
      CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
      CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
      CALL REDIS_PW_OVER_BANDS(WDES, W)

      CALL KPAR_SYNC_ALL(WDES,W)

      CALL STOP_TIMING("G",TIU6,'ORTHCH')
!-----------------------------------------------------------------------
! set  INFO%LPOTOK to .F. this requires a recalculation of the local pot.
!-----------------------------------------------------------------------
      INFO%LPOTOK=.FALSE.
!-----------------------------------------------------------------------
! if prediction of wavefunctions was performed and
! diagonalization of sub-space-matrix is selected then
! )  POTLOK: calculate potential according to predicted charge-density
! )  SETDIJ: recalculate the energy of the augmentation charges
! )  then perform a sub-space-diagonal. and generate new Fermi-weights
! )  set INFO%LPOTOK to true because potential is OK
! )  recalculate total energy
! if charge density not constant during band minimization
! )  calculate charge-density according to new wavefunctions
!     and set LPOTOK to false (requires recalculation of loc. potential)
! in all other cases the predicted charge density is used in the next
!   step of ELM
!-----------------------------------------------------------------------
!   wavefunctions are not diagonal, so if they are written to the file
!   or if we do not diagonalize before the optimization
!   rotate them now
  pre_subrot: IF (.NOT.INFO%LDIAG .OR. INFO%LONESW .OR. &
     &     INFO%LSTOP .OR. (MOD(NSTEP,10)==0 &
     &     .AND. PRED%IWAVPR >= 1 .AND.  PRED%IWAVPR < 10) ) THEN

      CALL START_TIMING("G")
      IF (IO%IU0>=0) WRITE(TIU0,*)'wavefunctions rotated'
      CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES, &
           INFO,XC_DATA,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF, &
           HAMILTONIAN%COULOMB_POT)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL SET_DD_PAW(WDES, P, XC_DATA, T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ,  RHOLM, CRHODE, &
          E, LASPH =INFO%LASPH , LCOREL=.FALSE.)

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,XC_DATA,IO%IU6)

      CALL STOP_TIMING("G",TIU6,'POTLOK')
      INFO%LPOTOK=.TRUE.

      IF (INFO%LDIAG) THEN
        IFLAG=3
      ELSE
        IFLAG=4
      ENDIF
      CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF)

      CALL REDIS_PW_OVER_BANDS(WDES, W)

      CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, &
               INFO%NELECT, INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

     ! for the selfconsistent update set W_F%CELTOT and TOTEN now
      IF (INFO%LONESW) W_F%CELTOT=W%CELTOT
      E%EBANDSTR= BANDSTRUCTURE_ENERGY(WDES, W)
      ! These are the correct energy terms but no necessity to write to TOTEN here.
      ! TOTEN = E%EBANDSTR + E%DENC + E%XCENC + E%TEWEN + E%PSCENC + E%EENTROPY &
      !       + E%PAWPS + E%PAWAE + INFO%EALLAT + E%EXHF + E%EDOTP + Ediel_sol

      CALL STOP_TIMING("G",IO%IU6,'EDDIAG')

    ENDIF pre_subrot

    IF (INFO%LONESW .AND. INFO%NELMDL == 0 ) THEN
      CALL START_TIMING("G")

      CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

      CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
           XC_DATA%LTAU,T_INFO%NIONS,W,WDES,KINEDEN)

      INFO%LPOTOK=.FALSE.

      CALL STOP_TIMING("G",IO%IU6,'CHARGE')

    ENDIF

  ELSE prepare_next_step
! in any case we have to call RSPHER at this point even if the ions do not move
! since the force routine overwrites the required arrays
      IF (INFO%LREAL) THEN
         CALL RSPHER(GRID,NONLR_S,LATT_CUR)
      ENDIF
  ENDIF prepare_next_step
!=======================================================================
!  update file WAVECAR if INFO%LSTOP = .TRUE.
!  or if wavefunctions on file TMPCAR are rotated
!=======================================================================
! after 10 steps rotate the wavefunctions on the file
      LTMP= MOD(NSTEP,10) == 0 .AND. PRED%IWAVPR >= 2 .AND. PRED%IWAVPR < 10
  wrtwave: IF (.NOT.IO%LH5 .AND. IO%LWAVE .AND. ( INFO%LSTOP .OR. LTMP .OR. LHFCALC)) THEN

      CALL OUTWAV(IO, WDES, W, LATT_INI,  EFERMI)
!-----------------------------------------------------------------------
! rotate wavefunctions on file (gives a better prediction)
!------------------------------------------------------------------------
      IF (LTMP) THEN
        PRED%INIPRE=10

        CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)

        IF (IO%IU0>=0) &
             WRITE(TIU0,*)'wavefunctions on file TMPCAR rotated'
! and read in wavefunctions  (destroyed by WAVPRE)
        CALL OPENWAV(IO, COMM)
        CALL INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, INFO%ISTART,  EFERMI )
        CALL CLOSEWAV

        CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
        CALL REDIS_PW_OVER_BANDS(WDES, W)
      ENDIF
   ENDIF wrtwave


#ifdef VASP_HDF5
   IF (IO%LH5 .OR. IO%WRT_DENSITY%LCHARGH5) THEN
      IF (IO%LH5 .AND. IO%LWAVE .AND. ( INFO%LSTOP .OR. LTMP .OR. LHFCALC)) THEN
         CALL VH5_WRITE_WAVEFUNCTIONS(IH5WAVEFILEID, IO, WDES, W, LATT_INI, EFERMI)
      END IF
   ENDIF
#endif

    END IF !LJ_ONLY AND LDO_AB_INITIO

!=======================================================================
! next electronic energy minimisation
      
    ! ML and related already stopped the LOOP+ timer above.
    IF (.NOT. (ML_LMLFF .OR. INFO%PLUGIN%FORCE_FIELD)) THEN
      CALL STOP_TIMING("LOOP+",IO%IU6,XMLTAG='totalsc')
    ENDIF

      ! for the final ionic step the calculation tag is closed later,
      ! as postprocessing might still take place
      IF (.NOT. INFO%LSTOP .AND. LDO_AB_INITIO .AND. (.NOT. LJ_ONLY))  THEN
         CALL XML_CLOSE_TAG("calculation")
         CALL XML_FLUSH
      ENDIF

! If INFO%LSTOP=.FALSE., overwrite LDO_AB_INITIO for the next MD step.
#if defined(MPI) || defined(MPI_CHAIN)
      IF (.NOT.INFO%LSTOP.AND.ML_LMLFF) THEN
         LDO_AB_INITIO=LDO_AB_INITIO_NEW
      ENDIF
#endif

      IF (LDO_AB_INITIO) THEN
         IF (isMetal(W%FERTOT).AND.(KPOINTS%ISMEAR == -4.OR.KPOINTS%ISMEAR == -5)) THEN
            CALL vtutor%alert("Tetrahedron method does not include variations&
               & of the Fermi occupations, so forces and stress will be&
               & inaccurate. We suggest using a different smearing scheme&
               & ISMEAR = 1 or 2 for metals in relaxations.")
         ENDIF
      ENDIF

      ENDDO ion
!!      PROFILING_STOP('ion_loop_main')
!=======================================================================
! here we are at the end of the required number of timesteps
!=======================================================================
      io_begin
      IF (IO%LOPEN) CALL WFORCE(IO%IU6)
      io_end

!=======================================================================
!  calculate electron-phonon interactions
!=======================================================================

      IF (ELPH_SETTINGS%DO_WSWQ) CALL ELPH_OVERLAP_WSWQ( &
         W, LATT_CUR, CQIJ, NONLR_S, NONL_S, P, T_INFO, &
         ELPH_SETTINGS%ALL_ELECTRON, IO)

      IF (CHECKPOINT_FD_SETTINGS%DOES_FINISH()) THEN
         IF (.NOT. ELPH_SETTINGS%DO_LOAD_DATA() .AND. ( PHON_SETTINGS%DO_INIT() .OR. ELPH_SETTINGS%DO_PERTURB() )) &
            CALL PHON_INIT(PRIM_CELL, FORCE_CONST, IO)

         CALL ELPH_IO(ELPH_INTERPOL, ELPH_SETTINGS, PRIM_CELL, T_INFO, POLAR_DATA, IO)

         IF (PHON_SETTINGS%DO_POLAR) CALL PHON_INIT_POLAR(PHON_SETTINGS, PRIM_CELL, POLAR_DATA)
         IF (ELPH_SETTINGS%DO_POLAR) CALL ELPH_INTERPOL%INIT_POLAR(ELPH_SETTINGS, POLAR_DATA)

         CALL PHONON_DRIVER(PRIM_CELL, T_INFO, PHON_SETTINGS, IO)

         IF (ELPH_SETTINGS%DO_CALC()) CALL ELPH_OUTPUT_STUFF(ELPH_INTERPOL, ELPH_SETTINGS, T_INFO, EFERMI, IO)
      ENDIF

!=======================================================================
!  if we are running an ab-initio calculation
!=======================================================================

    IF ((.NOT. LJ_ONLY) .AND. (LDO_AB_INITIO)) THEN
      IF (LWRT_CHGFIT()) THEN
         ! subtract the overlapping atomic charge density from CHTOT
         CALL ATOMIC_CHARGES(T_INFO,LATT_CUR,P,GRIDC,CSTRF,1.0_q,-1.0_q,CHTOT)
         ! write CHGFIT
         CALL WRITE_CHARGE_RC(INFO,T_INFO,LATT_CUR,GRIDC,CHTOT,IO,99)
         ! restore CHTOT
         CALL ATOMIC_CHARGES(T_INFO,LATT_CUR,P,GRIDC,CSTRF,1.0_q,1.0_q,CHTOT)
      ENDIF
      ! calculate stockholder charge partitioning
      CALL STOCKHOLDER_ANALYSIS(INFO,T_INFO,LATT_CUR,P,GRIDC,NCDIJ,CHTOT,DENCOR,IO)

!=======================================================================
!  write out some additional information
!  create the File NMRCAR
!=======================================================================
      IF (LNMRCAR.AND.(EGRAD_CALC_EFG().OR.LCHIMAG)) THEN         
         OPEN(UNIT=99, FILE="NMRCAR.magres", STATUS='UNKNOWN')
     
         CALL WRITE_MAGRES_HEADER(P,LATT_CUR,INFO,T_INFO,XC_DATA,KPOINTS,99)
         WRITE(99, '(A8)') "[magres]"

         ! working with routines from the modules as the required arrays are private
         IF (LCHIMAG) THEN
            CALL NMR_MAGRES_TENS(P,T_INFO, 99)
         ENDIF
         IF (EGRAD_CALC_EFG()) THEN
            CALL EFG_MAGRES_TENS(P,T_INFO,99)
         ENDIF

         WRITE(99, '(A9)') "[/magres]"
         CLOSE(99)
      ENDIF

!=======================================================================
!  write out some additional information
!=======================================================================
!---- write charge density, etc. to CHGCAR, etc. and/or vaspwave.h5
!     if we are interested in any density, write it here:
      CALL WRITE_DENSITIES(IO,XC_DATA,INFO,T_INFO,LATT_CUR,DYN,P,SYMM,&
            GRID,GRID_SOFT,GRIDC,SOFT_TO_C,WDES,W,KINEDEN,CHTOT,DENCOR,RHOLM,CHDEN)
!-----if we are interested in the total (local) potential write it here:
      CALL WRITE_POTENTIAL(IO, KINEDEN, GRID, GRIDC, GRID_SOFT, WDES, &
            INFO, XC_DATA, P, T_INFO, HAMILTONIAN%COULOMB_POT, E, LATT_CUR, CHDEN, CHTOT, CSTRF, &
            CVTOT, DENCOR, SV, HAMILTONIAN, SOFT_TO_C, XCSIF, DYN, &
            LMDIM, N_MIX_PAW, CDIJ, CQIJ)

      IF (IO%LWAP) THEN
         call getwave(EFERMI, W, CDIJ, CQIJ, GRIDC, SV, CVTOT, T_INFO, P, NONL_S, LATT_CUR, IO)
      ENDIF

#ifdef VASP_HDF5
!-----Add additional corrections to the potential here
      IF ( DIP%LCOR_DIP .and. DIP%IDIPCO < 4 ) THEN
         CALL VH5_WRITE_AVERAGE_POTENTIAL(IH5OUTFILEID, LATT_CUR, DIP%IDIPCO,&
         DIP%POTLIN, DIP%VACUUM, GRP_RESULTS, GRP_POTENTIAL) 
      END IF
#endif

      CALL PRINT_CMBJ_AUX(GRIDC,LATT_CUR,XC_DATA)

!=======================================================================
!  Write out the Eigenvalues
!=======================================================================

#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         CALL KPAR_SYNC_CELTOT(WDES,W)
      END IF
#endif

      io_begin
      DO NK=1,KPOINTS%NKPTS
        WRITE(22,*)
        WRITE(22,'(4E15.7)') WDES%VKPT(1,NK),WDES%VKPT(2,NK),WDES%VKPT(3,NK),KPOINTS%WTKPT(NK)
        DO N=1,WDES%NB_TOT
          IF (INFO%ISPIN==1) WRITE(22,852) N,REAL( W%CELTOT(N,NK,1) ,KIND=q), W%FERTOT(N,NK,1)
          IF (INFO%ISPIN==2) &
            WRITE(22,8852) N,REAL( W%CELTOT(N,NK,1) ,KIND=q) ,REAL( W%CELTOT(N,NK,INFO%ISPIN) ,KIND=q), W%FERTOT(N,NK,1), W%FERTOT(N,NK,INFO%ISPIN)
        ENDDO
      ENDDO
      IF (IO%LOPEN) CALL WFORCE(22)
#ifdef VASP_HDF5
      CALL VH5_WRITE_EIGENVAL(IH5OUTFILEID, WDES, W, KPOINTS, KPOINTS_FULL)
#endif
      io_end
      CALL XML_EIGENVALUES(W%CELTOT, W%FERTOT, WDES%NB_TOT, KPOINTS%NKPTS, INFO%ISPIN)

  852 FORMAT(1X,I6,4X,F14.6,2X,F9.6)
 8852 FORMAT(1X,I6,4X,F14.6,2X,F14.6,2X,F9.6,2X,F9.6)
!=======================================================================
!  get current density and write output
!=======================================================================
      CALL CURRENT( W, GRID_SOFT, GRIDC, GRIDUS, C_TO_US, SOFT_TO_C, P, LATT_CUR, &
          HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT, CHTOT, NONLR_S, NONL_S, &
          T_INFO, LMDIM, CRHODE, CDIJ, CQIJ, IRDMAX, IO%IU6, IO%IU0,  IO%NWRITE)

      CALL WRITE_ORBITALMAGOUT(IO%IU6)
      CALL XML_WRITE_ORBITALMAGOUT

!!! test_
!!      CALL NMR_CONVERSE_MAIN( &
!!         HAMILTONIAN, P, KPOINTS, SYMM, INFO, GRID, GRIDC, GRIDUS, C_TO_US, &
!!         WDES, W, NONL_S, NONLR_S, IRDMAX, T_INFO, LATT_CUR, LATT_INI, &
!!         LMDIM, CDIJ, CQIJ, CRHODE, SV, EFERMI, IO)
!!      STOP
!!! test_
!=======================================================================
!  calculate optical matrix elements
!=======================================================================
      !> call compression of bands
      !> this is usefull when doing GW in steps 
      !> however avoid compressing more than once 
      IF (  (MSDGW_NS > 0 .OR. MSDGW_F > 0) .AND. .NOT. MSD_COMPRESSED() )  THEN
         CALL MSD_COMPRESSION( WDES, W, LATT_CUR, INFO, T_INFO, IO )
         CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)
         !> overwrite WAVECAR if WAVECAR.psd has not been written 
         IF( .NOT. MSDGW_LWAVE ) THEN
         IF( .NOT. IO%LH5 .AND. IO%LWAVE ) THEN
            CALL OUTWAV(IO, WDES, W, LATT_INI,  EFERMI)
#ifdef VASP_HDF5
         ELSE IF (IO%LH5 .OR. IO%WRT_DENSITY%LCHARGH5) THEN
            IF (IO%LH5 .AND. IO%LWAVE .AND. ( INFO%LSTOP .OR. LTMP .OR. LHFCALC)) THEN
               CALL VH5_WRITE_WAVEFUNCTIONS(IH5WAVEFILEID, IO, WDES, W, LATT_INI, EFERMI)
            END IF
#endif
         ENDIF
         ENDIF
      ENDIF


      IF (IO%LOPTICS .AND. .NOT. LR_OPTIC_DONE()) THEN
      ! avoid call if already done

        CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)

        CALL START_TIMING("G")
        ! VASP onboard optics
        CALL LR_OPTIC( &
             KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
             T_INFO,INFO,XC_DATA,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
             GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,&
             CHTOT,DENCOR,CVTOT,CSTRF, &
             CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
             CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, &
             LPOT= EXXOEP==0 .AND. INFO%INICHG/=4 .AND.  .NOT. USE_OEP_IN_GW() .AND. .NOT. XC_DATA%LDOMETAGGA)
         ! potential must not be updated for OEP methods
         ! if potential was read from a file INFO%INICHG/=4
         ! or if a meta-GGA is used
        IF (NPAR ==1 .AND. KPAR==1 .AND. LPAW .AND. (.NOT.WDES%LNONCOLLINEAR)) THEN
           ! offboard optics by Juergen Furthmueller
           ALLOCATE(NABIJ(WDES%NB_TOT,WDES%NB_TOT))
           CALL CALC_NABIJ(NABIJ,W,WDES,P,KPOINTS,GRID_SOFT,LATT_CUR, &
                IO,INFO,T_INFO,COMM,IU0,55)
           DEALLOCATE(NABIJ)
        ENDIF
        CALL STOP_TIMING("G",IO%IU6,'OPTICS')
      ENDIF

!=======================================================================
!  Optical matrix elements and dielectric function for core electrons
!=======================================================================

      IF (CH_LSPEC) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*) "Starting Core_Hole_Spec"
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT-Z_CL, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
         CALL CALC_CORE_CON_MAT(T_INFO,IO,KPOINTS,LATT_CUR, &
                                P,WDES,W,SYMM%ISYM)
         IF (IO%IU0>=0) WRITE(TIU0,*) "Ending Core_Hole_Spec"
      ENDIF

!=======================================================================
! calculate electron-phonon interactions using plane-waves
!=======================================================================
      IF (ELPH_PW_SETTINGS%RUN) THEN
         CALL ELECTRON_PHONON_DRIVER(ELPH_PW_SETTINGS, HAMILTONIAN, KPOINTS, GRID, LATT_CUR, LATT_INI, &
                 T_INFO, NONLR_S, NONL_S, W, LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, &
                 GRIDC, GRIDUS, C_TO_US, IRDMAX)
      ENDIF

!=======================================================================
!  Write phelel_params.hdf5 for electron-phonon code
!=======================================================================

      IF (ELPH_DERIV_SETTINGS%GENERATE .AND. CHECKPOINT_FD_SETTINGS%DOES_FINISH()) THEN
         IF (XC_DATA%LTAU) THEN
             CALL ELPHON_DERIV_IO(IO, ELPH_DERIV_SETTINGS, PRIM_CELL, INFO, T_INFO, WDES, CDIJ, CQIJ, &
                                  GRIDC, CVTOT_DERIV, CDIJ_DERIV, FORCE_CONST, POLAR_DATA, MU_DERIV)
         ELSE
             CALL ELPHON_DERIV_IO(IO, ELPH_DERIV_SETTINGS, PRIM_CELL, INFO, T_INFO, WDES, CDIJ, CQIJ, &
                                  GRIDC, CVTOT_DERIV, CDIJ_DERIV, FORCE_CONST, POLAR_DATA)
         ENDIF
      ENDIF

!=======================================================================
!  calculate four orbital integrals
!=======================================================================
      CALL START_TIMING("G")

      CALL TWOELECTRON4O_MAIN( &
      P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
      T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM )

      IF (FOURORBIT/=0) THEN
         ! Set all DFT exchange contributions to zero
         XC_DATA_TEMP=XC_DATA
         XC_DATA_TEMP%ALDAX=1._q
         XC_DATA_TEMP%BLDAX=0._q
         XC_DATA_TEMP%ALDAC=1._q
         XC_DATA_TEMP%AGGAX=0._q
         XC_DATA_TEMP%BGGAX=0._q
         XC_DATA_TEMP%AGGAC=1._q
         XC_DATA_TEMP%AMGGAX=0._q
         XC_DATA_TEMP%AMGGAC=1._q
         XC_DATA_TEMP%AEXX=0._q
         XC_DATA_TEMP%BEXX=1._q
         CALL PUSH_XC_TYPE(XC_DATA,XC_DATA_TEMP)
         ! Initialize xc tables
         CALL SETUP_LDA_XC(XC_DATA,1,IO%IU6,IO%IU0,IO%IDIOT)
         ! calculate correlation contributions on PW grid
         CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,  &
              INFO,XC_DATA,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF,&
              HAMILTONIAN%COULOMB_POT)
         io_begin
            WRITE(*,'(A,F14.7)') 'LDA correlation energy (PW grid):',E%EXCG
         io_end
         ! Restore the original situation
         CALL POP_XC_TYPE(XC_DATA)
         ! reset the xc tables
         IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
            CALL SETUP_LDA_XC(XC_DATA,2,-1,-1,IO%IDIOT)
         ELSE
            CALL SETUP_LDA_XC(XC_DATA,1,-1,-1,IO%IDIOT)
         ENDIF
         ! recalculate plane wave contributions
         CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,  &
              INFO,XC_DATA,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF,&
              HAMILTONIAN%COULOMB_POT)
      ENDIF
      CALL STOP_TIMING("G",IO%IU6,'4ORBIT')

!=======================================================================
! MM: Stuff for Janos Angyan: write the AE-charge density to AECCAR2
!=======================================================================
      IF (LWRT_AECHG()) THEN
         CALL AUGCHG(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        &             LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
        &              LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX,.TRUE.,.FALSE.,IO%IU0)

#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
         OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'AECCAR2',STATUS='UNKNOWN')
         io_begin
         ! write header
         CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
        &             .FALSE., DYN%POSION)
         io_end
         ! write AE charge density
         CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
         CLOSE(99)
         ! write Fourier transform of AE charge density
         IF (LWRTSTRF()) THEN
            OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'STRFAC',STATUS='UNKNOWN')
            CALL WRTSTRF(GRIDC,LATT_CUR,CHTOT,99)
            CLOSE(99)
         ENDIF
!        OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'RADCHG',STATUS='UNKNOWN')
!        CALL WRT_RHO_RAD(WDES,P,T_INFO,INFO%LOVERL,LMDIM,CRHODE,99)
!        CLOSE(99)

         END IF

         ! and restore the total charge density
         CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                  LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
                  LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX)
#ifdef MPI
      ENDIF
#endif
! MM: end of addition
! embedding__
!=======================================================================
! write energy derivative w.r.t. external potential to DERIV
!=======================================================================
      IF (EXTPT_DEXTPOT()) THEN
         CALL EXTPT_DE_DVRN(WDES, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, &
                 P, T_INFO, LMDIM, CRHODE, CHTOT, CHDEN, IRDMAX )
#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
            OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'DERIV',STATUS='UNKNOWN')
            io_begin
            ! write header
            CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
        &             .FALSE., DYN%POSION)
            io_end
            ! write AE charge density
            CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
            CLOSE(99)
#ifdef MPI
         ENDIF
#endif
         ! and restore the total charge density
         CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                  LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
                  LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX)
    ENDIF
! embedding__

!=======================================================================
!  calculate ELF
!=======================================================================
      IF (IO%LELF) THEN
      ALLOCATE(CWORK(GRID_SOFT%MPLWV,WDES%NCDIJ))

      CALL ELF(GRID,GRID_SOFT,LATT_CUR,SYMM,NIOND, W,WDES,  &
               CHDEN,CWORK)

#ifdef MPI
     IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
! write ELF to file ELFCAR
      io_begin
      OPEN(UNIT=53,FILE='ELFCAR',STATUS='UNKNOWN')
      CALL OUTPOS(53,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
      io_end

      DO ISP=1,WDES%NCDIJ
         CALL OUTCHG(GRID_SOFT,53,.FALSE.,CWORK(1,ISP))
      ENDDO

      DEALLOCATE(CWORK)

      do_io CLOSE(53)

      END IF
#ifdef MPI
      ENDIF
#endif
#ifdef MPI
      IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
      IF (WRITE_MOMENTS()) CALL WR_MOMENTS(GRID,LATT_CUR,P,T_INFO,W,WDES,.TRUE.)
      IF (WRITE_DENSITY()) CALL WR_PROJ_CHARG(GRID,P,LATT_CUR,T_INFO,WDES)
#ifdef MPI
      END IF
#endif
      IF (LCALC_ORBITAL_MOMENT().AND.WDES%LNONCOLLINEAR) CALL WRITE_ORBITAL_MOMENT(WDES,T_INFO%NIONS,IO%IU6)
      IF (WDES%LSORBIT) CALL WRITE_SPINORB_MATRIX_ELEMENTS(WDES,T_INFO,IO)

!=======================================================================
!-WAH Write augmentation charges and projectors for differential waves
!=======================================================================

      IF (INFO%LOVERL) THEN
       IF (STM(4) > 0) THEN
#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
          do_io WRITE(*,*) "writing IETS projectors"
          CALL WRT_IETS(LMDIM, WDES%NIONS, WDES%NRSPINORS, CQIJ, WDES, W)
          do_io WRITE(*,*) "IETS projectors  written, exiting"
         END IF
#ifdef MPI
       END IF
#endif
      END IF
!=======================================================================
!  possibly decompose into Bloch states of a given primitive cell
!=======================================================================
!     CALL  KPROJ(IO%IU5, IO%IU0, IO%IU6, GRID, NONL_S, T_INFO, SYMM, P, LATT_CUR, KPOINTS, W)
      CALL  KPROJ(IO%IU5, IO%IU0, IO%IU6, GRID, LATT_CUR, W, SYMM, CQIJ)

!=======================================================================
!  total DOS, calculate ion and lm decomposed occupancies and dos
!=======================================================================
!     IF (JOBPAR/=0  .AND. IO%LORBIT<10 .AND. NPAR /=1) THEN
!        CALL vtutor%write(isAlert, PartialDOS)
!     ELSE IF (JOBPAR/=0 .OR. IO%LORBIT>=10  ) THEN
      IF (JOBPAR/=0 .OR. IO%LORBIT>=10  ) THEN
         DEALLOCATE(PAR,DOSPAR)

         IF (IO%LORBIT==1 .OR. IO%LORBIT==2 .OR. (IO%LORBIT>=11 .AND. IO%LORBIT<=14)) THEN
            LPAR=LMDIMP
         ELSE
            LPAR=LDIMP
         ENDIF

         ALLOCATE(PAR(WDES%NB_TOT,WDES%NKPTS,LPAR,T_INFO%NIONP,WDES%NCDIJ))

         IF (IO%LORBIT>=10) THEN
            CALL SPHPRO_FAST( &
                 GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
                 INFO%LOVERL,LMDIM,CQIJ, LPAR, LDIMP, LMDIMP, FINAL_STEP, IO%LORBIT,PAR, &
                 EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, SYMM)
         ELSE
            CALL SPHPRO( &
                 GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
                 INFO%LOVERL,LMDIM,CQIJ, LPAR, LDIMP, LMDIMP, LTRUNC, IO%LORBIT,PAR)
         ENDIF

         !  get and write partial / projected DOS ...

         ! make space free so that DOSPAR can take this space
         DEALLOCATE(W%CW)         ! make space free so that DOSPAR can take this space
         ALLOCATE (DOSPAR(NEDOS,LPAR,T_INFO%NIONP,WDES%NCDIJ))
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
              NEDOS, LPAR, T_INFO%NIONP, DOS, DOSI, PAR, DOSPAR)

         CALL SYMPDOS(WDES,T_INFO,SYMM,LATT_CUR,LDIMP,LMDIMP,NEDOS,LPAR,T_INFO%NIONP,DOSPAR)

         io_begin
         DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(NEDOS-1)

         WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
         DO I=1,NEDOS
            EN=KPOINTS%EMIN+DELTAE*(I-1)
            WRITE(16,7062) EN,(DOS(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN),(DOSI(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN)
         ENDDO

         DO NI=1,T_INFO%NIONP
            WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
            DO I=1,NEDOS
               EN=KPOINTS%EMIN+DELTAE*(I-1)
               WRITE(16,'(3X,F8.3,64E12.4)') &
                    &            EN,((DOSPAR(I,LPRO,NI,ISP),ISP=1,WDES%NCDIJ),LPRO=1,LPAR)
            ENDDO
         ENDDO

#ifdef VASP_HDF5
         io_begin
         CALL VH5_WRITE_DOS(IH5OUTFILEID, WDES, KPOINTS, DOS, DOSI, DOSPAR, EFERMI, T_INFO%NIONP, LPAR)
         io_end
#endif

         CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .TRUE., &
              DOS, DOSI, DOSPAR, NEDOS, LPAR, T_INFO%NIONP, WDES%NCDIJ)
         CALL XML_PROCAR(PAR, W%CELTOT, W%FERTOT, WDES%NB_TOT, WDES%NKPTS, LPAR ,T_INFO%NIONP,WDES%NCDIJ)
         io_end
      ELSE IF (DYN%IBRION/=1) THEN
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

         io_begin
         DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(NEDOS-1)

         WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
         DO I=1,NEDOS
            EN=KPOINTS%EMIN+DELTAE*(I-1)
            WRITE(16,7062) EN,(DOS(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN),(DOSI(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN)
         ENDDO

#ifdef VASP_HDF5
         CALL VH5_WRITE_DOS(IH5OUTFILEID, WDES, KPOINTS, DOS, DOSI, DOSPAR, EFERMI, T_INFO%NIONP, -1)
#endif
         CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .FALSE., &
              DOS, DOSI, DOSPAR, NEDOS, LPAR, T_INFO%NIONP, WDES%NCDIJ)
         io_end
      ENDIF

      CALL XML_CLOSE_TAG("calculation")
    END IF !LJ_ONLY AND LDO_AB_INITIO

      CALL XML_TAG("structure","finalpos")
      CALL XML_CRYSTAL(LATT_CUR%A, LATT_CUR%B, LATT_CUR%OMEGA)
      CALL XML_POSITIONS(T_INFO%NIONS, DYN%POSION)
      IF (T_INFO%LSDYN) CALL XML_LSDYN(T_INFO%NIONS,T_INFO%LSFOR(1,1))
      IF (DYN%IBRION<=0 .AND. DYN%NSW>0 ) CALL XML_VEL(T_INFO%NIONS, DYN%VEL)
      IF (T_INFO%LSDYN) CALL XML_NOSE(DYN%SMASS)
      CALL XML_CLOSE_TAG("structure")
!=======================================================================
! breath a sigh of relief - you have finished
! this jump is just a jump to the END statement
!=======================================================================
 5100 CONTINUE

      IF (MIX%IMIX==4 .AND. INFO%IALGO.NE.-1) THEN
        CALL CLBROYD(MIX%IUBROY)
      ENDIF

      IF (INFO%LSOFT.OR.INFO%LHARD) THEN
         io_begin
         INQUIRE(FILE="STOPCAR",ERR=5111,EXIST=LTMP)
         IF (LTMP) THEN
            IF (IO%IU0>0) &
            WRITE(TIU0,*) 'deleting file STOPCAR'
            IF (IO%LOPEN) OPEN(99,FILE='STOPCAR',ERR=5111)
            CLOSE(99,STATUS='DELETE',ERR=5111)
         ENDIF
 5111    CONTINUE
         io_end
      ENDIF
#ifdef shmem_bcast_buffer
      IF (LFOCK_SHM_ALLOC) THEN
         CALL M_shmem_dealloc(FOCK_SHM_BCST1)
         CALL M_shmem_dealloc(FOCK_SHM_BCST2)
         CALL M_shmem_dealloc(FOCK_SHM_BCST3)
         LFOCK_SHM_ALLOC=.FALSE.
      ENDIF
#endif
      CALL DUMP_ALLOCATE(IO%IU6)
      CALL DUMP_FINAL_TIMING(IO%IU6)
      CALL STOP_XML

#ifdef ML_AVAILABLE
! Finalize machine-learning program
      IF(ML_LMLFF) THEN
         DO I = 1, ML_TOTNUM_INSTANCES
            IF (ML_SUPER_HANDLE_MAIN(1)%FF%LIB .AND. (ML_SUPER_HANDLE_MAIN(1)%FF%ISTART == 2)) THEN
#ifdef libvaspml
               CALL VASPML%finalize()
#endif
            ELSE
               CALL MACHINE_LEARNING_FINISH (ML_SUPER_HANDLE_MAIN(I),DIR_APP,DIR_LEN,NSTEP)
               CALL MPI_BARRIER(ML_SUPER_HANDLE_MAIN(I)%PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
            ENDIF
         ENDDO
      ENDIF
#endif

#ifdef PROFILING
      ! started in INIT_PROFILING
      PROFILING_STOP('total_time')
      IF (IO%IU0>=0) THEN
        POUT=0
        PTHRESHOLD=0.0_q
        CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'PTHRESHOLD', PTHRESHOLD, IERR, WRITEXMLINCAR)
        IF (PTHRESHOLD/=0._q) POUT=POUT+1
        PROUTINE=''
        CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'PROUTINE', PROUTINE, 40, IERR, WRITEXMLINCAR)
        CALL STRIP(PROUTINE,N,'A')
        CALL LOWER(PROUTINE)
        IF (PROUTINE/='') POUT=POUT+1
        PLEVEL=0
        CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'PLEVEL', PLEVEL, IERR, WRITEXMLINCAR)
        IF (PLEVEL/=0) POUT=POUT+1
        PFLAT=.FALSE.
        CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'PFLAT', PFLAT, IERR, WRITEXMLINCAR)

        IF (PTHRESHOLD/=0.0_q) CALL PRINT_FULL_PROFILE(IO%IU6,PTHRESHOLD=PTHRESHOLD,PFLAT=PFLAT)
        IF (PROUTINE/='') CALL PRINT_ROUTINE_PROFILE(IO%IU6,PROUTINE=PROUTINE,PFLAT=PFLAT)
        IF (PLEVEL/=0) CALL PRINT_LEVEL_PROFILE(IO%IU6,PLEVEL=PLEVEL,PFLAT=PFLAT)
        IF (POUT==0)  CALL PRINT_FULL_PROFILE(IO%IU6,PFLAT=PFLAT)
      ENDIF
#endif

#ifdef VASP_HDF5
5200  CALL STOP_H5(DIR_APP(1:DIR_LEN))
#endif


#if defined(ACC_OFFLOAD) && defined(CUSOLVERMP)
    CALL CUSOLVERMP_FREE()
#endif
#if defined(ACC_OFFLOAD) && defined(CUBLASMP)
    CALL CUBLASMP_FREE()
#endif

      CALL STOP_PLUGINS()

      CALL vtutor%stopCode()

      CONTAINS

!**********************************************************************
!
!>  internal subroutine to perform optimization of electronic
!>  degrees of freedom
!
!**********************************************************************

      SUBROUTINE ELECTRONIC_OPTIMIZATION

      IF ( EXXOEP==1 ) THEN
         CALL ELMIN_LHF( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)
      ELSEIF ( EXXOEP>1) THEN
         CALL ELMIN_OEP( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)
      ELSEIF ( INFO%LONESW ) THEN
         CALL ELMIN_ALL( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)

          ! convergence reached switch back to DIIS
          IF (.NOT. INFO%LABORT .AND. INFO%LONESW_AUTO) THEN
             INFO%LONESW=.FALSE.
          ENDIF

      ELSE
         CALL ELMIN( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)


         IF (INFO%LABORT .AND. INFO%LONESW_AUTO) THEN
          INFO%LONESW=.TRUE.
          W_F%CELTOT=W%CELTOT ! copy current weights to W_F%CELTOT
          MIX%HARD_RESET=.TRUE.
         ! no convergence
         ! try to switch to CG for all bands
          CALL ELMIN_ALL( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)

         ENDIF
      ENDIF

      IF (INFO%LABORT) THEN
          ! Notify the user that electronic convergence was not reached ...
          ! ... we don't stop here, because we still want to do some cleanup
          CALL vtutor%alert("The electronic self-consistency was not achieved in &
              &the given number of steps (NELM). The forces and other quantities &
              &evaluated might not be reliable so examine the results carefully. &
              &If you find spurious results, we suggest increasing NELM, if you &
              &were close to convergence or switching to a different ALGO or &
              &adjusting the density mixing parameters otherwise.", dont_stop=.true.)

          ! Depending on STOP_ON we may want to quit the ionic relaxation or MD
          ! when electronic convergence was not reached in this step  ...
          ! ... we will always quit when INFO%LHARD=T (LABORT=T was read from STOPCAR)
          IF (vtutor%stopOn >= isAlert .OR. INFO%LHARD) THEN
              ! set variables to kill ionic relaxation
              DYN%IBRION = -1
              DYN%NSW = 0
              NSW_ACT = 0
          END IF
      END IF

      END SUBROUTINE ELECTRONIC_OPTIMIZATION

!**********************************************************************
!
!> internal subroutine to perform RPA calculations
!> \todo there are still logical flaws in this routine:
!>    LRPAFORCE will not work if LNBANDS is .TRUE.
!>    since STORE_WDES_GROUNDSTATE is not called
!
!**********************************************************************

      SUBROUTINE DO_RPA
        LOGICAL :: USELDApU_store
        CALL START_TIMING("LOOPX")

        ! need to open XML tag otherwise main will complain upon return (since it closes calculation)
        CALL XML_TAG("calculation")

        ! push the KPOINTS_ORIG to a stack and point KPOINTS_ORIG to KPOINTS
        CALL PUSH_KPOINTS_ORIG
        USELDApU_store=USELDApU()

#ifdef vasp6
         ! if all-in-one mode is selected, NBANDS can be used to reduce
         ! number of bands written to file 
         IF ( LALL_IN_ONE .OR. .NOT. LNBANDS .OR. NBANDSEXACT > 0 ) THEN
            MIX%LRESET=.TRUE.
            ! the first step is to restore the original XC-functional (as read from POTCAR/INCAR)
            ! in order to use the appropriate Hamiltonian
            CALL POP_XC_TYPE(XC_DATA)
            IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
               CALL SETUP_LDA_XC(XC_DATA,2,-1,-1,IO%IDIOT)
            ELSE
               CALL SETUP_LDA_XC(XC_DATA,1,-1,-1,IO%IDIOT)
            ENDIF
            ! now update the PAW one-center terms to current functional
            CALL SET_PAW_ATOM_POT( P, XC_DATA, T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, IO%IU6 )
            ! and restore the convergence corrections
            CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)

            ! optimize orbitals
            CALL ELMIN( &
                 HAMILTONIAN,KINEDEN, &
                 P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
                 T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
                 GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
                 CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
                 CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
                 CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
                 NSTEP,LMDIM,IRDMAX,NEDOS, &
                 TOTEN,EFERMI,LDIMP,LMDIMP)

            ! now we restore the HF xc-correlation functional
            CALL PUSH_XC_TYPE_FOR_GW(XC_DATA)
            CALL SETLDApU(.FALSE.)  ! switch off LDApU
            IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
               CALL SETUP_LDA_XC(XC_DATA,2,-1,-1,IO%IDIOT)
            ELSE
               CALL SETUP_LDA_XC(XC_DATA,1,-1,-1,IO%IDIOT)
            ENDIF
            ! now update the PAW one-center terms to current functional
            CALL SET_PAW_ATOM_POT( P, XC_DATA, T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, IO%IU6 )
            ! and restore the convergence corrections
            CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)

            IF (LRPAFORCE) THEN
               CALL STORE_WDES_GROUNDSTATE(WDES)
            ENDIF
            ! increase number of bands and perform exact diagonalization
            CALL START_TIMING("G")
            CALL EDDIAG_EXACT_UPDATE_NBANDS(HAMILTONIAN,KINEDEN, &
                 GRID,GRID_SOFT,GRIDC,GRIDB,GRIDUS,C_TO_US,SOFT_TO_C,B_TO_C,E, &
                 CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
                 IRDMAX,CRHODE,MIX,N_MIX_PAW,RHOLM,RHOLM_LAST,CHDEN, &
                 LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
                 LMDIM,CDIJ,CQIJ,SV,VDW_SET_GL,T_INFO, DYN, P,IO, XC_DATA, INFO, &
                 XCSIF, EWSIF, TSIF, EWIFOR, TIFOR, PRESS, TOTEN, KPOINTS, EFERMI, NEDOS, USELDApU_store)
            CALL STOP_TIMING("G",IO%IU6,'EDDIAG')

         ELSE
         ! orbitals are not updated
            IF (.NOT. LUSEPEAD() .AND. .NOT. LNABLA .AND. .NOT. EXXOEP>=1 .AND. .NOT. USE_OEP_IN_GW()) THEN
            ! optical properties do not work for conventional GW
            ! because the potential is non-local
            ! OEP are fine however
               IO%LOPTICS=.FALSE.
            ENDIF
            IF ( LNBANDS .AND. LRPAFORCE ) THEN
               CALL VTUTOR%ERROR(" Specific NBANDS set in INCAR is currently not implemented for RPA Forces.")
            ELSE IF ( LNBANDS .AND. LOEP ) THEN
               CALL vtutor%advice("Specific NBANDS set in INCAR should be&
                                 & close to the maximal number of plane-&
                                 &waves for OEP")
            ENDIF


            IF (  (MSDGW_NS > 0 .OR. MSDGW_F > 0) .AND. .NOT. MSD_COMPRESSED() )  THEN
               CALL MSD_COMPRESSION( WDES, W, LATT_CUR, INFO, T_INFO, IO )
               CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)
            ENDIF 
         ENDIF


         ! NBANDS at this stage is the default number of bands used in a KS calculation
         ! (if user doesn't set NBANDS)
         ! For LALL_IN_ONE , the GW routines still use the complete set of possible bands 
         ! determined above 
         ! The follownig matrix shows how many bands are written to WAVECAR 
         !                            |------- number of bands in WAVECAR ------ |
         !  NBANDS_WAVE  LALL_IN_ONE     LNBANDS=F                     LNBANDS=T
         !      -1            T          NBANDSGW                       NBANDS 
         !      -1            F          NBANDSEXACT                    NBANDS 
         !      >0            T     MIN(NBANDS_WAVE,NBANDSEXACT)   MIN( NBANDS_WAVE, NBANDS )
         !      >0            F          NBANDS                    MIN( NBANDS_WAVE, NBANDS )
         ! @TODO: make RPA force calculations compatible with NBANDS setting in INCAR 

         IF (ICHIREAL==1) THEN
            CALL CALCULATE_XI_REAL_GWRK( &
                 KINEDEN, HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
                 GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,E,&
                 CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
                 T_INFO, DYN, INFO, XC_DATA, IO, KPOINTS, SYMM, MIX, &
                 LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, & 
                 EFERMI, TOTEN, NEDOS, DOS, DOSI )
         ELSEIF (ICHIREAL>=1) THEN
            CALL CALCULATE_XI_REAL( &
                 KINEDEN, HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
                 GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,E,&
                 CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
                 T_INFO, DYN, INFO, XC_DATA, IO, KPOINTS, SYMM, MIX, &
                 LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, & 
                 EFERMI, TOTEN, NEDOS, DOS, DOSI )
         ELSE
#endif
            CALL CALCULATE_XI( &
                KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,E,&
                CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
                T_INFO, DYN, INFO, XC_DATA, IO, KPOINTS, SYMM, MIX, &
                LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, &
                EFERMI, NEDOS, DOS, DOSI )
#ifdef vasp6
         ENDIF

         ! Constrained RPA: DMFT_MATRIX_ELEMENTS
         ! TODO put this into chi.F
         IF (LCRPA .AND. ICHIREAL < 1 ) THEN
            ! only real-frequency GW code requires this step
            ! the corresponding matrix elements in space-time GW code are computed
            ! inside the XI_REAL routines for the full imaginary frequency grid
            CALL DMFT_MATRIX_ELEMENTS( &
            &   WDES,W,KPOINTS,GRID,T_INFO,INFO,P, &
            &   NONL_S,NONLR_S,SYMM,LATT_CUR,CQIJ,LMDIM,IO)
         ENDIF
         ! restore the original DFT functional and recalculate CDIJ
         CALL UPDATE_CDIJ(HAMILTONIAN,KINEDEN,GRID,GRID_SOFT, &
              GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,E,CHTOT,DENCOR, &
              CVTOT,CSTRF,IRDMAX,CRHODE,N_MIX_PAW,RHOLM,CHDEN, &
              LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM,LMDIM,CDIJ, &
              CQIJ,SV,T_INFO,P,IO,XC_DATA,INFO,XCSIF,USELDApU_store)

         ! store the RPA density matrix
         CALL RPA_DENSITY( W, LMDIM, IRDMAX, T_INFO, LATT_CUR, P, SYMM, &
              CDIJ, CQIJ, NONLR_S, NONL_S, &
              GRID, GRID_SOFT, GRIDC, GRIDUS, SOFT_TO_C, C_TO_US)

         IF (.NOT. LNBANDS .OR. NBANDSEXACT > 0 ) THEN
            ! restore the original number of orbitals
            CALL RESTORE_NBANDS( WDES, W, LATT_CUR, INFO, LMDIM, T_INFO, P, XC_DATA, IO)
         ENDIF
#endif
         CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)

         CALL KPAR_SYNC_ALL(WDES,W)
         IF ( .NOT. LCRPA .AND. WANNIER90() ) THEN
            CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
         ENDIF

         CALL SEPERATOR_TIMING(IO%IU6)
         CALL STOP_TIMING("LOOPX",IO%IU6,'GWTOTAL')

         IF (LRPAFORCE) THEN
            CALL LR_SKELETON( &
                 KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                 T_INFO,DYN,INFO,XC_DATA,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
                 GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E,&
                 CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
                 CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
                 CHDEN,SV,VDW_SET_GL,DOS,DOSI,CHAM, &
                 8,LMDIM,IRDMAX,NEDOS, &
                 TOTEN,EFERMI, TIFOR, &
                 .FALSE.)

            CALL WRITE_RPA_FORCE( TIFOR, T_INFO, LATT_CUR, SYMM, DYN%POSIOC, DYN%IBRION, IO)
#ifdef VASP_HDF5
            CALL VH5_WRITE_STRESS_FORCES_HIST_STEP(IH5INTERMEDIATEGROUP_ID, GRP_HISTORY, NSTEP/LBLOCK_HELP, &
                                                     TSIF*EVTOJ*1E22_q/LATT_CUR%OMEGA, TIFOR)
#endif
            CALL SEPERATOR_TIMING(IO%IU6)
            CALL STOP_TIMING("LOOPX",IO%IU6,'LINEAR')
         ENDIF

         ! restore the original KPOINTS from stack
         CALL POP_KPOINTS_ORIG

         ! now we restore the HF xc-correlation functional (but leave LDApU at INCAR level)
         CALL PUSH_XC_TYPE_FOR_GW(XC_DATA)
         IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
            CALL SETUP_LDA_XC(XC_DATA,2,-1,-1,IO%IDIOT)
         ELSE
            CALL SETUP_LDA_XC(XC_DATA,1,-1,-1,IO%IDIOT)
         ENDIF
         ! now update the PAW one-center terms to current functional
         CALL SET_PAW_ATOM_POT( P, XC_DATA, T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, IO%IU6 )
         ! and restore the convergence corrections
         CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)

         !> free optical matrix elements to recalculate them, if requested, after GW step
         !> however, do this only if pead is selected 
         IF(ALLOCATED( CDER_BETWEEN_STATES ) .AND. ( IO%LOPTICS .AND. LUSEPEAD() )  ) THEN
            DEALLOCATE( CDER_BETWEEN_STATES ) 
         ELSEIF(ALLOCATED( CDER_BETWEEN_STATES ) .AND. ( IO%LOPTICS .AND. .NOT. LUSEPEAD() )  ) THEN
            IF( LGW ) CALL vtutor%write(isAlert, LOPTICSAndGW)
         ENDIF


    END SUBROUTINE DO_RPA


!**********************************************************************
!
!> internal subroutine to perform a few performance tests
!> Output is written to IUT
!
!**********************************************************************

    SUBROUTINE PERFORMANCE_TEST
      io_begin
      IUT=IO%IU0

      IF (IUT>0) WRITE(IUT,5001)
 5001 FORMAT(/ &
     & ' All results refer to a run over all bands and one k-point'/ &
     & ' VNLACC   non local part of H'/ &
     & ' PROJ     calculate projection of all bands (contains FFTWAV)'/ &
     & ' RACC     non local part of H in real space (contains FFTEXT)'/ &
     & ' RPRO     calculate projection of all bands in real space '/ &
     & '          both calls contain on FFT (to be subtracted)'/ &
     & ' FFTWAV   FFT of a wavefunction to real space'/ &
     & ' FFTEXT   FFT to real space'/ &
     & ' ECCP     internal information only (subtract FFTWAV)'/ &
     & ' POTLOK   update of local potential (including one FFT)'/ &
     & ' SETDIJ   calculate stregth of US PP'/ &
     & ' ORTHCH   gramm-schmidt orth.  applying Choleski decomp.'/ &
     & ' LINCOM   unitary transformation of wavefunctions'/ &
     & ' LINUP    upper triangle transformation of wavefunctions'/ &
     & ' ORTHON   orthogonalisation of one band to all others')
       io_end


! set the wavefunction descriptor
      ISP=1
      NK=1
      CALL SETWDES(WDES,WDES1,NK)

      INFO%ISPIN=1
      INFO%RSPIN=2

      NPL=WDES%NPLWKP(NK)
      ALLOCATE(CWORK1(GRID%MPLWV),CWORK2(GRID%MPLWV),CPROTM(LMDIM*NIOND))
      W1%CR=>CWORK1

      CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror ))
#ifdef MPI
      IF (WDES%COMM_INTER%NCPU/=1) THEN
      CALL START_TIMING("G")
      DO I=1,10
         CALL REDIS_PW(WDES1, WDES%NBANDS, W%CW   (1,1,NK,1))
      ENDDO
      CALL STOP_TIMING("G",IUT,'10xRED')
      ENDIF
#endif

      IF(INFO%TURBO==0)THEN
         CALL START_TIMING("G")
         CALL STUFAK(GRIDC,T_INFO,CSTRF)
         CALL STOP_TIMING("G",IUT,'STUFAK')
      ENDIF

      CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES, &
           INFO,XC_DATA,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF,&
           HAMILTONIAN%COULOMB_POT)
      CALL STOP_TIMING("G",IUT,'POTLOK ')

      CALL FORLOC(GRIDC,P,T_INFO,LATT_CUR, CHTOT,TIFOR)
      CALL STOP_TIMING("G",IUT,'FORLOC')

      CALL SOFT_CHARGE(GRID,GRID_SOFT,W,WDES, CHDEN(1,1))
      CALL STOP_TIMING("G",IUT,'CHSP  ')


      CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
               LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
               LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX)

      CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)
      CALL STOP_TIMING("G",IUT,'DEPLE ')

      IF (INFO%LREAL) THEN

        CALL RSPHER(GRID,NONLR_S,LATT_CUR)
        CALL STOP_TIMING("G",IUT,'RSPHER')
        CWORK2=0
        CALL START_TIMING("G")
        CALL RACCT(NONLR_S,WDES,W,GRID,CDIJ,CQIJ, ISP, LMDIM, NK)
        CALL STOP_TIMING("G",IUT,'RACC')

      ELSE

        CALL PHASE(WDES,NONL_S,NK)
        NPL=WDES%NPLWKP(NK)
        CALL START_TIMING("G")
        DO N=1,WDES%NBANDS
          EVALUE=W%CELEN(N,1,1)
          CALL SETWAV(W,W1,WDES1,N,1)  ! allocation for W1%CR done above
          CALL VNLACC(NONL_S,W1,CDIJ,CQIJ, ISP, EVALUE,  CWORK2)
        ENDDO
        CALL STOP_TIMING("G",IUT,'VNLACC')
      ENDIF


      CALL START_TIMING("G")
      IF (INFO%LREAL) THEN
        CALL START_TIMING("G")
        CALL RPRO(NONLR_S,WDES,W,GRID,NK)
        CALL STOP_TIMING("G",IUT,'RPRO  ')
      ELSE
        CALL PROJ(NONL_S,WDES,W,NK)
        CALL STOP_TIMING("G",IUT,'PROJ  ')
      ENDIF

      CALL START_TIMING("G")
      DO  N=1,WDES%NBANDS
        CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CWORK1,W%CPTWFP(1,N,NK,1),GRID)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFTWAV')

      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RC%NP,CWORK1)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFTINI')


      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RC%NP,CWORK1)
        CALL FFT3D(CWORK1,GRID,1)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFT3DF')

      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RC%NP,CWORK1)
        CALL FFT3D(CWORK1,GRID,1)
        CALL FFT3D(CWORK1,GRID,-1)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFTFB ')

      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RL%NP,CWORK1)
        CALL FFTEXT(NPL,WDES%NINDPW(1,NK),CWORK1,CWORK2,GRID,.FALSE.)
      ENDDO

      CALL STOP_TIMING("G",IUT,'FFTEXT ')

      DO N=1,WDES%NBANDS
          CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CWORK1,W%CPTWFP(1,N,NK,1),GRID)
          CALL SETWAV(W,W1,WDES1,N,1)  ! allocation for W1%CR done above
          IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
             CALL ECCP_TAU(WDES1,W1,W1,LMDIM,CDIJ,GRID,SV,LATT_CUR,HAMILTONIAN%MU,W1%CELEN)
          ELSE
             CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ,GRID,SV, W1%CELEN)
          ENDIF
      ENDDO

      CALL STOP_TIMING("G",IUT,'ECCP ')

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL STOP_TIMING("G",IUT,'SETDIJ ')

      CALL SET_DD_PAW(WDES, P, XC_DATA, T_INFO, INFO%LOVERL, &
         INFO%ISPIN, LMDIM, CDIJ,  RHOLM, CRHODE, &
          E, LASPH =INFO%LASPH , LCOREL=.FALSE.)

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,XC_DATA,IO%IU6)

      CALL STOP_TIMING("G",IUT,'SETPAW ')

      CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
      CALL REDIS_PW_OVER_BANDS(WDES, W)
      CALL STOP_TIMING("G",IUT,'ORTHCH')

      IF (INFO%LDIAG) THEN
        IFLAG=3
      ELSE
        IFLAG=4
      ENDIF
      CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF)

      CALL REDIS_PW_OVER_BANDS(WDES, W)
      CALL STOP_TIMING("G",IUT,'EDDIAG')

! avoid that MATMUL is too clever
      ALLOCATE(CMAT(WDES%NB_TOT,WDES%NB_TOT))
      DO N1=1,WDES%NB_TOT
      DO N2=1,WDES%NB_TOT
        IF (N1==N2)  THEN
         CMAT(N1,N2)=0.99999_q
        ELSE
         CMAT(N1,N2)=EXP((0.7_q,0.5_q)/100)
       ENDIF
      ENDDO
      ENDDO

      NPRO= WDES%NPRO
      CALLMPI( SET_NPL_NPRO(WDES1, NPL, NPRO))
#ifdef MPI
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
      NCPU=1
#endif
      NRPLWV_RED=WDES%NRPLWV/NCPU
      NPROD_RED =WDES%NPROD /NCPU

      CALL START_TIMING("G")
      CALL LINCOM('F',W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1),CMAT(1,1), &
       WDES%NB_TOT,WDES%NB_TOT,NPL,0,NRPLWV_RED,NPROD_RED,WDES%NB_TOT, &
       W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1))
      CALL STOP_TIMING("G",IUT,'LINCOM')

      CALL LINCOM('F',W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1),CMAT(1,1), &
       WDES%NB_TOT,WDES%NB_TOT,NPL,NPRO,NRPLWV_RED,NPROD_RED,WDES%NB_TOT, &
       W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1))
      CALL STOP_TIMING("G",IUT,'LINCOM')

      CALL LINCOM('U',W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1),CMAT(1,1), &
       WDES%NB_TOT,WDES%NB_TOT,NPL,NPRO,NRPLWV_RED,NPROD_RED,WDES%NB_TOT, &
       W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1))
      CALL STOP_TIMING("G",IUT,'LINUP ')
    END SUBROUTINE PERFORMANCE_TEST


#ifdef MinimaxTest
!**********************************************************************
!
!> Minimax integration test routine
!
!**********************************************************************
     SUBROUTINE MM_INT_TEST
        REAL(q) :: E_MIN, E_MAX
        CALL DENSTA_SPIN( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
           INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
           NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

        E_MIN = MINVAL( ABS( W%CELTOT(:,:,:)-EFERMI ) )
        E_MAX = MAXVAL( ABS( W%CELTOT(:,:,:)-EFERMI ) )
        IF ( OMEGATL < 0 ) OMEGATL = -1.5*E_MAX

        IF (IO%IU6>=0) THEN
           WRITE(IO%IU6,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
           & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
           & ' E-fermi : ', 2E20.16,/)" ) EFERMI
           WRITE(IO%IU0,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
           & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
           & ' E-fermi : ', 2F8.4,/)" ) EFERMI

           IF ( LFINITE_TEMPERATURE ) THEN
              WRITE(IO%IU0,"(' E-min, E-max, E-max*beta: ',2F12.5,/)") E_MIN,E_MAX, E_MAX/KPOINTS%SIGMA
           ELSE
              WRITE(IO%IU0,"(' E-min, E-max, R: ',3F12.5,/)") E_MIN,E_MAX*E_MIN, E_MAX
           ENDIF
           IF( OMEGATL < 0 ) &
           WRITE(IO%IU0,"(' WARNING OMEGATL set automatically to 1.5*E-max',2F12.5,/)") ABS( OMEGATL )
        ENDIF
        IF ( OMEGATL < 0 ) OMEGATL = ABS( OMEGATL )

        ! determine grid
        CALL SET_IMAG_GRID_HANDLE(OMEGAMIN, OMEGATL, NOMEGA, OMEGAGRID, &
            KPOINTS%SIGMA, IMAG_GRIDS, WDES%COMM, IO)
        ! perform Integration test
        CALL MINIMAX_INTEGRATION_TEST(IMAG_GRIDS, W, EFERMI, IO)
        ! deallocate handle
        CALL DEALLOCATE_IMAG_GRID_HANDLE( IMAG_GRIDS )
     END SUBROUTINE
#endif

#ifdef GDBDEBUG
!**********************************************************************
!
!> Used for debugging, a file PAUSE_VASP is written containing all 
!> running PIDs. While this file is present VASP waits for debugger 
!> to attach to these PIDs. Once attached, PAUSE_VASP can be deleted 
!> and execution of program can be continued inside debugger.
!> Probably best used with tmpi-debug 
!
!**********************************************************************
  SUBROUTINE WAIT_ON_DEBUGGER()
    USE main_mpi
    IMPLICIT NONE
    INTERFACE
    FUNCTION GETPID() BIND(C, NAME='getpid')
       USE ISO_C_BINDING
       INTEGER(C_INT) :: GETPID
    END FUNCTION
    END INTERFACE
    INTEGER :: I
    INTEGER :: NODE_ME = 1
    INTEGER :: NCPU = 1
#ifdef MPI
    INTEGER :: PIDS( COMM_WORLD%NCPU )
    NODE_ME = COMM_WORLD%NODE_ME
    NCPU = COMM_WORLD%NCPU
    PIDS = 0 
    PIDS(COMM_WORLD%NODE_ME ) = getpid()
    CALLMPI( M_sum_i( COMM_WORLD, PIDS(1), COMM_WORLD%NCPU ) )
#else
    INTEGER :: PIDS(1)
    PIDS(1) = getpid()
#endif
    IF (NODE_ME == 1) THEN
       OPEN(100, FILE="PAUSE_VASP", STATUS="REPLACE" )
       DO I = 1, NCPU
          WRITE(100,'( I12 )' )PIDS(I)
       ENDDO
       CLOSE(100)
       WRITE(*,*) "use tmpi-debug or attach debugger to following PIDs and remove PAUSE_VASP manually:", PIDS
    ENDIF
    !> periodically check if PAUSE_VASP is present 
    DO WHILE ( PAUSE_VASP_PRESENT() )
       CALL SLEEP(1)
    END DO
    ! halt program one last time and continue execution in debugger 
    IF (NODE_ME == 1) THEN
       WRITE(*,*) "continue execution in debugger"
    ENDIF
    CALL SLEEP(3)

  END SUBROUTINE WAIT_ON_DEBUGGER
 
  !> checks if file PAUSE_VASP is present
  FUNCTION PAUSE_VASP_PRESENT( ) 
     LOGICAL PAUSE_VASP_PRESENT
     LOGICAL LEXIST
     INQUIRE( FILE="PAUSE_VASP", EXIST = LEXIST )
     PAUSE_VASP_PRESENT = LEXIST
  END FUNCTION PAUSE_VASP_PRESENT

#endif

  END PROGRAM
