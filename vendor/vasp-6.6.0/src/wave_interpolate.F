#include "symbol.inc"

!************************ MODULE WAVE_INTERPOLATE **********************
!>
!> @brief
!> this module provides functions to interpolate the band structure
!> to an arbitrary k-point grid
!>
!> The function wave_interpolate::WAVE_INTERPOLATOR_INIT keeps a pointer to the original WF
!> and created a new one WAVE_INTERPOLATOR%W with double the size of the
!> original wavefunctions
!>
!> The function wave_interpolate::WAVE_INTERPOLATOR_KPOINTS can be called interpolating
!> the wavefunctions to a specified batch of k-points
!> only the new wavefunctions are optimized
!>
!> INTERPOLATE_WAVE_FREE restores frees WAVE_INTERPOLATOR%w and WAVE_INTERPOLATOR%WDES
!>
!> This routine was writted by Henrique Miranda
!> Based on the INTERPOLATE_BAND_STR routine by gK and
!> INTERPOLATE_BANDSTR_GRDSHFT by eFL
!> @authors Henrique Miranda, eFL, gK
!>
!***********************************************************************

MODULE wave_interpolate

  USE prec
  USE wave_high
  USE lattice
  USE poscar
  USE mpimy
  USE mgrid
  USE nonl_high
  USE base
  USE pseudo
  USE kpoints_change
  USE constant
  USE choleski
  USE subrot
  USE david
  USE rmm_diis
  USE hamil_struct_def
  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  TYPE wave_interpolator_settings
     !> Mode to use for optional k-point grid
     INTEGER :: KPOINTS_OPT_MODE=1
     !> Logical to activate use of optional k-point grid
     LOGICAL :: LKPOINTS_OPT=.FALSE.
     !> Size of the number of k-points in each batch
     INTEGER :: KPOINTS_OPT_NKBATCH=-1
     !> Kspacing for kpoints_opt
     REAL(q) :: KSPACING_OPT=0.5
  END TYPE wave_interpolator_settings

  TYPE wave_interpolator
     !> Number of spin components
     INTEGER :: ISPIN
     !> Number of orbitals
     INTEGER :: NB_TOT
     !> Number of k-points in the original wavefunction datastructure
     INTEGER :: NKPTS_ORIG
     !> Max number of k-points in every batch
     INTEGER :: MAX_NKPTS_BATCH
     !> Number of k-points in batch
     INTEGER :: NKPTS_BATCH
     !> Size of the first-dimensions of the wavefunction.
     !> This is useful to ensure the dimensions don't change when computing for
     !> different k-point. The default is a negative number which means the size
     !> is determined by GEN_LAYOUT and can change depending on which k-points are computed
     INTEGER :: NGDIM
     !> Coordinates of k-points in batch
     REAL(q),ALLOCATABLE :: KPOINTS_BATCH(:,:)
     !> index of k-points in batch
     INTEGER,ALLOCATABLE :: KPOINTS_INDEX(:)
     !> option to decide with method is used to find intial guess for the wavefunctions
     !>  =1 Use WF in close k-point as initial guess
     !>  =2 Use derivative of the wavefunction to interpolate to new k-points
     INTEGER :: OPTION
     !> original grid structure
     TYPE (grid_3d),POINTER :: GRID_ORIG
     TYPE (grid_3d) :: GRID
     !> original kpoints structure
     TYPE (kpoints_struct),POINTER :: KPOINTS_ORIG
     !> Kpoints structure
     TYPE (kpoints_struct) :: KPOINTS
     !> Wavefunction container for computation of derivatives of the WFs
     TYPE (wavespin) :: WDIFF3(3)
     !> Reference to the original WF
     TYPE (wavespin),pointer :: W_ORIG
     !> Reference to the original wavefunction descriptor
     TYPE (wavedes),pointer :: WDES_ORIG
     !> wavefunction container
     TYPE (wavespin) :: W
     !> wavefunction descriptor
     TYPE (wavedes) :: WDES
     !> Logical to know wether a wavefunction and wavefunction descriptor are allocated
     LOGICAL :: HAS_W=.FALSE.
     !> Logical to know whether we have the derivative of the orbitals w.r.t or not
     LOGICAL :: HAS_WAVEDER = .FALSE.
     !> Converge all the bands and not only the occupied states
     LOGICAL :: FULL_ACC = .FALSE.
  END TYPE wave_interpolator

  !> Map the k-points in a WF object
  TYPE wave_mapper
     !> Number of k-points in the original wavefunction datastructure
     INTEGER :: NKPTS_ORIG
     !> Coordinates of k-points in batch
     REAL(q),ALLOCATABLE :: KPOINTS_BATCH(:,:)
     !> index of k-points in batch
     INTEGER,ALLOCATABLE :: KPOINTS_INDEX(:)
     !> Number of k-points in batch
     INTEGER :: NKPTS_BATCH
  END TYPE wave_mapper

  PUBLIC :: WAVE_INTERPOLATOR_INIT
  PUBLIC :: WAVE_INTERPOLATOR_KPOINTS
  PUBLIC :: WAVE_INTERPOLATOR_FREE

CONTAINS

!***********************************************************************
!> Read input variables controling the behaviour of the wave_interpolator
!***********************************************************************
SUBROUTINE WAVE_INTERPOLATOR_READER(INCAR,WAVE_INTERP_SET,IO)
  USE main_mpi, only: COMM_INB
  USE incar_reader, only: process_incar, incar_file, count_elements
  USE pead, ONLY: PEAD_REQUEST
  USE fock, ONLY: LHFCALC

  TYPE(incar_file),INTENT(INOUT) :: INCAR
  TYPE(wave_interpolator_settings),INTENT(OUT) :: WAVE_INTERP_SET
  TYPE(in_struct),INTENT(IN) :: IO

  ! local
  INTEGER IERR
  REAL(q) RDUM
  COMPLEX(q) CDUM
  LOGICAL LOPEN,LDUM

  !
  ! interpolation in k using linear response routines
  !
  WAVE_INTERP_SET%KPOINTS_OPT_MODE=1
  CALL PROCESS_INCAR(INCAR, 'KINTEROPT', WAVE_INTERP_SET%KPOINTS_OPT_MODE, IERR)
  IF (IERR==0) call vtutor%alert('The KINTEROPT INCAR tag is deprecated, use KPOINTS_OPT_MODE instead.')
  CALL PROCESS_INCAR(INCAR, 'KPOINTS_OPT_MODE', WAVE_INTERP_SET%KPOINTS_OPT_MODE, IERR)
  IF (WAVE_INTERP_SET%KPOINTS_OPT_MODE/10==1.OR.WAVE_INTERP_SET%KPOINTS_OPT_MODE/10==2) CALL PEAD_REQUEST

  !
  ! additional finer grid
  !
  INQUIRE(FILE="KPOINTS_OPT",EXIST=WAVE_INTERP_SET%LKPOINTS_OPT)
  CALL PROCESS_INCAR(INCAR, 'LKPOINTS_OPT', WAVE_INTERP_SET%LKPOINTS_OPT, IERR)

  !
  ! set number of k-points in the additional k-points batch
  !
  WAVE_INTERP_SET%KPOINTS_OPT_NKBATCH=-1
  CALL PROCESS_INCAR(INCAR, 'KPOINTS_OPT_NKBATCH', WAVE_INTERP_SET%KPOINTS_OPT_NKBATCH, IERR)

  !
  ! read KSPACING_OPT
  !
  WAVE_INTERP_SET%KSPACING_OPT=0.5_q
  IF (COUNT_ELEMENTS(INCAR, 'KSPACING_OPT')>0) WAVE_INTERP_SET%LKPOINTS_OPT=.TRUE.
  CALL PROCESS_INCAR(INCAR, 'KSPACING_OPT', WAVE_INTERP_SET%KSPACING_OPT, IERR)

  ! stop here if hybrid functionals with NCORE/=1 was requested
#ifdef MPI
  IF (WAVE_INTERP_SET%LKPOINTS_OPT.AND.(COMM_INB%NCPU/=1).AND.LHFCALC) THEN
     call vtutor%error('A calculation of the band-structure with hybrid-functionals on the k-points read '//&
                       'from KPOINTS_OPT file was requested with NCORE/=1. This is currently not supported.')
  ENDIF
#endif

END SUBROUTINE WAVE_INTERPOLATOR_READER

!***********************************************************************
!> Write information about wave_interpolator to the OUTCAR file
!***********************************************************************
SUBROUTINE WAVE_INTERPOLATOR_WRITE_OUTCAR(WAVE_INTERP,IU6)
  TYPE (wave_interpolator_settings),INTENT(IN) :: WAVE_INTERP
  INTEGER,INTENT(IN) :: IU6               !> output unit

  IF (IU6>=0) THEN
     WRITE(IU6,*)          ' Optional k-point grid parameters'
     WRITE(IU6,'(A,L4,A)') '   LKPOINTS_OPT  =  ',WAVE_INTERP%LKPOINTS_OPT,     '    use optional k-point grid'
     WRITE(IU6,'(A,I6,A)') '   KPOINTS_OPT_MODE=',WAVE_INTERP%KPOINTS_OPT_MODE, '    mode for optional k-point grid'
     WRITE(IU6,*)
  ENDIF

END SUBROUTINE WAVE_INTERPOLATOR_WRITE_OUTCAR

!***********************************************************************
!> Write information about the optinal k-points file to the xml file
!***********************************************************************
SUBROUTINE XML_WRITE_WAVE_INTERPOLATOR(WAVE_INTERP_SET)
  USE vaspxml
  TYPE(wave_interpolator_settings),INTENT(IN) :: WAVE_INTERP_SET

  !local varaibles
  LOGICAL :: LDUM
  INTEGER :: IDUM
  REAL(q) :: RDUM
  COMPLEX(q)  :: CDUM
  CHARACTER(1) :: CHARAC

  CALL XML_TAG("separator","optional k-points parameters")
  CALL XML_INCAR('KPOINTS_OPT_MODE','I',WAVE_INTERP_SET%KPOINTS_OPT_MODE,RDUM,CDUM,LDUM,CHARAC,1)
  CALL XML_INCAR('LKPOINTS_OPT','L',IDUM,RDUM,CDUM,WAVE_INTERP_SET%LKPOINTS_OPT,CHARAC,1)
  CALL XML_CLOSE_TAG
END SUBROUTINE XML_WRITE_WAVE_INTERPOLATOR

!***********************************************************************
!>
!> @brief Compute the derivative of the wavefunction with respect to k
!> this can be later on used to improve the intial guess of the wavefunction
!>
!***********************************************************************
SUBROUTINE WAVE_INTERPOLATOR_WAVEDER(WAVE_INTERP,HAMILTONIAN,W,KPOINTS,GRID,LATT_INI,LATT_CUR,T_INFO,&
                  NONLR_S,NONL_S,LMDIM,GRIDC,GRIDUS,C_TO_US,IRDMAX,P,CQIJ,CDIJ,SV,SYMM,INFO,IO)
  USE lattice, ONLY: latt
  USE subrot_cluster
  USE hamil_struct_def
  USE mlr_optic, ONLY: FOCK_K_DER_ANALYT
  USE mlrf_main, ONLY: LRF_RPHI0
  TYPE (wave_interpolator) :: WAVE_INTERP
  TYPE (ham_handle)     :: HAMILTONIAN
  TYPE (wavespin)       :: W
  TYPE (kpoints_struct) :: KPOINTS
  TYPE (grid_3d)        :: GRID
  TYPE (latt)           :: LATT_INI, LATT_CUR
  TYPE (potcar)         :: P(:)
  TYPE (grid_3d)        :: GRIDC
  TYPE (grid_3d)        :: GRIDUS
  TYPE (transit)        :: C_TO_US
  TYPE (info_struct)    :: INFO
  TYPE (type_info)      :: T_INFO
  TYPE (nonlr_struct)   :: NONLR_S
  TYPE (nonl_struct)    :: NONL_S
  TYPE (in_struct)      :: IO
  TYPE (symmetry)       :: SYMM
  INTEGER               :: LMDIM, IRDMAX
  OVERLAP :: CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
  OVERLAP :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
  RGRID   :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)
  ! local variables
  TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
  INTEGER :: IDIR
  LOGICAL :: LDONE

  ! Find degenerate energies
  NULLIFY(DEG_CLUSTER)
  CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER)

  ! Compute derivative of WFs
  DO IDIR=1,3
     CALL ALLOCW(W%WDES,WAVE_INTERP%WDIFF3(IDIR))
     WAVE_INTERP%WDIFF3(IDIR)%CW   =0
     WAVE_INTERP%WDIFF3(IDIR)%CPROJ=0

     ! We use linear response only (to allow arbitrary grids)
     ! TODO: The derivatives for the original points are always recomputed
     CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, &
          T_INFO,  NONLR_S, NONL_S, W, WAVE_INTERP%WDIFF3(IDIR), LMDIM, P, CQIJ, &
          SYMM, IDIR, LDONE, IO%IU0, IO%IU6)

     CALL LRF_RPHI0(HAMILTONIAN,P,NONLR_S,NONL_S,W,LATT_CUR, &
          T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
          CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, WAVE_INTERP%WDIFF3(IDIR), LDONE)
  ENDDO
  CALL FREE_DEG_CLUSTERS(W%WDES, DEG_CLUSTER)
  WAVE_INTERP%HAS_WAVEDER = .TRUE.

END SUBROUTINE WAVE_INTERPOLATOR_WAVEDER

!***********************************************************************
!>
!> @brief Initialize the wavefunction interpolator structure
!> increase the size of the wavefunction
!> in the old part the original wavefunctions are stored
!> the new part of the wavefunction has the size of the batches
!> of k-points that will be used
!>
!***********************************************************************

SUBROUTINE WAVE_INTERPOLATOR_INIT(WAVE_INTERP_SET, WAVE_INTERP, KPOINTS, VKPT_INTER, GRID, &
                  LATT_CUR, LATT_INI, T_INFO, NONL_S, W, P, &
                  SYMM, IO, MAX_NKPTS_BATCH, NGDIM, NB_TOT, USE_ORIGINAL_GRID)

  use fock, only: LHFCALC
  TYPE (wave_interpolator_settings),INTENT(IN) :: WAVE_INTERP_SET
  TYPE (wave_interpolator),TARGET :: WAVE_INTERP
  TYPE (grid_3d),TARGET, INTENT(IN) :: GRID
  TYPE (kpoints_struct),TARGET :: KPOINTS !< original set of kpoints (normally a regular grid)
  REAL(q) :: VKPT_INTER(:,:) !< coordinates of the k-points where we will want to interpolate quantities
  TYPE (latt)        LATT_CUR
  TYPE (latt)        LATT_INI
  TYPE (type_info)   T_INFO
  TYPE (nonl_struct) NONL_S
  TYPE (wavespin), TARGET, INTENT(IN) :: W !< original WF (we should not change anything here)
  TYPE (symmetry) :: SYMM
  TYPE (potcar) :: P(:)
  TYPE (wavespin) :: WDIFF3(3)
  TYPE (in_struct) :: IO
  INTEGER, OPTIONAL :: MAX_NKPTS_BATCH !< Select the size of the k-points batches
  INTEGER, OPTIONAL :: NGDIM !< Enforce the first dimensions of the wavefunction to take this value
  INTEGER, OPTIONAL :: NB_TOT !< Number of bands to interpolate
  LOGICAL, OPTIONAL :: USE_ORIGINAL_GRID !< Logical to control wether to use the original grid in the wave interpolator subroutines.
                                         !< Care should be taken when NCORE/=1 and the requested k-points are different from the
                                         !< ones used to generate the original GRID.
  ! local variables
  INTEGER IK
  INTEGER NKPTS_ORIG
  LOGICAL MY_USE_ORIGINAL_GRID
  REAL(q),ALLOCATABLE :: WTKPT_INTER(:)

  IF (IO%IU0>=0) WRITE(IO%IU0,*) "Start KPOINTS_OPT (optional k-point list driver)"

  ! save reference to the original WF
  WAVE_INTERP%W_ORIG => W
  ! save reference to the original WF descriptor
  WAVE_INTERP%WDES_ORIG => W%WDES
  ! save reference to the original KPOINTS
  WAVE_INTERP%KPOINTS_ORIG => KPOINTS
  ! save reference to the original GRID
  WAVE_INTERP%GRID_ORIG => GRID

  ! save some dimensions to the wave_interp structure.
  ! these must be consistent between the original and the new wavefunction object
  WAVE_INTERP%ISPIN = W%WDES%ISPIN
  WAVE_INTERP%NB_TOT = W%WDES%NB_TOT
  IF (PRESENT(NB_TOT)) WAVE_INTERP%NB_TOT=NB_TOT

  ! original number of k-points
  NKPTS_ORIG = 0
  IF (LHFCALC) NKPTS_ORIG = W%WDES%NKPTS
  ! can change this to any number to configure the size of the batches
  WAVE_INTERP%MAX_NKPTS_BATCH = W%WDES%NKPTS
  IF (WAVE_INTERP_SET%KPOINTS_OPT_NKBATCH>0) WAVE_INTERP%MAX_NKPTS_BATCH = WAVE_INTERP_SET%KPOINTS_OPT_NKBATCH
  IF (PRESENT(MAX_NKPTS_BATCH)) WAVE_INTERP%MAX_NKPTS_BATCH = MAX_NKPTS_BATCH
  WAVE_INTERP%NGDIM=-1; IF (PRESENT(NGDIM)) WAVE_INTERP%NGDIM=NGDIM

  ! Resize k-point structure
  CALL COPY_KPOINTS(KPOINTS,WAVE_INTERP%KPOINTS)
  WAVE_INTERP%KPOINTS%NKDIM = NKPTS_ORIG+WAVE_INTERP%MAX_NKPTS_BATCH
  WAVE_INTERP%KPOINTS%NKPTS = NKPTS_ORIG+WAVE_INTERP%MAX_NKPTS_BATCH
  DEALLOCATE(WAVE_INTERP%KPOINTS%WTKPT)
  DEALLOCATE(WAVE_INTERP%KPOINTS%VKPT)
  ALLOCATE(WAVE_INTERP%KPOINTS%VKPT(3,WAVE_INTERP%KPOINTS%NKPTS))
  ALLOCATE(WAVE_INTERP%KPOINTS%WTKPT(WAVE_INTERP%KPOINTS%NKPTS))
  WAVE_INTERP%KPOINTS%VKPT=0.0_q
  WAVE_INTERP%KPOINTS%WTKPT=0.0_q
  WAVE_INTERP%KPOINTS%VKPT(:,:NKPTS_ORIG)=WAVE_INTERP%KPOINTS_ORIG%VKPT(:,:NKPTS_ORIG)
  WAVE_INTERP%KPOINTS%WTKPT(:NKPTS_ORIG)=WAVE_INTERP%KPOINTS_ORIG%WTKPT(:NKPTS_ORIG)

  ! create a new WF with the original k-points + the new k-points that we want to compute
  CALL WDES_FROM_WDES(WAVE_INTERP%WDES,WAVE_INTERP%W_ORIG%WDES)

  ! create a new grid compatible with all the k-points we want to interpolate 
  WAVE_INTERP%GRID = GRID
  ALLOCATE(WTKPT_INTER(SIZE(VKPT_INTER,2)))
  WTKPT_INTER = 0.0_q
  CALL WDES_SET_KPOINTS(WAVE_INTERP%WDES,VKPT_INTER,WTKPT_INTER)
  DEALLOCATE(WTKPT_INTER)

  MY_USE_ORIGINAL_GRID=.FALSE.; IF (PRESENT(USE_ORIGINAL_GRID)) MY_USE_ORIGINAL_GRID = USE_ORIGINAL_GRID
  IF (MY_USE_ORIGINAL_GRID) THEN
      !In this case the GRID was generated outside this routine and is used here
      WAVE_INTERP%WDES%GRID => WAVE_INTERP%GRID
      WAVE_INTERP%WDES%NCOL = WAVE_INTERP%GRID%RC%NCOL
  ELSE
      CALL GEN_LAYOUT(WAVE_INTERP%GRID, WAVE_INTERP%WDES, LATT_CUR%B, LATT_INI%B, -1, .TRUE.)
#ifdef MPI
      CALL MAPSET(WAVE_INTERP%GRID)
#endif
      WAVE_INTERP%WDES%GRID => WAVE_INTERP%GRID
      WAVE_INTERP%WDES%NCOL = WAVE_INTERP%GRID%RC%NCOL
      !IF (WAVE_INTERP%GRID_ORIG%RC%NCOL /= WAVE_INTERP%GRID%RC%NCOL) THEN
      !   CALL vtutor%bug('The number of columns in the grid changed. original: '//str(WAVE_INTERP%GRID_ORIG%RC%NCOL)//&
      !                   ' new: '//str(WAVE_INTERP%GRID%RC%NCOL),__FILE__,__LINE__)
      !ENDIF
  ENDIF

  ! Update the mode of generating the KPOINTS list to explicit list.
  ! This is important because of the tratment of the singularity in the coulomb kernel in SET_FSG_STORE
  WAVE_INTERP%KPOINTS%MODE = ExplicitList

  IF (LHFCALC) THEN
      CALL SETUP_FULL_KPOINTS(WAVE_INTERP%KPOINTS,LATT_CUR, T_INFO%NIONS, SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, -1,IO%IU0)
  ENDIF

  ! save important dimensions in the data structure 
  WAVE_INTERP%NKPTS_ORIG      = NKPTS_ORIG
  IF (WAVE_INTERP_SET%KPOINTS_OPT_MODE/10==2.AND..NOT.WAVE_INTERP%HAS_WAVEDER) THEN
     CALL vtutor%bug('To use WAVEDER to obtain better estimatives one must first call WAVE_INTERPOLATOR_WAVEDER',__FILE__,__LINE__)
  ENDIF
  WAVE_INTERP%OPTION          = WAVE_INTERP_SET%KPOINTS_OPT_MODE

  ALLOCATE(WAVE_INTERP%KPOINTS_BATCH(3,WAVE_INTERP%MAX_NKPTS_BATCH))
  ALLOCATE(WAVE_INTERP%KPOINTS_INDEX(WAVE_INTERP%MAX_NKPTS_BATCH))
  WAVE_INTERP%NKPTS_BATCH = 0

  CALL KPAR_SYNC_ALL(W%WDES,W)

END SUBROUTINE WAVE_INTERPOLATOR_INIT


!***********************************************************************
!>
!> @brief Interpolate the wavefunctions to a new set of k-points
!> The k-point are passed in the WAVE_INTERP datastructure:
!> WAVE_INTERP%NKPTS_BATCH        - the number of k-points in the batch
!> WAVE_INTERP%KPOINTS_BATCH(3,:) - the coordinates of the k-points
!> and must be set before calling this routine
!>
!***********************************************************************

SUBROUTINE WAVE_INTERPOLATOR_KPOINTS(WAVE_INTERP, HAMILTONIAN, &
                  LATT_CUR, LATT_INI, T_INFO, NONLR_S, NONL_S, LMDIM, P, SV, &
                  CQIJ, CDIJ, SYMM, INFO, IU6, IU0, WDES_COMM)

  USE density_of_states, ONLY: DENINI
  USE pead, ONLY : WDES_FULL_PEAD, PEAD_WA_ROTATE, PEAD_DPSI_DK_1K
  USE fock_ace, ONLY : LFOCK_ACE,FOCK_ACE_ALLOCW,WACE,FOCK_ACE_CONSTRUCT,FOCK_ACE_DEALLOCW
  USE fock, ONLY: LHFCALC, WDES_FOCK, RESETUP_FOCK_WDES, SET_FOCK_KPOINTS
#ifdef fock_dblbuf
  USE fock_dbl, ONLY : FOCK_ALL_DBLBUF
#endif

  ! arguments 
  TYPE (wave_interpolator),TARGET :: WAVE_INTERP
  TYPE (ham_handle)  :: HAMILTONIAN
  TYPE (latt)        :: LATT_CUR
  TYPE (latt)        :: LATT_INI
  TYPE (type_info)   :: T_INFO
  TYPE (nonlr_struct):: NONLR_S
  TYPE (nonl_struct) :: NONL_S
  TYPE (potcar)      :: P(:)
  RGRID           :: SV(:,:)  !< DFT potential (used for the diagonalization of the Hamitlonian)
  OVERLAP         :: CQIJ (:,:,:,:) !< PAW one-center charge
  OVERLAP         :: CDIJ (:,:,:,:) !< PAW strength
  TYPE (symmetry) :: SYMM 
  TYPE (info_struct) :: INFO
  INTEGER IU6, IU0
  INTEGER LMDIM
  TYPE(wavedes),OPTIONAL :: WDES_COMM
  ! local
  TYPE (energy)   :: E
  TYPE (wavespin) :: DW_NEW(3)

  TYPE(wavefuna)  :: WK_FBZ, WK_IBZ, WK, DWK_FBZ(3)
  TYPE(wavedes1)  :: WDESK_FBZ, WDESK_IBZ, WDESK
  TYPE(wavedes)   :: WDES_FBZ
  COMPLEX(q),ALLOCATABLE :: CR(:)

  INTEGER :: NGDIM, NPAR
  INTEGER :: IDIR, IDIR_FBZ
  INTEGER :: IFLAG
  INTEGER :: NELM
  INTEGER :: I,J,K,L,IDX
  INTEGER :: ISP, N, M, IK, IKFULL, IKIBZ, CLOSEST_IKFULL
  INTEGER :: NKPTS_ORIG, IELM, NSIM, NB, ICOUEV
  INTEGER :: NRPLWV_
  REAL(q) :: DIST_CART, MIN_DIST, EXHF, DESUM1, RMS, OLD_TOTEN, ENINI
  INTEGER,ALLOCATABLE :: K_MAP_FROM_IBZ(:,:)
  REAL(q) :: VKPT1(3),VKPT2(3),DVKPT(3),DVKPT_TMP(3)
  REAL(q) :: S(3,3)
  CHARACTER(LEN=5) :: ALGO_STR

  COMPLEX(q),ALLOCATABLE :: CWL(:), CWR(:)
  COMPLEX(q),ALLOCATABLE :: CW_SAVE(:,:)
  GDEF,ALLOCATABLE :: CPROJ_SAVE(:,:)

  NKPTS_ORIG = WAVE_INTERP%NKPTS_ORIG

  ! if WF is allocated we deallocate
  IF (WAVE_INTERP%HAS_W) THEN
    CALL DEALLOCW(WAVE_INTERP%W)
    CALL DEALLOCWDES(WAVE_INTERP%WDES,.TRUE.)
    WAVE_INTERP%HAS_W=.FALSE.
  ENDIF

  ! create a new WF with the original k-points + the new k-points that we want to compute
  CALL WDES_FROM_WDES(WAVE_INTERP%WDES,WAVE_INTERP%W_ORIG%WDES)

  ! set communicators from another WF descriptor
  IF (PRESENT(WDES_COMM)) THEN
      CALL WDES_SET_COMMS(WAVE_INTERP%WDES,WDES_COMM%COMM,WDES_COMM%COMM_KINTER,&
                          WDES_COMM%COMM_KIN,WDES_COMM%COMM_INTER,WDES_COMM%COMM_INB)
  ENDIF

  CALL WDES_SET_NPRO(WAVE_INTERP%WDES,T_INFO,P,INFO%LOVERL)
  
  ! Replace the k-points with original kpoints + the k-points that we want to interpolate to
  DO IK=1,WAVE_INTERP%NKPTS_BATCH
     WAVE_INTERP%KPOINTS%VKPT(:,NKPTS_ORIG+IK)=WAVE_INTERP%KPOINTS_BATCH(:,IK)
  ENDDO
  CALL WDES_SET_KPOINTS(WAVE_INTERP%WDES,WAVE_INTERP%KPOINTS%VKPT,WAVE_INTERP%KPOINTS%WTKPT) 

  ! set first dimension of WF array
  IF (WAVE_INTERP%NGDIM>0) THEN
    CALL WDES_SET_NRPLWV(WAVE_INTERP%WDES,WAVE_INTERP%NGDIM)
  ELSE
    NGDIM = WDES_COUNT_INDPW(WAVE_INTERP%WDES,WAVE_INTERP%GRID,LATT_CUR%B,LATT_INI%B) 
    CALL WDES_SET_NRPLWV(WAVE_INTERP%WDES,NGDIM)
  ENDIF

  ! set number of bands
  IF (WAVE_INTERP%NB_TOT/=WAVE_INTERP%WDES_ORIG%NB_TOT) THEN
    NPAR = WAVE_INTERP%WDES%NB_PAR
    WAVE_INTERP%WDES%NB_TOT = ((WAVE_INTERP%NB_TOT+NPAR-1)/NPAR)*NPAR
    WAVE_INTERP%WDES%NBANDS=WAVE_INTERP%WDES%NB_TOT/WAVE_INTERP%WDES%NB_PAR
  ENDIF
 
  IF (LHFCALC) THEN
    ! now reset NKDIM and NKPTS_FOR_GEN_LAYOUT
    CALL SET_FULL_KPOINTS(WAVE_INTERP%WDES%NKPTS_FOR_GEN_LAYOUT,WAVE_INTERP%WDES%VKPT)
    CALL SET_FOCK_KPOINTS(WAVE_INTERP%WDES%NKDIM)
  ENDIF

  ! Replace the k-points with original kpoints + the k-points that we want to interpolate to
  DO IK=1,WAVE_INTERP%NKPTS_BATCH
     WAVE_INTERP%WDES%VKPT(:,NKPTS_ORIG+IK)=WAVE_INTERP%KPOINTS_BATCH(:,IK)
  ENDDO

  WAVE_INTERP%WDES%NCOL = WAVE_INTERP%GRID%RC%NCOL
  CALL ALLOCWDES(WAVE_INTERP%WDES,.TRUE.)
  WAVE_INTERP%WDES%GRID=>WAVE_INTERP%GRID
  WAVE_INTERP%WDES%NPLWKP=0
  WAVE_INTERP%WDES%LUSEINV=.FALSE.
  WAVE_INTERP%WDES%NINDPW=0
  !CALL GEN_LAYOUT(WAVE_INTERP%GRID, WAVE_INTERP%WDES, LATT_CUR%B, LATT_INI%B, -1, .TRUE.)
  CALL GEN_INDEX(WAVE_INTERP%GRID, WAVE_INTERP%WDES, LATT_CUR%B, LATT_INI%B, -1, -1, .TRUE.) 
  IF (LHFCALC) THEN
    CALL RESETUP_FOCK_WDES( WAVE_INTERP%WDES, LATT_CUR, LATT_INI, IU6)
  ENDIF
 
  ! allocate wf
  CALL ALLOCW(WAVE_INTERP%WDES,WAVE_INTERP%W)
  WAVE_INTERP%HAS_W = .TRUE.

  ! copy wavefunction coeficients and projectors
  DO ISP=1,WAVE_INTERP%ISPIN
     DO IK=1,NKPTS_ORIG
        if (SIZE(WAVE_INTERP%W_ORIG%WDES%NINDPW,1)>SIZE(WAVE_INTERP%WDES%NINDPW,1)) then
            call  vtutor%bug('dimensions are incompatible.'//&
                             ' original: '//str(SIZE(WAVE_INTERP%W_ORIG%WDES%NINDPW,1))//&
                             ' new: '//str(SIZE(WAVE_INTERP%WDES%NINDPW,1)),__FILE__,__LINE__)
        endif
        if (WAVE_INTERP%W_ORIG%WDES%NGVECTOR(IK)/=WAVE_INTERP%WDES%NGVECTOR(IK)) then
            call  vtutor%bug('',__FILE__,__LINE__)
        endif
        if (any(WAVE_INTERP%W_ORIG%WDES%NINDPW(:WAVE_INTERP%W_ORIG%WDES%NGVECTOR(IK),IK)/=&
                       WAVE_INTERP%WDES%NINDPW(:WAVE_INTERP%WDES%NGVECTOR(IK),IK))) then
            call  vtutor%bug('',__FILE__,__LINE__)
        endif

        DO N=1,WAVE_INTERP%NB_TOT
           WAVE_INTERP%W%CELTOT(N,IK,ISP) = WAVE_INTERP%W_ORIG%CELTOT(N,IK,ISP)
           WAVE_INTERP%W%FERTOT(N,IK,ISP) = WAVE_INTERP%W_ORIG%FERTOT(N,IK,ISP)
           WAVE_INTERP%W%AUXTOT(N,IK,ISP) = WAVE_INTERP%W_ORIG%AUXTOT(N,IK,ISP)
        ENDDO
     ENDDO
  ENDDO

  ! For now NCORE must be 1
  IF (WAVE_INTERP%WDES%COMM_INB%NCPU/=1.AND.NKPTS_ORIG>0) THEN
     call vtutor%bug('Copying orbitals from the original set is not yet implemented for NCORE/=1',__FILE__,__LINE__)
  ENDIF

  ! I check if the number of local bands is the same between W_ORIG AND W
  !IF (WAVE_INTERP%WDES_ORIG%NBANDS /= WAVE_INTERP%WDES%NBANDS) THEN
  !   call vtutor%bug('Number of bands in the original wavefunction '//str(WAVE_INTERP%WDES_ORIG%NBANDS)//&
  !                   ' is different from the one where we want to interpolate to '//str(WAVE_INTERP%WDES%NBANDS),__FILE__,__LINE__)
  !ENDIF

  NRPLWV_ = WAVE_INTERP%W_ORIG%WDES%NRPLWV
  DO ISP=1,WAVE_INTERP%ISPIN
     DO IK=1,NKPTS_ORIG
        DO N=1,MIN(WAVE_INTERP%WDES_ORIG%NBANDS,WAVE_INTERP%WDES%NBANDS)
           WAVE_INTERP%W%CW(:NRPLWV_,N,IK,ISP)   = WAVE_INTERP%W_ORIG%CW(:NRPLWV_,N,IK,ISP)
           WAVE_INTERP%W%CPROJ(:,N,IK,ISP)= WAVE_INTERP%W_ORIG%CPROJ(:,N,IK,ISP)
        ENDDO
     ENDDO
  ENDDO

  IF (NONL_S%LRECIP) THEN
     CALL NONL_DEALLOC(NONL_S)
     CALL NONL_ALLOC (NONL_S,T_INFO,P,WAVE_INTERP%W%WDES, .FALSE.)

     CALL SPHER(WAVE_INTERP%GRID, NONL_S, P, WAVE_INTERP%W%WDES, LATT_CUR, 1)
     CALL PHASE(WAVE_INTERP%W%WDES,NONL_S,0)
  ENDIF

  ! Initilize the orbitals
  SELECT CASE (WAVE_INTERP%OPTION/10)
     ! Use random wavefunctions
     CASE(0)
        DO ISP=1,WAVE_INTERP%ISPIN
            CALL DENINI(WAVE_INTERP%W%FERTOT(:,:,ISP),WAVE_INTERP%WDES%NB_TOT,WAVE_INTERP%WDES%NKPTS,INFO%NELECT,WAVE_INTERP%WDES%LNONCOLLINEAR, NKSTART=NKPTS_ORIG+1)
        ENDDO
        ! initialize with random numbers
        ENINI = INFO%ENINI !INFO%ENMAX
        CALL WFINIT(WAVE_INTERP%WDES, WAVE_INTERP%W, ENINI, NKSTART=NKPTS_ORIG+1)

     ! Use a closest WF as intial guess
     CASE(1,2)
        CALL vtutor%error('KPOINTS_OPT_MODE '//str(WAVE_INTERP%OPTION)//' is not implemented')
#if 0
        ! To initialize the new wavefunctions with a good initial guess
        ! we find the closest k-point from the coarse grid in the FBZ
        ! Loop over all the points to be interpolated
        ALLOCATE(K_MAP_FROM_IBZ(2,WAVE_INTERP%NKPTS_BATCH))
        DO IK=1,WAVE_INTERP%NKPTS_BATCH
           VKPT1 = WAVE_INTERP%KPOINTS_BATCH(:,IK)
           ! Find closest point in full grid using naive algorithm
           ! TODO: fast nearest-neighbour search
           MIN_DIST = 1e8
           CLOSEST_IKFULL = -1
           DO IKFULL=1,KPOINTS_FULL_ORIG%NKPTS
              VKPT2 = KPOINTS_FULL_ORIG%VKPT(:,IKFULL)
              DVKPT = VKPT1-VKPT2
              CALL DIRKAR(1, DVKPT(1), LATT_CUR%B)
              DIST_CART = DVKPT(1)**2+DVKPT(2)**2+DVKPT(3)**2
              IF (DIST_CART>MIN_DIST) CYCLE
              MIN_DIST = DIST_CART
              CLOSEST_IKFULL = IKFULL
           END DO
           K_MAP_FROM_IBZ(1,IK) = CLOSEST_IKFULL
           K_MAP_FROM_IBZ(2,IK) = KPOINTS_FULL_ORIG%NEQUIV(CLOSEST_IKFULL)
        END DO

        ! Replace the k-points with k-points that we want to interpolate to
        DO IK=1,WAVE_INTERP%NKPTS_BATCH
           W%WDES%VKPT(:,NKPTS_ORIG+IK)=WAVE_INTERP%KPOINTS_BATCH(:,IK)
        ENDDO
        CALL KPAR_SYNC_ALL(W%WDES,W)
        CALL RE_GEN_LAYOUT(GRID, W%WDES, WAVE_INTERP%KPOINTS, LATT_CUR, LATT_INI, IU6, IU0, NGDIM=WAVE_INTERP%NGDIM)
        CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
        IF (WAVE_INTERP%OPTION/10==2) THEN
           CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(1), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
           CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(2), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
           CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(3), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
        ENDIF

        ! to map from the adjacent k-point in the FBZ to the current point
        ! we use Fourier tranformation in real space
        ALLOCATE(CR(W%WDES%GRID%RL%NP))
        W_NEW=W
        W_NEW%WDES=>WDES_FULL_PEAD
        IF (WAVE_INTERP%OPTION/10==2) THEN
           ALLOCATE(CWR(GRID%NPLWV))
           ALLOCATE(CWL(W%WDES%NRPLWV))
           ALLOCATE(CW_SAVE(SIZE(W_NEW%CW,1),SIZE(W_NEW%CW,2)))
           ALLOCATE(CPROJ_SAVE(SIZE(W_NEW%CPROJ,1),SIZE(W_NEW%CPROJ,2)))
           DO IDIR=1,3
              DW_NEW(IDIR)=WAVE_INTERP%WDIFF3(IDIR)
              DW_NEW(IDIR)%WDES=>WDES_FULL_PEAD
           ENDDO
        ENDIF
        DO ISP=1,W%WDES%ISPIN
           DO IK=1,WAVE_INTERP%NKPTS_BATCH
              IKFULL = K_MAP_FROM_IBZ(1,IK)
              IKIBZ  = K_MAP_FROM_IBZ(2,IK)
              ! Create WDESK
              CALL SETWDES(W%WDES,WDESK,NKPTS_ORIG+IK)
              ! Create WDESK_FBZ
              CALL SETWDES(W_NEW%WDES,WDESK_FBZ,IKFULL)
              WDESK_FBZ%NRPLWV = W%WDES%NRPLWV
              IF (IKFULL==IKIBZ) THEN
                 ! k is a kpoint in the IBZ
                 WK_FBZ=ELEMENTS(W_NEW,WDESK_FBZ,ISP)
              ELSE
                 ! k is not a kpoint in the IBZ
                 CALL NEWWAVA(WK_FBZ,WDESK_FBZ,WDESK_FBZ%NBANDS)
                 CALL PEAD_WA_ROTATE(W_NEW,P,LATT_CUR,ISP,WK_FBZ)
              ENDIF
              IF (WAVE_INTERP%OPTION/10==2) THEN
                 ! Rotate all directions from IBZ to FBZ
                 DO IDIR=1,3
                    IF (IKFULL==IKIBZ) THEN
                       ! k is a kpoint in the IBZ
                       DWK_FBZ(IDIR)=ELEMENTS(DW_NEW(IDIR),WDESK_FBZ,ISP)
                    ELSE
                       ! k is not a kpoint in the IBZ
                       CALL NEWWAVA(DWK_FBZ(IDIR),WDESK_FBZ,WDESK_FBZ%NBANDS)
                       CALL PEAD_WA_ROTATE(DW_NEW(IDIR),P,LATT_CUR,ISP,DWK_FBZ(IDIR))
                    ENDIF
                 ENDDO
                 ! determine transformation matrix in real space when going from IBZ to FBZ k-point
                 S=0
                 DO L=1,3
                    DO K=1,3
                       DO J=1,3
                          DO I=1,3
                             S(L,I)=S(L,I)+LATT_CUR%A(L,K)*&
                                    KPOINTS_FULL%ISYMOP(J,K,IKFULL)*LATT_CUR%B(I,J)
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
                 ! k-point difference
                 DVKPT = WDESK%VKPT - WDESK_FBZ%VKPT
                 CALL KARDIR(1,DVKPT,LATT_CUR%B)
                 ! transform the k-point difference
                 DVKPT_TMP = DVKPT
                 DVKPT = MATMUL(S,DVKPT_TMP)
              ENDIF
              ! for all the bands
              DO NB=1,W%WDES%NBANDS
                 IF (WAVE_INTERP%OPTION/10==2) THEN
                    !CWL = WK_FBZ%CW(:,NB)
                    CWL = WK_FBZ%CW(:,NB) + CITPI*( &
                                            DVKPT(1)*DWK_FBZ(1)%CW(:,NB) + &
                                            DVKPT(2)*DWK_FBZ(2)%CW(:,NB) + &
                                            DVKPT(3)*DWK_FBZ(3)%CW(:,NB) )
                    ! rec (fbz point) -> real
                    !CALL FFTWAV(WDESK_FBZ,CR,WK_FBZ%CW(:,NB))
                    CALL FFTWAV(WDESK_FBZ%NGVECTOR,WDESK_FBZ%NINDPW,CR,CWL,WDESK_FBZ%GRID)
                    ! real -> rec (point that we want)
                    CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW,CR,W%CW(:,NB,NKPTS_ORIG+IK,ISP),WDESK%GRID,.FALSE.)
                    ! For the projectors we just copy from the FBZ
                    !W%CPROJ(:,NB,NKPTS_ORIG+IK,ISP) = WK_FBZ%CPROJ(:,NB)
                    W%CPROJ(:,NB,NKPTS_ORIG+IK,ISP) = WK_FBZ%CPROJ(:,NB) + CITPI*( &
                                                      DVKPT(1)*DWK_FBZ(1)%CPROJ(:,NB) + &
                                                      DVKPT(2)*DWK_FBZ(2)%CPROJ(:,NB) + &
                                                      DVKPT(3)*DWK_FBZ(3)%CPROJ(:,NB) )
                    W%CELTOT(NB,NKPTS_ORIG+IK,ISP) = W%CELTOT(NB,IKIBZ,ISP) + CITPI*( &
                                                     DVKPT(1)*DW_NEW(1)%CELTOT(NB,IKIBZ,ISP) + &
                                                     DVKPT(2)*DW_NEW(2)%CELTOT(NB,IKIBZ,ISP) + &
                                                     DVKPT(3)*DW_NEW(3)%CELTOT(NB,IKIBZ,ISP) )
                 ELSE
                    ! rec (fbz point) -> real
                    CALL FFTWAV(WDESK_FBZ%NGVECTOR,WDESK_FBZ%NINDPW,CR,WK_FBZ%CW(:,NB),WDESK%GRID)
                    ! real -> rec (point that we want)
                    CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW,CR,W%CW(:,NB,NKPTS_ORIG+IK,ISP),WDESK%GRID,.FALSE.)
                    ! For the projectors we just copy from the FBZ
                    W%CPROJ(:,NB,NKPTS_ORIG+IK,ISP)= WK_FBZ%CPROJ(:,NB)
                    W%CELTOT(NB,NKPTS_ORIG+IK,ISP) = W%CELTOT(NB,IKIBZ,ISP)
                 ENDIF
                 ! Copy the rest of the information
                 W%FERTOT(NB,NKPTS_ORIG+IK,ISP) = W%FERTOT(NB,IKIBZ,ISP)
                 W%AUXTOT(NB,NKPTS_ORIG+IK,ISP) = W%AUXTOT(NB,IKIBZ,ISP)
              ENDDO
              IF (IK/=KPOINTS_FULL%NEQUIV(IK)) THEN
                 CALL DELWAVA(WK_FBZ)
                 IF (WAVE_INTERP%OPTION/10==2) THEN
                    CALL DELWAVA(DWK_FBZ(1))
                    CALL DELWAVA(DWK_FBZ(2))
                    CALL DELWAVA(DWK_FBZ(3))
                 ENDIF
              ENDIF
           ENDDO ! kpoints batch
        ENDDO ! spin index
        DEALLOCATE(CR)
        DEALLOCATE(K_MAP_FROM_IBZ)

        DO ISP=1,W%WDES%ISPIN
            CALL DENINI(W%FERTOT(:,:,ISP),W%WDES%NB_TOT,W%WDES%NKPTS,INFO%NELECT,W%WDES%LNONCOLLINEAR, NKSTART=NKPTS_ORIG+1)
        ENDDO

        ! initialize with random numbers beyond batch
        ENINI = INFO%ENINI !INFO%ENMAX
        CALL WFINIT(W%WDES, W, ENINI, NKSTART=NKPTS_ORIG+WAVE_INTERP%NKPTS_BATCH+1)
#endif
  END SELECT

  ! Prepare call to diag by orthonormalizing the orbitals
  CALL SET_DATAKE(WAVE_INTERP%WDES, LATT_CUR%B)
  IF (NONLR_S%LREAL) THEN
     CALL RSPHER(WAVE_INTERP%GRID,NONLR_S,LATT_CUR)
  ELSE
     CALL SPHER(WAVE_INTERP%GRID,NONL_S,P,WAVE_INTERP%WDES,LATT_CUR, 1)
  ENDIF
  CALL PROALL(WAVE_INTERP%GRID, LATT_CUR, NONLR_S, NONL_S, WAVE_INTERP%W)
  CALL ORTHCH(WAVE_INTERP%WDES, WAVE_INTERP%W, WAVE_INTERP%WDES%LOVERL, LMDIM, CQIJ, NKSTART=NKPTS_ORIG+1)

  IF (LFOCK_ACE()) THEN
    CALL FOCK_ACE_ALLOCW(WAVE_INTERP%WDES)
  ENDIF

  ! Call diagonalization
  NSIM=WAVE_INTERP%WDES%NSIM*2
#ifdef MPI
  NSIM=((WAVE_INTERP%WDES%NSIM*2+WAVE_INTERP%WDES%COMM_INTER%NCPU-1)/ &
         WAVE_INTERP%WDES%COMM_INTER%NCPU)*WAVE_INTERP%WDES%COMM_INTER%NCPU
#endif

  OLD_TOTEN =BANDSTRUCTURE_ENERGY_NOWTKPT(WAVE_INTERP%WDES, WAVE_INTERP%W, NKPTS_ORIG+1)
  IF (IU0>=0) THEN
    WRITE(17,*)  '      N       E                     dE             ncg'
    WRITE(IU0,*) '      N       E                     dE             ncg'
  END IF

  NELM=INFO%NELM; IF (MOD(WAVE_INTERP%OPTION,10)==2) NELM=1
  DO IELM=1,NELM

#ifdef fock_dblbuf
     IF (LFOCK_ACE()) THEN
       CALL FOCK_ALL_DBLBUF(WAVE_INTERP%WDES,WAVE_INTERP%W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
                            NKSTART=NKPTS_ORIG+1,EX=E%EXHF,EX_ACFDT=E%EXHF_ACFDT,XI=WACE)
  !    CALL FOCK_ACE_ACC_PAW(W%WDES,W,LATT_CUR,NONL_S,NONLR_S)
       CALL FOCK_ACE_CONSTRUCT(WAVE_INTERP%WDES,WAVE_INTERP%W,NKSTART=NKPTS_ORIG+1)
     ENDIF
#endif

     SELECT CASE (MOD(WAVE_INTERP%OPTION,10))
        CASE (0)
          ! Lazy mode: do nothing!
        CASE (1)
             CALL EDDAV(HAMILTONIAN, P, WAVE_INTERP%GRID, INFO, LATT_CUR, NONLR_S, NONL_S, &
                       WAVE_INTERP%W, WAVE_INTERP%WDES, NSIM, LMDIM, CDIJ, CQIJ, RMS, DESUM1, ICOUEV, SV, &
                       E%EXHF, IU6,IU0, LDELAY=.FALSE., LSUBROTI=INFO%LDIAG, LHF=.TRUE.,&
                       LEMPTY=.FALSE.,NKSTART=NKPTS_ORIG+1,NKSTOP=NKPTS_ORIG+WAVE_INTERP%NKPTS_BATCH)
             ALGO_STR='DAV: '
        CASE (2)
             IFLAG = 3 ! exact diagonalization
             CALL EDDIAG_EXACT(HAMILTONIAN,WAVE_INTERP%GRID,LATT_CUR,NONLR_S,NONL_S,WAVE_INTERP%W,WAVE_INTERP%WDES,SYMM, &
                     LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,IU6,E%EXHF,E%EXHF_ACFDT,&
                     NKSTART=NKPTS_ORIG+1)
             ALGO_STR='DIA: '
     END SELECT

     E%EBANDSTR=BANDSTRUCTURE_ENERGY_NOWTKPT(WAVE_INTERP%WDES, WAVE_INTERP%W, NKPTS_ORIG+1)

     IF (IU0>=0) THEN
        WRITE(17,  '(A,I3,A,E20.12,A,E12.5,A,I6)') ALGO_STR,IELM,'   ',E%EBANDSTR,'   ',E%EBANDSTR-OLD_TOTEN,'   ',ICOUEV
        WRITE(IU0, '(A,I3,A,E20.12,A,E12.5,A,I6)') ALGO_STR,IELM,'   ',E%EBANDSTR,'   ',E%EBANDSTR-OLD_TOTEN,'   ',ICOUEV
     END IF

     IF (ABS(OLD_TOTEN-E%EBANDSTR) < ABS(INFO%EDIFF) ) EXIT
     OLD_TOTEN = E%EBANDSTR

  ENDDO

  IF (LFOCK_ACE()) CALL FOCK_ACE_DEALLOCW

  CONTAINS
  !> Compute the eigenvalues contribution to the total energy associated to the band-structure
  !>
  !> The difference in this routine w.r.t the hamil.F::BANDSTRUCTURE_ENERGY is that
  !> we don't use the weights of the k-points
  !> to compute the change of the total energy. This is necessary as the correct
  !> weights of the original points should be mantained for Hybrid calculations
  !> @todo Check whether to remove also FETOT as we might want to also converge empty states
  FUNCTION BANDSTRUCTURE_ENERGY_NOWTKPT(WDES, W, NKSTART) RESULT(E)
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W      ! wavefunction
    REAL(q)            E
    INTEGER            NKSTART
    INTEGER            ISP, NK, NB

    E=0
    DO ISP=1,WDES%ISPIN
       DO NK=NKSTART,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          IF (WAVE_INTERP%FULL_ACC) THEN
             DO NB=1,WDES%NB_TOT
                E=E+WDES%RSPIN* REAL( W%CELTOT(NB,NK,ISP) ,KIND=q) /WDES%ISPIN
             ENDDO
          ELSE
             DO NB=1,WDES%NB_TOT
                E=E+WDES%RSPIN* REAL( W%CELTOT(NB,NK,ISP) ,KIND=q) *W%FERTOT(NB,NK,ISP)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    CALLMPI( M_sum_d(WDES%COMM_KINTER, E, 1))
  END FUNCTION BANDSTRUCTURE_ENERGY_NOWTKPT

END SUBROUTINE WAVE_INTERPOLATOR_KPOINTS


!***********************************************************************
!>
!> @brief Restore wavefunctions to previous state and free interpolator structure
!>
!***********************************************************************

SUBROUTINE WAVE_INTERPOLATOR_FREE(WAVE_INTERP, KPOINTS, GRID, &
                  LATT_CUR, LATT_INI, T_INFO, NONL_S, W, P, SYMM) 

  use fock, only: LHFCALC, RESETUP_FOCK_WDES, SET_FOCK_KPOINTS
  TYPE (wave_interpolator) WAVE_INTERP
  TYPE (kpoints_struct) KPOINTS
  TYPE (grid_3d)     GRID
  TYPE (latt)        LATT_CUR
  TYPE (latt)        LATT_INI
  TYPE (type_info)   T_INFO
  TYPE (nonl_struct) NONL_S
  TYPE (potcar)      P(:)
  TYPE (wavespin), TARGET :: W
  TYPE (symmetry) SYMM
  ! local varaibles
  INTEGER :: IDIR

  CALL DEALLOCWDES(WAVE_INTERP%WDES,.TRUE.)
  DEALLOCATE(WAVE_INTERP%KPOINTS_BATCH)
  DEALLOCATE(WAVE_INTERP%KPOINTS_INDEX)
  IF (WAVE_INTERP%OPTION/10==2) THEN
     DO IDIR=1,3
         CALL DEALLOCW(WAVE_INTERP%WDIFF3(IDIR))
     ENDDO
  ENDIF

  IF (LHFCALC) THEN
    CALL RE_GEN_LAYOUT(GRID, W%WDES, KPOINTS, LATT_CUR, LATT_INI, -1,-1)
  ENDIF

  IF (NONL_S%LRECIP) THEN
     CALL NONL_DEALLOC(NONL_S)
     CALL NONL_ALLOC (NONL_S,T_INFO,P,W%WDES, .FALSE.)

     CALL SPHER(GRID, NONL_S, P, W%WDES, LATT_CUR, 1)
     CALL PHASE(W%WDES,NONL_S,0)
  ENDIF

END SUBROUTINE WAVE_INTERPOLATOR_FREE

!> Initialize wave mapper datastructure
SUBROUTINE WAVE_MAPPER_INIT(SELF,NKPTS_ORIG,NKPTS_BATCH,KPOINTS_BATCH,KPOINTS_INDEX)
  TYPE (wave_mapper) :: SELF
  INTEGER :: NKPTS_ORIG
  INTEGER :: NKPTS_BATCH
  real(q) :: KPOINTS_BATCH(:,:)
  INTEGER :: KPOINTS_INDEX(:)
  SELF%NKPTS_ORIG = NKPTS_ORIG
  SELF%NKPTS_BATCH = NKPTS_BATCH
  ALLOCATE(SELF%KPOINTS_BATCH(3,NKPTS_BATCH))
  SELF%KPOINTS_BATCH(:,:NKPTS_BATCH) = KPOINTS_BATCH(:,:NKPTS_BATCH)
  ALLOCATE(SELF%KPOINTS_INDEX(NKPTS_BATCH))
  SELF%KPOINTS_INDEX(:NKPTS_BATCH) = KPOINTS_INDEX(:NKPTS_BATCH)
END SUBROUTINE WAVE_MAPPER_INIT

!> Create a wave mapper datastructure from a wave_interp object
SUBROUTINE WAVE_MAPPER_FROM_WAVE_INTERP(SELF,WAVE_INTERP)
  TYPE (wave_mapper) :: SELF
  TYPE (wave_interpolator) :: WAVE_INTERP
  CALL WAVE_MAPPER_INIT(SELF,WAVE_INTERP%NKPTS_ORIG,WAVE_INTERP%NKPTS_BATCH,WAVE_INTERP%KPOINTS_BATCH,WAVE_INTERP%KPOINTS_INDEX)
END SUBROUTINE WAVE_MAPPER_FROM_WAVE_INTERP

!> Free wave mapper datastructure
SUBROUTINE WAVE_MAPPER_FREE(SELF)
  TYPE (wave_mapper) :: SELF
  DEALLOCATE(SELF%KPOINTS_BATCH)
  DEALLOCATE(SELF%KPOINTS_INDEX)
END SUBROUTINE WAVE_MAPPER_FREE

END MODULE wave_interpolate
