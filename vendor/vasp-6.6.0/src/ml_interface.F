#include "symbol.inc"

!**********************************************************************
!
! this module is interfacing VASP to a machine learning code
! it is part of VASP and coverged by the general VASP 
! copy right
!
!**********************************************************************
MODULE ml_interface
  USE mpimy
  USE prec
  USE tutor, ONLY: vtutor, isError, RandomSeed, argument
#ifdef ML_AVAILABLE
  USE ml_main_subroutines
#endif
  IMPLICIT NONE
! this variable specifies whether the MLFF calculation is executed or not.
     LOGICAL :: ML_LMLFF
! this variable determines whether ab initio calculation is executed.
     LOGICAL :: LDO_AB_INITIO=.TRUE.
     LOGICAL :: LRUN_VASP_TO_ML=.FALSE.
! this variable specifies the results of the judgement of the necessity of ab initio calculations on new configurations.
! if ML_LMLFF=.FALSE., LDO_AB_INITIO_NEW=LDO_AB_INITIO=.TRUE.
     LOGICAL :: LDO_AB_INITIO_NEW=.TRUE.
! this variable determines whether the calculations are executed using only the
! force field or not.
     LOGICAL :: ML_FF_LMLONLY = .FALSE.

  CONTAINS

!=======================================================================
!
! This subroutine reads variables necessary for the machine-learning force field
! generations from INCAR file.
!
!=======================================================================

   SUBROUTINE MACHINE_LEARNING_GENERAL_AND_INSTANCE_READER(INCAR)

      USE incar_reader, only: INCAR_FILE, PROCESS_INCAR, COUNT_ELEMENTS
#ifdef ML_AVAILABLE
      ! Warning: "tag" states are not set in this routine because the
      ! corresponding instance "ML_SUPER_TYPE%TAG_LIST" is not available here.
      USE ML_FF_TAGLIST, ONLY: TAGLIST_GET_ALT_TAGS
      USE ML_FF_STRUCT, ONLY: ML_FILESTREAM_GLOBAL, ML_TOTNUM_INSTANCES
#endif
      IMPLICIT NONE
      TYPE(INCAR_FILE), INTENT(INOUT) :: INCAR
      ! Local variables
      CHARACTER(len=:), ALLOCATABLE :: ALT_TAGS(:)
      INTEGER                       :: I
      ! Set ML_LMLFF as .FALSE.
      ML_LMLFF=.FALSE.
      ! Read ML_LMLFF.
      !   ML_LMLFF = .FALSE.: No machine-learning force field generations
      !   ML_LMLFF = .TRUE. : Execute machine-learning force field
      !   generations
      ! Need to set alternatives manually because the taglist module is
      ! not available here.
      ALT_TAGS = ['ML_LMLFF   ', 'ML_FF_LMLFF']
      DO I = 1, SIZE(ALT_TAGS)
         IF (COUNT_ELEMENTS(INCAR, TRIM(ALT_TAGS(I))) > 0) THEN
            CALL PROCESS_INCAR(INCAR, TRIM(ALT_TAGS(I)), ML_LMLFF)
            ! Skip alternative tags if one is already found.
            EXIT
         END IF
      END DO
#ifdef ML_AVAILABLE
      ! By default the number of instances is 1
      ML_TOTNUM_INSTANCES = 1
      ! At this point there is no need for multiple instances so we disable
      ! reading of this tag.
      !ALT_TAGS = TAGLIST_GET_ALT_TAGS('ML_TOTNUM_INSTANCES')
      !DO I = 1, SIZE(ALT_TAGS)
      !   IF (COUNT_ELEMENTS(INCAR, TRIM(ALT_TAGS(I))) > 0) THEN
      !      CALL PROCESS_INCAR(INCAR, TRIM(ALT_TAGS(I)), ML_TOTNUM_INSTANCES)
      !      ! Skip alternative tags if one is already found.
      !      EXIT
      !   END IF
      !END DO

      ! Global filestream variable
      ML_FILESTREAM_GLOBAL=500
      ALT_TAGS = TAGLIST_GET_ALT_TAGS('ML_FILESTREAM_START')
      DO I = 1, SIZE(ALT_TAGS)
         IF (COUNT_ELEMENTS(INCAR, TRIM(ALT_TAGS(I))) > 0) THEN
            CALL PROCESS_INCAR(INCAR, TRIM(ALT_TAGS(I)), ML_FILESTREAM_GLOBAL)
            ! Skip alternative tags if one is already found.
            EXIT
         END IF
      END DO
#else
      IF (ML_LMLFF) THEN
         CALL vtutor%error("This VASP executable does not support machine &
                           &learning (ML_LMLFF=.TRUE.)! Try recompiling with &
                           &-DML_AVAILABLE")
      ENDIF
#endif
   END SUBROUTINE MACHINE_LEARNING_GENERAL_AND_INSTANCE_READER

#ifdef ML_AVAILABLE
#if defined(MPI) || defined(MPI_CHAIN)
!=======================================================================
!
! Routines important for the handling of global machine learning handles.
!
!=======================================================================

     SUBROUTINE MACHINE_LEARNING_HANDLE_ALLOCATE(NUM_INSTANCES)
       USE ML_FF_STRUCT, ONLY: ML_SUPER_HANDLE_MAIN
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NUM_INSTANCES
       IF (ALLOCATED(ML_SUPER_HANDLE_MAIN)) THEN
          DEALLOCATE(ML_SUPER_HANDLE_MAIN)
       ENDIF
       ALLOCATE(ML_SUPER_HANDLE_MAIN(NUM_INSTANCES))
     END SUBROUTINE MACHINE_LEARNING_HANDLE_ALLOCATE

     SUBROUTINE MACHINE_LEARNING_HANDLE_DEALLOCATE
       USE ML_FF_STRUCT, ONLY: ML_SUPER_HANDLE_MAIN
       IMPLICIT NONE
       IF (ALLOCATED(ML_SUPER_HANDLE_MAIN)) THEN
          DEALLOCATE(ML_SUPER_HANDLE_MAIN)
       ENDIF
     END SUBROUTINE MACHINE_LEARNING_HANDLE_DEALLOCATE

!=======================================================================
!
! the following routine is called once, before initialization is called
! to set a global COMM_WORLD communicator. This communicator should be
! only used very sparsely for emergency aborts and etc.
!
!=======================================================================
    SUBROUTINE MACHINE_LEARNING_SET_COMM_WORLD_GLOBAL(VASP_COMMUNICATOR)
      USE mpimy, ONLY : communic
      USE ML_FF_STRUCT, ONLY : ML_MPI_PAR, ML_COMM_WORLD_GLOBAL, ML_IO_WRITE
      IMPLICIT NONE
      TYPE (communic)   :: VASP_COMMUNICATOR
      ML_COMM_WORLD_GLOBAL%MPI_COMM = VASP_COMMUNICATOR%MPI_COMM
      ML_COMM_WORLD_GLOBAL%NCPU     = VASP_COMMUNICATOR%NCPU
      ML_COMM_WORLD_GLOBAL%NODE_ME     = VASP_COMMUNICATOR%NODE_ME
      ! set here which rank is allowed to write (currently only the first rank)
      IF (ML_COMM_WORLD_GLOBAL%NODE_ME.EQ.1) THEN
         ML_IO_WRITE = .TRUE.
      ELSE
         ML_IO_WRITE = .FALSE.
      ENDIF
    END SUBROUTINE MACHINE_LEARNING_SET_COMM_WORLD_GLOBAL

#ifndef ML_LEARNING_NOT_ACTIVE
    SUBROUTINE INQUIRE_ML_AB_FILE(DIR_APP, DIR_LEN, ML_SUPER_HANDLE)
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE IOHANDLE
      IMPLICIT NONE
      CHARACTER(LEN=10), INTENT(IN) ::  DIR_APP
      INTEGER, INTENT(IN)           ::  DIR_LEN
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      CHARACTER(LEN=4)     :: FILE_POSTFIX
      ML_SUPER_HANDLE%FF%LMLABEXIST=.FALSE.
      CALL ML_SET_POSTFIX(ML_SUPER_HANDLE%INSTANCE,FILE_POSTFIX) 
      IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
         INQUIRE( FILE=DIR_APP(1:DIR_LEN) // 'ML_AB', EXIST=ML_SUPER_HANDLE%FF%LMLABEXIST )
      ELSE
         INQUIRE( FILE=DIR_APP(1:DIR_LEN) // 'ML_AB' // FILE_POSTFIX, EXIST=ML_SUPER_HANDLE%FF%LMLABEXIST )
      ENDIF
    END SUBROUTINE INQUIRE_ML_AB_FILE
#endif

    SUBROUTINE INQUIRE_ML_FF_FILE(DIR_APP, DIR_LEN, ML_SUPER_HANDLE)
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE IOHANDLE
      IMPLICIT NONE
      CHARACTER(LEN=10), INTENT(IN) ::  DIR_APP
      INTEGER, INTENT(IN)           ::  DIR_LEN
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      CHARACTER(LEN=4)     :: FILE_POSTFIX
      ML_SUPER_HANDLE%FF%LMLFFEXIST=.FALSE.
      CALL ML_SET_POSTFIX(ML_SUPER_HANDLE%INSTANCE,FILE_POSTFIX) 
      IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
         INQUIRE( FILE=DIR_APP(1:DIR_LEN) // 'ML_FF', EXIST=ML_SUPER_HANDLE%FF%LMLFFEXIST )
      ELSE
         INQUIRE( FILE=DIR_APP(1:DIR_LEN) // 'ML_FF' // FILE_POSTFIX, EXIST=ML_SUPER_HANDLE%FF%LMLFFEXIST )
      ENDIF
    END SUBROUTINE INQUIRE_ML_FF_FILE

!======================================================================
!> Handle modifications resulting from VCAIMAGES synchronization.
!>
!> In case of VCAIMAGES some output-related variables may have changed.
!> This subroutine propagates these changes back inside the ML
!> structures and handles any additional clean-up necessary.
!=======================================================================
    SUBROUTINE MACHINE_LEARNING_POST_SYNC_VCAIMAGES(ML_SUPER_HANDLE, OUTBLOCK_MAIN)
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE ML_FF_TAGLIST
      USE ml_ff_string, only: str
      IMPLICIT NONE
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      INTEGER, INTENT(IN)  :: OUTBLOCK_MAIN

      ! Early exit if no change.
      IF (OUTBLOCK_MAIN == ML_SUPER_HANDLE%FF%OUTBLOCK) RETURN

      ! If ML_OUTBLOCK was set in INCAR, set the tag state to "override".
      IF (ML_SUPER_HANDLE%TAG_LIST%GET_STATE("ML_OUTBLOCK") == TL_INCAR .OR. &
          ML_SUPER_HANDLE%TAG_LIST%GET_STATE("ML_OUTBLOCK") == TL_INCAR_ALT ) THEN
         CALL vtutor%alert("The ML_OUTBLOCK value from the INCAR file (" // &
            str(ML_SUPER_HANDLE%FF%OUTBLOCK) // ") was automatically reset to the minimum value (" &
            // str(OUTBLOCK_MAIN) // ") over all images used in this run (VCAIMAGES).")
         CALL ML_SUPER_HANDLE%TAG_LIST%SET_STATE("ML_OUTBLOCK", TL_OVERRIDE)
      ENDIF
      ! Update internal ML_OUTBLOCK.
      ML_SUPER_HANDLE%FF%OUTBLOCK = OUTBLOCK_MAIN

      ! Exit if ML_ESTBLOCK is zero, no need to change, whatever happened to ML_OUTBLOCK.
      IF (ML_SUPER_HANDLE%FF%IERR == 0) RETURN

      ! ML_OUTBLOCK cannot be set individually in each window. Hence, it could be either
      ! (1) == 1 in any mixed DFT, ML setup
      ! (2) >= 1 in a pure ML setup
      ! If the user has set a value for ML_ESTBLOCK in an image, it is still compatible here in both
      ! cases because compatibility was already checked in the ML reader. Hence, only in case (1)
      ! where one of the images set ML_ESTBLOCK to ML_OUTBLOCK by default, the value needs to be
      ! updated here.
      IF (ML_SUPER_HANDLE%TAG_LIST%GET_STATE("ML_ESTBLOCK") /= TL_INCAR .AND. &
          ML_SUPER_HANDLE%TAG_LIST%GET_STATE("ML_ESTBLOCK") /= TL_INCAR_ALT ) THEN
         ML_SUPER_HANDLE%FF%IERR = OUTBLOCK_MAIN
      ENDIF

    END SUBROUTINE MACHINE_LEARNING_POST_SYNC_VCAIMAGES

!=======================================================================
!
! the following routine is called once, after VASP has read
! the ionic positions + the Bravais lattice
! at this point the machine learning code is supposed to perform
! the basic setup
!
!=======================================================================

    SUBROUTINE MACHINE_LEARNING_INIT(ML_SUPER_HANDLE, &
                                     DIR_APP,         &
                                     DIR_LEN,         &
                                     DYN,             &
                                     INFO,            &
                                     LATT_CUR,        &
                                     T_INFO,          &
                                     IO,              &
                                     VASP_COMM)
      USE base
      USE lattice
      USE poscar
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_IO_WRITE
      USE ML_FF_CONSTANT, ONLY: AUTOA, PI, MUNIT, TUNIT
      USE ml_ff_tutor, only: ml_tutor
      IMPLICIT NONE
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR             ! Bravais lattice information
      TYPE (in_struct)   IO                   ! File units
      TYPE (type_info)   T_INFO               ! type and position information
      CHARACTER(LEN=10), INTENT(IN) ::  DIR_APP
      INTEGER, INTENT(IN)           ::  DIR_LEN
      INTEGER, INTENT(IN)           ::  VASP_COMM
      INTEGER                       ::  INIONS
      INTEGER                       ::  ITYPE

! FOr ISIF.LE.1 the stress tensor is not calculated, which is not allowed
      IF (ML_SUPER_HANDLE%FF%ISTART.LE.1 .AND. DYN%ISIF.LE.1) THEN
         CALL vtutor%error("ML_INTERFACE: ERROR: ML_MODE = train is not possible with ISIF = 0 or 1 &
                           &because the stress tensor is required, exiting ...")
      ENDIF

! One has to provide proper structure name in the POSCAR file otherwise the structure cannot be identified
! when reading the ML_AB file
      IF (T_INFO%SZNAM2.EQ."") THEN
         CALL vtutor%error("ML_INTERFACE: ERROR: No structure name was provided in the POSCAR file, exiting...")
      ENDIF
! Stop code if NSW>1 and ML_ISTART=3
! Give some variables from VASP to MLFF and allocate some necessary arrays
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A = LATT_CUR%A/AUTOA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B = LATT_CUR%B*AUTOA*2.0_q*PI
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A = LATT_CUR%A
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B = LATT_CUR%B
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP = DIR_APP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN = DIR_LEN
      ML_SUPER_HANDLE%FF%DIR_APP = DIR_APP
      ML_SUPER_HANDLE%FF%DIR_LEN = DIR_LEN
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%LABORT = INFO%LABORT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%COMM_VASP = VASP_COMM
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS = T_INFO%NIONS
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2 = TRIM(ADJUSTL(T_INFO%SZNAM2))
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA/AUTOA**3
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM = DYN%POTIM/TUNIT
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM = DYN%POTIM
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP = T_INFO%NTYP
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP))
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP))
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP))
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP = T_INFO%NITYP

      ! Calculate maximum number of atoms per species
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%MITYP = ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP(1)
      DO ITYPE=1, ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%MITYP = MAX(ML_SUPER_HANDLE%ML_INPUT_HANDLE%MITYP, &
                                                 ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP(ITYPE)) 
      ENDDO 

      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS = T_INFO%POMASS/MUNIT
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS = T_INFO%POMASS
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE = T_INFO%TYPE
! Set TIFOR_ML,TSIF_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)
      ENDIF
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS))
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)
      ENDIF
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3))
! Set VEL_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)
      ENDIF
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS))
! Set ionic positions
      IF (ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POSION)) DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POSION)
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POSION(1:3,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS))
      DO INIONS=1, ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%POSION(1:3,INIONS)= &
                  T_INFO%POSION(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                  T_INFO%POSION(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                  T_INFO%POSION(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
      ENDDO
      ! Initialize tutor (only one processor writes messages to screen).
      ml_tutor%discardOutputToNegativeUnit = .TRUE.
      IF (ML_IO_WRITE) THEN
         ml_tutor%unitOut = IO%IU6
         ml_tutor%unitErr = IO%IU0
      ELSE
         ml_tutor%unitOut = -1
         ml_tutor%unitErr = -1
      END IF
! Send variables to machine-learning library
      CALL MACHINE_LEARNING_INIT_LIB (ML_SUPER_HANDLE)
! We need to set this if ISTART=3
      IF (ML_SUPER_HANDLE%FF%ISTART.EQ.3) THEN
         LDO_AB_INITIO=.FALSE.
         LDO_AB_INITIO_NEW=.FALSE.
      ENDIF
      IF (ML_SUPER_HANDLE%FF%ISTART.EQ.4) THEN
         LDO_AB_INITIO=.FALSE.
         LDO_AB_INITIO_NEW=.FALSE.
      ENDIF
    END SUBROUTINE MACHINE_LEARNING_INIT

!=======================================================================
!
! the following routine whenever a new set of positions
! is available to VASP the ML code, might be already able
! to calculate forces etc. and in that case returns the energy
! forces and the stress tensor, as well as an error estimate
! if LDO_AB_INITIO_NEW is set to be true, VASP will do the 
! calculations for the present positions, otherwise
! we will use the predicted data
!
!=======================================================================

    SUBROUTINE ML_TO_VASP_MACHINE_LEARNING(ML_SUPER_HANDLE, DYN, LATT_CUR, &
               NIOND, NSTEP, T_INFO, TOTEN, TIFOR, TSIF, VASP_COMM)
      USE base
      USE classicfields
      USE constant
      USE ini
      USE lattice
      USE poscar
      USE mgrid
      USE mpimy
      USE msymmetry
      USE poscar
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_IO_WRITE
      USE ML_FF_CONSTANT, ONLY: MUNIT, TUNIT
      IMPLICIT NONE
!Input and output variables
      TYPE (ML_SUPER_TYPE)          :: ML_SUPER_HANDLE
      TYPE (dynamics)               :: DYN
      TYPE (latt), INTENT(IN)       :: LATT_CUR          ! Bravais lattice information
      TYPE (type_info), INTENT(IN)  :: T_INFO            ! type and position information
      INTEGER, INTENT(IN)           :: NIOND             ! The number of ions.
      INTEGER, INTENT(IN)           :: NSTEP             ! MD step
      REAL (q), INTENT(OUT)         :: TOTEN             ! total energy
      REAL (q), INTENT(OUT)         :: TIFOR(:,:)        !(1:3,1:NIOND) ! forces on all atoms
      REAL (q), INTENT(OUT)         :: TSIF(:,:)         !(1:3,1:3)      ! stress tensor
      INTEGER, INTENT(IN)           :: VASP_COMM
! Local variables related to symmetry.
      REAL (q)         :: TOTEN_TMP             ! total energy
      REAL (q)         :: TIFOR_TMP(1:3,1:NIOND)        !(1:3,1:NIOND) ! forces on all atoms
      REAL (q)         :: TSIF_TMP(1:3,1:3)         !(1:3,1:3)      ! stress tensor
      LOGICAL :: LDO_AB_INITIO_NEW_TMP
! Local variables related to bias forces.
      INTEGER   :: INIONS
      REAL(q)   :: POSION_LOC(1:3,T_INFO%NIONS)

! Do nothing and return if this flag is set to true (for instance ML_ISTART=3)
      IF (ML_SUPER_HANDLE%SKIP_ML_TO_VASP) THEN
         TOTEN=0.0_q
         TIFOR=0.0_q
         TSIF=0.0_q
         RETURN
      ENDIF
! Give some variables from VASP to MLFF
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A = LATT_CUR%A/AUTOA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B = LATT_CUR%B*AUTOA*2.0_q*PI
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%COMM_VASP = VASP_COMM
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS = T_INFO%NIONS
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP = T_INFO%NITYP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NSTEP = NSTEP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP = T_INFO%NTYP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS = T_INFO%POMASS/MUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM = DYN%POTIM/TUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2 = TRIM(ADJUSTL(T_INFO%SZNAM2))
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE = T_INFO%TYPE
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA/AUTOA**3
! Calculate positions in cartesian coordinates
! Set velocity. At NSTEP>1, velocity at the next step is estimated by adding 
! current accerelation to current velocity.
      DO INIONS=1, T_INFO%NIONS
         POSION_LOC(1:3,INIONS)=T_INFO%POSION(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                                T_INFO%POSION(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                                T_INFO%POSION(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,INIONS)=DYN%VEL(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                                                            DYN%VEL(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                                                            DYN%VEL(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,INIONS)=ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,INIONS) / &
                                                           ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM
      ENDDO
! Predict energy, force and stress tensor by using machine-learning force field.
      IF (ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP) THEN
         CALL ML_TO_VASP_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:,1:),ML_SUPER_HANDLE%ML_INPUT_HANDLE%B(1:,1:), &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NSTEP, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA,POSION_LOC,ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2,TOTEN_TMP,TIFOR_TMP,TSIF_TMP,LDO_AB_INITIO_NEW_TMP)
      ELSE
         CALL ML_TO_VASP_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:,1:),ML_SUPER_HANDLE%ML_INPUT_HANDLE%B(1:,1:), &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NSTEP, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA,POSION_LOC,ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2,TOTEN,TIFOR,TSIF,LDO_AB_INITIO_NEW)
      ENDIF 
    END SUBROUTINE ML_TO_VASP_MACHINE_LEARNING

#ifndef ML_LEARNING_NOT_ACTIVE
!=======================================================================
!
! the following routine passes the energies and forces
! calculated by VASP to the ML code
! note that if LDO_AB_INITIO is false, these are just the
! very same forces and energies as calculated by the ML code
! Not TOTEN, TIFOR, TSIF but TOTEN_ML, TIFOR_ML, TSIF_ML must be
! provided to machine-learning library.
!
!=======================================================================

    SUBROUTINE VASP_TO_ML_MACHINE_LEARNING( ML_SUPER_HANDLE, DIR_APP,DIR_LEN,INFO, LATT_CUR, NSTEP, T_INFO, DYN, IO, GRIDC, IONS, TOTEN, TIFOR, TSIF )
      USE constant
      USE base
      USE mgrid
      USE ini
      USE lattice
      USE poscar
      USE mpimy
      USE classicfields
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE ML_FF_CONSTANT, ONLY: MUNIT, TUNIT, EUNIT, FUNIT, SUNIT
      IMPLICIT NONE
! Input and output variables
      TYPE (ML_SUPER_TYPE)          :: ML_SUPER_HANDLE
      CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
      INTEGER          , INTENT(IN) :: DIR_LEN
      TYPE (info_struct)            :: INFO                      ! Information
      TYPE (latt)      , INTENT(IN) :: LATT_CUR    ! Bravais lattice information
      INTEGER          , INTENT(IN) :: NSTEP       ! MD step
      TYPE (type_info) , INTENT(IN) :: T_INFO      ! type and position information
      TYPE (dynamics)  , INTENT(IN) :: DYN
      TYPE (in_struct) , INTENT(IN) :: IO
      TYPE (grid_3d)   , INTENT(IN) :: GRIDC
      INTEGER                       :: IONS
      REAL(q)                       :: TOTEN
      REAL(q)                       :: TIFOR(:,:) !(1:3,1:IONS)
      REAL(q)                       :: TSIF(:,:) !(1:3,1:3)
! LOCAL variables
      REAL(q)                       :: OFIELD_E                        ! energy contribution from orderfield
      REAL(q)                       :: OFIELD_FOR (1:3,T_INFO%NIONS)    ! forces from order field
      REAL(q)                       :: FAKT
      INTEGER                       :: INIONS
      LOGICAL                       :: LDO_AB_INITIO_HELP
      REAL(q)                       :: POSION_LOC(1:3,T_INFO%NIONS)
      REAL(q)                       :: TOTEN_TMP
      REAL(q)                       :: TIFOR_TMP(1:3,1:T_INFO%NIONS)
      REAL(q)                       :: TSIF_TMP(1:3,1:3)
      REAL(q)                       :: VEL_TMP(1:3,1:T_INFO%NIONS)

! Do nothing and return if this flag is set to true (for instance ML_ISTART=3)
      IF (ML_SUPER_HANDLE%SKIP_VASP_TO_ML) RETURN

! Initialize helping arrays
      VEL_TMP=0.0_q

      ! Conversions needed for input to INIT_LIB
      FAKT=EVTOJ*1.0E+22_q/LATT_CUR%OMEGA !/AUTOA**3

      ! subtract contribution of OFIELD that was added in force.F
      CALL OFIELD(GRIDC%COMM, IO, LATT_CUR, DYN, T_INFO, OFIELD_FOR, OFIELD_E) 
      IF (DYN%PSTRESS/=0) THEN
         ! subtract stress contribution of energy that was added in main.F
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML=TOTEN - DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA - OFIELD_E
      ELSE
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML=TOTEN - OFIELD_E
      ENDIF
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:IONS)=TIFOR(1:3,1:IONS)-OFIELD_FOR(1:3,1:IONS)
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3)=TSIF(1:3,1:3)

! We need to do unit conversions here
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML=ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML/EUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:IONS)=ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:IONS)/FUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3)=ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3)*FAKT/SUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1:3) = LATT_CUR%A(1:3,1:3)/AUTOA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B(1:3,1:3) = LATT_CUR%B(1:3,1:3)*AUTOA*2.0_q*PI
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA/AUTOA**3
      DO INIONS=1, T_INFO%NIONS
         POSION_LOC(1:3,INIONS)=T_INFO%POSION(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                                T_INFO%POSION(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                                T_INFO%POSION(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
      ENDDO
! Do the actual calculations here
      CALL VASP_TO_ML_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE, ML_SUPER_HANDLE%ML_INPUT_HANDLE%A, ML_SUPER_HANDLE%ML_INPUT_HANDLE%B, &
           DIR_APP,DIR_LEN,INFO%LABORT, NSTEP, ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA, POSION_LOC, &
           TRIM(ADJUSTL(T_INFO%SZNAM2)), ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML, ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML, ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML )
    END SUBROUTINE VASP_TO_ML_MACHINE_LEARNING
#endif
! defined(ML_LEARNING_NOT_ACTIVE)

    LOGICAL FUNCTION ML_REQUESTS_EXIT(ML_SUPER_HANDLE)

      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, FF_PAR
      IMPLICIT NONE
      TYPE (ML_SUPER_TYPE), INTENT(IN) :: ML_SUPER_HANDLE

      ML_REQUESTS_EXIT = ML_SUPER_HANDLE%FF%SF_LIMIT_EXCEEDED

    END FUNCTION ML_REQUESTS_EXIT

!=======================================================================
!
! the following routine allows the ML code to terminate and
! write all data
!
!=======================================================================

    SUBROUTINE MACHINE_LEARNING_FINISH (ML_SUPER_HANDLE,DIR_APP,DIR_LEN,NSTEP)
      USE mpimy
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      IMPLICIT NONE
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
      INTEGER          , INTENT(IN) :: DIR_LEN
      INTEGER          , INTENT(IN) :: NSTEP ! Final MD step
      CALL MACHINE_LEARNING_FINISH_LIB (ML_SUPER_HANDLE,DIR_APP,DIR_LEN,NSTEP)
! Deallocate NITYP, POMASS and TYPE
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP)
      ENDIF
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS)
      ENDIF
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE)
      ENDIF
! Deallocate POSEIDON
      IF (ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POSION)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POSION)
      ENDIF
! Deallocate TIFOR_ML and TSIF_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)
      ENDIF
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)
      ENDIF
! Deallocate VEL_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)
      ENDIF
! Deallocate arrays storing the energies of isolated atoms.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%EATOM_VASP)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%EATOM_VASP)
      ENDIF
! Deallocate ML_ICOUPLE that stores the list of atoms, where the coupling parameter is introduced for the calculation of the chemical potential.
      IF(ALLOCATED(ML_SUPER_HANDLE%FFM%ICOUPLE)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%FFM%ICOUPLE)
      ENDIF
    END SUBROUTINE MACHINE_LEARNING_FINISH

#endif
! defined(MPI) || defined(MPI_CHAIN)
#endif
! ML_AVAILABLE
END MODULE ml_interface
