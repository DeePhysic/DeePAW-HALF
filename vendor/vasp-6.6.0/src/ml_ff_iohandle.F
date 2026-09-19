#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

!****************************************************************************************************
! Module to handle files used in machine-learning calculations
! 1) ABCAR is an input file storing all the necessary ab initio data for generating force field.
! 2) INCAR is an input file that contains force field parameters.
! 3) ABNCAR is an output file that should contain all the renewed ab initio data.
! 4) REGCAR is an output file that contains the results of the regressions.
!****************************************************************************************************

      MODULE IOHANDLE
        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        USE LATTICE_ML
        USE MPI_DATA
        USE TIME_CHECK
        IMPLICIT NONE

        CONTAINS

!****************************************************************************************************
! Increment global file_stream counter 
!****************************************************************************************************

        SUBROUTINE ML_FILE_STREAM_INCREMENT(FILE_STREAM)
          IMPLICIT NONE
          INTEGER, INTENT(INOUT) :: FILE_STREAM
          FILE_STREAM = FILE_STREAM + 1
        END SUBROUTINE ML_FILE_STREAM_INCREMENT

!****************************************************************************************************
! Set ending of filenames for given machine learning instance
!****************************************************************************************************

        SUBROUTINE ML_SET_POSTFIX(INSTANCE,FILE_POSTFIX)
          IMPLICIT NONE
          INTEGER, INTENT(IN)           :: INSTANCE
          CHARACTER(LEN=4), INTENT(OUT) :: FILE_POSTFIX
          CHARACTER(LEN=1)              :: H1
          CHARACTER(LEN=2)              :: H2
          CHARACTER(LEN=3)              :: H3
          IF (INSTANCE.EQ.1) THEN
             FILE_POSTFIX=""
          ELSE
             IF (INSTANCE.LT.10) THEN
                WRITE(H1,"(I1)") INSTANCE
                FILE_POSTFIX="_" // TRIM(H1)
             ELSE IF (INSTANCE.GT.9.AND.INSTANCE.LT.100) THEN
                WRITE(H2,"(I2)") INSTANCE
                FILE_POSTFIX="_" // TRIM(H2)
             ELSE IF (INSTANCE.GT.99.AND.INSTANCE.LT.1000) THEN
                WRITE(H3,"(I3)") INSTANCE
                FILE_POSTFIX="_" // TRIM(H3)
             ENDIF
          ENDIF
        END SUBROUTINE ML_SET_POSTFIX

!****************************************************************************************************
! Set file name in code
!****************************************************************************************************

        SUBROUTINE SET_FILE (CFILE_IN,IFILE_IN,CFILE_OUT,IFILE_OUT)
          IMPLICIT NONE
          CHARACTER(len=*), INTENT(IN)  :: CFILE_IN
          CHARACTER(len=*), INTENT(OUT) :: CFILE_OUT
          INTEGER IFILE_IN
          INTEGER IFILE_OUT
          CFILE_OUT=CFILE_IN
          IFILE_OUT=IFILE_IN
        END SUBROUTINE SET_FILE

        !======================================================================
        !> Open ML_FF, check version and skip header (if present).
        !======================================================================
        SUBROUTINE ROPEN_ML_FF(FF)

           USE ML_FF_CONSTANT
           USE ml_ff_tutor, ONLY: ml_tutor
           USE ml_ff_string, ONLY: str
           USE ML_FF_STRUCT, ONLY: FF_PAR
           IMPLICIT NONE
           TYPE(FF_PAR), INTENT(INOUT) :: FF

           ! Local variables
           CHARACTER(LEN=255) :: STMP
           INTEGER            :: IERR
           INTEGER            :: I

           STMP = ""

           OPEN(UNIT=FF%IU1,                               &
                FILE=FF%DIR_APP(1:FF%DIR_LEN) // FF%FFCAR, &
                STATUS='OLD',                              &
                FORM='UNFORMATTED',                        &
                ACCESS='STREAM',                           &
                IOSTAT=IERR)
           IF (IERR.NE.0) THEN
              CALL ml_tutor%error("Cannot open ML_FF file for reading.")
           ENDIF

           ! Check if header exists (version 0.1.0 had no header, for newer
           ! versions the user may have removed it).
           FF%ML_FF_HAS_HEADER = .FALSE.
           ! Read first 6 bytes.
           READ(FF%IU1) STMP(1:6)
           ! Check first six characters.
           IF (STMP(1:6) == "ML_FF ") FF%ML_FF_HAS_HEADER = .TRUE.

           IF (FF%ML_FF_HAS_HEADER) THEN
              ! Continue reading file to determine version.
              STMP = ""
              CALL GET_NEXT_STRING(FF%IU1, STMP, I, ".")
              READ(STMP(1:I-1), *) FF%ML_FF_VERSION_READ(1)
              CALL GET_NEXT_STRING(FF%IU1, STMP, I, ".")
              READ(STMP(1:I-1), *) FF%ML_FF_VERSION_READ(2)
              CALL GET_NEXT_STRING(FF%IU1, STMP, I, " ")
              READ(STMP(1:I-1), *) FF%ML_FF_VERSION_READ(3)

              ! Check if this ML_FF version is supported.
              CALL CHECK_ML_FF_VERSION_COMPATIBILITY(FF%ML_FF_VERSION_READ, &
                                                     ML_FF_VERSION_WRITE)

              ! Check if main file contents are binary or text.
              STMP = ""
              CALL GET_NEXT_STRING(FF%IU1, STMP, I, " ")
              IF      (STMP(1:I-1) == "binary") THEN
                 FF%ML_FF_BINARY = .TRUE.
              !ELSE IF (STMP(1:I-1) == "txt") THEN
              !   FF%ML_FF_BINARY = .FALSE.
              ELSE
                 CALL ml_tutor%error("Unknown ML_FF file mode or corrupt &
                                     &header.")
              END IF

              ! Skip everything up to next colon, STMP should then be '{ "date" :'.
              STMP = ""
              CALL GET_NEXT_STRING(FF%IU1, STMP, I, ":")

              ! Search for next comma, STMP should contain the date value (with a space and quote
              ! before and a quote after it).
              STMP = ""
              CALL GET_NEXT_STRING(FF%IU1, STMP, I, ",")
              FF%ML_FF_CREATION_DATETIME = STMP(3:I-2)
           END IF

           ! Go back to beginning of file.
           REWIND(FF%IU1)

           ! If header is present, skip it.
           IF (FF%ML_FF_HAS_HEADER) THEN
              DO I = 1, ML_FF_HEADER_SIZE
                 READ(FF%IU1) STMP(1:1)
              END DO
           END IF

           CONTAINS

           SUBROUTINE GET_NEXT_STRING(IU, STRING, I, DELIMITER)

              IMPLICIT NONE
              INTEGER         , INTENT(IN)    :: IU
              CHARACTER(LEN=*), INTENT(INOUT) :: STRING
              INTEGER         , INTENT(INOUT) :: I
              CHARACTER(LEN=1), INTENT(IN)    :: DELIMITER

              DO I = 1, LEN(STRING)
                 READ(IU) STRING(I:I)
                 IF (STMP(I:I) == DELIMITER) EXIT
              END DO

              IF (I == 1 .OR. I >= LEN(STRING)) THEN
                 CALL ml_tutor%error("Cannot find next string, maybe the ML_FF &
                                     &file header is corrupted.")
              END IF

           END SUBROUTINE GET_NEXT_STRING

        END SUBROUTINE ROPEN_ML_FF

        !======================================================================
        !> Check Input file types.
        !>
        !> Checks the input file types and gives error if two types have the
        !> same name (e.g. O O is not allowed, but O1 O2 would work)
        !======================================================================
        SUBROUTINE CHECK_TYPES_INPUT_STRUCTURE(TYPE_NAME,NTYP)
           USE ml_ff_tutor, ONLY: ml_tutor
           USE ml_ff_string, ONLY: str
           IMPLICIT NONE
           ! Input variables
           CHARACTER(LEN=2), ALLOCATABLE :: TYPE_NAME(:)
           INTEGER                       :: NTYP
           ! Local variables
           INTEGER                       :: ITYP, JTYP
           LOGICAL                       :: IS_DA_SAME
             
           IS_DA_SAME=.FALSE.
           DO ITYP=1, NTYP
              DO JTYP=1, NTYP
                 IF (ITYP.NE.JTYP) THEN
                    IF (TYPE_NAME(ITYP).EQ.TYPE_NAME(JTYP)) THEN
                       IS_DA_SAME=.TRUE.
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
           
           IF (IS_DA_SAME) THEN
              CALL ml_tutor%error("Using the same type name repeatedly in your &
                                  &structure input file (POSCAR) is not &
                                  &allowed and the code needs to exit. Please &
                                  &change the input such that each type has a &
                                  &unique name (limited to two characters). &
                                  &E.g. rename 'O O' to 'O1 O2', or group &
                                  &elements with the same type&
                                  & together. ")
           ENDIF
        END SUBROUTINE CHECK_TYPES_INPUT_STRUCTURE

        !======================================================================
        !> Check ML_FF version compatibility.
        !>
        !> For version history see documentation of ML_FF_VERSION_WRITE in ml_ff_constants.F
        !======================================================================
        SUBROUTINE CHECK_ML_FF_VERSION_COMPATIBILITY(VFILE, VMAX)

           USE ml_ff_tutor, ONLY: ml_tutor
           USE ml_ff_string, ONLY: str
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: VFILE(3)
           INTEGER, INTENT(IN) :: VMAX(3)

           LOGICAL :: SUPPORTED

           IF (ANY(VFILE < 0)) THEN
              CALL ml_tutor%error("Negative ML_FF version number encountered, &
                                  &maybe the file is corrupted.")
           END IF

           SUPPORTED = .TRUE.
           IF (VFILE(1) > VMAX(1)) THEN
              SUPPORTED = .FALSE.
           ELSE
              IF (VFILE(1) == VMAX(1) .AND. VFILE(2) > VMAX(2)) THEN
                 SUPPORTED = .FALSE.
              ELSE
                 IF (VFILE(1) == VMAX(1) .AND. &
                     VFILE(2) == VMAX(2) .AND. &
                     VFILE(3)  > VMAX(3)) SUPPORTED = .FALSE.
              END IF
           END IF
           IF (.NOT. SUPPORTED) THEN
              CALL ml_tutor%error("Found ML_FF version " //                   &
                 str(VFILE(1)) // "." //                                      &
                 str(VFILE(2)) // "." //                                      &
                 str(VFILE(3)) // " which is not supported by this executable &
                 &(maximum is " //                                            &
                 str(VMAX(1)) // "." //                                       &
                 str(VMAX(2)) // "." //                                       &
                 str(VMAX(3)) // "), please consider updating VASP.")
           END IF

        END SUBROUTINE CHECK_ML_FF_VERSION_COMPATIBILITY

!****************************************************************************************************
!> This routine determines whether the ML_FF file is in binary format or formatted string format.
!****************************************************************************************************
        SUBROUTINE CHECK_ML_FF_FORMAT(FF)

           USE ML_FF_STRUCT, ONLY: FF_PAR
           USE ml_ff_tutor, ONLY: ml_tutor
           IMPLICIT NONE
           ! Input and output variables
           TYPE (FF_PAR) :: FF

           ! Local variables
           CHARACTER :: HELP_CHAR
           INTEGER   :: IERR

           ! Assume binary format for ML_FF file by default.
           FF%ML_FF_BINARY = .TRUE.
           ! Open ML_FF file as unformatted.
           OPEN(UNIT=FF%IU1,                               &
                FILE=FF%DIR_APP(1:FF%DIR_LEN) // FF%FFCAR, &
                STATUS='OLD',                              &
                FORM='UNFORMATTED',                        &
                ACCESS='STREAM',                           &
                IOSTAT=IERR)
           ! IF there is a problem opening ML_FF then we need to abort with error message
           IF (IERR.NE.0) THEN
              CALL ml_tutor%error("Cannot open ML_FF file for reading.")
           ENDIF
           ! READ first character
           READ(FF%IU1) HELP_CHAR
           ! Non-binary ML_FF file starts with a line full of "*".
           IF (HELP_CHAR == "*") THEN
              FF%ML_FF_BINARY = .FALSE.
              ! Text version without header is defined as version 0.0.1.
              FF%ML_FF_VERSION_READ = [0, 0, 1]
           ENDIF

           CLOSE(FF%IU1)

        END SUBROUTINE CHECK_ML_FF_FORMAT

!****************************************************************************************************
! Read element number written in FFCAR file.
!****************************************************************************************************

        SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_1(FF,MTYP_OLD)
          USE ML_FF_CONSTANT, ONLY: ML_FF_VERSION_WRITE
          USE ML_FF_STRUCT, ONLY: FF_PAR
          IMPLICIT NONE
          ! Input and output variables
          TYPE (FF_PAR)                   :: FF
! Output variables.
          INTEGER         , INTENT(OUT)   :: MTYP_OLD
! Local variables
          INTEGER                         :: IHELP(1:3)
          INTEGER                         :: IDUMMY
          LOGICAL                         :: DO_READ_DESC_TYPE
! Open file.
          CALL ROPEN_ML_FF(FF)
! Do actual reading
          READ(FF%IU1) IHELP(1), IHELP(2), IHELP(3)
          ! Check version support.
          CALL CHECK_ML_FF_VERSION_COMPATIBILITY(IHELP, ML_FF_VERSION_WRITE)
          DO_READ_DESC_TYPE=.TRUE.
          IF (IHELP(1).EQ.0) THEN
             IF (IHELP(2).EQ.2) THEN
                IF (IHELP(3).EQ.0) THEN
                   DO_READ_DESC_TYPE=.FALSE.
                ENDIF
             ELSE IF (IHELP(2).LT.2) THEN
                DO_READ_DESC_TYPE=.FALSE.
             ENDIF
          ENDIF
          IF (DO_READ_DESC_TYPE) THEN
             ! Read descriptor type
             READ(FF%IU1) IDUMMY
          ENDIF
          ! Read number of elements
          READ(FF%IU1) MTYP_OLD
! Close file.
          CLOSE(FF%IU1)
        END SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_1

!****************************************************************************************************
! Read element data written in FFCAR file.
!****************************************************************************************************
        SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_2(FF,EATOM_OLD,MTYP_OLD,POMASS_OLD,TYPE_OLD)
          USE ML_FF_CONSTANT, ONLY: ML_FF_VERSION_WRITE
          USE ML_FF_STRUCT, ONLY: FF_PAR
          IMPLICIT NONE
! Input variables.
          TYPE (FF_PAR)                   :: FF
          INTEGER         , INTENT(IN)    :: MTYP_OLD
! Input and ouput variables.
          REAL(q)         , INTENT(OUT)   :: EATOM_OLD(:) !(1:MTYP_OLD)
          REAL(q)         , INTENT(OUT)   :: POMASS_OLD(:) !(1:MTYP_OLD)
          CHARACTER(LEN=2), INTENT(OUT)   :: TYPE_OLD(:) !(1:MTYP_OLD)
! Local variables
          INTEGER                         :: ICOL
          INTEGER                         :: IROW
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NWRITE
          INTEGER                         :: IDUMMY
          INTEGER                         :: IHELP(1:3)
          LOGICAL                         :: DO_READ_DESC_TYPE
! Set constant
          NWRITE = 3
! Open file.
          CALL ROPEN_ML_FF(FF)
! Do actual reading
          READ(FF%IU1) IHELP(1), IHELP(2), IHELP(3)
          ! Check version support.
          CALL CHECK_ML_FF_VERSION_COMPATIBILITY(IHELP, ML_FF_VERSION_WRITE)
          DO_READ_DESC_TYPE=.TRUE.
          IF (IHELP(1).EQ.0) THEN
             IF (IHELP(2).EQ.2) THEN
                IF (IHELP(3).EQ.0) THEN
                   DO_READ_DESC_TYPE=.FALSE.
                ENDIF
             ELSE IF (IHELP(2).LT.2) THEN
                DO_READ_DESC_TYPE=.FALSE.
             ENDIF
          ENDIF
          IF (DO_READ_DESC_TYPE) THEN
             ! Read descriptor type
             READ(FF%IU1) IDUMMY
          ENDIF
          ! Read number of elements
          READ(FF%IU1) IDUMMY
          NCOL=MTYP_OLD/NWRITE
          NMOD=MOD(MTYP_OLD,NWRITE)
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (TYPE_OLD(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (TYPE_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (EATOM_OLD(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (EATOM_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO ICOL=1, NCOL
             READ(FF%IU1) (POMASS_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(FF%IU1) (POMASS_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
! Unit conversion
          EATOM_OLD(1:MTYP_OLD)=EATOM_OLD(1:MTYP_OLD)/EUNIT
          POMASS_OLD(1:MTYP_OLD)=POMASS_OLD(1:MTYP_OLD)/MUNIT
! Close file.
          CLOSE(FF%IU1)
        END SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_2

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Read the number of configurations in the ABCAR file.
!****************************************************************************************************
         SUBROUTINE READ_ABCAR_HELP (ABCAR,IU,FF,AB_HELP)
            USE ML_FF_STRUCT, ONLY: FF_PAR, ABINITIO
            USE ml_ff_tutor, ONLY: ml_tutor
            IMPLICIT NONE
! Input variables.
            CHARACTER(LEN=9) , INTENT(IN)    :: ABCAR
            INTEGER          , INTENT(IN)    :: IU
            TYPE (FF_PAR)                    :: FF
! Output variables
            TYPE (ABINITIO)                  :: AB_HELP
! Local variables
            INTEGER                          :: IB
            INTEGER                          :: ICOL
            INTEGER                          :: ICONF
            INTEGER                          :: IDUMMY(1:2)
            INTEGER                          :: IERR
            INTEGER                          :: INIONS
            INTEGER                          :: INTYP
            INTEGER                          :: IROW
            INTEGER                          :: IXYZ
            INTEGER                          :: JXYZ
            INTEGER                          :: NATOM_HELP
            INTEGER                          :: NCOL
            INTEGER                          :: NMOD
            INTEGER                          :: NWRITE
            CHARACTER(LEN=40)                :: READ_HELPER
            REAL(q)                          :: VERSION
! Set constants.
            NWRITE=3
! Open file.
            OPEN(UNIT=IU,FILE=FF%DIR_APP(1:FF%DIR_LEN)//ABCAR, &
                 STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
            IF (IERR.NE.0) THEN
               CALL ml_tutor%error("Cannot open ML_AB file for reading.")
            ENDIF
! Start reading
            ! Version number
            READ(IU,*) AB_HELP%VERSION
            ! Read: The number of configurations
            CALL READ_THREE_EMPTY_LINES(IU)
            READ(IU,*) AB_HELP%NCONF
            ! Set: AB_HELP%MCONF
            AB_HELP%MCONF=AB_HELP%NCONF
            ! Read: The maximum number of atom type
            CALL READ_THREE_EMPTY_LINES(IU)
            READ(IU,*) AB_HELP%MTYP
            ! Set: Line modifier for reading
            NCOL=AB_HELP%MTYP/NWRITE
            NMOD=MOD(AB_HELP%MTYP,NWRITE)
            ! Allocate: The atom types in the data file
            ALLOCATE(AB_HELP%TYPE(AB_HELP%MTYP))
            ! Read: The atom types in the data file
            CALL READ_THREE_EMPTY_LINES(IU)
            DO ICOL=1, NCOL
               READ(IU,*) (AB_HELP%TYPE(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
            ENDDO
            IF(NMOD.NE.0) THEN
               READ(IU,*) (AB_HELP%TYPE(NWRITE*NCOL+IROW),IROW=1, NMOD)
            ENDIF
            ! Read: The maximum number of atoms per system
            CALL READ_THREE_EMPTY_LINES(IU)
            READ(IU,*) AB_HELP%MIONS
            ! Read: The maximum number of atoms per atom type
            CALL READ_THREE_EMPTY_LINES(IU)
            READ(IU,*) AB_HELP%MITYP
            ! Allocate: Reference atomic energy (eV)
            ALLOCATE(AB_HELP%EATOM(AB_HELP%MTYP))
            ! Read: Reference atomic energy (eV)
            CALL READ_THREE_EMPTY_LINES(IU)
            DO ICOL=1, NCOL
               READ(IU,*) (AB_HELP%EATOM(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
            ENDDO
            IF(NMOD.NE.0) THEN
               READ(IU,*) (AB_HELP%EATOM(NWRITE*NCOL+IROW),IROW=1, NMOD)
            ENDIF
            ! Allocate: Atomic mass
            ALLOCATE(AB_HELP%POMASS(AB_HELP%MTYP))
            ! Read: Atomic mass
            CALL READ_THREE_EMPTY_LINES(IU)
            DO ICOL=1, NCOL
               READ(IU,*) (AB_HELP%POMASS(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
            ENDDO
            IF(NMOD.NE.0) THEN
               READ(IU,*) (AB_HELP%POMASS(NWRITE*NCOL+IROW),IROW=1, NMOD)
            ENDIF
            ! Allocate: The numbers of basis sets per atom type
            ALLOCATE(AB_HELP%NBASIS(AB_HELP%MTYP))
            ! Read: The numbers of basis sets per atom type
            CALL READ_THREE_EMPTY_LINES(IU)
            DO ICOL=1, NCOL
               READ(IU,*) (AB_HELP%NBASIS(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
            ENDDO
            IF(NMOD.NE.0) THEN
               READ(IU,*) (AB_HELP%NBASIS(NWRITE*NCOL+IROW),IROW=1, NMOD)
            ENDIF
            ! Allocate: Basis set for
            AB_HELP%NB_MAX=AB_HELP%NBASIS(1)
            DO INTYP=1, AB_HELP%MTYP
               AB_HELP%NB_MAX=MAX(AB_HELP%NB_MAX,AB_HELP%NBASIS(INTYP))
            ENDDO
            ALLOCATE(AB_HELP%LNCONF(AB_HELP%NB_MAX,AB_HELP%MTYP))
            ALLOCATE(AB_HELP%LNIONS_FF(AB_HELP%NB_MAX,AB_HELP%MTYP))
            ! Read: Basis set for
            DO INTYP=1, AB_HELP%MTYP
               CALL READ_THREE_EMPTY_LINES(IU)
               IF (FF%ISTART .EQ. 3) THEN
                  ! In this mode we do only dummy reading
                  DO IB=1, AB_HELP%NBASIS(INTYP)
                     READ(IU,*) IDUMMY(1), IDUMMY(2)
                     AB_HELP%LNCONF(IB,INTYP)=1
                     AB_HELP%LNIONS_FF(IB,INTYP)=1
                  ENDDO
               ELSE
                  DO IB=1, AB_HELP%NBASIS(INTYP)
                     READ(IU,*) AB_HELP%LNCONF(IB,INTYP), AB_HELP%LNIONS_FF(IB,INTYP)
                     ! Exit if ML_AB has too many elements
                     IF (AB_HELP%LNCONF(IB,INTYP) .GT. AB_HELP%NCONF) THEN
                        CALL ml_tutor%error("ERROR: Trainings &
                              &configuration index within list of &
                              &local reference configurations in &
                              &ML_AB file is out of range &
                              &of the maximum number of training &
                              &configurations. Please correct &
                              &that in the &
                              &ML_AB file. Exiting...")
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO 
            ! Allocate: Arrays for training configurations
            ALLOCATE(AB_HELP%A(1:3,1:3,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%B(1:3,1:3,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%CTIFOR_ALLCONF_OLD(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%MITYP_CONF(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%NITYP(1:AB_HELP%MTYP,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%NTYP(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%NIONS(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%OMEGA(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%POSION(1:3,1:AB_HELP%MIONS,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%SZNAM2_INPUT(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%SZNAM2(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%TIFOR(1:3,1:AB_HELP%MIONS,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%TOTEN(1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%TSIF(1:3,1:3,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%TYPE_CONF(1:AB_HELP%MTYP,1:AB_HELP%MCONF))
            ALLOCATE(AB_HELP%CTIFOR_PRESENT(1:AB_HELP%MCONF))
            ! Read: Arrays for training configurations
            DO ICONF=1, AB_HELP%NCONF
               ! Read: Configuration num.
               READ(IU,*)
               READ(IU,*)
               ! Read: System name
               CALL READ_THREE_EMPTY_LINES(IU)
               READ(IU,'(A40)') AB_HELP%SZNAM2_INPUT(ICONF)
               AB_HELP%SZNAM2_INPUT(ICONF)=TRIM(ADJUSTL(AB_HELP%SZNAM2_INPUT(ICONF)))
               ! Read: The number of atom types
               CALL READ_THREE_EMPTY_LINES(IU)
               READ(IU,*) AB_HELP%NTYP(ICONF)
               ! Read: The number of atoms
               CALL READ_THREE_EMPTY_LINES(IU)
               READ(IU,*) AB_HELP%NIONS(ICONF)
               ! Read: Atom types and atom numbers
               CALL READ_THREE_EMPTY_LINES(IU)
               NATOM_HELP = 0
               DO INTYP=1, AB_HELP%NTYP(ICONF)
                  READ(IU,*) AB_HELP%TYPE_CONF(INTYP,ICONF), AB_HELP%NITYP(INTYP,ICONF)
                  ! Gracefully exit with error if number of atoms is larger than AB_HELP%MITYP
                  IF (AB_HELP%NITYP(INTYP,ICONF) .GT. AB_HELP%MITYP) THEN
                     CALL ml_tutor%error("ERROR: Number of atoms &
                          &for a given type &
                          &within training configurations is larger &
                          &than the maximum number of atoms &
                          &per atom type specified at the &
                          &beginning of ML_AB file. Please &
                          &adjust that. Exiting...")
                  ENDIF
                  NATOM_HELP = NATOM_HELP + AB_HELP%NITYP(INTYP,ICONF)
               ENDDO
               ! Check if number of atoms in this configuration exceeds max number of atoms
               IF (NATOM_HELP .GT. AB_HELP%MIONS) THEN
                  CALL ml_tutor%error("ERROR: Total number of atoms &
                       &within training configurations is larger &
                       &than the maximum number of atoms per system &
                       &specified at the &
                       &beginning of ML_AB file. Please &
                       &adjust that. Exiting...")
               ENDIF
               ! Read: Optional reading of CTIFOR
               READ(IU,*)
               READ(IU,*) READ_HELPER
               READ(IU,*)
               READ_HELPER=ADJUSTL(TRIM(READ_HELPER))
               IF (READ_HELPER(1:1).EQ.'C') THEN
                  AB_HELP%CTIFOR_PRESENT(ICONF)=.TRUE.
                  READ(IU,*) AB_HELP%CTIFOR_ALLCONF_OLD(ICONF)
                  CALL READ_THREE_EMPTY_LINES(IU)
               ELSE
                  AB_HELP%CTIFOR_PRESENT(ICONF)=.FALSE.
               ENDIF
               ! Read: Primitive lattice vectors (ang.)
               DO IXYZ=1, 3
                  READ(IU,*) (AB_HELP%A(JXYZ,IXYZ,ICONF),JXYZ=1, 3)
               ENDDO
               ! Read: Atomic positions (ang.)
               CALL READ_THREE_EMPTY_LINES(IU)
               DO INIONS=1, AB_HELP%NIONS(ICONF)
                  READ(IU,*) (AB_HELP%POSION(IXYZ,INIONS,ICONF),IXYZ=1, 3)
               ENDDO

               ! End optional actions for CTIFOR
               !Obtain Bravais Lattice
               DO IXYZ=1, 3
                  DO JXYZ=1, 3
                     AB_HELP%A(JXYZ,IXYZ,ICONF)=AB_HELP%A(JXYZ,IXYZ,ICONF)/AUTOA
                  ENDDO
               ENDDO
               !! shift atoms back into primitve cell
               ! compute direct lattice to shift atoms back to primitive cell
               CALL DIRECT_LATTICE( AB_HELP%A( 1:, 1:, ICONF ), &
                                    AB_HELP%B( 1:, 1:, ICONF ), &
                                    AB_HELP%OMEGA(ICONF ) )
               ! compute direct cordinates of positions
               CALL KARDIR( AB_HELP%NIONS( ICONF ), AB_HELP%POSION( :, :, ICONF ),&
                            AB_HELP%B( 1:, 1:, ICONF ) )
               ! shift direct positions back into the primitive box
               CALL TOPRIM( AB_HELP%NIONS( ICONF ), AB_HELP%POSION( :, :, ICONF ) )
               ! rescale positions back to cartesian coordinates
               CALL DIRKAR( AB_HELP%NIONS( ICONF ), AB_HELP%POSION( :, :, ICONF ),&
                            AB_HELP%A( 1:, 1:, ICONF ) )
               ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! Obtain reciprocal lattice
               CALL RECIPROCAL_LATTICE( AB_HELP%A( 1:, 1:, ICONF), &
                                        AB_HELP%B( 1:, 1:, ICONF), &
                                        AB_HELP%OMEGA( ICONF ) )

               ! Determine the largest volume amongst training structures
               IF (ICONF.EQ.1) THEN
                  AB_HELP%OMEGA_MAX = AB_HELP%OMEGA(ICONF)
               ELSE
                  AB_HELP%OMEGA_MAX = MAX(AB_HELP%OMEGA_MAX, AB_HELP%OMEGA(ICONF))
               ENDIF

               ! Read: Total energy (eV)
               CALL READ_THREE_EMPTY_LINES(IU)
               READ(IU,*) AB_HELP%TOTEN(ICONF)
               ! Read: Forces (eV ang.^-1)
               CALL READ_THREE_EMPTY_LINES(IU)
               DO INIONS=1, AB_HELP%NIONS(ICONF)
                  READ(IU,*) (AB_HELP%TIFOR(IXYZ,INIONS,ICONF),IXYZ=1, 3)
               ENDDO
               ! Read: Stress (kbar)
               READ(IU,*)
               READ(IU,*)
               ! Read: XX YY ZZ
               CALL READ_THREE_EMPTY_LINES(IU)
               READ(IU,*) AB_HELP%TSIF(1,1,ICONF),AB_HELP%TSIF(2,2,ICONF),AB_HELP%TSIF(3,3,ICONF) 
               ! Read: XY YZ ZX
               CALL READ_THREE_EMPTY_LINES(IU)
               READ(IU,*) AB_HELP%TSIF(2,1,ICONF),AB_HELP%TSIF(3,2,ICONF),AB_HELP%TSIF(3,1,ICONF) 
               ! Read: Possibly skip the lines for charges in very old versions (0.0) of ML_AB file
               IF (AB_HELP%VERSION.EQ.(0.0_q)) THEN
                  CALL READ_THREE_EMPTY_LINES(IU)
                  DO INIONS=1, AB_HELP%NIONS(ICONF)
                     READ(IU,*)
                  ENDDO
               ENDIF
            ENDDO             
            ! Close the file handle
            CLOSE(IU)
         END SUBROUTINE READ_ABCAR_HELP

         SUBROUTINE READ_THREE_EMPTY_LINES(IU)
            IMPLICIT NONE
            INTEGER, INTENT(IN)    :: IU
            INTEGER                :: I
            DO I=1,3
               READ(IU,*)
            ENDDO
         END SUBROUTINE READ_THREE_EMPTY_LINES

!****************************************************************************************************
! Read element number written in ABCAR file.
!****************************************************************************************************
         SUBROUTINE CHECK_ABCAR_VALIDITY (AB,FF)
            USE ML_FF_STRUCT, ONLY: ABINITIO, FF_PAR, ML_IO_WRITE
            IMPLICIT NONE
! Input and output variables
            TYPE (ABINITIO)                  :: AB
            TYPE (FF_PAR)                    :: FF
! Functions used
! Local variables
            INTEGER                          :: ICONF
            INTEGER                          :: IFLAG
            INTEGER                          :: JCONF
            INTEGER                          :: NSYS
            INTEGER                          :: SYS_HELP(1:AB%NCONF)

            ! Set first system
            NSYS=1
            SYS_HELP(1)=NSYS
            ! Set second system
            IF (AB%NCONF.GE.2) THEN
               IF (.NOT.CHECK_TWO_SYSTEMS(2,1,AB,FF)) THEN
                  NSYS=NSYS+1
               ENDIF
               SYS_HELP(2)=NSYS
            ENDIF
            ! Set remaining systems
            IF (AB%NCONF.GE.3) THEN
               DO ICONF=3, AB%NCONF
                  IFLAG=1
                  DO JCONF=1, ICONF-1
                     IF (CHECK_TWO_SYSTEMS(ICONF,JCONF,AB,FF)) THEN
                        SYS_HELP(ICONF)=SYS_HELP(JCONF)
                        IFLAG=0
                        EXIT 
                     ENDIF
                  ENDDO
                  IF (IFLAG.EQ.1) THEN
                     NSYS=NSYS+1
                     SYS_HELP(ICONF)=NSYS
                  ENDIF
               ENDDO
               ! Check systems
               IFLAG=1
               DO ICONF=3,AB%NCONF
                  IF (SYS_HELP(ICONF).NE.SYS_HELP(ICONF-1)) THEN
                     DO JCONF=1, ICONF-2
                        IF (SYS_HELP(ICONF).EQ.SYS_HELP(JCONF)) THEN
                           IFLAG=0
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               IF (IFLAG.EQ.0) THEN
                  IF (ML_IO_WRITE) THEN
                      WRITE(*,*) "WARNING: Your ML_AB file &
                        will be rearranged. Some of the structures &
                        are not ordered in the right way."
                   ENDIF
               ENDIF
            ENDIF
         END SUBROUTINE CHECK_ABCAR_VALIDITY

!****************************************************************************************************
! This function checks if two traning structures belong to the same system
! .FALSE.=different systems
! .TRUE.=same system
!****************************************************************************************************

         FUNCTION CHECK_TWO_SYSTEMS(INDEX_NEW,INDEX_PREV,AB,FF)
            USE ML_FF_STRUCT, ONLY: ABINITIO, FF_PAR
            IMPLICIT NONE
            LOGICAL                          :: CHECK_TWO_SYSTEMS
            ! Input variables
            INTEGER                          :: INDEX_NEW
            INTEGER                          :: INDEX_PREV
            TYPE (ABINITIO)                  :: AB
            TYPE (FF_PAR)                    :: FF
            ! Local variables
            INTEGER                          :: IFLAG
            INTEGER                          :: INTYP
            
            ! Check number of element types
            IF (AB%NTYP(INDEX_NEW).NE.AB%NTYP(INDEX_PREV)) THEN 
               CHECK_TWO_SYSTEMS=.FALSE.
               RETURN
            ENDIF
            ! Check number of atoms 
            IF (AB%NIONS(INDEX_NEW).NE.AB%NIONS(INDEX_PREV)) THEN
               CHECK_TWO_SYSTEMS=.FALSE.
               RETURN
            ENDIF
            ! Check individual type related stuff within structures
            DO INTYP=1, AB%NTYP(INDEX_NEW) 
               ! If types are not in the same order they are treated as different systems
               IF (AB%TYPE_CONF(INTYP,INDEX_NEW) .NE. AB%TYPE_CONF(INTYP,INDEX_PREV)) THEN
                  CHECK_TWO_SYSTEMS=.FALSE.
                  RETURN
               ENDIF
               ! If number of atoms per type are different, then they are different structure 
               IF (AB%NITYP(INTYP,INDEX_NEW).NE.AB%NITYP(INTYP,INDEX_PREV)) THEN
                  CHECK_TWO_SYSTEMS=.FALSE.
                  RETURN
               ENDIF
            ENDDO
            ! Check structure names
            IF (FF%LUSE_NAMES) THEN
               IF (AB%SZNAM2_INPUT(INDEX_NEW).NE.AB%SZNAM2_INPUT(INDEX_PREV)) THEN
                  CHECK_TWO_SYSTEMS=.FALSE.
                  RETURN
               ENDIF
            ENDIF
            ! Otherwise the two systems are the same
            CHECK_TWO_SYSTEMS=.TRUE.
            RETURN
         END FUNCTION CHECK_TWO_SYSTEMS

!****************************************************************************************************
! This function checks if two traning structures belong to the same system
! .FALSE.=different systems
! .TRUE.=same system
!****************************************************************************************************
         SUBROUTINE REARRANGE_ML_AB(ML_INPUT_HANDLE,AB,FF)
            USE ML_FF_STRUCT, ONLY: ML_FF_INPUT, ABINITIO, FF_PAR, ML_IO_WRITE 
            IMPLICIT NONE
! Input and output variables
            TYPE (ML_FF_INPUT)               :: ML_INPUT_HANDLE
            TYPE (ABINITIO)                  :: AB
            TYPE (FF_PAR)                    :: FF
! Local variables
            TYPE (ABINITIO)                  :: AB_TEMP
            INTEGER                          :: HELP_GROUP(1:AB%NCONF)
            INTEGER                          :: HELP_GROUP_CONF(1:AB%NCONF)
            INTEGER                          :: ICHECK
            INTEGER                          :: ICONF
            INTEGER                          :: IFIRST
            INTEGER                          :: IFLAG
            INTEGER                          :: ILAST
            INTEGER                          :: INDEX_CURRENT
            INTEGER                          :: INTYP
            INTEGER                          :: ISYS
            INTEGER                          :: JCONF
            INTEGER                          :: JFLAG
            INTEGER                          :: JNTYP
            INTEGER                          :: KCONF
            INTEGER                          :: KFLAG
            INTEGER                          :: KNTYP
            INTEGER                          :: NUM_STRUCT
            INTEGER                          :: INDEX_HELP(1:AB%NCONF)
            INTEGER                          :: INDEX_HELP_INVERSE(1:AB%NCONF)
            CHARACTER(LEN=40)                :: STRING_HELPER

            DO ICONF=1, AB%NCONF
               INDEX_HELP(ICONF)=ICONF
            ENDDO

            ! Type dependent arrays
            ALLOCATE(AB_TEMP%EATOM(1:AB%MTYP))
            ALLOCATE(AB_TEMP%LNCONF(1:AB%NB_MAX,1:AB%MTYP))
            ALLOCATE(AB_TEMP%LNIONS_FF(1:AB%NB_MAX,1:AB%MTYP))
            ALLOCATE(AB_TEMP%NBASIS(1:AB%MTYP))
            ALLOCATE(AB_TEMP%POMASS(1:AB%MTYP))
            ALLOCATE(AB_TEMP%TYPE(1:AB%MTYP))

            ! Structure dependent arrays
            ALLOCATE(AB_TEMP%A(1:3,1:3,1:1))
            ALLOCATE(AB_TEMP%B(1:3,1:3,1:1))
            ALLOCATE(AB_TEMP%CTIFOR_ALLCONF_OLD(1:1))
            ALLOCATE(AB_TEMP%MITYP_CONF(1:1))
            ALLOCATE(AB_TEMP%NITYP(1:AB%MTYP,1:1))
            ALLOCATE(AB_TEMP%NTYP(1:1))
            ALLOCATE(AB_TEMP%NIONS(1:1))
            ALLOCATE(AB_TEMP%OMEGA(1:1))
            ALLOCATE(AB_TEMP%POSION(1:3,1:AB%MIONS,1:1))
            ALLOCATE(AB_TEMP%SZNAM2(1:1))
            ALLOCATE(AB_TEMP%SZNAM2_INPUT(1:1))
            ALLOCATE(AB_TEMP%TIFOR(1:3,1:AB%MIONS,1:1))
            ALLOCATE(AB_TEMP%TOTEN(1:1))
            ALLOCATE(AB_TEMP%TSIF(1:3,1:3,1:1))
            ALLOCATE(AB_TEMP%TYPE_CONF(1:AB%MTYP,1:1))
            ALLOCATE(AB_TEMP%CTIFOR_PRESENT(1:1))
            
            ! First group structures roughly
            DO ICONF=2, AB%NCONF
               DO JCONF=ICONF-1,1,-1
                  IF (CHECK_TWO_SYSTEMS(JCONF,ICONF,AB,FF)) THEN
                     ! Copy the current structure (ICONF) to helping type
                     CALL COPY_INTO_HELPING_STRUCTURE(AB_TEMP,AB,1,ICONF)
                     ! Increment all structures by 1 from the (JCONF+1)th to the (ICONF-1)th structure
                     DO KCONF=ICONF-1,JCONF+1,-1
                        CALL COPY_INTO_HELPING_STRUCTURE_WITHIN_SAME_OBJECT(AB,(KCONF+1),KCONF)
                     ENDDO
                     ! Copy from the helping type to the (JCONF+1)th position
                     CALL COPY_INTO_HELPING_STRUCTURE(AB,AB_TEMP,(JCONF+1),1)
                     AB%SZNAM2(ICONF)=AB%SZNAM2(JCONF)
                     ! Update structure index list
                     CALL UPDATE_STRUCTURE_INDEX_LIST(INDEX_HELP,ICONF,JCONF)
                     ! Important exit loop after first occurrence of same structure
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
           
            ! Make a grouping before considering names
            ! This is needed to know how many structures are together
            HELP_GROUP_CONF=0
            NUM_STRUCT=1
            HELP_GROUP(1)=NUM_STRUCT
            HELP_GROUP_CONF(NUM_STRUCT)=HELP_GROUP_CONF(NUM_STRUCT)+1
            DO ICONF=2, AB%NCONF
               IF (CHECK_TWO_SYSTEMS(ICONF,ICONF-1,AB,FF)) THEN
                  HELP_GROUP(ICONF)=NUM_STRUCT
                  HELP_GROUP_CONF(NUM_STRUCT)=HELP_GROUP_CONF(NUM_STRUCT)+1
               ELSE
                  NUM_STRUCT=NUM_STRUCT+1
                  HELP_GROUP(ICONF)=NUM_STRUCT
                  HELP_GROUP_CONF(NUM_STRUCT)=HELP_GROUP_CONF(NUM_STRUCT)+1
               ENDIF
            ENDDO
            
            ! Here the structures are grouped by types, number of atoms
            ! and number of atoms per types. Inside the group we want to use
            ! the structure names from the ML_AB file to group them further
            INDEX_CURRENT = 1
            DO ISYS=1, NUM_STRUCT
               IF (HELP_GROUP_CONF(ISYS).GT.1) THEN
                  DO ICONF=(INDEX_CURRENT+1),(INDEX_CURRENT+HELP_GROUP_CONF(ISYS)-1)
                     DO JCONF=ICONF-1,INDEX_CURRENT,-1
                        IF (AB%SZNAM2_INPUT(ICONF).EQ.AB%SZNAM2_INPUT(JCONF)) THEN
                           CALL COPY_INTO_HELPING_STRUCTURE(AB_TEMP,AB,1,ICONF)
                           ! Increment all structures by 1 from the (JCONF+1)th to the (ICONF-1)th structure
                           DO KCONF=ICONF-1,JCONF+1,-1
                              CALL COPY_INTO_HELPING_STRUCTURE_WITHIN_SAME_OBJECT(AB,(KCONF+1),KCONF)
                           ENDDO
                           ! Copy from the helping type to the (JCONF+1)th position
                           CALL COPY_INTO_HELPING_STRUCTURE(AB,AB_TEMP,(JCONF+1),1)
                           AB%SZNAM2(ICONF)=AB%SZNAM2(JCONF)
                           ! Update structure index list
                           CALL UPDATE_STRUCTURE_INDEX_LIST(INDEX_HELP,ICONF,JCONF)
                           ! Important exit loop after first occurrence of same structure
                           EXIT
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
               INDEX_CURRENT=INDEX_CURRENT+HELP_GROUP_CONF(ISYS)
            ENDDO

            ! Give internal names for the structures
            NUM_STRUCT=1
            WRITE(AB%SZNAM2(1),*) NUM_STRUCT
            AB%SZNAM2(1)=TRIM(ADJUSTL(AB%SZNAM2(1)))
            DO ICONF=2, AB%NCONF
               IF (CHECK_TWO_SYSTEMS(ICONF,ICONF-1,AB,FF)) THEN
                  AB%SZNAM2(ICONF)=AB%SZNAM2(ICONF-1)
               ELSE
                  NUM_STRUCT=NUM_STRUCT+1
                  WRITE(AB%SZNAM2(ICONF),*) NUM_STRUCT
               ENDIF
               AB%SZNAM2(ICONF)=TRIM(ADJUSTL(AB%SZNAM2(ICONF)))
            ENDDO

            ICHECK=1
            ! Here we also need to configure the internal name for the new input structure
            DO ICONF=1, AB%NCONF
               IF (CHECK_NEW_STRUCTURE(ICONF)) THEN
                  ML_INPUT_HANDLE%SZNAM2=AB%SZNAM2(ICONF)
                  ICHECK=0
                  EXIT
               ENDIF
            ENDDO
            ! If new structure is not contained within old structure we set the new name here
            IF (ICHECK.EQ.1) THEN
               NUM_STRUCT=NUM_STRUCT+1
               WRITE(ML_INPUT_HANDLE%SZNAM2,*) NUM_STRUCT
            ENDIF

            ! For a continuation run ISTART=1 if the new structure has the same name as any of
            ! the old structures, we need to move the structures to the bottom
            ! and permutate such that the structures with the same name as the new structure
            ! align with the new structure
            IF (FF%ISTART.NE.3 .AND. ICHECK.EQ.0) THEN
               ! First check the first index of the structure with the same internal name as the input structure
               IFIRST=1
               DO ICONF=1,AB%NCONF
                  IF (AB%SZNAM2(ICONF).EQ.ML_INPUT_HANDLE%SZNAM2) THEN
                     IFIRST=MAX(IFIRST,ICONF)
                     EXIT
                  ENDIF
               ENDDO
               ! Second check the last index of the structure with the same internal name as the input structure
               ILAST=AB%NCONF
               DO ICONF=AB%NCONF,1,-1
                  IF (AB%SZNAM2(ICONF).EQ.ML_INPUT_HANDLE%SZNAM2) THEN
                     ILAST=MAX(IFIRST,ICONF)
                     EXIT
                  ENDIF
               ENDDO
               ! Copy everything between these two loop indices to the end of data and possibly update index list
               CALL REARRANGE_SAME_NAMES(AB,AB_TEMP,IFIRST,ILAST,INDEX_HELP)
               ! Determine the first entry of the structures with the same input name as the input structure
               IFIRST=1
               DO ICONF=1,AB%NCONF
                  IF (AB%SZNAM2_INPUT(ICONF).EQ.ML_INPUT_HANDLE%SZNAM2_INPUT) THEN
                     IFIRST=MAX(IFIRST,ICONF)
                     EXIT
                  ENDIF
               ENDDO
               ! Second check the last index of the structure with the same internal name as the input structure
               ILAST=AB%NCONF
               DO ICONF=AB%NCONF,1,-1
                  IF (AB%SZNAM2_INPUT(ICONF).EQ.ML_INPUT_HANDLE%SZNAM2_INPUT) THEN
                     ILAST=MAX(IFIRST,ICONF)
                     EXIT
                  ENDIF
               ENDDO
               ! Copy everything between these two loop indices to the end of data
               CALL REARRANGE_SAME_NAMES(AB,AB_TEMP,IFIRST,ILAST,INDEX_HELP)
            ENDIF
            

            ! Everything related to the order of atom types has to come later
            ! Copy type information to helping arrays
            AB_TEMP%EATOM(1:AB%MTYP)=AB%EATOM(1:AB%MTYP)
            AB_TEMP%LNCONF(1:AB%NB_MAX,1:AB%MTYP)=AB%LNCONF(1:AB%NB_MAX,1:AB%MTYP)
            AB_TEMP%LNIONS_FF(1:AB%NB_MAX,1:AB%MTYP)=AB%LNIONS_FF(1:AB%NB_MAX,1:AB%MTYP)
            AB_TEMP%NBASIS(1:AB%MTYP)=AB%NBASIS(1:AB%MTYP)
            AB_TEMP%POMASS(1:AB%MTYP)=AB%POMASS(1:AB%MTYP)
            AB_TEMP%TYPE(1:AB%MTYP)=AB%TYPE(1:AB%MTYP)

            ! Rearrange the order of atom types
            ! The atom types must occur in the order they appear in the training structures
            DO INTYP=1, AB%MTYP
               IF (INTYP.EQ.1) THEN
                  DO JNTYP=1, AB%MTYP
                     IF (AB_TEMP%TYPE(JNTYP).EQ.AB%TYPE_CONF(1,1)) THEN
                        CALL COPY_TYPE_INFORMATION(AB,AB_TEMP,INTYP,JNTYP)
                     ENDIF
                  ENDDO
               ELSE
                  IFLAG=1
                  DO ICONF=1, AB%NCONF
                     JFLAG=1
                     DO JNTYP=1, AB%NTYP(ICONF)
                        KFLAG=1
                        ! Check all types in the order if they match the first 
                        ! type on the current structure in the loop
                        DO KNTYP=1, INTYP-1
                           IF (AB%TYPE(KNTYP).EQ.AB%TYPE_CONF(JNTYP,ICONF)) THEN
                              KFLAG=0
                              EXIT
                           ENDIF
                        ENDDO
                        IF (KFLAG.EQ.1) THEN
                           DO KNTYP=1, AB%MTYP
                              IF (AB_TEMP%TYPE(KNTYP).EQ.AB%TYPE_CONF(JNTYP,ICONF)) THEN
                                 CALL COPY_TYPE_INFORMATION(AB,AB_TEMP,INTYP,KNTYP)
                              ENDIF
                           ENDDO
                           JFLAG=0
                           IFLAG=0
                           EXIT
                        ENDIF
                        IF (JFLAG.EQ.0) THEN
                           EXIT
                        ENDIF
                     ENDDO
                     IF (IFLAG.EQ.0) THEN
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF 
            ENDDO

            ! This subroutine makes the inverse list of structures for local reference configurations
            CALL MAKE_INVERSE_LIST(AB,INDEX_HELP,INDEX_HELP_INVERSE)
            ! Transform index lists of local reference configurations to new list
            CALL UPDATE_LOCAL_REFERENCE_CONFIGURATION_ENTRIES_ACCORDING_TO_LIST(AB,AB_TEMP,INDEX_HELP_INVERSE)

            DEALLOCATE(AB_TEMP%EATOM)
            DEALLOCATE(AB_TEMP%LNCONF)
            DEALLOCATE(AB_TEMP%LNIONS_FF)
            DEALLOCATE(AB_TEMP%NBASIS)
            DEALLOCATE(AB_TEMP%POMASS)
            DEALLOCATE(AB_TEMP%TYPE)
            DEALLOCATE(AB_TEMP%A)
            DEALLOCATE(AB_TEMP%B)
            DEALLOCATE(AB_TEMP%CTIFOR_ALLCONF_OLD)
            DEALLOCATE(AB_TEMP%MITYP_CONF)
            DEALLOCATE(AB_TEMP%NITYP)
            DEALLOCATE(AB_TEMP%NTYP)
            DEALLOCATE(AB_TEMP%NIONS)
            DEALLOCATE(AB_TEMP%OMEGA)
            DEALLOCATE(AB_TEMP%POSION)
            DEALLOCATE(AB_TEMP%SZNAM2)
            DEALLOCATE(AB_TEMP%SZNAM2_INPUT)
            DEALLOCATE(AB_TEMP%TIFOR)
            DEALLOCATE(AB_TEMP%TOTEN)
            DEALLOCATE(AB_TEMP%TSIF)
            DEALLOCATE(AB_TEMP%TYPE_CONF)
            DEALLOCATE(AB_TEMP%CTIFOR_PRESENT)

            CONTAINS

               ! The structure index of local reference configuration due to possible permutations is updated here
               SUBROUTINE UPDATE_LOCAL_REFERENCE_CONFIGURATION_ENTRIES_ACCORDING_TO_LIST(AB,AB_TEMP,INDEX_HELP_INVERSE)
                  IMPLICIT NONE
                  TYPE (ABINITIO)                  :: AB
                  TYPE (ABINITIO)                  :: AB_TEMP
                  INTEGER, INTENT(IN)              :: INDEX_HELP_INVERSE(:)
                  INTEGER                          :: IB, IB_NEW, INTYP0
                  ! Copy local reference information of structure number into helping array
                  AB_TEMP%LNCONF(1:AB%NB_MAX,1:AB%MTYP)=AB%LNCONF(1:AB%NB_MAX,1:AB%MTYP)
                  ! Loop over all structure numbers for local reference configuration and update them according
                  ! to the new list
                  DO INTYP0=1, AB%MTYP
                     DO IB=1, AB%NBASIS(INTYP0)
                        AB%LNCONF(IB,INTYP0)=INDEX_HELP_INVERSE(AB_TEMP%LNCONF(IB,INTYP0))
                     ENDDO
                  ENDDO
               END SUBROUTINE UPDATE_LOCAL_REFERENCE_CONFIGURATION_ENTRIES_ACCORDING_TO_LIST
             
               ! This subroutine does forward rotation of structure indices from JCONF+1 to ICONF
               SUBROUTINE UPDATE_STRUCTURE_INDEX_LIST(INDEX_HELP,ICONF,JCONF)
                  IMPLICIT NONE
                  INTEGER, INTENT(INOUT)           :: INDEX_HELP(:)
                  INTEGER                          :: ICONF, JCONF, KCONF
                  INTEGER                          :: IHELP
                  IHELP=INDEX_HELP(ICONF)
                  DO KCONF=ICONF-1,JCONF+1,-1
                     INDEX_HELP(KCONF+1)=INDEX_HELP(KCONF)
                  ENDDO
                  INDEX_HELP(JCONF+1)=IHELP
               END SUBROUTINE UPDATE_STRUCTURE_INDEX_LIST 

               ! This subroutine makes the inverse list of structures for local reference configurations
               ! so that e.g. 10th structure in the old can be mapped to new index via INDEX_HELP_INVERSE(10)
               SUBROUTINE MAKE_INVERSE_LIST(AB,INDEX_HELP,INDEX_HELP_INVERSE)
                  IMPLICIT NONE
                  TYPE (ABINITIO)                  :: AB
                  INTEGER, INTENT(IN)              :: INDEX_HELP(:)
                  INTEGER, INTENT(OUT)             :: INDEX_HELP_INVERSE(:)
                  INTEGER                          :: I, J
                  DO I=1, AB%NCONF
                     DO J=1, AB%NCONF
                        IF (I.EQ.INDEX_HELP(J)) THEN
                           INDEX_HELP_INVERSE(I)=J
                           EXIT
                        ENDIF
                     ENDDO
                  ENDDO
               END SUBROUTINE MAKE_INVERSE_LIST 

               ! This subroutine copies type information from OBJ2(I2) to OBJ1(I1)
               SUBROUTINE COPY_TYPE_INFORMATION(OBJ1,OBJ2,I1,I2)
                  IMPLICIT NONE
                  TYPE (ABINITIO)                  :: OBJ1
                  TYPE (ABINITIO)                  :: OBJ2
                  INTEGER                          :: I1
                  INTEGER                          :: I2
                  OBJ1%EATOM(I1)=OBJ2%EATOM(I2)
                  OBJ1%TYPE(I1)=OBJ2%TYPE(I2)
                  OBJ1%POMASS(I1)=OBJ2%POMASS(I2)
                  OBJ1%NBASIS(I1)=OBJ2%NBASIS(I2)
                  OBJ1%LNCONF(1:AB%NB_MAX,I1)=OBJ2%LNCONF(1:AB%NB_MAX,I2)
                  OBJ1%LNIONS_FF(1:AB%NB_MAX,I1)=OBJ2%LNIONS_FF(1:AB%NB_MAX,I2)
               END SUBROUTINE COPY_TYPE_INFORMATION

               ! This subroutine copies configurational information from OBJ2(I2) to OBJ1(I1)
               SUBROUTINE COPY_INTO_HELPING_STRUCTURE(OBJ1,OBJ2,I1,I2)
                  IMPLICIT NONE
                  TYPE (ABINITIO)                  :: OBJ1
                  TYPE (ABINITIO)                  :: OBJ2
                  INTEGER                          :: I1
                  INTEGER                          :: I2        
                  OBJ1%A(1:3,1:3,I1)=OBJ2%A(1:3,1:3,I2)
                  OBJ1%B(1:3,1:3,I1)=OBJ2%B(1:3,1:3,I2)
                  OBJ1%CTIFOR_ALLCONF_OLD(I1)=OBJ2%CTIFOR_ALLCONF_OLD(I2)
                  OBJ1%MITYP_CONF(I1)=OBJ2%MITYP_CONF(I2)
                  OBJ1%NITYP(1:AB%MTYP,I1)=OBJ2%NITYP(1:AB%MTYP,I2)
                  OBJ1%NTYP(I1)=OBJ2%NTYP(I2)
                  OBJ1%NIONS(I1)=OBJ2%NIONS(I2)
                  OBJ1%OMEGA(I1)=OBJ2%OMEGA(I2)
                  OBJ1%POSION(1:3,1:AB%MIONS,I1)=OBJ2%POSION(1:3,1:AB%MIONS,I2)
                  OBJ1%SZNAM2(I1)=OBJ2%SZNAM2(I2)
                  OBJ1%SZNAM2_INPUT(I1)=OBJ2%SZNAM2_INPUT(I2)
                  OBJ1%TIFOR(1:3,1:AB%MIONS,I1)=OBJ2%TIFOR(1:3,1:AB%MIONS,I2)
                  OBJ1%TOTEN(I1)=OBJ2%TOTEN(I2)
                  OBJ1%TSIF(1:3,1:3,I1)=OBJ2%TSIF(1:3,1:3,I2)
                  OBJ1%TYPE_CONF(1:AB%MTYP,I1)=OBJ2%TYPE_CONF(1:AB%MTYP,I2)
                  OBJ1%CTIFOR_PRESENT(I1)=OBJ2%CTIFOR_PRESENT(I2)
               END SUBROUTINE COPY_INTO_HELPING_STRUCTURE

               ! This subroutine copies configurational information from OBJ1(I2) to OBJ1(I1)
               ! Mainly because of NEC compiler
               SUBROUTINE COPY_INTO_HELPING_STRUCTURE_WITHIN_SAME_OBJECT(OBJ1,I1,I2)
                  IMPLICIT NONE
                  TYPE (ABINITIO)                  :: OBJ1
                  INTEGER                          :: I1
                  INTEGER                          :: I2
                  OBJ1%A(1:3,1:3,I1)=OBJ1%A(1:3,1:3,I2)
                  OBJ1%B(1:3,1:3,I1)=OBJ1%B(1:3,1:3,I2)
                  OBJ1%CTIFOR_ALLCONF_OLD(I1)=OBJ1%CTIFOR_ALLCONF_OLD(I2)
                  OBJ1%MITYP_CONF(I1)=OBJ1%MITYP_CONF(I2)
                  OBJ1%NITYP(1:AB%MTYP,I1)=OBJ1%NITYP(1:AB%MTYP,I2)
                  OBJ1%NTYP(I1)=OBJ1%NTYP(I2)
                  OBJ1%NIONS(I1)=OBJ1%NIONS(I2)
                  OBJ1%OMEGA(I1)=OBJ1%OMEGA(I2)
                  OBJ1%POSION(1:3,1:AB%MIONS,I1)=OBJ1%POSION(1:3,1:AB%MIONS,I2)
                  OBJ1%SZNAM2(I1)=OBJ1%SZNAM2(I2)
                  OBJ1%SZNAM2_INPUT(I1)=OBJ1%SZNAM2_INPUT(I2)
                  OBJ1%TIFOR(1:3,1:AB%MIONS,I1)=OBJ1%TIFOR(1:3,1:AB%MIONS,I2)
                  OBJ1%TOTEN(I1)=OBJ1%TOTEN(I2)
                  OBJ1%TSIF(1:3,1:3,I1)=OBJ1%TSIF(1:3,1:3,I2)
                  OBJ1%TYPE_CONF(1:AB%MTYP,I1)=OBJ1%TYPE_CONF(1:AB%MTYP,I2)
                  OBJ1%CTIFOR_PRESENT(I1)=OBJ1%CTIFOR_PRESENT(I2)
               END SUBROUTINE COPY_INTO_HELPING_STRUCTURE_WITHIN_SAME_OBJECT

               ! Sets the internal name of the new input structure
               FUNCTION CHECK_NEW_STRUCTURE(ICONF)
                  IMPLICIT NONE
                  LOGICAL                          :: CHECK_NEW_STRUCTURE
                  INTEGER                          :: ICONF
                  INTEGER                          :: IFLAG
                
                  ! Check number of atom types
                  IF (ML_INPUT_HANDLE%NTYP.NE.AB%NTYP(ICONF)) THEN
                     CHECK_NEW_STRUCTURE=.FALSE.
                     RETURN
                  ENDIF
                  ! Check number of atoms 
                  IF (ML_INPUT_HANDLE%NIONS.NE.AB%NIONS(ICONF)) THEN
                     CHECK_NEW_STRUCTURE=.FALSE.
                     RETURN
                  ENDIF
                  ! Check individual types within structures
                  DO INTYP=1, ML_INPUT_HANDLE%NTYP
                     IFLAG=1
                     ! Check type strings
                     DO JNTYP=1, ML_INPUT_HANDLE%NTYP
                        IF (ML_INPUT_HANDLE%TYPE(INTYP).EQ.AB%TYPE_CONF(JNTYP,ICONF)) THEN
                           IFLAG=0
                           EXIT
                        ENDIF
                     ENDDO
                     ! If same types are not present, then different structure, 
                     ! else check if number of atoms of each species is the same
                     IF (IFLAG.EQ.1) THEN
                        CHECK_NEW_STRUCTURE=.FALSE.
                        RETURN
                     ELSE
                        ! The loop counter is taken here due to the last exit
                        IF (ML_INPUT_HANDLE%NITYP(INTYP).NE.AB%NITYP(JNTYP,ICONF)) THEN
                           CHECK_NEW_STRUCTURE=.FALSE.
                           RETURN
                        ENDIF
                     ENDIF
                  ENDDO
                  ! Check structure names
                  IF (FF%LUSE_NAMES) THEN
                     IF (ML_INPUT_HANDLE%SZNAM2_INPUT.NE.AB%SZNAM2_INPUT(ICONF)) THEN
                        CHECK_NEW_STRUCTURE=.FALSE.
                        RETURN
                     ENDIF
                  ENDIF
                  ! Otherwise both systems are the same
                  CHECK_NEW_STRUCTURE=.TRUE.
                  RETURN
               END FUNCTION CHECK_NEW_STRUCTURE

               SUBROUTINE REARRANGE_SAME_NAMES(AB,AB_TEMP,IFIRST,ILAST,INDEX_HELP)
                  IMPLICIT NONE
                  TYPE (ABINITIO)                  :: AB
                  TYPE (ABINITIO)                  :: AB_TEMP
                  INTEGER                          :: IFIRST
                  INTEGER                          :: IHELP
                  INTEGER                          :: ILAST
                  INTEGER                          :: ICONF
                  INTEGER                          :: JCONF
                  INTEGER                          :: INDEX_HELP(:)
                  ! Copy everything between these two loops to the end
                  DO ICONF=IFIRST,ILAST
                     ! Copy at first index to helping type
                     CALL COPY_INTO_HELPING_STRUCTURE(AB_TEMP,AB,1,IFIRST)
                     ! Move each entry back by one postion
                     DO JCONF=IFIRST+1,AB%NCONF
                        ! Copy each position JCONF into the previous position JCONF-1
                        CALL COPY_INTO_HELPING_STRUCTURE_WITHIN_SAME_OBJECT(AB,(JCONF-1),JCONF)
                     ENDDO
                     ! Copy data from helping structure into last postion
                     CALL COPY_INTO_HELPING_STRUCTURE(AB,AB_TEMP,AB%NCONF,1)
                  ENDDO
                  ! Update helping index list backwards here (in principle the same operation as on the type variables)
                  DO ICONF=IFIRST, ILAST
                     IHELP=INDEX_HELP(IFIRST)
                     DO JCONF=IFIRST+1,AB%NCONF
                        INDEX_HELP(JCONF-1)=INDEX_HELP(JCONF)
                     ENDDO
                     INDEX_HELP(AB%NCONF)=IHELP
                  ENDDO
               END SUBROUTINE REARRANGE_SAME_NAMES
            ! End CONTAINS 
         END SUBROUTINE REARRANGE_ML_AB

!****************************************************************************************************
! Reset ML_MB if it is too small for ML_AB file
!****************************************************************************************************
        SUBROUTINE RESET_ML_MB(AB,FFM)
           USE ML_FF_STRUCT, ONLY: ABINITIO, FFM_PAR, ML_IO_WRITE 
           IMPLICIT NONE
! Input and output variables
           TYPE (ABINITIO)                  :: AB
           TYPE (FFM_PAR)                   :: FFM
           INTEGER                          :: ITYPE, MAX_MB_AB
           ! Set the target value for ML_MB here
           MAX_MB_AB=AB%NBASIS(1)
           DO ITYPE=1,AB%MTYP
              IF (AB%NBASIS(ITYPE).GT.MAX_MB_AB) THEN
                 MAX_MB_AB=AB%NBASIS(ITYPE)
              ENDIF
           ENDDO
           IF (MAX_MB_AB.GT.FFM%MB) THEN
              FFM%MB=MAX_MB_AB
           ENDIF
        END SUBROUTINE RESET_ML_MB


!****************************************************************************************************
! Read element number written in ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_ELEMENT_ORDER1 (AB_HELP,MTYP_OLD)
          USE ML_FF_STRUCT, ONLY: ABINITIO
          IMPLICIT NONE
! Input and ouput variables.
          TYPE (ABINITIO)                 :: AB_HELP
          INTEGER         , INTENT(INOUT) :: MTYP_OLD
          MTYP_OLD=AB_HELP%MTYP
        END SUBROUTINE READ_ABCAR_ELEMENT_ORDER1

!****************************************************************************************************
! Read element data written in ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_ELEMENT_ORDER2 (AB_HELP,MB_INPUT,MCONF_INPUT,EATOM_OLD,MIONS_OLD,MITYP_OLD,MTYP_OLD,POMASS_OLD,TYPE_OLD)
          USE ML_FF_STRUCT, ONLY: ABINITIO
          USE ml_ff_tutor, ONLY: ml_tutor
          USE ml_ff_string, ONLY: str
          IMPLICIT NONE
! Input variables.
          TYPE (ABINITIO)                 :: AB_HELP
          INTEGER         , INTENT(IN)    :: MB_INPUT
          INTEGER         , INTENT(IN)    :: MCONF_INPUT
          INTEGER         , INTENT(IN)    :: MTYP_OLD
! Input and ouput variables.
          REAL(q)         , INTENT(INOUT) :: EATOM_OLD(:) !(1:MTYP_OLD)
          INTEGER         , INTENT(INOUT) :: MIONS_OLD
          INTEGER         , INTENT(INOUT) :: MITYP_OLD
          REAL(q)         , INTENT(INOUT) :: POMASS_OLD(:) !(1:MTYP_OLD)
          CHARACTER(LEN=2), INTENT(INOUT) :: TYPE_OLD(:) !(1:MTYP_OLD)
! Local variables
          INTEGER                         :: I
          INTEGER                         :: MCONF_AB
          INTEGER                         :: NUM_BASIS_OLD(1:MTYP_OLD)
          MCONF_AB=AB_HELP%MCONF
          IF (MCONF_AB.GT.MCONF_INPUT) THEN
             CALL ml_tutor%error("ML_MCONF is smaller than the number &
                &of structures in the ML_AB file (" // str(MCONF_AB) // &
                "), please increase ML_MCONF.")
          ENDIF
          DO I=1, MTYP_OLD
             TYPE_OLD(I)=AB_HELP%TYPE(I)
          ENDDO
          MIONS_OLD=AB_HELP%MIONS
          MITYP_OLD=AB_HELP%MITYP
          DO I=1, MTYP_OLD
             EATOM_OLD(I)=AB_HELP%EATOM(I)
             POMASS_OLD(I)=AB_HELP%POMASS(I)
             NUM_BASIS_OLD(I)=AB_HELP%NBASIS(I)
          ENDDO
          DO I=1, MTYP_OLD
             IF (NUM_BASIS_OLD(I).GT.MB_INPUT) THEN
                CALL ml_tutor%error("ML_MB is smaller than the number &
                   &of local reference configurations in the ML_AB file &
                   &(atom type " // str(I) // ": " // str(NUM_BASIS_OLD(I)) // &
                   "), please increase ML_MB.")
             ENDIF
          ENDDO
! Convert units already immediately after reading in
          EATOM_OLD(1:MTYP_OLD)=EATOM_OLD(1:MTYP_OLD)/EUNIT
          POMASS_OLD(1:MTYP_OLD)=POMASS_OLD(1:MTYP_OLD)/MUNIT
        END SUBROUTINE READ_ABCAR_ELEMENT_ORDER2

!****************************************************************************************************
! Get the number of configurations in the ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_SYSTEM_ORDER1 (NCONF,AB_HELP)
           USE ML_FF_STRUCT, ONLY: ABINITIO
           IMPLICIT NONE
! Input and output variables.
           INTEGER         , INTENT(INOUT) :: NCONF
           TYPE (ABINITIO)                 :: AB_HELP
           NCONF = AB_HELP%NCONF           
        END SUBROUTINE READ_ABCAR_SYSTEM_ORDER1

!****************************************************************************************************
! Determine the number of different systems in the ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_SYSTEM_ORDER2 (NCONF,NSYS,SZNAM2_HELP2,SZNAM2_HELP2_INPUT,AB_HELP)
          USE ML_FF_STRUCT, ONLY: ABINITIO
          IMPLICIT NONE
! Input variables.
          INTEGER          , INTENT(IN)    :: NCONF
          TYPE (ABINITIO)                  :: AB_HELP
! Input and output variables.
          INTEGER          , INTENT(INOUT) :: NSYS
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_HELP2(:) !(1:NCONF+1)
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_HELP2_INPUT(:) !(1:NCONF+1)
! Local variables.
          INTEGER                          :: ICONF
          INTEGER                          :: IFLAG
          INTEGER                          :: INTYP
          INTEGER                          :: ISYS
          INTEGER                          :: MTYP
          INTEGER          , ALLOCATABLE   :: NB(:)
          CHARACTER(LEN=40)                :: SZNAM2_HELP1
! Get the system data written in the ABCAR file.
          NSYS=0
          MTYP=AB_HELP%MTYP
          ALLOCATE(NB(1:MTYP))
          DO INTYP=1, MTYP
             NB(INTYP)=AB_HELP%NBASIS(INTYP)
          ENDDO
          DO ICONF=1, NCONF
             SZNAM2_HELP1=AB_HELP%SZNAM2(ICONF)
             SZNAM2_HELP1=TRIM(ADJUSTL(SZNAM2_HELP1))
             IF(NSYS.EQ.0) THEN
                NSYS=NSYS+1
                SZNAM2_HELP2(NSYS)=SZNAM2_HELP1
                SZNAM2_HELP2_INPUT(NSYS)=AB_HELP%SZNAM2_INPUT(ICONF)
             ELSE
                IFLAG=1
                DO ISYS=1, NSYS
                   IF(SZNAM2_HELP2(ISYS).EQ.SZNAM2_HELP1) THEN
                      IFLAG=0
                      EXIT
                   ENDIF
                ENDDO
                IF(IFLAG.EQ.1) THEN
                   NSYS=NSYS+1
                   SZNAM2_HELP2(NSYS)=SZNAM2_HELP1
                   SZNAM2_HELP2_INPUT(NSYS)=AB_HELP%SZNAM2_INPUT(ICONF)
                ENDIF
             ENDIF
          ENDDO
! Deallocate array.
          DEALLOCATE(NB)
        END SUbROUTINE READ_ABCAR_SYSTEM_ORDER2

!****************************************************************************************************
! Read size written in ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_SIZE1 (AB_HELP,LMLMB_FF,MDATA_FFM,MTYP,NCONF)
          USE ML_FF_STRUCT, ONLY: ABINITIO
          IMPLICIT NONE
! Input variables.
          TYPE (ABINITIO)                 :: AB_HELP
          LOGICAL         , INTENT(IN)    :: LMLMB_FF
          INTEGER         , INTENT(IN)    :: MTYP
! In put and output variables.
          INTEGER         , INTENT(INOUT) :: MDATA_FFM
          INTEGER         , INTENT(INOUT) :: NCONF
! Local variables.
          INTEGER                         :: ICONF
          MDATA_FFM=0
          NCONF=AB_HELP%MCONF
          DO ICONF=1, NCONF
             IF (LMLMB_FF) THEN
                MDATA_FFM=MDATA_FFM+3*AB_HELP%NIONS(ICONF)+7
             ENDIF
          ENDDO
        END SUBROUTINE READ_ABCAR_SIZE1

!****************************************************************************************************
! Read size data in ABCAR file
! CTIFOR_FF, IERR_FF are overwritten.
!****************************************************************************************************
        SUBROUTINE READ_ABCAR_SIZE2 (AB_HELP,ML_INPUT_HANDLE,AB,FF,MITYP_NEW)
          USE ML_FF_STRUCT, ONLY : ML_FF_INPUT,ABINITIO,FF_PAR
          IMPLICIT NONE
! Input and output variables
          TYPE (ABINITIO)                 :: AB_HELP
          TYPE (ML_FF_INPUT)              :: ML_INPUT_HANDLE
          TYPE (ABINITIO)                 :: AB
          TYPE (FF_PAR)                   :: FF
! Input variables.
          INTEGER         , INTENT(IN)    :: MITYP_NEW
! Local variables.
          REAL(q),ALLOCATABLE             :: EATOM_TMP(:)
          INTEGER                         :: I
          INTEGER                         :: ICONF
          INTEGER                         :: IFLAG
          INTEGER                         :: INTYP
          INTEGER                         :: ISYS
          INTEGER                         :: JNTYP
          INTEGER                         :: MIONS_TMP
          INTEGER                         :: MITYP_TMP
          INTEGER                         :: MTYP_TMP
          INTEGER,ALLOCATABLE             :: NB_TMP(:)
          INTEGER,ALLOCATABLE             :: NITYP_TMP2(:)
          REAL(q),ALLOCATABLE             :: POMASS_TMP(:)
          CHARACTER(LEN=40)               :: SZNAM2_HELP1
          CHARACTER(LEN=2),ALLOCATABLE    :: TYPE_TMP(:)
! Set constants
          AB%NSYS=0
          AB%NCONF=AB_HELP%NCONF
          AB%MTYP=AB_HELP%MTYP
! Allocate helping arrays
          ALLOCATE(TYPE_TMP(AB%MTYP))
          ALLOCATE(EATOM_TMP(AB%MTYP))
          ALLOCATE(NB_TMP(AB%MTYP))
          ALLOCATE(NITYP_TMP2(AB%MTYP))
          ALLOCATE(POMASS_TMP(AB%MTYP))
          DO I=1, AB%MTYP
             TYPE_TMP(I)=AB_HELP%TYPE(I)
          ENDDO 
          MIONS_TMP=AB_HELP%MIONS
          MITYP_TMP=AB_HELP%MITYP
          DO I=1, AB%MTYP
             EATOM_TMP(I)=AB_HELP%EATOM(I)
             POMASS_TMP(I)=AB_HELP%POMASS(I)
             NB_TMP(I)=AB_HELP%NBASIS(I)
          ENDDO
! Define sizes
          AB%NSYS=FF%NSYS
          AB%MCONF=FF%MCONF
          IF (FF%ISTART.NE.3) THEN
             AB%MIONS=MAX(MIONS_TMP,ML_INPUT_HANDLE%NIONS)
             AB%MITYP=MAX(MITYP_TMP,MITYP_NEW)
! Check if current systems types are among all types
             MTYP_TMP=AB%MTYP
             DO INTYP=1, ML_INPUT_HANDLE%NTYP
                IFLAG=1
                DO JNTYP=1, AB%MTYP
                   IF (TYPE_TMP(JNTYP).EQ.ML_INPUT_HANDLE%TYPE(INTYP)) THEN
                      IFLAG=0
                      EXIT
                   ENDIF
                ENDDO
                IF (IFLAG.EQ.1) THEN
                   MTYP_TMP=MTYP_TMP+1
                ENDIF
             ENDDO
             AB%MTYP=MTYP_TMP
! Obtain maximum total number of ions in systemfrom current structure and ML_AB file
             AB%MIONS=MAX(AB%MIONS,ML_INPUT_HANDLE%NIONS)
! Obtain maximum number of atoms per atom type
             DO INTYP=1, ML_INPUT_HANDLE%NTYP
                MITYP_TMP=MAX(MITYP_TMP,ML_INPUT_HANDLE%NITYP(INTYP))
             ENDDO
             AB%MITYP=MITYP_TMP
          ELSE
             AB%MIONS=MIONS_TMP
             AB%MITYP=MITYP_TMP
          ENDIF
! Deallocate helping arrays
          DEALLOCATE(TYPE_TMP)
          DEALLOCATE(EATOM_TMP)
          DEALLOCATE(NB_TMP)
          DEALLOCATE(NITYP_TMP2)
          DEALLOCATE(POMASS_TMP)
        END SUBROUTINE READ_ABCAR_SIZE2

        SUBROUTINE READ_ABCAR_RESET_SIZES(AB_HELP,AB,FF,FFM,ML_INPUT_HANDLE,TAG_LIST)
          USE ML_FF_STRUCT, ONLY: ABINITIO, FF_PAR, FFM_PAR, ML_FF_INPUT, ML_MB_DEFAULT
          USE ML_FF_TAGLIST
          IMPLICIT NONE
! Input and output variables
          TYPE (ABINITIO)               :: AB_HELP
          TYPE (ABINITIO)               :: AB
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
          TYPE (ML_FF_INPUT)            :: ML_INPUT_HANDLE
          TYPE (TAGLIST), INTENT(INOUT) :: TAG_LIST
! Local variables
          INTEGER              :: I
          INTEGER              :: MCONF_HELP
          INTEGER              :: MTYP_TMP
          INTEGER, ALLOCATABLE :: NB_TMP(:) 
          INTEGER              :: NB_LARGEST_TMP
          MCONF_HELP=AB_HELP%MCONF
          MTYP_TMP=AB_HELP%MTYP
          ! Do allocation of helping array for local configuratins
          ALLOCATE(NB_TMP(1:MTYP_TMP))
          DO I=1,MTYP_TMP
             NB_TMP(I)=AB_HELP%NBASIS(I)
          ENDDO
          ! If ML_MB was not set in the INCAR file then 
          ! change the default here
          IF (TAG_LIST%GET_STATE("ML_MB") == TL_DEFAULT) THEN
             CALL TAG_LIST%SET_STATE("ML_MB", TL_OVERRIDE)
             IF (FF%ISTART.EQ.3) THEN
                FFM%MB=MIN(AB_HELP%MCONF*AB_HELP%MITYP,ML_MB_DEFAULT)
             ELSE
                NB_LARGEST_TMP=NB_TMP(1)
                DO I=1, MTYP_TMP
                   NB_LARGEST_TMP=MAX(NB_LARGEST_TMP,NB_TMP(I))
                ENDDO
                FFM%MB = NB_LARGEST_TMP!+MIN(MAX(FF%NSW,1),FF%MCONF_NEW)*AB_HELP%MITYP
             ENDIF
             FFM%MB_ALLOCATE=FFM%MB
             FFM%MB_INPUT=FFM%MB
          ENDIF
          ! If ML_MCONF was not set in the INCAR file then
          ! change the default here if needed
          IF (TAG_LIST%GET_STATE("ML_MCONF") == TL_DEFAULT) THEN
             CALL TAG_LIST%SET_STATE("ML_MCONF", TL_OVERRIDE)
             FF%MCONF=MCONF_HELP!+500+FF%MCONF_NEW
          ENDIF
          DEALLOCATE(NB_TMP)
        END SUBROUTINE READ_ABCAR_RESET_SIZES
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Read data in ABCAR file
!****************************************************************************************************

        SUBROUTINE READ_ABCAR (AB_HELP,AB,FF,A_AB,B_AB,EATOM,ISTART_FF,ISYS_AB, &
                   LADD_AB,LB_AB,LB_FFM,LBHEAD_FFM, &
                   LNCONF_FF,LNIONS_AB,LNIONS_FF,LNITYP_AB,LNTYP_AB,LSYS_AB, &
                   MTYP_AB,MTYP_FF, &
                   NB_FFM,NBMAX_TMP,NBTOT_FFM,NCONF_AB,NCONF_SYS_AB,NIONS_AB,NITYP_AB,NITYP_ABN, &
                   NSYS_AB,NSYS_FF,NTYP_AB,NTYP_NEW, &
                   OMEGA_AB,POMASS_AB,POMASS_FF,POSION_AB,SZNAM2_AB,SZNAM2_ADD_AB,SZNAM2_FF, &
                   TIFOR_AB,TOTEN_AB,TSIF_AB,TYPE_AB,TYPE_FF)
          USE ML_FF_STRUCT, ONLY: ABINITIO, FF_PAR
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
! Input variables.
          TYPE (ABINITIO)                  :: AB_HELP
          TYPE (ABINITIO)                  :: AB
          TYPE (FF_PAR)                    :: FF
          INTEGER          , INTENT(IN)    :: ISTART_FF
          INTEGER          , INTENT(IN)    :: NCONF_AB
          INTEGER          , INTENT(IN)    :: MTYP_AB
          INTEGER          , INTENT(IN)    :: MTYP_FF
          INTEGER          , INTENT(IN)    :: NITYP_ABN(:) !(1:MTYP_ABN)
          INTEGER          , INTENT(IN)    :: NSYS_AB
          INTEGER          , INTENT(IN)    :: NSYS_FF
          INTEGER          , INTENT(IN)    :: NTYP_NEW
          REAL(q)          , INTENT(IN)    :: POMASS_FF(:) !(1:MTYP_FF)
          CHARACTER(LEN=40), INTENT(IN)    :: SZNAM2_FF(:) !(1:NSYS_FF)
          CHARACTER(LEN=2) , INTENT(IN)    :: TYPE_FF(:) !(1:MTYP_FF)
! Input and output variables.
          REAL(q)          , INTENT(INOUT) :: A_AB(:,:,:) !(1:3,1:3,1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: B_AB(:,:,:) !(1:3,1:3,1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: EATOM(:)
          INTEGER          , INTENT(INOUT) :: ISYS_AB(:) !(1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LADD_AB(:,:) !(1:MTYP_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LB_AB(:,:) !(1:MIONS_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LB_FFM(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LBHEAD_FFM(:) !(1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LNCONF_FF(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LNIONS_FF(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LNIONS_AB(:,:,:) !(1:MITYP_AB,1:MTYP_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LNITYP_AB(:,:) !(1:MIONS_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LNTYP_AB(:,:) !(1:MIONS_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LSYS_AB(:) !(1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: NB_FFM(:) !(1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: NBMAX_TMP(:) !(1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: NBTOT_FFM
          INTEGER          , INTENT(INOUT) :: NCONF_SYS_AB(:) !(1:NSYS_AB)
          INTEGER          , INTENT(INOUT) :: NIONS_AB(:) !(1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: NITYP_AB(:,:) !(1:MTYP_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: NTYP_AB(:) !(1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: OMEGA_AB(:) !(1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: POMASS_AB(:) !(1:MTYP_AB)
          REAL(q)          , INTENT(INOUT) :: POSION_AB(:,:,:) !(1:3,1:MIONS_AB,1:MCONF_AB)
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_AB(:) !(1:MCONF_AB)
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_ADD_AB(:) !(1:NSYS_AB)
          REAL(q)          , INTENT(INOUT) :: TIFOR_AB(:,:,:) !(1:3,1:MIONS_AB,1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: TOTEN_AB(:) !(1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: TSIF_AB(:,:,:) !(1:3,1:3,1:MCONF_AB)
          CHARACTER(LEN=2) , INTENT(INOUT) :: TYPE_AB(:) !(1:MTYP_AB)
! Local variables.
          REAL(q)                          :: EATOM_TMP(1:MTYP_FF)
          INTEGER                          :: I
          INTEGER                          :: IB
          INTEGER                          :: ICONF
          INTEGER                          :: INIONS
          INTEGER                          :: INITYP
          INTEGER                          :: INTYP
          INTEGER                          :: INTYP0
          INTEGER                          :: ISYS
          INTEGER                          :: IXYZ
          INTEGER                          :: JNTYP0
          INTEGER                          :: JXYZ
          INTEGER                          :: MIONS_TMP
          INTEGER                          :: MITYP_TMP
          INTEGER                          :: MTYP_TMP
          INTEGER                          :: NB_TMP(1:MTYP_FF)
          INTEGER                          :: NITYP_TMP2(1:MTYP_FF)
          REAL(q)                          :: POMASS_TMP(1:MTYP_FF)
          CHARACTER(LEN=2)                 :: TYPE_TMP(1:MTYP_FF)
          EATOM_TMP = 0.0_q
! Initialize logical parameters
          NBMAX_TMP=0
          MTYP_TMP=AB_HELP%MTYP
          DO I=1, AB_HELP%MTYP
             TYPE_TMP(I)=AB_HELP%TYPE(I)
          ENDDO
          MIONS_TMP=AB_HELP%MIONS
          MITYP_TMP=AB_HELP%MITYP
          DO I=1, MTYP_TMP
             EATOM_TMP(I)=AB_HELP%EATOM(I)
             POMASS_TMP(I)=AB_HELP%POMASS(I)
          ENDDO
! Configure the model system name.
          DO ISYS=1, NSYS_FF
             SZNAM2_ADD_AB(ISYS)=SZNAM2_FF(ISYS)
          ENDDO
! Set TYPE and POMASS.
          DO INTYP0=1, MTYP_AB
             POMASS_AB(INTYP0)=POMASS_FF(INTYP0)
             TYPE_AB(INTYP0)=TYPE_FF(INTYP0)
          ENDDO
! Configure basis sets
          DO I=1, MTYP_TMP
             NB_TMP(I)=AB_HELP%NBASIS(I)
          ENDDO
          NBTOT_FFM=0
          DO INTYP=1, MTYP_TMP
             DO JNTYP0=1, MTYP_FF
                IF(TYPE_TMP(INTYP).EQ.TYPE_FF(JNTYP0)) THEN
                   NB_FFM(JNTYP0)=NB_TMP(INTYP)
                   NBTOT_FFM=NBTOT_FFM+NB_TMP(INTYP)
                   DO IB=1, NB_TMP(INTYP)
                      LNCONF_FF(IB,JNTYP0)=AB_HELP%LNCONF(IB,INTYP)
                      LNIONS_FF(IB,JNTYP0)=AB_HELP%LNIONS_FF(IB,INTYP)
                      LB_AB(LNIONS_FF(IB,JNTYP0),LNCONF_FF(IB,JNTYP0))=IB
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          LBHEAD_FFM(1)=0
          DO INTYP0=2, MTYP_FF
             LBHEAD_FFM(INTYP0)=LBHEAD_FFM(INTYP0-1)+NB_FFM(INTYP0-1)
          ENDDO
          NBTOT_FFM=0
          DO INTYP0=1, MTYP_FF
             DO IB=1, NB_FFM(INTYP0)
                NBTOT_FFM=NBTOT_FFM+1
                LB_FFM(IB,INTYP0)=NBTOT_FFM
             ENDDO
          ENDDO
          ! Loop over configurations
          DO ICONF=1, NCONF_AB
             SZNAM2_AB(ICONF)=AB_HELP%SZNAM2(ICONF)
             SZNAM2_AB(ICONF)=TRIM(ADJUSTL(SZNAM2_AB(ICONF)))
             AB%SZNAM2_INPUT(ICONF)=TRIM(ADJUSTL(AB_HELP%SZNAM2_INPUT(ICONF)))
             DO ISYS=1, NSYS_AB
                IF(SZNAM2_AB(ICONF).EQ.SZNAM2_ADD_AB(ISYS)) THEN
                   ISYS_AB(ICONF)=ISYS
                   EXIT
                ENDIF
             ENDDO
             NTYP_AB(ICONF)=AB_HELP%NTYP(ICONF)
             NIONS_AB(ICONF)=AB_HELP%NIONS(ICONF)
             DO INTYP=1, NTYP_AB(ICONF)
                TYPE_TMP(INTYP)=AB_HELP%TYPE_CONF(INTYP,ICONF)
                NITYP_TMP2(INTYP)=AB_HELP%NITYP(INTYP,ICONF)
                DO JNTYP0=1, MTYP_FF
                   IF(TYPE_TMP(INTYP) .EQ. TYPE_FF(JNTYP0)) THEN
                      LADD_AB(INTYP,ICONF)=JNTYP0
                      NITYP_AB(JNTYP0,ICONF)=NITYP_TMP2(INTYP)
                      NBMAX_TMP(JNTYP0)=NBMAX_TMP(JNTYP0)+NITYP_TMP2(INTYP)
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
             INIONS=0
             DO INTYP=1, NTYP_AB(ICONF)
                INTYP0=LADD_AB(INTYP,ICONF)
                DO INITYP=1, NITYP_AB(INTYP0,ICONF)
                   INIONS=INIONS+1
                   LNITYP_AB(INIONS,ICONF)=INITYP
                   LNIONS_AB(INITYP,INTYP0,ICONF)=INIONS
                   LNTYP_AB(INIONS,ICONF)=INTYP0
                ENDDO
             ENDDO
             ! Optional actions for CTIFOR
             AB%CTIFOR_PRESENT(ICONF)=AB_HELP%CTIFOR_PRESENT(ICONF)
             IF (AB%CTIFOR_PRESENT(ICONF)) THEN
                AB%CTIFOR_ALLCONF_OLD(ICONF)=AB_HELP%CTIFOR_ALLCONF_OLD(ICONF)
                AB%CTIFOR_ALLCONF_OLD(ICONF)=AB%CTIFOR_ALLCONF_OLD(ICONF)/FUNIT
             ENDIF
             ! End optional actions for CTIFOR
             !Obtain Bravais Lattice
             DO IXYZ=1, 3
                DO JXYZ=1, 3
                !A_AB(JXYZ,IXYZ,ICONF)=AB_HELP%A(JXYZ,IXYZ,ICONF)/AUTOA
                A_AB(JXYZ,IXYZ,ICONF)=AB_HELP%A(JXYZ,IXYZ,ICONF)
                B_AB(JXYZ,IXYZ,ICONF)=AB_HELP%B(JXYZ,IXYZ,ICONF)
                ENDDO
             ENDDO

             ! Obtain Wycoff positions
             DO INIONS=1, NIONS_AB(ICONF)
                DO IXYZ=1, 3
                   POSION_AB(IXYZ,INIONS,ICONF)=AB_HELP%POSION(IXYZ,INIONS,ICONF)/AUTOA
                ENDDO
             ENDDO
             
             ! Obtain cell volumes
             OMEGA_AB(ICONF) = AB_HELP%OMEGA(ICONF)
             !!! shift atoms back into primitve cell
             !! compute direct lattice to shift atoms back to primitive cell
             !CALL DIRECT_LATTICE( A_AB( 1:, 1:, ICONF ), &
             !                     B_AB( 1:, 1:, ICONF ), &
             !                     OMEGA_AB(ICONF ) )
             ! compute direct cordinates of positions
             !CALL KARDIR( NIONS_AB( ICONF ), POSION_AB( :, :, ICONF ),&
             !             B_AB( 1:, 1:, ICONF ) )
             !! shift direct positions back into the primitive box
             !CALL TOPRIM( NIONS_AB( ICONF ), POSION_AB( :, :, ICONF ) )
             !! rescale positions back to cartesian coordinates
             !CALL DIRKAR( NIONS_AB( ICONF ), POSION_AB( :, :, ICONF ),&
             !             A_AB( 1:, 1:, ICONF ) )
             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             !! Obtain reciprocal lattice
             !CALL RECIPROCAL_LATTICE( A_AB( 1:, 1:, ICONF), &
             !                         B_AB( 1:, 1:, ICONF), &
             !                         OMEGA_AB( ICONF ) )

             ! Obtain total energy
             TOTEN_AB(ICONF)=AB_HELP%TOTEN(ICONF)/EUNIT
             ! Obtain the forces
             DO INIONS=1, NIONS_AB(ICONF)
                DO IXYZ=1, 3
                   TIFOR_AB(IXYZ,INIONS,ICONF)=AB_HELP%TIFOR(IXYZ,INIONS,ICONF)/FUNIT
                ENDDO
             ENDDO
             ! Obtain stresses
             TSIF_AB(1,1,ICONF)=AB_HELP%TSIF(1,1,ICONF)
             TSIF_AB(2,2,ICONF)=AB_HELP%TSIF(2,2,ICONF)
             TSIF_AB(3,3,ICONF)=AB_HELP%TSIF(3,3,ICONF)
             TSIF_AB(2,1,ICONF)=AB_HELP%TSIF(2,1,ICONF)
             TSIF_AB(3,2,ICONF)=AB_HELP%TSIF(3,2,ICONF)
             TSIF_AB(3,1,ICONF)=AB_HELP%TSIF(3,1,ICONF)
             DO IXYZ=1, 2
                DO JXYZ=IXYZ+1, 3
                   TSIF_AB(IXYZ,JXYZ,ICONF)=TSIF_AB(JXYZ,IXYZ,ICONF)
                ENDDO
             ENDDO
             DO IXYZ=1, 3
                DO JXYZ=1, 3
                   TSIF_AB(JXYZ,IXYZ,ICONF)=TSIF_AB(JXYZ,IXYZ,ICONF)/SUNIT
                ENDDO
             ENDDO
          ENDDO
          ! If all structures in the ML_AB file contain a CTIFOR value (this
          ! is usually the case if ML_AB files are generated by VASP) then set
          ! corresponding flag to true. Later, this flag is used to determine
          ! whether the CTIFOR sections are written or if ML_ICRITERIA=3 can
          ! be used. ML_AB files with mixed CTIFOR presence are treated as if
          ! no CTIFOR values were present.
          IF (ALL(AB%CTIFOR_PRESENT(1:NCONF_AB))) THEN
             FF%LCTIFOR_PRESENT_IN_MLAB = .TRUE.
             IF (FF%ICRITERIA == 3) THEN
                CALL ml_tutor%advice( &
                   "Using Bayesian error threshold values from ML_AB file. &
                   &Sometimes better fitting results can be obtained with a &
                   &different threshold updating method (ML_ICRITERIA tag).")
             ENDIF
          ELSE
             FF%LCTIFOR_PRESENT_IN_MLAB = .FALSE.
             IF (FF%ICRITERIA == 3) THEN
                CALL ml_tutor%error( &
                   "Using Bayesian error threshold values from ML_AB file was &
                   &selected (ML_ICRITERIA=3) but not all structures contain &
                   &a CTIFOR value. Please use a different method to obtain &
                   &Bayesian error thresholds (ML_ICRITERIA tag).")
             ENDIF
             IF (ANY(AB%CTIFOR_PRESENT(1:NCONF_AB))) THEN
                CALL ml_tutor%alert( &
                   "Mixed CTIFOR presence in ML_AB file, i.e., some &
                   &structures contain CTIFOR values, some not. Present &
                   &CTIFOR values will be dropped in the output ML_ABN file.")
             ENDIF
          ENDIF
! Check whether the system stored in the database is the same as the current system
! This is only relevant for ISTART_FF.NE.3
          IF (ISTART_FF.NE.3) THEN
             DO ICONF=1, NCONF_AB
                LSYS_AB(ICONF)=1
                IF(NTYP_AB(ICONF).NE.NTYP_NEW) THEN
                   LSYS_AB(ICONF)=0
                ENDIF
                IF(LSYS_AB(ICONF).EQ.1) THEN
                   DO INTYP=1, NTYP_AB(ICONF)
                      INTYP0=LADD_AB(INTYP,ICONF)
                      IF(NITYP_AB(INTYP0,ICONF).NE.NITYP_ABN(INTYP0)) THEN
                         LSYS_AB(ICONF)=0
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDIF
! Check the number of configuration per model system.
          NCONF_SYS_AB=0
          DO ICONF=1, NCONF_AB
             NCONF_SYS_AB(ISYS_AB(ICONF))=NCONF_SYS_AB(ISYS_AB(ICONF))+1
          ENDDO
! Set the atomic reference energies read from ML_AB file
          DO INTYP=1, MTYP_FF
             EATOM(INTYP)=EATOM_TMP(INTYP)/EUNIT
          ENDDO
! Copying maximum volume amongst training structures
          AB%OMEGA_MAX = AB_HELP%OMEGA_MAX
        END SUBROUTINE READ_ABCAR

#endif

!****************************************************************************************************
! Some information needs to be obtained from the ML_FF file before allocating some arrays.
! FF%DESC_TYPE will be set here.
!****************************************************************************************************

        SUBROUTINE PRE_READ_FFCAR (FF,FFM)
           USE ML_FF_CONSTANT, ONLY: ML_FF_VERSION_WRITE
           USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR, ML_IO_WRITE
           USE ml_ff_tutor, ONLY: ml_tutor
           USE ml_ff_string, ONLY: str
           IMPLICIT NONE
           ! Input and output variables.
           LOGICAL                       :: DO_READ_DESC_TYPE
           TYPE (FF_PAR)                 :: FF
           TYPE (FFM_PAR)                :: FFM
           ! Local variables
           INTEGER                         :: ICOL
           INTEGER                         :: IROW
           INTEGER                         :: MTYP_LOC
           INTEGER                         :: NCOL
           INTEGER                         :: NMOD
           INTEGER, PARAMETER              :: NWRITE = 3
           INTEGER                         :: IDUMMY(1:NWRITE)
           REAL(q)                         :: RDUMMY(1:NWRITE)
           CHARACTER(LEN=2)                :: TYPE_LOC(1:NWRITE)
           LOGICAL                         :: LMLMB

           ! Open FFCAR file.
           CALL ROPEN_ML_FF(FF)
           ! Read the parameters.
           READ(FF%IU1) FF%ML_FF_VERSION_READ(1), FF%ML_FF_VERSION_READ(2), FF%ML_FF_VERSION_READ(3)
           ! Check version support.
           CALL CHECK_ML_FF_VERSION_COMPATIBILITY(FF%ML_FF_VERSION_READ, &
                                                  ML_FF_VERSION_WRITE)
           ! Descriptor type needs to be read here, so that it can be prefiously reopened.
           DO_READ_DESC_TYPE=.TRUE.
           IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
              IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
                 IF (FF%ML_FF_VERSION_READ(3).EQ.0) THEN
                    DO_READ_DESC_TYPE=.FALSE.
                 ENDIF
              ELSE IF (FF%ML_FF_VERSION_READ(2).LT.2) THEN
                 DO_READ_DESC_TYPE=.FALSE.
              ENDIF
           ENDIF 
           IF (DO_READ_DESC_TYPE) THEN
              ! Read the descriptor type
              READ(FF%IU1) FFM%DESC_TYPE
           ELSE
              ! Need to set DESC_TYPE dummywise here for older versions 
              FFM%DESC_TYPE=0
           ENDIF

           ! The number of element type
           READ(FF%IU1) MTYP_LOC
           NCOL=MTYP_LOC/NWRITE
           NMOD=MOD(MTYP_LOC,NWRITE)

           ! The element type
           DO ICOL=1, NCOL, 1
              READ(FF%IU1) (TYPE_LOC(IROW),IROW=1, NWRITE)
           ENDDO
           IF(NMOD .NE. 0) THEN
              READ(FF%IU1) (TYPE_LOC(IROW),IROW=1, NMOD)
           ENDIF

           ! Reference atomic energies (eV)
           DO ICOL=1, NCOL, 1
              READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
           ENDDO
           IF(NMOD .NE. 0) THEN
              READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
           ENDIF

           ! Atomic mass
           DO ICOL=1, NCOL, 1
              READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
           ENDDO
           IF(NMOD .NE. 0) THEN
              READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
           ENDIF

           ! Type of weighting
           READ(FF%IU1) IDUMMY(1)

           ! Weights
           READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)

           ! Mean square root error in energy (eV atom^-1)
           READ(FF%IU1) RDUMMY(1)

           ! Mean square root error in force (eV Angst^-1)
           READ(FF%IU1) RDUMMY(1)

           ! Mean square root error in stress tensor (kB)
           READ(FF%IU1) RDUMMY(1)

           ! Read LMLMB
           READ(FF%IU1) LMLMB
 
           ! Self-interaction correction?
           READ(FF%IU1) FFM%LSIC
          
           ! Close ML_FF file here
           CLOSE(FF%IU1)
        END SUBROUTINE PRE_READ_FFCAR

!****************************************************************************************************
! Read FFCAR file and overwrite the FF parameters by the parameters in this file.
!****************************************************************************************************

        SUBROUTINE OVERWRITE_FF_PAR_BY_PAR_IN_FFCAR (TAG_LIST,FF,FFM)
          USE ML_FF_CONSTANT, ONLY: ML_FF_VERSION_WRITE
          USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR, ML_IO_WRITE
          USE ML_FF_TAGLIST
          USE ml_ff_tutor, ONLY: ml_tutor
          USE ml_ff_string, ONLY: str
          IMPLICIT NONE
! Input and output variables.
          TYPE (TAGLIST), INTENT(INOUT) :: TAG_LIST
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
! Local variables.
          INTEGER                       :: I
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: IERR
          INTEGER                       :: INTYP
          INTEGER                       :: IROW
          INTEGER                       :: JNTYP
          INTEGER                       :: MTYP_LOC
          INTEGER                       :: NB_FFM(1:FF%MTYP)
          INTEGER                       :: NBTOT_FFM
          INTEGER                       :: NCOL
          INTEGER                       :: NMOD
          INTEGER, PARAMETER            :: NWRITE = 3
          CHARACTER(LEN=2)              :: TYPE_LOC(1:FF%MTYP)
          INTEGER                       :: IDUMMY(1:NWRITE)
          REAL(q)                       :: RDUMMY(1:NWRITE)
          LOGICAL                       :: DO_READ_DESC_TYPE
          INTEGER                       :: I_DESC_TYPE
! Initialize variables
          I_DESC_TYPE = 0
! Open FFCAR file.
          CALL ROPEN_ML_FF(FF)
! Read the parameters.
          ! Read the version number
          READ(FF%IU1) FF%ML_FF_VERSION_READ(1), FF%ML_FF_VERSION_READ(2), FF%ML_FF_VERSION_READ(3)
          ! Check version support.
          CALL CHECK_ML_FF_VERSION_COMPATIBILITY(FF%ML_FF_VERSION_READ, &
                                                 ML_FF_VERSION_WRITE)
          ! Check if version can read descriptor type
          DO_READ_DESC_TYPE=.TRUE.
          IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
             IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
                IF (FF%ML_FF_VERSION_READ(3).EQ.0) THEN
                   DO_READ_DESC_TYPE=.FALSE.
                ENDIF
             ELSE IF (FF%ML_FF_VERSION_READ(2).LT.2) THEN
                DO_READ_DESC_TYPE=.FALSE.
             ENDIF
          ENDIF
          IF (DO_READ_DESC_TYPE) THEN
             ! Read the descriptor type
             READ(FF%IU1) I_DESC_TYPE
          ENDIF
             
          ! The number of element type
          READ(FF%IU1) MTYP_LOC
          IF(FF%MTYP.NE.MTYP_LOC) THEN
             CALL ml_tutor%error("Atom types present in ML_FF &
                                  &inconsistent with &
                                  &prediction-only run (ML_ISTART = 2) &
                                  &with the given structure.")
          ENDIF
          NCOL=FF%MTYP/NWRITE
          NMOD=MOD(FF%MTYP,NWRITE)
          ! The element type
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (TYPE_LOC(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (TYPE_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, FF%MTYP
             IF(TYPE_LOC(I).NE.FF%TYPE(I)) THEN
                CALL ml_tutor%error("Inconsistent atom types in ML_FF for &
                   &prediction-only run (ML_ISTART = 2).")
             ENDIF
          ENDDO
          ! Reference atomic energies (eV)
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             FF%EATOM_REF(INTYP)=FF%EATOM_REF(INTYP)/EUNIT
          ENDDO
          ! Atomic mass
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
          ENDIF
          ! Type of weighting
          READ(FF%IU1) FF%IWEIGHT
          ! Weights
          READ(FF%IU1) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
          ! Mean square root error in energy (eV atom^-1)
          READ(FF%IU1) FF%STOTEN
          FF%STOTEN=FF%STOTEN/EUNIT
          ! Mean square root error in force (eV Angst^-1)
          READ(FF%IU1) FF%STIFOR
          FF%STIFOR=FF%STIFOR/FUNIT
          ! Mean square root error in stress tensor (kB)
          READ(FF%IU1) FF%STSIF
          FF%STSIF=FF%STSIF/SUNIT
! Read many-body parameters
          ! Many-body: true or false? - always true
          READ(FF%IU1) FF%LMLMB
          IF(FF%LMLMB) THEN
             ! Self-interaction correction?
             READ(FF%IU1) FFM%LSIC
             ! Type of similarity measure
             READ(FF%IU1) FFM%LSUPERVEC
             ! Weights for radial and angular descriptors
             READ(FF%IU1) FFM%W1,FFM%W2
             ! Type of cutoff function
             READ(FF%IU1) FFM%ICUT1,FFM%ICUT2
             ! Cutoff radius (Angst)
             READ(FF%IU1) FFM%RCUT1,FFM%RCUT2
             FFM%RCUT1=FFM%RCUT1/AUTOA
             FFM%RCUT2=FFM%RCUT2/AUTOA
             ! Broadening type
             READ(FF%IU1) FFM%IBROAD1,FFM%IBROAD2
             ! Broadening parameter
             READ(FF%IU1) FFM%SION1,FFM%SION2
             FFM%SION1=FFM%SION1/AUTOA
             FFM%SION2=FFM%SION2/AUTOA
             ! The maximum number of radial basis set
             READ(FF%IU1) FFM%MRB1,FFM%MRB2
             ! Read basis-function type
             FFM%BASIS_TYPE2=0
             IF (FF%ML_FF_VERSION_READ(1).GE.0 .AND. FF%ML_FF_VERSION_READ(2).GE.2 .AND. &
                 FF%ML_FF_VERSION_READ(3).GE.2) THEN
                READ(FF%IU1) FFM%BASIS_TYPE2
             ENDIF
             ! Read number of elements for summed basis
             IF (FFM%BASIS_TYPE2.EQ.4 .OR. FFM%BASIS_TYPE2.EQ.5) THEN
                READ(FF%IU1) FFM%MRB2_MAX
             ENDIF
             ! Optional read of reduced number of radial functions
             IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                READ(FF%IU1) FFM%RV%MRB2_RED
                READ(FF%IU1) FFM%RV%BASIS_TYPE2_RED
                IF (FFM%RV%BASIS_TYPE2_RED.EQ.4 .OR. FFM%RV%BASIS_TYPE2_RED.EQ.5) THEN
                   READ(FF%IU1) FFM%RV%MRB2_MAX_RED
                ENDIF
             ENDIF
             ! The number of grids for integrations
             READ(FF%IU1) FFM%NR1,FFM%NR2
             ! The number of spline grids 
             READ(FF%IU1) FFM%MSPL1,FFM%MSPL2
             FFM%NSPL1=FFM%MSPL1
             FFM%NSPL2=FFM%MSPL2
             ! Lmax for spherical harmonics
             READ(FF%IU1) FFM%LMAX1,FFM%LMAX2
             ! Optional read of reduced angular quantum number
             IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                READ(FF%IU1) FFM%RV%LMAX2_RED
             ENDIF
             ! Read ratio of high-resolution and low-resolution descriptor
             IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                READ(FF%IU1) FFM%RV%RATIO_DESC
                READ(FF%IU1) FFM%RV%RATIO_DESC_LOW_RES
                READ(FF%IU1) FFM%RV%FACTOR_DESC
             ENDIF
             IF(FFM%W1.GT.0.0_q) THEN
                FFM%MMVAR1=FFM%MRB1*MTYP_LOC
                FFM%NNVAR1=FFM%MMVAR1
             ENDIF
             ! Polynomial parameters for SOAP
             READ(FF%IU1) FFM%NHYP1,FFM%NHYP2
             ! Normalization
             READ(FF%IU1) FFM%LNORM1,FFM%LNORM2
             ! Radial window
             READ(FF%IU1) FFM%LWINDOW1,FFM%LWINDOW2
             IF(FFM%LWINDOW1) THEN
                ! Radial descriptor window
                READ(FF%IU1) FFM%IWINDOW1
             ENDIF
             IF(FFM%LWINDOW2) THEN
                ! Angular descriptor window
                READ(FF%IU1) FFM%IWINDOW2
             ENDIF
             ! Angular filtering
             READ(FF%IU1) FFM%LAFILT2
             IF(FFM%LAFILT2) THEN
                ! Type of angular filtering
                READ(FF%IU1) FFM%IAFILT2
                IF(FFM%IAFILT2.EQ.2) THEN
                   ! Damping parameter
                   READ(FF%IU1) FFM%AFILT2
                ENDIF
             ENDIF
             ! Radial metric
             READ(FF%IU1) FFM%LMETRIC1,FFM%LMETRIC2
             IF(FFM%LMETRIC1) THEN
                ! Radial descriptor metric
                READ(FF%IU1) FFM%NMETRIC1,FFM%RMETRIC1
             ENDIF
             IF(FFM%LMETRIC2) THEN
                ! Angular descriptor metric
                READ(FF%IU1) FFM%NMETRIC2,FFM%RMETRIC2
             ENDIF
             ! Variable transformation
             READ(FF%IU1) FFM%LVARTRAN1,FFM%LVARTRAN2
             IF(FFM%LVARTRAN1) THEN
                ! Radial variable transformation
                READ(FF%IU1) FFM%NVARTRAN1
             ENDIF
             IF(FFM%LVARTRAN2) THEN
                ! Angular variable transformation
                READ(FF%IU1) FFM%NVARTRAN2
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
! The number of radial basis sets for each angular component.
                NCOL=(FFM%LMAX2+1)/NWRITE
                NMOD=MOD(FFM%LMAX2+1,NWRITE)
                ! The number of radial basis sets for each angular channel
                DO ICOL=1, NCOL
                   READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                ENDIF
                ! Optional read for the number of reduced radial basis sets
                IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                   NCOL=(FFM%RV%LMAX2_RED+1)/NWRITE
                   NMOD=MOD(FFM%RV%LMAX2_RED+1,NWRITE)
                   ! The number of radial basis sets for each angular channel
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                   ENDDO
                   IF (NMOD .NE. 0) THEN
                      READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                   ENDIF
                ENDIF
! Lists of active descriptors.
                DO INTYP=1, MTYP_LOC
                   ! The number of active descriptors for XXX
                   READ(FF%IU1) FFM%NNVAR2(INTYP)
                   NCOL=(FFM%NNVAR2(INTYP))/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   ! List of active descriptors
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                   ENDIF
                   ! Optionale reducere descriptore
                   IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                      ! The number of active descriptors for XXX
                      READ(FF%IU1) FFM%RV%NNVAR2_RED(INTYP)
                      NCOL=(FFM%RV%NNVAR2_RED(INTYP))/NWRITE
                      NMOD=MOD(FFM%RV%NNVAR2_RED(INTYP),NWRITE)
                      ! List of active descriptors
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDIF
                ENDDO
! Set MMVAR2
                FFM%MMVAR2=FFM%NNVAR2(1)
                DO INTYP=2, MTYP_LOC, 1
                   FFM%MMVAR2=MAX(FFM%MMVAR2,FFM%NNVAR2(INTYP))
                ENDDO
! Optional set MMVAR2_RED
                IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                   FFM%RV%MMVAR2_RED=FFM%RV%NNVAR2_RED(1)
                   DO INTYP=2, MTYP_LOC, 1
                      FFM%RV%MMVAR2_RED=MAX(FFM%RV%MMVAR2_RED,FFM%RV%NNVAR2_RED(INTYP))
                   ENDDO
                ENDIF
! Lists of active descriptors for self-interaction correction if necessary.
                IF(FFM%LSIC) THEN
                   DO INTYP=1, MTYP_LOC
                      ! The number of active SIC descriptors for XXX 
                      READ(FF%IU1) FFM%NNVAR2_SIC(INTYP)
                      NCOL=(FFM%NNVAR2_SIC(INTYP))/NWRITE
                      NMOD=MOD(FFM%NNVAR2_SIC(INTYP),NWRITE)
                      ! List of active descriptors for SIC
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                   FFM%MMVAR2_SIC=FFM%NNVAR2_SIC(1)
                   DO INTYP=2, MTYP_LOC, 1
                      FFM%MMVAR2_SIC=MAX(FFM%MMVAR2_SIC,FFM%NNVAR2_SIC(INTYP))
                   ENDDO
                ENDIF
             ENDIF
             NCOL=MTYP_LOC/NWRITE
             NMOD=MOD(MTYP_LOC,NWRITE)
             ! The numbers of basis sets
             DO ICOL=1, NCOL
                READ(FF%IU1) (NB_FFM(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                READ(FF%IU1) (NB_FFM(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
! Set the required size automatically to the largest required if
! ML_ISTART=2
! If ML_MB was not read in then set to the first to
! first entry in ML_FF file and then to the largest 
             IF (TAG_LIST%GET_STATE("ML_MB") == TL_DEFAULT) THEN
                CALL TAG_LIST%SET_STATE("ML_MB", TL_OVERRIDE)
                FFM%MB=NB_FFM(1)+1
                DO INTYP=1, MTYP_LOC
                   IF (NB_FFM(INTYP) .GT. FFM%MB) THEN
                      FFM%MB=NB_FFM(INTYP)+1
                   ENDIF
                ENDDO
             ELSE
                DO INTYP=1, MTYP_LOC
                   IF (NB_FFM(INTYP) .GT. FFM%MB) THEN
                      CALL ml_tutor%error("ML_MB set in INCAR (" // str(FFM%MB) // &
                         ") is smaller than number of local reference configurations &
                         &in ML_FF (atom type " // str(INTYP) // ": " // &
                         str(NB_FFM(INTYP)) // "), please increase value or &
                         &remove tag from INCAR.")
                   ENDIF
                ENDDO
             ENDIF
! Further reading
             ! Type of the method to scale the energy data
             READ(FF%IU1) FFM%ISCALE_TOTEN
             ! Average energy per atom (eV atom^-1)
             READ(FF%IU1) FFM%TOTENAV
             FFM%TOTENAV=FFM%TOTENAV/EUNIT
             DO INTYP=1, MTYP_LOC
                NCOL=NB_FFM(INTYP)/NWRITE
                NMOD=MOD(NB_FFM(INTYP),NWRITE)
                ! Regression coefficients (eV atom^-1)
                ! For XXX
                DO ICOL=1, NCOL
                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                ENDIF
             ENDDO
             IF(FFM%W1 .GT. 0.0_q) THEN
                NCOL=FFM%MMVAR1/NWRITE
                NMOD=MOD(FFM%MMVAR1,NWRITE)
                DO INTYP=1, MTYP_LOC
                   ! C00 (radial descriptor)
                   ! For XXX
                   DO IB=1, NB_FFM(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             IF(FFM%W2 .GT. 0.0_q) THEN
                DO INTYP=1, MTYP_LOC
                   NCOL=FFM%NNVAR2(INTYP)/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   ! PS (angular descriptor) (-)
                   DO IB=1, NB_FFM(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDDO
                IF (I_DESC_TYPE.EQ.2 .OR. I_DESC_TYPE.EQ.4 .OR. I_DESC_TYPE.EQ.7 .OR. I_DESC_TYPE.EQ.9) THEN
                   DO INTYP=1, MTYP_LOC
                      NCOL=FFM%RV%NNVAR2_RED(INTYP)/NWRITE
                      NMOD=MOD(FFM%RV%NNVAR2_RED(INTYP),NWRITE)
                      ! PS (angular descriptor) (-)
                      DO IB=1, NB_FFM(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                         ENDDO
                         IF (NMOD .NE. 0) THEN
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
             ! Parameters for Bayesian error estimation
             ! Normalized noise parameter (-)     
             READ(FF%IU1) FFM%SIGV
             ! Normalized uncertainity parameter (-)        
             READ(FF%IU1) FFM%SIGW
             ! Variance in energies per atom in training data (eV atom^-1)     
             READ(FF%IU1) RDUMMY(1)
             ! Variance in forces in training data (eV Angst^-1)
             READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)
             ! Variance in stress tensors in training data (kB)      
             READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)
             READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)
             NBTOT_FFM=0
             DO INTYP=1, MTYP_LOC
                NBTOT_FFM=NBTOT_FFM+NB_FFM(INTYP)
             ENDDO
             NCOL=NBTOT_FFM/NWRITE
             NMOD=MOD(NBTOT_FFM,NWRITE)
             ! Read whether fitting was done with LFAST or not
             FF%LFAST_FROM_FILE = .FALSE.
             IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
                IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
                   IF (FF%ML_FF_VERSION_READ(3).GE.4) THEN
                      FF%LFAST_FROM_FILE = .TRUE.
                   ENDIF
                ELSEIF (FF%ML_FF_VERSION_READ(2).GT.2) THEN
                   FF%LFAST_FROM_FILE = .TRUE.
                ENDIF
             ELSE
                FF%LFAST_FROM_FILE = .TRUE.
             ENDIF
             IF (FF%LFAST_FROM_FILE) THEN
                READ(FF%IU1) FF%LFAST
                FF%LFAST_READ = FF%LFAST
             ENDIF

             ! FF%LFAST needs to be decided here
             IF (FFM%LSIC) THEN
                FF%LFAST = .FALSE.
             ELSE 
                IF (FF%LFORCESLOW) THEN
                   FF%LFAST = .FALSE.
                ELSE
                   ! Only reset here automatically if LFAST was not read (before ML_LOGFILE version 0.2.4)
                   IF (.NOT.FF%LFAST_FROM_FILE) THEN
                      FF%LFAST = .TRUE.
                   ENDIF
                ENDIF
             ENDIF
!             ! Covariance matrix (-)
!             DO IB=1, NBTOT_FFM
!                DO ICOL=1, NCOL
!                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
!                ENDDO
!                IF(NMOD.NE.0) THEN
!                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
!                ENDIF
!             ENDDO
          ENDIF
! Read charge parameters
! Keep this in for the moment to be compatible with old format (for LMLPC)
!          DO I=1, 3
!             READ(FF%IU1)
!          ENDDO
!          READ(FF%IU1) 
! Close FFCAR file.
          CLOSE(FF%IU1)
        END SUBROUTINE OVERWRITE_FF_PAR_BY_PAR_IN_FFCAR


!****************************************************************************************************
! Read FFCAR file.
!****************************************************************************************************

        SUBROUTINE READ_FFCAR (FF,FFM,PAR_SUP_HANDLE,ISYS_NEW)
          USE ML_FF_CONSTANT, ONLY: ML_FF_VERSION_WRITE
          USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR, PARALLEL_SUPER
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
! Input variables.
          TYPE (FF_PAR)                   :: FF
          TYPE (FFM_PAR)                  :: FFM
          TYPE (PARALLEL_SUPER)           :: PAR_SUP_HANDLE
          INTEGER         , INTENT(IN)    :: ISYS_NEW
! Local variables.
          INTEGER                         :: I
          INTEGER                         :: IB
          INTEGER                         :: ICOL
          INTEGER                         :: IERR
          INTEGER                         :: INTYP
          INTEGER                         :: INTYP0
          INTEGER                         :: IRB
          INTEGER                         :: IROW
          INTEGER                         :: IVAR
          INTEGER                         :: JNTYP
          INTEGER                         :: JNTYP0
          INTEGER                         :: JJNTYP0
          INTEGER                         :: JRB
          INTEGER                         :: L
          INTEGER                         :: LVAR_HELP(1:FFM%MMVAR2,1:FF%MTYP)
          INTEGER                         :: LVAR_HELP_RED(1:FFM%RV%MMVAR2_RED,1:FF%MTYP)
          INTEGER                         :: LVAR_SIC_HELP(1:FFM%MMVAR2_SIC,1:FF%MTYP)
          INTEGER                         :: MTYP_LOC
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NNVAR2_HELP1
          INTEGER                         :: NNVAR2_HELP2
          INTEGER                         :: NNVAR2_SIC_HELP1
          INTEGER                         :: NNVAR2_SIC_HELP2
          INTEGER                         :: NRB2(0:FFM%LMAX2)
          INTEGER                         :: NROW
          INTEGER, PARAMETER              :: NWRITE = 3
          CHARACTER(LEN=2)                :: TYPE_LOC(1:FF%MTYP)
          INTEGER                         :: IDUMMY(1:NWRITE)
          REAL(q)                         :: RDUMMY(1:NWRITE)
          LOGICAL                         :: LDUMMY
          LOGICAL                         :: DO_READ_DESC_TYPE
!          LOGICAL                         :: LREAD_SINV
          LOGICAL                         :: DO_READ_WRITECMAT
! Start profiling
          PROFILING_START('read_ffcar')
! Open FF%FFCAR file.
          CALL ROPEN_ML_FF(FF)
! Read the parameters.
          ! Read the version number
          READ(FF%IU1) FF%ML_FF_VERSION_READ(1), FF%ML_FF_VERSION_READ(2), FF%ML_FF_VERSION_READ(3)
          ! Check version support.
          CALL CHECK_ML_FF_VERSION_COMPATIBILITY(FF%ML_FF_VERSION_READ, &
                                                 ML_FF_VERSION_WRITE)
          DO_READ_DESC_TYPE=.TRUE.
          IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
             IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
                IF (FF%ML_FF_VERSION_READ(3).EQ.0) THEN
                   DO_READ_DESC_TYPE=.FALSE.
                ENDIF
             ELSE IF (FF%ML_FF_VERSION_READ(2).LT.2) THEN
                DO_READ_DESC_TYPE=.FALSE.
             ENDIF
          ENDIF
          IF (DO_READ_DESC_TYPE) THEN
             READ(FF%IU1) FFM%DESC_TYPE
          ELSE
             FFM%DESC_TYPE=0
          ENDIF
          ! The number of element type
          READ(FF%IU1) MTYP_LOC
          NCOL=MTYP_LOC/NWRITE
          NMOD=MOD(MTYP_LOC,NWRITE)
          ! The element type
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (TYPE_LOC(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (TYPE_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          ! Reference atomic energies (eV)
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             FF%EATOM_REF(INTYP)=FF%EATOM_REF(INTYP)/EUNIT
          ENDDO
          ! Atomic mass
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE) 
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
          ENDIF
          ! Type of weighting
          READ(FF%IU1) FF%IWEIGHT
          ! Weights
          READ(FF%IU1) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
          IF(FF%IWEIGHT.EQ.1) THEN
             FF%WTIFOR=FF%WTIFOR/FUNIT
             FF%WTOTEN=FF%WTOTEN/EUNIT
             FF%WTSIF=FF%WTSIF/SUNIT
          ENDIF
          ! Mean square root error in energy (eV atom^-1)
          READ(FF%IU1) FF%STOTEN
          FF%STOTEN=FF%STOTEN/EUNIT
          ! Mean square root error in force (eV Angst^-1)
          READ(FF%IU1) FF%STIFOR
          FF%STIFOR=FF%STIFOR/FUNIT
          ! Mean square root error in stress tensor (kB)
          READ(FF%IU1) FF%STSIF
          FF%STSIF=FF%STSIF/SUNIT
! Read many-body parameters
          ! Many-body: true or false?
          READ(FF%IU1) FF%LMLMB
          IF(FF%LMLMB) THEN
             ! Self-interaction correction?
             READ(FF%IU1) LDUMMY
             ! Type of similarity measure
             READ(FF%IU1) FFM%LSUPERVEC
             ! Weights for radial and angular descriptors
             READ(FF%IU1) FFM%W1,FFM%W2
             ! Type of cutoff function
             READ(FF%IU1) FFM%ICUT1,FFM%ICUT2
             ! Cutoff radius (Angst)
             READ(FF%IU1) FFM%RCUT1,FFM%RCUT2
             FFM%RCUT1=FFM%RCUT1/AUTOA
             FFM%RCUT2=FFM%RCUT2/AUTOA
             ! Broadening type
             READ(FF%IU1) FFM%IBROAD1,FFM%IBROAD2
             ! Broadening parameter
             READ(FF%IU1) FFM%SION1,FFM%SION2
             FFM%SION1=FFM%SION1/AUTOA
             FFM%SION2=FFM%SION2/AUTOA
            ! MRB1, MRB2: The maximum number of radial basis set 
             READ(FF%IU1) IDUMMY(1),IDUMMY(2)
            ! BASIS_TYPE2: Read in higher versions  
            IF (FF%ML_FF_VERSION_READ(1).GE.0 .AND. FF%ML_FF_VERSION_READ(2).GE.2 .AND. &
                FF%ML_FF_VERSION_READ(3).GE.2) THEN
               READ(FF%IU1) IDUMMY(1)
            ENDIF
            ! Optional read of MRB2_MAX
            IF (FFM%BASIS_TYPE2.EQ.4 .OR. FFM%BASIS_TYPE2.EQ.5) THEN
               READ(FF%IU1) IDUMMY(1)
            ENDIF
            ! Optional read of MRB_RED
             IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                READ(FF%IU1) IDUMMY(1)
                ! BASIS_TYPE2_RED: Read basis type for reduced descriptor
                READ(FF%IU1) IDUMMY(1)
                ! MRB2_MAX_RED
                IF (FFM%RV%BASIS_TYPE2_RED.EQ.4 .OR. FFM%RV%BASIS_TYPE2_RED.EQ.5) THEN
                   READ(FF%IU1) IDUMMY(1)
                ENDIF 
             ENDIF
            ! The number of grids for integrations
             READ(FF%IU1) FFM%NR1,FFM%NR2
            ! The number of spline grids
             READ(FF%IU1) FFM%MSPL1,FFM%MSPL2
            ! Lmax for spherical harmonics
             READ(FF%IU1) IDUMMY(1),IDUMMY(2)
            ! Optional read of LMAX2_RED
             IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                READ(FF%IU1) IDUMMY(1)
             ENDIF
            ! Optional read of desc ratio
             IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                READ(FF%IU1) RDUMMY(1)
                READ(FF%IU1) RDUMMY(1)
                READ(FF%IU1) RDUMMY(1)
             ENDIF
            ! Polynomial parameters for SOAP
             READ(FF%IU1) FFM%NHYP1,FFM%NHYP2
            ! Normalization 
             READ(FF%IU1) FFM%LNORM1,FFM%LNORM2
            ! Radial window
             READ(FF%IU1) FFM%LWINDOW1,FFM%LWINDOW2
             IF(FFM%LWINDOW1) THEN
               READ(FF%IU1) FFM%IWINDOW1
             ENDIF
             IF(FFM%LWINDOW2) THEN
               READ(FF%IU1) FFM%IWINDOW2
             ENDIF
            ! Angular filtering
             READ(FF%IU1) FFM%LAFILT2
             IF(FFM%LAFILT2) THEN
               ! Type of angular filtering
                READ(FF%IU1) FFM%IAFILT2
                IF(FFM%IAFILT2.EQ.2) THEN
                  ! Damping parameter
                   READ(FF%IU1) FFM%AFILT2
                ENDIF
             ENDIF
            ! Radial metric
             READ(FF%IU1) FFM%LMETRIC1,FFM%LMETRIC2
             IF(FFM%LMETRIC1) THEN
               READ(FF%IU1) FFM%NMETRIC1,FFM%RMETRIC1
             ENDIF
             IF(FFM%LMETRIC2) THEN
               READ(FF%IU1) FFM%NMETRIC2,FFM%RMETRIC2
             ENDIF
            ! Variable transformation
             READ(FF%IU1) FFM%LVARTRAN1,FFM%LVARTRAN2
             IF(FFM%LVARTRAN1) THEN
               READ(FF%IU1) FFM%NVARTRAN1
             ENDIF
             IF(FFM%LVARTRAN2) THEN
               READ(FF%IU1) FFM%NVARTRAN2
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
! Read the number of radial basis sets for each angular components.
                NCOL=(FFM%LMAX2+1)/NWRITE
                NMOD=MOD(FFM%LMAX2+1,NWRITE)
               ! The number of radial basis sets for each angular channel
                DO ICOL=1, NCOL, 1
                   READ(FF%IU1) (NRB2(NWRITE*(ICOL-1)+IROW),IROW=0, NWRITE-1)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (NRB2(NWRITE*NCOL+IROW),IROW=0, NMOD-1)
                ENDIF
! Optional read of the number of reduced radial basis sets
                IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                   NCOL=(FFM%RV%LMAX2_RED+1)/NWRITE
                   NMOD=MOD(FFM%RV%LMAX2_RED+1,NWRITE)
                   ! The number of radial basis sets for each angular channel
                   DO ICOL=1, NCOL, 1
                      READ(FF%IU1) (FFM%RV%NRB2_RED(NWRITE*(ICOL-1)+IROW),IROW=0, NWRITE-1)
                   ENDDO
                   IF (NMOD .NE. 0) THEN
                      READ(FF%IU1) (FFM%RV%NRB2_RED(NWRITE*NCOL+IROW),IROW=0, NMOD-1)
                   ENDIF
                ENDIF
! Read the list of active descriptors.
                DO INTYP=1, MTYP_LOC
                  !  The number of active descriptors for XXX
                   READ(FF%IU1) FFM%NNVAR2(INTYP)
                   NCOL=(FFM%NNVAR2(INTYP))/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                  ! List of active descriptors
                   DO ICOL=1, NCOL, 1
                      READ(FF%IU1) (LVAR_HELP(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (LVAR_HELP(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                   ENDIF
                  ! Optional read of list of active descriptors for reduced list
                   IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                      ! The number of active descriptors
                      READ(FF%IU1) FFM%RV%NNVAR2_RED(INTYP)
                      NCOL=(FFM%RV%NNVAR2_RED(INTYP))/NWRITE
                      NMOD=MOD(FFM%RV%NNVAR2_RED(INTYP),NWRITE)
                      ! List of active descriptors
                      DO ICOL=1, NCOL, 1
                         READ(FF%IU1) (LVAR_HELP_RED(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (LVAR_HELP_RED(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDIF
                ENDDO
! Read the list of active descriptors for SIC.
                IF(FFM%LSIC) THEN
                   DO INTYP=1, MTYP_LOC
                     ! The number of active SIC descriptors for XXX
                      READ(FF%IU1) FFM%NNVAR2_SIC(INTYP)
                      NCOL=(FFM%NNVAR2_SIC(INTYP))/NWRITE
                      NMOD=MOD(FFM%NNVAR2_SIC(INTYP),NWRITE)
                     ! List of active descriptors for SIC
                      DO ICOL=1, NCOL, 1
                         READ(FF%IU1) (LVAR_SIC_HELP(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (LVAR_SIC_HELP(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDIF
! Set lists.
! Determination of the numbers of variables related to the angular descriptor.
                FFM%LVAR=0
                FFM%LFLAG_VAR=.FALSE.
                FFM%WVAR=0.0_q
                IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4) THEN
                   FFM%RV%LVAR_RED=0
                   FFM%RV%LFLAG_VAR_RED=.FALSE.
                   FFM%RV%WVAR_RED=0.0_q
                ENDIF
                IF(FFM%LSIC) THEN
                   FFM%LVAR_SIC=0
                   FFM%LFLAG_VAR_SIC=.FALSE.
                ENDIF
                DO INTYP0=1, FF%MTYP
                   NNVAR2_HELP1=0
                   NNVAR2_HELP2=0
                   IF(FFM%LSIC) THEN
                      NNVAR2_SIC_HELP1=0
                      NNVAR2_SIC_HELP2=0
                   ENDIF
                   IF (FFM%DESC_TYPE.EQ.0) THEN 
                      DO JNTYP0=1, FF%MTYP
                         DO JJNTYP0=1, FF%MTYP
                            DO L=0, FFM%LMAX2
                               DO IRB=1, NRB2(L)
                                  DO JRB=IRB, NRB2(L)
                                     NNVAR2_HELP1=NNVAR2_HELP1+1
                                     DO IVAR=1, FFM%NNVAR2(INTYP0)
                                        IF (NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                           FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                           IF (JRB.EQ.IRB) THEN
                                              FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                           ELSE
                                              FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                           ENDIF
                                           EXIT
                                        ENDIF
                                     ENDDO
                                     IF (FFM%LSIC) THEN
                                        IF (JJNTYP0.EQ.JNTYP0) THEN
                                           NNVAR2_SIC_HELP1=NNVAR2_SIC_HELP1+1
                                           DO IVAR=1, FFM%NNVAR2_SIC(INTYP0)
                                              IF (NNVAR2_SIC_HELP1.EQ.LVAR_SIC_HELP(IVAR,INTYP0)) THEN
                                                 NNVAR2_SIC_HELP2=NNVAR2_SIC_HELP2+1
                                                 FFM%LVAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=NNVAR2_SIC_HELP2
                                                 FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=.TRUE.
                                                 EXIT
                                              ENDIF
                                           ENDDO
                                        ENDIF
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE IF (FFM%DESC_TYPE.EQ.1) THEN 
                      DO JNTYP0=1, FF%MTYP
                         DO L=0, FFM%LMAX2
                            DO IRB=1, NRB2(L)
                               DO JRB=IRB, NRB2(L)
                                  NNVAR2_HELP1=NNVAR2_HELP1+1
                                  DO IVAR=1, FFM%NNVAR2(INTYP0)
                                     IF (NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                        NNVAR2_HELP2=NNVAR2_HELP2+1
                                        FFM%LVAR(JRB,IRB,L,1,JNTYP0,INTYP0)=NNVAR2_HELP2
                                        FFM%LFLAG_VAR(JRB,IRB,L,1,JNTYP0,INTYP0)=.TRUE.
                                        IF (JRB.EQ.IRB) THEN
                                           FFM%WVAR(JRB,IRB,L,1,JNTYP0,INTYP0)=1.0_q
                                        ELSE
                                           FFM%WVAR(JRB,IRB,L,1,JNTYP0,INTYP0)=SQRT(2.0_q)
                                        ENDIF
                                        EXIT
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE IF (FFM%DESC_TYPE.EQ.2) THEN 
                      ! high resolution type-independent part
                      DO L=0, FFM%LMAX2
                         DO IRB=1, NRB2(L)
                            DO JRB=IRB, NRB2(L)
                               NNVAR2_HELP1=NNVAR2_HELP1+1
                               DO IVAR=1, FFM%NNVAR2(INTYP0)
                                  IF (NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                     NNVAR2_HELP2=NNVAR2_HELP2+1
                                     FFM%LVAR(JRB,IRB,L,1,1,INTYP0)=NNVAR2_HELP2
                                     FFM%LFLAG_VAR(JRB,IRB,L,1,1,INTYP0)=.TRUE.
                                     IF(JRB.EQ.IRB) THEN
                                        FFM%WVAR(JRB,IRB,L,1,1,INTYP0)=1.0_q
                                     ELSE
                                        FFM%WVAR(JRB,IRB,L,1,1,INTYP0)=SQRT(2.0_q)
                                     ENDIF
                                     EXIT
                                  ENDIF
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                      ! low resolution type-dependent part
                      NNVAR2_HELP1=0
                      NNVAR2_HELP2=0
                      DO JNTYP0=1, FF%MTYP
                         DO JJNTYP0=1, FF%MTYP
                            DO L=0, FFM%RV%LMAX2_RED
                               DO IRB=1, FFM%RV%NRB2_RED(L)
                                  DO JRB=IRB, FFM%RV%NRB2_RED(L)
                                     NNVAR2_HELP1=NNVAR2_HELP1+1
                                     DO IVAR=1, FFM%RV%NNVAR2_RED(INTYP0)
                                        IF (NNVAR2_HELP1.EQ.LVAR_HELP_RED(IVAR,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           FFM%RV%LVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                           FFM%RV%LFLAG_VAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                           IF (JRB.EQ.IRB) THEN
                                              FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                           ELSE
                                              FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                           ENDIF
                                           EXIT
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE IF (FFM%DESC_TYPE.EQ.3) THEN 
                      DO JNTYP0=1, FF%MTYP
                         DO JJNTYP0=1, FF%MTYP
                            DO L=0, FFM%LMAX2
                               DO IRB=1, NRB2(L)
                                  NNVAR2_HELP1=NNVAR2_HELP1+1
                                  DO IVAR=1, FFM%NNVAR2(INTYP0)
                                     IF (NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                        NNVAR2_HELP2=NNVAR2_HELP2+1
                                        FFM%LVAR(1,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                        FFM%LFLAG_VAR(1,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                        DO JRB=IRB, NRB2(L)
                                           IF (JRB.EQ.IRB) THEN
                                              FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                           ELSE
                                              FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                           ENDIF
                                        ENDDO
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE IF (FFM%DESC_TYPE.EQ.4) THEN 
                      ! First descriptor
                      DO JNTYP0=1, FF%MTYP
                         DO JJNTYP0=1, FF%MTYP
                            DO L=0, FFM%LMAX2
                               DO IRB=1, NRB2(L)
                                  DO JRB=IRB, NRB2(L)
                                     NNVAR2_HELP1=NNVAR2_HELP1+1
                                     DO IVAR=1, FFM%NNVAR2(INTYP0)
                                        IF (NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                           FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                           IF (JRB.EQ.IRB) THEN
                                              FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                           ELSE
                                              FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                           ENDIF
                                           EXIT
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                      ! Second descriptor
                      NNVAR2_HELP1=0
                      NNVAR2_HELP2=0
                      DO JNTYP0=1, FF%MTYP
                         DO JJNTYP0=1, FF%MTYP
                            DO L=0, FFM%RV%LMAX2_RED
                               DO IRB=1, FFM%RV%NRB2_RED(L)
                                  DO JRB=IRB, FFM%RV%NRB2_RED(L)
                                     NNVAR2_HELP1=NNVAR2_HELP1+1
                                     DO IVAR=1, FFM%RV%NNVAR2_RED(INTYP0)
                                        IF (NNVAR2_HELP1.EQ.LVAR_HELP_RED(IVAR,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           FFM%RV%LVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                           FFM%RV%LFLAG_VAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                           IF (JRB.EQ.IRB) THEN
                                              FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                           ELSE
                                              FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                           ENDIF
                                           EXIT
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF

! Read basis sets.
             NCOL=MTYP_LOC/NWRITE
             NMOD=MOD(MTYP_LOC,NWRITE)
            ! The numbers of basis sets
             DO ICOL=1, NCOL, 1
                READ(FF%IU1) (FFM%NB(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                READ(FF%IU1) (FFM%NB(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
            ! Type of the method to scale the energy data            
             READ(FF%IU1) FFM%ISCALE_TOTEN
            ! Average energy per atom (eV atom^-1)
             READ(FF%IU1) FFM%TOTENAV
             FFM%TOTENAV=FFM%TOTENAV/EUNIT
            DO INTYP=1, MTYP_LOC
                NCOL=FFM%NB(INTYP)/NWRITE
                NMOD=MOD(FFM%NB(INTYP),NWRITE)
               ! Regression coefficients (eV atom^-1)
                ! For XXX
#ifdef scaLAPACK
                CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM)
                IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,1),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,1),IROW=1, NMOD)
                   ENDIF
                ELSE
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD.NE.0) THEN
                      READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                   ENDIF
                ENDIF
                CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM,FFM%WMAT(1:,1:,INTYP), &
                     PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%NB(INTYP),1)
                CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                DO ICOL=1, NCOL
                   READ(FF%IU1) (FFM%WMAT(NWRITE*(ICOL-1)+IROW,1,INTYP),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (FFM%WMAT(NWRITE*NCOL+IROW,1,INTYP),IROW=1, NMOD)
                ENDIF
#endif
             ENDDO
             FFM%WMAT=FFM%WMAT/EUNIT
             IF(FFM%W1.GT.0.0_q) THEN
               NCOL=FFM%MMVAR1/NWRITE
                NMOD=MOD(FFM%MMVAR1,NWRITE)
                DO INTYP=1, MTYP_LOC
                  ! C00 (radial descriptor) (-)
                   ! For XXX
#ifdef scaLAPACK
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM,FFM%C00(1:,1:,INTYP), &
                        PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%MMVAR1,FFM%NB(INTYP))
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NB(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (FFM%C00(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (FFM%C00(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
                ENDDO
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
               DO INTYP=1, MTYP_LOC
                   NCOL=FFM%NNVAR2(INTYP)/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                  ! PS (angular descriptor) (-)
                   ! For XXX
#ifdef scaLAPACK
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM,FFM%PS(1:,1:,INTYP), &
                        PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%NNVAR2(INTYP),FFM%NB(INTYP))
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NB(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (FFM%PS(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (FFM%PS(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
               ENDDO
! Optional read of reduced angular descriptor components PS_RED
               IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                  DO INTYP=1, MTYP_LOC
                     NCOL=FFM%RV%NNVAR2_RED(INTYP)/NWRITE
                     NMOD=MOD(FFM%RV%NNVAR2_RED(INTYP),NWRITE)
                    ! PS (angular descriptor) (-)
#ifdef scaLAPACK
                     CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM_RED)
                     IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                        DO IB=1, FFM%NB(INTYP)
                           DO ICOL=1, NCOL
                              READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                           ENDDO
                           IF (NMOD .NE. 0) THEN
                               READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                           ENDIF
                        ENDDO
                     ELSE
                        DO IB=1, FFM%NB(INTYP)
                           DO ICOL=1, NCOL
                              READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                           ENDDO
                           IF (NMOD .NE. 0) THEN
                              READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                           ENDIF
                        ENDDO
                     ENDIF
                     CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM_RED,FFM%RV%PS_RED(1:,1:,INTYP), &
                          PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%RV%NNVAR2_RED(INTYP),FFM%NB(INTYP))
                     CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                     DO IB=1, FFM%NB(INTYP)
                        DO ICOL=1, NCOL
                           READ(FF%IU1) (FFM%RV%PS_RED(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                        ENDDO
                        IF (NMOD .NE. 0) THEN
                           READ(FF%IU1) (FFM%RV%PS_RED(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                        ENDIF
                     ENDDO
#endif
                  ENDDO
               ENDIF
             ENDIF
! Parameters for Bayesian error estimation
            ! Normalized noise parameter (-)
             READ(FF%IU1) FFM%SIGV
            ! Normalized uncertainity parameter (-)
             READ(FF%IU1) FFM%SIGW
! Here, FFM%SIG is assumed to be kBT. On the basis of the statistical
! dynamics, the variance should be given as sqrt(NIONS)*kBT.
            ! Variance in energies per atom in training data (eV atom^-1)
             READ(FF%IU1) FFM%SIG(1)
             FFM%SIG(1)=FFM%SIG(1)/EUNIT
            ! Variance in forces in training data (eV Angst^-1)
             READ(FF%IU1) FFM%SIG(2),FFM%SIG(3),FFM%SIG(4)
             FFM%SIG(2:4)=FFM%SIG(2:4)/FUNIT
            ! Variance in stress tensors in training data (kB)
             READ(FF%IU1) FFM%SIG(5),FFM%SIG(6),FFM%SIG(7)
             READ(FF%IU1) FFM%SIG(8),FFM%SIG(9),FFM%SIG(10)
             FFM%SIG(5:10)=FFM%SIG(5:10)/SUNIT
             FFM%SIG_SYS(1:10,ISYS_NEW)=FFM%SIG(1:10)
             FFM%NBTOT=0
             DO INTYP=1, MTYP_LOC
                FFM%NBTOT=FFM%NBTOT+FFM%NB(INTYP)
             ENDDO
             FFM%LBHEAD(1)=0
             DO INTYP=2, MTYP_LOC
                FFM%LBHEAD(INTYP)=FFM%LBHEAD(INTYP-1)+FFM%NB(INTYP-1)
             ENDDO

             ! From version 0.2.4 and higher LFAST needs to be dummy read here
             IF (FF%LFAST_FROM_FILE) THEN
                READ(FF%IU1) LDUMMY
             ENDIF

             ! Check version for reading of WRITECMAT 
             DO_READ_WRITECMAT = .TRUE.
             IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
                IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
                   IF (FF%ML_FF_VERSION_READ(3).LT.5) THEN
                      DO_READ_WRITECMAT = .FALSE.
                   ENDIF
                ELSE IF (FF%ML_FF_VERSION_READ(2).LT.2) THEN
                   DO_READ_WRITECMAT = .FALSE.
                ENDIF
             ENDIF

             ! Read WRITECMAT for ML_FF version >0.2.5 or set to .TRUE. for older versions
             ! except for 0.2.3 and 0.2.4 where it can be .FALSE. if LFAST was .TRUE.
             IF (DO_READ_WRITECMAT) THEN
                READ(FF%IU1) FFM%WRITECMAT
             ELSE
                FFM%WRITECMAT = .TRUE. 
                IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
                   IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
                      IF ((FF%ML_FF_VERSION_READ(3).EQ.3) .OR. (FF%ML_FF_VERSION_READ(3).EQ.4)) THEN
                         IF (FF%LFAST_FROM_FILE) THEN
                            IF (FF%LFAST_READ) FFM%WRITECMAT = .FALSE.
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF

!             ! Read last part of ML_FF containing the covariance matrix or inverse kernel matrix
!             ! We need to check the cases where SINV is written to the ML_FF file
!
!             LREAD_SINV = .TRUE.
!             IF (FF%ML_FF_VERSION_READ(1).EQ.0) THEN
!                IF (FF%ML_FF_VERSION_READ(2).EQ.2) THEN
!                   IF (FF%ML_FF_VERSION_READ(3).LT.3) THEN
!                      LREAD_SINV = .FALSE.
!                   ENDIF
!                ELSE IF (FF%ML_FF_VERSION_READ(2).LT.2) THEN
!                   LREAD_SINV = .FALSE.
!                ENDIF
!             ENDIF
!             IF (FF%LFAST_FROM_FILE) THEN
!                IF (.NOT.FF%LFAST_READ .AND. FFM%WRITECMAT) THEN
!                   LREAD_SINV = .FALSE.
!                ENDIF
!             ENDIF
!             IF (LREAD_SINV) THEN

             IF (.NOT.FFM%WRITECMAT) THEN 
                FF%LREAD_SINV_FROM_FILE = .TRUE.
#ifdef use_shmem
                IF (PAR_SUP_HANDLE%COMM_SHMEM%NODE_ME.EQ.0) THEN
#endif
                   DO INTYP0 = 1, FF%MTYP
                      NMOD=MOD(FFM%NB(INTYP0),NWRITE)
                      NCOL = FFM%NB(INTYP0) / NWRITE
                      DO IB = 1, FFM%NB(INTYP0)
!#ifdef scaLAPACK
                         DO ICOL = 1, NCOL
                            READ(FF%IU1) (FFM%SINV_ALL(NWRITE*(ICOL-1)+IROW,IB,INTYP0), IROW=1, NWRITE)
                         ENDDO
                         IF (NMOD.NE.0) THEN
                            READ(FF%IU1) (FFM%SINV_ALL(NWRITE*NCOL+IROW,IB,INTYP0), IROW=1, NMOD)
                         ENDIF
!#else
!                         DO ICOL = 1, NCOL
!                            READ(FF%IU1) (FFM%SINV(NWRITE*(ICOL-1)+IROW,IB,INTYP0), IROW=1, NWRITE)
!                         ENDDO
!                         IF (NMOD.NE.0) THEN
!                            READ(FF%IU1) (FFM%SINV(NWRITE*NCOL+IROW,IB,INTYP0), IROW=1, NMOD)
!                         ENDIF
!#endif
                      ENDDO
                   ENDDO
#ifdef use_shmem
                ENDIF
                CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
#endif
             ELSE
             ! Covariance matrix (-)
                IF (.NOT.FF%LFAST) THEN
                   NMOD=MOD(FFM%NBTOT,NWRITE)
                   NCOL=FFM%NBTOT/NWRITE
#ifdef scaLAPACK 
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NBTOT
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NBTOT
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM,FFM%CMAT(1:,1:), &
                        PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%NBTOT,FFM%NBTOT)
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NBTOT
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (FFM%CMAT(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (FFM%CMAT(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
                ENDIF
             ENDIF
             ! Be careful inserting reading here, since if LFAST=.TRUE. and LREAD_SINV=.TRUE.
             ! then only some cores reach properly to this point
          ENDIF
          CLOSE(FF%IU1)
! Stop profiling
          PROFILING_STOP('read_ffcar')
        END SUBROUTINE READ_FFCAR
!****************************************************************************************************


#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output regression results
!****************************************************************************************************

        SUBROUTINE OUTPUT_REG (COMM_WORLD,DIR_APP,DIR_LEN,IU3,MCONF,MIONS,NCONF,NIONS,NSTEP,REGCAR, &
                   TIFOR,TIFOR_MB,TOTEN,TOTEN_MB, &
                   TSIF,TSIF_MB)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
! Input variables
          TYPE (ML_MPI_PAR)                :: COMM_WORLD
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP     ! Location of directory passed from VASP.
          INTEGER, INTENT(IN)           :: DIR_LEN     ! Length of DIR_APP
          INTEGER, INTENT(IN)           :: IU3
          INTEGER, INTENT(IN)           :: MCONF
          INTEGER, INTENT(IN)           :: MIONS
          INTEGER, INTENT(IN)           :: NCONF
          INTEGER, INTENT(IN)           :: NIONS(:) !(1:NCONF)
          INTEGER, INTENT(IN)           :: NSTEP
          CHARACTER(LEN=10), INTENT(IN) :: REGCAR
          REAL(q), INTENT(IN)           :: TIFOR(:,:,:) !(1:3,1:MIONS,1:MCONF)
          REAL(q), INTENT(IN)           :: TIFOR_MB(:,:,:) !(1:3,1:MIONS,1:MCONF)
          REAL(q), INTENT(IN)           :: TOTEN(:) !(1:MCONF)
          REAL(q), INTENT(IN)           :: TOTEN_MB(:) !(1:MCONF)
          REAL(q), INTENT(IN)           :: TSIF(:,:,:) !(1:3,1:3,1:MCONF)
          REAL(q), INTENT(IN)           :: TSIF_MB(:,:,:) !(1:3,1:3,1:MCONF)
! Local variables
          INTEGER ICONF
          INTEGER IERR
          INTEGER INIONS
          INTEGER IXYZ
          INTEGER JXYZ
! Output regression results
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
             OPEN(UNIT=IU3,FILE=DIR_APP(1:DIR_LEN)// &
                  REGCAR,STATUS='UNKNOWN',IOSTAT=IERR)
! Check if error occured in opening file
             IF (IERR.NE.0) THEN
                CALL ml_tutor%error("Cannot open ML_REG file for writing.")
             ENDIF
             WRITE(IU3,1) NSTEP
             WRITE(IU3,2)
             DO ICONF=1, NCONF
                WRITE(IU3,3) TOTEN(ICONF)*EUNIT,(TOTEN_MB(ICONF))*EUNIT
             ENDDO
             WRITE(IU3,4)
             DO ICONF=1, NCONF
                DO INIONS=1, NIONS(ICONF)
                   DO IXYZ=1, 3
                      WRITE(IU3,3) TIFOR(IXYZ,INIONS,ICONF)*FUNIT, &
                                    (TIFOR_MB(IXYZ,INIONS,ICONF))*FUNIT
                   ENDDO
                ENDDO
             ENDDO
             WRITE(IU3,5)
             DO ICONF=1, NCONF
                DO IXYZ=1, 3
                   DO JXYZ=IXYZ, 3
                   WRITE(IU3,3) TSIF(JXYZ,IXYZ,ICONF)*SUNIT, &
                                 (TSIF_MB(JXYZ,IXYZ,ICONF))*SUNIT
                   ENDDO
                ENDDO
             ENDDO
             CLOSE(IU3)
          ENDIF
! Format
1         FORMAT (100("*")/5X,'Results at'1X,I10,1X,'MD step.')
2         FORMAT (100("=")/5X,'Total energies (eV)'/100("-"))
3         FORMAT (5X,E15.6,1X,E15.6)
4         FORMAT (100("=")/5X,'Forces (eV ang.^-1)'/100("-"))
5         FORMAT (100("=")/5X,'Stress (kbar)'/100("-"))
6         FORMAT (100("=")/5X,'Charge (e)'/100("-"))
        END SUBROUTINE OUTPUT_REG
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output new ab initio data
!****************************************************************************************************

        SUBROUTINE OUTPUT_ABNCAR (AB,FF,COMM_WORLD,A,ABNCAR,DIR_APP,DIR_LEN,EATOM_REF_FF,IU2,LADD,LNCONF,LNIONS,MB,MCONF,MITYP,MIONS,MTYP, &
                   NB,NCONF,NITYP,NIONS,NTYP,POMASS,POSION,SZNAM2,TIFOR,TOTEN,TSIF,TYPE)
          USE ML_FF_STRUCT, ONLY : ABINITIO, FF_PAR, ML_MPI_PAR
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
! Input variables
          TYPE (ABINITIO)               :: AB
          TYPE (FF_PAR)                 :: FF
          TYPE (ML_MPI_PAR)             :: COMM_WORLD
          REAL(q)          , INTENT(IN) :: A(:,:,:) !(1:3,1:3,1:MCONF)
          CHARACTER (LEN=10), INTENT(IN) :: ABNCAR
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
          INTEGER          , INTENT(IN) :: DIR_LEN
          REAL(q)          , INTENT(IN) :: EATOM_REF_FF(:) !(1:MTYP)
          INTEGER          , INTENT(IN) :: IU2
          INTEGER          , INTENT(IN) :: LADD(:,:) !(1:MTYP,1:MCONF)
          INTEGER          , INTENT(IN) :: LNCONF(:,:) !(1:MB,1:MTYP)
          INTEGER          , INTENT(IN) :: LNIONS(:,:) !(1:MB,1:MTYP)
          INTEGER          , INTENT(IN) :: MB
          INTEGER          , INTENT(IN) :: MCONF
          INTEGER          , INTENT(IN) :: MITYP
          INTEGER          , INTENT(IN) :: MIONS
          INTEGER          , INTENT(IN) :: MTYP
          INTEGER          , INTENT(IN) :: NB(:) !(1:MTYP)
          INTEGER          , INTENT(IN) :: NCONF
          INTEGER          , INTENT(IN) :: NITYP(:,:) !(1:MTYP,1:MCONF)
          INTEGER          , INTENT(IN) :: NIONS(:) !(1:MCONF)
          INTEGER          , INTENT(IN) :: NTYP(:) !(1:MCONF)
          REAL(q)          , INTENT(IN) :: POMASS(:) !(1:MTYP)
          REAL(q)          , INTENT(IN) :: POSION(:,:,:) !(1:3,1:MIONS,1:MCONF)
          CHARACTER(LEN=40), INTENT(IN) :: SZNAM2(:) !(1:MCONF)
          REAL(q)          , INTENT(IN) :: TIFOR(:,:,:) !(1:3,1:MIONS,1:MCONF)
          REAL(q)          , INTENT(IN) :: TOTEN(:) !(1:MCONF)
          REAL(q)          , INTENT(IN) :: TSIF(:,:,:) !(1:3,1:3,1:MCONF)
          CHARACTER (LEN=2), INTENT(IN) :: TYPE(:) !(1:MTYP)
! Local variables
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: ICONF
          INTEGER                       :: IERR
          INTEGER                       :: INIONS
          INTEGER                       :: INTYP,INTYP0
          INTEGER                       :: INITYP
          INTEGER                       :: IROW
          INTEGER                       :: IXYZ,JXYZ
          INTEGER                       :: LNCONF_LOC(1:MB,1:MTYP)
          INTEGER                       :: LNIONS_LOC(1:MB,1:MTYP)
          INTEGER                       :: NB_LOC(1:MTYP)
          INTEGER                       :: NCOL
          INTEGER                       :: NMOD
          INTEGER                       :: NWRITE
! Set constant.
          NWRITE=3
          NCOL=MTYP/NWRITE
          NMOD=MOD(MTYP,NWRITE)
! Set basis set lists
          LNCONF_LOC(1:MB,1:MTYP)=LNCONF(1:MB,1:MTYP)
          LNIONS_LOC(1:MB,1:MTYP)=LNIONS(1:MB,1:MTYP)
          NB_LOC(1:MTYP)=NB(1:MTYP)
! Output all the collected ab initio data
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
             OPEN(UNIT=IU2,FILE=DIR_APP(1:DIR_LEN)//ABNCAR,STATUS='UNKNOWN',IOSTAT=IERR)
! Check if error occured in opening file
             IF (IERR.NE.0) THEN
                CALL ml_tutor%error("Cannot open ML_ABN file for writing.")
             ENDIF
             WRITE(IU2,*) "1.0 Version"
             WRITE(IU2,5) NCONF
             WRITE(IU2,6) MTYP
             WRITE(IU2,7)
             DO ICOL=1, NCOL
                WRITE(IU2,8) (TYPE(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,8) (TYPE(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
             WRITE(IU2,9) MIONS
             WRITE(IU2,10) MITYP
             WRITE(IU2,11)
             DO ICOL=1, NCOL
                WRITE(IU2,*) (EATOM_REF_FF(NWRITE*(ICOL-1)+IROW)*EUNIT, IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,*) (EATOM_REF_FF(NWRITE*NCOL+IROW)*EUNIT,IROW=1, NMOD)
             ENDIF
             WRITE(IU2,12)
             DO ICOL=1, NCOL
                WRITE(IU2,*) (POMASS(NWRITE*(ICOL-1)+IROW)*MUNIT, IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,*) (POMASS(NWRITE*NCOL+IROW)*MUNIT,IROW=1, NMOD)
             ENDIF
             WRITE(IU2,13)
             DO ICOL=1, NCOL
                WRITE(IU2,14) (NB_LOC(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,14) (NB_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
             DO INTYP0=1, MTYP
                WRITE(IU2,15) TYPE(INTYP0)
                DO IB=1, NB_LOC(INTYP0)
                   WRITE(IU2,16) LNCONF_LOC(IB,INTYP0),LNIONS_LOC(IB,INTYP0)
                ENDDO
             ENDDO
             DO ICONF=1, NCONF
                WRITE(IU2,17) ICONF
                WRITE(IU2,18) SZNAM2(ICONF)
                WRITE(IU2,19) NTYP(ICONF)
                WRITE(IU2,20) NIONS(ICONF)
                WRITE(IU2,21)
                DO INTYP=1, NTYP(ICONF)
                   INTYP0=LADD(INTYP,ICONF)
                   WRITE(IU2,22) TYPE(INTYP0),NITYP(INTYP0,ICONF)
                ENDDO
                ! CTIFOR section is written under any of these conditions:
                ! * ISTART = 0   ... start from scratch.
                ! * ISTART = 3   ... MODE=select, generates CTIFOR corresponding
                !                    to ML_ICRITERIA setting.
                ! * ISTART = 1,4 ... only if CTIFOR is available for all
                !                    structures in ML_AB file.
                IF (FF%ISTART == 0 .OR. &
                    FF%ISTART == 3 .OR. &
                    FF%LCTIFOR_PRESENT_IN_MLAB) THEN
                   WRITE(IU2,30) FF%CALGO
                   WRITE(IU2,*) AB%CTIFOR_ALLCONF_NEW(ICONF)*FUNIT
                ENDIF
                WRITE(IU2,23)
                DO IXYZ=1, 3
                   WRITE(IU2,*) (A(JXYZ,IXYZ,ICONF)*AUTOA,JXYZ=1, 3)
                ENDDO
                WRITE(IU2,24)
                DO INIONS=1, NIONS(ICONF)
                   WRITE(IU2,*) (POSION(IXYZ,INIONS,ICONF)*AUTOA,IXYZ=1, 3)
                ENDDO
                WRITE(IU2,25)
                WRITE(IU2,*) TOTEN(ICONF)*EUNIT
                WRITE(IU2,26)
                DO INIONS=1, NIONS(ICONF)
                   WRITE(IU2,*) (TIFOR(IXYZ,INIONS,ICONF)*FUNIT,IXYZ=1, 3)
                ENDDO
                WRITE(IU2,27)
                WRITE(IU2,*) TSIF(1,1,ICONF)*SUNIT,TSIF(2,2,ICONF)*SUNIT,TSIF(3,3,ICONF)*SUNIT
                WRITE(IU2,28)
                WRITE(IU2,*) TSIF(1,2,ICONF)*SUNIT,TSIF(3,2,ICONF)*SUNIT,TSIF(3,1,ICONF)*SUNIT
             ENDDO
             CLOSE(IU2)
          ENDIF
5         FORMAT (50("*")/5X,'The number of configurations'/50("-")/5X,I6)
6         FORMAT (50("*")/5X,'The maximum number of atom type'/50("-")/5X,I3)
7         FORMAT (50("*")/5X,'The atom types in the data file'/50("-"))
8         FORMAT (5X,A2,1X,A2,1X,A2)
9         FORMAT (50("*")/5X,'The maximum number of atoms per system'/50("-")/5X,I10)
10        FORMAT (50("*")/5X,'The maximum number of atoms per atom type'/50("-")/5X,I10)
11        FORMAT (50("*")/5X,'Reference atomic energy (eV)'/50("-"))
12        FORMAT (50("*")/5X,'Atomic mass'/50("-"))
13        FORMAT (50("*")/5X,'The numbers of basis sets per atom type'/50("-"))
14        FORMAT (5X,I5,1X,I5,1X,I5)
15        FORMAT (50("*")/5X,'Basis set for'1X,A2/50("-"))
16        FORMAT (5X,I6,1X,I6)
17        FORMAT (50("*")/5X,'Configuration num.'1X,I6)
18        FORMAT (50("=")/5X,'System name'/50("-")/5X,A40)
19        FORMAT (50("=")/5X,'The number of atom types'/50("-")/5X,I3)
20        FORMAT (50("=")/5X,'The number of atoms'/50("-")/5X,I6)
21        FORMAT (50("*")/5X,'Atom types and atom numbers'/50("-"))
22        FORMAT (5X,A2,1X,I6)
23        FORMAT (50("=")/5X,'Primitive lattice vectors (ang.)'/50("-"))
24        FORMAT (50("=")/5X,'Atomic positions (ang.)'/50("-"))
25        FORMAT (50("=")/5X,'Total energy (eV)'/50("-"))
26        FORMAT (50("=")/5X,'Forces (eV ang.^-1)'/50("-"))
27        FORMAT (50("=")/5X,'Stress (kbar)'/50("-")/ &
                          5X,'XX YY ZZ'/50("-"))
28        FORMAT (50("-")/5X,'XY YZ ZX'/50("-"))
29        FORMAT (50("=")/5X,'Charges (e)'/50("-"))
30        FORMAT (50("=")/5X,'CTIFOR ( CALGO = ',I3,' )'/50("-"))
        END SUBROUTINE OUTPUT_ABNCAR
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
!> Output new FFCAR (FFNCAR) file.
!>
!> For version history see documentation of ML_FF_VERSION_WRITE in ml_ff_constants.F
!****************************************************************************************************
        SUBROUTINE OUT_FFNCAR (FF,FFM,NCONF,PAR_SUP_HANDLE,DIR_APP,DIR_LEN)
          USE ML_FF_CONSTANT, ONLY: ML_FF_VERSION_WRITE
          USE ML_FF_STRUCT, ONLY : FF_PAR,FFM_PAR,PARALLEL_SUPER
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
! Input variables
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
          TYPE (PARALLEL_SUPER)         :: PAR_SUP_HANDLE
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
          INTEGER          , INTENT(IN) :: DIR_LEN
          INTEGER          , INTENT(IN) :: NCONF !< Number of training structures.
! Local variables
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: IERR
          INTEGER                       :: INTYP0
          INTEGER                       :: JNTYP0
          INTEGER                       :: JJNTYP0
          INTEGER                       :: IRB
          INTEGER                       :: IROW
          INTEGER                       :: IVAR
          INTEGER                       :: JRB
          INTEGER                       :: L
          INTEGER                       :: LVAR_HELP(1:FFM%MMVAR2)
          INTEGER                       :: LVAR_HELP_RED(1:FFM%RV%MMVAR2_RED)
          INTEGER                       :: LVAR_SIC_HELP(1:FFM%MMVAR2_SIC)
          INTEGER                       :: NCOL
          INTEGER,PARAMETER             :: NROW = 3
          INTEGER                       :: NNVAR2_HELP1
          INTEGER                       :: NNVAR2_HELP2
          INTEGER                       :: NNVAR2_HELP1_RED
          INTEGER                       :: NNVAR2_HELP2_RED
! If FF%IFF=.TRUE., which means the presence of MLFF, output the FF parameters
          IF(FF%IFF) THEN
             PROFILING_START('out_ffncar')
             ! Store the creation date and time for log file and ML_FFN header.
             FF%ML_FF_CREATION_DATETIME = NOW()
             IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                ! Write ASCII header.
                CALL WRITE_FFN_HEADER(.TRUE., FF, FFM, NCONF, DIR_APP, DIR_LEN)
                ! Open file for binary writing.
                OPEN(FF%IU2,FILE=DIR_APP(1:DIR_LEN)//FF%FFNCAR, &
                     STATUS='OLD',FORM='UNFORMATTED',           &
                     ACCESS='STREAM',POSITION='APPEND',IOSTAT=IERR)
                ! Check if error occured in opening file
                IF (IERR.NE.0) THEN
                   CALL ml_tutor%error("Cannot open ML_FFN file for writing.")
                ENDIF

! Output the general parameters.
                ! Write the version number
                WRITE(FF%IU2) ML_FF_VERSION_WRITE(1), ML_FF_VERSION_WRITE(2), ML_FF_VERSION_WRITE(3)
                ! Write the type of angular descriptor, even if there is 
                ! no angular descriptor, just default it to 0
                WRITE(FF%IU2) (FFM%DESC_TYPE)
                ! The number of element type
                WRITE(FF%IU2) FF%MTYP 
                NCOL=FF%MTYP/NROW
                ! The element type
                DO ICOL=1, NCOL, 1
                   IROW=NROW*(ICOL-1)
                   WRITE(FF%IU2) (FF%TYPE(INTYP0), INTYP0=IROW+1,IROW+NROW)
                ENDDO
                IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                   WRITE(FF%IU2) (FF%TYPE(INTYP0), INTYP0=NCOL*NROW+1,NCOL*NROW+MOD(FF%MTYP,NROW))
                ENDIF
                ! Reference atomic energies (eV)
                NCOL=FF%MTYP/NROW
                DO ICOL=1, NCOL, 1
                   IROW=NROW*(ICOL-1)
                   WRITE(FF%IU2) (FF%EATOM_REF(INTYP0)*EUNIT, INTYP0=IROW+1,IROW+NROW)
                ENDDO
                IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                   WRITE(FF%IU2) (FF%EATOM_REF(INTYP0)*EUNIT, INTYP0=NCOL*NROW+1,NCOL*NROW+MOD(FF%MTYP,NROW))
                ENDIF
                ! Atomic mass
                NCOL=FF%MTYP/NROW
                DO ICOL=1, NCOL, 1
                   IROW=NROW*(ICOL-1)
                   WRITE(FF%IU2) (FF%POMASS(INTYP0)*MUNIT, INTYP0=IROW+1,IROW+NROW)
                ENDDO
                IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                   WRITE(FF%IU2) (FF%POMASS(INTYP0)*MUNIT, INTYP0=NCOL*NROW+1,NCOL*NROW+MOD(FF%MTYP,NROW))
                ENDIF
                ! Type of weighting
                WRITE(FF%IU2) FF%IWEIGHT
                IF(FF%IWEIGHT.EQ.1) THEN
                   WRITE(FF%IU2) FF%WTOTEN*EUNIT,FF%WTIFOR*FUNIT,FF%WTSIF*SUNIT
                ELSE
                   WRITE(FF%IU2) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
                ENDIF
                ! Mean square root error in energy (eV atom^-1)
                WRITE(FF%IU2) FF%STOTEN*EUNIT
                ! Mean square root error in force (eV Angst^-1)
                WRITE(FF%IU2) FF%STIFOR*FUNIT
                ! Mean square root error in stress tensor (kB)
                WRITE(FF%IU2) FF%STSIF*SUNIT
! Write Many-body parameters.
                ! Many-body: true or false? - always true
                WRITE(FF%IU2) FF%LMLMB
             ENDIF
! Output the parameters for many-body interactions.
             IF(FF%LMLMB) THEN
                IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                   ! Self-interaction correction?
                   WRITE(FF%IU2) FFM%LSIC
                   ! Type of similarity measure
                   WRITE(FF%IU2) FFM%LSUPERVEC
                   ! Weights for radial and angular descriptors
                   WRITE(FF%IU2) FFM%W1,FFM%W2
                   ! Type of cutoff function
                   WRITE(FF%IU2) FFM%ICUT1,FFM%ICUT2
                   ! Cutoff radius (Angst)
                   WRITE(FF%IU2) FFM%RCUT1*AUTOA,FFM%RCUT2*AUTOA
                   ! Broadening type
                   WRITE(FF%IU2) FFM%IBROAD1,FFM%IBROAD2
                   ! Broadening parameter
                   WRITE(FF%IU2) FFM%SION1*AUTOA,FFM%SION2*AUTOA
                   ! The maximum number of radial basis set
                   WRITE(FF%IU2) FFM%MRB1,FFM%MRB2
                   ! Basis-function type
                   WRITE(FF%IU2) FFM%BASIS_TYPE2
                   ! Optional write of Maximum MRB2 distance for summation
                   IF (FFM%BASIS_TYPE2.EQ.4 .OR. FFM%BASIS_TYPE2.EQ.5) THEN
                      WRITE(FF%IU2) FFM%MRB2_MAX
                   ENDIF 
                   ! Optional write for reduced descriptore
                   IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                      WRITE(FF%IU2) FFM%RV%MRB2_RED
                      ! Basis type
                      WRITE(FF%IU2) FFM%RV%BASIS_TYPE2_RED
                      ! Maximum distance for summation
                      IF (FFM%RV%BASIS_TYPE2_RED.EQ.4 .OR. FFM%RV%BASIS_TYPE2_RED.EQ.5) THEN
                         WRITE(FF%IU2) FFM%RV%MRB2_MAX_RED
                      ENDIF
                   ENDIF 
                   ! The number of grids for integrations
                   WRITE(FF%IU2) FFM%NR1,FFM%NR2
                   ! The number of spline grids
                   WRITE(FF%IU2) FFM%NSPL1,FFM%NSPL2
                   ! Lmax for spherical harmonics
                   WRITE(FF%IU2) FFM%LMAX1,FFM%LMAX2
                   ! Optional write for reduced descriptore
                   IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                      WRITE(FF%IU2) FFM%RV%LMAX2_RED
                   ENDIF
                   ! Optional write of desc ratio
                   IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                      WRITE(FF%IU2) FFM%RV%RATIO_DESC
                      WRITE(FF%IU2) FFM%RV%RATIO_DESC_LOW_RES
                      WRITE(FF%IU2) FFM%RV%FACTOR_DESC
                   ENDIF
                   ! Polynomial parameters for SOAP
                   WRITE(FF%IU2) FFM%NHYP1,FFM%NHYP2
                   ! Normalization
                   WRITE(FF%IU2) FFM%LNORM1,FFM%LNORM2
                   ! Radial window
                   WRITE(FF%IU2) FFM%LWINDOW1,FFM%LWINDOW2
                   IF(FFM%LWINDOW1) THEN
                      ! Radial descriptor window
                      WRITE(FF%IU2) FFM%IWINDOW1
                   ENDIF
                   IF(FFM%LWINDOW2) THEN
                      ! Angular descriptor window
                      WRITE(FF%IU2) FFM%IWINDOW2
                   ENDIF
                   ! Angular filtering
                   WRITE(FF%IU2) FFM%LAFILT2
                   IF(FFM%LAFILT2) THEN
                      ! Type of angular filtering
                      WRITE(FF%IU2) FFM%IAFILT2
                      IF(FFM%IAFILT2.EQ.2) THEN
                         ! Damping parameter
                         WRITE(FF%IU2) FFM%AFILT2
                      ENDIF
                   ENDIF
                   ! Radial metric
                   WRITE(FF%IU2) FFM%LMETRIC1,FFM%LMETRIC2
                   IF(FFM%LMETRIC1) THEN
                      ! Radial descriptor metric
                      WRITE(FF%IU2) FFM%NMETRIC1,FFM%RMETRIC1*AUTOA
                   ENDIF
                   IF(FFM%LMETRIC2) THEN
                      ! Angular descriptor metric
                      WRITE(FF%IU2) FFM%NMETRIC2,FFM%RMETRIC2*AUTOA
                   ENDIF
                   ! Variable transformation
                   WRITE(FF%IU2) FFM%LVARTRAN1,FFM%LVARTRAN2
                   IF(FFM%LVARTRAN1) THEN
                      ! Radial variable transformation
                      WRITE(FF%IU2) FFM%NVARTRAN1
                   ENDIF
                   IF(FFM%LVARTRAN2) THEN
                      ! Angular variable transformation
                      WRITE(FF%IU2) FFM%NVARTRAN2
                   ENDIF
                   IF(FFM%W2.GT.0.0_q) THEN
                      NCOL=(FFM%LMAX2+1)/NROW
                      ! The number of radial basis sets for each angular channel
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (FFM%NRB2(L), L=IROW,IROW+NROW-1)
                      ENDDO
                      IF(MOD(FFM%LMAX2+1,NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (FFM%NRB2(L), L=NCOL*NROW,FFM%LMAX2)
                      ENDIF
                      ! Optional write for reduced descriptore
                      IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                         NCOL=(FFM%RV%LMAX2_RED+1)/NROW
                         ! The number of radial basis sets for each angular channel
                         DO ICOL=1, NCOL, 1
                            IROW=NROW*(ICOL-1)
                            WRITE(FF%IU2) (FFM%RV%NRB2_RED(L), L=IROW,IROW+NROW-1)
                         ENDDO
                         IF (MOD(FFM%RV%LMAX2_RED+1,NROW) .GT. 0) THEN
                            WRITE(FF%IU2) (FFM%RV%NRB2_RED(L), L=NCOL*NROW,FFM%RV%LMAX2_RED)
                         ENDIF
                      ENDIF
                      DO INTYP0=1, FF%MTYP
                         LVAR_HELP=0
                         NNVAR2_HELP1=0
                         NNVAR2_HELP2=0
                         IF (FFM%DESC_TYPE.EQ.0) THEN
                            DO JNTYP0=1, FF%MTYP
                               DO JJNTYP0=1, FF%MTYP
                                  DO L=0, FFM%LMAX2
                                     DO IRB=1, FFM%NRB2(L)
                                        DO JRB=IRB, FFM%NRB2(L)
                                           NNVAR2_HELP1=NNVAR2_HELP1+1
                                           IF(FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                              NNVAR2_HELP2=NNVAR2_HELP2+1
                                              LVAR_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                           ENDIF
                                        ENDDO
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
                            DO JNTYP0=1, FF%MTYP
                               DO L=0, FFM%LMAX2
                                  DO IRB=1, FFM%NRB2(L)
                                     DO JRB=IRB, FFM%NRB2(L)
                                        NNVAR2_HELP1=NNVAR2_HELP1+1
                                        IF(FFM%LFLAG_VAR(JRB,IRB,L,1,JNTYP0,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           LVAR_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ELSE IF (FFM%DESC_TYPE.EQ.2) THEN
                            ! high resolution type-independent descriptor
                            DO L=0, FFM%LMAX2
                               DO IRB=1, FFM%NRB2(L)
                                  DO JRB=IRB, FFM%NRB2(L)
                                     NNVAR2_HELP1=NNVAR2_HELP1+1
                                     IF(FFM%LFLAG_VAR(JRB,IRB,L,1,1,INTYP0)) THEN
                                        NNVAR2_HELP2=NNVAR2_HELP2+1
                                        LVAR_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                            ! low resolution type-dependent descriptor
                            NNVAR2_HELP1_RED=0
                            NNVAR2_HELP2_RED=0
                            DO JNTYP0=1, FF%MTYP
                               DO JJNTYP0=1, FF%MTYP
                                  DO L=0, FFM%RV%LMAX2_RED
                                     DO IRB=1, FFM%RV%NRB2_RED(L)
                                        DO JRB=IRB, FFM%RV%NRB2_RED(L)
                                           NNVAR2_HELP1_RED=NNVAR2_HELP1_RED+1
                                           IF(FFM%RV%LFLAG_VAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                              NNVAR2_HELP2_RED=NNVAR2_HELP2_RED+1
                                              LVAR_HELP_RED(NNVAR2_HELP2_RED)=NNVAR2_HELP1_RED
                                           ENDIF
                                        ENDDO
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ELSE IF (FFM%DESC_TYPE.EQ.3) THEN
                            DO JNTYP0=1, FF%MTYP
                               DO JJNTYP0=1, FF%MTYP
                                  DO L=0, FFM%LMAX2
                                     DO IRB=1, FFM%NRB2(L)
                                        NNVAR2_HELP1=NNVAR2_HELP1+1
                                        IF(FFM%LFLAG_VAR(1,IRB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           LVAR_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ELSE IF (FFM%DESC_TYPE.EQ.4) THEN
                            ! First descriptor
                            DO JNTYP0=1, FF%MTYP
                               DO JJNTYP0=1, FF%MTYP
                                  DO L=0, FFM%LMAX2
                                     DO IRB=1, FFM%NRB2(L)
                                        DO JRB=IRB, FFM%NRB2(L)
                                           NNVAR2_HELP1=NNVAR2_HELP1+1
                                           IF(FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                              NNVAR2_HELP2=NNVAR2_HELP2+1
                                              LVAR_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                           ENDIF
                                        ENDDO
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                            ! Second descriptor
                            NNVAR2_HELP1_RED=0
                            NNVAR2_HELP2_RED=0
                            DO JNTYP0=1, FF%MTYP
                               DO JJNTYP0=1, FF%MTYP
                                  DO L=0, FFM%RV%LMAX2_RED
                                     DO IRB=1, FFM%RV%NRB2_RED(L)
                                        DO JRB=IRB, FFM%RV%NRB2_RED(L)
                                           NNVAR2_HELP1_RED=NNVAR2_HELP1_RED+1
                                           IF(FFM%RV%LFLAG_VAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                              NNVAR2_HELP2_RED=NNVAR2_HELP2_RED+1
                                              LVAR_HELP_RED(NNVAR2_HELP2_RED)=NNVAR2_HELP1_RED
                                           ENDIF
                                        ENDDO
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDIF 
                         NCOL=NNVAR2_HELP2/NROW
                         ! The number of active descriptors for XXX
                         WRITE(FF%IU2) NNVAR2_HELP2
                         ! List of active descriptors
                         DO ICOL=1, NCOL, 1
                            IROW=NROW*(ICOL-1)
                            WRITE(FF%IU2) (LVAR_HELP(IVAR), IVAR=IROW+1,IROW+NROW)
                         ENDDO
                         IF(MOD(NNVAR2_HELP2,NROW) .GT. 0) THEN
                            WRITE(FF%IU2) (LVAR_HELP(IVAR), IVAR=NCOL*NROW+1,NNVAR2_HELP2)
                         ENDIF
                         IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                            NCOL=NNVAR2_HELP2_RED/NROW
                            ! The number of active descriptors for XXX
                            WRITE(FF%IU2) NNVAR2_HELP2_RED
                            ! List of active descriptors
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (LVAR_HELP_RED(IVAR), IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF (MOD(NNVAR2_HELP2_RED,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (LVAR_HELP_RED(IVAR), IVAR=NCOL*NROW+1,NNVAR2_HELP2_RED)
                            ENDIF
                         ENDIF
                      ENDDO
                      ! Lists of active descriptors for self-interaction correction if necessary.
                      IF(FFM%LSIC) THEN
                         DO INTYP0=1, FF%MTYP
                            LVAR_SIC_HELP=0
                            NNVAR2_HELP1=0
                            NNVAR2_HELP2=0
                            DO JNTYP0=1, FF%MTYP
                               DO L=0, FFM%LMAX2
                                  DO IRB=1, FFM%NRB2(L)
                                     DO JRB=IRB, FFM%NRB2(L)
                                        NNVAR2_HELP1=NNVAR2_HELP1+1
                                        IF(FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           LVAR_SIC_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                            NCOL=NNVAR2_HELP2/NROW
                           ! The number of active SIC descriptors for XXX
                            WRITE(FF%IU2) NNVAR2_HELP2
                            ! List of active descriptors for SIC
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (LVAR_SIC_HELP(IVAR), IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(NNVAR2_HELP2,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (LVAR_SIC_HELP(IVAR), IVAR=NCOL*NROW+1,NNVAR2_HELP2)
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDIF
                   ! The numbers of basis sets
                   NCOL=FF%MTYP/NROW
                   DO ICOL=1, NCOL, 1
                      IROW=NROW*(ICOL-1)
                      WRITE(FF%IU2) (FFM%NB(INTYP0), INTYP0=IROW+1,IROW+NROW)
                   ENDDO
                   IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                      WRITE(FF%IU2) (FFM%NB(INTYP0), INTYP0=NCOL*NROW+1,FF%MTYP)
                   ENDIF
                   ! Type of the method to scale the energy data
                   WRITE(FF%IU2) FFM%ISCALE_TOTEN
                   ! Average energy per atom (eV atom^-1)
                   WRITE(FF%IU2) FFM%TOTENAV*EUNIT
                ENDIF
! Merge FFM%WMAT arrays if necessary
                DO INTYP0=1, FF%MTYP
#ifdef scaLAPACK
                   !Regression coefficients (eV atom^-1)
                   ! For XXX
                   CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM,FFM%WMAT(1:,1:,INTYP0), &
                        PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%NB(INTYP0),1)
                   IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      NCOL=FFM%NB(INTYP0)/NROW
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IB,1)*EUNIT,IB=IROW+1,IROW+NROW)
                      ENDDO
                      IF(MOD(FFM%NB(INTYP0),NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IB,1)*EUNIT,IB=NCOL*NROW+1,FFM%NB(INTYP0))
                      ENDIF
                   ENDIF
                   CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                   IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      NCOL=FFM%NB(INTYP0)/NROW
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (FFM%WMAT(IB,1,INTYP0)*EUNIT,IB=IROW+1,IROW+NROW)
                      ENDDO
                      IF(MOD(FFM%NB(INTYP0),NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (FFM%WMAT(IB,1,INTYP0)*EUNIT,IB=NCOL*NROW+1,FFM%NB(INTYP0))
                      ENDIF
                   ENDIF
#endif
                ENDDO
                IF(FFM%W1.GT.0.0_q) THEN
                   NCOL=FFM%MMVAR1/NROW
                   DO INTYP0=1, FF%MTYP
                     ! C00 (radial descriptor) (-)
                      ! For XXX
! Merge FFM%C00 arrays if necessary
#ifdef scaLAPACK
                      CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM,FFM%C00(1:,1:,INTYP0), &
                           PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%MMVAR1,FFM%NB(INTYP0))
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%MMVAR1,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%MMVAR1)
                            ENDIF
                         ENDDO
                      ENDIF
                      CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (FFM%C00(IVAR,IB,INTYP0),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%MMVAR1,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (FFM%C00(IVAR,IB,INTYP0),IVAR=NCOL*NROW+1,FFM%MMVAR1)
                            ENDIF
                         ENDDO
                      ENDIF
#endif
                   ENDDO
                ENDIF
! Output the angular descriptors.
! PS_SIC is not necessary because PS includes it.
                IF(FFM%W2.GT.0.0_q) THEN
                   DO INTYP0=1, FF%MTYP
                      NCOL=FFM%NNVAR2(INTYP0)/NROW
                      ! PS (angular descriptor) (-)
                      ! For XXX
!Merge FFM%PS if necessary
#ifdef scaLAPACK
                      CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM,FFM%PS(1:,1:,INTYP0), &
                           PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%NNVAR2(INTYP0),FFM%NB(INTYP0))
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%NNVAR2(INTYP0),NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%NNVAR2(INTYP0))
                            ENDIF
                         ENDDO
                      ENDIF
                      CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (FFM%PS(IVAR,IB,INTYP0),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%NNVAR2(INTYP0),NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (FFM%PS(IVAR,IB,INTYP0),IVAR=NCOL*NROW+1,FFM%NNVAR2(INTYP0))
                            ENDIF
                         ENDDO
                      ENDIF
#endif
                   ENDDO
                   IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                      DO INTYP0=1, FF%MTYP
                         NCOL=FFM%RV%NNVAR2_RED(INTYP0)/NROW
                         ! PS (angular descriptor) (-)
#ifdef scaLAPACK
                         CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM_RED,FFM%RV%PS_RED(1:,1:,INTYP0), &
                              PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%RV%NNVAR2_RED(INTYP0),FFM%NB(INTYP0))
                         IF (PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                            DO IB=1, FFM%NB(INTYP0)
                               DO ICOL=1, NCOL, 1
                                  IROW=NROW*(ICOL-1)
                                  WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                               ENDDO
                               IF (MOD(FFM%RV%NNVAR2_RED(INTYP0),NROW) .GT. 0) THEN
                                  WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%RV%NNVAR2_RED(INTYP0))
                               ENDIF
                            ENDDO
                         ENDIF
                         CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                         IF (PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                            DO IB=1, FFM%NB(INTYP0)
                               DO ICOL=1, NCOL, 1
                                  IROW=NROW*(ICOL-1)
                                  WRITE(FF%IU2) (FFM%RV%PS_RED(IVAR,IB,INTYP0),IVAR=IROW+1,IROW+NROW)
                               ENDDO
                               IF (MOD(FFM%RV%MMVAR2_RED,NROW) .GT. 0) THEN
                                  WRITE(FF%IU2) (FFM%RV%PS_RED(IVAR,IB,INTYP0),IVAR=NCOL*NROW+1,FFM%RV%NNVAR2_RED(INTYP0))
                               ENDIF
                            ENDDO
                         ENDIF
#endif
                      ENDDO
                   ENDIF
                ENDIF
                IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                   ! Normalized noise parameter (-)
                   WRITE(FF%IU2) FFM%SIGV
                   ! Normalized uncertainity parameter (-)
                   WRITE(FF%IU2) FFM%SIGW
! In this output, we assume that SIGtot/sqrt(NIONS) is constant, where
! SIGtot denotes the variance of the total energy of the systm (not per
! atom).
! This assumption is on the basis of the fact that the variance of the total energy
! should be proportional to sqrt(NIONS)*kBT. This fact indicates that
! SIGtot/sqrt(NIONS)=kBT=const. Then, const.=
! (SIGatom*NIONS)/sqrt(NIONS)=SIGatom*sqrt(NIONS), where SIGatom
! corresponds to the FFM%SIG.
! But I believe that this treatment should be done only for energy.
! For force and stress tensors, I believe that their variances do not
! depend on the system size. Hence, I do not scale them.
                   ! Variance in energies per atom in training data (eV atom^-1)
                   WRITE(FF%IU2) FFM%SIG(1)*EUNIT
                   ! Variance in forces in training data (eV Angst^-1)
                   WRITE(FF%IU2) FFM%SIG(2)*FUNIT,FFM%SIG(3)*FUNIT,FFM%SIG(4)*FUNIT
                   ! Variance in stress tensors in training data (kB)
                   WRITE(FF%IU2) FFM%SIG(5)*SUNIT,FFM%SIG(6)*SUNIT,FFM%SIG(7)*SUNIT 
                   WRITE(FF%IU2) FFM%SIG(8)*SUNIT,FFM%SIG(9)*SUNIT,FFM%SIG(10)*SUNIT
                   WRITE(FF%IU2) FF%LFAST
                   WRITE(FF%IU2) FFM%WRITECMAT
                ENDIF
                ! Control variable FF%LFAST is inserted here so it does only need
                ! to be read twice
! Writing the last bit of data depending on fitting mode
                ! Inverse of Kernel matrix 
                IF (.NOT.FFM%WRITECMAT) THEN
                   IF (PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      DO INTYP0 = 1, FF%MTYP
                         NCOL = FFM%NB(INTYP0) / NROW
                         DO IB = 1, FFM%NB(INTYP0)
#ifdef scaLAPACK
                            DO ICOL = 1, NCOL
                               IROW = NROW * (ICOL-1)
                               WRITE(FF%IU2) (FFM%SINV_ALL(IVAR,IB,INTYP0), IVAR=IROW+1, IROW+NROW)
                            ENDDO
                            IF (MOD(FFM%NB(INTYP0),NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (FFM%SINV_ALL(IVAR,IB,INTYP0), IVAR=NCOL*NROW+1, FFM%NB(INTYP0))
                            ENDIF
#else
                            DO ICOL = 1, NCOL
                               IROW = NROW * (ICOL-1)
                               WRITE(FF%IU2) (FFM%SINV(IVAR,IB,INTYP0), IVAR=IROW+1, IROW+NROW)
                            ENDDO
                            IF (MOD(FFM%NB(INTYP0),NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (FFM%SINV(IVAR,IB,INTYP0), IVAR=NCOL*NROW+1, FFM%NB(INTYP0))
                            ENDIF
#endif
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE
                ! Covariance matrix (-)
                   NCOL=FFM%NBTOT/NROW
#ifdef scaLAPACK
                   CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM,FFM%CMAT(1:,1:), &
                        PAR_SUP_HANDLE%CONTEXT_WORLD,FFM%NBTOT,FFM%NBTOT)
                   IF (PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      DO IB=1, FFM%NBTOT
                         DO ICOL=1, NCOL, 1
                            IROW=NROW*(ICOL-1)
                            WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                         ENDDO
                         IF (MOD(FFM%NBTOT,NROW) .GT. 0) THEN
                            WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%NBTOT)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                   IF (PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      DO IB=1, FFM%NBTOT
                         DO ICOL=1, NCOL, 1
                            IROW=NROW*(ICOL-1)
                            WRITE(FF%IU2) (FFM%CMAT(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                         ENDDO
                         IF (MOD(FFM%NBTOT,NROW) .GT. 0) THEN
                            WRITE(FF%IU2) (FFM%CMAT(IVAR,IB),IVAR=NCOL*NROW+1,FFM%NBTOT)
                         ENDIF
                      ENDDO
                   ENDIF
#endif
                ENDIF
             ENDIF
             PROFILING_STOP('out_ffncar')
          ENDIF
! Close file
          IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) CLOSE(FF%IU2)
        END SUBROUTINE OUT_FFNCAR

        !======================================================================
        !> Write ASCII header to ML_FFN file.
        !======================================================================
        SUBROUTINE WRITE_FFN_HEADER(BINARY, FF, FFM, NCONF, DIR_APP, DIR_LEN)

          USE ML_FF_CONSTANT
          USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
          LOGICAL          , INTENT(IN) :: BINARY
          TYPE (FF_PAR)    , INTENT(IN) :: FF
          TYPE (FFM_PAR)   , INTENT(IN) :: FFM
          INTEGER          , INTENT(IN) :: NCONF !< Number of training structure.
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
          INTEGER          , INTENT(IN) :: DIR_LEN

          ! Local variables.
          INTEGER            :: BYTES_LEFT
          INTEGER            :: IERR
          INTEGER            :: I
          INTEGER            :: ITMP
          REAL(q)            :: DTMP
          CHARACTER(LEN=255) :: STMP
          LOGICAL            :: IS_TRUNCATED

          BYTES_LEFT = ML_FF_HEADER_SIZE
          STMP = ""
          IS_TRUNCATED = .FALSE.

          ! Open file in ASCII mode.
          OPEN(FF%IU2,                             &
               FILE=DIR_APP(1:DIR_LEN)//FF%FFNCAR, &
               STATUS='REPLACE',                   &
               FORM='FORMATTED',                   &
               ACCESS='STREAM',                    &
               IOSTAT=IERR)
          IF (IERR.NE.0) THEN
             CALL ml_tutor%error("Cannot open ML_FFN file for writing.")
          ENDIF

          STMP = 'ML_FF'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          WRITE(STMP, '(1X,I0,("."),I0,("."),I0)') ML_FF_VERSION_WRITE
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          IF (BINARY) THEN
             WRITE(STMP, '(1X,A)') "binary"
          ELSE
             WRITE(STMP, '(1X,A)') "txt"
          END IF
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ' { "date" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          STMP = ' "' // FF%ML_FF_CREATION_DATETIME // '"'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_LFAST" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          
          STMP = " false"
          IF (FF%LFAST) STMP = " true"
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_DESC_TYPE" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          WRITE(STMP, '(1X,I3)') FFM%DESC_TYPE
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "types" : ['
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          DO I = 1, FF%MTYP
             WRITE(STMP, '(1X,A)') '"' // TRIM(FF%TYPE(I)) // '"'
             IF (I > 1) STMP = ',' // STMP(1:LEN(STMP)-1)
             CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          END DO

          STMP = ' ], "training_structures" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,I0)') NCONF
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "local_reference_cfgs" : ['
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          DO I = 1, FF%MTYP
             WRITE(STMP, '(1X,I0)') FFM%NB(I)
             IF (I > 1) STMP = ',' // STMP(1:LEN(STMP)-1)
             CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          END DO

          STMP = ' ], "descriptors" : ['
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          DO I = 1, FF%MTYP
             ITMP = 0
             IF (FFM%W1 > 0.0_q) ITMP = ITMP + FFM%NNVAR1
             IF (FFM%W2 > 0.0_q) ITMP = ITMP + FFM%NNVAR2(I)
             WRITE(STMP, '(1X,I0)') ITMP
             IF (I > 1) STMP = ',' // STMP(1:LEN(STMP)-1)
             CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          END DO

          STMP = ' ], "ML_IALGO_LINREG" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,I0)') FF%IALGO_LINREG
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_RCUT1" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,ES11.4)') FFM%RCUT1 * AUTOA
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_RCUT2" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,ES11.4)') FFM%RCUT2 * AUTOA
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_W1" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,ES11.4)') FFM%W1
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_SION1" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,ES11.4)') FFM%SION1 * AUTOA
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_SION2" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,ES11.4)') FFM%SION2 * AUTOA
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_LMAX2" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,I0)') FFM%LMAX2
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_MRB1" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,I0)') FFM%MRB1
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_MRB2" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,I0)') FFM%MRB2
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_IWEIGHT" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          WRITE(STMP, '(1X,I0)') FF%IWEIGHT
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_WTOTEN" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          DTMP = 1.0; IF(FF%IWEIGHT.EQ.1) DTMP = EUNIT
          WRITE(STMP, '(1X,ES11.4)') FF%WTOTEN * DTMP
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_WTIFOR" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          DTMP = 1.0; IF(FF%IWEIGHT.EQ.1) DTMP = FUNIT
          WRITE(STMP, '(1X,ES11.4)') FF%WTIFOR * DTMP
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ', "ML_WTSIF" :'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)
          DTMP = 1.0; IF(FF%IWEIGHT.EQ.1) DTMP = SUNIT
          WRITE(STMP, '(1X,ES11.4)') FF%WTSIF * DTMP
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          STMP = ' }'
          CALL WRITE_FFN_HEADER_ITEM(STMP, FF%IU2, BYTES_LEFT, IS_TRUNCATED)

          ! Fill the remaining bytes with spaces.
          DO I = BYTES_LEFT, 2, -1
             WRITE(FF%IU2, '(A)', ADVANCE='NO') " "
          END DO
          ! Finally write newline and close file.
          WRITE(FF%IU2, '(A)') ""
          CLOSE(FF%IU2)

        END SUBROUTINE WRITE_FFN_HEADER

        !======================================================================
        !> Writes a single item to ML_FFN header.
        !>
        !> If size of header would be exceeded, does not write input string.
        !> Instead sets IS_TRUNCATED to .TRUE. which prevents further output.
        !======================================================================
        SUBROUTINE WRITE_FFN_HEADER_ITEM(STRING, FILE_UNIT, BYTES_LEFT, IS_TRUNCATED)

          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
          CHARACTER(LEN=*), INTENT(IN)    :: STRING
          INTEGER         , INTENT(IN)    :: FILE_UNIT
          INTEGER         , INTENT(INOUT) :: BYTES_LEFT
          LOGICAL         , INTENT(INOUT) :: IS_TRUNCATED

          ! Local variables
          INTEGER :: STRING_LENGTH

          STRING_LENGTH = LEN(TRIM(STRING))

          ! If header was already truncated before, return immediately.
          IF (IS_TRUNCATED) RETURN

          ! Check if writing input string would exceed header size.
          IF (BYTES_LEFT - STRING_LENGTH < 2) THEN
             CALL ml_tutor%alert("ML_FFN header size is too small to collect &
                                 &all information, header is truncated.")
             IS_TRUNCATED = .TRUE.
          ! There is enough space left in header, write string now.
          ELSE
             BYTES_LEFT = BYTES_LEFT - STRING_LENGTH
             WRITE(FILE_UNIT, '(A)', ADVANCE='NO') TRIM(STRING)
          END IF

        END SUBROUTINE WRITE_FFN_HEADER_ITEM

        !======================================================================
        !> Returns string containing the current date and time in the format
        !> "2023-01-10T22:17:13.021".
        !======================================================================
        FUNCTION NOW()

           CHARACTER(LEN=23) :: NOW
           CHARACTER(LEN=8)  :: DATE
           CHARACTER(LEN=10) :: TIME

           CALL DATE_AND_TIME(DATE, TIME)

           NOW = DATE(1:4) // "-" // DATE(5:6) // "-" // DATE(7:8) // "T" // &
                 TIME(1:2) // ":" // TIME(3:4) // ":" // TIME(5:10)

        END FUNCTION NOW

#endif


#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output histogram
!****************************************************************************************************

        SUBROUTINE OUT_HISTOGRAM (COMM_WORLD,AV,DIR_APP,DIR_LEN,HISCAR,IU4,NX,SIG,X,Y)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
          USE ml_ff_tutor, ONLY: ml_tutor
          IMPLICIT NONE
! Input variables
          TYPE (ML_MPI_PAR)                :: COMM_WORLD
          REAL(q)          , INTENT(IN) :: AV(:) !(1:10)
          CHARACTER(LEN=10), INTENT(IN) ::  DIR_APP
          INTEGER          , INTENT(IN) ::  DIR_LEN
          CHARACTER (LEN=10), INTENT(IN) :: HISCAR
          INTEGER          , INTENT(IN) :: IU4
          INTEGER          , INTENT(IN) :: NX
          REAL(q)          , INTENT(IN) :: SIG(:) !(1:10)
          REAL(q)          , INTENT(IN) :: X(:) !(1:(NX+1))
          REAL(q)          , INTENT(IN) :: Y(:,:) !(1:(NX+1),1:10)
! Local variables
          INTEGER                       :: I
          INTEGER                       :: IERR
          INTEGER                       :: IX
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
! Open file
             OPEN(UNIT=IU4,FILE=DIR_APP(1:DIR_LEN)// &
                 HISCAR,STATUS='UNKNOWN',IOSTAT=IERR)
! Check if error occured in opening file
             IF (IERR.NE.0) THEN
                CALL ml_tutor%error("Cannot open ML_HIS file for writing.")
             ENDIF
! Output the calculated histograms
             WRITE(IU4,1) AV(1)*EUNIT
             WRITE(IU4,2) AV(2)*FUNIT,AV(3)*FUNIT,AV(4)*FUNIT
             WRITE(IU4,3) AV(5)*SUNIT/10.0_q,AV(8)*SUNIT/10.0_q,AV(10)*SUNIT/10.0_q, &
                          AV(6)*SUNIT/10.0_q,AV(7)*SUNIT/10.0_q,AV(9) *SUNIT/10.0_q
             WRITE(IU4,4) SIG(1)*EUNIT
             WRITE(IU4,5) SIG(2)*FUNIT,SIG(3)*FUNIT,SIG(4)*FUNIT
             WRITE(IU4,6) SIG(5)*SUNIT/10.0_q,SIG(8)*SUNIT/10.0_q,SIG(10)*SUNIT/10.0_q, &
                          SIG(6)*SUNIT/10.0_q,SIG(7)*SUNIT/10.0_q,SIG(9) *SUNIT/10.0_q
             WRITE(IU4,7)
             DO IX=1, NX
                WRITE(IU4,8) X(IX),Y(IX,1),(Y(IX,I),I=2, 4),Y(IX,5), &
                      Y(IX,8),Y(IX,10),Y(IX,6),Y(IX,7),Y(IX,9)
             ENDDO
! Close file
             CLOSE(IU4)
          ENDIF
! Format
1         FORMAT(100("*")/5X,'Average of energy (eV atom^-1)' &
                              /100("-")/5X,F16.8)
2         FORMAT(100("*")/5X,'Average of force (X Y Z) (eV ang^-1 &
                              atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8)
3         FORMAT(100("*")/5X,'Average of stress (XX YY ZZ XY XZ YZ) &
                        (GPa atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)
4         FORMAT(100("*")/5X,'Variance of energy (eV atom^-1)' &
                              /100("-")/5X,F16.8)
5         FORMAT(100("*")/5X,'Variance of force (X Y Z) (eV ang^-1 &
                              atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8)
6         FORMAT(100("*")/5X,'Variance of stress (XX YY ZZ XY XZ YZ) &
                  (GPa atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)
7         FORMAT(100("*")/5X,'Normalized histogram'/ &
                          5X,'X, HIS_E, HIS_F, HIS_S'/100("-"))
8         FORMAT(F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8, &
                F16.8,F16.8,F16.8,F16.8)
        END SUBROUTINE OUT_HISTOGRAM
#endif

!**********************************************************************
! Output local atomic energies
!***********************************************************************

        SUBROUTINE OUT_EATOM (COMM_WORLD,A,EATOM_KIN,EATOM_POT,IU,LADD,LNIONS,MIONS,MITYP,MTYP,NIONS,NITYP,NSTEP,NTYP,POSION,TYPE)
           USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
           IMPLICIT NONE
! Global variables.
           TYPE (ML_MPI_PAR)               :: COMM_WORLD
           REAL(q)         , INTENT(IN) :: A(:,:) !(1:3,1:3)
           REAL(q)         , INTENT(IN) :: EATOM_KIN(:) !(1:MIONS)
           REAL(q)         , INTENT(IN) :: EATOM_POT(:) !(1:MIONS)
           INTEGER         , INTENT(IN) :: IU
           INTEGER         , INTENT(IN) :: LADD(:) !(1:MTYP)
           INTEGER         , INTENT(IN) :: LNIONS(:,:) !(1:MITYP,1:MTYP)
           INTEGER         , INTENT(IN) :: MIONS
           INTEGER         , INTENT(IN) :: MITYP
           INTEGER         , INTENT(IN) :: MTYP
           INTEGER         , INTENT(IN) :: NIONS
           INTEGER         , INTENT(IN) :: NITYP(:) !(1:MTYP)
           INTEGER         , INTENT(IN) :: NSTEP
           INTEGER         , INTENT(IN) :: NTYP
           REAL(q)         , INTENT(IN) :: POSION(:,:) !(1:3,1:MIONS)
           CHARACTER(LEN=2), INTENT(IN) :: TYPE(:) !(1:MTYP)
! Local variables.
           INTEGER                      :: INIONS
           INTEGER                      :: INITYP
           INTEGER                      :: INTYP
           INTEGER                      :: INTYP0
           INTEGER                      :: IXYZ
           INTEGER                      :: JXYZ
           IF(COMM_WORLD%NODE_ME.EQ.0) THEN
              WRITE(IU,1) NSTEP
              WRITE(IU,2)
              DO IXYZ=1, 3
                 WRITE(IU,3) (A(JXYZ,IXYZ)*AUTOA,JXYZ=1, 3)
              ENDDO
              WRITE(IU,4)
              DO INTYP=1, NTYP
                 INTYP0=LADD(INTYP)
                 DO INITYP=1, NITYP(INTYP0)
                    INIONS=LNIONS(INITYP,INTYP0)
                    WRITE(IU,5) TYPE(INTYP0),(POSION(IXYZ,INIONS)*AUTOA,IXYZ=1, 3),EATOM_KIN(INIONS)*EUNIT,EATOM_POT(INIONS)*EUNIT
                 ENDDO
              ENDDO
           ENDIF
1          FORMAT (100("*")/5X,'NSTEP=',1X,I10)
2          FORMAT (100("=")/5X,'Lattice constant (Angst)'/100("-"))
3          FORMAT (5X,F16.8,1X,F16.8,1X,F16.8)
4          FORMAT (100("=")/5X,'Atomic positions (Angst) and its & 
                  energy (eV)'/100("-"))
5          FORMAT (5X,A2,1X,F16.8,F16.8,F16.8,F16.8,F16.8)
        END SUBROUTINE OUT_EATOM

!****************************************************************************************************
! Output the heat flux.
!****************************************************************************************************

        SUBROUTINE OUT_HEAT (COMM_WORLD,IU,NSTEP,QHEAT)
           USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
           IMPLICIT NONE
! Input variables.
           TYPE (ML_MPI_PAR)               :: COMM_WORLD
           INTEGER         , INTENT(IN) :: IU
           INTEGER         , INTENT(IN) :: NSTEP
           REAL(q)         , INTENT(IN) :: QHEAT(:) !(1:3)
! Local variables.
           INTEGER                      :: IXYZ
           IF(COMM_WORLD%NODE_ME.EQ.0) THEN
              WRITE(IU,1) NSTEP,(QHEAT(IXYZ)*EUNIT*AUTOA/TUNIT,IXYZ=1, 3)
           ENDIF
1          FORMAT (1X,'NSTEP=',I10,1X,'QXYZ=',E16.8e2,E16.8e2,E16.8e2)
        END SUBROUTINE OUT_HEAT

      END MODULE IOHANDLE
#endif
