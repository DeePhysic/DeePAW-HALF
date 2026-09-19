#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
!****************************************************************************************************
!> Module containing subroutines and functions necessary for output information on display.
!>
!> For version history see documentation of ML_LOGFILE_VERSION in ml_ff_constants.F
!****************************************************************************************************
      MODULE LOGFILE
        USE ML_FF_CONSTANT
        USE MPI_DATA
        USE ML_FF_PREC
        USE ML_FF_STRUCT, ONLY: ML_IO_WRITE

        IMPLICIT NONE

        !**********************************************************************
        !> "MSG" log line types.
        !**********************************************************************
        ENUM, BIND(C)
           !> Dummy state for usage as kind function argument.
           ENUMERATOR :: MSG_TYPE    = 0
           !> General information.
           ENUMERATOR :: MSG_INFO    = 1
           !> Warning message.
           ENUMERATOR :: MSG_WARNING = 2
        END ENUM

        !> Width of first column keyword.
        INTEGER, PARAMETER :: C1W = 15

        PRIVATE :: STRLJ,             &
                   C1W,               &
                   PRINT_TAG,         &
                   PRINT_TAG_REAL,    &
                   PRINT_TAG_INTEGER, &
                   PRINT_TAG_LOGICAL

        INTERFACE PRINT_TAG
           PROCEDURE PRINT_TAG_REAL
           PROCEDURE PRINT_TAG_INTEGER
           PROCEDURE PRINT_TAG_LOGICAL
           PROCEDURE PRINT_TAG_CHARACTER
        END INTERFACE PRINT_TAG

        CONTAINS

         SUBROUTINE PRINT_TAG_REAL(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            REAL(q),          INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",ES13.5,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),            &
               VAL,                                 &
               TAG_LIST%SYMBOL(TAG),                &
               TAG

         END SUBROUTINE PRINT_TAG_REAL

         SUBROUTINE PRINT_TAG_INTEGER(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            INTEGER,          INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",I13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),         &
               VAL,                              &
               TAG_LIST%SYMBOL(TAG),             &
               TAG

         END SUBROUTINE PRINT_TAG_INTEGER

         SUBROUTINE PRINT_TAG_LOGICAL(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            LOGICAL,          INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",L13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),         &
               VAL,                              &
               TAG_LIST%SYMBOL(TAG),             &
               TAG

         END SUBROUTINE PRINT_TAG_LOGICAL

         SUBROUTINE PRINT_TAG_CHARACTER(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            CHARACTER(LEN=*), INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",A13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),         &
               VAL,                              &
               TAG_LIST%SYMBOL(TAG),             &
               TAG

         END SUBROUTINE PRINT_TAG_CHARACTER

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "STATUS" line (+ preceding separator line).
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_STATUS(FF, LOGF, LDO_AB_INITIO)

            USE ML_FF_STRUCT, ONLY: FF_PAR, ML_MPI_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            LOGICAL,          INTENT(IN) :: LDO_AB_INITIO !< Ab initio T/F?

            ENUM, BIND(C)
               ENUMERATOR :: UNKNOWN    = -1
               ENUMERATOR :: STATE_ENUM = 0
               ENUMERATOR :: ACCURATE   = 1
               ENUMERATOR :: THRESHOLD  = 2
               ENUMERATOR :: LEARNING   = 3
               ENUMERATOR :: CRITICAL   = 4
               ENUMERATOR :: PREDICT    = 5
               ENUMERATOR :: PREDFAST   = 6
               ENUMERATOR :: INTERVAL   = 7
            END ENUM

            CHARACTER(LEN=*), PARAMETER :: STATE_STRINGS(UNKNOWN:INTERVAL) = [&
               "  unknown", &
               "?????????", &
               " accurate", &
               "threshold", &
               " learning", &
               " critical", &
               "  predict", &
               " predfast", &
               " interval" ]

            INTEGER(KIND(STATE_ENUM)) :: STATE !< State descriptor.

            STATE = STATE_ENUM

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! Separator line.
            CALL LOGLINE_SEPARATOR(LOGF)

            ! Decide current state.
            IF (FF%LMLONLY) THEN
               STATE = PREDICT
            ELSE
               IF (LDO_AB_INITIO) THEN
                  IF (FF%LDOUB) THEN
                     STATE = CRITICAL
                  ELSE
                     IF (FF%LGENFF) THEN
                        STATE = LEARNING
                     ELSE
                        IF (FF%LSAMPLE) THEN
                           STATE = THRESHOLD
                        ELSE
                           STATE = UNKNOWN ! This should never happen.
                        END IF
                     END IF
                  END IF
               ELSE
                  IF (FF%LTHRESHOLD) THEN
                     STATE = INTERVAL
                  ELSE
                     STATE = ACCURATE
                  END IF
               END IF
            END IF

            ! STATUS line.
            WRITE(LOGF%IU, '(A,1X,I9,1X,A,1X,I2,1X,L6,1X,L6,1X,I9,1X,I9)') &
               STRLJ("STATUS", C1W),                                       &
               FF%NSTEP,                                                   &
               STATE_STRINGS(STATE),                                       &
               STATE,                                                      &
               LDO_AB_INITIO,                                              &
               FF%IFF,                                                     &
               FF%NSTEP_FROM_SAMPLING,                                     &
               FF%NSTEP_FROM_GENFF

         END SUBROUTINE LOGLINE_STATUS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "LCONF" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_LCONF(FF, FFM, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! LCONF line.
            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ("LCONF", C1W),                   &
               FF%NSTEP
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU, '(1X,A2,2(1X,I9))', ADVANCE='NO') &
                  FF%TYPE(FF%LADD(I)),                          &
                  FFM%NB(FF%LADD(I)),                           &
                  FFM%NB(FF%LADD(I)) + FFM%NB_NEW(FF%LADD(I))
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_LCONF

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "STDAB" and "ERR" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_STDAB_ERR(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! STDAB line.
            ASSOCIATE(S => FF%SIG)
            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))')                &
               STRLJ("STDAB", C1W),                                 &
               FF%NSTEP,                                            &
               S(1) * EUNIT,                                        &
               SQRT((S(2)**2 + S(3)**2 + S(4)**2) / 3.0_q) * FUNIT, &
               SQRT((S(5)**2 + S(6)**2 + S(7)**2 +                  &
                     S(8)**2 + S(9)**2 + S(10)**2) / 6.0_q) * SUNIT
            END ASSOCIATE
            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))') &
               STRLJ("ERR", C1W),                    &
               FF%NSTEP,                             &
               FF%STOTEN * EUNIT,                    &
               FF%STIFOR * FUNIT,                    &
               FF%STSIF * SUNIT

         END SUBROUTINE LOGLINE_STDAB_ERR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "NORME" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_NORME(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            WRITE(LOGF%IU, '(A,1X,I9,3(1X,F16.3))') &
               STRLJ("NORME", C1W),                 &
               FF%NSTEP,                            &
               FF%RELTOTEN,                         &
               FF%RELTIFOR,                         &
               FF%RELTSIF

         END SUBROUTINE LOGLINE_NORME

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "ERPS" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_ERPS(FF, ABN, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, ABINITIONEW, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN)  :: FF !< FF setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            TYPE(LOG_PAR),    INTENT(IN)  :: LOGF !< Log file setup.
            ! Local variables
            INTEGER                       :: INTYP

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ("ERPS", C1W),                    &
               FF%NSTEP
            DO INTYP = 1, FF%MTYP_ALLOCATE
               IF (INTYP > FF%MTYP) THEN
                  IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,1(1X,A16))', ADVANCE='NO') &
                      FF%TYPE(INTYP), "NaN"
               ELSE
                  IF (ABN%NITYP(INTYP) <= 0) THEN
                     IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,1(1X,A16))', ADVANCE='NO') &
                        FF%TYPE(INTYP), "NaN"
                  ELSE
                     WRITE(LOGF%IU, '(1X,A2,1(1X,ES16.8))', ADVANCE='NO') &
                     FF%TYPE(INTYP),                                       &
                     FF%STIFOR_TYPE(INTYP) * FUNIT
                  ENDIF
               ENDIF
            ENDDO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_ERPS


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "NEPS" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_NEPS(FF, ABN, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, ABINITIONEW, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN)  :: FF !< FF setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            TYPE(LOG_PAR),    INTENT(IN)  :: LOGF !< Log file setup.
            ! Local variables
            INTEGER                       :: INTYP

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ("NEPS", C1W),                    &
               FF%NSTEP

            DO INTYP = 1, FF%MTYP_ALLOCATE
               IF (INTYP > FF%MTYP) THEN
                  IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,1(1X,A16))', ADVANCE='NO') &
                      FF%TYPE(INTYP), "NaN"
               ELSE
                  IF (ABN%NITYP(INTYP) <= 0) THEN
                     IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,1(1X,A16))', ADVANCE='NO') &
                        FF%TYPE(INTYP), "NaN"
                  ELSE
                     WRITE(LOGF%IU, '(1X,A2,1(1X,F16.3))', ADVANCE='NO') &
                     FF%TYPE(INTYP),                                      &
                     FF%RELTIFOR_TYPE(INTYP)
                  ENDIF
               ENDIF
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_NEPS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "LASTE" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_LAST_ERR(FF, LOGF)
            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))') &
               STRLJ("LASTE", C1W),                  &
               FF%NSTEP,                             &
               FF%STOTEN_LAST * EUNIT,               &
               FF%STIFOR_LAST * FUNIT,               &
               FF%STSIF_LAST  * SUNIT
           
         END SUBROUTINE LOGLINE_LAST_ERR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "LASTNE" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_LASTNE(FF, LOGF)
            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            WRITE(LOGF%IU, '(A,1X,I9,3(1X,F16.3))') &
               STRLJ("LASTNE", C1W),                &
               FF%NSTEP,                            &
               FF%RELTOTEN_LAST,                    &
               FF%RELTIFOR_LAST,                    &
               FF%RELTSIF_LAST  
           
         END SUBROUTINE LOGLINE_LASTNE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "BEE" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_BEE(FF, LOGF, ABN, FINAL_CALL)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(FF_PAR),      INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            !> Whether this is the final call in this time step.
            LOGICAL,           INTENT(IN) :: FINAL_CALL

            CHARACTER(LEN=:), ALLOCATABLE :: KEYWORD

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            KEYWORD = "BEE"
            IF (FINAL_CALL) KEYWORD = "BEEF"

            WRITE(LOGF%IU, '(A,1X,I9,6(1X,ES16.8))') &
               STRLJ(KEYWORD, C1W),                  &
               FF%NSTEP,                             &
               ABN%BEMAX(1, ABN%NCONF) * EUNIT,      &
               ABN%BEMAX(2, ABN%NCONF) * FUNIT,      &
               ABN%BEAV(2, ABN%NCONF) * FUNIT,       &
               FF%CTIFOR * FUNIT,                    &
               ABN%BEMAX(3, ABN%NCONF) * SUNIT,      &
               ABN%BEAV(3, ABN%NCONF) * SUNIT

         END SUBROUTINE LOGLINE_BEE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "BEPS" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_BEPS(FF, LOGF, ABN, FINAL_CALL)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(FF_PAR),      INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            !> Whether this is the final call in this time step.
            LOGICAL,           INTENT(IN) :: FINAL_CALL

            ! Local variables
            CHARACTER(LEN=:), ALLOCATABLE :: KEYWORD
            INTEGER                       :: INTYP

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            KEYWORD = "BEPS"
            IF (FINAL_CALL) KEYWORD = "BEFPS"

            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ(KEYWORD, C1W),                   &
               FF%NSTEP
            DO INTYP = 1, FF%MTYP_ALLOCATE
               IF (INTYP > FF%MTYP) THEN
                  IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,2(1X,A16))', ADVANCE='NO') &
                      FF%TYPE(INTYP), "NaN", "NaN"
               ELSE 
                  IF (ABN%NITYP(INTYP) <= 0) THEN
                     IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,2(1X,A16))', ADVANCE='NO') &
                        FF%TYPE(INTYP), "NaN", "NaN"
                  ELSE
                     WRITE(LOGF%IU, '(1X,A2,2(1X,ES16.8))', ADVANCE='NO') &
                           FF%TYPE(INTYP),                                   &
                     ABN%BEMAX_TYPE(2, INTYP, ABN%NCONF) * FUNIT,      &
                     ABN%BEAV_TYPE(2, INTYP, ABN%NCONF) * FUNIT
                  END IF 
               END IF
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_BEPS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "SF" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_SF(FF, LOGF, ABN, FINAL_CALL)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(FF_PAR),      INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            !> Whether this is the final call in this time step.
            LOGICAL,           INTENT(IN) :: FINAL_CALL
            CHARACTER(LEN=:), ALLOCATABLE :: KEYWORD

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            KEYWORD = "SF"
            IF (FINAL_CALL) KEYWORD = "SFF"

            WRITE(LOGF%IU, '(A,1X,I9,5(1X,ES16.8))') &
                  STRLJ(KEYWORD, C1W),               &
                  FF%NSTEP,                          &
                  ABN%SFMAX(ABN%NCONF),              &
                  ABN%SFMIN(ABN%NCONF),              &
                  ABN%SFMEAN(ABN%NCONF),             &
                  ABN%SFVAR(ABN%NCONF),              &
                  FF%CTIFOR * FUNIT

         END SUBROUTINE LOGLINE_SF

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "SPFPS" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_SPFPS(FF, LOGF, ABN, FINAL_CALL)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(FF_PAR),      INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            !> Whether this is the final call in this time step.
            LOGICAL,           INTENT(IN) :: FINAL_CALL
            CHARACTER(LEN=:), ALLOCATABLE :: KEYWORD
            INTEGER                       :: INTYP

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            KEYWORD = "SPFPS"
            IF (FINAL_CALL) KEYWORD = "SPFFPS"

            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ(KEYWORD, C1W),                   &
               FF%NSTEP
            DO INTYP = 1, FF%MTYP_ALLOCATE
               IF (INTYP > FF%MTYP) THEN
                  IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,2(1X,A16))', ADVANCE='NO') &
                      FF%TYPE(INTYP), "NaN", "NaN"
               ELSE
                  IF (ABN%NITYP(INTYP) <= 0) THEN
                     IF (FF%ISTART == 3) WRITE(LOGF%IU, '(1X,A2,2(1X,A16))', ADVANCE='NO') &
                        FF%TYPE(INTYP), "NaN", "NaN"
                  ELSE
                     WRITE(LOGF%IU, '(1X,A2,2(1X,ES16.8))', ADVANCE='NO') &
                        FF%TYPE(INTYP),                                   &
                        ABN%SFMAX_TYPE(INTYP, ABN%NCONF),                 &
                        ABN%SFMEAN_TYPE(INTYP, ABN%NCONF)
                  ENDIF
               END IF
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_SPFPS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "CFE" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_CFE(LOGF, NSTEP, FRMSE, FMAE, FMAX)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,          INTENT(IN) :: NSTEP !< Current time step.
            REAL(q),          INTENT(IN) :: FRMSE !< RMSE of forces.
            REAL(q),          INTENT(IN) :: FMAE !< MAE of forces.
            REAL(q),          INTENT(IN) :: FMAX !< Maximum error of forces.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! CFE line.
            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))') &
               STRLJ("CFE", C1W),                    &
               NSTEP,                                &
               FRMSE * FUNIT,                        &
               FMAE * FUNIT,                         &
               FMAX * FUNIT

         END SUBROUTINE LOGLINE_CFE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "THRHIST" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_THRHIST(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! THRHIST line.
            DO I = 1, FF%NHIS
               WRITE(LOGF%IU, '(A,1X,I9,1X,I4,1X,ES16.8)') &
                  STRLJ("THRHIST", C1W),                   &
                  FF%NSTEP,                                &
                  I,                                       &
                  FF%BEHIS(I) * FUNIT
            END DO

         END SUBROUTINE LOGLINE_THRHIST

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "THRUPD" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_THRUPD(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! THRUPD line.
            WRITE(LOGF%IU, '(A,1X,I9,4(1X,ES16.8))') &
               STRLJ("THRUPD", C1W),                 &
               FF%NSTEP,                             &
               FF%CTIFOR_PREV * FUNIT,               &
               FF%CTIFOR * FUNIT,                    &
               FF%CSIG_LOG,                          &
               FF%CSLOPE_LOG

         END SUBROUTINE LOGLINE_THRUPD

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "SPRSC" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_SPRSC(FF, FFM, ABN, LOGF, NB0, NCONF0, NCONF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, LOG_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,       INTENT(IN) :: NB0(:) !< LRC before spars.
            INTEGER,       INTENT(IN) :: NCONF0 !< Conf. before spars.
            INTEGER,       INTENT(IN) :: NCONF !< Conf. after spars.

            INTEGER :: INTYP

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! SPRSC line.
            WRITE(LOGF%IU, '(A,1X,I9,1x,I9,1x,I9)', ADVANCE='NO') &
               STRLJ("SPRSC", C1W),                               &
               FF%NSTEP,                                          &
               NCONF0,                                            &
               NCONF
            DO INTYP = 1, FF%MTYP
               WRITE(LOGF%IU, '(1X,A2,2(1X,I9))', ADVANCE='NO') &
               FF%TYPE(INTYP),                                   &
               NB0(INTYP),                                       &
               FFM%NB(INTYP)
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_SPRSC

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "NDESC" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_NDESC(FF, FFM, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! NDESC line.
            WRITE(LOGF%IU, '(A,1X,I9,1x,I9)', ADVANCE='NO') &
               STRLJ("NDESC", C1W),                         &
               FF%NSTEP,                                    &
               FFM%NNVAR1                                      
            IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN 
               DO I = 1, FF%MTYP
                  WRITE(LOGF%IU, '(1X,A2,1X,I9)', ADVANCE='NO') &
                     FF%TYPE(I),                                &
                     FFM%NNVAR2(I)+FFM%RV%NNVAR2_RED(I)                
               END DO
            ELSE
               DO I = 1, FF%MTYP
                  WRITE(LOGF%IU, '(1X,A2,1X,I9)', ADVANCE='NO') &
                     FF%TYPE(I),                                &
                     FFM%NNVAR2(I)
               END DO
            ENDIF
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_NDESC


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "NDESC_SIC" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_NDESC_SIC(FF, FFM, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! NDESC_SIC line.
            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ("NDESC_SIC", C1W),               &
               FF%NSTEP
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU, '(1X,A2,1X,I9)', ADVANCE='NO') &
                  FF%TYPE(I),                                &
                  FFM%NNVAR2_SIC(I)
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_NDESC_SIC

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "REGR" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_REGR(LOGF,     &
                                 NSTEP,    &
                                 LIN_REG,  &
                                 ITE,      &
                                 SIGW2,    &
                                 SIGV2,    &
                                 REGRATIO, &
                                 WRATIO,   &
                                 WMAX,     &
                                 ZTR,      &
                                 LEXIT)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,       INTENT(IN) :: NSTEP !< Current time step.
            INTEGER,       INTENT(IN) :: LIN_REG !< BLR, QR or SVD
            INTEGER,       INTENT(IN) :: ITE !< Regression loop index.
            REAL(q),       INTENT(IN) :: SIGW2 !< sigma_W^2.
            REAL(q),       INTENT(IN) :: SIGV2 !< sigma_V^2.
            REAL(q),       INTENT(IN) :: REGRATIO !< (SIGW2 / SIGV2) / WMAX.
            REAL(q),       INTENT(IN) :: ZTR !< Error of force field
            REAL(q),       INTENT(IN) :: WRATIO !< WMAX / WMIN.
            REAL(q),       INTENT(IN) :: WMAX !< highest eigenvalue/singular value
            LOGICAL,       INTENT(IN) :: LEXIT !< If last loop iteration.

            CHARACTER(LEN=C1W) :: COL1

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! Last iteration uses slightly different keyword ("F" for final).
            COL1 = "REGR"
            IF (LEXIT) COL1 = "REGRF"

            ! REGR line.
            WRITE(LOGF%IU, '(A,1X,I9,1x,I4,1X,I4,4(1X,ES16.8))', ADVANCE='NO') &
               STRLJ(COL1, C1W),                                               &
               NSTEP,                                                          &
               LIN_REG,                                                        &
               ITE,                                                            &
               SIGW2,                                                          &
               SIGV2,                                                          &
               REGRATIO,                                                       &
               ZTR

            ! Write WRATIO and WMAX only in last line because it does not
            ! change over iterations.
            IF (LEXIT) WRITE(LOGF%IU, '(1X,ES16.8)', ADVANCE='NO') WRATIO
            IF (LEXIT) WRITE(LOGF%IU, '(1X,ES16.8)', ADVANCE='NO') WMAX

            ! Newline.
            WRITE(LOGF%IU, '()')

         END SUBROUTINE LOGLINE_REGR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "DCOUPLE" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_DCOUPLE(LOGF, NSTEP, DHDLAMBDA)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,          INTENT(IN) :: NSTEP !< Current time step.
            REAL(q),          INTENT(IN) :: DHDLAMBDA !< dH/dlambda.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! DCOUPLE line.
            WRITE(LOGF%IU, '(A,1X,I9,1X,ES16.8)') &
               STRLJ("DCOUPLE", C1W),             &
               NSTEP,                             &
               DHDLAMBDA * EUNIT

         END SUBROUTINE LOGLINE_DCOUPLE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "FFOUT" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_FFOUT(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            ! Also return early if skipping of log lines is activated.
            IF (.NOT. ML_IO_WRITE .OR. LOGF%SKIP_LOG_LINES) RETURN

            ! FFOUT line.
            WRITE(LOGF%IU, '(A,1X,I9,1X,A23)') &
               STRLJ("FFOUT", C1W),            &
               FF%NSTEP,                       &
               FF%ML_FF_CREATION_DATETIME

         END SUBROUTINE LOGLINE_FFOUT

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "MSG" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_MSG(FF, LOGF, MESSAGE_TYPE, MESSAGE)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),            INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),           INTENT(IN) :: LOGF !< Log file setup.
            INTEGER(KIND(MSG_TYPE)), INTENT(IN) :: MESSAGE_TYPE
            CHARACTER(LEN=*),        INTENT(IN) :: MESSAGE

            CHARACTER(LEN=*), PARAMETER :: &
               TYPE_STRINGS(MSG_TYPE:MSG_WARNING) = ["?????????", &
                                                     "     info", &
                                                     "  warning"]

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! MSG line.
            WRITE(LOGF%IU, '(A,1X,I9,1X,A9,1X,A)') &
               STRLJ("MSG", C1W),                 &
               FF%NSTEP,                          &
               TYPE_STRINGS(MESSAGE_TYPE),        &
               MESSAGE

         END SUBROUTINE LOGLINE_MSG

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Separator for log lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_SEPARATOR(LOGF)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, '(80("-"))')

         END SUBROUTINE LOGLINE_SEPARATOR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about input parameters.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_ML_FF_PARAMETERS(LOGF,FF,FFM)
            USE ML_FF_CONSTANT
            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, FFM_PAR
            USE ml_ff_string, ONLY: str
            IMPLICIT NONE

            ! Input variables
            TYPE(LOG_PAR), INTENT(IN)    :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN)    :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN)    :: FFM !< FFM setup.
            ! Local variables
            INTEGER                      :: ITMP
            INTEGER                      :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* ML_FF INPUT FILE INFORMATION ************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            WRITE(LOGF%IU,'(A,A)')                                   &
               "ML_FF file version                              : ", &
               str(FF%ML_FF_VERSION_READ(1)) // "." //               &
               str(FF%ML_FF_VERSION_READ(2)) // "." //               &
               str(FF%ML_FF_VERSION_READ(3))

            WRITE(LOGF%IU,'(A)', ADVANCE='NO') &
               "ML_FF file mode                                 : "
            IF (FF%ML_FF_BINARY) THEN
               WRITE(LOGF%IU,'(A)') "binary"
            ELSE
               WRITE(LOGF%IU,'(A)') "text"
            ENDIF

            WRITE(LOGF%IU,'(A,L1)')                                  &
               "ML_FF header present                            : ", &
               FF%ML_FF_HAS_HEADER

            IF (FF%ML_FF_HAS_HEADER) THEN
               WRITE(LOGF%IU,'(A,A)')                                   &
                  "ML_FF creation timestamp                        : ", &
                  FF%ML_FF_CREATION_DATETIME
            ENDIF

            WRITE(LOGF%IU,'(A,L1)')                                  &
               "ML_FF supports fast prediction mode             : ", &
               FF%LFAST

            WRITE(LOGF%IU,'(A,ES16.8)')                              &
               "RMSE energies (eV atom^-1)                      : ", &
               FF%STOTEN * EUNIT

            WRITE(LOGF%IU,'(A,ES16.8)')                              &
               "RMSE forces (eV Angst^-1)                       : ", &
               FF%STIFOR * FUNIT

            WRITE(LOGF%IU,'(A,ES16.8)')                              &
               "RMSE stresses (kB)                              : ", &
               FF%STSIF * SUNIT

            WRITE(LOGF%IU,'(A,I0)')                                  &
               "Number of atom types                            : ", &
               FF%MTYP

            WRITE(LOGF%IU,'(A)', ADVANCE='NO') &
               "Atom types                                      :"
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU,'(1X,A12)', ADVANCE='NO') TRIM(ADJUSTL(FF%TYPE(I)))
            END DO
            WRITE(LOGF%IU,'(A)') ""

            WRITE(LOGF%IU,'(A)', ADVANCE='NO') &
               "Atomic reference energies                       :"
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU,'(1X,ES12.4)', ADVANCE='NO') FF%EATOM_REF(I) * EUNIT
            ENDDO
            WRITE(LOGF%IU,'(A)') ""

            WRITE(LOGF%IU,'(A)', ADVANCE='NO') &
               "Number of local reference configurations        :"
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU,'(1X,I12)', ADVANCE='NO') FFM%NB(I)
            ENDDO
            WRITE(LOGF%IU,'(A)') ""

            WRITE(LOGF%IU,'(A)', ADVANCE='NO') &
               "Number of descriptors                           :"
            DO I = 1, FF%MTYP
               ITMP = 0
               IF (FFM%W1 > 0.0_q) ITMP = ITMP + FFM%NNVAR1
               IF (FFM%W2 > 0.0_q) ITMP = ITMP + FFM%NNVAR2(I)
               WRITE(LOGF%IU,'(1X,I12)', ADVANCE='NO') ITMP
            ENDDO
            WRITE(LOGF%IU,'(A)') ""

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_ML_FF_PARAMETERS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about input parameters.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_PARAMETERS(TAG_LIST, LOGF, FF, FFM, EMP)

            USE ML_FF_CONSTANT
            USE ML_FF_TAGLIST
            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, FFM_PAR, EMP_POT_TYPE
            IMPLICIT NONE

            TYPE(TAGLIST), INTENT(INOUT) :: TAG_LIST !< List with tag occurrences.
            TYPE(LOG_PAR), INTENT(IN)    :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN)    :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN)    :: FFM !< FFM setup.
            TYPE(EMP_POT_TYPE), INTENT(IN) :: EMP

            INTEGER                      :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* MACHINE LEARNING SETTINGS ***************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            WRITE(LOGF%IU, '(A)') 'This section lists the available &
               &machine-learning related settings with a short description, &
               &their'
            WRITE(LOGF%IU, '(A)') 'selected values and the INCAR tags. The &
               &column between the value and the INCAR tag may contain a'
            WRITE(LOGF%IU, '(A)') '"state indicator" highlighting the origin &
               &of the value. Here is a list of possible indicators:'
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')                          &
               TAGLIST_STATE_SYMBOL(TL_DEFAULT),                          &
               "(empty) Tag was not provided in the INCAR file, a default &
               &value was chosen automatically."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)') &
               TAGLIST_STATE_SYMBOL(TL_INCAR),   &
               "Value was provided in the INCAR file."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')   &
               TAGLIST_STATE_SYMBOL(TL_INCAR_ALT), &
               "Value was provided in the INCAR file, deprecated tag."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')                       &
               TAGLIST_STATE_SYMBOL(TL_OVERRIDE),                      &
               "A value found in the INCAR file was overwritten by the &
               &contents of the ML_FF or ML_AB file."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')                            &
               TAGLIST_STATE_SYMBOL(TL_UNSET),                              &
               "The value for this tag was never set (please report this to &
               &the VASP developers)."
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') 'Tag values with associated units are given &
               &here in Angstrom/eV, if not specified otherwise.'
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') 'Please refer to the VASP online manual for &
               &a detailed description of available INCAR tags.'

            CALL SECTION_HEADER("General settings")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MODE                 ", FF%MODE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ISTART               ", FF%ISTART)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_TYPE                 ", FF%ML_TYPE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_OUTPUT_MODE          ", FF%OUTPUT_MODE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_OUTBLOCK             ", FF%OUTBLOCK)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LFAST                ", FF%LFAST)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LIB                  ", FF%LIB)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_PALGO                ", FF%PALGO)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NCSHMEM              ", FF%NCSHMEM)

            CALL SECTION_HEADER("Descriptor settings")
            WRITE(LOGF%IU, '(A)') "Radial descriptors:"
            WRITE(LOGF%IU, '(A)') "-------------------"
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RCUT1                ", FFM%RCUT1 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SION1                ", FFM%SION1 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB1                 ", FFM%MRB1)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "Angular descriptors:"
            WRITE(LOGF%IU, '(A)') "--------------------"
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_DESC_TYPE            ", FFM%DESC_TYPE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RCUT2                ", FFM%RCUT2 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SION2                ", FFM%SION2 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB2                 ", FFM%MRB2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LMAX2                ", FFM%LMAX2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LAFILT2              ", FFM%LAFILT2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_AFILT2               ", FFM%AFILT2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IAFILT2              ", FFM%IAFILT2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LSPARSDES            ", FFM%LSPARSDES)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NRANK_SPARSDES       ", FFM%NRANK_SPARSDES)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RDES_SPARSDES        ", FFM%RDES_SPARSDES)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LSIC                 ", FFM%LSIC)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_BASIS_TYPE2          ", FFM%BASIS_TYPE2)
            IF (FFM%BASIS_TYPE2.EQ.4 .OR. FFM%BASIS_TYPE2.EQ.5) THEN
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB2_MAX             ", FFM%MRB2_MAX)
            ENDIF
            IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LMAX2_RED             ", FFM%RV%LMAX2_RED)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB2_RED              ", FFM%RV%MRB2_RED)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_DESC_RATIO_DUAL     ", FFM%RV%RATIO_DESC)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_DESC_FACTORE_TESTE    ", FFM%RV%FACTOR_DESC)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_BASIS_TYPE2_RED       ", FFM%RV%BASIS_TYPE2_RED)
               IF (FFM%RV%BASIS_TYPE2_RED.EQ.4 .OR. FFM%RV%BASIS_TYPE2_RED.EQ.5) THEN
                  CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB2_MAX_RED          ", FFM%RV%MRB2_MAX_RED)
               ENDIF
            ENDIF

            CALL SECTION_HEADER("Kernel settings")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_W1                   ", FFM%W1)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NHYP                 ", FFM%NHYP1)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LSUPERVEC            ", FFM%LSUPERVEC)

            CALL SECTION_HEADER("Error estimation")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ICRITERIA            ", FF%ICRITERIA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IUPDATE_CRITERIA     ", FF%IUPDATE_CRITERIA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CTIFOR               ", FF%CTIFOR * FUNIT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SCLC_CTIFOR          ", FF%SCALE_CTIFOR)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CSIG                 ", FF%CSIG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CSLOPE               ", FF%CSLOPE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CX                   ", FF%CX)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CDOUB                ", FF%CDOUB)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CALGO                ", FF%CALGO)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ESTBLOCK             ", FF%IERR)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SF_LIMIT             ", FF%SF_LIMIT)

            CALL SECTION_HEADER("Sparsification and regression")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_EPS_LOW              ", FF%ML_EPS_LOW)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_EPS_REG              ", FF%ML_EPS_REG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LBASIS_DISCARD       ", FF%LBASIS_DISCARD)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IALGO_LINREG         ", FF%IALGO_LINREG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ISVD                 ", FF%ISVD)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IREG                 ", FFM%IREG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SIGV0                ", FFM%SIGV0)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SIGW0                ", FFM%SIGW0)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SAVECMAT             ", FFM%SAVECMAT)

            CALL SECTION_HEADER("Weighting and sampling")
            ! We need to take care of the different units if IWEIGHT=1
            IF (FF%IWEIGHT.EQ.1) THEN
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTOTEN               ", FF%WTOTEN * EUNIT)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTIFOR               ", FF%WTIFOR * FUNIT)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTSIF                ", FF%WTSIF * SUNIT)
            ELSE
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTOTEN               ", FF%WTOTEN)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTIFOR               ", FF%WTIFOR)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTSIF                ", FF%WTSIF)
            ENDIF
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MHIS                 ", FF%MHIS)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NMDINT               ", FF%NMDINT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IWEIGHT              ", FF%IWEIGHT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LUSE_NAMES           ", FF%LUSE_NAMES)

            CALL SECTION_HEADER("Static array sizes")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MB                   ", FFM%MB_TARGET)
            WRITE(LOGF%IU, '(A," : ",I13)')                                &
               "Maximum number of local reference configurations in memory &
               &(max. buffer size before sparsification)  ",               &
               FFM%MB
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MB_MIN               ", FFM%MB_MIN)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MCONF                ", FF%MCONF)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MCONF_NEW            ", FF%MCONF_NEW)

            CALL SECTION_HEADER("Special features")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LHEAT                ", FFM%LHEAT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LCOUPLE              ", FFM%LCOUPLE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NATOM_COUPLED        ", FFM%NATOM_COUPLED)
            WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION("ML_ICOUPLE"),     &
               SIZE(FFM%ICOUPLE),                     &
               TAG_LIST%SYMBOL("ML_ICOUPLE"),         &
               "SIZE( ML_ICOUPLE )"
            BLOCK
            CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH()) :: TMP_STRING
            DO I = 1, SIZE(FFM%ICOUPLE)
               TMP_STRING = "Ion index"
               WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A,I0,A)') &
                  ADJUSTR(TMP_STRING),                        &
                  FFM%ICOUPLE(I),                             &
                  TAG_LIST%SYMBOL("ML_ICOUPLE"),              &
                  "ML_ICOUPLE(",                              &
                  I,                                          &
                  ")"
            ENDDO
            END BLOCK
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RCOUPLE              ", FFM%RCOUPLE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LEMPPOT              ", FFM%LEMPPOT)
            IF (FFM%LEMPPOT) THEN
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_EMPPOT_RCUT          ", EMP%EMPPOT_RCUT * AUTOA)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SRPOT_B0             ", EMP%SRPOT_B0 * EUNIT)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SRPOT_S0             ", EMP%SRPOT_S0 * AUTOA)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SRPOT_N0             ", EMP%SRPOT_N0)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MOPOT_NM             ", EMP%MOPOT_NM)
               IF (EMP%MOPOT_NM.GT.0) THEN
                  WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
                     TAGLIST_DESCRIPTION("ML_MOPOT_DM"),    &
                     SIZE(EMP%MOPOT_DM),                    &
                     TAG_LIST%SYMBOL("ML_MOPOT_DM"),        &
                     "SIZE( ML_MOPOT_DM )"
                  BLOCK
                  CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH()) :: TMP_STRING
                  DO I = 1, SIZE(EMP%MOPOT_DM)
                     TMP_STRING = "ML_MOPOT_DM"
                     WRITE(LOGF%IU, '(A," : ",ES13.4,1X,A,1X,A,I0,A)') &
                        ADJUSTR(TMP_STRING),                           &
                        EMP%MOPOT_DM(I) * EUNIT,                       &
                        TAG_LIST%SYMBOL("ML_MOPOT_DM"),                &
                        "ML_MOPOT_DM(",                                &
                        I,                                             &
                        ")"
                  ENDDO
                  END BLOCK

                  WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
                     TAGLIST_DESCRIPTION("ML_MOPOT_RM"),    &
                     SIZE(EMP%MOPOT_RM),                    &
                     TAG_LIST%SYMBOL("ML_MOPOT_RM"),        &
                     "SIZE( ML_MOPOT_RM )"
                  BLOCK
                  CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH()) :: TMP_STRING
                  DO I = 1, SIZE(EMP%MOPOT_RM)
                     TMP_STRING = "ML_MOPOT_RM"
                     WRITE(LOGF%IU, '(A," : ",ES13.4,1X,A,1X,A,I0,A)') &
                        ADJUSTR(TMP_STRING),                           &
                        EMP%MOPOT_RM(I) * AUTOA,                       &
                        TAG_LIST%SYMBOL("ML_MOPOT_RM"),                &
                        "ML_MOPOT_RM(",                                &
                        I,                                             &
                        ")"
                  ENDDO
                  END BLOCK

                  WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
                     TAGLIST_DESCRIPTION("ML_MOPOT_RKM"),   &
                     SIZE(EMP%MOPOT_RKM),                   &
                     TAG_LIST%SYMBOL("ML_MOPOT_RKM"),       &
                     "SIZE( ML_MOPOT_RKM )"
                  BLOCK
                  CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH()) :: TMP_STRING
                  DO I = 1, SIZE(EMP%MOPOT_RKM)
                     TMP_STRING = "ML_MOPOT_RKM"
                     WRITE(LOGF%IU, '(A," : ",ES13.4,1X,A,1X,A,I0,A)')                    &
                        ADJUSTR(TMP_STRING),                                              &
                        EMP%MOPOT_RKM(I)**2 * 2.0_q * EMP%MOPOT_DM(I) * EUNIT / AUTOA**2, &
                        TAG_LIST%SYMBOL("ML_MOPOT_RKM"),                                  &
                        "ML_MOPOT_RKM(",                                                  &
                        I,                                                                &
                        ")"
                  ENDDO
                  END BLOCK

                  WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
                     TAGLIST_DESCRIPTION("ML_MOPOT_IJM"),   &
                     SIZE(EMP%MOPOT_IJM),                   &
                     TAG_LIST%SYMBOL("ML_MOPOT_IJM"),       &
                     "SIZE( ML_MOPOT_IJM )"
                  BLOCK
                  CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH()) :: TMP_STRING
                  DO I = 1, SIZE(EMP%MOPOT_IJM)
                     TMP_STRING = "ML_MOPOT_IJM"
                     WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A,I0,A)') &
                        ADJUSTR(TMP_STRING),                        &
                        EMP%MOPOT_IJM(I),                           &
                        TAG_LIST%SYMBOL("ML_MOPOT_IJM"),            &
                        "ML_MOPOT_IJM(",                            &
                        I,                                          &
                        ")"
                  ENDDO
                  END BLOCK
               ENDIF
            ENDIF

            CALL SECTION_HEADER("Reference energies")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ISCALE_TOTEN         ", FFM%ISCALE_TOTEN)
            WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION("ML_EATOM_REF"),   &
               SIZE(FF%EATOM_REF),                    &
               TAG_LIST%SYMBOL("ML_EATOM_REF"),       &
               "SIZE( ML_EATOM_REF )"
            BLOCK
            CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH() - 3) :: TMP_STRING
            DO I = 1, SIZE(FF%EATOM_REF)
               TMP_STRING = "Element"
               WRITE(LOGF%IU, '(A," ",A," : ",ES13.5,1X,A,1X,A,I0,A)') &
                  ADJUSTR(TMP_STRING),                                 &
                  FF%TYPE(I),                                          &
                  FF%EATOM_REF(I) * EUNIT,                             &
                  TAG_LIST%SYMBOL("ML_EATOM_REF"),                     &
                  "ML_EATOM_REF(",                                     &
                  I,                                                   &
                  ")"
            ENDDO
            END BLOCK

            CALL SECTION_HEADER("Input/output options")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LEATOM               ", FFM%LEATOM)

            ! Write input tag deprecation warnings.
            BLOCK
            !CHARACTER(LEN=:), ALLOCATABLE :: DEPRECATED(:,:)
            CHARACTER(LEN=TAGLIST_GET_TAG_WIDTH()), &
               ALLOCATABLE :: DEPRECATED(:,:)
            CALL TAG_LIST%GET_DEPRECATED(DEPRECATED)
            IF (SIZE(DEPRECATED, 1) > 0) THEN
               WRITE(LOGF%IU, *)
               WRITE(LOGF%IU, '(A)') 'WARNING: Deprecated input tags were &
                  &found, please consider the replacements:'
               DO I = 1, SIZE(DEPRECATED, 1)
                  WRITE(LOGF%IU, '(("* "),A,(" => "),A)') &
                     ADJUSTL(TRIM(DEPRECATED(I, 1))),     &
                     ADJUSTL(TRIM(DEPRECATED(I, 2)))
               END DO
            END IF
            END BLOCK

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

            CONTAINS

            SUBROUTINE SECTION_HEADER(HEADER)

               CHARACTER(LEN=*), INTENT(IN) :: HEADER

               WRITE(LOGF%IU, *)
               WRITE(LOGF%IU, *)
               WRITE(LOGF%IU, '(A)') HEADER
               WRITE(LOGF%IU, '(140("-"))')

            END SUBROUTINE SECTION_HEADER

         END SUBROUTINE LOGSECTION_PARAMETERS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about memory consumption.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_MEMORY(LOGF, LTYPE, MEM)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, MEMORY
            IMPLICIT NONE

            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            !> Whether this output contains estimated (0) or actual (1) memory
            !> consumption.
            INTEGER,          INTENT(IN) :: LTYPE
            !> Structure containing memory consumption information.
            TYPE(MEMORY),     INTENT(IN) :: MEM

            INTEGER                     :: I
            ! Width of memory string, match with MF string below!
            INTEGER, PARAMETER          :: W = 45
            ! Memory line output format.
            CHARACTER(LEN=*), PARAMETER :: MF = '(A45," : ",F8.1)'
            REAL(q)                     :: RMEM_DYN

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! Maximum dynamical memory.
            RMEM_DYN = MEM%RMEM(5)
            DO I = 6, 11
               RMEM_DYN = MAX(RMEM_DYN, MEM%RMEM(I))
            END DO

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* MEMORY INFORMATION **********************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            IF (LTYPE == 0) THEN
               WRITE(LOGF%IU, '(A)') "Estimated memory consumption for ML &
                                     &force field generation (MB):"
            ELSE IF (LTYPE == 1) THEN
               WRITE(LOGF%IU, '(A)') "Actual memory consumption for ML &
                                     &force field generation (MB):"
            END IF
            WRITE(LOGF%IU, *)
            !WRITE(LOGF%IU, MF) LJ("Static  Coulomb variables"), MEM%RMEM(1)
            !WRITE(LOGF%IU, MF) LJ("Coulomb  descriptors"), MEM%RMEM(6)
            !WRITE(LOGF%IU, MF) LJ("Coulomb  regression"), MEM%RMEM(8)
            !WRITE(LOGF%IU, MF) LJ("Coulomb  prediction"), MEM%RMEM(10)
            WRITE(LOGF%IU, MF) LJ("Persistent allocations for force field"), &
                               MEM%RMEM(2) +              &
                               MEM%RMEM_PS_ALL + MEM%RMEM_SOAP_ALL +         &
                               MEM%RMEM_WMAT_ALL
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- CMAT for basis"), &
                               MEM%RMEM_CMAT + MEM%RMEM_CMAT_ALL
            WRITE(LOGF%IU, MF) LJ("|-- SINV for basis"), &
                               MEM%RMEM_SINV + MEM%RMEM_SINV_ALL
            WRITE(LOGF%IU, MF) LJ("|-- FMAT for basis"), MEM%RMEM_FMAT
            WRITE(LOGF%IU, MF) LJ("|-- DESC for basis"), &
                               MEM%RMEM_PS + MEM%RMEM_PS_ALL
            WRITE(LOGF%IU, MF) LJ("|-- DESC product matrix"), &
                               MEM%RMEM_SOAP + MEM%RMEM_SOAP_ALL
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Persistent allocations for ab initio &
                                  &data"), &
                               MEM%RMEM(3) + MEM%RMEM(4)
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- Ab initio data"), MEM%RMEM(3)
            WRITE(LOGF%IU, MF) LJ("|-- Ab initio data (new)"), MEM%RMEM(4)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Temporary allocations for &
                                  &sparsification"), &
                               MEM%RMEM(5)
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- SVD matrices"), MAX(MEM%RMEM_SVD,MEM%RMEM_RST)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Other temporary allocations"), &
                               MEM%RMEM(7) + MEM%RMEM(9) + MEM%RMEM(11)
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- Descriptors"), MEM%RMEM(7)
            WRITE(LOGF%IU, MF) LJ("|-- Regression"), MEM%RMEM(9)
            WRITE(LOGF%IU, MF) LJ("|-- Prediction"), MEM%RMEM(11)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Total memory consumption"), &
                               MEM%RMEM(2) + MEM%RMEM_CMAT_ALL +             &
                               MEM%RMEM_PS_ALL + MEM%RMEM_SOAP_ALL +         &
                               MEM%RMEM_WMAT_ALL +                           &
                               MEM%RMEM(3) + MEM%RMEM(4) +                   &
                               MEM%RMEM(5) +                                 &
                               MEM%RMEM(7) + MEM%RMEM(9) + MEM%RMEM(11)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

            CONTAINS

            FUNCTION LJ(STRING)
               CHARACTER(LEN=*), INTENT(IN) :: STRING !< Input string.
               CHARACTER(LEN=W)             :: LJ !< Output string.
               LJ = STRLJ(STRING, W)
            END FUNCTION LJ

         END SUBROUTINE LOGSECTION_MEMORY

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Main loop header with log line descriptions.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_LOOP_HEADER(FF, ABN, LOGF, BASICS_ONLY)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            TYPE(FF_PAR),      INTENT(IN) :: FF !< FF setup.
            !> If true, print only "STATUS" and "MSG" in header.
            LOGICAL,           INTENT(IN) :: BASICS_ONLY

            INTEGER :: I
            !> If the number of types in the force field does not match with the ones present in the
            !> structure 
            !> To avoid many "NaN" columns we only write columns for all types in case of ML_MODE =
            !> select.
            INTEGER :: EXPECTED_NTYP
            !> An automatically allocatable string for first column keyword.
            CHARACTER(LEN=:), ALLOCATABLE :: TMPK
            !> An automatically allocatable string for format specifiers.
            CHARACTER(LEN=:), ALLOCATABLE :: TMPF
            !> An automatically allocatable string for separator lines.
            CHARACTER(LEN=:), ALLOCATABLE :: TMPS
            !> A very long string to capture the column list.
            CHARACTER(LEN=1000) :: TMPL
            !> Another very long string to capture the column list.
            CHARACTER(LEN=1000) :: TMPL2

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            IF (FF%ISTART == 3) THEN
               EXPECTED_NTYP = FF%MTYP
            ELSE
               EXPECTED_NTYP = ABN%NTYP
            END IF

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* MAIN LOOP *******************************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            TMPK = "STATUS"
            TMPF = '(A,1X,A9,1X,A9,1X,A2,1X,A6,1X,A6,1X,A9,1X,A9)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "state",                  &
               "is",                     &
               "doabin",                 &
               "iff",                    &
               "nsample",                &
               "ngenff"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line describes the overall status of each step.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ..... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'state ..... One-word description of step action'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "accurate"  (1) : Errors are low, force field is used'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "threshold" (2) : Errors exceeded threshold, structure is sampled from ab initio'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "learning"  (3) : Stored configurations are used for training force field'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "critical"  (4) : Errors are high, ab initio sampling and learning is enforced'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "predict"   (5) : Force field is used in prediction mode only, no error checking'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "predfast"  (6) : Prediction mode only, no error checking, no error estimates'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "interval"  (7) : Errors exceeded threshold, but too close to last sampling (ML_NMDINT)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'is ........ Integer representation of above one-word description (integer in parenthesis)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'doabin .... Perform ab initio calculation (T/F)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'iff ....... Force field available (T/F, False after startup hints to possible convergence problems)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nsample ... Number of steps since last reference structure collection (sample = T)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ngenff .... Number of steps since last force field generation (genff = T)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "MSG"
            TMPF = '(A,1X,A9,1X,A9,1X,A104)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "mtype",                  &
               "message"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows text messages regarding noticeable events in the current step.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ..... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'mtype ..... Message type, one of the following:'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "info"    : General information, not necessarily requiring any user action.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "warning" : Notes about unexpected events, please check log file carefully.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'message ... Actual message.'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 4)
            WRITE(LOGF%IU, *)

            ! Early exit option in case more detailed log lines are not
            ! supported/required in the given mode of operation.
            IF (BASICS_ONLY .AND. FF%IERR.EQ.0) THEN
               WRITE(LOGF%IU, '(140("*"))')
               RETURN
            END IF

            TMPK = "STDAB"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "std_energy",             &
               "std_force",              &
               "std_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the standard deviation of the collected ab initio reference data.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ........ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_energy ... Standard deviation in energy (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_force .... Standard deviation in forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_stress ... Standard deviation in stress (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "ERR"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "rmse_energy",            &
               "rmse_force",             &
               "rmse_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the RMSEs of the predictions with respect to ab initio results for the training data.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_energy ... RMSE of energies (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_force .... RMSE of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_stress ... RMSE of stress (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "NORME"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "nrmse_energy",           &
               "nrmse_force",            &
               "nrmse_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'Contains the normalized root mean square error (NRMSE) (standard'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'deviation in the the training set data is used for normalization).'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'Expressed in percent, i.e., NRMSE = 100% * RMSE / STDDEV.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_energy ... NRMSE of energies (%)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_force .... NRMSE of forces (%)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_stress ... NRMSE of stress (%)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "ERPS"
            TMPF = '(A,1X,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,1(1X,A16)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,1(1X,A16))') TRIM(TMPL2), &
               "el",                                           &
               "rmse_force"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains contributions to the force RMSE for each atomic species.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'NOTE: If this quantity cannot be computed for a particular species at this'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      step the corresponding columns will be filled with "NaN". This occurs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      if the force field contains more species than the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ........ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ........... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_force ... RMSE of forces for specific element (eV Angst^-1)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 2 * FF%MTYP)
            WRITE(LOGF%IU, *)

            TMPK = "NEPS"
            TMPF = '(A,1X,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,1(1X,A16)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,1(1X,A16))') TRIM(TMPL2), &
               "el",                                           &
               "nrmse_force"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'Contains the normalized root mean square error (NRMSE) of forces per atomic'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'species (standard deviation in the the training set data is used for'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'normalization). Expressed in percent, i.e., NRMSE = 100% * RMSE / STDDEV.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'NOTE: If this quantity cannot be computed for a particular species at this'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      step the corresponding columns will be filled with "NaN". This occurs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      if the force field contains more species than the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ............ Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_force ... Force NRMSE for specific element (%)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 2 * FF%MTYP)
            WRITE(LOGF%IU, *)

            TMPK = "LASTE"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "last_energy",            &
               "last_force",             &
               "last_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the RMSEs of the predictions with respect to ab initio results for the last configuration in the training data.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'last_energy ... RMSE of energy for last configuration (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'last_force .... RMSE of forces for last configuration (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'last_stress ... RMSE of stress for last configuration (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "LASTNE"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "nrsme_energy",           &
               "nrsme_force",            &
               "nrsme_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'Contains the normalized root mean square error (NRMSE) of the most recently'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'added structure (standard deviation in the the training set data is used'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'for normalization). Expressed in percent, i.e., NRMSE = 100% * RMSE / STDDEV.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_energy ... Mean relative error of energy'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_force .... Mean relative error of force'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nrmse_stress ... Mean relative error of stress'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "CFE"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "cur_for_mae",            &
               "cur_for_rmse",           &
               "cur_for_maxe"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains errors on forces with respect to ab initio results for the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'cur_for_mae .... Mean absolute error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'cur_for_rmse ... Root mean square error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'cur_for_maxe ... Maximum error of forces (eV Angst^-1)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "THRHIST"
            TMPF = '(A,1X,A9,1X,A4,1X,A16)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "nhis",                   &
               "behis"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the last ML_MHIS entries of the Bayesian errors of the forces for training structures.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'Entries are listed from newest to oldest, i.e. the most recent addition is listed first.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nhis .... Number of the Bayesian error of the force in the history'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'behis ... Bayesian error of the force (eV Angst^-1)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 4)
            WRITE(LOGF%IU, *)

            TMPK = "THRUPD"
            TMPF = '(A,1X,A9,4(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "ctifor_prev",            &
               "ctifor_new",             &
               "std_sig",                &
               "slope_sig"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the new and old threshold for the maximum Bayesian error of forces.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ctifor_prev ... Previous threshold for the maximum Bayesian error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ctifor_new .... New threshold for the maximum Bayesian error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_sig ....... Standard deviation of the collected Bayesion errors of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'slope_sig ..... Slope of the collected Bayesian errors of forces'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 6)
            WRITE(LOGF%IU, *)

            TMPK = "REGR"
            TMPF = '(A,1X,A9,1X,A4,1X,A4,4(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "linreg",                 &
               "ite",                    &
               "sigma_w_sqr",            &
               "sigma_v_sqr",            &
               "regratio",               &
               "error"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains results for each step of linear regression before the last step.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'linreg ........ Linear regression method (ML_IALGO_LINREG)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ite ........... Number of iteration'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_w_sqr ... Square of sigma_w (precision parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_v_sqr ... Square of sigma_v (noise parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'regratio ...... Ratio between (sigma_v_sqr/sigma_w_sqr) and the highest eigenvalue (wmax)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'error.......... 0.5 times squared norm of vector containing combined errors of force field'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "REGRF"
            TMPF = '(A,1X,A9,1X,A4,1x,A4,6(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "linreg",                 &
               "ite",                    &
               "sigma_w_sqr",            &
               "sigma_v_sqr",            &
               "regratio",               &
               "error",                  &
               "wratio",                 &
               "wmax"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains results of the last step of linear regression.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'linreg ........ Linear regression method (ML_IALGO_LINREG)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ite ........... Number of iteration'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_w_sqr ... Square of sigma_w (precision parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_v_sqr ... Square of sigma_v (noise parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'regratio ...... Ratio between (sigma_v_sqr/sigma_w_sqr) and the highest eigenvalue (wmax)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'error.......... 0.5 times squared norm of vector containing combined errors of force field'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'wratio ........ Ratio between the highest (wmax) and lowest (wmin) eigenvalue/singular value (wmax/wmin)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'wmax .......... Highest eigenvalue/singular value (wmax)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 10)
            WRITE(LOGF%IU, *)

            TMPK = "BEE"
            TMPF = '(A,1X,A9,6(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "bee_energy",             &
               "bee_max_force",          &
               "bee_ave_force",          &
               "threshold",              &
               "bee_max_stress",         &
               "bee_ave_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the Bayesian error estimations and the current threshold criterion,'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'immediately before re-training the force field.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_energy ....... Bayesian error estimate of energy per atom (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_force .... Bayesian error estimate of forces (max) (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_force .... Bayesian error estimate of forces (average) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'threshold ........ Current value of threshold criterion (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_stress ... Bayesian error estimate of stresses (max) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_stress ... Bayesian error estimate of stresses (average) (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "BEEF"
            TMPF = '(A,1X,A9,6(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "bee_energy",             &
               "bee_max_force",          &
               "bee_ave_force",          &
               "threshold",              &
               "bee_max_stress",         &
               "bee_ave_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the Bayesian error estimations and the current threshold criterion.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_energy ....... Bayesian error estimate of energy per atom (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_force .... Bayesian error estimate of forces (max) (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_force .... Bayesian error estimate of forces (average) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'threshold ........ Current value of threshold criterion (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_stress ... Bayesian error estimate of stresses (max) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_stress ... Bayesian error estimate of stresses (average) (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "BEPS"
            TMPF = '(A,1X,A9'
            DO I = 1, EXPECTED_NTYP
               TMPF = TMPF // ',1X,A2,2(1X,A16)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, EXPECTED_NTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A16))') TRIM(TMPL2), &
               "el",                                           &
               "bee_max_force",                                &
               "bee_ave_force"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the contributions to Bayesian force error estimation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'for each atomic species, immediately before re-training the force field.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'NOTE: If this quantity cannot be computed for a particular species at this'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      step the corresponding columns will be filled with "NaN". This occurs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      if the force field contains more species than the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ............... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_force .... Bayesian error estimate of forces (max) for specific element (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_force .... Bayesian error estimate of forces (average) for specific element (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 3 * EXPECTED_NTYP)
            WRITE(LOGF%IU, *)

            TMPK = "BEFPS"
            TMPF = '(A,1X,A9'
            DO I = 1, EXPECTED_NTYP
               TMPF = TMPF // ',1X,A2,2(1X,A16)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, EXPECTED_NTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A16))') TRIM(TMPL2), &
               "el",                                           &
               "bee_max_force",                                &
               "bee_ave_force"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the contributions to Bayesian force error estimation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'for each atomic species.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'NOTE: If this quantity cannot be computed for a particular species at this'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      step the corresponding columns will be filled with "NaN". This occurs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      if the force field contains more species than the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ............... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_force .... Bayesian error estimate of forces (max) for specific element (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_force .... Bayesian error estimate of forces (average) for specific element (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 3 * EXPECTED_NTYP)
            WRITE(LOGF%IU, *)

            TMPK = "SF"
            TMPF = '(A,1X,A9,5(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "sfmax",                  &
               "sfmin",                  &
               "sfmean",                 &
               "sfvar",                  &
               "threshold"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the spilling factor,'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'immediately before re-training the force field.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ....... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmax ....... Maximal spilling factor among atoms'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmin ....... Minimal spilling factor among atoms'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmean ...... Mean value of spilling factor'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfvar ....... Variance of spilling factor'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'threshold ... Value of threshold criterion for learning'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 7)
            WRITE(LOGF%IU, *)

            TMPK = "SFF"
            TMPF = '(A,1X,A9,5(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "sfmax",                  &
               "sfmin",                  &
               "sfmean",                 &
               "sfvar",                  &
               "threshold"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the spilling factor.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ....... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmax ....... Maximal spilling factor among atoms'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmin ....... Minimal spilling factor among atoms'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmean ...... Mean value of spilling factor'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfvar ....... Variance of spilling factor'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'threshold ... Value of threshold criterion for learning'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 7)
            WRITE(LOGF%IU, *)

            TMPK = "SPFPS"
            TMPF = '(A,1X,A9'
            DO I = 1, EXPECTED_NTYP
               TMPF = TMPF // ',1X,A2,2(1X,A16)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, EXPECTED_NTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A16))') TRIM(TMPL2), &
               "el",                                           &
               "sfmax",                                        &
               "sfmean"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows contributions to the spilling factor for each'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'atomic species, immediately before re-training the force field.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'NOTE: If this quantity cannot be computed for a particular species at this'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      step the corresponding columns will be filled with "NaN". This occurs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      if the force field contains more species than the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ....... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmax .... Spilling factor (max) over all atoms of specific element'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmean ... Spilling factor (mean) over all atoms of specific element'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 3 * EXPECTED_NTYP)
            WRITE(LOGF%IU, *)

            TMPK = "SPFFPS"
            TMPF = '(A,1X,A9'
            DO I = 1, EXPECTED_NTYP
               TMPF = TMPF // ',1X,A2,2(1X,A16)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, EXPECTED_NTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A16))') TRIM(TMPL2), &
               "el",                                           &
               "sfmax",                                        &
               "sfmean"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows contributions to the spilling factor for each'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'atomic species.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'NOTE: If this quantity cannot be computed for a particular species at this'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      step the corresponding columns will be filled with "NaN". This occurs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '      if the force field contains more species than the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ....... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmax .... Spilling factor (max) over all atoms of specific element'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sfmean ... Spilling factor (mean) over all atoms of specific element'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 3 * EXPECTED_NTYP)
            WRITE(LOGF%IU, *)

            TMPK = "DCOUPLE"
            TMPF = '(A,1X,A9,1X,A16)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "der_H_lambda"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the derivative of the Hamiltonian with respect to coupling constant (dH/dlambda).'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'der_H_lambda ... dH/dlambda'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 3)
            WRITE(LOGF%IU, *)

            TMPK = "FFOUT"
            TMPF = '(A,1X,A9,1X,A23)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "creation_datetime"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the date and time of creation of the ML_FFN file.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'The same string is written to the header of the ML_FFN file, which'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'simplifies the task of identifying later which ML_LOGFILE belongs'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'to which ML_FF(N).'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'creation_datetime ... Date and time of ML_FFN creation'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 3)
            WRITE(LOGF%IU, *)

            TMPK = "LCONF"
            TMPF = '(A,1X,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,2(1X,A9)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A9))') TRIM(TMPL2), &
               "el",                                          &
               "nlrc_old",                                    &
               "nlrc_new"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the number of local configurations'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'which were sampled from ab initio reference calculations.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ...... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ......... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_old ... Previous number of local reference configurations for this element'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_new ... Current number of local reference configurations for this element'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 3 * FF%MTYP)
            WRITE(LOGF%IU, *)

            TMPK = "SPRSC"
            TMPF = '(A,1X,A9,2(1X,A9)'
            WRITE(TMPL, TMPF // ')')     &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "nstr_prev",              &
               "nstr_spar"
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,2(1X,A9)'
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A9))') TRIM(TMPL2), &
               "el",                                        &
               "nlrc_prev",                                 &
               "nlrc_spar"
            END DO
            TMPF = TMPF // ')'
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the results of sparsification regarding the number'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'of reference structures and local reference configurations.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ....... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstr_prev ... Number of reference structures before sparsification'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstr_spar ... Number of reference structures after sparsification'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el .......... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_prev ... Number of local reference configurations before sparsification for this element'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_spar ... Number of local reference configurations after sparsification for this element'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 4 + 3 * FF%MTYP)
            WRITE(LOGF%IU, *)

            TMPK = "NDESC"
            TMPF = '(A,1X,A9,1x,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,1X,A9'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            WRITE(TMPL, '(A,1x,A9)') TRIM(TMPL), "descrad"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1x,A2,1x,A9)') TRIM(TMPL2), &
               "el",                                          &
               "descang"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the number of radial and angular descriptors for each element.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ....... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'descrad ..... Number of radial descriptors'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el .......... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'descang ..... Number of angular descriptors for a given element type'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 3 + 2 * FF%MTYP)
            WRITE(LOGF%IU, *)

            TMPK = "NDESC_SIC"
            TMPF = '(A,1X,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,1X,A9'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1x,A2,1x,A9)') TRIM(TMPL2), &
               "el",                                          &
               "descsic"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the number of self-interaction  '
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'correction (SIC) terms for the angular descriptor.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ....... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el .......... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'descsic ..... Number of SIC descriptors for a given element type'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 2 * FF%MTYP)
            WRITE(LOGF%IU, *)

            WRITE(LOGF%IU, '(140("*"))')

            CONTAINS

            FUNCTION CREATE_SEPARATOR(KEYWORD, LINE) RESULT(SEPARATOR)

               CHARACTER(LEN=*), INTENT(IN) :: KEYWORD
               CHARACTER(LEN=*), INTENT(IN) :: LINE
               CHARACTER(LEN=:), ALLOCATABLE :: SEPARATOR

               SEPARATOR = "# " // KEYWORD // " " // &
                           REPEAT("#", LEN(TRIM(LINE)) - LEN(KEYWORD) - 3)

            END FUNCTION CREATE_SEPARATOR

            SUBROUTINE CREATE_HEADER(KEYWORD, LINE, SEPARATOR, FRMT, NCOL)

               CHARACTER(LEN=*), INTENT(IN) :: KEYWORD
               CHARACTER(LEN=*), INTENT(IN) :: LINE
               CHARACTER(LEN=*), INTENT(IN) :: SEPARATOR
               CHARACTER(LEN=*), INTENT(IN) :: FRMT
               INTEGER,          INTENT(IN) :: NCOL

               ! String representation of column numbers (up to 99, 2 chars).
               CHARACTER(LEN=2) :: COL_NUMBERS(1:99)

               DO I = 1, SIZE(COL_NUMBERS)
                  WRITE(COL_NUMBERS(I), '(I2)') I
               END DO

               WRITE(LOGF%IU, '(A)') SEPARATOR
               WRITE(LOGF%IU, '(A)') TRIM(LINE)
               IF (NCOL <= SIZE(COL_NUMBERS)) THEN
                  WRITE(LOGF%IU, FRMT)            &
                     STRLJ("# " // KEYWORD, C1W), &
                     COL_NUMBERS(2:NCOL)
               END IF
               WRITE(LOGF%IU, '(A)') SEPARATOR

            END SUBROUTINE CREATE_HEADER

         END SUBROUTINE LOGSECTION_LOOP_HEADER

         SUBROUTINE LOGSECTION_AB_INITIO_DATA(ML_INPUT_HANDLE,FF,LOGF, AB, ABN)

            USE ML_FF_STRUCT, ONLY: ML_FF_INPUT, FF_PAR, LOG_PAR, ABINITIO, ABINITIONEW
            IMPLICIT NONE

            TYPE (ML_FF_INPUT)            :: ML_INPUT_HANDLE
            TYPE(FF_PAR),      INTENT(IN) :: FF
            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIO),    INTENT(IN) :: AB
            TYPE(ABINITIONEW), INTENT(IN) :: ABN

            INTEGER :: I
            INTEGER :: J
            INTEGER :: K
            INTEGER :: NSYS_HELP
            INTEGER :: TMP_NITYP
            CHARACTER(LEN=40) :: STRING_HELP

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* AVAILABLE AB INITIO DATA ****************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            WRITE(LOGF%IU, '(A,L)')                               &
               "Bayesian error threshold (CTIFOR) present for all &
               &structures: ",                                    &
               FF%LCTIFOR_PRESENT_IN_MLAB
            WRITE(LOGF%IU, '(A,I9,A,I9,A)')                         &
               "Number of stored (maximum) ab initio structures: ", &
               AB%NCONF,                                            &
               " (",                                                &
               AB%MCONF - ABN%MCONF, ")"
            DO I = 1, FF%NSYS_AB_ONLY
               WRITE(LOGF%IU, '(A,I3,A,I9,A)',ADVANCE='NO') &
                  " * System ",                    &
                  I,                               &
                  " : ",                           &
                  AB%NCONF_SYS(I),                 &
                  ' , Elements: '                 
               DO J=1, AB%NCONF
                  IF (TRIM(ADJUSTL(AB%SZNAM2(J))).EQ.TRIM(ADJUSTL(AB%SZNAM2_ADD(I)))) THEN
                     DO K=1, AB%NTYP(J)
                        WRITE(STRING_HELP,*) AB%NITYP(AB%LADD(K,J),J)
                        WRITE (LOGF%IU, '(A)',ADVANCE='NO') AB%TYPE(AB%LADD(K,J))//" "//TRIM(ADJUSTL(STRING_HELP))//" "
                     ENDDO
                     IF (FF%LUSE_NAMES) THEN 
                        WRITE (LOGF%IU, '(A)',ADVANCE='NO') "Name: "//TRIM(ADJUSTL(AB%SZNAM2_INPUT(J)))
                     ENDIF
                     EXIT
                  ENDIF
               ENDDO
               WRITE(LOGF%IU,*)
            ENDDO
!            WRITE(LOGF%IU, '(A)') "Maximum number of atoms per element:"
!            DO I = 1, AB%MTYP
!               TMP_NITYP = 0
!               DO J = 1, AB%NCONF
!                  TMP_NITYP = MAX(AB%NITYP(I, J), TMP_NITYP)
!               END DO
!               WRITE(LOGF%IU, '(A,A2,A,I9)') &
!                  " * Element ",             &
!                  AB%TYPE(I),                &
!                  " : ",                     &
!                  TMP_NITYP
!            END DO
            IF (FF%ISTART.NE.3) THEN
               WRITE(LOGF%IU, '(A)',ADVANCE='NO') &
                  "New System: "
               DO I=1, ML_INPUT_HANDLE%NTYP
                  WRITE(STRING_HELP,*) ML_INPUT_HANDLE%NITYP(I)
                  WRITE(LOGF%IU,'(A)',ADVANCE='NO') ML_INPUT_HANDLE%TYPE(I)//" "//TRIM(ADJUSTL(STRING_HELP))//" "
                  IF (FF%LUSE_NAMES) THEN 
                     WRITE (LOGF%IU, '(A)',ADVANCE='NO') "Name: "//TRIM(ADJUSTL(ML_INPUT_HANDLE%SZNAM2_INPUT))
                  ENDIF
               ENDDO
            ENDIF
            WRITE(LOGF%IU,*)

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_AB_INITIO_DATA

         SUBROUTINE LOGSECTION_SPARS_DES(LOGF, FF, FFM, NNVAR2_OLD)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, FFM_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            INTEGER,       INTENT(IN) :: NNVAR2_OLD(:)

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* INITIAL SPARSIFICATION OF DESCRIPTORS ***&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(13X,A)') "Number of descriptors"
            WRITE(LOGF%IU, '(A9,1X,A9,1X,A9,1X,A)') &
               "Element",                           &
               "before",                            &
               "after",                             &
               "sparsification"
            DO I = 1, FF%MTYP
               IF (FFM%NB(I) > 0) THEN
                  WRITE(LOGF%IU, '(A9,1X,I9,1X,I9)') &
                     FF%TYPE(I),                     &
                     NNVAR2_OLD(I),                  &
                     FFM%NNVAR2(I)
               END IF
            END DO
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_SPARS_DES

         SUBROUTINE LOGSECTION_TIME(LOGF, CPUTIM, VPUTIM)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN)    :: LOGF !< Log file setup.
            REAL(q),       INTENT(INOUT) :: CPUTIM(:) ! (1:12)
            REAL(q),       INTENT(INOUT) :: VPUTIM(:) ! (1:12)

            INTEGER :: I
            CHARACTER(LEN=:), ALLOCATABLE :: TFMT
            CHARACTER(LEN=50) :: STRTIM(12)

            STRTIM( 1) = "Setup (file I/O, parameters,...)"
            STRTIM( 2) = "Coulomb descriptor and design matrix"
            STRTIM( 3) = "pairwise descriptor and design matrix"
            STRTIM( 4) = "Descriptor and design matrix" ! formerly: "Many-body"
            STRTIM( 5) = "Sparsification of configurations"
            STRTIM( 6) = "Regression for Coulomb interactions"
            STRTIM( 7) = "Regression for pairwise interactions"
            STRTIM( 8) = "Regression" ! formerly: "for many-body interactions"
            STRTIM( 9) = "Prediction for Coulomb interactions"
            STRTIM(10) = "Prediction for pairwise interactions"
            STRTIM(11) = "Prediction" ! formerly: "for many-body interactions
            STRTIM(12) = "TOTAL"

            CPUTIM(12) = 0.0_q
            VPUTIM(12) = 0.0_q
            DO I = 1, 11
               CPUTIM(12) = CPUTIM(12) + CPUTIM(I)
               VPUTIM(12) = VPUTIM(12) + VPUTIM(I)
            ENDDO
            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* TIMING INFORMATION **********************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A50,3X,A18,3X,A18)') &
               STRLJ("Program part", 50),         &
               "system clock (sec)",              &
               "cpu time (sec)"
            WRITE(LOGF%IU, '(A50,("-|-"),A18,("-|-"),A18)') &
               REPEAT("-", 50),                             &
               REPEAT("-", 18),                             &
               REPEAT("-", 18)
            TFMT = '(A50,(" | "),F18.3,(" | "),F18.3)'
            WRITE(LOGF%IU, TFMT) STRTIM( 1), VPUTIM( 1), CPUTIM( 1)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 2), VPUTIM( 2), CPUTIM( 2)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 3), VPUTIM( 3), CPUTIM( 3)
            WRITE(LOGF%IU, TFMT) STRTIM( 4), VPUTIM( 4), CPUTIM( 4)
            WRITE(LOGF%IU, TFMT) STRTIM( 5), VPUTIM( 5), CPUTIM( 5)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 6), VPUTIM( 6), CPUTIM( 6)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 7), VPUTIM( 7), CPUTIM( 7)
            WRITE(LOGF%IU, TFMT) STRTIM( 8), VPUTIM( 8), CPUTIM( 8)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 9), VPUTIM( 9), CPUTIM( 9)
            ! WRITE(LOGF%IU, TFMT) STRTIM(10), VPUTIM(10), CPUTIM(10)
            WRITE(LOGF%IU, TFMT) STRTIM(11), VPUTIM(11), CPUTIM(11)
            WRITE(LOGF%IU, '(A50,("-|-"),A18,("-|-"),A18)') &
               REPEAT("-", 50),                             &
               REPEAT("-", 18),                             &
               REPEAT("-", 18)
            WRITE(LOGF%IU, TFMT) STRTIM(12), VPUTIM(12), CPUTIM(12)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_TIME

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output sparsification result
!****************************************************************************************************

        SUBROUTINE OUT_SPARS_DES_CUR (COMM_WORLD,LOGF,NNVAR2_NEW,NNVAR2_OLD,NRANK,TYPE,WRATIO_NEW,WRATIO_OLD)
          USE ML_FF_STRUCT, ONLY: ML_MPI_PAR, LOG_PAR
          IMPLICIT NONE
          TYPE (ML_MPI_PAR)               :: COMM_WORLD
          TYPE (LOG_PAR)               :: LOGF
          INTEGER         , INTENT(IN) :: NNVAR2_NEW
          INTEGER         , INTENT(IN) :: NNVAR2_OLD
          INTEGER         , INTENT(IN) :: NRANK
          CHARACTER(LEN=2), INTENT(IN) :: TYPE
          REAL(q)         , INTENT(IN) :: WRATIO_NEW
          REAL(q)         , INTENT(IN) :: WRATIO_OLD
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
             WRITE(LOGF%IU,1)
             WRITE(LOGF%IU,2) TYPE,NNVAR2_OLD,NNVAR2_NEW
          ENDIF
1         FORMAT (100("/")/5X,'Spars. des. results')
2         FORMAT (100("/")/5X,'Spars. results on',1X,A2/100("-")/ &
                          10X,'Descriptors from',1X,I5,1X,'to',1X,I5)
        END SUBROUTINE OUT_SPARS_DES_CUR
#endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Left-justify a given string with desired output width.
         !>
         !> If width > len(string) the output string is truncated (cut off at
         !> the right side).
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         PURE FUNCTION STRLJ(STRING, WIDTH)
            CHARACTER(LEN=*), INTENT(IN) :: STRING !< Input string.
            INTEGER, INTENT(IN)          :: WIDTH !< Desired width of output.
            CHARACTER(LEN=WIDTH)         :: STRLJ !< Output string.
            IF (WIDTH < LEN(STRING)) THEN
               STRLJ = STRING(1:WIDTH)
            ELSE
               STRLJ = ADJUSTL(STRING)
            END IF
         END FUNCTION STRLJ

      END MODULE LOGFILE
#endif
