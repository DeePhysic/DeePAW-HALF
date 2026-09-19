#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
      MODULE FORCE_FIELD_ADDITIONAL2
        USE ML_FF_CONSTANT
        USE IOHANDLE
        USE LOGFILE
        USE LAPACK_DATA
        USE MEMORY_USAGE
        USE MPI_DATA
        USE ML_FF_PREC
        USE TIME_CHECK
        USE SOAP_KERNEL
#ifdef use_shmem
        USE ML_FF_MPI_SHMEM
#endif
        IMPLICIT NONE

        CONTAINS
        
        SUBROUTINE COMPUTE_D0CLM(ASA, FF, FFM, MLSHM, PAR_SUP_HANDLE, &
                                      LADD, LNIONS_EST1, LNIONS_EST2, &
                                      LNTYP, MROW_EST1, MROW_EST2, MTYP, &
                                      NIONS, NNEIB_EST1, NNEIB_EST2, &
                                      NTYP, RMEM_HELP1, RMEM_HELP2, XYZ_EST1, XYZ_EST2, CALL_TYPE )
           USE ML_FF_STRUCT, ONLY: ASA_STUFF, PARALLEL_SUPER, ML_SHMEM, FF_PAR, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (ASA_STUFF)      :: ASA
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:)
           INTEGER, INTENT(IN)   :: LNTYP(:)
           INTEGER, INTENT(IN)   :: MROW_EST1
           INTEGER, INTENT(IN)   :: MROW_EST2
           INTEGER, INTENT(IN)   :: MTYP
           INTEGER, INTENT(IN)   :: NIONS
           INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:)
           INTEGER, INTENT(IN)   :: NTYP
           REAL(q), INTENT(IN)   :: XYZ_EST1(:,:,:,:)
           REAL(q), INTENT(IN)   :: XYZ_EST2(:,:,:,:)
           REAL(q), INTENT(INOUT) :: RMEM_HELP1
           REAL(q), INTENT(INOUT) :: RMEM_HELP2
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Local variables
           LOGICAL               :: COUPLE_HELP
           INTEGER               :: IERR
           INTEGER               :: INIONS
           INTEGER               :: INIONS_TMP
           INTEGER               :: INTYP0
           INTEGER               :: LOCVAR_COL_CLM1
           INTEGER               :: LOCVAR_COL_CLM1_COUPLE
           INTEGER               :: LOCVAR_COL_CLM2
           INTEGER               :: LOCVAR_COL_CLM2_COUPLE
           INTEGER               :: LOCVAR_COL_EST1
           INTEGER               :: LOCVAR_COL_EST2
           INTEGER               :: LOCVAR_COL_PS_SIC
           INTEGER               :: LOCVAR_COL_PS_SIC_COUPLE
! MPI parallelization
           PROFILING_START('compute_d0clm')
           IF (CALL_TYPE .EQ. "PRED") THEN
              COUPLE_HELP=FFM%LCOUPLE
           ELSE
              COUPLE_HELP=.FALSE.
           ENDIF
           DO INIONS_TMP=1, NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
              INIONS=INIONS_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
              IF (INIONS.LE.NIONS) THEN
                 INTYP0=LNTYP(INIONS)
                 IF (FFM%W1.GT.0.0_q) THEN
                    LOCVAR_COL_CLM1=INIONS
                    LOCVAR_COL_EST1=CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                    IF (FFM%LCOUPLE) THEN
                       LOCVAR_COL_CLM1_COUPLE=INIONS
                    ELSE
                       LOCVAR_COL_CLM1_COUPLE=1
                    ENDIF
                    LOCVAR_COL_PS_SIC=1
                    LOCVAR_COL_PS_SIC_COUPLE=1
                    CALL D0CLM_SOAP (MLSHM%CLM1(:,LOCVAR_COL_CLM1),MLSHM%CLM_COUPLE_HELP1_MANYBODY(:,LOCVAR_COL_CLM1_COUPLE), &
                                     FFM%DRSPL1,FFM%FNL1,FFM%FWIN1,INTYP0,LADD(1:),COUPLE_HELP,FFM%LFLAG_VAR, &
                                     FFM%LFLAG_VAR_SIC,FFM%LIONS_COUPLE,FFM%LMAX1,LNIONS_EST1(1:,1:,LOCVAR_COL_EST1), &
                                     .FALSE.,FFM%LVAR,FFM%LVAR_SIC,FFM%LWINDOW1,FF%MIONS,MROW_EST1,FFM%MRB1,FFM%MSPL1, &
                                     MTYP,MLSHM%MROW_CLM1_GLOBAL,MLSHM%MROW_CLM_COUPLE_HELP1_MANYBODY,MLSHM%MROW_PS_SIC_GLOBAL, &
                                     MLSHM%MROW_PS_SIC_COUPLE_GLOBAL,NNEIB_EST1(1:,LOCVAR_COL_EST1),FFM%NRB1(0:),FFM%NSPL1, &
                                     NTYP,MLSHM%PS_SIC(:,LOCVAR_COL_PS_SIC),MLSHM%PS_SIC_COUPLE(:,LOCVAR_COL_PS_SIC_COUPLE), &
                                     FFM%RCOUPLE_IONS,RMEM_HELP2,XYZ_EST1(1:,1:,1:,LOCVAR_COL_EST1), &
                                     ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
                    RMEM_HELP1=MAX(RMEM_HELP1,RMEM_HELP2)
                 ENDIF
                 IF (FFM%W2.GT.0.0_q) THEN
                    LOCVAR_COL_CLM2=INIONS
                    LOCVAR_COL_EST2=CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                    IF (FFM%LCOUPLE) THEN
                       LOCVAR_COL_CLM2_COUPLE=INIONS
                    ELSE
                       LOCVAR_COL_CLM2_COUPLE=1
                    ENDIF
                    IF (FFM%LSIC) THEN
                       LOCVAR_COL_PS_SIC=INIONS
                       IF (FFM%LCOUPLE) THEN
                          LOCVAR_COL_PS_SIC_COUPLE=INIONS
                       ELSE
                          LOCVAR_COL_PS_SIC_COUPLE=1
                       ENDIF
                    ELSE
                       LOCVAR_COL_PS_SIC=1
                       LOCVAR_COL_PS_SIC_COUPLE=1
                    ENDIF
                    IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7) THEN
                       CALL D0CLM_SOAP_NOSIC_DUAL (MLSHM%CLM2_RED(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                            FFM%DRSPL2,FFM%FNL2,FFM%RV%FNL2_RED,FFM%FWIN2, &
                            INTYP0,LADD(1:), &
                            FFM%LMAX2,FFM%RV%LMAX2_RED,LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                            FFM%LWINDOW2,FF%MIONS,MROW_EST2,FFM%MRB2,FFM%RV%MRB2_RED, &
                            FFM%MSPL2,MTYP, &
                            MLSHM%MROW_CLM2_GLOBAL, &
                            NNEIB_EST2(1:,LOCVAR_COL_EST2),FFM%NRB2(0:),FFM%RV%NRB2_RED(0:),FFM%NSPL2,NTYP, &
                            FFM%RCOUPLE_IONS,RMEM_HELP2,XYZ_EST2(1:,1:,1:,LOCVAR_COL_EST2), &
                            ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
                    ELSE
                       CALL D0CLM_SOAP (MLSHM%CLM2(:,LOCVAR_COL_CLM2),MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                                        FFM%DRSPL2,FFM%FNL2,FFM%FWIN2,INTYP0,LADD(1:),COUPLE_HELP,FFM%LFLAG_VAR, &
                                        FFM%LFLAG_VAR_SIC,FFM%LIONS_COUPLE,FFM%LMAX2,LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                                        FFM%LSIC,FFM%LVAR,FFM%LVAR_SIC,FFM%LWINDOW2,FF%MIONS,MROW_EST2,FFM%MRB2,FFM%MSPL2, &
                                        MTYP,MLSHM%MROW_CLM2_GLOBAL,MLSHM%MROW_CLM_COUPLE_HELP2_MANYBODY,MLSHM%MROW_PS_SIC_GLOBAL, &
                                        MLSHM%MROW_PS_SIC_COUPLE_GLOBAL,NNEIB_EST2(1:,LOCVAR_COL_EST2),FFM%NRB2(0:),FFM%NSPL2, &
                                        NTYP,MLSHM%PS_SIC(:,LOCVAR_COL_PS_SIC),MLSHM%PS_SIC_COUPLE(:,LOCVAR_COL_PS_SIC_COUPLE), &
                                        FFM%RCOUPLE_IONS,RMEM_HELP2,XYZ_EST2(1:,1:,1:,LOCVAR_COL_EST2), &
                                        ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
                    ENDIF
                    RMEM_HELP1=MAX(RMEM_HELP1,RMEM_HELP2)
                 ENDIF
              ENDIF
           ENDDO
           
#ifdef use_shmem
           CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
           !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#endif
! Reduction of expansion coefficients so that each process possesses the whole array CLM and PS_SIC.
#ifdef use_shmem
           IF (FFM%W1.GT.0.0_q) THEN
              IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM1,MLSHM%MROW_CLM1_GLOBAL*NIONS,MPI_REAL8,MPI_SUM, &
                      PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
              ENDIF
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF(FFM%LCOUPLE) THEN
                    IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                       CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM_COUPLE_HELP1_MANYBODY, &
                                          MLSHM%MROW_CLM_COUPLE_HELP1_MANYBODY*NIONS,MPI_REAL8,MPI_SUM, &
                            PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
           IF (FFM%W2.GT.0.0_q) THEN
              IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM2,MLSHM%MROW_CLM2_GLOBAL*NIONS,MPI_REAL8,MPI_SUM, &
                      PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                 IF (FFM%LSIC) THEN
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%PS_SIC,MLSHM%MROW_PS_SIC_GLOBAL*NIONS,MPI_REAL8,MPI_SUM, &
                         PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                 ENDIF
                 IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM2_RED,MLSHM%MROW_CLM2_GLOBAL_RED*NIONS,MPI_REAL8,MPI_SUM, &
                         PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                 ENDIF
              ENDIF
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF (FFM%LCOUPLE) THEN
                    IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                       CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM_COUPLE_HELP2_MANYBODY, &
                                          MLSHM%MROW_CLM_COUPLE_HELP2_MANYBODY*NIONS,MPI_REAL8,MPI_SUM, &
                                          PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                       IF (FFM%LSIC) THEN
                          CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%PS_SIC_COUPLE, &
                                             MLSHM%MROW_PS_SIC_COUPLE_GLOBAL*NIONS,MPI_REAL8,MPI_SUM, &
                                             PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
           CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
           !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#else
           IF (FFM%W1.GT.0.0_q) THEN
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM1,MLSHM%MROW_CLM1_GLOBAL*NIONS, &
                                 MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF (FFM%LCOUPLE) THEN
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM_COUPLE_HELP1_MANYBODY, &
                                       MLSHM%MROW_CLM_COUPLE_HELP1_MANYBODY*NIONS, &
                                       MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
                 ENDIF
              ENDIF
           ENDIF
           IF (FFM%W2.GT.0.0_q) THEN
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM2,MLSHM%MROW_CLM2_GLOBAL*NIONS, &
                                 MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              IF(FFM%LSIC) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%PS_SIC,MLSHM%MROW_PS_SIC_GLOBAL*NIONS, &
                                    MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              ENDIF
              IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM2_RED,MLSHM%MROW_CLM2_GLOBAL_RED*NIONS, &
                                    MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              ENDIF
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF (FFM%LCOUPLE) THEN
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%CLM_COUPLE_HELP2_MANYBODY, &
                                       MLSHM%MROW_CLM_COUPLE_HELP2_MANYBODY*NIONS, &
                                       MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
                    IF (FFM%LSIC) THEN
                       CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%PS_SIC_COUPLE, &
                                          MLSHM%MROW_PS_SIC_COUPLE_GLOBAL*NIONS,MPI_REAL8, &
                                          MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
                    ENDIF
                 END IF
              ENDIF
           ENDIF
#endif
           PROFILING_STOP('compute_d0clm')
        END SUBROUTINE COMPUTE_D0CLM


! Subroutine to calculate radial descriptors by using fast method.
        SUBROUTINE COMPUTE_D0C00(FFM, PAR_SUP_HANDLE, MLSHM, C00, C00_COUPLE, &
                                      INTYP0, LADD, LNIONS, MROW_C00, MROW_C00_COUPLE, MTYP, &
                                      NITYP, NTYP, CALL_TYPE) 
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, ML_SHMEM, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: C00(:,:,:)
           REAL(q), INTENT(INOUT) :: C00_COUPLE(:,:)
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS(:,:)
           INTEGER, INTENT(IN)   :: MROW_C00
           INTEGER, INTENT(IN)   :: MROW_C00_COUPLE
           INTEGER, INTENT(IN)   :: MTYP
           INTEGER, INTENT(IN)   :: NITYP(:)
           INTEGER, INTENT(IN)   :: NTYP
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Local variables
           LOGICAL               :: COUPLE_HELP
           INTEGER               :: INIONS
           INTEGER               :: INITYP
           INTEGER               :: INITYP_TMP
           INTEGER               :: LOCVAR_COL_C00
           INTEGER               :: LOCVAR_COL_C00_COUPLE
           INTEGER               :: LOCVAR_COL_CLM1
           INTEGER               :: LOCVAR_COL_CLM1_COUPLE
           PROFILING_START('compute_d0c00')
           IF (FFM%W1.GT.0.0_q) THEN
              IF (CALL_TYPE .EQ. "PRED") THEN
                 COUPLE_HELP=FFM%LCOUPLE
              ELSE
                 COUPLE_HELP=.FALSE.
              ENDIF
! MPI-parallel over INITYP loop.
              DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                 INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                 IF (INITYP.LE.NITYP(INTYP0)) THEN
                    INIONS=LNIONS(INITYP,INTYP0)
                    LOCVAR_COL_C00=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                    LOCVAR_COL_CLM1=INIONS
                    IF(FFM%LCOUPLE .AND. CALL_TYPE.EQ."PRED") THEN
                       LOCVAR_COL_C00_COUPLE=LOCVAR_COL_C00
                       LOCVAR_COL_CLM1_COUPLE=INIONS
                    ELSE
                       LOCVAR_COL_C00_COUPLE=1
                       LOCVAR_COL_CLM1_COUPLE=1
                    ENDIF
                    CALL D0C00_SOAP (C00(1:,LOCVAR_COL_C00,INTYP0),C00_COUPLE(1:,LOCVAR_COL_C00_COUPLE), &
                                     MLSHM%CLM1(:,LOCVAR_COL_CLM1), &
                                     MLSHM%CLM_COUPLE_HELP1_MANYBODY(:,LOCVAR_COL_CLM1_COUPLE), &
                                     LADD(1:),COUPLE_HELP,FFM%LMAX1,FFM%LNORM1, &
                                     MROW_C00,MROW_C00_COUPLE,FFM%MRB1,MTYP,MLSHM%MROW_CLM1_GLOBAL, &
                                     MLSHM%MROW_CLM_COUPLE_HELP1_MANYBODY,FFM%NRB1(0:),NTYP)
                 ENDIF
              ENDDO
           ENDIF
           PROFILING_STOP('compute_d0c00')
        END SUBROUTINE COMPUTE_D0C00

! Subroutine to calculate angular descriptors.

        SUBROUTINE COMPUTE_D0PS(FFM, MLSHM, PAR_SUP_HANDLE, INTYP0, &
                                     LADD, LNIONS, MROW_PS, MROW_PS_COUPLE, &
                                     MTYP, NITYP, NTYP, PS, PS_COUPLE, CALL_TYPE)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, ML_SHMEM, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: PS(:,:,:)
           REAL(q), INTENT(INOUT) :: PS_COUPLE(:,:)
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS(:,:)
           INTEGER, INTENT(IN)   :: MROW_PS
           INTEGER, INTENT(IN)   :: MROW_PS_COUPLE
           INTEGER, INTENT(IN)   :: MTYP
           INTEGER, INTENT(IN)   :: NITYP(:)
           INTEGER, INTENT(IN)   :: NTYP
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Local variables
           LOGICAL               :: COUPLE_HELP
           INTEGER               :: INIONS
           INTEGER               :: INTYP
           INTEGER               :: INITYP
           INTEGER               :: INITYP_TMP
           INTEGER               :: LOCVAR_COL_CLM2
           INTEGER               :: LOCVAR_COL_CLM2_COUPLE
           INTEGER               :: LOCVAR_COL_PS
           INTEGER               :: LOCVAR_COL_PS_COUPLE
           INTEGER               :: LOCVAR_COL_PS_SIC
           INTEGER               :: LOCVAR_COL_PS_SIC_COUPLE
           PROFILING_START('compute_d0ps')
           IF (FFM%W2.GT.0.0_q) THEN
              IF (CALL_TYPE .EQ. "PRED") THEN
                 COUPLE_HELP=FFM%LCOUPLE
              ELSE
                 COUPLE_HELP=.FALSE.
              ENDIF
! MPI-parallel.
              DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                 INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                 IF (INITYP.LE.NITYP(INTYP0)) THEN
                    INIONS=LNIONS(INITYP,INTYP0)
                    LOCVAR_COL_PS=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                    LOCVAR_COL_CLM2=INIONS
                    IF(FFM%LCOUPLE .AND. CALL_TYPE.EQ."PRED") THEN
                       LOCVAR_COL_CLM2_COUPLE=INIONS
                       LOCVAR_COL_PS_COUPLE=LOCVAR_COL_PS
                    ELSE
                       LOCVAR_COL_PS_COUPLE=1
                       LOCVAR_COL_CLM2_COUPLE=1
                    ENDIF
                    IF(FFM%LSIC) THEN
                       LOCVAR_COL_PS_SIC=INIONS
                       IF(FFM%LCOUPLE .AND. CALL_TYPE.EQ."PRED") THEN
                          LOCVAR_COL_PS_SIC_COUPLE=INIONS
                       ELSE
                          LOCVAR_COL_PS_SIC_COUPLE=1
                       ENDIF
                    ELSE
                       LOCVAR_COL_PS_SIC=1
                       LOCVAR_COL_PS_SIC_COUPLE=1
                    ENDIF
                    IF (FFM%DESC_TYPE.EQ.0) THEN
                       IF (FFM%LSIC) THEN
                          CALL D0PS_SOAP (FFM%AFILT2, MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                               MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                               FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                               FFM%LFLAG_VAR(:,:,0:,:,:,:), FFM%LFLAG_VAR_SIC(:,:,0:,:,:), &
                               FFM%L_LNRB2, FFM%LMAX2, FFM%LSIC, &
                               FFM%LVAR(:,:,0:,:,:,:), FFM%LVAR_SIC(:,:,0:,:,:), &
                               FFM%MLNRB2, &
                               FFM%MRB2, MTYP, &
                               FFM%NLNRB2, FFM%NRB2(0:), FFM%NRB_LNRB2, NTYP, &
                               PS(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                               MLSHM%PS_SIC(:,LOCVAR_COL_PS_SIC), &
                               MLSHM%PS_SIC_COUPLE(:,LOCVAR_COL_PS_SIC_COUPLE), FFM%WVAR(:,:,0:,:,:,:))
                       ELSE
                          CALL D0PS_SOAP_NOSIC (FFM%AFILT2, MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                               MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                               FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                               FFM%LFLAG_VAR(:,:,0:,:,:,:), &
                               FFM%L_LNRB2, FFM%LMAX2, &
                               FFM%LVAR(:,:,0:,:,:,:), &
                               FFM%MLNRB2, &
                               FFM%MRB2, MTYP, &
                               FFM%NLNRB2, FFM%NRB2(0:), FFM%NRB_LNRB2, NTYP, &
                               PS(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                               FFM%WVAR(:,:,0:,:,:,:), &
                               FFM%LMAX2)
                       ENDIF
                    ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
                       CALL D0PS_SOAP_NOSIC_LIN_ELEMENT (FFM%AFILT2, MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                            FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                            FFM%LFLAG_VAR(:,:,0:,:,:,:), &
                            FFM%L_LNRB2, FFM%LMAX2, &
                            FFM%LVAR(:,:,0:,:,:,:), &
                            FFM%MLNRB2, &
                            FFM%MRB2, MTYP, &
                            FFM%NLNRB2, FFM%NRB2(0:), FFM%NRB_LNRB2, NTYP, &
                            PS(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                            FFM%WVAR(:,:,0:,:,:,:))
                    ELSE IF (FFM%DESC_TYPE.EQ.2) THEN
                       ! High resolution type-independent descriptor
                       CALL D0PS_SOAP_NOSIC_NO_ELEMENT (FFM%AFILT2, MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                            FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                            FFM%LFLAG_VAR(:,:,0:,:,:,:), &
                            FFM%L_LNRB2, FFM%LMAX2, &
                            FFM%LVAR(:,:,0:,:,:,:), &
                            FFM%MLNRB2, &
                            FFM%MRB2, MTYP, &
                            FFM%NLNRB2, FFM%NRB2(0:), FFM%NRB_LNRB2, NTYP, &
                            PS(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                            FFM%WVAR(:,:,0:,:,:,:))
                       !  Low resolution type-dependent descriptor
                       CALL D0PS_SOAP_NOSIC (FFM%AFILT2, MLSHM%CLM2_RED(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                            FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                            FFM%RV%LFLAG_VAR_RED(:,:,0:,:,:,:), &
                            FFM%RV%L_LNRB2_RED, FFM%RV%LMAX2_RED, &
                            FFM%RV%LVAR_RED(:,:,0:,:,:,:), &
                            FFM%RV%MLNRB2_RED, &
                            FFM%RV%MRB2_RED, MTYP, &
                            FFM%RV%NLNRB2_RED, FFM%RV%NRB2_RED(0:), FFM%RV%NRB_LNRB2_RED, NTYP, &
                            FFM%RV%PS_UGANG_RED(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                            FFM%RV%WVAR_RED(:,:,0:,:,:,:), &
                            FFM%RV%LMAX2_RED)
                    ELSE IF (FFM%DESC_TYPE.EQ.3) THEN
                       CALL D0PS_SOAP_NOSIC_LIN_MRB (FFM%AFILT2, MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                            FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                            FFM%LFLAG_VAR(:,:,0:,:,:,:), &
                            FFM%L_LNRB2, FFM%LMAX2, &
                            FFM%LVAR(:,:,0:,:,:,:), &
                            FFM%MLNRB2, &
                            FFM%MRB2, MTYP, &
                            FFM%NLNRB2, FFM%NRB2(0:), FFM%NRB_LNRB2, NTYP, &
                            PS(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                            FFM%WVAR(:,:,0:,:,:,:), &
                            FFM%LMAX2)
                    ELSE IF (FFM%DESC_TYPE.EQ.4) THEN
                       ! First descriptor
                       CALL D0PS_SOAP_NOSIC (FFM%AFILT2, MLSHM%CLM2(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                            FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                            FFM%LFLAG_VAR(:,:,0:,:,:,:), &
                            FFM%L_LNRB2, FFM%LMAX2, &
                            FFM%LVAR(:,:,0:,:,:,:), &
                            FFM%MLNRB2, &
                            FFM%MRB2, MTYP, &
                            FFM%NLNRB2, FFM%NRB2(0:), FFM%NRB_LNRB2, NTYP, &
                            PS(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                            FFM%WVAR(:,:,0:,:,:,:), &
                            FFM%LMAX2)
                       ! Second descriptor
                       CALL D0PS_SOAP_NOSIC (FFM%AFILT2, MLSHM%CLM2_RED(:,LOCVAR_COL_CLM2), &
                            MLSHM%CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE), &
                            FFM%IAFILT2, INTYP0, FFM%LAFILT2, LADD(1:), COUPLE_HELP, &
                            FFM%RV%LFLAG_VAR_RED(:,:,0:,:,:,:), &
                            FFM%RV%L_LNRB2_RED, FFM%RV%LMAX2_RED, &
                            FFM%RV%LVAR_RED(:,:,0:,:,:,:), &
                            FFM%RV%MLNRB2_RED, &
                            FFM%RV%MRB2_RED, MTYP, &
                            FFM%RV%NLNRB2_RED, FFM%RV%NRB2_RED(0:), FFM%RV%NRB_LNRB2_RED, NTYP, &
                            FFM%RV%PS_UGANG_RED(1:,LOCVAR_COL_PS,INTYP0), PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                            FFM%RV%WVAR_RED(:,:,0:,:,:,:), &
                            FFM%RV%LMAX2_RED)
                    ENDIF
                 ENDIF
              ENDDO
           ENDIF
           PROFILING_STOP('compute_d0ps')
        END SUBROUTINE COMPUTE_D0PS

! Normalization.

        SUBROUTINE COMPUTE_NORMD0(FFM, PAR_SUP_HANDLE, C00, C00_COUPLE, INTYP0, &
                                                 LNIONS, MROW_C00, MROW_C00_COUPLE, MROW_PS, MROW_PS_COUPLE, &
                                                 NITYP, PS, PS_COUPLE, RNORM_C00, RNORM_PS, RNORM_ALL, &
                                                 CALL_TYPE)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (FFM_PAR)        :: FFM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: C00(:,:,:)
           REAL(q), INTENT(INOUT) :: C00_COUPLE(:,:)
           REAL(q), INTENT(INOUT) :: PS(:,:,:)
           REAL(q), INTENT(INOUT) :: PS_COUPLE(:,:)
           REAL(q), INTENT(INOUT) :: RNORM_C00(:,:)
           REAL(q), INTENT(INOUT) :: RNORM_PS(:,:)
           REAL(q), INTENT(INOUT) :: RNORM_ALL(:,:)
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: LNIONS(:,:)
           INTEGER, INTENT(IN)   :: MROW_C00
           INTEGER, INTENT(IN)   :: MROW_C00_COUPLE
           INTEGER, INTENT(IN)   :: MROW_PS
           INTEGER, INTENT(IN)   :: MROW_PS_COUPLE
           INTEGER, INTENT(IN)   :: NITYP(:)
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Local variables
           LOGICAL               :: COUPLE_HELP
           INTEGER               :: IERR
           INTEGER               :: INIONS 
           INTEGER               :: INITYP
           INTEGER               :: INITYP_TMP
           INTEGER               :: LOCVAR_COL_C00
           INTEGER               :: LOCVAR_COL_C00_COUPLE
           INTEGER               :: LOCVAR_COL_PS
           INTEGER               :: LOCVAR_COL_PS_RED
           INTEGER               :: LOCVAR_COL_PS_COUPLE
           INTEGER               :: LOCVAR_COL_RNORM_ALL
           INTEGER               :: LOCVAR_COL_RNORM_C00
           INTEGER               :: LOCVAR_COL_RNORM_PS
           INTEGER               :: LOCVAR_COL_RNORM_PS_RED
           PROFILING_START('compute_normd0')
! Initializations of norm.
           IF (.NOT.FFM%LSUPERVEC) THEN
              IF (FFM%W1.GT.0.0_q) THEN
                 RNORM_C00(1:NITYP(INTYP0),INTYP0)=0.0_q
              ENDIF
              IF (FFM%W2.GT.0.0_q) THEN
                 RNORM_PS(1:NITYP(INTYP0),INTYP0)=0.0_q
                 IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                    FFM%RV%RNORM_PS_RED(1:NITYP(INTYP0),INTYP0)=0.0_q
                 ENDIF
              ENDIF
           ELSE
              RNORM_ALL(1:NITYP(INTYP0),INTYP0)=0.0_q
           ENDIF
           IF (CALL_TYPE .EQ. "PRED") THEN
              COUPLE_HELP=FFM%LCOUPLE
           ELSE
              COUPLE_HELP=.FALSE.
           ENDIF
! MPI-parallel.
           DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
              INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
              IF (INITYP.LE.NITYP(INTYP0)) THEN
                 INIONS=LNIONS(INITYP,INTYP0)
                 IF(.NOT.FFM%LSUPERVEC) THEN
                    LOCVAR_COL_C00_COUPLE=1
                    IF(FFM%W1.GT.0.0_q) THEN
                       LOCVAR_COL_C00=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                       LOCVAR_COL_RNORM_C00=INITYP
                       IF (CALL_TYPE .EQ. "PRED") THEN 
                          IF (FFM%LCOUPLE) THEN
                             LOCVAR_COL_C00_COUPLE=LOCVAR_COL_C00
                          ENDIF
                       ENDIF
                    ELSE
                       LOCVAR_COL_C00=1
                       LOCVAR_COL_RNORM_C00=1
                    ENDIF
                    LOCVAR_COL_PS_COUPLE=1
                    IF(FFM%W2.GT.0.0_q) THEN
                       LOCVAR_COL_PS=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                       LOCVAR_COL_RNORM_PS=INITYP
                       IF (CALL_TYPE .EQ. "PRED") THEN 
                          IF (FFM%LCOUPLE) THEN
                             LOCVAR_COL_PS_COUPLE=LOCVAR_COL_PS
                          ENDIF
                       ENDIF
                       IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                          LOCVAR_COL_RNORM_PS_RED=LOCVAR_COL_RNORM_PS
                          LOCVAR_COL_PS_RED=LOCVAR_COL_PS
                       ELSE
                          LOCVAR_COL_RNORM_PS_RED=1
                          LOCVAR_COL_PS_RED=1
                       ENDIF
                    ELSE
                       LOCVAR_COL_PS=1
                       LOCVAR_COL_RNORM_PS=1
                       LOCVAR_COL_RNORM_PS_RED=1
                       LOCVAR_COL_PS_RED=1
                    ENDIF
                    LOCVAR_COL_RNORM_ALL=1
                 ELSE
                    LOCVAR_COL_C00_COUPLE=1
                    IF(FFM%W1.GT.0.0_q) THEN
                       LOCVAR_COL_C00=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                       IF (CALL_TYPE .EQ. "PRED") THEN 
                          IF (FFM%LCOUPLE) THEN
                             LOCVAR_COL_C00_COUPLE=LOCVAR_COL_C00
                          ENDIF
                       ENDIF
                    ELSE
                       LOCVAR_COL_C00=1
                    ENDIF
                    LOCVAR_COL_PS_COUPLE=1
                    IF(FFM%W2.GT.0.0_q) THEN
                       LOCVAR_COL_PS=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                       IF (CALL_TYPE .EQ. "PRED") THEN 
                          IF (FFM%LCOUPLE) THEN
                             LOCVAR_COL_PS_COUPLE=LOCVAR_COL_PS
                          ENDIF
                       ENDIF
                       IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                          LOCVAR_COL_PS_RED=LOCVAR_COL_PS
                       ELSE
                          LOCVAR_COL_PS_RED=1
                       ENDIF
                    ELSE
                       LOCVAR_COL_PS=1
                       LOCVAR_COL_PS_RED=1
                    ENDIF
                    LOCVAR_COL_RNORM_C00=1
                    LOCVAR_COL_RNORM_PS=1
                    LOCVAR_COL_RNORM_PS_RED=1
                    LOCVAR_COL_RNORM_ALL=INITYP
                 ENDIF
                 ! For DESC_TYPE=2 and DESC_TYPE=4 the reduced variables are called from FFM inside the subroutine
                 CALL NORMALIZATION_D0 (FFM,C00(1:,LOCVAR_COL_C00,INTYP0),C00_COUPLE(1:,LOCVAR_COL_C00_COUPLE), &
                                        FFM%LSUPERVEC,COUPLE_HELP,FFM%LNORM1,FFM%LNORM2, &
                                        FFM%NNVAR1, &
                                        FFM%NNVAR2(INTYP0),FFM%RV%NNVAR2_RED(INTYP0), &
                                        PS(1:,LOCVAR_COL_PS,INTYP0), &
                                        FFM%RV%PS_UGANG_RED(1:,LOCVAR_COL_PS_RED,INTYP0), &
                                        PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE), &
                                        RNORM_ALL(LOCVAR_COL_RNORM_ALL,INTYP0), &
                                        RNORM_C00(LOCVAR_COL_RNORM_C00,INTYP0), &
                                        RNORM_PS(LOCVAR_COL_RNORM_PS,INTYP0), &
                                        FFM%RV%RNORM_PS_RED(LOCVAR_COL_RNORM_PS_RED,INTYP0), &
                                        FFM%W1,FFM%W2)
              ENDIF
           ENDDO
           IF(.NOT.FFM%LSUPERVEC) THEN
              IF(FFM%W1.GT.0.0_q) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,RNORM_C00(1,INTYP0),NITYP(INTYP0), &
                                    MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              ENDIF
              IF(FFM%W2.GT.0.0_q) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,RNORM_PS(1,INTYP0),NITYP(INTYP0), &
                                    MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
                 IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,FFM%RV%RNORM_PS_RED(1,INTYP0),NITYP(INTYP0), &
                                       MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
                 ENDIF
              ENDIF
           ELSE
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,RNORM_ALL(1,INTYP0),NITYP(INTYP0), &
                                 MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
           ENDIF
           PROFILING_STOP('compute_normd0')
        END SUBROUTINE COMPUTE_NORMD0


! Subroutine to calculate SOAP similarity kernel using fast method.

        SUBROUTINE COMPUTE_SOAP(FF, FFM, MLSHM, PAR_SUP_HANDLE, C00, &
                                            INTYP, INTYP0, LHEAD, LNIONS, MCOL_SOAP, MROW_C00, &
                                            MROW_PS, MROW_SOAP, NB_MAX, NCOL_NITYP, NIONS, NITYP, &
                                            NTYP, PS, SOAP1, SOAP2, CALL_TYPE, &
                                            MROW_C00_COUPLE, MROW_PS_COUPLE, &
                                            C00_COUPLE, PS_COUPLE, &
                                            MROW_SOAP_COUPLE, MCOL_SOAP_COUPLE, &
                                            SOAP1_COUPLE, SOAP2_COUPLE)
           USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, ML_SHMEM, PARALLEL_SUPER, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: C00(:,:,:)
           REAL(q), INTENT(INOUT) :: PS(:,:,:)
           REAL(q), INTENT(INOUT) :: SOAP1(:,:)
           REAL(q), INTENT(INOUT) :: SOAP2(:,:)
           INTEGER, INTENT(IN)   :: INTYP
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: LHEAD(:)
           INTEGER, INTENT(IN)   :: LNIONS(:,:)
           INTEGER, INTENT(IN)   :: MCOL_SOAP
           INTEGER, INTENT(IN)   :: MROW_C00
           INTEGER, INTENT(IN)   :: MROW_PS
           INTEGER, INTENT(IN)   :: MROW_SOAP
           INTEGER, INTENT(IN)   :: NB_MAX
           INTEGER, INTENT(IN)   :: NCOL_NITYP
           INTEGER, INTENT(IN)   :: NIONS
           INTEGER, INTENT(IN)   :: NITYP(:)
           INTEGER, INTENT(IN)   :: NTYP
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
           INTEGER, INTENT(IN)   :: MROW_C00_COUPLE
           INTEGER, INTENT(IN)   :: MROW_PS_COUPLE
           REAL(q), INTENT(INOUT) :: C00_COUPLE(:,:)
           REAL(q), INTENT(INOUT) :: PS_COUPLE(:,:)
           INTEGER, INTENT(IN)   :: MROW_SOAP_COUPLE
           INTEGER, INTENT(IN)   :: MCOL_SOAP_COUPLE
           REAL(q), INTENT(INOUT) :: SOAP1_COUPLE(:,:)
           REAL(q), INTENT(INOUT) :: SOAP2_COUPLE(:,:)
! Local variables
           INTEGER               :: IB
           INTEGER               :: IERR
           INTEGER               :: IHEAD
           INTEGER               :: INIONS
           INTEGER               :: INITYP
           INTEGER               :: INITYP_TMP
           INTEGER               :: IPOLY
           INTEGER               :: LOCVAR_COL
           INTEGER               :: NATOM_NUM_HELP
           INTEGER               :: NB_USE
           REAL(q)               :: PREFAC_HELP
           REAL(q)               :: PREFAC_HELP_RED

           PROFILING_START('compute_soap')
! Set funny integer helping variable for checking if radial descriptor was used (see later)
           IPOLY=0
! Set descriptor weightfactor
           IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
              PREFAC_HELP=FFM%RV%FACTOR_DESC * FFM%RV%RATIO_DESC * FFM%W2
              PREFAC_HELP_RED=FFM%RV%FACTOR_DESC * FFM%RV%RATIO_DESC_LOW_RES * FFM%W2
           ELSE
              PREFAC_HELP=FFM%W2
           ENDIF
           ! Set header in the array.
           IF (CALL_TYPE .EQ. "PUTA") THEN
#ifdef scaLAPACK
              IHEAD=FFM%NB(INTYP0)+LHEAD(INTYP0)
#else
              IHEAD=FFM%NB(INTYP0)+1
#endif
           ELSE
#ifdef scaLAPACK
              IHEAD=LHEAD(INTYP0)
#else
              IHEAD=1
#endif
           ENDIF
           ! Set atom divisor
           IF (CALL_TYPE .EQ. "PRED") THEN
              NATOM_NUM_HELP = 1
           ELSE IF (CALL_TYPE .EQ. "PUTA") THEN
              NATOM_NUM_HELP = NIONS
           ELSE
              IF (FF%LTOTEN_SYSTEM) THEN
                 NATOM_NUM_HELP = 1
              ELSE
                 NATOM_NUM_HELP = NIONS
              ENDIF
           ENDIF
           ! Set the right NB here
           IF (CALL_TYPE .EQ. "PUTA") THEN
              NB_USE=FFM%NB_NEW(INTYP0)
           ELSE
              NB_USE=FFM%NB(INTYP0)
           ENDIF
! Calculate contribution from 1-point correlation function
           IF (FFM%W1.GT.0.0_q) THEN
              IF (NB_USE.GT.0.AND.NCOL_NITYP.GT.0) THEN
#ifdef scaLAPACK
                 CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR1,1.0_q, &
                             FFM%C00_ALL(1,IHEAD),FFM%MMVAR1,C00(1,1,INTYP0), &
                             MROW_C00,0.0_q,SOAP1(1,1),MROW_SOAP)
                 IF (CALL_TYPE .EQ. "PRED") THEN
                    IF (FFM%LCOUPLE) THEN
                       CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR1,1.0_q, &
                            FFM%C00_ALL(1,IHEAD),FFM%MMVAR1, &
                            C00_COUPLE(1,1),MROW_C00_COUPLE, &
                            0.0_q,SOAP1_COUPLE(1,1),MROW_SOAP)
                    ENDIF
                 ENDIF
#else
                 CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR1,1.0_q, &
                        FFM%C00(1,IHEAD,INTYP0),FFM%MMVAR1, &
                        C00(1,1,INTYP0),MROW_C00, &
                        0.0_q,SOAP1(1,1),MROW_SOAP)
                 IF (CALL_TYPE .EQ. "PRED") THEN
                    IF (FFM%LCOUPLE) THEN
                       CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR1,1.0_q, &
                            FFM%C00(1,IHEAD,INTYP0),FFM%MMVAR1, &
                            C00_COUPLE(1,1),MROW_C00_COUPLE, &
                            0.0_q,SOAP1_COUPLE(1,1),MROW_SOAP)
                    ENDIF
                 ENDIF
#endif
               
                 CALL POLY_SOAP1 (IPOLY,FFM%LSUPERVEC,MROW_SOAP,MCOL_SOAP,NATOM_NUM_HELP, &
                                  NB_USE,NCOL_NITYP,FFM%NHYP1,SOAP1(1:,1:),SOAP2(1:,1:),FFM%W1)
                 IF (CALL_TYPE .EQ. "PRED") THEN
                    IF (FFM%LCOUPLE) THEN
                       CALL POLY_SOAP_COUPLE1 (IPOLY,FFM%LSUPERVEC,MROW_SOAP,MCOL_SOAP,MROW_SOAP_COUPLE, &
                                               MCOL_SOAP_COUPLE,NB_USE,NCOL_NITYP,FFM%NHYP1, &
                                               SOAP1(1:,1:),SOAP1_COUPLE(1:,1:),SOAP2_COUPLE(1:,1:),FFM%W1)
                    ENDIF
                 ENDIF
                 IPOLY=1
              ENDIF
! Store the calculated SOAP1 to SOAP_C00.
! MPI-parallel.
              IF(.NOT.FFM%LSUPERVEC) THEN
                 IF (NB_USE.GT.0.AND.NCOL_NITYP.GT.0) THEN
                    DO IB=1, NB_USE
                       DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                          INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                          IF(INITYP.LE.NITYP(INTYP0)) THEN
                             INIONS=LNIONS(INITYP,INTYP0)
                             LOCVAR_COL=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                             MLSHM%SOAP_C00(INIONS,IB)=FFM%W1*REAL(FFM%NHYP1,q)*(SOAP1(IB,LOCVAR_COL))**(FFM%NHYP1-1)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDIF
! wait so that all shared memory arrays finish writing
#ifdef use_shmem
                 CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
                 !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#endif
! Then, all-reduce SOAP_C00 at the end of INTYP loop.
                 IF (INTYP.EQ.NTYP) THEN
#ifdef use_shmem
                    IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN 
                       CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_C00,MLSHM%MROW_SOAP_C00_GLOBAL*NB_MAX,MPI_REAL8,MPI_SUM, &
                            PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                    ENDIF
                    CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
                    !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#else
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_C00,MLSHM%MROW_SOAP_C00_GLOBAL*NB_MAX, &
                                       MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
#endif
                 ENDIF
              ENDIF
           ENDIF
! Calculate contribution from 2-point correlation function
           IF (FFM%W2.GT.0.0_q) THEN
              IF (NB_USE.GT.0.AND.NCOL_NITYP.GT.0) THEN
#ifdef scaLAPACK
                 CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR2(INTYP0),1.0_q, &
                        FFM%PS_ALL(1,IHEAD),FFM%MMVAR2, &
                        PS(1,1,INTYP0),MROW_PS, &
                        0.0_q,SOAP1(1,1),MROW_SOAP)
                 IF (CALL_TYPE .EQ. "PRED") THEN
                    IF (FFM%LCOUPLE) THEN
                       CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR2(INTYP0),1.0_q, &
                            FFM%PS_ALL(1,IHEAD),FFM%MMVAR2, &
                            PS_COUPLE(1,1),MROW_PS_COUPLE, &
                            0.0_q,SOAP1_COUPLE(1,1),MROW_SOAP)
                    ENDIF
                 ENDIF
#else
                 CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR2(INTYP0),1.0_q, &
                        FFM%PS(1,IHEAD,INTYP0),FFM%MMVAR2, &
                        PS(1,1,INTYP0),MROW_PS, &
                        0.0_q,SOAP1(1,1),MROW_SOAP)
                 IF (CALL_TYPE .EQ. "PRED") THEN
                    IF (FFM%LCOUPLE) THEN
                       CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%NNVAR2(INTYP0),1.0_q, &
                                   FFM%PS(1,IHEAD,INTYP0),FFM%MMVAR2, &
                                   PS_COUPLE(1,1),MROW_PS_COUPLE, &
                                   0.0_q,SOAP1_COUPLE(1,1),MROW_SOAP)
                    ENDIF
                 ENDIF 
#endif
                 CALL POLY_SOAP1 (IPOLY,FFM%LSUPERVEC,MROW_SOAP,MCOL_SOAP,NATOM_NUM_HELP, &
                                  NB_USE,NCOL_NITYP,FFM%NHYP2,SOAP1(1:,1:),SOAP2(1:,1:),PREFAC_HELP)
                 IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                    IPOLY=1
                 ENDIF
                 IF (CALL_TYPE .EQ. "PRED") THEN
                    IF (FFM%LCOUPLE) THEN
                       CALL POLY_SOAP_COUPLE1 (IPOLY,FFM%LSUPERVEC,MROW_SOAP,MCOL_SOAP, &
                                               MROW_SOAP_COUPLE,MCOL_SOAP_COUPLE,NB_USE, &
                                               NCOL_NITYP,FFM%NHYP2,SOAP1(1:,1:), &
                                               SOAP1_COUPLE(1:,1:),SOAP2_COUPLE(1:,1:),FFM%W2)
                    ENDIF
                 ENDIF
              ENDIF
! Store the calculated SOAP1 to SOAP_PS.
! MPI-parallel.
              IF (.NOT.FFM%LSUPERVEC) THEN
                 DO IB=1, NB_USE
                    DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                       INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                       IF(INITYP.LE.NITYP(INTYP0)) THEN
                          INIONS=LNIONS(INITYP,INTYP0)
                          LOCVAR_COL=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                          MLSHM%SOAP_PS(INIONS,IB)=PREFAC_HELP*REAL(FFM%NHYP2,q)*(SOAP1(IB,LOCVAR_COL))**(FFM%NHYP2-1)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF
              ! Optional calculation of reduced descriptor
              IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
                 IF (NB_USE.GT.0.AND.NCOL_NITYP.GT.0) THEN
#ifdef scaLAPACK
                    CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%RV%NNVAR2_RED(INTYP0),1.0_q, &
                           FFM%RV%PS_ALL_RED(1,IHEAD),FFM%RV%MMVAR2_RED, &
                           FFM%RV%PS_UGANG_RED(1,1,INTYP0),FFM%RV%MMVAR2_RED, &
                           0.0_q,SOAP1(1,1),MROW_SOAP)
#else
                    CALL DGEMM ('T','N',NB_USE,NCOL_NITYP,FFM%RV%NNVAR2_RED(INTYP0),1.0_q, &
                           FFM%RV%PS_RED(1,1,INTYP0),FFM%RV%MMVAR2_RED, &
                           FFM%RV%PS_UGANG_RED(1,1,INTYP0),FFM%RV%MMVAR2_RED, &
                           0.0_q,SOAP1(1,1),MROW_SOAP)
#endif
                    CALL POLY_SOAP1 (IPOLY,FFM%LSUPERVEC,MROW_SOAP,MCOL_SOAP,NATOM_NUM_HELP, &
                                     NB_USE,NCOL_NITYP,FFM%NHYP2,SOAP1(1:,1:),SOAP2(1:,1:),PREFAC_HELP_RED)
                 ENDIF
                 IF (.NOT.FFM%LSUPERVEC) THEN
                    DO IB=1, NB_USE
                       DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                          INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                          IF (INITYP.LE.NITYP(INTYP0)) THEN
                             INIONS=LNIONS(INITYP,INTYP0)
                             LOCVAR_COL=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                             MLSHM%SOAP_PS_RED(INIONS,IB)=PREFAC_HELP_RED*REAL(FFM%NHYP2,q)*(SOAP1(IB,LOCVAR_COL))**(FFM%NHYP2-1)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDIF
              ENDIF
              IF (.NOT.FFM%LSUPERVEC) THEN
#ifdef use_shmem
                 CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
                 !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#endif
! Then, all-reduce SOAP_PS at the end of INTYP loop.
                 IF (INTYP.EQ.NTYP) THEN
#ifdef use_shmem
!                    IF (PAR_SUP_HANDLE%COMM_SHMEM%NODE_ME.EQ.0) THEN
                    IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                       CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_PS,MLSHM%MROW_SOAP_PS_GLOBAL*NB_MAX,MPI_REAL8, &
                            MPI_SUM,PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                    ENDIF
                    CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
                    !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#else
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_PS,MLSHM%MROW_SOAP_PS_GLOBAL*NB_MAX, &
                                       MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
#endif
                    ! Optionally do the same for reduced descriptor
                    IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
#ifdef use_shmem
                       !IF (PAR_SUP_HANDLE%COMM_SHMEM%NODE_ME.EQ.0) THEN
                       IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                          CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_PS_RED,MLSHM%MROW_SOAP_PS_GLOBAL*NB_MAX,MPI_REAL8, &
                               MPI_SUM,PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                       ENDIF
                       CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
                       !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#else  
                       CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_PS_RED,MLSHM%MROW_SOAP_PS_GLOBAL*NB_MAX, &
                                          MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
#endif
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF

           IF (FFM%LSUPERVEC) THEN
! Store the calculated SOAP1 to SOAP_ALL.
! MPI-parallel.
              IF (NB_USE.GT.0.AND.NCOL_NITYP.GT.0) THEN
                 DO IB=1, NB_USE
                    DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                       INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                       IF (INITYP.LE.NITYP(INTYP0)) THEN
                          INIONS=LNIONS(INITYP,INTYP0)
                          LOCVAR_COL=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                          MLSHM%SOAP_ALL(INIONS,IB)=REAL(FFM%NHYP2,q)*(SOAP2(IB,LOCVAR_COL))**(FFM%NHYP2-1)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF
#ifdef use_shmem
              CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
              !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#endif
! Then, all-reduce SOAP_PS at the end of INTYP loop.
              IF (INTYP.EQ.NTYP) THEN
#ifdef use_shmem
                 !IF (PAR_SUP_HANDLE%COMM_SHMEM%NODE_ME.EQ.0) THEN
                 IF (PAR_SUP_HANDLE%COMM_ZERO%COLOR .EQ. 0) THEN
                    CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_ALL,MLSHM%MROW_SOAP_ALL_GLOBAL*NB_MAX,MPI_REAL8, &
                         MPI_SUM,PAR_SUP_HANDLE%COMM_ZERO%MPI_COMM,IERR)
                 ENDIF
                 CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
                 !CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_world%MPI_COMM,IERR)
#else
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,MLSHM%SOAP_ALL,MLSHM%MROW_SOAP_ALL_GLOBAL*NB_MAX, &
                                    MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
#endif
              ENDIF
! Calculate SAP similarity measure for super vector.
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF (FFM%LCOUPLE) THEN
                    CALL POLY_SOAP_COUPLE2 (MROW_SOAP,MCOL_SOAP,MROW_SOAP_COUPLE,MCOL_SOAP_COUPLE, &
                                            NB_USE,NCOL_NITYP,FFM%NHYP2,SOAP2(1:,1:),SOAP2_COUPLE(1:,1:))
                 ENDIF
              ENDIF
              CALL POLY_SOAP2 (MROW_SOAP,MCOL_SOAP,NATOM_NUM_HELP,NB_USE,NCOL_NITYP,FFM%NHYP2,SOAP2(1:,1:))
           ENDIF
           PROFILING_STOP('compute_soap')
        END SUBROUTINE COMPUTE_SOAP



! Subroutine to calculate energy components in the design-matrix and copy them to global arrays.

        SUBROUTINE COMPUTE_ENERGY_HELP(FF, FFM, PAR_SUP_HANDLE, ICONF, INTYP0, ENERGY_HELP, &
                                       LNIONS, MROW_SOAP, NCOL_NITYP, NIONS, &
                                       NITYP, SOAP2, VEC_LOC_ENERGY, CALL_TYPE, &
                                       MCOL_EATOM, MROW_EATOM, EATOM_HELP, &
                                       SOAP2_COUPLE, ENERGY_COUPLE_HELP)
           USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, PARALLEL_SUPER, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: SOAP2(:,:)
           REAL(q), INTENT(INOUT) :: ENERGY_HELP(:)
           REAL(q), INTENT(INOUT) :: VEC_LOC_ENERGY(:,:)
           INTEGER, INTENT(IN)   :: ICONF
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: LNIONS(:,:)
           INTEGER, INTENT(IN)   :: MROW_SOAP
           INTEGER, INTENT(IN)   :: NCOL_NITYP
           INTEGER, INTENT(IN)   :: NIONS
           INTEGER, INTENT(IN)   :: NITYP(:)
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Input stuff only needed for FMAT_NEW
           REAL(q), INTENT(INOUT) :: EATOM_HELP(:,:)
           INTEGER, INTENT(IN)   :: MCOL_EATOM
           INTEGER, INTENT(IN)   :: MROW_EATOM
! Thermodynamic integration stuff
           REAL(q), INTENT(INOUT) :: SOAP2_COUPLE(:,:)
           REAL(q), INTENT(INOUT) :: ENERGY_COUPLE_HELP(:)
! Local variables
           INTEGER               :: IB
           INTEGER               :: ICOL
           INTEGER               :: IERR
           INTEGER               :: INIONS
           INTEGER               :: INITYP
           INTEGER               :: INITYP_TMP
           INTEGER               :: IROW
           INTEGER               :: LOCVAR_COL
           INTEGER               :: LOCVAR_ROW
           INTEGER               :: TEST_COL
           INTEGER               :: TEST_ROW
           INTEGER               :: NB_USE

           PROFILING_START('compute_energy')
           ! Set the right NB here
           IF (CALL_TYPE .EQ. "PUTA") THEN
              NB_USE=FFM%NB_NEW(INTYP0)
           ELSE
              NB_USE=FFM%NB(INTYP0)
           ENDIF
           IF (NB_USE.GT.0.AND.NCOL_NITYP.GT.0) THEN
              IF (CALL_TYPE .EQ. "PRED") THEN
                 DO INITYP_TMP=1, NITYP(INTYP0), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                     INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                     IF (INITYP.LE.NITYP(INTYP0)) THEN
                        INIONS=LNIONS(INITYP,INTYP0)
                        LOCVAR_COL=CEILING(REAL(INITYP,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                        IF (FF%LTOTEN_SYSTEM) THEN
                           VEC_LOC_ENERGY(LOCVAR_COL,1)=FFM%RCOUPLE_IONS(INIONS)
                        ELSE
                           VEC_LOC_ENERGY(LOCVAR_COL,1)=FFM%RCOUPLE_IONS(INIONS)/REAL(NIONS,q)
                        ENDIF
                     ENDIF
                  ENDDO
              ELSE
                 VEC_LOC_ENERGY=1.0_q
              ENDIF
              CALL DGEMV ('N',NB_USE,NCOL_NITYP,1.0_q, &
                    SOAP2,MROW_SOAP,VEC_LOC_ENERGY(1,1), &
                    1,0.0_q,ENERGY_HELP(1),1)
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF (FFM%LCOUPLE) THEN
                    CALL DGEMV ('N',NB_USE,NCOL_NITYP,1.0_q, &
                         SOAP2_COUPLE,MROW_SOAP,VEC_LOC_ENERGY(1,1), &
                         1,0.0_q,ENERGY_COUPLE_HELP(1),1)
                 ENDIF
              ENDIF

           ELSE
              ENERGY_HELP=0.0_q
              IF (CALL_TYPE .EQ. "PRED") THEN
                 IF(FFM%LCOUPLE) THEN
                    ENERGY_COUPLE_HELP(1:NB_USE)=0.0_q
                 ENDIF
              ENDIF
           ENDIF
! All-reduce ENERGY_HELP.
           CALL MPI_ALLREDUCE(MPI_IN_PLACE,ENERGY_HELP,NB_USE,MPI_REAL8,MPI_SUM, &
                              PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
           IF (CALL_TYPE .EQ. "PRED") THEN
              IF (FFM%LCOUPLE) THEN
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,ENERGY_COUPLE_HELP(1),NB_USE, &
                                    MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              ENDIF
           ENDIF
           IF (CALL_TYPE .EQ. "INIT" .OR. CALL_TYPE .EQ. "PUTA") THEN
              ! Store the design matrix elements
              IF (NB_USE.GT.0) THEN
#ifdef scaLAPACK
                 TEST_COL=MOD((FFM%LTOTEN(ICONF)-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 IF (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.TEST_COL) THEN
                    LOCVAR_COL=GLOBAL2LOCAL(FFM%LTOTEN(ICONF),PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
                    LOCVAR_COL=FFM%LTOTEN(ICONF)
#endif
                    DO IB=1, NB_USE
#ifdef scaLAPACK
                       IF (CALL_TYPE .EQ. "INIT") THEN
                          TEST_ROW=MOD((IB-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                       ELSE IF (CALL_TYPE .EQ. "PUTA") THEN
                          TEST_ROW=MOD((IB+FFM%NB(INTYP0)-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                       ENDIF
                       IF (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                          CYCLE
                       ENDIF
                       IF (CALL_TYPE .EQ. "INIT") THEN
                          LOCVAR_ROW=GLOBAL2LOCAL(IB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
                       ELSE IF (CALL_TYPE .EQ. "PUTA") THEN
                          LOCVAR_ROW=GLOBAL2LOCAL(IB+FFM%NB(INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
                       ENDIF
#else
                       IF (CALL_TYPE .EQ. "INIT") THEN
                          LOCVAR_ROW=IB
                       ELSE IF (CALL_TYPE .EQ. "PUTA") THEN
                          LOCVAR_ROW=IB+FFM%NB(INTYP0)
                       ENDIF
#endif
                       FFM%FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)=ENERGY_HELP(IB)
                    ENDDO
#ifdef scaLAPACK
                 ENDIF
#endif
              ENDIF
           ELSE IF (CALL_TYPE .EQ. "PRED") THEN
              ! Store the atomic energy elements.
              IF (FFM%LEATOM.OR.FFM%LHEAT.OR.FFM%LCOUPLE) THEN
                 DO ICOL=1, MCOL_EATOM
                    DO IROW=1, MROW_EATOM
                       EATOM_HELP(IROW,ICOL)=SOAP2(IROW,ICOL)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
           PROFILING_STOP('compute_energy')
        END SUBROUTINE COMPUTE_ENERGY_HELP


! Subroutine to calculate derivatives of expansion coefficients by using fast method.

        SUBROUTINE COMPUTE_D1CLM(ASA, FF, FFM, PAR_SUP_HANDLE, DCLM1, DCLM2, &
                                      DPS_SIC, KNIONS, KNTYP0, &
                                      LADD, LNIONS_EST1, LNIONS_EST2, &
                                      MCOL_DCLM1, MCOL_DCLM2, MCOL_DPS_SIC, &
                                      MROW_DCLM1, MROW_DCLM2, &
                                      MROW_DPS_SIC, MROW_EST1, MROW_EST2, &
                                      MTYP, NNEIB_EST1, NNEIB_EST2, &
                                      NTYP, RMEM_HELP1, RMEM_HELP2, &
                                      XYZ_EST1, XYZ_EST2)
           USE ML_FF_STRUCT, ONLY: ASA_STUFF, PARALLEL_SUPER, ML_SHMEM, FF_PAR, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (ASA_STUFF)      :: ASA
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: DCLM1(:,:,:,:)
           REAL(q), INTENT(INOUT) :: DCLM2(:,:,:,:)
           REAL(q), INTENT(INOUT) :: DPS_SIC(:,:,:,:)
           REAL(q), INTENT(INOUT) :: RMEM_HELP1
           REAL(q), INTENT(INOUT) :: RMEM_HELP2
           INTEGER, INTENT(IN)   :: KNIONS
           INTEGER, INTENT(IN)   :: KNTYP0
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:)
           INTEGER, INTENT(IN)   :: MCOL_DCLM1
           INTEGER, INTENT(IN)   :: MCOL_DCLM2
           INTEGER, INTENT(IN)   :: MCOL_DPS_SIC
           INTEGER, INTENT(IN)   :: MROW_DCLM1
           INTEGER, INTENT(IN)   :: MROW_DCLM2
           INTEGER, INTENT(IN)   :: MROW_DPS_SIC
           INTEGER, INTENT(IN)   :: MROW_EST1
           INTEGER, INTENT(IN)   :: MROW_EST2
           INTEGER, INTENT(IN)   :: MTYP
           INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:)
           INTEGER, INTENT(IN)   :: NTYP
           REAL(q), INTENT(IN)   :: XYZ_EST1(:,:,:,:)
           REAL(q), INTENT(IN)   :: XYZ_EST2(:,:,:,:)
! Local variables
           INTEGER               :: LOCVAR_COL_EST1
           INTEGER               :: LOCVAR_COL_EST2
           PROFILING_START('compute_d1clm')
           IF(FFM%W1.GT.0.0_q) THEN
              LOCVAR_COL_EST1=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              CALL D1CLM_SOAP (DCLM1,FFM%DRSPL1,FFM%FNL1,FFM%FWIN1,KNIONS,KNTYP0, &
                   LADD(1:),FFM%LFLAG_VAR,FFM%LFLAG_VAR_SIC,FFM%LMAX1,LNIONS_EST1(1:,1:,LOCVAR_COL_EST1), &
                   .FALSE.,FFM%LVAR,FFM%LVAR_SIC,FFM%LWINDOW1, &
                   MCOL_DCLM1,MCOL_DPS_SIC,FF%MIONS,FFM%MRB1, &
                   MROW_DCLM1,MROW_DPS_SIC,MROW_EST1,FFM%MSPL1,MTYP, &
                   NNEIB_EST1(1:,LOCVAR_COL_EST1),FFM%NRB1(0:),FFM%NSPL1,NTYP, &
                   DPS_SIC,FFM%RCOUPLE_IONS,RMEM_HELP2,XYZ_EST1(1:,1:,1:,LOCVAR_COL_EST1), &
                   ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
              RMEM_HELP1=MAX(RMEM_HELP1,RMEM_HELP2)
           ENDIF
           IF(FFM%W2.GT.0.0_q) THEN
              LOCVAR_COL_EST2=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7) THEN
                 CALL D1CLM_SOAP_NOSIC_DUAL (FFM%RV%DCLM2_RED,DCLM2,FFM%DRSPL2, &
                      FFM%FNL2,FFM%RV%FNL2_RED,FFM%FWIN2,KNIONS,KNTYP0, &
                      LADD(1:),FFM%LMAX2,FFM%RV%LMAX2_RED,LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      FFM%LWINDOW2, &
                      MCOL_DCLM2,FF%MIONS,FFM%MRB2,FFM%RV%MRB2_RED, &
                      MROW_DCLM2,MROW_EST2,FFM%MSPL2,MTYP, &
                      NNEIB_EST2(1:,LOCVAR_COL_EST2),FFM%NRB2(0:),FFM%RV%NRB2_RED(0:),FFM%NSPL2,NTYP, &
                      FFM%RCOUPLE_IONS,RMEM_HELP2,XYZ_EST2(1:,1:,1:,LOCVAR_COL_EST2), &
                      ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
              ELSE
                 CALL D1CLM_SOAP (DCLM2,FFM%DRSPL2,FFM%FNL2,FFM%FWIN2,KNIONS,KNTYP0, &
                      LADD(1:),FFM%LFLAG_VAR,FFM%LFLAG_VAR_SIC, &
                      FFM%LMAX2,LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      FFM%LSIC,FFM%LVAR,FFM%LVAR_SIC,FFM%LWINDOW2, &
                      MCOL_DCLM2,MCOL_DPS_SIC,FF%MIONS,FFM%MRB2, &
                      MROW_DCLM2,MROW_DPS_SIC,MROW_EST2,FFM%MSPL2,MTYP, &
                      NNEIB_EST2(1:,LOCVAR_COL_EST2),FFM%NRB2(0:),FFM%NSPL2,NTYP, &
                      DPS_SIC,FFM%RCOUPLE_IONS,RMEM_HELP2,XYZ_EST2(1:,1:,1:,LOCVAR_COL_EST2), &
                      ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
              ENDIF
              RMEM_HELP1=MAX(RMEM_HELP1,RMEM_HELP2)
           ENDIF
           PROFILING_STOP('compute_d1clm')
        END SUBROUTINE COMPUTE_D1CLM


        SUBROUTINE COMPUTE_D1C00(FFM, MLSHM, PAR_SUP_HANDLE, DC00, DCLM1, &
                                      KNIONS, KNTYP0, LADD, LNIONS_EST1, &
                                      LNIONS_EST_MAX, MCOL_DC00, MCOL_DCLM1, &
                                      MROW_DC00, MROW_DCLM1, MROW_EST1, MROW_EST_MAX, &
                                      MTYP, MTYP_ALLOCATE, &
                                      NNEIB_EST1, NNEIB_EST_MAX, NTYP)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, ML_SHMEM, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: DC00(:,:,0:,:)
           REAL(q), INTENT(INOUT) :: DCLM1(:,:,:,:)
           INTEGER, INTENT(IN)   :: KNIONS
           INTEGER, INTENT(IN)   :: KNTYP0
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST_MAX(:,:,:)
           INTEGER, INTENT(IN)   :: MCOL_DC00
           INTEGER, INTENT(IN)   :: MCOL_DCLM1
           INTEGER, INTENT(IN)   :: MROW_DC00
           INTEGER, INTENT(IN)   :: MROW_DCLM1
           INTEGER, INTENT(IN)   :: MROW_EST1
           INTEGER, INTENT(IN)   :: MROW_EST_MAX
           INTEGER, INTENT(IN)   :: MTYP
           INTEGER, INTENT(IN)   :: MTYP_ALLOCATE
           INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST_MAX(:,:)
           INTEGER, INTENT(IN)   :: NTYP
! Local variables
           INTEGER :: LOCVAR_COL_EST1
           INTEGER :: LOCVAR_COL_EST_MAX
           PROFILING_START('compute_d1c00')
           IF (FFM%W1 .GT. 0.0_q) THEN
              LOCVAR_COL_EST1=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              IF (.NOT.FFM%LSUPERVEC) THEN
                 LOCVAR_COL_EST_MAX=1
              ELSE
                 LOCVAR_COL_EST_MAX=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              ENDIF
              CALL D1C00_SOAP (MLSHM%CLM1,DC00,DCLM1,FFM%LSUPERVEC,KNTYP0, &
                   LADD(1:),FFM%LMAX1,LNIONS_EST1(1:,1:,LOCVAR_COL_EST1), &
                   LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                   MLSHM%MCOL_CLM1_GLOBAL,MCOL_DC00,MCOL_DCLM1,FFM%MRB1,MROW_EST1, &
                   MROW_EST_MAX,MLSHM%MROW_CLM1_GLOBAL,MROW_DC00,MROW_DCLM1, &
                   MTYP,MTYP_ALLOCATE,NNEIB_EST1(1:,LOCVAR_COL_EST1), &
                   NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX),FFM%NRB1(0:),NTYP)
           ENDIF
           PROFILING_STOP('compute_d1c00')
        END SUBROUTINE COMPUTE_D1C00


! Subroutine to calculate derivatives of angular descriptors.

        SUBROUTINE COMPUTE_D1PS(FFM, MLSHM, PAR_SUP_HANDLE, DCLM2, &
                                     DPS, DPS_SIC, &
                                     KNIONS, KNTYP0, LADD, LNIONS_EST2, LNIONS_EST_MAX, &
                                     MCOL_DCLM2, MCOL_DPS, MCOL_DPS_SIC, &
                                     MROW_EST2, MROW_EST_MAX, MROW_DCLM2, &
                                     MROW_DPS, MROW_DPS_SIC, &
                                     MTYP, MTYP_ALLOCATE, NIONS, &
                                     NNEIB_EST2, NNEIB_EST_MAX, NTYP)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, ML_SHMEM, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: DCLM2(:,:,:,:)
           REAL(q), INTENT(INOUT) :: DPS(:,:,0:,:)
           REAL(q), INTENT(INOUT) :: DPS_SIC(:,:,:,:)
           INTEGER, INTENT(IN)   :: KNIONS
           INTEGER, INTENT(IN)   :: KNTYP0
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST_MAX(:,:,:)
           INTEGER, INTENT(IN)   :: MCOL_DCLM2
           INTEGER, INTENT(IN)   :: MCOL_DPS
           INTEGER, INTENT(IN)   :: MCOL_DPS_SIC
           INTEGER, INTENT(IN)   :: MROW_EST2
           INTEGER, INTENT(IN)   :: MROW_EST_MAX
           INTEGER, INTENT(IN)   :: MROW_DCLM2
           INTEGER, INTENT(IN)   :: MROW_DPS
           INTEGER, INTENT(IN)   :: MROW_DPS_SIC
           INTEGER, INTENT(IN)   :: MTYP
           INTEGER, INTENT(IN)   :: MTYP_ALLOCATE
           INTEGER, INTENT(IN)   :: NIONS
           INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST_MAX(:,:)
           INTEGER, INTENT(IN)   :: NTYP
! Local variables
           INTEGER :: LOCVAR_COL_EST2
           INTEGER :: LOCVAR_COL_EST_MAX
           PROFILING_START('compute_d1ps')
           IF (FFM%W2 .GT. 0.0_q) THEN
              LOCVAR_COL_EST2=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              IF (.NOT.FFM%LSUPERVEC) THEN
                 LOCVAR_COL_EST_MAX=1
              ELSE
                 LOCVAR_COL_EST_MAX=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              ENDIF
              IF (FFM%DESC_TYPE.EQ.0) THEN
                 CALL D1PS_SOAP (FFM%AFILT2, MLSHM%CLM2, DCLM2, DPS,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%LFLAG_VAR, &
                      FFM%LFLAG_VAR_SIC,FFM%L_LNRB2,FFM%LMAX2, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      FFM%LSIC,FFM%LVAR,FFM%LVAR_SIC, &
                      FFM%MLNRB2,FFM%MRB2, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%NLNRB2,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%NRB2(0:),FFM%NRB_LNRB2,NTYP,MLSHM%PS_SIC,FFM%WVAR, &
                      FFM%LMAX2)
              ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
                 CALL D1PS_SOAP_NOSIC_LIN_ELEM (FFM%AFILT2, MLSHM%CLM2, DCLM2, DPS,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%LFLAG_VAR, &
                      FFM%LFLAG_VAR_SIC,FFM%L_LNRB2,FFM%LMAX2, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      FFM%LSIC,FFM%LVAR,FFM%LVAR_SIC, &
                      MLSHM%MCOL_CLM2_GLOBAL,MCOL_DCLM2,MCOL_DPS, &
                      MCOL_DPS_SIC,MLSHM%MCOL_PS_SIC_GLOBAL,FFM%MLNRB2,FFM%MRB2, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MROW_DPS,MROW_DPS_SIC,MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%NLNRB2,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%NRB2(0:),FFM%NRB_LNRB2,NTYP,MLSHM%PS_SIC,FFM%WVAR)
              ELSE IF (FFM%DESC_TYPE.EQ.2) THEN
                 ! High resolution type-independent descriptor
                 CALL D1PS_SOAP_NOSIC_NO_ELEM (FFM%AFILT2, MLSHM%CLM2, DCLM2, DPS,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%LFLAG_VAR, &
                      FFM%LFLAG_VAR_SIC,FFM%L_LNRB2,FFM%LMAX2, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      .FALSE.,FFM%LVAR,FFM%LVAR_SIC, &
                      MLSHM%MCOL_CLM2_GLOBAL,MCOL_DCLM2,MCOL_DPS, &
                      MCOL_DPS_SIC,MLSHM%MCOL_PS_SIC_GLOBAL,FFM%MLNRB2,FFM%MRB2, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MROW_DPS,MROW_DPS_SIC,MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%NLNRB2,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%NRB2(0:),FFM%NRB_LNRB2,NTYP,MLSHM%PS_SIC,FFM%WVAR)
                 ! Low resolution type-dependent descriptor
                 CALL D1PS_SOAP (FFM%AFILT2, MLSHM%CLM2_RED, FFM%RV%DCLM2_RED, FFM%RV%DPS_UGANG_RED,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%RV%LFLAG_VAR_RED, &
                      FFM%LFLAG_VAR_SIC,FFM%RV%L_LNRB2_RED,FFM%RV%LMAX2_RED, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      .FALSE.,FFM%RV%LVAR_RED,FFM%LVAR_SIC, &
                      FFM%RV%MLNRB2_RED,FFM%RV%MRB2_RED, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%RV%NLNRB2_RED,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%RV%NRB2_RED(0:),FFM%RV%NRB_LNRB2_RED,NTYP,MLSHM%PS_SIC,FFM%RV%WVAR_RED, &
                      FFM%RV%LMAX2_RED)
              ELSE IF (FFM%DESC_TYPE.EQ.3) THEN
                 CALL D1PS_SOAP_LIN_MRB (FFM%AFILT2, MLSHM%CLM2, DCLM2, DPS,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%LFLAG_VAR, &
                      FFM%LFLAG_VAR_SIC,FFM%L_LNRB2,FFM%LMAX2, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      FFM%LSIC,FFM%LVAR,FFM%LVAR_SIC, &
                      FFM%MLNRB2,FFM%MRB2, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%NLNRB2,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%NRB2(0:),FFM%NRB_LNRB2,NTYP,MLSHM%PS_SIC,FFM%WVAR, &
                      FFM%LMAX2)
              ELSE IF (FFM%DESC_TYPE.EQ.4) THEN
                 ! First descriptor
                 CALL D1PS_SOAP (FFM%AFILT2, MLSHM%CLM2, DCLM2, DPS,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%LFLAG_VAR, &
                      FFM%LFLAG_VAR_SIC,FFM%L_LNRB2,FFM%LMAX2, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      FFM%LSIC,FFM%LVAR,FFM%LVAR_SIC, &
                      FFM%MLNRB2,FFM%MRB2, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%NLNRB2,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%NRB2(0:),FFM%NRB_LNRB2,NTYP,MLSHM%PS_SIC,FFM%WVAR, &
                      FFM%LMAX2)
                 ! Second descriptor
                 CALL D1PS_SOAP (FFM%AFILT2, MLSHM%CLM2_RED, FFM%RV%DCLM2_RED, FFM%RV%DPS_UGANG_RED,DPS_SIC, &
                      FFM%IAFILT2,FFM%LSUPERVEC,KNTYP0, &
                      FFM%LAFILT2,LADD(1:),FFM%RV%LFLAG_VAR_RED, &
                      FFM%LFLAG_VAR_SIC,FFM%RV%L_LNRB2_RED,FFM%RV%LMAX2_RED, &
                      LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                      LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                      .FALSE.,FFM%RV%LVAR_RED,FFM%LVAR_SIC, &
                      FFM%RV%MLNRB2_RED,FFM%RV%MRB2_RED, &
                      MROW_EST2,MROW_EST_MAX,MLSHM%MROW_CLM2_GLOBAL,MROW_DCLM2, &
                      MLSHM%MROW_PS_SIC_GLOBAL,MTYP,MTYP_ALLOCATE, &
                      NIONS,FFM%RV%NLNRB2_RED,NNEIB_EST2(1:,LOCVAR_COL_EST2),NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                      FFM%RV%NRB2_RED(0:),FFM%RV%NRB_LNRB2_RED,NTYP,MLSHM%PS_SIC,FFM%RV%WVAR_RED, &
                      FFM%RV%LMAX2_RED)
              ENDIF
           ENDIF
           PROFILING_STOP('compute_d1ps')
        END SUBROUTINE COMPUTE_D1PS

! Normalization.

        SUBROUTINE COMPUTE_NORMD1(FF, FFM, PAR_SUP_HANDLE, DC00, DPS, KNIONS, &
                                                 KNTYP0, LADD, LNIONS_EST1, LNIONS_EST2, LNIONS_EST_MAX, &
                                                 LNITYP, MCOL_DC00, MCOL_DPS, MCOL_RNORM_ALL, &
                                                 MCOL_RNORM_C00, MCOL_RNORM_PS, MIONS, MITYP, &
                                                 MROW_EST1, MROW_EST2, MROW_EST_MAX, MROW_DC00, &
                                                 MROW_DPS, MROW_RNORM_C00, MROW_RNORM_PS, &
                                                 NNEIB_EST1, NNEIB_EST2, NNEIB_EST_MAX, NTYP, &
                                                 RNORM_ALL, RNORM_C00, RNORM_PS)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FF_PAR, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
! Input and output varables
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: DC00(:,:,0:,:)
           REAL(q), INTENT(INOUT) :: DPS(:,:,0:,:)
           REAL(q), INTENT(INOUT) :: RNORM_ALL(:,:)
           REAL(q), INTENT(INOUT) :: RNORM_C00(:,:)
           REAL(q), INTENT(INOUT) :: RNORM_PS(:,:)
           INTEGER, INTENT(IN)   :: KNIONS
           INTEGER, INTENT(IN)   :: KNTYP0
           INTEGER, INTENT(IN)   :: LADD(:)
           INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST_MAX(:,:,:)
           INTEGER, INTENT(IN)   :: LNITYP(:)
           INTEGER, INTENT(IN)   :: MCOL_DC00
           INTEGER, INTENT(IN)   :: MCOL_DPS
           INTEGER, INTENT(IN)   :: MCOL_RNORM_ALL
           INTEGER, INTENT(IN)   :: MCOL_RNORM_C00
           INTEGER, INTENT(IN)   :: MCOL_RNORM_PS
           INTEGER, INTENT(IN)   :: MIONS
           INTEGER, INTENT(IN)   :: MITYP
           INTEGER, INTENT(IN)   :: MROW_EST1
           INTEGER, INTENT(IN)   :: MROW_EST2
           INTEGER, INTENT(IN)   :: MROW_EST_MAX
           INTEGER, INTENT(IN)   :: MROW_DC00
           INTEGER, INTENT(IN)   :: MROW_DPS
           INTEGER, INTENT(IN)   :: MROW_RNORM_C00
           INTEGER, INTENT(IN)   :: MROW_RNORM_PS
           INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST_MAX(:,:)
           INTEGER, INTENT(IN)   :: NTYP
! Local variables
           INTEGER :: LOCVAR_COL_EST1 
           INTEGER :: LOCVAR_COL_EST2
           INTEGER :: LOCVAR_COL_EST_MAX
           PROFILING_START('compute_normd1')
           IF (.NOT.FFM%LSUPERVEC) THEN
              IF (FFM%W1.GT.0.0_q) THEN
                 LOCVAR_COL_EST1=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              ELSE
                 LOCVAR_COL_EST1=1
              ENDIF
              IF (FFM%W2.GT.0.0_q) THEN
                 LOCVAR_COL_EST2=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              ELSE
                 LOCVAR_COL_EST2=1
              ENDIF
              LOCVAR_COL_EST_MAX=1
           ELSE
              IF (FFM%W1.GT.0.0_q) THEN
                 LOCVAR_COL_EST1=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              ELSE
                 LOCVAR_COL_EST1=1
              ENDIF
              IF (FFM%W2.GT.0.0_q) THEN
                 LOCVAR_COL_EST2=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
              ELSE
                 LOCVAR_COL_EST2=1
              ENDIF
              LOCVAR_COL_EST_MAX=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
           ENDIF
           CALL NORMALIZATION_D1 (FFM,DC00,DPS,FFM%RV%DPS_UGANG_RED, &
                                  FFM%LSUPERVEC,KNTYP0,LADD(1:), &
                                  LNIONS_EST1(1:,1:,LOCVAR_COL_EST1), &
                                  LNIONS_EST2(1:,1:,LOCVAR_COL_EST2), &
                                  LNIONS_EST_MAX(1:,1:,LOCVAR_COL_EST_MAX), &
                                  LNITYP(1:),FFM%LNORM1,FFM%LNORM2, &
                                  MROW_EST1,MROW_EST2,FF%MTYP_ALLOCATE, &
                                  NNEIB_EST1(1:,LOCVAR_COL_EST1),NNEIB_EST2(1:,LOCVAR_COL_EST2), &
                                  NNEIB_EST_MAX(1:,LOCVAR_COL_EST_MAX), &
                                  FFM%NNVAR1,FFM%NNVAR2, &
                                  FFM%RV%NNVAR2_RED,NTYP, &
                                  FFM%RCUT1,FFM%RCUT2, &
                                  RNORM_ALL(1:,1:),RNORM_C00(1:,1:), &
                                  RNORM_PS(1:,1:),FFM%RV%RNORM_PS_RED(1:,1:), &
                                  FFM%W1,FFM%W2)
           PROFILING_STOP('compute_normd1')
        END SUBROUTINE COMPUTE_NORMD1


! Subroutine to calculate derivatives of SOAP similarity kernel.

        SUBROUTINE COMPUTE_DSOAP(FFM, MLSHM, PAR_SUP_HANDLE, DC00, DPS, &
                                             DSOAP11, DSOAP12, DSOAP21, DSOAP22, &
                                             INTYP0, KNIONS, KNTYP0, LHEAD, &
                                             LNIONS_EST1, LNIONS_EST2, LNIONS_EST_MAX, &
                                             MCOL_DSOAP1, MCOL_DSOAP2, &
                                             MROW_DC00, MROW_DPS, MROW_EST1, MROW_EST2, &
                                             MROW_EST_MAX, MROW_DSOAP1, MROW_DSOAP2, &
                                             NNEIB_EST1, NNEIB_EST2, NNEIB_EST_MAX, &
                                             CALL_TYPE)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, ML_SHMEM, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: DC00(:,:,0:,:)
           REAL(q), INTENT(INOUT) :: DPS(:,:,0:,:)
           REAL(q), INTENT(INOUT) :: DSOAP11(:,:)
           REAL(q), INTENT(INOUT) :: DSOAP12(:,:)
           REAL(q), INTENT(INOUT) :: DSOAP21(:,:,:)
           REAL(q), INTENT(INOUT) :: DSOAP22(:,:,:)
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: KNIONS
           INTEGER, INTENT(IN)   :: KNTYP0
           INTEGER, INTENT(IN)   :: LHEAD(:)
           INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:)
           INTEGER, INTENT(IN)   :: LNIONS_EST_MAX(:,:,:)
           INTEGER, INTENT(IN)   :: MCOL_DSOAP1
           INTEGER, INTENT(IN)   :: MCOL_DSOAP2
           INTEGER, INTENT(IN)   :: MROW_DC00
           INTEGER, INTENT(IN)   :: MROW_DPS
           INTEGER, INTENT(IN)   :: MROW_EST1
           INTEGER, INTENT(IN)   :: MROW_EST2
           INTEGER, INTENT(IN)   :: MROW_EST_MAX
           INTEGER, INTENT(IN)   :: MROW_DSOAP1
           INTEGER, INTENT(IN)   :: MROW_DSOAP2
           INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST_MAX(:,:)
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Local variables
           INTEGER               :: IB
           INTEGER               :: IHEAD
           INTEGER               :: IXYZ
           INTEGER               :: LOCVAR_COL
           INTEGER               :: NNEIB_EST_DUMMY
           INTEGER               :: NB_USE
           PROFILING_START('compute_dsoap')
! Calculate header.
           IF (CALL_TYPE .EQ. "PUTA") THEN
#ifdef scaLAPACK
              IHEAD=FFM%NB(INTYP0)+LHEAD(INTYP0)
#else
              IHEAD=FFM%NB(INTYP0)+1
#endif
           ELSE
#ifdef scaLAPACK
              IHEAD=LHEAD(INTYP0)
#else
              IHEAD=1
#endif
           ENDIF
           ! Determine size of used local basis configs
           IF (CALL_TYPE .EQ. "PUTA") THEN
              NB_USE=FFM%NB_NEW(INTYP0)
           ELSE
              NB_USE=FFM%NB(INTYP0)
           ENDIF
! Calculate LOCVAR_COL.
           LOCVAR_COL=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
! Calculate derivatives of 1-point correlation kernel
           IF (FFM%W1 .GT. 0.0_q) THEN
              IF (.NOT.FFM%LSUPERVEC) THEN
                 NNEIB_EST_DUMMY=NNEIB_EST1(INTYP0,LOCVAR_COL)
              ELSE
                 NNEIB_EST_DUMMY=NNEIB_EST_MAX(INTYP0,LOCVAR_COL)
              ENDIF
              IF (NNEIB_EST_DUMMY+1.GT.0.AND.NB_USE.GT.0) THEN
                 DO IXYZ=1, 3
#ifdef scaLAPACK
                    CALL DGEMM ('T','N',NNEIB_EST_DUMMY+1,NB_USE,FFM%NNVAR1,1.0_q, &
                         DC00(1,1,IXYZ,INTYP0),MROW_DC00,FFM%C00_ALL(1,IHEAD),FFM%MMVAR1, &
                         0.0_q,DSOAP11(1,1),MROW_DSOAP1)
#else
                    CALL DGEMM ('T','N',NNEIB_EST_DUMMY+1,NB_USE,FFM%NNVAR1,1.0_q, &
                         DC00(1,1,IXYZ,INTYP0),MROW_DC00,FFM%C00(1,IHEAD,INTYP0),FFM%MMVAR1, &
                         0.0_q,DSOAP11(1,1),MROW_DSOAP1)
#endif
                    IF (.NOT.FFM%LSUPERVEC) THEN
                       CALL POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST1(1:,INTYP0,LOCVAR_COL), &
                                        NNEIB_EST1(INTYP0,LOCVAR_COL)+1, &
                                        NB_USE,DSOAP11(1:,1:),DSOAP21(1:,1:,IXYZ), &
                                        MLSHM%SOAP_C00)
                    ELSE
                       CALL POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST_MAX(1:,INTYP0,LOCVAR_COL), &
                                        NNEIB_EST_MAX(INTYP0,LOCVAR_COL)+1, &
                                        NB_USE,DSOAP11(1:,1:),DSOAP21(1:,1:,IXYZ), &
                                        MLSHM%SOAP_ALL)
                    ENDIF
                 ENDDO
              ENDIF
           ENDIF
! Calculate derivatives of 2-point correlation kernel
           IF (FFM%W2 .GT. 0.0_q) THEN
              IF (.NOT.FFM%LSUPERVEC) THEN
                 NNEIB_EST_DUMMY=NNEIB_EST2(INTYP0,LOCVAR_COL)
              ELSE
                 NNEIB_EST_DUMMY=NNEIB_EST_MAX(INTYP0,LOCVAR_COL)
              ENDIF
              IF (NNEIB_EST_DUMMY+1.GT.0.AND.NB_USE.GT.0) THEN
                 DO IXYZ=1, 3
#ifdef scaLAPACK
                    CALL DGEMM ('T','N',NNEIB_EST_DUMMY+1,NB_USE,FFM%NNVAR2(INTYP0),1.0_q, &
                         DPS(1,1,IXYZ,INTYP0),MROW_DPS,FFM%PS_ALL(1,IHEAD),FFM%MMVAR2, &
                         0.0_q,DSOAP12(1,1),MROW_DSOAP2)
#else
                    CALL DGEMM ('T','N',NNEIB_EST_DUMMY+1,NB_USE,FFM%NNVAR2(INTYP0),1.0_q, &
                         DPS(1,1,IXYZ,INTYP0),MROW_DPS,FFM%PS(1,IHEAD,INTYP0),FFM%MMVAR2, &
                         0.0_q,DSOAP12(1,1),MROW_DSOAP2)
#endif
                    IF (.NOT.FFM%LSUPERVEC) THEN
                       CALL POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST2(1:,INTYP0,LOCVAR_COL), &
                                        NNEIB_EST2(INTYP0,LOCVAR_COL)+1, &
                                        NB_USE,DSOAP12(1:,1:),DSOAP22(1:,1:,IXYZ), &
                                        MLSHM%SOAP_PS)
                    ELSE
                       CALL POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST_MAX(1:,INTYP0,LOCVAR_COL), &
                                        NNEIB_EST_MAX(INTYP0,LOCVAR_COL)+1, &
                                        NB_USE,DSOAP12(1:,1:),DSOAP22(1:,1:,IXYZ), &
                                        MLSHM%SOAP_ALL)
                    ENDIF
                    ! Optional reduced descriptor
                    IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
#ifdef scaLAPACK
                       CALL DGEMM ('T','N',NNEIB_EST_DUMMY+1,NB_USE,FFM%RV%NNVAR2_RED(INTYP0),1.0_q, &
                            FFM%RV%DPS_UGANG_RED(1,1,IXYZ,INTYP0),FFM%RV%MMVAR2_RED, &
                            FFM%RV%PS_ALL_RED(1,IHEAD),FFM%RV%MMVAR2_RED, &
                            0.0_q,DSOAP12(1,1),MROW_DSOAP2)
#else
                       CALL DGEMM ('T','N',NNEIB_EST_DUMMY+1,NB_USE,FFM%RV%NNVAR2_RED(INTYP0),1.0_q, &
                            FFM%RV%DPS_UGANG_RED(1,1,IXYZ,INTYP0),FFM%RV%MMVAR2_RED, &
                            FFM%RV%PS_RED(1,1,INTYP0),FFM%RV%MMVAR2_RED, &
                            0.0_q,DSOAP12(1,1),MROW_DSOAP2)
#endif
                       IF (.NOT.FFM%LSUPERVEC) THEN
                          ! Why POLY_DSOAP_RED instead of POLY_DSOAP ????
                          !CALL POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST2(1:,INTYP0,LOCVAR_COL), &
                          CALL POLY_DSOAP_RED (INTYP0,KNIONS,KNTYP0,LNIONS_EST2(1:,INTYP0,LOCVAR_COL), &
                                           NNEIB_EST2(INTYP0,LOCVAR_COL)+1, &
                                           NB_USE,DSOAP12(1:,1:),DSOAP22(1:,1:,IXYZ), &
                                           MLSHM%SOAP_PS_RED)
                       ELSE
                          !CALL POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST_MAX(1:,INTYP0,LOCVAR_COL), &
                          CALL POLY_DSOAP_RED (INTYP0,KNIONS,KNTYP0,LNIONS_EST_MAX(1:,INTYP0,LOCVAR_COL), &
                                           NNEIB_EST_MAX(INTYP0,LOCVAR_COL)+1, &
                                           NB_USE,DSOAP12(1:,1:),DSOAP22(1:,1:,IXYZ), &
                                           MLSHM%SOAP_ALL)
                       ENDIF
                    ENDIF
                 ENDDO
              ENDIF
           ENDIF
           PROFILING_STOP('compute_dsoap')
        END SUBROUTINE COMPUTE_DSOAP

! Subroutine to calculate force and stress components in the design matrix by using fast method.
! Copying subroutine is seperately written later.

        SUBROUTINE COMPUTE_FORCE_STRESS(FFM, PAR_SUP_HANDLE, DSOAP21, DSOAP22, &
                                                 FORCE_HELP1, INTYP0, KNIONS, KNTYP0, &
                                                 MROW_DSOAP1, MROW_DSOAP2, &
                                                 NNEIB_EST1, NNEIB_EST2, NNEIB_EST_MAX, &
                                                 OMEGA, STRESS_HELP, &
                                                 VEC_LOC_FORCE1, VEC_LOC_FORCE2, &
                                                 VEC_LOC_STRESS1, VEC_LOC_STRESS2, &
                                                 XYZ_EST1, XYZ_EST2, XYZ_EST_MAX, &
                                                 CALL_TYPE)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (FFM_PAR)        :: FFM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: DSOAP21(:,:,:)
           REAL(q), INTENT(INOUT) :: DSOAP22(:,:,:)
           REAL(q), INTENT(INOUT) :: FORCE_HELP1(:,:,:,:)
           REAL(q), INTENT(INOUT) :: STRESS_HELP(:,:,:,:)
           REAL(q), INTENT(INOUT) :: VEC_LOC_FORCE1(:,:)
           REAL(q), INTENT(INOUT) :: VEC_LOC_FORCE2(:,:)
           REAL(q), INTENT(INOUT) :: VEC_LOC_STRESS1(:,:)
           REAL(q), INTENT(INOUT) :: VEC_LOC_STRESS2(:,:)
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: KNIONS
           INTEGER, INTENT(IN)   :: KNTYP0
           INTEGER, INTENT(IN)   :: MROW_DSOAP1
           INTEGER, INTENT(IN)   :: MROW_DSOAP2
           INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:)
           INTEGER, INTENT(IN)   :: NNEIB_EST_MAX(:,:)
           REAL(q), INTENT(IN)   :: OMEGA
           REAL(q), INTENT(IN)   :: XYZ_EST1(:,:,:,:)
           REAL(q), INTENT(IN)   :: XYZ_EST2(:,:,:,:)
           REAL(q), INTENT(IN)   :: XYZ_EST_MAX(:,:,:,:)
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE            
! Local variables
           INTEGER               :: INEIB
           INTEGER               :: IXYZ
           INTEGER               :: JXYZ
           INTEGER               :: LOCVAR_COL
           INTEGER               :: NNEIB_EST_DUMMY
           INTEGER               :: NB_USE
           PROFILING_START('compute_force_stress')
           ! Determine size of used local basis configs
           IF (CALL_TYPE .EQ. "PUTA") THEN
              NB_USE=FFM%NB_NEW(INTYP0)
           ELSE
              NB_USE=FFM%NB(INTYP0)
           ENDIF
           LOCVAR_COL=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
! Contributions from radial descriptor.
           IF (FFM%W1.GT.0.0_q) THEN
              IF(.NOT.FFM%LSUPERVEC) THEN
                 NNEIB_EST_DUMMY=NNEIB_EST1(INTYP0,LOCVAR_COL)
              ELSE
                 NNEIB_EST_DUMMY=NNEIB_EST_MAX(INTYP0,LOCVAR_COL)
              ENDIF
              IF (NNEIB_EST_DUMMY+1.GT.0.AND.NB_USE.GT.0) THEN
                 DO IXYZ=1, 3
! Calculate force components.
                    DO INEIB=1, NNEIB_EST_DUMMY+1
                       IF ((INEIB.EQ.1).AND.(INTYP0.NE.KNTYP0)) THEN
                          VEC_LOC_FORCE1(INEIB,1)=0.0_q
                       ELSE
                          VEC_LOC_FORCE1(INEIB,1)=1.0_q
                       ENDIF
                    ENDDO
                    CALL DGEMV ('T',NNEIB_EST_DUMMY+1,NB_USE,1.0_q,DSOAP21(1,1,IXYZ),MROW_DSOAP1, &
                         VEC_LOC_FORCE1(1,1),1,1.0_q,FORCE_HELP1(1,IXYZ,LOCVAR_COL,INTYP0),1)
! Calculate stress components.
                    DO JXYZ=IXYZ, 3
                       DO INEIB=1, NNEIB_EST_DUMMY+1
                          IF ((INEIB.EQ.1).AND.(INTYP0.NE.KNTYP0)) THEN
                             VEC_LOC_STRESS1(INEIB,1)=0.0_q
                          ELSE
                             IF (.NOT.FFM%LSUPERVEC) THEN
                                VEC_LOC_STRESS1(INEIB,1)=-XYZ_EST1(JXYZ,INEIB,INTYP0,LOCVAR_COL)/OMEGA
                             ELSE
                                VEC_LOC_STRESS1(INEIB,1)=-XYZ_EST_MAX(JXYZ,INEIB,INTYP0,LOCVAR_COL)/OMEGA
                             ENDIF
                          ENDIF
                       ENDDO
                       CALL DGEMV ('T',NNEIB_EST_DUMMY+1,NB_USE,1.0_q,DSOAP21(1,1,IXYZ),MROW_DSOAP1, &
                            VEC_LOC_STRESS1(1,1),1,1.0_q,STRESS_HELP(1,JXYZ,IXYZ,INTYP0),1)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
! Contributions from angular descriptors.
           IF (FFM%W2.GT.0.0_q) THEN
              IF (.NOT.FFM%LSUPERVEC) THEN
                 NNEIB_EST_DUMMY=NNEIB_EST2(INTYP0,LOCVAR_COL)
              ELSE
                 NNEIB_EST_DUMMY=NNEIB_EST_MAX(INTYP0,LOCVAR_COL)
              ENDIF
              IF (NNEIB_EST_DUMMY+1.GT.0.AND.NB_USE.GT.0) THEN
                 DO IXYZ=1, 3
! Calculate force components.
                    DO INEIB=1, NNEIB_EST_DUMMY+1
                       IF ((INEIB.EQ.1).AND.(INTYP0.NE.KNTYP0)) THEN
                          VEC_LOC_FORCE2(INEIB,1)=0.0_q
                       ELSE
                          VEC_LOC_FORCE2(INEIB,1)=1.0_q
                       ENDIF
                    ENDDO
                    CALL DGEMV ('T',NNEIB_EST_DUMMY+1,NB_USE,1.0_q,DSOAP22(1,1,IXYZ),MROW_DSOAP2, &
                         VEC_LOC_FORCE2(1,1),1,1.0_q,FORCE_HELP1(1,IXYZ,LOCVAR_COL,INTYP0),1)
! Calculate stress components.
                    DO JXYZ=IXYZ, 3
                       DO INEIB=1, NNEIB_EST_DUMMY+1
                          IF ((INEIB.EQ.1).AND.(INTYP0.NE.KNTYP0)) THEN
                             VEC_LOC_STRESS2(INEIB,1)=0.0_q
                          ELSE
                             IF (.NOT.FFM%LSUPERVEC) THEN
                                VEC_LOC_STRESS2(INEIB,1)=-XYZ_EST2(JXYZ,INEIB,INTYP0,LOCVAR_COL)/OMEGA
                             ELSE
                                VEC_LOC_STRESS2(INEIB,1)=-XYZ_EST_MAX(JXYZ,INEIB,INTYP0,LOCVAR_COL)/OMEGA
                             ENDIF
                          ENDIF
                       ENDDO
                       CALL DGEMV ('T',NNEIB_EST_DUMMY+1,NB_USE,1.0_q,DSOAP22(1,1,IXYZ),MROW_DSOAP2, &
                            VEC_LOC_STRESS2(1,1),1,1.0_q,STRESS_HELP(1,JXYZ,IXYZ,INTYP0),1)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
           PROFILING_STOP('compute_force_stress')
        END SUBROUTINE COMPUTE_FORCE_STRESS


! Subroutine to copy the calculated force and stress components in the design-matrix to global arrays.

        SUBROUTINE COPY_FORCE_STRESS(FFM, PAR_SUP_HANDLE, ICONF, INTYP0, &
                                                      FORCE_HELP1, FORCE_HELP2, &
                                                      NIONS, STRESS_HELP, CALL_TYPE)
           USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (FFM_PAR)        :: FFM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           REAL(q), INTENT(INOUT) :: FORCE_HELP1(:,:,:,:)
           REAL(q), INTENT(INOUT) :: FORCE_HELP2(:)
           REAL(q), INTENT(INOUT) :: STRESS_HELP(:,:,:,:)
           INTEGER, INTENT(IN)   :: ICONF
           INTEGER, INTENT(IN)   :: INTYP0
           INTEGER, INTENT(IN)   :: NIONS
           CHARACTER(LEN=4), INTENT(IN) :: CALL_TYPE
! Local variables
           INTEGER               :: IB
           INTEGER               :: IERR
           INTEGER               :: IPROCESS
           INTEGER               :: IXYZ
           INTEGER               :: JXYZ
           INTEGER               :: KNIONS
           INTEGER               :: LOCVAR_COL
           INTEGER               :: LOCVAR_ROW
           INTEGER               :: TEST_COL
           INTEGER               :: TEST_ROW
           INTEGER               :: NB_USE
           PROFILING_START('copy_force_stress')
           ! Determine size of used local basis configs
           IF (CALL_TYPE .EQ. "PUTA") THEN
              NB_USE=FFM%NB_NEW(INTYP0)
           ELSE
              NB_USE=FFM%NB(INTYP0)
           ENDIF
           DO KNIONS=1, NIONS
              IPROCESS=MOD(KNIONS-1,PAR_SUP_HANDLE%COMM_WORLD%NCPU)
              DO IXYZ=1, 3
! Transfer the force elements to all processes.
                 IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.IPROCESS) THEN
                    LOCVAR_COL=CEILING(REAL(KNIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                    FORCE_HELP2(1:NB_USE)=FORCE_HELP1(1:NB_USE,IXYZ,LOCVAR_COL,INTYP0)
                 ENDIF
                  CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,FORCE_HELP2(1:),NB_USE,IPROCESS)
! Store the design matrix elements
#ifdef scaLAPACK
                 TEST_COL=MOD((FFM%LTIFOR(IXYZ,KNIONS,ICONF)-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.TEST_COL) THEN
                    LOCVAR_COL=GLOBAL2LOCAL(FFM%LTIFOR(IXYZ,KNIONS,ICONF),PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
                    LOCVAR_COL=FFM%LTIFOR(IXYZ,KNIONS,ICONF)
#endif
                    DO IB=1, NB_USE
#ifdef scaLAPACK
                       IF (CALL_TYPE .EQ. "PUTA") THEN
                          TEST_ROW=MOD((IB+FFM%NB(INTYP0)-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                       ELSE
                          TEST_ROW=MOD((IB-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                       ENDIF
                       IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                          CYCLE
                       ENDIF
                       IF (CALL_TYPE .EQ. "PUTA") THEN
                          LOCVAR_ROW=GLOBAL2LOCAL(IB+FFM%NB(INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
                       ELSE
                          LOCVAR_ROW=GLOBAL2LOCAL(IB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
                       ENDIF
#else
                       IF (CALL_TYPE .EQ. "PUTA") THEN
                          LOCVAR_ROW=IB+FFM%NB(INTYP0)
                       ELSE
                          LOCVAR_ROW=IB
                       ENDIF
#endif
                       FFM%FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)=FORCE_HELP2(IB)
                    ENDDO
#ifdef scaLAPACK
                 ENDIF
#endif
              ENDDO
           ENDDO
! Transfer and copy the calculated stress components in design matrix into the distributed FFM%FMAT array.
           DO IXYZ=1, 3
              DO JXYZ=IXYZ, 3
! All-reduce STRESS_HELP.
                 CALL MPI_ALLREDUCE(MPI_IN_PLACE,STRESS_HELP(1,JXYZ,IXYZ,INTYP0),NB_USE,MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
! Store the design matrix elements
#ifdef scaLAPACK
                 TEST_COL=MOD((FFM%LTSIF(JXYZ,IXYZ,ICONF)-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.TEST_COL) THEN
                    LOCVAR_COL=GLOBAL2LOCAL(FFM%LTSIF(JXYZ,IXYZ,ICONF),PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
                    LOCVAR_COL=FFM%LTSIF(JXYZ,IXYZ,ICONF)
#endif
                    DO IB=1, NB_USE
#ifdef scaLAPACK
                       IF (CALL_TYPE .EQ. "PUTA") THEN
                          TEST_ROW=MOD((IB+FFM%NB(INTYP0)-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                       ELSE
                          TEST_ROW=MOD((IB-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                       ENDIF
                       IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                          CYCLE
                       ENDIF
                       IF (CALL_TYPE .EQ. "PUTA") THEN
                          LOCVAR_ROW=GLOBAL2LOCAL(IB+FFM%NB(INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
                       ELSE
                          LOCVAR_ROW=GLOBAL2LOCAL(IB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
                       ENDIF
#else
                       IF (CALL_TYPE .EQ. "PUTA") THEN
                          LOCVAR_ROW=IB+FFM%NB(INTYP0)
                       ELSE
                          LOCVAR_ROW=IB
                       ENDIF 
#endif
                       FFM%FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)=STRESS_HELP(IB,JXYZ,IXYZ,INTYP0)
                    ENDDO
#ifdef scaLAPACK
                 ENDIF
#endif
              ENDDO
           ENDDO
           PROFILING_STOP('copy_force_stress')
        END SUBROUTINE COPY_FORCE_STRESS


        SUBROUTINE SET_DESCRIPTOR_COMBOS(FFM,FF)
           USE ML_FF_STRUCT, ONLY: FFM_PAR, FF_PAR, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (FFM_PAR)        :: FFM
           TYPE (FF_PAR)         :: FF
           INTEGER               :: INTYP0, JNTYP0, JJNTYP0
           INTEGER               :: IRB, JRB, L
           ! Initialize everything first
           FFM%NNVAR2=0
           FFM%LVAR=0
           FFM%LFLAG_VAR=.FALSE.
           FFM%WVAR=0.0_q
           IF(FFM%LSIC) THEN
              FFM%NNVAR2_SIC=0
              FFM%LVAR_SIC=0
              FFM%LFLAG_VAR_SIC=.FALSE.
           ENDIF
           IF (FFM%DESC_TYPE.EQ.2 .OR. FFM%DESC_TYPE.EQ.4 .OR. FFM%DESC_TYPE.EQ.7 .OR. FFM%DESC_TYPE.EQ.9) THEN
              FFM%RV%NNVAR2_RED=0
              FFM%RV%LVAR_RED=0
              FFM%RV%LFLAG_VAR_RED=.FALSE.
              FFM%RV%WVAR_RED=0.0_q
           ENDIF
           ! Now calculate combos
           IF (FFM%DESC_TYPE.EQ.0) THEN
              DO INTYP0=1, FF%MTYP
                 DO JNTYP0=1, FF%MTYP
                    DO JJNTYP0=1, FF%MTYP
                       DO L=0, FFM%LMAX2
                          DO IRB=1, FFM%NRB2(L)
                             DO JRB=IRB, FFM%NRB2(L)
                                FFM%NNVAR2(INTYP0)=FFM%NNVAR2(INTYP0)+1
                                FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=FFM%NNVAR2(INTYP0)
                                FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                IF(JRB.EQ.IRB) THEN
                                   FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                ELSE
                                   FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                ENDIF
                                IF(FFM%LSIC) THEN
                                   IF(JJNTYP0.EQ.JNTYP0) THEN
                                      FFM%NNVAR2_SIC(INTYP0)=FFM%NNVAR2_SIC(INTYP0)+1
                                      FFM%LVAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=FFM%NNVAR2_SIC(INTYP0)
                                      FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=.TRUE.
                                   ENDIF
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
              DO INTYP0=1, FF%MTYP
                 DO JNTYP0=1, FF%MTYP
                    DO L=0, FFM%LMAX2
                       DO IRB=1, FFM%NRB2(L)
                          DO JRB=IRB, FFM%NRB2(L)
                             FFM%NNVAR2(INTYP0)=FFM%NNVAR2(INTYP0)+1
                             FFM%LVAR(JRB,IRB,L,1,JNTYP0,INTYP0)=FFM%NNVAR2(INTYP0)
                             FFM%LFLAG_VAR(JRB,IRB,L,1,JNTYP0,INTYP0)=.TRUE.
                             IF(JRB.EQ.IRB) THEN
                                FFM%WVAR(JRB,IRB,L,1,JNTYP0,INTYP0)=1.0_q
                             ELSE
                                FFM%WVAR(JRB,IRB,L,1,JNTYP0,INTYP0)=SQRT(2.0_q)
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ELSE IF (FFM%DESC_TYPE.EQ.2) THEN
              ! First make the high resolution type-independent descriptor
              DO INTYP0=1, FF%MTYP
                 DO L=0, FFM%LMAX2
                    DO IRB=1, FFM%NRB2(L)
                       DO JRB=IRB, FFM%NRB2(L)
                          FFM%NNVAR2(INTYP0)=FFM%NNVAR2(INTYP0)+1
                          FFM%LVAR(JRB,IRB,L,1,1,INTYP0)=FFM%NNVAR2(INTYP0)
                          FFM%LFLAG_VAR(JRB,IRB,L,1,1,INTYP0)=.TRUE.
                          IF(JRB.EQ.IRB) THEN
                             FFM%WVAR(JRB,IRB,L,1,1,INTYP0)=1.0_q
                          ELSE
                             FFM%WVAR(JRB,IRB,L,1,1,INTYP0)=SQRT(2.0_q)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              ! Second descriptor
              DO INTYP0=1, FF%MTYP
                 DO JNTYP0=1, FF%MTYP
                    DO JJNTYP0=1, FF%MTYP
                       DO L=0, FFM%RV%LMAX2_RED
                          DO IRB=1, FFM%RV%NRB2_RED(L)
                             DO JRB=IRB, FFM%RV%NRB2_RED(L)
                                FFM%RV%NNVAR2_RED(INTYP0)=FFM%RV%NNVAR2_RED(INTYP0)+1
                                FFM%RV%LVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=FFM%RV%NNVAR2_RED(INTYP0)
                                FFM%RV%LFLAG_VAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                IF(JRB.EQ.IRB) THEN
                                   FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                ELSE
                                   FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ELSE IF (FFM%DESC_TYPE.EQ.3) THEN
              DO INTYP0=1, FF%MTYP
                 DO JNTYP0=1, FF%MTYP
                    DO JJNTYP0=1, FF%MTYP
                       DO L=0, FFM%LMAX2
                          DO IRB=1, FFM%NRB2(L)
                             FFM%NNVAR2(INTYP0)=FFM%NNVAR2(INTYP0)+1
                             FFM%LVAR(1,IRB,L,JJNTYP0,JNTYP0,INTYP0)=FFM%NNVAR2(INTYP0)
                             FFM%LFLAG_VAR(1,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                             DO JRB=IRB, FFM%NRB2(L)
                                IF(JRB.EQ.IRB) THEN
                                   FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                ELSE
                                   FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ELSE IF (FFM%DESC_TYPE.EQ.4) THEN
              ! First descriptor
              DO INTYP0=1, FF%MTYP
                 DO JNTYP0=1, FF%MTYP
                    DO JJNTYP0=1, FF%MTYP
                       DO L=0, FFM%LMAX2
                          DO IRB=1, FFM%NRB2(L)
                             DO JRB=IRB, FFM%NRB2(L)
                                FFM%NNVAR2(INTYP0)=FFM%NNVAR2(INTYP0)+1
                                FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=FFM%NNVAR2(INTYP0)
                                FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                IF(JRB.EQ.IRB) THEN
                                   FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                ELSE
                                   FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              ! Second descriptor
              DO INTYP0=1, FF%MTYP
                 DO JNTYP0=1, FF%MTYP
                    DO JJNTYP0=1, FF%MTYP
                       DO L=0, FFM%RV%LMAX2_RED
                          DO IRB=1, FFM%RV%NRB2_RED(L)
                             DO JRB=IRB, FFM%RV%NRB2_RED(L)
                                FFM%RV%NNVAR2_RED(INTYP0)=FFM%RV%NNVAR2_RED(INTYP0)+1
                                FFM%RV%LVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=FFM%RV%NNVAR2_RED(INTYP0)
                                FFM%RV%LFLAG_VAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                IF(JRB.EQ.IRB) THEN
                                   FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                ELSE
                                   FFM%RV%WVAR_RED(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        END SUBROUTINE SET_DESCRIPTOR_COMBOS

      END MODULE FORCE_FIELD_ADDITIONAL2
#endif
