#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
      MODULE FORCE_FIELD_ADDITIONAL
        USE ML_FF_CONSTANT
        USE IOHANDLE
        USE LOGFILE
        USE LAPACK_DATA
        USE MEMORY_USAGE
        USE MPI_DATA
        USE ML_FF_PREC
        USE TIME_CHECK
#ifdef use_shmem
        USE ML_FF_MPI_SHMEM
#endif
        IMPLICIT NONE

        CONTAINS

!------------------------------------------------------------------------------------------------------------------
! Calculate design matrix and spilling factor on new configurations.
! This subroutine calculates design-matrix elements on the new configuration
! using radial and angular descriptors, FFM%C00_ALL and FFM%PS_ALL, owned by all processes.
!------------------------------------------------------------------------------------------------------------------
        SUBROUTINE FMAT_NEW_MB3(ASA,AB,ABN,PAR_SUP_HANDLE,MLSHM,FF,FFM,LHEAD,LNIONS_EST1,LNIONS_EST2, &
                   LNIONS_EST_MAX,MCOL_EATOM, &
                   MNEIB_EST1,MNEIB_EST2,MNEIB_EST_MAX, &
                   MROW_C00,MROW_EATOM, &
                   MROW_EST1,MROW_EST2,MROW_EST_MAX, &
                   MROW_PS,MROW_SF,NB_MAX,NNEIB_EST1,NNEIB_EST2,NNEIB_EST_MAX,NORMALIZATION_XYZ1,NORMALIZATION_XYZ2, &
                   XYZ_EST1,XYZ_EST2,XYZ_EST1_NORM,XYZ_EST2_NORM,XYZ_EST_MAX,EATOM_HELP, &
                   ENERGY_COUPLE_HELP,ENERGY_HELP,HEAT_HELP1,HEAT_HELP2, &
                   RMEM_HELP1,RMEM_HELP2,SF,TIFOR_MB,TSIF_MB)
        USE SOAP_KERNEL
        USE ML_FF_STRUCT, ONLY: ASA_STUFF, ABINITIO, ABINITIONEW, PARALLEL_SUPER, ML_SHMEM, FF_PAR, FFM_PAR, ML_IO_WRITE
          IMPLICIT NONE
! Input and output variables
          TYPE (ASA_STUFF)      :: ASA
          TYPE (ABINITIO)       :: AB
          TYPE (ABINITIONEW)    :: ABN
          TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
          TYPE (ML_SHMEM)       :: MLSHM
          TYPE (FF_PAR)         :: FF
          TYPE (FFM_PAR)        :: FFM
          INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:) !(1:MROW_EST1,1:MTYP,1:MCOL_EST1)
          INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:) !(1:MROW_EST2,1:MTYP,1:MCOL_EST2)
          INTEGER, INTENT(IN)   :: LNIONS_EST_MAX(:,:,:) !(1:MROW_EST_MAX,1:MTYP,1:MCOL_EST_MAX)
          INTEGER, INTENT(IN)   :: MCOL_EATOM
          INTEGER, INTENT(IN)   :: MNEIB_EST1
          INTEGER, INTENT(IN)   :: MNEIB_EST2
          INTEGER, INTENT(IN)   :: MNEIB_EST_MAX
          INTEGER, INTENT(IN)   :: MROW_C00
          INTEGER, INTENT(IN)   :: MROW_EATOM
          INTEGER, INTENT(IN)   :: MROW_EST1
          INTEGER, INTENT(IN)   :: MROW_EST2
          INTEGER, INTENT(IN)   :: MROW_EST_MAX
          INTEGER, INTENT(IN)   :: MROW_PS
          INTEGER, INTENT(IN)   :: MROW_SF
          INTEGER, INTENT(IN)   :: NB_MAX  
          INTEGER, INTENT(IN)   :: NNEIB_EST1(:,:) !(1:MTYP,1:MCOL_EST1)
          INTEGER, INTENT(IN)   :: NNEIB_EST2(:,:) !(1:MTYP,1:MCOL_EST2)
          INTEGER, INTENT(IN)   :: NNEIB_EST_MAX(:,:) !(1:MTYP,1:MCOL_EST_MAX)
          REAL(q), INTENT(IN)   :: NORMALIZATION_XYZ1(:,:,:)
          REAL(q), INTENT(IN)   :: NORMALIZATION_XYZ2(:,:,:)
          REAL(q), INTENT(IN)   :: XYZ_EST1(:,:,:,:) !(1:MROW_EST1,1:3,1:MTYP,1:MCOL_EST1)
          REAL(q), INTENT(IN)   :: XYZ_EST2(:,:,:,:) !(1:MROW_EST2,1:3,1:MTYP,1:MCOL_EST2)
          REAL(q), INTENT(IN)   :: XYZ_EST1_NORM(:,:,:,:) !(1:3,1:MROW_EST1,1:MTYP,1:MCOL_EST1)
          REAL(q), INTENT(IN)   :: XYZ_EST2_NORM(:,:,:,:) !(1:3,1:MROW_EST2,1:MTYP,1:MCOL_EST2)
          REAL(q), INTENT(IN)   :: XYZ_EST_MAX(:,:,:,:) !(1:3,1:MROW_EST_MAX,1:MTYP,1:MCOL_EST_MAX)
          REAL(q), INTENT(INOUT)  :: EATOM_HELP(:,:,:) !(1:MROW_EATOM,1:MCOL_EATOM,1:FF%MTYP)
          REAL(q), INTENT(INOUT)  :: ENERGY_COUPLE_HELP(:,:,:) !(1:MROW_ENERGY_COUPLE_HELP,1:MCOL_ENERGY_COUPLE_HELP,1:MTYP)
          REAL(q), INTENT(INOUT)  :: ENERGY_HELP(:,:,:) !(1:MROW_ENERGY_HELP,1:MCOL_ENERGY_HELP,1:MTYP)
          REAL(q), INTENT(INOUT)  :: HEAT_HELP1(:,:,:,:) !(1:MROW_HEAT1,1:MCOL_HEAT1,1:FF%MTYP,1:MCOL_MIONS,1:3)
          REAL(q), INTENT(INOUT)  :: HEAT_HELP2(:,:,:,:) !(1:MROW_HEAT2,1:MCOL_HEAT2,1:FF%MTYP,1:MCOL_MIONS,1:3)
          INTEGER, INTENT(INOUT)  :: LHEAD(:)
          REAL(q), INTENT(INOUT)  :: RMEM_HELP1
          REAL(q), INTENT(INOUT)  :: RMEM_HELP2
          REAL(q), INTENT(INOUT)  :: SF(:) !(1:MROW_SF)
          REAL(q), INTENT(INOUT)  :: TIFOR_MB(:,:)
          REAL(q), INTENT(INOUT)  :: TSIF_MB(:,:)
! Local variables
          REAL(q), ALLOCATABLE  :: CLM_COUPLE_HELP1_MANYBODY(:,:,:)
          REAL(q), ALLOCATABLE  :: CLM_COUPLE_HELP2_MANYBODY(:,:,:)
          REAL(q), ALLOCATABLE  :: C00(:,:,:)
          REAL(q), ALLOCATABLE  :: C00_COUPLE(:,:,:)                   
          REAL(q), ALLOCATABLE  :: PS_COUPLE(:,:,:)                    
          REAL(q), ALLOCATABLE  :: PRODUCT_L_DPDC1(:,:,:)
          REAL(q), ALLOCATABLE  :: PRODUCT_L_DPDC2(:,:,:)
          REAL(q), ALLOCATABLE  :: L_MAT1(:,:,:)
          REAL(q), ALLOCATABLE  :: L_MAT2(:,:,:)
          REAL(q), ALLOCATABLE  :: PS(:,:,:)
          REAL(q), ALLOCATABLE  :: SOAP_ALL1(:,:,:)
          REAL(q), ALLOCATABLE  :: SOAP_ALL2(:,:,:)
          REAL(q), ALLOCATABLE  :: SOAP1(:,:,:)                          
          REAL(q), ALLOCATABLE  :: SOAP1_COUPLE(:,:,:)                   
          REAL(q), ALLOCATABLE  :: SOAP2(:,:,:)                          
          REAL(q), ALLOCATABLE  :: SOAP2_COUPLE(:,:,:)                   
          REAL(q), ALLOCATABLE  :: SOAP3(:,:)
          REAL(q), ALLOCATABLE  :: RNORM_ALL(:,:)
          REAL(q), ALLOCATABLE  :: RNORM_C00(:,:)
          REAL(q), ALLOCATABLE  :: RNORM_PS(:,:)
          INTEGER               :: MCOL_DCLM1
          INTEGER               :: MCOL_DCLM2
          INTEGER               :: MCOL_DSOAP1
          INTEGER               :: MCOL_DSOAP2
          INTEGER               :: MCOL_RNORM_ALL
          INTEGER               :: MCOL_RNORM_C00
          INTEGER               :: MCOL_RNORM_PS
          INTEGER               :: MCOL_SOAP
          INTEGER               :: MROW_CLM1  
          INTEGER               :: MROW_CLM2
          INTEGER               :: MROW_DCLM1
          INTEGER               :: MROW_DCLM2
          INTEGER               :: MROW_DSOAP1
          INTEGER               :: MROW_DSOAP2
          INTEGER               :: MROW_RNORM_ALL
          INTEGER               :: MROW_RNORM_C00
          INTEGER               :: MROW_RNORM_PS
          INTEGER               :: MROW_SOAP
          REAL(q)               :: OMEGA_INV
          REAL(q)               :: RMEM_HELP3
          REAL(q),ALLOCATABLE   :: SOAP_C001(:,:,:)
          REAL(q),ALLOCATABLE   :: SOAP_C002(:,:,:)
          REAL(q),ALLOCATABLE   :: SOAP_PS1(:,:,:)
          REAL(q),ALLOCATABLE   :: SOAP_PS2(:,:,:)
          INTEGER               :: INIONS
          INTEGER               :: INITYP
          INTEGER               :: LOCVAR_MAP_INIONS(1:ABN%MITYP,1:ABN%NTYP)
          INTEGER               :: LOCVAR_MAP_INITYP(1:ABN%MITYP,1:ABN%NTYP)
          INTEGER               :: MITYP_LocAL
          INTEGER               :: NITYP_LOCAL(1:ABN%NTYP)
          INTEGER               :: INI_LOC, INTYP, INTYP0, INIONS_TMP
          INTEGER               :: INTYP0_FF
          INTEGER               :: IHEAD

! Start profiling.
          PROFILING_START ('fmat_new_mb')
! Make index list for parallelization
          NITYP_LOCAL=0
          MITYP_LOCAL=0
          DO INIONS_TMP=1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
             INIONS=INIONS_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
             IF(INIONS.LE.ABN%NIONS) THEN
                INTYP0=ABN%TYPE_OF_EACH_ATOM_IN_ABN(INIONS)
                NITYP_LOCAL(INTYP0)=NITYP_LOCAL(INTYP0) + 1
                MITYP_LOCAL=MAX(MITYP_LOCAL,NITYP_LOCAL(INTYP0))
             ENDIF
          ENDDO 
! Initialization of memory usage.
          RMEM_HELP2=0.0_q
! Set some helping variables
          ! Define inverse volume
          OMEGA_INV=-1.0_q/ABN%OMEGA(ABN%NCONF)
! Make list for locations of beginings of basis sets in arrays.
! In this subroutine, this must be the same as FFM%LB(1,INTYP0)
          LHEAD=1
          DO INTYP0_FF=2, FF%MTYP, 1
             LHEAD(INTYP0_FF)=LHEAD(INTYP0_FF-1)+FFM%NB(INTYP0_FF-1)
          ENDDO
! Calculate array sizes for local arrays
          CALL CALCULATE_ARRAY_SIZES_FMAT_NEW_MB
! Allocate local arrays 
          CALL ALLOCATE_FMAT_NEW_MB
! Initialization of necessary arrays
          CALL INIT_VARIABLES_NEW_MB
! -------------------------------------------------------------------------
! Element loop for index i starts here
          NITYP_LOCAL=0
          DO INIONS_TMP=1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
             INIONS=INIONS_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
             IF(INIONS.LE.ABN%NIONS) THEN
                INTYP0=ABN%TYPE_OF_EACH_ATOM_IN_ABN(INIONS)
                ! Calculate element type in ff for local reference configuration type index
                INTYP0_FF=ABN%LADD(INTYP0)
                INITYP=ABN%LNITYP(INIONS)
                NITYP_LOCAL(INTYP0)=NITYP_LOCAL(INTYP0) + 1
                LOCVAR_MAP_INIONS(NITYP_LOCAL(INTYP0),INTYP0)=INIONS
                LOCVAR_MAP_INITYP(NITYP_LOCAL(INTYP0),INTYP0)=INITYP
                INI_LOC=NITYP_LOCAL(INTYP0)
                ! --------------------------------------------------------------------
                ! Calculate expansion coefficients and their derivatives
                CALL D0_D1_CLM_FMAT_NEW_MB(INIONS,INI_LOC,INTYP0)
                ! Calculate 1-point correlation (radial) descriptor.
                CALL D0C00_FMAT_NEW_MB
                ! Calculate 2-points correlation (angular) descriptor.
                CALL D0PS_FMAT_NEW_MB
                ! Normalization.
                CALL NORMALIZATION_D0_FMAT_NEW_MB
             ENDIF
          ENDDO
! Initialize the spilling factor.
          IF (FF%IERR.GT.0) THEN
             IF (MOD(FF%NSTEP,FF%IERR).EQ.0) THEN
                CALL INIT_SF_FMAT_NEW_MB
             ENDIF
          ENDIF
! -------------------------------------------------------------------------
          DO INTYP0 = 1, ABN%NTYP
             INTYP0_FF=ABN%LADD(INTYP0)
             IHEAD=LHEAD(INTYP0_FF)
             IF (NITYP_LOCAL(INTYP0).GT.0) THEN
                ! --------------------------------------------------------------------
                ! ENERGY:
                ! Calculate SOAP kernel.
                CALL SOAP_KERNEL_FMAT_NEW_MB
                ! Calculate energy components.
                CALL ENERGY_FMAT_NEW_MB
                ! Calculate spilling factor.
                IF (FF%IERR.GT.0) THEN
                   IF (MOD(FF%NSTEP,FF%IERR).EQ.0) THEN
                      CALL SF_FMAT_NEW_MB
                   ENDIF
                ENDIF
                ! --------------------------------------------------------------------
                ! FORCES and STRESS:
                ! Calculate dK/dX_i times w_i
                IF (FFM%LSUPERVEC) THEN
                   CALL CALCULATE_L_MAT_SUPERVEC
                ELSE
                   CALL CALCULATE_L_MAT
                ENDIF
                ! Calculate first tensor product L_i times dp/dC
                CALL LI_TIMES_DER_P
                ! Obtain force and stress components
                CALL CALCULATE_FORCE_STRESS_COMPONENTS_FAST
                ! --------------------------------------------------------------------
                ! Heat flux:
                IF (FFM%LHEAT) THEN
                   CALL CALCULATE_HEAT_COMPONENTS_FAST
                ENDIF
             ENDIF
          ENDDO
! Allreduce the spilling factor
          IF (FF%IERR.GT.0) THEN
             IF (MOD(FF%NSTEP,FF%IERR).EQ.0) THEN
                CALL ALLREDUCE_SF_FMAT_NEW_MB
             ENDIF
          ENDIF
! -------------------------------------------------------------------------
! All-reduce the stress components.
          CALL ALLREDUCE_STRESS_AND_FORCES_FMAT_NEW_MB
! Deallocate local arrays
          CALL DEALLOCATE_FMAT_NEW_MB
! Stop profiling.
          PROFILING_STOP ('fmat_new_mb')

          CONTAINS

! Calculate sizes of local arrays
          SUBROUTINE CALCULATE_ARRAY_SIZES_FMAT_NEW_MB
             IMPLICIT NONE
! Sizes for structural list.
             IF(.NOT.FFM%LSUPERVEC) THEN
                MROW_RNORM_ALL=1
                MCOL_RNORM_ALL=ABN%MTYP
             ELSE
                MROW_RNORM_ALL=ABN%MITYP
                MCOL_RNORM_ALL=ABN%MTYP
             ENDIF
! Size of SOAP
             MROW_SOAP=FFM%MB
             MCOL_SOAP=MITYP_LOCAL
! Radial descriptor
             IF(FFM%W1.GT.0.0_q) THEN
! Size of CLM1 and DCLM1.
                MROW_CLM1=(FFM%LMAX1+1)**2*FF%MTYP_ALLOCATE*FFM%MRB1
                MROW_DCLM1=(FFM%LMAX1+1)**2*FF%MTYP_ALLOCATE*FFM%MRB1
                MCOL_DCLM1=MNEIB_EST1+1
! Size of DSOAP1.
                IF(.NOT.FFM%LSUPERVEC) THEN
                   MROW_DSOAP1=MNEIB_EST1+1
                ELSE
                   MROW_DSOAP1=MNEIB_EST_MAX+1
                ENDIF
                MCOL_DSOAP1=FFM%MB
! Size of RNORM_C00
                IF(.NOT.FFM%LSUPERVEC) THEN
                   MROW_RNORM_C00=ABN%MITYP
                   MCOL_RNORM_C00=ABN%MTYP
                ELSE
                   MROW_RNORM_C00=1
                   MCOL_RNORM_C00=ABN%MTYP
                ENDIF
             ELSE
                MROW_CLM1=1
                MROW_DCLM1=1
                MCOL_DCLM1=1
                MROW_DSOAP1=1
                MCOL_DSOAP1=1
                MROW_RNORM_C00=1
                MCOL_RNORM_C00=ABN%MTYP
             ENDIF
! Angular descriptor
             IF(FFM%W2.GT.0.0_q) THEN
! Size of CLM2 and DCLM2.
                MROW_CLM2=(FFM%LMAX2+1)**2*FF%MTYP_ALLOCATE*FFM%MRB2
                MROW_DCLM2=(FFM%LMAX2+1)**2*FF%MTYP_ALLOCATE*FFM%MRB2
                MCOL_DCLM2=MNEIB_EST2+1
! Size of DSOAP2.
                IF(.NOT.FFM%LSUPERVEC) THEN
                   MROW_DSOAP2=(MNEIB_EST2+1)
                ELSE
                   MROW_DSOAP2=(MNEIB_EST_MAX+1)
                ENDIF
                MCOL_DSOAP2=FFM%MB
! Size of RNORM_PS
                IF(.NOT.FFM%LSUPERVEC) THEN
                   MROW_RNORM_PS=ABN%MITYP
                   MCOL_RNORM_PS=ABN%MTYP
                ELSE
                   MROW_RNORM_PS=1
                   MCOL_RNORM_PS=ABN%MTYP
                ENDIF
             ELSE
                MROW_CLM2=1
                MROW_DCLM2=1
                MCOL_DCLM2=1
                MROW_DSOAP2=1
                MCOL_DSOAP2=1
                MROW_RNORM_PS=1
                MCOL_RNORM_PS=ABN%MTYP
             ENDIF
          END SUBROUTINE CALCULATE_ARRAY_SIZES_FMAT_NEW_MB


! Allocate arrays used locally in FMAT_NEW_MB
          SUBROUTINE ALLOCATE_FMAT_NEW_MB
             IMPLICIT NONE
             LOGICAL           :: DO_ALLOCATE
             REAL(q)           :: PERCENT

             PERCENT=0.1_q

             ! Kernel related stuff
             IF (FFM%W1.GT.0.0_q) THEN
                IF (ALLOCATED(L_MAT1)) DEALLOCATE(L_MAT1)
                ALLOCATE(L_MAT1(1:FFM%MMVAR1,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(FFM%MMVAR1,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E+06_q


                IF (ALLOCATED(PRODUCT_L_DPDC1)) DEALLOCATE(PRODUCT_L_DPDC1)
                ALLOCATE(PRODUCT_L_DPDC1(1:MLSHM%MROW_CLM1_GLOBAL,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(MLSHM%MROW_CLM1_GLOBAL,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E+06_q

                IF (ALLOCATED(SOAP_C001)) DEALLOCATE(SOAP_C001)
                ALLOCATE(SOAP_C001(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP))
                IF (ALLOCATED(SOAP_C002)) DEALLOCATE(SOAP_C002)
                ALLOCATE(SOAP_C002(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+2.0_q*REAL(FFM%MB,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
             ELSE
                IF (ALLOCATED(L_MAT1)) DEALLOCATE(L_MAT1)
                ALLOCATE(L_MAT1(1:1,1:1,1:1))
                IF (ALLOCATED(PRODUCT_L_DPDC1)) DEALLOCATE(PRODUCT_L_DPDC1)
                ALLOCATE(PRODUCT_L_DPDC1(1:1,1:1,1:1))
                IF (ALLOCATED(SOAP_C001)) DEALLOCATE(SOAP_C001)
                ALLOCATE(SOAP_C001(1:1,1:1,1:1))
                IF (ALLOCATED(SOAP_C002)) DEALLOCATE(SOAP_C002)
                ALLOCATE(SOAP_C002(1:1,1:1,1:1))
             ENDIF
             IF (FFM%W2.GT.0.0_q) THEN
                IF (ALLOCATED(L_MAT2)) DEALLOCATE(L_MAT2)
                ALLOCATE(L_MAT2(1:FFM%MMVAR2,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(FFM%MMVAR2,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E+06_q

                IF (ALLOCATED(PRODUCT_L_DPDC2)) DEALLOCATE(PRODUCT_L_DPDC2)
                ALLOCATE(PRODUCT_L_DPDC2(1:MLSHM%MROW_CLM2_GLOBAL,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(MLSHM%MROW_CLM2_GLOBAL,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E+06_q

                IF (ALLOCATED(SOAP_PS1)) DEALLOCATE(SOAP_PS1)
                ALLOCATE(SOAP_PS1(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP))
                IF (ALLOCATED(SOAP_PS2)) DEALLOCATE(SOAP_PS2)
                ALLOCATE(SOAP_PS2(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+2.0_q*REAL(FFM%MB,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
             ELSE 
                IF (ALLOCATED(L_MAT2)) DEALLOCATE(L_MAT2)
                ALLOCATE(L_MAT2(1:1,1:1,1:1))
                IF (ALLOCATED(PRODUCT_L_DPDC2)) DEALLOCATE(PRODUCT_L_DPDC2)
                ALLOCATE(PRODUCT_L_DPDC2(1:1,1:1,1:1))
                IF (ALLOCATED(SOAP_PS1)) DEALLOCATE(SOAP_PS1)
                ALLOCATE(SOAP_PS1(1:1,1:1,1:1))
                IF (ALLOCATED(SOAP_PS2)) DEALLOCATE(SOAP_PS2)
                ALLOCATE(SOAP_PS2(1:1,1:1,1:1))
             ENDIF

             ! Coefficients, derivatives and descritpors
             IF (FFM%W1.GT.0.0_q) THEN
                ! Do allocation of CLM1 only in the first step for ML_ISTART=2
                IF (FF%ISTART.EQ.2) THEN
                   IF (FF%NSTEP.LE.1) THEN
                      DO_ALLOCATE=.TRUE.
                   ELSE
                      DO_ALLOCATE=.FALSE.
                   ENDIF
                ELSE
                   DO_ALLOCATE=.TRUE.
                ENDIF
                
                IF (DO_ALLOCATE) THEN
                   IF (ALLOCATED(FFM%CLM1)) DEALLOCATE(FFM%CLM1)
                   ALLOCATE (FFM%CLM1(1:MROW_CLM1,1:MITYP_LOCAL,1:ABN%NTYP)) 
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_CLM1,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
                ENDIF

                IF (ALLOCATED(C00)) DEALLOCATE(C00)
                ALLOCATE(C00(1:FFM%MMVAR1,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(FFM%MMVAR1,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q

     
                IF (FF%ISTART.EQ.2) THEN
                   IF (ALLOCATED(FFM%DCLM1)) THEN
                      IF (SIZE(FFM%DCLM1,3) .LT. MCOL_DCLM1) THEN
                         DO_ALLOCATE=.TRUE.
                      ELSE
                         DO_ALLOCATE=.FALSE.
                      ENDIF
                   ELSE
                      DO_ALLOCATE=.TRUE.
                   ENDIF
                   ! Make MCOL_DCLM1 PERCENT percent larger so that the number of reallocations is reduced
                   IF (DO_ALLOCATE) THEN
                      MCOL_DCLM1=MCOL_DCLM1 + NINT(PERCENT*MCOL_DCLM1)
                   ENDIF
                ELSE
                   DO_ALLOCATE=.TRUE.
                ENDIF
                IF (DO_ALLOCATE) THEN
                   IF (ALLOCATED(FFM%DCLM1)) DEALLOCATE(FFM%DCLM1)
                   ALLOCATE(FFM%DCLM1(1:MROW_DCLM1,1:3,1:MCOL_DCLM1,1:MITYP_LOCAL,1:ABN%NTYP))
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_DCLM1,q)*REAL(MCOL_DCLM1,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*3.0_q*8.0_q/1.0E+06_q

                   IF (ALLOCATED(FFM%DCLM1_HEAD)) DEALLOCATE(FFM%DCLM1_HEAD)
                   ALLOCATE(FFM%DCLM1_HEAD(1:MROW_DCLM1,1:3,1:FF%MTYP_ALLOCATE,1:MITYP_LOCAL,1:ABN%NTYP))
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_DCLM1,q)*REAL(FF%MTYP_ALLOCATE,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*3.0_q*8.0_q/1.0E+06_q
                ENDIF

                IF (ALLOCATED(RNORM_C00)) DEALLOCATE(RNORM_C00)
                ALLOCATE(RNORM_C00(1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q

                IF (ALLOCATED(CLM_COUPLE_HELP1_MANYBODY)) DEALLOCATE(CLM_COUPLE_HELP1_MANYBODY)
                IF (ALLOCATED(C00_COUPLE)) DEALLOCATE(C00_COUPLE)
                IF (FFM%LCOUPLE) THEN
                   ALLOCATE (CLM_COUPLE_HELP1_MANYBODY(1:MROW_CLM1,1:MITYP_LOCAL,1:ABN%NTYP)) 
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_CLM1,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
                   ALLOCATE (C00_COUPLE(1:FFM%MMVAR1,1:MITYP_LOCAL,1:ABN%NTYP)) 
                   RMEM_HELP1=RMEM_HELP1+REAL(FFM%MMVAR1,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
                ELSE
                   ALLOCATE (CLM_COUPLE_HELP1_MANYBODY(1:1,1:1,1:ABN%NTYP)) 
                   ALLOCATE (C00_COUPLE(1:1,1:1,1:ABN%NTYP)) 
                ENDIF
             ELSE 
                IF (ALLOCATED(FFM%CLM1)) DEALLOCATE(FFM%CLM1)
                ALLOCATE (FFM%CLM1(1:1,1:1,1:1)) 
                IF (ALLOCATED(C00)) DEALLOCATE(C00)
                ALLOCATE(C00(1:1,1:1,1:1))
                IF (ALLOCATED(FFM%DCLM1)) DEALLOCATE(FFM%DCLM1)
                ALLOCATE(FFM%DCLM1(1:1,1:1,1:1,1:1,1:1))
                IF (ALLOCATED(FFM%DCLM1_HEAD)) DEALLOCATE(FFM%DCLM1_HEAD)
                ALLOCATE(FFM%DCLM1_HEAD(1:1,1:1,1:1,1:1,1:1))
                IF (ALLOCATED(RNORM_C00)) DEALLOCATE(RNORM_C00)
                ALLOCATE(RNORM_C00(1:1,1:1))
                IF (ALLOCATED(CLM_COUPLE_HELP1_MANYBODY)) DEALLOCATE(CLM_COUPLE_HELP1_MANYBODY)
                ALLOCATE (CLM_COUPLE_HELP1_MANYBODY(1:1,1:1,1:ABN%NTYP)) 
                IF (ALLOCATED(C00_COUPLE)) DEALLOCATE(C00_COUPLE)
                ALLOCATE (C00_COUPLE(1:1,1:1,1:ABN%NTYP)) 
             ENDIF

             IF (FFM%W2.GT.0.0_q) THEN
                ! Do allocation of FFM%CLM2 only in the first step for ML_ISTART=2
                IF (FF%ISTART.EQ.2) THEN
                   IF (FF%NSTEP.LE.1) THEN
                      DO_ALLOCATE=.TRUE.
                   ELSE
                      DO_ALLOCATE=.FALSE.
                   ENDIF
                ELSE
                   DO_ALLOCATE=.TRUE.
                ENDIF
                IF (DO_ALLOCATE) THEN
                   IF (ALLOCATED(FFM%CLM2)) DEALLOCATE(FFM%CLM2)
                   ALLOCATE (FFM%CLM2(1:MROW_CLM2,1:MITYP_LOCAL,1:ABN%NTYP)) 
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_CLM2,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
                ENDIF

                IF (ALLOCATED(PS)) DEALLOCATE(PS)
                ALLOCATE(PS(1:FFM%MMVAR2,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(FFM%MMVAR2,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q

                IF (FF%ISTART.EQ.2) THEN
                   IF (ALLOCATED(FFM%DCLM2)) THEN
                      IF (SIZE(FFM%DCLM2,3) .LT. MCOL_DCLM2) THEN
                         DO_ALLOCATE=.TRUE.
                      ELSE
                         DO_ALLOCATE=.FALSE.
                      ENDIF
                   ELSE
                      DO_ALLOCATE=.TRUE.
                   ENDIF
                   ! Make MCOL_DCLM2 PERCENT percent larger so that the number of reallocations is reduced
                   IF (DO_ALLOCATE) THEN
                      MCOL_DCLM2=MCOL_DCLM2 + NINT(PERCENT*MCOL_DCLM2)
                   ENDIF
                ELSE
                   DO_ALLOCATE=.TRUE.
                ENDIF
                IF (DO_ALLOCATE) THEN
                   IF (ALLOCATED(FFM%DCLM2)) DEALLOCATE(FFM%DCLM2)
                   ALLOCATE(FFM%DCLM2(1:MROW_DCLM2,1:3,1:MCOL_DCLM2,1:MITYP_LOCAL,1:ABN%NTYP))
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_DCLM2,q)*REAL(MCOL_DCLM2,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*3.0_q*8.0_q/1.0E+06_q

                   IF (ALLOCATED(FFM%DCLM2_HEAD)) DEALLOCATE(FFM%DCLM2_HEAD)
                   ALLOCATE(FFM%DCLM2_HEAD(1:MROW_DCLM2,1:3,1:FF%MTYP_ALLOCATE,1:MITYP_LOCAL,1:ABN%NTYP))
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_DCLM2,q)*REAL(FF%MTYP_ALLOCATE,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*3.0_q*8.0_q/1.0E+06_q
                ENDIF

                IF (ALLOCATED(RNORM_PS)) DEALLOCATE(RNORM_PS)
                ALLOCATE(RNORM_PS(1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q

                IF (ALLOCATED(CLM_COUPLE_HELP2_MANYBODY)) DEALLOCATE(CLM_COUPLE_HELP2_MANYBODY)
                IF (ALLOCATED(PS_COUPLE)) DEALLOCATE(PS_COUPLE)
                IF (FFM%LCOUPLE) THEN
                   ALLOCATE (CLM_COUPLE_HELP2_MANYBODY(1:MROW_CLM2,1:MITYP_LOCAL,1:ABN%NTYP)) 
                   RMEM_HELP1=RMEM_HELP1+REAL(MROW_CLM2,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
                   ALLOCATE(PS_COUPLE(1:FFM%MMVAR2,1:MITYP_LOCAL,1:ABN%NTYP))
                   RMEM_HELP1=RMEM_HELP1+REAL(FFM%MMVAR2,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*4.0_q*8.0_q/1.0E6_q
                ELSE
                   ALLOCATE (CLM_COUPLE_HELP2_MANYBODY(1:1,1:1,1:ABN%NTYP)) 
                   ALLOCATE (PS_COUPLE(1:1,1:1,1:ABN%NTYP))
                ENDIF
             ELSE
                IF (ALLOCATED(FFM%CLM2)) DEALLOCATE(FFM%CLM2)
                ALLOCATE (FFM%CLM2(1:1,1:1,1:1)) 
                IF (ALLOCATED(PS)) DEALLOCATE(PS)
                ALLOCATE(PS(1:1,1:1,1:1))
                IF (ALLOCATED(FFM%DCLM2)) DEALLOCATE(FFM%DCLM2)
                ALLOCATE(FFM%DCLM2(1:1,1:1,1:1,1:1,1:1))
                IF (ALLOCATED(FFM%DCLM2_HEAD)) DEALLOCATE(FFM%DCLM2_HEAD)
                ALLOCATE(FFM%DCLM2_HEAD(1:1,1:1,1:1,1:1,1:1))
                IF (ALLOCATED(RNORM_PS)) DEALLOCATE(RNORM_PS)
                ALLOCATE(RNORM_PS(1:1,1:1))
                IF (ALLOCATED(CLM_COUPLE_HELP2_MANYBODY)) DEALLOCATE(CLM_COUPLE_HELP2_MANYBODY)
                ALLOCATE (CLM_COUPLE_HELP2_MANYBODY(1:1,1:1,1:ABN%NTYP)) 
                IF (ALLOCATED(PS_COUPLE)) DEALLOCATE(PS_COUPLE)
                ALLOCATE (PS_COUPLE(1:1,1:1,1:ABN%NTYP))
             ENDIF

             ! Allocate SOAP helping arrays, combine their memory requirement
             IF (ALLOCATED(SOAP1)) DEALLOCATE(SOAP1)
             ALLOCATE (SOAP1(1:MROW_SOAP,1:MITYP_LOCAL,1:ABN%NTYP))
             IF (ALLOCATED(SOAP2)) DEALLOCATE(SOAP2)
             ALLOCATE (SOAP2(1:MROW_SOAP,1:MITYP_LOCAL,1:ABN%NTYP))
             RMEM_HELP1=RMEM_HELP1+2.0_q*REAL(MROW_SOAP,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E+06_q
             IF (ALLOCATED(SOAP3)) DEALLOCATE(SOAP3)
             ALLOCATE (SOAP3(1:MROW_SOAP,1:MITYP_LOCAL))
             RMEM_HELP1=RMEM_HELP1+REAL(MROW_SOAP*MITYP_LOCAL*ABN%NTYP*2*8,q)/1.0E+06_q

             IF (ALLOCATED(SOAP1_COUPLE)) DEALLOCATE(SOAP1_COUPLE)
             IF (ALLOCATED(SOAP2_COUPLE)) DEALLOCATE(SOAP2_COUPLE)
             IF (FFM%LCOUPLE) THEN
                ALLOCATE (SOAP1_COUPLE(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP))
                ALLOCATE (SOAP2_COUPLE(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP))
                RMEM_HELP1=RMEM_HELP1+2.0_q*REAL(FFM%MB,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E+06_q
             ELSE
                ALLOCATE (SOAP1_COUPLE(1:1,1:1,1:ABN%MTYP_ALLOCATE))
                ALLOCATE (SOAP2_COUPLE(1:1,1:1,1:1))
             ENDIF

             ! Supervector stuff
             IF (ALLOCATED(SOAP_ALL1)) DEALLOCATE(SOAP_ALL1)
             ALLOCATE (SOAP_ALL1(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP)) 
             RMEM_HELP1=RMEM_HELP1+REAL(FFM%MB,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q
             IF (ALLOCATED(SOAP_ALL2)) DEALLOCATE(SOAP_ALL2)
             ALLOCATE (SOAP_ALL2(1:FFM%MB,1:MITYP_LOCAL,1:ABN%NTYP)) 
             RMEM_HELP1=RMEM_HELP1+REAL(FFM%MB,q)*REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q

             IF (ALLOCATED(RNORM_ALL)) DEALLOCATE(RNORM_ALL)
             ALLOCATE(RNORM_ALL(1:MITYP_LOCAL,1:ABN%NTYP))
             RMEM_HELP1=RMEM_HELP1+REAL(MITYP_LOCAL,q)*REAL(ABN%NTYP,q)*8.0_q/1.0E6_q

          END SUBROUTINE ALLOCATE_FMAT_NEW_MB

! Deallocate neccessary arrays for FMAT_NEW_MB
          SUBROUTINE DEALLOCATE_FMAT_NEW_MB
             IMPLICIT NONE
             LOGICAL           :: DO_DEALLOCATE
             !begin: new arrays instead of shared memory arrays
             IF (ALLOCATED(L_MAT1)) DEALLOCATE(L_MAT1)
             IF (ALLOCATED(PRODUCT_L_DPDC1)) DEALLOCATE(PRODUCT_L_DPDC1)
             IF (ALLOCATED(SOAP_C001)) DEALLOCATE(SOAP_C001)
             IF (ALLOCATED(SOAP_C002)) DEALLOCATE(SOAP_C002)
             IF (FF%ISTART.EQ.2) THEN
                DO_DEALLOCATE = .FALSE.
             ELSE
                DO_DEALLOCATE = .TRUE.
             ENDIF
             IF (DO_DEALLOCATE) THEN
                IF (ALLOCATED(FFM%CLM1)) DEALLOCATE(FFM%CLM1)
                IF (ALLOCATED(FFM%CLM2)) DEALLOCATE(FFM%CLM2)
                IF (ALLOCATED(FFM%DCLM1)) DEALLOCATE(FFM%DCLM1)
                IF (ALLOCATED(FFM%DCLM1_HEAD)) DEALLOCATE(FFM%DCLM1_HEAD)
                IF (ALLOCATED(FFM%DCLM2)) DEALLOCATE(FFM%DCLM2)
                IF (ALLOCATED(FFM%DCLM2_HEAD)) DEALLOCATE(FFM%DCLM2_HEAD)
             ENDIF
             IF (ALLOCATED(CLM_COUPLE_HELP1_MANYBODY)) DEALLOCATE(CLM_COUPLE_HELP1_MANYBODY)
             IF (ALLOCATED(C00)) DEALLOCATE(C00)
             IF (ALLOCATED(C00_COUPLE)) DEALLOCATE(C00_COUPLE)
             IF (ALLOCATED(RNORM_C00)) DEALLOCATE(RNORM_C00)
             IF (ALLOCATED(L_MAT2)) DEALLOCATE(L_MAT2)
             IF (ALLOCATED(PRODUCT_L_DPDC2)) DEALLOCATE(PRODUCT_L_DPDC2)
             IF (ALLOCATED(SOAP_PS1)) DEALLOCATE(SOAP_PS1)
             IF (ALLOCATED(SOAP_PS2)) DEALLOCATE(SOAP_PS2)
             IF (ALLOCATED(PS)) DEALLOCATE(PS)
             IF (ALLOCATED(RNORM_PS)) DEALLOCATE(RNORM_PS)
             IF (ALLOCATED(CLM_COUPLE_HELP2_MANYBODY)) DEALLOCATE(CLM_COUPLE_HELP2_MANYBODY)
             IF (ALLOCATED(PS_COUPLE)) DEALLOCATE(PS_COUPLE)
             IF (ALLOCATED(PS_COUPLE)) DEALLOCATE(PS_COUPLE)
             IF (ALLOCATED(SOAP_ALL1)) DEALLOCATE(SOAP_ALL1)
             IF (ALLOCATED(SOAP_ALL2)) DEALLOCATE(SOAP_ALL2)
             IF (ALLOCATED(RNORM_ALL)) DEALLOCATE(RNORM_ALL)
             IF (ALLOCATED(SOAP1)) DEALLOCATE(SOAP1)
             IF (ALLOCATED(SOAP2)) DEALLOCATE(SOAP2)
             IF (ALLOCATED(SOAP1_COUPLE)) DEALLOCATE(SOAP1_COUPLE)
             IF (ALLOCATED(SOAP2_COUPLE)) DEALLOCATE(SOAP2_COUPLE)
             IF (ALLOCATED(SOAP3)) DEALLOCATE(SOAP3)
          END SUBROUTINE DEALLOCATE_FMAT_NEW_MB

! Subroutine to initialize helping design-matrix arrays.
          SUBROUTINE INIT_VARIABLES_NEW_MB
             IMPLICIT NONE
             ! Initialization of energy components.
             ENERGY_HELP=0.0_q
             ! Initialize force stuff
             TIFOR_MB=0.0_q
             ! Initialize stress stuff
             TSIF_MB=0.0_q
             ! Initialization of heat components, if necessary.
             IF(FFM%LHEAT) THEN
                IF(FFM%W1.GT.0.0_q) THEN
                   HEAT_HELP1=0.0_q
                ENDIF
                IF(FFM%W2.GT.0.0_q) THEN
                   HEAT_HELP2=0.0_q
                ENDIF
             ENDIF
             ! Initialization of thermodynamic coupling stuff
             IF (FFM%LCOUPLE) ENERGY_COUPLE_HELP=0.0_q
             IF (FFM%W1.GT.0.0_q) THEN
                FFM%CLM1=0.0_q
                FFM%DCLM1_HEAD(1:MROW_DCLM1,1:3,1:ABN%NTYP,1:MITYP_LOCAL,1:ABN%NTYP)=0.0_q
                IF (FFM%LCOUPLE) CLM_COUPLE_HELP1_MANYBODY=0.0_q
             ENDIF
             IF (FFM%W2.GT.0.0_q) THEN
                FFM%CLM2=0.0_q
                FFM%DCLM2_HEAD(1:MROW_DCLM2,1:3,1:ABN%NTYP,1:MITYP_LOCAL,1:ABN%NTYP)=0.0_q
                IF (FFM%LCOUPLE) CLM_COUPLE_HELP2_MANYBODY=0.0_q
             ENDIF
             ! Initiliaze potential energy calculation
             IF (FFM%LEATOM.OR.FFM%LHEAT.OR.FFM%LCOUPLE) THEN
                ABN%EATOM_POT=0.0_q
             ENDIF
             ! Initialize SOAP arrays
             IF (.NOT.FFM%LSUPERVEC) THEN
                IF (FFM%W1.GT.0.0_q) THEN
                   SOAP_C001=0.0_q
                   SOAP_C002=0.0_q
                ENDIF
                IF (FFM%W2.GT.0.0_q) THEN
                   SOAP_PS1=0.0_q
                   SOAP_PS2=0.0_q
                ENDIF
             ELSE
                SOAP_ALL1=0.0_q
                SOAP_ALL2=0.0_q
             ENDIF
             ! Initialize helping arrays for derivatives 
             IF (FFM%W1.GT.0.0_q) THEN
                L_MAT1=0.0_q
                PRODUCT_L_DPDC1=0.0_q
             ENDIF
             IF (FFM%W2.GT.0.0_q) THEN
                L_MAT2=0.0_q
                PRODUCT_L_DPDC2=0.0_q
             ENDIF
             !Initialize norms
             IF(.NOT.FFM%LSUPERVEC) THEN
                IF (FFM%W1.GT.0.0_q) RNORM_C00=0.0_q
                IF (FFM%W2.GT.0.0_q) RNORM_PS=0.0_q
             ELSE
                RNORM_ALL=0.0_q
             ENDIF
          END SUBROUTINE INIT_VARIABLES_NEW_MB


          SUBROUTINE D0_D1_CLM_FMAT_NEW_MB(INIONS,INI_LOC,INTYP0)
             IMPLICIT NONE
             INTEGER               :: IERR
             INTEGER               :: INIONS,INI_LOC,INTYP0
             INTEGER               :: LOCVAR_COL
             INTEGER               :: LOC_HELP_COUPLE_CLM1, LOC_HELP_COUPLE_CLM2
             PROFILING_START('ml_fast_d0_d1_clm_fmat_new_mb')
             !Radial part
             LOCVAR_COL=CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
             IF (FFM%LCOUPLE) THEN
                LOC_HELP_COUPLE_CLM1=INI_LOC
                LOC_HELP_COUPLE_CLM2=INI_LOC
             ELSE
                LOC_HELP_COUPLE_CLM1=1
                LOC_HELP_COUPLE_CLM2=1
             ENDIF
             IF(FFM%W1.GT.0.0_q) THEN
!                CALL D0_D1_CLM_SOAP_TEST(FFM%DRSPL1,FFM%FNL1,FFM%RCOUPLE_IONS, &
!                     XYZ_EST1_NORM(1:,1:,1:,LOCVAR_COL),NORMALIZATION_XYZ1(1:,1:,LOCVAR_COL), &
!                     INIONS,INTYP0,INTYP0_FF,FFM%LMAX1,FFM%MRB1,ABN%MTYP,FFM%NSPL1,ABN%NTYP,ABN%LADD(1:), &
!                     LNIONS_EST1(1:,1:,LOCVAR_COL),NNEIB_EST1(1:,LOCVAR_COL),FFM%NRB1(0:),FFM%LCOUPLE,FFM%LIONS_COUPLE, &
!                     FFM%DCLM1(:,:,:,INI_LOC,INTYP0),FFM%DCLM1_HEAD(:,:,:,INI_LOC,INTYP0),FFM%CLM1(:,INI_LOC,INTYP0), &
!                     CLM_COUPLE_HELP1_MANYBODY(:,LOC_HELP_COUPLE_CLM1,INTYP0),FFM%IVAR_OFFSET1, &
!                     ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
                CALL D0_D1_CLM_SOAP_FAST(FFM%DRSPL1,FFM%FNL1,FFM%RCOUPLE_IONS, &
                     XYZ_EST1_NORM(1:,1:,1:,LOCVAR_COL),NORMALIZATION_XYZ1(1:,1:,LOCVAR_COL), &
                     INIONS,INTYP0,INTYP0_FF,FFM%LMAX1,FFM%MRB1,ABN%MTYP,FFM%NSPL1,ABN%NTYP,ABN%LADD(1:), &
                     LNIONS_EST1(1:,1:,LOCVAR_COL),NNEIB_EST1(1:,LOCVAR_COL),FFM%NRB1(0:),FFM%LCOUPLE,FFM%LIONS_COUPLE, &
                     FFM%DCLM1(:,:,:,INI_LOC,INTYP0),FFM%DCLM1_HEAD(:,:,:,INI_LOC,INTYP0),FFM%CLM1(:,INI_LOC,INTYP0), &
                     CLM_COUPLE_HELP1_MANYBODY(:,LOC_HELP_COUPLE_CLM1,INTYP0),FFM%IVAR_OFFSET1, &
                     ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
!                CALL D0_D1_CLM_SOAP_TEST(FFM%DRSPL2,FFM%FNL2,FFM%RCOUPLE_IONS, &
!                     XYZ_EST2_NORM(1:,1:,1:,LOCVAR_COL),NORMALIZATION_XYZ2(1:,1:,LOCVAR_COL), &
!                     INIONS,INTYP0,INTYP0_FF,FFM%LMAX2,FFM%MRB2,ABN%MTYP,FFM%NSPL2,ABN%NTYP,ABN%LADD(1:), &
!                     LNIONS_EST2(1:,1:,LOCVAR_COL),NNEIB_EST2(1:,LOCVAR_COL),FFM%NRB2(0:),FFM%LCOUPLE,FFM%LIONS_COUPLE, &
!                     FFM%DCLM2(:,:,:,INI_LOC,INTYP0),FFM%DCLM2_HEAD(:,:,:,INI_LOC,INTYP0),FFM%CLM2(:,INI_LOC,INTYP0), &
!                     CLM_COUPLE_HELP2_MANYBODY(:,LOC_HELP_COUPLE_CLM2,INTYP0),FFM%IVAR_OFFSET2, &
!                     ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
                CALL D0_D1_CLM_SOAP_FAST(FFM%DRSPL2,FFM%FNL2,FFM%RCOUPLE_IONS, &
                     XYZ_EST2_NORM(1:,1:,1:,LOCVAR_COL),NORMALIZATION_XYZ2(1:,1:,LOCVAR_COL), &
                     INIONS,INTYP0,INTYP0_FF,FFM%LMAX2,FFM%MRB2,ABN%MTYP,FFM%NSPL2,ABN%NTYP,ABN%LADD(1:), &
                     LNIONS_EST2(1:,1:,LOCVAR_COL),NNEIB_EST2(1:,LOCVAR_COL),FFM%NRB2(0:),FFM%LCOUPLE,FFM%LIONS_COUPLE, &
                     FFM%DCLM2(:,:,:,INI_LOC,INTYP0),FFM%DCLM2_HEAD(:,:,:,INI_LOC,INTYP0),FFM%CLM2(:,INI_LOC,INTYP0), &
                     CLM_COUPLE_HELP2_MANYBODY(:,LOC_HELP_COUPLE_CLM2,INTYP0),FFM%IVAR_OFFSET2, &
                     ASA%YLM3I,ASA%JS,ASA%INDCG,ASA%LMAXCG,ASA%YLM3LOOKUP_TABLE)
             ENDIF
             PROFILING_STOP('ml_fast_d0_d1_clm_fmat_new_mb')
          END SUBROUTINE D0_D1_CLM_FMAT_NEW_MB

! Calculate radial descriptors within the fast version method.
          SUBROUTINE D0C00_FMAT_NEW_MB
             IMPLICIT NONE
! Local variables
             INTEGER               :: LOCVAR_COL_C00_COUPLE
             INTEGER               :: LOCVAR_COL_CLM1_COUPLE
             PROFILING_START('ml_fast_d0c00_soap_fmat_new_mb')
             IF(FFM%W1.GT.0.0_q) THEN
                IF(FFM%LCOUPLE) THEN
                   LOCVAR_COL_C00_COUPLE=INI_LOC
                   LOCVAR_COL_CLM1_COUPLE=INI_LOC
                ELSE
                   LOCVAR_COL_C00_COUPLE=1
                   LOCVAR_COL_CLM1_COUPLE=1
                ENDIF
                CALL D0C00_SOAP_NEW (C00(1:,INI_LOC,INTYP0),C00_COUPLE(1:,LOCVAR_COL_C00_COUPLE,INTYP0), &
                     FFM%CLM1(:,INI_LOC,INTYP0),CLM_COUPLE_HELP1_MANYBODY(:,LOCVAR_COL_CLM1_COUPLE,INTYP0), &
                     ABN%LADD(1:),FFM%LCOUPLE,FFM%LMAX1,FFM%LNORM1, &
                     MROW_C00,FFM%MMVAR1,FFM%MRB1,ABN%MTYP,MLSHM%MROW_CLM1_GLOBAL, &
                     MLSHM%MROW_CLM_COUPLE_HELP1_MANYBODY,FFM%NRB1(0:),ABN%NTYP, &
                     FFM%IVAR_OFFSET1,FFM%OFFSET_MAP1(:,:,0:))
             ENDIF
             PROFILING_STOP('ml_fast_d0c00_soap_fmat_new_mb')
          END SUBROUTINE D0C00_FMAT_NEW_MB

! Calculate angular descriptor within the fast version method.
          SUBROUTINE D0PS_FMAT_NEW_MB
             IMPLICIT NONE
! Local variables
             INTEGER               :: LOCVAR_COL_CLM2_COUPLE
             INTEGER               :: LOCVAR_COL_PS_COUPLE
             PROFILING_START('ml_fast_d0ps_soap_fmat_new_mb')
             IF(FFM%W2.GT.0.0_q) THEN
                IF(FFM%LCOUPLE) THEN
                   LOCVAR_COL_PS_COUPLE=INI_LOC
                   LOCVAR_COL_CLM2_COUPLE=INI_LOC
                ELSE
                   LOCVAR_COL_PS_COUPLE=1
                   LOCVAR_COL_CLM2_COUPLE=1
                ENDIF
                IF (FFM%DESC_TYPE.EQ.0) THEN 
!                   CALL D0PS_SOAP_NOSIC_NEW (FFM%AFILT2,FFM%CLM2(:,INI_LOC,INTYP0), &
!                        CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE,INTYP0),FFM%IAFILT2, &
!                        INTYP0,INTYP0_FF,FFM%LAFILT2,ABN%LADD(1:),FFM%LCOUPLE,FFM%LFLAG_VAR, &
!                        FFM%L_LNRB2, &
!                        FFM%LMAX2,FFM%LNORM2,FFM%LVAR, &
!                        FFM%MLNRB2,MROW_PS,FFM%MMVAR2, &
!                        FFM%MRB2,ABN%MTYP,MLSHM%MROW_CLM2_GLOBAL,MLSHM%MROW_CLM_COUPLE_HELP2_MANYBODY, &
!                        FFM%NLNRB2,FFM%NRB2(0:),FFM%NRB_LNRB2,ABN%NTYP, &
!                        PS(1:,INI_LOC,INTYP0),PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE,INTYP0),FFM%WVAR, &
!                        FFM%IVAR_OFFSET2,FFM%OFFSET_MAP2)
                   CALL D0PS_SOAP_NOSIC_NEW (FFM%AFILT2,FFM%CLM2(:,INI_LOC,INTYP0), &
                        CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE,INTYP0),FFM%IAFILT2, &
                        INTYP0,INTYP0_FF,FFM%LAFILT2,ABN%LADD(1:),FFM%LCOUPLE,FFM%LFLAG_VAR, &
                        FFM%L_LNRB2, &
                        FFM%LMAX2,FFM%LVAR, &
                        FFM%MLNRB2, &
                        FFM%MRB2, &
                        FFM%NLNRB2,FFM%NRB2(0:),FFM%NRB_LNRB2,ABN%NTYP, &
                        PS(1:,INI_LOC,INTYP0),PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE,INTYP0),FFM%WVAR, &
                        FFM%IVAR_OFFSET2,FFM%OFFSET_MAP2)
                ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
!                   CALL D0PS_SOAP_NOSIC_LIN_ELEMENT_NEW (FFM%AFILT2,FFM%CLM2(:,INI_LOC,INTYP0), &
!                        CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE,INTYP0),FFM%IAFILT2, &
!                        INTYP0,INTYP0_FF,FFM%LAFILT2,ABN%LADD(1:),FFM%LCOUPLE,FFM%LFLAG_VAR, &
!                        FFM%L_LNRB2, &
!                        FFM%LMAX2,FFM%LNORM2,FFM%LVAR, &
!                        FFM%MLNRB2,MROW_PS,FFM%MMVAR2, &
!                        FFM%MRB2,ABN%MTYP,MLSHM%MROW_CLM2_GLOBAL,MLSHM%MROW_CLM_COUPLE_HELP2_MANYBODY, &
!                        FFM%NLNRB2,FFM%NRB2(0:),FFM%NRB_LNRB2,ABN%NTYP, &
!                        PS(1:,INI_LOC,INTYP0),PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE,INTYP0),FFM%WVAR, &
!                        FFM%IVAR_OFFSET2,FFM%OFFSET_MAP2)
                   CALL D0PS_SOAP_NOSIC_LIN_ELEMENT_FAST (FFM%AFILT2,FFM%CLM2(:,INI_LOC,INTYP0), &
                        CLM_COUPLE_HELP2_MANYBODY(:,LOCVAR_COL_CLM2_COUPLE,INTYP0),FFM%IAFILT2, &
                        INTYP0,INTYP0_FF,FFM%LAFILT2,ABN%LADD(1:),FFM%LCOUPLE,FFM%LFLAG_VAR, &
                        FFM%L_LNRB2, &
                        FFM%LMAX2,FFM%LVAR, &
                        FFM%MLNRB2, &
                        FFM%MRB2, &
                        FFM%NLNRB2,FFM%NRB2(0:),FFM%NRB_LNRB2,ABN%NTYP, &
                        PS(1:,INI_LOC,INTYP0),PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE,INTYP0),FFM%WVAR, &
                        FFM%IVAR_OFFSET2,FFM%OFFSET_MAP2)
                ENDIF
             ENDIF
             PROFILING_STOP('ml_fast_d0ps_soap_fmat_new_mb')
          END SUBROUTINE D0PS_FMAT_NEW_MB

! Normalization.
          SUBROUTINE NORMALIZATION_D0_FMAT_NEW_MB
             IMPLICIT NONE
! Local variables
             INTEGER               :: IERR
             INTEGER               :: LOCVAR_COL_C00_COUPLE
             INTEGER               :: LOCVAR_COL_PS_COUPLE
             INTEGER               :: INI_LOC1
             INTEGER               :: INI_LOC2
             INTEGER               :: INTYP0_1
             INTEGER               :: INTYP0_2
             INTEGER  :: I1
             PROFILING_START('ml_fast_normalization_d0_fmat_new_mb')
! Initialization of norm.
             IF(.NOT.FFM%LSUPERVEC) THEN
                IF(FFM%LCOUPLE) THEN
                   LOCVAR_COL_C00_COUPLE=INI_LOC
                   LOCVAR_COL_PS_COUPLE=INI_LOC
                ELSE
                   LOCVAR_COL_C00_COUPLE=1
                   LOCVAR_COL_PS_COUPLE=1
                ENDIF
             ELSE
                IF(FFM%LCOUPLE) THEN
                   LOCVAR_COL_C00_COUPLE=INI_LOC
                   LOCVAR_COL_PS_COUPLE=INI_LOC
                ELSE
                   LOCVAR_COL_C00_COUPLE=1
                   LOCVAR_COL_PS_COUPLE=1
                ENDIF
             ENDIF
             IF ( FFM%W1 .GT. 0.0_q ) THEN
                INI_LOC1 =  INI_LOC
                INTYP0_1 =  INTYP0
             ELSE
                INI_LOC1 =  1 
                INTYP0_1 =  1 
             END IF
             IF ( FFM%W2 .GT. 0.0_q ) THEN
                INI_LOC2 =  INI_LOC
                INTYP0_2 =  INTYP0
             ELSE
                INI_LOC2 =  1 
                INTYP0_2 =  1 
             END IF
             CALL NORMALIZATION_D0_FAST (C00(1:,INI_LOC1,INTYP0_1), &
                  C00_COUPLE(1:,LOCVAR_COL_C00_COUPLE,INTYP0_1),FFM%LSUPERVEC, &
                  FFM%LCOUPLE,FFM%LNORM1,FFM%LNORM2, &
                  MROW_C00,MROW_PS,FFM%NNVAR1, &
                  FFM%NNVAR2(INTYP0_FF),PS(1:,INI_LOC2,INTYP0_2), &
                  PS_COUPLE(1:,LOCVAR_COL_PS_COUPLE,INTYP0_2), &
                  RNORM_ALL(INI_LOC,INTYP0),RNORM_C00(INI_LOC1,INTYP0_1), &
                  RNORM_PS(INI_LOC2,INTYP0_2),FFM%W1,FFM%W2)
             PROFILING_STOP('ml_fast_normalization_d0_fmat_new_mb')
          END SUBROUTINE NORMALIZATION_D0_FMAT_NEW_MB

! Calculate SOAP kernel by using the fast version method.
          SUBROUTINE SOAP_KERNEL_FMAT_NEW_MB
             IMPLICIT NONE
             ! Local variables
             INTEGER               :: IB
             INTEGER               :: IERR
             INTEGER               :: ILOC
             INTEGER               :: LOCVAR_COL
             REAL(q)               :: FACTOR_HELP
             REAL(q)               :: C00_HELP, PS_HELP
             REAL(q)               :: RNORM_ALL_HELP
             REAL(q)               :: RNORM_C00_HELP
             REAL(q)               :: RNORM_PS_HELP
             REAL(q)               :: SOAP_ALL_HELP
             REAL(q)               :: EPS_TOL
             ! Set helping index for atoms
             LOCVAR_COL=1
             ! Set tolerance for minimal value of denominator of normalization
             EPS_TOL=1E-10_q
             IF (FFM%NB(INTYP0_FF).GT.0) THEN
                ! Initialize SOAP2 matrix
                DO ILOC=1, NITYP_LOCAL(INTYP0)
                   DO IB=1, FFM%NB(INTYP0_FF)
                      SOAP2(IB,ILOC,INTYP0)=0.0_q
                   ENDDO
                ENDDO
                IF (FFM%LCOUPLE) THEN
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      DO IB=1, FFM%NB(INTYP0_FF)
                         SOAP2_COUPLE(IB,ILOC,INTYP0)=0.0_q
                      ENDDO
                   ENDDO
                ENDIF
                ! Calculate contribution from 1-point correlation function (radial descriptor)
                IF(FFM%W1.GT.0.0_q) THEN
                   PROFILING_START('ml_fast_dgemm_soap1')
#ifdef scaLAPACK
                   CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR1,1.0_q, &
                               FFM%C00_ALL(1,IHEAD),FFM%MMVAR1, &
                               C00(1,1,INTYP0),MROW_C00, &
                               0.0_q,SOAP1(1,1,INTYP0),MROW_SOAP)
                   IF(FFM%LCOUPLE) THEN
                      CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR1,1.0_q, &
                                  FFM%C00_ALL(1,IHEAD),FFM%MMVAR1, &
                                  C00_COUPLE(1,1,INTYP0),FFM%MMVAR1, &
                                  0.0_q,SOAP1_COUPLE(1,1,INTYP0),MROW_SOAP)
                   ENDIF
#else
                   CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR1,1.0_q, &
                               FFM%C00(1,1,INTYP0_FF),FFM%MMVAR1, &
                               C00(1,1,INTYP0),MROW_C00, &
                               0.0_q,SOAP1(1,1,INTYP0),MROW_SOAP)
                   IF(FFM%LCOUPLE) THEN
                      CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR1,1.0_q, &
                                  FFM%C00(1,1,INTYP0_FF),FFM%MMVAR1, &
                                  C00_COUPLE(1,1,INTYP0),FFM%MMVAR1, &
                                  0.0_q,SOAP1_COUPLE(1,1,INTYP0),MROW_SOAP)
                   ENDIF
#endif
                   PROFILING_STOP('ml_fast_dgemm_soap1')

                   PROFILING_START('ml_fast_gather_soap_c00')
                   IF(.NOT.FFM%LSUPERVEC) THEN
                   ! Store the calculated SOAP1 to SOAP_C001,2.
                      FACTOR_HELP=FFM%W1*REAL(FFM%NHYP1,q)
                      DO ILOC=1, NITYP_LOCAL(INTYP0)
                         IF (RNORM_C00(ILOC,INTYP0).GT.EPS_TOL) THEN
                            RNORM_C00_HELP=FACTOR_HELP/RNORM_C00(ILOC,INTYP0)
                         ELSE
                            RNORM_C00_HELP=0.0_q
                         ENDIF
                         DO IB=1, FFM%NB(INTYP0_FF)
                            C00_HELP=SOAP1(IB,ILOC,INTYP0)**(FFM%NHYP1-1)
                            SOAP_C001(IB,ILOC,INTYP0)=C00_HELP*RNORM_C00_HELP
                            SOAP_C002(IB,ILOC,INTYP0)=SOAP_C001(IB,ILOC,INTYP0)*SOAP1(IB,ILOC,INTYP0)
                            SOAP2(IB,ILOC,INTYP0)=SOAP2(IB,ILOC,INTYP0)+C00_HELP*SOAP1(IB,ILOC,INTYP0)*FFM%W1
                         ENDDO
                      ENDDO
                      IF(FFM%LCOUPLE) THEN
                         DO ILOC=1, NITYP_LOCAL(INTYP0)
                            DO IB=1, FFM%NB(INTYP0_FF)
                               SOAP2_COUPLE(IB,ILOC,INTYP0)=SOAP2_COUPLE(IB,ILOC,INTYP0) + &
                                            FACTOR_HELP * SOAP1_COUPLE(IB,ILOC,INTYP0) * &
                                            SOAP1(IB,ILOC,INTYP0)**(FFM%NHYP1-1)
                            ENDDO
                         ENDDO
                      ENDIF
                   ELSE
                      DO ILOC=1, NITYP_LOCAL(INTYP0)
                         DO IB=1, FFM%NB(INTYP0_FF)
                            SOAP2(IB,ILOC,INTYP0)=SOAP2(IB,ILOC,INTYP0) + SOAP1(IB,ILOC,INTYP0)
                         ENDDO
                      ENDDO
                      IF (FFM%LCOUPLE) THEN
                         DO ILOC=1, NITYP_LOCAL(INTYP0)
                            DO IB=1, FFM%NB(INTYP0_FF)
                               SOAP2_COUPLE(IB,ILOC,INTYP0)= SOAP2_COUPLE(IB,ILOC,INTYP0) + &
                                                             SOAP1_COUPLE(IB,ILOC,INTYP0)
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDIF
                   PROFILING_STOP('ml_fast_gather_soap_c00')
                ENDIF
 
                ! Calculate contribution from 2-point correlation function (angular descriptore)
                IF (FFM%W2.GT.0.0_q) THEN
                   PROFILING_START('ml_fast_dgemm_soap2')
#ifdef scaLAPACK
                   CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR2(INTYP0_FF),1.0_q, &
                               FFM%PS_ALL(1,IHEAD),FFM%MMVAR2, &
                               PS(1,1,INTYP0),MROW_PS, &
                               0.0_q,SOAP1(1,1,INTYP0),MROW_SOAP)
                   ! Thermodynamic integration
                   IF (FFM%LCOUPLE) THEN
                      CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR2(INTYP0_FF),1.0_q, &
                                  FFM%PS_ALL(1,IHEAD),FFM%MMVAR2, &
                                  PS_COUPLE(1,1,INTYP0),FFM%MMVAR2, &
                                  0.0_q,SOAP1_COUPLE(1,1,INTYP0),MROW_SOAP)
                   ENDIF
#else
                   CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR2(INTYP0_FF),1.0_q, &
                               FFM%PS(1,1,INTYP0_FF),FFM%MMVAR2, &
                               PS(1,1,INTYP0),MROW_PS, &
                               0.0_q,SOAP1(1,1,INTYP0),MROW_SOAP)
                   ! Thermodynamic integration
                   IF (FFM%LCOUPLE) THEN
                      CALL DGEMM ('T','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NNVAR2(INTYP0_FF),1.0_q, &
                                  FFM%PS(1,1,INTYP0_FF),FFM%MMVAR2, &
                                  PS_COUPLE(1,1,INTYP0),FFM%MMVAR2, &
                                  0.0_q,SOAP1_COUPLE(1,1,INTYP0),MROW_SOAP)
                   ENDIF
#endif
                   PROFILING_STOP('ml_fast_dgemm_soap2')

                   ! Store the calculated SOAP1 to SOAP_PS1,2.
                   PROFILING_START('ml_fast_gather_soap_ps')
                   IF (.NOT.FFM%LSUPERVEC) THEN
                      FACTOR_HELP=FFM%W2*REAL(FFM%NHYP2,q)
                      DO ILOC=1, NITYP_LOCAL(INTYP0)
                         IF (RNORM_PS(ILOC,INTYP0).GT.EPS_TOL) THEN
                            RNORM_PS_HELP=FACTOR_HELP/RNORM_PS(ILOC,INTYP0)
                         ELSE
                            RNORM_PS_HELP=0.0_q
                         ENDIF
                         DO IB=1, FFM%NB(INTYP0_FF)
                            PS_HELP=(SOAP1(IB,ILOC,INTYP0))**(FFM%NHYP2-1)
                            SOAP_PS1(IB,ILOC,INTYP0)=PS_HELP*RNORM_PS_HELP
                            SOAP_PS2(IB,ILOC,INTYP0)=SOAP_PS1(IB,ILOC,INTYP0)*SOAP1(IB,ILOC,INTYP0)
                            SOAP2(IB,ILOC,INTYP0)=SOAP2(IB,ILOC,INTYP0)+PS_HELP*SOAP1(IB,ILOC,INTYP0)*FFM%W2
                         ENDDO
                      ENDDO
                      IF(FFM%LCOUPLE) THEN
                         DO ILOC=1, NITYP_LOCAL(INTYP0)
                            DO IB=1, FFM%NB(INTYP0_FF)
                               SOAP2_COUPLE(IB,ILOC,INTYP0)=SOAP2_COUPLE(IB,ILOC,INTYP0) + &
                                            FACTOR_HELP * SOAP1_COUPLE(IB,ILOC,INTYP0) * &
                                            SOAP1(IB,ILOC,INTYP0)**(FFM%NHYP2-1)
                            ENDDO
                         ENDDO
                      ENDIF
                   ELSE
                      DO ILOC=1, NITYP_LOCAL(INTYP0)
                         DO IB=1, FFM%NB(INTYP0_FF)
                            SOAP2(IB,ILOC,INTYP0)=SOAP2(IB,ILOC,INTYP0) + SOAP1(IB,ILOC,INTYP0)
                         ENDDO
                      ENDDO
                      IF (FFM%LCOUPLE) THEN
                         DO ILOC=1, NITYP_LOCAL(INTYP0)
                            DO IB=1, FFM%NB(INTYP0_FF)
                               SOAP2_COUPLE(IB,ILOC,INTYP0)= SOAP2_COUPLE(IB,ILOC,INTYP0) + &
                                                             SOAP1_COUPLE(IB,ILOC,INTYP0)
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDIF
                   PROFILING_STOP('ml_fast_gather_soap_ps')
                ENDIF
                ! This next part only needs to be done for the supervector
                IF (FFM%LSUPERVEC) THEN
                   PROFILING_START('ml_fast_gather_soap_all')
                   IF (FFM%LCOUPLE) THEN
                      FACTOR_HELP=REAL(FFM%NHYP2,q)
                      DO ILOC=1, NITYP_LOCAL(INTYP0)
                         IF (RNORM_ALL(ILOC,INTYP0).GT.EPS_TOL) THEN
                            RNORM_ALL_HELP=REAL(FFM%NHYP2,q)/RNORM_ALL(ILOC,INTYP0) 
                         ELSE
                            RNORM_ALL_HELP=0.0_q
                         ENDIF
                         DO IB=1, FFM%NB(INTYP0_FF)
                            SOAP_ALL_HELP=(SOAP2(IB,ILOC,INTYP0))**(FFM%NHYP2-1)
                            SOAP_ALL1(IB,ILOC,INTYP0)=SOAP_ALL_HELP*RNORM_ALL_HELP
                            SOAP2(IB,ILOC,INTYP0)=SOAP_ALL_HELP*SOAP2(IB,ILOC,INTYP0)
                            SOAP_ALL2(IB,ILOC,INTYP0)=SOAP2(IB,ILOC,INTYP0)*RNORM_ALL_HELP
                            SOAP2_COUPLE(IB,ILOC,INTYP0)=FACTOR_HELP* &
                                  SOAP2_COUPLE(IB,ILOC,INTYP0) * SOAP_ALL_HELP
                         ENDDO
                      ENDDO
                   ELSE
                      DO ILOC=1, NITYP_LOCAL(INTYP0)
                         IF (RNORM_ALL(ILOC,INTYP0).GT.EPS_TOL) THEN
                            RNORM_ALL_HELP=REAL(FFM%NHYP2,q)/RNORM_ALL(ILOC,INTYP0) 
                         ELSE
                            RNORM_ALL_HELP=0.0_q
                         ENDIF
                         DO IB=1, FFM%NB(INTYP0_FF)
                            SOAP_ALL_HELP=(SOAP2(IB,ILOC,INTYP0))**(FFM%NHYP2-1)
                            SOAP_ALL1(IB,ILOC,INTYP0)=SOAP_ALL_HELP*RNORM_ALL_HELP
                            SOAP2(IB,ILOC,INTYP0)=SOAP_ALL_HELP*SOAP2(IB,ILOC,INTYP0)
                            SOAP_ALL2(IB,ILOC,INTYP0)=SOAP2(IB,ILOC,INTYP0)*RNORM_ALL_HELP
                         ENDDO
                      ENDDO
                   ENDIF
                   PROFILING_STOP('ml_fast_gather_soap_all')
                ENDIF
                  
             ENDIF
          END SUBROUTINE SOAP_KERNEL_FMAT_NEW_MB

! Calculate energy components.
          SUBROUTINE ENERGY_FMAT_NEW_MB
            IMPLICIT NONE
            ! Local variables
            INTEGER               :: IB, ICOL, IERR, IROW
            INTEGER               :: ILOC, INIONS_LOC
            INTEGER               :: LOCVAR_COL
            REAL(q)               :: FACTOR_HELP
            ! Lapack routines
            REAL(q)               :: DDOT
! Calculate energy matrix elements
            PROFILING_START('ml_fast_energy')
            IF(FFM%NB(INTYP0).GT.0) THEN
               DO ILOC=1, NITYP_LOCAL(INTYP0)
                  INIONS_LOC = LOCVAR_MAP_INIONS(ILOC,INTYP0)
                  IF (FF%LTOTEN_SYSTEM) THEN
                     FACTOR_HELP=FFM%RCOUPLE_IONS(INIONS_LOC)
                  ELSE
                     FACTOR_HELP=FFM%RCOUPLE_IONS(INIONS_LOC)/REAL(ABN%NIONS,q)
                  ENDIF
                  CALL DAXPY(FFM%NB(INTYP0_FF),FACTOR_HELP,SOAP2(1,ILOC,INTYP0),1,ENERGY_HELP(1,1,INTYP0),1)
                  IF (FFM%LCOUPLE) THEN
                     CALL DAXPY(FFM%NB(INTYP0_FF),FACTOR_HELP,SOAP2_COUPLE(1,ILOC,INTYP0), &
                                1,ENERGY_COUPLE_HELP(1,1,INTYP0),1)
                  ENDIF
! Store the atomic energy elements.
                  IF(FFM%LEATOM.OR.FFM%LHEAT.OR.FFM%LCOUPLE) THEN
#ifdef scaLAPACK
                     ABN%EATOM_POT(INIONS_LOC)=ABN%EATOM_POT(INIONS_LOC) + &
                                               DDOT(FFM%NB(INTYP0_FF),SOAP2(1,ILOC,INTYP0), &
                                               1,FFM%WMAT_ALL(IHEAD),1)
#else
                     ABN%EATOM_POT(INIONS_LOC)=ABN%EATOM_POT(INIONS_LOC) + &
                                               DDOT(FFM%NB(INTYP0_FF),SOAP2(1,ILOC,INTYP0), &
                                               1,FFM%WMAT(1,1,INTYP0_FF),1)
#endif
                     IF (FFM%ISCALE_TOTEN.EQ.1) THEN
                        ABN%EATOM_POT(INIONS_LOC)=ABN%EATOM_POT(INIONS_LOC) + &
                            FF%EATOM_REF(INTYP0_FF)
                     ELSE IF (FFM%ISCALE_TOTEN.EQ.2) THEN
                        ABN%EATOM_POT(INIONS_LOC)=ABN%EATOM_POT(INIONS_LOC) + &
                            FF%EATOM_REF(INTYP0_FF) + FFM%TOTENAV
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            PROFILING_STOP('ml_fast_energy')
          END SUBROUTINE ENERGY_FMAT_NEW_MB

! Subroutine to initialize spilling factor.

          SUBROUTINE INIT_SF_FMAT_NEW_MB
            IMPLICIT NONE
            INTEGER :: INTYP,INTYP0
            INTEGER :: INIONS
            INTEGER :: INITYP
            DO INTYP=1, ABN%NTYP
               INTYP0=ABN%LADD(INTYP)
! If there is no basis sets, set spilling factor as unity.
               IF(FFM%NB(INTYP0).EQ.0.OR..NOT.FF%IFF) THEN
                  DO INITYP=1, ABN%NITYP(INTYP0)
                     INIONS=ABN%LNIONS(INITYP,INTYP0)
                     SF(INIONS)=1.0_q
                  ENDDO
               ELSE
! If FF%LMLONLY=.FALSE. or FF%LERR=.TRUE, then, calculate the spilling factor.
! Initialization of spilling factor. This is necessary for all processes
! because calculations are distributed to each process.
                  DO INITYP=1, ABN%NITYP(INTYP0)
                     INIONS=ABN%LNIONS(INITYP,INTYP0)
                     SF(INIONS)=0.0_q
                  ENDDO
               ENDIF
            ENDDO
          END SUBROUTINE INIT_SF_FMAT_NEW_MB

! Subroutine to calculate spilling factor.

          SUBROUTINE SF_FMAT_NEW_MB
            IMPLICIT NONE
            INTEGER :: IB, JB
            INTEGER :: IERR
            INTEGER :: ILOC
            INTEGER :: INIONS
            INTEGER :: INITYP
            INTEGER :: LOCVAR_COL
            REAL(q) :: SF_SELF
! LAPACK routine.
            REAL(q) :: DDOT
            INTEGER :: I1, I2
! Calculations of spilling factor.
            PROFILING_START ('sf_fmat_new_mb')
            IF(FFM%NB(INTYP0_FF).NE.0 .AND. FF%IFF) THEN
               CALL DGEMM ('N','N',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NB(INTYP0_FF),1.0_q,FFM%SINV_ALL(1,1,INTYP0_FF), &
                    FFM%MB,SOAP2(1,1,INTYP0),MROW_SOAP,0.0_q,SOAP3(1,1),MROW_SOAP)
               DO ILOC=1, NITYP_LOCAL(INTYP0)
! Inityp.
                  INITYP=LOCVAR_MAP_INITYP(ILOC,INTYP0)
! Inions.
                  INIONS=LOCVAR_MAP_INIONS(ILOC,INTYP0)
! Self term.
                  SF_SELF=0.0_q
                  IF(.NOT.FFM%LSUPERVEC) THEN
                     IF(FFM%W1.GT.0.0_q) THEN
                        IF(FFM%LNORM1) THEN
                           SF_SELF=SF_SELF+FFM%W1
                        ELSE
                           SF_SELF=SF_SELF+FFM%W1*RNORM_C00(INITYP,INTYP0)**FFM%NHYP1
                        ENDIF
                     ENDIF
                     IF(FFM%W2.GT.0.0_q) THEN
                        IF(FFM%LNORM2) THEN
                           SF_SELF=SF_SELF+FFM%W2
                        ELSE
                           SF_SELF=SF_SELF+FFM%W2*RNORM_PS(INITYP,INTYP0)**FFM%NHYP2
                        ENDIF
                     ENDIF
                  ELSE
! The kernel using the super vector must be always normalized!.
                     SF_SELF=SF_SELF+1.0_q
                  ENDIF
! Calculate normalized spilling factor.
                  SF(INIONS)=MIN(1.0_q,ABS(1.0_q-DDOT(FFM%NB(INTYP0_FF),SOAP2(1,ILOC,INTYP0),1,SOAP3(1,ILOC),1)/SF_SELF))
               ENDDO
            ENDIF
            PROFILING_STOP ('sf_fmat_new_mb')
          END SUBROUTINE SF_FMAT_NEW_MB

          SUBROUTINE ALLREDUCE_SF_FMAT_NEW_MB
            IMPLICIT NONE
            INTEGER :: IERR
            INTEGER :: INIONS
! All-reduce if INTYP=NTYP.
            PROFILING_START ('allred_sf_fmat_new_mb')
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,SF(1),ABN%NIONS,MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
            PROFILING_STOP ('allred_sf_fmat_new_mb')
          END SUBROUTINE ALLREDUCE_SF_FMAT_NEW_MB

          SUBROUTINE CALCULATE_L_MAT
             IMPLICIT NONE
             ! Local variables
             INTEGER                :: ILOC
             REAL(q)                :: HELPVEC(1:MITYP_LOCAL)
             ! Lapack routines
             REAL(q)                :: DDOT

             PROFILING_START('ml_fast_calculate_l_i')

             IF (FFM%NB(INTYP0_FF).GT.0) THEN
                ! Radial part
                IF (FFM%W1.GT.0.0_q) THEN
                   ! Calculate first half of equation for L_i
                   CALL DGEMM('N','N',FFM%NNVAR1,NITYP_LOCAL(INTYP0),FFM%NB(INTYP0_FF), &
                              1.0_q,FFM%C00_W_ALL(1,IHEAD),FFM%MMVAR1,SOAP_C001(1,1,INTYP0), &
                              FFM%MB_ALLOCATE,0.0_q,L_MAT1(1,1,INTYP0),FFM%MMVAR1)
                   ! Caclulate second half of equation for L_i
#ifdef scaLAPACK
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-1.0_q, &
                              SOAP_C002(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT_ALL(IHEAD),1,0.0_q,HELPVEC(1),1)
#else
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-1.0_q, &
                              SOAP_C002(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT(1,1,INTYP0_FF),1,0.0_q,HELPVEC(1),1)
#endif
                   ! Calculate elementwise product of X_i^d*HELPSCALAR_i (where d is descriptor dimension)
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      CALL DAXPY(FFM%NNVAR1,HELPVEC(ILOC),C00(1,ILOC,INTYP0),1,L_MAT1(1,ILOC,INTYP0),1)
                   ENDDO
                ENDIF
                ! Angular part
                IF (FFM%W2.GT.0.0_q) THEN
                   ! Calculate first half of equation for L_i
                   CALL DGEMM('N','N',FFM%NNVAR2(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NB(INTYP0_FF), &
                              1.0_q,FFM%PS_W_ALL(1,IHEAD),FFM%MMVAR2,SOAP_PS1(1,1,INTYP0), &
                              FFM%MB_ALLOCATE,0.0_q,L_MAT2(1,1,INTYP0),FFM%MMVAR2)
                   ! Caclulate second half of equation for L_i
#ifdef scaLAPACK
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-1.0_q, &
                              SOAP_PS2(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT_ALL(IHEAD),1,0.0_q,HELPVEC(1),1)
#else
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-1.0_q, &
                              SOAP_PS2(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT(1,1,INTYP0_FF),1,0.0_q,HELPVEC(1),1)
#endif
                   ! Calculate elementwise product of X_i^d*HELPSCALAR_i (where d is descriptor dimension)
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      CALL DAXPY(FFM%NNVAR2(INTYP0_FF),HELPVEC(ILOC),PS(1,ILOC,INTYP0),1,L_MAT2(1,ILOC,INTYP0),1)
                   ENDDO
                ENDIF
             ENDIF
             PROFILING_STOP('ml_fast_calculate_l_i')
          END SUBROUTINE CALCULATE_L_MAT



          SUBROUTINE CALCULATE_L_MAT_SUPERVEC
             IMPLICIT NONE
             ! Local variables
             INTEGER                :: ILOC
             REAL(q)                :: HELPVEC(1:MITYP_LOCAL)
             ! Lapack routines
             REAL(q)                :: DDOT

             PROFILING_START('ml_fast_calculate_l_i')

             IF (FFM%NB(INTYP0_FF).GT.0) THEN
                ! Radial part
                IF (FFM%W1.GT.0.0_q) THEN
                   ! Calculate first half of equation for L_i
                   CALL DGEMM('N','N',FFM%NNVAR1,NITYP_LOCAL(INTYP0),FFM%NB(INTYP0_FF), &
                              SQRT(FFM%W1),FFM%C00_W_ALL(1,IHEAD),FFM%MMVAR1,SOAP_ALL1(1,1,INTYP0), &
                              FFM%MB_ALLOCATE,0.0_q,L_MAT1(1,1,INTYP0),FFM%MMVAR1)
                   ! Caclulate second half of equation for L_i
#ifdef scaLAPACK
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-SQRT(FFM%W1), &
                              SOAP_ALL2(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT_ALL(IHEAD),1,0.0_q,HELPVEC(1),1)
#else
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-SQRT(FFM%W1), &
                              SOAP_ALL2(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT(1,1,INTYP0_FF),1,0.0_q,HELPVEC(1),1)
#endif
                   ! Calculate elementwise product of X_i^d*HELPSCALAR_i (where d is descriptor dimension)
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      CALL DAXPY(FFM%NNVAR1,HELPVEC(ILOC),C00(1,ILOC,INTYP0),1,L_MAT1(1,ILOC,INTYP0),1)
                   ENDDO
                ENDIF
                ! Angular part
                IF (FFM%W2.GT.0.0_q) THEN
                   ! Calculate first half of equation for L_i
                   CALL DGEMM('N','N',FFM%NNVAR2(INTYP0_FF),NITYP_LOCAL(INTYP0),FFM%NB(INTYP0_FF), &
                              SQRT(FFM%W2),FFM%PS_W_ALL(1,IHEAD),FFM%MMVAR2,SOAP_ALL1(1,1,INTYP0), &
                              FFM%MB_ALLOCATE,0.0_q,L_MAT2(1,1,INTYP0),FFM%MMVAR2)
                   ! Caclulate second half of equation for L_i
#ifdef scaLAPACK
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-SQRT(FFM%W2), &
                              SOAP_ALL2(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT_ALL(IHEAD),1,0.0_q,HELPVEC(1),1)
#else
                   CALL DGEMV('T',FFM%NB(INTYP0_FF),NITYP_LOCAL(INTYP0),-SQRT(FFM%W2), &
                              SOAP_ALL2(1,1,INTYP0),FFM%MB_ALLOCATE, &
                              FFM%WMAT(1,1,INTYP0_FF),1,0.0_q,HELPVEC(1),1)
#endif
                   ! Calculate elementwise product of X_i^d*HELPSCALAR_i (where d is descriptor dimension)
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      CALL DAXPY(FFM%NNVAR2(INTYP0_FF),HELPVEC(ILOC),PS(1,ILOC,INTYP0),1,L_MAT2(1,ILOC,INTYP0),1)
                   ENDDO
                ENDIF
             ENDIF
             PROFILING_STOP('ml_fast_calculate_l_i')
          END SUBROUTINE CALCULATE_L_MAT_SUPERVEC

!------------------------------------------------------------------------------------------------------------------
! This subroutine calculates dK/dX*dX/dp*dp/dC
! (dX/dp): this is one or zero and is not directly implemented
!------------------------------------------------------------------------------------------------------------------
          SUBROUTINE LI_TIMES_DER_P
             IMPLICIT NONE
             ! Local variables
             INTEGER               :: IVAR1, IVAR2 
             INTEGER               :: ILOC
             INTEGER               :: JJNTYP
             REAL(q)               :: DERIVATIVE_P1(1:FFM%MMVAR1,1:MLSHM%MROW_CLM1_GLOBAL)
             REAL(q)               :: DERIVATIVE_P2(1:FFM%MMVAR2,1:MLSHM%MROW_CLM2_GLOBAL)
             REAL(Q)               :: SUM_HELP
             ! Start timing
             PROFILING_START('ml_fast_fmat_l_times_der_p')

             ! Radial descriptor
             IF (FFM%W1.GT.0.0_q) THEN
                ! Calculate dp_nn'l^iJJ'/dc_nlm^i
                DO ILOC=1, NITYP_LOCAL(INTYP0)
                   DO IVAR1=1, MLSHM%MROW_CLM1_GLOBAL
                      ! Calculate dp/dC for a single nlm
                      DO IVAR2=1, FFM%NNVAR1
                         DERIVATIVE_P1(IVAR2,IVAR1)=FFM%PREFACTOR_HELP1(IVAR2,IVAR1,INTYP0)
                      ENDDO
                   ENDDO
                   CALL DGEMV('T',FFM%NNVAR1,MLSHM%MROW_CLM1_GLOBAL, &
                              1.0_q,DERIVATIVE_P1(1,1),FFM%MMVAR1,L_MAT1(1,ILOC,INTYP0), &
                              1,0.0_q,PRODUCT_L_DPDC1(1,ILOC,INTYP0),1)
                ENDDO
             ENDIF

             ! Angular descriptor
             IF (FFM%W2.GT.0.0_q) THEN
                IF (FFM%DESC_TYPE.EQ.0) THEN
                   ! Calculate dp_nn'l^iJJ'/dc_nlm^i
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      DO IVAR1=1, MLSHM%MROW_CLM2_GLOBAL
                         DO IVAR2=1, FFM%INDEX_NUM2(IVAR1,INTYP0)
                            PRODUCT_L_DPDC2(IVAR1,ILOC,INTYP0)   = &
                            PRODUCT_L_DPDC2(IVAR1,ILOC,INTYP0)   + &
                            FFM%CLM2(FFM%INDEX_PAIR_HELP2(IVAR2,IVAR1,INTYP0),ILOC,INTYP0) * &
                            FFM%PREFACTOR_HELP2(IVAR2,IVAR1,INTYP0) * &
                            L_MAT2(FFM%INDEX_MAP_TO_PS(IVAR2,IVAR1,INTYP0),ILOC,INTYP0)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
                   ! Calculate dp_nn'l^iJJ'/dc_nlm^i
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      DO IVAR1=1, MLSHM%MROW_CLM2_GLOBAL
                         ! Do first loop
                         DO IVAR2=1, FFM%INDEX_NUM2_B(IVAR1,INTYP0)
                            SUM_HELP=0.0_q
                            DO JJNTYP=1, ABN%NTYP
                               SUM_HELP=SUM_HELP + FFM%CLM2(FFM%INDEX_PAIR_HELP2_B(JJNTYP,IVAR2,IVAR1,INTYP0),ILOC,INTYP0)
                            ENDDO
                            PRODUCT_L_DPDC2(IVAR1,ILOC,INTYP0)   = &
                            PRODUCT_L_DPDC2(IVAR1,ILOC,INTYP0)   + &
                            SUM_HELP * FFM%PREFACTOR_HELP2_B(IVAR2,IVAR1,INTYP0) * &
                            L_MAT2(FFM%INDEX_MAP_TO_PS_B(IVAR2,IVAR1,INTYP0),ILOC,INTYP0)
                         ENDDO
                         ! Do second loop
                         DO IVAR2=1, FFM%INDEX_NUM2(IVAR1,INTYP0)
                            PRODUCT_L_DPDC2(IVAR1,ILOC,INTYP0)   = &
                            PRODUCT_L_DPDC2(IVAR1,ILOC,INTYP0)   + &
                            FFM%CLM2(FFM%INDEX_PAIR_HELP2(IVAR2,IVAR1,INTYP0),ILOC,INTYP0) * &
                            FFM%PREFACTOR_HELP2(IVAR2,IVAR1,INTYP0) * &
                            L_MAT2(FFM%INDEX_MAP_TO_PS(IVAR2,IVAR1,INTYP0),ILOC,INTYP0)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

             ! Stop timing
             PROFILING_STOP('ml_fast_fmat_l_times_der_p')
          END SUBROUTINE LI_TIMES_DER_P


!------------------------------------------------------------------------------------------------------------------
! This subroutine calculates directly the force and stress components via:
! (dK/dX)*(dX/dp)*(dp/dc)*(dc/dr_k)
! (dX/dp): this is one or zero and is not directly implemented
! This routine is parallelized and each core contains only some components that are non zero.
! Each core contains the whole force and stress arrays.
! These arrays are later allreduced so that each core contains the entire forces and stress.
!------------------------------------------------------------------------------------------------------------------
          SUBROUTINE CALCULATE_FORCE_STRESS_COMPONENTS_FAST
             IMPLICIT NONE
! Local variables
             INTEGER :: IXYZ, JXYZ
             INTEGER :: ILOC, INIONS_LOC
             INTEGER :: KNTYP, KNTYP0, KNEIB, KNIONS
             INTEGER :: LOCVAR_COL
             INTEGER :: IVAR
             INTEGER :: OFFSET
             INTEGER :: KNTYP0_FF
             REAL(q) :: PREFACTOR_STRESS, PREFACTOR_FORCE
             REAL(q) :: INPRODUCT_HELP
! Lapack routines
             REAL(q) :: DDOT
             PROFILING_START('ml_fast_fmat_calculate_force_stress')
             ! Radial descriptor
             IF (FFM%W1.GT.0.0_q) THEN
                DO ILOC=1, NITYP_LOCAL(INTYP0)
                   INIONS_LOC=LOCVAR_MAP_INIONS(ILOC,INTYP0)
                   LOCVAR_COL=CEILING(REAL(INIONS_LOC,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                   DO KNTYP0=1, ABN%NTYP
                      KNTYP0_FF=ABN%LADD(KNTYP0)
                      OFFSET=FFM%IVAR_OFFSET1*(KNTYP0-1)+1
                      IF (INTYP0.EQ.KNTYP0) THEN
                         KNIONS=LNIONS_EST1(1,KNTYP0_FF,LOCVAR_COL)
                         PREFACTOR_FORCE=FFM%RCOUPLE_IONS(KNIONS)
                         PREFACTOR_STRESS=OMEGA_INV * FFM%RCOUPLE_IONS(KNIONS)
                         DO IXYZ=1, 3
                            INPRODUCT_HELP=DDOT(MLSHM%MROW_CLM1_GLOBAL,PRODUCT_L_DPDC1(1,ILOC,INTYP0), &
                                                1,FFM%DCLM1_HEAD(1,IXYZ,KNTYP0,ILOC,INTYP0),1)
                            ! Forces
                            TIFOR_MB(IXYZ,KNIONS) = TIFOR_MB(IXYZ,KNIONS) - PREFACTOR_FORCE * &
                                                    INPRODUCT_HELP
                            !! Stresses are 0 for first element in neighbor dimension because
                            !! XYZ_EST2 is zero for that element 
                         ENDDO 
                      ENDIF
                      DO KNEIB=2, NNEIB_EST1(KNTYP0_FF,LOCVAR_COL)+1 
                         KNIONS=LNIONS_EST1(KNEIB,KNTYP0_FF,LOCVAR_COL)
                         PREFACTOR_FORCE=FFM%RCOUPLE_IONS(KNIONS)
                         PREFACTOR_STRESS=OMEGA_INV * FFM%RCOUPLE_IONS(KNIONS)
                         DO IXYZ=1, 3
                            INPRODUCT_HELP=DDOT(FFM%IVAR_OFFSET1,PRODUCT_L_DPDC1(OFFSET,ILOC,INTYP0), &
                                                1,FFM%DCLM1(OFFSET,IXYZ,KNEIB,ILOC,INTYP0),1)
                            ! Forces
                            TIFOR_MB(IXYZ,KNIONS) = TIFOR_MB(IXYZ,KNIONS) - PREFACTOR_FORCE * &
                                                    INPRODUCT_HELP
                            ! Stresses
                            DO JXYZ=IXYZ, 3
                               TSIF_MB(JXYZ,IXYZ) = TSIF_MB(JXYZ,IXYZ) + PREFACTOR_STRESS * &
                                                    XYZ_EST1(KNEIB,JXYZ,KNTYP0_FF,LOCVAR_COL) * &
                                                    INPRODUCT_HELP
                            ENDDO
                         ENDDO 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             ! Angular descriptor
             IF (FFM%W2.GT.0.0_q) THEN
                DO ILOC=1, NITYP_LOCAL(INTYP0)
                   INIONS_LOC=LOCVAR_MAP_INIONS(ILOC,INTYP0)
                   LOCVAR_COL=CEILING(REAL(INIONS_LOC,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                   DO KNTYP0=1, ABN%NTYP
                      KNTYP0_FF=ABN%LADD(KNTYP0)
                      OFFSET=FFM%IVAR_OFFSET2*(KNTYP0-1)+1
                      IF (INTYP0.EQ.KNTYP0) THEN
                         KNIONS=LNIONS_EST2(1,KNTYP0_FF,LOCVAR_COL)
                         PREFACTOR_FORCE=FFM%RCOUPLE_IONS(KNIONS)
                         PREFACTOR_STRESS=OMEGA_INV * FFM%RCOUPLE_IONS(KNIONS)
                         DO IXYZ=1, 3
                            INPRODUCT_HELP=DDOT(MLSHM%MROW_CLM2_GLOBAL,PRODUCT_L_DPDC2(1,ILOC,INTYP0), &
                                                1,FFM%DCLM2_HEAD(1,IXYZ,KNTYP0,ILOC,INTYP0),1)
                            ! Forces
                            TIFOR_MB(IXYZ,KNIONS) = TIFOR_MB(IXYZ,KNIONS) - PREFACTOR_FORCE * &
                                                    INPRODUCT_HELP
                            !! Stresses are 0 for first element in neighbor dimension because
                            !! XYZ_EST2 is zero for that element 
                         ENDDO 
                      ENDIF
                      DO KNEIB=2, NNEIB_EST2(KNTYP0_FF,LOCVAR_COL)+1 
                         KNIONS=LNIONS_EST2(KNEIB,KNTYP0_FF,LOCVAR_COL)
                         PREFACTOR_FORCE=FFM%RCOUPLE_IONS(KNIONS)
                         PREFACTOR_STRESS=OMEGA_INV * FFM%RCOUPLE_IONS(KNIONS)
                         DO IXYZ=1, 3
                            INPRODUCT_HELP=DDOT(FFM%IVAR_OFFSET2,PRODUCT_L_DPDC2(OFFSET,ILOC,INTYP0), &
                                                1,FFM%DCLM2(OFFSET,IXYZ,KNEIB,ILOC,INTYP0),1)
                            ! Forces
                            TIFOR_MB(IXYZ,KNIONS) = TIFOR_MB(IXYZ,KNIONS) - PREFACTOR_FORCE * &
                                                    INPRODUCT_HELP
                            ! Stresses
                            DO JXYZ=IXYZ, 3
                               TSIF_MB(JXYZ,IXYZ) = TSIF_MB(JXYZ,IXYZ) + PREFACTOR_STRESS * &
                                                    XYZ_EST2(KNEIB,JXYZ,KNTYP0_FF,LOCVAR_COL) * &
                                                    INPRODUCT_HELP
                            ENDDO
                         ENDDO 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             PROFILING_STOP('ml_fast_fmat_calculate_force_stress')
          END SUBROUTINE CALCULATE_FORCE_STRESS_COMPONENTS_FAST

          SUBROUTINE CALCULATE_HEAT_COMPONENTS_FAST
             IMPLICIT NONE
             ! Local variables
             INTEGER :: IXYZ, JXYZ
             INTEGER :: ILOC, INIONS_LOC
             INTEGER :: KNTYP, KNTYP0, KNEIB
             INTEGER :: LOCVAR_COL
             INTEGER :: IVAR
             INTEGER :: OFFSET
             INTEGER :: KNTYP0_FF
             INTEGER :: NNEIB_EST_DUMMY
             REAL(q) :: PREFACTOR_STRESS
             REAL(q) :: INPRODUCT_HELP
             ! Lapack routines
             REAL(q) :: DDOT
             PROFILING_START('ml_fast_heat_help')
             ! radial part
             IF (FFM%W1.GT.0.0_q) THEN
                DO IXYZ=1, 3
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      INIONS_LOC=LOCVAR_MAP_INIONS(ILOC,INTYP0)
                      LOCVAR_COL=CEILING(REAL(INIONS_LOC,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                      DO KNTYP0=1, ABN%NTYP
                         KNTYP0_FF=ABN%LADD(KNTYP0)
                         OFFSET=FFM%IVAR_OFFSET1*(KNTYP0-1)+1
                         IF (INTYP0.EQ.KNTYP0) THEN
                            INPRODUCT_HELP=DDOT(MLSHM%MROW_CLM1_GLOBAL,PRODUCT_L_DPDC1(1,ILOC,INTYP0), &
                                                1,FFM%DCLM1_HEAD(1,IXYZ,KNTYP0,ILOC,INTYP0),1)
                            HEAT_HELP1(1,KNTYP0_FF,LOCVAR_COL,IXYZ)=-INPRODUCT_HELP
                         ENDIF
                         DO KNEIB=2, NNEIB_EST1(KNTYP0_FF,LOCVAR_COL)+1 
                            INPRODUCT_HELP=DDOT(FFM%IVAR_OFFSET1,PRODUCT_L_DPDC1(OFFSET,ILOC,INTYP0), &
                                                1,FFM%DCLM1(OFFSET,IXYZ,KNEIB,ILOC,INTYP0),1)
                            HEAT_HELP1(KNEIB,KNTYP0_FF,LOCVAR_COL,IXYZ)=-INPRODUCT_HELP
                         ENDDO 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             ! angular part
             IF (FFM%W2.GT.0.0_q) THEN
                DO IXYZ=1, 3
                   DO ILOC=1, NITYP_LOCAL(INTYP0)
                      INIONS_LOC=LOCVAR_MAP_INIONS(ILOC,INTYP0)
                      LOCVAR_COL=CEILING(REAL(INIONS_LOC,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                      DO KNTYP0=1, ABN%NTYP
                         KNTYP0_FF=ABN%LADD(KNTYP0)
                         OFFSET=FFM%IVAR_OFFSET2*(KNTYP0-1)+1
                         IF (INTYP0.EQ.KNTYP0) THEN
                            INPRODUCT_HELP=DDOT(MLSHM%MROW_CLM2_GLOBAL,PRODUCT_L_DPDC2(1,ILOC,INTYP0), &
                                                1,FFM%DCLM2_HEAD(1,IXYZ,KNTYP0,ILOC,INTYP0),1)
                            HEAT_HELP2(1,KNTYP0_FF,LOCVAR_COL,IXYZ)=-INPRODUCT_HELP
                         ENDIF
                         DO KNEIB=2, NNEIB_EST2(KNTYP0_FF,LOCVAR_COL)+1 
                            INPRODUCT_HELP=DDOT(FFM%IVAR_OFFSET2,PRODUCT_L_DPDC2(OFFSET,ILOC,INTYP0), &
                                                1,FFM%DCLM2(OFFSET,IXYZ,KNEIB,ILOC,INTYP0),1)
                            HEAT_HELP2(KNEIB,KNTYP0_FF,LOCVAR_COL,IXYZ)=-INPRODUCT_HELP
                         ENDDO 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             PROFILING_STOP('ml_fast_heat_help')
          END SUBROUTINE CALCULATE_HEAT_COMPONENTS_FAST


! Reduction of help array containing stress tensor elements over atom types
          SUBROUTINE ALLREDUCE_STRESS_AND_FORCES_FMAT_NEW_MB
             IMPLICIT NONE
! Local variables
             INTEGER               :: IERR
             INTEGER               :: INTYP
             INTEGER               :: INTYP0
             INTEGER               :: INTYP0_FF
             INTEGER               :: IXYZ, JXYZ
             PROFILING_START('ml_fast_allred_en_fo_str')
             ! Allreduce energy help energy coupling help
             DO INTYP0=1,ABN%NTYP
                INTYP0_FF=ABN%LADD(INTYP0)
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,ENERGY_HELP(1,1,INTYP0),FFM%NB(INTYP0_FF),MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
             ENDDO
             IF(FFM%LCOUPLE) THEN
                DO INTYP0=1,ABN%NTYP
                   INTYP0_FF=ABN%LADD(INTYP0)
                   CALL MPI_ALLREDUCE(MPI_IN_PLACE,ENERGY_COUPLE_HELP(1,1,INTYP0),FFM%NB(INTYP0_FF),MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
                ENDDO
             ENDIF
             ! Allreduce potential energies
             IF (FFM%LEATOM.OR.FFM%LHEAT.OR.FFM%LCOUPLE) THEN
                CALL MPI_ALLREDUCE (MPI_IN_PLACE,ABN%EATOM_POT,ABN%NIONS,MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
             ENDIF 
             ! Allreduce forces
             CALL MPI_ALLREDUCE(MPI_IN_PLACE,TIFOR_MB(1,1),3*SIZE(TIFOR_MB,2),MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
             ! Allreduce stress
             CALL MPI_ALLREDUCE(MPI_IN_PLACE,TSIF_MB(1,1),9,MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
             ! Symmetrization of stress tensor.
             DO IXYZ=1, 2
                DO JXYZ=IXYZ+1, 3
                   TSIF_MB(IXYZ,JXYZ)=TSIF_MB(JXYZ,IXYZ)
                ENDDO
             ENDDO
             PROFILING_STOP('ml_fast_allred_en_fo_str')
          END SUBROUTINE ALLREDUCE_STRESS_AND_FORCES_FMAT_NEW_MB

        END SUBROUTINE FMAT_NEW_MB3


        SUBROUTINE CALCULATE_HELPING_MAPS(ABN,FF,FFM,MLSHM,PAR_SUP_HANDLE)
           USE ML_FF_STRUCT, ONLY: ABINITIONEW, FF_PAR, FFM_PAR, ML_SHMEM, PARALLEL_SUPER
           IMPLICIT NONE
! Input and output variables
           TYPE (ABINITIONEW)    :: ABN
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
! Local variableso   
           INTEGER               :: INTYP, INTYP0, INITYP
           INTEGER               :: JNTYP, JNTYP0, JJNTYP, JJNTYP0, ILNRB
           INTEGER               :: INTYP0_FF, JNTYP0_CLM_FF, JNTYP0_FF, JJNTYP0_FF
           INTEGER               :: L, IB, JB, M
           INTEGER               :: IVAR1, IVAR2, IVAR3, IVAR4
           INTEGER               :: JNTYP_CLM, JNTYP0_CLM, L_CLM, IB_CLM, M_CLM
           REAL(q)               :: PI8, FAC_HELP
           REAL(q)               :: HELP_PF(0:FFM%LMAX2)

           PROFILING_START('ml_fast_calculate_helping_maps')

           ! Set constant
           PI8=8.0_q*PI**2

           ! Radial descriptor
           IF (FFM%W1.GT.0) THEN
              ! Initialization
              FFM%PREFACTOR_HELP1=0.0_q
              ! Central atom
              DO INTYP0=1, ABN%NTYP
                 ! derivative C00
                 DO JNTYP0_CLM=1, ABN%NTYP
                    JNTYP0_CLM_FF=ABN%LADD(JNTYP0_CLM)
                    DO IB_CLM=1, FFM%NRB1(0)
                       IVAR2=FFM%MRB1*(JNTYP0_CLM_FF-1)+IB_CLM
                       IVAR1=FFM%IVAR_OFFSET1*(JNTYP0_CLM-1) + FFM%OFFSET_MAP1(IB_CLM,1,0)
                       FFM%PREFACTOR_HELP1(IVAR2,IVAR1,INTYP0)=1.0_q
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF

           ! Angular descriptor
           IF (FFM%W2.GT.0) THEN
              ! Define helping factor
              FAC_HELP=1.0_q/(PI8*REAL(FFM%NRB2(0),q)*REAL(FFM%LMAX2+1,q))**2 
              !Calculate l dependent prefactor
              DO L=0, FFM%LMAX2
                 IF(FFM%LAFILT2) THEN
                    IF(FFM%IAFILT2.EQ.1) THEN
                       HELP_PF(L)=SQRT(PI8/REAL((2*L+1),q))/SQRT(REAL((2*L+1),q))
                    ELSE IF(FFM%IAFILT2.EQ.2) THEN
                       HELP_PF(L)=SQRT(PI8/REAL((2*L+1),q))/(1.0_q+FFM%AFILT2*REAL((L*(L+1)),q)**2)**2
                    ENDIF
                 ELSE
                    HELP_PF(L)=SQRT(PI8/REAL((2*L+1),q))
                 ENDIF
              ENDDO
              ! Initialize descriptor index counter 
              FFM%INDEX_NUM2=0
              ! Initialize target arrays (Normal and SIC)
              FFM%INDEX_PAIR_HELP2=1
              FFM%PREFACTOR_HELP2=0.0_q
              ! Set defaults for values if linear-scaling descriptor is used
              IF (FFM%DESC_TYPE.EQ.1) THEN
                 FFM%INDEX_NUM2_B=0
                 FFM%INDEX_PAIR_HELP2_B=1
                 FFM%PREFACTOR_HELP2_B=0.0_q
              ENDIF
              IF (FFM%DESC_TYPE.EQ.0) THEN
                 ! Loop over descriptors, we run it for all central atom types,
                 ! possibly loop is enough for just one type
                 DO INTYP0=1, ABN%NTYP
                    INTYP0_FF=ABN%LADD(INTYP0)
                    DO JNTYP0_CLM=1, ABN%NTYP
                       JNTYP0_CLM_FF=ABN%LADD(JNTYP0_CLM)
                       !IVAR1=FFM%IVAR_OFFSET2*(JNTYP0_CLM-1)
                       DO L_CLM=0,FFM%LMAX2
                          DO M_CLM=1, 2*L_CLM+1
                             DO IB_CLM=1, FFM%NRB2(L_CLM)
                                !IVAR1=IVAR1+1
                                IVAR1=FFM%IVAR_OFFSET2*(JNTYP0_CLM-1) + FFM%OFFSET_MAP2(IB_CLM,M_CLM,L_CLM)
                                ! Loop over descriptoors and CLMs (non-derivatives)
                                DO JNTYP0=1, ABN%NTYP
                                   JNTYP0_FF=ABN%LADD(JNTYP0)
                                   DO JJNTYP0=1, ABN%NTYP
                                      JJNTYP0_FF=ABN%LADD(JJNTYP0)
                                      DO ILNRB=1, FFM%NLNRB2
                                         L=FFM%L_LNRB2(ILNRB)
                                         IB=FFM%NRB_LNRB2(ILNRB)
                                         IF (L_CLM.EQ.L) THEN
                                            DO JB=IB, FFM%NRB2(L)
                                               IF(FFM%LFLAG_VAR(JB,IB,L,JJNTYP0_FF,JNTYP0_FF,INTYP0_FF)) THEN
                                                  IVAR3=FFM%LVAR(JB,IB,L,JJNTYP0_FF,JNTYP0_FF,INTYP0_FF)
                                                  IVAR2=FFM%IVAR_OFFSET2*(JJNTYP0-1) + FFM%OFFSET_MAP2(JB,M_CLM,L_CLM)
                                                  IVAR4=FFM%IVAR_OFFSET2*(JNTYP0-1) + FFM%OFFSET_MAP2(IB,M_CLM,L_CLM)
                                                  IF (IVAR1.EQ.IVAR4 .OR. IVAR1.EQ.IVAR2) THEN
                                                     ! Increase number of considered indices for a given IVAR1 
                                                     FFM%INDEX_NUM2(IVAR1,INTYP0)=FFM%INDEX_NUM2(IVAR1,INTYP0) + 1
                                                     ! Have the index back to the whole powerspectrum
                                                     FFM%INDEX_MAP_TO_PS(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR3
                                                     ! Set the factor used for a given index
                                                     IF (IB.EQ.JB .AND. JNTYP0.EQ.JJNTYP0) THEN
                                                        FFM%PREFACTOR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=2.0_q*HELP_PF(L_CLM)
                                                     ELSE
                                                        IF (IB.EQ.JB .AND. JNTYP0.NE.JJNTYP0) THEN
                                                           FFM%PREFACTOR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=HELP_PF(L_CLM)
                                                        ELSE
                                                           FFM%PREFACTOR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=SQRT(2.0_q)*HELP_PF(L_CLM)
                                                        ENDIF
                                                     ENDIF
                                                     ! Actual index pairs are calculated here 
                                                     IF (IVAR1.EQ.IVAR4) THEN
                                                        FFM%INDEX_PAIR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR2
                                                     ELSE
                                                        FFM%INDEX_PAIR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR4
                                                     ENDIF
                                                  ENDIF
                                               ENDIF
                                            ENDDO
                                         ENDIF
                                      ENDDO
                                   ENDDO
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ELSE IF (FFM%DESC_TYPE.EQ.1) THEN
                 ! Loop over descriptors, we run it for all central atom types,
                 ! possibly loop is enough for just one type
                 DO INTYP0=1, ABN%NTYP
                    INTYP0_FF=ABN%LADD(INTYP0)
                    DO JNTYP0_CLM=1, ABN%NTYP
                       JNTYP0_CLM_FF=ABN%LADD(JNTYP0_CLM)
                       DO L_CLM=0,FFM%LMAX2
                          DO M_CLM=1, 2*L_CLM+1
                             DO IB_CLM=1, FFM%NRB2(L_CLM)
                                IVAR1=FFM%IVAR_OFFSET2*(JNTYP0_CLM-1) + FFM%OFFSET_MAP2(IB_CLM,M_CLM,L_CLM)
                                ! Loop over descriptoors and CLMs (non-derivatives)
                                DO JNTYP0=1, ABN%NTYP
                                   JNTYP0_FF=ABN%LADD(JNTYP0)
                                   DO ILNRB=1, FFM%NLNRB2
                                      L=FFM%L_LNRB2(ILNRB)
                                      IB=FFM%NRB_LNRB2(ILNRB)
                                      IF (L_CLM.EQ.L) THEN
                                         DO JB=IB, FFM%NRB2(L)
                                            IF(FFM%LFLAG_VAR(JB,IB,L,1,JNTYP0_FF,INTYP0_FF)) THEN
                                               IVAR4=FFM%IVAR_OFFSET2*(JNTYP0-1) + FFM%OFFSET_MAP2(IB,M_CLM,L_CLM)
                                               IVAR3=FFM%LVAR(JB,IB,L,1,JNTYP0_FF,INTYP0_FF)
                                               IF (IVAR1.EQ.IVAR4) THEN
                                                  ! Increase number of considered indices for a given IVAR1 
                                                  FFM%INDEX_NUM2_B(IVAR1,INTYP0)=FFM%INDEX_NUM2_B(IVAR1,INTYP0) + 1
                                                  FFM%INDEX_MAP_TO_PS_B(FFM%INDEX_NUM2_B(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR3
                                                  !Prefactors are calculated here, factor SQRT(2) due to upper triangle in radial basis functions
                                                  IF (IB.EQ.JB) THEN
                                                     FFM%PREFACTOR_HELP2_B(FFM%INDEX_NUM2_B(IVAR1,INTYP0),IVAR1,INTYP0)=HELP_PF(L_CLM)
                                                  ELSE
                                                     FFM%PREFACTOR_HELP2_B(FFM%INDEX_NUM2_B(IVAR1,INTYP0),IVAR1,INTYP0)=SQRT(2.0_q)*HELP_PF(L_CLM)
                                                  ENDIF
                                                  DO JJNTYP0=1, ABN%NTYP
                                                     IVAR2=FFM%IVAR_OFFSET2*(JJNTYP0-1) + FFM%OFFSET_MAP2(JB,M_CLM,L_CLM)
                                                     ! Actual index pairs are calculated here 
                                                     FFM%INDEX_PAIR_HELP2_B(JJNTYP0,FFM%INDEX_NUM2_B(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR2
                                                  ENDDO
                                               ENDIF
                                                 
                                               DO JJNTYP0=1, ABN%NTYP
                                                  IVAR2=FFM%IVAR_OFFSET2*(JJNTYP0-1) + FFM%OFFSET_MAP2(JB,M_CLM,L_CLM)
                                                  IF (IVAR1.EQ.IVAR2) THEN
                                                     ! Also increase index on the right hand side
                                                     FFM%INDEX_NUM2(IVAR1,INTYP0)=FFM%INDEX_NUM2(IVAR1,INTYP0) + 1
                                                     FFM%INDEX_MAP_TO_PS(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR3
                                                     IF (IB.EQ.JB) THEN
                                                        FFM%PREFACTOR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=HELP_PF(L_CLM)
                                                     ELSE
                                                        FFM%PREFACTOR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=SQRT(2.0_q)*HELP_PF(L_CLM)
                                                     ENDIF
                                                     ! Actual index pair is calculated here
                                                     FFM%INDEX_PAIR_HELP2(FFM%INDEX_NUM2(IVAR1,INTYP0),IVAR1,INTYP0)=IVAR4
                                                     ! Exit JJNTYP loop
                                                     EXIT
                                                  ENDIF
                                               ENDDO
                                            ENDIF
                                         ENDDO
                                      ENDIF
                                   ENDDO
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF

           PROFILING_STOP('ml_fast_calculate_helping_maps')

        END SUBROUTINE CALCULATE_HELPING_MAPS

        SUBROUTINE PS_ALL_MULTIPLY_WITH_WEIGHTS(FF,FFM,MLSHM,PAR_SUP_HANDLE)
           USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, ML_SHMEM, PARALLEL_SUPER
           IMPLICIT NONE
! Input and output variables
           TYPE (FF_PAR)         :: FF
           TYPE (FFM_PAR)        :: FFM
           TYPE (ML_SHMEM)       :: MLSHM
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
! Local variables
           INTEGER               :: IERR
           INTEGER               :: IHEAD
           INTEGER               :: ILOCCONF
           INTEGER               :: INTYP0
           INTEGER               :: IVAR
           INTEGER               :: LHEAD(1:FF%MTYP)
           
           PROFILING_START('ml_fast_pslocconf_tims_weights')
           ! Calculate offset for local reference numbering
           LHEAD=1
           DO INTYP0=2, FF%MTYP, 1
              LHEAD(INTYP0)=LHEAD(INTYP0-1)+FFM%NB(INTYP0-1)
           ENDDO
#ifdef use_shmem
           IF (PAR_SUP_HANDLE%COMM_SHMEM%NODE_ME.EQ.0) THEN
#endif
              ! Radial part
              IF (FFM%W1.GT.0.0_q) THEN
                 DO INTYP0=1, FF%MTYP
                    IHEAD=LHEAD(INTYP0)
                    DO ILOCCONF=1, FFM%NB(INTYP0)
                       DO IVAR=1, FFM%NNVAR1
#ifdef scaLAPACK
                          FFM%C00_W_ALL(IVAR,IHEAD+ILOCCONF-1)=FFM%C00_ALL(IVAR,IHEAD+ILOCCONF-1)*FFM%WMAT_ALL(IHEAD+ILOCCONF-1)
#else
                          FFM%C00_W_ALL(IVAR,IHEAD+ILOCCONF-1)=FFM%C00(IVAR,ILOCCONF,INTYP0)*FFM%WMAT(ILOCCONF,1,INTYP0)
#endif
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              ! Angular part
              IF (FFM%W2.GT.0.0_q) THEN
                 DO INTYP0=1, FF%MTYP
                    IHEAD=LHEAD(INTYP0)
                    DO ILOCCONF=1, FFM%NB(INTYP0)
                       DO IVAR=1, FFM%NNVAR2(INTYP0)
#ifdef scaLAPACK
                          FFM%PS_W_ALL(IVAR,IHEAD+ILOCCONF-1)=FFM%PS_ALL(IVAR,IHEAD+ILOCCONF-1)*FFM%WMAT_ALL(IHEAD+ILOCCONF-1)
#else
                          FFM%PS_W_ALL(IVAR,IHEAD+ILOCCONF-1)=FFM%PS(IVAR,ILOCCONF,INTYP0)*FFM%WMAT(ILOCCONF,1,INTYP0)
#endif
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
#ifdef use_shmem
           ENDIF
           CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_SHMEM%MPI_COMM,IERR)
#endif
           PROFILING_STOP('ml_fast_pslocconf_tims_weights')
        END SUBROUTINE PS_ALL_MULTIPLY_WITH_WEIGHTS

!------------------------------------------------------------------------------------------------------------------
! Calculate energy, force and stress tensor within fast version method.
! In this method, regression coefficients stored in all processes are 
! used to calculate energy, force and stress tensor without using scaLAPACK.
!------------------------------------------------------------------------------------------------------------------
        SUBROUTINE ZMAT_NEW_MB_NEW(ABN,PAR_SUP_HANDLE,FF,FFM,LHEAD, &
                   LNIONS_EST1,LNIONS_EST2,LNIONS_EST_MAX,  &
                   MNEIB_EST1,MNEIB_EST2,MNEIB_EST_MAX,XYZ_EST1,XYZ_EST2,XYZ_EST_MAX,EATOM_HELP, &
                   ENERGY_COUPLE_HELP,ENERGY_HELP,HEAT_HELP1,HEAT_HELP2,NNEIB_EST1,NNEIB_EST2, &
                   NNEIB_EST_MAX,QHEAT,RMEM,TIFOR_MB,TOTEN_MB,VEL)
          USE ML_FF_STRUCT, ONLY: ABINITIONEW, PARALLEL_SUPER, ML_SHMEM, FF_PAR, FFM_PAR
          IMPLICIT NONE
! Input and output variables
          TYPE (ABINITIONEW)    :: ABN
          TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
          TYPE (FF_PAR)         :: FF
          TYPE (FFM_PAR)        :: FFM
          INTEGER, INTENT(IN)   :: LHEAD(:)
          INTEGER, INTENT(IN)   :: LNIONS_EST1(:,:,:) !(1:MROW_EST1,1:MTYP,1:MCOL_EST1)
          INTEGER, INTENT(IN)   :: LNIONS_EST2(:,:,:) !(1:MROW_EST2,1:MTYP,1:MCOL_EST2)
          INTEGER, INTENT(IN)   :: LNIONS_EST_MAX(:,:,:)
          INTEGER, INTENT(IN)   :: MNEIB_EST1
          INTEGER, INTENT(IN)   :: MNEIB_EST2
          INTEGER, INTENT(IN)   :: MNEIB_EST_MAX
          REAL(q), INTENT(IN)   :: XYZ_EST1(:,:,:,:) !(1:MROW_EST1,1:3,1:MTYP,1:MCOL_EST1)
          REAL(q), INTENT(IN)   :: XYZ_EST2(:,:,:,:) !(1:MROW_EST2,1:3,1:MTYP,1:MCOL_EST2)
          REAL(q), INTENT(IN)   :: XYZ_EST_MAX(:,:,:,:) !(1:3,1:MROW_EST_MAX,1:MTYP,1:MCOL_EST_MAX)
          REAL(q), INTENT(INOUT)  :: EATOM_HELP(:,:,:) !(1:MROW_EATOM,1:MCOL_EATOM,1:FF%MTYP)
          REAL(q), INTENT(INOUT)  :: ENERGY_COUPLE_HELP(:,:,:) !(1:MROW_ENERGY_COUPLE_HELP,1:MCOL_ENERGY_COUPLE_HELP,1:MTYP)
          REAL(q), INTENT(INOUT)  :: ENERGY_HELP(:,:,:) !(1:MROW_ENERGY_HELP,1:MCOL_ENERGY_HELP,1:MTYP)
          REAL(q), INTENT(INOUT)  :: HEAT_HELP1(:,:,:,:) !(1:MROW_HEAT1,1:MCOL_HEAT1,1:FF%MTYP,1:MCOL_MIONS,1:3)
          REAL(q), INTENT(INOUT)  :: HEAT_HELP2(:,:,:,:) !(1:MROW_HEAT2,1:MCOL_HEAT2,1:FF%MTYP,1:MCOL_MIONS,1:3)
          INTEGER, INTENT(INOUT)  :: NNEIB_EST1(:,:) !(1:MTYP,1:MCOL_EST1)
          INTEGER, INTENT(INOUT)  :: NNEIB_EST2(:,:) !(1:MTYP,1:MCOL_EST2)
          INTEGER, INTENT(INOUT)  :: NNEIB_EST_MAX(:,:) !(1:MTYP,1:MCOL_EST_MAX)
          REAL(q), INTENT(INOUT)  :: QHEAT(:) !(1:3)
          REAL(q), INTENT(INOUT)  :: RMEM
          REAL(q), INTENT(INOUT)  :: TIFOR_MB(:,:) !(1:3,1:MIONS)
          REAL(q), INTENT(INOUT)  :: TOTEN_MB
          REAL(q), INTENT(INOUT)  :: VEL(:,:) !(1:3,1:MIONS)
! Local variables
          REAL(q)               :: DRIFT(1:3)
          INTEGER               :: IERR
          INTEGER               :: IHEAD
          INTEGER               :: INEIB
          INTEGER               :: INIONS, INIONS_TMP
          INTEGER               :: INITYP, INITYP_TMP
          INTEGER               :: INTYP, INTYP0, INTYP0_FF
          INTEGER               :: IXYZ, JXYZ, KXYZ
          INTEGER               :: KNEIB
          INTEGER               :: KNIONS, KNIONS_TMP
          INTEGER               :: KNTYP, KNTYP0, KNTYP0_FF
          INTEGER               :: LOCVAR_COL
          INTEGER               :: MCOL_VEC_LOC_HEAT
          INTEGER               :: MROW_VEC_LOC_HEAT
          INTEGER               :: NNEIB_EST_DUMMY
          REAL(q)               :: TOTEN_ATOM_REF
          REAL(q), ALLOCATABLE  :: VEC_LOC_HEAT(:,:)
! LAPACK routine.
          REAL(q) :: DDOT
! Start profiling
          PROFILING_START ('zmat_new_mb')
          IF(FF%IFF) THEN
! Allocate necessary arrays in this subroutine
             IF(FFM%LHEAT) THEN
! Size of VEC_LOC_HEAT.
                MROW_VEC_LOC_HEAT=1
                MCOL_VEC_LOC_HEAT=1
                IF(FFM%W1.GT.0.0_q) THEN
                   IF(.NOT.FFM%LSUPERVEC) THEN
                      MROW_VEC_LOC_HEAT=MAX(MNEIB_EST1+1,MROW_VEC_LOC_HEAT)
                   ELSE
                      MROW_VEC_LOC_HEAT=MAX(MNEIB_EST_MAX+1,MROW_VEC_LOC_HEAT)
                   ENDIF
                ENDIF
                IF(FFM%W2.GT.0.0_q) THEN
                   IF(.NOT.FFM%LSUPERVEC) THEN
                      MROW_VEC_LOC_HEAT=MAX(MNEIB_EST2+1,MROW_VEC_LOC_HEAT)
                   ELSE
                      MROW_VEC_LOC_HEAT=MAX(MNEIB_EST_MAX+1,MROW_VEC_LOC_HEAT)
                   ENDIF
                ENDIF
                IF(ALLOCATED(VEC_LOC_HEAT)) THEN
                   DEALLOCATE(VEC_LOC_HEAT)
                ENDIF
                ALLOCATE (VEC_LOC_HEAT(1:MROW_VEC_LOC_HEAT,1:MCOL_VEC_LOC_HEAT))
                RMEM=RMEM+REAL(MROW_VEC_LOC_HEAT,q)*REAL(MCOL_VEC_LOC_HEAT,q)*8.0_q/1.0E6_q
             ENDIF
                
! Initialization of energy, force and stress
             TOTEN_MB=0.0_q
             IF(FFM%LCOUPLE) THEN
                FFM%DHDRCOUPLE=0.0_q
             ENDIF
! Calculate reference atomic energy.
             TOTEN_ATOM_REF=0.0_q
             DO INTYP=1, ABN%NTYP
                INTYP0_FF=ABN%LADD(INTYP)
                DO INITYP=1, ABN%NITYP(INTYP0_FF)
                   TOTEN_ATOM_REF=TOTEN_ATOM_REF+FF%EATOM_REF(INTYP0_FF)
                ENDDO
             ENDDO

! Calculate energy, force and stress tensor.
! Calculate energy.
! IF statement is necessary since the saved arrays for the fast linear
! descriptor are different from the usual descriptors.
             DO INTYP0=1, ABN%NTYP
                INTYP0_FF=ABN%LADD(INTYP0)
                IHEAD=LHEAD(INTYP0_FF)
#ifdef scaLAPACK
                IF (FF%LTOTEN_SYSTEM) THEN
                   TOTEN_MB=TOTEN_MB+DDOT(FFM%NB(INTYP0_FF),ENERGY_HELP(1,1,INTYP0),1,FFM%WMAT_ALL(IHEAD),1)
                ELSE
                   TOTEN_MB=TOTEN_MB+REAL(ABN%NIONS,q)*DDOT(FFM%NB(INTYP0_FF),ENERGY_HELP(1,1,INTYP0),1,FFM%WMAT_ALL(IHEAD),1)
                ENDIF
#else
                IF (FF%LTOTEN_SYSTEM) THEN
                   TOTEN_MB=TOTEN_MB+DDOT(FFM%NB(INTYP0_FF),ENERGY_HELP(1,1,INTYP0),1,FFM%WMAT(1,1,INTYP0_FF),1)
                ELSE
                   TOTEN_MB=TOTEN_MB+REAL(ABN%NIONS,q)*DDOT(FFM%NB(INTYP0_FF),ENERGY_HELP(1,1,INTYP0),1,FFM%WMAT(1,1,INTYP0_FF),1)
                ENDIF
#endif
! Calculate the derivative of the energy with respect to the coupling constant.
                IF(FFM%LCOUPLE) THEN
#ifdef scaLAPACK
                   IF (FF%LTOTEN_SYSTEM) THEN
                      FFM%DHDRCOUPLE=FFM%DHDRCOUPLE+DDOT(FFM%NB(INTYP0_FF),ENERGY_COUPLE_HELP(1,1,INTYP0),1,FFM%WMAT_ALL(IHEAD),1)
                   ELSE
                      FFM%DHDRCOUPLE=FFM%DHDRCOUPLE+REAL(ABN%NIONS,q)*DDOT(FFM%NB(INTYP0_FF),ENERGY_COUPLE_HELP(1,1,INTYP0),1,FFM%WMAT_ALL(IHEAD),1)
                   ENDIF
#else
                   IF (FF%LTOTEN_SYSTEM) THEN
                      FFM%DHDRCOUPLE=FFM%DHDRCOUPLE+DDOT(FFM%NB(INTYP0_FF),ENERGY_COUPLE_HELP(1,1,INTYP0),1,FFM%WMAT(1,1,INTYP0_FF),1)
                   ELSE
                      FFM%DHDRCOUPLE=FFM%DHDRCOUPLE+REAL(ABN%NIONS,q)*DDOT(FFM%NB(INTYP0_FF),ENERGY_COUPLE_HELP(1,1,INTYP0),1,FFM%WMAT(1,1,INTYP0_FF),1)
                   ENDIF
#endif
                ENDIF
             ENDDO
! If FFM%ISCALE_TOTEN=1, the energy data is scaled to the energy of isolated atoms. Thus, this should be added.
             IF(FFM%ISCALE_TOTEN.EQ.1) THEN
                TOTEN_MB=TOTEN_MB+TOTEN_ATOM_REF
! If FFM%ISCALE_TOTEN=2, the energy data is scaled to the average of training energy data. Thus, the average energy must be also added.
             ELSE IF (FFM%ISCALE_TOTEN.EQ.2) THEN
                DO INTYP=1, ABN%NTYP
                   INTYP0_FF=ABN%LADD(INTYP)
                   DO INITYP=1, ABN%NITYP(INTYP0_FF)
                      INIONS=ABN%LNIONS(INITYP,INTYP0_FF)
                      TOTEN_MB=TOTEN_MB+FF%EATOM_REF(INTYP0_FF)+FFM%RCOUPLE_IONS(INIONS)*FFM%TOTENAV
                   ENDDO
                ENDDO
             ENDIF
! Calculate atomic energies, if necessary.
!TEST
!             IF(FFM%LHEAT.OR.FFM%LEATOM.OR.FFM%LCOUPLE) THEN
             IF(FFM%LHEAT.OR.FFM%LEATOM) THEN
!TEST
! Initialization.
                ABN%EATOM_KIN(1:ABN%NIONS)=0.0_q
! Drift term in forces.
                DRIFT(1:3)=0.0_q
                DO INIONS=1, ABN%NIONS
                   DO IXYZ=1, 3
                      DRIFT(IXYZ)=DRIFT(IXYZ)+TIFOR_MB(IXYZ,INIONS)
                   ENDDO
                ENDDO
                DRIFT(1:3)=-DRIFT(1:3)/REAL(ABN%NIONS,q)
! Renew the velocities by using the calculated forces.
                DO INTYP=1, ABN%NTYP
                   INTYP0_FF=ABN%LADD(INTYP)
                   DO INITYP=1, ABN%NITYP(INTYP0_FF)
                      INIONS=ABN%LNIONS(INITYP,INTYP0_FF)
                      DO IXYZ=1, 3
                         VEL(IXYZ,INIONS)=VEL(IXYZ,INIONS)+5E-01_q*(TIFOR_MB(IXYZ,INIONS)+DRIFT(IXYZ))/ABN%POMASS(INTYP0_FF)*ABN%POTIM
                      ENDDO
                   ENDDO
                ENDDO
! Calculate the atomic energies.
                DO INTYP=1, ABN%NTYP
                   INTYP0_FF=ABN%LADD(INTYP)
                   IHEAD=LHEAD(INTYP0_FF)
                   DO INITYP_TMP=1, ABN%NITYP(INTYP0_FF), PAR_SUP_HANDLE%COMM_WORLD%NCPU
                      INITYP=INITYP_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                      IF(INITYP.LE.ABN%NITYP(INTYP0_FF)) THEN
                         INIONS=ABN%LNIONS(INITYP,INTYP0_FF)
                         DO IXYZ=1, 3
                            ABN%EATOM_KIN(INIONS)=ABN%EATOM_KIN(INIONS)+5E-01_q*ABN%POMASS(INTYP0_FF)*VEL(IXYZ,INIONS)**2
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
                CALL MPI_ALLREDUCE (MPI_IN_PLACE,ABN%EATOM_KIN,ABN%NIONS,MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
             ENDIF
! If FFM%LCOUPLE=.TRUE., calculate the derivative of the Hamiltonian with respect to the coupling parameter.
             IF(FFM%LCOUPLE) THEN
                DO INIONS=1, FFM%NATOM_COUPLED
                   INTYP0_FF=ABN%LNTYP (FFM%ICOUPLE(INIONS))
                   FFM%DHDRCOUPLE=FFM%DHDRCOUPLE+ABN%EATOM_POT(FFM%ICOUPLE(INIONS))-FF%EATOM_REF(INTYP0_FF)
                ENDDO
             ENDIF
! Calculate heat flux, if necessary.
             IF(FFM%LHEAT) THEN
! Initialization of heat flux.
                QHEAT=0.0_q
! Calculations of interaction heat flux.
                IF(FFM%W1.GT.0.0_q) THEN
                   DO IXYZ=1, 3
                      DO INIONS_TMP=1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
                         INIONS=INIONS_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                         IF(INIONS.LE.ABN%NIONS) THEN
                            LOCVAR_COL=CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                            INTYP0=ABN%TYPE_OF_EACH_ATOM_IN_ABN(INIONS)
                            DO KNTYP0=1, ABN%NTYP
                               KNTYP0_FF=ABN%LADD(KNTYP0)
                               DO JXYZ=1, 3
                                  NNEIB_EST_DUMMY=NNEIB_EST1(KNTYP0_FF,LOCVAR_COL)
                                  DO KNEIB=1, NNEIB_EST_DUMMY+1
                                     KNIONS=LNIONS_EST1(KNEIB,KNTYP0_FF,LOCVAR_COL)
                                     IF((KNEIB.EQ.1).AND.(KNTYP0.NE.INTYP0)) THEN
                                        VEC_LOC_HEAT(KNEIB,1)=0.0_q
                                     ELSE
                                        VEC_LOC_HEAT(KNEIB,1)=XYZ_EST1(KNEIB,IXYZ,KNTYP0_FF,LOCVAR_COL) * &
                                                     VEL(JXYZ,KNIONS) * FFM%RCOUPLE_IONS(KNIONS)
                                     ENDIF
                                  ENDDO
                                  QHEAT(IXYZ)=QHEAT(IXYZ)+DDOT(NNEIB_EST_DUMMY+1,HEAT_HELP1(1,KNTYP0_FF,LOCVAR_COL,JXYZ),1,VEC_LOC_HEAT(1,1),1)
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
                IF(FFM%W2.GT.0.0_q) THEN
                   DO IXYZ=1, 3
                      DO INIONS_TMP=1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
                         INIONS=INIONS_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                         IF(INIONS.LE.ABN%NIONS) THEN
                            LOCVAR_COL=CEILING(REAL(INIONS,q)/REAL(PAR_SUP_HANDLE%COMM_WORLD%NCPU,q))
                            INTYP0=ABN%TYPE_OF_EACH_ATOM_IN_ABN(INIONS)
                            DO KNTYP0=1, ABN%NTYP
                               KNTYP0_FF=ABN%LADD(KNTYP0)
                               DO JXYZ=1, 3
                                  NNEIB_EST_DUMMY=NNEIB_EST2(KNTYP0_FF,LOCVAR_COL)
                                  DO KNEIB=1, NNEIB_EST_DUMMY+1
                                     KNIONS=LNIONS_EST2(KNEIB,KNTYP0_FF,LOCVAR_COL)
                                     IF((KNEIB.EQ.1).AND.(KNTYP0.NE.INTYP0)) THEN
                                        VEC_LOC_HEAT(KNEIB,1)=0.0_q
                                     ELSE
                                        VEC_LOC_HEAT(KNEIB,1)=XYZ_EST2(KNEIB,IXYZ,KNTYP0_FF,LOCVAR_COL) * &
                                                     VEL(JXYZ,KNIONS) * FFM%RCOUPLE_IONS(KNIONS)
                                     ENDIF
                                  ENDDO
                                  QHEAT(IXYZ)=QHEAT(IXYZ)+DDOT(NNEIB_EST_DUMMY+1,HEAT_HELP2(1,KNTYP0_FF,LOCVAR_COL,JXYZ),1,VEC_LOC_HEAT(1,1),1)
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
! Calculations of convection heat flux.
                DO IXYZ=1, 3
                   DO KNIONS_TMP=1, ABN%NIONS, PAR_SUP_HANDLE%COMM_WORLD%NCPU
                      KNIONS=KNIONS_TMP+PAR_SUP_HANDLE%COMM_WORLD%NODE_ME
                      IF(KNIONS.LE.ABN%NIONS) THEN
                         QHEAT(IXYZ)=QHEAT(IXYZ) +(ABN%EATOM_KIN(KNIONS)+ABN%EATOM_POT(KNIONS))*VEL(IXYZ,KNIONS)
                      ENDIF
                   ENDDO
                ENDDO
! Then, allreduce the flux.
                CALL MPI_ALLREDUCE (MPI_IN_PLACE,QHEAT,3,MPI_REAL8,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
             ENDIF
! Deallocate helping arrays for heat flux
             IF(FFM%LHEAT) THEN
                DEALLOCATE(VEC_LOC_HEAT)
             ENDIF
          ELSE
             TOTEN_MB=0.0_q
             IF(FFM%LEATOM.OR.FFM%LHEAT.OR.FFM%LCOUPLE) THEN
                ABN%EATOM_KIN=0.0_q
             ENDIF
             IF(FFM%LHEAT) THEN
                QHEAT=0.0_q
             ENDIF
          ENDIF
! Stop profiling
          PROFILING_STOP ('zmat_new_mb')
        END SUBROUTINE ZMAT_NEW_MB_NEW

        SUBROUTINE SET_IVAR_OFFSET(FFM)
           USE ML_FF_STRUCT, ONLY: FFM_PAR, ML_IO_WRITE
           IMPLICIT NONE
           ! Input and output variables
           TYPE (FFM_PAR)        :: FFM
           ! Local variables
           INTEGER               :: L, M, IB
           IF (FFM%W1.GT.0.0_q) THEN
              IF (ALLOCATED(FFM%OFFSET_MAP1)) DEALLOCATE(FFM%OFFSET_MAP1)
              ALLOCATE(FFM%OFFSET_MAP1(1:FFM%MRB1,1:(2*FFM%LMAX1+1),0:FFM%LMAX1))
              FFM%IVAR_OFFSET1=0
              DO L=0, FFM%LMAX1
                 DO M=1, 2*L+1
                    DO IB=1, FFM%NRB1(L)
                       FFM%IVAR_OFFSET1=FFM%IVAR_OFFSET1+1
                       FFM%OFFSET_MAP1(IB,M,L)=FFM%IVAR_OFFSET1
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
              
           IF (FFM%W2.GT.0.0_q) THEN
              IF (ALLOCATED(FFM%OFFSET_MAP2)) DEALLOCATE(FFM%OFFSET_MAP2)
              ALLOCATE(FFM%OFFSET_MAP2(1:FFM%MRB2,1:(2*FFM%LMAX2+1),0:FFM%LMAX2))
              FFM%IVAR_OFFSET2=0
              DO L=0, FFM%LMAX2
                 DO M=1, 2*L+1
                    DO IB=1, FFM%NRB2(L)
                       FFM%IVAR_OFFSET2=FFM%IVAR_OFFSET2+1
                       FFM%OFFSET_MAP2(IB,M,L)=FFM%IVAR_OFFSET2
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        END SUBROUTINE SET_IVAR_OFFSET

        SUBROUTINE DEALLOCATE_OFFSET_MAP(FFM)
           USE ML_FF_STRUCT, ONLY: FFM_PAR
           IMPLICIT NONE
           ! Input and output variables
           TYPE (FFM_PAR)        :: FFM
           IF (ALLOCATED(FFM%OFFSET_MAP1)) DEALLOCATE(FFM%OFFSET_MAP1)
           IF (ALLOCATED(FFM%OFFSET_MAP2)) DEALLOCATE(FFM%OFFSET_MAP2)
        END SUBROUTINE DEALLOCATE_OFFSET_MAP

        LOGICAL FUNCTION CHECK_SF_LIMIT(FF, LOGF, ABN)

           USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR, ABINITIONEW
           USE LOGFILE, ONLY: LOGLINE_MSG, LOGLINE_SEPARATOR, MSG_WARNING
           USE ml_ff_string, ONLY: str
           USE ml_ff_tutor, ONLY: ml_tutor

           IMPLICIT NONE
           TYPE(FF_PAR), INTENT(IN)      :: FF
           TYPE(LOG_PAR), INTENT(IN)     :: LOGF !< Log file setup.
           TYPE(ABINITIONEW), INTENT(IN) :: ABN

           INTEGER :: NATOMS_OVER_SF_LIMIT

           ! If limit already exceeded, do not check again.
           IF (FF%SF_LIMIT_EXCEEDED) THEN
              CHECK_SF_LIMIT = .TRUE.
              RETURN
           END IF

           NATOMS_OVER_SF_LIMIT = COUNT(ABN%SF(1:, ABN%NCONF) > FF%SF_LIMIT)
           IF (NATOMS_OVER_SF_LIMIT > 0) THEN
              CALL ml_tutor%alert("Spilling factor limit (" // str(FF%SF_LIMIT) //                &
                                  ") was exceeded in ionic step " // str(FF%NSTEP) //             &
                                  " (" // str(NATOMS_OVER_SF_LIMIT) // " atom(s)). " //           &
                                  "This suggests that the machine-learned force field is strongly &
                                  &extrapolating and therefore results are unreliable. To avoid &
                                  &waste of compute resources the simulation will be gracefully &
                                  &exited after the next block of ionic steps.") 
              CALL LOGLINE_MSG(FF,                                        &
                               LOGF,                                      &
                               MSG_WARNING,                               &
                               "Spilling factor limit was exceeded  (" // &
                               str(NATOMS_OVER_SF_LIMIT) //               &
                               " atom(s)), triggering early exit of ionic loop.")
              CHECK_SF_LIMIT = .TRUE.
           ELSE
              CHECK_SF_LIMIT = .FALSE.
           END IF

        END FUNCTION CHECK_SF_LIMIT

      END MODULE FORCE_FIELD_ADDITIONAL

#endif
