#include "symbol.inc"
#define boundary_check
#define MPI_barrier_after_bcast
!***********************************************************************
!
!> this module implements high level routines to operate on wavefunctions
!>
!> care is taken to avoid indexing pointer arrays since
!> this incures performance penalties on all Intel compilers
!> (but particularly on the efc compiler)
!> (constructs such as W1%CW(I, ..) are avoided)
!> instead low level F77 routines are used, which are implemented
!> at the end of the routine
!
!***********************************************************************

MODULE wave_high
  USE prec
  USE wave
  USE wave_mpi
!***********************************************************************
!
!> interfaces for functions that accept only F77 style arrays
!> these F77 routines would cause a problem if the first element
!> of a pointer array is passed to the F77 routines
!
!***********************************************************************

  INTERFACE
     SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       GDEF, OPTIONAL :: CPROJ
     END SUBROUTINE ARRAY_TO_W1
  END INTERFACE

  INTERFACE
     SUBROUTINE W1_TO_ARRAY( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       GDEF, OPTIONAL :: CPROJ
     END SUBROUTINE W1_TO_ARRAY
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
GPACC !$ACC ROUTINE VECTOR
       USE prec
GPOMP !$OMP DECLARE TARGET(ECCP_NL)
       GDEF      CNL
       INTEGER LMDIM, LMMAXC
       OVERLAP CDIJ
       GDEF CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       LOGICAL LOVERL
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_SIMPLE(WDES1, LMDIM, CDIJ, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_SIMPLE
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       COMPLEX(q) EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_C
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_CCDIJ
  END INTERFACE

  INTERFACE ELEMENT
     MODULE PROCEDURE W1_FROM_W
     MODULE PROCEDURE W1_FROM_WA
     MODULE PROCEDURE W1_FROM_WA2
  END INTERFACE

  INTERFACE ELEMENTS
     MODULE PROCEDURE WA_FROM_W
     MODULE PROCEDURE WA_FROM_WA
     MODULE PROCEDURE WA_FROM_WA2
  END INTERFACE

  INTERFACE REDISTRIBUTE_PROJ
     MODULE PROCEDURE W_REDIS_PROJ
     MODULE PROCEDURE WA_REDIS_PROJ
     MODULE PROCEDURE W1_REDIS_PROJ
  END INTERFACE

  INTERFACE REDISTRIBUTE_PW
     MODULE PROCEDURE W_REDIS_PW
     MODULE PROCEDURE WA_REDIS_PW
     MODULE PROCEDURE W1_REDIS_PW
  END INTERFACE

  INTERFACE ASSIGNMENT
     MODULE PROCEDURE W1_COPY_REVERSE_ARG
  END INTERFACE

  INTERFACE CONNECT
     MODULE PROCEDURE WA_CONNECT2_W1
     MODULE PROCEDURE WA_CONNECT2_AR
     MODULE PROCEDURE W1_CONNECT2_AR
  END INTERFACE

  INTERFACE DISCONNECT
     MODULE PROCEDURE WA_NULLIFY
     MODULE PROCEDURE W1_NULLIFY
  END INTERFACE

  CONTAINS

!***********************************************************************
!
!> FFT of a wavefunction to real space
!
!***********************************************************************

    SUBROUTINE FFTWAV_W1( W1)
      IMPLICIT NONE
      INTEGER ISPINOR
      TYPE (wavefun1)    W1


      PROFILING_START('fftwav_w1')
      IF (W1%WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_W1: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)
      ENDIF
      IF (W1%WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(W1%WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(W1%WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTWAV_W1: FFTSCA is not associated and LUSEINV is &
               &set " // str(ASSOCIATED(W1%WDES1%FFTSCA)) // " " // str(ASSOCIATED(W1%WDES1%NINDPW_INV)), &
               __FILE__, __LINE__)
         ENDIF
         DO ISPINOR=0,W1%WDES1%NRSPINORS-1
            CALL FFTWAV_USEINV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), W1%WDES1%NINDPW_INV(1), W1%WDES1%FFTSCA(1,2), &
                 W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), &
                 W1%CW(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
         ENDDO
      ELSE
         DO ISPINOR=0,W1%WDES1%NRSPINORS-1
            CALL FFTWAV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), &
                 W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), &
                 W1%CW(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
         ENDDO
      ENDIF

      PROFILING_STOP('fftwav_w1')
    END SUBROUTINE FFTWAV_W1

!
!> FFT of a wavefunction to real space
!> similar version as above, however, receives WDES1 and explicit "WORK" arrays CR and CW
!> this routine is new in vasp.6 and should consistently replace
!> FFTWAV
!

    SUBROUTINE FFTWAV_WDES1( WDES1, CR, CW)
      IMPLICIT NONE
      !> descriptor for one k-point
      TYPE (wavedes1)    WDES1
      !> input:  orbitals in real space
      COMPLEX(q):: CR(WDES1%GRID%RL%NP)
      !> output: orbitals in reciprocal space
      COMPLEX(q):: CW(WDES1%NRPLWV)

      IF (WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_WDES1: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)
      ENDIF
      IF (WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTWAV_WDES1: FFTSCA is not associated and LUSEINV &
               &is set " // str(ASSOCIATED(WDES1%FFTSCA)) // " " // str(ASSOCIATED(WDES1%NINDPW_INV)), &
               __FILE__, __LINE__)
         ENDIF
         CALL FFTWAV_USEINV(WDES1%NGVECTOR, WDES1%NINDPW(1), WDES1%NINDPW_INV(1), WDES1%FFTSCA(1,2), &
              CR(1),CW(1),WDES1%GRID)
      ELSE
         CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
              CR(1),CW(1),WDES1%GRID)
      ENDIF
    END SUBROUTINE FFTWAV_WDES1


!***********************************************************************
!
!> FFT of a wavefunction to real space
!> a WDES1 must be supplied
!> as well as the input and output array
!> this routine is new in vasp.6 and should consistently replace
!> FFTEXT
!
!***********************************************************************

    SUBROUTINE FFTEXT_WDES1(WDES1, CR, CW, LADD)
      IMPLICIT NONE
      INTEGER ISPINOR
      !> descriptor for one k-point
      TYPE (wavedes1)    WDES1
      !> input:  orbitals in real space
      COMPLEX(q):: CR(WDES1%GRID%RL%NP)
      !> output: orbitals in reciprocal space
      COMPLEX(q):: CW(WDES1%NRPLWV)
      !> add results to CW and does not clear CW beforehand
      LOGICAL    :: LADD

      IF (WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_WDES1: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)
      ENDIF
      IF (WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTEXT_WDES1: FFTSCA is not associated and LUSEINV &
               &is set " // str(ASSOCIATED(WDES1%FFTSCA)) // " " // str(ASSOCIATED(WDES1%NINDPW_INV)), &
               __FILE__, __LINE__)
         ENDIF
         CALL FFTEXT_USEINV(WDES1%NGVECTOR, WDES1%NINDPW(1), WDES1%FFTSCA(1,1), &
              CR(1), CW(1), WDES1%GRID, LADD)
      ELSE
         CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), &
              CR(1), CW(1), WDES1%GRID, LADD)
      ENDIF
    END SUBROUTINE FFTEXT_WDES1


!***********************************************************************
!
!> This routine accomplishes the same thing as the assignment:
!> ~~~
!>   W1_LHS = W1_RHS
!> ~~~
!> would, but allows for the use of "fast"-memory (in the sense that
!> the fm_*arrays are excluded from the assignment).
!
!***********************************************************************

  SUBROUTINE W1_ASSIGN(W1_LHS,W1_RHS)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1_RHS
    TYPE (wavefun1) :: W1_LHS

    W1_LHS%WDES1  =>W1_RHS%WDES1
    W1_LHS%FERWE  = W1_RHS%FERWE
    W1_LHS%AUX    = W1_RHS%AUX
    W1_LHS%CELEN  = W1_RHS%CELEN
    W1_LHS%NB     = W1_RHS%NB
    W1_LHS%ISP    = W1_RHS%ISP
    W1_LHS%LDO    = W1_RHS%LDO

    W1_LHS%CPTWFP =>W1_RHS%CPTWFP
    W1_LHS%CPROJ  =>W1_RHS%CPROJ
    W1_LHS%CR     =>W1_RHS%CR

  END SUBROUTINE W1_ASSIGN


!***********************************************************************
!
!> copy a W1 structure
!>  ~~~
!>  W2 = W1  (W1 -> W2)
!>  ~~~
!> the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
!> for syntactic sugar the reverse operation is also present
!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_COPY( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2
    INTEGER I

    PROFILING_START('w1_copy')

#if defined(OFFLOADING) || ! defined(_OPENMP)
#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL WZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif

    IF (W1%WDES1%LGAMMA) THEN
       CALL WDCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1, W2%CPROJ(1), 1)
    ELSE
       CALL WZCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is not correct",__FILE__,__LINE__)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)),__FILE__,__LINE__)
       ENDIF
       CALL WZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    ENDIF
#else
!$OMP PARALLEL SHARED(W1,W2) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NRPLWV
       W2%CW(I)=W1%CW(I)
    ENDDO
!$OMP END DO
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NPROD
       W2%CPROJ(I)=W1%CPROJ(I)
    ENDDO
!$OMP END DO
    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is not correct",__FILE__,__LINE__)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)),__FILE__,__LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
          W2%CR(I)=W1%CR(I)
       ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL
#endif
    PROFILING_STOP('w1_copy')

  END SUBROUTINE W1_COPY

  SUBROUTINE W1_COPY_NOCR( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_NOCR

  SUBROUTINE W1_COPY_CPROJ( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_CPROJ

  SUBROUTINE W1_COPY_REVERSE_ARG( W2, W1 )
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)), __FILE__, __LINE__)
       ENDIF

       CALL ZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    ENDIF

  END SUBROUTINE W1_COPY_REVERSE_ARG


!***********************************************************************
!
!> copy a WA structure
!> ~~~
!>  W2 = W1  (W1 -> W2)
!> ~~~
!> the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
!> the redistributed wavefunctions pointers in the destination
!> are also properly set
!
!***********************************************************************

  SUBROUTINE WA_COPY( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CW) /= SIZE(W2%CW)) THEN
       CALL vtutor%bug("WA_COPY: size mismatch in CW " // str(SIZE(W1%CW)) // &
          " " // str(SIZE(W2%CW)), __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       CALL vtutor%bug("WA_COPY: size mismatch in CPROJ " // str(SIZE(W1%CPROJ)) &
          // " " // str(SIZE(W2%CPROJ)), __FILE__, __LINE__)
    ENDIF

    CALL ZCOPY( SIZE(W1%CW), W1%CW(1,1), 1, W2%CW(1,1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ENDIF
    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CW(1,1))
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CW_RED=>W2%CPTWFP
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF
  END SUBROUTINE WA_COPY

  SUBROUTINE WA_COPY_CW( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CW) /= SIZE(W2%CW)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CW " // str(SIZE(W1%CW)) // &
          " " // str(SIZE(W2%CW)), __FILE__, __LINE__)
    ENDIF

    CALL ZCOPY( SIZE(W1%CW), W1%CW(1,1), 1, W2%CW(1,1), 1)

    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CW(1,1))
    ELSE
       W2%CW_RED=>W2%CPTWFP
    ENDIF
  END SUBROUTINE WA_COPY_CW

  SUBROUTINE WA_COPY_CPROJ( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CPROJ " // str(SIZE(W1%CPROJ)) &
          // " " // str(SIZE(W2%CPROJ)), __FILE__, __LINE__)
    ENDIF

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
      ENDIF
    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF
  END SUBROUTINE WA_COPY_CPROJ


!***********************************************************************
!
!> redistribute the wavefunction character for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PROJ( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PROJ(W1%WDES1, 1, W1%CPROJ(1))
  END SUBROUTINE W1_REDIS_PROJ

  SUBROUTINE WA_REDIS_PROJ( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PROJ(WA%WDES1, SIZE(WA%CPROJ,2), WA%CPROJ(1,1))
  END SUBROUTINE WA_REDIS_PROJ

  SUBROUTINE W_REDIS_PROJ( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS
#ifdef MPI
          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PROJ( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PROJ


!***********************************************************************
!
!> redistribute the plane wave coefficients for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PW( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PW(W1%WDES1, 1, W1%CW(1))
  END SUBROUTINE W1_REDIS_PW

  SUBROUTINE WA_REDIS_PW( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PW(WA%WDES1, SIZE(WA%CW,2), WA%CW(1,1))
  END SUBROUTINE WA_REDIS_PW


  SUBROUTINE W_REDIS_PW( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS
#ifdef MPI
          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PW( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PW


!***********************************************************************
!
!> update of vector
!>  ~~~
!>  W2 = W1*a + W2
!>  ~~~
!>
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_DAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1, W2
!$  INTEGER I

    PROFILING_START('w1_daxpy')

#if defined(OFFLOADING) || ! defined( _OPENMP)
    CALL WDAXPY( W1%WDES1%NPL*2, SCALE, W1%CW(1), 1, W2%CW(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL WDAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL WDAXPY( W1%WDES1%NPRO*2, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF

       ! real space wavefunction complex, SCALE real
       CALL WDAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
    ENDIF
#else
!$OMP PARALLEL SHARED(W1,W2,SCALE) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NPL
       W2%CW(I)=W2%CW(I)+SCALE*W1%CW(I)
    ENDDO
!$OMP END DO NOWAIT
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NPRO
       W2%CPROJ(I)=W2%CPROJ(I)+SCALE*W1%CPROJ(I)
    ENDDO
!$OMP END DO NOWAIT
    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
       W2%CR(I)=W2%CR(I)+SCALE*W1%CR(I)
    ENDDO
!$OMP END DO NOWAIT
    ENDIF
!$OMP END PARALLEL
#endif

    PROFILING_STOP('w1_daxpy')

  END SUBROUTINE W1_DAXPY

!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
  SUBROUTINE W1_GAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    GDEF SCALE
    TYPE (wavefun1)    W1, W2
!$  INTEGER I

    PROFILING_START('w1_gaxpy')

#if defined(OFFLOADING) || ! defined( _OPENMP)
    IF (W1%WDES1%LGAMMA) THEN
       ! wavefunction complex, SCALE real
       CALL WDAXPY( W1%WDES1%NPL*2, SCALE, W1%CW(1), 1, W2%CW(1), 1)
       CALL WDAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       ! wavefunction complex, SCALE complex
       CALL WZAXPY( W1%WDES1%NPL, SCALE, W1%CW(1), 1, W2%CW(1), 1)
       CALL WZAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF

       IF (W1%WDES1%LGAMMA) THEN
          ! real space wavefunction complex, SCALE real
          CALL WDAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ELSE
          ! real space wavefunction complex, SCALE complex
          CALL WZAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ENDIF
    ENDIF
#else
!$OMP PARALLEL SHARED(W1,W2,SCALE) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%NPL
          W2%CW(I)=W2%CW(I)+SCALE*W1%CW(I)
       ENDDO
!$OMP END DO NOWAIT
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%NPRO
          W2%CPROJ(I)=W2%CPROJ(I)+SCALE*W1%CPROJ(I)
       ENDDO
!$OMP END DO NOWAIT
    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
          W2%CR(I)=W2%CR(I)+SCALE*W1%CR(I)
       ENDDO
!$OMP END DO NOWAIT
    ENDIF
!$OMP END PARALLEL
#endif

    PROFILING_STOP('w1_gaxpy')

  END SUBROUTINE W1_GAXPY


!***********************************************************************
!
!> update of vector
!>  ~~~
!>  W1 = W1*a
!>  ~~~
!>
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_DSCAL( W1, SCALE)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1
!$  INTEGER I

    PROFILING_START('w1_dscal')

    ! since this function can be used to zero out an array we operate
    ! on all elements (dimension) and not only on those that are
    ! actually used
#if defined(OFFLOADING) || ! defined(_OPENMP)
    CALL WDSCAL( W1%WDES1%NRPLWV*2, SCALE, W1%CW(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL WDSCAL( W1%WDES1%NPROD, SCALE, W1%CPROJ(1), 1)
    ELSE
       CALL WDSCAL( W1%WDES1%NPROD*2, SCALE, W1%CPROJ(1),  1)
    ENDIF

    IF (ASSOCIATED(W1%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
       ! real space wavefunction complex, SCALE real
       CALL WDSCAL( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1)
    ENDIF
#else
!$OMP PARALLEL SHARED(W1,SCALE) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
     DO I=1,W1%WDES1%NRPLWV
        W1%CW(I)=SCALE*W1%CW(I)
     ENDDO
!$OMP END DO NOWAIT
!DIR$ IVDEP
!$OMP DO
     DO I=1,W1%WDES1%NPROD
        W1%CPROJ(I)=SCALE*W1%CPROJ(I)
     ENDDO
!$OMP END DO NOWAIT
    IF (ASSOCIATED(W1%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
     DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
        W1%CR(I)=SCALE*W1%CR(I)
     ENDDO
!$OMP END DO NOWAIT
    ENDIF
!$OMP END PARALLEL
#endif

    PROFILING_STOP('w1_dscal')

  END SUBROUTINE W1_DSCAL


!***********************************************************************
!
!> calculate the dot product between two wavefunctions
!> ~~~
!>  C=   W1^* x W2
!> ~~~
!> this is a substitue for the routine CINDPROD but mind
!> the W1 and W2 are interchanged
!>
!> @details @ref openmp :
!> a ZDOTC call is replaced by an explicit OMP PARALLEL loop
!> (for some reason BLAS1 call thread really badly).
!> The nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\", that is distributed over all available
!> threads.
!
!***********************************************************************

  FUNCTION W1_DOT( W1, W2, CQIJ) RESULT (C)

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefun1)   :: W1, W2
    GDEF              :: C
    OVERLAP, OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
  ! local
    GDEF    :: CNL, CTMP
    INTEGER :: LMDIM, NPRO, NPRO_, ISPINOR, ISPINOR_, LMMAXC, NT, NI

    REAL(q),    EXTERNAL :: WDDOT
    COMPLEX(q), EXTERNAL :: WZDOTC

!$  INTEGER :: I
!$ACC ROUTINE(ECCP_NL) VECTOR

    PROFILING_START('w1_dot')

#if defined(OFFLOADING) || ! defined(_OPENMP)
    IF (W1%WDES1%LGAMMA) THEN
       C=WDDOT(2*W1%WDES1%NPL,W1%CW(1),1,W2%CW(1),1)
    ELSE
       C=WZDOTC( W1%WDES1%NPL,W1%CW(1),1,W2%CW(1),1)
    ENDIF
#else
    C=0
!$OMP PARALLEL DO PRIVATE(I) SHARED(W1,W2) REDUCTION(+:C)
    DO I=1,W1%WDES1%NPL
       C=C+CONJG(W1%CW(I))*W2%CW(I)
    ENDDO
!$OMP END PARALLEL DO
#endif
    IF (PRESENT(CQIJ) .AND. W1%WDES1%LOVERL .AND.  W1%WDES1%NPROD>0 ) THEN
       CNL=0
       LMDIM=SIZE(CQIJ,1)

GPACC !$ACC ENTER DATA COPYIN(W1,W2) __IF_ASYNC__
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(CQIJ,W1,W1%WDES1,W1%CPROJ,W2,W2%CPROJ) REDUCTION(+:CNL) &
GPACC !$ACC PRIVATE(NT,LMMAXC,NPRO,NPRO_,CTMP) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:CQIJ,W1,W1%WDES1,W1%CPROJ,W2,W2%CPROJ) &
GPOMP !$OMP PRIVATE(NT,LMMAXC,NPRO,NPRO_,CTMP) REDUCTION(+:CNL) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP SHARED(LMDIM,W1,W2,CQIJ) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,CTMP) &
__CPU !$OMP REDUCTION(+:CNL)

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
       DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
          DO NI=1,W1%WDES1%NIONS
             NT=W1%WDES1%ITYP(NI)
             LMMAXC=W1%WDES1%LMMAX(NT)
             IF (LMMAXC==0) CYCLE
             NPRO =ISPINOR *(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)

             CTMP=0; CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W2%CPROJ(NPRO_+1),W1%CPROJ(NPRO+1),CTMP)
             CNL=CNL+CTMP
          ENDDO
       ENDDO
       ENDDO spinor

GPACC !$ACC EXIT DATA DELETE(W1,W2) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

__CPU !$OMP END PARALLEL DO

       C=C+CNL
    ENDIF
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    CALLMPI( M_sum_g(W1%WDES1%COMM_INB, C, 1))

    PROFILING_STOP('w1_dot')

  END FUNCTION W1_DOT


!***********************************************************************
!
!> calculate the inproduct between one wavefunction and
!> a set of wavefunctions with possible scales and add the result
!> to a third vector
!>  ~~~
!>  C=   WA^* x W1*SCALEA+C * SCALEC
!>  ~~~
!
!***********************************************************************

  SUBROUTINE W1_GEMV(WA, W1, C, CQIJ)

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefuna)    WA
    TYPE (wavefun1)    W1
    GDEF:: C(*)
    !> optional overlap operator
    OVERLAP, OPTIONAL:: CQIJ(:,:,:,:)
  ! local
    INTEGER :: LMDIM, N
    GDEF CRESUL(WA%WDES1%NPRO)

    PROFILING_START('w1_gemv')

    N=SIZE(WA%CW,2)

    IF (W1%WDES1%LGAMMA) THEN
       CALL WDGEMV( dtrans, 2* WA%WDES1%NPL, N, done , WA%CW(1,1) , &
                        2* WA%WDES1%NRPLWV, W1%CW(1) , 1 , dzero,  C(1), 1)
    ELSE
       CALL WZGEMV( ztrans,  WA%WDES1%NPL, N, zone , WA%CW(1,1) , &
                         WA%WDES1%NRPLWV, W1%CW(1) , 1 , zzero,  C(1), 1)
    ENDIF

    IF (PRESENT(CQIJ) .AND. WA%WDES1%LOVERL .AND.  WA%WDES1%NPROD>0 ) THEN
       LMDIM = SIZE(CQIJ,1)
GPACC !$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CRESUL) IF(OFFLOAD_ON)
       CALL OVERL1_SIMPLE(WA%WDES1, LMDIM, CQIJ(1,1,1,1), W1%CPROJ(1), CRESUL(1))
       IF (W1%WDES1%LGAMMA) THEN
          CALL WDGEMV( dtrans, WA%WDES1%NPRO,  N, done ,WA%CPROJ(1,1), &
                           WA%WDES1%NPROD, CRESUL , 1 , done,  C(1), 1)
       ELSE
          CALL WZGEMV( ztrans, WA%WDES1%NPRO,  N, zone ,WA%CPROJ(1,1), &
                           WA%WDES1%NPROD, CRESUL , 1 , zone,  C(1), 1)
       ENDIF
GPACC !$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CRESUL) IF(OFFLOAD_ON)
    ENDIF

    CALLMPI( M_sum_g(WA%WDES1%COMM_INB, C, N))

    PROFILING_STOP('w1_gemv')

  END SUBROUTINE W1_GEMV


!***********************************************************************
!
!> W1 descriptor from W array
!>
!> alternative to SETWAV (returns a W1 descriptor)
!> it is somewhat slimmed down to optimize performance
!> and returns only the wavefunction and wavefunction character pointers
!> the full version remains SETWAV
!
!***********************************************************************

  FUNCTION W1_FROM_W( W, WDES1, NB, ISP) RESULT (W1)
    IMPLICIT NONE
    INTEGER NB, ISP
    TYPE (wavespin) W
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    NK=WDES1%NK
#ifdef boundary_check
    IF (NB<=0 .OR. NB> SIZE(W%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_W: bounds exceed " // str(NB) // " " // &
          str(SIZE(W%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
    W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
    W1%WDES1 => WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_W


!***********************************************************************
!
!> W1 descriptor from WA descriptor
!
!***********************************************************************

  FUNCTION W1_FROM_WA( WA, N1) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1)
    W1%CPROJ =>WA%CPROJ(:,N1)
    W1%WDES1 =>WA%WDES1

    IF (ASSOCIATED(WA%CR)) THEN
       W1%CR =>WA%CR(:,N1)
    ELSE
       NULLIFY(W1%CR)
    ENDIF

    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA


!***********************************************************************
!
!> subindex a WA array
!> return WA(N1:N2)
!
!***********************************************************************

  FUNCTION WA_FROM_WA( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in WA_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
    IF (N2<N1 .OR. N2> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in WA_FROM_WA: bounds exceed " // str(N1) // " " // str(N2) &
          // " " // str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1:N2)
    W1%CPROJ =>WA%CPROJ(:,N1:N2)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA


!***********************************************************************
!
!> one element from a WA array using two indices
!> return WA(N1, N2)
!
!***********************************************************************

  FUNCTION W1_FROM_WA2( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA2


!***********************************************************************
!
!> subindex a WA array using two indices
!> return WA(N1:N12, N2)
!
!***********************************************************************

  FUNCTION WA_FROM_WA2( WA, N1, N12, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N12, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
    IF (N12<N1 .OR. N12> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // str(N12) &
          // " " // str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N12+ (N2-1)*WA%FIRST_DIM<N1+ (N2-1)*WA%FIRST_DIM .OR. N12+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N12+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF

#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA2


!***********************************************************************
!
!> get a WA structure from a wavefunction array
!> this is equivalent to the low level routine SETWAVA
!
!***********************************************************************

  FUNCTION WA_FROM_W( W, WDES1, ISP ) RESULT (WA)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1), TARGET ::  WDES1
    INTEGER ISP
    TYPE (wavefuna) ::  WA

    WA%CPTWFP=>W%CPTWFP(:,:,WDES1%NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,WDES1%NK,ISP)
    WA%FERWE =>W%FERWE(:,WDES1%NK,ISP)
    WA%AUX   =>W%AUX  (:,WDES1%NK,ISP)
    WA%CELEN =>W%CELEN(:,WDES1%NK,ISP)
    WA%WDES1 =>WDES1
    ! remember spin index
    WA%ISP =ISP
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
  END FUNCTION WA_FROM_W


!************************* SUBROUTINE NEWWAVA **************************
!
!>  create storage for a wavefunction array WA
!
!***********************************************************************

  SUBROUTINE NEWWAVA(WA, WDES1, NDIM, NDIM2, ALLOC_REAL)

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefuna) WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER :: NDIM
    INTEGER, OPTIONAL :: NDIM2
    LOGICAL, OPTIONAL :: ALLOC_REAL
    ! local variables
    INTEGER :: NDIMTOT

    NDIMTOT=NDIM
    WA%FIRST_DIM=0

    IF (PRESENT(NDIM2)) THEN
       NDIMTOT=NDIM*NDIM2
       WA%FIRST_DIM=NDIM
    ENDIF

    ALLOCATE(WA%CPTWFP(WDES1%NRPLWV,NDIMTOT), WA%CPROJ(WDES1%NPROD,NDIMTOT), &
             WA%CELEN(NDIMTOT), WA%FERWE(NDIMTOT), WA%AUX(NDIMTOT))
GPACC !$ACC ENTER DATA CREATE(WA%CPTWFP,WA%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
! Always, alloc needs to be set for openMP port if the pointer member is already on the device and the dope vector is updated here
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: WA%CPTWFP,WA%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))

DOESI WA%CPTWFP=0; WA%CPROJ=0 ! Breaks ???.

    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, NDIM*WDES1%NB_PAR, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED   =>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
GPACC !$ACC ENTER DATA CREATE(WA%CW_RED,WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(to: WA, WA%CW_RED, WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: WA%CW_RED, WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))

    ! optionally allocate space for the orbitals in real space
    IF (PRESENT(ALLOC_REAL)) THEN
       IF (ALLOC_REAL) THEN
          ALLOCATE(WA%CR(WDES1%GRID%MPLWV*WDES1%NRSPINORS,NDIMTOT))
GPACC !$ACC ENTER DATA CREATE(WA%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
! Always, alloc needs to be set for openMP port if the pointer member is already on the device and the dope vector is updated here
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: WA%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
       ENDIF
    ENDIF

    WA%ISP=-1

    WA%WDES1=>WDES1
GPACC !$ACC UPDATE DEVICE(WA%FIRST_DIM,WA%ISP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
!!      CALL COPYIN_TYPED_VAR(WA%WDES1)
!!      CALL ACC_ATTACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(TO: WA, WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
! set also the pointer
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
GPOMP !$OMP TARGET UPDATE TO(WA%FIRST_DIM,WA%ISP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))

  END SUBROUTINE NEWWAVA


!************************* SUBROUTINE NEWWAVA_PROJ *********************
!
!>  create storage for one wavefunction W array
!>  to store non local part only
!
!***********************************************************************

  SUBROUTINE NEWWAVA_PROJ(WA, WDES1, NDIM)

__GPU USE moffload

    USE prec
    USE wave_mpi
    IMPLICIT NONE
    TYPE (wavefuna) WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER, OPTIONAL :: NDIM
    IF (PRESENT(NDIM)) THEN
       ALLOCATE(WA%CPROJ(WDES1%NPROD,NDIM))
    ELSE
       ALLOCATE(WA%CPROJ(WDES1%NPROD,WDES1%NBANDS))
    ENDIF
GPACC !$ACC ENTER DATA CREATE(WA%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: WA%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
    WA%FIRST_DIM=0
    WA%ISP=-1
    NULLIFY(WA%CPTWFP)
    NULLIFY(WA%CW_RED)
    NULLIFY(WA%FERWE)
    NULLIFY(WA%AUX  )
    NULLIFY(WA%CELEN)
    ! set redistributed projection indices
    IF (WDES1%DO_REDIS) THEN
       IF (PRESENT(NDIM)) THEN
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
       ELSE
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  WDES1%NB_TOT, WA%CPROJ(1,1))
       ENDIF
    ELSE
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
GPACC !$ACC ENTER DATA CREATE(WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
    WA%WDES1=>WDES1
GPACC !$ACC UPDATE DEVICE(WA%FIRST_DIM,WA%ISP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
!!      CALL COPYIN_TYPED_VAR(WA%WDES1)
!!      CALL ACC_ATTACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE TO(WA%FIRST_DIM,WA%ISP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
GPOMP !$OMP TARGET ENTER DATA MAP(TO:WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA))
  END SUBROUTINE NEWWAVA_PROJ


!************************* SUBROUTINE NEWWAV_ARRAY *********************
!
!>  create storage for an array of single orbital structures W1(:)
!
!***********************************************************************

  SUBROUTINE NEWWAV_ARRAY(W1,WDES1,ALLOC_REAL,WA)

GPACC USE moffload

    IMPLICIT NONE

    TYPE (wavefun1) :: W1(:)
    TYPE (wavedes1) :: WDES1
    TYPE (wavefuna), OPTIONAL :: WA

    LOGICAL :: ALLOC_REAL

    ! local variables
    INTEGER :: NDIM,I

    NDIM=SIZE(W1)

    IF (PRESENT(WA)) THEN
!
! Allocate a wavefuna structure WA to create storage space for the pointer
! members of the array of single orbital structures W1(:).
! This guarantees that the respective orbital data arrays W1(:)%CPTWFP and
! W1(:)%CR will be laid out contiguously in memory.
!
!$ACC ENTER DATA CREATE(WA) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(:),1)) ASYNC(ACC_ASYNC_Q)
       ! allocate a wavefuna structure WA of size NDIM
       CALL NEWWAVA(WA, WDES1, NDIM, ALLOC_REAL=ALLOC_REAL)
       ! connect the W1 members to WA
       DO I=1,NDIM
          ! connect pointers
          W1(I) = ELEMENT(WA, I)
!$ACC ENTER DATA COPYIN(W1(I)%CPTWFP,W1(I)%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(W1(I)%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1).AND.ALLOC_REAL) ASYNC(ACC_ASYNC_Q)
!!!$ACC ENTER DATA CREATE(W1(I)%CPTWFP,W1(I)%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1)) ASYNC(ACC_ASYNC_Q)
!!!$ACC ENTER DATA CREATE(W1(I)%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1).AND.ALLOC_REAL) ASYNC(ACC_ASYNC_Q)

          ! set other members of W1 (see NEWWAV)
          W1(I)%FERWE= 0
          W1(I)%CELEN= 0
          W1(I)%NB   =-1
          W1(I)%ISP  =-1
          W1(I)%LDO  = .TRUE.

#ifdef ACC_OFFLOAD
          IF (OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1)) THEN
!$ACC UPDATE DEVICE(W1(I)%FERWE,W1(I)%CELEN,W1(I)%NB,W1(I)%ISP,W1(I)%LDO) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(W1(I)%WDES1) ASYNC(ACC_ASYNC_Q)
          ENDIF
#endif
       ENDDO
    ELSE
!
! In case WA is not present we fall back to the use of NEWWAV.
!
       DO I=1,NDIM
          CALL NEWWAV(W1(I), WDES1, ALLOC_REAL)
       ENDDO
    ENDIF
  END SUBROUTINE NEWWAV_ARRAY


!************************* SUBROUTINE DELWAV_ARRAY *********************
!
!>  destroy storage for an array of single orbital structures W1(:)
!>  that was possibly created using a wavefuna structure
!
!***********************************************************************

  SUBROUTINE DELWAV_ARRAY(W1,DEALLOC_REAL,WA)

GPACC USE moffload

    IMPLICIT NONE

    TYPE (wavefun1) :: W1(:)
    TYPE (wavefuna), OPTIONAL :: WA

    LOGICAL :: DEALLOC_REAL

    ! local variables
    INTEGER :: NDIM,I

    NDIM=SIZE(W1)

    IF (PRESENT(WA)) THEN
       IF (NDIM/=SIZE(WA%CPTWFP,2)) CALL vtutor%bug("DELWAV_ARRAY: partial deallocation is not supported",__FILE__,__LINE__)
       ! deallocate WA structure
       CALL DELWAVA(WA)
!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
       ! nullify W1 pointer members
       DO I=1,NDIM
!$ACC EXIT DATA DELETE(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1)) ASYNC(ACC_ASYNC_Q)
          NULLIFY(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1)
          IF (DEALLOC_REAL) THEN
!$ACC EXIT DATA DELETE(W1(I)%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1)) ASYNC(ACC_ASYNC_Q)
             NULLIFY(W1(I)%CR)
          ENDIF
       ENDDO
    ELSE
       DO I=1,NDIM
          CALL DELWAV(W1(I), DEALLOC_REAL)
       ENDDO
    ENDIF
  END SUBROUTINE DELWAV_ARRAY


!************************* FUNCTION WA_CONNECT2_W1 *********************
!
!>  Connect the pointer members of a wavefuna structure to the members
!>  of an array of single orbital structures
!>
!>  Mind: the latter has to have been created using NEWWAV_ARRAY, to
!>  ensure that things like W1(:)%CPTWFP are contiguous in memory
!>
!>  WA_CONNECT2_W1 will try to make sure that the members of WA and W1(:)
!>  are properly mapped
!
!***********************************************************************

  SUBROUTINE WA_CONNECT2_W1(W1,WA)

GPACC USE moffload

    USE iso_c_binding
    IMPLICIT NONE
    TYPE (wavefun1) :: W1(:)
    TYPE (wavefuna) :: WA

    ! local variables
    INTEGER :: I

    CALL c_f_pointer(c_loc(W1(1)%CPTWFP(1)),WA%CPTWFP,[SIZE(W1(1)%CPTWFP),SIZE(W1)])
    CALL c_f_pointer(c_loc(W1(1)%CPROJ (1)),WA%CPROJ ,[SIZE(W1(1)%CPROJ ),SIZE(W1)])

    IF (ASSOCIATED(W1(1)%CR)) CALL c_f_pointer(c_loc(W1(1)%CR(1)),WA%CR,[SIZE(W1(1)%CR),SIZE(W1)])

    WA%WDES1 => W1(1)%WDES1

    ! test consistency
    DO I=1,SIZE(W1)
       IF (.NOT.c_associated(c_loc(WA%CPTWFP(1,I)),c_loc(W1(I)%CPTWFP(1)))) &
          CALL vtutor%bug("WA_CONNECT2_W1: target member CPTWFP appears to be non-contiguous",__FILE__,__LINE__)
       IF (.NOT.c_associated(c_loc(WA%CPROJ(1,I)),c_loc(W1(I)%CPROJ(1)))) &
          CALL vtutor%bug("WA_CONNECT2_W1: target member CPROJ appears to be non-contiguous", __FILE__,__LINE__)
        IF (ASSOCIATED(WA%CR) .AND. .NOT.c_associated(c_loc(WA%CR(1,I)),c_loc(W1(I)%CR(1)))) &
          CALL vtutor%bug("WA_CONNECT2_W1: target member CR appears to be non-contiguous",    __FILE__,__LINE__)
    ENDDO

!$ACC ENTER DATA COPYIN(WA) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1).AND.ASSOCIATED(WA%CR)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE WA_CONNECT2_W1


!************************* FUNCTION WA_CONNECT2_AR *********************
!
!>  Allocate a wavefun1 array and connect its pointer members to
!>  two 2d-arrays, CW and CPROJ, and a wavedes1 structure.
!>  Mind the leading dimensions of CW and CPROJ must be larger than
!>  or equal to WDES1%NRPLWV and WDES1%NPROD, respectively.
!>  The second dimensions of CW and CPROJ must be equal.
!
!***********************************************************************

  SUBROUTINE WA_CONNECT2_AR(CW,CPROJ,WDES1,WA)

GPACC USE moffload

    IMPLICIT NONE
    COMPLEX(q), TARGET :: CW(:,:)
    GDEF,       TARGET :: CPROJ(:,:)

    TYPE (wavedes1), TARGET :: WDES1
    TYPE (wavefuna) :: WA

    ! local
    INTEGER :: NDIM2, I

#ifdef ACC_OFFLOAD
    LOGICAL :: ACC_ACTIVE
    ACC_ACTIVE=OFFLOAD_ON.AND.IS_ON_DEVICE(CW,1).AND.IS_ON_DEVICE(CPROJ,1)
#endif

    IF (SIZE(CW,2)/=SIZE(CPROJ,2)) &
       CALL vtutor%bug("WA_CONNECT2_AR: array sizes inconsistent",__FILE__,__LINE__)

    IF (SIZE(CW,1)<WDES1%NRPLWV) &
       CALL vtutor%bug("WA_CONNECT2_AR: first dimension of CW too small "// &
          str(SIZE(CW,1))//" "//str(WDES1%NRPLWV),__FILE__,__LINE__)

    IF (SIZE(CPROJ,1)<WDES1%NPROD) &
       CALL vtutor%bug("WA_CONNECT2_AR: first dimension of CPROJ too small "// &
          str(SIZE(CPROJ,1))//" "//str(WDES1%NPROD),__FILE__,__LINE__)

    WA%CPTWFP=>CW
    WA%CPROJ =>CPROJ
    WA%WDES1 =>WDES1

!$ACC ENTER DATA CREATE(WA) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE WA_CONNECT2_AR


!************************* SUBROUTINE WA_NULLIFY ***********************
!
!>  nullify the pointer members of a wavefuna structure
!
!***********************************************************************

  SUBROUTINE WA_NULLIFY(WA)

GPACC USE moffload

    IMPLICIT NONE
    TYPE (wavefuna) :: WA
!$ACC EXIT DATA DELETE(WA%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1).AND.ASSOCIATED(WA%CR)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1))  ASYNC(ACC_ASYNC_Q)
    NULLIFY(WA%CPTWFP,WA%CPROJ,WA%WDES1); IF (ASSOCIATED(WA%CR)) NULLIFY(WA%CR)
  END SUBROUTINE WA_NULLIFY


!************************* FUNCTION W1_CONNECT2_AR *********************
!
!>  Allocate a wavefun1 array and connect its pointer members to
!>  two 2d-arrays, CW and CPROJ, and a wavedes1 structure.
!>
!>  Mind the leading dimensions of CW and CPROJ must be larger than
!>  or equal to WDES1%NRPLWV and WDES1%NPROD, respectively.
!>  The second dimensions of CW and CPROJ must be equal.
!
!***********************************************************************

  SUBROUTINE W1_CONNECT2_AR(CW,CPROJ,WDES1,W1)

GPACC USE moffload

    IMPLICIT NONE
    COMPLEX(q), TARGET :: CW(:,:)
    GDEF,       TARGET :: CPROJ(:,:)

    TYPE (wavedes1), TARGET  :: WDES1
    TYPE (wavefun1), POINTER :: W1(:)

    ! local
    INTEGER :: NDIM2, I

#ifdef ACC_OFFLOAD
    LOGICAL :: ACC_ACTIVE
    ACC_ACTIVE=OFFLOAD_ON.AND.IS_ON_DEVICE(CW,1).AND.IS_ON_DEVICE(CPROJ,1)
#endif

    IF (SIZE(CW,2)/=SIZE(CPROJ,2)) &
       CALL vtutor%bug("W1_CONNECT2_AR: array sizes inconsistent",__FILE__,__LINE__)

    IF (SIZE(CW,1)<WDES1%NRPLWV) &
       CALL vtutor%bug("W1_CONNECT2_AR: first dimension of CW too small "// &
          str(SIZE(CW,1))//" "//str(WDES1%NRPLWV),__FILE__,__LINE__)

    IF (SIZE(CPROJ,1)<WDES1%NPROD) &
       CALL vtutor%bug("W1_CONNECT2_AR: first dimension of CPROJ too small "// &
          str(SIZE(CPROJ,1))//" "//str(WDES1%NPROD),__FILE__,__LINE__)

    NDIM2=SIZE(CW,2); ALLOCATE(W1(NDIM2))
!$ACC ENTER DATA CREATE(W1(:)) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
    DO I=1,NDIM2
       W1(I)%CPTWFP => CW(:,I)
       W1(I)%CPROJ  => CPROJ(:,I)
       W1(I)%WDES1  => WDES1
       W1(I)%LDO    =  .TRUE.
!$ACC ENTER DATA COPYIN(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!$ACC UPDATE DEVICE(W1(I)%LDO) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
    ENDDO
  END SUBROUTINE W1_CONNECT2_AR


!************************* SUBROUTINE W1_NULLIFY ***********************
!
!>  nullify the pointer members of a wavefun1 array
!
!***********************************************************************

  SUBROUTINE W1_NULLIFY(W1)

GPACC USE moffload

    IMPLICIT NONE
    TYPE (wavefun1), POINTER :: W1(:)
    ! local
    INTEGER :: I
    DO I=1,SIZE(W1)
!$ACC EXIT DATA DELETE(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1(I),1)) ASYNC(ACC_ASYNC_Q)
       NULLIFY(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1)
    ENDDO
!$ACC EXIT DATA DELETE(W1(:)) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(W1)
  END SUBROUTINE W1_NULLIFY


!************************* SUBROUTINE FFTWAV_WA ************************
!
!>  batched FFT (q -> r) of a group of orbitals
!
!***********************************************************************

  SUBROUTINE FFTWAV_WA(WA)
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
    ! local variables
    INTEGER :: ISPINOR

    IF (WA%WDES1%NK==0) &
       CALL vtutor%bug("FFTWAV_WA: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)

    IF (.NOT.ASSOCIATED(WA%CR)) &
       CALL vtutor%bug("FFTWAV_WA: can not perform FFT: CR is not allocated", __FILE__, __LINE__)

    IF (WA%WDES1%LUSEINV) THEN
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTWAV_USEINV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1), &
                                WA%WDES1%NINDPW_INV(1),WA%WDES1%FFTSCA(1,2), &
                                WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                                WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                                WA%WDES1%GRID)
       ENDDO
    ELSE
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTWAV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1), &
                         WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                         WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                         WA%WDES1%GRID)
       ENDDO
    ENDIF
  END SUBROUTINE FFTWAV_WA


!************************* SUBROUTINE FFTEXT_WA ************************
!
!>  batched FFT (r -> q) of a group of orbitals
!
!***********************************************************************

  SUBROUTINE FFTEXT_WA(WA,LADD)
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
    LOGICAL :: LADD
    ! local variables
    INTEGER :: ISPINOR

    IF (WA%WDES1%NK==0) &
       CALL vtutor%bug("FFTEXT_WA: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)

    IF (.NOT.ASSOCIATED(WA%CR)) &
       CALL vtutor%bug("FFTEXT_WA: can not perform FFT: CR is not allocated", __FILE__, __LINE__)

    IF (WA%WDES1%LUSEINV) THEN
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTEXT_USEINV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1),WA%WDES1%FFTSCA(1,1), &
                                WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                                WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                                WA%WDES1%GRID,LADD)
       ENDDO
    ELSE
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTEXT_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1), &
                         WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                         WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                         WA%WDES1%GRID,LADD)
       ENDDO
    ENDIF
  END SUBROUTINE FFTEXT_WA


!************************* SUBROUTINE SETWAVA **************************
!
!>  set one single wavefunction array (WA) from an array of wavefunctions
!
!***********************************************************************

  SUBROUTINE SETWAVA(W, WA, WDES1, ISP)
    USE prec
    IMPLICIT NONE
    INTEGER NB,ISP
    TYPE (wavespin) W
    TYPE (wavefuna) WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    NK=WDES1%NK

    WA%CPTWFP=>W%CPTWFP(:,:,NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,NK,ISP)
    WA%FERWE =>W%FERWE(:,NK,ISP)
    WA%AUX   =>W%AUX  (:,NK,ISP)
    WA%CELEN =>W%CELEN(:,NK,ISP)
    WA%WDES1 =>WDES1
    ! one dimensional indexing assumed
    WA%FIRST_DIM=0
    ! remember spin index
    WA%ISP =ISP
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WA%WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WA%WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF

  END SUBROUTINE SETWAVA


!************************* SUBROUTINE DELWAVA **************************
!
!>  destroy storage for a wavefunctionarray WA
!
!***********************************************************************

      SUBROUTINE DELWAVA(WA)

__GPU USE moffload

      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (.NOT. ASSOCIATED(WA%CPTWFP) .OR. .NOT. ASSOCIATED(WA%FERWE) .OR. &
          .NOT. ASSOCIATED(WA%CELEN)  .OR. .NOT. ASSOCIATED(WA%AUX)) THEN
         CALL vtutor%bug("internal error in DELWAVA: not all enities are associated, try DELWAVA_PROJ", __FILE__, __LINE__)
      ENDIF

GPACC !$ACC EXIT DATA DELETE(WA%CPTWFP,WA%CPROJ,WA%CW_RED,WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%CW_RED, WA%CPROJ_RED) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA, WA%CW_RED, WA%CPROJ_RED) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%CPTWFP,WA%CPROJ) IF(OFFLOAD_ON)
      DEALLOCATE(WA%CPTWFP,WA%CPROJ,WA%CELEN, WA%FERWE, WA%AUX)
      NULLIFY(WA%CW_RED,WA%CPROJ_RED)

      IF (ASSOCIATED(WA%CR)) THEN
GPACC !$ACC EXIT DATA DELETE(WA%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%CR) IF(OFFLOAD_ON)
         DEALLOCATE(WA%CR)
      ENDIF

GPACC !$ACC EXIT DATA DELETE(WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
!!        CALL DELETE_TYPED_VAR(WA%WDES1)
!!        CALL ACC_DETACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%WDES1) IF(OFFLOAD_ON)
! two times for pointer
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA,WA%WDES1) IF(OFFLOAD_ON)
      NULLIFY(WA%WDES1)
!!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
      END SUBROUTINE DELWAVA


!************************* SUBROUTINE DELWAVA_PROJ *********************
!
!>  destroy storage for a wavefunctionarray WA
!>  wavefunction character  only
!
!***********************************************************************

      SUBROUTINE DELWAVA_PROJ(WA)

__GPU USE moffload

      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (ASSOCIATED(WA%CPTWFP) .OR. ASSOCIATED(WA%FERWE) .OR. &
          ASSOCIATED(WA%CELEN) .OR. ASSOCIATED(WA%AUX)) THEN
         CALL vtutor%bug("internal error in DELWAVA_PROJ: enities are associated, try DELWAVA", __FILE__, __LINE__)
      ENDIF
GPACC !$ACC EXIT DATA DELETE(WA%CPROJ,WA%CPROJ_RED) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%CPROJ) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%CPROJ_RED) IF(OFFLOAD_ON)
      DEALLOCATE(WA%CPROJ)
      NULLIFY(WA%CPROJ_RED)
GPACC !$ACC EXIT DATA DELETE(WA%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
!!        CALL DELETE_TYPED_VAR(WA%WDES1)
!!        CALL ACC_DETACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%WDES1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%WDES1) IF(OFFLOAD_ON)

      NULLIFY(WA%WDES1)
!!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(WA,1)) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE DELWAVA_PROJ


!************************* FUNCTION PCW ********************************
!
!> index the wavefunction array or character array in WA
!
!***********************************************************************

  FUNCTION PCW( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    COMPLEX(q), POINTER :: PCW(:)
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    PCW=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCW

  FUNCTION PCPROJ( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    GDEF, POINTER :: PCPROJ(:)
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    PCPROJ=>WA%CPROJ(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCPROJ


!************************* SUBROUTINE ORTHON ***************************
!
!> orthogonalize a wavefunction W1 to all other bands
!> including the current band
!>
!> the subroutine uses BLAS 3 calls,
!
!***********************************************************************

  SUBROUTINE ORTHON(NK, W, W1, CQIJ, ISP)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1
    LOGICAL LOVERL
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP
  ! local
    GDEF :: CPRO(W%WDES%NBANDS),CWORK(W%WDES%NPRO)
    REAL(q) :: WFMAG
    INTEGER :: I

#ifdef MPI
    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       CALL vtutor%bug("internal error: ORTHON does not support band-par.", __FILE__, __LINE__)
    ENDIF
#endif
    IF (W1%WDES1%LOVERL) THEN
       CALL OVERL1_SIMPLE(W1%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP), W1%CPROJ(1),CWORK(1))
    ENDIF

    CALL GGEMV( trans , m_ W%WDES%NPLWKP(NK) , W%WDES%NBANDS ,one , W%CW(1,1,NK,ISP), &
                       m_ W%WDES%NRPLWV, W1%CW(1) , 1 , zero ,  CPRO(1), 1)

    IF (W1%WDES1%LOVERL) THEN
       IF (W%WDES%NPRO /= 0) &
            CALL GGEMV( trans ,  W%WDES%NPRO , W%WDES%NBANDS ,one , W%CPROJ(1,1,NK,ISP) , &
            W%WDES%NPROD, CWORK(1), 1 , one ,  CPRO(1), 1)
    ENDIF

    CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPRO(1), W%WDES%NBANDS))

    CALL GGEMM( 'N', 'N' , m_ W%WDES%NPLWKP(NK) , 1 , W%WDES%NBANDS , -one , &
         W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV , CPRO(1) , W%WDES%NBANDS , &
         one , W1%CW(1) , m_ W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0) &
         CALL GGEMM( 'N', 'N' ,  W%WDES%NPRO , 1 , W%WDES%NBANDS  , -one , &
         W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1) , W%WDES%NBANDS , &
         one , W1%CPROJ(1) ,  W%WDES%NPROD  )

  END SUBROUTINE ORTHON

  SUBROUTINE ORTHON1P(NK, W, W1, CQIJ, ISP, NB)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1(:)
    LOGICAL LOVERL
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP,NB(:)
  ! local
    GDEF :: CPRO(W%WDES%NBANDS,SIZE(W1)),CWORK(W%WDES%NPRO,SIZE(W1))
    GDEF,ALLOCATABLE :: CW1(:,:),CPROJ1(:,:)
    REAL(q) :: WFMAG
    INTEGER :: I,BW

    BW=128

#ifdef MPI
    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       CALL vtutor%bug("internal error: ORTHON does not support band-par.", __FILE__, __LINE__)
    ENDIF
#endif

    ALLOCATE(CW1(SIZE(W1(1)%CW,1),SIZE(W1)))
    DO I=1,SIZE(W1)
       CW1(:,I)=W1(I)%CW
    ENDDO

    DO I=1,SIZE(W1)
       IF (W1(I)%WDES1%LOVERL)THEN
          CALL OVERL1_SIMPLE(W1(I)%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP), W1(I)%CPROJ(1),CWORK(1,I))
       ELSE
          CWORK(:,I)=0.0_q
       ENDIF
    ENDDO

    CALL GGEMM( trans , 'N', W%WDES%NBANDS, size(w1), m_ W%WDES%NPLWKP(NK), &
           one , W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV, CW1(1,1) , m_ W%WDES%NRPLWV, zero ,  CPRO(1,1), W%WDES%NBANDS)


    IF (W%WDES%NPRO /= 0) &
         CALL GGEMM( trans , 'N', W%WDES%NBANDS , size(w1), W%WDES%NPRO , one ,&
           W%CPROJ(1,1,NK,ISP) , W%WDES%NPROD, CWORK(1,1), W%WDES%NPRO, one,  CPRO(1,1), W%WDES%NBANDS)

    CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPRO(1,1), W%WDES%NBANDS*SIZE(W1)))

    DO I=1,SIZE(W1)
       CPRO(NB(I),I)=0.d0
       IF(NB(I)-BW>=1)CPRO(1:NB(I)-BW,I)=0.d0
       IF(NB(I)+BW<=SIZE(CPRO,1))CPRO(NB(I)+BW:SIZE(CPRO,1),I)=0.d0
    ENDDO

    CALL GGEMM( 'N', 'N' , m_ W%WDES%NPLWKP(NK) , SIZE(W1) , W%WDES%NBANDS , -one , &
         W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV , CPRO(1,1) , W%WDES%NBANDS , &
         one , CW1(1,1) , m_ W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0)THEN
       DO I=1,SIZE(W1)
          CWORK(:,i)=W1(i)%CPROJ(:)
       ENDDO
       CALL GGEMM( 'N', 'N' ,  W%WDES%NPRO , SIZE(W1) , W%WDES%NBANDS  , -one , &
            W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1,1) , W%WDES%NBANDS , &
            one , CWORK(1,1) ,  W%WDES%NPROD  )
    ENDIF

    DO I=1,SIZE(W1)
       W1(I)%CW=CW1(:,I)
       W1(I)%CPROJ=CWORK(:,I)
    ENDDO

    DEALLOCATE(CW1)

  END SUBROUTINE ORTHON1P


!************************* SUBROUTINE CNORMN ***************************
!
!> this subroutine normalises a wavefunction
!>
!> subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE CNORMN(W, CQIJ, ISP, WSCAL)

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefun1) W

    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP
    REAL(q) WSCAL
  ! local
    GDEF    CP
    REAL(q) WFMAG
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, NPRO2, NPRO2_

    COMPLEX(q), EXTERNAL ::  WZDOTC
!$ACC ROUTINE(ECCP_NL) VECTOR

    PROFILING_START('cnormn')

    WFMAG=WZDOTC(W%WDES1%NPL,W%CW(1),1,W%CW(1),1)
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN
       CP  =0

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(W,CQIJ) REDUCTION(+:CP) &
GPACC !$ACC PRIVATE(NPRO,NPRO_,NPRO2,NPRO2_,NT,LMMAXC) DEFAULT(none) __IF_ASYNC__

! no need to limit the number of threads here
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) &
GPOMP !$OMP PRIVATE(NPRO,NPRO_,NPRO2,NPRO2_,NT,LMMAXC) &
#ifdef IFXWAR3
GPOMP !$OMP MAP(PRESENT,TO:W) REDUCTION(+:CP) IF(OFFLOAD_ON)
#else
GPOMP !$OMP MAP(PRESENT,TO:W,CQIJ) REDUCTION(+:CP) IF(OFFLOAD_ON)
#endif

       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1

__CPU        NPRO =ISPINOR *(W%WDES1%NPRO/2)
__CPU        NPRO_=ISPINOR_*(W%WDES1%NPRO/2)

             DO NI=1,W%WDES1%NIONS
__GPU           NPRO =ISPINOR *(W%WDES1%NPRO/2)
__GPU           NPRO_=ISPINOR_*(W%WDES1%NPRO/2)
                NT=W%WDES1%ITYP(NI)
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE

                NPRO2 =W%WDES1%LMBASE(NI)+NPRO
                NPRO2_=W%WDES1%LMBASE(NI)+NPRO_
                CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO2_+1),W%CPROJ(NPRO2+1),CP)
             ENDDO
          ENDDO
       ENDDO spinor

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       WFMAG=WFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W%WDES1%COMM_INB, WFMAG, 1))

    !-----check that it is non-zero
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller zero write a warning
!=======================================================================
#ifdef MPI
       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN
#endif
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
#ifdef MPI
       ENDIF
#endif
       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF
    CALL WZDSCAL( W%WDES1%NPL ,WSCAL,W%CW(1),1)
    CALL WGDSCAL( W%WDES1%NPRO,WSCAL,W%CPROJ(1),1)

    PROFILING_STOP('cnormn')

  END SUBROUTINE CNORMN


!************************* SUBROUTINE CNORMN_REAL **********************
!
!> performs operations on real space part of wavefunction
!> after a call to CPROJCN and CNORMN
!
!***********************************************************************

  SUBROUTINE CNORMN_REAL(W, W1, ISP, WSCAL, CSCPD )

__GPU USE moffload_struct_def

    IMPLICIT NONE
    TYPE (wavefun1)    W, W1

    INTEGER   ISP
    REAL(q) WSCAL
    COMPLEX(q) :: CSCPD
  ! local
    INTEGER ISPINOR, K, KK

    PROFILING_START('cnormn_real')

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W,W%WDES1,W%WDES1%GRID,W%WDES1%GRID%RL,W1) PRIVATE(KK) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(KK) &
GPOMP !$OMP MAP(PRESENT,TO:W,W%WDES1,W%WDES1%GRID,W%WDES1%GRID%RL,W1) IF(OFFLOAD_ON)
    DO ISPINOR=0,W%WDES1%NRSPINORS-1
       DO K=1,W%WDES1%GRID%RL%NP
          KK=K+ISPINOR*W%WDES1%GRID%MPLWV
          W%CR(KK)=(W%CR(KK)-CSCPD*W1%CR(KK))*WSCAL
       ENDDO
    ENDDO

    PROFILING_STOP('cnormn_real')

  END SUBROUTINE CNORMN_REAL


!************************* SUBROUTINE CNORMA ***************************
!
!> this subroutine calculates the norm of a wavefunction
!
!> @details @ref openmp :
!> a ZDOTC call is replaced by an explicit OMP PARALLEL loop
!> (for some reason BLAS1 call thread really badly).
!> The nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\", that is distributed over all available
!> threads.
!
!***********************************************************************

  SUBROUTINE CNORMA(W, CQIJ, ISP, WSCAL)

__GPU USE moffload

    IMPLICIT NONE

    TYPE (wavefun1) W
    OVERLAP :: CQIJ(:,:,:,:)
    INTEGER :: ISP
    REAL(q) :: WSCAL
  ! local
    GDEF    :: CP, CTMP
    REAL(q) :: WFMAG
    INTEGER :: ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NI, LMMAXC, LMDIM

    REAL(q), EXTERNAL ::  WDDOT

!$  INTEGER I
!$ACC ROUTINE(ECCP_NL) VECTOR

    PROFILING_START('cnorma')

#if defined(OFFLOADING) || ! defined(_OPENMP)
!      WFMAG=ZDOTC(W%WDES1%NPL,W%CW(1),1,W%CW(1),1)
       WFMAG=WDDOT(W%WDES1%NPL*2,W%CW(1),1,W%CW(1),1)
#else
       WFMAG=0
!$OMP PARALLEL DO SIMD SHARED(W) REDUCTION(+:WFMAG)
       DO I=1,W%WDES1%NPL
          WFMAG=WFMAG+CONJG(W%CW(I))*W%CW(I)
       ENDDO
!$OMP END PARALLEL DO SIMD
#endif
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN

       CP=0
       LMDIM=SIZE(CQIJ,1)

       PROFILING_START('eccp_nl')

GPACC !$ACC ENTER DATA COPYIN(W) __IF_ASYNC__
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(CQIJ,W) REDUCTION(+:CP) &
GPACC !$ACC PRIVATE(NT,LMMAXC,NPRO,NPRO_,CTMP) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(TO:W) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) &
GPOMP !$OMP PRIVATE(NT,LMMAXC,NPRO,NPRO_,CTMP) &
#ifdef IFXWAR3
GPOMP !$OMP MAP(PRESENT,TO:W) REDUCTION(+:CP) IF(OFFLOAD_ON)
#else
GPOMP !$OMP MAP(PRESENT,TO:W,CQIJ) REDUCTION(+:CP) IF(OFFLOAD_ON)
#endif

__CPU !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,CTMP) &
__CPU !$OMP SHARED(W,CQIJ,LMDIM,ISP) REDUCTION(+:CP)

       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1
             DO NI=1,W%WDES1%NIONS
                NT=W%WDES1%ITYP(NI)
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE
                NPRO =ISPINOR *(W%WDES1%NPRO/2)+W%WDES1%LMBASE(NI)
                NPRO_=ISPINOR_*(W%WDES1%NPRO/2)+W%WDES1%LMBASE(NI)

                CTMP=0; CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO_+1),W%CPROJ(NPRO+1),CTMP)
                CP=CP+CTMP
             ENDDO
          ENDDO
       ENDDO spinor

GPACC !$ACC EXIT DATA DELETE(W) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W) IF(OFFLOAD_ON)

__CPU !$OMP END PARALLEL DO

       PROFILING_STOP('eccp_nl')

       WFMAG=WFMAG+CP
    ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    CALLMPI( M_sum_d(W%WDES1%COMM_INB, WFMAG, 1))

    !-----check that it is non-zero
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller zero write a warning
!=======================================================================
#ifdef MPI
       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN
#endif
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
#ifdef MPI
       ENDIF
#endif
       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF

    PROFILING_STOP('cnorma')

  END SUBROUTINE CNORMA


!************************* SUBROUTINE CINPROD **************************
!
!> this subroutine calculates the inproduct between two wavefunctions
!> ~~~
!>  <W2 | S | W1>      = W2(G)*  W1(G) +  \sum_i W2_i*Q_ij W1_j
!> ~~~
!***********************************************************************

  SUBROUTINE CINPROD(W1,W2,CQIJ,CWFMAG)
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (wavefun1)    W2

    OVERLAP   CQIJ(:,:,:,:)
    COMPLEX(q) CWFMAG
  ! local
    GDEF      CP
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT

    CWFMAG=ZDOTC(W1%WDES1%NPL,W2%CW(1),1,W1%CW(1),1)
    CP  =0
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       NPRO=0

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W1%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)

             NIS =1
             DO NT=1,W1%WDES1%NTYP
                LMMAXC=W1%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CP)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W1%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       CWFMAG=CWFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W1%WDES1%COMM_INB, CWFMAG, 2))
  END SUBROUTINE CINPROD


!************************* SUBROUTINE PROJCN ***************************
!
!> this subroutine projects out from one wavefunction
!> CF another wavefunction  CPRO
!>
!> subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE PROJCN(W1, W2, CQIJ, ISP, CSCPD)

__GPU USE moffload

    IMPLICIT NONE

    TYPE (wavefun1)    W1,W2
    GDEF    CADD
    OVERLAP CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CSCPD
    ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL ::  WZDOTC
GPACC !$ACC ROUTINE(ECCP_NL) VECTOR

    CSCPD= GREAL(WZDOTC(W1%WDES1%NPL,W2%CW(1),1,W1%CW(1),1))
!=======================================================================
! if necessary caclulate <p| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       CADD=0
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(W1,W1%WDES1,W2,CQIJ) REDUCTION(+:CADD) &
GPACC !$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) &
GPOMP !$OMP MAP(PRESENT,TO:W1,W1%WDES1,W2,CQIJ) REDUCTION(+:CADD) &
GPOMP !$OMP PRIVATE(NPRO,NPRO_,NT,LMMAXC) IF(OFFLOAD_ON)
       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
             DO NI=1,W1%WDES1%NIONS
                   NT=W1%WDES1%ITYP(NI)
                   LMMAXC=W1%WDES1%LMMAX(NT)
                   IF (LMMAXC==0) CYCLE

                   NPRO =ISPINOR *(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
                   NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CADD)
             ENDDO
          ENDDO
       ENDDO spinor

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CSCPD=GREAL(CSCPD+CADD)
    ENDIF
!=======================================================================
! performe orthogonalisations
!=======================================================================
    CALLMPI( M_sum_z(W1%WDES1%COMM_INB, CSCPD, 1))

    CALL WZAXPY(W1%WDES1%NPL ,-CSCPD,W2%CW(1)   ,1,W1%CW(1)   ,1)
    IF (W1%WDES1%LGAMMA) THEN
       CALL WDAXPY(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ELSE
       CALL WZAXPY(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ENDIF

    RETURN
  END SUBROUTINE PROJCN


!************************ SUBROUTINE W1_GATHER ************************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!
!**********************************************************************

  SUBROUTINE W1_GATHER( W, NB1, NB2, ISP, W1)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=(NB2-NB1+1)*W%WDES%NB_PAR

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CW(1), &
               SIZE(W1(NN)%CW),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
          ENDIF
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif
    ENDIF
#endif
  END SUBROUTINE W1_GATHER


!************************ SUBROUTINE W1_GATHER_N **********************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the routine only up to NLOC bands are collected
!
!**********************************************************************

  SUBROUTINE W1_GATHER_N( W, NB1, NB2, ISP, W1, NLOC)
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> total number of bands to be collected
    INTEGER:: NLOC

  ! local
    INTEGER :: NN, N, NCPU

    PROFILING_START('w1_gather_n')

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       IF (NN>NLOC) EXIT
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CW(1), &
               SIZE(W1(NN)%CW),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
          ENDIF
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_n')

  END SUBROUTINE W1_GATHER_N


!************************ SUBROUTINE W1_GATHER_GLB ********************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied
!
!**********************************************************************


  SUBROUTINE W1_GATHER_GLB( W, NB1, NB2, ISP, W1)

    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror

    PROFILING_START('w1_gather_glb')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
          CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
#ifdef USENCCL
       IF ( W%WDES%COMM_INTER%LUSENCCL ) CALL nccl_group_start()
#endif
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
              CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                   W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
              CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                   W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
          ENDIF
       ENDDO
#ifdef USENCCL
       IF ( W%WDES%COMM_INTER%LUSENCCL ) CALL nccl_group_end()
#endif

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb')

  END SUBROUTINE W1_GATHER_GLB

!************************ SUBROUTINE W1_GATHER_GLB ********************
!
! This subroutine gathers a set of wavefunctions starting
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band
! indices are supplied
!
!**********************************************************************


  SUBROUTINE W1_GATHER_GLB_( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

  ! local
    INTEGER :: N_INTO_TOT, N, NB_LOCAL

    PROFILING_START('w1_gather_glb_')

    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(N_INTO_TOT)%WDES1,NB_LOCAL,ISP),W1(N_INTO_TOT))
          CALL FFTWAV_W1(W1(N_INTO_TOT))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CW(1), &
               SIZE(W1(N_INTO_TOT)%CW),MOD(N-1,W%WDES%NB_PAR)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N-1,W%WDES%NB_PAR)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
#endif
          ENDIF
       ENDDO
#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier(W%WDES%COMM_INTER))
#endif
    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_')

  END SUBROUTINE W1_GATHER_GLB_


!************************ SUBROUTINE W1_GATHER_GLB_NOFFT ********************
!> identical to W1_GATHER_GLB but does not FFT to real space
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOFFT( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror

    PROFILING_START('w1_gather_glb_nofft')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_nofft')

  END SUBROUTINE W1_GATHER_GLB_NOFFT


!************************ SUBROUTINE W1_IGATHER_GLB *******************
!
!> Gathers a set of wavefunctions starting from band NB1 until NB2
!> (global band indices) to all nodes using non-blocking bcast_from.
!>
!> In case the code is compiled with-Dshmem_bcast_buffer communication
!> will be between ranks within COMM_inter_node, i.e., only those ranks
!> within COMM_INTER that are NOT on the same physical node will talk to
!> eachother. MPI communication between ranks within COMM_INTER that
!> reside on the same node is not necessary since they access a common
!> shared memory segment.
!>
!> In case the code is NOT compiled with-Dshmem_bcast_buffer, the
!> communication will be between all ranks in COMM_INTER.
!
!**********************************************************************
#ifdef shmem_bcast_buffer
  SUBROUTINE W1_IGATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

  ! local
    INTEGER :: NI, N, NB_LOCAL

    INTEGER :: IDO(NB2-NB1+1)

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))

    PROFILING_START('w1_igather_glb')

    CALL M_barrier(W%WDES%COMM_intra_node)

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(NI)%WDES1,NB_LOCAL,ISP),W1(NI))
          CALL FFTWAV_W1(W1(NI))
          IDO(NI)=W%WDES%COMM_inter_node%NODE_ME
       ELSE
          IDO(NI)=0
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_inter_node%NCPU>1) THEN

       PROFILING_START('stage_1')

       CALLMPI( M_sum_i(W%WDES%COMM_inter_node,IDO,NB2-NB1+1))

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1
          IF (IDO(NI)>0) THEN
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CW(1), &
                  SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CR(1), &
                  SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
             IF (W%WDES%LOVERL) THEN
                nrequests=nrequests+1
#ifndef gammareal
                CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#else
                CALL M_ibcast_d_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#endif
             ENDIF
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

       PROFILING_STOP('stage_1')

    ENDIF

    CALL M_barrier(W%WDES%COMM_intra_node)
!!#ifndef shmem_bcast_buffer
!!!   CALL M_barrier(W%WDES%COMM_INTER)
!!
!!    IF (W%WDES%COMM_intra_node%NCPU>1) THEN
!!       DO N=1,NB2-NB1+1
!!          IF (IDO(N)>0) IDO(N)=W%WDES%COMM_intra_node%NODE_ME
!!       ENDDO
!!
!!       CALLMPI( M_sum_i(W%WDES%COMM_intra_node,IDO,NB2-NB1+1))
!!
!!       PROFILING_START('stage_2')
!!
!!       nrequests=0
!!       DO N=NB1,NB2
!!          NI=N-NB1+1
!!          IF (IDO(NI)>0) THEN
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CW(1), &
!!                  SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CR(1), &
!!                  SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
!!             IF (W%WDES%LOVERL) THEN
!!                nrequests=nrequests+1
!!#ifndef gammareal
!!                CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#else
!!                CALL M_ibcast_d_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#endif
!!             ENDIF
!!          ENDIF
!!       ENDDO
!!
!!       CALL M_waitall(nrequests,requests(1))
!!
!!       PROFILING_STOP('stage_2')
!!
!!    ENDIF
!!#endif
#endif
    PROFILING_STOP('w1_igather_glb')

  END SUBROUTINE W1_IGATHER_GLB
#else
  SUBROUTINE W1_IGATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

  ! local
    INTEGER :: NI, N, NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))

    PROFILING_START('w1_igather_glb')

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(NI)%WDES1,NB_LOCAL,ISP),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CW(1), &
               SIZE(W1(NI)%CW),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
               SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))

          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ibcast_d_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_igather_glb')

  END SUBROUTINE W1_IGATHER_GLB
#endif

!************************ SUBROUTINE W1_REDUCE_GLB ********************
!
!> Take the sum of W1(i)%CW, W1(i)%CPROJ, for i = 1, NB2-NB1+1 over all
!> ranks in COMM_inter_node, using reduce_to.
!>
!> The result is reduced onto the MPI-rank in COMM_inter_node that would
!> normally own (part of) the band with the global index N = NB1, ...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_REDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1):: WDES1
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2

  ! local
    INTEGER :: N_INTO_TOT, N
    INTEGER :: ierror

    INTEGER :: IDO(NB2-NB1+1)

    PROFILING_START('w1_reduce_glb')
#ifdef MPI
    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,WDES1%NB_PAR)+1==WDES1%NB_LOW) THEN
          IDO(N_INTO_TOT)=WDES1%COMM_inter_node%NODE_ME
       ELSE
          IDO(N_INTO_TOT)=0
       ENDIF
    ENDDO

    CALLMPI( M_sum_i(WDES1%COMM_inter_node,IDO,NB2-NB1+1))

    ! reduce W1
    IF (WDES1%DO_REDIS) THEN

       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          IF (IDO(N_INTO_TOT)>0) THEN
             CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CW(1), &
                  SIZE(W1(N_INTO_TOT)%CW),IDO(N_INTO_TOT))

             IF (WDES1%LOVERL) THEN
#ifndef gammareal
                CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
                     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
#else
                CALL M_reduce_d_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
                     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
#endif
             ENDIF
          ENDIF
       ENDDO

    ENDIF
#endif
    PROFILING_STOP('w1_reduce_glb')

  END SUBROUTINE W1_REDUCE_GLB


!************************ SUBROUTINE W1_IREDUCE_GLB *******************
!
!> Take the sum of W1(i)%CW, W1(i)%CPROJ, for i = 1, NB2-NB1+1 over all
!> ranks in COMM_INTER, using non-blocking reduce_to.
!>
!> The result is reduced onto the MPI-rank in COMM_INTER that would
!> normally own (part of) the band with the global index N = NB1, ...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_IREDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2

  ! local
    INTEGER :: NI, N
    INTEGER :: requests(2*(NB2-NB1+1))
    INTEGER :: nrequests

    PROFILING_START('w1_ireduce_glb')

    ! reduce W1
#ifdef MPI
    IF (WDES1%DO_REDIS) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CW(1), &
               SIZE(W1(NI)%CW),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))

          IF (WDES1%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ireduce_d_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_ireduce_glb')

  END SUBROUTINE W1_IREDUCE_GLB


!************************ SUBROUTINE W1_GATHER_GLB_NOCR ***************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied and the real space part is not
!> communicated
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU

    PROFILING_START('w1_gather_glb_nocr')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
          ENDIF
       ENDDO
    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_nocr')

  END SUBROUTINE W1_GATHER_GLB_NOCR

#ifdef use_shmem
!************************ SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM *********
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM(W,NB1,NB2,ISP,W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    ! local variables
    COMPLEX(q) :: CWORK(W1(1)%WDES1%NRPLWV)
    INTEGER :: N_INTO_TOT,N,NCPU

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
!         CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          ! the node that owns this band stores it onto CWORK
          IF (MOD(N_INTO_TOT-1,NCPU)+1==W%WDES%COMM_INTER%NODE_ME) CWORK=W1(N)%CW
          ! and broadcasts CWORK to all other nodes
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CWORK(1), &
               SIZE(CWORK),MOD(N_INTO_TOT-1,NCPU)+1)
          ! but only the root nodes in COMM_SHMEM need to store CWORK back into W1
          IF (MOD(N_INTO_TOT-1,NCPU)+1/=W%WDES%COMM_INTER%NODE_ME &
             .AND.W%WDES%COMM_SHMEM%NODE_ME==1) W1(N)%CW=CWORK

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
          ENDIF
       ENDDO
    ENDIF
#endif

  END SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM


!************************ SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM *********
!
!> Gathers a set of wavefunctions including band NB1 to NB2 over
!> all k-points to all nodes using shared memory
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: WF(:,:)
    ! local variables
    TYPE(wavefun1):: WAUX
    TYPE(wavedes1), TARGET :: WDESAUX
    COMPLEX(q) :: CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
    GDEF CPROJ(W%WDES%NPROD)
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS

    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
             ! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
    ! distribute WF to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             ! the node that owns this band stores it onto CWORK
             IF (MOD(N_INTO_TOT-1,NCPU)+1==W%WDES%COMM_INTER%NODE_ME) THEN
                CWORK(1:SIZE(WF(N,IK)%CR))=WF(N,IK)%CR(:); CPROJ(1:SIZE(WF(N,IK)%CPROJ))=WF(N,IK)%CPROJ(:)
             ENDIF
             ! and broadcasts CWORK to all other nodes
             CALL M_bcast_z_from(W%WDES%COMM_INTER,CWORK(1),SIZE(CWORK),MOD(N_INTO_TOT-1,NCPU)+1)
             IF (W%WDES%LOVERL) THEN
#ifdef gammareal
                CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1),SIZE(CPROJ),MOD(N_INTO_TOT-1,NCPU)+1)
#else
                CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1),SIZE(CPROJ),MOD(N_INTO_TOT-1,NCPU)+1)
#endif
             ENDIF
             ! but only the root nodes in COMM_SHMEM need to store CWORK back into WF
             IF (MOD(N_INTO_TOT-1,NCPU)+1/=W%WDES%COMM_INTER%NODE_ME.AND.W%WDES%COMM_SHMEM%NODE_ME==1) THEN
                WF(N,IK)%CR=CWORK(1:SIZE(WF(N,IK)%CR)); WF(N,IK)%CPROJ(:)=CPROJ(1:SIZE(WF(N,IK)%CPROJ))
             ENDIF
          ENDDO
       ENDDO
#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif
    ENDIF
#endif
  END SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM
#endif

!************************ SUBROUTINE W1_GATHER_GLB_ALLK ***************
!
!> Gathers a set of wavefunctions including band NB1 to NB2 over
!> all k-points to all nodes
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: WF(:,:)
    ! local
    TYPE (wavefun1):: WAUX
    TYPE(wavedes1), TARGET :: WDESAUX
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS


    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

    ! W1 contains bands NB1 to NB2 for all k-points
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
             ! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CR(1), &
                  W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)

             IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                     W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
                CALL M_bcast_d_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                     W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
             ENDIF
          ENDDO
       ENDDO

#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
  END SUBROUTINE W1_GATHER_GLB_ALLK


!************************ SUBROUTINE W1_GATHER_ARRAY ******************
!
!> This subroutine gathers a set of orbitals
!> from band NB1 until NB2 to all nodes
!>
!> compared to previous versions the set in collected
!> into work arrays
!
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY( W, NB1, NB2, ISP, W1, CR, CPROJ)
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array for FFT
    TYPE (wavefun1):: W1
    !> collected real space orbitals
    GDEF    :: CR(:,:)
    !> collected projected orbitals
    GDEF    :: CPROJ(:,:)

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )
       CALL FFTWAV_W1( W1)

       IF (NN> SIZE(CR,2) .OR. NN >SIZE(CPROJ,2)) THEN
          CALL vtutor%bug("internal error in W1_GATHER_ARRAY: bound exceed " // str(NN) // " " // &
             str(SIZE(CR,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
       ENDIF
       IF (SIZE(CR,1)/=W1%WDES1%GRID%NPLWV) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY: size mismatch ',SIZE(CR,1),W1%WDES1%GRID%NPLWV
       ENDIF

       CR(1:W1%WDES1%GRID%RL%NP, NN)=W1%CR(1:W1%WDES1%GRID%RL%NP)
       ! pad CR with zeros (just in case)
       CR(W1%WDES1%GRID%RL%NP+1:SIZE(CR,1), NN)=0
       CPROJ(:,NN)=W1%CPROJ(:)

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
        DO NN=1,NCPU
#ifndef gammareal
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
#else
          CALL M_bcast_d_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
#endif
          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
            ENDIF
         ENDDO
       ENDIF
#endif
      ENDDO

  END SUBROUTINE W1_GATHER_ARRAY


!************************ SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL *******
!
!> This subroutine gathers a set of orbitals
!> from band NB1 until NB2 to all nodes
!>
!> this version collects the place wave coefficients and the
!> CPROJ coefficients
!
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL( W, NB1, NB2, ISP, W1, CG, CPROJ)

GPACC USE moffload_struct_def

    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array for FFT
    TYPE (wavefun1):: W1
    !> collected real space orbitals
    COMPLEX(q):: CG(:,:)
    !> collected projected orbitals
    GDEF       :: CPROJ(:,:)

    COMPLEX(q) :: CWBUFF(SIZE(CG,1))       ! send buffer
    GDEF       :: CPROJBUFF(SIZE(CPROJ,1)) ! send buffer for projectors

  ! local
    INTEGER :: NN, N, NCPU

    PROFILING_START('w1_gather_array_reciprocal')

!$ACC ENTER DATA CREATE(CWBUFF,CPROJBUFF) __IF_ASYNC__

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), __FILE__, __LINE__)
    ENDIF

    NN=(NB2-NB1)*W%WDES%NB_PAR+W%WDES%NB_LOW
    IF (NN>SIZE(CG,2) .OR. NN>SIZE(CPROJ,2)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: bound exceeded " // str(NN) &
          // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    ! NB1 and NB2 are local indices for collecting the bands
    DO N=NB1,NB2
       !NB_LOW is the off-set of the local node
       !NN is index into CG array
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       ! copy data over to return arrays: CG and CPROJ
!$ACC KERNELS PRESENT(CG,CPROJ,W1) __IF_ASYNC__
       CG(1:W1%WDES1%NPL ,NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1:,NN)=0.0_q                 ! pad with zero

       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
!$ACC END KERNELS
#ifdef MPI
      ! distribute W1 to all nodes, CG(1,(N-NB)*NCPU)
       IF (W%WDES%DO_REDIS) THEN
          ! MPI_IN_PLACE was tested but found to be much slower (code has been removed)
!$ACC KERNELS PRESENT(CWBUFF,CPROJBUFF,CG,CPROJ) __IF_ASYNC__
          CPROJBUFF(:)=CPROJ(:,NN)
          CWBUFF(:)   =CG(:,NN)
!$ACC END KERNELS
          CALL M_allgathero_z(W%WDES%COMM_INTER,SIZE(CWBUFF),CWBUFF(1),CG(1,(N-NB1)*NCPU+1))

          IF (W%WDES%LOVERL) THEN
#ifdef gammareal
             CALL M_allgathero_d(W%WDES%COMM_INTER,SIZE(CPROJBUFF),CPROJBUFF(1),CPROJ(1,(N-NB1)*NCPU+1))
#else
             CALL M_allgathero_z(W%WDES%COMM_INTER,SIZE(CPROJBUFF),CPROJBUFF(1),CPROJ(1,(N-NB1)*NCPU+1))
#endif
          ENDIF
       ENDIF
#endif
    ENDDO

!$ACC EXIT DATA DELETE(CWBUFF,CPROJBUFF) __IF_ASYNC__

    PROFILING_STOP('w1_gather_array_reciprocal')

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL

!
!> A version of W1_GATHER_ARRAY_RECIPROCAL using non-blocking allgather operations
!
  SUBROUTINE W1_IGATHER_ARRAY_RECIPROCAL( W, NB1, NB2, ISP, W1, CG, CPROJ)
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array for FFT
    TYPE (wavefun1):: W1
    !> collected real space orbitals
    COMPLEX(q):: CG(:,:)
    !> collected projected orbitals
    GDEF       :: CPROJ(:,:)

    COMPLEX(q) :: CWBUFF(SIZE(CG,1))       ! send buffer
    GDEF       :: CPROJBUFF(SIZE(CPROJ,1)) ! send buffer for projectors

  ! local
    INTEGER :: NN, N, NCPU

    INTEGER :: nrequests, requests(2*(NB2-NB1+1))

    PROFILING_START('w1_igather_array_reciprocal')

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_ARRAY_RECIPROCAL: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("W1_IGATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), __FILE__, __LINE__)
    ENDIF

    NN=(NB2-NB1)*W%WDES%NB_PAR+W%WDES%NB_LOW
    IF (NN>SIZE(CG,2) .OR. NN>SIZE(CPROJ,2)) THEN
       CALL vtutor%bug("W1_IGATHER_ARRAY_RECIPROCAL: bound exceeded " // str(NN) &
          // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    nrequests=0

    ! NB1 and NB2 are local indices for collecting the bands
    DO N=NB1,NB2
       !NB_LOW is the off-set of the local node
       !NN is index into CG array
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       ! copy data over to return arrays: CG and CPROJ
       CG(1:W1%WDES1%NPL ,NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1:,NN)=0.0_q                 ! pad with zero

       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
#ifdef MPI
      ! distribute W1 to all nodes, CG(1,(N-NB)*NCPU)
       IF (W%WDES%DO_REDIS) THEN
          ! MPI_IN_PLACE was tested but found to be much slower (code has been removed)
          CPROJBUFF(:)=CPROJ(:,NN)
          CWBUFF(:)   =CG(:,NN)

          nrequests=nrequests+1
          CALL M_iallgathero_z(W%WDES%COMM_INTER,SIZE(CWBUFF),CWBUFF(2),CG(1,(N-NB1)*NCPU+1),requests(nrequests))

          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifdef gammareal
             CALL M_iallgathero_d(W%WDES%COMM_INTER,SIZE(CPROJBUFF),CPROJBUFF(1),CPROJ(1,(N-NB1)*NCPU+1),requests(nrequests))
#else
             CALL M_iallgathero_z(W%WDES%COMM_INTER,SIZE(CPROJBUFF),CPROJBUFF(1),CPROJ(1,(N-NB1)*NCPU+1),requests(nrequests))
#endif
          ENDIF
       ENDIF
#endif
    ENDDO

    CALLMPI( M_waitall(nrequests,requests(1)))

    PROFILING_STOP('w1_igather_array_reciprocal')

  END SUBROUTINE W1_IGATHER_ARRAY_RECIPROCAL

!
!> old version, uses bcast and is much slower
!
  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD( W, NB1, NB2, ISP, W1, CG, CPROJ)
    TYPE (wavespin) W        !< wavefunction
    INTEGER :: NB1           !< starting band
    INTEGER :: NB2           !< final band
    INTEGER :: ISP           !< spin
    TYPE (wavefun1):: W1     !< array for FFT
    COMPLEX(q) :: CG(:,:)    !< collected real space orbitals
    GDEF       :: CPROJ(:,:) !< collected projected orbitals

  ! local
    INTEGER :: NN, N, NCPU

    PROFILING_START('w1_gather_array_reciprocal_old')

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY_RECIPROCAL: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), __FILE__, __LINE__)
    ENDIF

    NN=(NB2-NB1)*W%WDES%NB_PAR+W%WDES%NB_LOW
    IF (NN>SIZE(CG,2) .OR. NN>SIZE(CPROJ,2)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: bound exceeded " // str(NN) &
          // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       !NB_LOW will be the off-set of the local node
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       !copy the local data to the position with global index in the strip
       CG(1:W1%WDES1%NPL, NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1: SIZE(CG,1), NN)=0

       !the first index CPROJ size can be smaller than the second one (less cpus/tau than total number)
       !but this copies only the data for valid indices in the first array (or?)
       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))

       ! distribute CG and CPROJ to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
          DO NN=1,NCPU
             !so broadcast the data from node holding index NN to all other nodes
             !I wonder how efficient this is, might be all right 10000 coefficients is 160kB of data,
             !a bit on the low side but fine
             CALL M_bcast_z_from(W%WDES%COMM_INTER,CG(1,NN+(N-NB1)*NCPU), &
                  SIZE(CG,1),MOD(NN-1,NCPU)+1)

             IF (W%WDES%LOVERL) THEN
                ! OK, the point is that the W%WDES%NPROD array can be longer than the size of CPROJ
                ! since both are a multiple of NCPU but once the total number and in the second case per tau
#ifndef gammareal
                CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
                     SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)
#else
                CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
                     SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)
#endif
             ENDIF
          ENDDO
       ENDIF
#endif
    ENDDO

    PROFILING_STOP('w1_gather_array_reciprocal_old')

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD


!************************ SUBROUTINE W1_GATHER_DISTR ******************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied
!
!**********************************************************************

  SUBROUTINE W1_GATHER_DISTR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    ! local variables
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: NBMIN, NBMAX
    TYPE (wavefun1):: WTMP

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .FALSE.)

    ! establish global band index interval
    NBMIN=NB1
    NBMAX=NB2
#ifdef MPI
    DO N=1,NCPU
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMIN, 1, n)
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMAX, 1, n)
       IF (NB1<NBMIN) NBMIN=NB1
       IF (NB2>NBMAX) NBMAX=NB2
    ENDDO
#endif

    DO N_INTO_TOT=NBMIN,NBMAX

       ! if band resides on this node copy it to WTMP
       IF (MOD(N_INTO_TOT-W%WDES%NB_LOW,NCPU)==0) THEN
          N=(N_INTO_TOT-W%WDES%NB_LOW)/NCPU+1
          CALL W1_COPY( ELEMENT( W, WTMP%WDES1, N, ISP), WTMP )
       ENDIF

#ifdef MPI
       ! broadcast WTMP from the node where it resides
       CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CW(1), &
            SIZE(WTMP%CW),MOD(N_INTO_TOT-1,NCPU)+1)

       IF (W%WDES%LOVERL) THEN
#ifndef gammareal
          CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          CALL M_bcast_d_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDIF
#endif
       ! if this band is targeted to reside on this node
       ! copy WTMP to W1(N_INTO_TOT-NB1+1)
       IF (N_INTO_TOT>=NB1 .AND. N_INTO_TOT<=NB2) THEN
          N=N_INTO_TOT-NB1+1
          CALL W1_COPY(WTMP,W1(N))
       ENDIF
    ENDDO

    ! FFT to real space
    DO N=1,(NB2-NB1)+1
       CALL FFTWAV_W1(W1(N))
    ENDDO

    CALL DELWAV(WTMP, .FALSE.)

  END SUBROUTINE W1_GATHER_DISTR

!************************ SUBROUTINE W1_GATHER_KSEL_NOFFT *************
!
! This subroutine gathers a set of wavefunctions starting
! from band NB1 until NB2 and distributes the data over k in a round
! robin fashion
! the data distribution is based on an index (k-point index) supplied as
! the last argument
!
!**********************************************************************


  SUBROUTINE W1_GATHER_KSEL_NOFFT( W, NB1, NB2, ISP, W1, NODE, NK)

    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NODE          ! node to recieve k-point
    INTEGER :: NK  ! k-points intex

   ! local
    INTEGER :: NN, N, NLOC, NCPU
    INTEGER :: NCP

    PROFILING_START('w1_gather_ksel_nofft')

    NCPU=W%WDES%NB_PAR
    MPLWV = SIZE(W%CPTWFP(:,1,NK,ISP))
    NPROD = W%WDES%NPROD

#ifdef USENCCL
    CALL nccl_group_start()
#endif
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          NCP = N
       ENDIF
#ifdef MPI
       IF (MOD(NODE-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN /= W%WDES%NB_LOW) THEN
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , &
                        W1(N_INTO_TOT-NB1+1)%CW(1), SIZE(W1(N_INTO_TOT-NB1+1)%CW))

                   IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                      CALL M_recv_z(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                      CALL M_recv_d(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
             CALL M_send_z(W%WDES%COMM_INTER, MOD(NODE-1,NCPU)+1, &
                      W%CPTWFP(:,NCP,NK,ISP),SIZE(W%CPTWFP(:,NCP,NK,ISP)))

             IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                CALL M_send_z(W%WDES%COMM_INTER, MOD(NODE-1,NCPU)+1, &
                     W%CPROJ(:,NCP,NK,ISP), W%WDES%NPROD)
#else
                CALL M_send_d(W%WDES%COMM_INTER, MOD(NODE-1,NCPU)+1, &
                     W%CPROJ(:,NCP,NK,ISP), W%WDES%NPROD)
#endif
             ENDIF
          ENDIF
       ENDIF
#endif
    ENDDO
#ifdef USENCCL
    CALL nccl_group_end()
#endif

    IF (MOD(NODE-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE(N_INTO_TOT,NCP) PRESENT(W1,W) __IF_ASYNC__
       DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
__CPU     N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
__CPU     IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
__CPU        NCP = N
__CPU        DO NN = 1, MPLWV
__CPU           W1(N_INTO_TOT-NB1+1)%CW(NN) = W%CPTWFP(NN,NCP,NK,ISP)
__CPU        ENDDO
__CPU     ENDIF
__GPU     DO NN = 1, MPLWV
__GPU        N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
__GPU        IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
__GPU           NCP = N
__GPU           W1(N_INTO_TOT-NB1+1)%CW(NN) = W%CPTWFP(NN,NCP,NK,ISP)
__GPU        ENDIF
__GPU     ENDDO
       ENDDO
!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE(N_INTO_TOT,NCP) PRESENT(W1,W) __IF_ASYNC__
       DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
__CPU     N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
__CPU     IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
__CPU        NCP = N
__CPU        DO NN = 1, NPROD
__CPU           W1(N_INTO_TOT-NB1+1)%CPROJ(NN) = W%CPROJ(NN,NCP,NK,ISP)
__CPU        ENDDO
__CPU     ENDIF
__GPU     DO NN = 1, NPROD
__GPU        N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
__GPU        IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
__GPU           NCP = N
__GPU           W1(N_INTO_TOT-NB1+1)%CPROJ(NN) = W%CPROJ(NN,NCP,NK,ISP)
__GPU        ENDIF
__GPU     ENDDO
       ENDDO
    ENDIF

    PROFILING_STOP('w1_gather_ksel_nofft')

  END SUBROUTINE W1_GATHER_KSEL_NOFFT

!************************ SUBROUTINE W1_GATHER_KSEL *******************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 and distributes the data over k in a round
!> robin fashion
!>
!> the data distribution is based on an index (k-point index) supplied as
!> the last argument
!
!**********************************************************************


  SUBROUTINE W1_GATHER_KSEL( W, NB1, NB2, ISP, W1, NK)

    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> k-points intex
    INTEGER:: NK

   ! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    PROFILING_START('w1_gather_ksel')

    NCPU=W%WDES%NB_PAR

    !$ACC ENTER DATA CREATE(WTMP) __IF_ASYNC__
    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF
#ifdef MPI
#ifdef USENCCL
       CALL nccl_group_start()
#endif
       IF (MOD(NK-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
          PROFILING_START('w1_gather_ksel_recv')
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
!                   WRITE(*,*) W%WDES%NB_LOW, 'receive', N_INTO_TOT-NB1+1, 'from' , NN
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , &
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

                   IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                      CALL M_recv_z(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                      CALL M_recv_d(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
          PROFILING_STOP('w1_gather_ksel_recv')
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
            PROFILING_START('w1_gather_ksel_send')
!             WRITE(*,*) W%WDES%NB_LOW,'send', N_INTO_TOT, 'to' , MOD(NK-1,NCPU)+1
             CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, &
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

             IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, &
                     WTMP%CPROJ(1), W%WDES%NPROD)
#else
                CALL M_send_d(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, &
                     WTMP%CPROJ(1), W%WDES%NPROD)
#endif
             ENDIF
            PROFILING_STOP('w1_gather_ksel_send')
          ENDIF
       ENDIF
#ifdef USENCCL
       CALL nccl_group_end()
#endif
#endif
    ENDDO

    CALL DELWAV(WTMP, .TRUE.)
    !$ACC EXIT DATA DELETE(WTMP) __IF_ASYNC__

    PROFILING_STOP('w1_gather_ksel')

  END SUBROUTINE W1_GATHER_KSEL



!************************ SUBROUTINE W1_GATHER_KNODESEL ***************
!
!> This subroutine is similar to #W1_GATHER_KSEL,
!> but at variance with
!> it, one explicitly select the node on which you gather the wave
!> fucntion (instead of distributing k-points in round robin fashion).
!
!**********************************************************************

  SUBROUTINE W1_GATHER_KNODESEL( W, NB1, NB2, ISP, W1, NODE)
    !> wavefunction
    TYPE (wavespin) W
    !> starting band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2
    !> spin
    INTEGER:: ISP
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> which node receives
    INTEGER:: NODE

   ! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF
#ifdef MPI
       IF (NODE==W%WDES%NB_LOW) THEN
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , &
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

                   IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                      CALL M_recv_z(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                      CALL M_recv_d(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
             CALL M_send_z(W%WDES%COMM_INTER,NODE, &
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

             IF (W%WDES%LOVERL) THEN
#ifndef gammareal
                CALL M_send_z(W%WDES%COMM_INTER,NODE, &
                     WTMP%CPROJ(1), W%WDES%NPROD)
#else
                CALL M_send_d(W%WDES%COMM_INTER,NODE, &
                     WTMP%CPROJ(1), W%WDES%NPROD)
#endif
             ENDIF
          ENDIF
       ENDIF
#endif
    ENDDO
    CALL DELWAV(WTMP, .TRUE.)
  END SUBROUTINE W1_GATHER_KNODESEL


!************************ SUBROUTINE W1_GATHER_W1 *********************
!
!> This subroutine gathers from a wavefunction array
!> instead of W
!
!**********************************************************************

  SUBROUTINE W1_GATHER_W1( W, NB2, W1_ORIG, W1)
    !> wavefunction
    TYPE (wavespin):: W
    !> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> final band
    INTEGER:: NB2

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)< NB2) THEN
       CALL vtutor%bug("internal error in W1_GATHER_W1: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    DO N=1,NB2
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=NB2*W%WDES%NB_PAR

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
          ENDIF
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif

  END SUBROUTINE W1_GATHER_W1


!************************ SUBROUTINE W1_GATHER_STRIP ******************
!
!> This subroutine gathers from a wavefunction array instead of W.
!> Collected are the bands W1_ORIG(N): NB1 <= NB_LOW+(N-1)*NB_PAR <= NB2,
!> into W1(1:NB2-NB1+1)
!
!**********************************************************************

  SUBROUTINE W1_GATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    !> wavefunction
    TYPE (wavespin):: W
    !> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> first band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2

  ! local
    INTEGER :: NN, N, NCPU

    PROFILING_START('w1_gather_strip')

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)*NCPU<NB2) THEN
       CALL vtutor%bug("internal error in W1_GATHER_STRIP: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_STRIP: W1 is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    DO N=1,SIZE(W1_ORIG)
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NN<NB1.OR.NN>NB2) CYCLE
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN-NB1+1) )
       CALL FFTWAV_W1( W1(NN-NB1+1))
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=NB1,NB2
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN
#ifndef gammareal
             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
             CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
          ENDIF
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_strip')

  END SUBROUTINE W1_GATHER_STRIP


!************************ SUBROUTINE W1_IGATHER_STRIP *****************
!
!> Gathers a set of wavefunctions starting from band NB1 until NB2
!> (global band indices) to all nodes using non-blocking bcast_from.
!> This subroutine gathers from a wavefunction array W1_ORIG instead
!> of W.
!>
!> In case the code is compiled with-Dshmem_bcast_buffer communication
!> will be between ranks within COMM_inter_node, i.e., only those ranks
!> within COMM_INTER that are NOT on the same physical node will talk to
!> eachother. MPI communication between ranks within COMM_INTER that
!> reside on the same node is not necessary since they access a common
!> shared memory segment.
!>
!> In case the code is NOT compiled with-Dshmem_bcast_buffer, the
!> communication will be between all ranks in COMM_INTER.
!
!**********************************************************************
#ifdef shmem_bcast_buffer
  SUBROUTINE W1_IGATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    !> wavefunction
    TYPE (wavespin):: W
    !> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> first band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2

  ! local
    INTEGER :: NI,N,NB_LOCAL
    INTEGER :: ierror

    INTEGER :: IDO(NB2-NB1+1)

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))


    PROFILING_START('w1_igather_strip')

    IF (SIZE(W1_ORIG)*W%WDES%NB_PAR<NB2) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1 is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    CALL M_barrier(W%WDES%COMM_intra_node)

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY_NOCR(W1_ORIG(NB_LOCAL),W1(NI)); CALL FFTWAV_W1(W1(NI))
          IDO(NI)=W%WDES%COMM_inter_node%NODE_ME
       ELSE
          IDO(NI)=0
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_inter_node%NCPU>1) THEN

       CALLMPI( M_sum_i(W%WDES%COMM_inter_node,IDO,NB2-NB1+1))

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1
          IF (IDO(NI)>0) THEN
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CW(1), &
                  SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CR(1), &
                  SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))

             IF (W%WDES%LOVERL) THEN
                nrequests=nrequests+1
#ifndef gammareal
                CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#else
                CALL M_ibcast_d_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#endif
             ENDIF
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF

    CALL M_barrier(W%WDES%COMM_intra_node)
!!#ifndef shmem_bcast_buffer
!!    IF (W%WDES%COMM_intra_node%NCPU>1) THEN
!!       DO N=1,NB2-NB1+1
!!          IF (IDO(N)>0) IDO(N)=W%WDES%COMM_intra_node%NODE_ME
!!       ENDDO
!!
!!       CALLMPI( M_sum_i(W%WDES%COMM_intra_node,IDO,NB2-NB1+1))
!!
!!       nrequests=0
!!       DO N=NB1,NB2
!!          NI=N-NB1+1
!!          IF (IDO(NI)>0) THEN
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CW(1), &
!!                  SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CR(1), &
!!                  SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
!!             IF (W%WDES%LOVERL) THEN
!!                nrequests=nrequests+1
!!#ifndef gammareal
!!                CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#else
!!                CALL M_ibcast_d_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!                     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#endif
!!             ENDIF
!!          ENDIF
!!       ENDDO
!!
!!       CALL M_waitall(nrequests,requests(1))
!!
!!    ENDIF
!!#endif
#endif
    PROFILING_STOP('w1_igather_strip')

  END SUBROUTINE W1_IGATHER_STRIP
#else
  SUBROUTINE W1_IGATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    !> wavefunction
    TYPE (wavespin):: W
    !> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
    !> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
    !> first band
    INTEGER:: NB1
    !> final band
    INTEGER:: NB2

  ! local
    INTEGER :: NI,N,NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))


    PROFILING_START('w1_igather_strip')

    IF (SIZE(W1_ORIG)*W%WDES%NB_PAR<NB2) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1 is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY_NOCR(W1_ORIG(NB_LOCAL),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CW(1), &
               SIZE(W1(NI)%CW),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
               SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))

          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ibcast_d_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_igather_strip')

  END SUBROUTINE W1_IGATHER_STRIP
#endif

!************************* SUBROUTINE WVREAL_PRECISE *******************
!
!> this subroutine forces the wavefunction to be real at the Gamma-point
!> it is required for the gamma point only mode
!> to avoid that small non real components develop
!> this version is exact and works through an FFT to real space
!> and then forcing the wavefunction to be real
!> the routine is required only if subspace rotations are performed
!> since the subspace rotation routine can not force wavefunctions
!> to become real
!>
!> if LORBITALREAL is .TRUE., the routine also tries to make the
!> orbitals real (in real space) by calculating the orbital
!> ~~~
!>   phi(r) = u_k(r) e^ikr
!> ~~~
!> taking the real part, storing the coefficients back
!> this is possible at special k-points such as points at the BZ boundary
!> ~~~
!> e.g. k=(+-0.5, +-0.5, +-0.5)
!>      k=(+-0.5, +-0.5,  0)
!>      k=(+-0.5,  0  ,  0)
!> ~~~
!>and permutations thereof
!>
!> the routine will not work properly at any other k-points
!> (and is in fact by-passed in this case)
!> ideally no subspace rotation should be called after calling this
!> routine (though orthogonalization is fine)
!
!***********************************************************************

  SUBROUTINE WVREAL_PRECISE(W)

GPACC USE moffload

    USE constant
    USE random_seeded, ONLY: RANE
    IMPLICIT NONE
    TYPE (wavespin) W
  ! local
    INTEGER :: NK, ISP, NB, ISPINOR, N
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    LOGICAL :: LPHASE
    COMPLEX(q) :: CPHASE(W%WDES%GRID%MPLWV)
    COMPLEX(q) :: AVERAGE_PHASE

! for spinors it is not (always) possible to make wave functions real-valued
! (for that reason we need to use the "_ncl" of VASP even for Gamma-only ...)
     IF (W%WDES%LNONCOLLINEAR) RETURN

#ifdef gammareal
!$ACC ENTER DATA CREATE(W1,WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
    !
    ! search for position of G=0 zero element
    !
    DO NK  =1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (ABS(SUM(W%WDES%VKPT(:,NK)*W%WDES%VKPT(:,NK)))< 1E-12) THEN
          CALL SETWDES(W%WDES,WDES1,NK)
          DO ISP=1,W%WDES%ISPIN
             DO NB=1,W%WDES%NBANDS
                CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
                CALL FFTWAV_W1(W1)
                DO ISPINOR =0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV) __IF_ASYNC__
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV),q) &
                          *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                   CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CW(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
#ifdef ACC_OFFLOAD
!$ACC EXIT DATA DELETE(W1) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
    CALL DELETE_TYPED_VAR(WDES1)
#endif
#else
    IF (W%WDES%LORBITALREAL) THEN
!$ACC ENTER DATA CREATE(W1,WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
!
! force orbitals to be real for non-gamma point only method
!
    DO NK  =1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (ABS(MOD(W%WDES%VKPT(1,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(2,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(3,NK)*2+100,1._q))< 1E-6) THEN
       CALL SETPHASE_WVREAL(W%WDES%VKPT(:,NK), W%WDES%GRID, CPHASE, LPHASE)
!$ACC ENTER DATA COPYIN(CPHASE) IF(OFFLOAD_ON .AND. LPHASE) ASYNC(ACC_ASYNC_Q)

       CALL SETWDES(W%WDES,WDES1,NK)
       DO ISP=1,W%WDES%ISPIN
          DO NB=1,W%WDES%NBANDS
             CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
             CALL FFTWAV_W1(W1)
             ! chose random phase in 1st and 4th quadrant, avoid proximity to +-i
             AVERAGE_PHASE=EXP(CITPI*(RANE()-0.5)/2.1)
             CALLMPI (M_bcast_z( WDES1%COMM_INB, AVERAGE_PHASE, 1))

             DO ISPINOR =0,WDES1%NRSPINORS-1
                IF (LPHASE) THEN
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV,CPHASE) &
!$ACC __IF_ASYNC__
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=CONJG(CPHASE(N))*REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*CPHASE(N)*AVERAGE_PHASE,q) &
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                ELSE
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV) __IF_ASYNC__
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*AVERAGE_PHASE,q) &
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                ENDIF

                CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CW(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
             ENDDO
          ENDDO
       ENDDO
!$ACC EXIT DATA DELETE(CPHASE) IF(OFFLOAD_ON .AND. LPHASE) ASYNC(ACC_ASYNC_Q)
       ENDIF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
#ifdef ACC_OFFLOAD
!$ACC EXIT DATA DELETE(W1) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
    CALL DELETE_TYPED_VAR(WDES1)
#endif
    ENDIF
#endif

    RETURN
  END SUBROUTINE WVREAL_PRECISE


  SUBROUTINE SETPHASE_WVREAL(VKPT, GRID, CPHASE, LPHASE)
    USE constant

    IMPLICIT NONE

    REAL(q) :: VKPT(3)
    TYPE (grid_3d) GRID
    COMPLEX(q) :: CPHASE(GRID%MPLWV)
    LOGICAL LPHASE
  ! local
    REAL(q),PARAMETER :: TINY=1E-6_q
    REAL(q) F1, F2, F3
    INTEGER NC, N, IND
    COMPLEX(q) C, CD, CSUM
    CSUM=0
#ifdef gammareal
    LPHASE=.FALSE.
    RETURN
#endif

    IF (ABS(VKPT(1))>TINY .OR. ABS(VKPT(2))>TINY .OR. ABS(VKPT(3))>TINY) THEN
       LPHASE=.TRUE.
       F1=TPI/GRID%NGX*VKPT(1)
       F2=TPI/GRID%NGY*VKPT(2)
       F3=TPI/GRID%NGZ*VKPT(3)

       IF (GRID%RL%NFAST==3) THEN
          CD=EXP(CMPLX(0,F3,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F1*(GRID%RL%I2(NC)-1)+F2*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                C=C*CD
             ENDDO
          ENDDO
       ELSE
          CD=EXP(CMPLX(0,F1,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F2*(GRID%RL%I2(NC)-1)+F3*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                CSUM=CSUM+C
                C=C*CD
             ENDDO
          ENDDO
       ENDIF
       LPHASE=.TRUE.
    ELSE
       LPHASE=.FALSE.
    ENDIF
  END SUBROUTINE SETPHASE_WVREAL

END MODULE wave_high

!***********************************************************************
!
!> assign a W1 structure an initial value from a wavefunction
!> array
!>
!> the plane wave coefficients and the wave function character (optional)
!> must be supplied
!> for performance reason no runtime checking on anything is performed
!
!***********************************************************************

  SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q):: C(*)
    GDEF, OPTIONAL :: CPROJ(*)

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,C, W1%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, C, 1, W1%CW(1), 1)
#endif
    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ENDIF
    ENDIF

  END SUBROUTINE ARRAY_TO_W1


!***********************************************************************
!
!> assign a W1 structure to an array (reverse of the previous operation)
!
!***********************************************************************

  SUBROUTINE W1_TO_ARRAY( W1, C,  CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q) :: C(*)
    GDEF, OPTIONAL :: CPROJ(*)

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 , W1%CW(1), C)
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, C, 1)
#endif
    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ENDIF
    ENDIF

  END SUBROUTINE W1_TO_ARRAY


!************************* SUBROUTINE ECCP_NL **************************
!
!> this subroutine calculates the expectation value of < c|H|cp>
!> where c and cp are two wavefunctions;
!>
!> non local part only for one ion. I have put this in a
!> separate routine because optimization
!> is than easier
!
!***********************************************************************

  SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
GPACC !$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET (ECCP_NL)
    GDEF      CNL
    INTEGER LMDIM, LMMAXC
    OVERLAP CDIJ(LMDIM,LMDIM)
    GDEF CPROJ1(LMMAXC),CPROJ2(LMMAXC)
 ! local
    INTEGER L, LP

!   PROFILING_START('eccp_nl')

__CPU !DIR$ IVDEP
__CPU !OCL NOVREL
GPACC !$ACC LOOP VECTOR COLLAPSE(2) REDUCTION(+:CNL)
GPOMP !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:CNL)
    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+CDIJ(LP,L)*CPROJ1(LP)*GCONJG(CPROJ2(L))
       ENDDO
    ENDDO

!   PROFILING_STOP('eccp_nl')

  END SUBROUTINE ECCP_NL

!***********************************************************************
! subroutines OVERL, OVERL1, OVERL_MU, et al:
!
!  OVERL calculates the action of the overlap operator onto a set of
!   wave functions characters (passed as a contiguous 2d array).
!
!  OVERL1 calculates the action of the Hamiltonian and overlap operator
!   (Dij - e Qij) onto a single wave function character (passed as wavefun1).
!
!  OVERL_MU is a batched version of OVERL1 (array of wavefun1).
!
!  OVERL1_SIMPLE and OVERL_MU_SIMPLE are simplified versions of OVERl1
!   and OVERL_MU where the term  -e Qij is omitted.
!
!  OVERL1_C is a version of OVERL1 for complex eigenvalues.
!  OVERL1_CCDIJ is a version of OVERL1 where CDIJ and CQIJ are always complex.
!
!***********************************************************************

!************************* SUBROUTINE OVERL ***************************
!
!> calculate the result of the overlap-operator acting onto a set of
!> wave function characters; F77 low level routine
!> ~~~
!>  CRESUL^n_N, nlm = sum_n'l'm' Q_N, n'l'm',nlm CPROF^n_N, n'l'm'
!>  CRESUL^n_N, i   = sum_j Q_N, j, i CPROF^n_N, j
!> ~~~
!>
!> n is the band index, N the ion index, nlm = i the index for the
!> one centre partial waves
!>
!> @details @ref openmp :
!> the loop over the bands owned locally by a particular MPI-rank
!> is distributed over all available OpenMP threads.
!
!**********************************************************************

  SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER OVERL
#endif
__GPU USE moffload
    USE wave

    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    LOGICAL LOVERL
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    GDEF CRESUL(WDES1%NPROD,WDES1%NBANDS),CPROF(WDES1%NPROD,WDES1%NBANDS)
  ! local
    GDEF CTMP
    INTEGER NB,NP,ISPINOR,ISPINOR_,NPRO,NPRO_,NPRO2,NPRO2_,NT,NI,LMMAXC,L,LP,NBANDS

    PROFILING_START('overl')

    NBANDS = WDES1%NBANDS

    IF (LOVERL) THEN

GPACC !$ACC KERNELS PRESENT(CRESUL)  __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CRESUL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
       CRESUL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

       spinor1: DO ISPINOR = 0, WDES1%NRSPINORS-1
          spinor2: DO ISPINOR_ = 0, WDES1%NRSPINORS-1
             NPRO  = ISPINOR  * (WDES1%NPRO / 2)
             NPRO_ = ISPINOR_ * (WDES1%NPRO / 2)

__CPU !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
__CPU !$OMP SHARED(WDES1,CRESUL,CQIJ,CPROF,NBANDS,NPRO,NPRO_,ISPINOR,ISPINOR_) &
__CPU !$OMP PRIVATE(NB,NPRO2,NPRO2_,NI,NT,LMMAXC,L,LP,CTMP)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) THREAD_LIMIT(64) &
GPOMP !$OMP PRIVATE(NB,NI,NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP) &
GPOMP !$OMP SHARED(WDES1,CQIJ,CPROF,CRESUL) IF(OFFLOAD_ON)
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,CQIJ,CPROF,CRESUL) COLLAPSE(2) VECTOR_LENGTH(32) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(NI,NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)
             bands: DO NB=1,NBANDS
                ions: DO NI=1,WDES1%NIONS
                   NT=WDES1%ITYP(NI)
                   LMMAXC=WDES1%LMMAX(NT)
                   NPRO2 =WDES1%LMBASE(NI) + NPRO
                   NPRO2_=WDES1%LMBASE(NI) + NPRO_

GPACC !$ACC LOOP VECTOR PRIVATE(CTMP,LP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP,LP) IF(OFFLOAD_ON)
!NEC$ select_vector
                   DO L =1,LMMAXC
#ifndef __NEC__
                      CTMP=0
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
                      DO LP=1,LMMAXC
                         CTMP=CTMP+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO2_,NB)
                      ENDDO

                      CRESUL(L+NPRO2,NB)=CRESUL(L+NPRO2,NB)+CTMP
#else // __NEC__
!NEC$ outerloop_unroll(8)
                      DO LP=1,LMMAXC
                         CRESUL(L+NPRO2,NB)=CRESUL(L+NPRO2,NB)+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO2_,NB)
                      ENDDO
#endif // __NEC__
                   ENDDO
                ENDDO ions
             ENDDO bands
          ENDDO spinor2
       ENDDO spinor1

    ENDIF

    PROFILING_STOP('overl')

  END SUBROUTINE OVERL

!************************* SUBROUTINE OVERL_MU ************************
!
!> This subroutine calculates the action of the Hamiltonian and
!> overlap-operator onto the wave function characters of an array
!> of wave functions (similar to what #OVERL1 does for a single wave
!> function)
!
!**********************************************************************
  SUBROUTINE OVERL_MU(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, W1, NSIM, CRESUL)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER OVERL_MU
#endif
     __GPU USE moffload_struct_def
     USE wave
     USE tutor, ONLY : vtutor
     IMPLICIT NONE

     TYPE(wavedes1) :: WDES1
     TYPE(wavefun1) :: W1(NSIM)
     INTEGER        :: LMDIM, NSIM
     OVERLAP        :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
     OVERLAP        :: CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
     REAL(q)        :: EVALUE(NSIM)
     GDEF           :: CRESUL(WDES1%NPROD,NSIM)

     ! Local variables
     LOGICAL        :: LDO(NSIM)
     INTEGER        :: ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NI, L, LP, NP
     GDEF           :: CTMP

     PROFILING_START('overl_mu')

     !----------------------------------------
     ! Initialize result
     !----------------------------------------
GPACC !$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CRESUL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
     CRESUL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

     LDO = W1%LDO

     !----------------------------------------
     ! Main computation
     !----------------------------------------
     spinor: DO ISPINOR = 0, WDES1%NRSPINORS-1
       DO ISPINOR_ = 0, WDES1%NRSPINORS-1

         NPRO  = ISPINOR  * (WDES1%NPRO / 2)
         NPRO_ = ISPINOR_ * (WDES1%NPRO / 2)

__CPU !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
__CPU !$OMP SHARED(WDES1,CRESUL,CQIJ,CDIJ,EVALUE,LDO,W1,NPRO,NPRO_,ISPINOR,ISPINOR_,NSIM) &
__CPU !$OMP PRIVATE(NP,NPRO2,NPRO2_,NI,L,LP,CTMP)

GPACC !$ACC PARALLEL LOOP GANG COLLAPSE(2) VECTOR_LENGTH(32) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(NPRO2,NPRO2_,CTMP,LP) PRESENT(WDES1,CDIJ,CRESUL,W1,CQIJ,EVALUE)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) THREAD_LIMIT(64) &
GPOMP !$OMP MAP(PRESENT,TO:WDES1,CDIJ,W1,CRESUL,CQIJ,EVALUE) &
GPOMP !$OMP PRIVATE(NPRO2,NPRO2_,CTMP,LP) IF(OFFLOAD_ON)

         DO NP = 1, NSIM
           DO NI = 1, WDES1%NIONS
             IF (.NOT.LDO(NP)) CYCLE
             NPRO2  = WDES1%LMBASE(NI) + NPRO
             NPRO2_ = WDES1%LMBASE(NI) + NPRO_

__CPU !$OMP PARALLEL DO SIMD PRIVATE(CTMP,LP)
GPACC !$ACC LOOP VECTOR PRIVATE(CTMP,LP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP,LP) IF(OFFLOAD_ON)
             DO L = 1, WDES1%LMMAX(WDES1%ITYP(NI))
               CTMP = 0
               DO LP = 1, WDES1%LMMAX(WDES1%ITYP(NI))
                 CTMP = CTMP + (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1) - &
                                EVALUE(NP)*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * &
                                W1(NP)%CPROJ(LP+NPRO2_)
               ENDDO
               CRESUL(L+NPRO2, NP) = CRESUL(L+NPRO2, NP) + CTMP
             ENDDO

           ENDDO
         ENDDO
       ENDDO
     ENDDO spinor

     PROFILING_STOP('overl_mu')

  END SUBROUTINE OVERL_MU

!************************* SUBROUTINE OVERL_MU_SIMPLE ******************
!
!> same as #OVERL_MU but without augmentation charge contributions
!
!**********************************************************************
  SUBROUTINE OVERL_MU_SIMPLE(WDES1, LMDIM, CDIJ, W1, NSIM, CRESUL)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER OVERL_MU_SIMPLE
#endif
     __GPU USE moffload_struct_def
     USE wave
     USE tutor, ONLY : vtutor
     IMPLICIT NONE

     TYPE(wavedes1) :: WDES1
     TYPE(wavefun1) :: W1(NSIM)
     INTEGER        :: LMDIM, NSIM
     OVERLAP        :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
     GDEF           :: CRESUL(WDES1%NPROD,NSIM)

     ! Local variables
     LOGICAL        :: LDO(NSIM)
     INTEGER        :: ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NI, L, LP, NP
     GDEF           :: CTMP

     PROFILING_START('overl_mu_simple')

     !----------------------------------------
     ! Initialize result
     !----------------------------------------
GPACC !$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CRESUL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
     CRESUL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

     LDO = W1%LDO

     !----------------------------------------
     ! Main computation
     !----------------------------------------
     spinor: DO ISPINOR = 0, WDES1%NRSPINORS-1
       DO ISPINOR_ = 0, WDES1%NRSPINORS-1

         NPRO  = ISPINOR  * (WDES1%NPRO / 2)
         NPRO_ = ISPINOR_ * (WDES1%NPRO / 2)

__CPU !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
__CPU !$OMP SHARED(WDES1,CRESUL,CDIJ,LDO,W1,NPRO,NPRO_,ISPINOR,ISPINOR_,NSIM) &
__CPU !$OMP PRIVATE(NP,NPRO2,NPRO2_,NI,L,LP,CTMP)

GPACC !$ACC PARALLEL LOOP GANG COLLAPSE(2) VECTOR_LENGTH(32) __IF_ASYNC__ &
GPACC !$ACC PRIVATE(NPRO2,NPRO2_,CTMP,LP) PRESENT(WDES1,CDIJ,CRESUL,W1)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) THREAD_LIMIT(64) &
GPOMP !$OMP MAP(PRESENT,TO:WDES1,CDIJ,W1,CRESUL) &
GPOMP !$OMP PRIVATE(NPRO2,NPRO2_,CTMP,LP) IF(OFFLOAD_ON)

         DO NP = 1, NSIM
           DO NI = 1, WDES1%NIONS
             IF (.NOT.LDO(NP)) CYCLE
             NPRO2  = WDES1%LMBASE(NI) + NPRO
             NPRO2_ = WDES1%LMBASE(NI) + NPRO_

__CPU !$OMP PARALLEL DO SIMD PRIVATE(CTMP,LP)
GPACC !$ACC LOOP VECTOR PRIVATE(CTMP,LP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP,LP) IF(OFFLOAD_ON)
             DO L = 1, WDES1%LMMAX(WDES1%ITYP(NI))
               CTMP = 0
               DO LP = 1, WDES1%LMMAX(WDES1%ITYP(NI))
                 CTMP = CTMP + CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1) * W1(NP)%CPROJ(LP+NPRO2_)
               ENDDO
               CRESUL(L+NPRO2, NP) = CRESUL(L+NPRO2, NP) + CTMP
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDDO spinor

     PROFILING_STOP('overl_mu_simple')

  END SUBROUTINE OVERL_MU_SIMPLE

!************************* SUBROUTINE OVERL1 ****************************
!
!> calculate the result of the Hamiltonian and overlap-operator acting
!> onto one wave function character;  F77 low level routine
!> ~~~
!>  CRESUL_N, nlm = sum_n'l'm' (D_N, n'l'm',nlm-e Q_N, n'l'm',nlm) CPROF_N, n'l'm'
!>  CRESUL_N, i   = sum_j (D_N, j, i-e Q_N, j, i) CPROF_N, j
!> ~~~
!>
!> N the ion index, nlm = i the index for the one centre partial waves
!>
!> @details @ref openmp :
!> the nested loops over \"types\" + \"ions-of-type\" are replaced
!> by a single loop over \"all ions\" that is distributed over
!> all available threads.
!
!**********************************************************************
  SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER OVERL1
#endif
__GPU USE moffload_struct_def

    USE wave
    IMPLICIT NONE

    TYPE (wavedes1)  :: WDES1
    INTEGER          :: LMDIM
    OVERLAP          :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    OVERLAP          :: CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q)          :: EVALUE
    GDEF             :: CPROF(WDES1%NPRO)
    GDEF             :: CRESUL(WDES1%NPRO)

    ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NI, L, LP
    GDEF CTMP

    PROFILING_START('overl1')

GPACC !$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CRESUL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CRESUL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC ENTER DATA COPYIN(EVALUE) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:EVALUE) IF(OFFLOAD_ON)

    spinor1: DO ISPINOR=0,WDES1%NRSPINORS-1
    spinor2: DO ISPINOR_=0,WDES1%NRSPINORS-1
       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

GPACC !$ACC PARALLEL LOOP GANG __IF_ASYNC__ &
GPACC !$ACC PRESENT(WDES1,CDIJ,CQIJ,CPROF,CRESUL,EVALUE) PRIVATE(NPRO2,NPRO2_,L,LP,CTMP)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WDES1,CDIJ,CQIJ,CPROF,CRESUL,EVALUE) &
GPOMP !$OMP PRIVATE(NPRO2,NPRO2_,L,LP,CTMP) THREAD_LIMIT(64) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP SHARED(NPRO,NPRO_,WDES1,EVALUE,ISPINOR,ISPINOR_,CDIJ,CQIJ,CPROF,CRESUL) &
__CPU !$OMP PRIVATE(NI,NPRO2,NPRO2_,L,LP,CTMP)
!NEC$ novector
       DO NI=1,WDES1%NIONS
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
GPACC !$ACC LOOP VECTOR PRIVATE(CTMP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP,LP) IF(OFFLOAD_ON)
!NEC$ select_vector
             DO L =1,WDES1%LMMAX(WDES1%ITYP(NI))
#ifndef __NEC__
                CTMP=0
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,WDES1%LMMAX(WDES1%ITYP(NI))
                   CTMP=CTMP+ &
                       (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                        EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
                ENDDO
!GPOMP !$OMP END PARALLEL DO
                CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
#else // __NEC__
!NEC$ outerloop_unroll(8)
                DO LP=1,WDES1%LMMAX(WDES1%ITYP(NI))
                   CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+ &
                       (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                        EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
                ENDDO
#endif // __NEC__
             ENDDO
       ENDDO
__CPU !$OMP END PARALLEL DO

    ENDDO spinor2
    ENDDO spinor1

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(EVALUE) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:EVALUE) IF(OFFLOAD_ON)
#endif

    PROFILING_STOP('overl1')

  END SUBROUTINE OVERL1

!************************* SUBROUTINE OVERL1_SIMPLE *******************
!
!> Same as OVERL1 but without augmentation charge contribution
!
!**********************************************************************
  SUBROUTINE OVERL1_SIMPLE(WDES1, LMDIM, CDIJ, CPROF, CRESUL)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER OVERL1_SIMPLE
#endif
__GPU USE moffload_struct_def

    USE wave
    IMPLICIT NONE

    TYPE (wavedes1)  :: WDES1
    INTEGER          :: LMDIM
    OVERLAP          :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    GDEF             :: CPROF(WDES1%NPRO)
    GDEF             :: CRESUL(WDES1%NPRO)

    ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NI, L, LP
    GDEF CTMP

    PROFILING_START('overl1_simple')

GPACC !$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CRESUL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CRESUL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    spinor1: DO ISPINOR=0,WDES1%NRSPINORS-1
    spinor2: DO ISPINOR_=0,WDES1%NRSPINORS-1
       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

       ! =====================================================
       ! Simpler kernel without EVALUE * CQIJ
       ! =====================================================
GPACC !$ACC PARALLEL LOOP GANG __IF_ASYNC__ &
GPACC !$ACC PRESENT(WDES1,CDIJ,CPROF,CRESUL) PRIVATE(NPRO2,NPRO2_,L,LP,CTMP)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WDES1,CDIJ,CPROF,CRESUL) &
GPOMP !$OMP PRIVATE(NPRO2,NPRO2_,L,LP,CTMP) THREAD_LIMIT(64) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP SHARED(NPRO,NPRO_,WDES1,ISPINOR,ISPINOR_,CDIJ,CPROF,CRESUL) &
__CPU !$OMP PRIVATE(NI,NPRO2,NPRO2_,L,LP,CTMP)
!NEC$ novector
       DO NI=1,WDES1%NIONS
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
GPACC !$ACC LOOP VECTOR PRIVATE(CTMP,LP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP,LP) IF(OFFLOAD_ON)
!NEC$ select_vector
             DO L =1,WDES1%LMMAX(WDES1%ITYP(NI))
#ifndef __NEC__
                CTMP=0
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,WDES1%LMMAX(WDES1%ITYP(NI))
                   CTMP=CTMP + &
                        CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO2_)
                ENDDO
                CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
#else // __NEC__
!NEC$ outerloop_unroll(8)
                DO LP=1,WDES1%LMMAX(WDES1%ITYP(NI))
                   CRESUL(L+NPRO2)=CRESUL(L+NPRO2) + &
                        CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO2_)
                ENDDO
#endif // __NEC__
             ENDDO
       ENDDO
__CPU !$OMP END PARALLEL DO

    ENDDO spinor2
    ENDDO spinor1

    PROFILING_STOP('overl1_simple')

  END SUBROUTINE OVERL1_SIMPLE

!************************* SUBROUTINE OVERL1_C ************************
!
!> identical to #OVERL1 but with complex EVALUE
!
!> @details @ref openmp :
!> the nested loops over \"types\" + \"ions-of-type\" are replaced
!> by a single loop over \"all ions\" that is distributed over
!> all available threads.
!
!**********************************************************************
  SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)

GPACC USE moffload_struct_def

    USE wave
    IMPLICIT NONE

    TYPE (wavedes1)  :: WDES1
    INTEGER          :: LMDIM
    OVERLAP          :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    OVERLAP          :: CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q)       :: EVALUE
    GDEF             :: CPROF(WDES1%NPRO)
    GDEF             :: CRESUL(WDES1%NPRO)

  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NI, L, LP
    GDEF CTMP

    PROFILING_START('overl1_c')

!$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
    CRESUL=0
!$ACC END KERNELS

!$ACC ENTER DATA COPYIN(EVALUE) __IF_ASYNC__

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
#ifdef ACC_OFFLOAD
!$ACC PARALLEL LOOP GANG __IF_ASYNC__ &
!$ACC PRESENT(WDES1,CDIJ,CQIJ,CPROF,CRESUL,EVALUE) PRIVATE(NPRO2,NPRO2_,L,LP,CTMP)
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(NPRO,NPRO_,WDES1,EVALUE,ISPINOR,ISPINOR_,CDIJ,CQIJ,CPROF,CRESUL) &
!$OMP PRIVATE(NI,NPRO2,NPRO2_,L,LP,CTMP)
#endif
       DO NI=1,WDES1%NIONS
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
          DO L =1,WDES1%LMMAX(WDES1%ITYP(NI))
             CTMP=0
!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
             DO LP=1,WDES1%LMMAX(WDES1%ITYP(NI))
                CTMP=CTMP+(CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                     EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
             ENDDO
             CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
          ENDDO
       ENDDO
#ifndef ACC_OFFLOAD
!$OMP END PARALLEL DO
#endif
    ENDDO
    ENDDO spinor

!$ACC EXIT DATA DELETE(EVALUE) __IF_ASYNC__

    PROFILING_STOP('overl1_c')

  END SUBROUTINE OVERL1_C


!************************* SUBROUTINE OVERL1_CCDIJ ********************
!
!> Same as #OVERL1 but with CDIJ and CQIJ always complex
!> for the ncl build this subroutine is identical to #OVERL1
!
!**********************************************************************
  SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1)  :: WDES1
    INTEGER          :: LMDIM
    COMPLEX(q)       :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q)       :: CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q)          :: EVALUE
    GDEF             :: CPROF(WDES1%NPRO)
    GDEF             :: CRESUL(WDES1%NPRO)

  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, L, LP

    PROFILING_START('overl1_ccdij')

    CRESUL=0
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

       NIS =1
       DO NT=1,WDES1%NTYP
          IF (WDES1%LMMAX(NT)/=0) THEN
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                IF (EVALUE==0) THEN
                   DO L =1,WDES1%LMMAX(NT)
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,WDES1%LMMAX(NT)
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ELSE
                   DO L =1,WDES1%LMMAX(NT)
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,WDES1%LMMAX(NT)
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                           EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ENDIF
                NPRO = WDES1%LMMAX(NT)+NPRO
                NPRO_= WDES1%LMMAX(NT)+NPRO_
             ENDDO
          ENDIF
          NIS = NIS+WDES1%NITYP(NT)
       ENDDO
    ENDDO
    ENDDO spinor
    PROFILING_STOP('overl1_ccdij')

  END SUBROUTINE OVERL1_CCDIJ
