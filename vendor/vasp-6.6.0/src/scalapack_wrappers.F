#include "symbol.inc"
#ifdef scaLAPACK
!************************* SUBROUTINE WPDTRTRI *************************
!
!>
!
!***********************************************************************
#undef WPDTRTRI
      SUBROUTINE WPDTRTRI(UPLO,DIAG, N,A,IA,JA,DESCA,INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDTRTRI_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO, DIAG

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDTRTRI_OFFLOAD(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PDTRTRI(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPDTRTRI


!************************* SUBROUTINE WPZTRTRI *************************
!
!>
!
!***********************************************************************
#undef WPZTRTRI
      SUBROUTINE WPZTRTRI(UPLO,DIAG, N,A,IA,JA,DESCA,INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZTRTRI_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO, DIAG

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZTRTRI_OFFLOAD(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PZTRTRI(UPLO,DIAG,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPZTRTRI

!************************* SUBROUTINE WPDPOTRF *************************
!
!>
!
!***********************************************************************
#undef WPDPOTRF
      SUBROUTINE WPDPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDPOTRF_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDPOTRF_OFFLOAD(UPLO,N,A,IA,JA,DESCA,INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PDPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPDPOTRF

!************************* SUBROUTINE WPZPOTRF *************************
!
!>
!
!***********************************************************************
#undef WPZPOTRF
      SUBROUTINE WPZPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZPOTRF_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*)
      INTEGER :: DESCA(*)
      INTEGER :: IA, INFO, JA, N
      CHARACTER(1) :: UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZPOTRF_OFFLOAD(UPLO,N,A,IA,JA,DESCA,INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PZPOTRF(UPLO,N,A,IA,JA,DESCA,INFO)

      RETURN
      END SUBROUTINE WPZPOTRF


!************************* SUBROUTINE WPSSYEVX *************************
!
!>
!
!***********************************************************************
#undef WPSSYEVX
      SUBROUTINE WPSSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                          IFAIL, ICLUSTR, GAP, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PSSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), GAP(*), W(*), WORK(*), Z(*)
      REAL(qs) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PSSYEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
GPACC                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
GPACC                         Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
GPACC                         IFAIL, ICLUSTR, GAP, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PSSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                   IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPSSYEVX


!************************* SUBROUTINE WPDSYEVX *************************
!
!>
!
!***********************************************************************
#undef WPDSYEVX
      SUBROUTINE WPDSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                          IFAIL, ICLUSTR, GAP, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), GAP(*), W(*), WORK(*), Z(*)
      REAL(q) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDSYEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
GPACC                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
GPACC                         Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
GPACC                         IFAIL, ICLUSTR, GAP, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PDSYEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, &
                   IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPDSYEVX


!************************* SUBROUTINE WPCHEEVX *************************
!
!>
!
!***********************************************************************
#undef WPCHEEVX
      SUBROUTINE WPCHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                          IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PCHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), WORK(*), Z(*)
      REAL(qs) :: GAP(*), RWORK(*), W(*)
      REAL(qs) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PCHEEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
GPACC                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
GPACC                         Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
GPACC                         IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PCHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                   IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPCHEEVX


!************************* SUBROUTINE WPZHEEVX *************************
!
!>
!
!***********************************************************************
#undef WPZHEEVX
      SUBROUTINE WPZHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                          VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                          Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                          IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: GAP(*), RWORK(*), W(*)
      REAL(q) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: DESCA(*), DESCZ(*), ICLUSTR(*), IFAIL(*), IWORK(*)
      INTEGER :: IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
      CHARACTER(1) :: JOBZ, RANGE, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZHEEVX_OFFLOAD(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
GPACC                         VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
GPACC                         Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
GPACC                         IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PZHEEVX(JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, &
                   VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, &
                   Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, &
                   IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO)

      RETURN
      END SUBROUTINE WPZHEEVX


!************************* SUBROUTINE WPSSYEVD *************************
!
!>
!
!***********************************************************************
#undef WPSSYEVD
      SUBROUTINE WPSSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                           WORK, LWORK, IWORK, LIWORK, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PSSYEVD_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), W(*), WORK(*), Z(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PSSYEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
GPACC                         WORK, LWORK, IWORK, LIWORK, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PSSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                   WORK, LWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPSSYEVD


!************************* SUBROUTINE WPDSYEVD *************************
!
!>
!
!***********************************************************************
#undef WPDSYEVD
      SUBROUTINE WPDSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                           WORK, LWORK, IWORK, LIWORK, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDSYEVD_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), W(*), WORK(*), Z(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDSYEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
GPACC                         WORK, LWORK, IWORK, LIWORK, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PDSYEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ,&
                   WORK, LWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPDSYEVD


!************************* SUBROUTINE WPCHEEVD *************************
!
!>
!
!***********************************************************************
#undef WPCHEEVD
      SUBROUTINE WPCHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PCHEEVD_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), WORK(*), Z(*)
      REAL(qs) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PCHEEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
GPACC                         WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PCHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                   WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPCHEEVD


!************************* SUBROUTINE WPZHEEVD *************************
!
!>
!
!***********************************************************************
#undef WPZHEEVD
      SUBROUTINE WPZHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZHEEVD_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*), IWORK(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZHEEVD_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
GPACC                         WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PZHEEVD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                   WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

      RETURN
      END SUBROUTINE WPZHEEVD

!************************* SUBROUTINE WPDSYEV **************************
!
!>
!
!***********************************************************************
#undef WPDSYEV
      SUBROUTINE WPDSYEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                         WORK, LWORK, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDSYEV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: W(*)
      INTEGER :: DESCA(*), DESCZ(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDSYEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
GPACC                        WORK, LWORK, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PDSYEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                  WORK, LWORK, INFO)

      RETURN
      END SUBROUTINE WPDSYEV


!************************* SUBROUTINE WPZHEEV **************************
!
!>
!
!***********************************************************************
#undef WPZHEEV
      SUBROUTINE WPZHEEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                         WORK, LWORK, RWORK, LRWORK, INFO)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZHEEV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), WORK(*), Z(*)
      REAL(q) :: RWORK(*), W(*)
      INTEGER :: DESCA(*), DESCZ(*)
      INTEGER :: IA, INFO, IZ, JA, JZ, LRWORK, LWORK, N
      CHARACTER(1) :: JOBZ, UPLO

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZHEEV_OFFLOAD(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
GPACC                        WORK, LWORK, RWORK, LRWORK, INFO)
GPACC    RETURN
GPACC ENDIF

      CALL PZHEEV(JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, &
                  WORK, LWORK, RWORK, LRWORK, INFO)

      RETURN
      END SUBROUTINE WPZHEEV


!************************* SUBROUTINE WPDTRMR2D ************************
!
!>
!
!***********************************************************************
#undef WPDTRMR2D
      SUBROUTINE WPDTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDTRMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PDTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPDTRMR2D


!************************* SUBROUTINE WPSTRMR2D ************************
!
!>
!
!***********************************************************************
#undef WPSTRMR2D
      SUBROUTINE WPSTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PSTRMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PSTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PSTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPSTRMR2D


!************************* SUBROUTINE WPCTRMR2D ************************
!
!>
!
!***********************************************************************
#undef WPCTRMR2D
      SUBROUTINE WPCTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PCTRMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PCTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PCTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPCTRMR2D


!************************* SUBROUTINE WPZTRMR2D ************************
!
!>
!
!***********************************************************************
#undef WPZTRMR2D
      SUBROUTINE WPZTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZTRMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT
      CHARACTER(1) :: UPLO, DIAG

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZTRMR2D_OFFLOAD(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PZTRMR2D(UPLO, DIAG, M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPZTRMR2D


!************************* SUBROUTINE WPSGEMR2D ************************
!
!>
!
!***********************************************************************
#undef WPSGEMR2D
      SUBROUTINE WPSGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PSGEMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PSGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PSGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPSGEMR2D


!************************* SUBROUTINE WPDGEMR2D ************************
!
!>
!
!***********************************************************************
#undef WPDGEMR2D
      SUBROUTINE WPDGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PDGEMR2D_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PDGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PDGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPDGEMR2D


!************************* SUBROUTINE WPCGEMR2D ************************
!
!>
!
!***********************************************************************
#undef WPCGEMR2D
      SUBROUTINE WPCGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PCGEMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PCGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PCGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPCGEMR2D


!************************* SUBROUTINE WPZGEMR2D ************************
!
!>
!
!***********************************************************************
#undef WPZGEMR2D
      SUBROUTINE WPZGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      USE prec
GPACC USE moffload_struct_def
GPACC USE moffload_scalapack, ONLY : PZGEMR2D_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(*), B(*)
      INTEGER :: DESCA(*), DESCB(*)
      INTEGER :: IA, JA, IB, JB, M, N, ICTXT

GPACC IF (OFFLOAD_ON) THEN
GPACC    CALL PZGEMR2D_OFFLOAD(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)
GPACC    RETURN
GPACC ENDIF

      CALL PZGEMR2D(M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICTXT)

      RETURN
      END SUBROUTINE WPZGEMR2D

#endif // scaLAPACK
