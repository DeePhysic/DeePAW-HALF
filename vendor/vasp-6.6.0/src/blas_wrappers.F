#include "symbol.inc"

!************************* SUBROUTINE WDCOPY ***************************
!
!> Wrapper for BLAS routine DCOPY. (DY := DX)
!> DX and DY must be present on the device if offloading is enabled. 
!
!***********************************************************************
#undef WDCOPY
      SUBROUTINE WDCOPY(N,DX,INCX,DY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DCOPY_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL DCOPY(N,DX,INCX,DY,INCY)

      RETURN
      END SUBROUTINE WDCOPY


!************************* SUBROUTINE WSCOPY ***************************
!
!> Wrapper for BLAS routine SCOPY. (DY := DX)
!> DX and DY must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WSCOPY
      SUBROUTINE WSCOPY(N,DX,INCX,DY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : SCOPY_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL SCOPY_OFFLOAD(N,DX,INCX,DY,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL SCOPY(N,DX,INCX,DY,INCY)

      RETURN
      END SUBROUTINE WSCOPY


!************************* SUBROUTINE WCCOPY ***************************
!
!> Wrapper for BLAS routine CCOPY. (ZY := ZX)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WCCOPY
      SUBROUTINE WCCOPY(N,ZX,INCX,ZY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : CCOPY_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL CCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL CCOPY(N,ZX,INCX,ZY,INCY)

      RETURN
      END SUBROUTINE WCCOPY


!************************* SUBROUTINE WZCOPY ***************************
!
!> Wrapper for BLAS routine ZCOPY. (ZY := ZX)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WZCOPY
      SUBROUTINE WZCOPY(N,ZX,INCX,ZY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZCOPY_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL ZCOPY(N,ZX,INCX,ZY,INCY)

      RETURN
      END SUBROUTINE WZCOPY


!************************* SUBROUTINE WDAXPY ***************************
!
!> Wrapper for BLAS routine DAXPY. (DY := DA*DX + DY)
!> DX and DY must be present on the device if offloading is enabled.
!> The scalar DA may be passed from the host.
!
!***********************************************************************
#undef WDAXPY
      SUBROUTINE WDAXPY(N,DA,DX,INCX,DY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DAXPY_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*),DY(*)
      REAL(q) :: DA
      INTEGER :: N,INCX,INCY

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DAXPY_OFFLOAD(N,DA,DX,INCX,DY,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL DAXPY(N,DA,DX,INCX,DY,INCY)

      RETURN
      END SUBROUTINE WDAXPY


!************************* SUBROUTINE WZAXPY ***************************
!
!> Wrapper for BLAS routine ZAXPY. (ZY := ZA*ZX + ZY)
!> ZX and ZY must be present on the device if offloading is enabled.
!> The scalar ZA may be passed from the host.
!
!***********************************************************************
#undef WZAXPY
      SUBROUTINE WZAXPY(N,ZA,ZX,INCX,ZY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZAXPY_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*),ZY(*)
      COMPLEX(q) :: ZA
      INTEGER :: N,INCX,INCY

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZAXPY_OFFLOAD(N,ZA,ZX,INCX,ZY,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL ZAXPY(N,ZA,ZX,INCX,ZY,INCY)

      RETURN
      END SUBROUTINE WZAXPY


!************************* SUBROUTINE WDSCAL ***************************
!
!> Wrapper for BLAS routine DSCAL. (DX := DA*DX)
!> DX must be present on the device if offloading is enabled.
!> The scalar DA may be passed from the host.
!
!***********************************************************************
#undef WDSCAL
      SUBROUTINE WDSCAL(N,DA,DX,INCX)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DSCAL_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*)
      REAL(q) :: DA
      INTEGER :: N,INCX

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DSCAL_OFFLOAD(N,DA,DX,INCX)
__GPU    RETURN
__GPU ENDIF

      CALL DSCAL(N,DA,DX,INCX)

      RETURN
      END SUBROUTINE WDSCAL


!************************* SUBROUTINE WZSCAL ***************************
!
!> Wrapper for BLAS routine ZSCAL. (ZX := ZA*ZX)
!> ZX must be present on the device if offloading is enabled.
!> The scalar ZA may be passed from the host.
!
!***********************************************************************
#undef WZSCAL
      SUBROUTINE WZSCAL(N,ZA,ZX,INCX)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZSCAL_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*)
      COMPLEX(q) :: ZA
      INTEGER :: N,INCX

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZSCAL_OFFLOAD(N,ZA,ZX,INCX)
__GPU    RETURN
__GPU ENDIF

      CALL ZSCAL(N,ZA,ZX,INCX)

      RETURN
      END SUBROUTINE WZSCAL


!************************* SUBROUTINE WZDSCAL **************************
!
!> Wrapper for BLAS routine ZDSCAL. (ZX := DA*ZX)
!> ZX must be present on the device if offloading is enabled.
!> The scalar DA may be passed from the host.
!
!***********************************************************************
#undef WZDSCAL
      SUBROUTINE WZDSCAL(N,DA,ZX,INCX)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZDSCAL_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*)
      REAL(q) :: DA
      INTEGER :: N,INCX

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZDSCAL_OFFLOAD(N,DA,ZX,INCX)
__GPU    RETURN
__GPU ENDIF

      CALL ZDSCAL(N,DA,ZX,INCX)

      RETURN
      END SUBROUTINE WZDSCAL


!************************* SUBROUTINE WDDOT ****************************
!
!> Wrapper for BLAS function DDOT. (WDDOT := DX . DY)
!> DX and DY must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WDDOT
      REAL(q) FUNCTION WDDOT(N,DX,INCX,DY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DDOT_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: DX(*),DY(*)
      INTEGER :: N,INCX,INCY
      ! local variables
      REAL(q), EXTERNAL :: DDOT

__GPU IF (OFFLOAD_ON) THEN
__GPU    WDDOT = DDOT_OFFLOAD(N,DX,INCX,DY,INCY)
__GPU    RETURN
__GPU ENDIF

      WDDOT = DDOT(N,DX,INCX,DY,INCY)

      RETURN
      END FUNCTION WDDOT


!************************* SUBROUTINE WZDOTC ***************************
!
!> Wrapper for BLAS function ZDOTC. (WZDOTC := ZX . ZY)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WZDOTC
      COMPLEX(q) FUNCTION WZDOTC(N,ZX,INCX,ZY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZDOTC_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY
      ! local variables
      COMPLEX(q), EXTERNAL :: ZDOTC

__GPU IF (OFFLOAD_ON) THEN
__GPU    WZDOTC = ZDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
__GPU    RETURN
__GPU ENDIF

      WZDOTC = ZDOTC(N,ZX,INCX,ZY,INCY)

      RETURN
      END FUNCTION WZDOTC


!************************* SUBROUTINE WCDOTC ***************************
!
!> Wrapper for BLAS function CDOTC. (WCDOTC := ZX . ZY)
!> ZX and ZY must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WCDOTC
      COMPLEX(qs) FUNCTION WCDOTC(N,ZX,INCX,ZY,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : CDOTC_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: ZX(*),ZY(*)
      INTEGER :: N,INCX,INCY
      ! local variables
      COMPLEX(qs), EXTERNAL :: CDOTC

__GPU IF (OFFLOAD_ON) THEN
__GPU    WCDOTC = CDOTC_OFFLOAD(N,ZX,INCX,ZY,INCY)
__GPU    RETURN
__GPU ENDIF

      WCDOTC = CDOTC(N,ZX,INCX,ZY,INCY)

      RETURN
      END FUNCTION WCDOTC


!************************* SUBROUTINE WDGEMV ***************************
!
!> Wrapper for BLAS routine DGEMV. 
!> (Y := ALPHA*A*X + BETA*Y, if MODE='N';
!> Y := ALPHA*A'*X + BETA*Y, if MODE='T';
!> Y := ALPHA*conjg(A')*X + BETA*Y, if MODE='C')
!>
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> A and X must be present on the device.
!> Y must be copied to the device if BETA /= 0.
!> Y must be created on the device if BETA = 0.
!> Y must be copied to the host at the end.
!
!***********************************************************************
#undef WDGEMV
      SUBROUTINE WDGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DGEMV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),X(*),Y(*)
      REAL(q) :: ALPHA,BETA
      INTEGER :: M,N,LDA,INCX,INCY
      CHARACTER(1) :: MODE

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL DGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      RETURN
      END SUBROUTINE WDGEMV


!************************* SUBROUTINE WZGEMV ***************************
!
!> Wrapper for BLAS routine ZGEMV.
!> (Y := ALPHA*A*X + BETA*Y, if MODE='N';
!> Y := ALPHA*A'*X + BETA*Y, if MODE='T';
!> Y := ALPHA*conjg(A')*X + BETA*Y, if MODE='C')
!>
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> A and X must be present on the device.
!> Y must be copied to the device if BETA /= 0.
!> Y must be created on the device if BETA = 0.
!> Y must be copied to the host at the end.
!
!***********************************************************************
#undef WZGEMV
      SUBROUTINE WZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZGEMV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),X(*),Y(*)
      COMPLEX(q) :: ALPHA,BETA
      INTEGER :: M,N,LDA,INCX,INCY
      CHARACTER(1) :: MODE

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZGEMV_OFFLOAD(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
__GPU    RETURN
__GPU ENDIF

      CALL ZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      RETURN
      END SUBROUTINE WZGEMV


!************************* SUBROUTINE WDGEMM ***************************
!
!> Wrapper for BLAS routine DGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WDGEMM
      SUBROUTINE WDGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DGEMM_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
      REAL(q) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
__GPU    RETURN
__GPU ENDIF

      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WDGEMM


!************************* SUBROUTINE WSGEMM ***************************
!
!> Wrapper for BLAS routine SGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WSGEMM
      SUBROUTINE WSGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : SGEMM_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
      REAL(qs) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
__GPU    RETURN
__GPU ENDIF

      CALL SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WSGEMM


!************************* SUBROUTINE WZGEMM ***************************
!
!> Wrapper for BLAS routine ZGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WZGEMM
      SUBROUTINE WZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZGEMM_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(q) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
__GPU    RETURN
__GPU ENDIF

      CALL ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WZGEMM


!************************* SUBROUTINE WCGEMM ***************************
!
!> Wrapper for BLAS routine CGEMM. (C := ALPHA*op(A)*op(B) + BETA*C)
!>    if TRANSA/B = 'N' or 'n', then op(A/B) = A;
!>    if TRANSA/B = 'T' or 't', then op(A/B) = A/B';
!>    if TRANSA/B = 'C' or 'c', then op(A/B) = conjg(A/B');
!> A, B and C must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WCGEMM
      SUBROUTINE WCGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : CGEMM_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(qs) :: ALPHA,BETA
      INTEGER :: M,N,K,LDA,LDB,LDC
      CHARACTER(1) :: TRANSA,TRANSB

      IF (N<1) RETURN

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL CGEMM_OFFLOAD(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
__GPU    RETURN
__GPU ENDIF

      CALL CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      RETURN
      END SUBROUTINE WCGEMM


!************************* SUBROUTINE WDTRMM ***************************
!
!> Wrapper for BLAS routine DTRMM.
!>    if SIDE = 'L' or 'l', then B := ALPHA*op(A)*B;
!>    if SIDE = 'R' or 'r', then B := ALPHA*B*op(A);
!>    if TRANSA = 'N' or 'n', then op(A) = A;
!>    if TRANSA = 'T' or 't', then op(A) = A';
!>    if TRANSA = 'C' or 'c', then op(A) = conjg(A');
!> A and B must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WDTRMM
      SUBROUTINE WDTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DTRMM_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*)
      REAL(q) :: ALPHA
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
__GPU    RETURN
__GPU ENDIF

      CALL DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      RETURN
      END SUBROUTINE WDTRMM


!************************* SUBROUTINE WZTRMM ***************************
!
!> Wrapper for BLAS routine ZTRMM.
!>    if SIDE = 'L' or 'l', then B := ALPHA*op(A)*B;
!>    if SIDE = 'R' or 'r', then B := ALPHA*B*op(A);
!>    if TRANSA = 'N' or 'n', then op(A) = A;
!>    if TRANSA = 'T' or 't', then op(A) = A';
!>    if TRANSA = 'C' or 'c', then op(A) = conjg(A');
!> A and B must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WZTRMM
      SUBROUTINE WZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : ZTRMM_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*)
      COMPLEX(q) :: ALPHA
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZTRMM_OFFLOAD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
__GPU    RETURN
__GPU ENDIF

      CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      RETURN
      END SUBROUTINE WZTRMM


!************************* SUBROUTINE WDLACPY **************************
!
!> Wrapper for BLAS routine DLACPY. (B := A)
!> Only a full copy via dcopy is implemented for GPU offloading.
!> This means UPLO = F is forced!
!> A and B must be present on the device if offloading is enabled.
!
!***********************************************************************
#undef WDLACPY
      SUBROUTINE WDLACPY(UPLO,M,N,A,LDA,B,LDB)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_blas, ONLY : DLACPY_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: LDA,LDB,M,N
      CHARACTER(1) :: UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DLACPY_OFFLOAD(UPLO,M,N,A,LDA,B,LDB)
__GPU    RETURN
__GPU ENDIF

      CALL DLACPY(UPLO,M,N,A,LDA,B,LDB)

      RETURN
      END SUBROUTINE WDLACPY
