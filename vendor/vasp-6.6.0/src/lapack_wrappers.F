#include "symbol.inc"

!***********************************************************************
! LAPACK
!***********************************************************************

!************************* SUBROUTINE WDPOTRF **************************
!
!> Wrapper for the LAPACK routine DPOTRF. (Cholesky factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************
#undef WDPOTRF
      SUBROUTINE WDPOTRF(UPLO,N,A,LDA,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DPOTRF_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DPOTRF(UPLO,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WDPOTRF


!************************* SUBROUTINE WZPOTRF **************************
!
!> Wrapper for the LAPACK routine ZPOTRF. (Cholesky factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************
#undef WZPOTRF
      SUBROUTINE WZPOTRF(UPLO,N,A,LDA,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZPOTRF_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZPOTRF_OFFLOAD(UPLO,N,A,LDA,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZPOTRF(UPLO,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WZPOTRF


!************************* SUBROUTINE WDGETRF **************************
!
!> Wrapper for the LAPACK routine DGETRF. (LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************
#undef WDGETRF
      SUBROUTINE WDGETRF(M,N,A,LDA,IPIV,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DGETRF_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*)
      INTEGER :: M,N,LDA,INFO,IPIV(*)

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DGETRF(M,N,A,LDA,IPIV,INFO)

      RETURN
      END SUBROUTINE WDGETRF


!************************* SUBROUTINE WZGETRF **************************
!
!> Wrapper for the LAPACK routine ZGETRF. (LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************
#undef WZGETRF
      SUBROUTINE WZGETRF(M,N,A,LDA,IPIV,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZGETRF_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*)
      INTEGER :: M,N,LDA,INFO,IPIV(*)

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZGETRF_OFFLOAD(M,N,A,LDA,IPIV,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZGETRF(M,N,A,LDA,IPIV,INFO)

      RETURN
      END SUBROUTINE WZGETRF


!************************* SUBROUTINE WDGETRS **************************
!
!> Wrapper for the LAPACK routine DGETRS. (Solves a system of linear
!> equations with a general matrix using the LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************
#undef WDGETRS
      SUBROUTINE WDGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DGETRS_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
      CHARACTER(1) :: TRANS

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      RETURN
      END SUBROUTINE WDGETRS


!************************* SUBROUTINE WZGETRS **************************
!
!> Wrapper for the LAPACK routine ZGETRS. (Solves a system of linear
!> equations with a general matrix using the LU factorization)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!> 3) INFO = INFO_ * 100 (+ error codes returned by the LAPACK routine
!>    if applicable in the implementation).
!
!***********************************************************************
#undef WZGETRS
      SUBROUTINE WZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZGETRS_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*)
      INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
      CHARACTER(1) :: TRANS

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZGETRS_OFFLOAD(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)

      RETURN
      END SUBROUTINE WZGETRS


!************************* SUBROUTINE WDTRTRI **************************
!
!> Wrapper for the LAPACK routine DTRTRI. (Inverts a triangular matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!
!***********************************************************************
#undef WDTRTRI
      SUBROUTINE WDTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DTRTRI_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO,DIAG

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WDTRTRI


!************************* SUBROUTINE WZTRTRI **************************
!
!> Wrapper for the LAPACK routine ZTRTRI. (Inverts a triangular matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!
!***********************************************************************
#undef WZTRTRI
      SUBROUTINE WZTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZTRTRI_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*)
      INTEGER :: N,LDA,INFO
      CHARACTER(1) :: UPLO,DIAG

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZTRTRI_OFFLOAD(UPLO,DIAG,N,A,LDA,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)

      RETURN
      END SUBROUTINE WZTRTRI


!************************* SUBROUTINE WDSYGV ***************************
!
!> Wrapper for the LAPACK routine DSYGV. (Solves the generalized
!> eigenvalue problem for a real symmetric-definite matrix pair)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************
#undef WDSYGV
      SUBROUTINE WDSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DSYGV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),B(LDB,*),W(*),WORK(*)
      INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DSYGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)

      RETURN
      END SUBROUTINE WDSYGV


!************************* SUBROUTINE WZHEGV ***************************
!
!> Wrapper for the LAPACK routine ZHEGV. (Solves the generalized
!> eigenvalue problem for a complex Hermitian-definite matrix pair)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and B must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************
#undef WZHEGV
      SUBROUTINE WZHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZHEGV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),B(LDB,*),WORK(*)
      REAL(q) :: W(*),RWORK(*)
      INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZHEGV_OFFLOAD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)

      RETURN
      END SUBROUTINE WZHEGV


!************************* SUBROUTINE WDSYEV ***************************
!
!> Wrapper for the LAPACK routine DSYEV. (Computes all eigenvalues and,
!> optionally, eigenvectors of a real symmetric matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************
#undef WDSYEV
      SUBROUTINE WDSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DSYEV_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),W(*),WORK(*)
      INTEGER :: N,LDA,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DSYEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)

      RETURN
      END SUBROUTINE WDSYEV


!************************* SUBROUTINE WZHEEV ***************************
!
!> Wrapper for the LAPACK routine ZHEEV. (Computes all eigenvalues and,
!> optionally, eigenvectors of a complex Hermitian matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) and are copied back
!>    to the host if W was not on the device before the call.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************
#undef WZHEEV
      SUBROUTINE WZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZHEEV_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),WORK(*)
      REAL(q) :: W(*),RWORK(*)
      INTEGER :: N,LDA,LWORK,INFO
      CHARACTER(1) :: JOBZ,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZHEEV_OFFLOAD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)

      RETURN
      END SUBROUTINE WZHEEV


!************************* SUBROUTINE WDSYEVX **************************
!
!> Wrapper for the LAPACK routine DSYEVX. (Computes selected eigenvalues
!> and, optionally, eigenvectors of a real symmetric matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************
#undef WDSYEVX
      SUBROUTINE WDSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : DSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(q) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*),VL,VU,ABSTOL
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL DSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL DSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WDSYEVX


!************************* SUBROUTINE WCHEEVX **************************
!
!>
!
!***********************************************************************
#undef WCHEEVX
      SUBROUTINE WCHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : CHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(qs) :: A(LDA,*),Z(LDZ,*),WORK(*)
      REAL(qs) :: W(*),RWORK(*),VL,VU,ABSTOL
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL CHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL CHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WCHEEVX


!************************* SUBROUTINE WZHEEVX **************************
!
!> Wrapper for the LAPACK routine ZHEEVX. (Computes selected eigenvalues
!> and, optionally, eigenvectors of a complex Hermitian matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrices A and Z must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.
!
!***********************************************************************
#undef WZHEEVX
      SUBROUTINE WZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : ZHEEVX_OFFLOAD

      IMPLICIT NONE

      COMPLEX(q) :: A(LDA,*),Z(LDZ,*),WORK(*)
      REAL(q) :: W(*),RWORK(*),VL,VU,ABSTOL
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL ZHEEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL ZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WZHEEVX


!************************* SUBROUTINE WSSYEVX **************************
!
!> Wrapper for the LAPACK routine SSYEVX. (Computes selected eigenvalues
!> and, optionally, eigenvectors of a real symmetric matrix)
!
!> For GPU offloading, the following conditions must be satisfied in the
!> implementation:
!> 1) Matrix A must be available on the device.
!> 2) Implement a workspace query within the implementation to determine
!>    the optimal workspace size for the LAPACK routine.
!> 3) Array W(1:N) must be created on the device.
!> 4) W (contains eigenvalues after LAPACK call) must be copied back
!>    to the host.
!> 5) INFO_ must be declared, created on the device, used in the
!>    subroutine call, and copied back to the host.      
!
!***********************************************************************
#undef WSSYEVX
      SUBROUTINE WSSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      USE prec
__GPU USE moffload_struct_def
__GPU USE moffload_lapack, ONLY : SSYEVX_OFFLOAD

      IMPLICIT NONE

      REAL(qs) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*)
      REAL(qs) :: ABSTOL,VL,VU
      INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
      CHARACTER(1) :: JOBZ,RANGE,UPLO

__GPU IF (OFFLOAD_ON) THEN
__GPU    CALL SSYEVX_OFFLOAD(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
__GPU    RETURN
__GPU ENDIF

      CALL SSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)

      RETURN
      END SUBROUTINE WSSYEVX


!***********************************************************************
! scaLAPACK
!***********************************************************************
