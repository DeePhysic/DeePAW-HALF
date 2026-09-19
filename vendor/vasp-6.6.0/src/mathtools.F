#include "symbol.inc"
MODULE mathtools
   USE prec
   USE string, ONLY: str
   USE tutor,  ONLY: vtutor
   IMPLICIT NONE 

   PRIVATE

   PUBLIC INVERT_REAL_MATRIX, det_3x3, matrix_inverse_3x3, det_3x3i, matrix_inverse_3x3i, SV_DECOMPOSER, RRQR_DECOMPOSER, erf_pop, erf_gaussian_pop

!                                _______________
!_______________________________/ SV_DECOMPOSER \_______________________________
!
!> @brief Wrapper for DGESVD/ZGESVD LAPACK routines
!>
!> This type facilitates the functionality for performing a singular-value
!> decomposition multiple times without reallocating any of the buffers or work
!> arrays (for example inside a loop).
!>
!> The public type variables are used as storage for the output arrays.
!>
!> @var U_MAT           Left matrix in the SVD
!> @var V_DAGGER_MAT    Right matrix in the SVD
!> @var SINGULAR_VALS   Singular values
!_______________________________________________________________________________
   TYPE SV_DECOMPOSER
      GDEF, ALLOCATABLE,      PUBLIC   :: U_MAT(:, :)
      GDEF, ALLOCATABLE,      PUBLIC   :: V_DAGGER_MAT(:, :)
      REAL(q), ALLOCATABLE,   PUBLIC   :: SINGULAR_VALS(:)

      CHARACTER(1),  PRIVATE  :: JOBU  = " "
      CHARACTER(1),  PRIVATE  :: JOBVT = " "
      INTEGER,       PRIVATE  :: DIM_1 = 0
      INTEGER,       PRIVATE  :: DIM_2 = 0
      INTEGER,       PRIVATE  :: MIN_DIM = 0
      INTEGER,       PRIVATE  :: LDU = 0
      INTEGER,       PRIVATE  :: LDVD = 0

      INTEGER,                PRIVATE  :: WORK_SIZE
      GDEF, ALLOCATABLE,      PRIVATE  :: WORK(:)
      REAL(q), ALLOCATABLE,   PRIVATE  :: RWORK(:)

      CONTAINS

      PROCEDURE,  PUBLIC   :: INIT => SV_DECOMPOSER_INIT
      PROCEDURE,  PUBLIC   :: CALC => SV_DECOMPOSER_CALC
      PROCEDURE,  PUBLIC   :: DEALLOC => SV_DECOMPOSER_DEALLOC
   END TYPE


!                               _________________
!______________________________/ RRQR_DECOMPOSER \______________________________
!
!> @brief Wrapper for xGEQP3 LAPACK routines
!>
!> This type facilitates the functionality for performing a rank-revealing QR
!> decomposition multiple times without reallocating any of the buffers or work
!> arrays (for example inside a loop).
!>
!> The public type variables are used as storage for the output arrays.
!>
!> @var PIVOT  Pivot values
!> @var TAU    Scalar factors of the elementary reflectors
!_______________________________________________________________________________
   TYPE RRQR_DECOMPOSER
      INTEGER, ALLOCATABLE,   PUBLIC   :: PIVOT(:)
      GDEF, ALLOCATABLE,      PUBLIC   :: TAU(:)

      INTEGER,       PRIVATE  :: DIM_1 = 0
      INTEGER,       PRIVATE  :: DIM_2 = 0
      INTEGER,       PRIVATE  :: LDA = 0

      INTEGER,                PRIVATE  :: WORK_SIZE
      GDEF, ALLOCATABLE,      PRIVATE  :: WORK(:)
      REAL(q), ALLOCATABLE,   PRIVATE  :: RWORK(:)

      CONTAINS

      PROCEDURE,  PUBLIC   :: INIT => RRQR_DECOMPOSER_INIT
      PROCEDURE,  PUBLIC   :: CALC => RRQR_DECOMPOSER_CALC
      PROCEDURE,  PUBLIC   :: DEALLOC => RRQR_DECOMPOSER_DEALLOC
   END TYPE


   CONTAINS

!                             ____________________
!____________________________/ SV_DECOMPOSER_INIT \_____________________________
!
!> @brief Class constructor
!>
!> Internal buffers and work arrays are set up based on the supplied matrix
!> dimensions of a matrix A that is about to be singular-value decomposed by
!> this calculator.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in] JOBU      Same as `JOBU` argument to DGESVD or ZGESVD
!> @param[in] JOBVT     Same as `JOBVT` argument to DGESVD or ZGESVD
!> @param[in] DIM_1     SIZE(A, 1)
!> @param[in] DIM_2     SIZE(A, 2)
!> @param[out] INFO     `INFO` parameter of DGESVD/ZGESVD query
!_______________________________________________________________________________
   SUBROUTINE SV_DECOMPOSER_INIT(THIS, JOBU, JOBVT, DIM_1, DIM_2, INFO)
      USE string, ONLY: uppercase
      IMPLICIT NONE

      CLASS(SV_DECOMPOSER),   INTENT(INOUT)  :: THIS
      CHARACTER(1),           INTENT(IN)     :: JOBU
      CHARACTER(1),           INTENT(IN)     :: JOBVT
      INTEGER,                INTENT(IN)     :: DIM_1
      INTEGER,                INTENT(IN)     :: DIM_2
      INTEGER,                INTENT(OUT)    :: INFO

      GDEF  :: ASK_DIM
      GDEF  :: DUMMY_G

      CALL THIS%DEALLOC()

      THIS%JOBU = uppercase(JOBU)
      THIS%JOBVT = uppercase(JOBVT)
      THIS%DIM_1 = DIM_1
      THIS%DIM_2 = DIM_2
      THIS%MIN_DIM = MIN(DIM_1, DIM_2)

      SELECT CASE (THIS%JOBU)
      CASE("A")
         ALLOCATE(THIS%U_MAT(DIM_1, DIM_1))
         THIS%LDU = DIM_1
      CASE("S")
         ALLOCATE(THIS%U_MAT(DIM_1, THIS%MIN_DIM))
         THIS%LDU = DIM_1
      CASE("O")
         IF (DIM_1 >= DIM_2) THEN
            ALLOCATE(THIS%U_MAT(1, 1))
            THIS%LDU = 1
         ELSE
            ALLOCATE(THIS%U_MAT(DIM_1, DIM_1))
            THIS%LDU = DIM_1
         ENDIF
      CASE("N")
         ALLOCATE(THIS%U_MAT(1, 1))
         THIS%LDU = 1
      CASE DEFAULT
         CALL vtutor%bug("SV_DECOMPOSER_INIT: JOBU '" // THIS%JOBU // "' not implemented!", __FILE__, __LINE__)
      END SELECT

      SELECT CASE (THIS%JOBVT)
      CASE("A")
         ALLOCATE(THIS%V_DAGGER_MAT(DIM_2, DIM_2))
         THIS%LDVD = DIM_2
      CASE("S")
         ALLOCATE(THIS%V_DAGGER_MAT(THIS%MIN_DIM, DIM_2))
         THIS%LDVD = THIS%MIN_DIM
      CASE("O")
         IF (DIM_1 >= DIM_2) THEN
            ALLOCATE(THIS%V_DAGGER_MAT(DIM_2, DIM_2))
            THIS%LDVD = DIM_2
         ELSE
            ALLOCATE(THIS%V_DAGGER_MAT(1, 1))
            THIS%LDVD = 1
         ENDIF
      CASE("N")
         ALLOCATE(THIS%V_DAGGER_MAT(1, 1))
         THIS%LDVD = 1
      CASE DEFAULT
         CALL vtutor%bug("SV_DECOMPOSER_INIT: JOBVT '" // THIS%JOBVT // "' not implemented!", __FILE__, __LINE__)
      END SELECT

      ALLOCATE(THIS%SINGULAR_VALS(THIS%MIN_DIM))

      ! Query optimal size of work array
#ifdef gammareal
      CALL DGESVD(THIS%JOBU, THIS%JOBVT, DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  ASK_DIM, -1, INFO)
#else
      ALLOCATE(THIS%RWORK(THIS%MIN_DIM * MAX(5 * (THIS%MIN_DIM + 1), 2 * (MAX(DIM_1, DIM_2) + THIS%MIN_DIM) + 1)))
      CALL ZGESVD(THIS%JOBU, THIS%JOBVT, DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  ASK_DIM, -1, THIS%RWORK, INFO)
#endif

      IF (INFO /= 0) THEN
         CALL THIS%DEALLOC()
         RETURN
      ENDIF

      THIS%WORK_SIZE = INT(ASK_DIM)
      ALLOCATE(THIS%WORK(THIS%WORK_SIZE))
   END SUBROUTINE


!                             ____________________
!____________________________/ SV_DECOMPOSER_CALC \_____________________________
!
!> @brief Performs singular-value decomposition of matrix A.
!>
!> The return values are either returned in A or are stored in the public type
!> variables SINGULAR_VALS, U_MAT or V_DAGGER_MAT, in accordance with the
!> chosen mode. This is consistent with how the output is defined in DGESVD and
!> ZGESDD.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in,out] A     Matrix A
!> @param[out] INFO     `INFO` parameter of DGESVD/ZGESVD
!_______________________________________________________________________________
   SUBROUTINE SV_DECOMPOSER_CALC(THIS, A, INFO)
      IMPLICIT NONE
   
      CLASS(SV_DECOMPOSER),   INTENT(INOUT)  :: THIS
      GDEF,                   INTENT(INOUT)  :: A(:, :) !! Content is destroyed on exit
      INTEGER,                INTENT(OUT)    :: INFO
      
      ! Perform singular value decomposition of A
#ifdef gammareal
      CALL DGESVD(THIS%JOBU, THIS%JOBVT, THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  THIS%WORK, THIS%WORK_SIZE, INFO)
#else
      CALL ZGESVD(THIS%JOBU, THIS%JOBVT, THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  THIS%WORK, THIS%WORK_SIZE, THIS%RWORK, INFO)
#endif
   END SUBROUTINE


!                            _______________________
!___________________________/ SV_DECOMPOSER_DEALLOC \___________________________
!
!> @brief Deallocates the buffer and work arrays.
!>
!> This can be used in case the allocated arrays are too large in between
!> operations. However, the object needs to be reinitialized for the next
!> operation. Deallocation happens at the end of scope regardless, so normally,
!> there shouldn't be a need to call this.
!>
!> @param[in,out] THIS  Instance reference
!_______________________________________________________________________________
   SUBROUTINE SV_DECOMPOSER_DEALLOC(THIS)
      IMPLICIT NONE
      CLASS(SV_DECOMPOSER),   INTENT(INOUT)  :: THIS

      IF (ALLOCATED(THIS%U_MAT)) DEALLOCATE(THIS%U_MAT)
      IF (ALLOCATED(THIS%V_DAGGER_MAT)) DEALLOCATE(THIS%V_DAGGER_MAT)
      IF (ALLOCATED(THIS%SINGULAR_VALS)) DEALLOCATE(THIS%SINGULAR_VALS)
      IF (ALLOCATED(THIS%WORK)) DEALLOCATE(THIS%WORK)
      IF (ALLOCATED(THIS%RWORK)) DEALLOCATE(THIS%RWORK)
   END SUBROUTINE


!                            ______________________
!___________________________/ RRQR_DECOMPOSER_INIT \____________________________
!
!> @brief @brief Class constructor
!>
!> Internal buffers and work arrays are set up based on the supplied matrix
!> dimensions of a matrix A that is about to be QR decomposed by this
!> calculator.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in] DIM_1     SIZE(A, 1)
!> @param[in] DIM_2     SIZE(A, 2)
!> @param[out] INFO     `INFO` parameter of xGEQP3 query
!_______________________________________________________________________________
   SUBROUTINE RRQR_DECOMPOSER_INIT(THIS, DIM_1, DIM_2, INFO)
      IMPLICIT NONE
      
      CLASS(RRQR_DECOMPOSER), INTENT(INOUT)  :: THIS
      INTEGER,                INTENT(IN)     :: DIM_1
      INTEGER,                INTENT(IN)     :: DIM_2
      INTEGER,                INTENT(OUT)    :: INFO

      GDEF  :: ASK_DIM
      GDEF  :: DUMMY_G

      CALL THIS%DEALLOC()

      THIS%DIM_1 = DIM_1
      THIS%DIM_2 = DIM_2

      ALLOCATE(THIS%PIVOT(DIM_2))
      THIS%PIVOT = 0

      ALLOCATE(THIS%TAU(MIN(DIM_1, DIM_2)))

      ! Query optimal size of work arrays
#ifdef gammareal
      CALL DGEQP3(DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%PIVOT, THIS%TAU, ASK_DIM, -1, INFO)
#else
      ALLOCATE(THIS%RWORK(2 * DIM_2))
      CALL ZGEQP3(DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%PIVOT, THIS%TAU, ASK_DIM, -1, THIS%RWORK, INFO)
#endif

      IF (INFO /= 0) THEN
         CALL THIS%DEALLOC()
         RETURN
      ENDIF

      THIS%WORK_SIZE = INT(ASK_DIM)
      ALLOCATE(THIS%WORK(THIS%WORK_SIZE))
   END SUBROUTINE


!                            ______________________
!___________________________/ RRQR_DECOMPOSER_CALC \____________________________
!
!> @brief Performs rank-revealing QR decomposition of matrix A.
!>
!> The return values are either returned in A or are stored in the public type
!> variables PIVOT and TAU. This is consistent with how the output is defined
!> in xGEQP3.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in,out] A     Matrix A
!> @param[out] INFO     `INFO` parameter of xGEQP3
!_______________________________________________________________________________
   SUBROUTINE RRQR_DECOMPOSER_CALC(THIS, A, INFO)
      IMPLICIT NONE

      CLASS(RRQR_DECOMPOSER), INTENT(INOUT)  :: THIS
      GDEF,                   INTENT(INOUT)  :: A(:, :) !! Content is destroyed on exit
      INTEGER,                INTENT(OUT)    :: INFO

      ! Perform rank-revealing QR decomposition
#ifdef gammareal
      CALL DGEQP3(THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%PIVOT, THIS%TAU, &
                  THIS%WORK, THIS%WORK_SIZE, INFO)
#else
      CALL ZGEQP3(THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%PIVOT, THIS%TAU, &
                  THIS%WORK, THIS%WORK_SIZE, THIS%RWORK, INFO)
#endif
   END SUBROUTINE


!                           _________________________
!__________________________/ RRQR_DECOMPOSER_DEALLOC \__________________________
!
!> @brief Deallocates the buffer and work arrays.
!>
!> This can be used in case the allocated arrays are too large in between
!> operations. However, the object needs to be reinitialized for the next
!> operation. Deallocation happens at the end of scope regardless, so normally,
!> there shouldn't be a need to call this.
!>
!> @param[in,out] THIS
!_______________________________________________________________________________
   SUBROUTINE RRQR_DECOMPOSER_DEALLOC(THIS)
      IMPLICIT NONE
      CLASS(RRQR_DECOMPOSER), INTENT(INOUT)  :: THIS

      IF (ALLOCATED(THIS%PIVOT)) DEALLOCATE(THIS%PIVOT)
      IF (ALLOCATED(THIS%TAU)) DEALLOCATE(THIS%TAU)
      IF (ALLOCATED(THIS%WORK)) DEALLOCATE(THIS%WORK)
      IF (ALLOCATED(THIS%RWORK)) DEALLOCATE(THIS%RWORK)
   END SUBROUTINE


!____________________________________________________________________________
!!inverts a real matrix using LU decomposition
!____________________________________________________________________________
 
   SUBROUTINE INVERT_REAL_MATRIX( MATRIX, IU0 )
      REAL(q)      :: MATRIX(:,:)
      INTEGER      :: IU0
      !local 
      INTEGER             :: N
      REAL(q),ALLOCATABLE :: A(:,:)
      INTEGER,ALLOCATABLE :: IPIV(:)
      REAL(q),ALLOCATABLE :: WORK(:)
      INTEGER             :: IWORK
      INTEGER             :: INFO

      PROFILING_START('invert_real_matrix')

      !obtain size of matrix
      N=SIZE(MATRIX,1)

      !consistency check
      IF( N/= SIZE(MATRIX,2)) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports matrix not square " // &
            str(N) // " " // str(SIZE(MATRIX,2)), __FILE__, __LINE__)
      ENDIF

      !store MATRIX to A
      ALLOCATE(A(N,N))
      A(1:N,1:N)=MATRIX(1:N,1:N)

      !allocation of IPIV
      ALLOCATE(IPIV(N))

      !LU decomposition
      CALL DGETRF(N,N,A(1,1),N,IPIV,INFO)
      IF ( INFO/=0) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports " // str(INFO), __FILE__, __LINE__)
      ENDIF

      !obtain optimal working array sizes
      ALLOCATE(WORK(1))
      CALL DGETRI(N,A(1,1),N,IPIV,WORK,-1,INFO)
      IWORK=WORK(1)
      IF ( IWORK <= 0 ) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports illegal size from " &
            // "DGETRI " // str(IWORK), __FILE__, __LINE__)
      ENDIF
      !allocate working array with optimal size
      DEALLOCATE(WORK)
      ALLOCATE(WORK(IWORK))

      !actual inversion comes now
      CALL DGETRI(N,A(1,1),N,IPIV,WORK,IWORK,INFO)
      IF ( INFO/=0) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports DGETRI reports: " // &
            str(INFO), __FILE__, __LINE__)
      ENDIF

      !store result into MATRIX, deallocate and return
      MATRIX(1:N,1:N)=A(1:N,1:N)

      PROFILING_STOP('invert_real_matrix')
   END SUBROUTINE INVERT_REAL_MATRIX

   pure function det_3x3i(a) result(det)

      integer, intent(in)  :: a(3, 3)
      integer              :: det

      det   = a(1, 1) * a(2, 2) * a(3, 3) &
            + a(2, 1) * a(3, 2) * a(1, 3) &
            + a(1, 2) * a(2, 3) * a(3, 1) &
            - a(1, 3) * a(2, 2) * a(3, 1) &
            - a(2, 1) * a(1, 2) * a(3, 3) &
            - a(3, 2) * a(2, 3) * a(1, 1)
   end function det_3x3i

   pure function matrix_inverse_3x3i(a) result(b)

      integer, intent(in)  :: a(3, 3)
      integer              :: b(3, 3)
      integer              :: det_inv

      ! calculate inverse determinant.
      det_inv = 1._q / det_3x3i(a)

      ! calculate inverse of the matrix.
      b(1, 1) = +det_inv * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
      b(2, 1) = -det_inv * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
      b(3, 1) = +det_inv * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
      b(1, 2) = -det_inv * (a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
      b(2, 2) = +det_inv * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
      b(3, 2) = -det_inv * (a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
      b(1, 3) = +det_inv * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
      b(2, 3) = -det_inv * (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
      b(3, 3) = +det_inv * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
   end function matrix_inverse_3x3i

   pure function det_3x3(a) result(det)

      real(q), intent(in)  :: a(3, 3)
      real(q)              :: det

      det   = a(1, 1) * a(2, 2) * a(3, 3) &
            + a(2, 1) * a(3, 2) * a(1, 3) &
            + a(1, 2) * a(2, 3) * a(3, 1) &
            - a(1, 3) * a(2, 2) * a(3, 1) &
            - a(2, 1) * a(1, 2) * a(3, 3) &
            - a(3, 2) * a(2, 3) * a(1, 1)
   end function det_3x3

   pure function matrix_inverse_3x3(a) result(b)

      real(q), intent(in)  :: a(3, 3)
      real(q)              :: b(3, 3)
      real(q)              :: det_inv

      ! calculate inverse determinant.
      det_inv = 1._q / det_3x3(a)

      ! calculate inverse of the matrix.
      b(1, 1) = +det_inv * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
      b(2, 1) = -det_inv * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
      b(3, 1) = +det_inv * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
      b(1, 2) = -det_inv * (a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
      b(2, 2) = +det_inv * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
      b(3, 2) = -det_inv * (a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
      b(1, 3) = +det_inv * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
      b(2, 3) = -det_inv * (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
      b(3, 3) = +det_inv * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
   end function matrix_inverse_3x3


!> A modified version of algorithm 680, rewritten in Fortran 2008.
!> G.P.M. Poppe, C.M.J. Wijers, More efficient computation of
!> the complex error-function, ACM Trans. Math. Software 16:38-46, 1990.
!>  and
!> G.P.M. Poppe, C.M.J. Wijers, Algorithm 680, Evaluation of the
!> complex error function, ACM Trans. Math. Software 16:47, 1990.
!
!> Given a complex number z, this function computes
!> the value of the Faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
!> where erfc is the complex complementary error-function and i
!> means sqrt(-1).  The accuracy of the algorithm for z in the 1st
!> and 2nd quadrant is 14 significant digits; in the 3rd and 4th
!> it is 13 significant digits outside a circular region with radius
!> 0.126 around a zero of the function.
!>
!> This code was taken from: https://sourceforge.net/projects/erf/files/
!> This is distributed with the BSD licence
!>
!
!$Id: LICENSE 4 2016-08-11 15:34:59Z mexas $
!
!Copyright (c) 2016, Anton Shterenlikht, The University of Bristol, UK
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are
!met:
!
!1. Redistributions of source code must retain the above copyright
!notice, this list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright
!notice, this list of conditions and the following disclaimer in the
!documentation and/or other materials provided with the distribution.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
!TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    function wpop( z )
GPACC !$ACC ROUTINE SEQ
GPOMP !$OMP DECLARE TARGET(wpop)
      complex( kind=q ) :: wpop
      complex( kind=q ) :: z
      real( kind=q ), parameter :: pi = 3.141592653589793238_q
      real( kind=q ), parameter :: factor = 2/sqrt(pi)
      logical :: a, b

      real( kind=q ) :: xabs, yabs, x, y, qrho
      real( kind=q ) :: xabsq, yquad, c, daux
      real( kind=q ) :: h, h2, qlambda, rx, ry, sx, sy, tx, ty, u1, u2, v1, v2, w1, xaux
      real( kind=q ) :: xquad, xsum, ysum, xi, yi, u, v
      integer :: i, j, kapn, n, np1, nu

      ! To avoid the complier uninitialised varning
      h2 = 0.0_q

      xi = real( z, kind=q )
      yi = aimag( z )
      xabs = abs( xi )
      yabs = abs( yi )
      x = xabs / 6.3_q
      y = yabs / 4.4_q
      qrho = x**2 + y**2
      xabsq = xabs**2
      xquad = xabsq - yabs**2
      yquad = 2*xabs*yabs

      a = qrho .lt. 0.085264_q

      branch1: if ( a ) then

        ! If ( qrho .lt. 0.085264 ) then the Faddeeva-function is evaluated
        ! using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
        ! n is the minimum number of terms needed to obtain the required
        ! accuracy
        qrho  = ( 1.0_q - 0.85_q * y ) * sqrt( qrho )
        n = nint( 6.0_q + 72.0_q * qrho )
        j = 2 * n + 1
        xsum = 1.0_q / real( j, kind=q )
        ysum = 0.0_q

        do i = n, 1, -1
           j = j - 2
           xaux = ( xsum * xquad - ysum * yquad) / real( i, kind=q )
           ysum = ( xsum * yquad + ysum * xquad) / real( i, kind=q )
           xsum = xaux + 1.0_q / real( j, kind=q )
        end do
        u1 = - factor * ( xsum * yabs + ysum * xabs ) + 1.0_q
        v1 = factor * ( xsum * xabs - ysum * yabs )
        daux = exp( -xquad )
        u2 = daux * cos( yquad )
        v2 = -daux * sin( yquad )
        u = u1 * u2 - v1 * v2
        v = u1 * v2 + v1 * u2
      else
        bran2: if ( qrho .gt. 1.0_q ) then
           ! If ( qrho .gt. 1) then w(z) is evaluated using the laplace
           ! continued fraction. nu is the minimum number of terms needed
           ! to obtain the required accuracy.
           h = 0.0_q
           kapn = 0
           qrho = sqrt( qrho )
           nu = int( 3.0_q + (1442.0_q / ( 26.0_q * qrho + 77.0_q )))
        else
           ! If ( qrho .ge. 0.085264 .and. qrho .le. 1.0_q ) then
           ! w(z) is evaluated by a truncated Taylor expansion,
           ! where the Laplace continued fraction is used to calculate
           ! the derivatives of w(z). KAPN is the minimum number of terms
           ! in the Taylor expansion needed to obtain the required accuracy.
           ! NU is the minimum number of terms of the continued fraction
           ! needed to calculate the derivatives with the required accuracy.
           qrho = ( 1.0_q - y ) * sqrt( 1.0_q - qrho )
           h = 1.88_q * qrho
           h2 = 2.0_q * h
           kapn = nint( 7.0_q + 34.0_q * qrho )
           nu   = nint( 16.0_q + 26.0_q * qrho )
        end if bran2

        b = h .gt. 0.0_q

        ! To avoid uninitialise compiler warning. qlambda is used
        ! only if (b), so can define to any value otherwise.
        qlambda = 0.0_q
        if ( b ) qlambda = h2**kapn

        rx = 0.0_q
        ry = 0.0_q
        sx = 0.0_q
        sy = 0.0_q

        do n = nu, 0, -1
           np1 = n + 1
           tx = yabs + h + np1 * rx
           ty = xabs - np1 * ry
           c = 0.5_q / (tx**2 + ty**2)
           rx = c * tx
           ry = c * ty
           if ( b .and. n .le. kapn ) then
              tx = qlambda + sx
              sx = rx*tx - ry*sy
              sy = ry*tx + rx*sy
              qlambda = qlambda / h2
           end if
        end do
           if ( h .eq. 0.0_q ) then
              u = factor * rx
              v = factor * ry
           else
              u = factor * sx
              v = factor * sy
           end if
           if ( yabs .eq. 0.0_q ) u = exp( -xabs**2 )
      end if branch1

      ! Evaluation of w(z) in the other quadrants
      if ( yi .lt. 0.0_q ) then
         if ( a ) then
            u2 = 2.0_q * u2
            v2 = 2.0_q * v2
         else
            xquad = -xquad
            w1 = 2.0_q * exp( xquad )
            u2 = w1 * cos( yquad )
            v2 = -w1 * sin( yquad )
         end if
         u = u2 - u
         v = v2 - v
         if ( xi .gt. 0.0_q ) v = -v
      else
         if ( xi .lt. 0.0_q ) v = -v
      end if
      wpop = cmplx( u, v, kind=q )

    end function wpop

    !> This is an error function of a complex argument, which uses wpop(z).
    function erf_pop( z )
GPACC !$ACC ROUTINE SEQ
GPOMP !$OMP DECLARE TARGET(erf_pop)
       complex(q) :: erf_pop
       complex(q) :: z
       ! local variables
       complex(q), parameter :: ci = cmplx( 0.0_q, 1.0_q, kind=q )
GPACC !$ACC ROUTINE(WPOP) SEQ
       erf_pop = 1.0_q - wpop( ci * z ) * exp( - z**2 )
    end function erf_pop

    !> This is an error function times a gaussian function of a complex argument, which uses wpop(z)
    !> erf_gaussian_pop(a+ib) = e^{-b^2}erf(a+ib)
    function erf_gaussian_pop( z )
GPACC !$ACC ROUTINE SEQ
GPOMP !$OMP DECLARE TARGET(erf_gaussian_pop)
       complex(q) :: erf_gaussian_pop
       complex(q) :: z
       ! local variables
       complex(q) :: c,d
       real(q) :: a,b
GPACC !$ACC ROUTINE(WPOP) SEQ
       a = real( z, kind=q )
       b = aimag( z )
       c = cmplx( -b, a, kind=q )
       d = cmplx( 0._q, -2._q*a*b, kind=q )
       erf_gaussian_pop = exp( - b*b ) - wpop( c ) * exp( - a*a ) * exp( d )
    end function erf_gaussian_pop

END MODULE mathtools
