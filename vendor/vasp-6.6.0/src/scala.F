#include "symbol.inc"
!=======================================================================
! RCS:  $Id: scala.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!=======================================================================
!> Module containing wrapper for scaLAPACK
!>
!> written by Gilles de Wijs (gD) and Georg Kresse (gK),
!> modified to run on any number of nodes by Dario Alfe,
!> BG_ routines were added by IBM fellow (bG)
!>
!> @note This is a common note for routines which are called with descriptor
!> DESCA (subroutines *_DESC and *_DESC2).
!> The routines can be called from any context or subgroup
!> of MPI processed as long as DESCA is properly set up.
!> They also do not refer to GSD. In that sense they are much safer to use, but
!> matter of fact the caller must "know" scalapack.
!> For the time being the check however that the passed matrix dimensions are
!> consistent with those set up in GSD by calling CHECK_scala. This should be
!> removed after more extensive testing.
!>
!=======================================================================

#ifdef scaLAPACK
  MODULE scala
    USE scala_helper
    USE prec
    USE mpimy
#ifdef ELPA
    USE elpa, ONLY : elpa_t, ELPA_INIT, ELPA_ALLOCATE, ELPA_SOLVER_2STAGE, ELPA_OK, ELPA_DEALLOCATE, ELPA_UNINIT
#endif

    IMPLICIT NONE

    ! customization of matrix diagonalization
    !> maximum cluster of eigenvector
    !>
    !>  see manpage for pSSYEVX pZHEEVX:
    !>    NCLUST is the number of eigenvalues in the largest cluster,
    !>    where a cluster is defined as a set of close eigenvalues:
    !>    {W(K),...,W(K+NCLUST-1)|W(J+1)<= W(J)+orfac*norm(A)}
    INTEGER, PRIVATE :: NCLUST=48
    REAL(q), PRIVATE :: ABSTOL_DEF=1E-10_q  !< specifies eigenvector orthogonality tolerance
    REAL(q), PRIVATE :: ORFAC=-1._q         !< controls reorthogonalisation of eigenvectors
    REAL(qs),PRIVATE :: ORFAC_single=-1._qs !< controls reorthogonalisation of eigenvectors

! Many routines in this module use a global scalapack descriptor GSD
! which is set up by INIT_SCALA (scala_struct.F). INIT_SCALA can be called
! repeatedly if the matrix dimensions or the processor grid change
! however, all subroutines that have an argument DESCA
! are save and they can be used in different contexts

! Set SET_EVX_WORKSPACE to bypass the PDSYEVX and PZHEEVX workspace queries in BG_PDSYEVX_ZHEEVX
! and PDSYEVX_ZHEEVX_DESC, and call the workspace query routine PDSYEVX_PZHEEVX_WORKSPACE instead.
! The aforementioned workspace queries are buggy in MKL.
#define SET_EVX_WORKSPACE

  CONTAINS
!=======================================================================
!
! the following routines rely on the global GSD communicator
! they however re-initialize GSD as well, so should be safe
! against changes of the matrix dimension
!
!=======================================================================
!
!> this subroutine determines U^-1: A= U+ U and invert U
!>
!> 1. calls INIT_scala
!> 2. calls the setup of distributed the distributed matrix
!> 3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!> 4. calls reconstruction of distributed data into patched matrix
!>
!> @note The sum over processors of the patched matrix is not realised
!>       in this module
!
!=======================================================================


    SUBROUTINE PPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           !< leading dimension of AMATIN
      INTEGER N              !< NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N) !< input/output matrix
! local variables
      INTEGER INFO

      GDEF, ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      PROFILING_START('ppotrf_trtri')

      CALL INIT_scala(COMM, N)

      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,GSD%DESCA)

        INFO=0
#ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, GSD%DESCA, INFO )
#else
        CALL PZPOTRF( 'U', N, A, 1, 1, GSD%DESCA, INFO )
#endif

        IF (INFO.NE.0) THEN
          CALL vtutor%error("PPOTRF_TRTRI, POTRF, INFO: " // str(INFO) // "\n STOP")
        ENDIF
#ifdef gammareal
         CALL PDTRTRI &
#else
         CALL PZTRTRI &
#endif
     &    ('U', 'N', N, A, 1, 1, GSD%DESCA, INFO)
        IF (INFO.NE.0) THEN
          CALL vtutor%error("PPOTRF_TRTRI, TRTRI, INFO: " // str(INFO))
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,GSD%DESCA)
   ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      PROFILING_STOP('ppotrf_trtri')

      RETURN
    END SUBROUTINE PPOTRF_TRTRI


!=======================================================================
!
!> this subroutine determines A^-1: A= U+ U and invert (U+ U)
!>
!> 1. calls INIT_scala
!> 2. calls the setup of distributed the distributed matrix
!> 3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!> 4. calls reconstruction of distributed data into patched matrix
!>
!> @note The sum over processors of the patched matrix is not realised
!>       in this module
!
!=======================================================================


    SUBROUTINE PPOTRF_POTRI (COMM, AMATIN, NDIM, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           ! leading dimension of AMATIN
      INTEGER N              ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N) ! input/output matrix
! local variables
      INTEGER INFO

      GDEF, ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      PROFILING_START('ppotrf_potri')

      CALL INIT_scala(COMM, N)
!-----------------------------------------------------------------------
! allocate workarray
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,GSD%DESCA)

        INFO=0
#ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, GSD%DESCA, INFO )
#else
        CALL PZPOTRF( 'U', N, A, 1, 1, GSD%DESCA, INFO )
#endif

        IF (INFO.NE.0) THEN
          CALL vtutor%error("PPOTRF_TRTRI, POTRF, INFO: " // str(INFO) // "\n STOP")
        ENDIF
#ifdef gammareal
         CALL PDPOTRI &
#else
         CALL PZPOTRI &
#endif
     &    ('U', N, A, 1, 1, GSD%DESCA, INFO)
        IF (INFO.NE.0) THEN
          CALL vtutor%error("PPOTRF_POTRI, POTRI, INFO: " // str(INFO))
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,GSD%DESCA)
   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      PROFILING_STOP('ppotrf_potri')

      RETURN
    END SUBROUTINE PPOTRF_POTRI


!=======================================================================
!
!> Call to pZHEEVX respectively pDSYEVX
!>
!> i.e. diagonalization of a symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE PDSYEVX_ZHEEVX(COMM, AMATIN, W, NDIM, N, ABSTOL_, AMATIN_SCALA)
#ifdef ELPA
      USE iso_c_binding
#endif
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          !< NxN matrix to be distributed
      INTEGER NDIM                       !< leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             !< input/output matrix
      REAL(q) W(N)                       !< eigenvalues
      REAL(q), OPTIONAL      ::  ABSTOL_ !< tolerance
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) !< alternative input/output matrix (scaLAPACK format)
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2
      REAL(q) :: ABSTOL

      GDEF, ALLOCATABLE, TARGET :: A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, ALLOCATABLE, TARGET :: Z(:) ! Z eigenvector matrix
      GDEF, POINTER        ::     AP(:) ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array

      INTEGER :: LWWORK, LIWORK, LRWORK ! work array sizes

      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE

#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS

      GDEF, POINTER :: EP(:,:), EZ(:,:)
#endif

      PROFILING_START('pdsyevx_zheevx')

      CALL INIT_scala(COMM, N)

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("PDSYEVX_ZHEEVX: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         IF (SIZE(AMATIN_SCALA)<GSD%MALLOCPQ) CALL vtutor%bug("PDSYEVX_ZHEEVX: AMATIN_SCALA too small: " &
               // str(SIZE(AMATIN_SCALA)) // " " // str(GSD%MALLOCPQ),__FILE__,__LINE__)
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF

      ALLOCATE(Z(GSD%MALLOCPQ))
DOESI Z = 0

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A, GSD%DESCA)
         ENDIF
         INFO=0

#ifdef ELPA
         ! call ELPA 2-stage solver
         IF (LELPA) THEN
            IF (ELPA_INIT(20171201) /= ELPA_OK) &
               CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

            CALL c_f_pointer(c_loc(AP),EP,[GSD%NP,GSD%NQ])
            CALL c_f_pointer(c_loc( Z),EZ,[GSD%NP,GSD%NQ])

            ET => ELPA_ALLOCATE()

            CALL ET%SET("na",  N, ISUCCESS)
            CALL ET%SET("nev", N, ISUCCESS)

            CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
            CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

            CALL ET%SET("nblk", GSD%DESCA(MB_), ISUCCESS)

            CALL ET%SET("mpi_comm_parent", COMM%MPI_COMM, ISUCCESS)
            CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
            CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

            ISUCCESS = ET%SETUP()

            !CALL ET%SET("nvidia-gpu", 1, ISUCCESS)
            CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)

            CALL ET%EIGENVECTORS(EP, W, EZ, ISUCCESS)

            CALL ELPA_DEALLOCATE(ET)
            CALL ELPA_UNINIT()

            NULLIFY(EP,EZ)
         ENDIF
#endif
         ! call scaLAPACK routine
         IF (.NOT. LELPA) THEN

            ! first call PDSYEVX/PZHEEVX to determine work array sizes
            LWWORK=-1 ; LIWORK=-1 ; LRWORK=-1
            ALLOCATE(RWORK(3), WWORK(3), IWORK(1))
DOESI       RWORK = 0; WWORK = 0; IWORK = 0

#ifdef gammareal
            CALL PDSYEVX( 'V', 'A', 'U', N, AP(1), 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                       -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, WWORK, &
                       LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )

            RWORK(1)=1
            WWORK(1)=WWORK(1)+NCLUST*N
#else
            CALL PZHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                       -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, &
                       WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
                       IFAIL, ICLUSTR, GAP, INFO )

            RWORK(1)=RWORK(1)+NCLUST*N
#endif

            LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
            LIWORK= MAX(1,IWORK(1))
            LRWORK= MAX(1,NINT(RWORK(1)))

            ! reallocate work arrays with optimal sizes
            DEALLOCATE(RWORK, WWORK, IWORK)
            ALLOCATE(RWORK(MAX(3,LRWORK)), WWORK(MAX(3,LWWORK)), IWORK(LIWORK))
DOESI       RWORK = 0.0_q; WWORK = 0; IWORK = 0

            ! and now the actual PDSYEVX/PZHEEVX calls to get the eigenvalues and eigenvectors of the matrix AP
#ifdef gammareal
            CALL PDSYEVX( 'V', 'A', 'U', N, AP(1), 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                       -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, WWORK, &
                       LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
            CALL PZHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                       -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, &
                       WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
                       IFAIL, ICLUSTR, GAP, INFO )
#endif

            DEALLOCATE(IWORK, WWORK, RWORK)

            IF (M.NE.N) THEN
               CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // &
                  str(LIWORK) // " " // str(GSD%DESCA(M_)) // "\nERROR in subspace rotation PDSYEVX/ &
                  &PZHEEVX: not enough eigenvalues found " // str(M) // " " // str(N))
            ENDIF
            IF (NZ.NE.N) THEN
               CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough eigenvectors &
                  &computed " // str(M) // " " // str(N))
            ENDIF
            DO I2=1,N
               IF (IFAIL(I2).NE.0) THEN
                  CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: I2,IFAIL= " // str(I2) &
                     // " " // str(IFAIL(I2)))
               ENDIF
            ENDDO
            IF (INFO<0) THEN
               CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ELSE IF (INFO>0) THEN
               IF (MOD(INFO,2).NE.0) THEN
                  CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
               ELSEIF (MOD(INFO/2,2).NE.0) THEN
                  ! this error condition happens very often, but does not imply that
                  ! the result is not usefull
                  ! IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
               ELSE
                  CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
               ENDIF
            ENDIF
         ENDIF

         IF (PRESENT(AMATIN_SCALA)) THEN
            AMATIN_SCALA(1:GSD%MALLOCPQ)=Z(1:GSD%MALLOCPQ)
         ELSE
            CALL RECON(AMATIN, NDIM, N, Z, GSD%DESCA)
         ENDIF

      ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF
      DEALLOCATE(Z)

      PROFILING_STOP('pdsyevx_zheevx')

      RETURN
    END SUBROUTINE PDSYEVX_ZHEEVX

!=======================================================================
!
!> Call to pZHEEVD respectively pDSYEVD
!>
!> i.e. diagonalization of a symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pDSYEV_ZHEEVD(COMM, AMATIN, W, NDIM, N, AMATIN_SCALA)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          !< NxN matrix to be distributed
      INTEGER NDIM                       !< leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             !< input/output matrix
      REAL(q) W(N)                       !< eigenvalues
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) !< alternative input/output matrix (scaLAPACK format)
! local variables
      GDEF, ALLOCATABLE,TARGET ::  A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, POINTER  ::  AP(:)          ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      INTEGER NODE_ME,IONODE

      CALL INIT_scala(COMM, N)
      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays (AP=>A)
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A(1), GSD%DESCA)
         ENDIF

         CALL PDSYEV_ZHEEVD_DESC(AP(:), W, N, GSD%DESCA, COMM)
! if desired this can be also replaced by the standard matrix diagonalization
!         CALL PDSYEV_ZHEEV_DESC(AP(:), W, N, GSD%DESCA, COMM)

         IF (PRESENT(AMATIN_SCALA)) THEN
            AMATIN_SCALA=AP
         ELSE
            CALL RECON(AMATIN, NDIM, N, A(1), GSD%DESCA)
         ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF

      RETURN
    END SUBROUTINE

!=======================================================================
!
!> Call to pZHEEV respectively pDSYEV
!>
!> i.e. diagonalization of a symmetric or hermitian matrix
!>  using standard alorithm
!
!=======================================================================

    SUBROUTINE pDSYEV_ZHEEV(COMM, AMATIN, W, NDIM, N, AMATIN_SCALA)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          !< NxN matrix to be distributed
      INTEGER NDIM                       !< leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             !< input/output matrix
      REAL(q) W(N)                       !< eigenvalues
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) !< alternative input/output matrix (scaLAPACK format)
! local variables
      GDEF, ALLOCATABLE,TARGET ::  A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, POINTER  ::  AP(:)          ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      INTEGER NODE_ME,IONODE

      CALL INIT_scala(COMM, N)
      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays (AP=>A)
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A(1), GSD%DESCA)
         ENDIF

         CALL PDSYEV_ZHEEV_DESC(AP(:), W, N, GSD%DESCA, COMM)

         IF (PRESENT(AMATIN_SCALA)) THEN
            AMATIN_SCALA=AP
         ELSE
            CALL RECON(AMATIN, NDIM, N, A(1), GSD%DESCA)
         ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF

      RETURN
    END SUBROUTINE

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version is virtually identical to PDSYEVX_ZHEEVX
!> when called with the additional AMATIN_SCALA argument,
!> i.e., `BG_PDSYEVX_ZHEEVX(COMM, A, W, N)` could be replaced by
!> `PDSYEVX_ZHEEVX(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)`.
!> A is a two dimension array here, however.
!
!=======================================================================

    SUBROUTINE BG_PDSYEVX_ZHEEVX(COMM, A, W, N)

__GPU USE moffload_struct_def
__GPU USE moffload_control

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               !< NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))  !< input/output matrix
      REAL(q) W(N)            !< eigenvalues
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2

      GDEF    Z(MAX(GSD%NP,1),MAX(GSD%NQ,1))
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
      INTEGER :: LWWORK, LIWORK, LRWORK ! work array sizes
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE
#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS
#endif

      PROFILING_START('bg_pdsyevx_zheevx')

      CALL CHECK_scala (N, 'BG_PDSYEVX_ZHEEVX' )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        CALL vtutor%error("PDSYEVX_ZHEEVX: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " " &
           // str(COMM%NCPU))
      ENDIF

!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

#if defined(CUSOLVERMP) && ! PGI_BEFORE_25_9
         IF (.NOT. LUSENCCL) THEN
GPACC !$ACC UPDATE SELF(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
            CALL OFFLOAD_ON_PUSH(.FALSE.)
         ENDIF
#endif

#ifdef CUSOLVERMP
GPACC !$ACC ENTER DATA CREATE(Z,W) __IF_ASYNC__
#else
GPACC !$ACC UPDATE SELF(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(A) IF(OFFLOAD_ON)
__GPU    CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

         INFO=0

#ifdef ELPA
         ! call ELPA 2-stage solver
         IF (LELPA) THEN
            IF (ELPA_INIT(20171201) /= ELPA_OK) &
               CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

            ET => ELPA_ALLOCATE()

            CALL ET%SET("na",  N, ISUCCESS)
            CALL ET%SET("nev", N, ISUCCESS)

            CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
            CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

            CALL ET%SET("nblk", GSD%DESCA(MB_), ISUCCESS)

            CALL ET%SET("mpi_comm_parent", COMM%MPI_COMM, ISUCCESS)
            CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
            CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

            ISUCCESS = ET%SETUP()

GPACC       IF (OFFLOAD_ON) THEN
GPACC          CALL ET%SET("nvidia-gpu", 1, ISUCCESS)
GPACC       ELSE
               CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)
GPACC       ENDIF

            CALL ET%EIGENVECTORS(A, W, Z, ISUCCESS)

            CALL ELPA_DEALLOCATE(ET)
            CALL ELPA_UNINIT()
         ENDIF
#endif
         ! call scaLAPACK routine
         IF (.NOT. LELPA) THEN

            ! first call PDSYEVX/PZHEEVX to determine work array sizes
            LWWORK=-1 ; LIWORK=-1 ; LRWORK=-1
            ALLOCATE(RWORK(3), WWORK(3), IWORK(1))
DOESI       RWORK = 0; WWORK = 0; IWORK = 0

#ifdef SET_EVX_WORKSPACE
            CALL PDSYEVX_PZHEEVX_WORKSPACE( GSD%DESCA, NCLUST, N, IWORK(1), RWORK, WWORK )
#else
#ifdef gammareal
            CALL WPDSYEVX( 'V', 'A', 'U', N, A, 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                        -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, WWORK, &
                        LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )

            RWORK(1)=1
            WWORK(1)=WWORK(1)+NCLUST*N
#else
            CALL WPZHEEVX( 'V', 'A', 'U', N, A, 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                        -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, &
                        WWORK, LWWORK, RWORK, LRWORK, IWORK, LIWORK, &
                        IFAIL, ICLUSTR, GAP, INFO )

            RWORK(1)=RWORK(1)+NCLUST*N
#endif
#endif
            LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
            LIWORK= MAX(1,IWORK(1))
            LRWORK= MAX(1,NINT(RWORK(1)))

            ! reallocate work arrays with optimal sizes
            DEALLOCATE(RWORK, WWORK, IWORK)
            ALLOCATE(RWORK(MAX(3,LRWORK)), WWORK(MAX(3,LWWORK)), IWORK(LIWORK))
DOESI       RWORK = 0.0_q; WWORK = 0; IWORK = 0

            ! and now the actual PDSYEVX/PZHEEVX calls to get the eigenvalues and eigenvectors of the matrix AP
#ifdef gammareal
            CALL WPDSYEVX( 'V', 'A', 'U', N, A, 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                        -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, WWORK, &
                        LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
            CALL WPZHEEVX( 'V', 'A', 'U', N, A, 1, 1, GSD%DESCA, ZERO, ZERO, 13,  &
                        -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, GSD%DESCA, &
                        WWORK, LWWORK, RWORK, LRWORK, IWORK, LIWORK, &
                        IFAIL, ICLUSTR, GAP, INFO )
#endif

             DEALLOCATE(IWORK, WWORK, RWORK)

#ifndef CUSOLVERMP
            IF (M.NE.N) THEN
              CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // &
                 str(LIWORK) // " " // str(GSD%DESCA(M_)) // "\nERROR in subspace rotation PDSYEVX: not &
                 &enough eigenvalues found " // str(M) // " " // str(N))
            ENDIF
            IF (NZ.NE.N) THEN
              CALL vtutor%error("ERROR in subspace rotation PDSYEVX: not enough eigenvectors computed " &
                 // str(M) // " " // str(N))
            ENDIF
            DO I2=1,N
              IF (IFAIL(I2).NE.0) THEN
                CALL vtutor%error("ERROR in subspace rotation PDSYEVX: I2,IFAIL= " // str(I2) // " " // &
                   str(IFAIL(I2)))
              ENDIF
            ENDDO
            IF (INFO<0) THEN
               CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ELSE IF (INFO>0) THEN
               IF (MOD(INFO,2).NE.0) THEN
                  CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
               ELSEIF (MOD(INFO/2,2).NE.0) THEN
                  ! this error condition happens very often, but does not imply that
                  ! the result is not usefull
                  ! IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
               ELSE
                  CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
               ENDIF
            ENDIF
#endif
         ENDIF

#ifdef CUSOLVERMP
GPACC !$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
#endif
         A=Z
#ifdef CUSOLVERMP
GPACC !$ACC END KERNELS
GPACC !$ACC EXIT DATA DELETE(Z) COPYOUT(W) __IF_ASYNC__
#else
__GPU    CALL OFFLOAD_ON_POP()
GPACC !$ACC UPDATE DEVICE(A) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(A) IF(OFFLOAD_ON)
#endif

#if defined(CUSOLVERMP) && ! PGI_BEFORE_25_9
         IF (.NOT. LUSENCCL) THEN
            CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         ENDIF
#endif

      ENDIF calc


      PROFILING_STOP('bg_pdsyevx_zheevx')

      RETURN
    END SUBROUTINE BG_PDSYEVX_ZHEEVX

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version is virtually identical to pDSSYED_ZHEEVD
!> when called with the additional AMATIN_SCALA argument,
!> i.e., `BG_pDSSYED_ZHEEVD(COMM, A, W, N)` could be replaced by
!> `pDSSYED_ZHEEVD(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)`
!> A is a two dimension array here, however.
!
!=======================================================================

    SUBROUTINE BG_pDSYEV_ZHEEVD(COMM, A, W, N)

GPACC USE moffload_struct_def

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               !< NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))  !< input/output matrix
      REAL(q) W(N)            !< eigenvalues
! local variables
      INTEGER NODE_ME,IONODE

      PROFILING_START('bg_pdsyev_zheevd')

!$ACC UPDATE SELF(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

      CALL CHECK_scala (N, 'BG_pDSYEV_ZHEEVD' )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " " &
           // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

        CALL PDSYEV_ZHEEVD_DESC2(A, W, N, GSD%DESCA, COMM)
      ENDIF calc

#ifdef CUSOLVERMP
      !$ACC UPDATE SELF(A) __IF_ASYNC__
#else
      !$ACC UPDATE DEVICE(A) __IF_ASYNC__
#endif

      PROFILING_STOP('bg_pdsyev_zheevd')
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      RETURN
    END SUBROUTINE

!=======================================================================
!
!> Matrix LU decompositions scaLAPACK aware version
!>
!> Determine A = U+ U (POTRF) and invert U (TRTRI)
!
!=======================================================================

    SUBROUTINE BG_PPOTRF_TRTRI(A, N, INFO)

GPACC USE moffload_struct_def
GPACC USE moffload_control

      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      INTEGER N                              !< NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1)) !< input/output matrix
      INTEGER INFO

#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS
#endif

      PROFILING_START('bg_ppotrf_trtri')

      INFO=0

      CALL CHECK_scala (N, 'BG_PPOTRF_TRTRI' )

#ifdef ELPA
      ! call ELPA cholesky and invert
      IF (LELPA) THEN
         IF (ELPA_INIT(20171201) /= ELPA_OK) &
            CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

         ET => ELPA_ALLOCATE()

         CALL ET%SET("na",  N, ISUCCESS)
         CALL ET%SET("nev", N, ISUCCESS)

         CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
         CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

         CALL ET%SET("nblk", GSD%DESCA(MB_), ISUCCESS)

         CALL ET%SET("mpi_comm_parent", GSD%ELPA_MPI_COMM_ALL, ISUCCESS)
         CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
         CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

         ISUCCESS = ET%SETUP()

GPACC    IF (OFFLOAD_ON) CALL ET%SET("nvidia-gpu", 1, ISUCCESS)

         CALL ET%CHOLESKY(A, ISUCCESS)
         CALL ET%INVERT_TRIANGULAR(A, ISUCCESS)

         CALL ELPA_DEALLOCATE(ET)
         CALL ELPA_UNINIT()
      ENDIF
#endif

      IF (.NOT. LELPA) THEN
#if defined(CUSOLVERMP) && ! PGI_BEFORE_25_9
         IF (.NOT. LUSENCCL) THEN
GPACC !$ACC UPDATE SELF(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
           CALL OFFLOAD_ON_PUSH(.FALSE.)
         ENDIF
#endif

#if ! defined(CUSOLVERMP)
GPACC !$ACC UPDATE SELF(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC    CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif

#ifdef gammareal
         CALL WPDPOTRF &
#else
         CALL WPZPOTRF &
#endif
           & ('U',N, A(1,1),1,1,GSD%DESCA,INFO)

#if defined(CUSOLVERMP) || defined(CUBLASMP)
#if ! defined(CUSOLVERMP)
!$ACC UPDATE DEVICE(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q)
         CALL OFFLOAD_ON_POP
#endif
#if ! defined(CUBLASMP)
!$ACC UPDATE SELF(A) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q)
         CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif
#endif

         IF (INFO/=0) RETURN

#ifdef gammareal
         CALL WPDTRTRI &
#else
         CALL WPZTRTRI &
#endif
             & ('U','N',N,A(1,1),1,1,GSD%DESCA,INFO)

#if ! defined(CUBLASMP)
GPACC    CALL OFFLOAD_ON_POP
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC UPDATE DEVICE(A) __IF_ASYNC__
#endif

#if defined(CUSOLVERMP) && ! PGI_BEFORE_25_9
         IF (.NOT. LUSENCCL) THEN
            CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         ENDIF
#endif
      ENDIF

      PROFILING_STOP('bg_ppotrf_trtri')

      END SUBROUTINE BG_PPOTRF_TRTRI

!=======================================================================
!
!> Matrix LU decompositions scaLAPACK aware version
!>
!> Determine A = U+ U (POTRF) and invert A = U+U (POTRI)
!
!=======================================================================

    SUBROUTINE BG_PPOTRF_POTRI(A, N, INFO)
      IMPLICIT NONE

      INTEGER N                              ! NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1)) ! input/output matrix
      INTEGER INFO

      PROFILING_START('bg_ppotrf_potri')

      CALL CHECK_scala (N, 'BG_PPOTRF_POTRI' )
#ifdef gammareal
         CALL PDPOTRF &
#else
         CALL PZPOTRF &
#endif
          & ('U',N, A(1,1),1,1,GSD%DESCA,INFO)
#ifdef gammareal
         CALL PDPOTRI &
#else
         CALL PZPOTRI &
#endif
          & ('U',N, A(1,1), 1, 1, GSD%DESCA, INFO)

      PROFILING_STOP('bg_ppotrf_potri')

     END SUBROUTINE BG_PPOTRF_POTRI

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version expects a proper scaLAPACK descriptor.
!> @note This routine does stop when eigenvectors cannot be
!> determined
!
!=======================================================================

     SUBROUTINE PDSYEVX_ZHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN, ABSTOL_)
#if defined ELPA
       USE iso_c_binding
#endif
GPACC  USE moffload_struct_def
GPACC  USE moffload_control

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF, TARGET :: A(:)    !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
       REAL(q), OPTIONAL :: ABSTOL_
! local variables
       REAL(q) :: ABSTOL
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(q) :: GAP(COMM%NCPU)
       REAL(q) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2

       GDEF, ALLOCATABLE, TARGET ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

#ifdef ELPA
       CLASS(elpa_t), POINTER :: ET
       INTEGER :: ISUCCESS

       GDEF, POINTER :: EP(:,:), EZ(:,:)
       INTEGER NP,NQ
       INTEGER MYROW, MYCOL, NPROW, NPCOL
       INTEGER    NUMROC
       EXTERNAL   NUMROC,BLACS_GRIDINFO
#endif

       PROFILING_START('pdsyevx_zheevx_desc')

       JOB='V'  ! default is eigenvalues and eigenvectors

       IF (PRESENT(JOBIN)) JOB=JOBIN

#if defined(CUSOLVERMP) && ! PGI_BEFORE_25_9
       IF (.NOT. LUSENCCL) THEN
GPACC !$ACC UPDATE SELF(A) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          CALL OFFLOAD_ON_PUSH(.FALSE.)
       ENDIF
#endif

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(3), WWORK(3), IWORK(1))
DOESI  RWORK = 0.0_q; WWORK = 0; IWORK = 0
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
DOESI  Z = 0
       !$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF

#ifdef ELPA
       CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
       NP = NUMROC(DESCA_(M_),DESCA(MB_),MYROW,0,NPROW)
       NQ = NUMROC(DESCA_(N_),DESCA(NB_),MYCOL,0,NPCOL)

       ! call ELPA 2-stage solver
       IF (LELPA) THEN
          IF (ELPA_INIT(20171201) /= ELPA_OK) &
            CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

          CALL c_f_pointer(c_loc(A(1)),EP,[NP,NQ])
          CALL c_f_pointer(c_loc(Z),EZ,[NP,NQ])

          ET => ELPA_ALLOCATE(ISUCCESS)

          CALL ET%SET("na",  N, ISUCCESS)
          CALL ET%SET("nev", N, ISUCCESS)

          CALL ET%SET("local_nrows", NP, ISUCCESS)
          CALL ET%SET("local_ncols", NQ, ISUCCESS)

          CALL ET%SET("nblk", DESCA(MB_), ISUCCESS)

          CALL ET%SET("mpi_comm_parent", COMM%MPI_COMM, ISUCCESS)
          CALL ET%SET("process_row", MYROW, ISUCCESS)
          CALL ET%SET("process_col", MYCOL, ISUCCESS)

!          CALL ET%SET("gpu", 1, ISUCCESS)
!          CALL ET%SET("debug", 1, ISUCCESS)
!          CALL ET%SET("verbose", 1, ISUCCESS)

!          CALL ET%SET("gpu", 1, ISUCCESS)
!          call ET%SET("use_gpu_id", int(gpuID,kind=c_int), ISUCCESS)
!          CALL ET%SET("gpu_cholesky",1, ISUCCESS)
!          CALL ET%SET("nvidia-gpu", 1, ISUCCESS)
!          CALL ET%SET("gpu_hermitian_multiply",1, ISUCCESS)

          CALL ET%SET("debug", 1, ISUCCESS)
          CALL ET%SET("verbose", 1, ISUCCESS)

          ISUCCESS = ET%SETUP()

!         CALL ET%PRINT_SETTINGS(ISUCCESS)

          CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)
          CALL ET%EIGENVECTORS(EP, W, EZ, ISUCCESS)

          CALL ELPA_DEALLOCATE(ET)
          CALL ELPA_UNINIT()

          NULLIFY(EP,EZ)
       ENDIF
#endif // ELPA
       ! call scaLAPACK routine
       IF (.NOT. LELPA) THEN
#ifdef SET_EVX_WORKSPACE
            CALL PDSYEVX_PZHEEVX_WORKSPACE( DESCA, NCLUST, N, IWORK(1), RWORK, WWORK )
#else
#ifdef gammareal
         CALL WPDSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, WWORK, &
              LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )

         RWORK(1)=1
         WWORK(1)=WWORK(1)+NCLUST*N
#else
         CALL WPZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
              WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )
         RWORK(1)=RWORK(1)+NCLUST*N
#endif
#endif
         LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
         LIWORK= MAX(1,IWORK(1))
         LRWORK= MAX(1,NINT(RWORK(1)))

         DEALLOCATE(RWORK, WWORK, IWORK)
         ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))
DOESI    RWORK = 0.0_q; WWORK = 0; IWORK = 0


#ifdef gammareal
         CALL WPDSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, WWORK, &
              LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
         CALL WPZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
              WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )
#endif
#ifndef ACC_OFFLOAD
         IF (M.NE.N) THEN
            CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // str(LIWORK) // &
               " " // str(DESCA(M_)) // "\nERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough &
               &eigenvalues found " // str(M) // " " // str(N))
         ENDIF
         IF (NZ.NE.N) THEN
            CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough eigenvectors &
               &computed " // str(M) // " " // str(N))
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: I2,IFAIL= " // str(I2) &
                  // " " // str(IFAIL(I2)))
            ENDIF
         ENDDO
         IF (INFO<0) THEN
            CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
         ELSE IF (INFO>0) THEN
            IF (MOD(INFO,2).NE.0) THEN
               CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ELSEIF (MOD(INFO/2,2).NE.0) THEN
            ! this error condition happens very often, but does not imply that
            ! the result is not usefull
            !  IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
            ELSE
               CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ENDIF
         ENDIF
#endif
       ENDIF ! IF (.NOT. LELPA) THEN

GPACC !$ACC WAIT IF(OFFLOAD_ON)
GPACC !$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
             A=Z
GPACC !$ACC END KERNELS
GPACC !$ACC EXIT DATA DELETE(Z) __IF_ASYNC__

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

#if defined(CUSOLVERMP) && ! PGI_BEFORE_25_9
       IF (.NOT. LUSENCCL) THEN
          CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(A,W) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       ENDIF
#endif

       PROFILING_STOP('pdsyevx_zheevx_desc')

     END SUBROUTINE PDSYEVX_ZHEEVX_DESC

!=======================================================================
!
!> Matrix diagonalization scaLAPACK aware version
!>
!> This version expects a proper scaLAPACK descriptor.
!> @note This routine does stop when eigenvectors cannot be
!> determined
!>
!
!=======================================================================

     SUBROUTINE PSSYEVX_CHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN, ABSTOL_)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEFS    A(:)           !< input/output matrix
       REAL(qs) W(N)           !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
       REAL(q), OPTIONAL :: ABSTOL_
! local variables
       REAL(qs) :: ABSTOL
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(qs) :: GAP(COMM%NCPU)
       REAL(qs) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2

       GDEFS, ALLOCATABLE ::    Z(:)
       REAL(qs), ALLOCATABLE::  RWORK(:) ! work array
       GDEFS, ALLOCATABLE   ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors

       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(3), WWORK(3), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
      !$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPSSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)+NCLUST*N
#else
       CALL WPCHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=RWORK(1)+NCLUST*N
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL WPSSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
       CALL WPCHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
#endif
#ifndef ACC_OFFLOAD
       IF (M.NE.N) THEN
          CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // str(LIWORK) // &
             " " // str(DESCA(M_)) // "\nERROR in subspace rotation PSSYEVX/ PCHEEVX: not enough &
             &eigenvalues found " // str(M) // " " // str(N))
       ENDIF
       IF (NZ.NE.N) THEN
          CALL vtutor%error("ERROR in subspace rotation PSSYEVX/ PCHEEVX: not enough eigenvectors &
             &computed " // str(M) // " " // str(N))
       ENDIF
       DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
             CALL vtutor%error("ERROR in subspace rotation PSSYEVX/ PCHEEVX: I2,IFAIL= " // str(I2) &
                // " " // str(IFAIL(I2)))
          ENDIF
       ENDDO
       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (MOD(INFO,2).NE.0) THEN
             CALL vtutor%error("ERROR eigenvector not converged PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
          ELSEIF (MOD(INFO/2,2).NE.0) THEN
          ! this error condition happens very often, but does not imply that
          ! the result is not usefull
!             IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PSSYEVX/ PCHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
          ELSE
             CALL vtutor%error("ERROR in diagonalization PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
          ENDIF
       ENDIF
#endif

       !$ACC WAIT IF(OFFLOAD_ON)
       !$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
       !$ACC END KERNELS
       !$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

     END SUBROUTINE PSSYEVX_CHEEVX_DESC


!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> This version expects a proper scaLAPACK descriptor.
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:)            !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ )  !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PDSYEV_ZHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PDSYEV_ZHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
!$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPDSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))
#ifdef gammareal
       CALL WPDSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSYEV_ZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PDSYEV_ZHEEVD_DESC' )
    END SUBROUTINE PDSYEV_ZHEEVD_DESC

!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> Same as #PDSYEV_ZHEEVD_DESC, but complex version regardless of gammareal
!
!=======================================================================

     SUBROUTINE PZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       COMPLEX(q) A(:)            !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ )  !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PZHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PZHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PZHEEVD failed using fallback PZHEEV"
          CALL PZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PZHEEVD_DESC' )
    END SUBROUTINE PZHEEVD_DESC

!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> This version uses a two-dimensional input array.
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEVD_DESC2(A, W, N, DESCA_, COMM, JOBIN)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:,:)          !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEF, ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

       !$ACC ENTER DATA COPYIN(Z,W) __IF_ASYNC__

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL WPZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

#ifdef CUSOLVERMP
      !$ACC KERNELS PRESENT(A, Z) __IF_ASYNC__
#endif
      A=Z
#ifdef CUSOLVERMP
      !$ACC END KERNELS
      !$ACC EXIT DATA DELETE(Z) COPYOUT(W) __IF_ASYNC__
#endif
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSYEV_ZHEEV_DESC2(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored in
! single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))


    END SUBROUTINE PDSYEV_ZHEEVD_DESC2

!=======================================================================
!
!> Matrix diagonalization using divide and conquer scaLAPACK
!>
!> Single precision version.
!
!=======================================================================

     SUBROUTINE PSSYEV_CHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
GPACC  USE moffload_struct_def
       IMPLICIT NONE

       INTEGER N                !< dimension of matrix
       GDEFS    A(:)            !< input/output matrix
       REAL(qs) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ )  !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEFS, ALLOCATABLE ::     Z(:)
       REAL(qs), ALLOCATABLE ::  RWORK(:) ! work array
       GDEFS, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PSSYEV_CHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PSSYEV_CHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
!$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPSSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL WPCHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       ! allocate optimal work space
       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL WPSSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL WPCHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEVD/ PCHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PSSYEVD/ PCHEEVD failed using fallback PSSYEV/ PCHEEV"
          CALL PSSYEV_CHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PSSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_s_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PSSYEV_CHEEVD_DESC' )
    END SUBROUTINE PSSYEV_CHEEVD_DESC

!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> Usually not very fast but robust
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)

GPACC  USE moffload_struct_def

       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:)            !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START('pdsyev_zheev_desc')

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
!$ACC ENTER DATA CREATE(Z) __IF_ASYNC__

#ifdef gammareal
       CALL WPDSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL WPZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL WPDSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL WPZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= " // str(INFO))
       ENDIF

!$ACC WAIT IF(OFFLOAD_ON)
!$ACC KERNELS PRESENT(A,Z) __IF_ASYNC__
       A=Z
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(Z) __IF_ASYNC__
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

       PROFILING_STOP('pdsyev_zheev_desc')

    END SUBROUTINE PDSYEV_ZHEEV_DESC

!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> Complex version regardless of gammareal
!
!=======================================================================

    SUBROUTINE PZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       COMPLEX(q)    A(:)      !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PZHEEV_DESC

!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> This version uses a two-dimensional input array.
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEV_DESC2(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:,:)          !< input/output matrix
       REAL(q) W(N)            !< eigenvalues
       INTEGER DESCA_( DLEN_ ) !< distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEF, ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PDSYEV_ZHEEV_DESC2


!=======================================================================
!
!> Matrix diagonalization using standard scaLAPACK diagonalization
!>
!> Single precision version
!
!=======================================================================

     SUBROUTINE PSSYEV_CHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEFS   A(:)            !< input/output matrix
       REAL(qs) W(N)           !< eigenvalues
       INTEGER :: DESCA_( DLEN_ )
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEFS, ALLOCATABLE ::     Z(:)
       REAL(qs), ALLOCATABLE :: RWORK(:) ! work array
       GDEFS, ALLOCATABLE   ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

#ifdef gammareal
       CALL PSSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PCHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PSSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PCHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEV/ PCHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PSSYEV_CHEEV_DESC


!=======================================================================
!
!> Makes an (approximately unitary) transformation matrix
!> exactly unitary (or orthogonal)
!>
!> If the orbitals (observing `<n|m>= delta_mn` )
!> are transformed by a matrix A, they will no longer be
!> exactly orthonormal
!> ~~~
!> |n'> = |n> A_nn'  (see LINCOM)
!> ~~~
!> but rather observe
!> ~~~
!> <n'| m'> =  A^+_n'm A_mm' = (A+ A)_n'm'
!> ~~~
!> The routine calculate `S=  A+ A` then determines
!>    `S= U+ U`  (U upper triangular matrix)
!> inverts U and updates the transformation matrix to
!>    `A <-  A U^-1`
!> one can easily show that
!> ~~~
!>(A U^-1)^+ (A U^-1) = 1
!> ~~~
!
!=======================================================================

     SUBROUTINE MAKE_UNITARY_DESC(A, N, DESCA )
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               !< dimension of matrix
       GDEF    A(:)            !< input/output matrix
       INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
       INTEGER :: INFO
! local variables
       GDEF, ALLOCATABLE ::     Z(:)

       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       ! calculate S= A+ A
       CALL PGGEMM( trans, 'N', N, N, N, one,  &
          A(1), 1, 1, DESCA, &
          A(1), 1, 1, DESCA, zero,  &
          Z(1), 1, 1, DESCA )

       ! determine S= U+ U  (U upper triangular matrix)
#ifdef gammareal
       CALL PDPOTRF &
#else
       CALL PZPOTRF &
#endif
          & ('U', N, Z(1), 1, 1, DESCA, INFO)

       IF (INFO/=0) THEN
          CALL vtutor%bug("internal error in VASP: MAKE_UNITARY_DESC PZPOTRF failed " // str(INFO), __FILE__, __LINE__)
       ENDIF
       ! invert U <- U^-1
#ifdef gammareal
       CALL PDTRTRI &
#else
       CALL PZTRTRI &
#endif
          & ('U','N', N, Z(1), 1, 1, DESCA, INFO)
       IF (INFO/=0) THEN
          CALL vtutor%bug("internal error in VASP: MAKE_UNITARY_DESC PZTRTRI failed " // str(INFO), __FILE__, __LINE__)
       ENDIF
! triangular update of matrix A <- A U^-1
! note that the lower part of Z is set so DGEMM does not work here
#ifdef gammareal
       CALL PDTRMM &
#else
       CALL PZTRMM &
#endif
     &                ('R', 'U', 'N', 'N' , N, N, one, &
     &                 Z, 1, 1, DESCA, A, 1, 1, DESCA)

       DEALLOCATE(Z)

     END SUBROUTINE MAKE_UNITARY_DESC

!=======================================================================
!
!> Performs a pseudo inverse of CHI1D
!>
!> CHI1D is passed as pointer to 1D array.
!>
!> The following lines can be used to reshape a 2D array into 1D form:
!> ~~~
!> USE, INTRINSIC :: iso_c_binding
!> GDEFCPTR, ALLOCATABLE   :: CHI1D(:)
!> CHI1D = RESHAPE( CHI_INV, (/MY_NROWS*MY_NCOLS/) )
!> ~~~
!
!=======================================================================
   SUBROUTINE PSEUDO_INVERSE_DESC( CHI1D, NP, DESC, THRESHHOLD, IU6 )
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF, POINTER        :: CHI1D(:)
      INTEGER              :: NP             !< dimension of submatrix
      INTEGER              :: DESC( DLEN_ )  !< distributed matrix descriptor array
      REAL(q)              :: THRESHHOLD
      INTEGER              :: IU6
      ! local
      INTEGER              :: J
      INTEGER              :: PINFO                !scaLAPACK routine info variable
      INTEGER              :: LWORK, LRWORK        !for scaLAPACK routine
      REAL(q)              :: E(NP)
      GDEF   ,ALLOCATABLE  :: EV(:)
      GDEF   ,ALLOCATABLE  :: V(:)
      GDEF, ALLOCATABLE    :: WORK(:)
#ifndef gammareal
      COMPLEX(q), ALLOCATABLE :: CWORK(:)             !needed for complex computation
#endif
      ALLOCATE(EV(SIZE(CHI1D)))
      ALLOCATE(V(SIZE(CHI1D)))
      ! initialize eigenvalues and eigenvectors
      E = zero
      EV = zero
      !
      !first get the optimal working arrays
      !
      ALLOCATE(WORK(1))
#ifdef gammareal
      !for real values only
      CALL PDSYEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), -1, PINFO)
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP=" // str(NP) // " 1st call of &
            &PDSYEV returns " // str(PINFO))
      ENDIF
#else
      !for a hermitian matrix
      ALLOCATE(CWORK(1))
      CALL PZHEEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1), &
           EV(1), 1, 1, DESC, WORK(1), -1, CWORK(1), -1, PINFO )
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP=" // str(NP) // " 1st call of &
            &PZHEEV returns " // str(PINFO))
      ENDIF
#endif

      !LWORK is now the optimal size of WORK
      LWORK=WORK(1)
      IF ( WORK(1)  /= 0 ) THEN   !if it's not zero reallocate WORK array
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
      ENDIF

      !diagonalize
#ifdef gammareal
      !real version
      CALL PDSYEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), LWORK , PINFO)
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP " // str(NP) // " 2nd call of &
            &PDSYEV returns " // str(PINFO))
      ENDIF
#else
      !in the complex version do the same with the CWORK array
      LRWORK=CWORK(1)
      IF ( LRWORK  /= 0 ) THEN   !if it's not zero reallocate WORK array
        DEALLOCATE(CWORK)
        ALLOCATE(CWORK(LRWORK))
      ENDIF
      !hermitian version
      CALL PZHEEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), LWORK, CWORK(1), LRWORK, PINFO )
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP " // str(NP) // " 2nd call of &
            &PZHEEV returns " // str(PINFO))
      ENDIF
      DEALLOCATE( CWORK )
#endif

      IF ( PINFO > 0 ) THEN
         WRITE(*,'(" Eigenvector #",I4," may not be converged" )') PINFO
      ENDIF

      IF (IU6>=0) WRITE(IU6,'(A)') ' eigenvalues of the response function'
      IF (IU6>=0) WRITE(IU6,'(2X,10F12.6)') E(1:NP)

      ! invert eigenvalues based on passed threshhold
      DO J=1,NP
         IF (ABS(E(J))<THRESHHOLD) THEN
            E(J)=0
         ELSE
            ! uncomment the following line to represent original matrix
            E(J)=1/E(J)
         ENDIF
      ENDDO

      ! form diagonal matrix with eigenvalues
      V = zero
      CALL ADD_TO_DIAGONALE_REAL(NP, V, E,  DESC )
      ! re-shape back to 2D array
       ! EIG * V,  multiply with rotation matrix
      CALL PGGEMM( 'N', 'N', NP, NP, NP, one,  &
         EV(1), 1, 1, DESC,              &
         V(1), 1, 1, DESC,         &
         zero, CHI1D(1), 1, 1, DESC )
      V = CHI1D
       ! 1/EIG * V,  multiply with rotation matrix
      CALL PGGEMM( 'N', trans, NP, NP, NP, one,  &
         EV(1), 1, 1, DESC,              &
         V(1), 1, 1, DESC,         &
         zero, CHI1D(1), 1, 1, DESC )

      DEALLOCATE( WORK )
      DEALLOCATE( EV )
      DEALLOCATE( V )
   END SUBROUTINE PSEUDO_INVERSE_DESC

!=======================================================================
!
!> Performs a Cholesky decomposition using the scalapack Descriptor DESC
!>
!=======================================================================
   SUBROUTINE PDPOTRF_PZPOTRF_DESC( CHI, NP,  DESC, JOBIN)
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF                 :: CHI(*)
      INTEGER              :: NP               !< dimension of submatrix
      INTEGER              :: DESC( DLEN_ )  !< distributed matrix descriptor array
      CHARACTER (LEN=1), OPTIONAL :: JOBIN
      ! local
      INTEGER              :: PINFO
      CHARACTER (LEN=1)    :: JOB

      ! default job is upper triangular part
      JOB = 'U'
      ! optionally pass lower triangular part
      IF ( PRESENT( JOBIN) ) JOB = JOBIN
      !call Cholesky for gamma or complex version
      CALL PGPOTRF( JOB, NP, CHI, 1, 1, DESC, PINFO)

      IF ( PINFO /=0 ) THEN
         CALL vtutor%error(" PGPOTRF failed for NP=" // str(NP) // " with error "// str(PINFO) )
      ENDIF

   END SUBROUTINE PDPOTRF_PZPOTRF_DESC

!=======================================================================
!
!> Determines the Frobenius Norm of a distributed matrix of dimension
!> N x N and returns the result in R
!> requires a MPI sum afterwards !!
!=======================================================================
   SUBROUTINE MATRIX_NORM_DESC( A, N,  DESC, R, NORMIN )
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF                   :: A(*)           !< distributed matrix
      INTEGER, INTENT(IN)    :: N              !< dimension of submatrix
      INTEGER                :: DESC( DLEN_ )  !< distributed matrix descriptor array
      REAL(q), INTENT(INOUT) :: R
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORMIN
      ! local
      CHARACTER(LEN=1)       :: NORM
      CHARACTER(LEN=1)       :: UPLO
      REAL(q)                :: WORK(1)
      REAL(q), EXTERNAL      :: PGLANXX

      NORM = 'F'
      IF( PRESENT( NORMIN) ) THEN
         IF( NORMIN == 'e' .OR. NORMIN == 'E' .OR. NORMIN == 'f' .OR. NORMIN == 'F' ) THEN
            NORM = NORMIN
         ELSE
            CALL VTUTOR%BUG( "requested norm "//NORMIN//" not implemented, only NORM = e,f supported", __FILE__, __LINE__ )
         ENDIF
      ENDIF

      UPLO='U'
      !frobenius norm of symmetric/hermitian real-valued matrix
      R = PGLANXX( NORM, UPLO, N, A, 1, 1, DESC, WORK )
      R = R**2
   END SUBROUTINE MATRIX_NORM_DESC

!=======================================================================
!
!> Performs a matrix inversion of CHI
!>
!=======================================================================
   SUBROUTINE PDGETRI_PZGETRI_DESC( INV_EPS, NP,  DESC)
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF                 :: INV_EPS(*)
      INTEGER, INTENT(IN)  :: NP               !< dimension of submatrix
      INTEGER, INTENT(IN)  :: DESC( DLEN_ )  !< distributed matrix descriptor array
      ! local
      INTEGER              :: PINFO
      GDEF, ALLOCATABLE    :: WORK(:)              !working array and eigenvectors
      GDEF, ALLOCATABLE    :: RWORK
      INTEGER,ALLOCATABLE  :: IPIV(:),IWORK(:)     !for scaLAPACK
      INTEGER              :: LWORK, LIWORK        !for scaLAPACK routine
      INTEGER              :: NPROW, NPCOL, MYROW, MYCOL
      INTEGER              :: MY_NROWS
      INTEGER              :: MY_NCOLS
      INTEGER, EXTERNAL    :: NUMROC

      PROFILING_START( 'pdgetri_pzgetri_desc')

      CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      !obtain local matrix sizes
      MY_NROWS = NUMROC(DESC(M_),DESC(MB_),MYROW,0,NPROW)
      MY_NCOLS = NUMROC(DESC(N_),DESC(NB_),MYCOL,0,NPCOL)

      !allocate IPIV
      ALLOCATE( IPIV(MY_NROWS+DESC(MB_)))

      !LU decomposition of matrix
      CALL PGGETRF( NP, NP, INV_EPS, 1, 1, DESC, IPIV, PINFO)
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("first call to PGGETRF returns " // str(PINFO) )
      ENDIF

      !obtain optimal size of working arrays
      ALLOCATE(WORK(1),IWORK(1))
      CALL PGGETRI( NP, INV_EPS, 1, 1, DESC, IPIV(1), WORK(1), -1 ,&
         IWORK(1), -1, PINFO)

      !allocate working arrays with optimum size
      RWORK = WORK(1)
      LWORK=WORK(1)
      LIWORK=IWORK(1)
      IF ( LWORK == 0 .OR. LIWORK ==0 ) THEN
         CALL vtutor%error("LWORK or LIWORK estimated from "//str([RWORK])//" is zero:" // str([LWORK, LIWORK]) )
      ENDIF

      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      ALLOCATE(WORK(LWORK))
      ALLOCATE(IWORK(LIWORK))

      !invert now
      CALL PGGETRI( NP, INV_EPS, 1, 1, DESC, IPIV(1), WORK(1), LWORK ,&
         IWORK(1), LIWORK, PINFO)
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error( "PGGETRI returns "//str(PINFO) )
      ENDIF

      !deallocate scaLAPACK's  axuillary arrays
      IF ( ALLOCATED( IPIV  )) DEALLOCATE( IPIV  )
      IF ( ALLOCATED( WORK  )) DEALLOCATE( WORK  )
      IF ( ALLOCATED( IWORK )) DEALLOCATE( IWORK )

      PROFILING_STOP( 'pdgetri_pzgetri_desc')
   END SUBROUTINE PDGETRI_PZGETRI_DESC

  END MODULE

!========================== ORTH1_DISTRI ===============================
!
!> ORTH1_DISTRI is essentially identical to #ORTH1
!>
!> It calculates a slice or stripe of the matrix
!> ~~~
!> O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
!> ~~~
!> and **then distributes the slice** in a scaLAPACK compatible manner
!> using #DISTRI_SLICE.
!> #ORTH1 requires the storage of a global NB_TOT time NB_TOT matrix
!> whereas ORTH1_DISTRI only the storage of the distributed matrix.
!>
!> Arguments are identical to #ORTH1, except for additional communicator
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)

__GPU   USE moffload_struct_def

        USE prec
        USE scala
        USE tutor, ONLY: vtutor
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        PROFILING_START('orth1_distri')

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI' )

GPACC !$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:COVL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:COVL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        COVL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC!$ACC END KERNELS
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        PROFILING_START('m_sum_g')
        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))
        PROFILING_STOP('m_sum_g')

#ifdef ELPA
        IF (LELPA) THEN
           IF (CSEL=='U') THEN
              CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),GSD%DESCA,NPOS,NPOS+NSTRIP-1)
           ELSE
              CALL vtutor%bug("internal error in VASP: ORTH1_DISTRI is only tested for CSEL=U try " &
                 // "DISTRI_SLICE_HERM_L", __FILE__, __LINE__)
           ENDIF
        ENDIF
#endif
        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),GSD%DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF

GPACC !$ACC EXIT DATA DELETE(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:COVL) IF(OFFLOAD_ON)

        PROFILING_STOP('orth1_distri')

      END SUBROUTINE ORTH1_DISTRI

!
!> Same as #ORTH1_DISTRI_HERM, this version always adds the Hermitian conjugated elements
!> to the matrix
!>
!> Costs little extra time, and is more secure, in the sense that
!> we can use the constructed matrix in places where the lower triangle is required as well
!

      SUBROUTINE ORTH1_DISTRI_HERM(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_HERM' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),GSD%DESCA,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),GSD%DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI_HERM


!========================== ORTH1_DISTRI_HERM_DESC =====================
!
!> Same as #ORTH1_DISTRI_HERM, this version requires an additional argument, the descriptor DESCA
!> of the matrix
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_HERM_DESC(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK, DESCA)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL
        INTEGER DESCA(*)

        PROFILING_START('orth1_distri_herm_desc')

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF

        PROFILING_STOP('orth1_distri_herm_desc')

      END SUBROUTINE ORTH1_DISTRI_HERM_DESC


!========================== ORTH1_DISTRI_DAVIDSON ======================
!
!> ORTH1_DISTRI_DAVIDSON is identical to #ORTH1_DISTRI
!> but also changes the intermediate COVL matrix as required in
!> the davidson.F routine (requires vector R)
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK, I
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_DAVIDSON' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        COVL(NPOS:NPOS+NSTRIP-1,1:NSTRIP)=0
        IF (MY_COMM%NODE_ME==MY_COMM%IONODE) THEN
           DO I=1,NSTRIP
              COVL(NPOS-1+I,I)=R(I)
           ENDDO
        ENDIF

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

#ifdef ELPA
        IF (LELPA) THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),GSD%DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF
#endif
        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),GSD%DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF

      END SUBROUTINE ORTH1_DISTRI_DAVIDSON

!========================== LINCOM_DISTRI ==============================
!
!> Variant of dfast::LINCOM accepting distributed matrix
!>
!> LINCOM_DISTRI does essentially the same as dfast::LINCOM but
!> accepts a distributed matrix COVL2 instead of a global
!> NB_TOT times NB_TOT matrix.
!> The required matrix COVL is reconstructed on the fly in memory
!> slice by slice.
!> The disadvantage is that the transformed wavefunctions need
!> to be stored in a temporary array.
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN,MY_COMM,NBLK)

__GPU   USE moffload_struct_def

        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q)  CF(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN)
        GDEF        COVL(*)

      ! local variables
        COMPLEX(q)  CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF_RESULT(NPROD,NIN)
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

DOESI   CF_RESULT=(0.0_q, 0.0_q); CPROF_RESULT=0; COVL_GLOBAL=0

        PROFILING_START('lincom_distri')

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI' )

GPACC !$ACC ENTER DATA CREATE(CF_RESULT,CPROF_RESULT,COVL_GLOBAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CF_RESULT,CPROF_RESULT,COVL_GLOBAL) IF(OFFLOAD_ON)
      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK

          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),GSD%DESCA,NPOS,NPOS+NOUT-1)

          PROFILING_START('m_sum_g')
          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))
          PROFILING_STOP('m_sum_g')

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),GSD%DESCA,NPOS,NPOS+NOUT-1)

        PROFILING_START('m_sum_g')
        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))
        PROFILING_STOP('m_sum_g')

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

GPACC !$ACC KERNELS PRESENT(CF,CPROF,CF_RESULT,CPROF_RESULT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CF,CPROF,CF_RESULT,CPROF_RESULT) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        CF   =CF_RESULT
        CPROF=CPROF_RESULT
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC EXIT DATA DELETE(CF_RESULT,CPROF_RESULT,COVL_GLOBAL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CF_RESULT,CPROF_RESULT,COVL_GLOBAL) IF(OFFLOAD_ON)

        PROFILING_STOP('lincom_distri')

      END SUBROUTINE LINCOM_DISTRI

!
!> similar to #LINCOM_DISTRI but add the Hermitian conjugated elements
!> in the lower triangle
!>
!> Can be used for "full" mode (MODE="F") if the transformation
!> is Hermitian and only upper triangle is set (ZPOTRI does that).
!
      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

        PROFILING_START('lincom_distri_herm')

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI_HERM' )

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),GSD%DESCA,NPOS,NPOS+NOUT-1)

          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),GSD%DESCA,NPOS,NPOS+NOUT-1)

        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

        PROFILING_STOP('lincom_distri_herm')

      END SUBROUTINE LINCOM_DISTRI_HERM


!========================== LINCOM_DISTRI_DESC =========================
!
!> LINCOM_DISTRI_DESC essentially the same as #LINCOM_DISTRI
!> but additionally requires the distributed matrix descriptor
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI_DESC(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK, DESCA)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
        INTEGER DESCA(*)                               !< distributed matrix descriptor array
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

      END SUBROUTINE LINCOM_DISTRI_DESC

!=======================================================================
!
!> Routine to call #DISTRI_SLICE without explicitly checking
!> the interface
!
!=======================================================================
      SUBROUTINE DISTRI_SLICE_NOINT(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
        USE prec
        USE scala
        IMPLICIT NONE

        INTEGER NDIM                                   !< first dimension of matrix
        INTEGER N                                      !< rank of matrix
        INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
        GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
        GDEF    A(*)                                   !< output distributed matrix
        INTEGER DESCA(*)                               !< distributed matrix descriptor array

        CALL DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

      END SUBROUTINE DISTRI_SLICE_NOINT


      SUBROUTINE BG_PDSYEVX_ZHEEVX_NOINT(COMM,A,W,N)
      USE prec
      USE scala
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               !< NxN matrix to be distributed
      GDEF    A(*)            !< input/output matrix
      REAL(q) W(N)            !< eigenvalues

      CALL BG_PDSYEVX_ZHEEVX(COMM,A(1),W,N)

      END SUBROUTINE BG_PDSYEVX_ZHEEVX_NOINT


      SUBROUTINE SCA_SHIFT_BANDS_LOW_NOINT(N, NLOW, A, SHIFT, DESCA )
      USE scala
      IMPLICIT NONE

      INTEGER N              !< actual dimension of the matrix
      INTEGER NLOW           !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA(*)       !< distributed matrix descriptor array

      CALL SCA_SHIFT_BANDS_LOW(N, NLOW, A, SHIFT, DESCA )
      END SUBROUTINE

      SUBROUTINE SCA_SHIFT_BANDS_HIGH_NOINT(N, NHIGH, A, SHIFT, DESCA )
      USE scala
      IMPLICIT NONE

      INTEGER N              !< actual dimension of the matrix
      INTEGER NHIGH          !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA(*)       !< distributed matrix descriptor array

      CALL SCA_SHIFT_BANDS_HIGH(N, NHIGH, A, SHIFT, DESCA )
      END SUBROUTINE

#else // scaLAPACK

  MODULE scala

    USE scala_struct_def
    USE prec
    USE mpimy

    IMPLICIT NONE

    CONTAINS
      SUBROUTINE PPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
        TYPE (communic) COMM
        INTEGER NDIM,N            !< NxN matrix to be distributed
        GDEF    AMATIN(NDIM,N)    !< input/output matrix
      END SUBROUTINE PPOTRF_TRTRI

      SUBROUTINE PPOTRF_POTRI (COMM, AMATIN, NDIM, N)
        TYPE (communic) COMM
        INTEGER NDIM,N            !< NxN matrix to be distributed
        GDEF    AMATIN(NDIM,N)    !< input/output matrix
      END SUBROUTINE PPOTRF_POTRI

      SUBROUTINE PDSYEVX_ZHEEVX(COMM,AMATIN,W, NDIM, N, ABSTOL, AMATIN_SCALA)
        TYPE (communic) COMM
        INTEGER NDIM,N                        !< NxN matrix to be distributed
        GDEF    AMATIN(NDIM,N)                !< input/output matrix
        REAL(q) W(N)                          !< eigenvalues
        REAL(q), OPTIONAL ::  ABSTOL          !< tolerance
        GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(N*N) !< alternative input/output matrix (scaLAPACK format)
      END SUBROUTINE PDSYEVX_ZHEEVX

      SUBROUTINE BG_PDSYEVX_ZHEEVX(COMM,A, W, N)
        TYPE (communic) COMM
        INTEGER N                             !< NxN matrix to be distributed
        GDEF    A(1,1)                        !< input/output matrix
        REAL(q) W(N)                          !< eigenvalues
      END SUBROUTINE BG_PDSYEVX_ZHEEVX

      SUBROUTINE BG_PPOTRF_TRTRI(A, N, INFO)
        TYPE (communic) COMM
        INTEGER N          !< NxN matrix to be distributed
        GDEF    A(1,1)     !< input/output matrix
        INTEGER INFO
      END SUBROUTINE BG_PPOTRF_TRTRI

      SUBROUTINE BG_PPOTRF_POTRI(A, N, INFO)
        TYPE (communic) COMM
        INTEGER N          !< NxN matrix to be distributed
        GDEF    A(1,1)     !< input/output matrix
        INTEGER INFO
      END SUBROUTINE BG_PPOTRF_POTRI

      SUBROUTINE BG_CHANGE_DIAGONALE(N,A,IU0)
        INTEGER N,IU0
        GDEF    A(*)
      END SUBROUTINE BG_CHANGE_DIAGONALE

      SUBROUTINE INIT_scala(COMM, N)
        IMPLICIT NONE
        TYPE (communic) :: COMM
        INTEGER N
      END SUBROUTINE INIT_scala

  END MODULE

      SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM, NBANDSK)
        USE prec
        USE scala

        IMPLICIT NONE

        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        TYPE(communic) :: MY_COMM
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI is not available"

      END SUBROUTINE ORTH1_DISTRI

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK)
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI_DAVIDSON is not available"

      END SUBROUTINE

      SUBROUTINE ORTH1_DISTRI_DAVMAT(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2, &
           CHAM_CORRECTION,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        GDEF      CHAM_CORRECTION(NSTRIP,NSTRIP)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI_DAVIDSON is not available"

      END SUBROUTINE ORTH1_DISTRI_DAVMAT

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(1,1)

        WRITE(0,*) "internal error in scala.F:  LINCOM_DISTRI is not available"

      END SUBROUTINE LINCOM_DISTRI

      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(1,1)

        WRITE(0,*) "internal error in scala.F:  LINCOM_DISTRI is not available"

      END SUBROUTINE LINCOM_DISTRI_HERM
#endif // scaLAPACK
