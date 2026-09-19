#include "symbol.inc"
!
!> This module the structure of the scaLAPACK descriptors
!> and the variables that control the use of scaLAPACK
!
  MODULE scala_struct_def

    USE mpimy, ONLY : communic

    IMPLICIT NONE

    !> In most cases DTYPE_ = BLOCK_CYCLIC_2D = 1
    INTEGER, PARAMETER :: DTYPE_=1
    !> The BLACS context handle, indicating the BLACS process grid A is
    !> distributed over. The context itself is global, but the handle (the
    !> integer value) may vary.
    INTEGER, PARAMETER :: CTXT_=2
    !> The number of rows in the global array A.
    INTEGER, PARAMETER :: M_=3
    !> The number of columns in the global array A.
    INTEGER, PARAMETER :: N_=4
    !> The blocking factor used to distribute the rows of the array.
    INTEGER, PARAMETER :: MB_=5
    !> The blocking factor used to distribute the columns of the array.
    INTEGER, PARAMETER :: NB_=6
    !> The process row over which the first row of the array A is distributed.
    INTEGER, PARAMETER :: RSRC_=7
    !> The process column over which the first column of the array A is distributed.
    INTEGER, PARAMETER :: CSRC_ =8
    !> The leading dimension of the local array.  LLD_A >= MAX(1,LOCr(M_A)).
    INTEGER, PARAMETER :: LLD_=9
    !> dimension of DESCA
    INTEGER, PARAMETER :: DLEN_=9

#ifdef scaLAPACK
    !> use scaLAPACK
    LOGICAL :: LscaLAPACK = .TRUE.
    !> perform parallel LU decomposition
    LOGICAL :: LscaLU = .FALSE.
    !> store subspace matrix block cyclically or not:
    !> VASP may use a "light-weighted" implementation of scaLAPACK, in the
    !> sense that the matrices are stored in the conventional (serial) manner
    !> and only when the scaLAPACK routines are called, the matrices are
    !> distributed and then treated by scaLAPACK.
    !> When LscaAWARE = T, VASP will use distributed matrices throughout.
    !> LscaAWARE is set by INIT_SCALAAWARE based on a heuristical criterion.
    !> but can be overridden by the user in the INCAR file.
    LOGICAL :: LscaAWARE = .FALSE.
    !> In case LscaAWARE was set from the INCAR file LscaAWARE_read = T
    LOGICAL :: LscaAWARE_read = .FALSE.
    !> use ELPA
#ifdef ELPA
    LOGICAL :: LELPA = .TRUE.
#else
    LOGICAL :: LELPA = .FALSE.
#endif

    TYPE scalapack_des
       !> descriptor array
       INTEGER :: DESCA(DLEN_)
       !> tells whether the descriptor is initialized
       INTEGER :: MAGIC=0
       !> dimensions of matrix for which the descriptor was created (mostly assumed square: N x N)
       INTEGER :: M, N
       !> block sizes used for the distributed matrix (mostly only NB is used)
       INTEGER :: MB, NB
       !> mpi communicator used for scaLAPACK
       INTEGER :: MPI_COMM
       !> number of PEs
       INTEGER :: NPROCS
       !> context handle of grid (usually MPI context)
       INTEGER :: ICTXT
       !> distributed matrix size on this process: number of rows
       INTEGER :: NP
       !> distributed matrix size on this process: number of cols
       INTEGER :: NQ
       !> processor grid dimensions
       INTEGER :: NPROW, NPCOL
       !> processor coordinates in procs. grid
       INTEGER :: MYROW, MYCOL
       !> allocation size of distributed matrix (local)
       INTEGER :: MALLOCPQ
#ifdef ELPA
       !> ELPA communicator
       !> (same as MPI_COMM when INIT_scala_DESC is used to set up the descriptor)
       INTEGER :: ELPA_MPI_COMM_ALL
       !> we may use ELPA only when all processes have some data
       LOGICAL :: LELPA
#endif
    END TYPE scalapack_des

    !> global scalapack descriptor used as "default" in VASP
    TYPE (scalapack_des) :: GSD

    !> default block sizes for the distributed matrices
    INTEGER, PRIVATE :: NB = 16

    !> describes a distributed (fermionic) Green's function
    !> or self-energy matrix in orbital basis
    TYPE greens_mat_des
       !> communicator for communication between all nodes
       TYPE (communic) :: COMM
       !> communicator for nodes sharing one time/frequency
       TYPE (communic) :: COMM_INTRA
       !> communicator between time/frequency points
       TYPE (communic) :: COMM_INTER
       !> maximal matrix size (required for allocation)
       INTEGER :: NROWS
       !> local row matrix size
       INTEGER :: MY_NROWS
       !> local column matrix size
       INTEGER :: MY_NCOLS
       !> BLACS context defined on top of COMM_INTAU
       INTEGER :: ICTXT
       !> matrix descriptor
       INTEGER :: DESC(DLEN_)
       !> processor grid dimensions
       INTEGER :: NPROW, NPCOL
       !> processor coordinates in procs. grid
       INTEGER :: MYROW, MYCOL
    END TYPE greens_mat_des

    CONTAINS
!
!> Routine to query the GSD%NP variable
!
    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=MAX(GSD%NP,1)
    END FUNCTION SCALA_NP

!
!> Routine to query the GSD%NQ variable
!
    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=MAX(GSD%NQ,1)
    END FUNCTION SCALA_NQ

!
!> Routine to query the NB variable (blocking factor for distributed matrices)
!
    FUNCTION GET_SCALA_NB()
      INTEGER GET_SCALA_NB
      GET_SCALA_NB=NB
    END FUNCTION GET_SCALA_NB

!
!> Routine to set the NB variable (blocking factor for distributed matrices)
!
    SUBROUTINE SET_SCALA_NB(IBLOCK)
      INTEGER IBLOCK
      NB=IBLOCK
    END SUBROUTINE SET_SCALA_NB

!
!> this subroutine sets LscaAWARE to a default value if not read from
!> INCAR and sets the block size NB for the distributed matrices
!
    SUBROUTINE INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      USE prec
      USE incar_reader, ONLY: PROCESS_INCAR
      USE reader_tags, ONLY: INCAR_F
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
      ! local variables
      REAL(q) :: STORAGE_HAM, STORAGE_WAVE

      CHARACTER(LEN=256) :: ENVVAR_VALUE
      INTEGER :: ENVVAR_STAT, IERR

#ifdef gammareal
      STORAGE_HAM =8._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU
#else
      STORAGE_HAM =16._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU
#endif
      IF (.NOT. LscaAWARE_read ) THEN
         IF (STORAGE_HAM>STORAGE_WAVE) THEN
            LscaAWARE=LscaLAPACK
         ELSE
            LscaAWARE=.FALSE.
         ENDIF
      ENDIF

      ! default block-cyclic distribution block size
#ifdef CUSOLVERMP
      NB=256
#else
      NB=16
#endif

      ! check and set block size NB from environment
      CALL GET_ENVIRONMENT_VARIABLE("VASP_SCALA_NB", ENVVAR_VALUE, STATUS=ENVVAR_STAT)

      IF (ENVVAR_STAT == 0) THEN
          READ(ENVVAR_VALUE,*,iostat=ENVVAR_STAT) NB
          IF (ENVVAR_STAT /= 0) THEN
              CALL vtutor%error("ERROR: INIT_SCALAAWARE: could not parse VASP_SCALA_NB value from env.")
          ENDIF
      ENDIF

      ! read block size from INCAR
      CALL PROCESS_INCAR(INCAR_F, "NBSCALA", NB, IERR)

    END SUBROUTINE INIT_SCALAAWARE

!
!> INIT_scala(myCOMM, N )
!> sets up a global structure for scaLAPACK communication
!> this is sufficiently flexible for most parts of VASP
!>
!> However in the more complicated cases, such as GW, it might
!> be necessary to use different contextes for instance
!> subdividing the nodes into different communication universes
!> in that case
!> INIT_scala_DESC(myCOMM, N , DESCA, GS)
!> should be used with local DESCA and GS (scalapack_des) descriptors
!
    SUBROUTINE INIT_scala(COMM, N)
      TYPE(communic) :: COMM
      INTEGER N
      ! local variables
      INTEGER :: DESCA(DLEN_)
      INTEGER :: NB

      CALL INIT_scala_DESC(COMM, N , DESCA, GSD)

    END SUBROUTINE INIT_scala

!
!> here is the context safe version (of #INIT_scala)
!
    SUBROUTINE INIT_scala_DESC(COMM, N, DESCA, GS)
      TYPE (scalapack_des) :: GS !< descriptor that should be initialized by the routine
      TYPE(communic) :: COMM     !< communicator
      INTEGER :: DESCA( DLEN_ )  !< distributed matrix descriptor array
      INTEGER :: N               !< size of the (square) matrix to be diagonalized

      ! local variables
      INTEGER :: NB

      IF (GS%MAGIC == 1789231) THEN
         IF (COMM%MPI_COMM/=GS%MPI_COMM .OR. GS%N/=N) THEN
            ! exit the BLACS context first
            CALL BLACS_GRIDEXIT(GS%ICTXT)
            ! require re-initialization
            GS%MAGIC=0
         ENDIF
      ENDIF

      IF (GS%MAGIC /= 1789231) THEN
         NB = GET_SCALA_NB()
         CALL INIT_SCALA_DES(COMM, 0, 0, N, N, NB, NB, GS)
         GS%MAGIC = 1789231
      ENDIF

      DESCA = GS%DESCA
#ifdef ELPA
! Note: when for whatever reason some communicator can not be used with ELPA,
!       the use of ELPA will be switched off completely. This is the legacy
!       behavior. At some point this should be changed to a more flexible handling.
      LELPA = LELPA.AND.GS%LELPA
#endif
    END SUBROUTINE INIT_scala_DESC


    SUBROUTINE INIT_SCALA_DES(COMM, NPROW, NPCOL, M, N, MB, NB, SDES)
      USE tutor, ONLY: vtutor
      USE string, ONLY: str

      IMPLICIT NONE

      !> communicator
      TYPE(communic) :: COMM
      !> number of processor rows and columns
      INTEGER :: NPROW, NPCOL
      !> matrix dimensions (M x N)
      INTEGER :: M, N
      !> blocking factors for block-cyclic distribution
      INTEGER :: MB, NB
      !> output: scalapack descriptor
      TYPE(scalapack_des) :: SDES

      ! local variables
      INTEGER :: INFO
      INTEGER, EXTERNAL :: NUMROC

      ! copy communicator info
      SDES%MPI_COMM = COMM%MPI_COMM
      SDES%NPROCS = COMM%NCPU

      IF (NPROW==0 .OR. NPCOL==0) THEN
         ! automatic determination of processor grid
         CALL FERMAT_RAZOR( SDES%NPROCS, SDES%NPROW, SDES%NPCOL )
      ELSE
         ! set processor grid as given (first check consistency)
         IF (NPROW*NPCOL /= SDES%NPROCS) &
            CALL vtutor%bug("INIT_SCALA_DES: NPROW * NPCOL /= NCPU: " // &
                 str(NPROW) // " " // str(NPCOL) // " " // str(SDES%NPROCS), __FILE__, __LINE__)
         SDES%NPROW = NPROW
         SDES%NPCOL = NPCOL
      ENDIF

      ! create processor grid context
      CALL M_barrier(COMM)
      CALL PROCMAP(COMM, SDES%ICTXT, 2, SDES%NPROW, SDES%NPCOL)
      ! get my row and column in processor grid
      CALL BLACS_GRIDINFO( SDES%ICTXT, SDES%NPROW, SDES%NPCOL, SDES%MYROW, SDES%MYCOL )

      ! size of the matrix to be distributed
      SDES%M = M ; SDES%N = N

      ! block-cyclic distribution with blocking factors MB and NB
      SDES%MB = MB ; SDES%NB = NB

      ! local matrix sizes of distributed matrix
      SDES%NP = NUMROC(SDES%M, SDES%MB, SDES%MYROW, 0, SDES%NPROW)   ! get number of rows on proc
      SDES%NQ = NUMROC(SDES%N, SDES%NB, SDES%MYCOL, 0, SDES%NPCOL)   ! get number of cols on proc
      ! local allocation size
      SDES%MALLOCPQ=MAX(SDES%NP*SDES%NQ,1)

      ! set up the scalapack descriptor array
      CALL DESCINIT(SDES%DESCA, SDES%M, SDES%N, SDES%MB, SDES%NB, 0, 0, SDES%ICTXT, MAX(1, SDES%NP), INFO)
      IF (INFO/=0) CALL vtutor%bug("INIT_SCALA_DES: DESCINIT INFO: " // str(INFO), __FILE__, __LINE__)

#ifdef ELPA
      ! if one of the ranks in COMM has no data we can not use ELPA
      INFO = 0; IF (SDES%NP==0 .OR. SDES%NQ==0) INFO = 1
      CALLMPI(M_sum_i(COMM, INFO, 1))

      IF (INFO==0) THEN
         SDES%LELPA = .TRUE.
         SDES%ELPA_MPI_COMM_ALL = COMM%MPI_COMM
      ELSE
         SDES%LELPA = .FALSE.
      ENDIF
#endif

    END SUBROUTINE INIT_SCALA_DES

!
!> FERMAT_RAZOR factorizes NPROCS into NPROW*NPCOL
!> with restriction that NPROW ~ NPCOL and NPROW >= NPCOL
!
    SUBROUTINE FERMAT_RAZOR( NPROCS, NPROW, NPCOL )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NPROCS
      INTEGER, INTENT(OUT) :: NPROW, NPCOL
      ! local
      INTEGER              :: NPRIMES
      INTEGER, POINTER     :: PRIMES(:) => NULL()
      INTEGER, POINTER     :: ALL_PRIMES(:) => NULL()
      INTEGER, POINTER     :: NEXP(:) => NULL()
      INTEGER              :: I,IP
      INTEGER              :: NODD
      INTEGER              :: A, B
      PROFILING_START( 'fermat_razor' )

      IF ( NPROCS == 1 ) THEN
         NPROW = 1
         NPCOL = 1
         PROFILING_STOP( 'fermat_razor' )
         RETURN
      ENDIF

      ! obtain all primes up to SQRT( NPROCS  )
      CALL ERATOSTHENES_SIEVE( NPROCS, PRIMES, ALL_PRIMES )
      NPRIMES = SIZE( PRIMES )
      ! obtain prime factorization of NPROCS
      CALL PRIME_FACTORIZATION( NPROCS, NEXP, PRIMES )

      ! first check if NPROCS is a perfect square
      IF ( IS_PERFECT_SQUARE( NPROCS ) ) THEN
         NPROW = PERFECT_SQUARE( NPROCS )
         NPCOL = PERFECT_SQUARE( NPROCS )
      ELSE
      ! the idea is to factorize NPROCS = A*B
      ! with restriction |A - B| -> 0
      ! fermat sieve works well if N is odd
      NODD = NPROCS/(PRIMES(1)**NEXP(1))

      ! in case NPROCS is even
      IF ( NODD == 1 ) THEN

         NPROW = PRIMES(1)**( NEXP(1)/2 )
         NPCOL = PRIMES(1)**( NEXP(1)/2+MOD(NEXP(1),2) )
      ! in case NPROCS is not even
      ! use fermats sieve on NODD
      ELSE

         ! perfect square nothing to do
         IF ( IS_PERFECT_SQUARE( NODD ) ) THEN
            NPROW = PERFECT_SQUARE( NODD )
            NPCOL = PERFECT_SQUARE( NODD )
         ELSE
            IF ( IS_PRIME( NODD, ALL_PRIMES ) ) THEN
               NPROW = 1
               NPCOL = NODD
            ELSE
               ! Fermat's way to factorize an odd integer:
               ! NODD = A**2- B**2 = ( A - B ) * ( A + B )
               ! and write A**2 - N = B**2
               ! and increase A by 1 until A**2 - N is a perfect square
               ! this gives B and finally the factorization
               A = CEILING( SQRT( REAL(NODD) ) )
               B = A*A - NODD
               DO
                  IF ( IS_PERFECT_SQUARE( B ) ) EXIT
                  A = A+1
                  B = A*A-NODD
               ENDDO
               NPROW = A+PERFECT_SQUARE( B )
               NPCOL = A-PERFECT_SQUARE( B )
            ENDIF
         ENDIF
         ! distribute remaining factors of 2
         ! such that NPROW is roughly NPCOL
         DO I = 1, NEXP(1)
            IF ( NPROW < NPCOL ) THEN
               NPROW = NPROW*2
            ELSE
               NPCOL = NPCOL*2
            ENDIF
         ENDDO
      ENDIF

      ENDIF ! perfect square

#ifdef CUSOLVERMP
      ! enforce NPCOL >= NPROW via check and swap (needed for e.g. NPROCS = 240 or 480)
      IF(NPROW > NPCOL) THEN
         A = NPCOL
         NPCOL = NPROW
         NPROW = A
      ENDIF
#endif

      IF ( NPROW*NPCOL /= NPROCS ) THEN
         CALL vtutor%error("ERROR, FERMAT_RAZOR failed " // str(NPROCS) // " " // str(NPROW) // " "&
             // str(NPCOL))
      ENDIF

      DEALLOCATE( PRIMES )
      DEALLOCATE( ALL_PRIMES )
      DEALLOCATE( NEXP )
      NULLIFY( PRIMES )
      NULLIFY( ALL_PRIMES )
      NULLIFY( NEXP )

      PROFILING_STOP( 'fermat_razor' )

      CONTAINS

      !> Eratosthenes_Sieve to obtain primes up to N
      SUBROUTINE ERATOSTHENES_SIEVE(N, PRIM_SQRT, PRIMES)
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: N
         INTEGER, POINTER     :: PRIM_SQRT(:)
         INTEGER              :: I, J, MAXPRIME, STAT
         LOGICAL              :: A(N)
         INTEGER, POINTER     :: PRIMES(:)

         MAXPRIME = FLOOR(SQRT(REAL(N)))
         A = .TRUE.

         DO I=2,MAXPRIME
            J = I*I
            DO
               A(J) = .FALSE.
               J = J + I ; IF ( J .GT. N ) EXIT
            ENDDO
         ENDDO !I

         ALLOCATE( PRIMES( COUNT(A)-1 ) )
         J = 1
         DO I=2,N ! SKIP 1
            IF ( .NOT. A(I) ) CYCLE
            PRIMES( J ) = I
            J = J + 1 ; IF ( J > SIZE(PRIMES) ) EXIT
         ENDDO


         J = 0
         DO I = 1, SIZE( PRIMES )
            IF ( PRIMES( I )*PRIMES( I ) > N ) THEN
               J = J+1
               EXIT
            ENDIF
            J = J+1
         ENDDO

         ALLOCATE( PRIM_SQRT( J ) )
         DO I = 1, J
            PRIM_SQRT( I ) = PRIMES( I )
         ENDDO
      END SUBROUTINE ERATOSTHENES_SIEVE

      !> obtains prime factoriztion of N
      SUBROUTINE PRIME_FACTORIZATION( N, NEXP, PRIMES )
         IMPLICIT NONE
         INTEGER          :: N
         INTEGER, POINTER :: NEXP(:)
         INTEGER          :: PRIMES(:)
         ! local
         INTEGER          :: NPRIMES
         INTEGER          :: IP
         INTEGER          :: M

         NPRIMES=SIZE( PRIMES )
         ALLOCATE( NEXP( NPRIMES ) )
         NEXP = 0
         ! 1 is not considered as prime here
         M = N
         !loop over all primes in PRIMES
         prime: DO IP = 1 , NPRIMES
            ! find exponent
            DO
               ! if this prime is a divisor
               IF( MOD( M, PRIMES( IP ) ) == 0 ) THEN
                  M = M/PRIMES( IP )
                  ! exponent of factor
                  NEXP( IP ) = NEXP( IP ) + 1
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO prime
      END SUBROUTINE PRIME_FACTORIZATION

      !> returns true if B is prime
      FUNCTION IS_PRIME( B, PRIMES )
         LOGICAL :: IS_PRIME
         INTEGER :: PRIMES(:)
         INTEGER :: B
         ! local
         INTEGER :: I
         INTEGER :: ILIMIT

         ILIMIT = SIZE( PRIMES )

         DO I = 1, ILIMIT
            IF ( B == PRIMES( I ) ) EXIT
         ENDDO

         IF ( I > ILIMIT ) THEN
            IS_PRIME = .FALSE.
         ELSE
            IS_PRIME = .TRUE.
         ENDIF
      END FUNCTION IS_PRIME

      !> returns true if B is perfect square
      FUNCTION IS_PERFECT_SQUARE( B )
         LOGICAL :: IS_PERFECT_SQUARE
         INTEGER :: B
         ! local
         INTEGER :: I
         INTEGER :: ILIMIT

         ILIMIT =  CEILING( SQRT( REAL(B) ) )
         DO I = 1, ILIMIT
            IF ( I*I == B ) EXIT
         ENDDO
         IF ( I > ILIMIT ) THEN
            IS_PERFECT_SQUARE = .FALSE.
         ELSE
            IS_PERFECT_SQUARE = .TRUE.
         ENDIF
      END FUNCTION IS_PERFECT_SQUARE

      !> returns perfect square of B
      FUNCTION PERFECT_SQUARE( B )
         INTEGER :: PERFECT_SQUARE
         INTEGER :: B
         ! local
         INTEGER :: I
         INTEGER :: ILIMIT

         ILIMIT =  CEILING( SQRT( REAL(B) ) )
         DO I = 1, ILIMIT
            IF ( I*I == B ) EXIT
         ENDDO
         PERFECT_SQUARE = I
      END FUNCTION PERFECT_SQUARE

    END SUBROUTINE FERMAT_RAZOR

!
!> PROCMAP maps processors onto a grid
!> and allows to perform scaLAPACK operations on a subset of nodes
!> taken from BLACS example code (written by Clint Whaley 7/26/94)
!> modified by gD and gK
!>
!> PROCMAP maps NPROW*NPCOL onto an VASP-Communicator
!> It is a nice system independent implementation
!>
!> @param[inout] CONTEXT (output) This integer is used by the BLACS to indicate
!>               a context. A context is a universe where messages exist and do not
!>              interact with other contexts messages. The context includes
!>              the definition of a grid, and each processs coordinates in it.
!>
!> @param[inout] MAPPING (input) Way to map processes to grid. Choices are:
!>               - `1` : row-major natural ordering
!>               - `2` : column-major natural ordering
!>
!> @param[inout] NPROW (input) The number of process rows the created grid
!>               should have.
!>
!> @param[inout] NPCOL (input) The number of process columns the created grid
!>               should have.
!
    SUBROUTINE PROCMAP(VCOMM, CONTEXT, MAPPING, NPROW, NPCOL)
      USE main_mpi, ONLY : communic
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(communic) :: VCOMM

!     .. Scalar Arguments ..
      INTEGER :: CONTEXT, MAPPING, NPROW, NPCOL
      INTEGER :: IMAP(NPROW, NPCOL)
!
      INTEGER,EXTERNAL :: BLACS_PNUM
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP
      INTEGER TMPCONTXT, NPROCS, I, J, K
      INTEGER BLACS_IDS(VCOMM%NCPU), MY_BLACS_ID
      INTEGER NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP
!
!     See how many processes there are in the system
!     NPROCS here refers to all processors being used

      NPROCS=VCOMM%NCPU
!
!     Temporarily map all processes into 1 x NPROCS grid
!
      TMPCONTXT=VCOMM%MPI_COMM
      CALL M_barrier( VCOMM )
      CALL BLACS_GRIDINIT( TMPCONTXT, 'Row', 1, NPROCS )
      CALL M_barrier( VCOMM )

!
!     If we want a row-major natural ordering
!
      CALL BLACS_GRIDINFO( TMPCONTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)

      MY_BLACS_ID= BLACS_PNUM(TMPCONTXT, MYROW_TMP, MYCOL_TMP)
      BLACS_IDS(VCOMM%NODE_ME) = MY_BLACS_ID
      CALL M_allgather_i(VCOMM, 1, BLACS_IDS)

      K=1
      IF (MAPPING .EQ. 1) THEN
         DO I = 1, NPROW
            DO J = 1, NPCOL
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
!
!     If we want a column-major natural ordering
!
      ELSE IF (MAPPING .EQ. 2) THEN
         DO J = 1, NPCOL
            DO I = 1, NPROW
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
      ELSE
         CALL vtutor%error("INIT_SCALA, PROCMAP: unknown mapping, STOP")
      END IF
!
!     Free temporary context
!
      CALL M_barrier( VCOMM )
      CALL BLACS_GRIDEXIT(TMPCONTXT)
!      WRITE (*,*) 'MAPPING',NPCOL,NPROW
!
!     Apply the new mapping to form desired context
!
      CONTEXT=VCOMM%MPI_COMM
      CALL BLACS_GRIDMAP( CONTEXT, IMAP, NPROW, NPROW, NPCOL )

!      WRITE(*,'("local map",16I3)') ((IMAP(I,J),I=1,NPROW),J=1,NPCOL)

      RETURN
    END SUBROUTINE PROCMAP


#else // scaLAPACK
!
! when scaLAPACK is not used some variables must still be defined
!
    !> use scaLAPACK
    LOGICAL :: LscaLAPACK = .FALSE.
    !> perform parallel LU decomposition
    LOGICAL :: LscaLU = .FALSE.
    !> store subspace matrix block cyclically
    LOGICAL :: LscaAWARE = .FALSE.
    !> was LscaAWARE read from file
    LOGICAL :: LscaAWARE_read = .FALSE.
    !> use ELPA
    LOGICAL :: LELPA = .FALSE.

    CONTAINS

    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=1
    END FUNCTION SCALA_NP

    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=1
    END FUNCTION SCALA_NQ

    FUNCTION GET_SCALA_NB()
      INTEGER GET_SCALA_NB
      GET_SCALA_NB=1
    END FUNCTION GET_SCALA_NB

    SUBROUTINE SET_SCALA_NB(IBLOCK)
      INTEGER IBLOCK
    END SUBROUTINE SET_SCALA_NB

    SUBROUTINE  INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      USE mpimy, ONLY : communic
      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
    END SUBROUTINE INIT_SCALAAWARE


#endif // scaLAPACK

  END MODULE scala_struct_def

#ifdef scaLAPACK
!
!> This module holds scaLAPACK related helper functions that do
!> not directly call actual scaLAPACK routines. Most concern the
!> manipulation of block-cyclic distributed matrices.
!
  MODULE scala_helper

    USE prec
    USE scala_struct_def

    CONTAINS

!
!> CHECK_scala aborts if the GSD%N is incompatible
!> to the supplied value N
!

    SUBROUTINE CHECK_scala (N, STRNG)
      USE tutor, ONLY : vtutor
      USE string, ONLY : str
      INTEGER N
      CHARACTER (LEN=*) :: STRNG
      IF (N /= GSD%N) THEN
         CALL vtutor%bug("CHECK_scala: GSD was not set up for matrix size: "//str(N) &
            // "\n CHECK_scala was called by " // STRNG, __FILE__, __LINE__)
      ENDIF
    END SUBROUTINE CHECK_scala

!
!> this routine can be used to set the workspace array for PDSYEVX and
!> PZHEEVX.
!

    SUBROUTINE PDSYEVX_PZHEEVX_WORKSPACE( DESC, NCLUST, N, LIWORK, RWORK, WWORK )
       INTEGER DESC( DLEN_ ) !< distributed matrix descriptor array
       INTEGER :: NCLUST
       INTEGER :: N
       INTEGER :: LIWORK
       REAL(q) :: RWORK(*)
       GDEF    :: WWORK(*)
       INTEGER :: NN, NB, NP0, NQ0
       INTEGER, EXTERNAL  :: NUMROC, ICEIL
       INTEGER :: NPROW, NPCOL, MYROW, MYCOL
       EXTERNAL :: BLACS_GRIDINFO

       CALL BLACS_GRIDINFO(DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)
       LIWORK=6*MAX(N,NPROW*NPCOL+1,4)

       NB = DESC(NB_)
       NN = MAX( N, NB, 2 )
       NP0 = NUMROC( NN, NB, 0, 0, NPROW )
       NQ0 = MAX( NUMROC( N, NB, 0, 0, NPCOL ), NB )
#ifdef gammareal
       RWORK(1)=1
       WWORK(1)=5*N+MAX(5*NN,NP0*NQ0)+ICEIL(N,NPROW*NPCOL)*NN+2*NB*NB+ &
          NCLUST*N
#else
       RWORK(1)=4*N+MAX( 5*NN, NP0 * NQ0 )+ICEIL( N, NPROW*NPCOL)*NN+ &
         NCLUST*N
       WWORK(1)=N + MAX((NP0 +NQ0 +NB ) * NB, 3)
#endif
    END SUBROUTINE PDSYEVX_PZHEEVX_WORKSPACE

!
!> Clean the matrix elements along the diagonal (force them to be real),
!> as it is done in the serial version
!
#define testslice
#ifndef testslice
   SUBROUTINE BG_CHANGE_DIAGONALE(N, A, IU0, NSTART, NSTOP)
      IMPLICIT NONE
      GDEF    A(*)                !< distributed matrix
      INTEGER N                   !< dimension of the matrix
      INTEGER IU0                 !< IO unit for error report
      INTEGER, OPTIONAL :: NSTART !< start with global column NSTART
      INTEGER, OPTIONAL :: NSTOP  !< stop at global column NSTOP
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER IROW, ICOL, ICOL1, ICOL2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
#ifndef gammareal
      CALL CHECK_scala (N, 'BG_CHANGE_DIAGONALE' )

      PROFILING_START('bg_change_diagonale')

      CALL BLACS_GRIDINFO(GSD%DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,GSD%DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,GSD%DESCA(NB_),MYCOL,0,NPCOL)

      ICOL1 = 1 ; IF (PRESENT(NSTART)) ICOL1 = NSTART
      ICOL2 = N ; IF (PRESENT(NSTOP) ) ICOL2 = NSTOP

      JCOL = 0

      ! loop over local columns
      DO J1=1,NQ,GSD%DESCA(NB_)
        J1RES=MIN(GSD%DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES

          IROW_JCOL = JCOL

          ! global column index
          ICOL = GSD%DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (ICOL < ICOL1 .OR. ICOL > ICOL2) GOTO 100

          ! loop over local rows
          DO I1=1,NP,GSD%DESCA(MB_)
            I1RES=MIN(GSD%DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES

              IROW_JCOL = IROW_JCOL+1

              ! global row index
              IROW = GSD%DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

              ! are we on the diagonal?
              IF( IROW == ICOL ) THEN
                 ! check whether the matrix is hermitian
                 IF (ABS(AIMAG(A(IROW_JCOL)))>1E-2_q .AND. IU0>=0) THEN
                    WRITE(IU0,*)'BG_CHANGE_DIAGONALE: WARNING: sub-space matrix is not hermitian', &
                    &                 AIMAG(A(IROW_JCOL)),GSD%DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                 ENDIF
                 ! make sure the diagonal is strictly real
                 A(IROW_JCOL) = REAL( A(IROW_JCOL) ,KIND=q)

                 ! since we have found the right row we may now move
                 ! on to the next column
                 GOTO 100
              ENDIF

            ENDDO
          ENDDO
100       JCOL = JCOL + GSD%DESCA(LLD_)
        ENDDO
      ENDDO

      PROFILING_STOP('bg_change_diagonale')
#endif //gammareal
      RETURN
    END SUBROUTINE BG_CHANGE_DIAGONALE
#else // testslice
    SUBROUTINE BG_CHANGE_DIAGONALE(N, A, IU0, NSTART, NSTOP)

__GPU USE moffload_struct_def

      IMPLICIT NONE

      GDEF    A(*)                !< distributed matrix
      INTEGER N                   !< dimension of the matrix
      INTEGER IU0                 !< IO unit for error report
      INTEGER, OPTIONAL :: NSTART !< start with global column NSTART
      INTEGER, OPTIONAL :: NSTOP  !< stop at global column NSTOP

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST, K
      INTEGER ROW, COL, ICOL1, ICOL2, A_MAX

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

#ifndef gammareal
      PROFILING_START('bg_change_diagonale')

      CALL CHECK_scala (N, 'BG_CHANGE_DIAGONALE' )

      CALL BLACS_GRIDINFO(GSD%DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,GSD%DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,GSD%DESCA(NB_),MYCOL,0,NPCOL)

      MB = GSD%DESCA(MB_ )
      NB = GSD%DESCA(NB_ )
      LD = GSD%DESCA(LLD_)

      A_MAX = NP+(NQ-1)*LD

      ICOL1 = 1 ; IF (PRESENT(NSTART)) ICOL1 = NSTART
      ICOL2 = N ; IF (PRESENT(NSTOP) ) ICOL2 = NSTOP

GPACC !$ACC PARALLEL LOOP PRESENT(A) PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:A(1:A_MAX)) &
GPOMP !$OMP PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) IF(OFFLOAD_ON)

      cols: DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         IF (COL < ICOL1 .OR. COL > ICOL2) CYCLE cols

         rows: DO I = 1, NP
            IBLCK = (I - 1) / MB
            IBASE = IBLCK * MB
            IREST = I - IBASE

            ! global row index
            ROW = MB * MYROW + IBASE * NPROW + IREST

            IF (COL == ROW) THEN
               K = I + (J - 1) * LD
#ifndef OFFLOADING
               IF ((ABS(AIMAG(A(K))) > 1E-2_q) .AND. (IU0 >= 0)) THEN
                  WRITE(IU0,*) 'BG_CHANGE_DIAGONALE: WARNING: subspace matrix is not hermitian', AIMAG(A(K)),COL
               ENDIF
#endif
               A(K) = REAL(A(K), KIND=q)
               ! move on to the next column, since we have now found
               ! the element on the diagonal in the current one
               EXIT rows
            ENDIF

         ENDDO rows
      ENDDO cols

      PROFILING_STOP('bg_change_diagonale')
#endif //gammareal
      RETURN
    END SUBROUTINE BG_CHANGE_DIAGONALE
#endif // testslice


!=======================================================================
!
! The following routines are called with descriptors DESCA.
! The routines can be called from any context or subgroup
! of MPI processed as long as DESCA is properly set up.
! They also do not refer to GSD. In that sense they are much safer to use,
! but matter of fact the caller must "know" scalapack.
! For the time being the check however that the passed matrix dimensions are
! consistent with those set up in GSD by calling CHECK_scala. This should be
! removed after more extensive testing.
!
!=======================================================================

!
!> Add diagonal elements AD (usually eigenvalues) to the supplied
!> distributed matrix
!

    SUBROUTINE ADD_TO_DIAGONALE(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(:)        !< distributed matrix
      COMPLEX(q) AD(:)       !< diagonal elements to be added to A
      INTEGER N              !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE

!
!> Similary to #SET_DIAGONALE for a real valued A
!> Also, A can be of any shape as long as compatible with DESCA
!

    SUBROUTINE SET_DIAGONALE_D(N, A, AD, DESCA )
GPACC USE moffload_struct_def
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      REAL(q) A(*)              !< distributed matrix
      REAL(q) AD(*)          !< diagonal elements to be added to A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST, K
      INTEGER MB, NB, LD, ROW, COL

      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      MB = GSD%DESCA(MB_ )
      NB = GSD%DESCA(NB_ )
      LD = GSD%DESCA(LLD_)

!$ACC PARALLEL LOOP PRESENT(A,AD) PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) __IF_ASYNC__
      cols: DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         rows: DO I = 1, NP
            IBLCK = (I - 1) / MB
            IBASE = IBLCK * MB
            IREST = I - IBASE

            ! global row index
            ROW = MB * MYROW + IBASE * NPROW + IREST

            IF (COL == ROW) THEN
               K = I + (J - 1) * LD
               A(K) = AD(ROW)
               ! move on to the next column, since we have now found
               ! the element on the diagonal in the current one
               EXIT rows
            ENDIF

         ENDDO rows
      ENDDO cols

      RETURN
    END SUBROUTINE SET_DIAGONALE_D

!
!> Similary to #SET_DIAGONALE for a complex valued A
!> Also, A can be of any shape as long as compatible with DESCA
!

    SUBROUTINE SET_DIAGONALE_Z(N, A, AD, DESCA )
GPACC USE moffload_struct_def
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(*)              !< distributed matrix
      REAL(q) AD(*)          !< diagonal elements to be added to A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST, K
      INTEGER MB, NB, LD, ROW, COL

      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      MB = GSD%DESCA(MB_ )
      NB = GSD%DESCA(NB_ )
      LD = GSD%DESCA(LLD_)

!$ACC PARALLEL LOOP PRESENT(A,AD) PRIVATE(JBLCK,JBASE,JREST,COL,I,IBLCK,IBASE,IREST,ROW,K) __IF_ASYNC__
      cols: DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         rows: DO I = 1, NP
            IBLCK = (I - 1) / MB
            IBASE = IBLCK * MB
            IREST = I - IBASE

            ! global row index
            ROW = MB * MYROW + IBASE * NPROW + IREST

            IF (COL == ROW) THEN
               K = I + (J - 1) * LD
               A(K) = AD(ROW)
               ! move on to the next column, since we have now found
               ! the element on the diagonal in the current one
               EXIT rows
            ENDIF

         ENDDO rows
      ENDDO cols

      RETURN
    END SUBROUTINE SET_DIAGONALE_Z

!
!> Same as #ADD_TO_DIAGONALE for a real valued function AD
!

    SUBROUTINE ADD_TO_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)       !< distributed matrix
      REAL(q) AD(:)      !< diagonal elements to be added to A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE_REAL

!
!> Similary to #ADD_TO_DIAGONALE for a GDEF valued function AD
!> Also, A can be of any shape as long as compatible with DESCA
!
       SUBROUTINE ADD_TO_DIAGONALE_GDEF(N, A, AD, DESCA )
         IMPLICIT NONE
         INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
         GDEF A(*)              !< distributed matrix
         REAL(q) AD(*)          !< diagonal elements to be added to A
         INTEGER N              !< dimension of the matrix
         ! local
         INTEGER NP,NQ
         INTEGER I1RES,J1RES,JCOL
         INTEGER MYROW, MYCOL, NPROW, NPCOL
         INTEGER I1,I2,J1,J2

         INTRINSIC  MIN
         INTEGER    NUMROC,IROW_JCOL
         EXTERNAL   NUMROC,BLACS_GRIDINFO

         CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
         NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
         NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

         JCOL = 0
         ! loop over local columns
         DO J1=1,NQ,DESCA(NB_)
           J1RES=MIN(DESCA(NB_),NQ-J1+1)
           DO J2=1,J1RES
             IROW_JCOL = JCOL

             ! loop over local rows
             DO I1=1,NP,DESCA(MB_)
               I1RES=MIN(DESCA(MB_),NP-I1+1)
               DO I2=1,I1RES
                 IROW_JCOL = IROW_JCOL+1
                 IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                      (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                   A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                   GOTO 100 ! next column, since we have found the right row
                 ENDIF
               ENDDO
             ENDDO
100          JCOL = JCOL + DESCA(LLD_)
           ENDDO
         ENDDO

         RETURN
       END SUBROUTINE ADD_TO_DIAGONALE_GDEF

!
!> Same as #ADD_TO_DIAGONALE_REAL, but set diagonal part
!

    SUBROUTINE SET_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)       !< distributed matrix
      REAL(q) AD(:)      !< diagonal elements to be subtracted from A
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2 <= SIZE(AD) ) THEN
                   A(IROW_JCOL) = AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                ENDIF
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_DIAGONALE_REAL

!
!> Determine diagonal elements AD of the supplied distributed matrix A
!> @attention Global sum is required after the call!
!

    SUBROUTINE DETERMINE_DIAGONALE(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(:)    !< distributed matrix
      COMPLEX(q) AD(:)   !< on exit: diagonal elements
      INTEGER N          !< dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DETERMINE_DIAGONALE

!
!> Same as #DETERMINE_DIAGONALE but for a GDEF input array (`REAL(q)` for gammaonly)
!

    SUBROUTINE DETERMINE_DIAGONALE_GDEF(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)       ! distributed matrix
      GDEF   AD(:)       ! on exit: diagonal elements
      INTEGER N          ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
      PROFILING_START( 'determine_diagonale_gdef' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      PROFILING_STOP( 'determine_diagonale_gdef' )
    END SUBROUTINE DETERMINE_DIAGONALE_GDEF

!
!> Similary to #DETERMINE_DIAGONALE but for a GDEF input array and real output array AD
!> should be called for hermitian matrices A only
!> also A can have any shape that is compatible with DESCA
!

    SUBROUTINE DETERMINE_DIAGONALE_GDEF_REAL(N, A, AD, DESCA)
      IMPLICIT NONE
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF A(*)    !< distributed matrix
      REAL(q) AD(:)   !< on exit: diagonal elements
      INTEGER N          !< dimension of the matrix
      !local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

    END SUBROUTINE DETERMINE_DIAGONALE_GDEF_REAL

    SUBROUTINE DISTRI(AMATIN, NDIM, N, A, DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEF    AMATIN(NDIM,N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("DISTRI: leading dimension of matrix too small", __FILE__, __LINE__)
      END IF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI

    !> same as #DISTRI using a single precision input array

    SUBROUTINE DISTRI_SINGLE(AMATIN, NDIM, N, A,DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM,N)
      GDEF    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("internal error in scala.F: leading dimension of matrix too small", __FILE__, __LINE__)
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE

    !> Same as #DISTRI using a single precision input array

    SUBROUTINE DISTRI_SINGLE_SINGLE(AMATIN, NDIM, N, A,DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM,N)
      GDEFS   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("internal error in scala.F: leading dimension of matrix too small", __FILE__, __LINE__)
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE_SINGLE


!=======================================================================
!
!> Merge distributed matrix into one large matrix
!
!=======================================================================

    SUBROUTINE RECON(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEF    AMATIN(NDIM, N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0._q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON

    !> Same as #RECON using a single precision output array

    SUBROUTINE RECON_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM, N)
      GDEF    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE


    !> Same as #RECON using a single precision output array

    SUBROUTINE RECON_SINGLE_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM, N)
      GDEFS   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE_SINGLE

!=======================================================================
!
!> Divide distributed matrix by eigenvalue difference
!>
!> ~~~
!>   (       |       )         (                |                )
!>   (       |       )         ( A_ij*(f_i-f_j) | A_ia*(f_i-f_a) )
!>   (  A_ij | A_ia  )         ( -------------- | -------------- )
!>   (       |       )         (    e_i -e_j    |    e_i -e_a    )
!>   (_______|_______)  ----\  (________________|________________)
!>   (       |       )  ----/  (                |                )
!>   (       |       )         ( A_ai*(f_a-f_i) | A_ab*(f_a-f_b) )
!>   ( A_ai  | A_ab  )         ( -------------- | -------------- )
!>   (       |       )         (    e_a -e_i    |    e_a -e_b    )
!>   (       |       )         (                |                )
!> ~~~
!>
!> This is usually required in standard perturbation theory.
!> The operation is only done in the occupied-unoccupied part
!> of the matrix.
!> The first version is used to construct the density matrix
!> from a given perturbation.
!> The second version is used to construct an approximately
!> unitary rotation matrix.
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF(N, A, CELTOT, FERTOT, DESCA, &
      LFINITE_T, AUXTOT )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )!< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      COMPLEX(q) CELTOT(:)  !< eigenvalues
      REAL(q)    FERTOT(:)  !< occupancies
      INTEGER N             !< dimension of the matrix
      LOGICAL :: LFINITE_T  !< finite temperature mode
      REAL(q)    AUXTOT(:)  !< derivatives of occupancies
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1

               IF ( LFINITE_T ) THEN
                  CALL SET_A_RS_AT_T( A(IROW_JCOL),        &
                     FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),        &
                     REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),q),&
                     REAL(CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q),&
                     AUXTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     AUXTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
               ELSE
                  CALL SET_A_RS_AT_T0( A(IROW_JCOL),       &
                     FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),        &
                     REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),q),&
                     REAL(CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q) )

               ENDIF

            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO
      RETURN

      CONTAINS

      ! ------------------------------------------------------------------------
      !> Gives matrix element at T=0
      !>
      !> F_R and F_S has integer values 0 and 1 only
      ! ------------------------------------------------------------------------
      SUBROUTINE SET_A_RS_AT_T0( A, FI, FJ, EI, EJ )
         GDEF, INTENT(INOUT) :: A
         REAL(q),INTENT(IN)  :: FI
         REAL(q),INTENT(IN)  :: FJ
         REAL(q),INTENT(IN)  :: EI
         REAL(q),INTENT(IN)  :: EJ
         ! local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q
         REAL(q) :: DIFF

         ! identical occupancies, clear rotation matrix
         DIFF=( FI- FJ )

         IF (ABS(DIFF)<DIFMAX) THEN
            A = 0
         ELSE
            ! (e_i - e_j)/(f_i-f_j)
            DIFF = ( EI - EJ ) / DIFF

            ! huge rotation element, should never occur except for tiny
            ! smearing parameters of the order of DIFMAX
            IF (ABS(DIFF)<DIFMAX) THEN
               A = 0
            ELSE
               A = A /DIFF
            ENDIF
         END IF

      END SUBROUTINE SET_A_RS_AT_T0

      ! ------------------------------------------------------------------------
      !> Gives matrix element at T>0
      !>
      !> FI and FJ between 0 and 1
      ! ------------------------------------------------------------------------
      SUBROUTINE SET_A_RS_AT_T( A, FI, FJ, EI, EJ, DFI, DFJ )
         GDEF, INTENT(INOUT) :: A
         REAL(q),INTENT(IN)  :: FI
         REAL(q),INTENT(IN)  :: FJ
         REAL(q),INTENT(IN)  :: EI
         REAL(q),INTENT(IN)  :: EJ
         REAL(q),INTENT(IN)  :: DFI
         REAL(q),INTENT(IN)  :: DFJ
         ! local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q

         IF (ABS( EI - EJ )<DIFMAX) THEN
            ! this was tested for small transition energies against equation below
            A = A * (DFI + DFJ)/2
         ELSE
            A =  A* ( FI- FJ ) / ( EI - EJ )
         END IF

      END SUBROUTINE SET_A_RS_AT_T

    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF

!=======================================================================
!
!> Divide matrix A by eigenvalue difference times eigenvalue (uocc_occ block
!> version)
!>
!> Divide matrix A by eigenvalue difference times eigenvalue in the
!> uocc_occ block ( i,j \in occ ; a,b \in uocc), remove uocc_uocc block
!> ~~~
!>   (       |       )         (          |           )
!>   (       |       )         (          |  e_i*A_ia )
!>   (  A_ij | A_ia  )         (   A_ij   |  -------- )
!>   (       |       )         (          |  e_i-e_a  )
!>   (_______|_______)  ----\  (__________|___________)  ---    A'
!>   (       |       )  ----/  (          |           )  ---
!>   (       |       )         ( e_i*A_ai |     0     )
!>   ( A_ai  | A_ab  )         ( -------- |    0 0    )
!>   (       |       )         ( e_i-e_a  |     0     )
!>   (       |       )         (          |           )
!> ~~~
!> In index notation:
!> ~~~
!>       A' = A_ij * ( e_i f_i - e_j f_j) / (e_i - e_j)
!>          = A_ij/2 [ f_i+f_j + (e_i+e_j) (f_i-f_j) /(e_i - e_j)]
!> ~~~
!>
!> @note At finite temperature one has 0 < f_i, f_j < 1.
!> So careful approximation of the fraction has to be done.
!> Here f is assumed to be a Fermi function.
!>
!> And add it to matrix B:
!> ~~~
!> B = B + A'
!> ~~~
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL(N, A, B, CELTOT, FERTOT, DESCA,&
      LFINITE_T, AUXTOT)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed Hamilton matrix
      GDEF    B(:)          !< distributed result matrix
      COMPLEX(q) CELTOT(:)  !< eigenvalues
      REAL(q)    FERTOT(:)  !< occupancies
      REAL(q)    AUXTOT(:)  !< occupancies
      INTEGER N             !< dimension of the matrix
      LOGICAL LFINITE_T
      REAL(q) :: EFERMI
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTRINSIC  MIN
      INTEGER    I, J
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1
               ! global column index is DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
               ! global row index DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

               J = DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
               I = DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

               ! compute matrix elements
               IF ( LFINITE_T ) THEN
                  CALL MATRIX_ELEMENT_AT_T( B(IROW_JCOL), A(IROW_JCOL),&
                      FERTOT(J), FERTOT(I), REAL(CELTOT(J),q), REAL(CELTOT(I),q), &
                      AUXTOT(I), AUXTOT(J))
               ELSE
                  CALL MATRIX_ELEMENT_AT_T0( B(IROW_JCOL), A(IROW_JCOL),&
                       FERTOT(J), FERTOT(I), REAL(CELTOT(J),q), REAL(CELTOT(I),q)  )
               ENDIF

            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      CONTAINS
      ! ------------------------------------------------------------------------
      !> Gives matrix element at T=0
      !>
      !> FI and FJ have integer values 0 and 1 only
      ! ------------------------------------------------------------------------
      SUBROUTINE MATRIX_ELEMENT_AT_T0( B, A, FJ, FI, EJ, EI )
         GDEF,INTENT(INOUT) :: B
         GDEF,INTENT(IN) :: A
         REAL(q), INTENT(IN) :: FJ, FI
         REAL(q), INTENT(IN) :: EJ, EI
         !local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q

         IF (ABS(EJ-EI)<DIFMAX) THEN
            ! not the exact solution but a sensible simple approximation at T=0
            B = B +  A/2 * ( FI +FJ )
         ELSE
            B = B + A * ( FJ*EJ - FI*EI ) / (EJ-EI)
         END IF
      END SUBROUTINE MATRIX_ELEMENT_AT_T0

      ! ------------------------------------------------------------------------
      !> Gives matrix element at T>0
      !>
      !> The matrix element at T>0 is obtained from
      !> ~~~
      !> sum_n i w_n g( iw_n) g( iw_n ) e^{-iw_n eta} = 1/2 ( sign( eta) -1)
      !> + (f_n e_n - f_m e_m ) / ( e_n - e_m )
      !> ~~~
      !> the actual Matsubara sum is
      !> ~~~
      !> 1/beta \sum_n ( iw_n ) G(iw_n + mu ) F G( iw_n + mu) e^{-i w_n eta)
      !>    fj ej - fi ei
      !> = ----------------
      !>       ej - ei
      !> ~~~
      !> This result follows from the evaluation of the series for ~ -w_n^2
      !> and the series for ~ -I w_n^3 Exp[-iw_n eta] (eta>0)
      !> The term proportional to -w_n^2 gives the result -1/2, while the term
      !> proportional to -I w_n^3 Exp[-iw_n 0-] yiels +1/2.
      !>
      !> For terms on the diagonal ei ~ ej
      !> make a Taylor exansion of fj ej  around ei
      !> ~~~
      !> fj ej ~ [fi + fi' (ej-ei)  ][ ei + (ej-ei) ]
      !>       ~ fi ei + fi' (ej -ei ) ei + fi (ej-ei)
      !> ~~~
      !> and insert this in the first term of the result giving
      !> ~~~
      !>    fi ei + fi' (ej -ei ) ei + fi (ej-ei) - fi ei
      !> = -----------------------------------------------  = fi' ei + fi
      !>                        ej - ei
      !> ~~~
      !> symmetrizing yields
      !> ~~~
      !>    fi + fj    fi' ei + fj' ej
      !> =  ------- +  ---------------
      !>       2              2
      !> ~~~
      !> to include off-diagonal terms the seond term is symmetrized as follows
      !> not this is valid since ei ~ ej holds true
      !> ~~~
      !>    fi + fj    fi' ei + fj' ej + fi' ej + fj' ei
      !> =  ------- +  ---------------------------------
      !>       2                       4
      !> ~~~
      ! ------------------------------------------------------------------------
      SUBROUTINE MATRIX_ELEMENT_AT_T( B, A, FJ, FI, EJ, EI, DFJ, DFI )
         GDEF,INTENT(INOUT) :: B
         GDEF,INTENT(IN) :: A
         REAL(q), INTENT(IN) :: FJ, FI
         REAL(q), INTENT(IN) :: EJ, EI
         REAL(q), INTENT(IN) :: DFJ, DFI
         !local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q
         ! energies close to each other
         IF (ABS(EJ-EI)<DIFMAX) THEN
            !B = B +  A* ( FI+FJ + EI * DFI + EJ*DFJ ) / 2
            ! taylor expansion if occupancies are close (tested against below)
            B = B +  A* ( FI +FJ + (EI +EJ) * (DFI + DFJ)/2) / 2
         ELSE
            !B = B + A*( FJ*EJ - FI*EI )/(EJ-EI)
            B = B + A*(( FJ*EJ - FI*EI )/(EJ-EI))
         ENDIF
       END SUBROUTINE MATRIX_ELEMENT_AT_T

    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL

!=======================================================================
!
!> Divide distribute matrix by eigenvalue difference
!>
!> This is usually required in standard perturbation theory.
!> The operation is only done in the occupied-unoccupied part
!> of the matrix.
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED(N, A, CELTOT, DESCA, LFINITE_T )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      GDEF  CELTOT(:)       !< eigenvalues
      INTEGER N             !< dimension of the matrix
      LOGICAL :: LFINITE_T
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1

              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- &
                           CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)<DIFMAX) THEN
                 A(IROW_JCOL) = 0
              ELSE
                 ! the equations above also involve
                 ! the difference between the occupancies, for the density matrix they would drop out
                 ! in other words: rotations in the (un)occupied-(un)occupied  block are irrelevant
                 ! for energies
                 ! tests indicate that inclusion of DIFCELL makes indeed little difference for energies
                 A(IROW_JCOL) = A(IROW_JCOL)/DIFCELL
                 IF (ABS(A(IROW_JCOL))>0.4_q) THEN
                    ! damp rotation matrix as in Loewdin case
                    A(IROW_JCOL)=A(IROW_JCOL)/ABS(ABS(A(IROW_JCOL)))*0.4_q
                 ENDIF
              ENDIF
              ! set diagonal to 1
              IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2==DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) THEN
                 A(IROW_JCOL)=1
              ENDIF


            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED

!=======================================================================
!
!> calculate |H_ai|^2 / (epsilon_a-epsilon_i) *(f_a-f_i)
!>
!> Both the matrix H_ai and the eigenvalues, and the occupancies are passed
!> down.
!> The calculation is only performed for |epsilon_a - epsilon_i| >DIFFMAX.
!
!=======================================================================

    SUBROUTINE SUM_SECOND_ORDER(N, A, CELTOT, FERTOT, CSUM, DESCA, LFINITE_T )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      REAL(q) CELTOT(:)     !< eigenvalues epsilon
      REAL(q) FERTOT(:)     !< occupancies
      COMPLEX(q) CSUM       !< final result |H_ai|^2 / (epsilon_a-epsilon_i)
      INTEGER N             !< dimension of the matrix
      LOGICAL :: LFINITE_T  !< finite temperature mode
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CSUM=0
      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- &
                          CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)>=DIFMAX) THEN
                 CSUM=CSUM+A(IROW_JCOL)*GCONJG(A(IROW_JCOL))/DIFCELL* &
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SUM_SECOND_ORDER

!=======================================================================
!
!> Clear occupied-occupied and unoccupied-unoccupied block of
!> a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
                     ! do nothing
               ELSE
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC

!=======================================================================
!
!> Clear occupied-unoccupied and unoccupied-occupied blocks of
!> a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC

!=======================================================================
!
!> Clear occupied-unoccupied and unoccupied-occupied blocks of
!> a square matrix using supplied fermi-weights
!>
!> Specifically the supplied matrix is multiplied by
!> ~~~
!> ferwe(ncol)*ferwe(nrow) + (1-ferwe(ncol))*(1-ferwe(nrow))
!> ~~~
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC_FER(N, A, FERTOT, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      COMPLEX(q) A(:)        !< distributed matrix
      INTEGER N              !< dimension of the matrix
      REAL(q)    FERTOT(:)   !< occupancies
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              A(IROW_JCOL)=A(IROW_JCOL)* &
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)* FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)+ &
                  (1-FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2))*(1-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)))
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC_FER

!=======================================================================
!
!> Set occupied-occupied block of a square matrix to identity matrix
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_IDENTITY(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=1._q
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_IDENTITY

!=======================================================================
!
!> Set occupied-occupied block to a square matrix
!>
!> ~~~
!> n+1        0             0           0
!> 0          n             0           0
!> 0          0             n-1         0
!> 0          0             0           1
!> ~~~
!> This will preserve the order of the states during diagonalization
!> of the density matrix.
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_INCREASING(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=LAST_OCCUPIED+1-(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_INCREASING

!=======================================================================
!
!> Conjugate the upper triangle of a matrix
!>
!> This routine constructs from the change of the density matrix
!> the unitary transformation for the orbitals.
!> One can also simply take the lower triangle and change the sign
!> in the upper.
!
!=======================================================================

    SUBROUTINE MINUS_UPPER(N, A, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array
      GDEF    A(:)          !< distributed matrix
      INTEGER N             !< dimension of the matrix
      INTEGER LAST_OCCUPIED !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) < &
                  (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                   A(IROW_JCOL)=-(A(IROW_JCOL))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MINUS_UPPER

!=======================================================================
!
!> Multiply distribute matrix by eigenvalue difference
!>
!> ~~~
!>   1/(iw - epsilon_i) and  1/(iw - epsilon_a)
!> ~~~
!> from the left and right
!
!=======================================================================

    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) CELTOTINV(:) !< multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE

!
!> Same as #MULTIPLY_BY_IVEIGENVALUE, but use only real part of transformation
!
    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) CELTOTINV(:) !< multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*REAL(CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL

!
!> Similar to #MULTIPLY_BY_IVEIGENVALUE, but use GDEF of A and real CELTOTINV
!> also A can have any shape that is compatible with DESCA
!
    SUBROUTINE MULTIPLY_A_BY_CELTOTINV_GDEF(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      GDEF A(*)               !< distributed matrix
      REAL(q) CELTOTINV(:)    !< multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

    END SUBROUTINE MULTIPLY_A_BY_CELTOTINV_GDEF

!=======================================================================
!
!> Right multiply matrix by a vector
!>
!> (most likely a BLACS routines exists to do this)
!> ~~~
!> A <- A R    or  A_ij = A_ij R_j
!> ~~~
!
!=======================================================================
!
    SUBROUTINE RIGHT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      GDEF A(:)               !< distributed matrix
      REAL(q) R(:)            !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_R

    SUBROUTINE RIGHT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) R(:)         !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_C

    SUBROUTINE LEFT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      GDEF A(:)               !< distributed matrix
      REAL(q) R(:)            !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE LEFT_MULTIPLY_BY_R

    SUBROUTINE LEFT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  !< distributed matrix descriptor array
      COMPLEX(q) A(:)         !< distributed matrix
      COMPLEX(q) R(:)         !< multiplication factors R
      INTEGER N               !< dimension of the matrix
      INTEGER LAST_OCCUPIED   !< last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE LEFT_MULTIPLY_BY_C

!=======================================================================
!
!> Distribute slice of the global matrix AMATIN to a distributed matrix A
!>
!> As input the routine expects several columns of the matrix
!> to be distributed.
!> The columns of the input matrix are COLUMN_LOW up to COLUMN_HIGH.
!> All rows are passed to the routine.
!
!=======================================================================
#ifndef testslice
    SUBROUTINE DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array
    ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('distri_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN

             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   A(IROW+(JCOL-1)*DESCA(LLD_))= &
                   AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      PROFILING_STOP('distri_slice')

      RETURN
    END SUBROUTINE DISTRI_SLICE
#else // testslice
    SUBROUTINE DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

__GPU USE moffload_struct_def

      IMPLICIT NONE

      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST
      INTEGER ROW, COL, COL_, A_MAX

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

      PROFILING_START('distri_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

      MB = DESCA(MB_ )
      NB = DESCA(NB_ )
      LD = DESCA(LLD_)
      A_MAX = NP+(NQ-1) * LD

GPACC !$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:A(1:A_MAX),AMATIN) &
GPOMP !$OMP PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) IF(OFFLOAD_ON)

      DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         ! storage column index
         COL_= COL - COLUMN_LOW + 1

         ! does this column match any of the columns of amatin,
         ! the current stripe of the global matrix
         IF (COL >= COLUMN_LOW .AND. COL <= COLUMN_HIGH) THEN

            DO I = 1, NP
               IBLCK = (I - 1) / MB
               IBASE = IBLCK * MB
               IREST = I - IBASE

               ! global row index
               ROW = MB * MYROW + IBASE * NPROW + IREST

               A(I + (J - 1) * LD) = AMATIN(ROW, COL_)
            ENDDO

         ENDIF
      ENDDO

      PROFILING_STOP('distri_slice')

      RETURN
    END SUBROUTINE DISTRI_SLICE
#endif // testslice

!
!> Distribute matrix and add hermitian conjugated elements
!
!
#ifndef testslice
    SUBROUTINE DISTRI_SLICE_HERM(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('distri_slice_herm')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part first
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                       A(IROW+(JCOL-1)*DESCA(LLD_))=GCONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                     A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      PROFILING_STOP('distri_slice_herm')

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM
#else // testslice
    SUBROUTINE DISTRI_SLICE_HERM(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

__GPU USE moffload_struct_def

      IMPLICIT NONE

      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST
      INTEGER ROW, COL, ROW_, COL_

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

      PROFILING_START('distri_slice_herm')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

      MB = DESCA(MB_ )
      NB = DESCA(NB_ )
      LD = DESCA(LLD_)

      ! conjugated part
GPACC !$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(IBLCK,IBASE,IREST,ROW,ROW_,J,JBLCK,JBASE,JREST,COL) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:A(1:NP+(NQ-1)*LD),AMATIN) &
GPOMP !$OMP PRIVATE(IBLCK,IBASE,IREST,ROW,ROW_,J,JBLCK,JBASE,JREST,COL) IF(OFFLOAD_ON)
      DO I = 1, NP
         IBLCK = (I - 1) / MB
         IBASE = IBLCK * MB
         IREST = I - IBASE

         ! global row index
         ROW = MB * MYROW + IBASE * NPROW + IREST

         ! storage row index
         ROW_= ROW - COLUMN_LOW + 1

         ! does this row match any of the columns of amatin,
         ! the current stripe (amatin) of the global matrix
         IF (ROW >= COLUMN_LOW .AND. ROW <= COLUMN_HIGH) THEN

            DO J = 1, NQ
               JBLCK = (J - 1) / NB
               JBASE = JBLCK * NB
               JREST = J - JBASE

               ! global column index
               COL = NB * MYCOL + JBASE * NPCOL + JREST

               A(I + (J - 1) * LD) = GCONJG( AMATIN(COL, ROW_) )
            ENDDO

         ENDIF
      ENDDO

      ! non-conjugated part
GPACC !$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:A(1:NP+(NQ-1)*LD),AMATIN) &
GPOMP !$OMP PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) IF(OFFLOAD_ON)
      DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         ! storage column index
         COL_= COL - COLUMN_LOW + 1

         ! does this column match any of the columns of amatin,
         ! the current stripe of the global matrix
         IF (COL >= COLUMN_LOW .AND. COL <= COLUMN_HIGH) THEN

            DO I = 1, NP
               IBLCK = (I -1) / MB
               IBASE = IBLCK * MB
               IREST = I - IBASE

               ! global row index
               ROW = MB * MYROW + IBASE * NPROW + IREST

               A(I + (J - 1) * LD) = AMATIN(ROW, COL_)
            ENDDO

         ENDIF
      ENDDO

      PROFILING_STOP('distri_slice_herm')

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM
#endif // testslice

!
!> Distribute matrix and add hermitian conjugated elements
!>
!> Similar to #DISTRI_SLICE_HERM. In this version it is assumed that the lower
!> triangle of the input matrix is set. This requires some extra precaution that
!> previously properly set elements are not overwritten. Only the lower part
!> of the input matrix is used as source.
!
    SUBROUTINE DISTRI_SLICE_HERM_L(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   !< first dimension of matrix
      INTEGER N                                      !< rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                !< starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  !< input matrix (globally known)
      GDEF    A(*)                                   !< on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         !< distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (COL>=ROW) A(IROW+(JCOL-1)*DESCA(LLD_))=GCONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW>=COL) A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM_L

!=======================================================================
!
!> Reconstruct a slice of the originally distributed matrix AMATIN
!
!=======================================================================
#ifndef testslice
    SUBROUTINE RECON_SLICE(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       !< distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('recon_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      PROFILING_STOP('recon_slice')

      RETURN
    END SUBROUTINE RECON_SLICE
#else // testslice
    SUBROUTINE RECON_SLICE(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

__GPU USE moffload_struct_def

      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       !< distributed matrix

      ! local variables
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER NP, NQ, MB, NB, LD
      INTEGER I, IBLCK, IBASE, IREST, J, JBLCK, JBASE, JREST
      INTEGER ROW, COL, COL_, A_MAX

      INTEGER NUMROC
      EXTERNAL NUMROC, BLACS_GRIDINFO

      PROFILING_START('recon_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      MB = DESCA(MB_ )
      NB = DESCA(NB_ )
      LD = DESCA(LLD_)

      A_MAX = NP+(NQ-1) * LD

GPACC !$ACC KERNELS PRESENT(AMATIN) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:AMATIN) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      AMATIN=0.0_q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      ! fill stripe of global matrix
GPACC !$ACC PARALLEL LOOP PRESENT(A,AMATIN) PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:A(1:A_MAX),AMATIN) &
GPOMP !$OMP PRIVATE(JBLCK,JBASE,JREST,COL,COL_,I,IBLCK,IBASE,IREST,ROW) IF(OFFLOAD_ON)

      DO J = 1, NQ
         JBLCK = (J - 1) / NB
         JBASE = JBLCK * NB
         JREST = J - JBASE

         ! global column index
         COL = NB * MYCOL + JBASE * NPCOL + JREST

         ! storage column index
         COL_= COL - COLUMN_LOW + 1

         ! does this column match any of the columns of amatin,
         ! the current stripe of the global matrix
         IF (COL >= COLUMN_LOW .AND. COL <= COLUMN_HIGH) THEN

            DO I = 1, NP
               IBLCK = (I - 1) / MB
               IBASE = IBLCK * MB
               IREST = I - IBASE

               ! global row index
               ROW = MB * MYROW + IBASE * NPROW + IREST

               AMATIN(ROW,COL_) = A(I + (J - 1) * LD)
            ENDDO

         ENDIF
      ENDDO

      PROFILING_STOP('recon_slice')

      RETURN
    END SUBROUTINE RECON_SLICE
#endif // testslice

!
!> Complex version of #RECON_SLICE
!

    SUBROUTINE RECON_SLICE_C(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      COMPLEX(q) AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      COMPLEX(q) A(*)                    !< distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_C

!
!> Similar to #RECON_SLICE_C but add complex conjugated elements
!> in the lower triangle in AMATIN
!
    SUBROUTINE RECON_SLICE_HERM(AMATIN,NDIM,N,A,DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_)              !< distributed matrix descriptor array
      INTEGER COLUMN_LOW, COLUMN_HIGH    !< first and last column of the matrix to be merged
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         !< on return: matrix stored in serial layout between column
      GDEF    A(*)                       !< distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON_SLICE_HERM' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     set conjugated part first restrict to lower triangle
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (ROW<COL) &
                     AMATIN(COL, ROW-COLUMN_LOW+1)= &
                        GCONJG(A(IROW+(JCOL-1)*DESCA(LLD_)))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     now the non conjugated part but restrict to upper triangle
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW<=COL) &
                      AMATIN(ROW,COL-COLUMN_LOW+1)= &
                            A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_HERM

!
!> Collect specific columns of the distribted matrix
!
    SUBROUTINE RECON_SLICE_REORDER(AMATIN, NDIM, N, A, DESCA, ICOL, NCOL)

      INTEGER DESCA( DLEN_ )             !< distributed matrix descriptor
      INTEGER NDIM                       !< first dimension of matrix AMATIN
      INTEGER N                          !< rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,NCOL)
                                         !< on return: matrix stored in serial layout between column
                                         !<   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       !< distributed matrix
      INTEGER NCOL                       !< number of columns to collect
      INTEGER ICOL(NCOL)                 !< global indices of the columns to be collected
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2,IC
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('recon_slice_reorder')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          DO IC=1,NCOL
            IF (ICOL(IC)==COL) THEN
               IROW=0
               DO I1=1,NP,DESCA(MB_)
                  I1RES=MIN(DESCA(MB_),NP-I1+1)
                  DO I2=1,I1RES
                     IROW=IROW+1
                     ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                     AMATIN(ROW,IC)= A(IROW+(JCOL-1)*DESCA(LLD_))
                  ENDDO
               ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      PROFILING_STOP('recon_slice_reorder')

      RETURN
    END SUBROUTINE RECON_SLICE_REORDER

!=======================================================================
!
!> Shift clear off diagonal components below NLOW
!
!=======================================================================

    SUBROUTINE SCA_SHIFT_BANDS_LOW(N, NLOW, A, SHIFT, DESCA )
      IMPLICIT NONE

      INTEGER N              !< actual dimension of the matrix
      INTEGER NLOW           !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array

 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2, JCOL_GLOBAL, IROW_GLOBAL

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      ! loop over local columns
      JCOL = 0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL = JCOL + 1
          JCOL_GLOBAL=(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)

          ! loop over local rows
          IROW = 0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW = IROW + 1
               IROW_GLOBAL=(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)

               IF (IROW_GLOBAL < NLOW .OR. JCOL_GLOBAL < NLOW ) THEN
                  IF (IROW_GLOBAL== JCOL_GLOBAL) THEN
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=A(IROW+ (JCOL-1)*DESCA(LLD_) )+SHIFT
                  ELSE
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=0.0_q
                  ENDIF
               ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SCA_SHIFT_BANDS_LOW

    SUBROUTINE SCA_SHIFT_BANDS_HIGH(N, NHIGH, A, SHIFT, DESCA )
      IMPLICIT NONE

      INTEGER N              !< dimension of the matrix
      INTEGER NHIGH          !< first band index not to be shifted
      GDEF    A(*)           !< distributed matrix
      REAL(q) SHIFT          !< magnitude of shift
      INTEGER DESCA( DLEN_ ) !< distributed matrix descriptor array

 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2, JCOL_GLOBAL, IROW_GLOBAL

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      ! loop over local columns
      JCOL = 0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL = JCOL + 1
          JCOL_GLOBAL=(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)

          ! loop over local rows
          IROW = 0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW = IROW + 1
               IROW_GLOBAL=(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)

               IF (IROW_GLOBAL > NHIGH .OR. JCOL_GLOBAL > NHIGH ) THEN
                  IF (IROW_GLOBAL== JCOL_GLOBAL) THEN
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=A(IROW+ (JCOL-1)*DESCA(LLD_) )+SHIFT
                  ELSE
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=0.0_q
                  ENDIF
               ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SCA_SHIFT_BANDS_HIGH

!
!> allocate an array for the Green function stored in the orbital
!> basis
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis
!> @param[in]  NKPTS_IRZ number of k-points
!> @param[in]  ISPIN     number of spin degrees of freedom
!> @param[in]  NOMEGA    number of frequency points (or 1)
!
  SUBROUTINE ALLOCATE_GREENS_MAT(GDES_MAT, CHAM_MAT,  NKPTS_IRZ, ISPIN, NOMEGA)
     USE ini
     USE tutor, ONLY: vtutor
     USE string, ONLY: str
     TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
     COMPLEX(q), POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis
     INTEGER             :: NKPTS_IRZ           ! in: number of k-points
     INTEGER             :: ISPIN               ! in: number of spin degrees of freedom
     INTEGER             :: NOMEGA              ! in: number of frequency points (or 1)
   ! local
     INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
     INTEGER,EXTERNAL :: NUMROC
     INTEGER          :: ISTAT

 ! calculate local size of matrices

     CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

     NP = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(MB_), MYROW, 0, NPROW)
     NQ = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(NB_), MYCOL, 0, NPCOL)

     ALLOCATE(CHAM_MAT(NP*NQ, NKPTS_IRZ, ISPIN, NOMEGA), STAT = ISTAT )
     IF ( ISTAT/=0 ) THEN
        CALL vtutor%error( "ALLOCATE_GREENS_MAT is not able to allocate "//&
          str( 16._q*NP*NQ*NKPTS_IRZ*ISPIN*NOMEGA/1024 ) //&
          " kB of data on MPI rank 0." )
     ENDIF
     CALL REGISTER_ALLOCATE(16._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")

  END SUBROUTINE ALLOCATE_GREENS_MAT

!
!> deallocate an array for the Green function stored in the orbital
!> basis
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis, that is deallocated
!
  SUBROUTINE DEALLOCATE_GREENS_MAT(GDES_MAT, CHAM_MAT)
    USE ini
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    COMPLEX(q), POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis

    IF (ASSOCIATED(CHAM_MAT)) THEN
       CALL DEREGISTER_ALLOCATE(16._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
       DEALLOCATE(CHAM_MAT)
       NULLIFY(CHAM_MAT)
    ELSE
       WRITE(*,*) 'internal warning in DEALLOCATE_GREENS_MAT: matrix is not allocated'
    ENDIF

  END SUBROUTINE DEALLOCATE_GREENS_MAT

!
!> same routines as ALLOCATE_GREENS_MAT
!> in imaginary time the matrix is real valued (and symmetric) at the Gamma point
!> hence once can use GDEF instead of COMPLEX(q)
!
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis
!> @param[in]  NKPTS_IRZ number of k-points
!> @param[in]  ISPIN     number of spin degrees of freedom
!> @param[in]  NOMEGA    number of frequency points (or 1)
!
  SUBROUTINE ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT,  NKPTS_IRZ, ISPIN, NOMEGA)
    USE ini
    USE tutor, ONLY: vtutor
    USE string, ONLY: str
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    GDEF, POINTER       :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis
    INTEGER             :: NKPTS_IRZ           ! in: number of k-points
    INTEGER             :: ISPIN               ! in: number of spin degrees of freedom
    INTEGER             :: NOMEGA              ! in: number of frequency points (or 1)
  ! local
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER          :: ISTAT

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(MB_), MYROW, 0, NPROW)
    NQ = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(NB_), MYCOL, 0, NPCOL)

    ALLOCATE(CHAM_MAT(NP*NQ, NKPTS_IRZ, ISPIN, NOMEGA), STAT = ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_GREENS_MAT_GDEF is not able to allocate "//&
         str( ndata*8._q*NP*NQ*NKPTS_IRZ*ISPIN*NOMEGA/1024 ) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CHAM_MAT=0
    CALL REGISTER_ALLOCATE(ndata*8._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")

  END SUBROUTINE ALLOCATE_GREENS_MAT_GDEF


!
!> same routines as DEALLOCATE_GREENS_MAT
!> in imaginary time the matrix is real valued (and symmetric) at the Gamma point
!> hence once can use GDEF instead of COMPLEX(q)
!> @param[in]  GDES_MAT  descriptor
!> @param[out] CHAM_MAT  matrix in orbital basis, that is deallocated
!
  SUBROUTINE DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT)
    USE ini
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    GDEF, POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis

    IF (ASSOCIATED(CHAM_MAT)) THEN
       CALL DEREGISTER_ALLOCATE(ndata*8._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
       DEALLOCATE(CHAM_MAT)
       NULLIFY(CHAM_MAT)
    ELSE
       WRITE(*,*) 'internal warning in DEALLOCATE_GREENS_MAT: matrix is not allocated'
    ENDIF

  END SUBROUTINE DEALLOCATE_GREENS_MAT_GDEF

!
!> destroys distributed response descriptor
!
   SUBROUTINE DESTROY_GDES_MAT( RDES )
      TYPE(greens_mat_des), POINTER :: RDES

      IF ( .NOT. ASSOCIATED( RDES ) ) RETURN

      CALL BLACS_GRIDEXIT( RDES%ICTXT )
      DEALLOCATE( RDES )
      NULLIFY( RDES )
   END SUBROUTINE DESTROY_GDES_MAT

!
!>  initalize column major distrubuted scalapack matrix descriptor
!
   SUBROUTINE RESPONSE_DES_INIT_COLMAJ( RDES, LOOPD, GDES )
      USE chi_glb
      USE minimax_struct
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      TYPE(greens_mat_des), POINTER :: RDES ! distributed response
      TYPE(loop_des)   :: LOOPD
      TYPE(greensfdes) :: GDES
      ! local
      INTEGER        :: NPROW0
      INTEGER        :: NPCOL0
      INTEGER        :: MYROW0
      INTEGER        :: MYCOL0
      INTEGER        :: LLD
      INTEGER        :: INFO
      INTEGER        :: NROWS
      INTEGER        :: NCOLS
      INTEGER, EXTERNAL :: NUMROC

      ! mark descriptor as initialized
      ALLOCATE( RDES )

      ! set communicators
      RDES%COMM_INTRA=LOOPD%COMM_IN_GROUP
      RDES%COMM_INTER=LOOPD%COMM_BETWEEN_GROUPS
      RDES%COMM=LOOPD%COMM

      !---------------------------------------------------------------------------------------
      !initialzie 1D column-major processor grid
      ! row major ordering (but column-major ordering is identical if there is only 1 column)
      NPROW0=1
      NPCOL0=LOOPD%COMM_IN_GROUP%NCPU
      ! actual size of local matrix
      NROWS  = GDES%RES_NRPLWV_ROW_DATA_POINTS
      !number of cols ( in parallel blocking factor )
      NCOLS  = GDES%RES_NRPLWV_COL_DATA_POINTS
      !
      RDES%NROWS = NROWS

      CALL PROCMAP( RDES%COMM_INTRA, RDES%ICTXT, 2, NPROW0, NPCOL0)
      CALL BLACS_GRIDINFO( RDES%ICTXT, NPROW0, NPCOL0, MYROW0, MYCOL0) !get the ids for scaLAPACK
      !define corresponding array descritpor DESCA
      !since CHI_WORK is distributed blockwisely,
      !the leading dimension of rows is simply NROWS,
      !check this
      LLD = MAX( 1, NUMROC(NROWS, NCOLS, MYROW0, 0, NPROW0) )
      IF ( LLD /= NROWS ) THEN
         CALL vtutor%bug("RESPONSE_DES_INIT_COLMAJ LLD has wrong size for NODE " // str(RDES%COMM%NODE_ME) // " " // str(LLD) // " " // str(NROWS), __FILE__, __LINE__)
      ENDIF
      CALL DESCINIT(RDES%DESC, NROWS , NROWS, NCOLS, NCOLS, 0, 0, RDES%ICTXT, &
         LLD, INFO )
      !Eventually take care for errors
      IF ( INFO /= 0 ) THEN
         CALL vtutor%bug("Error in 1D col-major initialization " // str(INFO), __FILE__,__LINE__)
      ENDIF

      !check if CHI_WORK has proper size
      RDES%MY_NROWS = NUMROC(RDES%DESC(M_),RDES%DESC(MB_),MYROW0,0,NPROW0)
      RDES%MY_NCOLS = NUMROC(RDES%DESC(N_),RDES%DESC(NB_),MYCOL0,0,NPCOL0)

   END SUBROUTINE RESPONSE_DES_INIT_COLMAJ

!
!>  initalize block cyclic distributed scalapack matrix descriptor
!
   SUBROUTINE RESPONSE_DES_INIT_BC( RDES, LOOPD, GDES )
      USE chi_glb
      USE minimax_struct
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      TYPE(greens_mat_des), POINTER :: RDES ! distributed response
      TYPE(loop_des)   :: LOOPD
      TYPE(greensfdes) :: GDES
      ! local
      INTEGER        :: LLD
      INTEGER        :: INFO
      INTEGER, EXTERNAL :: NUMROC
      INTEGER :: INU=-1
      INTEGER, PARAMETER :: ITRY = 6
      INTEGER :: NB
      ! mark descriptor as initialized
      ALLOCATE( RDES )

      ! set communicators
      RDES%COMM_INTRA=LOOPD%COMM_IN_GROUP
      RDES%COMM_INTER=LOOPD%COMM_BETWEEN_GROUPS
      RDES%COMM=LOOPD%COMM

      NB = GET_SCALA_NB()

      !---------------------------------------------------------------------------------------
      !initialzie block-cyclic 2D processor grid
      CALL FERMAT_RAZOR( RDES%COMM_INTRA%NCPU, RDES%NPROW, RDES%NPCOL )
      IF ( INU == 1 ) THEN
         WRITE(*,'(" Using a ",I3,"x",I3," processor grid for a ",I5," x ",I3,"*",I5," matrix")')&
            RDES%NPROW, RDES%NPCOL, GDES%RES_NRPLWV_ROW_DATA_POINTS,  RDES%COMM_INTRA%NCPU, GDES%RES_NRPLWV_COL_DATA_POINTS
      ENDIF

      ! update blocksize and processor grid, might be to large
      IF ( INU >= 0 ) THEN
         WRITE(*,'(" distribution is done in blocks of ",I3," x ", I3)') NB, NB
      ENDIF

      ! actual size of local matrix
      RDES%NROWS = GDES%RES_NRPLWV_ROW_DATA_POINTS

      CALL PROCMAP( RDES%COMM_INTRA, RDES%ICTXT, 2, RDES%NPROW, RDES%NPCOL)
      CALL BLACS_GRIDINFO( RDES%ICTXT, RDES%NPROW, RDES%NPCOL, RDES%MYROW, RDES%MYCOL) !get the ids for scaLAPACK

      !define corresponding array descritpor DESCA
      !since CHI_WORK is distributed blockwisely,
      !the leading dimension of rows is simply NROWS,
      !check this
      !leading dimension of local array for block cyclic distribution
      LLD = MAX( 1 , NUMROC(RDES%NROWS, NB, RDES%MYROW, 0, RDES%NPROW) )
      CALL DESCINIT(RDES%DESC, RDES%NROWS , RDES%NROWS, NB, NB, 0, 0, RDES%ICTXT, LLD, INFO )
      !Eventually take care for errors
      IF ( INFO /= 0 ) THEN
         CALL vtutor%error("Error in 2D block-cyclic initialization " // str(INFO))
      ENDIF

      !check if CHI_WORK has proper size
      RDES%MY_NROWS = NUMROC(RDES%DESC(M_),RDES%DESC(MB_),RDES%MYROW,0,RDES%NPROW)
      RDES%MY_NCOLS = NUMROC(RDES%DESC(N_),RDES%DESC(NB_),RDES%MYCOL,0,RDES%NPCOL)
   END SUBROUTINE RESPONSE_DES_INIT_BC

!
!> initialize Green's function descriptor in orbital representation
!> @param[out] GDES_MAT  pointer to generated descriptor
!> @param[in] IMAG_GRIDS time and freuency grid handle
!> @param[in] NB_TOT     maximum number of orbitals
!
  SUBROUTINE GREENS_MAT_DES_INIT( GDES_MAT, IGRID_DES, NB_TOT)
    USE chi_glb
    USE minimax_struct
    USE tutor, ONLY: vtutor
    USE string, ONLY: str

    TYPE (greens_mat_des), POINTER :: GDES_MAT   ! out: pointer to generated descriptor
    TYPE (loop_des)                :: IGRID_DES
    INTEGER                        :: NB_TOT     ! in: maximum number of orbitals
  ! local
    INTEGER :: NPROCS, NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER :: IERROR
    INTEGER :: NB_SCALA

    ALLOCATE(GDES_MAT)

    GDES_MAT%COMM=IGRID_DES%COMM
    GDES_MAT%COMM_INTRA=IGRID_DES%COMM_IN_GROUP
    GDES_MAT%COMM_INTER=IGRID_DES%COMM_BETWEEN_GROUPS

    ! few sanity tests
    IF (.NOT. LscaLAPACK) THEN
       CALL vtutor%bug("GREENS_MAT_DES_INIT requires scaLAPACK", __FILE__, __LINE__)
    ENDIF

    IF (GDES_MAT%COMM%NCPU /= GDES_MAT%COMM_INTRA%NCPU * GDES_MAT%COMM_INTER%NCPU) THEN
       CALL vtutor%bug("core number incompatible " // str(GDES_MAT%COMM%NCPU) // " " // str(GDES_MAT%COMM_INTRA%NCPU) // " " // str(GDES_MAT%COMM_INTER%NCPU), __FILE__, __LINE__)
    ENDIF

    GDES_MAT%NROWS=NB_TOT
    GDES_MAT%ICTXT=0
    NPROCS=GDES_MAT%COMM_INTRA%NCPU
!#define distribute_columns
#ifdef distribute_columns
    ! we chose to distribute the columns in a round robin fashion
    ! at least for the time being
    NPROW=1
    NPCOL=NPROCS

    CALL M_barrier( GDES_MAT%COMM_INTRA )
    CALL PROCMAP( GDES_MAT%COMM_INTRA, GDES_MAT%ICTXT, 2, NPROW, NPCOL)

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(NB_TOT,NB_TOT,MYROW,0,NPROW)   ! get number of rows on proc
    NQ = NUMROC(NB_TOT,1, MYCOL,0, NPCOL)      ! get number of cols on proc

    CALL DESCINIT(GDES_MAT%DESC, NB_TOT,NB_TOT, NB_TOT,1, 0,0, GDES_MAT%ICTXT,MAX( 1, NP ), IERROR ) ! setup descriptor
    IF (IERROR.NE.0) THEN
      CALL vtutor%bug("DESC, DESCINIT, INFO: " // str(IERROR), __FILE__, __LINE__)
    ENDIF
#else
! determine processor grid, according to naive guess
    CALL FERMAT_RAZOR( NPROCS, NPROW, NPCOL )

    CALL M_barrier( GDES_MAT%COMM_INTRA )
    CALL PROCMAP( GDES_MAT%COMM_INTRA, GDES_MAT%ICTXT, 2, NPROW, NPCOL)

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    ! set blocking factor compatible to scala.F
    ! this allows to call INIT_scala and use simple functions
    ! defined in scala.F
    NB_SCALA = GET_SCALA_NB()

    NP = NUMROC(NB_TOT,NB_SCALA, MYROW,0,NPROW)   ! get number of rows on proc
    NQ = NUMROC(NB_TOT,NB_SCALA, MYCOL,0,NPCOL)   ! get number of cols on proc

    CALL DESCINIT(GDES_MAT%DESC, NB_TOT,NB_TOT, NB_SCALA,NB_SCALA, 0,0, GDES_MAT%ICTXT, MAX( 1, NP ), IERROR ) ! setup descriptor
    IF (IERROR.NE.0) THEN
      CALL vtutor%bug("DESC, DESCINIT, INFO: " // str(IERROR) // " " // str(NB_TOT) // " " // str(NB_SCALA) // " " // str(NP),__FILE__, __LINE__)
    ENDIF
#endif
    GDES_MAT%MY_NROWS = NP
    GDES_MAT%MY_NCOLS = NQ
  END SUBROUTINE GREENS_MAT_DES_INIT

!
!> small routine to dump a distributed matrix
!> the descriptor in DESCA must properly describe the matrix
!> since RECON_SLICE calls check (GDEF version)
!
!> @param[in]   STRING          string to dump
!> @param[in]   CHAM_DISTRI     distributed matrix
!> @param[in]   GDES_MAT        descriptor
!> @param[in]   NB_TOT          leading dimension of distributed matrix
!> @param[in]   IU_             unit to write to (not dump for IU<0)
!> @param[in]   IUNIT           unit to write to (not dump for IU<0)
!
  SUBROUTINE DUMP_GREENS_HAM_GDEF( STRING, CHAM_DISTRI, NB_TOT, GDES_MAT, IU_,IUNIT)
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    GDEF              :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des):: GDES_MAT       ! in: descriptor
    INTEGER           :: NB_TOT
    INTEGER           :: IU_                ! unit to write to (not dump for IU<0)
    INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
  ! local
    INTEGER :: NDUMP=16
    GDEF,ALLOCATABLE   :: CHAM(:,:)
    INTEGER IU
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW
    LOGICAL :: LGAMMA_PRINT = .FALSE.

    IF ( PRESENT( IUNIT ) ) THEN
       IU = IUNIT
       NDUMP =  NB_TOT
    ELSE
       IU = IU_
    ENDIF
    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM)))


    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
    IF ( PRESENT( IUNIT ) ) THEN
       DO N2 = 1, NDUMP
          WRITE(IU,10)N1+COLUMN_LOW-1, N2+COLUMN_LOW-1, CHAM(N1+COLUMN_LOW-1,N2)
       ENDDO
    ELSE
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1+COLUMN_LOW-1,N2) ,KIND=q) ,N2=1,NDUMP)
    ENDIF
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
IF ( LGAMMA_PRINT ) THEN
    DO N1=1,NDUMP
    IF ( .NOT. PRESENT( IUNIT ) ) THEN
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1+COLUMN_LOW-1,N2)),N2=1,NDUMP)
    ENDIF
    ENDDO
    WRITE(IU,*)
ENDIF
#endif
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(1I3,3X,40F10.4)
2   FORMAT(1I3,3X,40F10.4)
10   FORMAT(2I4,2F20.10)

  END SUBROUTINE DUMP_GREENS_HAM_GDEF

!
!> small routine to dump a distributed matrix
!> the descriptor in DESCA must properly describe the matrix
!> since RECON_SLICE calls check
!> @param[in]   STRING          string to dump
!> @param[in]   CHAM_DISTRI     distributed matrix
!> @param[in]   GDES_MAT        descriptor
!> @param[in]   NB_TOT          leading dimension of distributed matrix
!> @param[in]   IU_             unit to write to (not dump for IU<0)
!> @param[in]   IUNIT           unit to write to (not dump for IU<0)
!
  SUBROUTINE DUMP_GREENS_HAM( STRING, CHAM_DISTRI, NB_TOT, GDES_MAT, IU_,IUNIT)
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    COMPLEX(q)        :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des):: GDES_MAT       ! in: descriptor
    INTEGER           :: NB_TOT
    INTEGER           :: IU_               ! unit to write to (not dump for IU<0)
    INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
  ! local
    INTEGER           :: IU                ! unit to write to (not dump for IU<0)
    INTEGER, PARAMETER :: NDUMP=8
    COMPLEX(q),ALLOCATABLE  :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW

    IF ( PRESENT( IUNIT ) ) THEN
       IU = IUNIT
    ELSE
       IU = IU_
    ENDIF
    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE_C(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_z(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM)))

    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1,N2+COLUMN_LOW-1) ,KIND=q) ,N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
    DO N1=1,NDUMP
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1,N2+COLUMN_LOW-1)),N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#endif
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(1I2,3X,40F16.8)
2   FORMAT(1I2,3X,40F16.8)

  END SUBROUTINE DUMP_GREENS_HAM

!
!> Determines the Frobenius Norm of a distributed matrix of dimension
!> N x N and returns the result in R
!> requires a MPI sum afterwards !!
!> @param[in]   A               distributed matrix
!> @param[in]   N               dimension of sub matrix of A for which norm is calculated
!> @param[in]   RDES2D          descriptor
!> @param[out]  R               squared Frobenius norm
!> @param[in]   NORMIN          optional argument that specifies which norm is calculated
!
    SUBROUTINE MATRIX_NORM_DESC_OLD( A, N,  RDES2D, R, NORMIN )
       USE tutor, ONLY: vtutor
       USE string, ONLY: str
       GDEF, INTENT(IN)       :: A(:)
       INTEGER, INTENT(IN)    :: N               !< dimension of submatrix
       TYPE( greens_mat_des ) :: RDES2D
       REAL(q), INTENT(INOUT) :: R
       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORMIN
       INTEGER               :: I, J
       INTEGER               :: NCOLS, NROWS
       INTEGER, EXTERNAL     :: NUMROC
       ! local
       ! padestrian computation of squared Frobenius norm of A:
       !CALL BLACS_GRIDINFO( RDES2D%ICTXT, RDES2D%NPROW, RDES2D%NPCOL, RDES2D%MYROW, RDES2D%MYCOL )
       NROWS = NUMROC(N,RDES2D%DESC(MB_), RDES2D%MYROW, 0, RDES2D%NPROW )
       NCOLS = NUMROC(N,RDES2D%DESC(NB_), RDES2D%MYCOL, 0, RDES2D%NPCOL )

       R = 0
       DO J = 1, NCOLS
          DO I = 1, NROWS
             R = R + A( I + (J-1)*RDES2D%MY_NROWS )*GCONJG( A( I + (J-1)*RDES2D%MY_NROWS ) )
          ENDDO
       ENDDO
       ! accumulate results from each rank
       CALLMPI( M_sum_d(RDES2D%COMM_INTRA, R , 1 ))
    END SUBROUTINE MATRIX_NORM_DESC_OLD

  END MODULE scala_helper
#endif // scaLAPACK
