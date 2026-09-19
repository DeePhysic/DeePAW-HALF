#include "symbol.inc"
#undef STOP

#if defined(MPI) || defined(MPI_CHAIN)

!=======================================================================
!
! in most MPI implementations the collective communcation
! is slow, hence by default we avoid them
! if you want to use them define use_collective (here or in the makefile)
!
! use_collective use MPI_alltoall
! avoid_async    tries to post syncronised send and read operations
!                such that collisions are avoided, this is usually slower
! PROC_GROUP     does communication is block wise in a group of
!                roughly PROC_GROUP processors
!                this reduces collisions
!
! in addition our own implementation of the collective communication
! routines allow
!=======================================================================

!#define use_collective
!#define avoid_async

#ifdef use_collective
#define use_collective_sum
#endif

#ifndef MPI_BLOCK
#define MPI_BLOCK   1000
!
! blocking over group of processors
! on a Gigabit ethernet this can improve performance substantially
! without hardware flowcontrol
! when hardware flowcontrol is enabled PROC_GROUP is usually not
! required
!
!#define PROC_GROUP 4
#endif
!
! OPENMPI seems to have a bcast bug
! whenever a node initiates a bcast the master node seems
! to return immediately, whereas all the other nodes seem
! to wait until the master node initiates ANOTHER MPI call
! a simple work around is to initiate a barrier command after
! every single bcast
! to do this set MPI_bcast_with_barrier the makefile
! or undocument this line
#define MPI_bcast_with_barrier
!#define MPI_avoid_bcast
!
! alternatively the main VASP code can initiate a barrier after
! a block of bcast calls
! this might be more efficient since less barrier are initiated
! this requires to define MPI_barrier_after_bcast
! in the makefile.
! The corresponding barriers have been inserted in
! wave_high.F and wave_mpi.F

!
! To fall back to blocking bcast in M_ibcast_z/d_from
#define MPI_avoid_ibcast_from

!
! To fall back to blocking reduce in M_ireduce_z/d_to
#define MPI_avoid_ireduce_to

!
! In CUDA-aware mode MPI_Ireduce is not supported yet
#define no_cuda_aware_ireduce

!=======================================================================
!
!> MPI communication routines for VASP
!>
!> all communication should be done using this interface to allow
!> adaption of other communication routines
!> routines were entirely rewritten by Kresse Georg,
!> but functionallity is similar to a module written by Peter
!> Lockey at Daresbury
!
!======================================================================
      MODULE mpimy

__GPU USE moffload_struct_def

      USE prec

#ifdef USENCCL
GPACC USE nccl, ONLY : ncclComm
GPOMP USE iso_c_binding, ONLY : ncclComm => c_ptr
#endif
#ifdef qd_emulate
      USE qdmodule, ONLY : qd_real, qd_complex
#endif

! Nbranch is the number of branches at each node in the gsum routines
! two should be always fine
!
      INTEGER Nbranch
      PARAMETER( Nbranch=2 )

      TYPE communic
!only COMM
        INTEGER MPI_COMM          !< MPI_Communicator
        INTEGER NODE_ME           !< node id starting from 1 ... NCPU
        INTEGER IONODE            !< node which has to do IO (set to 0 for no IO)
        INTEGER NCPU              !< total number of proc in this communicator
#ifdef USENCCL
        TYPE (ncclComm) NCCL_COMM
        LOGICAL :: LUSENCCL = .FALSE.
        LOGICAL :: LGRPNCCL = .TRUE.
#endif
      END TYPE

! Standard MPI include file.
! I would like to have everything in the header but "freaking" SGI
! compiler can not handle this, thus I have to use an include file
      INCLUDE "mpif.h"

#ifdef qd_emulate
!> create some emulated MPI_real16 types
      INTEGER :: M_qd_type
#endif
!> create operation for collective summation of quadruples
      INTEGER :: M_sum_qd_op

#ifdef MPI_INPLACE
! FIXME: check what M_cycle_d does with the intermediate array.
! REDUCEALL routines no longer use intermediate arrays
      INTEGER,PARAMETER ::  NDTMP=MPI_BLOCK
!> workspace for real operations
      REAL(q),SAVE    :: DTMP_m(NDTMP)
#else
! There are no global local sum routines in MPI, thus some workspace
! is required to store the results of the global sum
      INTEGER,PARAMETER :: NZTMP=MPI_BLOCK/2
      INTEGER,PARAMETER :: NDTMP=MPI_BLOCK
      INTEGER,PARAMETER :: NITMP=MPI_BLOCK
      INTEGER,PARAMETER :: NI8TMP=MPI_BLOCK
      INTEGER,PARAMETER :: NLTMP=MPI_BLOCK
! workspace for integer, complex, and real
      COMPLEX(q),SAVE   :: ZTMP_m(NZTMP)
      REAL(q),SAVE      :: DTMP_m(NDTMP)
      INTEGER,SAVE      :: ITMP_m(NITMP)
      INTEGER(qi8),SAVE :: I8TMP_m(NI8TMP)
      LOGICAL,SAVE      :: LTMP_m(NLTMP)

#ifndef IFC
! TODO: If the code is ever compiled with 4-byte floats, the 8-byte integer work array I8TMP_m would
! exceed the ZTMP_m buffer. So should we not switch ZTMP_m here with I8TMP_m? Or, alternatively, set
! NI8TMP = MPI_BLOCK / 2 instead?
      EQUIVALENCE (DTMP_m,ZTMP_m)
      EQUIVALENCE (ITMP_m,ZTMP_m)
      EQUIVALENCE (I8TMP_m,ZTMP_m)
      EQUIVALENCE (LTMP_m,ZTMP_m)
#endif
#endif

      CONTAINS

!*******************************************************************
!>computes local index of a cyclic distribution of an array
!>distributed among NCPU ranks
!>   @param[in]     NCPU number of MPI ranks
!>   @param[in]     NODE_ME rank id
!>   @param[out]    M index position of rank NODE_ME
!*******************************************************************

   FUNCTION LOCAL_INDEX_CYCLIC( NCPU, NODE_ME, M )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER  :: NODE_ME
      INTEGER  :: M, L
      INTEGER  :: NCPU
      INTEGER  :: IND
      INTEGER  :: LOCAL_INDEX_CYCLIC
      INTEGER  :: IBLOCK
      IF ( NCPU == 1 ) THEN
         LOCAL_INDEX_CYCLIC = M
         RETURN
      ENDIF
      ! number of ranks have to be at least 1
      IF ( NCPU < 1  ) THEN
         CALL vtutor%bug('LOCAL_INDEX_CYCLIC reports negative number of MPI ranks: NCPU < 1 ' &
            // str(NCPU), __FILE__, __LINE__)
      ENDIF
      ! in case there are more ranks than data available
      IBLOCK = MOD( M-1, NCPU )+1        ! local index array
      IND = FLOOR( (M-1)/REAL(NCPU) )+1  ! process id handling M
      ! if out of bounds
      IF ( IBLOCK /= NODE_ME ) THEN
         IND = -1000
      ENDIF
      LOCAL_INDEX_CYCLIC = IND
   END FUNCTION LOCAL_INDEX_CYCLIC


!
!> Common error handling for all MPI calls.
!>
!> The macro wraps the subroutine with the same name, to make it work with the
!> preprocessor, we choose a different capitalization for the subroutine than
!> for the macro.
!
#define CHECK_MPI_ERROR(ierror,vasp_func,mpi_func) CALL check_mpi_error(ierror, vasp_func, mpi_func, __FILE__, __LINE__)
      SUBROUTINE check_mpi_error(IERROR, VASP_FUNC, MPI_FUNC, FILENAME, LINE)
         USE string, ONLY: str
         USE tutor, ONLY: vtutor_bug
         INTEGER, INTENT(IN) :: IERROR, LINE
         CHARACTER(LEN=*), INTENT(IN) :: VASP_FUNC, MPI_FUNC, FILENAME
         IF (IERROR /= MPI_success) &
            CALL vtutor_bug(VASP_FUNC // ": " // MPI_FUNC // " returns: " // str(IERROR), FILENAME, LINE)
      END SUBROUTINE check_mpi_error

!
!> define an MPI operation for emulated quad type with the same
!> effect as the intrinsically defined MPI_SUM
!
      SUBROUTINE M_sum_qd_op_function( INVEC, INOUTVEC, LEN, DATATYPE )
#ifdef qd_emulate
         USE qdmodule, ONLY : operator(+)
#endif

         QDPREAL :: INVEC( LEN ), INOUTVEC( LEN )
         INTEGER :: LEN
         INTEGER :: DATATYPE
         ! local
         INTEGER :: I
#ifdef qd_emulate
         QDPREAL :: BUFF( LEN )

         DO I = 1, LEN
            BUFF( I ) = INOUTVEC( I )
         ENDDO
         DO I = 1, LEN
            BUFF( I ) = BUFF( I ) + INVEC( I )
         ENDDO
         DO I = 1, LEN
            INOUTVEC( I ) = BUFF( I )
         ENDDO
#else
         DO I = 1, LEN
            INOUTVEC( I ) = INVEC( I ) + INOUTVEC( I )
         ENDDO
#endif
      END SUBROUTINE M_sum_qd_op_function

!----------------------------------------------------------------------
!
!> initialise the basic communications (number of nodes, determine ionode)
!
!----------------------------------------------------------------------

      SUBROUTINE M_init( COMM )
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror

#ifndef _OPENMP
      CALL MPI_init( ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_init')
#else // _OPENMP
#ifdef OMP_OFFLOAD
      CALL MPI_init_thread( MPI_THREAD_MULTIPLE, i, ierror )
      IF ( i /= MPI_THREAD_MULTIPLE .OR. ierror /= MPI_success ) THEN
#else
      CALL MPI_init_thread( MPI_THREAD_SINGLE, i, ierror )
      IF ( i /= MPI_THREAD_SINGLE .OR. ierror /= MPI_success ) THEN
#endif
         CALL vtutor%bug('M_init: MPI_init_thread returns: ', __FILE__, __LINE__)
      ENDIF
#endif // _OPENMP
!
! initial communicator is world wide
! set only NCPU, NODE_ME and IONODE
! no internal setup done at this point
!
      COMM%MPI_COMM= MPI_comm_world

      CALL MPI_comm_rank( COMM%MPI_COMM, COMM%NODE_ME, ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_comm_rank')
      COMM%NODE_ME= COMM%NODE_ME+1

      CALL MPI_comm_size( COMM%MPI_COMM, COMM%NCPU , ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_comm_size')

      COMM%IONODE = 1

#ifdef qd_emulate
      ! define a MPI_sum operation for qd_type
      ! explain to MPI how type M_qd_type is defined
      ! a qd_type contains 4 contiguous double-precision arrays (real*8)
      CALL MPI_type_contiguous( 4, MPI_real8, M_qd_type, ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_type_contigous')
      ! commit M_qd_type
      CALL MPI_type_commit( M_qd_type, ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_type_commit')
      ! now create the qd_type-summation user operation
      CALL MPI_op_create( M_sum_qd_op_function, .TRUE., M_sum_qd_op, ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_op_create')
#else
      ! my own MPI summation of quadruple precision variables
      CALL MPI_op_create( M_sum_qd_op_function, .TRUE., M_sum_qd_op, ierror )
      CHECK_MPI_ERROR(ierror,'M_init','MPI_op_create')
#endif
      END SUBROUTINE M_init


!----------------------------------------------------------------------
!
!> creates a 2 dimensional cartesian topology and a communicator along
!> rows and columns of the process matrix
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide( COMM, NPAR, COMM_INTER, COMM_INB, reorder)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM, COMM_INTER, COMM_INB, COMM_CART
      INTEGER NPAR,NPAR_2
      INTEGER, PARAMETER :: ndims=2
      INTEGER :: dims(ndims)
      LOGICAL :: periods(ndims), reorder, remain_dims(ndims)
      INTEGER :: ierror

      IF (NPAR >= COMM%NCPU) NPAR=COMM%NCPU
      dims(1)       = NPAR
      dims(2)       = COMM%NCPU/ NPAR
      IF (dims(1)*dims(2) /= COMM%NCPU ) THEN
         WRITE(0,*) 'M_divide: can not subdivide ',COMM%NCPU,'nodes by',NPAR
      ENDIF

      periods(ndims)=.FALSE.

      CALL MPI_cart_create( COMM%MPI_COMM , ndims, dims, periods, reorder, &
                COMM_CART%MPI_COMM , ierror)
      CHECK_MPI_ERROR(ierror,'M_divide','MPI_cart_create')
! create the in-band communicator
      remain_dims(1)= .FALSE.
      remain_dims(2)= .TRUE.

      CALL MPI_cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INB%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_divide','MPI_cart_sub (1)')

! create the inter-band communicator
      remain_dims(1)= .TRUE.
      remain_dims(2)= .FALSE.

      CALL MPI_cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INTER%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_divide','MPI_cart_sub (2)')
! overwrite initial communicator by new one
      COMM=COMM_CART

      END SUBROUTINE M_divide


!----------------------------------------------------------------------
!
!> subdivides the communicator into two "images" (groups)
!>
!> the first one includes cores 1...NCORE, the second one all the remaining cores
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide2( COMM, NCORE, COMM_INTER, COMM_INB)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM, COMM_INTER, COMM_INB
      INTEGER NCORE
      INTEGER :: color, key
      INTEGER, PARAMETER :: ndims=2
      LOGICAL :: remain_dims(ndims)
      INTEGER :: i,ierror

      IF (NCORE >= COMM%NCPU) NCORE=COMM%NCPU
! first create the communicator COMM_INB which communicates inside one image
      IF (COMM%NODE_ME<=NCORE) THEN
         color=1
         key=COMM%NODE_ME
      ELSE
         color=2
         key=COMM%NODE_ME-NCORE
      ENDIF

      CALL MPI_comm_split(COMM%MPI_COMM, color, key, COMM_INB%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_divide2','MPI_comm_split (1)')

      CALL M_initc( COMM_INB)
!
! create the  communicator COMM_INTER
! this is only defined for the first node in each image
      color=MPI_UNDEFINED
      key=COMM%NODE_ME
      IF (COMM_INB%NODE_ME==1) THEN
         color=1
      ENDIF

      CALL MPI_comm_split(COMM%MPI_COMM, color, key, COMM_INTER%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_divide2','MPI_comm_split (2)')

      IF (COMM_INTER%MPI_COMM /= MPI_COMM_NULL) THEN
         CALL M_initc( COMM_INTER)
      ELSE
         COMM_INTER%NCPU=0
         COMM_INTER%IONODE=1
         ! for VASP internal reasons it is important that
         ! COMM%NODE_ME is equivalent to the "image"
         IF (COMM%NODE_ME<=NCORE) THEN
            COMM_INTER%NODE_ME=1
         ELSE
            COMM_INTER%NODE_ME=2
         ENDIF
      ENDIF

      END SUBROUTINE M_divide2


!----------------------------------------------------------------------
!> splits the communicator COMM into NTAUPAR groups
!>
!> splitting is cartesian and INTER is a proper communicator
!> that allows for collective communication among groups
!----------------------------------------------------------------------

      SUBROUTINE M_divide_general( COMM, NTAUPAR, COMM_INTER, COMM_INTRA, MAP )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      INTEGER NTAUPAR
      INTEGER MAP(NTAUPAR)  !< distribution matrix
      TYPE(communic) COMM, COMM_INTRA, COMM_INTER, COMM_CART
      INTEGER :: NDONE , NREMAIN, NGROUPS
      INTEGER :: ID, I
      INTEGER :: ierr

      !----------------------------------------------------------------------------------------
      !set up node IDs between groups (colums of processor grid)
      !note that the # of tau groups, i.e. COMM_INTER%NCPU is given
      !by NTAUPAR.
      COMM_INTER%NCPU = NTAUPAR
      COMM_INTER%NODE_ME = 0

      IF ( MOD( COMM%NCPU , NTAUPAR ) == 0 ) THEN
         !if NTAUPAR is a divisor of NCPU, there are COMM%NCPU / NTAUPAR cores in each group
         !hence the gropu id is basically NODE_ME/ (COMM%NCPU / NTAUPAR) + 1

         !care must be taken if COMM%NCPU / NTAUPAR  is a divisor of the ID.
         !in this case the group id is simply COMM%NODE_ME / ( COMM%NCPU / NTAUPAR )
         IF ( MOD( COMM%NODE_ME , COMM%NCPU / NTAUPAR ) == 0 ) THEN
            COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR )
         ELSE
            COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR )  + 1
         ENDIF

      ELSE
         IF ( COMM%NODE_ME <=  MOD( COMM%NCPU , NTAUPAR )*(COMM%NCPU / NTAUPAR + 1 )) THEN
            !in the case when NTAUPAR is not a divisor of NCPU the distribution is such that
            !the first MOD( NCPU , NTAUPAR ) groups have NCPU / NTAUPAR + 1 cores
            !so the first MOD( NCPU , NTAUPAR )*( NCPU / NTAUPAR + 1 ) cores have
            !following group IDs
            IF ( MOD( COMM%NODE_ME , COMM%NCPU / NTAUPAR + 1 ) == 0 ) THEN
               COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR + 1 )
            ELSE
               COMM_INTER%NODE_ME = COMM%NODE_ME / ( COMM%NCPU / NTAUPAR + 1 )  + 1
            ENDIF

         ELSE
            !the remaining NCPU - MOD( NCPU , NTAUPAR )*( NCPU / NTAUPAR + 1 ) cores
            !are distributed in groups of NCPU / NTAUPAR cores per group
            !so the groups MOD( NCPU , NTAUPAR ) + 1 ... NTAUPAR contain
            !NCPU / NTAUPAR cores each

            !NDONE CPUs done
            NDONE = MOD( COMM%NCPU , NTAUPAR )*( COMM%NCPU / NTAUPAR + 1 )

            !NREMAIN remaining CPUs
            NREMAIN = COMM%NCPU - NDONE

            !they need to be grouped into NGROUPS
            NGROUPS = NTAUPAR - MOD( COMM%NCPU , NTAUPAR )

            !where each of the NGROUPS groups contain NREMAIN / NGROUPS CPUs
            !check if any cpu is left out
            IF ( MOD( NREMAIN , NGROUPS ) /= 0 ) THEN
               CALL vtutor%bug("M_divide_general: some nodes are missing " // str(NREMAIN) // &
                  " " // str(NGROUPS) // " " // str(NREMAIN / NGROUPS), __FILE__, __LINE__)
            ENDIF

            IF ( .NOT.( (1 <= COMM%NODE_ME - NDONE).AND.(COMM%NODE_ME - NDONE <=  NREMAIN) ) ) THEN
               CALL vtutor%error("M_divide_general: NODE " // str(COMM%NODE_ME) // "has no group id, choose another NTAUPAR")
            ENDIF

            !distribute them according
             IF ( MOD( COMM%NODE_ME - NDONE , NREMAIN / NGROUPS ) == 0 ) THEN
                COMM_INTER%NODE_ME = ( COMM%NODE_ME - NDONE ) / ( NREMAIN / NGROUPS ) + &
                   MOD( COMM%NCPU , NTAUPAR )
             ELSE
                COMM_INTER%NODE_ME = ( COMM%NODE_ME - NDONE ) / ( NREMAIN / NGROUPS ) + 1 +&
                   MOD( COMM%NCPU , NTAUPAR )
             ENDIF
         ENDIF
      ENDIF

      !----------------------------------------------------------------------------------------
      !set up node IDs in groups (rows of processor grid)

      !the distribution matrix MAP(I)
      !gives us the number of nodes in the Ith tau group
      !as determined above there are COMM_INTER%NCPU tau groups
      COMM_INTRA%NCPU = MAP(COMM_INTER%NODE_ME)

      !the IDs inside each tau group can be set elegantly
      !using the distribution matrix
      IF ( COMM_INTER%NODE_ME == 1 ) THEN
         !in the first group the global ids conincide with the intau-ids
         COMM_INTRA%NODE_ME = COMM%NODE_ME
      ELSE
         !in the Ith group the total number of CPU from group 1 to I-1
         !are subtracted from the global node
         COMM_INTRA%NODE_ME = COMM%NODE_ME - SUM( MAP( 1 : COMM_INTER%NODE_ME - 1 ) )
      ENDIF

      !----------------------------------------------------------------------------------------
      !divide the communicator into columns, i.e. each tau group shares a separate communicator
      CALL MPI_comm_split( COMM%MPI_COMM , COMM_INTER%NODE_ME , COMM_INTRA%NODE_ME , &
         COMM_INTRA%MPI_COMM, ierr )
      !check if everything went ok
      IF ( ierr /= 0 ) THEN
         CALL vtutor%bug("M_divide_general: initialization of communicator in tau groups " &
            // "failed " // str(ierr), __FILE__, __LINE__)
      ENDIF

      !----------------------------------------------------------------------------------------
      !divide the communicator into rows, i.e. between tau groups
      !this means all nodes with the same COMM_INTRA%NODE_ME are in the new
      !communicator COMM_BEWWEENTAU%MPI_COMM
      CALL MPI_comm_split( COMM%MPI_COMM , COMM_INTRA%NODE_ME , COMM_INTER%NODE_ME , &
         COMM_INTER%MPI_COMM, ierr )
      !check if everything went ok
      IF ( ierr /= 0 ) THEN
         CALL vtutor%bug("M_divide_general: initialization of communicator between tau groups " &
            // "failed " // str(ierr), __FILE__, __LINE__)
      ENDIF

      !initialize the communicators, just to be sure
      !(although this should have happen already in M_DIVIDE_GENERAL)
      CALLMPI( M_initc( COMM_INTER ) )
      CALLMPI( M_initc( COMM_INTRA ) )

      END SUBROUTINE M_divide_general


!----------------------------------------------------------------------
!> splits communicator into NGROUPS groups
!>
!> if splitting is not cartesian so COMM_INTER is not a proper communicator
!> collective communication among groups using COMM_INTER will crash!
!> cartesian splitting satisfies MOD( COMM%NCPU , NGROUPS ) == 0
!----------------------------------------------------------------------

      SUBROUTINE M_divide_into_groups( COMM, COMM_INTER, COMM_INTRA, NGROUPS )
      IMPLICIT NONE
      TYPE(communic) :: COMM
      TYPE(communic) :: COMM_INTER
      TYPE(communic) :: COMM_INTRA
      INTEGER        :: NGROUPS
      ! local
      INTEGER :: ierror
      INTEGER :: MAP(COMM%NCPU,2)
      INTEGER :: MAP_(COMM%NCPU)
      INTEGER :: LMAP(NGROUPS)
      INTEGER :: ICOLOR, I, NCPU

      COMM_INTER%NCPU = 0
      COMM_INTRA%NCPU = 0

      ! map contains layout how COMM is subdivided into groups
      ! because we allow for non-cartesian splitting,
      ! there is no proper inter-communicator COMM_INTER
      ! which would allow for collective inter-group communications
      ! since one or more groups may contain more ranks than the other

      ! how may ranks are in each group?
      LMAP = 0
      IF ( MOD( COMM%NCPU , NGROUPS ) == 0 ) THEN
         LMAP(1:NGROUPS) = COMM%NCPU/NGROUPS
      ELSE
         LMAP( 1:MOD(COMM%NCPU,NGROUPS) ) = COMM%NCPU/NGROUPS +1
         LMAP( MOD(COMM%NCPU,NGROUPS)+1:NGROUPS ) = COMM%NCPU/NGROUPS
      ENDIF

      ! MAP stores the mapping of global communicator into groups
      MAP = 0

      ! first row of MAP tells us the new group id of each COMM%NODE_ME
      MAP_ = 0
      IF ( LMAP(1)+1 > COMM%NODE_ME ) MAP_( COMM%NODE_ME ) = 1
      DO ICOLOR = 1, NGROUPS-1
         IF ( SUM( LMAP( 1:ICOLOR ) ) < COMM%NODE_ME .AND.&
              SUM( LMAP( 1:ICOLOR+1 ))+1 > COMM%NODE_ME ) THEN
            MAP_( COMM%NODE_ME ) = ICOLOR+1
         ENDIF
      ENDDO
      CALL M_sum_i( COMM, MAP_, SIZE( MAP_ ) )

      ! second row of MAP tells us the new local id of each COMM%NODE_ME
      MAP( COMM%NODE_ME, 2 ) = SUM( LMAP(1:MAP_(COMM%NODE_ME)) )-COMM%NODE_ME+1
      CALL M_sum_i( COMM, MAP, SIZE( MAP ) )
      MAP( :, 1 ) = MAP_(:)

      ! initialize colors
      COMM_INTER%NODE_ME = MAP( COMM%NODE_ME, 1)
      ! initialize node ids in each color
      COMM_INTRA%NODE_ME = MAP( COMM%NODE_ME, 2)

      CALL MPI_comm_split( COMM%MPI_COMM, COMM_INTER%NODE_ME, &
         COMM_INTRA%NODE_ME, COMM_INTRA%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_divide_into_groups','MPI_comm_split (1)')

      CALL M_initc( COMM_INTRA )
      COMM_INTRA%NCPU = LMAP( COMM_INTER%NODE_ME )

      ! in case splitting is cartesian initialize inter communicator as well
      IF ( MOD( COMM%NCPU, NGROUPS ) == 0 ) THEN
         ICOLOR=COMM_INTRA%NODE_ME
      ELSE
         ! this is only defined for the first nodes in each group
         ICOLOR=MPI_UNDEFINED
         NCPU=LMAP( 1 )
         DO I = 1, NGROUPS-1
            NCPU = MIN( NCPU, LMAP( I+1 ) )
         ENDDO
         IF ( COMM_INTRA%NODE_ME <= NCPU ) THEN
            ICOLOR=COMM_INTRA%NODE_ME
         ENDIF
      ENDIF
      CALL MPI_comm_split( COMM%MPI_COMM, ICOLOR, &
         COMM_INTER%NODE_ME, COMM_INTER%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_divide_into_groups','MPI_comm_split (2)')

      IF (COMM_INTER%MPI_COMM /= MPI_COMM_NULL) THEN
         CALL M_initc( COMM_INTER )
      ELSE
         COMM_INTER%NCPU=0
      ENDIF

      END SUBROUTINE M_divide_into_groups


!----------------------------------------------------------------------
!
! M_divide_shmem:
!
!----------------------------------------------------------------------
#ifdef use_shmem
      SUBROUTINE M_divide_shmem(COMM_WORLD,COMM,N,COMM_SHMEM)
      USE c2f_interface
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"
      TYPE (communic) COMM_WORLD,COMM,COMM_SHMEM
      INTEGER N
      ! local variables
      CHARACTER*(MPI_MAX_PROCESSOR_NAME) myname,pname
      INTEGER :: COMM_group,COMM_SHMEM_group
      INTEGER :: group(0:COMM%NCPU-1)
      INTEGER :: resultlen,id_in_group
      INTEGER :: ierror
      INTEGER :: I,IGRP
      ! shmem variables for sanity check
      INTEGER(c_int)    :: shmid
      TYPE(c_ptr)       :: address
      INTEGER(c_size_t) :: k

      INTEGER*4, POINTER :: IDSHMEM(:)
      INTEGER :: ID

      IF (N>0) THEN
      ! explicitly specified division of COMM
         I=COMM%NCPU/N
         IF (I*N/=COMM%NCPU) THEN
            CALL vtutor%bug("M_divide_shmem: can not subdivide " // str(COMM%NCPU) // " nodes " &
               // "by " // str(N), __FILE__, __LINE__)
         ENDIF
         IGRP=0
         DO I=1,COMM%NCPU
            IF ((I-1)/N==(COMM%NODE_ME-1)/N) THEN
               group(IGRP)=I-1
               IGRP=IGRP+1
            ENDIF
         ENDDO
      ENDIF
      IF (N<0) THEN
      ! attempt automatic division of COMM
         IGRP=0
         CALL MPI_Get_processor_name(myname,resultlen,ierror)
         DO I=1,COMM%NCPU
            IF (I==COMM%NODE_ME) pname=myname
            CALL MPI_bcast(pname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,I-1,COMM%MPI_COMM,ierror)
            CHECK_MPI_ERROR(ierror,'M_divide_shmem','MPI_bcast')
#ifdef MPI_bcast_with_barrier
            CALL MPI_barrier(COMM%MPI_COMM,ierror)
#endif
            IF (pname(1:resultlen)==myname(1:resultlen)) THEN
               group(IGRP)=I-1
               IGRP=IGRP+1
            ENDIF
         ENDDO
      ENDIF

      CALL MPI_comm_group(COMM%MPI_COMM,COMM_group,ierror)
      CALL MPI_group_incl(COMM_group,IGRP,group,COMM_SHMEM_group,ierror)
      CALL MPI_comm_create(COMM%MPI_COMM,COMM_SHMEM_group,COMM_SHMEM%MPI_COMM,ierror)

      CALL M_initc(COMM_SHMEM)

      ! sanity check
      k=1
      IF (COMM_SHMEM%NODE_ME==1) CALL getshmem(INT(8*k,KIND=c_size_t),shmid)
      CALL M_bcast_i(COMM_SHMEM,shmid,1)
      CALL attachshmem(shmid,address)
      CALL c_f_pointer(address,IDSHMEM,[k])

      IF (COMM_SHMEM%NODE_ME==1) THEN
         ID=COMM_WORLD%NODE_ME; IDSHMEM=ID
      ENDIF
      CALL M_bcast_i(COMM_SHMEM,ID,1)

      ierror=0
      IF (ID/=IDSHMEM(1)) ierror=1
      CALL M_sum_i(COMM_WORLD,ierror,1)

      CALL detachshmem(address)
      IF (COMM_SHMEM%NODE_ME==1) CALL destroyshmem(shmid)

      IF (ierror>0) THEN
         CALL vtutor%bug("M_divide_shmem: not all procs in COMM_SHMEM seem to be on the same " &
            // "physical node.", __FILE__, __LINE__)
      ENDIF

      N=COMM_SHMEM%NCPU

!     DO I=1,COMM_WORLD%NCPU
!        IF (COMM_WORLD%NODE_ME==I) THEN
!           WRITE(*,'(A,I4,X,A,I4)') 'global id:',comm_world%node_me,'shmem id:',comm_shmem%node_me
!        ENDIF
!        CALL MPI_barrier(COMM_WORLD%MPI_COMM,ierror)
!     ENDDO

      END SUBROUTINE M_divide_shmem
#endif

!----------------------------------------------------------------------
!
! M_divide_intra_inter_node:
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide_intra_inter_node(COMM_WORLD,COMM,COMM_intra,COMM_inter)
      USE c2f_interface
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"
      TYPE (communic) COMM_WORLD,COMM,COMM_intra,COMM_inter,COMM_test
      ! local variables
      CHARACTER*(MPI_MAX_PROCESSOR_NAME) myname,pname
      INTEGER :: COMM_group,COMM_intra_group,COMM_inter_group
      INTEGER :: group(0:COMM%NCPU-1)
      INTEGER :: resultlen,myid
      INTEGER :: ierror
      INTEGER :: I,IGRP
      ! shmem variables for sanity check
      INTEGER(c_int)    :: shmid
      TYPE(c_ptr)       :: address
      INTEGER(c_size_t) :: k

      INTEGER*4, POINTER :: IDSHMEM(:)
      INTEGER :: ID

      ! attempt automatic division of COMM
#define use_comm_split_type
#ifndef use_comm_split_type
      IGRP=0
      CALL MPI_Get_processor_name(myname,resultlen,ierror)
      DO I=1,COMM%NCPU
         IF (I==COMM%NODE_ME) pname=myname
         CALL MPI_bcast(pname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,I-1,COMM%MPI_COMM,ierror)
         CHECK_MPI_ERROR(ierror,'M_divide_intra_inter_node','MPI_bcast (1)')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier(COMM%MPI_COMM,ierror)
#endif
         IF (pname(1:resultlen)==myname(1:resultlen)) THEN
            group(IGRP)=I-1
            IGRP=IGRP+1
         ENDIF
      ENDDO

      CALL MPI_comm_group(COMM%MPI_COMM,COMM_group,ierror)
      CALL MPI_group_incl(COMM_group,IGRP,group,COMM_intra_group,ierror)
      CALL MPI_comm_create(COMM%MPI_COMM,COMM_intra_group,COMM_intra%MPI_COMM,ierror)
#else
      CALL MPI_comm_split_type(COMM%MPI_COMM,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,COMM_intra%MPI_COMM,ierror)
      CALL MPI_comm_group(COMM_intra%MPI_COMM,COMM_intra_group,ierror)
      CALL MPI_comm_group(COMM%MPI_COMM,COMM_group,ierror)
#endif
      CALL M_initc(COMM_intra)

#ifdef use_shmem
      ! sanity check
      k=1
      IF (COMM_intra%NODE_ME==1) CALL getshmem(INT(8*k,KIND=c_size_t),shmid)
      CALL M_bcast_i(COMM_intra,shmid,1)
      CALL attachshmem(shmid,address)
      CALL c_f_pointer(address,IDSHMEM,[k])

      IF (COMM_intra%NODE_ME==1) THEN
         ID=COMM_WORLD%NODE_ME; IDSHMEM=ID
      ENDIF
      CALL M_bcast_i(COMM_intra,ID,1)

      ierror=0
      IF (ID/=IDSHMEM(1)) ierror=1
      CALL M_sum_i(COMM_WORLD,ierror,1)

      CALL detachshmem(address)
      IF (COMM_intra%NODE_ME==1) CALL destroyshmem(shmid)

      IF (ierror>0) THEN
         CALL vtutor%bug("M_divide_intra_inter_node: not all procs in COMM_intra seem to be on the " &
            // "same physical node.", __FILE__, __LINE__)
      ENDIF
#endif

      IGRP=0
      CALL MPI_group_rank(COMM_intra_group,myid,ierror)
      DO I=1,COMM%NCPU
         IF (I==COMM%NODE_ME) ID=myid
         CALL MPI_bcast(ID,1,MPI_INTEGER,I-1,COMM%MPI_COMM,ierror)
         CHECK_MPI_ERROR(ierror,'M_divide_intra_inter_node','MPI_bcast (2)')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier(COMM%MPI_COMM,ierror)
#endif
         IF (ID==myid) THEN
            group(IGRP)=I-1
            IGRP=IGRP+1
         ENDIF
      ENDDO

      CALL MPI_group_incl(COMM_group,IGRP,group,COMM_inter_group,ierror)
      CALL MPI_comm_create(COMM%MPI_COMM,COMM_inter_group,COMM_inter%MPI_COMM,ierror)

      CALL M_initc(COMM_inter)

!      DO I=1,COMM_WORLD%NCPU
!         IF (COMM_WORLD%NODE_ME==I) THEN
!            WRITE(*,'(A,I4,X,A,I4,X,A,I4)') 'global id:',comm_world%node_me,'intra node id:',comm_intra%node_me,'inter node id:',comm_inter%node_me
!         ENDIF
!         CALL MPI_barrier(COMM_WORLD%MPI_COMM,ierror)
!      ENDDO

      END SUBROUTINE M_divide_intra_inter_node


!----------------------------------------------------------------------
!
!> initialise a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_initc( COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror
      INTEGER id_in_group

      CALL MPI_comm_rank( COMM%MPI_COMM, id_in_group, ierror )
      CHECK_MPI_ERROR(ierror,'M_initc','MPI_comm_rank')

      CALL MPI_comm_size( COMM%MPI_COMM, COMM%NCPU, ierror )
      CHECK_MPI_ERROR(ierror,'M_initc','MPI_comm_size')

      COMM%NODE_ME = id_in_group + 1

      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_initc','MPI_barrier')
      COMM%IONODE =1

      END SUBROUTINE M_initc


!----------------------------------------------------------------------
!
!> frees a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_freec( COMM )
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror
      INTEGER id_in_group

      CALL MPI_comm_free( COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_freec','MPI_comm_free')

      COMM%NODE_ME = 0
      COMM%IONODE = 0
      COMM%NCPU = 0

      END SUBROUTINE M_freec

      END MODULE mpimy


!======================================================================
!
! all other routines are often called with either
! real or complex arrays, vectors or scalars, so I can not put
! them into the F90 module
!
!======================================================================
!----------------------------------------------------------------------
!
!> exit MPI, very simple just exit MPI
!
!----------------------------------------------------------------------

      SUBROUTINE M_exit()
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      INTEGER ierror

!      CALL MPI_barrier(MPI_comm_world, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('M_exit: MPI_barrier returns: ',ierror)

      ! freq quadruple summation operator
      CALL MPI_op_free( M_sum_qd_op, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('M_exit: MPI_op_free returns: ',ierror)

      CALL MPI_finalize( ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('M_exit: MPI_finalize returns: ',ierror)
      STOP

      END SUBROUTINE M_exit


!----------------------------------------------------------------------
!
!> exits MPI and program because of error, a message is printed
!
!----------------------------------------------------------------------

      SUBROUTINE M_stop(message)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER (LEN=*) message
      INTEGER request, status, ierror
      LOGICAL flag
      INTEGER i

      CALL MPI_ibarrier(MPI_comm_world, request, ierror )

      WRITE (*,*) message

      DO i=1,30
         CALL MPI_Test(request, flag, status, ierror )
         IF (flag) CALL M_exit()
         CALL SLEEP(1)
      ENDDO

      CALL MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      END SUBROUTINE M_stop


      SUBROUTINE M_stop_ierr(message, ierror)
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER (LEN=*) message
      INTEGER ierror

      WRITE (*,*) message, ierror

      CALL MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      END SUBROUTINE M_stop_ierr


!======================================================================
!
! Send and Receive routines, map directly onto MPI
!
!======================================================================

!----------------------------------------------------------------------
!
!> send n integers stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_i (COMM, node, ivec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)

      INTEGER status(MPI_status_size), ierror

      CALL MPI_send( ivec(1), n, MPI_integer, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_send_i','MPI_send')

      END SUBROUTINE M_send_i


!----------------------------------------------------------------------
!
!> send n 8-byte integers stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_i8 (COMM, node, ivec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER(qi8) ivec(n)

      INTEGER status(MPI_status_size), ierror

      CALL MPI_send( ivec(1), n, MPI_integer8, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_send_i8','MPI_send')

      END SUBROUTINE M_send_i8


!----------------------------------------------------------------------
!
!> receive n integers into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_i(COMM, node, ivec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)
      INTEGER status(MPI_status_size), ierror

      CALL MPI_recv( ivec(1), n, MPI_integer , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CHECK_MPI_ERROR(ierror,'M_recv_i','MPI_recv')

      END SUBROUTINE M_recv_i


!----------------------------------------------------------------------
!
!> receive n 8-byte integers into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_i8(COMM, node, ivec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER(qi8) ivec(n)
      INTEGER status(MPI_status_size), ierror

      CALL MPI_recv( ivec(1), n, MPI_integer8 , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CHECK_MPI_ERROR(ierror,'M_recv_i8','MPI_recv')

      END SUBROUTINE M_recv_i8


!----------------------------------------------------------------------
!
!> send n double complex stored in zvec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_z (COMM, node, zvec, n)

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      COMPLEX(q) :: zvec(n)

      INTEGER status(MPI_status_size), ierror

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(zvec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_send_d(zvec(1), n*2, node-1, COMM, ACC_ASYNC_Q)
            RETURN
         ENDIF
#endif // USENCCL
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(zvec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(zvec)
         CALL MPI_send( zvec(1), n, MPI_double_complex, node-1, 200, &
        &               COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_send_z','MPI_send')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_send( zvec(1), n, MPI_double_complex, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_send_z','MPI_send')

      END SUBROUTINE M_send_z


!----------------------------------------------------------------------
!
!> receive n double complex into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_z(COMM, node, zvec, n )

      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      COMPLEX(q) :: zvec(n)
      INTEGER status(MPI_status_size), ierror

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(zvec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_recv_d(zvec(1), n*2, node-1, COMM, ACC_ASYNC_Q)
            RETURN
         ENDIF
#endif // USENCCL
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(zvec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(zvec)
         CALL MPI_recv( zvec(1), n, MPI_double_complex , node-1, 200, &
        &               COMM%MPI_COMM, status, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_recv_z','MPI_recv')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_recv( zvec(1), n, MPI_double_complex , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CHECK_MPI_ERROR(ierror,'M_recv_z','MPI_recv')

      END SUBROUTINE M_recv_z


!----------------------------------------------------------------------
!
!> send n double stored in ivec to node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_d (COMM, node, dvec, n)

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      REAL(q) :: dvec(n)

      INTEGER status(MPI_status_size), ierror

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(dvec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_send_d(dvec(1), n, node-1, COMM, ACC_ASYNC_Q)
            RETURN
         ENDIF
#endif // USENCCL
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(dvec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(dvec)
         CALL MPI_send( dvec(1), n, MPI_double_precision, node-1, 200, &
        &               COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_send_d','MPI_send')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_send( dvec(1), n, MPI_double_precision, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_send_d','MPI_send')

      END SUBROUTINE M_send_d


!----------------------------------------------------------------------
!
!> receive n double  into array ivec from node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_d(COMM, node, dvec, n )

      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      REAL(q) :: dvec(n)
      INTEGER status(MPI_status_size), ierror

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(dvec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_recv_d(dvec(1), n, node-1, COMM, ACC_ASYNC_Q)
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(dvec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(dvec)
         CALL MPI_recv( dvec(1), n, MPI_double_precision , node-1, 200, &
        &               COMM%MPI_COMM, status, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_recv_d','MPI_recv')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_recv( dvec(1), n, MPI_double_precision , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      CHECK_MPI_ERROR(ierror,'M_recv_d','MPI_recv')

      END SUBROUTINE M_recv_d


!======================================================================
!
! global sum and maximum routines
!
!======================================================================

#ifdef MPI_INPLACE
! employ the IN_PLACE versions of MPI_allreduce and do not split

!----------------------------------------------------------------------
!
!> Check if the logicals are true on all cores
!
!----------------------------------------------------------------------

      SUBROUTINE M_and(COMM, lvec, n)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(INOUT) :: lvec(n)

      CALL M_op(COMM, MPI_LAND, 'M_and', lvec, n)

      END SUBROUTINE M_and

!----------------------------------------------------------------------
!
!> Check if any of the logicals are true on all cores
!
!----------------------------------------------------------------------

      SUBROUTINE M_or(COMM, lvec, n)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(INOUT) :: lvec(n)

      CALL M_op(COMM, MPI_LOR, 'M_or', lvec, n)

      END SUBROUTINE M_or


!----------------------------------------------------------------------
!
!> Reduce logicals across MPI ranks using some MPI operation
!
!----------------------------------------------------------------------
      SUBROUTINE M_op(COMM, operation, op_str, lvec, n)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: operation
      CHARACTER(LEN=*), INTENT(IN) :: op_str
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(INOUT) :: lvec(n)

      INTEGER ierror

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug(trim(op_str)//": invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce(MPI_IN_PLACE, lvec(1), n, MPI_LOGICAL, operation, COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,op_str,'MPI_allreduce')

      END SUBROUTINE M_op

!----------------------------------------------------------------------
!
!> Performs an in-place global sum on n integers in vector ivec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i(COMM, ivec, n )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: ivec(n)

      INTEGER :: ierror, status(MPI_status_size), ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_INTEGER, &
         &                MPI_SUM, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_sum_i','MPI_allreduce')

      END SUBROUTINE M_sum_i

!----------------------------------------------------------------------
!
!> Performs an in-place global sum on n 8-byte integers in vector ivec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i8(COMM, ivec, n )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      INTEGER(qi8), INTENT(INOUT) :: ivec(n)

      INTEGER :: ierror, status(MPI_status_size), ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_INTEGER8, &
         &                MPI_SUM, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_sum_i8','MPI_allreduce')

      END SUBROUTINE M_sum_i8


!----------------------------------------------------------------------
!
!> Performs an in-place global min on n integers in vector ivec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_i(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: ivec(n)

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_min_i')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_min_i(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_min_i')
            RETURN
         ENDIF
#endif // USENCCL

! invoke CUDA-aware in-place MPI_allreduce
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec)
         CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_integer, &
        &                    MPI_min, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_min_i','MPI_allreduce')
         PROFILING_STOP('m_min_i')
         RETURN
      ENDIF
#endif // OFFLOADING

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_INTEGER, &
         &                MPI_MIN, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_min_i','MPI_allreduce')

      PROFILING_STOP('m_min_i')

      END SUBROUTINE M_min_i


!----------------------------------------------------------------------
!
!> Performs an in-place global min on n 8-byte integers in vector ivec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_i8(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      INTEGER(qi8), INTENT(INOUT) :: ivec(n)

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_min_i8')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_min_i8(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_min_i8')
            RETURN
         ENDIF
#endif // USENCCL

! invoke CUDA-aware in-place MPI_allreduce
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec)
         CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_integer8, &
        &                    MPI_min, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_min_i8','MPI_allreduce')
         PROFILING_STOP('m_min_i8')
         RETURN
      ENDIF
#endif // OFFLOADING

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_INTEGER, &
         &                MPI_MIN, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_min_i8','MPI_allreduce')

      PROFILING_STOP('m_min_i8')

      END SUBROUTINE M_min_i8


!----------------------------------------------------------------------
!
!> Performs an in-place global maximum on n integers in vector ivec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: ivec(n)

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_max_i')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_max_i(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_max_i')
            RETURN
         ENDIF
#endif // USENCCL

! invoke CUDA-aware in-place MPI_allreduce
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec)
         CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_integer, &
        &                    MPI_max, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_max_i','MPI_allreduce')
         PROFILING_STOP('m_max_i')
         RETURN
      ENDIF
#endif // OFFLOADING

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_INTEGER, &
         &                MPI_MAX, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_max_i','MPI_allreduce')

      PROFILING_STOP('m_max_i')

      END SUBROUTINE M_max_i


!----------------------------------------------------------------------
!
!> Performs an in-place global maximum on n 8-byte integers in vector ivec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i8(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      INTEGER(qi8), INTENT(INOUT) :: ivec(n)

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_max_i8')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_max_i8(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_max_i8')
            RETURN
         ENDIF
#endif // USENCCL

! invoke CUDA-aware in-place MPI_allreduce
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec)
         CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_integer8, &
        &                    MPI_max, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_max_i8','MPI_allreduce')
         PROFILING_STOP('m_max_i8')
         RETURN
      ENDIF
#endif // OFFLOADING

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, ivec(1), n, MPI_INTEGER8, &
         &                MPI_MAX, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_max_i8','MPI_allreduce')

      PROFILING_STOP('m_max_i8')

      END SUBROUTINE M_max_i8


!----------------------------------------------------------------------
!
!> Performs an in-place global maximum on n doubles in vector vec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      REAL(q), INTENT(INOUT) :: vec(n)

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_DOUBLE_PRECISION, &
         &                MPI_MAX, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_max_d','MPI_allreduce')

      END SUBROUTINE M_max_d


!----------------------------------------------------------------------
!
!> Performs an in-place global minimum on n doubles in vector vec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      REAL(q), INTENT(INOUT) :: vec(n)

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_DOUBLE_PRECISION, &
         &                MPI_MIN, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_min_d','MPI_allreduce')

      END SUBROUTINE M_min_d


!----------------------------------------------------------------------
!
!> Performs an in-place global sum on n doubles in vector vec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      REAL(q), INTENT(INOUT) :: vec(n)

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_sumb_d')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_sum_d(vec(1), vec(1), n, COMM, ACC_ASYNC_Q)
!note(sm): this wait basically defeats the purpose of using NCCL asynchronously
GPACC !$ACC WAIT(ACC_ASYNC_Q)
            PROFILING_STOP('m_sumb_d')
            RETURN
         ENDIF
#endif // USENCCL

! invoke CUDA-aware in-place MPI_allreduce
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)

         CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_DOUBLE_PRECISION, &
        &                    MPI_SUM, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_sumb_d','MPI_allreduce')
         PROFILING_STOP('m_sumb_d')
         RETURN
      ENDIF
#endif // OFFLOADING

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_DOUBLE_PRECISION, &
         &                MPI_SUM, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_sumb_d','MPI_allreduce')

      PROFILING_STOP('m_sumb_d')

      END SUBROUTINE M_sumb_d


!----------------------------------------------------------------------
!
!> Performs an in-place global maximum on n singles in vector vec,
!> invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_s(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      REAL(qs), INTENT(INOUT) :: vec(n)

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_s: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_sumb_s')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_sum_s(vec(1), vec(1), n, COMM, ACC_ASYNC_Q)
!note(sm): this wait basically defeats the purpose of using NCCL asynchronously
GPACC !$ACC WAIT(ACC_ASYNC_Q)
            PROFILING_STOP('m_sumb_s')
            RETURN
         ENDIF
#endif // USENCCL

! invoke CUDA-aware in-place MPI_allreduce
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_REAL, &
        &                    MPI_SUM, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_sumb_s','MPI_allreduce')
         PROFILING_STOP('m_sumb_s')
         RETURN
      ENDIF
#endif // OFFLOADING

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_REAL, &
         &                MPI_SUM, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_sumb_s','MPI_allreduce')

      PROFILING_STOP('m_sumb_s')

      END SUBROUTINE M_sumb_s


!----------------------------------------------------------------------
!
!> Performs an in-place global sum on n double precision complex items
!> in the vector vec, invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_z(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      COMPLEX(q), INTENT(INOUT) :: vec(n)

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_z: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_DOUBLE_COMPLEX, &
         &                MPI_SUM, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_sumb_z','MPI_allreduce')

      END SUBROUTINE M_sumb_z


!----------------------------------------------------------------------
!
!> Performs an in-place global product on n double precision complex
!> items in the vector vec, invoking MPI_allreduce.
!
!----------------------------------------------------------------------

      SUBROUTINE M_prodb_z(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      COMPLEX(q), INTENT(INOUT) :: vec(n)

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_prodb_z: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_DOUBLE_COMPLEX, &
         &                MPI_PROD, COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_prodb_z','MPI_allreduce')

      END SUBROUTINE M_prodb_z


!----------------------------------------------------------------------
!
!> Performs an in-place global sum on n doubles in vector vec
!>
!> The array is split into smaller pieces and MPI_allreduce is invoked
!> without the MPI_INPLACE flag.
!> This routine is for testing purposes against the non splitting
!> version.
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_d_splitting(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_d_splitting: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sumb_d_splitting','MPI_allreduce')

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)

      ENDDO

      END SUBROUTINE M_sumb_d_splitting

#else // MPI_INPLACE

! split the arrays and copy for MPI_allreduce

!----------------------------------------------------------------------
!
!> performs a global and on n logicals in vector lvec
!
!----------------------------------------------------------------------

      SUBROUTINE M_and(COMM, lvec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      LOGICAL lvec(n)

      INTEGER ierror, ichunk

      CALL M_op(COMM, MPI_LAND, 'M_and', lvec, n)

      END SUBROUTINE M_and

!----------------------------------------------------------------------
!
!> performs a global or on n logicals in vector lvec
!
!----------------------------------------------------------------------

      SUBROUTINE M_or(COMM, lvec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      LOGICAL lvec(n)

      INTEGER ierror, ichunk

      CALL M_op(COMM, MPI_LOR, 'M_or', lvec, n)

      END SUBROUTINE M_or

!----------------------------------------------------------------------
!
!> Reduce logicals across MPI ranks using some MPI operation
!
!----------------------------------------------------------------------
      SUBROUTINE M_op(COMM, operation, op_str, lvec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: operation
      CHARACTER(LEN=*), INTENT(IN) :: op_str
      INTEGER,INTENT(IN) ::  n
      LOGICAL,INTENT(INOUT) :: lvec(n)

      INTEGER j,k
      INTEGER ierror, ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug(trim(op_str)//": invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NLTMP
         ichunk = MIN( n-j+1 , NLTMP)

         CALL MPI_allreduce(lvec(j), LTMP_m(1), ichunk, MPI_logical, &
     &                      operation, COMM%MPI_COMM, ierror)
         CHECK_MPI_ERROR(ierror, op_str, 'MPI_allreduce')

         DO k = 0, ichunk-1
            lvec(j+k) = LTMP_m(k+1)
         ENDDO

      ENDDO

      END SUBROUTINE M_op

!----------------------------------------------------------------------
!
!> performs a global sum on n integers in vector ivec
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
     &                       MPI_sum, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sum_i','MPI_allreduce')

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO

      ENDDO

      END SUBROUTINE M_sum_i


!----------------------------------------------------------------------
!
!> performs a global sum on n integers in vector ivec
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i8(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER(qi8) ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NI8TMP
         ichunk = MIN( n-j+1 , NI8TMP)

         CALL MPI_allreduce( ivec(j), I8TMP_m(1), ichunk, MPI_integer8, &
     &                       MPI_sum, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sum_i8','MPI_allreduce')

         DO k = 0, ichunk-1
            ivec(j+k) = I8TMP_m(k+1)
         ENDDO

      ENDDO

      END SUBROUTINE M_sum_i8


!----------------------------------------------------------------------
!
!> performs a global min on n integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_i(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_min_i')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_min_i(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_min_i')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC ENTER DATA CREATE(ITMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:ITMP_m)
         DO j = 1, n, NITMP
            ichunk = MIN( n-j+1 , NITMP)

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec,ITMP_m)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec,ITMP_m)
            CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                                MPI_min, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CHECK_MPI_ERROR(ierror,'M_min_i','MPI_allreduce')

GPACC !$ACC PARALLEL LOOP PRESENT(ivec,ITMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD PRIVATE(K) IF(OFFLOAD_ON)
            DO k = 0, ichunk-1
               ivec(j+k) = ITMP_m(k+1)
            ENDDO
         ENDDO
GPACC !$ACC EXIT DATA DELETE(ITMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:ITMP_m)
         PROFILING_STOP('m_min_i')
         RETURN
      ENDIF
#endif // OFFLOADING

!  there is no inplace global min in MPI, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                             MPI_min, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_min_i','MPI_allreduce')

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO
      ENDDO

      PROFILING_STOP('m_min_i')

      END SUBROUTINE M_min_i


!----------------------------------------------------------------------
!
!> performs a global min on n 8-byte integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_i8(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER(qi8) ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_min_i8')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_min_i8(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_min_i8')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC ENTER DATA CREATE(I8TMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:I8TMP_m)
         DO j = 1, n, NI8TMP
            ichunk = MIN( n-j+1 , NI8TMP)

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec,I8TMP_m)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec,I8TMP_m)
            CALL MPI_allreduce( ivec(j), I8TMP_m(1), ichunk, MPI_integer8, &
                                MPI_min, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CHECK_MPI_ERROR(ierror,'M_min_i8','MPI_allreduce')

GPACC !$ACC PARALLEL LOOP PRESENT(ivec,I8TMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD PRIVATE(K) IF(OFFLOAD_ON)
            DO k = 0, ichunk-1
               ivec(j+k) = I8TMP_m(k+1)
            ENDDO
         ENDDO
GPACC !$ACC EXIT DATA DELETE(I8TMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:I8TMP_m)
         PROFILING_STOP('m_min_i8')
         RETURN
      ENDIF
#endif // OFFLOADING

!  there is no inplace global min in MPI, thus we have to use
!  a work array

      DO j = 1, n, NI8TMP
         ichunk = MIN( n-j+1 , NI8TMP)

         CALL MPI_allreduce( ivec(j), I8TMP_m(1), ichunk, MPI_integer8, &
                             MPI_min, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_min_i8','MPI_allreduce')

         DO k = 0, ichunk-1
            ivec(j+k) = I8TMP_m(k+1)
         ENDDO
      ENDDO

      PROFILING_STOP('m_min_i8')

      END SUBROUTINE M_min_i8


!----------------------------------------------------------------------
!
! M_max_i: performs a global max on n integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_max_i')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_max_i(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_max_i')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC ENTER DATA CREATE(ITMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:ITMP_m)
         DO j = 1, n, NITMP
            ichunk = MIN( n-j+1 , NITMP)

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec,ITMP_m)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec,ITMP_m)
            CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
           &                    MPI_max, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CHECK_MPI_ERROR(ierror,'M_max_i','MPI_allreduce')

GPACC !$ACC PARALLEL LOOP PRESENT(ivec,ITMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD PRIVATE(K) IF(OFFLOAD_ON)
            DO k = 0, ichunk-1
               ivec(j+k) = ITMP_m(k+1)
            ENDDO
         ENDDO
GPACC !$ACC EXIT DATA DELETE(ITMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:ITMP_m)
         PROFILING_STOP('m_max_i')
         RETURN
      ENDIF
#endif // OFFLOADING

!  there is no inplace global max in MPI, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         CALL MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                             MPI_max, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_max_i','MPI_allreduce')

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO
      ENDDO

      PROFILING_STOP('m_max_i')

      END SUBROUTINE M_max_i


!----------------------------------------------------------------------
!
! M_max_i8: performs a global max on n 8-byte integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i8(COMM, ivec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER(qi8) ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_max_i8')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(ivec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_max_i8(ivec(1), ivec(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_max_i8')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC ENTER DATA CREATE(I8TMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:I8TMP_m)
         DO j = 1, n, NI8TMP
            ichunk = MIN( n-j+1 , NI8TMP)

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(ivec,I8TMP_m)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(ivec,I8TMP_m)
            CALL MPI_allreduce( ivec(j), I8TMP_m(1), ichunk, MPI_integer, &
           &                    MPI_max, COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CHECK_MPI_ERROR(ierror,'M_max_i8','MPI_allreduce')

GPACC !$ACC PARALLEL LOOP PRESENT(ivec,I8TMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD PRIVATE(K) IF(OFFLOAD_ON)
            DO k = 0, ichunk-1
               ivec(j+k) = I8TMP_m(k+1)
            ENDDO
         ENDDO
GPACC !$ACC EXIT DATA DELETE(I8TMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:I8TMP_m)
         PROFILING_STOP('m_max_i8')
         RETURN
      ENDIF
#endif // OFFLOADING

!  there is no inplace global max in MPI, thus we have to use
!  a work array

      DO j = 1, n, NI8TMP
         ichunk = MIN( n-j+1 , NI8TMP)

         CALL MPI_allreduce( ivec(j), I8TMP_m(1), ichunk, MPI_integer8, &
                             MPI_max, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_max_i8','MPI_allreduce')

         DO k = 0, ichunk-1
            ivec(j+k) = I8TMP_m(k+1)
         ENDDO
      ENDDO

      PROFILING_STOP('m_max_i8')

      END SUBROUTINE M_max_i8


!----------------------------------------------------------------------
!
!> performs a global max search on n doubles in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_max_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global max in MPI, thus we have to use
!  a work array
      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_max, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_max_d','MPI_allreduce')

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_max_d


!----------------------------------------------------------------------
!
!> performs a global min search on n doubles in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_min_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_min_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global max in MPI, thus we have to use
!  a work array
      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_min, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_min_d','MPI_allreduce')

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_min_d


!----------------------------------------------------------------------
!
!> performs a global sum on n doubles in vector vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!>  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_sumb_d')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
GPACC !$ACC ENTER DATA CREATE(DTMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:DTMP_m)
         DO j = 1, n, NDTMP
            ichunk = MIN( n-j+1 , NDTMP)
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec,DTMP_m)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec,DTMP_m)
            CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                                MPI_double_precision, MPI_sum, &
                                COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CALL WDCOPY(ichunk, DTMP_m(1), 1,  vec(j), 1)
         ENDDO
GPACC !$ACC EXIT DATA DELETE(DTMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:DTMP_m)
         PROFILING_STOP('m_sumb_d')
         RETURN
      ENDIF
#endif // OFFLOADING

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sumb_d','MPI_allreduce')

         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      PROFILING_STOP('m_sumb_d')

      END SUBROUTINE M_sumb_d


!----------------------------------------------------------------------
!
!> performs a global sum on n singles in vector vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!>  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_s(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_s: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_sumb_s')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
GPACC !$ACC ENTER DATA CREATE(DTMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:DTMP_m)
         DO j = 1, n, NDTMP
            ichunk = MIN( n-j+1 , NDTMP)
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec,DTMP_m)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec,DTMP_m)
            CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                                MPI_REAL, MPI_sum, &
                                COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CALL WSCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
         ENDDO
GPACC !$ACC EXIT DATA DELETE(DTMP_m) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(DELETE:DTMP_m)
         PROFILING_STOP('m_sumb_s')
         RETURN
      ENDIF
#endif // OFFLOADING

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_REAL, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sumb_s','MPI_allreduce')

         CALL SCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      PROFILING_STOP('m_sumb_s')

      END SUBROUTINE M_sumb_s


!----------------------------------------------------------------------
!
!> performs a global sum on n complex items in vector vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!>  faster alternative routines can be found below
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_z(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sumb_z: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array
      DO j = 1, n, NZTMP
         ichunk = MIN( n-j+1 , NZTMP)

         CALL MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sumb_z','MPI_allreduce')

         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_sumb_z


!----------------------------------------------------------------------
!
!> performs a global product on n complex numbers in vec
!>
!>  uses MPI_allreduce which is usually very inefficient
!
!----------------------------------------------------------------------

      SUBROUTINE M_prodb_z(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_prodb_z: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         CALL MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_prod, &
                             COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_prodb_z','MPI_allreduce')

         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)
      ENDDO

      END SUBROUTINE M_prodb_z

#endif // MPI_INPLACE

!----------------------------------------------------------------------
!
!> to make live easier, a global sum for 2 scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_2(COMM, v1, v2)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      REAL(q) vec(2),v1,v2

      vec(1)=v1
      vec(2)=v2

      CALL M_sumb_d(COMM, vec, 2)

      v1=vec(1)
      v2=vec(2)

      END SUBROUTINE M_sum_2


!----------------------------------------------------------------------
!
!> a global sum for 3 scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_3(COMM, v1, v2, v3)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      REAL(q) vec(3),v1,v2,v3

      vec(1)=v1
      vec(2)=v2
      vec(3)=v3

      CALL M_sumb_d(COMM, vec, 3)

      v1=vec(1)
      v2=vec(2)
      v3=vec(3)

      END SUBROUTINE M_sum_3


!----------------------------------------------------------------------
!
!> a global sum for 4 scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_4(COMM, v1, v2, v3, v4)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      REAL(q) vec(4),v1,v2,v3,v4

      vec(1)=v1
      vec(2)=v2
      vec(3)=v3
      vec(4)=v4

      CALL M_sumb_d(COMM, vec, 4)

      v1=vec(1)
      v2=vec(2)
      v3=vec(3)
      v4=vec(4)

      END SUBROUTINE M_sum_4


!======================================================================
!
!> Global barrier routine
!
!======================================================================

      SUBROUTINE M_barrier(COMM )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_barrier','MPI_barrier')

      END SUBROUTINE M_barrier

!======================================================================
!
! Global Copy Routines
!
!======================================================================

!----------------------------------------------------------------------
!
!> copy n characters of KIND=c_char from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_c_char(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE iso_c_binding, ONLY : c_char
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      CHARACTER(KIND=c_char) :: vec(n)
      INTEGER n

      ! local
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_c_char: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      CALL MPI_Bcast(vec(1), n, MPI_CHAR, COMM%IONODE-1, COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_bcast_c_char','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
      END SUBROUTINE M_bcast_c_char


!----------------------------------------------------------------------
!
!> copy n integers from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_i: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF
#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n, MPI_integer, COMM%IONODE-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else
      CALL MPI_bcast( vec(1), n, MPI_integer, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif
      END SUBROUTINE M_bcast_i


!----------------------------------------------------------------------
!
!> copy n 8-byte integers from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i8(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER(qi8) vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_i8: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF
#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n, MPI_integer8, COMM%IONODE-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i8','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else
      CALL MPI_bcast( vec(1), n, MPI_integer8, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i8','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif
      END SUBROUTINE M_bcast_i8


!----------------------------------------------------------------------
!
!> copy n logical from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_l(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      LOGICAL vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_l: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF
#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n, MPI_logical, COMM%IONODE-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_l','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else
      CALL MPI_bcast( vec(1), n, MPI_logical, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_l','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif
      END SUBROUTINE M_bcast_l


!----------------------------------------------------------------------
!
!> copy n integers from node inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i_from(COMM, vec, n , inode)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER inode
      INTEGER vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_i_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF
#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n, MPI_integer, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else
      CALL MPI_bcast( vec(1), n, MPI_integer, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif
      END SUBROUTINE M_bcast_i_from


!----------------------------------------------------------------------
!
!> copy n 8-byte integers from node inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i8_from(COMM, vec, n , inode)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER inode
      INTEGER(qi8) vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_i8_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF
#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n, MPI_integer8, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i8_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else
      CALL MPI_bcast( vec(1), n, MPI_integer8, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_i8_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif
      END SUBROUTINE M_bcast_i8_from


!----------------------------------------------------------------------
!
!> copy n double precision from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_d: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_bcast_d')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_d(vec(1), n, COMM%IONODE-1, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_bcast_d')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_double_precision, COMM%IONODE-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_bcast_d','MPI_bcast')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         PROFILING_STOP('m_bcast_d')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_double_precision, COMM%IONODE-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_d','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_double_precision, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_d','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast

      PROFILING_STOP('m_bcast_d')

      END SUBROUTINE M_bcast_d


!----------------------------------------------------------------------
!
!> copy n single precision from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_s(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_s: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_bcast_s')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_real, COMM%IONODE-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_bcast_s','MPI_bcast')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         PROFILING_STOP('m_bcast_s')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_real, COMM%IONODE-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_s','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_real, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_s','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast
      PROFILING_STOP('m_bcast_s')

      END SUBROUTINE M_bcast_s


!----------------------------------------------------------------------
!
!> copy n double precision complex from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z(COMM, vec, n )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_bcast_z')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_d(vec(1), n*2, COMM%IONODE-1, COMM, ACC_ASYNC_Q)
!note(sm): this wait basically defeats the purpose of using NCCL asynchronously
GPACC !$ACC WAIT(ACC_ASYNC_Q)
            PROFILING_STOP('m_bcast_z')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_double_complex, COMM%IONODE-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_bcast_z','MPI_bcast')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         PROFILING_STOP('m_bcast_z')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_double_complex, COMM%IONODE-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_z','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_double_complex, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_z','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast
      PROFILING_STOP('m_bcast_z')

      END SUBROUTINE M_bcast_z


!----------------------------------------------------------------------
!
!> copy n single precision complex from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_c_from(COMM, vec, n, inode )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(qs) vec(n)

      INTEGER inode, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_c_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_s(vec(1), n*2, inode-1, COMM, ACC_ASYNC_Q)
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_complex, inode-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_bcast_c_from','MPI_bcast')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_complex, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_c_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_complex, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_c_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast
      END SUBROUTINE M_bcast_c_from


!----------------------------------------------------------------------
!
!> copy n double precision complex from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z_from(COMM, vec, n, inode )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER inode, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_bcast_z_from')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_d(vec(1), n*2, inode-1, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_bcast_z_from')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_bcast_z_from','MPI_bcast')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         PROFILING_STOP('m_bcast_z_from')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_z_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_z_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast

      PROFILING_STOP('m_bcast_z_from')

      END SUBROUTINE M_bcast_z_from


!----------------------------------------------------------------------
!
!> copy n double precision complex from inode to all nodes the number
!> of items is INTEGER (KIND=8)
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z8_from(COMM, vec, n, inode )

      USE mpimy
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER (KIND=selected_int_kind(15)) :: n
      COMPLEX(q) vec(n)
! local
      INTEGER (KIND=selected_int_kind(15)) :: npos, nstep=2**24 ! nmax about 268 Mbyte bcast
      INTEGER :: nsend

      INTEGER inode, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_z8_from: invalid vector size", __FILE__, __LINE__)
      END IF

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         DO npos=1,n,nstep

            nsend=MIN(n-npos+1,nstep)

            CALL MPI_bcast( vec(npos), nsend,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
           &                ierror )
            CHECK_MPI_ERROR(ierror,'M_bcast_z8_from','MPI_bcast')
         ENDDO
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      DO npos=1,n,nstep

         nsend=MIN(n-npos+1,nstep)

         CALL MPI_ibcast( vec(npos), nsend,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
              &          request, ierror )
         CHECK_MPI_ERROR(ierror,'M_bcast_z8_from','MPI_ibcast')

         CALL MPI_wait( request, status, ierror )
      ENDDO
#else // MPI_avoid_bcast
      DO npos=1,n,nstep

         nsend=MIN(n-npos+1,nstep)

         CALL MPI_bcast( vec(npos), nsend,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
              &        ierror )
         CHECK_MPI_ERROR(ierror,'M_bcast_z8_from','MPI_bcast')
      ENDDO
#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast
      END SUBROUTINE M_bcast_z8_from


!----------------------------------------------------------------------
!
!> copy n single precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_s_from(COMM, vec, n, inode )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)

      INTEGER inode, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_s_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_real, inode-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_real, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_s_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_real, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_s_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast
      END SUBROUTINE M_bcast_s_from


!----------------------------------------------------------------------
!
!> copy n double precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d_from(COMM, vec, n, inode )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER inode, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_d_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_bcast_d_from')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_d(vec(1), n, inode-1, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_bcast_d_from')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         CALL MPI_bcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
        &                ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA

         CHECK_MPI_ERROR(ierror,'M_bcast_d_from','MPI_bcast')
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         PROFILING_STOP('m_bcast_d_from')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      CALL MPI_ibcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_d_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else // MPI_avoid_bcast
      CALL MPI_bcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_d_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast

      PROFILING_STOP('m_bcast_d_from')

      END SUBROUTINE M_bcast_d_from


!-----------------------------------------------------------------------
!
!> copy n double precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d8_from(COMM, vec, n, inode )

      USE mpimy
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER (KIND=selected_int_kind(15)) :: n
      REAL(q) vec(n)
! local
      INTEGER (KIND=selected_int_kind(15)) :: npos, nstep=2**24 ! nmax about 268 Mbyte bcast
      INTEGER :: nsend

      INTEGER inode, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_d8_from: invalid vector size", __FILE__, __LINE__)
      END IF

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
         DO npos=1,n,nstep

            nsend=MIN(n-npos+1,nstep)

            CALL MPI_bcast( vec(npos), nsend,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
           &                ierror )
            CHECK_MPI_ERROR(ierror,'M_bcast_d8_from','MPI_bcast')
         ENDDO
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
#ifdef MPI_bcast_with_barrier
         CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
         RETURN
      ENDIF
#endif // OFFLOADING

#ifdef MPI_avoid_bcast
      DO npos=1,n,nstep

         nsend=MIN(n-npos+1,nstep)

         CALL MPI_ibcast( vec(npos), nsend,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
              &          request, ierror )
         CHECK_MPI_ERROR(ierror,'M_bcast_d8_from','MPI_ibcast')

         CALL MPI_wait( request, status, ierror )
      ENDDO
#else // MPI_avoid_bcast
      DO npos=1,n,nstep

         nsend=MIN(n-npos+1,nstep)

         CALL MPI_bcast( vec(npos), nsend,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
              &                ierror )
         CHECK_MPI_ERROR(ierror,'M_bcast_d8_from','MPI_bcast')
      ENDDO
#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif // MPI_avoid_bcast
      END SUBROUTINE M_bcast_d8_from


!----------------------------------------------------------------------
!
!> copy n double precision complex from inode to all nodes (non-blocking)
!
!----------------------------------------------------------------------

      SUBROUTINE M_ibcast_z_from(COMM, vec, n, inode, request)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER inode
      INTEGER request

      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_ibcast_z_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_ibcast_z_from')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_d(vec(1), n*2, inode-1, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_ibcast_z_from')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
#ifndef MPI_avoid_ibcast_from
         CALL MPI_ibcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
        &                 request, ierror )
#else
         CALL MPI_bcast ( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
        &                          ierror )
         request = MPI_REQUEST_NULL
#endif
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_ibcast_z_from','MPI_ibcast')
         PROFILING_STOP('m_ibcast_z_from')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifndef MPI_avoid_ibcast_from
      CALL MPI_ibcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                 request, ierror )
#else
      CALL MPI_bcast ( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                          ierror )
      request = MPI_REQUEST_NULL
#endif
      CHECK_MPI_ERROR(ierror,'M_ibcast_z_from','MPI_ibcast')

      PROFILING_STOP('m_ibcast_z_from')

      END SUBROUTINE M_ibcast_z_from


!----------------------------------------------------------------------
!
!> copy n double precision from inode to all nodes (non-blocking)
!
!----------------------------------------------------------------------

      SUBROUTINE M_ibcast_d_from(COMM, vec, n, inode, request)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER inode
      INTEGER request

      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_ibcast_d_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_ibcast_d_from')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_bcast_d(vec(1), n, inode-1, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_ibcast_d_from')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
#ifndef MPI_avoid_ibcast_from
         CALL MPI_ibcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
        &                 request, ierror )
#else
         CALL MPI_bcast ( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
        &                          ierror )
         request = MPI_REQUEST_NULL
#endif
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_ibcast_d_from','MPI_ibcast')
         PROFILING_STOP('m_ibcast_d_from')
         RETURN
      ENDIF
#endif // OFFLOADING

#ifndef MPI_avoid_ibcast_from
      CALL MPI_ibcast( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                 request, ierror )
#else
      CALL MPI_bcast ( vec(1), n,  MPI_double_precision, inode-1, COMM%MPI_COMM, &
     &                          ierror )
      request = MPI_REQUEST_NULL
#endif
      CHECK_MPI_ERROR(ierror,'M_ibcast_d_from','MPI_ibcast')

      PROFILING_STOP('m_ibcast_d_from')

      END SUBROUTINE M_ibcast_d_from


!-----------------------------------------------------------------------
!
!> copy n quad precision from inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_qd_from(COMM, vec, n, inode )

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n, inode
      QDPREAL vec(n)

      INTEGER i, ierror
#ifdef MPI_avoid_bcast
      INTEGER request, status(MPI_status_size)
#endif
      PROFILING_START('m_bcast_qd_from' )
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_bcast_qd_from: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF
#ifdef MPI_avoid_bcast
#ifndef qd_emulate
      CALL MPI_ibcast( vec(1), n, MPI_real16, inode-1, COMM%MPI_COMM, &
     &                request, ierror )
      CHECK_MPI_ERROR(ierror,'M_bcast_qd_from','MPI_ibcast')

      CALL MPI_wait( request, status, ierror )
#else
      ! in case quadruple precision is emulated with qdmodule:
      ! vec(n) is an array of type qd_real and the type qd_real
      ! is a sequence of real*8  4-dimensional arrays re(4)
      ! so we broad cast following
      DO i = 1, n
         CALL MPI_ibcast( vec(i)%re(1), 4, MPI_real8, inode-1, COMM%MPI_COMM, &
        &                request, ierror )
         CHECK_MPI_ERROR(ierror,'M_bcast_qd_from','MPI_ibcast')

         CALL MPI_wait( request, status, ierror )
      ENDDO
#endif
#else
#ifndef qd_emulate
      CALL MPI_bcast( vec(1), n, MPI_real16, inode-1, COMM%MPI_COMM, &
     &                ierror )
#else
      ! in case quadruple precision is emulated with qdmodule:
      ! vec(n) is an array of type qd_real and the type qd_real
      ! is a sequence of real*8  4-dimensional arrays re(4)
      ! so we broad cast following
      DO i = 1, n
         CALL MPI_bcast( vec(i)%re(1), 4, MPI_real8, inode-1, COMM%MPI_COMM, &
        &                ierror )
      ENDDO
#endif
      CHECK_MPI_ERROR(ierror,'M_bcast_qd_from','MPI_bcast')

#ifdef MPI_bcast_with_barrier
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
#endif
#endif
      PROFILING_STOP('m_bcast_qd_from' )

      END SUBROUTINE M_bcast_qd_from

!======================================================================
!
! Global Exchange Routine
!
!======================================================================

!----------------------------------------------------------------------
!
!> complex global exchange routine
!>
!> This maps directly onto MPI_alltoallv which is usually very slow;
!> therefore an alternative implementation optimised for clusters exists
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_z(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv)

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      COMPLEX(q) xsnd(*)         !< send buffer
      COMPLEX(q) xrcv(*)         !< receive buffer

      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
      INTEGER nsnd(COMM%NCPU+1) !< number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) !< number of data recv from each node

!----------------------------------------------------------------------
#ifdef use_collective
!----------------------------------------------------------------------
      INTEGER r, ierror

      PROFILING_START('m_alltoallv_z')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(xsnd).AND.IS_ON_DEVICE(xrcv).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (COMM%LGRPNCCL) CALL nccl_group_start()
            DO r=1,COMM%NCPU
               IF(nsnd(r)>0) CALL nccl_send_d(xsnd(psnd(r)+1), nsnd(r)*2, r-1, COMM, ACC_ASYNC_Q)
               IF(nrcv(r)>0) CALL nccl_recv_d(xrcv(prcv(r)+1), nrcv(r)*2, r-1, COMM, ACC_ASYNC_Q)
            ENDDO
            IF (COMM%LGRPNCCL) CALL nccl_group_end()
            PROFILING_STOP('m_alltoallv_z')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(xsnd,xrcv)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(xsnd,xrcv)
         CALL MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_double_complex, &
                             xrcv(1), nrcv(1), prcv(1), MPI_double_complex, &
                             COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_alltoallv_z','MPI_alltoallv')
         PROFILING_STOP('m_alltoallv_z')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_double_complex, &
                          xrcv(1), nrcv(1), prcv(1), MPI_double_complex, &
                          COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_alltoallv_z','MPI_alltoallv')

!----------------------------------------------------------------------
#else // use_collective
!----------------------------------------------------------------------
      INTEGER ierror,sndcount,rcvcount,prcv_,psnd_,i,in
      INTEGER :: tag=201
      INTEGER :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=MPI_BLOCK/2
      INTEGER       :: block, p, nstat
      INTEGER       :: maxsnd_rcvcount

      PROFILING_START('m_alltoallv_z')

      maxsnd_rcvcount=MAX(MAXVAL(nsnd(1:COMM%NCPU)),MAXVAL(nrcv(1:COMM%NCPU)))

    DO block = 0, maxsnd_rcvcount-1, max_
      p        = 1 + block   ! pointer to the current block base address
      nstat    = 0

      ! initiate the receive and send on all nodes
      ! local copy is done below
      ! handle remaining  NCPU-1 packages
      DO in = 1, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         rcvcount=MIN(max_, nrcv(i+1)-block)
         prcv_   =prcv(i+1)

         IF (rcvcount>0) THEN
            nstat=nstat+1
            CALL MPI_irecv( xrcv(prcv_+p), rcvcount, MPI_double_complex, &
                       i, tag,  COMM%MPI_COMM, request(nstat), ierror )
            CHECK_MPI_ERROR(ierror,'M_alltoallv_z','MPI_irecv')
         ENDIF
      ENDDO

      ! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         sndcount=MIN(max_, nsnd(i+1)-block)
         psnd_   =psnd(i+1)

         IF (sndcount>0) THEN
            nstat=nstat+1
            CALL MPI_isend( xsnd(psnd_+p), sndcount, MPI_double_complex, &
     &                  i, tag,  COMM%MPI_COMM, request(nstat), ierror )
            CHECK_MPI_ERROR(ierror,'M_alltoallv_z','M_isend')
         ENDIF
      ENDDO

     ! local memory copy for data kept on the local node

      sndcount=MIN(max_, nsnd(COMM%NODE_ME)-block)
      prcv_   =prcv(COMM%NODE_ME)
      psnd_   =psnd(COMM%NODE_ME)

      IF (sndcount>0) &
         CALL ZCOPY( sndcount, xsnd( psnd_+p), 1, xrcv( prcv_+p) , 1 )

      CALL MPI_waitall(nstat , request, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoallv_z','MPI_waitall')

    ENDDO
!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoallv_z','MPI_barrier')
!----------------------------------------------------------------------
#endif // use_collective
!----------------------------------------------------------------------

      PROFILING_STOP('m_alltoallv_z')

      END SUBROUTINE M_alltoallv_z


!----------------------------------------------------------------------
!
!> real cyclic exchange routine which maps directly onto
!
!----------------------------------------------------------------------

      SUBROUTINE M_cycle_d(COMM, xsnd, nsnd)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"


      TYPE(communic) COMM
      REAL(q) xsnd(nsnd)         !< send/ receive buffer
      INTEGER nsnd

      INTEGER ierror, i, j, ichunk
      INTEGER :: tag=202
      INTEGER :: request(2)
      INTEGER :: A_OF_STATUSES(2)
      INTEGER :: NDTMP_=3

      DO j = 1, nsnd, NDTMP_
         ichunk = MIN( nsnd-j+1 , NDTMP_)

         ! initiate the receive from node-1 (note zero based)
         i = MOD(COMM%NODE_ME-1-1+COMM%NCPU , COMM%NCPU)
         CALL MPI_irecv( DTMP_m(1), ichunk, MPI_double_precision, &
              i, tag,  COMM%MPI_COMM, request(1), ierror )
         CHECK_MPI_ERROR(ierror,'M_cycle_d','MPI_irecv')

         ! initiate the send
         i = MOD(COMM%NODE_ME-1+1 , COMM%NCPU)  ! i zero based
         CALL MPI_isend( xsnd(j), ichunk, MPI_double_precision, &
              i, tag,  COMM%MPI_COMM, request(2), ierror )
         CHECK_MPI_ERROR(ierror,'M_cycle_d','MPI_isend')

         ! wait for send and receive to finish
         CALL MPI_waitall(2 , request, A_OF_STATUSES, ierror)
         CHECK_MPI_ERROR(ierror,'M_cycle_d','MPI_waitall')

         ! copy result back
         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  xsnd(j) , 1)
      END DO

      END SUBROUTINE M_cycle_d


!----------------------------------------------------------------------
!
!> integer routine which maps directly onto MPI_alltoallv
!>
!> this is used only once by VASP
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_i(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER xsnd(*)           !< send buffer
      INTEGER xrcv(*)           !< receive buffer

      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
      INTEGER nsnd(COMM%NCPU+1) !< number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) !< number of data recv from each node

! local data
      INTEGER ierror

      CALL MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_integer, &
                          xrcv(1), nrcv(1), prcv(1), MPI_integer, &
                          COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_alltoall_i','MPI_alltoallv')

      END SUBROUTINE M_alltoall_i


!----------------------------------------------------------------------
!
!> 8-byte integer routine which maps directly onto MPI_alltoallv
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_i8(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER(qi8) xsnd(*)           !< send buffer
      INTEGER(qi8) xrcv(*)           !< receive buffer

      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
      INTEGER nsnd(COMM%NCPU+1) !< number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) !< number of data recv from each node

! local data
      INTEGER ierror

      CALL MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_integer8, &
                          xrcv(1), nrcv(1), prcv(1), MPI_integer8, &
                          COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_alltoall_i8','MPI_alltoallv')

      END SUBROUTINE M_alltoall_i8


!----------------------------------------------------------------------
!
!> requires as only input the number of data nsnd send from each node to each other node
!>
!> it assumes a continous data arrangement on sender and receiver
!> and sets up the arrays which are required for MPI_alltoall
!> nrcv, psnd, prcv
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_simple(COMM, nsnd, nrcv, psnd, prcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
! input
      INTEGER nsnd(COMM%NCPU)   !< number of data send to each node
! output
      INTEGER psnd(COMM%NCPU+1) !< location of data in send buffer (0 based)
      INTEGER nrcv(COMM%NCPU)   !< number of data recv from each node
      INTEGER prcv(COMM%NCPU+1) !< location of data in recv buffer (0 based)
! local variable
      INTEGER ierror,i

      ! only one simple MPI_alltoall is required
      ! to find number of received data on each node
      CALL MPI_alltoall( nsnd(1),  1, MPI_integer, &
                         nrcv(1),  1, MPI_integer, &
                         COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_alltoallv_simple','MPI_alltoall')

      ! now set the locations assuming linear arrangement of data
      psnd(1)=0
      prcv(1)=0

      DO i=1,COMM%NCPU
        psnd(i+1)=psnd(i)+nsnd(i)
        prcv(i+1)=prcv(i)+nrcv(i)
      ENDDO

      END SUBROUTINE M_alltoallv_simple


!----------------------------------------------------------------------
!
!> real global exchange routine
!>
!> redistributes an array from distribution over bands to
!> distribution over coefficient (or vice versa)
!>      original distribution             final distribution
!>     |  1  |  2  |  3  |  4  |       |  1  |  1  |  1  |  1  |
!>     |  1  |  2  |  3  |  4  |       |  2  |  2  |  2  |  2  |
!>     |  1  |  2  |  3  |  4  |  <->  |  3  |  3  |  3  |  3  |
!>     |  1  |  2  |  3  |  4  |       |  4  |  4  |  4  |  4  |
!>
!> mind that only (n/NCPU) *NCPU data are exchanged
!> it is the responsability of the user to guarantee that n is
!> correct
!>
!> \param xsnd the array to be redistributed (having n elements)
!> \param xrcv the result array with n/NCPU elements received from each processor
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d(COMM, n, xsnd, xrcv )

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=MPI_BLOCK
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
                       proc_group, group_base, i_in_group, irequests
      INTEGER       :: r

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_alltoall_d')

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------
#if defined(use_collective) || defined(OFFLOADING)
!----------------------------------------------------------------------
#ifdef OFFLOADING
      IF (IS_ON_DEVICE(xsnd).AND.IS_ON_DEVICE(xrcv).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (COMM%LGRPNCCL) CALL nccl_group_start()
            DO r=0,COMM%NCPU-1
               IF(sndcount>0) CALL nccl_send_d(xsnd(r*sndcount+1), sndcount, r, COMM, ACC_ASYNC_Q)
               IF(rcvcount>0) CALL nccl_recv_d(xrcv(r*rcvcount+1), rcvcount, r, COMM, ACC_ASYNC_Q)
            ENDDO
            IF (COMM%LGRPNCCL) CALL nccl_group_end()
            PROFILING_STOP('m_alltoall_d')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(xsnd,xrcv)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(xsnd,xrcv)
         CALL MPI_alltoall( xsnd(1), sndcount, MPI_double_precision, &
        &                   xrcv(1), rcvcount, MPI_double_precision, &
        &                   COMM%MPI_COMM, ierror )
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_alltoall')
         PROFILING_STOP('m_alltoall_d')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_alltoall( xsnd(1), sndcount, MPI_double_precision, &
     &                   xrcv(1), rcvcount, MPI_double_precision, &
     &                   COMM%MPI_COMM, ierror )

      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_alltoall')

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_barrier')
!----------------------------------------------------------------------
#else // use_collective || OFFLOADING
#ifndef avoid_async
#ifndef PROC_GROUP
!
!  this version  throws all data in one go onto all nodes
!  then waits for all sends and receives to finish
!----------------------------------------------------------------------
      ! initiate the receive on all nodes
      ! local copy has already been done on each node handle remaining  NCPU-1 packages
      DO block = 0, sndcount-1, max_
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

      DO in = 1, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         CALL MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(in), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_irecv')
      ENDDO

      ! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         CALL MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(in+ COMM%NCPU-1), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_isend')
      ENDDO

      ! local memory copy for data kept on the local node
      ! overlaps with communication
      CALL DCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )

      CALL MPI_waitall((COMM%NCPU-1)*2, request, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_waitall')
      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_barrier')
!----------------------------------------------------------------------
#else // PROC_GROUP
!
!  this version  groups the processors into groups
!  with PROC_GROUP members each
!  communication is first done within one group and
!  then between any two groups
!----------------------------------------------------------------------
      actual_proc_group=MIN(PROC_GROUP,COMM%NCPU)
      ! comm_proc_group is allways larger than COMM%NCPU
      com_proc_group =((COMM%NCPU+actual_proc_group-1)/ actual_proc_group)*actual_proc_group

!      WRITE(*,*) COMM%NODE_ME, actual_proc_group,com_proc_group

      !
      ! blocking on the nodes
      !
      DO proc_group=0, com_proc_group -1, actual_proc_group
      !
      ! blocking on the data
      !
      DO block = 0, sndcount-1, max_
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

         irequests=0 ! counts the number of issued requests

         DO in = 0, actual_proc_group-1
           ! base index of the present group
           group_base = ((COMM%NODE_ME-1)/actual_proc_group)*actual_proc_group

           ! within the group send to node i_in_group
           i_in_group = MOD(in+COMM%NODE_ME-1 , actual_proc_group)

           ! actual node from which to receive (i is zero based)

           i = MOD(i_in_group+group_base+proc_group , com_proc_group)
!           WRITE(*,*) 'receive',i,COMM%NODE_ME-1

           ! local copies are handled below, and take care of sends from smaller groups
           IF ( i/=COMM%NODE_ME-1 .AND. i<COMM%NCPU ) THEN
              irequests=irequests+1
              CALL MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_double_precision, &
                   i, tag,  COMM%MPI_COMM, request(irequests), ierror )
              CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_irecv')
           ENDIF
         ENDDO

      ! initiate the send on all nodes
         DO in = 0, actual_proc_group-1
           ! base index of the present group
           group_base = ((COMM%NODE_ME-1)/actual_proc_group)*actual_proc_group

           ! send to node i_in_group
           i_in_group = MOD(in+COMM%NODE_ME-1 , actual_proc_group)

           ! actual node to which to send
           i = MOD(i_in_group+group_base-proc_group+com_proc_group , com_proc_group)

!           WRITE(*,*) 'send',COMM%NODE_ME-1,i

           ! local receives are handled below, also take care of receives from smaller groups
           IF ( i/=COMM%NODE_ME-1 .AND. i<COMM%NCPU ) THEN
              irequests=irequests+1
              CALL MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_double_precision, &
                   i, tag,  COMM%MPI_COMM, request(irequests), ierror )
              CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_isend')
           ENDIF
         ENDDO

      CALL MPI_waitall(irequests, request, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_waitall')

      ! local memory copy for data kept on the local node
      ! overlaps with communication
      CALL DCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )

      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_barrier')

      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_barrier')
!----------------------------------------------------------------------
#endif // PROC_GROUP
#else // avoid_async
! more moderate version copies data between two nodes and then blocks
! not so bad if we have full duplex point to point connections
!----------------------------------------------------------------------
      ! local memory copy for data kept on the local node
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

      ! initiate the receive and send on all nodes
      ! local copy has already been done on each node
      ! handle remaining  NCPU-1 packages
      DO in = 1, COMM%NCPU-1
         ! receive from: (own node id) - in
         i = MOD(-in+COMM%NODE_ME-1 +COMM%NCPU, COMM%NCPU)  ! i zero based
         CALL MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(1), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_irecv')

         !
         ! sent to: (own node id) + in
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based
         CALL MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(2), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_isend')

         CALL MPI_waitall(2, request, A_OF_STATUSES, ierror)
         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_waitall')

!         CALL MPI_barrier( COMM%MPI_COMM, ierror )
!         CHECK_MPI_ERROR(ierror,'M_alltoall_d','MPI_barrier')

      ENDDO
!----------------------------------------------------------------------
#endif // avoid_async
#endif // use_collective || OFFLOADING
!----------------------------------------------------------------------

      PROFILING_STOP('m_alltoall_d')

      END SUBROUTINE M_alltoall_d


      SUBROUTINE M_alltoall_d_omp(COMM, n, xsnd, xrcv )
      USE mpimy
      USE mopenmp_struct_def, ONLY : omp_nthreads_alltoall
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i

      INTEGER, SAVE :: tag=201
      INTEGER       :: in

      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
! test_
!     INTEGER, PARAMETER :: max_=MPI_BLOCK
      INTEGER       :: max_
! test_
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

!$    INTEGER       :: omp_id
!$    INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU
! test_
      max_=(sndcount+omp_nthreads_alltoall-1)/omp_nthreads_alltoall
! test_

!----------------------------------------------------------------------
!  this version  throws all data in one go onto all nodes
!  then waits for all sends and receives to finish
!----------------------------------------------------------------------
      ! initiate the receive on all nodes
      ! local copy has already been done on each node handle remaining  NCPU-1 packages
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(omp_nthreads_alltoall) &
!$OMP SHARED(sndcount,max_,COMM,xrcv,xsnd) PRIVATE(block,omp_id,sndcount_,p,in,i,request,ierror,A_OF_STATUSES)
      DO block = 0, sndcount-1, max_
!$    omp_id   = OMP_GET_THREAD_NUM()
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

      DO in = 1, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         CALL MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, __omp_id,  COMM%MPI_COMM, request(in), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_d_omp','MPI_irecv')
      ENDDO

      ! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         CALL MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, __omp_id,  COMM%MPI_COMM, request(in+ COMM%NCPU-1), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_d_omp','MPI_isend')
      ENDDO

      ! local memory copy for data kept on the local node
      ! overlaps with communication
      CALL DCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )

      CALL MPI_waitall((COMM%NCPU-1)*2, request, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_d_omp','MPI_waitall')
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE M_alltoall_d_omp


!----------------------------------------------------------------------
!
!> real global exchange routine (single precision)
!>
!> redistributes an array from distribution over bands to
!> distribution over coefficient (or vice versa)
!>      original distribution             final distribution
!>     |  1  |  2  |  3  |  4  |       |  1  |  1  |  1  |  1  |
!>     |  1  |  2  |  3  |  4  |       |  2  |  2  |  2  |  2  |
!>     |  1  |  2  |  3  |  4  |  <->  |  3  |  3  |  3  |  3  |
!>     |  1  |  2  |  3  |  4  |       |  4  |  4  |  4  |  4  |
!>
!> mind that only (n/NCPU) *NCPU data are exchanged
!> it is the responsability of the user to guarantee that n is
!> correct
!>
!> \param xsnd the array to be redistributed (having n elements)
!> \param xrcv the result array with n/NCPU elements received from each processor
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_s(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=MPI_BLOCK
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------
#ifdef use_collective
!----------------------------------------------------------------------
      CALL MPI_alltoall( xsnd(1), sndcount, MPI_real, &
     &                   xrcv(1), rcvcount, MPI_real, &
     &                   COMM%MPI_COMM, ierror )
      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_alltoall')

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_barrier')

!----------------------------------------------------------------------
#else
#ifndef avoid_async
#ifndef PROC_GROUP
!
!  this version  throws all data in one go onto all nodes
!  then waits for all sends and receives to finish
!----------------------------------------------------------------------
      ! initiate the receive on all nodes
      ! local copy has already been done on each node handle remaining  NCPU-1 packages
      DO block = 0, sndcount-1, max_
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

      DO in = 1, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         CALL MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_real, &
     &                  i, tag,  COMM%MPI_COMM, request(in), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_irecv')
      ENDDO

      ! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         CALL MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_real, &
     &                  i, tag,  COMM%MPI_COMM, request(in+ COMM%NCPU-1), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_isend')
      ENDDO

      ! local memory copy for data kept on the local node
      ! overlaps with communication
      CALL SCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )
      CALL MPI_waitall((COMM%NCPU-1)*2, request, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_waitall')
      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_barrier')
!----------------------------------------------------------------------
#else
!
!  this version  groups the processors into groups
!  with PROC_GROUP members each
!  communication is first done within one group and
!  then between any two groups
!----------------------------------------------------------------------
      actual_proc_group=MIN(PROC_GROUP,COMM%NCPU)
      ! comm_proc_group is allways larger than COMM%NCPU
      com_proc_group =((COMM%NCPU+actual_proc_group-1)/ actual_proc_group)*actual_proc_group

!      WRITE(*,*) COMM%NODE_ME, actual_proc_group,com_proc_group

      !
      ! blocking on the nodes
      !
      DO proc_group=0, com_proc_group -1, actual_proc_group
      !
      ! blocking on the data
      !
      DO block = 0, sndcount-1, max_
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

         irequests=0 ! counts the number of issued requests

         DO in = 0, actual_proc_group-1
           ! base index of the present group
           group_base = ((COMM%NODE_ME-1)/actual_proc_group)*actual_proc_group

           ! within the group send to node i_in_group
           i_in_group = MOD(in+COMM%NODE_ME-1 , actual_proc_group)

           ! actual node from which to receive (i is zero based)

           i = MOD(i_in_group+group_base+proc_group , com_proc_group)
!           WRITE(*,*) 'receive',i,COMM%NODE_ME-1

           ! local copies are handled below, and take care of sends from smaller groups
           IF ( i/=COMM%NODE_ME-1 .AND. i<COMM%NCPU ) THEN
              irequests=irequests+1
              CALL MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_real, &
                   i, tag,  COMM%MPI_COMM, request(irequests), ierror )
              CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_irecv')
           ENDIF
         ENDDO

      ! initiate the send on all nodes
         DO in = 0, actual_proc_group-1
           ! base index of the present group
           group_base = ((COMM%NODE_ME-1)/actual_proc_group)*actual_proc_group

           ! send to node i_in_group
           i_in_group = MOD(in+COMM%NODE_ME-1 , actual_proc_group)

           ! actual node to which to send
           i = MOD(i_in_group+group_base-proc_group+com_proc_group , com_proc_group)

!           WRITE(*,*) 'send',COMM%NODE_ME-1,i

           ! local receives are handled below, also take care of receives from smaller groups
           IF ( i/=COMM%NODE_ME-1 .AND. i<COMM%NCPU ) THEN
              irequests=irequests+1
              CALL MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_real, &
                   i, tag,  COMM%MPI_COMM, request(irequests), ierror )
              CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_isend')
           ENDIF
         ENDDO

      CALL MPI_waitall(irequests, request, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_waitall')

      ! local memory copy for data kept on the local node
      ! overlaps with communication
      CALL SCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )

      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_barrier')

      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_barrier')
!----------------------------------------------------------------------
#endif
#else
! more moderate version copies data between two nodes and then blocks
! not so bad if we have full duplex point to point connections
!----------------------------------------------------------------------
      ! local memory copy for data kept on the local node
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

      ! initiate the receive and send on all nodes
      ! local copy has already been done on each node
      ! handle remaining  NCPU-1 packages
      DO in = 1, COMM%NCPU-1
         ! receive from: (own node id) - in
         i = MOD(-in+COMM%NODE_ME-1 +COMM%NCPU, COMM%NCPU)  ! i zero based
         CALL MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_real, &
     &                  i, tag,  COMM%MPI_COMM, request(1), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_irecv')

         !
         ! sent to: (own node id) + in
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based
         CALL MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_real, &
     &                  i, tag,  COMM%MPI_COMM, request(2), ierror )
         CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_isend')

         CALL MPI_waitall(2, request, A_OF_STATUSES, ierror)
         CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_waitall')

!         CALL MPI_barrier( COMM%MPI_COMM, ierror )
!         CHECK_MPI_ERROR(ierror,'M_alltoall_s','MPI_barrier')

      ENDDO
!----------------------------------------------------------------------
#endif
#endif
!----------------------------------------------------------------------

      END SUBROUTINE M_alltoall_s


!----------------------------------------------------------------------
!
!> complex global exchange routine
!>
!> uses #M_alltoall_d with twice as many elements
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_z(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) xsnd(n), xrcv(n)

      CALL M_alltoall_d(COMM, n*2, xsnd, xrcv)

      END SUBROUTINE M_alltoall_z


!----------------------------------------------------------------------
!
!> real global exchange routine (nonblocking)
!>
!> redistributes an array from distribution over bands to
!> distribution over coefficient (or vice versa)
!>      original distribution             final distribution
!>     |  1  |  2  |  3  |  4  |       |  1  |  1  |  1  |  1  |
!>     |  1  |  2  |  3  |  4  |       |  2  |  2  |  2  |  2  |
!>     |  1  |  2  |  3  |  4  |  <->  |  3  |  3  |  3  |  3  |
!>     |  1  |  2  |  3  |  4  |       |  4  |  4  |  4  |  4  |
!>
!> mind that only (n/NCPU) *NCPU data are exchanged
!> it is the responsability of the user to guarantee that n is
!> correct
!>
!> \param xsnd the array to be redistributed (having n elements)
!> \param xrcv the result array with n/NCPU elements received from each processor
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d_async(COMM, n, xsnd, xrcv, tag, srequest, rrequest )

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER tag
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror
      INTEGER i,j,in

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_alltoall_d_async')

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

      ! local memory copy for data kept on the local node
GPACC !$ACC PARALLEL LOOP PRESENT(xrcv,xsnd) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:xrcv,xsnd) IF(OFFLOAD_ON)
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

      ! initiate send and receive on all nodes
      ! local copy has already been done each node send NCPU-1 packages
      j=1
#ifdef OFFLOADING
      IF (IS_ON_DEVICE(xsnd).AND.IS_ON_DEVICE(xrcv).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (COMM%LGRPNCCL) CALL nccl_group_start()
            DO in = 0, COMM%NCPU-1
               ! send to node in + own node id
               ! such a construct should allow for an efficient use of the network
               i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)
               IF ( COMM%NODE_ME-1 /= i ) THEN
                  IF(sndcount>0) CALL nccl_send_d(xsnd(i*sndcount+1), sndcount, i, COMM, ACC_ASYNC_Q)
                  IF(rcvcount>0) CALL nccl_recv_d(xrcv(i*rcvcount+1), rcvcount, i, COMM, ACC_ASYNC_Q)
                  srequest(j) = ACC_ASYNC_Q ; rrequest(j) = ACC_ASYNC_Q ; j = j + 1
               ENDIF
            ENDDO
            IF (COMM%LGRPNCCL) CALL nccl_group_end()
            PROFILING_STOP('m_alltoall_d_async')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(xsnd,xrcv)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(xsnd,xrcv)
            DO in = 0, COMM%NCPU-1
               ! send to node in + own node id
               ! such a construct should allow for an efficient use of the network
               i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)
               IF ( COMM%NODE_ME-1 /= i) THEN
                  CALL MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
           &                  i, tag,  COMM%MPI_COMM, srequest(j), ierror )
                  CHECK_MPI_ERROR(ierror,'M_alltoall_d_async','MPI_isend')
                  CALL MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
           &                  i, tag,  COMM%MPI_COMM, rrequest(j), ierror )
                  CHECK_MPI_ERROR(ierror,'M_alltoall_d_async','MPI_irecv')
                  j=j+1
               ENDIF
            ENDDO
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_alltoall_d_async')
         RETURN
      ENDIF
#endif // OFFLOADING

      DO in = 0, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network

         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)
         IF ( COMM%NODE_ME-1 /= i) THEN
            CALL MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, srequest(j), ierror )
            CHECK_MPI_ERROR(ierror,'M_alltoall_d_async','MPI_isend')
            CALL MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, rrequest(j), ierror )
            CHECK_MPI_ERROR(ierror,'M_alltoall_d_async','MPI_irecv')
            j=j+1
         ENDIF
      ENDDO

      PROFILING_STOP('m_alltoall_d_async')

      END SUBROUTINE M_alltoall_d_async


      SUBROUTINE M_alltoall_wait(COMM, srequest, rrequest )

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      INTEGER ierror
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,COMM%NCPU)

      ! wait for the NCPU-1 outstanding packages

#if defined(OFFLOADING) && defined(USENCCL)
      IF (OFFLOAD_ON .AND. COMM%LUSENCCL) THEN
         ! We assume that all send and receives associated
         ! with srequest and rrequest were put into the same
         ! asynchronous execution queue. If this is not the
         ! case the single wait below will not suffice ...
GPACC !$ACC WAIT(rrequest(1))
         RETURN
      ENDIF
#endif
      CALL MPI_waitall(COMM%NCPU-1, srequest, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_wait','MPI_waitall (1)')
      CALL MPI_waitall(COMM%NCPU-1, rrequest, A_OF_STATUSES, ierror)
      CHECK_MPI_ERROR(ierror,'M_alltoall_wait','MPI_waitall (2)')

      END SUBROUTINE M_alltoall_wait


!----------------------------------------------------------------------
!
!> performs a fast global sum on n doubles in vector vec
!  (algorithm by Kresse Georg)
!>
!> uses complete interchange algorithm
!  (my own invention, but I guess some people must know it)
!> exchange data between nodes, sum locally and
!> interchange back, this algorithm is faster than typical MPI based
!> algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_d(COMM, vec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(q) vec(n)

      REAL(q), ALLOCATABLE :: vec_inter(:)
    ! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=MPI_BLOCK

    ! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
     ! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
        ! sum localy
         DO i=2, COMM%NCPU
           CALL DAXPY(ncount, 1.0_q, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
     ! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
     ! backward exchange
         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL DCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

     ! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_d: " // str(n_) // " " // str(nsummed), __FILE__, __LINE__)
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_d(COMM, vec(nsummed+1), n-nsummed)

      DEALLOCATE(vec_inter)

      END SUBROUTINE M_sumf_d


!----------------------------------------------------------------------
!
!> performs a fast global sum on n doubles in vector vec
! (algorithm by Kresse Georg)
!>
!> This a special version for giant arrays exceeding an element count n
!> of 2^31-1 (the maximum size that can be handled with 32-bit INTEGER).
!> Mainly, the major difference to routine #M_sumf_d is that the last
!> argument ("n") is now declared as 64-bit INTEGER -- and that
!> in addition some intermediate integer operation need also to
!> be done with 64-bit INTEGER temporary variables
!  (adapted by jF)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_d8(COMM, vec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
! this must all declared as 64-bit INTEGER
      INTEGER(qi8) n,nsummed,ndo,i,j,n_,mmax
! this must stay to be 32-bit INTEGER and a new "n4" is needed for MPI calls
      INTEGER ncount,info,n4
      REAL(q) vec(n)

      REAL(q), ALLOCATABLE :: vec_inter(:)
    ! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=MPI_BLOCK

    ! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n/COMM%NCPU < 2147483648_qi8) THEN
! here the case that n/COMM%CPU still fits into the 32-bit INTEGER range ...
         n4=n/COMM%NCPU
         mmax=MIN(n4,max_)
      ELSE
! ... if not: (since max_ is of type 32-bit INTEGER) max_ *must* be the minimum
         mmax=max_
      ENDIF
      ALLOCATE(vec_inter(mmax*COMM%NCPU))
GPACC !$ACC ENTER DATA CREATE(vec_inter) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:vec_inter) IF(OFFLOAD_ON)

      nsummed=0_qi8
      n_=n/COMM%NCPU

      DO ndo=0_qi8,n_-1_qi8,mmax
     ! forward exchange
     ! since the maximum value of mmax is finally determined by max_ (MPI_BLOCK)
     ! which is a 32-bit INTEGER (and hence can never violate the allowed range
     ! of 32-bit INTEGER numbers) I could safely assume that "ncount" may be of
     ! type "INTEGER" -- this was essential since DAXPY and DCOPY as well as
     ! all MPI calls require 32-bit INTEGER arguments (for a standard BLAS);
     ! be warned that the situation might change if COMM%NCPU goes into the
     ! millions (since actually not ncount alone but ncount*COMM%NCPU is handed
     ! over as argument!) but then not only here but in ALL other subroutines we
     ! would run into trouble -- or needed BLAS/MPI routines that can accept
     ! 64-bit INTEGERS (still rare but already existing ...); at the moment I
     ! could just recommend to set a small enough value for max_ (MPI_BLOCK)
     ! which must be simply smaller than INT(2147483647/NCPU) ...
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
        ! sum localy
         DO i=2, COMM%NCPU
           CALL WDAXPY(ncount, 1.0_q, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
     ! replicate data (will be send to each proc)
GPACC !$ACC PARALLEL LOOP PRESENT(vec,vec_inter) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:vec,vec_inter) IF(OFFLOAD_ON)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
     ! backward exchange
         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL WDCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

     ! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_d8: " // str(n_) // " " // str(nsummed), __FILE__, __LINE__)
      ENDIF

! Now the remaining few elements (we can use M_sumb_d for this task since the
! number of elements left cannot be larger than "max_" [being a 32-bit INTEGER])
      IF (n-nsummed /= 0_qi8 ) THEN
        n4=n-nsummed
        CALL M_sumb_d(COMM, vec(nsummed+1), n4)
      ENDIF

GPACC !$ACC EXIT DATA DELETE(vec_inter) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:vec_inter) IF(OFFLOAD_ON)
      DEALLOCATE(vec_inter)

      END SUBROUTINE M_sumf_d8


!----------------------------------------------------------------------
!
!> performs a fast global sum on n singles in vector vec
!  (algorithm by Kresse Georg)
!>
!> uses complete interchange algorithm
!  (my own invention, but I guess some people must know it)
!> exchange data between nodes, sum locally and
!> interchange back, this algorithm is faster than typical MPI based
!> algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_s(COMM, vec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(qs) vec(n)
      REAL(qs), ALLOCATABLE :: vec_inter(:)
    ! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=MPI_BLOCK

    ! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
     ! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
        ! sum localy
         DO i=2, COMM%NCPU
           CALL SAXPY(ncount, 1.0, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
     ! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
     ! backward exchange
         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL SCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

     ! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_s: " // str(n_) // " " // str(nsummed), __FILE__, __LINE__)
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_s(COMM, vec(nsummed+1), n-nsummed)

      DEALLOCATE(vec_inter)
      END SUBROUTINE M_sumf_s


!----------------------------------------------------------------------
!
!> performs a fast global sum on n singles in vector vec
!  (algorithm by Kresse Georg)
!> This a special version for giant arrays exceeding an element count n of
!> 2^31-1 (the maximum size that can be handled with 32-bit INTEGER).
!> Mainly, the major difference to routine M_sumf_s is that the last
!> argument ("n") is now declared as 64-bit INTEGER -- and that
!> in addition some intermediate integer operation need also to
!> be done with 64-bit INTEGER temporary variables
!  (adapted by jF)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_s8(COMM, vec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(communic) COMM
! this must all declared as 64-bit INTEGER
      INTEGER(qi8) n,nsummed,ndo,i,j,n_,mmax
! this must stay to be 32-bit INTEGER and a new "n4" is needed for MPI calls
      INTEGER ncount,info,n4
      REAL(qs) vec(n)
      REAL(qs), ALLOCATABLE :: vec_inter(:)
    ! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=MPI_BLOCK

    ! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n/COMM%NCPU < 2147483648_qi8) THEN
! here the case that n/COMM%CPU still fits into the 32-bit INTEGER range ...
         n4=n/COMM%NCPU
         mmax=MIN(n4,max_)
      ELSE
! ... if not: (since max_ is of type 32-bit INTEGER) max_ *must* be the minimum
         mmax=max_
      ENDIF
      ALLOCATE(vec_inter(mmax*COMM%NCPU))

      nsummed=0_qi8
      n_=n/COMM%NCPU

      DO ndo=0_qi8,n_-1_qi8,mmax
     ! forward exchange
! same story as in M_sumf_d8: we can again safely keep ncount to be INTEGER
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
        ! sum localy
         DO i=2, COMM%NCPU
           CALL SAXPY(ncount, 1.0, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
     ! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
     ! backward exchange
         CALL M_alltoall_s(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
         CALL SCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
      ENDDO

     ! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         CALL vtutor%bug("M_sumf_s8: " // str(n_) // " " // str(nsummed), __FILE__, __LINE__)
      ENDIF

! Now the remaining few elements (we can use M_sumb_s for this task since the
! number of elements left cannot be larger than "max_" [being a 32-bit INTEGER])
      IF (n-nsummed /= 0_qi8 ) THEN
        n4=n-nsummed
        CALL M_sumb_s(COMM, vec(nsummed+1), n4)
      ENDIF

      DEALLOCATE(vec_inter)
      END SUBROUTINE M_sumf_s8


!----------------------------------------------------------------------
!
!> performs a fast global sum on n (double) complex in vector 'vec' (see #M_sumf_d)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(2*n)
      CALL M_sumf_d(COMM, vec, 2*n)
      END SUBROUTINE M_sumf_z


!----------------------------------------------------------------------
!
! performs a fast global sum on n (single) complex in vector 'vec' (see #M_sumf_s)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_c(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(2*n)
      CALL M_sumf_s(COMM, vec, 2*n)
      END SUBROUTINE M_sumf_c


!----------------------------------------------------------------------
!
!> performs a sum on n double complex numbers; it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(2*n)

#ifdef use_collective_sum
      CALL M_sumb_d(COMM, vec, 2*n)
#else
      IF ( 2*n>MPI_BLOCK) THEN
         CALL M_sumf_d(COMM, vec, 2*n)
      ELSE
         CALL M_sumb_d(COMM, vec, 2*n)
      ENDIF
#endif
      END SUBROUTINE M_sum_z


!----------------------------------------------------------------------
!
!> performs a sum on n double complex numbers
!>
!> This is for giant arrays with giant dimensions (parameter 2*n potentially
!> exceeding the range of 32-bit INTEGER numbers); this uses the special
!> version M_sum_d8 of routine M_sum_d and just uses 2*n as the
!> corresponding number of doubles to be communicated (like in M_sumf_z)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_z8(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER(qi8) n
      REAL(q) vec(2_qi8*n)

      CALL M_sum_d8(COMM, vec, 2_qi8*n)

      END SUBROUTINE M_sum_z8


!----------------------------------------------------------------------
!
!> performs a sum on n single complex numbers; it uses either sumb_s or sumf_s
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_c(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(2*n)

#ifdef use_collective_sum
      CALL M_sumb_s(COMM, vec, 2*n)
#else
      IF ( 2*n>MPI_BLOCK) THEN
         CALL M_sumf_s(COMM, vec, 2*n)
      ELSE
         CALL M_sumb_s(COMM, vec, 2*n)
      ENDIF
#endif
      END SUBROUTINE M_sum_c


!----------------------------------------------------------------------
!
!> performs a sum on n quadruple complex numbers
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_qd(COMM, vec, n)

      USE mpimy
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic), INTENT(IN) :: COMM
      INTEGER, INTENT(IN) :: n
      QDPREAL             :: vec(n)

      INTEGER  ierror

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         CALL vtutor%bug("M_sum_qd: invalid vector size n " // str(n), __FILE__, __LINE__)
      END IF

      PROFILING_START('m_sum_qd')

#ifndef qd_emulate
! invoke in-place version of MPI_allreduce, with own operation
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, MPI_real16, &
         &                M_sum_qd_op, COMM%MPI_COMM, ierror )
#else
! invoke in-place version of MPI_allreduce
      CALL MPI_allreduce( MPI_IN_PLACE, vec(1), n, M_qd_type, &
         &                M_sum_qd_op, COMM%MPI_COMM, ierror )
#endif
      CHECK_MPI_ERROR(ierror,'M_sum_qd','MPI_allreduce')

      PROFILING_STOP('m_sum_qd')

      END SUBROUTINE M_sum_qd


!----------------------------------------------------------------------
!
!> performs a sum on n double prec numbers; it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

#ifdef use_collective_sum
      CALL M_sumb_d(COMM, vec, n)
#else
      IF ( n>MPI_BLOCK) THEN
         CALL M_sumf_d(COMM, vec, n)
      ELSE
         CALL M_sumb_d(COMM, vec, n)
      ENDIF
#endif
      END SUBROUTINE M_sum_d


!----------------------------------------------------------------------
!
!> performs a sum on n double precision numbers
!>
!> This is for giant arrays with giant dimensions (parameter n potentially
!> exceeding the range of 32-bit INTEGER numbers); this mainly uses special
!> version M_sumf_d_giant of routine M_sumf_d (use of M_sumb_d not always
!> possible for "use_collective_sum" -- only second call always possible)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_d8(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER(qi8) n
      REAL(q) vec(n)
! we need additional local variables
      INTEGER(qi8) max_
      INTEGER n4

! return ranks that are not in the group
      IF ( COMM%MPI_COMM == MPI_COMM_NULL ) THEN
         RETURN
      ENDIF

#ifdef use_collective_sum
! this is now a bit more complicated than in routine M_sum_d ...
      max_=2147483647_qi8          ! 2^31-1 = maximum possible 32-bit integer
      IF ( n>max_) THEN
! here we have no choice; we cannot use M_sumb_d (MPI limitations to 32-bit
! integers!) and therefore we have to use M_sumf_d8 instead ... !!
         CALL M_sumf_d8(COMM, vec, n)
      ELSE
! only "n" smaller than 2^31 (limit for 32-bit integers) allows use of M_sumb_d
         n4=n
         CALL M_sumb_d(COMM, vec, n4)
      ENDIF
#else
! ... and also here we need to take a bit care ...
      max_=MPI_BLOCK
      IF (n>max_) THEN
! ... for safety we also have to call the "8" version here ...
         CALL M_sumf_d8(COMM, vec, n)
      ELSE
! ... and since MPI_BLOCK (which is a 32-bit INTEGER) is the possible maximum
! of n inside this ELSE branch we can keep the call to M_sumb_d here ...
         n4=n
         CALL M_sumb_d(COMM, vec, n4)
      ENDIF
#endif
      END SUBROUTINE M_sum_d8


!----------------------------------------------------------------------
!
!> performs a sum on n single prec numbers; it uses either sumb_s or sumf_s
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_single(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(qs) vec(n)

#ifdef use_collective_sum
      CALL M_sumb_s(COMM, vec, n)
#else
      IF ( n>MPI_BLOCK) THEN
         CALL M_sumf_s(COMM, vec, n)
      ELSE
         CALL M_sumb_s(COMM, vec, n)
      ENDIF
#endif
      END SUBROUTINE M_sum_single


!----------------------------------------------------------------------
!
!> performs a sum on n single prec numbers
!>
!> This is for giant arrays with giant dimensions (parameter n potentially
!> exceeding the range of 32-bit INTEGER numbers); this mainly uses special
!> version M_sumf_s8 of routine M_sumf_s (use of M_sumb_s not always
!> possible for "use_collective_sum" -- only second call always possible);
!> the structure is basically the same as in M_sum_d8 ("same story")
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_single8(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER(qi8) n
      REAL(qs) vec(n)
! local
      INTEGER(qi8) max_
      INTEGER n4

#ifdef use_collective_sum
      max_=2147483647_qi8
      IF ( n>max_) THEN
         CALL M_sumf_s8(COMM, vec, n)
      ELSE
         n4=n
         CALL M_sumb_s(COMM, vec, n4)
      ENDIF
#else
      max_=MPI_BLOCK
      IF ( n>MPI_BLOCK) THEN
         CALL M_sumf_s8(COMM, vec, n)
      ELSE
         n4=n
         CALL M_sumb_s(COMM, vec, n4)
      ENDIF
#endif
      END SUBROUTINE M_sum_single8


!----------------------------------------------------------------------
!
!> performs a sum on n double prec numbers to master
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_master_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER :: I,thisdata,ierror

      IF (COMM%NCPU == 1 ) RETURN

      DO i=1,n,MPI_BLOCK
         thisdata=min(n-i+1,MPI_BLOCK)

         CALL MPI_reduce( vec(i), DTMP_m(1), thisdata, &
                             MPI_double_precision, MPI_sum, &
                             0, COMM%MPI_COMM, ierror )
         CHECK_MPI_ERROR(ierror,'M_sum_master_d','MPI_reduce')

         IF (COMM%NODE_ME==1) THEN
            vec(i:i+thisdata-1)=DTMP_m(1:thisdata)
         ENDIF
      END DO

      END SUBROUTINE M_sum_master_d


!----------------------------------------------------------------------
!
!> performs a sum on n double complex numbers to master
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_master_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      CALL M_sum_master_d(COMM,vec,2*n)

      END SUBROUTINE M_sum_master_z

!----------------------------------------------------------------------
!
!> copy nrcv double complex from node i=1,COMM%NCPU to x((i-1)*nrcv+1:i*nrcv)
!> on all nodes using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgather_z(COMM, nrcv, x)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      COMPLEX(q) x(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF( nrcv > 0 ) CALL nccl_allgather_d(x(nrcv*(COMM%NODE_ME-1)+1), x(1), nrcv*2, COMM, ACC_ASYNC_Q)
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x)
            CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_double_complex, &
                             x(1), nrcv, MPI_double_complex, &
                             COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
            CHECK_MPI_ERROR(ierror,'M_allgather_z','MPI_allgather')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_double_complex, &
                         x(1), nrcv, MPI_double_complex, &
                         COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgather_z','MPI_allgather')

      END SUBROUTINE M_allgather_z


!----------------------------------------------------------------------
!
!> copy nrcv double precision from node i=1,COMM%NCPU to x((i-1)*nrcv+1:i*nrcv)
!> on all nodes using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgather_d(COMM, nrcv, x)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      REAL(q) x(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_double_precision, &
                         x(1), nrcv, MPI_double_precision, &
                         COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgather_d','MPI_allgather')

      END SUBROUTINE M_allgather_d

      SUBROUTINE M_allgather_i(COMM, nrcv, x)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      INTEGER x(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_integer, &
                         x(1), nrcv, MPI_integer, &
                         COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgather_i','MPI_allgather')

      END SUBROUTINE M_allgather_i

      SUBROUTINE M_allgather_i8(COMM, nrcv, x)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv

      INTEGER(qi8) x(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_allgather(MPI_IN_PLACE, 0, MPI_integer8, &
                         x(1), nrcv, MPI_integer, &
                         COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgather_i8','MPI_allgather')

      END SUBROUTINE M_allgather_i8


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgather_z
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgather_z(COMM, nrcv, x, request)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv
      INTEGER request

      COMPLEX(q) x(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1 ) RETURN

      CALL MPI_iallgather(MPI_IN_PLACE, 0, MPI_double_complex, &
                          x(1), nrcv, MPI_double_complex, &
                          COMM%MPI_COMM, request, ierror)
      CHECK_MPI_ERROR(ierror,'M_iallgather_z','MPI_iallgather')

      END SUBROUTINE M_iallgather_z


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgather_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgather_d(COMM, nrcv, x, request)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv
      INTEGER request

      REAL(q) x(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      CALL MPI_iallgather(MPI_IN_PLACE, 0, MPI_double_precision, &
                          x(1), nrcv, MPI_double_precision, &
                          COMM%MPI_COMM, request, ierror)
      CHECK_MPI_ERROR(ierror,'M_iallgather_d','MPI_iallgather')

      END SUBROUTINE M_iallgather_d


!----------------------------------------------------------------------
!
!> non-blocking version of #M_gather_z , gathers x from ranks in comm
!> into y
!
!----------------------------------------------------------------------

      SUBROUTINE M_gathero_z(COMM, n, x, y)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      COMPLEX(q) x(*),y(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_gathero_z')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.IS_ON_DEVICE(y).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_allgather_d(x(1), y(1), n*2, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_gathero_z')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x,y)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
         CALL MPI_gather( x(1), n, MPI_double_complex, &
                          y(1), n, MPI_double_complex, &
                          COMM%IONODE-1, COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_gathero_z')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_gather( x(1), n, MPI_double_complex, &
                       y(1), n, MPI_double_complex, &
                       COMM%IONODE-1, COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_gathero_z','MPI_gather')

      PROFILING_STOP('m_gathero_z')

      END SUBROUTINE M_gathero_z


!----------------------------------------------------------------------
!
!> non-blocking version of #M_gather_d , gathers x from ranks in comm
!> into y
!
!----------------------------------------------------------------------

      SUBROUTINE M_gathero_d(COMM, n, x, y)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      REAL(q) x(*),y(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_gathero_d')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.IS_ON_DEVICE(y).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF ( n > 0 ) CALL nccl_allgather_d(x(1), y(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_gathero_d')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x,y)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
            CALL MPI_gather( x(1), n, MPI_double_precision, &
                             y(1), n, MPI_double_precision, &
                             COMM%IONODE-1, COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_gathero_d')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_gather( x(1), n, MPI_double_precision, &
                       y(1), n, MPI_double_precision, &
                       COMM%IONODE-1, COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_gathero_d','MPI_gather')

      PROFILING_STOP('m_gathero_d')

      END SUBROUTINE M_gathero_d


!----------------------------------------------------------------------
!
!> copy n double complex x(1:n) from nodes i=1,COMM%NCPU to y((i-1)*n+1:i*n)
!> on all nodes using out-of-place communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgathero_z(COMM, n, x, y)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      COMPLEX(q) x(*),y(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_allgathero_z')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.IS_ON_DEVICE(y).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF (COMM%LUSENCCL) THEN
            IF ( n > 0 ) CALL nccl_allgather_d(x(1), y(1), n*2, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_allgathero_z')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x,y)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
            CALL MPI_allgather( x(1), n, MPI_double_complex, &
                                y(1), n, MPI_double_complex, &
                                COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_allgathero_z')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_allgather( x(1), n, MPI_double_complex, &
                          y(1), n, MPI_double_complex, &
                          COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgathero_z','MPI_allgather')

      PROFILING_STOP('m_allgathero_z')

      END SUBROUTINE M_allgathero_z


!----------------------------------------------------------------------
!
!> copy n double precision x(1:n) from nodes i=1,COMM%NCPU to y((i-1)*n+1:i*n)
!> on all nodes using out-of-place communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgathero_d(COMM, n, x, y)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n

      REAL(q) x(*),y(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_allgathero_d')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.IS_ON_DEVICE(y).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF (COMM%LUSENCCL) THEN
            IF ( n > 0 ) CALL nccl_allgather_d(x(1), y(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_allgathero_d')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x,y)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
         CALL MPI_allgather( x(1), n, MPI_double_precision, &
                             y(1), n, MPI_double_precision, &
                             COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_allgathero_d')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_allgather( x(1), n, MPI_double_precision, &
                          y(1), n, MPI_double_precision, &
                          COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgathero_d','MPI_allgather')

      PROFILING_STOP('m_allgathero_d')

      END SUBROUTINE M_allgathero_d


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgathero_z
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgathero_z(COMM, n, x, y, request)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      INTEGER request

      COMPLEX(q) x(*),y(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_iallgathero_z')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.IS_ON_DEVICE(y).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF (COMM%LUSENCCL) THEN
            IF ( n > 0 ) CALL nccl_allgather_d(x(1), y(1), n*2, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_iallgathero_z')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x,y)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
         CALL MPI_iallgather( x(1), n, MPI_double_complex, &
                              y(1), n, MPI_double_complex, &
                              COMM%MPI_COMM, request, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_iallgathero_z')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_iallgather( x(1), n, MPI_double_complex, &
                           y(1), n, MPI_double_complex, &
                           COMM%MPI_COMM, request, ierror)
      CHECK_MPI_ERROR(ierror,'M_iallgathero_z','MPI_iallgather')

      PROFILING_STOP('m_iallgathero_z')

      END SUBROUTINE M_iallgathero_z


!----------------------------------------------------------------------
!
!> non-blocking version of #M_allgathero_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_iallgathero_d(COMM, n, x, y, request)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      INTEGER request

      REAL(q) x(*),y(*)

      ! local variables
      INTEGER ierror

      IF (COMM%NCPU == 1) RETURN

      PROFILING_START('m_iallgathero_d')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(x).AND.IS_ON_DEVICE(y).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF (COMM%LUSENCCL) THEN
            IF ( n > 0 ) CALL nccl_allgather_d(x(1), y(1), n, COMM, ACC_ASYNC_Q)
            PROFILING_STOP('m_iallgathero_d')
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(x,y)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
         CALL MPI_iallgather( x(1), n, MPI_double_precision, &
                              y(1), n, MPI_double_precision, &
                              COMM%MPI_COMM, request, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         PROFILING_STOP('m_iallgathero_d')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_iallgather( x(1), n, MPI_double_precision, &
                           y(1), n, MPI_double_precision, &
                           COMM%MPI_COMM, request, ierror)
      CHECK_MPI_ERROR(ierror,'M_iallgathero_d','MPI_iallgather')

      PROFILING_STOP('m_iallgathero_d')

      END SUBROUTINE M_iallgathero_d


!----------------------------------------------------------------------
!
!> copy nrcv(i) n double complex from node i=1,COMM%NCPU to
!> c(prcv(i)+1:prcv(i)+nrcv(i)) on all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgatherv_z(COMM, x, c, nrcv, prcv)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv(COMM%NCPU)
      INTEGER prcv(COMM%NCPU)

      COMPLEX(q) x(*), c(*)

      ! local variables
      INTEGER i, r, ierror

      IF (COMM%NCPU == 1 ) THEN
#ifndef CRAYHIP
GPACC !$ACC KERNELS PRESENT(x,c) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:x(1:nrcv(1)),c(prcv(1)+1:prcv(1)+nrcv(1))) IF(OFFLOAD_ON)
         c(prcv(1)+1:prcv(1)+nrcv(1)) = x(1:nrcv(1))
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#else // CRAYHIP
         CALL WZCOPY(nrcv(1), x(1), 1, c(prcv(1)+1), 1)
#endif // CRAYHIP
         RETURN
      ENDIF

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(c).AND.IS_ON_DEVICE(x).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (COMM%LGRPNCCL) CALL nccl_group_start()
            DO r=1,COMM%NCPU
               IF(nrcv(COMM%NODE_ME)>0) CALL nccl_send_d(x(1), nrcv(COMM%NODE_ME)*2, r-1, COMM, ACC_ASYNC_Q)
               IF(nrcv(r)>0)            CALL nccl_recv_d(c(prcv(r)+1),    nrcv(r)*2, r-1, COMM, ACC_ASYNC_Q)
            ENDDO
            IF (COMM%LGRPNCCL) CALL nccl_group_end()
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(X,C)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(X,C)
         CALL MPI_allgatherv(x(1), nrcv(COMM%NODE_ME), MPI_double_complex, &
                             c(1), nrcv(1), prcv(1), MPI_double_complex, &
                             COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_allgatherv_z','MPI_allgatherv')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_allgatherv(x(1), nrcv(COMM%NODE_ME), MPI_double_complex, &
                          c(1), nrcv(1), prcv(1), MPI_double_complex, &
                          COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgatherv_z','MPI_allgatherv')

      END SUBROUTINE M_allgatherv_z


!----------------------------------------------------------------------
!
!> copy nrcv(i) n double precision from node i=1,COMM%NCPU to
!> c(prcv(i)+1:prcv(i)+nrcv(i)) on all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_allgatherv_d(COMM, x, c, nrcv, prcv)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER nrcv(COMM%NCPU)
      INTEGER prcv(COMM%NCPU)

      REAL(q) x(*), c(*)

      ! local variables
      INTEGER i, r, ierror

      IF (COMM%NCPU == 1 ) THEN
#ifndef CRAYHIP
GPACC !$ACC KERNELS PRESENT(x,c) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:x(1:nrcv(1)),c(prcv(1)+1:prcv(1)+nrcv(1))) IF(OFFLOAD_ON)
         c(prcv(1)+1:prcv(1)+nrcv(1)) = x(1:nrcv(1))
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#else // CRAYHIP
         CALL WDCOPY(nrcv(1), x(1), 1, c(prcv(1)+1), 1)
#endif // CRAYHIP
         RETURN
      ENDIF

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(c).AND.IS_ON_DEVICE(x).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (COMM%LGRPNCCL) CALL nccl_group_start()
            DO r=1,COMM%NCPU
               IF(nrcv(COMM%NODE_ME)>0) CALL nccl_send_d(x(1), nrcv(COMM%NODE_ME), r-1, COMM, ACC_ASYNC_Q)
               IF(nrcv(r)>0)            CALL nccl_recv_d(c(prcv(r)+1),    nrcv(r), r-1, COMM, ACC_ASYNC_Q)
            ENDDO
            IF (COMM%LGRPNCCL) CALL nccl_group_end()
            RETURN
         ENDIF
#endif // USENCCL

GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(X,C)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(X,C)
         CALL MPI_allgatherv(x(1), nrcv(COMM%NODE_ME), MPI_double_precision, &
                             c(1), nrcv(1), prcv(1), MPI_double_precision, &
                             COMM%MPI_COMM, ierror)
GPOMP !$OMP END TARGET DATA
GPACC !$ACC END HOST_DATA
         CHECK_MPI_ERROR(ierror,'M_allgatherv_d','MPI_allgatherv')
         RETURN
      ENDIF
#endif // OFFLOADING

      CALL MPI_allgatherv(x(1), nrcv(COMM%NODE_ME), MPI_double_precision, &
                           c(1), nrcv(1), prcv(1), MPI_double_precision, &
                           COMM%MPI_COMM, ierror)
      CHECK_MPI_ERROR(ierror,'M_allgatherv_d','MPI_allgatherv')

      END SUBROUTINE M_allgatherv_d


!----------------------------------------------------------------------
!
!> reduce (MPI_sum) n complex doubles to node i from all nodes to node i
!> using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_reduce_z_to(COMM, vec, n, inode)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER inode

      ! local variables
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      PROFILING_START('m_reduce_z_to')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (n>0) CALL nccl_reduce_d(vec(1), vec(1), n*2, inode-1, COMM, ACC_ASYNC_Q)
         ELSE
#endif
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
            IF (COMM%NODE_ME==inode) THEN
               CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
            ELSE
               CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_complex, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
            ENDIF
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
            CHECK_MPI_ERROR(ierror,'M_reduce_z_to','MPI_reduce')
#ifdef USENCCL
         ENDIF
#endif
         PROFILING_STOP('m_reduce_z_to')
         RETURN
      ENDIF
#endif

      IF (COMM%NODE_ME==inode) THEN
         CALL MPI_reduce( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ELSE
         CALL MPI_reduce( vec(1),       vec(1), n, MPI_double_complex, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ENDIF

      CHECK_MPI_ERROR(ierror,'M_reduce_z_to','MPI_reduce')

      PROFILING_STOP('m_reduce_z_to')

      END SUBROUTINE M_reduce_z_to


!----------------------------------------------------------------------
!
!> reduce (MPI_sum) n doubles to node i from all nodes to node i using
!> inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_reduce_d_to(COMM, vec, n, inode)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER inode

      ! local variables
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      PROFILING_START('m_reduce_d_to')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (n>0) CALL nccl_reduce_d(vec(1), vec(1), n, inode-1, COMM, ACC_ASYNC_Q)
         ELSE
#endif
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
            IF (COMM%NODE_ME==inode) THEN
               CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
            ELSE
               CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_precision, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
            ENDIF
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
            CHECK_MPI_ERROR(ierror,'M_reduce_d_to','MPI_reduce')
#ifdef USENCCL
         ENDIF
#endif
         PROFILING_STOP('m_reduce_d_to')
         RETURN
      ENDIF
#endif

      IF (COMM%NODE_ME==inode) THEN
         CALL MPI_reduce( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ELSE
         CALL MPI_reduce( vec(1),       vec(1), n, MPI_double_precision, &
                          MPI_sum, inode-1, COMM%MPI_COMM, ierror)
      ENDIF

      CHECK_MPI_ERROR(ierror,'M_reduce_d_to','MPI_reduce')

      PROFILING_STOP('m_reduce_d_to')

      END SUBROUTINE M_reduce_d_to


!----------------------------------------------------------------------
!
!> non-blocking reduce (MPI_sum) n complex doubles from all nodes to node i
!> using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_ireduce_z_to(COMM, vec, n, inode, request)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER inode
      INTEGER request

      ! local variables
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      PROFILING_START('m_ireduce_z_to')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (n>0) CALL nccl_reduce_d(vec(1), vec(1), n*2, inode-1, COMM, ACC_ASYNC_Q)
         ELSE
#endif
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
            IF (COMM%NODE_ME==inode) THEN
#if !defined(no_cuda_aware_ireduce) && !defined(MPI_avoid_ireduce_to)
               CALL MPI_ireduce( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                                 MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
               CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
            ELSE
#if !defined(no_cuda_aware_ireduce) && !defined(MPI_avoid_ireduce_to)
               CALL MPI_ireduce( vec(1),       vec(1), n, MPI_double_complex, &
                                 MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
               CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_complex, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
            ENDIF
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
            CHECK_MPI_ERROR(ierror,'M_ireduce_z_to','MPI_ireduce')
#if defined(no_cuda_aware_ireduce) || defined(MPI_avoid_ireduce_to)
            request = MPI_REQUEST_NULL
#endif
#ifdef USENCCL
         ENDIF
#endif
         PROFILING_STOP('m_ireduce_z_to')
         RETURN
      ENDIF
#endif

      IF (COMM%NODE_ME==inode) THEN
#ifndef MPI_avoid_ireduce_to
         CALL MPI_ireduce( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                           MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
         CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_complex, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
      ELSE
#ifndef MPI_avoid_ireduce_to
         CALL MPI_ireduce( vec(1),       vec(1), n, MPI_double_complex, &
                           MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
         CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_complex, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
      ENDIF

      CHECK_MPI_ERROR(ierror,'M_ireduce_z_to','MPI_ireduce')
#ifdef MPI_avoid_ireduce_to
      request = MPI_REQUEST_NULL
#endif
      PROFILING_STOP('m_ireduce_z_to')

      END SUBROUTINE M_ireduce_z_to


!----------------------------------------------------------------------
!
!> non-blocking reduce (MPI_sum) n doubles to node i from all nodes to node i
!> using inplace communication
!
!----------------------------------------------------------------------

      SUBROUTINE M_ireduce_d_to(COMM, vec, n, inode, request)

      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER inode
      INTEGER request

      ! local variables
      INTEGER ierror

! quick return if possible
      IF (COMM%NCPU == 1 ) RETURN

      PROFILING_START('m_ireduce_d_to')

#ifdef OFFLOADING
      IF (IS_ON_DEVICE(vec).AND.OFFLOAD_ON) THEN
#ifdef USENCCL
         IF ( COMM%LUSENCCL ) THEN
            IF (n>0) CALL nccl_reduce_d(vec(1), vec(1), n, inode-1, COMM, ACC_ASYNC_Q)
         ELSE
#endif
GPACC !$ACC WAIT(ACC_ASYNC_Q)
GPACC !$ACC HOST_DATA USE_DEVICE(vec)
GPOMP !$OMP TARGET DATA USE_DEVICE_ADDR(vec)
            IF (COMM%NODE_ME==inode) THEN
#if !defined(no_cuda_aware_ireduce) && !defined(MPI_avoid_ireduce_to)
               CALL MPI_ireduce( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                                 MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
               CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
            ELSE
#if !defined(no_cuda_aware_ireduce) && !defined(MPI_avoid_ireduce_to)
               CALL MPI_ireduce( vec(1),       vec(1), n, MPI_double_precision, &
                                 MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
               CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_precision, &
                                 MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
            ENDIF
GPACC !$ACC END HOST_DATA
GPOMP !$OMP END TARGET DATA
            CHECK_MPI_ERROR(ierror,'M_ireduce_d_to','MPI_ireduce')
#if defined(no_cuda_aware_ireduce) || defined(MPI_avoid_ireduce_to)
            request = MPI_REQUEST_NULL
#endif
#ifdef USENCCL
         ENDIF
#endif
         PROFILING_STOP('m_ireduce_d_to')
         RETURN
      ENDIF
#endif

      IF (COMM%NODE_ME==inode) THEN
#ifndef MPI_avoid_ireduce_to
         CALL MPI_ireduce( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                           MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
         CALL MPI_reduce ( MPI_IN_PLACE, vec(1), n, MPI_double_precision, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
      ELSE
#ifndef MPI_avoid_ireduce_to
         CALL MPI_ireduce( vec(1),       vec(1), n, MPI_double_precision, &
                           MPI_sum, inode-1, COMM%MPI_COMM, request, ierror)
#else
         CALL MPI_reduce ( vec(1),       vec(1), n, MPI_double_precision, &
                           MPI_sum, inode-1, COMM%MPI_COMM,          ierror)
#endif
      ENDIF

      CHECK_MPI_ERROR(ierror,'M_ireduce_d_to','MPI_ireduce')
#ifdef MPI_avoid_ireduce_to
      request = MPI_REQUEST_NULL
#endif
      PROFILING_STOP('m_ireduce_d_to')

      END SUBROUTINE M_ireduce_d_to


!----------------------------------------------------------------------
!
!> just a wrapper of MPI_waitall
!
!----------------------------------------------------------------------

      SUBROUTINE M_waitall(n,requests)
      USE mpimy
      IMPLICIT NONE
      INTEGER n,requests(n)
      ! local variables
      INTEGER ierror

      CALL MPI_waitall(n, requests, MPI_STATUSES_IGNORE, ierror)

      CHECK_MPI_ERROR(ierror,'M_waitall','MPI_waitall')

      END SUBROUTINE M_waitall


!----------------------------------------------------------------------
!
!> Subroutine to compare integer values among ranks within a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_compare_i(COMM, TARGET_VALUE, LSAME)
         USE mpimy
         IMPLICIT NONE
         TYPE(communic), INTENT(IN)     :: COMM
         INTEGER, INTENT(IN)            :: TARGET_VALUE
         LOGICAL, INTENT(OUT)           :: LSAME
         ! Local variables
         INTEGER :: ARRAY_VALUE(1:COMM%NCPU)
         INTEGER :: ICORE

         ! Initialize arrays holding values and write into position of process number
         ARRAY_VALUE = 0
         ARRAY_VALUE(COMM%NODE_ME) = TARGET_VALUE

         ! After this each core holds the value for each core in the array
         CALLMPI ( M_sum_i(COMM, ARRAY_VALUE(1), COMM%NCPU) )

         ! Initialize break condition
         LSAME = .TRUE.

         ! Everybody calculates I_THE_SAME, but only one core
         ! broadcasts at the end to make sure every process is synced
         DO ICORE = 1, COMM%NCPU
            IF (ARRAY_VALUE(ICORE) .NE. ARRAY_VALUE(1)) THEN
               LSAME = .FALSE.
            ENDIF
         ENDDO
         CALLMPI( M_bcast_l( COMM, LSAME, 1))
      END SUBROUTINE M_compare_i

!----------------------------------------------------------------------
!
!> Subroutine to compare integer values among ranks within a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_compare_d(COMM, TARGET_VALUE, TOLERANCE, LSAME)
         USE mpimy
         IMPLICIT NONE
         TYPE(communic), INTENT(IN)     :: COMM
         REAL(q), INTENT(IN)            :: TARGET_VALUE
         REAL(q), INTENT(IN)            :: TOLERANCE
         LOGICAL, INTENT(OUT)           :: LSAME
         ! Local variables
         REAL(q) :: ARRAY_VALUE(1:COMM%NCPU)
         INTEGER :: ICORE

         ! Initialize arrays holding values and write into position of process number
         ARRAY_VALUE = 0.0_q
         ARRAY_VALUE(COMM%NODE_ME) = TARGET_VALUE

         ! After this each core holds the value for each core in the array
         CALLMPI ( M_sum_d(COMM, ARRAY_VALUE(1), COMM%NCPU) )

         ! Initialize boolean
         LSAME = .TRUE.

         ! Everybody calculates I_THE_SAME, but only one core
         ! broadcasts at the end to make sure every process is synced
         DO ICORE = 1, COMM%NCPU
            IF (ABS(ARRAY_VALUE(ICORE)-ARRAY_VALUE(1)) .GT. TOLERANCE) THEN
               LSAME = .FALSE.
            ENDIF
         ENDDO
         CALLMPI( M_bcast_l( COMM, LSAME, 1))
      END SUBROUTINE M_compare_d

#else // defined(MPI) || defined(MPI_CHAIN)

!======================================================================
! RCS:  $Id: mpi.F,v 1.6 2003/06/27 13:22:20 kresse Exp kresse $
!
!> dummy module if MPI is not used
!>
!> a few files will not compile with this dummy module
!> i.e. fftmpi.F fftmpi_map.F
!
!======================================================================
      MODULE mpimy
      !> dummy communicator for serial execution
      TYPE communic
        INTEGER :: NODE_ME = 1
        INTEGER :: IONODE = 1
        INTEGER :: NCPU = 1
      END TYPE
      CONTAINS

!*******************************************************************
!>computes local index of a cyclic distribution of an array
!>distributed among NCPU ranks
!>   @param[in]     NCPU number of MPI ranks
!>   @param[in]     NODE_ME rank id
!>   @param[out]    M index position of rank NODE_ME
!*******************************************************************
   FUNCTION LOCAL_INDEX_CYCLIC( NCPU, NODE_ME, M )
      IMPLICIT NONE
      INTEGER  :: NODE_ME
      INTEGER  :: M
      INTEGER  :: NCPU
      INTEGER  :: LOCAL_INDEX_CYCLIC
      LOCAL_INDEX_CYCLIC = M
   END FUNCTION LOCAL_INDEX_CYCLIC

      SUBROUTINE mpi_dummy
      WRITE(*,*)'Im a DEC compiler so I need this line'
      END SUBROUTINE
      END MODULE

#endif // defined(MPI) || defined(MPI_CHAIN)
