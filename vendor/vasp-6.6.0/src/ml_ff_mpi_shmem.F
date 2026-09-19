#include "symbol.inc"
#ifdef use_shmem
!**********************************************************************
!> @brief
!! This module provides the means of allocating and deallocating
!! shared memory segments by means of the MPI-3 shared memory windows
!! functionality.
!! Furthermore this module provides semaphore functionality based on an
!! MPI-3 remote memory access (RMA) implementation of the (MCS) locks
!! of: J. M. Mellor-Crummey, M. L. Scott, ACM Trans. Comput. Syst. 9,
!! 21-65, 1991.
!!
!! Shared memory:
!!
!!  mpi_shmem::m_shmem
!!  mpi_shmem::m_shmem_alloc
!!  mpi_shmem::m_shmem_dealloc
!!
!! Semaphores (aka locks):
!!
!!  mpi_shmem::m_shmem_lock
!!  mpi_shmem::m_getsem
!!  mpi_shmem::m_locksem
!!  mpi_shmem::m_unlocksem
!!  mpi_shmem::m_destroysem
!!
!! @note
!!  To use the functionality in this module VASP needs to be
!!  compiled with:
!!  <tt>-Duse_shmem</tt>,
!!  <tt>-Dshmem_rproj</tt>, or
!!  <tt>-Dshmem_bcast_buffer</tt>
!!
!! @note
!!  When compiled with <tt>-Dsysv</tt> the routines in this module
!!  will fall back to the use of SYSV shared memory segments and
!!  semaphores, instead of using the MPI-3 shmem functionality and
!!  MCS locks.
!!
!**********************************************************************
      MODULE ml_ff_mpi_shmem

      USE MPI_HELP, ONLY: M_BCAST
      USE ML_FF_PREC
      USE ML_FF_STRUCT, ONLY: ML_MPI_PAR, M_shmem

#ifndef sysv
#define SIZE_KIND MPI_ADDRESS_KIND
#define GET_SIZE(x) CALL MPI_Sizeof(x,sizeoftype,ierror)

      USE iso_c_binding

      USE mpi_f08, ONLY : MPI_ADDRESS_KIND, &
                          MPI_INFO_NULL, &
                          MPI_MODE_NOSUCCEED, &
                          MPI_KEYVAL_INVALID, &
                          MPI_WIN_NULL_COPY_FN, &
                          MPI_WIN_NULL_DELETE_FN, &
                          MPI_WIN_BASE, &
                          MPI_REPLACE, &
                          MPI_INTEGER, &
                          MPI_Win, &
                          MPI_Comm, &
                          MPI_Sizeof, &
                          MPI_Win_allocate_shared, &
                          MPI_Win_shared_query, &
                          MPI_Win_fence, &
                          MPI_Win_free, &
                          MPI_Win_allocate, &
                          MPI_Win_lock_all, &
                          MPI_Win_unlock_all, &
                          MPI_Win_flush, &
                          MPI_Fetch_and_op, &
                          MPI_Accumulate, &
                          MPI_Compare_and_swap


#else
#define SIZE_KIND c_size_t
#define GET_SIZE(x) sizeoftype=SIZEOF(x)

      USE ml_ff_c2f_interface
#endif

      IMPLICIT NONE

      PUBLIC :: M_shmem,M_shmem_alloc,M_shmem_dealloc, &
                M_shmem_lock,M_getsem,M_locksem, &
                M_unlocksem,M_destroysem

#ifdef sysv
      PUBLIC :: M_shmem_read_shmmax,M_shmem_check_max_allocation_vs_shmmax
#endif
      PRIVATE

!****************** SUBROUTINE M_shmem_alloc **************************
!> @details
!! mpi_shmem::m_shmem_alloc allocates a 1d, 2d, or 3d, REAL(q) or
!! COMPLEX(q) shared memory segment shared by the MPI ranks in comm
!! and returns pointers to it.
!!
!! @param[in]  comm    Communicator
!!                       (mpimy::ML_MPI_PAR)
!! @param[in]  dim1    The first dimension of the shmem segment
!!                       (integer)
!! @param[in]  dim2    The second dimension of the shmem segment
!!                       (integer, optional)
!! @param[in]  dim3    The third dimension of the shmem segment
!!                       (integer, optional)
!! @param[out] segment Shared memory object handle
!!                       (mpi_shmem::m_shmem)
!! @param[out] ptr     Pointer to shared memory segment
!!                       (1d, 2d, or 3d pointer of type REAL(q)
!!                        or COMPLEX(q))
!! @param[out] istat   Exit status, istat=0 on succes 
!!                       (integer, optional)
!! @remark
!!  Actual allocation is only done by MPI rank 0
!!
!! @warning
!!  All MPI ranks in comm must call this routine
!! @warning
!!  The MPI ranks in comm must share direct access to the same
!!  physical memory (i.e., they must reside on the same physical
!!  compute node).
!!
!**********************************************************************

      INTERFACE M_shmem_alloc
         MODULE PROCEDURE shmem_alloc_c_1d
         MODULE PROCEDURE shmem_alloc_c_2d
         MODULE PROCEDURE shmem_alloc_c_3d
         MODULE PROCEDURE shmem_alloc_r_1d
         MODULE PROCEDURE shmem_alloc_r_2d
         MODULE PROCEDURE shmem_alloc_r_3d
      END INTERFACE

#ifndef sysv
      TYPE MCS_lock
         TYPE(MPI_Win)           :: lock_win
         TYPE(ML_MPI_PAR), POINTER :: comm => NULL()
         INTEGER, POINTER        :: lmem(:)
         INTEGER                 :: iownlock
      END TYPE MCS_lock
#endif

!****************** SUBROUTINE M_shmem_lock ***************************
!> @brief
!! Lock object handle
!!
!! The lock object handle is basically an array of mpi_shmem::mcs_lock
!!
!**********************************************************************

      TYPE M_shmem_lock
#ifndef sysv
         TYPE(MCS_lock), POINTER :: lock(:) => NULL()
         INTEGER                 :: ilock0
#else
         INTEGER(KIND=c_int)     :: lock
         TYPE(ML_MPI_PAR), POINTER :: comm => NULL()
#endif
      END TYPE M_shmem_lock

#ifndef sysv
      ! MCS lock related internal variables
      INTEGER :: nextRank=1, blocked=2, lockTail=3
#else
      CHARACTER(LEN=200)        :: FILE_READ_SHMMAX
      INTEGER                   :: IERR_READ_SHMMAX
      INTEGER                   :: IERR_OPEN_SHMMAX
      INTEGER(KIND=c_size_t)    :: AVAILABLE_MEM
      CHARACTER(LEN=200)        :: AVAILABLE_MEM_STRING
      LOGICAL                   :: OVER18
#endif
      CONTAINS

#ifdef sysv
!****************** SUBROUTINE M_shmem_read_shmmax ********************
!
! This routine opens the file where the maximum size 
! of shared memory segments is stored. 
! The file in most systems is at /proc/sys/kernel/shmmax. 
! If the file cannot be opened due to differences in the operation 
! system, the code gives a warning that the file could not be read,
! but the code continues without checks of the shared memory. 
!
!**********************************************************************

      SUBROUTINE M_shmem_read_shmmax(comm_world,comm)
         USE tutor
         TYPE (ML_MPI_PAR), TARGET :: comm_world
         TYPE (ML_MPI_PAR)         :: comm
         ! Local variables
         INTEGER                 :: FUNIT
         INTEGER                 :: IERR_READ_TMP
         INTEGER                 :: myid, ierror
         INCLUDE "mpif.h"

         myid=comm%node_me+1
         
         FUNIT = 99999 
         FILE_READ_SHMMAX = "/proc/sys/kernel/shmmax"


         OVER18 = .FALSE.
         IERR_READ_SHMMAX = 0
         IERR_OPEN_SHMMAX = 0
         IERR_READ_TMP = 0
         IF (myid==1) THEN
            OPEN(UNIT=FUNIT,FILE=TRIM(ADJUSTL(FILE_READ_SHMMAX)),STATUS='OLD',IOSTAT=IERR_OPEN_SHMMAX)

            ! First we read into a string because we have no clue how big this number is.
            ! On some computer this number is set to the highest number of an unsigned
            ! 8 byte integer, which is not avaiable on all Fortran machines (Fortran standard
            ! requires signed 8 byte integer). 
            IF (IERR_OPEN_SHMMAX .EQ. 0) THEN
               READ(FUNIT,*,IOSTAT=IERR_READ_SHMMAX) AVAILABLE_MEM_STRING
               IERR_READ_TMP = IERR_READ_SHMMAX

               ! So here we need to check if size of available memory has more than 18 digits
               ! because 8 bit signed integer (2*63) has 19 digits. With 19 digits the size could
               ! already above 8 bit signed integer. That is anyway quasi infinite amount 
               ! of segment size for shared memory. If it's below 19 digits we need to read it
               ! into a signed integer array and check if we have enough.
               IF (LEN(TRIM(ADJUSTL(AVAILABLE_MEM_STRING))) .GT. 18) THEN
                  OVER18 = .TRUE.
               ELSE
                  REWIND(99999)
                  READ(99999,*,IOSTAT=IERR_READ_SHMMAX) AVAILABLE_MEM
                  IERR_READ_SHMMAX = IERR_READ_SHMMAX + IERR_READ_TMP
               ENDIF
            ENDIF
            CLOSE(99999)
         ENDIF

         ! If desired file could not be open then no memory checking can be done
         ! and an warning message is given. Nevertheless we need to continue
         ! because the location of the maximum shared memory segment size
         ! can be system dependent.
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, IERR_OPEN_SHMMAX, &
                            1, MPI_INTEGER, &
                            MPI_SUM, comm_world, ierror)
         IF (IERR_OPEN_SHMMAX .GT. 0 ) THEN
            CALL vtutor%alert("Could not open: " //TRIM(ADJUSTL(FILE_READ_SHMMAX))// ". &
                 &Checking maximum system V shmem segment size is disabled. & 
                 &Continuing with the calculation.")
         ENDIF

         ! If the system V file exists but there was a problem with 
         ! reading it, the code woul give an error and exit.
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, IERR_READ_SHMMAX, &
                            1, MPI_INTEGER, &
                            MPI_SUM, comm_world, ierror)
         IF (IERR_READ_SHMMAX .GT. 0 ) THEN
            CALL vtutor%alert("Could not read: " //TRIM(ADJUSTL(FILE_READ_SHMMAX))// ". &
                 &Checking maximum system V shmem segment size is disabled. & 
                 &Continuing with the calculation.")
         ENDIF
      END SUBROUTINE M_shmem_read_shmmax

!****************** M_shmem_check_max_allocation_vs_shmmax*************
!
! The input value max_segmentsize is the maximal segmentsize that
! is going to be allocated later. If that size is larger than the
! maximally available shared memory segment (AVAILABLE_MEM) than
! the code exits with an error. Otherwise codes continues without
! any additional messages.
!
!**********************************************************************
      SUBROUTINE M_shmem_check_max_allocation_vs_shmmax(comm_world, comm, max_segmentsize)
         USE tutor
         USE string, only: str
         TYPE (ML_MPI_PAR), TARGET :: comm_world
         TYPE (ML_MPI_PAR)         :: comm
         INTEGER(kind=8), INTENT(IN)       :: max_segmentsize
         ! Local variables
         INTEGER                   :: MEMORY_CHECK
         INTEGER                   :: myid, ierror
         INTEGER(KIND=c_size_t)    :: segmenttypesize 
         INTEGER(KIND=c_size_t)    :: PRODUCT_HELP
         REAL(q)                   :: rdummy
         INTEGER                   :: sizeoftype
         
         INCLUDE "mpif.h"

         myid=comm%node_me+1

         MEMORY_CHECK = 0

         ! Get size of 
         GET_SIZE(rdummy)
         segmenttypesize = INT(sizeoftype,KIND=c_size_t)
         PRODUCT_HELP = segmenttypesize * INT(max_segmentsize,KIND=c_size_t)

         ! Do checking of available memory if file containing
         ! shared memory segments could be read succesfully.
         IF (IERR_OPEN_SHMMAX.EQ.0 .AND. IERR_READ_SHMMAX.EQ.0) THEN
            IF (myid==1) THEN
               IF (OVER18) THEN
                  MEMORY_CHECK = 0
               ELSE
                  IF (AVAILABLE_MEM .LT. PRODUCT_HELP) THEN
                     MEMORY_CHECK = 1
                  ELSE
                     MEMORY_CHECK = 0
                  ENDIF
               ENDIF
            ENDIF

            CALL MPI_ALLREDUCE(MPI_IN_PLACE, MEMORY_CHECK,  &
                               1, MPI_INTEGER, &
                               MPI_SUM, comm_world, ierror)

            IF (MEMORY_CHECK .GT. 0 ) THEN
               CALL vtutor%error("Maximum system V shmem segment size &
                  &(see " //TRIM(ADJUSTL(FILE_READ_SHMMAX))// ") &
                  &is too small: " //str(AVAILABLE_MEM)// ". &
                  &The maximum required memory is "&
                  &// str(PRODUCT_HELP) // ". Please increase &
                  &the size or switch to shared memory without &
                  &system V (compile without -Dsysv).")
            ENDIF
         ENDIF
      END SUBROUTINE M_shmem_check_max_allocation_vs_shmmax
#endif


!****************** SUBROUTINE M_shmem_alloc_c_1d *********************
!
! Allocates a 1d complex(q) memory array of size (dim1)
! shared by the MPI ranks in comm and returns a pointer to it:
!
!    segment%ftr(1:dim1)
!
! Actual allocation is only done by MPI rank 0
!
! N.B.: All MPI ranks in comm must call this routine
!
! N.B.II: the MPI ranks in comm must share direct access to the same
!         physical memory (i.e., they must reside on the same physical
!         compute node).
!
!**********************************************************************

      SUBROUTINE shmem_alloc_c_1d(comm,dim1_help,segment,ptr,istat)

      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem)          :: segment
      COMPLEX(q), POINTER     :: ptr(:)
      INTEGER                 :: dim1_help
      INTEGER, OPTIONAL       :: istat

      ! local variables
      INTEGER(KIND=SIZE_KIND) :: dim1
      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype

      INTEGER                 :: ierror

      COMPLEX(q)              :: cdummy

      dim1 = INT(dim1_help,KIND=SIZE_KIND)

      CALL M_barrier(comm)

      IF (PRESENT(istat)) istat=0

      segment%desc%comm => comm

      IF (segment%desc%comm%ncpu>1) THEN
         GET_SIZE(cdummy)
         segmentsize=sizeoftype*MAX(dim1,1)

         CALL M_getshmem(segment%desc%comm, segmentsize, sizeoftype, &
                         segment%desc%shmem_win, segment%desc%shmem_ptr, ierror )
      ELSE
         ALLOCATE(segment%fptr_c(1:dim1),STAT=ierror)
      ENDIF

      IF (PRESENT(istat)) THEN
         IF (ierror/=0) THEN
            istat=ierror; RETURN
         ENDIF
      ENDIF

      IF (segment%desc%comm%ncpu>1) CALL c_f_pointer(segment%desc%shmem_ptr,segment%fptr_c,[dim1])

      ptr(1:dim1) => segment%fptr_c(1:dim1)

      segment%desc%dim1=dim1

      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE shmem_alloc_c_1d


!****************** SUBROUTINE M_shmem_alloc_c_2d *********************
!
! Allocates a 2d complex(q) memory array of size (dim1,dim2)
! shared by the MPI ranks in comm and returns a pointer to it:
!
!    segment%ftr(1:dim1,1:dim2)
!
! Actual allocation is only done by MPI rank 0
!
! N.B.: All MPI ranks in comm must call this routine
!
! N.B.II: the MPI ranks in comm must share direct access to the same
!         physical memory (i.e., they must reside on the same physical
!         compute node).
!
!**********************************************************************

      SUBROUTINE shmem_alloc_c_2d(comm,dim1_help,dim2_help,segment,ptr,istat)

      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem)          :: segment
      COMPLEX(q), POINTER     :: ptr(:,:)
      INTEGER                 :: dim1_help,dim2_help
      INTEGER, OPTIONAL       :: istat

      ! local variables
      INTEGER(KIND=SIZE_KIND) :: dim1,dim2
      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype

      INTEGER                 :: ierror

      COMPLEX(q)              :: cdummy

      dim1 = INT(dim1_help,KIND=SIZE_KIND)
      dim2 = INT(dim2_help,KIND=SIZE_KIND)

      CALL M_barrier(comm)

      IF (PRESENT(istat)) istat=0

      segment%desc%comm => comm

      IF (segment%desc%comm%ncpu>1) THEN
         GET_SIZE(cdummy)
         segmentsize=sizeoftype*MAX(dim1*dim2,1)

         CALL M_getshmem(segment%desc%comm, segmentsize, sizeoftype, &
                         segment%desc%shmem_win, segment%desc%shmem_ptr, ierror )
      ELSE
         ALLOCATE(segment%fptr_c(1:dim1*dim2),STAT=ierror)
      ENDIF

      IF (PRESENT(istat)) THEN
         IF (ierror/=0) THEN
            istat=ierror; RETURN
         ENDIF
      ENDIF

      IF (segment%desc%comm%ncpu>1) CALL c_f_pointer(segment%desc%shmem_ptr,segment%fptr_c,[dim1*dim2])

      ptr(1:dim1,1:dim2) => segment%fptr_c(1:dim1*dim2)

      segment%desc%dim1=dim1; segment%desc%dim2=dim2

      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE shmem_alloc_c_2d


!****************** SUBROUTINE M_shmem_alloc_c_3d *********************
!
! Allocates a 3d complex(q) memory array of size (dim1,dim2,dim3)
! shared by the MPI ranks in comm and returns a pointer to it:
!
!    segment%ftr(1:dim1,1:dim2,1:dim3)
!
! Actual allocation is only done by MPI rank 0
!
! N.B.: All MPI ranks in comm must call this routine
!
! N.B.II: the MPI ranks in comm must share direct access to the same
!         physical memory (i.e., they must reside on the same physical
!         compute node).
!
!**********************************************************************

      SUBROUTINE shmem_alloc_c_3d(comm,dim1_help,dim2_help,dim3_help,segment,ptr,istat)

      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem)          :: segment
      COMPLEX(q), POINTER     :: ptr(:,:,:)
      INTEGER                 :: dim1_help,dim2_help,dim3_help
      INTEGER, OPTIONAL       :: istat

      ! local variables
      INTEGER(KIND=SIZE_KIND) :: dim1,dim2,dim3
      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype

      INTEGER                 :: ierror

      COMPLEX(q)              :: cdummy

      dim1 = INT(dim1_help,KIND=SIZE_KIND)
      dim2 = INT(dim2_help,KIND=SIZE_KIND)
      dim3 = INT(dim3_help,KIND=SIZE_KIND)

      CALL M_barrier(comm)

      IF (PRESENT(istat)) istat=0

      segment%desc%comm => comm

      IF (segment%desc%comm%ncpu>1) THEN
         GET_SIZE(cdummy)
         segmentsize=sizeoftype*MAX(dim1*dim2*dim3,1)

         CALL M_getshmem(segment%desc%comm, segmentsize, sizeoftype, &
                         segment%desc%shmem_win, segment%desc%shmem_ptr, ierror )
      ELSE
         ALLOCATE(segment%fptr_c(1:dim1*dim2*dim3),STAT=ierror)
      ENDIF

      IF (PRESENT(istat)) THEN
         IF (ierror/=0) THEN
            istat=ierror; RETURN
         ENDIF
      ENDIF

      IF (segment%desc%comm%ncpu>1) CALL c_f_pointer(segment%desc%shmem_ptr,segment%fptr_c,[dim1*dim2*dim3])

      ptr(1:dim1,1:dim2,1:dim3) => segment%fptr_c(1:dim1*dim2*dim3)

      segment%desc%dim1=dim1; segment%desc%dim2=dim2; segment%desc%dim3=dim3

      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE shmem_alloc_c_3d


!****************** SUBROUTINE M_shmem_alloc_r_1d *********************
!
! Allocates a 1d real(q) memory array of size (dim1)
! shared by the MPI ranks in comm and returns a pointer to it:
!
!    segment%ftr(1:dim1)
!
! Actual allocation is only done by MPI rank 0
!
! N.B.: All MPI ranks in comm must call this routine
!
! N.B.II: the MPI ranks in comm must share direct access to the same
!         physical memory (i.e., they must reside on the same physical
!         compute node).
!
!**********************************************************************

      SUBROUTINE shmem_alloc_r_1d(comm,dim1_help,segment,ptr,istat)

      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem)          :: segment
      REAL(q), POINTER        :: ptr(:)
      INTEGER                 :: dim1_help
      INTEGER, OPTIONAL       :: istat

      ! local variables
      INTEGER(KIND=SIZE_KIND) :: dim1
      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype

      INTEGER                 :: ierror

      REAL(q)                 :: rdummy

      dim1 = INT(dim1_help,KIND=SIZE_KIND)

      CALL M_barrier(comm)

      IF (PRESENT(istat)) istat=0

      segment%desc%comm => comm

      IF (segment%desc%comm%ncpu>1) THEN
         GET_SIZE(rdummy)
         segmentsize=sizeoftype*MAX(dim1,1)

         CALL M_getshmem(segment%desc%comm, segmentsize, sizeoftype, &
                         segment%desc%shmem_win, segment%desc%shmem_ptr, ierror )
      ELSE
         ALLOCATE(segment%fptr_r(1:dim1),STAT=ierror)
      ENDIF

      IF (PRESENT(istat)) THEN
         IF (ierror/=0) THEN
            istat=ierror; RETURN
         ENDIF
      ENDIF

      IF (segment%desc%comm%ncpu>1) CALL c_f_pointer(segment%desc%shmem_ptr,segment%fptr_r,[dim1])

      ptr(1:dim1) => segment%fptr_r(1:dim1)

      segment%desc%dim1=dim1

      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE shmem_alloc_r_1d


!****************** SUBROUTINE M_shmem_alloc_r_2d *********************
!
! Allocates a 2d real(q) memory array of size (dim1,dim2)
! shared by the MPI ranks in comm and returns a pointer to it:
!
!    segment%ftr(1:dim1,1:dim2)
!
! Actual allocation is only done by MPI rank 0
!
! N.B.: All MPI ranks in comm must call this routine
!
! N.B.II: the MPI ranks in comm must share direct access to the same
!         physical memory (i.e., they must reside on the same physical
!         compute node).
!
!**********************************************************************

      SUBROUTINE shmem_alloc_r_2d(comm,dim1_help,dim2_help,segment,ptr,istat)

      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem)          :: segment
      REAL(q), POINTER        :: ptr(:,:)
      INTEGER                 :: dim1_help,dim2_help
      INTEGER, OPTIONAL       :: istat

      ! local variables
      INTEGER(KIND=SIZE_KIND) :: dim1,dim2
      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype

      INTEGER                 :: ierror

      REAL(q)                 :: rdummy

      dim1 = INT(dim1_help,KIND=SIZE_KIND)
      dim2 = INT(dim2_help,KIND=SIZE_KIND)

      CALL M_barrier(comm)

      IF (PRESENT(istat)) istat=0

      segment%desc%comm => comm

      IF (segment%desc%comm%ncpu>1) THEN
         GET_SIZE(rdummy)
         segmentsize=sizeoftype*MAX(dim1*dim2,1)

         CALL M_getshmem(segment%desc%comm, segmentsize, sizeoftype, &
                         segment%desc%shmem_win, segment%desc%shmem_ptr, ierror )
      ELSE
         ALLOCATE(segment%fptr_r(1:dim1*dim2),STAT=ierror)
      ENDIF

      IF (PRESENT(istat)) THEN
         IF (ierror/=0) THEN
            istat=ierror; RETURN
         ENDIF
      ENDIF

      IF (segment%desc%comm%ncpu>1) CALL c_f_pointer(segment%desc%shmem_ptr,segment%fptr_r,[dim1*dim2])

      ptr(1:dim1,1:dim2) => segment%fptr_r(1:dim1*dim2)

      segment%desc%dim1=dim1; segment%desc%dim2=dim2

      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE shmem_alloc_r_2d


!****************** SUBROUTINE M_shmem_alloc_r_3d *********************
!
! Allocates a 3d real(q) memory array of size (dim1,dim2,dim3)
! shared by the MPI ranks in comm and returns a pointer to it:
!
!    segment%ftr(1:dim1,1:dim2,1:dim3)
!
! Actual allocation is only done by MPI rank 0
!
! N.B.: All MPI ranks in comm must call this routine
!
! N.B.II: the MPI ranks in comm must share direct access to the same
!         physical memory (i.e., they must reside on the same physical
!         compute node).
!
!**********************************************************************

      SUBROUTINE shmem_alloc_r_3d(comm,dim1_help,dim2_help,dim3_help,segment,ptr,istat)

      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem)          :: segment
      REAL(q), POINTER        :: ptr(:,:,:)
      INTEGER                 :: dim1_help,dim2_help,dim3_help
      INTEGER, OPTIONAL       :: istat

      ! local variables
      INTEGER(KIND=SIZE_KIND) :: dim1,dim2,dim3
      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype

      INTEGER                 :: ierror

      REAL(q)                 :: rdummy

      dim1 = INT(dim1_help,KIND=SIZE_KIND)
      dim2 = INT(dim2_help,KIND=SIZE_KIND)
      dim3 = INT(dim3_help,KIND=SIZE_KIND)

      CALL M_barrier(comm)

      IF (PRESENT(istat)) istat=0

      segment%desc%comm => comm

      IF (segment%desc%comm%ncpu>1) THEN
         GET_SIZE(rdummy)
         segmentsize=sizeoftype*MAX(dim1*dim2*dim3,1)

         CALL M_getshmem(segment%desc%comm, segmentsize, sizeoftype, &
                         segment%desc%shmem_win, segment%desc%shmem_ptr, ierror )
      ELSE
         ALLOCATE(segment%fptr_r(1:dim1*dim2*dim3),STAT=ierror)
      ENDIF

      IF (PRESENT(istat)) THEN
         IF (ierror/=0) THEN
            istat=ierror; RETURN
         ENDIF
      ENDIF

      IF (segment%desc%comm%ncpu>1) CALL c_f_pointer(segment%desc%shmem_ptr,segment%fptr_r,[dim1*dim2*dim3])

      ptr(1:dim1,1:dim2,1:dim3) => segment%fptr_r(1:dim1*dim2*dim3)

      segment%desc%dim1=dim1; segment%desc%dim2=dim2; segment%desc%dim3=dim3

      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE shmem_alloc_r_3d


!****************** SUBROUTINE M_getshmem *****************************
!
!**********************************************************************

      SUBROUTINE M_getshmem(comm,segmentsize,sizeoftype,shmem_win,shmem_ptr,istat)
#ifndef sysv

      TYPE (ML_MPI_PAR)         :: comm
      TYPE(MPI_Win)           :: shmem_win
      TYPE (c_ptr)            :: shmem_ptr

      INTEGER(KIND=SIZE_KIND) :: segmentsize
      INTEGER                 :: sizeoftype
      INTEGER                 :: istat
      ! local variables
      TYPE (MPI_comm)         :: mpi_comm_f08
      INTEGER(KIND=SIZE_KIND) :: ssize
      INTEGER                 :: disp_unit
      INTEGER                 :: myid,ierror

      myid=comm%node_me+1
      ssize=segmentsize
      IF (myid/=1) ssize=0

      mpi_comm_f08%mpi_val=comm%mpi_comm
      CALL MPI_Win_allocate_shared(ssize, sizeoftype, MPI_INFO_NULL, mpi_comm_f08, &
                                   shmem_ptr, shmem_win, istat)

      IF (myid/=1) CALL MPI_Win_shared_query(shmem_win, 0, &
                                             ssize, disp_unit, shmem_ptr, ierror)

      RETURN

#else
      TYPE (ML_MPI_PAR)      :: comm
      INTEGER(KIND=c_int)    :: shmem_win
      TYPE (c_ptr)           :: shmem_ptr

      INTEGER(KIND=c_size_t) :: segmentsize
      INTEGER                :: sizeoftype
      INTEGER                :: istat

      ! local variables
      INTEGER(KIND=c_int)     :: shmem_id
      INTEGER                 :: myid
      INTEGER                 :: ierror
      INCLUDE "mpif.h"

      myid=comm%node_me+1

      istat=0

      IF (myid==1) THEN
         CALL getshmem_error(segmentsize,shmem_id)
         IF (shmem_id==-1) istat=istat+1
      ENDIF

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, istat, 1, MPI_INTEGER, &
                         MPI_SUM, comm, ierror)
      IF (istat/=0) RETURN

      CALL M_BCAST(comm,shmem_id,0)
      shmem_win=shmem_id

      CALL attachshmem(shmem_id,shmem_ptr)

      IF (myid==1) CALL destroyshmem(shmem_id)

      CALL M_barrier(comm)

      RETURN

#endif
      END SUBROUTINE M_getshmem


!****************** SUBROUTINE M_shmem_dealloc ************************
!> @brief
!! mpi_shmem::m_shmem_dealloc deallocates a shared memory segment
!!
!! @param segment Shared memory object handle (mpi_shmem::m_shmem)
!!
!! @warning
!!  All MPI ranks in segment\%desc\%comm must call this routine
!!
!**********************************************************************

      SUBROUTINE M_shmem_dealloc(segment)
      TYPE (M_shmem) :: segment
      ! local variables
      INTEGER :: ierror

      IF (.NOT.ASSOCIATED(segment%desc%comm)) RETURN

      IF (segment%desc%comm%ncpu>1) THEN
         CALL M_barrier(segment%desc%comm)
#ifndef sysv
         CALL MPI_Win_fence(MPI_MODE_NOSUCCEED, segment%desc%shmem_win, ierror)
         CALL MPI_Win_free(segment%desc%shmem_win, ierror)
#else
         CALL detachshmem(segment%desc%shmem_ptr)
#endif
         CALL M_barrier(segment%desc%comm)
      ELSE
         IF (ASSOCIATED(segment%fptr_r)) DEALLOCATE(segment%fptr_r)
         IF (ASSOCIATED(segment%fptr_c)) DEALLOCATE(segment%fptr_c)
      ENDIF

      NULLIFY(segment%desc%comm,segment%fptr_r,segment%fptr_c)

      segment%desc%dim1=0; segment%desc%dim2=0; segment%desc%dim3=0

      RETURN
      END SUBROUTINE M_shmem_dealloc


!****************** SUBROUTINE M_getsem *******************************
!> @brief
!! mpi_shmem::m_getsem allocates, creates, and initialises
!! locks\%lock(0:n). These locks are shared by the MPI ranks in comm
!!
!! @param[in]  comm  Communicator (mpimy::ML_MPI_PAR)
!! @param[in]  n     Number of locks requested (integer)
!! @param[out] locks Lock object handle (mpi_shmem::m_shmem_lock)
!!
!! @warning
!!  All MPI ranks in comm must call this routine
!!
!**********************************************************************
!!#define debug
      SUBROUTINE M_getsem(comm,n,locks)
#ifndef sysv

      TYPE (ML_MPI_PAR)     :: comm
      TYPE (M_shmem_lock) :: locks
      INTEGER :: n
      ! local variables
      INTEGER :: i

      IF (comm%ncpu==1) RETURN

      ALLOCATE(locks%lock(0:n))

      DO i=0,n
         CALL MCS_lockCreate(comm,locks%lock(i))
      ENDDO

      locks%ilock0=0

      RETURN

#else
      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (M_shmem_lock)     :: locks
      INTEGER :: n
      INTEGER :: myid

      IF (comm%ncpu==1) RETURN

      locks%comm => comm
      myid=comm%node_me+1
      IF (comm%node_me==1) CALL getsem(n,locks%lock)

      CALL M_bcast_i(comm,locks%lock,1)

      CALL M_barrier(comm)

      RETURN

#endif
      END SUBROUTINE M_getsem


!****************** SUBROUTINE M_locksem ******************************
!> @brief
!! mpi_shmem::m_locksem acquires a lock
!!
!! @param[in,out] locks Lock object handle (mpi_shmem::m_shmem_lock)
!! @param[in]     n     number of the lock (integer)
!!
!! @details
!! For n>0 mpi_shmem::m_locksem will acquire locks\%lock(n)
!!
!! For n=0 mpi_shmem::m_locksem will acquire the "super-lock"
!!        locks\%lock(0):
!!
!! @warning
!!  mpi_shmem::m_locksem will raise an internal error in two cases:
!!     -# when an MPI rank attempts to acquire lock(0) before all
!!        other locks are freed
!!     -# when an MPI rank attempts to acquire any lock(n/=0) while
!!        lock(0) is locked
!! @warning
!!  The n=0 case involves an mpimy::m_barrier call: all MPI ranks in
!!  locks\%lock(0)\%comm have to call mpi_shmem::m_locksem(locks,0)
!!
!**********************************************************************

      SUBROUTINE M_locksem(locks,n)
      USE tutor, ONLY: vtutor
#ifndef sysv

      TYPE (M_shmem_lock) :: locks
      INTEGER :: n

      ! local variables
      INTEGER :: i,isum

      IF (.NOT.ASSOCIATED(locks%lock)) RETURN

      IF (n>0) THEN
         IF (locks%ilock0==1) THEN
            CALL vtutor%bug('M_locksem: internal error: lock0 has not been freed yet', __FILE__, __LINE__)
         ENDIF
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(n)%comm%node_me,'   requesting lock:',n
#endif
         CALL MCS_LockAcquire(locks%lock(n))
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(n)%comm%node_me,' has acquired lock:',n
#endif
         RETURN
      ENDIF

      ! check whether all locks are free
      isum=SUM(locks%lock(:)%iownlock)
      CALL M_sum_i(locks%lock(0)%comm,isum,1)
      IF (isum/=0) THEN
         CALL vtutor%bug('M_locksem: internal error: not all locks were released', __FILE__, __LINE__)
      ENDIF

      ! first rank0 acquires lock0
      IF (locks%lock(0)%comm%node_me==0) THEN
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(i)%comm%node_me,'   requesting lock:',i
#endif
         CALL MCS_LockAcquire(locks%lock(0))
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(i)%comm%node_me,' has acquired lock:',i
#endif
         locks%ilock0=1
      ENDIF

      CALL M_bcast_i(locks%lock(0)%comm,locks%ilock0,1)
      CALL M_barrier(locks%lock(0)%comm)
#ifdef debug
      WRITE(*,'(2(a,i4))') 'rank:',locks%lock(0)%comm%node_me,' ilock0:',locks%ilock0
#endif

      ! now the other ranks try to acquire lock0 as well
      IF (locks%lock(0)%comm%node_me/=0) THEN
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(0)%comm%node_me,'   requesting lock:',n
#endif
         CALL MCS_LockAcquire(locks%lock(0))
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(0)%comm%node_me,' has acquired lock:',n
#endif
      ENDIF

      RETURN

#else

      TYPE (M_shmem_lock) :: locks
      INTEGER :: n

      IF (.NOT.ASSOCIATED(locks%comm)) RETURN

      CALL locksem(locks%lock,n)

      RETURN

#endif
      END SUBROUTINE M_locksem


!****************** SUBROUTINE M_unlocksem ****************************
!> @brief
!! mpi_shmem::m_unlocksem releases a lock
!!
!! @param[in,out] locks Lock object handle (mpi_shmem::m_shmem_lock)
!! @param[in]     n     number of the lock (integer)
!!
!! @details
!! For n>0 mpi_shmem::m_locksem will free locks\%lock(n)
!!
!! For n=0 mpi_shmem::m_locksem will free the "super-lock"
!!        locks\%lock(0):
!!
!! @warning
!!  The only difference between the two cases is that freeing lock(0)
!!  is followed by an mpimy::m_barrier call: all MPI ranks in
!!  locks\%lock(0)\%comm have to call mpi_shmem::m_unlocksem(locks,0)
!!
!**********************************************************************

      SUBROUTINE M_unlocksem(locks,n)
#ifndef sysv

      TYPE (M_shmem_lock) :: locks
      INTEGER :: n

      ! local variables
      INTEGER :: i

      IF (.NOT.ASSOCIATED(locks%lock)) RETURN

      IF (n>0) THEN
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(n)%comm%node_me,' will release lock:',n
#endif
         CALL MCS_LockRelease(locks%lock(n))
#ifdef debug
         WRITE(*,'(2(a,i4))') 'rank:',locks%lock(n)%comm%node_me,' has released lock:',n
#endif
         RETURN
      ENDIF

#ifdef debug
      WRITE(*,'(2(a,i4))') 'rank:',locks%lock(0)%comm%node_me,' will release lock:',n
#endif
      CALL MCS_LockRelease(locks%lock(0))
#ifdef debug
      WRITE(*,'(2(a,i4))') 'rank:',locks%lock(0)%comm%node_me,' has released lock:',n
#endif

      IF (locks%lock(0)%comm%node_me==0) locks%ilock0=0
      CALL M_bcast_i(locks%lock(0)%comm,locks%ilock0,1)
      CALL M_barrier(locks%lock(0)%comm)
#ifdef debug
      WRITE(*,'(2(a,i4))') 'rank:',locks%lock(0)%comm%node_me,' ilock0:',locks%ilock0
#endif
      RETURN

#else

      TYPE (M_shmem_lock) :: locks
      INTEGER :: n

      IF (.NOT.ASSOCIATED(locks%comm)) RETURN

      CALL unlocksem(locks%lock,n)

      RETURN

#endif
      END SUBROUTINE M_unlocksem


!****************** SUBROUTINE M_destroysem ***************************
!> @brief
!! mpi_shmem::m_destroysem destroys and deallocates locks\%lock(*)
!!
!! @param locks Lock object handle (mpi_shmem::m_shmem_lock)
!!
!! @warning
!!  Calling mpi_shmem::m_destroysem before all locks\%lock(*) are
!!  freed raises an internal error
!!
!**********************************************************************

      SUBROUTINE M_destroysem(locks)
#ifndef sysv
      USE tutor, ONLY: vtutor

      TYPE (M_shmem_lock) :: locks

      ! local variables
      INTEGER :: i,isum

      IF (.NOT.ASSOCIATED(locks%lock)) RETURN

      ! check whether all locks are free
      isum=SUM(locks%lock(:)%iownlock)
      CALL M_sum_i(locks%lock(0)%comm,isum,1)
      IF (isum/=0) THEN
         CALL vtutor%bug('M_destroysem: internal error: not all locks were released', __FILE__, __LINE__)
      ENDIF

      DO i=0,SIZE(locks%lock)-1
         CALL MCS_LockDestroy(locks%lock(i))
      ENDDO

      DEALLOCATE(locks%lock)
      NULLIFY(locks%lock)

      RETURN

#else

      TYPE (M_shmem_lock) :: locks

      IF (.NOT.ASSOCIATED(locks%comm)) RETURN

      CALL M_barrier(locks%comm)

      IF (locks%comm%node_me==0) CALL destroysem(locks%lock)

      NULLIFY(locks%comm)

      RETURN

#endif
      END SUBROUTINE M_destroysem

#ifndef sysv
!****************** SUBROUTINE MCS_LockCreate *************************
!
!**********************************************************************

      SUBROUTINE MCS_LockCreate(comm,lock)
      TYPE (ML_MPI_PAR), TARGET :: comm
      TYPE (MCS_lock) :: lock
      ! local variables
      TYPE (MPI_comm)                :: mpi_comm_f08
      TYPE (c_ptr)                   :: baseptr
      INTEGER(KIND=MPI_ADDRESS_KIND) :: winsize
      INTEGER                        :: myid,idummy,sizeoftype,ndim,ierror

      lock%comm => comm

      myid=lock%comm%node_me+1

      CALL MPI_Sizeof(idummy, sizeoftype, ierror)
      ndim=2; IF (myid==1) ndim=ndim+1
      winsize=ndim*sizeoftype

      mpi_comm_f08%mpi_val=lock%comm%mpi_comm
      CALL MPI_Win_allocate(winsize, sizeoftype, MPI_INFO_NULL, mpi_comm_f08, baseptr, lock%lock_win, ierror)

      CALL MPI_Win_lock_all(0, lock%lock_win, ierror)

      CALL c_f_pointer(baseptr,lock%lmem,[ndim])

      lock%lmem(nextRank)=-1; lock%lmem(blocked)=0; lock%iownlock=0
      IF (myid==1) lock%lmem(lockTail)=-1
#ifdef debug
      WRITE(*,'(a,i4,a)')  'rank:',lock%comm%node_me,' has created the lock'
      WRITE(*,'(3(a,i4))') 'rank:',lock%comm%node_me,' blocked:',lock%lmem(blocked),' next:',lock%lmem(nextRank)
#endif
      CALL M_barrier(comm)

      RETURN
      END SUBROUTINE MCS_LockCreate


!****************** SUBROUTINE MCS_LockAcquire ************************
!
!**********************************************************************

      SUBROUTINE MCS_LockAcquire(lock)
      TYPE (MCS_lock) :: lock
      ! local variables
      INTEGER :: predecessor,ierror

      ! when this rank already owns the lock
      IF (lock%iownlock==1) RETURN

      lock%lmem(blocked)=1

      ! Replace (MPI_REPLACE) the number of the rank that should give me the lock, i.e.,
      ! the value of lock%lmem(lockTail) on rank "0", by "lock%comm%node_me" and return
      ! the number of the rank that should pass me the lock in the variable predecessor.
      ! This will signify to any other rank that they will have to obtain this lock from me ...
      CALL MPI_Fetch_and_op(lock%comm%node_me, predecessor, MPI_INTEGER, &
                            0, INT(lockTail-1,KIND=MPI_ADDRESS_KIND), MPI_REPLACE, &
                            lock%lock_win, ierror)

      ! this forces the previous MPI_Fetch_and_op to complete
      CALL MPI_Win_flush(0, lock%lock_win, ierror)

      ! in case some other rank is supposed to pass on this lock
      IF (predecessor/=-1) THEN
         ! Tell rank predecessor that it should pass the lock on to me, i.e.,
         ! replace (MPI_REPLACE) lock%lmem(nextRank) on rank "predecessor" by
         ! with my rank number (lock%comm%node_me).
         CALL MPI_Accumulate(lock%comm%node_me, 1, MPI_INTEGER, &
                             predecessor-1, INT(nextRank-1,KIND=MPI_ADDRESS_KIND), 1, MPI_INTEGER, MPI_REPLACE, &
                             lock%lock_win, ierror)

         ! this forces the previous MPI_Accumulate to complete
         CALL MPI_Win_flush(predecessor-1, lock%lock_win, ierror)
#ifdef debug
         WRITE(*,'(a,i4,a,i4,a)') 'rank:',lock%comm%node_me,' waiting for rank:',predecessor,' to release the lock'
#endif
         ! now wait for rank predecessor to pass on the lock ...
         waiting: DO
            CALL MPI_Win_sync(lock%lock_win, ierror)
            IF (lock%lmem(blocked)==0) EXIT waiting
         ENDDO waiting
      ENDIF
#ifdef debug
      WRITE(*,'(a,i4,a)') 'rank:',lock%comm%node_me,' has acquired the lock'
#endif
      lock%lmem(blocked)=0 ; lock%iownlock=1
#ifdef debug
      WRITE(*,'(3(a,i4))') 'rank:',lock%comm%node_me,' blocked:',lock%lmem(blocked),' next:',lock%lmem(nextRank)
#endif

      RETURN
      END SUBROUTINE MCS_LockAcquire


!****************** SUBROUTINE MCS_LockRelease ************************
!
!**********************************************************************

      SUBROUTINE MCS_LockRelease(lock)
      TYPE (MCS_lock) :: lock
      ! local variables
      INTEGER :: nullrank=-1,izero=0,curtail,ierror

      ! when this rank does not own the lock
      IF (lock%iownlock==0) RETURN

      IF (lock%lmem(nextRank)==-1) THEN
         ! Check whether I (lock%comm%node_me) am the last rank to have requested this
         ! lock, i.e., check whether lock%lmem(lockTail) at rank 0 equals lock%comm%node_me.
         ! If this is the case we set lock%lmem(lockTail)=nullrank (=-1). The original value
         ! of lock%lmem(lockTail) is returned in the varianle curtail.
         CALL MPI_Compare_and_swap(nullrank, lock%comm%node_me, curtail, MPI_INTEGER, &
                                   0, INT(lockTail-1,KIND=MPI_ADDRESS_KIND), &
                                   lock%lock_win, ierror)

         ! this forces the previous MPI_Compare_and_swap to complete
         CALL MPI_Win_flush(0, lock%lock_win, ierror)

         ! when I am the last to have requested this lock
         IF (curtail==lock%comm%node_me) THEN
#ifdef debug
            WRITE(*,'(a,i4,a)') 'rank:',lock%comm%node_me,' has released the lock'
#endif
            lock%iownlock=0 ; RETURN
         ENDIF

         ! in case I am not the last rank to have requested this lock I will
         ! have to wait for this rank to tell me its number ...
         waiting: DO
            CALL MPI_Win_sync(lock%lock_win, ierror)
            IF (lock%lmem(nextRank)/=-1) EXIT waiting
         ENDDO waiting
      ENDIF

      ! pass on the lock to rank lock%lmem(nextRank), i.e., set the variable lock%lmem(blocked)
      ! on rank lock%lmem(nextRank) to izero (=0). This tells the target rank that it now holds the lock.
      ! We use MPI_Accumulate with MPI_REPLACE  because that constitutes an "atomic" put operation.
      CALL MPI_Accumulate(izero, 1, MPI_INTEGER, &
                          lock%lmem(nextRank)-1, INT(blocked-1,KIND=MPI_ADDRESS_KIND), 1, MPI_INTEGER, MPI_REPLACE, &
                          lock%lock_win, ierror)

      ! this forces the previous MPI_Accumulate to complete
      CALL MPI_Win_flush(lock%lmem(nextRank)-1, lock%lock_win, ierror)

#ifdef debug
      WRITE(*,'(2(a,i4))') 'rank:',lock%comm%node_me,' has released to lock to rank:',lock%lmem(nextRank)
#endif
      lock%lmem(nextRank)=-1 ; lock%iownlock=0
#ifdef debug
      WRITE(*,'(3(a,i4))') 'rank:',lock%comm%node_me,' blocked:',lock%lmem(blocked),' next:',lock%lmem(nextRank)
#endif

      RETURN
      END SUBROUTINE MCS_LockRelease


!****************** SUBROUTINE MCS_LockDestroy ************************
!
!**********************************************************************

      SUBROUTINE MCS_LockDestroy(lock)
      TYPE (MCS_lock) :: lock
      ! local variables
      INTEGER :: ierror

      CALL M_barrier(lock%comm)

#ifdef debug
      WRITE(*,'(a,i4,a)')  'rank:',lock%comm%node_me,' will destroy the lock'
      WRITE(*,'(3(a,i4))') 'rank:',lock%comm%node_me,' blocked:',lock%lmem(blocked),' next:',lock%lmem(nextRank)
      if (lock%comm%node_me==1) write(*,'(a,i4)') 'tail:',lock%lmem(lockTail)
#endif
      CALL MPI_Win_unlock_all(lock%lock_win, ierror)

      CALL MPI_Win_free(lock%lock_win, ierror)

      CALL M_barrier(lock%comm)

      RETURN
      END SUBROUTINE MCS_LockDestroy
#endif
      END MODULE ml_ff_mpi_shmem
#else
      MODULE ml_ff_mpi_shmem
      CONTAINS
      SUBROUTINE mpi_shmem_dummy
      WRITE(*,*)'Im a DEC compiler so I need this line'
      END SUBROUTINE
      END MODULE ml_ff_mpi_shmem
#endif
