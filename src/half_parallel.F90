module half_parallel
  use half_kinds, only: dp, i64
#ifdef HALF_HAVE_MPI
  use mpi, only: MPI_Init, MPI_Finalize, MPI_Initialized, MPI_Finalized, MPI_Comm_rank, MPI_Comm_size, &
    MPI_Abort, MPI_Allreduce, MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_MIN, MPI_DOUBLE_PRECISION, MPI_INTEGER8
#endif
  implicit none
  private
  public :: parallel_initialize, parallel_finalize, parallel_abort, parallel_rank, parallel_size, &
    parallel_root, parallel_owns, parallel_sum, parallel_min

  integer, save :: saved_rank = 0, saved_size = 1
  logical, save :: started_here = .false.

  interface parallel_sum
    module procedure parallel_sum_real_1d, parallel_sum_real_2d, parallel_sum_int64_1d
  end interface

contains
  subroutine parallel_initialize()
#ifdef HALF_HAVE_MPI
    logical :: initialized
    integer :: ierr
    call MPI_Initialized(initialized,ierr)
    if(.not.initialized)then
      call MPI_Init(ierr)
      started_here=.true.
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD,saved_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,saved_size,ierr)
#endif
  end subroutine

  subroutine parallel_finalize()
#ifdef HALF_HAVE_MPI
    logical :: finalized
    integer :: ierr
    call MPI_Finalized(finalized,ierr)
    if(started_here.and..not.finalized)call MPI_Finalize(ierr)
#endif
  end subroutine

  subroutine parallel_abort(code)
    integer,intent(in)::code
#ifdef HALF_HAVE_MPI
    integer::ierr
    if(saved_size>1)call MPI_Abort(MPI_COMM_WORLD,code,ierr)
#endif
  end subroutine

  pure integer function parallel_rank()
    parallel_rank=saved_rank
  end function

  pure integer function parallel_size()
    parallel_size=saved_size
  end function

  pure logical function parallel_root()
    parallel_root=saved_rank==0
  end function

  pure logical function parallel_owns(index)
    integer,intent(in)::index
    parallel_owns=mod(index-1,saved_size)==saved_rank
  end function

  subroutine parallel_sum_real_1d(values)
    real(dp),intent(inout)::values(:)
#ifdef HALF_HAVE_MPI
    integer::ierr
    if(saved_size>1)call MPI_Allreduce(MPI_IN_PLACE,values,size(values),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
  end subroutine

  subroutine parallel_sum_real_2d(values)
    real(dp),intent(inout)::values(:,:)
#ifdef HALF_HAVE_MPI
    integer::ierr
    if(saved_size>1)call MPI_Allreduce(MPI_IN_PLACE,values,size(values),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
  end subroutine

  subroutine parallel_sum_int64_1d(values)
    integer(i64),intent(inout)::values(:)
#ifdef HALF_HAVE_MPI
    integer::ierr
    if(saved_size>1)call MPI_Allreduce(MPI_IN_PLACE,values,size(values),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
  end subroutine

  subroutine parallel_min(value)
    real(dp),intent(inout)::value
#ifdef HALF_HAVE_MPI
    integer::ierr
    if(saved_size>1)call MPI_Allreduce(MPI_IN_PLACE,value,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
#endif
  end subroutine
end module half_parallel
