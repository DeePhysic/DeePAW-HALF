#ifdef VASP_HDF5
#include "symbol.inc"
#define BAIL(e)  if (e/=0) return

!> Low level subroutines for writting data to the hdf5 files
module vhdf5_base
  use vhdf5_struct_def
  use mpimy, only: communic
  use hdf5

  implicit none

  !> Interface to read fortran datatypes to the hdf5 file
  interface vh5_read
     module procedure vh5_read_logical_scalar
     module procedure vh5_read_logical_array_1d
     module procedure vh5_read_logical_array_2d
     module procedure vh5_read_logical_array_3d
     module procedure vh5_read_integer_scalar
     module procedure vh5_read_integer_array_1d
     module procedure vh5_read_integer_array_2d
     module procedure vh5_read_integer_array_3d
     module procedure vh5_read_real_scalar
     module procedure vh5_read_real_array_1d
     module procedure vh5_read_real_array_2d
     module procedure vh5_read_real_array_3d
     module procedure vh5_read_real_array_4d
     module procedure vh5_read_double_scalar
     module procedure vh5_read_double_array_1d
     module procedure vh5_read_double_array_2d
     module procedure vh5_read_double_array_3d
     module procedure vh5_read_double_array_4d
     module procedure vh5_read_double_array_5d
     module procedure vh5_read_double_array_6d
     module procedure vh5_read_double_array_7d
     module procedure vh5_read_complex_array_1d
     module procedure vh5_read_complex_array_2d
     module procedure vh5_read_complex_array_3d
     module procedure vh5_read_double_complex_array_1d
     module procedure vh5_read_double_complex_array_2d
     module procedure vh5_read_double_complex_array_3d
     module procedure vh5_read_double_complex_array_4d
     module procedure vh5_read_double_complex_array_5d
     module procedure vh5_read_double_complex_array_6d
     module procedure vh5_read_string_scalar
     module procedure vh5_read_string_array_1d
     !> Interface for parallel reading: ionode reads, then broadcasts to all ranks
     module procedure vh5_read_par_logical_scalar
     module procedure vh5_read_par_integer_scalar
     module procedure vh5_read_par_integer_array_1d
     module procedure vh5_read_par_integer_array_2d
     module procedure vh5_read_par_integer_array_3d
     module procedure vh5_read_par_real_scalar
     module procedure vh5_read_par_real_array_1d
     module procedure vh5_read_par_real_array_2d
     module procedure vh5_read_par_real_array_3d
     module procedure vh5_read_par_real_array_4d
     module procedure vh5_read_par_double_scalar
     module procedure vh5_read_par_double_array_1d
     module procedure vh5_read_par_double_array_2d
     module procedure vh5_read_par_double_array_3d
     module procedure vh5_read_par_double_array_4d
     module procedure vh5_read_par_double_array_5d
     module procedure vh5_read_par_double_array_6d
     module procedure vh5_read_par_double_array_7d
     module procedure vh5_read_par_double_complex_array_1d
     module procedure vh5_read_par_double_complex_array_2d
     module procedure vh5_read_par_double_complex_array_3d
     module procedure vh5_read_par_double_complex_array_4d
     module procedure vh5_read_par_double_complex_array_5d
     module procedure vh5_read_par_double_complex_array_6d
  end interface vh5_read

  !> Interface to write fortran datatypes to the hdf5 file
  interface vh5_write
     module procedure vh5_write_logical_scalar
     module procedure vh5_write_logical_array_1d
     module procedure vh5_write_logical_array_2d
     module procedure vh5_write_logical_array_3d
     module procedure vh5_write_integer_scalar
     module procedure vh5_write_integer_array_1d
     module procedure vh5_write_integer_array_2d
     module procedure vh5_write_integer_array_3d
     module procedure vh5_write_integer_array_4d
     module procedure vh5_write_real_scalar
     module procedure vh5_write_real_array_1d
     module procedure vh5_write_real_array_2d
     module procedure vh5_write_real_array_3d
     module procedure vh5_write_real_array_4d
     module procedure vh5_write_real_array_5d
     module procedure vh5_write_double_scalar
     module procedure vh5_write_double_array_1d
     module procedure vh5_write_double_array_2d
     module procedure vh5_write_double_array_3d
     module procedure vh5_write_double_array_4d
     module procedure vh5_write_double_array_5d
     module procedure vh5_write_double_array_6d
     module procedure vh5_write_double_array_7d
     module procedure vh5_write_double_complex_array_1d
     module procedure vh5_write_double_complex_array_2d
     module procedure vh5_write_double_complex_array_3d
     module procedure vh5_write_double_complex_array_4d
     module procedure vh5_write_double_complex_array_5d
     module procedure vh5_write_double_complex_array_6d
     module procedure vh5_write_complex_array_1d
     module procedure vh5_write_complex_array_2d
     module procedure vh5_write_complex_array_3d
     module procedure vh5_write_complex_array_4d
     module procedure vh5_write_complex_array_5d
     module procedure vh5_write_string_scalar
     module procedure vh5_write_string_array_1d
  end interface vh5_write

  interface vh5_group_open
     module procedure vh5_group_open
     module procedure vh5_group_par_open
  end interface vh5_group_open

  interface vh5_group_close
     module procedure vh5_group_close
     module procedure vh5_group_par_close
  end interface vh5_group_close

  interface vh5_exists
     module procedure vh5_exists
     module procedure vh5_exists_par
  end interface vh5_exists

  interface vh5_get_dimensions
     module procedure vh5_get_dimensions
     module procedure vh5_get_dimensions_par
  end interface vh5_get_dimensions

  public :: vh5_dump_original_to_h5
  public :: read_potcar_from_potcarh5
  public :: read_potcar_from_vaspin

contains

  !======================================================================================
  ! non-hdf5 utility functions
  !======================================================================================
  elemental function logical_to_integer(logical_value) result(logasint)
    logical, intent(in) :: logical_value
    integer:: logasint
    !--------------------------------------------
    if (logical_value) then
       logasint = 1
    else
       logasint = 0
    end if
  end function logical_to_integer

  elemental function integer_to_logical(integer_value) result(intaslog)
    integer, intent(in) :: integer_value
    logical :: intaslog
    !--------------------------------------------
    if (integer_value == 1) then
       intaslog = .true.
    else
       intaslog = .false.
    end if
  end function integer_to_logical

  integer function check_rank(expected_rank, actual_rank) result(ierr)
    integer, intent(in) :: expected_rank, actual_rank
    if (expected_rank /= actual_rank) then
      ierr = actual_rank
    else
      ierr = 0
    end if
  end function check_rank

  !======================================================================================
  ! Initializing HDF5
  !======================================================================================

  !> @brief Initialize hdf5 library
  !>
  !> Needs to be called at the start of the application, before using any hdf5
  !> functionality
  !>
  !> @returns hdf5 lib error code
  function vh5_start() result(ierr)
    integer:: ierr
    !-------------------------
    call h5open_f(ierr); BAIL(ierr)
    call h5pcreate_f(H5P_LINK_CREATE_F, create_inter_prop, ierr); BAIL(ierr)
    call h5pset_create_inter_group_f(create_inter_prop, 1, ierr)
  end function vh5_start

  !> @brief close hdf5 library
  !>
  !> Needs to be called at the end of the application, to free global data structures
  !> allocated by vh5start
  !> @returns hdf5 lib error code
  function vh5_end() result(ierr)
    integer :: ierr
    !--------------------------
    call h5close_f(ierr)
  end function vh5_end

  !======================================================================================
  ! HDF5 file interface: open, close, create and delete h5 files
  !======================================================================================

  !> @brief open vasp hdf5 file for reading
  !>
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_open_read(filename, fileid) result(ierr)
    integer:: ierr
    logical:: file_exists
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    !------------------------------------------------------
    ! open file
    ierr = 0
    inquire(file=filename, exist=file_exists)
    if (.not.file_exists) ierr=1; BAIL(ierr)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, fileid, ierr)
  end function vh5_file_open_read

  !> @brief open vasp hdf5 file for read and write access
  !>
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_open_readwrite(filename, fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    !------------------------------------------------------
    ! open file
    ierr = 0; if (.not.lwriteh5) return
    call h5fopen_f(filename, H5F_ACC_RDWR_F, fileid, ierr)
  end function vh5_file_open_readwrite

  !> @brief close hdf5 file
  !>
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_close(fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    !------------------------------------------------------
!    call h5fflush_f(fileid, ierr)
    ierr = 0
    if( fileid >= 0 ) call h5fclose_f(fileid, ierr)
    fileid = -1
  end function vh5_file_close

  !> @brief create vasp hdf5 file
  !>
  !> create hdf5 file, fail if file already exists
  !>
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @param use_swmr should the file be opened in single-writer-multiple-reader mode
  !> @returns hdf5 lib error code
  function vh5_file_create(filename, fileid, use_swmr, comm) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    logical, optional, intent(in) :: use_swmr
    integer, optional, intent(in) :: comm
    ! local variables
    integer(HID_T) :: plist_id
    integer ierr1
    !------------------------------------------------------
    ! create file
    ierr = 0
    if (.not.lwriteh5) return
    plist_id = vh5_file_access_property(use_swmr, comm)
    CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, fileid, ierr, access_prp = plist_id)
    call h5pclose_f(plist_id, ierr1); VH5_CHECK(ierr1)
  end function vh5_file_create

  !> @brief create vasp hdf5 file
  !>
  !> create hdf5 file, overwrite if file already exists
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_create_or_overwrite(filename, fileid, use_swmr, comm) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    logical, optional, intent(in) :: use_swmr
    integer, optional, intent(in) :: comm
    ! local variables
    integer(HID_T) :: plist_id
    integer :: ierr1
    !------------------------------------------------------
    ! create file
    ierr = 0
    if (.not.lwriteh5) return
    plist_id = vh5_file_access_property(use_swmr, comm)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileid, ierr, access_prp = plist_id)
    call h5pclose_f(plist_id, ierr1); VH5_CHECK(ierr1)
  end function vh5_file_create_or_overwrite

  !> create a property for the appropriate file access
  integer(HID_T) function vh5_file_access_property(use_swmr, comm) result(plist_id)
    logical, intent(in), optional :: use_swmr
    integer, intent(in), optional :: comm
    integer ierr
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr); VH5_CHECK(ierr)
    if (present(comm)) &
        call vh5_file_access_mpiio(plist_id, comm)
    if (present(use_swmr)) &
        call vh5_file_access_swmr(plist_id, use_swmr)
  end function vh5_file_access_property

  !> @brief set access property for MPI-IO
  subroutine vh5_file_access_mpiio(plist_id, comm)
     use tutor, only: vtutor
#ifdef MPI
     include "mpif.h"
#endif
     integer,intent(in) :: comm
     integer(HID_T), intent(inout) :: plist_id
     ! local variables
     integer :: ierr
     integer :: info
     !
     ! Setup file access property list with parallel I/O access.
     !
#if defined(VASP_HDF5_MPIIO) && defined(MPI)
     info = MPI_INFO_NULL
     call h5pset_fapl_mpio_f(plist_id, comm, info, ierr); VH5_CHECK(ierr)
#else
     call vtutor%error('A request for opening the HDF5 file in parallel was made &
                       &but the code was not compiled with parallel HDF5 support. &
                       &You should ensure that the HDF5 library has parallel HDF5 support:\n&
                       &  $ h5fc -showconfig\n&
                       & and recompile vhdf5_base.F with -DVASP_HDF5_MPIIO option.')
#endif
  end subroutine vh5_file_access_mpiio

  !> @brief set access property for SWMR
  subroutine vh5_file_access_swmr(plist_id, use_swmr)
      integer(HID_T), intent(inout) :: plist_id
      logical, intent(in) :: use_swmr
      integer ierr
      if (.not.use_swmr) return
      call h5pset_file_locking_f(plist_id, use_file_locking=.false., &
          ignore_disabled_locks=.true., hdferr=ierr); VH5_CHECK(ierr)
  end subroutine vh5_file_access_swmr


  !> @brief create access property list for MPI-IO
  function vh5_data_access_mpiio() result(plist_id)
     use tutor, only: vtutor
     integer(HID_T) :: plist_id
     ! local variables
     integer :: ierr
     integer :: info
#if defined(VASP_HDF5_MPIIO) && defined(MPI)
     call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr); VH5_CHECK(ierr)
     call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr); VH5_CHECK(ierr)
#else
     call vtutor%error('A request for opening the HDF5 file in parallel was made &
                       &but the code was not compiled with parallel HDF5 support. &
                       &You should ensure that the HDF5 library has parallel HDF5 support:\n&
                       &  $ h5fc -showconfig\n&
                       & and recompile vhdf5_base.F with -DVASP_HDF5_MPIIO option.')
#endif
  end function vh5_data_access_mpiio

  !> @brief close a property list
  function vh5_plist_close(plist) result(ierr)
    integer :: ierr
    integer(HID_T), intent(in) :: plist
    call h5pclose_f(plist, ierr)
  end function vh5_plist_close

  !> @brief close a vasp hdf5 file
  function vh5_file_close_writing(fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    !------------------------------------------------------
!    call h5fflush_f(fileid, ierr)
    ierr = 0; if (.not.lwriteh5) return
    if( fileid >= 0 ) then
       call h5fclose_f(fileid, ierr)
    endif
    fileid = -1
  end function vh5_file_close_writing

  !> @brief delete file
  !>
  !> deletes file (any file, not just hdf5). works by opening the
  !> file with status='old', then closing it with status='delete'.
  !> I.e.: you need to have the proper access rights to open the file,
  !> and you need to pass a unit used for opening.
  !>
  !> @param filename name of file to delete
  !> @param unit io unit to use for deleting the file
  !> @returns error code. 0 if the file could be deleted, otherwise
  !>     the error status of the last open or close command
  function vh5_file_delete(filename, iounit) result(ierr)
    integer:: ierr
    integer:: iounit
    logical:: file_exists
    character(len=*), intent(in) :: filename
    !----------------------------------------
    integer:: istat
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      open(unit=iounit, iostat=istat, file=filename, status='old')
      if (istat == 0) then
        close(iounit, iostat=ierr, status='delete')
      else
        ierr = istat
      endif
    else
      ierr = 0
    endif
  end function vh5_file_delete

  function vh5_file_flush(fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) :: fileid
    !------------------------------------------------------
    ierr = 0; if (.not.lwriteh5 .or. .not.lsynch5) return
    call h5fflush_f(fileid, H5F_SCOPE_LOCAL_F, ierr)
  end function vh5_file_flush

  subroutine vh5_error(ierr, srcfile, srcline)
    use string
    use tutor, only: vtutor
    integer,intent(in) :: ierr, srcline
    character(len=*), intent(in) :: srcfile
    if (ierr==0) return
    call vtutor%bug('HDF5 call in '//srcfile//':'//str(srcline)//' produced error: ' &
       // str(ierr), srcfile, srcline)
  end subroutine vh5_error

  !======================================================================================
  ! HDF5 group interface: create, open, close,... groups in H5 files
  !======================================================================================

  !> @brief open group, fail if it does not exist
  !>
  !> @close locid hdf5 handle (HID_T) to file or group
  !> @param groupname name of group to create
  !> @param groupid handle of the created group
  !> @returns hdf5 lib error code
  function vh5_group_open(locid, groupname, groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   locid
    character(len=*), intent(in) :: groupname
    integer(HID_T), intent(out) ::  groupid
    !----------------------------------------
    call h5gopen_f(locid, groupname, groupid, ierr)
  end function vh5_group_open

  !> @brief close group, fail if it does not exist
  !>
  !> @param groupid hdf5 handle (HID_T) to group
  !> @returns hdf5 lib error code
  function vh5_group_close(groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   groupid
    !----------------------------------------
    call h5gclose_f(groupid, ierr)
  end function vh5_group_close

  !> @brief create group in file, fail if group exists
  !>
  !> @param locid hdf5 handle (HID_T) to file or group
  !> @param groupname name of group to create
  !> @param groupid handle of the created group
  !> @returns hdf5 lib error code
  function vh5_group_create(locid, groupname, groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   locid
    character(len=*), intent(in) :: groupname
    integer(HID_T), intent(out) ::  groupid
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Create a property list and set it to allow creation of intermediate groups
    call h5gcreate_f(locid, groupname, groupid, ierr, lcpl_id = create_inter_prop)
  end function vh5_group_create

  !> @brief create group in file, fail if group exists
  !>
  !> @param locid hdf5 handle (HID_T) to file or group
  !> @param groupname name of group to create
  !> @param groupid handle of the created group
  !> @returns hdf5 lib error code
  function vh5_group_open_or_create(locid, groupname, groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   locid
    character(len=*), intent(in) :: groupname
    integer(HID_T), intent(out) ::  groupid
    !----------------------------------------
    logical :: link_exists
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    link_exists = vh5_exists(locid, groupname); BAIL(ierr)
    if (link_exists) then
       ! group already exists, open it
       ierr = vh5_group_open(locid, groupname, groupid); BAIL(ierr);
    else
       ! group does not exist, create
       ierr = vh5_group_create(locid, groupname, groupid); BAIL(ierr);
    end if
  end function vh5_group_open_or_create

  !> @brief close group if it was opened for writing, fail if it does not exist
  !>
  !> @param groupid hdf5 handle (HID_T) to group
  !> @returns hdf5 lib error code
  function vh5_group_close_writing(groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   groupid
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    call h5gclose_f(groupid, ierr)
  end function vh5_group_close_writing

  !======================================================================================
  ! HDF5 links interface
  !======================================================================================

  !> @brief check for the existence of a link
  !>
  !> @param locid hdf5 handle (HID_T) to group, dataset or attribute
  !> @param link_name name of the link
  !> @returns true if the link is found, false otherwise, in case of internal error stop the code
  recursive function vh5_exists(locid, link_name) result(exists)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: link_name
    logical :: exists
    integer :: ierr
    integer :: sep_loc
    ! Search for '/' from the right
    sep_loc = scan(link_name, '/', .true.)
    if (sep_loc <= 1) then
       ! Reached root/top-level name
       call h5lexists_f(locid, link_name, exists, ierr); VH5_CHECK(ierr)
    else
       ! Recursively call to get parent name
       exists = vh5_exists(locid, link_name(1 : sep_loc - 1))
       if (exists) then
          call h5lexists_f(locid, link_name, exists, ierr); VH5_CHECK(ierr)
       endif
    endif
  end function vh5_exists

  !> @brief delete link if it exists
  !>
  !> @param locid hdf5 handle (HID_T) to group, dataset or attribute
  !> @param link_name name of the link
  !> @return hdf5 lib error code
  function vh5_delete(locid, link_name) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: link_name
    integer :: ierr
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    if (vh5_exists(locid, link_name)) then
       call h5ldelete_f(locid, link_name, ierr); BAIL(ierr)
    endif
  end function vh5_delete

  !======================================================================================
  ! HDF5 attributes interface
  !======================================================================================
  !> @brief Add a simple attribute to a dataset
  !>
  function vh5_add_attribute_dataset(groupid, dataset_name, attr_name, attr_val, skip_present) result(ierr)
    integer(HID_T), intent(in) :: groupid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(in) :: attr_name
    character(LEN=*), intent(in) :: attr_val
    logical, optional, intent(in) :: skip_present
    !---------------------------------------------------------
    integer(HID_T) :: dset_id, attr_id, aspace, string_type
    integer(HSIZE_T) :: dimensions(1)
    integer(HSIZE_T) :: strlen
    logical :: my_skip_present, exists
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create string datatype
    strlen = len(attr_val)
    if (strlen == 0) ierr = 1
    BAIL(ierr);
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    ! Open dataspace
    call h5dopen_f(groupid, dataset_name, dset_id, ierr); BAIL(ierr)
    call h5screate_f(H5S_SCALAR_F, aspace, ierr); BAIL(ierr)
    call h5aexists_f(dset_id, attr_name, exists, ierr); BAIL(ierr)
    if (exists) then
       if (my_skip_present) return
       call h5aopen_f(dset_id, attr_name, attr_id, ierr); BAIL(ierr)
    else
       call h5acreate_f(dset_id, attr_name, string_type, aspace, attr_id, ierr); BAIL(ierr)
    endif
    ! Write the data
    call h5awrite_f(attr_id, string_type, attr_val, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5tclose_f(string_type, ierr); BAIL(ierr)
    call h5aclose_f(attr_id, ierr); BAIL(ierr)
    call h5sclose_f(aspace, ierr); BAIL(ierr)
    call h5dclose_f(dset_id, ierr)
  end function vh5_add_attribute_dataset

  !======================================================================================
  ! HDF5 dataset interface: read, write datasets
  !======================================================================================

  !> @brief if link to dataset exists open it
  !>
  !> @param dataset_name is the path to the dataset to open
  !> @param if the link is not found ierr=3
  !> @returns hdf5 handle (HID_T) to dataset
  function vh5_open_dataset_ifexists(locid, dataset_name, dset) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer(HID_T), intent(out) :: dset
    !----------------------------------------
    integer :: ierr
    !----------------------------------------
    if (vh5_exists(locid, dataset_name)) then
        call h5dopen_f(locid, dataset_name, dset, ierr)
    else
        ierr = 3
    endif
  end function vh5_open_dataset_ifexists

  !> @brief if link to dataset exists open it otherwise create it
  function vh5_create_open_dataset(locid, dataset_name, dset) result(ierr)
    integer(HID_T), intent(in) :: locid          !> hdf5 handle (HID_T) to root group
    character(LEN=*), intent(in) :: dataset_name !> string with path to the dataset to open
    integer(HID_T), intent(out) :: dset          !> hdf5 handle (HID_T) to dataset
    !------------------------------------------
    integer :: ierr
    integer(HID_T) :: dspace
    !------------------------------------------
    if (vh5_exists(locid, dataset_name)) then
       ! Open dataset
       call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
       ! Create dataspace
       call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
       ! Create the dataset
       call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
       ! Close dataspace
       call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
  end function vh5_create_open_dataset

  !---------------------------------------------------------
  ! single precision
  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_real_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: x
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! read the data
    call h5dread_f(dset, H5T_NATIVE_REAL, x, dimensions, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_real_scalar

  function vh5_read_real_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.)), intent(out) :: array(rank)
    !---------------------------------------------------------
    integer :: ierr
    integer(HID_T) :: dspace, memspace, dset
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_REAL, array, count, ierr, memspace, dspace); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_real_subarray_nd

  function vh5_read_real_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    real(kind(0.)), intent(out) :: array(*)
    !---------------------------------------------------------
    integer :: ierr, rank_or_error
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: maxdims(rank)
    !---------------------------------------------------------
    !
    dimensions(:rank) = 0
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_REAL, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_real_array_nd

  function vh5_read_real_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    !---------------------------------------------------------
    integer :: ierr
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 1, my_dimensions, array)
    dimensions = my_dimensions(1)
  end function vh5_read_real_array_1d

  function vh5_read_real_array_2d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:,:)
    integer, optional, intent(out) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(2)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 2, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_real_array_2d

  function vh5_read_real_array_3d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:,:,:)
    integer, optional, intent(out) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(3)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 3, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_real_array_3d

  function vh5_read_real_array_4d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:,:,:,:)
    integer, optional, intent(out) :: dimensions(4)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(4)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 4, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_real_array_4d

  !---------------------------------------------------------
  ! single precision
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_real_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_REAL, x, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_real_scalar

  function vh5_create_real_array_nd(locid, dataset_name, rank, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Create dataspace
    call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
    ! Create the dataset
    call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_create_real_array_nd

  function vh5_write_real_subarray_nd(locid, dataset_name, rank, start, count, array, extend) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.)), intent(in) :: array(*)
    integer(HSIZE_T), optional, intent(in) :: extend(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i, rank_or_error
    integer(HID_T) :: dspace, memspace, dset
    integer(HSIZE_T) :: dimensions(rank), maxdims(rank), new_dimensions(rank)
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    !
    if (present(extend)) then
        ! Get dimensions
        call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
        ierr = check_rank(rank, rank_or_error); BAIL(ierr)
        new_dimensions = dimensions
        do i=1,rank
            if ((extend(i)>0).and.(start(i)>dimensions(i))) new_dimensions(i) = new_dimensions(i) + max(extend(i),start(i)-dimensions(i))
        end do
        if (any(new_dimensions /= dimensions)) then
            call h5dset_extent_f(dset, new_dimensions, ierr); BAIL(ierr)
            call h5sclose_f(dspace, ierr); BAIL(ierr)
            call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
        endif
    endif
    !
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_REAL, array, count, ierr, memspace, dspace); BAIL(ierr)
    !
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_real_subarray_nd

  function vh5_write_real_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    real(kind(0.)), intent(in) :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Create dataspace
    call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
    ! Create the dataset
    call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_REAL, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_real_array_nd

  function vh5_write_real_array_1d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name).and.my_skip_present) return
    dimensions(1) = size(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 1, dimensions, array)
  end function vh5_write_real_array_1d

  function vh5_write_real_array_2d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(2)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 2, dimensions, array)
  end function vh5_write_real_array_2d

  function vh5_write_real_array_3d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(3)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 3, dimensions, array)
  end function vh5_write_real_array_3d

  function vh5_write_real_array_4d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(4)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 4, dimensions, array)
  end function vh5_write_real_array_4d

  function vh5_write_real_array_5d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(5)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 5, dimensions, array)
  end function vh5_write_real_array_5d

  function vh5_write_real_array_6d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(6)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 6, dimensions, array)
  end function vh5_write_real_array_6d

  !---------------------------------------------------------
  ! double precision
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_double_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if (my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, x, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_write_double_scalar

  function vh5_create_double_array_nd(locid, dataset_name, rank, dimensions, max_dimensions, chunk_size) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer(HSIZE_T), optional, intent(in):: max_dimensions(rank)
    integer(HSIZE_T), optional, intent(in):: chunk_size(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i
    integer(HSIZE_T) :: my_chunk_size(rank)
    integer(HID_T) :: dspace, dset, crp_list
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Create dataspace
    if (present(max_dimensions)) then
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr); BAIL(ierr)
        my_chunk_size = 1
        do i=1,rank
            if (max_dimensions(i) /= H5S_UNLIMITED_F) cycle
            my_chunk_size(i) = vh5_default_chunk_size
        end do
        if (present(chunk_size)) my_chunk_size = chunk_size
        call h5pset_chunk_f(crp_list, rank, my_chunk_size, ierr); BAIL(ierr)
        !create extensible dataset
        call h5screate_simple_f(rank, dimensions, dspace, ierr, max_dimensions); BAIL(ierr)
    else
        call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
    endif
    ! Create the dataset
    if (present(max_dimensions)) then
        call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, crp_list, lcpl_id = create_inter_prop); BAIL(ierr)
        call h5pclose_f(crp_list, ierr); BAIL(ierr)
    else
        call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
    endif
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_create_double_array_nd

  function vh5_write_double_subarray_nd(locid, dataset_name, rank, start, count, array, extend, xfer_prp) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.d0)), intent(in) :: array(*)
    integer(HID_T), optional, intent(in) :: xfer_prp
    integer(HSIZE_T), optional, intent(in) :: extend(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i, rank_or_error
    integer(HID_T) :: dspace, memspace, dset
    integer(HSIZE_T) :: dimensions(rank), maxdims(rank), new_dimensions(rank)
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    !
    if (present(extend)) then
        ! Get dimensions
        call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
        ierr = check_rank(rank, rank_or_error); BAIL(ierr)
        new_dimensions = dimensions
        do i=1,rank
            if ((extend(i)>0).and.(start(i)>dimensions(i))) new_dimensions(i) = new_dimensions(i) + max(extend(i),start(i)-dimensions(i))
        end do
        if (any(new_dimensions /= dimensions)) then
            call h5dset_extent_f(dset, new_dimensions, ierr); BAIL(ierr)
            call h5sclose_f(dspace, ierr); BAIL(ierr)
            call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
        endif
    endif
    !
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, array, count, ierr, memspace, dspace, xfer_prp=xfer_prp); BAIL(ierr)
    !
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_double_subarray_nd

  function vh5_write_integer_subarray_nd(locid, dataset_name, rank, start, count, array, extend) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    integer, intent(in) :: array(*)
    integer(HSIZE_T), optional, intent(in) :: extend(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i, rank_or_error
    integer(HID_T) :: dspace, memspace, dset
    integer(HSIZE_T) :: dimensions(rank), maxdims(rank), new_dimensions(rank)
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    !
    if (present(extend)) then
        ! Get dimensions
        call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
        ierr = check_rank(rank, rank_or_error); BAIL(ierr)
        new_dimensions = dimensions
        do i=1,rank
            if ((extend(i)>0).and.(start(i)>dimensions(i))) new_dimensions(i) = new_dimensions(i) + max(extend(i),start(i)-dimensions(i))
        end do
        if (any(new_dimensions /= dimensions)) then
            call h5dset_extent_f(dset, new_dimensions, ierr); BAIL(ierr)
            call h5sclose_f(dspace, ierr); BAIL(ierr)
            call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
        endif
    endif
    !
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, array, count, ierr, memspace, dspace); BAIL(ierr)
    !
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_integer_subarray_nd

  function vh5_write_double_array_nd(locid, dataset_name, rank, dimensions, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    real(kind(0.d0)), intent(in) :: array(*)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_write_double_array_nd

  function vh5_write_double_array_1d(locid, dataset_name, array, numb, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    integer, optional :: numb
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions(1) = size(array); if (present(numb)) dimensions(1) = numb
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ierr = vh5_write_double_array_nd(locid, dataset_name, 1, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_1d

  function vh5_write_double_array_2d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(2)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 2, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_2d

  function vh5_write_double_array_3d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(3)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 3, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_3d

  function vh5_write_double_array_4d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(4)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 4, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_4d

  function vh5_write_double_array_5d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(5)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 5, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_5d

  function vh5_write_double_array_6d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(6)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 6, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_6d

  function vh5_write_double_array_7d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(7)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 7, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_7d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_double_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: x
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! read the data to the dataset
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, x, dimensions, ierr); BAIL(ierr)
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_read_double_scalar

  function vh5_read_double_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.d0)), intent(out) :: array(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, memspace, dset
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, array, count, ierr, memspace, dspace); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_read_double_subarray_nd

  function vh5_read_double_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    real(kind(0.d0)), intent(out) :: array(*)
    integer :: ierr, rank_or_error
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: maxdims(rank)
    !---------------------------------------------------------
    !
    dimensions(:rank)=0
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_double_array_nd

  function vh5_read_double_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 1, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions(1)
  end function vh5_read_double_array_1d

  function vh5_read_double_array_2d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:)
    integer, optional, intent(out) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(2)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 2, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_2d

  function vh5_read_double_array_3d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:)
    integer, optional, intent(out) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(3)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 3, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_3d

  function vh5_read_double_array_4d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:)
    integer, optional, intent(out) :: dimensions(4)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(4)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 4, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_4d

  function vh5_read_double_array_5d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:,:)
    integer, optional, intent(out) :: dimensions(5)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(5)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 5, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_5d

  function vh5_read_double_array_6d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:,:,:)
    integer, optional, intent(out) :: dimensions(6)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(6)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 6, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_6d

  function vh5_read_double_array_7d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:,:,:,:)
    integer, optional, intent(out) :: dimensions(7)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(7)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 7, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_7d

  !---------------------------------------------------------
  ! complex (double precision)
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_create_double_complex_array_nd(locid, dataset_name, rank, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    ierr = vh5_create_double_array_nd(locid, dataset_name, rank+1, rdimensions); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_create_double_complex_array_nd

  function vh5_write_double_complex_subarray_nd(locid, dataset_name, rank, start, count, array, xfer_prp) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.d0)), intent(in), target :: array(*)
    integer(HID_T), optional, intent(in) :: xfer_prp
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_write_double_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray, xfer_prp = xfer_prp)
  end function vh5_write_double_complex_subarray_nd

 !> @brief  write a complex array of ndimensions
 !>
 !> The strategy here is to cast the nD array to a 1D real array
 !> while still keeping the information about the dimensions
 !>
 function vh5_write_double_complex_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    complex(kind(0.d0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rdimensions)])
    ierr = vh5_write_double_array_nd(locid, dataset_name, rank+1, rdimensions, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_nd

  function vh5_write_double_complex_array_1d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:)
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_2d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_1d

  function vh5_write_double_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:)
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_3d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_2d

  function vh5_write_double_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:)
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_4d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_3d

  function vh5_write_double_complex_array_4d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(5) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_5d(locid, dataset_name, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_4d

  function vh5_write_double_complex_array_5d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(6) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_6d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_5d

  function vh5_write_double_complex_array_6d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(7) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_7d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_6d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_double_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.d0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_read_double_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_read_double_complex_subarray_nd

  function vh5_read_double_complex_array_1d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_2d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_1d

  function vh5_read_double_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_3d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_2d

  function vh5_read_double_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_4d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_3d

  function vh5_read_double_complex_array_4d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(5) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_5d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_4d

  function vh5_read_double_complex_array_5d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(6) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_6d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_5d

  function vh5_read_double_complex_array_6d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(7) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_7d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_6d

  !---------------------------------------------------------
  ! complex
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
   function vh5_create_complex_array_nd(locid, dataset_name, rank, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    ierr = vh5_create_real_array_nd(locid, dataset_name, rank+1, rdimensions); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_create_complex_array_nd

  function vh5_write_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_write_real_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_write_complex_subarray_nd

 !> @brief  write a complex array of ndimensions
 !>
 !> The strategy here is to cast the nD array to a 1D real array
 !> while still keeping the information about the dimensions
 !>
 function vh5_write_complex_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    complex(kind(0.0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rdimensions)])
    ierr = vh5_write_real_array_nd(locid, dataset_name, rank+1, rdimensions, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_nd

  function vh5_write_complex_array_1d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:)
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_2d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_1d

  function vh5_write_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:)
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_3d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_2d

  function vh5_write_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:,:)
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_4d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_3d

  function vh5_write_complex_array_4d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(5) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_5d(locid, dataset_name, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_4d

  function vh5_write_complex_array_5d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(6) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_6d(locid, dataset_name, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_5d


  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_read_real_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_read_complex_subarray_nd

  function vh5_read_complex_array_1d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(out), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    real(kind(0.)), pointer :: rarray(:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_real_array_2d(locid, dataset_name, rarray)
  end function vh5_read_complex_array_1d

  function vh5_read_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(out), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    real(kind(0.)), pointer :: rarray(:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_real_array_3d(locid, dataset_name, rarray)
  end function vh5_read_complex_array_2d

  function vh5_read_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(out), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    real(kind(0.)), pointer :: rarray(:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_real_array_4d(locid, dataset_name, rarray)
  end function vh5_read_complex_array_3d

  !---------------------------------------------------------
  ! integer
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_integer_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_STD_I32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, x, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_integer_scalar

  function vh5_write_integer_array_nd(locid, dataset_name, rank, dimensions, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    integer, intent(in) :: array(*)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_STD_I32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_integer_array_nd

  function vh5_write_integer_array_1d(locid, dataset_name, array, numb, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    integer, optional :: numb
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Check if explicit dimensions are passed
    dimensions(1) = size(array); if (present(numb)) dimensions(1) = numb
    ! Finally write
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 1, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_1d

  function vh5_write_integer_array_2d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(2)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 2, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_2d

  function vh5_write_integer_array_3d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(3)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 3, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_3d

  function vh5_write_integer_array_4d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(4)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 4, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_4d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_integer_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: x
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    !---------------------------------------------------------
    !
    ! Open the dataset
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! Read the data
    call h5dread_f(dset, H5T_NATIVE_INTEGER, x, dimensions, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_integer_scalar

  function vh5_read_integer_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(rank)
    integer, intent(out) :: array(*)
    integer :: ierr, rank_or_error
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: maxdims(rank) ! Array to store max dimension sizes
    !---------------------------------------------------------
    !
    dimensions = 0
    ! Open the dataset
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data
    call h5dread_f(dset, H5T_NATIVE_INTEGER, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_read_integer_array_nd

  function vh5_read_integer_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_integer_array_nd(locid, dataset_name, 1, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions(1)
  end function vh5_read_integer_array_1d

  function vh5_read_integer_array_2d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: array(:,:)
    integer, optional, intent(out) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(2)
    !---------------------------------------------------------
    ierr = vh5_read_integer_array_nd(locid, dataset_name, 2, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_integer_array_2d

  function vh5_read_integer_array_3d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: array(:,:,:)
    integer, optional, intent(out) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(3)
    !---------------------------------------------------------
    ierr = vh5_read_integer_array_nd(locid, dataset_name, 3, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_integer_array_3d

  !---------------------------------------------------------
  ! logical
  !
  ! Note that we can cast integer to logical as the fortran standard
  ! implies that the size of logical is equal to the size of integer
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_logical_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !-------------------------------------------------------
    logical :: my_skip_present
    !-------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! write stuff
    ierr = vh5_write_integer_scalar(locid, dataset_name, logical_to_integer(x), skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_scalar

  function vh5_write_logical_array_1d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create pointer
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ! Write stuff
    ierr = vh5_write_integer_array_1d(locid, dataset_name, iarray, skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_array_1d

  function vh5_write_logical_array_2d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(in) :: array(:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create pointer
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ! Write stuff
    ierr = vh5_write_integer_array_2d(locid, dataset_name, iarray, skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_array_2d

  function vh5_write_logical_array_3d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(in) :: array(:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:,:)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create pointer
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ! Write stuff
    ierr = vh5_write_integer_array_3d(locid, dataset_name, iarray, skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_array_3d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_logical_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, intent(out) :: x
    integer :: y, ierr
    !-------------------------------------------------------
    ierr = vh5_read_integer_scalar(locid, dataset_name, y)
    x = integer_to_logical(y)
  end function vh5_read_logical_scalar

  function vh5_read_logical_array_1d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(out) :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:)
    !---------------------------------------------------------
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ierr = vh5_read_integer_array_1d(locid, dataset_name, iarray)
  end function vh5_read_logical_array_1d

  function vh5_read_logical_array_2d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(out) :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:)
    !---------------------------------------------------------
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ierr = vh5_read_integer_array_2d(locid, dataset_name, iarray)
  end function vh5_read_logical_array_2d

  function vh5_read_logical_array_3d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(out) :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:,:)
    !---------------------------------------------------------
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ierr = vh5_read_integer_array_3d(locid, dataset_name, iarray)
  end function vh5_read_logical_array_3d

  !---------------------------------------------------------
  ! string (aka character array)
  !
  ! Note: I didn't implement writing string arrays, because without
  ! jagged arrays, I don't think they are really used in Fortran
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_string_scalar(locid, dataset_name, str, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(in) :: str
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset, string_type
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    integer(HSIZE_T) :: strlen
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    strlen = len(str)
    if (strlen == 0) ierr = 1
    BAIL(ierr);
    ! Create datatype
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if (my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, string_type, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, string_type, str, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5tclose_f(string_type, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_string_scalar

  function vh5_write_string_array_nd(locid, dataset_name, rank, dimensions, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    character(LEN=*), intent(in) :: array(*)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset, string_type
    integer(HSIZE_T) :: strlen
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    strlen = len(array(1))
    if (strlen == 0) ierr = 1
    BAIL(ierr);
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if (my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace dataset
      call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, string_type, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, string_type, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5tclose_f(string_type, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_string_array_nd

  function vh5_write_string_array_1d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(in) :: array(:)
    integer :: ierr
    logical, optional, intent(in) :: skip_present
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name).and.my_skip_present) return
    ! Get dimensions and write data
    dimensions(1) = size(array)
    ierr = vh5_write_string_array_nd(locid, dataset_name, 1, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_string_array_1d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_string_scalar(locid, dataset_name, str) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=:), allocatable, intent(inout) :: str
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset, file_type, string_type
    integer :: data_class, string_class
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    integer(HSIZE_T) :: strlen
    !---------------------------------------------------------
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! Read the data to the dataset
    call H5dget_type_f(dset, file_type, ierr); BAIL(ierr)
    call h5tget_class_f(file_type, data_class, ierr); BAIL(ierr)
    call h5tget_class_f(H5T_FORTRAN_S1, string_class, ierr); BAIL(ierr)
    if (data_class /= string_class) then
      ierr = 5
    else
      call H5tget_size_f(file_type, strlen, ierr); BAIL(ierr)
      if (allocated(str)) deallocate(str)
      allocate(character(len=strlen) :: str)
      call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
      call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
      !
      call h5dread_f(dset, string_type, str, dimensions, ierr); BAIL(ierr)
    endif
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5tclose_f(string_type, ierr); BAIL(ierr)
  end function vh5_read_string_scalar

  function vh5_read_string_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    character(LEN=*), intent(out) :: array(*)
    integer :: ierr, rank_or_error
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset, string_type
    integer(HSIZE_T) :: strlen
    integer(HSIZE_T) :: maxdims(rank)
    !---------------------------------------------------------
    dimensions(:rank) = 0
    strlen = len(array(1))
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data
    call h5dread_f(dset, string_type, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5tclose_f(string_type, ierr); BAIL(ierr)
  end function vh5_read_string_array_nd

  function vh5_read_string_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_string_array_nd(locid, dataset_name, 1, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions(1)
  end function vh5_read_string_array_1d

  !================================================
  !  Write file contents in HDF5
  !================================================

  !> @brief dump a file in text format to hdf5
  function vh5_dump_original_to_h5(groupid, subgroup, filename) result(ierr)
    use string, only: string_from_file
    integer(HID_T), intent(in)       :: groupid
    character(LEN=*), intent(in)    :: subgroup
    character(LEN=*), intent(in)    :: filename
    character(len=:), allocatable :: content
    integer(HID_T) :: subgroupid
    integer :: ierr
    content = string_from_file(filename, ierr); BAIL(ierr)
    VH5_CHECK(vh5_group_open_or_create(groupid, subgroup, subgroupid))
    VH5_CHECK(vh5_write(subgroupid, 'content', content))
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    deallocate(content)
  end function vh5_dump_original_to_h5

  ! @brief function for retrieving potcar from h5
  function read_potcar_from_potcarh5(filenameh5, potentialtype, potcarname, sha256, content) result(ierr)

    character(len=*), intent(in) :: potentialtype
    character(len=*), intent(in) :: potcarname
    character(len=*), intent(inout) :: sha256
    character(len=:), allocatable :: content

    character(len=*) :: filenameh5
    ! needs to be larger by 1 --> 65 otherwise last character is truncated
    character(len=65) :: name_buffer ! buffer to hold object's name
    integer(HID_T) :: fileid, groupid, potid, contentid, datatype, string_type
    integer(SIZE_T) :: datasize, buf_size
    integer(HSIZE_T), dimension(1) :: dimensions ! dummy, not used
    integer :: n, i
    integer :: ih5err
    integer :: ierr
    integer :: nmembers, npotcars_found
    integer :: type
    logical :: exists
    integer :: length_user_sha256
    logical :: found_sha256

    BAIL(vh5_file_open_read(filenameh5, fileid))
    BAIL(vh5_group_open(fileid, potentialtype, groupid))
!    write(*,*) 'Search for '//potcarname
    call h5lexists_f(groupid, potcarname, exists, ih5err); VH5_CHECK(ih5err)
    if (.not.exists) then
      VH5_CHECK(vh5_group_close(groupid))
      VH5_CHECK(vh5_file_close(fileid))
      ierr = -1
      return
    end if
    call h5gopen_f(groupid, potcarname, potid, ih5err); VH5_CHECK(ih5err)
    ! check for SHA256
    !write(*,*)'CREATE_POTCAR, sha256: ', sha256
    if (sha256=='') then
      ! get latest attribute
      call h5lexists_f(potid, 'latest', exists, ih5err); VH5_CHECK(ih5err)
      if (.not.exists) then
        VH5_CHECK(vh5_group_close(potid))
        VH5_CHECK(vh5_group_close(groupid))
        VH5_CHECK(vh5_file_close(fileid))
        ierr = -1
        return
      end if
      call h5dopen_f(potid, 'latest', contentid, ih5err); VH5_CHECK(ih5err)
      buf_size = len(sha256)
      call h5gget_linkval_f(potid, 'latest', buf_size, sha256, ih5err); VH5_CHECK(ih5err)
    else
      found_sha256=.false.
      call h5gn_members_f(groupid, potcarname, nmembers, ih5err); VH5_CHECK(ih5err)
      npotcars_found = 0
      do i=0,nmembers-1
        call h5gget_obj_info_idx_f(potid, '.' , i, name_buffer, type, ih5err); VH5_CHECK(ih5err)
        length_user_sha256 = index(sha256, ' ') - 1
        if (sha256(1:length_user_sha256)==name_buffer(1:length_user_sha256)) then
          call h5dopen_f(potid, name_buffer, contentid, ih5err); VH5_CHECK(ih5err)
          npotcars_found = npotcars_found + 1
          found_sha256 = .true.
          sha256 = name_buffer(1:len(sha256)-1)
        end if
      enddo
      if (npotcars_found > 1) then
        VH5_CHECK(vh5_group_close(groupid))
        VH5_CHECK(vh5_file_close(fileid))
        ierr = -2
        return
      end if
      if (.not.found_sha256) then
        VH5_CHECK(vh5_group_close(groupid))
        VH5_CHECK(vh5_file_close(fileid))
        ierr = -1
        return
      end if
    end if

    call h5dget_type_f(contentid, datatype, ih5err); VH5_CHECK(ih5err)
    call h5tget_size_f(datatype, datasize, ih5err); VH5_CHECK(ih5err)
    if (allocated(content)) deallocate(content)
    allocate(character(datasize) :: content)
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ih5err); VH5_CHECK(ih5err)
    call h5tset_size_f(string_type, datasize, ih5err); VH5_CHECK(ih5err)
    call h5dread_f(contentid, string_type, content, dimensions, ih5err); VH5_CHECK(ih5err)
    call h5dclose_f(contentid, ih5err); VH5_CHECK(ih5err)
    ! close groups and file
    call h5gclose_f(potid, ih5err); VH5_CHECK(ih5err)
    call h5gclose_f(groupid, ih5err); VH5_CHECK(ih5err)
    ierr = vh5_file_close_writing(fileid)
  end function read_potcar_from_potcarh5

  function read_potcar_from_vaspin(content) result(ierr)
    character(len=:), allocatable :: content
    integer(HID_T) fileid, groupid, potid, contentid, datatype, string_type, dspace
    integer(SIZE_T) datasize
    integer(HSIZE_T), dimension(1) :: dimensions ! dummy, not used
    integer :: n, i
    integer :: ih5err
    integer :: ierr
    integer :: nmembers, npotcars_found
    integer :: type
    logical :: exists
    integer :: length_user_sha256
    logical :: found_sha256

    VH5_CHECK(vh5_group_open(ih5ininputgroup_id, subgrp_potcar, groupid))
    call h5dopen_f(groupid, 'content', contentid, ih5err); VH5_CHECK(ih5err)
    call h5dget_type_f(contentid, datatype, ih5err); VH5_CHECK(ih5err)
    call h5tget_size_f(datatype, datasize, ih5err); VH5_CHECK(ih5err)
    if (allocated(content)) deallocate(content)
    allocate(character(datasize) :: content)
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ih5err); VH5_CHECK(ih5err)
    call h5tset_size_f(string_type, datasize, ih5err); VH5_CHECK(ih5err)
    call h5dread_f(contentid, string_type, content, dimensions, ih5err); VH5_CHECK(ih5err)
    call h5dclose_f(contentid, ih5err); VH5_CHECK(ih5err)
    call h5gclose_f(groupid, ih5err); VH5_CHECK(ih5err)

  end function read_potcar_from_vaspin

  integer function vh5_number_elements(locid, dataset_name, number_elements) result (ierr)
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: dataset_name
    integer, intent(out) :: number_elements
    integer(HSIZE_T) npoints
    integer(HID_T) dspace, dset
    number_elements = 0
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_npoints_f(dspace, npoints, ierr); BAIL(ierr)
    number_elements = npoints
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_number_elements

  !> Get dataset dimensions
  !> returns rank and dimensions of dataset as integer array
  integer function vh5_get_dimensions(locid, dataset_name, dimensions, rank) result (ierr)
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: dataset_name
    integer,intent(out) :: dimensions(:)
    integer,intent(out) :: rank
    ! local variables
    integer(HSIZE_T) :: maxdims(VH5_MAX_RANK)
    integer(HSIZE_T) :: my_dimensions(VH5_MAX_RANK)
    integer(HID_T) :: dspace, dset
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, my_dimensions, maxdims, rank) !TODO
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    dimensions(:rank) = my_dimensions(:rank)
  end function vh5_get_dimensions

!=======================================================================
!
!> Parallel read wrappers: ionode reads from HDF5, then broadcasts
!>
!> These functions provide a clean interface for parallel I/O where
!> only the ionode reads from the HDF5 file and then broadcasts the
!> data to all other ranks in the communicator.
!
!=======================================================================

  !> Get dataset dimensions in parallel
  !> only ionode reads dimensions, then broadcasts to all ranks
  integer function vh5_get_dimensions_par(locid, dataset_name, dimensions, rank, COMM) result (ierr)
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: dataset_name
    integer,intent(out) :: dimensions(:)
    integer,intent(out) :: rank
    TYPE (communic), INTENT(IN) :: COMM
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(VH5_MAX_RANK)
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_get_dimensions(locid, dataset_name, dimensions, rank)
    endif
#ifdef MPI
    call M_bcast_i(COMM, ierr, 1); BAIL(ierr)
    call M_bcast_i(COMM, rank, 1)
    call M_bcast_i(COMM, dimensions, rank)
#endif
  end function vh5_get_dimensions_par


  !> check if link exists in parallel, only ionode checks
  !> but communicates the result to all ranks
  function vh5_exists_par(locid, link_name, comm) result(exists)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: link_name
    type(communic), intent(in) :: comm
    logical :: exists
    integer :: ierr
    !----------------------------------------
    if (comm%node_me == comm%ionode) then
       exists = vh5_exists(locid, link_name)
    end if
    CALLMPI( M_bcast_l(comm,exists,1) )
  end function vh5_exists_par

  !> open group in parallel, only ionode opens the group
  !> but communicates the error code to all ranks
  function vh5_group_par_open(locid, groupname, groupid, comm) result(ierr)
      integer:: ierr
      integer(HID_T), intent(in) ::   locid
      character(len=*), intent(in) :: groupname
      integer(HID_T), intent(out) ::  groupid
      type(communic), intent(in) :: comm
      !----------------------------------------
      if (comm%node_me == comm%ionode) then
       ierr = vh5_group_open(locid, groupname, groupid)
      end if
      CALLMPI( M_bcast_i(comm,ierr,1) )
  end function vh5_group_par_open

  !> close group in parallel, only ionode closes the group
  !> but communicates the error code to all ranks
  function vh5_group_par_close(groupid,comm) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   groupid
    type(communic), intent(in) :: comm
    !----------------------------------------
    if (comm%node_me == comm%ionode) then
       ierr = vh5_group_close(groupid)
    end if
    CALLMPI( M_bcast_i(comm,ierr,1) )
  end function vh5_group_par_close

  !> Read scalar real in parallel
  function vh5_read_par_real_scalar(locid, dataset_name, x, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(inout) :: x
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, x)
    endif
    CALLMPI( M_bcast_s(COMM, x, 1) )
    CALLMPI( M_bcast_i(COMM, ierr, 1) )
  end function vh5_read_par_real_scalar

  !> Read 1D real array in parallel
  function vh5_read_par_real_array_1d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(inout) :: array(:)
    integer, optional, intent(inout) :: dimensions
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 1) )
    CALLMPI( M_bcast_s(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_real_array_1d

  !> Read 2D real array in parallel
  function vh5_read_par_real_array_2d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(inout) :: array(:,:)
    integer, optional, intent(inout) :: dimensions(2)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 2) )
    CALLMPI( M_bcast_s(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_real_array_2d

  !> Read 3D real array in parallel
  function vh5_read_par_real_array_3d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(inout) :: array(:,:,:)
    integer, optional, intent(inout) :: dimensions(3)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 3) )
    CALLMPI( M_bcast_s(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_real_array_3d

  !> Read 4D real array in parallel
  function vh5_read_par_real_array_4d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(inout) :: array(:,:,:,:)
    integer, optional, intent(inout) :: dimensions(4)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 4) )
    CALLMPI( M_bcast_s(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_real_array_4d

  !> Read scalar double in parallel
  function vh5_read_par_double_scalar(locid, dataset_name, x, COMM) result(ierr)
    use mpi
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: x
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, x)
    endif
    CALLMPI( M_bcast_d(COMM, x, 1) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_scalar

  !> Read 1D double array in parallel
  function vh5_read_par_double_array_1d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:)
    integer, optional, intent(inout) :: dimensions
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 1) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_1d

  !> Read 2D double array in parallel
  function vh5_read_par_double_array_2d(locid, dataset_name, array, COMM, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer, optional, intent(inout) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 2) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_2d

  !> Read 3D double array in parallel
  function vh5_read_par_double_array_3d(locid, dataset_name, array, comm, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer, optional, intent(inout) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 3) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_3d

  !> Read 4D double array in parallel
  function vh5_read_par_double_array_4d(locid, dataset_name, array, COMM, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:,:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer, optional, intent(inout) :: dimensions(4)
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 4) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_4d

  !> Read 5D double array in parallel
  function vh5_read_par_double_array_5d(locid, dataset_name, array, COMM, dimensions) result(ierr)
    use mpi
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:,:,:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer, optional, intent(inout) :: dimensions(5)
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 5) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_5d

  !> Read 6D double array in parallel
  function vh5_read_par_double_array_6d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:,:,:,:,:,:)
    integer, optional, intent(inout) :: dimensions(6)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 6) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_6d

  !> Read 7D double array in parallel
  function vh5_read_par_double_array_7d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(inout) :: array(:,:,:,:,:,:,:)
    integer, optional, intent(inout) :: dimensions(7)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 7) )
    CALLMPI( M_bcast_d(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_array_7d

  !> Read 1D double complex array in parallel
  function vh5_read_par_double_complex_array_1d(locid, dataset_name, array, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(inout) :: array(:)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array)
    endif
    CALLMPI( M_bcast_z(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_complex_array_1d

  !> Read 2D double complex array in parallel
  function vh5_read_par_double_complex_array_2d(locid, dataset_name, array, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(inout) :: array(:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array)
    endif
    CALLMPI( M_bcast_z(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_complex_array_2d

  !> Read 3D double complex array in parallel
  function vh5_read_par_double_complex_array_3d(locid, dataset_name, array, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(inout) :: array(:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array)
    endif
    CALLMPI( M_bcast_z(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_complex_array_3d

  !> Read 4D double complex array in parallel
  function vh5_read_par_double_complex_array_4d(locid, dataset_name, array, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(inout) :: array(:,:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array)
    endif
    CALLMPI( M_bcast_z(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_complex_array_4d

  !> Read 5D double complex array in parallel
  function vh5_read_par_double_complex_array_5d(locid, dataset_name, array, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(inout) :: array(:,:,:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array)
    endif
    CALLMPI( M_bcast_z(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_complex_array_5d

  !> Read 6D double complex array in parallel
  function vh5_read_par_double_complex_array_6d(locid, dataset_name, array, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(inout) :: array(:,:,:,:,:,:)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array)
    endif
    CALLMPI( M_bcast_z(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_double_complex_array_6d

  !> Read scalar integer in parallel
  function vh5_read_par_integer_scalar(locid, dataset_name, x, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(inout) :: x
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, x)
    endif
    CALLMPI( M_bcast_i(COMM, x, 1) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_integer_scalar

  !> Read 1D integer array in parallel
  function vh5_read_par_integer_array_1d(locid, dataset_name, array, comm, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(inout) :: array(:)
    TYPE (communic), INTENT(IN) :: COMM
    integer, optional, intent(inout) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 1) )
    CALLMPI( M_bcast_i(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_integer_array_1d

  !> Read 2D integer array in parallel
  function vh5_read_par_integer_array_2d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(inout) :: array(:,:)
    integer, optional, intent(inout) :: dimensions(2)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 2) )
    CALLMPI( M_bcast_i(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_integer_array_2d

  !> Read 3D integer array in parallel
  function vh5_read_par_integer_array_3d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(inout) :: array(:,:,:)
    integer, optional, intent(inout) :: dimensions(3)
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 3) )
    CALLMPI( M_bcast_i(COMM, array, size(array)) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_integer_array_3d

  !> Read scalar logical in parallel
  function vh5_read_par_logical_scalar(locid, dataset_name, x, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, intent(inout) :: x
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, x)
    endif
    CALLMPI( M_bcast_l(comm, x, 1) )
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_logical_scalar

  !> Read 1D string array in parallel
  function vh5_read_par_string_array_1d(locid, dataset_name, array, dimensions, COMM) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(inout) :: array(:)
    integer, optional, intent(inout) :: dimensions
    TYPE (communic), INTENT(IN) :: COMM
    integer :: ierr
    integer :: i
    !---------------------------------------------------------
    ierr = 0
    if (COMM%NODE_ME == COMM%IONODE) then
       ierr = vh5_read(locid, dataset_name, array, dimensions)
    endif
    if (present(dimensions)) CALLMPI( M_bcast_i(COMM, dimensions, 1) )
    do i = 1, size(array)
       CALLMPI( M_bcast_c_char(COMM, array(i), len(array(i))) )
    enddo
    CALLMPI( M_bcast_i(comm, ierr, 1) )
  end function vh5_read_par_string_array_1d

end module vhdf5_base
#else
!--------------------------------------------------------------------
! Dummy module
!--------------------------------------------------------------------
module vhdf5_base
   contains
   subroutine vhdf5_base_dummy
   write(*,*)'Im a DEC compiler so I need this line'
   end subroutine
end module vhdf5_base
#endif
