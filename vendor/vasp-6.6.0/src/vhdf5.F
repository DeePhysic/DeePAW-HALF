#include "symbol.inc"
#ifdef VASP_HDF5

!> High level subroutines for writing data to the hdf5 files
!> The routines here should only return an error if that error should be handled by the
!> calling code. There are very few such cases, most errors should be handled inside the
!> routines by calling vtutor%error.
module vhdf5
  use prec
  use poscar, only : type_info, dynamics, nh_chains
  use base, only : info_struct
  use wave, only: wavedes, wavespin
  use lattice, only : latt, dirkar
  use mkpoints, only : kpoints_struct
  use mkpoints_struct_def, only : skpoints_full
  use locproj_struct, only : LPRJ_function
  use vhdf5_base
  use hdf5, only : HID_T, HSIZE_T, SIZE_T

  implicit none
  private

  integer,parameter :: MAX_LEN_GROUP=50
  public :: MAX_LEN_GROUP

  ! Types re-exported from hdf5
  public :: HID_T
  public :: HSIZE_T
  public :: SIZE_T

  ! Filenames, IDs, and group names
  public :: vaspin
  public :: vaspwpot
  public :: vaspwpotfull
  public :: ih5infileid
  public :: ih5outfileid
  public :: ih5wavefileid
  public :: ih5wpotfileid
  public :: ih5wpotfullfileid
  public :: ih5intermediategroup_id
  public :: grp_charge
  public :: grp_elphon
  public :: grp_history
  public :: grp_intermediate
  public :: grp_locpot
  public :: grp_pair_correlation
  public :: grp_potential
  public :: grp_results
  public :: grp_tau
  public :: grp_wanproj

  ! Global state
  public :: hdf5_found
  public :: incar_found
  public :: kpoints_found
  public :: lwriteh5
  public :: poscar_found
  public :: vaspwave_found
  public :: wavecar_found
  public :: wpot_found
  public :: wpotfull_found

  ! HDF5 lifecycle
  public :: start_h5
  public :: stop_h5
  public :: open_vaspin
  public :: open_vaspout
  public :: close_vaspout
  public :: open_vaspwave
  public :: close_vaspwave

  ! Low-level helpers and error handling
  public :: vh5_data_access_mpiio
  public :: vh5_delete
  public :: vh5_error
  public :: vh5_exists
  public :: vh5_file_close
  public :: vh5_file_close_writing
  public :: vh5_file_create_or_overwrite
  public :: vh5_file_delete
  public :: vh5_file_open_read
  public :: vh5_file_open_readwrite
  public :: vh5_get_dimensions
  public :: vh5_group_close
  public :: vh5_group_close_writing
  public :: vh5_group_open_or_create
  public :: vh5_plist_close
  public :: vh5_read
  public :: vh5_write
  public :: vh5_create_double_array_nd
  public :: vh5_create_double_complex_array_nd
  public :: vh5_read_double_complex_subarray_nd
  public :: vh5_write_double_subarray_nd
  public :: vh5_write_double_complex_subarray_nd

  ! High-level read APIs
  public :: vh5_read_force_constants
  public :: vh5_read_pot
  public :: vh5_read_rho_paw
  public :: vh5_read_wanproj
  public :: vh5_read_waveder
  public :: vh5_read_waveder_header
  public :: vh5_read_wavefunction_header
  public :: vh5_read_wavefunctions
  public :: vh5_read_wpot
  public :: vh5_read_wpot_header

  ! High-level write APIs
  public :: kpoints_along_path_hdf5write
  public :: vh5_dump_original_to_h5
  public :: vh5_init_parcharge_output
  public :: vh5_write_atommom
  public :: vh5_write_average_potential
  public :: vh5_write_band_velocities
  public :: vh5_write_bandgap
  public :: vh5_write_born_charges
  public :: vh5_write_bsebands
  public :: vh5_write_bsefatband
  public :: vh5_write_bseindex
  public :: vh5_write_charge
  public :: vh5_write_conductivity_dynamic
  public :: vh5_write_dielectric_dynamic
  public :: vh5_write_dielectric_energies
  public :: vh5_write_dielectric_scalar_dynamic
  public :: vh5_write_dielectric_static
  public :: vh5_write_dielectric_vector
  public :: vh5_write_dipole_moment
  public :: vh5_write_dos
  public :: vh5_write_dosheader
  public :: vh5_write_dynmat
  public :: vh5_write_efermi
  public :: vh5_write_eigenval
  public :: vh5_write_eigenvalheader
  public :: vh5_write_elastic_modulus
  public :: vh5_write_energies_hist_start
  public :: vh5_write_energies_hist_step
  public :: vh5_write_exciton_start
  public :: vh5_write_exciton_step
  public :: vh5_write_force_constants
  public :: vh5_write_grid_quantity
  public :: vh5_write_internal_strain
  public :: vh5_write_lattice_ion
  public :: vh5_write_lattice_ion_hist_start
  public :: vh5_write_lattice_ion_hist_step
  public :: vh5_write_locproj
  public :: vh5_write_moments_hist_start
  public :: vh5_write_moments_hist_step
  public :: vh5_write_nics_for_posnics
  public :: vh5_write_nmrcur_meta
  public :: vh5_write_optical_transitions
  public :: vh5_write_oszicar_end_elmin_step
  public :: vh5_write_oszicar_hist_step
  public :: vh5_write_oszicar_start
  public :: vh5_write_pair_correlation_start
  public :: vh5_write_pair_correlation_step
  public :: vh5_write_paw_dij_qij
  public :: vh5_write_paw_occupancies
  public :: vh5_write_phonon_dos
  public :: vh5_write_phonons
  public :: vh5_write_phonons_simple
  public :: vh5_write_piezoelectric_tensor
  public :: vh5_write_positions_trailer
  public :: vh5_write_potential
  public :: vh5_write_projectors
  public :: vh5_write_self_energy
  public :: vh5_write_spin_resolved_dipole_moment
  public :: vh5_write_stress_forces_hist_start
  public :: vh5_write_stress_forces_hist_step
  public :: vh5_write_to_parcharge_array
  public :: vh5_write_velocity_hist_start
  public :: vh5_write_velocity_hist_step
  public :: vh5_write_version
  public :: vh5_write_wanproj
  public :: vh5_write_waveder
  public :: vh5_write_wavefunctions
  public :: vh5_write_wpot
  public :: vh5_write_xijkl
  public :: vh5_write_xijkl_header

  interface vh5_write_dielectric_energies
    module procedure vh5_write_dielectric_energies_grid, vh5_write_dielectric_energies_mesh, &
                     vh5_write_dielectric_energies_grid_offset
  end interface vh5_write_dielectric_energies

  interface vh5_write_xijkl
    module procedure vh5_write_xijkl_2d, vh5_write_xijkl_2d_real, vh5_write_xijkl_4d, &
      vh5_write_xijkl_5d, vh5_write_xijkl_6d, vh5_write_xijkl_7d
  end interface vh5_write_xijkl

contains

  !> @brief dumps warning if group not present in file
  function vh5_check_group_presence(fileid, comm, iu0, group, parent_group ) result( ierr )
    USE wave
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type( communic ), intent(in) :: comm
    integer, intent(in) :: iu0
    character(len=*), intent(in) :: group !< group name
    character(len=*), intent(in), optional :: parent_group !< parent group name
    integer :: ierr
    !> local variables
    logical :: lexist

    ierr = 0
    ! check if group is present:
    lexist  = vh5_exists(fileid,group,comm)
    if( .NOT. lexist ) then
       if( present( parent_group ) ) then
          if( iu0>=0 ) write(*,*)"vaspwave.h5 does not contain a dataset with name "//trim(parent_group)//"/"//trim(group)
       else
          if( iu0>=0 ) write(*,*)"vaspwave.h5 does not contain a dataset with name "//trim(group)
       endif
       ierr = -2
    endif
  end function vh5_check_group_presence

  !> Initialize the HDF5 library and open vaspin.h5
  subroutine start_h5(input_file)
    character(len=*), intent(in), optional :: input_file  !< overwrite the default filename vaspin.h5
    ! vh5_start() initializes the hdf5 library and needs to be called before using any
    ! hdf5 related functionality
    VH5_CHECK(vh5_start())
    if (hdf5_found) call open_vaspin(input_file)
  end subroutine start_h5

  !> Open vaspin.h5 (or optionally the file specified by input_file)
  subroutine open_vaspin(input_file)
    character(len=*), intent(in), optional :: input_file  !< overwrite the default filename vaspin.h5
    character(len=:), allocatable :: filename
    if (present(input_file)) then
      filename = input_file
    else
      filename = VASPIN
    end if
    VH5_CHECK(vh5_file_open_read(filename, ih5infileid))
    VH5_CHECK(vh5_group_open(ih5infileid, grp_input, ih5ininputgroup_id))
  end subroutine open_vaspin

  !> Create the vaspout.h5 file for the HDF5 output
  subroutine open_vaspout(lsynch5_loc, subdir, output_file)
    logical, intent(in) :: lsynch5_loc  !< synchronize the HDF5 output to allow access with other processes
    character(len=*), intent(in), optional :: subdir  !< subdirectory in which the file is created
    character(len=*), intent(in), optional :: output_file  !< overwrite the default filename vaspout.h5
    character(len=:), allocatable :: filename
    if (present(output_file)) then
      filename = output_file
    else
      filename = VASPOUT
    end if
    if (present(subdir)) &
      filename = subdir // filename
    lsynch5 = lsynch5_loc
    VH5_CHECK(vh5_file_create_or_overwrite(filename, ih5outfileid, lsynch5))
    VH5_CHECK(vh5_group_open_or_create(ih5outfileid, grp_input, ih5outinputgroup_id))
    VH5_CHECK(vh5_group_open_or_create(ih5outfileid, grp_intermediate, ih5intermediategroup_id))
    call vh5_write_version(ih5outfileid)
  end subroutine open_vaspout

  !> Create the vaspwave.h5 file for the HDF5 output
  subroutine open_vaspwave(io, comm, subdir)
    use base
    use mpimy, only: communic
    type(in_struct)      :: io
    type(communic), intent(in) :: comm
    character(len=*), intent(in), optional :: subdir  !< subdirectory in which the file is created
    character(len=:), allocatable :: filename
    logical :: lexist

    if (.not.(io%lh5 .or. io%wrt_density%lchargh5)) return
    ! open vaspwave
    filename = "vaspwave.h5"
    if (present(subdir)) &
      filename = subdir // filename
    inquire( file=filename, exist = lexist )
    if( lexist) then
       VH5_CHECK(vh5_file_open_readwrite(filename,ih5wavefileid))
    else
       VH5_CHECK(vh5_file_delete(filename, 112))
       VH5_CHECK(vh5_file_create_or_overwrite(filename, ih5wavefileid))
       call vh5_write_version(ih5wavefileid)
    endif

  end subroutine open_vaspwave

  !> Close the file opened by open_vaspout
  subroutine close_vaspout
    VH5_CHECK(vh5_group_close_writing(ih5outinputgroup_id))
    VH5_CHECK(vh5_group_close_writing(ih5intermediategroup_id))
    VH5_CHECK(vh5_file_close_writing(ih5outfileid))
  end subroutine close_vaspout

  !> Close vaspwave opened if opened as rw
  subroutine close_vaspwave()
    VH5_CHECK(vh5_file_close_writing(ih5wavefileid))
  end subroutine close_vaspwave

  subroutine stop_h5(subdir)
    use mpimy, only: communic
    character(len=*), intent(in), optional :: subdir
    character(len=:), allocatable :: subdir_
    logical image_incar
    logical image_kpoints
    !
    ! if no INCAR, try to close vaspin.h5
    !
    if (.not. incar_found .and. hdf5_found) then
      VH5_CHECK(vh5_group_close(ih5ininputgroup_id))
      VH5_CHECK(vh5_file_close(ih5infileid))
    endif
    !
    ! write original  input data to vaspout.h5
    !
    if (present(subdir)) then
      subdir_ = subdir
    else
      subdir_ = ''
    end if
    if (poscar_found.or.incar_found.or.kpoints_found) then
      VH5_CHECK(vh5_group_open_or_create(ih5outfileid, grp_original, ih5outoriginalgroup_id))
    endif
    if (poscar_found) then
      VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_poscar, subdir_ // "POSCAR"))
    endif
    if (kpoints_found) then
      inquire(file=subdir_ // 'KPOINTS', exist=image_kpoints)
      if (image_kpoints) then
        VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_kpoints, subdir_ // "KPOINTS"))
      else
        VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_kpoints, "KPOINTS"))
      end if
    endif
    if (incar_found) then
      inquire(file=subdir_ // 'INCAR', exist=image_incar)
      if (image_incar) then
        VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_incar, subdir_ // "INCAR"))
      else
        VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_incar, "INCAR"))
      end if
    endif
    if (poscar_found.or.incar_found.or.kpoints_found) then
      VH5_CHECK(vh5_group_close_writing(ih5outoriginalgroup_id))
    endif
    !
    ! close vaspout.h5
    !
    call close_vaspout()
    call close_vaspwave()
    ! vh5_end() frees memory allocated by vh5_start() and must be called at the end
    ! of any program using hdf5
    VH5_CHECK(vh5_end())
  end subroutine stop_h5

  !> @brief write Vasp version into hdf5 file
  subroutine vh5_write_version(fileid)
    use version, only: major, minor, patch
    integer(HID_T), intent(in) :: fileid
    integer(HID_T) :: groupid
    VH5_CHECK(vh5_group_open_or_create(fileid, grp_version, groupid))
    VH5_CHECK(vh5_write(groupid, "major", major))
    VH5_CHECK(vh5_write(groupid, "minor", minor))
    VH5_CHECK(vh5_write(groupid, "patch", patch))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_version

  !> @brief write header of doscar files.
  subroutine vh5_write_dosheader(fileid, typeinfo, jobpar, mdynamics, wave, lattice, info, kpoints, group, subgroup)
    integer, intent(in)              :: jobpar !< if io%lorbit > 10 then jobpar is 1
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(type_info), intent(in)      :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    type(dynamics), intent(in)       :: mdynamics !< MD data structure (atom positions, velocities, etc.)
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(latt), intent(in)           :: lattice !< lattice description
    type(info_struct), intent(in)    :: info !< info structure
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
    !-------------------------------------------------------------------------------------
    real(q) :: aomega
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    aomega = lattice%omega/typeinfo%nions
    VH5_CHECK(vh5_write(subgroupid, "nions", typeinfo%nions))
    VH5_CHECK(vh5_write(subgroupid, "jobpar", jobpar))
    VH5_CHECK(vh5_write(subgroupid, "aomega", aomega))
    VH5_CHECK(vh5_write(subgroupid, "anorm", lattice%anorm*1E-10_q))
    VH5_CHECK(vh5_write(subgroupid, "potim", mdynamics%potim*1e-15_q))
    VH5_CHECK(vh5_write(subgroupid, "temperature", mdynamics%temp))
    VH5_CHECK(vh5_write(subgroupid, "system", info%sznam1))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dosheader

  !> @brief write header of eigenval files.
  subroutine vh5_write_eigenvalheader(fileid, typeinfo, mdynamics, wave, lattice, info, kpoints, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(type_info), intent(in)      :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    type(dynamics), intent(in)       :: mdynamics !< MD data structure (atom positions, velocities, etc.)
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(latt), intent(in)           :: lattice !< lattice description
    type(info_struct), intent(in)    :: info !< info structure
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    real(q) :: aomega
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    aomega = lattice%omega/typeinfo%nions
    VH5_CHECK(vh5_write(subgroupid, "nions", typeinfo%nions))
    VH5_CHECK(vh5_write(subgroupid, "nblocks", mdynamics%nblock * mdynamics%kblock))
    VH5_CHECK(vh5_write(subgroupid, "ispin", wave%ispin))
    VH5_CHECK(vh5_write(subgroupid, "aomega", aomega))
    VH5_CHECK(vh5_write(subgroupid, "anorm", lattice%anorm*1E-10_q))
    VH5_CHECK(vh5_write(subgroupid, "potim", mdynamics%potim*1e-15_q))
    VH5_CHECK(vh5_write(subgroupid, "temperature", mdynamics%temp))
    VH5_CHECK(vh5_write(subgroupid, "system", info%sznam1))
    VH5_CHECK(vh5_write(subgroupid, "nelectrons", info%nelect))
    VH5_CHECK(vh5_write(subgroupid, "kpoints", kpoints%nkpts))
    VH5_CHECK(vh5_write(subgroupid, "nb_tot", wave%nb_tot))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_eigenvalheader

  !> @brief write electronic energies
  subroutine vh5_write_eigenval(fileid, wave, wavesp, kpoints, kpoints_full, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(wavespin), intent(in)       :: wavesp !< wavefunctions including band index and spin
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    type(skpoints_full), pointer, intent(in), optional :: kpoints_full
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    logical :: write_kpoints_full
    !-------------------------------------------------------------------------------------
    write_kpoints_full = .false.
    if (present(kpoints_full)) then
      ! make sure kpoints_full has been set up
      if (associated(kpoints_full)) write_kpoints_full = .true.
    endif
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", kpoints%vkpt))
    VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_weight", kpoints%wtkpt))
    VH5_CHECK(vh5_write(subgroupid, "eigenvalues", real(wavesp%celtot)))
    VH5_CHECK(vh5_write(subgroupid, "fermiweights", wavesp%fertot))
    ! write tetrahedra information if available
    if (kpoints%ltet) then
       VH5_CHECK(vh5_write(subgroupid, 'num_tetrahedra', kpoints%ntet))
       VH5_CHECK(vh5_write(subgroupid, 'volume_weight_tetrahedra', kpoints%volwgt))
       VH5_CHECK(vh5_write(subgroupid, 'coordinate_id_tetrahedra', kpoints%idtet(:,1:kpoints%ntet)))
     endif
    ! write full kpoint mesh
    if (write_kpoints_full) then
      VH5_CHECK(vh5_write(subgroupid, "kpoint_coords_full", kpoints_full%vkpt))
      VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_weight_full", kpoints_full%wtkpt))
      VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_symop", kpoints_full%isymop))
      VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_mapping", kpoints_full%nequiv))
    endif

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_eigenval

  !> @brief write electronic energies
  subroutine vh5_write_band_velocities(fileid, velocities, group, subgroup)
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    real(q), intent(in) :: velocities(:,:,:,:) !< band velocities
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "velocities", velocities))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_band_velocities

  !> @brief write phonon frequencies and eigenvectors
  subroutine vh5_write_phonons(fileid, prim_lat, prim_pos, prim_nityp, prim_ityp, q_mesh, frequ, evecs, group, subgroup)
    use tutor, ONLY: vtutor
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    real(q), intent(in)               :: prim_lat(3,3) !< lattice vectors of the primitive cell
    real(q), intent(in)               :: prim_pos(:,:) !< primitive positions within the primitive cell
    integer, intent(in)               :: prim_nityp(:) !< number of atoms of each type in the primitive cell
    character(len=*), intent(in)      :: prim_ityp(:)  !< string identifier of each of the atomic types
    type (kpoints_struct), intent(in) :: q_mesh
    real(q), intent(in) :: frequ(:, :)
    complex(q), intent(in) :: evecs(:, :, :)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: nmodes, nqpoints, nions
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "primitive/scale", 1.0_q))
    VH5_CHECK(vh5_write(subgroupid, "primitive/lattice_vectors", prim_lat))
    VH5_CHECK(vh5_write(subgroupid, "primitive/position_ions", prim_pos))
    VH5_CHECK(vh5_write(subgroupid, "primitive/number_ion_types", prim_nityp))
    VH5_CHECK(vh5_write(subgroupid, "primitive/ion_types", prim_ityp))

    VH5_CHECK(vh5_write(subgroupid, "qpoint_coords", q_mesh%vkpt))
    VH5_CHECK(vh5_write(subgroupid, "qpoints_symmetry_weight", q_mesh%wtkpt))

    if (q_mesh%nkpts/=size(frequ,2)) call vtutor%bug('Inconsistent dimension for frequencies array',__FILE__,__LINE__)
    if (q_mesh%nkpts/=size(evecs,3)) call vtutor%bug('Inconsistent dimension for frequencies array',__FILE__,__LINE__)
    nmodes   = size(frequ,1)
    nqpoints = size(frequ,2)
    nions = nmodes/3
    VH5_CHECK(vh5_write(subgroupid, "nions", nions))
    VH5_CHECK(vh5_write(subgroupid, "frequencies", frequ))
    VH5_CHECK(vh5_write(subgroupid, "eigenvectors", reshape(evecs,[3,nions,nmodes,nqpoints])))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_phonons

  !> Write phonon density of states
  subroutine vh5_write_phonon_dos(fileid, energies, phonon_dos, group, subgroup)
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    real(q), intent(in) :: energies(:)
    real(q), intent(in) :: phonon_dos(:, :, :)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: nmodes
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS; if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "dos_mesh", energies))
    VH5_CHECK(vh5_write(subgroupid, "dos", sum(sum(phonon_dos,dim=2),dim=2)))
    VH5_CHECK(vh5_write(subgroupid, "dospar", phonon_dos))
    VH5_CHECK(vh5_write(subgroupid, "directions", ["x","y","z"]))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_phonon_dos

  !> @brief Same as above but simpler since it does not use the q_data structure
  subroutine vh5_write_phonons_simple(fileid, vkpt, eigenvalues, eigenvectors,k_path_length, group, subgroup)
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    real(q),intent(in) :: vkpt(:,:)
    real(q),intent(in) :: eigenvalues(:,:)
    complex(q),intent(in) :: eigenvectors(:,:,:)
    real(q),intent(in) :: k_path_length(:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: i, nmodes, nqpoints
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    integer(SIZE_T) :: start2(2), count2(2), start3(3), count3(3)
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    nqpoints = size(vkpt,2)
    nmodes = size(eigenvectors,1)
    VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", vkpt))
    VH5_CHECK(vh5_write(subgroupid, "k_path_length", k_path_length))

    VH5_CHECK(vh5_write(subgroupid, "nions", nmodes/3))
    VH5_CHECK(vh5_write(subgroupid, "eigenvalues", eigenvalues))
    VH5_CHECK(vh5_write(subgroupid, "eigenvectors", eigenvectors))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_phonons_simple

  !> @brief write electronic density of states
  subroutine vh5_write_dos(fileid, wave, kpoints, dos, dosi, dospar, efermi, nionp, lpar, group, subgroup, scale)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    real(q), intent(in) :: dos(:,:) !< density of states
    real(q), intent(in) :: dosi(:,:) !< integrated dos
    real(q), intent(in) :: dospar(:,:,:,:) !< partial density of states
    real(q), intent(in) :: efermi !< fermi energy
    integer, intent(in) :: nionp !< number of ions
    integer, intent(in) :: lpar
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
    real(q), intent(in),optional :: scale !< scale dos and idos before writing
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: i, nedos
    real(q) :: my_scale, deltae
    real(q) :: energies(size(dos,1))
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    my_scale    = 1.d0; if (present(scale)) my_scale = scale
    nedos = size(dos,1)

    VH5_CHECK(vh5_write(subgroupid, "nionp", nionp))
    VH5_CHECK(vh5_write(subgroupid, "ncdij", wave%ncdij))
    VH5_CHECK(vh5_write(subgroupid, "lpar", lpar))

    VH5_CHECK(vh5_write(subgroupid, "emax", kpoints%emax))
    VH5_CHECK(vh5_write(subgroupid, "emin", kpoints%emin))
    VH5_CHECK(vh5_write(subgroupid, "nedos", nedos))
    VH5_CHECK(vh5_write(subgroupid, "efermi", efermi))

    deltae=(kpoints%emax-kpoints%emin)/(nedos-1)
    do i=1,nedos
        energies(i)=kpoints%emin+deltae*(i-1)
    end do
    VH5_CHECK(vh5_write(subgroupid, "energies", energies))
    VH5_CHECK(vh5_write(subgroupid, "dos", dos*my_scale))
    VH5_CHECK(vh5_write(subgroupid, "dosi", dosi*my_scale))
    VH5_CHECK(vh5_write(subgroupid, "dospar", dospar))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dos

  !> @brief: small helper to only update the fermi level in the h5 archive
  subroutine vh5_write_efermi(fileid, efermi, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    real(q), intent(in) :: efermi !< fermi energy
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "efermi", efermi))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_efermi

  !> @brief: write the directions for the NMR response current
  subroutine vh5_write_nmrcur_meta(fileid,wrt_nmrcur,group, subgroup)
    use string, only: str
    use tutor, only: vtutor
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    integer, intent(in) :: wrt_nmrcur
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_NMRCUR)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    character(len=2) :: d(3)
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS; if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_NMRCUR;  if (present(subgroup)) my_subgroup = trim(subgroup)

    d=["bx","by","bz"]

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    if (wrt_nmrcur==4) then
      VH5_CHECK(vh5_write(subgroupid, "perturbations", d))
    elseif (wrt_nmrcur>0.and.wrt_nmrcur<4) then
      VH5_CHECK(vh5_write(subgroupid, "perturbations", [d(wrt_nmrcur)]))
    else
      call vtutor%bug("unknown setting for wrt_nmrcur. "//str(wrt_nmrcur),__FILE__, __LINE__)
    endif
  end subroutine vh5_write_nmrcur_meta

  !> @brief write lattice parameters and positions
  subroutine vh5_write_lattice_ion(fileid, sznam, typeinfo, scale, amat, selective_dynamics, posions, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    character(LEN=*), intent(in)  :: sznam !< string with system name read from poscar
    type(type_info), intent(in)   :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    real(q), intent(in) :: scale !< scale for the lattice vectors amat
    real(q), intent(in) :: amat(3,3) !< lattice vectors
    logical :: selective_dynamics !< logical determining whether we are running selective dynamics
    real(q), dimension(3,typeinfo%nions), intent(in) :: posions !< ion positions
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_POSITIONS)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: posgroupid
    integer(HID_T) :: unitgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;      if (present(group))    my_group = group
    my_subgroup = GRP_POSITIONS; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), unitgroupid))
    VH5_CHECK(vh5_group_open_or_create(unitgroupid, trim(my_subgroup), posgroupid))

    ! lattice and system information
    VH5_CHECK(vh5_write(posgroupid, "system", sznam))
    VH5_CHECK(vh5_write(posgroupid, "scale", scale))
    VH5_CHECK(vh5_write(posgroupid, "lattice_vectors", amat/scale))
    VH5_CHECK(vh5_write(posgroupid, "direct_coordinates", 1))
    VH5_CHECK(vh5_write(posgroupid, "ion_types", typeinfo%type(1:typeinfo%ntyp)))
    VH5_CHECK(vh5_write(posgroupid, "ion_sha256", typeinfo%sha256(1:typeinfo%ntyp)))
    VH5_CHECK(vh5_write(posgroupid, "number_ion_types", typeinfo%nityp(1:typeinfo%ntyp)))

    !positions of the ions
    VH5_CHECK(vh5_write(posgroupid, "position_ions", posions))

    !selective dynamics
    if (selective_dynamics) then
       VH5_CHECK(vh5_write(posgroupid, "selective_dynamics", selective_dynamics))
       VH5_CHECK(vh5_write(posgroupid, "selective_dynamics_ions", typeinfo%lsfor))
    end if

    VH5_CHECK(vh5_group_close_writing(posgroupid))
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_lattice_ion

  subroutine vh5_write_exciton_start(fileid, subgroup, ngx, ngy, ngz, ncv)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: ngx, ngy, ngz, ncv
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(4), maxdims(4), chunk_size(4)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    dims = [ngx,ngy,ngz,1]
    maxdims = [integer(SIZE_T)::ngx,ngy,ngz,ncv]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "exciton_charge", 4, dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_exciton_start

  subroutine vh5_write_exciton_step(fileid, subgroup, ngx, ngy, ngz, lambda, charge)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: ngx, ngy, ngz, lambda
    complex(q), intent(in) :: charge(:)
    !-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(4), count(4), extend(4)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    start = [1,1,1,lambda]
    count = [ngx,ngy,ngz,1]
    extend = [0,0,0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "exciton_charge", 4, start, count, real(charge), extend))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_exciton_step

  !> @brief start writing of lattice vectors and ion positions history to hdf5 file
  subroutine vh5_write_lattice_ion_hist_start(fileid, subgroup, nions, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    type(dynamics), intent(in) :: dyn
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3), chunk_size(3)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    dims = [3,nions,1]
    maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    chunk_size = [1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "position_ions", 3, dims, maxdims, chunk_size))

    dims = [3,3,1]
    maxdims = [integer(SIZE_T)::3,3,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "lattice_vectors", 3, dims, maxdims))

    VH5_CHECK(vh5_write(groupid, "scale", 1.0_q))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_lattice_ion_hist_start

  !> @brief write the lattice vectors and ion positions for one ionic step to hdf5 file
  subroutine vh5_write_lattice_ion_hist_step(fileid, subgroup, nstep, amat, posions)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: amat(3,3)
    real(q), intent(in) :: posions(:,:)
    !-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    nions = size(posions,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "position_ions", 3, start, count, posions, extend))

    count = [3,3,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "lattice_vectors", 3, start, count, amat, extend))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_lattice_ion_hist_step

  !> create entry in HDF5 file to store the ion velocities
  subroutine vh5_write_velocity_hist_start(fileid, subgroup, nions, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    type(dynamics), intent(in) :: dyn
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3), chunk_size(3)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    dims = [3,nions,1]
    maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    chunk_size = [1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "ion_velocities", 3, dims, maxdims, chunk_size))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_velocity_hist_start

  !> write the ion velocities for one ionic step to hdf5 file
  subroutine vh5_write_velocity_hist_step(fileid, subgroup, nstep, lattice, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    type(latt), intent(in) :: lattice
    type(dynamics), intent(in) :: dyn
    !-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
    real(q), allocatable :: velocities(:,:)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    ! convert velocities to A/fs
    allocate(velocities, source=dyn%vel)
    if (dyn%ibrion /= 44 .and. dyn%ibrion /= 40) then
       do nions = 1, size(velocities, 2)
          velocities(:,nions) = velocities(:,nions) / dyn%potim
          call dirkar(1,velocities(:,nions),lattice%a)
       end do
    end if

    nions = size(velocities,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "ion_velocities", 3, start, count, velocities, extend))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_velocity_hist_step

  !> @brief start writing of enegies history to hdf5 file
  subroutine vh5_write_energies_hist_start(fileid, subgroup, ibrion)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in) :: ibrion
    !-------------------------------------------------------------
    integer :: extent
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(2), maxdims(2)
    character(LEN=24) :: energies_tags(7)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    if (ibrion == 0 .or. ibrion == 4) then
        energies_tags = [ 'ion-electron   TOTEN   ',&
                          'kinetic energy EKIN    ',&
                          'kin. lattice   EKIN_LAT',&
                          'temperature    TEIN    ',&
                          'nose potential ES      ',&
                          'nose kinetic   EPS     ',&
                          'total energy   ETOTAL  ' ]
        VH5_CHECK(vh5_write(groupid, "energies_tags", energies_tags))
        extent = 7
    else
        energies_tags(1:3) = [ 'free energy    TOTEN   ',&
                               'energy without entropy ',&
                               'energy(sigma->0)       ']
        VH5_CHECK(vh5_write(groupid, "energies_tags", energies_tags(1:3)))
        extent = 3
    endif

    ! Create an extensible array that will contain all the energies
    dims = [extent,1]
    maxdims = [integer(SIZE_T)::extent,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "energies", 2, dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_energies_hist_start

  !> @brief write the energies for one ionic step to hdf5 file
  subroutine vh5_write_energies_hist_step(fileid, subgroup, nstep, energies)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: energies(:)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    start = [1,nstep]
    count = [size(energies),1]
    extend = [0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "energies", 2, start, count, energies, extend))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_energies_hist_step

  !> @brief write bse index
  subroutine vh5_write_bseindex(fileid, bseindex, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    integer,intent(in) :: bseindex(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "bse_index", bseindex))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_bseindex

  !> @brief write bands included in bse and the Fermi energy
  subroutine vh5_write_bsebands(fileid, bsebands, efermi, group, subgroup)
    use bse_struct, only : banddesc
    integer(HID_T), intent(in)  :: fileid
    type(banddesc), intent(in) :: bsebands(:)
    real(q), intent(in), optional :: efermi
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    integer :: isp
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "bse_vbmax", (/ (bsebands(isp)%vbmax, isp=1,size(bsebands)) /)))
    VH5_CHECK(vh5_write(subgroupid, "bse_vbmin", (/ (bsebands(isp)%vbmin, isp=1,size(bsebands)) /)))
    VH5_CHECK(vh5_write(subgroupid, "bse_cbmin", (/ (bsebands(isp)%cbmin, isp=1,size(bsebands)) /)))
    VH5_CHECK(vh5_write(subgroupid, "bse_cbmax", (/ (bsebands(isp)%cbmax, isp=1,size(bsebands)) /)))
    if (present(efermi)) then
       VH5_CHECK(vh5_write(subgroupid, "efermi", efermi))
    end if

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_bsebands

  !> @brief write bse fatbands
  subroutine vh5_write_bsefatband(fileid, fatbands, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    complex(q),intent(in) :: fatbands(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "bse_fatbands", fatbands))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_bsefatband

  !> @brief start writing of stress and forces history to hdf5 file
  subroutine vh5_write_stress_forces_hist_start(fileid, subgroup, nions, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    type(dynamics), intent(in) :: dyn
    !-------------------------------------------------------------
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: dims(3), maxdims(3), chunk_size(3)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))

    dims = [3,nions,1]
    maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    chunk_size = [1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    VH5_CHECK(vh5_create_double_array_nd(posgroupid, "forces", 3, dims, maxdims, chunk_size))

    dims = [3,3,1]
    maxdims = [integer(SIZE_T)::3,3,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(posgroupid, "stress", 3, dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(posgroupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_stress_forces_hist_start

  !> @brief write the stress and the forces for one ionic step to hdf5 file
  subroutine vh5_write_stress_forces_hist_step(fileid, subgroup, nstep, stress, forces)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: stress(3,3)
    real(q), intent(in) :: forces(:,:)
    !-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))

    nions = size(forces,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(posgroupid, "forces", 3, start, count, forces, extend))

    count = [3,3,1]
    VH5_CHECK(vh5_write_double_subarray_nd(posgroupid, "stress", 3, start, count, stress, extend))

    VH5_CHECK(vh5_group_close_writing(posgroupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_stress_forces_hist_step

  !> @brief start writing output of OSZICAR to hdf5 file and store constants
  subroutine vh5_write_oszicar_start(fileid, subgroup)
      integer(HID_T), intent(in) :: fileid
      character(len=*), intent(in) :: subgroup
!
      integer(HID_T) :: posgroupid
      integer(SIZE_T) :: dims(2), maxdims(2)
      character(len=6), allocatable :: label(:)

      VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))
      dims = [7, 1]
      maxdims = [integer(SIZE_T)::7, H5S_UNLIMITED_F]
      VH5_CHECK(vh5_create_double_array_nd(posgroupid, "oszicar", 2, dims, maxdims))
      VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))
      ! Used to store if the calculation converged electronically
      dims = [1, 1]
      maxdims = [integer(SIZE_T)::1, H5S_UNLIMITED_F]
      VH5_CHECK(vh5_create_double_array_nd(posgroupid, "electronic_step_converged", 2, dims, maxdims))
      VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))

      allocate(label(7))
      label = [character(len=6) :: "N", "E", "dE", "deps", "ncg", "rms", "rms(c)"]
      VH5_CHECK(vh5_write(posgroupid, "oszicar_label", label))
      VH5_CHECK(vh5_write(posgroupid, "electronic_step_converged_dtype", "bool"))
  end subroutine vh5_write_oszicar_start

  !> @brief write the output data of OSZICAR for one electronic step
  subroutine vh5_write_oszicar_hist_step(fileid, subgroup, convergence_data)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    real(q), intent(in) :: convergence_data(:)
    !
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
    integer, save :: completed_steps = 1

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))
    start = [1,completed_steps]
    count = [7,1]
    extend = [0,1]
    completed_steps = completed_steps + 1
    VH5_CHECK(vh5_write_double_subarray_nd(posgroupid, "oszicar", 2, start, count, convergence_data, extend))
    VH5_CHECK(vh5_group_close_writing(posgroupid))
  end subroutine vh5_write_oszicar_hist_step

  !> @brief write the output of the OSZICAR at the end of the electronic step
  subroutine vh5_write_oszicar_end_elmin_step(fileid, subgroup, is_elmin_converged)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    logical, intent(in)          :: is_elmin_converged

    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
    integer, save :: completed_steps = 1
    integer :: converged(1)

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))
    start = [1, completed_steps]
    count = [1, 1]
    extend = [0, 1]

    completed_steps = completed_steps + 1
    if (is_elmin_converged) then
        converged = 0
    else
        converged = 1
    end if
    VH5_CHECK(vh5_write_integer_subarray_nd(posgroupid, "electronic_step_converged", 2, start, count, converged, extend))
    VH5_CHECK(vh5_group_close_writing(posgroupid))
  end subroutine vh5_write_oszicar_end_elmin_step

  !> @brief start writing of charges and magnetic, or orbital-moments history to hdf5 file
  subroutine vh5_write_moments_hist_start(fileid, subgroup, shape_, prefix, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: shape_(3)  !< [n_orbitals,n_ions,n_components]
    character(LEN=*), intent(in) :: prefix     !< select "spin" or "orbital" moments
    type(dynamics), intent(in)   :: dyn
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(4), maxdims(4), chunk_size(4)
    character(len=6), allocatable :: component_tags(:)
    character(len=1), allocatable :: orbital_tags(:)
    !-------------------------------------------------------------
    if (prefix == "spin") then
       allocate(component_tags(4),orbital_tags(4))
       component_tags = ['charge', 'x     ', 'y     ', 'z     ']
       orbital_tags = ['s', 'p', 'd', 'f']
    else if (prefix == "orbital") then
       allocate(component_tags(3),orbital_tags(3))
       component_tags = ['x     ', 'y     ', 'z     ']
       orbital_tags = ['p', 'd', 'f']
    end if
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))
    VH5_CHECK(vh5_write(groupid, "magnetism/" // prefix // "_moments/components", component_tags(:shape_(3))))
    VH5_CHECK(vh5_write(groupid, "magnetism/" // prefix // "_moments/orbitals", orbital_tags(:shape_(1))))
    ! Create an extensible array that will contain all the charges and magnetic moments
    dims = [shape_, 1]
    maxdims = [int(shape_, kind=SIZE_T), H5S_UNLIMITED_F]
    chunk_size = [1,1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "magnetism/" // prefix // "_moments/values", size(dims), dims, maxdims, chunk_size))

    VH5_CHECK(vh5_group_close_writing(groupid))
    deallocate(component_tags,orbital_tags)
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_moments_hist_start

  !> @brief write the charges and magnetic, or orbital-moments for one ionic step to hdf5 file
  subroutine vh5_write_moments_hist_step(locid, subgroup, nstep, moments, prefix)
    integer(HID_T), intent(in)   :: locid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep          !< ionic step
    real(q), intent(in)          :: moments(:,:,:) !< values of the moments in shape [n_orbitals,n_ions,n_components]
    character(LEN=*), intent(in) :: prefix         !< select "spin" or "orbital" moments
    !-------------------------------------------------------------------------------------
    integer ierr
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(4), count_(4), extend(4)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(locid, trim(subgroup), groupid))
    start = [1,1,1,nstep]
    count_ = [shape(moments),1]
    extend = [0,0,0,1]
    ierr = vh5_write_double_subarray_nd(groupid, "magnetism/" // prefix // "_moments/values", size(start), &
        start, count_, moments, extend); VH5_CHECK(ierr)

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(locid))
  end subroutine vh5_write_moments_hist_step

  subroutine vh5_write_pair_correlation_start(locid, group, typeinfo, dyn, pair_correlation)
    use base, only: paco_struct
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: group
    type(type_info), intent(in) :: typeinfo
    type(dynamics), intent(in) :: dyn
    type(paco_struct), intent(in) :: pair_correlation
    !
    integer point, num_points, num_pairs, max_steps
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3)
    real(q) step
    integer :: itype, jtype, index_
    real(q), allocatable :: distances(:)
    character(len=5), allocatable :: labels(:)
    !
    VH5_CHECK(vh5_group_open_or_create(locid, trim(group), groupid))
    !
    num_points = pair_correlation%npaco + 1
    num_pairs = typeinfo%ntyp * (typeinfo%ntyp + 1) / 2 + 1  ! add 1 for the total pair correlation function
    max_steps = dyn%nsw / (dyn%nblock * dyn%kblock)
    dims = [num_points, num_pairs, 0]
    maxdims = [num_points, num_pairs, max_steps]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "function", size(dims), dims, maxdims))
    !
    allocate(distances(num_points))
    step = pair_correlation%apaco / pair_correlation%npaco
    distances = [(step * (point - 0.5), point=1, num_points)]
    VH5_CHECK(vh5_write(groupid, "distances", distances))
    !
    allocate(labels(num_pairs))
    index_ = 1
    labels(index_) = "total"
    do itype = 1, typeinfo%ntyp
       do jtype = itype, typeinfo%ntyp
          index_ = index_ + 1
          labels(index_) = trim(typeinfo%type(itype)) // "~" // trim(typeinfo%type(jtype))
       end do
    end do
    VH5_CHECK(vh5_write(groupid, "labels", labels))
    !
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(locid))
    !
  end subroutine vh5_write_pair_correlation_start

  subroutine vh5_write_pair_correlation_step(locid, group, step, dyn, pair_correlation)
    use base, only: paco_struct
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: group
    integer, intent(in) :: step
    type(dynamics), intent(in) :: dyn
    type(paco_struct), intent(in) :: pair_correlation
    !
    integer shape_(2), sample
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count_(3), extend(3)
    !
    shape_ = shape(pair_correlation%sipaco)
    sample = step / (dyn%nblock * dyn%kblock)
    start = [1, 1, sample]
    count_ = [shape_, 1]
    extend = [0, 0, 1]
    VH5_CHECK(vh5_group_open_or_create(locid, trim(group), groupid))
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "function", 3, start, count_, pair_correlation%sipaco, extend))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(locid))
    !
  end subroutine vh5_write_pair_correlation_step


  subroutine vh5_write_bandgap(locid, group, method, bandgap)
    use bandgap_struct
    use tutor, only: vtutor
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: group
    integer, intent(in) :: method
    type(bandgap_info), intent(in) :: bandgap(:)
    !
    character(len=23), parameter :: labels(14) = [ &
        "valence band maximum   ", &
        "conduction band minimum", &
        "direct gap bottom      ", &
        "direct gap top         ", &
        "Fermi energy           ", &
        "kx (VBM)               ", &
        "ky (VBM)               ", &
        "kz (VBM)               ", &
        "kx (CBM)               ", &
        "ky (CBM)               ", &
        "kz (CBM)               ", &
        "kx (direct)            ", &
        "ky (direct)            ", &
        "kz (direct)            "]
    character(len=:), allocatable :: dataset
    real(q) data_(size(labels))
    integer(HID_T) groupid
    integer(SIZE_T) dims(3), maxdims(3), start(3), extend(3)
    integer component
    integer, save :: step(2) = 1
    !
    select case (method)
    case (use_fermi_weights)
        dataset = "band/gap_from_weight"
    case (use_local_fermi_energy)
        dataset = "band/gap_from_kpoint"
    case default
        call vtutor%bug("Bandgap method to write to hdf5 not implemented", __FILE__, __LINE__)
    end select
    !
    dims = [size(labels), 1, 1]
    maxdims = [size(labels, kind=SIZE_T), size(bandgap, kind=SIZE_T), H5S_UNLIMITED_F]
    extend = [0, 1, 1]
    VH5_CHECK(vh5_group_open_or_create(locid, trim(group) // "/electron", groupid))
    if (all(step == 1)) &
      VH5_CHECK(vh5_write(groupid, "band/labels", labels))
    if (step(method) == 1) &
      VH5_CHECK(vh5_create_double_array_nd(groupid, dataset, size(dims), dims, maxdims))
    do component = 1, size(bandgap)
        associate(gap => bandgap(component))
            data_ = [ &
                gap%fundamental_valence%eigenvalue, &
                gap%fundamental_conduction%eigenvalue, &
                gap%direct_valence%eigenvalue, &
                gap%direct_conduction%eigenvalue, &
                gap%fermi_energy, &
                gap%fundamental_valence%kvector, &
                gap%fundamental_conduction%kvector, &
                gap%direct_valence%kvector &
            ]
            start = [1, component, step(method)]
            VH5_CHECK(vh5_write_double_subarray_nd(groupid, dataset, size(dims), start, dims, data_, extend))
        end associate
    end do
    VH5_CHECK(vh5_group_close_writing(groupid))
    step(method) = step(method) + 1
    VH5_CHECK(vh5_file_flush(locid))
    !
  end subroutine vh5_write_bandgap


  !> @brief write optical transitions
  subroutine vh5_write_optical_transitions(fileid, opticaltransitions, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: opticaltransitions(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "opticaltransitions", opticaltransitions))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_optical_transitions

  !> @brief write contcar trailer
  !>
  !> writes the informations contained in the contcar trailer in the "traditional"
  !> VASP output files to the hdf5 file. Writes lattice velocities and velocities
  !>
  !> @param fileid hdf5 file handle
  !> @param lattice lattice description (def. in lattice.F)
  !> @param typeinfo typeinfo data structure (atom types, number of ions, etc. defined in
  !>        poscar.f)
  !> @param mdynamics MD data structure (atom positions, velocities, etc. defined in poscar.f)
  subroutine vh5_write_positions_trailer(fileid, lattice, typeinfo, mdynamics,nhchains, group, subgroup)
    use poscar, only: has_lattice_velocities
    integer(HID_T), intent(in) :: fileid
    type(latt), intent(in):: lattice
    type(type_info), intent(in):: typeinfo
    type(dynamics), intent(in):: mdynamics
    type(nh_chains), intent(in)      :: nhchains !< nose hoover chains thermostat data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_POSITIONS)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    integer:: ni
    real(q), dimension(3,typeinfo%nions) :: velocities_out
    real(q), dimension(3) :: tmp
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;      if (present(group))    my_group = group
    my_subgroup = GRP_POSITIONS; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    if (has_lattice_velocities(mdynamics, lattice)) then
       VH5_CHECK(vh5_write(subgroupid, "lattice_velocities", lattice%avel/mdynamics%potim))
    end if

    if (typeinfo%nionp > typeinfo%nions) then
       VH5_CHECK(vh5_write(subgroupid, "number_empty_sphere_types", typeinfo%nityp(typeinfo%ntyp+1:typeinfo%ntypp)))
       VH5_CHECK(vh5_write(subgroupid, "position_empty_spheres", mdynamics%posion(:,typeinfo%nions+1:typeinfo%nionp)))
    endif

    if (mdynamics%ibrion/=44 .and. mdynamics%ibrion/=40) then
       do ni = 1, typeinfo%nions
          tmp = mdynamics%vel(:,ni)/mdynamics%potim
          call dirkar(1,tmp,lattice%a)
          velocities_out(:,ni) = tmp
       end do
    else
       velocities_out = mdynamics%vel
    end if
    VH5_CHECK(vh5_write(subgroupid, "ion_velocities", velocities_out))
    VH5_CHECK(vh5_write(subgroupid, "direct_coordinates_velocities", 0))
    if (mdynamics%init==1) then
       VH5_CHECK(vh5_write(subgroupid, "dyn_init", mdynamics%init))
       VH5_CHECK(vh5_write(subgroupid, "potim", mdynamics%potim))
       IF (nhchains%LINIT) THEN
         VH5_CHECK(vh5_write(subgroupid, "nhc_thermostat_x", nhchains%x(1:nhchains%NCHAINSMAX)))
         VH5_CHECK(vh5_write(subgroupid, "nhc_thermostat_p", nhchains%p(1:nhchains%NCHAINSMAX)))
       ELSE
         VH5_CHECK(vh5_write(subgroupid, "nose_thermostat", mdynamics%snose))
       ENDIF
       VH5_CHECK(vh5_write(subgroupid, "predictor_coordinates", mdynamics%posion(:,1:typeinfo%nions)))
       VH5_CHECK(vh5_write(subgroupid, "predictor_coordinates_2", mdynamics%d2))
       VH5_CHECK(vh5_write(subgroupid, "predictor_coordinates_3", mdynamics%d3))
    end if

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_positions_trailer

  !> @brief write force constants
  subroutine vh5_write_force_constants(fileid, force_constants, group, subgroup)
    integer(HID_T), intent(in)       :: fileid
    real(q) :: force_constants(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;           if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "force_constants", force_constants))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_force_constants

  !> @brief write dynamical matrix
  subroutine vh5_write_dynmat(fileid, dynmat, eigenvalues, eigenvectors, group, subgroup)
    use constant
    integer(HID_T), intent(in) :: fileid
    real(q), intent(in) :: dynmat(:,:), eigenvalues(:), eigenvectors(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid, ssubgroupid
    complex(q), allocatable :: sqrt_eigenvalues(:)
    character(LEN=40) :: my_group, my_subgroup
    integer num_ions
    !-------------------------------------------------------------------------------------
    allocate(sqrt_eigenvalues(size(eigenvalues)))
    ! eigenvalues of dynmat are -w^2
    sqrt_eigenvalues = convert_frequency(sqrt(abs(eigenvalues)), unit_internal, unit_eV)
    where (eigenvalues > 0)
        sqrt_eigenvalues = sqrt_eigenvalues * cmplx(0.0, 1.0, kind=q)
    end where
    num_ions = size(eigenvalues) / 3
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))
    VH5_CHECK(vh5_write(subgroupid, "hessian", dynmat))

    VH5_CHECK(vh5_group_open_or_create(subgroupid, "dynmat", ssubgroupid))
    VH5_CHECK(vh5_write(ssubgroupid, "eigenvalues", sqrt_eigenvalues))
    VH5_CHECK(vh5_write(ssubgroupid, "eigenvectors", reshape(eigenvectors, [3, num_ions, size(eigenvalues)])))
    VH5_CHECK(vh5_group_close_writing(ssubgroupid))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dynmat

  !> @brief write born effective charges
  subroutine vh5_write_elastic_modulus(fileid, prefix, elastic_modulus, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: elastic_modulus(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_elastic_modulus", elastic_modulus))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_elastic_modulus

  !> @brief write born effective charges
  subroutine vh5_write_born_charges(fileid, born_charges, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: born_charges(:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "born_charges", born_charges))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_born_charges

  !> @brief write internal strain tensor
  subroutine vh5_write_internal_strain(fileid, internal_strain, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: internal_strain(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "internal_strain", internal_strain))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_internal_strain

  !> @brief write piezoelectric tensor
  subroutine vh5_write_piezoelectric_tensor(fileid, prefix, piezoelectric_tensor, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: piezoelectric_tensor(:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_piezoelectric_tensor", piezoelectric_tensor))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_piezoelectric_tensor

  !> @brief write dipole moment
  subroutine vh5_write_dipole_moment(fileid, prefix, dipole_moment, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: dipole_moment(:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_dipole_moment", dipole_moment))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dipole_moment

  !> @brief write spin resolved dipole moment
  subroutine vh5_write_spin_resolved_dipole_moment(fileid, dipole_moment, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: dipole_moment(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "spin_resolved_dipole_moment", dipole_moment))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_spin_resolved_dipole_moment

  !> @brief write dielectric tensor
  subroutine vh5_write_dielectric_static(fileid, dielectric_tag, dielectric_tensor, method, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(LEN=*),intent(in) :: dielectric_tag
    real(q),intent(in) :: dielectric_tensor(:,:)
    character(LEN=*), intent(in), optional :: method
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DIELECTRIC)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, dielectric_tag // "_dielectric_tensor", dielectric_tensor))
    if (present(method)) &
       VH5_CHECK(vh5_write(subgroupid, "method_dielectric_tensor", method))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_static

  subroutine vh5_write_dielectric_energies_grid(fileid, energy_step, num_energies, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q), intent(in) :: energy_step
    integer, intent(in) :: num_energies
    character(len=*), intent(in), optional :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    real(q) :: energies(num_energies)
    integer :: ii
    energies = [(energy_step * (ii - 1), ii = 1, num_energies)]
    call vh5_write_dielectric_energies_mesh(fileid, energies, prefix, group, subgroup)
  end subroutine vh5_write_dielectric_energies_grid

  subroutine vh5_write_dielectric_energies_grid_offset(fileid, energy_step, num_energies, energy_offset, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q), intent(in) :: energy_step
    integer, intent(in) :: num_energies
    real(q), intent(in) :: energy_offset
    character(len=*), intent(in), optional :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    real(q) :: energies(num_energies)
    integer :: ii
    energies = [(energy_offset + energy_step * (ii - 1), ii = 1, num_energies)]
    call vh5_write_dielectric_energies_mesh(fileid, energies, prefix, group, subgroup)
  end subroutine vh5_write_dielectric_energies_grid_offset

  subroutine vh5_write_dielectric_energies_mesh(fileid, energies, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q), intent(in) :: energies(:)
    character(len=*), intent(in), optional :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    character(LEN=:), allocatable :: dataset
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    dataset = "energies_dielectric_function"
    if (present(prefix)) dataset = prefix // "_" // dataset
    VH5_CHECK(vh5_write(subgroupid, dataset, energies))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_energies_mesh

  !> @brief write energy dependent scalar dielectric function
  subroutine vh5_write_dielectric_scalar_dynamic(fileid, dielectric_function, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    complex(q), intent(in) :: dielectric_function(:)
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_dielectric_function", dielectric_function))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_scalar_dynamic

  !> @brief write energy dependent dielectric function
  subroutine vh5_write_dielectric_dynamic(fileid, dielectric_function, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    complex(q), intent(in) :: dielectric_function(:,:,:)
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_dielectric_function", dielectric_function))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_dynamic

  !> @brief write energy dependent conductivity
  subroutine vh5_write_conductivity_dynamic(fileid, nedos, edos, conductivity, efermi, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q),intent(in) :: efermi
    real(q),intent(in) :: edos(nedos), conductivity(nedos,3,3)
    integer,intent(in) :: nedos
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_CONDUCTIVITY)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    real(q) :: energies(nedos)
    integer :: i
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;         if (present(group))    my_group = group
    my_subgroup = GRP_CONDUCTIVITY; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, my_subgroup, subgroupid))

    ! write stuff
    VH5_CHECK(vh5_write(subgroupid, "energies", edos))
    VH5_CHECK(vh5_write(subgroupid, "conductivity_dynamic", conductivity))
    VH5_CHECK(vh5_write(subgroupid, "efermi", efermi))

    ! close/deallocate
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_conductivity_dynamic

  !> @brief write contents of waveder to h5 file, into same group as wavefunction
  subroutine vh5_write_waveder(fileid, wdes, nbands_cder, nodes_in_dielectric_function, wplasmon, &
        cder_between_states)
    use string, only: str
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    type(wavedes), intent(in) :: wdes !< wave function descriptor
    integer, intent(in) :: nbands_cder !< number of occupied bands rounded to MPI ranks, for T>0 all bands
    real(q), intent(in) :: nodes_in_dielectric_function !< position of innermst node in dielectric function
    real(q), intent(in) :: wplasmon(:,:) !< plasmon frequency
    GDEFS,intent(in) :: cder_between_states(:,:,:,:,:)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, k_groupid, spin_groupid
    character(LEN=256) :: spingroupname, kgroupname
    integer :: isp, k
    GDEFS,allocatable :: cbuf(:,:,:)

    if (vh5_exists(fileid, GRP_WAVEDER, wdes%comm)) then
       VH5_CHECK(vh5_delete(fileid, GRP_WAVEDER))
    endif
    ! store to same wavefunctions group
    VH5_CHECK(vh5_group_open_or_create(fileid, GRP_WAVEDER, groupid))

    ! write stuff
    VH5_CHECK(vh5_write(groupid, "nb_tot", wdes%nb_tot ) )
    VH5_CHECK(vh5_write(groupid, "nbands_cder", nbands_cder ) )
    VH5_CHECK(vh5_write(groupid, "nkpts", wdes%nkpts ) )
    VH5_CHECK(vh5_write(groupid, "ispin", wdes%ispin ) )
    VH5_CHECK(vh5_write(groupid, "nodes_in_dielectric_function", nodes_in_dielectric_function) )
    VH5_CHECK(vh5_write(groupid, "wplasmon", wplasmon ) )
    allocate( cbuf( wdes%nb_tot, nbands_cder, 3 ) )
    spin: do isp=1,wdes%ispin
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       VH5_CHECK(vh5_group_open_or_create(groupid, trim(spingroupname), spin_groupid))
       kpoints: do k=1,wdes%nkpts
          write(kgroupname, '(2A)') "kpoint_"//trim(str(k))
          VH5_CHECK(vh5_group_open_or_create(spin_groupid, trim(kgroupname), k_groupid))

         ! restore this in the following way, in order to be able to downsample w.r.t. k-points when reading it
         ! cder_between_states(1:nb_tot,1:nbands_cder,1:nkpts,1:ispin,1:3)
          cbuf = cder_between_states( :, :, k, isp, : )
          VH5_CHECK(vh5_write(k_groupid, "cder_between_states", cbuf ) )
          VH5_CHECK(vh5_group_close_writing(k_groupid))
       enddo kpoints
       VH5_CHECK(vh5_group_close_writing(spin_groupid))
    enddo spin
    deallocate( cbuf )
    ! close/deallocate
    VH5_CHECK(vh5_group_close_writing(groupid))
    !VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_waveder

  !> @brief reads header of waveder
  function vh5_read_waveder_header(fileid, wdes, nbands_cder, nodes_in_dielectric_function, wplasmon, iu0, group) result (ierr)
    use string, only: str
    USE base
    USE tutor, ONLY: vtutor, isAlert, BSEoptics, newArgument
    USE chi_glb, ONLY: LBSE
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    type(wavedes), intent(in) :: wdes !< wave function descriptor
    integer, intent(inout) :: nbands_cder !< number of occupied bands rounded to MPI ranks, for T>0 all bands
    real(q), intent(inout) :: nodes_in_dielectric_function !< position of innermst node in dielectric function
    real(q), intent(inout) :: wplasmon(:,:) !< plasmon frequency
    integer, intent(in) :: iu0
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_WAVEDER)
    integer :: ierr
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    character(LEN=256) :: spingroupname
    integer :: nb_totf, nkptsf, ispinf, nbands_cderf
    integer :: nb_tot, nkpts, ispin
    integer :: nk
    integer :: nkpts_non_zero
    character(len=MAX_LEN_GROUP) :: my_group
    logical :: lkpoints_shifted
    my_group = GRP_WAVEDER;     if (present(group))    my_group = group


    nb_tot=wdes%nb_tot
    nkpts=wdes%nkpts
    ispin=wdes%ispin

    nkpts_non_zero=0
    do nk=1,wdes%nkpts
       if (wdes%wtkpt(nk)==0) exit
       nkpts_non_zero=nkpts_non_zero+1
    enddo

    ! check if group is present:
    ierr= vh5_check_group_presence( fileid, wdes%comm, iu0, trim( my_group ) )
    if( ierr /= 0 ) return
    ! read from same wavefunctions group
    VH5_CHECK(vh5_group_open(fileid, trim(my_group), groupid, wdes%comm))
    ! read stuff
    ierr = vh5_read(groupid, "nodes_in_dielectric_function", nodes_in_dielectric_function, wdes%comm)
    if( ierr /=0 ) then
       if( iu0>=0 ) write(*,*)"vaspwave.h5 does not contain information about wavefunction derivatives."
       if( iu0>=0 ) write(*,*)"Did you forget to set LOPTICS = .TRUE. when generating this file?"
       ierr = -1
       goto 100
    endif

    ierr = vh5_read(groupid, "wplasmon", wplasmon, wdes%comm )
    if( ierr /=0 ) then
       if( iu0>=0 ) write(*,*)"vaspwave.h5 does not contain plasmon frequency."
       if( iu0>=0 ) write(*,*)"Did you forget to set LOPTICS = .TRUE. when generating this file?"
       ierr = 5
       goto 100
    endif

    VH5_CHECK(vh5_read(groupid, "nb_tot", nb_totf,wdes%comm ) )
    IF( nb_totf < nb_tot ) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined with too few bands: "//&
          str( [ nb_totf, nb_tot])
       ierr = 1
       goto 100
    ENDIF

    if( nb_totf /= nb_tot ) then
       if (lbse) then
          CALL vtutor%write(isAlert, BSEoptics, newArgument(ival=[nb_totf, nb_totf, min(nb_totf,nb_tot)]))
       endif
    endif

    VH5_CHECK(vh5_read(groupid, "nbands_cder", nbands_cder,wdes%comm ) )
    IF( nbands_cder < 1  ) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined with too few bands(2): "//&
          str( [ nbands_cder])
       ierr = 2
       goto 100
    ENDIF

    VH5_CHECK(vh5_read(groupid, "nkpts", nkptsf, wdes%comm ) )
    ! make sure we do not have shifted k-point grids
    if( nkptsf < nkpts .and. nkptsf /= nkpts_non_zero ) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined with too few k-points: "//&
          str( [ nkptsf, nkpts])
       ierr = 3
       goto 100
    ENDIF

    VH5_CHECK(vh5_read(groupid, "ispin", ispinf, wdes%comm ) )
    if( ispinf < ispin ) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined without spin! "
       ierr = 4
       goto 100
    endif

100 continue
    ! close/deallocate
    VH5_CHECK(vh5_group_close(groupid, wdes%comm))
  end function vh5_read_waveder_header

  !> @brief reads contents of waveder from h5 file
  function vh5_read_waveder(fileid, wdes, nbands_cder, nodes_in_dielectric_function, wplasmon, &
        cder_between_states, iu0, group) result (ierr)
    use string, only: str
    USE base
    USE mkpoints
    USE full_kpoints , ONLY:IKPT_SELECT
    USE tutor
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    type(wavedes), intent(in) :: wdes !< wave function descriptor
    integer, intent(inout) :: nbands_cder !< number of occupied bands rounded to MPI ranks, for T>0 all bands
    real(q), intent(inout) :: nodes_in_dielectric_function !< position of innermst node in dielectric function
    real(q), intent(inout) :: wplasmon(:,:) !< plasmon frequency
    GDEFS,intent(inout) :: cder_between_states(:,:,:,:,:)
    integer, intent(in) :: iu0
    integer :: ierr
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_WAVE)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, k_groupid, spin_groupid
    character(LEN=256) :: spingroupname, kgroupname
    integer :: isp, k, kf
    integer :: nb_totf, nkptsf, ispinf, nbands_cderf
    integer :: nb_tot, nkpts, ispin
    integer :: kfound
    integer :: nkpts_non_zero
    integer :: n, m, i
    GDEFS,allocatable:: cbuf(:,:,:)
    character(len=MAX_LEN_GROUP) :: my_group
    integer :: node_me
    integer :: ionode
    !-------------------------------------------------------------------------------------
    node_me=0
    ionode=0
#ifdef MPI
    node_me=wdes%comm%node_me
    ionode =wdes%comm%ionode
#endif
    my_group = GRP_WAVEDER;     if (present(group))    my_group = group


    ! check if group is present:
    ierr= vh5_check_group_presence( fileid, wdes%comm, iu0, trim( my_group ) )
    if( ierr /= 0 ) return
    nb_tot= size( cder_between_states, 1 )
    nbands_cder = size( cder_between_states, 2 )
    nkpts = size( cder_between_states, 3 )
    ispin = size( cder_between_states, 4 )
    ! error code
    ierr=0
             
    nkpts_non_zero=0
    do k=1,wdes%nkpts
       if (wdes%wtkpt(k)==0) exit
       nkpts_non_zero=nkpts_non_zero+1
    enddo

    ! read from same wavefunctions group
    VH5_CHECK(vh5_group_open(fileid, trim(my_group), groupid, wdes%comm))

    ! read stuff
    ierr = vh5_read(groupid, "wplasmon", wplasmon, wdes%comm )
    if( ierr /=0 ) then
       if( iu0>=0 ) write(*,*)"vaspwave.h5 does not contain wavefunction derivative (waveder) information."
       if( iu0>=0 ) write(*,*)"Did you forget to set LOPTICS = .TRUE. when generating this file?"
       ierr = -1
       VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    endif

    VH5_CHECK(vh5_read(groupid, "nb_tot", nb_totf, wdes%comm ) )
    IF( nb_totf < nb_tot ) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined with too few bands: "//&
          str( [ nb_totf, nb_tot])
       ierr = 1
       VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    ENDIF
    VH5_CHECK(vh5_read(groupid, "nbands_cder", nbands_cderf, wdes%comm ) )
    IF( nbands_cderf < nbands_cder ) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined with too few virtual bands: "//&
          str( [ nbands_cderf, nbands_cder])
       ierr = 2
       VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    ENDIF
    VH5_CHECK(vh5_read(groupid, "nkpts", nkptsf, wdes%comm ) )
    
    if( nkptsf < nkpts .and. nkptsf /= nkpts_non_zero) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined with too few k-points: "//&
          str( [ nkptsf, nkpts])
       ierr = 3
       VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    endif
    VH5_CHECK(vh5_read(groupid, "ispin", ispinf, wdes%comm ) )
    if( ispinf < ispin) then
       if( iu0>=0 ) write(*,*) "optical transition elements in vaspwave.h5 have been determined without spin! "
       ierr = 4
       VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    endif
    VH5_CHECK(vh5_read(groupid, "nodes_in_dielectric_function", nodes_in_dielectric_function, wdes%comm) )

    !> print info about downsampling

    if( iu0>=0 ) write(*,*)"reading wavefunction derivatives from vaspwave.h5"
    i=0
    if( nkptsf > nkpts ) then
       i=i+1
       if( iu0>=0 ) write(*,*) "vaspwave.h5 stores derivatives for denser k-point grid than presently selected."
    endif
    if( nb_totf> nb_tot ) then
       i=i+1
       if( iu0>=0 ) write(*,*) "vaspwave.h5 stores derivatives for more bands than presently selected."
    endif
    if( ispinf > ispin ) then
       i=i+1
       if( iu0>=0 ) write(*,*) "vaspwave.h5 stores derivatives for more spin channels than presently selected."
    endif

    if( i > 0 ) then
       if( iu0>=0 ) write(*,'(a)', advance='no') " trying downsample derivative data..."
    endif

    io_begin
    ! read optical transition elements
    allocate( cbuf( nb_totf, nbands_cderf, 3))
    io_end

    ierr =0
    spin: do isp=1,ispinf
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       VH5_CHECK(vh5_group_open(groupid, trim(spingroupname), spin_groupid, wdes%comm))
       kfound=wdes%nkdim - nkpts_non_zero
       kpoints_loop: do kf=1, nkptsf
          if ( allocated( ikpt_select ) ) then
             if ( ikpt_select(kf) < 1 ) cycle kpoints_loop
             k=ikpt_select(kf)
          else
             k=kf
          endif
          write(kgroupname, '(2A)') "kpoint_"//trim(str(kf))
          cder_between_states(:,:,k,isp,:) = 0

          VH5_CHECK(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid, wdes%comm))
         ! cder_between_states(1:nb_tot,1:nbands_cder,1:nkpts,1:ispin,1:3)
          io_begin
          ierr = vh5_read(k_groupid, "cder_between_states", cbuf)
          do i = 1, 3
             do m = 1, min( nbands_cderf, nbands_cder )
                do n = 1, min( nb_totf, nb_tot)
                   cder_between_states(n,m,k,isp,i) = cbuf(n,m,i)
                enddo
             enddo
          enddo
          io_end
          CALLMPI( M_sum_i( wdes%COMM, ierr, 1 ) )
          if (ierr/=0) return

          !> increase k-point number
          kfound=kfound+1

          VH5_CHECK(vh5_group_close(k_groupid, wdes%comm))
       enddo kpoints_loop
       VH5_CHECK(vh5_group_close(spin_groupid, wdes%comm))

       if( kfound < nkpts ) then
          if( iu0>=0 ) write(*,*) "failed"
          if( iu0>=0 ) write(*,*) "not all k-points of optical transition elements in vaspwave.h5 have been read: &
             &"//str( [ kfound, nkpts ])
          ierr = 5
          return
       endif

    enddo spin
    io_begin
    deallocate( cbuf )
    io_end

    CALLMPI( M_bcast_i(wdes%COMM, ierr, 1 ) )

    !> synchronize data
#ifdef gammareal
    CALLMPI( M_bcast_s(wdes%COMM, cder_between_states(1,1,1,1,1), size(cder_between_states) ))
#else
    CALLMPI( M_bcast_s(wdes%COMM, cder_between_states(1,1,1,1,1), size(cder_between_states)*2 ))
#endif
    ! close/deallocate
    VH5_CHECK(vh5_group_close(groupid, wdes%comm))
    if( ierr == 0 ) then
       if( iu0>=0 ) write(*,*)" succeeded."
    else
       if( iu0>=0 ) write(*,*)" FAILED!"
       cder_between_states = 0
       call vtutor%alert("Please regenerate the wavefunction derivative with LOPTICS, if you intended to use it for&
           & the long-wave limit of the polarizability. This job will continue without this information." )
    endif
  end function vh5_read_waveder

  !> @brief write charge (single spin or all spins)
  !>
  !> @param fileid hdf5 file handle
  !> @param subgroup subgroup to store the position data in
  !> @param grid grid_3d grid description type (from mgrid.F)
  !> @param chargetot charge array. shape is (ispin, nodes)
  subroutine vh5_write_charge(fileid, subgroup, grid, chargetot)
    use mgrid_struct_def, only: grid_3d
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    type(grid_3d) :: grid
    complex(q) :: chargetot(:,:)
    !-------------------------------------------------------------------------------------
    complex(q), allocatable :: work(:,:)
    integer spin
    !-------------------------------------------------------------------------------------
    ! bring density to real space
    allocate(work, source=chargetot)
    do spin = 1, size(work, 2)
      call FFT3D(work(:,spin), grid, 1)
    end do
    !
    call vh5_write_grid_quantity(fileid, subgroup, "charge", grid, work)
    !
  end subroutine vh5_write_charge

  !> ---------------------------- SUBROUTINE VH5 WRITE POTENTIAL---------------------------
  !> @brief write local potential on the plane wave grid
  !>  for spin-unpolarized calculations it writes one dataset
  !>  for spin-polarized calculations (ISPIN=2) it writes spin up and spin down potential
  !>  for noncollinear calculations it writes the spinor representation (4 datasets)
  !> -_------------------------------------------------------------------------------------
  subroutine vh5_write_potential(fileid, key, gridc, cvtot, group, subgroup)
    use mgrid_struct_def, only: grid_3d
    integer(HID_T), intent(in)             :: fileid    !< id of h5 file
    character(len=*), intent(in) :: key                 !< key used to write the potential to the file
    character(len=*), intent(in), optional :: group     !< group in h5 file, default: GRP_RESULTS
    character(len=*), intent(in), optional :: subgroup  !< subgroup in h5 file, default: GRP_POTENTIAL
    type(grid_3d), intent(in) :: gridc                  !< real space grid
    complex(q), intent(in)    :: cvtot(:,:)             !< potential computed by POTLOK
    !-------------------------------------------------------------------------------------
    character(len=:), allocatable :: my_group
    !-------------------------------------------------------------------------------------
    if (present(group)) then
      my_group = group
    else
      my_group = GRP_RESULTS
    end if
    if (present(subgroup)) then
      my_group = my_group // "/" // subgroup
    else
      my_group = my_group // "/" // GRP_POTENTIAL
    end if
    !
    call vh5_write_grid_quantity(fileid, my_group, key, gridc, cvtot)
    !
  end subroutine vh5_write_potential

  !*******************************************************************
  !> write a quantity on a grid to the hdf5 file
  !>
  !> this function assumes the following conditions are fulfilled
  !> - the quantity is given in real space and is real (imaginary part is discarded
  !> - the grid is consistent with other quantities in the same group
  !*******************************************************************
  subroutine vh5_write_grid_quantity(fileid, group, key, grid, quantity)
    use prec, only: q
    use mgrid, only: grid_3d
    use string, only: str
    use tutor, only: vtutor
    integer(HID_T), intent(in) ::  fileid   !< id of h5 file
    character(len=*), intent(in) :: group   !< name of the group in the HDF5 file
    character(len=*), intent(in) :: key     !< name of the dataset in the HDF5 file
    type(grid_3d), intent(in) :: grid       !< grid on which the quantity is defined
    complex(q), intent(in) :: quantity(:,:) !< quantity which is written to the HDF5 file
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(HSIZE_T) :: start(4), count_(4), dimensions(4)
    real(q), allocatable :: work(:)
    integer size_work, status_, spin, plane
    !-------------------------------------------------------------------------------------
    ! check for consistency of grid data
    if (grid%nplwv /= grid%ngx * grid%ngy * grid%ngz) then
       call vtutor%bug('internal ERROR: vh5_write_grid_quantity gridc%nplwv = ' // &
         str(grid%nplwv) // 'is not compatible with NGX = ' // str(grid%ngx) // &
         ', NGY = ' // str(grid%ngy) // ', NGZ = ' // str(grid%ngz), __FILE__, __LINE__)
    end if
    !
    ! allocate work arrays
    size_work = grid%ngx * grid%ngy
    allocate(work(size_work), stat=status_)
    if (status_ /= 0) then
       write (*,*) "WARNING: insufficient memory to write " // key // " to hdf5 file"
       return
    end if
    !
    VH5_CHECK(vh5_group_open_or_create(fileid, group, groupid))
    !
    ! write grid dimensions
    VH5_CHECK(vh5_write(groupid, "grid", [grid%ngx, grid%ngy, grid%ngz]))
    !
    ! write potential
    dimensions = [grid%ngx, grid%ngy, grid%ngz, size(quantity, 2)]
    VH5_CHECK(vh5_create_double_array_nd(groupid, key, 4, dimensions))
    do spin = 1, size(quantity, 2)
       ! write one plane of the grid quantity
       do plane = 1, grid%ngz
          call mrg_grid_rl_plane(grid, work, quantity(:,spin), plane)
          start = [1, 1, plane, spin]
          count_ = [grid%ngx, grid%ngy, 1, 1]
          VH5_CHECK(vh5_write_double_subarray_nd(groupid, key, 4, start, count_, work))
       end do
    end do
    !
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
    !
  end subroutine vh5_write_grid_quantity

  !> write self-energy of low-scaling job
  subroutine vh5_write_self_energy(fileid, w, w_prev, sigmau, sigmao, nkpts_irz, imag_grids, nbandsgw, io, subgroup)
    USE constant
    USE base
    USE wave
    USE vaspxml
    USE minimax_struct
    USE chi_glb, ONLY: map_self_energy
    integer(HID_T), intent(in)         :: fileid
    type (wavespin), intent(in)        :: w               !wave function
    type (wavespin), intent(in)        :: w_prev          ! quasi-particle previous iteration
    GDEF, intent(in)                   :: sigmao(:,:,:,:)  !diagonal of sin transformed gf in bloch domain
    GDEF, intent(in)                   :: sigmau(:,:,:,:)  !diagonal of cos transformed gf in bloch domain
    integer, intent(in)                :: nkpts_irz
    type(imag_grid_handle), intent(in) :: imag_grids        !time and frequency grid handle
    integer, intent(in)                :: nbandsgw        !number of bands to be written to file
    type(in_struct), intent(in)        :: io
    character(LEN=*), intent(in), optional :: subgroup
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid, mygroupid
    integer(SIZE_T) :: dims4(4), start4(4)
    integer(SIZE_T) :: dims3(3), start3(3)
    integer(SIZE_T) :: dims2(2), start2(2)
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup

    integer                :: i
    integer                :: nb,nk,isp
    complex(q), allocatable :: sigma_c(:,:,:,:)
    complex(q), allocatable :: eps(:,:,:)
    real(q), allocatable    :: occ(:,:,:)
    complex(q), allocatable :: complex_frequencies(:,:)
    !-------------------------------------------------------------------------------------
    my_group    = GRP_INTERMEDIATE
    my_subgroup = GRP_SELFENERGY
    !-------------------------------------------------------------
    ! open the groups
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    if( present( subgroup ) ) then
       VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), mygroupid))
       VH5_CHECK(vh5_group_open_or_create(mygroupid, trim(subgroup), subgroupid))
    else
       VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))
    endif

    allocate(sigma_c(imag_grids%nomega, nbandsgw, nkpts_irz, w%wdes%ispin))
    allocate(eps(nbandsgw, nkpts_irz, w%wdes%ispin))
    allocate(occ(nbandsgw, nkpts_irz, w%wdes%ispin))
    allocate(complex_frequencies(imag_grids%nomega, w%wdes%ispin))

    dims2 = [ imag_grids%nomega, w%wdes%ispin ]
    start2 = [1,1]

    dims3 = [ nbandsgw, nkpts_irz, w%wdes%ispin ]
    start3 = [1,1,1]

    dims4 = [ imag_grids%nomega, nbandsgw, nkpts_irz, w%wdes%ispin ]
    start4 = [1,1,1,1]

    VH5_CHECK(vh5_write(subgroupid, "dims", [nbandsgw, nkpts_irz, w%wdes%ispin]))

    do isp = 1, w%wdes%ispin
       do i = 1, imag_grids%nomega
          complex_frequencies(i,isp) = cmplx(w_prev%efermi(isp), imag_grids%fer_re(i), q )
       enddo
    enddo
    VH5_CHECK(vh5_write(subgroupid, "complex_frequencies", complex_frequencies ))

    ! dump one-electron energies (eigenvalues of HF hamiltonian used for basis representation)
    ! required to find the poles of the Green's function as follows
    ! eps + w + efermi - Sigma(w) = 0
    eps(1:nbandsgw,1:nkpts_irz, 1:w%wdes%ispin) = w%celtot(1:nbandsgw, 1:nkpts_irz, 1:w%wdes%ispin)
    VH5_CHECK(vh5_write(subgroupid, "hf_energies", eps ))

    ! dump previous quasiparticle energies
    eps(1:nbandsgw,1:nkpts_irz, 1:w%wdes%ispin) = w_prev%celtot(1:nbandsgw, 1:nkpts_irz, 1:w%wdes%ispin)
    VH5_CHECK(vh5_write(subgroupid, "prev_qp_energies", eps ))

    ! dump previous quasiparticle occupations
    occ(1:nbandsgw,1:nkpts_irz, 1:w%wdes%ispin) = w_prev%fertot(1:nbandsgw, 1:nkpts_irz, 1:w%wdes%ispin)
    VH5_CHECK(vh5_write(subgroupid, "prev_qp_occupations", occ ))

    ! map from time representation of correlation part of self-energy to frequency domain
    do isp = 1, w%wdes%ispin
       do nk = 1, nkpts_irz
           call map_self_energy( sigmau(:,:,nk,isp), sigmao(:,:,nk,isp), imag_grids, nbandsgw, sigma_c(:,:,nk,isp) )
       enddo
    enddo
    VH5_CHECK(vh5_write(subgroupid, "sigma_c", sigma_c))

    deallocate( complex_frequencies )
    deallocate( occ )
    deallocate( eps )
    deallocate( sigma_c )

    if( present( subgroup ) ) then
       VH5_CHECK(vh5_group_close_writing(subgroupid))
       VH5_CHECK(vh5_group_close_writing(mygroupid))
    else
       VH5_CHECK(vh5_group_close_writing(subgroupid))
    endif
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_self_energy

  !*******************************************************************
  !> initialize the array in the hdf5 file for writing partial charges
  !>
  !> this function creates a 6D array where partial charges are written
  !> to. Two 1D integer arrays for band and kpoint numbers are also
  !> written to the file. If they just contain a 0, this means that
  !> there will be an average over all selected bands or kpoints.
  !*******************************************************************
  subroutine vh5_init_parcharge_output(fileid, grid, nspins, bands, kpoints, group, subgroup)
    use mgrid, only : grid_3d
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(grid_3d), intent(in) :: grid      ! grid for potential / charge
    integer, intent(in) :: nspins, bands(:), kpoints(:)  ! nr. of spin channels and list of bands and kpoints
    character(len=*), intent(in), optional :: group !< group name (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PARCHG)
    !-------------------------------------------------------------------------------------
    character(len=20), parameter :: dataset_name = "parchg"
    integer, parameter :: rank = 6
    integer(HID_T) :: groupid, subgroupid
    integer(HSIZE_T) :: dimensions(rank)
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup

    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = SUBGRP_PARCHG; if (present(subgroup)) my_subgroup = trim(subgroup)

    ! open the groups
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    ! write grid data
    VH5_CHECK(vh5_write(subgroupid, "grid", [grid%ngx, grid%ngy, grid%ngz]))

    ! write lists of bands and kpoints
    VH5_CHECK(vh5_write(subgroupid, "bands", bands))
    VH5_CHECK(vh5_write(subgroupid, "kpoints", kpoints))

    ! prepare the array dimensions and the array
    dimensions = [grid%ngx, grid%ngy, grid%ngz, nspins, size(bands), size(kpoints)]
    VH5_CHECK(vh5_create_double_array_nd(subgroupid, dataset_name, rank, dimensions))

    ! close the groups again
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
 end subroutine vh5_init_parcharge_output

 !*******************************************************************
 !> write the partial charge density to the hdf5 file
 !>
 !> This function writes a parchg array to the hdf5 file, eiter for
 !> a single spin channel or for two of them. The array is written
 !> to the given band and kpoint index.
 !*******************************************************************

 subroutine vh5_write_to_parcharge_array(fileid, grid, chgtot, iband, ikpt, group, subgroup)
    use mgrid, only : grid_3d
    use string, only: str
    USE tutor, only: vtutor
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(grid_3d), intent(in) :: grid      ! grid for potential / charge
    complex(q), intent(in) :: chgtot(:,:) !<The (partial) charge density
    integer, intent(in) :: iband !< counting band index
    integer, intent(in) :: ikpt !< counting kpoint index
    character(len=*), intent(in), optional :: group !< group name (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PARCHG)
    !-------------------------------------------------------------------------------------
    complex(q), allocatable :: work_fft(:,:) !< work array for fft3d
    real(q), allocatable :: work(:) !< work array for real space charge writing
    character(len=20), parameter :: dataset_name = "parchg"
    integer, parameter :: rank = 6
    integer(HSIZE_T) :: start(6), count_(6)
    integer(HID_T) :: groupid, subgroupid, kpts_id, band_id
    integer :: spin, size_work, status_, plane
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = SUBGRP_PARCHG; if (present(subgroup)) my_subgroup = trim(subgroup)

    ! check for consistency of grid data
    if (grid%nplwv /= grid%ngx * grid%ngy * grid%ngz) then
       call vtutor%bug('internal ERROR: vh5_write_to_parcharge_array gridc%nplwv = ' // &
       str(grid%nplwv) // 'is not compatible with NGX = ' // str(grid%ngx) // &
       ', NGY = ' // str(grid%ngy) // ', NGZ = ' // str(grid%ngz), __FILE__, __LINE__)
    end if

    allocate(work_fft, source=chgtot)
    do spin = 1, size(work_fft, 2)
       call FFT3D(work_fft(:,spin), grid, 1)
    end do

    !
    ! allocate work arrays
    size_work = grid%ngx * grid%ngy
    allocate(work(size_work), stat=status_)
    if (status_ /= 0) then
       write (*,*) "WARNING: insufficient memory to write " // dataset_name // " to hdf5 file"
       return
    end if

    ! open the groups
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    ! write the partial charge density
    do spin = 1, size(work_fft, 2)
       ! write one plane of the grid quantity
       do plane = 1, grid%ngz
          call mrg_grid_rl_plane(grid, work, work_fft(:,spin), plane)
          start = [1, 1, plane, spin, iband, ikpt]
          count_ = [grid%ngx, grid%ngy, 1, 1, 1, 1]
          VH5_CHECK(vh5_write_double_subarray_nd(subgroupid, dataset_name, rank, start, count_, work))
       end do
    end do

    ! close the groups again
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
 end subroutine vh5_write_to_parcharge_array

  subroutine vh5_write_average_potential(fileid, latt_cur, direction, avg_potential, vacuum_potential, group, subgroup)
    use prec
    use LATTICE
    integer(HID_T), intent(in) :: fileid
    integer :: direction
    type(latt)      :: latt_cur
    real(q), intent(in) ::  avg_potential(:), vacuum_potential(2)
    character(len=*), intent(in), optional :: group     !< group in h5 file, default: GRP_RESULTS
    character(len=*), intent(in), optional :: subgroup  !< subgroup in h5 file, default: GRP_POTENTIAL
    integer(HID_T)       :: groupid, subgroupid
    character(len=20)    :: my_group, my_subgroup
    real(q), allocatable :: distance(:)
    integer :: i

    my_group = GRP_RESULTS;      if (present(group))    my_group = trim(group)
    my_subgroup = GRP_POTENTIAL; if (present(subgroup)) my_subgroup = trim(subgroup)

    ! Generate grid information for the averaged potential
    allocate(distance(size(avg_potential,1)))
    do i = 1, size(avg_potential,1)
      distance(i) = real(i,q) / size(avg_potential,1,q) * latt_cur%ANORM(direction)
    end do

    VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, my_subgroup, subgroupid))
    VH5_CHECK(vh5_write(subgroupid, "vacuum_potential", vacuum_potential))
    VH5_CHECK(vh5_write(subgroupid, "average_potential_along_IDIPOL", avg_potential))
    VH5_CHECK(vh5_write(subgroupid, "distance_along_IDIPOL", distance))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_average_potential

  !*******************************************************************
  !> write the PAW occupancies to a hdf5 file
  !*******************************************************************
  subroutine vh5_write_paw_occupancies(fileid, subgroup, p, t_info, loverl, rholm_store, comm, ispin)
    use pseudo
    use poscar
    use wave
    use constant
    use paw

    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type (type_info)   t_info
    type (potcar),target::  p(t_info%ntyp)
    logical  loverl          ! overlap matrix used ?
    real(q) rholm_store(:)   ! storage for the channel occupancies
    type(communic) :: comm
    integer :: ispin

  ! local variables
    integer(HID_T) :: unitgroupid
    type (potcar),pointer :: pp
    integer nt, ni, i
    integer ibase, iadd, lymax, lmmax
    integer, allocatable::  nelements(:)
    integer, external :: maxl_aug
    character(LEN=40) :: datasetname

    real(q), allocatable::  buffer(:,:)

    allocate(nelements(t_info%nions))

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid))
!=======================================================================
! quick return if possible
!=======================================================================
    if (.not.loverl .or. mimic_us ) return

    lymax =maxl_aug(t_info%ntyp,p)
    lmmax=(lymax+1)**2

    allocate( buffer(lmmax*lmmax, t_info%nions))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
    ibase=1

    ion: do ni=1,t_info%nions
       buffer(:,ni)=0
       iadd  =0

       nt=t_info%ityp(ni)
       pp=> p(nt)

       nelements(ni) = 0
       if (do_local(ni)) then
          call retrieve_rholm( buffer(:,ni), rholm_store(ibase:), &
               metric(ibase:), iadd, pp, .true., nelements(ni))

          if (nelements(ni) > lmmax*lmmax) then
             call vtutor%bug('internal error: wrt_rho_paw running out of buffer', __FILE__, __LINE__)
          endif
          ibase=ibase+iadd
       endif

       CALLMPI( m_sum_i(comm, nelements(ni), 1))
       CALLMPI( m_sum_d(comm, buffer(:,ni), nelements(ni)))

    enddo ion
    write(datasetname,'(A,I0.2)') "aug_occupancies_ions_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, nelements))
    write(datasetname,'(A,I0.2)') "aug_occupancies_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, buffer))

#ifdef rholm_complex
    ibase=1

    again: do ni=1,t_info%nions
       buffer(:,ni)=0
       iadd  =0

       nt=t_info%ityp(ni)
       pp=> p(nt)

       nelements =0
       if (do_local(ni)) then
          call retrieve_rholm( buffer(:,ni), rholm_store(ibase+ubound(rholm_store,1)/2:), &
               metric(ibase+ubound(rholm_store,1)/2:), iadd, pp, .true., nelements(ni))

          if (nelements(ni) > lmmax*lmmax) then
             call vtutor%bug('internal error: wrt_rho_paw running out of buffer', __FILE__, __LINE__)
          endif

          ibase=ibase+iadd
       endif

       CALLMPI( m_sum_i(comm, nelements(ni), 1))
       CALLMPI( m_sum_d(comm, buffer(:,ni), nelements(ni)))

    enddo again
    write(datasetname,'(A,I0.2)') "imag_aug_occupancies_ions_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, nelements))
    write(datasetname,'(A,I0.2)') "imag_aug_occupancies_ions_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, buffer))

#endif
    deallocate(nelements)
    deallocate(buffer)
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_paw_occupancies

  !> write PAW strength (Dij) and overlap (Qij) to the hdf5 file
  !> note that qij and dij are distributed over NCORE so we must gather and write
  subroutine vh5_write_paw_dij_qij(fileid,wdes,t_info,p,lmdim,cdij,cqij,group)
    use wave, only: ni_global
    use pseudo, only: potcar
    integer(HID_T), intent(in)   :: fileid
    type (wavedes), intent(in)   :: wdes
    type (type_info), intent(in) :: t_info
    type (potcar), intent(in) :: p(:)
    integer, intent(in) :: lmdim
    OVERLAP,intent(in) :: cdij(:,:,:,:)
    OVERLAP,intent(in) :: cqij(:,:,:,:)
    character(len=*),optional,intent(in) :: group
    !local variables
    integer(HID_T) :: groupid
    character(:),allocatable :: my_group
    integer :: i, j, isp, at, ni_local, ni, lmax
    complex(q),allocatable :: cdij_tot(:,:,:,:)
    complex(q),allocatable :: cqij_tot(:,:,:,:)
    integer,allocatable :: lps(:,:)

    ! collect the data among the mpi ranks
    allocate(cdij_tot(lmdim, lmdim, t_info%nions, wdes%ncdij))
    allocate(cqij_tot(lmdim, lmdim, t_info%nions, wdes%ncdij))
    cdij_tot = cmplx(0.0_q,0.0_q,q)
    cqij_tot = cmplx(0.0_q,0.0_q,q)
    do isp = 1, wdes%ncdij
       do ni_local = 1, wdes%nions
          ni = ni_global(ni_local, wdes%comm_inb)
          cqij_tot(:, :, ni, isp) = cmplx(cqij(:, :, ni_local, isp),kind=q)
          cdij_tot(:, :, ni, isp) = cmplx(cdij(:, :, ni_local, isp),kind=q)
       enddo
    enddo
    CALLMPI( M_sum_z(wdes%comm_inb, cdij_tot, size(cdij_tot)) )
    CALLMPI( M_sum_z(wdes%comm_inb, cqij_tot, size(cqij_tot)) )

    allocate(lps(maxval(p(:)%lmax),size(p)))
    lps(:,:)=-1
    do i=1,size(p)
        lmax = p(i)%lmax
        lps(:lmax,i) = p(i)%lps(:lmax)
    enddo

    ! write the data to hdf5 file
    my_group = GRP_RESULTS//'/'//GRP_PAW
    if (present(group)) my_group = trim(my_group)
    VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))

    ! write data
    VH5_CHECK(vh5_write(groupid,'dij',cdij_tot))
    VH5_CHECK(vh5_write(groupid,'qij',cqij_tot))
    VH5_CHECK(vh5_write(groupid,'lmdim',lmdim))
    VH5_CHECK(vh5_write(groupid,'ldim',p(:)%ldim))
    VH5_CHECK(vh5_write(groupid,'lmax',p(:)%lmax))
    VH5_CHECK(vh5_write(groupid,'lreal',p(:)%lreal))
    VH5_CHECK(vh5_write(groupid,'lps',lps))

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_paw_dij_qij

  subroutine vh5_write_atommom(fileid, subgroup, t_info)
    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type (type_info)   t_info

    integer(HID_T) :: unitgroupid

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid))
    VH5_CHECK(vh5_write(unitgroupid, "atomic_moments", t_info%atomom))
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))

  end subroutine vh5_write_atommom

  subroutine vh5_write_dielectric_vector(fileid, subgroup, wavevector, prefix)
    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    real(q), intent(in)             :: wavevector(:)
    character(LEN=*), intent(in) :: prefix     !< select "wave" or "qpoint"

    integer(HID_T) :: unitgroupid

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid))
    VH5_CHECK(vh5_write(unitgroupid, prefix // "_vector", wavevector))
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))

  end subroutine vh5_write_dielectric_vector


  !> @brief write wavefunctions
  !>
  !> @param fileid hdf5 file handle
  !> @param io io handle
  !> @param wave_des wave function descriptor (see wave.f)
  !> @param wave_spin wave functions (see wave.f)
  !> @param latt_ini lattice definition
  !> @param efermi fermi energy
  !> @param nbands_dump
  subroutine vh5_write_wavefunctions(fileid, io, wave_des, wave_spin, latt_ini, efermi, nbands_dump)
    USE base
    USE string
    USE wave
    USE lattice
    USE main_mpi
    USE tutor, only: vtutor

    integer(HID_T), intent(in)       :: fileid
    type(in_struct) :: io
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(latt)      :: latt_ini
    real (q)        :: efermi
    integer, optional :: nbands_dump
    !----------------------------------------------------------------------
    ! local variables
    integer(HID_T) :: groupid, spin_groupid, k_groupid
    integer :: ierr, i, j, k, npl, irec, ireclw_old, isp
    integer(HSIZE_T) :: dimensions(2), start(2), count(2),dim_nb_tot(1)
!    integer :: node_me, ionode
    integer npl_tot, iu0
    integer nb_tot_old
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl
    logical, pointer :: lconjg(:)=>null()
    integer, pointer :: index(:)=>null()
    character(LEN=256) :: spingroupname, kgroupname
    !----------------------------------------------------------------------
    ! local work arrays
    type(wavedes1) wdes1
    ! FIXME CW is replaced by CPTFWP via macro in symbol.inc. still need to find out
    ! why, and what the exact effects are
    complex(q), allocatable :: CW(:)
    real(q):: rdum_old, rispin_old, rnkpts_old, rnb_tot_old, efermi_old, enmax_old, a_old(3,3)
    real(q):: rnkpx_old, rnkpy_old, rnkpz_old
    real(q), parameter :: tiny=1e-6_q
    integer :: node_me =1
    integer :: ionode
    character(len=255) :: msg
    !----------------------------------------------------------------------
    !
     node_me=0
     ionode=0
#ifdef MPI
     node_me=wave_des%comm%node_me
     ionode =wave_des%comm%ionode
#endif

    ! backup old number of bands
    nb_tot_old = wave_des%nb_tot
    if ( present( nbands_dump ) ) then
       if( nbands_dump > 1 ) then
          wave_des%nb_tot = nbands_dump
       else
          call vtutor%bug(" silly number of bands passed "//str(nbands_dump),__FILE__, __LINE__ )
       endif
    endif
    !
    ! report to  stdout
    !

    if (io%iu6>=0) write(io%iu6,*)'writing wavefunctions to vaspwave.h5'
    if (io%iu0>=0) write(io%iu0,*)'writing wavefunctions to vaspwave.h5'
    !
    ! set up constants
    !
    npl_tot = maxval(wave_des%nplwkp_tot)
#ifdef WAVECAR_double
    io%ireclw=max(npl_tot,7)*io%icmplx
    ireclw_old=max(max(npl_tot,6),((wave_des%nb_tot*3+1)/2+2))*io%icmplx
#else
    io%ireclw=max((npl_tot+1)/2,7)*io%icmplx
    ireclw_old=max(max((npl_tot+1)/2,6),((wave_des%nb_tot*3+1)/2+2))*io%icmplx
#endif
    !
    ! allocate work arrays
    !
    allocate(CW(npl_tot))

    ! in order to increase exchangeability of WAVECAR files across IEEE platforms
    ! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
    rispin=wave_des%ispin
    rdum  =io%ireclw
    rnkpts =wave_des%nkpts
    rnb_tot=wave_des%nb_tot
    dim_nb_tot = [wave_des%nb_tot]

    !
    ! write constants
    !
    VH5_CHECK(vh5_group_open_or_create(fileid, GRP_WAVE, groupid))

    rdum_old = rdum
    rispin_old = rispin
    rnkpts_old = rnkpts
    rnb_tot_old = rnb_tot
    enmax_old = wave_des%enmax
    efermi_old = efermi
    a_old = latt_ini%a
    ! read old entries
    if( vaspwave_found ) then
       read_block: block
          ierr = vh5_read( groupid, "rdum" , rdum_old, wave_des%COMM );      if (ierr /= 0) exit read_block
          ierr = vh5_read( groupid, "rispin" , rispin_old, wave_des%COMM );  if (ierr /= 0) exit read_block
          ierr = vh5_read( groupid, "rnkpts", rnkpts_old, wave_des%COMM );   if (ierr /= 0) exit read_block
          ierr = vh5_read( groupid, "rnb_tot", rnb_tot_old, wave_des%COMM ); if (ierr /= 0) exit read_block
          ierr = vh5_read( groupid, "efermi", efermi_old, wave_des%COMM );if (ierr /= 0) exit read_block
          ierr = vh5_read( groupid, "enmax", enmax_old, wave_des%COMM );  if (ierr /= 0) exit read_block
          ierr = vh5_read( groupid, "amat", a_old, wave_des%COMM );     if (ierr /= 0) exit read_block
       end block read_block

       ! in case all of these entries are present in file and have been read in
       if ( ierr == 0 ) then
          irec=0
          ! check if they have changed
          msg='The following parameters are overwritten:'//new_line('a')
          if ( abs( rdum - rdum_old ) > tiny ) then
             ierr = ierr +1
             irec= irec + 1
             msg=trim(msg)//"   Plane-wave coefficients: "//str( int(rdum) )//" <- "//str(int(rdum_old))//new_line('a')
          endif
          if( abs( rispin - rispin_old ) > tiny ) then
             ierr = ierr +1
             irec= irec + 1
             msg=trim(msg)//"   ISPIN: "//str( int(rispin) )//" <- "//str(int(rispin_old))//new_line('a')
          endif
          if( abs( rnkpts - rnkpts_old ) > tiny) then
             ierr = ierr +1
             irec= irec + 1
             msg=trim(msg)//"   k-points: "//str( int(rnkpts) )//" <- "//str(int(rnkpts_old))//new_line('a')
          endif
          if( abs( rnb_tot - rnb_tot_old ) > tiny ) then
             ierr = ierr +1
             irec= irec + 1
             msg=trim(msg)//"   NBANDS: "//str( int(rnb_tot) )//" <- "//str(int(rnb_tot_old))//new_line('a')
          endif
          if( abs( efermi - efermi_old ) > tiny ) then
             irec= irec + 1
             msg=trim(msg)//"   E-fermi: "//str( efermi )//" <- "//str(efermi_old)//", diff: "//str(abs(efermi- efermi_old))//new_line('a')
          endif
          if( abs( wave_des%enmax - enmax_old ) > tiny ) then
             irec= irec + 1
             msg=trim(msg)//"   ENCUT: "//str( wave_des%enmax )//" <- "//str(enmax_old)//new_line('a')
          endif
          if( norm2( latt_ini%a - a_old ) > tiny ) then
             irec= irec + 1
             msg=trim(msg)//"   Lattice vectors |A-Anew|: "//str( norm2( latt_ini%a-a_old) ) //new_line('a')
          endif
          if( irec > 0 .and. io%iu0>=0 ) write(*,*) trim(msg)
          if( irec > 0 .and. io%iu6>=0 ) write(io%iu6,*) trim( msg )
          if( ierr > 0 ) then
             call vtutor%alert( "vaspwave.h5 might store redundant data. To reduce &
                &its filesize a call to h5repack could be appropriate!" )
          endif
       endif
    endif

    ! check if vaspwave group already exists
    if (vh5_exists(fileid, GRP_WAVE, wave_des%comm)) then
       VH5_CHECK(vh5_group_close_writing(groupid))
       VH5_CHECK(vh5_delete(fileid, GRP_WAVE))
       VH5_CHECK(vh5_group_open_or_create(fileid, GRP_WAVE, groupid))
    endif

    VH5_CHECK(vh5_write(groupid, "rdum", rdum))
    VH5_CHECK(vh5_write(groupid, "rispin", rispin))
    !
    ! write some more parameters
    !
    VH5_CHECK(vh5_write(groupid, "rnkpts", rnkpts))
    VH5_CHECK(vh5_write(groupid, "rnb_tot", rnb_tot))
    VH5_CHECK(vh5_write(groupid, "efermi", efermi))
    VH5_CHECK(vh5_write(groupid, "enmax", wave_des%enmax))
    VH5_CHECK(vh5_write(groupid, "amat", latt_ini%a))
    !
    !
    !
    !
    ! loop over spin, kpoints, bands
    !
    spin: DO ISP=1,wave_des%ISPIN
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       VH5_CHECK(vh5_group_open_or_create(groupid, trim(spingroupname), spin_groupid))
       kpoints: DO K=1,wave_des%NKPTS
          write(kgroupname, '(2A)') "kpoint_"//trim(str(k))
          VH5_CHECK(vh5_group_open_or_create(spin_groupid, trim(kgroupname), k_groupid))
#ifdef MPI
          IF (wave_des%COMM_KINTER%NCPU > 1) THEN
             IF (wave_des%COMM_KINTER%NODE_ME /= 1) THEN
                IF (MOD(K-1,wave_des%COMM_KINTER%NCPU).NE.wave_des%COMM_KINTER%NODE_ME-1) THEN
                   CYCLE
                ELSE
                   CALLMPI( M_send_z(wave_des%COMM_KINTER, 1, wave_spin%CW(1,1,K,ISP), SIZE(wave_spin%CW,1)*SIZE(wave_spin%CW,2)) )
                END IF
             ELSE
                IF (MOD(K-1,wave_des%COMM_KINTER%NCPU).NE.wave_des%COMM_KINTER%NODE_ME-1) THEN
                   CALLMPI( M_recv_z(wave_des%COMM_KINTER, MOD(K-1,wave_des%COMM_KINTER%NCPU)+1, wave_spin%CW(1,1,K,ISP), SIZE(wave_spin%CW,1)*SIZE(wave_spin%CW,2)) )
                ENDIF
             ENDIF
          END IF
#endif
          CALL SETWDES(wave_des,WDES1,K)
#ifdef gammareal
          CALL GEN_INDEX_GAMMA(WDES1, latt_ini%B, wave_des%ENMAX, INDEX, LCONJG)
#endif
          npl=wave_des%nplwkp_tot(k)
          rnpl=npl
          ! write number of plane waves, k-point coordinates and all eigenvalues and
          ! occupation numbers for current k-point K
          ! write eigenvalues in real format
          ! io_begin
          VH5_CHECK(vh5_write(k_groupid, "num_planewaves", npl))
          VH5_CHECK(vh5_write(k_groupid, "vkpt", wave_des%vkpt(:,k)))
          VH5_CHECK(vh5_write(k_groupid, "celtot", wave_spin%celtot(1:wave_des%nb_tot,k,isp)))
          VH5_CHECK(vh5_write(k_groupid, "fertot", wave_spin%fertot(1:wave_des%nb_tot,k,isp)))
          ! io_end


          dimensions = [npl, wave_des%nb_tot]
#ifdef WAVECAR_double
          VH5_CHECK(vh5_create_double_complex_array_nd(k_groupid, "wave", 2, dimensions))
#else
          VH5_CHECK(vh5_create_complex_array_nd(k_groupid, "wave", 2, dimensions))
#endif

          DO J=1,wave_des%nb_tot
#ifdef gammareal
             CALL MRG_PW_BAND_GAMMA(WDES1, J, CW, wave_spin%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
             CALL MRG_PW_BAND(WDES1, J, CW, wave_spin%CW(1,1,K,ISP))
#endif
             start = [1,J]
             count = [npl,1]
             VH5_CHECK(vh5_write_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, CW))
          enddo

#ifdef gammareal
          call free_index_gamma(index, lconjg)
#endif
          VH5_CHECK(vh5_group_close_writing(k_groupid))
       enddo kpoints
       VH5_CHECK(vh5_group_close_writing(spin_groupid))
    enddo spin

    deallocate(CW)
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
    ! reset number of bands
    wave_des%nb_tot =  nb_tot_old
  end subroutine vh5_write_wavefunctions

  !> read force constants from hdf5 file
  subroutine vh5_read_force_constants(fileid,force_constants)
    use tutor, only: vtutor
    use string, only: str
    integer(HID_T),intent(in) :: fileid
    real(q),allocatable,intent(inout) :: force_constants(:,:)
    ! local variables
    integer :: dims(2)
    integer :: rank
    real(q) :: fc(1,1)
    integer(HID_T) :: ih5inresultsgroup_id
    integer :: ierr

    ierr = vh5_get_dimensions(fileid, '/results/linear_response/force_constants', dims, rank)
    if (ierr/=0) then
        call vtutor%error("Could not find '/results/linear_response/force_constants' in " // VASPIN)
    endif
    if (rank/=2) then
        call vtutor%error("The force_constants array in the HDF5 file has rank "//str(rank)//&
                          ", but rank 2 is expected.")
    endif
    ! check the dimensions
    if (any(shape(force_constants)/=dims)) then
        call vtutor%error("The dimensions of the force_constants array in the HDF5 file "//str(dims(1))//&
                          " are incompatible with the current VASP run "//str(size(force_constants,1))//&
                          ". Please check if the POSCAR is the same you used for the supercell calculation")
    endif
    VH5_CHECK(vh5_read(fileid, '/results/linear_response/force_constants', force_constants))
  end subroutine vh5_read_force_constants

  !> @brief read wave function header
  !>
  !> @param fileid hdf5 file handle
  !> @param wave_des wave function descriptor (see wave.f)
  !> @param latt_ini lattice definition
  !> @param latt_cur lattice definition
  !> @param enmaxi
  !> @param istart
  !> @param iu0
  function vh5_read_wavefunction_header(fileid, wave_des, latt_ini, latt_cur, enmaxi, istart, iu0, group ) result(ierr)
    USE base
    USE wave
    USE lattice
    USE main_mpi

    integer(HID_T), intent(in) :: fileid
    type(wavedes), intent(in) :: wave_des
    type(latt), intent(inout) :: latt_ini
    type(latt), intent(in) :: latt_cur
    real (q), intent(out) :: enmaxi
    integer, intent(inout) :: istart
    integer         :: iu0
    character(len=*),optional,intent(in) :: group
    integer :: ierr
    !-------------------------------------------------------------
    ! local variables
    integer(HID_T) :: groupid
    integer :: i,j
    logical :: ldiff
    real(q) :: rnb_tot, rnkpts
    integer :: nbandf, nkptsf
    character(len=MAX_LEN_GROUP) :: my_group
    my_group = GRP_WAVE;     if (present(group))    my_group = group
    !-------------------------------------------------------------

    ! check if group is present:
    ierr= vh5_check_group_presence( fileid, wave_des%comm, iu0, trim( my_group ) )
    if( ierr /= 0 ) then
       if( trim( my_group ) == GRP_WAVE ) then
          ! tell vasp to initialize wavefunction later on in main.F
          istart = 0
       endif
       return
    endif
    ! open group
    ierr = vh5_group_open(fileid, trim(my_group), groupid, wave_des%comm)
    if(ierr == 0) then
       if (iu0 >= 0) then
          write(iu0,*) 'found vaspwave.h5, reading the header of group: ', trim(my_group)
       endif

       VH5_CHECK(vh5_read(groupid, "rnkpts", rnkpts, wave_des%comm))
       VH5_CHECK(vh5_read(groupid, "rnb_tot", rnb_tot, wave_des%comm))
       VH5_CHECK(vh5_read(groupid, "enmax", enmaxi, wave_des%comm))
       VH5_CHECK(vh5_read(groupid, "amat", latt_ini%a, wave_des%comm))

       NKPTSF=NINT(rnkpts)
       NBANDF=NINT(rnb_tot)

       CALL LATTIC(LATT_INI)

       IF (ISTART==2 .AND. ENMAXI /= wave_des%ENMAX) THEN
          CALL vtutor%error('ERROR: ENMAX changed please set ISTART to 1')
       ENDIF

       IF (.NOT. (NBANDF==wave_des%NB_TOT .OR. (wave_des%NRSPINORS==2 .AND.  NBANDF*wave_des%NRSPINORS==wave_des%NB_TOT))) THEN
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'number of bands has changed, file:',NBANDF,' present:',wave_des%NB_TOT
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'trying to continue reading vaspwave.h5, but it might fail'
       ENDIF
       IF (NKPTSF/=wave_des%NKPTS) THEN
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'number of k-points has changed, file:',NKPTSF,' present:',wave_des%NKPTS
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'trying to continue reading vaspwave.h5, but it might fail'
       ENDIF
       IF (ISTART ==1) THEN

          LDIFF=.FALSE.
          DO I=1,3
             DO J=1,3
                IF (ABS(LATT_INI%A(I,J)-LATT_CUR%A(I,J)) > 1E-4) LDIFF=.TRUE.
             ENDDO
          ENDDO
          IF (ENMAXI /= wave_des%ENMAX) LDIFF=.TRUE.
          IF (LDIFF) THEN
             IF (IU0>=0) &
                  WRITE(IU0,*)'WAVECAR: different cutoff or change in lattice found'
          ENDIF

       ENDIF
    else
       latt_ini%a = latt_cur%a
       call lattic(latt_ini)
       if (istart==3) then
          if (iu0 >= 0) then
             write(iu0,*) "ERROR: can't restart (ISTART=3) with wavefunctions on file"
          endif
       else
          istart = 0
       endif
    endif
    ! close group
    VH5_CHECK( vh5_group_close(groupid, wave_des%comm) )
  end function vh5_read_wavefunction_header


  !> @brief read wave functions
  !>
  !> @param fileid hdf5 file handle
  !> @param io io handle
  !> @param wave_des wave function descriptor (see wave.f)
  !> @param symm contains information about used symmetry in 1. BZ
  !> @param wave_spin wave functions (see wave.f)
  !> @param latt_cur lattice definition
  !> @param latt_ini lattice definition
  !> @param istart
  !> @param efermi fermi energy
  function vh5_read_wavefunctions(fileid, io, wave_des, wave_spin, symm, kpoints, grid, latt_cur, latt_ini, istart, efermi, group ) result(ierr)
    USE base
    USE string
    USE wave_high
    use mgrid
    USE lattice
    USE main_mpi
    USE full_kpoints , ONLY:IKPT_SELECT, FIND_KPOINTS_IN_WAVECAR
    USE mkpoints

    integer(HID_T), intent(in)       :: fileid
    type(in_struct) :: io
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(symmetry) :: symm
    type(kpoints_struct) :: kpoints
    type(grid_3d)   :: grid
    type(latt)      :: latt_ini
    type(latt)      :: latt_cur
    integer         :: istart
    real (q)        :: efermi
    character(len=*),optional,intent(in) :: group
    !-------------------------------------------------------------
    ! local variables
    integer(HID_T) :: groupid, spin_groupid, k_groupid
    integer :: ierr, i, j, k, n, ireclw_old, isp, kpl
    integer :: npl, nmaxf, ifail, irecl, is, ispinread, iu0
    integer :: npl_tot, nbandf, nkptsf, malloc, ncomm, nplread
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl, enmax, rbandf, rkptsf, rtag, enmaxf
    logical :: ldiff
    logical, pointer :: lconjg(:)=>null()
    integer, pointer :: index(:)=>null()
    character(LEN=256) :: spingroupname, kgroupname
    type(wavedes1) wave_des1
    ! FIXME CW is replaced by CPTFWP via macro in symbol.inc. still need to find out
    ! why, and what the exact effects are
    complex(q), allocatable :: CW(:)
    REAL(q) VKPT(3)
    INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
    COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
    CHARACTER(:),ALLOCATABLE :: MSG
    LOGICAL ::  SINGLE_PREC, SPAN_RECORDS
    REAL(q), ALLOCATABLE :: rbuf(:)
    COMPLEX(q), ALLOCATABLE :: cbuf(:)
    integer(HSIZE_T) :: start(2), count(2)
    INTEGER :: nb_min
    INTEGER :: NKPT_COUNT
    character(len=MAX_LEN_GROUP) :: my_group
    integer :: node_me
    integer :: ionode, ii(1)
    integer, allocatable :: ikpt_set(:)
    !----------------------------------------------------------------------
    node_me=0
    ionode=0
#ifdef MPI
    node_me=comm%node_me
    ionode =comm%ionode
#endif

    my_group = GRP_WAVE;     if (present(group))    my_group = group

    iu0 = io%iu0
    ! check if group is present:
    ierr= vh5_check_group_presence( fileid, wave_des%comm, iu0, trim( my_group ) )
    if( ierr /= 0 ) then
       if( trim( my_group ) == GRP_WAVE ) then
       else
          if( iu0>=0 ) write(*,*)"You can create one as follows:"
          if( iu0>=0 ) write(*,*)"h5copy -i vaspwave.h5 -o vaspwave.h5 -s /wave -d "//trim(my_group)//" -p"
          if( io%iu6>=0 ) write(io%iu6,*)"Your can create one as follows:"
          if( io%iu6>=0 ) write(io%iu6,*)"h5copy -i vaspwave.h5 -o vaspwave.h5 -s /wave -d "//trim(my_group)//" -p"
       endif
       return
    endif

    if (iu0 >= 0) then
       write(iu0,*) 'found vaspwave.h5, reading group = ', trim(my_group)
    endif

    VH5_CHECK( vh5_group_open(fileid, trim(my_group), groupid, wave_des%comm) )

    RTAG=0
    VH5_CHECK(vh5_read(groupid, "rispin", rispin, wave_des%comm))

    ISPINREAD=NINT(RISPIN)
    VH5_CHECK(vh5_read(groupid, "rnkpts", rkptsf, wave_des%comm))
    VH5_CHECK(vh5_read(groupid, "rnb_tot", rbandf, wave_des%comm))
    VH5_CHECK(vh5_read(groupid, "enmax", enmaxf, wave_des%comm))
    VH5_CHECK(vh5_read(groupid, "amat", latt_ini%a, wave_des%comm))
    VH5_CHECK(vh5_read(groupid, "efermi", efermi, wave_des%comm))
    wave_spin%efermi = efermi

    nkptsf=nint(rkptsf)
    nbandf=nint(rbandf)
    if( allocated(ikpt_select) ) then
       call vtutor%bug("ikpt_select is allocated, this should not be happending!", __FILE__, __LINE__)
    endif
!
! if number of k-points in WAVECAR > then current number of k-points, we need to select
! appropriate k-points
!
    if ( io%ldownsample .and. nkptsf > wave_des%nkpts ) then
        ! this selects the proper k-point from file
        ALLOCATE(ikpt_select( nkptsf ))
        ALLOCATE(ikpt_set( nkptsf ) )
        ikpt_set = -1

        if ( nkptsf > wave_des%nkpts ) then
          if ( iu0>=0 )write(iu0,'(a,i4,a)',advance='no')&
             ' wavecar contains more than',wave_des%nkpts,&
             ' k-points, downsampling...'
        endif

        call find_kpoints_in_wavecar(wave_des,latt_cur,symm,kpoints,nkptsf,ikpt_select,io%iu0)

        !> ikpt_select maps from file -> wave_des
        !> ikpt_set maps from wave_des -> file
        do k = 1, wave_des%nkpts
           ii = findloc( ikpt_select, k )
           ikpt_set( k) = ii(1)
        enddo
    else
        ! legacy way of reading more plane-wave coefficients that stored on file:
        allocate(ikpt_set( wave_des%nkpts ))

        ! read the first k-points and set them too
        do k = 1, wave_des%nkpts
           ikpt_set(k) = k
        enddo
        if ( nkptsf > wave_des%nkpts ) then
           ! file contains more k-points than job
           ! ignore the remaining ones
        else
           ! file contains same or less k-points
           ! continue reading last k-point
           do k = nkptsf, size( ikpt_set )
              ikpt_set(k) = nkptsf
           enddo
        endif
    endif

    !> the number of bands read in might be smaller
    nb_min = min( nbandf, wave_des%nb_tot )
    allocate( cbuf(1:nbandf))
    allocate( rbuf(1:nbandf))

    CALL LATTIC(LATT_INI)
    !=======================================================================
    ! read WAVECAR file, number of bands agree
    !=======================================================================
    IF (.NOT.(wave_des%NRSPINORS==2 .AND. NBANDF*wave_des%NRSPINORS == wave_des%NB_TOT)) THEN

       spin:    DO ISP=1, MIN(wave_des%ISPIN, ISPINREAD)
          write(spingroupname, '(2A)') "spin_"//trim(str(isp))
          VH5_CHECK(vh5_group_open(groupid, trim(spingroupname), spin_groupid, wave_des%comm))
          !NKPT_COUNT= KPOINTS%NKDIM - KPOINTS%NKPTS_NON_ZERO
          NKPT_COUNT= 0
          kpts: DO K = 1, WAVE_DES%NKPTS

             if ( ikpt_set(k) < 1 ) then
                 if(io%iu0>=0 .and. io%nwrite>2) write(iu0,'(a,i5,a)')' kpoint,',k,' skipped.'
                 cycle kpts
             endif
             ! read correct k-point from file
             write(kgroupname, '(2a)') "kpoint_"//trim(str( ikpt_set(k) ))
             if( io%nwrite>2 .and. io%iu0 >=0 ) write(*,*)' mapping k-point: '//str(k)//' -> '//str( ikpt_set(k) )

             NKPT_COUNT = NKPT_COUNT + 1
             VH5_CHECK(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid, wave_des%comm))

             CALL SETWDES(wave_des,wave_des1,K)
#ifdef gammareal
             CALL GEN_INDEX_GAMMA(wave_des1, LATT_CUR%B, wave_des%ENMAX, INDEX, LCONJG)
#endif

             io_begin
             VH5_CHECK(vh5_read(k_groupid, "vkpt", vkpt))
             VH5_CHECK(vh5_read(k_groupid, "celtot", cbuf))
             wave_spin%CELTOT(1:nb_min,K,ISP) = cbuf(1:nb_min)
             VH5_CHECK(vh5_read(k_groupid, "fertot", rbuf))
             wave_spin%FERTOT(1:nb_min,K,ISP) = rbuf(1:nb_min)
             io_end

             VH5_CHECK(vh5_read(k_groupid, "num_planewaves", nplread, wave_des%comm))
             NPL=wave_des%NPLWKP_TOT(K)

             MALLOC=MAX(NPL, NPLREAD)
             ALLOCATE(CW1(MALLOC),CW2(MALLOC),IND(MALLOC),INDI(MALLOC))

             io_begin
             ! create an index to allow for change of cutoff or cell size
             IF (ISTART==2) THEN
                ! gK: fix ISTART=2 in INWAV_FAST
                CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                     wave_des%ENMAX, ENMAXF, NPL/wave_des%NRSPINORS, NPLREAD/wave_des%NRSPINORS, IND, INDI, MALLOC, IFAIL )
             ELSE
                CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                     wave_des%ENMAX, ENMAXF, NPL/wave_des%NRSPINORS, NPLREAD/wave_des%NRSPINORS, IND, INDI, MALLOC, IFAIL )
             ENDIF
             io_end

             CALLMPI( M_bcast_i( wave_des%COMM, IFAIL, 1))
             IF (IFAIL /=0) THEN
                IF ( IU0>=0 )WRITE(IU0,*)' '
                VH5_CHECK(vh5_group_close(k_groupid, wave_des%comm))
                VH5_CHECK(vh5_group_close(spin_groupid, wave_des%comm))
                GOTO 220
             ENDIF

             band: DO J=1,NBANDF
                IF (J>wave_des%NB_TOT) CYCLE
                start = [1,J]
                count = [nplread,1]

                io_begin
                VH5_CHECK(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2))

                CW1=0
                DO IS=1,wave_des%NRSPINORS
                   ! store the wave function coefficients according to new cutoff
                   CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                        CW1((IS-1)*NPL/wave_des%NRSPINORS+1), CW2((IS-1)*NPLREAD/wave_des%NRSPINORS+1))
                ENDDO
                io_end
#ifdef MPI
                IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                   CALLMPI( M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1)) )
                END IF
#endif
#ifdef gammareal
                ! gamma only go from seriel layout to parallel layout if required
                CALL DIS_PW_BAND_GAMMA(wave_des1, J, CW1, wave_spin%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
                CALL DIS_PW_BAND(wave_des1, J, CW1, wave_spin%CW(1,1,K,ISP))
#endif
             ENDDO band

#ifdef gammareal
             CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
             DEALLOCATE(CW1,CW2,IND,INDI)
             VH5_CHECK(vh5_group_close(k_groupid, wave_des%comm))
          ENDDO kpts

          ! copy eigenvalues and weights to all nodes
          NCOMM=wave_des%NB_TOT*wave_des%NKPTS
          CALLMPI( M_bcast_z(wave_des%COMM, wave_spin%CELTOT(1,1,ISP),NCOMM ))
          CALLMPI( M_bcast_d(wave_des%COMM, wave_spin%FERTOT(1,1,ISP),NCOMM ))

          VH5_CHECK(vh5_group_close(spin_groupid, wave_des%comm))


          IF ( NKPT_COUNT == wave_des%NKPTS ) THEN
             IF(IU0>=0)WRITE(IU0,*)'succeeded.'
          ELSE
             call vtutor%alert("K-points found ("//str( NKPT_COUNT)//") inconsistent with requested points ("//&
                str(wave_des%NKPTS)//"). You might want to regenerate the wavefunction and restart this job!")
          ENDIF
       ENDDO spin

       IF (ISPINREAD > wave_des%ISPIN .AND. IU0>=0 ) THEN
          WRITE(IU0,*) 'down-spin wavefunctions not read'
       ENDIF

       IF (NBANDF<wave_des%NB_TOT) THEN
          IF (IU0>=0) WRITE(IU0,*) 'random initialization beyond band ',NBANDF
          CALL WFINIT(wave_des, wave_spin, 1E10_q, NBANDF+1) ! ENINI=1E10 not cutoff restriction
       ENDIF
       IF (IU0 >= 0) WRITE(IU0,*) 'the vaspwave.h5 file was read successfully'

       IF (wave_des%ISPIN<=ISPINREAD) THEN
          GOTO 260 ! close group, free buffers and return
       ENDIF
       !
       !  spin down is missing
       !
       IF (IU0>=0) &
            WRITE(IU0,*) 'No down-spin wavefunctions found', &
            &             ' --> setting down-spin equal up-spin ...'

       DO K=1,wave_des%NKPTS
          wave_spin%CELTOT(1:wave_des%NB_TOT,K,2)=wave_spin%CELTOT(1:wave_des%NB_TOT,K,1)
          wave_spin%FERTOT(1:wave_des%NB_TOT,K,2)=wave_spin%FERTOT(1:wave_des%NB_TOT,K,1)
          NPL=wave_des%NPLWKP(K)
          wave_spin%CW(1:NPL,1:wave_des%NBANDS,K,2)=wave_spin%CW(1:NPL,1:wave_des%NBANDS,K,1)
       ENDDO

       GOTO 260 ! close group, free buffers and return
       !=======================================================================
       ! read collinear WAVECAR file for a non collinear run
       !=======================================================================
    ELSE
       wave_spin%CW=0

       spin2:    DO ISP=1, ISPINREAD
          IF (IU0>=0.AND. ISP==1) &
               WRITE(IU0,*) 'reading wavefunctions of collinear run, up'
          IF (IU0>=0.AND. ISP==2) &
               WRITE(IU0,*) 'reading wavefunctions of collinear run, down'
          write(spingroupname, '(2A)') "spin_"//trim(str(isp))
          VH5_CHECK(vh5_group_open(groupid, trim(spingroupname), spin_groupid, wave_des%comm))
          NKPT_COUNT= KPOINTS%NKDIM - KPOINTS%NKPTS_NON_ZERO
          kpts2: do k = 1,wave_des%nkpts
             if ( ikpt_set(k) < 1 ) then
                 if(io%iu0>=0 .and. io%nwrite>2) write(iu0,'(a,i5,a)')' kpoint,',k,' skipped.'
                 cycle  kpts2
             endif
             write(kgroupname, '(2A)') "kpoint_"//trim(str( ikpt_set(k)))
             if( io%nwrite>2 .and. io%iu0 >=0 ) write(*,*)' mapping k-point: '//str(k)//' -> '//str( ikpt_set(k) )

             NKPT_COUNT = NKPT_COUNT + 1
             VH5_CHECK(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid, wave_des%comm))

             CALL SETWDES(wave_des,wave_des1,K)

             io_begin
             VH5_CHECK(vh5_read_double_array_1d(k_groupid, "vkpt", vkpt))
             VH5_CHECK(vh5_read_double_complex_array_1d(k_groupid, "celtot",cbuf ))
             VH5_CHECK(vh5_read_double_array_1d(k_groupid, "fertot", rbuf))
             do j=1, nb_min
                wave_spin%CELTOT(j+wave_des%nb_tot*(isp-1),K,1) = cbuf(j)
                wave_spin%FERTOT(j+wave_des%nb_tot*(isp-1),K,1) = rbuf(j)
             enddo
             io_end

             VH5_CHECK(vh5_read(k_groupid, "num_planewaves", nplread, wave_des%comm))
             NPL=wave_des%NPLWKP_TOT(K)/2

             MALLOC=MAX(NPL, NPLREAD)
             ALLOCATE(CW1(2*MALLOC),CW2(2*MALLOC),IND(MALLOC),INDI(MALLOC))

             io_begin
             IF (ISTART==2) THEN
                ! gK: fix ISTART=2 in INWAV_FAST
                CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                     wave_des%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
             ELSE
                CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                     wave_des%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
             ENDIF ! istart == 2
             io_end

             CALLMPI( M_bcast_i( wave_des%COMM, IFAIL, 1))
             IF (IFAIL /=0) THEN
                IF ( IU0>=0 )WRITE(IU0,*)' '
                VH5_CHECK(vh5_group_close(k_groupid, wave_des%comm))
                VH5_CHECK(vh5_group_close(spin_groupid, wave_des%comm))
                GOTO 230
             ENDIF

             IF (ISP==1) THEN
                DO J=1,NBANDF
                   IF (J>wave_des%NB_TOT) CYCLE
                   start = [1,J]
                   count = [nplread,1]

                   io_begin
                   VH5_CHECK(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2))

                   CW1=0
                   CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1, CW2)
                   io_end

#ifdef MPI
                   IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                      CALLMPI( M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1)) )
                   END IF ! wave_des
#endif
                   CALL DIS_PW_BAND(wave_des1, J, CW1, wave_spin%CW(1,1,K,ISP))

                   ! copy immediately to second panel
                   wave_spin%CELTOT(J+nb_min,K,1)=wave_spin%CELTOT(J,K,1)
                   wave_spin%FERTOT(J+nb_min,K,1)=wave_spin%FERTOT(J,K,1)

                   CW1(NPL+1:2*NPL)=CW1(1:NPL)
                   CW1(1:NPL)=0

                   CALL DIS_PW_BAND(wave_des1, J+nb_min, CW1, wave_spin%CW(1,1,K,1))
                ENDDO
             ELSE ! ISP==1
                DO J=1,NBANDF
                   IF (J>wave_des%NB_TOT) CYCLE
                   start = [1,J]
                   count = [nplread,1]

                   io_begin
                   VH5_CHECK(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2))

                   CW1=0
                   CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1(NPL+1), CW2)
                   io_end
#ifdef MPI
                   IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                      CALLMPI( M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1)) )
                   END IF
#endif
                   CALL DIS_PW_BAND(wave_des1, J+nb_min, CW1, wave_spin%CW(1,1,K,1))
                ENDDO

             ENDIF  ! ISP == 1

             DEALLOCATE(CW1,CW2,IND,INDI)

             VH5_CHECK(vh5_group_close(k_groupid, wave_des%comm))
          ENDDO kpts2

          IF ( NKPT_COUNT == wave_des%NKPTS .AND. ISP == MIN(wave_des%ISPIN, ISPINREAD) ) THEN
             IF(IU0>=0)WRITE(IU0,*)'succeeded'
          ELSE
             IF( NKPT_COUNT < wave_des%NKPTS ) THEN
                IF ( IU0>=0 )WRITE(IU0,*)' FAILED! '
                VH5_CHECK(vh5_group_close(spin_groupid, wave_des%comm))
                GOTO 250
             ENDIF
          ENDIF

          VH5_CHECK(vh5_group_close(spin_groupid, wave_des%comm))
       ENDDO spin2

       ! copy eigenvalues and weights to all nodes
       NCOMM=wave_des%NB_TOT*wave_des%NKPTS
       CALLMPI( M_bcast_z(wave_des%COMM, wave_spin%CELTOT(1,1,1),NCOMM ))
       CALLMPI( M_bcast_d(wave_des%COMM, wave_spin%FERTOT(1,1,1),NCOMM ))

       IF (IU0 >= 0) WRITE(IU0,*) 'the vaspwave.h5 file was read successfully'
       GOTO 260 ! close group, free buffers and return
    ENDIF ! num bands agree

    MSG = 'ERROR: while reading vaspwave.h5, file is incompatible'
    IF (wave_des%ENMAX /= ENMAXF) &
       MSG = MSG // ' the energy cutoff has changed (new,old) ' // &
          str(wave_des%ENMAX) // " " // str(ENMAXF)
    IF (NBANDF /= wave_des%NB_TOT) &
       MSG = MSG // ' the number of bands has changed (new,old) ' // &
          str(wave_des%NB_TOT) // " " // str(NBANDF)
    IF (NKPTSF /= wave_des%NKPTS) &
       MSG = MSG // ' the number of k-points has changed (new,old) ' // &
          str(wave_des%NKPTS) // " " // str(NKPTSF)

    CALL vtutor%error(MSG)

    ! can not do anything with WAVECAR
    ! hard stop, pull all breaks
    !=======================================================================
220 CONTINUE

#if defined( wNGZhalf ) || defined( wNGXhalf) || defined( NGZhalf ) || defined( NGXhalf)
    CALL vtutor%error('ERROR: while reading vaspwave.h5, plane wave coefficients changed ' &
       // str(NPL) // ' ' // str(NPLREAD))
#else
      CALL vtutor%error("ERROR: while reading vaspwave.h5, plane wave coefficients changed.&
              & If you intend to read the wavefunction from a vasp_std job,&
              & double NBANDS in INCAR to continue.")
#endif

230 CONTINUE

#if defined( wNGZhalf ) || defined( wNGXhalf) || defined( NGZhalf ) || defined( NGXhalf)
#else
      CALL vtutor%error("ERROR2: while reading vaspwave.h5, plane wave coefficients changed.&
              & If you intend to read the wavefunction from a vasp_std job,&
              & double NBANDS in INCAR to continue.")
#endif

250 CONTINUE

#if defined( wNGZhalf ) || defined( wNGXhalf) || defined( NGZhalf ) || defined( NGXhalf)
#else
      CALL vtutor%error("ERROR, while reading vaspwave.h5 of vasp_std: " // str(NKPT_COUNT) // " out of " // &
         str(wave_des%NKPTS)//" k-points found. Using a wavefunction from a vasp_std job with ISYM=0(-1) might be successful here.")
#endif

260 CONTINUE

    ! close group
    VH5_CHECK( vh5_group_close(groupid, wave_des%comm) )
    deallocate( cbuf, rbuf)
    !> get rid of IKPT_SELECT if number of k-points same on file
    IF ( NKPTSF==wave_des%NKPTS.and.ALLOCATED(IKPT_SELECT)) DEALLOCATE(IKPT_SELECT)
    RETURN

  end function vh5_read_wavefunctions

  !> @brief write projectors to hdf5 file
  subroutine vh5_write_projectors(fileid, typeinfo, wave_des, wave_spin, lorbit, &
       lpar, par, lchar, phase, group, subgroup)
    integer(HID_T), intent(in)   :: fileid !< hdf5 file handle
    type(type_info), intent(in)  :: typeinfo
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    integer :: lorbit, lpar
    real(q), dimension(wave_des%nb_tot, wave_des%nkpts, lpar, typeinfo%nionp, wave_des%ncdij) :: par
    character(len=*), dimension(:) :: lchar
    complex(q), dimension(:,:,:,:,:) :: phase
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PROJECTORS)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=256) :: lcharout
    integer :: nl
    integer(HSIZE_T), dimension(5) :: par_shape;
    integer(HSIZE_T), dimension(5) :: phase_shape;
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group    = GRP_RESULTS;    if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PROJECTORS; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    ! write procar type
    if (lorbit == 11) then
       VH5_CHECK(vh5_write(subgroupid, "procar_type", 'PROCAR lm decomposed'))
    else if (lorbit == 12) then
       VH5_CHECK(vh5_write(subgroupid, "procar_type", 'PROCAR lm decomposed + phase'))
    else
       VH5_CHECK(vh5_write(subgroupid, "procar_type", 'PROCAR new format'))
    end if
    !
    ! write parameters
    ! FIXME these parameters are also written elsewhere. we need to discuss how to structure
    ! a file where these things are only written once. ATM, I'm writing them out because
    ! I'm not sure whether they correspond to the same state of the program as the other
    ! instances
    VH5_CHECK(vh5_write(subgroupid, "nb_tot", wave_des%nb_tot))
    VH5_CHECK(vh5_write(subgroupid, "nionp", typeinfo%nionp))
    VH5_CHECK(vh5_write(subgroupid, "nkpoints", wave_des%nkpts))
    VH5_CHECK(vh5_write(subgroupid, "kpoints", wave_des%vkpt))
    VH5_CHECK(vh5_write(subgroupid, "celtot", wave_spin%celtot))
    !
    ! write lchars
    !
    VH5_CHECK(vh5_write(subgroupid, "lchar", lchar))
    !
    ! write par
    !
    par_shape = shape(par)
    VH5_CHECK(vh5_write_double_array_nd(subgroupid, "par", 5, par_shape, par))
    !
    ! write phase
    !
    if (lorbit == 12) then
       phase_shape = shape(phase)
       VH5_CHECK(vh5_write_double_complex_array_nd(subgroupid, "phase", 5, phase_shape, phase))
    end if
    !
    ! calculate and write parsum and sumion
    !

    !
    ! close/deallocate
    !
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_projectors

  subroutine vh5_write_locproj(fileid, lprj_functions, lprj_covl, lnoncol, group, subgroup)
    use string, only: str

    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    TYPE (LPRJ_function), ALLOCATABLE :: lprj_functions(:)
    GDEF, ALLOCATABLE :: LPRJ_COVL(:,:,:,:)
    integer(HID_T) :: nproj !nr. of proj.
    integer(HID_T):: ip, groupid, subgroupid, paramgroupid
    logical, intent(in), optional :: lnoncol
    character(len=*), intent(in), optional :: group
    character(len=*), intent(in), optional :: subgroup
    character(len=15), allocatable :: radial_type_list(:)
    character(len=10), allocatable :: ang_type_list(:)
    REAL(q), allocatable :: lprj_coord_list(:,:)           !< position of the localized orbital in reduced coordinates
    CHARACTER(LEN=15) :: PRKIND(3)=(/"PAW projector  ","PS partial wave","Hydrogen-like  "/)
    CHARACTER(LEN=10) :: LABEL(1:7,-5:3)
    DATA LABEL(1:7,-5) &
     &/"  sp3d2-1 ","  sp3d2-2 ","  sp3d2-3 ","  sp3d2-4 ","  sp3d2-5 ","  sp3d2-6 ","          "/
      DATA LABEL(1:7,-4) &
     &/"   sp3d-1 ","   sp3d-2 ","   sp3d-3 ","   sp3d-4 ","   sp3d-5 ","          ","          "/
      DATA LABEL(1:7,-3) &
     &/"   sp3-1  ","   sp3-2  ","   sp3-3  ","   sp3-4  ","          ","          ","          "/
      DATA LABEL(1:7,-2) &
     &/"   sp2-1  ","   sp2-2  ","   sp2-3  ","          ","          ","          ","          "/
      DATA LABEL(1:7,-1) &
     &/"    sp-1  ","    sp-2  ","          ","          ","          ","          ","          "/
      DATA LABEL(1:7, 0) &
     &/"      s   ","          ","          ","          ","          ","          ","          "/
      DATA LABEL(1:7, 1) &
     &/"     py   ","     pz   ","     px   ","          ","          ","          ","          "/
      DATA LABEL(1:7, 2) &
     &/"    dxy   ","    dyz   ","    dz2   ","    dxz   ","   dx2-y2 ","          ","          "/
      DATA LABEL(1:7, 3) &
     &/"fy(3x2-y2)","    fxyz  ","    fyz2  ","    fz3   ","    fxz2  "," fz(x2-y2)","fx(x2-3y2)"/
    character(LEN=40) :: my_group, my_subgroup

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_LOCPROJ; if (present(subgroup)) my_subgroup = subgroup

    ! create / open groups
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))
    VH5_CHECK(vh5_group_open_or_create(subgroupid, trim("parameters"), paramgroupid))

    ! allocate lists of strings
    nproj=SIZE(LPRJ_COVL,4)
    allocate(radial_type_list(nproj))
    allocate(ang_type_list(nproj))
    allocate(lprj_coord_list(3,nproj))
    ! write string info in lists first
    do ip=1,nproj
      radial_type_list(ip) = trim(adjustl(prkind(lprj_functions(ip)%radial_type)))
      ang_type_list(ip) = trim(adjustl(LABEL(lprj_functions(ip)%m,lprj_functions(ip)%l)))
      lprj_coord_list(:,ip) = lprj_functions(ip)%r
    enddo

    ! write locproj parameters
    VH5_CHECK(vh5_write(paramgroupid, "site", lprj_functions(:)%poscar_site))
    VH5_CHECK(vh5_write(paramgroupid, "coordinates", lprj_coord_list))
    VH5_CHECK(vh5_write(paramgroupid, "rad_type", radial_type_list))
    VH5_CHECK(vh5_write(paramgroupid, "ang_type", ang_type_list))
    VH5_CHECK(vh5_write(paramgroupid, "orbital", lprj_functions(:)%n))
    VH5_CHECK(vh5_write(paramgroupid, "alpha", lprj_functions(:)%za))
    VH5_CHECK(vh5_write(paramgroupid, "species_number", lprj_functions(:)%species))
    if (present(lnoncol)) then
      VH5_CHECK(vh5_write(paramgroupid, "lnoncollinear", lnoncol))
    endif

    ! the locproj data itself and a small format header
    VH5_CHECK(vh5_write(subgroupid, "data", GCONJG(lprj_covl)))
    VH5_CHECK(vh5_write(subgroupid, "format", "[proj_index, spin_index, kpts_index, band_index]"))

    VH5_CHECK(vh5_group_close_writing(paramgroupid))
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine

  !> write the k-points along path to a file
  subroutine kpoints_along_path_hdf5write(fileid,kpoints_ikibz,kpoints_dists,kpoints_coord,group,subgroup)
     integer(HID_T), intent(in) :: fileid !< hdf5 file handle
     real(q) :: kpoints_dists(:)
     integer :: kpoints_ikibz(:)
     real(q) :: kpoints_coord(:,:)
     character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
     character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PROJECTORS)
     !----------------------------------------------------------------------
     integer(HID_T) :: groupid, subgroupid
     character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
     !----------------------------------------------------------------------
     my_group    = GRP_RESULTS;         if (present(group))    my_group    = trim(group)
     my_subgroup = SUBGRP_KPOINTS_LINE; if (present(subgroup)) my_subgroup = trim(subgroup)

     VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
     VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

     VH5_CHECK(vh5_write(subgroupid, "kpoints_ikibz", kpoints_ikibz))
     VH5_CHECK(vh5_write(subgroupid, "kpoints_dists", kpoints_dists))
     VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", kpoints_coord))
  end subroutine kpoints_along_path_hdf5write

  !> @brief writes header info of four-point coulomb integrals and related quanties
  subroutine vh5_write_xijkl_header(fileid, nwtot, ispin, comega, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    integer, intent(in) :: nwtot !< number of basis functions
    integer, intent(in) :: ispin !< number of spin indices
    complex(q), intent(in) :: comega(:) !< array of complex frequency points considered
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, 'nwtot', nwtot))
    VH5_CHECK(vh5_write(subgroupid, 'ispin', ispin))
    VH5_CHECK(vh5_write(subgroupid, 'nomega', size( comega )))
    VH5_CHECK(vh5_write(subgroupid, 'comega', comega))
    if (ispin == 1) then
       VH5_CHECK(vh5_write(subgroupid, 'spin_labels', ['total']))
    else
       VH5_CHECK(vh5_write(subgroupid, 'spin_labels', ['up~up    ', 'down~down', 'up~down  ']))
    end if
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_header

  !> @brief writes wannier position in two-centre coulomb integrals
  subroutine vh5_write_xijkl_2d(fileid, tag, xijkl, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    character(len=*), intent(in) :: tag !< tag of matrix elements, vijkl, wijkl, uijkl or pijkl
    complex(q), intent(in) :: xijkl(:,:) !< positions of wannier in two-centre coulomb integrals
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, tag, xijkl))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_2d

  !> @brief writes wannier position in two-centre coulomb integrals
  subroutine vh5_write_xijkl_2d_real(fileid, tag, xijkl, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    character(len=*), intent(in) :: tag !< tag of matrix elements, vijkl, wijkl, uijkl or pijkl
    real(q), intent(in) :: xijkl(:,:) !< positions of wannier in two-centre coulomb integrals
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, tag, xijkl))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_2d_real

  !> @brief writes bare coulomb integrals
  subroutine vh5_write_xijkl_4d(fileid, tag, xijkl, group, subgroup, compact_dim)
    use tutor, only: vtutor
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    character(len=*), intent(in) :: tag !< tag of matrix elements, vijkl, wijkl, uijkl or pijkl
    GDEF, intent(in) :: xijkl(:,:,:,:) !< four-point bare coulomb integrals
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    logical, intent(in), optional :: compact_dim !< Should we compress the first 4 dimensions into 1 dimension
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    GDEF, allocatable :: work(:), work3d(:,:,:)
    logical my_compact_dim
    integer n1, n2, n3, n4, nglobal, nspin, index_
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    my_compact_dim = .false.;  if (present(compact_dim)) my_compact_dim = compact_dim
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    if (my_compact_dim) then
       nglobal = size(xijkl, 1)
       if ((size(xijkl, 2) /= nglobal).or.(size(xijkl, 3) /= nglobal).or.(size(xijkl, 4) /= nglobal)) &
          call vtutor%bug("xijkl should have the same dimension for all indices", __FILE__, __LINE__)
       allocate(work(size(xijkl)))
       do n1 = 1, nglobal
          do n2 = 1, nglobal
             do n3 = 1, nglobal
                do n4 = 1, nglobal
                   index_ = ((((n1 - 1) * nglobal) + (n2 - 1) * nglobal) + (n3 - 1)) * nglobal + n4
                   work(index_) = xijkl(n1,n4,n3,n2)
                end do
             end do
          end do
       end do
       VH5_CHECK(vh5_write(subgroupid, tag, work))
    else
       nspin = size(xijkl, 3)
       if (size(xijkl, 4) /= nspin) &
          call vtutor%bug("last two dimensions of xijkl should be spin dimensions", __FILE__, __LINE__)
       if (nspin == 1) then
          VH5_CHECK(vh5_write(subgroupid, tag, xijkl(:,:,:,1)))
       else
          ! we only compute one of the off-diagonal spin components due to spin symmetry
          allocate(work3d(size(xijkl,1), size(xijkl,2), 3))
          work3d(:,:,1) = xijkl(:,:,1,1)
          work3d(:,:,2) = xijkl(:,:,2,2)
          work3d(:,:,3) = xijkl(:,:,1,2)
          VH5_CHECK(vh5_write(subgroupid, tag, work3d))
       end if
    end if

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_4d

  !> @brief writes bare coulomb integrals
  subroutine vh5_write_xijkl_5d(fileid, tag, xijkl, group, subgroup)
    use tutor, only: vtutor
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    character(len=*), intent(in) :: tag !< tag of matrix elements, vijkl, wijkl, uijkl or pijkl
    GDEF, intent(in) :: xijkl(:,:,:,:,:) !< four-point bare coulomb integrals with spin up-up, up-dn, dn-dn
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    GDEF, allocatable :: work(:,:,:,:)
    integer nspin
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    nspin = size(xijkl, 3)
    if (size(xijkl, 4) /= nspin) &
       call vtutor%bug("last two dimensions of xijkl should be spin dimensions", __FILE__, __LINE__)
    if (nspin == 1) then
       VH5_CHECK(vh5_write(subgroupid, tag, xijkl(:,:,:,1,:)))
    else
       ! we only compute one of the off-diagonal spin components due to spin symmetry
       allocate(work(size(xijkl,1), size(xijkl,2), 3, size(xijkl,5)))
       work(:,:,1,:) = xijkl(:,:,1,1,:)
       work(:,:,2,:) = xijkl(:,:,2,2,:)
       work(:,:,3,:) = xijkl(:,:,1,2,:)
       VH5_CHECK(vh5_write(subgroupid, tag, work))
    end if

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_5d

  !> @brief writes frequency dependent coulomb integrals
  subroutine vh5_write_xijkl_6d(fileid, tag, xijkl, group, subgroup)
    use tutor, only: vtutor
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    character(len=*), intent(in) :: tag !< tag of matrix elements, vijkl, wijkl, uijkl or pijkl
    GDEF, intent(in) :: xijkl(:,:,:,:,:,:) !< four-point coulomb integrals with spin up-up, up-dn, dn-dn and frequency
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    GDEF, allocatable :: work(:,:)
    integer n1, n2, n3, n4, nglobal, nspin, index_
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    nglobal = size(xijkl, 1)
    nspin = size(xijkl, 5)
    if ((size(xijkl, 2) /= nglobal).or.(size(xijkl, 3) /= nglobal).or.(size(xijkl, 4) /= nglobal)) &
       call vtutor%bug("xijkl should have the same dimension for all indices", __FILE__, __LINE__)
    if (size(xijkl, 6) /= nspin) &
       call vtutor%bug("last two dimensions of xijkl should be spin dimensions", __FILE__, __LINE__)
    allocate(work(nglobal**4, 2*nspin-1))
    do n1 = 1, nglobal
       do n2 = 1, nglobal
          do n3 = 1, nglobal
             do n4 = 1, nglobal
                index_ = (((n1 - 1) * nglobal + n2 - 1) * nglobal + n3 - 1) * nglobal + n4
                work(index_,1) = xijkl(n1,n4,n3,n2,1,1)
                ! we only compute one of the off-diagonal spin components due to spin symmetry
                if (nspin == 2) then
                   work(index_,2) = xijkl(n1,n4,n3,n2,2,2)
                   work(index_,3) = xijkl(n1,n4,n3,n2,1,2)
                end if
             end do
          end do
       end do
    end do
    VH5_CHECK(vh5_write(subgroupid, tag, work))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_6d

  !> @brief writes frequency dependent coulomb integrals
  subroutine vh5_write_xijkl_7d(fileid, tag, xijkl, group, subgroup)
    use tutor, only: vtutor
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    character(len=*), intent(in) :: tag !< tag of matrix elements, vijkl, wijkl, uijkl or pijkl
    GDEF, intent(in) :: xijkl(:,:,:,:,:,:,:) !< four-point coulomb integrals with spin up-up, up-dn, dn-dn, frequency and points
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    GDEF, allocatable :: work(:,:,:)
    integer n1, n2, n3, n4, nglobal, nspin, index_
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = SUBGRP_CRPA; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    nglobal = size(xijkl, 1)
    nspin = size(xijkl, 5)
    if ((size(xijkl, 2) /= nglobal).or.(size(xijkl, 3) /= nglobal).or.(size(xijkl, 4) /= nglobal)) &
       call vtutor%bug("xijkl should have the same dimension for all indices", __FILE__, __LINE__)
    if (size(xijkl, 6) /= nspin) &
       call vtutor%bug("two spin dimensions of xijkl should be the same size", __FILE__, __LINE__)
    allocate(work(nglobal**4, 2*nspin-1, size(xijkl, 7)))
    do n1 = 1, nglobal
       do n2 = 1, nglobal
          do n3 = 1, nglobal
             do n4 = 1, nglobal
                index_ = (((n1 - 1) * nglobal + n2 - 1) * nglobal + n3 - 1) * nglobal + n4
                work(index_,1,:) = xijkl(n1,n4,n3,n2,1,1,:)
                ! we only compute one of the off-diagonal spin components due to spin symmetry
                if (nspin == 2) then
                   work(index_,2,:) = xijkl(n1,n4,n3,n2,2,2,:)
                   work(index_,3,:) = xijkl(n1,n4,n3,n2,1,2,:)
                end if
             end do
          end do
       end do
    end do
    VH5_CHECK(vh5_write(subgroupid, tag, work))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_xijkl_7d

  !> @brief writes screended potential to file
  subroutine vh5_write_wpot(fileid, nq, np, lrealstore, wpot, w, nomega, lfull, ldiag, group, subgroup)
    use chi_glb, only: responsefunction
    use string, only: str
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    integer, intent(in) :: nq  !< q-point
    integer, intent(in) :: np  !< number of plane waves
    logical, intent(in) :: lrealstore !< are we storing real arrays?
    complex(q), intent(in) :: w(:,:) !< screened potential
    type(responsefunction), intent(in) :: wpot !< response function type
    logical, intent(in) :: lfull !< store full W (screened Coulomb) matrix
    logical, intent(in) :: ldiag !< store diag of W (screened Coulomb) matrix
    integer, intent(in), optional :: nomega !< frequency point to be stored
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: SUBGRP_CRPA)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    integer :: nomega_store
    integer :: i
    complex(q) :: head(3,3)
    complex(q), allocatable :: wdumpdiag(:)
    complex(q), allocatable :: wdump(:,:)
    complex(q), allocatable :: wing(:,:)
    real(q), allocatable :: wdumpdiagr(:)
    real(q), allocatable :: wdumpr(:,:)
    real(q), allocatable :: wingr(:,:)
    !-------------------------------------------------------------------------------------

    my_group = GRP_WPOT;     if (present(group))    my_group = group
    my_subgroup = 'qpoint_'//str(nq); if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    nomega_store = 1
    if( present( nomega ) ) nomega_store = nomega

    VH5_CHECK(vh5_write(subgroupid, 'np', np))
    head =  wpot%head(:,:,nomega_store)
    VH5_CHECK(vh5_write(subgroupid, 'head', head))
    VH5_CHECK(vh5_write(subgroupid, 'lrealstore', lrealstore))
    VH5_CHECK(vh5_write(subgroupid, 'lfull', lfull))
    VH5_CHECK(vh5_write(subgroupid, 'ldiag', ldiag))
    VH5_CHECK(vh5_write(subgroupid, 'vqpt', wpot%vkpt))
    VH5_CHECK(vh5_write(subgroupid, 'omega', wpot%comega))
    if( lrealstore ) then
       VH5_CHECK(vh5_write(subgroupid, 'wfullr', wdumpr))
       allocate( wingr( np,3) )
       wingr(1:np, 1:3 ) = wpot%wingr(1:np,:,nomega_store)
       VH5_CHECK(vh5_write(subgroupid, 'wingr', wingr))
       wingr(1:np, 1:3 ) = wpot%cwingr(1:np,:,nomega_store)
       VH5_CHECK(vh5_write(subgroupid, 'cwingr', wingr))
       deallocate( wingr )
       if ( lfull ) then
           allocate( wdumpr(np,np) )
           wdumpr(:,:) = real(w(1:np,1:np),q)
           VH5_CHECK(vh5_write(subgroupid, 'wfullr', wdumpr))
           deallocate( wdumpr )
       endif
       if ( ldiag ) then
           allocate( wdumpdiagr(np) )
           do i =1 , np
              wdumpdiagr(i) = real(w(i,i),q)
           enddo
           VH5_CHECK(vh5_write(subgroupid, 'wr', wdumpdiagr))
           deallocate( wdumpdiagr )
       endif
    else
       allocate( wing( np,3) )
       wing(1:np, 1:3 ) = wpot%wing(1:np,:,nomega_store)
       VH5_CHECK(vh5_write(subgroupid, 'wing', wing) )
       wing(1:np, 1:3 ) = wpot%cwing(1:np,:,nomega_store)
       VH5_CHECK(vh5_write(subgroupid, 'cwing', wing) )
       deallocate( wing )
       if ( lfull ) then
           allocate( wdump(np,np) )
           wdump(:,:) = w(1:np,1:np)
           VH5_CHECK(vh5_write(subgroupid, 'wfull', wdump))
           deallocate( wdump )
       endif
       if ( ldiag ) then
           allocate( wdumpdiag(np) )
           do i =1 , np
              wdumpdiag(i) = w(i,i)
           enddo
           VH5_CHECK(vh5_write(subgroupid, 'w', wdumpdiag))
           deallocate( wdumpdiag )
       endif
    endif

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    !VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_wpot

  !> Read the q-point header from vaspwpot.h5 or vaspwpotfull.h5
  function vh5_read_wpot_header(fileid, nq, vqpt, group) result(ierr)
    use string, only: str
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    integer, intent(in) :: nq  !< q-point
    real(q), intent(out) :: vqpt(:) !< q-points read from file
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    integer :: ierr
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    my_group = GRP_WPOT;     if (present(group))    my_group = group
    my_subgroup = 'qpoint_'//str(nq)

    ierr = vh5_group_open(fileid, trim(my_group), groupid)
    if (ierr/=0) return

    ! check existence of dataset 
    ! only try to open dataset if it is present in file, otherwise drop error
    ! this prevents HDF5 warnings that polute the stdout
    if (vh5_exists(groupid, trim(my_subgroup) )) then
       ierr = vh5_group_open(groupid, trim(my_subgroup), subgroupid)
    else
       vqpt = 10 
       ierr = -1
    endif
    if (ierr/=0) return

    ierr = vh5_read(subgroupid, "vqpt", vqpt)
    if (ierr/=0) return

    VH5_CHECK(vh5_group_close(subgroupid))
    VH5_CHECK(vh5_group_close(groupid))
  end function vh5_read_wpot_header

  !> Read the screened coulomb potential from vaspwpot.h5
  function vh5_read_wpot(fileid, wpot, wgwq, nq, io, lreadfull, group) result(ierr)
    use chi_glb, only: responsefunction
    use string, only: str
    use wave
    use base
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    integer, intent(in) :: nq  !< q-point
    type(wavedes1), intent(in) :: wgwq
    type(responsefunction), intent(inout) :: wpot !< response function type
    type(in_struct), intent(in) :: io
    logical, intent(in) :: lreadfull !< read full W (screened Coulomb) matrix or diagonal only
    integer :: ierr
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    ! read variables
    logical :: lfull, lrealstore
    integer :: np, i, np_wg
    integer :: nomega_store = 1
    character(len=:),allocatable :: filename
    real(q), allocatable :: wr(:),wingr(:,:),wfullr(:,:)
    complex(q), allocatable :: w(:),wing(:,:),wfull(:,:)
    complex(q) :: head(3,3)

    !-------------------------------------------------------------------------------------
    my_group = GRP_WPOT;     if (present(group))    my_group = group
    my_subgroup = 'qpoint_'//str(nq)
    filename = 'unknown file'
    if (fileid==ih5wpotfileid) filename = VASPWPOT
    if (fileid==ih5wpotfullfileid) filename = VASPWPOTFULL

    ierr = vh5_group_open(fileid, trim(my_group), groupid)
    if (ierr/=0) return
    VH5_CHECK( vh5_group_open(groupid, trim(my_subgroup), subgroupid) )
    if (ierr/=0) return

    if (io%iu0 >= 0 .and. io%nwrite > 2) write(io%iu0,*) GRP_WPOT//' found in '//trim(filename)//', reading '//trim(my_subgroup)
    ! legacy routine READ_WPOT checks for np1 and np2 though np1=np2 is always written to WFULL
    VH5_CHECK(vh5_read(subgroupid, "np", np))

    np_wg=wgwq%ngvector
    if (wgwq%lgamma) np_wg=np_wg*2
    if (np_wg /= np) then
          call vtutor%error("number of g-vectors in "//trim(my_group)//"/"//trim(my_subgroup)//" of "//trim(filename)//" (" &
              // str(np) //") does not match " &
              // "the basis set of the response function ("// str(np_wg) // "). \n" &
              // "make sure that ENCUTGW is the same in both calculations.")
    endif

    VH5_CHECK(vh5_read(subgroupid, 'lfull', lfull))
    if (lfull .neqv. lreadfull) then
       ierr = 1
       goto 100
    endif

    VH5_CHECK(vh5_read(subgroupid, 'lrealstore', lrealstore))
    ! sanity check
    if (lrealstore .neqv. wpot%lrealstore) then
       call vtutor%error("Inconsistent lrealstore="//str(lrealstore)//" flag in wpot "//&
                         "group of "//filename//"and response function lrealstore="//str(wpot%lrealstore)//".")
    endif

    ! head is always complex
    VH5_CHECK(vh5_read(subgroupid, 'head', head) )
    wpot%head(:,:,nomega_store) = head

    ! now decide if gamma point is read
    if( wpot%lrealstore)  then
       wpot%responser(:,:,nomega_store)=0
       wpot%wingr(:,:,nomega_store)=0
       wpot%cwingr(:,:,nomega_store)=0

       allocate( wingr(np,3) )
       wingr=0
       VH5_CHECK( vh5_read(subgroupid, 'wingr', wingr) )
       wpot%wingr(1:np,:,nomega_store) = wingr

       wingr=0
       VH5_CHECK( vh5_read(subgroupid, 'cwingr', wingr) )
       wpot%cwingr(1:np,:,nomega_store) = wingr
       deallocate( wingr )

       if( wpot%llead ) then
          if ( lfull ) then
             allocate( wfullr(np,np) )
             VH5_CHECK( vh5_read(subgroupid, 'wfullr', wfullr) )
             do i = 1, np
                wpot%responser(:np,i,nomega_store) = wfullr(:,i)
             enddo
             deallocate(wfullr)
          else
             allocate( wr(np) )
             VH5_CHECK( vh5_read(subgroupid, 'wr', wr) )
             do i = 1, np
                wpot%responser(i,i,nomega_store) = wr(i)
             enddo
             deallocate(wr)
          endif
       endif
    else
       wpot%responsefun(:,:,nomega_store)=0
       wpot%wing(:,:,nomega_store)=0
       wpot%cwing(:,:,nomega_store)=0

       allocate( wing(np,3) )
       wing=0
       VH5_CHECK( vh5_read(subgroupid, 'wing', wing) )
       wpot%wing(1:np,:,nomega_store) = wing

       wing=0
       VH5_CHECK( vh5_read(subgroupid, 'cwing', wing) )
       wpot%cwing(1:np,:,nomega_store) = wing
       deallocate( wing )

       if( wpot%llead ) then
          if ( lfull ) then
             allocate( wfull(np,np) )
             VH5_CHECK( vh5_read(subgroupid, 'wfull', wfull) )
             do i = 1, np
                wpot%responsefun(:np,i,nomega_store) = wfull(:,i)
             enddo
             deallocate(wfull)
          else
             allocate( w(np) )
             VH5_CHECK( vh5_read(subgroupid, 'w', w) )
             do i = 1, np
                wpot%responsefun(i,i,nomega_store) = w(i)
             enddo
             deallocate(w)
          endif
       endif
    endif

    if( trim( my_group ) == GRP_WPOT ) wpot_found = .false.

100 continue

#ifdef use_shmem
    ! force the other ranks  to wait for root core
    if (wpot%lshmem)  then
       CALLMPI( M_barrier(wpot%comm_shmem ) )
    endif
#endif

    VH5_CHECK(vh5_group_close(subgroupid))
    VH5_CHECK(vh5_group_close(groupid))
  end function vh5_read_wpot

  !> @brief writes wanproj to h5
  subroutine vh5_write_wanproj(fileid, nb_tot, ispin, nkpts, nw, vkpts, exclude_bands, lwindow, comm, u, io)
    use mpimy
    use base
    use string
    use tutor, ONLY: vtutor
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    integer, intent(in) :: nb_tot !< number of bands
    integer, intent(in) :: ispin !< number of spin indices
    integer, intent(in) :: nkpts !< number of k-points
    integer, intent(in) :: nw !< number of wannier functions
    real(q), intent(in) :: vkpts(:,:) !< k-point list
    logical, intent(in) :: exclude_bands(:) !< info about excluded bands in projection
    logical, intent(in) :: lwindow(:,:,:) !< energy window
    type(communic), intent(in) :: comm
    GDEF, allocatable, intent(in) :: u(:,:,:,:) !< rotation matrix for each kpoint and spin
    type(in_struct), intent(in) :: io
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    !-------------------------------------------------------------------------------------
    ! local variables
    integer i,j, ip, ipp
    integer isp,nk
    integer ispin_old
    integer nkpts_old
    integer nb_tot_old
    integer nw_old
    integer irec
    integer, allocatable :: included_bands(:)
    integer(HSIZE_T) :: dim_incl_bands(1), count(2)
    character(LEN=256) :: spingroupname, kgroupname
    integer(HID_T) :: spin_groupid, k_groupid
    GDEF, allocatable :: uij(:,:)
    integer :: ierr
    !----------------------------------------------------------------------

    ierr = vh5_check_group_presence(fileid, comm, io%iu0, GRP_WANPROJ )
    if( ierr==0 ) then
       if( io%iu0>=0 ) write(*,*)"dataset "//trim(GRP_WANPROJ)//" found in vaspwave.h5, replacing it ..."
       call vtutor%alert( "vaspwave.h5 might store redundant data. To reduce &
          &its filesize a call to h5repack could be appropriate!" )
       ! delete the entire dataset
       VH5_CHECK(vh5_delete(fileid, GRP_WANPROJ))
    endif

    VH5_CHECK(vh5_group_open_or_create(fileid, GRP_WANPROJ, groupid))

    VH5_CHECK(vh5_write(groupid, 'ispin', ispin))
    VH5_CHECK(vh5_write(groupid, 'nkpts', nkpts))
    VH5_CHECK(vh5_write(groupid, 'nb_tot', nb_tot))
    VH5_CHECK(vh5_write(groupid, 'nw', nw))
    count=[ 3, nkpts ]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "vkpts", 2, count))
    VH5_CHECK(vh5_write(groupid, 'vkpts', vkpts) )

    !> dump Wannier projection matrices to h5
    !> store only part that is non-zero
    do isp = 1, ispin
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       VH5_CHECK(vh5_group_open_or_create(groupid, trim(spingroupname), spin_groupid))
       do nk=1, nkpts
          write(kgroupname, '(2A)') "kpoint_"//trim(str(nk))
          VH5_CHECK(vh5_group_open_or_create(spin_groupid, trim(kgroupname), k_groupid))

          !> count number of relevant bands for this k-point
          ip=0; ipp=0
          do i=1,nb_tot
             if (exclude_bands(i)) cycle
             ip=ip+1
             if (.not.lwindow(ip,nk,isp)) cycle
             ipp=ipp+1
          enddo
          allocate( included_bands(ipp) )
          VH5_CHECK(vh5_write(k_groupid, 'nbands', ipp))
          dim_incl_bands=[ ipp ]

          !> and store them in array included_bands
          allocate( uij(ipp,nw) )
          ip=0; ipp=0
          do i=1,nb_tot
             if (exclude_bands(i)) cycle
             ip=ip+1
             if (.not.lwindow(ip,nk,isp)) cycle
             ipp=ipp+1
             included_bands(ipp) = i
             uij( ipp, 1:nw ) = u( i, 1:nw, nk, isp )
          enddo
          !VH5_CHECK(vh5_create_integer_array_nd(k_groupid, "incl_bands", 1, dim_incl_bands))
          VH5_CHECK(vh5_write(k_groupid, "incl_bands", included_bands) )
          !> store row by row of projection matrix
          count = [ ipp, nw ]
#ifdef gammareal
          VH5_CHECK(vh5_create_double_array_nd(k_groupid, "u", 2, count))
#else
          VH5_CHECK(vh5_create_double_complex_array_nd(k_groupid, "u", 2, count))
#endif
          VH5_CHECK(vh5_write( k_groupid, 'u', uij ) )

          deallocate( uij )
          deallocate( included_bands )
          VH5_CHECK(vh5_group_close_writing(k_groupid))
       enddo ! nkpts
       VH5_CHECK(vh5_group_close_writing(spin_groupid))
    enddo !isp

    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_wanproj

  !> @brief reads wanproj
  function vh5_read_wanproj(fileid, wdes, io, nw, uwan, group) result( ierr)
    USE base
    USE string
    USE wave
    USE full_kpoints, ONLY : KPOINTS_FULL, LIDENTICAL_KPOINT
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(wavedes), intent(in) :: wdes !< wave function descriptor
    type(in_struct) :: io
    integer, intent(out) :: nw !< number of wannier states read from file
    !> u_{nm}^k tranformation matrix between bloch and wannier functions
    GDEF, allocatable, intent(out) :: uwan(:,:,:,:)
    character(len=*),optional,intent(in) :: group
    integer :: ierr
    ! local variables
    character(len=MAX_LEN_GROUP) :: my_group
    integer :: iu0
    integer i,j, k, ip, ipp
    integer isp,nk
    ! variables that are read
    integer(HID_T) :: groupid, spingroupid, spin_groupid, k_groupid
    character(LEN=256) :: spingroupname, kgroupname
    integer            ispin, nkpts, nb_tot, nbands
    real(q)            vkpt(3)
    integer, allocatable :: included_bands(:)
    integer, allocatable :: kpt(:)
    GDEF, allocatable :: uij(:,:)
    real(q), allocatable :: vkpts(:,:)
    integer :: node_me, ionode
    !----------------------------------------------------------------------
    node_me=0
    ionode=0
#ifdef MPI
    node_me=wdes%comm%node_me
    ionode =wdes%comm%ionode
#endif

    my_group = GRP_WANPROJ;     if (present(group))    my_group = group

    iu0 = io%iu0

    ierr= vh5_check_group_presence( fileid, wdes%comm, iu0, trim( my_group ) )
    if( ierr /= 0 ) return
    VH5_CHECK(vh5_group_open(fileid, trim(my_group), groupid, wdes%comm))

    if (iu0 >= 0) then
       write(iu0,*) 'found vaspwave.h5, reading file, group = ', trim(my_group)
    endif

    VH5_CHECK(vh5_read(groupid, "ispin", ispin, wdes%comm))
    VH5_CHECK(vh5_read(groupid, "nkpts", nkpts, wdes%comm))
    VH5_CHECK(vh5_read(groupid, 'nb_tot', nb_tot, wdes%comm))
    VH5_CHECK(vh5_read(groupid, 'nw', nw, wdes%comm))

    !maybe we need to kill the job here
    if (nkpts < wdes%nkpts) then
       call vtutor%error('inconsistent nkpts in '//trim(my_group)//str(nkpts)//' '//str(wdes%nkpts))
    endif
    if (ispin < wdes%ispin) then
       call vtutor%error('inconsistent nkpts in '//trim(my_group)//str(ispin)//' '//str(wdes%ispin))
    endif
    if (nb_tot <  wdes%nb_tot) then
       call vtutor%alert('wannier projection for bands beyond no. '//str(nb_tot)//&
          ' not found in wanproj. i will ignore those bands in the following.')
    endif

    allocate(uwan(wdes%nb_tot,nw,wdes%nkpts,wdes%ispin))
    uwan = 0

    !write info on band reduction
    if ( nb_tot > wdes%nb_tot .and. iu0>=0 ) &
       write(*,'(" downsampling from ",i4," -> ",i4, " bands")')nb_tot, wdes%nb_tot
    !write info on k-points downsampling
    if ( nkpts > wdes%nkpts  .and. iu0>=0 ) &
       write(*,'(" downsampling from ",i4," -> ",i4, " k-points")')nkpts, wdes%nkpts
    !write info on spin reduction
    if ( ispin > wdes%ispin  .and. iu0>=0 ) &
       write(*,'(" downsampling from ",i4," -> ",i4, " spins")')ispin, wdes%ispin

    !> filter relevant k-points from list stored in file
    allocate( kpt( wdes%nkpts ) )
    kpt=0
    io_begin
       k=0
       allocate( vkpts( 3, nkpts ) )
       VH5_CHECK(vh5_read(groupid, "vkpts", vkpts))
       ! search correct kpoint list
       do i=1, nkpts
          do nk = 1, wdes%nkpts
                !check if current k-point of list corresponds to k-point of file
                IF (.NOT.LIDENTICAL_KPOINT(WDES%VKPT(:,NK),VKPTS(:,i))) CYCLE
                k=k+1
                kpt( k ) = i
          enddo
       enddo
       deallocate( vkpts )
    io_end
    CALLMPI( M_bcast_i( wdes%COMM, k,1) )
    CALLMPI( M_bcast_i( wdes%COMM, kpt, wdes%nkpts))

    if( k /= wdes%nkpts ) then
       call vtutor%error('vaspwave.h5 dataset '//GRP_WANPROJ//' does not contain all k-points: '//&
          str( [ k, wdes%nkpts]) )
    endif

    !> read relevant projection matrices
    do isp = 1, min( ispin, wdes%ispin )
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       VH5_CHECK(vh5_group_open(groupid, trim(spingroupname), spin_groupid, wdes%comm))

       do k = 1, wdes%nkpts
             ! nk is the k-point in file that stores requested k-point
             nk = kpt( k )
             write(kgroupname, '(2A)') "kpoint_"//trim(str(nk))
             VH5_CHECK(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid, wdes%comm))

             ! restrict reading to ionode
             io_begin
             ! read how many bands are projected
             VH5_CHECK(vh5_read(k_groupid, "nbands", nbands))
             ! this array will store the relevant bands that are projected
             allocate( included_bands( nbands ) )
             ! this is the non-zero part of the projection matrix from bloch to wannier basis
             allocate( uij( nbands, nw ) )
             ! read relevant entries on file
             VH5_CHECK(vh5_read(k_groupid, "incl_bands", included_bands))
             VH5_CHECK(vh5_read(k_groupid, "u", uij)  )
             ! store these entries in the actual projection matrix
             ! that is used in LINCOM
             do j = 1, nw
                do i = 1, nbands
                    uwan( included_bands( i ), j, k, isp ) = uij(i,j)
                enddo
             enddo
             ! remove temporary storage from heap
             deallocate( uij )
             deallocate( included_bands )
             io_end
             ! all ranks can now close this group
             VH5_CHECK(vh5_group_close(k_groupid, wdes%comm))
       enddo
       VH5_CHECK(vh5_group_close(spin_groupid, wdes%comm))
    enddo
    deallocate( kpt )
    VH5_CHECK(vh5_group_close(groupid, wdes%comm))

#ifdef gammareal
    CALLMPI(M_sum_d8(wdes%comm,uwan,size(uwan,kind=qi8)))
#else
    CALLMPI(M_sum_g8(wdes%comm,uwan,size(uwan,kind=qi8)))
#endif

  end function vh5_read_wanproj

  !> @brief reads PAW occupancies from a dataset group
  function vh5_read_rho_paw(groupid, my_group, p, t_info, loverl, rholm_store, ispin, comm, io) result( ierr)
    USE base
    USE string
    USE pseudo
    USE poscar
    USE wave
    USE paw
    USE constant
    integer(HID_T), intent(in) :: groupid !< hdf5 file handle
    character(len=*), intent(in) :: my_group !< group name
    type(type_info), intent(in) :: t_info
    type(potcar), intent(in),target  ::  p (t_info%ntyp)
    logical, intent(in) :: loverl
    real(q), intent(inout) :: rholm_store(:,:)   !< storage for the channel occupancies
    integer, intent(in) :: ispin !< selected spin
    type(communic), intent(in) :: comm
    type(in_struct), intent(in) :: io
    integer :: ierr
    !> local variables
    integer :: iu0
    character(LEN=40) :: datasetname
    integer, allocatable :: nelements(:)
    integer :: ni
    integer :: lymax, lmmax
    real(q), allocatable :: buffer(:,:)
    type(potcar), pointer :: pp
    integer :: ibase, iadd, i_error
    integer :: nelements_read
    integer, external :: maxl_aug
    !-------------------------------------------------------------------------------------

    iu0 = io%iu0
    if (.not.loverl .or. mimic_us ) return

    lymax = maxl_aug( t_info%ntyp, p )
    lmmax = (lymax+1)**2

    write(datasetname,'(A,I0.2)') "aug_occupancies_ions_s",ispin
    ierr= vh5_check_group_presence( groupid, comm, iu0, trim(datasetname), parent_group = my_group )
    if( ierr /= 0 ) return

    ! read number of elements and number of species
    allocate(nelements(t_info%nions))
    ierr = vh5_read( groupid, trim(datasetname), nelements, comm )
    if( ierr /=0 ) then
       deallocate( nelements )
       call vtutor%error( "Reading "//trim(my_group)//"/"//trim(datasetname)//" failed with error code: "//str(ierr) )
    endif

    write(datasetname,'(A,I0.2)') "aug_occupancies_s",ispin
    ierr= vh5_check_group_presence( groupid, comm, iu0, trim(datasetname), parent_group = my_group )
    if( ierr /= 0 ) return

    allocate( buffer( lmmax*lmmax, t_info%nions ) )
    buffer = 0
    VH5_CHECK(vh5_read( groupid, trim(datasetname), buffer, comm ))

    ! loop over all ions
    ibase = 1
    ion: do ni=1,t_info%nions
         iadd = 0

         pp => p( t_info%ityp( ni ) )

         if (nelements(ni)> lmmax*lmmax) then
            CALL vtutor%bug("running out of buffer "//str([ nelements(ni), lmmax**2, ni] ), __FILE__, __LINE__)
         endif

         nelements_read = 0
         if( do_local( ni ) ) then
            call store_rholm( buffer(:,ni), rholm_store(ibase:,ispin), &
                 metric(ibase:), iadd, pp, .true. , nelements_read)
            ibase = ibase + iadd
         endif

         CALLMPI( MPI_barrier( comm%mpi_comm, i_error ))
         CALLMPI( M_sum_i(comm, nelements_read, 1))

         IF (nelements(ni) /= nelements_read)  THEN
           ierr=1
         ENDIF

         if ( ierr /=0 ) then
            if( iu0>= 0 ) write(*,*) "Corrupt data in "//trim(my_group)//"/"//trim(datasetname)//" found for ion "//str(ni)
            deallocate( nelements )
            deallocate( buffer )
            return
         endif
    enddo ion


    ! now try to find out whether or not the imaginary part of the occupancies were written to h5
    write(datasetname,'(A,I0.2)') "imag_aug_occupancies_ions_s",ispin
    ierr= vh5_check_group_presence( groupid, comm, -1, trim(datasetname), parent_group = my_group )
    if( ierr /= 0 ) then
       deallocate( nelements )
       deallocate( buffer )
       ierr = 0
       return
    endif

    !
    ! continue with imaginary part
    !
    if( ierr == 0 ) then
       if( iu0>=0 ) write(*,*)" reading imaginary PAW occupancies from "//trim(my_group)//"/"//trim(datasetname)
    endif

    ! read number of elements and number of species
    ierr = vh5_read( groupid, trim(datasetname), nelements, comm )
    if( ierr /=0 ) then
       call vtutor%error( "Reading "//trim(my_group)//"/"//trim(datasetname)//" failed with error code: "//str(ierr) )
    endif

    write(datasetname,'(A,I0.2)') "aug_occupancies_s",ispin
    ierr= vh5_check_group_presence( groupid, comm, iu0, trim(datasetname), parent_group = my_group )
    if( ierr /= 0 ) return

    ibase =1
    ion_imag: do ni=1,t_info%nions
         iadd = 0

         pp => p( t_info%ityp( ni ) )

         if (nelements(ni)> lmmax*lmmax) then
            CALL vtutor%bug("running out of buffer "//str([ nelements(ni), lmmax**2, ni] ), __FILE__, __LINE__)
         endif

         CALLMPI( M_bcast_d( comm, buffer(1,ni), nelements(ni) ) )

         nelements_read = 0
         if( do_local( ni ) ) then
            call store_rholm( buffer(:,ni), rholm_store(ibase+ubound(rholm_store,1)/2:,ispin), &
                 metric(ibase+ubound(rholm_store,1)/2:), iadd, pp, .true. , nelements_read)
            ibase = ibase + iadd
         endif

         CALLMPI( MPI_barrier( comm%mpi_comm, i_error ))
         CALLMPI( M_sum_i(comm, nelements_read, 1))

         IF (nelements(ni) /= nelements_read)  THEN
           ierr=1
         ENDIF

         if ( ierr /=0 ) then
            if( iu0>= 0 ) write(*,*) "Corrupt data in "//trim(my_group)//"/"//trim(datasetname)//" found for ion "//str(ni)
            deallocate( buffer )
            return
         endif
    enddo ion_imag


    deallocate( nelements )
    deallocate( buffer )

  end function vh5_read_rho_paw

  !> @brief reads charge  dataset (or oep)
  function vh5_inchg(fileid, gridc, ng, cvtot, isp, iu0, group, comm) result(ierr)
    USE base
    USE string
    USE wave
    USE pseudo
    USE paw
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(grid_3d), intent(in)  :: gridc
    integer, intent(in) :: ng(3) !< grid dimensions
    complex(q), intent(out) ::  cvtot(:,:)
    integer, intent(in) :: isp !< spin channel to read
    integer, intent(in) :: iu0 !< unit number for stdout
    character(len=*), intent(in), optional :: group !< group name (default: GRP_POTENTIAL)
    type(communic), intent(in) :: comm
    ! local variables
    integer :: ierr
    integer :: nz
    real(q), pointer  :: work(:)
    character(len=MAX_LEN_GROUP) :: my_group, datasetname
    integer :: dims(4)
    integer :: rank
    integer(HSIZE_T) :: start(4), count(4)
    integer(HID_T) :: groupid
    real(q), allocatable :: pot_plane(:,:)
    !-------------------------------------------------------------------------------------
    my_group    = GRP_POT;  if (present(group))    my_group    = trim(group)

    ierr=vh5_group_open_or_create(fileid, trim(my_group), groupid)
    if( ierr /=0 ) then
       if( iu0>=0 ) write(iu0,*) "Opening group "//trim(my_group)//" failed!"
       return
    endif

    ierr=vh5_get_dimensions_par(groupid, "total", dims, rank, comm)
    if( ierr /=0 ) then
       if( iu0>=0 ) write(iu0,*) "Reading dataset dimensions of "//trim(my_group)//"/total failed!"
       ierr = 1
       goto 100
    endif

    ! sanity check
    if (rank /= 4) then
       if( iu0>=0 ) write(iu0,*) "Dataset "//trim(my_group)//" has wrong rank "//str(rank)//" (expected 4)!"
       ierr = 1
       goto 100
    endif
    if (ng(1)==0 .and. ng(2)==0 .and. ng(3)==0) then
       if( iu0>=0 ) write(iu0,*) "Dataset "//trim(my_group)//" contains no points!"
       ierr = 1
       goto 100
    else if (ng(1)/=gridc%ngx .or. ng(2)/=gridc%ngy .or. ng(3)/=gridc%ngz ) then
       if( iu0>=0 ) write(iu0,*) "Dataset "//trim(my_group)//" sampled with other grid:"//str(ng)//" vs "&
          //str( [gridc%ngx, gridc%ngy, gridc%ngz ] )
       ierr = 1
       goto 100
    else if( ng(1) /= dims( 1) .or. ng(2) /= dims(2) .or. ng(3) /= dims(3) ) then
       if( iu0>=0 ) write(iu0,*) "Dataset "//trim(my_group)//" sampled with other grid:"//str(ng)//" vs "&
          //str( dims(1:3) )
       ierr = 1
       goto 100
    endif

    allocate( pot_plane(gridc%ngx, gridc%ngy) )

    do nz = 1, gridc%ngz
       if (comm%node_me == comm%ionode)  then
           start = [1, 1, nz, isp]
           count = [gridc%ngx, gridc%ngy, 1, 1]
           ! read potential
           ierr =  vh5_read_double_subarray_nd(groupid, "total", 4, start, count, pot_plane )
           if( ierr /=0 ) then
              if( iu0>=0 ) write(iu0,*) "Reading dataset "//trim(my_group)//"/total failed!"
              goto 100
           endif
       endif
       ! distribute ngx x ngz plane wave coefficients
       call dis_grid_rl_plane(gridc, pot_plane, cvtot(1,isp), .true., nz )
    enddo
    ! fft from real to reciprocal space
    call fft_rc_scale(cvtot(1,isp),cvtot(1,isp),gridc)
    call setunb_compat(cvtot(1,isp),gridc)

100 continue
    if (vh5_group_close_writing( groupid )/=0) then
       if( iu0>=0 ) write(iu0,*) "closing group "//trim(my_group)//" failed!"
    endif
    if (allocated(pot_plane)) deallocate( pot_plane )
  end function vh5_inchg

  !> @brief reads pot  dataset (or oep)
  function vh5_read_pot(fileid, gridc, t_info, cvtot, wdes, latt_cur, p, lmdim, cdij, n_mix_paw, io, group) result( ierr)
    USE base
    USE string
    USE wave
    USE pseudo
    USE paw
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(grid_3d), intent(in)  :: gridc
    type(type_info), intent(in) :: t_info
    type(potcar), intent(in) ::  p (t_info%ntyp)
    type(latt), intent(in) ::  latt_cur
    type(wavedes), intent(in) :: wdes
    complex(q), intent(out) ::  cvtot(gridc%mplwv,wdes%ncdij)
    integer, intent(in) :: lmdim
    OVERLAP, intent(out) ::  cdij(lmdim,lmdim,wdes%nions,wdes%ncdij)
    integer, intent(in) ::  n_mix_paw
    type(in_struct) :: io
    character(len=*), intent(in), optional :: group !< group name (default: GRP_POTENTIAL)
    integer :: ierr
    ! local variables
    real(q) :: scale, amat(3,3)
    integer :: ng(3), nz, isp
    real(q), pointer  :: work(:)
    real(q)         :: dlm(n_mix_paw,wdes%ncdij)
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group
    integer :: iu0
    real(q), parameter :: tiny=1e-6_q
    integer :: dims(3)

    !-------------------------------------------------------------------------------------
    my_group    = GRP_POT;  if (present(group))    my_group    = trim(group)

    iu0 = io%iu0
    ! check if group is present:
    ierr= vh5_check_group_presence( fileid, wdes%comm, iu0, trim( my_group ) )
    if( ierr /= 0 ) return

    if( iu0>=0 ) write(*,*)'reading vaspwave.h5, group = '//trim(my_group )
    ! open group
    VH5_CHECK(vh5_group_open(fileid, trim(my_group), groupid, wdes%comm))

    ! this is kind of a sanity check if dataset is consistent with present setting
    VH5_CHECK(vh5_group_open(groupid, 'position', subgroupid, wdes%comm))
    VH5_CHECK(vh5_read(subgroupid, "scale", scale, wdes%comm))
    VH5_CHECK(vh5_read(subgroupid, "lattice_vectors", amat, wdes%comm))
    VH5_CHECK(vh5_group_close(subgroupid, wdes%comm))

    if( norm2( latt_cur%a - amat*scale ) > tiny ) then
       if( iu0>=0 ) write(*,*)"WARNING: lattice in dataset "//trim(my_group)&
          //" seems to have changed by: "//str(norm2( latt_cur%a - amat*scale ))
    endif

    ! read grid dimensions
    VH5_CHECK(vh5_read(groupid, "grid", ng, wdes%comm))
    ierr = vh5_inchg(fileid, gridc, ng, cvtot, 1, iu0, trim(my_group), wdes%comm)
    if (ierr /=0 ) then
       VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    endif

    CALL FFT3D(CVTOT(1,1),GRIDC,1)

    ! just in case initialize the magnetization density to 0
    do isp=2,wdes%ncdij
      call rc_add(cvtot(1,isp),0.0_q,cvtot(1,isp),0.0_q,cvtot(1,isp),gridc)
    enddo

    dlm(:,1)=0
    ! read PAW occupancies like in RD_RHO_PAW
    ierr = vh5_read_rho_paw( groupid, trim(my_group), p, t_info, wdes%loverl, dlm, 1, wdes%comm, io )

    if( ierr /=0 ) then
       if( iu0>=0 ) write(*,*)"WARNING: paw-occupancies not read correctly for isp=1"
          VH5_CHECK(vh5_group_close(groupid, wdes%comm))
       return
    endif

    ! read spin density
    spin: do isp = 2, wdes%ncdij
       ierr = vh5_inchg(fileid, gridc, ng, cvtot, isp, iu0, trim(my_group), wdes%comm)
       if (ierr /=0 ) then
          VH5_CHECK(vh5_group_close(groupid, wdes%comm))
          return
       endif
       CALL FFT3D(CVTOT(1,isp),GRIDC,1)
       ! read PAW occupancies like in RD_RHO_PAW
       ierr = vh5_read_rho_paw( groupid, trim(my_group), p, t_info, wdes%loverl, dlm, isp, wdes%comm, io )
       if( ierr /=0 ) then
          if( isp==2 ) then
             ! we could not read any entry in the magnetisation density return ICHARG=-1
             ! (magnetisation according to overlapping atoms)
             ierr = -1
             if( iu0>=0 ) write(*,*)"WARNING: failed to read magnetisation density from vaspwave.h5"
             VH5_CHECK(vh5_group_close(groupid, wdes%comm))
             return
          else
             call rc_add(cvtot(1,2),1.0_q,cvtot(1,2),0.0_q,cvtot(1,4),gridc) ! copy
             call rc_add(cvtot(1,2),0.0_q,cvtot(1,2),0.0_q,cvtot(1,2),gridc) ! =0
             call rc_add(cvtot(1,3),0.0_q,cvtot(1,3),0.0_q,cvtot(1,3),gridc) ! =0

             dlm(:,4)=dlm(:,2)
             dlm(:,2)=0
             dlm(:,3)=0
             exit spin
          endif
       endif
    enddo spin

    ! one centre corrections are passed through mixer and added now
    call retrive_rho_paw(wdes, p, t_info, wdes%loverl, lmdim, cdij, dlm)

    ! successful exit
    ierr = 0
    VH5_CHECK(vh5_group_close(groupid, wdes%comm))
  end function vh5_read_pot

  !*******************************************************************
  !> write NICS and POSNICS data to the hdf5 file
  !*******************************************************************
  subroutine vh5_write_nics_for_posnics(fileid, subgroup, nicsdata_sym, taunics, npos_nics, lnmrshield)
    use prec, only: q
    use string, only: str
    use tutor, only: vtutor
    integer(HID_T), intent(in) ::  fileid   !< id of h5 file
    real(q), intent(in) :: nicsdata_sym(:,:,:)
    real(q), intent(in) :: taunics(:,:)
    integer, intent(in) :: npos_nics
    logical, intent(in) :: lnmrshield
    character(len=*), intent(in) :: subgroup     !< name of the dataset in the HDF5 file
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: subgroupid
    character(len=44) :: label_
    !-------------------------------------------------------------------------------------
    ! check for consistency
    if (size(taunics,1) /= npos_nics) then
       call vtutor%bug('internal ERROR: vh5_write_nics_for_posnics size(taunics,1) = ' // &
         str(size(taunics,1)) // ' is not compatible with npos_nics = ' // str(npos_nics) // &
         '.', __FILE__, __LINE__)
    end if
    if (size(nicsdata_sym,3) < npos_nics) then
       call vtutor%bug('internal ERROR: vh5_write_nics_for_posnics',__FILE__, __LINE__)
    end if

    ! open groups
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), subgroupid))

    ! write positions provided in POSNICS file
    VH5_CHECK(vh5_write(subgroupid, "positions", taunics))

    ! write NICS
    VH5_CHECK(vh5_write(subgroupid, "values", nicsdata_sym(:,:,1:npos_nics)))

    if (lnmrshield) then
      label_="Nuclear independent chemical shielding (ppm)"
    else
      label_="Nuclear independent chemical shift (ppm)"
    endif
    VH5_CHECK(vh5_write(subgroupid, "label", trim(label_)))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_file_flush(fileid))

  end subroutine vh5_write_nics_for_posnics

end module vhdf5
#else
!--------------------------------------------------------------------
! Dummy module
!--------------------------------------------------------------------
module vhdf5
   use vhdf5_base
   contains
   subroutine vasphdf5_dummy
      write(*,*)'Im a DEC compiler so I need this line'
   end subroutine vasphdf5_dummy
end module vhdf5
#endif
