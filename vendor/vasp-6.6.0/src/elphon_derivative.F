#include "symbol.inc"

module elphon_derivative
   use prec
   use base, only: in_struct, info_struct
   use msymmetry, only: prim_cell_t
   use poscar_struct_def, only: type_info
   use wave_struct_def, only: wavedes
   use mgrid_struct_def, only: grid_3d, transit
   use tutor, only: vtutor

#ifdef VASP_HDF5
   use hdf5, only: hid_t
   use vhdf5_base, only: &
      vh5_file_create_or_overwrite, &
      vh5_file_close_writing, &
      vh5_write, vh5_read, &
      vh5_error
#endif

   implicit none
   private

   public elphon_deriv_settings, elphon_deriv_reader
   public elphon_deriv_io

   public calc_elphon_derivatives, write_elphon_checkpoint, read_elphon_checkpoint

   type elphon_deriv_settings
      logical :: generate = .false.
      real(q) :: prim_latt(3, 3)
      logical :: use_prim_latt = .false.
      integer :: fft_mesh(3) = [0, 0, 0]
   contains
      procedure, public :: check => elphon_deriv_check
   end type

   type grid_interpolator_t
      type(grid_3d), public :: source_grid
      type(grid_3d), public :: target_grid

      real(q), private :: recip_target_lattice(3, 3)
      real(q), private :: recip_source_lattice(3, 3)
      complex(q), allocatable, private :: buffer(:)
      integer, private :: target_to_source(3, 3)
      logical, private :: is_aligned
      type(transit), private :: aligned_transit
      integer, allocatable, private :: general_transit(:)
      integer, allocatable, private :: plane_wave_weight(:)
   contains
      procedure, public :: init => grid_interpolator_init
      procedure, public :: interpolate => grid_interpolator_interpolate
      procedure, private :: set_general_transit => grid_interpolator_set_general_transit

      ! _test_mani write finalizer
   end type

   type min_image_mapper_t
      real(q) :: cart_to_source(3, 3)
      real(q) :: standard_to_cart(3, 3)
      real(q) :: source_to_standard(3, 3)
      real(q) :: shift_vectors(3, 8)
   contains
      procedure, public :: init => min_image_mapper_init
      procedure, public :: get_shift => min_image_mapper_get_shift
   end type

contains

   subroutine min_image_mapper_init(this, source_lattice)
      class(min_image_mapper_t), intent(out) :: this
      real(q), intent(in) :: source_lattice(3, 3)
      real(q) :: std_latt(3, 3)
      integer :: shift_idx
      integer :: x, y, z
      integer :: i_dum
      real(q) :: r_dum(6)

      std_latt = source_lattice
      call lattyp(std_latt(:, 1), std_latt(:, 2), std_latt(:, 3), i_dum, r_dum, -1)
      this%standard_to_cart = std_latt
      this%cart_to_source = matrix_inverse_3x3_r(source_lattice)
      this%source_to_standard = matmul(matrix_inverse_3x3_r(std_latt), source_lattice)

      shift_idx = 0
      do z = -1, 0
         do y = -1, 0
            do x = -1, 0
               shift_idx = shift_idx + 1
               this%shift_vectors(:, shift_idx) = matmul(std_latt, [x, y, z])
            enddo
         enddo
      enddo
   end subroutine

   function min_image_mapper_get_shift(this, point, square_dist) result(shift)
      integer :: shift(3)
      class(min_image_mapper_t), intent(in) :: this
      real(q), intent(in) :: point(3)
      real(q), optional, intent(out) :: square_dist

      real(q) :: trafo_point(3)
      real(q) :: shifted_point(3)
      real(q) :: min_dist, cur_dist
      integer :: shift_idx
      integer :: min_idx

      trafo_point = matmul(this%source_to_standard, point)
      trafo_point = modulo(trafo_point, 1._q)
      trafo_point = matmul(this%standard_to_cart, trafo_point)

      min_dist = huge(1._q)
      do shift_idx = 1, 8
         shifted_point = trafo_point + this%shift_vectors(:, shift_idx)
         cur_dist = dot_product(shifted_point, shifted_point)
         if (cur_dist < min_dist) then
            min_dist = cur_dist
            min_idx = shift_idx
         endif
      enddo

      shifted_point = trafo_point + this%shift_vectors(:, min_idx)
      shifted_point = matmul(this%cart_to_source, shifted_point)
      shift = nint(shifted_point - point)

      if (present(square_dist)) square_dist = min_dist
   end function

   subroutine elphon_deriv_reader(incar, settings)
      use incar_reader, only: incar_file, process_incar, read_type_error

      type(incar_file), intent(inout) :: incar
      type(elphon_deriv_settings), intent(out) :: settings

      integer :: ierr

      call process_incar(incar, "elph_pot_generate", settings%generate)
      call process_incar(incar, "elph_pot_fft_mesh", settings%fft_mesh)
      call process_incar(incar, "elph_pot_lattice", settings%prim_latt, ierr)

      if (ierr == 0) then
         settings%use_prim_latt = .true.
      else if (ierr == read_type_error) then
         call vtutor%error("Error reading ELPH_POT_LATTICE: " // &
            "expected 9 real numbers for the primitive lattice vectors.")
      endif

#ifndef VASP_HDF5
      if (settings%generate) call vtutor%error( &
         "In order to use ELPH_POT_GENERATE=True, VASP needs to be compiled with HDF5 support.")
#endif
   end subroutine

   subroutine elphon_deriv_check(this, ibrion, nfree, prim_cell)
      class(elphon_deriv_settings), intent(in) :: this
      integer, intent(in) :: ibrion
      integer, intent(in) :: nfree
      type(prim_cell_t), intent(in) :: prim_cell

      type(prim_cell_t) :: tmp_prim_cell
      character(len=:), allocatable :: error_msg

      if (.not. this%generate) return

      if (ibrion /= 6) call vtutor%error( &
         "ELPH_POT_GENERATE=True requires IBRION=6 for finite differences at the moment. " // &
         "Please set IBRION=6 in your INCAR file.")

      if (nfree /= 0 .and. nfree /= 2) call vtutor%error( &
         "ELPH_POT_GENERATE=True requires NFREE=2 for central differences at the moment. " // &
         "Please omit NFREE in your INCAR file or set NFREE=2.")

      if (this%use_prim_latt) then
         call tmp_prim_cell%init_from_latt(prim_cell, this%prim_latt, error_msg)
         if (allocated(error_msg)) call vtutor%error( &
            "The primitive cell provided by ELPH_POT_LATTICE is invalid: " // error_msg)
      endif
   end subroutine

   subroutine elphon_deriv_io(io, settings, prim_cell, info, t_info, wdes, &
      Dij, Qij, V_grid, dV_du, dDij_du, force_const, polar_data, dmu_du)
      use elphon_util, only: polar_data_t
      use poscar, only: outpos
#ifdef VASP_HDF5
      use vhdf5, only: ih5outfileid, vh5_write_lattice_ion
#endif

      type(in_struct), intent(in) :: io
      type(elphon_deriv_settings), intent(in) :: settings
      type(prim_cell_t), target, intent(in) :: prim_cell
      type(info_struct), intent(in) :: info
      type(type_info), intent(in) :: t_info
      type(wavedes), intent(in) :: wdes
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      OVERLAP, intent(in) :: Qij(:, :, :, :)
      type(grid_3d), intent(in) :: V_grid
      complex(q), intent(in) :: dV_du(:, :, :, :)
      OVERLAP, intent(in) :: dDij_du(:, :, :, :, :, :)
      real(q), intent(in) :: force_const(:, :)
      type(polar_data_t), intent(in) :: polar_data
      complex(q), optional, intent(in) :: dmu_du(:, :, :, :)

      type(prim_cell_t), pointer :: prim_cell_ptr
      type(type_info) :: t_info_prim
      integer :: prim_fft_mesh(3)
      real(q), allocatable, target :: prim_positions(:, :)
      integer, allocatable, target :: prim_nityp(:)
      integer :: tmp
      integer :: unit

      ! Generate small headers for stdout and OUTCAR
      if (io%iu0 >= 0) write(io%iu0, *) "Generating phelel_params.hdf5"
      if (io%iu6 >= 0) then
         write(io%iu6, *)
         write(io%iu6, *) "Generation of phelel_params.hdf5"
         write(io%iu6, *) "================================"
         write(io%iu6, *)
      endif

      if (settings%use_prim_latt) then
         ! Initialize primitive cell from provided lattice
         allocate(prim_cell_ptr)
         call prim_cell_ptr%init_from_latt(prim_cell, settings%prim_latt)
      else
         ! Use internal primitive cell
         if (io%iu0 >= 0) then
            write(io%iu0, "(A)") "  No primitive cell provided"
            write(io%iu0, "(A)") "  Using primitive cell determined by VASP"
         endif
         prim_cell_ptr => prim_cell
      endif

      if (any(settings%fft_mesh == 0)) then
         ! Determine FFT grid dimensions of primitive cell from ENCUT
         call find_fft_dims(info, prim_cell_ptr%prim_latt, prim_fft_mesh)
         if (io%iu0 >= 0) then
            write(io%iu0, "(A)") "  No primitive FFT mesh was provided"
            write(io%iu0, "(A, 3I6)") "  Primitive FFT mesh automatically set to ", prim_fft_mesh
         endif
      else
         ! Use provided FFT grid dimensions
         prim_fft_mesh = settings%fft_mesh
         if (io%iu0 >= 0) write(io%iu0, "(A, 3I6)") "  Primitive FFT mesh was provided as ", prim_fft_mesh
      endif

      ! Write primitive-cell information to OUTCAR
      if (io%iu6 >= 0) then
         write(io%iu6, "(A)") "Primitive cell"
         write(io%iu6, "(A, 3F14.8)") "   a1 = ", prim_cell_ptr%prim_latt%a(:, 1)
         write(io%iu6, "(A, 3F14.8)") "   a2 = ", prim_cell_ptr%prim_latt%a(:, 2)
         write(io%iu6, "(A, 3F14.8)") "   a3 = ", prim_cell_ptr%prim_latt%a(:, 3)
         write(io%iu6, *)
         write(io%iu6, "(A, 3I6)") "Primitive FFT mesh = ", prim_fft_mesh
         write(io%iu6, *)
      endif

      ! Prepare output structures and arrays
      call prim_cell_ptr%get_prim_posion(prim_positions)
      call prim_cell_ptr%get_prim_nityp(t_info, prim_nityp)
      t_info_prim = t_info
      t_info_prim%nions = prim_cell_ptr%num_prim_atoms
      t_info_prim%posion => prim_positions
      t_info_prim%nityp => prim_nityp

      ! Write primitive-cell structure to POSCAR format
      if (io%iu0 >= 0) then
         open(newunit=unit, file="CONTCAR_ELPH", status="replace")
         call outpos(unit, .true., info%sznam1, t_info_prim, &
            prim_cell_ptr%prim_latt%scale, prim_cell_ptr%prim_latt%a, .false., &
            prim_positions, .false.)
         close(unit)
      endif

#ifdef VASP_HDF5
      ! Write primitive-cell structure to HDF5 output
      call vh5_write_lattice_ion(ih5outfileid, info%sznam1, t_info_prim, &
         prim_cell_ptr%prim_latt%scale, prim_cell_ptr%prim_latt%a, .false., &
         prim_positions, subgroup="elph_potential/primitive_positions")

      ! Write primitive-cell FFT grid dimensions to HDF5 output
      VH5_CHECK(vh5_write(ih5outfileid, "results/elph_potential/primitive_fft_mesh", prim_fft_mesh))

      call write_phelel_params( &
         io, prim_cell_ptr, t_info, wdes, prim_fft_mesh, Dij, Qij, &
         V_grid, dV_du, dDij_du, force_const, polar_data, dmu_du=dmu_du)
#endif

      if (.not. associated(prim_cell_ptr, prim_cell)) deallocate(prim_cell_ptr)
   end subroutine

#ifdef VASP_HDF5
   subroutine write_phelel_params(io, prim_cell, t_info, wdes, fft_mesh, &
      Dij, Qij, V_grid, dV_du, dDij_du, force_const, polar_data, filename, dmu_du)
      use elphon_util, only: polar_data_t

      type(in_struct), intent(in) :: io
      type(prim_cell_t), intent(in) :: prim_cell
      type(type_info), intent(in) :: t_info
      type(wavedes), intent(in) :: wdes
      integer, intent(in) :: fft_mesh(3)
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      OVERLAP, intent(in) :: Qij(:, :, :, :)
      type(grid_3d), intent(in) :: V_grid
      complex(q), intent(in) :: dV_du(:, :, :, :)
      OVERLAP, intent(in) :: dDij_du(:, :, :, :, :, :)
      real(q), intent(in) :: force_const(:, :)
      type(polar_data_t), intent(in) :: polar_data
      character(*), optional, intent(in) :: filename
      complex(q), optional, intent(in) :: dmu_du(:, :, :, :)

      integer(hid_t) :: file_id
      character(:), allocatable :: act_filename
      integer :: supercell_matrix(3, 3)
      integer :: super_target_to_source(3, 3)
      integer :: prim_source_to_target(3, 3)
      type(grid_interpolator_t) :: grid_interpol
      integer, allocatable :: lattice_points(:, :)
      real(q), allocatable :: grid_points(:, :)
      integer :: interpol_fft_mesh(3)

      act_filename = "phelel_params.hdf5"
      if (present(filename)) act_filename = filename

      supercell_matrix = nint(prim_cell%super_to_prim)

      call find_target_grid( &
         supercell_matrix, fft_mesh, interpol_fft_mesh, &
         super_target_to_source, prim_source_to_target)

      call grid_interpol%init(V_grid, prim_cell%super_latt%a, super_target_to_source, interpol_fft_mesh)

      call generate_lattice_points(supercell_matrix, lattice_points)
      call generate_grid_points(lattice_points, fft_mesh, grid_points)

      ! _test_mani warn about FFT grid change

      VH5_CHECK(vh5_file_create_or_overwrite(act_filename, file_id))

      call write_primitive_cell_info(io, file_id, prim_cell, t_info)
      call write_paw_overlaps(io, file_id, prim_cell, wdes, Dij, Qij, dDij_du)
      call write_local_pot_derivative(io, file_id, &
         grid_interpol, grid_points, supercell_matrix, super_target_to_source, dV_du, "dVdu")
      if (present(dmu_du)) then
         call write_local_pot_derivative(io, file_id, &
            grid_interpol, grid_points, supercell_matrix, super_target_to_source, dmu_du, "dmudu")
      endif
      call write_force_constants(io, file_id, force_const)
      call write_polar_data(io, file_id, polar_data)
      call write_grid_dimensions(io, file_id, &
         lattice_points, grid_points, fft_mesh, prim_cell%prim_to_super)

      VH5_CHECK(vh5_file_close_writing(file_id))
   end subroutine

   subroutine write_primitive_cell_info(io, group_id, prim_cell, t_info)
      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(prim_cell_t), intent(in) :: prim_cell
      type(type_info), intent(in) :: t_info

      integer, allocatable :: prim_indices(:)
      real(q), allocatable :: prim_positions(:, :)

      if (io%iu0 < 0) return

      allocate(prim_indices(prim_cell%num_prim_atoms))
      prim_indices = prim_cell%super_idx(1, :)

      VH5_CHECK(vh5_write(group_id, "primitive_lattice", transpose(prim_cell%prim_latt%a)))
      VH5_CHECK(vh5_write(group_id, "supercell_lattice", transpose(prim_cell%super_latt%a)))
      VH5_CHECK(vh5_write(group_id, "unitcell_lattice", transpose(prim_cell%prim_latt%a)))
      VH5_CHECK(vh5_write(group_id, "supercell_matrix", transpose(nint(prim_cell%super_to_prim))))
      VH5_CHECK(vh5_write(group_id, "supercell_positions", prim_cell%posion_super))
      VH5_CHECK(vh5_write(group_id, "p2s_map", prim_indices-1))
      VH5_CHECK(vh5_write(group_id, "atom_indices_in_derivatives", prim_indices-1))
      VH5_CHECK(vh5_write(group_id, "primitive_masses", t_info%pomass(t_info%ityp(prim_indices))))

      allocate(prim_positions(3, prim_cell%num_prim_atoms))
      prim_positions = matmul(prim_cell%super_to_prim, prim_cell%posion_super(:, prim_indices))
      VH5_CHECK(vh5_write(group_id, "primitive_positions", prim_positions))
   end subroutine

   subroutine write_paw_overlaps(io, group_id, prim_cell, wdes, Dij, Qij, dDij_du)
      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(prim_cell_t), intent(in) :: prim_cell
      type(wavedes), intent(in) :: wdes
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      OVERLAP, intent(in) :: Qij(:, :, :, :)
      OVERLAP, intent(in) :: dDij_du(:, :, :, :, :, :)

      integer :: num_chan, ncdij
      integer :: num_prim_atoms, num_displacements
      integer :: total_size
      complex(q), allocatable :: overlap_out(:, :, :, :)
      complex(q), allocatable :: overlap_deriv_out(:, :, :, :, :, :)

      num_chan = size(Dij, 1)
      ncdij = size(Dij, 4)
      num_displacements = size(dDij_du, 5)
      num_prim_atoms = size(dDij_du, 6)

      allocate(overlap_out(num_chan, num_chan, prim_cell%num_prim_atoms, wdes%ncdij))
      total_size = size(overlap_out)

      call gather_overlap(Dij)
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "Dij", overlap_out))
      call gather_overlap(Qij)
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "qij", overlap_out))

      deallocate(overlap_out)
      allocate(overlap_deriv_out( &
         num_chan, num_chan, prim_cell%num_super_atoms, &
         num_displacements, num_prim_atoms, wdes%ncdij))

      call gather_overlap_derivative(dDij_du)
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "dDijdu", overlap_deriv_out))
      overlap_deriv_out = 0
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "dqijdu", overlap_deriv_out))

   contains

      subroutine gather_overlap(local_overlap)
         use wave, only: ni_global
         OVERLAP, target, intent(in) :: local_overlap(:, :, :, :)

         integer :: local_idx, global_idx
         integer :: prim_idx

         overlap_out = 0
         do prim_idx = 1, prim_cell%num_prim_atoms
            do local_idx = 1, wdes%nions
               global_idx = ni_global(local_idx, wdes%comm_inb)
               if (global_idx == prim_cell%super_idx(1, prim_idx)) then
                  overlap_out(:, :, prim_idx, :) = local_overlap(:, :, local_idx, :)
               endif
            enddo
         enddo
         CALLMPI(M_sum_z(wdes%comm_inb, overlap_out, size(overlap_out)))
      end subroutine

      subroutine gather_overlap_derivative(local_overlap)
         use wave, only: ni_global
         OVERLAP, target, intent(in) :: local_overlap(:, :, :, :, :, :)

         integer :: local_idx, global_idx
         integer :: icdij

         overlap_deriv_out = 0
         do icdij = 1, wdes%ncdij
            do local_idx = 1, wdes%nions
               global_idx = ni_global(local_idx, wdes%comm_inb)
               overlap_deriv_out(:, :, global_idx, :, :, icdij) = &
                  local_overlap(:, :, local_idx, icdij, :, :)
            enddo
         enddo
         CALLMPI(M_sum_z(wdes%comm_inb, overlap_deriv_out, size(overlap_deriv_out)))
      end subroutine
   end subroutine

   subroutine write_local_pot_derivative(io, group_id, &
      grid_interpol, grid_points, supercell_matrix, super_target_to_source, dV_du, groupname)

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(grid_interpolator_t), intent(inout) :: grid_interpol
      real(q), intent(in) :: grid_points(:, :)
      integer, intent(in) :: supercell_matrix(3, 3)
      integer, intent(in) :: super_target_to_source(3, 3)
      complex(q), intent(in) :: dV_du(:, :, :, :)
      character(len=*),  intent(in) :: groupname

      integer :: num_prim_atoms, ncdij
      integer :: prim_idx, icdij, cart_idx
      complex(q), allocatable :: interpol_buffer(:)
      complex(q), allocatable :: dV_du_out(:, :, :, :)

      num_prim_atoms = size(dV_du, 4)
      ncdij = size(dV_du, 2)

      allocate(interpol_buffer(grid_interpol%target_grid%mplwv))
      allocate(dV_du_out(size(grid_points, 2), 3, num_prim_atoms, ncdij))

      ! Iterate over all degrees of freedom
      do prim_idx = 1, num_prim_atoms
         do cart_idx = 1, 3
            do icdij = 1, ncdij
               ! Interpolate the data to the interpolation grid
               call grid_interpol%interpolate( &
                  dV_du(:, icdij, cart_idx, prim_idx), interpol_buffer)
               ! Transfer potential to output grid at relevant grid points
               call transfer_to_grid_points( &
                  grid_interpol%target_grid, grid_points, supercell_matrix, super_target_to_source, &
                  interpol_buffer, dV_du_out(:, cart_idx, prim_idx, icdij))
            enddo
         enddo
      enddo

      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, groupname, dV_du_out))
   end subroutine

   subroutine write_force_constants(io, group_id, force_const)
      use string, only: str

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      real(q), intent(in) :: force_const(:, :)

      real(q), allocatable :: force_const_blocks(:, :, :, :)
      integer :: array_size
      integer :: num_atoms
      integer :: atom_idx_1, atom_idx_2

      array_size = size(force_const, 1)
      num_atoms = array_size / 3
      if (num_atoms * 3 /= array_size .or. array_size /= size(force_const, 2)) &
         call vtutor%bug("IFC matrix has wrong size: " // &
         str(shape(force_const)) // " " // str(num_atoms), __FILE__, __LINE__)

      allocate(force_const_blocks(3, 3, num_atoms, num_atoms))
      do atom_idx_2 = 1, num_atoms
         do atom_idx_1 = 1, num_atoms
            force_const_blocks(:, :, atom_idx_1, atom_idx_2) = &
               force_const(atom_idx_1*3-2 : atom_idx_1*3, atom_idx_2*3-2 : atom_idx_2*3)
         enddo
      enddo

      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "force_constants", force_const_blocks))
   end subroutine

   subroutine write_polar_data(io, group_id, polar_data)
      use elphon_util, only: polar_data_t

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(polar_data_t), intent(in) :: polar_data

      real(q), allocatable :: buffer(:, :, :)
      integer :: num_prim, prim_idx

      if (allocated(polar_data%dielectric_tensor)) then
         VH5_CHECK(vh5_write(group_id, "dielectric_constant", transpose(polar_data%dielectric_tensor)))
      endif
      if (allocated(polar_data%born_eff_charges)) then
         num_prim = size(polar_data%born_eff_charges, 3)
         allocate(buffer(3, 3, num_prim))
         do prim_idx = 1, num_prim
            buffer(:, :, prim_idx) = transpose(polar_data%born_eff_charges(:, :, prim_idx))
         enddo
         VH5_CHECK(vh5_write(group_id, "born_effective_charges", buffer))
      endif
   end subroutine

   subroutine write_grid_dimensions(io, group_id, &
      lattice_points, grid_points, prim_fft_mesh, prim_to_super)

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      integer, intent(in) :: lattice_points(:, :)
      real(q), intent(in) :: grid_points(:, :)
      integer, intent(in) :: prim_fft_mesh(3)
      real(q), intent(in) :: prim_to_super(3, 3)

      real(q), allocatable :: grid_points_out(:, :)

      allocate(grid_points_out(3, size(grid_points, 2)))

      grid_points_out = matmul(prim_to_super, grid_points)

      VH5_CHECK(vh5_write(group_id, "lattice_point", lattice_points))
      VH5_CHECK(vh5_write(group_id, "grid_point", grid_points_out))
      VH5_CHECK(vh5_write(group_id, "FFT_mesh", prim_fft_mesh))
   end subroutine
#endif

   subroutine find_fft_dims(info, latt_cur, fft_mesh)
      use lattice, only: latt
      use constant, only: tpi, autoa, rytoev

      type(info_struct), intent(in) :: info
      type(latt), intent(in) :: latt_cur
      integer, intent(out) :: fft_mesh(3)

      real(q) :: cutoff(3)
      real(q) :: weight

      cutoff = sqrt(info%enmax / rytoev) / (tpi / (latt_cur%anorm(:) / autoa))

      if (info%szprec(1:7)=='singlen') then
         weight = 3
      else if (info%szprec(1:1)=='h' .or. info%szprec(1:1)=='a' .or. info%szprec(1:1)=='s' ) then
         weight = 4
      else
         weight = 3
      endif
      fft_mesh = cutoff * weight + 0.5_q
      call fftchk(fft_mesh)
   end subroutine

   subroutine find_target_grid(supercell_matrix, prim_fft_mesh, &
      target_fft_mesh, super_target_to_source, prim_source_to_target)
      use mkpoints, only: SmithNormalForm

      integer, intent(in) :: supercell_matrix(3, 3)
      integer, intent(in) :: prim_fft_mesh(3)
      integer, intent(out) :: target_fft_mesh(3)
      integer, intent(out) :: super_target_to_source(3, 3)
      integer, intent(out) :: prim_source_to_target(3, 3)

      integer :: i
      integer :: matrix(3, 3)
      integer :: normal_form(3, 3)

      if (is_diagonal_and_positive(supercell_matrix)) then
         super_target_to_source = 0
         prim_source_to_target = 0
         do i = 1, 3
            target_fft_mesh(i) = supercell_matrix(i, i) * prim_fft_mesh(i)
            super_target_to_source(i, i) = 1
            prim_source_to_target(i, i) = 1
         enddo
      else
         do i = 1, 3
            matrix(:, i) = prim_fft_mesh(:) * supercell_matrix(:, i)
         enddo
         call SmithNormalForm(matrix, prim_source_to_target, normal_form, super_target_to_source)
         do i = 1, 3
            target_fft_mesh(i) = normal_form(i, i)
         enddo
      endif

   contains
      pure logical function is_diagonal_and_positive(m)
         integer, intent(in) :: m(3, 3)
         integer :: a, b
         is_diagonal_and_positive = .false.
         do a = 1, 3
            do b = 1, 3
               if (a == b) then
                  if (m(a, b) < 0) return
               else
                  if (m(a, b) /= 0) return
               endif
            enddo
         enddo
         is_diagonal_and_positive = .true.
      end function
   end subroutine

   subroutine grid_interpolator_init(this, initial_grid, initial_lattice, target_to_initial, target_fft_mesh)
      use mgrid, only: inilgrd, gen_rc_grid, gen_rc_sub_grid

      class(grid_interpolator_t), intent(out) :: this
      type(grid_3d), intent(in) :: initial_grid
      real(q), intent(in) :: initial_lattice(3, 3)
      integer, intent(in) :: target_to_initial(3, 3)
      integer, intent(in) :: target_fft_mesh(3)

      integer, parameter :: identity(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

      this%source_grid = initial_grid
      this%recip_target_lattice = transpose(matrix_inverse_3x3_r( &
         matmul(initial_lattice, target_to_initial)))
      this%recip_source_lattice = transpose(matrix_inverse_3x3_r(initial_lattice))

      associate(tgrid => this%target_grid)
         ! Allocate storage for the FFT.
         allocate(this%buffer(initial_grid%mplwv))
         this%target_to_source = target_to_initial

         ! Set target FFT mesh.
         tgrid%ngptar = target_fft_mesh
         ! call FFTCHK(tgrid%ngptar)

         ! Initialize target grid.
         call inilgrd(tgrid%ngptar(1), tgrid%ngptar(2), tgrid%ngptar(3), tgrid)

         this%is_aligned = all(target_to_initial == identity)

         if (this%is_aligned) then
            ! Point to the same communicators as the initial grid.
            tgrid%comm => initial_grid%comm
            tgrid%comm_kin => initial_grid%comm_kin
            tgrid%comm_kinter => initial_grid%comm_kinter
            ! Create target grid as sub grid of the initial grid.
            ! This also sets the transition table between the two.
            call gen_rc_sub_grid(tgrid, initial_grid, this%aligned_transit, .true., .true.)
         else
            ! Create new communicator for target grid.
            call init_grid_comm(tgrid)
            ! Create target grid as a separate grid.
            call gen_rc_grid(tgrid)
            ! Create map between the two grids.
            call this%set_general_transit()
         endif

         ! Setup communication part of target grid.
         CALLMPI(mapset(tgrid))
      end associate

   contains

      subroutine init_grid_comm(grid)
#if defined(MPI) || defined(MPI_CHAIN)
         use mpi, only: MPI_COMM_SELF
         use mpimy, only: M_initc
#endif
         type(grid_3d), intent(inout) :: grid

         ! Initialize communicator of target grid to be fully serial.
         ! _test_mani this can be improved in the future
         allocate(grid%comm)
#if defined(MPI) || defined(MPI_CHAIN)
         grid%comm%mpi_comm = MPI_COMM_SELF
         CALL M_initc(grid%comm)
#endif
         ! Just point everything to the same serial communicator.
         ! _test_mani this can be improved in the future
         grid%comm_kin => grid%comm
         grid%comm_kinter => grid%comm
      end subroutine
   end subroutine

   subroutine grid_interpolator_interpolate(this, initial_values, interpol_values)
      class(grid_interpolator_t), intent(inout) :: this
      complex(q), intent(in) :: initial_values(:)
      complex(q), intent(out) :: interpol_values(:)

      integer :: pw_idx, mapped_pw_idx
      integer :: start_idx
      complex(q) :: coeff

      ! Copy initial values into FFT buffer
      this%buffer(:) = initial_values(:)
      ! FFT values into reciprocal space
      call FFT3D(this%buffer, this%source_grid, -1)
      ! Scale with number of points
      this%buffer = this%buffer / this%source_grid%nplwv

      ! Initialize plane-wave components to zero
      interpol_values = 0
      ! Copy plane-wave components into buffer of target FFT grid
      if (this%is_aligned) then
         ! Simple case, just copy plane-wave components
         call cp_grid(this%source_grid, this%target_grid, &
            this%aligned_transit, this%buffer, interpol_values)
      else
         start_idx = 1
         if (this%source_grid%real2cplx) start_idx = -this%source_grid%rc%np
         do pw_idx = start_idx, this%source_grid%rc%np
            mapped_pw_idx = this%general_transit(pw_idx)
            if (mapped_pw_idx /= 0) then
               coeff = this%buffer(abs(pw_idx))
               if (pw_idx < 0 .neqv. mapped_pw_idx < 0) coeff = conjg(coeff)
               interpol_values(abs(mapped_pw_idx)) = interpol_values(abs(mapped_pw_idx)) + coeff
            endif
         enddo

         do pw_idx = 1, this%target_grid%rc%np
            if (this%plane_wave_weight(pw_idx) > 1) then
               interpol_values(pw_idx) = interpol_values(pw_idx) / this%plane_wave_weight(pw_idx)
            endif
         enddo

         ! Accumulate all plane-wave components onto one rank
         ! _test_mani right now this is done for all ranks
         CALLMPI(M_sum_z(this%source_grid%comm, interpol_values(1), this%target_grid%mplwv))
      endif

      ! FFT values into real space
      call FFT3D(interpol_values, this%target_grid, 1)
   end subroutine

   subroutine grid_interpolator_set_general_transit(this)
      class(grid_interpolator_t), intent(inout) :: this

      integer, allocatable :: g_map(:, :, :)
      integer :: lower_bounds_1(3)
      integer :: lower_bounds_2(3)
      integer :: upper_bounds_2(3)
      integer :: column_idx
      integer :: x_idx, y_idx, z_idx
      integer :: g_vec(3)
      integer :: shift(3)
      integer :: pw_idx
      integer :: gz
      real(q) :: min_dist
      real(q), allocatable :: distances(:)
      type(min_image_mapper_t) :: min_image_mapper

      associate (gr1 => this%source_grid, gr2 => this%target_grid)
         ! Get min and max G-vector elements of target grid in each direction.
         lower_bounds_1 = [minval(gr1%lpctx), minval(gr1%lpcty), minval(gr1%lpctz)]
         lower_bounds_2 = [minval(gr2%lpctx), minval(gr2%lpcty), minval(gr2%lpctz)]
         upper_bounds_2 = [maxval(gr2%lpctx), maxval(gr2%lpcty), maxval(gr2%lpctz)]

         ! This array maps all G-vectors of the target grid to the corresponding
         ! plane-wave indices in the source grid. The indices of the map array
         ! are the G-vector components in direct coordinates.
         ! This mapping is done locally, i.e., no MPI communication takes place.
         ! _test_mani likely this has to change in the more general case.
         allocate(g_map( &
            lower_bounds_2(1) : upper_bounds_2(1), &
            lower_bounds_2(2) : upper_bounds_2(2), &
            lower_bounds_2(3) : upper_bounds_2(3)))
         ! G-vectors that do not exist in either grid are mapped to the value 0.
         g_map = 0

         allocate(distances(gr2%rc%np))
         distances = huge(1._q)

         ! We set this up purely to get minimum Cartesian distances for the G-vectors
         call min_image_mapper%init(this%recip_target_lattice * spread(gr2%ngptar, 1, 3))

         ! Loop over all G-vectors of the target grid
         pw_idx = 0
         do column_idx = 1, gr2%rc%ncol
            y_idx = gr2%rc%i2(column_idx)
            z_idx = gr2%rc%i3(column_idx)
            do x_idx = 1, gr2%rc%nrow
               ! Get plane-wave index of target grid
               pw_idx = pw_idx + 1
               ! Get G-vector of target grid in direct coordinates
               g_vec = [gr2%lpctx(x_idx), gr2%lpcty(y_idx), gr2%lpctz(z_idx)]
               ! Map G-vector to plane-wave index
               g_map(g_vec(1), g_vec(2), g_vec(3)) = pw_idx
               ! Set minimum distance for this G-vector
               shift = min_image_mapper%get_shift(real(g_vec, q) / gr2%ngptar, distances(pw_idx))

               ! Map -G explicitly if half-grid mode is used
               if (gr2%real2cplx .and. g_vec(3) > 0) then
                  g_vec = -g_vec
                  g_vec = modulo(g_vec - lower_bounds_2, gr2%ngptar) + lower_bounds_2
                  ! Sometimes we map back to positive z so we need to check
                  if (g_vec(3) < 0) g_map(g_vec(1), g_vec(2), g_vec(3)) = -pw_idx
               endif
            enddo
         enddo

         allocate(this%plane_wave_weight(gr2%rc%np))
         this%plane_wave_weight = 0

         if (gr1%real2cplx) then
            allocate(this%general_transit(-gr1%rc%np : gr1%rc%np))
         else
            allocate(this%general_transit(gr1%rc%np))
         endif
         this%general_transit = 0

         call min_image_mapper%init(this%recip_source_lattice * spread(gr1%ngptar, 1, 3))

         ! Loop over all G-vectors of the source grid
         pw_idx = 0
         do column_idx = 1, gr1%rc%ncol
            y_idx = gr1%rc%i2(column_idx)
            z_idx = gr1%rc%i3(column_idx)
            do x_idx = 1, gr1%rc%nrow
               ! Get plane-wave index of source grid
               pw_idx = pw_idx + 1
               ! Get G-vector of source grid in direct coordinates
               g_vec = [gr1%lpctx(x_idx), gr1%lpcty(y_idx), gr1%lpctz(z_idx)]
               ! Shift G-vector to have smallest distance from origin in Cartesian coordinates
               shift = min_image_mapper%get_shift(real(g_vec, q) / gr1%ngptar, min_dist)
               g_vec = g_vec + shift * gr1%ngptar
               ! Get corresponding G-vector in target grid
               g_vec = matmul(g_vec, this%target_to_source)
               ! Map plane-wave index of source grid to plane-wave index of target grid
               this%general_transit(pw_idx) = map_g_vec(g_vec)

               ! Map -G explicitly if half-grid mode is used
               if (gr1%real2cplx .and. gr1%lpctz(z_idx) > 0) then
                  gz = -gr1%lpctz(z_idx)
                  gz = modulo(gz - lower_bounds_1(3), gr1%ngz) + lower_bounds_1(3)
                  ! Sometimes we map back to positive z so we need to check
                  if (gz < 0) this%general_transit(-pw_idx) = map_g_vec(-g_vec)
               endif
            enddo
         enddo
      end associate

      CALLMPI(M_sum_i(this%source_grid%comm, this%plane_wave_weight, this%target_grid%rc%np))

   contains
      integer function map_g_vec(g_in) result(pw_idx)
         integer, intent(in) :: g_in(3)
         integer :: g(3)
         integer :: abs_idx

         g = modulo(g_in - lower_bounds_2, this%target_grid%ngptar) + lower_bounds_2
         pw_idx = g_map(g(1), g(2), g(3))

         abs_idx = abs(pw_idx)
         if (min_dist < distances(abs_idx) + 1e-7_q) then
            this%plane_wave_weight(abs_idx) = this%plane_wave_weight(abs_idx) + 1
         else
            pw_idx = 0
         endif
      end function
   end subroutine

   pure integer function det_3x3_i(a) result(det)
      integer, intent(in) :: a(3, 3)

      det   = a(1, 1) * a(2, 2) * a(3, 3) &
            + a(2, 1) * a(3, 2) * a(1, 3) &
            + a(1, 2) * a(2, 3) * a(3, 1) &
            - a(1, 3) * a(2, 2) * a(3, 1) &
            - a(2, 1) * a(1, 2) * a(3, 3) &
            - a(3, 2) * a(2, 3) * a(1, 1)
   end function

   pure function matrix_inverse_3x3_i(a) result(b)
      integer, intent(in) :: a(3, 3)
      integer :: b(3, 3)
      integer :: det_inv

      ! calculate inverse determinant.
      det_inv = 1 / det_3x3_i(a)

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
   end function

   pure real(q) function det_3x3_r(a) result(det)
      real(q), intent(in) :: a(3, 3)

      det   = a(1, 1) * a(2, 2) * a(3, 3) &
            + a(2, 1) * a(3, 2) * a(1, 3) &
            + a(1, 2) * a(2, 3) * a(3, 1) &
            - a(1, 3) * a(2, 2) * a(3, 1) &
            - a(2, 1) * a(1, 2) * a(3, 3) &
            - a(3, 2) * a(2, 3) * a(1, 1)
   end function

   pure function matrix_inverse_3x3_r(a) result(b)
      real(q), intent(in) :: a(3, 3)
      real(q) :: b(3, 3)
      real(q) :: det_inv

      ! calculate inverse determinant.
      det_inv = 1._q / det_3x3_r(a)

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
   end function

   subroutine generate_lattice_points(supercell_matrix, lattice_points)
      use mkpoints, only: HermiteNormalForm
      integer, intent(in) :: supercell_matrix(3, 3)
      integer, allocatable, intent(out) :: lattice_points(:, :)

      integer :: normal_form(3, 3), tmp(3, 3)
      integer :: supercell_size(3)
      integer :: point_idx
      integer :: x, y, z
      integer :: i

      call HermiteNormalForm(supercell_matrix, normal_form, tmp)
      do i = 1, 3
         supercell_size(i) = normal_form(i, i)
      enddo

      allocate(lattice_points(3, product(supercell_size)))

      point_idx = 0
      do x = 0, supercell_size(1) - 1
         do y = 0, supercell_size(2) - 1
            do z = 0, supercell_size(3) - 1
               point_idx = point_idx + 1
               lattice_points(:, point_idx) = [x, y, z]
            enddo
         enddo
      enddo
   end subroutine

   subroutine generate_grid_points(lattice_points, fft_mesh, grid_points)
      integer, intent(in) :: lattice_points(:, :)
      integer, intent(in) :: fft_mesh(3)
      real(q), allocatable, intent(out) :: grid_points(:, :)

      integer :: i1, i2, i3
      integer :: point_idx
      integer :: lattice_point_count
      integer :: lattice_point_idx

      lattice_point_count = size(lattice_points, 2)

      allocate(grid_points(3, lattice_point_count * product(fft_mesh)))

      point_idx = 0
      do lattice_point_idx = 1, lattice_point_count
         do i3 = 0, fft_mesh(3) - 1
            do i2 = 0, fft_mesh(2) - 1
               do i1 = 0, fft_mesh(1) - 1
                  point_idx = point_idx + 1
                  grid_points(:, point_idx) = &
                     lattice_points(:, lattice_point_idx) + &
                     [real(q):: i1, i2, i3] / fft_mesh
               enddo
            enddo
         enddo
      enddo
   end subroutine

   subroutine transfer_to_grid_points(grid, grid_points, &
      supercell_matrix, super_target_to_source, pot_in, pot_out)

      type(grid_3d), intent(in) :: grid
      real(q), intent(in) :: grid_points(:, :)
      integer, intent(in) :: supercell_matrix(3, 3)
      integer, intent(in) :: super_target_to_source(3, 3)
      complex(q), intent(in) :: pot_in(:)
      complex(q), intent(out) :: pot_out(:)

      integer :: num_out_points, point_idx
      integer :: fast_idx, slow_idx_1, slow_idx_2
      integer :: column_idx, pot_idx
      real(q) :: vec(3)
      integer :: idx_vec(3)
      real(q) :: trafo(3, 3)

      num_out_points = size(pot_out)

      trafo = matrix_inverse_3x3_r(real(matmul(supercell_matrix, super_target_to_source), q))

      fast_idx = grid%rl%nfast
      slow_idx_1 = modulo(fast_idx, 3) + 1
      slow_idx_2 = modulo(slow_idx_1, 3) + 1

      pot_out = 0

      do point_idx = 1, num_out_points
         vec = matmul(trafo, grid_points(:, point_idx))
         idx_vec = modulo(nint(vec * grid%ngptar), grid%ngptar)

         column_idx = grid%rl%index(idx_vec(slow_idx_1), idx_vec(slow_idx_2))
         if (column_idx == 0) cycle
         pot_idx = idx_vec(fast_idx) + (column_idx - 1) * grid%ngptar(fast_idx) + 1
         pot_out(point_idx) = get_val(pot_in, pot_idx)
      enddo

      CALLMPI(M_sum_z(grid%comm, pot_out(1), num_out_points))

   contains

      pure RGRID function get_val(arr, complex_idx)
         complex(q), intent(in) :: arr(:)
         integer, intent(in) :: complex_idx
#ifdef realmode
         integer :: array_idx

         array_idx = (complex_idx + 1) / 2
         if (mod(complex_idx, 2) == 1) then
            get_val = real(arr(array_idx), q)
         else
            get_val = aimag(arr(array_idx))
         endif
#else
         get_val = arr(complex_idx)
#endif
      end function
   end subroutine

   ! Scatter full grid from ionode to distributed grid
   ! Uses DIS_GRID_RL_PLANE to distribute each z-plane from ionode
   subroutine read_v_from_hdf5(grid, group_id, dataset_name, local_data)
#ifdef VASP_HDF5
      use vhdf5_base, only: hsize_t, vh5_read_double_subarray_nd

      type(grid_3d), intent(in) :: grid
      integer(hid_t), intent(in) :: group_id
      character(*), intent(in) :: dataset_name
      complex(q), intent(out) :: local_data(:, :)

      integer :: grid_read(3)
      integer :: nspin, ispin, nz
      real(q), allocatable :: plane(:, :)
      integer(hsize_t) :: start(4), amount(4)

      VH5_CHECK(vh5_read(group_id, dataset_name//"/grid", grid_read))
      if (any(grid_read /= [grid%ngx, grid%ngy, grid%ngz])) call vtutor%error( &
         "Could not read grid from HDF5: grid dimensions do not match")

      nspin = size(local_data, 2)
      allocate(plane(grid%ngx, grid%ngy))

      do ispin = 1, nspin
         do nz = 1, grid%ngz
            ! Read plane from HDF5 on ionode
            if (grid%comm%node_me == grid%comm%ionode) then
               start = [1, 1, nz, ispin]
               amount = [grid%ngx, grid%ngy, 1, 1]
               VH5_CHECK(vh5_read_double_subarray_nd(group_id, dataset_name//"/values", 4, start, amount, plane))
            endif

            ! Broadcast plane to all ranks
            call DIS_GRID_RL_PLANE(grid, plane, local_data(1, ispin), .true., nz)
         enddo
      enddo
#else
      type(grid_3d), intent(in) :: grid
      integer, intent(in) :: group_id
      character(*), intent(in) :: dataset_name
      complex(q), intent(out) :: local_data(:, :)
#endif
   end subroutine

   ! Write Dij to HDF5 using global sum within HDF5 group
   subroutine write_dij_to_hdf5(wdes, group_id, dataset_name, local_dij)
#ifdef VASP_HDF5
      use wave, only: ni_global
      use vhdf5_base, only: vh5_create_complex_array_nd, vh5_write_complex_subarray_nd, lwriteh5

      type(wavedes), intent(in) :: wdes
      integer(hid_t), intent(in) :: group_id
      character(*), intent(in) :: dataset_name
      OVERLAP, intent(in) :: local_dij(:, :, :, :)

      integer :: num_chan, ncdij, nions_global
      integer :: local_idx, global_idx
      integer :: writer_group
      OVERLAP, allocatable :: dij_full(:, :, :, :)

      num_chan = size(local_dij, 1)
      ncdij = size(local_dij, 4)

      ! Determine which group has the HDF5 writer
      writer_group = 0
      if (lwriteh5) writer_group = 1
      CALLMPI(M_max_i(wdes%comm_inb, writer_group, 1))
      
      ! Only the writer group participates
      if (writer_group == 0) return

      ! Compute global number of ions via reduction
      nions_global = wdes%nions
      CALLMPI(M_sum_i(wdes%comm_inb, nions_global, 1))

      ! Allocate full array (on all ranks in writer group)
      allocate(dij_full(num_chan, num_chan, nions_global, ncdij))
      dij_full = 0

      ! Each rank fills in its local ions
      do local_idx = 1, wdes%nions
         global_idx = ni_global(local_idx, wdes%comm_inb)
         dij_full(:, :, global_idx, :) = local_dij(:, :, local_idx, :)
      enddo

      ! Global sum across comm_inb (only within writer group)
#ifdef realmode
      CALLMPI(M_sum_d(wdes%comm_inb, dij_full, size(dij_full)))
#else
      CALLMPI(M_sum_z(wdes%comm_inb, dij_full, size(dij_full)))
#endif

      ! Write to HDF5 (only writer rank)
      VH5_CHECK(vh5_write(group_id, dataset_name, dij_full))
#else
      type(wavedes), intent(in) :: wdes
      integer, intent(in) :: group_id
      character(*), intent(in) :: dataset_name
      OVERLAP, intent(in) :: local_dij(:, :, :, :)
#endif
   end subroutine

   ! Read Dij from HDF5 ion-by-ion
   subroutine read_dij_from_hdf5(wdes, group_id, dataset_name, local_dij)
#ifdef VASP_HDF5
      use wave, only: ni_global
      use string, only: str
      use vhdf5_base, only: hsize_t, VH5_MAX_RANK, vh5_get_dimensions, &
         vh5_read_double_subarray_nd, vh5_read_double_complex_subarray_nd

      type(wavedes), intent(in) :: wdes
      integer(hid_t), intent(in) :: group_id
      character(*), intent(in) :: dataset_name
      OVERLAP, intent(out) :: local_dij(:, :, :, :)

      integer :: local_idx, global_idx
      integer :: dataset_rank
      logical :: dataset_is_complex
      integer :: dataset_dims(VH5_MAX_RANK)
      integer(hsize_t) :: start(4), amount(4)
#ifndef realmode
      real(q), allocatable :: real_buffer(:, :, :, :)
#endif

      dataset_is_complex = .false.

      VH5_CHECK(vh5_get_dimensions(group_id, dataset_name, dataset_dims, dataset_rank))
      if (dataset_rank == 4) then
         dataset_is_complex = .false.
      else if (dataset_rank == 5 .and. int(dataset_dims(1)) == 2) then
         dataset_is_complex = .true.
      else
         call vtutor%error( &
            "Dataset " // trim(dataset_name) // " has unsupported rank for Dij data: " // &
            str(dataset_rank))
      endif

#ifdef realmode
      if (dataset_is_complex) call vtutor%error( &
         "Cannot read complex-valued " // trim(dataset_name) // " in a real-only executable.")
#else
      if (.not. dataset_is_complex) then
         allocate(real_buffer(size(local_dij, 1), size(local_dij, 2), 1, size(local_dij, 4)))
      endif
#endif

      ! Process each local ion on this rank (no communication needed)
      do local_idx = 1, wdes%nions
         global_idx = ni_global(local_idx, wdes%comm_inb)
         
         ! Read this ion's data directly from HDF5 (each rank reads independently)
         start = [1, 1, global_idx, 1]
         amount = [size(local_dij, 1), size(local_dij, 2), 1, size(local_dij, 4)]
#ifdef realmode
         VH5_CHECK(vh5_read_double_subarray_nd(group_id, dataset_name, 4, start, amount, local_dij(:, :, local_idx, :)))
#else
         if (dataset_is_complex) then
            VH5_CHECK(vh5_read_double_complex_subarray_nd(group_id, dataset_name, 4, start, amount, local_dij(:, :, local_idx, :)))
         else
            VH5_CHECK(vh5_read_double_subarray_nd(group_id, dataset_name, 4, start, amount, real_buffer))
            local_dij(:, :, local_idx, :) = cmplx(real_buffer(:, :, 1, :), 0._q, kind=q)
         endif
#endif
      enddo
#else
      type(wavedes), intent(in) :: wdes
      integer, intent(in) :: group_id
      character(*), intent(in) :: dataset_name
      OVERLAP, intent(out) :: local_dij(:, :, :, :)
#endif
   end subroutine

   ! _test_mani this should be removed when the new finite-difference code is here
   subroutine calc_elphon_derivatives(progress,  &
      v, dv_du, dij, ddij_du, mu, dmu_du)
      use findiff_struct_def, only: findiff_progress_t

      type(findiff_progress_t), intent(in) :: progress

      complex(q), intent(in) :: v(:, :)
      complex(q), intent(inout) :: dv_du(:, :, :, :)
      OVERLAP, intent(in) :: dij(:, :, :, :)
      OVERLAP, intent(inout) :: ddij_du(:, :, :, :, :, :)
      complex(q), optional, intent(in) :: mu(:, :)
      complex(q), optional, intent(inout) :: dmu_du(:, :, :, :)

      if (present(mu) .neqv. present(dmu_du)) call vtutor%bug( &
         "mu and dmu_du must be both present or both absent", __FILE__, __LINE__)

      call calc_deriv_v(progress%dof_index, progress%displ_index, progress%step, v, dv_du)
      call calc_deriv_dij(progress%dof_index, progress%displ_index, progress%step, dij, ddij_du)
      if (present(mu)) then
         call calc_deriv_v(progress%dof_index, progress%displ_index, progress%step, mu, dmu_du)
      endif
   end subroutine

   subroutine write_elphon_checkpoint(checkpoint, gridc, wdes, v, dij, mu)
      use checkpointing, only: checkpoint_struct

      type(checkpoint_struct), intent(inout) :: checkpoint
      type(grid_3d), intent(in) :: gridc
      type(wavedes), intent(in) :: wdes
      complex(q), intent(in) :: v(:, :)
      OVERLAP, intent(in) :: dij(:, :, :, :)
      complex(q), optional, intent(in) :: mu(:, :)

      call checkpoint%write_data("elphon", write_data)

   contains
      subroutine write_data(group_id)
#ifdef VASP_HDF5
         use vhdf5, only: vh5_write_grid_quantity
         integer(hid_t), intent(in) :: group_id

         call vh5_write_grid_quantity(group_id, 'V', 'values', gridc, v)
         if (present(mu)) call vh5_write_grid_quantity(group_id, 'mu', 'values', gridc, mu)
         call write_dij_to_hdf5(wdes, group_id, 'Dij', dij)
#else
         integer, intent(in) :: group_id
#endif
      end subroutine
   end subroutine

   subroutine read_elphon_checkpoint(checkpoint, displ_count, step, gridc, wdes, dv_du, ddij_du, dmu_du)
      use checkpointing, only: checkpoint_struct

      type(checkpoint_struct), intent(in) :: checkpoint
      integer, intent(in) :: displ_count
      real(q), intent(in) :: step
      type(grid_3d), intent(in) :: gridc
      type(wavedes), intent(in) :: wdes
      complex(q), intent(inout) :: dv_du(:, :, :, :)
      OVERLAP, intent(inout) :: ddij_du(:, :, :, :, :, :)
      complex(q), optional, intent(inout) :: dmu_du(:, :, :, :)

      call checkpoint%read_data("elphon", read_data)

   contains
      subroutine read_data(group_id, idx)
#ifdef VASP_HDF5
         integer(hid_t), intent(in) :: group_id
         integer, intent(in) :: idx

         integer :: dof_index, displ_index
         complex(q), allocatable :: v_read(:, :)
         OVERLAP, allocatable :: dij_read(:, :, :, :)

         ! Infer DOF and displacement indices from idx using progress structure
         dof_index = (idx - 1) / displ_count + 1
         displ_index = modulo(idx - 1, displ_count) + 1

         ! Allocate temporary arrays for reading
         allocate(v_read(size(dv_du, 1), size(dv_du, 2)))
         allocate(dij_read(size(ddij_du, 1), size(ddij_du, 2), size(ddij_du, 3), size(ddij_du, 4)))

         ! Read and scatter V from full grid on ionode to distributed grid
         call read_v_from_hdf5(gridc, group_id, 'V', v_read)

         ! Read global Dij and scatter to distributed array ion-by-ion (NCORE > 1)
         call read_dij_from_hdf5(wdes, group_id, 'Dij', dij_read)

         ! Apply finite difference formulas using existing helper routines
         ! Use progress%step as it's stored in the progress structure
         call calc_deriv_v(dof_index, displ_index, step, v_read, dv_du)
         call calc_deriv_dij(dof_index, displ_index, step, dij_read, ddij_du)

         ! Handle mu if present - reuse v_read memory
         if (present(dmu_du)) then
            deallocate(v_read)
            allocate(v_read(size(dmu_du, 1), size(dmu_du, 2)))
            call read_v_from_hdf5(gridc, group_id, 'mu', v_read)
            call calc_deriv_v(dof_index, displ_index, step, v_read, dmu_du)
         endif
#else
         integer, intent(in) :: group_id
         integer, intent(in) :: idx
#endif
      end subroutine
   end subroutine

   subroutine calc_deriv_v(dof_idx, displ_idx, step, quant, deriv)
      integer, intent(in) :: dof_idx
      integer, intent(in) :: displ_idx
      real(q), intent(in) :: step
      complex(q), intent(in) :: quant(:, :)
      complex(q), intent(inout) :: deriv(:, :, :, :)

      integer :: prim_idx, cart_idx

      prim_idx = (dof_idx + 2) / 3
      cart_idx = modulo(dof_idx - 1, 3) + 1

      if (displ_idx == 1) then
         deriv(:, :, cart_idx, prim_idx) = quant(:, :)
      else
         deriv(:, :, cart_idx, prim_idx) = &
            (deriv(:, :, cart_idx, prim_idx) - quant(:, :)) / &
            (2 * step)
      endif
   end subroutine

   subroutine calc_deriv_dij(dof_idx, displ_idx, step, quant, deriv)
      integer, intent(in) :: dof_idx
      integer, intent(in) :: displ_idx
      real(q), intent(in) :: step
      OVERLAP, intent(in) :: quant(:, :, :, :)
      OVERLAP, intent(inout) :: deriv(:, :, :, :, :, :)

      integer :: prim_idx, cart_idx

      prim_idx = (dof_idx + 2) / 3
      cart_idx = modulo(dof_idx - 1, 3) + 1

      if (displ_idx == 1) then
         deriv(:, :, :, :, cart_idx, prim_idx) = quant(:, :, :, :)
      else
         deriv(:, :, :, :, cart_idx, prim_idx) = &
            (deriv(:, :, :, :, cart_idx, prim_idx) - quant(:, :, :, :)) / &
            (2 * step)
      endif
   end subroutine

end module
