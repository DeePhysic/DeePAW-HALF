#include "symbol.inc"

#ifdef libvaspml

module mvaspml

   use iso_c_binding, only: c_associated, &
                            c_char,       &
                            c_double,     &
                            c_int,        &
                            c_int8_t,     &
                            c_null_ptr,   &
                            c_ptr,        &
                            C_NULL_CHAR,  &
                            C_NEW_LINE
   use ml_ff_prec
   use tutor, only: vtutor
   use mpi_help, only : ml_mpi_par
   use verlet_cell_list_mod, only : verlet_cell_list_type

   implicit none

   private
   public :: vaspml_wrapper

   !================================================================================================
   !> Wrapper around VASPml library.
   !>
   !> This derived type stores a C pointer which works as a handle to a VASPml object. This object
   !> is responsible for carrying out machine-learned force field related work, e.g. computing
   !> predictions with a given force field. The object is created when the reader() member function
   !> is called and destroyed with finalize(). Member functions are used to pass data between VASP
   !> and VASPml.
   !================================================================================================
   type vaspml_wrapper

      private
      !> C-pointer to instance of VASPml object.
      type(c_ptr)                   :: handle = c_null_ptr
      type(verlet_cell_list_type)   :: nearest_neighbor_radial
      type(verlet_cell_list_type)   :: nearest_neighbor_angular
      real(kind=q),allocatable      :: posion(:, :)
      real(kind=q)                  :: lattice(1:3, 1:3) = 0.0_q
      real(kind=q)                  :: inverse_lattice(1:3, 1:3) = 0.0_q
      integer                       :: nions = 0
      integer                       :: ntype = 0
      integer, allocatable          :: itype(:)
      type(ml_mpi_par)              :: communicator
      real(kind=q)                  :: w1 = 0.0_q
      real(kind=q)                  :: w2 = 0.0_q
      real(kind=q)                  :: rcut1 = 0.0_q
      real(kind=q)                  :: rcut2 = 0.0_q
      integer                       :: number_atoms_per_node = 0
      integer, allocatable          :: number_atoms_per_type_per_node(:)
      character(len=1), allocatable :: types(:)
      real(kind=q)                  :: omega = 0.0_q

      contains

      procedure, public, pass(this) :: reader => vaspml_reader
      procedure, public, pass(this) :: init => vaspml_init
      procedure, public, pass(this) :: post_sync_vcaimages => vaspml_post_sync_vcaimages
      procedure, public, pass(this) :: compute => vaspml_compute
      procedure, public, pass(this) :: finalize => vaspml_finalize
      procedure, pass(this) :: fill_neighbor_list_arrays => vaspml_fill_neighbor_list_arrays
      procedure, pass(this) :: rescale_atoms => vaspml_rescale_atoms
      procedure, pass(this) :: define_elements_per_thread => vaspml_define_elements_per_thread
      procedure, pass(this) :: fill_force_array => vaspml_fill_force_array
      procedure, pass(this) :: get_energy => vaspml_get_energy
      procedure, pass(this) :: get_stress_tensor => vaspml_get_stress_tensor

   end type vaspml_wrapper

   interface

      !=============================================================================================
      !> Interface to C function which creates a VASPml object handle.
      !=============================================================================================
      function create_handle() result(ptr) bind(c, name="vaspml_createHandle")
         use iso_c_binding, only: c_ptr
         type(c_ptr) :: ptr
      end function create_handle

      !=============================================================================================
      !> Interface to C function which prepares the setup from basic run information.
      !=============================================================================================
      subroutine prepare_setup(ptr,        &
                               ntyp,       &
                               nsw,        &
                               ibrion,     &
                               lmlabexist, &
                               lmlffexist) bind(c, name="vaspml_prepareSetup")
         use iso_c_binding, only: c_ptr, c_int, c_int8_t
         type(c_ptr), intent(inout)         :: ptr
         integer(kind=c_int), intent(in)    :: ntyp
         integer(kind=c_int), intent(in)    :: nsw
         integer(kind=c_int), intent(in)    :: ibrion
         integer(kind=c_int8_t), intent(in) :: lmlabexist
         integer(kind=c_int8_t), intent(in) :: lmlffexist
      end subroutine prepare_setup

      !=============================================================================================
      !> Interface to C function which transfers and processes the INCAR file.
      !>
      !> @note This does not yet trigger the full INCAR setup but only the first part
      !> (io::readIncar()). The setup part (io::setupFromIncar()) will be called in the init()
      !> function when the atom type names are known.
      !=============================================================================================
      subroutine read_incar(ptr, incar) bind(c, name="vaspml_readIncar")
         use iso_c_binding, only: c_ptr, c_char
         type(c_ptr), intent(inout)                       :: ptr
         character(kind=c_char), dimension(*), intent(in) :: incar
      end subroutine read_incar

      !=============================================================================================
      !> Interface to C function which finalizes setup (io::setupFromIncar()).
      !=============================================================================================
      subroutine setup_from_incar(ptr, types) bind(c, name="vaspml_setupFromIncar")
         use iso_c_binding, only: c_ptr, c_char
         type(c_ptr), intent(inout)                       :: ptr
         character(kind=c_char), dimension(*), intent(in) :: types
      end subroutine setup_from_incar

      !=============================================================================================
      !> Interface to C function which updates setup after VCAIMAGES were synchronized.
      !=============================================================================================
      subroutine update_setup(ptr, outblockSync) bind(c, name="vaspml_updateSetup")
         use iso_c_binding, only: c_ptr, c_int
         type(c_ptr), intent(inout)      :: ptr
         integer(kind=c_int), intent(in) :: outblockSync
      end subroutine update_setup

      !=============================================================================================
      !> Interface to C function which passes communicator.
      !=============================================================================================
      subroutine setup_mpi(ptr, mpi_comm) bind(c, name="vaspml_setupMpi")
         use iso_c_binding, only: c_ptr, c_int
         type(c_ptr), intent(inout) :: ptr
         integer(c_int), intent(in) :: mpi_comm
      end subroutine setup_mpi

      !=============================================================================================
      !> Interface to C function which reads in the force field data.
      !=============================================================================================
      subroutine setup_force_field(ptr) bind(c, name="vaspml_setupForceField")
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(inout) :: ptr
      end subroutine setup_force_field

      !=============================================================================================
      !> Interface to C function which returns ML_W1.
      !=============================================================================================
      function get_w1(ptr) bind(c, name="vaspml_get_W1") result(w1)
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout) :: ptr
         real(kind=c_double)        :: w1
      end function get_w1

      !=============================================================================================
      !> Interface to C function which returns ML_W2.
      !=============================================================================================
      function get_w2(ptr) bind(c, name="vaspml_get_W2") result(w2)
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout) :: ptr
         real(kind=c_double)        :: w2
      end function get_w2

      !=============================================================================================
      !> Interface to C function which returns ML_RCUT1.
      !=============================================================================================
      function get_rcut1(ptr) bind(c, name="vaspml_get_RCUT1") result(rcut1)
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout) :: ptr
         real(kind=c_double)        :: rcut1
      end function get_rcut1

      !=============================================================================================
      !> Interface to C function which returns ML_RCUT2.
      !=============================================================================================
      function get_rcut2(ptr) bind(c, name="vaspml_get_RCUT2") result(rcut2)
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout) :: ptr
         real(kind=c_double)        :: rcut2
      end function get_rcut2

      subroutine reset_neighbor_arrays(ptr, key) bind(c, name="vaspml_resetNeighborArrays")
         use iso_c_binding, only: c_ptr, c_char
         type(c_ptr), intent(inout)                       :: ptr
         character(kind=c_char), dimension(*), intent(in) :: key
      end subroutine reset_neighbor_arrays

      !=============================================================================================
      !> Interface to C function which assigns number of atoms per core to c++ neighbor list.
      !=============================================================================================
      subroutine resize_neighbor_arrays(ptr, nions, key) bind(c, name="vaspml_resizeNeighborArrays")
         use iso_c_binding, only: c_ptr, c_int, c_char
         type(c_ptr), intent(inout)                       :: ptr
         integer(kind=c_int)                              :: nions
         character(kind=c_char), dimension(*), intent(in) :: key
      end subroutine resize_neighbor_arrays

      !=============================================================================================
      !> Interface to C function which destroys VASPml object.
      !=============================================================================================
      subroutine destroy_handle(ptr) bind(c, name="vaspml_destroyHandle")
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(inout) :: ptr
      end subroutine destroy_handle

      !=============================================================================================
      !> Interface to C function which sets the number of atoms per type per node for the
      !> c++ neighbor list
      !>
      !> ntypes number of types determines length of number_atoms_per_type_per_node
      !> number_atoms_per_type_per_node stores the number of elements for certain type for
      !> a single thread
      !=============================================================================================
      subroutine set_types_per_atom(ptr,                            &
                                    ntypes,                         &
                                    number_atoms_per_type_per_node, &
                                    key) bind(c, name="vaspml_set_nAtomsType")
         use iso_c_binding, only: c_ptr, c_int, c_char
         type(c_ptr), intent(inout)      :: ptr
         integer(kind=c_int), intent(in) :: ntypes
         integer(kind=c_int), intent(in) :: number_atoms_per_type_per_node(1:ntypes)
         character(kind=c_char), dimension(*), intent(in) :: key
      end subroutine set_types_per_atom

      subroutine fill_neighbor_list_single_atom(ptr,              &
                                                number_neighbors, &
                                                number_element,   &
                                                central_type,     &
                                                neighbors,        &
                                                neighbors_types,  &
                                                neighbor_dist,    &
                                                neighbor_connect, &
                                                key)              &
                                                bind(c, name="vaspml_fillNeighborArraysSingleAtom")
         use iso_c_binding, only: c_ptr, c_int, c_double, c_char
         type(c_ptr), intent(inout)                       :: ptr
         integer(kind=c_int),intent(in)                   :: number_neighbors
         integer(kind=c_int),intent(in)                   :: number_element
         integer(kind=c_int),intent(in)                   :: central_type
         integer(kind=c_int),intent(in)                   :: neighbors(1:number_neighbors)
         integer(kind=c_int),intent(in)                   :: neighbors_types(1:number_neighbors)
         real(kind=c_double),intent(in)                   :: neighbor_dist(1:number_neighbors)
         real(kind=c_double),intent(in)                   :: neighbor_connect(1:3, number_neighbors)
         character(kind=c_char), dimension(*), intent(in) :: key
      end subroutine fill_neighbor_list_single_atom

      !=============================================================================================
      !> Update the VASPml object such that forces are computed. These
      !> forces can then be sent to vasp for ionic update
      !=============================================================================================
      subroutine update_vaspml(ptr, omega) bind(c, name="vaspml_update")
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout)      :: ptr
         real(kind=c_double), intent(in) :: omega
      end subroutine update_vaspml

      !=============================================================================================
      !> set up the type map which converts structure types to force field types in the vaspml
      !> types atom types by name as defined in the POSCAR file
      !=============================================================================================
      subroutine set_type_map(ptr, types) bind(c, name="vaspml_set_typeMap")
         use iso_c_binding, only: c_ptr, c_char
         type(c_ptr), intent(inout)                       :: ptr
         character(kind=c_char), dimension(*), intent(in) :: types
      end subroutine set_type_map

      !=============================================================================================
      !> compute total force acting on a single atom
      !> ptr c pionter to vaspml interface
      !> nions total number of ions (not number of ions per core or mpi rank)
      !> ion_index ion index of atom for which force is computed
      !> local_index Ion index in the c++ part of the code
      !> tifor array in which forces are stored.
      !=============================================================================================
      subroutine fill_force_single_atom(ptr,         &
                                        nions,       &
                                        ion_index,   &
                                        local_index, &
                                        tifor) bind(c, name="vaspml_fillForceSingleAtom")
         use iso_c_binding, only: c_ptr, c_int, c_double
         type(c_ptr), intent(inout)        :: ptr
         integer(kind=c_int), intent(in)   :: nions
         integer(kind=c_int), intent(in)   :: ion_index
         integer(kind=c_int), intent(in)   :: local_index
         real(kind=c_double),intent(inout) :: tifor(1:3, 1:nions)
      end subroutine fill_force_single_atom

      subroutine compute_forces(ptr, nions, tifor) bind(c, name="vaspml_computeForces")
         use iso_c_binding, only: c_ptr, c_int, c_double
         type(c_ptr), intent(inout)        :: ptr
         integer(kind=c_int), intent(in)   :: nions
         real(kind=c_double),intent(inout) :: tifor(1:3, 1:nions)
      end subroutine compute_forces

      !=============================================================================================
      !> return potential energy of system from c++ code
      !=============================================================================================
      subroutine get_potential_energy(ptr, toten) bind(c, name="vaspml_getPotentialEnergy")
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout)        :: ptr
         real(kind=c_double),intent(inout) :: toten
      end subroutine get_potential_energy

      !=============================================================================================
      !> return stress tensor of system from c++ code
      !=============================================================================================
      subroutine get_stress_tensor_vaspml(ptr, tsif) bind(c, name="vaspml_getStressTensor")
         use iso_c_binding, only: c_ptr, c_double
         type(c_ptr), intent(inout)        :: ptr
         real(kind=c_double),intent(inout) :: tsif(1:3, 1:3)
      end subroutine get_stress_tensor_vaspml

      subroutine make_global_ion_index(ptr, nions) bind(c, name="vaspml_makeGlobalIonIndex")
         use iso_c_binding, only: c_ptr, c_int
         type(c_ptr), intent(inout)      :: ptr
         integer(kind=c_int), intent(in) :: nions
      end subroutine make_global_ion_index

      subroutine copy_neighbor_arrays(ptr) bind(c, name="vaspml_copyNeighborArrays")
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(inout) :: ptr
      end subroutine

      subroutine copy_force_array_to_cpu(ptr) bind(c, name="vaspml_copyForceArrayToCPU")
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(inout) :: ptr
      end subroutine

      subroutine copy_force_array_to_gpu(ptr) bind(c, name="vaspml_copyForceArrayToGPU")
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(inout) :: ptr
      end subroutine

   end interface

   contains

   !================================================================================================
   !> Create VASPml object and transfer INCAR data (key/value - pairs).
   !================================================================================================
   subroutine vaspml_reader(this, incar, ntyp, nsw, ibrion, lmlabexist, lmlffexist)

      use incar_reader, only: incar_file, str
      class(vaspml_wrapper), intent(inout) :: this
      type(incar_file), intent(inout)      :: incar
      integer, intent(in)                  :: ntyp
      integer, intent(in)                  :: nsw
      integer, intent(in)                  :: ibrion
      logical, intent(in)                  :: lmlabexist
      logical, intent(in)                  :: lmlffexist

      integer :: i
      integer :: n
      integer(kind=c_int) :: c_ntyp
      integer(kind=c_int) :: c_nsw
      integer(kind=c_int) :: c_ibrion
      integer(kind=c_int8_t) :: c_lmlabexist
      integer(kind=c_int8_t) :: c_lmlffexist
      character(len=:), allocatable :: incar_string
      character(kind=c_char), allocatable :: incar_to_c(:)

      if (c_associated(this%handle)) then
         call vtutor%bug("vaspml_wrapper: Trying to create another VASPml object &
                         &although pointer is already associated.", &
                         __FILE__, __LINE__)
      end if
      this%handle = create_handle()
      if (.not. c_associated(this%handle)) then
         call vtutor%bug("vaspml_wrapper: Could not create VASPml object.", &
                         __FILE__, __LINE__)
      end if

      ! Combine all INCAR tags to a single string using the overloaded str() method from the
      ! incar_reader module. Then convert deferred-length character string to one-dimensional
      ! character array which can be handed over to C.
      incar_string = str(incar)
      incar_to_c = f2c_dlcstring(incar_string)
      c_ntyp = ntyp
      c_nsw = nsw
      c_ibrion = ibrion
      c_lmlabexist = f2c_logical(lmlabexist)
      c_lmlffexist = f2c_logical(lmlffexist)

      call prepare_setup(this%handle, c_ntyp, c_nsw, c_ibrion, c_lmlabexist, c_lmlffexist)
      call read_incar(this%handle, incar_to_c)

   end subroutine vaspml_reader

   !================================================================================================
   !> Initialize VASPml object and set up force field.
   !>
   !> communicator mpi world communicator from vasp
   !================================================================================================
   subroutine vaspml_init(this,            &
                          communicator,    &
                          posion,          &
                          lattice,         &
                          inverse_lattice, &
                          nions,           &
                          ntype,           &
                          itype,           &
                          types)

      use ml_ff_constant, only : autoa, pi
      use ml_main_subroutines, only : initialize_cell_list_algorithm
      class(vaspml_wrapper), intent(inout) :: this
      type(ml_mpi_par), intent(in)         :: communicator
      real(kind=q), intent(in)             :: posion(:, :)
      real(kind=q), intent(in)             :: lattice(:, :)
      real(kind=q), intent(in)             :: inverse_lattice(:, :)
      integer, intent(in)                  :: nions
      integer, intent(in)                  :: ntype
      integer, intent(in)                  :: itype(:)
      character(len=2),intent(in)          :: types(:)

      this%communicator         = communicator
      this%communicator%node_me = this%communicator%node_me - 1
      this%nions                = nions
      if (nions < this%communicator%ncpu) then
         call vtutor%error("Cannot initialize VASPml because there are more M" // "PI ranks than &
                           &atoms. Please reduce the number of M" // "PI processes.")
      end if
      this%ntype           = ntype
      this%lattice         = lattice / autoa
      this%inverse_lattice = inverse_lattice * autoa * 2.0_q * pi

      if (allocated(this%posion)) deallocate(this%posion)
      allocate(this%posion(1:3, nions))
      if (allocated(this%itype)) deallocate(this%itype)
      allocate(this%itype(1:nions))
      this%itype = itype
      this%types = f2c_character1d(types, ";")

      call setup_from_incar(this%handle, this%types)
      call setup_mpi(this%handle,this%communicator%mpi_comm)
      call setup_force_field(this%handle)

      ! rescale positions to cartesian coordinates
      call this%rescale_atoms(posion, this%lattice)
      this%w1    =  get_w1(this%handle)
      this%w2    =  get_w2(this%handle)
      this%rcut1 =  get_rcut1(this%handle)
      this%rcut2 =  get_rcut2(this%handle)
      call initialize_cell_list_algorithm(this%nearest_neighbor_radial,  &
                                          this%nearest_neighbor_angular, &
                                          this%rcut1,                    &
                                          this%w1,                       &
                                          this%rcut2,                    &
                                          this%w2,                       &
                                          this%nions,                    &
                                          this%lattice,                  &
                                          this%inverse_lattice,          &
                                          this%itype,                    &
                                          this%posion,                   &
                                          this%ntype,                    &
                                          this%communicator%ncpu,        &
                                          this%communicator%node_me)
      call make_global_ion_index(this%handle, nions)

   end subroutine vaspml_init

   !================================================================================================
   !> Forward potentially updated tag values (VCAIMAGES) to VASPml setup.
   !================================================================================================
   subroutine vaspml_post_sync_vcaimages(this, outblock)

      class(vaspml_wrapper), intent(inout) :: this
      integer, intent(in)                  :: outblock

      integer(kind=c_int) :: c_outblock

      c_outblock = outblock
      call update_setup(this%handle, c_outblock)

   end subroutine vaspml_post_sync_vcaimages

   !================================================================================================
   !> Compute and transfer neighbor lists, execute force field prediction, transfer back data.
   !================================================================================================
   subroutine vaspml_compute(this,            &
                             posion,          &
                             lattice,         &
                             inverse_lattice, &
                             nions,           &
                             ntype,           &
                             itype,           &
                             types,           &
                             omega,           &
                             toten,           &
                             tifor,           &
                             tsif)

      use ml_ff_constant, only : autoa, pi
      use force_field, only : compute_nearest_neighbor_arrays
      class(vaspml_wrapper), intent(inout) :: this
      real(kind=q), intent(in)             :: posion(:, :)
      real(kind=q), intent(in)             :: lattice(:, :)
      real(kind=q), intent(in)             :: inverse_lattice(:, :)
      integer, intent(in)                  :: nions
      integer, intent(in)                  :: ntype
      integer, intent(in)                  :: itype(:)
      character(len=2),intent(in)          :: types(:)
      real(kind=q), intent(in)             :: omega
      real(kind=q), intent(out)            :: toten
      real(kind=q), intent(out)            :: tifor(:, :)
      real(kind=q), intent(out)            :: tsif(:, :)

      integer :: mneib_est1_dummy
      integer :: mneib_est2_dummy

      PROFILING_START('VASPML-UPDATE')
      PROFILING_START('PREPARE_VASPML')
      this%nions = nions
      this%ntype = ntype
      this%lattice = lattice / autoa
      this%inverse_lattice = inverse_lattice * autoa * 2.0_q * pi
      this%omega = omega / autoa**3
      if (allocated(this%posion)) deallocate(this%posion)
      allocate(this%posion(1:3, nions))
      if (allocated(this%itype)) deallocate(this%itype)
      allocate(this%itype(1:nions))
      this%nions = nions
      this%itype = itype
      this%types = f2c_character1d(types, ";")
      call set_type_map( this%handle, this%types )
      ! rescale positions to cartesian coordinates
      call this%rescale_atoms(posion, this%lattice)
      PROFILING_STOP('PREPARE_VASPML')

      PROFILING_START('compute_nearest_neighbor_arrays')
      call compute_nearest_neighbor_arrays(this%nearest_neighbor_radial,  &
                                           this%nearest_neighbor_angular, &
                                           this%lattice,                  &
                                           this%inverse_lattice,          &
                                           this%itype,                    &
                                           this%ntype,                    &
                                           this%nions,                    &
                                           this%nions,                    &
                                           this%posion,                   &
                                           this%communicator,             &
                                           ! radial variables
                                           mneib_est1_dummy,              &
                                           this%w1,                       &
                                           this%rcut1,                    &
                                           ! angular variables
                                           mneib_est2_dummy,              &
                                           this%w2,                       &
                                           this%rcut2)
      PROFILING_STOP('compute_nearest_neighbor_arrays')
      PROFILING_START('fill_neighbor_list_arrays')
      !! fill neighbor arrays
      call this%fill_neighbor_list_arrays()
      PROFILING_STOP('fill_neighbor_list_arrays')
      !!! call interface update
      PROFILING_START('update_vaspml')
      call update_vaspml(this%handle, this%omega)
      PROFILING_STOP('update_vaspml')
      !!! fill toten
      PROFILING_START('get_energy')
      call this%get_energy(toten)
      PROFILING_STOP('get_energy')
      !!! fill tifor
      PROFILING_START('fill_force_array')
      call this%fill_force_array(tifor)
      PROFILING_STOP('fill_force_array')
      !!! fill tisif
      PROFILING_START('get_stress_tensor')
      call this%get_stress_tensor(tsif)
      PROFILING_STOP('get_stress_tensor')
      PROFILING_STOP('VASPML-UPDATE')

   end subroutine vaspml_compute

   !================================================================================================
   !> Destroy VASPml object.
   !================================================================================================
   subroutine vaspml_finalize(this)

      class(vaspml_wrapper), intent(inout) :: this

      if (.not. c_associated(this%handle)) then
         call vtutor%bug("vaspml_wrapper: Trying to destroy VASPml object &
                         &although pointer is not yet associated.", &
                         __FILE__, __LINE__)
      end if
      call destroy_handle(this%handle)
      this%handle = c_null_ptr

   end subroutine vaspml_finalize


   !===============================!
   ! Private type-bound procedures !
   !===============================!


   subroutine vaspml_fill_neighbor_list_arrays(this)

      class(vaspml_wrapper), intent(inout) :: this

      integer :: inions
      integer :: inions_tmp
      integer :: number_element

      call this%define_elements_per_thread()
      if (this%w1 > 0) then
         !call reset_neighbor_arrays(this%handle, "2-body"//C_NULL_CHAR)
         call resize_neighbor_arrays(this%handle,                &
                                     this%number_atoms_per_node, &
                                     "2-body"//C_NULL_CHAR)
         call set_types_per_atom(this%handle,                         &
                                 this%ntype,                          &
                                 this%number_atoms_per_type_per_node, &
                                 "2-body"//C_NULL_CHAR)
      end if
      if (this%w2 > 0) then
         !call reset_neighbor_arrays(this%handle, "3-body"//C_NULL_CHAR)
         call resize_neighbor_arrays(this%handle,                &
                                     this%number_atoms_per_node, &
                                     "3-body"//C_NULL_CHAR)
         call set_types_per_atom(this%handle,                         &
                                 this%ntype,                          &
                                 this%number_atoms_per_type_per_node, &
                                 "3-body"//C_NULL_CHAR)
      end if

      number_element = 0
      do inions_tmp = 1, this%nions, this%communicator%ncpu
         inions = inions_tmp + this%communicator%node_me
         if (inions <= this%nions) then
            ! fill 2-body neighbor list
            if (this%w1 > 0 ) then
               ! Mind: counting starts at zero in C++ -> minus 1
               call fill_neighbor_list_single_atom(                                 &
                  this%handle,                                                      &
                  size(this%nearest_neighbor_radial%nn_list(inions)%neighbors),     &
                  number_element,                                                   &
                  this%itype(inions) - 1,                                           &
                  this%nearest_neighbor_radial%nn_list(inions)%neighbors - 1,       &
                  this%nearest_neighbor_radial%nn_list(inions)%neighbors_types - 1, &
                  this%nearest_neighbor_radial%nn_list(inions)%neighbor_dist,       &
                  -this%nearest_neighbor_radial%nn_list(inions)%neighbor_connect,   &
                  "2-body"//C_NULL_CHAR)
            end if
            ! fill 3-body neighbor list
            if (this%w2 > 0 ) then
               ! Mind: counting starts at zero in C++ -> minus 1
               call fill_neighbor_list_single_atom(                                  &
                  this%handle,                                                       &
                  size(this%nearest_neighbor_angular%nn_list(inions)%neighbors),     &
                  number_element,                                                    &
                  this%itype(inions) - 1,                                            &
                  this%nearest_neighbor_angular%nn_list(inions)%neighbors - 1,       &
                  this%nearest_neighbor_angular%nn_list(inions)%neighbors_types - 1, &
                  this%nearest_neighbor_angular%nn_list(inions)%neighbor_dist,       &
                  -this%nearest_neighbor_angular%nn_list(inions)%neighbor_connect,   &
                  "3-body"//C_NULL_CHAR)
            end if
            number_element = number_element + 1
         end if
      end do
      !call copy_neighbor_arrays(this%handle)

   end subroutine vaspml_fill_neighbor_list_arrays

   !================================================================================================
   !> rescale 3xN atoms array by a 3x3 matrix
   !> can be used to convert atomic positions from direct to cartesian coordinates and vice versa
   !>
   !> posion atomic positions stored as (xyz, number ions)
   !> lattice defines the matrix by which the positions are multiplied
   !================================================================================================
   subroutine vaspml_rescale_atoms(this, posion, lattice)

      class(vaspml_wrapper), intent(inout) :: this
      real(kind=q), intent(in)             :: posion(:, :)
      real(kind=q), intent(in)             :: lattice(:, :)

      integer :: inions

      do inions = 1, this%nions
         this%posion(1:3, inions) =               &
            posion(1, inions) * lattice(1:3, 1) + &
            posion(2, inions) * lattice(1:3, 2) + &
            posion(3, inions) * lattice(1:3, 3)
      enddo

   end subroutine vaspml_rescale_atoms

   !================================================================================================
   !> determine how many atoms each thread on the c++ side will get. Also determine how many
   !> atoms per type per thread will be passed to the c++ code.
   !>
   !> computes
   !> number_atoms_per_node storing the total number of atoms after round robin distribution
   !> is done
   !> number_atoms_per_type_per_node the total number of atoms per type after round robin
   !> distribution of atoms.
   !================================================================================================
   subroutine vaspml_define_elements_per_thread(this)

      class(vaspml_wrapper), intent(inout) :: this

      integer :: inions_tmp
      integer :: inions

      if (allocated(this%number_atoms_per_type_per_node)) then
         deallocate(this%number_atoms_per_type_per_node)
      end if
      allocate(this%number_atoms_per_type_per_node(1:this%ntype))
      this%number_atoms_per_type_per_node = 0
      this%number_atoms_per_node = 0
      do inions_tmp = 1, this%nions, this%communicator%ncpu
         inions = inions_tmp + this%communicator%node_me
         if (inions <= this%nions) then
            this%number_atoms_per_type_per_node(this%itype(inions)) = &
               this%number_atoms_per_type_per_node(this%itype(inions)) + 1
            this%number_atoms_per_node = this%number_atoms_per_node + 1
         end if
      end do

   end subroutine vaspml_define_elements_per_thread

   subroutine vaspml_fill_force_array(this, tifor)

      use ml_ff_constant, only : funit
      use mpi_help
      class(vaspml_wrapper), intent(inout) :: this
      real(kind=q), intent(inout)          :: tifor(:, :)

      integer :: local_index
      integer :: ion_index_tmp
      integer :: ion_index
      integer :: ierror

      tifor = 0.0_q
      local_index = 0
      do ion_index_tmp = 1, this%nions, this%communicator%ncpu
         ion_index = ion_index_tmp + this%communicator%node_me
         if (ion_index <= this%nions) then
            call fill_force_single_atom(this%handle,   &
                                        this%nions,    &
                                        ion_index - 1, &
                                        local_index,   &
                                        tifor(:, :))
            local_index = local_index + 1
         end if
      end do
      tifor(:, :) = tifor(:, :) * funit
      !call compute_forces(this%handle, this%nions, tifor)
      ! TODO Remove direct calls to MPI, use wrappers!
      call MPI_Allreduce(MPI_IN_PLACE,               &
                         tifor(1, 1),                &
                         3 * this%nions,             &
                         MPI_DOUBLE,                 &
                         MPI_SUM,                    &
                         this%communicator%mpi_comm, &
                         ierror)

   end subroutine vaspml_fill_force_array

   subroutine vaspml_get_energy(this, toten)

      use ml_ff_constant, only : eunit
      use mpi_help
      class(vaspml_wrapper), intent(inout) :: this
      real(kind=q), intent(inout)          :: toten

      integer :: ierror

      toten = 0.0_q
      call get_potential_energy(this%handle, toten)
      toten = toten * eunit
      ! TODO Remove direct calls to MPI, use wrappers!
      call MPI_Allreduce(MPI_IN_PLACE,               &
                         toten,                      &
                         1,                          &
                         MPI_DOUBLE,                 &
                         MPI_SUM,                    &
                         this%communicator%mpi_comm, &
                         ierror)

   end subroutine vaspml_get_energy

   subroutine vaspml_get_stress_tensor(this, tsif)

      use ml_ff_constant, only : sunit, evtoj
      use mpi_help
      class(vaspml_wrapper), intent(inout) :: this
      real(kind=q), intent(inout)          :: tsif(1:3, 1:3)

      integer :: ierror
      real(kind=q) :: fac

      tsif = 0.0_q
      call get_stress_tensor_vaspml(this%handle, tsif)
      fac = evtoj * 1d22 / this%omega / autoa**3
      tsif = tsif * sunit / fac
      ! TODO Remove direct calls to MPI, use wrappers!
      call MPI_Allreduce(MPI_IN_PLACE,               &
                         tsif(1, 1),                 &
                         9,                          &
                         MPI_DOUBLE,                 &
                         MPI_SUM,                    &
                         this%communicator%mpi_comm, &
                         ierror)

   end subroutine vaspml_get_stress_tensor

   !================================================================================================
   !> Convert a given logical to an integer.
   !>
   !> Unfortunately, c_bool is unusable because it translates to C99 type _Bool which is generally
   !> not identical to C++ bool. Instead we use c_int8_t here with the mapping:
   !> * .TRUE.  <---> 1
   !> * .FALSE. <---> 0
   !> With int8_t instead of int on the C++ side we avoid mixing up integers and logicals, i.e.,
   !> int8_t in the wrappers signals that the variable is actually a logical on the Fortran side.
   !================================================================================================
   function f2c_logical(input) result(output)

      integer(kind=c_int8_t) :: output
      logical, intent(in)    :: input

      if (input) then
         output = 1
      else
         output = 0
      end if

   end function f2c_logical

   !================================================================================================
   !> Convert deferred-length character string to one-dimensional character array (with len=1, to
   !> pass to C).
   !>
   !> Replace Fortran new_line() character with C_NEW_LINE (not sure whether there could be a
   !> difference). Also adds null terminator character at the end.
   !================================================================================================
   function f2c_dlcstring(input) result(output)

      character(kind=c_char), allocatable       :: output(:)
      character(len=:), intent(in), allocatable :: input

      integer :: i

      allocate(output(len(input) + 1))
      do i = 1, len(input)
         if (input(i:i) == new_line("n")) then
            output(i) = C_NEW_LINE
         else
            output(i) = input(i:i)
         end if
      end do
      output(len(input) + 1) = C_NULL_CHAR

   end function f2c_dlcstring

   !================================================================================================
   !> Convert one-dimensional array of fixed-length character type to one-dimensional
   !> character array (with len=1, to pass to C) with specified separator.
   !>
   !> Replace Fortran new_line() character with C_NEW_LINE (not sure whether there could be a
   !> difference). Also adds null terminator character at the end.
   !================================================================================================
   function f2c_character1d(input, separator) result(output)

      character(kind=c_char), allocatable :: output(:)
      character(len=*), intent(in)        :: input(:)
      character(len=*), intent(in)        :: separator

      integer :: i, j, k
      integer :: nchars
      integer :: nitems
      integer :: nalloc
      integer :: nsep

      nchars = len(input)
      nitems = size(input)
      nsep = len(separator)

      nalloc = nitems * nchars + 1
      if (nitems > 1) nalloc = nalloc + (nitems - 1) * nsep
      allocate(output(nalloc))

      k = 1
      do i = 1, nitems
         do j = 1, nchars
            if (input(i)(j:j) == new_line("n")) then
               output(k) = C_NEW_LINE
            else
               output(k) = input(i)(j:j)
            end if
            k = k + 1
         end do
         if (i < nitems) then
            do j = 1, nsep
               if (separator(j:j) == new_line("n")) then
                  output(k) = C_NEW_LINE
               else
                  output(k) = separator(j:j)
               end if
               k = k + 1
            end do
         end if
      end do
      output(nalloc) = C_NULL_CHAR

   end function f2c_character1d

end module mvaspml

#endif
