#include "symbol.inc"

!> Wrapper around the DFTD4 and simple-DFTD3 libraries
!>
!> Translates the VASP structural data into a format compatible with the DFTD4
!> and simple-DFTD3 libraries. DFTD4 and simple-DFTD3 are external libraries and need to be
!> installed separately. Please check the
!> [documentation](https://github.com/dftd4/dftd4) and [documentation](https://dftd3.readthedocs.io)
!> for more information.
module vdw_sd3_d4

   use constant, only: AuToA, RyToEv
   use prec, only: q
   use base, only: in_struct
   use tutor, only: vtutor
#ifdef DFTD4
   use dftd4, only: rational_damping_param
#endif
#ifdef SDFTD3
   use dftd3, only: d3_param, &
                    rational_d3_damping_param=>rational_damping_param, &
                    zero_d3_damping_param=>zero_damping_param, &
                    mzero_d3_damping_param=>mzero_damping_param, &
                    optimizedpower_d3_damping_param=>optimizedpower_damping_param
   use mctc_env, only : error_type
#endif

   implicit none

   integer, parameter :: unknown_functional = 1, unknown_atom = -1, vdw_dftd4 = 13, vdw_sdftd3 = 15

#ifdef DFTD4
   logical, parameter :: dftd4_available = .true.
   character(len=50) :: dftd4_model='d4'
#else
   logical, parameter :: dftd4_available = .false.
   type rational_damping_param
   end type rational_damping_param
#endif
#ifdef SDFTD3
   logical, parameter :: sdftd3_available = .true.
   character(len=50) :: sdftd3_damping=""
#else
   logical, parameter :: sdftd3_available = .false.
   type rational_d3_damping_param
   end type rational_d3_damping_param
   type zero_d3_damping_param
   end type zero_d3_damping_param
   type mzero_d3_damping_param
   end type mzero_d3_damping_param
   type optimizedpower_d3_damping_param
   end type optimizedpower_d3_damping_param
#endif

   type vdw_settings
      integer ivdw
      real(q) vdw_cnradius
      real(q) vdw_radius
      character(len=50) :: xc_functional
      logical l_user_defined
      type(in_struct) io
      type(rational_damping_param) parameters
      type(rational_d3_damping_param) parameters_rational_d3
      type(zero_d3_damping_param) parameters_zero_d3
      type(mzero_d3_damping_param) parameters_mzero_d3
      type(optimizedpower_d3_damping_param) parameters_optimizedpower_d3
   end type vdw_settings

   type vdw_structure
      real(q), pointer :: cell(:,:) => null()
      integer, pointer :: types(:) => null()
      character(len=2), pointer :: elements(:) => null()
      real(q), pointer :: positions(:,:) => null()
   end type vdw_structure

   type vdw_interaction
      real(q) energy, stress(3,3)
      real(q), allocatable :: forces(:,:)
   end type vdw_interaction

   interface vdw_structure
      module procedure make_structure
   end interface vdw_structure

   private
   public default_settings, dftd4_available, sdftd3_available, unknown_functional, &
      vdw_settings, vdw_dftd4, vdw_sdftd3, vdw_forces_d4, vdw_forces_sd3, &
      vdw_interaction, vdw_structure, vdw_reader, vdw_writer

contains

   !> Merge information from various types used by VASP into a single type used
   !> by the vdW routine.
   type(vdw_structure) function make_structure(latt_cur, dyn, t_info, elem) result (structure)
      use lattice, only: latt
      use poscar, only: dynamics, type_info
      type(latt), target, intent(in) :: latt_cur
      type(type_info), target, intent(in) :: t_info
      type(dynamics), target, intent(in) :: dyn
      character(len=2), target, intent(in) :: elem(:)
      structure%cell => latt_cur%a
      structure%types => t_info%iTyp
      structure%elements => elem
      structure%positions => dyn%posIon
   end function make_structure

   !> Compute the vdw interaction with DFT-D4
   type(vdw_interaction) function vdw_forces_d4(io, settings, structure) result (results)
      use base, only: in_struct
#ifdef DFTD4
      use dftd4, only: structure_type, d4_model, realspace_cutoff, new_d4_model, get_dispersion
#ifndef DFTD4_API_V3
      use dftd4, only: d4s_model, new_d4s_model
      use mctc_env, only : error_type
#endif
#endif
      type(in_struct), intent(in) :: io
      type(vdw_settings), intent(in) :: settings
      type(vdw_structure), intent(in) :: structure
#ifdef DFTD4
      type(structure_type) struct
      type(d4_model) model_d4
#ifndef DFTD4_API_V3
      type(d4s_model) model_d4s
      type(error_type), allocatable :: error
#endif
      !
      allocate(results%forces(3, size(structure%types)))
      struct = convert_structure_to_dftd4(structure)
#ifdef DFTD4_API_V3
      if (dftd4_model=='d4') then
         call new_d4_model(model_d4, struct)
      else
         call vtutor%error("The DFTD4_MODEL tag has a wrong value. The possible value is 'd4'.")
      endif
      call get_dispersion(struct, model_d4, settings%parameters, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
         results%energy, results%forces, results%stress)
#else
      if (dftd4_model=='d4') then
         call new_d4_model(error, model_d4, struct)
      elseif (dftd4_model=='d4s') then
         call new_d4s_model(error, model_d4s, struct)
      else
         call vtutor%error("The DFTD4_MODEL tag has a wrong value. The possible values are 'd4' and 'd4s'.")
      endif
      if (allocated(error)) then
         call vtutor%error("DFTD4 init error:"//trim(error%message))
      endif
      if (dftd4_model=='d4') then
         call get_dispersion(struct, model_d4, settings%parameters, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
            results%energy, results%forces, results%stress)
      elseif (dftd4_model=='d4s') then
         call get_dispersion(struct, model_d4s, settings%parameters, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
            results%energy, results%forces, results%stress)
      endif
#endif
      call convert_units_of_results(results)

  1001 FORMAT(/,"  DFT-D4 for vdW energy calculation",&
              /,"  ---------------------------------------------------------------------------------------------")
  1002 FORMAT("  Edisp (eV):",F11.5)
      if (io%iu6>=0) then
         write(io%iu6,1001)
         write(io%iu6,1002) results%energy
      endif

#else
      call vtutor%error("The code was not compiled with DFTD4 support.")
#endif
   end function vdw_forces_d4

   !> Compute the vdw interaction with simple-DFTD3
   type(vdw_interaction) function vdw_forces_sd3(io, settings, structure) result (results)
      use base, only: in_struct
#ifdef SDFTD3
      use dftd3, only : d3_model, realspace_cutoff, new_d3_model, get_dispersion
      use mctc_io, only : structure_type
#endif
      type(in_struct), intent(in) :: io
      type(vdw_settings), intent(in) :: settings
      type(vdw_structure), intent(in) :: structure
#ifdef SDFTD3
      type(structure_type) struct
      type(d3_model) model
      !
      allocate(results%forces(3, size(structure%types)))
      struct = convert_structure_to_sdftd3(structure)
      call new_d3_model(model, struct)

      if ((sdftd3_damping=='rational').OR.(sdftd3_damping=='mrational')) then
         call get_dispersion(struct, model, settings%parameters_rational_d3, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
            results%energy, results%forces, results%stress)
      elseif (sdftd3_damping=='zero') then
         call get_dispersion(struct, model, settings%parameters_zero_d3, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
            results%energy, results%forces, results%stress)
      elseif (sdftd3_damping=='mzero') then
         call get_dispersion(struct, model, settings%parameters_mzero_d3, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
            results%energy, results%forces, results%stress)
      elseif (sdftd3_damping=='optimizedpower') then
         call get_dispersion(struct, model, settings%parameters_optimizedpower_d3, realspace_cutoff(cn=settings%vdw_cnradius/AuToA,disp2=settings%vdw_radius/AuToA), &
            results%energy, results%forces, results%stress)
      else
         call vtutor%error("The SDFTD3_DAMPING tag has a wrong value. The possible values are 'rational', 'zero', 'mrational', 'mzero' and 'optimizedpower'.")
      endif
      call convert_units_of_results(results)

  1001 FORMAT(/,"  simple-DFT-D3 for vdW energy calculation",&
              /,"  ---------------------------------------------------------------------------------------------")
  1002 FORMAT("  Edisp (eV):",F11.5)
      if (io%iu6>=0) then
         write(io%iu6,1001)
         write(io%iu6,1002) results%energy
      endif

#else
      call vtutor%error("The code was not compiled with simple-DFTD3 support.")
#endif
   end function vdw_forces_sd3

#ifdef DFTD4
   !> Convert VASP data to the format used by DFTD4
   type (structure_type) function convert_structure_to_dftd4(structure) result (struct)
      use dftd4, only: structure_type, new
      type(vdw_structure), intent(in) :: structure
      integer, allocatable :: atomic_numbers(:)
      real(q), allocatable :: positions(:,:)
      real(q) lattice(3,3)
      integer num_atom
      !
      num_atom = size(structure%types)
      allocate(atomic_numbers(num_atom), positions(3, num_atom))
      atomic_numbers = atomic_number(structure%elements(structure%types))
      if (any(atomic_numbers == unknown_atom)) &
         call vtutor%bug("Could not convert the elements to atomic numbers", __FILE__, __LINE__)
      lattice = structure%cell / AuToA
      positions = matmul(lattice, structure%positions)
      call new(struct, atomic_numbers, positions, lattice=lattice)
   end function convert_structure_to_dftd4
#endif

#ifdef SDFTD3
   !> Convert VASP data to the format used by SDFTD3
   type (structure_type) function convert_structure_to_sdftd3(structure) result (struct)
      use mctc_io, only : structure_type, new
      type(vdw_structure), intent(in) :: structure
      integer, allocatable :: atomic_numbers(:)
      real(q), allocatable :: positions(:,:)
      real(q) lattice(3,3)
      integer num_atom
      !
      num_atom = size(structure%types)
      allocate(atomic_numbers(num_atom), positions(3, num_atom))
      atomic_numbers = atomic_number(structure%elements(structure%types))
      if (any(atomic_numbers == unknown_atom)) &
         call vtutor%bug("Could not convert the elements to atomic numbers", __FILE__, __LINE__)
      lattice = structure%cell / AuToA
      positions = matmul(lattice, structure%positions)
      call new(struct, atomic_numbers, positions, lattice=lattice)
   end function convert_structure_to_sdftd3
#endif

#if defined(DFTD4) || defined(SDFTD3)
   !> Map an ion type onto the atomic number
   elemental integer function atomic_number(element)
      use string, only: lowercase
      character(len=2), intent(in) :: element
      character(len=2), parameter :: elements(103) = ['h ', 'he', &
         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne', &
         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', &
         'k ', 'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', &
            'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', &
         'rb', 'sr', 'y ', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', &
            'cd', 'in', 'sn', 'sb', 'te', 'i ', 'xe', &
         'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', &
            'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w ', 're', 'os', 'ir', 'pt', &
            'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', &
         'fr', 'ra', 'ac', 'th', 'pa', 'u ', 'np', 'pu', 'am', 'cm', 'bk', 'cf', &
            'es', 'fm', 'md', 'no', 'lr']
      character(len=2) element_
      integer i
      element_ = lowercase(element)
      atomic_number = unknown_atom
      do i = 1, size(elements)
         if (element_ == elements(i)) then
            atomic_number = i
            exit
         endif
      end do
   end function atomic_number

   !> Convert between atomic units used by DFTD4 and simple-DFTD3 and Angstroem-eV units used by VASP.
   subroutine convert_units_of_results(results)
      type(vdw_interaction), intent(inout) :: results
      real(q), parameter :: HtrToEv = 2 * RyToEv
      results%energy = results%energy * HtrToEv
      results%forces = -results%forces * HtrToEv / AuToA
      results%stress = -results%stress * HtrToEv
   end subroutine convert_units_of_results
#endif

   !> Read all settings in the INCAR file related to the DFTD4 or simple-DFTD3 library
   subroutine vdw_reader(io, xc_functional, settings)
      use base, only: in_struct
      use reader_tags, only: writexmlincar, process_incar, open_incar_if_found, close_incar_if_found
      type(in_struct), intent(in) :: io
      character(len=*), intent(in) :: xc_functional
      type(vdw_settings), intent(out) :: settings
      integer ierr, ivdw
      logical lopen
      !
      ivdw = 0
      call open_incar_if_found(io%iu5, lopen)
      call process_incar(lopen, io%iu0, io%iu5, 'IVDW', ivdw, ierr, writexmlincar)
      settings%l_user_defined = .false.
      if (ivdw == 13) then
         call read_dftd4_tags(settings,xc_functional)
      elseif (ivdw == 15) then
        call read_sdftd3_tags(settings,xc_functional)
      endif
      call close_incar_if_found(io%iu5)
      settings%ivdw = ivdw
      settings%io = io
      !
   contains
      !
      subroutine read_dftd4_tags(settings,xc_functional_in)
         type(vdw_settings),intent(out) :: settings
         character(len=*),intent(in) :: xc_functional_in
#ifdef DFTD4
         character(len=50) :: xc_functional
         character(len=:), allocatable :: user_defined
         real(q) rdum
         logical unknown
         !
         xc_functional = xc_functional_in
         call process_incar(lopen, io%iu0, io%iu5, 'DFTD4_XC', xc_functional, 50, ierr, writexmlincar)
         if (ierr==0) call lower(xc_functional)
         settings = default_settings(xc_functional, ierr)
         unknown = ierr == unknown_functional
         user_defined = ''
         !
         call process_incar(lopen, io%iu0, io%iu5, 'DFTD4_MODEL', dftd4_model, 50, ierr, writexmlincar)
         call lower(dftd4_model)
         !
         associate (param => settings%parameters)
            call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
            if (ierr == 0) then
               call check_interval('VDW_S6', [0.0_q, 10.0_q], rdum, param%s6, user_defined)
            elseif (unknown) then
               call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
            endif
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
            if (ierr == 0) then
               call check_interval('VDW_S8', [-10.0_q, 10.0_q], rdum, param%s8, user_defined)
            elseif (unknown) then
               call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
            endif
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_A1', rdum, ierr, writexmlincar)
            if (ierr == 0) then
               call check_interval('VDW_A1', [-10.0_q, 10.0_q], rdum, param%a1, user_defined)
            elseif (unknown) then
               call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
            endif
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_A2', rdum, ierr, writexmlincar)
            if (ierr == 0) then
               call check_interval('VDW_A2', [0.0_q, 20.0_q], rdum, param%a2, user_defined)
            elseif (unknown) then
               call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
            endif
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_S9', rdum, ierr, writexmlincar)
            if (ierr == 0) call check_interval('VDW_S9', [0.0_q, 1.0_q], rdum, param%s9, user_defined)
         end associate
         !
         ! The default values in the DFTD4 library are vdw_cnradius=30 Bohr and vdw_radius=60 Bohr.
         settings%vdw_cnradius = 30._q*AuToA
         settings%vdw_radius = 60._q*AuToA
         call process_incar(lopen, io%iu0, io%iu5, 'VDW_CNRADIUS', settings%vdw_cnradius, ierr, writexmlincar)
         call process_incar(lopen, io%iu0, io%iu5, 'VDW_RADIUS', settings%vdw_radius, ierr, writexmlincar)
         !
         if (len_trim(user_defined) > 0) settings%l_user_defined = .true.
#else
         call vtutor%error("The code was not compiled with DFTD4 support. Please &
            recompile the code adding -DDFTD4 to the CPP_OPTIONS and add the &
            necessary include and link parameters.")
#endif
      end subroutine read_dftd4_tags

      subroutine read_sdftd3_tags(settings,xc_functional_in)
         type(vdw_settings),intent(out) :: settings
         character(len=*),intent(in) :: xc_functional_in
#ifdef SDFTD3
         character(len=50) :: xc_functional
         character(len=:), allocatable :: user_defined
         real(q) rdum
         logical unknown
         !
         xc_functional = xc_functional_in
         call process_incar(lopen, io%iu0, io%iu5, 'SDFTD3_XC', xc_functional, 50, ierr, writexmlincar)
         if (ierr==0) call lower(xc_functional)
         user_defined = ''
         !
         call process_incar(lopen, io%iu0, io%iu5, 'SDFTD3_DAMPING', sdftd3_damping, 50, ierr, writexmlincar)
         if (ierr==0) then
            call lower(sdftd3_damping)
         else
            call vtutor%error("SDFTD3_DAMPING has to be specified. The possible values are 'rational', 'zero', 'mrational', 'mzero' and 'optimizedpower'.")
         endif

         if (sdftd3_damping=='rational') then

            settings = default_settings_rational_d3(xc_functional, ierr)
            unknown = ierr == unknown_functional
            associate (param => settings%parameters_rational_d3)
               call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S6', [0.0_q, 10.0_q], rdum, param%s6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S8', [-10.0_q, 10.0_q], rdum, param%s8, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_A1', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_A1', [-10.0_q, 10.0_q], rdum, param%a1, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_A2', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_A2', [0.0_q, 20.0_q], rdum, param%a2, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S9', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_S9', [0.0_q, 1.0_q], rdum, param%s9, user_defined)
            end associate

         elseif (sdftd3_damping=='zero') then

            settings = default_settings_zero_d3(xc_functional, ierr)
            unknown = ierr == unknown_functional
            associate (param => settings%parameters_zero_d3)
               call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S6', [0.0_q, 10.0_q], rdum, param%s6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8 and VDW_SR have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S8', [-10.0_q, 10.0_q], rdum, param%s8, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8 and VDW_SR have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_SR', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_SR', [0.0_q, 10.0_q], rdum, param%rs6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8 and VDW_SR have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_SR8', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_SR8', [0.0_q, 10.0_q], rdum, param%rs8, user_defined)
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S9', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_S9', [0.0_q, 1.0_q], rdum, param%s9, user_defined)
            end associate

         elseif (sdftd3_damping=='mrational') then

            settings = default_settings_mrational_d3(xc_functional, ierr)
            unknown = ierr == unknown_functional
            associate (param => settings%parameters_rational_d3)
               call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S6', [0.0_q, 10.0_q], rdum, param%s6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S8', [-10.0_q, 10.0_q], rdum, param%s8, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_A1', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_A1', [-10.0_q, 10.0_q], rdum, param%a1, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_A2', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_A2', [0.0_q, 20.0_q], rdum, param%a2, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1 and VDW_A2 have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S9', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_S9', [0.0_q, 1.0_q], rdum, param%s9, user_defined)
            end associate

         elseif (sdftd3_damping=='mzero') then

            settings = default_settings_mzero_d3(xc_functional, ierr)
            unknown = ierr == unknown_functional
            associate (param => settings%parameters_mzero_d3)
               call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S6', [0.0_q, 10.0_q], rdum, param%s6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_SR and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S8', [-10.0_q, 10.0_q], rdum, param%s8, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_SR and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_SR', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_SR', [0.0_q, 10.0_q], rdum, param%rs6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_SR and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_BETA', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_BETA', [0.0_q, 20.0_q], rdum, param%bet, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_SR and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_SR8', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_SR8', [0.0_q, 10.0_q], rdum, param%rs8, user_defined)
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S9', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_S9', [0.0_q, 1.0_q], rdum, param%s9, user_defined)
            end associate

         elseif (sdftd3_damping=='optimizedpower') then

            settings = default_settings_optimizedpower_d3(xc_functional, ierr)
            unknown = ierr == unknown_functional
            associate (param => settings%parameters_optimizedpower_d3)
               call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S6', [0.0_q, 10.0_q], rdum, param%s6, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1, VDW_A2 and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_S8', [-10.0_q, 10.0_q], rdum, param%s8, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1, VDW_A2 and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_A1', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_A1', [-10.0_q, 10.0_q], rdum, param%a1, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1, VDW_A2 and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_A2', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_A2', [0.0_q, 20.0_q], rdum, param%a2, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1, VDW_A2 and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_BETA', rdum, ierr, writexmlincar)
               if (ierr == 0) then
                  call check_interval('VDW_BETA', [0.0_q, 20.0_q], rdum, param%bet, user_defined)
               elseif (unknown) then
                  call vtutor%error('The xc functional is not supported. VDW_S6, VDW_S8, VDW_A1, VDW_A2 and VDW_BETA have to be defined in INCAR.')
               endif
               call process_incar(lopen, io%iu6, io%iu5, 'VDW_S9', rdum, ierr, writexmlincar)
               if (ierr == 0) call check_interval('VDW_S9', [0.0_q, 1.0_q], rdum, param%s9, user_defined)
            end associate

         else
            call vtutor%error("The SDFTD3_DAMPING tag is not specified.")
         endif
         !
         ! The default values in the simple-DFTD3 library are vdw_cnradius=40 Bohr and vdw_radius=60 Bohr.
         settings%vdw_cnradius = 40._q*AuToA
         settings%vdw_radius = sqrt(9000._q)*AuToA
         call process_incar(lopen, io%iu0, io%iu5, 'VDW_CNRADIUS', settings%vdw_cnradius, ierr, writexmlincar)
         call process_incar(lopen, io%iu0, io%iu5, 'VDW_RADIUS', settings%vdw_radius, ierr, writexmlincar)
         !
         if (len_trim(user_defined) > 0) settings%l_user_defined = .true.
#else
         call vtutor%error("The code was not compiled with simple-DFTD3 support. Please &
            recompile the code adding -DSDFTD3 to the CPP_OPTIONS and add the &
            necessary include and link parameters.")
#endif
      end subroutine read_sdftd3_tags
      !
      subroutine check_interval(tag, sane_interval, new_value, use_value, user_defined)
         use string, only: str
         character(len=*), intent(in) :: tag
         character(len=10) :: sane_interval_min, sane_interval_max
         real(q), intent(in) :: sane_interval(2), new_value
         real(q), intent(inout) :: use_value
         character(len=:), allocatable, intent(inout) :: user_defined
         if ((minval(sane_interval) <= new_value).and.(new_value <= maxval(sane_interval))) then
            user_defined = user_defined // '\n' // tag // ' = ' // str(new_value)
            use_value = new_value
         else
            write(sane_interval_min, "(F4.0)") minval(sane_interval)
            write(sane_interval_max, "(F4.0)") maxval(sane_interval)
            call vtutor%error(tag // ' is not reasonable. The value should be in the' // '\n' // &
               'interval [' // trim(sane_interval_min) // ',' // trim(sane_interval_max) // ']')
         endif
      end subroutine check_interval
      !
   end subroutine vdw_reader

   subroutine vdw_writer(settings)
       type(vdw_settings), intent(in) :: settings
       character(len=*), parameter :: fmt_int = '(3x,a,t13,i6,4x,a)'
       integer unit_
       unit_ = settings%io%iu6
       if (unit_ < 0) &
           return
       write(unit_,*) 'Van der Waals corrections'
       write(unit_,fmt_int) 'IVDW    =', settings%ivdw, 'specifies the selected vdW correction'
       if (settings%ivdw == 13) then
#ifdef DFTD4
          call write_dftd4_info(unit_, settings)
#else
          write(unit_,'(3x,a)') 'DFT-D4  = not available'
#endif
       elseif (settings%ivdw == 15) then
#ifdef SDFTD3
          if ((sdftd3_damping=='rational').OR.(sdftd3_damping=='mrational')) then
             call write_sdftd3_rational_info(unit_, settings)
          elseif (sdftd3_damping=='zero') then
             call write_sdftd3_zero_info(unit_, settings)
          elseif (sdftd3_damping=='mzero') then
             call write_sdftd3_mzero_info(unit_, settings)
          elseif (sdftd3_damping=='optimizedpower') then
             call write_sdftd3_optimizedpower_info(unit_, settings)
          else
             call vtutor%error("The SDFTD3_DAMPING tag is not specified.")
          endif
#else
          write(unit_,'(3x,a)') 'simple-DFT-D3  = not available'
#endif
       endif
       write(unit_,'(a)')
   end subroutine vdw_writer

#ifdef DFTD4
   subroutine write_dftd4_info(unit_, settings)
       use dftd4, only: get_dftd4_version
       integer, intent(in) :: unit_
       type(vdw_settings), intent(in) :: settings
       character(len=:), allocatable :: version
       character(len=*), parameter :: fmt_real = '(3x,a,t20,f12.8,2x,a)'
       call get_dftd4_version(string=version)
       write(unit_,'(3x,a)') 'Using the DFT-D4 (v' // version // ') external library'
       write(unit_,'(3x,a)') 'DFTD4_MODEL    = ' // trim(dftd4_model)
       if (settings%l_user_defined) then
          write(unit_,'(3x,a)') 'The DFT-D4 parameters are defined by the user'
       else
          write(unit_,'(3x,a)') 'The DFT-D4 parameters are for DFTD4_XC = ' // trim(settings%xc_functional) // ' (can be modified):'
       endif
       write(unit_,fmt_real) 'VDW_S6         =', settings%parameters%s6, 's6, scaling of the dipole-dipole dispersion'
       write(unit_,fmt_real) 'VDW_S8         =', settings%parameters%s8, 's8, scaling of the dipole-quadrupole dispersion'
       write(unit_,fmt_real) 'VDW_A1         =', settings%parameters%a1, 'a1, scaling of the critical radii'
       write(unit_,fmt_real) 'VDW_A2         =', settings%parameters%a2, 'a2, offset of the critical radii'
       write(unit_,fmt_real) 'VDW_S9         =', settings%parameters%s9, 's9, scaling of the three-body dispersion'
       write(unit_,fmt_real) 'VDW_RADIUS     =', settings%vdw_radius,    'two-body interaction cutoff (in Angstroms)'
       write(unit_,fmt_real) 'VDW_CNRADIUS   =', settings%vdw_cnradius,  'CN, coordination number cutoff (in Angstroms)'
       write(unit_,'(3x,a)') 'Fixed parameters (cannot be modified):'
       write(unit_,fmt_real) 'alpha          =', settings%parameters%alp,'alpha, exponent of the zero damping (three-body term)'
   end subroutine write_dftd4_info
#endif

#ifdef SDFTD3
   subroutine write_sdftd3_rational_info(unit_, settings)
       use dftd3, only: get_dftd3_version
       integer, intent(in) :: unit_
       type(vdw_settings), intent(in) :: settings
       character(len=:), allocatable :: version
       character(len=*), parameter :: fmt_real = '(3x,a,t20,f12.8,2x,a)'
       call get_dftd3_version(string=version)
       write(unit_,'(3x,a)') 'Using the simple-DFT-D3 (v' // version // ') external library'
       write(unit_,'(3x,a)') 'SDFTD3_DAMPING = ' // trim(sdftd3_damping)
       if (settings%l_user_defined) then
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are defined by the user'
       else
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are for SDFTD3_XC = ' // trim(settings%xc_functional) // ' (can be modified):'
       endif
       write(unit_,fmt_real) 'VDW_S6         =', settings%parameters_rational_d3%s6,   's6, scaling of the dipole-dipole dispersion'
       write(unit_,fmt_real) 'VDW_S8         =', settings%parameters_rational_d3%s8,   's8, scaling of the dipole-quadrupole dispersion'
       write(unit_,fmt_real) 'VDW_A1         =', settings%parameters_rational_d3%a1,   'a1, scaling of the critical radii'
       write(unit_,fmt_real) 'VDW_A2         =', settings%parameters_rational_d3%a2,   'a2, offset of the critical radii'
       write(unit_,fmt_real) 'VDW_S9         =', settings%parameters_rational_d3%s9,   's9, scaling of the three-body dispersion'
       write(unit_,fmt_real) 'VDW_RADIUS     =', settings%vdw_radius,                  'two-body interaction cutoff (in Angstroms)'
       write(unit_,fmt_real) 'VDW_CNRADIUS   =', settings%vdw_cnradius,                'CN, coordination number cutoff (in Angstroms)'
       write(unit_,'(3x,a)') 'Fixed parameters (cannot be modified):'
       write(unit_,fmt_real) 'alpha          =', settings%parameters_rational_d3%alp+2,'alpha, exponent of the zero damping (three-body term)'
   end subroutine write_sdftd3_rational_info

   subroutine write_sdftd3_zero_info(unit_, settings)
       use dftd3, only: get_dftd3_version
       integer, intent(in) :: unit_
       type(vdw_settings), intent(in) :: settings
       character(len=:), allocatable :: version
       character(len=*), parameter :: fmt_real = '(3x,a,t20,f12.8,2x,a)'
       call get_dftd3_version(string=version)
       write(unit_,'(3x,a)') 'Using the simple-DFT-D3 (v' // version // ') external library'
       write(unit_,'(3x,a)') 'SDFTD3_DAMPING = ' // trim(sdftd3_damping)
       if (settings%l_user_defined) then
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are defined by the user'
       else
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are for SDFTD3_XC = ' // trim(settings%xc_functional) // ' (can be modified):'
       endif
       write(unit_,fmt_real) 'VDW_S6         =', settings%parameters_zero_d3%s6,   's6, scaling of the dipole-dipole dispersion'
       write(unit_,fmt_real) 'VDW_S8         =', settings%parameters_zero_d3%s8,   's8, scaling of the dipole-quadrupole dispersion'
       write(unit_,fmt_real) 'VDW_SR         =', settings%parameters_zero_d3%rs6,  'sr6, scaling of the dipole-dipole damping'
       write(unit_,fmt_real) 'VDW_SR8        =', settings%parameters_zero_d3%rs8,  'sr8, scaling of the dipole-quadrupole damping'
       write(unit_,fmt_real) 'VDW_S9         =', settings%parameters_zero_d3%s9,   's9, scaling of the three-body dispersion'
       write(unit_,fmt_real) 'VDW_RADIUS     =', settings%vdw_radius,              'two-body interaction cutoff (in Angstroms)'
       write(unit_,fmt_real) 'VDW_CNRADIUS   =', settings%vdw_cnradius,            'CN, coordination number cutoff (in Angstroms)'
       write(unit_,'(3x,a)') 'Fixed parameters (cannot be modified):'
       write(unit_,fmt_real) 'alpha6         =', settings%parameters_zero_d3%alp,  'alpha6, exponent of the zero damping (two-body term)'
       write(unit_,fmt_real) 'alpha          =', settings%parameters_zero_d3%alp+2,'alpha, exponent of the zero damping (three-body term)'
   end subroutine write_sdftd3_zero_info

   subroutine write_sdftd3_mzero_info(unit_, settings)
       use dftd3, only: get_dftd3_version
       integer, intent(in) :: unit_
       type(vdw_settings), intent(in) :: settings
       character(len=:), allocatable :: version
       character(len=*), parameter :: fmt_real = '(3x,a,t20,f12.8,2x,a)'
       call get_dftd3_version(string=version)
       write(unit_,'(3x,a)') 'Using the simple-DFT-D3 (v' // version // ') external library'
       write(unit_,'(3x,a)') 'SDFTD3_DAMPING = ' // trim(sdftd3_damping)
       if (settings%l_user_defined) then
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are defined by the user'
       else
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are for SDFTD3_XC = ' // trim(settings%xc_functional) // ' (can be modified):'
       endif
       write(unit_,fmt_real) 'VDW_S6         =', settings%parameters_mzero_d3%s6,   's6, scaling of the dipole-dipole dispersion'
       write(unit_,fmt_real) 'VDW_S8         =', settings%parameters_mzero_d3%s8,   's8, scaling of the dipole-quadrupole dispersion'
       write(unit_,fmt_real) 'VDW_SR         =', settings%parameters_mzero_d3%rs6,  'sr6, scaling of the dipole-dipole damping'
       write(unit_,fmt_real) 'VDW_SR8        =', settings%parameters_mzero_d3%rs8,  'sr8, scaling of the dipole-quadrupole damping'
       write(unit_,fmt_real) 'VDW_BETA       =', settings%parameters_mzero_d3%bet,  'beta, offset for damping radius'
       write(unit_,fmt_real) 'VDW_S9         =', settings%parameters_mzero_d3%s9,   's9, scaling of the three-body dispersion'
       write(unit_,fmt_real) 'VDW_RADIUS     =', settings%vdw_radius,               'two-body interaction cutoff (in Angstroms)'
       write(unit_,fmt_real) 'VDW_CNRADIUS   =', settings%vdw_cnradius,             'CN, coordination number cutoff (in Angstroms)'
       write(unit_,'(3x,a)') 'Fixed parameters (cannot be modified):'
       write(unit_,fmt_real) 'alpha6         =', settings%parameters_mzero_d3%alp,  'alpha6, exponent of the zero damping (two-body term)'
       write(unit_,fmt_real) 'alpha          =', settings%parameters_mzero_d3%alp+2,'alpha, exponent of the zero damping (three-body term)'
   end subroutine write_sdftd3_mzero_info

   subroutine write_sdftd3_optimizedpower_info(unit_, settings)
       use dftd3, only: get_dftd3_version
       integer, intent(in) :: unit_
       type(vdw_settings), intent(in) :: settings
       character(len=:), allocatable :: version
       character(len=*), parameter :: fmt_real = '(3x,a,t20,f12.8,2x,a)'
       call get_dftd3_version(string=version)
       write(unit_,'(3x,a)') 'Using the simple-DFT-D3 (v' // version // ') external library'
       write(unit_,'(3x,a)') 'SDFTD3_DAMPING = ' // trim(sdftd3_damping)
       if (settings%l_user_defined) then
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are defined by the user'
       else
          write(unit_,'(3x,a)') 'The simple-DFT-D3 parameters are for SDFTD3_XC = ' // trim(settings%xc_functional) // ' (can be modified):'
       endif
       write(unit_,fmt_real) 'VDW_S6         =', settings%parameters_optimizedpower_d3%s6,   's6, scaling of the dipole-dipole dispersion'
       write(unit_,fmt_real) 'VDW_S8         =', settings%parameters_optimizedpower_d3%s8,   's8, scaling of the dipole-quadrupole dispersion'
       write(unit_,fmt_real) 'VDW_A1         =', settings%parameters_optimizedpower_d3%a1,   'a1, scaling of the critical radii'
       write(unit_,fmt_real) 'VDW_A2         =', settings%parameters_optimizedpower_d3%a2,   'a2, offset of the critical radii'
       write(unit_,fmt_real) 'VDW_BETA       =', settings%parameters_optimizedpower_d3%bet,  'beta, power for the zero-damping component'
       write(unit_,fmt_real) 'VDW_S9         =', settings%parameters_optimizedpower_d3%s9,   's9, scaling of the three-body dispersion'
       write(unit_,fmt_real) 'VDW_RADIUS     =', settings%vdw_radius,                        'two-body interaction cutoff (in Angstroms)'
       write(unit_,fmt_real) 'VDW_CNRADIUS   =', settings%vdw_cnradius,                      'CN, coordination number cutoff (in Angstroms)'
       write(unit_,'(3x,a)') 'Fixed parameters (cannot be modified):'
       write(unit_,fmt_real) 'alpha          =', settings%parameters_optimizedpower_d3%alp+2,'alpha, exponent of the zero damping (three-body term)'
   end subroutine write_sdftd3_optimizedpower_info
#endif

   type(vdw_settings) function default_settings(xc_functional, ierr) result (settings)
#ifdef DFTD4
      use dftd4, only: get_rational_damping, damping_param
      class(damping_param), allocatable :: param
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef DFTD4
      settings%xc_functional = xc_functional
      call get_rational_damping(xc_functional, param)
      if (.not.allocated(param)) then
          ierr = unknown_functional
          call get_rational_damping('pbe', param)
      endif
      select type (param)
      type is (rational_damping_param)
          settings%parameters = param
      class default
          call vtutor%bug("These damping parameters are not implemented", __FILE__, __LINE__)
      end select
#endif
   end function default_settings

   type(vdw_settings) function default_settings_rational_d3(xc_functional, ierr) result (settings)
#ifdef SDFTD3
      use dftd3, only : get_rational_damping, new_rational_damping
      type(rational_d3_damping_param) :: param
      type(d3_param) :: inp
      type(error_type), allocatable :: error
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef SDFTD3
      settings%xc_functional = xc_functional
      call get_rational_damping(inp,xc_functional,error,s9=0.0_q)
      if (allocated(error)) then
         ierr = unknown_functional
         call get_rational_damping(inp,'pbe',error,s9=0.0_q)
      endif
      call new_rational_damping(param, inp)
      settings%parameters_rational_d3 = param
#endif
   end function default_settings_rational_d3

   type(vdw_settings) function default_settings_zero_d3(xc_functional, ierr) result (settings)
#ifdef SDFTD3
      use dftd3, only : get_zero_damping, new_zero_damping
      type(zero_d3_damping_param) :: param
      type(d3_param) :: inp
      type(error_type), allocatable :: error
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef SDFTD3
      settings%xc_functional = xc_functional
      call get_zero_damping(inp,xc_functional,error,s9=0.0_q)
      if (allocated(error)) then
         ierr = unknown_functional
         call get_zero_damping(inp,'pbe',error,s9=0.0_q)
      endif
      call new_zero_damping(param, inp)
      settings%parameters_zero_d3 = param
#endif
   end function default_settings_zero_d3

   type(vdw_settings) function default_settings_mrational_d3(xc_functional, ierr) result (settings)
#ifdef SDFTD3
      use dftd3, only : get_mrational_damping, new_rational_damping
      type(rational_d3_damping_param) :: param
      type(d3_param) :: inp
      type(error_type), allocatable :: error
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef SDFTD3
      settings%xc_functional = xc_functional
      call get_mrational_damping(inp,xc_functional,error,s9=0.0_q)
      if (allocated(error)) then
         ierr = unknown_functional
         call get_mrational_damping(inp,'pbe',error,s9=0.0_q)
      endif
      call new_rational_damping(param, inp)
      settings%parameters_rational_d3 = param
#endif
   end function default_settings_mrational_d3

   type(vdw_settings) function default_settings_mzero_d3(xc_functional, ierr) result (settings)
#ifdef SDFTD3
      use dftd3, only : get_mzero_damping, new_mzero_damping
      type(mzero_d3_damping_param) :: param
      type(d3_param) :: inp
      type(error_type), allocatable :: error
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef SDFTD3
      settings%xc_functional = xc_functional
      call get_mzero_damping(inp,xc_functional,error,s9=0.0_q)
      if (allocated(error)) then
         ierr = unknown_functional
         call get_mzero_damping(inp,'pbe',error,s9=0.0_q)
      endif
      call new_mzero_damping(param, inp)
      settings%parameters_mzero_d3 = param
#endif
   end function default_settings_mzero_d3

   type(vdw_settings) function default_settings_optimizedpower_d3(xc_functional, ierr) result (settings)
#ifdef SDFTD3
      use dftd3, only : get_optimizedpower_damping, new_optimizedpower_damping
      type(optimizedpower_d3_damping_param) :: param
      type(d3_param) :: inp
      type(error_type), allocatable :: error
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef SDFTD3
      settings%xc_functional = xc_functional
      call get_optimizedpower_damping(inp,xc_functional,error,s9=0.0_q)
      if (allocated(error)) then
         ierr = unknown_functional
         call get_optimizedpower_damping(inp,'pbe',error,s9=0.0_q)
      endif
      call new_optimizedpower_damping(param, inp)
      settings%parameters_optimizedpower_d3 = param
#endif
   end function default_settings_optimizedpower_d3

end module vdw_sd3_d4
