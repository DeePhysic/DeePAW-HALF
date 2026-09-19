#include "symbol.inc"

!> Wrapper around the libMBD library
!>
!> Translates the VASP structural data into a format compatible with the libMBD
!> library. libMBD is an external library and needs to be installed separately.
!> Please check the [documentation](https://github.com/libmbd/libmbd) for more
!> information.
!>
!> The RPA-related options are disabled (see commented lines below) since
!> it is not clear how to use them.
module libmbd

   use prec, only: q
#ifdef LIBMBD
   use mbd, only : mbd_input_t, mbd_calc_t
#endif
   implicit none

   type libmbd_structure
      integer, pointer :: nions => null()
      real(q), pointer :: cell(:,:) => null()
      real(q), pointer :: cell_recip(:,:) => null()
      integer, pointer :: types(:) => null()
      character(len=2), pointer :: elements(:) => null()
      real(q), pointer :: positions(:,:) => null()
   end type libmbd_structure

   type libmbd_interaction
      real(q) energy, stress(3,3)
      real(q), allocatable :: forces(:,:)
!      real(q), allocatable :: spectrum(:), modes(:,:)
   end type libmbd_interaction

   contains

   type(libmbd_structure) function structure_settings_libmbd(latt_cur,dyn,t_info,elem) result (structure)
      use lattice, only: latt
      use poscar, only: dynamics, type_info
      type(latt), target, intent(in) :: latt_cur
      type(dynamics), target, intent(in) :: dyn
      type(type_info), target, intent(in) :: t_info
      character(len=2), target, intent(in) :: elem(:)
      structure%nions => t_info%nions
      structure%cell => latt_cur%a
      structure%cell_recip => latt_cur%b
      structure%types => t_info%iTyp
      structure%elements => elem
      structure%positions => dyn%posIon
   end function structure_settings_libmbd

#ifdef LIBMBD
   subroutine input_settings_libmbd(io,xc,t_info,structure,kpoints,input,ivdw)
      use constant, only: AuToA
      use poscar, only: dynamics, type_info
      use base, only: in_struct
      use reader_tags, only: writexmlincar, process_incar, open_incar_if_found, close_incar_if_found
      use setexm_struct_def, only : xc_info
      use setexm, only: verbalize_dftmethod
      use mkpoints, only: kpoints_struct
      use tutor, only: vtutor
      type(in_struct), intent(in) :: io
      type(xc_info), intent(in) :: xc
      type(type_info), target, intent(in) :: t_info
      type(libmbd_structure), intent(in) :: structure
      type(kpoints_struct), intent(in) :: kpoints
      type(mbd_input_t) :: input
      integer, intent(in) :: ivdw

      integer i, ierr, idum(3)
      real(q) :: rdum, rdumvec(1:t_info%nions)
      character(10) :: dftmethod
      character(40) :: cdum
      logical :: lopen, ldum

      if (.not.allocated(input%atom_types)) allocate(input%atom_types(structure%nions))
      if (.not.allocated(input%coords)) allocate(input%coords(3,structure%nions))
      if (.not.allocated(input%lattice_vectors)) allocate(input%lattice_vectors(3,3))

      call open_incar_if_found(io%iu5, lopen)

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_METHOD', cdum, 40, ierr, writexmlincar)
      if (ierr==0) then
         call lower(cdum)
         input%method=cdum
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_MAX_ATOMS_PER_BLOCK', idum(1), ierr, writexmlincar)
      if (ierr==0) input%max_atoms_per_block=idum(1)

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_LOG_LEVEL', idum(1), ierr, writexmlincar)
      if (ierr==0) input%log_level=idum(1)

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_CALCULATE_FORCES', ldum, ierr, writexmlincar)
      if (ierr==0) input%calculate_forces=ldum

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_CALCULATE_VDW_PARAMS_GRADIENTS', ldum, ierr, writexmlincar)
      if (ierr==0) input%calculate_vdw_params_gradients=ldum

!      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_CALCULATE_SPECTRUM', ldum, ierr, writexmlincar)
!      if (ierr==0) input%calculate_spectrum=ldum

!      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_DO_RPA', ldum, ierr, writexmlincar)
!      if (ierr==0) input%do_rpa=ldum

!      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_RPA_ORDERS', ldum, ierr, writexmlincar)
!      if (ierr==0) input%rpa_orders=ldum

!      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_RPA_RESCALE_EIGS', ldum, ierr, writexmlincar)
!      if (ierr==0) input%rpa_rescale_eigs=ldum

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_N_OMEGA_GRID', idum(1), ierr, writexmlincar)
      if (ierr==0) input%n_omega_grid=idum(1)

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_K_GRID_SHIFT', rdum, ierr, writexmlincar)
      if (ierr==0) input%k_grid_shift=rdum

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_ZERO_NEGATIVE_EIGVALS', ldum, ierr, writexmlincar)
      if (ierr==0) input%zero_negative_eigvals=ldum

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_XC', cdum, 40, ierr, writexmlincar)
      if (ierr==0) then
         call lower(cdum)
         select case (cdum)
            case ('none')
               input%xc=''
            case default
               input%xc=cdum
         end select
      else
         dftmethod=verbalize_dftmethod(xc)
         select case (dftmethod)
            case ('pbe')
               input%xc='pbe'
            case ('pbe0')
               input%xc='pbe0'
            case ('hse03','hse06')
               input%xc='hse'
            case ('am05')
               input%xc='am05'
            case ('revpbe')
               input%xc='revpbe'
         end select
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_TS_D', rdum, ierr, writexmlincar)
      if (ierr==0) then
         if (input%xc=='') then
            input%ts_d=rdum
         else
            call vtutor%error("Error: LIBMBD_TS_D can be specified only for LIBMBD_XC=NONE.")
         endif
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_TS_SR', rdum, ierr, writexmlincar)
      if (ierr==0) then
         if (input%xc=='') then
            input%ts_sr=rdum
         else
            call vtutor%error("Error: LIBMBD_TS_SR can be specified only for LIBMBD_XC=NONE.")
         endif
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_MBD_A', rdum, ierr, writexmlincar)
      if (ierr==0) then
         if (input%xc=='') then
            input%mbd_a=rdum
         else
            call vtutor%error("Error: LIBMBD_MBD_A can be specified only for LIBMBD_XC=NONE.")
         endif
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_MBD_BETA', rdum, ierr, writexmlincar)
      if (ierr==0) then
         if (input%xc=='') then
            input%mbd_beta=rdum
         else
            call vtutor%error("Error: LIBMBD_MBD_BETA can be specified only for LIBMBD_XC=NONE.")
         endif
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_VDW_PARAMS_KIND', cdum, 40, ierr, writexmlincar)
      if (ierr==0) then
         call lower(cdum)
         input%vdw_params_kind=cdum
      endif

      do i=1,structure%nions
         input%atom_types(i) = structure%elements(structure%types(i))
      enddo

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_ALPHA', rdumvec(1:t_info%nions), t_info%nions, ierr, writexmlincar)
      if (ierr==0) then
         if (.not.allocated(input%free_values)) allocate(input%free_values(3,1:t_info%nions))
         input%free_values(1,1:t_info%nions)=rdumvec(1:t_info%nions)
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_C6AU', rdumvec(1:t_info%nions), t_info%nions, ierr, writexmlincar)
      if (ierr==0) then
         if (.not.allocated(input%free_values)) allocate(input%free_values(3,1:t_info%nions))
         input%free_values(2,1:t_info%nions)=rdumvec(1:t_info%nions)
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_R0AU', rdumvec(1:t_info%nions), t_info%nions, ierr, writexmlincar)
      if (ierr==0) then
         if (.not.allocated(input%free_values)) allocate(input%free_values(3,1:t_info%nions))
         input%free_values(3,1:t_info%nions)=rdumvec(1:t_info%nions)
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_K_GRID', idum(1:3), 3, ierr, writexmlincar)
      if (ierr==0) then
         input%k_grid=idum(1:3)
      else
         input%k_grid = [kpoints%nkpx, kpoints%nkpy, kpoints%nkpz]
      endif

      call process_incar(lopen, io%iu0, io%iu5, 'LIBMBD_PARALLEL_MODE', cdum, 40, ierr, writexmlincar)
      if (ierr==0) then
         call lower(cdum)
         input%parallel_mode=cdum
      endif

      call close_incar_if_found(io%iu5)

  1001 FORMAT(/,"  libMBD for vdW energy calculation")
  1002 FORMAT("  IVDW = ",I3)
  1003 FORMAT("  LIBMBD_METHOD = ",A10)
  1004 FORMAT("  LIBMBD_N_OMEGA_GRID = ",I6)
  1005 FORMAT("  LIBMBD_K_GRID_SHIFT = ",F9.3)
  1006 FORMAT("  LIBMBD_XC = ",A10)
  1007 FORMAT("  LIBMBD_TS_D = ",F9.3)
  1008 FORMAT("  LIBMBD_TS_SR = ",F9.3)
  1009 FORMAT("  LIBMBD_MBD_A = ",F9.3)
  1010 FORMAT("  LIBMBD_MBD_BETA = ",F9.3)
  1011 FORMAT("  LIBMBD_VDW_PARAMS_KIND = ",A10)
  1012 FORMAT("  LIBMBD_K_GRID = ",I8,X,I8,X,I8)
  1013 FORMAT("  LIBMBD_ALPHA = ")
  1014 FORMAT("  LIBMBD_C6AU = ")
  1015 FORMAT("  LIBMBD_R0AU = ")

      if (io%iu6>=0) then
         write(io%iu6,1001)
         write(io%iu6,1002) ivdw
         write(io%iu6,1003) input%method
         write(io%iu6,1004) input%n_omega_grid
         write(io%iu6,1005) input%k_grid_shift
         write(io%iu6,1006) input%xc
         write(io%iu6,1007) input%ts_d
         write(io%iu6,1008) input%ts_sr
         write(io%iu6,1009) input%mbd_a
         write(io%iu6,1010) input%mbd_beta
         write(io%iu6,1011) input%vdw_params_kind
         write(io%iu6,1012) input%k_grid
         if (allocated(input%free_values)) then
            write(io%iu6,advance='no',fmt=1013)
            do i=1,structure%nions
               write(io%iu6,advance='no',fmt='(X,F9.3)') input%free_values(1,i)
            enddo
            write(io%iu6,fmt='(A3)') " au"
            write(io%iu6,advance='no',fmt=1014)
            do i=1,structure%nions
               write(io%iu6,advance='no',fmt='(X,F9.3)') input%free_values(2,i)
            enddo
            write(io%iu6,fmt='(A3)') " au"
            write(io%iu6,advance='no',fmt=1015)
            do i=1,structure%nions
               write(io%iu6,advance='no',fmt='(X,F9.3)') input%free_values(3,i)
            enddo
            write(io%iu6,fmt='(A3)') " au"
         endif
      endif

   end subroutine input_settings_libmbd
#endif

#ifdef LIBMBD
   type(libmbd_interaction) function vdw_forces_libmbd(io,t_info,ivdw,structure,input,relvol,comm) result (res)
      use base, only: in_struct
      use constant, only: RyToEv, AuToA
      use poscar, only: type_info
      use tutor, only: vtutor
      use mpimy
      type(libmbd_structure), intent(in) :: structure
      real(q) :: relvol(:)
      integer :: i, ivdw
      real(q), parameter :: HtrToEv=2*RyToEv
      type(in_struct), intent(in) :: io
      type(type_info), intent(in) :: t_info
      type(communic), intent(in) :: comm
      type(mbd_input_t) :: input
      type(mbd_calc_t) :: calc

      input%coords = matmul(structure%cell,structure%positions)/AuToA
      input%lattice_vectors = structure%cell/AuToA

      !input%comm = comm%mpi_comm
      call calc%init(input)
      call libmbd_error_handler(calc,'init',comm)

      call calc%update_vdw_params_from_ratios(relvol)
      call libmbd_error_handler(calc,'update_vdw_params_from_ratios',comm)

      call calc%evaluate_vdw_method(res%energy)
      call libmbd_error_handler(calc,'evaluate_vdw_method',comm)
      res%energy = res%energy*HtrToEv

      if (input%calculate_forces) then
         if (.not.allocated(res%forces)) allocate(res%forces(3,structure%nions))
         call calc%get_gradients(res%forces)
         call libmbd_error_handler(calc,'get_gradients',comm)
         res%forces = -res%forces*HtrToEv/AuToA

         call calc%get_lattice_stress(res%stress)
         call libmbd_error_handler(calc,'get_lattice_stress',comm)
         res%stress = matmul(res%stress,structure%cell_recip*AuToA)
         res%stress = matmul(res%stress,transpose(input%lattice_vectors))
         res%stress = -res%stress*HtrToEv
      endif

!      if (input%calculate_spectrum) then
!         if (.not.allocated(res%spectrum)) allocate(res%spectrum(3*structure%nions))
!         if (.not.allocated(res%modes)) allocate(res%modes(3*structure%nions,3*structure%nions))
!         call calc%get_spectrum_modes(res%spectrum,res%modes)
!      endif

  1001 FORMAT(/,"  libMBD for vdW energy calculation",&
              /,"  ---------------------------------------------------------------------------------------------")
  1002 FORMAT(/,"  Parameters of vdW forcefield:",&
              /,"            RELVOL",&
              /,"   -----------------------------")
  1003 FORMAT("   ",A2,5X,F12.8)
  1004 FORMAT("  Edisp (eV):",F11.5)

      if (io%iu6>=0) then
          write(io%iu6,1001)
          write(io%iu6,1002)
          do i=1,structure%nions
             write(io%iu6,1003) input%atom_types(i), relvol(i)
          enddo
          write(io%iu6,*)
          write(io%iu6,1004) res%energy
      endif

      call calc%destroy()

      contains

      !> call libmbd error handler
      !> analyzes the error code from libmbd and prints the error message if needed
      subroutine libmbd_error_handler(calc,routine,comm)
         type(mbd_calc_t) :: calc
         character(len=*), intent(in) :: routine
         type(communic), intent(in) :: comm
         ! local variables
         integer :: code
         character(200) :: origin, msg

         origin=''
         msg=''
         code=0
         call calc%get_exception(code,origin,msg)
         CALLMPI( M_sum_i(comm, code, 1) )
         if (code > 0) then
            call vtutor%error('libmdb '//trim(routine)//' error: '//trim(msg))
         end if
      end subroutine libmbd_error_handler
   end function vdw_forces_libmbd
#endif

end module libmbd
