module half_types
  use half_constants, only: twopi
  use half_kinds, only: dp, i32, i64
  use half_math, only: determinant3, inverse3
  implicit none
  private
  public :: crystal_t, charge_grid_t, plane_wave_basis_t, potcar_t

  type :: crystal_t
    character(len=:), allocatable :: system_name
    integer(i32) :: ntypes = 0
    integer(i32) :: nions = 0
    character(len=16), allocatable :: species(:)
    integer(i32), allocatable :: counts(:)
    real(dp) :: lattice(3,3) = 0.0_dp
    real(dp) :: reciprocal(3,3) = 0.0_dp
    real(dp) :: volume = 0.0_dp
    real(dp), allocatable :: positions(:,:)
  contains
    procedure :: update_geometry => crystal_update_geometry
  end type crystal_t

  type :: charge_grid_t
    integer(i32) :: shape(3) = 0
    real(dp), allocatable :: values(:)
  contains
    procedure :: size => charge_size
    procedure :: electron_count => charge_electron_count
    procedure :: normalize_electron_count => charge_normalize_electron_count
  end type charge_grid_t

  type :: plane_wave_basis_t
    real(dp) :: encut = 0.0_dp
    integer(i32) :: shape(3) = 0
    integer(i64) :: npw = 0
    integer(i32), allocatable :: gvectors(:,:)
    real(dp), allocatable :: qvectors(:,:)
    real(dp), allocatable :: kinetic(:)
    integer(i64), allocatable :: fft_index(:)
  end type plane_wave_basis_t

  type :: potcar_t
    character(len=:), allocatable :: header
    character(len=16) :: element = ""
    character(len=8) :: lexch = "PE"
    real(dp) :: zval = 0.0_dp, pomass = 0.0_dp, enmax = 0.0_dp
    real(dp) :: eaug = 0.0_dp, eatom = 0.0_dp, dexccore = 0.0_dp
    real(dp) :: psp_gmax = 0.0_dp, pspnl_gmax = 0.0_dp, paw_rmax = 0.0_dp
    logical :: has_core = .false.
    integer(i32) :: channels = 0, nmax = 0
    integer(i32), allocatable :: lps(:)
    real(dp), allocatable :: psp_local(:), psp_g_grid(:), pspcor(:), psprho(:)
    real(dp), allocatable :: dion(:,:), pspnl(:,:), psprnl(:,:), pspnl_rmax(:)
    real(dp), allocatable :: rgrid(:), potae(:), potps(:), potpsc(:)
    real(dp), allocatable :: rhoae(:), rhops(:), wae(:,:), wps(:,:)
    real(dp), allocatable :: qpaw(:,:), qpaw_l(:,:,:), qato(:,:)
    real(dp) :: qpaw_rebuild_max_delta = 0.0_dp
  end type potcar_t

contains
  subroutine crystal_update_geometry(self)
    class(crystal_t), intent(inout) :: self
    self%volume = abs(determinant3(self%lattice))
    self%reciprocal = twopi*transpose(inverse3(self%lattice))
  end subroutine crystal_update_geometry

  pure function charge_size(self) result(n)
    class(charge_grid_t), intent(in) :: self
    integer(i64) :: n
    n = int(self%shape(1),i64)*int(self%shape(2),i64)*int(self%shape(3),i64)
  end function charge_size

  pure function charge_electron_count(self) result(nelect)
    class(charge_grid_t), intent(in) :: self
    real(dp) :: nelect
    if (.not. allocated(self%values)) then
      nelect = 0.0_dp
    else
      nelect = sum(self%values)/real(size(self%values), dp)
    end if
  end function charge_electron_count

  subroutine charge_normalize_electron_count(self,target_electrons,original_electrons,scale)
    class(charge_grid_t),intent(inout)::self
    real(dp),intent(in)::target_electrons
    real(dp),intent(out),optional::original_electrons,scale
    real(dp)::current,factor
    if(.not.allocated(self%values))error stop 'HALF: cannot normalize an unallocated charge grid'
    if(target_electrons<=0.0_dp.or.target_electrons/=target_electrons) &
      error stop 'HALF: density normalization target must be positive and finite'
    current=self%electron_count()
    if(abs(current)<=tiny(1.0_dp).or.current/=current) &
      error stop 'HALF: cannot normalize a charge grid with zero or non-finite integral'
    factor=target_electrons/current;self%values=self%values*factor
    if(present(original_electrons))original_electrons=current
    if(present(scale))scale=factor
  end subroutine charge_normalize_electron_count
end module half_types
