module half_constants
  use half_kinds, only: dp
  implicit none
  private
  public :: pi, twopi, autoa, rytoev, hatoev, felect, edeps, hsqdtm
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
  real(dp), parameter :: twopi = 2.0_dp*pi
  real(dp), parameter :: autoa = 0.529177249_dp
  real(dp), parameter :: rytoev = 13.605826_dp
  real(dp), parameter :: hatoev = 2.0_dp*rytoev
  real(dp), parameter :: felect = 2.0_dp*rytoev*autoa
  real(dp), parameter :: edeps = 4.0_dp*pi*felect
  ! hbar^2 / (2 m_e), in eV Angstrom^2. Uses HAPPY/VASP constants.
  real(dp), parameter :: hsqdtm = rytoev*autoa*autoa
end module half_constants
