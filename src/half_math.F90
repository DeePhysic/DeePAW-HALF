module half_math
  use half_kinds, only: dp
  implicit none
  private
  public :: determinant3, inverse3
contains
  pure function determinant3(a) result(det)
    real(dp), intent(in) :: a(3,3)
    real(dp) :: det
    det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function determinant3

  function inverse3(a) result(inv)
    real(dp), intent(in) :: a(3,3)
    real(dp) :: inv(3,3), det
    det = determinant3(a)
    if (abs(det) < 1.0e-18_dp) error stop "HALF: singular lattice"
    inv(1,1) =  (a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
    inv(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))/det
    inv(1,3) =  (a(1,2)*a(2,3)-a(1,3)*a(2,2))/det
    inv(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))/det
    inv(2,2) =  (a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
    inv(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))/det
    inv(3,1) =  (a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
    inv(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))/det
    inv(3,3) =  (a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
  end function inverse3
end module half_math
