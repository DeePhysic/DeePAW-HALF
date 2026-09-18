module half_basis
  use half_constants, only: hsqdtm
  use half_kinds, only: dp, i32, i64
  use half_types, only: crystal_t, plane_wave_basis_t
  implicit none
  private
  public :: build_plane_wave_basis
contains
  subroutine build_plane_wave_basis(crystal, shape, encut, kpoint, basis)
    type(crystal_t), intent(in) :: crystal
    integer(i32), intent(in) :: shape(3)
    real(dp), intent(in) :: encut, kpoint(3)
    type(plane_wave_basis_t), intent(out) :: basis
    integer(i32) :: n1, n2, n3
    integer(i64) :: count, cursor, index0
    real(dp) :: q(3), kinetic

    if (encut <= 0.0_dp) error stop "HALF: ENCUT must be positive"
    count = 0
    do n1 = -shape(1)/2, (shape(1)-1)/2
      do n2 = -shape(2)/2, (shape(2)-1)/2
        do n3 = -shape(3)/2, (shape(3)-1)/2
          q = matmul([real(n1,dp)+kpoint(1), real(n2,dp)+kpoint(2), &
                      real(n3,dp)+kpoint(3)], crystal%reciprocal)
          if (hsqdtm*dot_product(q,q) < encut) count = count + 1
        end do
      end do
    end do

    basis%encut = encut
    basis%shape = shape
    basis%npw = count
    allocate(basis%gvectors(count,3), basis%qvectors(count,3), &
             basis%kinetic(count), basis%fft_index(count))
    cursor = 0
    do n1 = -shape(1)/2, (shape(1)-1)/2
      do n2 = -shape(2)/2, (shape(2)-1)/2
        do n3 = -shape(3)/2, (shape(3)-1)/2
          q = matmul([real(n1,dp)+kpoint(1), real(n2,dp)+kpoint(2), &
                      real(n3,dp)+kpoint(3)], crystal%reciprocal)
          kinetic = hsqdtm*dot_product(q,q)
          if (kinetic < encut) then
            cursor = cursor + 1
            basis%gvectors(cursor,:) = [n1,n2,n3]
            basis%qvectors(cursor,:) = q
            basis%kinetic(cursor) = kinetic
            index0 = int(modulo(n1,shape(1)),i64)*shape(2)*shape(3) &
                   + int(modulo(n2,shape(2)),i64)*shape(3) &
                   + int(modulo(n3,shape(3)),i64)
            basis%fft_index(cursor) = index0 + 1_i64
          end if
        end do
      end do
    end do
  end subroutine build_plane_wave_basis
end module half_basis
