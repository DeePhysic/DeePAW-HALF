program half_slice_benchmark
  use half_basis, only: build_plane_wave_basis
  use half_chgcar, only: read_chgcar
  use half_kinds, only: dp, i64
  use half_types, only: charge_grid_t, crystal_t, plane_wave_basis_t
  implicit none
  character(len=1024) :: path, argument
  integer :: argc, ios, iteration, iterations
  integer(i64) :: clock_begin, clock_end, clock_rate, final_npw
  real(dp) :: elapsed, final_electrons, final_kinetic_sum

  argc = command_argument_count()
  if (argc == 1) then
    call get_command_argument(1, argument)
    if (trim(argument) == "--help" .or. trim(argument) == "-h") then
      call usage()
      stop
    end if
  end if
  if (argc < 1 .or. argc > 2) then
    call usage()
    error stop 2
  end if
  call get_command_argument(1, path)
  iterations = 20
  if (argc == 2) then
    call get_command_argument(2, argument)
    read(argument, *, iostat=ios) iterations
    if (ios /= 0 .or. iterations < 1) error stop "HALF: invalid iteration count"
  end if

  call run_once(path, final_npw, final_electrons, final_kinetic_sum)
  call system_clock(clock_begin, clock_rate)
  do iteration = 1, iterations
    call run_once(path, final_npw, final_electrons, final_kinetic_sum)
  end do
  call system_clock(clock_end)
  elapsed = real(clock_end-clock_begin,dp)/real(clock_rate,dp)

  write(*,'(A)') "HALF Python-parity slice benchmark"
  write(*,'(A,I0)') "iterations: ", iterations
  write(*,'(A,F0.9)') "total_seconds: ", elapsed
  write(*,'(A,F0.9)') "seconds_per_iteration: ", elapsed/real(iterations,dp)
  write(*,'(A,I0)') "gamma_plane_waves: ", final_npw
  write(*,'(A,ES24.16)') "smooth_electrons: ", final_electrons
  write(*,'(A,ES24.16)') "kinetic_checksum_eV: ", final_kinetic_sum

contains
  subroutine run_once(input_path, npw, electrons, kinetic_sum)
    character(len=*), intent(in) :: input_path
    integer(i64), intent(out) :: npw
    real(dp), intent(out) :: electrons, kinetic_sum
    type(crystal_t) :: crystal
    type(charge_grid_t) :: charge
    type(plane_wave_basis_t) :: basis
    call read_chgcar(trim(input_path), crystal, charge)
    call build_plane_wave_basis(crystal, charge%shape, 400.0_dp, [0.0_dp,0.0_dp,0.0_dp], basis)
    npw = basis%npw
    electrons = charge%electron_count()
    kinetic_sum = sum(basis%kinetic)
  end subroutine run_once

  subroutine usage()
    write(*,'(A)') "Usage: half-slice-benchmark CHGCAR [ITERATIONS]"
  end subroutine usage
end program half_slice_benchmark
