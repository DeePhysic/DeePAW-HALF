program half_inspect
  use half_basis, only: build_plane_wave_basis
  use half_chgcar, only: read_chgcar
  use half_cuda, only: cuda_density_mean, cuda_runtime_summary
  use half_kinds, only: dp, i64
  use half_types, only: charge_grid_t, crystal_t, plane_wave_basis_t
  implicit none
  type(crystal_t) :: crystal
  type(charge_grid_t) :: charge
  type(plane_wave_basis_t) :: basis
  character(len=1024) :: path, argument, basis_output
  character(len=256) :: gpu_name
  integer :: argc, major, minor, ios
  integer(i64) :: memory_bytes
  real(dp) :: encut, cpu_count, gpu_count

  argc = command_argument_count()
  if (argc == 1) then
    call get_command_argument(1, argument)
    if (trim(argument) == "--help" .or. trim(argument) == "-h") then
      call usage()
      stop
    end if
  end if
  if (argc < 1 .or. argc > 3) then
    call usage()
    error stop 2
  end if
  call get_command_argument(1, path)
  encut = 400.0_dp
  if (argc >= 2) then
    call get_command_argument(2, argument)
    read(argument, *, iostat=ios) encut
    if (ios /= 0 .or. encut <= 0.0_dp) error stop "HALF: invalid ENCUT"
  end if
  if (argc >= 3) then
    call get_command_argument(3, basis_output)
  else
    basis_output = ""
  end if

  call read_chgcar(trim(path), crystal, charge)
  call build_plane_wave_basis(crystal, charge%shape, encut, [0.0_dp,0.0_dp,0.0_dp], basis)
  cpu_count = charge%electron_count()
  gpu_count = cuda_density_mean(charge%values)
  if (abs(cpu_count-gpu_count) > 5.0e-10_dp) error stop "HALF: CPU/GPU density mismatch"
  call cuda_runtime_summary(gpu_name, major, minor, memory_bytes)
  if (len_trim(basis_output) > 0) call write_basis_csv(trim(basis_output), basis)

  write(*, '(A)') "HALF inspect 0.2.0"
  write(*, '(A,A)') "system: ", trim(crystal%system_name)
  write(*, '(A,I0)') "ions: ", crystal%nions
  write(*, '(A,3(I0,1X))') "grid: ", charge%shape
  write(*, '(A,F0.10)') "volume_A3: ", crystal%volume
  write(*, '(A,F0.10)') "smooth_electrons: ", gpu_count
  write(*, '(A,F0.3)') "encut_eV: ", encut
  write(*, '(A,I0)') "gamma_plane_waves: ", basis%npw
  write(*, '(A,A)') "gpu: ", trim(gpu_name)
  write(*, '(A,I0,A,I0)') "compute_capability: ", major, ".", minor
  write(*, '(A,F0.3)') "gpu_memory_GiB: ", real(memory_bytes,dp)/(1024.0_dp**3)
contains
  subroutine usage()
    write(*, '(A)') "Usage: half-inspect CHGCAR [ENCUT_eV [basis.csv]]"
  end subroutine usage

  subroutine write_basis_csv(filename, pw)
    character(len=*), intent(in) :: filename
    type(plane_wave_basis_t), intent(in) :: pw
    integer :: unit, status
    integer(i64) :: i
    open(newunit=unit, file=filename, status="replace", action="write", iostat=status)
    if (status /= 0) error stop "HALF: cannot open basis CSV output"
    write(unit, '(A)') "n1,n2,n3,fft_index_1based,qx,qy,qz,kinetic_eV"
    do i = 1, pw%npw
      write(unit, '(I0,A,I0,A,I0,A,I0,4(A,ES25.16E3))') &
        pw%gvectors(i,1), ',', pw%gvectors(i,2), ',', pw%gvectors(i,3), ',', &
        pw%fft_index(i), ',', pw%qvectors(i,1), ',', pw%qvectors(i,2), ',', &
        pw%qvectors(i,3), ',', pw%kinetic(i)
    end do
    close(unit)
  end subroutine write_basis_csv
end program half_inspect
