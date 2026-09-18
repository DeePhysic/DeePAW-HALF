program half_inspect_cpu
  use half_basis, only: build_plane_wave_basis
  use half_chgcar, only: read_chgcar
  use half_cpu, only: cpu_density_mean
  use half_kinds, only: dp
  use half_types, only: charge_grid_t, crystal_t, plane_wave_basis_t
  implicit none
  type(crystal_t) :: crystal
  type(charge_grid_t) :: charge
  type(plane_wave_basis_t) :: basis
  character(len=1024) :: path, argument
  integer :: argc, ios
  real(dp) :: encut
  argc = command_argument_count()
  if (argc == 1) then
    call get_command_argument(1,argument)
    if (trim(argument) == "--help" .or. trim(argument) == "-h") then
      write(*,'(A)') "Usage: half-inspect CHGCAR [ENCUT_eV]"
      stop
    end if
  end if
  if (argc < 1 .or. argc > 2) then
    write(*,'(A)') "Usage: half-inspect CHGCAR [ENCUT_eV]"
    error stop 2
  end if
  call get_command_argument(1,path)
  encut = 400.0_dp
  if (argc == 2) then
    call get_command_argument(2,argument)
    read(argument,*,iostat=ios) encut
    if (ios /= 0 .or. encut <= 0.0_dp) error stop "HALF: invalid ENCUT"
  end if
  call read_chgcar(trim(path), crystal, charge)
  call build_plane_wave_basis(crystal, charge%shape, encut, [0.0_dp,0.0_dp,0.0_dp], basis)
  write(*,'(A)') "HALF inspect CPU 0.2.0"
  write(*,'(A,A)') "system: ", trim(crystal%system_name)
  write(*,'(A,I0)') "ions: ", crystal%nions
  write(*,'(A,3(I0,1X))') "grid: ", charge%shape
  write(*,'(A,F0.10)') "volume_A3: ", crystal%volume
  write(*,'(A,F0.10)') "smooth_electrons: ", cpu_density_mean(charge%values)
  write(*,'(A,F0.3)') "encut_eV: ", encut
  write(*,'(A,I0)') "gamma_plane_waves: ", basis%npw
  write(*,'(A)') "backend: CPU Fortran"
end program half_inspect_cpu
