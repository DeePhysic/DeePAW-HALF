program half_benchmark
  use half_chgcar, only: read_chgcar
  use half_cpu, only: cpu_density_benchmark, cpu_max_threads
  use half_cuda, only: cuda_density_benchmark, cuda_runtime_summary
  use half_kinds, only: dp, i64
  use half_types, only: charge_grid_t, crystal_t
  implicit none
  type(crystal_t) :: crystal
  type(charge_grid_t) :: charge
  character(len=1024) :: path, argument
  character(len=256) :: gpu_name
  integer :: argc, iterations, replicas, ios, major, minor, replica
  integer(i64) :: memory_bytes
  real(dp) :: cpu_seconds, gpu_seconds, cpu_checksum, gpu_checksum
  real(dp), allocatable :: workload(:)

  argc = command_argument_count()
  if (argc < 1 .or. argc > 3) then
    write(*,'(A)') "Usage: half-benchmark CHGCAR [iterations [grid_replicas]]"
    error stop 2
  end if
  call get_command_argument(1, path)
  iterations = 2000
  replicas = 1
  if (argc >= 2) then
    call get_command_argument(2, argument)
    read(argument,*,iostat=ios) iterations
    if (ios /= 0 .or. iterations <= 0) error stop "HALF: invalid iteration count"
  end if
  if (argc == 3) then
    call get_command_argument(3, argument)
    read(argument,*,iostat=ios) replicas
    if (ios /= 0 .or. replicas <= 0) error stop "HALF: invalid grid replica count"
  end if
  call read_chgcar(trim(path), crystal, charge)
  allocate(workload(size(charge%values)*replicas))
  do replica = 1, replicas
    workload((replica-1)*size(charge%values)+1:replica*size(charge%values)) = charge%values
  end do
  call cpu_density_benchmark(workload, iterations, cpu_seconds, cpu_checksum)
  call cuda_density_benchmark(workload, iterations, gpu_seconds, gpu_checksum)
  call cuda_runtime_summary(gpu_name, major, minor, memory_bytes)
  if (abs(cpu_checksum-gpu_checksum) > 1.0e-8_dp*abs(cpu_checksum)) &
    error stop "HALF: CPU/CUDA benchmark checksum mismatch"

  write(*,'(A)') "HALF density backend benchmark"
  write(*,'(A,I0)') "grid_values: ", size(workload,kind=i64)
  write(*,'(A,I0)') "grid_replicas: ", replicas
  write(*,'(A,I0)') "iterations: ", iterations
  write(*,'(A,I0)') "cpu_openmp_threads: ", cpu_max_threads()
  write(*,'(A,A)') "gpu: ", trim(gpu_name)
  write(*,'(A,I0,A,I0)') "compute_capability: ", major, ".", minor
  write(*,'(A,F0.9)') "cpu_seconds: ", cpu_seconds
  write(*,'(A,F0.9)') "cuda_kernel_seconds: ", gpu_seconds
  write(*,'(A,F0.4)') "kernel_speedup: ", cpu_seconds/gpu_seconds
  write(*,'(A,ES24.16E3)') "cpu_checksum: ", cpu_checksum
  write(*,'(A,ES24.16E3)') "cuda_checksum: ", gpu_checksum
end program half_benchmark
