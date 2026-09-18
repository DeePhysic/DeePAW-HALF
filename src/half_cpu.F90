module half_cpu
  use half_kinds, only: dp, i64
  use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: cpu_density_mean, cpu_density_benchmark, cpu_max_threads
contains
  pure function cpu_density_mean(values) result(mean_value)
    real(dp), intent(in) :: values(:)
    real(dp) :: mean_value
    mean_value = sum(values)/real(size(values,kind=i64),dp)
  end function cpu_density_mean

  subroutine cpu_density_benchmark(values, iterations, elapsed_seconds, checksum)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: iterations
    real(dp), intent(out) :: elapsed_seconds, checksum
    real(dp), allocatable :: output(:)
    integer(i64) :: tick_start, tick_stop, tick_rate, i
    integer :: iteration
    real(dp), parameter :: alpha = 0.999999_dp, beta = 0.000001_dp
    allocate(output(size(values)))
    output = 0.0_dp
    call system_clock(tick_start, tick_rate)
    !$omp parallel default(shared) private(iteration,i)
    do iteration = 1, iterations
      !$omp do schedule(static)
      do i = 1_i64, size(values,kind=i64)
        output(i) = alpha*values(i) + beta*output(i)
      end do
      !$omp end do
    end do
    !$omp end parallel
    call system_clock(tick_stop)
    elapsed_seconds = real(tick_stop-tick_start,dp)/real(tick_rate,dp)
    checksum = sum(output)
    deallocate(output)
  end subroutine cpu_density_benchmark

  function cpu_max_threads() result(count)
    integer :: count
    count = omp_get_max_threads()
  end function cpu_max_threads
end module half_cpu
