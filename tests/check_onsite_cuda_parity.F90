program check_onsite_cuda_parity
  use half_kinds,only:dp
  use half_types,only:potcar_t,crystal_t,charge_grid_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar,validate_potcar_structure
  use half_paw,only:onsite_species_correction_t
  use half_uspp,only:augmentation_occupancy_t
  use half_onsite_density,only:atomic_onsite_occupation
  use half_onsite_functional,only:evaluate_onsite_corrections
  use half_cuda_onsite_functional,only:evaluate_onsite_corrections_cuda
  implicit none
  type(potcar_t),allocatable::potcars(:)
  type(crystal_t)::crystal
  type(charge_grid_t)::charge
  type(augmentation_occupancy_t),allocatable::occupancy(:)
  type(onsite_species_correction_t),allocatable::cpu(:),gpu(:)
  real(dp),allocatable::reference(:,:)
  real(dp)::cpu_dc,cpu_energy,gpu_dc,gpu_energy,dij_error,amplitude
  integer::it,iat,j,n,lmax
  character(len=1024)::charge_path,potcar_path,arg
  if(command_argument_count()<3)error stop 'usage: half-onsite-cuda-parity-check CHGCAR POTCAR LMAX [AMPLITUDE]'
  call get_command_argument(1,charge_path)
  call get_command_argument(2,potcar_path)
  call get_command_argument(3,arg);read(arg,*)lmax
  amplitude=0.03_dp
  if(command_argument_count()>=4)then
    call get_command_argument(4,arg);read(arg,*)amplitude
  end if
  call read_chgcar(trim(charge_path),crystal,charge)
  call read_potcar(trim(potcar_path),potcars)
  call validate_potcar_structure(potcars,crystal)
  allocate(occupancy(size(potcars)))
  do it=1,size(potcars)
    call atomic_onsite_occupation(potcars(it),reference)
    n=size(reference,1)
    allocate(occupancy(it)%matrix(crystal%counts(it),n,n))
    do iat=1,crystal%counts(it)
      occupancy(it)%matrix(iat,:,:)=reference
      do j=1,n
        occupancy(it)%matrix(iat,j,j)=occupancy(it)%matrix(iat,j,j)+ &
          amplitude*real(j,dp)/real(n,dp)
      end do
      if(n>=2)then
        occupancy(it)%matrix(iat,1,2)=occupancy(it)%matrix(iat,1,2)+0.002_dp
        occupancy(it)%matrix(iat,2,1)=occupancy(it)%matrix(iat,2,1)+0.002_dp
      end if
    end do
  end do
  call evaluate_onsite_corrections(potcars,crystal,occupancy,lmax,cpu,cpu_dc,cpu_energy)
  call evaluate_onsite_corrections_cuda(potcars,crystal,occupancy,lmax,gpu,gpu_dc,gpu_energy)
  dij_error=0.0_dp
  do it=1,size(potcars)
    dij_error=max(dij_error,maxval(abs(cpu(it)%dij-gpu(it)%dij)))
  end do
  write(*,'(A,ES16.8)')'energy_delta_eV=',gpu_energy-cpu_energy
  write(*,'(A,ES16.8)')'double_counting_delta_eV=',gpu_dc-cpu_dc
  write(*,'(A,ES16.8)')'max_Dij_delta_eV=',dij_error
  if(abs(gpu_energy-cpu_energy)>1.0e-5_dp.or.abs(gpu_dc-cpu_dc)>1.0e-5_dp.or. &
    dij_error>1.0e-4_dp)error stop 'CUDA onsite functional parity failed'
end program
