program half_force_check
  use half_kinds,only:dp
  use half_types,only:crystal_t,charge_grid_t,potcar_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar,validate_potcar_structure
  use half_potential,only:build_veff_pbe
  use half_local_forces,only:nlcc_forces
  use half_math,only:inverse3
  implicit none
  type(crystal_t)::crystal,shifted
  type(charge_grid_t)::charge
  type(potcar_t),allocatable::potcars(:)
  real(dp),allocatable::veff(:),vxc(:),analytic(:,:),numeric(:,:)
  real(dp)::eh,exc,exv,ep,em,h,invlat(3,3),delta(3),dummy_h,dummy_v
  character(len=1024)::chgcar,potcar
  integer::iat,alpha
  if(command_argument_count()/=2)error stop 'usage: half-force-check CHGCAR POTCAR'
  call get_command_argument(1,chgcar);call get_command_argument(2,potcar)
  call read_chgcar(trim(chgcar),crystal,charge);call read_potcar(trim(potcar),potcars)
  call validate_potcar_structure(potcars,crystal)
  call build_veff_pbe(charge,potcars,crystal,veff,eh,exc,exv,vxc)
  allocate(analytic(crystal%nions,3),numeric(crystal%nions,3));call nlcc_forces(vxc,potcars,crystal,charge%shape,analytic)
  invlat=inverse3(crystal%lattice);h=1.0e-5_dp;numeric=0.0_dp
  do iat=1,crystal%nions;do alpha=1,3
    delta=0.0_dp;delta(alpha)=h;shifted=crystal
    shifted%positions(iat,:)=shifted%positions(iat,:)+matmul(delta,invlat)
    call build_veff_pbe(charge,potcars,shifted,veff,dummy_h,ep,dummy_v)
    shifted=crystal;shifted%positions(iat,:)=shifted%positions(iat,:)-matmul(delta,invlat)
    call build_veff_pbe(charge,potcars,shifted,veff,dummy_h,em,dummy_v)
    numeric(iat,alpha)=-(ep-em)/(2*h)
  end do;end do
  write(*,'(a,es24.16)')'nlcc_force_max_error_eV_per_A=',maxval(abs(analytic-numeric))
  write(*,'(a,es24.16)')'nlcc_force_max_analytic_eV_per_A=',maxval(abs(analytic))
  if(maxval(abs(analytic-numeric))>2.0e-5_dp)error stop 'HALF: analytic NLCC force regression'
end program half_force_check
