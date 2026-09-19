program half_uspp_inspect
  use half_kinds,only:dp
  use half_types,only:crystal_t,charge_grid_t,plane_wave_basis_t,potcar_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar
  use half_basis,only:build_plane_wave_basis
  use half_paw,only:paw_species_t,build_paw_operators
  use half_potential,only:build_veff_pbe
  use half_uspp,only:build_uspp_dij_cpu
  implicit none
  type(crystal_t)::crystal
  type(charge_grid_t)::charge
  type(plane_wave_basis_t)::basis
  type(potcar_t),allocatable::potcars(:)
  type(paw_species_t),allocatable::paw(:)
  real(dp),allocatable::veff(:)
  real(dp)::encut,eh,exc
  character(len=1024)::charge_path,potcar_path,arg
  integer::it,iat
  if(command_argument_count()<2)then
    write(*,'(A)')'Usage: half-uspp-inspect CHGCAR POTCAR [ENCUT]';stop 2
  end if
  call get_command_argument(1,charge_path);call get_command_argument(2,potcar_path);encut=400.0_dp
  if(command_argument_count()>=3)then;call get_command_argument(3,arg);read(arg,*)encut;end if
  call read_chgcar(trim(charge_path),crystal,charge);call read_potcar(trim(potcar_path),potcars)
  call build_plane_wave_basis(crystal,charge%shape,encut,[0.0_dp,0.0_dp,0.0_dp],basis)
  call build_paw_operators(potcars,crystal,basis,paw);call build_veff_pbe(charge,potcars,crystal,veff,eh,exc)
  call build_uspp_dij_cpu(veff,charge%shape,potcars,crystal,paw)
  do it=1,size(paw);do iat=1,paw(it)%natoms
    write(*,'(A,I0,A,I0,A,3(ES24.16,1X))')'species=',it,' atom=',iat,' min/max/sum=', &
      minval(paw(it)%dij_atom(iat,:,:)),maxval(paw(it)%dij_atom(iat,:,:)),sum(paw(it)%dij_atom(iat,:,:))
  end do;end do
end program
