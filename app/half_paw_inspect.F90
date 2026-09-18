program half_paw_inspect
  use half_kinds,only:dp
  use half_types,only:crystal_t,charge_grid_t,plane_wave_basis_t,potcar_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar
  use half_basis,only:build_plane_wave_basis
  use half_paw,only:paw_species_t,build_paw_operators
  implicit none
  type(crystal_t)::c
  type(charge_grid_t)::rho
  type(plane_wave_basis_t)::b
  type(potcar_t),allocatable::p(:)
  type(paw_species_t),allocatable::paw(:)
  character(len=1024)::chg,pot,arg
  real(dp)::encut
  if(command_argument_count()/=3)then
    write(*,'(A)')'usage: half-paw-inspect CHGCAR POTCAR ENCUT'; stop 2
  end if
  call get_command_argument(1,chg); call get_command_argument(2,pot); call get_command_argument(3,arg); read(arg,*)encut
  call read_chgcar(trim(chg),c,rho); call read_potcar(trim(pot),p)
  call build_plane_wave_basis(c,rho%shape,encut,[0.0_dp,0.0_dp,0.0_dp],b)
  call build_paw_operators(p,c,b,paw)
  write(*,'(A,I0)')'npw=',b%npw; write(*,'(A,I0)')'nlm=',paw(1)%nlm
  write(*,'(A,ES24.16)')'projector_real_sum=',sum(real(paw(1)%projectors,dp))
  write(*,'(A,ES24.16)')'projector_imag_sum=',sum(aimag(paw(1)%projectors))
  write(*,'(A,ES24.16)')'projector_norm2=',sum(abs(paw(1)%projectors)**2)
  write(*,'(A,ES24.16)')'dij_sum=',sum(paw(1)%dij)
  write(*,'(A,ES24.16)')'qij_sum=',sum(paw(1)%qij)
end program
