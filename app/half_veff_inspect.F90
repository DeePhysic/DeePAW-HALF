program half_veff_inspect
  use half_kinds,only:dp
  use half_types,only:crystal_t,charge_grid_t,potcar_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar
  use half_potential,only:build_veff_lda
  implicit none
  type(crystal_t)::c
  type(charge_grid_t)::rho
  type(potcar_t),allocatable::p(:)
  real(dp),allocatable::veff(:)
  real(dp)::eh,exc
  character(len=1024)::chg,pot
  if(command_argument_count()/=2)then; write(*,'(A)')'usage: half-veff-inspect CHGCAR POTCAR'; stop 2; end if
  call get_command_argument(1,chg); call get_command_argument(2,pot)
  call read_chgcar(trim(chg),c,rho); call read_potcar(trim(pot),p); call build_veff_lda(rho,p,c,veff,eh,exc)
  write(*,'(A,ES24.16)')'hartree_eV=',eh; write(*,'(A,ES24.16)')'xc_eV=',exc
  write(*,'(A,ES24.16)')'veff_sum=',sum(veff); write(*,'(A,ES24.16)')'veff_norm2=',sum(veff*veff)
  write(*,'(A,ES24.16)')'veff_min=',minval(veff); write(*,'(A,ES24.16)')'veff_max=',maxval(veff)
end program
