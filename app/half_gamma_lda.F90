program half_gamma_lda
  use half_kinds,only:dp
  use half_types,only:crystal_t,charge_grid_t,plane_wave_basis_t,potcar_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar
  use half_basis,only:build_plane_wave_basis
  use half_paw,only:paw_species_t,build_paw_operators
  use half_potential,only:build_veff_lda,build_veff_pbe
  use half_dense_solver,only:solve_dense_gamma
  implicit none
  type(crystal_t)::c
  type(charge_grid_t)::rho
  type(plane_wave_basis_t)::basis
  type(potcar_t),allocatable::p(:)
  type(paw_species_t),allocatable::paw(:)
  real(dp),allocatable::veff(:),eig(:)
  real(dp)::eh,exc,smin,smax,encut
  integer::i,nbands
  character(len=1024)::chg,pot,arg,xc
  if(command_argument_count()/=5)then; write(*,'(A)')'usage: half-gamma CHGCAR POTCAR ENCUT NBANDS lda|pbe'; stop 2; end if
  call get_command_argument(1,chg); call get_command_argument(2,pot)
  call get_command_argument(3,arg); read(arg,*)encut; call get_command_argument(4,arg); read(arg,*)nbands
  call get_command_argument(5,xc)
  call read_chgcar(trim(chg),c,rho); call read_potcar(trim(pot),p)
  call build_plane_wave_basis(c,rho%shape,encut,[0.0_dp,0.0_dp,0.0_dp],basis)
  call build_paw_operators(p,c,basis,paw)
  if(trim(xc)=='pbe')then; call build_veff_pbe(rho,p,c,veff,eh,exc); else; call build_veff_lda(rho,p,c,veff,eh,exc); end if
  call solve_dense_gamma(veff,basis,paw,eig,smin,smax)
  write(*,'(A,I0)')'npw=',basis%npw; write(*,'(A,ES24.16)')'overlap_min=',smin
  write(*,'(A,ES24.16)')'overlap_max=',smax
  do i=1,min(nbands,size(eig)); write(*,'(A,I0,A,ES24.16)')'eigenvalue_',i,'=',eig(i); end do
end program
