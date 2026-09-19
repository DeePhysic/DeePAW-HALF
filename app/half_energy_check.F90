program half_energy_check
  use half_kinds,only:dp
  use half_types,only:crystal_t
  use half_energy,only:compute_occupations,ewald_energy
  implicit none
  type(crystal_t)::crystal
  real(dp)::eig(2,3),weight(2),mu,band,entropy,e1,e2
  real(dp),allocatable::occ(:,:)

  eig=0.0_dp;eig(:,1:2)=reshape([0.0_dp,1.0_dp,2.0_dp,3.0_dp],[2,2])
  weight=[0.5_dp,0.5_dp]
  call compute_occupations(eig(:,1:2),weight,2.0_dp,0.0_dp,occ,mu,band,entropy)
  if(maxval(abs(occ-reshape([2.0_dp,2.0_dp,0.0_dp,0.0_dp],[2,2])))>1e-14_dp) &
    error stop 'HALF: zero-temperature occupation regression'
  if(abs(mu-1.0_dp)>1e-14_dp.or.abs(band-1.0_dp)>1e-14_dp.or.abs(entropy)>1e-14_dp) &
    error stop 'HALF: zero-temperature energy regression'

  eig=reshape([-1.0_dp,-0.8_dp,0.5_dp,0.7_dp,2.0_dp,2.2_dp],[2,3])
  weight=[0.25_dp,0.75_dp]
  call compute_occupations(eig,weight,2.0_dp,0.1_dp,occ,mu,band,entropy)
  if(abs(sum(spread(weight,2,3)*occ)-2.0_dp)>1e-12_dp.or.entropy>=0.0_dp) &
    error stop 'HALF: finite-temperature occupation regression'

  crystal%nions=2;crystal%ntypes=1;crystal%lattice=0.0_dp
  crystal%lattice(1,1)=4.0_dp;crystal%lattice(2,2)=4.0_dp;crystal%lattice(3,3)=4.0_dp
  allocate(crystal%counts(1),crystal%species(1),crystal%positions(2,3))
  crystal%counts=2;crystal%species='X'
  crystal%positions=reshape([0.0_dp,0.5_dp,0.0_dp,0.5_dp,0.0_dp,0.5_dp],[2,3])
  call crystal%update_geometry()
  call ewald_energy(crystal,[1.0_dp,1.0_dp],e1,eta=0.35_dp,tolerance=1e-12_dp)
  call ewald_energy(crystal,[1.0_dp,1.0_dp],e2,eta=0.80_dp,tolerance=1e-12_dp)
  if(abs(e1-e2)>1e-10_dp)error stop 'HALF: Ewald split-parameter regression'
  write(*,'(a,es24.16)')'ewald_ev=',e1
  write(*,'(a,es24.16)')'finite_temperature_mu_ev=',mu
  write(*,'(a)')'HALF occupation and Ewald checks passed'
end program
