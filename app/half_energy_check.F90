program half_energy_check
  use half_kinds,only:dp
  use half_types,only:crystal_t,plane_wave_basis_t
  use half_math,only:inverse3
  use half_energy,only:compute_occupations,ewald_energy,ewald_forces
  use half_paw,only:paw_species_t
  use half_forces,only:add_nonlocal_paw_forces
  use half_potcar_interp,only:vasp_local_spline_second,vasp_local_spline_eval,vasp_four_point_eval
  use half_kpoints,only:kpoint_set_t,symmetrize_scalar_grid,symmetrize_atomic_vectors
  implicit none
  type(crystal_t)::crystal,shifted
  type(plane_wave_basis_t)::basis
  type(paw_species_t),allocatable::paw(:)
  type(kpoint_set_t)::symset
  real(dp)::eig(2,3),weight(2),mu,band,entropy,e1,e2,ep,em,h,phase,table(8),m2table(8),dx
  real(dp)::analytic(2,3),numeric(2,3),inverse_lattice(3,3),delta(3)
  real(dp)::nlforce(1,3),nlnumeric(3),eval(1),focc(1),r0(3),eplus,eminus
  real(dp)::symgrid(8),symvectors(2,3),expected_vectors(2,3)
  complex(dp)::waves(3,1),projector0(2,3),cc(2)
  integer::i,j,ig
  real(dp),allocatable::occ(:,:)

  eig=0.0_dp;eig(:,1:2)=reshape([0.0_dp,1.0_dp,2.0_dp,3.0_dp],[2,2])
  weight=[0.5_dp,0.5_dp]
  call compute_occupations(eig(:,1:2),weight,2.0_dp,0.0_dp,occ,mu,band,entropy)
  if(maxval(abs(occ-reshape([2.0_dp,2.0_dp,0.0_dp,0.0_dp],[2,2])))>1e-14_dp) &
    error stop 'HALF: zero-temperature occupation regression'
  if(abs(mu-1.0_dp)>1e-14_dp.or.abs(band-1.0_dp)>1e-14_dp.or.abs(entropy)>1e-14_dp) &
    error stop 'HALF: zero-temperature energy regression'

  do i=1,8;table(i)=real(i-1,dp)**3;end do
  if(abs(vasp_four_point_eval(table,8.0_dp,2.3_dp)-2.3_dp**3)>1e-13_dp) &
    error stop 'HALF: VASP four-point POTCAR interpolation regression'
  do i=1,8;table(i)=real(i-1,dp)**2;end do
  call vasp_local_spline_second(table,8.0_dp,m2table);dx=1.0e-7_dp
  if(abs(m2table(8))>1e-14_dp.or.abs((vasp_local_spline_eval(table,m2table,8.0_dp,dx)- &
      vasp_local_spline_eval(table,m2table,8.0_dp,0.0_dp))/dx)>1e-6_dp) &
    error stop 'HALF: VASP local-potential spline boundary regression'

  symset%nsym=2
  allocate(symset%rotations(3,3,2),symset%translations(3,2),symset%cart_rotations(3,3,2),symset%atom_map(2,2))
  symset%rotations=0;symset%translations=0.0_dp;symset%cart_rotations=0.0_dp
  do i=1,3;symset%rotations(i,i,:)=1;symset%cart_rotations(i,i,1)=1.0_dp;end do
  symset%translations(:,2)=[0.5_dp,0.0_dp,0.0_dp]
  symset%cart_rotations(1,1,2)=-1.0_dp;symset%cart_rotations(2,2,2)=-1.0_dp
  symset%cart_rotations(3,3,2)=1.0_dp;symset%atom_map(:,1)=[1,2];symset%atom_map(:,2)=[2,1]
  symgrid=[(real(i,dp),i=1,8)];call symmetrize_scalar_grid(symset,[2,2,2],symgrid)
  if(maxval(abs(symgrid-[1.5_dp,1.5_dp,3.5_dp,3.5_dp,5.5_dp,5.5_dp,7.5_dp,7.5_dp]))>1e-14_dp) &
    error stop 'HALF: scalar-grid space-group expansion regression'
  symvectors=reshape([1.0_dp,4.0_dp,2.0_dp,5.0_dp,3.0_dp,6.0_dp],[2,3])
  expected_vectors(1,:)=[-1.5_dp,-1.5_dp,4.5_dp];expected_vectors(2,:)=[1.5_dp,1.5_dp,4.5_dp]
  call symmetrize_atomic_vectors(symset,symvectors)
  if(maxval(abs(symvectors-expected_vectors))>1e-14_dp) &
    error stop 'HALF: atomic-vector space-group expansion regression'

  eig=reshape([-1.0_dp,-0.8_dp,0.5_dp,0.7_dp,2.0_dp,2.2_dp],[2,3])
  weight=[0.25_dp,0.75_dp]
  call compute_occupations(eig,weight,2.0_dp,0.1_dp,occ,mu,band,entropy)
  if(abs(sum(spread(weight,2,3)*occ)-2.0_dp)>1e-12_dp.or.entropy>=0.0_dp) &
    error stop 'HALF: finite-temperature occupation regression'

  crystal%nions=2;crystal%ntypes=1;crystal%lattice=0.0_dp
  crystal%lattice(1,1)=4.0_dp;crystal%lattice(2,2)=4.0_dp;crystal%lattice(3,3)=4.0_dp
  allocate(crystal%counts(1),crystal%species(1),crystal%positions(2,3))
  crystal%counts=2;crystal%species='X'
  crystal%positions(1,:)=[0.07_dp,0.11_dp,0.19_dp]
  crystal%positions(2,:)=[0.31_dp,0.46_dp,0.73_dp]
  call crystal%update_geometry()
  call ewald_energy(crystal,[1.0_dp,1.0_dp],e1,eta=0.35_dp,tolerance=1e-12_dp)
  call ewald_energy(crystal,[1.0_dp,1.0_dp],e2,eta=0.80_dp,tolerance=1e-12_dp)
  if(abs(e1-e2)>1e-10_dp)error stop 'HALF: Ewald split-parameter regression'
  call ewald_forces(crystal,[1.0_dp,1.0_dp],analytic,eta=0.35_dp,tolerance=1e-12_dp)
  inverse_lattice=inverse3(crystal%lattice);h=1e-5_dp;numeric=0.0_dp
  do i=1,2;do j=1,3
    delta=0.0_dp;delta(j)=h;shifted=crystal
    shifted%positions(i,:)=shifted%positions(i,:)+matmul(delta,inverse_lattice)
    call ewald_energy(shifted,[1.0_dp,1.0_dp],ep,eta=0.35_dp,tolerance=1e-12_dp)
    shifted=crystal;shifted%positions(i,:)=shifted%positions(i,:)-matmul(delta,inverse_lattice)
    call ewald_energy(shifted,[1.0_dp,1.0_dp],em,eta=0.35_dp,tolerance=1e-12_dp)
    numeric(i,j)=-(ep-em)/(2*h)
  end do;end do
  write(*,'(a,es24.16)')'ewald_force_max_error=',maxval(abs(analytic-numeric))
  if(maxval(abs(analytic-numeric))>2e-7_dp)error stop 'HALF: analytic Ewald force regression'
  basis%npw=3;allocate(basis%qvectors(3,3));basis%qvectors=reshape([ &
    0.2_dp,0.7_dp,-0.4_dp, -0.3_dp,0.5_dp,0.8_dp, 0.6_dp,-0.2_dp,0.9_dp],[3,3])
  allocate(paw(1));paw(1)%natoms=1;paw(1)%nlm=2
  allocate(paw(1)%projectors(1,2,3),paw(1)%dij(2,2),paw(1)%qij(2,2),paw(1)%dij_atom(1,2,2))
  projector0=reshape([cmplx(0.4_dp,0.1_dp,dp),cmplx(-0.2_dp,0.3_dp,dp), &
    cmplx(0.1_dp,-0.2_dp,dp),cmplx(0.5_dp,0.4_dp,dp),cmplx(-0.3_dp,0.2_dp,dp), &
    cmplx(0.2_dp,0.6_dp,dp)],[2,3])
  waves(:,1)=[cmplx(0.3_dp,-0.1_dp,dp),cmplx(-0.2_dp,0.4_dp,dp),cmplx(0.5_dp,0.2_dp,dp)]
  paw(1)%dij_atom(1,:,:)=reshape([1.2_dp,0.15_dp,0.15_dp,0.8_dp],[2,2])
  paw(1)%dij=paw(1)%dij_atom(1,:,:);paw(1)%qij=reshape([0.1_dp,0.02_dp,0.02_dp,0.07_dp],[2,2])
  eval=[0.37_dp];focc=[1.6_dp];r0=[0.13_dp,-0.21_dp,0.08_dp]
  do ig=1,3
    phase=dot_product(basis%qvectors(ig,:),r0)
    paw(1)%projectors(1,:,ig)=projector0(:,ig)*cmplx(cos(phase),-sin(phase),dp)
  end do
  nlforce=0.0_dp;call add_nonlocal_paw_forces(basis,paw,eval,focc,0.4_dp,waves,nlforce)
  do j=1,3
    delta=0.0_dp;delta(j)=h
    call synthetic_nonlocal_energy(r0+delta,eplus)
    call synthetic_nonlocal_energy(r0-delta,eminus)
    nlnumeric(j)=-(eplus-eminus)/(2*h)
  end do
  write(*,'(a,es24.16)')'nonlocal_force_max_error=',maxval(abs(nlforce(1,:)-nlnumeric))
  if(maxval(abs(nlforce(1,:)-nlnumeric))>2e-9_dp)error stop 'HALF: analytic non-local force regression'
  write(*,'(a,es24.16)')'ewald_ev=',e1
  write(*,'(a,es24.16)')'finite_temperature_mu_ev=',mu
  write(*,'(a)')'HALF occupation and Ewald checks passed'
contains
  subroutine synthetic_nonlocal_energy(position,energy)
    real(dp),intent(in)::position(3);real(dp),intent(out)::energy
    complex(dp)::metric(2,2)
    integer::g
    cc=(0.0_dp,0.0_dp)
    do g=1,3
      phase=dot_product(basis%qvectors(g,:),position)
      cc=cc+conjg(projector0(:,g)*cmplx(cos(phase),-sin(phase),dp))*waves(g,1)
    end do
    metric=cmplx(paw(1)%dij_atom(1,:,:)-eval(1)*paw(1)%qij,0.0_dp,dp)
    energy=0.4_dp*focc(1)*real(dot_product(cc,matmul(metric,cc)),dp)
  end subroutine synthetic_nonlocal_energy
end program
