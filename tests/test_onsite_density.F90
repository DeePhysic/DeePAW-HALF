program test_onsite_density
  use half_kinds,only:dp
  use half_constants,only:pi
  use half_types,only:potcar_t
  use half_onsite_density,only:build_onsite_radial_multipoles,atomic_onsite_occupation,onsite_hartree_energy
  use half_onsite_functional,only:evaluate_onsite_nonlinearity
  implicit none
  type(potcar_t)::p,p_projector
  real(dp)::occupation(1,1),y00
  real(dp)::p_occupation(3,3),g20,eh_ae,eh_ps
  real(dp),allocatable::ae(:,:),ps(:,:),compensation(:),atomic_occupation(:,:),dij_delta(:,:)
  real(dp)::double_counting_delta,energy_delta,e_plus,e_minus,d_fd
  p%channels=1
  allocate(p%lps(1),p%rgrid(3),p%wae(3,1),p%wps(3,1),p%rhoae(3),p%rhops(3),p%qpaw_l(1,1,3))
  p%lps=0;p%rgrid=[0.1_dp,0.2_dp,0.3_dp];p%paw_rmax=0.3_dp
  p%wae(:,1)=[1.0_dp,2.0_dp,3.0_dp]
  p%wps(:,1)=0.5_dp*p%wae(:,1)
  p%rhoae=0.1_dp;p%rhops=0.05_dp
  p%qpaw_l=0.0_dp;p%qpaw_l(1,1,1)=2.0_dp
  allocate(p%qato(1,1));p%qato=1.0_dp
  call atomic_onsite_occupation(p,atomic_occupation)
  if(size(atomic_occupation,1)/=1.or.abs(atomic_occupation(1,1)-1.0_dp)>1.0e-12_dp) &
    error stop 'atomic s occupation mismatch'
  occupation=1.0_dp;y00=1.0_dp/sqrt(4.0_dp*pi)
  call build_onsite_radial_multipoles(p,occupation,2,ae,ps,compensation)
  if(maxval(abs(ae(:,1)-(p%rhoae+y00*p%wae(:,1)**2)))>1.0e-12_dp)error stop 'AE monopole mismatch'
  if(maxval(abs(ps(:,1)-(p%rhops+y00*p%wps(:,1)**2)))>1.0e-12_dp)error stop 'PS monopole mismatch'
  if(abs(compensation(1)-2.0_dp*y00)>1.0e-12_dp)error stop 'compensation monopole mismatch'
  if(maxval(abs(ae(:,2:)))>1.0e-12_dp.or.maxval(abs(ps(:,2:)))>1.0e-12_dp) &
    error stop 'spurious nonspherical density from s projector'
  call onsite_hartree_energy(p,occupation,0,eh_ae,eh_ps)
  if(eh_ae<=0.0_dp.or.eh_ps<=0.0_dp)error stop 'onsite Hartree energies must be positive'
  call evaluate_onsite_nonlinearity(p,occupation,0,dij_delta,double_counting_delta,energy_delta)
  if(maxval(abs(dij_delta))>1.0e-8_dp.or.abs(double_counting_delta)>1.0e-8_dp) &
    error stop 'onsite nonlinear correction is not zero at atomic reference'
  occupation=1.1_dp
  call evaluate_onsite_nonlinearity(p,occupation,0,dij_delta,double_counting_delta,energy_delta)
  if(abs(dij_delta(1,1))<1.0e-9_dp.or.abs(energy_delta)<1.0e-9_dp) &
    error stop 'onsite nonlinear correction did not respond to occupation'
  d_fd=dij_delta(1,1)
  occupation=1.1001_dp
  call evaluate_onsite_nonlinearity(p,occupation,0,dij_delta,double_counting_delta,e_plus)
  occupation=1.0999_dp
  call evaluate_onsite_nonlinearity(p,occupation,0,dij_delta,double_counting_delta,e_minus)
  if(abs(d_fd-(e_plus-e_minus)/0.0002_dp)>1.0e-4_dp) &
    error stop 'onsite Hamiltonian is not the derivative of onsite energy'
  p_projector%channels=1
  allocate(p_projector%lps(1),p_projector%rgrid(3),p_projector%wae(3,1), &
    p_projector%wps(3,1),p_projector%rhoae(3),p_projector%rhops(3),p_projector%qpaw_l(1,1,3))
  p_projector%lps=1;p_projector%rgrid=p%rgrid;p_projector%paw_rmax=p%paw_rmax
  p_projector%wae=1.0_dp;p_projector%wps=0.0_dp
  p_projector%rhoae=0.0_dp;p_projector%rhops=0.0_dp
  p_projector%qpaw_l=0.0_dp;p_projector%qpaw_l(1,1,3)=3.0_dp
  allocate(p_projector%qato(1,1));p_projector%qato=0.25_dp
  call atomic_onsite_occupation(p_projector,atomic_occupation)
  if(size(atomic_occupation,1)/=3)error stop 'atomic p occupation dimension mismatch'
  if(abs(atomic_occupation(1,1)-0.25_dp)>1.0e-12_dp.or. &
    abs(atomic_occupation(2,2)-0.25_dp)>1.0e-12_dp.or. &
    abs(atomic_occupation(3,3)-0.25_dp)>1.0e-12_dp)error stop 'atomic p occupation mismatch'
  p_occupation=0.0_dp;p_occupation(2,2)=1.0_dp
  call build_onsite_radial_multipoles(p_projector,p_occupation,0,ae,ps,compensation)
  if(size(ae,2)/=1)error stop 'onsite lmax=0 did not truncate the angular channels'
  call build_onsite_radial_multipoles(p_projector,p_occupation,2,ae,ps,compensation)
  g20=sqrt(5.0_dp)/(5.0_dp*sqrt(pi))
  if(maxval(abs(ae(:,7)-g20))>1.0e-12_dp)error stop 'p-projector quadrupole mismatch'
  if(abs(compensation(7)-3.0_dp*g20)>1.0e-12_dp)error stop 'quadrupole compensation mismatch'
  if(maxval(abs(ae(:,2:4)))>1.0e-12_dp)error stop 'spurious odd multipole from p projector'
  call onsite_hartree_energy(p_projector,p_occupation,2,eh_ae,eh_ps)
  if(eh_ae<=0.0_dp.or.eh_ps<0.0_dp)error stop 'p-projector Hartree energies are invalid'
  print '(A)','onsite radial multipoles: PASS'
end program
