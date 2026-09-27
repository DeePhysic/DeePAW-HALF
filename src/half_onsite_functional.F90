module half_onsite_functional
  ! Incremental nonlinear PAW one-centre functional about POTCAR QATO.
  ! Linear kinetic/ionic terms remain in DION; fixed smooth-potential
  ! compensation coupling remains in QDEP.  The residual and its derivative
  ! are evaluated from the same discrete AE/PS Hartree + spherical PBE energy.
  use half_kinds,only:dp
  use half_constants,only:pi
  use half_types,only:potcar_t,crystal_t
  use half_onsite_density,only:atomic_onsite_occupation,build_onsite_radial_multipoles, &
    add_onsite_compensation,hartree_multipole
  use half_paw_atomic_energy,only:radial_pbe_energy
  use half_uspp,only:radial_weights,augmentation_occupancy_t,gaunt_numeric
  use half_paw,only:paw_species_t,onsite_species_correction_t
  implicit none
  private
  public::evaluate_onsite_nonlinearity,evaluate_onsite_corrections,apply_onsite_corrections_cpu
contains
  subroutine evaluate_onsite_corrections(potcars,crystal,occupancy,lmax,correction,double_counting_delta,energy_delta)
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    type(augmentation_occupancy_t),intent(in)::occupancy(:)
    integer,intent(in)::lmax
    type(onsite_species_correction_t),allocatable,intent(out)::correction(:)
    real(dp),intent(out)::double_counting_delta,energy_delta
    real(dp),allocatable::dij(:,:),reference(:,:),reference_gradient(:,:)
    real(dp)::dc,e,reference_energy
    integer::it,iat,n
    if(size(potcars)/=size(occupancy).or.size(potcars)/=size(crystal%counts)) &
      error stop 'HALF: onsite species count mismatch'
    allocate(correction(size(potcars)))
    double_counting_delta=0.0_dp;energy_delta=0.0_dp
    do it=1,size(potcars)
      call atomic_onsite_occupation(potcars(it),reference)
      call functional_energy(potcars(it),reference,lmax,reference_energy)
      call functional_gradient(potcars(it),reference,lmax,reference_gradient)
      n=size(occupancy(it)%matrix,2)
      if(size(occupancy(it)%matrix,1)/=crystal%counts(it)) &
        error stop 'HALF: onsite atom count mismatch'
      allocate(correction(it)%dij(crystal%counts(it),n,n))
      do iat=1,crystal%counts(it)
        call evaluate_onsite_nonlinearity(potcars(it),occupancy(it)%matrix(iat,:,:),lmax,dij,dc,e, &
          reference,reference_energy,reference_gradient)
        correction(it)%dij(iat,:,:)=dij
        double_counting_delta=double_counting_delta+dc
        energy_delta=energy_delta+e
      end do
    end do
  end subroutine

  subroutine apply_onsite_corrections_cpu(paw,correction)
    type(paw_species_t),intent(inout)::paw(:)
    type(onsite_species_correction_t),intent(in)::correction(:)
    integer::it
    if(size(paw)/=size(correction))error stop 'HALF: onsite correction species count mismatch'
    do it=1,size(paw)
      if(any(shape(paw(it)%dij_atom)/=shape(correction(it)%dij))) &
        error stop 'HALF: onsite correction Dij dimensions mismatch'
      paw(it)%dij_atom=paw(it)%dij_atom+correction(it)%dij
    end do
  end subroutine

  subroutine evaluate_onsite_nonlinearity(p,occupation,lmax,dij_delta,double_counting_delta,energy_delta, &
      reference_in,reference_energy_in,reference_gradient_in)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),allocatable,intent(out)::dij_delta(:,:)
    real(dp),intent(out)::double_counting_delta,energy_delta
    real(dp),intent(in),optional::reference_in(:,:),reference_energy_in,reference_gradient_in(:,:)
    real(dp),allocatable::reference(:,:),gradient(:,:),reference_gradient(:,:)
    real(dp)::current_energy,reference_energy
    if(lmax<0)error stop 'HALF: nonlinear onsite functional requires lmax >= 0'
    if(present(reference_in))then
      if(.not.present(reference_energy_in).or..not.present(reference_gradient_in)) &
        error stop 'HALF: incomplete cached onsite reference'
      reference=reference_in
    else
      call atomic_onsite_occupation(p,reference)
    end if
    if(any(shape(reference)/=shape(occupation)))error stop 'HALF: onsite functional dimension mismatch'
    call functional_energy(p,occupation,lmax,current_energy)
    if(present(reference_energy_in))then
      reference_energy=reference_energy_in
    else
      call functional_energy(p,reference,lmax,reference_energy)
    end if
    call functional_gradient(p,occupation,lmax,gradient)
    if(present(reference_gradient_in))then
      reference_gradient=reference_gradient_in
    else
      call functional_gradient(p,reference,lmax,reference_gradient)
    end if
    dij_delta=gradient-reference_gradient
    dij_delta=0.5_dp*(dij_delta+transpose(dij_delta))
    energy_delta=current_energy-reference_energy-sum(reference_gradient*(occupation-reference))
    double_counting_delta=energy_delta-sum(occupation*dij_delta)
  end subroutine

  subroutine functional_energy(p,occupation,lmax,energy)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),intent(out)::energy
    real(dp),allocatable::ae(:,:),ps(:,:),compensation(:)
    call build_onsite_radial_multipoles(p,occupation,lmax,ae,ps,compensation)
    call functional_energy_from_radial(p,ae,ps,compensation,energy)
  end subroutine

  subroutine functional_energy_from_radial(p,ae,ps,compensation,energy)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::ae(:,:),ps(:,:),compensation(:)
    real(dp),intent(out)::energy
    real(dp),allocatable::pseudo(:,:),weights(:),ae_density(:),ps_density(:)
    real(dp)::ae_hartree,ps_hartree,ae_xc,ps_xc
    integer::i,n,k,l
    pseudo=ps
    call add_onsite_compensation(p,compensation,pseudo)
    n=size(p%rgrid);allocate(weights(n),ae_density(n),ps_density(n))
    call radial_weights(p%rgrid,weights)
    ae_hartree=0.0_dp;ps_hartree=0.0_dp
    do k=1,size(compensation)
      l=int(sqrt(real(k-1,dp)))
      ae_hartree=ae_hartree+hartree_multipole(p%rgrid,weights,ae(:,k),l)
      ps_hartree=ps_hartree+hartree_multipole(p%rgrid,weights,pseudo(:,k),l)
    end do
    do i=1,n
      ae_density(i)=ae(i,1)/(sqrt(4.0_dp*pi)*p%rgrid(i)**2)
      ps_density(i)=pseudo(i,1)/(sqrt(4.0_dp*pi)*p%rgrid(i)**2)
    end do
    ae_xc=radial_pbe_energy(ae_density,p%rgrid,weights)
    ps_xc=radial_pbe_energy(ps_density,p%rgrid,weights)
    energy=(ae_hartree-ps_hartree)+(ae_xc-ps_xc)
  end subroutine

  subroutine functional_gradient(p,occupation,lmax,gradient)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),allocatable,intent(out)::gradient(:,:)
    real(dp),allocatable::ae(:,:),ps(:,:),compensation(:),ae_shift(:,:),ps_shift(:,:),q_shift(:)
    real(dp),allocatable::plus_ae(:,:),minus_ae(:,:),plus_ps(:,:),minus_ps(:,:),plus_q(:),minus_q(:)
    integer,allocatable::channel(:),angular(:),magnetic(:)
    real(dp)::step,plus,minus,factor,g,q
    integer::a,b,n,k,l,m,ich
    n=size(occupation,1);allocate(gradient(n,n),channel(n),angular(n),magnetic(n));gradient=0.0_dp
    call build_onsite_radial_multipoles(p,occupation,lmax,ae,ps,compensation)
    allocate(ae_shift(size(ae,1),size(ae,2)),ps_shift(size(ps,1),size(ps,2)),q_shift(size(compensation)))
    allocate(plus_ae(size(ae,1),size(ae,2)),minus_ae(size(ae,1),size(ae,2)), &
      plus_ps(size(ps,1),size(ps,2)),minus_ps(size(ps,1),size(ps,2)), &
      plus_q(size(compensation)),minus_q(size(compensation)))
    k=0
    do ich=1,p%channels
      l=p%lps(ich)
      do m=-l,l
        k=k+1;channel(k)=ich;angular(k)=l;magnetic(k)=m
      end do
    end do
    do b=1,n
      do a=1,b
        ! Large-core species can have a radial functional in the keV range;
        ! use a central step large enough to keep cancellation below the
        ! micro-eV-scale onsite fixed-point tolerance.
        step=1.0e-3_dp*max(1.0_dp,abs(occupation(a,b)))
        factor=1.0_dp;if(a/=b)factor=2.0_dp
        ae_shift=0.0_dp;ps_shift=0.0_dp;q_shift=0.0_dp;k=0
        do l=0,lmax
          do m=-l,l
            k=k+1
            if(abs(angular(a)-angular(b))>l.or.l>angular(a)+angular(b))cycle
            if(mod(angular(a)+angular(b)+l,2)/=0)cycle
            g=factor*gaunt_numeric(angular(a),magnetic(a),angular(b),magnetic(b),l,m)
            if(abs(g)<1.0e-15_dp)cycle
            ae_shift(:,k)=g*p%wae(:,channel(a))*p%wae(:,channel(b))
            ps_shift(:,k)=g*p%wps(:,channel(a))*p%wps(:,channel(b))
            if(l+1<=size(p%qpaw_l,3))then
              q=p%qpaw_l(channel(a),channel(b),l+1)
              q_shift(k)=g*q
            end if
          end do
        end do
        plus_ae=ae+step*ae_shift;minus_ae=ae-step*ae_shift
        plus_ps=ps+step*ps_shift;minus_ps=ps-step*ps_shift
        plus_q=compensation+step*q_shift;minus_q=compensation-step*q_shift
        call functional_energy_from_radial(p,plus_ae,plus_ps,plus_q,plus)
        call functional_energy_from_radial(p,minus_ae,minus_ps,minus_q,minus)
        factor=2.0_dp*step;if(a/=b)factor=4.0_dp*step
        gradient(a,b)=(plus-minus)/factor
        gradient(b,a)=gradient(a,b)
      end do
    end do
  end subroutine
end module
