module half_onsite_functional
  ! Incremental nonlinear PAW one-centre functional about POTCAR QATO.
  ! Linear kinetic/ionic terms remain in DION; fixed smooth-potential
  ! compensation coupling remains in QDEP.  The residual and its derivative
  ! are evaluated from the same discrete AE/PS Hartree + spherical PBE energy.
  use half_kinds,only:dp
  use half_constants,only:pi,autoa,hatoev
  use half_types,only:potcar_t,crystal_t
  use half_onsite_density,only:atomic_onsite_occupation,build_onsite_radial_multipoles, &
    add_onsite_compensation,hartree_multipole,onsite_angular_table_t,prepare_onsite_angular_table
  use half_paw_atomic_energy,only:radial_pbe_energy
  use half_xc,only:pbe_energy_density
  use half_uspp,only:radial_weights,augmentation_occupancy_t,gaunt_numeric
  use half_paw,only:paw_species_t,onsite_species_correction_t
  implicit none
  private
  public::evaluate_onsite_nonlinearity,evaluate_onsite_corrections,apply_onsite_corrections_cpu,angular_pbe_energy
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
    type(onsite_angular_table_t)::angular_table
    integer::it,iat,n
    if(size(potcars)/=size(occupancy).or.size(potcars)/=size(crystal%counts)) &
      error stop 'HALF: onsite species count mismatch'
    allocate(correction(size(potcars)))
    double_counting_delta=0.0_dp;energy_delta=0.0_dp
    do it=1,size(potcars)
      if(lmax>0)call prepare_onsite_angular_table(lmax,angular_table)
      call atomic_onsite_occupation(potcars(it),reference)
      call functional_energy(potcars(it),reference,lmax,reference_energy,angular_table)
      call functional_gradient(potcars(it),reference,lmax,reference_gradient,angular_table)
      n=size(occupancy(it)%matrix,2)
      if(size(occupancy(it)%matrix,1)/=crystal%counts(it)) &
        error stop 'HALF: onsite atom count mismatch'
      allocate(correction(it)%dij(crystal%counts(it),n,n))
      do iat=1,crystal%counts(it)
        call evaluate_onsite_nonlinearity(potcars(it),occupancy(it)%matrix(iat,:,:),lmax,dij,dc,e, &
          reference,reference_energy,reference_gradient,angular_table)
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
      reference_in,reference_energy_in,reference_gradient_in,angular_table_in)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),allocatable,intent(out)::dij_delta(:,:)
    real(dp),intent(out)::double_counting_delta,energy_delta
    real(dp),intent(in),optional::reference_in(:,:),reference_energy_in,reference_gradient_in(:,:)
    type(onsite_angular_table_t),intent(in),optional::angular_table_in
    real(dp),allocatable::reference(:,:),gradient(:,:),reference_gradient(:,:)
    type(onsite_angular_table_t)::angular_table
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
    if(present(angular_table_in))then
      angular_table=angular_table_in
    else if(lmax>0)then
      call prepare_onsite_angular_table(lmax,angular_table)
    end if
    call functional_energy(p,occupation,lmax,current_energy,angular_table)
    if(present(reference_energy_in))then
      reference_energy=reference_energy_in
    else
      call functional_energy(p,reference,lmax,reference_energy,angular_table)
    end if
    call functional_gradient(p,occupation,lmax,gradient,angular_table)
    if(present(reference_gradient_in))then
      reference_gradient=reference_gradient_in
    else
      call functional_gradient(p,reference,lmax,reference_gradient,angular_table)
    end if
    dij_delta=gradient-reference_gradient
    dij_delta=0.5_dp*(dij_delta+transpose(dij_delta))
    energy_delta=current_energy-reference_energy-sum(reference_gradient*(occupation-reference))
    double_counting_delta=energy_delta-sum(occupation*dij_delta)
  end subroutine

  subroutine functional_energy(p,occupation,lmax,energy,angular_table)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),intent(out)::energy
    type(onsite_angular_table_t),intent(in),optional::angular_table
    real(dp),allocatable::ae(:,:),ps(:,:),compensation(:)
    call build_onsite_radial_multipoles(p,occupation,lmax,ae,ps,compensation)
    call functional_energy_from_radial(p,ae,ps,compensation,energy,angular_table)
  end subroutine

  subroutine functional_energy_from_radial(p,ae,ps,compensation,energy,angular_table)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::ae(:,:),ps(:,:),compensation(:)
    real(dp),intent(out)::energy
    type(onsite_angular_table_t),intent(in),optional::angular_table
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
    if(size(compensation)>1.and.present(angular_table))then
      ae_xc=angular_pbe_energy(ae,p%rgrid,weights,angular_table)
      ps_xc=angular_pbe_energy(pseudo,p%rgrid,weights,angular_table)
    else
      ae_xc=radial_pbe_energy(ae_density,p%rgrid,weights)
      ps_xc=radial_pbe_energy(ps_density,p%rgrid,weights)
    end if
    energy=(ae_hartree-ps_hartree)+(ae_xc-ps_xc)
  end subroutine

  subroutine functional_gradient(p,occupation,lmax,gradient,angular_table)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),allocatable,intent(out)::gradient(:,:)
    type(onsite_angular_table_t),intent(in),optional::angular_table
    real(dp),allocatable::ae(:,:),ps(:,:),compensation(:),ae_shift(:,:),ps_shift(:,:),q_shift(:)
    real(dp),allocatable::plus_ae(:,:),minus_ae(:,:),plus_ps(:,:),minus_ps(:,:),plus_q(:),minus_q(:)
    real(dp),allocatable::pseudo(:,:),xc_ae_gradient(:,:),xc_ps_gradient(:,:),weights(:)
    integer,allocatable::channel(:),angular(:),magnetic(:)
    real(dp)::step,plus,minus,factor,g,q
    integer::a,b,n,k,l,m,ich
    n=size(occupation,1);allocate(gradient(n,n),channel(n),angular(n),magnetic(n));gradient=0.0_dp
    call build_onsite_radial_multipoles(p,occupation,lmax,ae,ps,compensation)
    if(lmax>0.and.present(angular_table))then
      pseudo=ps
      call add_onsite_compensation(p,compensation,pseudo)
      allocate(weights(size(p%rgrid)))
      call radial_weights(p%rgrid,weights)
      call angular_pbe_density_gradient(ae,p%rgrid,weights,angular_table,xc_ae_gradient)
      call angular_pbe_density_gradient(pseudo,p%rgrid,weights,angular_table,xc_ps_gradient)
    end if
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
        if(lmax>0.and.present(angular_table))then
          call add_onsite_compensation(p,q_shift,ps_shift)
          plus_ps=pseudo+step*ps_shift;minus_ps=pseudo-step*ps_shift
          call hartree_energy_from_radial(p%rgrid,weights,plus_ae,plus_ps,plus)
          call hartree_energy_from_radial(p%rgrid,weights,minus_ae,minus_ps,minus)
          plus=plus+2.0_dp*step*(sum(xc_ae_gradient*ae_shift)-sum(xc_ps_gradient*ps_shift))
        else
          plus_ps=ps+step*ps_shift;minus_ps=ps-step*ps_shift
          plus_q=compensation+step*q_shift;minus_q=compensation-step*q_shift
          call functional_energy_from_radial(p,plus_ae,plus_ps,plus_q,plus,angular_table)
          call functional_energy_from_radial(p,minus_ae,minus_ps,minus_q,minus,angular_table)
        end if
        factor=2.0_dp*step;if(a/=b)factor=4.0_dp*step
        gradient(a,b)=(plus-minus)/factor
        gradient(b,a)=gradient(a,b)
      end do
    end do
  end subroutine

  subroutine hartree_energy_from_radial(r,w,ae,ps,energy)
    real(dp),intent(in)::r(:),w(:),ae(:,:),ps(:,:)
    real(dp),intent(out)::energy
    integer::k,l
    energy=0.0_dp
    do k=1,size(ae,2)
      l=int(sqrt(real(k-1,dp)))
      energy=energy+hartree_multipole(r,w,ae(:,k),l)-hartree_multipole(r,w,ps(:,k),l)
    end do
  end subroutine

  subroutine angular_pbe_density_gradient(radial,r,w,table,gradient)
    real(dp),intent(in)::radial(:,:),r(:),w(:)
    type(onsite_angular_table_t),intent(in)::table
    real(dp),allocatable,intent(out)::gradient(:,:)
    real(dp),allocatable::rho(:,:),theta(:,:),phi(:,:)
    real(dp)::n_raw,n,dr,dt,df,rau,scale,sigma,hn,hs,fn,fs,radial_factor, &
      c_n,c_r,c_t,c_p,delta_r
    integer::i,j,k,lo,hi,nrad,naug
    nrad=size(r);naug=size(radial,2)
    rho=matmul(radial,table%y)
    theta=matmul(radial,table%dtheta)
    phi=matmul(radial,table%dphi)
    allocate(gradient(nrad,naug));gradient=0.0_dp
    do i=1,nrad
      lo=max(1,i-1);hi=min(nrad,i+1)
      if(i==1)hi=2
      if(i==nrad)lo=nrad-1
      rau=r(i)/autoa;delta_r=(r(hi)-r(lo))/autoa
      do j=1,table%npoints
        n_raw=rho(i,j)/(r(i)*r(i))*autoa**3
        n=max(n_raw,1.0e-30_dp)
        dr=(rho(hi,j)/r(hi)**2-rho(lo,j)/r(lo)**2)/delta_r*autoa**3
        dt=theta(i,j)/r(i)**2*autoa**3/rau
        df=phi(i,j)/r(i)**2*autoa**3/(rau*table%sin_theta(j))
        sigma=dr*dr+dt*dt+df*df
        hn=max(1.0e-12_dp,1.0e-6_dp*n)
        hs=max(1.0e-10_dp,1.0e-5_dp*sigma)
        fn=(pbe_energy_density(n+hn,sigma)-pbe_energy_density(max(1.0e-30_dp,n-hn),sigma))/ &
          (n+hn-max(1.0e-30_dp,n-hn))
        if(n_raw<=1.0e-30_dp)fn=0.0_dp
        fs=(pbe_energy_density(n,sigma+hs)-pbe_energy_density(n,max(0.0_dp,sigma-hs)))/ &
          (sigma+hs-max(0.0_dp,sigma-hs))
        scale=(w(i)/autoa)*rau*rau*hatoev*table%weight(j)
        radial_factor=autoa**3/r(i)**2
        c_n=scale*fn*radial_factor
        c_r=scale*2.0_dp*fs*dr*autoa**3/delta_r
        c_t=scale*2.0_dp*fs*dt*radial_factor/rau
        c_p=scale*2.0_dp*fs*df*radial_factor/(rau*table%sin_theta(j))
        do k=1,naug
          gradient(i,k)=gradient(i,k)+c_n*table%y(k,j)+c_t*table%dtheta(k,j)+c_p*table%dphi(k,j)
          gradient(hi,k)=gradient(hi,k)+c_r*table%y(k,j)/r(hi)**2
          gradient(lo,k)=gradient(lo,k)-c_r*table%y(k,j)/r(lo)**2
        end do
      end do
    end do
  end subroutine

  real(dp) function angular_pbe_energy(radial,r,w,table)result(energy)
    real(dp),intent(in)::radial(:,:),r(:),w(:)
    type(onsite_angular_table_t),intent(in)::table
    real(dp),allocatable::rho(:,:),theta(:,:),phi(:,:)
    real(dp)::na,dr,dt,df,rau,scale,sigma
    integer::i,j,lo,hi,nrad,naug
    nrad=size(r);naug=size(radial,2)
    if(size(table%y,1)/=naug)error stop 'HALF: angular XC multipole count mismatch'
    rho=matmul(radial,table%y)
    theta=matmul(radial,table%dtheta)
    phi=matmul(radial,table%dphi)
    energy=0.0_dp
    do i=1,nrad
      lo=max(1,i-1);hi=min(nrad,i+1)
      if(i==1)hi=2
      if(i==nrad)lo=nrad-1
      rau=r(i)/autoa;scale=(w(i)/autoa)*rau*rau*hatoev
      do j=1,table%npoints
        na=max(rho(i,j)/(r(i)*r(i))*autoa**3,1.0e-30_dp)
        dr=((rho(hi,j)/(r(hi)*r(hi)))-(rho(lo,j)/(r(lo)*r(lo))))/ &
          ((r(hi)-r(lo))/autoa)*autoa**3
        dt=theta(i,j)/(r(i)*r(i))*autoa**3/rau
        df=phi(i,j)/(r(i)*r(i))*autoa**3/rau
        sigma=dr*dr+dt*dt+(df/table%sin_theta(j))**2
        energy=energy+scale*table%weight(j)*pbe_energy_density(na,sigma)
      end do
    end do
  end function
end module
