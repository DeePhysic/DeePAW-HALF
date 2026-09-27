module half_onsite_density
  ! Occupation-dependent one-centre radial density for the CPU PAW functional.
  ! The CUDA radial path uses the same multipole and compensation conventions.
  use half_kinds,only:dp
  use half_types,only:potcar_t
  use half_constants,only:pi,felect
  use half_uspp,only:gaunt_numeric,radial_weights,bessel_root,compensation_coefficients,sph_bessel
  implicit none
  private
  public::build_onsite_radial_multipoles,atomic_onsite_occupation, &
    add_onsite_compensation,onsite_hartree_energy,hartree_multipole
contains
  subroutine atomic_onsite_occupation(p,occupation)
    type(potcar_t),intent(in)::p
    real(dp),allocatable,intent(out)::occupation(:,:)
    integer,allocatable::channel(:),angular(:),magnetic(:)
    integer::nlm,ich,l,m,a,b
    nlm=0
    do ich=1,p%channels
      nlm=nlm+2*p%lps(ich)+1
    end do
    allocate(channel(nlm),angular(nlm),magnetic(nlm),occupation(nlm,nlm))
    occupation=0.0_dp;a=0
    do ich=1,p%channels
      l=p%lps(ich)
      do m=-l,l
        a=a+1;channel(a)=ich;angular(a)=l;magnetic(a)=m
      end do
    end do
    do b=1,nlm
      do a=1,nlm
        if(angular(a)==angular(b).and.magnetic(a)==magnetic(b)) &
          occupation(a,b)=p%qato(channel(a),channel(b))
      end do
    end do
  end subroutine

  subroutine build_onsite_radial_multipoles(p,occupation,lmax,ae,ps,compensation)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),allocatable,intent(out)::ae(:,:),ps(:,:),compensation(:)
    integer,allocatable::channel(:),angular(:),magnetic(:)
    integer::nlm,nrad,naug,ich,m,a,b,l,k
    real(dp)::g,q
    if(lmax<0)error stop 'HALF: onsite radial multipoles require lmax >= 0'
    nlm=0
    do ich=1,p%channels
      nlm=nlm+2*p%lps(ich)+1
    end do
    if(size(occupation,1)/=nlm.or.size(occupation,2)/=nlm) &
      error stop 'HALF: onsite occupation dimensions do not match POTCAR projectors'
    if(maxval(abs(occupation-transpose(occupation)))>1.0e-9_dp) &
      error stop 'HALF: onsite occupation matrix must be symmetric'
    nrad=size(p%rgrid);naug=(lmax+1)*(lmax+1)
    allocate(channel(nlm),angular(nlm),magnetic(nlm),ae(nrad,naug),ps(nrad,naug),compensation(naug))
    ae=0.0_dp;ps=0.0_dp;compensation=0.0_dp
    a=0
    do ich=1,p%channels
      l=p%lps(ich)
      do m=-l,l
        a=a+1;channel(a)=ich;angular(a)=l;magnetic(a)=m
      end do
    end do
    k=0
    do l=0,lmax
      do m=-l,l
        k=k+1
        do b=1,nlm
          do a=1,nlm
            if(abs(angular(a)-angular(b))>l.or.l>angular(a)+angular(b))cycle
            if(mod(angular(a)+angular(b)+l,2)/=0)cycle
            g=occupation(a,b)*gaunt_numeric(angular(a),magnetic(a),angular(b),magnetic(b),l,m)
            if(abs(g)<1.0e-15_dp)cycle
            ae(:,k)=ae(:,k)+g*p%wae(:,channel(a))*p%wae(:,channel(b))
            ps(:,k)=ps(:,k)+g*p%wps(:,channel(a))*p%wps(:,channel(b))
            if(l+1<=size(p%qpaw_l,3))then
              q=p%qpaw_l(channel(a),channel(b),l+1)
              compensation(k)=compensation(k)+g*q
            end if
          end do
        end do
      end do
    end do
    ! POTCAR core densities are radial coefficients of the Y_00 channel.
    ae(:,1)=ae(:,1)+p%rhoae
    ps(:,1)=ps(:,1)+p%rhops
  end subroutine

  subroutine add_onsite_compensation(p,compensation,ps)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::compensation(:)
    real(dp),intent(inout)::ps(:,:)
    integer::l,m,k,ir,lmax
    real(dp)::z1,z2,c1,c2,r,shape_value
    if(any(shape(ps)/=[size(p%rgrid),size(compensation)])) &
      error stop 'HALF: onsite compensation dimensions mismatch'
    lmax=nint(sqrt(real(size(compensation),dp)))-1
    if((lmax+1)*(lmax+1)/=size(compensation))error stop 'HALF: invalid onsite multipole count'
    k=0
    do l=0,lmax
      call bessel_root(l,1,z1);call bessel_root(l,2,z2)
      call compensation_coefficients(l,p%paw_rmax,z1,z2,c1,c2)
      do m=-l,l
        k=k+1
        if(abs(compensation(k))<1.0e-18_dp)cycle
        do ir=1,size(p%rgrid)
          r=p%rgrid(ir)
          if(r>p%paw_rmax)cycle
          shape_value=c1*sph_bessel(l,z1*r/p%paw_rmax)+c2*sph_bessel(l,z2*r/p%paw_rmax)
          ps(ir,k)=ps(ir,k)+compensation(k)*shape_value*r*r
        end do
      end do
    end do
  end subroutine

  subroutine onsite_hartree_energy(p,occupation,lmax,ae_energy,ps_energy)
    type(potcar_t),intent(in)::p
    real(dp),intent(in)::occupation(:,:)
    integer,intent(in)::lmax
    real(dp),intent(out)::ae_energy,ps_energy
    real(dp),allocatable::ae(:,:),ps(:,:),compensation(:),weights(:)
    integer::k,l
    call build_onsite_radial_multipoles(p,occupation,lmax,ae,ps,compensation)
    call add_onsite_compensation(p,compensation,ps)
    allocate(weights(size(p%rgrid)))
    call radial_weights(p%rgrid,weights)
    ae_energy=0.0_dp;ps_energy=0.0_dp
    do k=1,size(compensation)
      l=int(sqrt(real(k-1,dp)))
      ae_energy=ae_energy+hartree_multipole(p%rgrid,weights,ae(:,k),l)
      ps_energy=ps_energy+hartree_multipole(p%rgrid,weights,ps(:,k),l)
    end do
  end subroutine

  real(dp) function hartree_multipole(r,w,q,l)result(energy)
    real(dp),intent(in)::r(:),w(:),q(:)
    integer,intent(in)::l
    real(dp),allocatable::inner(:),outer(:)
    real(dp)::potential
    integer::i,n
    n=size(r);allocate(inner(n),outer(n))
    inner(1)=w(1)*q(1)*r(1)**l
    do i=2,n
      inner(i)=inner(i-1)+w(i)*q(i)*r(i)**l
    end do
    outer(n)=0.0_dp
    do i=n-1,1,-1
      outer(i)=outer(i+1)+w(i+1)*q(i+1)/r(i+1)**(l+1)
    end do
    energy=0.0_dp
    do i=1,n
      potential=4.0_dp*pi/real(2*l+1,dp)* &
        (inner(i)/r(i)**(l+1)+r(i)**l*outer(i))
      energy=energy+0.5_dp*felect*w(i)*q(i)*potential
    end do
  end function
end module
