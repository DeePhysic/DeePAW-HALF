module half_paw_atomic_energy
  use half_kinds,only:dp
  use half_constants,only:pi,autoa,hatoev,felect
  use half_types,only:potcar_t
  use half_xc,only:pbe_energy_density
  implicit none
  private
  public::compute_paw_atomic_double_counting
contains
  subroutine compute_paw_atomic_double_counting(potcars,counts,use_pbe,total,ae_total,ps_total)
    type(potcar_t),intent(in)::potcars(:)
    integer,intent(in)::counts(:)
    logical,intent(in)::use_pbe
    real(dp),intent(out)::total
    real(dp),intent(out),optional::ae_total,ps_total
    real(dp)::ae,ps,ae_one,ps_one
    integer::it
    if(size(counts)/=size(potcars))error stop 'HALF: PAW atomic-energy species dimensions differ'
    ae=0.0_dp;ps=0.0_dp
    do it=1,size(potcars)
      call species_double_counting(potcars(it),use_pbe,ae_one,ps_one)
      ae=ae+real(counts(it),dp)*ae_one;ps=ps+real(counts(it),dp)*ps_one
    end do
    total=ae+ps
    if(present(ae_total))ae_total=ae
    if(present(ps_total))ps_total=ps
  end subroutine

  subroutine species_double_counting(p,use_pbe,double_ae,double_ps)
    type(potcar_t),intent(in)::p
    logical,intent(in)::use_pbe
    real(dp),intent(out)::double_ae,double_ps
    real(dp),allocatable::rho_ae(:),rho_ps(:),val_ae(:),val_ps(:),weights(:),shape0(:)
    real(dp)::y00,aug0,eh_ae,eh_ps,xc_ae,xc_ps,vxc_ae,vxc_ps,core_xc,qae,qps
    integer::a,b,l,n
    if(.not.use_pbe)error stop 'HALF: automatic PAW atomic double counting currently requires PBE'
    n=size(p%rgrid);allocate(rho_ae(n),rho_ps(n),val_ae(n),val_ps(n),weights(n),shape0(n))
    rho_ae=0.0_dp;rho_ps=0.0_dp;y00=1.0_dp/sqrt(4.0_dp*pi);aug0=0.0_dp
    do b=1,p%channels;do a=1,p%channels
      if(p%lps(a)/=p%lps(b))cycle;l=p%lps(a)
      rho_ae=rho_ae+y00*real(2*l+1,dp)*p%qato(a,b)*p%wae(:,a)*p%wae(:,b)
      rho_ps=rho_ps+y00*real(2*l+1,dp)*p%qato(a,b)*p%wps(:,a)*p%wps(:,b)
    end do;end do
    call radial_weights(p%rgrid,weights)
    aug0=0.0_dp
    do b=1,p%channels;do a=1,p%channels
      if(p%lps(a)/=p%lps(b))cycle;l=p%lps(a)
      aug0=aug0+y00*real(2*l+1,dp)*p%qato(a,b)* &
        sum(weights*(p%wae(:,a)*p%wae(:,b)-p%wps(:,a)*p%wps(:,b)))
    end do;end do
    call compensation_shape0(p%rgrid,p%paw_rmax,shape0)
    rho_ps=rho_ps+aug0*shape0*p%rgrid*p%rgrid
    call radial_physical_density(rho_ae,p%rgrid,val_ae)
    call radial_physical_density(rho_ps,p%rgrid,val_ps)
    qae=sum(weights*4.0_dp*pi*p%rgrid*p%rgrid*val_ae)
    qps=sum(weights*4.0_dp*pi*p%rgrid*p%rgrid*val_ps)
    eh_ae=spherical_hartree(val_ae,p%rgrid,weights)
    eh_ps=spherical_hartree(val_ps,p%rgrid,weights)
    call xc_double_counting(val_ae,p%rhoae,p%rgrid,weights,xc_ae,vxc_ae,core_xc)
    double_ae=-eh_ae+xc_ae-vxc_ae-(core_xc+p%dexccore)
    call xc_double_counting(val_ps,p%rhops,p%rgrid,weights,xc_ps,vxc_ps,core_xc)
    double_ps=eh_ps-xc_ps+vxc_ps
    write(0,'(A,A,5ES18.8)')'HALF PAW atomic species ',trim(p%element),qae,qps,aug0,double_ae,double_ps
  end subroutine

  subroutine radial_physical_density(radial,r,density)
    real(dp),intent(in)::radial(:),r(:)
    real(dp),intent(out)::density(:)
    integer::i
    do i=1,size(r)
      density(i)=radial(i)/(sqrt(4.0_dp*pi)*max(r(i)*r(i),1.0e-30_dp))
    end do
  end subroutine

  real(dp) function spherical_hartree(density,r,w)result(energy)
    real(dp),intent(in)::density(:),r(:),w(:)
    real(dp),allocatable::shell(:)
    integer::i,j
    allocate(shell(size(r)));shell=4.0_dp*pi*r*r*density
    energy=0.0_dp
    do j=1,size(r);do i=1,size(r)
      energy=energy+0.5_dp*felect*w(i)*w(j)*shell(i)*shell(j)/max(r(i),r(j))
    end do;end do
  end function

  subroutine xc_double_counting(valence,core_radial,r,w,exc,dvxc,core_exc)
    real(dp),intent(in)::valence(:),core_radial(:),r(:),w(:)
    real(dp),intent(out)::exc,dvxc,core_exc
    real(dp),allocatable::core(:)
    real(dp)::t,eplus,eminus
    allocate(core(size(r)));call radial_physical_density(core_radial,r,core)
    t=1.0e-5_dp
    exc=radial_pbe_energy(core+valence,r,w)
    eplus=radial_pbe_energy(core+(1.0_dp+t)*valence,r,w)
    eminus=radial_pbe_energy(core+(1.0_dp-t)*valence,r,w)
    dvxc=(eplus-eminus)/(2.0_dp*t)
    core_exc=radial_pbe_energy(core,r,w)
  end subroutine

  real(dp) function radial_pbe_energy(density,r,w)result(energy)
    real(dp),intent(in)::density(:),r(:),w(:)
    real(dp),allocatable::nau(:),gradient(:)
    real(dp)::sigma,rau
    integer::i,n
    n=size(r);allocate(nau(n),gradient(n));nau=max(density*autoa**3,1.0e-30_dp)
    gradient(1)=(nau(2)-nau(1))/((r(2)-r(1))/autoa)
    do i=2,n-1
      gradient(i)=(nau(i+1)-nau(i-1))/((r(i+1)-r(i-1))/autoa)
    end do
    gradient(n)=(nau(n)-nau(n-1))/((r(n)-r(n-1))/autoa)
    energy=0.0_dp
    do i=1,n
      sigma=gradient(i)*gradient(i);rau=r(i)/autoa
      energy=energy+(w(i)/autoa)*4.0_dp*pi*rau*rau*pbe_energy_density(nau(i),sigma)*hatoev
    end do
  end function

  subroutine compensation_shape0(r,rc,shape)
    real(dp),intent(in)::r(:),rc
    real(dp),intent(out)::shape(:)
    real(dp)::z1,z2,c1,c2
    integer::i
    call bessel_root(0,1,z1);call bessel_root(0,2,z2)
    call compensation_coefficients(rc,z1,z2,c1,c2)
    do i=1,size(r)
      if(r(i)<=rc)then
        shape(i)=c1*sph_bessel0(z1*r(i)/rc)+c2*sph_bessel0(z2*r(i)/rc)
      else
        shape(i)=0.0_dp
      end if
    end do
  end subroutine

  subroutine compensation_coefficients(rc,z1,z2,c1,c2)
    real(dp),intent(in)::rc,z1,z2
    real(dp),intent(out)::c1,c2
    real(dp)::a11,a12,a21,a22,det
    integer,parameter::n=4000
    real(dp)::h,x
    integer::i
    a11=z1/rc*sph_bessel0_deriv(z1);a12=z2/rc*sph_bessel0_deriv(z2)
    h=rc/real(n,dp);a21=0.0_dp;a22=0.0_dp
    do i=0,n
      x=real(i,dp)*h
      a21=a21+merge(1.0_dp,merge(2.0_dp,4.0_dp,mod(i,2)==0),i==0.or.i==n)*sph_bessel0(z1*x/rc)*x*x
      a22=a22+merge(1.0_dp,merge(2.0_dp,4.0_dp,mod(i,2)==0),i==0.or.i==n)*sph_bessel0(z2*x/rc)*x*x
    end do
    a21=a21*h/3.0_dp;a22=a22*h/3.0_dp;det=a11*a22-a12*a21;c1=-a12/det;c2=a11/det
  end subroutine

  pure real(dp) function sph_bessel0(x)result(v)
    real(dp),intent(in)::x
    if(abs(x)<1.0e-8_dp)then;v=1.0_dp-x*x/6.0_dp;else;v=sin(x)/x;end if
  end function
  pure real(dp) function sph_bessel0_deriv(x)result(v)
    real(dp),intent(in)::x
    v=(x*cos(x)-sin(x))/(x*x)
  end function
  subroutine bessel_root(l,which,root)
    integer,intent(in)::l,which
    real(dp),intent(out)::root
    if(l/=0)error stop 'HALF: only the spherical compensation root is requested'
    root=real(which,dp)*pi
  end subroutine

  subroutine radial_weights(r,w)
    real(dp),intent(in)::r(:);real(dp),intent(out)::w(:)
    integer::i,n,last;real(dp)::h0,h1,s,alpha,beta,eta
    n=size(r);w=0.0_dp;last=merge(n,n-1,mod(n,2)==1)
    do i=1,last-2,2
      h0=r(i+1)-r(i);h1=r(i+2)-r(i+1);s=(h0+h1)/6.0_dp
      w(i)=w(i)+s*(2.0_dp-h1/h0);w(i+1)=w(i+1)+s*(h0+h1)**2/(h0*h1);w(i+2)=w(i+2)+s*(2.0_dp-h0/h1)
    end do
    if(mod(n,2)==0)then
      h0=r(n-1)-r(n-2);h1=r(n)-r(n-1);alpha=(2*h1*h1+3*h0*h1)/(6*(h0+h1))
      beta=(h1*h1+3*h0*h1)/(6*h0);eta=h1**3/(6*h0*(h0+h1))
      w(n)=w(n)+alpha;w(n-1)=w(n-1)+beta;w(n-2)=w(n-2)-eta
    end if
  end subroutine
end module half_paw_atomic_energy
