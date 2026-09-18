module half_xc
  use half_kinds,only:dp
  use half_constants,only:pi,autoa,hatoev
  use half_types,only:crystal_t
  use half_fft,only:fft3_forward,fft3_backward
  implicit none
  private
  public::perdew_zunger_xc,pbe_xc
contains
  subroutine perdew_zunger_xc(rho,omega,vxc,exc)
    real(dp),intent(in)::rho(:),omega
    real(dp),intent(out)::vxc(:),exc
    real(dp)::density,rs,ex,vx,ec,vc,sr,denom,lnrs,rhoau,dvau
    integer::i
    if(size(vxc)/=size(rho))error stop 'HALF: XC size mismatch'
    vxc=0; exc=0; dvau=(omega/real(size(rho),dp))/(autoa**3)
    do i=1,size(rho)
      density=rho(i)/omega; rhoau=density*autoa**3
      if(rhoau<=1.0e-30_dp)cycle
      rs=(3.0_dp/(4*pi*rhoau))**(1.0_dp/3.0_dp)
      ex=-0.4582_dp/rs; vx=-(4.0_dp/3.0_dp)*0.4582_dp/rs
      if(rs>=1)then
        sr=sqrt(rs); denom=1+1.0529_dp*sr+0.3334_dp*rs
        ec=-0.1423_dp/denom
        vc=ec*(1+(7.0_dp/6.0_dp)*1.0529_dp*sr+(4.0_dp/3.0_dp)*0.3334_dp*rs)/denom
      else
        lnrs=log(rs); ec=0.0311_dp*lnrs-0.048_dp+0.002_dp*rs*lnrs-0.0116_dp*rs
        vc=0.0311_dp*lnrs+(-0.048_dp-0.0311_dp/3)+(2.0_dp/3)*0.002_dp*rs*lnrs+(-0.0232_dp-0.002_dp)*rs/3
      end if
      vxc(i)=(vx+vc)*hatoev; exc=exc+(ex+ec)*rhoau*dvau*hatoev
    end do
  end subroutine

  subroutine pbe_xc(rho,shape,crystal,vxc,exc)
    real(dp),intent(in)::rho(:)
    integer,intent(in)::shape(3)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(out)::vxc(:),exc
    complex(dp),allocatable,target::work(:),freq(:),deriv(:)
    real(dp),allocatable::density(:),sigma(:),grad(:,:),fs(:),fn(:),flux(:),edens(:)
    real(dp)::q(3),rhoau,dn,np,nm,ds,sp,sm,dvau
    integer::n,alpha,i,i1,i2,i3,n1,n2,n3,idx
    n=size(rho); allocate(work(n),freq(n),deriv(n),density(n),sigma(n),grad(3,n),fs(n),fn(n),flux(n),edens(n))
    density=max(rho/crystal%volume*autoa**3,1.0e-20_dp); sigma=0
    work=cmplx(rho/crystal%volume,0.0_dp,dp); call fft3_forward(shape,work,freq)
    do alpha=1,3
      idx=0
      do i1=0,shape(1)-1; n1=fft_integer(i1,shape(1))
        do i2=0,shape(2)-1; n2=fft_integer(i2,shape(2))
          do i3=0,shape(3)-1; n3=fft_integer(i3,shape(3)); idx=idx+1
            q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal)
            deriv(idx)=cmplx(0.0_dp,q(alpha),dp)*freq(idx)
          end do
        end do
      end do
      call fft3_backward(shape,deriv,work); grad(alpha,:)=real(work,dp)/real(n,dp)*autoa**4
      sigma=sigma+grad(alpha,:)**2
    end do
    do i=1,n
      if(rho(i)/crystal%volume*autoa**3<=1.0e-20_dp)then
        density(i)=1.0e-20_dp; sigma(i)=0
      end if
      dn=max(abs(density(i))*2.0e-5_dp,1.0e-12_dp); nm=max(density(i)-dn,1.0e-30_dp); np=density(i)+dn
      fn(i)=(pbe_energy_density(np,sigma(i))-pbe_energy_density(nm,sigma(i)))/(np-nm)
      ds=max(abs(sigma(i))*2.0e-5_dp,1.0e-14_dp); sm=max(sigma(i)-ds,0.0_dp); sp=sigma(i)+ds
      fs(i)=(pbe_energy_density(density(i),sp)-pbe_energy_density(density(i),sm))/(sp-sm)
      edens(i)=pbe_energy_density(density(i),sigma(i))
      if(rho(i)/crystal%volume*autoa**3<=1.0e-20_dp)fs(i)=0
    end do
    vxc=fn
    do alpha=1,3
      flux=2*fs*grad(alpha,:); work=cmplx(flux,0.0_dp,dp); call fft3_forward(shape,work,freq)
      idx=0
      do i1=0,shape(1)-1; n1=fft_integer(i1,shape(1))
        do i2=0,shape(2)-1; n2=fft_integer(i2,shape(2))
          do i3=0,shape(3)-1; n3=fft_integer(i3,shape(3)); idx=idx+1
            q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal)
            deriv(idx)=cmplx(0.0_dp,q(alpha),dp)*freq(idx)
          end do
        end do
      end do
      call fft3_backward(shape,deriv,work); vxc=vxc-autoa*real(work,dp)/real(n,dp)
    end do
    do i=1,n
      if(rho(i)/crystal%volume*autoa**3<=1.0e-20_dp)vxc(i)=0
    end do
    vxc=vxc*hatoev; dvau=(crystal%volume/real(n,dp))/autoa**3
    exc=0
    do i=1,n
      if(rho(i)/crystal%volume*autoa**3>1.0e-20_dp)exc=exc+edens(i)*dvau*hatoev
    end do
  end subroutine

  pure function pbe_energy_density(density,sigma)result(f)
    real(dp),intent(in)::density,sigma
    real(dp)::f,n,s,rs,kf,exlda,s2,fx,ec,de,ks2,t2,aa,at2,ratio,hh
    real(dp),parameter::mu=0.2195149727645171_dp,kappa=0.804_dp,beta=0.06672455060314922_dp
    real(dp),parameter::gamma=(1.0_dp-log(2.0_dp))/(pi*pi)
    n=max(density,1.0e-30_dp); s=max(sigma,0.0_dp); rs=(3/(4*pi*n))**(1.0_dp/3.0_dp)
    kf=(3*pi*pi*n)**(1.0_dp/3.0_dp); exlda=-3*kf/(4*pi); s2=s/(4*kf*kf*n*n)
    fx=1+kappa-kappa/(1+mu*s2/kappa); call pw92(rs,ec,de); ks2=4*kf/pi; t2=s/(4*ks2*n*n)
    aa=beta/gamma/(exp(-ec/gamma)-1); at2=aa*t2
    if(t2==0)then; ratio=0
    else if(at2<1.0e100_dp)then; ratio=t2*(1+at2)/(1+at2+at2*at2)
    else; ratio=(1/aa)*(1+1/at2)/(1+1/at2+(1/at2)**2); end if
    hh=gamma*log(1+beta/gamma*ratio); f=n*(exlda*fx+ec+hh)
  end function
  pure subroutine pw92(rs,ec,dec)
    real(dp),intent(in)::rs
    real(dp),intent(out)::ec,dec
    real(dp)::sr,q,dq,la
    real(dp),parameter::a=.031091_dp,a1=.21370_dp,b1=7.5957_dp,b2=3.5876_dp,b3=1.6382_dp,b4=.49294_dp
    sr=sqrt(rs); q=2*a*(b1*sr+b2*rs+b3*rs*sr+b4*rs*rs)
    dq=a*(b1/sr+2*b2+3*b3*sr+4*b4*rs); la=log(1+1/q)
    ec=-2*a*(1+a1*rs)*la; dec=-2*a*a1*la+2*a*(1+a1*rs)*dq/(q*q+q)
  end subroutine
  pure integer function fft_integer(index0,n)
    integer,intent(in)::index0,n
    if(index0<=(n-1)/2)then; fft_integer=index0; else; fft_integer=index0-n; end if
  end function
end module half_xc
