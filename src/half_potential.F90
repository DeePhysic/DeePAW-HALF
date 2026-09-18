module half_potential
  use half_kinds,only:dp,i32,i64
  use half_constants,only:pi,felect,edeps
  use half_types,only:crystal_t,charge_grid_t,potcar_t
  use half_fft,only:fft3_forward,fft3_backward
  use half_xc,only:perdew_zunger_xc,pbe_xc
  implicit none
  private
  public::build_veff_lda,build_veff_pbe
contains
  subroutine build_veff_lda(charge,potcars,crystal,veff,e_hartree,e_xc)
    type(charge_grid_t),intent(in)::charge
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    real(dp),allocatable,intent(out)::veff(:)
    real(dp),intent(out)::e_hartree,e_xc
    call build_veff(charge,potcars,crystal,.false.,veff,e_hartree,e_xc)
  end subroutine
  subroutine build_veff_pbe(charge,potcars,crystal,veff,e_hartree,e_xc)
    type(charge_grid_t),intent(in)::charge
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    real(dp),allocatable,intent(out)::veff(:)
    real(dp),intent(out)::e_hartree,e_xc
    call build_veff(charge,potcars,crystal,.true.,veff,e_hartree,e_xc)
  end subroutine
  subroutine build_veff(charge,potcars,crystal,use_pbe,veff,e_hartree,e_xc)
    type(charge_grid_t),intent(in)::charge
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    logical,intent(in)::use_pbe
    real(dp),allocatable,intent(out)::veff(:)
    real(dp),intent(out)::e_hartree,e_xc
    complex(dp),allocatable,target::rho_in(:),rho_g(:),local_g(:),local_c(:),core_g(:),core_c(:)
    real(dp),allocatable::local_r(:),core_r(:),vxc(:),m2local(:,:),m2core(:,:)
    real(dp)::q(3),g2,gn,radial,core_radial,phase,rpos(3)
    integer::i1,i2,i3,n1,n2,n3,idx,it,iat,ion0,n,nt
    n=size(charge%values); nt=size(potcars)
    allocate(rho_in(n),rho_g(n),local_g(n),local_c(n),core_g(n),core_c(n),local_r(n),core_r(n),vxc(n),veff(n))
    allocate(m2local(1000,nt),m2core(1000,nt)); m2core=0
    do it=1,nt
      call spline_second(potcars(it)%psp_local,potcars(it)%psp_gmax,m2local(:,it))
      if(potcars(it)%has_core)call spline_second(potcars(it)%pspcor,potcars(it)%psp_gmax,m2core(:,it))
    end do
    rho_in=cmplx(charge%values,0.0_dp,dp); call fft3_forward(charge%shape,rho_in,rho_g); rho_g=rho_g/real(n,dp)
    local_g=(0.0_dp,0.0_dp); core_g=(0.0_dp,0.0_dp); e_hartree=0
    idx=0
    do i1=0,charge%shape(1)-1
      n1=fft_integer(i1,charge%shape(1))
      do i2=0,charge%shape(2)-1
        n2=fft_integer(i2,charge%shape(2))
        do i3=0,charge%shape(3)-1
          n3=fft_integer(i3,charge%shape(3)); idx=idx+1
          q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal); g2=dot_product(q,q); gn=sqrt(g2)
          if(g2>1.0e-12_dp)then
            local_g(idx)=edeps*rho_g(idx)/(g2*crystal%volume)
            e_hartree=e_hartree+0.5_dp*real(conjg(rho_g(idx))*local_g(idx),dp)
          end if
          ion0=0
          do it=1,nt
            radial=spline_eval(potcars(it)%psp_local,m2local(:,it),potcars(it)%psp_gmax,gn)
            if(g2>1.0e-12_dp)radial=radial-4*pi*potcars(it)%zval*felect/g2
            if(g2<=1.0e-12_dp.or.gn>potcars(it)%psp_gmax-potcars(it)%psp_gmax/1000)radial=0
            core_radial=0
            if(potcars(it)%has_core.and.gn<=potcars(it)%psp_gmax-3*potcars(it)%psp_gmax/1000) &
              core_radial=spline_eval(potcars(it)%pspcor,m2core(:,it),potcars(it)%psp_gmax,gn)
            do iat=1,crystal%counts(it)
              rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice); phase=dot_product(q,rpos)
              local_g(idx)=local_g(idx)+radial*cmplx(cos(phase),-sin(phase),dp)/crystal%volume
              core_g(idx)=core_g(idx)+core_radial*cmplx(cos(phase),-sin(phase),dp)
            end do
            ion0=ion0+crystal%counts(it)
          end do
        end do
      end do
    end do
    local_c=local_g; call fft3_backward(charge%shape,local_c,rho_in); local_r=real(rho_in,dp)
    core_c=core_g; call fft3_backward(charge%shape,core_c,rho_in); core_r=real(rho_in,dp)
    if(use_pbe)then
      call pbe_xc(charge%values+core_r,charge%shape,crystal,vxc,e_xc)
    else
      call perdew_zunger_xc(charge%values+core_r,crystal%volume,vxc,e_xc)
    end if
    veff=local_r+vxc
  end subroutine

  pure integer function fft_integer(index0,n)
    integer,intent(in)::index0,n
    if(index0<=(n-1)/2)then; fft_integer=index0; else; fft_integer=index0-n; end if
  end function
  subroutine spline_second(y,xmax,m2)
    real(dp),intent(in)::y(:),xmax
    real(dp),intent(out)::m2(:)
    real(dp)::h,w
    real(dp),allocatable::diag(:),rhs(:),upper(:)
    integer::n,i
    n=size(y); h=xmax/real(n,dp); allocate(diag(n),rhs(n),upper(n)); diag=0; rhs=0; upper=0
    diag(2)=1; rhs(2)=((y(3)-y(2))/h-(y(2)-y(1))/h)/h
    do i=3,n-2
      diag(i)=4*h; upper(i)=h; rhs(i)=6*((y(i+1)-y(i))/h-(y(i)-y(i-1))/h)
    end do
    diag(n-1)=1; rhs(n-1)=((y(n)-y(n-1))/h-(y(n-1)-y(n-2))/h)/h
    do i=3,n-1
      w=merge(h,0.0_dp,i<=n-2)/diag(i-1); diag(i)=diag(i)-w*upper(i-1); rhs(i)=rhs(i)-w*rhs(i-1)
    end do
    m2(n-1)=rhs(n-1)/diag(n-1)
    do i=n-2,2,-1; m2(i)=(rhs(i)-upper(i)*m2(i+1))/diag(i); end do
    m2(1)=2*m2(2)-m2(3); m2(n)=2*m2(n-1)-m2(n-2)
  end subroutine
  pure function spline_eval(y,m2,xmax,x)result(value)
    real(dp),intent(in)::y(:),m2(:),xmax,x
    real(dp)::value,h,t
    integer::n,k
    n=size(y); h=xmax/real(n,dp); k=max(1,min(n-1,int(floor(x/h))+1)); t=(x-real(k-1,dp)*h)/h
    value=(1-t)*y(k)+t*y(k+1)+h*h*(((1-t)**3-(1-t))*m2(k)+(t**3-t)*m2(k+1))/6
  end function
end module half_potential
