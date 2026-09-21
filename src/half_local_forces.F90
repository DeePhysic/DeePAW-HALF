module half_local_forces
  use half_kinds,only:dp
  use half_constants,only:pi,felect
  use half_types,only:crystal_t,charge_grid_t,potcar_t
  use half_fft,only:fft3_forward
  implicit none
  private
  public::local_ionic_forces,nlcc_forces
contains
  subroutine local_ionic_forces(charge,potcars,crystal,forces,energy)
    ! VASP FORLOC in a full complex FFT representation.  The input density
    ! is transformed once; differentiation acts only on exp[-i G.R_I].
    type(charge_grid_t),intent(in)::charge
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(out)::forces(:,:)
    real(dp),intent(out),optional::energy
    complex(dp),allocatable,target::rho_r(:),rho_g(:)
    real(dp),allocatable::charge_c(:),m2(:)
    real(dp)::q(3),g2,gn,radial,phase,rpos(3),term,total_energy
    integer::i1,i2,i3,n1,n2,n3,idx,it,iat,ion0,n
    if(size(forces,1)/=crystal%nions.or.size(forces,2)/=3) &
      error stop 'HALF: local force dimensions do not match the crystal'
    n=size(charge%values);allocate(charge_c(n),rho_r(n),rho_g(n));call reorder(charge%values,charge%shape,charge_c)
    rho_r=cmplx(charge_c,0.0_dp,dp);call fft3_forward(charge%shape,rho_r,rho_g);rho_g=rho_g/real(n,dp)
    forces=0.0_dp;total_energy=0.0_dp
    do it=1,size(potcars)
      allocate(m2(size(potcars(it)%psp_local)));call spline_second(potcars(it)%psp_local,potcars(it)%psp_gmax,m2)
      ion0=sum(crystal%counts(:it-1))
      do iat=1,crystal%counts(it)
        rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice);idx=0
        do i1=0,charge%shape(1)-1
          n1=fft_integer(i1,charge%shape(1))
          do i2=0,charge%shape(2)-1
            n2=fft_integer(i2,charge%shape(2))
            do i3=0,charge%shape(3)-1
              n3=fft_integer(i3,charge%shape(3));idx=idx+1
              q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal);g2=dot_product(q,q)
              if(g2<=1.0e-12_dp)cycle
              gn=sqrt(g2);if(gn>potcars(it)%psp_gmax-potcars(it)%psp_gmax/real(size(m2),dp))cycle
              radial=spline_eval(potcars(it)%psp_local,m2,potcars(it)%psp_gmax,gn)- &
                4*pi*potcars(it)%zval*felect/g2
              phase=dot_product(q,rpos)
              term=radial/crystal%volume*aimag(conjg(rho_g(idx))*cmplx(cos(phase),-sin(phase),dp))
              forces(ion0+iat,:)=forces(ion0+iat,:)-q*term
              total_energy=total_energy+radial/crystal%volume*real(conjg(rho_g(idx))* &
                cmplx(cos(phase),-sin(phase),dp),dp)
            end do
          end do
        end do
      end do
      deallocate(m2)
    end do
    if(present(energy))energy=total_energy
  end subroutine local_ionic_forces

  subroutine nlcc_forces(vxc,potcars,crystal,shape,forces)
    ! VASP FORCOR/FORHAR(LPAR=.TRUE.) in a full complex FFT layout:
    ! F_I = - integral v_xc(r) d n_core,I(r-R_I)/dR_I dr.
    ! CHGCAR/core arrays store Omega*n, so the integral is a grid average.
    real(dp),intent(in)::vxc(:)
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    integer,intent(in)::shape(3)
    real(dp),intent(out)::forces(:,:)
    complex(dp),allocatable,target::work(:),vxc_g(:)
    real(dp),allocatable::m2(:)
    real(dp)::q(3),g2,gn,radial,phase,rpos(3),term
    integer::n,idx,i1,i2,i3,n1,n2,n3,it,iat,ion0
    n=size(vxc)
    if(n/=product(shape).or.size(forces,1)/=crystal%nions.or.size(forces,2)/=3) &
      error stop 'HALF: NLCC force dimensions do not match the FFT grid/crystal'
    allocate(work(n),vxc_g(n));work=cmplx(vxc,0.0_dp,dp);call fft3_forward(shape,work,vxc_g)
    vxc_g=vxc_g/real(n,dp);forces=0.0_dp
    do it=1,size(potcars)
      ion0=sum(crystal%counts(:it-1));if(.not.potcars(it)%has_core)cycle
      allocate(m2(size(potcars(it)%pspcor)));call spline_second(potcars(it)%pspcor,potcars(it)%psp_gmax,m2)
      do iat=1,crystal%counts(it)
        rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice);idx=0
        do i1=0,shape(1)-1
          n1=fft_integer(i1,shape(1))
          do i2=0,shape(2)-1
            n2=fft_integer(i2,shape(2))
            do i3=0,shape(3)-1
              n3=fft_integer(i3,shape(3));idx=idx+1
              q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal);g2=dot_product(q,q)
              if(g2<=1.0e-12_dp)cycle
              gn=sqrt(g2);if(gn>potcars(it)%psp_gmax-3*potcars(it)%psp_gmax/real(size(m2),dp))cycle
              radial=spline_eval(potcars(it)%pspcor,m2,potcars(it)%psp_gmax,gn)
              phase=dot_product(q,rpos)
              term=radial*aimag(conjg(vxc_g(idx))*cmplx(cos(phase),-sin(phase),dp))
              forces(ion0+iat,:)=forces(ion0+iat,:)-q*term
            end do
          end do
        end do
      end do
      deallocate(m2)
    end do
  end subroutine nlcc_forces

  subroutine reorder(input,shape,output)
    real(dp),intent(in)::input(:);integer,intent(in)::shape(3);real(dp),intent(out)::output(:)
    integer::i1,i2,i3,source,destination
    do i1=0,shape(1)-1;do i2=0,shape(2)-1;do i3=0,shape(3)-1
      source=i1+shape(1)*(i2+shape(2)*i3)+1
      destination=i1*shape(2)*shape(3)+i2*shape(3)+i3+1;output(destination)=input(source)
    end do;end do;end do
  end subroutine reorder
  pure integer function fft_integer(index0,n)result(value)
    integer,intent(in)::index0,n
    if(index0<=(n-1)/2)then;value=index0;else;value=index0-n;end if
  end function fft_integer
  subroutine spline_second(y,xmax,m2)
    real(dp),intent(in)::y(:),xmax;real(dp),intent(out)::m2(:)
    real(dp)::h,w;real(dp),allocatable::diag(:),rhs(:),upper(:);integer::n,i
    n=size(y);h=xmax/real(n,dp);allocate(diag(n),rhs(n),upper(n));diag=0;rhs=0;upper=0
    diag(2)=1;rhs(2)=((y(3)-y(2))/h-(y(2)-y(1))/h)/h
    do i=3,n-2;diag(i)=4*h;upper(i)=h;rhs(i)=6*((y(i+1)-y(i))/h-(y(i)-y(i-1))/h);end do
    diag(n-1)=1;rhs(n-1)=((y(n)-y(n-1))/h-(y(n-1)-y(n-2))/h)/h
    do i=3,n-1;w=merge(h,0.0_dp,i<=n-2)/diag(i-1);diag(i)=diag(i)-w*upper(i-1);rhs(i)=rhs(i)-w*rhs(i-1);end do
    m2(n-1)=rhs(n-1)/diag(n-1);do i=n-2,2,-1;m2(i)=(rhs(i)-upper(i)*m2(i+1))/diag(i);end do
    m2(1)=2*m2(2)-m2(3);m2(n)=2*m2(n-1)-m2(n-2)
  end subroutine spline_second
  pure real(dp) function spline_eval(y,m2,xmax,x)result(value)
    real(dp),intent(in)::y(:),m2(:),xmax,x;real(dp)::h,t;integer::n,k
    n=size(y);h=xmax/real(n,dp);k=max(1,min(n-1,int(floor(x/h))+1));t=(x-real(k-1,dp)*h)/h
    value=(1-t)*y(k)+t*y(k+1)+h*h*(((1-t)**3-(1-t))*m2(k)+(t**3-t)*m2(k+1))/6
  end function spline_eval
end module half_local_forces
