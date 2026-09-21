module half_paw
  use half_kinds, only: dp, i32
  use half_constants, only: pi
  use half_types, only: crystal_t, plane_wave_basis_t, potcar_t
  implicit none
  private
  public :: paw_species_t, build_paw_operators
  type :: paw_species_t
    integer(i32) :: natoms=0, nlm=0
    complex(dp), allocatable :: projectors(:,:,:)
    real(dp), allocatable :: dij(:,:), qij(:,:)
    ! DION is species-wide, but the MIMIC_US correction DeltaD depends on
    ! Veff around each ion.  Keep the assembled D matrix per atom.
    real(dp), allocatable :: dij_atom(:,:,:)
    real(dp), allocatable :: ddij_atom(:,:,:,:)
  end type
contains
  subroutine build_paw_operators(potcars,crystal,basis,paw)
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    type(plane_wave_basis_t),intent(in)::basis
    type(paw_species_t),allocatable,intent(out)::paw(:)
    integer::it,ich,m,ilm,jlm,jch,mj,iat,ig,ion0,nlm,l
    real(dp)::q(3),qn,theta,phi,rval,ylm,phase,rpos(3),pref
    real(dp),allocatable::spline_m2(:,:)
    complex(dp)::ilfac
    allocate(paw(size(potcars))); ion0=0; pref=1.0_dp/sqrt(crystal%volume)
    do it=1,size(potcars)
      nlm=0
      do ich=1,potcars(it)%channels; nlm=nlm+2*potcars(it)%lps(ich)+1; end do
      paw(it)%natoms=crystal%counts(it); paw(it)%nlm=nlm
      allocate(paw(it)%projectors(crystal%counts(it),nlm,basis%npw),paw(it)%dij(nlm,nlm),paw(it)%qij(nlm,nlm), &
        paw(it)%dij_atom(crystal%counts(it),nlm,nlm),paw(it)%ddij_atom(crystal%counts(it),nlm,nlm,3))
      allocate(spline_m2(size(potcars(it)%pspnl,1),potcars(it)%channels))
      do ich=1,potcars(it)%channels
        call spline_second(potcars(it)%pspnl(:,ich),potcars(it)%pspnl_gmax,spline_m2(:,ich))
      end do
      paw(it)%projectors=(0.0_dp,0.0_dp); paw(it)%dij=0; paw(it)%qij=0;paw(it)%ddij_atom=0; ilm=0
      do ich=1,potcars(it)%channels
        l=potcars(it)%lps(ich)
        do m=-l,l
          ilm=ilm+1; jlm=0
          do jch=1,potcars(it)%channels
            do mj=-potcars(it)%lps(jch),potcars(it)%lps(jch)
              jlm=jlm+1
              if(l==potcars(it)%lps(jch).and.m==mj)then
                paw(it)%dij(ilm,jlm)=potcars(it)%dion(ich,jch)
                paw(it)%qij(ilm,jlm)=potcars(it)%qpaw(ich,jch)
              end if
            end do
          end do
          ilfac=(-cmplx(0.0_dp,1.0_dp,dp))**l
          do iat=1,crystal%counts(it)
            rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice)
            do ig=1,basis%npw
              q=basis%qvectors(ig,:); qn=sqrt(dot_product(q,q))
              if(qn>1.0e-10_dp)then
                theta=acos(max(-1.0_dp,min(1.0_dp,q(3)/qn))); phi=atan2(q(2),q(1))
              else
                theta=0; phi=0
              end if
              rval=spline_eval(potcars(it)%pspnl(:,ich),spline_m2(:,ich),potcars(it)%pspnl_gmax,qn)
              if(qn>0.99_dp*potcars(it)%pspnl_gmax)rval=0
              ylm=real_ylm(l,m,theta,phi); phase=dot_product(q,rpos)
              paw(it)%projectors(iat,ilm,ig)=pref*ilfac*rval*ylm*cmplx(cos(phase),-sin(phase),dp)
            end do
          end do
        end do
      end do
      do iat=1,paw(it)%natoms
        paw(it)%dij_atom(iat,:,:)=paw(it)%dij
      end do
      deallocate(spline_m2)
      ion0=ion0+crystal%counts(it)
    end do
  end subroutine

  subroutine spline_second(y,xmax,m2)
    real(dp),intent(in)::y(:),xmax
    real(dp),intent(out)::m2(:)
    real(dp)::h,w
    real(dp),allocatable::diag(:),rhs(:),upper(:)
    integer::n,i
    n=size(y); h=xmax/real(n,dp); allocate(diag(n),rhs(n),upper(n)); diag=0; rhs=0; upper=0
    ! Uniform-grid not-a-knot conditions eliminate the endpoint second
    ! derivatives, leaving a tridiagonal system for M(2:n-1).
    diag(2)=1; rhs(2)=((y(3)-y(2))/h-(y(2)-y(1))/h)/h
    do i=3,n-2
      diag(i)=4*h; upper(i)=h; rhs(i)=6*((y(i+1)-y(i))/h-(y(i)-y(i-1))/h)
    end do
    diag(n-1)=1; rhs(n-1)=((y(n)-y(n-1))/h-(y(n-1)-y(n-2))/h)/h
    do i=3,n-1
      w=merge(h,0.0_dp,i<=n-2)/diag(i-1)
      diag(i)=diag(i)-w*upper(i-1); rhs(i)=rhs(i)-w*rhs(i-1)
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

  pure function real_ylm(l,m,theta,phi)result(y)
    integer,intent(in)::l,m
    real(dp),intent(in)::theta,phi
    real(dp)::y,p,norm,x
    integer::am
    am=abs(m); x=cos(theta); p=assoc_legendre(l,am,x)
    norm=sqrt(real(2*l+1,dp)/(4*pi)*factorial(l-am)/factorial(l+am))
    if(m>0)then
      y=sqrt(2.0_dp)*(-1.0_dp)**m*norm*p*cos(real(am,dp)*phi)
    else if(m<0)then
      y=sqrt(2.0_dp)*(-1.0_dp)**m*norm*p*sin(real(am,dp)*phi)
    else
      y=norm*p
    end if
  end function
  pure function assoc_legendre(l,m,x)result(p)
    integer,intent(in)::l,m
    real(dp),intent(in)::x
    real(dp)::p,pmm,pmmp1,pll,fact
    integer::i,ll
    pmm=1
    if(m>0)then
      fact=1; do i=1,m; pmm=-pmm*fact*sqrt(max(0.0_dp,(1-x)*(1+x))); fact=fact+2; end do
    end if
    if(l==m)then; p=pmm; return; end if
    pmmp1=x*real(2*m+1,dp)*pmm; if(l==m+1)then; p=pmmp1; return; end if
    do ll=m+2,l; pll=(real(2*ll-1,dp)*x*pmmp1-real(ll+m-1,dp)*pmm)/real(ll-m,dp); pmm=pmmp1; pmmp1=pll; end do
    p=pmmp1
  end function
  pure function factorial(n)result(v)
    integer,intent(in)::n
    real(dp)::v
    integer::i
    v=1; do i=2,n; v=v*i; end do
  end function
end module half_paw
