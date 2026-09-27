module half_potcar_interp
  use half_kinds,only:dp
  implicit none
  private
  public::vasp_local_spline_second,vasp_local_spline_eval,vasp_four_point_eval
contains
  subroutine vasp_local_spline_second(y,xmax,m2)
    ! Equivalent to VASP SPLCOF(...,Y1P=0): the reciprocal local
    ! pseudopotential has zero first derivative at G=0 and a natural
    ! boundary at the final table point.
    real(dp),intent(in)::y(:),xmax
    real(dp),intent(out)::m2(:)
    real(dp),allocatable::lower(:),diag(:),upper(:),rhs(:)
    real(dp)::h,w
    integer::n,i
    n=size(y)
    if(size(m2)/=n.or.n<3)error stop 'HALF: invalid local-potential spline table'
    h=xmax/real(n,dp)
    allocate(lower(n),diag(n),upper(n),rhs(n));lower=0.0_dp;upper=0.0_dp;rhs=0.0_dp
    diag(1)=2.0_dp*h;upper(1)=h;rhs(1)=6.0_dp*(y(2)-y(1))/h
    do i=2,n-1
      lower(i)=h;diag(i)=4.0_dp*h;upper(i)=h
      rhs(i)=6.0_dp*((y(i+1)-y(i))/h-(y(i)-y(i-1))/h)
    end do
    diag(n)=1.0_dp;rhs(n)=0.0_dp
    do i=2,n
      w=lower(i)/diag(i-1);diag(i)=diag(i)-w*upper(i-1);rhs(i)=rhs(i)-w*rhs(i-1)
    end do
    m2(n)=rhs(n)/diag(n)
    do i=n-1,1,-1;m2(i)=(rhs(i)-upper(i)*m2(i+1))/diag(i);end do
  end subroutine vasp_local_spline_second

  pure real(dp) function vasp_local_spline_eval(y,m2,xmax,x)result(value)
    real(dp),intent(in)::y(:),m2(:),xmax,x
    real(dp)::h,t
    integer::n,k
    n=size(y);h=xmax/real(n,dp);k=max(1,min(n-1,int(floor(x/h))+1));t=(x-real(k-1,dp)*h)/h
    value=(1.0_dp-t)*y(k)+t*y(k+1)+h*h*(((1.0_dp-t)**3-(1.0_dp-t))*m2(k)+ &
      (t**3-t)*m2(k+1))/6.0_dp
  end function vasp_local_spline_eval

  pure real(dp) function vasp_four_point_eval(y,xmax,x)result(value)
    ! VASP RHOATO/FORHAR interpolation for PSPCOR and PSPRHO.
    real(dp),intent(in)::y(:),xmax,x
    real(dp)::arg,rem,v1,v2,v3,v4,t0,t1,t2,t3
    integer::n,k
    n=size(y);arg=x*real(n,dp)/xmax+1.0_dp;k=max(2,min(n-2,int(arg)));rem=arg-real(k,dp)
    v1=y(k-1);v2=y(k);v3=y(k+1);v4=y(k+2)
    t0=v2;t1=(6.0_dp*v3-2.0_dp*v1-3.0_dp*v2-v4)/6.0_dp
    t2=(v1+v3-2.0_dp*v2)/2.0_dp;t3=(v4-v1+3.0_dp*(v2-v3))/6.0_dp
    value=t0+rem*(t1+rem*(t2+rem*t3))
  end function vasp_four_point_eval
end module half_potcar_interp
