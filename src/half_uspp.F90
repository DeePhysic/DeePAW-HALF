module half_uspp
  use half_kinds,only:dp,i32
  use half_constants,only:pi
  use half_types,only:potcar_t,crystal_t
  use half_paw,only:paw_species_t
  use half_math,only:inverse3
  use half_fft,only:fft3_forward,fft3_backward
  implicit none
  private
  public::build_uspp_dij_cpu,augmentation_occupancy_t,initialize_augmentation_occupancy, &
    accumulate_augmentation_occupancy,add_augmentation_density
  type::augmentation_occupancy_t
    real(dp),allocatable::matrix(:,:,:)
  end type
contains
  subroutine initialize_augmentation_occupancy(paw,occupancy)
    type(paw_species_t),intent(in)::paw(:)
    type(augmentation_occupancy_t),allocatable,intent(out)::occupancy(:)
    integer::it
    allocate(occupancy(size(paw)))
    do it=1,size(paw)
      allocate(occupancy(it)%matrix(paw(it)%natoms,paw(it)%nlm,paw(it)%nlm));occupancy(it)%matrix=0.0_dp
    end do
  end subroutine

  subroutine accumulate_augmentation_occupancy(paw,eigenvectors,occupations,kweight,occupancy)
    type(paw_species_t),intent(in)::paw(:)
    complex(dp),intent(in)::eigenvectors(:,:)
    real(dp),intent(in)::occupations(:),kweight
    type(augmentation_occupancy_t),intent(inout)::occupancy(:)
    complex(dp),allocatable::c(:)
    integer::it,iat,ib,a,b,ig,nlm
    do it=1,size(paw)
      nlm=paw(it)%nlm;allocate(c(nlm))
      do iat=1,paw(it)%natoms;do ib=1,size(occupations)
        if(occupations(ib)==0.0_dp)cycle;c=(0.0_dp,0.0_dp)
        do a=1,nlm;do ig=1,size(eigenvectors,1)
          c(a)=c(a)+conjg(paw(it)%projectors(iat,a,ig))*eigenvectors(ig,ib)
        end do;end do
        do b=1,nlm;do a=1,nlm
          occupancy(it)%matrix(iat,a,b)=occupancy(it)%matrix(iat,a,b)+ &
            kweight*occupations(ib)*real(conjg(c(a))*c(b),dp)
        end do;end do
      end do;end do
      deallocate(c)
    end do
  end subroutine

  subroutine add_augmentation_density(potcars,crystal,shape,occupancy,density_f,integrated_charge)
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    integer(i32),intent(in)::shape(3)
    type(augmentation_occupancy_t),intent(in)::occupancy(:)
    real(dp),intent(inout)::density_f(:)
    real(dp),intent(out),optional::integrated_charge
    integer::it,iat,nlm,lmax,naug,nrad,pairs,ich,l,m,a,b,k,ir,i1,i2,i3,index,ion0
    integer,allocatable::chan(:),lv(:),mv(:),laug(:),maug(:)
    real(dp),allocatable::rw(:),multipole(:,:),coeff(:),root1(:),root2(:),c1(:),c2(:)
    real(dp)::moment,frac(3),df(3),cart(3),r,theta,phi,g,value,total_before,total_after
    if(size(density_f)/=product(shape).or.size(occupancy)/=size(potcars)) &
      error stop 'HALF: augmentation-density dimensions are inconsistent'
    total_before=sum(density_f)/real(size(density_f),dp);ion0=0
    do it=1,size(potcars)
      nlm=0;lmax=0
      do ich=1,potcars(it)%channels;nlm=nlm+2*potcars(it)%lps(ich)+1;lmax=max(lmax,2*potcars(it)%lps(ich));end do
      lmax=active_qdep_lmax(lmax)
      naug=(lmax+1)*(lmax+1);nrad=size(potcars(it)%rgrid);pairs=nlm*nlm
      allocate(chan(nlm),lv(nlm),mv(nlm),laug(naug),maug(naug),rw(nrad),multipole(pairs,naug), &
        coeff(naug),root1(0:lmax),root2(0:lmax),c1(0:lmax),c2(0:lmax))
      a=0;do ich=1,potcars(it)%channels;do m=-potcars(it)%lps(ich),potcars(it)%lps(ich)
        a=a+1;chan(a)=ich;lv(a)=potcars(it)%lps(ich);mv(a)=m
      end do;end do
      call radial_weights(potcars(it)%rgrid,rw);multipole=0.0_dp;k=0
      do l=0,lmax
        call bessel_root(l,1,root1(l));call bessel_root(l,2,root2(l))
        call compensation_coefficients(l,potcars(it)%paw_rmax,root1(l),root2(l),c1(l),c2(l))
        do m=-l,l;k=k+1;laug(k)=l;maug(k)=m;end do
      end do
      do b=1,nlm;do a=1,nlm;k=0
        do l=0,lmax;do m=-l,l;k=k+1
          if(abs(lv(a)-lv(b))<=l.and.l<=lv(a)+lv(b).and.mod(lv(a)+lv(b)+l,2)==0)then
            moment=potcars(it)%qpaw_l(chan(a),chan(b),l+1)
            multipole((b-1)*nlm+a,k)=gaunt_numeric(lv(a),mv(a),lv(b),mv(b),l,m)*moment
          end if
        end do;end do
      end do;end do
      do iat=1,crystal%counts(it)
        coeff=0.0_dp
        do k=1,naug;do b=1,nlm;do a=1,nlm
          coeff(k)=coeff(k)+occupancy(it)%matrix(iat,a,b)*multipole((b-1)*nlm+a,k)
        end do;end do;end do
        do i3=0,shape(3)-1;do i2=0,shape(2)-1;do i1=0,shape(1)-1
          frac=[real(i1,dp)/shape(1),real(i2,dp)/shape(2),real(i3,dp)/shape(3)]
          df=frac-crystal%positions(ion0+iat,:);df=df-anint(df);cart=matmul(df,crystal%lattice);r=sqrt(sum(cart*cart))
          if(r<=potcars(it)%paw_rmax)then
            if(r>1.0e-14_dp)then;theta=acos(max(-1.0_dp,min(1.0_dp,cart(3)/r)));phi=atan2(cart(2),cart(1))
            else;theta=0.0_dp;phi=0.0_dp;end if
            value=0.0_dp
            do k=1,naug;l=laug(k);g=c1(l)*sph_bessel(l,root1(l)*r/potcars(it)%paw_rmax)+ &
                c2(l)*sph_bessel(l,root2(l)*r/potcars(it)%paw_rmax)
              value=value+coeff(k)*g*ylm(l,maug(k),cos(theta),phi)
            end do
            index=i1+shape(1)*(i2+shape(2)*i3)+1;density_f(index)=density_f(index)+crystal%volume*value
          end if
        end do;end do;end do
      end do
      ion0=ion0+crystal%counts(it);deallocate(chan,lv,mv,laug,maug,rw,multipole,coeff,root1,root2,c1,c2)
    end do
    total_after=sum(density_f)/real(size(density_f),dp)
    if(present(integrated_charge))integrated_charge=total_after-total_before
  end subroutine add_augmentation_density
  subroutine build_uspp_dij_cpu(veff,shape,potcars,crystal,paw)
    real(dp),intent(in)::veff(:)
    integer(i32),intent(in)::shape(3)
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    type(paw_species_t),intent(inout)::paw(:)
    complex(dp),allocatable,target::work(:),freq(:)
    real(dp),allocatable::coeff(:)
    integer::n,i,t,i1,i2,i3
    real(dp)::b1,b2,b3
    n=size(veff)
#ifdef HALF_QDEP_DIRECT_GRID
    do i=1,size(potcars)
      call build_species_direct(veff,shape,potcars(i),crystal,i,paw(i))
    end do
    return
#else
    allocate(work(n),freq(n),coeff(n));work=cmplx(veff,0.0_dp,dp)
    call fft3_forward(shape,work,freq)
    do i=1,n
      t=i-1;i3=mod(t,shape(3));t=t/shape(3);i2=mod(t,shape(2));i1=t/shape(2)
      b1=(2.0_dp+cos(2.0_dp*pi*real(i1,dp)/real(shape(1),dp)))/3.0_dp
      b2=(2.0_dp+cos(2.0_dp*pi*real(i2,dp)/real(shape(2),dp)))/3.0_dp
      b3=(2.0_dp+cos(2.0_dp*pi*real(i3,dp)/real(shape(3),dp)))/3.0_dp
      freq(i)=freq(i)/(b1*b2*b3)
    end do
    call fft3_backward(shape,freq,work);coeff=real(work,dp)/real(n,dp)
    do i=1,size(potcars)
      call build_species(coeff,shape,potcars(i),crystal,i,paw(i))
    end do
#endif
  end subroutine

  subroutine build_species_direct(veff,shape,potcar,crystal,itype,paw)
    ! VASP SETDIJ construction for PAW/MIMIC_US: put the normalized
    ! compensation functions on the actual FFT points around each ion and
    ! contract them directly with Veff.  This avoids the additional
    ! spline-to-radial-shell and angular-quadrature approximations.
    real(dp),intent(in)::veff(:)
    integer(i32),intent(in)::shape(3)
    type(potcar_t),intent(in)::potcar
    type(crystal_t),intent(in)::crystal
    integer,intent(in)::itype
    type(paw_species_t),intent(inout)::paw
    integer::nlm,lmax,naug,nrad,natoms,pairs,ich,l,m,a,b,k,iat,ion0
    integer::i1,i2,i3,index,alpha
    integer,allocatable::chan(:),lv(:),mv(:),laug(:),maug(:)
    real(dp),allocatable::rw(:),multipole(:,:),root1(:),root2(:),c1(:),c2(:),radial_v(:),radial_grad(:,:)
    real(dp)::moment,frac(3),df(3),cart(3),cart_minus(3),cart_plus(3),kernel,scale,h
    nlm=0;lmax=0
    do ich=1,potcar%channels
      nlm=nlm+2*potcar%lps(ich)+1;lmax=max(lmax,2*potcar%lps(ich))
    end do
    lmax=active_qdep_lmax(lmax)
    if(nlm/=paw%nlm)error stop 'HALF: direct-grid QDEP channel mismatch'
    naug=(lmax+1)*(lmax+1);nrad=size(potcar%rgrid);natoms=crystal%counts(itype);pairs=nlm*nlm
    allocate(chan(nlm),lv(nlm),mv(nlm),laug(naug),maug(naug),rw(nrad),multipole(pairs,naug), &
      root1(0:lmax),root2(0:lmax),c1(0:lmax),c2(0:lmax),radial_v(naug),radial_grad(naug,3))
    a=0
    do ich=1,potcar%channels;do m=-potcar%lps(ich),potcar%lps(ich)
      a=a+1;chan(a)=ich;lv(a)=potcar%lps(ich);mv(a)=m
    end do;end do
    call radial_weights(potcar%rgrid,rw);k=0
    do l=0,lmax
      call bessel_root(l,1,root1(l));call bessel_root(l,2,root2(l))
      call compensation_coefficients(l,potcar%paw_rmax,root1(l),root2(l),c1(l),c2(l))
      do m=-l,l;k=k+1;laug(k)=l;maug(k)=m;end do
    end do
    multipole=0.0_dp
    do b=1,nlm;do a=1,nlm;k=0
      do l=0,lmax;do m=-l,l;k=k+1
        if(abs(lv(a)-lv(b))<=l.and.l<=lv(a)+lv(b).and.mod(lv(a)+lv(b)+l,2)==0)then
          moment=potcar%qpaw_l(chan(a),chan(b),l+1)
          multipole((b-1)*nlm+a,k)=gaunt_numeric(lv(a),mv(a),lv(b),mv(b),l,m)*moment
        end if
      end do;end do
    end do;end do
    ion0=sum(crystal%counts(:itype-1));scale=crystal%volume/real(product(shape),dp)
    ! This differentiates the compact QDEP kernel, not a displaced atom or
    ! a repeated eigensolve.  The small Cartesian stencil is independent of
    ! the ionic finite-difference validation path.
    h=1.0e-5_dp
    do iat=1,natoms
      radial_v=0.0_dp;radial_grad=0.0_dp
      do i1=0,shape(1)-1;do i2=0,shape(2)-1;do i3=0,shape(3)-1
        index=i1*shape(2)*shape(3)+i2*shape(3)+i3+1
        frac=[real(i1,dp)/real(shape(1),dp),real(i2,dp)/real(shape(2),dp),real(i3,dp)/real(shape(3),dp)]
        df=frac-crystal%positions(ion0+iat,:);df=df-anint(df);cart=matmul(df,crystal%lattice)
        if(sum(cart*cart)<=potcar%paw_rmax*potcar%paw_rmax)then
          do k=1,naug
            kernel=qdep_kernel(laug(k),maug(k),cart,potcar%paw_rmax,root1(laug(k)),root2(laug(k)), &
              c1(laug(k)),c2(laug(k)))
            radial_v(k)=radial_v(k)+scale*veff(index)*kernel
            do alpha=1,3
              cart_minus=cart;cart_plus=cart;cart_minus(alpha)=cart_minus(alpha)-h;cart_plus(alpha)=cart_plus(alpha)+h
              radial_grad(k,alpha)=radial_grad(k,alpha)+scale*veff(index)* &
                (qdep_kernel(laug(k),maug(k),cart_minus,potcar%paw_rmax,root1(laug(k)),root2(laug(k)), &
                   c1(laug(k)),c2(laug(k)))- &
                 qdep_kernel(laug(k),maug(k),cart_plus,potcar%paw_rmax,root1(laug(k)),root2(laug(k)), &
                   c1(laug(k)),c2(laug(k))))/(2.0_dp*h)
            end do
          end do
        end if
      end do;end do;end do
      paw%dij_atom(iat,:,:)=paw%dij
      do b=1,nlm;do a=1,nlm
        paw%dij_atom(iat,a,b)=paw%dij_atom(iat,a,b)+sum(multipole((b-1)*nlm+a,:)*radial_v)
        do alpha=1,3
          paw%ddij_atom(iat,a,b,alpha)=sum(multipole((b-1)*nlm+a,:)*radial_grad(:,alpha))
        end do
      end do;end do
      paw%dij_atom(iat,:,:)=0.5_dp*(paw%dij_atom(iat,:,:)+transpose(paw%dij_atom(iat,:,:)))
      do alpha=1,3
        paw%ddij_atom(iat,:,:,alpha)=0.5_dp*(paw%ddij_atom(iat,:,:,alpha)+transpose(paw%ddij_atom(iat,:,:,alpha)))
      end do
    end do
  end subroutine build_species_direct

  real(dp) function qdep_kernel(l,m,cart,rc,z1,z2,c1,c2)result(value)
    integer,intent(in)::l,m
    real(dp),intent(in)::cart(3),rc,z1,z2,c1,c2
    real(dp)::r,x,phi,g
    r=sqrt(sum(cart*cart))
    if(r>rc)then;value=0.0_dp;return;end if
    g=c1*sph_bessel(l,z1*r/rc)+c2*sph_bessel(l,z2*r/rc)
    if(r<1.0e-14_dp)then
      if(l==0.and.m==0)then;value=g/sqrt(4.0_dp*pi);else;value=0.0_dp;end if
    else
      x=max(-1.0_dp,min(1.0_dp,cart(3)/r));phi=atan2(cart(2),cart(1));value=g*ylm(l,m,x,phi)
    end if
  end function qdep_kernel

  subroutine build_species(coeff,shape,potcar,crystal,itype,paw)
    real(dp),intent(in)::coeff(:)
    integer(i32),intent(in)::shape(3)
    type(potcar_t),intent(in)::potcar
    type(crystal_t),intent(in)::crystal
    integer,intent(in)::itype
    type(paw_species_t),intent(inout)::paw
    integer,parameter::ntheta=HALF_QDEP_NTHETA,nphi=HALF_QDEP_NPHI
    integer::nlm,lmax,naug,nrad,natoms,pairs,nang,ich,l,m,a,b,k,ir,iat,ip,it,ion0,alpha
    integer,allocatable::chan(:),lv(:),mv(:),laug(:)
    real(dp),allocatable::rw(:),gshape(:,:),multipole(:,:),dirs(:,:),yweight(:,:),xg(:),wg(:), &
      vlm(:,:),radial_v(:),radial_grad(:,:),pos(:,:),invlat(:,:)
    real(dp)::root1,root2,c1,c2,moment,cart(3),frac(3),value,grad_value(3)
    nlm=0;lmax=0
    do ich=1,potcar%channels
      nlm=nlm+2*potcar%lps(ich)+1;lmax=max(lmax,2*potcar%lps(ich))
    end do
    lmax=active_qdep_lmax(lmax)
    if(nlm/=paw%nlm)error stop 'HALF: CPU QDEP channel mismatch'
    naug=(lmax+1)*(lmax+1);nrad=size(potcar%rgrid);natoms=crystal%counts(itype)
    pairs=nlm*nlm;nang=ntheta*nphi
    allocate(chan(nlm),lv(nlm),mv(nlm));a=0
    do ich=1,potcar%channels
      do m=-potcar%lps(ich),potcar%lps(ich)
        a=a+1;chan(a)=ich;lv(a)=potcar%lps(ich);mv(a)=m
      end do
    end do
    allocate(rw(nrad),gshape(nrad,lmax+1));call radial_weights(potcar%rgrid,rw)
    do l=0,lmax
      call bessel_root(l,1,root1);call bessel_root(l,2,root2)
      call compensation_coefficients(l,potcar%paw_rmax,root1,root2,c1,c2)
      do ir=1,nrad
        if(potcar%rgrid(ir)<=potcar%paw_rmax)then
          gshape(ir,l+1)=c1*sph_bessel(l,root1*potcar%rgrid(ir)/potcar%paw_rmax)+ &
            c2*sph_bessel(l,root2*potcar%rgrid(ir)/potcar%paw_rmax)
        else
          gshape(ir,l+1)=0.0_dp
        end if
      end do
    end do
    allocate(multipole(pairs,naug));multipole=0.0_dp
    do b=1,nlm;do a=1,nlm
      k=0
      do l=0,lmax;do m=-l,l
        k=k+1
        if(abs(lv(a)-lv(b))<=l.and.l<=lv(a)+lv(b).and.mod(lv(a)+lv(b)+l,2)==0)then
          moment=potcar%qpaw_l(chan(a),chan(b),l+1)
          multipole((b-1)*nlm+a,k)=gaunt_numeric(lv(a),mv(a),lv(b),mv(b),l,m)*moment
        end if
      end do;end do
    end do;end do
    allocate(laug(naug),dirs(nang,3),yweight(nang,naug),xg(ntheta),wg(ntheta))
    call gauss_legendre(ntheta,xg,wg);k=0
    do l=0,lmax;do m=-l,l;k=k+1;laug(k)=l;end do;end do
    ip=0
    do it=1,ntheta;do m=0,nphi-1
      ip=ip+1
      dirs(ip,1)=sqrt(max(0.0_dp,1.0_dp-xg(it)*xg(it)))*cos(2*pi*real(m,dp)/real(nphi,dp))
      dirs(ip,2)=sqrt(max(0.0_dp,1.0_dp-xg(it)*xg(it)))*sin(2*pi*real(m,dp)/real(nphi,dp))
      dirs(ip,3)=xg(it);k=0
      do l=0,lmax;do b=-l,l
        k=k+1;yweight(ip,k)=ylm(l,b,xg(it),2*pi*real(m,dp)/real(nphi,dp))*wg(it)*2*pi/real(nphi,dp)
      end do;end do
    end do;end do
    ion0=sum(crystal%counts(:itype-1));allocate(pos(natoms,3),invlat(3,3),vlm(naug,nrad),radial_v(naug),radial_grad(naug,3))
    do iat=1,natoms;pos(iat,:)=matmul(crystal%positions(ion0+iat,:),crystal%lattice);end do
    invlat=inverse3(crystal%lattice)
    do iat=1,natoms
      vlm=0.0_dp
      do ir=1,nrad
        do ip=1,nang
          cart=pos(iat,:)+potcar%rgrid(ir)*dirs(ip,:);frac=modulo(matmul(cart,invlat),1.0_dp)
          value=cubic_sample(coeff,shape,frac)
          do k=1,naug;vlm(k,ir)=vlm(k,ir)+value*yweight(ip,k);end do
        end do
      end do
      do k=1,naug
        radial_v(k)=sum(rw*potcar%rgrid*potcar%rgrid*gshape(:,laug(k)+1)*vlm(k,:))
      end do
      do alpha=1,3
        vlm=0.0_dp
        do ir=1,nrad
          do ip=1,nang
            cart=pos(iat,:)+potcar%rgrid(ir)*dirs(ip,:);frac=modulo(matmul(cart,invlat),1.0_dp)
            call cubic_sample_gradient(coeff,shape,frac,invlat,grad_value);value=grad_value(alpha)
            do k=1,naug;vlm(k,ir)=vlm(k,ir)+value*yweight(ip,k);end do
          end do
        end do
        do k=1,naug
          radial_grad(k,alpha)=sum(rw*potcar%rgrid*potcar%rgrid*gshape(:,laug(k)+1)*vlm(k,:))
        end do
      end do
      paw%dij_atom(iat,:,:)=paw%dij
      do b=1,nlm;do a=1,nlm
        paw%dij_atom(iat,a,b)=paw%dij_atom(iat,a,b)+ &
          sum(multipole((b-1)*nlm+a,:)*radial_v)
      end do;end do
      paw%dij_atom(iat,:,:)=0.5_dp*(paw%dij_atom(iat,:,:)+transpose(paw%dij_atom(iat,:,:)))
      do alpha=1,3;do b=1,nlm;do a=1,nlm
        paw%ddij_atom(iat,a,b,alpha)=sum(multipole((b-1)*nlm+a,:)*radial_grad(:,alpha))
      end do;end do
      paw%ddij_atom(iat,:,:,alpha)=0.5_dp*(paw%ddij_atom(iat,:,:,alpha)+transpose(paw%ddij_atom(iat,:,:,alpha)))
      end do
    end do
  end subroutine

  pure integer function active_qdep_lmax(potcar_lmax) result(lmax)
    integer,intent(in)::potcar_lmax
    lmax=potcar_lmax
#if HALF_QDEP_LMAX >= 0
    lmax=min(lmax,HALF_QDEP_LMAX)
#endif
  end function active_qdep_lmax

  real(dp) function cubic_sample(coeff,shape,frac)result(value)
    real(dp),intent(in)::coeff(:),frac(3)
    integer(i32),intent(in)::shape(3)
    integer::ix,iy,iz,a,b,c,ia,ib,ic,index
    real(dp)::coord(3),w1(4),w2(4),w3(4)
    coord=frac*real(shape,dp);ix=floor(coord(1));iy=floor(coord(2));iz=floor(coord(3))
    call cubic_weights(coord(1)-ix,w1);call cubic_weights(coord(2)-iy,w2);call cubic_weights(coord(3)-iz,w3)
    value=0.0_dp
    do a=1,4;ia=modulo(ix+a-2,shape(1))
      do b=1,4;ib=modulo(iy+b-2,shape(2))
        do c=1,4;ic=modulo(iz+c-2,shape(3));index=ia*shape(2)*shape(3)+ib*shape(3)+ic+1
          value=value+w1(a)*w2(b)*w3(c)*coeff(index)
        end do
      end do
    end do
  end function

  subroutine cubic_sample_gradient(coeff,shape,frac,invlat,gradient)
    real(dp),intent(in)::coeff(:),frac(3),invlat(3,3)
    integer(i32),intent(in)::shape(3)
    real(dp),intent(out)::gradient(3)
    integer::ix,iy,iz,a,b,c,ia,ib,ic,index,alpha
    real(dp)::coord(3),w1(4),w2(4),w3(4),dw1(4),dw2(4),dw3(4),df(3),value
    coord=frac*real(shape,dp);ix=floor(coord(1));iy=floor(coord(2));iz=floor(coord(3))
    call cubic_weights(coord(1)-ix,w1);call cubic_weights(coord(2)-iy,w2);call cubic_weights(coord(3)-iz,w3)
    call cubic_weights_derivative(coord(1)-ix,dw1);call cubic_weights_derivative(coord(2)-iy,dw2)
    call cubic_weights_derivative(coord(3)-iz,dw3);df=0.0_dp
    do a=1,4;ia=modulo(ix+a-2,shape(1));do b=1,4;ib=modulo(iy+b-2,shape(2));do c=1,4
      ic=modulo(iz+c-2,shape(3));index=ia*shape(2)*shape(3)+ib*shape(3)+ic+1;value=coeff(index)
      df(1)=df(1)+dw1(a)*w2(b)*w3(c)*value*real(shape(1),dp)
      df(2)=df(2)+w1(a)*dw2(b)*w3(c)*value*real(shape(2),dp)
      df(3)=df(3)+w1(a)*w2(b)*dw3(c)*value*real(shape(3),dp)
    end do;end do;end do
    do alpha=1,3;gradient(alpha)=sum(df*invlat(alpha,:));end do
  end subroutine cubic_sample_gradient

  pure subroutine cubic_weights(t,w)
    real(dp),intent(in)::t;real(dp),intent(out)::w(4)
    w(1)=(1-t)**3/6;w(2)=(3*t**3-6*t*t+4)/6
    w(3)=(-3*t**3+3*t*t+3*t+1)/6;w(4)=t**3/6
  end subroutine
  pure subroutine cubic_weights_derivative(t,w)
    real(dp),intent(in)::t;real(dp),intent(out)::w(4)
    w(1)=-0.5_dp*(1-t)**2
    w(2)=1.5_dp*t*t-2.0_dp*t
    w(3)=-1.5_dp*t*t+t+0.5_dp
    w(4)=0.5_dp*t*t
  end subroutine cubic_weights_derivative

  subroutine radial_weights(r,w)
    real(dp),intent(in)::r(:);real(dp),intent(out)::w(:)
    integer::i,n,last;real(dp)::h0,h1,s,alpha,beta,eta,hlog,spread
    n=size(r)
    if(n>=3.and.all(r>0.0_dp))then
      hlog=log(r(2)/r(1));spread=maxval(abs(log(r(2:n)/r(1:n-1))-hlog))
      if(spread<1.0e-10_dp*max(1.0_dp,abs(hlog)))then
        ! VASP SET_SIMP: Simpson integration in x=log(r), hence dr=r dx.
        w=0.0_dp
        do i=3,n,2
          w(i)=w(i)+r(i)*hlog/3.0_dp
          w(i-1)=4.0_dp*r(i-1)*hlog/3.0_dp
          w(i-2)=w(i-2)+r(i-2)*hlog/3.0_dp
        end do
        return
      end if
    end if
    w=0.0_dp;last=merge(n,n-1,mod(n,2)==1)
    do i=1,last-2,2
      h0=r(i+1)-r(i);h1=r(i+2)-r(i+1);s=(h0+h1)/6
      w(i)=w(i)+s*(2-h1/h0);w(i+1)=w(i+1)+s*(h0+h1)**2/(h0*h1);w(i+2)=w(i+2)+s*(2-h0/h1)
    end do
    if(mod(n,2)==0)then
      h0=r(n-1)-r(n-2);h1=r(n)-r(n-1);alpha=(2*h1*h1+3*h0*h1)/(6*(h0+h1))
      beta=(h1*h1+3*h0*h1)/(6*h0);eta=h1**3/(6*h0*(h0+h1))
      w(n)=w(n)+alpha;w(n-1)=w(n-1)+beta;w(n-2)=w(n-2)-eta
    end if
  end subroutine

  pure real(dp) function sph_bessel(l,x)result(v)
    integer,intent(in)::l;real(dp),intent(in)::x
    real(dp)::jm,j,jn,den,x2,term;integer::n,i,k
    if(abs(x)<real(l+2,dp))then
      den=1;do i=1,l;den=den*real(2*i+1,dp);end do;x2=x*x;term=x**l/den;v=term
      do k=0,30;term=-term*x2/(2*real(k+1,dp)*real(2*l+2*k+3,dp));v=v+term
        if(abs(term)<1e-16_dp*max(1.0_dp,abs(v)))exit
      end do;return
    end if
    jm=sin(x)/x;if(l==0)then;v=jm;return;end if;j=sin(x)/(x*x)-cos(x)/x
    if(l==1)then;v=j;return;end if
    do n=1,l-1;jn=real(2*n+1,dp)*j/x-jm;jm=j;j=jn;end do;v=j
  end function

  subroutine bessel_root(l,which,root)
    integer,intent(in)::l,which;real(dp),intent(out)::root
    real(dp)::x0,x1,f0,f1;integer::found
    x0=1e-4_dp;f0=sph_bessel(l,x0);found=0
    do while(found<which)
      x1=x0+0.05_dp;f1=sph_bessel(l,x1)
      if(f0*f1<0)then;found=found+1;if(found==which)exit;end if;x0=x1;f0=f1
    end do
    do while(x1-x0>1e-13_dp)
      root=0.5_dp*(x0+x1)
      if(f0*sph_bessel(l,root)<=0)then;x1=root;else;x0=root;f0=sph_bessel(l,x0);end if
    end do;root=0.5_dp*(x0+x1)
  end subroutine

  subroutine compensation_coefficients(l,rc,z1,z2,c1,c2)
    integer,intent(in)::l;real(dp),intent(in)::rc,z1,z2;real(dp),intent(out)::c1,c2
    real(dp)::a11,a12,a21,a22,det,xg(32),wg(32),rg(32)
    a11=z1/rc*sph_bessel_deriv(l,z1);a12=z2/rc*sph_bessel_deriv(l,z2)
    call gauss_legendre(32,xg,wg);rg=0.5_dp*rc*(xg+1);wg=0.5_dp*rc*wg
    a21=sum(wg*sph_bessel_vec(l,z1*rg/rc)*rg**(l+2));a22=sum(wg*sph_bessel_vec(l,z2*rg/rc)*rg**(l+2))
    det=a11*a22-a12*a21;c1=-a12/det;c2=a11/det
  end subroutine

  pure real(dp) function sph_bessel_deriv(l,x)result(v)
    integer,intent(in)::l;real(dp),intent(in)::x
    if(l==0)then;v=-sph_bessel(1,x);else;v=sph_bessel(l-1,x)-real(l+1,dp)/x*sph_bessel(l,x);end if
  end function
  pure function sph_bessel_vec(l,x)result(v)
    integer,intent(in)::l;real(dp),intent(in)::x(:);real(dp)::v(size(x));integer::i
    do i=1,size(x);v(i)=sph_bessel(l,x(i));end do
  end function

  real(dp) function gaunt_numeric(l1,m1,l2,m2,l3,m3)result(g)
    integer,intent(in)::l1,m1,l2,m2,l3,m3;integer,parameter::nt=12,np=24
    real(dp)::xg(nt),wg(nt),phi;integer::it,ip
    call gauss_legendre(nt,xg,wg);g=0
    do it=1,nt;do ip=0,np-1
      phi=2*pi*real(ip,dp)/real(np,dp)
      g=g+wg(it)*2*pi/real(np,dp)*ylm(l1,m1,xg(it),phi)*ylm(l2,m2,xg(it),phi)*ylm(l3,m3,xg(it),phi)
    end do;end do
  end function

  subroutine gauss_legendre(n,x,w)
    integer,intent(in)::n;real(dp),intent(out)::x(n),w(n)
    integer::i,j,m;real(dp)::z,z1,p1,p2,p3,pp
    m=(n+1)/2
    do i=1,m
      z=cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
      do
        p1=1;p2=0;do j=1,n;p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/real(j,dp);end do
        pp=real(n,dp)*(z*p1-p2)/(z*z-1);z1=z;z=z1-p1/pp;if(abs(z-z1)<1e-15_dp)exit
      end do
      x(i)=-z;x(n+1-i)=z;w(i)=2/((1-z*z)*pp*pp);w(n+1-i)=w(i)
    end do
  end subroutine

  pure real(dp) function ylm(l,m,x,phi)result(y)
    integer,intent(in)::l,m;real(dp),intent(in)::x,phi
    real(dp)::p,norm;integer::am
    am=abs(m);p=legendre(l,am,x);norm=sqrt(real(2*l+1,dp)/(4*pi)*factorial(l-am)/factorial(l+am))
    if(m>0)then;y=sqrt(2.0_dp)*(-1.0_dp)**m*norm*p*cos(real(am,dp)*phi)
    else if(m<0)then;y=sqrt(2.0_dp)*(-1.0_dp)**m*norm*p*sin(real(am,dp)*phi)
    else;y=norm*p;end if
  end function

  pure real(dp) function legendre(l,m,x)result(p)
    integer,intent(in)::l,m;real(dp),intent(in)::x
    real(dp)::p0,p1,p2,f;integer::i,n
    p0=1;f=1;do i=1,m;p0=-p0*f*sqrt(max(0.0_dp,1-x*x));f=f+2;end do
    if(l==m)then;p=p0;return;end if;p1=x*real(2*m+1,dp)*p0
    if(l==m+1)then;p=p1;return;end if
    do n=m+2,l;p2=(real(2*n-1,dp)*x*p1-real(n+m-1,dp)*p0)/real(n-m,dp);p0=p1;p1=p2;end do;p=p1
  end function
  pure real(dp) function factorial(n)result(v)
    integer,intent(in)::n;integer::i;v=1;do i=2,n;v=v*i;end do
  end function
end module half_uspp
