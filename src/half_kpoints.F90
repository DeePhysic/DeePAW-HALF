module half_kpoints
  use iso_c_binding,only:c_int,c_double
  use half_kinds,only:dp,i32
  use half_types,only:crystal_t
  implicit none
  private
  public::kpoint_set_t,read_explicit_kpoints,gamma_centered_mesh,gamma_centered_irreducible_mesh,generate_cubic_band_path
  type::kpoint_set_t
    integer(i32)::nk=0
    integer(i32)::divisions(3)=0
    integer(i32)::full_count=0
    real(dp),allocatable::points(:,:),weights(:)
    integer(i32),allocatable::multiplicities(:)
  end type
#ifdef HALF_HAVE_SPGLIB
  interface
    integer(c_int) function spg_get_ir_reciprocal_mesh(grid_address,mapping,mesh,is_shift,time_reversal, &
        lattice,positions,types,natoms,symprec)bind(C,name='spg_get_ir_reciprocal_mesh')
      import::c_int,c_double
      integer(c_int)::grid_address(*),mapping(*)
      integer(c_int),intent(in)::mesh(*),is_shift(*)
      integer(c_int),value::time_reversal,natoms
      real(c_double),intent(in)::lattice(*),positions(*)
      integer(c_int),intent(in)::types(*)
      real(c_double),value::symprec
    end function
  end interface
#endif
contains
  subroutine read_explicit_kpoints(path,set)
    character(len=*),intent(in)::path
    type(kpoint_set_t),intent(out)::set
    integer::unit,ios,i,n
    character(len=1024)::line
    real(dp)::total
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'HALF: cannot open KPOINTS'
    read(unit,'(A)',iostat=ios)line;if(ios/=0)error stop 'HALF: KPOINTS missing title'
    read(unit,*,iostat=ios)n;if(ios/=0.or.n<1)error stop 'HALF: explicit KPOINTS requires a positive count'
    read(unit,'(A)',iostat=ios)line
    if(ios/=0.or.(lower(line(1:1))/='r'.and.lower(line(1:1))/='k')) &
      error stop 'HALF: explicit KPOINTS must use reciprocal coordinates'
    set%nk=n;set%full_count=n;allocate(set%points(n,3),set%weights(n),set%multiplicities(n))
    do i=1,n
      read(unit,*,iostat=ios)set%points(i,:),set%weights(i)
      if(ios/=0.or.set%weights(i)<=0.or.any(set%points(i,:)/=set%points(i,:))) &
        error stop 'HALF: invalid explicit KPOINTS row'
    end do
    close(unit);total=sum(set%weights);set%weights=set%weights/total;set%multiplicities=1
  end subroutine

  subroutine gamma_centered_mesh(crystal,kspacing,set)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(in)::kspacing
    type(kpoint_set_t),intent(out)::set
    integer::i1,i2,i3,k,n1,n2,n3
    if(kspacing<=0.or.kspacing/=kspacing)error stop 'HALF: KSPACING must be positive and finite'
    do k=1,3
      set%divisions(k)=max(1,ceiling(sqrt(sum(crystal%reciprocal(k,:)**2))/kspacing))
    end do
    n1=set%divisions(1);n2=set%divisions(2);n3=set%divisions(3)
    set%nk=n1*n2*n3;set%full_count=set%nk
    allocate(set%points(set%nk,3),set%weights(set%nk),set%multiplicities(set%nk));k=0
    do i1=0,n1-1;do i2=0,n2-1;do i3=0,n3-1
      k=k+1;set%points(k,:)=[wrapped(i1,n1),wrapped(i2,n2),wrapped(i3,n3)]
    end do;end do;end do
    set%weights=1.0_dp/real(set%nk,dp);set%multiplicities=1
  end subroutine

  subroutine gamma_centered_irreducible_mesh(crystal,kspacing,set,symprec,time_reversal)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(in)::kspacing
    type(kpoint_set_t),intent(out)::set
    real(dp),intent(in),optional::symprec
    logical,intent(in),optional::time_reversal
#ifdef HALF_HAVE_SPGLIB
    type(kpoint_set_t)::full
    integer(c_int),allocatable::address(:,:),mapping(:),types(:)
    integer(c_int)::mesh(3),shift(3),use_tr,nat,nir
    real(c_double),allocatable::positions(:,:)
    real(c_double)::lattice(3,3),tol
    integer::i,it,iat,ion,k,mult
    call gamma_centered_mesh(crystal,kspacing,full);mesh=full%divisions;shift=0;nat=crystal%nions
    tol=1e-5_dp;if(present(symprec))tol=symprec
    use_tr=1;if(present(time_reversal))use_tr=merge(1,0,time_reversal)
    allocate(address(3,full%full_count),mapping(full%full_count),positions(3,nat),types(nat))
    lattice=transpose(crystal%lattice);ion=0
    do it=1,crystal%ntypes;do iat=1,crystal%counts(it)
      ion=ion+1;positions(:,ion)=crystal%positions(ion,:);types(ion)=it
    end do;end do
    nir=spg_get_ir_reciprocal_mesh(address,mapping,mesh,shift,use_tr,lattice,positions,types,nat,tol)
    if(nir<=0)error stop 'HALF: spglib irreducible k-mesh failed'
    set%nk=nir;set%full_count=full%full_count;set%divisions=full%divisions
    allocate(set%points(nir,3),set%weights(nir),set%multiplicities(nir));k=0
    do i=0,full%full_count-1
      mult=count(mapping==i)
      if(mult>0)then
        k=k+1;set%points(k,:)=real(address(:,i+1),dp)/real(mesh,dp)
        set%multiplicities(k)=mult;set%weights(k)=real(mult,dp)/real(full%full_count,dp)
      end if
    end do
    if(k/=nir)error stop 'HALF: inconsistent spglib mapping'
#else
    error stop 'HALF: irreducible k mesh requires spglib; set HALF_SPGLIB_ROOT and rebuild'
#endif
  end subroutine

  subroutine generate_cubic_band_path(crystal,npoints,requested_path,set,path_used)
    type(crystal_t),intent(in)::crystal
    integer,intent(in)::npoints
    character(len=*),intent(in)::requested_path
    type(kpoint_set_t),intent(out)::set
    character(len=*),intent(out)::path_used
    real(dp)::gram(3,3),diag_mean,ratio,coordinates(6,3),label_points(len_trim(requested_path)+32,3)
    real(dp)::lengths(len_trim(requested_path)+31),total,x0,remaining,t
    real(dp),allocatable::buffer(:,:)
    character(len=1)::labels(6),c
    character(len=256)::path
    logical::disconnect(len_trim(requested_path)+31),pending
    integer::kind,nlabel,i,j,k,nsegment,n
    if(npoints<2)error stop 'HALF: band-path npoints must be at least two'
    gram=matmul(crystal%lattice,transpose(crystal%lattice));diag_mean=sum([gram(1,1),gram(2,2),gram(3,3)])/3
    if(maxval(abs([gram(1,1),gram(2,2),gram(3,3)]-diag_mean))>1e-6_dp*diag_mean) &
      error stop 'HALF: automatic band paths currently require a cubic primitive cell; use explicit KPOINTS'
    ratio=(gram(1,2)+gram(1,3)+gram(2,3))/(3*diag_mean)
    if(abs(ratio)<1e-6_dp)then
      kind=1;labels(1:4)=['G','X','M','R'];coordinates(1,:)=[0._dp,0._dp,0._dp];coordinates(2,:)=[.5_dp,0._dp,0._dp]
      coordinates(3,:)=[.5_dp,.5_dp,0._dp];coordinates(4,:)=[.5_dp,.5_dp,.5_dp];path='GXMGRX,MR'
    else if(abs(ratio-.5_dp)<1e-6_dp)then
      kind=2;labels=['G','K','L','U','W','X'];coordinates(1,:)=[0._dp,0._dp,0._dp];coordinates(2,:)=[.375_dp,.375_dp,.75_dp]
      coordinates(3,:)=[.5_dp,.5_dp,.5_dp];coordinates(4,:)=[.625_dp,.25_dp,.625_dp]
      coordinates(5,:)=[.5_dp,.25_dp,.75_dp];coordinates(6,:)=[.5_dp,0._dp,.5_dp];path='GXWKGLUWLK,UX'
    else if(abs(ratio+1._dp/3._dp)<1e-6_dp)then
      kind=3;labels(1:4)=['G','H','P','N'];coordinates(1,:)=[0._dp,0._dp,0._dp];coordinates(2,:)=[.5_dp,-.5_dp,.5_dp]
      coordinates(3,:)=[.25_dp,.25_dp,.25_dp];coordinates(4,:)=[0._dp,.5_dp,0._dp];path='GHNGPH,PN'
    else;error stop 'HALF: unsupported cubic primitive metric; use explicit KPOINTS';end if
    if(len_trim(requested_path)>0)path=trim(requested_path);path_used=trim(path)
    nlabel=0;pending=.false.;disconnect=.false.
    do i=1,len_trim(path)
      c=path(i:i)
      if(c==',')then;pending=.true.;cycle;end if
      nlabel=nlabel+1
      if(nlabel>1)disconnect(nlabel-1)=pending
      pending=.false.;k=0
      do j=1,merge(4,6,kind/=2);if(labels(j)==c)then;k=j;exit;end if;end do
      if(k==0)error stop 'HALF: unknown label in cubic band path'
      label_points(nlabel,:)=coordinates(k,:)
    end do
    if(nlabel<2)error stop 'HALF: band path needs at least two labels'
    nsegment=nlabel-1;total=0
    do i=1,nsegment
      if(disconnect(i))then;lengths(i)=0
      else;lengths(i)=sqrt(sum(matmul(label_points(i+1,:)-label_points(i,:),crystal%reciprocal)**2));end if
      total=total+lengths(i)
    end do
    allocate(buffer(npoints+2*nlabel,3));k=0;x0=0
    do i=1,nsegment
      remaining=total-x0
      if(abs(remaining)<1e-12_dp)then;n=0
      else;n=max(2,nint(lengths(i)*real(npoints-k,dp)/remaining));end if
      do j=0,n-2
        t=real(j,dp)/real(n-1,dp);k=k+1;buffer(k,:)=label_points(i,:)+t*(label_points(i+1,:)-label_points(i,:))
      end do
      x0=x0+lengths(i)
    end do
    k=k+1;buffer(k,:)=label_points(nlabel,:);set%nk=k;set%full_count=k
    allocate(set%points(k,3),set%weights(k),set%multiplicities(k));set%points=buffer(:k,:)
    set%weights=1.0_dp/real(k,dp);set%multiplicities=1
  end subroutine

  pure real(dp) function wrapped(index,n)result(value)
    integer,intent(in)::index,n
    value=real(index,dp)/real(n,dp);if(value>0.5_dp)value=value-1.0_dp
  end function
  pure character function lower(c)result(out)
    character,intent(in)::c;integer::v
    v=iachar(c);if(v>=iachar('A').and.v<=iachar('Z'))then;out=achar(v+32);else;out=c;end if
  end function
end module half_kpoints
