module half_kpoints
  use iso_c_binding,only:c_int,c_double
  use half_kinds,only:dp,i32
  use half_types,only:crystal_t
  use half_math,only:inverse3
  implicit none
  private
  public::kpoint_set_t,read_explicit_kpoints,gamma_centered_mesh,gamma_centered_irreducible_mesh,generate_cubic_band_path, &
    symmetrize_scalar_grid,symmetrize_atomic_vectors
  type::kpoint_set_t
    integer(i32)::nk=0
    integer(i32)::divisions(3)=0
    integer(i32)::full_count=0
    real(dp),allocatable::points(:,:),weights(:)
    integer(i32),allocatable::multiplicities(:)
    integer(i32)::nsym=1
    integer(i32),allocatable::rotations(:,:,:),atom_map(:,:)
    real(dp),allocatable::translations(:,:),cart_rotations(:,:,:)
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
    integer(c_int) function spg_get_symmetry(rotation,translation,max_size,lattice,positions,types,natoms,symprec) &
        bind(C,name='spg_get_symmetry')
      import::c_int,c_double
      integer(c_int)::rotation(*)
      real(c_double)::translation(*)
      integer(c_int),value::max_size,natoms
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
    call identity_symmetry(set,0)
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
    call identity_symmetry(set,crystal%nions)
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
    call load_space_group(crystal,tol,set)
#else
    error stop 'HALF: irreducible k mesh requires spglib; set HALF_SPGLIB_ROOT and rebuild'
#endif
  end subroutine

  subroutine identity_symmetry(set,natoms)
    type(kpoint_set_t),intent(inout)::set
    integer,intent(in)::natoms
    integer::i
    set%nsym=1
    allocate(set%rotations(3,3,1),set%translations(3,1),set%cart_rotations(3,3,1),set%atom_map(max(0,natoms),1))
    set%rotations=0;set%translations=0.0_dp;set%cart_rotations=0.0_dp
    do i=1,3;set%rotations(i,i,1)=1;set%cart_rotations(i,i,1)=1.0_dp;end do
    do i=1,natoms;set%atom_map(i,1)=i;end do
  end subroutine identity_symmetry

#ifdef HALF_HAVE_SPGLIB
  subroutine load_space_group(crystal,tol,set)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(in)::tol
    type(kpoint_set_t),intent(inout)::set
    integer,parameter::maxsym=384
    integer(c_int),allocatable::raw_rotation(:,:,:),types(:)
    real(c_double),allocatable::raw_translation(:,:),positions(:,:)
    real(c_double)::lattice(3,3)
    real(dp)::rmat(3,3),cart(3,3),target(3),delta(3),best
    integer(c_int)::nsym
    integer::op,i,j,it,iat,ion,bestj
    allocate(raw_rotation(3,3,maxsym),raw_translation(3,maxsym),positions(3,crystal%nions),types(crystal%nions))
    lattice=transpose(crystal%lattice);ion=0
    do it=1,crystal%ntypes;do iat=1,crystal%counts(it)
      ion=ion+1;positions(:,ion)=crystal%positions(ion,:);types(ion)=it
    end do;end do
    nsym=spg_get_symmetry(raw_rotation,raw_translation,maxsym,lattice,positions,types,crystal%nions,tol)
    if(nsym<=0.or.nsym>maxsym)error stop 'HALF: spglib symmetry-operation enumeration failed'
    set%nsym=nsym
    allocate(set%rotations(3,3,nsym),set%translations(3,nsym),set%cart_rotations(3,3,nsym), &
      set%atom_map(crystal%nions,nsym))
    do op=1,nsym
      ! C stores rotation[operation][row][column]; the first two indices are
      ! reversed when the same memory is viewed as a Fortran array.
      set%rotations(:,:,op)=transpose(raw_rotation(:,:,op))
      set%translations(:,op)=raw_translation(:,op)
      rmat=real(set%rotations(:,:,op),dp)
      cart=matmul(transpose(crystal%lattice),matmul(rmat,inverse3(transpose(crystal%lattice))))
      set%cart_rotations(:,:,op)=cart
      do i=1,crystal%nions
        target=matmul(rmat,crystal%positions(i,:))+set%translations(:,op)
        best=huge(1.0_dp);bestj=0
        do j=1,crystal%nions
          if(types(j)/=types(i))cycle
          delta=target-crystal%positions(j,:);delta=delta-anint(delta)
          if(sqrt(sum(matmul(delta,crystal%lattice)**2))<best)then
            best=sqrt(sum(matmul(delta,crystal%lattice)**2));bestj=j
          end if
        end do
        if(bestj==0.or.best>10.0_dp*tol)error stop 'HALF: symmetry operation does not map the atomic structure'
        set%atom_map(i,op)=bestj
      end do
    end do
  end subroutine load_space_group
#endif

  subroutine symmetrize_scalar_grid(set,shape,values)
    type(kpoint_set_t),intent(in)::set
    integer,intent(in)::shape(3)
    real(dp),intent(inout)::values(:)
    real(dp),allocatable::result(:)
    real(dp)::frac(3),target(3),scaled(3)
    integer::op,i1,i2,i3,j1,j2,j3,source,destination
    if(size(values)/=product(shape))error stop 'HALF: symmetry grid size mismatch'
    if(set%nsym<=1)return
    allocate(result(size(values)));result=0.0_dp
    do op=1,set%nsym
      do i3=0,shape(3)-1;do i2=0,shape(2)-1;do i1=0,shape(1)-1
        frac=[real(i1,dp)/shape(1),real(i2,dp)/shape(2),real(i3,dp)/shape(3)]
        target=matmul(real(set%rotations(:,:,op),dp),frac)+set%translations(:,op)
        target=target-floor(target);scaled=target*real(shape,dp)
        if(maxval(abs(scaled-anint(scaled)))>1.0e-7_dp) &
          error stop 'HALF: FFT grid is incompatible with a crystal symmetry translation'
        j1=modulo(nint(scaled(1)),shape(1));j2=modulo(nint(scaled(2)),shape(2));j3=modulo(nint(scaled(3)),shape(3))
        source=i1+shape(1)*(i2+shape(2)*i3)+1
        destination=j1+shape(1)*(j2+shape(2)*j3)+1
        result(destination)=result(destination)+values(source)/real(set%nsym,dp)
      end do;end do;end do
    end do
    values=result
  end subroutine symmetrize_scalar_grid

  subroutine symmetrize_atomic_vectors(set,vectors)
    type(kpoint_set_t),intent(in)::set
    real(dp),intent(inout)::vectors(:,:)
    real(dp),allocatable::result(:,:)
    integer::op,i,j
    if(set%nsym<=1)return
    if(size(vectors,2)/=3.or.size(vectors,1)/=size(set%atom_map,1)) &
      error stop 'HALF: symmetry atomic-vector size mismatch'
    allocate(result(size(vectors,1),3));result=0.0_dp
    do op=1,set%nsym;do i=1,size(vectors,1)
      j=set%atom_map(i,op)
      result(j,:)=result(j,:)+matmul(set%cart_rotations(:,:,op),vectors(i,:))/real(set%nsym,dp)
    end do;end do
    vectors=result
  end subroutine symmetrize_atomic_vectors

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
