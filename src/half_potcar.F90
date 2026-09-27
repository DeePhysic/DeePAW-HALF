module half_potcar
  use half_kinds, only: dp, i32
  use half_types, only: potcar_t,crystal_t
#ifdef HALF_HAVE_HDF5
  use half_vaspwave,only:is_hdf5_file,extract_hdf5_potcar
#endif
  implicit none
  private
  public :: read_potcar,validate_potcar_structure
  integer, parameter :: npspts=1000, npsnl=100, npsrnl=100
contains
  subroutine validate_potcar_structure(datasets,crystal)
    type(potcar_t),intent(in)::datasets(:)
    type(crystal_t),intent(in)::crystal
    integer::i
    if(size(datasets)/=crystal%ntypes)error stop 'HALF: POTCAR dataset count does not match structure species count'
    do i=1,crystal%ntypes
      if(index(lower(trim(crystal%species(i))),'type')==1)cycle
      if(trim(lower(datasets(i)%element))/=trim(lower(crystal%species(i)))) &
        error stop 'HALF: POTCAR dataset order does not match structure species order'
    end do
  end subroutine

  subroutine read_potcar(path,datasets)
    character(len=*),intent(in)::path
    type(potcar_t),allocatable,intent(out)::datasets(:)
    type(potcar_t),allocatable::work(:)
    type(potcar_t)::item
    character(len=1024)::line,read_path
    integer::unit,ios,n
    logical::temporary
    read_path=path;temporary=.false.
#ifdef HALF_HAVE_HDF5
    if(is_hdf5_file(path))then;call extract_hdf5_potcar(path,read_path);temporary=.true.;end if
#endif
    open(newunit=unit,file=trim(read_path),status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'HALF: cannot open POTCAR'
    allocate(work(32)); n=0
    do
      call next_nonblank(unit,line,ios)
      if(ios/=0)exit
      if(index(line,'End of Dataset')>0)cycle
      call read_dataset(unit,trim(adjustl(line)),item)
      n=n+1; if(n>size(work))error stop 'HALF: too many POTCAR datasets'
      work(n)=item
    end do
    if(temporary)then;close(unit,status='delete');else;close(unit);end if
    allocate(datasets(n)); if(n>0)datasets=work(:n)
  end subroutine

  subroutine read_dataset(unit,header,p)
    integer,intent(in)::unit
    character(len=*),intent(in)::header
    type(potcar_t),intent(out)::p
    character(len=1024)::line
    integer::ios
    p%header=header; call element_from_header(header,p%element)
    read(unit,*,iostat=ios)p%zval; call require(ios,'POTCAR ZVAL')
    call read_psctr(unit,p)
    call next_nonblank(unit,line,ios); call require(ios,'POTCAR local header')
    if(index(lower(line),'local part')==0)error stop 'HALF: expected POTCAR local part'
    call read_local(unit,p); call read_projectors_and_paw(unit,p)
  end subroutine

  subroutine read_psctr(unit,p)
    integer,intent(in)::unit
    type(potcar_t),intent(inout)::p
    character(len=1024)::line
    integer::ios
    do
      read(unit,'(A)',iostat=ios)line; call require(ios,'POTCAR PSCTR')
      if(index(lower(line),'end of psctr')>0)exit
      call real_tag(line,'ENMAX',p%enmax); call real_tag(line,'POMASS',p%pomass)
      call real_tag(line,'ZVAL',p%zval); call real_tag(line,'EAUG',p%eaug)
      call real_tag(line,'EATOM',p%eatom); call real_tag(line,'DEXC',p%dexccore)
      call string_tag(line,'LEXCH',p%lexch); call string_tag(line,'VRHFIN',p%element)
    end do
  end subroutine

  subroutine read_local(unit,p)
    integer,intent(in)::unit
    type(potcar_t),intent(inout)::p
    character(len=1024)::line
    integer::ios,i
    real(dp),allocatable::scratch(:)
    read(unit,*,iostat=ios)p%psp_gmax; call require(ios,'POTCAR local grid')
    allocate(p%psp_local(npspts),p%psp_g_grid(npspts))
    read(unit,*,iostat=ios)p%psp_local; call require(ios,'POTCAR local potential')
    do i=1,npspts; p%psp_g_grid(i)=p%psp_gmax*real(i-1,dp)/real(npspts,dp); end do
    do
      call next_nonblank(unit,line,ios); call require(ios,'POTCAR local optional data')
      if(index(lower(line),'gradient')>0)then
        read(unit,'(A)',iostat=ios)line; call require(ios,'POTCAR gradient selector')
      else if(index(lower(line),'core charge-density (partial)')>0)then
        allocate(p%pspcor(npspts)); read(unit,*,iostat=ios)p%pspcor
        call require(ios,'POTCAR partial core charge'); p%has_core=.true.
      else if(index(lower(line),'kinetic energy density')>0)then
        allocate(scratch(npspts)); read(unit,*,iostat=ios)scratch
        call require(ios,'POTCAR kinetic density'); deallocate(scratch)
      else if(index(lower(line),'atomic pseudo charge-density')>0)then
        allocate(p%psprho(npspts)); read(unit,*,iostat=ios)p%psprho
        call require(ios,'POTCAR atomic pseudo charge'); exit
      else
        error stop 'HALF: unknown POTCAR local-data selector'
      end if
    end do
  end subroutine

  subroutine read_projectors_and_paw(unit,p)
    integer,intent(in)::unit
    type(potcar_t),intent(inout)::p
    character(len=1024)::line
    integer::ios,lval,npro,i,j,oldn,newn
    real(dp)::rmax
    integer(i32),allocatable::ltmp(:)
    real(dp),allocatable::dtmp(:,:),gtmp(:,:),rtmp(:,:),maxtmp(:),flat(:)
    read(unit,*,iostat=ios)p%pspnl_gmax; call require(ios,'POTCAR projector grid')
    allocate(p%lps(0),p%dion(0,0),p%pspnl(npsnl,0),p%psprnl(npsrnl,0),p%pspnl_rmax(0))
    do
      call next_nonblank(unit,line,ios); call require(ios,'POTCAR projector header')
      if(index(lower(line),'paw radial sets')>0)exit
      if(index(lower(line),'non local part')==0)error stop 'HALF: expected non-local projector block'
      read(unit,*,iostat=ios)lval,npro,rmax; call require(ios,'POTCAR projector channels')
      oldn=p%channels; newn=oldn+npro
      allocate(ltmp(newn),dtmp(newn,newn),gtmp(npsnl,newn),rtmp(npsrnl,newn),maxtmp(newn)); dtmp=0
      if(oldn>0)then
        ltmp(:oldn)=p%lps; dtmp(:oldn,:oldn)=p%dion; gtmp(:,:oldn)=p%pspnl
        rtmp(:,:oldn)=p%psprnl; maxtmp(:oldn)=p%pspnl_rmax
      end if
      ltmp(oldn+1:newn)=lval; maxtmp(oldn+1:newn)=rmax
      allocate(flat(npro*npro)); read(unit,*,iostat=ios)flat; call require(ios,'POTCAR DION')
      do i=1,npro; do j=1,npro; dtmp(oldn+i,oldn+j)=flat((i-1)*npro+j); end do; end do
      deallocate(flat)
      do i=oldn+1,newn
        call next_nonblank(unit,line,ios); call require(ios,'POTCAR reciprocal header')
        if(index(lower(line),'reciprocal space part')==0)error stop 'HALF: expected reciprocal projector'
        read(unit,*,iostat=ios)gtmp(:,i); call require(ios,'POTCAR reciprocal projector')
        call next_nonblank(unit,line,ios); call require(ios,'POTCAR real header')
        if(index(lower(line),'real space part')==0)error stop 'HALF: expected real projector'
        read(unit,*,iostat=ios)rtmp(:,i); call require(ios,'POTCAR real projector')
      end do
      call move_alloc(ltmp,p%lps); call move_alloc(dtmp,p%dion); call move_alloc(gtmp,p%pspnl)
      call move_alloc(rtmp,p%psprnl); call move_alloc(maxtmp,p%pspnl_rmax); p%channels=newn
    end do
    call read_paw(unit,p)
  end subroutine

  subroutine read_paw(unit,p)
    integer,intent(in)::unit
    type(potcar_t),intent(inout)::p
    character(len=1024)::line,key
    integer::ios,i,n
    real(dp),allocatable::flat(:),scratch(:)
    read(unit,*,iostat=ios)p%nmax,p%paw_rmax; call require(ios,'POTCAR PAW grid header')
    read(unit,'(A)',iostat=ios)line; call require(ios,'POTCAR PAW format')
    call next_nonblank(unit,line,ios); call require(ios,'POTCAR augmentation header')
    if(index(lower(line),'augmentation charges')==0)error stop 'HALF: expected augmentation charges'
    n=p%channels; allocate(flat(n*n),p%qpaw(n,n),p%qato(n,n)); p%qato=0
    read(unit,*,iostat=ios)flat; call require(ios,'POTCAR QPAW'); call row_major(flat,p%qpaw); deallocate(flat)
    allocate(p%rgrid(p%nmax),p%potae(p%nmax),p%potps(p%nmax),p%potpsc(p%nmax))
    allocate(p%rhoae(p%nmax),p%rhops(p%nmax),p%wae(p%nmax,n),p%wps(p%nmax,n))
    p%rgrid=0; p%potae=0; p%potps=0; p%potpsc=0; p%rhoae=0; p%rhops=0; p%wae=0; p%wps=0
    do
      call next_nonblank(unit,line,ios); call require(ios,'POTCAR PAW selector'); key=lower(adjustl(line))
      if(index(key,'occupancies in atom')>0.or.index(key,'uccopancies in atom')>0)then
        allocate(flat(n*n)); read(unit,*,iostat=ios)flat; call require(ios,'POTCAR QATO')
        call row_major(flat,p%qato); deallocate(flat)
      else if(trim(key)=='grid')then
        read(unit,*,iostat=ios)p%rgrid; call require(ios,'POTCAR radial grid')
      else if(trim(key)=='aepotential')then
        read(unit,*,iostat=ios)p%potae; call require(ios,'POTCAR AE potential')
      else if(index(key,'core charge-density (pseudized)')>0)then
        read(unit,*,iostat=ios)p%rhops; call require(ios,'POTCAR pseudo core charge')
      else if(index(key,'core charge-density')>0)then
        read(unit,*,iostat=ios)p%rhoae; call require(ios,'POTCAR AE core charge')
      else if(index(key,'local pseudopotential core')>0)then
        read(unit,*,iostat=ios)p%potpsc; call require(ios,'POTCAR core potential')
      else if(index(key,'pspotential valence only')>0)then
        read(unit,*,iostat=ios)p%potps; call require(ios,'POTCAR valence potential')
      else if(index(key,'kinetic energy-density')>0.or.index(key,'mkinetic energy-density')>0)then
        allocate(scratch(p%nmax)); read(unit,*,iostat=ios)scratch
        call require(ios,'POTCAR kinetic radial data'); deallocate(scratch)
      else if(index(key,'pseudo wavefunction')>0)then
        read(unit,*,iostat=ios)p%wps(:,1); call require(ios,'POTCAR pseudo wavefunction')
        call next_nonblank(unit,line,ios); call require(ios,'POTCAR AE wavefunction header')
        read(unit,*,iostat=ios)p%wae(:,1); call require(ios,'POTCAR AE wavefunction')
        do i=2,n
          call next_nonblank(unit,line,ios); call require(ios,'POTCAR pseudo wavefunction header')
          read(unit,*,iostat=ios)p%wps(:,i); call require(ios,'POTCAR pseudo wavefunction')
          call next_nonblank(unit,line,ios); call require(ios,'POTCAR AE wavefunction header')
          read(unit,*,iostat=ios)p%wae(:,i); call require(ios,'POTCAR AE wavefunction')
        end do
        exit
      else if(key(1:1)=='t')then
        allocate(scratch(n*n)); read(unit,*,iostat=ios)scratch
        call require(ios,'POTCAR total charge'); deallocate(scratch)
      else
        error stop 'HALF: unknown POTCAR PAW selector'
      end if
    end do
    do
      read(unit,'(A)',iostat=ios)line
      if(ios/=0.or.index(line,'End of Dataset')>0)exit
    end do
    call rebuild_paw_qion(p)
  end subroutine

  subroutine rebuild_paw_qion(p)
    ! VASP SET_PAW_AUG does not retain the tabulated PAW overlap matrix as
    ! the final QION.  It rebuilds QPAW(i,j,L) from the AE/PS partial waves
    ! with the POTCAR logarithmic-grid Simpson weights and assigns
    ! QION(i,j)=QPAW(i,j,0).  Keep the same one-centre definition here so
    ! H, S, augmentation density, and forces share identical L=0 moments.
    type(potcar_t),intent(inout)::p
    real(dp),allocatable::w(:),rebuilt(:,:)
    real(dp)::hlog,spread
    integer::i,j,k,l,lmin,lmax,n
    n=size(p%rgrid);allocate(w(n),rebuilt(p%channels,p%channels));w=0.0_dp;rebuilt=0.0_dp
    if(n<3.or.any(p%rgrid<=0.0_dp))error stop 'HALF: invalid PAW logarithmic radial grid'
    hlog=log(p%rgrid(2)/p%rgrid(1));spread=maxval(abs(log(p%rgrid(2:n)/p%rgrid(1:n-1))-hlog))
    if(spread>1.0e-10_dp*max(1.0_dp,abs(hlog)))error stop 'HALF: PAW radial grid is not logarithmic'
    do k=3,n,2
      w(k)=w(k)+p%rgrid(k)*hlog/3.0_dp
      w(k-1)=4.0_dp*p%rgrid(k-1)*hlog/3.0_dp
      w(k-2)=w(k-2)+p%rgrid(k-2)*hlog/3.0_dp
    end do
    p%qpaw_rebuild_max_delta=0.0_dp
    lmax=2*maxval(p%lps);allocate(p%qpaw_l(p%channels,p%channels,lmax+1));p%qpaw_l=0.0_dp
    do j=1,p%channels;do i=1,p%channels
      lmin=abs(p%lps(i)-p%lps(j))
      do l=lmin,p%lps(i)+p%lps(j),2
        p%qpaw_l(i,j,l+1)=sum(w*(p%wae(:,i)*p%wae(:,j)-p%wps(:,i)*p%wps(:,j))*p%rgrid**l)
      end do
      if(p%lps(i)==p%lps(j))then
        rebuilt(i,j)=p%qpaw_l(i,j,1)
        p%qpaw_rebuild_max_delta=max(p%qpaw_rebuild_max_delta,abs(rebuilt(i,j)-p%qpaw(i,j)))
      end if
    end do;end do
    p%qpaw=rebuilt
  end subroutine rebuild_paw_qion

  subroutine row_major(flat,matrix)
    real(dp),intent(in)::flat(:)
    real(dp),intent(out)::matrix(:,:)
    integer::i,j,nc
    nc=size(matrix,2)
    do i=1,size(matrix,1); do j=1,nc; matrix(i,j)=flat((i-1)*nc+j); end do; end do
  end subroutine
  subroutine next_nonblank(unit,line,ios)
    integer,intent(in)::unit
    character(len=*),intent(out)::line
    integer,intent(out)::ios
    do; read(unit,'(A)',iostat=ios)line; if(ios/=0.or.len_trim(line)>0)return; end do
  end subroutine
  subroutine real_tag(line,tag,value)
    character(len=*),intent(in)::line,tag
    real(dp),intent(inout)::value
    character(len=1024)::part
    integer::at,eq,semi,ios
    at=index(line,tag); if(at==0)return; eq=index(line(at:),'='); if(eq==0)return
    part=adjustl(line(at+eq:)); semi=index(part,';'); if(semi>0)part=part(:semi-1)
    read(part,*,iostat=ios)value
  end subroutine
  subroutine string_tag(line,tag,value)
    character(len=*),intent(in)::line,tag
    character(len=*),intent(inout)::value
    character(len=1024)::part
    integer::at,eq,last
    at=index(line,tag); if(at==0)return; eq=index(line(at:),'='); if(eq==0)return
    part=adjustl(line(at+eq:)); last=scan(part,' :;'//achar(9)); if(last==0)last=len_trim(part)+1
    value=part(:last-1)
  end subroutine
  subroutine element_from_header(line,element)
    character(len=*),intent(in)::line
    character(len=*),intent(out)::element
    character(len=64)::a,b
    integer::ios
    read(line,*,iostat=ios)a,b; if(ios==0)then; element=b; else; element=''; end if
  end subroutine
  pure function lower(text)result(out)
    character(len=*),intent(in)::text
    character(len=len(text))::out
    integer::i,c
    do i=1,len(text); c=iachar(text(i:i)); if(c>=65.and.c<=90)then; out(i:i)=achar(c+32); else; out(i:i)=text(i:i); end if; end do
  end function
  subroutine require(ios,what)
    integer,intent(in)::ios
    character(len=*),intent(in)::what
    if(ios/=0)then; write(*,'(A)')'HALF: failed to read '//trim(what); error stop; end if
  end subroutine
end module half_potcar
