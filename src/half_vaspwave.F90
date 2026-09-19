module half_vaspwave
  use iso_c_binding,only:c_int,c_int64_t,c_double,c_float,c_char,c_null_char
  use half_kinds,only:dp
  use half_types,only:crystal_t,charge_grid_t,plane_wave_basis_t
  use half_math,only:inverse3
  implicit none
  private
  public::wave_block_t,write_vaspwave_h5,is_hdf5_file,read_vaspwave_h5,extract_hdf5_potcar
  type::wave_block_t
    complex(dp),allocatable::coefficients(:,:)
  end type
  interface
    integer(c_int) function is_hdf5_c(filename)bind(C,name='half_is_hdf5_c')
      import::c_int,c_char
      character(c_char),intent(in)::filename(*)
    end function
    integer(c_int) function extract_potcar_c(filename,output,output_size)bind(C,name='half_extract_potcar_c')
      import::c_int,c_char
      character(c_char),intent(in)::filename(*)
      character(c_char),intent(out)::output(*)
      integer(c_int),value::output_size
    end function
    integer(c_int) function probe_c(filename,ntypes,nions,grid)bind(C,name='half_probe_vaspwave_c')
      import::c_int,c_char
      character(c_char),intent(in)::filename(*)
      integer(c_int),intent(out)::ntypes,nions,grid(3)
    end function
    integer(c_int) function read_c(filename,ntypes,nions,grid,system,system_size,species,counts,lattice,positions, &
        direct,charge)bind(C,name='half_read_vaspwave_c')
      import::c_int,c_double,c_char
      character(c_char),intent(in)::filename(*)
      integer(c_int),value::ntypes,nions,system_size
      integer(c_int),intent(in)::grid(3)
      character(c_char),intent(out)::system(*),species(*)
      integer(c_int),intent(out)::counts(*),direct
      real(c_double),intent(out)::lattice(*),positions(*),charge(*)
    end function
    integer(c_int) function write_c(filename,system,ntypes,species,counts,nions,lattice,positions,grid,charge, &
        encut,fermi,nk,nb,kpoints,eigenvalues,occupations,npws,offsets,coefficients)bind(C,name='half_write_vaspwave_h5_c')
      import::c_int,c_int64_t,c_double,c_float,c_char
      character(c_char),intent(in)::filename(*),system(*),species(*)
      integer(c_int),value::ntypes,nions,nk,nb
      integer(c_int),intent(in)::counts(*),grid(*),npws(*)
      integer(c_int64_t),intent(in)::offsets(*)
      real(c_double),intent(in)::lattice(*),positions(*),charge(*),kpoints(*),eigenvalues(*),occupations(*)
      real(c_double),value::encut,fermi
      real(c_float),intent(in)::coefficients(*)
    end function
  end interface
contains
  logical function is_hdf5_file(filename)result(value)
    character(len=*),intent(in)::filename
    character(c_char),allocatable::path(:)
    call c_string(filename,path);value=is_hdf5_c(path)/=0
  end function

  subroutine extract_hdf5_potcar(filename,temporary_path)
    character(len=*),intent(in)::filename
    character(len=*),intent(out)::temporary_path
    character(c_char),allocatable::path(:),output(:)
    integer::i,n,status
    call c_string(filename,path);allocate(output(len(temporary_path)+1));output=c_null_char
    status=extract_potcar_c(path,output,int(size(output),c_int));if(status/=0)error stop 'HALF: HDF5 input has no /input/potcar/content'
    n=0;do while(n<size(output).and.output(n+1)/=c_null_char);n=n+1;end do
    temporary_path='';do i=1,min(n,len(temporary_path));temporary_path(i:i)=output(i);end do
  end subroutine

  subroutine read_vaspwave_h5(filename,crystal,charge)
    character(len=*),intent(in)::filename
    type(crystal_t),intent(out)::crystal
    type(charge_grid_t),intent(out)::charge
    character(c_char),allocatable::path(:),system(:),species(:)
    integer(c_int)::ntypes,nions,grid(3),direct,status
    integer(c_int),allocatable::counts(:)
    real(c_double),allocatable::lattice(:,:),positions(:,:)
    integer::i,j,n
    call c_string(filename,path);status=probe_c(path,ntypes,nions,grid)
    if(status/=0)error stop 'HALF: invalid vaspwave.h5 charge/structure input'
    allocate(system(1024),species(16*ntypes),counts(ntypes),lattice(3,3),positions(3,nions))
    charge%shape=int(grid);allocate(charge%values(int(grid(1))*int(grid(2))*int(grid(3))))
    status=read_c(path,ntypes,nions,grid,system,1024_c_int,species,counts,lattice,positions,direct,charge%values)
    if(status/=0)error stop 'HALF: cannot read vaspwave.h5 charge/structure input'
    n=0;do while(n<size(system).and.system(n+1)/=c_null_char);n=n+1;end do
    allocate(character(len=max(1,n))::crystal%system_name);do i=1,n;crystal%system_name(i:i)=system(i);end do
    crystal%ntypes=ntypes;crystal%nions=nions;allocate(crystal%species(ntypes),crystal%counts(ntypes),crystal%positions(nions,3))
    crystal%counts=int(counts);crystal%lattice=transpose(lattice);crystal%positions=transpose(positions)
    do i=1,ntypes
      crystal%species(i)='';do j=1,16;crystal%species(i)(j:j)=species(16*(i-1)+j);end do
    end do
    call crystal%update_geometry()
    if(direct==0)crystal%positions=matmul(crystal%positions,inverse3(crystal%lattice))
  end subroutine

  subroutine write_vaspwave_h5(filename,crystal,charge,encut,kpoints,eigenvalues,occupations,fermi,bases,waves)
    character(len=*),intent(in)::filename
    type(crystal_t),intent(in)::crystal
    type(charge_grid_t),intent(in)::charge
    real(dp),intent(in)::encut,kpoints(:,:),eigenvalues(:,:),occupations(:,:),fermi
    type(plane_wave_basis_t),intent(in)::bases(:)
    type(wave_block_t),intent(in)::waves(:)
    character(c_char),allocatable::cfilename(:),csystem(:),species(:)
    integer(c_int),allocatable::counts(:),grid(:),npws(:)
    integer(c_int64_t),allocatable::offsets(:)
    integer,allocatable::order(:)
    real(c_double),allocatable::lattice(:,:),positions(:,:),kp(:,:),eig(:,:),occ(:,:)
    real(c_float),allocatable::packed(:)
    integer::nk,nb,ik,ib,ig,npw,total,cursor,j,status
    nk=size(kpoints,1);nb=size(eigenvalues,2)
    if(size(kpoints,2)/=3.or.size(eigenvalues,1)/=nk.or.any(shape(occupations)/=shape(eigenvalues))) &
      error stop 'HALF: HDF5 band array mismatch'
    if(size(bases)/=nk.or.size(waves)/=nk)error stop 'HALF: HDF5 wave block count mismatch'
    if(any(occupations< -1e-10_dp).or.any(occupations>2.0_dp+1e-10_dp))error stop 'HALF: invalid HDF5 occupations'
    call c_string(filename,cfilename);call c_string(crystal%system_name,csystem)
    allocate(species(16*crystal%ntypes));species=' '
    do ik=1,crystal%ntypes;do j=1,min(16,len_trim(crystal%species(ik)))
      species(16*(ik-1)+j)=crystal%species(ik)(j:j)
    end do;end do
    counts=int(crystal%counts,c_int);grid=int(charge%shape,c_int);allocate(npws(nk),offsets(nk))
    total=0
    do ik=1,nk;npws(ik)=int(bases(ik)%npw,c_int);offsets(ik)=int(total,c_int64_t);total=total+2*npws(ik)*nb;end do
    allocate(packed(total));cursor=0
    do ik=1,nk
      npw=npws(ik)
      if(.not.allocated(waves(ik)%coefficients))error stop 'HALF: missing HDF5 coefficients'
      if(any(shape(waves(ik)%coefficients)/=[npw,nb]))error stop 'HALF: HDF5 coefficient shape mismatch'
      call vasp_permutation(bases(ik),order)
      do ib=1,nb;do ig=1,npw
        cursor=cursor+1;packed(cursor)=real(waves(ik)%coefficients(order(ig),ib),c_float)
        cursor=cursor+1;packed(cursor)=real(aimag(waves(ik)%coefficients(order(ig),ib)),c_float)
      end do;end do
      deallocate(order)
    end do
    allocate(lattice(3,3),positions(3,crystal%nions),kp(3,nk),eig(nb,nk),occ(nb,nk))
    lattice=transpose(crystal%lattice);positions=transpose(crystal%positions);kp=transpose(kpoints)
    eig=transpose(eigenvalues);occ=transpose(occupations)
    status=write_c(cfilename,csystem,int(crystal%ntypes,c_int),species,counts,int(crystal%nions,c_int), &
      lattice,positions,grid,charge%values,encut,fermi,int(nk,c_int),int(nb,c_int),kp,eig,occ,npws,offsets,packed)
    if(status/=0)error stop 'HALF: failed to write vaspwave.h5'
  end subroutine

  subroutine vasp_permutation(basis,order)
    type(plane_wave_basis_t),intent(in)::basis
    integer,allocatable,intent(out)::order(:)
    integer,allocatable::lookup(:)
    integer::gx,gy,gz,nx,ny,nz,index0,k,i
    nx=basis%shape(1);ny=basis%shape(2);nz=basis%shape(3);allocate(lookup(0:nx*ny*nz-1));lookup=0
    do i=1,basis%npw;lookup(basis%fft_index(i)-1)=i;end do
    allocate(order(basis%npw));k=0
    do gz=0,nz-1;do gy=0,ny-1;do gx=0,nx-1
      index0=gx*ny*nz+gy*nz+gz
      if(lookup(index0)>0)then;k=k+1;order(k)=lookup(index0);end if
    end do;end do;end do
    if(k/=basis%npw)error stop 'HALF: incomplete VASP plane-wave permutation'
  end subroutine

  subroutine c_string(input,output)
    character(len=*),intent(in)::input
    character(c_char),allocatable,intent(out)::output(:)
    integer::i,n
    n=len_trim(input);allocate(output(n+1));do i=1,n;output(i)=input(i:i);end do;output(n+1)=c_null_char
  end subroutine
end module half_vaspwave
