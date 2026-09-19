module half_c_api
  use iso_c_binding,only:c_int,c_int32_t,c_int64_t,c_double,c_double_complex,c_char,c_null_char,c_ptr,c_loc, &
    c_f_pointer,c_associated
  use half_kinds,only:dp
  use half_types,only:plane_wave_basis_t
  use half_library,only:half_context_t,LIB_SUCCESS=>HALF_SUCCESS,LIB_INVALID=>HALF_ERROR_INVALID_ARGUMENT, &
    HALF_BACKEND_AUTO
  implicit none
  private
  integer(c_int),parameter::HALF_SUCCESS=0,HALF_INVALID_ARGUMENT=1,HALF_IO=2,HALF_INVALID_HANDLE=3, &
    HALF_CAPACITY=4,HALF_UNAVAILABLE=5,HALF_INTERNAL=6
  integer(c_int),parameter::MAX_CONTEXTS=32
  type(half_context_t),save::contexts(MAX_CONTEXTS)
  logical,save::used(MAX_CONTEXTS)=.false.
  integer(c_int64_t),save::generation(MAX_CONTEXTS)=0_c_int64_t
  character(c_char),target,save::version(6)=[character(c_char)::'0','.', '5','.', '0',c_null_char]
contains
  integer(c_int) function half_get_abi_version()bind(C,name='half_get_abi_version')
    half_get_abi_version=1
  end function

  type(c_ptr) function half_get_version_string()bind(C,name='half_get_version_string')
    half_get_version_string=c_loc(version(1))
  end function

  integer(c_int) function half_get_capabilities()bind(C,name='half_get_capabilities')
    half_get_capabilities=0
#ifdef HALF_HAVE_MKL
    half_get_capabilities=ior(half_get_capabilities,1)
#endif
#ifdef HALF_HAVE_CUDA
    half_get_capabilities=ior(half_get_capabilities,2)
#endif
#ifdef HALF_HAVE_MPI
    half_get_capabilities=ior(half_get_capabilities,4)
#endif
#ifdef HALF_HAVE_HDF5
    half_get_capabilities=ior(half_get_capabilities,8)
#endif
#ifdef HALF_HAVE_SPGLIB
    half_get_capabilities=ior(half_get_capabilities,16)
#endif
  end function

  integer(c_int) function half_create_from_files(charge_c,potential_c,encut,xc,use_uspp,handle,error,error_capacity) &
      bind(C,name='half_create_from_files')
    character(c_char),intent(in)::charge_c(*),potential_c(*)
    real(c_double),value::encut
    integer(c_int),value::xc,use_uspp,error_capacity
    integer(c_int64_t),intent(out)::handle
    character(c_char),intent(out)::error(*)
    character(len=:),allocatable::charge,potential
    character(len=512)::message
    character(len=3)::xc_name
    integer::slot,status
    handle=0;call clear_error(error,error_capacity)
    call import_string(charge_c,charge);call import_string(potential_c,potential)
    if(len(charge)==0.or.len(potential)==0.or.(xc/=1.and.xc/=2))then
      half_create_from_files=HALF_INVALID_ARGUMENT;call export_error('invalid path or XC selector',error,error_capacity);return
    end if
    slot=0
    do status=1,MAX_CONTEXTS;if(.not.used(status))then;slot=status;exit;end if;end do
    if(slot==0)then;half_create_from_files=HALF_INTERNAL;call export_error('HALF context registry is full',error,error_capacity);return;end if
    xc_name=merge('pbe','lda',xc==2)
    call contexts(slot)%initialize_files(charge,potential,real(encut,dp),xc_name,use_uspp/=0,status,message,HALF_BACKEND_AUTO,'evd')
    if(status/=LIB_SUCCESS)then
      half_create_from_files=int(status,c_int);call export_error(trim(message),error,error_capacity);return
    end if
    used(slot)=.true.;generation(slot)=generation(slot)+1
    handle=generation(slot)*int(MAX_CONTEXTS,c_int64_t)+int(slot,c_int64_t)
    half_create_from_files=HALF_SUCCESS
  end function

  integer(c_int) function half_create_from_files_backend(charge_c,potential_c,encut,xc,use_uspp,backend,solver, &
      handle,error,error_capacity)bind(C,name='half_create_from_files_backend')
    character(c_char),intent(in)::charge_c(*),potential_c(*)
    real(c_double),value::encut
    integer(c_int),value::xc,use_uspp,backend,solver,error_capacity
    integer(c_int64_t),intent(out)::handle
    character(c_char),intent(out)::error(*)
    character(len=:),allocatable::charge,potential
    character(len=512)::message
    character(len=3)::xc_name,solver_name
    integer::slot,status,i
    handle=0;call clear_error(error,error_capacity)
    call import_string(charge_c,charge);call import_string(potential_c,potential)
    if(len(charge)==0.or.len(potential)==0.or.(xc/=1.and.xc/=2).or.backend<0.or.backend>2.or.solver<1.or.solver>2)then
      half_create_from_files_backend=HALF_INVALID_ARGUMENT
      call export_error('invalid path, XC, backend, or solver selector',error,error_capacity);return
    end if
    slot=0;do i=1,MAX_CONTEXTS;if(.not.used(i))then;slot=i;exit;end if;end do
    if(slot==0)then
      half_create_from_files_backend=HALF_INTERNAL;call export_error('HALF context registry is full',error,error_capacity);return
    end if
    xc_name=merge('pbe','lda',xc==2);solver_name=merge('evj','evd',solver==2)
    call contexts(slot)%initialize_files(charge,potential,real(encut,dp),xc_name,use_uspp/=0,status,message,backend,solver_name)
    if(status/=LIB_SUCCESS)then
      half_create_from_files_backend=int(status,c_int);call export_error(trim(message),error,error_capacity);return
    end if
    used(slot)=.true.;generation(slot)=generation(slot)+1
    handle=generation(slot)*int(MAX_CONTEXTS,c_int64_t)+int(slot,c_int64_t)
    half_create_from_files_backend=HALF_SUCCESS
  end function

  integer(c_int) function half_destroy(handle,error,error_capacity)bind(C,name='half_destroy')
    integer(c_int64_t),value::handle
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    integer::slot
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_destroy=HALF_INVALID_HANDLE;call export_error('invalid HALF context handle',error,error_capacity);return;end if
    call contexts(slot)%clear();used(slot)=.false.;half_destroy=HALF_SUCCESS
  end function

  integer(c_int) function half_get_system_info(handle,nions,ntypes,grid,lattice,volume,error,error_capacity) &
      bind(C,name='half_get_system_info')
    integer(c_int64_t),value::handle
    integer(c_int),intent(out)::nions,ntypes,grid(3)
    real(c_double),intent(out)::lattice(9),volume
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    integer::slot,i,j,k
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_get_system_info=bad_handle(error,error_capacity);return;end if
    nions=contexts(slot)%crystal%nions;ntypes=contexts(slot)%crystal%ntypes;grid=contexts(slot)%charge%shape
    k=0;do i=1,3;do j=1,3;k=k+1;lattice(k)=contexts(slot)%crystal%lattice(i,j);end do;end do
    volume=contexts(slot)%crystal%volume;half_get_system_info=HALF_SUCCESS
  end function

  integer(c_int) function half_get_basis_size(handle,kpoint,npw,error,error_capacity)bind(C,name='half_get_basis_size')
    integer(c_int64_t),value::handle
    real(c_double),intent(in)::kpoint(3)
    integer(c_int64_t),intent(out)::npw
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    type(plane_wave_basis_t)::basis
    character(len=512)::message
    integer::slot,status
    call clear_error(error,error_capacity);npw=0;slot=context_slot(handle)
    if(slot==0)then;half_get_basis_size=bad_handle(error,error_capacity);return;end if
    call contexts(slot)%make_basis(real(kpoint,dp),basis,status,message)
    if(status/=LIB_SUCCESS)then;half_get_basis_size=status;call export_error(trim(message),error,error_capacity);return;end if
    npw=basis%npw;half_get_basis_size=HALF_SUCCESS
  end function

  integer(c_int) function half_get_basis(handle,kpoint,capacity,gvec,kinetic,error,error_capacity)bind(C,name='half_get_basis')
    integer(c_int64_t),value::handle,capacity
    real(c_double),intent(in)::kpoint(3)
    integer(c_int32_t),intent(out)::gvec(*)
    real(c_double),intent(out)::kinetic(*)
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    type(plane_wave_basis_t)::basis
    character(len=512)::message
    integer::slot,status,i,j
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_get_basis=bad_handle(error,error_capacity);return;end if
    call contexts(slot)%make_basis(real(kpoint,dp),basis,status,message)
    if(status/=LIB_SUCCESS)then;half_get_basis=status;call export_error(trim(message),error,error_capacity);return;end if
    if(capacity<basis%npw)then;half_get_basis=HALF_CAPACITY;call export_error('basis output capacity is too small',error,error_capacity);return;end if
    do i=1,basis%npw
      do j=1,3;gvec(3*(i-1)+j)=basis%gvectors(i,j);end do
      kinetic(i)=basis%kinetic(i)
    end do
    half_get_basis=HALF_SUCCESS
  end function

  integer(c_int) function half_assemble_hs(handle,kpoint,ld,h_c,s_c,error,error_capacity)bind(C,name='half_assemble_hs')
    integer(c_int64_t),value::handle,ld
    real(c_double),intent(in)::kpoint(3)
    type(c_ptr),value::h_c,s_c
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    complex(c_double_complex),pointer::h_out(:,:),s_out(:,:)
    complex(dp),allocatable::h(:,:),s(:,:)
    character(len=512)::message
    integer::slot,status,n
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_assemble_hs=bad_handle(error,error_capacity);return;end if
    if(.not.c_associated(h_c).or..not.c_associated(s_c))then
      half_assemble_hs=HALF_INVALID_ARGUMENT;call export_error('H and S output pointers are required',error,error_capacity);return
    end if
    call contexts(slot)%assemble_hs(real(kpoint,dp),h,s,status,message)
    if(status/=LIB_SUCCESS)then;half_assemble_hs=status;call export_error(trim(message),error,error_capacity);return;end if
    n=size(h,1)
    if(ld<n)then;half_assemble_hs=HALF_CAPACITY;call export_error('matrix leading dimension is too small',error,error_capacity);return;end if
    call c_f_pointer(h_c,h_out,[int(ld),n]);call c_f_pointer(s_c,s_out,[int(ld),n])
    h_out(1:n,:)=h;s_out(1:n,:)=s;half_assemble_hs=HALF_SUCCESS
  end function

  integer(c_int) function half_apply_hs(handle,kpoint,nstates,ld,psi_c,hpsi_c,spsi_c,error,error_capacity) &
      bind(C,name='half_apply_hs')
    integer(c_int64_t),value::handle,ld
    integer(c_int),value::nstates,error_capacity
    real(c_double),intent(in)::kpoint(3)
    type(c_ptr),value::psi_c,hpsi_c,spsi_c
    character(c_char),intent(out)::error(*)
    complex(c_double_complex),pointer::psi(:,:),hpsi(:,:),spsi(:,:)
    complex(dp),allocatable::ho(:,:),so(:,:)
    integer(c_int64_t)::npw
    character(len=512)::message
    integer::slot,status
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_apply_hs=bad_handle(error,error_capacity);return;end if
    if(nstates<1.or..not.c_associated(psi_c).or..not.c_associated(hpsi_c).or..not.c_associated(spsi_c))then
      half_apply_hs=HALF_INVALID_ARGUMENT;call export_error('invalid state block arguments',error,error_capacity);return
    end if
    status=half_get_basis_size(handle,kpoint,npw,error,error_capacity);if(status/=HALF_SUCCESS)then;half_apply_hs=status;return;end if
    if(ld<npw)then;half_apply_hs=HALF_CAPACITY;call export_error('state leading dimension is too small',error,error_capacity);return;end if
    call c_f_pointer(psi_c,psi,[int(ld),nstates]);call c_f_pointer(hpsi_c,hpsi,[int(ld),nstates]);call c_f_pointer(spsi_c,spsi,[int(ld),nstates])
    call contexts(slot)%apply_hs(real(kpoint,dp),psi(1:npw,:),ho,so,status,message)
    if(status/=LIB_SUCCESS)then;half_apply_hs=status;call export_error(trim(message),error,error_capacity);return;end if
    hpsi(1:npw,:)=ho;spsi(1:npw,:)=so;half_apply_hs=HALF_SUCCESS
  end function

  integer(c_int) function half_solve_kpoint(handle,kpoint,nbands,eigenvalues,eigenvectors_c,ld,overlap_min,overlap_max, &
      error,error_capacity)bind(C,name='half_solve_kpoint')
    integer(c_int64_t),value::handle,ld
    integer(c_int),value::nbands,error_capacity
    real(c_double),intent(in)::kpoint(3)
    real(c_double),intent(out)::eigenvalues(*),overlap_min,overlap_max
    type(c_ptr),value::eigenvectors_c
    character(c_char),intent(out)::error(*)
    complex(c_double_complex),pointer::vectors_out(:,:)
    real(dp),allocatable::values(:)
    complex(dp),allocatable::vectors(:,:)
    type(plane_wave_basis_t)::basis
    character(len=512)::message
    integer::slot,status
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_solve_kpoint=bad_handle(error,error_capacity);return;end if
    if(nbands<1)then;half_solve_kpoint=HALF_INVALID_ARGUMENT;call export_error('band count must be positive',error,error_capacity);return;end if
    if(c_associated(eigenvectors_c))then
      call contexts(slot)%solve_kpoint(real(kpoint,dp),nbands,values,vectors,overlap_min,overlap_max,status,message,basis)
    else
      call contexts(slot)%solve_kpoint(real(kpoint,dp),nbands,values,overlap_min=overlap_min,overlap_max=overlap_max, &
        status=status,message=message,basis_out=basis)
    end if
    if(status/=LIB_SUCCESS)then;half_solve_kpoint=status;call export_error(trim(message),error,error_capacity);return;end if
    eigenvalues(:nbands)=values
    if(c_associated(eigenvectors_c))then
      if(ld<basis%npw)then;half_solve_kpoint=HALF_CAPACITY;call export_error('eigenvector leading dimension is too small',error,error_capacity);return;end if
      call c_f_pointer(eigenvectors_c,vectors_out,[int(ld),nbands]);vectors_out(1:basis%npw,:)=vectors
    end if
    half_solve_kpoint=HALF_SUCCESS
  end function

  integer(c_int) function half_solve_kpoint_mapped(handle,kpoint,nbands,host_npw,host_gvec,eigenvalues, &
      eigenvectors_c,ld,overlap_min,overlap_max,error,error_capacity)bind(C,name='half_solve_kpoint_mapped')
    integer(c_int64_t),value::handle,host_npw,ld
    integer(c_int),value::nbands,error_capacity
    real(c_double),intent(in)::kpoint(3)
    integer(c_int32_t),intent(in)::host_gvec(*)
    real(c_double),intent(out)::eigenvalues(*),overlap_min,overlap_max
    type(c_ptr),value::eigenvectors_c
    character(c_char),intent(out)::error(*)
    complex(c_double_complex),pointer::vectors_out(:,:)
    real(dp),allocatable::values(:)
    complex(dp),allocatable::vectors(:,:)
    type(plane_wave_basis_t)::basis
    integer,allocatable::table(:)
    character(len=512)::message
    integer::slot,status,n,table_size,i,j,p,host_i
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_solve_kpoint_mapped=bad_handle(error,error_capacity);return;end if
    if(nbands<1.or.host_npw<1.or.ld<host_npw.or..not.c_associated(eigenvectors_c))then
      half_solve_kpoint_mapped=HALF_INVALID_ARGUMENT
      call export_error('invalid mapped-solve band count, basis size, leading dimension, or output pointer',error,error_capacity)
      return
    end if
    call contexts(slot)%solve_kpoint(real(kpoint,dp),nbands,values,vectors,overlap_min,overlap_max,status,message,basis)
    if(status/=LIB_SUCCESS)then
      half_solve_kpoint_mapped=status;call export_error(trim(message),error,error_capacity);return
    end if
    n=basis%npw
    if(host_npw/=n)then
      half_solve_kpoint_mapped=HALF_INVALID_ARGUMENT
      call export_error('host and HALF plane-wave bases have different sizes',error,error_capacity);return
    end if
    table_size=1
    do while(table_size<2*n);table_size=table_size*2;end do
    allocate(table(table_size));table=0
    do host_i=1,n
      p=g_hash(host_gvec(3*host_i-2),host_gvec(3*host_i-1),host_gvec(3*host_i),table_size)
      do while(table(p)/=0)
        i=table(p)
        if(host_gvec(3*i-2)==host_gvec(3*host_i-2).and.host_gvec(3*i-1)==host_gvec(3*host_i-1).and. &
            host_gvec(3*i)==host_gvec(3*host_i))then
          half_solve_kpoint_mapped=HALF_INVALID_ARGUMENT
          call export_error('host G-vector list contains a duplicate',error,error_capacity);return
        end if
        p=mod(p,table_size)+1
      end do
      table(p)=host_i
    end do
    call c_f_pointer(eigenvectors_c,vectors_out,[int(ld),nbands]);vectors_out=cmplx(0.0_c_double,0.0_c_double,c_double_complex)
    do i=1,n
      p=g_hash(int(basis%gvectors(i,1),c_int32_t),int(basis%gvectors(i,2),c_int32_t), &
        int(basis%gvectors(i,3),c_int32_t),table_size)
      host_i=0
      do while(table(p)/=0)
        j=table(p)
        if(host_gvec(3*j-2)==basis%gvectors(i,1).and.host_gvec(3*j-1)==basis%gvectors(i,2).and. &
            host_gvec(3*j)==basis%gvectors(i,3))then;host_i=j;exit;end if
        p=mod(p,table_size)+1
      end do
      if(host_i==0)then
        half_solve_kpoint_mapped=HALF_INVALID_ARGUMENT
        call export_error('host and HALF plane-wave bases contain different G vectors',error,error_capacity);return
      end if
      vectors_out(host_i,1:nbands)=vectors(i,1:nbands)
    end do
    eigenvalues(:nbands)=values
    half_solve_kpoint_mapped=HALF_SUCCESS
  end function

  integer function g_hash(gx,gy,gz,table_size)result(position)
    integer(c_int32_t),intent(in)::gx,gy,gz
    integer,intent(in)::table_size
    integer(c_int64_t)::value
    value=int(gx,c_int64_t)*73856093_c_int64_t+int(gy,c_int64_t)*19349663_c_int64_t+ &
      int(gz,c_int64_t)*83492791_c_int64_t
    position=int(modulo(value,int(table_size,c_int64_t)))+1
  end function

  integer function context_slot(handle)result(slot)
    integer(c_int64_t),intent(in)::handle
    integer(c_int64_t)::g
    slot=int(modulo(handle,int(MAX_CONTEXTS,c_int64_t)))
    if(slot==0)slot=MAX_CONTEXTS
    g=(handle-int(slot,c_int64_t))/int(MAX_CONTEXTS,c_int64_t)
    if(handle<=0.or.slot<1.or.slot>MAX_CONTEXTS.or..not.used(slot).or.g/=generation(slot))slot=0
  end function

  integer(c_int) function bad_handle(error,capacity)
    character(c_char),intent(out)::error(*)
    integer(c_int),intent(in)::capacity
    call export_error('invalid HALF context handle',error,capacity);bad_handle=HALF_INVALID_HANDLE
  end function

  subroutine import_string(input,output)
    character(c_char),intent(in)::input(*)
    character(len=:),allocatable,intent(out)::output
    integer::n,i
    n=0;do while(n<4096);if(input(n+1)==c_null_char)exit;n=n+1;end do
    allocate(character(len=n)::output);do i=1,n;output(i:i)=input(i);end do
  end subroutine

  subroutine clear_error(error,capacity)
    character(c_char),intent(out)::error(*)
    integer(c_int),intent(in)::capacity
    if(capacity>0)error(1)=c_null_char
  end subroutine

  subroutine export_error(message,error,capacity)
    character(len=*),intent(in)::message
    character(c_char),intent(out)::error(*)
    integer(c_int),intent(in)::capacity
    integer::i,n
    if(capacity<=0)return
    n=min(len_trim(message),int(capacity)-1)
    do i=1,n;error(i)=message(i:i);end do;error(n+1)=c_null_char
  end subroutine
end module half_c_api
