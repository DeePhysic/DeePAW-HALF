module half_c_api
  use iso_c_binding,only:c_int,c_int32_t,c_int64_t,c_float,c_double,c_double_complex,c_char,c_null_char,c_ptr,c_loc, &
    c_f_pointer,c_associated,c_sizeof,c_null_ptr
  use half_kinds,only:dp
  use half_types,only:plane_wave_basis_t,potcar_t
  use half_potcar,only:read_potcar
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
  character(c_char),target,save::version(6)=[character(c_char)::'0','.', '6','.', '0',c_null_char]
  type,bind(C)::half_request_geometry_v1_c
    integer(c_int32_t)::struct_size,flags,nions,ntypes,grid(3),reserved_i32
    real(c_double)::lattice(9)
    type(c_ptr)::species,positions_fractional
    integer(c_int64_t)::reserved(8)
  end type
  type,bind(C)::half_escn_request_v1_c
    integer(c_int32_t)::struct_size,flags,nions,grid(3),reserved_i32
    real(c_double)::cell(9)
    type(c_ptr)::atomic_numbers,positions_cartesian
    integer(c_int64_t)::reserved(8)
  end type
  type,bind(C)::half_escn_result_v1_c
    integer(c_int32_t)::struct_size,flags,grid(3),reserved_i32
    integer(c_int64_t)::capacity
    type(c_ptr)::density,nu,alpha,beta,risk
    real(c_double)::elapsed_seconds
    integer(c_int64_t)::reserved(8)
  end type
  interface
    integer(c_int) function half_escn_predict_c(base_url,request,result,timeout,error,error_capacity) &
        bind(C,name='half_escn_predict')
      import::c_char,c_int,c_ptr
      character(c_char),intent(in)::base_url(*)
      type(c_ptr),value::request,result
      integer(c_int),value::timeout,error_capacity
      character(c_char),intent(out)::error(*)
    end function
  end interface
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
#ifdef HALF_HAVE_ESCN_API
    half_get_capabilities=ior(half_get_capabilities,32)
#endif
  end function

  integer(c_int) function half_create_from_escn(endpoint_c,potential_c,geometry_c,encut,xc,use_uspp,backend,solver, &
      normalization,timeout,handle,error,error_capacity)bind(C,name='half_create_from_escn')
    character(c_char),intent(in)::endpoint_c(*),potential_c(*)
    type(c_ptr),value::geometry_c
    real(c_double),value::encut
    integer(c_int),value::xc,use_uspp,backend,solver,normalization,timeout,error_capacity
    integer(c_int64_t),intent(out)::handle
    character(c_char),intent(out)::error(*)
    type(half_request_geometry_v1_c),pointer::geometry
    type(half_request_geometry_v1_c)::geometry_layout
    type(half_escn_request_v1_c),target::api_request
    type(half_escn_result_v1_c),target::api_result
    integer(c_int32_t),pointer::species_c(:)
    real(c_double),pointer::positions_c(:)
    integer(c_int32_t),allocatable,target::atomic_numbers(:),species(:)
    real(c_double),allocatable,target::positions_cartesian(:)
    real(c_float),allocatable,target::api_density(:)
    real(dp),allocatable::positions(:,:),density_c(:),density_f(:)
    type(potcar_t),allocatable::potcars(:)
    real(dp)::lattice(3,3),target_electrons,model_electrons,scale
    character(len=:),allocatable::endpoint,potential
    character(len=512)::message
    character(len=3)::xc_name,solver_name
    integer(c_int)::remote_status
    integer(c_int64_t)::npoints64
    integer::slot,status,i,j,i1,i2,i3,source,destination,npoints,z
    handle=0_c_int64_t;call clear_error(error,error_capacity)
    call import_string(endpoint_c,endpoint);call import_string(potential_c,potential)
    if(len(endpoint)==0.or.len(potential)==0.or..not.c_associated(geometry_c).or. &
        (xc/=1.and.xc/=2).or.backend<0.or.backend>2.or.solver<1.or.solver>4.or. &
        (normalization/=0.and.normalization/=1))then
      half_create_from_escn=HALF_INVALID_ARGUMENT
      call export_error('invalid eSCN URL, path, geometry, execution policy, or normalization',error,error_capacity);return
    end if
    call c_f_pointer(geometry_c,geometry)
    if(geometry%struct_size<int(c_sizeof(geometry_layout),c_int32_t).or.geometry%flags/=0.or.geometry%nions<1.or. &
        geometry%nions>10000.or. &
        geometry%ntypes<1.or.any(geometry%grid<1).or..not.c_associated(geometry%species).or. &
        .not.c_associated(geometry%positions_fractional))then
      half_create_from_escn=HALF_INVALID_ARGUMENT
      call export_error('invalid request geometry structure',error,error_capacity);return
    end if
    npoints64=int(geometry%grid(1),c_int64_t)*int(geometry%grid(2),c_int64_t)*int(geometry%grid(3),c_int64_t)
    if(npoints64<1_c_int64_t.or.npoints64>2000000_c_int64_t)then
      half_create_from_escn=HALF_INVALID_ARGUMENT
      call export_error('eSCN grid must contain at most 2000000 points',error,error_capacity);return
    end if
    npoints=int(npoints64)
    slot=0
    do i=1,MAX_CONTEXTS;if(.not.used(i))then;slot=i;exit;end if;end do
    if(slot==0)then
      half_create_from_escn=HALF_INTERNAL;call export_error('HALF context registry is full',error,error_capacity);return
    end if
    call c_f_pointer(geometry%species,species_c,[int(geometry%nions)])
    call c_f_pointer(geometry%positions_fractional,positions_c,[3*int(geometry%nions)])
    allocate(species(geometry%nions),positions(geometry%nions,3),atomic_numbers(geometry%nions), &
      positions_cartesian(3*geometry%nions),api_density(npoints),density_c(npoints),density_f(npoints))
    species=species_c
    do i=1,geometry%nions
      do j=1,3;positions(i,j)=positions_c(3*(i-1)+j);end do
    end do
    do i=1,3
      do j=1,3;lattice(i,j)=geometry%lattice(3*(i-1)+j);end do
    end do
    if(any(species<1).or.any(species>geometry%ntypes))then
      half_create_from_escn=HALF_INVALID_ARGUMENT
      call export_error('request species indices are outside the POTCAR type range',error,error_capacity);return
    end if
    call read_potcar(potential,potcars)
    if(size(potcars)/=geometry%ntypes)then
      half_create_from_escn=HALF_INVALID_ARGUMENT
      call export_error('POTCAR dataset count does not match request type count',error,error_capacity);return
    end if
    do i=1,geometry%nions
      z=element_atomic_number(potcars(species(i))%element)
      if(z==0)then
        half_create_from_escn=HALF_INVALID_ARGUMENT
        call export_error('cannot map a POTCAR element to an atomic number',error,error_capacity);return
      end if
      atomic_numbers(i)=z
      do j=1,3
        positions_cartesian(3*(i-1)+j)=dot_product(positions(i,:),lattice(:,j))
      end do
    end do
    api_request%struct_size=int(c_sizeof(api_request),c_int32_t);api_request%flags=0_c_int32_t
    api_request%nions=geometry%nions;api_request%grid=geometry%grid;api_request%reserved_i32=0_c_int32_t
    do i=1,3
      do j=1,3;api_request%cell(3*(i-1)+j)=lattice(i,j);end do
    end do
    api_request%atomic_numbers=c_loc(atomic_numbers);api_request%positions_cartesian=c_loc(positions_cartesian)
    api_request%reserved=0_c_int64_t
    api_result%struct_size=int(c_sizeof(api_result),c_int32_t);api_result%flags=0_c_int32_t
    api_result%grid=0_c_int32_t;api_result%reserved_i32=0_c_int32_t;api_result%capacity=int(npoints,c_int64_t)
    api_result%density=c_loc(api_density);api_result%nu=c_null_ptr;api_result%alpha=c_null_ptr
    api_result%beta=c_null_ptr;api_result%risk=c_null_ptr;api_result%elapsed_seconds=0.0_c_double
    api_result%reserved=0_c_int64_t
    remote_status=half_escn_predict_c(endpoint_c,c_loc(api_request),c_loc(api_result),timeout,error,error_capacity)
    if(remote_status/=HALF_SUCCESS)then;half_create_from_escn=remote_status;return;end if
    density_c=real(api_density,dp)
    if(normalization==1)then
      target_electrons=0.0_dp
      do i=1,geometry%nions;target_electrons=target_electrons+potcars(species(i))%zval;end do
      model_electrons=sum(density_c)/real(npoints,dp)
      if(abs(model_electrons)<=tiny(1.0_dp))then
        half_create_from_escn=HALF_INVALID_ARGUMENT
        call export_error('cannot valence-normalize an eSCN density with zero integral',error,error_capacity);return
      end if
      scale=target_electrons/model_electrons;density_c=density_c*scale
    end if
    do i1=0,geometry%grid(1)-1
      do i2=0,geometry%grid(2)-1
        do i3=0,geometry%grid(3)-1
          source=i1*geometry%grid(2)*geometry%grid(3)+i2*geometry%grid(3)+i3+1
          destination=i1+geometry%grid(1)*(i2+geometry%grid(2)*i3)+1
          density_f(destination)=density_c(source)
        end do
      end do
    end do
    xc_name=merge('pbe','lda',xc==2)
    select case(solver)
    case(1);solver_name='evd'
    case(2);solver_name='evj'
    case(3);solver_name='evx'
    case(4);solver_name='acc'
    end select
    call contexts(slot)%initialize_density(potential,real(encut,dp),xc_name,use_uspp/=0,geometry%nions,geometry%ntypes, &
      geometry%grid,lattice,species,positions,density_f,status,message,backend,solver_name)
    if(status/=LIB_SUCCESS)then
      half_create_from_escn=int(status,c_int);call export_error(trim(message),error,error_capacity);return
    end if
    used(slot)=.true.;generation(slot)=generation(slot)+1
    handle=generation(slot)*int(MAX_CONTEXTS,c_int64_t)+int(slot,c_int64_t)
    half_create_from_escn=HALF_SUCCESS
  end function half_create_from_escn

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
    if(len(charge)==0.or.len(potential)==0.or.(xc/=1.and.xc/=2).or.backend<0.or.backend>2.or.solver<1.or.solver>4)then
      half_create_from_files_backend=HALF_INVALID_ARGUMENT
      call export_error('invalid path, XC, backend, or solver selector',error,error_capacity);return
    end if
    slot=0;do i=1,MAX_CONTEXTS;if(.not.used(i))then;slot=i;exit;end if;end do
    if(slot==0)then
      half_create_from_files_backend=HALF_INTERNAL;call export_error('HALF context registry is full',error,error_capacity);return
    end if
    xc_name=merge('pbe','lda',xc==2)
    select case(solver)
    case(1);solver_name='evd'
    case(2);solver_name='evj'
    case(3);solver_name='evx'
    case(4);solver_name='acc'
    end select
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

  integer(c_int) function half_set_request_geometry(handle,geometry_c,error,error_capacity) &
      bind(C,name='half_set_request_geometry')
    integer(c_int64_t),value::handle
    type(c_ptr),value::geometry_c
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    type(half_request_geometry_v1_c),pointer::geometry
    type(half_request_geometry_v1_c)::geometry_layout
    integer(c_int32_t),pointer::species_c(:)
    real(c_double),pointer::positions_c(:)
    integer(c_int32_t),allocatable::species(:)
    real(dp),allocatable::positions(:,:)
    real(dp)::lattice(3,3)
    character(len=512)::message
    integer::slot,status,i,j
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_set_request_geometry=bad_handle(error,error_capacity);return;end if
    if(.not.c_associated(geometry_c))then
      half_set_request_geometry=HALF_INVALID_ARGUMENT
      call export_error('request geometry pointer is required',error,error_capacity);return
    end if
    call c_f_pointer(geometry_c,geometry)
    if(geometry%struct_size<int(c_sizeof(geometry_layout),c_int32_t).or.geometry%flags/=0)then
      half_set_request_geometry=HALF_INVALID_ARGUMENT
      call export_error('unsupported request geometry structure size or flags',error,error_capacity);return
    end if
    if(geometry%nions<1.or.geometry%ntypes<1.or..not.c_associated(geometry%species).or. &
        .not.c_associated(geometry%positions_fractional))then
      half_set_request_geometry=HALF_INVALID_ARGUMENT
      call export_error('request geometry arrays and positive dimensions are required',error,error_capacity);return
    end if
    call c_f_pointer(geometry%species,species_c,[int(geometry%nions)])
    call c_f_pointer(geometry%positions_fractional,positions_c,[3*int(geometry%nions)])
    allocate(species(geometry%nions),positions(geometry%nions,3))
    species=species_c
    do i=1,geometry%nions
      do j=1,3;positions(i,j)=positions_c(3*(i-1)+j);end do
    end do
    do i=1,3
      do j=1,3;lattice(i,j)=geometry%lattice(3*(i-1)+j);end do
    end do
    call contexts(slot)%set_request_geometry(geometry%nions,geometry%ntypes,geometry%grid,lattice,species,positions,status,message)
    if(status/=LIB_SUCCESS)then
      half_set_request_geometry=int(status,c_int);call export_error(trim(message),error,error_capacity);return
    end if
    half_set_request_geometry=HALF_SUCCESS
  end function

  integer(c_int) function half_get_request_geometry(handle,nions_c,ntypes_c,grid_c,lattice_c,species_capacity, &
      species_c,positions_capacity,positions_c,error,error_capacity)bind(C,name='half_get_request_geometry')
    integer(c_int64_t),value::handle,species_capacity,positions_capacity
    type(c_ptr),value::nions_c,ntypes_c,grid_c,lattice_c,species_c,positions_c
    integer(c_int),value::error_capacity
    character(c_char),intent(out)::error(*)
    integer(c_int32_t),pointer::nions,ntypes,grid(:),species(:)
    real(c_double),pointer::lattice(:),positions(:)
    integer::slot,i,j
    call clear_error(error,error_capacity);slot=context_slot(handle)
    if(slot==0)then;half_get_request_geometry=bad_handle(error,error_capacity);return;end if
    if(.not.contexts(slot)%has_request_geometry)then
      half_get_request_geometry=HALF_UNAVAILABLE
      call export_error('request geometry has not been set',error,error_capacity);return
    end if
    if(.not.c_associated(nions_c).or..not.c_associated(ntypes_c))then
      half_get_request_geometry=HALF_INVALID_ARGUMENT
      call export_error('nions and ntypes output pointers are required',error,error_capacity);return
    end if
    call c_f_pointer(nions_c,nions);call c_f_pointer(ntypes_c,ntypes)
    nions=contexts(slot)%request_nions;ntypes=contexts(slot)%request_ntypes
    if(c_associated(grid_c))then
      call c_f_pointer(grid_c,grid,[3]);grid=contexts(slot)%request_grid
    end if
    if(c_associated(lattice_c))then
      call c_f_pointer(lattice_c,lattice,[9])
      do i=1,3
        do j=1,3;lattice(3*(i-1)+j)=contexts(slot)%request_lattice(i,j);end do
      end do
    end if
    if(c_associated(species_c))then
      if(species_capacity<contexts(slot)%request_nions)then
        half_get_request_geometry=HALF_CAPACITY
        call export_error('request species output capacity is too small',error,error_capacity);return
      end if
      call c_f_pointer(species_c,species,[int(species_capacity)])
      species(1:contexts(slot)%request_nions)=contexts(slot)%request_species
    end if
    if(c_associated(positions_c))then
      if(positions_capacity<3_c_int64_t*contexts(slot)%request_nions)then
        half_get_request_geometry=HALF_CAPACITY
        call export_error('request positions output capacity is too small',error,error_capacity);return
      end if
      call c_f_pointer(positions_c,positions,[int(positions_capacity)])
      do i=1,contexts(slot)%request_nions
        do j=1,3;positions(3*(i-1)+j)=contexts(slot)%request_positions(i,j);end do
      end do
    end if
    half_get_request_geometry=HALF_SUCCESS
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

  integer function element_atomic_number(element)result(number)
    character(len=*),intent(in)::element
    character(len=3),parameter::symbols(118)=[character(len=3):: &
      'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', &
      'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr', &
      'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg', &
      'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
      'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']
    character(len=16)::wanted
    integer::i
    wanted=lower_ascii(adjustl(trim(element)));number=0
    do i=1,size(symbols)
      if(trim(wanted)==trim(lower_ascii(symbols(i))))then;number=i;return;end if
    end do
  end function element_atomic_number

  pure function lower_ascii(input)result(output)
    character(len=*),intent(in)::input
    character(len=len(input))::output
    integer::i,code
    output=input
    do i=1,len(input)
      code=iachar(input(i:i))
      if(code>=iachar('A').and.code<=iachar('Z'))output(i:i)=achar(code+32)
    end do
  end function lower_ascii

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
