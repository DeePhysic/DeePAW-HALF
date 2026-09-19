module half_api
  use iso_c_binding
  implicit none
  integer(c_int),parameter::HALF_SUCCESS=0,HALF_ERROR_INVALID_ARGUMENT=1,HALF_ERROR_IO=2, &
    HALF_ERROR_INVALID_HANDLE=3,HALF_ERROR_CAPACITY=4,HALF_ERROR_UNAVAILABLE=5,HALF_ERROR_INTERNAL=6
  integer(c_int),parameter::HALF_XC_LDA=1,HALF_XC_PBE=2
  integer(c_int),parameter::HALF_BACKEND_AUTO=0,HALF_BACKEND_CPU=1,HALF_BACKEND_CUDA=2
  integer(c_int),parameter::HALF_SOLVER_EVD=1,HALF_SOLVER_EVJ=2
  integer(c_int),parameter::HALF_CAP_CPU=1,HALF_CAP_CUDA=2,HALF_CAP_MPI=4,HALF_CAP_HDF5=8,HALF_CAP_SPGLIB=16
  interface
    integer(c_int) function half_get_abi_version()bind(C)
      import::c_int
    end function
    type(c_ptr) function half_get_version_string()bind(C)
      import::c_ptr
    end function
    integer(c_int) function half_create_from_files_backend(charge,potential,encut,xc,use_uspp,backend,solver,handle,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t
      character(c_char),intent(in)::charge(*),potential(*)
      real(c_double),value::encut
      integer(c_int),value::xc,use_uspp,backend,solver,error_capacity
      integer(c_int64_t),intent(out)::handle
      character(c_char),intent(out)::error(*)
    end function
    integer(c_int) function half_get_capabilities()bind(C)
      import::c_int
    end function
    integer(c_int) function half_create_from_files(charge,potential,encut,xc,use_uspp,handle,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t
      character(c_char),intent(in)::charge(*),potential(*)
      real(c_double),value::encut
      integer(c_int),value::xc,use_uspp,error_capacity
      integer(c_int64_t),intent(out)::handle
      character(c_char),intent(out)::error(*)
    end function
    integer(c_int) function half_destroy(handle,error,error_capacity)bind(C)
      import::c_char,c_int,c_int64_t
      integer(c_int64_t),value::handle
      character(c_char),intent(out)::error(*)
      integer(c_int),value::error_capacity
    end function
    integer(c_int) function half_get_system_info(handle,nions,ntypes,grid,lattice,volume,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t
      integer(c_int64_t),value::handle
      integer(c_int),intent(out)::nions,ntypes,grid(3)
      real(c_double),intent(out)::lattice(9),volume
      character(c_char),intent(out)::error(*)
      integer(c_int),value::error_capacity
    end function
    integer(c_int) function half_get_basis_size(handle,kpoint,npw,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t
      integer(c_int64_t),value::handle
      real(c_double),intent(in)::kpoint(3)
      integer(c_int64_t),intent(out)::npw
      character(c_char),intent(out)::error(*)
      integer(c_int),value::error_capacity
    end function
    integer(c_int) function half_get_basis(handle,kpoint,capacity,gvec,kinetic,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int32_t,c_int64_t
      integer(c_int64_t),value::handle,capacity
      real(c_double),intent(in)::kpoint(3)
      integer(c_int32_t),intent(out)::gvec(*)
      real(c_double),intent(out)::kinetic(*)
      character(c_char),intent(out)::error(*)
      integer(c_int),value::error_capacity
    end function
    integer(c_int) function half_assemble_hs(handle,kpoint,ld,h,s,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t,c_ptr
      integer(c_int64_t),value::handle,ld
      real(c_double),intent(in)::kpoint(3)
      type(c_ptr),value::h,s
      character(c_char),intent(out)::error(*)
      integer(c_int),value::error_capacity
    end function
    integer(c_int) function half_apply_hs(handle,kpoint,nstates,ld,psi,hpsi,spsi,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t,c_ptr
      integer(c_int64_t),value::handle,ld
      integer(c_int),value::nstates,error_capacity
      real(c_double),intent(in)::kpoint(3)
      type(c_ptr),value::psi,hpsi,spsi
      character(c_char),intent(out)::error(*)
    end function
    integer(c_int) function half_solve_kpoint(handle,kpoint,nbands,eigenvalues,eigenvectors,ld,smin,smax,error,error_capacity)bind(C)
      import::c_char,c_double,c_int,c_int64_t,c_ptr
      integer(c_int64_t),value::handle,ld
      integer(c_int),value::nbands,error_capacity
      real(c_double),intent(in)::kpoint(3)
      real(c_double),intent(out)::eigenvalues(*),smin,smax
      type(c_ptr),value::eigenvectors
      character(c_char),intent(out)::error(*)
    end function
  end interface
end module half_api
