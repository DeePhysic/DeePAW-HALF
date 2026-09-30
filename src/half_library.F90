module half_library
  use, intrinsic::ieee_arithmetic,only:ieee_is_finite
  use half_kinds,only:dp,i32
  use half_types,only:crystal_t,charge_grid_t,plane_wave_basis_t,potcar_t
  use half_chgcar,only:read_chgcar
  use half_potcar,only:read_potcar,validate_potcar_structure
  use half_basis,only:build_plane_wave_basis
  use half_paw,only:paw_species_t,build_paw_operators
#ifdef HALF_HAVE_MKL
  use half_fft,only:fft3_forward,fft3_backward
  use half_uspp,only:build_uspp_dij_cpu
  use half_potential,only:build_veff_lda,build_veff_pbe
  use half_dense_solver,only:assemble_dense_gamma,solve_dense_gamma
#endif
#ifdef HALF_HAVE_CUDA
  use cudafor
  use half_cuda_solver,only:solve_dense_gamma_cuda_full
  use half_cuda_potential,only:build_veff_cuda
  use half_cuda_assembly,only:cuda_hs_operator_t,prepare_hs_operator_cuda,apply_hs_operator_cuda,destroy_hs_operator_cuda
#endif
  implicit none
  private
  public::half_context_t,HALF_SUCCESS,HALF_ERROR_INVALID_ARGUMENT,HALF_ERROR_IO,HALF_ERROR_UNAVAILABLE, &
    HALF_BACKEND_AUTO,HALF_BACKEND_CPU,HALF_BACKEND_CUDA

  integer,parameter::HALF_SUCCESS=0,HALF_ERROR_INVALID_ARGUMENT=1,HALF_ERROR_IO=2,HALF_ERROR_UNAVAILABLE=5
  integer,parameter::HALF_BACKEND_AUTO=0,HALF_BACKEND_CPU=1,HALF_BACKEND_CUDA=2

  type::half_context_t
    type(crystal_t)::crystal
    type(charge_grid_t)::charge
    type(potcar_t),allocatable::potcars(:)
    real(dp),allocatable::veff(:)
    integer(i32),allocatable::request_species(:)
    real(dp),allocatable::request_positions(:,:)
    integer(i32)::request_grid(3)=0_i32,request_nions=0_i32,request_ntypes=0_i32
    real(dp)::request_lattice(3,3)=0.0_dp
    logical::has_request_geometry=.false.
    type(plane_wave_basis_t)::apply_basis
    type(paw_species_t),allocatable::apply_paw(:)
    real(dp)::apply_kpoint(3)=0.0_dp
    logical::apply_ready=.false.
#ifdef HALF_HAVE_CUDA
    type(cuda_hs_operator_t)::gpu_apply_operator
#endif
    real(dp)::encut=0.0_dp,e_hartree=0.0_dp,e_xc=0.0_dp,e_xc_potential=0.0_dp
    integer::backend=HALF_BACKEND_AUTO
    character(len=3)::solver='evd'
    logical::use_pbe=.true.,use_uspp=.true.,ready=.false.
  contains
    procedure::initialize_files=>context_initialize_files
    procedure::initialize_density=>context_initialize_density
    procedure::clear=>context_clear
    procedure::make_basis=>context_make_basis
    procedure::set_request_geometry=>context_set_request_geometry
    procedure::assemble_hs=>context_assemble_hs
    procedure::apply_hs=>context_apply_hs
    procedure::solve_kpoint=>context_solve_kpoint
  end type
contains
  subroutine context_initialize_files(self,charge_path,potential_path,encut,xc,use_uspp,status,message,backend,solver)
    class(half_context_t),intent(inout)::self
    character(len=*),intent(in)::charge_path,potential_path,xc
    character(len=*),intent(in),optional::solver
    real(dp),intent(in)::encut
    logical,intent(in)::use_uspp
    integer,intent(out)::status
    integer,intent(in),optional::backend
    character(len=*),intent(out)::message
    logical::exists
    call self%clear();status=HALF_SUCCESS;message=''
    if(encut<=0.0_dp.or.encut/=encut)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='ENCUT must be positive and finite';return
    end if
    if(trim(xc)/='lda'.and.trim(xc)/='pbe')then
      status=HALF_ERROR_INVALID_ARGUMENT;message='XC must be lda or pbe';return
    end if
    self%backend=HALF_BACKEND_AUTO;if(present(backend))self%backend=backend
    if(self%backend==HALF_BACKEND_AUTO)then
#ifdef HALF_HAVE_CUDA
      self%backend=HALF_BACKEND_CUDA
#else
      self%backend=HALF_BACKEND_CPU
#endif
    end if
    if(self%backend/=HALF_BACKEND_CPU.and.self%backend/=HALF_BACKEND_CUDA)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='backend must be AUTO, CPU, or CUDA';return
    end if
#ifndef HALF_HAVE_MKL
    if(self%backend==HALF_BACKEND_CPU)then;status=HALF_ERROR_UNAVAILABLE;message='CPU API backend requires oneMKL';return;end if
#endif
#ifndef HALF_HAVE_CUDA
    if(self%backend==HALF_BACKEND_CUDA)then;status=HALF_ERROR_UNAVAILABLE;message='CUDA API backend is not compiled';return;end if
#endif
    self%solver='evd';if(present(solver))self%solver=solver
    if(self%solver/='evd'.and.self%solver/='evj'.and.self%solver/='evx'.and.self%solver/='acc')then
      status=HALF_ERROR_INVALID_ARGUMENT;message='solver must be evd, evj, evx, or acc';return
    end if
    if(self%backend==HALF_BACKEND_CPU.and.self%solver/='evd')then
      status=HALF_ERROR_UNAVAILABLE;message='EVJ, EVX, and ACC are available only with CUDA';return
    end if
    inquire(file=trim(charge_path),exist=exists)
    if(.not.exists)then;status=HALF_ERROR_IO;message='charge input does not exist';return;end if
    inquire(file=trim(potential_path),exist=exists)
    if(.not.exists)then;status=HALF_ERROR_IO;message='potential input does not exist';return;end if
    call read_chgcar(trim(charge_path),self%crystal,self%charge)
    call read_potcar(trim(potential_path),self%potcars)
    call validate_potcar_structure(self%potcars,self%crystal)
    self%encut=encut;self%use_pbe=trim(xc)=='pbe';self%use_uspp=use_uspp
    if(self%backend==HALF_BACKEND_CPU)then
#ifdef HALF_HAVE_MKL
    if(self%use_pbe)then
      call build_veff_pbe(self%charge,self%potcars,self%crystal,self%veff,self%e_hartree,self%e_xc,self%e_xc_potential)
    else
      call build_veff_lda(self%charge,self%potcars,self%crystal,self%veff,self%e_hartree,self%e_xc,self%e_xc_potential)
    end if
#endif
    end if
    self%ready=.true.
  end subroutine

  subroutine context_initialize_density(self,potential_path,encut,xc,use_uspp,nions,ntypes,grid,lattice, &
      species,positions,density,status,message,backend,solver)
    class(half_context_t),intent(inout)::self
    character(len=*),intent(in)::potential_path,xc
    character(len=*),intent(in),optional::solver
    real(dp),intent(in)::encut,lattice(3,3),positions(:,:),density(:)
    integer(i32),intent(in)::nions,ntypes,grid(3),species(:)
    logical,intent(in)::use_uspp
    integer,intent(out)::status
    integer,intent(in),optional::backend
    character(len=*),intent(out)::message
    logical::exists
    integer::it,ion,next
    call self%clear();status=HALF_SUCCESS;message=''
    if(encut<=0.0_dp.or..not.ieee_is_finite(encut))then
      status=HALF_ERROR_INVALID_ARGUMENT;message='ENCUT must be positive and finite';return
    end if
    if(trim(xc)/='lda'.and.trim(xc)/='pbe')then
      status=HALF_ERROR_INVALID_ARGUMENT;message='XC must be lda or pbe';return
    end if
    if(nions<1.or.ntypes<1.or.any(grid<1).or.size(species)/=nions.or.size(positions,1)/=nions.or. &
        size(positions,2)/=3.or.size(density)/=product(grid))then
      status=HALF_ERROR_INVALID_ARGUMENT;message='invalid in-memory density dimensions';return
    end if
    if(any(species<1).or.any(species>ntypes).or.any(.not.ieee_is_finite(lattice)).or. &
        any(.not.ieee_is_finite(positions)).or.any(.not.ieee_is_finite(density)))then
      status=HALF_ERROR_INVALID_ARGUMENT;message='invalid in-memory geometry or density values';return
    end if
    if(abs(determinant3(lattice))<1.0e-12_dp)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='density lattice is singular';return
    end if
    inquire(file=trim(potential_path),exist=exists)
    if(.not.exists)then;status=HALF_ERROR_IO;message='potential input does not exist';return;end if

    self%backend=HALF_BACKEND_AUTO;if(present(backend))self%backend=backend
    if(self%backend==HALF_BACKEND_AUTO)then
#ifdef HALF_HAVE_CUDA
      self%backend=HALF_BACKEND_CUDA
#else
      self%backend=HALF_BACKEND_CPU
#endif
    end if
    if(self%backend/=HALF_BACKEND_CPU.and.self%backend/=HALF_BACKEND_CUDA)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='backend must be AUTO, CPU, or CUDA';return
    end if
#ifndef HALF_HAVE_MKL
    if(self%backend==HALF_BACKEND_CPU)then;status=HALF_ERROR_UNAVAILABLE;message='CPU API backend requires oneMKL';return;end if
#endif
#ifndef HALF_HAVE_CUDA
    if(self%backend==HALF_BACKEND_CUDA)then;status=HALF_ERROR_UNAVAILABLE;message='CUDA API backend is not compiled';return;end if
#endif
    self%solver='evd';if(present(solver))self%solver=solver
    if(self%solver/='evd'.and.self%solver/='evj'.and.self%solver/='evx'.and.self%solver/='acc')then
      status=HALF_ERROR_INVALID_ARGUMENT;message='solver must be evd, evj, evx, or acc';return
    end if
    if(self%backend==HALF_BACKEND_CPU.and.self%solver/='evd')then
      status=HALF_ERROR_UNAVAILABLE;message='EVJ, EVX, and ACC are available only with CUDA';return
    end if

    call read_potcar(trim(potential_path),self%potcars)
    if(size(self%potcars)/=ntypes)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='POTCAR dataset count does not match request type count';call self%clear();return
    end if
    self%crystal%system_name='DeePAW-eSCN remote density'
    self%crystal%nions=nions;self%crystal%ntypes=ntypes;self%crystal%lattice=lattice
    allocate(self%crystal%species(ntypes),self%crystal%counts(ntypes),self%crystal%positions(nions,3))
    self%crystal%counts=0
    do it=1,ntypes
      self%crystal%species(it)=self%potcars(it)%element
      self%crystal%counts(it)=count(species==it)
      if(self%crystal%counts(it)==0)then
        status=HALF_ERROR_INVALID_ARGUMENT;message='each request type must be used by at least one atom';call self%clear();return
      end if
    end do
    next=0
    do it=1,ntypes
      do ion=1,nions
        if(species(ion)==it)then;next=next+1;self%crystal%positions(next,:)=positions(ion,:);end if
      end do
    end do
    call self%crystal%update_geometry()
    call validate_potcar_structure(self%potcars,self%crystal)
    self%charge%shape=grid;self%charge%values=density
    self%encut=encut;self%use_pbe=trim(xc)=='pbe';self%use_uspp=use_uspp
    self%request_nions=nions;self%request_ntypes=ntypes;self%request_grid=grid
    self%request_lattice=lattice;self%request_species=species;self%request_positions=positions
    self%has_request_geometry=.true.
    if(self%backend==HALF_BACKEND_CPU)then
#ifdef HALF_HAVE_MKL
      if(self%use_pbe)then
        call build_veff_pbe(self%charge,self%potcars,self%crystal,self%veff,self%e_hartree,self%e_xc,self%e_xc_potential)
      else
        call build_veff_lda(self%charge,self%potcars,self%crystal,self%veff,self%e_hartree,self%e_xc,self%e_xc_potential)
      end if
#endif
    end if
    self%ready=.true.
  end subroutine context_initialize_density

  subroutine context_clear(self)
    class(half_context_t),intent(inout)::self
#ifdef HALF_HAVE_CUDA
    if(self%gpu_apply_operator%n>0)call destroy_hs_operator_cuda(self%gpu_apply_operator)
#endif
    if(allocated(self%apply_paw))deallocate(self%apply_paw)
    self%apply_ready=.false.;self%apply_kpoint=0.0_dp
    if(allocated(self%potcars))deallocate(self%potcars)
    if(allocated(self%veff))deallocate(self%veff)
    if(allocated(self%request_species))deallocate(self%request_species)
    if(allocated(self%request_positions))deallocate(self%request_positions)
    if(allocated(self%charge%values))deallocate(self%charge%values)
    if(allocated(self%crystal%system_name))deallocate(self%crystal%system_name)
    if(allocated(self%crystal%species))deallocate(self%crystal%species)
    if(allocated(self%crystal%counts))deallocate(self%crystal%counts)
    if(allocated(self%crystal%positions))deallocate(self%crystal%positions)
    self%ready=.false.;self%encut=0.0_dp;self%backend=HALF_BACKEND_AUTO;self%solver='evd'
    self%has_request_geometry=.false.;self%request_grid=0;self%request_nions=0;self%request_ntypes=0
    self%request_lattice=0.0_dp
  end subroutine

  subroutine context_set_request_geometry(self,nions,ntypes,grid,lattice,species,positions,status,message)
    class(half_context_t),intent(inout)::self
    integer(i32),intent(in)::nions,ntypes,grid(3),species(:)
    real(dp),intent(in)::lattice(3,3),positions(:,:)
    integer,intent(out)::status
    character(len=*),intent(out)::message
    status=HALF_SUCCESS;message=''
    if(.not.self%ready)then;status=HALF_ERROR_INVALID_ARGUMENT;message='HALF context is not initialized';return;end if
    if(nions<1.or.ntypes<1.or.any(grid<1).or.size(species)/=nions.or.size(positions,1)/=nions.or.size(positions,2)/=3)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='invalid request geometry dimensions';return
    end if
    if(any(species<1).or.any(species>ntypes).or.any(lattice/=lattice).or.any(positions/=positions))then
      status=HALF_ERROR_INVALID_ARGUMENT;message='invalid request geometry values';return
    end if
    if(abs(determinant3(lattice))<1.0e-12_dp)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='request lattice is singular';return
    end if
    self%request_nions=nions;self%request_ntypes=ntypes;self%request_grid=grid
    self%request_lattice=lattice;self%request_species=species;self%request_positions=positions
    self%has_request_geometry=.true.
  end subroutine

  pure real(dp) function determinant3(a)result(value)
    real(dp),intent(in)::a(3,3)
    value=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+ &
      a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function

  subroutine context_make_basis(self,kpoint,basis,status,message)
    class(half_context_t),intent(in)::self
    real(dp),intent(in)::kpoint(3)
    type(plane_wave_basis_t),intent(out)::basis
    integer,intent(out)::status
    character(len=*),intent(out)::message
    status=HALF_SUCCESS;message=''
    if(.not.self%ready)then;status=HALF_ERROR_INVALID_ARGUMENT;message='HALF context is not initialized';return;end if
    if(any(kpoint/=kpoint))then;status=HALF_ERROR_INVALID_ARGUMENT;message='k point must be finite';return;end if
    call build_plane_wave_basis(self%crystal,self%charge%shape,self%encut,kpoint,basis)
  end subroutine

  subroutine prepare(self,kpoint,basis,paw,status,message)
    class(half_context_t),intent(in)::self
    real(dp),intent(in)::kpoint(3)
    type(plane_wave_basis_t),intent(out)::basis
    type(paw_species_t),allocatable,intent(out)::paw(:)
    integer,intent(out)::status
    character(len=*),intent(out)::message
    call self%make_basis(kpoint,basis,status,message);if(status/=HALF_SUCCESS)return
    call build_paw_operators(self%potcars,self%crystal,basis,paw)
#ifdef HALF_HAVE_MKL
    if(self%use_uspp)call build_uspp_dij_cpu(self%veff,self%charge%shape,self%potcars,self%crystal,paw)
#endif
  end subroutine

  subroutine context_assemble_hs(self,kpoint,h,s,status,message,basis_out)
    class(half_context_t),intent(in)::self
    real(dp),intent(in)::kpoint(3)
    complex(dp),allocatable,intent(out)::h(:,:),s(:,:)
    integer,intent(out)::status
    character(len=*),intent(out)::message
    type(plane_wave_basis_t),intent(out),optional::basis_out
    type(plane_wave_basis_t)::basis
    type(paw_species_t),allocatable::paw(:)
    if(self%backend==HALF_BACKEND_CUDA)then
#ifdef HALF_HAVE_CUDA
      status=HALF_ERROR_UNAVAILABLE;message='CUDA host H/S export is not available in ABI v1; use half_solve_kpoint';return
#else
      status=HALF_ERROR_UNAVAILABLE;message='CUDA backend is not compiled';return
#endif
    end if
#ifdef HALF_HAVE_MKL
    call prepare(self,kpoint,basis,paw,status,message);if(status/=HALF_SUCCESS)return
    call assemble_dense_gamma(self%veff,basis,paw,h,s)
    if(present(basis_out))basis_out=basis
#else
    status=HALF_ERROR_UNAVAILABLE;message='H/S assembly requires oneMKL';return
#endif
  end subroutine

  subroutine context_apply_hs(self,kpoint,psi,hpsi,spsi,status,message)
    class(half_context_t),intent(inout)::self
    real(dp),intent(in)::kpoint(3)
    complex(dp),intent(in)::psi(:,:)
    complex(dp),allocatable,intent(out)::hpsi(:,:),spsi(:,:)
    integer,intent(out)::status
    character(len=*),intent(out)::message
    complex(dp),allocatable::grid(:),real_grid(:),weighted_grid(:),local_grid(:),coeff(:),dcoeff(:),qcoeff(:)
    integer::n,ngrid,ib,it,iat,a
#ifdef HALF_HAVE_CUDA
    real(dp),device,allocatable::dveff(:)
    complex(dp),device,allocatable::dpsi(:,:),dhpsi(:,:),dspsi(:,:)
    real(dp)::eh,exc
#endif
    if(self%backend==HALF_BACKEND_CUDA)then
#ifdef HALF_HAVE_CUDA
      if(.not.self%apply_ready.or.any(self%apply_kpoint/=kpoint))then
        if(self%gpu_apply_operator%n>0)call destroy_hs_operator_cuda(self%gpu_apply_operator)
        call self%make_basis(kpoint,self%apply_basis,status,message)
        if(status/=HALF_SUCCESS)return
        call build_veff_cuda(self%charge,self%potcars,self%crystal,self%use_pbe,dveff,eh,exc,.false.)
        call prepare_hs_operator_cuda(dveff,self%apply_basis,self%potcars,self%crystal, &
          self%gpu_apply_operator,self%use_uspp,projector_tile_atoms=8)
        deallocate(dveff)
        self%apply_kpoint=kpoint;self%apply_ready=.true.
      end if
      n=int(self%apply_basis%npw)
      if(size(psi,1)/=n)then
        status=HALF_ERROR_INVALID_ARGUMENT;message='state leading dimension does not match plane-wave basis';return
      end if
      allocate(dpsi(n,size(psi,2)),dhpsi(n,size(psi,2)),dspsi(n,size(psi,2)))
      dpsi=psi
      call apply_hs_operator_cuda(self%gpu_apply_operator,dpsi,dhpsi,dspsi)
      allocate(hpsi(n,size(psi,2)),spsi(n,size(psi,2)))
      hpsi=dhpsi;spsi=dspsi
      status=HALF_SUCCESS;message='';return
#else
      status=HALF_ERROR_UNAVAILABLE;message='CUDA backend is not compiled';return
#endif
    end if
#ifdef HALF_HAVE_MKL
    if(.not.self%apply_ready.or.any(self%apply_kpoint/=kpoint))then
      if(allocated(self%apply_paw))deallocate(self%apply_paw)
      call prepare(self,kpoint,self%apply_basis,self%apply_paw,status,message)
      if(status/=HALF_SUCCESS)return
      self%apply_kpoint=kpoint;self%apply_ready=.true.
    end if
    n=int(self%apply_basis%npw);ngrid=size(self%veff)
    if(size(psi,1)/=n)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='state leading dimension does not match plane-wave basis';return
    end if
    allocate(hpsi(n,size(psi,2)),spsi(n,size(psi,2)),grid(ngrid),real_grid(ngrid), &
      weighted_grid(ngrid),local_grid(ngrid))
    do ib=1,size(psi,2)
      grid=(0.0_dp,0.0_dp)
      grid(self%apply_basis%fft_index)=psi(:,ib)
      call fft3_backward(self%apply_basis%shape,grid,real_grid)
      weighted_grid=real_grid*self%veff
      call fft3_forward(self%apply_basis%shape,weighted_grid,local_grid)
      hpsi(:,ib)=self%apply_basis%kinetic*psi(:,ib)+ &
        local_grid(self%apply_basis%fft_index)/real(ngrid,dp)
      spsi(:,ib)=psi(:,ib)
      do it=1,size(self%apply_paw)
        allocate(coeff(self%apply_paw(it)%nlm),dcoeff(self%apply_paw(it)%nlm), &
          qcoeff(self%apply_paw(it)%nlm))
        do iat=1,self%apply_paw(it)%natoms
          do a=1,self%apply_paw(it)%nlm
            coeff(a)=dot_product(self%apply_paw(it)%projectors(iat,a,:),psi(:,ib))
          end do
          dcoeff=matmul(self%apply_paw(it)%dij_atom(iat,:,:),coeff)
          qcoeff=matmul(self%apply_paw(it)%qij,coeff)
          do a=1,self%apply_paw(it)%nlm
            hpsi(:,ib)=hpsi(:,ib)+self%apply_paw(it)%projectors(iat,a,:)*dcoeff(a)
            spsi(:,ib)=spsi(:,ib)+self%apply_paw(it)%projectors(iat,a,:)*qcoeff(a)
          end do
        end do
        deallocate(coeff,dcoeff,qcoeff)
      end do
    end do
    status=HALF_SUCCESS;message=''
#else
    status=HALF_ERROR_UNAVAILABLE;message='CPU matrix-free apply requires oneMKL'
#endif
  end subroutine

  subroutine context_solve_kpoint(self,kpoint,nbands,eigenvalues,eigenvectors,overlap_min,overlap_max,status,message,basis_out)
    class(half_context_t),intent(in)::self
    real(dp),intent(in)::kpoint(3)
    integer,intent(in)::nbands
    real(dp),allocatable,intent(out)::eigenvalues(:)
    complex(dp),allocatable,intent(out),optional::eigenvectors(:,:)
    real(dp),intent(out)::overlap_min,overlap_max
    integer,intent(out)::status
    character(len=*),intent(out)::message
    type(plane_wave_basis_t),intent(out),optional::basis_out
    type(plane_wave_basis_t)::basis
    type(paw_species_t),allocatable::paw(:)
    real(dp),allocatable::all_values(:)
    complex(dp),allocatable::all_vectors(:,:)
    if(self%backend==HALF_BACKEND_CUDA)then
#ifdef HALF_HAVE_CUDA
      call self%make_basis(kpoint,basis,status,message);if(status/=HALF_SUCCESS)return
      if(nbands<1.or.nbands>basis%npw)then
        status=HALF_ERROR_INVALID_ARGUMENT;message='invalid requested band count';return
      end if
      if(present(eigenvectors))then
        call solve_dense_gamma_cuda_full(self%charge,self%potcars,self%crystal,basis,self%use_pbe,all_values, &
          overlap_min,overlap_max,solver=self%solver,use_uspp=self%use_uspp,eigenvectors=all_vectors,target_bands=nbands)
        if(self%solver=='evx')then
          eigenvectors=all_vectors
        else
          eigenvectors=all_vectors(:,:nbands)
        end if
      else
        call solve_dense_gamma_cuda_full(self%charge,self%potcars,self%crystal,basis,self%use_pbe,all_values, &
          overlap_min,overlap_max,solver=self%solver,use_uspp=self%use_uspp,target_bands=nbands)
      end if
      if(self%solver=='evx')then
        eigenvalues=all_values
      else
        eigenvalues=all_values(:nbands)
      end if
      if(present(basis_out))basis_out=basis
      status=HALF_SUCCESS;message='';return
#else
      status=HALF_ERROR_UNAVAILABLE;message='CUDA backend is not compiled';return
#endif
    end if
#ifdef HALF_HAVE_MKL
    call prepare(self,kpoint,basis,paw,status,message);if(status/=HALF_SUCCESS)return
    if(nbands<1.or.nbands>basis%npw)then
      status=HALF_ERROR_INVALID_ARGUMENT;message='invalid requested band count';return
    end if
    if(present(eigenvectors))then
      call solve_dense_gamma(self%veff,basis,paw,all_values,overlap_min,overlap_max,all_vectors)
      eigenvectors=all_vectors(:,:nbands)
    else
      call solve_dense_gamma(self%veff,basis,paw,all_values,overlap_min,overlap_max)
    end if
    eigenvalues=all_values(:nbands)
    if(present(basis_out))basis_out=basis
#else
    status=HALF_ERROR_UNAVAILABLE;message='k-point solve requires oneMKL';return
#endif
  end subroutine
end module half_library
