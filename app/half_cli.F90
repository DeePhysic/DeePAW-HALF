program half_cli
  use half_kinds, only: dp, i64
  use half_types, only: crystal_t, charge_grid_t, plane_wave_basis_t, potcar_t
  use half_chgcar, only: read_chgcar
  use half_potcar, only: read_potcar,validate_potcar_structure
  use half_basis, only: build_plane_wave_basis
  use half_kpoints, only: kpoint_set_t,read_explicit_kpoints,gamma_centered_mesh,gamma_centered_irreducible_mesh,generate_cubic_band_path
  use half_paw, only: paw_species_t, build_paw_operators
  use half_uspp, only: build_uspp_dij_cpu
  use half_energy, only: compute_occupations,ewald_energy,read_vasp_eigenval
  use half_math,only:inverse3
  use half_cpu, only: cpu_density_mean
#ifdef HALF_CLI_HAVE_MKL
  use half_potential, only: build_veff_lda, build_veff_pbe
  use half_dense_solver, only: solve_dense_gamma
#endif
#ifdef HALF_CLI_CUDA
  use half_cuda_solver, only: solve_dense_gamma_cuda_full
#endif
#ifdef HALF_CLI_HAVE_HDF5
  use half_vaspwave,only:wave_block_t,write_vaspwave_h5
#endif
  implicit none
  character(len=1024) :: invocation, command
  integer :: argument_offset

  call get_command_argument(0, invocation)
  invocation = basename(trim(invocation))
  argument_offset = 1
  select case (trim(invocation))
  case ('half-validate-gamma')
    command = 'gamma'; argument_offset = 0
  case ('half-bands')
    command = 'bands'; argument_offset = 0
  case ('half-energy')
    command = 'energy'; argument_offset = 0
  case default
    if (command_argument_count() == 0) then
      call print_help(6); stop
    end if
    call get_command_argument(1, command)
  end select

  select case (trim(command))
  case ('help', '--help', '-h')
    call print_help(6)
  case ('version', '--version', '-V')
    write(*,'(A)') 'DeePAW-HALF 0.4.0'
  case ('gamma', 'validate-gamma')
    call command_gamma(argument_offset)
  case ('inspect')
    call command_inspect(argument_offset)
  case ('potcar')
    call command_potcar(argument_offset)
  case ('paw')
    call command_paw(argument_offset)
  case ('kpoints')
    call command_kpoints(argument_offset)
  case ('bands')
    call command_bands(argument_offset)
  case ('energy')
    call command_energy(argument_offset)
  case default
    call fail('unknown command: '//trim(command)//'; run half --help')
  end select

contains
  subroutine print_help(unit)
    integer,intent(in)::unit
    write(unit,'(A)') 'DeePAW-HALF - Harris fixed-density reconstruction in Fortran'
    write(unit,'(A)') ''
    write(unit,'(A)') 'Usage: half COMMAND [OPTIONS]'
    write(unit,'(A)') ''
    write(unit,'(A)') 'Commands:'
    write(unit,'(A)') '  gamma    validate the Gamma-point eigenspectrum (HAPPY-compatible)'
    write(unit,'(A)') '  inspect  inspect a CHGCAR and construct its Gamma plane-wave basis'
    write(unit,'(A)') '  potcar   inspect POTCAR datasets'
    write(unit,'(A)') '  paw      inspect reciprocal-space PAW operators'
    write(unit,'(A)') '  kpoints  generate a Gamma-centered full or irreducible KSPACING mesh'
    write(unit,'(A)') '  bands    reconstruct bands at arbitrary reciprocal-space k points'
    write(unit,'(A)') '  energy   evaluate the HAPPY-compatible fixed-density total energy'
    write(unit,'(A)') '  version  print the HALF version'
    write(unit,'(A)') ''
    write(unit,'(A)') 'Run half COMMAND --help for command-specific options.'
  end subroutine

  subroutine command_kpoints(offset)
    integer,intent(in)::offset
    type(crystal_t)::crystal
    type(charge_grid_t)::charge
    type(kpoint_set_t)::set
    character(len=1024)::path,arg,value
    real(dp)::kspacing,symprec
    integer::i,narg,ios,k
    logical::full_mesh
    narg=command_argument_count()
    if(narg<offset+1)call fail('kpoints requires CHGCAR')
    call get_command_argument(offset+1,path);kspacing=0.5_dp;symprec=1e-5_dp;full_mesh=.false.;i=offset+2
    do while(i<=narg)
      call get_command_argument(i,arg)
      select case(trim(arg))
      case('--kspacing')
        call option_value(i,narg,'--kspacing',value);read(value,*,iostat=ios)kspacing
        if(ios/=0.or.kspacing<=0)call fail('invalid --kspacing value')
      case('--symprec')
        call option_value(i,narg,'--symprec',value);read(value,*,iostat=ios)symprec
        if(ios/=0.or.symprec<=0)call fail('invalid --symprec value')
      case('--full')
        full_mesh=.true.
      case('--help','-h')
        write(*,'(A)')'Usage: half kpoints CHGCAR [--kspacing VALUE] [--symprec VALUE] [--full]';return
      case default
        call fail('unknown kpoints option: '//trim(arg))
      end select
      i=i+1
    end do
    call read_chgcar(trim(path),crystal,charge)
    if(full_mesh)then;call gamma_centered_mesh(crystal,kspacing,set)
    else;call gamma_centered_irreducible_mesh(crystal,kspacing,set,symprec,.true.);end if
    write(*,'(A)')'{'
    write(*,'(A,I0,A,I0,A,I0,A)')'  "divisions": [',set%divisions(1),', ', &
      set%divisions(2),', ',set%divisions(3),'],'
    write(*,'(A,I0,A)')'  "full_kpoint_count": ',set%full_count,','
    write(*,'(A,I0,A)')'  "kpoint_count": ',set%nk,','
    write(*,'(A)')'  "kpoints": ['
    do k=1,set%nk
      write(*,'(A,3(ES24.16,A),ES24.16,A)',advance='no')'    [',set%points(k,1),', ', &
        set%points(k,2),', ',set%points(k,3),', ',set%weights(k),']'
      if(k<set%nk)then;write(*,'(A)')',';else;write(*,'(A)')'';end if
    end do
    write(*,'(A)')'  ]';write(*,'(A)')'}'
  end subroutine

  subroutine print_gamma_help(unit)
    integer,intent(in)::unit
    write(unit,'(A)') 'Usage: half gamma CHARGE POTENTIAL [OPTIONS]'
    write(unit,'(A)') '       half-validate-gamma CHARGE POTENTIAL [OPTIONS]'
    write(unit,'(A)') ''
    write(unit,'(A)') 'Options:'
    write(unit,'(A)') '  --encut EV                 plane-wave cutoff (default: 400)'
    write(unit,'(A)') '  --bands N                  eigenvalues to report (default: 8)'
    write(unit,'(A)') '  --kpoint KX KY KZ          fractional reciprocal k point (default: Gamma)'
    write(unit,'(A)') '  --xc lda|pbe               exchange-correlation model (default: pbe)'
    write(unit,'(A)') '  --backend auto|cpu|cuda    compute backend (default: auto)'
    write(unit,'(A)') '  --solver evd|evj           CUDA eigensolver: divide-and-conquer or Jacobi (default: evd)'
    write(unit,'(A)') '  --uspp-dij                 add potential-dependent MIMIC_US QDEP correction'
    write(unit,'(A)') '  --reference-eigenval FILE  compare against a VASP EIGENVAL file'
    write(unit,'(A)') '  --output FILE              JSON result (default: gamma_validation.json)'
  end subroutine

  subroutine command_gamma(offset)
    integer,intent(in)::offset
    type(crystal_t)::crystal
    type(charge_grid_t)::rho
    type(plane_wave_basis_t)::basis
    type(potcar_t),allocatable::potcars(:)
    type(paw_species_t),allocatable::paw(:)
    real(dp),allocatable::veff(:),eigenvalues(:),reference(:)
    real(dp)::eh,exc,smin,smax,encut,elapsed,potential_seconds,assembly_seconds,gpu_seconds,kpoint(3)
    logical::use_uspp
    integer::nbands,i,narg,ios,tick0,tick1,rate,unit
    character(len=1024)::charge_path,potential_path,arg,value,output_path,reference_path
    character(len=16)::xc,backend,actual_backend,solver,actual_solver
    narg=command_argument_count()
    if(narg>=offset+1)then
      call get_command_argument(offset+1,arg)
      if(trim(arg)=='--help'.or.trim(arg)=='-h')then; call print_gamma_help(6); return; end if
    end if
#if !defined(HALF_CLI_HAVE_MKL) && !defined(HALF_CLI_CUDA)
    call fail('gamma requires either CUDA or oneMKL; enable a numerical backend and rebuild HALF')
#else
    if(narg<offset+2)then; call print_gamma_help(0); call fail('gamma requires CHARGE and POTENTIAL'); end if
    call get_command_argument(offset+1,charge_path)
    call get_command_argument(offset+2,potential_path)
    encut=400.0_dp; nbands=8; xc='pbe'; backend='auto'; solver='evd';use_uspp=.false.;kpoint=0.0_dp
    output_path='gamma_validation.json'; reference_path=''
    i=offset+3
    do while(i<=narg)
      call get_command_argument(i,arg)
      select case(trim(arg))
      case('--encut')
        call option_value(i,narg,'--encut',value); read(value,*,iostat=ios)encut
        if(ios/=0.or.encut<=0.0_dp)call fail('invalid --encut value')
      case('--bands')
        call option_value(i,narg,'--bands',value); read(value,*,iostat=ios)nbands
        if(ios/=0.or.nbands<1)call fail('invalid --bands value')
      case('--kpoint')
        if(i+3>narg)call fail('--kpoint requires KX KY KZ')
        call get_command_argument(i+1,value);read(value,*,iostat=ios)kpoint(1)
        if(ios/=0)call fail('invalid --kpoint KX value')
        call get_command_argument(i+2,value);read(value,*,iostat=ios)kpoint(2)
        if(ios/=0)call fail('invalid --kpoint KY value')
        call get_command_argument(i+3,value);read(value,*,iostat=ios)kpoint(3)
        if(ios/=0)call fail('invalid --kpoint KZ value');i=i+3
      case('--xc')
        call option_value(i,narg,'--xc',xc); xc=lower(trim(xc))
        if(trim(xc)/='lda'.and.trim(xc)/='pbe')call fail('--xc must be lda or pbe')
      case('--backend')
        call option_value(i,narg,'--backend',backend); backend=lower(trim(backend))
        if(trim(backend)/='auto'.and.trim(backend)/='cpu'.and.trim(backend)/='cuda') &
          call fail('--backend must be auto, cpu, or cuda')
      case('--solver')
        call option_value(i,narg,'--solver',solver); solver=lower(trim(solver))
        if(trim(solver)/='evd'.and.trim(solver)/='evj')call fail('--solver must be evd or evj')
      case('--output')
        call option_value(i,narg,'--output',output_path)
      case('--reference-eigenval')
        call option_value(i,narg,'--reference-eigenval',reference_path)
      case('--uspp-dij')
        use_uspp=.true.
      case('--help','-h')
        call print_gamma_help(6); return
      case default
        call fail('unknown gamma option: '//trim(arg))
      end select
      i=i+1
    end do

    call system_clock(tick0,rate)
    call read_chgcar(trim(charge_path),crystal,rho)
    call read_potcar(trim(potential_path),potcars)
    call validate_potcar_structure(potcars,crystal)
    call build_plane_wave_basis(crystal,rho%shape,encut,kpoint,basis)
    potential_seconds=0.0_dp;assembly_seconds=0.0_dp;gpu_seconds=0.0_dp
#ifdef HALF_CLI_CUDA
    if(trim(backend)=='auto'.or.trim(backend)=='cuda')then
      actual_backend='cuda'
      actual_solver=solver
      call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',eigenvalues,smin,smax, &
        potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp)
    else
      actual_backend='cpu'
      actual_solver='evd'
      if(trim(solver)/='evd')call fail('--solver evj is available only with --backend cuda')
#ifndef HALF_CLI_HAVE_MKL
      call fail('the CPU gamma backend requires oneMKL; use --backend cuda or rebuild with MKLROOT')
#else
      call build_paw_operators(potcars,crystal,basis,paw)
      if(trim(xc)=='pbe')then;call build_veff_pbe(rho,potcars,crystal,veff,eh,exc)
      else;call build_veff_lda(rho,potcars,crystal,veff,eh,exc);end if
      if(use_uspp)call build_uspp_dij_cpu(veff,rho%shape,potcars,crystal,paw)
      call solve_dense_gamma(veff,basis,paw,eigenvalues,smin,smax)
#endif
    end if
#else
    if(trim(backend)=='cuda')call fail('this HALF build has no CUDA backend')
    actual_backend='cpu'
    actual_solver='evd'
    if(trim(solver)/='evd')call fail('--solver evj is available only with --backend cuda')
    call build_paw_operators(potcars,crystal,basis,paw)
    if(trim(xc)=='pbe')then;call build_veff_pbe(rho,potcars,crystal,veff,eh,exc)
    else;call build_veff_lda(rho,potcars,crystal,veff,eh,exc);end if
    if(use_uspp)call build_uspp_dij_cpu(veff,rho%shape,potcars,crystal,paw)
    call solve_dense_gamma(veff,basis,paw,eigenvalues,smin,smax)
#endif
    call system_clock(tick1)
    elapsed=real(tick1-tick0,dp)/real(rate,dp)
    nbands=min(nbands,size(eigenvalues))
    if(len_trim(reference_path)>0)then
      call read_eigenval_gamma(trim(reference_path),nbands,reference)
    else
      allocate(reference(0))
    end if
    if(trim(output_path)/='-')then
      open(newunit=unit,file=trim(output_path),status='replace',action='write',iostat=ios)
      if(ios/=0)call fail('cannot write output: '//trim(output_path))
      call write_gamma_json(unit,charge_path,potential_path,xc,actual_backend,encut,basis%npw, &
        eigenvalues(:nbands),reference,smin,smax,elapsed,potential_seconds,assembly_seconds,gpu_seconds,actual_solver,use_uspp,kpoint)
      close(unit)
    end if
    call write_gamma_json(6,charge_path,potential_path,xc,actual_backend,encut,basis%npw, &
      eigenvalues(:nbands),reference,smin,smax,elapsed,potential_seconds,assembly_seconds,gpu_seconds,actual_solver,use_uspp,kpoint)
#endif
  end subroutine

  subroutine write_gamma_json(unit,charge,potential,xc,backend,encut,npw,eigenvalues,reference, &
      smin,smax,elapsed,potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp,kpoint)
    integer,intent(in)::unit
    integer(i64),intent(in)::npw
    character(len=*),intent(in)::charge,potential,xc,backend,solver
    logical,intent(in)::use_uspp
    real(dp),intent(in)::encut,eigenvalues(:),reference(:),smin,smax,elapsed,potential_seconds,assembly_seconds,gpu_seconds,kpoint(3)
    real(dp),allocatable::errors(:)
    integer::n
    n=min(size(eigenvalues),size(reference)); allocate(errors(n))
    if(n>0)errors=eigenvalues(:n)-reference(:n)
    write(unit,'(A)') '{'
    write(unit,'(A)') '  "implementation": "DeePAW-HALF 0.4.0",'
    write(unit,'(A,A,A)') '  "charge": "',trim(charge),'",'
    write(unit,'(A,A,A)') '  "potential": "',trim(potential),'",'
    write(unit,'(A,A,A)') '  "xc": "',trim(xc),'",'
    write(unit,'(A,A,A)') '  "backend": "',trim(backend),'",'
    write(unit,'(A,A,A)') '  "solver": "',trim(solver),'",'
    write(unit,'(A,A,A)') '  "uspp_dij": ',merge('true ','false',use_uspp),','
    write(unit,'(A,ES24.16,A,ES24.16,A,ES24.16,A)') '  "kpoint_fractional": [', &
      kpoint(1),', ',kpoint(2),', ',kpoint(3),'],'
    write(unit,'(A,ES24.16,A)') '  "encut_eV": ',encut,','
    write(unit,'(A,I0,A)') '  "plane_waves": ',npw,','
    call write_real_array(unit,'gamma_eigenvalues_eV',eigenvalues,.true.)
    call write_real_array(unit,'reference_gamma_eigenvalues_eV',reference,.true.)
    call write_real_array(unit,'gamma_errors_eV',errors,.true.)
    if(n>0)then
      write(unit,'(A,ES24.16,A)') '  "gamma_max_abs_error_eV": ',maxval(abs(errors)),','
      write(unit,'(A,ES24.16,A)') '  "gamma_rms_error_eV": ',sqrt(sum(errors**2)/real(n,dp)),','
    else
      write(unit,'(A)') '  "gamma_max_abs_error_eV": null,'
      write(unit,'(A)') '  "gamma_rms_error_eV": null,'
    end if
    write(unit,'(A,ES24.16,A)') '  "overlap_eigenvalue_min": ',smin,','
    write(unit,'(A,ES24.16,A)') '  "overlap_eigenvalue_max": ',smax,','
    write(unit,'(A,ES24.16,A)') '  "elapsed_seconds": ',elapsed,','
    write(unit,'(A,ES24.16,A)') '  "potential_seconds": ',potential_seconds,','
    write(unit,'(A,ES24.16,A)') '  "assembly_seconds": ',assembly_seconds,','
    write(unit,'(A,ES24.16)') '  "gpu_solver_seconds": ',gpu_seconds
    write(unit,'(A)') '}'
  end subroutine

  subroutine write_real_array(unit,name,values,comma)
    integer,intent(in)::unit
    character(len=*),intent(in)::name
    real(dp),intent(in)::values(:)
    logical,intent(in)::comma
    integer::i
    write(unit,'(A,A,A)',advance='no') '  "',trim(name),'": ['
    do i=1,size(values)
      if(i>1)write(unit,'(A)',advance='no') ', '
      write(unit,'(ES24.16)',advance='no') values(i)
    end do
    if(comma)then; write(unit,'(A)') '],'; else; write(unit,'(A)') ']'; end if
  end subroutine

  subroutine read_eigenval_gamma(path,wanted,values)
    character(len=*),intent(in)::path
    integer,intent(in)::wanted
    real(dp),allocatable,intent(out)::values(:)
    character(len=2048)::line
    integer::unit,ios,i,nelect,nkpts,nbands,index
    real(dp)::energy,occupation
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)call fail('cannot open reference EIGENVAL: '//trim(path))
    do i=1,5; read(unit,'(A)',iostat=ios)line; if(ios/=0)call fail('invalid EIGENVAL header'); end do
    read(unit,*,iostat=ios)nelect,nkpts,nbands
    if(ios/=0.or.nbands<1)call fail('invalid EIGENVAL dimensions')
    read(unit,'(A)',iostat=ios)line
    read(unit,'(A)',iostat=ios)line
    allocate(values(min(wanted,nbands)))
    do i=1,size(values)
      read(unit,*,iostat=ios)index,energy,occupation
      if(ios/=0)call fail('cannot read Gamma bands from EIGENVAL')
      values(i)=energy
    end do
    close(unit)
  end subroutine

  subroutine command_inspect(offset)
    integer,intent(in)::offset
    type(crystal_t)::crystal
    type(charge_grid_t)::charge
    type(plane_wave_basis_t)::basis
    character(len=1024)::path,arg,value
    integer::narg,i,ios
    real(dp)::encut
    narg=command_argument_count()
    if(narg>=offset+1)then
      call get_command_argument(offset+1,arg)
      if(trim(arg)=='--help'.or.trim(arg)=='-h')then
        write(*,'(A)')'Usage: half inspect CHARGE [--encut EV]'; return
      end if
    end if
    if(narg<offset+1)call fail('inspect requires CHARGE')
    call get_command_argument(offset+1,path); encut=400.0_dp; i=offset+2
    do while(i<=narg)
      call get_command_argument(i,arg)
      if(trim(arg)/='--encut')call fail('unknown inspect option: '//trim(arg))
      call option_value(i,narg,'--encut',value); read(value,*,iostat=ios)encut
      if(ios/=0.or.encut<=0.0_dp)call fail('invalid --encut value'); i=i+1
    end do
    call read_chgcar(trim(path),crystal,charge)
    call build_plane_wave_basis(crystal,charge%shape,encut,[0.0_dp,0.0_dp,0.0_dp],basis)
    write(*,'(A)') 'DeePAW-HALF inspect'
    write(*,'(A,A)') 'system: ',trim(crystal%system_name)
    write(*,'(A,I0)') 'ions: ',crystal%nions
    write(*,'(A,3(I0,1X))') 'grid: ',charge%shape
    write(*,'(A,F0.10)') 'volume_A3: ',crystal%volume
    write(*,'(A,F0.10)') 'smooth_electrons: ',cpu_density_mean(charge%values)
    write(*,'(A,F0.3)') 'encut_eV: ',encut
    write(*,'(A,I0)') 'gamma_plane_waves: ',basis%npw
  end subroutine

  subroutine command_potcar(offset)
    integer,intent(in)::offset
    type(potcar_t),allocatable::p(:)
    character(len=1024)::path,arg
    integer::i
    if(command_argument_count()<offset+1)call fail('potcar requires POTENTIAL')
    call get_command_argument(offset+1,arg)
    if(trim(arg)=='--help'.or.trim(arg)=='-h')then; write(*,'(A)')'Usage: half potcar POTENTIAL'; return; end if
    path=arg; call read_potcar(trim(path),p)
    write(*,'(A,I0)')'datasets: ',size(p)
    do i=1,size(p)
      write(*,'(A,I0,A,A,A,F0.6,A,F0.3,A,I0)')'dataset ',i,': element=',trim(p(i)%element), &
        ' zval=',p(i)%zval,' enmax_eV=',p(i)%enmax,' channels=',p(i)%channels
    end do
  end subroutine

  subroutine command_paw(offset)
    integer,intent(in)::offset
    type(crystal_t)::crystal
    type(charge_grid_t)::rho
    type(plane_wave_basis_t)::basis
    type(potcar_t),allocatable::p(:)
    type(paw_species_t),allocatable::paw(:)
    character(len=1024)::chg,pot,arg,value
    real(dp)::encut
    integer::narg,i,ios
    narg=command_argument_count()
    if(narg<offset+2)call fail('paw requires CHARGE and POTENTIAL')
    call get_command_argument(offset+1,arg)
    if(trim(arg)=='--help'.or.trim(arg)=='-h')then
      write(*,'(A)')'Usage: half paw CHARGE POTENTIAL [--encut EV]'; return
    end if
    chg=arg; call get_command_argument(offset+2,pot); encut=400.0_dp; i=offset+3
    do while(i<=narg)
      call get_command_argument(i,arg)
      if(trim(arg)/='--encut')call fail('unknown paw option: '//trim(arg))
      call option_value(i,narg,'--encut',value); read(value,*,iostat=ios)encut
      if(ios/=0.or.encut<=0.0_dp)call fail('invalid --encut value'); i=i+1
    end do
    call read_chgcar(trim(chg),crystal,rho); call read_potcar(trim(pot),p);call validate_potcar_structure(p,crystal)
    call build_plane_wave_basis(crystal,rho%shape,encut,[0.0_dp,0.0_dp,0.0_dp],basis)
    call build_paw_operators(p,crystal,basis,paw)
    write(*,'(A,I0)')'plane_waves: ',basis%npw
    do i=1,size(paw)
      write(*,'(A,I0,A,I0,A,ES16.8,A,ES16.8)')'species ',i,': projectors=',paw(i)%nlm, &
        ' projector_norm2=',sum(abs(paw(i)%projectors)**2),' qij_sum=',sum(paw(i)%qij)
    end do
  end subroutine

  subroutine command_bands(offset)
    integer,intent(in)::offset
    type(crystal_t)::crystal
    type(charge_grid_t)::rho
    type(plane_wave_basis_t)::basis
    type(potcar_t),allocatable::potcars(:)
    type(paw_species_t),allocatable::paw(:)
    type(kpoint_set_t)::set
    real(dp),allocatable::veff(:),values(:),eigenvalues(:,:)
    complex(dp),allocatable::vectors(:,:)
    integer(i64),allocatable::plane_waves(:)
    real(dp)::encut,eh,exc,smin,smax,potential_seconds,assembly_seconds,gpu_seconds
    integer::narg,i,ios,ik,nbands,unit,npoints
    character(len=1024)::charge_path,potential_path,kpoints_path,arg,value,output_path,hdf5_path,path_spec,path_used
    character(len=16)::xc,backend,solver,actual_backend
    logical::use_uspp
#ifdef HALF_CLI_HAVE_HDF5
    type(plane_wave_basis_t),allocatable::wave_bases(:)
    type(wave_block_t),allocatable::waves(:)
    real(dp),allocatable::wave_occ(:,:)
    real(dp)::wave_mu,wave_band,wave_entropy,nelect
#endif
    narg=command_argument_count()
    if(narg>=offset+1)then
      call get_command_argument(offset+1,arg)
      if(trim(arg)=='--help'.or.trim(arg)=='-h')then;call print_bands_help(6);return;end if
    end if
#if !defined(HALF_CLI_HAVE_MKL) && !defined(HALF_CLI_CUDA)
    call fail('bands requires either CUDA or oneMKL; enable a numerical backend and rebuild HALF')
#else
    if(narg<offset+2)then;call print_bands_help(0);call fail('bands requires CHARGE and POTENTIAL');end if
    call get_command_argument(offset+1,charge_path);call get_command_argument(offset+2,potential_path)
    kpoints_path='';path_spec='';path_used='';i=offset+3
    if(i<=narg)then;call get_command_argument(i,arg);if(arg(1:1)/='-')then;kpoints_path=arg;i=i+1;end if;end if
    encut=400.0_dp;nbands=8;npoints=60;xc='pbe';backend='auto';solver='evd';use_uspp=.true.;output_path='bands.json';hdf5_path=''
    do while(i<=narg)
      call get_command_argument(i,arg)
      select case(trim(arg))
      case('--encut');call option_value(i,narg,'--encut',value);read(value,*,iostat=ios)encut
        if(ios/=0.or.encut<=0)call fail('invalid --encut value')
      case('--bands');call option_value(i,narg,'--bands',value);read(value,*,iostat=ios)nbands
        if(ios/=0.or.nbands<1)call fail('invalid --bands value')
      case('--npoints');call option_value(i,narg,'--npoints',value);read(value,*,iostat=ios)npoints
        if(ios/=0.or.npoints<2)call fail('invalid --npoints value')
      case('--path');call option_value(i,narg,'--path',path_spec)
      case('--xc');call option_value(i,narg,'--xc',xc);xc=lower(trim(xc))
        if(trim(xc)/='lda'.and.trim(xc)/='pbe')call fail('--xc must be lda or pbe')
      case('--backend');call option_value(i,narg,'--backend',backend);backend=lower(trim(backend))
        if(trim(backend)/='auto'.and.trim(backend)/='cpu'.and.trim(backend)/='cuda')call fail('--backend must be auto, cpu, or cuda')
      case('--solver');call option_value(i,narg,'--solver',solver);solver=lower(trim(solver))
        if(trim(solver)/='evd'.and.trim(solver)/='evj')call fail('--solver must be evd or evj')
      case('--output');call option_value(i,narg,'--output',output_path)
      case('--output-prefix');call option_value(i,narg,'--output-prefix',value);output_path=trim(value)//'.json'
      case('--vaspwave-h5');call option_value(i,narg,'--vaspwave-h5',hdf5_path)
      case('--uspp-dij');use_uspp=.true.
      case('--no-uspp-dij');use_uspp=.false.
      case('--help','-h');call print_bands_help(6);return
      case default;call fail('unknown bands option: '//trim(arg))
      end select
      i=i+1
    end do
    call read_chgcar(trim(charge_path),crystal,rho);call read_potcar(trim(potential_path),potcars);call validate_potcar_structure(potcars,crystal)
    if(len_trim(kpoints_path)>0)then;call read_explicit_kpoints(trim(kpoints_path),set);path_used=kpoints_path
    else;call generate_cubic_band_path(crystal,npoints,path_spec,set,path_used);end if
#ifdef HALF_CLI_CUDA
    if(trim(backend)=='auto'.or.trim(backend)=='cuda')then;actual_backend='cuda'
    else;actual_backend='cpu';end if
#else
    if(trim(backend)=='cuda')call fail('this HALF build has no CUDA backend');actual_backend='cpu'
#endif
#ifndef HALF_CLI_HAVE_MKL
    if(trim(actual_backend)=='cpu')call fail('the CPU bands backend requires oneMKL')
#endif
    if(trim(actual_backend)=='cpu'.and.trim(solver)/='evd')call fail('--solver evj is available only with CUDA')
    if(len_trim(hdf5_path)>0.and.trim(solver)/='evd')call fail('vaspwave.h5 export requires --solver evd')
#ifndef HALF_CLI_HAVE_HDF5
    if(len_trim(hdf5_path)>0)call fail('this HALF build has no HDF5 support')
#endif
    allocate(eigenvalues(set%nk,nbands),plane_waves(set%nk))
#ifdef HALF_CLI_HAVE_HDF5
    if(len_trim(hdf5_path)>0)allocate(wave_bases(set%nk),waves(set%nk))
#endif
    do ik=1,set%nk
      call build_plane_wave_basis(crystal,rho%shape,encut,set%points(ik,:),basis);plane_waves(ik)=basis%npw
      if(nbands>basis%npw)call fail('requested bands exceed plane waves at a k point')
#ifdef HALF_CLI_CUDA
      if(trim(actual_backend)=='cuda')then
        if(len_trim(hdf5_path)>0)then
          call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax, &
            potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp,eigenvectors=vectors)
        else
          call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax, &
            potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp)
        end if
      else
#endif
#ifdef HALF_CLI_HAVE_MKL
        call build_paw_operators(potcars,crystal,basis,paw)
        if(trim(xc)=='pbe')then;call build_veff_pbe(rho,potcars,crystal,veff,eh,exc)
        else;call build_veff_lda(rho,potcars,crystal,veff,eh,exc);end if
        if(use_uspp)call build_uspp_dij_cpu(veff,rho%shape,potcars,crystal,paw)
        if(len_trim(hdf5_path)>0)then;call solve_dense_gamma(veff,basis,paw,values,smin,smax,vectors)
        else;call solve_dense_gamma(veff,basis,paw,values,smin,smax);end if
#endif
#ifdef HALF_CLI_CUDA
      end if
#endif
      eigenvalues(ik,:)=values(:nbands)
#ifdef HALF_CLI_HAVE_HDF5
      if(len_trim(hdf5_path)>0)then;wave_bases(ik)=basis;waves(ik)%coefficients=vectors(:,:nbands);end if
#endif
      write(0,'(A,I0,A,I0,A,ES14.6)')'HALF bands: k point ',ik,'/',set%nk,' E1_eV=',values(1)
    end do
#ifdef HALF_CLI_HAVE_HDF5
    if(len_trim(hdf5_path)>0)then
      nelect=0;do i=1,size(potcars);nelect=nelect+potcars(i)%zval*real(crystal%counts(i),dp);end do
      call compute_occupations(eigenvalues,set%weights,nelect,0.0_dp,wave_occ,wave_mu,wave_band,wave_entropy)
      call write_vaspwave_h5(trim(hdf5_path),crystal,rho,encut,set%points,eigenvalues,wave_occ,wave_mu,wave_bases,waves)
    end if
#endif
    if(trim(output_path)=='-')then;unit=6
    else;open(newunit=unit,file=trim(output_path),status='replace',action='write',iostat=ios)
      if(ios/=0)call fail('cannot write output: '//trim(output_path));end if
    call write_bands_json(unit,charge_path,potential_path,path_used,xc,actual_backend,solver,encut,set,eigenvalues,plane_waves,use_uspp)
    if(unit/=6)then;close(unit);call write_bands_json(6,charge_path,potential_path,path_used,xc,actual_backend,solver,encut,set,eigenvalues,plane_waves,use_uspp);end if
#endif
  end subroutine

  subroutine print_bands_help(unit)
    integer,intent(in)::unit
    write(unit,'(A)')'Usage: half bands CHARGE POTENTIAL [KPOINTS] [OPTIONS]'
    write(unit,'(A)')'       half-bands CHARGE POTENTIAL [KPOINTS] [OPTIONS]'
    write(unit,'(A)')'Options: --encut EV --bands N --npoints N --path LABELS --xc lda|pbe --backend auto|cpu|cuda'
    write(unit,'(A)')'         --solver evd|evj --no-uspp-dij --output FILE --output-prefix PREFIX --vaspwave-h5 FILE'
    write(unit,'(A)')'KPOINTS is a VASP explicit reciprocal-coordinate file; omit it for an automatic cubic band path.'
  end subroutine

  subroutine write_bands_json(unit,charge,potential,kfile,xc,backend,solver,encut,set,eigenvalues,plane_waves,use_uspp)
    integer,intent(in)::unit
    character(len=*),intent(in)::charge,potential,kfile,xc,backend,solver
    real(dp),intent(in)::encut,eigenvalues(:,:)
    type(kpoint_set_t),intent(in)::set
    integer(i64),intent(in)::plane_waves(:)
    logical,intent(in)::use_uspp
    integer::ik,ib
    write(unit,'(A)')'{';write(unit,'(A)')'  "implementation": "DeePAW-HALF 0.4.0",'
    write(unit,'(A,A,A)')'  "charge": "',trim(charge),'",';write(unit,'(A,A,A)')'  "potential": "',trim(potential),'",'
    write(unit,'(A,A,A)')'  "kpoints_source": "',trim(kfile),'",';write(unit,'(A,A,A)')'  "xc": "',trim(xc),'",'
    write(unit,'(A,A,A)')'  "backend": "',trim(backend),'",';write(unit,'(A,A,A)')'  "solver": "',trim(solver),'",'
    write(unit,'(A,A,A)')'  "uspp_dij": ',merge('true ','false',use_uspp),',';write(unit,'(A,ES24.16,A)')'  "encut_eV": ',encut,','
    write(unit,'(A)')'  "kpoints": ['
    do ik=1,set%nk
      write(unit,'(A,3(ES24.16,A),ES24.16,A)',advance='no')'    [',set%points(ik,1),', ',set%points(ik,2),', ',set%points(ik,3),', ',set%weights(ik),']'
      if(ik<set%nk)then;write(unit,'(A)')',';else;write(unit,'(A)')'';end if
    end do
    write(unit,'(A)')'  ],';write(unit,'(A)',advance='no')'  "plane_waves": ['
    do ik=1,size(plane_waves);if(ik>1)write(unit,'(A)',advance='no')', ';write(unit,'(I0)',advance='no')plane_waves(ik);end do
    write(unit,'(A)')'],';write(unit,'(A)')'  "eigenvalues_eV": ['
    do ik=1,size(eigenvalues,1)
      write(unit,'(A)',advance='no')'    ['
      do ib=1,size(eigenvalues,2);if(ib>1)write(unit,'(A)',advance='no')', ';write(unit,'(ES24.16)',advance='no')eigenvalues(ik,ib);end do
      if(ik<size(eigenvalues,1))then;write(unit,'(A)')'],';else;write(unit,'(A)')']';end if
    end do
    write(unit,'(A)')'  ]';write(unit,'(A)')'}'
  end subroutine

  subroutine command_energy(offset)
    integer,intent(in)::offset
    type(crystal_t)::crystal
    type(charge_grid_t)::rho
    type(plane_wave_basis_t)::basis
    type(potcar_t),allocatable::potcars(:)
    type(paw_species_t),allocatable::paw(:)
    type(kpoint_set_t)::set,reference_set
    real(dp),allocatable::veff(:),values(:),eigenvalues(:,:),occupation(:,:),charges(:),reference_eigenvalues(:,:)
    complex(dp),allocatable::vectors(:,:)
    integer(i64),allocatable::plane_waves(:)
    real(dp)::encut,kspacing,symprec,sigma,nelect,density_electrons,mu,band_energy,entropy_term,reference_nelect
    real(dp)::eh,exc,exv,smin,smax,potential_seconds,assembly_seconds,gpu_seconds
    real(dp)::ewald,ewald_real,ewald_recip,ewald_self,ewald_background,atomic_reference,local_g0
    real(dp)::paw_atomic,internal_energy,free_energy
    real(dp)::force_step,displaced_energies(2),inverse_lattice(3,3),delta_cart(3)
    real(dp),allocatable::forces(:,:)
    type(crystal_t)::displaced
    integer::narg,i,ios,ik,it,iat,ion,nbands,minimum_bands,unit
    character(len=1024)::charge_path,potential_path,kpoints_path,arg,value,output_path,hdf5_path,reference_path
    character(len=16)::xc,backend,solver,actual_backend
    logical::use_uspp,full_mesh,have_atomic_override,have_paw_atomic,do_forces,use_reference
#ifdef HALF_CLI_HAVE_HDF5
    type(plane_wave_basis_t),allocatable::wave_bases(:)
    type(wave_block_t),allocatable::waves(:)
#endif
    if(command_argument_count()>=offset+1)then
      call get_command_argument(offset+1,arg)
      if(trim(arg)=='--help'.or.trim(arg)=='-h')then
        call print_energy_help(6);return
      end if
    end if
#if !defined(HALF_CLI_HAVE_MKL) && !defined(HALF_CLI_CUDA)
    call fail('energy requires either CUDA or oneMKL; enable a numerical backend and rebuild HALF')
#else
    narg=command_argument_count();if(narg<offset+2)then;call print_energy_help(0);call fail('energy requires CHARGE and POTENTIAL');end if
    call get_command_argument(offset+1,charge_path);call get_command_argument(offset+2,potential_path)
    encut=400.0_dp;kspacing=0.5_dp;symprec=1e-5_dp;sigma=0.0_dp;nbands=0
    xc='pbe';backend='auto';solver='evd';kpoints_path='';output_path='energy.json';hdf5_path='';reference_path=''
    use_uspp=.true.;full_mesh=.false.;have_atomic_override=.false.;have_paw_atomic=.false.;do_forces=.false.
    paw_atomic=0;atomic_reference=0;force_step=1e-3_dp;i=offset+3
    do while(i<=narg)
      call get_command_argument(i,arg)
      select case(trim(arg))
      case('--encut');call option_value(i,narg,'--encut',value);read(value,*,iostat=ios)encut
        if(ios/=0.or.encut<=0)call fail('invalid --encut value')
      case('--kspacing');call option_value(i,narg,'--kspacing',value);read(value,*,iostat=ios)kspacing
        if(ios/=0.or.kspacing<=0)call fail('invalid --kspacing value')
      case('--symprec');call option_value(i,narg,'--symprec',value);read(value,*,iostat=ios)symprec
        if(ios/=0.or.symprec<=0)call fail('invalid --symprec value')
      case('--sigma');call option_value(i,narg,'--sigma',value);read(value,*,iostat=ios)sigma
        if(ios/=0.or.sigma<0)call fail('invalid --sigma value')
      case('--bands');call option_value(i,narg,'--bands',value);read(value,*,iostat=ios)nbands
        if(ios/=0.or.nbands<1)call fail('invalid --bands value')
      case('--kpoints-file');call option_value(i,narg,'--kpoints-file',kpoints_path)
      case('--reference-eigenval');call option_value(i,narg,'--reference-eigenval',reference_path)
      case('--no-kpoint-symmetry');full_mesh=.true.
      case('--xc');call option_value(i,narg,'--xc',xc);xc=lower(trim(xc))
        if(trim(xc)/='lda'.and.trim(xc)/='pbe')call fail('--xc must be lda or pbe')
      case('--backend');call option_value(i,narg,'--backend',backend);backend=lower(trim(backend))
        if(trim(backend)/='auto'.and.trim(backend)/='cpu'.and.trim(backend)/='cuda')call fail('--backend must be auto, cpu, or cuda')
      case('--solver');call option_value(i,narg,'--solver',solver);solver=lower(trim(solver))
        if(trim(solver)/='evd'.and.trim(solver)/='evj')call fail('--solver must be evd or evj')
      case('--no-uspp-dij');use_uspp=.false.
      case('--forces');do_forces=.true.
      case('--force-step');call option_value(i,narg,'--force-step',value);read(value,*,iostat=ios)force_step
        if(ios/=0.or.force_step<=0)call fail('invalid --force-step value')
      case('--paw-atomic-double-counting');call option_value(i,narg,'--paw-atomic-double-counting',value);read(value,*,iostat=ios)paw_atomic
        if(ios/=0)call fail('invalid --paw-atomic-double-counting value');have_paw_atomic=.true.
      case('--atomic-reference-energy');call option_value(i,narg,'--atomic-reference-energy',value);read(value,*,iostat=ios)atomic_reference
        if(ios/=0)call fail('invalid --atomic-reference-energy value');have_atomic_override=.true.
      case('--output');call option_value(i,narg,'--output',output_path)
      case('--vaspwave-h5');call option_value(i,narg,'--vaspwave-h5',hdf5_path)
      case('--help','-h');call print_energy_help(6);return
      case default;call fail('unknown energy option: '//trim(arg))
      end select
      i=i+1
    end do
    call read_chgcar(trim(charge_path),crystal,rho);call read_potcar(trim(potential_path),potcars);call validate_potcar_structure(potcars,crystal)
    nelect=0.0_dp
    do it=1,size(potcars);nelect=nelect+potcars(it)%zval*real(crystal%counts(it),dp);end do
    density_electrons=rho%electron_count()
    if(abs(density_electrons-nelect)>5e-4_dp)call fail('smooth CHGCAR electron count does not match POTCAR ZVAL')
    minimum_bands=ceiling(nelect/2.0_dp);if(nbands==0)nbands=minimum_bands+8
    if(nbands<minimum_bands)call fail('too few bands for POTCAR electron count')
    if(len_trim(kpoints_path)>0)then;call read_explicit_kpoints(trim(kpoints_path),set)
    else if(full_mesh)then;call gamma_centered_mesh(crystal,kspacing,set)
    else;call gamma_centered_irreducible_mesh(crystal,kspacing,set,symprec,.true.);end if
    use_reference=len_trim(reference_path)>0
    if(use_reference)then
      call read_vasp_eigenval(trim(reference_path),reference_set,reference_eigenvalues,reference_nelect)
      if(reference_set%nk/=set%nk.or.size(reference_eigenvalues,2)<nbands)call fail('EIGENVAL dimensions do not match requested k points/bands')
      if(maxval(abs(reference_set%points-set%points))>2e-8_dp.or.maxval(abs(reference_set%weights-set%weights))>2e-8_dp) &
        call fail('EIGENVAL coordinates or weights do not match k points')
      if(abs(reference_nelect-nelect)>5e-4_dp)call fail('EIGENVAL electron count does not match POTCAR')
      if(len_trim(hdf5_path)>0)call fail('vaspwave.h5 export requires reconstructed eigenvectors, not EIGENVAL')
      if(do_forces)call fail('forces require reconstructed eigenstates, not EIGENVAL')
    end if
#ifdef HALF_CLI_CUDA
    if(trim(backend)=='auto'.or.trim(backend)=='cuda')then;actual_backend='cuda';else;actual_backend='cpu';end if
#else
    if(trim(backend)=='cuda')call fail('this HALF build has no CUDA backend');actual_backend='cpu'
#endif
#ifndef HALF_CLI_HAVE_MKL
    if(trim(actual_backend)=='cpu')call fail('the CPU energy backend requires oneMKL')
#endif
    if(trim(actual_backend)=='cpu'.and.trim(solver)/='evd')call fail('--solver evj is available only with CUDA')
    if(len_trim(hdf5_path)>0.and.trim(solver)/='evd')call fail('vaspwave.h5 export requires --solver evd')
#ifndef HALF_CLI_HAVE_HDF5
    if(len_trim(hdf5_path)>0)call fail('this HALF build has no HDF5 support')
#endif
    allocate(eigenvalues(set%nk,nbands),plane_waves(set%nk));eh=0;exc=0;exv=0
    if(use_reference)eigenvalues=reference_eigenvalues(:,:nbands)
#ifdef HALF_CLI_HAVE_HDF5
    if(len_trim(hdf5_path)>0)allocate(wave_bases(set%nk),waves(set%nk))
#endif
    do ik=1,set%nk
      call build_plane_wave_basis(crystal,rho%shape,encut,set%points(ik,:),basis);plane_waves(ik)=basis%npw
      if(nbands>basis%npw)call fail('requested bands exceed plane waves at a k point')
      if(use_reference.and.ik>1)cycle
#ifdef HALF_CLI_CUDA
      if(trim(actual_backend)=='cuda')then
        if(ik==1)then
          if(len_trim(hdf5_path)>0)then
            call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax, &
              potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp,eh,exc,exv,vectors)
          else
            call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax, &
              potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp,eh,exc,exv)
          end if
        else
          if(len_trim(hdf5_path)>0)then
            call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax, &
              potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp,eigenvectors=vectors)
          else
            call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax, &
              potential_seconds,assembly_seconds,gpu_seconds,solver,use_uspp)
          end if
        end if
      else
#endif
#ifdef HALF_CLI_HAVE_MKL
        call build_paw_operators(potcars,crystal,basis,paw)
        if(trim(xc)=='pbe')then;call build_veff_pbe(rho,potcars,crystal,veff,eh,exc,exv)
        else;call build_veff_lda(rho,potcars,crystal,veff,eh,exc,exv);end if
        if(use_uspp)call build_uspp_dij_cpu(veff,rho%shape,potcars,crystal,paw)
        if(len_trim(hdf5_path)>0)then;call solve_dense_gamma(veff,basis,paw,values,smin,smax,vectors)
        else;call solve_dense_gamma(veff,basis,paw,values,smin,smax);end if
#endif
#ifdef HALF_CLI_CUDA
      end if
#endif
      if(.not.use_reference)eigenvalues(ik,:)=values(:nbands)
#ifdef HALF_CLI_HAVE_HDF5
      if(len_trim(hdf5_path)>0)then;wave_bases(ik)=basis;waves(ik)%coefficients=vectors(:,:nbands);end if
#endif
      write(0,'(A,I0,A,I0,A,ES14.6)')'HALF energy: k point ',ik,'/',set%nk,' E1_eV=',values(1)
    end do
    call compute_occupations(eigenvalues,set%weights,nelect,sigma,occupation,mu,band_energy,entropy_term)
    if(sigma>0.and.maxval(occupation(:,nbands))>1e-6_dp)call fail('highest computed band is occupied; increase --bands')
#ifdef HALF_CLI_HAVE_HDF5
    if(len_trim(hdf5_path)>0)call write_vaspwave_h5(trim(hdf5_path),crystal,rho,encut,set%points,eigenvalues,occupation,mu,wave_bases,waves)
#endif
    allocate(charges(crystal%nions));ion=0
    do it=1,size(potcars);do iat=1,crystal%counts(it);ion=ion+1;charges(ion)=potcars(it)%zval;end do;end do
    call ewald_energy(crystal,charges,ewald,ewald_real,ewald_recip,ewald_self,ewald_background)
    if(.not.have_atomic_override)then;atomic_reference=0;do it=1,size(potcars)
      atomic_reference=atomic_reference+potcars(it)%eatom*real(crystal%counts(it),dp);end do;end if
    local_g0=0;do it=1,size(potcars)
      if(allocated(potcars(it)%psp_local))local_g0=local_g0+potcars(it)%psp_local(1)*real(crystal%counts(it),dp)
    end do;local_g0=local_g0*nelect/crystal%volume
    internal_energy=band_energy-eh-exv+exc+ewald+local_g0+atomic_reference
    if(have_paw_atomic)internal_energy=internal_energy+paw_atomic
    free_energy=internal_energy+entropy_term
    if(do_forces)then
      allocate(forces(crystal%nions,3));inverse_lattice=inverse3(crystal%lattice)
      do iat=1,crystal%nions;do i=1,3;do ion=1,2
        displaced=crystal;delta_cart=0.0_dp;delta_cart(i)=merge(-force_step,force_step,ion==1)
        displaced%positions(iat,:)=modulo(displaced%positions(iat,:)+matmul(delta_cart,inverse_lattice),1.0_dp)
        call evaluate_free_energy(displaced,rho,potcars,encut,kspacing,kpoints_path,full_mesh,symprec,nbands,sigma,xc, &
          actual_backend,solver,use_uspp,atomic_reference,paw_atomic,have_paw_atomic,displaced_energies(ion))
      end do
      forces(iat,i)=-(displaced_energies(2)-displaced_energies(1))/(2.0_dp*force_step)
      write(0,'(A,I0,A,I0,A,ES16.8)')'HALF force: atom ',iat,' axis ',i,' force_eV_per_A=',forces(iat,i)
      end do;end do
    else
      allocate(forces(0,0))
    end if
    if(trim(output_path)=='-')then;unit=6
    else;open(newunit=unit,file=trim(output_path),status='replace',action='write',iostat=ios)
      if(ios/=0)call fail('cannot write output: '//trim(output_path));end if
    call write_energy_json(unit,xc,actual_backend,solver,set,nelect,density_electrons,mu,band_energy,entropy_term, &
      eh,exv,exc,ewald,ewald_real,ewald_recip,ewald_self,ewald_background,local_g0,atomic_reference,paw_atomic, &
      have_paw_atomic,internal_energy,free_energy,smin,plane_waves,use_uspp,forces,do_forces,force_step)
    if(unit/=6)then;close(unit);call write_energy_json(6,xc,actual_backend,solver,set,nelect,density_electrons,mu,band_energy,entropy_term, &
      eh,exv,exc,ewald,ewald_real,ewald_recip,ewald_self,ewald_background,local_g0,atomic_reference,paw_atomic, &
      have_paw_atomic,internal_energy,free_energy,smin,plane_waves,use_uspp,forces,do_forces,force_step);end if
#endif
  end subroutine

  subroutine print_energy_help(unit)
    integer,intent(in)::unit
    write(unit,'(A)')'Usage: half energy CHARGE POTENTIAL [OPTIONS]'
    write(unit,'(A)')'       half-energy CHARGE POTENTIAL [OPTIONS]'
    write(unit,'(A)')'Options: --encut EV --kspacing VALUE --kpoints-file FILE --no-kpoint-symmetry'
    write(unit,'(A)')'         --symprec VALUE --bands N --sigma EV --xc lda|pbe'
    write(unit,'(A)')'         --reference-eigenval EIGENVAL'
    write(unit,'(A)')'         --backend auto|cpu|cuda --solver evd|evj --no-uspp-dij --output FILE'
    write(unit,'(A)')'         --vaspwave-h5 FILE'
    write(unit,'(A)')'         --forces --force-step ANGSTROM'
    write(unit,'(A)')'         --atomic-reference-energy EV --paw-atomic-double-counting EV'
  end subroutine

  subroutine evaluate_free_energy(crystal,rho,potcars,encut,kspacing,kpoints_path,full_mesh,symprec,nbands,sigma,xc, &
      backend,solver,use_uspp,atomic_reference,paw_atomic,have_paw_atomic,free_energy)
    type(crystal_t),intent(in)::crystal
    type(charge_grid_t),intent(in)::rho
    type(potcar_t),intent(in)::potcars(:)
    real(dp),intent(in)::encut,kspacing,symprec,sigma,atomic_reference,paw_atomic
    integer,intent(in)::nbands
    character(len=*),intent(in)::kpoints_path,xc,backend,solver
    logical,intent(in)::full_mesh,use_uspp,have_paw_atomic
    real(dp),intent(out)::free_energy
    type(kpoint_set_t)::mesh
    type(plane_wave_basis_t)::basis
    type(paw_species_t),allocatable::paw(:)
    real(dp),allocatable::veff(:),values(:),eigenvalues(:,:),occupation(:,:),charges(:)
    real(dp)::eh,exc,exv,smin,smax,ps,as,gs,nelect,mu,band,entropy,ewald,local_g0,internal
    integer::ik,it,iat,ion
    if(len_trim(kpoints_path)>0)then;call read_explicit_kpoints(trim(kpoints_path),mesh)
    else if(full_mesh)then;call gamma_centered_mesh(crystal,kspacing,mesh)
    else;call gamma_centered_irreducible_mesh(crystal,kspacing,mesh,symprec,.true.);end if
    allocate(eigenvalues(mesh%nk,nbands));eh=0;exc=0;exv=0
    do ik=1,mesh%nk
      call build_plane_wave_basis(crystal,rho%shape,encut,mesh%points(ik,:),basis)
      if(nbands>basis%npw)call fail('requested bands exceed plane waves in displaced force evaluation')
#ifdef HALF_CLI_CUDA
      if(trim(backend)=='cuda')then
        if(ik==1)then
          call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax,ps,as,gs, &
            solver,use_uspp,eh,exc,exv)
        else
          call solve_dense_gamma_cuda_full(rho,potcars,crystal,basis,trim(xc)=='pbe',values,smin,smax,ps,as,gs,solver,use_uspp)
        end if
      else
#endif
#ifdef HALF_CLI_HAVE_MKL
        call build_paw_operators(potcars,crystal,basis,paw)
        if(trim(xc)=='pbe')then;call build_veff_pbe(rho,potcars,crystal,veff,eh,exc,exv)
        else;call build_veff_lda(rho,potcars,crystal,veff,eh,exc,exv);end if
        if(use_uspp)call build_uspp_dij_cpu(veff,rho%shape,potcars,crystal,paw)
        call solve_dense_gamma(veff,basis,paw,values,smin,smax)
#endif
#ifdef HALF_CLI_CUDA
      end if
#endif
      eigenvalues(ik,:)=values(:nbands)
    end do
    nelect=0;do it=1,size(potcars);nelect=nelect+potcars(it)%zval*real(crystal%counts(it),dp);end do
    call compute_occupations(eigenvalues,mesh%weights,nelect,sigma,occupation,mu,band,entropy)
    if(sigma>0.and.maxval(occupation(:,nbands))>1e-6_dp)call fail('highest force-evaluation band is occupied; increase --bands')
    allocate(charges(crystal%nions));ion=0
    do it=1,size(potcars);do iat=1,crystal%counts(it);ion=ion+1;charges(ion)=potcars(it)%zval;end do;end do
    call ewald_energy(crystal,charges,ewald)
    local_g0=0;do it=1,size(potcars)
      if(allocated(potcars(it)%psp_local))local_g0=local_g0+potcars(it)%psp_local(1)*real(crystal%counts(it),dp)
    end do;local_g0=local_g0*nelect/crystal%volume
    internal=band-eh-exv+exc+ewald+local_g0+atomic_reference
    if(have_paw_atomic)internal=internal+paw_atomic
    free_energy=internal+entropy
  end subroutine

  subroutine write_energy_json(unit,xc,backend,solver,set,nelect,density_electrons,mu,band,entropy,eh,exv,exc,ewald, &
      er,eg,es,eb,local_g0,atomic_reference,paw_atomic,have_paw_atomic,internal,free,smin,plane_waves,use_uspp, &
      forces,have_forces,force_step)
    integer,intent(in)::unit
    character(len=*),intent(in)::xc,backend,solver
    type(kpoint_set_t),intent(in)::set
    real(dp),intent(in)::nelect,density_electrons,mu,band,entropy,eh,exv,exc,ewald,er,eg,es,eb,local_g0
    real(dp),intent(in)::atomic_reference,paw_atomic,internal,free,smin
    integer(i64),intent(in)::plane_waves(:)
    real(dp),intent(in)::forces(:,:),force_step
    logical,intent(in)::have_paw_atomic,use_uspp,have_forces
    integer::iat
    write(unit,'(A)')'{';write(unit,'(A)')'  "implementation": "DeePAW-HALF 0.4.0",'
    write(unit,'(A,A,A)')'  "xc": "',trim(xc),'",';write(unit,'(A,A,A)')'  "backend": "',trim(backend),'",'
    write(unit,'(A,A,A)')'  "solver": "',trim(solver),'",';write(unit,'(A,A,A)')'  "uspp_dij": ',merge('true ','false',use_uspp),','
    write(unit,'(A,I0,A)')'  "full_kpoint_count": ',set%full_count,',';write(unit,'(A,I0,A)')'  "irreducible_kpoint_count": ',set%nk,','
    write(unit,'(A,I0,A,I0,A)')'  "plane_wave_range": [',minval(plane_waves),', ',maxval(plane_waves),'],'
    write(unit,'(A,ES24.16,A)')'  "electron_count": ',nelect,',';write(unit,'(A,ES24.16,A)')'  "smooth_density_electron_count": ',density_electrons,','
    write(unit,'(A,ES24.16,A)')'  "fermi_level_eV": ',mu,',';write(unit,'(A,ES24.16,A)')'  "band_energy_eV": ',band,','
    write(unit,'(A,ES24.16,A)')'  "entropy_minus_Ts_eV": ',entropy,',';write(unit,'(A,ES24.16,A)')'  "hartree_double_counting_eV": ',-eh,','
    write(unit,'(A,ES24.16,A)')'  "xc_potential_double_counting_eV": ',-exv,',';write(unit,'(A,ES24.16,A)')'  "xc_energy_eV": ',exc,','
    write(unit,'(A,ES24.16,A)')'  "ewald_energy_eV": ',ewald,',';write(unit,'(A,4(ES24.16,:,A))')'  "ewald_components_eV": [',er,', ',eg,', ',es,', ',eb,'],'
    write(unit,'(A,ES24.16,A)')'  "local_pseudopotential_g0_correction_eV": ',local_g0,','
    write(unit,'(A,ES24.16,A)')'  "atomic_reference_energy_eV": ',atomic_reference,','
    if(have_paw_atomic)then;write(unit,'(A,ES24.16,A)')'  "paw_atomic_double_counting_eV": ',paw_atomic,','
    else;write(unit,'(A)')'  "paw_atomic_double_counting_eV": null,';end if
    write(unit,'(A,ES24.16,A)')'  "minimum_overlap_eigenvalue": ',smin,','
    write(unit,'(A,ES24.16,A)')'  "internal_energy_eV": ',internal,',';write(unit,'(A,ES24.16,A)')'  "free_energy_eV": ',free,','
    if(have_forces)then
      write(unit,'(A,ES24.16,A)')'  "force_step_Angstrom": ',force_step,',';write(unit,'(A)')'  "forces_eV_per_Angstrom": ['
      do iat=1,size(forces,1)
        write(unit,'(A,3(ES24.16,:,A))',advance='no')'    [',forces(iat,1),', ',forces(iat,2),', ',forces(iat,3),']'
        if(iat<size(forces,1))then;write(unit,'(A)')',';else;write(unit,'(A)')'';end if
      end do
      write(unit,'(A)')'  ]'
    else
      write(unit,'(A)')'  "forces_eV_per_Angstrom": null'
    end if
    write(unit,'(A)')'}'
  end subroutine

  subroutine option_value(index,narg,name,value)
    integer,intent(inout)::index
    integer,intent(in)::narg
    character(len=*),intent(in)::name
    character(len=*),intent(out)::value
    if(index>=narg)call fail(trim(name)//' requires a value')
    index=index+1; call get_command_argument(index,value)
  end subroutine

  subroutine unavailable(name,reason)
    character(len=*),intent(in)::name,reason
    call fail(trim(name)//' is not available in HALF 0.4.0: '//trim(reason))
  end subroutine

  subroutine fail(message)
    character(len=*),intent(in)::message
    write(0,'(A)') 'half: error: '//trim(message)
    error stop 2
  end subroutine

  pure function lower(text)result(output)
    character(len=*),intent(in)::text
    character(len=len(text))::output
    integer::i,c
    output=text
    do i=1,len(text)
      c=iachar(output(i:i)); if(c>=iachar('A').and.c<=iachar('Z'))output(i:i)=achar(c+32)
    end do
  end function

  pure function basename(path)result(name)
    character(len=*),intent(in)::path
    character(len=len(path))::name
    integer::i
    i=scan(path,'/',back=.true.)
    if(i>0)then; name=path(i+1:); else; name=path; end if
  end function
end program half_cli
