program half_cli
  use half_kinds, only: dp, i64
  use half_types, only: crystal_t, charge_grid_t, plane_wave_basis_t, potcar_t
  use half_chgcar, only: read_chgcar
  use half_potcar, only: read_potcar
  use half_basis, only: build_plane_wave_basis
  use half_paw, only: paw_species_t, build_paw_operators
  use half_uspp, only: build_uspp_dij_cpu
  use half_cpu, only: cpu_density_mean
#ifdef HALF_CLI_HAVE_MKL
  use half_potential, only: build_veff_lda, build_veff_pbe
  use half_dense_solver, only: solve_dense_gamma
#endif
#ifdef HALF_CLI_CUDA
  use half_cuda_solver, only: solve_dense_gamma_cuda_full
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
    write(unit,'(A)') '  bands    reserved for HAPPY-compatible band reconstruction'
    write(unit,'(A)') '  energy   reserved for HAPPY-compatible total-energy evaluation'
    write(unit,'(A)') '  version  print the HALF version'
    write(unit,'(A)') ''
    write(unit,'(A)') 'Run half COMMAND --help for command-specific options.'
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
    call read_chgcar(trim(chg),crystal,rho); call read_potcar(trim(pot),p)
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
    character(len=1024)::arg
    if(command_argument_count()>=offset+1)then
      call get_command_argument(offset+1,arg)
      if(trim(arg)=='--help'.or.trim(arg)=='-h')then
        write(*,'(A)')'Usage: half bands CHARGE POTENTIAL [OPTIONS]'
        write(*,'(A)')'       half-bands CHARGE POTENTIAL [OPTIONS]'
        write(*,'(A)')''
        write(*,'(A)')'Planned HAPPY-compatible options: --output-prefix, --encut, --bands,'
        write(*,'(A)')'--npoints, --path, --xc, and --vaspwave-h5.'
        write(*,'(A)')''
        write(*,'(A)')'Status: arbitrary-k Hamiltonians are not available in HALF 0.4.0.'
        return
      end if
    end if
    call unavailable('bands','arbitrary-k Hamiltonians are still being ported')
  end subroutine

  subroutine command_energy(offset)
    integer,intent(in)::offset
    character(len=1024)::arg
    if(command_argument_count()>=offset+1)then
      call get_command_argument(offset+1,arg)
      if(trim(arg)=='--help'.or.trim(arg)=='-h')then
        write(*,'(A)')'Usage: half energy CHARGE POTENTIAL [OPTIONS]'
        write(*,'(A)')'       half-energy CHARGE POTENTIAL [OPTIONS]'
        write(*,'(A)')''
        write(*,'(A)')'Planned HAPPY-compatible options: --encut, --kspacing, --kpoints-file,'
        write(*,'(A)')'--bands, --sigma, --xc, --forces, --vaspwave-h5, and --output-prefix.'
        write(*,'(A)')''
        write(*,'(A)')'Status: total-energy and force terms are not available in HALF 0.4.0.'
        return
      end if
    end if
    call unavailable('energy','total-energy and force terms are still being ported')
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
