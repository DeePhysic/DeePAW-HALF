module half_local_forces
  use half_kinds,only:dp
  use half_constants,only:pi,felect
  use half_types,only:crystal_t,charge_grid_t,potcar_t,plane_wave_basis_t
  use half_fft,only:fft3_forward,fft3_backward
  use half_potential,only:build_veff_lda,build_veff_pbe
  use half_potcar_interp,only:vasp_local_spline_second,vasp_local_spline_eval,vasp_four_point_eval
  implicit none
  private
  public::local_ionic_forces,nlcc_forces,accumulate_smooth_density,harris_correction_forces
contains
  subroutine local_ionic_forces(charge,potcars,crystal,forces,energy)
    ! VASP FORLOC in a full complex FFT representation.  The input density
    ! is transformed once; differentiation acts only on exp[-i G.R_I].
    type(charge_grid_t),intent(in)::charge
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    real(dp),intent(out)::forces(:,:)
    real(dp),intent(out),optional::energy
    complex(dp),allocatable,target::rho_r(:),rho_g(:)
    real(dp),allocatable::charge_c(:),m2(:)
    real(dp)::q(3),g2,gn,radial,phase,rpos(3),term,total_energy
    integer::i1,i2,i3,n1,n2,n3,idx,it,iat,ion0,n
    if(size(forces,1)/=crystal%nions.or.size(forces,2)/=3) &
      error stop 'HALF: local force dimensions do not match the crystal'
    n=size(charge%values);allocate(charge_c(n),rho_r(n),rho_g(n));call reorder(charge%values,charge%shape,charge_c)
    rho_r=cmplx(charge_c,0.0_dp,dp);call fft3_forward(charge%shape,rho_r,rho_g);rho_g=rho_g/real(n,dp)
    forces=0.0_dp;total_energy=0.0_dp
    do it=1,size(potcars)
      allocate(m2(size(potcars(it)%psp_local)));call vasp_local_spline_second(potcars(it)%psp_local,potcars(it)%psp_gmax,m2)
      ion0=sum(crystal%counts(:it-1))
      do iat=1,crystal%counts(it)
        rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice);idx=0
        do i1=0,charge%shape(1)-1
          n1=fft_integer(i1,charge%shape(1))
          do i2=0,charge%shape(2)-1
            n2=fft_integer(i2,charge%shape(2))
            do i3=0,charge%shape(3)-1
              n3=fft_integer(i3,charge%shape(3));idx=idx+1
              q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal);g2=dot_product(q,q)
              if(g2<=1.0e-12_dp)cycle
              gn=sqrt(g2);if(gn>potcars(it)%psp_gmax-potcars(it)%psp_gmax/real(size(m2),dp))cycle
              radial=vasp_local_spline_eval(potcars(it)%psp_local,m2,potcars(it)%psp_gmax,gn)- &
                4*pi*potcars(it)%zval*felect/g2
              phase=dot_product(q,rpos)
              term=radial/crystal%volume*aimag(conjg(rho_g(idx))*cmplx(cos(phase),-sin(phase),dp))
              forces(ion0+iat,:)=forces(ion0+iat,:)-q*term
              total_energy=total_energy+radial/crystal%volume*real(conjg(rho_g(idx))* &
                cmplx(cos(phase),-sin(phase),dp),dp)
            end do
          end do
        end do
      end do
      deallocate(m2)
    end do
    if(present(energy))energy=total_energy
  end subroutine local_ionic_forces

  subroutine nlcc_forces(vxc,potcars,crystal,shape,forces)
    ! VASP FORCOR/FORHAR(LPAR=.TRUE.) in a full complex FFT layout:
    ! F_I = - integral v_xc(r) d n_core,I(r-R_I)/dR_I dr.
    ! CHGCAR/core arrays store Omega*n, so the integral is a grid average.
    real(dp),intent(in)::vxc(:)
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    integer,intent(in)::shape(3)
    real(dp),intent(out)::forces(:,:)
    complex(dp),allocatable,target::work(:),vxc_g(:)
    real(dp)::q(3),g2,gn,radial,phase,rpos(3),term
    integer::n,idx,i1,i2,i3,n1,n2,n3,it,iat,ion0
    n=size(vxc)
    if(n/=product(shape).or.size(forces,1)/=crystal%nions.or.size(forces,2)/=3) &
      error stop 'HALF: NLCC force dimensions do not match the FFT grid/crystal'
    allocate(work(n),vxc_g(n));work=cmplx(vxc,0.0_dp,dp);call fft3_forward(shape,work,vxc_g)
    vxc_g=vxc_g/real(n,dp);forces=0.0_dp
    do it=1,size(potcars)
      ion0=sum(crystal%counts(:it-1));if(.not.potcars(it)%has_core)cycle
      do iat=1,crystal%counts(it)
        rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice);idx=0
        do i1=0,shape(1)-1
          n1=fft_integer(i1,shape(1))
          do i2=0,shape(2)-1
            n2=fft_integer(i2,shape(2))
            do i3=0,shape(3)-1
              n3=fft_integer(i3,shape(3));idx=idx+1
              q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal);g2=dot_product(q,q)
              if(g2<=1.0e-12_dp)cycle
              gn=sqrt(g2);if(gn>potcars(it)%psp_gmax-3*potcars(it)%psp_gmax/real(size(potcars(it)%pspcor),dp))cycle
              radial=vasp_four_point_eval(potcars(it)%pspcor,potcars(it)%psp_gmax,gn)
              phase=dot_product(q,rpos)
              term=radial*aimag(conjg(vxc_g(idx))*cmplx(cos(phase),-sin(phase),dp))
              forces(ion0+iat,:)=forces(ion0+iat,:)-q*term
            end do
          end do
        end do
      end do
    end do
  end subroutine nlcc_forces

  subroutine accumulate_smooth_density(basis,eigenvectors,occupations,kweight,shape,density_f)
    type(plane_wave_basis_t),intent(in)::basis
    complex(dp),intent(in)::eigenvectors(:,:)
    real(dp),intent(in)::occupations(:),kweight
    integer,intent(in)::shape(3)
    real(dp),intent(inout)::density_f(:)
    complex(dp),allocatable,target::grid_g(:),grid_r(:)
    real(dp),allocatable::density_c(:)
    integer::n,ib,ig,i1,i2,i3,source,destination
    n=product(shape)
    if(size(density_f)/=n.or.size(eigenvectors,1)/=basis%npw.or.size(eigenvectors,2)<size(occupations)) &
      error stop 'HALF: wave-density dimensions are inconsistent'
    allocate(grid_g(n),grid_r(n),density_c(n));density_c=0.0_dp
    do ib=1,size(occupations)
      if(occupations(ib)==0.0_dp)cycle
      grid_g=(0.0_dp,0.0_dp)
      do ig=1,basis%npw;grid_g(basis%fft_index(ig))=eigenvectors(ig,ib);end do
      call fft3_backward(shape,grid_g,grid_r)
      density_c=density_c+kweight*occupations(ib)*abs(grid_r)**2
    end do
    do i1=0,shape(1)-1;do i2=0,shape(2)-1;do i3=0,shape(3)-1
      source=i1*shape(2)*shape(3)+i2*shape(3)+i3+1
      destination=i1+shape(1)*(i2+shape(2)*i3)+1
      density_f(destination)=density_f(destination)+density_c(source)
    end do;end do;end do
  end subroutine accumulate_smooth_density

  subroutine harris_correction_forces(input_density,output_density,potcars,crystal,use_pbe,forces)
    ! Frozen-input-density Harris response, equivalent to VASP's
    ! CHGGRA -> FORHAR(LPAR=.FALSE.) path:
    !
    !   g_Hxc = [v_H[n_out-n_in] + f_xc[n_in+n_core](n_out-n_in)]
    !   F_I^Harris = - integral g_Hxc(r) d n_atom,I(r-R_I)/dR_I dr.
    !
    ! The full effective-potential central difference supplies both the
    ! linear Hartree response and the XC kernel.  FORHAR couples this field
    ! to PSPRHO (the atomic valence reference density), not PSPCOR; PSPCOR
    ! belongs to the separate NLCC/FORCOR force above.
    type(charge_grid_t),intent(in)::input_density,output_density
    type(potcar_t),intent(in)::potcars(:)
    type(crystal_t),intent(in)::crystal
    logical,intent(in)::use_pbe
    real(dp),intent(out)::forces(:,:)
    type(charge_grid_t)::plus_density,minus_density
    real(dp),allocatable::delta_f(:),veff_plus(:),veff_minus(:),vxc(:),response(:)
    real(dp)::eh,exc,exv,t
    integer::n
    n=size(input_density%values)
    if(size(output_density%values)/=n.or.any(input_density%shape/=output_density%shape)) &
      error stop 'HALF: Harris input/output density grids differ'
    allocate(delta_f(n),response(n))
    delta_f=output_density%values-input_density%values;t=1.0e-4_dp
    plus_density%shape=input_density%shape;minus_density%shape=input_density%shape
    allocate(plus_density%values(n),minus_density%values(n))
    plus_density%values=input_density%values+t*delta_f;minus_density%values=input_density%values-t*delta_f
    if(use_pbe)then
      call build_veff_pbe(plus_density,potcars,crystal,veff_plus,eh,exc,exv,vxc)
      call build_veff_pbe(minus_density,potcars,crystal,veff_minus,eh,exc,exv,vxc)
    else
      call build_veff_lda(plus_density,potcars,crystal,veff_plus,eh,exc,exv,vxc)
      call build_veff_lda(minus_density,potcars,crystal,veff_minus,eh,exc,exv,vxc)
    end if
    response=(veff_plus-veff_minus)/(2*t)
    call atomic_radial_density_forces(response,potcars,crystal,input_density%shape,.false.,forces)
  end subroutine harris_correction_forces

  subroutine atomic_radial_density_forces(potential,potcars,crystal,shape,use_core,forces)
    real(dp),intent(in)::potential(:);type(potcar_t),intent(in)::potcars(:);type(crystal_t),intent(in)::crystal
    integer,intent(in)::shape(3);logical,intent(in)::use_core;real(dp),intent(out)::forces(:,:)
    complex(dp),allocatable,target::work(:),potential_g(:)
    real(dp),allocatable::radial_values(:)
    real(dp)::q(3),g2,gn,radial,phase,rpos(3),term;integer::n,idx,i1,i2,i3,n1,n2,n3,it,iat,ion0
    n=size(potential);allocate(work(n),potential_g(n));work=cmplx(potential,0.0_dp,dp)
    call fft3_forward(shape,work,potential_g);potential_g=potential_g/real(n,dp);forces=0.0_dp
    do it=1,size(potcars)
      ion0=sum(crystal%counts(:it-1))
      if(use_core)then
        if(.not.potcars(it)%has_core)cycle;radial_values=potcars(it)%pspcor
      else
        if(.not.allocated(potcars(it)%psprho))cycle;radial_values=potcars(it)%psprho
      end if
      do iat=1,crystal%counts(it);rpos=matmul(crystal%positions(ion0+iat,:),crystal%lattice);idx=0
        do i1=0,shape(1)-1;n1=fft_integer(i1,shape(1));do i2=0,shape(2)-1;n2=fft_integer(i2,shape(2))
          do i3=0,shape(3)-1;n3=fft_integer(i3,shape(3));idx=idx+1
            q=matmul([real(n1,dp),real(n2,dp),real(n3,dp)],crystal%reciprocal);g2=dot_product(q,q);if(g2<=1e-12_dp)cycle
            gn=sqrt(g2);if(gn>potcars(it)%psp_gmax-3*potcars(it)%psp_gmax/real(size(radial_values),dp))cycle
            radial=vasp_four_point_eval(radial_values,potcars(it)%psp_gmax,gn);phase=dot_product(q,rpos)
            term=radial*aimag(conjg(potential_g(idx))*cmplx(cos(phase),-sin(phase),dp))
            forces(ion0+iat,:)=forces(ion0+iat,:)-q*term
          end do;end do
        end do
      end do;deallocate(radial_values)
    end do
  end subroutine atomic_radial_density_forces

  subroutine reorder(input,shape,output)
    real(dp),intent(in)::input(:);integer,intent(in)::shape(3);real(dp),intent(out)::output(:)
    integer::i1,i2,i3,source,destination
    do i1=0,shape(1)-1;do i2=0,shape(2)-1;do i3=0,shape(3)-1
      source=i1+shape(1)*(i2+shape(2)*i3)+1
      destination=i1*shape(2)*shape(3)+i2*shape(3)+i3+1;output(destination)=input(source)
    end do;end do;end do
  end subroutine reorder
  subroutine reorder_c_to_f(input,shape,output)
    real(dp),intent(in)::input(:);integer,intent(in)::shape(3);real(dp),intent(out)::output(:)
    integer::i1,i2,i3,source,destination
    do i1=0,shape(1)-1;do i2=0,shape(2)-1;do i3=0,shape(3)-1
      source=i1*shape(2)*shape(3)+i2*shape(3)+i3+1
      destination=i1+shape(1)*(i2+shape(2)*i3)+1;output(destination)=input(source)
    end do;end do;end do
  end subroutine reorder_c_to_f
  pure integer function fft_integer(index0,n)result(value)
    integer,intent(in)::index0,n
    if(index0<=(n-1)/2)then;value=index0;else;value=index0-n;end if
  end function fft_integer
end module half_local_forces
