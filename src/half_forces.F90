module half_forces
  use half_kinds,only:dp
  use half_types,only:plane_wave_basis_t
  use half_paw,only:paw_species_t
  implicit none
  private
  public::add_nonlocal_paw_forces
contains
  subroutine add_nonlocal_paw_forces(basis,paw,eigenvalues,occupations,kweight,eigenvectors,forces)
    ! Analytic non-local PAW/USPP force for a generalized eigenproblem:
    !
    !   F_I = -sum_n f_n w_k <psi_n|dH_I-eps_n dS_I|psi_n>
    !       = -2 sum_n f_n w_k Re[c_n^H (D_I-eps_n Q) dc_n/dR_I].
    !
    ! Projectors have the HALF/VASP phase exp[-i(G+k).R_I], hence
    ! dc_a/dR_alpha = i sum_G q_G,alpha conj(P_Ga) psi_G.
    type(plane_wave_basis_t),intent(in)::basis
    type(paw_species_t),intent(in)::paw(:)
    real(dp),intent(in)::eigenvalues(:),occupations(:),kweight
    complex(dp),intent(in)::eigenvectors(:,:)
    real(dp),intent(inout)::forces(:,:)
    complex(dp),allocatable::c(:),dc(:,:),metric_c(:)
    real(dp),allocatable::metric(:,:)
    integer::it,iat,ion0,ib,a,b,ig,alpha,nlm
    if(size(eigenvectors,1)/=basis%npw.or.size(eigenvectors,2)<size(eigenvalues).or. &
       size(occupations)/=size(eigenvalues).or.size(forces,2)/=3) &
      error stop 'HALF: non-local force dimensions are inconsistent'
    ion0=0
    do it=1,size(paw)
      nlm=paw(it)%nlm
      allocate(c(nlm),dc(nlm,3),metric_c(nlm),metric(nlm,nlm))
      do iat=1,paw(it)%natoms
        if(ion0+iat>size(forces,1))error stop 'HALF: non-local force atom index overflow'
        do ib=1,size(eigenvalues)
          if(occupations(ib)==0.0_dp)cycle
          c=(0.0_dp,0.0_dp);dc=(0.0_dp,0.0_dp)
          do a=1,nlm
            do ig=1,basis%npw
              c(a)=c(a)+conjg(paw(it)%projectors(iat,a,ig))*eigenvectors(ig,ib)
              do alpha=1,3
                dc(a,alpha)=dc(a,alpha)+cmplx(0.0_dp,basis%qvectors(ig,alpha),dp)* &
                  conjg(paw(it)%projectors(iat,a,ig))*eigenvectors(ig,ib)
              end do
            end do
          end do
          metric=paw(it)%dij_atom(iat,:,:)-eigenvalues(ib)*paw(it)%qij
          metric_c=matmul(metric,c)
          do alpha=1,3
            forces(ion0+iat,alpha)=forces(ion0+iat,alpha)-2.0_dp*kweight*occupations(ib)* &
              real(dot_product(metric_c,dc(:,alpha)),dp)
          end do
        end do
      end do
      deallocate(c,dc,metric_c,metric);ion0=ion0+paw(it)%natoms
    end do
  end subroutine add_nonlocal_paw_forces
end module half_forces
