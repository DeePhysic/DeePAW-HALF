module half_dense_solver
  use iso_c_binding,only:c_int,c_char,c_double,c_double_complex
  use half_kinds,only:dp
  use half_types,only:plane_wave_basis_t
  use half_paw,only:paw_species_t
  use half_fft,only:fft3_forward
  implicit none
  private
  public::solve_dense_gamma,assemble_dense_gamma
  interface
    function lapacke_zheev(layout,jobz,uplo,n,a,lda,w)bind(C,name='LAPACKE_zheev')result(info)
      import c_int,c_char,c_double,c_double_complex
      integer(c_int),value::layout,n,lda
      character(c_char),value::jobz,uplo
      complex(c_double_complex)::a(*)
      real(c_double)::w(*)
      integer(c_int)::info
    end function
    function lapacke_zhegv(layout,itype,jobz,uplo,n,a,lda,b,ldb,w)bind(C,name='LAPACKE_zhegv')result(info)
      import c_int,c_char,c_double,c_double_complex
      integer(c_int),value::layout,itype,n,lda,ldb
      character(c_char),value::jobz,uplo
      complex(c_double_complex)::a(*),b(*)
      real(c_double)::w(*)
      integer(c_int)::info
    end function
  end interface
contains
  subroutine solve_dense_gamma(veff,basis,paw,eigenvalues,overlap_min,overlap_max)
    real(dp),intent(in)::veff(:)
    type(plane_wave_basis_t),intent(in)::basis
    type(paw_species_t),intent(in)::paw(:)
    real(dp),allocatable,intent(out)::eigenvalues(:)
    real(dp),intent(out)::overlap_min,overlap_max
    complex(dp),allocatable::h(:,:),s(:,:),scopy(:,:)
    integer::n,info
    n=basis%npw; allocate(scopy(n,n),eigenvalues(n))
    call assemble_dense_gamma(veff,basis,paw,h,s)
    scopy=s
    info=lapacke_zheev(102,'N','U',n,scopy,n,eigenvalues)
    if(info/=0)error stop 'HALF: overlap eigensolve failed'
    overlap_min=eigenvalues(1); overlap_max=eigenvalues(n)
    info=lapacke_zhegv(102,1,'N','U',n,h,n,s,n,eigenvalues)
    if(info/=0)error stop 'HALF: generalized eigensolve failed'
  end subroutine

  subroutine assemble_dense_gamma(veff,basis,paw,h,s)
    real(dp),intent(in)::veff(:)
    type(plane_wave_basis_t),intent(in)::basis
    type(paw_species_t),intent(in)::paw(:)
    complex(dp),allocatable,intent(out)::h(:,:),s(:,:)
    complex(dp),allocatable,target::vin(:),vg(:)
    integer::i,j,it,iat,a,b,idx,n
    complex(dp)::z
    n=basis%npw; allocate(vin(size(veff)),vg(size(veff)),h(n,n),s(n,n))
    vin=cmplx(veff,0.0_dp,dp); call fft3_forward(basis%shape,vin,vg); vg=vg/real(size(veff),dp)
    do j=1,n
      do i=1,n
        idx=fft_index_of_difference(basis%gvectors(i,:),basis%gvectors(j,:),basis%shape)
        h(i,j)=vg(idx); s(i,j)=(0.0_dp,0.0_dp)
      end do
      h(j,j)=h(j,j)+basis%kinetic(j); s(j,j)=1
    end do
    do it=1,size(paw)
      do iat=1,paw(it)%natoms
        do j=1,n
          do i=1,n
            z=(0.0_dp,0.0_dp)
            do a=1,paw(it)%nlm; do b=1,paw(it)%nlm
              z=z+paw(it)%projectors(iat,a,i)*paw(it)%dij_atom(iat,a,b)*conjg(paw(it)%projectors(iat,b,j))
            end do; end do
            h(i,j)=h(i,j)+z; z=(0.0_dp,0.0_dp)
            do a=1,paw(it)%nlm; do b=1,paw(it)%nlm
              z=z+paw(it)%projectors(iat,a,i)*paw(it)%qij(a,b)*conjg(paw(it)%projectors(iat,b,j))
            end do; end do
            s(i,j)=s(i,j)+z
          end do
        end do
      end do
    end do
    h=0.5_dp*(h+conjg(transpose(h))); s=0.5_dp*(s+conjg(transpose(s)))
  end subroutine
  pure integer function fft_index_of_difference(gi,gj,shape)result(index1)
    integer,intent(in)::gi(3),gj(3),shape(3)
    integer::a,b,c
    a=modulo(gi(1)-gj(1),shape(1)); b=modulo(gi(2)-gj(2),shape(2)); c=modulo(gi(3)-gj(3),shape(3))
    index1=a*shape(2)*shape(3)+b*shape(3)+c+1
  end function
end module half_dense_solver
