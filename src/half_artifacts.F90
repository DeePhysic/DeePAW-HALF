module half_artifacts
  use iso_c_binding,only:c_int,c_double,c_int64_t,c_char,c_null_char
  use half_kinds,only:dp,i64
  implicit none
  private
  public::write_bands_artifacts,write_energy_npz
  interface
    integer(c_int) function bands_npz(prefix,nk,nb,kpoints,x,eig,shifted,npw,smin,vbm)bind(C,name='half_write_bands_npz_c')
      import::c_char,c_int,c_double,c_int64_t
      character(c_char),intent(in)::prefix(*)
      integer(c_int),value::nk,nb
      real(c_double),intent(in)::kpoints(*),x(*),eig(*),shifted(*),smin(*);integer(c_int64_t),intent(in)::npw(*);real(c_double),value::vbm
    end function
    integer(c_int) function bands_png(prefix,nk,nb,x,shifted)bind(C,name='half_write_bands_png_c')
      import::c_char,c_int,c_double
      character(c_char),intent(in)::prefix(*);integer(c_int),value::nk,nb;real(c_double),intent(in)::x(*),shifted(*)
    end function
    integer(c_int) function energy_npz(prefix,nk,nb,kpoints,weights,eig,occ,nf,forces)bind(C,name='half_write_energy_npz_c')
      import::c_char,c_int,c_double
      character(c_char),intent(in)::prefix(*);integer(c_int),value::nk,nb,nf
      real(c_double),intent(in)::kpoints(*),weights(*),eig(*),occ(*),forces(*)
    end function
  end interface
contains
  subroutine write_bands_artifacts(prefix,kpoints,x,eig,shifted,npw,smin,vbm,status)
    character(len=*),intent(in)::prefix
    real(dp),intent(in)::kpoints(:,:),x(:),eig(:,:),shifted(:,:),smin(:),vbm
    integer(i64),intent(in)::npw(:);integer,intent(out)::status
    character(kind=c_char,len=:),allocatable::c
    c=trim(prefix)//c_null_char;status=bands_npz(c,size(eig,1),size(eig,2),kpoints,x,eig,shifted,npw,smin,vbm)
    if(status==0)status=bands_png(c,size(eig,1),size(eig,2),x,shifted)
  end subroutine
  subroutine write_energy_npz(prefix,kpoints,weights,eig,occ,forces,status)
    character(len=*),intent(in)::prefix
    real(dp),intent(in)::kpoints(:,:),weights(:),eig(:,:),occ(:,:),forces(:,:);integer,intent(out)::status
    character(kind=c_char,len=:),allocatable::c
    c=trim(prefix)//c_null_char;status=energy_npz(c,size(eig,1),size(eig,2),kpoints,weights,eig,occ,size(forces,1),forces)
  end subroutine
end module half_artifacts
