#include "symbol.inc"

!> Module to deal with scattering events and their mapping
module triplets
  use prec, only: q
  use lattice, only: latt
  use mkpoints_struct_def, only: kpoints_struct
  implicit none

  !> Contains all possible internal representations of one k-point
  type kpoint_type
    integer :: kindex !< integer index as computed by the hash table
    integer :: ibz !< integer index in the list of k-points
    integer :: kint(3) !< three intergers defining the k-point
    !TODO: integer :: g(3) !< integer part of the k-point
  end type kpoint_type

  !> Contains all the representations of a scattering event between 3 states
  type triplet_type
    type(kpoint_type) :: k
    type(kpoint_type) :: kp
    type(kpoint_type) :: q
  end type triplet_type
  contains

  !> initialize the elph_kpoint structure from kindex
  pure subroutine kpoint_init_kindex(self, kpoints, kindex)
    use mkpoints, only: kindex_vkpt, kindex_kint, kint_vkpt, ikibz_kindex
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    integer,intent(in) :: kindex
    !PROFILING_START('kpoint_init_kindex')
    self%kindex = kindex
    self%kint = kindex_kint(kpoints, self%kindex)
    self%ibz = kpoints%fbz2ibz(1, self%kindex)
    !PROFILING_STOP('kpoint_init_kindex')
  end subroutine kpoint_init_kindex

  !> initialize the elph_kpoint structure from kibz
  pure subroutine kpoint_init_kibz(self, kpoints, latt_cur, kibz)
    use mkpoints, only: kindex_vkpt, kindex_kint, kint_vkpt, kint_kindex, kpoint_integers
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: kibz
    !PROFILING_START('kpoint_init_kibz')
    self%ibz = kibz
    call kpoint_integers(kpoints,latt_cur,kpoints%vkpt(:,self%ibz),self%kint)
    self%kindex = kint_kindex(kpoints,self%kint)
    !PROFILING_STOP('kpoint_init_kibz')
  end subroutine kpoint_init_kibz

  !> initialize the elph_kpoint structure from vkpt
  pure subroutine kpoint_init_kint(self, kpoints, kint)
    use mkpoints, only: kpoint_integers, kindex_kint, kint_kindex, ikibz_kindex, kint_mod
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    integer,intent(in) :: kint(3)
    !PROFILING_START('kpoint_init_kint')
    self%kint = kint_mod(kpoints,kint)
    self%kindex = kint_kindex(kpoints,self%kint)
    self%ibz = kpoints%fbz2ibz(1,self%kindex)
    !PROFILING_STOP('kpoint_init_kint')
  end subroutine kpoint_init_kint

  !> initialize the elph_kpoint structure from vkpt
  pure subroutine kpoint_init_vkpt(self, kpoints, latt_cur, vkpt)
    use mkpoints, only: kpoint_integers, kindex_kint, kint_kindex, ikibz_kindex
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    real(q),intent(in) :: vkpt(3)
    !PROFILING_START('kpoint_init_vkpt')
    call kpoint_integers(kpoints,latt_cur,vkpt,self%kint)
    self%kindex = kint_kindex(kpoints,self%kint)
    self%ibz = kpoints%fbz2ibz(1,self%kindex)
    !PROFILING_STOP('kpoint_init_vkpt')
  end subroutine kpoint_init_vkpt

  !> intialize the triplet_type from two integer indexes
  subroutine triplet_init_k_and_kp(self, kpoints, latt_cur, ik_fbz, ikp_fbz)
    use mkpoints, only: kindex_from_k_and_kp
    type(triplet_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: ik_fbz !< integer index in the hash table of kp
    integer,intent(in) :: ikp_fbz !< integer index in the hash table of iq
    ! local variables
    integer :: iq_fbz
    PROFILING_START('triplet_init_k_and_kp')
    call kpoint_init_kindex(self%k, kpoints, ik_fbz)
    call kpoint_init_kindex(self%kp, kpoints, ikp_fbz)
    iq_fbz = kindex_from_k_and_kp(kpoints, latt_cur, ik_fbz, ikp_fbz)
    call kpoint_init_kindex(self%q, kpoints, iq_fbz)
    PROFILING_STOP('triplet_init_k_and_kp')
  end subroutine triplet_init_k_and_kp

  !> intialize the triplet_type from two integer indexes
  subroutine triplet_init_kp_and_q(self, kpoints, latt_cur, ikp_fbz, iq_fbz)
    use mkpoints, only: kindex_from_kp_and_q
    type(triplet_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: ikp_fbz !< integer index in the hash table of kp
    integer,intent(in) :: iq_fbz !< integer index in the hash table of iq
    ! local variables
    integer :: ik_fbz
    PROFILING_START('triplet_init_kp_and_q')
    call kpoint_init_kindex(self%kp, kpoints, ikp_fbz)
    call kpoint_init_kindex(self%q, kpoints, iq_fbz)
    ik_fbz = kindex_from_kp_and_q(kpoints, latt_cur, ikp_fbz, iq_fbz)
    call kpoint_init_kindex(self%k, kpoints, ik_fbz)
    PROFILING_STOP('triplet_init_kp_and_q')
  end subroutine triplet_init_kp_and_q

end module triplets
