#include "symbol.inc"

module pot_electrostat

__GPU use moffload_struct_def
__GPU use moffload
      use prec
      use mpimy
      use mgrid
      use lattice
      use constant
      use pot_struct_def
      use coulomb_cutoff
      use pseudo_struct_def, only: NPSPTS, potcar
      use poscar_struct_def, only: type_info
      use charge, only: rhoadd

      implicit none

      private

      public :: pothar, potion, potion_particle_mesh
      public :: potential_from_ion_charge, rescaling_factor

      ! needed for the unit tests only
      public :: electrostatic_energy, pscenc_energy

contains

!************************ SUBROUTINE POTHAR ****************************
!
!> This subroutine is just a wrapper that passes off to
!> POTHAR_ (legacy, in pot.F) or POTHAR_TRUNCATED (below)
!> The former is the default, the latter is called if the
!> coulomb kernel is truncated.
!
!***********************************************************************

    subroutine pothar(gridc, latt_cur, chtot, cvd, denc, coulomb_pot)

      type(grid_3d) :: gridc
      type(latt) :: latt_cur
      complex(q) :: chtot(*)
      complex(q) :: cvd(*)
      real(q) :: denc

      type(coulomb_potential), optional :: coulomb_pot

      PROFILING_START('pothar')

      if (present(coulomb_pot)) then
         if (coulomb_pot%kernel_truncate%active) then
            call pothar_truncated(gridc, latt_cur, coulomb_pot, chtot, cvd, denc)
         else
            call pothar_(gridc, latt_cur, chtot, cvd, denc)
         endif
      else
         call pothar_(gridc, latt_cur, chtot, cvd, denc)
      endif

      PROFILING_STOP('pothar')

      return
    end subroutine pothar


!************************ SUBROUTINE POTION ****************************
!
!> This subroutine is just a wrapper that passes off to
!> POTION_ (legacy, in pot.F) or POTION_TRUNCATED (below)
!> The former is the default, the latter is called if the
!> coulomb kernel is truncated.
!
!***********************************************************************

    subroutine potion(gridc, p, latt_cur, t_info, cvps, cdvps, cstrf, pscenc, coulomb_pot)

      type(grid_3d) :: gridc
      type(latt) :: latt_cur
      type(potcar) :: p(:)
      type(type_info) :: t_info
      complex(q) :: cvps(*), cdvps(:)
      complex(q) :: cstrf(:,:)
      real(q) :: pscenc

      type(coulomb_potential), optional :: coulomb_pot

      PROFILING_START('potion')

      if (present(coulomb_pot)) then
         if (coulomb_pot%kernel_truncate%active) then
            call potion_truncated(gridc, p, latt_cur, t_info, coulomb_pot, cvps, cstrf, pscenc)
            cdvps = 0
         else
            call potion_(gridc, p, latt_cur, t_info, cvps, cdvps, cstrf, pscenc)
         endif
      else
         call potion_(gridc, p, latt_cur, t_info, cvps, cdvps, cstrf, pscenc)
      endif

      PROFILING_STOP('potion')

      return
    end subroutine potion


!************************ SUBROUTINE POTION_PARTICLE_MESH **************
!
!> This subroutine is just a wrapper that passes off to
!> POTION_PARTICLE_MESH_ (legacy) or POTION_PARTICLE_MESH_TRUNCATED
!> The former is the default, the latter is called if the
!> coulomb kernel is truncated.
!
!***********************************************************************

    subroutine potion_particle_mesh(gridc, p, latt_cur, t_info, cvps, pscenc, tewen, for, coulomb_pot)

      type (grid_3d) gridc
      type (type_info) t_info
      type (potcar) p(:)
      type (latt) latt_cur
      complex(q) cvps(*)
      real(q) pscenc,tewen
      real(q), optional :: for(:,:)

      type(coulomb_potential), optional :: coulomb_pot

      PROFILING_START('potion_particle_mesh')

      if (present(coulomb_pot)) then
         if (coulomb_pot%kernel_truncate%active) then
            call potion_particle_mesh_truncated(gridc, p, latt_cur, t_info, coulomb_pot, cvps, pscenc, tewen, for)
         else
            call potion_particle_mesh_(gridc, p, latt_cur, t_info, cvps, pscenc, tewen, for)
         endif
      else
         call potion_particle_mesh_(gridc, p, latt_cur, t_info, cvps, pscenc, tewen, for)
      endif

      PROFILING_STOP('potion_particle_mesh')

      return
    end subroutine potion_particle_mesh


    subroutine pothar_truncated(gridc, latt_cur, coulomb_pot, charge_density, &
                                potential, denc)
        type(grid_3d), intent(in)            :: gridc
        type(latt), intent(in)               :: latt_cur
        type(coulomb_potential), intent(in)  :: coulomb_pot
        complex(q)                           :: charge_density(gridc%rc%np)
        complex(q), intent(out)              :: potential(gridc%rc%np)
        real(q), intent(out)                 :: denc

        PROFILING_START('pothar_truncated')
        !
        if (coulomb_pot%kernel_truncate%active) then
            !
            call cutoff_hartree_potential(gridc, latt_cur, &
                                          coulomb_pot, &
                                          charge_density, potential)
        else
            !
            call potential_from_charge_density(gridc, latt_cur, &
                                               coulomb_pot, &
                                               charge_density, &
                                               potential)
            !
            call setunb(potential, gridc)
        end if
        !
        call electrostatic_energy(gridc, charge_density, potential, &
                                  denc)
        !
        PROFILING_STOP('pothar_truncated')
    end subroutine

    subroutine potion_truncated(gridc, p, latt_cur, t_info, coulomb_pot, &
                                ion_potential, structure_factor, pscenc, &
                                over_ions)
      type (grid_3d), intent(in)           :: gridc
      type (type_info), intent(in)         :: t_info
      type (potcar), intent(in)            :: p (t_info%ntyp)
      type (latt), intent(in)              :: latt_cur
      type (coulomb_potential), intent(in) :: coulomb_pot
      complex(q), intent(in)               :: structure_factor(gridc%mplwv,t_info%ntyp)
      complex(q), intent(out)              :: ion_potential(gridc%rc%np)
      real(q), intent(out)                 :: pscenc       ! non-electrostatic contributions to energy
      integer, optional                    :: over_ions(:) ! if less than the total number of ions is required
      PROFILING_START('potion_truncated')

      call pscenc_energy(p,t_info,latt_cur,coulomb_pot,pscenc)
      !
GPACC !$acc kernels present(ion_potential) __IF_ASYNC__
GPOMP !$omp target map(present,to: ion_potential) if(offload_on)
      ion_potential = 0._q
GPOMP !$omp end target
GPACC !$acc end kernels
      !
      if (coulomb_pot%kernel_truncate%active) then
         !
         call cutoff_ion_potential(p, t_info, gridc, latt_cur, &
                                   coulomb_pot, structure_factor, &
                                   ion_potential, &
                                   use_grad_kernels=.false., &
                                   over_ions=over_ions)
      else
         !
         call potential_from_ion_charge(p, t_info, gridc, latt_cur, &
                                        coulomb_pot, structure_factor, &
                                        ion_potential, &
                                        use_grad_kernels=.false., &
                                        over_ions=over_ions)
         !
         call setunb(ion_potential, gridc)
      end if

      PROFILING_STOP('potion_truncated')
    end subroutine

    subroutine grad_potion(gridc, p, latt_cur, t_info, coulomb_pot, &
                           grad_ion_potential, structure_factor)
      type (grid_3d), intent(in)           :: gridc
      type (type_info), intent(in)         :: t_info
      type (potcar), intent(in)            :: p (t_info%ntyp)
      type (latt), intent(in)              :: latt_cur
      type (coulomb_potential), intent(in) :: coulomb_pot
      complex(q), intent(in)               :: structure_factor(gridc%mplwv,t_info%ntyp)
      complex(q), intent(out)              :: grad_ion_potential(gridc%rc%np)

      !
      !$acc kernels present(grad_ion_potential) __IF_ASYNC__
      grad_ion_potential = 0._q
      !$acc end kernels
      !
      if (coulomb_pot%kernel_truncate%active) then
         !
         call cutoff_ion_potential(p, t_info, gridc, latt_cur, &
                                   coulomb_pot, structure_factor, &
                                   grad_ion_potential, &
                                   use_grad_kernels=.true.)
      else
         !
         call potential_from_ion_charge(p, t_info, gridc, latt_cur, &
                                   coulomb_pot, structure_factor, &
                                   grad_ion_potential, &
                                   use_grad_kernels=.true.)
         !
         call setunb(grad_ion_potential, gridc)
      end if
      !
      end subroutine

      subroutine potential_from_charge_density(gridc, latt_cur, &
                                               coulomb_pot, &
                                               charge_density, &
                                               potential)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER potential_from_charge_density
#endif
          type (grid_3d), intent(in)           :: gridc
          type (latt), intent(in)              :: latt_cur
          type(coulomb_potential), intent(in)  :: coulomb_pot
          complex(q), intent(in)               :: charge_density(gridc%rc%np)
          complex(q), intent(inout)            :: potential(gridc%rc%np)
          integer                              :: idx(3)
          real(q)                              :: g_vector(3)
          real(q)                              :: kernel
          integer                              :: n1, n2, n3, ni, nc

GPACC !$acc routine (index_to_gvector) seq
GPACC !$acc routine (kernel_for_gvector) seq
GPACC !$acc routine (potential_from_kernel_for_gvector) seq
          !
GPACC !$acc parallel loop collapse(2) private(n2,n3,ni,idx,g_vector,kernel) &
GPACC !$acc present(gridc,latt_cur,coulomb_pot,charge_density,potential) __IF_ASYNC__
          !
GPOMP !$omp target teams distribute parallel do collapse(2) private(n2,n3,ni,idx,g_vector,kernel) &
GPOMP !$omp map(present,to: gridc,latt_cur,coulomb_pot,charge_density,potential) if(offload_on)
          !
          !
          do nc=1,gridc%rc%ncol
__CPU     n2= gridc%rc%i2(nc)
__CPU     n3= gridc%rc%i3(nc)
          do n1=1,gridc%rc%nrow
            !
__GPU       n2= gridc%rc%i2(nc)
__GPU       n3= gridc%rc%i3(nc)
            !
            ni = (nc - 1) * gridc%rc%nrow + n1
            !
            idx = [gridc%lpctx(n1), gridc%lpcty(n2), gridc%lpctz(n3)]
            !
            call index_to_gvector(latt_cur, idx, g_vector)
            !
            call kernel_for_gvector(coulomb_pot, &
                                    g_vector, kernel, all(idx==0))
            !
            call potential_from_kernel_for_gvector(kernel, &
                                                   charge_density(ni), &
                                                   latt_cur%omega, &
                                                   potential(ni))
          enddo
          enddo
      end subroutine

      subroutine potential_from_kernel_for_gvector(kernel, &
                                                        charge_density,&
                                                        volume, &
                                                        potential)
GPACC !$acc routine seq
GPOMP !$omp declare target(potential_from_kernel_for_gvector)
         real(q), intent(in)     :: kernel
         complex(q), intent(in)  :: charge_density
         real(q), intent(in)     :: volume
         complex(q), intent(out) :: potential

         potential = kernel * charge_density / volume
      end subroutine

      subroutine electrostatic_energy(gridc, charge_density, potential,&
                                      energy)
         type(grid_3d), intent(in) :: gridc
         complex(q), intent(in)    :: charge_density(gridc%rc%np)
         complex(q), intent(in)    :: potential(gridc%rc%np)
         real(q), intent(out)      :: energy
         integer :: nc, n1, n2, n3, ni, ncol, nrow
         real(q) :: factm, dum

         energy = 0._q
         ! cast loop bounds in local variables for faster GPU execution
         ncol = gridc%rc%ncol
         nrow = gridc%rc%nrow
         !
GPACC    !$acc parallel loop collapse(2) reduction(+:energy) private(n2,n3,ni,dum,factm) &
GPACC    !$acc present(gridc,potential,charge_density) __IF_ASYNC__
         !
GPOMP    !$omp target teams distribute parallel do collapse(2) reduction(+:energy) &
GPOMP    !$omp private(n2,n3,ni,dum,factm) map(present,to:gridc,potential,charge_density) if(offload_on)
         !
         do nc = 1, ncol
__CPU    n2 = gridc%rc%i2(nc)
__CPU    n3 = gridc%rc%i3(nc)
         do n1 = 1, nrow
__GPU      n2 = gridc%rc%i2(nc)
__GPU      n3 = gridc%rc%i3(nc)
           !
           ni = (nc - 1) * nrow + n1
           SETFACT1
           SETFACT
           !
           dum = MULFACT potential(ni) * conjg(charge_density(ni))
           energy = energy + dum
         enddo
         enddo
         !
         energy = -energy / 2._q
         CALLMPI( M_sum_d(gridc%comm, energy, 1))
      end subroutine

      subroutine iteration_information(iterating_over_ions, &
                                       iterate_over_array, &
                                       total_number_types, &
                                       over_ions)
         logical, intent(out)              :: iterating_over_ions
         integer, allocatable, intent(out) :: iterate_over_array(:)
         integer, intent(in)               :: total_number_types
         integer, optional, intent(in)     :: over_ions(:)
         integer                           :: nt

         if (present(over_ions)) then
            !
            iterating_over_ions = .true.
         else
            !
            iterating_over_ions = .false.
         end if
         !
         if (iterating_over_ions) then
             !
             allocate(iterate_over_array, source=over_ions)
         else
             !
             allocate(iterate_over_array(total_number_types))
             !
             do nt = 1, total_number_types
                !
                iterate_over_array(nt) = nt
             end do
         end if
      end subroutine

      subroutine potential_from_ion_charge(p, t_info, gridc, latt_cur, &
                                           coulomb_pot, &
                                           structure_factor, &
                                           potential, &
                                           use_grad_kernels, &
                                           over_ions)
         type(grid_3d), intent(in)           :: gridc
         type(type_info), intent(in)         :: t_info
         type(potcar), intent(in)            :: p (t_info%ntyp)
         type (latt), intent(in)             :: latt_cur
         type(coulomb_potential), intent(in) :: coulomb_pot
         complex(q), intent(in)              :: structure_factor(gridc%mplwv,t_info%ntyp)
         complex(q), intent(out)             :: potential(gridc%rc%np)
         logical, intent(in)                 :: use_grad_kernels
         integer, optional                   :: over_ions(:)
         integer                             :: nc, ni, n1, n2, n3
         integer                             :: j, idx_type
         real(q)                             :: g_vector(3)
         real(q)                             :: g_norm
         real(q)                             :: psgma2, zz
         real(q)                             :: kernel, ps_kernel
         real(q)                             :: summed_kernel
         integer, allocatable                :: iterate_over(:)
         integer                             :: idx(3)
         logical                             :: iterating_over_ions
         complex(q)                          :: charge
         complex(q)                          :: potential_
         logical                             :: is_g0
GPACC !$acc routine (index_to_gvector) seq
GPACC !$acc routine (kernel_for_gvector) seq
GPACC !$acc routine (kernel_gradient_for_gvector) seq
GPACC !$acc routine (pseudo_kernel_for_gvector) seq
GPACC !$acc routine (pseudo_kernel_gradient_for_gvector) seq
GPACC !$acc routine (phase_factor_for_index) seq
GPACC !$acc routine (potential_from_kernel_for_gvector) seq

GPACC    !$acc enter data copyin(t_info) __IF_ASYNC__
GPACC    !$acc enter data copyin(t_info%vca,t_info%posion) __IF_ASYNC__
         !
GPOMP    !$omp target enter data map(to: t_info) if(offload_on)
GPOMP    !$omp target enter data map(to: t_info%vca,t_info%posion) if(offload_on)
IFXW2 GPOMP !$omp target enter data map(always,alloc: t_info%vca,t_info%posion) if(offload_on)
         !
GPACC    !$acc kernels present(potential) __IF_ASYNC__
GPOMP    !$omp target map(present,to: potential) if(offload_on)
IFXW8 GPOMP !$omp teams workdistribute if(offload_on)
         potential = 0._q
IFXW8 GPOMP !$omp end teams workdistribute
GPOMP    !$omp end target
GPACC    !$acc end kernels
         !
         call iteration_information(iterating_over_ions, iterate_over, &
                                    t_info%ntyp, over_ions)
         !
         do j = 1, size(iterate_over)
           !
           if (iterating_over_ions) then
               !
               idx_type = t_info%ityp(iterate_over(j))
           else
               !
               idx_type = iterate_over(j)
           end if
           !
           psgma2 = p(idx_type)%psgmax - p(idx_type)%psgmax / npspts
           !
           zz = -p(idx_type)%zvalf
           !
           !
GPACC !$acc parallel loop collapse(2) private(n2,n3,ni,idx,g_vector,g_norm,is_g0,charge,kernel,ps_kernel,summed_kernel,potential_) &
GPACC !$acc& present(gridc,p,t_info,latt_cur,coulomb_pot,structure_factor,potential) __IF_ASYNC__
           !
! for OpenMP present checks are performed with MAP(PRESENT,TO:...), so for conditional data,
! MAP(ALLOC:...) is used to avoid potentially copying the data for the parallel loop.
GPOMP !$omp target teams distribute parallel do collapse(2) &
GPOMP !$omp private(n2,n3,ni,idx,g_vector,g_norm,is_g0,charge,kernel,ps_kernel,summed_kernel,potential_) &
GPOMP !$omp map(present,to: gridc,p,t_info,latt_cur,coulomb_pot,structure_factor,potential) if(offload_on)
           !
           !
           do nc = 1, gridc%rc%ncol
__CPU      n2 = gridc%rc%i2(nc)
__CPU      n3 = gridc%rc%i3(nc)
           do n1 = 1, gridc%rc%nrow
__GPU        n2 = gridc%rc%i2(nc)
__GPU        n3 = gridc%rc%i3(nc)
             !
             ni = (nc - 1) * gridc%rc%nrow + n1
             !
             idx = [gridc%lpctx(n1), gridc%lpcty(n2), gridc%lpctz(n3)]
             !
             is_g0 = all(idx==0)
             !
             call index_to_gvector(latt_cur, idx, g_vector)
             !
             g_norm = sqrt(  g_vector(1)**2._q &
                           + g_vector(2)**2._q &
                           + g_vector(3)**2._q )
             !
             if ( g_norm < psgma2 )  then
                 !
                 if (use_grad_kernels) then
                    !
                     call kernel_gradient_for_gvector(&
                                     coulomb_pot, &
                                     g_vector, kernel, is_g0)
                     !
                     call pseudo_kernel_gradient_for_gvector(&
                                     p, t_info, coulomb_pot, &
                                     g_vector, idx_type, ps_kernel, &
                                     is_g0)
                else
                     !
                     call kernel_for_gvector(&
                                     coulomb_pot, &
                                     g_vector, kernel, is_g0)
                     !
                     call pseudo_kernel_for_gvector(&
                                     p, t_info, coulomb_pot, &
                                     g_vector, idx_type, ps_kernel, &
                                     is_g0)
                end if
                !
                if (iterating_over_ions) then
                   !
                   call phase_factor_for_index(idx, &
                                   t_info%posion(:,iterate_over(j)), &
                                   t_info%vca(idx_type), &
                                   charge)
                else
                  !
                   charge = structure_factor(ni, idx_type)
                end if
                !
                summed_kernel = ps_kernel + zz * kernel
                !
                call potential_from_kernel_for_gvector(summed_kernel, &
                                                       charge, &
                                                       latt_cur%omega,&
                                                       potential_)
                potential(ni) = potential(ni) + potential_
             else
                potential(ni) = 0._q
             endif
           enddo
           enddo
           enddo
GPACC      !$acc exit data delete(t_info) __IF_ASYNC__
GPACC      !$acc exit data delete(t_info%vca,t_info%posion) __IF_ASYNC__
           !
GPOMP      !$omp target exit data map(release: t_info%vca,t_info%posion) if(offload_on)
IFXW2 GPOMP !$omp target exit data map(release: t_info%vca,t_info%posion) if(offload_on)
GPOMP      !$omp target exit data map(release: t_info) if(offload_on)
           !
      end subroutine

!> calculate the contribution to the total energy from the non-coulomb
!> part of the g=0 component of the pseudopotential and the force on the
!> unit cell due to the change in this energy as the size of the cell
!> changes
      subroutine pscenc_energy(p, t_info, latt_cur, coulomb_pot, pscenc)
          type(type_info), intent(in)         :: t_info
          type(potcar), intent(in)            :: p (t_info%ntyp)
          type(latt), intent(in)              :: latt_cur
          type(coulomb_potential), intent(in) :: coulomb_pot
          real(q), intent(out)                :: pscenc
          real(q)                             :: zvsum
          integer :: nt

          if (coulomb_pot%kernel_truncate%active) then
              !
              pscenc = 0._q
          else
              !
              zvsum=0
              do nt=1,t_info%ntyp
                 !
                 zvsum=zvsum+p(nt)%zvalf*t_info%nityp(nt)*t_info%vca(nt)
              enddo
              !
              pscenc = 0._q
              do nt=1,t_info%ntyp
                 !
                 pscenc=pscenc + p(nt)%pscore * t_info%vca(nt) &
                                              * t_info%nityp(nt) &
                                              * (zvsum / latt_cur%omega)
              enddo
          end if
      end subroutine

     subroutine cutoff_hartree_potential(gridc, latt_cur, coulomb_pot, &
                                         charge_density, potential)
        type(grid_3d), intent(in)            :: gridc
        type(latt), intent(in)               :: latt_cur
        type(coulomb_potential), intent(in)  :: coulomb_pot
        complex(q)                           :: charge_density(gridc%rc%np)
        complex(q), intent(out)              :: potential(gridc%rc%np)
        type(latt)                           :: latt_cur_padded
        type(coulomb_potential)              :: periodic_coulomb_pot
        complex(q), allocatable              :: padded_charge_density(:)
        complex(q), allocatable              :: coarse_charge_density(:)
        complex(q), allocatable              :: padded_potential(:)
        complex(q), allocatable              :: coarse_open_potential(:)
        complex(q), allocatable              :: coarse_periodic_potential(:)
        complex(q), allocatable              :: periodic_potential(:)
        real(q)                              :: rescale

        call rescaling_factor(coulomb_pot%kernel_truncate, rescale)
        !
__GPU   call copyin_typed_var(periodic_coulomb_pot)
        !
        allocate(padded_charge_density(coulomb_pot%kernel_truncate%fine_grid%mplwv))
        allocate(padded_potential(coulomb_pot%kernel_truncate%fine_grid%mplwv))
        !
        padded_charge_density = 0._q
        padded_potential = 0._q
        !
GPACC   !$acc enter data copyin(padded_charge_density) __IF_ASYNC__
GPACC   !$acc enter data copyin(padded_potential) __IF_ASYNC__
        !
GPOMP   !$omp target enter data map(to: padded_charge_density) if(offload_on)
GPOMP   !$omp target enter data map(to: padded_potential) if(offload_on)
        !
        call generate_padded_lattice(coulomb_pot%kernel_truncate, &
                                     latt_cur, latt_cur_padded)
        !
GPACC   !$acc enter data copyin(latt_cur_padded) __IF_ASYNC__
GPOMP   !$omp target enter data map(to: latt_cur_padded) if(offload_on)
        !
        call setunb_compat(charge_density, gridc)
        !
        if (coulomb_pot%kernel_truncate%coarsen_before_pad) then
            !
GPACC       !$acc kernels present(potential) __IF_ASYNC__
GPOMP       !$omp target map(present,to: potential) if(offload_on)
IFXW8 GPOMP !$omp teams workdistribute if(offload_on)
            potential = 0._q
IFXW8 GPOMP !$omp end teams workdistribute
GPOMP       !$omp end target
GPACC       !$acc end kernels
            !
            allocate(coarse_charge_density(coulomb_pot%kernel_truncate%coarse_grid%mplwv))
            coarse_charge_density = 0._q
GPACC       !$acc enter data copyin(coarse_charge_density) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: coarse_charge_density) if(offload_on)
            !
            allocate(coarse_open_potential(coulomb_pot%kernel_truncate%coarse_grid%mplwv))
            coarse_open_potential = 0._q
GPACC       !$acc enter data copyin(coarse_open_potential) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: coarse_open_potential) if(offload_on)
            !
            allocate(coarse_periodic_potential(coulomb_pot%kernel_truncate%coarse_grid%mplwv))
            coarse_periodic_potential = 0._q
GPACC       !$acc enter data copyin(coarse_periodic_potential) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: coarse_periodic_potential) if(offload_on)
            !
            allocate(periodic_potential(gridc%mplwv))
            periodic_potential = 0._q
            !
GPACC       !$acc enter data copyin(periodic_potential) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: periodic_potential) if(offload_on)
            call cp_grid(coulomb_pot%kernel_truncate%fine_grid, &
                         coulomb_pot%kernel_truncate%coarse_grid, &
                         coulomb_pot%kernel_truncate%translation_table, &
                         charge_density, coarse_charge_density)
            !
            call setunb(coarse_charge_density, &
                        coulomb_pot%kernel_truncate%coarse_grid)
            !
            call pad3d(coulomb_pot%kernel_truncate, &
                       coarse_charge_density, &
                       padded_charge_density, 1, rescale)
            !
            call setunb_compat(padded_charge_density, &
                        coulomb_pot%kernel_truncate%fine_grid)
            !
            call potential_from_charge_density(coulomb_pot%kernel_truncate%fine_grid,&
                                latt_cur_padded,coulomb_pot,&
                                padded_charge_density, padded_potential)
            !
            call setunb_compat(padded_potential, coulomb_pot%kernel_truncate%fine_grid)
            !
            call pad3d(coulomb_pot%kernel_truncate, coarse_open_potential, &
                       padded_potential, -1)
            !
            call setunb_compat(coarse_open_potential, coulomb_pot%kernel_truncate%coarse_grid)
            !
            call potential_from_charge_density(coulomb_pot%kernel_truncate%coarse_grid, latt_cur, &
                                periodic_coulomb_pot, &
                                coarse_charge_density, &
                                coarse_periodic_potential)
            !
            call setunb(coarse_periodic_potential, coulomb_pot%kernel_truncate%coarse_grid)
            !
            call rc_add(coarse_open_potential, 1._q, &
                        coarse_periodic_potential, -1._q, &
                        coarse_open_potential, coulomb_pot%kernel_truncate%coarse_grid)
            !
            call cpb_grid(coulomb_pot%kernel_truncate%fine_grid, coulomb_pot%kernel_truncate%coarse_grid, coulomb_pot%kernel_truncate%translation_table, &
                          coarse_open_potential, potential)
            !
            call setunb_compat(potential, coulomb_pot%kernel_truncate%fine_grid)
            !
            call potential_from_charge_density(gridc, latt_cur, &
                                periodic_coulomb_pot, charge_density, &
                                periodic_potential)
            !
            call setunb_compat(periodic_potential, gridc)
            !
            call rc_add(potential, 1._q, &
                        periodic_potential, 1._q, &
                        potential, gridc)
            !
GPACC       !$acc exit data delete(coarse_charge_density) __IF_ASYNC__
GPACC       !$acc exit data delete(coarse_open_potential) __IF_ASYNC__
GPACC       !$acc exit data delete(coarse_periodic_potential) __IF_ASYNC__
GPACC       !$acc exit data delete(periodic_potential) __IF_ASYNC__
            !
GPOMP       !$omp target exit data map(release: coarse_charge_density) if(offload_on)
GPOMP       !$omp target exit data map(release: coarse_open_potential) if(offload_on)
GPOMP       !$omp target exit data map(release: coarse_periodic_potential) if(offload_on)
GPOMP       !$omp target exit data map(release: periodic_potential) if(offload_on)
            !
        else
            !
            call pad3d(coulomb_pot%kernel_truncate, charge_density, &
                       padded_charge_density, 1, rescale)
            !
            call setunb_compat(padded_charge_density, coulomb_pot%kernel_truncate%fine_grid)
            !
            call potential_from_charge_density(coulomb_pot%kernel_truncate%fine_grid, &
                                latt_cur_padded,coulomb_pot, &
                                padded_charge_density,padded_potential)
            !
            call setunb_compat(padded_potential, coulomb_pot%kernel_truncate%fine_grid)
            !
            call pad3d(coulomb_pot%kernel_truncate, potential, padded_potential, -1)
        end if
        !
GPACC   !$acc exit data delete(padded_charge_density) __IF_ASYNC__
GPACC   !$acc exit data delete(padded_potential) __IF_ASYNC__
GPACC   !$acc exit data delete(latt_cur_padded) __IF_ASYNC__
        !
GPOMP   !$omp target exit data map(release: padded_charge_density) if(offload_on)
GPOMP   !$omp target exit data map(release: padded_potential) if(offload_on)
GPOMP   !$omp target exit data map(release: latt_cur_padded) if(offload_on)
        !
        call setunb_compat(potential, gridc)
__GPU   call delete_typed_var(periodic_coulomb_pot)
        !
        return
     end subroutine

     subroutine cutoff_ion_potential(p, t_info, gridc, latt_cur, &
                                     coulomb_pot, structure_factor, &
                                     potential, use_grad_kernels, &
                                     over_ions)
        type (grid_3d), intent(in)           :: gridc
        type (type_info), intent(in)         :: t_info
        type (potcar), intent(in)            :: p (t_info%ntyp)
        type (latt), intent(in)              :: latt_cur
        type (coulomb_potential), intent(in) :: coulomb_pot
        complex(q), intent(in)               :: structure_factor(gridc%mplwv,t_info%ntyp)
        complex(q), intent(inout)            :: potential(gridc%rc%np)
        logical, intent(in)                  :: use_grad_kernels
        integer, optional                    :: over_ions(:)
        type(type_info)                      :: t_info_pad
        type(latt)                           :: latt_cur_pad
        complex(q), allocatable              :: structure_factor_pad(:,:)
        complex(q), allocatable              :: structure_factor_coarse(:,:)
        complex(q), allocatable              :: potential_pad(:)
        complex(q), allocatable              :: potential_coarse(:)
        complex(q), allocatable              :: potential_coarse_periodic(:)
        complex(q), allocatable              :: potential_periodic(:)
        type(coulomb_potential)              :: coulomb_pot_periodic
        real(q), allocatable                 :: padded_posion(:,:)
        integer :: nt

__GPU   call copyin_typed_var(coulomb_pot_periodic)
IFXW2 GPOMP !$omp target enter data map(always,alloc: coulomb_pot_periodic) if(offload_on)
        !
        allocate(structure_factor_pad(coulomb_pot%kernel_truncate%fine_grid%mplwv, t_info%ntyp))
        allocate(structure_factor_coarse(coulomb_pot%kernel_truncate%coarse_grid%mplwv, t_info%ntyp))
        allocate(potential_pad(coulomb_pot%kernel_truncate%fine_grid%mplwv))
        allocate(potential_coarse(coulomb_pot%kernel_truncate%coarse_grid%mplwv))
        allocate(potential_coarse_periodic(coulomb_pot%kernel_truncate%coarse_grid%mplwv))
        allocate(potential_periodic(gridc%mplwv))
        !
        do nt = 1, t_info%ntyp
            !
            call setunb_compat(structure_factor(:,nt), gridc)
        end do
        !
        call generate_padded_lattice(coulomb_pot%kernel_truncate, latt_cur, &
                                     latt_cur_pad)
GPACC   !$acc enter data copyin(latt_cur_pad) __IF_ASYNC__
GPOMP   !$omp target enter data map(to: latt_cur_pad) if(offload_on)
        !
        call generate_padded_t_info(coulomb_pot%kernel_truncate, t_info, &
                                    t_info_pad, padded_posion)
        !
        structure_factor_pad = 0._q
        !
        call stufak(coulomb_pot%kernel_truncate%fine_grid, t_info_pad, structure_factor_pad)
        !
GPACC   !$acc enter data copyin(structure_factor_pad) __IF_ASYNC__
GPOMP   !$omp target enter data map(to: structure_factor_pad) if(offload_on)
        !
        do nt = 1, t_info%ntyp
            !
            call setunb_compat(structure_factor_pad(:,nt), coulomb_pot%kernel_truncate%fine_grid)
        end do
        !
        potential_pad = 0._q
GPACC   !$acc enter data copyin(potential_pad) __IF_ASYNC__
GPOMP   !$omp target enter data map(to: potential_pad) if(offload_on)
        call potential_from_ion_charge(p, t_info_pad,&
                                       coulomb_pot%kernel_truncate%fine_grid,latt_cur_pad, &
                                       coulomb_pot, &
                                       structure_factor_pad, &
                                       potential_pad, &
                                       use_grad_kernels=.false., &
                                       over_ions=over_ions)
        !
        call setunb(potential_pad, coulomb_pot%kernel_truncate%fine_grid)
        !
        if (coulomb_pot%kernel_truncate%coarsen_before_pad) then
            !
            potential_coarse = 0._q
GPACC       !$acc enter data copyin(potential_coarse) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: potential_coarse) if(offload_on)
            call pad3d(coulomb_pot%kernel_truncate, potential_coarse, &
                       potential_pad, -1)
            !
            call setunb_compat(potential_coarse, coulomb_pot%kernel_truncate%coarse_grid)
            !
            structure_factor_coarse = 0._q
GPACC       !$acc enter data copyin(structure_factor_coarse) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: structure_factor_coarse) if(offload_on)
            !
            do nt = 1, t_info%ntyp
                !
                call cp_grid(coulomb_pot%kernel_truncate%fine_grid, coulomb_pot%kernel_truncate%coarse_grid, &
                             coulomb_pot%kernel_truncate%translation_table, &
                             structure_factor(:,nt), &
                             structure_factor_coarse(:,nt))
                !
                call setunb(structure_factor_coarse(:,nt), &
                            coulomb_pot%kernel_truncate%coarse_grid)
            end do
            !
            potential_coarse_periodic = 0._q
GPACC       !$acc enter data copyin(potential_coarse_periodic) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: potential_coarse_periodic) if(offload_on)
            !
            call potential_from_ion_charge(p, t_info, coulomb_pot%kernel_truncate%coarse_grid, &
                                   latt_cur, &
                                   coulomb_pot_periodic, &
                                   structure_factor_coarse,&
                                   potential_coarse_periodic, &
                                   use_grad_kernels=use_grad_kernels, &
                                   over_ions=over_ions)
            !
            call setunb(potential_coarse_periodic, coulomb_pot%kernel_truncate%coarse_grid)
            !
            call rc_add(potential_coarse, 1.0_q, &
                        potential_coarse_periodic, -1.0_q, &
                        potential_coarse, coulomb_pot%kernel_truncate%coarse_grid)
GPACC       !$acc kernels present(potential) __IF_ASYNC__
GPOMP       !$omp target map(present,to: potential) if(offload_on)
            potential = 0._q
GPOMP       !$omp end target
GPACC       !$acc end kernels
            !
            call cpb_grid(coulomb_pot%kernel_truncate%fine_grid, coulomb_pot%kernel_truncate%coarse_grid, coulomb_pot%kernel_truncate%translation_table, &
                          potential_coarse, potential)
            !
            call setunb_compat(potential, coulomb_pot%kernel_truncate%fine_grid)
            !
            potential_periodic = 0._q
GPACC       !$acc enter data copyin(potential_periodic) __IF_ASYNC__
GPOMP       !$omp target enter data map(to: potential_periodic) if(offload_on)
            !
            call potential_from_ion_charge(p, t_info, gridc, latt_cur,&
                                   coulomb_pot_periodic, &
                                   structure_factor,&
                                   potential_periodic, &
                                   use_grad_kernels=use_grad_kernels, &
                                   over_ions=over_ions)
            !
            call setunb(potential_periodic, gridc)
            !
            call rc_add(potential, 1.0_q, &
                        potential_periodic, 1.0_q, &
                        potential, gridc)
GPACC       !$acc exit data delete(potential_coarse) __IF_ASYNC__
GPACC       !$acc exit data delete(structure_factor_coarse) __IF_ASYNC__
GPACC       !$acc exit data delete(potential_coarse_periodic) __IF_ASYNC__
GPACC       !$acc exit data delete(potential_periodic) __IF_ASYNC__
            !
GPOMP       !$omp target exit data map(release: potential_coarse) if(offload_on)
GPOMP       !$omp target exit data map(release: structure_factor_coarse) if(offload_on)
GPOMP       !$omp target exit data map(release: potential_coarse_periodic) if(offload_on)
GPOMP       !$omp target exit data map(release: potential_periodic) if(offload_on)
            !
        else
            !
            call pad3d(coulomb_pot%kernel_truncate, potential, potential_pad, -1)
        end if
        !
GPACC   !$acc exit data delete(structure_factor_pad) __IF_ASYNC__
GPACC   !$acc exit data delete(potential_pad) __IF_ASYNC__
GPACC   !$acc exit data delete(latt_cur_pad) __IF_ASYNC__
        !
GPOMP   !$omp target exit data map(release: structure_factor_pad) if(offload_on)
GPOMP   !$omp target exit data map(release: potential_pad) if(offload_on)
GPOMP   !$omp target exit data map(release: latt_cur_pad) if(offload_on)
        !
IFXW2 GPOMP !$omp target exit data map(release: coulomb_pot_periodic) if(offload_on)
__GPU   call delete_typed_var(coulomb_pot_periodic)
        !
        call setunb(potential,gridc)
   end subroutine

   subroutine potion_particle_mesh_truncated(gridc, p, latt_cur, t_info, &
                                             coulomb_pot, potential, &
                                             pscenc, energy, force)
       type(grid_3d), intent(in)           :: gridc
       type(type_info), intent(in)         :: t_info
       type(potcar)                        :: p(t_info%ntyp)
       type(latt), intent(in)              :: latt_cur
       type(coulomb_potential), intent(in) :: coulomb_pot
       complex(q), intent(out)             :: potential(gridc%rc%np)
       real(q), intent(out)                :: pscenc
       real(q), intent(out)                :: energy
       real(q), optional                   :: force(:,:)

       real(q), allocatable                :: qtot(:)
       complex(q), allocatable             :: charge(:)
       integer                             :: n_type
       real(q) :: zvsum,g,gx,gy,gz,tewen0
       integer :: nt,nis,ni,n,n1,n2,n3,nc,factm

       PROFILING_START('potion_particle_mesh_truncated')

       allocate(charge(gridc%MPLWV))
       allocate(qtot(t_info%nions))
       !
       call pscenc_energy(p, t_info, latt_cur, coulomb_pot, pscenc)
       !
       charge = 0._q; potential = 0._q; qtot = 0._q
       !
       do n_type = 1, t_info%ntyp
          p(n_type)%usespl => p(n_type)%rhospl
          p(n_type)%usez = -p(n_type)%zvalf
          p(n_type)%usecut = p(n_type)%rcutrho
       end do
       !
       call rhoadd(t_info, latt_cur, p, gridc, charge, qtot)
       !
       call pothar_truncated(gridc, latt_cur, coulomb_pot, charge, potential, &
                   energy)
       energy = -1_q * energy - pscenc
       do nt=1,t_info%ntyp
         energy = energy - p(nt)%eself*t_info%vca(nt)*t_info%nityp(nt)
       enddo
       !
       call FFT3D(potential, gridc, 1)
       !
       ! calculate force on each ion
       !
       ! (dE/dR_i,a) = e^2 \int d^3 r (drho_i(|r-R_i|)/dR_i,a) V(r)
       ! i=1..NION
       ! a=1..3 (x,y,z)
       !
       ! drho_i(|r-R_i|)/dR_i,a = -rho'(|r-R_i|) (r_a-R_i,a)/|r-R_i|
       !
       ! Evaluate f_i(r)=d rho/dr*(r-R_i)/|r-R_i| on the grid for each ion
       ! and integrate F_i = \int dr f_i(r) V(r)
       !
       nis = 1
       !
       if (present(force)) then
           type: do nt = 1, t_info%ntyp
              if (.not.associated(p(nt)%usespl)) then
                 nis = nis+t_info%nityp(nt); cycle type
              endif
              ions: do ni=nis,t_info%nityp(nt)+nis-1
                 call rhoder(t_info,latt_cur,p(nt),ni,gridc,qtot,potential,force(1,ni))
              enddo ions
              nis = nis+t_info%nityp(nt)
           enddo type
           force=force*latt_cur%omega
       end if
       ! take the potential back to reciprocal space
       call rl_add(potential,1._q/gridc%nplwv,potential,0.0_q,potential,gridc)
       call FFT3D(potential,gridc,-1)
       call setunb_compat(potential,gridc)

       PROFILING_STOP('potion_particle_mesh_truncated')
   end subroutine

   subroutine rescaling_factor(kernel_truncate, rescale)
       type(kernel_truncation) :: kernel_truncate
       real(q) :: rescale

       if (kernel_truncate%dimensionality == 0) then
           !
           rescale = real(kernel_truncate%padding_factor,q)**3.0_q
       else if (kernel_truncate%dimensionality == 2) then
           !
           rescale = real(kernel_truncate%padding_factor,q)
       end if
   end subroutine


!************************ SUBROUTINE POTION_PARTICLE_MESH **************
!
!***********************************************************************

      subroutine potion_particle_mesh_(gridc,p, latt_cur, t_info, cvps, pscenc, tewen, for)
      use ini
      use prec
      use mgrid
      use pseudo
      use lattice
      use poscar
      use constant
      use charge
      implicit none
      type (grid_3d) gridc
      type (type_info) t_info
      type (potcar) p(t_info%ntyp)
      type (latt) latt_cur
      complex(q) cvps(gridc%rc%np)
      real(q) pscenc,tewen
      real(q), optional :: for(:,:)
      ! local variables
      real(q) zvsum,g,gx,gy,gz,tewen0
      integer nt,nis,ni,n,n1,n2,n3,nc,factm
      real(q), allocatable :: qtot(:)
      complex(q),allocatable :: chg(:)

      PROFILING_START('potion_particle_mesh_')

      zvsum=0
      do nt=1,t_info%ntyp
         zvsum=zvsum+p(nt)%zvalf*t_info%nityp(nt)*t_info%vca(nt)
      enddo

      pscenc=0
      do nt=1,t_info%ntyp
         pscenc=pscenc+p(nt)%pscore*t_info%vca(nt)*t_info%nityp(nt)*(zvsum/latt_cur%omega)
      enddo

      tewen=-pscenc
      do nt=1,t_info%ntyp
         tewen=tewen-p(nt)%eself*t_info%vca(nt)*t_info%nityp(nt)
      enddo

      cvps =0._q
      allocate(chg(gridc%mplwv))
      chg=0

      ! add pseudo chg.dens. and fft to recipr.space
      do nt=1,t_info%ntyp
         p(nt)%usespl=>p(nt)%rhospl
         p(nt)%usez=-p(nt)%zvalf
         p(nt)%usecut=p(nt)%rcutrho
      end do

      allocate(qtot(t_info%nions))

      call rhoadd(t_info,latt_cur,p,gridc,chg,qtot)

      ! multiply by 4pi/g^2
      ! loop over recipr.lattice points
      n=0
      tewen0=0._q
      col: do nc=1,gridc%rc%ncol
         n2= gridc%rc%i2(nc)
         n3= gridc%rc%i3(nc)
         row: do n1=1,gridc%rc%nrow
            n=n+1
            SETFACT1
            SETFACT

            !=======================================================================
            ! calculate the magnitude of the reciprocal lattice vector
            !=======================================================================
            gx= gridc%lpctx(n1)*latt_cur%b(1,1)+gridc%lpcty(n2)*latt_cur%b(1,2)+gridc%lpctz(n3)*latt_cur%b(1,3)
            gy= gridc%lpctx(n1)*latt_cur%b(2,1)+gridc%lpcty(n2)*latt_cur%b(2,2)+gridc%lpctz(n3)*latt_cur%b(2,3)
            gz= gridc%lpctx(n1)*latt_cur%b(3,1)+gridc%lpcty(n2)*latt_cur%b(3,2)+gridc%lpctz(n3)*latt_cur%b(3,3)

            g=sqrt(gx**2+gy**2+gz**2)*2*pi
            if ( (gridc%lpctx(n1)/=0) .or. (gridc%lpcty(n2)/=0) .or. (gridc%lpctz(n3)/=0) ) then
!              cvps(n) = cvps(n) + ((edeps/latt_cur%omega)*chg(n)/g**2)
               cvps(n)=((edeps/latt_cur%omega)*chg(n)/g**2)
            else
               cvps(n)=0._q
            endif
            tewen0=tewen0+MULFACT cvps(n)*conjg(chg(n))*0.5_q

         enddo row
      enddo col
      CALLMPI( m_sum_d(gridc%comm,tewen0,1))
      tewen=tewen+tewen0

      call setunb(cvps,gridc)
      deallocate(chg)

      if (present(for)) then
         ! bring the potential to real space
         call fft3d(cvps(1),gridc,1)
         ! calculate force on each ion
         !
         ! (de/dr_i,a) = e^2 \int d^3 r (drho_i(|r-r_i|)/dr_i,a) v(r)
         ! i=1..nion
         ! a=1..3 (x,y,z)
         !
         ! drho_i(|r-r_i|)/dr_i,a = -rho'(|r-r_i|) (r_a-r_i,a)/|r-r_i|
         !
         ! evaluate f_i(r)=d rho/dr*(r-r_i)/|r-r_i| on the grid for each ion
         ! and integrate f_i = \int dr f_i(r) v(r)
         nis=1
         type: do nt=1,t_info%ntyp
            if (.not.associated(p(nt)%usespl)) then
               nis = nis+t_info%nityp(nt); cycle type
            endif
            ions: do ni=nis,t_info%nityp(nt)+nis-1
               call rhoder(t_info,latt_cur,p(nt),ni,gridc,qtot,cvps,for(1,ni))
            enddo ions
            nis = nis+t_info%nityp(nt)
         enddo type
         for=for*latt_cur%omega
         ! take the potential back to reciprocal space
         call rl_add(cvps(1),1._q/gridc%nplwv,cvps(1),0.0_q,cvps(1),gridc)
         call fft3d(cvps(1),gridc,-1)
         call setunb_compat(cvps(1),gridc)
      endif

      deallocate(qtot)

      PROFILING_STOP('potion_particle_mesh_')

      return
      end subroutine potion_particle_mesh_

end module

!************************ SUBROUTINE POTHAR ****************************
!
!> This subroutine calculates the hartree potential from the electronic
!> charge density.
!>
!> The correction to the
!> total energy due to overcounting the hartree energy on summing the
!> electronic eigenvalues is also computed (hartree contribution to the
!> total energy = 0.5*sum (vh(g)*rho(-g)). The sum of eigenvalues gives
!> sum (vh(g)*rho(-g)) where vh(g) is the hartree potential at wavevector
!> g and rho(g) is the charge density at wavevector g).
!
!***********************************************************************

      SUBROUTINE POTHAR_(GRIDC,LATT_CUR,CHTOT,CVD,DENC)

__GPU USE moffload_struct_def

      USE prec
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) CVD(GRIDC%RC%NP),CHTOT(GRIDC%RC%NP)

      PROFILING_START('pothar_')

      DENC=0._q
!=======================================================================
! scale the hartree potential by edeps divided by the volume of the unit
! cell
!=======================================================================
      SCALE=EDEPS/LATT_CUR%OMEGA/TPI**2
!=======================================================================
! calculate the hartree potential on the grid of reciprocal lattice
! vectors and the correction to the total energy
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NI,FACTM,GX,GY,GZ,GSQU) &
GPACC !$ACC PRESENT(GRIDC,CVD,CHTOT,LATT_CUR) __IF_ASYNC__ DEFAULT(none)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP PRIVATE(N2,N3,NI,FACTM,GX,GY,GZ,GSQU) MAP(PRESENT,TO:GRIDC,CVD,CHTOT,LATT_CUR) IF(OFFLOAD_ON)
      col: DO NC=1,GRIDC%RC%NCOL
__CPU N2= GRIDC%RC%I2(NC)
__CPU N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
__GPU N2= GRIDC%RC%I2(NC)
__GPU N3= GRIDC%RC%I3(NC)

        NI=(NC-1)*GRIDC%RC%NROW+N1
        SETFACT1
        SETFACT

        GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
        GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
        GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

        GSQU=GX**2+GY**2+GZ**2
!=======================================================================
! since the G=0 coulomb contributions to the hartree, ewald and
! electron-ion energies are individually divergent but together sum to
! zero, set the hartree potential at G=0 to zero.
!=======================================================================
        IF ((GRIDC%LPCTX(N1)==0).AND.(GRIDC%LPCTY(N2)==0).AND.(GRIDC%LPCTZ(N3)==0)) &
     & THEN
          CVD(NI)=(0.0_q,0.0_q)
        ELSE
          CVD(NI)=CHTOT(NI)/GSQU*SCALE
        ENDIF
      ENDDO row
      ENDDO col
      CALL SETUNB(CVD,GRIDC)
!=======================================================================
! calculate the correction to the total energy
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(force:2) REDUCTION(+:DENC) PRIVATE(DUM,N2,N3,NI,FACTM) &
GPACC !$ACC PRESENT(GRIDC,CVD,CHTOT) __IF_ASYNC__ DEFAULT(none)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) REDUCTION(+:DENC) &
GPOMP !$OMP PRIVATE(DUM,N2,N3,NI,FACTM) MAP(PRESENT,TO:GRIDC,CVD,CHTOT) IF(OFFLOAD_ON)
      col2: DO NC=1,GRIDC%RC%NCOL
__CPU N2= GRIDC%RC%I2(NC)
__CPU N3= GRIDC%RC%I3(NC)
      row2: DO N1=1,GRIDC%RC%NROW
__GPU N2= GRIDC%RC%I2(NC)
__GPU N3= GRIDC%RC%I3(NC)

        NI=(NC-1)*GRIDC%RC%NROW+N1
        SETFACT1
        SETFACT

        DUM=MULFACT CVD(NI)*CONJG(CHTOT(NI))
        DENC=DENC+DUM
      ENDDO row2
      ENDDO col2
      DENC=-DENC/2
      CALLMPI( M_sum_d(GRIDC%COMM,DENC,1))

      PROFILING_STOP('pothar_')

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE POTION ****************************
!
!> This subroutine calculates the pseudopotential and its derivatives
!> multiplied by the partial structure-factors
!> on the grid of reciprocal lattice vectors.
!
!***********************************************************************

      SUBROUTINE POTION_(GRIDC,P,LATT_CUR,T_INFO,CVPS,CDVPS,CSTRF,PSCENC)

__GPU USE moffload_struct_def

      USE prec
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      COMPLEX(q) CVPS(GRIDC%RC%NP),CDVPS(GRIDC%RC%NP)

      PROFILING_START('potion_')

!=======================================================================
! calculate the contribution to the total energy from the non-coulomb
! part of the g=0 component of the pseudopotential and the force on the
! unit cell due to the change in this energy as the size of the cell
! changes
!=======================================================================
      ZVSUM=0
      DO NT=1,T_INFO%NTYP
         ZVSUM=ZVSUM+P(NT)%ZVALF*T_INFO%NITYP(NT)*T_INFO%VCA(NT)
      ENDDO

      PSCENC=0
      DO NT=1,T_INFO%NTYP
         PSCENC=PSCENC+P(NT)%PSCORE*T_INFO%VCA(NT)*T_INFO%NITYP(NT)*(ZVSUM/LATT_CUR%OMEGA)
      ENDDO

GPACC !$ACC KERNELS PRESENT(CVPS,CDVPS) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO: CVPS,CDVPS) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CVPS =0
      CDVPS=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!=======================================================================
! loop over all types of atoms
! multiply structur factor by local pseudopotential
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP

      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS
      ZZ=  -4*PI*P(NT)%ZVALF*FELECT

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,N,GX,GY,GZ,G,I,REM,VPST,DVPST) &
GPACC !$ACC PRESENT(GRIDC,CVPS,CDVPS,LATT_CUR,P,CSTRF) __IF_ASYNC__ DEFAULT(none)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP PRIVATE(N2,N3,N,GX,GY,GZ,G,I,REM,VPST,DVPST) &
GPOMP !$OMP MAP(PRESENT,TO:GRIDC,CVPS,CDVPS,LATT_CUR,P,CSTRF) IF(OFFLOAD_ON)

      col: DO NC=1,GRIDC%RC%NCOL
__CPU N2= GRIDC%RC%I2(NC)
__CPU N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
__GPU N2= GRIDC%RC%I2(NC)
__GPU N3= GRIDC%RC%I3(NC)
        N=(NC-1)*GRIDC%RC%NROW+N1
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*2*PI
        IF ( ( (GRIDC%LPCTX(N1)/=0) .OR. (GRIDC%LPCTY(N2)/=0) .OR. &
     &         (GRIDC%LPCTZ(N3)/=0) ) .AND. (G<PSGMA2) ) THEN
!=======================================================================
! convert the magnitude of the reciprocal lattice vector to a position
! in the pseudopotential arrays and interpolate the pseudopotential and
! its derivative
!=======================================================================
        I  =INT(G*ARGSC)+1
        REM=G-P(NT)%PSP(I,1)
        VPST =(P(NT)%PSP(I,2)+REM*(P(NT)%PSP(I,3)+ &
     &                     REM*(P(NT)%PSP(I,4)  +REM*P(NT)%PSP(I,5))))
        DVPST= P(NT)%PSP(I,3)+REM*(P(NT)%PSP(I,4)*2+REM*P(NT)%PSP(I,5)*3)

        CVPS (N)=CVPS (N)+( VPST+ ZZ / G**2)  /LATT_CUR%OMEGA*CSTRF(N,NT)
        CDVPS(N)=CDVPS(N)+( DVPST- 2*ZZ /G**3)/LATT_CUR%OMEGA*CSTRF(N,NT)
        ELSE
        CVPS (N)=0._q
        CDVPS(N)=0._q
        ENDIF

      ENDDO row
      ENDDO col
      ENDDO typ
! set local potential to zero
!      CVPS=0
!      CDVPS=0
      CALL SETUNB(CVPS,GRIDC)
      CALL SETUNB(CDVPS,GRIDC)

      PROFILING_STOP('potion_')

      RETURN
      END SUBROUTINE
