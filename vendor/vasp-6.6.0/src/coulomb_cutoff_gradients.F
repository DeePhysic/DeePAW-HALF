#include "symbol.inc"

module coulomb_cutoff_gradients

__GPU use moffload_struct_def
__GPU use moffload
      use prec, only: q
      use mgrid, only: grid_3d, inilgrd, gen_rc_sub_grid, gen_rc_grid
      use lattice, only: latt, lattic, dirkar
      use pot_electrostat
      use pot_struct_def
      use pseudo_struct_def, only: potcar
      use poscar_struct_def, only: type_info
      use mgrid_struct_def, only: transit
      use base, only: info_struct
      use constant, only: EDEPS, PI, TPI, FELECT, CITPI
      use coulomb_cutoff

      implicit none

contains

    !> This routine computes the local component of the forces for the
    !> coulomb truncation method. It does so by taking into account
    !> the interaction of the entire padded cell when computing the
    !> derivatives. That is, the force is computed as
    !> Sum_g(padded) ig_padded exp(-ig_padded r_padded) rho(g_padded)
    !> V(g_padded)
    subroutine truncate_forloc(gridc, p, t_info, latt_cur, &
                               coulomb_pot, charge_density, &
                               structure_factor, forces)
       type(grid_3d), intent(in)            :: gridc
       type(latt), intent(in)               :: latt_cur
       type(latt)                           :: latt_cur_pad
       type(type_info), intent(in)          :: t_info
       type(type_info)                      :: t_info_pad
       type(potcar), intent(in)             :: p(t_info%ntyp)
       type(coulomb_potential), intent(in)  :: coulomb_pot
       complex(q)                           :: charge_density(:)
       complex(q)                           :: structure_factor(:,:)
       real(q), intent(out)                 :: forces(:,:)
       integer                              :: i, nis, nt, n1, n2, n3
       integer                              :: nc, ni, nion, ni_add
       integer                              :: cg_ncol, cg_nrow
       complex(q)                           :: for
       real(q)                              :: for1, for2, for3
       real(q)                              :: dummy, factm
       complex(q)                           :: charge
       real(q), allocatable                 :: padded_posion(:,:)
       complex(q), allocatable              :: structure_factor_pad(:,:)
       complex(q), allocatable              :: potential_pad(:)
       complex(q), allocatable              :: padded_charge_density(:)
       complex(q), allocatable              :: coarse_charge_density(:)
       real(q)                              :: rescale

       PROFILING_START('truncate_forloc')
       !
       allocate(structure_factor_pad(coulomb_pot%kernel_truncate%fine_grid%mplwv, t_info%ntyp))
       allocate(potential_pad(coulomb_pot%kernel_truncate%fine_grid%mplwv))
       allocate(padded_charge_density(coulomb_pot%kernel_truncate%fine_grid%mplwv))
       !
       structure_factor_pad = 0._q
       potential_pad = 0._q
       padded_charge_density = 0._q
       !
       GPACC !$acc enter data copyin(structure_factor_pad) __IF_ASYNC__
       GPACC !$acc enter data copyin(potential_pad) __IF_ASYNC__
       GPACC !$acc enter data copyin(padded_charge_density) __IF_ASYNC__
       GPOMP !$omp target enter data map(to:structure_factor_pad,potential_pad,padded_charge_density) IF(OFFLOAD_ON)
       !
       if (coulomb_pot%kernel_truncate%coarsen_before_pad) then
           !
           allocate(coarse_charge_density(coulomb_pot%kernel_truncate%coarse_grid%mplwv))
           !
           coarse_charge_density = 0._q
           GPACC !$acc enter data copyin(coarse_charge_density) __IF_ASYNC__
           GPOMP !$omp target enter data map(to:coarse_charge_density) IF(OFFLOAD_ON)
           !
       end if
       !
       call rescaling_factor(coulomb_pot%kernel_truncate, rescale)
       !
       call generate_padded_lattice(coulomb_pot%kernel_truncate, latt_cur, &
                                    latt_cur_pad)
       !
       GPACC !$acc enter data copyin(latt_cur_pad) __IF_ASYNC__
       GPOMP !$omp target enter data map(to:latt_cur_pad) IF(OFFLOAD_ON)
       !
       call generate_padded_t_info(coulomb_pot%kernel_truncate, t_info, &
                                   t_info_pad, padded_posion)
       !
       call stufak(coulomb_pot%kernel_truncate%fine_grid, t_info_pad, structure_factor_pad)
       !
       do nt = 1, t_info%ntyp
           !
           call setunb(structure_factor_pad(:,nt), coulomb_pot%kernel_truncate%fine_grid)
       end do
       !
       call setunb(charge_density, gridc)
       !
       if (coulomb_pot%kernel_truncate%coarsen_before_pad) then
           !
           call cp_grid(coulomb_pot%kernel_truncate%fine_grid, coulomb_pot%kernel_truncate%coarse_grid, &
                        coulomb_pot%kernel_truncate%translation_table, &
                        charge_density, coarse_charge_density)
           !
           call setunb(coarse_charge_density, coulomb_pot%kernel_truncate%coarse_grid)
           !
           call pad3d(coulomb_pot%kernel_truncate, &
                      coarse_charge_density, &
                      padded_charge_density, 1, rescale)
       else
           call pad3d(coulomb_pot%kernel_truncate, charge_density, &
                      padded_charge_density, 1, rescale)
       end if
       !
       call setunb(padded_charge_density, coulomb_pot%kernel_truncate%fine_grid)
       !
       nis = 1
       !
       do nt = 1, t_info%ntyp
       !
       ni_add = t_info%nityp(nt)
       !
       ! cast loop bounds in local variables for faster GPU execution
       cg_ncol = coulomb_pot%kernel_truncate%fine_grid%RC%NCOL
       cg_nrow = coulomb_pot%kernel_truncate%fine_grid%RC%NROW
       !
       do nion = nis, ni_add + nis - 1
           !acc kernels present(potential_pad) __IF_ASYNC__
           potential_pad = 0._q
           !acc end kernels
           call potential_from_ion_charge(p, t_info_pad,&
                                          coulomb_pot%kernel_truncate%fine_grid,latt_cur_pad, &
                                          coulomb_pot, &
                                          structure_factor_pad, &
                                          potential_pad, &
                                          use_grad_kernels=.false., &
                                          over_ions=[nion])
           !
           call setunb(potential_pad, coulomb_pot%kernel_truncate%fine_grid)
           !
           for1 = 0.0_q; for2 = 0.0_q; for3 = 0.0_q
           !
           GPACC !$acc parallel loop reduction(+:for1,for2,for3) &
           GPACC !$acc private(n2,n3,ni,factm,for) &
           GPACC !$acc present(coulomb_pot%kernel_truncate%fine_grid,charge_density)
           GPOMP !$omp target teams distribute parallel do reduction(+:for1,for2,for3) &
           GPOMP !$omp private(n2,n3,ni,factm,for) map(present,to:coulomb_pot%kernel_truncate%fine_grid, &
           GPOMP !$omp charge_density) IF(OFFLOAD_ON)
           do nc = 1, cg_ncol
__CPU      n2 = coulomb_pot%kernel_truncate%fine_grid%RC%I2(nc)
__CPU      n3 = coulomb_pot%kernel_truncate%fine_grid%RC%I3(nc)
           do n1 = 1, cg_nrow
__GPU          n2 = coulomb_pot%kernel_truncate%fine_grid%RC%I2(nc)
__GPU          n3 = coulomb_pot%kernel_truncate%fine_grid%RC%I3(nc)
               !
               ni = (nc-1) * coulomb_pot%kernel_truncate%fine_grid%RC%NROW + n1
               !
               SETFACT1
               SETFACT
               !
               for = cmplx(0.0_q, -1.0_q) * potential_pad(ni) &
                   * MULFACT conjg(padded_charge_density(ni))
               !
               for1 = for1 - coulomb_pot%kernel_truncate%fine_grid%LPCTX_(n1) * real(for,q)
               for2 = for2 - coulomb_pot%kernel_truncate%fine_grid%LPCTY_(n2) * real(for,q)
               for3 = for3 - coulomb_pot%kernel_truncate%fine_grid%LPCTZ_(n3) * real(for,q)
           end do
           end do
           !
           forces(1,nion) = for1 * TPI
           forces(2,nion) = for2 * TPI
           forces(3,nion) = for3 * TPI
           enddo
           !
           nis = nis + ni_add
           !
       enddo
       !
       CALLMPI( M_sum_d(coulomb_pot%kernel_truncate%fine_grid%comm, forces(1,1), t_info_pad%nions*3))
       !
       call dirkar(t_info_pad%nions, forces, latt_cur_pad%b)
       !
       GPACC !$acc exit data delete(structure_factor_pad) __IF_ASYNC__
       GPACC !$acc exit data delete(potential_pad) __IF_ASYNC__
       GPACC !$acc exit data delete(padded_charge_density) __IF_ASYNC__
       GPACC !$acc exit data delete(latt_cur_pad) __IF_ASYNC__
       GPOMP !$omp target exit data map(release:structure_factor_pad,&
       GPOMP !$omp potential_pad,padded_charge_density,latt_cur_pad) IF(OFFLOAD_ON)
       !
       if (coulomb_pot%kernel_truncate%coarsen_before_pad) then
           !
           GPACC !$acc exit data delete(coarse_charge_density) __IF_ASYNC__
           GPOMP !$omp target exit data map(release:coarse_charge_density) IF(OFFLOAD_ON)
           !
       end if
       !
       PROFILING_STOP('truncate_forloc')
    end subroutine
end module
