#include "symbol.inc"
!> Implements improved tetrahedron method for BZ integration described in
!> Bloechl, et al., Phys. Rev. B, 49, 16223 (1994)
!>
!> Usage: First create a tetrahedron calculator using the new_calculator
!> method passing all the relevant quantities from VASP. Then use the
!> member functions of the calculator object to evaluate dos, number of
!> electrons, Fermi energy, or Fermi weights as needed.
!>
!> Internally the calculator object creates a list of all tetrahedra to
!> loop over which is then used by all the other routines. You can also
!> provide a smearing width, which will expand the energy on a grid
!> around the Fermi energy.
module tet

   use fermi_energy, only: smearing_base
   use mpimy, only: communic
   use prec, only: q

   implicit none

   enum, bind(c)
       enumerator :: no_degeneracy
       enumerator :: lowest_3_energies_degenerate
       enumerator :: highest_3_energies_degenerate
       enumerator :: all_energies_degenerate
   end enum
   real(q), parameter :: degeneracy_threshold = 1e-8_q
   real(q), parameter :: converged_num_electron = 1e-10_q
   real(q), parameter :: minimal_energy_difference = 1e-14_q

   real(q), parameter :: grid_points(21) = [ &
      -5.765539551081341_q, &
      -4.094770202368816_q, &
      -3.179348013945789_q, &
      -2.5367996436014075_q, &
      -2.033013883417398_q, &
      -1.6111385714884052_q, &
      -1.2414103917154518_q, &
      -0.9059517905313856_q, &
      -0.5928144852051971_q, &
      -0.2932063844006723_q, &
      0.0_q, &
      0.2932063844006724_q, &
      0.5928144852051971_q, &
      0.9059517905313855_q, &
      1.2414103917154518_q, &
      1.6111385714884052_q, &
      2.033013883417398_q, &
      2.536799643601407_q, &
      3.1793480139457895_q, &
      4.094770202368816_q, &
      5.765539551081341_q &
   ]
   real(q), parameter :: grid_weights(21) = [ &
      0.016017228257774137_q, &
      0.03695378977085292_q, &
      0.057134425426857156_q, &
      0.07610011362837935_q, &
      0.09344442345603382_q, &
      0.10879729916714831_q, &
      0.12183141605372842_q, &
      0.13226893863333739_q, &
      0.13988739479107312_q, &
      0.14452440398997007_q, &
      0.14608113364969047_q, &
      0.14452440398997007_q, &
      0.13988739479107312_q, &
      0.13226893863333739_q, &
      0.12183141605372842_q, &
      0.10879729916714831_q, &
      0.09344442345603382_q, &
      0.07610011362837935_q, &
      0.057134425426857156_q, &
      0.03695378977085292_q, &
      0.016017228257774137_q &
   ]

   type tetrahedron
       real(q) energy(4), volume_weight
       integer :: spin = 0, kpoints(4), band = 0
   contains
       procedure :: e => tetrahedron_energy_difference
   end type tetrahedron

   type, private :: energy_grid
       real(q), allocatable :: energies(:), weights(:)
   end type energy_grid

   type, extends(smearing_base) :: tetrahedron_calculator
       type(tetrahedron), allocatable :: tetra(:)
       real(q), allocatable :: volume_weight_tetra(:)
       type(energy_grid) :: grid
       type(communic) comm_tet
       logical :: use_deltas = .false.
       real(q), private :: lower_bound, energy_window
       real(q) :: num_electrons = -1
   contains
       procedure :: dos => calculate_dos
       procedure :: dos_with_weight
       procedure :: fermi_weights => calculate_fermi_weights
       procedure :: number_electrons => calculate_number_electrons
   end type tetrahedron_calculator

   type fermi_weights_state
       real(q) :: fermi_energy
       real(q), allocatable :: weight_kpoints(:)
       logical :: use_bloechl = .true.
       integer :: unit_out = -1
   end type fermi_weights_state

   type, private :: degeneracy_state
       logical :: warned_about_issues = .false.
       type(tetrahedron), allocatable :: degenerate_tetra(:)
   end type degeneracy_state

   type, private :: weight_sums
       real(q) :: eigenvalues = 0.0_q, bloechl_corrections = 0.0_q, electrons = 0.0_q
   end type weight_sums

   type dos_helper
      ! integration boundaries
      ! first: first index above the lower end of the range spanned by the tetrahedron
      ! last: last index below the upper end of the range spanned by tetrahedron
      integer :: first = 1, last = 0
      ! flat band correction is required if no energy point within tetrahedron
      real(q) :: flat_band_correction = 0.0_q
      ! constant in range E1 < E <= E2
      real(q) :: c21 = dzero
      ! constants in range E2 < E <= E3
      real(q) :: c32(4) = dzero
      ! constant in range E3 < E <= E4
      real(q) :: c43 = dzero
   end type dos_helper

contains

   !> Create a new calculator instance splitting the band structure into individual
   !> tetrahedra.
   type(tetrahedron_calculator) function new_calculator( &
         band_structure, spin_degeneracy, volume_weights, kpoint_indices, &
         comm, num_electrons, smearing, use_deltas, npoints) &
         result (calculator)
      use mpimy
      real(q), intent(in) :: band_structure(:,:,:)  !< dimensions are (bands,kpoints,spin)
      real(q), intent(in) :: spin_degeneracy  !< either 1 or 2 depending how many electrons are in each band
      real(q), intent(in) :: volume_weights(:)  !< volume and weight of every tetrahedron
      integer, intent(in) :: kpoint_indices(:,:)  !< maps all tetrahedron to their respective k-points
      type(communic), intent(in), optional :: comm  !< parallelize over this communicator
      real(q), intent(in), optional :: num_electrons  !< the number of electrons in the system
      real(q), intent(in), optional :: smearing  !< specifies the smearing width
      logical, intent(in), optional :: use_deltas  !< use the exact derivatives
      integer, intent(in), optional :: npoints !< number of grid points for the Gauss-Legendre integration grid
      integer spin, ii, band, bands_for_checksum, tetra
      real(q) checksum
      !
      call calculator%initialize(band_structure, spin_degeneracy)
      calculator%grid = new_grid(smearing,npoints)
      calculator%lower_bound = minval(band_structure)
      calculator%energy_window = maxval(band_structure) - calculator%lower_bound
      if (present(comm)) then
         calculator%comm_tet = comm
      else
#ifdef MPI
         calculator%comm_tet%mpi_comm = MPI_COMM_SELF
         call M_initc(calculator%comm_tet)
#endif
      end if
      if (present(num_electrons)) then
         calculator%num_electrons = num_electrons
         bands_for_checksum = ceiling(0.5_q * num_electrons)
      else
         bands_for_checksum = 0
      end if
      if (present(use_deltas)) calculator%use_deltas = use_deltas
      allocate(calculator%tetra(size(band_structure, 1) * size(kpoint_indices, 2) * size(band_structure, 3)))
      checksum = 0
      tetra = 0
      do spin = 1, size(band_structure, 3)
         do ii = 1, size(kpoint_indices, 2)
            do band = 1, size(band_structure, 1)
               tetra = tetra + 1
               calculator%tetra(tetra) = new_tetrahedron( &
                  spin_degeneracy * volume_weights(ii), &
                  band_structure(band, kpoint_indices(:,ii), spin), &
                  kpoint_indices(:,ii) &
               )
               calculator%tetra(tetra)%band = band
               calculator%tetra(tetra)%spin = spin
               if (band <= bands_for_checksum) &
                  checksum = checksum + spin_degeneracy * volume_weights(ii)
            end do
         end do
      end do
      calculator%attainable_convergence_estimate = abs(nint(checksum) - checksum)
      !
   end function new_calculator

   type(energy_grid) function new_grid(smearing,npoints)
      use gauss_quad, only: gauss_legendre
      integer, intent(in), optional :: npoints !< number of grid points for the Gauss-Legendre integration grid
      real(q), intent(in), optional :: smearing
      integer :: i
      integer npoints_
      real(q) smearing_
      npoints_ = 21
      if (present(npoints)) npoints_ = npoints
      smearing_ = 0.0_q
      if (present(smearing)) smearing_ = smearing
      if (smearing_ > 0.0_q) then
         allocate(new_grid%energies(npoints_))
         allocate(new_grid%weights(npoints_))
         call gauss_legendre(-1.0_q,1.0_q,new_grid%energies,new_grid%weights,npoints_,eps=1e-12_q)
         new_grid%energies = smearing * log((1+new_grid%energies)/(1-new_grid%energies))
         new_grid%weights = 0.5_q * new_grid%weights
      else
         allocate(new_grid%energies(1), new_grid%weights(1))
         new_grid%energies = 0.0_q
         new_grid%weights = 1.0_q
      end if
   end function new_grid

   !> Calculate the DOS with the tetrahedron method.
   subroutine calculate_dos(calculator, energies, dos, integrated_dos)
      use tutor, only: vtutor
      class(tetrahedron_calculator), intent(in) :: calculator
      real(q), intent(in) :: energies(:) !< energies at which to evaluate the dos
      real(q), intent(out) :: dos(:,:) !< dimensions are (size(energy),spin)
      real(q), intent(out) :: integrated_dos(:,:) !< dimensions are (size(energy),spin)
      type(tetrahedron) tetra
      type(dos_helper) helper
      type(communic) comm
      integer ii
      !
      if (size(dos, 1) /= size(energies)) &
         call vtutor%bug("grid and dos incompatible", __FILE__, __LINE__)
      if (size(integrated_dos, 1) /= size(energies)) &
         call vtutor%bug("grid and integrated_dos incompatible", __FILE__, __LINE__)
      !
      dos = dzero
      integrated_dos = dzero
      comm = calculator%comm_tet
      ! loop over tetrahedrons
      do ii = comm%node_me, size(calculator%tetra), comm%ncpu
         tetra = calculator%tetra(ii)
         helper = new_dos_helper(tetra, energies)
         call add_tetrahedron_dos(tetra, helper, energies, dos(:,tetra%spin))
         call add_tetrahedron_integrated_dos(tetra, helper, energies, integrated_dos(:,tetra%spin))
      end do
      !
      CALLMPI( M_sum_d(comm, dos, size(dos)) )
      CALLMPI( M_sum_d(comm, integrated_dos, size(integrated_dos)) )
      !
   end subroutine calculate_dos

   pure subroutine add_tetrahedron_dos(tetra, helper, energies, dos)
      type(tetrahedron), intent(in) :: tetra
      type(dos_helper), intent(in) :: helper
      real(q), intent(in) :: energies(:)
      real(q), intent(inout) :: dos(:)
      integer ii, lower
      !
      do ii = helper%first, helper%last
         dos(ii) = dos(ii) + evaluate_dos(helper, tetra%energy, energies(ii)) &
            + helper%flat_band_correction
      end do
   end subroutine add_tetrahedron_dos

   pure type(dos_helper) function new_dos_helper(tetra, energies) result(helper)
      type(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: energies(:)
      call find_integration_boundaries(tetra, energies, helper)
      if (helper%flat_band_correction > 0.0_q) return
      call calculate_integration_terms(tetra, helper)
   end function new_dos_helper

   pure subroutine find_integration_boundaries(tetra, energies, helper)
      type(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: energies(:)
      type(dos_helper), intent(inout) :: helper
      integer min_, max_, mid
      helper%first = bisect(energies, tetra%energy(1)) + 1
      helper%last = bisect(energies, tetra%energy(4))
      ! if tetrahedron is between two energy points its weight is assigned to the next point
      if (max(1, helper%last) < min(helper%first, size(energies))) then
         helper%flat_band_correction = tetra%volume_weight / (energies(helper%first) - energies(helper%last))
         helper%last = min(helper%first, size(energies))
      end if
   end subroutine find_integration_boundaries

   pure integer function bisect(energies, search) result (res)
      real(q), intent(in) :: energies(:), search
      integer min_, max_, mid
      if (search < energies(1)) then
         res = 0
      else if (search > energies(size(energies))) then
         res = size(energies)
      else
         min_ = 1
         max_ = size(energies)
         do while (max_ - min_ > 1)
            mid = (min_ + max_) / 2
            if (energies(mid) > search) then
               max_ = mid
            else
               min_ = mid
            end if
         end do
         res = min_
      end if
   end function bisect

   pure subroutine calculate_integration_terms(tetra, helper)
      type(tetrahedron), intent(in) :: tetra
      type(dos_helper), intent(inout) :: helper
      associate(vol => tetra%volume_weight, e21 => tetra%energy(2) - tetra%energy(1), &
         e31 => tetra%energy(3) - tetra%energy(1), e41 => tetra%energy(4) - tetra%energy(1), &
         e32 => tetra%energy(3) - tetra%energy(2), e42 => tetra%energy(4) - tetra%energy(2), &
         e43 => tetra%energy(4) - tetra%energy(3))
         if (e21 > dzero) then
            helper%c21 = vol / (e21 * e31 * e41)
         end if
         if (e32 > dzero) then
            helper%c32(4) = -vol * (e31 + e42) / (e31 * e41 * e32 * e42)
            helper%c32(3) = 3.0_q * vol / (e31 * e41)
         end if
         helper%c32(2) = helper%c32(3) * e21
         helper%c32(1) = helper%c32(2) * e21 / 3.0_q
         if (e43 > dzero) then
            helper%c43 = -vol / (e43 * e42 * e41)
         end if
      end associate
   end subroutine calculate_integration_terms

   pure real(q) function evaluate_dos(helper, energies, energy) result (dos)
      type(dos_helper), intent(in) :: helper
      real(q), intent(in) :: energies(4), energy
      dos = dzero
      associate(e1 => energy - energies(1), e2 => energy - energies(2), &
         e3 => energy - energies(3), e4 => energy - energies(4))
         if (e4 > dzero) then
            dos = 0.0_q
         else if (e3 > dzero) then
            ! Case energy between energySorted(3) and energySorted(4):
            dos = -3.0_q * helper%c43 * e4**2
         else if (e2 > dzero) then
            ! Case energy between energySorted(2) and energySorted(3):
            dos = helper%c32(2) + e2 * (2.0_q * helper%c32(3) + 3.0_q * e2 * helper%c32(4))
         else if (e1 > dzero) then
            ! Case energy between energySorted(1) and energySorted(2):
            dos = 3.0_q * helper%c21 * e1**2
         else
            dos = 0.0_q
         end if
      end associate
   end function evaluate_dos

   !> Compute the weigthed DOS
   !>
   !> The weights have two dimensions that can be used to obtain weighted DOS w.r.t
   !> a projection and with spin resolution
   subroutine dos_with_weight(calculator, energies, weights, weighted_dos)
      class(tetrahedron_calculator), intent(in) :: calculator
      real(q), intent(in) :: energies(:) !< energies at which to evaluate the DOS
      real(q), intent(in) :: weights(:,:,:,:) !< dimensions are (bands,nkpoints,projection,spin)
      real(q), intent(out) :: weighted_dos(:,:,:) !< dimensions are (size(energies),projection,spin)
      !
      real(q), allocatable :: deltas(:,:), tetra_weights(:,:)
      type(tetrahedron) tetra
      type(dos_helper) helper
      type(communic) comm
      integer ii, num_nonzero, num_deltas, num_proj
      !
      if (calculator%use_deltas) then
         num_deltas = 4
      else
         num_deltas = 1
      end if
      num_proj = size(weights, 3)
      allocate(deltas(size(energies), num_deltas), tetra_weights(num_deltas, num_proj))
      weighted_dos = dzero
      comm = calculator%comm_tet
      !
      ! loop over tetrahedrons
      do ii = comm%node_me, size(calculator%tetra), comm%ncpu
         tetra = calculator%tetra(ii)
         helper = new_dos_helper(tetra, energies)
         num_nonzero = helper%last - helper%first + 1
         if (num_nonzero == 0) &
             cycle
         call tetrahedron_weights(tetra, weights, tetra_weights)
         call tetrahedron_deltas(tetra, helper, energies, deltas)
         call dgemm('n', 'n', num_nonzero, num_proj, num_deltas, 1.0_q, deltas(helper%first, 1), size(energies), &
            tetra_weights, num_deltas, 1.0_q, weighted_dos(helper%first, 1, tetra%spin), size(energies))
      end do
      CALLMPI( M_sum_d(comm, weighted_dos, size(weighted_dos)) )
   end subroutine dos_with_weight

   !> Evaluate the occupations at every k point
   !>
   !> Based on the determined Fermi level, assign weights to each individual band to
   !> describing by how many electrons it is occupied. Special care is required when
   !> the bands are degenerate at the Fermi energy, because a tetrahedron may be partially
   !> occupied if all corners lie on this degenerate band. Otherwise the algorithm is
   !> numerically unstable, because there is a step function in the occupations.
   subroutine calculate_fermi_weights(calculator, state, fermi_weights)
      class(tetrahedron_calculator), intent(in) :: calculator
      type(fermi_weights_state), intent(in) :: state  !< describes the Fermi energy and some settings used to determine the weights
      real(q), intent(out) :: fermi_weights(:,:,:) !< dimensions are (band,kpoints,spin)
      type(degeneracy_state) degeneracy
      real(q) integrated_weights(8), error_num_electron
      real(q), allocatable :: energies(:), weights(:,:)
      type(tetrahedron) tetra
      type(weight_sums) sum_
      integer ii, jj, kk
      !
      allocate(energies, mold=calculator%grid%energies)
      energies = state%fermi_energy + calculator%grid%energies
      allocate(weights(8, size(energies)))
      !
      fermi_weights = dzero
      !
      do ii = 1, size(calculator%tetra)
         tetra = calculator%tetra(ii)
         call update_degeneracy(tetra, state, degeneracy)
         weights = tetrahedron_fermi_weights(tetra, energies, state%use_bloechl)
         integrated_weights = matmul(weights, calculator%grid%weights)
         call update_sums(integrated_weights, tetra%energy, sum_)
         do jj = 1, size(tetra%kpoints)
            kk = tetra%kpoints(jj)
            fermi_weights(tetra%band, kk, tetra%spin) = fermi_weights(tetra%band, kk, tetra%spin) &
                + integrated_weights(jj) + integrated_weights(jj + size(tetra%kpoints))
         end do
      end do
      !
      error_num_electron = calculate_error_num_electron(calculator%num_electrons, sum_%electrons)
      if (abs(error_num_electron) > converged_num_electron .and. allocated(degeneracy%degenerate_tetra)) &
         call correct_degeneracy(error_num_electron, degeneracy%degenerate_tetra, fermi_weights)
      if (allocated(state%weight_kpoints)) &
         call normalize_for_spin_and_kpoint_weight(calculator%spin_degeneracy, &
            state%weight_kpoints, fermi_weights)
      call log_sums(state, sum_)
      !
   end subroutine calculate_fermi_weights

   pure real(q) function calculate_error_num_electron(expected_electrons, actual_electrons) result (error)
      real(q), intent(in) :: expected_electrons, actual_electrons
      error = 0
      if (expected_electrons <= 0) &
         ! only correct if we know how many electrons are in the system
         return
      error = actual_electrons - expected_electrons
   end function calculate_error_num_electron

   subroutine update_degeneracy(tetra, state, degeneracy)
       use tutor, only: vtutor
       type(tetrahedron), intent(in) :: tetra
       type(fermi_weights_state), intent(in) :: state
       type(degeneracy_state), intent(inout) :: degeneracy
       logical issue_with_bloechl_correction
       integer ierr
       ierr = degenerate_with_energy(tetra, state%fermi_energy)
       if (ierr == all_energies_degenerate) &
           call add_to_list(tetra, degeneracy%degenerate_tetra)
       !
       issue_with_bloechl_correction = ((ierr == lowest_3_energies_degenerate).or.(ierr == highest_3_energies_degenerate)) &
           .and. state%use_bloechl
       if (issue_with_bloechl_correction .and. .not.degeneracy%warned_about_issues) then
           call vtutor%warning("Fermi level nearly degenerate with energies of tetrahedron.&
               & This may result in numerical instabilities because the density changes by&
               & tiny changes of the Fermi energy. If the calculation converges, the total&
               & energy will still be reliable. If you notice convergence issues, you can&
               & try to\n&
               &*) Use a Gamma centered k-point grid\n&
               &*) Use tetrahedron method without corrections ISMEAR = -4\n&
               &*) Don't use tetrahedron method\n&
               &   ISMEAR = 0 for semiconductors, 1 or 2 for metals")
           degeneracy%warned_about_issues = .true.
       end if
   end subroutine update_degeneracy

   subroutine add_to_list(new_tetra, list)
      type(tetrahedron), intent(in) :: new_tetra
      type(tetrahedron), intent(inout), allocatable :: list(:)
      type(tetrahedron), allocatable :: copy(:)
      if (allocated(list)) then
         allocate(copy(size(list) + 1))
         copy(:size(list)) = list
         copy(size(copy)) = new_tetra
         call move_alloc(copy, list)
      else
         allocate(list(1))
         list(1) = new_tetra
      end if
   end subroutine add_to_list

   subroutine update_sums(weights, energies, sum_)
      real(q), intent(in) :: weights(8), energies(4)
      type(weight_sums), intent(inout) :: sum_
      sum_%eigenvalues = sum_%eigenvalues + dot_product(weights(1:4), energies)
      sum_%bloechl_corrections = sum_%bloechl_corrections + dot_product(weights(5:8), energies)
      sum_%electrons = sum_%electrons + sum(weights(1:4))
   end subroutine update_sums

   subroutine correct_degeneracy(error_num_electron, degenerate_tetra, fermi_weights)
      real(q), intent(in) :: error_num_electron
      type(tetrahedron), intent(in) :: degenerate_tetra(:)
      real(q), intent(inout) :: fermi_weights(:,:,:)
      real(q) excess_fraction
      integer ii, jj
      excess_fraction = 0.25_q * error_num_electron / sum(degenerate_tetra%volume_weight)
         ! factor 0.25, because weight is distributed over 4 k-points
      do ii = 1, size(degenerate_tetra)
         associate(tetra => degenerate_tetra(ii))
            do jj = 1, size(tetra%kpoints)
               fermi_weights(tetra%band, tetra%kpoints(jj), tetra%spin) = &
                  fermi_weights(tetra%band, tetra%kpoints(jj), tetra%spin) - excess_fraction * tetra%volume_weight
            end do
         end associate
      end do
   end subroutine correct_degeneracy

   subroutine normalize_for_spin_and_kpoint_weight(weight_spin, weight_kpoints, fermi_weights)
      use string, only: str
      use tutor, only: vtutor
      real(q), intent(in) :: weight_spin, weight_kpoints(:)
      real(q), intent(inout) :: fermi_weights(:,:,:)
      integer ikpt
      if (size(weight_kpoints) < 4) &
         call vtutor%error('VERY BAD NEWS! internal error in subroutine BZINTS: &
            &Tetrahedron method fails (number of k-points < 4) ' // str(size(weight_kpoints)))
      do concurrent (ikpt = 1:size(fermi_weights,2))
         fermi_weights(:,ikpt,:) = fermi_weights(:,ikpt,:) / &
            (weight_spin * weight_kpoints(ikpt))
      end do
   end subroutine normalize_for_spin_and_kpoint_weight

   subroutine log_sums(state, sum_)
      type(fermi_weights_state), intent(in) :: state
      type(weight_sums), intent(in) :: sum_
      if (state%unit_out >= 0.and.state%unit_out <= 99) then
         write(state%unit_out, '(1x,a,f12.6,a,f12.6,a)') &
            'BZINTS: Fermi energy:', state%fermi_energy, ';', sum_%electrons, ' electrons'
         write(state%unit_out, '(9x,a,e13.6,a,2x,a,f12.6)') &
            'Band energy:', sum_%eigenvalues + sum_%bloechl_corrections, ';', &
            'BLOECHL correction:', sum_%bloechl_corrections
      end if
   end subroutine log_sums


   !> Evaluate how many electrons are in all tetrahedra based on the Fermi energy
   real(q) function calculate_number_electrons(this, fermi_energy) result (num_elec)
      class(tetrahedron_calculator), intent(in) :: this
      real(q), intent(in) :: fermi_energy
      type(tetrahedron) tetra
      type(dos_helper) helper
      type(communic) comm
      real(q), allocatable :: energies(:), int_dos(:), degenerate_volume(:)
      real(q) sum_degen_vol, error_num_elec
      integer ii
      !
      allocate(energies, int_dos, degenerate_volume, mold=this%grid%energies)
      energies = fermi_energy + this%grid%energies
      int_dos = 0
      degenerate_volume = 0
      comm = this%comm_tet
      ! loop over tetrahedrons
      do ii = comm%node_me, size(this%tetra), comm%ncpu
         tetra = this%tetra(ii)
         helper = new_dos_helper(tetra, energies)
         where (degenerate_with_energy(tetra, energies) == all_energies_degenerate)
             degenerate_volume = degenerate_volume + tetra%volume_weight
         else where
             int_dos = int_dos + tetrahedron_integrated_dos(tetra, helper, energies)
         end where
      end do
      num_elec = sum(this%grid%weights * int_dos)
      sum_degen_vol = sum(this%grid%weights * degenerate_volume)
      CALLMPI( M_sum_d(comm, num_elec, 1) )
      CALLMPI( M_sum_d(comm, sum_degen_vol, 1) )
      !
      error_num_elec = this%num_electrons - num_elec
      if (0 < error_num_elec .and. error_num_elec < sum_degen_vol) &
          num_elec = this%num_electrons
      !
   end function calculate_number_electrons

   type(tetrahedron) function new_tetrahedron(volume_weight, energy, kpoints) result (tetra)
      real(q), intent(in) :: volume_weight, energy(4)
      integer, intent(in) :: kpoints(4)
      integer order(4)
      tetra%volume_weight = volume_weight
      call sort(energy, tetra%energy, order)
      tetra%kpoints = kpoints(order)
   end function new_tetrahedron

   subroutine sort(unsorted, sorted, order)
      real(q), intent(in) :: unsorted(4)
      real(q), intent(out) :: sorted(4)
      integer, intent(out), optional :: order(4)
      sorted = unsorted
      if (present(order)) order = [1, 2, 3, 4]
      ! optimal sorting network for vector of length 4
      call reorder(1, 3)
      call reorder(2, 4)
      call reorder(1, 2)
      call reorder(3, 4)
      call reorder(2, 3)
   contains
      subroutine reorder(ii, jj)
         integer, intent(in) :: ii, jj
         integer indx
         real(q) val
         if (sorted(ii) > sorted(jj)) then
            val = sorted(ii)
            sorted(ii) = sorted(jj)
            sorted(jj) = val
            if (present(order)) then
               indx = order(ii)
               order(ii) = order(jj)
               order(jj) = indx
            end if
         end if
      end subroutine reorder
   end subroutine sort


   pure subroutine add_tetrahedron_integrated_dos(tetra, helper, energies, integrated_dos)
      ! This subroutine adds up the contributions to the number
      ! of states for one single tetrahedron
      !
      ! properties of the current tetrahedron
      class(tetrahedron), intent(in) :: tetra
      type(dos_helper), intent(in) :: helper
      ! energies at which the integrated DOS is evaluated
      real(q), intent(in) :: energies(:)
      real(q), intent(inout) :: integrated_dos(:)
      integer ii
      !
      do ii = helper%first, helper%last
         integrated_dos(ii) = integrated_dos(ii) + evaluate_integrated_dos(tetra, helper, energies(ii))
      end do
      integrated_dos(helper%last + 1:) = integrated_dos(helper%last + 1:) + tetra%volume_weight
   end subroutine add_tetrahedron_integrated_dos

   pure function tetrahedron_integrated_dos(tetra, helper, energies) result (integrated_dos)
      ! This subroutine adds up the contributions to the number
      ! of states for one single tetrahedron
      !
      ! properties of the current tetrahedron
      class(tetrahedron), intent(in) :: tetra
      type(dos_helper), intent(in) :: helper
      ! energies at which the integrated DOS is evaluated
      real(q), intent(in) :: energies(:)
      real(q) integrated_dos(size(energies))
      !
      integrated_dos = 0.0_q
      call add_tetrahedron_integrated_dos(tetra, helper, energies, integrated_dos)
   end function tetrahedron_integrated_dos

   pure real(q) function evaluate_integrated_dos(tetra, helper, energy) result (int_dos)
      type(tetrahedron), intent(in) :: tetra
      type(dos_helper), intent(in) :: helper
      real(q), intent(in) :: energy
      associate(e1 => energy - tetra%energy(1), e2 => energy - tetra%energy(2), &
         e3 => energy - tetra%energy(3), e4 => energy - tetra%energy(4))
         if (e4 > dzero) then
            ! Case energy greater than tetra%energy(4) (might probably happen for last element):
            int_dos = tetra%volume_weight
         else if (e3 > dzero) then
            ! Case energy between tetra%energy(3) and tetra%energy(4):
            int_dos = tetra%volume_weight - helper%c43 * e4**3
         else if (e2 > dzero) then
            ! Case energy between tetra%energy(2) and tetra%energy(3):
            int_dos = helper%c32(1) + e2 * (helper%c32(2) + e2 * (helper%c32(3) + e2 * helper%c32(4)))
         else if (e1 > dzero) then
            ! Case energy between tetra%energy(1) and tetra%energy(2):
            int_dos = helper%c21 * e1**3
         else
            int_dos = 0.0_q
         end if
      end associate
   end function evaluate_integrated_dos

   elemental real(q) function tetrahedron_energy_difference(tetra, higher, lower) result (diff)
      class(tetrahedron), intent(in) :: tetra
      integer, intent(in) :: higher, lower
      diff = tetra%energy(higher) - tetra%energy(lower)
   end function tetrahedron_energy_difference

   pure function tetrahedron_fermi_weights(tetra, energies, use_bloechl) result (weights)
      ! This routine makes the weights for integration up to fermi_energy for
      ! one single tetrahedron ('contributions to occupation numbers').
      !
      ! properties of the current tetrahedron
      class(tetrahedron), intent(in) :: tetra
      ! energies: possible choices for the Fermi energy or integration grid
      real(q), intent(in) :: energies(:)
      ! whether Bloechl corrections should be applied to the weights
      logical, intent(in) :: use_bloechl
      !
      ! weights(1-4,:): 'normal' weights
      ! weights(5-8,:): Bloechl-corrections to the weights
      real(q) :: weights(8, size(energies))
      integer ii
      !
      weights = dzero
      !
      do ii = 1, size(energies)
         ! Lowest energy still >= fermi_energy --> no contributions to weights
         if (tetra%energy(1) >= energies(ii)) cycle
         if (energies(ii) >= tetra%energy(4)) then
            ! Highest energy still <= fermi_energy --> just add up full weight
            weights(1:4,ii) = 0.25_q * tetra%volume_weight
         else if (energies(ii) > tetra%energy(3)) then
            weights(1:4,ii) = fermi_weights_between_34(tetra, energies(ii))
         else if (energies(ii) > tetra%energy(2)) then
            weights(1:4,ii) = fermi_weights_between_23(tetra, energies(ii))
         else
            weights(1:4,ii) = fermi_weights_between_12(tetra, energies(ii))
         end if
         !
         if (use_bloechl) then
            weights(5:8,ii) = bloechl_corrections(tetra, energies(ii))
         end if
      end do
   end function tetrahedron_fermi_weights

   pure function fermi_weights_between_34(tetra, fermi_energy) result (weights)
      type(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: fermi_energy
      real(q) weights(4), cc
      cc = 0.25_q * tetra%volume_weight * (tetra%energy(4) - fermi_energy)**3 / &
         (tetra%e(4,1) * tetra%e(4,2) * tetra%e(4,3))
      weights(1) = 0.25_q * tetra%volume_weight - cc * (tetra%energy(4) - fermi_energy) / tetra%e(4,1)
      weights(2) = 0.25_q * tetra%volume_weight - cc * (tetra%energy(4) - fermi_energy) / tetra%e(4,2)
      weights(3) = 0.25_q * tetra%volume_weight - cc * (tetra%energy(4) - fermi_energy) / tetra%e(4,3)
      weights(4) = 0.25_q * tetra%volume_weight - cc * (4.0_q - (tetra%energy(4) - fermi_energy) * &
         (done / tetra%e(4,1) + done / tetra%e(4,2) + done / tetra%e(4,3)))
   end function fermi_weights_between_34

   pure function fermi_weights_between_23(tetra, fermi_energy) result (weights)
      type(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: fermi_energy
      real(q) weights(4), c1, c2, c3
      c1 = 0.25_q * tetra%volume_weight * (fermi_energy - tetra%energy(1))**2 / (tetra%e(4,1) * tetra%e(3,1))
      c2 = 0.25_q * tetra%volume_weight * (fermi_energy - tetra%energy(1)) * &
         (fermi_energy - tetra%energy(2)) * (tetra%energy(3) - fermi_energy) / &
         (tetra%e(4,1) * tetra%e(3,2) * tetra%e(3,1))
      c3 = 0.25_q * tetra%volume_weight * (fermi_energy - tetra%energy(2))**2 * (tetra%energy(4) - fermi_energy) / &
         (tetra%e(4,2) * tetra%e(3,2) * tetra%e(4,1))
      weights(1) = c1 + (c1 + c2) * (tetra%energy(3) - fermi_energy) / tetra%e(3,1) + &
         (c1 + c2 + c3) * (tetra%energy(4) - fermi_energy) / tetra%e(4,1)
      weights(2) = c1 + c2 + c3 + (c2 + c3) * (tetra%energy(3) - fermi_energy) / tetra%e(3,2) + &
         c3 * (tetra%energy(4) - fermi_energy) / tetra%e(4,2)
      weights(3) = (c1 + c2) * (fermi_energy - tetra%energy(1)) / tetra%e(3,1) + &
         (c2 + c3) * (fermi_energy - tetra%energy(2)) / tetra%e(3,2)
      weights(4) = (c1 + c2 + c3) * (fermi_energy - tetra%energy(1)) / tetra%e(4,1) + &
         c3 * (fermi_energy - tetra%energy(2)) / tetra%e(4,2)
   end function fermi_weights_between_23

   pure function fermi_weights_between_12(tetra, fermi_energy) result (weights)
      type(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: fermi_energy
      real(q) weights(4), cc
      cc = 0.25_q * tetra%volume_weight * (fermi_energy - tetra%energy(1))**3 / &
         (tetra%e(2,1) * tetra%e(3,1) * tetra%e(4,1))
      weights(1) = cc * (4.0_q - (fermi_energy - tetra%energy(1)) * &
         (done / tetra%e(2,1) + done / tetra%e(3,1) + done / tetra%e(4,1)))
      weights(2) = cc * (fermi_energy - tetra%energy(1)) / tetra%e(2,1)
      weights(3) = cc * (fermi_energy - tetra%energy(1)) / tetra%e(3,1)
      weights(4) = cc * (fermi_energy - tetra%energy(1)) / tetra%e(4,1)
   end function fermi_weights_between_12

   pure function bloechl_corrections(tetra, fermi_energy) result (weights)
      type(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: fermi_energy
      type(dos_helper) helper
      real(q) weights(4), dos_fermi(1), energies(1)
      integer ii
      energies = fermi_energy
      helper = new_dos_helper(tetra, energies)
      dos_fermi = 0.0_q
      call add_tetrahedron_dos(tetra, helper, energies, dos_fermi)
      do ii = 1, 4
         weights(ii) = 0.025_q * sum(tetra%energy(:) - tetra%energy(ii)) * dos_fermi(1)
      end do
   end function bloechl_corrections

   pure subroutine tetrahedron_weights(tetra, weights, tetra_weights)
      class(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: weights(:,:,:,:)
      real(q), intent(out) :: tetra_weights(:,:)
      !
      if (size(tetra_weights, 1) == 4) then
         tetra_weights = weights(tetra%band, tetra%kpoints, :, tetra%spin)
      else
         ! average over corners of tetrahedron
         tetra_weights(1,:) = sum(weights(tetra%band, tetra%kpoints, :, tetra%spin), dim=1) / 4
      end if
   end subroutine tetrahedron_weights

   pure subroutine tetrahedron_deltas(tetra, helper, energies, deltas)
      ! This routine makes the weights for the density of states at energy for
      ! one single tetrahedron.
      ! properties of the current tetrahedron
      class(tetrahedron), intent(in) :: tetra
      ! energy at which the deltas are evaluated
      real(q), intent(in) :: energies(:)
      ! helper type for dos calculation
      type(dos_helper), intent(in) :: helper
      ! derivatives of the Fermi weights, note: only computed between helper%first and helper%last
      real(q), intent(out) :: deltas(:,:)
      integer ii
      !

      if (size(deltas, 2) == 4) then
         do ii = helper%first, helper%last
            call deltas_single_energy_tetra(tetra, energies(ii), deltas(ii,:))
         end do
      else
         deltas(helper%first:helper%last,1) = 0.0_q
         call add_tetrahedron_dos(tetra, helper, energies, deltas(:,1))
      end if
   end subroutine tetrahedron_deltas

   pure subroutine deltas_single_energy_tetra(tetra, energy, deltas)
       type(tetrahedron), intent(in) :: tetra
       real(q), intent(in) :: energy
       real(q), intent(out) :: deltas(4)
       ! energy out of bounds --> no contributions to deltas
       call deltas_single_energy(tetra%volume_weight, energy, tetra%energy, deltas)
   end subroutine deltas_single_energy_tetra

   pure subroutine deltas_single_energy(volume_weight, energy, energies, deltas)
       real(q), intent(in) :: volume_weight
       real(q), intent(in) :: energy
       real(q), intent(in) :: energies(4)
       real(q), intent(out) :: deltas(4)
       ! energy out of bounds --> no contributions to deltas
       if ((energy <= energies(1)).or.(energy >= energies(4))) then
          deltas = dzero
       else if (energy > energies(3)) then
          call deltaWeights43(volume_weight, energies, energy, deltas)
       else if (energy > energies(2)) then
          call deltaWeights32(volume_weight, energies, energy, deltas)
       else
          call deltaWeights21(volume_weight, energies, energy, deltas)
       end if
   end subroutine deltas_single_energy

   pure subroutine deltaWeights21(volumeWeightTetra, energySorted, Energy, weights)
      ! Case Energy between energySorted(1) and energySorted(2)
      real(q), intent(in) :: volumeWeightTetra, energySorted(4), Energy
      real(q), intent(out) :: weights(4)
      real(q) a21, a31, a41, prod
      real(q) d21, d31, d41, dprod
      associate(vol => 0.25_q * volumeWeightTetra, e1 => energySorted(1), &
         e2 => energySorted(2), e3 => energySorted(3), e4 => energySorted(4))
         d21 = done / (e2 - e1)
         d31 = done / (e3 - e1)
         d41 = done / (e4 - e1)
         a21 = (Energy - e1) * d21
         a31 = (Energy - e1) * d31
         a41 = (Energy - e1) * d41
         prod = a21 * a31 * a41 * vol
         dprod = ( d21 * a31 * a41 +&
                   a21 * d31 * a41 +&
                   a21 * a31 * d41 ) * vol
         !weights(1) = prod * (4.0_q - a21 - a31 - a41)
         weights(1) = prod * (- d21 - d31 - d41) + dprod * (4.0_q - a21 - a31 - a41)
         !weights(2) = prod * a21
         weights(2) = prod * d21 + dprod * a21
         !weights(3) = prod * a31
         weights(3) = prod * d31 + dprod * a31
         !weights(4) = prod * a41
         weights(4) = prod * d41 + dprod * a41
      end associate
   end subroutine deltaWeights21

   pure subroutine deltaWeights32(volumeWeightTetra, energySorted, Energy, weights)
      ! Case Energy between energySorted(2) and energySorted(3)
      real(q), intent(in) :: volumeWeightTetra, energySorted(4), Energy
      real(q), intent(out) :: weights(4)
      real(q) a31, a41, a32, a42, v1, v2, v3
      real(q) d31, d41, d32, d42, dv1, dv2, dv3
      associate(vol => 0.25_q * volumeWeightTetra, e1 => energySorted(1), &
         e2 => energySorted(2), e3 => energySorted(3), e4 => energySorted(4))
         d31 = done / (e3 - e1)
         d41 = done / (e4 - e1)
         d32 = done / (e3 - e2)
         d42 = done / (e4 - e2)
         a31 = (Energy - e1) * d31
         a41 = (Energy - e1) * d41
         a32 = (Energy - e2) * d32
         a42 = (Energy - e2) * d42
         v1 = a31 * a41
         v2 = a31 * a42 * (done - a41)
         v3 = a42 * a32 * (done - a31)
         dv1 = d31 * a41 + a31 * d41
         dv2 = d31 * a42 * (done - a41) +&
               a31 * d42 * (done - a41) +&
               a31 * a42 * (- d41)
         dv3 = d42 * a32 * (done - a31) +&
               a42 * d32 * (done - a31) +&
               a42 * a32 * ( - d31)
         !weights(1) = (v1 * (3.0_q - a31 - a41) +&
         !              v2 * (2.0_q - a31 - a41) +&
         !              v3 * (done - a31)) * vol
         weights(1) = (dv1 * (3.0_q - a31 - a41) + v1 * (- d31 - d41) +&
                       dv2 * (2.0_q - a31 - a41) + v2 * (- d31 - d41) +&
                       dv3 * (done - a31) + v3 * (- d31)) * vol
         !weights(2) = (v1 +&
         !              v2 * (2.0_q - a42) +&
         !              v3 * (3.0_q - a32 - a42)) * vol
         weights(2) = (dv1 +&
                       dv2 * (2.0_q - a42) + v2 * (- d42) +&
                       dv3 * (3.0_q - a32 - a42) + v3 * (- d32 - d42)) * vol
         !weights(3) = (v1 * a31 +&
         !              v2 * a31 +&
         !              v3 * (a31 + a32)) * vol
         weights(3) = (dv1 * a31 + v1 * d31 +&
                       dv2 * a31 + v2 * d31 +&
                       dv3 * (a31 + a32) + v3 * (d31 + d32)) * vol
         !weights(4) = (v1 * a41 +&
         !              v2 * (a41 + a42) +&
         !              v3 * a42) * vol
         weights(4) = (dv1 * a41 + v1 * d41 +&
                       dv2 * (a41 + a42) + v2 * (d41 + d42) +&
                       dv3 * a42 + v3 * d42) * vol
      end associate
   end subroutine deltaWeights32

   pure subroutine deltaWeights43(volumeWeightTetra, energySorted, Energy, weights)
      ! Case Energy between energySorted(3) and energySorted(4)
      real(q), intent(in) :: volumeWeightTetra, energySorted(4), Energy
      real(q), intent(out) :: weights(4)
      real(q) a14, a24, a34, prod
      real(q) d14, d24, d34, dprod
      associate(vol => 0.25_q * volumeWeightTetra, e1 => energySorted(1), &
         e2 => energySorted(2), e3 => energySorted(3), e4 => energySorted(4))
         d14 = done / (e1 - e4)
         d24 = done / (e2 - e4)
         d34 = done / (e3 - e4)
         a14 = (Energy - e4) * d14
         a24 = (Energy - e4) * d24
         a34 = (Energy - e4) * d34
         prod = a14 * a24 * a34 * vol
         dprod = (d14 * a24 * a34+&
                  a14 * d24 * a34+&
                  a14 * a24 * d34) * vol
         !weights(1) = vol - prod * a14
         weights(1) = - (dprod * a14 + prod * d14)
         !weights(2) = vol - prod * a24
         weights(2) = - (dprod * a24 + prod * d24)
         !weights(3) = vol - prod * a34
         weights(3) = - (dprod * a34 + prod * d34)
         !weights(4) = vol - prod * (4.0_q - a14 - a24 - a34)
         weights(4) = - ( dprod * (4.0_q - a14 - a24 - a34) + prod * (- d14 - d24 - d34) )
      end associate
   end subroutine deltaWeights43

   elemental integer function degenerate_with_energy(tetra, energy) result (res)
      class(tetrahedron), intent(in) :: tetra
      real(q), intent(in) :: energy
      logical degeneracy(4)
      degeneracy(1:2) = tetra%energy(1:2) + 0.5_q * degeneracy_threshold > energy
      degeneracy(3:4) = tetra%energy(3:4) - 0.5_q * degeneracy_threshold < energy
      if (degeneracy(1).and.degeneracy(4)) then
         res = all_energies_degenerate
      else if (degeneracy(1).and.degeneracy(3)) then
         res = lowest_3_energies_degenerate
      else if (degeneracy(2).and.degeneracy(4)) then
         res = highest_3_energies_degenerate
      else
         res = no_degeneracy
      end if
   end function degenerate_with_energy

end module tet

!> Implements the weights for the tetrahedron method as described in the paper
!> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
!> Journal of Physics C: Solid State Physics 12, 2991 (1979).
!> These same expressions are implemented in (spglib,kspclib) by Atsushi Togo
!> The implementation here is made by Henrique Miranda
module tet_macdonald
  use prec
  implicit none

  public :: tetralite_delta
  public :: tetralite_delta1
  public :: tetralite_delta2
  public :: tetralite_delta3
  public :: tetralite_delta4
  public :: tetralite_theta

  contains
  pure real(q) function f(n,m,w,eig)
    integer,intent(in) :: n,m
    real(q),intent(in) :: w,eig(4)
    f = ((w - eig(m)) / (eig(n) - eig(m)))
  end function

  pure real(q) function g1(w,eig)
    real(q),intent(in) :: w,eig(4)
    g1 = (3.0_q *&
          f(2, 1, w, eig) *&
          f(3, 1, w, eig) /&
          (eig(4) - eig(1)))
    !g1 = (3.0_q *&
    !      n1(w, eig) /&
    !      (w - eig(1)))
  end function

  pure real(q) function g2(w,eig)
    real(q),intent(in) :: w,eig(4)
    g2 = (3.0_q /&
          (eig(4) - eig(1)) *&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))
  end function

  pure real(q) function g3(w,eig)
    real(q),intent(in) :: w,eig(4)
    g3 = (3.0_q *&
            f(2, 4, w, eig) *&
            f(3, 4, w, eig) /&
            (eig(4) - eig(1)))
    !g3 = (3.0_q *&
    !        (1.0_q - n3(w,eig)) /&
    !        (eig(4) - w))
  end function

  pure real(q) function I10(w,eig)
    real(q),intent(in) :: w,eig(4)
    I10 = (f(1, 2, w, eig) +&
           f(1, 3, w, eig) +&
           f(1, 4, w, eig))/3.0_q
  end function

  pure real(q) function I11(w,eig)
    real(q),intent(in) :: w,eig(4)
    I11 = f(2, 1, w, eig)/3.0_q
  end function

  pure real(q) function I12(w,eig)
    real(q),intent(in) :: w,eig(4)
    I12 = f(3, 1, w, eig)/3.0_q
  end function

  pure real(q) function I13(w,eig)
    real(q),intent(in) :: w,eig(4)
    I13 = f(4, 1, w, eig)/3.0_q
  end function

  pure real(q) function I20(w,eig)
    real(q),intent(in) :: w,eig(4)
    I20 = (f(1, 4, w, eig) +&
           f(1, 3, w, eig) *&
           f(3, 1, w, eig) *&
           f(2, 3, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))/3.0_q
  end function

  pure real(q) function I21(w,eig)
    real(q),intent(in) :: w,eig(4)
    I21 = (f(2, 3, w, eig) +&
           f(2, 4, w, eig) *&
           f(2, 4, w, eig) *&
           f(3, 2, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))/3.0_q
  end function

  pure real(q) function I22(w,eig)
    real(q),intent(in) :: w,eig(4)
    I22 = (f(3, 2, w, eig) +&
           f(3, 1, w, eig) *&
           f(3, 1, w, eig) *&
           f(2, 3, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))/3.0_q
  end function

  pure real(q) function I23(w,eig)
    real(q),intent(in) :: w,eig(4)
    I23 = (f(4, 1, w, eig) +&
           f(4, 2, w, eig) *&
           f(2, 4, w, eig) *&
           f(3, 2, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))/3.0_q
  end function

  pure real(q) function I30(w,eig)
    real(q),intent(in) :: w,eig(4)
    I30 = f(1, 4, w, eig)/3.0_q
  end function

  pure real(q) function I31(w,eig)
    real(q),intent(in) :: w,eig(4)
    I31 = f(2, 4, w, eig)/3.0_q
  end function

  pure real(q) function I32(w,eig)
    real(q),intent(in) :: w,eig(4)
    I32 = f(3, 4, w, eig)/3.0_q
  end function

  pure real(q) function I33(w,eig)
    real(q),intent(in) :: w,eig(4)
    I33 = (f(4, 1, w, eig) +&
           f(4, 2, w, eig) +&
           f(4, 3, w, eig))/3.0_q
  end function

  ! w1 < w < w2
  pure real(q) function n1(w,eig)
    real(q),intent(in) :: w,eig(4)
    n1 = (f(2, 1, w, eig) *&
          f(3, 1, w, eig) *&
          f(4, 1, w, eig))
  end function

  ! w2 < w < w3
  pure real(q) function n2(w,eig)
    real(q),intent(in) :: w,eig(4)
    n2 = (f(4, 2, w, eig) *&
          f(3, 2, w, eig) +&
          f(4, 1, w, eig) *&
          f(2, 4, w, eig) *&
          f(3, 2, w, eig) +&
          f(4, 1, w, eig) *&
          f(3, 1, w, eig) *&
          f(2, 3, w, eig))
  end function

  ! w3 < w < w4
  pure real(q) function n3(w,eig)
    real(q),intent(in) :: w,eig(4)
    n3 = (1.0_q -&
          f(1, 4, w, eig) *&
          f(2, 4, w, eig) *&
          f(3, 4, w, eig))
  end function

  pure real(q) function J11(w,eig)
    real(q),intent(in) :: w,eig(4)
    integer :: f12,f13,f14
    J11 = (1.0_q +&
           f(1, 2, w, eig) +&
           f(1, 3, w, eig) +&
           f(1, 4, w, eig))/4.0_q
  end function

  pure real(q) function J12(w,eig)
    real(q),intent(in) :: w,eig(4)
    J12 = f(2, 1, w, eig)/4.0_q
  end function

  pure real(q) function J13(w,eig)
    real(q),intent(in) :: w,eig(4)
    J13 = f(3, 1, w, eig)/4.0_q
  end function

  pure real(q) function J14(w,eig)
    real(q),intent(in) :: w,eig(4)
    J14 = f(4, 1, w, eig)/4.0_q
  end function

  pure real(q) function vt1(w,eig)
    real(q),intent(in) :: w,eig(4)
    vt1 = f(4, 2, w, eig)*f(3, 2, w, eig)
  end function

  pure real(q) function vt2(w,eig)
    real(q),intent(in) :: w,eig(4)
    vt2 = f(4, 1, w, eig)*f(2, 4, w, eig)*f(3, 2, w, eig)
  end function

  pure real(q) function vt3(w,eig)
    real(q),intent(in) :: w,eig(4)
    vt3 = f(4, 1, w, eig)*f(3, 1, w, eig)*f(2, 3, w, eig)
  end function

  pure real(q) function J21(w,eig)
    real(q),intent(in) :: w,eig(4)
    !J21 = (f(4, 2, w, eig) * f(3, 2, w, eig) +&
    !       f(4, 1, w, eig) * f(2, 4, w, eig) * f(3, 2, w, eig) *&
    !      (1.0_q + f(1, 4, w, eig)) +&
    !       f(4, 1, w, eig) * f(3, 1, w, eig) * f(2, 3, w, eig) *&
    !      (1.0_q + f(1, 4, w, eig) + f(1, 3, w, eig))) / 4.0_q / n2(w, eig)
    J21 = (vt1(w, eig) +&
           vt2(w, eig) * (1.0_q + f(1, 4, w, eig)) +&
           vt3(w, eig) * (1.0_q + f(1, 4, w, eig) + f(1, 3, w, eig))) / 4.0_q / n2(w, eig)
  end function

  pure real(q) function J22(w,eig)
    real(q),intent(in) :: w,eig(4)
    !J22 = (f(4, 2, w, eig) * f(3, 2, w, eig) *&
    !       (1.0 + f(2, 4, w, eig) + f(2, 3, w, eig)) +&
    !       f(4, 1, w, eig) * f(2, 4, w, eig) * f(3, 2, w, eig) *&
    !       (f(2, 4, w, eig) + f(2, 3, w, eig)) +&
    !       f(4, 1, w, eig) * f(3, 1, w, eig) * f(2, 3, w, eig) * f(2, 3, w, eig)) / 4.0_q / n2(w, eig)
    J22 = (vt1(w, eig) * (1.0 + f(2, 4, w, eig) + f(2, 3, w, eig)) +&
           vt2(w, eig) * (f(2, 4, w, eig) + f(2, 3, w, eig)) +&
           vt3(w, eig) * f(2, 3, w, eig)) / 4.0_q / n2(w, eig)
  end function

  pure real(q) function J23(w,eig)
    real(q),intent(in) :: w,eig(4)
    !J23 = (f(4, 2, w, eig) * f(3, 2, w, eig) * f(3, 2, w, eig) +&
    !       f(4, 1, w, eig) * f(2, 4, w, eig) * f(3, 2, w, eig) * f(3, 2, w, eig) +&
    !       f(4, 1, w, eig) * f(3, 1, w, eig) * f(2, 3, w, eig) *&
    !       (f(3, 2, w, eig) + f(3, 1, w, eig))) / 4.0_q / n2(w,eig)
    J23 = (vt1(w, eig) * f(3, 2, w, eig) +&
           vt2(w, eig) * f(3, 2, w, eig) +&
           vt3(w, eig) * (f(3, 2, w, eig) + f(3, 1, w, eig))) / 4.0_q / n2(w,eig)
  end function

  pure real(q) function J24(w,eig)
    real(q),intent(in) :: w,eig(4)
    !J24 = (f(4, 2, w, eig) * f(3, 2, w, eig) *&
    !       f(4, 2, w, eig) +&
    !       f(4, 1, w, eig) * f(2, 4, w, eig) * f(3, 2, w, eig) *&
    !       (f(4, 2, w, eig) +&
    !        f(4, 1, w, eig)) +&
    !       f(4, 1, w, eig) * f(3, 1, w, eig) * f(2, 3, w, eig) *&
    !       f(4, 1, w, eig)) / 4.0_q / n2(w,eig)
    J24 = (vt1(w, eig) * f(4, 2, w, eig) +&
           vt2(w, eig) * (f(4, 2, w, eig) + f(4, 1, w, eig)) +&
           vt3(w, eig) * f(4, 1, w, eig)) / 4.0_q / n2(w,eig)
  end function

  pure real(q) function J31(w,eig)
    real(q),intent(in) :: w,eig(4)
    J31 = (1.0_q -&
            f(1, 4, w, eig) *&
            f(1, 4, w, eig) *&
            f(2, 4, w, eig) *&
            f(3, 4, w, eig)) / 4.0_q / n3(w, eig)
  end function

  pure real(q) function J32(w,eig)
    real(q),intent(in) :: w,eig(4)
    J32 = (1.0_q -&
            f(1, 4, w, eig) *&
            f(2, 4, w, eig) *&
            f(2, 4, w, eig) *&
            f(3, 4, w, eig)) / 4.0_q / n3(w, eig)
  end function

  pure real(q) function J33(w,eig)
    real(q),intent(in) :: w,eig(4)
    J33 = (1.0_q -&
            f(1, 4, w, eig) *&
            f(2, 4, w, eig) *&
            f(3, 4, w, eig) *&
            f(3, 4, w, eig)) / 4.0_q / n3(w, eig)
  end function

  pure real(q) function J34(w,eig)
    real(q),intent(in) :: w,eig(4)
    J34 = (1.0_q -&
            f(1, 4, w, eig) *&
            f(2, 4, w, eig) *&
            f(3, 4, w, eig) *&
            (1.0_q +&
             f(4, 1, w, eig) +&
             f(4, 2, w, eig) +&
             f(4, 3, w, eig))) / 4.0_q / n3(w, eig)
  end function

  !> Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  subroutine tetralite_theta(volume_weight,eig,nw,wvals,alpha,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(in) :: alpha(4) !< weight of the matrix element
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,n

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 tetrahedron completly empty
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        n = volume_weight*n1(w,eig)
        integral(iw) = integral(iw)+n*J11(w,eig)*alpha(1)
        integral(iw) = integral(iw)+n*J12(w,eig)*alpha(2)
        integral(iw) = integral(iw)+n*J13(w,eig)*alpha(3)
        integral(iw) = integral(iw)+n*J14(w,eig)*alpha(4)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        n = volume_weight*n2(w,eig)
        integral(iw) = integral(iw)+n*J21(w,eig)*alpha(1)
        integral(iw) = integral(iw)+n*J22(w,eig)*alpha(2)
        integral(iw) = integral(iw)+n*J23(w,eig)*alpha(3)
        integral(iw) = integral(iw)+n*J24(w,eig)*alpha(4)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        n = volume_weight*n3(w,eig)
        integral(iw) = integral(iw)+n*J31(w,eig)*alpha(1)
        integral(iw) = integral(iw)+n*J32(w,eig)*alpha(2)
        integral(iw) = integral(iw)+n*J33(w,eig)*alpha(3)
        integral(iw) = integral(iw)+n*J34(w,eig)*alpha(4)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) then
        n = volume_weight/4.0
        integral(iw) = integral(iw)+n*alpha(1)
        integral(iw) = integral(iw)+n*alpha(2)
        integral(iw) = integral(iw)+n*alpha(3)
        integral(iw) = integral(iw)+n*alpha(4)
        cycle
      endif

    end do

  end subroutine tetralite_theta

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  subroutine tetralite_delta(volume_weight,eig,nw,wvals,alpha,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(in) :: alpha(4) !< weight of the matrix element
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,g

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 nothing to do
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        g = volume_weight*g1(w,eig)
        integral(iw) = integral(iw)+g*I10(w,eig)*alpha(1)
        integral(iw) = integral(iw)+g*I11(w,eig)*alpha(2)
        integral(iw) = integral(iw)+g*I12(w,eig)*alpha(3)
        integral(iw) = integral(iw)+g*I13(w,eig)*alpha(4)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        g = volume_weight*g2(w,eig)
        integral(iw) = integral(iw)+g*I20(w,eig)*alpha(1)
        integral(iw) = integral(iw)+g*I21(w,eig)*alpha(2)
        integral(iw) = integral(iw)+g*I22(w,eig)*alpha(3)
        integral(iw) = integral(iw)+g*I23(w,eig)*alpha(4)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        g = volume_weight*g3(w,eig)
        integral(iw) = integral(iw)+g*I30(w,eig)*alpha(1)
        integral(iw) = integral(iw)+g*I31(w,eig)*alpha(2)
        integral(iw) = integral(iw)+g*I32(w,eig)*alpha(3)
        integral(iw) = integral(iw)+g*I33(w,eig)*alpha(4)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) exit

    end do

  end subroutine tetralite_delta

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  !> Here only the weight of the first corner is used
  pure subroutine tetralite_delta1(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,g

    do iw=1,nw
      w = wvals(iw)
      ! w < e1 nothing to do
      if (w < eig(1)) cycle
      ! e1 < w < e2
      if (w < eig(2)) then
        g = volume_weight*g1(w,eig)
        integral(iw) = integral(iw)+g*I10(w,eig)
        cycle
      endif
      ! e2 < eps < e3
      if (w < eig(3)) then
        g = volume_weight*g2(w,eig)
        integral(iw) = integral(iw)+g*I20(w,eig)
        cycle
      endif
      ! e3 < eps
      if (w < eig(4)) then
        g = volume_weight*g3(w,eig)
        integral(iw) = integral(iw)+g*I30(w,eig)
        cycle
      endif
      ! e4 < eps
      if (eig(4) < w) exit
    end do
  end subroutine tetralite_delta1

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  !> Here only the weight of the second corner is used
  pure subroutine tetralite_delta2(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,g

    do iw=1,nw
      w = wvals(iw)
      ! w < e1 nothing to do
      if (w < eig(1)) cycle
      ! e1 < w < e2
      if (w < eig(2)) then
        g = volume_weight*g1(w,eig)
        integral(iw) = integral(iw)+g*I11(w,eig)
        cycle
      endif
      ! e2 < eps < e3
      if (w < eig(3)) then
        g = volume_weight*g2(w,eig)
        integral(iw) = integral(iw)+g*I21(w,eig)
        cycle
      endif
      ! e3 < eps
      if (w < eig(4)) then
        g = volume_weight*g3(w,eig)
        integral(iw) = integral(iw)+g*I31(w,eig)
        cycle
      endif
      ! e4 < eps
      if (eig(4) < w) exit
    end do
  end subroutine tetralite_delta2

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  !> Here only the weight of the third corner is used
  pure subroutine tetralite_delta3(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,g

    do iw=1,nw
      w = wvals(iw)
      ! w < e1 nothing to do
      if (w < eig(1)) cycle
      ! e1 < w < e2
      if (w < eig(2)) then
        g = volume_weight*g1(w,eig)
        integral(iw) = integral(iw)+g*I12(w,eig)
        cycle
      endif
      ! e2 < eps < e3
      if (w < eig(3)) then
        g = volume_weight*g2(w,eig)
        integral(iw) = integral(iw)+g*I22(w,eig)
        cycle
      endif
      ! e3 < eps
      if (w < eig(4)) then
        g = volume_weight*g3(w,eig)
        integral(iw) = integral(iw)+g*I32(w,eig)
        cycle
      endif
      ! e4 < eps
      if (eig(4) < w) exit
    end do
  end subroutine tetralite_delta3

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  !> Here only the weight of the fourth corner is used
  pure subroutine tetralite_delta4(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,g

    do iw=1,nw
      w = wvals(iw)
      ! w < e1 nothing to do
      if (w < eig(1)) cycle
      ! e1 < w < e2
      if (w < eig(2)) then
        g = volume_weight*g1(w,eig)
        integral(iw) = integral(iw)+g*I13(w,eig)
        cycle
      endif
      ! e2 < eps < e3
      if (w < eig(3)) then
        g = volume_weight*g2(w,eig)
        integral(iw) = integral(iw)+g*I23(w,eig)
        cycle
      endif
      ! e3 < eps
      if (w < eig(4)) then
        g = volume_weight*g3(w,eig)
        integral(iw) = integral(iw)+g*I33(w,eig)
        cycle
      endif
      ! e4 < eps
      if (eig(4) < w) exit
    end do
  end subroutine tetralite_delta4

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  pure subroutine tetralite_theta1(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,n

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 tetrahedron completly empty
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        n = volume_weight*n1(w,eig)
        integral(iw) = integral(iw)+n*J11(w,eig)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        n = volume_weight*n2(w,eig)
        integral(iw) = integral(iw)+n*J21(w,eig)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        n = volume_weight*n3(w,eig)
        integral(iw) = integral(iw)+n*J31(w,eig)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) then
        n = volume_weight/4.0
        integral(iw) = integral(iw)+n
        cycle
      endif

    end do

  end subroutine tetralite_theta1

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  pure subroutine tetralite_theta2(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,n

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 tetrahedron completly empty
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        n = volume_weight*n1(w,eig)
        integral(iw) = integral(iw)+n*J12(w,eig)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        n = volume_weight*n2(w,eig)
        integral(iw) = integral(iw)+n*J22(w,eig)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        n = volume_weight*n3(w,eig)
        integral(iw) = integral(iw)+n*J32(w,eig)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) then
        n = volume_weight/4.0
        integral(iw) = integral(iw)+n
        cycle
      endif

    end do

  end subroutine tetralite_theta2

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  pure subroutine tetralite_theta3(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,n

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 tetrahedron completly empty
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        n = volume_weight*n1(w,eig)
        integral(iw) = integral(iw)+n*J13(w,eig)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        n = volume_weight*n2(w,eig)
        integral(iw) = integral(iw)+n*J23(w,eig)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        n = volume_weight*n3(w,eig)
        integral(iw) = integral(iw)+n*J33(w,eig)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) then
        n = volume_weight/4.0
        integral(iw) = integral(iw)+n
        cycle
      endif

    end do

  end subroutine tetralite_theta3

  !> @brief Routine to evaluate the delta, the expressions are taken from
  !> A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  !> Journal of Physics C: Solid State Physics 12, 2991 (1979).
  pure subroutine tetralite_theta4(volume_weight,eig,nw,wvals,integral)
    real(q),intent(in) :: volume_weight !< weight associated to one tetrahedron
    real(q),intent(in) :: eig(4) !< eigenvalues at the corners of the tetrahedron (must be sorted)
    integer,intent(in) :: nw !< number of energies
    real(q),intent(in) :: wvals(nw) !< energies at which to evaluate integrated DOS
    real(q),intent(inout) :: integral(nw) !< value of the delta integral for each energy
    ! local variables
    integer :: iw
    real(q) :: w,n

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 tetrahedron completly empty
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        n = volume_weight*n1(w,eig)
        integral(iw) = integral(iw)+n*J14(w,eig)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        n = volume_weight*n2(w,eig)
        integral(iw) = integral(iw)+n*J24(w,eig)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        n = volume_weight*n3(w,eig)
        integral(iw) = integral(iw)+n*J34(w,eig)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) then
        n = volume_weight/4.0
        integral(iw) = integral(iw)+n
        cycle
      endif

    end do

  end subroutine tetralite_theta4

  !> @brief Compute the contibutions from 24 tetrahedra to one k-point
  function weight_24tetra_delta(volume_weight,origin,nw,w,energies_tetra)
      real(q),intent(in) :: volume_weight
      integer,intent(in) :: origin(24)
      integer,intent(in) :: nw
      real(q),intent(in) :: w(nw)
      real(q),intent(in) :: energies_tetra(4,24)
      real(q) :: weight_24tetra_delta(nw)
      ! local variables
      integer :: itet, ic, iw
      integer :: order(4)
      integer :: inv_order(4)
      real(q) :: etetra(4)
      real(q) :: weights(4)
      real(q) :: weights_ordered(4)

      weight_24tetra_delta = 0.0_q
      do itet=1,24
          ! each tetrahedra has 4 energies
          etetra = energies_tetra(:,itet)
          ! sort according to those energies
          order = [1,2,3,4]
          call sort_etetra(etetra,order)
          do ic=1,4
              inv_order(order(ic)) = ic
          enddo
          select case(inv_order(origin(itet)))
              case(1)
                   call tetralite_delta1(volume_weight,etetra,size(w,1),w,weight_24tetra_delta)
              case(2)
                   call tetralite_delta2(volume_weight,etetra,size(w,1),w,weight_24tetra_delta)
              case(3)
                   call tetralite_delta3(volume_weight,etetra,size(w,1),w,weight_24tetra_delta)
              case(4)
                   call tetralite_delta4(volume_weight,etetra,size(w,1),w,weight_24tetra_delta)
          end select
      enddo
  end function weight_24tetra_delta

  !> @brief Compute the contibutions from 24 tetrahedra to one k-point
  function weight_24tetra_theta(volume_weight,origin,nw,w,energies_tetra)
      real(q),intent(in) :: volume_weight
      integer,intent(in) :: origin(24)
      integer,intent(in) :: nw
      real(q),intent(in) :: w(nw)
      real(q),intent(in) :: energies_tetra(4,24)
      real(q) :: weight_24tetra_theta(nw)
      ! local variables
      integer :: itet, ic, iw
      integer :: order(4)
      integer :: inv_order(4)
      real(q) :: etetra(4)
      real(q) :: weights(4)
      real(q) :: weights_ordered(4)

      weight_24tetra_theta = 0.0_q
      do itet=1,24
          ! each tetrahedra has 4 energies
          etetra = energies_tetra(:,itet)
          ! sort according to those energies
          order = [1,2,3,4]
          call sort_etetra(etetra,order)
          do ic=1,4
              inv_order(order(ic)) = ic
          enddo
          select case(inv_order(origin(itet)))
              case(1)
                   call tetralite_theta1(volume_weight,etetra,size(w,1),w,weight_24tetra_theta)
              case(2)
                   call tetralite_theta2(volume_weight,etetra,size(w,1),w,weight_24tetra_theta)
              case(3)
                   call tetralite_theta3(volume_weight,etetra,size(w,1),w,weight_24tetra_theta)
              case(4)
                   call tetralite_theta4(volume_weight,etetra,size(w,1),w,weight_24tetra_theta)
          end select
      enddo
  end function weight_24tetra_theta

  !> @brief sort the tetrahedra according to the indexes representing each corner
  pure subroutine sort_etetra(list,perm)

      integer, intent(inout) :: perm(4)
      real(q), intent(inout) :: list(4)

      integer :: ia,ib,ic,id
      integer :: ilow1,ilow2,ihigh1,ihigh2
      integer :: ilowest,ihighest
      integer :: imiddle1,imiddle2
      real(q) :: va,vb,vc,vd
      real(q) :: vlow1,vlow2,vhigh1,vhigh2
      real(q) :: vlowest,vhighest
      real(q) :: vmiddle1,vmiddle2

      va = list(1); ia = perm(1)
      vb = list(2); ib = perm(2)
      vc = list(3); ic = perm(3)
      vd = list(4); id = perm(4)

      if (va < vb) then
          vlow1 = va; vhigh1 = vb
          ilow1 = ia; ihigh1 = ib
      else
          vlow1 = vb; vhigh1 = va
          ilow1 = ib; ihigh1 = ia
      endif

      if (vc < vd) then
          vlow2 = vc; vhigh2 = vd
          ilow2 = ic; ihigh2 = id
      else
          vlow2 = vd; vhigh2 = vc
          ilow2 = id; ihigh2 = ic
      endif

      if (vlow1 < vlow2) then
          vlowest  = vlow1; vmiddle1 = vlow2
          ilowest  = ilow1; imiddle1 = ilow2
      else
          vlowest  = vlow2; vmiddle1 = vlow1
          ilowest  = ilow2; imiddle1 = ilow1
      endif

      if (vhigh1 > vhigh2) then
          vhighest = vhigh1; vmiddle2 = vhigh2
          ihighest = ihigh1; imiddle2 = ihigh2
      else
          vhighest = vhigh2; vmiddle2 = vhigh1
          ihighest = ihigh2; imiddle2 = ihigh1
      endif

      if (vmiddle1 < vmiddle2) then
          list = [vlowest,vmiddle1,vmiddle2,vhighest]
          perm = [ilowest,imiddle1,imiddle2,ihighest]
      else
          list = [vlowest,vmiddle2,vmiddle1,vhighest]
          perm = [ilowest,imiddle2,imiddle1,ihighest]
      endif

  end subroutine sort_etetra

end module tet_macdonald

!> @brief This module implements an hash table to keep track of unique tetrahedra
module m_tetra_hash
    use prec
    implicit none
    !> A list of tetrahedra. Each tetrahedra is represented by 4 indexes of the corners in the IBZ
    type :: tetra_bucket
        integer :: nstep !< number of elements to add each time the list is increased
        integer :: ntetra !< number of tetrahedra
        integer,allocatable :: tetras(:,:) !< a list of tetrahedra
    end type

    !> An hash table to keep track of unique tetrahedra
    type :: tetra_hash
        integer :: nkibz
        type(tetra_bucket),allocatable :: tetrab(:) !< a list of tetrahedron buckets
    end type
    contains

    !> @brief Initialize a bucket of tetrahedra
    subroutine tetra_bucket_init(self,nstep)
        type(tetra_bucket) :: self
        integer,intent(in) :: nstep
        self%nstep = nstep
        self%ntetra = 0
        allocate(self%tetras(0:4,self%nstep))
    end subroutine tetra_bucket_init

    !> @brief Add an element to the bucket
    subroutine tetra_bucket_add(self,tetra)
        type(tetra_bucket) :: self
        integer,intent(in) :: tetra(4)
        ! local variables
        integer :: nsize
        integer,allocatable :: tmp_tetras(:,:)
        self%ntetra=self%ntetra+1
        nsize = size(self%tetras,2)
        ! if the list is not large enough increase it
        if (nsize<self%ntetra) then
            allocate(tmp_tetras(0:4,nsize+self%nstep))
            tmp_tetras(:,:nsize)=self%tetras
            call move_alloc(tmp_tetras,self%tetras)
        endif
        self%tetras(1:4,self%ntetra)=tetra
        self%tetras(0,self%ntetra)=1
    end subroutine tetra_bucket_add

    !> @brief Find an element in the bucket
    subroutine tetra_bucket_find(self,tetra,idx)
        type(tetra_bucket) :: self
        integer,intent(in) :: tetra(4)
        integer,intent(out) :: idx
        do idx=1,self%ntetra
            if (all(tetra==self%tetras(1:4,idx))) return
        enddo
        idx = -1
    end subroutine

    !> @brief Free the bucket of tetrahedra
    subroutine tetra_bucket_free(self)
        type(tetra_bucket) :: self
        deallocate(self%tetras)
    end subroutine tetra_bucket_free

    !> @brief Initialize the tetrahedron hash table
    subroutine tetra_hash_init(self,nkibz,nstep)
        type(tetra_hash) :: self
        integer,intent(in) :: nkibz
        integer,intent(in) :: nstep
        integer :: ikibz
        self%nkibz = nkibz
        allocate(self%tetrab(self%nkibz))
        do ikibz=1,self%nkibz
            call tetra_bucket_init(self%tetrab(ikibz),nstep)
        enddo
    end subroutine tetra_hash_init

    !> @brief Add a tetrahedron to the hash table or increate the multiplicity of an existing one
    subroutine tetra_hash_add(self,tetra)
        type(tetra_hash) :: self
        integer,intent(in) :: tetra(4)
        ! Local variables
        integer :: ikibz, idx
        !take the first element of the tetrahedron (lowest) and use it to get the bucket
        ikibz = tetra(1)
        call tetra_bucket_find(self%tetrab(ikibz),tetra,idx)
        ! if tetrahedron not found then add it
        if (idx==-1) then
            call tetra_bucket_add(self%tetrab(ikibz),tetra)
        ! if found then increase the multiplicity
        else
            self%tetrab(ikibz)%tetras(0,idx)=self%tetrab(ikibz)%tetras(0,idx)+1
        endif
    end subroutine tetra_hash_add

    !> @brief Return all the tetrahedra on a list
    subroutine tetra_hash_condense(self,ntet,idtet)
        type(tetra_hash) :: self
        integer,intent(out) :: ntet
        integer,allocatable,intent(out) :: idtet(:,:)
        ! local variables
        integer :: itet
        integer :: ikibz
        ! get tetrahedra count
        ntet=0
        do ikibz=1,self%nkibz
            ntet = ntet + self%tetrab(ikibz)%ntetra
        enddo
        allocate(idtet(0:4,ntet))
        ! copy tetrahedrons
        ntet=0
        do ikibz=1,self%nkibz
            do itet=1,self%tetrab(ikibz)%ntetra
                ntet=ntet+1
                idtet(:,ntet)=self%tetrab(ikibz)%tetras(:,itet)
            enddo
        enddo
    end subroutine tetra_hash_condense

    !> @brief Free memory for the hash table
    subroutine tetra_hash_free(self)
        type(tetra_hash) :: self
        integer :: ikibz
        do ikibz=1,self%nkibz
            call tetra_bucket_free(self%tetrab(ikibz))
        end do
        deallocate(self%tetrab)
    end subroutine tetra_hash_free

    !> @brief compute irreducible tetrahera using a hashtable
    subroutine tetirr_hash(kpoints,b,ntet,idtet,iu6)
      use mkpoints_struct_def, only: kpoints_struct
      use mkpoints, only: kint_mod, kint_kindex, sort_tetra, get_tetdisp, kpoints_is_regular
      use tutor, only: vtutor
      type(kpoints_struct),intent(in) :: kpoints
      real(q),intent(in) :: b(3,3) !< reciprocal lattice vectors
      integer,intent(out) :: ntet
      integer,allocatable,intent(out) :: idtet(:,:)
      integer,intent(in) :: iu6 !< unit to write output
      ! local variables
      type(tetra_hash) :: tetrah
      integer,parameter :: nstep=24
      integer :: ix,iy,iz
      integer :: itet,ic,ikfbz
      integer :: kint(3), kint_full(3)
      integer :: tetra_ibz(4), order(4)
      integer :: kcut(3,4,6)

      PROFILING_START('tetirr_hash')

      ! sanity check
      if (.not.kpoints_is_regular(kpoints)) then
         call vtutor%bug('tetirr_hash does not work with non-regular grids',__FILE__,__LINE__)
      endif

      call get_tetdisp(kpoints%b,kcut)

      ! transform kcut
      do itet=1,6
         do ic=1,4
            kcut(:,ic,itet) = matmul(kpoints%ilsnf,kcut(:,ic,itet))
         enddo
      enddo

      call tetra_hash_init(tetrah,kpoints%nkpts,nstep)
      ! loop over all the points in the fbz
      do ix=1,kpoints%isnf(1,1)
        kint_full(1)=ix
        do iy=1,kpoints%isnf(2,2)
          kint_full(2)=iy
          do iz=1,kpoints%isnf(3,3)
            kint_full(3)=iz
            ! loop over all the tetrahedra that split this microzone
            do itet=1,6
              ! loop over the corners of each tetrahedra
              do ic=1,4
                kint = kint_mod(kpoints,kcut(:,ic,itet)+kint_full)
                ikfbz = kint_kindex(kpoints,kint)
                ! store the ibz index point corresponding to this corner
                tetra_ibz(ic) = kpoints%fbz2ibz(1,ikfbz)
              enddo
              ! sort the indexes of corners of the tetrahedra
              call sort_tetra(tetra_ibz,order)
              ! store this tetrahedra or increate multiplicity
              call tetra_hash_add(tetrah,tetra_ibz)
            enddo
          enddo
        enddo
      enddo

      call tetra_hash_condense(tetrah,ntet,idtet)
      call tetra_hash_free(tetrah)
      if (iu6>=0) then
        write(iu6,'(1x,a,i6,a,i8)') 'tetirr_hash: found ',ntet,' inequivalent tetrahedra from ',&
                                     6*kpoints%nkpx*kpoints%nkpy*kpoints%nkpz
      endif
      PROFILING_STOP('tetirr_hash')
    end subroutine tetirr_hash

END MODULE m_tetra_hash
