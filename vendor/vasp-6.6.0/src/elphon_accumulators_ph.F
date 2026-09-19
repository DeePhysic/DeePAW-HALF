#include "symbol.inc"

module elphon_accumulators_ph
    use prec
    use constant
    use bandgap_struct, only: bandgap_info
    use elphon_base, only: elph_selfen_approx_t
    use elphon_accumulators, only: elph_accumulator_selfen_ph
    implicit none

    real(q), parameter :: PHON_TOL=1e-4 !< phonon tolerance to skip the calculation (eV)
    real(q), parameter :: EL_TOL=4e-5 !< electron tolerance to detect degenerate eigenvalues (eV)
    real(q), parameter :: PLANCK_SI = 6.626070040e-34_q
    real(q), parameter :: THZTOEV = 1e12 * PLANCK_SI / EVTOJ / 2 / PI

   contains

    !> @brief Initialize the accumulator object for the self-energy
    subroutine elph_accumulator_selfen_init_ph(self, nbands, nbands_sum, nkpoints, nspin, rspin, nw, enwin, delta, &
                                               ntemps, temps, band_start, band_stop, approx, broad_tol, transport, &
                                               natoms, eigenvalues, phonon_freqs_ibz, nb_totk, band_idx, efermi, carrier_per_cell, nelect, &
                                               fan, dfan, dw, id_idx, id_size, id_name)
        use constant
        use ini, only: register_allocate
        use elphon_base, only: celtot_max
        use elphon_base, only: elph_selfen_approx_t

        type(elph_accumulator_selfen_ph) :: self
        integer, intent(in) :: nw
        integer, intent(in) :: nbands
        integer, intent(in) :: nbands_sum
        integer, intent(in) :: nkpoints
        integer, intent(in) :: nspin
        real(q), intent(in) :: rspin
        integer, intent(in) :: natoms
        real(q), intent(in) :: enwin
        real(q), intent(in) :: delta
        real(q), intent(in) :: efermi(ntemps)
        real(q), intent(in) :: carrier_per_cell(ntemps) !< number of carriers per cell (with doping)
        real(q), intent(in) :: nelect                   !< number of carriers per cell (without doping)
        integer, intent(in) :: ntemps
        integer, intent(in) :: band_start
        integer, intent(in) :: band_stop
        type(elph_selfen_approx_t), intent(in) :: approx !< use static approximation for the electron self-energy
        real(q), intent(in) :: broad_tol
        logical, intent(in) :: transport
        logical, intent(in) :: fan !< store the Fan contribution to the self-energy
        logical, intent(in) :: dfan !< store the derivative of the Fan contribution to the self-energy
        logical, intent(in) :: dw !< store the Debye-Waller contribution to the self-energy
        integer, intent(in) :: id_idx(:) !< integer indentifying this electron self-energy accumulator
        integer, intent(in) :: id_size(:) !< size of the integer indentifying this electron self-energy accumulator
        character(len=*), intent(in) :: id_name(:) !< list of names indentifying this electron self-energy accumulator
        !> temperatures at which the self-energy is to be evaluated (in K)
        real(q), intent(in) :: temps(ntemps)
        complex(q), intent(in) :: eigenvalues(:, :, :) !< eigenvalues of the electronic structure
        real(q), intent(in) :: phonon_freqs_ibz(:, :) !< phonon frequencies in the irreducible BZ
        integer, intent(in) :: nb_totk(:, :)
        integer, intent(in) :: band_idx(:) !< index of the global band that is present locally

        self%nbands = nbands
        self%nbands_sum = nbands; if (nbands_sum > 0 .and. nbands_sum < nbands) self%nbands_sum = nbands_sum
        self%approx = approx
        self%id_idx = id_idx
        self%id_size = id_size
        self%id_name = id_name
        self%broad_tol = broad_tol
        self%transport = transport
        self%tetrahedron = delta <= 0
        self%kcut24_initialized = .false.
        self%nkpoints = nkpoints
        self%nspin = nspin
        self%spindeg = nint(rspin)
        if (nw > 0) then
            self%nw = nw + mod(nw + 1, 2)
            self%energy_grid_mode = 0
        else
            self%energy_grid_mode = 1
            self%nw = max(2, abs(nw))
        end if
        self%emax = maxval(phonon_freqs_ibz) + 5.0_q
        self%emin = minval(phonon_freqs_ibz) - 5.0_q
        self%delta = abs(delta)
        self%enwin = enwin
        self%band_start = 1; if (band_start > 0 .and. band_start < nbands) self%band_start = band_start
        self%band_stop = nbands; if (band_stop > 0 .and. band_stop < nbands) self%band_stop = band_stop
        self%ntemps = ntemps
        allocate (self%tempsk(self%ntemps))
        self%tempsk = temps
        self%natoms = natoms
        self%mode_start = 1
        self%mode_stop = natoms*3
        allocate (self%efermi(self%ntemps))
        self%efermi = efermi
        allocate (self%carrier_per_cell(self%ntemps))
        self%carrier_per_cell = carrier_per_cell
        self%nelect = nelect
        allocate(self%band_idx(size(band_idx)))
        self%band_idx=band_idx

        ! allocate compute_mask array with the information of whether a particular
        ! k-point should be computed
        allocate (self%compute_mask(self%mode_stop, self%nkpoints))
        call register_allocate(VMEM_ARR_SIZE(self%compute_mask), 'elph_accumulator_selfen')
        self%compute_mask = .false.
    end subroutine elph_accumulator_selfen_init_ph

    !> Choose the k-points to compute explicitely
    subroutine elph_accumulator_selfen_selector_iqpt_ph(self, iqpts)
        type(elph_accumulator_selfen_ph) :: self
        integer :: iqpts(:)
        ! local variables
        integer :: iqpt, iqibz, isp
        ! loop over ikpts and set compute_mask
        do iqpt = 1, size(iqpts)
            iqibz = iqpts(iqpt)
            self%compute_mask(:, iqibz) = .true.
        end do
    end subroutine elph_accumulator_selfen_selector_iqpt_ph

    !> Select modes and q-points for which to compute the phonon self-energy
    !> based on the interval from elph_selfen_energy_window(1) up to to elph_selfen_energy_window(2)
    !> On output compute_mask is true for the states for which we want to compute the self-energy
    subroutine elph_accumulator_selfen_selector_energy_window_ph(self, phonon_freqs_ibz, select_energy_window)
        use wave_struct_def, only: wavespin
        type(elph_accumulator_selfen_ph) :: self
        real(q), intent(in) :: phonon_freqs_ibz(:,:)
        real(q), intent(in) :: select_energy_window(2)
        ! local variables
        integer :: iqibz, imode
        real(q) :: omegamq
        self%select_energy_window = select_energy_window

        self%compute_mask = .false.
        ! spin unpolarized and non-collinear spin case
        do iqibz = 1, self%nkpoints
            do imode = self%mode_start, self%mode_stop
                omegamq = phonon_freqs_ibz(imode, iqibz)
                self%compute_mask(imode, iqibz) = &
                    self%compute_mask(imode, iqibz) .or. &
                    (omegamq < self%select_energy_window(2) .and. &
                     omegamq > self%select_energy_window(1))
            end do ! bands
        end do ! kpoints
    end subroutine elph_accumulator_selfen_selector_energy_window_ph

    !> Use information of compute_mask to only store the states for which
    !> the bubble self-energy will be computed
    subroutine elph_accumulator_selfen_sparsify_ph(self, phonon_freqs_ibz)
        use ini, only: register_allocate
        use elphon_accumulators, only: thztoev
        type(elph_accumulator_selfen_ph) :: self
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        ! local variables
        integer :: iqibz, imode, imq
        integer :: nmq, iw, n
        real(q) :: omegamq

        ! set index of each element
        nmq = 0
        allocate (self%mq_idx(self%mode_start:self%mode_stop, self%nkpoints))
        self%mq_idx = 0
        do iqibz = 1, self%nkpoints
            do imode = self%mode_start, self%mode_stop
                if (.not. self%compute_mask(imode, iqibz)) cycle
                nmq = nmq + 1
                self%mq_idx(imode, iqibz) = nmq
            end do
        end do

        ! create energy grids
        self%nmq = nmq
        allocate (self%energies(self%nw, self%nmq))
        do iqibz = 1, self%nkpoints
            do imode = self%mode_start, self%mode_stop
                omegamq = phonon_freqs_ibz(imode, iqibz)*thztoev
                imq = self%mq_idx(imode, iqibz)
                if (imq < 1) cycle
                do iw = 1, self%nw
                    if (self%energy_grid_mode == 0) then
                        self%energies(iw, imq) = &
                            omegamq + ((real(iw, q) - 0.5_q)/self%nw - 0.5_q)*self%enwin
                    else
                        self%energies(iw, imq) = &
                            self%emin + (self%emax - self%emin)/(self%nw - 1)*(iw - 1)
                    end if
                end do ! frequencies
            end do ! bands
        end do ! kpoints

        ! allocate arrays
        allocate (self%selfen_bubble(self%ntemps, self%nw, self%nmq))
        allocate (self%selfen_ph(self%ntemps, self%nw, self%nkpoints, self%natoms*3, self%nspin))
        call register_allocate(VMEM_ARR_SIZE(self%selfen_ph) + &
                               VMEM_ARR_SIZE(self%mq_idx) + &
                               VMEM_ARR_SIZE(self%selfen_bubble),'elph_accumulator_selfen_ph')
        self%selfen_bubble(:,:,:) = cmplx(0.0_q, 0.0_q, q)
        self%selfen_ph(:, :, :, :, :) = cmplx(0.0_q, 0.0_q, q)
    end subroutine elph_accumulator_selfen_sparsify_ph

    !> @brief Write information about the computation of the self-enrgy
    subroutine elph_accumulator_selfen_print_ph(self, naccumulator, iu)
        use string, only: str
        use elphon_base, only: elph_selfen_approx_str
        type(elph_accumulator_selfen_ph), intent(in) :: self
        integer, intent(in) :: naccumulator !< index of the accumulator
        integer, intent(in) :: iu
        ! local variables
        integer :: ikpt, nkpoints_compute, ib, itemp
        logical, allocatable :: compute_band(:)
        if (iu < 0) return
        write (iu, '(A,I3)') 'Phonon self-energy accumulator N =', naccumulator
        write (iu, '(A)') '----------------------------------'
        write (iu, '(A,I8,A,I8,A)') ' Mode range:               [', self%mode_start, ':', self%mode_stop, ']'
        write (iu, '(A,I6)') ' Number of bands to sum over:', self%nbands_sum
        write (iu, '(A,A)') ' Scattering approximation: ', elph_selfen_approx_str(self%approx, verbose=.true.)
        write (iu, '(A,L)') ' Static self-energy: ', self%approx%static
        !write(iu,'(A,F6.3)')      ' Valence band maximum:   ', self%vbm
        !write(iu,'(A,F6.3)')      ' Conduction band minimum:', self%cbm
        write (iu, '(A,F6.3)') ' Complex imaginary shift (delta):', self%delta
        write (iu, '(A)') ' Chemical potentials mu(T):'
        write (iu, '(2A20)') 'Temperature (K)', 'mu (eV)'
        do itemp = 1, self%ntemps
            write (iu, '(2F20.8)') self%tempsk(itemp), self%efermi(itemp)
        end do
        if (self%transport) then
            write (iu, '(A,F8.3,A,F8.3,A,F8.3,A)') ' Transport energy range:   [', self%transport_emin, ':', &
                self%transport_emax, ']  wich corresponds to ', &
                (self%transport_emax - self%transport_emin), ' eV'
        end if
        if (.not. allocated(self%compute_mask)) return
        nkpoints_compute = 0
        do ikpt = 1, self%nkpoints
            if (any(self%compute_mask(:, ikpt))) then
                nkpoints_compute = nkpoints_compute + 1
            end if
        end do
        write (iu, '(A,I6,A,I6,A)') ' Number of selected k-points at which to compute the self-energy: [', &
            nkpoints_compute, ' / ', self%nkpoints, ']'
        write (iu, '(A)') ' Selected bands and k-points at which to compute the self-energy:'

        write (iu, '(A)', advance='no') '      '
        allocate (compute_band(lbound(self%compute_mask, 1):ubound(self%compute_mask, 1)))
        compute_band = .false.
        do ib = lbound(self%compute_mask, 1), ubound(self%compute_mask, 1)
            if (.not. any(self%compute_mask(ib, :))) cycle
            compute_band(ib) = .true.
            write (iu, '(I4)', advance='no') ib
        end do
        write (iu, *)

        do ikpt = 1, self%nkpoints
            if (.not. any(self%compute_mask(:, ikpt))) cycle
            write (iu, '(I6)', advance='no') ikpt
            do ib = lbound(self%compute_mask, 1), ubound(self%compute_mask, 1)
                if (.not. compute_band(ib)) cycle
                write (iu, '(L4)', advance='no') self%compute_mask(ib, ikpt)
            end do
            write (iu, *)
        end do
        write (iu, *)
    end subroutine elph_accumulator_selfen_print_ph

#ifdef VASP_HDF5
    !> @brief Write the information of the phonon self-energy accumulator to a file
    subroutine elph_accumulator_selfen_hdf5write_ph(self, phonon_freqs_ibz, fileid, idx, group, subgroup)
        use vhdf5
        use string, only: str
        implicit none

        type(elph_accumulator_selfen_ph) :: self
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        integer, intent(in) :: idx
        integer(HID_T), intent(in) :: fileid
        character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
        character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
        ! local variables
        integer(HID_T) :: groupid, subgroupid
        character(len=MAX_LEN_GROUP) :: my_group
        character(len=MAX_LEN_GROUP) :: my_subgroup

        my_group = GRP_RESULTS; if (present(group)) my_group = trim(group)
        my_subgroup = GRP_ELPHON; if (present(subgroup)) my_subgroup = trim(subgroup)

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group)//'/'//trim(my_subgroup)//'/phonons', groupid))
        VH5_CHECK(vh5_group_open_or_create(groupid, 'self_energy_'//str(idx), subgroupid))

        ! write important metadata
        VH5_CHECK(vh5_write(subgroupid, 'nbands', self%nbands))
        VH5_CHECK(vh5_write(subgroupid, 'nbands_sum', self%nbands_sum))
        VH5_CHECK(vh5_write(subgroupid, 'delta', self%delta))
        VH5_CHECK(vh5_write(subgroupid, 'band_start', self%band_start))
        VH5_CHECK(vh5_write(subgroupid, 'band_stop', self%band_stop))
        VH5_CHECK(vh5_write(subgroupid, 'nw', self%nw))
        VH5_CHECK(vh5_write(subgroupid, 'enwin', self%enwin))
        VH5_CHECK(vh5_write(subgroupid, 'efermi', self%efermi))
        VH5_CHECK(vh5_write(subgroupid, 'tetrahedron', self%tetrahedron))
        VH5_CHECK(vh5_write(subgroupid, 'select_energy_window', self%select_energy_window))
        VH5_CHECK(vh5_write(subgroupid, 'mq_idx', self%mq_idx))

        ! write main data
        VH5_CHECK(vh5_write(subgroupid, 'selfen_ph', self%selfen_ph))
        VH5_CHECK(vh5_write(subgroupid, 'selfen_bubble', self%selfen_bubble))
        VH5_CHECK(vh5_write(subgroupid, 'phonon_freqs_ibz', phonon_freqs_ibz))
        !VH5_CHECK(vh5_write(subgroupid, "kpoint_generating_vectors", kpoints%B))
        !VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", kpoints%vkpt))
        !VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_weight", kpoints%wtkpt))
        VH5_CHECK(vh5_write(subgroupid, 'temps', self%tempsk))
        VH5_CHECK(vh5_write(subgroupid, 'energies', self%energies))
        !VH5_CHECK(vh5_write(subgroupid, 'frequency_grid', self%w))
        !VH5_CHECK(vh5_write(subgroupid, 'a2F', self%a2F))
        !VH5_CHECK(vh5_write(subgroupid, 'phonon_dos', self%phonon_dos))
        !if (self%nspin == 1 .and. self%spindeg == 2) then   ! non spin polarized
        !    VH5_CHECK(vh5_write(subgroupid, 'dos_at_ef', self%DOSEF*2))
        !else   ! spin polarized and non collinear
        !    VH5_CHECK(vh5_write(subgroupid, 'dos_at_ef', self%DOSEF))
        !end if

        ! close group
        VH5_CHECK(vh5_group_close_writing(subgroupid))
        VH5_CHECK(vh5_group_close_writing(groupid))

    end subroutine elph_accumulator_selfen_hdf5write_ph
#endif

    !> @brief Free the accumulator structure
    subroutine elph_accumulator_selfen_free_ph(self)
        use ini, only: deregister_allocate
        type(elph_accumulator_selfen_ph) :: self
        deallocate (self%tempsk)
        deallocate (self%energies)
        deallocate (self%efermi)
        deallocate (self%carrier_per_cell)
        call deregister_allocate(VMEM_ARR_SIZE(self%selfen_ph) + &
                                 VMEM_ARR_SIZE(self%mq_idx) + &
                                 VMEM_ARR_SIZE(self%compute_mask), 'elph_accumulator_selfen_ph')
        deallocate (self%selfen_ph)
        deallocate (self%compute_mask)
        deallocate (self%mq_idx)
    end subroutine elph_accumulator_selfen_free_ph

    !> @brief Sum results across the nodes
    subroutine elph_accumulator_selfen_sum_ph(self, comm_batch)
        use mpimy, only: communic
        type(elph_accumulator_selfen_ph) :: self
        type(communic) :: comm_batch
        ! local variables
        integer :: itemp, iw, iq, imode, ispin, imq
        PROFILING_START('selfen_sum_ph')
        CALLMPI(m_sum_z(comm_batch, self%selfen_bubble, size(self%selfen_bubble)))
        PROFILING_STOP('selfen_sum_ph')
        ! fill in values in selfen_ph from selfen_bubble
        do imode = 1, self%natoms*3
            do iq = 1, self%nkpoints
                imq = self%mq_idx(imode, iq)
                if (imq < 1) cycle
                do iw = 1, self%nw
                   do itemp = 1, self%ntemps
                       self%selfen_ph(itemp, iw, iq, imode, :) = self%selfen_bubble(itemp, iw, imq)
                   end do
                enddo
            end do
        end do
    end subroutine elph_accumulator_selfen_sum_ph

    subroutine elph_accumulator_selfen_ph_get_lifetimes(self,phonon_lifetimes)
        use tutor, only: vtutor
        type(elph_accumulator_selfen_ph) :: self
        real(q), allocatable, intent(out) :: phonon_lifetimes(:, :, :, :) ! (itemp, iq, imode, ispin)
        ! local variables
        integer :: itemp, iq, imode, ispin, imq
        if (self%nw/=1) then
            call vtutor%error('Phonon self energy computed for more than one frequency (nw) this is not alowed')
        end if
        allocate(phonon_lifetimes(self%ntemps,self%nkpoints,self%natoms*3,self%nspin))
        ! fill in values in selfen_ph from selfen_bubble
        do imode = 1, self%natoms*3
            do iq = 1, self%nkpoints
                imq = self%mq_idx(imode, iq)
                do itemp = 1, self%ntemps
                    phonon_lifetimes(itemp, iq, imode, :) = -2.0_q*aimag(self%selfen_bubble(itemp, 1, imq))
                end do
            end do
        end do
    end subroutine elph_accumulator_selfen_ph_get_lifetimes

    !> find degenerate phonon eigenvalues and average the phonon self-energies
    subroutine elph_accumulator_selfen_ph_degavg(self,phonon_freqs_ibz)
       use wave_struct_def, only: wavespin
       type(elph_accumulator_selfen_ph) :: self
       real(q),intent(in) :: phonon_freqs_ibz(:,:)
       ! local variables
       integer :: ib,iq
       integer :: ibstart,ibstop,ibks
       integer :: ndeg
       complex(q) :: bubble_avg(self%ntemps,self%nw,self%nspin)
       PROFILING_START('selfen_ph_degavg')
       ! find degenerate states
       do iq=1,self%nkpoints
           ib = 1
           do while (ib<self%natoms*3)
               ibstart = ib
               ibstop  = ib
               do while( abs(phonon_freqs_ibz(ibstop+1,iq)-phonon_freqs_ibz(ibstart,iq))<PHON_TOL )
                  ibstop=ibstop+1
                  if (ibstop>=self%natoms*3) exit
               enddo

               ! debug
               !write(*,'(a,i3,a)',advance='no') 'ndeg:',ndeg, ' nbands:'
               !do ib=ibstart,ibstop
               !   write(*,'(i3)',advance='no') ib
               !enddo
               !write(*,*)

               ! sum all the elements
               bubble_avg=0
               ndeg = 0
               do ib=ibstart,ibstop
                   ndeg=ndeg+1
                   bubble_avg(:,:,:) = bubble_avg(:,:,:) + self%selfen_ph(:,:,iq,ib,:)
               enddo
               ! and set them
               do ib=ibstart,ibstop
                   self%selfen_ph(:,:,iq,ib,:) = bubble_avg(:,:,:)/ndeg
               enddo
               ! move to next eigenvalue
               ib=ibstop+1
           enddo
       enddo
       PROFILING_STOP('selfen_ph_degavg')
    end subroutine elph_accumulator_selfen_ph_degavg

    !> Write the results from the calculation of the self-energy to the outcar
    subroutine elph_accumulator_selfen_print_result_ph(self, i, phonon_freqs_ibz, iu)
        use string, only: str
        type(elph_accumulator_selfen_ph) :: self
        integer, intent(in) :: i
        real(q), intent(in) :: phonon_freqs_ibz(:, :)  ! (band, ir-qpts)
        integer, intent(in) :: iu
        ! local variables
        integer :: ikpt, ib, ispin
        integer :: itemp
        real(q) :: reph, imph

        if (iu < 0) return

        WRITE (iu, '(A,I3)') 'Phonon self-energy accumulator N=', i
        if (self%nw > 1) then
            write (iu, *) 'nw > 1 so the self-energy is not written'
            return
        end if

        do itemp = 1, self%ntemps
            write (iu, '(A,F8.0,A)') ' T= ', self%tempsk(itemp), ' K'
            write (iu, '(2A6,3A12)') 'iqpt', 'iband', 'freq THz', 're(ph) eV', 'im(ph) eV'
            do ikpt = 1, self%nkpoints
                do ib = 1, self%natoms*3
                    do ispin = 1, self%nspin
                        reph = real(self%selfen_ph(itemp, 1, ikpt, ib, ispin), q)
                        imph = aimag(self%selfen_ph(itemp, 1, ikpt, ib, ispin))
                        write (iu, '(3I6,3F12.6)') ikpt, ib, ispin, &
                            phonon_freqs_ibz(ib, ikpt), &
                            reph, imph
                    end do
                end do
            end do !kpoints
            write (iu, *)
        end do !itemp
        write (iu, *)
    end subroutine elph_accumulator_selfen_print_result_ph

    !> Check wether this q-point contributes to the imaginary part of the self-energy
    !> of this k-point
    function elph_accumulator_selfen_tet_contributes_low_ph(self, kpoints, triplet, &
                                                            eleig_kp_celen, eleig_k_celtot, eph_thz &
                                                            ) result(contributes)
        use mkpoints, only: kpoints_struct, get_indexes
        use mpimy, only: communic
        use elphon_accumulators, only: get_energies_complex
        use elphon_accumulators, only: tetrahedron_kcontributes
        use triplets, only: triplet_type

        implicit none

        type(elph_accumulator_selfen_ph), intent(in) :: self
        type(kpoints_struct), intent(in) :: kpoints !< kpoints structure
        type(triplet_type), intent(in) :: triplet
        real(q), intent(in) :: eph_thz(:)
        complex(q), intent(in) :: eleig_kp_celen(:, :) !< electron eigenvalues in the IBZ
        complex(q), intent(in) :: eleig_k_celtot(:, :) !< electron eigenvalues for all bands in the IBZ
        logical contributes

        ! local variables
        integer :: imode, ibfbz, ibibz
        integer :: iqibz !< index of the q in the ibz
        integer :: ikibz !< index of the k in the ibz
        integer :: ikpibz !< index of the k' in the ibz
        integer :: idx_k(4, 24)
        integer :: idx_kp(4, 24)
        real(q) :: eleig_tetra_k(4, 24)
        real(q) :: eleig_tetra_kp(4, 24)
        real(q) :: eph(self%natoms*3)

        eph = eph_thz*thztoev
        ikibz = triplet%k%ibz
        ikpibz = triplet%kp%ibz
        iqibz = triplet%q%ibz

        call get_indexes(kpoints, triplet%k%kint, self%kcut24, idx_k)
        call get_indexes(kpoints, triplet%kp%kint, self%kcut24, idx_kp)

        ! loop over bands at k'
        do ibfbz = 1, size(eleig_kp_celen, 1)
            call get_energies_complex(idx_kp, ibfbz, eleig_kp_celen, eleig_tetra_kp)
            ! loop over bands at k
            do ibibz = self%band_start, self%band_stop
                call get_energies_complex(idx_k, ibibz, eleig_k_celtot, eleig_tetra_k)
                if (tetrahedron_kcontributes(eleig_tetra_k - eleig_tetra_kp, eph)) then
                    contributes = .true.
                    return
                end if
            end do
        end do !bands at k'
        contributes = .false.
    end function elph_accumulator_selfen_tet_contributes_low_ph

    !> Accumulate the contribution to phonon self-energy.
    !> This is a high lever wrapper for the lorentzian (not implemented) and tetrahedron method.
    subroutine elph_accumulator_selfen_tet_only_ph(self, kpoints, little_group_orders, ispin, triplet, &
                                                   eleig_kp_celen, eleig_k_celtot, enfbz_idx, eph_thz, band_start_k, g, &
                                                   rspin, comm)
        use mkpoints, only: kpoints_struct
        use lattice, only: latt
        use mpimy, only: communic
        use triplets, only: triplet_type
        implicit none

        type(elph_accumulator_selfen_ph), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints !< kpoints structure
        real(q), intent(in) :: little_group_orders(:)  !< orders or little groups at ir-kpoints
        integer, intent(in) :: ispin !< spin index
        type(triplet_type), intent(in) :: triplet
        integer, intent(in) :: enfbz_idx(:) !< index of the bands
        real(q), intent(in) :: eph_thz(:) !< phonon frequencies in THz with 2pi of an ir-qpoint
        complex(q), intent(in) :: eleig_kp_celen(:, :) !< electron eigenvalues in the IBZ
        complex(q), intent(in) :: eleig_k_celtot(:, :) !< electron eigenvalues for all bands in the IBZ
        integer, intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        complex(q), intent(in) :: g(:, :, :) ! (band at k, band at k', phonon mode)
        real(q), intent(in) :: rspin

        type(communic) :: comm

        if (self%tetrahedron) then
            call elph_accumulator_selfen_tet_ph(self, kpoints, little_group_orders, ispin, triplet, &
                                                eleig_kp_celen, eleig_k_celtot, enfbz_idx, eph_thz, band_start_k, &
                                                g, rspin, comm)
        end if
    end subroutine elph_accumulator_selfen_tet_only_ph

    !> @brief Accumulate the contribution to the phonon self-energy using the tetrahedron method.
    !> Contributions from q-points in FBZ are accumurated to those of corresponding ir-q-points.
    !> This works because only ir-k-points come in this subroutine. The idea is:
    !> (ir-k, kp, q) is mapped to (R.k, R.kp, R.q=ir-q)
    !> implicitly under the invariance of g for rotation.
    !> The phonon lifetime is obtained from the imaginary part of the SE as
    !> hbar/tau_qj= 2 Im SE
    subroutine elph_accumulator_selfen_tet_ph(self, kpoints, little_group_orders, ispin, triplet, &
                                              eleig_kp_celen, eleig_k_celtot, enfbz_idx, eph_thz, &
                                              band_start_k, g, rspin, comm)
        use constant, only: bolkev
        use mkpoints, only: kpoints_struct, get_indexes
        use tet_macdonald, only: weight_24tetra_delta
        use lattice, only: latt
        use tutor, only: vtutor
        use elphon_base, only: fermi_dirac
        use elphon_accumulators, only: phon_tol, thztoev, get_energies_complex
        use mpimy, only: communic
        use triplets, only: triplet_type

        implicit none

        type(elph_accumulator_selfen_ph), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints !< kpoints structure
        real(q), intent(in) :: little_group_orders(:)  !< orders or little groups at ir-kpoints
        integer, intent(in) :: ispin !< spin index
        type(triplet_type), intent(in) :: triplet
        integer, intent(in) :: enfbz_idx(:) !< index of the bands
        real(q), intent(in) :: eph_thz(:)
        complex(q), intent(in) :: eleig_kp_celen(:, :) !< electron eigenvalues in the IBZ
        complex(q), intent(in) :: eleig_k_celtot(:, :) !< electron eigenvalues for all bands in the IBZ
        integer, intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        complex(q), intent(in) :: g(:, :, :) ! (band at k, band at k', phonon mode)
        real(q), intent(in) :: rspin
        type(communic) :: comm !< communicator to paralelize the calculation

        ! local variables
        integer :: imode, iw, ibfbz, ibibz, imq, itemp, irot
        integer :: iqibz !< index of the q in the ibz
        integer :: ikibz !< index of the k in the ibz
        integer :: ikpibz !< index of the k' in the ibz
        integer :: idx_k(4, 24)
        integer :: idx_kp(4, 24)
        real(q) :: eleig_tetra_k(4, 24)
        real(q) :: eleig_tetra_kp(4, 24)
        real(q) :: deltamw(self%nw)
        real(q) :: fk(self%ntemps)
        real(q) :: fkp(self%ntemps)
        real(q) :: temps(self%ntemps)
        real(q) :: eph(self%natoms*3)
        complex(q) :: selfen_ph(self%ntemps, self%nw)

        PROFILING_START('elph_accumulator_selfen_ph_tet')
        eph = eph_thz*thztoev
        temps = self%tempsk*bolkev ! convert to eV

        ikibz = triplet%k%ibz
        ikpibz = triplet%kp%ibz
        iqibz = triplet%q%ibz

        ! TODO: check if any q-point is selected

        call get_indexes(kpoints, triplet%k%kint, self%kcut24, idx_k)
        call get_indexes(kpoints, triplet%kp%kint, self%kcut24, idx_kp)

        ! loop over bands at k'
        do ibfbz = 1, size(eleig_kp_celen, 1)
            ! when ncore/=1 we can still use this comunicator to paralelize the self-energy computation
            if (mod(ibfbz, comm%ncpu) /= comm%node_me - 1) cycle
            if (enfbz_idx(ibfbz) > self%nbands_sum) cycle

            call get_energies_complex(idx_kp, ibfbz, eleig_kp_celen, eleig_tetra_kp)
            do itemp = 1, self%ntemps
                fkp(itemp) = fermi_dirac(real(eleig_kp_celen(ibfbz, ikpibz), q) - self%efermi(itemp), temps(itemp))
            end do

            ! loop over bands at k
            do ibibz = self%band_start, self%band_stop
                call get_energies_complex(idx_k, ibibz, eleig_k_celtot, eleig_tetra_k)
                do itemp = 1, self%ntemps
                    fk(itemp) = fermi_dirac(real(eleig_k_celtot(ibibz, ikibz), q) - self%efermi(itemp), temps(itemp))
                end do

                ! loop over phonon modes
                do imode = 1, self%natoms*3
                    if (abs(eph(imode)) < phon_tol) cycle
                    ! check if this mode is selected
                    imq = self%mq_idx(imode, iqibz)
                    if (imq < 1) cycle

                    ! compute the delta function using the tetrahedron method
                    deltamw = weight_24tetra_delta(kpoints%volwgt, self%kcut24_origin, &
                                                   self%nw, self%energies(:,imq), eleig_tetra_k - eleig_tetra_kp)


                    ! loop over energies (this loop will be removed since nw=1 for phonon selfenergy.)
                    selfen_ph(:, :) = cmplx(0.0_q, 0.0_q, q)
                    do iw = 1, self%nw
                        ! First dim of selfen_ph counts temperatures.
                        selfen_ph(:, iw) = g(ibibz - band_start_k + 1, ibfbz, imode) &
                                           *conjg(g(ibibz - band_start_k + 1, ibfbz, imode)) &
                                           *(fk - fkp) &
                                           *deltamw(iw)*cmplx(0.0_q, pi, q) &
                                           /little_group_orders(triplet%k%ibz)*little_group_orders(triplet%q%ibz)
                    end do

                    self%selfen_bubble(:, :, imq) = self%selfen_bubble(:, :, imq) + selfen_ph(:, :)
                end do ! modes
            end do ! k bands
        end do ! k' bands
        PROFILING_STOP('elph_accumulator_selfen_ph_tet')
    end subroutine elph_accumulator_selfen_tet_ph


end module elphon_accumulators_ph
