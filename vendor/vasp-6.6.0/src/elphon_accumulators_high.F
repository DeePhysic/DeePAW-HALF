#include "symbol.inc"
!> This module implements the accumulators objects which are collections of accumulators
!> these are used to collect electron-phonon quantities such as the electron self-energy,
!> phonon self-energy and electron-phonon matrix elements
!> Written by Henrique Miranda
module elphon_accumulators_high
    use prec
    use constant
    use mpimy, only: communic
    use bandgap_struct, only: bandgap_info
    use triplets, only: triplet_type
    use elphon_accumulators, only: elph_accumulator_selfen, elph_accumulator_selfen_ph, elph_accumulator_g, elph_accumulator_g_hdf5mpiio, &
                                   elph_accumulator_selfen_shallido, elph_accumulator_selfen_fan_contributes, &
                                   elph_accumulator_selfen_fan, elph_accumulator_selfen_dw_g, elph_accumulator_g_hdf5mpiio_free, &
                                   elph_accumulator_g_accumulate, elph_accumulator_g_init, &
                                   elph_accumulator_selfen_dw_t, elph_accumulator_selfen_dw_sum_t, elph_accumulator_selfen_free, &
                                   elph_accumulator_selfen_sum, elph_accumulator_selfen_kk, elph_accumulator_selfen_degavg, &
                                   elph_accumulator_selfen_gaps_compute, elph_accumulator_selfen_gaps_write, &
                                   elph_accumulator_selfen_print_result, elph_accumulator_g_hdf5mpiio_accumulate, elph_accumulator_g_hdf5mpiio_init
    use elphon_accumulators_ph, only: elph_accumulator_selfen_tet_contributes_low_ph, elph_accumulator_selfen_free_ph, &
                                      elph_accumulator_selfen_tet_only_ph, elph_accumulator_selfen_sum_ph, elph_accumulator_selfen_ph_degavg, &
                                      elph_accumulator_selfen_print_result_ph
#ifdef VASP_HDF5
    use elphon_accumulators_ph, only: elph_accumulator_selfen_hdf5write_ph
    use elphon_accumulators, only: elph_accumulator_selfen_hdf5write, elph_accumulator_g_vh5write_free
    use vhdf5
#endif
#ifdef _OPENMP
    USE omp_lib_kinds, ONLY : OMP_LOCK_KIND
#endif
    USE elphon_base, ONLY : elph_selfen_approx_t
    implicit none

    !> contains a list of accumulators for electron self-energy
    type elph_accumulators_selfen
        !> number of accumulators_el
        integer :: naccumulators_el = 0
        !> number of accumulators_ph
        integer :: naccumulators_ph = 0
        !> number of accumulators_mels
        integer :: naccumulators_mels = 0
        !> name of the id of each accumulator
        character(len=15) :: id_name(4)
        !> number of possible elements for each id
        integer :: id_size(4)
        !> electron self-energy accumulators
        type(elph_accumulator_selfen),allocatable :: selfen_el(:)
        !> phonon self-energy accumulators
        type(elph_accumulator_selfen_ph),allocatable :: selfen_ph(:)
    end type elph_accumulators_selfen

    !> contains a list of accumulators for g matrix elements
    type elph_accumulators_g
        !> logical to determine wether to use hdf5 with mpiio support
        logical :: lhdf5_mpiio
        !> electron-phonon matrix element accumulator
        type(elph_accumulator_g) :: g_store
        !> electron-phonon matrix element accumulator to hdf5mpiio
        type(elph_accumulator_g_hdf5mpiio) :: g_hdf5mpiio
    end type elph_accumulators_g

    real(q), parameter :: PHON_TOL=1e-4 !< phonon tolerance to skip the calculation (eV)
    real(q), parameter :: EL_TOL=4e-5 !< electron tolerance to detect degenerate eigenvalues (eV)
    real(q), parameter :: PLANCK_SI = 6.626070040e-34_q
    real(q), parameter :: THZTOEV = 1e12 * PLANCK_SI / EVTOJ / 2 / PI
    contains

    ! from here on its about routines that deal with a set of accumulators
    subroutine elph_accumulators_init(self,elph_set)
        use elphon_base, only: elph_settings, elph_driver_el, elph_driver_ph, elph_driver_mels
        use wave_struct_def, only: wavedes
        type(elph_accumulators_selfen) :: self
        type(elph_settings),intent(in) :: elph_set
        self%naccumulators_el = 0
        self%naccumulators_ph = 0
        self%naccumulators_mels = 0
        select case(elph_set%driver)
           case(elph_driver_el)
             ! electron self-energy
             self%naccumulators_el = elph_set%selfen_naccumulators
             allocate(self%selfen_el(self%naccumulators_el))
           case(elph_driver_ph)
             ! phonon self-energy
             self%naccumulators_ph = elph_set%selfen_naccumulators
             allocate (self%selfen_ph(self%naccumulators_ph))
           case(elph_driver_mels)
             ! matrix elements only, no accumulators needed but should mark all transitions as needed
             ! TODO: actually rename elph_accumulators_selfen to elph_accumulators and move the g accumulators there as well
             self%naccumulators_mels = 1
        end select
    end subroutine elph_accumulators_init

    subroutine elph_accumulators_free(self)
        type(elph_accumulators_selfen) :: self
        ! local variables
        integer :: i
        do i=1,self%naccumulators_el
            call elph_accumulator_selfen_free(self%selfen_el(i))
        enddo
        do i = 1, self%naccumulators_ph
            call elph_accumulator_selfen_free_ph(self%selfen_ph(i))
        end do
    end subroutine elph_accumulators_free

    !> checck if the k, q and kp scattering event contributes to the any of the self-energy accumulators
    function elph_accumulators_contributes(self,elph_set,w_dense_ibz,phonon_freqs_ibz,kpoints,&
                                           triplet,isp) result(contributes)
        use elphon_base, only: elph_settings
        use lattice, only: latt
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct
        use mkpoints, only: kindex_vkpt, ikibz_kindex
        type(elph_accumulators_selfen) :: self
        type(elph_settings) :: elph_set
        type(wavespin) :: w_dense_ibz
        real(q) :: phonon_freqs_ibz(:,:)
        type(kpoints_struct) :: kpoints
        type(triplet_type) :: triplet
        integer :: isp
        logical :: contributes

        PROFILING_START('elph_contributes')
        contributes = .true.
        if (.not.elph_accumulators_shallido(self,triplet%k%ibz,isp)) then
            contributes = .false.
            PROFILING_STOP('elph_contributes')
            return
        endif
        if (elph_accumulators_delta_is_zero(self,elph_set,w_dense_ibz,phonon_freqs_ibz,&
                                            kpoints,triplet,isp)) then
            contributes = .false.
            PROFILING_STOP('elph_contributes')
            return
        endif
        PROFILING_STOP('elph_contributes')
    end function elph_accumulators_contributes

    function elph_accumulators_shallido(self,ik_ibz,isp) result(shallido)
        type(elph_accumulators_selfen) :: self
        integer,intent(in) :: ik_ibz
        integer,intent(in) :: isp
        logical :: shallido
        ! local variables
        integer :: i
        shallido=.false.
        ! first check is wether this k-point is marked for computation
        do i=1,self%naccumulators_el
           if (elph_accumulator_selfen_shallido(self%selfen_el(i),ik_ibz,isp)) then
              shallido=.true.
              return
           endif
        enddo
        do i=1,self%naccumulators_ph
           shallido=.true.
        enddo
        do i=1,self%naccumulators_mels
           shallido=.true.
        enddo
    end function elph_accumulators_shallido

    ! Check if any of the delta functions for the accumulators is not zero
    function elph_accumulators_delta_is_zero(self,elph_set,w_dense_ibz,phonon_freqs_ibz,kpoints_inter,&
                                             triplet,isp) result(delta_is_zero)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct
        type(elph_accumulators_selfen),intent(in) :: self
        type(elph_settings),intent(in) :: elph_set
        type(wavespin),intent(in) :: w_dense_ibz
        real(q),intent(in) :: phonon_freqs_ibz(:,:)
        type(kpoints_struct),intent(in) :: kpoints_inter
        type(triplet_type), intent(in) :: triplet
        integer,intent(in) :: isp
        logical :: delta_is_zero
        ! local variables
        integer :: i
        delta_is_zero = .false.
        ! next check is specific for the case where we are only interested in the imaginary part
        ! of the self-energy
        if (elph_set%selfen_imag_skip) then
           ! check for electron self-energy
           do i=1,self%naccumulators_el
              if (elph_accumulator_selfen_fan_contributes(self%selfen_el(i),kpoints_inter,&
                         triplet,isp,w_dense_ibz%celen(:,:,isp),phonon_freqs_ibz)) then
                 return
              endif
           enddo
           ! check for phonon self-energy
           do i = 1, self%naccumulators_ph
               if (elph_accumulator_selfen_tet_contributes_low_ph( &
                   self%selfen_ph(i), kpoints_inter, &
                   triplet, w_dense_ibz%celen(:, :, isp), w_dense_ibz%celtot(:, :, isp), &
                   phonon_freqs_ibz(:, triplet%q%ibz))) then
                   return
               end if
           end do
           delta_is_zero = .true.
        endif
    end function elph_accumulators_delta_is_zero

    !> Receive the matrix elements and loop over all the accumulators storing what is required
    subroutine elph_accumulators_accumulate(self,elph_set,band_start_k,band_start_kp,triplet,isp,&
                                            kpoints,waverot,w_dense_ibz,wdes_dense_ibz,phonon_freqs_ibz,&
                                            latt_cur,gmode_kkp,gcart_kkp,velocity)
        use lattice, only: latt
        use elphon_base, only: elph_settings
        use mkpoints_struct_def, only: kpoints_struct
        use mkpoints, only: kindex_vkpt, kpoints_fbz_to_ibz
        use wave_struct_def, only: wavespin, wavedes
        use wave_rotate, only: wave_rotator, rotate_vkpt_irot
        type(elph_accumulators_selfen) :: self
        type(elph_settings) :: elph_set
        integer,intent(in) :: band_start_k
        integer,intent(in) :: band_start_kp
        type(triplet_type),intent(in) :: triplet
        integer,intent(in) :: isp
        type(kpoints_struct) :: kpoints
        type(wave_rotator) :: waverot
        type(wavespin) :: w_dense_ibz
        type(wavedes) :: wdes_dense_ibz
        real(q) :: phonon_freqs_ibz(:,:)
        type(latt) :: latt_cur
        complex(q),intent(in) :: gmode_kkp(:,:,:)
        complex(q),intent(in) :: gcart_kkp(:,:,:,:)
        real(q), allocatable, intent(in) :: velocity(:,:,:,:)
        ! local variables
        integer :: i, ikibz, ikpibz, iqibz, irotp
        real(q) :: vk(3), vkp(3), vqpt(3)

        PROFILING_START('elph_accumulators_accumulate')
        ikibz = triplet%k%ibz
        ikpibz= triplet%kp%ibz
        iqibz = triplet%q%ibz

        ikibz = triplet%k%ibz
        vk = kpoints%vkpt(:,ikibz)
        call kpoints_fbz_to_ibz(kpoints,triplet%kp%kindex,ikpibz,irotp)
        vkp = rotate_vkpt_irot(kpoints%vkpt(:,ikpibz),waverot%igrpop,irotp)
        vqpt = vk-vkp

        if (elph_set%selfen_fan) then
           do i=1,self%naccumulators_el
              call elph_accumulator_selfen_fan(self%selfen_el(i),kpoints,latt_cur,&
                           triplet,isp,&
                           w_dense_ibz%celen(:,:,isp),phonon_freqs_ibz,&
                           real(w_dense_ibz%celen(:,ikpibz,isp),q),phonon_freqs_ibz(:,iqibz),&
                           band_start_k,band_start_kp,gmode_kkp,1.0_q/kpoints%nkfull,wdes_dense_ibz%comm_inb,velocity,waverot%igrpop)
           enddo
        endif

        ! check if q=0
        if (ikibz==ikpibz.and.irotp==1) then
            if (elph_set%selfen_dw) then
                do i=1,self%naccumulators_el
                   call elph_accumulator_selfen_dw_g(self%selfen_el(i),ikibz,isp,&
                                                     real(w_dense_ibz%celtot(:,ikpibz,isp),q),&
                                                     real(w_dense_ibz%celen(:,ikpibz,isp),q),band_start_k,band_start_kp,gcart_kkp,wdes_dense_ibz%comm_inb)
                enddo
            endif
        endif
        PROFILING_STOP('elph_accumulators_accumulate')
    end subroutine elph_accumulators_accumulate


    subroutine elph_accumulators_accumulate_ph(self, elph_set, band_start_k, triplet, isp, &
                                               kpoints, little_group_orders, w_dense_ibz, wdes_dense_ibz, &
                                               phonon_freqs_ibz, gmode_kkp)
        use lattice, only: latt
        use elphon_base, only: elph_settings
        use mkpoints_struct_def, only: kpoints_struct
        use wave_struct_def, only: wavespin, wavedes
        use wave_rotate, only: wave_rotator
        use triplets, only: triplet_type
        implicit none

        type(elph_accumulators_selfen), intent(inout) :: self
        type(elph_settings), intent(in) :: elph_set
        integer, intent(in) :: band_start_k
        type(triplet_type), intent(in) :: triplet
        integer, intent(in) :: isp
        type(kpoints_struct), intent(in) :: kpoints
        real(q), intent(in) :: little_group_orders(:)  !< orders or little groups at ir-kpoints
        type(wavespin), intent(in) :: w_dense_ibz
        type(wavedes), intent(in) :: wdes_dense_ibz
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        complex(q), intent(in) :: gmode_kkp(:, :, :)
        ! local variables
        integer :: i

        PROFILING_START('elph_accumulators_accumulate_ph')
        do i = 1, self%naccumulators_ph
            call elph_accumulator_selfen_tet_only_ph(self%selfen_ph(i), kpoints, little_group_orders, isp, triplet, &
                                                     w_dense_ibz%celen(:, :, isp), w_dense_ibz%celtot(:, :, isp), &
                                                     self%selfen_ph(i)%band_idx, phonon_freqs_ibz(:, triplet%q%ibz), &
                                                     band_start_k, gmode_kkp, wdes_dense_ibz%rspin, wdes_dense_ibz%comm_inb)
        end do
        PROFILING_STOP('elph_accumulators_accumulate_ph')
    end subroutine elph_accumulators_accumulate_ph


    !> Prepare the computation of the Debye-Waller term
    subroutine elph_accumulators_prepare_dw(self,elph_set,fbz2ibz,elph_ifc,wdes_dense_ibz,comm_global)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavedes
        use wave_rotate, only: kpoint_star, fbz2ibz_stars, kstar_from_fbz2ibz
        use elphon_potential, only: phonon_ifc, phonon_ifc_get_phonons_zero_disp
        type(elph_accumulators_selfen) :: self
        type(elph_settings) :: elph_set
        type(fbz2ibz_stars) :: fbz2ibz
        type(phonon_ifc) :: elph_ifc
        type(wavedes) :: wdes_dense_ibz
        type(communic) :: comm_global
        ! local variables
        type(kpoint_star) :: kstar
        integer :: ikp_ibz, ikp_star, i
        real(q) :: vkp(3), vk(3), vqpt(3)
        real(q),allocatable :: w_q(:)
        complex(q),allocatable :: e_q(:,:)
        ! Compute phonon modes in the FBZ and accumulate for computation of the Debye-Waller term
        if (elph_set%selfen_dw.or.elph_set%selfen_dws) then
           allocate(w_q(1:3*elph_ifc%natoms_pc))
           allocate(e_q(1:3*elph_ifc%natoms_pc,1:3*elph_ifc%natoms_pc))

           do ikp_ibz=1,wdes_dense_ibz%nkpts
#ifdef MPI
              if (mod(ikp_ibz-1,comm_global%ncpu).ne.comm_global%node_me-1) cycle
#endif
              call kstar_from_fbz2ibz(kstar, wdes_dense_ibz%vkpt(:,ikp_ibz), &
                                                      fbz2ibz, ikp_ibz)
              do ikp_star=1,kstar%nkpts
                 vk = [0.0_q,0.0_q,0.0_q]
                 vkp = kstar%vkpt(:,ikp_star)
                 vqpt=vk-vkp
                 call phonon_ifc_get_phonons_zero_disp(elph_ifc,vqpt,w_q,e_q)
                 do i=1,elph_set%selfen_naccumulators
                    call elph_accumulator_selfen_dw_t(self%selfen_el(i),&
                                                      w_q,e_q,1.0_q/fbz2ibz%nkfbz)
                 enddo
              end do !loop over k' star
           end do !loop over k'
           ! sum t matrix
           do i=1,elph_set%selfen_naccumulators
              call elph_accumulator_selfen_dw_sum_t(self%selfen_el(i),comm_global)
           enddo
        endif
    end subroutine elph_accumulators_prepare_dw

    !> Check if any of the accumulators requires the calculaiton of the velocity
    !> which is sometimes needed for ERTA or MRTA approximations of the scattering
    logical function elph_accumulators_require_velocity(self)
        use elphon_base, only: elph_scattering_approx_mrta_tau, elph_scattering_approx_mrta_lambda, &
                               elph_scattering_approx_erta_tau, elph_scattering_approx_erta_lambda
        type(elph_accumulators_selfen) :: self
        ! local variables
        integer :: i
        elph_accumulators_require_velocity = .false.
        do i=1,self%naccumulators_el
            select case(self%selfen_el(i)%approx%scattering_approx)
                case(elph_scattering_approx_mrta_lambda,&
                     elph_scattering_approx_mrta_tau,&
                     elph_scattering_approx_erta_lambda,&
                     elph_scattering_approx_erta_tau)
                     elph_accumulators_require_velocity = .true.
                    return
                case default
            end select
        enddo
    end function

    !> Finalize calculation of the self-energy
    subroutine elph_accumulators_finalize(self,celtot,phonon_freqs_ibz,comm_global)
        use wave_struct_def, only: wavespin
        type(elph_accumulators_selfen) :: self
        complex(q),intent(in) :: celtot(:,:,:)
        real(q),intent(in) :: phonon_freqs_ibz(:,:)
        type(communic),intent(in) :: comm_global
        ! local variables
        integer :: i
        ! finalize electron self-energy accumulators
        do i=1,self%naccumulators_el
            call elph_accumulator_selfen_sum(self%selfen_el(i),comm_global)
            call elph_accumulator_selfen_kk(self%selfen_el(i))
            call elph_accumulator_selfen_degavg(self%selfen_el(i),celtot)
        enddo
        ! finalize phonon self-energy accumulators
        do i = 1, self%naccumulators_ph
            call elph_accumulator_selfen_sum_ph(self%selfen_ph(i), comm_global)
            call elph_accumulator_selfen_ph_degavg(self%selfen_ph(i), phonon_freqs_ibz)
        end do
    end subroutine elph_accumulators_finalize

#ifdef VASP_HDF5
    !> @brief Write the information of the electron self-energy accumulator to a file
    subroutine elph_accumulators_selfen_hdf5write(self,fileid,group)
        use vhdf5
        use string, only: str
        use elphon_base, only: elph_selfen_approx_hdf5write
        type(elph_accumulators_selfen) :: self
        integer(HID_T),intent(in) :: fileid
        character(len=*), intent(in), optional :: group
        !local variables
        integer(HID_T) :: groupid
        character(len=MAX_LEN_GROUP) :: my_group

        my_group    = GRP_RESULTS//'/'//GRP_ELPHON//'/electrons/self_energy_meta'
        if (present(group))    my_group    = trim(group)

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))

        VH5_CHECK(vh5_write(groupid,'ncalculators',product(self%id_size)))
        VH5_CHECK(vh5_write(groupid,'id_size',self%id_size))
        VH5_CHECK(vh5_write(groupid,'id_name',self%id_name))

        ! close group
        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine elph_accumulators_selfen_hdf5write
#endif

    !> Write results to hdf5 file
    subroutine elph_accumulators_print(self,w_dense_ibz,phonon_freqs_ibz,io)
        use base, only: in_struct
        use wave_struct_def, only: wavespin
#ifdef VASP_HDF5
        USE vhdf5, ONLY: IH5OUTFILEID
#endif
        type(elph_accumulators_selfen) :: self
        type(wavespin),intent(in) :: w_dense_ibz
        real(q),intent(in) :: phonon_freqs_ibz(:,:)
        type(in_struct),intent(in) :: io
        ! local variables
        integer :: i

        !write metadata about electron-phonon self-energy accumulators
#ifdef VASP_HDF5
        call elph_accumulators_selfen_hdf5write(self,ih5outfileid)
#endif

        ! electron self-energy
        do i=1,self%naccumulators_el
            call elph_accumulator_selfen_print_result(self%selfen_el(i),i,w_dense_ibz%celtot,io%iu6)
            call elph_accumulator_selfen_gaps_compute(self%selfen_el(i))
            call elph_accumulator_selfen_gaps_write(self%selfen_el(i),i,io%iu6)
            call elph_accumulator_selfen_gaps_write(self%selfen_el(i),i,io%iu0)
#ifdef VASP_HDF5
            call elph_accumulator_selfen_hdf5write(self%selfen_el(i),ih5outfileid,i)
#endif
        enddo

        ! phonon self-energy
        do i = 1, self%naccumulators_ph
#ifdef VASP_HDF5
            call elph_accumulator_selfen_hdf5write_ph(self%selfen_ph(i), phonon_freqs_ibz, ih5outfileid, i)
#endif
            call elph_accumulator_selfen_print_result_ph(self%selfen_ph(i), i, phonon_freqs_ibz, io%iu6)
        end do

    end subroutine elph_accumulators_print

    !> Initialize accumulators for matrix elements
    subroutine elph_accumulators_g_init(self,band_start_k,band_stop_k,nbands_k,&
                                        band_start_kp,band_stop_kp,nbands_kp,nkpts_k,nkpts_k_local,nkpts_kp,nkpts_kp_local,ispin,natoms_pc,band_idx,lhdf5_mpiio,comm_global)
        type(elph_accumulators_g) :: self
        integer,intent(in) :: band_start_k !< Starting band index for k
        integer,intent(in) :: band_stop_k !< Ending band index for k
        integer,intent(in) :: nbands_k !< Number of bands for k
        integer,intent(in) :: band_start_kp !< Starting band index for k'
        integer,intent(in) :: band_stop_kp !< Ending band index for k'
        integer,intent(in) :: nbands_kp !< Number of bands for k'
        integer,intent(in) :: nkpts_k !< Number of k-points
        integer,intent(in) :: nkpts_k_local !< Number of k-points in the current mpi rank
        integer,intent(in) :: nkpts_kp !< Number of k'-points
        integer,intent(in) :: nkpts_kp_local !< Number of k'-points in the current mpi rank
        integer,intent(in) :: ispin !< Spin channel
        integer,intent(in) :: natoms_pc !< Number of atoms in the primitive cell
        integer,intent(in) :: band_idx(:)
        logical,intent(in) :: lhdf5_mpiio
        type(communic) :: comm_global
        self%lhdf5_mpiio = lhdf5_mpiio
        if (lhdf5_mpiio) then
            ! open an HDF5 file in parallel
            call elph_accumulator_g_hdf5mpiio_init(self%g_hdf5mpiio,nbands_k,band_start_k,nbands_k,band_start_kp,nbands_kp,&
                                      natoms_pc,nkpts_k,&
                                      nkpts_kp,ispin,comm_global)
        else
            ! store electron-phonon matrix elements in memory and then write to file
            call elph_accumulator_g_init(self%g_store,nkpts_k,nkpts_kp,nkpts_k_local,nkpts_kp_local,natoms_pc,ispin,&
                                band_start_k,band_stop_k,nbands_k,band_start_kp,band_stop_kp,nbands_kp,band_idx,comm_global)
        endif
    end subroutine elph_accumulators_g_init

    !> Accumulate matrix elements in the g matrix elements accumlators
    subroutine elph_accumulators_g_accumulate(self,wdes1,ik,ikp_ibz,isp,gmode_kkp,w_q)
        use wave_struct_def, only: wavedes1
        type(elph_accumulators_g) :: self
        type(wavedes1) :: wdes1
        integer,intent(in) :: ik
        integer,intent(in) :: ikp_ibz
        integer,intent(in) :: isp
        complex(q),intent(in) :: gmode_kkp(:,:,:)
        real(q),intent(in) :: w_q(:)
        if (self%lhdf5_mpiio) then
            ! write the matrix elements to hdf5 file
            call elph_accumulator_g_hdf5mpiio_accumulate(self%g_hdf5mpiio, ik, &
                          wdes1, ikp_ibz, isp, gmode_kkp, w_q)
        else
            call elph_accumulator_g_accumulate(self%g_store,ik,ikp_ibz,isp,gmode_kkp,w_q)
        endif
    end subroutine elph_accumulators_g_accumulate

    !> Free electron-phonon g matirx elements accumulators
    subroutine elph_accumulators_g_free(self)
        type(elph_accumulators_g) :: self
#ifdef VASP_HDF5
        if (self%lhdf5_mpiio) then
            call elph_accumulator_g_hdf5mpiio_free(self%g_hdf5mpiio)
        else
            call elph_accumulator_g_vh5write_free(self%g_store)
        endif
#endif
    end subroutine elph_accumulators_g_free

end module
