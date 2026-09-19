#include "symbol.inc"

! Unified driver to compute electron and phonon self-energies due to electron-phonon coupling
! elphon_mels_sandwich routine should have a caching mechanism for the potential
! two consecutive calls with the same q-point skip the interpolation
module elphon_driver
 use prec
 use mkpoints_struct_def, only : kpoints_struct
 use elphon_base, only: elph_settings
 use elphon_mels, only: elph_mels_type
 use elphon_accumulators_high, only: elph_accumulators_selfen, elph_accumulators_g
 implicit none

 contains
 subroutine electron_phonon_driver(ELPH_SET, HAMILTONIAN, KPOINTS, GRID, LATT_CUR, LATT_INI, &
            T_INFO, NONLR_S, NONL_S, W, LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, &
            GRIDC, GRIDUS, C_TO_US, IRDMAX)
__GPU USE moffload_struct_def
    USE base, ONLY: info_struct, in_struct, symmetry
    USE hamil_struct_def, ONLY: ham_handle
    USE tutor, ONLY: vtutor
    USE string, ONLY: str
    USE mpimy, ONLY : communic
    USE mgrid_struct_def, ONLY: grid_3d, transit
    USE wave, ONLY: wdes_from_wdes, wdes_set_npro
    USE lattice, ONLY: latt
    USE nonl_struct_def, ONLY: nonl_struct
    USE nonlr_struct_def, ONLY: nonlr_struct
    USE poscar, ONLY: type_info
    USE pseudo, ONLY: potcar
    USE wave_struct_def, ONLY: wavespin, wavedes
    USE ini, ONLY: dump_allocate
    USE elphon_kgrid, ONLY: elph_kgrid
    USE elphon_base, ONLY: elph_driver_el, elph_driver_ph, elph_driver_mels, elph_header_write, elph_footer_write, &
                           elph_write_settings, elph_scattering_approx_mrta_tau, elph_scattering_approx_mrta_lambda, &
                           elph_scattering_approx_erta_tau, elph_scattering_approx_erta_lambda, elph_scattering_approx_crta, &
                           elph_mels_use_pw, elph_mels_use_wann, elph_read_ifc_pot, scattering_approx_all_crta
    USE elphon_mels, ONLY: elph_mels_init_pw, elph_mels_init_wannier, elph_mels_compute_velocity, elph_mels_ifc_init, elph_mels_g_init, elph_mels_free, elph_mels_barrier
    USE elphon_common, ONLY: elphon_comm_init, elphon_get_kpoints_inter, elphon_prepare_selfenergy_accumulators, elphon_set_nbands, &
                             transport_elphon, elphon_prepare_ikpt
    USE elphon_accumulators_high, ONLY: elph_accumulators_init,elph_accumulators_finalize, elph_accumulators_print, elph_accumulators_free, &
                                   thztoev, elph_accumulators_require_velocity
#ifdef no_elph_superconductivity
#else
    use elphon_superconductivity, only: superconductivity_elphon
#endif
    use electron_transport, only: chemical_potential
    use bandgap_struct, only: bandgap_info
    TYPE (elph_settings)  :: ELPH_SET !< information for the electron-phonon calculation
    TYPE (ham_handle)     :: HAMILTONIAN
    TYPE (kpoints_struct) :: KPOINTS
    TYPE (grid_3d),TARGET :: GRID
    TYPE (latt)           :: LATT_CUR
    TYPE (latt)           :: LATT_INI
    TYPE (type_info)      :: T_INFO
    TYPE (nonlr_struct)   :: NONLR_S
    TYPE (nonl_struct)    :: NONL_S
    TYPE (potcar)         :: P(:)
    TYPE (grid_3d)        :: GRIDC
    TYPE (grid_3d)        :: GRIDUS
    TYPE (transit)        :: C_TO_US
    TYPE (info_struct)    :: INFO
    TYPE (in_struct)      :: IO
    TYPE (symmetry)       :: SYMM
    INTEGER               :: LMDIM, IRDMAX
    TYPE (wavespin), TARGET :: W
    OVERLAP :: CQIJ (:,:,:,:)
    OVERLAP :: CDIJ (:,:,:,:)
    RGRID   :: SV(:,:)
    ! local variables
    type(communic) :: comm_kbatch !< Communicator for each batch of kpoints
    type(communic) :: comm_global !< Communicator for each batch of kpoints
    type(elph_mels_type) :: elph_mels
    type(kpoints_struct) :: kpoints_inter
    type(elph_accumulators_selfen) :: elph_accumulators !< accumulators for electron self-energy
    type(wavedes)  :: wdes_dense_ibz
    type(elph_kgrid) :: kgrid ! kspclib related variables
    logical :: read_ifc
    logical :: qpoints_exists
    real(q), allocatable :: velocity(:,:,:,:) !< electron velocity matrix elements
    type(chemical_potential) :: chempot !< chemical potentials computed at different temperatures and doping levels
    type(bandgap_info),allocatable :: bandgap(:)

#ifdef OFFLOADING
    LOGICAL :: OFFLOAD_ALLOW_SAVE
!$ACC WAIT IF(OFFLOAD_ON)
    OFFLOAD_ALLOW_SAVE = OFFLOAD_ALLOW
    OFFLOAD_ALLOW = .FALSE.
#endif

    call elph_header_write(io)
    call elph_write_settings(elph_set,io%iu6)

    PROFILING_START('elph_driver')
    call wdes_from_wdes(wdes_dense_ibz,w%wdes)
    comm_global = w%wdes%comm
    if (elph_set%ncore/=w%wdes%comm_inb%ncpu) then
       call vtutor%alert('The NCORE='//str(elph_set%ncore)//' setting in the INCAR was not satisfied, try to continue using NCORE='//str(w%wdes%comm_inb%ncpu))
       elph_set%ncore = w%wdes%comm_inb%ncpu
    endif
    call elphon_comm_init(wdes_dense_ibz, elph_set%kpar, elph_set%ncore, comm_kbatch, comm_global, io)
    ! Check if the NCORE setting is compatible with the one in the SCF calculation
    if (wdes_dense_ibz%comm_inb%ncpu/=w%wdes%comm_inb%ncpu) then
       call vtutor%error('The NCORE='//str(w%wdes%comm_inb%ncpu)//' setting in the INCAR cannot be used in the electron-phonon calculation. '//&
                         'You might need to change ELPH_KPAR <= Number of MPI ranks / NCORE.')
    endif

    ! set number of bands based on the input by user or using number of plane-waves
    call elphon_set_nbands(wdes_dense_ibz,elph_set,info%nelect,t_info,maxval(w%wdes%nplwkp_tot))
    call wdes_set_npro(wdes_dense_ibz,t_info,p,info%loverl)

    ! read kpoint_dense file
    call elphon_get_kpoints_inter(kpoints_inter, kgrid, kpoints, elph_set, latt_cur, t_info, wdes_dense_ibz, comm_kbatch, symm%isym, io)

    select case(elph_set%mels_use)
        case(elph_mels_use_pw)
          call elph_mels_init_pw(elph_mels,elph_set,latt_ini,latt_cur,kpoints,kpoints_inter,kgrid,w,wdes_dense_ibz,lmdim,&
                                 hamiltonian,grid,nonl_s,nonlr_s,info,t_info,symm,p,sv,cqij,cdij,chempot,bandgap,io,comm_kbatch,comm_global)
        case(elph_mels_use_wann)
          call elph_mels_init_wannier(elph_mels,elph_set,latt_ini,latt_cur,kpoints,kpoints_inter,w,wdes_dense_ibz,lmdim,irdmax,&
                                 hamiltonian,grid,gridc,gridus,nonl_s,nonlr_s,info,t_info,symm,p,sv,cqij,cdij,chempot,bandgap,io,comm_kbatch,comm_global)
    end select

    ! initialize the force constants
    if (elph_read_ifc_pot(elph_set)) then
        call elph_mels_ifc_init(elph_mels,elph_set,kpoints_inter,qpoints_exists)
        if (qpoints_exists) then
            PROFILING_STOP('elph_driver')
            return
        endif
    endif

    call elphon_prepare_ikpt(elph_set,kpoints_inter,latt_cur,elph_mels%wdes_dense_ibz,io)
    call elph_accumulators_init(elph_accumulators,elph_set)
    ! prepare self-energy accumulators
    call elphon_prepare_selfenergy_accumulators(elph_accumulators, latt_cur, kpoints_inter, elph_set, &
                                                elph_mels%wdes_dense_ibz, &
                                                elph_mels%w_dense_ibz, elph_mels%phonon_freqs_ibz, t_info, info, chempot, bandgap, io)

    if (elph_set%transport.or.elph_accumulators_require_velocity(elph_accumulators)) then
        call elph_mels_compute_velocity(elph_mels,velocity,hamiltonian,elph_set,kpoints_inter,w,latt_ini,latt_cur,elph_mels%wdes_dense_ibz%grid,gridc,gridus,c_to_us,p,lmdim,&
                                        nonlr_s,nonl_s,info,t_info,irdmax,sv,cqij,cdij,symm,io,comm_kbatch)
    endif

    ! if only CRTA is requested we can skip the computation of the electron self-energy altogether
    if (elph_set%transport.and.scattering_approx_all_crta(elph_set)) then
        call transport_elphon(elph_set,elph_accumulators%selfen_el,velocity,latt_cur,kgrid,kpoints_inter,&
                              elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,bandgap,symm,info,io,comm_global)
        PROFILING_STOP('elph_driver')
        return
    endif

    ! initialize the ingredients to compute electron-phonon matrix elements
    call elph_mels_g_init(elph_mels,elph_set,kpoints_inter,w,elph_accumulators,qpoints_exists)
    if (qpoints_exists) then
        PROFILING_STOP('elph_driver')
        return
    endif

    ! report memory usage
    call dump_allocate(io%iu6)
    call wforce(io%iu6)

    ! main loop where electron-phonon matrix elements are computed and accumulated
    select case (elph_set%driver)
        case(elph_driver_el)
            call electron_phonon_selfen_el_driver(elph_set,elph_mels,elph_accumulators,kpoints_inter,velocity)
        case(elph_driver_ph)
            call electron_phonon_selfen_ph_driver(elph_set,elph_mels,elph_accumulators,kpoints_inter)
        case(elph_driver_mels)
            call electron_phonon_mels_driver(elph_set,elph_mels,kpoints_inter)
    end select

    ! mpi barrier + free mpi windows and caches
    call elph_mels_barrier(elph_mels,comm_global)

    ! finalize accumulators
    call elph_accumulators_finalize(elph_accumulators,elph_mels%w_dense_ibz%celtot,elph_mels%phonon_freqs_ibz,comm_global)

#ifdef no_elph_superconductivity
#else
    ! superconductivity stuff
    if (elph_set%superconductivity) then
        call superconductivity_elphon(elph_set, elph_accumulators, elph_mels%phonon_freqs_ibz, elph_mels%w_dense_ibz%celtot, kpoints_inter, latt_cur, io, comm_global)
    end if
#endif

    call elph_accumulators_print(elph_accumulators,elph_mels%w_dense_ibz,elph_mels%phonon_freqs_ibz,io)

    ! transport calculation
    if (elph_set%transport) then
        call transport_elphon(elph_set,elph_accumulators%selfen_el,velocity,latt_cur,kgrid,kpoints_inter,&
                              elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,bandgap,symm,info,io,comm_global)
    endif

    ! free accumulator objects
    call elph_accumulators_free(elph_accumulators)

    call elph_mels_free(elph_mels,kpoints,w,nonl_s)
    PROFILING_STOP('elph_driver')
    call elph_footer_write(io)

#ifdef OFFLOADING
    OFFLOAD_ALLOW = OFFLOAD_ALLOW_SAVE
#endif

 end subroutine electron_phonon_driver

 !> Main loop for the computation of the electron-phonon self-energy
 subroutine electron_phonon_selfen_el_driver(elph_set,elph_mels,elph_accumulators,kpoints,velocity)
    USE elphon_base, ONLY: find_value, elph_progress
    USE elphon_mels, ONLY: elph_mels_compute_g, elph_mels_set_ikp_fbz_list_pointer
    USE elphon_accumulators_high, ONLY: elph_accumulators_accumulate, elph_accumulators_contributes, elph_accumulators_prepare_dw,&
                                   elph_accumulators_shallido, elph_accumulators_delta_is_zero
    USE triplets, ONLY: triplet_type, triplet_init_kp_and_q, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_vkpt,&
                        kpoint_init_kint
    USE mkpoints, ONLY: kindex_from_kp_and_q, ikfbz_kindex, vkpt_kindex, kindex_kint, kint_mod, kint_kindex
    type(elph_settings) :: elph_set
    type(elph_mels_type) :: elph_mels
    type(elph_accumulators_selfen) :: elph_accumulators
    type(kpoints_struct) :: kpoints
    real(q), allocatable, intent(in) :: velocity(:,:,:,:) !< electron velocity matrix elements
    ! local variables
    type(triplet_type) :: triplet
    logical :: contributes
    logical,parameter :: new_loop = .true. !< if true, loop over q points, otherwise loop over k'
    integer,pointer :: ikp_fbz_list(:) !< indexes of the k' points present locally
    integer :: isp
    integer :: iq_fbz,ikp_fbz
    integer :: ik_ibz
    integer :: ik
    integer :: ibz2kindex(kpoints%nkpts)
    complex(q),allocatable :: gcart_kkp(:,:,:,:)
    complex(q),allocatable :: gmode_kkp(:,:,:)

    if (elph_set%dry_run) return

    call elph_mels_set_ikp_fbz_list_pointer(elph_mels,ikp_fbz_list)

    ! create mapping from the ibz to kindex
    do ik_ibz=1,kpoints%nkpts
        ibz2kindex(ik_ibz) = vkpt_kindex(kpoints,elph_mels%latt_cur,kpoints%vkpt(:,ik_ibz))
    enddo

    allocate(gcart_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3,elph_mels%natoms_pc))
    allocate(gmode_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3*elph_mels%natoms_pc))

    ! prepare computation of the debye waller term
    call elph_accumulators_prepare_dw(elph_accumulators,elph_set,elph_mels%fbz2ibz,elph_mels%elph_ifc,elph_mels%wdes_dense_ibz,elph_mels%comm_global)
    PROFILING_START('elph_el_loop')
    ! loop for electron self-energy
    if (new_loop) then
    ! loop over the q points in the fbz
    do iq_fbz=1,kpoints%nkfull
        call elph_progress('fbz qpoint',iq_fbz,kpoints%nkfull,elph_mels%io%iu0)
        call kpoint_init_kindex(triplet%q, kpoints, iq_fbz)

        ! loop over k'
        do ik=1,size(ikp_fbz_list)
            ikp_fbz = ikp_fbz_list(ik)
            call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)

            call kpoint_init_kint(triplet%k, kpoints, triplet%q%kint+triplet%kp%kint)

            ! skip if not ibz
            if (triplet%k%kindex/=ibz2kindex(triplet%k%ibz)) cycle

            do isp=1,elph_mels%ispin
                ! check if delta(ek-ek'+wq) or delta(ek-ek'-wq) function is zero
                contributes = elph_accumulators_contributes(elph_accumulators,elph_set,elph_mels%w_dense_ibz,elph_mels%phonon_freqs_ibz,&
                                                            kpoints,triplet,isp)
                if (.not.contributes) cycle
                ! compute the electron-phonon matrix elements
                call elph_mels_compute_g(elph_mels,elph_set,kpoints,triplet%k%kindex,ikp_fbz,iq_fbz,isp,gcart_kkp,gmode_kkp)
                ! use the electron-phonon matrix elements to compute the electron self-energy
                call elph_accumulators_accumulate(elph_accumulators,elph_set,elph_mels%band_start_k,elph_mels%band_start_kp,triplet,isp,kpoints,elph_mels%waverot,&
                                                  elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,elph_mels%phonon_freqs_ibz,elph_mels%latt_cur,gmode_kkp,gcart_kkp, &
                                                  velocity)
            enddo !loop over spin
        enddo !loop over k'
    enddo !loop over q

    else

    do ik=1,size(ikp_fbz_list)
        call elph_progress('fbz kpoint',ik,size(ikp_fbz_list),elph_mels%io%iu0)
        ikp_fbz = ikp_fbz_list(ik)
        call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)

        do ik_ibz=1,kpoints%nkpts
            if (.not.elph_accumulators_shallido(elph_accumulators,ik_ibz,isp)) cycle

            call kpoint_init_kibz(triplet%k, kpoints, elph_mels%latt_cur, ik_ibz)
            ! compute q=k-kp
            call kpoint_init_kint(triplet%q, kpoints, triplet%k%kint-triplet%kp%kint)

            do isp=1,elph_mels%ispin
                ! check if delta(ek-ek'+wq) or delta(ek-ek'-wq) function is zero
                contributes = elph_accumulators_contributes(elph_accumulators,elph_set,elph_mels%w_dense_ibz,elph_mels%phonon_freqs_ibz,&
                                                            kpoints,triplet,isp)
                if (.not.contributes) cycle
                ! compute the electron-phonon matrix elements
                call elph_mels_compute_g(elph_mels,elph_set,kpoints,triplet%k%kindex,ikp_fbz,triplet%q%kindex,isp,gcart_kkp,gmode_kkp)
                ! use the electron-phonon matrix elements to compute the electron self-energy
                call elph_accumulators_accumulate(elph_accumulators,elph_set,elph_mels%band_start_k,elph_mels%band_start_kp,triplet,isp,kpoints,elph_mels%waverot,&
                                                  elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,elph_mels%phonon_freqs_ibz,elph_mels%latt_cur,gmode_kkp,gcart_kkp, &
                                                  velocity)
            enddo !loop over spin
        enddo !loop over k
    enddo !loop over k'
    endif ! new_loop
    PROFILING_STOP('elph_el_loop')
 end subroutine electron_phonon_selfen_el_driver

 !> Main loop for the computation of the phonon self-energy
 subroutine electron_phonon_selfen_ph_driver(elph_set, elph_mels, elph_accumulators, kpoints)
     use elphon_base, only: elph_settings, elph_progress
     use elphon_mels_pw, only: elph_mels_pw_type
     use elphon_mels, only: elph_mels_type, elph_mels_compute_g
     use elphon_accumulators_high, only: elph_accumulators_contributes, elph_accumulators_accumulate_ph
     use triplets, only: triplet_type, kpoint_init_kindex, kpoint_init_kint
     use mkpoints, only: vkpt_kindex, kpoints_struct
     implicit none

     type(elph_settings), intent(in) :: elph_set
     type(elph_mels_type), target, intent(in) :: elph_mels
     type(elph_accumulators_selfen), intent(inout) :: elph_accumulators
     type(kpoints_struct), intent(in) :: kpoints
     ! local variables
     type(elph_mels_pw_type),pointer :: elph_mels_pw
     type(triplet_type) :: triplet
     logical :: contributes, delta_is_zero
     integer :: isp
     integer :: ikp_fbz, iq_fbz
     integer :: ik_ibz
     integer :: ikp
     real(q), allocatable :: little_group_orders(:)
     integer :: ibz2kindex(kpoints%nkpts)
     complex(q), allocatable :: gcart_kkp(:, :, :, :)
     complex(q), allocatable :: gmode_kkp(:, :, :)

     if (elph_set%dry_run) return
     ! temporary
     elph_mels_pw => elph_mels%pw

     allocate (little_group_orders(kpoints%nkpts))
     do ik_ibz = 1, kpoints%nkpts
         little_group_orders(ik_ibz) = elph_mels_pw%waverot%nrotk/(kpoints%nkfull*kpoints%wtkpt(ik_ibz))
     end do
     !write(*, *) "little_group_orders", little_group_orders

     ! create mapping from the ibz to kindex
     do ik_ibz = 1, kpoints%nkpts
         ibz2kindex(ik_ibz) = vkpt_kindex(kpoints, elph_mels_pw%latt_cur, kpoints%vkpt(:, ik_ibz))
     end do

     allocate (gcart_kkp(elph_mels_pw%nbands_k, elph_mels_pw%nbands, 3, elph_mels_pw%natoms_pc))
     allocate (gmode_kkp(elph_mels_pw%nbands_k, elph_mels_pw%nbands, 3*elph_mels_pw%natoms_pc))

     PROFILING_START('elph_ph_loop')

     ! loop over the q points in the fbz.
     ! Those of (q, k', ir-k) are selected.
     do iq_fbz = 1, kpoints%nkfull
         call elph_progress('fbz qpoint', iq_fbz, kpoints%nkfull, elph_mels_pw%io%iu0)
         call kpoint_init_kindex(triplet%q, kpoints, iq_fbz)

         ! loop over k' in each distribution bin
         do ikp = 1, size(elph_mels_pw%ikp_fbz_list)
             ikp_fbz = elph_mels_pw%ikp_fbz_list(ikp)
             call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)
             call kpoint_init_kint(triplet%k, kpoints, triplet%q%kint + triplet%kp%kint)

             ! Pass through only k in irreducible k-points.
             if (triplet%k%kindex /= ibz2kindex(triplet%k%ibz)) cycle

             do isp=1,elph_mels%ispin
                 contributes = elph_accumulators_contributes( &
                               elph_accumulators, elph_set, elph_mels_pw%w_dense_ibz, elph_mels_pw%phonon_freqs_ibz, &
                               kpoints, triplet, isp)
                 if (.not. contributes) cycle

                 ! compute the electron-phonon matrix elements
                 call elph_mels_compute_g(elph_mels, elph_set, kpoints, triplet%k%kindex, ikp_fbz, iq_fbz, isp, gcart_kkp, gmode_kkp)

                 ! use the electron-phonon matrix elements to compute the electron self-energy
                 call elph_accumulators_accumulate_ph( &
                     elph_accumulators, elph_set, elph_mels_pw%band_start_k, &
                     triplet, isp, kpoints, little_group_orders, &
                     elph_mels_pw%w_dense_ibz, elph_mels_pw%wdes_dense_ibz, &
                     elph_mels_pw%phonon_freqs_ibz, gmode_kkp)
             enddo
         end do !loop over k'
     end do !loop over q
     PROFILING_STOP('elph_ph_loop')

 end subroutine electron_phonon_selfen_ph_driver


 !> Main loop for the computation of electron-phonon matrix elements
 subroutine electron_phonon_mels_driver(elph_set,elph_mels,kpoints)
#ifdef VASP_HDF5
    USE vhdf5, ONLY: vh5_write_lattice_ion
    USE elphon_base, ONLY: vh5_write_elph_kpoints_regular, vh5_write_elph_kpoints_list, vh5_write_elph_eigenvalues_list
#endif
    USE elphon_base, ONLY: find_value, &
                           get_elph_kpoints_eigenvals_list, get_elph_nkpoints_local, &
                           elph_progress, wdes_get_global_band_index
    USE elphon_mels, ONLY: elph_mels_compute_g, elph_mels_compute_g_ibz, elph_mels_set_wavemap_pointer, elph_mels_set_ikp_fbz_list_pointer
    USE elphon_accumulators, ONLY: elph_accumulator_g
    USE elphon_accumulators_high, ONLY: elph_accumulators_g_init, elph_accumulators_g_accumulate, &
                                   elph_accumulators_g_free
    USE triplets, ONLY: triplet_type, triplet_init_kp_and_q, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_vkpt,&
                        kpoint_init_kint
    USE mkpoints_struct_def, ONLY: LineMode, ExplicitList
    USE mkpoints, ONLY: kindex_from_kp_and_q, ikfbz_kindex, vkpt_kindex, kindex_kint, kint_mod, kint_kindex
    USE wave_rotate, ONLY: kpoint_star, wave_rotator_get_star
    USE wave_interpolate, ONLY: wave_mapper
    USE wave_struct_def, ONLY: wavedes1
    USE wave, ONLY: setwdes
    USE vhdf5
    type(elph_settings) :: elph_set
    type(elph_mels_type),target :: elph_mels
    type(kpoints_struct) :: kpoints
    ! local variables
    type(triplet_type) :: triplet
    type(wave_mapper),pointer :: wave_map
    type(elph_accumulators_g) :: elph_accumulators
    type(wavedes1) :: wdes1
    type(kpoint_star) :: star
#ifdef VASP_HDF5
    integer(HID_T) :: fileid
#endif
    logical :: lhdf5_mpiio
    integer,pointer :: ikp_fbz_list(:) !< indexes of the k' points present locally
    integer :: isp, ikp_progress
    integer :: nkpts_k_local, nkpts_kp_local
    integer :: ik_ibz
    integer :: ikp_ibz, ikp_write
    integer :: istar, ik, ikp
    integer :: nkpts_k, nkpts_kp
    integer,allocatable :: kpoint_index(:) ! indexes of the k-points to compute
    integer,allocatable :: band_idx(:)
    real(q),allocatable :: wtkpt_k(:)
    real(q),allocatable :: vkpt_k(:,:)
    real(q),allocatable :: vkpt_kp(:,:)
    real(q),allocatable :: eig_k(:,:,:)
    real(q),allocatable :: eig_kp(:,:,:)
    complex(q),allocatable,target :: gcart_kkp(:,:,:,:)
    complex(q),allocatable,target :: gmode_kkp(:,:,:)
    complex(q),pointer :: g_kkp_write(:,:,:) ! this is just for writing matrix elements to hdf5

    allocate(gcart_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3,elph_mels%natoms_pc))
    allocate(gmode_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3*elph_mels%natoms_pc))

    call elph_mels_set_wavemap_pointer(elph_mels,wave_map)
    call elph_mels_set_ikp_fbz_list_pointer(elph_mels,ikp_fbz_list)

    ! select k-points
    if (allocated(elph_set%selfen_ikpt)) then
        allocate(kpoint_index(size(elph_set%selfen_ikpt)))
        kpoint_index = elph_set%selfen_ikpt
    endif

#ifdef VASP_HDF5
    call get_elph_kpoints_eigenvals_list(kpoints,elph_mels%fbz2ibz,kpoint_index,elph_mels%w_dense_ibz%celtot,&
                                         nkpts_k,nkpts_kp,wtkpt_k,vkpt_k,vkpt_kp,eig_k,eig_kp,lwriteh5)
#endif
    call get_elph_nkpoints_local(wave_map,elph_mels%waverot,elph_mels%fbz2ibz,kpoints,nkpts_k,nkpts_k_local,nkpts_kp_local)
    call wdes_get_global_band_index(elph_mels%wdes_dense_ibz,band_idx)

    lhdf5_mpiio = elph_set%usehdf5mpiio
    call elph_accumulators_g_init(elph_accumulators,elph_mels%band_start_k,elph_mels%band_stop_k,elph_mels%nbands_k,&
                                                    elph_mels%band_start_kp,elph_mels%band_stop_kp,elph_mels%nbands_kp,&
                                                    nkpts_k,nkpts_k_local,nkpts_kp,nkpts_kp_local,elph_mels%ispin,elph_mels%natoms_pc,band_idx,lhdf5_mpiio,elph_mels%comm_global)

    ! need wdes1 to write matrix elements block
    call setwdes(elph_mels%wdes_dense_ibz,wdes1,1)

    PROFILING_START('elph_mels_pw_loop')
    if (kpoints%mode==LineMode.or.kpoints%mode==ExplicitList) then

        ! fix k and k' is k-points

        ! loop over k' points in the local ibz
        do ikp=1,wave_map%nkpts_batch
            call elph_progress("k'-point",ikp,wave_map%nkpts_batch,elph_mels%io%iu0)
            ikp_ibz = wave_map%kpoints_index(ikp+wave_map%nkpts_orig)
            do ik=1,size(kpoint_index)
                ik_ibz = kpoint_index(ik)
                ! loop for computing electron-phonon matrix elements along a path
                do isp=1,elph_mels%ispin
                    ! compute the electron-phonon matrix elements
                    call elph_mels_compute_g_ibz(elph_mels,elph_set,kpoints,ik_ibz,ikp_ibz,isp,gcart_kkp,gmode_kkp)
                    ! accumulate matrix elements
                    call elph_accumulators_g_accumulate(elph_accumulators,wdes1,ik,ikp_ibz,isp,gmode_kkp,elph_mels%w_q)
                enddo ! loop over spin
            enddo !loop over k
        enddo !loop over k'

        ! TODO fix q k is kpoints and k'

    else
        ! loop over k' points in the local ibz
        !do ik=1,size(ikp_fbz_list)
        !    ikp_fbz = ikp_fbz_list(ik)
        !    call kpoint_init_kindex(triplet%kp, kpoints, elph_mels%latt_cur, ikp_fbz)

        ! loop over k' points in the local ibz
        ikp_progress = 0
        do ikp=1,wave_map%nkpts_batch
            ikp_ibz = wave_map%kpoints_index(ikp+wave_map%nkpts_orig)
            ! loop over points in the star of k'
            call wave_rotator_get_star(elph_mels%waverot,kpoints%vkpt(:,ikp_ibz),star)
            do istar=1,star%nkpts
                call kpoint_init_vkpt(triplet%kp, kpoints, elph_mels%latt_cur, star%vkpt(:,istar))
                ikp_write = sum(elph_mels%fbz2ibz%star_size(1:ikp_ibz-1)) + istar
                ikp_progress = ikp_progress + 1
                call elph_progress("k' point",ikp_progress,nkpts_kp_local,elph_mels%io%iu0)

                ! loop over k in the ibz
                do ik=1,size(kpoint_index)
                    ik_ibz = kpoint_index(ik)
                    call kpoint_init_kibz(triplet%k, kpoints, elph_mels%latt_cur, ik_ibz)
                    ! compute q=k-kp
                    call kpoint_init_kint(triplet%q, kpoints, triplet%k%kint-triplet%kp%kint)

                    ! loop for computing electron-phonon matrix elements on a regular k-point grid
                    do isp=1,elph_mels%ispin
                        ! compute the electron-phonon matrix elements
                        call elph_mels_compute_g(elph_mels,elph_set,kpoints,triplet%k%kindex,triplet%kp%kindex,triplet%q%kindex,isp,gcart_kkp,gmode_kkp)

                        ! accumulate matrix elements
                        if (elph_set%mels_cart) then
                            g_kkp_write(1:elph_mels%nbands_k, 1:elph_mels%nbands_local, 1:3*elph_mels%natoms_pc) => gcart_kkp
                        else
                            g_kkp_write => gmode_kkp
                        end if
                        ! write matrix elements in the mode basis
                        call elph_accumulators_g_accumulate(elph_accumulators,wdes1,ik_ibz,ikp_write,isp,g_kkp_write,elph_mels%w_q)
                    enddo ! loop over spin
                enddo !loop over k
            enddo
        enddo !loop over k'
    endif
    PROFILING_STOP('elph_mels_pw_loop')

#ifdef VASP_HDF5
    call elph_accumulators_g_free(elph_accumulators)

    ! write data to HDF5 file
    VH5_CHECK(vh5_file_open_readwrite('vaspelph.h5',fileid))
    call vh5_write_lattice_ion(fileid, elph_mels%info%sznam1, elph_mels%t_info, elph_mels%latt_cur%scale, elph_mels%latt_cur%a,&
                               .false., elph_mels%t_info%posion)
    call vh5_write_elph_kpoints_regular(fileid,elph_mels%waverot%nrotk,elph_mels%waverot%igrpop,&
                                elph_mels%fbz2ibz%irot_fbz2ibz,elph_mels%fbz2ibz%indx_fbz2ibz)
    call vh5_write_elph_kpoints_list(fileid,wtkpt_k,vkpt_k,vkpt_kp)
    if (allocated(vkpt_k))  deallocate(vkpt_k)
    if (allocated(vkpt_kp)) deallocate(vkpt_kp)

    ! write eigenvalues
    call vh5_write_elph_eigenvalues_list(fileid,eig_k,eig_kp)
    if (allocated(eig_k))  deallocate(eig_k)
    if (allocated(eig_kp)) deallocate(eig_kp)
    if (lwriteh5) then
        VH5_CHECK(vh5_file_close(fileid))
    endif
#endif

 end subroutine electron_phonon_mels_driver

end module elphon_driver
