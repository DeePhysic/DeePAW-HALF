#include "symbol.inc"

!> This module is a high-level wrapper for the routines wich compute
!> electron-phonon matrix elements
module elphon_mels
    use prec
    use base, ONLY: info_struct, in_struct, symmetry
    use poscar, ONLY: type_info
    use elphon_base, only: elph_settings, elph_mels_use_pw, elph_mels_use_wann
    use elphon_mels_pw, only: elph_mels_pw_type
    use elphon_mels_wann, only: elph_mels_wann_type
    use elphon_potential_struct, only : phonon_ifc
    use nonl_struct_def, ONLY: nonl_struct
    use nonlr_struct_def, ONLY: nonlr_struct
    use mgrid_struct_def, ONLY: grid_3d, transit
    use mkpoints_struct_def, ONLY: kpoints_struct
    use wave, ONLY: wavespin, wavedes
    use wave_rotate, only: wave_rotator, fbz2ibz_stars
    use mpimy, ONLY: communic
    use lattice, ONLY: latt
    use pseudo, ONLY: potcar
    implicit none

    !> Structure for computing electron-phonon matrix elements using different approaches
    !> Currently support for plane-waves and wannier basis
    type elph_mels_type
        integer :: elph_mels_use !< select which mode to use to compute electron-phonon matrix elements
        type(elph_mels_pw_type) :: pw !< using plane-waves
        type(elph_mels_wann_type) :: wann !< using wannier basis
        ! dimeinsions for the computation of the electron-phonon matrix elements
        integer :: ispin !< number of spins
        integer :: natoms_pc !< number of atoms in the primitive cell
        ! information to compute electron-phonon matrix elements
        integer :: nbands !< total number of bands
        integer :: nbands_local !< number of bands locally (according to the distribution in wavespin)
        integer :: nbands_k !< total number of bands at k
        integer :: band_start_k !< index of the first band at k
        integer :: band_stop_k !< index of the last band at k
        integer :: nbands_kp !< total number of bands at k'
        integer :: band_start_kp !< index of the first band at k'
        integer :: band_stop_kp !< index of the last band at k'
        ! communicators
        type(communic) :: comm_kbatch
        type(communic) :: comm_global
        ! common structures
        type(latt) :: latt_ini !< initial lattice structure
        type(latt) :: latt_cur !< current lattice structure
        type(in_struct) :: io !< input/output structure
        type(type_info) :: t_info !< type information
        type(info_struct) :: info !< information structure
        ! phonons
        integer :: phonon_det_p2s !< vol of supercell/vol primitive cell for phonon calculation
        integer :: phonon_natoms_sc !< number of atoms in the supercell
        real(q),pointer :: phonon_freqs_ibz(:,:) !< store the phonon frequencies in the IBZ
        complex(q),allocatable :: e_q(:,:) !< CACHES phonon eigenvectors at iq_cache (3*natoms,nmodes)
        real(q),allocatable :: w_q(:) !< CACHES phonon frequencies at iq_cache (nmodes)
        type (phonon_ifc) :: elph_ifc
        ! pointers to pw or wan objects
        type(wavespin),pointer :: w_dense_ibz
        type(wavedes),pointer :: wdes_dense_ibz
        type(fbz2ibz_stars),pointer :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
        type(wave_rotator),pointer :: waverot !< dataype to rotate wfs
        !TODO: add supercell wannier structure here
    end type elph_mels_type

    contains
    !> initialize elph_mels structure to compute electron-phonon
    !> matrix elemetns using the plane-wave basis
    subroutine elph_mels_init_pw(self,elph_set,latt_ini,latt_cur,kpoints,kpoints_inter,kgrid,w,wdes_dense_ibz,lmdim,&
                                 hamiltonian,grid,nonl_s,nonlr_s,info,t_info,symm,p,sv,cqij,cdij,chempot,bandgap,io,comm_kbatch,comm_global)
        use elphon_kgrid, ONLY: elph_kgrid
        use elphon_base, ONLY: elph_driver_el, scattering_approx_all_crta, elph_read_ifc_pot
        use elphon_mels_pw, ONLY: elph_mels_pw_read_dimensions, elph_mels_pw_init, elph_mels_pw_get_w_dense, elph_mels_pw_compute_wfs
        use hamil_struct_def, ONLY: ham_handle
        use wave_struct_def, ONLY: wavespin, wavedes
        use electron_transport, only: chemical_potential
        use bandgap_struct, only: bandgap_info
        type(elph_mels_type),target :: self
        type(elph_settings),intent(in) :: elph_set
        type(ham_handle),intent(in) :: hamiltonian
        type(grid_3d),target :: grid
        type(latt),intent(in) :: latt_cur
        type(latt),intent(in) :: latt_ini
        type(kpoints_struct),intent(in) :: kpoints
        type(kpoints_struct) :: kpoints_inter
        type(wavespin),intent(in) :: w
        type(wavedes),intent(in) :: wdes_dense_ibz
        type(elph_kgrid),intent(in) :: kgrid ! kspclib related variables
        type(nonlr_struct),intent(inout) :: nonlr_s
        type(nonl_struct),intent(inout) :: nonl_s
        type(info_struct),intent(in) :: info
        type(in_struct),intent(in) :: io
        type(type_info),intent(in) :: t_info
        type(symmetry),intent(in) :: symm
        integer,intent(in) :: lmdim
        type(potcar),intent(in) :: p(:)
        type(chemical_potential),intent(inout) :: chempot
        type(bandgap_info),allocatable,intent(inout) :: bandgap(:)
        OVERLAP,intent(in) :: cqij(:,:,:,:)
        OVERLAP,intent(in) :: cdij(:,:,:,:)
        RGRID,intent(in) :: sv(:,:)
        type(communic),intent(in) :: comm_kbatch
        type(communic),intent(in) :: comm_global

        self%elph_mels_use = elph_mels_use_pw

        ! initialize elphon object
        call elph_mels_pw_init(self%pw,elph_set,wdes_dense_ibz,lmdim,grid,w%wdes%ngdim,latt_ini,latt_cur,info,p,symm,t_info,io,comm_kbatch,comm_global)

        ! read dimensions from phelel_params.hdf5
        if (elph_read_ifc_pot(elph_set)) then
            call elph_mels_pw_read_dimensions(self%pw,self%phonon_natoms_sc,self%phonon_det_p2s,filename='phelel_params.hdf5')
        endif

        ! prepare WF descriptor in the dense IBZ as well as wave_rotator and fbz2ibz
        call elph_mels_pw_get_w_dense(self%pw,self%pw%wdes_dense_ibz,self%pw%w_dense_ibz,kpoints_inter)
        self%waverot => self%pw%waverot
        self%fbz2ibz => self%pw%fbz2ibz

        ! compute KS states on dense grid this also produces the wave_map
        call elph_mels_pw_compute_wfs(self%pw,kpoints,kpoints_inter,kgrid,w,hamiltonian,nonlr_s,nonl_s,&
                                      cqij,cdij,sv,chempot,bandgap)

        ! set basic dimensions
        self%natoms_pc = self%pw%natoms_pc
        self%ispin = self%pw%ispin
        ! point things
        self%wdes_dense_ibz => self%pw%wdes_dense_ibz
        self%w_dense_ibz => self%pw%w_dense_ibz
        ! copy things
        self%latt_cur = latt_cur
        self%latt_ini = latt_ini
        self%comm_global = comm_global
        self%comm_kbatch = comm_kbatch
        self%io = io
        self%t_info = t_info
        self%info = info
    end subroutine elph_mels_init_pw

    !> initialize elph_mels structure to compute electron-phonon
    !> matrix elements using the wannier basis
    subroutine elph_mels_init_wannier(self,elph_set,latt_ini,latt_cur,kpoints,kpoints_inter,w,wdes_dense_ibz,lmdim,irdmax,&
                                      hamiltonian,grid,gridc,gridus,nonl_s,nonlr_s,info,t_info,symm,&
                                      p,sv,cqij,cdij,chempot,bandgap,io,comm_kbatch,comm_global)
        use hamil_struct_def, ONLY: ham_handle
        use elphon_mels_wann, only: elph_mels_wann_init, elph_mels_wann_get_w_dense, &
                                    elph_mels_wann_compute_wfs
        use wannier_interpolation, only: wannier_interpol_setup_kpoints, wannier_interpol_from_wanproj, wannier_interpol_from_mlwf, &
                                         wannier_interpol_setup_ham, wannier_interpol_setup_cder, wannier_interpolate, &
                                         wannier_interpol_free
        use electron_transport, only: chemical_potential
        use elphon_common, only: elphon_run_precalculation, elphon_run_chempot
        use elphon_base, ONLY: elph_driver_el, scattering_approx_all_crta, elph_read_ifc_pot
        use elphon_mels_pw, ONLY: elph_mels_pw_read_dimensions
        use wap, ONLY: read_phelel_phonon_dimensions_hdf5
        use bandgap_struct, only: bandgap_info
        use full_kpoints, only: lfull_kpoints
        use mlwf, only: mlwf_global
        type(elph_mels_type),target :: self
        type(elph_settings),intent(in) :: elph_set
        type(ham_handle),intent(in) :: hamiltonian
        type(grid_3d),intent(in) :: grid
        type(grid_3d),intent(in) :: gridc
        type(grid_3d),intent(in) :: gridus
        type(latt),intent(in) :: latt_cur
        type(latt),intent(in) :: latt_ini
        type(kpoints_struct),intent(in) :: kpoints
        type(kpoints_struct) :: kpoints_inter
        type(wavespin),intent(in) :: w
        type(wavedes) :: wdes_dense_ibz
        type(nonlr_struct),intent(inout) :: nonlr_s
        type(nonl_struct),intent(inout) :: nonl_s
        type(info_struct),intent(in) :: info
        type(in_struct),intent(in) :: io
        type(type_info),intent(in) :: t_info
        type(symmetry),intent(in) :: symm
        integer,intent(in) :: lmdim
        integer,intent(in) :: irdmax
        type(potcar),intent(in) :: p(:)
        type(chemical_potential),intent(inout) :: chempot
        type(bandgap_info),allocatable,intent(inout) :: bandgap(:)
        OVERLAP,intent(in) :: cqij(:,:,:,:)
        OVERLAP,intent(in) :: cdij(:,:,:,:)
        RGRID,intent(in) :: sv(:,:)
        type(communic),intent(in) :: comm_kbatch
        type(communic),intent(in) :: comm_global
        !local variables
        integer :: i
        logical :: lfull_kpoints_save

        self%elph_mels_use = elph_mels_use_wann

        ! read dimensions from phelel_params.hdf5
        if (elph_read_ifc_pot(elph_set)) then
            call read_phelel_phonon_dimensions_hdf5(self%natoms_pc,self%phonon_natoms_sc,self%phonon_det_p2s,filename='phelel_params.hdf5')
        endif

        ! initialize wann mels structure
        call elph_mels_wann_init(self%wann,wdes_dense_ibz,grid,latt_ini,latt_cur,symm,w%wdes%ngdim)

        ! initialize wannier_interpol struct
        call wannier_interpol_setup_kpoints(self%wann%wann_interp,kpoints,latt_cur,symm,w)
        !call wannier_interpol_from_wanproj(self%wann%wann_interp,io)
        call wannier_interpol_from_mlwf(self%wann%wann_interp,mlwf_global,io)
        call wannier_interpol_setup_ham(self%wann%wann_interp,w,io)

        !shenenigans
        lfull_kpoints_save = lfull_kpoints
        lfull_kpoints = .false.
        ! create wdes_dense_ibz and w_dense_ibz as well as wave_rotator and fbz2ibz
        call elph_mels_wann_get_w_dense(self%wann,self%wann%wdes_dense_ibz,self%wann%w_dense_ibz,kpoints_inter)
        lfull_kpoints = lfull_kpoints_save
        self%waverot => self%wann%waverot
        self%fbz2ibz => self%wann%fbz2ibz

        call elph_mels_wann_compute_wfs(self%wann,kpoints_inter,symm,latt_ini,latt_cur,t_info,info,io,comm_kbatch)

        ! prepare WF descriptor in the dense IBZ
        wdes_dense_ibz%nkpts = kpoints_inter%nkpts

        ! compute chemical potential and bandgap
        call elphon_run_chempot(self%wann%w_dense_ibz, chempot, bandgap, &
                                latt_cur, latt_ini, t_info, &
                                symm, &
                                info, elph_set, self%wann%wdes_dense_ibz, comm_kbatch, comm_global, &
                                kpoints_inter, io)

        ! set basic dimensions
        self%natoms_pc = sum(t_info%nityp)
        self%ispin = self%wann%ispin
        ! point things
        self%wdes_dense_ibz => self%wann%wdes_dense_ibz
        self%w_dense_ibz => self%wann%w_dense_ibz
        ! copy things
        self%latt_cur = latt_cur
        self%latt_ini = latt_ini
        self%comm_global = comm_global
        self%comm_kbatch = comm_kbatch
        self%io = io
        self%t_info = t_info
        self%info = info
    end subroutine elph_mels_init_wannier

    !> compute velocity matrix elements using plane-waves or wannier basis
    subroutine elph_mels_compute_velocity(self,velocity,hamiltonian,elph_set,kpoints_inter,w,latt_ini,latt_cur,grid,gridc,gridus,c_to_us,p,lmdim,&
                                          nonlr_s,nonl_s,info,t_info,irdmax,sv,cqij,cdij,symm,io,comm_kbatch)
        use base, ONLY: info_struct, in_struct, symmetry
        use hamil_struct_def, ONLY: ham_handle
        use lattice, ONLY: latt
        use pseudo, ONLY: potcar
        use elphon_common, ONLY: elphon_get_velocity, elphon_write_velocity
        use elphon_mels_wann, ONLY: elph_mels_wann_get_velocity
        use wannier_interpolation, only: wannier_interpol_setup_cder, wannier_interpolate_cder
        type(elph_mels_type) :: self
        real(q),allocatable,intent(out) :: velocity(:,:,:,:)
        type(ham_handle),intent(in) :: hamiltonian !< handle to the hamiltonian
        type(elph_settings),intent(in) :: elph_set
        type(kpoints_struct),intent(in) :: kpoints_inter
        type(wavespin),intent(in) :: w
        type(latt),intent(in) :: latt_cur
        type(latt),intent(in) :: latt_ini
        type(grid_3d),target :: grid
        type(grid_3d),target :: gridc
        type(grid_3d),target :: gridus
        type(transit) :: c_to_us
        type(nonlr_struct),intent(inout) :: nonlr_s
        type(nonl_struct),intent(inout) :: nonl_s
        integer,intent(in) :: lmdim, irdmax
        type(potcar),intent(in) :: p(:)
        OVERLAP,intent(in) :: cqij(:,:,:,:)
        OVERLAP,intent(in) :: cdij(:,:,:,:)
        RGRID,intent(in) :: sv(:,:)
        type(info_struct),intent(in) :: info
        type(type_info),intent(in) :: t_info
        type(symmetry),intent(in) :: symm
        type(in_struct),intent(in) :: io
        type(communic),intent(in) :: comm_kbatch
        !local variables

        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                call elphon_get_velocity(velocity, hamiltonian, p, lmdim, self%pw%nbands, self%pw%wf_ncdij, latt_cur, latt_ini, symm, nonlr_s, nonl_s, &
                         t_info, elph_set, self%pw%wdes_dense_ibz, self%pw%w_dense_ibz, &
                         self%pw%wave_map, self%pw%wdes_kbatch, self%pw%w_kbatch, sv, cqij, cdij, &
                         grid, gridc, gridus, c_to_us, irdmax, info, io, comm_kbatch)
                call elphon_write_velocity(elph_set, velocity, self%pw%w_dense_ibz%celtot, io)
                !TODO: this should go in some other place?
                deallocate(self%pw%w_kbatch%CPROJ)
            case(elph_mels_use_wann)
                call elph_mels_wann_get_velocity(self%wann,hamiltonian,w,grid,gridc,gridus,c_to_us,irdmax,lmdim,t_info,&
                                          nonlr_s,nonl_s,p,sv,cqij,cdij,latt_ini,latt_cur,symm,velocity,info,io,comm_kbatch)
                call elphon_write_velocity(elph_set, velocity, self%wann%w_dense_ibz%celtot, io)
        end select
    end subroutine elph_mels_compute_velocity

    !> initialize the elphon_mels structure to compute the phonon dispersion from
    !> the interatomic force constants
    subroutine elph_mels_ifc_init(self,elph_set,kpoints_inter,qpoints_exists)
        use elphon_potential, ONLY: phonon_ifc_init, phonon_ifc_from_hdf5
        use elphon_common, ONLY: elphon_prepare_phonon
        type(elph_mels_type),target :: self
        type(elph_settings),intent(in) :: elph_set !< elphon settings
        type(kpoints_struct),intent(in) :: kpoints_inter
        logical :: qpoints_exists
        ! local variables
        logical :: make_hermitian

        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                make_hermitian = elph_set%ifc_asr /= -2
                call phonon_ifc_init(self%elph_ifc,self%phonon_det_p2s,self%natoms_pc,self%phonon_natoms_sc,make_hermitian)
                call phonon_ifc_from_hdf5(self%elph_ifc,filename='phelel_params.hdf5')
                call elphon_prepare_phonon(self%elph_ifc, elph_set, self%latt_cur, self%pw%grid, self%pw%t_info, self%pw%io, self%pw%comm_global)
                call elph_mels_ifc_compute_dispersion(self,qpoints_exists,filename='QPOINTS')
                call elph_mels_ifc_compute_ibz_phonons(self,elph_set,self%pw%phonon_freqs_ibz,kpoints_inter)
                self%phonon_freqs_ibz => self%pw%phonon_freqs_ibz
            case(elph_mels_use_wann)
                make_hermitian = elph_set%ifc_asr /= -2
                call phonon_ifc_init(self%elph_ifc,self%phonon_det_p2s,self%natoms_pc,self%phonon_natoms_sc,make_hermitian)
                call phonon_ifc_from_hdf5(self%elph_ifc,filename='phelel_params.hdf5')
                call elphon_prepare_phonon(self%elph_ifc, elph_set, self%latt_cur, self%wann%grid, self%t_info, self%io, self%comm_global)
                call elph_mels_ifc_compute_dispersion(self,qpoints_exists,filename='QPOINTS')
                call elph_mels_ifc_compute_ibz_phonons(self,elph_set,self%wann%phonon_freqs_ibz,kpoints_inter)
                self%phonon_freqs_ibz => self%wann%phonon_freqs_ibz
        end select
        allocate(self%w_q(3*self%natoms_pc))
        allocate(self%e_q(3*self%natoms_pc,3*self%natoms_pc))
    end subroutine elph_mels_ifc_init

    !> compute phonon dispersion on a path specified in the QPOINTS file
    subroutine elph_mels_ifc_compute_dispersion(self,qpoints_exists,filename)
        USE mkpoints, ONLY: rd_kpoints
        USE elphon_potential, ONLY: phonon_ifc_dispersion
        type(elph_mels_type) :: self
        logical :: qpoints_exists
        character(len=*),intent(in) :: filename
        !local variables
        type(kpoints_struct) :: qpoints
        ! if compute phonon dispersion
        inquire(file=filename, exist=qpoints_exists)
        if (qpoints_exists) then
            qpoints%ltet = .false.
            call rd_kpoints(qpoints, self%latt_cur, .false., .false., self%io%iu6, self%io%iu0, file_name = filename)

            ! we use the d convention for the dynamical matrix of phonon: exp(iq.l)
            call phonon_ifc_dispersion(self%elph_ifc,qpoints%nkpts,qpoints%vkpt)
        endif
    end subroutine elph_mels_ifc_compute_dispersion

    !> Compute phonon frequencies in the IBZ
    subroutine elph_mels_ifc_compute_ibz_phonons(self,elph_set,phonon_freqs_ibz,kpoints_inter)
        use ini, only: register_allocate
        use elphon_base, only: elph_settings
        use elphon_common, only: elphon_compute_ibz_phonons
        type(elph_mels_type) :: self
        type(elph_settings),intent(in) :: ELPH_SET
        real(q),allocatable,intent(out) :: phonon_freqs_ibz(:,:) !< phonon frequencies in the IBZ
        type(kpoints_struct) :: kpoints_inter !< kpoints in the IBZ
        ! allocate phonons in the IBZ
        allocate(phonon_freqs_ibz(self%elph_ifc%natoms_pc*3,kpoints_inter%nkpts))
        call register_allocate(VMEM_ARR_SIZE(phonon_freqs_ibz),"elph_phonon_freqs")
        ! compute phonon frequencies in the IBZ
        call elphon_compute_ibz_phonons(self%elph_ifc,elph_set,kpoints_inter,phonon_freqs_ibz)
    end subroutine elph_mels_ifc_compute_ibz_phonons

    !> Initialize the elphon_mels structure to compute the g matrix elements
    !> using either plane-waves or wannier basis
    subroutine elph_mels_g_init(self,elph_set,kpoints_inter,w,elph_accumulators,qpoints_exists)
        use elphon_base, only: elph_driver_el, elph_driver_ph, elph_driver_mels
        use elphon_mels_pw, ONLY: elph_mels_pw_pot_init, elph_mels_pw_pot_read_hdf5, elph_mels_pw_pot_setup, &
                                  elph_mels_pw_pot_dispersion, elph_mels_pw_set_nbandsk, elph_mels_pw_set_nbandskp, elph_mels_pw_wf_write, &
                                  elph_mels_pw_cache_init, elph_mels_pw_windows_init, elph_mels_pw_wave_caches_init, &
                                  elph_mels_pw_wf_redistribute, elph_mels_pw_cache_prefill, &
                                  wave_map_get_ikpfbz
        use elphon_accumulators_high, only: elph_accumulators_selfen
        type(elph_mels_type) :: self
        type(elph_settings),intent(in) :: elph_set
        type(kpoints_struct),intent(in) :: kpoints_inter
        type(wavespin),intent(in) :: w
        type(elph_accumulators_selfen),intent(in) :: elph_accumulators
        logical,intent(out) :: qpoints_exists !< whether QPOINTS_POT file exists
        ! local variables
        integer :: num_wann !< number of wannier functions
        integer :: start_k, stop_k

        qpoints_exists = .false.
        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                ! initialize the electron-phonon potential from phelel_params.hdf5
                call elph_mels_pw_pot_init(self%pw,w%wdes%COMM_intra_node)
                call elph_mels_pw_pot_read_hdf5(self%pw,filename='phelel_params.hdf5')
                call elph_mels_pw_pot_setup(self%pw)
                call elph_mels_pw_pot_dispersion(self%pw,qpoints_exists,filename='QPOINTS_POT')
                if (qpoints_exists) return

                ! set the number of bands to be used at k and k' depending on the min and max values among all the accumulators
                select case (elph_set%driver)
                    case (elph_driver_ph)
                      call elph_mels_pw_set_nbandsk(self%pw,minval(elph_accumulators%selfen_ph(:)%band_start),maxval(elph_accumulators%selfen_ph(:)%band_stop))
                      call elph_mels_pw_set_nbandskp(self%pw,1,self%pw%nbands)
                    case (elph_driver_el)
                      call elph_mels_pw_set_nbandsk(self%pw,minval(elph_accumulators%selfen_el(:)%band_start),maxval(elph_accumulators%selfen_el(:)%band_stop))
                      call elph_mels_pw_set_nbandskp(self%pw,minval(elph_accumulators%selfen_el(:)%band_start_kp),maxval(elph_accumulators%selfen_el(:)%band_stop_kp))
                      !Uncomment the line bellow to use all the bands
                      !call elph_mels_pw_set_nbandskp(self%pw,1,self%pw%nbands)
                    case (elph_driver_mels)
                        start_k=1
                        if (elph_set%selfen_band_start>0.and.elph_set%selfen_band_start<self%pw%nbands) start_k=elph_set%selfen_band_start
                        stop_k=self%pw%nbands
                        if (elph_set%selfen_band_stop>0.and.elph_set%selfen_band_start<self%pw%nbands) stop_k=elph_set%selfen_band_stop
                        call elph_mels_pw_set_nbandsk(self%pw,start_k,stop_k)
                        call elph_mels_pw_set_nbandskp(self%pw,1,self%pw%nbands)
                end select
                ! set number of bands in elph_mels
                self%nbands = self%pw%nbands
                self%nbands_local = self%pw%nbands_local
                self%nbands_k = self%pw%nbands_k
                self%band_start_k = self%pw%band_start_k
                self%band_stop_k = self%pw%band_stop_k
                self%nbands_kp = self%pw%nbands_kp
                self%band_start_kp = self%pw%band_start_kp
                self%band_stop_kp = self%pw%band_stop_kp

                ! write information about the wavefunctions to screen
                call elph_mels_pw_wf_write(self%pw)
                call elph_mels_pw_cache_init(self%pw)

                ! initialize windows for one-sided communication of the wavefunctions
                call elph_mels_pw_windows_init(self%pw)
                ! initialize the wavefunction caches to store the distributed wavefunctions
                call elph_mels_pw_wave_caches_init(self%pw)

                ! redistribute the wavefunctions to the processes that need them while balancing the workload
                ! and prefill the cache with the wavefunctions
                call elph_mels_pw_wf_redistribute(self%pw,elph_accumulators,kpoints_inter)
                call elph_mels_pw_cache_prefill(self%pw,elph_accumulators)

                call wave_map_get_ikpfbz(self%pw%wave_map,self%pw%waverot,self%latt_cur,kpoints_inter,self%pw%ikp_fbz_list)

            case(elph_mels_use_wann)
                ! nothing to do here
                num_wann = self%wann%wann_interp%nw
                self%nbands = num_wann
                self%nbands_local = num_wann
                self%nbands_k = num_wann
                self%band_start_k = 1
                self%band_stop_k = num_wann
                self%nbands_kp = num_wann
                self%band_start_kp = 1
                self%band_stop_kp = num_wann
                call wave_map_get_ikpfbz(self%wann%wave_map,self%wann%waverot,self%latt_cur,kpoints_inter,self%wann%ikp_fbz_list)
        end select
    end subroutine elph_mels_g_init

    !> Return a pointer to the wave_map object in elph_mels
    subroutine elph_mels_set_wavemap_pointer(self,wavemap)
        use wave_interpolate, only: wave_mapper
        type(elph_mels_type),target :: self
        type(wave_mapper),pointer,intent(out) :: wavemap

        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                wavemap => self%pw%wave_map
            case(elph_mels_use_wann)
                wavemap => self%wann%wave_map
        end select
    end subroutine elph_mels_set_wavemap_pointer

    subroutine elph_mels_set_ikp_fbz_list_pointer(self,ikp_fbz_list)
        use wave_interpolate, only: wave_mapper
        type(elph_mels_type),target :: self
        integer,pointer,intent(out) :: ikp_fbz_list(:)

        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                ikp_fbz_list => self%pw%ikp_fbz_list
            case(elph_mels_use_wann)
                ikp_fbz_list => self%wann%ikp_fbz_list
        end select
    end subroutine elph_mels_set_ikp_fbz_list_pointer

    !> compute g matrix element at a specific k and k' pair
    !> where k and k' are in the ibz
    subroutine elph_mels_compute_g_ibz(self,elph_set,kpoints,ikfbz,ikpfbz,isp,gcart_kkp,gmode_kkp)
        use elphon_mels_pw, ONLY: elph_mels_pw_compute_g_ibz
        use elphon_mels_wann, ONLY: elph_mels_wann_compute_g_ibz
        use elphon_base, only: get_g_kkp_modes
        use elphon_potential, only: phonon_ifc_get_phonons_zero_disp
        type(elph_mels_type) :: self
        type(elph_settings),intent(in) :: elph_set
        type(kpoints_struct),intent(in) :: kpoints
        integer,intent(in) :: ikfbz ! k-point index in the full BZ
        integer,intent(in) :: ikpfbz ! k' point index in the full BZ
        integer,intent(in) :: isp ! spin index
        complex(q),intent(out) :: gcart_kkp(:,:,:,:) ! g matrix element in the cartesian coordinates
        complex(q),intent(out) :: gmode_kkp(:,:,:) ! g matrix element in the mode coordinates
        ! local variables
        real(q) :: vqpt(3) ! q-point in the IBZ
        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                call elph_mels_pw_compute_g_ibz(self%pw,elph_set,kpoints,ikfbz,ikpfbz,isp,vqpt,gcart_kkp)
                ! Interpolate the dynamical matrix and compute phonon frequencies and eigenvectors
                call phonon_ifc_get_phonons_zero_disp(self%elph_ifc,vqpt,self%w_q,self%e_q)
                ! Transform the electron-phonon matrix elements from cartesian to mode basis
                call get_g_kkp_modes(self%natoms_pc,gcart_kkp,self%e_q,gmode_kkp)
            case(elph_mels_use_wann)
                call elph_mels_wann_compute_g_ibz(self%wann,elph_set,kpoints,ikfbz,ikpfbz,isp,vqpt,gcart_kkp)
                ! Interpolate the dynamical matrix and compute phonon frequencies and eigenvectors
                call phonon_ifc_get_phonons_zero_disp(self%elph_ifc,vqpt,self%w_q,self%e_q)
                ! Transform the electron-phonon matrix elements from cartesian to mode basis
                call get_g_kkp_modes(self%natoms_pc,gcart_kkp,self%e_q,gmode_kkp)
        end select
    end subroutine elph_mels_compute_g_ibz

    !> compute g matrix element at a specific k and k' pair
    !> where k and k' can be in the full BZ
    subroutine elph_mels_compute_g(self,elph_set,kpoints,ikfbz,ikpfbz,iqfbz,isp,gcart_kkp,gmode_kkp)
        use elphon_mels_pw, ONLY: elph_mels_pw_compute_g
        use elphon_mels_wann, ONLY: elph_mels_wann_compute_g
        use elphon_base, only: get_g_kkp_modes
        use elphon_potential, only: phonon_ifc_get_phonons_zero_disp
        type(elph_mels_type) :: self
        type(elph_settings),intent(in) :: elph_set
        type(kpoints_struct),intent(in) :: kpoints
        integer,intent(in) :: ikfbz ! k-point index in the full BZ
        integer,intent(in) :: ikpfbz ! k' point index in the full BZ
        integer,intent(in) :: iqfbz ! q-point index in the full BZ
        integer,intent(in) :: isp ! spin index
        complex(q),intent(out) :: gcart_kkp(:,:,:,:) ! g matrix element in the cartesian coordinates
        complex(q),intent(out) :: gmode_kkp(:,:,:) ! g matrix element in the mode coordinates
        ! local variables
        real(q) :: vqpt(3) ! q-point in the IBZ
        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                call elph_mels_pw_compute_g(self%pw,elph_set,kpoints,ikfbz,ikpfbz,iqfbz,isp,vqpt,gcart_kkp)
                ! Interpolate the dynamical matrix and compute phonon frequencies and eigenvectors
                call phonon_ifc_get_phonons_zero_disp(self%elph_ifc,vqpt,self%w_q,self%e_q)
                ! Transform the electron-phonon matrix elements from cartesian to mode basis
                call get_g_kkp_modes(self%natoms_pc,gcart_kkp,self%e_q,gmode_kkp)
            case(elph_mels_use_wann)
                call elph_mels_wann_compute_g(self%wann,elph_set,kpoints,ikfbz,ikpfbz,iqfbz,isp,vqpt,gcart_kkp)
                ! Interpolate the dynamical matrix and compute phonon frequencies and eigenvectors
                call phonon_ifc_get_phonons_zero_disp(self%elph_ifc,vqpt,self%w_q,self%e_q)
                ! Transform the electron-phonon matrix elements from cartesian to mode basis
                call get_g_kkp_modes(self%natoms_pc,gcart_kkp,self%e_q,gmode_kkp)
        end select
    end subroutine elph_mels_compute_g

    !> synchronize all the processes and free wfs and one-sided comm windows
    subroutine elph_mels_barrier(self,comm_global)
        USE elphon_mels_pw, ONLY: elph_mels_pw_windows_free, &
                                  elph_mels_pw_wave_caches_free
        type(elph_mels_type) :: self
        type(communic),intent(in) :: comm_global
        PROFILING_START('elph_inbalance')
        CALLMPI( M_barrier( comm_global) )
        PROFILING_STOP('elph_inbalance')
        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                call elph_mels_pw_windows_free(self%pw)
                call elph_mels_pw_wave_caches_free(self%pw)
            case(elph_mels_use_wann)
                ! nothing to do here
        end select
    end subroutine elph_mels_barrier

    !> free memory
    subroutine elph_mels_free(self,kpoints,w,nonl_s)
        use wave_struct_def, ONLY: wavespin
        use mkpoints_struct_def, ONLY: kpoints_struct
        use nonl_struct_def, ONLY: nonl_struct
        use elphon_potential, ONLY: phonon_ifc_free
        use elphon_mels_pw, ONLY: elph_mels_pw_cache_free, elph_mels_pw_pot_free, elph_mels_pw_wfs_free
        use wannier_interpolation, only: wannier_interpol_free
        type(elph_mels_type) :: self
        type(kpoints_struct) :: kpoints
        type(wavespin) :: w
        type(nonl_struct) :: nonl_s
        select case(self%elph_mels_use)
            case(elph_mels_use_pw)
                ! free pw objects
                call elph_mels_pw_cache_free(self%pw)
                call elph_mels_pw_pot_free(self%pw)
                call elph_mels_pw_wfs_free(self%pw,kpoints,w,nonl_s)
                call phonon_ifc_free(self%elph_ifc)
            case(elph_mels_use_wann)
                ! free wannier objects
                call wannier_interpol_free(self%wann%wann_interp)
        end select
    end subroutine elph_mels_free

end module elphon_mels
