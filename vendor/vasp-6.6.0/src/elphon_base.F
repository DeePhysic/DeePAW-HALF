#include "symbol.inc"

!> @brief This module provides low-level routines for handling electron-phonon calculations
module elphon_base
    use base
    use vhdf5
    use tutor, only: vtutor
    use string, only: str
    use incar_reader, only: process_incar, incar_file, count_elements
    use wave_interpolate, only: wave_interpolator_settings
    implicit none

    !> enumerate all possible scattering approximations
    enum,bind(c)
         enumerator :: elph_scattering_approx_crta,&
                       elph_scattering_approx_serta,&
                       elph_scattering_approx_erta_lambda,&
                       elph_scattering_approx_erta_tau,&
                       elph_scattering_approx_mrta_lambda,&
                       elph_scattering_approx_mrta_tau
    end enum

    !> options for compution the workload for each MPI rank
    enum,bind(c)
         enumerator :: elph_workload_mode_ibz,&
                       elph_workload_mode_fbz
    end enum

    !> options for wf redistribution
    enum,bind(c)
         enumerator :: elph_comm_twosided,&
                       elph_comm_onesided
    end enum

    !> type containing the approximations to be used in the computation of the electron self-energy due to electron-phonon coupling
    type elph_selfen_approx_t
        !> use static formulation of the self-energy
        logical :: static
        !> enum for the different scattering approximations to be used to solve the Boltzmann equation
        integer :: scattering_approx
    end type


    !> types of driver for the electron-phonon calculation
    enum, bind(c)
        enumerator :: elph_driver_el,&
                      elph_driver_ph,&
                      elph_driver_mels
    end enum

    !> modes to compute electron-phonon matrix elements
    enum,bind(c)
        enumerator :: elph_mels_use_pw,&
                      elph_mels_use_wann
    end enum

    !> @brief Container for all the input variables related to the electron-phonon calculation
    type elph_settings
        ! dimensions
        integer :: lmdim
        integer :: ngridpoints !< number of grid points in supercell = PROD(FFT_MESH)*DET_P2S
        integer :: det_p2s !< vol supercell/vol primitive cell
        integer :: fft_mesh(3) !< fft mesh = grid mesh in the primitive cell
        integer :: natoms_sc !< number of atoms in the supercell
        integer :: natoms_pc !< number of atoms in the primitive cell
        integer :: natoms_pc_local !< number of atoms in the primitive cell
        integer :: pot_ncdij !< number of components of the potential (1,2 or 4)
        integer :: wf_ncdij !< number of components of the wf
                            !< ncdij=1 unpolarized calculation
                            !< ncdij=2 spin-polarized calculation and
                            !< ncdij=4 noncollinear calculation
        integer :: npw !< total number of plane-waves
        integer :: npw_local !< storage needed for in-place FFT (grid%mplwv)
        integer :: npr_local !< local number of real-space grid points for fft = ELPH_SET%NP=GRID%RL%NP
        integer :: ispin !< number of spins
        integer :: nproj !< local number of projectors
        integer :: nproj_tot !< total number of projectors
        integer :: nrspinors !< number of spinors (1 for collinear, 2 for non collinear)
        integer :: ntyp
        integer,allocatable :: nityp(:)
        integer,allocatable :: lmmax(:)
        ! dimensions for phonons
        integer :: phonon_det_p2s !< vol supercell/vol primitive cell
        integer :: phonon_natoms_sc !< number of atoms in the phonon supercell
        ! options for the run
        integer :: nbands !< number of bands to be computed in the electron-phonon driver
        integer :: nbands_k !< number of bands for k
        integer :: band_start_k !< index of the first band for k
        integer :: band_stop_k !< index of the last band for k
        integer :: band_start_kp !< index of the first band for k'
        integer :: band_stop_kp !< index of the last band for k'
        type(wave_interpolator_settings) :: wave_interp_set !< options for computing the orbitals NSCF in the KPOINTS_ELPH kgrid
        integer,allocatable :: nbands_sum(:) !< number of bands to sum over
        integer,allocatable :: selfen_ikpt(:) !< index of the ibz points for which to compute the self-energy
        real(q),allocatable :: selfen_kpts(:,:) !< coordinates of the ibz kpoints for which to compute the self-energy
        integer :: selfen_band_start !< first band for which to compute the self-energy
        integer :: selfen_band_stop  !< last band for which to compute the self-energy
        integer :: selfen_band_start_kp !< first band for which to compute the self-energy
        integer :: selfen_band_stop_kp !< last band for which to compute the self-energy
        integer :: selfen_nw !< number of frequency points at which to compute the self-energy
        logical :: run !< whether to run the electron-phonon driver
        logical :: prepare !< whether to prepare the electron-phonon calculation
        logical :: metagga !< use metagga functional
        integer :: kpar !< number of kpoints per core
        integer :: ncore !< number of cores to use per-band
        integer :: driver !< logical controlling wether to run the electron or phonon self-energy
        type(elph_selfen_approx_t),allocatable :: selfen_approx(:) !< use static formulation of the self-energy
        logical :: selfen_fan !< compute fan self-energy
        logical :: selfen_dfan !< compute fan self-energy
        logical :: selfen_dw !< compute dw self-energy
        logical :: selfen_fans !< compute fan self-energy
        logical :: selfen_dws !< compute dw self-energy
        logical :: transport !< compute transport properties
        integer :: transport_driver !< choose which driver to use to compute transport properties
        integer :: fermi_nedos !< number of energy points for fermi energy computation
        integer :: transport_nedos !< number of energy points for transport computation
        integer :: transport_nedos_plot !< number of energy points for plotting the transport function
        real(q) :: transport_emin_plot !< min of energy grid for plotting of transport function
        real(q) :: transport_emax_plot !< max of energy grid for plotting of transport function
        real(q) :: transport_emin !< min of energy grid for computing the transport function
        real(q) :: transport_emax !< max of energy grid for computing the transport function
        logical :: superconductivity !< compute superconductivity
        integer :: nedos !< number of energy points at which to compute DOS
        real(q) :: kspacing !< k-spacing used to generate generalized regular mesh
        real(q) :: selfen_wrange !< energy window at which to compute the self-energy
        real(q),allocatable :: selfen_delta(:) !< imaginary energy shift to compute the self-energy
        integer :: selfen_ntemps !< number of temperatures at which to compute the self-energy
        real(q),allocatable :: selfen_temps(:) !< temperatures at which to compute the self-energy
        real(q) :: wf_cache_mb !< size fo the cache for the WFs in MB (default=0)
        logical :: wf_cache_prefill !< choose whether to fill in the cache before the start of the calculation
        logical :: use_pot_cache !< choose whether to cache the results of the interpolation of the electron-phonon potential
        logical :: wf_redistribute !< redistribute the orbitals before starting the computation
        integer :: wf_redistribute_opt !< option for redistributing the orbitals before starting the computation
        integer :: wf_comm_opt !< option to communicate wfs between nodes
        logical :: writemels !< write electron-phonon matrix elements to hdf5 file
        logical :: usehdf5mpiio !< whether or not to use hdf5 with MPIIO support
        logical :: useblas !< use blas implementation (default = .true.)
        logical :: userecip !< compute the electron-phonon matrix element using g-space
        logical :: useshmem !< use shared memory
        logical :: writepot !< write the electron-phonon potential to a file for visualization
        logical :: rotateprojectors !< rotate the projectors from the IBZ instead of recomputing them
        integer :: lr_sft !< use slow fourier transform to compute the long-range part of the potential in the supercell
                          !< -1 - no
                          !<  0 - only use if FFT requires too much storage (heuristic)
                          !<  1 - yes
        integer :: mels_use !< which mode to use to compute electron-phonon matrix elements
        logical :: lgaddv !< add the contribution from g_v
        logical :: lgaddt !< add the contribution from g_t
        logical :: lgaddd !< add the contribution from g_d
        logical :: lgaddr !< add the contribution from g_r
        logical :: lgaddp !< add the contribution from g_p
        logical :: lgaddq !< add the contribution from g_q
        logical :: lgadds !< add the contribution from g_q
        logical :: dry_run !< skip the main loop over k and k'
        integer :: velocity_mode    !< select mode to compute the velocity matrix elements.
                                    !< Currently two approaches are implemented:
                                    !< 1 - Laurent routines which use a AE approach
                                    !< 2 - Vasp routines set to use PS approach (currently might use more memory)
        logical :: write_hdf5vel !< write the electron group velocities to the hdf5 file
        logical :: write_textvel !< write the electron group velocities to a text file
        real(q) :: qval1 !< quadrupoles for atom 1 (this is only temporary)
        real(q) :: qval2 !< quadrupoles for atom 2 (this is only temporary)
        integer :: elph_lr !< treatment of the long-range part of the potential
        integer :: ifc_lr !< treatment of the long-range part of the force-constants
        real(q) :: encutlr !< g-space cutoff for the long-range contribution
        real(q) :: rcut !< radial cutof for coulomb potential
        integer :: ifc_asr  !< impose acoustic sum rule (ASR) to the force-constants.
                            !< if positive represents number of iterations to perform,
                            !< if negative the ASR is not imposed
        logical :: ignore_imag_phonons !< continue the calculation even if imag frequencies are present
        logical :: mels_cart !< write the electron-phonon matrix elements in basis of cartesian displacements instead of mode displacements
        integer :: selfen_ncarrier_den !< number of carrier densities at which to compute the self-energy
        real(q),allocatable :: selfen_carrier_den(:) !< carrier densities at which to compute the self-energy in cm^-3
        integer :: selfen_ncarrier_per_cell !< number of doping  electron concentration per cell at which to compute the self-energy
        real(q),allocatable :: selfen_carrier_per_cell(:) !< doping  electron concentration per cell at which to compute the self-energy
        integer :: selfen_nmu !< number chemical potential with respect to the Fermi level
        real(q),allocatable :: selfen_mu(:) ! chemical potential with respect to the Fermi level
        integer :: ismear !< smearing option to use to determine the fermi level in the electron_phonon driver
        integer :: selfen_naccumulators !< total number of accumulators to be used
        logical :: selfen_imag_skip !< skip the computation of q and k conbinations for which the delta is zero (requires tetrahedron method)
        !> set the energy window beyond which the delta function obtained from
        !> the imaginary part of the electron self-energy is considered to be zero.
        !> The energy window is computed such that a fraction (between 0 and 1) of the integral
        !> of the broadening function is excluded.
        !> This is only used when ELPH_SELFEN_IMAG_SKIP=.TRUE. and ELPH_SELFEN_DELTA>0
        real(q) :: selfen_broad_tol
        logical :: selfen_g_skip !< skip the computation of the explicit computation of the g matrix elements
        !> Set the the energy window beyond which the derivative of the fermi dirac
        !> distribution is considered to be zero. The energy window is computed such that a fraction (between 0 and 1)
        !> of the integral of the derivative of the fermi function is excluded from the integrals of the transport function.
        !> This is only used when ELPH_TRANSPORT_DRIVER = 1.
        real(q) :: transport_dfermi_tol
        real(q) :: transport_relaxation_time !< when > 0 use constant relaxation time to compute transport quantities
        real(q) :: selfen_energy_window(2) !< Energy cutoff (eV) for choosing which states to compute.
        logical :: selfen_gaps !< Find the direct and indirect gaps and select those states to be computed
        character(len=:),allocatable :: mode !< set up the default for the different modes = (transport, zpr, spectral, ...)
        !> For phonon self-energy
        integer :: nedos_ph ! number of energy grid points for phonon DOS
        real(q) :: emax_ph ! max phonon energy (eV)
        real(q) :: max_matsubara ! maximum matsubara frequency, in unit of the maximum phonon frequency
        integer :: matsubara_ntemps !< number of temperatures at which to compute matsubara freq
        real(q),allocatable :: matsubara_temps(:) !< temperatures at which to compute matsubara freq
        real(q) :: mu_star !! Morel Anderson mu*
    end type elph_settings

    !> @brief Container for arrays that can be re-used during the electron-phonon computation
    type elph_caches_t
        complex(q), allocatable :: v_cproj(:,:)
        complex(q), allocatable :: v_cr(:,:,:)
        complex(q), allocatable :: ccqij(:,:,:,:)
        complex(q), allocatable :: ccdij(:,:,:,:)
        complex(q), allocatable :: ccdrij(:,:,:,:,:)
    end type

    real(q), parameter :: FERMI_MAX_ARG=200   !< Maximum argument to pass to fermi-dirac function
                                                !< otherwise the limits are considered
    real(q), parameter :: BOSE_MAX_ARG=200 !< Maximum argument to pass to bose-einstein function
                                            !< otherwise the limits are considered

    contains

    subroutine elph_selfen_approx_init(self,static,scattering_approx)
        type(elph_selfen_approx_t),intent(out) :: self
        logical :: static !< use static approximation or not
        integer :: scattering_approx
        self%static = static
        self%scattering_approx = scattering_approx
    end subroutine

    function elph_selfen_approx_str(self,verbose) result(res)
        use string, only: str
        type(elph_selfen_approx_t),intent(in) :: self
        character(len=:), allocatable :: res
        logical,optional :: verbose
        !local variables
        logical :: my_verbose
        my_verbose=.false.
        if (present(verbose)) my_verbose = verbose
        if (my_verbose) then
            select case(self%scattering_approx)
                case(elph_scattering_approx_crta)
                    res = "constant relaxation-time approximation (CRTA)"
                case(elph_scattering_approx_serta)
                    res = "self-energy relaxation-time approximation (SERTA)"
                case(elph_scattering_approx_erta_lambda)
                    res = "energy relaxation-time approximation (ERTA_LAMBDA)"
                case(elph_scattering_approx_erta_tau)
                    res = "energy relaxation-time approximation (ERTA_TAU)"
                case(elph_scattering_approx_mrta_lambda)
                    res = "momentum relaxation-time approximation (MRTA_LAMBDA)"
                case(elph_scattering_approx_mrta_tau)
                    res = "momentum relaxation-time approximation (MRTA_TAU)"
            end select
        else
            select case(self%scattering_approx)
                case(elph_scattering_approx_crta)
                    res = "CRTA"
                case(elph_scattering_approx_serta)
                    res = "SERTA"
                case(elph_scattering_approx_erta_lambda)
                    res = "ERTA_LAMBDA"
                case(elph_scattering_approx_erta_tau)
                    res = "ERTA_TAU"
                case(elph_scattering_approx_mrta_lambda)
                    res = "MRTA_LAMBDA"
                case(elph_scattering_approx_mrta_tau)
                    res = "MRTA_TAU"
            end select
        endif
    end function

#ifdef VASP_HDF5
    subroutine elph_selfen_approx_hdf5write(self,groupid)
        type(elph_selfen_approx_t),intent(in) :: self
        integer(HID_T) :: groupid
        VH5_CHECK(vh5_write(groupid,'static',self%static))
        VH5_CHECK(vh5_write(groupid,'scattering_approximation',elph_selfen_approx_str(self)))
    end subroutine
#endif

    !> @brief Allocate some caches needed for electron-phonon calculations
    !> This is to avoid allocating/deallocating inside loops or
    !> spamming the main thread with variables
    subroutine elph_caches_init(self,nproj_tot,nbands,npr_local,nbands_k,nrspinors)
        use ini, only: register_allocate
        type(elph_caches_t) :: self
        integer,intent(in) :: nproj_tot
        integer,intent(in) :: nbands
        integer,intent(in) :: npr_local !< number of local real space points grid%rl%np
        integer,intent(in) :: nbands_k
        integer,intent(in) :: nrspinors
        allocate(self%v_cproj(nproj_tot,nbands))
        allocate(self%v_cr(npr_local,nrspinors,nbands_k))
        call register_allocate(VMEM_ARR_SIZE(self%v_cproj)+&
                               VMEM_ARR_SIZE(self%v_cr),"elph_cache")
    end subroutine elph_caches_init

    !> @brief Dellocate some caches needed for electron-phonon calculations
    subroutine elph_caches_free(self)
        use ini, only: deregister_allocate
        type(elph_caches_t) :: self
        call deregister_allocate(   VMEM_ARR_SIZE(self%v_cproj)+&
                                    VMEM_ARR_SIZE(self%v_cr),"elph_cache")
        deallocate(self%v_cproj)
        deallocate(self%v_cr)
        if (allocated(self%ccdij)) deallocate(self%ccdij)
        if (allocated(self%ccqij)) deallocate(self%ccqij)
        if (allocated(self%ccdrij)) deallocate(self%ccdrij)
    end subroutine elph_caches_free

    !> @brief Set band_start_k and band_stop_k and keep consistency with nbands_k
    subroutine elph_settings_set_nbandsk(self,start_k,stop_k)
        type(elph_settings) :: self
        integer,intent(in) :: start_k
        integer,intent(in) :: stop_k
        self%band_start_k=start_k
        self%band_stop_k=stop_k
        self%nbands_k=self%band_stop_k-self%band_start_k+1
    end subroutine elph_settings_set_nbandsk

    !> @brief Set defaults for general settings
    subroutine electron_phonon_reader_general_default(elphin,xc,wave_interp_set)
        use setexm, only: xc_info
        type(elph_settings), intent(inout) :: elphin
        type(xc_info),intent(in) :: xc
        type(wave_interpolator_settings), intent(inout) :: wave_interp_set

        elphin%run=.false.
        elphin%prepare=.false.
        ! we want to check for metagga and not the deorbitalized versions
        elphin%metagga=(xc%ldometagga.and.(xc%ltau.or.xc%lmu))
        elphin%driver=elph_driver_el
        elphin%dry_run=.false.
        elphin%nbands=-1
        elphin%nedos=5001
        elphin%kspacing=0.5
        elphin%encutlr=50._q
        elphin%rcut=0._q
        elphin%writepot=.false.
        elphin%qval1=0._q
        elphin%qval2=0._q
        elphin%wave_interp_set=wave_interp_set
    end subroutine

    !> @brief Read general settings
    subroutine electron_phonon_reader_general_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io

        integer :: ierr

        call process_incar(incar, 'elph_prepare', elphin%prepare, ierr)
        call process_incar(incar, 'elph_run', elphin%run, ierr)
        call process_incar(incar, 'elph_dry_run', elphin%dry_run, ierr)
        call process_incar(incar, 'elph_nbands', elphin%nbands, ierr)
        call process_incar(incar, 'elph_nedos', elphin%nedos, ierr)
        call process_incar(incar, 'elph_kspacing', elphin%kspacing, ierr)
        call process_incar(incar, 'encutlr', elphin%encutlr, ierr)
        call process_incar(incar, 'elph_writepot', elphin%writepot, ierr)
        call process_incar(incar, 'elph_rcut', elphin%rcut, ierr)
        call process_incar(incar, 'qval', elphin%qval1, ierr)
        elphin%qval2=-elphin%qval1
        call process_incar(incar, 'qval2', elphin%qval2, ierr)
    end subroutine

    !> @brief Set defaults for phonon-related settings
    subroutine electron_phonon_reader_phonon_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        elphin%ignore_imag_phonons=.false.
        elphin%mels_cart=.false.
        elphin%ifc_asr=1
        elphin%ifc_lr=1
    end subroutine

    !> @brief Read phonon-related settings
    subroutine electron_phonon_reader_phonon_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io

        logical :: lopen
        integer :: ierr

        call process_incar(incar, 'elph_ignore_imag_phonons', elphin%ignore_imag_phonons, ierr)
        call process_incar(incar, 'elph_mels_cart', elphin%mels_cart, ierr)
        call process_incar(incar, 'ifc_asr', elphin%ifc_asr, ierr)
        call process_incar(incar, 'ifc_lr', elphin%ifc_lr, ierr)
    end subroutine

    !> @brief Set defaults for settings related to matrix elements
    subroutine electron_phonon_reader_matrix_elements_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        elphin%mels_use = elph_mels_use_pw
        elphin%lgaddv = .true.
        elphin%lgaddt = elphin%metagga
        elphin%lgaddd = .true.
        elphin%lgaddr = .true.
        elphin%lgaddp = .true.
        elphin%lgaddq = .false.
        elphin%lgadds = .false.

        elphin%elph_lr=1
        elphin%writemels=.false.
        elphin%usehdf5mpiio=.false.
    end subroutine

    !> @brief Read settings related to matrix elements
    subroutine electron_phonon_reader_matrix_elements_read(incar, elphin, io)
        use string, only: lowercase
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io

        logical :: lopen
        integer :: ierr
        character(len=:),allocatable :: elph_decompose
        character(len=:),allocatable :: elph_mels
        character(len=:),allocatable :: elph_mels_kind
        character(len=:),allocatable :: writemels

        elph_mels="pw"
        call process_incar(incar, 'elph_mels', elph_mels, ierr)
        if (lowercase(elph_mels)=='wann') then
            elphin%mels_use = elph_mels_use_wann
        else if (lowercase(elph_mels)=='pw') then
            elphin%mels_use = elph_mels_use_pw
        endif

        elph_mels_kind=""
        call process_incar(incar, 'elph_mels_kind', elph_mels_kind, ierr)
        select case(lowercase(elph_mels_kind))
            case("ae")
                elphin%lgaddv = .true.
                elphin%lgaddt = elphin%metagga
                elphin%lgaddd = .true.
                elphin%lgaddr = .true.
                elphin%lgaddp = .true.
                elphin%lgaddq = .false.
                elphin%lgadds = .false.
            case("ps")
                elphin%lgaddv = .true.
                elphin%lgaddt = elphin%metagga
                elphin%lgaddd = .true.
                elphin%lgaddr = .false.
                elphin%lgaddp = .false.
                elphin%lgaddq = .true.
                elphin%lgadds = .false.
            case("")
                 ! do nothing, keep the defaults
            case default
                call vtutor%alert("unrecognized option for ELPH_MELS_KIND. Using default or ELPH_DECOMPOSE.")
        end select

        elph_decompose=""
        call process_incar(incar, 'elph_decompose', elph_decompose, ierr)
        elph_decompose=lowercase(elph_decompose)
        if (ierr==0) then
            elphin%lgaddv = index(elph_decompose, 'v') /= 0
            elphin%lgaddt = index(elph_decompose, 't') /= 0
            elphin%lgaddd = index(elph_decompose, 'd') /= 0
            elphin%lgaddr = index(elph_decompose, 'r') /= 0
            elphin%lgaddp = index(elph_decompose, 'p') /= 0
            elphin%lgaddq = index(elph_decompose, 'q') /= 0
            elphin%lgadds = index(elph_decompose, 's') /= 0
        endif

        call process_incar(incar, 'elph_lr', elphin%elph_lr, ierr)
        ! try to read if writemels=hdf5mpiio
        call process_incar(incar, 'elph_writemels', writemels, ierr)
        if (ierr==0) then
            writemels=lowercase(writemels)
            ! do we have hdf5mpiio?
            if (writemels=="hdf5mpiio") then
                elphin%usehdf5mpiio=.true.
                elphin%writemels=.true.
            endif
            ! do we have writemels?
            if (index(writemels, 't') /= 0) then
                elphin%writemels=.true.
            endif
        endif

    end subroutine

    !> @brief Set defaults for performance-related settings
    subroutine electron_phonon_reader_performance_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        elphin%wf_cache_prefill=.true.
        elphin%wf_cache_mb=1000
        elphin%use_pot_cache=.true.
        elphin%wf_redistribute=.false.
        elphin%wf_redistribute_opt=elph_workload_mode_ibz
        elphin%wf_comm_opt=elph_comm_twosided
        elphin%kpar=-1
        elphin%ncore=1
        elphin%useblas=.true.
        elphin%userecip=elphin%metagga
        elphin%useshmem=.false.
        elphin%rotateprojectors=.true.
        elphin%lr_sft=0
    end subroutine

    !> @brief Read performance-related settings
    subroutine electron_phonon_reader_performance_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io

        logical :: lopen
        integer :: ierr

        call process_incar(incar, 'elph_wf_cache_prefill', elphin%wf_cache_prefill, ierr)
        call process_incar(incar, 'elph_wf_cache_mb', elphin%wf_cache_mb, ierr)
        call process_incar(incar, 'elph_use_pot_cache', elphin%use_pot_cache, ierr)
        call process_incar(incar, 'elph_wf_redistribute', elphin%wf_redistribute, ierr)
        call process_incar(incar, 'elph_wf_redistribute_opt', elphin%wf_redistribute_opt, ierr)
        call process_incar(incar, 'elph_wf_comm_opt', elphin%wf_comm_opt, ierr)
        call process_incar(incar, 'elph_kpar', elphin%kpar, ierr)
        call process_incar(incar, 'ncore', elphin%ncore, ierr)
        call process_incar(incar, 'elph_useblas', elphin%useblas, ierr)
        call process_incar(incar, 'elph_userecip', elphin%userecip, ierr)
        call process_incar(incar, 'elph_useshmem', elphin%useshmem, ierr)
        call process_incar(incar, 'elph_rotateprojectors', elphin%rotateprojectors, ierr)
        call process_incar(incar, 'elph_lr_sft', elphin%lr_sft, ierr)
    end subroutine

    !> @brief Set defaults for settings related to chemical potential
    subroutine electron_phonon_reader_chemical_potential_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        elphin%ismear = 0
        elphin%fermi_nedos = 501 ! only used when ELPH_ISMEAR=-24

        elphin%selfen_ntemps = 6
        if (allocated(elphin%selfen_temps)) deallocate(elphin%selfen_temps)
        allocate(elphin%selfen_temps(elphin%selfen_ntemps))
        elphin%selfen_temps = [0, 100, 200, 300, 400, 500]

        ! carrier densities at which to compute the electron-phonon self-energy
        ! We have 3 tags to gives this information. Nb carrier per cm^3, nb carrier per cell, chemical potential.
        ! The precedence is defined by chemical potential > Nb carrier per cm^3 > nb carrier per cell
        elphin%selfen_ncarrier_per_cell = 1
        if (allocated(elphin%selfen_carrier_per_cell)) deallocate(elphin%selfen_carrier_per_cell)
        allocate(elphin%selfen_carrier_per_cell(elphin%selfen_ncarrier_per_cell))
        elphin%selfen_carrier_per_cell(1) = 0.0_q

        elphin%selfen_ncarrier_den=-1
        if (allocated(elphin%selfen_carrier_den)) deallocate(elphin%selfen_carrier_den)

        elphin%selfen_nmu=-1
        if (allocated(elphin%selfen_mu)) deallocate(elphin%selfen_mu)
    end subroutine

    !> @brief Read settings related to chemical potential
    subroutine electron_phonon_reader_chemical_potential_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io
        ! Local variables
        logical :: lopen
        integer :: ierr
        integer :: i, irange, ishift, j, n
        integer :: num, num_d, num_c, num_m, num_dr, num_cr, num_mr
        real(q) :: carrier_sign
        real(q) :: valr
        real(q) :: delta
        real(q),allocatable :: selfen_carrier_den_range(:)
        real(q),allocatable :: selfen_carrier_per_cell_range(:)
        real(q) :: selfen_mu_range(3)
        real(q) :: selfen_temps_range(3)

        call process_incar(incar, 'ismear', elphin%ismear, ierr)
        call process_incar(incar, 'elph_ismear', elphin%ismear, ierr)
        call process_incar(incar, 'elph_fermi_nedos',elphin%fermi_nedos, ierr)

        ! temperatures at which to compute the self-energy
        num = count_elements(incar, "elph_selfen_temps")
        if (num .ne. 0) then
            elphin%selfen_ntemps = num
            deallocate(elphin%selfen_temps)
            allocate(elphin%selfen_temps(num))
            call process_incar(incar, 'elph_selfen_temps', elphin%selfen_temps, ierr)
        endif
        num_dr = count_elements(incar, 'elph_selfen_temps_range')
        if (num_dr==3) then
            call process_incar(incar, 'elph_selfen_temps_range', selfen_temps_range, ierr)
            num_d = int(selfen_temps_range(3))
            elphin%selfen_ntemps = num_d
            if (allocated(elphin%selfen_temps)) deallocate(elphin%selfen_temps)
            allocate(elphin%selfen_temps(elphin%selfen_ntemps))
            delta =(selfen_temps_range(2)-selfen_temps_range(1))/(elphin%selfen_ntemps-1)
            do i=1,elphin%selfen_ntemps
                elphin%selfen_temps(i) = (i-1)*delta+selfen_temps_range(1)
            enddo
        endif

        ! carrier densities at which to compute the electron-phonon self-energy
        num_d = count_elements(incar, 'elph_selfen_carrier_den')
        if (num_d .ne. 0) then
            elphin%selfen_ncarrier_den = num_d
            if (allocated(elphin%selfen_carrier_den)) deallocate(elphin%selfen_carrier_den)
            allocate(elphin%selfen_carrier_den(num_d))
            call process_incar(incar, 'elph_selfen_carrier_den', elphin%selfen_carrier_den, ierr)
        endif
        num_dr = count_elements(incar, 'elph_selfen_carrier_den_range')
        if (mod(num_dr,3)==0.and.num_dr/=0) then
            allocate(selfen_carrier_den_range(num_dr))
            call process_incar(incar, 'elph_selfen_carrier_den_range', selfen_carrier_den_range, ierr)
            num_d = 0
            ! count number of points in all the ranges
            do irange=1,num_dr/3
                ishift = (irange-1)*3
                n = int(selfen_carrier_den_range(ishift+3))
                num_d = num_d + n
            enddo
            elphin%selfen_ncarrier_den = num_d
            if (allocated(elphin%selfen_carrier_den)) deallocate(elphin%selfen_carrier_den)
            allocate(elphin%selfen_carrier_den(elphin%selfen_ncarrier_den))
            ! generate points for all the ranges
            j=1
            do irange=1,num_dr/3
                ishift = (irange-1)*3
                n = int(selfen_carrier_den_range(ishift+3))
                !check if the sign is the same
                if (sign(1.0_q,selfen_carrier_den_range(ishift+1)) /= sign(1.0_q,selfen_carrier_den_range(ishift+2))) then
                    call vtutor%error('The lower and upper bounds in elph_selfen_carrier_den_range '//str(selfen_carrier_den_range(ishift+1))//&
                                      ' and ' //str(selfen_carrier_den_range(ishift+2))//' have different signs. This is not allowed.')
                endif
                carrier_sign = sign(1.0_q,selfen_carrier_den_range(ishift+1))
                selfen_carrier_den_range(ishift+1) = log10(abs(selfen_carrier_den_range(ishift+1)))
                selfen_carrier_den_range(ishift+2) = log10(abs(selfen_carrier_den_range(ishift+2)))
                elphin%selfen_carrier_den(j:j+n-1) = carrier_sign*10**(linspace(selfen_carrier_den_range(ishift+1),selfen_carrier_den_range(ishift+2),n))
                j=j+n
            enddo
        endif

        num_c = count_elements(incar, 'elph_selfen_carrier_per_cell')
        if (num_c .ne. 0) then
            elphin%selfen_ncarrier_per_cell = num_c
            if (allocated(elphin%selfen_carrier_per_cell)) deallocate(elphin%selfen_carrier_per_cell)
            allocate(elphin%selfen_carrier_per_cell(num_c))
            call process_incar(incar, 'elph_selfen_carrier_per_cell', elphin%selfen_carrier_per_cell, ierr)
        endif
        num_cr = count_elements(incar, 'elph_selfen_carrier_per_cell_range')
        if (mod(num_cr,3)==0.and.num_cr/=0) then
            allocate(selfen_carrier_per_cell_range(num_cr))
            call process_incar(incar, 'elph_selfen_carrier_per_cell_range', selfen_carrier_per_cell_range, ierr)
            num_c = 0
            ! count number of points in all the ranges
            do irange=1,num_cr/3
                ishift = (irange-1)*3
                n = int(selfen_carrier_per_cell_range(ishift+3))
                num_c = num_c + n
            enddo
            elphin%selfen_ncarrier_per_cell = num_c
            if (allocated(elphin%selfen_carrier_per_cell)) deallocate(elphin%selfen_carrier_per_cell)
            allocate(elphin%selfen_carrier_per_cell(elphin%selfen_ncarrier_per_cell))
            ! generate points for all the ranges
            j=1
            do irange=1,num_cr/3
                ishift = (irange-1)*3
                n = int(selfen_carrier_per_cell_range(ishift+3))
                elphin%selfen_carrier_per_cell(j:j+n-1) = linspace(selfen_carrier_per_cell_range(ishift+1),selfen_carrier_per_cell_range(ishift+2),n)
                j=j+n
            enddo
        endif

        num_m = count_elements(incar, 'elph_selfen_mu')
        if (num_m .ne. 0) then
            elphin%selfen_nmu = num_m
            if (allocated(elphin%selfen_mu)) deallocate(elphin%selfen_mu)
            allocate(elphin%selfen_mu(num_m))
            call process_incar(incar, 'elph_selfen_mu', elphin%selfen_mu, ierr)
        endif
        num_mr = count_elements(incar, 'elph_selfen_mu_range')
        if (num_mr==3) then
            call process_incar(incar, 'elph_selfen_mu_range', selfen_mu_range, ierr)
            num_m = int(selfen_mu_range(3))
            elphin%selfen_nmu = num_m
            if (allocated(elphin%selfen_mu)) deallocate(elphin%selfen_mu)
            allocate(elphin%selfen_mu(elphin%selfen_nmu))
            delta = (selfen_mu_range(2)-selfen_mu_range(1))/(elphin%selfen_nmu-1)
            do i=1,elphin%selfen_nmu
                elphin%selfen_mu(i) = (i-1)*delta+selfen_mu_range(1)
            enddo
        endif

        if (num_m > 0) then
            if (allocated(elphin%selfen_carrier_den))   deallocate(elphin%selfen_carrier_den)
            elphin%selfen_ncarrier_den=-1
            if (allocated(elphin%selfen_carrier_per_cell)) deallocate(elphin%selfen_carrier_per_cell)
            elphin%selfen_ncarrier_per_cell=-1
        elseif (num_d > 0) then
            if (allocated(elphin%selfen_carrier_per_cell)) deallocate(elphin%selfen_carrier_per_cell)
            elphin%selfen_ncarrier_per_cell=-1
            if (allocated(elphin%selfen_mu)) deallocate(elphin%selfen_mu)
            elphin%selfen_nmu=-1
        elseif (num_c > 0) then
            if (allocated(elphin%selfen_carrier_den)) deallocate(elphin%selfen_carrier_den)
            elphin%selfen_ncarrier_den=-1
            if (allocated(elphin%selfen_mu)) deallocate(elphin%selfen_mu)
            elphin%selfen_nmu=-1
        endif
    end subroutine

    !> @brief Set defaults for settings related to self energy
    subroutine electron_phonon_reader_self_energy_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        allocate(elphin%selfen_approx(1))
        call elph_selfen_approx_init(elphin%selfen_approx(1),static=.false.,scattering_approx=elph_scattering_approx_serta)
        elphin%selfen_fan=.false.
        elphin%selfen_dfan=.false.
        elphin%selfen_dw=.false.
        elphin%selfen_fans=.false.
        elphin%selfen_dws=.false.

        allocate(elphin%nbands_sum(1))
        elphin%nbands_sum=-1
        elphin%selfen_band_start=-1
        elphin%selfen_band_stop=-1
        elphin%selfen_band_start_kp=-1
        elphin%selfen_band_stop_kp=-1

        if (allocated(elphin%selfen_ikpt)) deallocate(elphin%selfen_ikpt)
        if (allocated(elphin%selfen_kpts)) deallocate(elphin%selfen_kpts)

        elphin%selfen_nw=1
        elphin%selfen_wrange=0
        elphin%selfen_gaps=.false.

        allocate(elphin%selfen_delta(1))
        elphin%selfen_delta=0.01_q
        elphin%selfen_energy_window=0.0_q
        elphin%selfen_g_skip=.false.
        elphin%selfen_imag_skip=.false.
        elphin%selfen_broad_tol=0.0
    end subroutine

    !> @brief Read settings related to self energy
    subroutine electron_phonon_reader_self_energy_read(incar, elphin, io)
        use string, only: lowercase
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io
        ! logical variables
        logical :: lopen
        logical :: static
        integer :: ierr
        integer :: val
        integer :: num
        real(q) :: valr
        logical :: crta_present, serta_present
        logical :: erta_lambda_present, erta_tau_present
        logical :: mrta_lambda_present, mrta_tau_present
        character(len=:),allocatable :: inplin

        static=.false.
        call process_incar(incar, 'elph_selfen_static', static, ierr)

        inplin = 'serta'

        call process_incar(incar, 'elph_scattering_approx', inplin, ierr)
        call scattering_approx_from_str(elphin%selfen_approx, static, inplin)

        call process_incar(incar, 'elph_selfen_fan', elphin%selfen_fan, ierr)
        call process_incar(incar, 'elph_selfen_dfan', elphin%selfen_dfan, ierr)
        call process_incar(incar, 'elph_selfen_dw', elphin%selfen_dw, ierr)
        call process_incar(incar, 'elph_selfen_fans', elphin%selfen_fans, ierr)
        call process_incar(incar, 'elph_selfen_dws', elphin%selfen_dws, ierr)

        num = count_elements(incar, "elph_nbands_sum")
        if (num .ne. 0) then
            if (allocated(elphin%nbands_sum)) deallocate(elphin%nbands_sum)
            allocate(elphin%nbands_sum(num))
            call process_incar(incar, 'elph_nbands_sum', elphin%nbands_sum, ierr)
        endif

        call process_incar(incar, 'elph_selfen_band_start', elphin%selfen_band_start, ierr)
        call process_incar(incar, 'elph_selfen_band_stop', elphin%selfen_band_stop, ierr)
        call process_incar(incar, 'elph_selfen_band_start_kp', elphin%selfen_band_start_kp, ierr)
        call process_incar(incar, 'elph_selfen_band_stop_kp', elphin%selfen_band_stop_kp, ierr)


        num = count_elements(incar, 'elph_selfen_ikpt')
        if (num/=0) then
            if (allocated(elphin%selfen_ikpt)) deallocate(elphin%selfen_ikpt)
            allocate(elphin%selfen_ikpt(num))
            call process_incar(incar, 'elph_selfen_ikpt', elphin%selfen_ikpt, ierr)
        endif

        num = count_elements(incar, 'elph_selfen_kpts')
        if (num/=0) then
            if (allocated(elphin%selfen_kpts)) deallocate(elphin%selfen_kpts)
            allocate(elphin%selfen_kpts(3,num/3))
            call process_incar(incar, 'elph_selfen_kpts', elphin%selfen_kpts, ierr)
        endif

        call process_incar(incar, 'elph_selfen_nw', elphin%selfen_nw, ierr)
        call process_incar(incar, 'elph_selfen_enwin', elphin%selfen_wrange, ierr)
        if (ierr==0) call vtutor%alert('the ELPH_SELFEN_ENWIN incar tag is deprecated, use ELPH_SELFEN_WRANGE instead.')
        call process_incar(incar, 'elph_selfen_wrange', elphin%selfen_wrange, ierr)
        call process_incar(incar, 'elph_selfen_gaps', elphin%selfen_gaps, ierr)

        num = count_elements(incar, 'elph_selfen_delta')
        if (num .ne. 0) then
            if (allocated(elphin%selfen_delta)) deallocate(elphin%selfen_delta)
            allocate(elphin%selfen_delta(num))
            call process_incar(incar, 'elph_selfen_delta', elphin%selfen_delta, ierr)
        endif

        call process_incar(incar, 'elph_selfen_energy_window', elphin%selfen_energy_window, ierr)
        call process_incar(incar, 'elph_selfen_g_skip', elphin%selfen_g_skip, ierr)
        call process_incar(incar, 'elph_selfen_imag_skip',elphin%selfen_imag_skip, ierr)
        call process_incar(incar, 'elph_selfen_broad_tol',elphin%selfen_broad_tol, ierr)
    end subroutine

    !> @brief Set defaults for settings related to self energy
    subroutine electron_phonon_reader_phonon_self_energy_default(elphin)
        type(elph_settings), intent(inout) :: elphin

    end subroutine electron_phonon_reader_phonon_self_energy_default

    !> @brief Set defaults for settings related to self energy
    !> These can still be overwritten by explicit INCAR tags
    subroutine electron_phonon_reader_phonon_self_energy_meta(elphin)
        type(elph_settings), intent(inout) :: elphin
        ! use tetrahedron method
        if (allocated(elphin%selfen_delta)) deallocate(elphin%selfen_delta)
        allocate(elphin%selfen_delta(1))
        elphin%selfen_delta=0.0_q
        ! computation of chemical potential
        elphin%ismear = -24
        elphin%fermi_nedos = 501 ! only used when ELPH_ISMEAR=-24

        ! with those inputs, the bound will be put to 1 and nbands by
        ! elph_accumulator_selfen_init
        if (allocated(elphin%nbands_sum)) deallocate(elphin%nbands_sum)
        allocate(elphin%nbands_sum(1))
        elphin%nbands_sum=-1
        elphin%selfen_band_start=-1
        elphin%selfen_band_stop=-1
        elphin%selfen_band_start_kp=-1
        elphin%selfen_band_stop_kp=-1

        elphin%wf_redistribute=.true.
    end subroutine electron_phonon_reader_phonon_self_energy_meta

   !>ead settings related to chemical potential
    subroutine electron_phonon_reader_phonon_self_energy_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io
        ! Local variables
        logical :: lopen
        integer :: ierr
        integer :: i
        integer :: num, num_d, num_dr

    end subroutine electron_phonon_reader_phonon_self_energy_read

    ! read a list of scattering approximation from a string
    subroutine scattering_approx_from_str(selfen_approx, static, inplin_in)
        use string, only: lowercase
        type(elph_selfen_approx_t),allocatable :: selfen_approx(:)
        logical,intent(in) :: static
        character(*),intent(in) :: inplin_in
        ! local variables
        integer :: num
        logical :: crta_present, serta_present
        logical :: erta_lambda_present, erta_tau_present
        logical :: mrta_lambda_present, mrta_tau_present
        character(len=256) :: inplin

        inplin = lowercase(inplin_in)
        ! for now just add one or two instances depending on wether mrta is present or not
        num = 0
        crta_present         = index(inplin, 'crta') > 0
        serta_present        = index(inplin, 'serta') > 0
        erta_lambda_present  = index(inplin, 'erta_lambda') > 0
        erta_tau_present     = index(inplin, 'erta_tau') > 0
        mrta_lambda_present  = index(inplin, 'mrta_lambda') > 0
        mrta_tau_present     = index(inplin, 'mrta_tau') > 0
        if (crta_present)        num=num+1
        if (serta_present)       num=num+1
        if (erta_lambda_present) num=num+1
        if (erta_tau_present)    num=num+1
        if (mrta_lambda_present) num=num+1
        if (mrta_tau_present)    num=num+1
        deallocate(selfen_approx)
        allocate(selfen_approx(num))
        num = 0
        if (crta_present) then
            num = num+1
            call elph_selfen_approx_init(selfen_approx(num),static=static,scattering_approx=elph_scattering_approx_crta)
        endif
        if (serta_present) then
            num = num+1
            call elph_selfen_approx_init(selfen_approx(num),static=static,scattering_approx=elph_scattering_approx_serta)
        endif
        if (erta_lambda_present) then
            num = num+1
            call elph_selfen_approx_init(selfen_approx(num),static=static,scattering_approx=elph_scattering_approx_erta_lambda)
        endif
        if (erta_tau_present) then
            num = num+1
            call elph_selfen_approx_init(selfen_approx(num),static=static,scattering_approx=elph_scattering_approx_erta_tau)
        endif
        if (mrta_lambda_present) then
            num = num+1
            call elph_selfen_approx_init(selfen_approx(num),static=static,scattering_approx=elph_scattering_approx_mrta_lambda)
        endif
        if (mrta_tau_present) then
            num = num+1
            call elph_selfen_approx_init(selfen_approx(num),static=static,scattering_approx=elph_scattering_approx_mrta_tau)
        endif
    end subroutine

    !> @brief Set defaults for transport-related settings
    subroutine electron_phonon_reader_transport_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        elphin%transport=.false.
        elphin%transport_driver=2
        elphin%transport_nedos=501
        elphin%transport_nedos_plot=-1
        elphin%transport_emin_plot=huge(elphin%transport_emin_plot)
        elphin%transport_emax_plot=huge(elphin%transport_emax_plot)
        elphin%transport_emin=huge(elphin%transport_emin)
        elphin%transport_emax=huge(elphin%transport_emax)
        elphin%transport_dfermi_tol=1e-6_q
        elphin%transport_relaxation_time=1e-14_q
        elphin%velocity_mode=2
        elphin%write_hdf5vel=.true.
        elphin%write_textvel=.false.
    end subroutine

    !> Read superconductivity-related settings
    subroutine electron_phonon_reader_superconductivity_default(elphin)
        type(elph_settings), intent(inout) :: elphin

        elphin%superconductivity=.false.

        ! phonon dos
        elphin%nedos_ph = 5001

        ! phonon maximum energy: if -1, later put to 101% of maximum
        ! frequency of the spectrum
        elphin%emax_ph = -1.0_q

        elphin%matsubara_ntemps = 6
        if (allocated(elphin%matsubara_temps)) deallocate(elphin%matsubara_temps)
        allocate(elphin%matsubara_temps(elphin%matsubara_ntemps))
        elphin%matsubara_temps = [0, 1, 2, 3, 4, 5]

        ! max matsubara frequency
        elphin%max_matsubara=10.0_q
        ! Morel Anderson pseudo potential
        elphin%mu_star=0.1_q
    end subroutine

    !> @brief Read transport-related settings
    subroutine electron_phonon_reader_transport_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io

        logical :: lopen
        integer :: ierr

        call process_incar(incar, 'elph_transport', elphin%transport, ierr)
        call process_incar(incar, 'elph_transport_driver',elphin%transport_driver, ierr)

        select case (elphin%transport_driver)
            case(1)
                elphin%transport_nedos=5001
            case(2)
                elphin%transport_nedos=501
            case default
                call vtutor%error('elph_transport_driver must be 1 or 2')
        end select

        call process_incar(incar, 'elph_transport_nedos', elphin%transport_nedos, ierr)
        call process_incar(incar, 'transport_nedos', elphin%transport_nedos, ierr)

        call process_incar(incar, 'elph_transport_nedos_plot', elphin%transport_nedos_plot, ierr)
        call process_incar(incar, 'transport_nedos_plot', elphin%transport_nedos_plot, ierr)

        call process_incar(incar, 'elph_transport_emin_plot', elphin%transport_emin_plot, ierr)
        call process_incar(incar, 'elph_transport_emax_plot', elphin%transport_emax_plot, ierr)

        call process_incar(incar, 'elph_transport_emin', elphin%transport_emin, ierr)
        call process_incar(incar, 'elph_transport_emax', elphin%transport_emax, ierr)

        call process_incar(incar, 'elph_transport_dfermi_tol', elphin%transport_dfermi_tol, ierr)
        call process_incar(incar, 'transport_dfermi_tol', elphin%transport_dfermi_tol, ierr)

        call process_incar(incar, 'transport_relaxation_time', elphin%transport_relaxation_time, ierr)
        call process_incar(incar, 'elph_velocity_mode', elphin%velocity_mode, ierr)
        call process_incar(incar, 'elph_write_hdf5vel', elphin%write_hdf5vel, ierr)
        call process_incar(incar, 'elph_write_textvel', elphin%write_textvel, ierr)
    end subroutine

    !> read inputs related to superconductivity
    subroutine electron_phonon_reader_superconductivity_read(incar, elphin, io)
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(inout) :: elphin
        type(in_struct),intent(in) :: io
        ! local variables
        integer :: i, ierr
        integer :: num, num_d, num_dr
        real(q) :: delta
        real(q) :: matsubara_temps_range(3)

        ! temperatures at which to compute the matsubara freq
        num = count_elements(incar, "elph_matsubara_temps")
        elphin%matsubara_ntemps = num
        if (num .ne. 0) then
            deallocate(elphin%matsubara_temps)
            allocate(elphin%matsubara_temps(num))
            call process_incar(incar, 'elph_matsubara_temps', elphin%matsubara_temps, ierr)
        endif
        num_dr = count_elements(incar, 'elph_matsubara_temps_range')
        if (num_dr==3) then
            call process_incar(incar, 'elph_matsubara_temps_range', matsubara_temps_range, ierr)
            num_d = int(matsubara_temps_range(3))
            elphin%matsubara_ntemps = num_d
            if (allocated(elphin%matsubara_temps)) deallocate(elphin%matsubara_temps)
            allocate(elphin%matsubara_temps(elphin%matsubara_ntemps))
            delta =(matsubara_temps_range(2)-matsubara_temps_range(1))/(elphin%matsubara_ntemps-1)
            do i=1,elphin%matsubara_ntemps
                elphin%matsubara_temps(i) = (i-1)*delta+matsubara_temps_range(1)
            enddo
        endif

       ! mu star
       call process_incar(incar, 'elph_mu_star', elphin%mu_star, ierr)
       ! maximum matsubara frequency
       call process_incar(incar, 'elph_max_matsubara', elphin%max_matsubara, ierr)
       ! number of frequency grid point to compute dos and eliashberg function
       call process_incar(incar, 'elph_nedos_ph', elphin%nedos_ph, ierr)
       ! maximum phonon energy (eV) to compute dos and eliashberg function
       call process_incar(incar, 'elph_emax_ph', elphin%emax_ph, ierr)
    end subroutine

    !> @brief Set default settings for the transport mode
    !> These can still be overwritten by explicit INCAR tags
    subroutine electron_phonon_reader_transport_meta(elphin)
        type(elph_settings),intent(inout) :: elphin
        ! local variable
        logical :: static

        elphin%run=.true.
        elphin%transport = .true.

        ! what is computed
        elphin%selfen_fan=.true.
        elphin%selfen_dfan=.false.
        elphin%selfen_dw=.false.
        elphin%selfen_fans=.false.
        elphin%selfen_dws=.false.

        elphin%selfen_imag_skip=.true.
        elphin%selfen_broad_tol=1e-4
        elphin%wf_redistribute=.true.

        ! select serta and mrta_lambda scattering approximation by default
        static=.false.
        call scattering_approx_from_str(elphin%selfen_approx, static, 'serta mrta_lambda')

        ! list of carrier densities for which to compute the scattering matrix elements
        elphin%selfen_ncarrier_den=13
        if (allocated(elphin%selfen_carrier_den)) deallocate(elphin%selfen_carrier_den)
        allocate(elphin%selfen_carrier_den(elphin%selfen_ncarrier_den))
        elphin%selfen_carrier_den(1:elphin%selfen_ncarrier_den)=[real(q):: &
            -1e+21_q, -1e+20_q, -1e+19_q, -1e+18_q, -1e+17_q, -1e+16_q,  &
            0.0_q, &
            1e+16_q, 1e+17_q, 1e+18_q, 1e+19_q, 1e+20_q, 1e+21_q]

        ! use tetrahedron method
        if (allocated(elphin%selfen_delta)) deallocate(elphin%selfen_delta)
        allocate(elphin%selfen_delta(1))
        elphin%selfen_delta=0.0_q
    end subroutine

    !> @brief Set default settings for the superconductivity mode
    subroutine electron_phonon_reader_superconductivity_meta(elphin)
        type(elph_settings),intent(inout) :: elphin

        elphin%superconductivity = .true.
    end subroutine

    !> Set default settings for the bandstructure renormalization mode
    !> These can still be overwritten by explicit INCAR tags
    subroutine electron_phonon_reader_renormalization_meta(elphin)
        type(elph_settings),intent(inout) :: elphin

        elphin%run=.true.
        ! what is computed
        elphin%selfen_fan=.true.
        elphin%selfen_dfan=.false.
        elphin%selfen_dw=.true.
        elphin%selfen_fans=.false.
        elphin%selfen_dws=.false.

        elphin%selfen_gaps=.true. ! select KS states that form the gap to be computed
        elphin%nbands = -2 ! tells VASP to compute all the bands (same as the number of plane-waves)

        ! use ia small finite broadening
        if (allocated(elphin%selfen_delta)) deallocate(elphin%selfen_delta)
        allocate(elphin%selfen_delta(1))
        elphin%selfen_delta=0.01_q
    end subroutine

    !> @brief Implements logic related to INCAR tags after all reading is done
    subroutine electron_phonon_reader_post_logic(elphin,xc,wave_interp_set,io)
        use string, only: str, lowercase
        use fock_glb, only: lhfcalc
        use wave_interpolate, only: wave_interpolator_settings
        use setexm, only: xc_info
        type(elph_settings), intent(inout) :: elphin
        type(xc_info),intent(in) :: xc
        type(wave_interpolator_settings), intent(inout) :: wave_interp_set
        type(in_struct) :: io

        if (elphin%run.and.lhfcalc) call vtutor%error( &
            'Electron-phonon calculations using hybrid functionals are currently not supported.')

#ifdef no_elph_superconductivity
        if (elphin%driver==elph_driver_ph) call vtutor%error( &
            'Phonon self-energy calculation due to electron-phonon coupling (ELPH_DRIVER=PH) is currently not supported')

        if (lowercase(elphin%mode)=='superconductivity') call vtutor%error( &
            'Superconductivity property calculations are currently not supported')
#else
        if (elphin%driver==elph_driver_ph) call vtutor%alert( &
            'Phonon self-energy calculation due to electron-phonon coupling (ELPH_DRIVER=PH) is currently not supported')

        if (lowercase(elphin%mode)=='superconductivity') call vtutor%alert( &
            'Superconductivity property calculations are currently not supported')
#endif

        if (elphin%driver==elph_driver_ph) then
            if ((elphin%selfen_band_start .gt. 0) .or. (elphin%selfen_band_stop .gt. 0)) then
                call vtutor%error( 'For a phonon SE calculation the values of' // &
                                   'ELPH_SELFEN_BAND_START, ELPH_SELFEN_BAND_STOP and ELPH_NBANDS_SUM must not be set in the INCAR.' // &
                                   ' Default values (1 and NBANDS) will be used')
            endif
            if ((size(elphin%nbands_sum) .gt. 1) .or. (elphin%nbands_sum(1) .gt. 0)) then
                call vtutor%error( 'For a phonon SE calculation the values of' // &
                                   'ELPH_SELFEN_BAND_START, ELPH_SELFEN_BAND_STOP and ELPH_NBANDS_SUM must not be set in the INCAR.' // &
                                   ' Default values (1 and NBANDS) will be used')
            endif
        endif

        if (elphin%prepare) then
            io%wrt_potential%total=.true.
            io%wrt_potential%paw=.true.
            if (xc%ldometagga.and.(xc%ltau.or.xc%lmu)) then
                io%wrt_potential%xcmu=.true.
            end if
        endif

        if (size(elphin%selfen_approx)==0) then
            call vtutor%error('ELPH_SCATTERING_APPROX was not specified or recognized. '//&
                              'Must be one or a combination of CRTA, SERTA, ERTA_TAU, ERTA_LAMBDA, MRTA_TAU or MRTA_LAMBDA')
        endif

        if (elphin%transport .and. any(elphin%selfen_temps < 1.0) .and. (elphin%ismear >= -1)) call vtutor%alert( &
            'You selected a smearing method (elph_ismear=' // str(elphin%ismear) // &
            ') and a temperature below 1K among the list of temperatures (elph_selfen_temps). ' // &
            'This will likely yield inaccurate results. We recommend setting elph_ismear=-15.')

        if ((elphin%transport_driver .ne. 1) .and. (elphin%transport_driver .ne. 2)) &
            call vtutor%error('ELPH_TRANSPORT_DRIVER must be 1 or 2')

        if (elphin%writemels) elphin%driver = elph_driver_mels

        if (elphin%driver==elph_driver_mels .and. elphin%wf_redistribute) then
            call vtutor%error('ELPH_WF_REDISTRIBUTE=.TRUE. is not supported in combination with ELPH_DRIVER=MELS')
        endif

        if (elphin%transport_dfermi_tol <= 0 .or. elphin%transport_dfermi_tol >= 1.0_q) &
            call vtutor%error('Invalid value for transport_dfermi_tol. It must be in the interval from ]0,1[ ')

        if (elphin%selfen_dfan.and.all(elphin%selfen_delta==0.0_q)) then
            call vtutor%error('Computation of the derivative of the Fan self-energy (ELPH_SELFEN_DFAN=.TRUE.) only implemented for smearing methods. '//&
                              'Please set ELPH_SELFEN_DELTA > 0.0')
        endif

        elphin%selfen_naccumulators = size(elphin%selfen_delta)*&
                                      size(elphin%nbands_sum)*&
                                      size(elphin%selfen_approx)
        elphin%selfen_naccumulators = elphin%selfen_naccumulators * &
                    max(1,elphin%selfen_ncarrier_per_cell)* &
                    max(1,elphin%selfen_ncarrier_den)* &
                    max(1,elphin%selfen_nmu)

        if (elphin%lgaddt.and.(.not.elphin%metagga)) then
            call vtutor%error('ELPH_DECOMPOSE=T is only implemented for meta-GGA functionals. Please remove T from ELPH_DECOMPOSE or use a meta-GGA functional.')
        endif

        ! HOTFIX for elph_nbands = 2
        ! TODO: Fix properly (better meta tags?)
        if (elphin%nbands == -2) then
            call vtutor%warning("ELPH_NBANDS=-2 requires exact diagonalization. " // &
                "Setting KPOINTS_OPT=2 internally.")
            elphin%wave_interp_set%kpoints_opt_mode = 2
        endif
   end subroutine

    !> @brief Read information related to the electron-phonon calculation
    subroutine electron_phonon_reader(incar,elphin,xc,wave_interp_set,io)
        use string, only: lowercase
        use setexm, only: xc_info
        use wave_interpolate, only: wave_interpolator_settings
        type(incar_file),intent(inout) :: incar
        type(elph_settings),intent(out) :: elphin
        type(xc_info),intent(in) :: xc
        type(wave_interpolator_settings),intent(inout) :: wave_interp_set
        type(in_struct),intent(inout) :: io

        logical :: lopen
        integer :: ierr
        character(len=:),allocatable :: elph_driver

        ! defaults values
        call electron_phonon_reader_general_default(elphin, xc, wave_interp_set)
        call electron_phonon_reader_phonon_default(elphin)
        call electron_phonon_reader_matrix_elements_default(elphin)
        call electron_phonon_reader_performance_default(elphin)
        call electron_phonon_reader_chemical_potential_default(elphin)
        call electron_phonon_reader_self_energy_default(elphin)
        call electron_phonon_reader_transport_default(elphin)
        call electron_phonon_reader_phonon_self_energy_default(elphin)
        call electron_phonon_reader_superconductivity_default(elphin)

        ! electron-phonon meta tag
        elph_driver = 'el'
        call process_incar(incar, 'elph_driver', elph_driver, ierr)
        select case (lowercase(elph_driver))
            case ('el')
                elphin%driver=elph_driver_el
            case ('ph')
                elphin%driver=elph_driver_ph
            case ('mels')
                elphin%driver=elph_driver_mels
            case default
                call vtutor%error('Unknown value '//trim(elph_driver)//'for elph_driver for electron-phonon computation. '//&
                                  'The possible values are: EL for electron self-energy and PH for phonon self-energy')
        end select

        ! transport meta tag
        elphin%mode=''
        call process_incar(incar, 'elph_mode', elphin%mode, ierr)
        select case (lowercase(elphin%mode))
            case('transport')
                call electron_phonon_reader_transport_meta(elphin) ! some of the self energy variables are changed
            case('renorm')
                call electron_phonon_reader_renormalization_meta(elphin) ! some of the self energy variables are changed
            case('superconductivity')
                call electron_phonon_reader_superconductivity_meta(elphin) ! some of the self energy variables are changed
            case('')
                ! do nothing, just use the defaults
            case default
                call vtutor%error('Unknown value '//trim(elphin%mode)//'for elph_mode for electron-phonon computation. '//&
                                  'The possible values are: TRANSPORT, RENORM or SUPERCONDUCTIVITY')
        end select
        ! phonon driver tag
        if (elphin%driver == elph_driver_ph) then
            call electron_phonon_reader_phonon_self_energy_meta(elphin)
        endif

        ! read in values
        call electron_phonon_reader_general_read(incar, elphin, io)
        call electron_phonon_reader_phonon_read(incar, elphin, io)
        call electron_phonon_reader_matrix_elements_read(incar, elphin, io)
        call electron_phonon_reader_performance_read(incar, elphin, io)
        call electron_phonon_reader_chemical_potential_read(incar, elphin, io)
        call electron_phonon_reader_self_energy_read(incar, elphin, io)
        call electron_phonon_reader_transport_read(incar, elphin, io)
        call electron_phonon_reader_phonon_self_energy_read(incar, elphin, io)
        call electron_phonon_reader_superconductivity_read(incar, elphin, io)

        ! apply checks and other logic on settings
        call electron_phonon_reader_post_logic(elphin, xc, wave_interp_set, io)
    end subroutine electron_phonon_reader

    !> Write the parsed settings to the OUTCAR.
    !> The INCAR values have interdependencies, i.e. one value might alter the value of another settings.
    !> The logic that implements this might change between different versions of vasp making it very difficult to verify which
    !> settings were actually used in the calculation. In this routine we write these settings to the OUTCAR similar to what
    !> is done for other INCAR variables in vasp.
    subroutine elph_write_settings(elphin,iu)
        type(elph_settings),intent(in) :: elphin
        integer,intent(in) :: iu
        !local variables
        integer :: i

        if (iu<0) return
        write(iu,'(a)')
        write(iu,'(a)') 'General:'
        write(iu,'(a,l)')    '  elph_run = ', elphin%run
        write(iu,'(a,i1)')   '  elph_driver = ', elphin%driver
        write(iu,'(a,l)')    '  elph_dry_run = ', elphin%dry_run
        write(iu,'(a,i6)')   '  elph_nbands = ', elphin%nbands
        write(iu,'(a,i6)')   '  elph_nedos = ', elphin%nedos
        write(iu,'(a,f8.3)') '  elph_kspacing = ', elphin%kspacing
        write(iu,'(a,f8.3)') '  encutlr = ', elphin%encutlr

        write(iu,'(a)')
        write(iu,'(a)') 'Phonon:'
        write(iu,'(a,l)') '  elph_ignore_imag_phonons = ', elphin%ignore_imag_phonons
        write(iu,'(a,i2)') '  ifc_asr = ', elphin%ifc_asr
        write(iu,'(a,i2)') '  ifc_lr = ', elphin%ifc_lr

        write(iu,'(a)')
        write(iu,'(a)') 'Matrix elements:'
        write(iu,'(a)') '                   v t d r p q s'
        write(iu,'(a,6(l,a))') '  elph_decompose = ', elphin%lgaddv,' ',&
                                                      elphin%lgaddt,' ',&
                                                      elphin%lgaddd,' ',&
                                                      elphin%lgaddr,' ',&
                                                      elphin%lgaddp,' ',&
                                                      elphin%lgaddq,' ',&
                                                      elphin%lgadds!, 'AE: vdpr  PS: vdq'
        write(iu,'(a,i2)') '  elph_lr = ', elphin%elph_lr
        write(iu,'(a)',advance='no') '  elph_writemels = '
        if (elphin%usehdf5mpiio) then
            write(iu,'(a)') 'hdf5mpiio'
        else
            write(iu,'(l)') elphin%writemels
        endif

        write(iu,'(a)')
        write(iu,'(a)') 'Performance:'
        write(iu,'(a,l)')    '  elph_wf_cache_prefill = ',elphin%wf_cache_prefill
        write(iu,'(a,f8.3)') '  elph_wf_cache_mb = ',elphin%wf_cache_mb
        write(iu,'(a,l)')    '  elph_use_pot_cache = ',elphin%use_pot_cache
        write(iu,'(a,l)')    '  elph_wf_redistribute = ',elphin%wf_redistribute
        write(iu,'(a,i4)')   '  elph_wf_redistribute_opt = ',elphin%wf_redistribute_opt
        write(iu,'(a,i4)')   '  elph_kpar = ',elphin%kpar
        write(iu,'(a,i4)')   '  elph_ncore = ',elphin%ncore
        write(iu,'(a,l)')    '  elph_useblas = ',elphin%useblas
        write(iu,'(a,l)')    '  elph_userecip = ',elphin%userecip
        write(iu,'(a,l)')    '  elph_rotateprojectors = ',elphin%rotateprojectors
        write(iu,'(a,i4)')   '  elph_wf_comm_opt = ',elphin%wf_comm_opt

        write(iu,'(a)')
        write(iu,'(a)') 'Chemical potential:'
        write(iu,'(a,i3)')     '  elph_ismear = ',elphin%ismear
        write(iu,'(a,i6)')     '  elph_fermi_nedos = ',elphin%fermi_nedos

        write(iu,'(a)') '  elph_selfen_temps ='
        do i=1,elphin%selfen_ntemps
            write(iu,'(a,f8.3)') '    ',elphin%selfen_temps(i)
        enddo

        write(iu,'(a)') '  elph_selfen_carrier_den ='
        do i=1,elphin%selfen_ncarrier_den
            write(iu,'(a,es10.3)') '    ',elphin%selfen_carrier_den(i)
        enddo

         write(iu,'(a)') '  elph_selfen_carrier_per_cell ='
        do i=1,elphin%selfen_ncarrier_per_cell
            write(iu,'(a,e10.3)') '    ',elphin%selfen_carrier_per_cell(i)
        enddo

        write(iu,'(a)') '  elph_selfen_mu ='
        do i=1,elphin%selfen_nmu
            write(iu,'(a,f8.3)') '    ',elphin%selfen_mu(i)
        enddo

        write(iu,'(a)')
        write(iu,'(a)') 'Self-energy:'
        write(iu,'(a,l)') '  elph_scattering_approx = '
        do i=1,size(elphin%selfen_approx)
            write(iu,'(i3,a,a)') i,' ',elph_selfen_approx_str(elphin%selfen_approx(i),verbose=.true.)
        enddo
        write(iu,'(a,l)') '  elph_selfen_fan = ',elphin%selfen_fan
        write(iu,'(a,l)') '  elph_selfen_dw = ',elphin%selfen_dw
        write(iu,'(a,l)') '  elph_selfen_fans = ',elphin%selfen_fans
        write(iu,'(a,l)') '  elph_selfen_dws = ',elphin%selfen_dws
        write(iu,'(a,i6)') '  elph_selfen_band_start = ',elphin%selfen_band_start
        write(iu,'(a,i6)') '  elph_selfen_band_stop = ',elphin%selfen_band_stop
        write(iu,'(a,i6)') '  elph_nbands_sum ='
        do i=1,size(elphin%nbands_sum,1)
            write(iu,'(i6)') elphin%nbands_sum(i)
        enddo
        write(iu,'(a)') '  elph_selfen_ikpt ='
        if (allocated(elphin%selfen_ikpt)) then
            do i=1,size(elphin%selfen_ikpt)
                write(iu,'(a,i6)') '    ',elphin%selfen_ikpt(i)
            enddo
        endif
        write(iu,'(a)') '  elph_selfen_kpts ='
        if (allocated(elphin%selfen_kpts)) then
            do i=1,size(elphin%selfen_kpts,2)
                write(iu,'(a,3f8.3)') '    ',elphin%selfen_kpts(:,i)
            enddo
        endif
        write(iu,'(a,i6)') '  elph_selfen_nw = ',elphin%selfen_nw
        write(iu,'(a,f8.3)') '  elph_selfen_wrange = ',elphin%selfen_wrange
        write(iu,'(a,l)') '  elph_selfen_gaps = ',elphin%selfen_gaps
        write(iu,'(a)') '  elph_selfen_delta ='
        do i=1,size(elphin%selfen_delta)
            write(iu,'(a,f8.3)') '    ',elphin%selfen_delta(i)
        enddo
        write(iu,'(a,f8.3,a,f8.3,a)') '  elph_selfen_energy_window = [',&
                                      elphin%selfen_energy_window(1),':',&
                                      elphin%selfen_energy_window(2),']'
        write(iu,'(a,l)') '  elph_selfen_g_skip = ',elphin%selfen_g_skip
        write(iu,'(a,l)') '  elph_selfen_imag_skip = ',elphin%selfen_imag_skip
        write(iu,'(a,e8.3)') '  elph_selfen_broad_tol = ',elphin%selfen_broad_tol


        write(iu,*)
        write(iu,'(a)') 'Transport:'
        write(iu,'(a,l)')    '  elph_transport = ',elphin%transport
        write(iu,'(a,i6)')   '  elph_transport_nedos_plot = ',elphin%transport_nedos_plot
        if (elphin%transport_emin_plot==huge(1.0_q)) then
            write(iu,'(a)') '  elph_transport_emin_plot = min(eigenvalues)'
        else
            write(iu,'(a,f8.3)') '  elph_transport_emin_plot = ',elphin%transport_emin_plot
        endif

        if (elphin%transport_emax_plot==huge(1.0_q)) then
            write(iu,'(a)') '  elph_transport_emax_plot = max(eigenvalues)'
        else
            write(iu,'(a,f8.3)') '  elph_transport_emax_plot = ',elphin%transport_emax_plot
        endif

        write(iu,'(a,e8.3)') '  elph_transport_dfermi_tol = ',elphin%transport_dfermi_tol
        write(iu,'(a,i3)')   '  elph_transport_driver = ',elphin%transport_driver
        write(iu,'(a,i6)')   '  elph_transport_nedos = ',elphin%transport_nedos
        write(iu,'(a,e8.3)') '  transport_relaxation_time = ',elphin%transport_relaxation_time
        write(iu,'(a,i3)')   '  elph_velocity_mode = ',elphin%velocity_mode

        write(iu,*)

    end subroutine

#ifdef VASP_HDF5
    !> @brief write information related to the k-points and symmetries
    subroutine vh5_write_elph_kpoints_regular(fileid,nrotk,igrpop,irot_fbz2ibz,indx_fbz2ibz)
        use vhdf5
        integer(HID_T) :: fileid
        integer :: nrotk !< Number of reciprocal symmetry operations
        integer :: igrpop(:,:,:) !< Integer reciprocal symmetry operations
        integer :: irot_fbz2ibz(:) !< Index of the rotation from IBZ to FBZ
        integer :: indx_fbz2ibz(:) !< Index from the FBZ to the IBZ
        ! local variables
        integer(HID_T) :: groupid

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, 'kpoints', groupid))

        VH5_CHECK(vh5_write(groupid,'nrotk',nrotk))
        VH5_CHECK(vh5_write(groupid,'irot_fbz2ibz',irot_fbz2ibz))
        VH5_CHECK(vh5_write(groupid,'indx_fbz2ibz',indx_fbz2ibz))
        VH5_CHECK(vh5_write(groupid,'igrpop',igrpop))

        ! close group
        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine vh5_write_elph_kpoints_regular

    !> @brief write lists of k-points
    subroutine vh5_write_elph_kpoints_list(fileid,wtkpt_k,vkpt_k,vkpt_kp)
        use vhdf5
        integer(HID_T) :: fileid
        real(q) :: wtkpt_k(:) !< weights of the k-points at k
        real(q) :: vkpt_k(:,:) !< Coordinates of the IBZ k-points at k
        real(q) :: vkpt_kp(:,:) !< Coordinates of the IBZ k-points at k'
        ! local variables
        integer(HID_T) :: groupid

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, 'kpoints', groupid))

        ! open group
        VH5_CHECK(vh5_write(groupid,'wtkpt_k',wtkpt_k))
        VH5_CHECK(vh5_write(groupid,'vkpt_k',vkpt_k))
        VH5_CHECK(vh5_write(groupid,'vkpt_kp',vkpt_kp))

        ! close group
        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine vh5_write_elph_kpoints_list

    subroutine vh5_write_elph_eigenvalues_list(fileid,eig_k,eig_kp)
        use vhdf5
        integer(HID_T) :: fileid
        real(q),intent(in) :: eig_k(:,:,:)
        real(q),intent(in) :: eig_kp(:,:,:)
        ! local variables
        integer(HID_T) :: groupid

        ! open group
        !VH5_CHECK(vh5_group_open_or_create(fileid, GRP_ELECTRON, groupid))
        ! keep compatiblity for now
        VH5_CHECK(vh5_group_open_or_create(fileid, 'matrix_elements', groupid))

        ! open group
        VH5_CHECK(vh5_write(groupid,'eigenvalues_k',eig_k))
        VH5_CHECK(vh5_write(groupid,'eigenvalues_kp',eig_kp))

        ! close group
        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine vh5_write_elph_eigenvalues_list
#endif

    !> Generate an explicit list of kpoints to write to the hdf5 file
    !> Only one node will allocate the eigenvalues and k-points arrays
    subroutine get_elph_kpoints_eigenvals_list(kpoints,fbz2ibz,kpoint_index,celtot,nkpts_k,nkpts_kp,wtkpt_k,vkpt_k,vkpt_kp,eig_k,eig_kp,lmaster)
        use mkpoints_struct_def, only: kpoints_struct, LineMode, ExplicitList
        use wave_rotate, only: fbz2ibz_stars, rotate_vkpt_irot
        type(kpoints_struct) :: kpoints !< structure containing a list of k-points
        type(fbz2ibz_stars) :: fbz2ibz
        integer,intent(in) :: kpoint_index(:) !< array with the k-points that are selected to be computed
        complex(q),intent(in) :: celtot(:,:,:)
        integer,intent(out) :: nkpts_k !< number of k-points at k (set for all ranks)
        integer,intent(out) :: nkpts_kp !< number of k-points at k' (set for all ranks)
        real(q),allocatable,intent(out) :: wtkpt_k(:) !< weights of the k-points at k (only set when lmaster=.true.)
        real(q),allocatable,intent(out) :: vkpt_k(:,:) !< k-point coordinates at k (only set when lmaster=.true.)
        real(q),allocatable,intent(out) :: vkpt_kp(:,:) !< k-point coordinates at k' (only set when lmaster=.true.)
        real(q),allocatable,intent(out) :: eig_k(:,:,:) !< eigenvalues list at k (only set when lmaster=.true.)
        real(q),allocatable,intent(out) :: eig_kp(:,:,:) !< eigenvalues list at k' (only set when lmaster=.true.)
        logical,intent(in) :: lmaster !< is this rank master?
        !local variables
        integer :: nbands, ispin, isp, indx, ik_ibz, ik
        integer :: irot
        nbands = size(celtot,1)
        ispin = size(celtot,3)

        if (kpoints%mode==LineMode.or.kpoints%mode==ExplicitList) then
            ! fix k and k' is from the kpoints structure
            nkpts_k = size(kpoint_index)
            nkpts_kp = kpoints%nkpts
            if (lmaster) then
                allocate(vkpt_kp(3,nkpts_kp))
                allocate(vkpt_k(3,nkpts_k))
                allocate(wtkpt_k(nkpts_k))
                vkpt_kp = kpoints%vkpt
                vkpt_k(:,:) = kpoints%vkpt(:,kpoint_index)
                wtkpt_k(:) = kpoints%wtkpt(kpoint_index)
                ! get eigenvalues for k
                allocate(eig_k(nbands,nkpts_k,ispin))
                do isp=1,ispin
                    do ik=1,size(kpoint_index)
                        ik_ibz = kpoint_index(ik)
                        eig_k(:,ik,isp) = real(celtot(:,ik_ibz,isp),q)
                    enddo
                enddo
                ! get eigenvalues for k'
                allocate(eig_kp(nbands,nkpts_kp,ispin))
                do isp=1,ispin
                    do ik=1,nkpts_kp
                        eig_kp(:,ik,isp) = real(celtot(:,ik,isp),q)
                    enddo
                enddo
            endif
            ! TODO: fix q, k is kpoints and k' is k-q
            !nkpts_k = kpoints%nkpts
            !nkpts_kp = kpoints%nkpts
            !allocate(vkpt_kp(3,nkpts_kp))
            !allocate(vkpt_k(3,nkpts_k))
            !vkpt_k = kpoints%vkpt
            !do ik=1,nkpts_kp
            !    vkpt_kp(:,ik) = vkpt_k(:,ik)-vqpt
            !enddo
        else
            ! regular k-point mesh
            nkpts_k = size(kpoint_index)
            nkpts_kp = fbz2ibz%nkfbz
            if (lmaster) then
                allocate(vkpt_kp(3,nkpts_kp))
                allocate(vkpt_k(3,nkpts_k))
                allocate(wtkpt_k(nkpts_k))
                vkpt_k(:,:) = kpoints%vkpt(:,kpoint_index)
                wtkpt_k(:) = kpoints%wtkpt(kpoint_index)
                do ik=1,fbz2ibz%nkfbz
                    irot = fbz2ibz%irot_fbz2ibz(ik)
                    indx = fbz2ibz%indx_fbz2ibz(ik)
                    vkpt_kp(:,ik) = rotate_vkpt_irot(kpoints%vkpt(:,indx),fbz2ibz%igrpop,irot)
                enddo
                ! get eigenvalues at k' by expanding from the IBZ to the FBZ
                allocate(eig_kp(nbands,nkpts_kp,ispin))
                do isp=1,ispin
                    do ik=1,fbz2ibz%nkfbz
                        !irot = fbz2ibz%irot_fbz2ibz(ik)
                        indx = fbz2ibz%indx_fbz2ibz(ik)
                        eig_kp(:,ik,isp) = real(celtot(:,indx,isp),q)
                    enddo
                enddo
                ! get eigenvalues at k by copying from the IBZ
                allocate(eig_k(nbands,nkpts_k,ispin))
                do isp=1,ispin
                    do ik=1,size(kpoint_index)
                        ik_ibz = kpoint_index(ik)
                        eig_k(:,ik,isp) = real(celtot(:,ik_ibz,isp),q)
                    enddo
                enddo
            endif
        endif
    end subroutine get_elph_kpoints_eigenvals_list

    !> Get number of k-points to be treated locally
    subroutine get_elph_nkpoints_local(wave_map,waverot,fbz2ibz,kpoints,nkpts_k,nkpts_k_local,nkpts_kp_local)
        use wave_rotate, only: wave_rotator, kpoint_star, fbz2ibz_stars, wave_rotator_get_star
        use wave_interpolate, only: wave_mapper
        use mkpoints_struct_def, only: kpoints_struct
        type(wave_mapper),intent(in) :: wave_map
        type(wave_rotator),intent(in) :: waverot
        type(fbz2ibz_stars),intent(in) :: fbz2ibz
        type(kpoints_struct),intent(in) :: kpoints
        integer,intent(in) :: nkpts_k
        integer,intent(out) :: nkpts_k_local
        integer,intent(out) :: nkpts_kp_local
        ! local variables
        integer :: ikp, ikp_ibz
        type(kpoint_star) :: star
        ! count star size of k-points that are owned locally to have correct progress bar
        nkpts_kp_local = 0
        do ikp=1,wave_map%nkpts_batch
            ikp_ibz = wave_map%kpoints_index(ikp+wave_map%nkpts_orig)
            ! loop over points in the star of k'
            call wave_rotator_get_star(waverot,kpoints%vkpt(:,ikp_ibz),star)
            nkpts_kp_local = nkpts_kp_local + fbz2ibz%star_size(ikp_ibz)
        enddo
        nkpts_k_local = nkpts_k
    end subroutine get_elph_nkpoints_local

#ifdef VASP_HDF5
    subroutine vh5_write_dvdu_avg(fileid,natoms,ncdij,dvdu_avg,dvdu_sr_avg,dvdu_lr_avg,&
                                    ddijdu_avg,ddijdu_sr_avg,ddijdu_lr_avg,nqpts,vqpts)
        use vhdf5
        integer(HID_T) :: fileid
        integer :: natoms !< Number of atoms
        integer :: ncdij
        integer :: nqpts !< Number of q-points
        complex(q),intent(in) :: dvdu_avg(3,natoms,ncdij,nqpts)
        complex(q),intent(in) :: dvdu_sr_avg(3,natoms,ncdij,nqpts)
        complex(q),intent(in) :: dvdu_lr_avg(3,natoms,ncdij,nqpts)
        complex(q),intent(in) :: ddijdu_avg(3,natoms,ncdij,nqpts)
        complex(q),intent(in) :: ddijdu_sr_avg(3,natoms,ncdij,nqpts)
        complex(q),intent(in) :: ddijdu_lr_avg(3,natoms,ncdij,nqpts)
        real(q) :: vqpts(3,nqpts)
        ! local variables
        integer(HID_T) :: groupid

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, 'results/'//trim(GRP_ELPHON)//'/dvdu', groupid))

        ! write some dimensions to file
        VH5_CHECK(vh5_write(groupid, 'nqpts',  nqpts))
        VH5_CHECK(vh5_write(groupid, 'natoms', natoms))

        ! write the coordinates of the qpoints
        VH5_CHECK(vh5_write(groupid, 'qpoints', vqpts))

        ! Write the average of the potential
        VH5_CHECK(vh5_write(groupid, 'dvdu_avg', dvdu_avg))
        VH5_CHECK(vh5_write(groupid, 'dvdu_sr_avg', dvdu_sr_avg))
        VH5_CHECK(vh5_write(groupid, 'dvdu_lr_avg', dvdu_lr_avg))
        VH5_CHECK(vh5_write(groupid, 'ddijdu_avg', ddijdu_avg))
        VH5_CHECK(vh5_write(groupid, 'ddijdu_sr_avg', ddijdu_sr_avg))
        VH5_CHECK(vh5_write(groupid, 'ddijdu_lr_avg', ddijdu_lr_avg))

        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine vh5_write_dvdu_avg
#endif

    !> Write header indicating the start of an electron-phonon calculation
    subroutine elph_header_write(io)
        type(in_struct),intent(in) :: io
        if (io%iu0>=0) write(io%iu0,*) 'electron-phonon calculation'
        if (io%iu6>=0) then
           write(io%iu6,*) ''
           write(io%iu6,*) 'electron-phonon calculation'
           write(io%iu6,*) '==========================='
        endif
    end subroutine

    !> Write footer indicating the end of an electron-phonon calculation
    subroutine elph_footer_write(io)
        type(in_struct),intent(in) :: io
        if (io%iu6>0) then
           write(io%iu6,'(a)') '------------------------ end of electron-phonon driver reached -----------------------------------------'
           write(io%iu6,*)
        endif
    end subroutine

    !> Write information about the electron-phonon computation
    subroutine elph_wf_write(fft_mesh,ispin,nrspinors,band_start_kp,band_stop_kp,nbands_kp,&
                             band_start_k,bands_stop_k,nbands_k,elph_nbands,iu)
        type(elph_settings) :: elph_set
        integer,intent(in) :: fft_mesh(3)
        integer,intent(in) :: ispin
        integer,intent(in) :: nrspinors
        integer,intent(in) :: band_start_kp
        integer,intent(in) :: band_stop_kp
        integer,intent(in) :: nbands_kp
        integer,intent(in) :: band_start_k
        integer,intent(in) :: bands_stop_k
        integer,intent(in) :: nbands_k
        integer,intent(in) :: elph_nbands
        integer :: iu
        if (iu>=0) then
            write(iu,*) "============== electron-phonon wfs ================="
            write(iu,*) "fft_mesh:      ",  fft_mesh
            write(iu,*) wf_spin_to_string(ispin,nrspinors)
            write(iu,*) "band_start_k': ", band_start_kp
            write(iu,*) "band_stop_k':  ", band_stop_kp
            write(iu,*) "nbands_k':     ", nbands_kp
            write(iu,*) "band_start_k:  ", band_start_k
            write(iu,*) "band_stop_k:   ", bands_stop_k
            write(iu,*) "nbands_k:      ", nbands_k
            write(iu,*) "elph_nbands:   ", elph_nbands
            write(iu,*) "===================================================="
        endif
    end subroutine elph_wf_write

    !> From the value of ispin and nrspinors get a string identifying wether its a
    !> spin-unpolarized, spin-polarized or noncollinear calculation
    function wf_spin_to_string(ispin,nrspinors) result(string)
        integer,intent(in) :: ispin
        integer,intent(in) :: nrspinors
        character(len=:), allocatable :: string
        string = "unknown"
        if (ispin==1.and.nrspinors==1) then
            string = pot_spin_to_string(1)
        endif
        if (ispin==2.and.nrspinors==1) then
            string = pot_spin_to_string(2)
        endif
        if (ispin==1.and.nrspinors==2) then
            string = pot_spin_to_string(4)
        endif
    end function wf_spin_to_string

    !> From the spinor dimension of the potential get a string identifying wether its a
    !> spin-unpolarized, spin-polarized or noncollinear calculation
    function pot_spin_to_string(ncdij) result(string)
        integer,intent(in) :: ncdij
        character(len=:), allocatable :: string
        select case(ncdij)
        case(1)
            string = "spin-unpolarized (ispin=1, nrspinors=1, ncdij=1)"
        case(2)
            string = "spin-polarized (ispin=2, nrspinors=1, ncdij=2)"
        case(4)
            string = "noncollinear (ispin=1, nrspinors=2, ncdij=4)"
        end select
    end function pot_spin_to_string

    !> write a string corresponding to the option to use slow fourier transform
    function lr_sft_string(lr_sft) result(string)
        integer,intent(in) :: lr_sft
        character(len=:), allocatable :: string
        select case(lr_sft)
        case(-1)
            string = "no"
        case(0)
            string = "heuristic"
        case(1)
            string = "yes"
        case default
            call vtutor%error('Unknown option for EPLH_LR_SFT. Must be -1, 0 or 1.')
        end select
    end function lr_sft_string

    !> logic to decide wether the force constants and potential should be read
    function elph_read_ifc_pot(elph_set)
        type(elph_settings),intent(in) :: elph_set
        logical :: elph_read_ifc_pot
        elph_read_ifc_pot = .true.
        if (elph_set%driver==elph_driver_el.and.scattering_approx_all_crta(elph_set)) elph_read_ifc_pot = .false.
    end function elph_read_ifc_pot

    !> return true if all scattering approximations are CRTA
    function scattering_approx_all_crta(elph_set) result(all_crta)
        type(elph_settings),intent(in) :: elph_set
        logical :: all_crta
        all_crta = all(elph_set%selfen_approx(:)%scattering_approx==elph_scattering_approx_crta)
    end function scattering_approx_all_crta

    !> show progress of an electron-phonon calculation
    subroutine elph_progress(tag,i,ni,iu,step)
        character(len=*) :: tag !< name of the progress bar
        integer,intent(in) :: i !< current element
        integer,intent(in) :: ni !< total number of elements
        integer,intent(in) :: iu !< unit in which to report progress
        integer,optional,intent(in) :: step !< report every step
        ! local variables
        integer :: my_step
        my_step = 100 ! by default report only every 100 steps
        if (present(step)) my_step = step
        if (iu<0) return
        if (i==1.or.mod(i,my_step)==0) then
            write(iu,'(a,a,i8,a,i8,a,i8,a)') tag,' [',i,' - ',min(i+my_step-1,ni),' /',ni,']'
        endif
    end subroutine elph_progress

    !> Get max energy from celtot array
    function celtot_max(celtot,nb_totk) result(emax)
        complex(q) :: celtot(:,:,:)
        integer :: nb_totk(:,:)
        real(q) :: emax
        !local variables
        integer :: ispin, ikibz
        integer :: nspin, nkpoints
        integer :: ibandmax
        nspin=size(celtot,3)
        nkpoints=size(celtot,2)
        emax=real(celtot(1,1,1),q)
        do ispin=1,nspin
            do ikibz=1,nkpoints
                ibandmax = nb_totk(ikibz,ispin)
                emax = max(emax,real(celtot(ibandmax,ikibz,ispin),q))
            enddo
        enddo
    end function celtot_max

    !> @brief Function to compute the Fermi-Dirac distribution
    pure function fermi_dirac(e,temp) result(f)
        real(q),intent(in) :: e !< energy of the state minus the fermi energy
        real(q),intent(in) :: temp !< temperature of the state
        real(q) :: f !< fermi occupation
        !local variables
        real(q) :: arg
        if (temp < 1e-12) then
            ! this is just a heaviside function
            if (e<0) then
                f=1.0_q
            else
                f=0.0_q
            endif
        else
            arg = e/temp
            if (arg>FERMI_MAX_ARG) then
                f=0.0_q
            else if (arg<-FERMI_MAX_ARG) then
                f=1.0_q
            else
                f=1.0_q/(exp(arg)+1.0_q)
            endif
        endif
    end function fermi_dirac

    !> @brief function to compute the derivative w.r.t. energy of the fermi-dirac distrubtion
    pure function dfermi_dirac(e,temp) result(df)
        real(q),intent(in) :: e !< energy of the state
        real(q),intent(in) :: temp !< temperature of the state
        real(q) :: df !< derivative of the occupation w.r.t energy  --- > chemical potential I guess, otherwise there is a minus ?
        !local variables
        real(q) :: arg, exparg
        if (temp < 1e-12) then
            ! the result is a delta function in this case but it should be handled
            ! depending on the context, not here
            df=0.0_q
        else
            arg = e/temp
            if (arg>FERMI_MAX_ARG.or.arg<-FERMI_MAX_ARG) then
                df=0.0_q
            else
                exparg=exp(arg)
                df=exparg/(exparg+1.0_q)**2/temp
            endif
        endif
    end function dfermi_dirac

    !> @brief Compute the Bose-Einstein distribution
    pure function bose_einstein(e,temp) result(n)
        real(q),intent(in) :: e !< energy of the state
        real(q),intent(in) :: temp !< temperature of the state
        real(q) :: n !< occupation
        !local variables
        real(q) :: arg
        if (temp < 1e-12) then
            ! in this case the occupation should be huge for 'e' close to zero
            ! but this limiting case will not be handled here
            if (abs(e)<1e-12) then
                n = 0.0_q
            else
                n = temp/e
            endif
        else
            arg = e/temp
            if (arg>BOSE_MAX_ARG) then
                n = 0.0_q
            else
                n = 1.0_q/(exp(arg)-1.0_q)
            endif
        endif
    end function bose_einstein

    !> Compute min and max energy range on which we should evaluate the transport function
    !> By default use dfermi_tol and smearing to choose the energy range (based on the meaning of dfermi_tol).
    !> If alternatively emin_in and emax_in are specified by the user then
    !> emin=emin_in
    !> emax=emax_in
    subroutine transport_emin_emax(mu,smearing,dfermi_tol, emin_in, emax_in, emin,emax)
        real(q),intent(in) :: mu !< chemical potential. center the energy grid arounf this value
        real(q),intent(in) :: smearing
        real(q),intent(in) :: dfermi_tol !< percentage of the integral that is excluded
        real(q),intent(in) :: emin_in !< min energy range specified be the user, by default is huge(1.0_q)
        real(q),intent(in) :: emax_in !< max energy range specified be the user, by default is huge(1.0_q)
        real(q),intent(out) :: emin !< min energy at which the transport function will be evaluated for integration
        real(q),intent(out) :: emax !< max energy at whcih the transport function will be evaluated for integration
        ! local variables
        real(q) :: y
        !in this definition dfermi_tol is in [0,1[ and the results improve as we approach 1
        !y = -(dfermi_tol+1)/(dfermi_tol-1)

        !in this definition dfermi_tol is in ]0,1] and the results improve as we approach 0
        !this definition is related to the one above by dfermi_tol' = 1-dfermi_tol
        y = (2-dfermi_tol)/(dfermi_tol)
        emin = mu-smearing*log(y)
        emax = mu+smearing*log(y)

        if (abs(emax_in-emin_in) > 1e-6_q) then
            emin=emin_in
            emax=emax_in
        endif
    end subroutine transport_emin_emax

    !> Generate a linear grid with n energy points between emin and emax
    function linspace(xmin,xmax,n) result(x)
        real(q),intent(in) :: xmin
        real(q),intent(in) :: xmax
        integer,intent(in) :: n
        real(q),allocatable :: x(:)
        !local variables
        real(q) :: dx
        integer :: i
        allocate(x(n))
        dx = (xmax-xmin)/(n-1)
        do i=1,n
            x(i) = xmin+dx*(i-1)
        enddo
    end function linspace

  !> @brief Find index in array with value
  !> In Fortran 2008 revision an intrinsic function was introduced,
  !> but older compilers (for examlple gcc7) do not support it.
  function find_value(array,val) result(idx)
       integer :: array(:)
       integer :: val
       ! local variables
       integer :: idx
       do idx=1,size(array,1)
           if(array(idx)==val) return
       enddo
       idx = 0
  end function find_value

  !> @brief Transform the g matrix elements from cartesian to mode basis
  subroutine get_g_kkp_modes(natoms_pc,gcart_kkp,e_q,gmode_kkp)
      integer,intent(in) :: natoms_pc
      complex(q), intent(in) :: gcart_kkp(:,:,:,:)
      complex(q), intent(in) :: e_q(:,:)
      complex(q) :: gmode_kkp(:,:,:)
      integer i0,i1,i2,ind

      gmode_kkp=cmplx(0.0_q,0.0_q,q)

      do i0=1,3*natoms_pc
         ind=0
         do i1=1,natoms_pc
            do i2=1,3
               ind=ind+1
               gmode_kkp(:,:,i0)=gmode_kkp(:,:,i0)+gcart_kkp(:,:,i2,i1)*e_q(ind,i0)
            enddo
         enddo
      enddo

  end subroutine get_g_kkp_modes

  !> get global index of the local bands from wavedes
  subroutine wdes_get_global_band_index(wdes,band_idx)
      use wave_struct_def, only: wavedes
      type(wavedes),intent(in) :: wdes
      integer,allocatable,intent(out) :: band_idx(:)
      call get_global_band_index(wdes%nbands,wdes%nb_tot,wdes%nb_par,wdes%nb_low,band_idx)
  end subroutine wdes_get_global_band_index

  !> get global index of local bands from low level elements of wavedes
  subroutine get_global_band_index(nbands,nb_tot,nb_par,nb_low,band_idx)
      use wave, only: nb_local_low
      integer,intent(in) :: nbands
      integer,intent(in) :: nb_tot
      integer,intent(in) :: nb_par
      integer,intent(in) :: nb_low
      integer,allocatable,intent(out) :: band_idx(:)
      ! local variables
      integer :: nb, nb_global
      allocate(band_idx(nbands))
      do nb_global=1,nb_tot
         nb = nb_local_low(nb_global,nb_par,nb_low)
         if (nb/=0) band_idx(nb)=nb_global
      enddo ! loop over bands
  end subroutine get_global_band_index

end module elphon_base
