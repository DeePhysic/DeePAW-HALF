#include "symbol.inc"

! This module contains the structure and routines required for the computation of the
! electron phonon matrix elements in the PAW basis
module elphon_mels_pw
  use prec
  use mpimy, ONLY: communic
  use base, ONLY : in_struct, info_struct, symmetry
  use mkpoints_struct_def, only : kpoints_struct
  use poscar_struct_def, only: latt, type_info
  use pseudo_struct_def, only: potcar
  use mgrid_struct_def, only: grid_3d
  use wave_struct_def, only : wavespin, wavefuna, wavedes, wavedes1
  use wave_interpolate, only : wave_mapper
  use wave_rotate, only: wave_rotator, fbz2ibz_stars
  use wave_windower, only : wave_window, wave_cache
  use wave_interpolate, only: wave_interpolator
  use elphon_potential_struct, only : elphon_pot, phonon_ifc
  use elphon_base, only : elph_caches_t, elph_settings
  implicit none

  !> Datastructure for the computation of electron-phonon matrix elements
  !> in the PAW basis
  type elph_mels_pw_type
      ! incar settings used to generate the matrix elements
      type(elph_settings) :: elph_set
      ! dimensions for the computation of the electron-phonon matrix elements
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
                          !< ncdij=4 non-collinear calculation
      integer :: npw !< local number of grid points for fft = ELPH_SET%NPW=GRID%RL%NP (ASK_HM)
      integer :: npw_local !< storage needed for in-place FFT (grid%mplwv)
      integer :: npr_local !< local number of real-space grid points for fft = ELPH_SET%NP=GRID%RL%NP
      integer :: ispin !< number of spins
      integer :: nproj !< number of projectors
      integer :: nproj_tot !< number of projectors
      integer :: nrspinors !< number of spinors (1 for collinear, 2 for non collinear)
      integer :: ntyp !< total number of atom types
      integer :: ngdim_orig !< number of g-vectors in the original wf
      integer,allocatable :: nityp(:) !< number of ions per type (ntyp)
      integer,allocatable :: lmmax(:) !< max lm channel per atom type (ntyp)
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
      type(potcar),allocatable :: p(:) ! set in init
      type(latt) :: latt_ini ! set in init
      type(latt) :: latt_cur ! set in init
      type(kpoints_struct) :: kpoints_orig
      type(fbz2ibz_stars) :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
      type(info_struct) :: info ! set in init
      type(type_info) :: t_info ! set in init
      type(in_struct) :: io
      type(grid_3d) :: grid_orig !< original grid comming from SCF calculation
      type(grid_3d) :: grid !< grid to be used in electron-phonon calculation
                            !< the number of local columns might change because k-point set also changes
      type(symmetry) :: symm
      ! wfs
      type(wavespin) :: w_kbatch !< container for the k-distributed wfs
      type(wavedes) :: wdes_kbatch !< descriptor of the k-distributed wfs
      type(wavespin) :: w_dense_ibz !< container for eigenvalues in the dense ibz (not distributed)
      type(wavedes) :: wdes_dense_ibz !< descriptor of the w_dense_ibz
      type(wave_interpolator) :: wave_interp !< instance of wave interpolator used to compute WFs at arbitrary k-points
      type(wave_rotator) :: waverot !< dataype to rotate wfs
      type(wave_window) :: wavwin_cw !< wave window for the plane-wave part of the orbitals
      type(wave_window) :: wavwin_cproj !< wave window for the projections part of the orbitals
      type(wave_cache) :: cache_cw !< wave cache for the plane-wave part of the orbitals
      type(wave_cache) :: cache_cproj !< wave cache for the projections part of the orbitals
      type(wave_mapper) :: wave_map !< mapping of the k-points in the k-point batch to the IBZ
      integer,allocatable :: ikp_fbz_list(:) !< indexed of the k' points present locally
      integer :: workload_mode !< mode to compute the workload of each cpu
      ! electron-phonon potential
      complex(q),allocatable :: dvdu_kkp(:,:,:,:) !< CACHES electron-phonon potential at iq_cache (npw_local,3,natoms_pc,ncdij)
      complex(q),allocatable :: dmudu_kkp(:,:,:,:) !< CACHES derivative of kinetic energy density at iq_cache (npw_local,3,natoms_pc,ncdij)
      complex(q),allocatable :: ddijdu_kkp(:,:,:,:,:,:) !< CACHES electron-phonon potential at iq_cache (lmmax,lmmax,natoms_sc,3,natoms_pc_local,ncdij)
      real(q) :: vqpt(3) !< CACHES coordinates of the q-point where the potential was interpolated
      type (elphon_pot) :: elph_pot
      ! phonons
      real(q),allocatable :: phonon_freqs_ibz(:,:) !< store the phonon frequencies in the IBZ
      ! projections
      GDEF,allocatable :: dprojk(:,:,:,:) !CACHES
      GDEF,allocatable :: dprojkp(:,:,:,:) !CACHES
      GDEF,allocatable :: dproj_kbatch(:,:,:,:,:,:) !< contains the projections and its derivatives for the local batch of k-points
      real(q),allocatable :: velocity(:,:,:,:)
      ! caches
      type (elph_caches_t) :: elph_caches

      integer :: iq_cache !< store the index of the q-point that was used last
  end type elph_mels_pw_type

  type wf_redistribution
      !> original wave_map
      type(wave_mapper) :: wave_map_orig
      !> workload before redistribution
      integer,allocatable :: workload_orig(:)
      !> original number of k-points per batch
      integer,allocatable :: nkpts_batch_orig(:)
      !> workload per k-point (does not change with redistribution)
      integer,allocatable :: workload_ikp(:)
      !> wave_map after redistributing the workload
      type(wave_mapper) :: wave_map
      !> workload after redistribution
      integer,allocatable :: workload(:)
      !> number of k-points per batch
      integer,allocatable :: nkpts_batch(:)
  end type wf_redistribution

  contains

  !> get information about the grid in which we want to compute the electron-phonon matrix elements
  subroutine elph_mels_pw_init(self,elph_set,wdes,lmdim,grid,ngdim_orig,latt_ini,latt_cur,info,p,symm,t_info,io,comm_kbatch,comm_global)
      use string, only: str
      use elphon_base, only: elph_caches_init
      type(elph_mels_pw_type) :: self
      type(elph_settings),intent(in) :: elph_set
      type(wavedes),intent(in) :: wdes
      integer,intent(in) :: lmdim
      type(grid_3d),intent(in) :: grid
      integer,intent(in) :: ngdim_orig
      type(latt),intent(in) :: latt_ini
      type(latt),intent(in) :: latt_cur
      type(info_struct),intent(in) :: info
      type(potcar),intent(in) :: p(:)
      type(symmetry),intent(in) :: symm
      type(type_info),intent(in) :: t_info
      type(in_struct),intent(in) :: io
      type(communic),intent(in) :: comm_kbatch
      type(communic),intent(in) :: comm_global
      ! local variables
      integer :: nt

      ! set some important dimensions based on the wavefunction descritor
      self%ispin=wdes%ispin
      self%lmdim=lmdim
      self%nproj=wdes%npro/wdes%nrspinors
      self%nproj_tot=wdes%npro_tot/wdes%nrspinors
      self%nbands=wdes%nb_tot
      self%nbands_local=wdes%nbands
      self%nbands_k=wdes%nb_tot
      self%nbands_kp=wdes%nb_tot
      self%wf_ncdij=wdes%ncdij
      self%nrspinors=wdes%nrspinors
      self%npw=grid%ngx*grid%ngy*grid%ngz !grid%rl%np
      self%npw_local=grid%mplwv
      self%npr_local=grid%rl%np
      self%fft_mesh=[grid%ngx,grid%ngy,grid%ngz]
      self%ngdim_orig = ngdim_orig

      self%ntyp = t_info%ntyp
      allocate(self%nityp(self%ntyp))
      self%nityp = t_info%nityp
      allocate(self%lmmax(self%ntyp))
      do nt=1,self%ntyp
         self%lmmax(nt) = p(nt)%lmmax
      enddo

      ! set communicators
      self%comm_kbatch = comm_kbatch
      self%comm_global = comm_global

      ! keep pointers to common structures
      self%elph_set = elph_set
      self%t_info = t_info
      self%io = io
      self%latt_cur = latt_cur
      self%latt_ini = latt_ini
      self%info = info
      self%grid_orig = grid
      self%grid = grid
      self%wdes_dense_ibz = wdes
      self%symm = symm
      allocate(self%p(size(p)))
      self%p(:) = p(:)
      ! I set it based on the number of atoms in the POSCAR file read by vasp
      ! but this must be consistent with phelel_params.hdf5 which will be checked later
      self%natoms_pc = sum(t_info%nityp)
      self%natoms_pc_local = sum(wdes%nityp)

      self%iq_cache=-1
  end subroutine elph_mels_pw_init

  !> Read dimnensions from hdf5 file
  subroutine elph_mels_pw_read_dimensions(self,phonon_natoms_sc,phonon_det_p2s,filename)
      use tutor, only: vtutor
      use string, only: str
      use wap, only: read_phelel_dimensions_hdf5
      type(elph_mels_pw_type),intent(inout) :: self
      integer,intent(out) :: phonon_natoms_sc !< number of atoms in the supercell for phonon calculation
      integer,intent(out) :: phonon_det_p2s !< det_p2s for phonon calculation
      character(len=*),intent(in) :: filename
      ! local variables
      integer :: natoms_pc_check
      real(q) :: lattice_pc(3,3) ! dummy not used

      ! Read the potential
      call read_phelel_dimensions_hdf5(self%ngridpoints,self%det_p2s,self%pot_ncdij,&
                                       self%lmdim,self%natoms_pc,self%natoms_sc,&
                                       self%fft_mesh,lattice_pc,&
                                       phonon_natoms_sc,phonon_det_p2s, filename)

      ! check if number of atoms in the primitive cell is consitent with the current calculation
      if (self%natoms_pc /= sum(self%t_info%nityp)) then
          call vtutor%error('Number of atoms in the primitive cell read from phelel_params.hdf5 '//&
                            'is not consistent with the POSCAR file.')
      endif

      ! MPI sanity check
      self%natoms_pc_local = sum(self%wdes_dense_ibz%nityp)
      natoms_pc_check = self%natoms_pc_local
      CALLMPI( m_sum_i( self%wdes_dense_ibz%comm_inb, natoms_pc_check, 1))
      if (natoms_pc_check/=self%natoms_pc) then
         call vtutor%bug('sum of local atoms '//str(natoms_pc_check)//&
              ' is different from the total number of atoms '//str(self%natoms_pc),__FILE__,__LINE__)
      endif

  end subroutine elph_mels_pw_read_dimensions

  !> Set range of bands for which to compute the electron-phonon matrix elements
  !> We can choose start and stop bands individually for k and k'
  subroutine elph_mels_pw_set_nbandsk(self,start_k,stop_k)
      type(elph_mels_pw_type) :: self
      integer,intent(in) :: start_k
      integer,intent(in) :: stop_k
      ! bands for KS states in k
      self%band_start_k=start_k
      self%band_stop_k=stop_k
      self%nbands_k=self%band_stop_k-self%band_start_k+1
  end subroutine elph_mels_pw_set_nbandsk

  !> Set range of bands for which to compute the electron-phonon matrix elements at k'
  subroutine elph_mels_pw_set_nbandskp(self,start_kp,stop_kp)
      type(elph_mels_pw_type) :: self
      integer,intent(in) :: start_kp
      integer,intent(in) :: stop_kp
      ! bands for KS states in k'
      self%band_start_kp=start_kp
      self%band_stop_kp=stop_kp
      self%nbands_kp=self%band_stop_kp-self%band_start_kp+1
  end subroutine elph_mels_pw_set_nbandskp

  !> write information about WFs
  subroutine elph_mels_pw_wf_write(self)
      use elphon_base, ONLY: elph_wf_write
      type(elph_mels_pw_type) :: self
      call elph_wf_write(self%fft_mesh,self%ispin,self%nrspinors,&
                         self%band_start_kp,self%band_stop_kp,self%nbands_kp,&
                         self%band_start_k,self%band_stop_k,self%nbands_k,&
                         self%elph_set%nbands,self%io%iu0)
      call elph_wf_write(self%fft_mesh,self%ispin,self%nrspinors,&
                         self%band_start_kp,self%band_stop_kp,self%nbands_kp,&
                         self%band_start_k,self%band_stop_k,self%nbands_k,&
                         self%elph_set%nbands,self%io%iu6)
  end subroutine elph_mels_pw_wf_write

  !> Initialize temporary work arrays needed for the computation of electron-phonon matrix elements
  subroutine elph_mels_pw_cache_init(self)
      use ini, only: register_allocate
      use elphon_base, only: elph_caches_init
      use elphon_potential, ONLY: elphon_pot_transfer_ncdij
      use tutor, only: vtutor
      type(elph_mels_pw_type) :: self
      ! local variables
      integer :: ierr

      call elph_caches_init(self%elph_caches,self%nproj_tot,self%nbands,self%npr_local,self%nbands_k,self%nrspinors)

      call elphon_pot_transfer_ncdij(self%elph_pot,self%wf_ncdij,self%elph_caches%ccdij,&
                                     self%elph_caches%ccqij,self%elph_caches%ccdrij)

      ! do allocations
      allocate(self%dprojk(self%nproj,self%nbands_k,0:3,self%nrspinors))
      allocate(self%dprojkp(self%nproj,self%nbands_kp,0:3,self%nrspinors))
      ! allocate the electron-phonon potential folded from supercell to primitive cell
      allocate(self%dvdu_kkp(self%npr_local,3,self%natoms_pc,self%wf_ncdij))
      allocate(self%dmudu_kkp(self%npr_local,3,self%natoms_pc,self%wf_ncdij))
      allocate(self%ddijdu_kkp(self%lmdim,self%lmdim,self%natoms_pc,3,self%natoms_pc_local,self%wf_ncdij))
  end subroutine elph_mels_pw_cache_init

  !> Free temporary arrays needed for the computaiton of the electron-phonon matrix elements
  subroutine elph_mels_pw_cache_free(self)
      use ini, only: deregister_allocate
      use elphon_base, only: elph_caches_free
      type(elph_mels_pw_type) :: self

      call elph_caches_free(self%elph_caches)

      deallocate(self%dvdu_kkp)
      deallocate(self%dmudu_kkp)
      deallocate(self%ddijdu_kkp)
  end subroutine elph_mels_pw_cache_free

  !> Initialize windows to access k-point information from other nodes
  subroutine elph_mels_pw_windows_init(self)
      use elphon_common, only: create_wave_windows
      type(elph_mels_pw_type) :: self
      call create_wave_windows(self%wavwin_cw,self%wavwin_cproj,&
                               self%nproj,self%nrspinors,self%wdes_dense_ibz%nkpts,&
                               self%wave_map,self%w_kbatch,self%dproj_kbatch,self%wdes_kbatch,self%comm_global)
  end subroutine elph_mels_pw_windows_init

  subroutine elph_mels_pw_wave_caches_init(self)
      use elphon_common, only: create_wave_caches
      type(elph_mels_pw_type) :: self
      call create_wave_caches(self%cache_cw,self%cache_cproj,&
                              self%elph_set%wf_cache_mb,self%nproj,self%nproj_tot,self%nrspinors,self%wdes_dense_ibz%nkpts,&
                              self%wdes_kbatch)
  end subroutine elph_mels_pw_wave_caches_init

  !> redistribute the k-points of the wfs according to the work that needs to be performed
  subroutine elph_mels_pw_wf_redistribute(self,elph_accumulators,kpoints)
      use elphon_base, only: elph_workload_mode_ibz, elph_workload_mode_fbz, &
                             elph_comm_twosided, elph_comm_onesided, elph_driver_ph
      use elphon_accumulators_high, only: elph_accumulators_selfen
      use tutor, only: vtutor
      use string, only: str
      type(elph_mels_pw_type) :: self
      type(elph_accumulators_selfen) :: elph_accumulators
      type(kpoints_struct) :: kpoints
      !local variables
      type(wf_redistribution) :: wf_redis

      if (.not.self%elph_set%wf_redistribute) return

      PROFILING_START('wf_redistribute_driver')
      call wf_redistribute_init(wf_redis,kpoints%nkpts,self%comm_kbatch)
      call wave_map_get_ikpfbz(self%wave_map,self%waverot,self%latt_cur,kpoints,self%ikp_fbz_list)
      select case (self%elph_set%wf_redistribute_opt)
        case (elph_workload_mode_ibz)
          call wf_redistribute_compute_workload_ibz( &
              wf_redis,self%elph_set,elph_accumulators,self%ikp_fbz_list,&
              self%fbz2ibz,self%w_dense_ibz,self%wdes_dense_ibz,&
              kpoints,self%phonon_freqs_ibz,self%latt_cur,self%wave_map,self%comm_kbatch,self%io)
        case (elph_workload_mode_fbz)
          call wf_redistribute_compute_workload( &
              wf_redis,self%elph_set,elph_accumulators,self%ikp_fbz_list,&
              self%fbz2ibz,self%w_dense_ibz,self%wdes_dense_ibz,&
              kpoints,self%phonon_freqs_ibz,self%latt_cur,self%wave_map,self%comm_kbatch,self%io)
        case default
          call vtutor%error('Invalid value of elph_wf_redistribute_opt '//str(self%elph_set%wf_redistribute_opt))
      end select
      call wf_redistribute_write_workload_cpu('current',wf_redis%workload_orig,wf_redis%nkpts_batch_orig,self%io%iu0)
      call wf_redistribute_write_workload_kpoint(wf_redis,self%io%iu0)
      call wf_redistribute_workload(wf_redis,self%wdes_dense_ibz,kpoints,self%comm_kbatch)
      call wf_redistribute_write_workload_cpu('new',wf_redis%workload,wf_redis%nkpts_batch,self%io%iu0)

      ! redistribute the wavefunctions according to the new workload
      ! wave_map, w_kbatch, wdes_kbatch, dproj_kbatch are updated
      select case (self%elph_set%wf_comm_opt)
        case(elph_comm_twosided)
          call wf_redistribute_wfs_twosided(wf_redis,self%nrspinors,self%nproj,self%wavwin_cw,self%wavwin_cproj,&
                                            self%latt_cur,self%latt_ini,self%wdes_dense_ibz,self%wave_map,self%wdes_kbatch,self%w_kbatch,self%dproj_kbatch,self%comm_global)
        case(elph_comm_onesided)
          call wf_redistribute_wfs_window(wf_redis,self%nrspinors,self%nproj,self%wavwin_cw,self%wavwin_cproj,&
                                          self%latt_cur,self%latt_ini,self%wdes_dense_ibz,self%wave_map,self%wdes_kbatch,self%w_kbatch,self%dproj_kbatch)
      end select
      call elph_mels_pw_windows_init(self)

      call wf_redistribute_free(wf_redis)
      PROFILING_STOP('wf_redistribute_driver')
  end subroutine elph_mels_pw_wf_redistribute

  subroutine wf_redistribute_init(self,nkpts,comm_kbatch)
      type(wf_redistribution) :: self
      integer,intent(in) :: nkpts
      type(communic),intent(in) :: comm_kbatch
      allocate(self%workload_orig(comm_kbatch%ncpu))
      self%workload_orig = 0
      allocate(self%nkpts_batch_orig(comm_kbatch%ncpu))
      self%nkpts_batch_orig = 0
      allocate(self%workload(comm_kbatch%ncpu))
      self%workload = 0
      allocate(self%nkpts_batch(comm_kbatch%ncpu))
      self%nkpts_batch = 0
      allocate(self%workload_ikp(nkpts))
      self%workload_ikp = 0
  end subroutine wf_redistribute_init

  !> compute workload given the current WF distribution
  !> this routine sets:
  !>  wf_redistribution%workload_orig - original workload for each mpi rank
  !>  wf_redistribution%nkpts_batch_orig - original number of k-points per batch
  !>  wf_redistribution%workload_ikp  - workload per k' point
  !> it considers the scattering between k in the ibz and k' in the fbz
  subroutine wf_redistribute_compute_workload(self,elph_set,elph_accumulators,ikp_fbz_list,fbz2ibz,&
                                              w_dense_ibz,wdes_dense_ibz,kpoints,phonon_freqs_ibz,latt_cur,wave_map,comm_kbatch,io)
      use mkpoints, only: ikibz_kindex, vkpt_kindex
      use triplets, only: triplet_type, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_vkpt, kpoint_init_kint
      use elphon_base, only: find_value, elph_progress
      use elphon_accumulators_high, only: elph_accumulators_selfen, elph_accumulators_shallido, elph_accumulators_delta_is_zero
      type(wf_redistribution) :: self
      type(elph_settings) :: elph_set
      type(elph_accumulators_selfen) :: elph_accumulators
      integer,target,intent(in) :: ikp_fbz_list(:)
      type(fbz2ibz_stars),intent(in) :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
      type(wavespin),intent(in) :: w_dense_ibz !< wavefunction storage in the dense ibz
      type(wavedes),intent(in) :: wdes_dense_ibz !< wavefunction descriptor for the dense ibz
      type(kpoints_struct),intent(in) :: kpoints
      real(q),intent(in) :: phonon_freqs_ibz(:,:)
      type(latt),intent(in) :: latt_cur
      type(wave_mapper),intent(in) :: wave_map !< map local WFs to points in the ibz
      type(communic),intent(in) :: comm_kbatch
      type(in_struct),intent(in) :: io
      ! local variables
      type(triplet_type) :: triplet
      !integer :: ik, ik_ibz, ikp_ibz, isp, idx, workload
      !integer :: ik_fbz, ikp_fbz, iq_fbz
      !logical :: contributes, delta_is_zero
      integer :: ik, ik_ibz, isp, workload
      integer :: ikp_fbz
      logical :: delta_is_zero

      PROFILING_START('compute_workload')
      !do_iqfbz = 0

      do isp=1,wdes_dense_ibz%ispin
          ! loop over all local fbz k' points
          do ik=1,size(ikp_fbz_list)
              call elph_progress('workload fbz kpoint',ik,size(ikp_fbz_list),io%iu0)
              ikp_fbz = ikp_fbz_list(ik)
              call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)
              ! loop over all the k points
              do ik_ibz=1,wdes_dense_ibz%nkpts
                  if (.not.elph_accumulators_shallido(elph_accumulators,ik_ibz,isp)) cycle

                  call kpoint_init_kibz(triplet%k, kpoints, latt_cur, ik_ibz)
                  call kpoint_init_kint(triplet%q, kpoints, triplet%k%kint-triplet%kp%kint)

                  delta_is_zero = elph_accumulators_delta_is_zero( &
                                    elph_accumulators,elph_set,w_dense_ibz,phonon_freqs_ibz,&
                                    kpoints,triplet,isp)
                  if (delta_is_zero) cycle

                  !do_iqfbz(triplet%q%kindex) = &
                  !do_iqfbz(triplet%q%kindex) + 1

                  workload = 1
                  self%workload_orig(comm_kbatch%node_me) = &
                  self%workload_orig(comm_kbatch%node_me) + workload
                  self%workload_ikp(triplet%kp%ibz) = &
                  self%workload_ikp(triplet%kp%ibz) + workload
              enddo !k
          enddo !k'
      enddo !isp

      !write(*,*) 'will do', count(do_iqfbz/=0), '/', size(do_iqfbz), &
      !           'max:', maxval(do_iqfbz), &
      !           'avg:', 1.0_q*sum(do_iqfbz)/count(do_iqfbz/=0), &
      !           'potential savings:', 100.0_q*count(do_iqfbz/=0)/sum(do_iqfbz)
      self%nkpts_batch_orig(comm_kbatch%node_me) = wave_map%nkpts_batch

      CALLMPI( M_sum_i( comm_kbatch, self%workload_orig, size(self%workload_orig)) )
      CALLMPI( M_sum_i( comm_kbatch, self%nkpts_batch_orig, size(self%nkpts_batch_orig)) )
      CALLMPI( M_sum_i( comm_kbatch, self%workload_ikp, size(self%workload_ikp)) )
      PROFILING_STOP('compute_workload')
  end subroutine wf_redistribute_compute_workload

  !> compute workload given the current WF distribution
  !> this routine sets:
  !>  wf_redistribution%workload_orig - original workload for each mpi rank
  !>  wf_redistribution%nkpts_batch_orig - original number of k-points per batch
  !>  wf_redistribution%workload_ikp  - workload per k' point
  !> it considers the scattering between k in the ibz and k' also in the ibz
  !> this is faster than the above
  subroutine wf_redistribute_compute_workload_ibz(self,elph_set,elph_accumulators,ikp_fbz_list,fbz2ibz,&
                                                   w_dense_ibz,wdes_dense_ibz,kpoints,phonon_freqs_ibz,latt_cur,wave_map,comm_kbatch,io)
      use mkpoints, only: ikibz_kindex, vkpt_kindex
      use triplets, only: triplet_type, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_vkpt, kpoint_init_kint
      use elphon_base, only: find_value, elph_progress
      use elphon_accumulators_high, only: elph_accumulators_delta_is_zero, elph_accumulators_selfen, elph_accumulators_shallido
      type(wf_redistribution) :: self
      type(elph_settings) :: elph_set
      type(elph_accumulators_selfen) :: elph_accumulators
      integer,target,intent(in) :: ikp_fbz_list(:)
      type(fbz2ibz_stars),intent(in) :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
      type(wavespin),intent(in) :: w_dense_ibz !< wavefunction storage in the dense ibz
      type(wavedes),intent(in) :: wdes_dense_ibz !< wavefunction descriptor for the dense ibz
      type(kpoints_struct),intent(in) :: kpoints
      real(q),intent(in) :: phonon_freqs_ibz(:,:)
      type(latt),intent(in) :: latt_cur
      type(wave_mapper),intent(in) :: wave_map !< map local WFs to points in the ibz
      type(communic),intent(in) :: comm_kbatch
      type(in_struct),intent(in) :: io
      ! local variables
      type(triplet_type) :: triplet
      integer :: ikp_ibz, isp, workload
      logical :: delta_is_zero
      integer :: ik, ik_ibz

      PROFILING_START('compute_workload_ibz')

      do isp=1,wdes_dense_ibz%ispin
          ! loop over all the local k' points
          do ik=1,wave_map%nkpts_batch
              call elph_progress('workload ibz kpoint',ik,wave_map%nkpts_batch,io%iu0)
              ikp_ibz = wave_map%kpoints_index(ik)
              call kpoint_init_kibz(triplet%kp, kpoints, latt_cur, ikp_ibz)
              ! loop over all the k points
              do ik_ibz=1,wdes_dense_ibz%nkpts
                  if (.not.elph_accumulators_shallido(elph_accumulators,ik_ibz,isp)) cycle

                  call kpoint_init_kibz(triplet%k, kpoints, latt_cur, ik_ibz)
                  call kpoint_init_kint(triplet%q, kpoints, triplet%k%kint-triplet%kp%kint)

                  delta_is_zero = elph_accumulators_delta_is_zero( &
                                    elph_accumulators,elph_set,w_dense_ibz,phonon_freqs_ibz,&
                                    kpoints,triplet,isp)
                  if (delta_is_zero) cycle

                  !do_iqfbz(triplet%q%kindex) = &
                  !do_iqfbz(triplet%q%kindex) + 1

                  workload = fbz2ibz%star_size(ikp_ibz)

                  self%workload_orig(comm_kbatch%node_me) = &
                  self%workload_orig(comm_kbatch%node_me) + workload
                  self%workload_ikp(triplet%kp%ibz) = &
                  self%workload_ikp(triplet%kp%ibz) + workload
              enddo !k
          enddo !k'
      enddo !isp

      !write(*,*) 'will do', count(do_iqfbz/=0), '/', size(do_iqfbz), &
      !           'max:', maxval(do_iqfbz), &
      !           'avg:', 1.0_q*sum(do_iqfbz)/count(do_iqfbz/=0), &
      !           'potential savings:', 100.0_q*count(do_iqfbz/=0)/sum(do_iqfbz)
      self%nkpts_batch_orig(comm_kbatch%node_me) = wave_map%nkpts_batch

      CALLMPI( M_sum_i( comm_kbatch, self%workload_orig, size(self%workload_orig)) )
      CALLMPI( M_sum_i( comm_kbatch, self%nkpts_batch_orig, size(self%nkpts_batch_orig)) )
      CALLMPI( M_sum_i( comm_kbatch, self%workload_ikp, size(self%workload_ikp)) )
      PROFILING_STOP('compute_workload_ibz')
  end subroutine wf_redistribute_compute_workload_ibz

  !> write workload per cpu
  subroutine wf_redistribute_write_workload_cpu(tag,workload,nkpts_batch,iu)
      use tutor, only: vtutor
      character(len=*) :: tag
      integer,intent(in) :: workload(:)
      integer,intent(in) :: nkpts_batch(:)
      integer :: iu
      ! local variables
      integer :: i

      if (iu<0) return

      if (size(workload) /= size(nkpts_batch)) call vtutor%bug('Incomaptible dimensions of workload and nkpts_batch array',__FILE__,__LINE__)
      write(iu,*) tag, ' workload distribution per cpu'
      write(iu,*) '  cpu  workload     nkpts'
      do i=1,size(workload)
         write(iu,'(i6,2i10)') i, workload(i), nkpts_batch(i)
      enddo
      write(iu,*) 'sum', sum(workload)
      write(iu,*)
  end subroutine wf_redistribute_write_workload_cpu

  !> write workload per kpoint
  subroutine wf_redistribute_write_workload_kpoint(self,iu)
      type(wf_redistribution),intent(in) :: self
      integer :: iu
      ! local variables
      integer :: ikp_ibz
      integer :: workload

      if (iu<0) return
      write(iu,*) 'current workload distribution per kpoint'
      write(iu,*) '  ikp_ibz  workload'
      do ikp_ibz=1,size(self%workload_ikp)
         ! don't print k-points for which there is no workload
         workload = self%workload_ikp(ikp_ibz)
         if (workload==0) cycle
         write(iu,'(i6,i10)') ikp_ibz, workload
      enddo
      write(iu,*) 'sum', sum(self%workload_ikp)
      write(iu,*)
  end subroutine wf_redistribute_write_workload_kpoint

  !> Compute k-point redistribution such that the workload is similar among all the cpus
  !> Will write the wave_map, workload and workload_ikp in self
  subroutine wf_redistribute_workload(self,wdes_dense_ibz,kpoints,comm_kbatch)
      use wave_interpolate, only: wave_mapper_init
      type(wf_redistribution),intent(inout) :: self
      type(wavedes),intent(in) :: wdes_dense_ibz
      type(kpoints_struct),intent(in) :: kpoints
      type(communic),intent(in) :: comm_kbatch
      ! local variables
      integer :: ik, ikp_ibz, nkpts, cpu_next
      integer,allocatable :: workload_ikp_sorted(:)
      integer,allocatable :: sortidx(:)
      integer,allocatable :: kpoints_index(:)
      real(q),allocatable :: kpoints_batch(:,:)

      ! sort k-points in the irreducible brillouin zone by workload
      allocate(workload_ikp_sorted(wdes_dense_ibz%nkpts))
      allocate(sortidx(wdes_dense_ibz%nkpts))
      workload_ikp_sorted = self%workload_ikp(:)
      do ikp_ibz=1,wdes_dense_ibz%nkpts
          sortidx(ikp_ibz) = ikp_ibz
      enddo
      call sort_redis(wdes_dense_ibz%nkpts,workload_ikp_sorted,sortidx)

      nkpts = 0
      allocate(kpoints_index(wdes_dense_ibz%nkpts))
      allocate(kpoints_batch(3,wdes_dense_ibz%nkpts))
      ! Now we know how much work each k' point requires
      ! so we distribute the work equally among CPUS
      do ik=1,wdes_dense_ibz%nkpts
         ! we loop first over the k' points that have the largest workload
         ikp_ibz = sortidx(ik)
         ! if there is no work to be done for this k-point then skip it
         if (self%workload_ikp(ikp_ibz)==0) cycle
         ! get kbatch  with least amount work and increment
         cpu_next = minloc(self%workload,1)
         self%workload(cpu_next) = &
         self%workload(cpu_next) + self%workload_ikp(ikp_ibz)
         ! if the cpu that gets the workload is not me then we are done here
         if (cpu_next/=comm_kbatch%node_me) cycle
         ! add this k-point to batch
         nkpts = nkpts + 1
         kpoints_batch(:,nkpts) = kpoints%vkpt(:,ikp_ibz)
         kpoints_index(nkpts) = ikp_ibz
      enddo
      self%nkpts_batch(comm_kbatch%node_me) = nkpts
      CALLMPI( M_sum_i( comm_kbatch, self%nkpts_batch, size(self%nkpts_batch)) )

      ! initialize wave_map structure with new redistribution
      call wave_mapper_init(self%wave_map,0,nkpts,kpoints_batch,kpoints_index)

  end subroutine wf_redistribute_workload

  !> Redistribute the wfs according to the wave_map structure in self
  subroutine wf_redistribute_wfs_window(self,nrspinors,nproj,wavwin_cw,wavwin_cproj,&
                                        latt_cur,latt_ini,wdes_dense_ibz,wave_map,wdes_kbatch,w_kbatch,dproj_kbatch)
      use elphon_base, only: find_value
      use wave_windower, only: wave_window_get, wave_window_free
      use wave_rotate, only: get_wdes
      use wave, only: allocw, deallocw
      use ini, only: register_allocate, deregister_allocate
      use iso_c_binding, only: c_loc, c_f_pointer
      use tutor, only: vtutor
      type(wf_redistribution),intent(in) :: self
      integer,intent(in) :: nrspinors
      integer,intent(in) :: nproj
      type(wave_window),intent(in) :: wavwin_cw !< wave window for the plane-wave part of the orbitals
      type(wave_window),intent(in) :: wavwin_cproj !< wave window for the projections part of the orbitals
      type(latt),intent(in) :: latt_cur
      type(latt),intent(in) :: latt_ini
      type(wavedes),intent(in) :: wdes_dense_ibz
      type(wave_mapper),intent(inout) :: wave_map !< [in] the old wavefunction map
                                                  !< [out] new wavefunction map from self%wave_map
      type(wavedes),intent(inout) :: wdes_kbatch !< [in] the wavefunction descriptor for the initial distribution of k-points
                                                 !< [out] the wavefunction descriptor for the new wave_map
      type(wavespin),intent(inout) :: w_kbatch !< [in] wfs for initial distribution of k-points
                                               !< [out] wfs for the for the new wave_map
      GDEF,allocatable,target,intent(inout) :: dproj_kbatch(:,:,:,:,:,:) !< [in] dproj for initial distribution of k-points
                                                             !< [out] dproj for the new wave_map
      ! local variables
      integer :: isp, ikp_ibz, ikp_local, ik, i, ib
      integer :: nbkp, nbkp_global
      logical :: has_wfkp_locally
      real(q),allocatable :: wtkpt(:)
      complex(q),pointer :: cproj_ptr(:,:,:,:)
      complex(q),allocatable :: cproj_tmp(:,:,:,:)
      complex(q),allocatable :: cw_tmp(:,:,:,:)

      PROFILING_START('wf_redistribute_wfs_window')
      ! Allocate memory to store PW and PROJ for the WFs in the new distribution
      allocate(cw_tmp(wdes_kbatch%nrplwv,wdes_kbatch%nbands,self%wave_map%nkpts_batch,wdes_kbatch%ispin))
      allocate(cproj_tmp(nproj*4*nrspinors,wdes_kbatch%nbands,self%wave_map%nkpts_batch,wdes_kbatch%ispin))
      ! Check if dimensions of dproj are correct
      if (size(dproj_kbatch,4)/=wdes_kbatch%nbands) call vtutor%bug('wrong 4 dim for dproj_kbatch',__FILE__,__LINE__)
      if (size(dproj_kbatch,5)/=wdes_kbatch%nkpts) call vtutor%bug('wrong 5 dim for dproj_kbatch',__FILE__,__LINE__)
      if (size(dproj_kbatch,6)/=wdes_kbatch%ispin) call vtutor%bug('wrong 6 dim for dproj_kbatch',__FILE__,__LINE__)
      ! Create a pointer to dproj_kbatch
      call c_f_pointer(c_loc(dproj_kbatch),cproj_ptr,&
           [nproj*4*nrspinors,wdes_kbatch%nbands,wdes_kbatch%nkpts,wdes_kbatch%ispin])

      ! report new and old distribution
      !write(*,'(a,i3,a)',advance='no') 'old redis', comm%node_me, ':'
      !do i=1,size(wave_map%kpoints_index)
      !  write(*,'(i3)',advance='no') wave_map%kpoints_index(i)
      !enddo
      !write(*,*)

      !write(*,'(a,i3,a)',advance='no') 'new redis', comm%node_me, ':'
      !do i=1,size(self%wave_map%kpoints_index)
      !  write(*,'(i3)',advance='no') self%wave_map%kpoints_index(i)
      !enddo
      !write(*,*)

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cw    old", comm%node_me, ib
      !do i=1,size(w_kbatch%CW,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(w_kbatch%CW(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cproj old", comm%node_me, ib
      !do i=1,size(cproj_ptr,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cproj_ptr(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      ! Fill in the data locally with the new distribution determined by self%wave_map
      do isp=1,wdes_kbatch%ispin
         do ik=1,self%wave_map%nkpts_batch
            ikp_ibz = self%wave_map%kpoints_index(ik)
            ! careful here! we need to use the old wave_map to get the ikp locally
            ikp_local = find_value(wave_map%kpoints_index,ikp_ibz)
            has_wfkp_locally = ikp_local/=0
            do nbkp_global=1,wdes_kbatch%nb_tot
               if ( mod(nbkp_global-1,wdes_kbatch%nb_par)+1 == wdes_kbatch%nb_low) then
                  nbkp=1+(nbkp_global-1)/wdes_kbatch%nb_par
                  if (has_wfkp_locally) then
                     cw_tmp(:,nbkp,ik,isp)  = w_kbatch%CW(:,nbkp,ikp_local,isp)
                     cproj_tmp(:,nbkp,ik,isp) = cproj_ptr(:,nbkp,ikp_local,isp)
                  else
                     call wave_window_get(wavwin_cw,   nbkp_global,ikp_ibz,isp,   cw_tmp(:,nbkp,ik,isp))
                     call wave_window_get(wavwin_cproj,nbkp_global,ikp_ibz,isp,cproj_tmp(:,nbkp,ik,isp))
                  endif
               endif
            enddo !bands
         enddo !kpoints
      enddo !spin
      call wave_window_free(wavwin_cw)
      call wave_window_free(wavwin_cproj)

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cw    new", comm%node_me, ib
      !do i=1,size(cw_tmp,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cw_tmp(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cproj new", comm%node_me, ib
      !do i=1,size(cw_tmp,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cproj_tmp(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      ! replace wave_map
      wave_map = self%wave_map

      ! deallocate old dproj array and allocate new
      call deregister_allocate(VMEM_ARR_SIZE(dproj_kbatch),"dproj")
      deallocate(dproj_kbatch)
      allocate(dproj_kbatch(nproj,0:3,nrspinors,&
                                  wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin))
      call register_allocate(VMEM_ARR_SIZE(dproj_kbatch),"dproj")
      call c_f_pointer(c_loc(dproj_kbatch),cproj_ptr,&
           [nproj*4*nrspinors,wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin])

      ! deallocate old wf
      call deallocw(w_kbatch)

      ! now create wf descriptor for new set of k-points
      deallocate(wdes_kbatch%vkpt)
      deallocate(wdes_kbatch%wtkpt)
      allocate(wtkpt(wave_map%nkpts_batch))
      wtkpt = 1.0
      call get_wdes(wdes_kbatch,wdes_dense_ibz,wave_map%nkpts_batch,wave_map%kpoints_batch,wtkpt,latt_cur,latt_ini)

      ! allocate new wf
      call allocw(wdes_kbatch,w_kbatch)

      ! copy data
      do isp=1,wdes_kbatch%ispin
         do ik=1,wdes_kbatch%nkpts
            do nbkp=1,wdes_kbatch%nbands
               w_kbatch%CW(:,nbkp,wave_map%nkpts_orig+ik,isp) =  cw_tmp(:,nbkp,ik,isp)
               cproj_ptr(:,nbkp,wave_map%nkpts_orig+ik,isp) = cproj_tmp(:,nbkp,ik,isp)
            enddo
         enddo
      enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cw    new", comm%node_me, ib
      !do i=1,size(w_kbatch%CW,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(w_kbatch%CW(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cproj new", comm%node_me, ib
      !do i=1,size(cproj_ptr,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cproj_ptr(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      deallocate(cw_tmp)
      deallocate(cproj_tmp)
      PROFILING_STOP('wf_redistribute_wfs_window')
  end subroutine wf_redistribute_wfs_window

  !> Redistribute the wfs according to the wave_map structure in self
  subroutine wf_redistribute_wfs_twosided(self,nrspinors,nproj,wavwin_cw,wavwin_cproj,&
                                          latt_cur,latt_ini,wdes_dense_ibz,wave_map,wdes_kbatch,w_kbatch,dproj_kbatch,comm)
      use elphon_base, only: find_value
      use wave_windower, only: wave_window_get, wave_window_free
      use wave_rotate, only: get_wdes
      use wave, only: allocw, deallocw
      use ini, only: register_allocate, deregister_allocate
      use iso_c_binding, only: c_loc, c_f_pointer
      use tutor, only: vtutor
#ifdef MPI
      use mpimy, only: mpi_double_complex, MPI_STATUS_SIZE, MPI_TAG_UB
#endif
     type(wf_redistribution),intent(in) :: self
      integer,intent(in) :: nrspinors
      integer,intent(in) :: nproj
      type(wave_window),intent(in) :: wavwin_cw !< wave window for the plane-wave part of the orbitals
      type(wave_window),intent(in) :: wavwin_cproj !< wave window for the projections part of the orbitals
      type(latt),intent(in) :: latt_cur
      type(latt),intent(in) :: latt_ini
      type(wavedes),intent(in) :: wdes_dense_ibz
      type(wave_mapper),intent(inout) :: wave_map !< [in] the old wavefunction map
                                                  !< [out] new wavefunction map from self%wave_map
      type(wavedes),intent(inout) :: wdes_kbatch !< [in] the wavefunction descriptor for the initial distribution of k-points
                                                 !< [out] the wavefunction descriptor for the new wave_map
      type(wavespin),intent(inout) :: w_kbatch !< [in] wfs for initial distribution of k-points
                                               !< [out] wfs for the for the new wave_map
      GDEF,allocatable,target,intent(inout) :: dproj_kbatch(:,:,:,:,:,:) !< [in] dproj for initial distribution of k-points
                                                             !< [out] dproj for the new wave_map
      type(communic),intent(in) :: comm !< global communicator
      ! local variables
      integer :: isp, ikp_ibz, ikp_local, ik, ib
      integer :: nbkp, nbkp_global, nmaxrequests
      integer :: i, tag, tag_cproj, irank, ireq_recv, ireq_send, icore, ncore, ierr, nkpoints
      integer,allocatable :: bks_rank_dest(:,:,:,:)
      integer,allocatable :: request_send_cw(:)
      integer,allocatable :: request_send_cproj(:)
      integer,allocatable :: request_recv_cw(:)
      integer,allocatable :: request_recv_cproj(:)
      integer,allocatable :: statuses(:,:)
      logical :: has_wfkp_locally
      real(q),allocatable :: wtkpt(:)
      complex(q),pointer :: cproj_ptr(:,:,:,:)
      complex(q),allocatable :: cproj_tmp(:,:,:,:)
      complex(q),allocatable :: cw_tmp(:,:,:,:)

#ifdef MPI
      PROFILING_START('wf_redistribute_wfs_twosided')
      ierr = 0
      ! Allocate memory to store PW and PROJ for the WFs in the new distribution
      allocate(cw_tmp(wdes_kbatch%nrplwv,wdes_kbatch%nbands,self%wave_map%nkpts_batch,wdes_kbatch%ispin))
      allocate(cproj_tmp(nproj*4*nrspinors,wdes_kbatch%nbands,self%wave_map%nkpts_batch,wdes_kbatch%ispin))
      ! Check if dimensions of dproj are correct
      if (size(dproj_kbatch,4)/=wdes_kbatch%nbands) call vtutor%bug('wrong 4 dim for dproj_kbatch',__FILE__,__LINE__)
      if (size(dproj_kbatch,5)/=wdes_kbatch%nkpts) call vtutor%bug('wrong 5 dim for dproj_kbatch',__FILE__,__LINE__)
      if (size(dproj_kbatch,6)/=wdes_kbatch%ispin) call vtutor%bug('wrong 6 dim for dproj_kbatch',__FILE__,__LINE__)
      ! Create a pointer to dproj_kbatch
      call c_f_pointer(c_loc(dproj_kbatch),cproj_ptr,&
           [nproj*4*nrspinors,wdes_kbatch%nbands,wdes_kbatch%nkpts,wdes_kbatch%ispin])

      ! report new and old distribution
      !write(*,'(a,i3,a)',advance='no') 'old redis', comm%node_me, ':'
      !do i=1,size(wave_map%kpoints_index)
      !  write(*,'(i3)',advance='no') wave_map%kpoints_index(i)
      !enddo
      !write(*,*)

      !write(*,'(a,i3,a)',advance='no') 'new redis', comm%node_me, ':'
      !do i=1,size(self%wave_map%kpoints_index)
      !  write(*,'(i3)',advance='no') self%wave_map%kpoints_index(i)
      !enddo
      !write(*,*)

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cw    old", comm%node_me, ib
      !do i=1,size(w_kbatch%CW,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(w_kbatch%CW(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cproj old", comm%node_me, ib
      !do i=1,size(cproj_ptr,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cproj_ptr(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      ! we know the ranks where the WFs are self%bks_rank(self%icore,ib,ikpt,ispin)
      ! so these will be the ranks passed to the receive statements
      ! we need to know where the WFs will go bks_rank_dest()
      ! so these will be the ranks passed to the send statements
      nkpoints = wavwin_cw%nkpoints
      ncore = wdes_kbatch%comm_inb%ncpu
      icore = wdes_kbatch%comm_inb%node_me
      allocate(bks_rank_dest(ncore,wdes_kbatch%nb_tot,nkpoints,wdes_kbatch%ispin))

      bks_rank_dest(:,:,:,:) = 0
      ! fill in destination array
      do isp=1,wdes_kbatch%ispin
         do ik=1,self%wave_map%nkpts_batch
            ikp_ibz = self%wave_map%kpoints_index(ik)
            ! careful here! we need to use the old wave_map to get the ikp locally
            ikp_local = find_value(wave_map%kpoints_index,ikp_ibz)
            has_wfkp_locally = ikp_local/=0
            do nbkp_global=1,wdes_kbatch%nb_tot
               if ( mod(nbkp_global-1,wdes_kbatch%nb_par)+1 == wdes_kbatch%nb_low) then
                  nbkp=1+(nbkp_global-1)/wdes_kbatch%nb_par
                  if (has_wfkp_locally) then
                     cw_tmp(:,nbkp,ik,isp)  = w_kbatch%CW(:,nbkp,ikp_local,isp)
                     cproj_tmp(:,nbkp,ik,isp) = cproj_ptr(:,nbkp,ikp_local,isp)
                  else
                     bks_rank_dest(icore,nbkp_global,ikp_ibz,isp) = comm%node_me
                  endif
               endif
            enddo !bands
         enddo !kpoints
      enddo !spin
      CALLMPI( M_sum_i(comm, bks_rank_dest, size(bks_rank_dest)) )

      ! TODO this can be smaller most of the times but I need to precalculate it
      nmaxrequests = wdes_kbatch%nb_tot*nkpoints*wdes_kbatch%ispin
      allocate(request_send_cw(nmaxrequests))
      allocate(request_send_cproj(nmaxrequests))
      allocate(request_recv_cw(nmaxrequests))
      allocate(request_recv_cproj(nmaxrequests))

      ireq_recv = 0
      ireq_send = 0
      ! send statements
      ! we loop over the current wave_map and send all the data that is required by other nodes
      do isp=1,wdes_kbatch%ispin
         do ik=1,wave_map%nkpts_batch
            ikp_ibz = wave_map%kpoints_index(ik)
            do nbkp_global=1,wdes_kbatch%nb_tot
               tag = (isp-1)*ncore*wdes_kbatch%nb_tot*nkpoints+&
                     (ikp_ibz-1)*ncore*wdes_kbatch%nb_tot+&
                     (nbkp_global-1)*ncore+&
                     icore
               tag_cproj = tag+ncore*wdes_kbatch%nb_tot*nkpoints*wdes_kbatch%ispin
               if ( mod(nbkp_global-1,wdes_kbatch%nb_par)+1 == wdes_kbatch%nb_low) then
                  nbkp=1+(nbkp_global-1)/wdes_kbatch%nb_par
                  irank = bks_rank_dest(icore,nbkp_global,ikp_ibz,isp)
                  if (irank==0) cycle
                  ireq_send = ireq_send+1
                  if (ireq_send>nmaxrequests) then
                     call vtutor%bug('number of requests exceeds the max number of requests',__FILE__,__LINE__)
                  endif
                  CALLMPI( MPI_isend( w_kbatch%CW(:,nbkp,ik,isp), size(w_kbatch%CW,1), MPI_double_complex, irank-1, tag, comm, request_send_cw(ireq_send), ierr) )
                  CALLMPI( MPI_isend( cproj_ptr(:,nbkp,ik,isp),   size(cproj_ptr,1),   MPI_double_complex, irank-1, tag_cproj, comm, request_send_cproj(ireq_send), ierr) )
               endif
            enddo !bands
         enddo !kpoints
      enddo !spin

      ! receive statements
      ! we loop over the new wave_map and receive all the data that we require
      do isp=1,wdes_kbatch%ispin
         do ik=1,self%wave_map%nkpts_batch
            ikp_ibz = self%wave_map%kpoints_index(ik)
            ! careful here! we need to use the old wave_map to get the ikp locally
            ikp_local = find_value(wave_map%kpoints_index,ikp_ibz)
            has_wfkp_locally = ikp_local/=0
            do nbkp_global=1,wdes_kbatch%nb_tot
               tag = (isp-1)*ncore*wdes_kbatch%nb_tot*nkpoints+&
                     (ikp_ibz-1)*ncore*wdes_kbatch%nb_tot+&
                     (nbkp_global-1)*ncore+&
                     icore
               tag_cproj = tag+ncore*wdes_kbatch%nb_tot*nkpoints*wdes_kbatch%ispin
               if ( mod(nbkp_global-1,wdes_kbatch%nb_par)+1 == wdes_kbatch%nb_low) then
                  nbkp=1+(nbkp_global-1)/wdes_kbatch%nb_par
                  if (has_wfkp_locally) then
                     cw_tmp(:,nbkp,ik,isp)  = w_kbatch%CW(:,nbkp,ikp_local,isp)
                     cproj_tmp(:,nbkp,ik,isp) = cproj_ptr(:,nbkp,ikp_local,isp)
                  else
                     !TODO: carefull here wavwin_cw%bks_rank might be different from wavwin_cproj%bks_rank but normally isn't
                     irank = wavwin_cw%bks_rank(icore,nbkp_global,ikp_ibz,isp)
                     ireq_recv = ireq_recv+1
                     if (ireq_recv>nmaxrequests) then
                        call vtutor%bug('number of requests exceeds the max number of requests',__FILE__,__LINE__)
                     endif
                     CALLMPI( MPI_irecv( cw_tmp(:,nbkp,ik,isp),    size(cw_tmp,1),    MPI_double_complex, irank-1, tag, comm, request_recv_cw(ireq_recv), ierr) )
                     CALLMPI( MPI_irecv( cproj_tmp(:,nbkp,ik,isp), size(cproj_tmp,1), MPI_double_complex, irank-1, tag_cproj, comm, request_recv_cproj(ireq_recv), ierr) )
                  endif
               endif
            enddo !bands
         enddo !kpoints
      enddo !spin
      !write(*,*) comm%node_me, icore, "send:", ireq_send, "recv:", ireq_recv

      ! wait for all requests
      allocate(statuses(MPI_STATUS_SIZE,ireq_send))
      CALL MPI_waitall(ireq_send, request_send_cw,    statuses, ierr)
      if (ierr/=0) call vtutor%bug('Bug sending CW',__FILE__,__LINE__)
      CALL MPI_waitall(ireq_send, request_send_cproj, statuses, ierr)
      if (ierr/=0) call vtutor%bug('Bug sending CPROJ',__FILE__,__LINE__)
      deallocate(statuses)

      allocate(statuses(MPI_STATUS_SIZE,ireq_recv))
      CALL MPI_waitall(ireq_recv, request_recv_cw,    statuses, ierr)
      if (ierr/=0) call vtutor%bug('Bug receiving CW',__FILE__,__LINE__)
      CALL MPI_waitall(ireq_recv, request_recv_cproj, statuses, ierr)
      if (ierr/=0) call vtutor%bug('Bug receiving CPROJ',__FILE__,__LINE__)
      deallocate(statuses)

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cw    new", comm%node_me, ib
      !do i=1,size(cw_tmp,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cw_tmp(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cproj new", comm%node_me, ib
      !do i=1,size(cw_tmp,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cproj_tmp(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      ! this is needed otherwise we start replacing the WFs on some nodes before they have been send/arrived on/to others
      call wave_window_free(wavwin_cw)
      call wave_window_free(wavwin_cproj)

      ! replace wave_map
      wave_map = self%wave_map

      ! deallocate old dproj array and allocate new
      call deregister_allocate(VMEM_ARR_SIZE(dproj_kbatch),"dproj")
      deallocate(dproj_kbatch)
      allocate(dproj_kbatch(nproj,0:3,nrspinors,&
                                  wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin))
      call register_allocate(VMEM_ARR_SIZE(dproj_kbatch),"dproj")
      call c_f_pointer(c_loc(dproj_kbatch),cproj_ptr,&
           [nproj*4*nrspinors,wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin])

      ! deallocate old wf
      call deallocw(w_kbatch)

      ! now create wf descriptor for new set of k-points
      deallocate(wdes_kbatch%vkpt)
      deallocate(wdes_kbatch%wtkpt)
      allocate(wtkpt(wave_map%nkpts_batch))
      wtkpt = 1.0
      call get_wdes(wdes_kbatch,wdes_dense_ibz,wave_map%nkpts_batch,wave_map%kpoints_batch,wtkpt,latt_cur,latt_ini)

      ! allocate new wf
      call allocw(wdes_kbatch,w_kbatch)

      ! copy data
      do isp=1,wdes_kbatch%ispin
         do ik=1,wdes_kbatch%nkpts
            do nbkp=1,wdes_kbatch%nbands
               w_kbatch%CW(:,nbkp,wave_map%nkpts_orig+ik,isp) =  cw_tmp(:,nbkp,ik,isp)
               cproj_ptr(:,nbkp,wave_map%nkpts_orig+ik,isp) = cproj_tmp(:,nbkp,ik,isp)
            enddo
         enddo
      enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cw    new", comm%node_me, ib
      !do i=1,size(w_kbatch%CW,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(w_kbatch%CW(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      !do ib=1,wdes_kbatch%nbands
      !write(*,'(a,2i3)',advance='no') "cproj new", comm%node_me, ib
      !do i=1,size(cproj_ptr,3)
      !  write(*,'(F10.4)',advance='no') sum(abs(cproj_ptr(:,ib,i,1)))
      !enddo
      !write(*,*)
      !enddo

      deallocate(cw_tmp)
      deallocate(cproj_tmp)
      PROFILING_STOP('wf_redistribute_wfs_twosided')
#endif
 end subroutine wf_redistribute_wfs_twosided

  subroutine wf_redistribute_free(self)
      type(wf_redistribution) :: self
      deallocate(self%workload_orig)
      deallocate(self%nkpts_batch_orig)
      deallocate(self%workload)
      deallocate(self%nkpts_batch)
      deallocate(self%workload_ikp)
  end subroutine wf_redistribute_free

  subroutine elph_mels_pw_cache_prefill(self,elph_accumulators)
      use elphon_base, only: elph_driver_ph
      use elphon_accumulators_high, only: elph_accumulators_selfen
      type(elph_mels_pw_type) :: self
      type(elph_accumulators_selfen),intent(in) :: elph_accumulators
      select case (self%elph_set%wf_comm_opt)
        case (0)
            call elph_mels_pw_cache_prefill_bcast(self,elph_accumulators)
        case (1)
            call elph_mels_pw_cache_prefill_window(self,elph_accumulators)
      end select
  end subroutine elph_mels_pw_cache_prefill

  !> Fill in the WF cache with the states that should be present on all the nodes
  subroutine elph_mels_pw_cache_prefill_window(self,elph_accumulators)
      use elphon_accumulators_high, only: elph_accumulators_selfen, elph_accumulators_shallido
      use wave_windower, only: wave_window_get, wave_cache_set, wave_window_flush_local_all
      use elphon_base, only: find_value
      type(elph_mels_pw_type) :: self
      type(elph_accumulators_selfen),intent(in) :: elph_accumulators
      ! local variables
      integer :: nbk_global,isp,ikp_fbz,ik_ibz,idx
      logical :: has_wfk_locally
      complex(q),allocatable :: cw1(:)
      complex(q),allocatable :: cproj1(:)

      if (.not.self%elph_set%wf_cache_prefill) return

      PROFILING_START('wf_cache_prefill')
      allocate(cw1(self%wdes_dense_ibz%nrplwv))
      allocate(cproj1(self%nproj*4*self%nrspinors))
      ! Gather all the required WFs
      do isp=1,self%wdes_dense_ibz%ispin
        ikp_fbz=1
        do ik_ibz=1,self%wdes_dense_ibz%nkpts
          if (.not.elph_accumulators_shallido(elph_accumulators,ik_ibz,isp)) cycle
          idx = find_value(self%wave_map%kpoints_index,ik_ibz)
          has_wfk_locally = idx/=0
          ! if I have these orbitals locally then skip
          ! otherwise get orbitals for this k-point
          do nbk_global=self%band_start_k,self%band_stop_k
            ! check if I have the bands locally otherwise get them from the node that has them
            if ( mod(nbk_global-1,self%wdes_dense_ibz%nb_par)+1 == self%wdes_dense_ibz%nb_low .and.&
                 has_wfk_locally) cycle
            call wave_window_get(self%wavwin_cw,nbk_global,ik_ibz,isp,cw1)
            call wave_cache_set(self%cache_cw,nbk_global,ik_ibz,isp,cw1)
            call wave_window_get(self%wavwin_cproj,nbk_global,ik_ibz,isp,cproj1)
            call wave_cache_set(self%cache_cproj,nbk_global,ik_ibz,isp,cproj1)
          enddo
        enddo
      enddo
      ! TODO: this is needed otherwise the communication does not happen
      call elph_mels_pw_windows_free(self)
      call elph_mels_pw_windows_init(self)
      !call wave_window_flush_local_all(self%wavwin_cw)
      !call wave_window_flush_local_all(self%wavwin_cproj)
      PROFILING_STOP('wf_cache_prefill')

  end subroutine elph_mels_pw_cache_prefill_window


  !> Fill in the WF cache with the states that should be present on all the nodes
  subroutine elph_mels_pw_cache_prefill_bcast(self,elph_accumulators)
      use elphon_accumulators_high, only: elph_accumulators_selfen, elph_accumulators_shallido
      use wave_windower, only: wave_window_get, wave_cache_set, wave_window_flush_local_all
      use elphon_base, only: find_value
#ifdef MPI
      use mpimy, only: MPI_double_complex
#endif
      use wave, only: nb_local_low
      use iso_c_binding, only: c_loc, c_f_pointer
      type(elph_mels_pw_type),target :: self
      type(elph_accumulators_selfen),intent(in) :: elph_accumulators
      ! local variables
      integer :: nbk_global,nbk_local,isp,ikp_fbz,ik_ibz,idx,icore,ncore,irank
      integer :: ierr,maxnrplwv,nrplwv,maxnproj,nprojdim
      logical :: has_wfk_locally
      complex(q),allocatable :: cw1(:)
      complex(q),allocatable :: cproj1(:)
      complex(q),pointer :: cproj_ptr(:,:,:,:)


      if (.not.self%elph_set%wf_cache_prefill) return

#ifdef MPI
      PROFILING_START('wf_cache_prefill_bcast')
      maxnrplwv = maxval(self%wavwin_cw%nplw_all)
      nrplwv = self%wdes_dense_ibz%nrplwv
      allocate(cw1(maxnrplwv))

      nprojdim = self%nproj*4*self%nrspinors
      maxnproj = nprojdim
      call M_max_i(self%comm_global,maxnproj,1)
      allocate(cproj1(maxnproj))
      ncore = self%wdes_dense_ibz%comm_inb%ncpu

      call c_f_pointer(c_loc(self%dproj_kbatch),cproj_ptr,&
           [self%nproj*4*self%nrspinors,self%wdes_kbatch%nbands,self%wave_map%nkpts_batch,self%wdes_kbatch%ispin])
      ! Gather all the required WFs
      do isp=1,self%wdes_dense_ibz%ispin
        ikp_fbz=1
        do ik_ibz=1,self%wdes_dense_ibz%nkpts
          if (.not.elph_accumulators_shallido(elph_accumulators,ik_ibz,isp)) cycle
          idx = find_value(self%wave_map%kpoints_index,ik_ibz)
          has_wfk_locally = idx/=0
          ! if I have these orbitals locally then skip
          ! otherwise get orbitals for this k-point
          do nbk_global=self%band_start_k,self%band_stop_k
            nbk_local = nb_local_low(nbk_global,self%wdes_dense_ibz%nb_par,self%wdes_dense_ibz%nb_low)

            ! we want to broadcast information within all the nodes
            do icore=1,ncore
                ! check if I have the bands locally and set them
                if ( mod(nbk_global-1,self%wdes_dense_ibz%nb_par)+1 == self%wdes_dense_ibz%nb_low .and.&
                     has_wfk_locally) then
                    cw1(:nrplwv) = self%w_kbatch%CW(:,nbk_local,self%wave_map%nkpts_orig+idx,isp)
                    cproj1(:nprojdim) = cproj_ptr(:,nbk_local,self%wave_map%nkpts_orig+idx,isp)
                endif

                ! now broadcast this the orbital to all the nodes
                irank = self%wavwin_cw%bks_rank(icore,nbk_global,ik_ibz,isp)
                call MPI_Bcast(cw1,    maxnrplwv, MPI_double_complex, irank-1, self%comm_global, ierr)
                irank = self%wavwin_cproj%bks_rank(icore,nbk_global,ik_ibz,isp)
                call MPI_Bcast(cproj1, maxnproj,  MPI_double_complex, irank-1, self%comm_global, ierr)

                ! whoever needs these orbitals (same icore) should keep them
                if (icore/=self%wdes_dense_ibz%comm_inb%node_me) cycle
                call wave_cache_set(self%cache_cw,nbk_global,ik_ibz,isp,cw1(:nrplwv))
                call wave_cache_set(self%cache_cproj,nbk_global,ik_ibz,isp,cproj1(:nprojdim))
            enddo
          enddo
        enddo
      enddo
      PROFILING_STOP('wf_cache_prefill_bcast')
#endif
  end subroutine elph_mels_pw_cache_prefill_bcast

  subroutine elph_mels_pw_windows_free(self)
      use wave_windower, only: wave_window_free, wave_cache_free
      type(elph_mels_pw_type) :: self
      call wave_window_free(self%wavwin_cw)
      call wave_window_free(self%wavwin_cproj)
  end subroutine elph_mels_pw_windows_free

  subroutine elph_mels_pw_wave_caches_free(self)
      use wave_windower, only: wave_window_free, wave_cache_free
      type(elph_mels_pw_type) :: self
      call wave_cache_free(self%cache_cw)
      call wave_cache_free(self%cache_cproj)
  end subroutine elph_mels_pw_wave_caches_free

  !> Compute dense wavefunction descriptor and allocate without plane-waves (only eigenvalues will be used)
  !> build wdes_dense_ibz and w_dense_ibz inside elph_mels_pw
  !> The waverot and fbz2ibz structure are also initialized here.
  subroutine elph_mels_pw_get_w_dense(self,wdes_dense_ibz,w_dense_ibz,kpoints)
      use tutor, only: vtutor
      use string, only: str
      use wave_rotate, only: wave_rotator_init, fbz2ibz_stars_init, fbz2ibz_stars_dummy_init, kpoint_star, kstar_from_fbz2ibz
      use mkpoints, only: kpoints_is_regular, LineMode, vkpt_kindex
      use elphon_common, only: elphon_get_waverot_fbz2ibz, elphon_get_wdes_w
      type(elph_mels_pw_type),intent(inout) :: self
      type(wavedes),intent(inout) :: wdes_dense_ibz
      type(wavespin),intent(inout) :: w_dense_ibz
      type(kpoints_struct),intent(inout) :: kpoints
      !local variables
      integer :: fft_mesh_wf(3)

      ! First check the consistency of the FFT mesh
      fft_mesh_wf = [self%grid_orig%ngx,self%grid_orig%ngy,self%grid_orig%ngz]
      if (any(fft_mesh_wf/=self%fft_mesh)) then
         call vtutor%error('error: fft grid from phelel '//str(self%fft_mesh)//&
                           ' different from wf '//str(fft_mesh_wf))
      endif

      ! intializes the waverot and fbz2ibz structures also update fbz2ibz mapping in kpoints
      call elphon_get_waverot_fbz2ibz(self%waverot,self%fbz2ibz,self%symm,self%latt_cur,self%latt_ini,&
                                      self%grid_orig,self%ngdim_orig,wdes_dense_ibz,kpoints)

      ! Initialize wavefunction descriptor for the dense IBZ
      ! BEWARE: this changes grid!
      call elphon_get_wdes_w(wdes_dense_ibz,w_dense_ibz,kpoints,self%fbz2ibz,self%grid,self%latt_cur,self%latt_ini,self%ngdim_orig)

  end subroutine elph_mels_pw_get_w_dense

  !> prepare WFs on the regular grid
  subroutine elph_mels_pw_compute_wfs(self,kpoints_orig,kpoints_inter,kgrid,w_orig,hamiltonian,nonlr_s,nonl_s,&
                                   cqij,cdij,sv,chempot,bandgap)
      use hamil_struct_def, only: ham_handle
      use elphon_kgrid, only: elph_kgrid
      use nonl_struct_def, only: nonl_struct
      use nonlr_struct_def, only: nonlr_struct
      use elphon_common, only: elphon_prepare_wave_interp, elphon_run_precalculation, elphon_run_chempot, elphon_compute_projections
      use bandgap_struct, only: bandgap_info
      use electron_transport, only: chemical_potential
      type(elph_mels_pw_type) :: self
      type(kpoints_struct) :: kpoints_orig
      type(kpoints_struct) :: kpoints_inter
      type(elph_kgrid) :: kgrid
      type(wavespin) :: w_orig
      type(ham_handle) :: hamiltonian
      type(nonlr_struct) :: nonlr_s
      type(nonl_struct)  :: nonl_s
      OVERLAP :: cqij(:,:,:,:)
      OVERLAP :: cdij(:,:,:,:)
      RGRID   :: sv(:,:)
      type(chemical_potential) :: chempot
      type(bandgap_info),allocatable :: bandgap(:)
      !
      ! Set up WAVE_INTERP (orbital interpolation).
      !
      call elphon_prepare_wave_interp(self%elph_set%wave_interp_set, self%wave_interp, self%fbz2ibz, kpoints_inter, self%comm_kbatch, &
                           kpoints_orig, self%grid, self%latt_cur, self%latt_ini, &
                           self%t_info, nonl_s, w_orig, self%p, self%symm, self%io, &
                           self%wdes_dense_ibz%ngdim, self%wdes_dense_ibz%nb_tot)

      !
      ! Run pre-calculation of physical properties.
      !
      ! orbital interpolation: w_dense_ibz, wave_interp
      call elphon_run_precalculation(self%w_dense_ibz, self%wave_interp, self%wave_map, &
                                     hamiltonian, self%latt_cur, self%latt_ini, self%t_info, &
                                     nonlr_s, nonl_s, self%w_kbatch, self%lmdim, self%p, sv, cqij, cdij, self%symm, &
                                     self%info, self%elph_set, self%wdes_dense_ibz, self%comm_kbatch, self%comm_global, &
                                     self%wdes_kbatch, kpoints_inter, self%io)
      ! chemical potential and number of electrons: chempot
      call elphon_run_chempot(self%w_dense_ibz, chempot, bandgap, &
                              self%latt_cur, self%latt_ini, self%t_info, &
                              self%symm, &
                              self%info, self%elph_set, self%wdes_dense_ibz, self%comm_kbatch, self%comm_global, &
                              kpoints_inter, self%io)
      ! displacement derivative of projectors calculation: dproj_kbatch
      call elphon_compute_projections(self%wdes_kbatch,self%w_kbatch,nonl_s,self%wave_map,self%dproj_kbatch,self%latt_cur)
  end subroutine elph_mels_pw_compute_wfs

  !> Read electron-phonon potential and initialize structure
  subroutine elph_mels_pw_pot_init(self,comm_shmem)
      use elphon_potential, only: elphon_pot_init
      type(elph_mels_pw_type) :: self
      type(communic),pointer :: comm_shmem

      !> set "dimension" in ELPH_POT (see elphon_pot type) from self, and allocate arrays in ELPH_POT
      call elphon_pot_init(self%elph_pot,self%grid,self%ngridpoints,self%npr_local,self%natoms_pc,self%natoms_pc_local,self%natoms_sc,&
                           self%pot_ncdij,self%det_p2s,self%lmdim,self%fft_mesh,self%elph_set%useshmem,self%elph_set%writepot,self%elph_set%metagga,comm_shmem)
  end subroutine elph_mels_pw_pot_init

  subroutine elph_mels_pw_pot_read_hdf5(self,filename)
      use elphon_potential, only: elphon_pot_from_hdf5
      type(elph_mels_pw_type) :: self
      character(len=*),optional,intent(in) :: filename
      ! read potentials from 'phelel_params.hdf5',
      call elphon_pot_from_hdf5(self%elph_pot,self%wdes_dense_ibz%comm_inb,self%comm_kbatch,filename)
  end subroutine elph_mels_pw_pot_read_hdf5

  subroutine elph_mels_pw_pot_setup(self)
      use elphon_common, only: elphon_prepare_elph_pot
      use elphon_potential, only: elphon_pot_transfer_ncdij
      type(elph_mels_pw_type) :: self
      call elphon_prepare_elph_pot(self%elph_pot, self%elph_set, self%latt_cur, self%grid, self%p, self%t_info, self%io)
      call elphon_pot_transfer_ncdij(self%elph_pot,self%wf_ncdij,self%elph_caches%ccdij,self%elph_caches%ccqij,self%elph_caches%ccdrij)
  end subroutine elph_mels_pw_pot_setup

  !> compute dipersion of the electron-phonon potential along a path
  subroutine elph_mels_pw_pot_dispersion(self,qpoints_exists,filename)
      USE elphon_potential, ONLY: elphon_test_lr, elphon_check_lr_addremove, elphon_write_lr, elphon_check_lr_slowfastft
      USE mkpoints, ONLY: rd_kpoints
      type(elph_mels_pw_type) :: self
      logical :: qpoints_exists
      character(len=*),intent(in) :: filename
      ! local variables
      type(kpoints_struct) :: qpoints
      inquire(file=filename, exist=qpoints_exists)
      if (qpoints_exists) then
          qpoints%ltet = .false.
          call rd_kpoints(qpoints, self%latt_cur, .false., .false., self%io%iu6, self%io%iu0, file_name = filename)

          ! be careful, these subroutines need to be corrected since I changed
          ! self%fft_mesh_scp to self%fft_mesh_sc
          select case (qpoints%sznamk)
            case ('slowfast')
                call elphon_check_lr_slowfastft(self%elph_pot,self%io%iu0)
            case ('writelr')
                call elphon_write_lr(self%elph_pot,self%io%iu0)
            case ('addremove')
                call elphon_check_lr_addremove(self%elph_pot,qpoints%nkpts,qpoints%vkpt)
            case default
                call elphon_test_lr(self%elph_pot,qpoints%nkpts,qpoints%vkpt)
          end select
      endif
  end subroutine elph_mels_pw_pot_dispersion

  !> Generate an integer list with all the k' points in the FBZ
  !> that need to be handled locally
  !> This is done by looping over the k' points in the IBZ
  !> getting the star of each k' point and then
  !> getting the k-point index in the IBZ for each point in the star
  subroutine wave_map_get_ikpfbz(wave_map,waverot,latt_cur,kpoints,ikp_fbzs)
      use wave_rotate, ONLY: wave_rotator_get_star, kpoint_star
      use mkpoints_struct_def, ONLY: LineMode
      use mkpoints, ONLY: vkpt_kindex
      use tutor, only: vtutor
      type(wave_mapper),intent(in) :: wave_map !< map k-points in w_kbatch to the dense IBZ
      type(wave_rotator),intent(in) :: waverot !< wave function rotator
      type(latt),intent(in) :: latt_cur !< current lattice
      type(kpoints_struct),intent(in) :: kpoints
      integer,allocatable,intent(out) :: ikp_fbzs(:)
      ! local variables
      type(kpoint_star) :: star
      integer,allocatable :: ikp_fbzs_tmp(:)
      integer :: ik, ik_ibz, nkpts, ik_count, istar

      PROFILING_START('elph_get_ikpfbz')
      ik_count = 0
      nkpts = wave_map%nkpts_batch-wave_map%nkpts_orig
      allocate(ikp_fbzs_tmp(nkpts*48))
      ! loop over k' points in the local ibz
      do ik=1,wave_map%nkpts_batch
          ik_ibz = wave_map%kpoints_index(ik+wave_map%nkpts_orig)
          ! loop over points in the star of k'
          call wave_rotator_get_star(waverot,kpoints%vkpt(:,ik_ibz),star)
          do istar=1,star%nkpts
              ik_count = ik_count+1
              if (ik_count>size(ikp_fbzs_tmp)) call vtutor%bug('Size of temporary array for k-point integers is insufficient',__FILE__,__LINE__)
              if (kpoints%mode==LineMode) then
                  ikp_fbzs_tmp(ik_count)=ik_count
              else
                  ikp_fbzs_tmp(ik_count)=vkpt_kindex(kpoints,latt_cur,star%vkpt(:,istar))
              endif
          enddo
      enddo
      ! allocate and copy filled data
      allocate(ikp_fbzs(ik_count))
      ikp_fbzs(:) = ikp_fbzs_tmp(:ik_count)
      deallocate(ikp_fbzs_tmp)
      PROFILING_STOP('elph_get_ikpfbz')
  end subroutine wave_map_get_ikpfbz

  !> Subroutine to get the WFs and projections for a particular k-point
  subroutine get_pws_ibz(full,band_start,band_stop,ikibz,isp,w_dense_ibz,w_kbatch,cache_cw,wavwin_cw,wave_map,wdesk,wk)
      use tutor, only: vtutor
      use string, only: str
      use elphon_base, only : find_value
      use wave_high, only : elements, setwdes, newwava, wa_copy_cw
      use wave_windower, only : wave_window_cache_get, wave_window_cache_getk
      logical,intent(in) :: full
      integer,intent(in) :: band_start !< starting band for which to get the plane-waves
      integer,intent(in) :: band_stop !< last band for which to get the plane-waves
      integer,intent(in) :: ikibz !< coordinates of the points in the ibz
      integer,intent(in) :: isp !< spin index
      type(wavespin),intent(in) :: w_dense_ibz
      type(wavespin),intent(in) :: w_kbatch !< local batch of wavefunctions
      type(wave_cache),intent(in) :: cache_cw
      type(wave_window),intent(in) :: wavwin_cw
      type(wave_mapper),intent(in) :: wave_map !< map k-points in w_kbatch to the dense IBZ
      type(wavedes1),intent(inout) :: wdesk
      type(wavefuna) :: wk
      ! local variables
      logical :: has_wfk_locally
      integer :: idx, iband, nbands, nbk, nbk_global

      ! check if we have this k-point locally
      idx = find_value(wave_map%kpoints_index,ikibz)
      has_wfk_locally = idx/=0

      if (full) then
          nbands = band_stop-band_start+1
          wdesk%nbands = nbands
          call newwava(wk,wdesk,nbands)
          !in this case we get the pw coefficients for all the bands between band_start and band_stop
          do nbk_global=band_start,band_stop
              iband = nbk_global-band_start+1
              ! Check if I have the bands locally otherwise get them from the node that has them
              if ( mod(nbk_global-1,wdesk%nb_par)+1 == wdesk%nb_low .and.&
                   has_wfk_locally) then
                   ! I have the WF locally at this index
                   nbk=1+(nbk_global-1)/wdesk%nb_par
                   wk%CW(:,iband)=w_kbatch%CW(:,nbk,wave_map%nkpts_orig+idx,isp)
              else
                   call wave_window_cache_get(wavwin_cw,cache_cw,nbk_global,ikibz,isp,wk%CW(:,iband))
              endif
              wk%celen(iband) = w_dense_ibz%celtot(nbk_global,ikibz,isp)
          enddo
      else
          wdesk%nbands=0
          do nbk_global=band_start,band_stop
              if ( mod(nbk_global-1,wdesk%nb_par)+1 /= wdesk%nb_low ) cycle
              wdesk%nbands = wdesk%nbands+1
          enddo
          ! allocate a container for the wavefunction
          call newwava(wk,wdesk,wdesk%nbands)
          ! request plane-wave part of the wavefunction
          ! Check if I have the bands locally otherwise get them from the node that has them
          do nbk_global=band_start,band_stop
              if ( mod(nbk_global-1,wdesk%nb_par)+1 == wdesk%nb_low .and.&
                   has_wfk_locally) then
                   ! I have the WF locally at this index
                   nbk=1+(nbk_global-1)/wdesk%nb_par
                   wk%CW(:,nbk-band_start+1)=w_kbatch%CW(:,nbk,wave_map%nkpts_orig+idx,isp)
              else
                   call wave_window_cache_get(wavwin_cw,cache_cw,nbk_global,ikibz,isp,wk%CW(:,nbk))
              endif
          enddo
          do nbk_global=band_start,band_stop
              if ( mod(nbk_global-1,wdesk%nb_par)+1 == wdesk%nb_low ) then
                nbk=1+(nbk_global-1)/wdesk%nb_par
                wk%celen(nbk-band_start+1) = w_dense_ibz%celtot(nbk_global,ikibz,isp)
              endif
          enddo
      endif
  end subroutine get_pws_ibz

  ! get projections from local allocation if available otherwise recompute from wf
  subroutine get_projections_ibz(full,band_start,band_stop,ikibz,isp,dproj_kbatch,wdes_kbatch,cache_cproj,wavwin_cproj,wave_map,dprojk)
      use iso_c_binding, only: c_loc, c_f_pointer
      use tutor, only: vtutor
      use string, only: str
      use elphon_base, only: find_value
      use wave_windower, only : wave_window_cache_get, wave_window_cache_getk
      logical,intent(in) :: full
      integer,intent(in) :: band_start, band_stop
      integer,intent(in) :: ikibz
      integer,intent(in) :: isp
      GDEF,intent(in) :: dproj_kbatch(:,0:,:,:,:,:)
      type(wavedes),intent(in) :: wdes_kbatch !< descriptor of the local batch of wavefunctions
      type(wave_cache),intent(in) :: cache_cproj
      type(wave_window),intent(in) :: wavwin_cproj
      type(wave_mapper),intent(in) :: wave_map !< map k-points in w_kbatch to the dense IBZ
      GDEF :: dprojk(:,:,0:,:) !< (nproj,nbands,0:3,nrspinors)
      !local variables
      integer :: ispinor, ik
      integer :: nbk, nbk_global
      integer :: nproj, nrspinors, ispin
      logical :: has_wfk_locally
      GDEF,allocatable,target :: dproj(:,:,:,:)
      complex(q),pointer :: cproj1_ptr(:)
      complex(q),pointer :: cproj1k_ptr(:,:)

      nproj = size(dproj_kbatch,1)
      nrspinors = size(dproj_kbatch,3)
      ispin = size(dproj_kbatch,6)
      if (size(dprojk,1) /= nproj) then
         call vtutor%bug('Wrong value for dimension 1 in dprojk expected: '//&
                         str(nproj)//' got: '//str(size(dprojk,1)),__FILE__,__LINE__)
      endif
      if (size(dprojk,2) /= band_stop-band_start+1) then
         call vtutor%bug('Wrong value for dimension 2 in dprojk expected: '//&
                         str(band_stop-band_start+1)//' got: '//str(size(dprojk,2)),__FILE__,__LINE__)
      endif
      if (size(dprojk,4) /= nrspinors) then
         call vtutor%bug('Wrong value for dimension 4 in dprojk expected: '//&
                         str(nrspinors)//' got: '//str(size(dprojk,4)),__FILE__,__LINE__)
      endif

      ik = find_value(wave_map%kpoints_index,ikibz)
      has_wfk_locally = ik/=0

      if (full) then
          allocate(dproj(nproj,0:3,nrspinors,1))
          call c_f_pointer(c_loc(dproj),cproj1_ptr,&
                           [nproj*4*nrspinors])
          do nbk_global=band_start,band_stop
            ! check if I have the bands locally otherwise get them from the node that has them
            if ( mod(nbk_global-1,wdes_kbatch%nb_par)+1 == wdes_kbatch%nb_low .and.&
                 has_wfk_locally) then
              ! compute local band index
              nbk=1+(nbk_global-1)/wdes_kbatch%nb_par
              ! I have the dproj locally at this index
              do ispinor=1,nrspinors
                  dprojk(:,nbk_global-band_start+1,0,ispinor) = dproj_kbatch(:,0,ispinor,nbk,ik,isp)
                  dprojk(:,nbk_global-band_start+1,1,ispinor) = dproj_kbatch(:,1,ispinor,nbk,ik,isp)
                  dprojk(:,nbk_global-band_start+1,2,ispinor) = dproj_kbatch(:,2,ispinor,nbk,ik,isp)
                  dprojk(:,nbk_global-band_start+1,3,ispinor) = dproj_kbatch(:,3,ispinor,nbk,ik,isp)
              enddo ! spinor
            else
              ! I have to look for the cproj in the cache
              call wave_window_cache_get(wavwin_cproj,cache_cproj,nbk_global,ikibz,isp,cproj1_ptr)
              do ispinor=1,nrspinors
                  dprojk(:,nbk_global-band_start+1,0,ispinor) = dproj(:,0,ispinor,1)
                  dprojk(:,nbk_global-band_start+1,1,ispinor) = dproj(:,1,ispinor,1)
                  dprojk(:,nbk_global-band_start+1,2,ispinor) = dproj(:,2,ispinor,1)
                  dprojk(:,nbk_global-band_start+1,3,ispinor) = dproj(:,3,ispinor,1)
              enddo ! spinor
            endif
          enddo ! loop over bands
          deallocate(dproj)
      else
          if (has_wfk_locally) then
              do nbk_global=band_start,band_stop
                if (mod(nbk_global-1,wdes_kbatch%nb_par)+1 /= wdes_kbatch%nb_low) cycle
                ! compute local band index
                nbk=1+(nbk_global-1)/wdes_kbatch%nb_par
                do ispinor=1,nrspinors
                  dprojk(:,nbk-band_start+1,0,ispinor)=dproj_kbatch(:,0,ispinor,nbk,ik,isp)
                  dprojk(:,nbk-band_start+1,1,ispinor)=dproj_kbatch(:,1,ispinor,nbk,ik,isp)
                  dprojk(:,nbk-band_start+1,2,ispinor)=dproj_kbatch(:,2,ispinor,nbk,ik,isp)
                  dprojk(:,nbk-band_start+1,3,ispinor)=dproj_kbatch(:,3,ispinor,nbk,ik,isp)
                enddo
              enddo
          else
              call vtutor%bug('Requested wf at '//str(ikibz)//' is not present locally',__FILE__,__LINE__)
              allocate(dproj(nproj,0:3,nrspinors,1))
              call c_f_pointer(c_loc(dproj),cproj1k_ptr,&
                               [nproj*4*nrspinors,wdes_kbatch%nbands])
              do nbk_global=band_start,band_stop
                  if (mod(nbk_global-1,wdes_kbatch%nb_par)+1 /= wdes_kbatch%nb_low) cycle
                  ! compute local band index
                  nbk=1+(nbk_global-1)/wdes_kbatch%nb_par
                  ! I have to look for the cproj in the cache
                  call wave_window_cache_get(wavwin_cproj,cache_cproj,nbk_global,ikibz,isp,cproj1_ptr)
                  do ispinor=1,nrspinors
                      dprojk(:,nbk-band_start+1,0,ispinor) = dproj(:,0,ispinor,1)
                      dprojk(:,nbk-band_start+1,1,ispinor) = dproj(:,1,ispinor,1)
                      dprojk(:,nbk-band_start+1,2,ispinor) = dproj(:,2,ispinor,1)
                      dprojk(:,nbk-band_start+1,3,ispinor) = dproj(:,3,ispinor,1)
                  enddo ! spinor
              enddo ! loop over bands
              deallocate(dproj)
          endif
      endif !if full
  end subroutine get_projections_ibz

  !> get pw and projections at a specific k-point between band_start and band_stop
  subroutine get_wfs_ibz(full,kpoints,band_start,band_stop,ikibz,isp,w_dense_ibz,w_kbatch,wdes_kbatch,dproj_kbatch,wave_map,&
                     wavwin_cw,cache_cw,wavwin_cproj,cache_cproj,latt_cur,latt_ini,wdes1_container,wkfbz,dproj_fbz)
      use wave_high, only: newwava, setwdes, element, delwava
      use wave_rotate, only: wave_rotator, wave_rotator_rotate_pw, wave_rotator_rotate_dproj, wavedes1_container, wdes1_container_init, rotate_vkpt_irot
      use mkpoints, only: kpoints_fbz_to_ibz, kindex_vkpt, ikibz_kindex
      use nonl_struct_def, only: nonl_struct
      use nonl_high, only: spher, phase, nonl_alloc, nonl_dealloc
      use wap, only: get_dproj
      logical,intent(in) :: full
      type(kpoints_struct),intent(in) :: kpoints
      integer :: band_start
      integer :: band_stop
      integer :: isp !< index of the spin channel
      integer :: ikibz !< index of the k-point we want to obtain
      GDEF,intent(in) :: dproj_kbatch(:,0:,:,:,:,:)
      type(wavespin),intent(in) :: w_dense_ibz !< wavefunctions in the IBZ
      type(wavespin),intent(in) :: w_kbatch !< local batch of wavefunctions
      type(wavedes),intent(in) :: wdes_kbatch !< descriptor of the local batch of wavefunctions
      type(wave_mapper) :: wave_map !< mapping of the k-points in the k-point batch to the IBZ
      type(wave_window),intent(in) :: wavwin_cw
      type(wave_cache),intent(in) :: cache_cw
      type(wave_window),intent(in) :: wavwin_cproj
      type(wave_cache),intent(in) :: cache_cproj
      type(latt) :: latt_cur
      type(latt) :: latt_ini
      type(wavedes1_container) :: wdes1_container
      type(wavefuna) :: wkfbz
      GDEF :: dproj_fbz(:,:,0:,:)
      ! local variables
      real(q) :: vk(3)

      PROFILING_START('get_wfs_ibz')
      vk = kpoints%vkpt(:,ikibz)
      call wdes1_container_init(wdes1_container,wdes_kbatch,vk,1.0_q,latt_cur,latt_ini)
      call get_pws_ibz(full,band_start,band_stop,ikibz,isp,w_dense_ibz,w_kbatch,&
                       cache_cw,wavwin_cw,wave_map,wdes1_container%wdes1,wkfbz)
      call get_projections_ibz(full,band_start,band_stop,ikibz,isp,&
                               dproj_kbatch,wdes_kbatch,cache_cproj,wavwin_cproj,wave_map,dproj_fbz)
      PROFILING_STOP('get_wfs_ibz')
  end subroutine get_wfs_ibz

  !> get pw and projections at a specific k-point between band_start and band_stop
  !> if the k-point is not in the IBZ then it is rotated from its image in the IBZ
  subroutine get_wfs(full,rotateprojectors,kpoints,band_start,band_stop,ikfbz,isp,w_dense_ibz,wdes_dense_ibz,w_kbatch,wdes_kbatch,dproj_kbatch,wave_map,&
                     wavwin_cw,cache_cw,wavwin_cproj,cache_cproj,waverot,latt_cur,latt_ini,p,t_info,wdes1_container,wkfbz,dproj_fbz)
      use wave_high, only: newwava, setwdes, element, delwava
      use wave_rotate, only: wave_rotator, wave_rotator_rotate_pw, wave_rotator_rotate_dproj, wavedes1_container, wdes1_container_init, rotate_vkpt_irot
      use mkpoints, only: kpoints_fbz_to_ibz, kindex_vkpt, ikibz_kindex
      use nonl_struct_def, only: nonl_struct
      use nonl_high, only: spher, phase, nonl_alloc, nonl_dealloc
      use wap, only: get_dproj
      use tutor, only: vtutor
      logical,intent(in) :: full !< wether to gather all the bands in the current MPI rank
      logical,intent(in) :: rotateprojectors !< wether to rotate the projectors and their derivatives from the IBZ to the FBZ
                                             !< alternatively their are simply recomputed
      type(kpoints_struct),intent(in) :: kpoints
      integer :: band_start !< index of the first band which we want to obtain
      integer :: band_stop !< index of the last band which we want to obtain
      integer :: isp !< index of the spin channel
      integer :: ikfbz !< index of the k-point in the FBZ we want to obtain
      GDEF,intent(in) :: dproj_kbatch(:,0:,:,:,:,:) !< container for the projections and their ionic derivaties in the IBZ
      type(wavespin),intent(in) :: w_dense_ibz !< wavefunctions in the IBZ
      type(wavedes),intent(in) :: wdes_dense_ibz !< descriptor for the WFs in the IBZ
      type(wavespin),intent(in) :: w_kbatch !< local batch of wavefunctions
      type(wavedes),intent(in) :: wdes_kbatch !< descriptor of the local batch of wavefunctions
      type(wave_mapper) :: wave_map !< mapping of the k-points in the k-point batch to the IBZ
      type(wave_window),intent(in) :: wavwin_cw !< wave window for the PW coefficients
      type(wave_cache),intent(in) :: cache_cw !< cache for the PW coefficients
      type(wave_window),intent(in) :: wavwin_cproj !< wave window for the projections
      type(wave_cache),intent(in) :: cache_cproj !< cahce for the projections
      type(wave_rotator),intent(in) :: waverot
      type(latt),intent(in) :: latt_cur
      type(latt),intent(in) :: latt_ini
      type(potcar),intent(in) :: p(:)
      type(type_info),intent(in) :: t_info
      type(wavedes1_container) :: wdes1_container !< container for the wavefunction descriptor at ikfbz
      type(wavefuna) :: wkfbz !< container for the PW coefficients in the FBZ point
      GDEF :: dproj_fbz(:,:,0:,:) !< container for the projections and their ionic directional derivatives in the FBZ point
      ! local variables
      type(nonl_struct) :: nonl_s
      type(wavefuna) :: wkibz
      type(wavedes1) :: wdesk
      integer :: ikibz, irot, nbk, nbk_global
      real(q) :: vkfbz_rotated(3)
      GDEF,allocatable :: dproj_ibz(:,:,:,:)
      GDEF,allocatable :: dproj(:,:,:)

      PROFILING_START('get_wfs')
      ! get ibz corresponding to fbz index
      call kpoints_fbz_to_ibz(kpoints,ikfbz,ikibz,irot)
      ! get wavefunction descriptor at a certain k-point coordinate in the fbz
      ! TODO: in principle I could use vkfbz or vkfbz_rotated but sometimes they are not the same
      !vkfbz = kindex_vkpt(kpoints,latt_cur,ikfbz)
      vkfbz_rotated = rotate_vkpt_irot(kpoints%vkpt(:,ikibz),waverot%igrpop,irot)

      call wdes1_container_init(wdes1_container,wdes_kbatch,vkfbz_rotated,1.0_q,latt_cur,latt_ini)
      if (ikfbz==ikibz_kindex(kpoints,latt_cur,ikibz)) then
          if (irot/=1) call vtutor%bug('wrong rotation index',__FILE__,__LINE__)
          ! if the current point is in the IBZ we just need to get the information
          call get_pws_ibz(full,band_start,band_stop,ikibz,isp,w_dense_ibz,w_kbatch,&
                           cache_cw,wavwin_cw,wave_map,wdes1_container%wdes1,wkfbz)
          call get_projections_ibz(full,band_start,band_stop,ikibz,isp,&
                                   dproj_kbatch,wdes_kbatch,cache_cproj,wavwin_cproj,wave_map,dproj_fbz)
      else
          ! if not then we get the wf in the IBZ
          call setwdes(wdes_dense_ibz,wdesk,ikibz)
          call get_pws_ibz(full,band_start,band_stop,ikibz,isp,w_dense_ibz,w_kbatch,&
                           cache_cw,wavwin_cw,wave_map,wdesk,wkibz)
          ! allocate wavefuna at ikfbz
          wdes1_container%wdes1%nbands = wdesk%nbands
          call newwava(wkfbz,wdes1_container%wdes1,wdesk%nbands)
          ! rotate it to the FBZ
          call wave_rotator_rotate_pw(waverot,wkibz,wkfbz,irot,latt_cur)
          wkfbz%celen=wkibz%celen
          ! there are two options to get the projections in the FBZ,
          ! one is to rotate them from the IBZ and the other one is to rotate the plane-wave part and recompute them.
          ! For the moment we recompute because rotating from the IBZ is not correctly implemented.
          if (rotateprojectors) then
              ! get projections and its derivatives in the IBZ ...
              allocate(dproj_ibz(size(dproj_fbz,1),size(dproj_fbz,2),0:3,size(dproj_fbz,4)))
              call get_projections_ibz(full,band_start,band_stop,ikibz,isp,&
                                       dproj_kbatch,wdes_kbatch,cache_cproj,wavwin_cproj,wave_map,dproj_ibz)
              ! ...and rotate to the FBZ
              call wave_rotator_rotate_dproj(waverot,wkibz,dproj_ibz,wkfbz,dproj_fbz,irot,latt_cur,p)
              deallocate(dproj_ibz)
          else
              allocate(dproj(size(dproj_fbz,1),0:3,size(dproj_fbz,4)))
              ! create non-local projector object
              !if (info%lreal) then
              !  call vtutor%error('lreal=.true. or lreal=auto not implemented')
              !  call nonlr_setup(nonlr_s,t_info,p,info%lreal,wdes_dense_ibz_star%lspiral)
              !  call nonlr_alloc(nonlr_s)
              !  call rspher(wdes_dense_ibz%grid,nonlr_s,latt_cur)
              !else
                call nonl_alloc(nonl_s,t_info,p,wdes1_container%wdes,.false.)
                call spher(wdes_dense_ibz%grid,nonl_s,p,wdes1_container%wdes,latt_cur, 1)
                call phase(wdes1_container%wdes,nonl_s,0)
              !endif

              ! Compute the derivatives of the projectors for k'
              call phase(wdes1_container%wdes,nonl_s,wdes1_container%wdes1%nk)
              !$omp parallel do private(nbk,dproj)
              do nbk_global=band_start,band_stop
                  if (mod(nbk_global-1,wdes1_container%wdes%nb_par)+1 /= wdes1_container%wdes%nb_low) cycle
                  ! compute local band index
                  nbk=1+(nbk_global-1)/wdes1_container%wdes%nb_par
                  call get_dproj(latt_cur,nonl_s,&
                                 element(wkfbz,nbk-band_start+1),wdes1_container%wdes1,dproj)
                  dproj_fbz(:,nbk-band_start+1,:,:) = dproj(:,:,:)
              enddo

              ! deallocate this for now (new allocations might be done inside the loop)
              call nonl_dealloc(nonl_s)
              deallocate(dproj)
          endif
          call delwava(wkibz)
      endif
      PROFILING_STOP('get_wfs')
  end subroutine get_wfs

  !> Fourier transform the WFs to real space
  subroutine fft_wfs(wdesk,nrspinors,cw,cr)
      use tutor, only: vtutor
      use string, only: str
      type(wavedes1),intent(in) :: wdesk
      integer,intent(in) :: nrspinors
      complex(q),intent(in) :: cw(:,:)
      complex(q) :: cr(:,:,:)
      !local variables
      integer :: nb,ispinor,ipwstart
      complex(q) :: cr_tmp(wdesk%grid%mplwv)
      if (size(cr,1)/=wdesk%grid%rl%np) then
          call vtutor%bug('Wrong dimensions of CR to store the real space WFs. '//&
                          'Got '//str(size(cr,1))//' expected '//str(wdesk%grid%rl%np),__FILE__,__LINE__)
      endif
      do ispinor=1,nrspinors
          ipwstart = 1+(ispinor-1)*wdesk%ngvector
          do nb=1,wdesk%nbands
              call FFTWAV(wdesk%ngvector, wdesk%nindpw, &
                          cr_tmp(:), &
                          cw(ipwstart,nb), wdesk%grid)
              !copy the data
              cr(:,ispinor,nb) = cr_tmp(:wdesk%grid%rl%np)
          enddo ! loop over bands
      enddo ! loop over spinors
  end subroutine fft_wfs

  !> Fourier transform the WFs to real space
  subroutine fft_nabla_wfs(wdesk,latt_cur,nrspinors,cw,cr)
    use constant, only: citpi
    use tutor, only: vtutor
    use string, only: str
    type(wavedes1),intent(in) :: wdesk
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: nrspinors
    complex(q),intent(in) :: cw(:,:)
    complex(q) :: cr(:,:,:,:)
    !local variables
    integer :: nb,ispinor,ipwstart,idir,ipw
    complex(q) :: cr_tmp(wdesk%grid%mplwv)
    complex(q) :: cwork(wdesk%ngvector)
    real(q) :: g

    if (size(cr,1)/=wdesk%grid%rl%np) then
        call vtutor%bug('Wrong dimensions of CR to store the real space WFs. '//&
                        'Got '//str(size(cr,1))//' expected '//str(wdesk%grid%rl%np),__FILE__,__LINE__)
    endif

    do idir=1,3
        do nb=1,wdesk%nbands
            do ispinor=1,nrspinors
               ipwstart = (ispinor-1)*wdesk%ngvector
               do ipw=1,wdesk%ngvector
                  !multiply by ctpi*g
                  g=((wdesk%igx(ipw)+wdesk%vkpt(1))*latt_cur%b(idir,1)+&
                     (wdesk%igy(ipw)+wdesk%vkpt(2))*latt_cur%b(idir,2)+&
                     (wdesk%igz(ipw)+wdesk%vkpt(3))*latt_cur%b(idir,3))
                  cwork(ipw) = cw(ipwstart+ipw,nb)*g*citpi
               enddo
               call FFTWAV(wdesk%ngvector, wdesk%nindpw, &
                           cr_tmp(:), &
                           cwork, wdesk%grid)
               !copy the data
               cr(:,ispinor,nb,idir) = cr_tmp(:wdesk%grid%rl%np)
            enddo !loop over spinors
        enddo ! loop over bands
    enddo ! loop over directions
  end subroutine fft_nabla_wfs

  !> Compute electron-phonon matrix elements in the plane-wave basis
  !> between specific k-points in the IBZ. This basically skips rotating the k-points
  subroutine elph_mels_pw_compute_g_ibz(self,elph_set,kpoints,ikibz,ikpibz,isp,vqpt,gcart_kkp)
      use mkpoints, only: kindex_vkpt, kpoint_integers, vkpt_kindex, kpoints_fbz_to_ibz
      use wave_high, only: delwava
      use wave_rotate, only: wavedes1_container, wdes1_container_free, rotate_vkpt_irot
      use elphon_potential, only: elphon_pot_get_dvdu_prim, elphon_pot_get_ddijdu_prim
      type(elph_mels_pw_type) :: self
      type(elph_settings),intent(in) :: elph_set !< electron-phonon set
      type(kpoints_struct),intent(in) :: kpoints
      integer,intent(in) :: ikibz !<index of the k-point in the fbz
      integer,intent(in) :: ikpibz !<index of the k'-point in the fbz
      integer,intent(in) :: isp
      real(q),intent(out) :: vqpt(3) !< q-point in the reciprocal space
      complex(q),intent(out) :: gcart_kkp(:,:,:,:)
      ! local variables
      real(q) :: vk(3), vkp(3)
      type(wavedes1_container) :: wdesk_container
      type(wavedes1_container) :: wdeskp_container
      type(wavefuna) :: wk_fbz
      type(wavefuna) :: wkp_fbz

      vk = kpoints%vkpt(:,ikibz)
      vkp = kpoints%vkpt(:,ikpibz)
      vqpt = vk-vkp

      ! early exit if computing g is not required
      if (elph_set%selfen_g_skip) then
        gcart_kkp = 1.0_q
        return
      endif

      PROFILING_START('elph_compute_g')
      ! get wf at ikibz
      call get_wfs_ibz(.true.,kpoints,self%band_start_k,self%band_stop_k,ikibz,isp,&
                       self%w_dense_ibz,self%w_kbatch,self%wdes_kbatch,self%dproj_kbatch,self%wave_map,&
                       self%wavwin_cw,self%cache_cw,self%wavwin_cproj,self%cache_cproj,self%latt_cur,self%latt_ini,&
                       wdesk_container,wk_fbz,self%dprojk)

      ! get wf at ikpibz
      call get_wfs_ibz(.false.,kpoints,self%band_start_kp,self%band_stop_kp,ikpibz,isp,&
                       self%w_dense_ibz,self%w_kbatch,self%wdes_kbatch,self%dproj_kbatch,self%wave_map,&
                       self%wavwin_cw,self%cache_cw,self%wavwin_cproj,self%cache_cproj,self%latt_cur,self%latt_ini,&
                       wdeskp_container,wkp_fbz,self%dprojkp)

      ! compute electron-phonon potential
      if (elph_set%lgaddv) &
          call elphon_pot_get_dvdu_prim(self%elph_pot,vqpt,self%elph_pot%has_lr_active,self%elph_pot%dvdu,self%dvdu_kkp)
      if (elph_set%lgaddt) &
          call elphon_pot_get_dvdu_prim(self%elph_pot,vqpt,.false.,self%elph_pot%dmudu,self%dmudu_kkp)
      if (elph_set%lgaddd) &
          call elphon_pot_get_ddijdu_prim(self%elph_pot,vqpt,self%ddijdu_kkp)

      ! compute electron-phonon matrix elements
      call elph_sandwich(elph_set,self%elph_pot,self%elph_caches,self%latt_cur,&
                         self%lmdim,self%nrspinors,self%nproj,self%nproj_tot,self%npw,self%ntyp,self%nityp,&
                         self%lmmax,self%nbands_k,self%nbands_kp,isp,&
                         wk_fbz,[0.0_q,0.0_q,0.0_q],self%w_dense_ibz%celtot(self%band_start_k:,ikibz,isp),self%dprojk,&
                         wkp_fbz,wkp_fbz%celen,self%dprojkp,&
                         vqpt,self%dvdu_kkp,self%dmudu_kkp,self%ddijdu_kkp,gcart_kkp)
      CALLMPI( M_sum_z(self%wdes_dense_ibz%comm_inb, gcart_kkp, size(gcart_kkp)) )

      ! cleanup
      call wdes1_container_free(wdesk_container)
      call wdes1_container_free(wdeskp_container)
      call delwava(wk_fbz)
      call delwava(wkp_fbz)
      PROFILING_STOP('elph_compute_g')

  end subroutine elph_mels_pw_compute_g_ibz

  !> Compute electron-phonon matrix elements for ik and ikp pair of orbitals
  subroutine elph_mels_pw_compute_g(self,elph_set,kpoints,ikfbz,ikpfbz,iqfbz,isp,vqpt,gcart_kkp)
      use tutor, only: vtutor
      use mkpoints, only: kindex_vkpt, kpoint_integers, vkpt_kindex, kpoints_fbz_to_ibz
      use wave_high, only: delwava
      use wave_rotate, only: wavedes1_container, wdes1_container_free, rotate_vkpt_irot
      use elphon_potential, only: elphon_pot_get_dvdu_prim, elphon_pot_get_ddijdu_prim
      use string, only: str
      type(elph_mels_pw_type) :: self
      type(elph_settings),intent(in) :: elph_set !< electron-phonon settings
      type(kpoints_struct),intent(in) :: kpoints
      integer,intent(in) :: ikfbz !<index of he k-point in the fbz
      integer,intent(in) :: ikpfbz !<index of the k'-point in the fbz
      integer,intent(in) :: iqfbz !<index of the q-point in the fbz
      integer,intent(in) :: isp !< index of the spin channel
      real(q) :: vqpt(3) !< q-point in the reciprocal space
      complex(q),intent(out) :: gcart_kkp(:,:,:,:) !< electron-phonon matrix elements in the cartesian basis
      ! local variables
      integer :: ikibz, ikpibz, iqfbz_local, irotp
      real(q) :: vqpt_pot(3), vk(3), vkp(3)
      real(q) :: qdiff(3)
      type(wavedes1_container) :: wdesk_container
      type(wavedes1_container) :: wdeskp_container
      type(wavefuna) :: wk_fbz
      type(wavefuna) :: wkp_fbz

      ! get indexes of kpoints
      ikibz = kpoints%fbz2ibz(1,ikfbz)
      vk = kpoints%vkpt(:,ikibz)
      call kpoints_fbz_to_ibz(kpoints,ikpfbz,ikpibz,irotp)
      vkp = rotate_vkpt_irot(kpoints%vkpt(:,ikpibz),self%waverot%igrpop,irotp)

      vqpt = vk-vkp
      iqfbz_local = vkpt_kindex(kpoints,self%latt_cur,vqpt)
      ! sanity check
      if (iqfbz/=iqfbz_local) then
        call vtutor%error('Wrong index of the q-point index '//str(iqfbz)//' instead of '//str(iqfbz_local)//&
                          ' for computation of electron-phonon matrix element, computed ')
      endif

      ! early exit if computing g is not required
      if (elph_set%selfen_g_skip) then
        gcart_kkp = 1.0_q
        return
      endif

      PROFILING_START('elph_compute_g')
      ! get wf at ikfbz
      call get_wfs(.true.,elph_set%rotateprojectors,kpoints,self%band_start_k,self%band_stop_k,ikfbz,isp,&
                   self%w_dense_ibz,self%wdes_dense_ibz,self%w_kbatch,self%wdes_kbatch,self%dproj_kbatch,self%wave_map,&
                   self%wavwin_cw,self%cache_cw,self%wavwin_cproj,self%cache_cproj,self%waverot,self%latt_cur,self%latt_ini,self%p,self%t_info,&
                   wdesk_container,wk_fbz,self%dprojk)

      ! get wf at ikpfbz
      call get_wfs(.false.,elph_set%rotateprojectors,kpoints,self%band_start_kp,self%band_stop_kp,ikpfbz,isp,&
                   self%w_dense_ibz,self%wdes_dense_ibz,self%w_kbatch,self%wdes_kbatch,self%dproj_kbatch,self%wave_map,&
                   self%wavwin_cw,self%cache_cw,self%wavwin_cproj,self%cache_cproj,self%waverot,self%latt_cur,self%latt_ini,self%p,self%t_info,&
                   wdeskp_container,wkp_fbz,self%dprojkp)

      qdiff = [0.0_q,0.0_q,0.0_q]
      if (elph_set%use_pot_cache) then
          ! Here we use a caching algorithm for computing the electron-phonon potential
          ! this works by checking if the q-point at which we need to interpolate
          ! the potential is the same one as before. If that is the case then we don't need
          ! to interpolate it again. This can lead to a speedup of 10x in calculations
          ! where the interpolation of the potential is the bottleneck.
          ! we use the fact that:
          !   q = k-k'
          ! and
          !   q'+ G = q
          ! this means that
          !   q'= (k-G)-k'
          ! if we computed the electron-phonon matrix elements using
          !   g(k,k',q)
          ! then if we can compute for q' as
          !   g(k-G,k',q')

          ! if we are dealing with a new q-point interpolate the potential again
          if (self%iq_cache/=iqfbz) then
              if (elph_set%lgaddv) &
                  call elphon_pot_get_dvdu_prim(self%elph_pot,vqpt,self%elph_pot%has_lr_active,self%elph_pot%dvdu,self%dvdu_kkp)
              if (elph_set%lgaddt) &
                  call elphon_pot_get_dvdu_prim(self%elph_pot,vqpt,.false.,self%elph_pot%dmudu,self%dmudu_kkp)
              if (elph_set%lgaddd) &
                  call elphon_pot_get_ddijdu_prim(self%elph_pot,vqpt,self%ddijdu_kkp)
              self%vqpt = vqpt
              self%iq_cache = iqfbz
          endif

          ! multiply cproj part of the WF at k by phase factor
          ! the plane-wave part is mutiplied inside the elph_sandwich routine
          qdiff = vqpt-self%vqpt
          call multiply_phase_cproj(wk_fbz%wdes1,self%elph_pot%centers,qdiff,self%dprojk)
          vqpt_pot = self%vqpt
      else
          if (elph_set%lgaddv) &
              call elphon_pot_get_dvdu_prim(self%elph_pot,vqpt,self%elph_pot%has_lr_active,self%elph_pot%dvdu,self%dvdu_kkp)
          if (elph_set%lgaddt) &
              call elphon_pot_get_dvdu_prim(self%elph_pot,vqpt,.false.,self%elph_pot%dmudu,self%dmudu_kkp)
          if (elph_set%lgaddd) &
              call elphon_pot_get_ddijdu_prim(self%elph_pot,vqpt,self%ddijdu_kkp)
          vqpt_pot = vqpt
      endif

      ! compute electron-phonon matrix elements
      call elph_sandwich(elph_set,self%elph_pot,self%elph_caches,self%latt_cur,&
                         self%lmdim,self%nrspinors,self%nproj,self%nproj_tot,self%npw,self%ntyp,self%nityp,&
                         self%lmmax,self%nbands_k,self%nbands_kp,isp,&
                         wk_fbz,qdiff,self%w_dense_ibz%celtot(self%band_start_k:,ikibz,isp),self%dprojk,&
                         wkp_fbz,wkp_fbz%celen,self%dprojkp,&
                         vqpt_pot,self%dvdu_kkp,self%dmudu_kkp,self%ddijdu_kkp,gcart_kkp)
      CALLMPI( M_sum_z(self%wdes_dense_ibz%comm_inb, gcart_kkp, size(gcart_kkp)) )

      ! cleanup
      call wdes1_container_free(wdesk_container)
      call wdes1_container_free(wdeskp_container)
      call delwava(wk_fbz)
      call delwava(wkp_fbz)
      PROFILING_STOP('elph_compute_g')

  end subroutine elph_mels_pw_compute_g

  subroutine get_phase_pw(grid,npw_local,vqpt,phase)
      use constant, only: citpi
      type(grid_3d) :: grid
      integer :: npw_local
      complex(q) :: phase(npw_local)
      real(q) :: vqpt(3)
      ! local variables
      real(q) :: r(3)
      integer :: nx,ny,nz,nc
      integer :: i,ibands,ispinor,nbands,nrspinor

      PROFILING_START('get_phase_pw')
      ! electron-phonon potential
      i = 0
      if (grid%rl%nfast==1) then
        do nc=1,grid%rl%ncol
           ny=grid%rl%i2(nc)-1
           nz=grid%rl%i3(nc)-1
           r(2) = 1.0_q*ny/grid%ngy
           r(3) = 1.0_q*nz/grid%ngz
           do nx=0,grid%rl%nrow-1
              r(1) = 1.0_q*nx/grid%ngx
              i=i+1
              phase(i) = exp(citpi*dot_product(vqpt,r))
           enddo
        enddo
      elseif (grid%rl%nfast==3) then
        do nc=1,grid%rl%ncol
           nx=grid%rl%i2(nc)-1
           ny=grid%rl%i3(nc)-1
           r(1) = 1.0_q*nx/grid%ngx
           r(2) = 1.0_q*ny/grid%ngy
           do nz=0,grid%rl%nrow-1
              r(3) = 1.0_q*nz/grid%ngz
              i=i+1
              phase(i) = exp(citpi*dot_product(vqpt,r))
           enddo
        enddo
      endif
      PROFILING_STOP('get_phase_pw')
  end subroutine get_phase_pw

  subroutine multiply_phase_pw(grid,npw_local,vqpt,cr)
      use constant, only: citpi
      type(grid_3d) :: grid
      integer :: npw_local
      real(q) :: vqpt(3)
      complex(q) :: cr(:,:,:)
      ! local variables
      complex(q) :: phase(npw_local)
      integer :: ibands,ispinor,nbands,nrspinor
      call get_phase_pw(grid,npw_local,vqpt,phase)

      PROFILING_START('multiply_phase_pw')
      nbands = size(cr,2)
      nrspinor = size(cr,3)
      do ispinor=1,nrspinor
        do ibands=1,nbands
          cr(:,ibands,ispinor) = cr(:,ibands,ispinor)*phase(:)
        enddo
      enddo
      PROFILING_STOP('multiply_phase_pw')
  end subroutine multiply_phase_pw

  subroutine multiply_phase_cproj(wdes1,centers,vqpt,dprojk)
      use constant, only: citpi
      use wave, only: ni_global
      type(wavedes1) :: wdes1
      real(q),intent(in) :: centers(:,:)
      real(q),intent(in) :: vqpt(3)
      GDEF :: dprojk(:,:,0:,:)
      ! local variables
      integer :: lmbase_at, lmmax, ibands, iproj, iatom_local, ii, ispinor, iatom_global
      integer :: nbands, nrspinor, nt
      real(q) :: angle_q
      complex(q) :: phase_p
      ! multiply projections by phase factors
      lmbase_at=0
      iatom_local=0
      nbands = size(dprojk,2)
      nrspinor = size(dprojk,4)
      ! loop over atoms
      do nt = 1, wdes1%ntyp
      do ii = 1, wdes1%nityp(nt)
         lmmax = wdes1%lmmax(nt)
         iatom_local = iatom_local+1
         iatom_global = ni_global(iatom_local,wdes1%comm_inb)

         angle_q=dot_product(vqpt,centers(:,iatom_global))
         phase_p=exp(citpi*angle_q)

         do ispinor=1,nrspinor
           do ibands=1,nbands
             do iproj=1,lmmax
               dprojk(lmbase_at+iproj,ibands,:,ispinor) = dprojk(lmbase_at+iproj,ibands,:,ispinor)*phase_p
             enddo
           enddo
         enddo
         lmbase_at = lmbase_at + lmmax
      enddo
      enddo
  end subroutine multiply_phase_cproj

  subroutine elph_mels_pw_pot_free(self)
      use elphon_potential, only: elphon_pot_free
      type(elph_mels_pw_type) :: self
      call elphon_pot_free(self%elph_pot)
  end subroutine elph_mels_pw_pot_free

  subroutine elph_mels_pw_wfs_free(self,kpoints,w,nonl_s)
      use nonl_struct_def, only: nonl_struct
      use wave_interpolate, only: wave_interpolator_free
      use wave_rotate, only: wave_rotator_free, fbz2ibz_stars_free
      use wave, only: deallocw, deallocw_noplanewave
      use ini, only: deregister_allocate
      type(elph_mels_pw_type) :: self
      type(kpoints_struct) :: kpoints
      type(wavespin) :: w
      type(nonl_struct) :: nonl_s
      call wave_interpolator_free(self%wave_interp, kpoints, self%grid, &
                    self%latt_cur, self%latt_ini, self%t_info, nonl_s, w, self%p, self%symm)
      call fbz2ibz_stars_free(self%fbz2ibz)
      call wave_rotator_free(self%waverot)
      call deallocw(self%w_kbatch)
      call deregister_allocate(VMEM_ARR_SIZE(self%dproj_kbatch), "dproj")
      deallocate(self%dproj_kbatch)
      call deregister_allocate(VMEM_ARR_SIZE(self%w_dense_ibz%celtot)+&
                               VMEM_ARR_SIZE(self%w_dense_ibz%fertot)+&
                               VMEM_ARR_SIZE(self%w_dense_ibz%auxtot),"elph_w_dense")
      call deallocw_noplanewave(self%w_dense_ibz)
  end subroutine elph_mels_pw_wfs_free

  !> High level routine to compute the sandwich.
  !> The contributions to the electron phonon matrix element that are computed
  !> are the ones selected in elph_set. It translates to a call to elph_sandwich_low.
  subroutine elph_sandwich(elph_set,elph_pot,elph_caches,latt_cur,&
                           lmdim,nrspinors,nproj,nproj_tot,npw,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                           wk_fbz,kphase,celenk,dprojk,&
                           wkp_fbz,celenkp,dprojkp,&
                           vqpt,dvdu_kkp,dmudu_kkp,ddijdu_kkp,gcart_kkp)
     use elphon_base, only: elph_settings, elph_caches_t
     type (elph_settings)  :: elph_set !< settings for the electron-phonon calculation
     type (elphon_pot) :: elph_pot !< electron-phonon potential
     type (elph_caches_t)  :: elph_caches !< caches for the electron-phonon calculation
     integer,intent(in) :: lmdim
     integer,intent(in) :: nrspinors
     integer,intent(in) :: nproj
     integer,intent(in) :: nproj_tot
     integer,intent(in) :: npw
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: lmmax_all(:)
     integer,intent(in) :: nbands_k !< number of bands at k
     integer,intent(in) :: nbands_kp !< number of bands at k'
     integer,intent(in) :: isp !< spin index
     type (latt) :: latt_cur !< lattice structure
     type (wavefuna),intent(in) :: wk_fbz  !< wavefunction at k
     type (wavefuna),intent(in) :: wkp_fbz  !< wavefunction at k
     real(q),intent(in) :: kphase(3) !< additional phase to mutiply at k
     complex(q),intent(in) :: celenk(:)
     complex(q),intent(in) :: celenkp(:)
     GDEF,intent(in) :: dprojk(:,:,0:,:)
     GDEF,intent(in) :: dprojkp(:,:,0:,:)
     real(q) :: vqpt(3)
     complex(q),intent(in) :: dvdu_kkp(:,:,:,:)
     complex(q),intent(in) :: dmudu_kkp(:,:,:,:)
     complex(q),intent(in) :: ddijdu_kkp(:,:,:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:)

     call elph_sandwich_low(elph_set,elph_pot,elph_caches,latt_cur,&
                            lmdim,nrspinors,nproj,nproj_tot,npw,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                            elph_set%lgaddr,elph_set%lgaddv,elph_set%lgaddt,elph_set%lgaddp,&
                            elph_set%lgadds,elph_set%lgaddq,elph_set%lgaddd,&
                            wk_fbz,kphase,celenk,dprojk,&
                            wkp_fbz,celenkp,dprojkp,&
                            vqpt,dvdu_kkp,dmudu_kkp,ddijdu_kkp,gcart_kkp)

  end subroutine elph_sandwich

  !> High level routine to compute the sandwich.
  !> Only the additional gs contributions are computed.
  !> It translates to a call to elph_sandwich_low.
  subroutine elph_sandwich_gs(elph_set,elph_pot,elph_caches,latt_cur,&
                              lmdim,nrspinors,nproj,nproj_tot,npw,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                              wk_fbz,kphase,celenk,dprojk,&
                              wkp_fbz,celenkp,dprojkp,&
                              vqpt,dvdu_kkp,dmudu_kkp,ddijdu_kkp,gcart_kkp)
     use elphon_base, only: elph_settings, elph_caches_t
     type (elph_settings)  :: elph_set !< settings for the electron-phonon calculation
     type (elphon_pot) :: elph_pot !< electron-phonon potential
     type (elph_caches_t)  :: elph_caches !< caches for the electron-phonon calculation
     integer,intent(in) :: lmdim
     integer,intent(in) :: nrspinors
     integer,intent(in) :: nproj
     integer,intent(in) :: nproj_tot
     integer,intent(in) :: npw
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: lmmax_all(:)
     integer,intent(in) :: nbands_k !< number of bands at k
     integer,intent(in) :: nbands_kp !< number of bands at k'
     integer,intent(in) :: isp !< spin index
     type (latt) :: latt_cur
     type (wavefuna),intent(in) :: wk_fbz  !< wavefunction at k
     type (wavefuna),intent(in) :: wkp_fbz  !< wavefunction at k
     real(q),intent(in) :: kphase(3)
     complex(q),intent(in) :: celenk(:)
     complex(q),intent(in) :: celenkp(:)
     GDEF,intent(in) :: dprojk(:,:,0:,:)
     GDEF,intent(in) :: dprojkp(:,:,0:,:)
     real(q),intent(in) :: vqpt(3)
     complex(q),intent(in) :: dvdu_kkp(:,:,:,:)
     complex(q),intent(in) :: dmudu_kkp(:,:,:,:)
     complex(q),intent(in) :: ddijdu_kkp(:,:,:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:)
     ! local variables
     logical :: lgaddr,lgaddt,lgaddv,lgaddp,lgadds,lgaddq,lgaddd
     lgaddr=.false.
     lgaddt=.false.
     lgaddv=.false.
     lgaddp=.false.
     lgadds=.true.
     lgaddq=.false.
     lgaddd=.false.
     call elph_sandwich_low(elph_set,elph_pot,elph_caches,latt_cur,&
                            lmdim,nrspinors,nproj,nproj_tot,npw,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                            lgaddr,lgaddv,lgaddt,lgaddp,lgadds,lgaddq,lgaddd,&
                            wk_fbz,kphase,celenk,dprojk,&
                            wkp_fbz,celenkp,dprojkp,&
                            vqpt,dvdu_kkp,dmudu_kkp,ddijdu_kkp,gcart_kkp)
  end subroutine elph_sandwich_gs

  !> Low level routine to compute the sandwich
  subroutine elph_sandwich_low(elph_set,elph_pot,elph_caches,latt_cur,&
                               lmdim,nrspinors,nproj,nproj_tot,npw,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                               lgaddr,lgaddv,lgaddt,lgaddp,lgadds,lgaddq,lgaddd,&
                               wk_fbz,kphase,celenk,dprojk,&
                               wkp_fbz,celenkp,dprojkp,&
                               vqpt,dvdu_kkp,dmudu_kkp,ddijdu_kkp,gcart_kkp)
     use elphon_base, only: elph_settings, elph_caches_t
     type (elph_settings)  :: elph_set !< settings for the electron-phonon calculation
     type (elphon_pot) :: elph_pot !< electron-phonon potential
     type (elph_caches_t)  :: elph_caches !< caches for the electron-phonon calculation
     integer,intent(in) :: lmdim
     integer,intent(in) :: nrspinors
     integer,intent(in) :: nproj
     integer,intent(in) :: nproj_tot
     integer,intent(in) :: npw
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: lmmax_all(:)
     integer,intent(in) :: nbands_k !< number of bands at k
     integer,intent(in) :: nbands_kp !< number of bands at k'
     integer,intent(in) :: isp !< spin index
     logical,intent(in) :: lgaddr,lgaddv,lgaddt,lgaddp,lgadds,lgaddq,lgaddd
     type (latt) :: latt_cur
     type (wavefuna) :: wk_fbz  !< wavefunction at k
     type (wavefuna) :: wkp_fbz  !< wavefunction at k
     real(q),intent(in) :: kphase(3)
     complex(q),intent(in) :: celenk(:)
     complex(q),intent(in) :: celenkp(:)
     GDEF,intent(in) :: dprojk(:,:,0:,:)
     GDEF,intent(in) :: dprojkp(:,:,0:,:)
     real(q) :: vqpt(3)
     complex(q),intent(in) :: dvdu_kkp(:,:,:,:)
     complex(q),intent(in) :: dmudu_kkp(:,:,:,:)
     complex(q),intent(in) :: ddijdu_kkp(:,:,:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:)
     !local variables
     integer :: idir

     if (elph_set%selfen_g_skip) then
         gcart_kkp = (1.0_q,0.0_q)
         return
     endif

     PROFILING_START('elph_sandwich')
     gcart_kkp = 0
     if (elph_set%userecip) then
         call elph_sandwich_pw_recip(lgaddv,lgaddt,latt_cur,&
                                 elph_caches,nrspinors,npw,ntyp,nityp,nbands_k,nbands_kp,isp,&
                                 wk_fbz,kphase,wkp_fbz,&
                                 vqpt,elph_pot%centers,elph_pot%primitive_positions,&
                                 dvdu_kkp,dmudu_kkp,gcart_kkp)
         call elph_sandwich_blas(lgaddr,lgaddp,lgadds,lgaddq,lgaddd,&
                                 elph_caches,lmdim,nrspinors,nproj,nproj_tot,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                                 wk_fbz,celenk,dprojk,celenkp,dprojkp,&
                                 vqpt,elph_pot%centers,elph_pot%primitive_positions,&
                                 ddijdu_kkp,elph_caches%ccdij,elph_caches%ccqij,&
                                 elph_caches%ccdrij,gcart_kkp)
     else if (elph_set%useblas) then
         call elph_sandwich_pw_blas(lgaddv,&
                                 elph_caches,nrspinors,npw,ntyp,nityp,nbands_k,nbands_kp,isp,&
                                 wk_fbz,kphase,wkp_fbz,&
                                 vqpt,elph_pot%centers,elph_pot%primitive_positions,&
                                 dvdu_kkp,gcart_kkp)
         call elph_sandwich_blas(lgaddr,lgaddp,lgadds,lgaddq,lgaddd,&
                                 elph_caches,lmdim,nrspinors,nproj,nproj_tot,ntyp,nityp,lmmax_all,nbands_k,nbands_kp,isp,&
                                 wk_fbz,celenk,dprojk,celenkp,dprojkp,&
                                 vqpt,elph_pot%centers,elph_pot%primitive_positions,&
                                 ddijdu_kkp,elph_caches%ccdij,elph_caches%ccqij,&
                                 elph_caches%ccdrij,gcart_kkp)
     else
         call elph_sandwich_pw_ref(lgaddv,&
                                nrspinors,npw,ntyp,nityp,nbands_k,nbands_kp,isp,&
                                wk_fbz,wkp_fbz,&
                                vqpt,kphase,elph_pot%centers,elph_pot%primitive_positions,&
                                dvdu_kkp,gcart_kkp)
         call elph_sandwich_ref(lgaddr,lgaddp,lgadds,lgaddq,lgaddd,&
                                lmdim,nrspinors,nproj,nbands_k,nbands_kp,isp,&
                                wk_fbz,celenk,dprojk,celenkp,dprojkp,&
                                vqpt,elph_pot%centers,elph_pot%primitive_positions,&
                                ddijdu_kkp,elph_caches%ccdij,elph_caches%ccqij,&
                                elph_caches%ccdrij,gcart_kkp)
     endif
     PROFILING_STOP('elph_sandwich')

  end subroutine elph_sandwich_low

  !> Compute the electron-phonon matrix elements using the WFs at two different k-points
  !> This implementation is not the most efficient, it is used for reference
  !> The efficient version of this routine is elph_sandwich_blas
  subroutine elph_sandwich_pw_ref(lgaddv,&
                               nrspinors,npw,ntyp,nityp,nbands_k,nbands_kp,isp,&
                               wk_fbz,wkp_fbz,&
                               vqpt,kphase,centers,primitive_positions,&
                               dvdu_kkp,gcart_kkp)
     use constant
     use elphon_base, only: elph_settings
     logical,intent(in) :: lgaddv
     integer,intent(in) :: nrspinors
     integer,intent(in) :: npw
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: nbands_k
     integer,intent(in) :: nbands_kp
     integer,intent(in) :: isp
     type (wavefuna) :: wk_fbz  !< wavefunction at k
     type (wavefuna) :: wkp_fbz  !< wavefunction at k'
     real(q),intent(in) :: vqpt(3)
     real(q),intent(in) :: kphase(3) !< additional phase to mutiply at k
     real(q) :: centers(:,:)
     real(q) :: primitive_positions(:,:)
     complex(q),intent(in) :: dvdu_kkp(:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:)
     ! Local variables
     integer :: lmmax
     integer :: nbk, nbkp
     integer :: alpha
     integer :: spinors, spinors_k, spinors_kp
     integer :: atoms_disp
     integer :: nt,ii ! looping over atoms in the VASP way
     integer :: lmbase_at
     integer :: npr_local
     real(q) :: angle_q
     complex(q) :: phase_q
     complex(q) :: elph_coupling
     complex(q) :: elph_coupling_v
     complex(q) :: cr(wk_fbz%wdes1%grid%rl%np,nrspinors,wk_fbz%wdes1%nbands)
     complex(q) :: crp(wk_fbz%wdes1%grid%rl%np,nrspinors,wkp_fbz%wdes1%nbands)

     npr_local = wk_fbz%wdes1%grid%rl%np
     ! fft WFs
     call fft_wfs(wk_fbz%wdes1,nrspinors,wk_fbz%CW,cr)
     call fft_wfs(wkp_fbz%wdes1,nrspinors,wkp_fbz%CW,crp)
     if (any(abs(kphase(:))>1e-8)) then
        call multiply_phase_pw(wk_fbz%wdes1%grid,npr_local,kphase,cr)
     endif

     ! loop over spinors
     do spinors_k=1,nrspinors
     do spinors_kp=1,nrspinors
        lmbase_at=1
        spinors= isp+(spinors_kp-1)+2*(spinors_k-1)  ! = 1, 2, 3, 4

           atoms_disp=0
           ! loop over atoms
           do nt = 1, ntyp
           do ii = 1, nityp(nt)
              lmmax = wk_fbz%wdes1%lmmax(nt)
              atoms_disp=atoms_disp+1

              ! use Bloch theorem to get the coupling derivative for the atoms in POSCAR
              ! may not always be wise, so it is not included in the above routines
              angle_q=dot_product(vqpt,centers(:,atoms_disp)-primitive_positions(:,atoms_disp))
              phase_q=exp(citpi*angle_q)

              ! loop over directions
              do alpha =1,3

                 ! loop over bands
                 do nbk=1,nbands_k
                 do nbkp=1,nbands_kp

                    elph_coupling = cmplx(0.0_q,0.0_q,q)
                    if (lgaddv) then
                        PROFILING_START('elph_sandwich_gv')
                        call elph_gv(npw, &
                                     cr(:,spinors_k,nbk),crp(:,spinors_kp,nbkp), &
                                     dvdu_kkp(:,alpha,atoms_disp,spinors), &
                                     elph_coupling_v)
                        elph_coupling = elph_coupling + elph_coupling_v*phase_q
                        PROFILING_STOP('elph_sandwich_gv')
                    ENDIF

                    gcart_kkp(nbk,nbkp,alpha,atoms_disp)=&
                    gcart_kkp(nbk,nbkp,alpha,atoms_disp)+elph_coupling
                 enddo ! nbands k'
                 enddo ! nbands k

              enddo ! directions

              lmbase_at = lmbase_at + lmmax
           enddo ! loop over all atoms of type
           enddo ! loop over atom types
     enddo ! spinor k'
     enddo ! spinor k
  end subroutine elph_sandwich_pw_ref

  !> Compute the electron-phonon matrix elements using the WFs at two different k-points
  !> This implementation is not the most efficient, it is used for reference
  !> The efficient version of this routine is elph_sandwich_blas
  subroutine elph_sandwich_ref(lgaddr,lgaddp,lgadds,lgaddq,lgaddd,&
                               lmdim,nrspinors,nproj,nbands_k,nbands_kp,isp,&
                               wk_ibz,celenk,dprojk,&
                               celenkp,dprojkp,&
                               vqpt,centers,primitive_positions,&
                               ddijdu_kkp,ccdij,ccqij,ccdrij,gcart_kkp)
     use constant
     use wave, only: ni_global
     use elphon_base, only: elph_settings
     logical,intent(in) :: lgaddr,lgaddp,lgadds,lgaddq,lgaddd
     integer,intent(in) :: lmdim
     integer,intent(in) :: nrspinors
     integer,intent(in) :: nproj
     integer,intent(in) :: nbands_k
     integer,intent(in) :: nbands_kp
     integer,intent(in) :: isp
     type (wavefuna) :: wk_ibz  !< wavefunction at k
     complex(q),intent(in) :: celenk(:)
     complex(q),intent(in) :: celenkp(:)
     GDEF,intent(in) :: dprojk(:,:,0:,:)
     GDEF,intent(in) :: dprojkp(:,:,0:,:)
     real(q) :: vqpt(3)
     real(q) :: centers(:,:)
     real(q) :: primitive_positions(:,:)
     complex(q),intent(in) :: ddijdu_kkp(:,:,:,:,:,:)
     complex(q),intent(in) :: ccdij(:,:,:,:)
     complex(q),intent(in) :: ccqij(:,:,:,:)
     complex(q),intent(in) :: ccdrij(:,:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:)
     ! Local variables
     integer :: lmmax
     integer :: nbk, nbkp
     integer :: alpha
     integer :: spinors, spinors_k, spinors_kp
     integer :: iatom_local, iatom_global
     integer :: nt,ii ! looping over atoms in the VASP way
     integer :: lmbase_at
     real(q) :: angle_q
     complex(q) :: phase_q,phase_p
     complex(q) :: elph_coupling
     complex(q) :: elph_coupling_p,elph_coupling_r,elph_coupling_d,elph_coupling_q,elph_coupling_s

     ! loop over spinors
     do spinors_k=1,nrspinors
     do spinors_kp=1,nrspinors
        lmbase_at=1
        spinors= isp+(spinors_kp-1)+2*(spinors_k-1)  ! = 1, 2, 3, 4

           iatom_local=0
           ! loop over atoms
           do nt = 1, wk_ibz%wdes1%ntyp
           do ii = 1, wk_ibz%wdes1%nityp(nt)
              lmmax = wk_ibz%wdes1%lmmax(nt)
              iatom_local = iatom_local+1
              iatom_global = ni_global(iatom_local,wk_ibz%wdes1%comm_inb)

              ! use Bloch theorem to get the coupling derivative for the atoms in POSCAR
              ! may not always be wise, so it is not included in the above routines
              angle_q=dot_product(vqpt,centers(:,iatom_global)-primitive_positions(:,iatom_global))
              phase_q=exp(citpi*angle_q)
              angle_q=dot_product(vqpt,-centers(:,iatom_global))
              phase_p=exp(citpi*angle_q)*phase_q

              ! loop over directions
              do alpha =1,3

                 ! loop over bands
                 do nbk=1,nbands_k
                 do nbkp=1,nbands_kp

                    elph_coupling = cmplx(0.0_q,0.0_q,q)
                    if (lgaddd) then
                        PROFILING_START('elph_sandwich_gd')
                        call elph_gd(ddijdu_kkp(:,:,:,alpha,iatom_global,spinors), &
                                     wk_ibz,nproj, &
                                     dprojk(:,nbk,0,spinors_k), &
                                     dprojkp(:,nbkp,0,spinors_kp), &
                                     elph_coupling_d)
                        elph_coupling = elph_coupling + elph_coupling_d*phase_q
                        PROFILING_STOP('elph_sandwich_gd')
                    endif

                    if (lgaddr) then
                        PROFILING_START('elph_sandwich_gr')
                        call elph_gr(wk_ibz,&
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,0,spinors_k), &
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,0,spinors_kp), &
                                     nt,lmdim,ccdrij(:,:,iatom_global,alpha,spinors),&
                                     celenk(nbk),celenkp(nbkp), &
                                     elph_coupling_r)
                        elph_coupling = elph_coupling + elph_coupling_r*phase_p
                        PROFILING_STOP('elph_sandwich_gr')
                    endif

                    if (lgaddp) then
                        PROFILING_START('elph_sandwich_gp')
                        call elph_gp(wk_ibz,&
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,0,spinors_k),  &
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,0,spinors_kp), &
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,alpha,spinors_k),&
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,alpha,spinors_kp),&
                                     nt, &
                                     ccdij(:,:,iatom_global,spinors), &
                                     ccqij(:,:,iatom_global,spinors), &
                                     celenk(nbk),celenkp(nbkp), &
                                     elph_coupling_p)
                        elph_coupling = elph_coupling + elph_coupling_p*phase_p
                        PROFILING_STOP('elph_sandwich_gp')
                    endif

                    if (lgadds) then
                        PROFILING_START('elph_sandwich_gs')
                        call elph_gs(wk_ibz,&
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,0,spinors_k),  &
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,0,spinors_kp), &
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,alpha,spinors_k),&
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,alpha,spinors_kp),&
                                     nt, &
                                     ccqij(:,:,iatom_global,spinors), &
                                     elph_coupling_s)
                        elph_coupling = elph_coupling + elph_coupling_s*phase_p
                        PROFILING_STOP('elph_sandwich_gs')
                    endif

                    if (lgaddq) then
                        PROFILING_START('elph_sandwich_gq')
                        call elph_gq(wk_ibz,&
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,0,spinors_k),  &
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,0,spinors_kp), &
                                     dprojk(lmbase_at:lmbase_at+lmmax-1,nbk,alpha,spinors_k),&
                                     dprojkp(lmbase_at:lmbase_at+lmmax-1,nbkp,alpha,spinors_kp),&
                                     nt, &
                                     ccdij(:,:,iatom_global,spinors), &
                                     ccqij(:,:,iatom_global,spinors), &
                                     celenk(nbk), &
                                     elph_coupling_q)
                        elph_coupling = elph_coupling + elph_coupling_q*phase_p
                        PROFILING_STOP('elph_sandwich_gq')
                    endif

                    gcart_kkp(nbk,nbkp,alpha,iatom_global)=&
                    gcart_kkp(nbk,nbkp,alpha,iatom_global)+elph_coupling
                 enddo ! nbands k'
                 enddo ! nbands k

              enddo ! directions

              lmbase_at = lmbase_at + lmmax
           enddo ! loop over all atoms of type
           enddo ! loop over atom types
     enddo ! spinor k'
     enddo ! spinor k

  end subroutine elph_sandwich_ref

  !> Compute gV contribution to the electron-phonon matrix element
  !> see Eq 76 in PRB 100, 174304 (2019)
  subroutine elph_gv(ngrid_prim,cr,crp,dvdu_kkp,elph_coupling)
     integer, intent(in) :: ngrid_prim
     complex(q), intent(in) :: cr(:)
     complex(q), intent(in) :: crp(:)
     complex(q), intent(in) :: dvdu_kkp(:)
     complex(q), intent(out) :: elph_coupling

     elph_coupling=dot_product(cr(:),dvdu_kkp(:)*crp(:))/ngrid_prim

  end subroutine elph_gv

  !> Compute gD contribution to the electron-phonon matrix element
  !> see Eq 77 in PRB 100, 174304 (2019)
  subroutine elph_gd(ddijdu_kkp,wk_ibz,nproj,proj_k,proj_kp,elph_coupling)
     use tutor, only: vtutor
     use string, only: str
     complex(q), intent(in) :: ddijdu_kkp(:,:,:)
     type (wavefuna) :: wk_ibz
     integer, intent(in) :: nproj
     GDEF, intent(in) :: proj_k(:)
     GDEF, intent(in) :: proj_kp(:)
     complex(q), intent(out) :: elph_coupling
     ! Local variables
     integer lm,lmp
     integer ni,nt,ii,lmbase

     elph_coupling=cmplx(0.0_q,0.0_q,q)

     ni=0
     lmbase=0
     do nt = 1,wk_ibz%wdes1%ntyp ! loop for atoms in PRIMITIVE_POSITION
     do ii = 1, wk_ibz%wdes1%nityp(nt)
        ni=ni+1
        do lm=1,wk_ibz%wdes1%lmmax(nt)
        do lmp=1,wk_ibz%wdes1%lmmax(nt) ! not very clean because I dont use WKP_FBZ
           ! I use LMP, LM, but I dont know if this is correct
           elph_coupling=elph_coupling+&
           GCONJG(proj_k(lmbase+lm))*ddijdu_kkp(lmp,lm,ni)*(proj_kp(lmbase+lmp))
        enddo
        enddo
        lmbase = lmbase + wk_ibz%wdes1%lmmax(nt)
     enddo
     enddo

     if (lmbase /= nproj) then
        call vtutor%error('error in elph_make_gd '//str(lmbase)//' /= '//str(nproj))
     endif
  end subroutine elph_gd

  !> Compute gR contribution to the electron-phonon matrix element
  !> see Eq 79 in PRB 100, 174304 (2019)
  subroutine elph_gr(wk_ibz,proj_k,proj_kp,nt,lmdim,ccdrij,e_k,e_kp,elph_coupling)
     type (wavefuna) :: wk_ibz
     integer,intent(in) :: nt
     GDEF, intent(in) :: proj_k(1:wk_ibz%wdes1%lmmax(nt))
     GDEF, intent(in) :: proj_kp(1:wk_ibz%wdes1%lmmax(nt))
     integer, intent(in) :: lmdim
     complex(q), intent(in) :: ccdrij(1:lmdim,1:lmdim)
     complex(q), intent(in) :: e_k
     complex(q), intent(in) :: e_kp
     complex(q), intent(out) :: elph_coupling
     integer :: lm,lmp

     elph_coupling=cmplx(0.0_q,0.0_q,q)

     do lm=1,wk_ibz%wdes1%lmmax(nt)
     do lmp=1,wk_ibz%wdes1%lmmax(nt) ! not very clean because I dont use WKP_FBZ
        elph_coupling=elph_coupling+GCONJG(proj_k(lm))*ccdrij(lmp,lm)*(proj_kp(lmp))
     enddo
     enddo

     elph_coupling=elph_coupling*(e_kp-e_k)

  end subroutine elph_gr

  !> Compute gP contribution to the electron-phonon matrix element
  !> see Eq 78 in PRB 100, 174304 (2019)
  subroutine elph_gp(wk_ibz,proj_k,proj_kp,dproj_k,dproj_kp,nt,&
                     ccdij,ccqij,e_k,e_kp,elph_coupling)
     type (wavefuna) :: wk_ibz
     integer,intent(in) :: nt
     GDEF, intent(in) :: proj_k(:)
     GDEF, intent(in) :: proj_kp(:)
     GDEF, intent(in) :: dproj_k(:)
     GDEF, intent(in) :: dproj_kp(:)
     complex(q), intent(in) :: ccdij(:,:)
     complex(q), intent(in) :: ccqij(:,:)
     complex(q), intent(in) :: e_k
     complex(q), intent(in) :: e_kp
     complex(q),  INTENT(OUT) :: ELPH_COUPLING
     ! Local variables
     integer lm,lmp

     elph_coupling=cmplx(0.0_q,0.0_q,q)

     do lm=1,wk_ibz%wdes1%lmmax(nt)
     do lmp=1,wk_ibz%wdes1%lmmax(nt) ! not very clean because I dont use WKP_FBZ
        elph_coupling=elph_coupling+ &
                      GCONJG(dproj_k(lm))*(ccdij(lmp,lm)-e_kp*ccqij(lmp,lm))*(proj_kp(lmp))+ &
                      GCONJG(proj_k(lm))*(ccdij(lmp,lm)-e_k*ccqij(lmp,lm))*(dproj_kp(lmp))
     enddo
     enddo

  end subroutine elph_gp

  !> Compute gQ contribution to the electron-phonon matrix element
  subroutine elph_gq(wk_ibz,proj_k,proj_kp,dproj_k,dproj_kp,nt,&
                     ccdij,ccqij,e_k,elph_coupling)
     type (wavefuna) :: wk_ibz
     integer,intent(in) :: nt
     GDEF, intent(in) :: proj_k(:)
     GDEF, intent(in) :: proj_kp(:)
     GDEF, intent(in) :: dproj_k(:)
     GDEF, intent(in) :: dproj_kp(:)
     complex(q), intent(in) :: ccdij(:,:)
     complex(q), intent(in) :: ccqij(:,:)
     complex(q), intent(in) :: e_k
     complex(q),  INTENT(OUT) :: ELPH_COUPLING
     ! Local variables
     integer lm,lmp

     elph_coupling=cmplx(0.0_q,0.0_q,q)

     do lm=1,wk_ibz%wdes1%lmmax(nt)
     do lmp=1,wk_ibz%wdes1%lmmax(nt) ! not very clean because I dont use WKP_FBZ
        elph_coupling=elph_coupling+ &
                      GCONJG(dproj_k(lm))*(ccdij(lmp,lm)-e_k*ccqij(lmp,lm))*(proj_kp(lmp))+ &
                      GCONJG(proj_k(lm))*(ccdij(lmp,lm)-e_k*ccqij(lmp,lm))*(dproj_kp(lmp))
     enddo
     enddo

  end subroutine elph_gq

  !> Compute gS contribution to the electron-phonon matrix element
  subroutine elph_gs(wk_ibz,proj_k,proj_kp,dproj_k,dproj_kp,nt,&
                     ccqij,elph_coupling)
     type (wavefuna) :: wk_ibz
     integer,intent(in) :: nt
     GDEF, intent(in) :: proj_k(:)
     GDEF, intent(in) :: proj_kp(:)
     GDEF, intent(in) :: dproj_k(:)
     GDEF, intent(in) :: dproj_kp(:)
     complex(q), intent(in) :: ccqij(:,:)
     complex(q),  INTENT(OUT) :: ELPH_COUPLING
     ! Local variables
     integer lm,lmp

     elph_coupling=cmplx(0.0_q,0.0_q,q)

     do lm=1,wk_ibz%wdes1%lmmax(nt)
     do lmp=1,wk_ibz%wdes1%lmmax(nt) ! not very clean because I dont use WKP_FBZ
        elph_coupling=elph_coupling+ &
                      GCONJG(dproj_k(lm))*(ccqij(lmp,lm))*(proj_kp(lmp))+ &
                      GCONJG(proj_k(lm))*(ccqij(lmp,lm))*(dproj_kp(lmp))
     enddo
     enddo

  end subroutine elph_gs

  !> Compute the sandwich of the planewave part in reciprocal space
  !> using blas
  !>
  !> This approach requires a bit more FFTs but smaller zgemm contractions
  subroutine elph_sandwich_pw_recip(lgaddv, lgaddt, latt_cur, &
                                   elph_caches, nrspinors, npw, ntyp, nityp, nbands_k, nbands_kp, isp,&
                                   wk_fbz, kphase, wkp_fbz, vqpt, centers, primitive_positions, &
                                   dvdu_kkp, dmudu_kkp, gcart_kkp)
     use constant
     use string, only: str
     use elphon_base, only: elph_settings, elph_caches_t
     use wave_high, only: newwava, delwava
     use mlr_optic, only: inprod_w
     use tutor, only: vtutor
     !> Data structure containing information for the electron-phonon calculation
     integer,intent(in) :: nrspinors
     integer,intent(in) :: npw
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: nbands_k
     integer,intent(in) :: nbands_kp
     integer,intent(in) :: isp
     logical,intent(in) :: lgaddv
     logical,intent(in) :: lgaddt
     type (latt),intent(in) :: latt_cur
     type (elph_caches_t)  :: elph_caches
     type (wavefuna) :: wk_fbz
     type (wavefuna) :: wkp_fbz
     real(q),intent(in) :: kphase(3)
     real(q) :: vqpt(3)
     real(q) :: centers(:,:)
     real(q) :: primitive_positions(:,:)
     complex(q),intent(in) :: dvdu_kkp(:,:,:,:)
     complex(q),intent(in) :: dmudu_kkp(:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:) !nbands_k,nbands_kp,alpha,iatom_global
     ! Local variables
     type(wavefuna) :: hwk_fbz
     integer :: ii, nt, m
     integer :: nbk, nbkp
     integer :: iatom_global, alpha, idir_nabla
     integer :: spinors
     integer :: spinors_k, spinors_kp
     integer :: ipwstart
     integer :: npr_local
     real(q) :: angle_q
     complex(q) :: phase_q,phase_p
     complex(q) :: v_cr(wk_fbz%wdes1%grid%rl%np)
     complex(q) :: cwork(wk_fbz%wdes1%grid%mplwv)
     complex(q) :: nabla_cw(wk_fbz%wdes1%grid%mplwv)
     real(q) :: g1,g2,g3
     complex(q) :: phase(wk_fbz%wdes1%grid%rl%np)

     if (.not.(lgaddv.or.lgaddt)) return

     CALL newwava(hwk_fbz,wk_fbz%wdes1,wk_fbz%wdes1%nbands)
     npr_local = wk_fbz%wdes1%grid%rl%np

     ! loop over atoms
     iatom_global=0
     do nt = 1, ntyp
     do ii = 1, nityp(nt)
        iatom_global=iatom_global+1

        ! use Bloch theorem to get the coupling derivative for the atoms in POSCAR
        ! may not always be wise, so it is not included in the above routines
        angle_q=dot_product(vqpt,centers(:,iatom_global)-primitive_positions(:,iatom_global))
        phase_q=exp(citpi*angle_q)
        angle_q=dot_product(vqpt,-centers(:,iatom_global))
        phase_p=exp(citpi*angle_q)*phase_q

        phase=1.0_q
        if (any(abs(kphase(:))>1e-8)) then
           call get_phase_pw(wk_fbz%wdes1%grid,npr_local,kphase,phase)
        endif

        do alpha=1,3
            ! TODO: switch order of spinor index
            !do spinors_k=1,nrspinors
            !    elph_caches%v(:,spinors_k) = conjg(dvdu_kkp(:,alpha,iatom_global,spinors_k))
            !enddo
            do nbk=1,nbands_k
                ! apply H | psi >
                ! fft | psi_nk > to real space
                do spinors_k=1,nrspinors
                    ipwstart = (spinors_k-1)*wk_fbz%wdes1%ngvector
                    ! compute V(r) |psi>
                    if (lgaddv) then
                        ! FFT |psi> to real space -> <r|psi>
                        CALL FFTWAV(wk_fbz%wdes1%ngvector,wk_fbz%wdes1%nindpw,cwork,wk_fbz%CW(ipwstart+1,nbk),wk_fbz%wdes1%grid)
                        ! apply V(r) < r | psi_nk > in real space
                        !TODO: I would like to use this routine or something similar but the potential must be real
                        ! which is not the case here where it must be complex
                        !CALL VHAMIL(wk_fbz%wdes1,wk_fbz%wdes1%grid,elph_caches%v,cwork,v_cr(:,:))
                        v_cr = 0.0_q
                        do spinors_kp=1,nrspinors
                           spinors= isp+(spinors_kp-1)+2*(spinors_k-1)
                           v_cr(:) = &
                           v_cr(:) + &
                           conjg(dvdu_kkp(:,alpha,iatom_global,spinors))*cwork(:npr_local)*phase(:)
                        enddo
                        ! FFT V(r) < r | psi_nk > to reciprocal space -> < G |\hat{V} psi >
                        cwork = v_cr(:)
                        CALL FFTEXT(wkp_fbz%wdes1%ngvector,wkp_fbz%wdes1%nindpw,cwork,hwk_fbz%CW(ipwstart+1,nbk),wkp_fbz%wdes1%grid,.false.)
                    else
                        hwk_fbz%CW(:,nbk) = 0.0_q
                    endif
                    ! compute nabla.(mu(r) nabla) |psi>
                    if (lgaddt) then
                       do idir_nabla=1,3
                           ! nabla |psi>
                           call calc_nabla_cw(wk_fbz%wdes1, wk_fbz%wdes1%vkpt+kphase, latt_cur%b, -citpi, idir_nabla, wk_fbz%CW(ipwstart+1:ipwstart+wk_fbz%wdes1%ngvector,nbk), nabla_cw)
                           ! FFT | nabla psi > to real space -> < r | nabla psi >
                           CALL FFTWAV(wk_fbz%wdes1%ngvector,wk_fbz%wdes1%nindpw,cwork,nabla_cw,wk_fbz%wdes1%grid)
                           ! mu(r) . < r | nabla psi >
                           v_cr = 0.0_q
                           do spinors_kp=1,nrspinors
                               spinors= isp+(spinors_kp-1)+2*(spinors_k-1)
                               v_cr(:) = &
                               v_cr(:) + &
                               conjg(dmudu_kkp(:,alpha,iatom_global,spinors))*cwork(:npr_local)*phase(:)
                           enddo
                           ! FFT to reciprocal space -> < G | \hat{mu} nabla psi >
                           cwork = v_cr(:)
                           CALL FFTEXT(wkp_fbz%wdes1%ngvector,wkp_fbz%wdes1%nindpw,cwork,nabla_cw,wkp_fbz%wdes1%grid,.false.)
                           ! nabla | \hat{mu} nabla psi>  citpi is negative because we want the conjugate of H
                           call calc_nabla_cw(wkp_fbz%wdes1, wkp_fbz%wdes1%vkpt, latt_cur%b, -citpi, idir_nabla, nabla_cw, nabla_cw)
                           ! add to hwk_fbz%cw
                           do m=1,wkp_fbz%wdes1%ngvector
                               hwk_fbz%CW(ipwstart+m,nbk)=hwk_fbz%CW(ipwstart+m,nbk)-nabla_cw(m)
                           enddo
                       enddo !idir_nabla
                    endif
                enddo !spinors
            enddo !bands at k
            PROFILING_START('elph_sandwich_blas_gv_gemm')
            call zgemm('C','N',nbands_k,nbands_kp,wkp_fbz%wdes1%ngvector*nrspinors,&
                       phase_q/npw,&
                       hwk_fbz%CW,size(hwk_fbz%CW,1),&
                       wkp_fbz%CW,size(wkp_fbz%CW,1),&
                       (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
            PROFILING_STOP('elph_sandwich_blas_gv_gemm')
        enddo
     enddo
     enddo
     call delwava(hwk_fbz)

  end subroutine elph_sandwich_pw_recip

  !> compute the gradient of a single orbital
  subroutine calc_nabla_cw(wdes1,vkpt,b,c,idir,cw,nabla_cw)
      type(wavedes1),intent(in) :: wdes1 !< wave function descriptor for the current state
      integer,intent(in) :: idir !< direction in which to compute nabla
      real(q),intent(in) :: vkpt(3)
      real(q),intent(in) :: b(:,:) !< reciprocal space vectors
      complex(q),intent(in) :: c !< factor to multiply with
      complex(q),intent(in) :: cw(:) !< plane-wave coefficient
      complex(q),intent(out) :: nabla_cw(:) !< nabla plane-wave coefficient
      !local variables
      integer :: m
      real(q) :: g1,g2,g3
      real(q) :: g

      do m=1,wdes1%ngvector
         g1=wdes1%igx(m)+vkpt(1)
         g2=wdes1%igy(m)+vkpt(2)
         g3=wdes1%igz(m)+vkpt(3)
         g=(g1*b(idir,1)+&
            g2*b(idir,2)+&
            g3*b(idir,3))
         nabla_cw(m) = c*g*cw(m)
      enddo
  end subroutine calc_nabla_cw

  !> Compute the plane-wave part of the sanwhich of the electron-phonon potential
  !> between two states k and k'
  subroutine elph_sandwich_pw_blas(lgadd, &
                                   elph_caches, nrspinors, npw, ntyp, nityp, nbands_k, nbands_kp, isp,&
                                   wk_fbz, kphase, wkp_fbz, vqpt, centers, primitive_positions, &
                                   dvdu_kkp, gcart_kkp)
     use constant
     use string, only: str
     use elphon_base, only: elph_settings, elph_caches_t
     use tutor, only: vtutor
     !> Data structure containing information for the electron-phonon calculation
     integer,intent(in) :: nrspinors
     integer,intent(in) :: npw
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: nbands_k
     integer,intent(in) :: nbands_kp
     integer,intent(in) :: isp
     logical,intent(in) :: lgadd
     type (wavefuna) :: wk_fbz
     type (wavefuna) :: wkp_fbz
     real(q),intent(in) :: kphase(3)
     type (elph_caches_t)  :: elph_caches
     real(q) :: vqpt(3)
     real(q) :: centers(:,:)
     real(q) :: primitive_positions(:,:)
     complex(q),intent(in) :: dvdu_kkp(:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:) !nbands_k,nbands_kp,alpha,iatom_global
     ! Local variables
     integer :: ii, nt
     integer :: npr_local
     integer :: nbk
     integer :: iatom_global, alpha
     integer :: spinors
     integer :: spinors_k, spinors_kp
     real(q) :: angle_q
     complex(q) :: phase_q,phase_p
     complex(q) :: cr(wk_fbz%wdes1%grid%rl%np,nrspinors,wk_fbz%wdes1%nbands)
     complex(q) :: crp(wkp_fbz%wdes1%grid%rl%np,nrspinors,wkp_fbz%wdes1%nbands)

     if (.not.lgadd) return

     npr_local = wk_fbz%wdes1%grid%rl%np
     ! fft WFs
     call fft_wfs(wk_fbz%wdes1,nrspinors,wk_fbz%CW,cr)
     call fft_wfs(wkp_fbz%wdes1,nrspinors,wkp_fbz%CW,crp)
     if (any(abs(kphase(:))>1e-8)) then
        call multiply_phase_pw(wk_fbz%wdes1%grid,npr_local,kphase,cr)
     endif

     ! consistency check
     if (size(elph_caches%v_cr,1) /= npr_local) then
        call vtutor%bug('Wrong first 1 dimension in elph_caces%v_cr: '//str(size(elph_caches%v_cr,1))//&
                        ' expected: '//str(npr_local),__FILE__,__LINE__)
     endif
     if (size(elph_caches%v_cr,3) /= nbands_k) then
        call vtutor%bug('Wrong second 3 dimension in elph_caces%v_cr: '//str(size(elph_caches%v_cr,2))//&
                        ' expected: '//str(nbands_k),__FILE__,__LINE__)
     endif
     if (size(dvdu_kkp,1) /= npr_local) then
        call vtutor%bug('Wrong dimension 1 in dvdu_kkp: '//str(size(dvdu_kkp,1))//&
                        ' expected: '//str(npr_local),__FILE__,__LINE__)
     endif
     if (size(gcart_kkp,1)/=nbands_k) then
        call vtutor%bug('Wrong dimension 1 in gcart_kkp: '//str(size(gcart_kkp,1))//&
                        ' expected: '//str(nbands_k),__FILE__,__LINE__)
     endif
     iatom_global=0
     ! loop over atoms
     do nt = 1, ntyp
     do ii = 1, nityp(nt)
        iatom_global=iatom_global+1

        ! use Bloch theorem to get the coupling derivative for the atoms in POSCAR
        ! may not always be wise, so it is not included in the above routines
        angle_q=dot_product(vqpt,centers(:,iatom_global)-primitive_positions(:,iatom_global))
        phase_q=exp(citpi*angle_q)
        angle_q=dot_product(vqpt,-centers(:,iatom_global))
        phase_p=exp(citpi*angle_q)*phase_q

        do alpha =1,3
           ! Compute gV
           PROFILING_START('elph_sandwich_blas_gv_mul')
           elph_caches%v_cr=cmplx(0.0_q,0.0_q,q)
           do spinors_k=1,nrspinors
           do spinors_kp=1,nrspinors
              spinors= isp+(spinors_kp-1)+2*(spinors_k-1)
              do nbk=1,nbands_k
                 elph_caches%v_cr(:,spinors_kp,nbk) = &
                 elph_caches%v_cr(:,spinors_kp,nbk) + &
                 conjg(dvdu_kkp(:,alpha,iatom_global,spinors))*cr(:,spinors_k,nbk)
              enddo
           enddo
           enddo
           PROFILING_STOP('elph_sandwich_blas_gv_mul')
           PROFILING_START('elph_sandwich_blas_gv_gemm')
           call zgemm('C','N',nbands_k,nbands_kp,npr_local*nrspinors,&
                      phase_q/npw,&
                      elph_caches%v_cr,npr_local*nrspinors,&
                      crp,npr_local*nrspinors,&
                      (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
           PROFILING_STOP('elph_sandwich_blas_gv_gemm')
        enddo ! alpha
     enddo !ii
     enddo !n
  end subroutine elph_sandwich_pw_blas

  ! HM: This version is about 10 times faster than the above because if uses BLAS.
  !     20s vs 2s for 4x4x4 Diamond with 24 bands
  !     The speedup might be more significant for larger number of bands (ZPR calculations)
  subroutine elph_sandwich_blas(lgaddr,lgaddp,lgadds,lgaddq,lgaddd,&
                                elph_caches, lmdim, nrspinors, nproj, nproj_tot, ntyp, nityp, lmmax_all, nbands_k, nbands_kp, isp,&
                                wk_ibz, celenk, dprojk,&
                                celenkp,dprojkp,&
                                vqpt,centers,primitive_positions,&
                                ddijdu_kkp,ccdij,ccqij,ccdrij,gcart_kkp)
     use constant
     use wave, only: ni_global
     use elphon_base, only: elph_settings, elph_caches_t
     !> Data structure containing information for the electron-phonon calculation
     integer,intent(in) :: lmdim
     integer,intent(in) :: nrspinors
     integer,intent(in) :: nproj
     integer,intent(in) :: nproj_tot
     integer,intent(in) :: ntyp
     integer,intent(in) :: nityp(:)
     integer,intent(in) :: lmmax_all(:) !< lmmax for all the atoms
     integer,intent(in) :: nbands_k
     integer,intent(in) :: nbands_kp
     integer,intent(in) :: isp
     logical,intent(in) :: lgaddr,lgaddp,lgadds,lgaddq,lgaddd
     type (elph_caches_t)  :: elph_caches
     type (wavefuna) :: wk_ibz
     complex(q),intent(in) :: celenk(:)
     complex(q),intent(in) :: celenkp(:)
     GDEF,target,intent(in) :: dprojk(:,:,0:,:)
     GDEF,target,intent(in) :: dprojkp(:,:,0:,:)
     real(q) :: vqpt(3)
     real(q) :: centers(:,:)
     real(q) :: primitive_positions(:,:)
     complex(q),intent(in) :: ddijdu_kkp(:,:,:,:,:,:)
     complex(q),intent(in) :: ccdij(:,:,:,:)
     complex(q),intent(in) :: ccqij(:,:,:,:)
     complex(q),intent(in) :: ccdrij(:,:,:,:,:)
     complex(q),intent(inout) :: gcart_kkp(:,:,:,:)
     ! Local variables
     integer :: ii, nt
     integer :: iip, ntp
     integer :: nip
     integer :: nbk, nbkp
     integer :: lmbase, lmmax, lmmaxp
     integer :: iatom_local, iatom_global, alpha
     integer :: spinors
     integer :: spinors_k, spinors_kp
     integer :: lmbase_at
     real(q) :: e_k, e_kp
     real(q) :: angle_q
     complex(q) :: phase_q,phase_p
     complex(q) :: gcart_kkp_tmp(nbands_k,nbands_kp)
     complex(q) :: deqij(lmdim,lmdim)
     GDEF,allocatable :: dprojkp_mrg(:,:)
     GDEF,allocatable :: dprojk_mrg(:,:)

     if (wk_ibz%wdes1%comm_inb%ncpu/=1.and.lgaddd) then
         allocate(dprojkp_mrg(wk_ibz%wdes1%npro_tot,nbands_kp))
         allocate(dprojk_mrg(wk_ibz%wdes1%npro_tot,nbands_k))
     endif

     ! It might be interesting to save the spinor-resolved electron-phonon matrix elements
     do spinors_k=1,nrspinors
     do spinors_kp=1,nrspinors
        if (wk_ibz%wdes1%comm_inb%ncpu/=1.and.lgaddd) then
            do nbkp=1,nbands_kp
               call mrg_proj(wk_ibz%wdes1,dprojkp_mrg(:,nbkp),dprojkp(:,nbkp,0,spinors_kp))
            enddo
            do nbk=1,nbands_k
               call mrg_proj(wk_ibz%wdes1,dprojk_mrg(:,nbk),dprojk(:,nbk,0,spinors_k))
            enddo
        endif

        lmbase_at=1
        spinors= isp+(spinors_kp-1)+2*(spinors_k-1)
        iatom_local=0
        ! loop over atoms
        do nt = 1, wk_ibz%wdes1%ntyp
        do ii = 1, wk_ibz%wdes1%nityp(nt)
           lmmax = wk_ibz%wdes1%lmmax(nt)
           if (lmmax==0) cycle
           iatom_local = iatom_local+1
           iatom_global = ni_global(iatom_local,wk_ibz%wdes1%comm_inb)

           ! use Bloch theorem to get the coupling derivative for the atoms in POSCAR
           ! may not always be wise, so it is not included in the above routines
           angle_q=dot_product(vqpt,centers(:,iatom_global)-primitive_positions(:,iatom_global))
           phase_q=exp(citpi*angle_q)
           angle_q=dot_product(vqpt,-centers(:,iatom_global))
           phase_p=exp(citpi*angle_q)*phase_q

           do alpha =1,3
              ! Compute gR
              if (lgaddr) then
                 PROFILING_START('elph_sandwich_blas_gr')
                 !$omp parallel do private(nbkp)
                 do nbkp=1,nbands_kp
                    ! This 'T' is needed to agree with the non-BLAS implementation
                    call zgemv('T',lmmax,lmmax,(1.0_q,0.0_q),&
                               ccdrij(:,:,iatom_global,alpha,spinors),size(ccdrij,1),&
                               dprojkp(lmbase_at,nbkp,0,spinors_kp),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbkp),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,dprojk(lmbase_at,1,0,spinors_k),size(dprojk,1),&
                            elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            (0.0_q,0.0_q),gcart_kkp_tmp(:,:),nbands_k)
                 ! I don't know how to do this using BLAS
                 !$omp parallel do private(nbkp,nbk,e_kp,e_k)
                 do nbkp=1,nbands_kp
                    e_kp = real(celenkp(nbkp),q)
                    do nbk=1,nbands_k
                       e_k = real(celenk(nbk),q)
                       gcart_kkp(nbk,nbkp,alpha,iatom_global)=&
                       gcart_kkp(nbk,nbkp,alpha,iatom_global)+&
                       gcart_kkp_tmp(nbk,nbkp)*(e_kp-e_k)
                    end do
                 end do
                 !$omp end parallel do
                 PROFILING_STOP('elph_sandwich_blas_gr')
              endif

              ! Compute gP
              if (lgaddp) then
                 PROFILING_START('elph_sandwich_blas_gp')
                 !$omp parallel do private(nbkp,e_kp,deqij)
                 do nbkp=1,nbands_kp
                    e_kp = real(celenkp(nbkp),q)
                    deqij = &
                    ccdij(:,:,iatom_global,spinors)-e_kp*ccqij(:,:,iatom_global,spinors)
                    ! I am not sure about 'T' here, the matrix is symmetric
                    call zgemv('T',lmmax,lmmax,(1.0_q,0.0_q),deqij,size(deqij,1),&
                               dprojkp(lmbase_at,nbkp,0,spinors_kp),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbkp),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,dprojk(lmbase_at,1,alpha,spinors_k),size(dprojk,1),&
                            elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)

                 !$omp parallel do private(nbkp,e_k,deqij)
                 do nbk=1,nbands_k
                    e_k = real(celenk(nbk),q)
                    deqij = &
                    conjg(ccdij(:,:,iatom_global,spinors)-e_k*ccqij(:,:,iatom_global,spinors))
                    call zgemv('N',lmmax,lmmax,(1.0_q,0.0_q),deqij,size(deqij,1),&
                               dprojk(lmbase_at,nbk,0,spinors_k),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbk),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            dprojkp(lmbase_at,1,alpha,spinors_kp),size(dprojkp,1),&
                            (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
                 PROFILING_STOP('elph_sandwich_blas_gp')
              endif

              ! Compute gS
              if (lgadds) then
                 PROFILING_START('elph_sandwich_blas_gs')
                 !$omp parallel do private(nbkp,deqij,e_kp)
                 do nbkp=1,nbands_kp
                    e_kp = real(celenkp(nbkp),q)
                    deqij = &
                    ccqij(:,:,iatom_global,spinors)
                    ! I am not sure about 'T' here, the matrix is symmetric
                    call zgemv('T',lmmax,lmmax,(1.0_q,0.0_q),deqij,size(deqij,1),&
                               dprojkp(lmbase_at,nbkp,0,spinors_kp),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbkp),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,dprojk(lmbase_at,1,alpha,spinors_k),size(dprojk,1),&
                            elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)

                 !$omp parallel do private(nbk,e_k,deqij)
                 do nbk=1,nbands_k
                    e_k = real(celenk(nbk),q)
                    deqij = &
                    conjg(ccqij(:,:,iatom_global,spinors))
                    ! I am not sure about 'T' here, the matrix is symmetric
                    call zgemv('N',lmmax,lmmax,(1.0_q,0.0_q),deqij,size(deqij,1),&
                               dprojk(lmbase_at,nbk,0,spinors_k),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbk),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            dprojkp(lmbase_at,1,alpha,spinors_kp),size(dprojkp,1),&
                            (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
                 PROFILING_STOP('elph_sandwich_blas_gs')
              endif

              ! Compute gQ
              if (lgaddq) then
                 PROFILING_START('elph_sandwich_blas_gq')
                 !$omp parallel do private(nbk,deqij,e_k)
                 do nbk=1,nbands_k
                    e_k = real(celenk(nbk),q)
                    deqij = &
                    conjg(ccdij(:,:,iatom_global,spinors)-e_k*ccqij(:,:,iatom_global,spinors))
                    call zgemv('N',lmmax,lmmax,(1.0_q,0.0_q),deqij,size(deqij,1),&
                               dprojk(lmbase_at,nbk,0,spinors_k),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbk),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            dprojkp(lmbase_at,1,alpha,spinors_kp),size(dprojkp,1),&
                            (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)

                 !$omp parallel do private(nbk,e_k,deqij)
                 do nbk=1,nbands_k
                    e_k = real(celenk(nbk),q)
                    deqij = &
                    conjg(ccdij(:,:,iatom_global,spinors)-e_k*ccqij(:,:,iatom_global,spinors))
                    call zgemv('N',lmmax,lmmax,(1.0_q,0.0_q),deqij,size(deqij,1),&
                               dprojk(lmbase_at,nbk,alpha,spinors_k),1,&
                               (0.0_q,0.0_q),elph_caches%v_cproj(1,nbk),1)
                 enddo
                 !$omp end parallel do
                 call zgemm('C','N',nbands_k,nbands_kp,lmmax,&
                            phase_p,elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                            dprojkp(lmbase_at,1,0,spinors_kp),size(dprojkp,1),&
                            (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
                 PROFILING_STOP('elph_sandwich_blas_gq')
              endif

              ! Compute gD
              if (lgaddd) then
                 PROFILING_START('elph_sandwich_blas_gd')
                 if (wk_ibz%wdes1%comm_inb%ncpu==1) then
                     ! ncore==1
                     !$omp parallel do private(nbkp,lmbase,nip,lmmaxp)
                     do nbkp=1,nbands_kp
                        nip=0; lmbase=1
                        do ntp = 1,wk_ibz%wdes1%ntyp ! loop for atoms in PRIMITIVE_POSITION
                        do iip = 1,wk_ibz%wdes1%nityp(ntp)
                           lmmaxp = wk_ibz%wdes1%lmmax(ntp)
                           nip=nip+1
                           ! I am not sure about 'T' here, the matrix is symmetric
                           call zgemv('T',lmmaxp,lmmaxp,(1.0_q,0.0_q),&
                                      ddijdu_kkp(:,:,nip,alpha,iatom_global,spinors),size(ddijdu_kkp,1),&
                                      dprojkp(lmbase,nbkp,0,spinors_kp),1,&
                                      (0.0_q,0.0_q),elph_caches%v_cproj(lmbase,nbkp),1)
                           lmbase=lmbase+lmmaxp
                        enddo !iip
                        enddo !ntp
                     enddo ! nbands_kp
                     !$omp end parallel do
                     call zgemm('C','N',nbands_k,nbands_kp,nproj,&
                                phase_q,dprojk(:,:,0,spinors_k),size(dprojk,1),&
                                elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                                (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
                 else
                     !ncore/=1
                     !$omp parallel do private(nbkp,lmbase,nip,lmmaxp)
                     do nbkp=1,nbands_kp
                        nip=0; lmbase=1
                        do ntp = 1,ntyp ! loop for atoms in PRIMITIVE_POSITION
                        do iip = 1,nityp(ntp)
                           lmmaxp = lmmax_all(ntp)
                           nip = nip+1
                           call zgemv('T',lmmaxp,lmmaxp,(1.0_q,0.0_q),&
                                      ddijdu_kkp(:,:,nip,alpha,iatom_local,spinors),size(ddijdu_kkp,1),&
                                      dprojkp_mrg(lmbase,nbkp),1,&
                                      (0.0_q,0.0_q),elph_caches%v_cproj(lmbase,nbkp),1)
                           lmbase=lmbase+lmmaxp
                        enddo !iip
                        enddo !ntp
                     enddo ! nbands_kp
                     !$omp end parallel do
                     call zgemm('C','N',nbands_k,nbands_kp,nproj_tot,&
                                phase_q,dprojk_mrg(:,:),size(dprojk_mrg,1),&
                                elph_caches%v_cproj,size(elph_caches%v_cproj,1),&
                                (1.0_q,0.0_q),gcart_kkp(:,:,alpha,iatom_global),nbands_k)
                 endif
                 PROFILING_STOP('elph_sandwich_blas_gd')
              endif
           enddo ! alpha
           lmbase_at = lmbase_at + lmmax
        enddo !ii
        enddo !n
        !!$omp end parallel
     enddo ! spinor_kp
     enddo ! spinor_k
     if (wk_ibz%wdes1%comm_inb%ncpu/=1.and.lgaddd) then
        deallocate(dprojkp_mrg)
        deallocate(dprojk_mrg)
     endif
  end subroutine elph_sandwich_blas

end module elphon_mels_pw
