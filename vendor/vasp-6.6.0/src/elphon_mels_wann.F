#include "symbol.inc"

module elphon_mels_wann
  use prec
  use base, only: symmetry
  use poscar_struct_def, only: latt
  use wannier_interpolation, only: wannier_interpolator
  use wave_interpolate, only: wave_mapper
  use wave, ONLY: wavespin, wavedes
  use mkpoints_struct_def, ONLY: kpoints_struct
  use mgrid_struct_def, ONLY: grid_3d
  use mpimy, ONLY: communic
  use wave_rotate, only: wave_rotator, fbz2ibz_stars

  !> Shim containing the wannier interpolator plus additional structures
  !> needed in the context of electron-phonon coupling calculations
  type elph_mels_wann_type
      type(wannier_interpolator) :: wann_interp
      type(wavedes) :: wdes_dense_ibz
      type(wavespin) :: w_dense_ibz
      complex(q),allocatable :: celtot(:,:,:) !< wavefunction coefficients in the wannier basis
      GDEF, allocatable :: ham_kpath(:,:,:,:) !< eigenvectors of the wannier hamiltonian on a regular grid
      type(grid_3d) :: grid
      type(grid_3d) :: grid_orig
      type(wave_rotator) :: waverot !< dataype to rotate wfs
      type(fbz2ibz_stars) :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
      type(symmetry) :: symm
      type(latt) :: latt_ini ! set in init
      type(latt) :: latt_cur ! set in init
      integer :: ngdim_orig !< number of g-vectors in the original wf
      integer :: ispin !< number of spins
      type(wave_mapper) :: wave_map !< mapping of the k-points in the k-point batch to the IBZ
      integer,allocatable :: ikp_fbz_list(:) !< indexed of the k' points present locally
      ! phonons
      real(q),allocatable :: phonon_freqs_ibz(:,:) !< store the phonon frequencies in the IBZ
  end type elph_mels_wann_type

  contains

  !> initialize elph_mels_wann structure
  subroutine elph_mels_wann_init(self,wdes,grid,latt_ini,latt_cur,symm,ngdim_orig)
      type(elph_mels_wann_type),intent(inout) :: self
      type(wavedes),intent(in) :: wdes
      type(grid_3d),intent(in) :: grid
      type(latt),intent(in) :: latt_ini
      type(latt),intent(in) :: latt_cur
      type(symmetry),intent(in) :: symm
      integer,intent(in) :: ngdim_orig

      self%ispin = wdes%ispin
      self%grid = grid
      self%grid_orig = grid
      self%latt_ini = latt_ini
      self%latt_cur = latt_cur
      self%symm = symm
      self%wdes_dense_ibz = wdes
      self%ngdim_orig = ngdim_orig
  end subroutine elph_mels_wann_init

  !> get wave function descriptor and wavefunction container in the dense IBZ
  !> The waverot and fbz2ibz structure are also initialized here.
  subroutine elph_mels_wann_get_w_dense(self,wdes_dense_ibz,w_dense_ibz,kpoints)
      use tutor, only: vtutor
      use string, only: str
      use wave_rotate, only: wave_rotator_init, fbz2ibz_stars_init, fbz2ibz_stars_dummy_init, kpoint_star, kstar_from_fbz2ibz
      use mkpoints, only: kpoints_is_regular, LineMode, vkpt_kindex
      use elphon_common, only: elphon_get_waverot_fbz2ibz, elphon_get_wdes_w
      type(elph_mels_wann_type),intent(inout) :: self
      type(wavedes),intent(inout) :: wdes_dense_ibz
      type(wavespin),intent(inout) :: w_dense_ibz
      type(kpoints_struct),intent(inout) :: kpoints

      ! intializes the waverot and fbz2ibz structures also update fbz2ibz mapping in kpoints
      call elphon_get_waverot_fbz2ibz(self%waverot,self%fbz2ibz,self%symm,self%latt_cur,self%latt_ini,&
                                      self%grid_orig,self%ngdim_orig,wdes_dense_ibz,kpoints)

      ! Initialize wavefunction descriptor for the dense IBZ
      ! BEWARE: this changes grid!
      call elphon_get_wdes_w(wdes_dense_ibz,w_dense_ibz,kpoints,self%fbz2ibz,self%grid,self%latt_cur,self%latt_ini,self%ngdim_orig)

  end subroutine elph_mels_wann_get_w_dense

  !> distribute the k-points over the processes, compute the hamiltonian
  !> in the wannier basis and then compute the eigenstates
  subroutine elph_mels_wann_compute_wfs(self,kpoints_inter,symm,latt_ini,latt_cur,t_info,info,io,comm_kbatch)
      use base, only: info_struct, in_struct
      use wannier_interpolation, only: wannier_interpolate_ham
      use wave_interpolate, only: wave_mapper_init
      use poscar, only: type_info
      type(elph_mels_wann_type),intent(inout) :: self
      type(kpoints_struct),intent(in) :: kpoints_inter !< k-points in the IBZ where the wannier hamiltonian is computed
      type(symmetry),intent(in) :: symm !< symmetry operations
      type(latt),intent(in) :: latt_ini !< initial lattice
      type(latt),intent(in) :: latt_cur !< current lattice
      type(type_info),intent(in) :: t_info !< timing information
      type(info_struct),intent(in) :: info !< information structure
      type(in_struct),intent(in) :: io !< IO structure
      type(communic),intent(in) :: comm_kbatch !< communicator for the k-point batch
      ! local variables
      integer :: ikpt,ikcoarse,ikdense
      integer :: isp
      integer :: nkpts_batch
      real(q) :: kpoints_batch(3,kpoints_inter%nkpts)
      integer :: kpoints_index(kpoints_inter%nkpts)

      ! distribute IBZ k-points over the processes
      nkpts_batch = 0
      do ikpt = 1, kpoints_inter%nkpts
#ifdef MPI
          if (mod(ikpt,comm_kbatch%ncpu)/=comm_kbatch%node_me-1) cycle
#endif
          nkpts_batch = nkpts_batch + 1
          kpoints_batch(:,nkpts_batch) = kpoints_inter%vkpt(:,ikpt)
          kpoints_index(nkpts_batch) = ikpt
      end do
      ! create a dummy wave_map (need to be carefull when we introduce paralelism)
      call wave_mapper_init(self%wave_map,0,nkpts_batch,kpoints_batch,kpoints_index)

      ! compute wannier eigenstates on the distributed k-point mesh
      allocate(self%celtot(self%wann_interp%nw,nkpts_batch,self%wdes_dense_ibz%ispin))
      call wannier_interpolate_ham(self%wann_interp,self%celtot,&
                                   self%ham_kpath,self%wave_map%nkpts_batch,self%wave_map%kpoints_batch,&
                                   symm,latt_ini,latt_cur,t_info,info,io)

      ! gather all the celtot to the root process
      self%w_dense_ibz%celtot = 0.0_q
      do isp=1,self%wdes_dense_ibz%ispin
         do ikpt=1,self%wave_map%nkpts_batch
            ikcoarse = self%wave_map%nkpts_orig+ikpt
            ikdense = self%wave_map%kpoints_index(ikpt)
            self%w_dense_ibz%celtot(:,ikdense,isp) = self%celtot(:,ikcoarse,isp)
         end do
      end do
      CALLMPI( M_sum_z( comm_kbatch, self%w_dense_ibz%celtot, size(self%w_dense_ibz%celtot)))
  end subroutine elph_mels_wann_compute_wfs

  !> compute the velocity operator in the wannier basis
  subroutine elph_mels_wann_get_velocity(self,hamiltonian,w,grid,gridc,gridus,c_to_us,irdmax,lmdim,t_info,&
                                          nonlr_s,nonl_s,p,sv,cqij,cdij,latt_ini,latt_cur,symm,velocity,info,io,comm_kbatch)
        use base, only: info_struct, in_struct
        use hamil_struct_def, only: ham_handle
        use wannier_interpolation, only: wannier_interpol_setup_cder, wannier_interpolate_cder
        use mkpoints_struct_def, only: kpoints_struct
        use mgrid_struct_def, only: grid_3d
        use nonlr_struct_def, only: nonlr_struct
        use nonl_struct_def, only: nonl_struct
        use mgrid, only: transit
        use poscar, only: type_info
        use pseudo, only: potcar
        use tutor, only: vtutor
        use string, only: str
        type(elph_mels_wann_type),intent(inout) :: self
        type(ham_handle),intent(in) :: hamiltonian !< hamiltonian structure
        type(wavespin),intent(in) :: w
        type(grid_3d),intent(in) :: grid !< grid structure
        type(grid_3d),intent(in) :: gridc !< coarse grid structure
        type(grid_3d),intent(in) :: gridus !< ultra-coarse grid structure
        type(transit),intent(in) :: c_to_us !< coarse to ultra-coarse grid mapping
        integer,intent(in) :: irdmax !< maximum number of irreducible directions
        integer,intent(in) :: lmdim !< dimension of the local mode
        type(type_info),intent(in) :: t_info !< timing information
        type(nonlr_struct),intent(in) :: nonlr_s
        type(nonl_struct),intent(in) :: nonl_s
        type(potcar),intent(in) :: p(:) !< number of irreducible directions
        RGRID,intent(in)   :: sv(DIMREAL(grid%mplwv),w%wdes%ncdij)
        OVERLAP,intent(in) :: cqij(:,:,:,:)
        OVERLAP,intent(in) :: cdij(:,:,:,:)
        type(latt),intent(in) :: latt_ini !< initial lattice
        type(latt),intent(in) :: latt_cur !< current lattice
        type(symmetry),intent(in) :: symm !< symmetry operations
        real(q),allocatable,intent(out) :: velocity(:,:,:,:) !< velocity operator in the wannier basis
        type(info_struct),intent(in) :: info !< information structure
        type(in_struct),intent(in) :: io !< IO structure
        type(communic),intent(in) :: comm_kbatch !< communicator for the k-point batch
        ! local variables
        type(kpoints_struct) :: kpoints_dummy
        GDEFS,allocatable :: chams_kpath(:,:,:,:,:) !< hamiltonian in the k-path
        real(q),allocatable :: velocity_kbatch(:,:,:,:) !< velocity operator in the k-point batch

        call wannier_interpol_setup_cder(self%wann_interp,hamiltonian,w,kpoints_dummy,grid,gridc,gridus,c_to_us,irdmax,lmdim,t_info,&
                                         nonlr_s,nonl_s,p,sv,cqij,cdij,latt_ini,latt_cur,symm,info,io)
        call wannier_interpolate_cder(self%wann_interp,self%ham_kpath,self%celtot,velocity_kbatch,&
                                      chams_kpath,self%wave_map%nkpts_batch,self%wave_map%kpoints_batch,symm,latt_ini,latt_cur,t_info,info,io)

        ! allocate the velocity operator
        allocate(velocity(self%wdes_dense_ibz%nb_tot,self%wdes_dense_ibz%nkpts,self%ispin,3)) ! 3 directions, nkpts_orig, spins, 3 components
        velocity = 0.0_q
        ! collect the velocity operator
        do idir=1,3
          do isp=1,self%ispin
             do ik=1,self%wave_map%nkpts_batch
                ikcoarse = self%wave_map%nkpts_orig+ik
                ikdense = self%wave_map%kpoints_index(ik)
                velocity(:,ikdense,isp,idir) = velocity_kbatch(:,ikcoarse,isp,idir)     ! in ev ang
             enddo !loop kpoints
          enddo !loop spins
       enddo ! loop directions
       CALLMPI( M_sum_d( comm_kbatch, velocity, size(velocity)))

  end subroutine elph_mels_wann_get_velocity

  !> compute g matrix for the wannier hamiltonian
  subroutine elph_mels_wann_compute_g_ibz(self,elph_set,kpoints,ikfbz,ikpfbz,isp,vqpt,gcart_kkp)
      use elphon_base, only: elph_settings
      type(elph_mels_wann_type),intent(inout) :: self
      type(elph_settings),intent(in) :: elph_set
      type(kpoints_struct),intent(in) :: kpoints
      integer,intent(in) :: ikfbz !< index of the k-point in the FBZ
      integer,intent(in) :: ikpfbz !< index of the k' point in the FBZ
      integer,intent(in) :: isp !< spin index
      real(q),intent(out) :: vqpt(3) !< phonon q-point in cartesian coordinates
      complex(q),intent(inout) :: gcart_kkp(3) !< g-vector in cartesian coordinates
      ! local variables
      integer :: ikpibz,ikibz,iqfbz
      integer :: irotp
      real(q),allocatable :: vk(:),vkp(:)

      ! get indexes of kpoints
      vk = kpoints%vkpt(:,ikibz)
      vkp = kpoints%vkpt(:,ikpibz)
      vqpt = vk-vkp

      ! early exit if computing g is not required
      if (elph_set%selfen_g_skip) then
        gcart_kkp = 1.0_q
        return
      endif

      ! TODO: implement this
      gcart_kkp = 1.0_q
  end subroutine elph_mels_wann_compute_g_ibz

  !> compute g matrix for the wannier hamiltonian
  subroutine elph_mels_wann_compute_g(self,elph_set,kpoints,ikfbz,ikpfbz,iqfbz,isp,vqpt,gcart_kkp)
      use tutor, only: vtutor
      use string, only: str
      use elphon_base, only: elph_settings
      use mkpoints, only: kindex_vkpt, kpoint_integers, vkpt_kindex, kpoints_fbz_to_ibz
      use wave_rotate, only: rotate_vkpt_irot
      type(elph_mels_wann_type),intent(inout) :: self
      type(elph_settings),intent(in) :: elph_set
      type(kpoints_struct),intent(in) :: kpoints
      integer,intent(in) :: ikfbz !< index of the k-point in the FBZ
      integer,intent(in) :: ikpfbz !< index of the k' point in the FBZ
      integer,intent(in) :: isp !< spin index
      real(q),intent(out) :: vqpt(3) !< phonon q-point in cartesian coordinates
      complex(q),intent(out) :: gcart_kkp(:,:,:,:) !< electron-phonon matrix elements in the cartesian basis
      ! local variables
      integer :: ikpibz,ikibz,iqfbz_local
      real(q),allocatable :: vk(:),vkp(:)

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

      ! TODO: implement this
      gcart_kkp = 1.0_q
  end subroutine elph_mels_wann_compute_g

  !> free all datastructures from elph_mels_wann structure
  subroutine elph_mels_wann_free(self)
      type(elph_mels_wann_type),intent(inout) :: self
  end subroutine elph_mels_wann_free

end module elphon_mels_wann
