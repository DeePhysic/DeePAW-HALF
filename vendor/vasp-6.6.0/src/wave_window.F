#include "symbol.inc"

!>
!> @ brief Helper routiens to store and access a distibuted WF using one-sided comunication
!> by Henrique Miranda
!>
!> Currently we only support distribution over bands and k-points
!> 1. (local) intialize the object (window_init)
!> 2. (local) call store with the index bands and k-points that it owns (window_store)
!> 3. (global) call start (window_start) which propagates the info to all CPUs
!> 4. (one-sided) each cpu can acess any part of the WF using window_get
!> 5. wait for the data to arrive (can be used to overlap comunication and computation)
!> 6. (global) window_free cleans up the mess
!>

module wave_windower
    use prec
    use mpimy
    use string, only: str
    use tutor, only: vtutor
    implicit none

#ifndef MPI
    integer,parameter :: MPI_ADDRESS_KIND = 8
    integer,parameter :: MPI_DOUBLE_COMPLEX = 8
    integer,parameter :: MPI_STATUS_SIZE = 8
    integer,parameter :: MPI_LOCK_SHARED = 8
    integer,parameter :: MPI_INFO_NULL = 0
#endif

    !> Block of cache for the WFs
    type wave_cache_block
        integer :: nstates !< number of states in this block
        complex(q),allocatable :: wf(:,:) !< contains a set of bands
    end type wave_cache_block

    !> Datatype to cache wavefunctions
    type wave_cache
        integer :: nbands !< number of bands
        integer :: nkpoints !< number of k-points in the caching structure
        integer :: nspin !< number of spin components (same as ISPIN)
        integer :: nplw !< number of plane-waves
        integer :: band_start !< index of the first band
        integer :: band_stop !< index of the last band
        real(q) :: wf_cache_mb !< max size of the cache in MB
        integer :: nblocks !< current number of blocks
        integer :: nblocks_max !< maximum number of blocks
        integer :: nstates !< number of states
        integer(qi8) :: nstates_per_block !< numbero of states per block
        integer(qi8) :: nstates_max !< maximum number of WF blocks
        !> Mapping where to find the wavefunctions
        !> First dimension contains index of block and position within the block where to find the WF
        !> Other dimensions are iband, ikpt and ispin
        integer,allocatable :: bks_idx(:,:,:,:)
        !> Container for the wavefunctions
        type(wave_cache_block),pointer :: wave_blocks(:)
    end type wave_cache

    !> Datatype used to get wavefunctions from different CPUs
    !> using mpi one-sided comunication
    type wave_window
        ! global data
        integer :: nbands !< number of bands
        integer :: nkpoints !< number of k-points
        integer :: nspin !< number of spin components
        integer :: ncore !< number of cores sharing the information of one KS state (comm_inb%ncpu)
        integer :: nplw !< number of plane waves
        integer :: band_start !< index of the first band
        integer :: band_stop !< index of the last band
        !> Number of plane-waves in every other rank
        integer,allocatable :: nplw_all(:)
        !> MPI communicator between the different nodes holding the WF data
        type(communic),pointer :: comm
        integer,allocatable :: requests(:,:,:)
        ! local data
        integer :: icore !< index of the node for ncore (comm_inb%node_me)
        integer :: nkpoints_local !< number of k-points in local array
        integer :: nbands_local !< number of bands in local array
        integer :: nspin_local !< number of spins in local array
        integer,allocatable :: bks_idx(:,:,:,:) !< index of kpoint and bands (3,b,k,s) in local array
        integer :: win !< index of the window
        ! global data
        !> rank that contains the desired k-point and band (b,k,s)
        integer,allocatable :: bks_rank(:,:,:,:)
        integer :: bks_rank_size
        !> displacement in memory for the desired k-point and band
        integer(kind=MPI_ADDRESS_KIND),allocatable :: bks_disp(:,:,:,:)
        integer :: bks_disp_size
    end type wave_window

    ! window
    public :: wave_window_init
    public :: wave_window_store
    public :: wave_window_info
    public :: wave_window_start
    public :: wave_window_get
    public :: wave_window_free

    ! cache
    public :: wave_cache_init
    public :: wave_cache_setk
    public :: wave_cache_hask
    public :: wave_cache_info
    public :: wave_cache_getk
    public :: wave_cache_free

    ! combined
    public :: wave_window_cache_get
    public :: wave_window_cache_getk
    contains

    !> @brief Initialize the wave_window object
    subroutine wave_window_init(self,nbands,nkpoints,nspin,nplw,comm,comm_inb,band_start)
        type(wave_window) :: self
        type(communic),target :: comm !< global communicator where to create the window
        type(communic) :: comm_inb !< intra-band comunicator
        integer,intent(in) :: nbands,nkpoints,nspin,nplw
        integer,optional,intent(in) :: band_start

        PROFILING_START('wave_window_init')
        self%nbands = nbands
        self%nkpoints = nkpoints
        self%nspin = nspin
        self%nplw = nplw
        self%ncore = comm_inb%ncpu
        self%icore = comm_inb%node_me
        self%band_start = 1; if (present(band_start)) self%band_start=band_start
        self%band_stop = self%band_start+self%nbands-1
        !allocate the global position of WFs rank index
        allocate(self%bks_rank(self%ncore,self%nbands,self%nkpoints,self%nspin))
        self%bks_rank = 0
        self%bks_rank_size = product(shape(self%bks_rank))
        allocate(self%bks_disp(self%ncore,self%nbands,self%nkpoints,self%nspin))
        self%bks_disp = 0_MPI_ADDRESS_KIND
        ! TODO: this is a very ugly workaround to be able to sum integers of a different type
        ! we should implement a new interface in mpi.F
        self%bks_disp_size = product(shape(self%bks_disp))*storage_size(self%bks_disp)/storage_size(self%bks_rank)
        !local variables 
        self%nkpoints_local = -1
        self%nbands_local = -1
        self%comm => comm
        PROFILING_STOP('wave_window_init')
    end subroutine wave_window_init

    !> @brief Get the bks array from a wavedes structure
    !> returns allocated bks_idx
    subroutine get_bks_wavedes(wdes,bks_idx,nkpts,kpoints_index,band_start)
        use wave_struct_def
        type(wavedes) :: wdes
        integer,allocatable :: bks_idx(:,:,:,:)
        integer :: nkpts
        integer :: kpoints_index(:)
        integer,optional :: band_start
        !local variables
        integer :: is, ib, ik
        integer :: nkpoints, nspin
        integer :: my_band_start
        my_band_start = 1; if (present(band_start)) my_band_start = band_start
        nkpoints = nkpts
        nspin = wdes%ispin
        !if (size(kpoints_index)/=nkpoints) then
        !    call vtutor%bug('kpoints_index array has the wrong size got '//str(size(kpoints_index))//' but expected '//str(nkpoints),__FILE__,__LINE__)
        !endif
        allocate(bks_idx(3,wdes%nbands,nkpoints,nspin))
        bks_idx = 0
        do is=1,nspin
            do ik=1,nkpoints
                do ib=1,wdes%nbands
                    bks_idx(1,ib,ik,is) = wdes%nb_low+(ib-1)*wdes%nb_par+my_band_start-1
                    bks_idx(2,ib,ik,is) = kpoints_index(ik)
                    bks_idx(3,ib,ik,is) = is
                end do ! loop over bands
            end do ! loop over kpoints
        end do ! loop over spins
    end subroutine get_bks_wavedes

    !> @brief Create a window to local WF and expose it to the other CPUs
    subroutine wave_window_store(self,bks_idx,wf)
        !use iso_c_binding
        type(wave_window) :: self
        integer :: bks_idx(:,:,:,:) !> indexes of the desired WFs (nbands,nkpoints)
        complex(q) :: wf(:,:,:,:) !> container for the WF (nplw,nbands,nkpoints)
        !local variables
        !complex(q),pointer :: local_wf(:,:,:,:)
        !type(c_ptr) :: baseptr
        integer :: ib,ik,is
        integer :: ib_idx,ik_idx,is_idx
        integer :: ierror, doublesize, lowerbound
        integer(kind=MPI_ADDRESS_KIND) :: displacement
        integer(kind=MPI_ADDRESS_KIND) :: winsize
        integer :: nbands, nspin, nkpoints

#ifdef MPI
        PROFILING_START('wave_window_store')
        ! store k-points and bands dimensions
        self%nbands_local = size(bks_idx,2)
        self%nkpoints_local = size(bks_idx,3)
        self%nspin_local = size(bks_idx,4)
        nbands = self%nbands_local
        nkpoints = self%nkpoints_local
        nspin = self%nspin_local
        ! allocate indexes
        allocate(self%bks_idx(size(bks_idx,1),self%nbands_local,self%nkpoints_local,self%nspin_local))
        self%bks_idx = bks_idx
        ! compute size
        call mpi_type_get_extent( MPI_DOUBLE_COMPLEX, lowerbound, doublesize, ierror )
        winsize = 1
        winsize = winsize * doublesize * self%nbands_local * self%nkpoints_local * self%nspin_local * self%nplw
        ! create window to my WF
        call mpi_win_create(wf, winsize, doublesize, MPI_INFO_NULL, self%comm%mpi_comm, self%win, ierror)
        !call mpi_win_create_dynamic(MPI_INFO_NULL, self%comm%mpi_comm, self%win, ierror)
        !call mpi_win_attach(self%win, wf, winsize, ierror)
        !call mpi_win_allocate(winsize, doublesize, MPI_INFO_NULL, self%comm%mpi_comm, baseptr, self%win, ierror)
        !call c_f_pointer(baseptr,local_wf,shape(wf))
        !local_wf = wf

        if (ierror/=0) call vtutor%error('Failed to create MPI window')
        ! store indexes and displacements in bks_rank
        do is=1,self%nspin_local
            do ik=1,self%nkpoints_local
                do ib=1,self%nbands_local
                    ib_idx = self%bks_idx(1,ib,ik,is)
                    ik_idx = self%bks_idx(2,ib,ik,is)
                    is_idx = self%bks_idx(3,ib,ik,is)
                    self%bks_rank(self%icore,ib_idx,ik_idx,is_idx) = self%comm%node_me
                    displacement = ((ib-1)+(ik-1)*nbands+(is-1)*nkpoints*nbands)
                    displacement = displacement*self%nplw
                    if ((displacement+self%nplw)*doublesize > winsize) call vtutor%bug('Displacement too large',__FILE__,__LINE__)
                    if (self%bks_disp(self%icore,ib_idx,ik_idx,is_idx)/=0) call vtutor%bug('Displacement already set for ib: '//str(ib_idx)//' ik: '//str(ik_idx)//' is: '//str(is_idx),__FILE__,__LINE__)
                    self%bks_disp(self%icore,ib_idx,ik_idx,is_idx) = displacement
                end do ! loop over bands
            end do ! loop over kpoints
        end do ! loop over spins
        PROFILING_STOP('wave_window_store')
#endif
    end subroutine wave_window_store

    !> write information about the contents of the cache
    subroutine wave_window_info(self,tag,iu)
        type(wave_window),intent(in) :: self
        integer,intent(in) :: iu
        character(len=*) :: tag
        ! local variables
        integer :: ispin, ikpt, iband
        integer :: ib_idx, ik_idx, is_idx
        if (iu<=0) return
        write(iu,'(2a)') 'wave window ',tag
        write(iu,'(a)')  '======================'
        write(iu,'(a,i6)') 'ncore:      ', self%ncore
        write(iu,'(a,i6)') 'icore:      ', self%icore
        write(iu,'(a,i6)') 'band_start: ', self%band_start
        write(iu,'(a,i6)') 'band_stop:  ', self%band_stop
        write(iu,'(a,i6)') 'nbands:     ', self%nbands
        write(iu,'(a,i6)') 'nkpoints:   ', self%nkpoints
        write(iu,'(a,i6)') 'nspin:      ', self%nspin

        do ispin=1,self%nspin_local
          do ikpt=1,self%nkpoints_local
            do iband=1,self%nbands_local
                ib_idx = self%bks_idx(1,iband,ikpt,ispin)
                ik_idx = self%bks_idx(2,iband,ikpt,ispin)
                is_idx = self%bks_idx(3,iband,ikpt,ispin)
                write(iu,*) iband, ikpt, ispin, ib_idx, ik_idx, is_idx, self%bks_rank(self%icore,ib_idx,ik_idx,is_idx)
            enddo
          enddo
        enddo
        write(iu,*)
    end subroutine wave_window_info

    !> @brief Comunicate among the CPUs which nodes contain which information and where
    subroutine wave_window_start(self)
        type(wave_window) :: self
        ! local variables
        integer :: group
        integer :: ierror

        PROFILING_START('wave_window_start')
        allocate(self%nplw_all(self%comm%ncpu))
        self%nplw_all = 0
        self%nplw_all(self%comm%node_me) = self%nplw
        ! transmit nplw among all the CPUS
        CALLMPI( M_sum_i(self%comm, self%nplw_all, self%comm%ncpu) )
        ! transmit the indexes among all the CPUS 
        CALLMPI( M_sum_i(self%comm, self%bks_rank, self%bks_rank_size) )
        ! transmit the indexes among all the CPUS 
        CALLMPI( M_sum_i(self%comm, self%bks_disp, self%bks_disp_size) )
        !call mpi_comm_group(self%comm%mpi_comm, self%comm_group, ierror)
        !call mpi_group_incl(self%comm_group, 1, self%comm%node_me-1, group, ierror)
        !call mpi_win_post(group, 0, self%win, ierror)
        !call mpi_win_fence(0,self%win,ierror)
        !if (ierror/=0) call vtutor%error('Failed at mpi_win_fence')
        PROFILING_STOP('wave_window_start')
    end subroutine wave_window_start

    !> @brief Get the WFs for a certain pair of k-points and bands indexes
    subroutine wave_window_getbks(self,bks_idx,wf)
        type(wave_window) :: self
        integer :: bks_idx(:,:,:,:) !> indexes of the desired WFs (nbands,nkpoints)
        complex(q) :: wf(:,:,:,:) !> container for the WF (nplw,nbands,nkpoints,nspin)
#ifdef MPI
        !local variables
        integer :: ib, ik, is
        integer :: ib_idx, ik_idx, is_idx
        integer :: ierror
        integer :: target_rank
        integer(kind=MPI_ADDRESS_KIND) :: target_disp
        integer :: nbands, nkpoints, nspin
        ! check if the structure is initialized
        if (.not.allocated(self%bks_rank)) then
            call vtutor%bug('The wave_window structure was not initialized, '//&
                            'please call wave_window_init before calling this function',__FILE__,__LINE__)
        endif
        ! first check which rank has the information we need
        nbands = size(bks_idx,2)
        nkpoints = size(bks_idx,3)
        nspin = size(bks_idx,4)
        ! allocate array for requests
        !if (allocated(self%requests)) then
        !    call vtutor%bug('Some pending one-sided mpi comm requests were found',__FILE__,__LINE__)
        !endif
        !if (allocated(self%requests)) call vtutor%bug('There are pending requests',__FILE__,__LINE__)
        !allocate(self%requests(nbands,nkpoints,nspin))
        !self%requests = MPI_REQUEST_NULL
        do is=1,nspin
            do ik=1,nkpoints
                do ib=1,nbands
                    is_idx = bks_idx(3,ib,ik,is)
                    if (is_idx<1) cycle
                    ik_idx = bks_idx(2,ib,ik,is)
                    if (ik_idx<1) cycle
                    ib_idx = bks_idx(1,ib,ik,is)
                    if (ib_idx<1) cycle
                    target_rank = self%bks_rank(self%icore,ib_idx,ik_idx,is_idx)
                    target_disp = self%bks_disp(self%icore,ib_idx,ik_idx,is_idx)
                    !write(*,*) self%comm%node_me, target_rank, target_disp, ib, ik, is
                    !  get information for this band and k-point pair
                    if (target_rank<0) then
                        call vtutor%bug('Invalid rank '//str(target_rank)//' for ib: '//str(ib_idx)//&
                                          ' ik: '//str(ik_idx)//' is: '//str(is_idx),__FILE__,__LINE__)
                    endif
                    if (target_disp<0) then
                        call vtutor%bug('Invalid disp for ib: '//str(ib_idx)//&
                                          ' ik: '//str(ik_idx)//' is: '//str(is_idx),__FILE__,__LINE__)
                    endif
                    if (self%nplw_all(target_rank)>size(wf,1)) then
                        call vtutor%bug('Incompatible size of WF '//str(self%nplw_all(target_rank))//' '//&
                              str(size(wf,1)),__FILE__,__LINE__)
                    endif
                    call mpi_get(wf(:,ib,ik,is), self%nplw_all(target_rank), MPI_DOUBLE_COMPLEX, target_rank-1, target_disp,&
                                 self%nplw_all(target_rank), MPI_DOUBLE_COMPLEX, self%win, ierror )
                    if (ierror/=0) call vtutor%error('Error in mpi_Rget in wave_window')
                end do !loop over bands
            end do ! loop over kpoints
        end do ! loop over spins
#endif
    end subroutine wave_window_getbks

    !> @brief Get the WFs for a certain k-point and spin
    subroutine wave_window_getk(self,ispin,ikpt,wf,band_start,band_stop)
        type(wave_window) :: self
        integer,intent(in) :: ispin
        integer,intent(in) :: ikpt
        complex(q),intent(inout) :: wf(:,:) !> container for the WF (nplw,nbands)
        integer,optional,intent(in) :: band_start
        integer,optional,intent(in) :: band_stop
#ifdef MPI
        !local variables
        integer :: ierror
        integer :: target_rank
        integer :: origin_count
        integer(kind=MPI_ADDRESS_KIND) :: target_disp
        integer :: nbands
        integer :: i
        integer :: my_band_start
        integer :: my_band_stop

        my_band_start = 1;          if (present(band_start)) my_band_start = band_start
        my_band_stop = self%nbands; if (present(band_stop))  my_band_stop = band_stop
        PROFILING_START('wave_window_getk')
        ! check if the structure is initialized
        if (.not.allocated(self%bks_rank)) then
            call vtutor%bug('The wave_window structure was not initialized, '//&
                            'please call wave_window_init before calling this function',__FILE__,__LINE__)
        endif
        ! first check which rank has the information we need
        target_rank = self%bks_rank(self%icore,my_band_start,ikpt,ispin)
        target_disp = self%bks_disp(self%icore,my_band_start,ikpt,ispin)
        !  get information for this band and k-point pair
        if (target_rank<=0) then
            call vtutor%bug('Invalid rank '//str(target_rank)//&
                              ' ik: '//str(ikpt)//' is: '//str(ispin),__FILE__,__LINE__)
        endif
        if (target_disp<0) then
            call vtutor%bug('Invalid disp for ik: '//str(ikpt)//' is: '//str(ispin),__FILE__,__LINE__)
        endif
        if (self%nbands/=size(wf,2)) then
            call vtutor%bug('Incompatible number of bands '//str(self%nbands)//' '//str(size(wf,2)),__FILE__,__LINE__)
        endif
        !The dimensions of the CPTWFP array in the wavespin type of W_KBATCH
        !might in some cases be inconsistent with the WDES_DENSE_IBZ.
        !This is because the first k-points of W_KBATCH are the ones from KPOINTS
        !(this is done so that in the future we can handle hybrid functionals calculations)
        !and if incommensurate grids are specified then the WF sizes can be different.
        !In these cases, when using `WAVE_WINDOW` to copy data we need to be careful
        !that the memory alignment is preserved.
        !I implemented a workaround that transfers the information band-by-band.
        if (self%nplw_all(target_rank)/=size(wf,1)) then
            ! Here I will transfer band-by-band to avoid the allocation of a large array.
            ! To avoid multiple comunication requests a single request to an array with the
            ! (self%nplw_all(target_rank),self%nbands) dimension could be made
            PROFILING_START('wave_window_lock')
            call mpi_win_lock(MPI_LOCK_SHARED,target_rank-1,0,self%win,ierror)
            PROFILING_STOP('wave_window_lock')
            origin_count = min(size(wf,1),self%nplw_all(target_rank))
            do i=my_band_start,my_band_stop
                target_disp = self%bks_disp(self%icore,i,ikpt,ispin)
                call mpi_get(wf(:,i), origin_count, MPI_DOUBLE_COMPLEX, &
                             target_rank-1, target_disp,&
                             origin_count, MPI_DOUBLE_COMPLEX, self%win, ierror)
                if (ierror/=0) call vtutor%error('Error in mpi_get')
            enddo
            PROFILING_START('wave_window_unlock')
            call mpi_win_unlock(target_rank-1,self%win,ierror)
            PROFILING_STOP('wave_window_unlock')
        else if (self%nplw_all(target_rank)==size(wf,1)) then
            ! In this case the arrays are compatible so the transfer can be done on a single step
            PROFILING_START('wave_window_lock')
            !write(*,*) 'cpu', self%comm%node_me, 'requests', target_rank
            call mpi_win_lock(MPI_LOCK_SHARED,target_rank-1,0,self%win,ierror)
            PROFILING_STOP('wave_window_lock')
            nbands = my_band_stop-my_band_start+1
            call mpi_get(wf, self%nplw_all(target_rank)*nbands, MPI_DOUBLE_COMPLEX, &
                             target_rank-1, target_disp, self%nplw_all(target_rank)*nbands, &
                             MPI_DOUBLE_COMPLEX, self%win, ierror)
            PROFILING_START('wave_window_unlock')
            call mpi_win_unlock(target_rank-1,self%win,ierror)
            PROFILING_STOP('wave_window_unlock')
        endif
        PROFILING_STOP('wave_window_getk')
#endif
    end subroutine wave_window_getk

    !> @brief Get the WFs for a certain k-point and spin
    subroutine wave_window_get(self,ib,ikpt,ispin,wf)
        type(wave_window) :: self
        integer,intent(in) :: ib
        integer,intent(in) :: ikpt
        integer,intent(in) :: ispin
        complex(q),intent(inout) :: wf(:) !> container for the WF for one band (nplw)
#ifdef MPI
        !local variables
        integer :: ierror
        integer :: target_rank
        integer :: origin_count
        integer(kind=MPI_ADDRESS_KIND) :: target_disp

        PROFILING_START('wave_window_get')
        ! check if the structure is initialized
        if (.not.allocated(self%bks_rank)) then
            call vtutor%bug('The wave_window structure was not initialized, '//&
                            'please call wave_window_init before calling this function',__FILE__,__LINE__)
        endif
        ! first check which rank has the information we need
        target_rank = self%bks_rank(self%icore,ib,ikpt,ispin)
        target_disp = self%bks_disp(self%icore,ib,ikpt,ispin)
        !  get information for this band and k-point pair
        if (target_rank<=0) then
            call vtutor%bug('Invalid rank '//str(target_rank)//&
                              ' ib: '//str(ib)//' ik: '//str(ikpt)//' is: '//str(ispin),__FILE__,__LINE__)
        endif
        if (target_disp<0) then
            call vtutor%bug('Invalid disp for ib: '//str(ib)//' ik: '//str(ikpt)//' is: '//str(ispin),__FILE__,__LINE__)
        endif
        PROFILING_START('wave_window_lock')
        call mpi_win_lock(MPI_LOCK_SHARED,target_rank-1,0,self%win,ierror)
        PROFILING_STOP('wave_window_lock')
        origin_count = min(size(wf,1),self%nplw_all(target_rank))
        call mpi_get(wf(:), origin_count, MPI_DOUBLE_COMPLEX, &
                     target_rank-1, target_disp,&
                     origin_count, MPI_DOUBLE_COMPLEX, self%win, ierror)
        if (ierror/=0) call vtutor%error('Error in mpi_get')
        PROFILING_START('wave_window_unlock')
        call mpi_win_unlock(target_rank-1,self%win,ierror)
        PROFILING_STOP('wave_window_unlock')
        PROFILING_STOP('wave_window_get')
#endif
    end subroutine wave_window_get

    !> Get wavefunction for a particular state
    !> First check if the cache already contains the point, if it does not then
    !> get it at store it in the cache
    subroutine wave_window_cache_get(self,wcache,iband,ikpt,ispin,wf)
        type(wave_window) :: self
        type(wave_cache) :: wcache
        integer,intent(in) :: ispin
        integer,intent(in) :: ikpt
        integer,intent(in) :: iband
        complex(q),intent(inout) :: wf(:) !> container for the WF (nplw,nbands)
        ! local variables
        logical :: wcache_has
        ! check if this wavefunction is stored in cache
        wcache_has = all(wcache%bks_idx(:,iband,ikpt,ispin)/=0)
        if (wcache_has.and.wcache%wf_cache_mb>0) then
            call wave_cache_get(wcache,iband,ikpt,ispin,wf)
            return
        endif
        ! use the wave_window to get the wavefunction
        call wave_window_get(self,iband,ikpt,ispin,wf)
        ! if this wavefunction is not in cache then store it
        if ((.not.wcache_has).and.wcache%wf_cache_mb>0) then
           call wave_cache_set(wcache,iband,ikpt,ispin,wf)
        endif
    end subroutine wave_window_cache_get

    !> Get wavefunction at a particular k-point using cache
    !> First check if the cache already contains the point, if it does not then
    !> get it at store it in the cache
    subroutine wave_window_cache_getk(self,wcache,ispin,ikpt,wf,band_start,band_stop)
        type(wave_window) :: self
        type(wave_cache) :: wcache
        integer,intent(in) :: ispin
        integer,intent(in) :: ikpt
        complex(q),intent(inout) :: wf(:,:) !> container for the WF (nplw,nbands)
        integer,optional,intent(in) :: band_start
        integer,optional,intent(in) :: band_stop
        ! local variables
        logical :: wcache_hask
        ! check if this wavefunction is stored in cache
        wcache_hask = wave_cache_hask(wcache,ispin,ikpt,band_start,band_stop)
        if (wcache_hask.and.wcache%wf_cache_mb>0) then
            call wave_cache_getk(wcache,ispin,ikpt,wf,band_start,band_stop)
            return
        endif
        ! use the wave_window to get the wavefunction
        call wave_window_getk(self,ispin,ikpt,wf,band_start,band_stop)
        ! if this wavefunction is not in cache then store it
        if ((.not.wcache_hask).and.wcache%wf_cache_mb>0) then
           call wave_cache_setk(wcache,ispin,ikpt,wf,band_start,band_stop)
        endif
    end subroutine wave_window_cache_getk

    !> Perform all the local pensing requests
    subroutine wave_window_flush_local_all(self)
        type(wave_window),target :: self
        !local variables
        integer :: ierror
        PROFILING_START('wave_window_flush_local_all')
        CALLMPI( mpi_win_flush_local_all(self%win,ierror) )
        PROFILING_STOP('wave_window_flush_local_all')
    end subroutine wave_window_flush_local_all

    !> Perform all the remote pending requests
    subroutine wave_window_flush_all(self)
        type(wave_window),target :: self
        !local variables
        integer :: ierror
        PROFILING_START('wave_window_flush_all')
        CALLMPI( mpi_win_flush_all(self%win,ierror) )
        PROFILING_STOP('wave_window_flush_all')
    end subroutine wave_window_flush_all

    !> @brief  Wait for all the pending requests to complete
    subroutine wave_window_wait(self)
        type(wave_window),target :: self
        integer,allocatable :: statuses(:,:)
        !local variables
        integer :: nbands, nkpoints, nspin
        integer :: ierror
        PROFILING_START('wave_window_wait')
        if (.not.allocated(self%requests)) then
            call vtutor%bug('There are no pending requests',__FILE__,__LINE__)
        endif
        nbands = size(self%requests,1)
        nkpoints = size(self%requests,2)
        nspin = size(self%requests,3)
        ! allocate statuses
        allocate(statuses(MPI_STATUS_SIZE,nbands*nkpoints*nspin))
        CALLMPI( mpi_waitall(nbands*nkpoints*nspin,self%requests,statuses,ierror) )
        if (ierror/=0) then
            call vtutor%error('At least one one-sided mpi comm request could not be completed in wave_window')
        endif
        deallocate(statuses)
        deallocate(self%requests)
        PROFILING_STOP('wave_window_wait')
    end subroutine wave_window_wait

    !> @brief Free the wave_window object and the memory
    subroutine wave_window_free(self)
        type(wave_window) :: self
        !local variables
        integer :: ierror
        deallocate(self%bks_idx)
        deallocate(self%bks_rank)
        deallocate(self%bks_disp)
        deallocate(self%nplw_all)
        !call mpi_win_wait(self%win, ierror)
        !call mpi_win_fence(0,self%win,ierror)
        !if (ierror/=0) call vtutor%error('error win_fence')
        CALLMPI( mpi_win_free(self%win,ierror) )
        if (ierror/=0) call vtutor%error('Error calling win_free in wave_window')
    end subroutine wave_window_free

    !> Initialize the wavefunction structure
    subroutine wave_cache_init(self,nbands,nkpoints,nspin,nplw,wf_cache_mb,band_start)
        type(wave_cache) :: self
        integer,intent(in) :: nbands
        integer,intent(in) :: nkpoints
        integer,intent(in) :: nspin
        integer,intent(in) :: nplw
        real(q),intent(in) :: wf_cache_mb
        integer,optional,intent(in) :: band_start
        type(wave_cache_block) :: wave_block ! dummy
        self%nbands = nbands
        self%nkpoints = nkpoints
        self%nspin = nspin
        self%nplw = nplw
        self%wf_cache_mb = wf_cache_mb
        self%band_start = 1; if (present(band_start)) self%band_start=band_start
        self%band_stop = self%band_start+self%nbands-1
        self%nstates_max = nint(wf_cache_mb*1024*1024/(VMEM_SIZE(wave_block%wf)*self%nplw),qi8)
        self%nblocks = 0
        self%nstates = 0
        self%nblocks_max = 100 ! this can be changed
        self%nstates_per_block = 1_qi8*(self%nstates_max-1)/self%nblocks_max+1
        allocate(self%wave_blocks(self%nblocks_max))
        allocate(self%bks_idx(2,self%band_start:self%band_start+self%nbands,self%nkpoints,self%nspin))
        self%bks_idx = 0
        call wave_cache_new_block(self)
    end subroutine wave_cache_init

    !> Create a new block
    subroutine wave_cache_new_block(self)
        use ini, only: register_allocate
        type(wave_cache) :: self
        self%nblocks=self%nblocks+1
        if (self%nblocks>self%nblocks_max) call vtutor%bug('The number of blocks is too large',__FILE__,__LINE__)
        allocate(self%wave_blocks(self%nblocks)%wf(self%nplw,self%nstates_per_block))
        call register_allocate(VMEM_ARR_SIZE(self%wave_blocks(self%nblocks)%wf),"wf_cache")
        self%wave_blocks(self%nblocks)%nstates=0
    end subroutine wave_cache_new_block

    !> Set the wavefunctions at a particular k-point in the wavefunction datastructure
    subroutine wave_cache_set(self,iband,ikpoint,ispin,wf)
        use ini, only: register_allocate
        type(wave_cache) :: self
        integer,intent(in) :: iband
        integer,intent(in) :: ikpoint
        integer,intent(in) :: ispin
        complex(q) :: wf(:) !< plane-waves for one state (nplw)
        ! local variables
        integer :: nplw
        nplw = size(wf,1)
        if (nplw/=self%nplw) then
            call vtutor%bug('wrong size of wfs in wave_cache '//&
                            str(nplw)//' '//str(self%nplw),__FILE__,__LINE__)
        endif
        ! if the memory is full then don't store anymore
        if (self%nstates==self%nstates_max) call vtutor%error('wf cache is full') !return
        ! if this block is full, create a new one
        if (self%wave_blocks(self%nblocks)%nstates>=self%nstates_per_block) then
            call wave_cache_new_block(self)
        endif
        ! increment states and store
        self%nstates=self%nstates+1
        call wave_cache_block_set(self%wave_blocks(self%nblocks),wf)
        ! set entries on index array
        self%bks_idx(1,iband-self%band_start+1,ikpoint,ispin) = self%nblocks
        self%bks_idx(2,iband-self%band_start+1,ikpoint,ispin) = self%wave_blocks(self%nblocks)%nstates
    end subroutine wave_cache_set

    !> Set a new state in the cache block
    subroutine wave_cache_block_set(self,wf)
        type(wave_cache_block) :: self
        complex(q) :: wf(:) !< plane-waves for one state (nplw)
        self%nstates = self%nstates + 1
        self%wf(:,self%nstates) = wf
    end subroutine wave_cache_block_set

    !> Set the wavefunctions at a particular k-point in the wavefunction datastructure
    subroutine wave_cache_setk(self,ispin,ikpt,wf,band_start,band_stop)
        use ini, only: register_allocate
        type(wave_cache) :: self
        integer, intent(in) :: ikpt
        integer, intent(in) :: ispin
        complex(q) :: wf(:,:) !< wavefunction array (nplw,nbands)
        integer,optional :: band_start
        integer,optional :: band_stop
        ! local variables
        integer :: nbands, nplw
        integer :: iband
        integer :: my_band_start
        integer :: my_band_stop
        my_band_start=1;          if (present(band_start)) my_band_start = band_start
        my_band_stop=self%nbands; if (present(band_stop))  my_band_stop = band_stop
        nbands = size(wf,2)
        if (nbands/=self%nbands) then
            call vtutor%bug('wrong number of bands in wave_cache',__FILE__,__LINE__)
        endif
        nplw = size(wf,1)
        if (nplw/=self%nplw) then
            call vtutor%bug('wrong size of wfs in wave_cache '//&
                            str(nplw)//' '//str(self%nplw),__FILE__,__LINE__)
        endif
        ! set entries on index array
        do iband=my_band_start,my_band_stop
            call wave_cache_set(self,iband,ikpt,ispin,wf(:,iband))
        enddo
    end subroutine wave_cache_setk

    !> Check if we have all the bands at a particular k-point
    function wave_cache_hask(self,ispin,ikpt,band_start,band_stop)
        type(wave_cache) :: self
        integer, intent(in) :: ikpt
        integer, intent(in) :: ispin
        integer, optional, intent(in) :: band_start
        integer, optional, intent(in) :: band_stop
        ! local variables
        integer :: my_band_start
        integer :: my_band_stop
        logical :: wave_cache_hask
        my_band_start = self%band_start; if (present(band_start)) my_band_start = band_start
        my_band_stop = self%band_stop;   if (present(band_stop)) my_band_stop = band_stop
        wave_cache_hask = all(self%bks_idx(:,my_band_start:my_band_stop,ikpt,ispin)/=0)
    end function wave_cache_hask

    !> write information about the contents of the cache
    subroutine wave_cache_info(self,tag,iu)
        type(wave_cache),intent(in) :: self
        integer,intent(in) :: iu
        character(len=*) :: tag
        ! local variables
        integer :: ispin, ikpt, iband
        if (iu<=0) return
        write(iu,'(2a)') 'wave cache ',tag
        write(iu,'(a)')  '======================'
        write(iu,'(a,i6)') 'band_start: ', self%band_start
        write(iu,'(a,i6)') 'band_stop:  ', self%band_stop
        write(iu,'(a,i6)') 'nkpoints:   ', self%nkpoints
        write(iu,'(a,i6)') 'nspin:      ', self%nspin
        write(iu,'(a,i6)') 'nstates:    ', self%nstates
        write(iu,'(a,i10,a)') 'memory max: ', self%wf_cache_mb*1000, ' (kBytes)'
        write(iu,'(a,i10,a)') 'memory used:', int(VMEM_ARR_SIZE(self%wave_blocks(1)%wf)*self%nblocks/1000), ' (kBytes)'

        do ispin=1,self%nspin
          do ikpt=1,self%nkpoints
            do iband=1,self%nbands
                if (all(self%bks_idx(:,iband,ikpt,ispin)==0)) cycle
                write(iu,'(3I3,A3,2I3)') iband,ikpt,ispin,' - ',self%bks_idx(:,iband,ikpt,ispin)
            enddo
          enddo
        enddo
        write(iu,*)
    end subroutine wave_cache_info

    !> Get a copy of a state
    subroutine wave_cache_get(self,iband,ikpt,ispin,wf)
        type(wave_cache) :: self
        integer, intent(in) :: ikpt
        integer, intent(in) :: ispin
        complex(q), intent(out) :: wf(:)
        ! local varibles
        integer :: nplw
        integer :: iband, ipos, iblock
        nplw = size(wf,1)
        if (nplw/=self%nplw) then
            call vtutor%bug('wrong size of wfs in wave_cache '//&
                            str(nplw)//' '//str(self%nplw),__FILE__,__LINE__)
        endif
        iblock=self%bks_idx(1,iband,ikpt,ispin)
        ipos=self%bks_idx(2,iband,ikpt,ispin)
        if (iblock==0) call vtutor%bug('The element with index ib: '//str(iband)//' ik: '//str(ikpt)//' is: '//str(ispin)//&
                                       ' was not found in wave_cache',__FILE__,__LINE__)
        if (ipos==0) call vtutor%bug('The element with index ib: '//str(iband)//' ik: '//str(ikpt)//' is: '//str(ispin)//&
                                     ' was not found in wave_cache',__FILE__,__LINE__)
        wf(:) = self%wave_blocks(iblock)%wf(:,ipos)
    end subroutine wave_cache_get

    !> Get pointer to a state
    subroutine wave_cache_get_ptr(self,iband,ikpt,ispin,wf)
        type(wave_cache),target :: self
        integer, intent(in) :: ikpt
        integer, intent(in) :: ispin
        complex(q), pointer, intent(out) :: wf(:)
        ! local varibles
        integer :: iband, ipos, iblock
        iblock=self%bks_idx(1,iband,ikpt,ispin)
        ipos=self%bks_idx(2,iband,ikpt,ispin)
        if (iblock==0) call vtutor%bug('The element with index ib: '//str(iband)//' ik: '//str(ikpt)//' is: '//str(ispin)//&
                                       ' was not found in wave_cache',__FILE__,__LINE__)
        if (ipos==0) call vtutor%bug('The element with index ib: '//str(iband)//' ik: '//str(ikpt)//' is: '//str(ispin)//&
                                     ' was not found in wave_cache',__FILE__,__LINE__)
        wf => self%wave_blocks(iblock)%wf(:,ipos)
    end subroutine wave_cache_get_ptr

    !> Get all the states at a particular k-point
    subroutine wave_cache_getk(self,ispin,ikpt,wf,band_start,band_stop)
        type(wave_cache) :: self
        integer, intent(in) :: ikpt
        integer, intent(in) :: ispin
        complex(q), intent(out) :: wf(:,:)
        integer,optional,intent(in) :: band_start
        integer,optional,intent(in) :: band_stop
        ! local varibles
        integer :: nplw
        integer :: iband
        integer :: my_band_start
        integer :: my_band_stop
        nplw = size(wf,1)
        if (nplw/=self%nplw) then
            call vtutor%bug('wrong size of wfs in wave_cache '//&
                            str(nplw)//' '//str(self%nplw),__FILE__,__LINE__)
        endif
        my_band_start=self%band_start
        if (present(band_start)) my_band_start=band_start
        my_band_stop=self%band_start+self%nbands-1
        if (present(band_stop)) my_band_stop=band_stop
        do iband=my_band_start,my_band_stop
            call wave_cache_get(self,iband,ikpt,ispin,wf(:,iband))
        enddo
    end subroutine wave_cache_getk

    !> Free wave cache datastructure
    subroutine wave_cache_free(self)
        use ini, only: deregister_allocate
        type(wave_cache) :: self
        ! local variables
        integer :: i
        do i=1,self%nblocks
            call deregister_allocate(VMEM_ARR_SIZE(self%wave_blocks(i)%wf),"wf_cache")
            deallocate(self%wave_blocks(i)%wf)
        enddo
        deallocate(self%wave_blocks)
        deallocate(self%bks_idx)
    end subroutine wave_cache_free
end module wave_windower
