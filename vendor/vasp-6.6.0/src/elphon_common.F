#include "symbol.inc"

module elphon_common
  use prec
  implicit none
  contains
  !>
  !> @brief Private subroutine to initialize communicators.
  !>
  !> This routine splits the COMM_IN communicator in the following order:
  !> 1. among NKBATCH k-point batches
  !> 2. among KPAR for vasp k-point paralelism (always dummy=1 because we don't use this VASP level)
  !> 3. among NPAR = NCPU/NCORE for band or plane-wave paralelism
  !> This is used by subroutine ELECTRON_PHONON.
  !>
  subroutine elphon_comm_init(wdes_kbatch, nkbatch, ncore, comm_kbatch, comm_in, io)
   use base, only : in_struct
   use mpimy, only : communic
#ifdef MPI
   use mpimy, only : m_divide, m_initc
#endif
   use wave_struct_def, only: wavedes
   use tutor, only: vtutor
   use string, only: str

   type(wavedes),intent(inout) :: wdes_kbatch
   integer,intent(inout) :: nkbatch !< number of k-point batches
   integer,intent(in) :: ncore !< number of cores treating the same orbital
   type(communic),intent(inout) :: comm_in !< input comunicator
   type(communic),intent(inout) :: comm_kbatch !< communicator for each batch of kpoints
   type(in_struct),intent(in) :: io
#ifdef MPI
   ! local variables
   type(communic),target :: comm  !< new wf communicator
   type(communic),target :: comm_kinter !< between k-points communicator
   type(communic),target :: comm_kin   !< in k-point communicator
   type(communic),target :: comm_inter !< between-band communicator
   type(communic),target :: comm_inb   !< in-band communicator

   integer :: kpar_dummy, npar

   ! check if the settings are valid
   if (mod(comm_in%ncpu,ncore)/=0) then
     call vtutor%error('Total number of MPI ranks '//str(comm_in%ncpu)//' is not divisible by NCORE='//str(ncore))
   endif
   ! check if the total number of cores is correct
   if (nkbatch<1) then
       nkbatch = comm_in%ncpu/ncore
   else
       if (mod(comm_in%ncpu,nkbatch)/=0) then
         call vtutor%error('Total number of MPI ranks '//str(comm_in%ncpu)//' is not divisible by ELPH_KPAR='//str(nkbatch))
       endif
   endif
   call m_divide( comm_in, nkbatch, comm_kbatch, comm, .true. )
   call m_initc( comm_in )
   call m_initc( comm )
   call m_initc( comm_kbatch )

   ! create new k-points comunicator
   kpar_dummy=1
   call m_divide( comm, kpar_dummy, comm_kinter, comm_kin, .false.)
   call m_initc( comm )
   call m_initc( comm_kinter )
   call m_initc( comm_kin )

   ! create inter and intra band comunicators
   npar=max(comm_kin%ncpu/ncore,1)
   call m_divide( comm_kin, npar, comm_inter, comm_inb, .false.)
   call m_initc( comm_kin )
   call m_initc( comm_inter )
   call m_initc( comm_inb )

   ! set comunicators
   wdes_kbatch%comm = comm
   wdes_kbatch%comm_kinter = comm_kinter
   wdes_kbatch%comm_kin = comm_kin
   wdes_kbatch%comm_inter = comm_inter
   wdes_kbatch%comm_inb = comm_inb

   ! set dimensions
   ! TODO: can change the number of nbands to be divisible by NCORE see logic in main.F
   wdes_kbatch%nb_low=comm_inter%node_me
   wdes_kbatch%nb_par=npar

   if (io%iu0>0) then
      write(io%iu0,'(A,I4,A,I4,A)') " elph_comm: each k-point on ", comm%ncpu, " cores, ",  comm_kbatch%ncpu, " groups"
      write(io%iu0,'(A,I4,A,I4,A)') " elph_comm: one band on ", comm_inb%ncpu, " cores, ", comm_inter%ncpu, " groups"
   endif
#endif
  end subroutine elphon_comm_init

  !> @brief Private subroutine to prepare elph_pot.
  !>   Set the dimensions and allocates arrays in ELPH_POT
  !>   Reads potentials arrays from phelel_params.hdf5
  !>   Prepare atoms centered grids
  !>   Prepare supercell lattice vectors for grids point and supercell atoms
  !>   for obtaining the phases in the mimimun images calculations
  subroutine elphon_prepare_elph_pot(elph_pot, elph_set, latt_cur, grid, p, t_info, io)
   use base, only: in_struct
   use lattice, only: latt
   use pseudo, only: potcar
   use poscar, only: type_info
   use mpimy, only: communic
   use mgrid_struct_def, only: grid_3d
#ifdef VASP_HDF5
   !use vhdf5, only: ih5outfileid
   !use elphon_potential, only: elphon_write_ccdrij
#endif
   use elphon_potential, only: elphon_pot_setup, elphon_pot_setup_rij, &
   elphon_quadrupoles_set, elphon_pot_setup_lr, elphon_pot_write_lr, &
   elphon_pot_write_structural, elphon_pot, &
   lattice_check, centers_check, elphon_write_pot, &
   elphon_pot_write_dielectric_born_quad
   use elphon_base, only: elph_settings

   type(elphon_pot),intent(inout) :: elph_pot
   type(elph_settings),intent(in) :: elph_set
   type(latt),intent(in) :: latt_cur
   type(grid_3d),target,intent(in) :: grid
   type(potcar),intent(in) :: p(:)
   type(type_info),intent(in):: t_info
   type(in_struct),intent(in) :: io

   PROFILING_START('elph_prepare_elph_pot')

   ! skip the computation of the g matrix elemetns
   elph_pot%skip = elph_set%selfen_g_skip

   !> check that lattice vector are the same than in POSCAR, and
   !> primitive_positions as well, including order, up to a lattice vector.
   call lattice_check(elph_pot%lattice_pc,elph_pot%natoms_pc,elph_pot%primitive_positions,latt_cur,t_info)

   !> enforce that primitive_positions are exactly those of t_info%posion (= POSCAR folded to [0,1[)
   elph_pot%primitive_positions(:,:)=t_info%posion(:,:)

   !>  This subroutine defines several arrays used for the interpolation of the
   !>  potential using the minimal image convention:
   !>      ELPH_POT%supercell_atoms, ELPH_POT%center, ELPH_POT%grid_point,
   !>      ELPH_POT%sp_lattice_points, ELPH_POT%sp_lattice_atoms
   call elphon_pot_setup(elph_pot,t_info,p)

   !> add a checking routine like phonon_ifc_setup_check

   !> compute ELPH_POT%ccdrij(LMP,LM,NI,ISP,DIR) = R_{LM LMP}(NI)
   call elphon_pot_setup_rij(elph_pot,p,t_info)

   !> check that center are POSCAR atoms, upto a lattice vector, including order
   call centers_check(elph_pot%centers,t_info)

   ! not checked
   if (elph_set%qval1/=0.or.elph_set%qval2/=0) then
      call elphon_quadrupoles_set(elph_pot%quadrupoles,elph_set%qval1,elph_set%qval2)
   endif

   !> lr is for long range, not linear response
   !> ELPH_SET%ENCUTLR = energy cut off, which give g cut off,
   !> for long range calculation

   if (elph_set%elph_lr==1) then
      !>compute vr= ELPH_POT%dvdu_longrange(self%ngridpoints,3,self%natoms_pc) on the supercell grid
      !> compute ddijdu_longrange(lm,lmp,iatom_sc,idir,iatom)
      ! xxx I beleive lm,lmp -> lmp,lm
      call elphon_pot_setup_lr(elph_pot,grid,elph_set%encutlr,elph_set%lr_sft,io%iu0)

      !> ELPH_POT%dvdu = ELPH_POT%dvdu - ELPH_POT%dvdu_longrange
      !> ELPH_POT%ddijdu = ELPH_POT%ddijdu - ELPH_POT%ddijdu_longrange
      ! ASK_HM where is has_lr_active put to .FALSE. ?
      !call elphon_pot_lr_set(ELPH_POT,.TRUE.)
      !deallocate(elph_pot%dvdu_longrange)
   ENDIF

   elph_pot%rcut = elph_set%rcut
#ifdef VASP_HDF5
   !call elphon_write_ccdrij(ELPH_POT,ih5outfileid)
#endif
   if (elph_set%writepot) then
        call elphon_write_pot(ELPH_POT,IO%IU0)
   endif
   call elphon_pot_write_lr(ELPH_POT,IO%IU0)
   call elphon_pot_write_lr(ELPH_POT,IO%IU6)
   call elphon_pot_write_dielectric_born_quad(ELPH_POT,IO%IU6)
   call elphon_pot_write_structural(ELPH_POT,IO%IU6)
   PROFILING_STOP('elph_prepare_elph_pot')
  end subroutine elphon_prepare_elph_pot

  !> @brief Private subroutine to prepare elph_ifc.
  subroutine elphon_prepare_phonon(ELPH_IFC, ELPH_SET, LATT_CUR, GRID, T_INFO, IO, COMM_GLOBAL)
   USE base, ONLY: in_struct
   USE lattice, ONLY: latt
   USE poscar, ONLY: type_info
   USE elphon_base, ONLY: elph_settings
   USE mgrid_struct_def, ONLY: grid_3d
   USE mpimy, only: communic
   USE elphon_potential, ONLY: phonon_ifc, &
         phonon_ifc_init, phonon_ifc_from_hdf5, phonon_ifc_setup, &
         phonon_ifc_lr_set, phonon_ifc_setup_lr, phonon_ifc_asr, phonon_ifc_write_lr, &
         lattice_check, centers_check, phonon_ifc_setup_check

   type(phonon_ifc),intent(INOUT) :: ELPH_IFC
   type(elph_settings),intent(in) :: ELPH_SET
   type(latt),intent(in) :: LATT_CUR
   type(grid_3d),target,intent(in) :: GRID
   type(type_info),intent(in):: T_INFO
   type(in_struct),intent(in) :: IO
   type(communic),intent(in) :: COMM_GLOBAL

   PROFILING_START('elph_prepare_phonon')

   !> check that lattice vector are then same than in POSCAR, and
   !> primitive_positions as well, including order, up to a lattice vector.
   !> Because the same is in elphon_prepare_elph_pot, the phonon and
   !> electron-phonon calculations must have the same primitive cell.
   call lattice_check(elph_ifc%lattice_pc,elph_ifc%natoms_pc,elph_ifc%primitive_positions,latt_cur,t_info)

   !> enforce that primitive_positions are exactly those of t_info%posion = POSCAR folded to [0,1[
   elph_ifc%primitive_positions(:,:)=t_info%posion(:,:)

   !>documentation can be found in elphon_pot_setup
   call phonon_ifc_setup(elph_ifc,t_info)

   !> check shortest vectors are the same than in phelel
   call phonon_ifc_setup_check(elph_ifc)

   !> check that center are POSCAR atoms, upto a lattice vector, including order
   call centers_check(elph_ifc%centers,t_info)

   !> accoutis sum rule, and also take the average with the transposed
   ! ASK_HM why dont we do that after LR part ?
   if (elph_set%ifc_asr>0) call phonon_ifc_asr(elph_ifc,elph_set%ifc_asr)

   if (elph_set%ifc_lr==1) then
      ! compute long range force constants_
      call phonon_ifc_setup_lr(elph_ifc,grid,elph_set%encutlr,COMM_GLOBAL)
      !  force_constants = force_constants - force_constants_lr
      call phonon_ifc_lr_set(ELPH_IFC,.TRUE.)
   end if
   call phonon_ifc_write_lr(elph_ifc,io%iu0,io%iu6)
   PROFILING_STOP('elph_prepare_phonon')
  end subroutine elphon_prepare_phonon

  !> @brief Private subroutine to prepare kpoints_inter.
  subroutine elphon_get_kpoints_inter(kpoints_inter, kgrid, kpoints, elph_set, latt_cur, t_info, &
                                      wdes, comm_kbatch, isym, io)
   use base, only: in_struct
   use constant, only : pi
   use string, only: str
   USE tutor, ONLY: vtutor
   USE mpimy, ONLY: communic
   USE lattice, ONLY: latt
   USE poscar, ONLY: type_info
   USE mkpoints, ONLY: kpoints_struct, copy_kpoints, rd_kpoints, kpoints_is_regular
   USE m_tetra_hash, ONLY: tetirr_hash
   USE wave_struct_def, ONLY: wavedes
   USE elphon_base, ONLY: elph_settings
   USE elphon_kgrid, ONLY: elph_kgrid

   type(kpoints_struct),intent(inout) :: kpoints_inter
   type(elph_kgrid),intent(inout) :: kgrid
   type(kpoints_struct),intent(in) :: kpoints
   type(elph_settings),intent(in) :: elph_set
   type(latt),intent(in) :: latt_cur
   type(type_info),intent(in):: t_info
   type(wavedes),target,intent(in) :: wdes
   type(communic),target,intent(in) :: comm_kbatch
   integer,intent(in) :: isym
   type(in_struct),intent(in) :: io
   ! local variables
   integer :: ik
   character(13) :: kpoints_filename
   logical :: kpoints_exists
   integer,allocatable :: idtet(:,:)

   PROFILING_START('elph_get_kpoints_inter')

   ! Every node reads the number of k-points
   call copy_kpoints(kpoints,kpoints_inter)

   !KPOINTS_INTER%LTET=ANY(ELPH_SET%SELFEN_DELTA<=0).OR.KPOINTS_INTER%ISMEAR<=-4.OR.ELPH_SET%TRANSPORT
   kpoints_inter%ltet=.false.
   kpoints_inter%lkmap=.true.
   kpoints_inter%ismear=elph_set%ismear

   kpoints_inter%spacing=elph_set%kspacing

   !try reading from KPOINTS_ELPH file
   kpoints_filename = "KPOINTS_ELPH"
   INQUIRE(file=KPOINTS_FILENAME, exist=kpoints_exists)
   !if KPOINTS_ELPH does not exist, try KPOINTS_DENSE and emit deprecation warning
   if (.not.kpoints_exists) then
        INQUIRE(file="KPOINTS_DENSE", exist=kpoints_exists)
        if (kpoints_exists) then
             call vtutor%alert('The use of KPOINTS_DENSE is deprecated, use a KPOINTS_ELPH file or ELPH_KSPACING in the INCAR instead.')
             kpoints_filename = "KPOINTS_DENSE"
        endif
   endif
#ifdef oldsym
   call rd_kpoints(kpoints_inter, latt_cur, &
        isym>=0.and..not. wdes%lsorbit.and..not. wdes%lspiral, &
        isym<0, -1, -1, file_name=kpoints_filename)
#else
   call rd_kpoints(kpoints_inter, latt_cur, &
        isym>=0.and..not. wdes%lnoncollinear, isym<0, &
        -1, -1, file_name=kpoints_filename)
#endif
   if (kpoints_is_regular(kpoints_inter)) then
      call tetirr_hash(kpoints_inter,latt_cur%b,kpoints_inter%ntet,idtet,-1)
      if (associated(kpoints_inter%idtet)) deallocate(kpoints_inter%idtet)
      allocate(kpoints_inter%idtet(0:4,kpoints_inter%ntet))
      kpoints_inter%idtet = idtet(0:4,:kpoints_inter%ntet)
      deallocate(idtet)
   endif

   ! write information about k-points
   if (io%iu6>=0) then
        write(io%iu6,*) kpoints_filename
        write(io%iu6,*) '--------------'
        write(io%iu6,*)
        write(io%iu6,*) 'kpoints in the ibz:',kpoints_inter%nkpts
        write(io%iu6,*) 'kpoints in the fbz:',kpoints_inter%isnf(1,1)*&
                                              kpoints_inter%isnf(2,2)*&
                                              kpoints_inter%isnf(3,3)
        write(io%iu6,*)
        write(io%iu6,*) 'Smith normal form'
        write(io%iu6,*) 'D = P R Q'
        write(io%iu6,*)

        write(io%iu6,*) 'D:'
        write(io%iu6,*) kpoints_inter%isnf(1,1)
        write(io%iu6,*) kpoints_inter%isnf(2,2)
        write(io%iu6,*) kpoints_inter%isnf(3,3)
        write(io%iu6,*)

        write(io%iu6,*) 'P:'
        write(io%iu6,*) kpoints_inter%ilsnf(1,:)
        write(io%iu6,*) kpoints_inter%ilsnf(2,:)
        write(io%iu6,*) kpoints_inter%ilsnf(3,:)
        write(io%iu6,*)
   endif

   ! safety check for more k-points than CPUs
   if (kpoints_inter%nkpts<comm_kbatch%ncpu) then
      call vtutor%alert('Number of k-points in the irreducible brillouin zone: '//str(kpoints_inter%nkpts)//&
                        ' is lower than ELPH_KPAR='//str(comm_kbatch%ncpu)//'. '//&
                        'which is inefficient. We recommend manually setting ELPH_KPAR to a value smaller than '//&
                        str(kpoints_inter%nkpts)//' and relaunch the calculation.')
   end if
   PROFILING_STOP('elph_get_kpoints_inter')

  end subroutine elphon_get_kpoints_inter

  subroutine elphon_prepare_ikpt(elph_set, kpoints_inter, latt_cur, wdes_dense_ibz, io)
   use base, only: in_struct
   use string, only: str
   use tutor, only: vtutor
   use lattice, only: latt
   use mkpoints, only: kpoints_struct, ExplicitList, LineMode, kpoint_integers, kint_kindex, kpoints_along_path, kpoints_along_path_write
   use wave_struct_def, only: wavedes
   use elphon_base, only: elph_settings
#ifdef VASP_HDF5
   use vhdf5, only: HID_T, GRP_RESULTS, ih5outfileid, kpoints_along_path_hdf5write, vh5_group_open_or_create, vh5_error
#endif
   type(elph_settings),intent(inout) :: elph_set
   type(kpoints_struct),intent(in) :: kpoints_inter
   type(latt),intent(in) :: latt_cur
   type(wavedes),intent(in) :: wdes_dense_ibz
   type(in_struct),intent(in) :: io
   ! local variables
   integer :: i, ikfbz, ikpt, kint(3)
   logical :: kpoints_line_exists
#ifdef VASP_HDF5
   integer(HID_T) :: groupid
#endif
   real(q),allocatable :: kpoints_dists(:)
   real(q),allocatable :: kpoints_vkpt(:,:)

   ! Check if ikpts is allocated, if not use all the points
   if (.not.allocated(elph_set%selfen_ikpt)) then
      allocate(elph_set%selfen_ikpt(wdes_dense_ibz%nkpts))
      elph_set%selfen_ikpt = [(i, i = 1, wdes_dense_ibz%nkpts, 1)]
   endif

   ! Check if ikpts values are within range
   if (any(elph_set%selfen_ikpt<=0.or.elph_set%selfen_ikpt>wdes_dense_ibz%nkpts)) then
      call vtutor%error('Invalid IKPT indexes in ELPH_SELFEN_IKPT: '//str(elph_set%selfen_ikpt)//&
                        '. Should be between 1 and '//str(wdes_dense_ibz%nkpts))
   ENDIF

   ! Now if ELPH_SELFEN_KPTS is set then we use these k-point coordinates to set the SELFEN_IKPT array
   if (allocated(elph_set%selfen_kpts)) then
      deallocate(elph_set%selfen_ikpt)
      allocate(elph_set%selfen_ikpt(size(elph_set%selfen_kpts,2)))
      do i=1,size(elph_set%selfen_kpts,2)
        if (kpoints_inter%mode==LineMode.or.kpoints_inter%mode==ExplicitList) then
            do ikpt=1,kpoints_inter%nkpts
                if (sqrt(sum((kpoints_inter%vkpt(:,ikpt)-elph_set%selfen_kpts(:,i))**2))<1e-3_q) then
                    elph_set%selfen_ikpt(i)=ikpt
                endif
                if (elph_set%selfen_ikpt(i)==-1) then
                    call vtutor%error('Could not find '//str(elph_set%selfen_kpts(:,i))//' in the k-point list')
                endif
            enddo
        else
            ! Since the coordinates are input by the user, use a lower threshold for comparing points
            call kpoint_integers(kpoints_inter,latt_cur,elph_set%selfen_kpts(:,i),kint,eps=1e-3_q)
            if (all(kint==-1)) call vtutor%error( 'Could not find '//str(elph_set%selfen_kpts(:,i))//' in the regular grid')
            ikfbz = kint_kindex(kpoints_inter,kint)
            elph_set%selfen_ikpt(i) = kpoints_inter%fbz2ibz(1,ikfbz)
        endif
      enddo
      ! some output
      if (io%iu6>0) then
          write(io%iu6,*) ' SELFEN_KPTS               IKPT  MATCHED TO'
          do i=1,size(elph_set%selfen_kpts,2)
            ikpt = elph_set%selfen_ikpt(i)
            write(io%iu6,'(3f8.4,a,i6,3f8.4)') elph_set%selfen_kpts(:,i), '->', ikpt, kpoints_inter%vkpt(:,ikpt)
          enddo
      endif
   endif

   inquire(file="KPOINTS_LINE", exist=kpoints_line_exists)
   if (kpoints_line_exists) then
       deallocate(elph_set%selfen_ikpt)
       call kpoints_along_path(kpoints_inter,latt_cur,kpoints_dists,kpoints_vkpt,elph_set%selfen_ikpt)
       call kpoints_along_path_write(elph_set%selfen_ikpt,kpoints_dists,kpoints_vkpt,io%iu6)
#ifdef VASP_HDF5
       VH5_CHECK(vh5_group_open_or_create(ih5outfileid,trim(GRP_RESULTS),groupid))
       call kpoints_along_path_hdf5write(groupid,elph_set%selfen_ikpt,kpoints_dists,kpoints_vkpt,group='electron_phonon')
#endif
   endif
  end subroutine elphon_prepare_ikpt

  !> @brief Private subroutine to initialize the self-energy accumulators
  subroutine elphon_prepare_selfenergy_accumulators(elph_accumulators, latt_cur, kpoints_inter, elph_set, wdes_dense_ibz, &
                                                    w_dense_ibz, phonon_freqs_ibz, t_info, info, chempot, bandgap, io)
   use base, only: in_struct,info_struct
   use tutor, only: vtutor
   use string, only: str
   use mkpoints_struct_def, only: kpoints_struct
   use lattice, only: latt
   use poscar, only: type_info
   use wave_struct_def, only: wavedes, wavespin
   use vhdf5
   use elphon_accumulators_high, only: elph_accumulators_selfen
   use elphon_base, only: elph_settings, wdes_get_global_band_index, elph_driver_el, elph_driver_ph
   use electron_transport, only: chemical_potential
   use bandgap_struct, only: bandgap_info
   use elphon_accumulators, only: elph_accumulator_selfen_print
   use elphon_accumulators_ph, only: elph_accumulator_selfen_print_ph

   type(elph_accumulators_selfen),intent(inout) :: elph_accumulators
   type(latt),intent(in) :: latt_cur
   type(kpoints_struct),intent(in) :: kpoints_inter
   type(elph_settings),intent(inout) :: elph_set
   type(wavedes),intent(in) :: wdes_dense_ibz
   type(wavespin),target,intent(inout) :: w_dense_ibz
   real(q),intent(in) :: phonon_freqs_ibz(:,:)
   type(type_info),intent(in) :: t_info
   type(info_struct),intent(in) :: info
   type(chemical_potential),intent(in) :: chempot
   type(bandgap_info),intent(in) :: bandgap(:)
   type(in_struct),intent(in) :: io

   ! local variables
   integer :: idx, i, idelta, inbands_sum, jcarrier, ikpt, iapprox, ikfbz
   integer,allocatable :: band_idx(:)

   PROFILING_START('elph_prepare_selfenergy_accumulators')

   ! keep track of the index of the bands
   call wdes_get_global_band_index(wdes_dense_ibz,band_idx)

   elph_accumulators%id_name = [character(len=15) :: "selfen_delta","nbands_sum","selfen_muij","selfen_approx"]
   ! store accumulators metadata to write to the hdf5 file
   elph_accumulators%id_size(1) = size(elph_set%selfen_delta)
   ! loop over number of bands to sum
   elph_accumulators%id_size(2) = size(elph_set%nbands_sum)
   ! loop over chemical potentials
   elph_accumulators%id_size(3) = size(chempot%selfen_muij,2)
   ! loop over approximations
   elph_accumulators%id_size(4) = size(elph_set%selfen_approx)

   idx = 0
   ! loop over broadening
   do idelta=1,size(elph_set%selfen_delta)
      ! loop over number of bands to sum
      do inbands_sum=1,size(elph_set%nbands_sum)
         ! loop over chemical potentials
         do jcarrier=1,size(chempot%selfen_muij,2)
         ! loop over approximations
         do iapprox=1,size(elph_set%selfen_approx)
            idx = idx+1
            select case (elph_set%driver)
               case (elph_driver_el)
                  if (idx>elph_accumulators%naccumulators_el) then
                     call vtutor%bug('Index is larger than number of electron accumulators',__FILE__,__LINE__)
                  endif
                  call elphon_prepare_selfenergy_accumulator_el(elph_accumulators%selfen_el(idx), idelta, inbands_sum, jcarrier, iapprox, &
                                                              elph_accumulators%id_size, elph_accumulators%id_name, elph_set, &
                                                              wdes_dense_ibz, w_dense_ibz, t_info, info, chempot, bandgap, &
                                                              phonon_freqs_ibz, band_idx, kpoints_inter)
               case (elph_driver_ph)
                  if (idx>elph_accumulators%naccumulators_ph) then
                     call vtutor%bug('Index is larger than number of phonon accumulators',__FILE__,__LINE__)
                  endif
                  call elphon_prepare_selfenergy_accumulator_ph(elph_accumulators%selfen_ph(idx), idelta, inbands_sum, jcarrier, iapprox, &
                                                              elph_accumulators%id_size, elph_accumulators%id_name, elph_set, &
                                                              wdes_dense_ibz, w_dense_ibz, t_info, info, chempot, bandgap, &
                                                              phonon_freqs_ibz, band_idx, kpoints_inter)
            end select
         enddo ! loop over approximations
         enddo ! loop over chemical potentials
      enddo ! loop over number of bands
   enddo ! loop over broadening

   ! write some information to the output
   do i=1,elph_accumulators%naccumulators_el
      call elph_accumulator_selfen_print(elph_accumulators%selfen_el(i),i,io%iu6)
   enddo
   do i=1,elph_accumulators%naccumulators_ph
      call elph_accumulator_selfen_print_ph(elph_accumulators%selfen_ph(i),i,io%iu6)
   end do
   PROFILING_STOP('elph_prepare_selfenergy_accumulators')

  end subroutine elphon_prepare_selfenergy_accumulators

  !> Private subroutine to initialize the self-energy accumulators for electrons
  subroutine elphon_prepare_selfenergy_accumulator_el(selfen, idelta, inbands_sum, jcarrier, iapprox, &
                                                      id_size, id_name, elph_set, wdes_dense_ibz, w_dense_ibz, t_info, &
                                                      info, chempot, bandgap, phonon_freqs_ibz, band_idx, &
                                                      kpoints_inter)
   use wave_struct_def, only: wavedes, wavespin
   use elphon_accumulators, only: elph_accumulator_selfen, elph_accumulator_selfen_init, &
        elph_accumulator_selfen_selector_ikpt, elph_accumulator_selfen_selector_energy_window, &
        elph_accumulator_selfen_selector_transport, elph_accumulator_selfen_selector_gaps, &
        elph_accumulator_selfen_sparsify, elph_accumulator_selfen_set_bandgap, &
        elph_accumulator_selfen_get_band_start_stop, elph_accumulator_selfen_get_band_start_stop_kp
   use elphon_base, only: elph_settings, elph_scattering_approx_crta
   use poscar, only: type_info
   use base, only: in_struct,info_struct
   use electron_transport, only: chemical_potential
   use bandgap_struct, only: bandgap_info
   use mkpoints_struct_def, only: kpoints_struct
   use mkpoints, only: kpoints_init_tetrahedron
   type(elph_accumulator_selfen),intent(inout) :: selfen
   integer,intent(in) :: idelta, inbands_sum, jcarrier, iapprox
   type(elph_settings),intent(in) :: elph_set
   integer,intent(in) :: id_size(4)
   character(len=15),intent(in) :: id_name(4)
   type(wavedes),intent(in) :: wdes_dense_ibz
   type(wavespin),target,intent(inout) :: w_dense_ibz
   type(type_info),intent(in) :: t_info
   type(info_struct),intent(in) :: info
   type(chemical_potential),intent(in) :: chempot
   type(bandgap_info),intent(in) :: bandgap(:)
   real(q),intent(in) :: phonon_freqs_ibz(:,:) !< phonon frequencies in the IBZ
   integer, intent(in) :: band_idx(:)
   type(kpoints_struct),intent(in) :: kpoints_inter

   call elph_accumulator_selfen_init(selfen,wdes_dense_ibz%nb_tot,&
                                    elph_set%nbands_sum(inbands_sum),wdes_dense_ibz%nkpts,&
                                    wdes_dense_ibz%ispin,wdes_dense_ibz%rspin,elph_set%selfen_nw,&
                                    elph_set%selfen_wrange,&
                                    elph_set%selfen_delta(idelta),elph_set%selfen_ntemps,&
                                    elph_set%selfen_temps,&
                                    elph_set%selfen_band_start,elph_set%selfen_band_stop,&
                                    elph_set%selfen_band_start_kp,elph_set%selfen_band_stop_kp,&
                                    elph_set%selfen_approx(iapprox),elph_set%selfen_broad_tol,&
                                    elph_set%transport,&
                                    t_info%nions,w_dense_ibz%celtot,wdes_dense_ibz%nb_totk,band_idx,&
                                    chempot%selfen_muij(:,jcarrier),chempot%selfen_nij(:,jcarrier),info%nelect,&
                                    elph_set%selfen_fan, elph_set%selfen_dfan, elph_set%selfen_dw,&
                                    [idelta,inbands_sum,jcarrier,iapprox],id_size,id_name)
   call elph_accumulator_selfen_set_bandgap(selfen,bandgap)
   call kpoints_init_tetrahedron(kpoints_inter,selfen%kcut24,selfen%kcut24_origin)
   selfen%kcut24_initialized = .true.
   if (size(elph_set%selfen_ikpt)>0) then
       call elph_accumulator_selfen_selector_ikpt(selfen,elph_set%selfen_ikpt)
   endif
   if (any(elph_set%selfen_energy_window>0)) then
       call elph_accumulator_selfen_selector_energy_window(selfen,w_dense_ibz%celtot,elph_set%selfen_energy_window)
   endif
   if (elph_set%transport) then
       ! The following variables are set:
       ! SELFEN%compute_mask SELFEN%band_start SELFEN%band_stop
       ! SELFEN%transport_emin, SELFEN%transport_emax, SELFEN%transport_nedos
       call elph_accumulator_selfen_selector_transport(selfen,w_dense_ibz%celtot,kpoints_inter, &
                                                       elph_set%transport_driver,elph_set%transport_emin_plot,elph_set%transport_emax_plot,&
                                                       elph_set%transport_nedos_plot>0,elph_set%transport_dfermi_tol,elph_set%transport_nedos,&
                                                       elph_set%transport_emin,elph_set%transport_emax)

       ! Now that we know for which nk states we need to compute the self-energy, we can use this to determine
       ! band_start_kp and band_stop_kp
       if ( elph_set%selfen_approx(iapprox)%scattering_approx == elph_scattering_approx_crta) then
           ! in the case of CRTA phonon_freqs_ibz are not allocated so we cannot determine the max phonon frequency
           ! so we cannot do anything here
       else
           call elph_accumulator_selfen_get_band_start_stop_kp(selfen,w_dense_ibz%celtot,phonon_freqs_ibz,elph_set%selfen_band_start_kp,elph_set%selfen_band_stop_kp)
       end if
   endif
   if (elph_set%selfen_gaps) then
       call elph_accumulator_selfen_selector_gaps(selfen,w_dense_ibz%celtot)
   endif
   call elph_accumulator_selfen_get_band_start_stop(selfen)
   call elph_accumulator_selfen_sparsify(selfen,w_dense_ibz%celtot)
  end subroutine elphon_prepare_selfenergy_accumulator_el

  !> Private subroutine to initialize the self-energy accumulators for phonons
  subroutine elphon_prepare_selfenergy_accumulator_ph(selfen, idelta, inbands_sum, jcarrier, iapprox, &
                                                      id_size, id_name, elph_set, wdes_dense_ibz, w_dense_ibz, t_info, &
                                                      info, chempot, bandgap, phonon_freqs_ibz, band_idx, &
                                                      kpoints_inter)
   use wave_struct_def, only: wavedes, wavespin
   use elphon_accumulators_ph, only: elph_accumulator_selfen_ph, elph_accumulator_selfen_init_ph, &
        elph_accumulator_selfen_selector_iqpt_ph, elph_accumulator_selfen_selector_energy_window_ph, &
        elph_accumulator_selfen_sparsify_ph
   use elphon_base, only: elph_settings
   use poscar, only: type_info
   use base, only: in_struct,info_struct
   use electron_transport, only: chemical_potential
   use bandgap_struct, only: bandgap_info
   use mkpoints_struct_def, only: kpoints_struct
   use mkpoints, only: kpoints_init_tetrahedron
   type(elph_accumulator_selfen_ph),intent(inout) :: selfen
   integer,intent(in) :: idelta, inbands_sum, jcarrier, iapprox
   type(elph_settings),intent(in) :: elph_set
   integer,intent(in) :: id_size(4)
   character(len=15),intent(in) :: id_name(4)
   type(wavedes),intent(in) :: wdes_dense_ibz
   type(wavespin),target,intent(inout) :: w_dense_ibz
   type(type_info),intent(in) :: t_info
   type(info_struct),intent(in) :: info
   type(chemical_potential),intent(in) :: chempot
   type(bandgap_info),intent(in) :: bandgap(:)
   real(q),intent(in) :: phonon_freqs_ibz(:,:)
   integer, intent(in) :: band_idx(:)
   type(kpoints_struct),intent(in) :: kpoints_inter

   call elph_accumulator_selfen_init_ph( &
       selfen, wdes_dense_ibz%nb_tot, &
       elph_set%nbands_sum(inbands_sum), wdes_dense_ibz%nkpts, &
       wdes_dense_ibz%ispin, wdes_dense_ibz%rspin, elph_set%selfen_nw, &
       elph_set%selfen_wrange, &
       elph_set%selfen_delta(idelta), elph_set%selfen_ntemps, &
       elph_set%selfen_temps, &
       elph_set%selfen_band_start, elph_set%selfen_band_stop, &
       elph_set%selfen_approx(iapprox), elph_set%selfen_broad_tol, &
       elph_set%transport, &
       t_info%nions, w_dense_ibz%celtot, phonon_freqs_ibz, wdes_dense_ibz%nb_totk, band_idx,&
       chempot%selfen_muij(:, jcarrier), chempot%selfen_nij(:, jcarrier), info%nelect, &
       elph_set%selfen_fan, elph_set%selfen_dfan, elph_set%selfen_dw, &
       [idelta, inbands_sum, jcarrier, iapprox], id_size, id_name)
   call kpoints_init_tetrahedron(kpoints_inter, selfen%kcut24, selfen%kcut24_origin)
   selfen%kcut24_initialized = .true.
   if (size(elph_set%selfen_ikpt) > 0) then
       call elph_accumulator_selfen_selector_iqpt_ph(selfen, elph_set%selfen_ikpt)
   end if
   if (any(elph_set%selfen_energy_window > 0)) then
       call elph_accumulator_selfen_selector_energy_window_ph( &
           selfen, phonon_freqs_ibz, elph_set%selfen_energy_window)
   end if
   call elph_accumulator_selfen_sparsify_ph(selfen, phonon_freqs_ibz)
  end subroutine elphon_prepare_selfenergy_accumulator_ph

  !> set number of bands to be computed non-self consistently for the
  !> electron-phonon calculation
  subroutine elphon_set_nbands(wdes,elph_set,nelect,t_info,max_nplwkp_tot)
   use wave_struct_def, only: wavedes
   use wave, only: get_nbands_serial
   use elphon_base, only: elph_settings
   use poscar, only: type_info
   use tutor, only: vtutor, isError, NumberOfElectrons, newArgument
   type(wavedes) :: wdes
   type(elph_settings) :: elph_set
   real(q),intent(in) :: nelect
   type(type_info),intent(in) :: t_info
   integer,intent(in) :: max_nplwkp_tot
   !local variables
   integer :: npar
   integer :: itut(2)
   character(20) :: ctut(1)

   ! set number of bands
   npar = wdes%nb_par
   if (elph_set%nbands>0) then
        elph_set%nbands=((elph_set%nbands+npar-1)/npar)*npar
        wdes%nb_tot=elph_set%nbands
   elseif (elph_set%nbands==-1) then
        ! Choose the same default that is used for the SCF calculation
        elph_set%nbands=get_nbands_serial(wdes%lnoncollinear,wdes%ispin,nelect,t_info)
        elph_set%nbands=((elph_set%nbands+npar-1)/npar)*npar
        wdes%nb_tot=elph_set%nbands
   elseif (elph_set%nbands==-2) then
        ! If NBANDS==-2 then the number of bands is set from the maximum number of plane-waves
        elph_set%nbands = max_nplwkp_tot
        elph_set%nbands = ((elph_set%nbands+npar-1)/npar)*npar
        wdes%nb_tot=elph_set%nbands
   endif
   wdes%nb_tot=((wdes%nb_tot+npar-1)/npar)*npar
   wdes%nbands=wdes%nb_tot/wdes%nb_par
   elph_set%nbands=wdes%nb_tot

   ! check if the number of bands is sufficient to accomodate all the electrons.
   ! if this is not he case, the routines that compute the fermi level simply hang.
   ! in some cases one might want to use a smaller number of bands just to compute
   ! for example electron-phonon matrix elements along a path but we will
   ! not allow this
   if (((.not. wdes%lnoncollinear).and.  nelect>real(wdes%nb_tot*2,kind=q)).or. &
           ((wdes%lnoncollinear).and.((nelect*2)>real(wdes%nb_tot*2,kind=q)))) then
      itut(1)=nelect ; itut(2)=wdes%nb_tot ; ctut(1)='ELPH_NBANDS'
      call vtutor%write(iserror, NumberOfElectrons, newArgument(ival = itut, cval = ctut))
   endif

  end subroutine elphon_set_nbands

  !> get wave function descriptor and wavefunction container in the dense IBZ
  !> This initializes the waverot and fbz2ibz structures
  !> TODO: the kpoints structure fbz2ibz map is updated here but it should be done somewhere else
  subroutine elphon_get_waverot_fbz2ibz(waverot,fbz2ibz,symm,latt_cur,latt_ini,grid_orig,ngdim_orig,wdes_dense_ibz,kpoints)
    use tutor, only: vtutor
    use string, only: str
    use base, only: symmetry
    use poscar_struct_def, only: latt
    use mgrid_struct_def, only: grid_3d
    use wave_struct_def, only : wavespin, wavedes
    use mkpoints_struct_def, only : kpoints_struct
    use wave_rotate, only: wave_rotator, wave_rotator_init, &
                           fbz2ibz_stars, fbz2ibz_stars_init, fbz2ibz_stars_dummy_init, kpoint_star, kstar_from_fbz2ibz
    use mkpoints, only: kpoints_is_regular, LineMode, vkpt_kindex
    type(wave_rotator),intent(out) :: waverot !< wave function rotator
    type(fbz2ibz_stars),intent(out) :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
    type(symmetry),intent(in) :: symm !< symmetry operations
    type(latt),intent(in) :: latt_cur !< current lattice
    type(latt),intent(in) :: latt_ini !< initial lattice
    type(grid_3d),intent(in) :: grid_orig !< original grid
    integer,intent(in) :: ngdim_orig !< minimum number of g-vectors to use when creating the storage for the WFs
    type(wavedes),intent(in) :: wdes_dense_ibz
    type(kpoints_struct),intent(inout) :: kpoints !< the fbz2ibz mapping in the kpoints structure will be updated
    !local variables
    type(kpoint_star) :: kstar
    integer :: ik,ik_star,ikfbz,irot
    logical :: timerev

    ! Intialize the wave function rotator
    timerev = symm%isym>=0.and..not.wdes_dense_ibz%lnoncollinear
    call wave_rotator_init(waverot,grid_orig,latt_cur,symm,timerev,wdes_dense_ibz%ispin)

    ! Generate mapping from the IBZ to the FBZ using stars
    if (kpoints%mode==LineMode) then
       call fbz2ibz_stars_dummy_init(fbz2ibz,symm,kpoints%nkpts,kpoints%vkpt)
    else
       call fbz2ibz_stars_init(fbz2ibz,symm,timerev,kpoints%nkpts,kpoints%vkpt)
    endif

    ! The following piece of code replaces the rotation operations that were used to
    ! expand the IBZ to the FBZ. This set of rotation operations is not unique so
    ! small numeric differences might appear depending on which rotation operation
    ! is used to generate the FBZ in particular for coarse k-point meshes.
    ! If the code bellow is executed we recover the reference resutls for the testsuite
    if (kpoints_is_regular(kpoints)) then
        ! use fbz2ibz_stars to store the index of the rotations that expand the FBZ
        kpoints%fbz2ibz(2,:) = 0
        if (allocated(kpoints%invmap)) then
           deallocate(kpoints%invmap)
        endif
        allocate(kpoints%invmap(48))
        kpoints%invmap = symm%invmap

        do ik=1,kpoints%nkpts
           call kstar_from_fbz2ibz(kstar,kpoints%vkpt(:,ik),fbz2ibz,ik)
           do ik_star=1,kstar%nkpts
              ikfbz = vkpt_kindex(kpoints,latt_cur,kstar%vkpt(:,ik_star))
              irot = kstar%irot(ik_star)
              kpoints%fbz2ibz(2,ikfbz) = sign(1,irot)*kpoints%invmap(abs(irot))
           enddo
        enddo
        if (any(kpoints%fbz2ibz(2,:)==0)) then
           call vtutor%error('Could not find the rotations to create the FBZ from the IBZ')
        endif
    endif

  end subroutine elphon_get_waverot_fbz2ibz

  !> Pinitialize wdes_dense_ibz and wdes_dense.
  !> Note that this subroutine will also change the layouts in GRID
  !> because of the kpoints in the wdes structure
  subroutine elphon_get_wdes_w(wdes_dense_ibz, w_dense_ibz, &
                               kpoints_inter, fbz2ibz, grid, latt_cur, latt_ini, ngdim)
   use lattice, only: latt
   use mkpoints, only: kpoints_struct
   use ini, only: register_allocate
   use wave, only: allocw_noplanewave, allocw, allocwdes, deallocwdes, wdes_set_nrplwv, wdes_count_indpw, &
                   count_indpw, wdes_set_kpoints, gen_layout
   use wave_struct_def, only: wavedes, wavespin
   use mgrid_struct_def, only: grid_3d
   use tutor, only: vtutor
   use wave_rotate, only: fbz2ibz_stars, kpoint_star, kstar_from_fbz2ibz

   type(wavedes),intent(inout) :: wdes_dense_ibz
   type(wavespin),target,intent(inout) :: w_dense_ibz
   type(kpoints_struct),intent(in) :: kpoints_inter
   type(fbz2ibz_stars),intent(in) :: fbz2ibz
   type(grid_3d),target,intent(inout) :: grid
   type(latt),intent(in) :: latt_cur
   type(latt),intent(in) :: latt_ini
   integer,intent(in) :: ngdim !< min number of g-vectors to use when creating wavefunction storage
   !local variables
   type(kpoint_star) :: kstar
   integer :: ik_ibz, ik_fbz, ik_star, nrplwv_max, nrplwv

   PROFILING_START('elph_get_wdes')

   ! Intialize wavefunction descriptor for the (dense) IBZ
   wdes_dense_ibz%grid => grid

   call wdes_set_kpoints(wdes_dense_ibz,kpoints_inter%vkpt,kpoints_inter%wtkpt)

#if 1
   ! TODO: Here I am generating the points in the FBZ such that I can call GEN_LAYOUT
   ! and find out what columns will be needed to perform the FFTs.
   ! There is a better way to do this by spliting parts of the GEN_LAYOUT routine
   ! that determine the number of columns and use the fbz2ibz structure to generate the
   ! FBZ on the fly.
   deallocate(wdes_dense_ibz%vkpt)
   allocate(wdes_dense_ibz%vkpt(3,fbz2ibz%nkfbz))
   ! loop over full kpoints and compute the max number of g-vectors
   ! This is a better version of the following line of code
   ! WDES_DENSE_IBZ%NGDIM = WDES_COUNT_INDPW(WDES_DENSE_IBZ,GRID,LATT_CUR%B, LATT_INI%B)
   ik_fbz=0
   do ik_ibz=1,wdes_dense_ibz%nkpts
       call kstar_from_fbz2ibz(kstar, kpoints_inter%vkpt(:,ik_ibz), &
                               fbz2ibz, ik_ibz)
       do ik_star=1,kstar%nkpts
          ik_fbz=ik_fbz+1
          wdes_dense_ibz%vkpt(:,ik_fbz) = kstar%vkpt(:,ik_star)
       enddo
   enddo
   if (ik_fbz/=fbz2ibz%nkfbz) call vtutor%bug('Expasion of the k-points from the IBZ to FBZ failed',__FILE__,__LINE__)
   wdes_dense_ibz%nkpts_for_gen_layout = fbz2ibz%nkfbz
#endif

   call gen_layout(grid, wdes_dense_ibz, latt_cur%b, latt_ini%b, -1, .true.)
   call gen_index(grid, wdes_dense_ibz, latt_cur%b, latt_ini%b, -1, -1, .true.)
#ifdef MPI
   call mapset(grid)
#endif

#if 1
   deallocate(wdes_dense_ibz%vkpt)
   allocate(wdes_dense_ibz%vkpt(3,wdes_dense_ibz%nkpts))
   wdes_dense_ibz%vkpt(:,:) = kpoints_inter%vkpt(:,:)
#endif

   if (size(wdes_dense_ibz%nplwkp,1)/=wdes_dense_ibz%nkpts) then
       call vtutor%bug('wrong dimensions of wdes_dense_ibz%nplwkp',__FILE__,__LINE__)
   endif

   ! loop over full kpoints and compute the max number of g-vectors
   ! This is a better version of the following line of code
   ! wdes_dense_ibz%ngdim = wdes_count_indpw(wdes_dense_ibz,grid,latt_cur%b, latt_ini%b)
   do ik_ibz=1,wdes_dense_ibz%nkpts
       call kstar_from_fbz2ibz(kstar, wdes_dense_ibz%vkpt(:,ik_ibz), &
                               fbz2ibz, ik_ibz)
       nrplwv_max=0
       do ik_star=1,kstar%nkpts
          nrplwv = count_indpw(grid, wdes_dense_ibz%enmax, latt_cur%b, latt_ini%b, kstar%vkpt(:,ik_star))
          nrplwv_max = max(nrplwv_max,nrplwv)
       enddo
       wdes_dense_ibz%nplwkp(ik_ibz) = nrplwv_max
   enddo
   wdes_dense_ibz%ngdim = maxval(wdes_dense_ibz%nplwkp)
   ! end loop over full kpoints

   wdes_dense_ibz%ngdim = max(ngdim, wdes_dense_ibz%ngdim)
   call wdes_set_nrplwv(wdes_dense_ibz,wdes_dense_ibz%ngdim)

   ! now we setup wdes again but with the correct dimensions for the IBZ
   call deallocwdes(wdes_dense_ibz,.true.)
   call allocwdes(wdes_dense_ibz,.true.)
   wdes_dense_ibz%nplwkp=0
   wdes_dense_ibz%luseinv=.false.
   wdes_dense_ibz%nindpw=0
   call gen_index(grid, wdes_dense_ibz, latt_cur%b, latt_ini%b, -1, -1, .true.)

   ! allocate a container for interpolated eigenvalues
   call allocw_noplanewave(wdes_dense_ibz, w_dense_ibz)
   call register_allocate(VMEM_ARR_SIZE(w_dense_ibz%celtot)+&
                          VMEM_ARR_SIZE(w_dense_ibz%fertot)+&
                          VMEM_ARR_SIZE(w_dense_ibz%auxtot),"elph_w_dense")

   PROFILING_STOP('elph_get_wdes')

  end subroutine elphon_get_wdes_w

  !> Generate wave function descriptor for the dense Brillouin zone
  !> Generate descriptor only for k-points for which we will need to compute the self-energy
  subroutine elphon_generate_wdes_dense(elph_set,selfen,latt_cur,latt_ini,&
                                        wdes_dense_ibz,kpoints_kidx,wdes_densek_ibz)
   use elphon_base, only: elph_settings
   use lattice, only: latt
   use elphon_accumulators, only: elph_accumulator_selfen, elph_accumulator_selfen_shallido
   use wave_struct_def, only: wavedes
   use wave, only: allocwdes
   type (elph_settings)  :: elph_set !< information for the electron-phonon calculation
   type (elph_accumulator_selfen) :: selfen(elph_set%selfen_naccumulators)
   type (latt) :: latt_cur
   type (latt) :: latt_ini
   type (wavedes) :: wdes_dense_ibz !< wavefunction descriptor for the dense ibz
   type (wavedes) :: wdes_densek_ibz !< wavefunction descriptor for the dense ibz
   integer,allocatable :: kpoints_kidx(:) !< mapping of the k-points included in the descriptor to the k-points in the ibz
   !local variables
   integer :: i, isp
   integer :: nkpts_k, ik_ibz
   logical :: shallido
   real(q) :: kpoints_vkpt(3,wdes_dense_ibz%nkpts)

   PROFILING_START('elph_wdes_dense')
   allocate(kpoints_kidx(wdes_dense_ibz%nkpts))
   kpoints_kidx=-1
   nkpts_k=0
   do ik_ibz=1,wdes_dense_ibz%nkpts
     ! check if need to do this k-point
     if ((elph_set%selfen_fan.or.elph_set%selfen_dw.or.elph_set%selfen_fans.or.elph_set%selfen_dws).and.(.not.elph_set%dry_run)) then
       shallido = .false.
       do isp=1,wdes_dense_ibz%ispin
           do i=1,elph_set%selfen_naccumulators
              shallido = shallido .or. &
                         elph_accumulator_selfen_shallido(selfen(i),ik_ibz,isp)
           enddo
       enddo
       if (.not.shallido) cycle
     endif
     nkpts_k=nkpts_k+1
     kpoints_vkpt(:,nkpts_k)=wdes_dense_ibz%vkpt(:,ik_ibz)
     kpoints_kidx(ik_ibz) = nkpts_k
   enddo

   wdes_densek_ibz=wdes_dense_ibz
   wdes_densek_ibz%grid=>wdes_dense_ibz%grid

   wdes_densek_ibz%nkdim= nkpts_k
   wdes_densek_ibz%nkpts= nkpts_k
   wdes_densek_ibz%nkpts_for_gen_layout= nkpts_k
   allocate(wdes_densek_ibz%vkpt(3,nkpts_k))
   wdes_densek_ibz%vkpt(:,:) = kpoints_vkpt(:,:nkpts_k)
   allocate(wdes_densek_ibz%wtkpt(nkpts_k))
   wdes_densek_ibz%wtkpt(:) = 1.0_q

   call allocwdes(wdes_densek_ibz,.true.)
   wdes_densek_ibz%nplwkp=0
   wdes_densek_ibz%luseinv=.false.
   call gen_index (wdes_densek_ibz%grid, wdes_densek_ibz, latt_cur%b, latt_ini%b, -1, -1, .true.)

   PROFILING_STOP('elph_wdes_dense')
  end subroutine elphon_generate_wdes_dense

  !>
  !> @brief Private subroutine to free communicators.
  !>
  !> THis is used by subroutine ELECTRON_PHONON.
  !>
  subroutine elphon_comm_free(comm_kbatch)
   use mpimy, only: communic
#ifdef MPI
   use mpimy, only: m_freec
#endif
   use wave_struct_def, only: wavedes

   !type(wavedes),intent(inout) :: wdes_kbatch
   type(communic),target,intent(inout) :: comm_kbatch !< communicator for each batch

   CALLMPI(m_freec(comm_kbatch))
   !CALLMPI(m_freec(wdes_kbatch%comm))
   !CALLMPI(m_freec(wdes_kbatch%comm_kinter))
   !CALLMPI(m_freec(wdes_kbatch%comm_kin))
   !CALLMPI(m_freec(wdes_kbatch%comm_inter))
   !CALLMPI(m_freec(wdes_kbatch%comm_inb))
  end subroutine elphon_comm_free

  !> @brief Private subroutine to prepare elph inputs.
  !>
  !> The following variable are set
  !>    + The "dimensions" in ELPH_SET, defined in elph_settings,
  !>      are read from phelel_params.hdf5 and from WDES_KBATCH
  !>
  !>    + lattice_pc in ELPH_POT from phelel_params.hdf5
  !>
  !>    + nbands in ELPH_SET = WDES_KBATCH%NB_TOT
  !>
  !>    + nbands_k in ELPH_SET, from INCAR,  = ELPH_SET%NBANDS or ELPH_SET%SELFEN_BAND_STOP
  subroutine elphon_prepare_elph_inputs(ELPH_SET, ELPH_POT, WDES_DENSE_IBZ, T_INFO, P, GRID, W)
   USE string, ONLY: str
   USE tutor, ONLY: vtutor
   USE elphon_base, ONLY: elph_settings
   USE wap, ONLY: read_phelel_dimensions_hdf5
   USE wave_struct_def, ONLY: wavedes, wavespin
   USE mgrid_struct_def, ONLY: grid_3d
   USE elphon_potential, ONLY: elphon_pot
   USE pseudo, ONLY: potcar
   USE poscar, ONLY: type_info
   type(elph_settings),intent(inout) :: ELPH_SET
   type(elphon_pot),intent(inout) :: ELPH_POT
   type(wavedes),intent(in) :: WDES_DENSE_IBZ
   type(potcar),intent(in) :: P(:)
   type(type_info),intent(in) :: T_INFO
   type(grid_3d),target :: GRID
   type(wavespin),target :: W
   !local
   integer :: FFT_MESH_WF(3)
   integer :: NATOMS_PC_CHECK
   integer :: nt

   ! Read the potential
   call read_phelel_dimensions_hdf5(ELPH_SET%NGRIDPOINTS,ELPH_SET%DET_P2S,ELPH_SET%POT_NCDIJ,&
                                    ELPH_SET%LMDIM,ELPH_SET%NATOMS_PC,ELPH_SET%NATOMS_SC,&
                                    ELPH_SET%FFT_MESH,ELPH_POT%LATTICE_PC,&
                                    ELPH_SET%PHONON_NATOMS_SC,ELPH_SET%PHONON_DET_P2S)

   ELPH_SET%ISPIN=WDES_DENSE_IBZ%ISPIN
   ELPH_SET%NPROJ=WDES_DENSE_IBZ%NPRO/WDES_DENSE_IBZ%NRSPINORS
   ELPH_SET%NPROJ_TOT=WDES_DENSE_IBZ%NPRO_TOT/WDES_DENSE_IBZ%NRSPINORS
   ELPH_SET%NBANDS=WDES_DENSE_IBZ%NB_TOT
   ELPH_SET%NBANDS_K=WDES_DENSE_IBZ%NB_TOT
   ELPH_SET%NPW=W%WDES%GRID%NGX*W%WDES%GRID%NGY*W%WDES%GRID%NGZ
   ELPH_SET%NPR_LOCAL=GRID%RL%NP
   ELPH_SET%NPW_LOCAL=GRID%MPLWV

   ELPH_SET%NTYP = T_INFO%NTYP
   ALLOCATE(ELPH_SET%NITYP(ELPH_SET%NTYP))
   ELPH_SET%NITYP = T_INFO%NITYP
   ALLOCATE(ELPH_SET%LMMAX(ELPH_SET%NTYP))
   DO NT=1,ELPH_SET%NTYP
      ELPH_SET%LMMAX(NT) = P(NT)%LMMAX
   ENDDO

   ELPH_SET%NATOMS_PC_LOCAL = SUM(WDES_DENSE_IBZ%NITYP)
   NATOMS_PC_CHECK = ELPH_SET%NATOMS_PC_LOCAL
   CALLMPI( M_sum_i( WDES_DENSE_IBZ%COMM_INB, NATOMS_PC_CHECK, 1))
   IF (NATOMS_PC_CHECK/=ELPH_SET%NATOMS_PC) THEN
      call vtutor%bug('Sum of local atoms '//str(NATOMS_PC_CHECK)//&
           ' is different from the total number of atoms '//str(ELPH_SET%NATOMS_PC),__FILE__,__LINE__)
   ENDIF

   ELPH_SET%NRSPINORS=WDES_DENSE_IBZ%NRSPINORS
   FFT_MESH_WF = [W%WDES%GRID%NGX,W%WDES%GRID%NGY,W%WDES%GRID%NGZ]
   if (ANY(FFT_MESH_WF/=ELPH_SET%FFT_MESH)) then
      call vtutor%error('Error: FFT grid from phelel '//str(ELPH_SET%FFT_MESH)//&
                        ' different from WF '//str(FFT_MESH_WF))
   endif
  end subroutine elphon_prepare_elph_inputs

  !> @brief Private subroutine to prepare wave_interp.
  subroutine elphon_prepare_wave_interp(WAVE_INTERP_SET, WAVE_INTERP, FBZ2IBZ, KPOINTS_INTER, COMM_KBATCH, &
                        KPOINTS, GRID, LATT_CUR, LATT_INI, &
                        T_INFO, NONL_S, W_ORIG, P, SYMM, IO, NGDIM, NB_TOT)
   USE base, ONLY: in_struct, symmetry
   USE mkpoints, ONLY: kpoints_struct
   USE mpimy, ONLY: communic
   USE lattice, ONLY: latt
   USE poscar, ONLY: type_info
   USE pseudo, ONLY: potcar
   USE mgrid_struct_def, ONLY: grid_3d
   USE nonl_high, ONLY: nonl_struct
   USE wave_struct_def, ONLY: wavespin
   USE wave_rotate, ONLY: fbz2ibz_stars
   USE wave_interpolate, ONLY: wave_interpolator, wave_interpolator_settings, wave_interpolator_init

   type(wave_interpolator_settings),intent(in) :: WAVE_INTERP_SET !< Wave interpolator settings
   type(wave_interpolator),intent(inout) :: WAVE_INTERP !< Wave interpolator
   type(fbz2ibz_stars),intent(in) :: FBZ2IBZ
   type(kpoints_struct),intent(in) :: KPOINTS_INTER
   type(communic),target,intent(in) :: COMM_KBATCH
   type(kpoints_struct),intent(in) :: KPOINTS
   type(grid_3d),target,intent(in) :: GRID
   type(latt),intent(in) :: LATT_CUR
   type(latt),intent(in) :: LATT_INI
   type(type_info),intent(in):: T_INFO
   type(nonl_struct),intent(in) :: NONL_S
   type(wavespin),target,intent(in) :: W_ORIG !< original WF from SCF run
   type(potcar),intent(in) :: P(:)
   type(symmetry),intent(in) :: SYMM
   type(in_struct),intent(in) :: IO
   integer,intent(in) :: NGDIM
   integer,intent(in) :: NB_TOT !number of bands to compute using the wave_interpolator
   integer :: IKPT
   integer :: NKPTS_BATCH
   integer :: STAR_SIZE, CPU_NEXT
   integer :: STAR_SIZE_COUNT(48)
   integer,allocatable :: WORKLOAD(:)
   real(q),allocatable :: KPOINTS_BATCH(:,:)
   integer,allocatable :: KPOINTS_INDEX(:)

   PROFILING_START('elph_prepare_wave_interp')

   ! Each node creates its own batch of k-points
   ! TODO: The distribution here is crucial to avoid load inbalance on high-symmetry systems!
   ! We should pre-compute the size of the star each k-point. The bigger the star
   ! the larger the amount of work needs to be done for that k-point.
   ! Sorting the k-points according to the size of the star and distributing in a
   ! round-bin fashion should lead to minimum load inbalance.
   ALLOCATE(WORKLOAD(COMM_KBATCH%NCPU))
   ALLOCATE(KPOINTS_BATCH(3,KPOINTS_INTER%NKPTS))
   ALLOCATE(KPOINTS_INDEX(KPOINTS_INTER%NKPTS))
   WORKLOAD=0
   NKPTS_BATCH=0
#if 0
   DO IKPT=1,KPOINTS_INTER%NKPTS
#ifdef MPI
      IF (MOD(IKPT,COMM_KBATCH%NCPU)/=COMM_KBATCH%NODE_ME-1) CYCLE
#endif
      ! If the batch is not complete add this k-point to batch
      NKPTS_BATCH = NKPTS_BATCH + 1
      KPOINTS_BATCH(:,NKPTS_BATCH) = KPOINTS_INTER%VKPT(:,IKPT)
      KPOINTS_INDEX(NKPTS_BATCH) = IKPT
      WORKLOAD(COMM_KBATCH%NODE_ME) = WORKLOAD(COMM_KBATCH%NODE_ME) + FBZ2IBZ%STAR_SIZE(IKPT)
   END DO
#else
   ! Count unique stars
   STAR_SIZE_COUNT = 0
   DO IKPT=1,KPOINTS_INTER%NKPTS
      STAR_SIZE = FBZ2IBZ%STAR_SIZE(IKPT)
      STAR_SIZE_COUNT(STAR_SIZE) = STAR_SIZE_COUNT(STAR_SIZE) + 1
   ENDDO

   ! Start by distributing the k-points with largest stars and slowly
   ! move to distributing the k-points with smaller and smaller ones
   ! each k-point is distributed to the CPU that has the lowest workload
   CPU_NEXT = 0
   DO STAR_SIZE=48,1,-1
      IF (STAR_SIZE_COUNT(STAR_SIZE)==0) CYCLE
      ! We have at least one start with this size
      DO IKPT=1,KPOINTS_INTER%NKPTS
         IF (FBZ2IBZ%STAR_SIZE(IKPT)/=STAR_SIZE) CYCLE
         ! Get CPU with minimum workload
         CPU_NEXT = MINLOC(WORKLOAD,1)
         WORKLOAD(CPU_NEXT) = WORKLOAD(CPU_NEXT) + FBZ2IBZ%STAR_SIZE(IKPT)
#ifdef MPI
         IF (CPU_NEXT/=COMM_KBATCH%NODE_ME) CYCLE
#endif
         ! If the batch is not complete add this k-point to batch
         NKPTS_BATCH = NKPTS_BATCH + 1
         KPOINTS_BATCH(:,NKPTS_BATCH) = KPOINTS_INTER%VKPT(:,IKPT)
         KPOINTS_INDEX(NKPTS_BATCH) = IKPT
      ENDDO
   ENDDO
#endif

   ! Initialize wave function interpolator
   CALL WAVE_INTERPOLATOR_INIT(WAVE_INTERP_SET, WAVE_INTERP, KPOINTS, KPOINTS_BATCH, GRID, LATT_CUR, LATT_INI, &
                          T_INFO, NONL_S, W_ORIG, P, SYMM, IO, &
                          NKPTS_BATCH, NGDIM=NGDIM, NB_TOT=NB_TOT, USE_ORIGINAL_GRID=.TRUE.)

   WAVE_INTERP%NKPTS_BATCH = NKPTS_BATCH
   WAVE_INTERP%KPOINTS_BATCH = KPOINTS_BATCH(:,:NKPTS_BATCH)
   WAVE_INTERP%KPOINTS_INDEX = KPOINTS_INDEX(:NKPTS_BATCH)
   ! The batch is complete so launch computation
   WRITE(*,'(A,I6,A,I6,A,I6,A,I6,A,I6,A,I6,A)') &
    'k-point batch [',WAVE_INTERP%NKPTS_BATCH,' / ',KPOINTS_INTER%NKPTS,&
      ' ] workload [',WORKLOAD(COMM_KBATCH%NODE_ME),' / ',SUM(FBZ2IBZ%STAR_SIZE),&
        '] on node [',COMM_KBATCH%NODE_ME,' / ',COMM_KBATCH%NCPU,']'
   DEALLOCATE(KPOINTS_INDEX)
   DEALLOCATE(KPOINTS_BATCH)
   DEALLOCATE(WORKLOAD)
   PROFILING_STOP('elph_prepare_wave_interp')

  end subroutine elphon_prepare_wave_interp

  !> Compute the electron velocity in eV . Ang
  !> Two possible approaches are used depending on elph_set%velocity_mode
  !> 1 - get_momentum routine in wap.F (computes derivatives of the AE orbitals, uses less memory)
  !> 2 - lrf_epsilon_mels in linear_optics.F (does not require derivatives of the AE orbitals uses more memory)
  subroutine elphon_get_velocity(VELOCITY, HAMILTONIAN, P, LMDIM, NBANDS, WF_NCDIJ, LATT_CUR, LATT_INI, SYMM, NONLR_S, NONL_S, &
                                 T_INFO, ELPH_SET, WDES_DENSE_IBZ, W_DENSE_IBZ, &
                                 WAVE_MAP, WDES_KBATCH, W_KBATCH, SV, CQIJ, CDIJ, &
                                 GRID, GRIDC, GRIDUS, C_TO_US, IRDMAX, INFO, IO, COMM_KBATCH)
   USE constant, only: felect
   USE base, ONLY: info_struct, symmetry, in_struct
   USE tutor, only: vtutor
   USE mpimy, ONLY: communic
   USE lattice, ONLY: latt
   USE pseudo, ONLY: potcar
   USE mkpoints, ONLY: kpoints_struct
   USE poscar, ONLY: type_info
   USE nonl_high, ONLY: nonl_struct, nonlr_struct
   USE wave_struct_def, ONLY: wavedes, wavedes1, wavespin, wavefuna
   USE wave_interpolate, ONLY: wave_mapper
   USE wave, ONLY: setwdes
   USE wave_high, ONLY: ELEMENTS
   USE elphon_base, ONLY: elph_settings
   USE wap, ONLY: get_momentum, getwave_rij, get_ccrij_a
   USE mlr_optic, ONLY: lrf_epsilon_mels
   USE subrot_cluster, ONLY: eigenf_cluster_pointer, find_deg_clusters, free_deg_clusters
   USE mgrid_struct_def, ONLY: grid_3d, transit
   USE hamil_struct_def, ONLY: ham_handle
   real(q),allocatable,intent(out) :: VELOCITY(:,:,:,:)
   type(ham_handle),intent(in) :: HAMILTONIAN
   type(potcar),intent(in) :: P(:)
   integer,intent(in) :: LMDIM
   integer,intent(in) :: NBANDS !< total number of bands for which the compute the velocity
   integer,intent(in) :: WF_NCDIJ !< spinor dimension of the potential according to the spinor dimension of the WFs
   type(latt),intent(in) :: LATT_CUR
   type(latt),intent(in) :: LATT_INI
   type(symmetry),intent(in) :: SYMM
   type(nonlr_struct),intent(in) :: NONLR_S
   type(nonl_struct),intent(in) :: NONL_S
   type(type_info),intent(in) :: T_INFO
   type(elph_settings),intent(in) :: ELPH_SET
   type(wavedes),intent(in) :: WDES_DENSE_IBZ
   type(wavespin),intent(in) :: W_DENSE_IBZ
   type(wave_mapper),intent(in) :: WAVE_MAP
   type(wavedes),intent(in) :: WDES_KBATCH
   type(wavespin),target,intent(in) :: W_KBATCH
   RGRID,intent(in) :: SV(:,:)
   OVERLAP :: CQIJ (:,:,:,:)
   OVERLAP :: CDIJ (:,:,:,:)
   type(grid_3d),intent(in) :: GRID
   type(grid_3d),intent(in) :: GRIDC
   type(grid_3d),intent(in) :: GRIDUS
   type(transit),intent(in) :: C_TO_US
   type(info_struct),intent(in) :: INFO
   type(in_struct),intent(in) :: IO
   type(communic),target,intent(in) :: COMM_KBATCH
   integer,intent(in) :: IRDMAX

   ! local variables
   real(q), parameter :: EL_TOL=4e-5 !< electron tolerance to detect degenerate eigenvalues (eV)
   LOGICAL :: LUSEPEAD_FALSE
   integer :: isp, ik, ib, nb_global, idir, ikcoarse, ikdense
   integer :: ibstart, ibstop, ndeg
   complex(q), allocatable :: CHAM(:,:,:)
   complex(q), allocatable :: CCRIJ(:,:,:,:,:)
   complex(q), allocatable :: CCRIJ_A(:,:,:,:,:)
   type(wavedes1) :: WDESK_IBZ
   type(wavefuna) :: WK_IBZ
   type(kpoints_struct) :: KPOINTS_DUMMY
   type(eigenf_cluster_pointer),pointer :: DEG_CLUSTER(:,:)
   GDEFS, allocatable :: CDER_BETWEEN_STATES(:,:,:,:,:)
   REAL(q), allocatable :: ENERGY_DER(:,:,:,:)
   REAL(q), allocatable :: VELOCITY_TMP(:,:)

   PROFILING_START('elph_get_velocity')
   IF (NBANDS /= WDES_DENSE_IBZ%NB_TOT) THEN
     call vtutor%bug('Incompatible number of bands of the wavefunciton',__FILE__,__LINE__)
   ENDIF
   IF (WF_NCDIJ /= WDES_DENSE_IBZ%NCDIJ) THEN
     call vtutor%bug('Incompatible spinor dimension of the wavefunction',__FILE__,__LINE__)
   ENDIF
   ALLOCATE(VELOCITY(NBANDS,WDES_DENSE_IBZ%NKPTS,WDES_DENSE_IBZ%ISPIN,3))
   VELOCITY=0.0_q

   SELECT CASE (ELPH_SET%VELOCITY_MODE)
     CASE(1)
       ALLOCATE(CHAM(NBANDS,NBANDS,3))
       ALLOCATE(CCRIJ( LMDIM, LMDIM, T_INFO%NIONS, WF_NCDIJ, 2))
       CCRIJ = CMPLX(0.0_q,0.0_q,q)
       call getwave_rij(WF_NCDIJ,LMDIM,P,T_INFO,CCRIJ )                ! get radial part of R matrix
       ! or CALL  getwave_RIJ_LOC(WDES_KBATCH,LMDIM,P,NONL_S, CCRIJ)              ! get radial part of R matrix 
       ALLOCATE(CCRIJ_A(LMDIM, LMDIM, T_INFO%NIONS, 3, WF_NCDIJ))
       CCRIJ_A = CMPLX(0.0_q,0.0_q,q)
       call get_ccrij_a(WF_NCDIJ,LMDIM,T_INFO,CCRIJ,P,CCRIJ_A)              ! multiply R matrix with angular part
       ! This should be avoided
       if (WF_NCDIJ==4) then
         CCRIJ_A(:,:,:,:,2)=0
         CCRIJ_A(:,:,:,:,3)=0
       endif
       ! or call get_CCRIJ_A_LOC(WDES_KBATCH,LMDIM,NONL_S,CCRIJ,P,CCRIJ_A)        ! multiply R matrix with angular part
       DO ISP=1,WDES_DENSE_IBZ%ISPIN
          DO IK=1,WAVE_MAP%NKPTS_BATCH
             IKCOARSE = WAVE_MAP%NKPTS_ORIG+IK
             IKDENSE = WAVE_MAP%KPOINTS_INDEX(IK)
             CALL SETWDES(WDES_KBATCH,WDESK_IBZ,IKCOARSE)
             WK_IBZ = ELEMENTS(W_KBATCH, WDESK_IBZ, ISP)
             CALL get_momentum(LATT_CUR,WK_IBZ,LMDIM,CCRIJ_A,CHAM(:,:,:))
             DO IDIR=1,3
                DO NB_GLOBAL=1,WDES_DENSE_IBZ%NB_TOT
                   VELOCITY(NB_GLOBAL,IKDENSE,ISP,IDIR) = REAL(CHAM(NB_GLOBAL,NB_GLOBAL,IDIR),q)*FELECT! keep only diagonal, but other are computed too
                ENDDO
             ENDDO ! loop directions
          ENDDO !loop kpoints
       ENDDO !loop spins
       DEALLOCATE(CHAM)
       DEALLOCATE(CCRIJ_A)
       DEALLOCATE(CCRIJ)
     CASE(2)
       LUSEPEAD_FALSE=.FALSE.
       ALLOCATE(CDER_BETWEEN_STATES(WDES_KBATCH%NB_TOT,WDES_KBATCH%NB_TOT,WDES_KBATCH%NKPTS,WDES_KBATCH%ISPIN,3))
       ALLOCATE(ENERGY_DER(WDES_KBATCH%NB_TOT, WDES_KBATCH%NKPTS, WDES_KBATCH%ISPIN, 3))
       NULLIFY(DEG_CLUSTER)
       CALL FIND_DEG_CLUSTERS(WDES_KBATCH, W_KBATCH, DEG_CLUSTER) ! DEGENERACY_THRESHOLD=2E-3_q
       CALL LRF_EPSILON_MELS(HAMILTONIAN, P, W_KBATCH%WDES, W_KBATCH, LATT_CUR, LATT_INI, T_INFO, INFO, IO, KPOINTS_DUMMY, SYMM, &
                              NONL_S, NONLR_S, IRDMAX, SV, DEG_CLUSTER, LUSEPEAD_FALSE, INFO%LOVERL, WDES_KBATCH%NB_TOT, &
                              LMDIM, CDIJ, CQIJ, GRID, GRIDC, GRIDUS, C_TO_US, ENERGY_DER, CDER_BETWEEN_STATES )
       CALL FREE_DEG_CLUSTERS(WDES_KBATCH, DEG_CLUSTER)
       DO IDIR=1,3
          DO ISP=1,WDES_DENSE_IBZ%ISPIN
             DO IK=1,WAVE_MAP%NKPTS_BATCH
                IKCOARSE = WAVE_MAP%NKPTS_ORIG+IK
                IKDENSE = WAVE_MAP%KPOINTS_INDEX(IK)
                ! prepare command to get non diagonal
                VELOCITY(:,IKDENSE,ISP,IDIR) = ENERGY_DER(:,IKCOARSE,ISP,IDIR)     ! in eV Ang
             ENDDO !loop kpoints
          ENDDO !loop spins
       ENDDO ! loop directions
       DEALLOCATE(CDER_BETWEEN_STATES)
       DEALLOCATE(ENERGY_DER)
   END SELECT
   CALLMPI( M_sum_d( COMM_KBATCH, VELOCITY, SIZE(VELOCITY)))

#if 1
   ! average over degenerate states
   ALLOCATE(VELOCITY_TMP(SIZE(VELOCITY,1),SIZE(VELOCITY,4)))
   DO ISP=1,WDES_DENSE_IBZ%ISPIN
      DO IK=1,WAVE_MAP%NKPTS_BATCH
        IB=1
        DO WHILE (IB<NBANDS)
            IBSTART = IB
            IBSTOP  = IB
            DO WHILE( abs(w_DENSE_IBZ%celtot(ibstop+1,ik,isp)-&
                          w_DENSE_IBZ%celtot(ibstart,ik,isp))<EL_TOL )
                IBSTOP = IBSTOP+1
                IF (IBSTOP>=NBANDS) EXIT
            ENDDO

            ! sum all the elements
            VELOCITY_TMP=0.0_q
            ndeg = 0
            do ib=ibstart,ibstop
                ndeg=ndeg+1
                VELOCITY_TMP(:,:) = VELOCITY_TMP(:,:) + VELOCITY(:,IK,ISP,:)
            enddo
            ! and set them
            do ib=ibstart,ibstop
                VELOCITY(:,IK,ISP,:) = VELOCITY_TMP(:,:)/ndeg
            enddo
            ! move to next eigenvalue
            ib=ibstop+1
        ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(VELOCITY_TMP)
#endif
   PROFILING_STOP('elph_get_velocity')

  end subroutine elphon_get_velocity

  !> Write the electron velocities
  subroutine elphon_write_velocity(elph_set,velocity,celtot,io)
   use elphon_base, only: elph_settings
   use base, only: in_struct
#ifdef VASP_HDF5
   use vhdf5, only: IH5OUTFILEID
#endif
   type(elph_settings),intent(in) :: ELPH_SET
   real(q),intent(in) :: velocity(:,:,:,:)
   complex(q),intent(in) :: celtot(:,:,:)
   type(in_struct), intent(in) :: io
   if (elph_set%write_textvel) call elphon_velocity_textwrite(velocity, celtot, io)
#ifdef VASP_HDF5
   if (elph_set%write_hdf5vel) call elphon_velocity_hdf5write(IH5OUTFILEID, velocity)
#endif
  end subroutine elphon_write_velocity

  !> Write out the electron velocities to a text file.
  !> Only master writes
  subroutine elphon_velocity_textwrite(velocity,celtot,io)
   use base, only: in_struct
   use tutor, only: vtutor
   use string, only: str
   real(q),intent(in) :: velocity(:,:,:,:)
   complex(q),intent(in) :: celtot(:,:,:)
   type(in_struct),intent(in) :: io
   ! local variables
   integer :: ispin, nkpoints, nbands
   integer :: idir,isp,ik,nb

   if (io%iu6<0) return

   !safety checks
   ispin = size(velocity,3)
   if (ispin/=size(celtot,3)) then
        call vtutor%bug('Incompatible number of spin channels for eigenvalues: '//str(size(celtot,3))//&
                        ' and band velocity matrix elements: '//str(ispin),__FILE__,__LINE__)
   endif

   nkpoints = size(velocity,2)
   if (nkpoints/=size(celtot,2)) then
        call vtutor%bug('Incompatible number of k-points for eigenvalues: '//str(size(celtot,2))//&
                        ' and band velocity matrix elements: '//str(nkpoints),__FILE__,__LINE__)
   endif

   nbands = size(velocity,1)
   if (nbands/=size(celtot,1)) then
        call vtutor%bug('Incompatible number of bands for eigenvalues: '//str(size(celtot,1))//&
                        ' and band velocity matrix elements: '//str(nbands),__FILE__,__LINE__)
   endif

   open(unit=99,file='velocity')
   do idir=1,3
      do isp=1,ispin
         do ik=1,nkpoints
            do nb=1,nbands
               write(99,'(4i3,2f12.6)') nb,ik,isp,idir,real(celtot(nb,ik,isp),q),&
                                          velocity(nb,ik,isp,idir)
            enddo
         enddo !ikpoints
      enddo !ispin
   enddo !idir
   close(99)
  end subroutine elphon_velocity_textwrite

#ifdef VASP_HDF5
  !> Write the electron group velocities to hdf5
  subroutine elphon_velocity_hdf5write(fileid,velocity,group)
   use base, only: in_struct
   use vhdf5
   integer(HID_T),intent(in) :: fileid
   real(q),intent(in) :: velocity(:,:,:,:)
   character(len=MAX_LEN_GROUP),optional,intent(in) :: group
   !local variables
   integer(HID_T) :: groupid
   character(len=MAX_LEN_GROUP) :: my_group

   my_group    = GRP_RESULTS//'/'//GRP_ELPHON//'/electrons'
   if (present(group))    my_group    = trim(group)

   ! open group
   VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))

   VH5_CHECK(vh5_write(groupid,'velocity',velocity))

   VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine elphon_velocity_hdf5write
#endif

  !> Compute the phonons in the IBZ
  subroutine elphon_compute_ibz_phonons(ELPH_IFC,ELPH_SET,KPOINTS,PHONON_FREQS_IBZ)
   USE elphon_base, ONLY: elph_settings
   USE tutor, ONLY: vtutor
   USE string, ONLY: str
   USE elphon_potential, ONLY: phonon_ifc, phonon_ifc_get_phonons
   USE mkpoints, ONLY: kpoints_struct
   USE elphon_accumulators, ONLY: PHON_TOL, THZTOEV
   type(phonon_ifc),intent(in) :: ELPH_IFC
   type(elph_settings),intent(in) :: ELPH_SET
   type(kpoints_struct),intent(in) :: KPOINTS
   real(q),intent(inout) :: PHONON_FREQS_IBZ(:,:)
   ! local variables
   integer :: ik
   complex(q),allocatable :: e_Q(:,:)
   ! Sanity check
   IF (KPOINTS%NKPTS/=SIZE(PHONON_FREQS_IBZ,2)) THEN
      call vtutor%error('Wrong dimensions of PHONON_FREQS_IBZ expected: '//&
                        str(SIZE(PHONON_FREQS_IBZ,2))//' got '//str(KPOINTS%NKPTS))
   ENDIF
   ! Compute phonon frequencies in IBZ, used for tetrahedron method to compute self-energy.
   ALLOCATE(e_Q(1:3*ELPH_IFC%NATOMS_PC,1:3*ELPH_IFC%NATOMS_PC))
   DO IK=1,KPOINTS%NKPTS
      CALL PHONON_IFC_GET_PHONONS(ELPH_IFC,KPOINTS%VKPT(:,IK),PHONON_FREQS_IBZ(:,IK),e_Q)
      IF (ANY(PHONON_FREQS_IBZ(:,IK)*THZTOEV<-PHON_TOL).AND..NOT.ELPH_SET%IGNORE_IMAG_PHONONS) THEN
          CALL vtutor%error('Negative frequencies found in the phonon dispersion. '//&
                 'At q-point = '//str(KPOINTS%VKPT(:,IK))//&
                 ' the 3 lowest freqs are '//str(PHONON_FREQS_IBZ(:3,IK))//&
                 '. You can silence this error by setting ELPH_IGNORE_IMAG_PHONONS '//&
                 'in the INCAR file.')
      ENDIF
   END DO
   DEALLOCATE(e_Q)
  end subroutine elphon_compute_ibz_phonons

  !> Subroutine to pre-calculate physical properties.
  !>
  !> Orbital interpolation: W_DENSE_IBZ, WAVE_INTERP
  subroutine elphon_run_precalculation(W_DENSE_IBZ, WAVE_INTERP, WAVE_MAP, &
                                       HAMILTONIAN, LATT_CUR, LATT_INI, T_INFO, &
                                       NONLR_S, NONL_S, W_KBATCH, LMDIM, P, SV, CQIJ, CDIJ, SYMM, &
                                       INFO, ELPH_SET, WDES_DENSE_IBZ, COMM_KBATCH, COMM_GLOBAL, &
                                       WDES_KBATCH, KPOINTS_INTER, IO)
   USE base, ONLY: symmetry, info_struct, in_struct
   USE lattice, ONLY: latt
   USE string, ONLY: str
   USE poscar, ONLY: type_info
   USE pseudo, ONLY: potcar
   USE mpimy, ONLY: communic
   USE mkpoints, ONLY: kpoints_struct
   USE ini, ONLY: register_allocate
   USE hamil_struct_def, ONLY: ham_handle
   USE nonl_high, ONLY: nonl_struct, nonlr_struct, phase
   USE wave_struct_def, ONLY: wavedes, wavespin, wavedes1
   USE wave_interpolate, ONLY: wave_interpolator, wave_interpolator_kpoints, wave_mapper, wave_mapper_from_wave_interp
   USE elphon_base, ONLY: elph_settings
#ifdef VASP_HDF5
   USE vhdf5, ONLY : HID_T
#endif
   USE elphon_kgrid, ONLY: elph_kgrid

   type(wavespin),target,intent(inout) :: W_DENSE_IBZ
   type(wave_interpolator),intent(inout) :: WAVE_INTERP
   type(wave_mapper),intent(inout) :: WAVE_MAP
   type(ham_handle),intent(in):: HAMILTONIAN
   type(latt),intent(in) :: LATT_CUR
   type(latt),intent(in) :: LATT_INI
   type(type_info),intent(in):: T_INFO
   type(nonlr_struct),intent(in) :: NONLR_S
   type(nonl_struct),intent(in) :: NONL_S
   type(wavespin),target,intent(inout) :: W_KBATCH
   integer,intent(in) :: LMDIM
   type(potcar),intent(in) :: P(:)
   RGRID,intent(in) :: SV(:,:)
   OVERLAP,intent(in) :: CQIJ(:,:,:,:)
   OVERLAP,intent(in) :: CDIJ(:,:,:,:)
   type(symmetry),intent(in) :: SYMM
   type(info_struct),intent(in) :: INFO
   type(in_struct),intent(in) :: IO
   type(elph_settings),intent(in) :: ELPH_SET
   type(wavedes),intent(in) :: WDES_DENSE_IBZ
   type(communic),target,intent(in) :: COMM_KBATCH
   type(communic),target,intent(in) :: COMM_GLOBAL
   type(wavedes),intent(inout) :: WDES_KBATCH
   type(kpoints_struct),intent(inout) :: KPOINTS_INTER
#ifdef VASP_HDF5
   integer(HID_T) :: groupid
#endif
   integer :: ierr

   ! local variables
   integer :: isp, ik, ikcoarse, ikdense, i, ibandmax
   real(q),allocatable :: ENERGY(:)
   real(q),allocatable :: DOS(:,:)
   real(q),allocatable :: DOSI(:,:)
   real(q),allocatable :: DOSPAR(:,:,:,:)
   real(q) :: EFERMI

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute orbitals on a dense k-point mesh !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   PROFILING_START('elph_wfs')
   !Deactivate IO when using wave interpolator IU6 = -1; IU0 = -1
   CALL WAVE_INTERPOLATOR_KPOINTS(WAVE_INTERP, HAMILTONIAN, &
              LATT_CUR, LATT_INI, T_INFO, NONLR_S, NONL_S, LMDIM, P, SV, &
              CQIJ, CDIJ, SYMM, INFO, -1, -1, WDES_COMM=WDES_DENSE_IBZ)
   PROFILING_STOP('elph_wfs')

   ! Now the orbitals and its descriptor live inside WAVE_INTERP
   W_KBATCH = WAVE_INTERP%W
   WDES_KBATCH = WAVE_INTERP%WDES

   ! Get mapping from k-points in the batch to the IBZ
   CALL WAVE_MAPPER_FROM_WAVE_INTERP(WAVE_MAP,WAVE_INTERP)

   ! Save eigenvalues
   W_DENSE_IBZ%CELTOT = 0.0_q
   WDES_DENSE_IBZ%NB_TOTK = 0
   DO ISP=1,WDES_DENSE_IBZ%ISPIN
      DO IK=1,WAVE_MAP%NKPTS_BATCH
         IKCOARSE = WAVE_MAP%NKPTS_ORIG+IK
         IKDENSE = WAVE_MAP%KPOINTS_INDEX(IK)
         W_DENSE_IBZ%CELTOT(:,IKDENSE,ISP) = W_KBATCH%CELTOT(:,IKCOARSE,ISP)
         WDES_DENSE_IBZ%NB_TOTK(IKDENSE,ISP) = WDES_KBATCH%NB_TOTK(IKCOARSE,ISP)
      END DO
   END DO
   PROFILING_START('elph_inbalance_wfs')
   CALLMPI( MPI_barrier( COMM_GLOBAL, IERR ) )
   PROFILING_STOP('elph_inbalance_wfs')
   ! Transmit all the eigenvalues across nodes
   CALLMPI( M_sum_z( COMM_KBATCH, W_DENSE_IBZ%CELTOT, SIZE(W_DENSE_IBZ%CELTOT)))
   CALLMPI( M_sum_i( COMM_KBATCH, WDES_DENSE_IBZ%NB_TOTK, SIZE(WDES_DENSE_IBZ%NB_TOTK)))


  end subroutine elphon_run_precalculation

  !> Chemical potential and number of electrons: CHEMPOT
  subroutine elphon_run_chempot(W_DENSE_IBZ, CHEMPOT, BANDGAP, &
                                LATT_CUR, LATT_INI, T_INFO, &
                                SYMM, &
                                INFO, ELPH_SET, WDES_DENSE_IBZ, COMM_KBATCH, COMM_GLOBAL, &
                                KPOINTS_INTER, IO)
   USE base, ONLY: symmetry, info_struct, in_struct
   USE lattice, ONLY: latt
   USE string, ONLY: str
   USE poscar, ONLY: type_info
   USE pseudo, ONLY: potcar
   USE mpimy, ONLY: communic
   USE mkpoints, ONLY: kpoints_struct
   USE ini, ONLY: register_allocate
   USE hamil_struct_def, ONLY: ham_handle
   USE nonl_high, ONLY: nonl_struct, nonlr_struct, phase
   USE wave_struct_def, ONLY: wavedes, wavespin, wavedes1
   USE wave_interpolate, ONLY: wave_interpolator, wave_interpolator_kpoints, wave_mapper, wave_mapper_from_wave_interp
   USE elphon_base, ONLY: elph_settings
#ifdef VASP_HDF5
   USE vhdf5, ONLY : HID_T, GRP_ELPHON, GRP_RESULTS, IH5OUTFILEID, vh5_error, VH5_GROUP_OPEN_OR_CREATE, vh5_group_close_writing, VH5_WRITE_DOS, VH5_WRITE_EIGENVAL
   USE electron_transport, ONLY: CHEMICAL_POTENTIAL_HDF5WRITE
#endif
   USE electron_transport, ONLY: chemical_potential,GET_Nij_MUij_DENSTA, WRITE_Nij_MUij_HEADER, WRITE_Nij_MUij_FINAL, &
                                 COMPUTE_DOSINTDOS_TETRAHEDRON
   USE elphon_kgrid, ONLY: elph_kgrid
   use bandgap_struct, only: bandgap_info, use_fermi_weights, use_local_fermi_energy
   use bandgap_tools, only: find_bandgap, write_to_unit_verbose

   type(wavespin),target,intent(inout) :: W_DENSE_IBZ
   type(chemical_potential) :: chempot
   type(bandgap_info),allocatable,intent(out) :: bandgap(:)
   type(latt),intent(in) :: LATT_CUR
   type(latt),intent(in) :: LATT_INI
   type(type_info),intent(in):: T_INFO
   type(symmetry),intent(in) :: SYMM
   type(info_struct),intent(in) :: INFO
   type(in_struct),intent(in) :: IO
   type(elph_settings),intent(in) :: ELPH_SET
   type(wavedes),intent(in) :: WDES_DENSE_IBZ
   type(communic),target,intent(in) :: COMM_KBATCH
   type(communic),target,intent(in) :: COMM_GLOBAL
   type(kpoints_struct),intent(inout) :: KPOINTS_INTER
#ifdef VASP_HDF5
   integer(HID_T) :: groupid
#endif
   integer :: ierr

   ! local variables
   integer :: isp, ik, ikcoarse, ikdense, i, ibandmax
   real(q),allocatable :: ENERGY(:)
   real(q),allocatable :: DOS(:,:)
   real(q),allocatable :: DOSI(:,:)
   real(q),allocatable :: DOSPAR(:,:,:,:)
   real(q) :: EFERMI

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Chemical potential and number of electrons !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   PROFILING_START('elph_carrier_fermi')
   CALL WRITE_Nij_MUij_header(ELPH_SET,IO%IU6)
   CALL GET_Nij_MUij_DENSTA(ELPH_SET, KPOINTS_INTER, LATT_CUR, WDES_DENSE_IBZ, W_DENSE_IBZ, &
                            INFO, CHEMPOT%SELFEN_Nij, CHEMPOT%SELFEN_MUij, EFERMI, COMM_GLOBAL)

   CALL WRITE_Nij_MUij_final(ELPH_SET%selfen_temps,CHEMPOT%SELFEN_NIJ,CHEMPOT%SELFEN_MUIJ,KPOINTS_INTER%SIGMA,EFERMI,IO%IU6)
   ! Find band-gap using band-tools
   !CALL FIND_BANDGAP(W_DENSE_IBZ, INFO%NELECT, use_fermi_weights, bandgap)
   CALL FIND_BANDGAP(W_DENSE_IBZ, INFO%NELECT, use_local_fermi_energy, bandgap)
   if (io%iu6>0) then
      call write_to_unit_verbose(io%iu6, bandgap)
   endif

#ifdef VASP_HDF5
   VH5_CHECK(VH5_GROUP_OPEN_OR_CREATE(IH5OUTFILEID, TRIM(GRP_RESULTS)//'/'//TRIM(GRP_ELPHON)//'/electrons/chemical_potential', GROUPID))
   CALL CHEMICAL_POTENTIAL_HDF5WRITE(ELPH_SET%selfen_temps,ELPH_SET%selfen_carrier_den,ELPH_SET%selfen_carrier_per_cell,ELPH_SET%selfen_mu,&
                                   EFERMI,CHEMPOT%SELFEN_NIJ,CHEMPOT%SELFEN_MUIJ,groupid)
   VH5_CHECK(vh5_group_close_writing(GROUPID))
#endif
   PROFILING_STOP('elph_carrier_fermi')

   ! Compute DOS and integral of DOS
   KPOINTS_INTER%EMAX = REAL(W_DENSE_IBZ%CELTOT(1,1,1),q)
   KPOINTS_INTER%EMIN = REAL(W_DENSE_IBZ%CELTOT(WDES_DENSE_IBZ%NB_TOTK(1,1),1,1),q)
   DO ISP=1,WDES_DENSE_IBZ%ISPIN
       DO IK=1,WDES_DENSE_IBZ%NKPTS
           IBANDMAX = WDES_DENSE_IBZ%NB_TOTK(IK,ISP)
           KPOINTS_INTER%EMAX = MAX(KPOINTS_INTER%EMAX,REAL(W_DENSE_IBZ%CELTOT(IBANDMAX,IK,ISP),q))
           KPOINTS_INTER%EMIN = MIN(KPOINTS_INTER%EMIN,REAL(W_DENSE_IBZ%CELTOT(1,IK,ISP),q))
       ENDDO
   ENDDO
   KPOINTS_INTER%EMAX = KPOINTS_INTER%EMAX+5.0_q
   KPOINTS_INTER%EMIN = KPOINTS_INTER%EMIN-5.0_q
   ALLOCATE(ENERGY(ELPH_SET%NEDOS))
   ALLOCATE(DOSI(ELPH_SET%NEDOS,WDES_DENSE_IBZ%ISPIN))
   ALLOCATE(DOS(ELPH_SET%NEDOS,WDES_DENSE_IBZ%ISPIN))
   ! TODO: compute this
   ALLOCATE(DOSPAR(1,1,1,1))
   DOSPAR = 0.0_q
   DO I=1,ELPH_SET%NEDOS
      ENERGY(I) = KPOINTS_INTER%EMIN + (KPOINTS_INTER%EMAX-KPOINTS_INTER%EMIN)/(ELPH_SET%NEDOS-1)*(I-1)
   END DO
   CALL COMPUTE_DOSINTDOS_TETRAHEDRON(KPOINTS_INTER, WDES_DENSE_IBZ, W_DENSE_IBZ%CELTOT, LATT_CUR, ENERGY, DOS, DOSI, COMM_GLOBAL)

   ! Write the eigenvalues and DOS to HDF5 file
#ifdef VASP_HDF5
   VH5_CHECK(VH5_GROUP_OPEN_OR_CREATE(IH5OUTFILEID, TRIM(GRP_RESULTS)//'/electron_phonon', GROUPID))
   CALL VH5_WRITE_DOS(GROUPID, WDES_DENSE_IBZ, KPOINTS_INTER, DOS, DOSI, DOSPAR, EFERMI, &
                      T_INFO%NIONP,-1,GROUP="electrons", SUBGROUP="dos")
   CALL VH5_WRITE_EIGENVAL(GROUPID, WDES_DENSE_IBZ, W_DENSE_IBZ, KPOINTS_INTER, &
                           GROUP="electrons", SUBGROUP="eigenvalues")
   VH5_CHECK(vh5_group_close_writing(GROUPID))
#endif
   DEALLOCATE(DOS,DOSI,DOSPAR)

  end subroutine elphon_run_chempot

  !> Compute phonon limited mobility using the electron-phonon lifetimes
  SUBROUTINE TRANSPORT_ELPHON(ELPH_SET,SELFEN,VELOCITY,LATT_CUR,KGRID,KPOINTS_INTER,W_DENSE_IBZ,WDES_DENSE_IBZ,BANDGAP,SYMM,INFO,IO,COMM)
    USE mpimy, ONLY: communic
    USE base, ONLY: energy, info_struct, in_struct, symmetry
    USE lattice, ONLY: latt
    USE mkpoints, ONLY: kpoints_struct, tetirr_star
    USE m_tetra_hash, ONLY: tetirr_hash
    USE elphon_base, ONLY: elph_settings, elph_scattering_approx_crta
    USE wave_struct_def, ONLY: wavedes, wavespin
    USE elphon_accumulators, ONLY: elph_accumulator_selfen, elph_accumulator_selfen_get
    USE tutor, ONLY: vtutor
    USE string, ONLY: str
#ifdef VASP_HDF5
    USE vhdf5, ONLY: IH5OUTFILEID
    USE electron_transport, ONLY: transport_hdf5write, transports_hdf5write
#endif
    USE elphon_kgrid, ONLY: elph_kgrid
    USE electron_transport, ONLY: transport_type, transport_init, transport_compute, transport_compute_gausslegendre, &
                                  transport_free, relaxation_time, transport_print_header, &
                                  transport_print, transport_print_results,&
                                  print_transport_results_outcar, transport_set_results_array_line,&
                                  transport_compute_tau_average
    USE bandgap_struct, ONLY: bandgap_info
    type(elph_settings),intent(in) :: ELPH_SET
    type(elph_accumulator_selfen),intent(in) :: SELFEN(:)
    real(q),intent(inout) :: VELOCITY(:,:,:,:) !<velocity matrix elements
    type(latt),intent(in) :: LATT_CUR
    type(elph_kgrid),intent(in) :: kgrid !< kspclib related variables
    type(kpoints_struct),intent(inout) :: KPOINTS_INTER
    type(wavedes),intent(in) :: WDES_DENSE_IBZ !< wavefunction descriptor for the dense IBZ
    type(wavespin),intent(in) :: W_DENSE_IBZ !< wavefunction for the dense IBZ
    type(bandgap_info),intent(in) :: bandgap(:)
    type(symmetry),intent(in) :: SYMM
    type(info_struct),intent(in) :: INFO
    type(in_struct),intent(in) :: IO
    type(communic),intent(in) :: COMM
    ! local variables
    INTEGER :: I, J
    INTEGER :: NTET
    TYPE (transport_type) :: TRANSPORT
    COMPLEX(q), ALLOCATABLE :: SELF_ENERGY(:,:,:)
    REAL(q), ALLOCATABLE :: TAU(:,:,:)
    REAL(q), ALLOCATABLE :: TRANSPORT_RESULTS(:,:,:)
    INTEGER,ALLOCATABLE,TARGET :: IDTET(:,:)
    REAL(q), PARAMETER :: HBAR=6.582119569e-16_q

    IF (.NOT.ELPH_SET%TRANSPORT) RETURN

    ! inputs: ELPH_SET, VELOCITY, WDES_DENSE_IBZ, WDES_DENSE_IBZ, LATT_CUR%A, kgrid,  and SELFEN if we use the elph interactionto compute TAU
    IF(IO%IU0>0) write(IO%IU0,*) 'transport calculation'
    IF(IO%IU6>0) THEN
        write(IO%IU6,*) ''
        write(IO%IU6,*) 'transport calculation'
        write(IO%IU6,*) '====================='
    ENDIF

    IF (WDES_DENSE_IBZ%NKPTS/=KPOINTS_INTER%NKPTS) THEN
        call vtutor%bug('Inconsistent number of k-points in KPOINTS NKPTS='//str(KPOINTS_INTER%NKPTS)//&
                        ' and WDES NKPTS='//str(WDES_DENSE_IBZ%NKPTS)//' structures',__FILE__,__LINE__)
    ENDIF
    PROFILING_START('transport_driver')

    ALLOCATE(TAU(WDES_DENSE_IBZ%NB_TOT,WDES_DENSE_IBZ%NKPTS,WDES_DENSE_IBZ%ISPIN))

    ! If the tetrahedra are not allocated then compute them
    IF (.NOT.ASSOCIATED(KPOINTS_INTER%IDTET)) THEN
        !CALL TETIRR_STAR(KPOINTS,SYMM,LATT_CUR%B,NTET,IDTET,-1)
        !WRITE(*,*) 'STAR', NTET, KPOINTS%NKPX, KPOINTS%NKPY, KPOINTS%NKPZ
        !DEALLOCATE(IDTET)
        CALL TETIRR_HASH(KPOINTS_INTER,LATT_CUR%B,NTET,IDTET,-1)
        !WRITE(*,*) 'HASH', NTET, KPOINTS_INTER%NKPX, KPOINTS_INTER%NKPY, KPOINTS_INTER%NKPZ
        !COPY DATA
        KPOINTS_INTER%NTET = NTET
        ALLOCATE(KPOINTS_INTER%IDTET(0:4,NTET))
        KPOINTS_INTER%IDTET = IDTET
        DEALLOCATE(IDTET)
    ENDIF
    IF(IO%IU6>0) THEN
        WRITE(IO%IU6,*)
        WRITE(IO%IU6,*) 'Using', KPOINTS_INTER%NTET, 'irreducible tetrahedra out of a total of', 6*KPOINTS_INTER%ISNF(1,1)*KPOINTS_INTER%ISNF(2,2)*KPOINTS_INTER%ISNF(3,3)
        WRITE(IO%IU6,*)
    ENDIF

    ! the number of temperatures must be the same on all accumualtors for now
    ! (this can can easily be changed if needed)
    IF (ANY(SELFEN(:)%NTEMPS/=SELFEN(1)%NTEMPS)) THEN
        call vtutor%bug('Number of temperatures must be the same on all accumulators',__FILE__,__LINE__)
    ENDIF
    ALLOCATE(TRANSPORT_RESULTS(7,SELFEN(1)%NTEMPS,SIZE(SELFEN)))

#ifdef VASP_HDF5
    CALL TRANSPORTS_HDF5WRITE(IH5OUTFILEID, SELFEN(1)%ID_NAME, SELFEN(1)%ID_SIZE)
#endif
    DO I=1,SIZE(SELFEN)
        ! set dimensions and basic variables in the TRANSPORT structure
        CALL TRANSPORT_INIT(TRANSPORT,LATT_CUR%A,SYMM%NROT,SYMM%ISYMOP,ELPH_SET%SELFEN_TEMPS(:), &
                            SELFEN(I)%CARRIER_PER_CELL(:), INFO%NELECT, SELFEN(I)%APPROX, WDES_DENSE_IBZ%RSPIN, SELFEN(I)%EFERMI(:), &
                            ELPH_SET%TRANSPORT_EMIN_PLOT,ELPH_SET%TRANSPORT_EMAX_PLOT, &
                            ELPH_SET%TRANSPORT_EMIN,ELPH_SET%TRANSPORT_EMAX,SELFEN(I)%TRANSPORT_NEDOS, &
                            ELPH_SET%TRANSPORT_NEDOS_PLOT, &
                            ELPH_SET%TRANSPORT_DFERMI_TOL, BANDGAP, SELFEN(I)%ID_IDX, SELFEN(I)%ID_NAME, TAU, VELOCITY)

        ! calculate transport quantities one temperature at a time
        DO J=1,SELFEN(I)%NTEMPS

            IF (SELFEN(I)%APPROX%SCATTERING_APPROX == ELPH_SCATTERING_APPROX_CRTA) THEN
                ! in eV^{-1} if ELPH_SET%TRANSPORT_RELAXATION_TIME in second in INCAR
                ! typical value is TAU=10^{-14} s
                TAU=ELPH_SET%TRANSPORT_RELAXATION_TIME/(HBAR)
            ELSE
                ALLOCATE(SELF_ENERGY(SELFEN(I)%BAND_START:SELFEN(I)%BAND_STOP,SELFEN(I)%NKPOINTS,SELFEN(I)%NSPIN))
                CALL ELPH_ACCUMULATOR_SELFEN_GET(SELFEN(I),J,1,SELF_ENERGY)  ! XXX why 1 ???

                ! get TAU = 1/(2 |Im Sigma|) in eV^{-1} (\tau/\hbar in fact)
                CALL RELAXATION_TIME(SELFEN(I)%BAND_START, SELFEN(I)%BAND_STOP,SELFEN(I)%NBANDS, &
                                     SELFEN(I)%NKPOINTS,SELFEN(I)%NSPIN, &
                                     SELF_ENERGY,TAU)
                DEALLOCATE(SELF_ENERGY)
            ENDIF

            CALL TRANSPORT_COMPUTE_TAU_AVERAGE(TRANSPORT,W_DENSE_IBZ%CELTOT,J)
            SELECT CASE (ELPH_SET%TRANSPORT_DRIVER)
               CASE (1)
                  CALL TRANSPORT_COMPUTE(TRANSPORT,J,WDES_DENSE_IBZ,W_DENSE_IBZ,KPOINTS_INTER,COMM)
               CASE (2)
                  CALL TRANSPORT_COMPUTE_GAUSSLEGENDRE(TRANSPORT,J,WDES_DENSE_IBZ,W_DENSE_IBZ,KPOINTS_INTER,COMM)
               CASE DEFAULT
                  CALL vtutor%error('Invalid option for ELPH_TRANSPORT_DRIVER, should be 1 or 2')
            END SELECT

            !CALL TRANSPORT_PRINT_RESULTS(TRANSPORT,1,1,IO%IU0)
            !CALL TRANSPORT_PRINT_RESULTS(TRANSPORT,1,1,IO%IU6)

            ! store the results on an array for output in the output and OUTCAR
            CALL TRANSPORT_SET_RESULTS_ARRAY_LINE(TRANSPORT,1,1,J,TRANSPORT_RESULTS(:,J,I))

        ENDDO ! temperatures
        !write information about transport calculation to outcar
        CALL TRANSPORT_PRINT(TRANSPORT,I,IO%IU6)

#ifdef VASP_HDF5
        CALL TRANSPORT_HDF5WRITE(TRANSPORT,IH5OUTFILEID,I)
#endif
        CALL TRANSPORT_FREE(TRANSPORT)
    ENDDO ! self-energy objects

    ! write data to outcar and stdout
    DO I=1,SIZE(SELFEN)
        CALL TRANSPORT_PRINT_HEADER(IO%IU0,I)
        CALL TRANSPORT_PRINT_HEADER(IO%IU6,I)
        SELECT CASE (ELPH_SET%TRANSPORT_DRIVER)
            CASE(1)
                CALL PRINT_TRANSPORT_RESULTS_OUTCAR(TRANSPORT_RESULTS(:,:,I),'Linear grids',IO%IU0)
                CALL PRINT_TRANSPORT_RESULTS_OUTCAR(TRANSPORT_RESULTS(:,:,I),'Linear grids',IO%IU6)
            CASE(2)
                CALL PRINT_TRANSPORT_RESULTS_OUTCAR(TRANSPORT_RESULTS(:,:,I),'Gauss-Legendre grids',IO%IU0)
                CALL PRINT_TRANSPORT_RESULTS_OUTCAR(TRANSPORT_RESULTS(:,:,I),'Gauss-Legendre grids',IO%IU6)
        END SELECT
    ENDDO

    DEALLOCATE(TRANSPORT_RESULTS)

    DEALLOCATE(TAU)

    IF (IO%IU6>0) THEN
        WRITE(IO%IU6,'(A)') '------------------------ end of transport driver reached -----------------------------------------------'
        WRITE(IO%IU6,*)
    ENDIF
    PROFILING_STOP('transport_driver')

  END SUBROUTINE TRANSPORT_ELPHON


  !> Displacement derivative of projectors calculation: DPROJ_KBATCH
  subroutine elphon_compute_projections(WDES_KBATCH,W_KBATCH,NONL_S,WAVE_MAP,DPROJ_KBATCH,LATT_CUR)
   use ini, only: register_allocate
   use wave_struct_def, only: wavedes, wavedes1, wavespin
   use wave_high, only: element
   USE nonl_struct_def, ONLY: nonl_struct
   USE nonl_high, ONLY: phase
   use wave_interpolate, only: wave_mapper
   USE wave, ONLY: setwdes
   use lattice, only: latt
   USE wap, ONLY: get_dproj
   type(wavedes),intent(in) :: WDES_KBATCH
   type(wavespin),intent(in) :: W_KBATCH
   type(nonl_struct),intent(in) :: NONL_S
   type(latt),intent(in) :: LATT_CUR
   type(wave_mapper),intent(in) :: WAVE_MAP
   GDEF,allocatable,intent(out) :: DPROJ_KBATCH(:,:,:,:,:,:)
   ! local variables
   GDEF,allocatable :: DPROJ(:,:,:)
   type(wavedes1) :: WDESKP_IBZ
   INTEGER :: IK, ISP, ISPIN, NBK, NPROJ, NRSPINORS

   NRSPINORS = WDES_KBATCH%NRSPINORS
   NPROJ = WDES_KBATCH%NPRO/NRSPINORS
   ISPIN = WDES_KBATCH%ISPIN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Displacement derivative of projectors !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ! compute projections for the WFs in this batch
   ALLOCATE(DPROJ(NPROJ,0:3,NRSPINORS))
   ALLOCATE(DPROJ_KBATCH(NPROJ,0:3,NRSPINORS,WDES_KBATCH%NBANDS,&
                               WAVE_MAP%NKPTS_BATCH,ISPIN))
   CALL REGISTER_ALLOCATE(VMEM_ARR_SIZE(DPROJ_KBATCH), "dproj")
   DO IK=1,WAVE_MAP%NKPTS_BATCH
      CALL SETWDES(WDES_KBATCH,WDESKP_IBZ,WAVE_MAP%NKPTS_ORIG+IK)
      CALL PHASE(WDES_KBATCH,NONL_S,WDESKP_IBZ%NK)
      DO ISP=1,ISPIN
         !$omp parallel do private(NBK,DPROJ)
         DO NBK=1,WDES_KBATCH%NBANDS
            CALL get_DPROJ(LATT_CUR,NONL_S,&
                                 ELEMENT(W_KBATCH,WDESKP_IBZ,NBK,ISP),WDESKP_IBZ,DPROJ)
            DPROJ_KBATCH(:,:,:,NBK,IK,ISP) = DPROJ
         ENDDO
         !$omp end parallel do
      ENDDO
   END DO
   DEALLOCATE(DPROJ)
  end subroutine elphon_compute_projections

  !> Create wave windows for the plane-waves and projections
  subroutine create_wave_windows(wavwin_cw,wavwin_cproj,nproj,nrspinors,nkpts,&
                                 wave_map,w_kbatch,dproj_kbatch,wdes_kbatch,comm_global)
      use wave_windower, only: wave_window, wave_window_init, wave_window_start, wave_window_store, get_bks_wavedes, wave_cache_init
      use wave_struct_def, only: wavedes, wavespin
      use wave_interpolate, only: wave_mapper
      use mpimy, only: communic
      use tutor, only: vtutor
      use iso_c_binding, only: c_loc, c_f_pointer
      type (wave_window) :: wavwin_cw !< wave window
      type (wave_window) :: wavwin_cproj !< wave window
      type (wave_mapper),intent(in) :: wave_map
      integer,intent(in) :: nproj
      integer,intent(in) :: nkpts !< total number of distributed k-points
      integer,intent(in) :: nrspinors
      GDEF,allocatable,target :: dproj_kbatch(:,:,:,:,:,:)
      type(wavespin),intent(in) :: w_kbatch
      type(wavedes),intent(in) :: wdes_kbatch
      type(communic) :: comm_global
      !local variables
      complex(q),pointer :: cw_ptr(:,:,:,:)
      complex(q),pointer :: cproj_ptr(:,:,:,:)
      integer,allocatable :: bks_idx(:,:,:,:)

      call wave_window_init(wavwin_cw,wdes_kbatch%nb_tot,nkpts,wdes_kbatch%ispin,&
                            wdes_kbatch%nrplwv,comm_global,wdes_kbatch%comm_inb)
      call get_bks_wavedes(wdes_kbatch,bks_idx,wave_map%nkpts_batch,wave_map%kpoints_index)
      ! Here is check if the number of k-points in the batch is larger than zero and
      ! create a pointer. If the number o k-points is zero then I get problems with
      ! the bounds check so I create a dummy pointer (needed to avoid segfault in wave_window_store)
      if (wave_map%nkpts_batch>0) then
        call c_f_pointer(c_loc(w_kbatch%CW(1,1,wave_map%nkpts_orig+1,1)),cw_ptr,&
                         [wdes_kbatch%nrplwv,wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin])
      else
        ! dummy pointer
        call c_f_pointer(c_loc(w_kbatch%CW),cw_ptr,&
                         [wdes_kbatch%nrplwv,wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin])
      endif
      if (size(w_kbatch%CW,1)/=wdes_kbatch%nrplwv) then
        call vtutor%bug('number of plane-waves',__FILE__,__LINE__)
      endif
      if (size(w_kbatch%CW,2)/=wdes_kbatch%nbands) then
        call vtutor%bug('number of bands',__FILE__,__LINE__)
      endif
      if (size(w_kbatch%CW,3)/=wave_map%nkpts_batch) then
        call vtutor%bug('number of k-points',__FILE__,__LINE__)
      endif
      call wave_window_store(wavwin_cw,bks_idx,cw_ptr)
      call wave_window_start(wavwin_cw)
      deallocate(bks_idx)

      call wave_window_init(wavwin_cproj,wdes_kbatch%nb_tot,nkpts,wdes_kbatch%ispin,&
                            nproj*4*nrspinors,comm_global,wdes_kbatch%comm_inb)
      call get_bks_wavedes(wdes_kbatch,bks_idx,wave_map%nkpts_batch,wave_map%kpoints_index)
      call c_f_pointer(c_loc(dproj_kbatch),cproj_ptr,&
           [nproj*4*nrspinors,wdes_kbatch%nbands,wave_map%nkpts_batch,wdes_kbatch%ispin])
      call wave_window_store(wavwin_cproj,bks_idx,cproj_ptr)
      call wave_window_start(wavwin_cproj)
      deallocate(bks_idx)
  end subroutine create_wave_windows

  subroutine create_wave_caches(cache_cw,cache_cproj,wf_cache_mb,nproj,nproj_tot,nrspinors,nkpts,wdes_kbatch)
      use wave_windower, only: wave_cache, wave_cache_init
      use wave_struct_def, only: wavedes
      type (wave_cache) :: cache_cw !< wave window
      type (wave_cache) :: cache_cproj !< wave window
      real(q),intent(in) :: wf_cache_mb
      integer,intent(in) :: nproj
      integer,intent(in) :: nproj_tot
      integer,intent(in) :: nkpts
      integer,intent(in) :: nrspinors
      type(wavedes),intent(in) :: wdes_kbatch
      !local variables
      real(q) :: cproj_cache_mb
      cproj_cache_mb = wf_cache_mb/wdes_kbatch%nrplwv * nproj_tot*4*nrspinors
      call wave_cache_init(cache_cw,wdes_kbatch%nb_tot,nkpts,wdes_kbatch%ispin,wdes_kbatch%nrplwv,wf_cache_mb)
      call wave_cache_init(cache_cproj,wdes_kbatch%nb_tot,nkpts,wdes_kbatch%ispin,nproj*4*nrspinors,cproj_cache_mb)
  end subroutine create_wave_caches
END MODULE

