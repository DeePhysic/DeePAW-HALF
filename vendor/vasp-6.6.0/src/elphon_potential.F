#include "symbol.inc"

!>
!> @brief Module to handle the electron-phonon potential
!> Most of the code here was originally written by Laurent Chaput.
!> Uses data produced by phelel written by Atsushi Togo.
!> Treatment of the long-range potential, refactoring,
!> hdf5 reading and some optimizations by Henrique Miranda.
!>
module elphon_potential
    use prec
    use poscar_struct_def
    use mgrid_struct_def
    use elphon_potential_struct
    use pseudo, only: potcar
    use lattice, only: lattic
#ifdef shmem_elphon
    use mpi_shmem
#endif

    implicit none

    private :: get_phonons_zero_disp
    complex(q),parameter :: ci = cmplx(0.0_q,1.0_q,q)
    real(q), parameter :: supression_threshold = log(1e-10_q) !< todo: check this value
    contains

    !> @brief initialize the phonon force constants
    subroutine phonon_ifc_init(self,det_p2s,natoms_pc,natoms_sc,make_hermitian)
        type(phonon_ifc) :: self
        integer,intent(in) :: det_p2s
        integer,intent(in) :: natoms_pc
        integer,intent(in) :: natoms_sc
        logical,optional,intent(in) :: make_hermitian

        self%make_hermitian=.true.
        if (present(make_hermitian)) self%make_hermitian=make_hermitian
        self%has_lr=.false.
        self%has_lr_active=.false.
        self%det_p2s = det_p2s
        self%natoms_pc = natoms_pc
        self%natoms_sc = natoms_sc
        allocate(self%lattice_point(3,self%det_p2s))
        allocate(self%centers(3,self%natoms_pc))
        allocate(self%primitive_positions(3,self%natoms_pc))
        allocate(self%atom_indices_in_derivatives(self%natoms_pc))
        allocate(self%supercell_positions(3,self%natoms_sc))
        allocate(self%born_eff_charges(3,3,self%natoms_pc))
        allocate(self%quadrupoles(3,3,3,self%natoms_pc))
        self%quadrupoles = 0.0_q
        allocate(self%phonon_masses(self%natoms_pc))
        allocate(self%force_constants(3,3,self%natoms_sc,self%natoms_sc))
        allocate(self%shortest_vectors(3,27,self%natoms_pc,self%natoms_sc))
        allocate(self%shortest_vector_multiplicities(self%natoms_pc,self%natoms_sc))
    end subroutine phonon_ifc_init

    !> Free memory of the ifc datastructure
    subroutine phonon_ifc_free(self)
        use ini, only: deregister_allocate
        type(phonon_ifc) :: self
        deallocate(self%lattice_point)
        deallocate(self%centers)
        deallocate(self%primitive_positions)
        deallocate(self%phonon_masses)
        deallocate(self%force_constants)

        if (allocated(self%force_constants_lr)) then
            deallocate(self%force_constants_lr)
        endif

        call deregister_allocate(VMEM_ARR_SIZE(self%sp_lattice_atoms)+&
                                 VMEM_ARR_SIZE(self%sp_lattice_atoms_size), 'ifc_lattice_atoms')
        if (allocated(self%sp_lattice_atoms)) deallocate(self%sp_lattice_atoms)
        if (allocated(self%sp_lattice_atoms_size)) deallocate(self%sp_lattice_atoms_size)
        if (allocated(self%centers)) deallocate(self%centers)
        if (allocated(self%supercell_atoms)) deallocate(self%supercell_atoms)
        if (allocated(self%supercell_atoms_map_to_phelel)) deallocate(self%supercell_atoms_map_to_phelel)
    end subroutine phonon_ifc_free

    !> @brief read the potential information from hdf5 file
    subroutine phonon_ifc_from_hdf5(self,filename)
        use vhdf5
        use tutor, only: vtutor, isError, ElectronPhononRequiresHDF5
        type(phonon_ifc) :: self
        character(len=*),optional,intent(in) :: filename
#ifdef VASP_HDF5
        ! local variables
        integer :: i,ierr
        real(q) :: lattice_transpose(3,3)
        integer :: int_matrix(3, 3)
        integer(hid_t):: fileid
        character(len=:),allocatable :: my_filename

        my_filename = 'phelel_params.hdf5'; if (present(filename)) my_filename = filename
        ierr = vh5_file_open_read(my_filename, fileid)
        if (ierr/=0) call vtutor%error('Failed to open '//trim(my_filename)//'.  Please check if the file is present in the current path.')

        ! information regarding the grids
        ierr = vh5_read(fileid,'lattice_point',self%lattice_point)
        if (ierr/=0) then
            VH5_CHECK(vh5_read(fileid,'phonon_lattice_point',self%lattice_point))
        endif

        ! information regarding the lattice parameters
        ierr = vh5_read(fileid,'phonon_supercell_matrix',int_matrix)
        if (ierr/=0) then
            VH5_CHECK(vh5_read(fileid,'supercell_matrix',int_matrix))
        endif
        self%supercell_matrix = transpose(int_matrix)

        ierr = vh5_read(fileid,'phonon_supercell_lattice',lattice_transpose)
        if (ierr/=0) then
            VH5_CHECK(vh5_read(fileid,'supercell_lattice',lattice_transpose))
        endif
        self%lattice_sp = transpose(lattice_transpose)

        ierr = vh5_read(fileid,'phonon_supercell_positions',self%supercell_positions)
        if (ierr/=0) then
            VH5_CHECK(vh5_read(fileid,'supercell_positions',self%supercell_positions))
        endif

        ierr = vh5_read(fileid,'phonon_p2s_map',self%atom_indices_in_derivatives)
        if (ierr/=0) then
            ierr = vh5_read(fileid,'atom_indices_in_derivatives',self%atom_indices_in_derivatives)
        endif
        if (ierr/=0) then
            VH5_CHECK(vh5_read(fileid,'p2s_map',self%atom_indices_in_derivatives))
        endif

        ierr = vh5_read(fileid,'primitive_matrix',lattice_transpose)
        if (ierr/=0) then
           self%primitive_matrix = 0
           self%primitive_matrix(1,1) = 1;
           self%primitive_matrix(2,2) = 1;
           self%primitive_matrix(3,3) = 1;
        else
           self%primitive_matrix = transpose(lattice_transpose)
        endif

        VH5_CHECK(vh5_read(fileid,'primitive_lattice',lattice_transpose))
        self%lattice_pc = transpose(lattice_transpose)

        VH5_CHECK(vh5_read(fileid,'primitive_positions',self%primitive_positions))


        ! information regarding the force constants
        ! no transpose ?
        VH5_CHECK(vh5_read(fileid,'force_constants',self%force_constants))
        VH5_CHECK(vh5_read(fileid,'primitive_masses',self%phonon_masses))

        ! shortest_vectors: used only to chekc the structure
        ierr = vh5_read(fileid,'phonon_shortest_vectors',self%shortest_vectors)
        if (ierr/=0) deallocate(self%shortest_vectors)
        ierr = vh5_read(fileid,'phonon_shortest_vector_multiplicities',self%shortest_vector_multiplicities)
        if (ierr/=0) deallocate(self%shortest_vector_multiplicities)


        ! read the born effective charges
        ierr = vh5_read(fileid,'born_effective_charges',self%born_eff_charges)
        if (ierr/=0) self%born_eff_charges = 0
        do i=1,self%natoms_pc
           self%born_eff_charges(:,:,i)= transpose(self%born_eff_charges(:,:,i))
        enddo

        ! read the dielectric tensor
        ierr = vh5_read(fileid,'dielectric_constant',self%dielectric_tensor)
        if (ierr/=0) then
            self%dielectric_tensor = 0
            self%dielectric_tensor(1,1) = 1
            self%dielectric_tensor(2,2) = 1
            self%dielectric_tensor(3,3) = 1
        endif
        self%dielectric_tensor = transpose(self%dielectric_tensor)

        VH5_CHECK(vh5_file_close(fileid))
#else
        call vtutor%write(isError,ElectronPhononRequiresHDF5)
#endif
    end subroutine phonon_ifc_from_hdf5

    !> @brief setup the the interpolation of dynamical matrix
    subroutine phonon_ifc_setup(self,t_info)
        use ini, only: register_allocate
        use tutor
        use mathtools, only: matrix_inverse_3x3
        type(phonon_ifc) :: self
        type (type_info),intent(in)  :: t_info
        ! local variables
        integer :: i1, i2, is, iatom
        integer i,is_phelel
        real(q) :: scinv(3,3)
        real(q), allocatable :: shift_sp_atoms(:,:)

        PROFILING_START('phonon_ifc_setup')
        ! The matrix used to transform the coordinates from primitive/current to supercell
        scinv = matrix_inverse_3x3(real(self%supercell_matrix,q))
        self%prim2sp=matmul(scinv,self%primitive_matrix)

        ! The matrix used to transform the coordinates from supercell to primitive/current
        self%sp2prim = matrix_inverse_3x3(self%prim2sp)

        !> compute reciprocal lattice vectors, stored as columns, without 2pi
        self%latt_pc%a = self%lattice_pc
        call lattic(self%latt_pc)

        !> compute sp cell reciprocal lattice vectors, stored as columns, without 2pi
        self%latt_sp%a = self%lattice_sp
        call lattic(self%latt_sp)

        ! from here the array we use is SUPERCELL_ATOMS and not SUPERCELL_POSITIONS.
        ! We want to remove SUPERCELL_ATOMS_MAP_TO_PHELEL and SUPERCELL_ATOMS_MAP_FROM_PHELEL,
        ! and to do that we just want tp reorder phelel outputs
        !> Define SUPERCELL_ATOMS, where sp atoms are in [0,1[_SP and in 'convenient' order (like grid points)
        !>
        !> Suggestion: - remove SUPERCELL_POSITIONS from the code --> DONE in make_supercell_atoms
        !>             - we want to remove SUPERCELL_ATOMS_MAP_TO_PHELEL
        !>               and  SUPERCELL_ATOMS_MAP_FROM_PHELEL
        call make_supercell_atoms(self%prim2sp,self%lattice_point,self%primitive_positions,&
             self%supercell_positions,self%supercell_atoms,self%supercell_atoms_map_to_phelel,&
             self%supercell_atoms_map_from_phelel,shift_sp_atoms)

        !> create ityp array for supercell
        allocate(self%ityp_pc(self%natoms_pc))
        allocate(self%ityp_sc(self%natoms_sc))
        if (size(t_info%ityp)/=self%natoms_pc) then
           call vtutor%error('elph_ifc_setup: Inconsistent size of the array with the atomic types')
        endif
        self%ityp_pc = t_info%ityp
        is = 0
        do i=1,self%det_p2s
           do i1=1,self%natoms_pc
              is=is+1
              is_phelel = self%supercell_atoms_map_to_phelel(is)
              self%ityp_sc(is_phelel) = self%ityp_pc(i1)
          enddo
        enddo

        ! position of the atoms for which we get the phonon eigenvectors. Here
        ! we need to think more to be consistent with ELPHON
        do is=1,self%natoms_pc
            ! +1 because we start with 0 in python
            iatom = self%supercell_atoms_map_from_phelel(self%atom_indices_in_derivatives(is)+1)
            self%centers(:,is) = self%supercell_atoms(:,iatom)
        enddo

        ! prepare supercell atoms for summations
        call atoms_centered_grids(self%lattice_sp,self%natoms_pc,self%natoms_sc,self%centers,&
                                  self%supercell_atoms,shift_sp_atoms,self%max_points_atoms,self%sp_lattice_atoms,&
                                  self%sp_lattice_atoms_size)
        deallocate(shift_sp_atoms)
        call register_allocate(VMEM_ARR_SIZE(self%sp_lattice_atoms)+&
                               VMEM_ARR_SIZE(self%sp_lattice_atoms_size), 'ifc_lattice_atoms')

        ! convert derivative centers to primitive cell reduced coordinates
        self%centers=matmul(self%sp2prim,self%centers)

        do i1=1,self%natoms_pc
           do i2=1,size(self%sp_lattice_atoms,2)
              self%sp_lattice_atoms(1:3,i2,i1)=matmul(self%sp2prim,self%sp_lattice_atoms(1:3,i2,i1))
           enddo
        enddo
        PROFILING_STOP('phonon_ifc_setup')

    end subroutine phonon_ifc_setup

    !> @brief check that the phonon setup is the same than in phelel
    subroutine phonon_ifc_setup_check(self)
        use tutor
        use string, only: str
        use sym_prec, only: tiny
        type (phonon_ifc),intent(in) :: self
        ! local vectors
        logical :: found
        integer :: is,is_phelel,sizep,shift,j,is_phelel0,i0,i1,i2,k
        real(q) :: diff(3)
        real(q) :: vec(3),vec_phelel(3)
        ! check structure using shortest vectors
        is=0
        do i2=1,self%det_p2s
           do i1=1,self%natoms_pc
              is=is+1
              is_phelel = self%supercell_atoms_map_to_phelel(is)

              do i0=1,self%natoms_pc
                 is_phelel0=self%atom_indices_in_derivatives(i0)+1

                 sizep = self%sp_lattice_atoms_size(1,is,i0)
                 shift = self%sp_lattice_atoms_size(2,is,i0)

                 if (allocated(self%shortest_vector_multiplicities)) then
                    if (self%shortest_vector_multiplicities(i0,is_phelel) /= sizep) then
                       call vtutor%error('error in the multiplicities of shortest vectors')
                    endif
                 endif

                 do j=1,sizep

                    vec=self%primitive_positions(:,i1)+self%lattice_point(:,i2)+&
                        self%sp_lattice_atoms(:,shift+j,i0)-self%centers(:,i0)

                    if (allocated(self%shortest_vectors)) then
                       ! loop over shortest vectors in phelel and find vec
                       found=.false.
                       do k=1,sizep
                           vec_phelel=self%shortest_vectors(:,k,i0,is_phelel)

                           diff = vec(:)-vec_phelel(:)
                           if (all(abs(diff)<tiny)) found=.true.
                       enddo
                       ! check if vec was found
                       if (.not.found) then
                           call vtutor%error('Shortest vectors are different from phelel. '//&
                                             'The following vector was not found:'//str(vec))
                       endif
                    endif

                 enddo

              enddo !i0 natoms_pc

          enddo !i1 natoms_pc
        enddo !is det_p2s

    end subroutine phonon_ifc_setup_check

    !> @brief enforce acoustic sum rule
    !> the algorithm implemented here is the same one as in phonopy by Atsushi Togo
    subroutine phonon_ifc_asr(self,asr_iter)
        type (phonon_ifc) :: self
        integer,intent(in) :: asr_iter
        ! local variables
        integer :: iter
        integer :: iatom1, iatom2, idir1, idir2
        real(q) :: fstore(3,3)
        real(q),allocatable :: fctmp(:,:,:,:)

        allocate(fctmp(3,3,self%natoms_sc,self%natoms_sc))

        do iter=1,asr_iter
            ! subtract drift along row
            do iatom1=1,self%natoms_sc
                fstore = 0
                do iatom2=1,self%natoms_sc
                    fstore = fstore + self%force_constants(:,:,iatom1,iatom2)
                enddo
                do iatom2=1,self%natoms_sc
                    self%force_constants(:,:,iatom1,iatom2) = &
                    self%force_constants(:,:,iatom1,iatom2) - fstore/self%natoms_sc
                enddo
            enddo
            ! subtract drift along column
            do iatom1=1,self%natoms_sc
                fstore = 0
                do iatom2=1,self%natoms_sc
                    fstore = fstore + self%force_constants(:,:,iatom2,iatom1)
                enddo
                do iatom2=1,self%natoms_sc
                    self%force_constants(:,:,iatom2,iatom1) = &
                    self%force_constants(:,:,iatom2,iatom1) - fstore/self%natoms_sc
                enddo
            enddo
            ! average with the transpose
            fctmp = self%force_constants
            do iatom1=1,self%natoms_sc
            do iatom2=1,self%natoms_sc
                do idir1=1,3
                do idir2=1,3
                    self%force_constants(idir1,idir2,iatom2,iatom1) = &
                    self%force_constants(idir1,idir2,iatom2,iatom1) + &
                    fctmp(idir2,idir1,iatom1,iatom2)
                enddo
                enddo
            enddo
            enddo
            self%force_constants = .5_q*self%force_constants
        enddo

        deallocate(fctmp)

        ! translational symmetry
        do iatom1=1,self%natoms_sc
            fstore=0
            do iatom2=1,self%natoms_sc
                if (iatom2==iatom1) cycle
                fstore = fstore + self%force_constants(:,:,iatom2,iatom1)
            enddo
            self%force_constants(:,:,iatom1,iatom1) = -0.5*(fstore+transpose(fstore))
        enddo
    end subroutine phonon_ifc_asr

    !> @brief Setup the computation fo the long-range force constants
    subroutine phonon_ifc_setup_lr(self,grid,encutlr,comm)
    !> force_constants_lr(alpha_1,alpha_2,phelel sp index1, phelel sp index2)
    !>        = sum(g) (4pi/V_S) e^2/(4pi epsilon_0)
    !>                 *(g*Z_atom1)_ap (g*Z_atom2)_a/(g.epsilon.g)
    !>                 * exp(-g.epsilon.g/(4 ewald_parameter^2))
    !>
    !>
    !>  set to zero if 1) exp(-g.epsilon.g/(4  ewald_parameter^2)) < exp(supression_threshold)
    !>                 2) g < 10^(-10)
    !>
        use constant
        use lattice
        use tutor
        type (phonon_ifc) :: self
        type (grid_3d) :: grid
        real(q),intent(in) :: encutlr !< G-space cutoff for the long-range interaction
        type(communic),intent(in) :: comm
        ! Local variables
        real(q) :: gvec(3), gvec_car(3), vqpt(3)
        real(q) :: vec_car(3)
        integer :: ngvecs
        integer :: iatom
        integer :: ig, i1, i2
        integer :: is1, is2
        integer :: alphap, alpha, iatom1, iatom2
        integer :: is_phelel, i0_phelel
        real(q) :: arg
        real(q) :: prefactor, dielectric_factor, exp_supression
        real(q) :: born_factor(3,self%natoms_pc)
        real(q) :: atoms_sc_cart(3,self%natoms_sc)
        complex(q) :: phase
        integer, allocatable :: nindpw(:)
        integer, allocatable :: gvecs(:,:)
        type(communic) :: comm_self
        type(grid_3d) :: grid_pc
        type(grid_3d) :: grid_sc
        integer :: NG(3,3)
        real(q) :: glattice_sp(3,3)
        integer :: fft_mesh(3)
        if ((.not.is_identity(3,self%dielectric_tensor)).or.&
            any(self%born_eff_charges/=0).or.&
            any(self%quadrupoles/=0)) self%has_lr = .true.
        if (.not.self%has_lr) return
        PROFILING_START('ifc_setup_lr')

        self%encutlr = encutlr
        ! setup ewald parameter
        self%ewald_param = get_ewald_parameter(self%dielectric_tensor,self%encutlr)

        ! cartesian components of the supercell atoms
        atoms_sc_cart=matmul(self%lattice_sp,self%supercell_atoms)

        ! bounds for looking for the sphere of g-vectors from the supercell
        NG=nint(transpose(self%sp2prim))
        fft_mesh(1) = abs(NG(1,1))*grid%ngx+abs(NG(1,2))*grid%ngy+abs(NG(1,3))*grid%ngz
        fft_mesh(2) = abs(NG(2,1))*grid%ngx+abs(NG(2,2))*grid%ngy+abs(NG(2,3))*grid%ngz
        fft_mesh(3) = abs(NG(3,1))*grid%ngx+abs(NG(3,2))*grid%ngy+abs(NG(3,3))*grid%ngz
        call fftchk(fft_mesh)

        ! create grid_sc
        CALLMPI(comm_init_dummy(comm_self))
        call grid3d_init(grid_sc,fft_mesh(1),fft_mesh(2),fft_mesh(3),&
                            comm_self,comm_self,comm_self,real2cplx=.false.,lreal=.false.,lparallel=.false.)

        ! get g-sphere
        vqpt = 0._q
        call get_gsphere(grid_sc,self%latt_sp%b,encutlr,vqpt,ngvecs,gvecs,nindpw)
        glattice_sp=self%latt_sp%b

        ! allocate force constants
        allocate(self%force_constants_lr(3,3,self%natoms_sc,self%natoms_sc))

        prefactor = 4._q * pi /self%latt_sp%omega * felect

        ! Loop over G vectors
        self%ngvecs_sc = 0
        self%force_constants_lr=0.0_q
        do ig=1,ngvecs
           gvec = gvecs(:,ig)
#ifdef MPI
           if (mod(ig-1,comm%ncpu)/=comm%node_me-1) cycle
#endif
           ! Only take non-zero G.
           if (all(abs(gvec) < 1e-10_q)) cycle

           gvec_car = tpi*matmul(glattice_sp,gvec)

           dielectric_factor = 1._q / dot_product(gvec_car,matmul(self%dielectric_tensor,gvec_car))

           ! Beyond a certain lower limit there is no need to account for these anymore
           arg = - 1._q / (dielectric_factor * 4._q * self%ewald_param * self%ewald_param)
           if (arg<supression_threshold) cycle
           self%ngvecs_sc = self%ngvecs_sc + 1
           exp_supression = exp(arg)

           ! Compute born factors
           do iatom=1,self%natoms_pc
              born_factor(:,iatom) = matmul(gvec_car,self%born_eff_charges(:,:,iatom))
           enddo

           ! Compute the long-range force constants
           is1=0
           do i1=1,self%det_p2s
              do iatom1=1,self%natoms_pc
                 is1=is1+1
                 is_phelel = self%supercell_atoms_map_to_phelel(is1)

                 is2=0
                 do i2=1,self%det_p2s
                    do iatom2=1,self%natoms_pc
                       is2=is2+1
                       i0_phelel = self%supercell_atoms_map_to_phelel(is2)

                       vec_car=atoms_sc_cart(:,is2)-atoms_sc_cart(:,is1)
                       arg=dot_product(gvec_car,vec_car)
                       phase=exp(-ci*arg)

                       do alpha=1,3
                          do alphap=1,3

                          self%force_constants_lr(alpha,alphap,is_phelel,i0_phelel) = &
                          self%force_constants_lr(alpha,alphap,is_phelel,i0_phelel) + &
                          real(born_factor(alpha,iatom1)*born_factor(alphap,iatom2) * &
                              prefactor*dielectric_factor*exp_supression*phase)

                          enddo
                       enddo


                    enddo ! iatom2
                 enddo ! det_p2s
              enddo ! iatom1
           enddo ! det_p2s

        enddo ! loop over g vectors
        CALLMPI( M_sum_d(comm, self%force_constants_lr, size(self%force_constants_lr)) )
        CALLMPI( M_sum_i(comm, self%ngvecs_sc, 1))

        !translational symmetry
        !do iatom1=1,self%natoms_sc
        !    fstore=0
        !    do iatom2=1,self%natoms_sc
        !       if (iatom2==iatom1) cycle
        !       fstore = fstore + self%force_constants_lr(:,:,iatom2,iatom1)
        !    enddo
        !    self%force_constants_lr(:,:,iatom1,iatom1) = -0.5*(fstore+transpose(fstore))
        !enddo

        ! prepare a list of g-vectors for computation of LR part in the primitive cell
        vqpt = 0._q
        ! ASK_HM : I think one has too, otherwise the dynamical matrix may not be real

        ! create grid_pc
        CALLMPI(comm_init_dummy(comm_self))
        call grid3d_init(grid_pc,grid%ngx,grid%ngy,grid%ngz,&
                         comm_self,comm_self,comm_self,real2cplx=.false.,lreal=.false.,lparallel=.false.)
        call get_gsphere(grid_pc,self%latt_pc%b,2*self%encutlr,vqpt,self%ngvecs_pc,self%gvecs,self%nindpw)

        PROFILING_STOP('ifc_setup_lr')

    end subroutine phonon_ifc_setup_lr

    !> @brief Set option for the treatment of the long-range interaction
    !> This routine modifies the has_lr_active in the phonon_ifc structure
    subroutine phonon_ifc_lr_set(self,lr_activate)
        type (phonon_ifc) :: self
        logical,intent(in) :: lr_activate !< this indicates the new value that has_lr_active should take
        ! If we don't have the long-range part information there is nothing to do
        if (.not.self%has_lr) return
        ! first restore V_t in self%dvdu
        if (self%has_lr_active) then
            self%force_constants = self%force_constants + self%force_constants_lr
        endif
        ! now compute V_sr
        self%has_lr_active = lr_activate
        if (self%has_lr_active) then
            self%force_constants = self%force_constants - self%force_constants_lr
        endif
    end subroutine phonon_ifc_lr_set

    !> @brief Compute the dynamical matrix by interpolating the force constants
    subroutine phonon_ifc_get_dynmat_sr(self,vq,dynmat)
        use constant
        type(phonon_ifc) :: self
        real(q), intent(in) :: vq(3)
        complex(q) ::  dynmat(3*self%natoms_pc,3*self%natoms_pc)
        ! local variables
        real(q) vec(3),angle
        complex(q) phase,weight
        integer i0,i1,i2,is,ind0,ind1
        integer j,alpha,alphap
        integer :: is_phelel, i0_phelel
        integer :: shift, sizep

        PROFILING_START('elph_dynmat_sr')
        is=0
        do i2=1,self%det_p2s
           do i1=1,self%natoms_pc
              is=is+1
              is_phelel = self%supercell_atoms_map_to_phelel(is)

              do i0=1,self%natoms_pc
                 i0_phelel = self%atom_indices_in_derivatives(i0)+1
                 weight=cmplx(0.0_q,0.0_q,q)
                 sizep = self%sp_lattice_atoms_size(1,is,i0)
                 shift = self%sp_lattice_atoms_size(2,is,i0)
                 do j=1,sizep
                    ! C convention: exp(i q . (R_{l' \tau'} - R_{\tau}) )
                    vec=self%primitive_positions(:,i1)+self%lattice_point(:,i2)+&
                        self%sp_lattice_atoms(:,shift+j,i0)-self%centers(:,i0)
                    angle=dot_product(vq,vec)
                    phase=exp(citpi*angle)
                    weight=weight+phase
                 enddo
                 weight=weight/sizep
                 do alpha =1,3
                    ind0=3*(i0-1)+alpha
                    do alphap=1,3
                       ind1=3*(i1-1)+alphap
                       dynmat(ind0,ind1)=dynmat(ind0,ind1)+&
                                         self%force_constants(alphap,alpha,is_phelel,i0_phelel) / &
                                         sqrt(self%phonon_masses(i0)*self%phonon_masses(i1))*weight

                    enddo
                 enddo

              enddo

           enddo
        enddo

        ! unitary transformation to get D convention exp (i q R_{l'})
        do i1=1,self%natoms_pc
           do i0=1,self%natoms_pc
              vec=self%primitive_positions(:,i0)-self%primitive_positions(:,i1)
              angle=dot_product(vq,vec)
              phase=exp(citpi*angle)

              do alpha =1,3
                 ind0=3*(i0-1)+alpha
                 do alphap=1,3
                    ind1=3*(i1-1)+alphap
                    dynmat(ind0,ind1)=dynmat(ind0,ind1)*phase
                 enddo
              enddo

           enddo
        enddo

        ! make hermitian
        if (self%make_hermitian) dynmat = .5_q*(dynmat + conjg(transpose(dynmat)))
        PROFILING_STOP('elph_dynmat_sr')
    end subroutine phonon_ifc_get_dynmat_sr

    !> @brief Compute long-range dynamical matrix in the D convention
    subroutine phonon_ifc_get_dynmat_lr(self,vq,dynmat)
        use constant
        type(phonon_ifc) :: self
        real(q), intent(in) :: vq(3)
        complex(q),intent(inout) :: dynmat(:,:)
        ! Local variables
        integer :: ig
        integer :: idx1, idx2
        integer :: iatom1, iatom2, idir1, idir2
        integer :: counter
        real(q) :: arg
        real(q) :: prefactor, dielectric_factor, exp_supression
        complex(q) :: weight
        real(q) :: vec(3), qpgvec(3), qpgvec_car(3)
        real(q) :: born_factor(3,self%natoms_pc)

        PROFILING_START('elph_dynmat_lr')
        prefactor = 4._q * pi / self%latt_pc%omega * felect

        dynmat = cmplx(0.0_q,0.0_q,q)
        counter = 0
        do ig=1,self%ngvecs_pc

           qpgvec = vq+self%gvecs(:,ig)
           qpgvec_car = tpi*matmul(self%latt_pc%b,qpgvec)

           ! Only take non-zero Q+G
           if (all(abs(qpgvec) < 1e-10_q)) cycle

           ! Pre-compute stuff.
           dielectric_factor = 1._q / dot_product(qpgvec_car,matmul(self%dielectric_tensor,qpgvec_car))
           arg = - 1._q / (dielectric_factor * 4._q * self%ewald_param * self%ewald_param)
           if (arg<supression_threshold) cycle
           counter = counter + 1
           exp_supression = exp(arg)

           do iatom1 = 1, self%natoms_pc
              born_factor(:, iatom1) = matmul(qpgvec_car,self%born_eff_charges(:, :, iatom1))
           enddo

           do iatom1=1,self%natoms_pc
              do iatom2=1,self%natoms_pc
                 vec=self%primitive_positions(:,iatom2)-self%primitive_positions(:,iatom1)
                 arg = dot_product(qpgvec,vec)
                 weight = exp(citpi*arg)
                 do idir1=1,3
                     idx1 = (iatom1-1)*3+idir1
                     do idir2=1,3
                         idx2 = (iatom2-1)*3+idir2
                         dynmat(idx2,idx1) = dynmat(idx2,idx1)+&
                               prefactor*dielectric_factor*exp_supression*&
                               born_factor(idir2,iatom2)*born_factor(idir1,iatom1)/&
                               sqrt(self%phonon_masses(iatom1)*self%phonon_masses(iatom2))*weight
                     enddo !idir2
                 enddo !idir1
              enddo !iatom2
           enddo !iatom1

        enddo ! loop over g


        ! if needed add gamma contribution
        if (all(abs(vq) < 1e-10_q)) then
           qpgvec = [1._q,0._q,0._q]
           qpgvec_car = tpi*matmul(self%latt_pc%b,qpgvec)

           do iatom1 = 1, self%natoms_pc
              born_factor(:, iatom1) = matmul(qpgvec_car,self%born_eff_charges(:, :, iatom1))
           enddo

           dielectric_factor = 1._q / dot_product(qpgvec_car,matmul(self%dielectric_tensor,qpgvec_car))
           do iatom1=1,self%natoms_pc
              do iatom2=1,self%natoms_pc
                 do idir1=1,3
                     idx1 = (iatom1-1)*3+idir1
                     do idir2=1,3
                         idx2 = (iatom2-1)*3+idir2
                         dynmat(idx2,idx1) = dynmat(idx2,idx1)+&
                               prefactor*dielectric_factor*&
                               born_factor(idir2,iatom2)*born_factor(idir1,iatom1)/&
                               sqrt(self%phonon_masses(iatom1)*self%phonon_masses(iatom2))
                     enddo
                  enddo
              enddo
           enddo
        endif

        ! make hermitian
        if (self%make_hermitian) dynmat = .5_q*(dynmat + conjg(transpose(dynmat)))
        PROFILING_STOP('elph_dynmat_lr')
    end subroutine phonon_ifc_get_dynmat_lr

    !> @brief Compute phonon dynamical matrix and diagonalize
    subroutine phonon_ifc_get_phonons(self,vq,w_q,e_q)
        use constant
        type(phonon_ifc) :: self
        real(q), intent(in) :: vq(3)
        real(q), intent(inout) ::  w_q(:)
        complex(q), intent(inout) ::  e_q(:,:)
        ! local variables
        complex(q) ::  dynmat_lr(3*self%natoms_pc,3*self%natoms_pc)
        complex(q), allocatable :: work(:)
        integer lwork,info
        real(q), allocatable :: rwork(:)
        integer i0
        real(q), parameter :: ang= 1.0e-10_q    ! Ang
        real(q), parameter :: fac_to_omega_thz= 1.0e-12_q*sqrt(evtoj/(amtokg*ang*ang))
        PROFILING_START('ifc_phonons')
        e_q=cmplx(0.0_q,0.0_q,q)

        ! compute sr part
        call phonon_ifc_get_dynmat_sr(self,vq,e_q)

        ! compute lr part
        if (self%has_lr_active) then
            ! Coulomb dynamical matrix: dynmat_lr
            call phonon_ifc_get_dynmat_lr(self,vq,dynmat_lr)
            e_q=e_q+dynmat_lr
        endif

        ! diagonalize matrix
        lwork=6*self%natoms_pc-1
        allocate(work(1:lwork))
        allocate(rwork(1:9*self%natoms_pc-2))
        call zheev('V', 'U', 3 * self%natoms_pc, e_q, 3 * self%natoms_pc, w_q, work, lwork, rwork, info)
        deallocate(work)
        deallocate(rwork)

        ! calculate sqrt of eigenvalues and convert unit
        do i0=1,3*self%natoms_pc
           if (w_q(i0) .gt. 0.0_q) then
               w_q(i0)= sqrt( w_q(i0))*fac_to_omega_thz
           else
               w_q(i0)=-sqrt(-w_q(i0))*fac_to_omega_thz
           endif
        enddo
        PROFILING_STOP('ifc_phonons')
    end subroutine phonon_ifc_get_phonons

    subroutine phonon_ifc_get_phonons_zero_disp(self,vq,w_q,e_q)
        use constant
        type(phonon_ifc) :: self
        real(q), intent(in) :: vq(3)
        real(q), intent(inout) ::  w_q(:)
        complex(q), intent(inout) ::  e_q(:,:)
        ! local variables
        integer :: ierr
        call phonon_ifc_get_phonons(self,vq,w_q,e_q)
        call get_phonons_zero_disp(self%natoms_pc,self%phonon_masses,w_q,e_q,ierr)
    end subroutine phonon_ifc_get_phonons_zero_disp

    !> @brief write information about the long-range part of the interatomic force constants
    subroutine phonon_ifc_write_lr(self,iu0,iu6)
        type(phonon_ifc),intent(in) :: self
        integer,intent(in) :: iu0
        integer,intent(in) :: iu6
        if (iu0>=0) then
            write(iu0,*) '========== interatomic force constants  ============'
            write(iu0,*) 'treatment of LR:', self%has_lr_active
            if (self%has_lr_active) then
                write(iu0,*) 'encutlr:', self%encutlr
                write(iu0,*) 'ewald parameter:', self%ewald_param
                write(iu0,*) 'number of g-vectors to evaluate LR part in the supercell:', self%ngvecs_sc
                write(iu0,*) 'number of g-vectors to evaluate LR part in the unit cell:', self%ngvecs_pc
            endif
            write(iu0,*) '===================================================='
        endif
        if (iu6>=0) then
            if (self%has_lr_active) then
                call write_dielectric_born_quad(self%natoms_pc,self%dielectric_tensor,&
                                                self%born_eff_charges,self%quadrupoles,iu6)
            endif
        endif
    end subroutine phonon_ifc_write_lr

    !> @brief Write stuctural information stored in the elphon_pot structure
    subroutine phonon_ifc_write_structural(self,iu0)
        type(phonon_ifc),intent(in) :: self
        integer,intent(in) :: iu0
        ! local variables
        integer :: is
        if (iu0<0) return
        write(iu0,*)
        write(iu0,*) "Structural information"
        write(iu0,*) "--------------------------"
        write(iu0,*)
        write(iu0,*) "Supercell lattice vectors (SP)"
        write(iu0,'(3f14.8,A)') self%lattice_sp(:,1), ' (A1)'
        write(iu0,'(3f14.8,A)') self%lattice_sp(:,2), ' (A2)'
        write(iu0,'(3f14.8,A)') self%lattice_sp(:,3), ' (A3)'
        write(iu0,*)
        write(iu0,*) "Supercell atoms, in reduced coordinates of the supercell"
        do is=1,self%natoms_sc
            write(iu0,'(I4,3f14.8,I3)') is, self%supercell_atoms(1:3,is), self%ityp_sc(is)
        enddo
        write(iu0,*)
        write(iu0,*) "Primitive lattice vectors (PRIM) = POSCAR vectors"
        write(iu0,'(3f14.8,A)') self%lattice_pc(:,1), ' (a1)'
        write(iu0,'(3f14.8,A)') self%lattice_pc(:,2), ' (a2)'
        write(iu0,'(3f14.8,A)') self%lattice_pc(:,3), ' (a3)'
        write(iu0,*)
        write(iu0,*) "Derivative centers, in reduced coordinates of the PRIMITIVE CELL = POSCAR atoms [l]"
        do is=1,self%natoms_pc
           write(iu0,'(I4,3f14.8,I3)') is, self%centers(1:3,is), self%ityp_pc(is)
        enddo
        write(iu0,*)
        write(iu0,*) "Coordinates transformation matrix M from primitive cell (PRIM) to supercell (SP)"
        write(iu0,*) "  (x1_sp,x2_sp,x3_sp) = M (x1_prim,x2_prim,x3_prim)"
        write(iu0,*) "           (A1 A2 A3) = (a1 a2 a3) M^{-1}"
        write(iu0,'(3f14.8)') self%prim2sp(1,:)
        write(iu0,'(3f14.8)') self%prim2sp(2,:)
        write(iu0,'(3f14.8)') self%prim2sp(3,:)
        write(iu0,*)
        write(iu0,*) "Coordinates transformation matrix M^{-1} from the supercell (SP) to primitive cell (PRIM)"
        write(iu0,*) "  (x1_prim,x2_prim,x3_prim) = M^{-1} (x1_sp,x2_sp,x3_sp)"
        write(iu0,*) "                 (a1 a2 a3) = (A1 A2 A3) M"
        write(iu0,'(3f14.8)') self%sp2prim(1,:)
        write(iu0,'(3f14.8)') self%sp2prim(2,:)
        write(iu0,'(3f14.8)') self%sp2prim(3,:)
        write(iu0,*)
        write(iu0,*) "  Notation:"
        write(iu0,*) "  (x1,x2,x3) = (x1  x2  x3)^T   (column vector)"
        write(iu0,*) "  (a1 a2 a3) = (a1_1 a2_1 a3_1)"
        write(iu0,*) "               (a1_2 a2_2 a3_2)"
        write(iu0,*) "               (a1_3 a2_3 a3_3)"
        write(iu0,*)
    end subroutine phonon_ifc_write_structural

    !> @brief write the dielectric constant, born effective charges, and quadrupoles
    subroutine write_dielectric_born_quad(natom,dielectric_tensor,born_eff_charges,quadrupoles,iu)
        integer,intent(in) :: natom
        real(q),intent(in) :: dielectric_tensor(3,3)
        real(q),intent(in) :: born_eff_charges(3,3,natom)
        real(q),intent(in) :: quadrupoles(3,3,3,natom)
        integer,intent(in) :: iu
        ! Local variables
        integer :: iatom

        if (iu<=0) return

        if (.not.is_identity(3,dielectric_tensor)) then
            write(iu,*) 'Dielectric tensor'
            write(iu,'(3f10.5)') dielectric_tensor(:,1)
            write(iu,'(3f10.5)') dielectric_tensor(:,2)
            write(iu,'(3f10.5)') dielectric_tensor(:,3)
            write(iu,*)
        endif

        if (any(born_eff_charges/=0)) then
            write(iu,*) 'Born effective charges'
            do iatom=1,natom
                write(iu,*) 'atom', iatom
                write(iu,'(3f10.5)') born_eff_charges(:,1,iatom)
                write(iu,'(3f10.5)') born_eff_charges(:,2,iatom)
                write(iu,'(3f10.5)') born_eff_charges(:,3,iatom)
            end do
            write(iu,*)
        else
            write(iu,*) 'Born effective charges (not found)'
            write(iu,*)
        endif

        if (any(quadrupoles/=0)) then
            write(iu,*) 'Quadrupole tensor'
            do iatom=1,natom
                write(iu,*) 'atom', iatom
                write(iu,*) 'dir 1'
                write(iu,'(3f10.5)') quadrupoles(:,1,1,iatom)
                write(iu,'(3f10.5)') quadrupoles(:,1,2,iatom)
                write(iu,'(3f10.5)') quadrupoles(:,1,3,iatom)
                write(iu,*) 'dir 2'
                write(iu,'(3f10.5)') quadrupoles(:,2,1,iatom)
                write(iu,'(3f10.5)') quadrupoles(:,2,2,iatom)
                write(iu,'(3f10.5)') quadrupoles(:,2,3,iatom)
                write(iu,*) 'dir 3'
                write(iu,'(3f10.5)') quadrupoles(:,3,1,iatom)
                write(iu,'(3f10.5)') quadrupoles(:,3,2,iatom)
                write(iu,'(3f10.5)') quadrupoles(:,3,3,iatom)
            enddo
            write(iu,*)
        else
            write(iu,*) 'Quadrupole tensor (not found)'
            write(iu,*)
        endif
    end subroutine write_dielectric_born_quad

    !> @brief Compute the phonon dispersion on a list of k-points and write to file
    subroutine phonon_ifc_dispersion(self,nkpts,vkpt)
        use vhdf5
        USE constant
        type(phonon_ifc) :: self
        integer,intent(in) :: nkpts
        real(q),intent(in) :: vkpt(3,nkpts)
        ! Local variables
        integer :: ikpt
        real(q) :: eigenvalues(self%natoms_pc*3,nkpts)
        complex(q) :: eigenvectors(self%natoms_pc*3,self%natoms_pc*3,nkpts)
        real(q) :: kvec(3,nkpts), length(nkpts)

        ! k vectors in cartesians components
        kvec=matmul(self%latt_pc%b,vkpt)*tpi

        ! Compute phonon dispersion
        do ikpt=1,nkpts
            if (ikpt .eq. 1) then
               length(ikpt)=0.0_q
            else
               length(ikpt)=length(ikpt-1)+norm2(kvec(:,ikpt)-kvec(:,ikpt-1)) ! cumulated length along the path
            endif
            call phonon_ifc_get_phonons(self,vkpt(:,ikpt),eigenvalues(:,ikpt),eigenvectors(:,:,ikpt))
        enddo
        ! Write to hdf5 file
#ifdef VASP_HDF5
        call vh5_write_phonons_simple(ih5outfileid,vkpt,eigenvalues,eigenvectors,length)
#endif
    end subroutine phonon_ifc_dispersion

    !> @brief Initialize the electron-phonon potential structure
    subroutine elphon_pot_init(self,grid_pc,ngridpoints,npr_local,natoms_pc,natoms_pc_local,natoms_sc,ncdij,det_p2s,lmdim,fft_mesh,&
                               useshmem,writepot,metagga,comm_shmem)
        use ini, only: register_allocate
        use tutor, only: vtutor
        type(elphon_pot),intent(inout) :: self
        type(grid_3d),target,intent(in) :: grid_pc !< primtive lattice grid
        integer,intent(in) :: ngridpoints !< number of grids points in the supercell
        integer,intent(in) :: npr_local !< local number of grid points in the primitive cell (depends on the value of ncore)
        integer,intent(in) :: natoms_pc
        integer,intent(in) :: natoms_pc_local
        integer,intent(in) :: natoms_sc
        integer,intent(in) :: ncdij
        integer,intent(in) :: det_p2s
        integer,intent(in) :: lmdim
        integer,intent(in) :: fft_mesh(3)
        logical,intent(in) :: useshmem
        logical,intent(in) :: writepot !< write the electron-phonon potential to file for visualization
        logical,intent(in) :: metagga !< wether we are using a metagga functional or not
        type(communic),pointer :: comm_shmem
        ! local variables
        integer :: istat

        self%skip=.false.
        self%has_lr=.false.
        self%has_lr_active=.false.
        self%metagga=metagga
        ! set dimensions
        self%ngridpoints = ngridpoints
        self%ngridpoints_local = det_p2s*npr_local
        self%natoms_pc = natoms_pc
        self%natoms_pc_local = natoms_pc_local
        self%natoms_sc = natoms_sc
        self%ncdij = ncdij
        self%det_p2s = det_p2s
        self%lmdim = lmdim
        self%fft_mesh = fft_mesh
        self%npw = product(fft_mesh)
        self%npr_local = npr_local
        ! TODO: this should be properly initilized
        self%grid_pc = grid_pc

        ! structural data
        allocate(self%lattice_point(3,self%det_p2s))
        allocate(self%primitive_positions(3,self%natoms_pc))
        allocate(self%centers(3,self%natoms_pc))
        allocate(self%atom_indices_in_derivatives(self%natoms_pc))
        allocate(self%supercell_positions(3,self%natoms_sc))
        allocate(self%born_eff_charges(3,3,self%natoms_pc))
        allocate(self%quadrupoles(3,3,3,self%natoms_pc))
        self%quadrupoles = 0.0_q
        allocate(self%phonon_masses(self%natoms_pc))

        ! Allocate memory (on all cpus)
        self%has_shmem = .false.
#ifdef shmem_elphon
        self%has_shmem = useshmem
#endif
        if (self%has_shmem) then
            self%comm_shmem = comm_shmem
            ! poor man's way o checking for NCORE/=1 and shared memory activated
            if (self%ngridpoints/=self%ngridpoints_local.and.self%comm_shmem%ncpu/=1) then
                call vtutor%error('Shared memory for the electron-phonon driver NCSHMEM/=1 '//&
                                  'and NCORE/=1 is not supported.')
            endif
            istat=1
#ifdef shmem_elphon
            call M_shmem_alloc(self%comm_shmem,int(self%ngridpoints,kind=c_size_t),3_c_size_t,int(self%natoms_pc,kind=c_size_t),int(self%ncdij,kind=c_size_t),self%dvdu_shmem,self%dvdu,istat)
#endif
            if (istat/=0) then
                call vtutor%error('Failed allocation of shared memory segment for electron-phonon potential (dvdu).')
            endif
        else
            self%comm_shmem%node_me=1
            self%comm_shmem%ncpu=1
            allocate(self%dvdu(self%ngridpoints_local,3,self%natoms_pc,self%ncdij))
        endif
        if (self%metagga) allocate(self%dmudu(self%ngridpoints_local,3,self%natoms_pc,self%ncdij))
        if (writepot) allocate(self%dvdu_longrange(self%ngridpoints_local,3,self%natoms_pc))
        allocate(self%grid_point(3,self%ngridpoints_local))
        allocate(self%ccdij(self%lmdim,self%lmdim,self%natoms_pc,self%ncdij))
        allocate(self%ccqij(self%lmdim,self%lmdim,self%natoms_pc,self%ncdij))
        allocate(self%ddijdu(self%lmdim,self%lmdim,self%natoms_sc,3,self%natoms_pc_local,self%ncdij))
        allocate(self%dqijdu(self%lmdim,self%lmdim,self%natoms_sc,3,self%natoms_pc,self%ncdij))
        allocate(self%ccdrij(self%lmdim,self%lmdim,self%natoms_pc,3,self%ncdij))
        call register_allocate(VMEM_ARR_SIZE(self%grid_point), 'pot_grid')
        call register_allocate(VMEM_ARR_SIZE(self%dvdu)/self%comm_shmem%ncpu, 'pot_dvdu')
    end subroutine elphon_pot_init

    !> @brief read the potential information from phelel_params.hdf5
    !>    Dij like arrays are obtained following VASP convention, for example
    !>        ELPH_POT%ccdij(L',L,...)= <L|h|L'>-<L_ps|h_ps|L'_ps>
    !>
    !>   The lattice related arrays contains vectors columnwise, according to the rules:
    !>         lattice_sp = lattice_conv  * supercell_matrix
    !>         lattice_pc=  lattice_conv  * primitive_matrix
    !>
    !>  The Born effective charges are
    !>  born_eff_charges(alpha,beta,atom_pc)= (Omega / q) * d P_alpha / d R_(atom_pc, beta)
    !>

    subroutine elphon_pot_from_hdf5(self,comm_inb,comm_kbatch,filename)
        use vhdf5
        use tutor, only: vtutor, isError, ElectronPhononRequiresHDF5
        use wave, only: ni_local
        type(elphon_pot),intent(inout) :: self
        type(communic),target,intent(in) :: comm_inb !< intra-band comunicator (responsible for ion paralelism for projections)
        type(communic),target,intent(in) :: comm_kbatch !< comunicator between k-points
        character(len=*),optional,intent(in) :: filename
#ifdef VASP_HDF5
        ! local variables
        integer :: shift_local, shift_global
        integer :: i,ispinor,iatom,iatom_local,idir,ierr
        integer :: nx,ny,nz,nc
        real(q) :: lattice_transpose(3,3)
        integer :: int_matrix(3, 3)
        real(q) :: lattice_conv(3,3)    !> conventional unit cell, not used
        integer(hid_t) :: fileid
        integer(hid_t) :: start(6), count(6)
        real(q),allocatable :: grid_point_tmp(:,:)
        complex(q),allocatable :: dvdu_tmp(:)
        complex(q),allocatable :: ddijdu_tmp(:,:,:,:)
        character(len=:),allocatable :: my_filename

        my_filename = 'phelel_params.hdf5'; if (present(filename)) my_filename = filename
        ierr = vh5_file_open_read(my_filename, fileid)
        if (ierr/=0) call vtutor%error('Failed to open '//trim(my_filename)//'.  Please check if the file is present in the current path.')

        associate(grid_pc => self%grid_pc)
        ! information regarding the grids
        if (grid_pc%comm%ncpu==1) then
            VH5_CHECK(vh5_read(fileid,'grid_point',self%grid_point))
        else
            allocate(grid_point_tmp(3,self%ngridpoints))
            VH5_CHECK(vh5_read(fileid,'grid_point',grid_point_tmp))
            do i=1,self%det_p2s
               shift_local = (i-1)*self%npr_local
               shift_global = (i-1)*self%npw
               if (grid_pc%rl%nfast==3) then
                  do nc=1,grid_pc%rl%ncol
                     nx = grid_pc%rl%i2(nc)
                     ny = grid_pc%rl%i3(nc)
                     do nz=1,grid_pc%rl%nrow
                        self%grid_point(:,nz+(nc-1)*grid_pc%rl%nrow+shift_local)=&
                        grid_point_tmp(:,nx+(ny-1)*(grid_pc%ngx)+(nz-1)*(grid_pc%ngx*grid_pc%ngy)+shift_global)
                     enddo
                  enddo
               else
                  do nc=1,grid_pc%rl%ncol
                     ny = grid_pc%rl%i2(nc)
                     nz = grid_pc%rl%i3(nc)
                     do nx=1,grid_pc%rl%nrow
                        self%grid_point(:,nx+(nc-1)*grid_pc%rl%nrow+shift_local)=&
                        grid_point_tmp(:,nx+(ny-1)*(grid_pc%ngx)+(nz-1)*(grid_pc%ngx*grid_pc%ngy)+shift_global)
                     enddo
                  enddo
               endif
            enddo
            deallocate(grid_point_tmp)
        endif
        VH5_CHECK(vh5_read(fileid,'lattice_point',self%lattice_point))

        ! information regarding the local potential
        if (self%comm_shmem%node_me==1) then
        if (grid_pc%comm%ncpu==1) then
            ierr = vh5_read(fileid,'dVdu',self%dvdu)
            if (ierr/=0) call vtutor%error('Failed to read the potential dVdu from the phelel_params.hdf5 file.')
        else
            ierr = elphon_pot_from_hdf5_ncore(self,grid_pc,self%dvdu,fileid,'dVdu')
            if (ierr/=0) call vtutor%error('Failed to read the potential dVdu from the phelel_params.hdf5 file.')
        endif
        endif

        if (self%metagga) then
        if (self%comm_shmem%node_me==1) then
        if (grid_pc%comm%ncpu==1) then
            ierr = vh5_read(fileid,'dmudu',self%dmudu)
            if (ierr/=0) call vtutor%error('Failed to read the derivatives of the kinetic energy dmudu from the phelel_params.hdf5 file.')
        else
            ierr = elphon_pot_from_hdf5_ncore(self,grid_pc,self%dmudu,fileid,'dmudu')
            if (ierr/=0) call vtutor%error('Failed to read the derivatives of the kinetic energy dmudu from the phelel_params.hdf5 file.')
        endif
        endif
        endif
        call elphon_pot_shmem_barrier(self)
        end associate

        ! information regarding the projectors
        VH5_CHECK(vh5_read(fileid,'Dij',self%ccdij))
        VH5_CHECK(vh5_read(fileid,'qij',self%ccqij))
        self%comm_inb = comm_inb
        self%comm_kbatch = comm_kbatch
        if (comm_inb%ncpu==1) then
            VH5_CHECK(vh5_read(fileid,'dDijdu',self%ddijdu))
        else
            allocate(ddijdu_tmp(self%lmdim,self%lmdim,self%natoms_sc,3))
            do ispinor=1,self%ncdij
                do iatom=1,self%natoms_pc
                    iatom_local = ni_local(iatom,comm_inb)
                    if (iatom_local==0) cycle
                    start = [1,1,1,1,iatom,ispinor]
                    count = [self%lmdim,self%lmdim,self%natoms_sc,3,1,1]
                    VH5_CHECK(vh5_read_double_complex_subarray_nd(fileid,'dDijdu',6,start,count,ddijdu_tmp))
                    self%ddijdu(:,:,:,:,iatom_local,ispinor)=ddijdu_tmp(:,:,:,:)
                enddo
            enddo
            deallocate(ddijdu_tmp)
        endif
        VH5_CHECK(vh5_read(fileid,'dqijdu',self%dqijdu))

        ! information regarding the lattice parameters
        VH5_CHECK(vh5_read(fileid,'unitcell_lattice',lattice_transpose))
        lattice_conv= transpose(lattice_transpose)

        VH5_CHECK(vh5_read(fileid,'supercell_matrix',int_matrix))
        self%supercell_matrix = transpose(int_matrix)

        VH5_CHECK(vh5_read(fileid,'supercell_lattice',lattice_transpose))
        self%lattice_sp = transpose(lattice_transpose)

        VH5_CHECK(vh5_read(fileid,'supercell_positions',self%supercell_positions))

        ierr = vh5_read(fileid,'atom_indices_in_derivatives',self%atom_indices_in_derivatives)
        if (ierr/=0) then
           VH5_CHECK(vh5_read(fileid,'p2s_map',self%atom_indices_in_derivatives))
        endif

        ierr = vh5_read(fileid,'primitive_matrix',lattice_transpose)
        if (ierr/=0) then
           self%primitive_matrix = 0
           self%primitive_matrix(1,1) = 1;
           self%primitive_matrix(2,2) = 1;
           self%primitive_matrix(3,3) = 1;
        else
           self%primitive_matrix = transpose(lattice_transpose)
        endif

        VH5_CHECK(vh5_read(fileid,'primitive_lattice',lattice_transpose))
        self%lattice_pc = transpose(lattice_transpose)

        VH5_CHECK(vh5_read(fileid,'primitive_positions',self%primitive_positions))

        ! read the born effective charges
        ierr = vh5_read(fileid,'born_effective_charges',self%born_eff_charges)
        if (ierr/=0) self%born_eff_charges = 0
        do i=1,self%natoms_pc
           self%born_eff_charges(:,:,i)= transpose( self%born_eff_charges(:,:,i))
        enddo

        ! read the dielectric tensor
        ierr = vh5_read(fileid,'dielectric_constant',self%dielectric_tensor)
        if (ierr/=0) then
            self%dielectric_tensor = 0
            self%dielectric_tensor(1,1) = 1
            self%dielectric_tensor(2,2) = 1
            self%dielectric_tensor(3,3) = 1
        endif
        self%dielectric_tensor = transpose(self%dielectric_tensor)

        VH5_CHECK(vh5_file_close(fileid))
#else
        call vtutor%write(isError,ElectronPhononRequiresHDF5)
#endif
    end subroutine elphon_pot_from_hdf5

#ifdef VASP_HDF5
    !> Read potential from phelel_params.hdf5 with x fast and distribute it for the case
    !> of ncore/=1 to either x or z fast layout
    function elphon_pot_from_hdf5_ncore(self,grid,pot,fileid,groupname) result(ierr)
        use vhdf5
        use tutor, only: vtutor
        type(elphon_pot) :: self
        type(grid_3d),intent(in) :: grid
        complex(q),intent(inout) :: pot(:,:,:,:)
        integer(hid_t),intent(in) :: fileid
        character(len=*),intent(in) :: groupname
        integer :: ierr
        ! local variables
        integer :: shift_local, shift_global
        integer :: i,ispinor,iatom,idir,ierr_local
        integer :: nx,ny,nz,nc
        integer(HID_T) :: start(4), count(4)
        real(q),allocatable :: grid_point_tmp(:,:)
        complex(q),allocatable :: pot_tmp(:)

        ierr = 0
        allocate(pot_tmp(self%ngridpoints))
        do ispinor=1,self%ncdij
            do iatom=1,self%natoms_pc
                do idir=1,3
                    start = [1,idir,iatom,ispinor]
                    count = [self%ngridpoints,1,1,1]
                    ierr_local = vh5_read_double_complex_subarray_nd(fileid,groupname,4,start,count,pot_tmp)
                    ! return error code if reading failed
                    if (ierr_local/=0) then
                       ierr = ierr_local
                       return
                    endif
                    ! transfer full potential from phelel_params.h5 (with x fast) to dvdu (either x or z fast)
                    do i=1,self%det_p2s
                        shift_local = (i-1)*self%npr_local
                        shift_global = (i-1)*self%npw
                        if (grid%rl%nfast==3) then
                           do nc=1,grid%rl%ncol
                              nx = grid%rl%i2(nc)
                              ny = grid%rl%i3(nc)
                              do nz=1,grid%rl%nrow
                                 pot(nz+(nc-1)*grid%rl%nrow+shift_local,idir,iatom,ispinor)=&
                                 pot_tmp(nx+(ny-1)*(grid%ngx)+(nz-1)*(grid%ngx*grid%ngy)+shift_global)
                              enddo
                           enddo
                        else
                           do nc=1,grid%rl%ncol
                              ny = grid%rl%i2(nc)
                              nz = grid%rl%i3(nc)
                              do nx=1,grid%rl%nrow
                                 pot(nx+(nc-1)*grid%rl%nrow+shift_local,idir,iatom,ispinor)=&
                                 pot_tmp(nx+(ny-1)*(grid%ngx)+(nz-1)*(grid%ngx*grid%ngy)+shift_global)
                              enddo
                           enddo
                        endif
                    enddo
                enddo
            enddo
        enddo
        deallocate(pot_tmp)

    end function elphon_pot_from_hdf5_ncore
#endif

    !> @brief This subroutine defines several arrays used for the interpolation of the
    !>        potential using the minimal image convention:
    !>  1) ELPH_POT%supercell_atoms = ELPH_POT%primitive_positions + ELPH_POT%lattice_point + shift_sp_atoms
    !>                              are positions of atoms in the supercell [0,1[,
    !>                              in reduced coordinates of the sp
    !>
    !>     ELPH_POT%primitive_positions + ELPH_POT%lattice_point are supercell positions from phelel)
    !>
    !>  2) ELPH_POT%centers are the positions of atoms for which the derivatice of the potential are known,
    !>     in reduced coordinates of the sp
    !>
    !>  3) ELPH_POT%grid_point= r + ELPH_POT%lattice_points+ shift_sp
    !>      where grid_point is in[0,1[, in reduced coordinates wrt sp lattice vectors
    !>           r+ELPH_POT%lattice_points are the grid points from phelel,
    !>         and  r=(i/NGX,j/NGY,k/NGZ)
    !>
    !>  4) ELPH_POT%sp_lattice_points = ArgMin_L |L+ELPH_POT%grid_point-ELPH_POT%centers|+shift_sp
    !>                                + r+ELPH_POT%lattice_points
    !>     where L are supercell lattice vectors. (in reduced coordinates of the primitive cell)
    !>     Those vectors are used in minimum images calculations when the sum over
    !>     grid points is performed
    !>
    !>     ELPH_POT%sp_lattice_points(1:3,1:ELPH_POT%max_points_lattice,ELPH_POT%natoms_pc)
    !>     max_points_lattice = sum_(sp grid points)  nb nearest L
    !>
    !>  5) ELPH_POT%sp_lattice_atoms= ArgMin_L |L+ELPH_POT%supercell_atoms-ELPH_POT%centers|+shift_sp_atoms,
    !>     where L are supercell lattice vectors. (in reduce coordinate of the primitive cell)
    !>
    !>     ELPH_POT%sp_lattice_atoms(1:3,ELPH_POT%max_points_atoms,ELPH_POT%natoms_pc)
    !>     max_points_atoms = sum_(natoms_sc)  nb  nereast L
    !>
    !>    Those vectors are used in minimum images calculations when sum over
    !>    supercell atoms is performed
    !>
    !>  6) Suggestion: maybe add ELPH_POT%primitive_positions + ELPH_POT%lattice_point to
    !>              ELPH_POT%sp_lattice_atoms

    subroutine elphon_pot_setup(self,t_info,p)
        use ini, only: register_allocate
        use tutor
        use mathtools, only: matrix_inverse_3x3
        use sym_prec, only: tiny
        type(elphon_pot) :: self
        type (type_info),intent(in) :: t_info
        type (potcar),target,intent(in) :: p(:)
        ! local variables
        integer :: is, i1, i2, iatom
        integer i,is_phelel
#ifdef shmem_elphon
        integer :: istat
#endif
        real(q) :: diff(3)
        real(q) :: scinv(3,3)
        real(q), allocatable :: shift_sp(:,:)
        real(q), allocatable :: shift_sp_atoms(:,:)

        PROFILING_START('elph_pot_setup')
        ! The matrix used to transform the coordinates from primitive/current to supercell
        scinv = matrix_inverse_3x3(real(self%supercell_matrix,q))
        self%prim2sp=matmul(scinv,self%primitive_matrix)

        ! The matrix used to transform the coordinates from supercell to primitive/current
        self%sp2prim = matrix_inverse_3x3(self%prim2sp)

        !> compute reciprocal lattice vectors, stored as columns, without 2pi
        self%latt_pc%a = self%lattice_pc
        call lattic(self%latt_pc)

        !> compute sp cell reciprocal lattice vectors, stored as columns, without 2pi
        self%latt_sp%a = self%lattice_sp
        call lattic(self%latt_sp)

        ! from here the array we use is SUPERCELL_ATOMS and not SUPERCELL_POSITIONS.
        ! We want to remove SUPERCELL_ATOMS_MAP_TO_PHELEL and SUPERCELL_ATOMS_MAP_FROM_PHELEL,
        ! and to do that we just want tp reorder phelel outputs
        !> Define SUPERCELL_ATOMS, where sp atoms are in [0,1[_SP and in 'convenient' order (like grid points)
        !>
        !> Suggestion: - remove SUPERCELL_POSITIONS from the code
        !>             - we want to remove SUPERCELL_ATOMS_MAP_TO_PHELEL
        !>               and  SUPERCELL_ATOMS_MAP_FROM_PHELEL
        call make_supercell_atoms(self%prim2sp,self%lattice_point,self%primitive_positions,&
             self%supercell_positions,self%supercell_atoms,self%supercell_atoms_map_to_phelel,&
             self%supercell_atoms_map_from_phelel,shift_sp_atoms)

        !> create ityp array for supercell
        allocate(self%ityp_pc(self%natoms_pc))
        allocate(self%ityp_sc(self%natoms_sc))
        if (size(t_info%ityp)/=self%natoms_pc) then
           call vtutor%error('elph_pot_setup: Inconsistent size of the array with the atomic types')
        endif
        self%ityp_pc = t_info%ityp
        is = 0
        do i=1,self%det_p2s
           do i1=1,self%natoms_pc
              is=is+1
              is_phelel = self%supercell_atoms_map_to_phelel(is)
              self%ityp_sc(is_phelel) = self%ityp_pc(i1)
          enddo
        enddo

        !> store pointer to the pseudopotential
        self%p => p

        !>  The atoms in the supercell [0,1[_SP for which the potential derivatives are known
        !>  They are in general different from the atoms in POSCAR, but can be mapped to them by a
        !>  primitive lattice translation.
        !>
        !> Suggestion: Add +1 to atom_indices_in_derivatives at read in time
        do is=1,self%natoms_pc
            ! +1 because we start with 0 in python
            iatom = self%supercell_atoms_map_from_phelel(self%atom_indices_in_derivatives(is)+1)
            self%centers(:,is) = self%supercell_atoms(:,iatom)
        enddo

        !>  put grid points in [0,1[_SP according to : GRID_POINT=r+l+SHIFT_SP, with r and l from phelel
        allocate(shift_sp(3,self%ngridpoints_local))
        call grid_in_01(self%grid_pc,self%fft_mesh,self%prim2sp,self%det_p2s,self%ngridpoints_local,&
                        self%lattice_point,self%grid_point,shift_sp)

        !>  compute SP_LATTICE_POINTS=LPOINTS+SHIFT_SP, with LPOINTS the L vectors which
        !>  gives equiv GRID_POINTS wrt the CENTERS.
        call atoms_centered_grids_count(self%lattice_sp,self%natoms_pc,self%ngridpoints_local,self%centers,self%grid_point,self%max_points_lattice)

#ifdef shmem_elphon
        call M_shmem_alloc(self%comm_shmem,3_c_size_t,int(self%max_points_lattice,kind=c_size_t),int(self%natoms_pc,kind=c_size_t),self%sp_lattice_points_shmem,self%sp_lattice_points,istat)
        if (istat/=0) then
            call vtutor%error('Failed allocation of shared memory segment for electron-phonon potential (sp_lattice_points).')
        endif
#else
        allocate(self%sp_lattice_points(3,self%max_points_lattice,self%natoms_pc))
#endif
        allocate(self%sp_lattice_points_size(2,self%ngridpoints_local,self%natoms_pc))
        self%sp_lattice_points_size=0
        if (self%comm_shmem%node_me==1) then
        call atoms_centered_grids_generate(self%lattice_sp,self%natoms_pc,self%ngridpoints_local,self%centers,self%grid_point,shift_sp, &
                                           self%sp_lattice_points,self%sp_lattice_points_size)
        endif
#ifdef shmem_elphon
        CALLMPI( M_sum_i(self%comm_shmem, self%sp_lattice_points_size,size(self%sp_lattice_points_size)))
#endif
        deallocate(shift_sp)
        call register_allocate(VMEM_ARR_SIZE(self%sp_lattice_points)/self%comm_shmem%ncpu+&
                               VMEM_ARR_SIZE(self%sp_lattice_points_size)/self%comm_shmem%ncpu, 'pot_lattice_points')

        !> compute sp_lattice_atoms= LPOINTS+shift_sp_atoms, with LPOINTS the L
        !> vectors which gives equiv supercell_atoms wrt the CENTERS.
        call atoms_centered_grids(self%lattice_sp,self%natoms_pc,self%natoms_sc,self%centers,&
                                  self%supercell_atoms,shift_sp_atoms,self%max_points_atoms,self%sp_lattice_atoms,&
                                  self%sp_lattice_atoms_size)
        deallocate(shift_sp_atoms)

        !> convert derivative centers to primitive cell reduced coordinates
        self%centers=matmul(self%sp2prim,self%centers)

        !> convert supercell translations to primitive cell reduced coordinates
        if (self%comm_shmem%node_me==1) then
        do i1=1,self%natoms_pc
           do i2=1,size(self%sp_lattice_points,2)
              self%sp_lattice_points(1:3,i2,i1)=matmul(self%sp2prim,self%sp_lattice_points(1:3,i2,i1))
           enddo
        enddo
        endif

        !> add r + lattice_point to sp_lattice_points
        if (self%comm_shmem%node_me==1) then
        call prepare_grid(self%grid_pc,self%natoms_pc,self%det_p2s,self%npr_local,self%fft_mesh,self%lattice_point,&
                          self%sp_lattice_points,self%sp_lattice_points_size)
        endif
        call elphon_pot_shmem_barrier(self)

        !> convert sp_lattice_atoms to primitive cell reduced coordinates
        do i1=1,self%natoms_pc
           do i2=1,size(self%sp_lattice_atoms,2)
              self%sp_lattice_atoms(1:3,i2,i1)=matmul(self%sp2prim,self%sp_lattice_atoms(1:3,i2,i1))
           enddo
        enddo

        !> we check that center derivatives are the atoms in POSCAR, up to a lattice translation
        !>  PRIMITIVE_POSITIONS have been MADE to be EXACTLY POSCAR (t_info%posion in fact) atoms.
        outer: do i1=1,self%natoms_pc
           diff = self%centers(:,i1)-self%primitive_positions(:,i1)
           if (all(abs(diff-nint(diff)) < tiny)) cycle outer
           call vtutor%error('derivative centers are not primitive cell atoms')
        enddo outer

        PROFILING_STOP('elph_pot_setup')
    end subroutine elphon_pot_setup

    subroutine elphon_pot_shmem_barrier(self)
        type(elphon_pot),intent(inout) :: self
        ! local variables
        integer :: ierr
#ifdef shmem_elphon
        CALL M_barrier(self%comm_shmem)
#endif
    end subroutine elphon_pot_shmem_barrier

    !> check consistentcy between a set of arrays with the atomic positions and the vasp datastructures
    subroutine lattice_check(lattice_pc,natoms_pc,primitive_positions,latt_cur,t_info)
        use tutor
        use string, only: str
        use sym_prec, only: tiny
        real(q),intent(in) :: lattice_pc(3,3)
        integer,intent(in) :: natoms_pc
        real(q),intent(in) :: primitive_positions(:,:)
        type (latt),intent(in)       :: latt_cur
        type (type_info),intent(in)  :: t_info
        ! local variables
        integer :: ni
        real(q) :: diff(3)

        ! check the primitive cell
        do ni=1,3
           if (any(abs(lattice_pc(:,ni)-latt_cur%a(:,ni)) > tiny)) then
              write(*,*) abs(lattice_pc(:,ni)-latt_cur%a(:,ni)), tiny
              call vtutor%error('The lattice in phelel and POSCAR are different')
           endif
        enddo

        ! check dimensions of arrays
        if (size(primitive_positions,2) /= t_info%nions .or.&
            natoms_pc /= t_info%nions) then
           call vtutor%error('Number of atoms in POSCAR '//str(t_info%nions)//&
                             ' are different from the ones in the phelel primitive cell '//str(size(primitive_positions,2)))
        endif

        ! check the positions
        do ni=1,t_info%nions
           diff = t_info%posion(:,ni)-primitive_positions(:,ni)
           if (any(abs(diff-nint(diff))>tiny)) then
              call vtutor%error('positions in phelel and in POSCAR are different atom:'//str(ni)//&
                                ' poscar:'//str(t_info%posion(:,ni))//&
                                ' phelel:'//str(primitive_positions(:,ni)))
           endif
        enddo
    end subroutine lattice_check

    subroutine centers_check(centers, t_info)
        use tutor
        use string, only: str
        use sym_prec, only: tiny
        real(q),intent(in) :: centers(:,:)
        type (type_info),intent(in)  :: t_info
        ! local variables
        integer :: ni
        real(q) :: diff(3)
        ! check number of atoms
        if (size(centers,2)/=t_info%nions) then
           call vtutor%error('Number of atoms in POSCAR '//str(t_info%nions)//&
                             ' are different from the ones in the phelel primitive cell '//str(size(centers,2)))
        endif
        ! check that centers are primitve cell atoms, and given in the same order
        do ni=1,t_info%nions
           diff = t_info%posion(:,ni)-centers(:,ni)
           if (any(abs(diff-nint(diff))>tiny)) then
              call vtutor%error('Derivative centers are not equivalent to the primitive cell atoms')
           endif
        enddo
    end subroutine centers_check

    !> @brief Compute the Rij matrix
    !>  compute ELPH_POT%ccdrij(LMP,LM,NI,ISP,DIR) = R_{LM LMP}(NI)
    !>                       = < \phi_{LM}(NI) | d/dR_{NI} \phi_{LMP}(NI) >
    !>                      - < \tilde{\phi}_{LM}(NI) | d/dR_{NI}  \tilde{\phi}_{LMP}(NI) >
    !> NI are atoms in POSCAR
    subroutine elphon_pot_setup_rij(self,p,t_info)
        USE pseudo
        USE poscar, ONLY: type_info
        USE wap
        type (elphon_pot),intent(inout) :: self
        type (potcar) :: p(:)
        type (type_info) :: t_info


        !> Compute R0 and R1 in Eqs D3 and 4 of PRB 100, 174304 (2019)
        allocate(self%ccdrij_radial(self%lmdim,self%lmdim,self%natoms_pc,self%ncdij,2))
        call getwave_rij(self%ncdij,self%lmdim,p,t_info,self%ccdrij_radial)

        !> Compute Eq 73 PRB 100, 174304 (2019)
        call get_ccrij_a(self%ncdij,self%lmdim,t_info,self%ccdrij_radial,p,self%ccdrij)

        ! This should be avoided
        if (self%ncdij==4) then
           self%ccdrij(:,:,:,:,2)=0
           self%ccdrij(:,:,:,:,3)=0
        endif
    end subroutine elphon_pot_setup_rij

    !> transfer the potential from representation stored in file to the
    !> representation to be used in this calculation
    !> The representation is determined by the integer ncdij
    !> ncdij = 1 - spin unpolarized
    !>         2 - spin polarized
    !>         4 - noncollinear spin
    subroutine elphon_pot_transfer_ncdij(self,wf_ncdij,ccdij,ccqij,ccdrij)
        use tutor, only: vtutor
        use elphon_base, only: pot_spin_to_string
        type (elphon_pot),target,intent(inout) :: self
        integer, intent(in) :: wf_ncdij
        complex(q), allocatable, target, intent(out) :: ccdij(:,:,:,:)
        complex(q), allocatable, target, intent(out) :: ccqij(:,:,:,:)
        complex(q), allocatable, target, intent(out) :: ccdrij(:,:,:,:,:)
        ! local variables
        integer :: n
        complex(q), pointer :: a(:,:)
        complex(q), pointer :: b(:,:)
        complex(q), pointer, contiguous :: ccdij_ptr(:,:,:,:)
        complex(q), pointer, contiguous :: ccqij_ptr(:,:,:,:)
        complex(q), pointer, contiguous :: ccdrij_ptr(:,:,:,:,:)

        if (wf_ncdij/=self%ncdij) then
            call vtutor%alert('The electron-phonon potential is '//&
                       pot_spin_to_string(self%ncdij)//&
                       ' while the current calculation is '//&
                       pot_spin_to_string(wf_ncdij)//&
                       '. VASP will transfer the known information to the current calculation.')
        endif

        allocate(ccdij(self%lmdim,self%lmdim,self%natoms_pc,wf_ncdij))
        n=product(shape(self%ccdij))/self%ncdij
        if (n/=self%lmdim*self%lmdim*self%natoms_pc) call vtutor%bug('',__FILE__,__LINE__)
        ccdij_ptr => self%ccdij
        a(1:n,1:self%ncdij) => ccdij_ptr
        b(1:n,1:wf_ncdij) => ccdij
        call transfer_ncdij(a,self%ncdij,b,wf_ncdij)

        allocate(ccqij(self%lmdim,self%lmdim,self%natoms_pc,wf_ncdij))
        n=product(shape(self%ccqij))/self%ncdij
        if (n/=self%lmdim*self%lmdim*self%natoms_pc) call vtutor%bug('',__FILE__,__LINE__)
        ccqij_ptr => self%ccqij
        a(1:n,1:self%ncdij) => ccqij_ptr
        b(1:n,1:wf_ncdij) => ccqij
        call transfer_ncdij(a,self%ncdij,b,wf_ncdij)

        allocate(ccdrij(self%lmdim,self%lmdim,self%natoms_pc,3,wf_ncdij))
        n=product(shape(self%ccdrij))/self%ncdij
        if (n/=self%lmdim*self%lmdim*self%natoms_pc*3) call vtutor%bug('',__FILE__,__LINE__)

        ccdrij_ptr => self%ccdrij
        a(1:n,1:self%ncdij) => ccdrij_ptr
        b(1:n,1:wf_ncdij) => ccdrij
        call transfer_ncdij(a,self%ncdij,b,wf_ncdij)
    end subroutine elphon_pot_transfer_ncdij

    !> @brief Set the quadrupole matrices for diamond-like structures
    subroutine elphon_quadrupoles_set(quadrupoles,qval1,qval2)
        use string, only: str
        use tutor, only: vtutor
        real(q),intent(inout) :: quadrupoles(:,:,:,:)
        real(q),intent(in) :: qval1
        real(q),intent(in) :: qval2
        if (size(quadrupoles,4) /= 2) then
           call vtutor%error('Number of atoms '//str(size(quadrupoles,4))//' in the unit cell'//&
                             ' not allowed to set quadrupoles from input file. Must be 2.')
        endif
        !atom 1
        !dir 1
        quadrupoles(:,1,1,1) = [ 0._q, 0._q, 0._q]
        quadrupoles(:,2,1,1) = [ 0._q, 0._q,qval1]
        quadrupoles(:,3,1,1) = [ 0._q,qval1, 0._q]
        !dir 2
        quadrupoles(:,1,2,1) = [ 0._q, 0._q,qval1]
        quadrupoles(:,2,2,1) = [ 0._q, 0._q, 0._q]
        quadrupoles(:,3,2,1) = [qval1, 0._q, 0._q]
        !dir 3
        quadrupoles(:,1,3,1) = [ 0._q,qval1, 0._q]
        quadrupoles(:,2,3,1) = [qval1, 0._q, 0._q]
        quadrupoles(:,3,3,1) = [ 0._q, 0._q, 0._q]

        !atom 2
        !dir 1
        quadrupoles(:,1,1,2) = [ 0._q, 0._q, 0._q]
        quadrupoles(:,2,1,2) = [ 0._q, 0._q,qval2]
        quadrupoles(:,3,1,2) = [ 0._q,qval2, 0._q]
        !dir 2
        quadrupoles(:,1,2,2) = [ 0._q, 0._q,qval2]
        quadrupoles(:,2,2,2) = [ 0._q, 0._q, 0._q]
        quadrupoles(:,3,2,2) = [qval2, 0._q, 0._q]
        !dir 3
        quadrupoles(:,1,3,2) = [ 0._q,qval2, 0._q]
        quadrupoles(:,2,3,2) = [qval2, 0._q, 0._q]
        quadrupoles(:,3,3,2) = [ 0._q, 0._q, 0._q]
    end subroutine elphon_quadrupoles_set

    !> Write information about the long-range part of the potential
    subroutine elphon_pot_write_lr(self,iu)
        use elphon_base, only: pot_spin_to_string, lr_sft_string
        type(elphon_pot),intent(in) :: self
        integer,intent(in) :: iu
        if (iu<0) return
        write(iu,*) '=========== electron-phonon potential  ============='
        write(iu,*) 'treatment of LR:', self%has_lr_active
        write(iu,*) 'shmem:',self%has_shmem
        write(iu,*) pot_spin_to_string(self%ncdij)
        if (self%has_shmem) then
             write(iu,*) 'shmem ranks on master:', self%comm_shmem%ncpu
        endif
        if (self%has_lr_active) then
             write(iu,*) 'encutlr:', self%encutlr
             write(iu,*) 'ewald parameter:', self%ewald_param
             !write(iu,*) 'coulomb truncation radius:', self%rcut
             write(iu,*) 'use sft: ', lr_sft_string(self%use_lr_sft)
             write(iu,*) 'fft_mesh_sc: ', self%fft_mesh_sc
             write(iu,*) 'fft_mesh_scp:', self%fft_mesh_scp
             write(iu,*) 'number of g-vectors to evaluate LR part in the supercell:', self%ngvecs_sc
             write(iu,*) 'number of g-vectors to evaluate LR part in the unit cell:', self%ngvecs_pc
        endif
        write(iu,*) '===================================================='
    end subroutine elphon_pot_write_lr

    !> Write dielectric tensor, born effective charges and quadrupoles to output
    subroutine elphon_pot_write_dielectric_born_quad(self,iu)
        type(elphon_pot),intent(in) :: self
        integer,intent(in) :: iu
        if (iu<0) return
        if (.not.self%has_lr_active) return
        call write_dielectric_born_quad(self%natoms_pc,self%dielectric_tensor,&
                                        self%born_eff_charges,self%quadrupoles,iu)
    end subroutine elphon_pot_write_dielectric_born_quad

#ifdef VASP_HDF5
    !> write rij matrix to hdf5 file
    subroutine elphon_write_ccdrij(self,fileid)
        use vhdf5
        type(elphon_pot),intent(in) :: self
        integer(HID_T) :: fileid
        ! local variables
        integer(HID_T) :: groupid, subgroupid

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, 'results/'//trim(GRP_ELPHON), groupid))
        VH5_CHECK(vh5_group_open_or_create(groupid, 'electrons', subgroupid))

        VH5_CHECK(vh5_write(subgroupid, 'drij_radial', self%ccdrij_radial))
        VH5_CHECK(vh5_write(subgroupid, 'drij', self%ccdrij))
    end subroutine elphon_write_ccdrij
#endif

    !> Write stuctural information stored in the elphon_pot structure
    subroutine elphon_pot_write_structural(self,iu)
        type(elphon_pot),intent(in) :: self
        integer,intent(in) :: iu
        ! local variables
        integer :: is
        if (iu<0) return
        write(iu,*)
        write(iu,*) "Structural information"
        write(iu,*) "--------------------------"
        write(iu,*)
        write(iu,*) "Supercell lattice vectors (SP)"
        write(iu,'(3f14.8,A)') self%lattice_sp(:,1), ' (A1)'
        write(iu,'(3f14.8,A)') self%lattice_sp(:,2), ' (A2)'
        write(iu,'(3f14.8,A)') self%lattice_sp(:,3), ' (A3)'
        write(iu,*)
        write(iu,*) "Supercell atoms, in reduced coordinates of the supercell"
        do is=1,self%natoms_sc
            write(iu,'(3f14.8,I3)') self%supercell_atoms(1:3,is), self%ityp_sc(is)
        enddo
        write(iu,*)
        write(iu,*) "Primitive lattice vectors (PRIM) = POSCAR vectors"
        write(iu,'(3f14.8,A)') self%lattice_pc(:,1), ' (a1)'
        write(iu,'(3f14.8,A)') self%lattice_pc(:,2), ' (a2)'
        write(iu,'(3f14.8,A)') self%lattice_pc(:,3), ' (a3)'
        write(iu,*)
        write(iu,*) "Derivative centers, in reduced coordinates of the PRIMITIVE CELL = POSCAR atoms [l]"
        do is=1,self%natoms_pc
           write(iu,'(3f14.8,I3)') self%centers(1:3,is), self%ityp_pc(is)
        enddo
        write(iu,*)
        write(iu,*) "Coordinates transformation matrix M from primitive cell (PRIM) to supercell (SP)"
        write(iu,*) "  (x1_sp,x2_sp,x3_sp) = M (x1_prim,x2_prim,x3_prim)"
        write(iu,*) "           (A1 A2 A3) = (a1 a2 a3) M^{-1}"
        write(iu,'(3f14.8)') self%prim2sp(1,:)
        write(iu,'(3f14.8)') self%prim2sp(2,:)
        write(iu,'(3f14.8)') self%prim2sp(3,:)
        write(iu,*)
        write(iu,*) "Coordinates transformation matrix M^{-1} from the supercell (SP) to primitive cell (PRIM)"
        write(iu,*) "  (x1_prim,x2_prim,x3_prim) = M^{-1} (x1_sp,x2_sp,x3_sp)"
        write(iu,*) "                 (a1 a2 a3) = (A1 A2 A3) M"
        write(iu,'(3f14.8)') self%sp2prim(1,:)
        write(iu,'(3f14.8)') self%sp2prim(2,:)
        write(iu,'(3f14.8)') self%sp2prim(3,:)
        write(iu,*)
        write(iu,*) "  Notation:"
        write(iu,*) "  (x1,x2,x3) = (x1  x2  x3)^T   (column vector)"
        write(iu,*) "  (a1 a2 a3) = (a1_1 a2_1 a3_1)"
        write(iu,*) "               (a1_2 a2_2 a3_2)"
        write(iu,*) "               (a1_3 a2_3 a3_3)"
        write(iu,*)
    end subroutine elphon_pot_write_structural

    !> @brief Setup the computation of the long-range part of the potential
    subroutine elphon_pot_setup_lr(self,grid,encutlr,use_lr_sft,iu)
        use ini, only: register_allocate
        use constant
        use lattice
        use pseudo, only: potcar
        use mgrid_struct_def
        use string, only: str
        use mpimy, only: communic
        type (elphon_pot) :: self
        type (grid_3d),target :: grid
        real(q),intent(in) :: encutlr !< G-space cutoff for the long-range interaction
        integer,intent(in) :: use_lr_sft
        integer,intent(in) :: iu !< output unit to repor status of the calculation
        ! Local variables
        integer :: ngvecs
        integer :: idir,iatom
        real(q) :: sp_at(3,self%natoms_pc)
        real(q) :: vqpt(3)
        complex(q), allocatable :: v_long_range_g(:,:,:)
        complex(q), allocatable :: dvdu_longrange(:)
        complex(q), allocatable :: ddijdu_longrange(:,:,:,:,:)
        integer, allocatable :: gvecs(:,:)
        integer, allocatable :: nindpw(:)
        logical, allocatable :: gmask(:)
        integer, allocatable :: map_phelel2vasp(:)
        type(communic) :: comm

        !> if one of epsilon, Z, Q are present, then LR contribution is computed
        if ((.not.is_identity(3,self%dielectric_tensor)).or.&
            any(self%born_eff_charges/=0).or.&
            any(self%quadrupoles/=0)) self%has_lr = .true.
        if (.not.self%has_lr) return
        PROFILING_START('elph_pot_setup_lr')

        !> store cutoff for the long-range part and compute the ewald parameter from it
        self%encutlr = encutlr
        self%ewald_param = get_ewald_parameter(self%dielectric_tensor,self%encutlr)

        call get_snf_scp(self)
        self%rcut = 0.0_q

        ! compute the smallest grid spacing in each direction
        ! self%grid_point are sp grid point in [0,1[
        ! self%fft_mesh_sc =nint ( [1/min(self%grid_point(1)) 1/min(self%grid_point(2)) 1/min(self%grid_point(3))])
        call estimate_mesh(self%ngridpoints_local,self%grid_point,self%fft_mesh_sc,grid%comm)

        ! sometimes to map the potential in the SC to the primitive cell requires
        ! performing FFTs on very large regular meshes. In that case we fallback to slow fourier interpolation
        self%use_lr_sft = use_lr_sft
        if (1.0_q*product(self%fft_mesh_sc)/self%ngridpoints>16.0_q.and.self%use_lr_sft==0) self%use_lr_sft = 1

        ! Gamma centered sphered of g vectors
        vqpt = 0._q

        ! create a dense grid_sc grid in the supercell, along the supercell lattice
        ! vectors, the mapping to the phelel grid, and the g vectors reciprocal to grid_sc
        CALLMPI(comm_init_dummy(comm))
        if (self%use_lr_sft>0) then
            call get_grid_in_supercell(self,comm,vqpt,self%grid_sc,ngvecs,gvecs,nindpw)
        else
            call get_grid_in_supercell(self,comm,vqpt,self%grid_sc,ngvecs,gvecs,nindpw,map_phelel2vasp)
        endif

        ! Compute reciprocal space components v(g), for the previously computed g's,
        !  of the long range potnetiel : v_long_range_g
        allocate(gmask(ngvecs))
        allocate(v_long_range_g(ngvecs,3,self%natoms_pc))
        sp_at=matmul(self%prim2sp,self%centers) ! in supercell (self%latt_sp) coordinates
        call dvdu_long_range_g(self%latt_sp,vqpt,self%natoms_pc,&
                        self%born_eff_charges,self%quadrupoles,self%dielectric_tensor,&
                        sp_at,self%rcut,self%ewald_param,ngvecs,gvecs,gmask,&
                        v_long_range_g)

        ! compute dvdu_longrange
        PROFILING_START('elph_dvdu_setup_lr')
        ! number of non zero coef in v(g)
        self%ngvecs_sc = count(gmask)

        ! the potential in real space, on the phelel sp grid
        allocate(dvdu_longrange(self%ngridpoints_local))
        !> perform the fft (sum_g) to get the real space points on grid_sc, but only keep the
        !> points that belong to phelel grid self%grid_point
        !>     vr = sum(g) vg exp(+igr) (VASP convention)
        !>        = sum(g) (4pi/V_sp) e^2/(4pi epsilon_0)
        !>                            * (i g*Z_(alpha,atom) + Q )/ (g.epsilon.g)
        !>                            * exp(ig.(r-R_atom))
        !>                            * exp(-g.epsilon.g/(4 ewald_parameter^2))
        !>
        !> V_sp is the volule of the supercell
        !> g is the set of g vectors defined by grid_sc,self%fft_mesh_sc(1),self%fft_mesh_sc(2),self%fft_mesh_sc(3)
        !>  dvdu_longrange(i,idir,iatom), with i the phelel grid
        !> with this method we sum over more g, since we invent a dense grid, than the corresponding true real grid. I
        !> am not sure how correct it is
        if (iu>0) write(iu,*) 'Computing LR part of the electron-phonon potential'
        do iatom=1,self%natoms_pc
            do idir=1,3
                if (self%comm_shmem%node_me/=1) cycle
                if (iu>0) write(iu,*) 'atom: ',iatom,'/',self%natoms_pc,' idir:',idir
                if (self%use_lr_sft>0) then
                    call get_lr_sft(self%ngridpoints_local,ngvecs,gvecs,gmask,self%grid_point,&
                                    v_long_range_g(:,idir,iatom),dvdu_longrange,self%comm_kbatch)
                else
                    call get_lr_fft(self%grid_sc,ngvecs,nindpw,self%ngridpoints_local,map_phelel2vasp,&
                                    v_long_range_g(:,idir,iatom),dvdu_longrange)
                endif
                if (self%ncdij==1) then
                    self%dvdu(:,idir,iatom,1) = self%dvdu(:,idir,iatom,1) - dvdu_longrange
                elseif (self%ncdij==2) then
                    self%dvdu(:,idir,iatom,1) = self%dvdu(:,idir,iatom,1) - dvdu_longrange
                    self%dvdu(:,idir,iatom,2) = self%dvdu(:,idir,iatom,2) - dvdu_longrange
                elseif (self%ncdij==4) then
                    self%dvdu(:,idir,iatom,1) = self%dvdu(:,idir,iatom,1) - dvdu_longrange
                    self%dvdu(:,idir,iatom,4) = self%dvdu(:,idir,iatom,4) - dvdu_longrange
                endif
                if (allocated(self%dvdu_longrange)) then
                    self%dvdu_longrange(:,idir,iatom) = dvdu_longrange
                endif
            enddo
        enddo
        call elphon_pot_shmem_barrier(self)
        PROFILING_STOP('elph_dvdu_setup_lr')

        ! compute ddijdu_longrange
        PROFILING_START('elph_ddijdu_setup_lr')
        allocate(ddijdu_longrange(self%lmdim,self%lmdim,self%natoms_sc,3,self%natoms_pc_local))
        call get_ddijdu_longrange(self,gmask,ngvecs,gvecs,v_long_range_g,ddijdu_longrange,self%comm_inb)
        if (self%ncdij==1) then
            self%ddijdu(:,:,:,:,:,1) = self%ddijdu(:,:,:,:,:,1) - ddijdu_longrange
        elseif (self%ncdij==2) then
            self%ddijdu(:,:,:,:,:,1) = self%ddijdu(:,:,:,:,:,1) - ddijdu_longrange
            self%ddijdu(:,:,:,:,:,2) = self%ddijdu(:,:,:,:,:,2) - ddijdu_longrange
        elseif (self%ncdij==4) then
            self%ddijdu(:,:,:,:,:,1) = self%ddijdu(:,:,:,:,:,1) - ddijdu_longrange
            self%ddijdu(:,:,:,:,:,4) = self%ddijdu(:,:,:,:,:,4) - ddijdu_longrange
        endif
        deallocate(ddijdu_longrange)
        PROFILING_STOP('elph_ddijdu_setup_lr')

        self%has_lr_active = .true.

        deallocate(gmask)
        deallocate(v_long_range_g)

        ! generate g-sphere for evaluation of the LR part of the potential
        vqpt = 0._q
        call get_gsphere(self%grid_pc,self%latt_pc%b,2*self%encutlr,vqpt,self%ngvecs_pc_local,self%gvecs_local,self%nindpw_local)

        self%comm = grid%comm
        block
        type(grid_3d) :: grid_tmp
        integer,allocatable :: nindpw_dummy(:)
        type(communic) :: comm_dummy
        CALLMPI(comm_init_dummy(comm_dummy))
        ! initialize a grid to perform the FFT
        call grid3d_init(grid_tmp,self%fft_mesh(1),self%fft_mesh(2),self%fft_mesh(3),&
                         comm_dummy,comm_dummy,comm_dummy,real2cplx=.false.,lreal=.false.,lparallel=.false.)
        call get_gsphere(grid_tmp,self%latt_pc%b,2*self%encutlr,vqpt,self%ngvecs_pc,self%gvecs,nindpw_dummy)
        allocate(self%gmask(self%ngvecs_pc))
        deallocate(nindpw_dummy)
        end block
        PROFILING_STOP('elph_pot_setup_lr')

    end subroutine elphon_pot_setup_lr

    !> Compute the long-range part of the ddij terms
    subroutine get_ddijdu_longrange(self,gmask,ngvecs,gvecs,v_long_range_g,ddijdu_longrange,comm_inb)
        use constant
        use mpimy, only: communic
        use wave, only: ni_local
        type (elphon_pot) :: self
        logical, allocatable, intent(in)  :: gmask(:)
        integer, intent(in) :: ngvecs
        integer, allocatable, intent(in) :: gvecs(:,:)
        complex(q), allocatable, intent(inout) :: v_long_range_g(:,:,:)
        complex(q), allocatable, intent(inout) :: ddijdu_longrange(:,:,:,:,:)
        type(communic), intent(in) :: comm_inb
        ! local variables
        integer :: ngvecs_filtered
        real(q),allocatable :: gvecs_tmp(:,:)
        complex(q) :: angle
        integer ig,j,iatom_local,iatom_global,iatom,idir,ntypes,iatom_sc,ityp,lm,lmp
        complex(q),allocatable :: qijb_ityp(:,:,:,:)
        complex(q),allocatable :: phases(:,:)

        ! Get g-vectors list and compress v_long_range_g
        ngvecs_filtered = count(gmask)
        allocate(gvecs_tmp(3,ngvecs_filtered))
        ig = 0
        do j=1,ngvecs
           if (.not.gmask(j)) cycle
           ig = ig+1
           gvecs_tmp(:,ig)=gvecs(:,j)
           ! compress storage of v_long_range_g (remove zero entries)
           do iatom_global=1,self%natoms_pc
              iatom_local=ni_local(iatom_global,comm_inb)
              if (iatom_local==0) cycle
              do idir=1,3
                 v_long_range_g(ig,idir,iatom_local)=v_long_range_g(j,idir,iatom_global)
              enddo
           enddo
        enddo

        ntypes = size(self%p)
        allocate(qijb_ityp(ngvecs_filtered,self%lmdim,self%lmdim,ntypes))
        ! xxx may be we have this ???
        ! xxx qijb_ityp(g,l_Im_I,l_Jm_J,ityp)=int d^3r (phi_I(r) phi_J(r)-
        ! \tilde{phi}_I(r) \tilde{phi}_J(r)) exp(+ i g r)
        ! xxx exp(+ i g r) and not exp(- i g r) to be compatible witht the phase
        ! coming next
        ! xxx but check i^L and not (-i)^L -> I believe OK
        call getqijb_ityp(ngvecs_filtered,gvecs_tmp,ntypes,self%p,&
                          self%lmdim,self%latt_sp,qijb_ityp)

        ! Compute atom phase factors
        ! xxx exp( i g R_i) with i in supercell
        allocate(phases(ngvecs_filtered,self%natoms_sc))
        do iatom_sc=1,self%natoms_sc
          iatom = self%supercell_atoms_map_from_phelel(iatom_sc)
          do ig=1,ngvecs_filtered
            angle = gvecs_tmp(1,ig)*self%supercell_atoms(1,iatom)+&
                    gvecs_tmp(2,ig)*self%supercell_atoms(2,iatom)+&
                    gvecs_tmp(3,ig)*self%supercell_atoms(3,iatom)
            phases(ig,iatom_sc) = exp(citpi*angle)
          enddo
        enddo
        deallocate(gvecs_tmp)

        ! Compute dijdu lr contribution
        ! xxx lm or lmp and is the phase correct ?
        ! xxx what formula do we want to code ?
        ! xxx ddijdu_longrange(lm,lmp,iatom_sc,idir,iatom)=
        ! xxx      sum(g)  QIJB(NGVECS,LMDIM,LMDIM,NTYPES)
        ! xxx                            *(4pi/V_sp) e^2/(4pi epsilon_0)
        ! xxx                            * (i Z_atom*g + Q )/ (g.epsilon.g)
        ! xxx                            * exp(ig.(R_i-R_atom))
        ! xxx                            * exp(-g.epsilon.g/(4 ewald_parameter^2))
        do iatom=1,self%natoms_pc_local
          do idir=1,3
            do iatom_sc=1,self%natoms_sc
              ityp = self%ityp_sc(iatom_sc)
              do lmp=1,self%lmdim
              do lm=1,self%lmdim
                ddijdu_longrange(lm,lmp,iatom_sc,idir,iatom) = &
                sum(qijb_ityp(:,lm,lmp,ityp)*v_long_range_g(:ngvecs_filtered,idir,iatom)*phases(:,iatom_sc))
              enddo
              enddo
            enddo
          enddo
        enddo
        deallocate(qijb_ityp)

    end subroutine get_ddijdu_longrange

    !> Build mapping from phelel supercell to vasp
    !> [A1' A2' A3' ]= [A1 A2 A3] Q
    !> qin=Q^{-1}
    !> grid_point = reduced coordinates along [A1 A2 A3] axes
    !> grid_sc = reduced coordinates along [A1' A2' A3']
    !> map_phelel2vasp(index in grid_point)=index in grid_sc
    !> map_regular2phelel(index in grid_sc )=index in grid_point
    !>                                       or -1 if grid_point is not in  grid_sc
    subroutine get_mapping(ngridpoints_local, grid_point, qinv, grid_sc, map_phelel2vasp, map_regular2phelel)
        use string, only: str
        use tutor, only: vtutor
        integer, intent(in) :: ngridpoints_local
        real(q), intent(in) :: grid_point(:,:)
        integer, intent(in) :: qinv(3,3)
        type(grid_3d), intent(in) :: grid_sc
        integer, allocatable, intent(out) :: map_phelel2vasp(:)
        integer, allocatable, optional, intent(out) :: map_regular2phelel(:)
        ! local variables
        integer :: i,ivasp,iphelel,ig,nx,ny,nz,nc
        integer :: fft_mesh_sc(3)
        integer,allocatable :: map(:,:,:,:)
        real(q) :: nr(3)

        fft_mesh_sc(1) = grid_sc%ngx
        fft_mesh_sc(2) = grid_sc%ngy
        fft_mesh_sc(3) = grid_sc%ngz

        ! create mapping from phelel2vasp
        allocate(map(2,fft_mesh_sc(1),fft_mesh_sc(2),fft_mesh_sc(3)))
        map = -1
        do i=1,ngridpoints_local
            ! integral coordinates along [A1' A2' A3']
            nr = matmul(qinv,grid_point(:,i))*fft_mesh_sc
            ! check if integer
            if (any(abs((nint(nr)-nr))>1e-6)) call vtutor%bug('grid_point not integer in this FFT grid: '//str(nr),__FILE__,__LINE__)
            ! ASK_HM
            nx = modulo(nint(nr(1)),fft_mesh_sc(1))+1
            ny = modulo(nint(nr(2)),fft_mesh_sc(2))+1
            nz = modulo(nint(nr(3)),fft_mesh_sc(3))+1
            map(1,nx,ny,nz)=i ! the point i in grid_point as integral coordinates nx,ny,nz along [A1' A2' A3']
        enddo
        ! index of vasp grid on regular grid
        if (grid_sc%rl%nfast==1) then
            do nc=1,grid_sc%rl%ncol
               ny=grid_sc%rl%i2(nc)
               nz=grid_sc%rl%i3(nc)
               do nx=1,grid_sc%rl%nrow
                  ivasp=nx+(ny-1)*fft_mesh_sc(1)+(nz-1)*fft_mesh_sc(1)*fft_mesh_sc(2)
                  map(2,nx,ny,nz)=ivasp
               enddo
            enddo
        elseif (grid_sc%rl%nfast==3) then
            do nc=1,grid_sc%rl%ncol
               nx=grid_sc%rl%i2(nc)
               ny=grid_sc%rl%i3(nc)
               do nz=1,grid_sc%rl%nrow
                  ivasp=nx+(ny-1)*fft_mesh_sc(1)+(nz-1)*fft_mesh_sc(1)*fft_mesh_sc(2)
                  map(2,nx,ny,nz)=ivasp
               enddo
            enddo
        endif
        ! allocate mapping
        allocate(map_phelel2vasp(ngridpoints_local))
        if (present(map_regular2phelel)) allocate(map_regular2phelel(grid_sc%ngx*grid_sc%ngy*grid_sc%ngz))
        map_phelel2vasp=-1
        if (present(map_regular2phelel)) map_regular2phelel=-1
        ig = 0
        do nz=1,fft_mesh_sc(3)
          do ny=1,fft_mesh_sc(2)
            do nx=1,fft_mesh_sc(1)
              ig = ig+1
              iphelel = map(1,nx,ny,nz)
              ivasp   = map(2,nx,ny,nz)
              if (iphelel>0) map_phelel2vasp(iphelel)=ivasp
              if (present(map_regular2phelel)) map_regular2phelel(ig)=iphelel
            enddo
          enddo
        enddo
        !if (any(map_phelel2vasp<1))    call vtutor%error('Failed to create mapping from phelel to vasp')
        !if (any(map_regular2phelel<1)) call vtutor%error('Failed to create mapping from regular to vasp')
        deallocate(map)
    end subroutine get_mapping

    !> get the lr part of the potential in real space using slow fourier transformation
    subroutine get_lr_sft(ngridpoints,ngvecs,gvecs,gmask,grid_point,v_long_range_g,dvdu_longrange,comm)
        use constant, only: tpi
        use tutor, only: vtutor
        use string, only: str
#ifdef MKL_VML
        include "mkl_vml.f90"
#endif
        integer,intent(in) :: ngridpoints
        integer,intent(in) :: ngvecs
        integer,intent(in) :: gvecs(:,:)
        logical,intent(in) :: gmask(:)
        real(q),intent(in) :: grid_point(:,:) !< local grid points
        complex(q),intent(in) :: v_long_range_g(:)
        complex(q) :: dvdu_longrange(:)
        type(communic),intent(in) :: comm
        !local variables
        integer :: i, j
        real(q) :: nr(3), gvec2pi(3)
        real(q),allocatable :: angles(:)
        complex(q),allocatable :: phase(:)
        if (ngridpoints/=size(dvdu_longrange)) then
            call vtutor%bug('Incompatible dimensions in the potential array wanted '//str(ngridpoints)//' found '//str(size(dvdu_longrange)),__FILE__,__LINE__)
        endif
        if (ngridpoints/=size(grid_point,2)) then
            call vtutor%bug('Incompatible dimensions in the grid_point array wanted '//str(ngridpoints)//' found '//str(size(grid_point,2)),__FILE__,__LINE__)
        endif
        allocate(angles(ngridpoints))
        allocate(phase(ngridpoints))
        dvdu_longrange = cmplx(0.0_q,0.0_q,q)
        ! Loop over G vectors
        do j=1,ngvecs
#ifdef MPI
           if (mod(j-1,comm%ncpu)/=comm%node_me-1) cycle
#endif
           if (.not.gmask(j)) cycle
           gvec2pi = tpi*gvecs(:,j)
           do i=1,ngridpoints
             nr = grid_point(:,i)
             angles(i) = gvec2pi(1)*nr(1)+&
                         gvec2pi(2)*nr(2)+&
                         gvec2pi(3)*nr(3)
           enddo
#ifdef MKL_VML
           call vmzcis( ngridpoints, angles, phase, VML_EP )
#else
           do i=1,ngridpoints
             phase(i) = exp(ci*angles(i))
           enddo
#endif
           dvdu_longrange(:) = &
           dvdu_longrange(:) + &
           v_long_range_g(j)*phase(:)
        enddo ! loop over g-vectors
        CALLMPI( M_sum_z(comm, dvdu_longrange, size(dvdu_longrange)) )
        deallocate(angles)
        deallocate(phase)
    end subroutine get_lr_sft

    !> get the lr part of the potential in real space using fast fourier transformation
    subroutine get_lr_fft(grid,ngvecs,nindpw,ngridpoints_local,map_phelel2vasp,v_long_range_g,dvdu_longrange)
        use string, only: str
        type(grid_3d) :: grid
        integer,intent(in) :: ngvecs
        integer,intent(in) :: nindpw(:)
        integer,intent(in) :: ngridpoints_local
        integer,intent(in) :: map_phelel2vasp(:)
        complex(q),intent(in) :: v_long_range_g(:) !> potential in reciprocal space
        complex(q),intent(inout) :: dvdu_longrange(:)
        ! local variables
        integer :: i, ivasp
        complex(q),allocatable :: dvdu_longrange_tmp(:)
        allocate(dvdu_longrange_tmp(grid%ngx*grid%ngy*grid%ngz))
        ! call FFTWAV to get V(r) on the auxiliary grid defined by VASP
        !>    vr= sum(g) vg exp(+igr)
        call FFTWAV(ngvecs,nindpw,dvdu_longrange_tmp(:),&
                    v_long_range_g(:),grid)

        ! set V(r) from the auxiliary grid to the phelel grid
        do i=1,ngridpoints_local
            ivasp = map_phelel2vasp(i)
            dvdu_longrange(i) = dvdu_longrange_tmp(ivasp)
        enddo
        deallocate(dvdu_longrange_tmp)
    end subroutine get_lr_fft

    !> write the lattice parameters and the atomic positions in the POSCAR format
    subroutine write_latt_atoms(iu,p,lattice,natoms,atom_positions,atom_types)
        integer :: iu
        type (potcar) :: p(:)
        real(q) :: lattice(3,3)
        integer,intent(in) :: natoms
        real(q),intent(in) :: atom_positions(3,natoms)
        integer,intent(in) :: atom_types(natoms)
        ! local variables
        integer :: i
        if (iu>0) then
            write(iu,*) 'electron-phonon potential'
            write(iu,*) 1.0
            write(iu,*) lattice(:,1)
            write(iu,*) lattice(:,2)
            write(iu,*) lattice(:,3)
            do i=1,size(p)
                write(iu,"(A,' ')",advance='no') p(i)%element
            enddo
            write(iu,*)
            do i=1,size(p)
                write(iu,"(I4,' ')",advance='no') count(atom_types==i)
            enddo
            write(iu,*)
            write(iu,*) 'Direct'
            do i=1,natoms
              write(iu,*) atom_positions(:,i)
            enddo
            write(iu,*)
        endif
    end subroutine write_latt_atoms

    !> Write the potential to a text file (same formar as LOCPOT file)
    subroutine write_pot(iu,ngx,ngy,ngz,ngridpoints,pot,mapping)
        integer,intent(in) :: iu
        integer,intent(in) :: ngx,ngy,ngz
        integer,intent(in) :: ngridpoints
        complex(q),intent(in) :: pot(ngridpoints)
        integer,optional,intent(in) :: mapping(ngridpoints)
        !local variables
        integer :: i, imap
        write(iu,'(3i5)') ngx,ngy,ngz
        do i=1,ngridpoints
          imap = i; if (present(mapping)) imap = mapping(i)
          if ( mod(i,5) == 0 ) then
            write(iu,'(1(1x,e17.11))') real(pot(imap))
          else
            write(iu,'(1(1x,e17.11))',advance='no') real(pot(imap))
          endif
        enddo
    end subroutine write_pot

    !> @brief Get the electron-phonon potential on the primitive cell
    subroutine elphon_pot_get_dvdu_prim(self,vqpt,has_lr_active,dvdu,dvdu_kkp)
        type(elphon_pot) :: self
        real(q),intent(in) :: vqpt(3)
        logical :: has_lr_active
        !> Potential in the supercell in real space
        complex(q),target,contiguous,intent(in) :: dvdu(:,:,:,:)
        !> Potential in the primitive cell at a particular q in real space
        complex(q),target,contiguous,intent(out) :: dvdu_kkp(:,:,:,:)
        ! Local variables
        integer :: n,wf_ncdij
        complex(q) :: dvdu_kkp_lr(self%npr_local,3,self%natoms_pc)
        complex(q),contiguous,pointer :: dvdu_kkp_ptr(:,:,:,:)
        complex(q),pointer :: a(:,:)
        complex(q),pointer :: b(:,:)

        ! early return if set
        if (self%skip) return

        ! check dimensions of the wanted potential
        wf_ncdij = size(dvdu_kkp,4)
        if (wf_ncdij==self%ncdij) then
            dvdu_kkp_ptr => dvdu_kkp
        else
            allocate(dvdu_kkp_ptr(self%npr_local,3,self%natoms_pc,self%ncdij))
        endif

        ! compute short and long-range parts
        call elphon_pot_get_dvdu_prim_sr(self,vqpt,dvdu,dvdu_kkp_ptr)
        if (has_lr_active) then
            call elphon_pot_get_dvdu_prim_lr(self,vqpt,dvdu_kkp_lr)
            select case(self%ncdij)
              case(1)
                dvdu_kkp_ptr(:,:,:,1) = dvdu_kkp_ptr(:,:,:,1) + dvdu_kkp_lr
              case(2)
                dvdu_kkp_ptr(:,:,:,1) = dvdu_kkp_ptr(:,:,:,1) + dvdu_kkp_lr
                dvdu_kkp_ptr(:,:,:,2) = dvdu_kkp_ptr(:,:,:,2) + dvdu_kkp_lr
              case(3)
                dvdu_kkp_ptr(:,:,:,1) = dvdu_kkp_ptr(:,:,:,1) + dvdu_kkp_lr
                dvdu_kkp_ptr(:,:,:,4) = dvdu_kkp_ptr(:,:,:,4) + dvdu_kkp_lr
            end select
        endif

        ! transfer the potential if necessary
        if (wf_ncdij/=self%ncdij) then
            n = self%npr_local*3*self%natoms_pc
            a(1:n,1:self%ncdij) => dvdu_kkp_ptr
            b(1:n,1:wf_ncdij) => dvdu_kkp
            call transfer_ncdij(a,self%ncdij,b,wf_ncdij)
            deallocate(dvdu_kkp_ptr)
        endif
    end subroutine elphon_pot_get_dvdu_prim

    subroutine transfer_ncdij(a,a_ncdij,b,b_ncdij)
        use tutor, only: vtutor
        ! local variables
        complex(q),intent(in)  :: a(:,:)
        complex(q),intent(out) :: b(:,:)
        integer,intent(in) :: a_ncdij
        integer,intent(in) :: b_ncdij
        select case(b_ncdij)
            case(1)
                select case(a_ncdij)
                    case(1)
                        b=a
                    case(2)
                        call vtutor%error('Transfering data from potential with ISPIN=2 to a spin unpolarized calculation is not implemented')
                        b(:,1)=a(:,1)+a(:,2)
                    case(4)
                        b(:,1)=a(:,1)+a(:,4)
                end select
            case(2)
                select case(a_ncdij)
                    case(1)
                        b(:,1)=a(:,1)
                        b(:,2)=a(:,1)
                    case(2)
                        b=a
                    case(4)
                        b(:,1)=a(:,1)
                        b(:,2)=a(:,4)
                end select
            case(4)
                select case(a_ncdij)
                    case(1)
                        b(:,1)=a(:,1)
                        b(:,2)=0.0_q
                        b(:,3)=0.0_q
                        b(:,4)=a(:,1)
                    case(2)
                        b(:,1)=a(:,1)
                        b(:,2)=0.0_q
                        b(:,3)=0.0_q
                        b(:,4)=a(:,2)
                    case(4)
                        b=a
                end select
        end select
    end subroutine transfer_ncdij

    !> @brief Test interpolation with the special treatment of the long-range part of the potential
    subroutine elphon_test_lr(self,nkpts,vkpt)
        use lattice
        use vhdf5_base
        use elphon_base
        type (elphon_pot) :: self
        integer,intent(in) :: nkpts
        real(q),intent(in) :: vkpt(3,nkpts)
        ! Local variables
        complex(q), allocatable :: dvdu_kkp(:,:,:,:)
        complex(q), allocatable :: dvdu_kkp_lr(:,:,:,:)
        complex(q), allocatable :: ddijdu_kkp(:,:,:,:,:,:)
        complex(q), allocatable :: ddijdu_kkp_lr(:,:,:,:,:,:)
        complex(q), allocatable :: dvdu_avg(:,:,:,:)
        complex(q), allocatable :: dvdu_sr_avg(:,:,:,:)
        complex(q), allocatable :: dvdu_lr_avg(:,:,:,:)
        complex(q), allocatable :: ddijdu_avg(:,:,:,:)
        complex(q), allocatable :: ddijdu_sr_avg(:,:,:,:)
        complex(q), allocatable :: ddijdu_lr_avg(:,:,:,:)
        integer :: iq, iatom, idir

        ! Compute the long-range part of the potential
        allocate(dvdu_kkp(self%npw,3,self%natoms_pc,self%ncdij))
        allocate(dvdu_kkp_lr(self%npw,3,self%natoms_pc,1))
        allocate(ddijdu_kkp(self%lmdim,self%lmdim,self%natoms_pc,3,self%natoms_pc,self%ncdij))
        allocate(ddijdu_kkp_lr(self%lmdim,self%lmdim,self%natoms_pc,3,self%natoms_pc,1))

        allocate(dvdu_avg(3,self%natoms_pc,self%ncdij,nkpts))
        allocate(dvdu_sr_avg(3,self%natoms_pc,self%ncdij,nkpts))
        allocate(dvdu_lr_avg(3,self%natoms_pc,self%ncdij,nkpts))
        allocate(ddijdu_avg(3,self%natoms_pc,self%ncdij,nkpts))
        allocate(ddijdu_sr_avg(3,self%natoms_pc,self%ncdij,nkpts))
        allocate(ddijdu_lr_avg(3,self%natoms_pc,self%ncdij,nkpts))

        ! initialize writting of the potential along a path
        do iq=1,nkpts
            write(*,*) '[',iq,'/',nkpts,']'
            call elphon_pot_get_dvdu_prim(self,vkpt(:,iq),self%has_lr_active,self%dvdu,dvdu_kkp)
            call elphon_pot_get_ddijdu_prim(self,vkpt(:,iq),ddijdu_kkp)
            do iatom=1,self%natoms_pc
                do idir=1,3
                    dvdu_avg(idir,iatom,1,iq) = sum(dvdu_kkp(:,idir,iatom,1))/self%npw
                    ddijdu_avg(idir,iatom,1,iq) = sum(ddijdu_kkp(:,:,:,idir,iatom,1))/self%natoms_pc
                enddo
            enddo
        enddo

        do iq=1,nkpts
            write(*,*) '[',iq,'/',nkpts,']'
            ! compute potential without long-range part
            call elphon_pot_get_dvdu_prim_sr(self,vkpt(:,iq),self%dvdu,dvdu_kkp)
            call elphon_pot_get_ddijdu_prim_sr(self,vkpt(:,iq),ddijdu_kkp)
            ! compute long-range potential
            if (self%has_lr_active) then
                call elphon_pot_get_dvdu_prim_lr(self,vkpt(:,iq),dvdu_kkp_lr(:,:,:,1))
                call elphon_pot_get_ddijdu_prim_lr(self,vkpt(:,iq),ddijdu_kkp_lr(:,:,:,:,:,1))
            endif

            ! save information for file
            do iatom=1,self%natoms_pc
                do idir=1,3
                    dvdu_sr_avg(idir,iatom,1,iq) = sum(dvdu_kkp(:,idir,iatom,1))/self%npw
                    dvdu_lr_avg(idir,iatom,1,iq) = sum(dvdu_kkp_lr(:,idir,iatom,1))/self%npw
                    ddijdu_sr_avg(idir,iatom,1,iq) = sum(ddijdu_kkp(:,:,:,idir,iatom,1))/self%natoms_pc
                    ddijdu_lr_avg(idir,iatom,1,iq) = sum(ddijdu_kkp_lr(:,:,:,idir,iatom,1))/self%natoms_pc
                enddo
            enddo
        enddo

        ! write information to hdf5 file
#ifdef VASP_HDF5
        call vh5_write_dvdu_avg(ih5outfileid,self%natoms_pc,self%ncdij,&
                                dvdu_avg,dvdu_sr_avg,dvdu_lr_avg,&
                                ddijdu_avg,ddijdu_sr_avg,ddijdu_lr_avg,nkpts,vkpt)
#endif
        deallocate(dvdu_avg)
        deallocate(dvdu_sr_avg)
        deallocate(dvdu_lr_avg)
        deallocate(ddijdu_avg)

        deallocate(dvdu_kkp)
        deallocate(dvdu_kkp_lr)
        deallocate(ddijdu_kkp)
        deallocate(ddijdu_kkp_lr)
    end subroutine elphon_test_lr

    !> Check the computation of the long-range part of the potential
    !> by comparing the results using slow and fast fourier transform.
    !>
    !> By default the long-range part is always computed using the fast fourier transformation
    !> only in this test we compare with the slow version.
    subroutine elphon_check_lr_slowfastft(self,iu0)
        use tutor, only: vtutor
        use string, only: str
        type (elphon_pot) :: self
        integer,intent(in) :: iu0
        ! local variables
        integer :: iatom, idir
        integer :: ngvecs
        real(q) :: vqpt(3)
        real(q) :: maxdiff, avgdiff
        real(q) :: sp_at(3,self%natoms_pc)
        complex(q), allocatable :: dvdu_longrange_sft(:)
        complex(q), allocatable :: dvdu_longrange_fft(:)
        complex(q), allocatable :: v_long_range_g(:,:,:)
        integer, allocatable :: map_phelel2vasp(:)
        integer, allocatable :: nindpw(:)
        logical, allocatable :: gmask(:)
        integer, allocatable :: gvecs(:,:)

        if (iu0>=0) write(iu0,*) 'Check LR with slow and fast Fourier transform'

        vqpt = 0._q
        call get_grid_in_supercell(self,self%comm,vqpt,self%grid_sc,ngvecs,gvecs,nindpw,map_phelel2vasp)

        allocate(gmask(ngvecs))
        allocate(v_long_range_g(ngvecs,3,self%natoms_pc))

        ! Compute the potential in reciprocal space
        sp_at=matmul(self%prim2sp,self%centers) ! in supercell (self%latt_sp) coordinates
        call dvdu_long_range_g(self%latt_sp,vqpt,self%natoms_pc,&
                               self%born_eff_charges,self%quadrupoles,self%dielectric_tensor,&
                               sp_at,self%rcut,self%ewald_param,ngvecs,gvecs,gmask,&
                               v_long_range_g)

        ! allocate and compute the potential using slow fourier transform
        allocate(dvdu_longrange_sft(self%ngridpoints))
        allocate(dvdu_longrange_fft(self%ngridpoints))

        ! compare the results
        do iatom=1,self%natoms_pc
            do idir=1,3
                call get_lr_sft(self%ngridpoints,ngvecs,gvecs,gmask,self%grid_point,&
                                v_long_range_g(:,idir,iatom),dvdu_longrange_sft,self%grid_pc%comm_kinter)

                call get_lr_fft(self%grid_sc,ngvecs,nindpw,self%ngridpoints,map_phelel2vasp,&
                                v_long_range_g(:,idir,iatom),dvdu_longrange_fft)

                write(*,*) 'iatom:',iatom, 'idir:',idir
                maxdiff = maxval(abs(dvdu_longrange_sft(:)-&
                                     dvdu_longrange_fft(:)))
                avgdiff = sum(abs(dvdu_longrange_sft(:)-&
                                  dvdu_longrange_fft(:)))/self%ngridpoints
                if (avgdiff>1e-7 .or. maxdiff>1e-6) then
                    call vtutor%error('Difference found between computation with fast and slow Fourier transform.'//&
                                      ' maxdiff: '//str(maxdiff)//' > 1e-8  avgdiff:'//str(avgdiff)//' > 1e-8')
                endif
                write(*,*) 'max diff:   ', maxdiff
                write(*,*) 'avg diff:   ', avgdiff
                write(*,*) 'max fft:    ', maxval(real(dvdu_longrange_fft(:))),&
                                           maxval(aimag(dvdu_longrange_fft(:)))
                write(*,*) 'max slow ft:', maxval(real(dvdu_longrange_sft(:))),&
                                           maxval(aimag(dvdu_longrange_sft(:)))
                write(*,*) 'min fft:    ', minval(real(dvdu_longrange_fft(:))),&
                                           minval(aimag(dvdu_longrange_fft(:)))
                write(*,*) 'min slow ft:', minval(real(dvdu_longrange_sft(:))),&
                                           minval(aimag(dvdu_longrange_sft(:)))
            enddo
        enddo
    end subroutine elphon_check_lr_slowfastft

    !> Check if removing the LR part in the supercell and adding it back in the primitive cell
    !> results in a potential close to zero
    subroutine elphon_check_lr_addremove(self,nkpts,vkpt)
        type (elphon_pot) :: self
        integer,intent(in) :: nkpts
        real(q),intent(in) :: vkpt(3,nkpts)
        ! local variables
        integer :: iq, iatom, idir
        integer :: ngvecs
        real(q) :: vqpt(3)
        real(q) :: sp_at(3,self%natoms_pc)
        complex(q), allocatable :: v_long_range_g(:,:,:)
        complex(q), allocatable :: dvdu_longrange_save(:,:,:,:)
        complex(q), allocatable :: dvdu_longrange_tmp(:)
        complex(q), allocatable :: ddijdu_longrange(:,:,:,:,:)
        complex(q), allocatable :: dvdu_kkp(:,:,:,:)
        complex(q), allocatable :: ddijdu_kkp(:,:,:,:,:,:)
        integer, allocatable :: gvecs(:,:)
        integer, allocatable :: nindpw(:)
        logical, allocatable :: gmask(:)
        integer, allocatable :: map_phelel2vasp(:)

        ! Gamma centered sphered of g vectors
        vqpt = 0._q

        !> create a dense grid_sc grid in the supercell, along the supercell lattice
        !> vectors, the mapping to the phelel grid, and the g vectors reciprocal to grid_sc
        call get_grid_in_supercell(self,self%comm,vqpt,self%grid_sc,ngvecs,gvecs,nindpw,map_phelel2vasp)

        ! Compute reciprocal space components v(g), for the previously computed g's,
        !  of the long range potnetiel : v_long_range_g
        allocate(gmask(ngvecs))
        allocate(v_long_range_g(ngvecs,3,self%natoms_pc))
        sp_at=matmul(self%prim2sp,self%centers) ! in supercell (self%latt_sp) coordinates
        call dvdu_long_range_g(self%latt_sp,vqpt,self%natoms_pc,&
                        self%born_eff_charges,self%quadrupoles,self%dielectric_tensor,&
                        sp_at,self%rcut,self%ewald_param,ngvecs,gvecs,gmask,&
                        v_long_range_g)

        ! set potential in the SC to be the LR part
        allocate(dvdu_longrange_tmp(self%ngridpoints))
        allocate(dvdu_longrange_save(self%ngridpoints,3,self%natoms_pc,self%ncdij))
        dvdu_longrange_save = self%dvdu

        ! set dvdu
        do iatom=1,self%natoms_pc
            do idir=1,3
                call get_lr_fft(self%grid_sc,ngvecs,nindpw,self%ngridpoints,map_phelel2vasp,&
                                v_long_range_g(:,idir,iatom),dvdu_longrange_tmp)
                self%dvdu(:,idir,iatom,1) = -dvdu_longrange_tmp
            enddo
        enddo
        ! set ddijju
        allocate(ddijdu_longrange(self%lmdim,self%lmdim,self%natoms_sc,3,self%natoms_pc))
        call get_ddijdu_longrange(self,gmask,ngvecs,gvecs,v_long_range_g,ddijdu_longrange,self%comm_inb)
        self%ddijdu(:,:,:,:,:,1) = -ddijdu_longrange
        deallocate(v_long_range_g)

        ! compute the potential in the primitive cell and check if its zero
        allocate(dvdu_kkp(self%npw,3,self%natoms_pc,self%ncdij))
        allocate(ddijdu_kkp(self%lmdim,self%lmdim,self%natoms_pc,3,self%natoms_pc,self%ncdij))
        do iq=1,nkpts
            call elphon_pot_get_dvdu_prim(self,vkpt(:,iq),self%has_lr_active,self%dvdu,dvdu_kkp)
            call elphon_pot_get_ddijdu_prim(self,vkpt(:,iq),ddijdu_kkp)
            write(*,*) vkpt(:,iq)
            write(*,*) 'dvdu:', sum(dvdu_kkp)/self%npw
            write(*,*) 'ddij:', sum(ddijdu_kkp)
        enddo

        ! restore
        self%dvdu = dvdu_longrange_save
    end subroutine elphon_check_lr_addremove

    !> write potentials in the supercell in the format X Y Z Re(V)  Im(V)
    !> where X Y Z are the cartesian coordinates of the grid points, and V the potential
    subroutine elphon_write_pot(self,iu0)
        use tutor, only: vtutor
        use string, only: str
        type (elphon_pot) :: self
        integer,intent(in) :: iu0
        integer :: iatom, idir,ind,ispinor
        integer :: i, ix, iy, iz
        integer, allocatable :: map(:)
        real(q) xsnf(3,self%ngridpoints),x_atoms_snf(3,self%natoms_sc)
        real(q) lattice_sp_snf(3,3)
        integer NG(3,3)
        type(snf) :: snf_fft
        ! SNF transformed supercell
        NG = 0
        NG(1,1) = self%fft_mesh(1)
        NG(2,2) = self%fft_mesh(2)
        NG(3,3) = self%fft_mesh(3)
        call snf_init(snf_fft,matmul(NG,nint(self%sp2prim)))

        lattice_sp_snf=matmul(self%lattice_sp,snf_fft%qm)

        x_atoms_snf=matmul(snf_fft%qinv,self%supercell_atoms)
        do i=1,self%natoms_sc
            x_atoms_snf(1,i)=modulo(x_atoms_snf(1,i),1.0_q)
            x_atoms_snf(2,i)=modulo(x_atoms_snf(2,i),1.0_q)
            x_atoms_snf(3,i)=modulo(x_atoms_snf(3,i),1.0_q)
        enddo

        xsnf=matmul(snf_fft%qinv,self%grid_point)
        allocate(map(snf_fft%d(1)*snf_fft%d(2)*snf_fft%d(3)))
        map=-1
        do i=1,self%ngridpoints
            xsnf(1,i)=xsnf(1,i)*snf_fft%d(1)
            xsnf(2,i)=xsnf(2,i)*snf_fft%d(2)
            xsnf(3,i)=xsnf(3,i)*snf_fft%d(3)

            ix=modulo(nint(xsnf(1,i)),snf_fft%d(1))+1
            iy=modulo(nint(xsnf(2,i)),snf_fft%d(2))+1
            iz=modulo(nint(xsnf(3,i)),snf_fft%d(3))+1
            ind=(iz-1)*snf_fft%d(1)*snf_fft%d(2)+&
                (iy-1)*snf_fft%d(1)+&
                ix
            map(ind)=i
        enddo
        if (any(map(:) < 0)) then
            call vtutor%error('Failed to create mapping in elphon_write_pot')
        endif
        if (iu0>=0)  then
        write(iu0,*) 'Write potentials derivatives in supercell to files'


        do iatom=1,self%natoms_pc
            do idir=1,3
                write(iu0,*) 'iatom:',iatom, 'idir:',idir
                ! write cell infos
                open(unit=99,file='pot_dvdu.'//str(iatom)//'.'//str(idir)//'.dat')
                call write_latt_atoms(99,self%p,lattice_sp_snf,self%natoms_sc,x_atoms_snf,self%ityp_sc)
                ! write potential
                do ispinor=1,self%ncdij
                    call write_pot(99,snf_fft%d(1),snf_fft%d(2),snf_fft%d(3),product(snf_fft%d),&
                                   self%dvdu(:,idir,iatom,ispinor),map)
                    write(99,*) ! to separate spin(or) components
                enddo
                close(99)

                if (.not. allocated(self%dvdu_longrange)) cycle

                ! write cell infos
                open(unit=99,file='pot_dvdu_lr.'//str(iatom)//'.'//str(idir)//'.dat')
                call write_latt_atoms(99,self%p,lattice_sp_snf,self%natoms_sc,x_atoms_snf,self%ityp_sc)
                ! write potential
                call write_pot(99,snf_fft%d(1),snf_fft%d(2),snf_fft%d(3),product(snf_fft%d),&
                               self%dvdu_longrange(:,idir,iatom),map)
                write(99,*) ! to separate spin(or) components
                close(99)
            enddo
        enddo

        endif
        deallocate(map)
    end subroutine elphon_write_pot

    !> Write the potential to a file in the disk using the LOCPOT format
    subroutine elphon_write_lr(self,iu0)
        use tutor, only: vtutor
        use string, only: str
        type (elphon_pot) :: self
        integer,intent(in) :: iu0
        ! local variables
        integer :: iatom, idir
        real(q) :: lattice_spp(3,3)
        integer,allocatable :: map_phelel2vasp(:)
        integer,allocatable :: map_regular2phelel(:)
        real(q),allocatable :: supercell_positions_poscar(:,:)
        type(grid_3d) :: grid_sc

        if (iu0>=0) write(iu0,*) 'write LR part of the potential'

        ! create supercell grid
        call grid3d_init(grid_sc,self%fft_mesh_sc(1),self%fft_mesh_sc(2),self%fft_mesh_sc(3),&
                         self%comm,self%comm,self%comm,real2cplx=.false.,lreal=.false.,lparallel=.false.)

        ! create mapping from regular grid to phelel: map_regular2phelel
        ! map_phelel2vasp is dummy (not used)
        ! note that the coordinates in grid_point have to be transformed by the Q matrix
        call get_mapping(self%ngridpoints_local, self%grid_point, self%snf_fft%qinv, grid_sc, map_phelel2vasp, map_regular2phelel)
        if (any(map_regular2phelel<1)) call vtutor%bug('Failed to create mapping from regular to vasp',__FILE__,__LINE__)

        ! reorder atomic supercell positions so as to have the same order as in a poscar file
        allocate(supercell_positions_poscar(3,self%natoms_sc))
        do iatom=1,self%natoms_sc
            supercell_positions_poscar(:,iatom) = matmul(self%snf_fft%qinv,self%supercell_atoms(:,self%supercell_atoms_map_from_phelel(iatom)))
        enddo

        ! Compute the lattice of the supercell prime
        lattice_spp = matmul(self%lattice_sp,self%snf_fft%qm)

        do iatom=1,self%natoms_pc
            do idir=1,3
                write(*,*) 'iatom:',iatom, 'idir:',idir
                ! write lr potential to file
                !open(unit=99,file='potlr.'//str(iatom)//'.'//str(idir)//'.dat')
                !call write_latt_atoms(99,self%p,lattice_spp,self%natoms_sc,supercell_positions_poscar,self%ityp_sc)
                !call write_pot(99,self%fft_mesh_sc(1),self%fft_mesh_sc(2),self%fft_mesh_sc(3),&
                !           self%ngridpoints,self%dvdu_longrange(:,idir,iatom),map_regular2phelel)
                !close(99)
                ! write potential to file
                open(unit=99,file='potsr.'//str(iatom)//'.'//str(idir)//'.dat')
                call write_latt_atoms(99,self%p,lattice_spp,self%natoms_sc,supercell_positions_poscar,self%ityp_sc)
                call write_pot(99,self%fft_mesh_sc(1),self%fft_mesh_sc(2),self%fft_mesh_sc(3),&
                           self%ngridpoints,self%dvdu(:,idir,iatom,1),map_regular2phelel)
                close(99)
            enddo
        enddo
    end subroutine elphon_write_lr

    !> @brief Get the short-range part of the potential in the primitive cell at a given q-point
    subroutine elphon_pot_get_dvdu_prim_sr(self,qpt,dvdu,dvdu_kkp)
        use constant
#ifdef MKL_VML
        include "mkl_vml.f90"
#endif
        type(elphon_pot) :: self
        real(q),intent(in) :: qpt(3)
        complex(q),intent(in) :: dvdu(:,:,:,:)
        complex(q),intent(out) :: dvdu_kkp(:,:,:,:)
        ! local variables
        integer :: idir, ispinor
        integer :: i, j, k, h
        integer :: sizep,shift,shifp
        integer :: atoms_disp
        real(q) :: vec(3),qpt2pi(3)
        complex(q) :: phase
        complex(q) :: phases(self%npr_local)
        real(q) :: angles_all(self%max_points_lattice)
        complex(q) :: phases_all(self%max_points_lattice)
        PROFILING_START('elph_potential_sr')

        dvdu_kkp=cmplx(0.0_q,0.0_q,q)
        qpt2pi=-tpi*qpt

        !$omp parallel default(none) private(atoms_disp,shift) shared(self,dvdu,qpt2pi,dvdu_kkp,angles_all,phases_all,phases)
        do atoms_disp=1,self%natoms_pc
           ! compute angles
           !$omp do private(h,vec) schedule(static)
           do h=1,self%max_points_lattice
              vec=self%sp_lattice_points(:,h,atoms_disp) ! D convention
              angles_all(h)=qpt2pi(1)*vec(1)+qpt2pi(2)*vec(2)+qpt2pi(3)*vec(3)
           enddo
           !$omp end do

           ! compute phases
#ifdef MKL_VML
           call vmzcis( self%max_points_lattice, angles_all, phases_all, VML_EP )
#else
           !$omp do private(h) schedule(static)
           do h=1,self%max_points_lattice
              phases_all(h) = exp(ci*angles_all(h))
           enddo
           !$omp end do
#endif

           do i=1,self%det_p2s
              shift = (i-1)*self%npr_local

              ! average phases
              !$omp do private(k,j,phase,sizep,shifp) schedule(static)
              do k=1,self%npr_local
                 sizep = self%sp_lattice_points_size(1,shift+k,atoms_disp)
                 shifp = self%sp_lattice_points_size(2,shift+k,atoms_disp)
                 phase = cmplx(0.0_q,0.0_q,q)
                 do j=1,sizep
                    phase=phase+phases_all(shifp+j)
                 enddo
                 phases(k) = phase/sizep
              enddo ! loop over gridpoints
              !$omp end do

              ! multiply phases
              !$omp do private(k,ispinor,idir) collapse(3) schedule(static)
              do ispinor=1,self%ncdij
                 do idir=1,3
                    do k=1,self%npr_local
                        dvdu_kkp(k,idir,atoms_disp,ispinor)=&
                        dvdu_kkp(k,idir,atoms_disp,ispinor)+&
                        phases(k)*dvdu(shift+k,idir,atoms_disp,ispinor)
                    enddo ! loop over points
                 enddo ! loop over directions
              enddo ! loop over spinors
              !$omp end do
           enddo ! loop over primitive cells in supercell

        enddo ! loop over atoms
        !$omp end parallel

        PROFILING_STOP('elph_potential_sr')

    end subroutine elphon_pot_get_dvdu_prim_sr

    !> @brief Compute the long-range part of the potential
    subroutine elphon_pot_get_dvdu_prim_lr(self,vqpt,v_lr_r)
        use tutor, only: vtutor
        use string, only: str
        type(elphon_pot) :: self
        real(q),intent(in) :: vqpt(3)
        complex(q),intent(out) :: v_lr_r(:,:,:) !< Long range potential in real space
        ! Local variables
        integer :: idir, iatom
        complex(q) :: v_lr_r_tmp(self%grid_pc%mplwv)
        complex(q) :: v_lr_g(self%ngvecs_pc_local,3,self%natoms_pc) !< Long range potential in reciprocal space

        PROFILING_START('elph_potential_lr')
        if (self%grid_pc%rl%np /= size(v_lr_r,1)) then
            call vtutor%bug('Wrong first dimension of v_lr_r array. found: '//str(size(v_lr_r,1))//&
                            ' expected: '//str(self%grid_pc%rl%np),__FILE__,__LINE__)
        endif
        call dvdu_long_range_g(self%latt_pc,vqpt,self%natoms_pc,&
                               self%born_eff_charges,self%quadrupoles,self%dielectric_tensor,&
                               self%centers,self%rcut,self%ewald_param,self%ngvecs_pc_local,self%gvecs_local,self%gmask,v_lr_g)
        !write(*,*) 'pot', count(self%gmask), self%ngvecs_pc

        ! Fourier transform to real space
        do iatom=1,self%natoms_pc
            do idir=1,3
                call FFTWAV(self%ngvecs_pc_local,self%nindpw_local,v_lr_r_tmp,&
                            v_lr_g(:,idir,iatom),self%grid_pc)
                v_lr_r(:,idir,iatom) = v_lr_r_tmp(:self%grid_pc%rl%np)
            enddo
        enddo
        PROFILING_STOP('elph_potential_lr')

    end subroutine elphon_pot_get_dvdu_prim_lr

    subroutine elphon_pot_get_ddijdu_prim(self,vqpt,ddijdu_kkp)
        use pseudo, only: potcar
        type(elphon_pot) :: self
        real(q),intent(in) :: vqpt(3)
        complex(q),target,contiguous,intent(inout) ::  ddijdu_kkp(:,:,:,:,:,:)
        ! Local variables
        integer :: n,wf_ncdij
        complex(q) :: ddijdu_kkp_lr(self%lmdim,self%lmdim,self%natoms_pc,3,self%natoms_pc_local)
        complex(q),contiguous,pointer :: ddijdu_kkp_ptr(:,:,:,:,:,:)
        complex(q),pointer :: a(:,:)
        complex(q),pointer :: b(:,:)

        ! early return if set
        if (self%skip) return

        ! check the dimensions of the wanted potential
        wf_ncdij = size(ddijdu_kkp,6)
        if (wf_ncdij==self%ncdij) then
            ddijdu_kkp_ptr => ddijdu_kkp
        else
            allocate(ddijdu_kkp_ptr(self%lmdim,self%lmdim,self%natoms_pc,3,self%natoms_pc,self%ncdij))
        endif

        ! compute short and long-range parts
        call elphon_pot_get_ddijdu_prim_sr(self,vqpt,ddijdu_kkp_ptr)
        if (self%has_lr_active) then
            call elphon_pot_get_ddijdu_prim_lr(self,vqpt,ddijdu_kkp_lr)
            select case (self%ncdij)
              case(1)
                ddijdu_kkp_ptr(:,:,:,:,:,1) = ddijdu_kkp_ptr(:,:,:,:,:,1) + ddijdu_kkp_lr
              case(2)
                ddijdu_kkp_ptr(:,:,:,:,:,1) = ddijdu_kkp_ptr(:,:,:,:,:,1) + ddijdu_kkp_lr
                ddijdu_kkp_ptr(:,:,:,:,:,2) = ddijdu_kkp_ptr(:,:,:,:,:,2) + ddijdu_kkp_lr
              case(4)
                ddijdu_kkp_ptr(:,:,:,:,:,1) = ddijdu_kkp_ptr(:,:,:,:,:,1) + ddijdu_kkp_lr
                ddijdu_kkp_ptr(:,:,:,:,:,4) = ddijdu_kkp_ptr(:,:,:,:,:,4) + ddijdu_kkp_lr
            end select
        endif

        !transfer the potential if necessary
        if (wf_ncdij/=self%ncdij) then
            n = self%lmdim*self%lmdim*self%natoms_pc*3*self%natoms_pc_local
            a(1:n,1:self%ncdij) => ddijdu_kkp_ptr
            b(1:n,1:wf_ncdij) => ddijdu_kkp
            call transfer_ncdij(a,self%ncdij,b,wf_ncdij)
            deallocate(ddijdu_kkp_ptr)
        endif
    end subroutine elphon_pot_get_ddijdu_prim

    !> @brief Get the derivative of dij in the primitive cell
    subroutine elphon_pot_get_ddijdu_prim_sr(self,qpt,ddijdu_kkp)
        use constant
        use wave, only : ni_local
        type(elphon_pot) :: self
        real(q),intent(in) :: qpt(3)
        complex(q), intent(inout) ::  ddijdu_kkp(:,:,:,:,:,:)
        ! local variables
        real(q) :: vec(3),angle
        complex(q) :: phase,weight
        integer :: shift, sizep
        integer :: atoms_disp,atoms_disp_local,i1,is,is_phelel
        integer :: i,j,alpha,spinors

        PROFILING_START('elph_strength_sr')

        ddijdu_kkp=0

        is=0
        do i=1,self%det_p2s
           do i1=1,self%natoms_pc
              is=is+1

              is_phelel = self%supercell_atoms_map_to_phelel(is)
              do atoms_disp=1,self%natoms_pc
                 atoms_disp_local = ni_local(atoms_disp,self%comm_inb)
                 if (atoms_disp_local==0) cycle
                 weight=0
                 sizep = self%sp_lattice_atoms_size(1,is,atoms_disp)
                 shift = self%sp_lattice_atoms_size(2,is,atoms_disp)
                 do j=1,sizep
                    !vec=self%primitive_positions(:,i1)+self%lattice_point(:,i)+&
                    !    self%sp_lattice_atoms(:,shift+j,atoms_disp)-&
                    !    self%centers(:,atoms_disp) ! C convention
                    vec=self%primitive_positions(:,i1)+self%lattice_point(:,i)+&
                        self%sp_lattice_atoms(:,shift+j,atoms_disp) ! D convention
                    angle=qpt(1)*vec(1)+qpt(2)*vec(2)+qpt(3)*vec(3)
                    phase=exp(-citpi*angle)
                    weight=weight+phase
                 enddo
                 weight=weight/sizep

                 do spinors=1,self%ncdij
                    do alpha =1,3
                       ddijdu_kkp(:,:,i1,alpha,atoms_disp_local,spinors)=&
                       ddijdu_kkp(:,:,i1,alpha,atoms_disp_local,spinors)+&
                       weight*self%ddijdu(:,:,is_phelel, alpha, atoms_disp_local, spinors)
                    enddo
                 enddo ! spinors
              enddo ! atoms_disp

           enddo
        enddo
        PROFILING_STOP('elph_strength_sr')

    end subroutine elphon_pot_get_ddijdu_prim_sr

    !> Compute the long-range part of the change of the PAW strength in the primtive cell
    !> This corresponds to the dDij/du which can be computed as:
    !> dDij/du = /int Q_{ij}^L(r) V_{LR}(r) dr
    !> V^{LR}_q(r) = i e/4pi \sum_{G/=q} [ (q+G).Zeff ]/[(q+G).eps.(q+G)] e^i(q+G).(r+R)
    subroutine elphon_pot_get_ddijdu_prim_lr(self,vqpt,ddijdu_kkp_lr)
        use constant
        use pseudo, only: potcar
        use wave, only: ni_local
        type(elphon_pot) :: self
        real(q),intent(in) :: vqpt(3)
        complex(q), intent(inout) ::  ddijdu_kkp_lr(:,:,:,:,:)
        ! local variables
        integer :: iatom_local, iatom_global, iatom_pc, idir, j
        integer :: ngvecs_filtered, ntypes
        integer :: lm, lmp, ig, ityp
        real(q) :: angle
        real(q),allocatable :: gvecs_tmp(:,:)
        real(q),allocatable :: qpgvecs_tmp(:,:)
        complex(q) :: v_lr_g(self%ngvecs_pc,3,self%natoms_pc) !< Long range potential in reciprocal space
        complex(q),allocatable :: qijb_ityp(:,:,:,:)
        complex(q),allocatable :: phases(:,:)

        PROFILING_START('elph_ddijdu_lr')
        call dvdu_long_range_g(self%latt_pc,vqpt,self%natoms_pc,&
                               self%born_eff_charges,self%quadrupoles,self%dielectric_tensor,&
                               self%centers,self%rcut,self%ewald_param,self%ngvecs_pc,self%gvecs,self%gmask,v_lr_g)

        ! Get g-vectors list and compress v_lr_g
        ngvecs_filtered = count(self%gmask)
        allocate(gvecs_tmp(3,ngvecs_filtered))
        allocate(qpgvecs_tmp(3,ngvecs_filtered))
        ig = 0
        do j=1,self%ngvecs_pc
            if (.not.self%gmask(j)) cycle
            ig = ig+1
            qpgvecs_tmp(:,ig)=vqpt+self%gvecs(:,j)
            gvecs_tmp(:,ig)=self%gvecs(:,j)
            ! compress storage of v_lr_g (remove zero entries)
            do iatom_global=1,self%natoms_pc
                do idir=1,3
                v_lr_g(ig,idir,iatom_global)=v_lr_g(j,idir,iatom_global)
                enddo
            enddo
        enddo

        ! Compute qijb for all the g-vectors and atom types
        ntypes = size(self%p)
        allocate(qijb_ityp(ngvecs_filtered,self%lmdim,self%lmdim,ntypes))
        call getqijb_ityp(ngvecs_filtered,qpgvecs_tmp,ntypes,self%p,&
                            self%lmdim,self%latt_pc,qijb_ityp)
        deallocate(qpgvecs_tmp)

        ! Compute atom phase factors
        allocate(phases(ngvecs_filtered,self%natoms_pc))
        do iatom_pc=1,self%natoms_pc
            do ig=1,ngvecs_filtered
            angle = gvecs_tmp(1,ig)*self%primitive_positions(1,iatom_pc)+&
                    gvecs_tmp(2,ig)*self%primitive_positions(2,iatom_pc)+&
                    gvecs_tmp(3,ig)*self%primitive_positions(3,iatom_pc)
            phases(ig,iatom_pc) = exp(citpi*angle)
            enddo
        enddo
        deallocate(gvecs_tmp)

        ! Compute dijdu lr contribution
        do iatom_global=1,self%natoms_pc
            iatom_local = ni_local(iatom_global, self%comm_inb)
            if (iatom_local==0) cycle
            do idir=1,3
            do iatom_pc=1,self%natoms_pc
                ityp = self%ityp_pc(iatom_pc)
                do lmp=1,self%lmdim
                do lm=1,self%lmdim
                ddijdu_kkp_lr(lm,lmp,iatom_pc,idir,iatom_local) = &
                sum(qijb_ityp(:,lm,lmp,ityp)*v_lr_g(:ngvecs_filtered,idir,iatom_global)*phases(:,iatom_pc))
                enddo
                enddo
            enddo
            enddo
        enddo
        deallocate(qijb_ityp)
        deallocate(phases)

        PROFILING_STOP('elph_ddijdu_lr')

    end subroutine elphon_pot_get_ddijdu_prim_lr

    !> @brief Free memory of elphon potential datastructure
    subroutine elphon_pot_free(self)
        use ini, only: deregister_allocate
        use mgrid, only: dealloc_grd
        type(elphon_pot) :: self
        deallocate(self%lattice_point)
        deallocate(self%centers)
        deallocate(self%primitive_positions)
        deallocate(self%atom_indices_in_derivatives)
        deallocate(self%born_eff_charges)
        deallocate(self%quadrupoles)
        deallocate(self%phonon_masses)

        call deregister_allocate(VMEM_ARR_SIZE(self%grid_point), 'pot_grid')
        call deregister_allocate(VMEM_ARR_SIZE(self%dvdu)/self%comm_shmem%ncpu, 'pot_dvdu')

        if (self%has_shmem) then
#ifdef shmem_elphon
            call M_shmem_dealloc(self%dvdu_shmem)
#endif
        else
            deallocate(self%dvdu)
        endif
        if (allocated(self%dmudu)) deallocate(self%dmudu)
        deallocate(self%grid_point)
        deallocate(self%ccdij)
        deallocate(self%ccqij)
        deallocate(self%ddijdu)
        deallocate(self%dqijdu)
        deallocate(self%ccdrij_radial)
        deallocate(self%ccdrij)

        call deregister_allocate(VMEM_ARR_SIZE(self%sp_lattice_points)/self%comm_shmem%ncpu+&
                                 VMEM_ARR_SIZE(self%sp_lattice_points_size)/self%comm_shmem%ncpu, 'pot_lattice_points')
        if (self%has_shmem) then
#ifdef shmem_elphon
            call M_shmem_dealloc(self%sp_lattice_points_shmem)
#endif
        else
            deallocate(self%sp_lattice_points)
        endif
        deallocate(self%sp_lattice_points_size)

        if (self%has_lr_active) then
            !call dealloc_grd(self%grid_pc)
            call dealloc_grd(self%grid_sc)
        endif
    end subroutine elphon_pot_free

    !> @brief Transform grid_point to interval [0,1[, and result the shift shift_sp
    !> grid_point = r + l + shift_sp, where r is in the primitive cell and l the
    !> lattice vectors from phelel
    subroutine grid_in_01(grid,fft_mesh,prim2sp,det_p2s,ngridpoints,lattice_point,grid_point,shift_sp)
        type(grid_3d), intent(in) :: grid
        integer, intent(in) :: fft_mesh(3)
        real(q), intent(in) :: prim2sp(3,3)
        integer, intent(in) :: det_p2s
        integer, intent(in) :: ngridpoints
        real(q), intent(in) :: lattice_point(3,det_p2s)
        real(q), intent(inout) :: grid_point(3,ngridpoints)
        !> the shift from (r+l) to grid_point. Must be integer.
        real(q), intent(out):: shift_sp(3,ngridpoints)
        integer i,nc,ix,iy,iz,ind

        ind=0
        if (grid%rl%nfast==3) then
           do i=1,det_p2s
           do nc=1,grid%rl%ncol
              ix = grid%rl%i2(nc)-1
              iy = grid%rl%i3(nc)-1
              do iz=0,grid%rl%nrow-1
                 ind=ind+1
                 shift_sp(:,ind)=(1.0_q*[ix,iy,iz])/fft_mesh(:)+lattice_point(:,i)
              enddo
           enddo
           enddo
        else
           do i=1,det_p2s
           do nc=1,grid%rl%ncol
              iy = grid%rl%i2(nc)-1
              iz = grid%rl%i3(nc)-1
              do ix=0,grid%rl%nrow-1
                 ind=ind+1
                 shift_sp(:,ind)=(1.0_q*[ix,iy,iz])/fft_mesh(:)+lattice_point(:,i)
              enddo
           enddo
           enddo
        endif

        do i=1,ngridpoints
           shift_sp(:,i)=matmul(prim2sp,shift_sp(:,i)) ! now wrt supercell lattice vectors
        enddo

        do i=1,ngridpoints
           grid_point(1,i)=modulo(grid_point(1,i),1.0_q)
           grid_point(2,i)=modulo(grid_point(2,i),1.0_q)
           grid_point(3,i)=modulo(grid_point(3,i),1.0_q)
        enddo

        shift_sp=grid_point-shift_sp

    end subroutine grid_in_01

    subroutine make_supercell_atoms(prim2sp,lattice_point,primitive_positions,supercell_positions,&
                                    supercell_atoms,supercell_atoms_map_to_phelel,&
                                    supercell_atoms_map_from_phelel,shift_sp)
        use tutor, only : vtutor
        use sym_prec, only: tiny
        real(q), intent(in) :: prim2sp(3,3)
        real(q),allocatable, intent(in) :: lattice_point(:,:)
        real(q),allocatable, intent(in) :: primitive_positions(:,:)
        real(q),allocatable, intent(inout) :: supercell_positions(:,:)
        real(q),allocatable, intent(out) :: supercell_atoms(:,:)
        integer,allocatable, intent(out) :: supercell_atoms_map_to_phelel(:)
        integer,allocatable, intent(out) :: supercell_atoms_map_from_phelel(:)
        real(q),allocatable, intent(out) :: shift_sp(:,:)
        ! local variables
        integer :: i1,i2,is
        integer :: natoms_pc,natoms_sc
        real(q) :: diff(3)

        natoms_pc=size(primitive_positions,2)
        natoms_sc=size(supercell_positions,2)

        ! ordering in phelel is not convenient yet
        allocate(supercell_atoms(1:3,1:natoms_sc))
        allocate(supercell_atoms_map_to_phelel(1:natoms_sc))
        allocate(supercell_atoms_map_from_phelel(1:natoms_sc))
        allocate(shift_sp(3,natoms_sc))

        ! make supercell
        is=0
        do i2=1,size(lattice_point,2)
           do i1=1,natoms_pc
              is=is+1
              shift_sp(:,is)= primitive_positions(:,i1)+lattice_point(:,i2)
           enddo
        enddo
        shift_sp=matmul(prim2sp,shift_sp)
        do is=1,natoms_sc
           do i1=1,3
              supercell_atoms(i1,is) = modulo(shift_sp(i1,is),1.0_q)
           enddo
        enddo
        shift_sp=supercell_atoms-shift_sp

        ! map to phelel ordereing: to be supressed, i hope
        supercell_atoms_map_to_phelel=-1
        do i1=1,natoms_sc
            do i2=1,natoms_sc
               diff = supercell_atoms(:,i1)-supercell_positions(:,i2)
               if (all(abs(diff-nint(diff)) < tiny)) then
                  if (supercell_atoms_map_to_phelel(i1) > 0) then
                      call vtutor%error('error in sp-atom mapping to phelel (1)')
                  else
                      supercell_atoms_map_to_phelel(i1)=i2
                      supercell_atoms_map_from_phelel(i2)=i1
                  endif
               endif
            enddo
            if (supercell_atoms_map_to_phelel(i1) < 0) then
               call vtutor%error('error in sp-atom mapping to phelel (2)')
            endif
        enddo

    deallocate(supercell_positions)
    end subroutine make_supercell_atoms

    !>@brief Find the r space coordinates of points closest to the atoms
    subroutine atoms_centered_grids(lattice_sp,natoms_pc,ngridpoints,centers,grid_point,shift_sp, &
                                    max_nlpoints,sp_lattice_points,sp_lattice_points_size)
        use ini, only: register_allocate
        real(q), intent(in) :: lattice_sp(3,3) !< supercell lattice vectors
        integer, intent(in) :: natoms_pc
        integer, intent(in) :: ngridpoints
        !> Reduced SP coordinates of the atoms where the potential derivatives are known
        real(q), intent(in) :: centers(:,:)
        real(q), intent(in) :: grid_point(:,:)
        real(q), intent(in) :: shift_sp(:,:)
        !> equivalent SP lattice points around a given atom
        real(q), allocatable, intent(out) :: sp_lattice_points(:,:,:)
        integer, allocatable, intent(out) :: sp_lattice_points_size(:,:,:)
        integer, intent(out) :: max_nlpoints

        PROFILING_START('atoms_centered_grids')

        ! Determine max_nlpoints
        call atoms_centered_grids_count(lattice_sp,natoms_pc,ngridpoints,centers,grid_point,max_nlpoints)

        ! Now that we know max_nlpoints we can allocate sp_lattice_points
        allocate(sp_lattice_points(3,max_nlpoints,natoms_pc))
        allocate(sp_lattice_points_size(2,ngridpoints,natoms_pc))

        ! fill in the values of the points in sp_lattice_points
        call atoms_centered_grids_generate(lattice_sp,natoms_pc,ngridpoints,centers,grid_point,shift_sp, &
                                           sp_lattice_points,sp_lattice_points_size)

        PROFILING_STOP('atoms_centered_grids')

    end subroutine atoms_centered_grids

    !>@brief Find the r space coordinates of points closest to the atoms
    subroutine atoms_centered_grids_count(lattice_sp,natoms_pc,ngridpoints,centers,grid_point, &
                                          max_nlpoints)
        use ini, only: register_allocate
        real(q), intent(in) :: lattice_sp(3,3) !< supercell lattice vectors
        integer, intent(in) :: natoms_pc
        integer, intent(in) :: ngridpoints
        !> Reduced SP coordinates of the atoms where the potential derivatives are known
        real(q), intent(in) :: centers(:,:)
        real(q), intent(in) :: grid_point(:,:)
        integer, intent(out) :: max_nlpoints
        !Local variables
        real(q) :: lattice_sp_perp(3,3) ! supercell lattice vectors
        real(q) :: lattice_sp_norm2(3)
        integer :: nlpoints ! number of closest SP lattice points
        integer :: total_nlpoints
        integer :: lpoints(3,12)  ! list of closest SP lattice points
        integer :: i,j
        real(q) :: x(3)

        PROFILING_START('atoms_centered_grids_count')

        call ddecomp_a_per(lattice_sp,lattice_sp_perp)
        lattice_sp_norm2(1) = dot_product(lattice_sp_perp(:,1),lattice_sp_perp(:,1))
        lattice_sp_norm2(2) = dot_product(lattice_sp_perp(:,2),lattice_sp_perp(:,2))
        lattice_sp_norm2(3) = dot_product(lattice_sp_perp(:,3),lattice_sp_perp(:,3))

        ! Determine max_nlpoints
        max_nlpoints=0
        do i=1,natoms_pc
           total_nlpoints=0
           do j=1,ngridpoints
              x=centers(:,i)-grid_point(:,j)
              call closest_lattice_sp(lattice_sp,lattice_sp_norm2,lattice_sp_perp,x,lpoints,nlpoints)
              total_nlpoints = total_nlpoints+nlpoints

           enddo
           max_nlpoints=max(total_nlpoints,max_nlpoints)
        enddo

        PROFILING_STOP('atoms_centered_grids_count')

    end subroutine atoms_centered_grids_count

    !>@brief Find the r space coordinates of points closest to the atoms
    subroutine atoms_centered_grids_generate(lattice_sp,natoms_pc,ngridpoints,centers,grid_point,shift_sp, &
                                             sp_lattice_points,sp_lattice_points_size)
        use ini, only: register_allocate
        real(q), intent(in) :: lattice_sp(3,3) !< supercell lattice vectors
        integer, intent(in) :: natoms_pc
        integer, intent(in) :: ngridpoints
        !> Reduced SP coordinates of the atoms where the potential derivatives are known
        real(q), intent(in) :: centers(:,:)
        real(q), intent(in) :: grid_point(:,:)
        real(q), intent(in) :: shift_sp(:,:)
        !> equivalent SP lattice points around a given atom
        real(q), intent(inout) :: sp_lattice_points(:,:,:)
        integer, intent(inout) :: sp_lattice_points_size(:,:,:)
        !Local variables
        real(q) :: lattice_sp_perp(3,3) ! supercell lattice vectors
        real(q) :: lattice_sp_norm2(3)
        integer :: nlpoints ! number of closest SP lattice points
        integer :: lpoints(3,12)  ! list of closest SP lattice points
        integer :: i,j,k,ind
        real(q) :: x(3)

        PROFILING_START('atoms_centered_grids_generate')

        call ddecomp_a_per(lattice_sp,lattice_sp_perp)
        lattice_sp_norm2(1) = dot_product(lattice_sp_perp(:,1),lattice_sp_perp(:,1))
        lattice_sp_norm2(2) = dot_product(lattice_sp_perp(:,2),lattice_sp_perp(:,2))
        lattice_sp_norm2(3) = dot_product(lattice_sp_perp(:,3),lattice_sp_perp(:,3))

        ! fill in the values of the points in sp_lattice_points
        do i=1,natoms_pc
           ind=0
           do j=1,ngridpoints
              x=centers(:,i)-grid_point(:,j)

              !lpoints= ArgMin_L |L+grid_point-center|, with L supercell lattice vectors
              call closest_lattice_sp(lattice_sp,lattice_sp_norm2,lattice_sp_perp,x,lpoints,nlpoints)

              ! TODO: check if using nlpoints=1 changes the final results
              sp_lattice_points_size(1,j,i)=nlpoints
              sp_lattice_points_size(2,j,i)=ind
              do k=1,nlpoints
                 ind=ind+1
                 ! coordinates wrt supercell lattice vectors
                 sp_lattice_points(:,ind,i)=1.0_q*lpoints(:,k)+shift_sp(:,j)
              enddo

           enddo
        enddo

        PROFILING_STOP('atoms_centered_grids_generate')

    end subroutine atoms_centered_grids_generate

    !> @brief Find the closest lattice point
    subroutine closest_lattice_sp(lattice_sp,lattice_sp_norm2,lattice_sp_perp,xin,lpoints,nlpoints)
        use sym_prec, only: tiny
        real(q), intent(in) :: lattice_sp(3,3)
        real(q), intent(in) :: lattice_sp_norm2(3)
        real(q), intent(in) :: lattice_sp_perp(3,3)
        real(q), intent(in):: xin(3)
        integer, intent(out) :: lpoints(3,12) ! closest points
        integer, intent(out) :: nlpoints
        ! local variables
        real(q) :: rho
        integer :: i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
        integer :: i1,i2,i3
        real(q) :: x(3),p(3)
        real(q) :: min_value
        real(q) :: lnorm, dp
        real(q) :: dist

        x=matmul(lattice_sp,xin)

        rho=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

        dp = x(1)*lattice_sp_perp(1,1)+&
             x(2)*lattice_sp_perp(2,1)+&
             x(3)*lattice_sp_perp(3,1)
        lnorm  = sqrt(lattice_sp_norm2(1))
        i1_min=  floor(dp/lattice_sp_norm2(1)-rho/lnorm)
        i1_max=ceiling(dp/lattice_sp_norm2(1)+rho/lnorm)

        dp = x(1)*lattice_sp_perp(1,2)+&
             x(2)*lattice_sp_perp(2,2)+&
             x(3)*lattice_sp_perp(3,2)
        lnorm  = sqrt(lattice_sp_norm2(2))
        i2_min=  floor(dp/lattice_sp_norm2(2)-rho/lnorm)
        i2_max=ceiling(dp/lattice_sp_norm2(2)+rho/lnorm)

        dp = x(1)*lattice_sp_perp(1,3)+&
             x(2)*lattice_sp_perp(2,3)+&
             x(3)*lattice_sp_perp(3,3)
        lnorm  = sqrt(lattice_sp_norm2(3))
        i3_min=  floor(dp/lattice_sp_norm2(3)-rho/lnorm)
        i3_max=ceiling(dp/lattice_sp_norm2(3)+rho/lnorm)

        ! HM: This version is faster and shorter
        ! xxx WRONG(?) : I dont think it is correct if the atom is "on" a grid point
        nlpoints=0
        min_value=1e8
        do i1=i1_min,i1_max
          do i2=i2_min,i2_max
              do i3=i3_min,i3_max
                 p=lattice_sp(:,1)*i1+lattice_sp(:,2)*i2+lattice_sp(:,3)*i3-x
                 dist=sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
                 ! get minimum distance
                 if (min_value>dist) then
                    ! restart generatorr if the distance is very different from previous distance
                    if (abs(dist-min_value)>tiny) nlpoints=0
                    min_value=dist
                 endif
                 ! if distance is minimum add point
                 if (abs(dist-min_value)<tiny) then
                    nlpoints=nlpoints+1
                    lpoints(1,nlpoints) = i1
                    lpoints(2,nlpoints) = i2
                    lpoints(3,nlpoints) = i3
                 endif
             enddo
           enddo
        enddo
    end subroutine closest_lattice_sp

    !> @brief A_per contain as column the vectors A1_per,A2_per,A3_per, which are
    !>  components of A_i perpendicular to the other 2 vectors
    subroutine ddecomp_a_per(a,a_per)
        real(q), intent(in) :: a(1:3,1:3)
        real(q), intent(out) :: a_per(1:3,1:3)
        real(q)  a1_par(3,1)
        real(q)  a1_per(3,1)
        real(q)  a2_par(3,1)
        real(q)  a2_per(3,1)
        real(q)  a3_par(3,1)
        real(q)  a3_per(3,1)
        real(q)  a12(3,2)
        real(q)  a13(3,2)
        real(q)  a23(3,2)
        real(q) vec(3,1)

        a12(:,1)=a(:,1)
        a12(:,2)=a(:,2)

        a13(:,1)=a(:,1)
        a13(:,2)=a(:,3)

        a23(:,1)=a(:,2)
        a23(:,2)=a(:,3)

        vec(:,1)=a(:,1)
        call ddecomp_par_perp(a23,vec,a1_par,a1_per)

        vec(:,1)=a(:,2)
        call ddecomp_par_perp(a13,vec,a2_par,a2_per)

        vec(:,1)=a(:,3)
        call ddecomp_par_perp(a12,vec,a3_par,a3_per)

        a_per(:,1)=a1_per(:,1)
        a_per(:,2)=a2_per(:,1)
        a_per(:,3)=a3_per(:,1)

    end subroutine ddecomp_a_per


    !> @brief decompose the vector X into a part that is parallel to the range of A,
    !> and a part that is perpendicular to it: X = X_par + X_per
    !> the vectors and matrices must be allocated before the call.
    subroutine ddecomp_par_perp(a,x,x_par,x_per)
       real(q), intent(in) :: a(3,2)
       real(q), intent(in) :: x(3,1)
       real(q),  intent(out) :: x_par(3,1)
       real(q), intent(inout) :: x_per(3,1)
       integer i,j
       real(q)  g(2,2)
       real(q)  y(2,1)
       real(q)  z(2,1)
       real(q), allocatable :: work(:)
       integer, allocatable ::  ipiv(:)
       integer :: info

       ! sub-metric tensor
       do i=1,2
          do j=1,2
             g(i,j)=dot_product(a(:,i),a(:,j))
           enddo
       enddo

       ! compute inverse
       allocate(ipiv(1:2))
       allocate(work(1:2))
       call dgetrf( 2, 2,g, 2, ipiv, info )     ! which version to use ???
       call dgetri( 2, g, 2, ipiv, work, 2, info ) ! which version to use ???
       deallocate(work)
       deallocate(ipiv)

       ! projection  of X on the range of A
       do i=1,2
           y(i,1)=dot_product(x(:,1),a(:,i))
       enddo

       ! components  of X in the range of A
       z=matmul(g,y)

       ! X_par

       x_par=0.0d0
       do i=1,2
          do j=1,3
            x_par(j,1)=x_par(j,1)+z(i,1)*a(j,i)
          enddo
       enddo

       ! x_per
       x_per=x-x_par

    end subroutine ddecomp_par_perp

    !>@brief Add r + lattice_point to sp_lattice_points
    subroutine prepare_grid(grid,natoms_pc,det_p2s,np,fft_mesh,lattice_point,&
                            sp_lattice_points,sp_lattice_points_size)
        type(grid_3d), intent(in) :: grid
        integer, intent(in) :: natoms_pc
        integer, intent(in) :: det_p2s
        integer, intent(in) :: np !< number of points in the fft grid
        integer, intent(in) :: fft_mesh(3)
        real(q), intent(in) :: lattice_point(3,det_p2s)
        real(q), intent(inout) :: sp_lattice_points(:,:,:)
        integer, intent(in) :: sp_lattice_points_size(:,:,:)
        ! local variables
        integer :: i, j, k
        integer :: ix,iy,iz,nc
        integer :: sizep,shift,shifp
        integer :: atoms_disp
        real(q) :: vec(3)

        do atoms_disp=1,natoms_pc
           do i=1,det_p2s
              shift = (i-1)*np
              k=0
              if (grid%rl%nfast==3) then
                  do nc=1,grid%rl%ncol
                     ix = grid%rl%i2(nc)-1
                     vec(1)=1.0_q*ix/fft_mesh(1)
                     iy = grid%rl%i3(nc)-1
                     vec(2)=1.0_q*iy/fft_mesh(2)
                     do iz=0,grid%rl%nrow-1
                        vec(3)=1.0_q*iz/fft_mesh(3)
                        k=k+1
                        sizep = sp_lattice_points_size(1,shift+k,atoms_disp)
                        shifp = sp_lattice_points_size(2,shift+k,atoms_disp)
                        do j=1,sizep
                           sp_lattice_points(:,shifp+j,atoms_disp)=&
                           sp_lattice_points(:,shifp+j,atoms_disp)+&
                           vec+lattice_point(:,i) ! D convention
                        enddo !sizep
                     enddo !nrow
                  enddo ! nc
              else
                  do nc=1,grid%rl%ncol
                     iy = grid%rl%i2(nc)-1
                     vec(2)=1.0_q*iy/fft_mesh(2)
                     iz = grid%rl%i3(nc)-1
                     vec(3)=1.0_q*iz/fft_mesh(3)
                     do ix=0,grid%rl%nrow-1
                        vec(1)=1.0_q*ix/fft_mesh(1)
                        k=k+1
                        sizep = sp_lattice_points_size(1,shift+k,atoms_disp)
                        shifp = sp_lattice_points_size(2,shift+k,atoms_disp)
                        do j=1,sizep
                           sp_lattice_points(:,shifp+j,atoms_disp)=&
                           sp_lattice_points(:,shifp+j,atoms_disp)+&
                           vec+lattice_point(:,i) ! D convention
                        enddo !sizep
                     enddo !nrow
                  enddo ! nc
              endif
           enddo ! loop over primitive cells in supercell
        enddo ! loop over atoms

    end subroutine prepare_grid

    !> @brief Obtian the EWALD parameter in the same way as in phonon.F
    function get_ewald_parameter(dielectric_tensor,encutlr) result(ewald_param)
        !> ewald parameter is lambda in Eq. 73 of PRB 55, 10355 (1997)
        !> If we ask for the exponential to be less than alpha=exp(supression_threshold)
        !> then ewald_param = sqrt(2m Ecut/hbar^2 epsilon/(-4 supression_threshold))
        use constant, only: hsqdtm
        real(q),intent(in) :: dielectric_tensor(3,3)
        real(q) :: encutlr
        real(q) :: ewald_param
        ! local variables
        integer :: cart_idx
        real(q)  :: eigenvals(3)
        real(q)  :: dielectric_tensor_tmp(3, 3)
        real(q)  :: work(9)
        integer  :: info

        ! Take maximum eigenvalue of dielectric tensor as measure
        dielectric_tensor_tmp = dielectric_tensor
        call dsyev("n", "u", 3, dielectric_tensor_tmp, 3, eigenvals, work, size(work), info)
        ewald_param=maxval(eigenvals)
        ! Fallback previous method if diagonalization fails
        if (info /= 0) then
            ewald_param = 0
            do cart_idx = 1, 3
               ewald_param = ewald_param + dielectric_tensor(cart_idx, cart_idx)
            enddo
            ewald_param=ewald_param/3
        endif

        ewald_param= sqrt(encutlr/hsqdtm *  ewald_param  / ( -4* supression_threshold ) )
    end function get_ewald_parameter

    !@brief compute the long-range part of the electron-phonon potential in reciprocal space
    subroutine dvdu_long_range_g(latt_cur,vqpt,natoms,born_eff_charges,quadrupoles,dielectric_tensor,&
                                 atomic_positions,rcut,ewald_parameter,ngvecs,gvecs,gmask,v_long_range_g)
    !>   v_long_range_g(g,alpha,atom) = (4pi/Omega) e^2/(4pi epsilon_0)
    !>                            * (i (q+g)Z_atom + Q )_(alpha, atom)/ ((q+g).epsilon.(q+g))
    !>                            * exp(-i(q+g).R_atom)
    !>                            * exp(-(q+g).epsilon.(q+g)/(4  ewald_parameter^2))
    !>
    !>  set to zero if 1) exp(-(q+g).epsilon.(q+g)/(4  ewald_parameter^2)) < exp(supression_threshold)
    !>                 2) q+G < 10^(-10)
    !>  gmask is set to true if 1) and 2) are fulfiled
    !>
    !> FELEC=e^2/(4pi epsilon_0) in eV AA
        use constant
        use poscar_struct_def
        type (latt) :: latt_cur
        real(q) :: vqpt(3) !> reduced coordinates of the q-point wrt to latt_cur%b
        integer :: natoms
        real(q) :: rcut !< radial coulomb truncation
        real(q) :: ewald_parameter
        real(q) :: born_eff_charges(3,3,natoms) !> born charges Z_atom
        real(q) :: quadrupoles(3,3,3,natoms)
        real(q) :: dielectric_tensor(3,3)     ! epsilon
        real(q) :: atomic_positions(3,natoms) !> reduced coordinates of atoms wrt to latt_cur%a, R_atom
        integer :: ngvecs !> number of g-vectors for wich to compute the long-range part
        integer :: gvecs(:,:)  !> reduced coordinates of the g vectors wrt to 2pi*latt_cur%b. dimension (ngvecs) 
        logical :: gmask(:) !> dimension (ngvecs)
        complex(q) :: v_long_range_g(:,:,:) !> dimension (ngvecs,3,natoms)
        ! Local variables
        logical :: has_born_eff_charges
        logical :: has_quadrupoles
        logical :: has_dielectric_tensor
        integer :: ind, iatom, idir, jdir, kdir
        real(q) :: prefactor, born_factor(3,natoms), dielectric_factor, quadrupole_factor, exp_supression
        real(q) :: qpgvec(3),qpgvec_car(3)
        real(q) :: angle,arg
        real(q) :: truncation
        complex(q) :: phase

        has_born_eff_charges = any(born_eff_charges/=0)
        has_quadrupoles = any(quadrupoles/=0)
        has_dielectric_tensor = .not.is_identity(3,dielectric_tensor)
        v_long_range_g = cmplx(0.0_q,0.0_q,q)
        gmask = .false.
        if (.not.has_born_eff_charges.and..not.has_quadrupoles.and..not.has_dielectric_tensor) return

        prefactor = 4._q * pi / latt_cur%omega * felect

        ! loop over sphere of G vectors
        do ind=1,ngvecs
            qpgvec = vqpt+gvecs(:,ind)

            ! Only take non-zero q + G.
            if (all(abs(qpgvec) < 1e-10_q)) cycle

            qpgvec_car = tpi*matmul(latt_cur%b,qpgvec)
            dielectric_factor = 1._q / dot_product(qpgvec_car,matmul(dielectric_tensor,qpgvec_car))

            ! Beyond a certain lower limit there is no need to account for these anymore
            arg = - 1._q / (dielectric_factor * 4._q * ewald_parameter * ewald_parameter)

            ! ASK_HM : for me this is strange. We fix ewald_parameter using
            ! g_max <=> ECUTLR. Here why don't we just say "if g < g_max cycle"
            if (arg<supression_threshold) cycle
            exp_supression = exp(arg)
            gmask(ind) = .true.

            ! Coulomb truncation
            truncation = 1.0_q
            if (rcut>0) truncation = (1.0_q-cos(norm2(qpgvec_car)*rcut))

            do iatom=1,natoms
                born_factor(:,iatom) = matmul(qpgvec_car,born_eff_charges(:,:,iatom))
            enddo

            ! We are here so we compute the long-range contribution
            do iatom=1,natoms
                angle = dot_product(qpgvec,atomic_positions(:,iatom))
                phase = exp(-citpi*angle)
                do idir=1,3
                    ! due to quadrupoles
                    quadrupole_factor = 0
                    if (has_quadrupoles) then
                        do jdir=1,3
                        do kdir=1,3
                            quadrupole_factor = quadrupole_factor+&
                                 quadrupoles(kdir,jdir,idir,iatom)*qpgvec_car(kdir)*qpgvec_car(jdir)
                        enddo
                        enddo
                    endif
                    ! total contribution
                    v_long_range_g(ind,idir,iatom) = prefactor*dielectric_factor*exp_supression*&
                             (cmplx(0_q,1_q)*born_factor(idir,iatom)+0.5_q*quadrupole_factor)*phase*truncation
               enddo
            enddo
        enddo

    end subroutine dvdu_long_range_g

    !> This is a wrapper to the functionality implemented in pead.F
    subroutine getqijb(b,p,natoms,lmdim,ityp,latt_cur,qijb)
      USE asa
      USE poscar
      USE pseudo
      USE lattice
      USE constant
      USE wave_high
      USE radial, ONLY: SPHERICAL_BESSEL_TRANSFORM_R2Q
      TYPE(potcar),TARGET :: P(:)
      TYPE(latt) :: LATT_CUR
      INTEGER :: NATOMS
      INTEGER :: LMDIM
      INTEGER :: ITYP(NATOMS)
      REAL(q) :: B(3)
      COMPLEX(q) :: QIJB(:,:,:)
      ! local variables
      TYPE (potcar), POINTER :: PP
      INTEGER NI,NT
      INTEGER Ch1,Ch2
      INTEGER LMAX,LMMAX
      INTEGER L,LL,LLP,LM,LMP,M,MP
      INTEGER LDEP_INDEX,LMINDX,ISTART,IEND,IC
      REAL(q) BX,BY,BZ,BB
      REAL(q) XYZPTS(1,3)

      REAL(q), ALLOCATABLE :: YLM(:,:)
      REAL(q), ALLOCATABLE :: FR(:),FG(:,:),G(:)

      COMPLEX(q) CSET,CFAK
      COMPLEX(q) QLLMM(LMDIM,LMDIM)

      QIJB=CMPLX(0.0_q,0.0_q,q)

      BX=(B(1)*LATT_CUR%B(1,1)+B(2)*LATT_CUR%B(1,2)+B(3)*LATT_CUR%B(1,3))*TPI
      BY=(B(1)*LATT_CUR%B(2,1)+B(2)*LATT_CUR%B(2,2)+B(3)*LATT_CUR%B(2,3))*TPI
      BZ=(B(1)*LATT_CUR%B(3,1)+B(2)*LATT_CUR%B(3,2)+B(3)*LATT_CUR%B(3,3))*TPI

      BB=MAX(SQRT(BX*BX+BY*BY+BZ*BZ),1E-10_q)

      LMAX=0
      DO NT=1,SIZE(P)
         LMAX=MAX(LMAX,MAXVAL(P(NT)%LPS))
      ENDDO
      LMAX=2*LMAX
      LMMAX=(LMAX+1)**2

      ALLOCATE(YLM(1,LMMAX))
      XYZPTS(1,1)=BX/BB; XYZPTS(1,2)=BY/BB; XYZPTS(1,3)=BZ/BB
      CALL SETYLM(LMAX,1,YLM,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))

      CSET=CMPLX(0._q,1.0_q,q)

      ion: DO NI=1,NATOMS

    ! Strictly speaking we should account for e^(-i b \cdot R_ion)
    ! (see Eq. B.7 of J. Phys.:Condens. Matter 19, 036215 (2007))
    !      CFAK=CITPI*(ATOMIC_POSITIONS(1,NI)*B(1)+&
    !                  ATOMIC_POSITIONS(2,NI)*B(2)+&
    !                  ATOMIC_POSITIONS(3,NI)*B(3))
    !      CFAK=EXP(CFAK)
    ! but in VASP this factor should not be added
    !(see remarks below Eq 10.35 of the thesis of gK)
          CFAK=CMPLX(1._q,0._q,q)

          NT=ITYP(NI); PP=>P(NT)
          IF (PP%PSDMAX==0) CYCLE ion

          ALLOCATE(FR(PP%R%NMAX),G(1),FG(1,LMAX+1))

          G(1)=BB

          QLLMM=0

          LDEP_INDEX=1
          ! loop over all channels
          LM=1
          chan1: DO Ch1=1,PP%LMAX
              LMP=LM
              chan2: DO Ch2=Ch1,PP%LMAX

                  ! angular moments l and lp of these channels
                  LL=PP%LPS(Ch1)
                  LLP=PP%LPS(Ch2)

                  CALL YLM3LOOKUP(LL,LLP,LMINDX)

                  FR(:)=(PP%WAE(:,Ch1)*PP%WAE(:,Ch2)-PP%WPS(:,Ch1)*PP%WPS(:,Ch2))/PP%R%R(:)/PP%R%R(:)

                  DO L=0,LL+LLP
                     CALL SPHERICAL_BESSEL_TRANSFORM_R2Q(L,PP%R,FR,G(:),FG(:,L+1))
                  ENDDO

                  DO M =1,2*LL+1
                  DO MP=1,2*LLP+1
                     LMINDX=LMINDX+1

                     ISTART=INDCG(LMINDX)
                     IEND  =INDCG(LMINDX+1)

                     DO IC=ISTART,IEND-1
                        QLLMM(LM+M-1,LMP+MP-1)=&
                        QLLMM(LM+M-1,LMP+MP-1)+ &
                        CFAK*(CSET**REAL(JL(IC),q))*YLM3(IC)*YLM(1,JS(IC))*FG(1,JL(IC)+1)
                     ENDDO
                  ENDDO
                  ENDDO

                  ! fill lower triangle
                  IF (Ch1/=Ch2) THEN
                     DO M =1,2*LL+1
                     DO MP=1,2*LLP+1
                        QLLMM(LMP+MP-1,LM+M-1)=QLLMM(LM+M-1,LMP+MP-1)
                     ENDDO
                     ENDDO
                  ENDIF
                  LMP=LMP+2*LLP+1
              ENDDO chan2
              LM =LM +2*LL +1
          ENDDO chan1
          QIJB(:,:,NI)=QLLMM(:,:)

          DEALLOCATE(FR,G,FG)

      ENDDO ion
      DEALLOCATE(YLM)
    end subroutine getqijb

    !> This is a wrapper to the functionality implemented in pead.F
    subroutine getqijb_ityp(ngvecs,gvecs,ntypes,p,lmdim,latt_cur,qijb)
      USE poscar
      USE pseudo, ONLY: potcar
      USE lattice
      USE constant
      USE asa, ONLY: indcg, jl, js, ylm3, setylm, ylm3lookup
      USE radial, ONLY: SPHERICAL_BESSEL_TRANSFORM_R2Q
      INTEGER,INTENT(IN) :: NTYPES
      TYPE(potcar),TARGET :: P(NTYPES)
      TYPE(latt),INTENT(IN) :: LATT_CUR
      INTEGER,INTENT(IN) :: LMDIM
      INTEGER,INTENT(IN) :: NGVECS
      REAL(q),INTENT(IN) :: GVECS(3,NGVECS)
      COMPLEX(q) :: QIJB(NGVECS,LMDIM,LMDIM,NTYPES)
      ! local variables
      TYPE (potcar), POINTER :: PP
      INTEGER IG
      INTEGER NT
      INTEGER Ch1,Ch2
      INTEGER LMAX,LMMAX
      INTEGER L,LL,LLP,LM,LMP,M,MP
      INTEGER LMINDX,ISTART,IEND,IC
      REAL(q) BX,BY,BZ
      REAL(q) B(3)
      REAL(q) BNORM(NGVECS)
      REAL(q) XYZPTS(NGVECS,3)

      REAL(q), ALLOCATABLE :: YLM(:,:)
      REAL(q), ALLOCATABLE :: FR(:),FG(:,:)

      COMPLEX(q),PARAMETER :: CSET=CMPLX(0._q,1._q,q)
      COMPLEX(q) :: QIJB_ELEMENTS(NGVECS)

      QIJB=CMPLX(0.0_q,0.0_q,q)

      ! xxx G_i/|G| if |G|> 1E-10_q
      ! xxx if not what is it ? shoundn't we put XYZPTS to zero ?
      ! xxx WRONG(?)

      ! Loop over g-vectors
      do ig=1,ngvecs
          B=GVECS(:,ig)

          BX=(B(1)*LATT_CUR%B(1,1)+B(2)*LATT_CUR%B(1,2)+B(3)*LATT_CUR%B(1,3))*TPI
          BY=(B(1)*LATT_CUR%B(2,1)+B(2)*LATT_CUR%B(2,2)+B(3)*LATT_CUR%B(2,3))*TPI
          BZ=(B(1)*LATT_CUR%B(3,1)+B(2)*LATT_CUR%B(3,2)+B(3)*LATT_CUR%B(3,3))*TPI

          BNORM(ig)=MAX(SQRT(BX*BX+BY*BY+BZ*BZ),1E-10_q)

          XYZPTS(ig,1)=BX/BNORM(ig)
          XYZPTS(ig,2)=BY/BNORM(ig)
          XYZPTS(ig,3)=BZ/BNORM(ig)
      end do ! g-vectors

      ! Compute spherical harmonics
      LMAX=0
      DO NT=1,SIZE(P)
         LMAX=MAX(LMAX,MAXVAL(P(NT)%LPS))
      ENDDO
      LMAX=2*LMAX
      LMMAX=(LMAX+1)**2

      ! xxx compute  real spherical harmonics Slm(\hat{G})
      ALLOCATE(YLM(NGVECS,LMMAX))
      CALL SETYLM(LMAX,NGVECS,YLM,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))

      ! Loop over atom types and compute
      ! xxx compute what ?
      ALLOCATE(FG(NGVECS,LMAX+1))
      DO NT=1,NTYPES
          PP=>P(NT)
          IF (PP%PSDMAX==0) CYCLE
          ALLOCATE(FR(PP%R%NMAX))

          ! loop over all channels
          LM=1
          chan1: DO Ch1=1,PP%LMAX
              LMP=LM
              chan2: DO Ch2=Ch1,PP%LMAX

                  ! angular moments l and lp of these channels
                  LL=PP%LPS(Ch1)
                  LLP=PP%LPS(Ch2)

                  CALL YLM3LOOKUP(LL,LLP,LMINDX)

                  ! xxx  [R_1(r)R_1(r)-\tilde{R}_1(r)\tilde{R}_2(r)]
                  FR(:)=(PP%WAE(:,Ch1)*PP%WAE(:,Ch2)-PP%WPS(:,Ch1)*PP%WPS(:,Ch2))/PP%R%R(:)/PP%R%R(:)

                  ! Compute radial integral
                  DO L=0,LL+LLP
                     ! xxx 4pi *int dr r^2 j_L(|G|r) [R_1(r)R_1(r)-\tilde{R}_1(r)\tilde{R}_2(r)]
                     CALL SPHERICAL_BESSEL_TRANSFORM_R2Q(L,PP%R,FR,BNORM(:),FG(:,L+1))
                  ENDDO

                  ! Compute angular part
                  DO M =1,2*LL+1
                  DO MP=1,2*LLP+1
                     LMINDX=LMINDX+1

                     ISTART=INDCG(LMINDX)
                     IEND  =INDCG(LMINDX+1)

                     QIJB_ELEMENTS = CMPLX(0.0_q,0.0_q,q)

                    ! xxx sum i^L ....
                    ! xxx JL and YLM3 are know ?
                     DO IC=ISTART,IEND-1
                        QIJB_ELEMENTS(:)=QIJB_ELEMENTS(:)+ &
                        (CSET**REAL(JL(IC),q))*YLM3(IC)*YLM(:,JS(IC))*FG(:,JL(IC)+1)
                     ENDDO

                     QIJB(:,LM+M-1,LMP+MP-1,NT)=QIJB_ELEMENTS
                  ENDDO
                  ENDDO

                  ! fill lower triangle
                  IF (Ch1/=Ch2) THEN
                     DO M =1,2*LL+1
                     DO MP=1,2*LLP+1
                        QIJB(:,LMP+MP-1,LM+M-1,NT)=QIJB(:,LM+M-1,LMP+MP-1,NT)
                     ENDDO
                     ENDDO
                  ENDIF
                  LMP=LMP+2*LLP+1
              ENDDO chan2
              LM =LM +2*LL +1
          ENDDO chan1

          DEALLOCATE(FR)
      ENDDO ! loop atom types
      DEALLOCATE(FG)
      DEALLOCATE(YLM)
    end subroutine getqijb_ityp

    !> @brief get a sphere of G-vectors within the cutoff
    subroutine get_gsphere(grid,b,enmax,vkpt,ngvecs,gvecs,nindpw)
        use constant, only: hsqdtm, tpi
        use mgrid_struct_def
        type (grid_3d),intent(in) :: grid
        real(q),intent(in) :: b(3,3) !< reciprocal lattice vectors
        real(q),intent(in) :: enmax !< cutoff energy in eV
        real(q),intent(in) :: vkpt(3) !< generate the g-phere centered around this k-point
        integer,allocatable,intent(out) :: gvecs(:,:)
        integer,allocatable,intent(out) :: nindpw(:)
        integer,intent(out) :: ngvecs
        call get_gsphere_ngvecs(grid,b,enmax,vkpt,ngvecs)
        call get_gsphere_gvecs_nindpw(grid,b,enmax,vkpt,ngvecs,gvecs,nindpw)
    end subroutine get_gsphere

    !> Get number of g-vectors inside a sphere with enmax cutoff
    subroutine get_gsphere_ngvecs(grid,b,enmax,vkpt,ngvecs)
        use constant, only: hsqdtm, tpi
        use mgrid_struct_def
        type (grid_3d),intent(in) :: grid
        real(q),intent(in) :: b(3,3) !< reciprocal lattice vectors
        real(q),intent(in) :: enmax !< cutoff energy in eV
        real(q),intent(in) :: vkpt(3) !< ASK_HM generate the g-phere centered around this k-point
        integer,intent(out) :: ngvecs
        !local variables
        integer :: maxgvecs, stepgvecs
        integer :: nlboxi
        integer :: n1,n2,n3,nc
        real(q) :: gx,gy,gz
        real(q) :: g1,g2,g3
        integer :: in_this_row
        real(q) :: energy
        integer,allocatable :: gvecs_tmp(:,:)
        integer,allocatable :: nindpw_tmp(:)

        ngvecs=0
        col: do nc=1,grid%rc%ncol
           n2=grid%rc%i2(nc) ; g2=(grid%lpcty(n2)+vkpt(2))
           n3=grid%rc%i3(nc) ; g3=(grid%lpctz(n3)+vkpt(3))

           row: do n1=1,grid%rc%nrow
              g1=(grid%lpctx(n1)+vkpt(1))
              gx= (g1*b(1,1)+g2*b(1,2)+g3*b(1,3)) *tpi
              gy= (g1*b(2,1)+g2*b(2,2)+g3*b(2,3)) *tpi
              gz= (g1*b(3,1)+g2*b(3,2)+g3*b(3,3)) *tpi
              energy=hsqdtm*((gx**2)+(gy**2)+(gz**2))
              if(energy<enmax) then
                 ngvecs=ngvecs+1
              endif
           enddo row
        enddo col
    end subroutine get_gsphere_ngvecs

    !> @brief get a sphere of G-vectors within the cutoff
    subroutine get_gsphere_gvecs_nindpw(grid,b,enmax,vkpt,ngvecs,gvecs,nindpw)
        use constant, only: hsqdtm, tpi
        use mgrid_struct_def
        use tutor, only: vtutor
        type (grid_3d),intent(in) :: grid
        real(q),intent(in) :: b(3,3) !< reciprocal lattice vectors
        real(q),intent(in) :: enmax !< cutoff energy in eV
        real(q),intent(in) :: vkpt(3) !< generate the g-phere centered around this k-point
        integer,intent(in) :: ngvecs
        integer,allocatable,intent(out) :: gvecs(:,:)
        integer,allocatable,intent(out) :: nindpw(:)
        !local variables
        integer :: nlboxi
        integer :: n1,n2,n3,nc
        real(q) :: gx,gy,gz
        real(q) :: g1,g2,g3
        integer :: in_this_row
        integer :: ind
        real(q) :: energy

        ! initial guess for the maximum size
        allocate(gvecs(3,ngvecs))
        allocate(nindpw(ngvecs))

        nlboxi=0
        ind=0
        col: do nc=1,grid%rc%ncol
           n2=grid%rc%i2(nc) ; g2=(grid%lpcty(n2)+vkpt(2))
           n3=grid%rc%i3(nc) ; g3=(grid%lpctz(n3)+vkpt(3))
           in_this_row=0

           row: do n1=1,grid%rc%nrow
              nlboxi=nlboxi+1

              g1=(grid%lpctx(n1)+vkpt(1))
              gx= (g1*b(1,1)+g2*b(1,2)+g3*b(1,3)) *tpi
              gy= (g1*b(2,1)+g2*b(2,2)+g3*b(2,3)) *tpi
              gz= (g1*b(3,1)+g2*b(3,2)+g3*b(3,3)) *tpi
              energy=hsqdtm*((gx**2)+(gy**2)+(gz**2))
              if(energy<enmax) then
                 in_this_row=in_this_row+1

                 ind=ind+1
                 ! if number of g-vectors is larger than storage stop the code
                 if (ind>ngvecs) then
                    call vtutor%bug('Index of g vector larger than storage space',__FILE__,__LINE__)
                 endif
                 gvecs(1,ind)=grid%lpctx(n1)
                 gvecs(2,ind)=grid%lpcty(n2)
                 gvecs(3,ind)=grid%lpctz(n3)
                 nindpw(ind)=nlboxi
              endif
           enddo row
        enddo col
    end subroutine get_gsphere_gvecs_nindpw

    !> intialize the information on a given grid3d datastructure
    subroutine grid3d_init(self,ngx,ngy,ngz,comm,comm_kin,comm_kinter,real2cplx,lreal,lparallel)
        use mgrid_struct_def, only: grid_3d
        use mgrid, only: inilgrd, rec_stdlay, real_stdlay
        type(grid_3d),target :: self
        integer,intent(in) :: ngx,ngy,ngz
        type(communic),target,intent(in) :: comm
        type(communic),target,intent(in) :: comm_kin
        type(communic),target,intent(in) :: comm_kinter
        logical,intent(in) :: real2cplx
        logical,intent(in) :: lreal
        logical,intent(in) :: lparallel
        call inilgrd(ngx,ngy,ngz,self)
        allocate(self%comm,self%comm_kin,self%comm_kinter)
        self%comm = comm
        self%comm_kin = comm_kin
        self%comm_kinter = comm_kinter
        if (real2cplx) then
            self%ngx_rd = ngx/2
            self%ngy_rd = ngy/2
            self%ngz_rd = ngz/2
        else
            self%ngx_rd = ngx
            self%ngy_rd = ngy
            self%ngz_rd = ngz
        endif
        self%real2cplx=real2cplx
        self%lreal=lreal
        ! initialize layouts
        self%rc%index=>null()
        self%in%i2=>null()
        self%in%i3=>null()
        self%in%index=>null()
        call rec_stdlay(self)
        call real_stdlay(self,lparallel=lparallel)
        self%rl_fft=>self%rl
    end subroutine grid3d_init

    !> @brief this routine scales the eigenvectors with the caracteristic length of the oscillator
    !> If we have negative phonon energies ierr becomes 1.
    subroutine get_phonons_zero_disp(natoms_pc,phonon_masses,w_q,e_q,ierr)
        use constant
        integer,intent(in) :: natoms_pc !< number of atoms in the primitive cell
        real(q), intent(in) :: phonon_masses(:) !< phonon masses in atomic units
        real(q), intent(in) ::  w_q(:) !< phonon frequencies
        complex(q), intent(inout) ::  e_q(:,:) !< phonon eigenvectors (modfied uppon exit)
        integer, intent(out) :: ierr
        ! local variables
        integer i0,i1,i2,ind
        real(q), parameter :: hbar = 1.0545718e-34_q  ! js
        real(q), parameter :: fac_to_ang= 10000.0_q*sqrt(hbar/amtokg)

        ierr = 0
        do i0=1,3*natoms_pc
           if (w_q(i0) > 0.0) then
              ind=0
              do i1=1,natoms_pc
                 do i2=1,3
                    ind=ind+1
                    e_q(ind,i0) = e_q(ind,i0)*fac_to_ang/sqrt(2.0_q*phonon_masses(i1)*w_q(i0)) !in ang
                 enddo
              enddo
           else
              ierr = 1
              e_q(:,i0)=cmplx(0.0_q,0.0_q,q)
           endif
        enddo
    end subroutine get_phonons_zero_disp

    !> check if a matrix is identity
    logical function is_identity(n,mat)
        integer,intent(in) :: n
        real(q),intent(in) :: mat(n,n)
        ! local variables
        integer :: i
        real(q) :: identity(n,n)
        identity=0.0_q
        do i=1,3
            identity(i,i)=1.0_q
        enddo
        is_identity = all(abs(identity - mat) < 1e-8)
    end function is_identity

    !> Compute Smith Normal form decomposition of a matrix
    subroutine snf_init(self,s)
        use mkpoints, only: SmithNormalForm
        use mathtools, only: det_3x3, matrix_inverse_3x3
        type(snf) :: self
        integer :: s(3,3)
        ! local variables
        integer :: d(3,3),i
        integer check, diag

        ! check if s is diagonal, then don't do SNF
        diag=0
        check=abs(s(1,2))+abs(s(1,3))+abs(s(2,1))+abs(s(2,3))+abs(s(3,1))+abs(s(3,2))
        if (check .ne. 0) diag=1
        if (s(1,1) .eq. 0) diag=1
        if (s(2,2) .eq. 0) diag=1
        if (s(3,3) .eq. 0) diag=1
        if (diag .eq. 0) then
            self%qm=0
            self%pm=0
            self%qinv=0
            self%pinv=0
            do i=1,3
                self%d(i) = s(i,i)
                self%qm(i,i)=1
                self%pm(i,i)=1
                self%qinv(i,i)=1
                self%pinv(i,i)=1
            enddo
        else
            ! Decompose the matrix S into S = PM D QM
            call SmithNormalForm(s,self%pm,d,self%qm)
            ! TODO: these should use an integer version of det_3x3
            self%qm=self%qm*det_3x3(real(self%qm,q)) ! to be sure det =1
            self%pm=self%pm*det_3x3(real(self%pm,q))

            self%d(1) = d(1,1)
            self%d(2) = d(2,2)
            self%d(3) = d(3,3)

            self%qinv = nint(matrix_inverse_3x3(real(self%qm,q)))
            self%pinv = nint(matrix_inverse_3x3(real(self%pm,q)))
        endif
    end subroutine snf_init

    !> From a set of regular points (in fractional coordinates) determine a mesh that contains them
    !> compute the smallest grid spacing in each direction
    !> grid_point are sp grid point in [0,1[
    !> mesh =nint ( [1/min(grid_point(1)) 1/min(grid_point(2)) 1/min(grid_point(3))])
    subroutine estimate_mesh(ngridpoints,grid_point,mesh,comm)
        integer,intent(in) :: ngridpoints
        real(q),intent(in) :: grid_point(:,:)
        integer,intent(out) :: mesh(:)
        type(communic) :: comm
        ! local variables
        integer :: i
        real(q) :: x,y,z
        real(q) :: small_grid_point(3)

        small_grid_point = [1e8,1e8,1e8]
        ! store smallest vector larger than zero along each direction
        do i=1,ngridpoints
            x = abs(mod(grid_point(1,i),1.0_q))
            y = abs(mod(grid_point(2,i),1.0_q))
            z = abs(mod(grid_point(3,i),1.0_q))
            if (x>1e-8.and.x<small_grid_point(1)) small_grid_point(1)=x
            if (y>1e-8.and.y<small_grid_point(2)) small_grid_point(2)=y
            if (z>1e-8.and.z<small_grid_point(3)) small_grid_point(3)=z
        enddo
        CALLMPI(M_min_d(comm,small_grid_point,3))
        ! compute mesh by inverting
        mesh = nint(1.0_q/small_grid_point)
    end subroutine estimate_mesh

    !> This subroutine return the grid grid_sc defined by
    !> self%fft_mesh_sc(1),self%fft_mesh_sc(2),self%fft_mesh_sc(3) along self%latt_sp%a
    !> as well as the corresponding g vectors in a sphere hbar^2(vqpt+g)^2/2m < self%<encutlr
    !> It also returns the mapping array  map_phelel2vasp
    !> between this grid and the one by phelel in self%grid_point
    subroutine get_grid_in_supercell(self,comm,vqpt,grid_sc,ngvecs,gvecs,nindpw,map_phelel2vasp)
        use tutor
        type (elphon_pot) :: self
        type(communic),target,intent(in) :: comm
        real(q), intent(in) ::  vqpt(3)
        type(grid_3d), intent(out) :: grid_sc
        integer, intent(out):: ngvecs
        integer, allocatable, intent(out) :: gvecs(:,:)
        integer, allocatable, intent(out) :: nindpw(:)
        integer, allocatable, optional, intent(out) :: map_phelel2vasp(:)
        ! localm variables
        integer :: identity(3,3)
        !> create the  grid_sc structure
        !> using the above smallest grid spacing, create a vasp grid, in the
        !> supercell, which contains the phelel grid
        call grid3d_init(grid_sc,self%fft_mesh_sc(1),self%fft_mesh_sc(2),self%fft_mesh_sc(3),&
                    comm,comm,comm,real2cplx=.false.,lreal=.false.,lparallel=.false.)

        ! identity matrix
        identity=0
        identity(1,1)=1
        identity(2,2)=1
        identity(3,3)=1

        !> create mapping arrays between the 2 grids
        !> map_phelel2vasp(iphelel)=ivasp
        !> map_regular2phelel(ivasp)=iphelel, or -1
        if (present(map_phelel2vasp)) then
            call get_mapping(self%ngridpoints_local, self%grid_point, identity, grid_sc, map_phelel2vasp)
            if (any(map_phelel2vasp<1))    call vtutor%error('Failed to create mapping from phelel to vasp')
        endif

        ! ASK_HM : so gvecs are reduced coordinates along the reciprocal vectors
        ! which INCLUDE 2pi (not VASP convention)
        !  gx= (g1*b(1,1)+g2*b(1,2)+g3*b(1,3)) *tpi
        !  gvecs(1,ind)=grid%lpctx(n1)

        !> build sphere of vqpt+G vectors
        !> keep g vectors g=i_1*latt_sp%b_1 + i_2*latt_sp%b_2 + i3*latt_sp%b_3
        !> with i_1=1..NGX=self%fft_mesh_sc(1), ect ect (roughly)
        !> such that hbar^2/2m|q+G|^2<encutlr
        ! ASK_HM: what isnindpw
        call get_gsphere(grid_sc,self%latt_sp%b,self%encutlr,vqpt,ngvecs,gvecs,nindpw)

    end subroutine get_grid_in_supercell

    subroutine get_snf_scp(self)
        type (elphon_pot) :: self
        integer :: N(3,3)

        !> Compute Smith Normal of the microzone to supercell matrix
        !> [A1 A2 A3]=[am1 am2 am3] S
        !> S=N * sp2prim
        N = 0
        N(1,1) = self%fft_mesh(1)
        N(2,2) = self%fft_mesh(2)
        N(3,3) = self%fft_mesh(3)

        !> S = self%snf_fft%pm * diag(self%snf_fft%d) * self%snf_fft%qm
        call snf_init(self%snf_fft,matmul(N,nint(self%sp2prim)))

        !> set fft mesh of the transformed supercell
        self%fft_mesh_scp = self%snf_fft%d
    end subroutine get_snf_scp

    subroutine get_fft_bounds(ng,lattice_pc,lattice_sp)
        use mathtools, only: matrix_inverse_3x3
        ! ng(:) : grid of the primitive cell
        ! ng(:) : grid of the snf surpercell
        integer, intent(in) :: ng(3)
        real(q), intent(in) :: lattice_pc(3,3)
        real(q), intent(in) :: lattice_sp(3,3)
        type(snf) :: snf_fft
        integer :: S(3,3)
        real(q) :: am(3,3),M(3,3)

        am(:,1)=lattice_pc(:,1)/ng(1)
        am(:,2)=lattice_pc(:,2)/ng(2)
        am(:,3)=lattice_pc(:,3)/ng(3)
        M=matrix_inverse_3x3(am)
        S=nint(matmul(lattice_sp,M))
        call snf_init(snf_fft,S)
    end subroutine get_fft_bounds

    !> Create a dummy mpi communicator
    !> This routine should be moved somewhere else where it can be used by more routines
    subroutine comm_init_dummy(self)
        type(communic) ::  self
#ifdef MPI
        self%mpi_comm = MPI_COMM_SELF
#endif
        self%ncpu=1
        self%node_me=1
        self%ionode=1
    end subroutine comm_init_dummy

end module elphon_potential
