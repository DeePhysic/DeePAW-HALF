#include "symbol.inc"

!>
!> @brief This module provides routines to rotate one wavefunction from one k-point to another
!> The code is scavenged from the PEAD routines but here we want a light-weight interface:
!> -> very small memory footprint
!> -> no global variables
!> by Henrique Miranda

module wave_rotate
    use prec
    use constant, only: tpi
    use mgrid, only: grid_3d
    use base, only: symmetry
    use lattice, only: latt
    use full_kpoints, only: mp
    use wave_struct_def, only: wavedes, wavedes1
    use mkpoints, only: kindex_snform, snform, rotation_handle
    use wave_struct_def, only: wavefuna
    use pseudo, only: potcar
    use tutor, only: vtutor
    implicit none

    !> @brief Datastructure with all the information for rotating wavefunctions
    type wave_rotator
        integer :: ispin
        real(q) :: b(3,3) !< vectors of the reciprocal lattice
        real(q) :: a(3,3)
        ! information about the grid
        integer :: ngx, ngy, ngz !< dimensions of the FFT box
        integer,allocatable :: igridind(:,:,:)
        ! information about rotation
        logical :: timerev !< wether to use time-reversal symmetry or not
        integer :: nrotk !< number of space group symmetry operations
        integer :: invmap(48) !< map symmetry operations (IGRPOP) to their inverse
        integer :: isymop(3,3,48) !< real space rotation symmetry operations
        integer :: igrpop(3,3,48) !< reciprocal space symmetry operations
        real(q) :: gtrans(3,48) !< translation part of the symmetry operation
        complex(q) :: rssymop(2,2,48) !< spin rotation matrices
        real(q) :: magrot(48,1) !< magnetic rotation
        integer,allocatable :: rotmap(:,:)
        ! temporary arrays to perform the rotation
        complex(q),allocatable :: cphase(:) !< additional phase for the wave function
        integer,allocatable :: nindpw(:) !< array for mapping G vectors
        type (rotation_handle):: rot_handle
    end type wave_rotator

    !> @brief structure to hold information about the star of a k-point
    !> A star of a k-point is obtained by applying all the possible symmetry
    !> operations to a k-point.
    !> Different symmetries might generate the same k-point in which case
    !> only the first symmetry that generates a new k-point is used
    type kpoint_star
        integer :: nkpts !< number of k-points in the star
        integer :: irot(48) !< index of the rotation (negative index means rotation + inversion symmetry)
        real(q) :: vkpt(3,48) !< coordinates of the kpoints
        real(q) :: wtkpt(48) !< weights of the kpoints (not really used, just to avoid problems with SETWDES)
    end type kpoint_star

    !> @brief Datastructure for going from the IBZ to the FBZ using the k-point stars
    type fbz2ibz_stars
        integer :: nkibz !< number of points in IBZ
        integer :: nkfbz !< number of points in FBZ
        integer,allocatable :: star_size(:) !< size of the star of each ibz point
        integer,allocatable :: irot_fbz2ibz(:) !< index of the rotation (negative index means rotation + inversion symmetry)
        integer,allocatable :: indx_fbz2ibz(:) !< index ffor going from the FBZ to the IBZ
        ! information about rotation
        logical :: timerev !< wether to use time-reversal symmetry or not
        integer :: nrotk !< number of space group symmetry operations
        integer :: isymop(3,3,48) !< real space rotation symmetry operations
        integer :: igrpop(3,3,48) !< reciprocal space symmetry operations
        real(q) :: gtrans(3,48) !< translation part of the symmetry operation
        integer,allocatable :: rotmap(:,:)
        real(q) :: magrot(48,1) !< magnetic rotation
    end type fbz2ibz_stars

    !> A container type for creating wavedes1 at an arbitrary k-point
    type wavedes1_container
        real(q) :: vkpt(3,1)
        real(q) :: wtkpt(1)
        type(wavedes) :: wdes
        type(wavedes1) :: wdes1
    end type wavedes1_container

    contains
    !> @brief Initialization of wave_rotator
    subroutine wave_rotator_init(self,grid,latt_cur,symm,timerev,ispin)
        use spinsym, only: rssymop
        type (wave_rotator) :: self
        type (grid_3d) :: grid
        type (latt) :: latt_cur
        type (symmetry) :: symm
        logical :: timerev
        integer :: ispin
        ! local variables
        integer :: i
        integer :: nions
        self%timerev = timerev
        self%b = latt_cur%b
        self%ispin = ispin
        ! store teh dimensions of the FFT grid
        self%ngx = grid%ngx
        self%ngy = grid%ngy
        self%ngz = grid%ngz
        allocate(self%igridind(self%ngx,self%ngy,self%ngz))
        allocate(self%nindpw(self%ngx*self%ngy*self%ngz))
        allocate(self%cphase(self%ngx*self%ngy*self%ngz))
        ! copy rotmap
        nions = size(symm%rotmap,1)
        allocate(self%rotmap(nions,symm%nrot))
        self%rotmap = symm%rotmap(:,:,1)
        ! copy rotation and translation matrices
        self%nrotk = symm%nrot
        self%isymop = symm%isymop
        self%rssymop(:,:,:) = rssymop(:,:,:)
        self%invmap = symm%invmap
        self%igrpop = symm%igrpop
        self%gtrans = symm%gtrans
        do i=1,size(symm%magrot,1)
          self%magrot(i,1) = symm%magrot(i,1)
        enddo

        self%rot_handle%lsetup = .false.
    end subroutine wave_rotator_init

    !> Apply all the symmetry transformations to a k-point and store the new unique k-points
    !> as well as the symmetry operations used and any translation vector
    subroutine wave_rotator_get_star(self,vkpt,star)
        type (wave_rotator),intent(in) :: self
        real(q),intent(in) :: vkpt(3)
        type (kpoint_star),intent(out) :: star
        call get_star(self%igrpop,self%nrotk,self%timerev,vkpt,star)
    end subroutine wave_rotator_get_star

    subroutine get_star(igrpop,nrotk,timerev,vkpt,star)
        integer,intent(in) :: igrpop(:,:,:)
        integer,intent(in) :: nrotk
        logical,intent(in) :: timerev
        type (kpoint_star),intent(out) :: star
        real(q),intent(in) :: vkpt(3)

        ! local variables
        integer :: nop, idx, indk
        integer :: iop(3,3)
        real(q) :: rop(3,3)
        real(q) :: vkpt_new(3), vkpt_car(3)
        star%nkpts = 0
        ! apply all the symmetry operations to k-point
        do nop=1,nrotk
            rop = igrpop(:,:,nop)
            vkpt_new = rotate_vkpt_rop(vkpt,rop)
            !check if current k-point already exists
            idx = find_coord(star%vkpt(:,:star%nkpts),vkpt_new)
            if (idx==0) then
                star%nkpts = star%nkpts + 1
                star%irot(star%nkpts) = nop
                star%vkpt(:,star%nkpts) = vkpt_new
            endif
        enddo
        ! now with inversion
        if (timerev) then
        do nop=1,nrotk
            ! apply inversion symmetry to form to get IOP
            !call sgrprd(invers,igrpop(:,:,nop),iop)
            !rop = iop
            rop = -1 * igrpop(:,:,nop)
            vkpt_new=rotate_vkpt_rop(vkpt,rop)
            !check if current k-point already exists
            idx = find_coord(star%vkpt(:,:star%nkpts),vkpt_new)
            if (idx==0) then
                star%nkpts = star%nkpts + 1
                star%irot(star%nkpts) = -nop
                star%vkpt(:,star%nkpts) = vkpt_new
            endif
        enddo
        endif

        contains
        function find_coord(vkpts,vkpt) result(idx)
             real(q) :: vkpts(:,:)
             real(q) :: vkpt(3)
             ! local variables
             integer :: idx
             do idx=1,size(vkpts,2)
                 if (all(abs(modulo(vkpts(:,idx)-vkpt+6.5,1._q)-0.5_q)<1e-6)) return
             enddo
             idx = 0
        end function find_coord
    end subroutine get_star

    !> @brief Rotate the wavefunctions from one k-point to another using the index of a rotation
    subroutine wave_rotator_rotate(self,wa_original,wa_rotated,irot,latt_cur,p)
        use kpoints_change, only: generate_rot_handle_low
        use string, only: str
        use wave_high, only: element
        use spinsym, only: rotate_wave_spin_recip_low, rotate_wave_character_spin
        use string, ONLY: str
        use wave_struct_def, ONLY: wavefun1
        type(wave_rotator) :: self
        type(wavefuna),intent(in) :: wa_original
        type(wavefuna) :: wa_rotated
        integer,intent(in) :: irot
        type (latt),intent(in) :: latt_cur
        type (potcar),intent(in) :: p(:)

        PROFILING_START('wave_rotator_rotate')
        call wave_rotator_rotate_pw(self,wa_original,wa_rotated,irot,latt_cur)
        call wave_rotator_rotate_proj(self,wa_original,wa_original%CPROJ,wa_rotated,wa_rotated%CPROJ,irot,latt_cur,p)
        PROFILING_STOP('wave_rotator_rotate')
    end subroutine wave_rotator_rotate

    !> Rotate plane-wave part of the wavefunction
    subroutine wave_rotator_rotate_pw(self,wa_original,wa_rotated,irot,latt_cur)
        use kpoints_change, only: generate_rot_handle_low
        use wave_high, only: element
        use string, ONLY: str
        use wave_struct_def, ONLY: wavefun1
        use spinsym, only: rotate_wave_spin_recip_low
        type(wave_rotator) :: self
        type(wavefuna),intent(in) :: wa_original
        type(wavefuna) :: wa_rotated
        integer,intent(in) :: irot
        type (latt),intent(in) :: latt_cur
        ! local variables
        logical :: linv, lshift
        integer :: ng1,ng2,ng3
        integer :: ng1i,ng2i,ng3i
        integer :: ind
        integer :: i, n, ni, ispinor, abs_irot
        integer :: ngvector_pos_original, ngvector_original, ngvector_tot_original, ngvector_tot_rotated
        integer :: ngvector_pos_rotated, ngvector_rotated
        integer :: rot(3,3)
        integer :: spinflip
        real(q) :: phase
        complex(q),allocatable :: cptwfp(:)
        complex(q),allocatable :: cptwfp_rotated(:)
        ! find which rotation operation we are going to use
        linv   = irot < 0
        abs_irot = abs(irot)
        lshift = any(abs(self%gtrans(:,abs_irot))>1e-6)
        ! create a mapping from g-vectors to the box
        ngvector_tot_original = wa_original%wdes1%ngvector
        ngvector_tot_rotated  = wa_rotated%wdes1%ngvector
        CALLMPI(M_sum_i(wa_original%wdes1%comm_inb,ngvector_tot_original,1))
        CALLMPI(M_sum_i(wa_rotated%wdes1%comm_inb,ngvector_tot_rotated,1))
        if (ngvector_tot_original/=ngvector_tot_rotated) then
            call vtutor%bug('Number of g-vectors is incompatible '//&
                            str(ngvector_tot_original)//' '//&
                            str(ngvector_tot_rotated),__FILE__,__LINE__)
        endif
        self%igridind = 0
        self%nindpw = 0
        self%cphase = cmplx(0_q,0_q)
        ! create a mapping from the g-vectors of g-sphere of the rotated wf to the fft box
        do ni=1,wa_rotated%wdes1%ngvector
            ng1=wa_rotated%wdes1%igx(ni)
            ng2=wa_rotated%wdes1%igy(ni)
            ng3=wa_rotated%wdes1%igz(ni)
            !write(*,*) 'p',ng1,ng2,ng3,ni
            self%igridind(mp(ng1,self%ngx),mp(ng2,self%ngy),mp(ng3,self%ngz))=ni+wa_rotated%wdes1%ngvector_pos-1
        end do
        CALLMPI( M_sum_i(wa_original%wdes1%comm_inb, self%igridind, size(self%igridind)) )
        ! rotate the g-vectors sphere
        if (linv) then
            rot = -1 * self%igrpop(:,:,abs_irot)
        else
            rot = self%igrpop(:,:,abs_irot)
        endif
        ! rotate the original g-sphere and find where the new indexes will be
        do ni=1,wa_original%wdes1%ngvector
            ng1i=wa_original%wdes1%igx(ni)
            ng2i=wa_original%wdes1%igy(ni)
            ng3i=wa_original%wdes1%igz(ni)
            ng1=ng1i*rot(1,1)+ng2i*rot(2,1)+ng3i*rot(3,1)
            ng2=ng1i*rot(1,2)+ng2i*rot(2,2)+ng3i*rot(3,2)
            ng3=ng1i*rot(1,3)+ng2i*rot(2,3)+ng3i*rot(3,3)
            ind = self%igridind(mp(ng1,self%ngx),mp(ng2,self%ngy),mp(ng3,self%ngz))
            self%nindpw(ni+wa_original%wdes1%ngvector_pos-1) = ind
            !TODO: add some more checks
            if (ind==0) then
                call vtutor%error("G vector "//str([ng1i,ng2i,ng3i])//" rotated to "//str([ng1,ng2,ng3])//&
                                  " not found when mapping the wavefunction from k to k'")
            endif
            ! compute phase factor due to rotation
            phase = -tpi*(ng1*self%gtrans(1,abs_irot)+ &
                          ng2*self%gtrans(2,abs_irot)+ &
                          ng3*self%gtrans(3,abs_irot))
            self%cphase(ni+wa_original%wdes1%ngvector_pos-1) = cmplx(cos(phase),sin(phase),kind=q)
        end do
        CALLMPI( M_sum_i(wa_original%wdes1%comm_inb, self%nindpw, ngvector_tot_original) )
        CALLMPI( M_sum_z(wa_original%wdes1%comm_inb, self%cphase, ngvector_tot_original) )
        ! now we are ready to rotate the PW part of the wavefunction
        ngvector_original = wa_original%wdes1%ngvector
        ngvector_pos_original = wa_original%wdes1%ngvector_pos
        ngvector_rotated = wa_rotated%wdes1%ngvector
        ngvector_pos_rotated = wa_rotated%wdes1%ngvector_pos
        ! ncore=1 version
        if (wa_original%wdes1%comm_inb%ncpu==1) then
            do n=1,wa_original%wdes1%nbands
                do ispinor=0,wa_original%wdes1%nrspinors-1
                    ! rotate
                    call rotate_wave( ngvector_original, &
                                      wa_rotated%cptwfp(1+ispinor*ngvector_rotated,n), &
                                      wa_original%cptwfp(1+ispinor*ngvector_original,n), &
                                      self%cphase, self%nindpw, linv, lshift )
                enddo
#ifndef oldsym
                spinflip=0
                if (self%ispin==2) then
                    if (self%magrot(abs_irot,1)==-1) then
                        spinflip = 1
                    endif
                endif
                call rotate_wave_spin_recip_low( wa_rotated%wdes1%lnoncollinear, self%rssymop(:,:,abs_irot), ngvector_rotated, &
                                                 spinflip, wa_rotated%cptwfp(:,n) )
#endif

            enddo
        else
            ! ncore/=1
            allocate(cptwfp(ngvector_tot_original))
            allocate(cptwfp_rotated(ngvector_tot_original))
            if (size(self%cphase,1)<ngvector_tot_original) call vtutor%bug('Wrong dimensions of cphase array',__FILE__,__LINE__)
            if (size(self%nindpw,1)<ngvector_tot_original) call vtutor%bug('Wrong dimensions of nindpw array',__FILE__,__LINE__)
            do n=1,wa_original%wdes1%nbands
                do ispinor=0,wa_original%wdes1%nrspinors-1
                    ! gather wf in one node
                    cptwfp=0.0_q
                    do i=1,ngvector_original
                        cptwfp(i+ngvector_pos_original-1)=wa_original%cptwfp(i+ispinor*ngvector_original,n)
                    enddo
                    CALLMPI( M_sum_z(wa_original%wdes1%comm_inb, cptwfp, size(cptwfp)) )
                    ! now rotate it like in the above
                    call rotate_wave( ngvector_tot_original, &
                                      cptwfp_rotated, &
                                      cptwfp, &
                                      self%cphase, self%nindpw, linv, lshift )
                    ! only keep the components that we will need
                    do i=1,ngvector_rotated
                        wa_rotated%cptwfp(i+ispinor*ngvector_rotated,n) = cptwfp_rotated(i+ngvector_pos_rotated-1)
                    enddo
                enddo
#ifndef oldsym
                spinflip=0
                if (self%ispin==2) then
                    if (self%magrot(abs_irot,1)==-1) then
                        spinflip = 1
                    endif
                endif
                call rotate_wave_spin_recip_low( wa_rotated%wdes1%lnoncollinear, self%rssymop(:,:,abs_irot), ngvector_rotated, &
                                                 spinflip, wa_rotated%cptwfp(:,n) )
#endif
           enddo
        endif
    end subroutine wave_rotator_rotate_pw

    !> Rotate projections part
    subroutine wave_rotator_rotate_proj(self,wa_original,cproj_original,wa_rotated,cproj_rotated,irot,latt_cur,p)
        use kpoints_change, only: generate_rot_handle_low
        use wave_high, only: element, ni_global
        use string, ONLY: str
        use wave_struct_def, ONLY: wavefun1
        use spinsym, only: rotate_wave_character_spin
        use full_kpoints, ONLY: rotate_wave_character_low, rotate_wave_character
        type(wave_rotator) :: self
        type(wavefuna),intent(in) :: wa_original
        GDEF,target,intent(in) :: cproj_original(:,:) !< (nproj,nbands)
        type(wavefuna) :: wa_rotated
        GDEF,target,intent(out) :: cproj_rotated(:,:) !< (nproj,nbands)
        integer,intent(in) :: irot
        type (latt),intent(in) :: latt_cur
        type (potcar),intent(in) :: p(:)
        !local variables
        type(wavefun1) :: w1_rotated
        type(wavefun1) :: w1_original
        logical :: linv
        integer :: n, ni, nis, nip, nt, npro, nprop, ispinor, abs_irot
        integer :: ni_glbl
        ! find which rotation operation we are going to use
        linv   = irot < 0
        abs_irot = abs(irot)
        do n=1,wa_original%wdes1%nbands
            w1_original = element(wa_original,n)
            w1_original%CPROJ=>cproj_original(:,n)
            w1_original%wdes1%nk=abs_irot
            w1_rotated = element(wa_rotated,n)
            w1_rotated%CPROJ=>cproj_rotated(:,n)
            call rotate_wave_character_low(self%rot_handle,p,latt_cur,w1_original%wdes1,w1_original,w1_rotated,&
                                           linv,self%isymop(:,:,abs_irot),self%rotmap(:,abs_irot),self%rssymop(:,:,abs_irot))
        enddo
    end subroutine wave_rotator_rotate_proj

    !> Rotate the projections and their derivatives with respect to atomic positions to a new k-point
    !> TODO: this is not working yet
    subroutine wave_rotator_rotate_dproj(self,wa_original,dproj_original,wa_rotated,dproj_rotated,irot,latt_cur,p)
        use kpoints_change, only: generate_rot_handle_low
        use wave_high, only: element
        use string, ONLY: str
        use wave_struct_def, ONLY: wavefun1
        type(wave_rotator) :: self
        type(wavefuna),intent(in) :: wa_original
        GDEF,intent(in) :: dproj_original(:,:,0:,:)
        type(wavefuna) :: wa_rotated
        GDEF :: dproj_rotated(:,:,0:,:)
        integer,intent(in) :: irot
        type (latt),intent(in) :: latt_cur
        type (potcar),intent(in) :: p(:)
        ! local variables
        integer :: i,j,k,l
        real(q) :: rop(3,3)
        GDEF :: vec(3),vec_rotated(3)
        GDEF,allocatable :: dproj_original_(:,:)
        GDEF,allocatable :: dproj_rotated_(:,:)
        logical :: linv
        integer :: ib,ispinor,iproj,idir,nbands,nproj,nrspinors

        nproj = size(dproj_original,1)
        nbands = size(dproj_original,2)
        nrspinors = size(dproj_original,4)
        !TODO: change this
        allocate(dproj_original_(nproj*nrspinors,nbands))
        allocate(dproj_rotated_(nproj*nrspinors,nbands))
        ! rotate projections
        call pack_to_vasp(0,dproj_original,dproj_original_)
        call wave_rotator_rotate_proj(self,wa_original,dproj_original_,wa_rotated,dproj_rotated_,irot,latt_cur,p)
        call unpack_from_vasp(0,dproj_rotated_,dproj_rotated)
        ! rotate derivatives of projections
        do idir=1,3
            call pack_to_vasp(idir,dproj_original,dproj_original_)
            call wave_rotator_rotate_proj(self,wa_original,dproj_original_,wa_rotated,dproj_rotated_,irot,latt_cur,p)
            call unpack_from_vasp(idir,dproj_rotated_,dproj_rotated)
        enddo
        ! rotate final derivatives of projections inplace
#if 0
        ! taken from KDER_WAVE
        rop=0
        do l=1,3
        do k=1,3
           do j=1,3
           do i=1,3
              rop(l,i)=rop(l,i)+latt_cur%a(l,k)*self%isymop(j,k,abs(irot))*latt_cur%b(i,j)
           enddo
           enddo
        enddo
        enddo
#else
        rop = matmul(latt_cur%a,matmul(transpose(self%isymop(:,:,abs(irot))),transpose(latt_cur%b)))
#endif
        linv = irot < 0
        ! in-place 3D rotation of the derivatives of the projections
        do ispinor=1,nrspinors
           do ib=1,nbands
              do iproj=1,nproj
                vec = dproj_rotated(iproj,ib,1:3,ispinor)
                vec_rotated = matmul(rop,vec)
                dproj_rotated(iproj,ib,1:3,ispinor)=vec_rotated
              enddo
           enddo
        enddo
        contains
        subroutine pack_to_vasp(idir,dproj,dproj_vasp)
            integer :: idir
            GDEF,intent(in) :: dproj(:,:,0:,:)
            GDEF,intent(out) :: dproj_vasp(:,:)
            do ib=1,nbands
               do ispinor=1,nrspinors
                  do iproj=1,nproj
                     dproj_vasp(iproj+(ispinor-1)*nproj,ib) = dproj(iproj,ib,idir,ispinor)
                  end do
               enddo
            enddo
        end subroutine pack_to_vasp

        subroutine unpack_from_vasp(idir,dproj_vasp,dproj)
            integer :: idir
            GDEF,intent(in) :: dproj_vasp(:,:)
            GDEF,intent(out) :: dproj(:,:,0:,:)
            do ib=1,nbands
               do ispinor=1,nrspinors
                  do iproj=1,nproj
                     dproj(iproj,ib,idir,ispinor) = dproj_vasp(iproj+(ispinor-1)*nproj,ib)
                  enddo
               enddo
            enddo
        end subroutine unpack_from_vasp
    end subroutine wave_rotator_rotate_dproj

    !> @brief Free the rotation object
    subroutine wave_rotator_free(self)
        use full_kpoints, only: deallocate_rot_handle
        type(wave_rotator) :: self
        deallocate(self%nindpw)
        deallocate(self%cphase)
        deallocate(self%igridind)
        deallocate(self%rotmap)
        call deallocate_rot_handle(self%rot_handle)
    end subroutine wave_rotator_free

    !> @brief generate fbz2ibz_stars datastructure with information
    !> of how to go from the IBZ to the FBZ using k-point stars
    subroutine fbz2ibz_stars_init(self,symm,timerev,nkibz,vkpt)
       type (fbz2ibz_stars) :: self
       type (symmetry) :: symm !< datastructure for rotating the WFs
       logical :: timerev !< flag governing wether we have time reversal symmetry
       integer :: nkibz !< Number of k-points in the IBZ
       real(q) :: vkpt(3,nkibz) !< Coordinates of the points in the IBZ
       ! local variables
       type (kpoint_star) :: kstar
       integer :: fbz_nkpt_count
       integer :: ik_ibz, ik_fbz, ik_star

       self%nkibz = nkibz
       ! save infromation related to the symmetry transformations
       call fbz2ibz_stars_set_rot(self,symm,timerev)

       ! Compute the size of the full brillouin zone
       fbz_nkpt_count = 0
       do ik_ibz=1,nkibz
           ! get star for k'
           call get_star(self%igrpop,self%nrotk,self%timerev,vkpt(:,ik_ibz),kstar)
           ! increament the size of the fbz
           fbz_nkpt_count = fbz_nkpt_count + kstar%nkpts
       enddo
       self%nkfbz = fbz_nkpt_count

       ! allocate main arrays
       allocate(self%irot_fbz2ibz(self%nkfbz))
       allocate(self%indx_fbz2ibz(self%nkfbz))
       allocate(self%star_size(self%nkibz))

       ! generate all the stars of k
       fbz_nkpt_count = 0
       do ik_ibz=1,nkibz

           ! get star for k'
           call get_star(self%igrpop,self%nrotk,self%timerev,vkpt(:,ik_ibz),kstar)

           ! store the size of the stars
           self%star_size(ik_ibz) = kstar%nkpts

           do ik_star=1,kstar%nkpts
               ik_fbz = fbz_nkpt_count+ik_star
               self%irot_fbz2ibz(ik_fbz) = kstar%irot(ik_star)
               self%indx_fbz2ibz(ik_fbz) = ik_ibz
           enddo
           fbz_nkpt_count = fbz_nkpt_count + kstar%nkpts
       enddo
    end subroutine fbz2ibz_stars_init

    !> save infromation related to the symmetry transformations
    subroutine fbz2ibz_stars_set_rot(self,symm,timerev)
       type (fbz2ibz_stars) :: self
       type (symmetry), intent(in) :: symm !< datastructure for rotating the WFs
       logical, intent(in) :: timerev
       !local
       integer :: i
       self%nrotk = symm%nrot
       self%timerev = timerev
       self%isymop = symm%isymop
       self%igrpop = symm%igrpop
       self%gtrans = symm%gtrans
       allocate(self%rotmap(size(symm%rotmap,1),size(symm%rotmap,2)))
       self%rotmap = symm%rotmap(:,:,1)
       do i=1,size(symm%magrot,1)
         self%magrot(i,1) = symm%magrot(i,1)
       enddo
    end subroutine fbz2ibz_stars_set_rot

    !> get the first start of the k-point with the largest multiplicity
    !> this is useful to estimate the memory usage
    subroutine fbz2ibz_stars_get_large_kstar(self,kstar,vkpt)
       type (fbz2ibz_stars) :: self
       type (kpoint_star) :: kstar
       real(q) :: vkpt(:,:)
       ! local variables
       integer :: ik_ibz, largest_star_size
       ! get size of the largest star
       largest_star_size = maxval(self%star_size)
       ! get index of the first ibz kpoints with the largest star
       do ik_ibz=1,self%nkibz
           if (self%star_size(ik_ibz)==largest_star_size) exit
       enddo
       ! return the star for this k-point
       call get_star(self%igrpop,self%nrotk,self%timerev,vkpt(:,ik_ibz),kstar)
    end subroutine fbz2ibz_stars_get_large_kstar

    !> @brief generate fbz2ibz_stars datastructure with information
    !> of how to go from the IBZ to the FBZ using k-point stars
    !> but in this case no symmetries are used at all
    subroutine fbz2ibz_stars_dummy_init(self,symm,nkibz,vkpt)
       type (fbz2ibz_stars) :: self
       type (symmetry) :: symm !< datastructure for rotating the WFs
       integer,intent(in) :: nkibz !< Number of k-points in the IBZ
       real(q),intent(in) :: vkpt(3,nkibz) !< Coordinates of the points in the IBZ
       ! local variables
       type (kpoint_star) :: kstar
       integer :: fbz_nkpt_count
       integer :: ik_ibz, ik_fbz, ik_star

       self%nkfbz = nkibz
       self%nkibz = nkibz
       ! save infromation related to the symmetry transformations
       call fbz2ibz_stars_set_rot(self,symm,timerev=.false.)
       self%nrotk = 1

       ! allocate main arrays
       allocate(self%irot_fbz2ibz(self%nkfbz))
       allocate(self%indx_fbz2ibz(self%nkfbz))
       allocate(self%star_size(self%nkibz))

       ! generate all the stars of k
       do ik_ibz=1,nkibz
           self%star_size(ik_ibz) = 1
           self%irot_fbz2ibz(ik_ibz) = 1
           self%indx_fbz2ibz(ik_ibz) = ik_ibz
       enddo
    end subroutine fbz2ibz_stars_dummy_init

    !> @brief free the fbz2ibz datastructure
    subroutine fbz2ibz_stars_free(self)
       type (fbz2ibz_stars) :: self
       deallocate(self%irot_fbz2ibz)
       deallocate(self%indx_fbz2ibz)
       deallocate(self%star_size)
    end subroutine fbz2ibz_stars_free

    !> @brief Apply symmetry operations stored in fbz2ibz and return arms of the star.
    subroutine kstar_from_fbz2ibz(kstar,vkpt,fbz2ibz,ikp_ibz)
        type(kpoint_star),intent(inout) :: kstar !< star structure
        real(q),intent(in) :: vkpt(3) !< k-point coordinate for which to get the star
        type(fbz2ibz_stars),intent(in) :: fbz2ibz !< object containing map of ibz to fbz using stars
        integer,intent(in) :: ikp_ibz !< index of the k-point in the irreducible Brillouin zone

        integer :: nop, ik_fbz, ik_star
        real(q) :: rop(3,3)
        real(q) :: vkpt_new(3)

        ik_fbz = sum(fbz2ibz%star_size(1:(ikp_ibz - 1)))

        kstar%nkpts = fbz2ibz%star_size(ikp_ibz)
        do ik_star = 1, kstar%nkpts
            nop = fbz2ibz%irot_fbz2ibz(ik_fbz + ik_star)
            vkpt_new = rotate_vkpt_irot(vkpt,fbz2ibz%igrpop,nop)
            kstar%irot(ik_star) = nop
            kstar%vkpt(:,ik_star) = vkpt_new(:)
        end do
    end subroutine kstar_from_fbz2ibz

    !> Rotate a k-point according to a symmetry operation in reciprocal space
    !> using the index of the rotation in reciprocal space
    !> If the index of the rotation is negative apply time reversal
    function rotate_vkpt_irot(vkpt,igrpop,irot) result(vkpt_rotated)
        real(q),intent(in) :: vkpt(3)
        integer,intent(in) :: igrpop(3,3,48)
        integer,intent(in) :: irot
        real(q) :: vkpt_rotated(3)
        !local variables
        real(q) :: rop(3,3)
        rop(:,:) = igrpop(:,:,abs(irot))
        if (irot < 0) then
            rop(:,:) = -1 * rop(:,:)
        end if
        vkpt_rotated=rotate_vkpt_rop(vkpt,rop)
    end function rotate_vkpt_irot

    !> Rotate a k-point according to a symmetry operation in reciprocal space
    !> using the rotation matrix in reciprocal space
    function rotate_vkpt_rop(vkpt,rop) result(vkpt_rotated)
        real(q),intent(in) :: vkpt(3)
        real(q),intent(in) :: rop(3,3)
        real(q) :: vkpt_rotated(3)
        vkpt_rotated(1)=vkpt(1)*rop(1,1)+vkpt(2)*rop(2,1)+vkpt(3)*rop(3,1)
        vkpt_rotated(2)=vkpt(1)*rop(1,2)+vkpt(2)*rop(2,2)+vkpt(3)*rop(3,2)
        vkpt_rotated(3)=vkpt(1)*rop(1,3)+vkpt(2)*rop(2,3)+vkpt(3)*rop(3,3)
    end function rotate_vkpt_rop

    !> Create a wavefunction descriptor for all the k-points in the star of k
    subroutine kstar_get_wdes(self,wdes,wdes_ref,latt_cur,latt_ini)
        use wave_struct_def, only: wavedes
        use wave, only: allocwdes
        type(kpoint_star),target,intent(inout) :: self !< star structure
        type(wavedes),intent(inout) :: wdes !< output wave descriptor
        type(wavedes),intent(in) :: wdes_ref !< reference wave descriptor
        type(latt),intent(in) :: latt_cur
        type(latt),intent(in) :: latt_ini
        call get_wdes(wdes,wdes_ref,self%nkpts,self%vkpt,self%wtkpt,latt_cur,latt_ini)
    end subroutine kstar_get_wdes

    !> create a wavefunction descriptor for a set of k-points
    subroutine get_wdes(wdes,wdes_ref,nkpts,vkpt,wtkpt,latt_cur,latt_ini)
        use wave_struct_def, only: wavedes
        use wave, only: allocwdes
        type(wavedes),intent(inout) :: wdes !< output wave descriptor
        type(wavedes),intent(in) :: wdes_ref !< reference wave descriptor
        integer,intent(in) :: nkpts
        real(q),target,intent(in) :: vkpt(:,:)
        real(q),target,intent(in) :: wtkpt(:)
        type(latt),intent(in) :: latt_cur
        type(latt),intent(in) :: latt_ini
        wdes = wdes_ref
        wdes%nkdim = nkpts
        wdes%nkpts = nkpts
        wdes%nkpts_for_gen_layout = nkpts
        wdes%vkpt=>vkpt
        wdes%wtkpt=>wtkpt
        call allocwdes(wdes,.true.)
        wdes%nplwkp=0
        wdes%luseinv=.false.
        wdes%nindpw=0
        call gen_index(wdes%grid, wdes, latt_cur%b, latt_ini%b, -1, -1, .true.)
    end subroutine get_wdes

    !> Initialize wavedes1_container from a reference wavefunction descriptor
    !> These simple routines are needed because I need a wavedes1 at an arbitrary k-point.
    !> To set the k-point coordinate I need to call gen_index, but that requires wavedes
    !> so in the end I need to call setwdes to get wavedes1 from wavedes.
    !> the structure also keeps track of the coordinate of the k-point and its weight
    !> to ensure there are no memory leaks (no allocated pointers in this structure)
    subroutine wdes1_container_init(self,wdes_ref,vkpt,wtkpt,latt_cur,latt_ini)
        use wave, only: setwdes
        type(wavedes1_container) :: self
        type(wavedes),intent(in) :: wdes_ref !< reference wave descriptor
        real(q),target,intent(in) :: wtkpt
        real(q),target,intent(in) :: vkpt(3)
        type(latt),intent(in) :: latt_cur
        type(latt),intent(in) :: latt_ini
        self%wtkpt(1) = wtkpt
        self%vkpt(:,1) = vkpt
        ! create wdes with 1 kpoint
        call get_wdes(self%wdes,wdes_ref,1,self%vkpt,self%wtkpt,latt_cur,latt_ini)
        ! point wdes1 to wdes
        call setwdes(self%wdes,self%wdes1,1)
    end subroutine wdes1_container_init

    !> free wavedes1_container
    subroutine wdes1_container_free(self)
        use wave, only: deallocwdes
        type(wavedes1_container) :: self
        ! deallocate wdes
        call deallocwdes(self%wdes,.true.)
    end subroutine wdes1_container_free

end module wave_rotate

!in principle for each k in the IBZ we need to sum over all the k' in the FBZ
!star of k-point -> reduces the amount of symmetries but some might still remain
!use the symmetries in the star of k to obtain a new IBZ
!compute electron-phonon for the the q-points in the star of k
