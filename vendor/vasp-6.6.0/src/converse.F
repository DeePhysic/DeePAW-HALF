#include "symbol.inc"

!***********************************************************************
!
!***********************************************************************

      MODULE mconverse

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: NMR_CONVERSE_MAIN

      CONTAINS

!************************ SUBROUTINE NMR_CONVERSE_MAIN *****************
!
!> Compute the orbital contribution to the magnetisation
!> using the approach of Zwanziger et al.: Eq. 36 of PRB 107, 165157 (2023)
!
!***********************************************************************

      SUBROUTINE NMR_CONVERSE_MAIN( &
         HAMILTONIAN, P, KPOINTS, SYMM, INFO, GRID, GRIDC, GRIDUS, C_TO_US, &
         WDES, W, NONL_S, NONLR_S, IRDMAX, T_INFO, LATT_CUR, LATT_INI, &
         LMDIM, CDIJ, CQIJ, CRHODE, SV, EFERMI, IO)

        USE prec
        USE base
        USE pseudo_struct_def
        USE mgrid_struct_def
        USE hamil_struct_def
        USE poscar_struct_def
        USE mkpoints_struct_def
        USE wave_struct_def
        USE pead
        USE lattice
        USE nonl_high
        USE subrot_cluster
        USE msymmetry
        USE mlr_optic
        USE mlrf_main
        USE mlr_main_nmr

        !> vector potential and meta-gga specific potential terms
        TYPE (ham_handle) :: HAMILTONIAN
        !> PAW datasets
        TYPE(potcar) :: P(:)
        !> k-point information
        TYPE(kpoints_struct) :: KPOINTS
        !> symmetry information
        TYPE(symmetry) :: SYMM
        !> Control variables
        TYPE(info_struct) :: INFO
        !> FFT grid for the orbitals
        TYPE(grid_3d) :: GRID
        !> grids for potentials and charge
        TYPE(grid_3d) :: GRIDC, GRIDUS
        !> map between GRIDC and GRIDUS
        TYPE(transit) :: C_TO_US
        !> wave function descriptor
        TYPE(wavedes), TARGET :: WDES
        !> the groundstate orbitals
        TYPE(wavespin) :: W
        !> real space projectors
        TYPE(nonlr_struct) :: NONLR_S
        !> reciprocal space projectors
        TYPE(nonl_struct) :: NONL_S
        !> definition of the crystalline basis
        TYPE(type_info) :: T_INFO
        !> definition of the unit cells
        TYPE(latt) :: LATT_CUR
        TYPE(latt) :: LATT_INI
        !> io units
        TYPE(in_struct) :: IO
        !> first dimension of CDIJ and CQIJ
        INTEGER :: LMDIM
        !> the PAW strength and overlap parameters
        OVERLAP :: CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
                   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
                 CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

        !> the local potential
        RGRID :: SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)

        !> the Fermi energy
        REAL(q) :: EFERMI

        INTEGER :: IRDMAX

        ! local variables
        TYPE(eigenf_cluster_pointer), POINTER :: DEG_CLUSTER(:,:)
        TYPE(wavespin) :: DWDK

        COMPLEX(q), ALLOCATABLE, TARGET :: RPHI(:,:,:,:,:)
        GDEF, ALLOCATABLE, TARGET :: RPHI_CPROJ(:,:,:,:,:)

        REAL(q) :: MORB(3)

        INTEGER :: IDIR

        LOGICAL :: LDUMMY

        ALLOCATE(RPHI(WDES%NRPLWV,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3), &
                 RPHI_CPROJ(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3))

        RPHI=0
        RPHI_CPROJ=0

        ! compute -i | d u_nk / dk > and -i < p_i | d u_nk / dk >
        IF (LUSEPEAD()) THEN
           LPEAD_RETURN_Q_CPROJ=.FALSE.
           ! ... from finite differences
           CALL PEAD_DPSI_DK_ALL(W, KPOINTS, P, CQIJ, LATT_CUR, T_INFO, RPHI_d=RPHI, RPHI_CPROJ_d=RPHI_CPROJ)
           LPEAD_RETURN_Q_CPROJ=.TRUE.
        ELSE
           NULLIFY(DEG_CLUSTER)
           CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)

           IF (IO%IU0>=0) THEN
              WRITE(IO%IU0,*) 'Linear response G [H, r] |phi>, progress :'
              WRITE(17,*)     'Linear response G [H, r] |phi>, progress :'
           ENDIF
           IF (IO%IU6>=0) &
              WRITE(IO%IU6,*) 'Linear response G [H, r] |phi>, progress :'

           DO IDIR = 1, 3
              IF (IO%IU0>=0) THEN
                 WRITE(IO%IU0,'(A,I3)') '  Direction: ',IDIR
                 WRITE(17,'(A,I3)')     '  Direction: ',IDIR
              ENDIF
              IF (IO%IU6>=0) &
                 WRITE(IO%IU6,'(A,I3)') '  Direction: ',IDIR

              DWDK%CW => RPHI(:,:,:,:,IDIR)
              DWDK%CPROJ => RPHI_CPROJ(:,:,:,:,IDIR)
              DWDK%WDES => WDES

              CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, T_INFO, NONLR_S, NONL_S, &
                 W, DWDK, LMDIM, P, CQIJ, SYMM, IDIR, LDUMMY, IO%IU0, IO%IU6)

              ! ... or from linear response
              CALL LRF_RPHI(HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, &
                 T_INFO, INFO, IO, GRID, GRIDC, GRIDUS, C_TO_US, IRDMAX, &
                 CDIJ, CQIJ, SV, LMDIM, DEG_CLUSTER, IDIR, .TRUE., &
                 RPHI_d=DWDK%CW, RPHI_CPROJ_d=DWDK%CPROJ)
           ENDDO

           CALL FREE_DEG_CLUSTERS(WDES, DEG_CLUSTER)
        ENDIF

        ! the routines above return -i | d u_nk / dk > and -i < p_i | d u_nk /dk >
        RPHI = RPHI * (0._q,1_q)
        RPHI_CPROJ = RPHI_CPROJ * (0._q,1._q)

        ! these need to be allocated, otherwise the code will crash in SETWAV (in the calltree of PROALL below)
        ALLOCATE(DWDK%FERWE(WDES%NBANDS,WDES%NKPTS,WDES%ISPIN), DWDK%CELEN(WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))

        DWDK%WDES => WDES
        DO IDIR = 1, 3
           DWDK%CW => RPHI(:,:,:,:,IDIR)
           DWDK%CPROJ => RPHI_CPROJ(:,:,:,:,IDIR)
           ! project | d u_nk / dk_idir > onto the empty subspace
           CALL PROJ_IN_EMPTY_SUBSPACE_BLAS(W, GRID, NONL_S, NONLR_S, LATT_CUR, LMDIM, CQIJ, INFO, DWDK)
           ! compute < p_i | d u_nk / dk_idir >
           CALL PROALL(GRID, LATT_CUR, NONLR_S, NONL_S, DWDK)
        ENDDO

        DEALLOCATE(DWDK%FERWE, DWDK%CELEN)

        CALL MORBCONV(HAMILTONIAN, P, INFO, GRID, WDES, W, RPHI, RPHI_CPROJ, &
           NONL_S, NONLR_S, T_INFO, LATT_CUR, LMDIM, CDIJ, CQIJ, CRHODE, SV, EFERMI, MORB)

        DEALLOCATE(RPHI, RPHI_CPROJ)

        RETURN
      END SUBROUTINE NMR_CONVERSE_MAIN


!************************ SUBROUTINE MORBCONV **************************
!
!> Compute the orbital contribution to the magnetisation
!> using the approach of Zwanziger et al.: Eq. 36 of PRB 107, 165157 (2023)
!
!***********************************************************************
#define debug
      SUBROUTINE MORBCONV(HAMILTONIAN, P, INFO, GRID, WDES, W, RPHI, RPHI_CPROJ, &
         NONL_S, NONLR_S, T_INFO, LATT_CUR, LMDIM, CDIJ, CQIJ, CRHODE, SV, EFERMI, &
         MORB)

        USE prec
        USE constant
        USE base
        USE pseudo_struct_def
        USE mgrid_struct_def
        USE hamil_struct_def
        USE poscar_struct_def
        USE paw
        USE lattice
        USE wave_high
        USE nonl_high
        USE hamil
        USE morbitalmag
        USE dfast, ONLY : NSTRIP_STANDARD
        USE fock_glb, ONLY: LHFCALC, AEXX
        USE fock, ONLY : FOCK_ACC
        USE mlr_main_nmr, ONLY : EIJK

        !> vector potential and meta-gga specific potential terms
        TYPE (ham_handle) :: HAMILTONIAN
        !> PAW datasets
        TYPE(potcar), TARGET :: P(:)
        !> Control variables
        TYPE(info_struct) :: INFO
        !> FFT grid for the orbitals
        TYPE(grid_3d) :: GRID
        !> wave function descriptor
        TYPE(wavedes) :: WDES
        !> the groundstate orbitals
        TYPE(wavespin) :: W

        !> the derivative of the orbitals w.r.t. bloch vector | d u_nk / dk >
        COMPLEX(q), TARGET :: RPHI(WDES%NRPLWV,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3)
        !> projections < p_i | d u_nk / dk >
        GDEF, TARGET :: RPHI_CPROJ(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3)

        !> real space projectors
        TYPE(nonlr_struct) :: NONLR_S
        !> reciprocal space projectors
        TYPE(nonl_struct) :: NONL_S
        !> definition of the crystalline basis
        TYPE(type_info) :: T_INFO
        !> definition of the unit cells
        TYPE(latt) :: LATT_CUR
        !> first dimension of CDIJ and CQIJ
        INTEGER :: LMDIM
        !> the PAW strength and overlap parameters
        OVERLAP :: CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
                   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
                 CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

        !> the local potential
        RGRID :: SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)

        !> the Fermi energy
        REAL(q) :: EFERMI

        !> the orbital magnetization
        REAL(q) :: MORB(3)

        ! local variables
        TYPE(wavedes1), TARGET :: WDES1
        TYPE(wavefun1) :: W1, WDK(3), WTMP

        TYPE(nonlr_struct) :: NONLR_D(3)
        TYPE(nonl_struct) :: NONL_D(3)

        TYPE(potcar), POINTER :: PP

        GDEF :: OVL1(WDES%NPROD,WDES%NBANDS), OVL2(WDES%NPROD,WDES%NBANDS,3)

        GDEF :: CPROJD(WDES%NPROD,WDES%NBANDS,3), COVL(WDES%NPROD)
        GDEF :: DSDK(WDES%NB_TOT,WDES%NB_TOT,3)
        GDEF, TARGET :: CDD(WDES%NPROD,3,3)

        COMPLEX(q) :: CTMP, CTMP1(3,3), CTMP2(3,3), CTMP3(3,3), CTMP4(3,3)
        COMPLEX(q) :: CHRN1(3,3), CHRN3(3,3)

        REAL(q) :: E0NK, MTOT(3), M1C(3), MP1C(3), MD1C(3), MLQ(3), WEIGHT

        INTEGER :: ISP, NK, NB, NOCC
        INTEGER :: IDIR, JDIR, KDIR
        INTEGER :: NI, NIP, NT

        TYPE(wavefun1), ALLOCATABLE :: WHF1(:)
        TYPE(wavefuna) :: WHFA
        COMPLEX(q) :: CDUMMY
        INTEGER :: NPOS, NSTRIP, NSTRIP_ACT
        LOGICAL :: LHFCALC_LOCAL
#ifdef debug
        COMPLEX(q) :: MTMP1(3), MTMP2(3), MTMP3(3), MTMP4(3), CHRN(3)
        REAL(q) :: MPARA(3,T_INFO%NIONS), MDIA(3)
#endif

        MTOT = 0 ; M1C = 0 ; MP1C = 0 ; MD1C = 0 ; MLQ = 0
        CHRN = 0

        CALL SETWDES(WDES,WDES1,0)

        ! workspace for the FFT of an orbital, and for its projection onto the projectors
        CALL NEWWAV_NOCW(W1, WDES1)
        DO IDIR = 1, 3
           CALL NEWWAV_R(WDK(IDIR), WDES1)
        ENDDO

        DO IDIR = 1, 3
           IF (NONLR_S%LREAL) THEN
              CALL NONLR_ASSIGN(NONLR_D(IDIR), NONLR_S)
              CALL NONLR_ALLOC_CRREXP(NONLR_D(IDIR))
           ELSE
              CALL NONL_ALLOC_DER(NONL_S, NONL_D(IDIR))
              ! returns -i < d p_ik / dk_idir | at all k
              CALL SPHER_DER(GRID, NONL_D(IDIR), P, WDES, LATT_CUR, IDIR)
           ENDIF
        ENDDO

        LHFCALC_LOCAL = LHFCALC .AND. AEXX/=0
        IF (LHFCALC_LOCAL) THEN
           NSTRIP = MIN(NSTRIP_STANDARD, WDES%NBANDS)
           ALLOCATE(WHF1(NSTRIP))
           CALL NEWWAVA(WHFA, WDES1, NSTRIP)
        ENDIF

        CTMP1 = 0 ; CTMP2 = 0 ; CTMP3 = 0 ; CTMP4 = 0
        CHRN1 = 0 ; CHRN3 = 0

!>=============================================================================
!>
!>  MTOT
!>
!>=============================================================================
        spin: DO ISP = 1,WDES%ISPIN
           kpoint: DO NK = 1, WDES%NKPTS
#ifdef MPI
              ! parallelization over k-points
              IF (MOD(NK-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE kpoint
#endif
              CALL SETWDES(WDES, WDES1, NK)

              ! OVL1(i, n) = \sum_j Q_ij < p_jk | u_nk >  for i=1,lmdim and n runs over all bands
              OVL1 = 0 ; CALL OVERL(WDES1, WDES1%LOVERL, LMDIM, CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP), OVL1(1,1))

              OVL2 = 0 ; CPROJD = 0

              DO IDIR = 1, 3
                 ! OVL2(i, n, idir) = \sum_j Q_ij < p_jk | d u_nk / dk_idir >  for i=1,lmdim and n runs over all bands
                 CALL OVERL(WDES1, WDES1%LOVERL, LMDIM, CQIJ(1,1,1,ISP), RPHI_CPROJ(1,1,NK,ISP,IDIR), OVL2(1,1,IDIR))
                 IF (NONLR_S%LREAL) THEN
                    ! returns < p_ik |r for k=NK
                    CALL PHASERR(GRID, LATT_CUR, NONLR_D(IDIR), NK, WDES, IDIR)
                 ELSE
                    NONL_D(IDIR)%NK = -1
                    CALL PHASE(WDES, NONL_D(IDIR), NK)
                 ENDIF

                 ! compute CPROJD(i, n, idir) = < d p_i / d k_idir | u_nk > for all bands at (NK,ISP)
                 DO NB = 1, WDES%NBANDS
                    W1%CW => W%CW(:,NB,NK,ISP)
                    W1%CPROJ = 0
                    IF (NONLR_S%LREAL) THEN
                       CALL FFTWAV_W1(W1)
                       CALL RPRO1(NONLR_D(IDIR), WDES1, W1)
                    ELSE
                       CALL PROJ1(NONL_D(IDIR), WDES1, W1)
                    ENDIF
                    ! account for the fact that we projected onto -i < d p_i / dk_idir |
                    CPROJD(1:WDES1%NPRO,NB,IDIR) = W1%CPROJ(1:WDES1%NPRO) * (0._q, 1._q)
                 ENDDO
              ENDDO

              bands: DO NB = 1, WDES%NBANDS
                 ! only occupied bands contribute
                 IF (isEmpty(W%FERWE(NB,NK,ISP))) CYCLE bands

                 WEIGHT=WDES1%RSPIN*W%FERWE(NB,NK,ISP)*WDES%WTKPT(NK)

                 DO IDIR = 1, 3
                    WDK(IDIR)%CW => RPHI(:,NB,NK,ISP,IDIR)
                    WDK(IDIR)%CPROJ => RPHI_CPROJ(:,NB,NK,ISP,IDIR)
                    ! bring | d psi_nk / dk_idir > to real space
                    CALL FFTWAV_W1(WDK(IDIR))
                 ENDDO

                 ! CTMP1(i, j) =  < d u_nk / dk_idir | H_k + (e(0)_nk - 2 e_fermi) S_k | d u_nk / dk_jdir >
                 ! (second term on the RHS of Eq. 36 of PRB 107, 165157 (2023))
                 DO IDIR = 1, 3
                    DO JDIR = 1, 3
                       ! we will not have contributions for IDIR=JDIR (Levi-Civita, below)
                       IF (IDIR==JDIR) CYCLE

                       !  compute < d u_nk / dk_idir | H_k  | d u_nk / dk_jdir >
                       CTMP = 0
#ifdef test_eccp_
                       CALL ECCP_(HAMILTONIAN, GRID, LATT_CUR, WDES1, WDK(JDIR), WDK(IDIR), LMDIM, CDIJ, SV, ISP, CTMP)
#else
                       IF (ASSOCIATED(HAMILTONIAN%AVEC)) THEN
                          CALL ECCP_VEC(WDES1, WDK(JDIR), WDK(IDIR), LMDIM, CDIJ(1,1,1,ISP), WDES1%GRID, SV(1,ISP), HAMILTONIAN%AVEC, CTMP)
                       ELSE
                          CALL ECCP(WDES1, WDK(JDIR), WDK(IDIR), LMDIM, CDIJ(1,1,1,ISP), WDES1%GRID, SV(1,ISP), CTMP)
                       ENDIF
#endif
                       ! add < d u_nk / dk_idir | H_k | d u_nk / dk_jdir > * i/2 * w_k
                       CTMP1(IDIR,JDIR) = CTMP1(IDIR,JDIR) + CTMP * (0._q, 0.5_q) * WEIGHT ! WDES%WTKPT(NK)
                       ! compute < d u_nk / dk_idir | S_k | d u_nk / dk_jdir >
                       CTMP = W1_DOT( WDK(IDIR), WDK(JDIR), CQIJ)
                       ! and add < d u_nk / dk_idir | (e(0)_nk - 2 e_fermi) S_k | d u_nk / dk_jdir > * i/2 * w_k  for the orbital magnetization
                       CTMP1(IDIR,JDIR) = CTMP1(IDIR,JDIR) + CTMP * (W%CELEN(NB,NK,ISP)-2._q*EFERMI) * (0._q, 0.5_q) * WEIGHT ! WDES%WTKPT(NK)
                       ! and < d u_nk / dk_idir | S_k | d u_nk / dk_jdir > * i * w_k  for the Chern number
                       CHRN1(IDIR,JDIR) = CHRN1(IDIR,JDIR) + CTMP * (0._q, 1.0_q) * WEIGHT ! WDES%WTKPT(NK)
                    ENDDO
                 ENDDO

                 ! CTMP2(i, j) = \sum_ij < u_nk | d p_ik / dk_idir > (D_ij - e(0)_nk Q_ij) < d p_jk / dk_jdir | u_nk >
                 ! (part of the first term on the RHS of Eq. 36 of PRB 107, 165157 (2023): see Eq. 42 of the aforementioned paper)
                 DO IDIR = 1, 3
                    DO JDIR = 1, 3
                       ! we will not have contributions for IDIR=JDIR (Levi-Civita, below)
                       IF (IDIR==JDIR) CYCLE

                       E0NK = REAL(W%CELEN(NB,NK,ISP))
                       ! COVL(i) = \sum_j (D_ij - e(0)_nk Q_ij) < d p_jk / dk_jdir | u_nk >
                       COVL = 0 ; CALL OVERL1(WDES1, LMDIM, CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), E0NK, CPROJD(1,NB,JDIR), COVL(1))
                       ! \sum_ij < u_nk | d p_ik / dk_idir > COVL(i)
                       CTMP2(IDIR,JDIR) = CTMP2(IDIR,JDIR) + SUM( GCONJG(CPROJD(1:WDES1%NPRO,NB,IDIR))*COVL(1:WDES1%NPRO) ) * (0._q,-0.5_q) * WEIGHT ! WDES%WTKPT(NK)
                    ENDDO
                 ENDDO

                 ! CDD(i, idir, jdir) = < d p_ik / dk_idir | d u_nk / dk_jdir >
                 CDD = 0
                 DO IDIR = 1, 3
                    DO JDIR = 1, 3
                       ! we will not have contributions for IDIR=JDIR (Levi-Civita, below)
                       IF (IDIR==JDIR) CYCLE

                       WTMP%CPROJ => CDD(:,IDIR,JDIR)
                       IF (NONLR_S%LREAL) THEN
                          WTMP%CR => WDK(JDIR)%CR
                          CALL RPRO1(NONLR_D(IDIR), WDES1, WTMP)
                       ELSE
                          WTMP%CW => WDK(JDIR)%CW
                          CALL PROJ1(NONL_D(IDIR), WDES1, WTMP)
                       ENDIF
                       ! account for the fact that we projected onto -i < d p_i / dk_idir |
                       WTMP%CPROJ = WTMP%CPROJ * (0._q, 1._q)
                    ENDDO
                 ENDDO


                 ! third term on the RHS of Eq. 36 of PRB 107, 165157 (2023)
                 DO IDIR = 1, 3
                    DO JDIR = 1, 3
                       ! we will not have contributions for IDIR=JDIR (Levi-Civita, below)
                       IF (IDIR==JDIR) CYCLE

                       CTMP = &
                            ! \sum_ij ( < u_nk | d p_ik / dk_idir > Q_ij < p_jk | d u_nk / dk_jdir > +
                          ( SUM( GCONJG(CPROJD(1:WDES1%NPRO,NB,IDIR))*OVL2(1:WDES1%NPRO,NB,JDIR) ) +&
                            ! < u_nk | p_ik > Q_ij < d p_jk / dk_idir | d u_nk / dk_jdir > +
                            SUM( GCONJG(OVL1(1:WDES1%NPRO,NB))*CDD(1:WDES1%NPRO,IDIR,JDIR) ) +&
                            ! < d u_nk / dk_idir | d p_ik / dk_jdir > Q_ij < p_jk | u_nk > +
                            SUM( GCONJG(CDD(1:WDES1%NPRO,JDIR,IDIR))*OVL1(1:WDES1%NPRO,NB) ) +&
                            ! < d u_nk / dk_idir | p_ik > Q_ij < d p_jk / dk_jdir | u_nk > ) *
                            SUM( GCONJG(OVL2(1:WDES1%NPRO,NB,IDIR))*CPROJD(1:WDES1%NPRO,NB,JDIR) ) )

                       ! add CTMP * (e(0)_nk - e_fermi) * i/2 * w_k for the orbital magnetization
                       CTMP3(IDIR,JDIR) = CTMP3(IDIR,JDIR) + CTMP * (W%CELEN(NB,NK,ISP)-EFERMI) * (0._q, 0.5_q) * WEIGHT ! WDES%WTKPT(NK)
                       ! and CTMP * i/2 * w_k for the Chern number
                       CHRN3(IDIR,JDIR) = CHRN3(IDIR,JDIR) + CTMP * (0._q, 0.5_q) * WEIGHT ! WDES%WTKPT(NK)
                    ENDDO
                 ENDDO

              ENDDO bands

              IF (LHFCALC_LOCAL) THEN
                 IF (NONLR_S%LREAL) THEN
                    CALL PHASER(GRID, LATT_CUR, NONLR_S, NK, WDES)
                 ELSE
                    CALL PHASE(WDES, NONL_S, NK)
                 ENDIF

                 DO JDIR = 1, 3
                    strips: DO NPOS = 1, WDES%NBANDS, NSTRIP
                       NSTRIP_ACT = MIN(WDES%NBANDS-NPOS+1, NSTRIP)
                       ! | WHF1_n > = | d u_nk / dk_jdir > for the strip of bands n = NPOS, .., NPOS+NSTRIP_ACT-1
                       DO NB = 1, NSTRIP_ACT
                          WHF1(NB)%CW => RPHI(:,NB+NPOS-1,NK,ISP,JDIR)
                          WHF1(NB)%CPROJ => RPHI_CPROJ(:,NB+NPOS-1,NK,ISP,JDIR)
                          WHF1(NB)%WDES1 => WDES1
                       ENDDO
                       ! | WHFA_n > = X(W) | WHF1_n > for n = 1, NSTRIP_ACT
                       ! X(W) is the Fock operator constructed from the groundstate orbitals W
                       CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, NONLR_S, NONL_S, NK, ISP, &
                          1, NSTRIP_ACT, WHFA%CW(:,:), P, CQIJ(1,1,1,1), CDUMMY, WHF1)
                       DO IDIR = 1, 3
                          IF (IDIR==JDIR) CYCLE
                          DO NB = 1, NSTRIP_ACT
                             WEIGHT=WDES1%RSPIN*W%FERWE(NB+NPOS-1,NK,ISP)*WDES%WTKPT(NK)
                             WHF1(NB)%CW => RPHI(:,NB+NPOS-1,NK,ISP,IDIR)
                             ! CTMP = < d_u_nk / dk_idir | X(W) | d_u_nk / dk_jdir >
                             CTMP = W1_DOT(WHF1(NB), ELEMENT(WHFA, NB))
                             ! add CTMP * 1/2 * w_k
                             CTMP1(IDIR,JDIR) = CTMP1(IDIR,JDIR) + CTMP * (0._q, 0.5_q) * WEIGHT
                          ENDDO
                       ENDDO
                    ENDDO strips
                 ENDDO
              ENDIF

              ! redistribute OVL1 and CPROJD to "over-planewave" distribution
              IF (WDES1%DO_REDIS) THEN
                 CALL REDIS_PROJ(WDES1, WDES1%NBANDS, OVL1(1,1))
                 CALL REDIS_PROJ(WDES1, WDES1%NBANDS*3, CPROJD(1,1,1))
              ENDIF

              NOCC = COUNT(.NOT.isEmpty(W%FERTOT(:,NK,ISP)))

              ! DSDK(n, n', IDIR) = \sum_ij < u_nk | d p_ik / dk_idir > Q_ij < p_jk | u_n'k > + < u_nk | p_ik > Q_ij < d p_jk / dk_idir | u_n'k >
              DO IDIR = 1, 3
                 ! compute \sum_ij < u_nk | d p_ik / dk_idir > Q_ij < p_jk | u_n'k >
                 CALL GGEMM(trans, 'N', NOCC, NOCC, WDES1%NPRO_O_RED, &
                            one, CPROJD(1,1,IDIR), WDES1%NPROD_RED, OVL1(1,1), WDES1%NPROD_RED, &
                            zero, DSDK(1,1,IDIR), WDES1%NB_TOT )
                 ! add < u_nk | p_ik > Q_ij < d p_jk / dk_idir | u_n'k >
                 CALL GGEMM(trans, 'N', NOCC, NOCC, WDES1%NPRO_O_RED, &
                            one, OVL1(1,1), WDES1%NPROD_RED, CPROJD(1,1,IDIR), WDES1%NPROD_RED, &
                            one, DSDK(1,1,IDIR), WDES1%NB_TOT )
              ENDDO

              ! sum over all ranks
              CALLMPI( M_sum_g(WDES1%COMM_KIN, DSDK(1,1,1), SIZE(DSDK)) )

              ! \sum_nn' \sum_ij < u_nk | d p_ik / dk_idir> Qij < p_jk | u_n'k > + < u_nk | p_ik > Qij < d p_jk / dk_idir | u_n'k >
              !                  < u_n'k | d p_ik / dk_jdir> Qij < p_jk | u_nk > + < u_n'k | p_ik > Qij < d p_jk / dk_jdir | u_nk >
              !                  (e(0)_nk - e_fermi)
              ! (fourth term on the RHS of Eq. 36 of PRB 107, 165157 (2023))
              DO IDIR = 1, 3
                 DO JDIR = 1, 3
                    ! we will not have contributions for IDIR=JDIR (Levi-Civita, below)
                    IF (IDIR==JDIR) CYCLE

                    DO NB = 1, NOCC
                       WEIGHT=WDES1%RSPIN*W%FERTOT(NB,NK,ISP)*WDES%WTKPT(NK)
                       CTMP4(IDIR,JDIR) = CTMP4(IDIR,JDIR) + SUM( GCONJG(DSDK(1:NOCC,NB,IDIR))*DSDK(1:NOCC,NB,JDIR) )*(W%CELTOT(NB,NK,ISP)-EFERMI) * (0._q,-0.5_q) * WEIGHT ! WDES%WTKPT(NK)
                    ENDDO
                 ENDDO
              ENDDO

           ENDDO kpoint
        ENDDO spin

        ! CTMP1 has already been reduced over COMM_INB, so here we only need to reduce over COMM_INTER
        CALLMPI( M_sum_z(WDES%COMM_INTER, CTMP1, SIZE(CTMP1)) )
        ! and over the ranks in COMM_KINTER
        CALLMPI( M_sum_z(WDES%COMM_KINTER, CTMP1, SIZE(CTMP1)) )

        ! CTMP2 and CTMP3 need to be reduced over all ranks
        CALLMPI( M_sum_z(WDES%COMM, CTMP2, SIZE(CTMP2)) )
        CALLMPI( M_sum_z(WDES%COMM, CTMP3, SIZE(CTMP3)) )

        ! CTMP4 has already been reduced over COMM_KIN, and only needs a reduction over COMM_KINTER
        CALLMPI( M_sum_z(WDES%COMM_KINTER, CTMP4, SIZE(CTMP4)) )

        ! Collect quantities related to the Chern number
        CALLMPI( M_sum_z(WDES%COMM_INTER, CHRN1, SIZE(CHRN1)) )
        CALLMPI( M_sum_z(WDES%COMM_KINTER, CHRN1, SIZE(CHRN1)) )
        CALLMPI( M_sum_z(WDES%COMM, CHRN3, SIZE(CHRN3)) )

#ifdef debug
        MTMP1 = 0 ; MTMP2 = 0 ; MTMP3 = 0 ; MTMP4 = 0
        CHRN = 0
        DO IDIR = 1, 3
           DO KDIR = 1, 3
              DO JDIR = 1, 3
                 IF (EIJK(IDIR,JDIR,KDIR)==0) CYCLE
                 MTMP1(IDIR) = MTMP1(IDIR) + CTMP1(JDIR,KDIR) * EIJK(IDIR,JDIR,KDIR)
                 MTMP2(IDIR) = MTMP2(IDIR) + CTMP2(JDIR,KDIR) * EIJK(IDIR,JDIR,KDIR)
                 MTMP3(IDIR) = MTMP3(IDIR) + CTMP3(JDIR,KDIR) * EIJK(IDIR,JDIR,KDIR)
                 MTMP4(IDIR) = MTMP4(IDIR) + CTMP4(JDIR,KDIR) * EIJK(IDIR,JDIR,KDIR)

                 CHRN(IDIR) = CHRN(IDIR) + (CHRN1(JDIR,KDIR) + CHRN3(JDIR,KDIR)) * EIJK(IDIR,JDIR,KDIR)
              ENDDO
           ENDDO
        ENDDO

        WRITE(*,'(a,3(2x,2f20.10))') 'MTMP1',MTMP1 * 1.E6_q/(2._q*RYTOEV)/AUTOA/AUTOA
        WRITE(*,'(a,3(2x,2f20.10))') 'MTMP2',MTMP2 * 1.E6_q/(2._q*RYTOEV)/AUTOA/AUTOA
        WRITE(*,'(a,3(2x,2f20.10))') 'MTMP3',MTMP3 * 1.E6_q/(2._q*RYTOEV)/AUTOA/AUTOA
        WRITE(*,'(a,3(2x,2f20.10))') 'MTMP4',MTMP4 * 1.E6_q/(2._q*RYTOEV)/AUTOA/AUTOA
        WRITE(*,'(a,3(2x,2f20.10))') 'CHRN ',CHRN  * TPI*TPI
#endif
        ! Levi-Civita
        DO IDIR = 1, 3
           DO KDIR = 1, 3
              DO JDIR = 1, 3
                 IF (EIJK(IDIR,JDIR,KDIR)==0) CYCLE
                 MTOT(IDIR) = MTOT(IDIR) + ( CTMP1(JDIR,KDIR) + CTMP2(JDIR,KDIR) + CTMP3(JDIR,KDIR) + CTMP4(JDIR,KDIR) ) * EIJK(IDIR,JDIR,KDIR)
              ENDDO
           ENDDO
        ENDDO
        ! unit ?
        MTOT = MTOT /(2._q*RYTOEV)/AUTOA/AUTOA
        MTOT = MTOT*2._q ! For consistency with Ceresoli converse.

!>=============================================================================
!>
!>  M1C
!>
!>=============================================================================
#ifdef debug
        MPARA = 0 ; MDIA = 0
#endif
        ! \sum_nk \sum_ij < u_nk | p_ik > (D^1_ij - \tilde D^1_ij) < p_jk | u_nk >
        ! (remainder of the first term on the RHS of Eq. 36 of PRB 107, 165157 (2023): see Eqs. 33-35 of the aforementioned paper)
        ! loop over all atoms (global)
        DO NI = 1, T_INFO%NIONS
           ! local index of atom NI
           NIP = NI_LOCAL(NI, WDES%COMM_INB)
           ! if atom NI is not treated by this rank, then NIP=0
           ! ... in that case we cycle
           IF (NIP==0) CYCLE
           ! the type of this atom
           NT = T_INFO%ITYP(NI); PP => P(NT)

           IF (ASSOCIATED(PP%QPAW)) THEN
              IF (DO_LOCAL(NI)) THEN
                 MLQ = 0 ; CALL CALC_MU_LQ(PP, CRHODE(:,:,NIP,1), MLQ)
! mM: unsure about the sign here ...
!!                 M1C = M1C + MLQ
                 M1C = M1C - MLQ
#ifdef debug
                 MPARA(:,NI) = MLQ
#endif
                 IF (NI==MAGATOM) THEN
                    CALL CALC_MU_DIA(PP, CRHODE(:,:,NIP,1), MAGDIPOL, MD1C)
!!                    CALL CALC_MU_MAGATOM(PP, CRHODE(:,:,NIP,1), MP1C, MD1C, .FALSE.)
! mM: unsure about the sign here as well ...
!!                    M1C = M1C + MD1C
                    M1C = M1C - MD1C
#ifdef debug
                    MDIA = MD1C
#endif
                 ENDIF
              ENDIF
           ENDIF
        ENDDO

#ifdef debug
        CALLMPI( M_sum_d(WDES%COMM_INTER, MPARA, 3*T_INFO%NIONS) )
        CALLMPI( M_sum_d(WDES%COMM_INTER, MDIA, 3) )

        WRITE(*,'(a)') 'MPARA'
        DO NI = 1, T_INFO%NIONS
           WRITE(*,'(i4,2x,3f14.7)') NI, MPARA(:,NI) * 1.E6_q
        ENDDO
        WRITE(*,'(a)') 'MDIA'
        WRITE(*,'(i4,2x,3f14.7)') MAGATOM, MDIA * 1.E6_q
#endif
        CALLMPI( M_sum_d(WDES%COMM_INTER, M1C, 3) )

!>=============================================================================
!>
!>  WRAP-UP/OUTPUT
!>
!>=============================================================================
        MORB = -(MTOT + M1C)

#ifdef debug
         WRITE(*,'(a,2x,3f14.7)') 'MTOT',-MTOT * 1.E6_q
         WRITE(*,'(a,2x,3f14.7)') 'M1C ',-M1C  * 1.E6_q
         WRITE(*,'(a,2x,3f14.7)') 'MORB', MORB * 1.E6_q
#endif
        ! cleanup
        IF (LHFCALC_LOCAL) THEN
           CALL DELWAVA(WHFA)
           DEALLOCATE(WHF1)
        ENDIF

        DO IDIR = 1, 3
           IF (NONLR_S%LREAL) THEN
              CALL NONLR_DEALLOC_CRREXP(NONLR_D(IDIR))
           ELSE
              CALL NONL_DEALLOC_DER(NONL_D(IDIR))
           ENDIF
           CALL DELWAV_R(WDK(IDIR))
        ENDDO
        CALL DELWAV_NOCW(W1)

        RETURN
      END SUBROUTINE MORBCONV


!************************ SUBROUTINE CALC_MU_DIA ***********************
!
!> Compute the total moment of the diamagnetic one-center current
!> ~~~
!>     M = \int d3r r x \sum_ab < a | ( m_s x r) | b > \rho_ab
!> ~~~
!> First we compute the LM decomposed diamagnetic current
!> ~~~
!>     \vec{j} = \sum_LM \vec{J}_LM(r) Y_LM(\hat r)
!> ~~~
!>  where
!> ~~~
!>     J_LM (r)_idir = \sum_ij ( R_i(r) R_j(r) - \tilde R_i(r) \tilde R_j(r) ) / r^3 *
!>                     \sum_{m_i, m_j} \sum_{lm_k} CG(lm_i, lm_k, LM) r < lm_k | (m x \hat r)_idir | lm_j > * \rho_(lm_i,lm_j)
!>
!>             CG(lm_i, lm_k, LM) = \int \hat{r} Y_lm_i(\hat r) Y_lm_k(\hat r) Y_LM(\hat r)
!> ~~~
!>  where i and j run over the different (n, l) channels in the PAW dataset,
!>  m_i and m_j run over all m_i = -l_i, ..., l_i, m_j = -l_j, ..., l_j,
!>  and lm_k = (l_k, m_k) runs over all l_k = 0, ..., LMAXNABLA, m_k = -l_k, ..., l_k
!>
!>  With J_LM we have
!> ~~~
!>     M = \int d3r \vec{r} x \vec{j} = sqrt(4\pi) \sum_LM < 00 | \hat r | LM > x \int dr r J_LM(r)
!> ~~~
!>  in the above we used
!> ~~~
!>     \vec{r} = sqrt(4\pi) r \hat{r} Y_00(\hat r)
!>
!>             Y_00(\hat r) = 1/sqrt(4\pi)
!> ~~~
!
!***********************************************************************

      SUBROUTINE CALC_MU_DIA(PP, CRHODE, MS, MDIA)

        USE prec
        USE constant
        USE pseudo_struct_def
        USE asa

        !> PAW dataset
        TYPE (potcar), POINTER :: PP
        !> one-center occupancy matrix
        OVERLAP :: CRHODE(:,:)
        !> external magnetic dipole
        REAL(q) :: MS(3)
        !> diamagnetic moment
        REAL(q) :: MDIA(3)

        ! local variables
        INTEGER :: NMAX
        INTEGER :: I, J, K, IDIR
        INTEGER :: L_i, L_j, L_k, L
        INTEGER :: M_i, M_j, M_k
        INTEGER :: LM_i, LM_j, LM_k, LM, LMMAX
        INTEGER :: ISTART, IEND, LMIND, IC

        INTEGER :: LMAXNABLA, LMAX

        REAL(q) :: CG
        REAL(q), ALLOCATABLE :: TMP(:)
        REAL(q), ALLOCATABLE :: YLM_MxR_YLM(:,:,:)

        REAL(q), ALLOCATABLE :: YLM_NABLA_YLM(:,:,:),YLM_X_YLM(:,:,:)

        COMPLEX(q), ALLOCATABLE :: JLM(:,:,:), JLMR(:,:)

        COMPLEX(q) :: MU(3)
        INTEGER, EXTERNAL :: MAXL1

        ! maximum l-quantum number of the channels in this PAW dataset
        LMAX = MAXL1(PP)
        ! \nable Y_lm has components up to l+1
        LMAXNABLA = MAXL1(PP)+1

        ALLOCATE(YLM_NABLA_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAXNABLA+2)*(2*LMAXNABLA+2),0:3))
        ALLOCATE(YLM_X_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAX+2)*(2*LMAX+2),0:3))

        YLM_NABLA_YLM = 0
        YLM_X_YLM = 0

! After the call to SETYLM_NABLA_YLM
!
!    YLM_NABLA_YLM(lm, l'm', i) = |r| < lm | \nabla_i | l'm' >
!
!        YLM_X_YLM(lm, l'm', i) = < lm | r_i / |r| | l'm' >
!
! in principle we should go up to LMAX+LMAXNABLA+1 = 2*LMAX + 2, but currently
! that does not work for f-electrons, so we'll stick to 2*LMAX+1 for now
        CALL SETYLM_NABLA_YLM(2*LMAX+1,YLM_NABLA_YLM,YLM_X_YLM)

! YLM_MxR_YLM(lm, l'm', i) = < lm | (m x \hat r)_i | l'm' >
        ALLOCATE(YLM_MxR_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAX+2)*(2*LMAX+2),0:3))

        YLM_MxR_YLM(:,:,1)=MS(2)*YLM_X_YLM(:,:,3)-MS(3)*YLM_X_YLM(:,:,2)
        YLM_MxR_YLM(:,:,2)=MS(3)*YLM_X_YLM(:,:,1)-MS(1)*YLM_X_YLM(:,:,3)
        YLM_MxR_YLM(:,:,3)=MS(1)*YLM_X_YLM(:,:,2)-MS(2)*YLM_X_YLM(:,:,1)

! N.B.: CALC_MU_DIA actually only needs YLM_MxR_YLM, so at this point we
! could in principle deallocate YLM_NABLA_YLM and YLM_X_YLM ...

        NMAX=PP%R%NMAX
        ALLOCATE(TMP(NMAX), JLM(NMAX,(LMAXNABLA+LMAX+1)*(LMAXNABLA+LMAX+1),3))

!       J_LM (r)_idir = \sum_ij ( R_i(r) R_j(r) - \tilde R_i(r) \tilde R_j(r) ) / r^3 *
!                       \sum_{m_i} \sum_{m_j} \sum_{lm_k} CG(lm_i, lm_k, LM) r < lm_k | (m x \hat r)_idir | lm_j > * \rho_(lm_i,lm_j)

        JLM = 0 ; LMMAX = 0
        ! loop over all channels of the bra < i |
        LM_i = 1 ; DO I = 1, PP%LMAX
           ! loop over all channels of the ket | j >
           LM_j = 1 ; DO J = 1, PP%LMAX

              L_i = PP%LPS(I) ! l_i
              L_j = PP%LPS(J) ! l_j

              ! loop over m_j
              DO M_j = 1, 2*L_j+1
                 ! loop over l_k
                 LM_k = 1 ; DO L_k = 0, LMAXNABLA

                    ! LMIND is the starting index for all non-zero Clebsch-Gordan coefficients
                    ! \int Y_{l_i, m_i} Y_{l_k, m_k} Y_{L, M} d\hat r  for a particular l_i and l_k
                    CALL YLM3LOOKUP(L_i, L_k, LMIND)

                    ! loop over m_i
                    DO M_i = 1, 2*L_i+1
                       ! loop over m
                       DO M_k = 1, 2*L_k+1

                          LMIND = LMIND+1
                          ISTART = INDCG(LMIND)
                          IEND = INDCG(LMIND+1)

                          ! loop over all LM for which \int Y_{l_i, m_i} Y_{l_k, m_k} Y_{L, M} d\hat r /= 0
                          DO IC = ISTART, IEND-1

                             ! CG(LM_i, LM_k, LM) = \int Y_{l_i, m_i} Y_{l_k, m_k} Y_{L, M} d\hat r
                             CG = YLM3(IC)

                             LM = JS(IC) ! LM
                             L  = JL(IC) ! L

                             LMMAX = MAX(LM, LMMAX)

                             ! loop over cartesian directions
                             DO IDIR = 1, 3
                                ! if < LM_k | (m x \hat r)_idir | LM_j > /= 0
                                IF (ABS(YLM_MxR_YLM(LM_k+M_k-1,L_j*L_j+M_j,IDIR))>1E-10_q) THEN

                                   ! TMP(r) = CG( LM_i, LM_k, LM) ( R_i(r) R_j(r) - \tilde R_i(r) \tilde R_j(r) )/ r^2 < LM_k | (m x \hat r)_idir | LM_j >
                                   TMP = 0
                                   DO K = 1, NMAX
                                      TMP(K) = CG* YLM_MxR_YLM(LM_k+M_k-1,L_j*L_j+M_j,IDIR) / PP%R%R(K)**2 * &
                                               ( PP%WAE(K,I) * PP%WAE(K,J) - PP%WPS(K,I) * PP%WPS(K,J) )
                                   ENDDO

                                   ! WRITE (*,'(3I3,2f16.10)') J,I,L2,PP%QPAW(I,J,L2)

                                   ! J_LM (r)_idir = J_LM (r)_idir  - TMP(r) * \rho_{lm_i, lm_j}
                                   DO K = 1, NMAX
                                      JLM(K,LM,IDIR) = JLM(K,LM,IDIR) - TMP(K)*CRHODE(LM_i+M_i-1,LM_j+M_j-1)*MOMTOMOM
                                   ENDDO
                                ENDIF
                             ENDDO  ! IDIR

                          ENDDO  ! IC1

                       ENDDO ! M_k
                    ENDDO ! M_i

                    LM_k=LM_k+2*L_k+1
                 ENDDO ! L_k

              ENDDO ! M_j

              LM_j=LM_j+2*L_j+1
           ENDDO ! J

           LM_i=LM_i+2*L_i+1
        ENDDO ! I


        ! \int dr r J_LM(r)
        ALLOCATE(JLMR(LMMAX,3)); JLMR = 0
        DO IDIR = 1, 3
           DO LM = 1, LMMAX
              DO K = 1, NMAX
                 JLMR(LM,IDIR) = JLMR(LM,IDIR) + JLM(K,LM,IDIR) *PP%R%R(K) *PP%R%SI(K)
              ENDDO
           ENDDO
        ENDDO

! calculate the diamagnetic moment
!
!       \int d3r \vec{r} x \vec{j} = sqrt(4\pi) \sum_LM < 00 | \hat r | LM > x \int dr J_LM(r)
!
!    in the above we used
!
!       \vec{r} = sqrt(4\pi) r \hat{r} Y_00(\hat r)  since Y_00(\hat r) = 1/sqrt(4\pi)
!
!       \vec{j} = \sum_LM Y_LM(\hat r) J_LM(r)

        MU = 0
        DO LM = 1, LMMAX
           MU(1) = MU(1) + (YLM_X_YLM(1,LM,2)*JLMR(LM,3)-YLM_X_YLM(1,LM,3)*JLMR(LM,2))
           MU(2) = MU(2) + (YLM_X_YLM(1,LM,3)*JLMR(LM,1)-YLM_X_YLM(1,LM,1)*JLMR(LM,3))
           MU(3) = MU(3) + (YLM_X_YLM(1,LM,1)*JLMR(LM,2)-YLM_X_YLM(1,LM,2)*JLMR(LM,1))
        ENDDO
        ! divide by Y_00
        MU = MU*SQRT(4._q*PI)

        ! WRITE (*,'("one center magnetic moment*1E6",6E16.7)') MU*1E6

        MDIA = REAL(MU, KIND=q)

        DEALLOCATE(YLM_NABLA_YLM, YLM_X_YLM, YLM_MxR_YLM)
        DEALLOCATE(TMP,JLM,JLMR)

        RETURN
      END SUBROUTINE CALC_MU_DIA


!************************ SUBROUTINE PROJ_IN_EMPTY_SUBSPACE_BLAS *******
!
!***********************************************************************

      SUBROUTINE PROJ_IN_EMPTY_SUBSPACE_BLAS( &
         W0, GRID, NONL_S, NONLR_S, LATT_CUR, LMDIM, CQIJ, INFO, W1, LEFTS)

        USE prec
        USE base
        USE dfast
        USE lattice
        USE mgrid
        USE wave_high
        USE nonl_high

        IMPLICIT NONE

        TYPE (wavespin) W0
        TYPE (wavespin) W1
        TYPE (grid_3d) GRID
        TYPE (nonl_struct) NONL_S
        TYPE (nonlr_struct) NONLR_S
        TYPE (latt) LATT_CUR
        TYPE (info_struct) INFO
        INTEGER LMDIM
        OVERLAP CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)

        LOGICAL, OPTIONAL :: LEFTS

        ! local variables
        TYPE (wavedes1) WDES1
        TYPE (wavefun1), ALLOCATABLE :: WTMP(:)
        ! redistributed plane wave coefficients
        COMPLEX(q), POINTER :: CW0_RED(:,:),CW1_RED(:,:),CH_RED(:,:)
        INTEGER NCPU,NRPLWV_RED,NPL,NPRO
        LOGICAL DO_REDIS

        INTEGER ISP,NK,ISPINOR
        INTEGER NB_TOT,NBANDS,NSTRIP,NSTRIP_ACT,NB_START,N
        INTEGER NPOSB1,NPOSB2,NSTRPB,NPOSPL,NSTRPL

        REAL(q), ALLOCATABLE :: DUMMY(:)
        LOGICAL, ALLOCATABLE :: LDO(:)

        COMPLEX(q), ALLOCATABLE, TARGET :: CH(:,:)
        COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)

        GDEF,ALLOCATABLE,TARGET :: COVL(:,:)

        LOGICAL :: LEFT

        LEFT = .FALSE. ; IF (PRESENT(LEFTS)) LEFT = LEFTS
#ifdef MPI
        ! number of procs involved in band distribution
        NCPU=W0%WDES%COMM_INTER%NCPU
#else
        NCPU=1
#endif
        IF (NCPU/=1) THEN
           DO_REDIS=.TRUE.
           NRPLWV_RED=W0%WDES%NRPLWV/NCPU
        ELSE
           DO_REDIS=.FALSE.
           NRPLWV_RED=W0%WDES%NRPLWV
        ENDIF

        NB_TOT=W0%WDES%NB_TOT
        NBANDS=W0%WDES%NBANDS

!!        ! set NSTRIP between [1 and 32]
!!        NSTRIP=NSTRIP_STANDARD_GLOBAL
        ! the above was much too large
        NSTRIP=MIN(4,NBANDS)

        CALL SETWDES(W0%WDES,WDES1,0)
        ALLOCATE(WTMP(NSTRIP))
        DO N=1,NSTRIP
           CALL NEWWAV(WTMP(N),WDES1,INFO%LREAL)
        ENDDO

        ALLOCATE(DUMMY(NSTRIP),LDO(NSTRIP))

        IF (INFO%LREAL) THEN
           ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,NSTRIP))
        ENDIF

        ALLOCATE(CH(W0%WDES%NRPLWV,NBANDS))
        ! get pointers for redistribution to "over-plane-wave coefficients"
        IF (NCPU/=1) THEN
           CALL SET_WPOINTER( CH_RED, NRPLWV_RED, NB_TOT, CH(1,1))
        ELSE
           CH_RED => CH(:,:)
        ENDIF

        ALLOCATE(COVL(NB_TOT,NB_TOT))
!=======================================================================
        spin:  DO ISP=1,W0%WDES%ISPIN
        kpoints: DO NK=1,W0%WDES%NKPTS
#ifdef MPI
           IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!=======================================================================
           CALL SETWDES(W0%WDES,WDES1,NK)

           IF (INFO%LREAL) THEN
              CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,W0%WDES)
           ELSE
              NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK)
           ENDIF

           ! compute S|phi(0)> and store it in CH
           CH=(0._q,0._q)
           DO NB_START=1,NBANDS,NSTRIP
              NSTRIP_ACT=MIN(NSTRIP,NBANDS-NB_START+1)

              DO N=1,NSTRIP_ACT
                 CALL W1_COPY(ELEMENT(W0,WDES1,NB_START+N-1,ISP),WTMP(N)); IF (INFO%LREAL) CALL FFTWAV_W1(WTMP(N))
                 CH(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB_START+N-1)=WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
              ENDDO

              IF (INFO%LREAL) THEN
                 CWORK1=0; LDO=.FALSE.; LDO(1:NSTRIP_ACT)=.TRUE.; DUMMY=0
                 CALL RACCMU(NONLR_S,WDES1,WTMP,LMDIM,CQIJ,CQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSTRIP,LDO)
                 DO N=1,NSTRIP_ACT
                    DO ISPINOR=0,WDES1%NRSPINORS-1
                       CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,N),CH(1+ISPINOR*WDES1%NGVECTOR,NB_START+N-1),GRID,.TRUE.)
                    ENDDO
                 ENDDO
              ELSE
                 DO N=1,NSTRIP_ACT
                    CALL VNLACC_ADD(NONL_S,WTMP(N),CQIJ,CQIJ,1,0._q,CH(:,NB_START+N-1))
                 ENDDO
              ENDIF
           ENDDO

           ! get pointers for redistribution to "over-plane-wave coefficients"
           IF (NCPU/=1) THEN
              CALL SET_WPOINTER(CW0_RED, NRPLWV_RED, NB_TOT,  W0%CW(1,1,NK,ISP))
              CALL SET_WPOINTER(CW1_RED, NRPLWV_RED, NB_TOT, W1%CW(1,1,NK,ISP))
           ELSE
              CW0_RED => W0%CW(:,:,NK,ISP)
              CW1_RED => W1%CW(:,:,NK,ISP)
           ENDIF

           ! set number of wavefunctions after redistribution
           NPL = WDES1%NPL
           NPRO= WDES1%NPRO
           CALLMPI(SET_NPL_NPRO(WDES1, NPL, NPRO))

           ! redistribute to "over-plane-wave coefficients"
           IF (DO_REDIS) THEN
              CALL REDIS_PW(WDES1, NBANDS, W0%CW(1,1,NK,ISP))
              CALL REDIS_PW(WDES1, NBANDS, W1%CW(1,1,NK,ISP))
              CALL REDIS_PW(WDES1, NBANDS, CH(1,1))
           ENDIF

           IF (LEFT) THEN
              ! calculate COVL(i,j)=<phi(0)_i|phi(1)_j>
              COVL=zero
              IF (NPL/=0) THEN
                 NSTRPB=NB_TOT
                 DO NPOSB1=1,NB_TOT,NSTRPB
                    DO NPOSB2=1,NB_TOT,NSTRPB
                       CALL GGEMM('C','N',MIN(NSTRPB,NB_TOT-NPOSB2+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),m_ NPL, &
                          one,CW0_RED(1,NPOSB2),m_ NRPLWV_RED,CW1_RED(1,NPOSB1),m_ NRPLWV_RED, &
                          zero,COVL(NPOSB2,NPOSB1),NB_TOT)
                    ENDDO
                 ENDDO
              ENDIF
           ELSE
              ! calculate COVL(i,j)=<phi(0)_i|S|phi(1)_j>
              COVL=zero
              IF (NPL/=0) THEN
                 NSTRPB=NB_TOT
                 DO NPOSB1=1,NB_TOT,NSTRPB
                    DO NPOSB2=1,NB_TOT,NSTRPB
                       CALL GGEMM('C','N',MIN(NSTRPB,NB_TOT-NPOSB2+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),m_ NPL, &
                          one,CH_RED(1,NPOSB2),m_ NRPLWV_RED,CW1_RED(1,NPOSB1),m_ NRPLWV_RED, &
                          zero,COVL(NPOSB2,NPOSB1),NB_TOT)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF

           CALLMPI( M_sum_g(W0%WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))

           ! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
           DO N=1,W0%WDES%NB_TOT
              IF (W0%FERTOT(N,NK,ISP)<1E-4_q) COVL(N,:)=zero
           ENDDO

           IF (LEFT) THEN
              ! subtract S|phi(0)_j><phi(0)_j|phi(1)_i> from |phi(1)_i>
              IF (NPL/=0) THEN
                 NSTRPB=NB_TOT; NSTRPL=m_ NPL
                 DO NPOSB1=1,NB_TOT,NSTRPB
                    DO NPOSPL=1,m_ NPL,NSTRPL
                       CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),NB_TOT, &
                         -one,CH_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,NPOSB1),NB_TOT, &
                          one,CW1_RED(NPOSPL,NPOSB1),m_ NRPLWV_RED)
                    ENDDO
                 ENDDO
              ENDIF
           ELSE
              ! subtract |phi(0)_j><phi(0)_j|S|phi(1)_i> from |phi(1)_i>
              IF (NPL/=0) THEN
                 NSTRPB=NB_TOT; NSTRPL=m_ NPL
                 DO NPOSB1=1,NB_TOT,NSTRPB
                    DO NPOSPL=1,m_ NPL,NSTRPL
                       CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),NB_TOT, &
                         -one,CW0_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,NPOSB1),NB_TOT, &
                          one,CW1_RED(NPOSPL,NPOSB1),m_ NRPLWV_RED)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF

           ! redistribute back to original data distribution
           IF (DO_REDIS) THEN
              CALL REDIS_PW(WDES1, NBANDS, W0%CW(1,1,NK,ISP))
              CALL REDIS_PW(WDES1, NBANDS, W1%CW(1,1,NK,ISP))
           ENDIF

!=======================================================================
        ENDDO kpoints
        ENDDO spin
!=======================================================================

        ! deallocation
        IF (NCPU/=1) NULLIFY(CW0_RED,CW1_RED,CH_RED)

        DO N=1,NSTRIP
           CALL DELWAV(WTMP(N),INFO%LREAL)
        ENDDO

        DEALLOCATE(WTMP,DUMMY,LDO,CH,COVL)

        IF (INFO%LREAL) DEALLOCATE(CWORK1)

        RETURN
      END SUBROUTINE PROJ_IN_EMPTY_SUBSPACE_BLAS

      END MODULE mconverse
