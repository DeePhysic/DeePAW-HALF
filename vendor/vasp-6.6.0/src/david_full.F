#include "symbol.inc"
! #define debug

#if defined(MPI) && defined(scaLAPACK)

! Forcibly add at least one residual per orbital to the reduced basis
#define always_compute_first_residual

! Add the diagonal elements (in reciprocal space) of the PAW strength part the Hamiltonian
!!#define non_local_contributions_in_preconditioner

! Use the fast reciprocal-space projector scheme (when LREAL = F)
#define wa_proj_cprojk
#define wa_vnlacc_cprojk
#if defined(wa_proj_cprojk) || defined(wa_vnlacc_cprojk)
#define cprojk
#endif

! Asynchronous redistribution in separate buffers using NCCL
!!#define nccl_async_dblbuf
#if defined(nccl_async_dblbuf) && defined(USENCCL)
#define GPACC2
#else
#define GPACC2 !!
#endif

! Use divide-and-conquer
#define use_pdsyevd_pzheevd

! Extra verbose
!!#define verbose

      MODULE david_full

      USE prec
      USE wave_struct_def
      USE wave_mpi, ONLY : redis_pw_ctr

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EDDAV_FULL

      TYPE redis_handle
        TYPE(redis_pw_ctr), POINTER :: H
        !> local index of the first orbital that is redistributed in this handle
        INTEGER :: ISTART
        !> number of orbitals being redistributed in this handle
        INTEGER :: NSTRIP
      END TYPE redis_handle

      TYPE buffer
         !> Use asynchronous redistribution or not?
         LOGICAL :: LASYNC
         !> Handles for asynchronous redistribution
         TYPE(redis_handle), ALLOCATABLE :: HANDLE(:)
         !> Descriptor for use with the asynchronous redsitribution handles
         TYPE(wavedes) :: WDES
      END TYPE buffer

      CONTAINS

!************************ SUBROUTINE EDDAV_FULL ************************
!
!***********************************************************************

      SUBROUTINE EDDAV_FULL(HAMILTONIAN, SV, CDIJ, CQIJ, NONLR_S, NONL_S, INFO, LATT_CUR, WDES, W, &
            DEX, ETHR, RMS, DESUM, EXHF, NRED, IU0, IU6, LDELAY, LRESIDUAL, LHF, NKSTART, NKSTOP)

__GPU   USE moffload_struct_def
__GPU   USE moffload, ONLY : COPYIN_TYPED_VAR, DELETE_TYPED_VAR

        USE base

        USE wave_struct_def
        USE hamil_struct_def
        USE nonl_struct_def
        USE nonlr_struct_def
        USE poscar_struct_def

        USE nonl_high, ONLY : PHASE, PHASER
        USE wave_high, ONLY : SETWDES, ELEMENTS, NEWWAVA, DELWAVA, REDISTRIBUTE_PW, REDISTRIBUTE_PROJ, determineNumberOccupied
        USE scala, ONLY : LscaAWARE, INIT_scala, SCALA_NP, SCALA_NQ, GSD, BG_CHANGE_DIAGONALE
        USE fock_glb, ONLY : LHFCALC, AEXX
        USE fock_ace, ONLY : LFOCK_ACE
        USE pead, ONLY : LUSEPEAD

        USE tutor, ONLY : vtutor
        USE string, ONLY : str

        USE iso_c_binding, ONLY : c_loc, c_f_pointer

#if defined(verbose) && defined(PROFILING)
        USE profiling, ONLY : PRINT_ROUTINE_PROFILE
#endif

        IMPLICIT NONE

        !> contains the magnetic vector potential and the derivative
        !> of the xc-energy density w.r.t. the kinetic energy density
        TYPE(ham_handle) :: HAMILTONIAN
        !> local potential
        RGRID :: SV(:,:)
        !> PAW strength parameters and augmentation charges
        OVERLAP :: CDIJ(:,:,:,:), CQIJ(:,:,:,:)
        !> real space projection operators
        TYPE (nonlr_struct) :: NONLR_S
        !> reciprocal space projection operators
        TYPE (nonl_struct) :: NONL_S
        !> a whole bunch of control variables
        TYPE (info_struct) :: INFO
        !> the lattice vectors
        TYPE (latt) :: LATT_CUR

        !> descriptor of wave functions
        TYPE(wavedes) :: WDES
        !> wave functions
        TYPE(wavespin) :: W

        !> An estimate of the expected change in the total energy
        REAL(q) :: DEX

        !> orbital energy threshold from previous call
        REAL(q) :: ETHR
        !> sum of the norm of the final residuals
        REAL(q) :: RMS
        !> sum of the absolute changes in the eigenvalues
        REAL(q) :: DESUM
        !> Fock contribution to the double counting energy
        REAL(q) :: EXHF
        !> number of residuals that were calculated
        INTEGER :: NRED

        !> output channels
        INTEGER :: IU0, IU6

        !> delay phase? If yes, IPREC = 8
        LOGICAL :: LDELAY

        !> enforce decreasing residuals
        LOGICAL :: LRESIDUAL

        !> add Fock contributions (ACE)
        LOGICAL :: LHF

        !> treat k-points NKSTART .. NKSTOP
        INTEGER, OPTIONAL :: NKSTART, NKSTOP

        ! local variables
        TYPE(wavedes) :: WDES_DAV

        ! Orbitals and action of Hamiltonian and overlap operator on the orbitals
        TYPE(wavefuna) :: WA, HWA, SWA
        ! Common descriptor for WA, HWA, and SWA
        TYPE(wavedes1) :: WDES1
        ! Hamiltonian and overlap in reduced representation
        GDEF, ALLOCATABLE :: CHAM(:,:), COVL(:,:)
        ! Eigenvectors
        GDEF, ALLOCATABLE :: EVEC(:,:)
        ! Eigenvalues
        REAL(q), ALLOCATABLE :: EV(:), EV_PREV(:), EV_INI(:)

        ! For the preconditioner
        REAL(q), ALLOCATABLE :: DATAKE_RED(:)
#ifdef non_local_contributions_in_preconditioner
        REAL(q), ALLOCATABLE :: DIIG_RED(:), QIIG_RED(:)
#endif
        REAL(q) :: SLOCAL, DEVMAX
        INTEGER :: IPREC

        REAL(q) :: E_THRESHOLD_DYNAMIC
        REAL(q) :: E_THRESHOLD_STRICT, E_THRESHOLD_LOOSE

        LOGICAL :: LSUCCESS

        INTEGER :: NKSTART_, NKSTOP_, NK
        INTEGER :: ISP
        INTEGER :: NB_MAX, NB_TOT, NB_TOT_PREV
        INTEGER :: NKRYLOV
        INTEGER :: ITER, MAXITER, I
        INTEGER :: NRESTART
        INTEGER :: NDEPTH, MAXDEPTH

        LOGICAL :: LFORCE_FIRST_RESIDUAL

        INTEGER, ALLOCATABLE :: NOCC(:)
        INTEGER :: NSTRICT

        REAL(q), ALLOCATABLE :: FNORM(:), FNORM_INI(:), RWORK(:)
        INTEGER, ALLOCATABLE :: IDO(:), IDEPTH(:)
        LOGICAL, ALLOCATABLE :: LDO(:)
        INTEGER :: NDO

        REAL(q) :: EDCHF

        REAL(q) :: TMP1, TMP2

        LOGICAL :: LHF_LOCAL
        LOGICAL :: LscaAWARE_LOCAL

        PROFILING_START('eddav_full')

GPACC !$ACC WAIT IF(OFFLOAD_ON)
        SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)

        LscaAWARE_LOCAL = LscaAWARE
__GPU   LscaAWARE_LOCAL = LscaAWARE_LOCAL .AND. (.NOT.OFFLOAD_ON)

        ! does the Hamiltonian include Fock contributions?
        LHF_LOCAL = (LHF .AND. LHFCALC .AND. AEXX/=0)

        ! for the moment Fock contributions can only be added through the ACE
        IF (LHF_LOCAL .AND. .NOT.LFOCK_ACE()) &
           CALL vtutor%error("EDDAV_FULL: the non-blocked Davidson minimizer can only add&
              & Fock contributions to the Hamiltonian by means of the adaptively compressed&
              & exchange (ACE), please set LFOCKACE = .TRUE. in your INCAR file")

!-----------------------------------------------------------------------
! Control variables
!-----------------------------------------------------------------------
        E_THRESHOLD_DYNAMIC = MAX(ABS(INFO%EDIFF), DEX)
        E_THRESHOLD_DYNAMIC = E_THRESHOLD_DYNAMIC/INFO%NELECT*WDES%RSPIN/10

        E_THRESHOLD_STRICT = MIN(E_THRESHOLD_DYNAMIC, ETHR)
        E_THRESHOLD_LOOSE = MAX(5*E_THRESHOLD_STRICT, 1E-4_q)

        ! set the number of orbitals that will be converged tightly
        ! (i.e., in accordance with E_THRESHOLD_STRICT), the rest will
        ! only loosely converged (cf. E_THRESHOLD_LOOSE).
        CALL determineNumberOccupied(W%FERTOT, NOCC)

        IF (INFO%NSTRICT==-1) THEN
           ! by default we only converge (partially) occupied orbitals tightly
           NSTRICT = MAXVAL(NOCC)
        ELSE
           ! unless the user has specified something else ... but even then
           ! at least all occupied orbitals will be treated strictly
           ! (NOCC <= NSTRICT <= NB_TOT)
           NSTRICT = MAX(MAXVAL(NOCC), MIN(INFO%NSTRICT, WDES%NB_TOT))
        ENDIF

        ! the maximum number of iterations in this optimization step
        MAXITER = INFO%DAVITER

        ! NKRYLOV determines the maximum size of the search space:
        ! it will not exceed NKRYLOV*WDES%NB_TOT
        NKRYLOV = MAX(INFO%NKRYLOV,2)

        ! this selects the preconditioner
        IPREC = INFO%IALGO
#ifdef non_local_contributions_in_preconditioner
        ! Taking non-local contributions into accounr in the precconditioner
        ! is only implemented for the reciprocal space projector scheme
        IF (IPREC==1 .AND. INFO%LREAL) IPREC = 9
#endif
        ! the maximum depth of the search space for each orbital
        NDEPTH = INFO%NDAV

        ! only treat orbitals for which the change in the initial eigenvalue
        ! exceeds the thresholds (E_THRESHOLD_STRICT or E_THRESHOLD_LOOSE, respectively)
        LFORCE_FIRST_RESIDUAL = .FALSE.

        ! the "delay" phase is a special case
        IF (LDELAY) THEN
           E_THRESHOLD_STRICT = .25_q
           E_THRESHOLD_LOOSE = 1.0_q

           ! During the delay phase we do two iterations on each orbital at most
           NKRYLOV = MIN(NKRYLOV,4) ; NDEPTH = 10

!!           ! and we default to a preconditioner that does not use the current eigenenergies
!!           IPREC = 8

           ! make sure at least one residual is calculated for every orbital
           LFORCE_FIRST_RESIDUAL = .TRUE.
        ENDIF

#ifdef always_compute_first_residual
        LFORCE_FIRST_RESIDUAL = .TRUE.
#endif

        ! return the threshold that is used in this step to the caller
        ETHR = E_THRESHOLD_STRICT

#ifdef verbose
        !
        IF (IU0>=0) THEN
           WRITE(IU0,'(A,2X,L2   )') 'LDELAY                :',LDELAY
           WRITE(IU0,'(A,2X,L2   )') 'LRESIDUAL             :',LRESIDUAL
           WRITE(IU0,'(A,2X,L2   )') 'LFORCE_FIRST_RESIDUAL :',LFORCE_FIRST_RESIDUAL

           WRITE(IU0,'(A,2X,E10.3)') 'EDIFF                 :',INFO%EDIFF
           WRITE(IU0,'(A,2X,E10.3)') 'DEX                   :',DEX
           WRITE(IU0,'(A,2X,I4   )') 'NKRYLOV               :',NKRYLOV
           WRITE(IU0,'(A,2X,I4   )') 'NDEPTH                :',NDEPTH
           WRITE(IU0,'(A,2X,I4   )') 'NSTRICT               :',NSTRICT

           WRITE(IU0,'(A,2X,E10.3)') 'E_THRESHOLD_STRICT    :',E_THRESHOLD_STRICT
           WRITE(IU0,'(A,2X,E10.3)') 'E_THRESHOLD_LOOSE     :',E_THRESHOLD_LOOSE
           WRITE(IU0,*)
        ENDIF
#endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        DESUM = 0._q
        RMS = 0._q
        LSUCCESS = .TRUE.
        NRESTART = 0
        MAXDEPTH = 0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

        NKSTART_ = 1
        IF (PRESENT(NKSTART)) NKSTART_ = NKSTART

        NKSTOP_ = WDES%NKPTS
        IF (PRESENT(NKSTOP)) NKSTOP_ = NKSTOP

        WDES_DAV = WDES
        ! increase the number of orbitals by a factor NKRYLOV
        WDES_DAV%NBANDS = WDES%NBANDS*NKRYLOV
        WDES_DAV%NB_TOT = WDES%NB_TOT*NKRYLOV
        ! some care need to be taken when changing pointer members ...
        ALLOCATE(WDES_DAV%NB_TOTK, SOURCE=WDES%NB_TOTK)
        WDES_DAV%NB_TOTK = WDES%NB_TOTK*NKRYLOV

__GPU   CALL COPYIN_TYPED_VAR(WDES_DAV)

GPACC !$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1) IF(OFFLOAD_ON)
        CALL SETWDES(WDES_DAV, WDES1, 0)

GPACC !$ACC ENTER DATA CREATE(WA,HWA,SWA) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WA,HWA,SWA) IF(OFFLOAD_ON)
        CALL NEWWAVA(WA,  WDES1, WDES1%NBANDS)
        CALL NEWWAVA(HWA, WDES1, WDES1%NBANDS)
        CALL NEWWAVA(SWA, WDES1, WDES1%NBANDS)

        ! the maximum size of the reduced basis
        NB_MAX = WDES_DAV%NB_TOT

        ! a bit of a hack ...
        ! we need to pass down the fermi-weights to the WA_FOCK_ACE_APPLY_ACC
        ! to be able to compute Fock contribution to the double counting energy
        DEALLOCATE(WA%FERWE); ALLOCATE(WA%FERWE(NB_MAX)) ; WA%FERWE = 0
GPACC !$ACC ENTER DATA COPYIN(WA%FERWE) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:WA%FERWE) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:WA%FERWE) IF(OFFLOAD_ON)

        IF (.NOT.LscaAWARE_LOCAL) THEN
           ALLOCATE(CHAM(NB_MAX, NB_MAX))
           ALLOCATE(COVL(NB_MAX, NB_MAX))
           ALLOCATE(EVEC(NB_MAX, NB_MAX))
GPACC !$ACC ENTER DATA CREATE(CHAM,COVL,EVEC) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CHAM,COVL,EVEC) IF(OFFLOAD_ON)
        ELSE
           CALL INIT_scala(WDES_DAV%COMM_KIN, NB_MAX)
           ALLOCATE(CHAM(SCALA_NP(), SCALA_NQ()))
           ALLOCATE(COVL(SCALA_NP(), SCALA_NQ()))
           ALLOCATE(EVEC(SCALA_NP(), SCALA_NQ()))
        ENDIF

        ALLOCATE(EV(NB_MAX), EV_PREV(WDES%NB_TOT), EV_INI(WDES%NB_TOT))
GPACC !$ACC ENTER DATA CREATE(EV) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:EV) IF(OFFLOAD_ON)

        ! these arrays keeps track of the orbitals that still need to be optimized
        ALLOCATE(LDO(WDES%NB_TOT), IDO(WDES%NB_TOT), IDEPTH(WDES%NB_TOT), &
                 FNORM(WDES%NB_TOT), FNORM_INI(WDES%NB_TOT), RWORK(WDES%NB_TOT))

        ! average local potential, needed by the preconditioner
GPACC !$ACC ENTER DATA CREATE(SLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:SLOCAL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(SLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:SLOCAL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        SLOCAL = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP PRESENT(WDES,SLOCAL,SV) REDUCTION(+:SLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES,SLOCAL,SV) REDUCTION(+:SLOCAL) IF(OFFLOAD_ON)
        DO I = 1, WDES%GRID%RL%NP
           SLOCAL = SLOCAL + SV(I,1)
        ENDDO
        CALLMPI( M_sum_d(WDES%COMM_INB, SLOCAL, 1))

GPACC !$ACC KERNELS PRESENT(SLOCAL,WDES) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:SLOCAL,WDES,WDES%GRID) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        SLOCAL = SLOCAL / WDES%GRID%NPLWV
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        ALLOCATE(DATAKE_RED(WDES%NRPLWV_RED))
GPACC !$ACC ENTER DATA CREATE(DATAKE_RED) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:DATAKE_RED) IF(OFFLOAD_ON)
#ifdef non_local_contributions_in_preconditioner
        ALLOCATE(DIIG_RED(WDES%NRPLWV_RED), QIIG_RED(WDES%NRPLWV_RED))
#endif

        EXHF = 0
        NRED = 0

        spin: DO ISP = 1, WDES%ISPIN
        kpoints: DO NK = NKSTART_, NKSTOP_
#ifdef MPI
           IF (MOD(NK-1, WDES_DAV%COMM_KINTER%NCPU) /= WDES_DAV%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!-----------------------------------------------------------------------
! Startup phase
!-----------------------------------------------------------------------

           ! set the descriptor for wave functions at k-point NK
           CALL SETWDES(WDES_DAV, WDES1, NK)

           ! set the k-point dependent parts of the projectors
           IF (INFO%LREAL) THEN
              CALL PHASER(WDES_DAV%GRID, LATT_CUR, NONLR_S, NK, WDES_DAV)
           ELSE
              CALL PHASE(WDES_DAV, NONL_S, NK)
#ifdef cprojk
              CALL CPROJK(NONL_S, WDES1)
#endif
           ENDIF

           ! copy the wave functions at k-point NK and spin channel ISP to WA
           CALL WA_COPY(ELEMENTS(W, WDES1, ISP), WA)

           ! copy the Fermi-weights for all states at this k-point to WA%FERWE
           ! WA_FOCK_ACE_APPLY_ACC needs them to compute the Fock contribution
           ! to the double counting energy
           NB_TOT = WDES%NB_TOT
GPACC !$ACC KERNELS PRESENT(WA%FERWE, W%FERTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:WA%FERWE,W%FERTOT) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
           WA%FERWE(1:NB_TOT) = W%FERTOT(1:NB_TOT,NK,ISP)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

           ! make sure these wave functions are in over-plane-wave distribution
           IF (.NOT.W%OVER_BAND) THEN
              CALL REDISTRIBUTE_PW(WA)
           ENDIF
           ! signify that we are in over-plane-wave distribution now
           WA%OVER_BAND = .TRUE.

           ! compute H | WA > and S | WA >
           CALL HSPSI(HAMILTONIAN, SV, CDIJ, CQIJ, ISP, NONLR_S, NONL_S, LATT_CUR, EDCHF, LHF_LOCAL, WDES1, WA, HWA, SWA, 1, WDES%NB_TOT)
           ! keep track of the number of times we applied the hamiltonian
           NRED = NRED + WDES%NB_TOT
           ! the Fock contribution to the double counting energy is calculated
           ! from the original set of orbitals in the HSPSI call above
           EXHF = EXHF + EDCHF

           ! compute < WA | H | WA >
           CALL BRAKET(WDES1, WA, HWA, ISP, CHAM, 1, WDES%NB_TOT)
#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(CHAM) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(CHAM) IF(OFFLOAD_ON)
           CALL DUMP_HAM_BLOCK("Hamilton matrix eddav_full", CHAM, LscaAWARE_LOCAL, GSD%DESCA, WDES_DAV%COMM_KIN, 1, WDES%NB_TOT, 1, WDES%NB_TOT, IU0)
#endif

           ! compute < WA | S | WA >
           CALL BRAKET(WDES1, WA, SWA, ISP, COVL, 1, WDES%NB_TOT)
  
#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(COVL) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(COVL) IF(OFFLOAD_ON)
           CALL DUMP_HAM_BLOCK("overlap matrix eddav_full", COVL, LscaAWARE_LOCAL, GSD%DESCA, WDES_DAV%COMM_KIN, 1, WDES%NB_TOT, 1, WDES%NB_TOT, IU0)
#endif

           ! solve the generalized eigenvalue problem
           CALL HWESW(WDES1%COMM_KIN, CHAM, COVL, EV, EVEC, WDES%NB_TOT, IU0)
GPACC !$ACC UPDATE SELF(EV) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(EV) IF(OFFLOAD_ON)
#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(EVEC) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(EVEC) IF(OFFLOAD_ON)
           CALL DUMP_HAM_BLOCK("eigenvectors eddav_full", EVEC, LscaAWARE_LOCAL, GSD%DESCA, WDES_DAV%COMM_KIN, 1, WDES%NB_TOT, 1, WDES%NB_TOT, IU0)
           IF (IU0>=0) WRITE(*,'(F14.7)') EV(1:WDES%NB_TOT)
#endif
           ! store the eigenvalues
           EV_PREV(1:WDES%NB_TOT) = EV(1:WDES%NB_TOT)
           EV_INI(1:WDES%NB_TOT) = EV(1:WDES%NB_TOT)

           IF (LFORCE_FIRST_RESIDUAL) THEN
              LDO = .TRUE.
           ELSE
              ! convergence checks
              LDO(1:NSTRICT) = ABS(EV(1:NSTRICT)-W%CELTOT(1:NSTRICT,NK,ISP)) > E_THRESHOLD_STRICT

              IF (NSTRICT<WDES%NB_TOT) THEN
                 LDO(NSTRICT+1:WDES%NB_TOT) = ABS(EV(NSTRICT+1:WDES%NB_TOT)-W%CELTOT(NSTRICT+1:WDES%NB_TOT,NK,ISP)) > E_THRESHOLD_LOOSE
              ENDIF
           ENDIF

           ! the current size of the reduced basis
           NB_TOT = WDES%NB_TOT

           ! simply taken from EDDAV (used in the preconditioner):
           ! the difference between the largest and smallest eigenenergy (divided by 4)
           DEVMAX = ABS(MAXVAL(EV(1:WDES%NB_TOT))-MINVAL(EV(1:WDES%NB_TOT)))/4
           ! needed by the preconditioner as well
           CALL SET_DATAKE_RED(WDES1, DATAKE_RED)
#ifdef non_local_contributions_in_preconditioner
           IF (IPREC==1) CALL SET_NL_RED(WDES1, NONL_S, CDIJ, CQIJ, DIIG_RED, QIIG_RED)
#endif

           IDEPTH = 0
           FNORM = 0._q

!-----------------------------------------------------------------------
! Iterative phase
!-----------------------------------------------------------------------
           idav: DO ITER = 1, MAXITER

              NB_TOT_PREV = NB_TOT

              NDO = 0
              DO I = 1, SIZE(LDO)
                 IF (LDO(I)) THEN
                    NDO = NDO+1
                    IDO(NDO) = I
                    IDEPTH(I) = IDEPTH(I)+1
                 ENDIF
              ENDDO

              ! expand the basis by adding the residuals of the orbitals that
              ! are not sufficiently converged yet
              CALL EXPAND(WDES1, WA, HWA, SWA, EV, EVEC, NB_TOT, NDO, IDO)

              ! precondition and normalize the newly added basis functions
#ifdef non_local_contributions_in_preconditioner
              CALL PRECONDITION(IPREC, DATAKE_RED, DIIG_RED, QIIG_RED, SLOCAL, DEVMAX, EV, WDES1, WA, NB_TOT_PREV+1, NB_TOT, RWORK)
#else
              CALL PRECONDITION(IPREC, DATAKE_RED, SLOCAL, DEVMAX, EV, WDES1, WA, NB_TOT_PREV+1, NB_TOT, RWORK)
#endif
              ! keep track of the norm of the residuals
              DO I = 1, NDO
                 FNORM(IDO(I)) = RWORK(I)
              ENDDO
              ! store the initial residuals, may be used to construct a
              ! break-off criterion
              IF (ITER==1) FNORM_INI = FNORM

              ! compute H | WA > and S | WA >
              CALL HSPSI(HAMILTONIAN, SV, CDIJ, CQIJ, ISP, NONLR_S, NONL_S, LATT_CUR, EDCHF, LHF_LOCAL, WDES1, WA, HWA, SWA, NB_TOT_PREV+1, NB_TOT)
              ! keep track of the number of times we have applied the hamiltonian
              NRED = NRED + NB_TOT-NB_TOT_PREV

              ! compute < WA | H | WA >
              CALL BRAKET(WDES1, WA, HWA, ISP, CHAM, NB_TOT_PREV+1, NB_TOT)
#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(CHAM) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(CHAM) IF(OFFLOAD_ON)
          CALL DUMP_HAM_BLOCK("Hamilton matrix eddav_full #it: "//str(ITER), CHAM, LscaAWARE_LOCAL, GSD%DESCA, WDES_DAV%COMM_KIN, 1, NB_TOT, 1, NB_TOT, IU0)
#endif

              ! compute < WA | S | WA >
              CALL BRAKET(WDES1, WA, SWA, ISP, COVL, NB_TOT_PREV+1, NB_TOT)

#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(COVL) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(COVL) IF(OFFLOAD_ON)
              CALL DUMP_HAM_BLOCK("overlap matrix eddav_full #it: "//str(ITER), COVL, LscaAWARE_LOCAL, GSD%DESCA, WDES_DAV%COMM_KIN, 1, NB_TOT, 1, NB_TOT, IU0)
#endif

              ! elements on the diagonal of H and S must be real
              IF (LscaAWARE_LOCAL) THEN
                 CALL BG_CHANGE_DIAGONALE(WDES1%NB_TOT, CHAM(1,1), IU0, 1, NB_TOT)
                 CALL BG_CHANGE_DIAGONALE(WDES1%NB_TOT, COVL(1,1), IU0, 1, NB_TOT)
              ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(CHAM,COVL) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CHAM,COVL) IF(OFFLOAD_ON)
                 DO I = 1, NB_TOT
                    CHAM(I,I) = CMPLX(REAL(CHAM(I,I), KIND=q), 0._q, KIND=q)
                    COVL(I,I) = CMPLX(REAL(COVL(I,I), KIND=q), 0._q, KIND=q)
                 ENDDO
              ENDIF

              ! solve the generalized eigenvalue problem
              CALL HWESW(WDES1%COMM_KIN, CHAM, COVL, EV, EVEC, NB_TOT, IU0)
GPACC !$ACC UPDATE SELF(EV) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(EV) IF(OFFLOAD_ON)
#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(EVEC) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(EVEC) IF(OFFLOAD_ON)
              CALL DUMP_HAM_BLOCK("eigenvectors eddav_full #it: "//str(ITER), EVEC, LscaAWARE_LOCAL, GSD%DESCA, WDES_DAV%COMM_KIN, 1, NB_TOT, 1, NB_TOT, IU0)
              IF (IU0>=0) WRITE(*,'(F14.7)') EV(1:WDES%NB_TOT)
#endif

              ! convergence checks
              LDO(1:NSTRICT) = ABS(EV(1:NSTRICT)-EV_PREV(1:NSTRICT)) > E_THRESHOLD_STRICT

              IF (NSTRICT<WDES%NB_TOT) THEN
                 LDO(NSTRICT+1:WDES%NB_TOT) = ABS(EV(NSTRICT+1:WDES%NB_TOT)-EV_PREV(NSTRICT+1:WDES%NB_TOT)) > E_THRESHOLD_LOOSE
              ENDIF

              ! force the residuals of the first NSTRICT orbitals to decrease
              ! w.r.t. the previous iteration (stored in W%AUXTOT)
              IF (LRESIDUAL) THEN
                 LDO(1:NSTRICT) = LDO(1:NSTRICT) .OR. &
                    ((FNORM(1:NSTRICT) > W%AUXTOT(1:NSTRICT,NK,ISP)*0.75_q) .AND. (FNORM(1:NSTRICT) > INFO%FBREAK) .AND. (IDEPTH(1:NSTRICT) < 4))
              ENDIF

              ! during the delay phase the norm of the residuals of the occupied orbitals
              ! is forced to decrease
              IF (LDELAY) THEN
                 LDO(1:NSTRICT) = LDO(1:NSTRICT) .OR. ((FNORM(1:NSTRICT) > FNORM_INI(1:NSTRICT)*0.75_q))
              ENDIF

              ! if the maximum iteration depth is reached for an orbital
              ! it will not be optimized further
              DO I = 1, WDES%NB_TOT
                 IF (LDO(I) .AND. IDEPTH(I)==NDEPTH) THEN
                    LDO(I) = .FALSE. ; LSUCCESS = .FALSE.
                 ENDIF
              ENDDO

              ! how many orbitals should be worked on in the nexr iteration?
              NDO = COUNT(LDO)
#ifdef debug
              IF (IU0>=0) THEN
                 DO I = 1, WDES%NB_TOT
                    WRITE(*,'(I6,3F14.7,E12.4,I4)') I, EV_PREV(I), EV(I), EV(I)-EV_PREV(I), FNORM(I), IDEPTH(I)
                 ENDDO
                 WRITE(*,'(A,I4,4(2X,A,I6))') 'ISP:',ISP,'NK:',NK,'ITER:',ITER,'NB:',NB_TOT,'NDO:',NDO
              ENDIF
#endif
              EV_PREV(1:WDES%NB_TOT) = EV(1:WDES%NB_TOT)

              ! simply taken from EDDAV (used in the preconditioner):
              ! the difference between the largest and smallest eigenenergy (divided by 4)
              DEVMAX = ABS(MAXVAL(EV(1:WDES%NB_TOT))-MINVAL(EV(1:WDES%NB_TOT)))/4

              IF ( NDO==0 .OR. ITER==MAXITER ) THEN
                 ! up till now the projections have always been in the default distribution
                 ! ... redistribute them now
                 CALL REDISTRIBUTE_PROJ(WA)

                 ! extract WDES%NB_TOT orbitals (plane wave coefficients and projections)
                 CALL EXTRACT(WDES1, WA, EVEC, NB_TOT, WDES%NB_TOT, LPROJ=.TRUE.)

                 ! redistribute the projections back into the default distribution
                 CALL REDISTRIBUTE_PROJ(WA)
                 ! redistribute the plane wave coefficients to the original distribution
                 IF (.NOT. W%OVER_BAND) CALL REDISTRIBUTE_PW(WA)

                 ! copy them back to W
                 CALL WZCOPY(WDES%NBANDS*WDES1%NRPLWV, WA%CW(1,1), 1, W%CW(1,1,NK,ISP), 1)
                 CALL WDCOPY(ndata*WDES%NBANDS*WDES1%NPROD, WA%CPROJ(1,1), 1, W%CPROJ(1,1,NK,ISP), 1)

                 ! Update the sum of changes in the eigenvalues and the
                 ! sum of the norm of the final residuals
                 TMP1 = 0._q ; TMP2 = 0._q
                 DO I = 1, WDES%NB_TOT
                    ! add absolute changes in the eigenvalues to DESUM
                    TMP1 = TMP1 + W%FERTOT(I,NK,ISP)*ABS(EV(I)-EV_INI(I))
!!                    TMP1 = TMP1 + W%FERTOT(I,NK,ISP)*ABS(EV(I)-W%CELTOT(I,NK,ISP))

                    ! and sum the norm of the final residuals
                    TMP2 = TMP2 + W%FERTOT(I,NK,ISP)*FNORM(I)
                 ENDDO
                 DESUM = DESUM + TMP1*WDES%RSPIN*WDES%WTKPT(NK)
                 RMS = RMS + TMP2*WDES%NRSPINORS*WDES%RSPIN*WDES%WTKPT(NK)/WDES%NB_TOT

                 ! and store the eigenenergies
                 W%CELTOT(1:WDES%NB_TOT,NK,ISP) = EV(1:WDES%NB_TOT)

                 ! keep track of the final residuals as well
                 W%AUXTOT(1:WDES%NB_TOT,NK,ISP) = FNORM(1:WDES%NB_TOT)

                 ! What was the maximum iteration depth?
                 MAXDEPTH = MAX(MAXDEPTH, MAXVAL(IDEPTH(1:WDES%NB_TOT)))

                 ! Did we converge all orbitals?
                 LSUCCESS = LSUCCESS .AND. (NDO==0)
#ifdef verbose
                 IF (IU0>=0) THEN
                    DO I = 1, WDES%NB_TOT
                       WRITE(IU0,'(I6,4F14.7,E12.4,I4)') I, EV(I), EV_INI(I), EV(I)-EV_INI(I), W%FERTOT(I,NK,ISP), FNORM(I), IDEPTH(I)
                    ENDDO
                    WRITE(IU0,'(A,I4,3(2X,A,I6),2(2X,A,E10.3))') 'ISP:',ISP,'NK:',NK,'ITER:',ITER,'NB:',NB_TOT,'DESUM:',DESUM,'RMS:',RMS
                 ENDIF
#endif
                ! stop iterating
                 EXIT idav
              ELSEIF ( NB_TOT+NDO > NB_MAX ) THEN
                 ! extract reduced search space and continue

                 ! up till now the projections have always been in the default distribution
                 ! ... redistribute them now
                 CALL REDISTRIBUTE_PROJ(WA)

                 CALL EXTRACT(WDES1, WA,  EVEC, NB_TOT, WDES%NB_TOT, LPROJ=.TRUE.)

                 ! redistribute the projections back into the default distribution
                 CALL REDISTRIBUTE_PROJ(WA)

                 CALL EXTRACT(WDES1, HWA, EVEC, NB_TOT, WDES%NB_TOT)
                 CALL EXTRACT(WDES1, SWA, EVEC, NB_TOT, WDES%NB_TOT)

                 ! now reduce basis set size
                 NB_TOT = WDES%NB_TOT

                 ! we do not need to rotate these orbitals anymore
                 ! (has been done on extraction), so EVEC is now the
                 ! identity matrix
                 IF (LscaAWARE_LOCAL) THEN
                    CALL SET_TO_IDENTITY(NB_TOT, EVEC, GSD%DESCA)
                 ELSE
GPACC !$ACC KERNELS PRESENT(EVEC) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:EVEC) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
                    EVEC(1:NB_TOT,1:NB_TOT) = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP PRESENT(EVEC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:EVEC) IF(OFFLOAD_ON)
                    DO I = 1, NB_TOT
                       EVEC(I,I) = 1._q
                    ENDDO
                 ENDIF

                 ! the reduced Hamiltonian is now diagonal as well
                 IF (LscaAWARE_LOCAL) THEN
                    CALL SET_DIAGONAL(EV(1:NB_TOT), CHAM, GSD%DESCA)
                 ELSE
GPACC !$ACC KERNELS PRESENT(CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHAM) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
                    CHAM(1:NB_TOT,1:NB_TOT) = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP PRESENT(CHAM,EV) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CHAM,EV) IF(OFFLOAD_ON)
                    DO I = 1, NB_TOT
                       CHAM(I,I) = EV(I)
                    ENDDO
                 ENDIF

                 ! and similarly for the overlap matrix
                 IF (LscaAWARE_LOCAL) THEN
                    CALL SET_TO_IDENTITY(NB_TOT, COVL, GSD%DESCA)
                 ELSE
GPACC !$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:COVL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
                    COVL(1:NB_TOT,1:NB_TOT) = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP PRESENT(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:COVL) IF(OFFLOAD_ON)
                    DO I = 1, NB_TOT
                       COVL(I,I) = 1._q
                    ENDDO
                 ENDIF

                 ! keep track of the number if times we restart
                 NRESTART = NRESTART+1
              ELSE
                 ! continue
              ENDIF

           ENDDO idav
        ENDDO kpoints
        ENDDO spin

        IF (LHF_LOCAL) THEN
           CALLMPI( M_sum_d(W%WDES%COMM_KIN,    EXHF, 1))
           CALLMPI( M_sum_d(W%WDES%COMM_KINTER, EXHF, 1))
        ENDIF

#ifdef MPI
        IF (WDES%COMM_KINTER%NCPU>1) THEN
           CALL KPAR_SYNC_CELTOT(WDES,W)
           IF (LHF_LOCAL .OR. LUSEPEAD()) CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
        ENDIF
#endif

!-----------------------------------------------------------------------
! Deallocation
!-----------------------------------------------------------------------
GPACC !$ACC EXIT DATA DELETE(CHAM,COVL,EVEC,EV,DATAKE_RED,SLOCAL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CHAM,COVL,EVEC,EV,DATAKE_RED,SLOCAL) IF(OFFLOAD_ON)
        DEALLOCATE(CHAM, COVL, EVEC, EV, EV_PREV, EV_INI, LDO, IDO, IDEPTH, &
                   FNORM, FNORM_INI, RWORK, DATAKE_RED, NOCC)

#ifdef non_local_contributions_in_preconditioner
        DEALLOCATE(DIIG_RED, QIIG_RED)
#endif

#ifdef cprojk
        IF (.NOT.INFO%LREAL .AND. ASSOCIATED(NONL_S%CPROJK)) THEN
GPACC !$ACC EXIT DATA DELETE(NONL_S%CPROJK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:NONL_S%CPROJK) IF(OFFLOAD_ON)
#ifdef IFXWAR2
   GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:NONL_S%CPROJK) IF(OFFLOAD_ON)
#endif
           DEALLOCATE(NONL_S%CPROJK); NULLIFY(NONL_S%CPROJK)
        ENDIF
#endif

GPACC !$ACC EXIT DATA DELETE(WA%FERWE) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%FERWE) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA%FERWE) IF(OFFLOAD_ON)
        CALL DELWAVA(WA)
        CALL DELWAVA(HWA)
        CALL DELWAVA(SWA)
GPACC !$ACC EXIT DATA DELETE(WA,HWA,SWA) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WA,HWA,SWA) IF(OFFLOAD_ON)

__GPU   CALL DELETE_TYPED_VAR(WDES1)
__GPU   CALL DELETE_TYPED_VAR(WDES_DAV)

        ! Needs to be deallocated after the DELETE_TYPED_VAR(WDES_DAV) call
        DEALLOCATE(WDES_DAV%NB_TOTK)

#ifdef verbose
        IF (IU0>=0) WRITE(*,'(A,E10.3,2X,A,L2,2X,A,I6,2(2X,A,I3),2(2X,A,E10.3))') &
          'E_THRESHOLD:',E_THRESHOLD_STRICT,'LSUCCESS:',LSUCCESS,'NRED:',NRED,'MAXDEPTH:',MAXDEPTH,'NRESTART:',NRESTART,'DESUM:',DESUM,'RMS:',RMS
#endif

GPACC !$ACC WAIT IF(OFFLOAD_ON)
        SET_ACC_ASYNC_Q(INT(ACC_ASYNC_SYNC))

        PROFILING_STOP('eddav_full')

#if defined(verbose) && defined(PROFILING)
        IF (IU0>=0) CALL PRINT_ROUTINE_PROFILE(IU0, 'eddav_full', PFLAT=.FALSE.)
#endif
        RETURN
      END SUBROUTINE EDDAV_FULL


!************************ SUBROUTINE WA_COPY ***************************
!
!> Copy one wavefuna structure to another
!> ~~~
!>  W2 = W1  (W1 -> W2)
!> ~~~
!> The argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
!> the redistributed wavefunctions pointers in the destination
!> are also properly set
!>
!> This version of WA_COPY is almost a literal copy of the one in
!> wave_high. The difference is that here W2 may be larger than W1
!> and W2%WDES1 is not set to point to W1%WDES1
!
!***********************************************************************

      SUBROUTINE WA_COPY(W1, W2)

__GPU   USE moffload_struct_def

        USE wave_struct_def
        USE tutor, ONLY : vtutor
        USE string, ONLY : str

        IMPLICIT NONE

        TYPE (wavefuna), INTENT(IN) :: W1
        TYPE (wavefuna) :: W2

        PROFILING_START('wa_copy')

        IF (SIZE(W1%CW) > SIZE(W2%CW)) THEN
           CALL vtutor%bug("WA_COPY: destination CW too small " // str(SIZE(W1%CW)) // &
              " " // str(SIZE(W2%CW)), __FILE__, __LINE__)
        ENDIF

        IF (SIZE(W1%CPROJ) > SIZE(W2%CPROJ)) THEN
           CALL vtutor%bug("WA_COPY: destination CPROJ too small " // str(SIZE(W1%CPROJ)) &
              // " " // str(SIZE(W2%CPROJ)), __FILE__, __LINE__)
        ENDIF

        CALL WZCOPY(SIZE(W1%CW), W1%CW(1,1), 1, W2%CW(1,1), 1)

        IF (W1%WDES1%LGAMMA) THEN
           CALL WDCOPY(SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
        ELSE
           CALL WZCOPY(SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
        ENDIF

        ! one dimensional indexing assumed
        W2%FIRST_DIM=0

        ! remember spin index
        W2%ISP =W1%ISP

!!        ! set redistributed wavefunction indices
!!        IF (W2%WDES1%DO_REDIS) THEN
!!           CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CW(1,1))
!!           CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
!!        ELSE
!!           W2%CW_RED=>W2%CPTWFP
!!           W2%CPROJ_RED=>W2%CPROJ
!!        ENDIF

GPACC !$ACC UPDATE DEVICE(W2%FIRST_DIM,W2%ISP) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(W2%FIRST_DIM,W2%ISP) IF(OFFLOAD_ON)

        PROFILING_STOP('wa_copy')

        RETURN
      END SUBROUTINE WA_COPY

#ifdef cprojk
!!!************************ SUBROUTINE CPROJK ****************************
!!!
!!!***********************************************************************
!!
!!      SUBROUTINE CPROJK(NONL_S, WDES1)
!!
!!        USE nonl_struct_def
!!        USE wave_struct_def
!!
!!        USE wave, ONLY : NI_GLOBAL
!!
!!        IMPLICIT NONE
!!
!!        TYPE(nonl_struct) :: NONL_S
!!        TYPE(wavedes1) :: WDES1
!!
!!        ! local variables
!!        INTEGER :: NK, ISPINOR, ISPIRAL
!!        INTEGER :: NI, NT, NI_GLB, NT_GLB, LMMAXC, LMBASE, LM
!!        INTEGER :: NRPLWV_RED, NGVECTOR, I, IBASE, IGLB, ILOC, IRANK
!!
!!        PROFILING_START('cprojk')
!!
!!        IF (.NOT.ASSOCIATED(NONL_S%CPROJK)) &
!!           ALLOCATE(NONL_S%CPROJK(WDES1%NRPLWV_RED, WDES1%NPRO_TOT))
!!
!!        NK = WDES1%NK
!!
!!        ! maximum number of plane-wave components the ranks in COMM_INTER own
!!        ! in over-plane-wave distribution ... the last rank may own less
!!        NRPLWV_RED = WDES1%NRPLWV_RED
!!        ! size of the plane-wave basis
!!        NGVECTOR = WDES1%NGVECTOR
!!
!!        NONL_S%CPROJK = 0
!!
!!        spinor: DO ISPINOR = 0, WDES1%NRSPINORS-1
!!           ISPIRAL = 1 ; IF (NONL_S%LSPIRAL) ISPIRAL = ISPINOR+1
!!
!!           IBASE = ISPINOR * NGVECTOR
!!
!!           ion: DO NI = 1, WDES1%NIONS
!!              NT = WDES1%ITYP(NI)
!!              LMMAXC = NONL_S%LMMAX(NT)
!!              IF (LMMAXC == 0) CYCLE ion
!!
!!              LMBASE = WDES1%NPRO_POS(NI) + ISPINOR*WDES1%NPRO_TOT/2
!!
!!              NT_GLB = WDES1%NT_GLOBAL(NT)
!!              NI_GLB = NI_GLOBAL(NI, WDES1%COMM_INB)
!!
!!              DO LM = 1, LMMAXC
!!
!!                 DO I = 1, NGVECTOR
!!                    IGLB = I + IBASE
!!                    IRANK = (IGLB-1)/NRPLWV_RED + 1
!!                    ILOC = MOD(IGLB-1, NRPLWV_RED) + 1
!!                    IF (IRANK == WDES1%COMM_INTER%NODE_ME) THEN
!!                       NONL_S%CPROJK(ILOC, LM+LMBASE) = &
!!                          CONJG( NONL_S%QPROJ(I, LM, NT_GLB, NK, ISPIRAL) * NONL_S%CREXP(I, NI_GLB) * NONL_S%CQFAK(LM,NT_GLB) )
!!                    ENDIF
!!                 ENDDO
!!              ENDDO
!!
!!           ENDDO ion
!!
!!        ENDDO spinor
!!
!!        CALLMPI( M_sum_z(WDES1%COMM_INB, NONL_S%CPROJK(1,1), SIZE(NONL_S%CPROJK)))
!!
!!        PROFILING_STOP('cprojk')
!!
!!        RETURN
!!      END SUBROUTINE CPROJK


!************************ SUBROUTINE CPROJK ****************************
!
!***********************************************************************

      SUBROUTINE CPROJK(NONL_S, WDES1)

__GPU   USE moffload_struct_def

        USE nonl_struct_def
        USE wave_struct_def

        IMPLICIT NONE

        TYPE(nonl_struct) :: NONL_S
        TYPE(wavedes1) :: WDES1

        ! local variables
        INTEGER :: NK, ISPINOR, ISPIRAL
        INTEGER :: NI, LMBASE, LM
        INTEGER :: NRPLWV_RED, NGVECTOR, I, IBASE, IGLB, ILOC, IRANK, NODE_ME

        PROFILING_START('cprojk')

        IF (.NOT.ASSOCIATED(NONL_S%CPROJK)) THEN
           ALLOCATE(NONL_S%CPROJK(WDES1%NRPLWV_RED, WDES1%NPRO_TOT))
GPACC !$ACC ENTER DATA CREATE(NONL_S%CPROJK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:NONL_S%CPROJK) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: NONL_S%CPROJK) IF(OFFLOAD_ON)
        ENDIF

        NK = WDES1%NK

        ! maximum number of plane-wave components the ranks in COMM_INTER own
        ! in over-plane-wave distribution ... the last rank may own less
        NRPLWV_RED = WDES1%NRPLWV_RED
        ! size of the plane-wave basis
        NGVECTOR = WDES1%NGVECTOR

        NODE_ME = WDES1%COMM_INTER%NODE_ME

GPACC !$ACC KERNELS PRESENT(NONL_S%CPROJK) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:NONL_S%CPROJK) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        NONL_S%CPROJK = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(NONL_S) PRIVATE(ISPIRAL,IBASE,LMBASE,LM,I,IGLB,IRANK,ILOC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:NONL_S) &
GPOMP !$OMP PRIVATE(ISPIRAL,IBASE,LMBASE,LM,I,IGLB,IRANK,ILOC) IF(OFFLOAD_ON)
        spinor: DO ISPINOR = 0, WDES1%NRSPINORS-1
__CPU      ISPIRAL = 1 ; IF (NONL_S%LSPIRAL) ISPIRAL = ISPINOR+1
__CPU      IBASE = ISPINOR * NGVECTOR

           ion: DO NI = 1, NONL_S%NIONS
__GPU         ISPIRAL = 1 ; IF (NONL_S%LSPIRAL) ISPIRAL = ISPINOR+1
__GPU         IBASE = ISPINOR * NGVECTOR

              LMBASE = NONL_S%LMBASE(NI) + ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)

GPACC !$ACC LOOP COLLAPSE(2)
GPOMP !$OMP PARALLEL DO COLLAPSE(2) IF(OFFLOAD_ON)
              DO LM = 1, NONL_S%LMMAX(NONL_S%ITYP(NI))
                 DO I = 1, NGVECTOR
                    IGLB = I + IBASE
                    IRANK = (IGLB-1)/NRPLWV_RED + 1
                    ILOC = MOD(IGLB-1, NRPLWV_RED) + 1
                    IF (IRANK == NODE_ME) THEN
                       NONL_S%CPROJK(ILOC, LM+LMBASE) = &
                          CONJG( NONL_S%QPROJ(I, LM, NONL_S%ITYP(NI), NK, ISPIRAL) * NONL_S%CREXP(I, NI) * NONL_S%CQFAK(LM, NONL_S%ITYP(NI)) )
                    ENDIF
                 ENDDO
              ENDDO

           ENDDO ion

        ENDDO spinor

        PROFILING_STOP('cprojk')

        RETURN
      END SUBROUTINE CPROJK
#endif // cprojk

!************************ SUBROUTINE HSPSI *****************************
!
!> Compute the action of the Hamiltonian and overlap operators on a
!> strip of wave functions W.
!> The result is stored in HW%CW and SW%CW, respectively.
!> On exit W%CW, HW%CW, and SW%CW will all be in over-plane-wave
!> distribution.
!
!***********************************************************************

      SUBROUTINE HSPSI(HAMILTONIAN, SV, CDIJ, CQIJ, ISP, NONLR_S, NONL_S, LATT_CUR, EDCHF, LHF, WDES1, W, HW, SW, NBSTART, NBSTOP)

__GPU   USE moffload_struct_def
GPACC   USE moffload, ONLY : ACC_ATTACH_ASYNC, ACC_DETACH_ASYNC, ACC_SET_ASYNC_Q, ACC_SYNC_ASYNC_Q

        USE wave_struct_def
        USE hamil_struct_def
        USE nonl_struct_def
        USE nonlr_struct_def
        USE poscar_struct_def

        USE wave_high, ONLY : FFTWAV_W1, ELEMENTS
        USE nonl_high, ONLY : W1_PROJALL
        USE hamil, ONLY : HAMILTMU
        USE wave_mpi, ONLY : LASYNC
        USE fock_ace, ONLY : LFOCK_ACE

        USE tutor, ONLY : vtutor

        IMPLICIT NONE

        !> Contains the magnetic vector potential and the derivative
        !> of the xc-energy density w.r.t. the kinetic energy density
        TYPE(ham_handle) :: HAMILTONIAN
        !> local potential
        RGRID :: SV(:,:)
        !> PAW strength parameters and augmentation charges
        OVERLAP :: CDIJ(:,:,:,:), CQIJ(:,:,:,:)
        !> Spin channel
        INTEGER :: ISP
         !> real space projection operators
        TYPE (nonlr_struct) :: NONLR_S
        !> reciprocal space projection operators
        TYPE (nonl_struct) :: NONL_S

        !> In case we use a meta-GGA we need the reciprocal lattive vectors
        TYPE(latt) :: LATT_CUR

        !> Orbitals and action of Hamiltonian and overlap operator on the orbitals
        TYPE(wavefuna) :: W, HW, SW
        !> Common descriptor for H, HW, and SW
        TYPE(wavedes1), TARGET :: WDES1

        !> Fock contribution to double counting energy
        REAL(q) :: EDCHF

        !> Do we need to add Fock contributions?
        LOGICAL :: LHF

        !> Consider bands [NBSTART .. NBSTOP] (global indices)
        INTEGER, OPTIONAL :: NBSTART, NBSTOP

        ! local variables
        TYPE(wavefun1), ALLOCATABLE :: W1(:)

        TYPE(buffer) :: BUF(2)

        REAL(q), ALLOCATABLE :: ZEROS(:)
        INTEGER :: NBSTART_LOCAL, NBSTOP_LOCAL, NBSTART_GLOBAL, NBSTOP_GLOBAL
        INTEGER :: NBANDS, NSTRIP, NSTRIPS, NB, NSTRIP_ACT, ISTRIP, I, MPLWV
        INTEGER :: IBUF, IBUF_PREV, IBUF_NEXT

        LOGICAL :: LASYNC_LOCAL

GPACC2  INTEGER :: IQ

        PROFILING_START('hspsi')

        LASYNC_LOCAL = LASYNC
__GPU   IF (OFFLOAD_ON) LASYNC_LOCAL = .FALSE.

        NBSTART_LOCAL = 1; NBSTART_GLOBAL = 1
        IF (PRESENT(NBSTART)) THEN
           NBSTART_GLOBAL = NBSTART
           NBSTART_LOCAL = (NBSTART_GLOBAL+WDES1%NB_PAR-1)/WDES1%NB_PAR
        ENDIF

        NBSTOP_LOCAL = WDES1%NBANDS ; NBSTOP_GLOBAL = WDES1%NB_TOT
        IF (PRESENT(NBSTOP)) THEN
           NBSTOP_GLOBAL = NBSTOP
           NBSTOP_LOCAL = (NBSTOP_GLOBAL+WDES1%NB_PAR-1)/WDES1%NB_PAR
        ENDIF

        ! quick return if possible
        IF (NBSTART_GLOBAL>NBSTOP_GLOBAL) THEN
           PROFILING_STOP('hspsi')
           RETURN
        ENDIF

        NBANDS = NBSTOP_LOCAL-NBSTART_LOCAL+1

        NSTRIP = WDES1%NSIM
        NSTRIP = MIN(NSTRIP, NBANDS)
!@TODO: choice of NSTRIP needs to be evaluated and optimized

        ! the total number of strips
        NSTRIPS = (NBANDS+NSTRIP-1)/NSTRIP

        ALLOCATE(W1(NSTRIP))
GPACC !$ACC ENTER DATA CREATE(W1(:)) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W1) IF(OFFLOAD_ON)
        MPLWV=WDES1%GRID%MPLWV*WDES1%NRSPINORS
        DO I = 1, NSTRIP
           W1(I)%WDES1 => WDES1
GPACC      IF (OFFLOAD_ON) CALL ACC_ATTACH_ASYNC(W1(I)%WDES1, ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:W1(I)%WDES1) IF(OFFLOAD_ON)
#ifdef OFFLOADING
#ifdef CRAYWAR3
!$OMP TARGET IF(OFFLOAD_ON)
          W1(I)%WDES1 => WDES1
!$OMP END TARGET
#endif // CRAYWAR3
#endif // OFFLOADING
           ALLOCATE(W1(I)%CR(MPLWV))
           NULLIFY(W1(I)%CPROJ)
           NULLIFY(W1(I)%CW)
GPACC !$ACC ENTER DATA CREATE(W1(I)%CR) __IF_ASYNC__
#ifdef INTELMKL
! For the IFX compiler, a derived type pointer's size on the device is only correct
! if the pointer is put on the device with ALWAYS,ALLOC mapping clause.
! Usually this is handled by NEWWAVA in wave_high.F
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:W1(I)%CR) IF(OFFLOAD_ON)
#else
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W1(I)%CR) IF(OFFLOAD_ON)
#endif
        ENDDO

#ifdef wa_proj_cprojk
        IF (.NOT.NONLR_S%LREAL) CALL WA_PROJ_CPROJK(NONL_S, W, NBSTART_GLOBAL, NBSTOP_GLOBAL)
#endif

        ! create the asynchronous communication buffers
        CALL CREATE_BUF(WDES1, NSTRIP, LASYNC_LOCAL, BUF(1), 4)
        CALL CREATE_BUF(WDES1, NSTRIP, LASYNC_LOCAL, BUF(2), 4)

        ! In case the wave functions are in over-plave-wave distribution
        ! we need to redistribute them (WDES1%OVER_BAND = .TRUE.)
        IF (W%OVER_BAND) THEN
           ! redistribute the wave functions in strip 1 (to default distribution)
           CALL REDIS(BUF(1), W, NBSTART_LOCAL, NBSTART_LOCAL+NSTRIP-1, 1)
           ! and wait for the redistribution to finish
           CALL WAIT(BUF(1), W, 1)
        ENDIF

        ALLOCATE(ZEROS(NSTRIP)) ; ZEROS = 0._q

        ! initially the action of the Hamltonian and overlap operator
        ! on the wave functions is in default dsitribution
        HW%OVER_BAND = .FALSE. ; SW%OVER_BAND = .FALSE.

GPACC2  IQ = NSTRIP + 1
GPACC2  IF (OFFLOAD_ON .AND. LUSENCCL) SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)

        blocks: DO NB = NBSTART_LOCAL, NBSTOP_LOCAL, NSTRIP
           ! the current strip number
           ISTRIP = (NB-NBSTART_LOCAL)/NSTRIP + 1
           ! the buffer for this strip
           IBUF = MODULO(ISTRIP+1,2) + 1
           ! the actual number of orbitals in this strip
           NSTRIP_ACT = MIN(NSTRIP,NBSTOP_LOCAL-NB+1)
         !   write(*,'(3(A,I3))') 'block: ',NB, ' NSTRIP_ACT: ',NSTRIP_ACT, ' NSTRIP - NSTRIP_ACT: ',NSTRIP-NSTRIP_ACT

GPACC2     IF (OFFLOAD_ON .AND. LUSENCCL) SET_ACC_ASYNC_Q(IQ)

           IF (W%OVER_BAND .AND. ISTRIP < NSTRIPS) THEN
              IBUF_NEXT = MODULO(ISTRIP+2,2) + 1
              ! redistribute the strip of orbitals ISTRIP+1 (to default distribution)
              CALL REDIS(BUF(IBUF_NEXT), W, NB+NSTRIP, MIN(NBSTOP_LOCAL, NB+2*NSTRIP-1), 1)
           ENDIF

           IF (ISTRIP > 1) THEN
              IBUF_PREV = MODULO(ISTRIP,2) + 1
              ! redistribute the strip of orbitals ISTRIP-1, and the action of the Hamiltonian
              ! and overlap on the aforeementioned (to over-plane-wave distribution)
              CALL REDIS(BUF(IBUF_PREV), W,  NB-NSTRIP, NB-1, 4)
              CALL REDIS(BUF(IBUF_PREV), HW, NB-NSTRIP, NB-1, 2)
              CALL REDIS(BUF(IBUF_PREV), SW, NB-NSTRIP, NB-1, 3)
           ENDIF

GPACC2     IF (OFFLOAD_ON .AND. LUSENCCL) SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)
           CALL SETSTRIP(W, NB, NB+NSTRIP_ACT-1, W1)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
           ! FFT strip of orbitals in buffer ibuf to real space
           DO I = 1, NSTRIP_ACT
              IF (W1(I)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(I)
                 CALL FFTWAV_W1(W1(I))
              ENDIF
           ENDDO

           ! compute the projections of the orbitals in the current strip on the PAW projectors
!@TODO: absorb into HPSIMU
#ifdef wa_proj_cprojk
           IF (NONLR_S%LREAL) THEN
#endif
GPACC         IF (OFFLOAD_ON) CALL ACC_SYNC_ASYNC_Q(SIZE(W1),W1(:)%LDO)
              SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)
              CALL W1_PROJALL(WDES1, W1, NONLR_S, NONL_S, NSTRIP_ACT)
#ifdef wa_proj_cprojk
           ENDIF
#endif
           ! compute the action of the Hamiltonian on orbitals strip ISTRIP in buffer IBUF
           IF (ASSOCIATED(HAMILTONIAN%AVEC)) THEN
!@TODO: add support for the vector potential
              CALL vtutor%error("HSPSI: the vector potential is not implemented &
                                &for the non-blocked Davidson minimizer yet, sorry!")
           ELSE IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
              CALL HPSIMU_TAU( &
                 WDES1, W1, NONLR_S, NONL_S, ZEROS, CDIJ, SV, ISP, ELEMENTS(HW, NB, NB+NSTRIP_ACT-1), &
                 LATT_CUR, HAMILTONIAN%MU)
           ELSE
              CALL HPSIMU( &
                 WDES1, W1, NONLR_S, NONL_S, ZEROS, CDIJ, SV, ISP, ELEMENTS(HW, NB, NB+NSTRIP_ACT-1))
           ENDIF

           ! compute the action of the overlap operator on orbital strip ISTRIP in buffer IBUF
#ifdef wa_vnlacc_cprojk
           IF (NONLR_S%LREAL) &
#endif
           CALL SPSIMU(WDES1, W1, NONLR_S, NONL_S, ZEROS, CQIJ, ELEMENTS(SW, NB, NB+NSTRIP_ACT-1))

           IF (W%OVER_BAND .AND. ISTRIP < NSTRIPS) THEN
              ! wait for the redistribution of the strip of orbitals ISTRIP+1
              ! (to default distribution) to finish
              CALL WAIT(BUF(IBUF_NEXT), W, 1)
           ENDIF

           IF (ISTRIP > 1) THEN
              ! wait for the redistribution of orbitals ISTRIP-1, and the action of the Hamiltonian
              ! and overlap operators on the aforementioned (to over-plane-wave distribution) to finish
              CALL WAIT(BUF(IBUF_PREV), W,  4)
              CALL WAIT(BUF(IBUF_PREV), HW, 2)
              CALL WAIT(BUF(IBUF_PREV), SW, 3)
           ENDIF

GPACC2 !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON .AND. LUSENCCL)
        ENDDO blocks

        ! redistribute the orbitals of the last strip, and the action of the Hamiltonian
        ! and overlap operators on the aforementioned (to over-plane-wave distribution)
        CALL REDIS(BUF(IBUF), W,  NBSTOP_LOCAL-NSTRIP_ACT+1, NBSTOP_LOCAL, 4)
        CALL REDIS(BUF(IBUF), HW, NBSTOP_LOCAL-NSTRIP_ACT+1, NBSTOP_LOCAL, 2)
        CALL REDIS(BUF(IBUF), SW, NBSTOP_LOCAL-NSTRIP_ACT+1, NBSTOP_LOCAL, 3)

        ! wait for all remaining redistribution in BUF(IBUF) to finish
        CALL WAIT(BUF(IBUF), W,  4)
        CALL WAIT(BUF(IBUF), HW, 2)
        CALL WAIT(BUF(IBUF), SW, 3)

!!GPACC2 !$ACC WAIT IF(OFFLOAD_ON .AND. LUSENCCL)

        ! throw away the buffers
        CALL DESTROY_BUF(BUF(1))
        CALL DESTROY_BUF(BUF(2))

        ! After this point everything (orbitals, Hamiltonian, and overlap)
        ! is in over-plane-wave distribution
        ! Paradoxically this means we need to set WDES1%OVER_BAND = .TRUE. ...
        W%OVER_BAND = .TRUE. ; HW%OVER_BAND = .TRUE. ; SW%OVER_BAND = .TRUE.

        ! free workspace
        DO I = 1, SIZE(W1)
           IF (ASSOCIATED(W1(I)%WDES1)) THEN
GPACC         IF (OFFLOAD_ON) CALL ACC_DETACH_ASYNC(W1(I)%WDES1, ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1(I)%WDES1) IF(OFFLOAD_ON)
              NULLIFY(W1(I)%WDES1)
           ENDIF
           IF (ASSOCIATED(W1(I)%CW)) THEN
GPACC         IF (OFFLOAD_ON) CALL ACC_DETACH_ASYNC(W1(I)%CW, ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1(I)%CW) IF(OFFLOAD_ON)
              NULLIFY(W1(I)%CW)
           ENDIF
           IF (ASSOCIATED(W1(I)%CPROJ)) THEN
GPACC         IF (OFFLOAD_ON) CALL ACC_DETACH_ASYNC(W1(I)%CPROJ, ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1(I)%CPROJ) IF(OFFLOAD_ON)
              NULLIFY(W1(I)%CPROJ)
           ENDIF
           IF (ASSOCIATED(W1(I)%CR)) THEN
GPACC !$ACC EXIT DATA DELETE(W1(I)%CR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1(I)%CR) IF(OFFLOAD_ON)
              DEALLOCATE(W1(I)%CR) ; NULLIFY(W1(I)%CR)
           ENDIF
        ENDDO
GPACC !$ACC EXIT DATA DELETE(W1(:)) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1) IF(OFFLOAD_ON)
        DEALLOCATE(W1)

        DEALLOCATE(ZEROS)

#ifdef wa_vnlacc_cprojk
        ! In the reciprocal-space-projector case we can compute the action of
        ! the non-local part of the Hamiltonian and overlap operators here
        IF (.NOT.NONLR_S%LREAL) THEN
           ! Beware the first orbital treated in the blocks-loop above has
           ! global index (NBSTART_LOCAL-1)*WDES1%NB_PAR+1, this may be smalles
           ! than NBSTART_GLOBAL ... just reset NBSTART_GLOBAL to take care of this
           NBSTART_GLOBAL = (NBSTART_LOCAL -1) * WDES1%NB_PAR + 1

           ! < G | HW_n > = < G | HW_n > + \sum_i < G | p_i > D_ij < p_j | \psi_n >
           CALL WA_VNLACC_CPROJK(NONL_S, W, HW, CDIJ, ISP, NBSTART_GLOBAL, NBSTOP_GLOBAL)
           ! < G | SW_n > =  < G | \psi_n >

GPACC !$ACC KERNELS PRESENT(SW%CW_RED) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:SW%CW_RED) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
           SW%CW_RED(:,NBSTART_GLOBAL:NBSTOP_GLOBAL) = W%CW_RED(:,NBSTART_GLOBAL:NBSTOP_GLOBAL)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
           ! < G | SW_n > = < G | SW_n > + \sum_i < G | p_i > Q_ij < p_j | \psi_n >
           CALL WA_VNLACC_CPROJK(NONL_S, W, SW, CQIJ, ISP, NBSTART_GLOBAL, NBSTOP_GLOBAL)
        ENDIF
#endif

        ! For hybrid functionals
        EDCHF = 0
        IF (LHF .AND. LFOCK_ACE()) THEN
           ! Add action of the ACE on the orbitals W, and compute the Fock contribution
           ! to the double counting energy
           CALL WA_FOCK_ACE_APPLY_ACC(W, HW, EDCHF, ISP, NBSTART_GLOBAL, NBSTOP_GLOBAL)
        ENDIF

        PROFILING_STOP('hspsi')

        CONTAINS

!************************ SUBROUTINE SETSTRIP **************************
!
!***********************************************************************

        SUBROUTINE SETSTRIP(W, ISTART, ISTOP, W1)

__GPU     USE moffload_struct_def
GPACC     USE moffload, ONLY : ACC_ATTACH_ASYNC
          USE wave_struct_def

          IMPLICIT NONE

          TYPE(wavefuna) :: W
          TYPE(wavefun1) :: W1(:)
          INTEGER :: ISTART,  ISTOP

          ! local variables
          INTEGER :: I, IGLOBAL

          PROFILING_START('setstrip')

          ! Set the single wave function pointers
          DO I = 1, ISTOP-ISTART+1
             IF (ASSOCIATED(W1(I)%CPROJ)) THEN
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1(I)%CPROJ) IF(OFFLOAD_ON)
               NULLIFY(W1(I)%CPROJ)
             END IF

             IF (ASSOCIATED(W1(I)%CW)) THEN
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1(I)%CW) IF(OFFLOAD_ON)
               NULLIFY(W1(I)%CW)
             END IF

             W1(I)%CW => W%CW(:,ISTART+I-1)
GPACC        IF (OFFLOAD_ON) CALL ACC_ATTACH_ASYNC(W1(I)%CW, ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:W1(I)%CW) IF(OFFLOAD_ON)          

             W1(I)%CPROJ => W%CPROJ(:,ISTART+I-1)
GPACC        IF (OFFLOAD_ON) CALL ACC_ATTACH_ASYNC(W1(I)%CPROJ, ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC:W1(I)%CPROJ) IF(OFFLOAD_ON)
             W1(I)%LDO = .TRUE.
#ifdef CRAYWAR3
!$OMP TARGET IF(OFFLOAD_ON)
              W1(I)%CW => W%CW(:,ISTART+I-1)
              W1(I)%CPROJ => W%CPROJ(:,ISTART+I-1)
!$OMP END TARGET
#endif // CRAYWAR3
          ENDDO

          ! The last strip might not be full so we need to mark those
          ! entries that should not be touched by HAMILTMU etc
          DO I = ISTOP-ISTART+2, SIZE(W1)
             W1(I)%LDO = .FALSE.
          ENDDO

          PROFILING_STOP('setstrip')

          RETURN
        END SUBROUTINE SETSTRIP

      END SUBROUTINE HSPSI


!************************ SUBROUTINE HPSIMU ****************************
!
!> Compute the action of the Hamiltonian on a strip of orbitals.
!> This is an almost literal copy (ughh, yes I know) of HAMILTMU:
!> the only difference is that the evaluation of the non-local part
!> of the Hamiltonian in reciprocal space may be skipped when the
!> code is compiled with -Dwa_vnlacc_cprojk
!
!@TODO: find a better solution than (almost) copying HAMILTMU
!@TODO: the code below can be simplified if support for
!       -Uwa_vnlacc_cprojk would be dropped
!
!***********************************************************************

      SUBROUTINE HPSIMU(WDES1, W1, NONLR_S, NONL_S, EVALUE, CDIJ, SV, ISP, HWA)

__GPU   USE moffload_struct_def
GPACC   USE moffload, ONLY : ACC_SET_ASYNC_Q,  ACC_SYNC_ASYNC_Q

        USE wave_struct_def
        USE nonlr_struct_def
        USE nonl_struct_def

        USE nonl_high, ONLY : RACCMU_, VNLACC

        IMPLICIT NONE

        TYPE(wavedes1) :: WDES1
        TYPE(wavefun1) :: W1(:)
        TYPE(wavefuna) :: HWA
        TYPE(nonlr_struct) :: NONLR_S
        TYPE(nonl_struct) :: NONL_S

        REAL(q) :: EVALUE(:)
        OVERLAP :: CDIJ(:,:,:,:)
        RGRID :: SV(:,:)

        INTEGER :: ISP

        ! local variables
        COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS, SIZE(W1))
        INTEGER NP, N

        PROFILING_START('hpsimu')

GPACC !$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

        DO NP = 1, SIZE(W1)
           IF (W1(NP)%LDO) THEN
GPACC         IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
              CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
           ENDIF
        ENDDO

        IF (NONLR_S%LREAL) THEN
           CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CDIJ, ISP, EVALUE, CWORK1)
           DO NP = 1, SIZE(W1)
              IF (W1(NP)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
                 CALL KINHAMIL(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                      WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), HWA%CW(1,NP))
              ENDIF
           ENDDO
        ELSE
           DO NP = 1, SIZE(W1)
              IF (W1(NP)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
#ifndef wa_vnlacc_cprojk
                 CALL VNLACC(NONL_S, W1(NP), CDIJ, CDIJ, ISP, EVALUE(NP), HWA%CW(:,NP))
                 CALL KINHAMIL(WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                      WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), HWA%CW(1,NP))
#else
                 CALL KINHAMIL(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                      WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CW(1), HWA%CW(1,NP))
#endif
              ENDIF
           ENDDO
        ENDIF

GPACC   IF (OFFLOAD_ON) CALL ACC_SYNC_ASYNC_Q(SIZE(W1),W1(:)%LDO)
        SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)
GPACC !$ACC EXIT DATA DELETE(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

        PROFILING_STOP('hpsimu')

        RETURN
      END SUBROUTINE HPSIMU


!************************ SUBROUTINE HPSIMU_TAU ************************
!
!> Compute the action of the Hamiltonian on a strip of orbitals.
!> This is an almost literal copy (ughh, yes I know) of HAMILTMU_TAU:
!> the only difference is that the evaluation of the non-local part
!> of the Hamiltonian in reciprocal space may be skipped when the
!> code is compiled with -Dwa_vnlacc_cprojk
!
!@TODO: find a better solution than (almost) copying HAMILTMU_TAU
!@TODO: the code below can be simplified if support for
!       -Uwa_vnlacc_cprojk would be dropped
!
!***********************************************************************

      SUBROUTINE HPSIMU_TAU(WDES1, W1, NONLR_S, NONL_S, EVALUE, CDIJ, SV, ISP, HWA, LATT_CUR, MU)

__GPU   USE moffload_struct_def
GPACC   USE moffload, ONLY : ACC_SET_ASYNC_Q,  ACC_SYNC_ASYNC_Q

        USE wave_struct_def
        USE nonlr_struct_def
        USE nonl_struct_def
        USE poscar_struct_def

        USE nonl_high, ONLY : RACCMU_, VNLACC

        IMPLICIT NONE

        TYPE(wavedes1) :: WDES1
        TYPE(wavefun1) :: W1(:)
        TYPE(wavefuna) :: HWA
        TYPE(nonlr_struct) :: NONLR_S
        TYPE(nonl_struct) :: NONL_S
        TYPE(latt) :: LATT_CUR

        REAL(q) :: EVALUE(:)
        OVERLAP :: CDIJ(:,:,:,:)
        RGRID :: SV(:,:)
        RGRID :: MU(:,:)

        INTEGER :: ISP

        ! local variables
        COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS, SIZE(W1))
#ifdef OFFLOADING
        COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV),SIZE(W1))
        COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV),SIZE(W1))
#else
        COMPLEX(q) :: CWORK2(__align_DIM(WDES1%NRPLWV)     __omp_arg(3))
        COMPLEX(q) :: CWORK3(__align_DIM(WDES1%GRID%MPLWV) __omp_arg(6))
#endif

        INTEGER NP, N

        PROFILING_START('hpsimu_tau')

GPACC !$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

        DO NP = 1, SIZE(W1)
           IF (W1(NP)%LDO) THEN
GPACC         IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
              CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
           ENDIF
        ENDDO

        IF (NONLR_S%LREAL) THEN
           CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CDIJ, ISP, EVALUE, CWORK1)
           DO NP = 1, SIZE(W1)
              IF (W1(NP)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
                 CALL KINHAMIL_TAU(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., .TRUE., &
                      WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef OFFLOADING
                      CWORK2(1,NP), CWORK3(1,NP), &
#else
                      CWORK2(1 __omp_arg(1)), CWORK3(1 __omp_arg(1)), &
#endif
                      EVALUE(NP), W1(NP)%CW(1), HWA%CW(1,NP))
              ENDIF
           ENDDO
        ELSE
           DO NP = 1, SIZE(W1)
              IF (W1(NP)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
#ifndef wa_vnlacc_cprojk
                 CALL VNLACC(NONL_S, W1(NP), CDIJ, CDIJ, ISP, EVALUE(NP), HWA%CW(:,NP))
                 CALL KINHAMIL_TAU(WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., .TRUE., &
                      WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef OFFLOADING
                      CWORK2(1,NP), CWORK3(1,NP), &
#else
                      CWORK2(1 __omp_arg(1)), CWORK3(1 __omp_arg(1)), &
#endif
                      EVALUE(NP), W1(NP)%CW(1), HWA%CW(1,NP))

#else
                 CALL KINHAMIL_TAU(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., .TRUE., &
                      WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
#ifdef OFFLOADING
                      CWORK2(1,NP), CWORK3(1,NP), &
#else
                      CWORK2(1 __omp_arg(1)), CWORK3(1 __omp_arg(1)), &
#endif
                      EVALUE(NP), W1(NP)%CW(1), HWA%CW(1,NP))

#endif
              ENDIF
           ENDDO
        ENDIF

GPACC   IF (OFFLOAD_ON) CALL ACC_SYNC_ASYNC_Q(SIZE(W1),W1(:)%LDO)
        SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)
GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK2,CWORK3,LATT_CUR%B,LATT_CUR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

        PROFILING_STOP('hpsimu_tau')

        RETURN
      END SUBROUTINE HPSIMU_TAU


!************************ SUBROUTINE SPSIMU ****************************
!
!> Compute the action of the overlap operator on a strip of orbitals
!
!***********************************************************************

      SUBROUTINE SPSIMU(WDES1, W1, NONLR_S, NONL_S, EVALUE, CQIJ, SWA)

__GPU   USE moffload_struct_def
GPACC   USE moffload, ONLY : ACC_SET_ASYNC_Q, ACC_SYNC_ASYNC_Q

        USE wave_struct_def
        USE nonlr_struct_def
        USE nonl_struct_def

        USE nonl_high, ONLY : RACCMU_, VNLACC_ADD

        IMPLICIT NONE

        TYPE(wavedes1) :: WDES1
        TYPE(wavefun1) :: W1(:)
        TYPE(wavefuna) :: SWA
        TYPE(nonlr_struct) :: NONLR_S
        TYPE(nonl_struct) :: NONL_S

        REAL(q) :: EVALUE(:)
        OVERLAP :: CQIJ(:,:,:,:)

        ! local variables
        COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS, SIZE(W1))
        INTEGER NP, N

        PROFILING_START('spsimu')

GPACC !$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK1) IF(OFFLOAD_ON)

        IF (NONLR_S%LREAL) THEN
GPACC !$ACC KERNELS PRESENT(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
           CWORK1 = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
           CALL RACCMU_(NONLR_S, WDES1, W1, CQIJ, CQIJ, 1, EVALUE, CWORK1)
           DO NP = 1, SIZE(W1)
              IF (W1(NP)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
                 CALL FFTHAMIL(WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., -1._q, W1(NP)%CW(1), SWA%CW(1,NP))
              ENDIF
           ENDDO
        ELSE
           DO NP = 1, SIZE(W1)
              IF (W1(NP)%LDO) THEN
GPACC            IF (OFFLOAD_ON) CALL ACC_SET_ASYNC_Q(NP)
                 CALL WZCOPY(WDES1%NRPLWV, W1(NP)%CW(1), 1, SWA%CW(1,NP), 1)
                 CALL VNLACC_ADD(NONL_S, W1(NP), CQIJ, CQIJ, 1, EVALUE(NP), SWA%CW(:,NP))
              ENDIF
           ENDDO
        ENDIF

GPACC   IF (OFFLOAD_ON) CALL ACC_SYNC_ASYNC_Q(SIZE(W1),W1(:)%LDO)
        SET_ACC_ASYNC_Q(ACC_ASYNC_ASYNC)
GPACC !$ACC EXIT DATA DELETE(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK1) IF(OFFLOAD_ON)

        PROFILING_STOP('spsimu')

        RETURN
      END SUBROUTINE SPSIMU

#ifdef wa_proj_cprojk
!************************ SUBROUTINE WA_PROJ_CPROJK ********************
!
!***********************************************************************

      SUBROUTINE WA_PROJ_CPROJK(NONL_S, WA, NBSTART, NBSTOP)

__GPU   USE moffload_struct_def

        USE nonl_struct_def
        USE wave_struct_def

        USE tutor, ONLY : vtutor

        IMPLICIT NONE

        !> reciprocal space projection operators
        TYPE(nonl_struct) :: NONL_S
        !> orbitals
        TYPE(wavefuna) :: WA
        !> project bands [NBSTART .. NBSTOP] (global indices) on the projectors
        INTEGER :: NBSTART, NBSTOP

        ! local variables
        COMPLEX(q), ALLOCATABLE :: CPROJ(:,:)
        INTEGER :: NPL_RED, NPRO_TOT, NRPLWV_RED, NSTRIP, NPOS, NSTRIP_ACT
#ifdef over_projections_distribution
        INTEGER :: NPROD_RED, NPRO_BASE, NPRO_RED
#else
        INTEGER :: N, NB_GLOBAL, NB_LOCAL, NPRO, NI, NT, LMMAXC, L, NPRO_POS
#endif

        PROFILING_START('wa_proj_cprojk')

        IF (.NOT.WA%OVER_BAND) &
           CALL vtutor%bug("WA_PROJ_CPROJK: WA is not distributed over plane-waves", __FILE__, __LINE__)

        ! the number of plane-wave components this rank owns in "over-plane-wave" distribution
        NPL_RED = WA%WDES1%NPL_RED
        ! the total number of projections
        NPRO_TOT = WA%WDES1%NPRO_TOT
        ! maximum number of plane waves any rank owns in "over-plane-wave" distribution
        NRPLWV_RED = WA%WDES1%NRPLWV_RED

#ifdef over_projections_distribution
        ! the number of projections this rank owns in "over-projections" distribution
        NPRO_RED = WA%WDES1%NPRO_RED

        NPROD_RED = WA%WDES1%NPROD_RED
        NPRO_BASE = (WA%WDES1%COMM_KIN%NODE_ME-1)*NPROD_RED
#endif

        NSTRIP = NBSTOP-NBSTART+1
!@TODO: strip size needs to be evaluated and optimized

        ALLOCATE(CPROJ(NPRO_TOT, NSTRIP))
GPACC !$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ) IF(OFFLOAD_ON)

        DO NPOS = NBSTART, NBSTOP, NSTRIP
           NSTRIP_ACT = MIN(NSTRIP, NBSTOP-NPOS+1)

           PROFILING_START('wa_proj_cprojk|gemm')

           CALL WZGEMM('C', 'N', NPRO_TOT, NSTRIP_ACT, NPL_RED, &
                      zone, NONL_S%CPROJK(1,1), NRPLWV_RED, WA%CW_RED(1,NPOS), NRPLWV_RED, &
                      zzero, CPROJ(1,1), NPRO_TOT)

           PROFILING_STOP('wa_proj_cprojk|gemm')

           IF (WA%WDES1%COMM_KIN%NCPU > 1) THEN
              CALLMPI( M_sum_z(WA%WDES1%COMM_KIN, CPROJ(1,1), NPRO_TOT*NSTRIP_ACT))
           ENDIF
#ifdef over_projections_distribution
           ! store directly in "over-projections" distribution
GPACC !$ACC KERNELS PRESENT(WA%CPROJ_RED,CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:WA%CPROJ_RED,CPROJ) IF(OFFLOAD_ON)
           IF (NPRO_RED > 0) &
              WA%CPROJ_RED(1:NPRO_RED, NPOS:NPOS+NSTRIP_ACT-1) = CPROJ(NPRO_BASE+1:NPRO_BASE+NPRO_RED, 1:NSTRIP_ACT)
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#else
           istrip: DO N = 1, NSTRIP_ACT
              NB_GLOBAL = NBSTART+N-1
              ! should this rank store the projection for orbital NB_GLOBAL?
              IF (MOD(NB_GLOBAL-1,WA%WDES1%NB_PAR)+1 == WA%WDES1%NB_LOW) THEN
                 ! this is the local band index
                 NB_LOCAL = (NB_GLOBAL-1)/WA%WDES1%NB_PAR + 1

                 ! in case the projections of orbital NB_LOCAL are fully owned by a single rank
                 IF (WA%WDES1%COMM_INB%NCPU==1) THEN
GPACC !$ACC KERNELS PRESENT(WA,CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:WA,CPROJ) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
                    WA%CPROJ(1:WA%WDES1%NPRO, NB_LOCAL) = CPROJ(1:WA%WDES1%NPRO, N)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
                    CYCLE istrip
                 ENDIF

                 ! in case the projections of orbital NB_LOCAL are distributed over COMM_INB
GPACC !$ACC PARALLEL LOOP PRESENT(WA,CPROJ) PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WA,CPROJ) &
GPOMP !$OMP PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) IF(OFFLOAD_ON)
                 DO NI = 1, WA%WDES1%NIONS
                    NT = WA%WDES1%ITYP(NI)
                    LMMAXC = WA%WDES1%LMMAX(NT)
                    IF (LMMAXC/=0) THEN
                       NPRO_POS = WA%WDES1%NPRO_POS(NI)
                       NPRO = WA%WDES1%LMBASE(NI)
GPACC !$ACC LOOP
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
                       DO L = 1, LMMAXC
                          WA%CPROJ(L+NPRO, NB_LOCAL) = CPROJ(L+NPRO_POS, N)
                       ENDDO
                    ENDIF
                 ENDDO

                 IF (.NOT.WA%WDES1%LNONCOLLINEAR) CYCLE istrip

GPACC !$ACC PARALLEL LOOP PRESENT(WA,CPROJ) PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WA,CPROJ) &
GPOMP !$OMP PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) IF(OFFLOAD_ON)
                 DO NI = 1, WA%WDES1%NIONS
                    NT = WA%WDES1%ITYP(NI)
                    LMMAXC = WA%WDES1%LMMAX(NT)
                    IF (LMMAXC/=0) THEN
                       NPRO_POS = WA%WDES1%NPRO_POS(NI) + NPRO_TOT/2
                       NPRO = WA%WDES1%LMBASE(NI) + WA%WDES1%LMBASE(WA%WDES1%NIONS+1)
GPACC !$ACC LOOP
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
                       DO L = 1, LMMAXC
                          WA%CPROJ(L+NPRO, NB_LOCAL) = CPROJ(L+NPRO_POS, N)
                       ENDDO
                    ENDIF
                 ENDDO

              ENDIF
           ENDDO istrip
#endif
        ENDDO

GPACC !$ACC EXIT DATA DELETE(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CPROJ) IF(OFFLOAD_ON)
        DEALLOCATE(CPROJ)

        PROFILING_STOP('wa_proj_cprojk')

        RETURN
      END SUBROUTINE WA_PROJ_CPROJK
#endif // wa_proj_cprojk
#ifdef wa_vnlacc_cprojk
!************************ SUBROUTINE WA_VNLACC_CPROJK ******************
!
!***********************************************************************

      SUBROUTINE WA_VNLACC_CPROJK(NONL_S, WA, WACC, CIJ, ISP, NBSTART, NBSTOP)

__GPU   USE moffload_struct_def

        USE nonl_struct_def
        USE wave_struct_def

        USE tutor, ONLY : vtutor

        IMPLICIT NONE

        !> reciprocal space projection operators
        TYPE(nonl_struct) :: NONL_S
        !> orbitals
        TYPE(wavefuna) :: WA
        !> action
        TYPE(wavefuna) :: WACC
        !> PAW strength or overlap matrix
        OVERLAP :: CIJ(:,:,:,:)
        !> spin component
        INTEGER :: ISP
        !> add action to bands [NBSTART .. NBSTOP] (global indices)
        INTEGER :: NBSTART, NBSTOP

        ! local variables
        COMPLEX(q), ALLOCATABLE :: CPROJ(:,:)
        GDEF :: GCIJP(WA%WDES1%NPRO)

        INTEGER :: NPOS, NSTRIP, NSTRIP_ACT
        INTEGER :: N, NB_PAR, NB_GLOBAL, NB_LOCAL
        INTEGER :: NPRO_TOT, NPRO, NPRO_POS, NI, NT, LMMAXC, L
        INTEGER :: NPL_RED, NRPLWV_RED

        PROFILING_START('wa_vnlacc_cprojk')

        IF (.NOT.WA%OVER_BAND) &
           CALL vtutor%bug("WA_VMLACC_CPROJK: WA is not distributed over plane-waves", __FILE__, __LINE__)

        IF (.NOT.WACC%OVER_BAND) &
           CALL vtutor%bug("WA_VMLACC_CPROJK: WACC is not distributed over plane-waves", __FILE__, __LINE__)

        NPRO = WA%WDES1%NPRO
        NPRO_TOT = WA%WDES1%NPRO_TOT

        NPL_RED = WA%WDES1%NPL_RED
        NRPLWV_RED = WA%WDES1%NRPLWV_RED

        NSTRIP = NBSTOP-NBSTART+1
!@TODO: strip size needs to be evaluated and optimized

        ALLOCATE(CPROJ(NPRO_TOT, NSTRIP))
GPACC !$ACC ENTER DATA CREATE(CPROJ,GCIJP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ,GCIJP) IF(OFFLOAD_ON)

        DO NPOS = NBSTART, NBSTOP, NSTRIP
           NSTRIP_ACT = MIN(NSTRIP, NBSTOP-NPOS+1)
GPACC !$ACC KERNELS PRESENT(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CPROJ) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
           CPROJ = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET  
GPACC !$ACC END KERNELS
           istrip: DO N = 1, NSTRIP_ACT
              NB_GLOBAL = NBSTART+N-1
              ! is this rank responsible for work on orbital NB_GLOBAL?
              IF (MOD(NB_GLOBAL-1, WA%WDES1%NB_PAR)+1 == WA%WDES1%NB_LOW ) THEN
                 ! this is the local band index
                 NB_LOCAL = (NB_GLOBAL-1)/WA%WDES1%NB_PAR + 1

                 ! calculate GCIJP(i) = \sum_j C_ij < p_j | NB_LOCAL >
                 CALL OVERL1_SIMPLE(WA%WDES1, SIZE(CIJ,1), CIJ(1,1,1,ISP), WA%CPROJ(1,NB_LOCAL), GCIJP(1))

                 ! and store GCIJP at the right place into CPROJ ...
                 ! in case the projections of orbital NB_LOCAL are owned by a single rank
                 IF (WA%WDES1%COMM_INB%NCPU==1) THEN
GPACC !$ACC KERNELS PRESENT(CPROJ,GCIJP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CPROJ,GCIJP) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
                    CPROJ(1:NPRO,N) = GCIJP(1:NPRO)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET  
GPACC !$ACC END KERNELS
                    CYCLE istrip
                 ENDIF

                 ! in case the projections of orbital NB_LOCAL are shared over COMM_INB
GPACC !$ACC PARALLEL LOOP PRESENT(WA%WDES1,CPROJ,GCIJP) PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WA%WDES1,CPROJ,GCIJP) &
GPOMP !$OMP PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) IF(OFFLOAD_ON)
                 DO NI = 1, WA%WDES1%NIONS
                    NT = WA%WDES1%ITYP(NI)
                    LMMAXC = WA%WDES1%LMMAX(NT)
                    IF (LMMAXC/=0) THEN
                       NPRO = WA%WDES1%LMBASE(NI)
                       NPRO_POS = WA%WDES1%NPRO_POS(NI)
GPACC !$ACC LOOP
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
                       DO L = 1, LMMAXC
                          CPROJ(L+NPRO_POS, N) = GCIJP(L+NPRO)
                       ENDDO
                    ENDIF
                 ENDDO

                 IF (.NOT.WA%WDES1%LNONCOLLINEAR) CYCLE istrip

GPACC !$ACC PARALLEL LOOP PRESENT(WA%WDES1,CPROJ,GCIJP) PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WA%WDES1,CPROJ,GCIJP) &
GPOMP !$OMP PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) IF(OFFLOAD_ON)
                 DO NI = 1, WA%WDES1%NIONS
                    NT = WA%WDES1%ITYP(NI)
                    LMMAXC = WA%WDES1%LMMAX(NT)
                    IF (LMMAXC/=0) THEN
                       NPRO = WA%WDES1%LMBASE(NI) + WA%WDES1%LMBASE(WA%WDES1%NIONS+1)
                       NPRO_POS = WA%WDES1%NPRO_POS(NI) + NPRO_TOT/2
GPACC !$ACC LOOP
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
                       DO L = 1, LMMAXC
                          CPROJ(L+NPRO_POS, N) = GCIJP(L+NPRO)
                       ENDDO
                    ENDIF
                 ENDDO

              ENDIF
           ENDDO istrip

           CALLMPI( M_sum_z(WA%WDES1%COMM_KIN, CPROJ(1,1), NPRO_TOT*NSTRIP_ACT))

           PROFILING_START('wa_vnlacc_cprojk|gemm')

           CALL WZGEMM('N', 'N', NPL_RED, NSTRIP_ACT, NPRO_TOT, &
                      zone, NONL_S%CPROJK(1,1), NRPLWV_RED, CPROJ(1,1), NPRO_TOT, &
                      zone, WACC%CW_RED(1,NPOS), NRPLWV_RED)

           PROFILING_STOP('wa_vnlacc_cprojk|gemm')
        ENDDO

GPACC !$ACC EXIT DATA DELETE(CPROJ,GCIJP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CPROJ,GCIJP) IF(OFFLOAD_ON)
        DEALLOCATE(CPROJ)

        PROFILING_STOP('wa_vnlacc_cprojk')

        RETURN
      END SUBROUTINE WA_VNLACC_CPROJK
#endif // wa_vnlacc_cprojk

!************************ SUBROUTINE WA_FOCK_ACE_APPLY_ACC *************
!
!> Apply the adaptively compressed exchange (ACE) to a set of orbitals.
!> The action of the ACE operator on the orbitals WA is added to WACC.
!
!***********************************************************************

      SUBROUTINE WA_FOCK_ACE_APPLY_ACC(WA, WACC, EDCHF, ISP, NBSTART, NBSTOP)

__GPU   USE moffload_struct_def
        USE wave_struct_def

        USE wave_high, ONLY : ELEMENTS
        USE fock_ace, ONLY : WACE

        IMPLICIT NONE

        !> orbitals
        TYPE(wavefuna) :: WA
        !> action
        TYPE(wavefuna) :: WACC
        !> contribution to Fock double counting energy
        REAL(q) :: EDCHF
        !> spin component
        INTEGER :: ISP
        !> add action to bands [NBSTART .. NBSTOP] (global indices)
        INTEGER :: NBSTART, NBSTOP

        ! local variables
        TYPE(wavefuna) :: WX
        COMPLEX(q), ALLOCATABLE :: CACC(:,:)
        GDEF, ALLOCATABLE :: COVL(:,:)
        COMPLEX(q) :: CDCHF
        REAL(q) :: WEIGHTK, WEIGHT
        INTEGER :: NB_ACE, NSTRIP, NRPLWV_RED
        INTEGER :: IB, IG, NG

        PROFILING_START('wa_fock_ace_apply_acc')

        WX = ELEMENTS(WACE, WA%WDES1, ISP)
        NB_ACE = SIZE(WX%CW_RED, 2)

        NSTRIP = NBSTOP-NBSTART+1
        NRPLWV_RED = WACC%WDES1%NRPLWV_RED

        ALLOCATE(COVL(NB_ACE,NSTRIP), CACC(NRPLWV_RED,NSTRIP))
GPACC !$ACC ENTER DATA CREATE(COVL,CACC) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:COVL,CACC) IF(OFFLOAD_ON)

        ! compute COVL(i,j)=< \tilde X_i | \tilde \psi_j >,
        ! for i \in [1,NB_ACE] and j \in [NBSTART,NBSTOP]
        CALL WGGEMM(trans,'N', NB_ACE, NSTRIP, m_ WA%WDES1%NPL_RED, &
       &   one, WX%CW_RED(1,1), m_ WX%WDES1%NRPLWV_RED, WA%CW_RED(1,NBSTART), m_ WA%WDES1%NRPLWV_RED, &
       &   zero, COVL(1,1), NB_ACE)

        CALLMPI( M_sum_g(WACE%WDES%COMM_KIN,COVL(1,1),NB_ACE*NSTRIP))

        ! CACC(:,j) = \sum_i COVL(i,j) WA%CW_RED(:,i)
        CALL WGGEMM('N','N', m_ WA%WDES1%NPL_RED, NSTRIP, NB_ACE, &
       &   one, WX%CW_RED(1,1), m_ WX%WDES1%NRPLWV_RED, COVL(1,1), NB_ACE, &
       &   zero, CACC(1,1), m_ NRPLWV_RED)

GPACC !$ACC EXIT DATA DELETE(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:COVL) IF(OFFLOAD_ON)
        DEALLOCATE(COVL)

        ! add the action of the ACE operator on WA to WACC and
        ! calculate its contribution to the double counting energy
        CDCHF = 0
GPACC !$ACC ENTER DATA COPYIN(CDCHF) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:CDCHF) IF(OFFLOAD_ON)

        WEIGHTK = 0.5_q*WA%WDES1%WTKPT*WA%WDES1%RSPIN
        NG = WA%WDES1%NPL_RED
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WA, CACC, WACC, CDCHF) PRIVATE(WEIGHT) REDUCTION(+:CDCHF) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(WEIGHT) &
GPOMP !$OMP MAP(PRESENT,TO:WA,CACC,WACC,CDCHF) REDUCTION(+:CDCHF) IF(OFFLOAD_ON)
__CPU !$OMP PARALLEL SHARED(NSTRIP, NBSTART, WEIGHTK, NG, WA, CACC, WACC) PRIVATE(IB, WEIGHT, IG) REDUCTION(+:CDCHF)
        DO IB = 1, NSTRIP
__CPU      WEIGHT = WEIGHTK*WA%FERWE(IB+NBSTART-1)
!!           WEIGHT = WEIGHTK !!*WA%FERWE(IB+NBSTART-1)
__CPU !$OMP DO
           DO IG = 1, NG
__GPU         WEIGHT = WEIGHTK*WA%FERWE(IB+NBSTART-1)
              WACC%CW_RED(IG,IB+NBSTART-1) = WACC%CW_RED(IG,IB+NBSTART-1) - CACC(IG,IB)
              CDCHF = CDCHF + CONJG(WA%CW_RED(IG,IB+NBSTART-1))*CACC(IG,IB)*WEIGHT
           ENDDO
__CPU !$OMP END DO
        ENDDO
__CPU !$OMP END PARALLEL

GPACC !$ACC EXIT DATA COPYOUT(CDCHF) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:CDCHF) IF(OFFLOAD_ON)

GPACC !$ACC EXIT DATA DELETE(CACC) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CACC) IF(OFFLOAD_ON)
        DEALLOCATE(CACC)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
        EDCHF = REAL(CDCHF, q)

        PROFILING_STOP('wa_fock_ace_apply_acc')

        RETURN
      END SUBROUTINE WA_FOCK_ACE_APPLY_ACC


!************************ SUBROUTINE BRAKET ****************************
!
!> Compute the upper triangular part of C_ij = < WA1_i | WA2_j >
!> where i, j = [1 .. WDES1%B_TOT]
!> Building this matrix is done in striped fashion, using j-strips of
!> size NSTRIP * NCPU.
!> In case the optional arguments NBSTART and NBSTOP are present only
!> those strips will be calculated that contain j = [NBSTART .. NBSTOP]
!> On exit WA1 and WA2 will be in over-plane-wave distribution.
!
!***********************************************************************

      SUBROUTINE BRAKET(WDES1, WA1, WA2, ISP, COVL, NBSTART, NBSTOP)

__GPU   USE moffload_struct_def

        USE wave_struct_def

        USE dfast, ONLY : NSTRIPD
        USE scala, ONLY : LscaAWARE, INIT_scala
        USE tutor, ONLY : vtutor
        USE string, ONLY : str

        IMPLICIT NONE

        !> Wave functions | WA1 > and | WA2 >
        TYPE(wavefuna) :: WA1, WA2
        !> Common descriptor for WA1 and WA2
        TYPE(wavedes1) :: WDES1
        !> Overlap matrix  COVL(i,j) = < WA1_i | WA2_j >
        GDEF :: COVL(:,:)
        !> Spin channel
        INTEGER :: ISP
        !> Consider WA2 bands [NBSTART .. NBSTOP] (global indices)
        INTEGER, OPTIONAL :: NBSTART, NBSTOP

        INTEGER :: NBANDS, NSTRIP, NB_TOT, NB_RED, NSTRIP_RED, NCPU, NBSTOP_RED, NBSTART_RED

        LOGICAL :: LscaAWARE_LOCAL

        PROFILING_START('braket')

        NBSTART_RED = 1
        IF (PRESENT(NBSTART)) NBSTART_RED = NBSTART

        NBSTOP_RED = WDES1%NB_TOT
        IF (PRESENT(NBSTOP)) NBSTOP_RED = NBSTOP

        ! quick return if possible
        IF (NBSTART_RED>NBSTOP_RED) THEN
           PROFILING_STOP('braket')
           RETURN
        ENDIF

        IF (.NOT.WA1%OVER_BAND .OR. .NOT.WA2%OVER_BAND) THEN
           CALL vtutor%bug("BRAKET: wave functions not in over-plane-wave distribution " &
            // str(.NOT.WA1%OVER_BAND) // " " // str(.NOT.WA2%OVER_BAND), __FILE__, __LINE__)
        ENDIF

        NCPU = WDES1%COMM_INTER%NCPU

        NB_TOT = WDES1%NB_TOT

        NBANDS = NBSTOP_RED-NBSTART_RED+1

!!        NSTRIP = NSTRIP_STANDARD*NCPU
!!        NSTRIP = MIN(NSTRIP, WDES1%NSIM)
        NSTRIP = NSTRIPD
        NSTRIP = MIN(NSTRIP, NBANDS)
!@TODO: the choice of NSTRIP needs to be evaluated and optimized

        LscaAWARE_LOCAL = LscaAWARE
__GPU   LscaAWARE_LOCAL = LscaAWARE_LOCAL .AND. (.NOT.OFFLOAD_ON)

        IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES1%COMM_KIN, WDES1%NB_TOTK(ISP))

        blocks: DO NB_RED = NBSTART_RED, NBSTOP_RED, NSTRIP
           ! the actual number of orbitals in this strip
           NSTRIP_RED = MIN(NSTRIP,NBSTOP_RED-NB_RED+1)

           ! contract
           IF (.NOT.LscaAWARE_LOCAL) THEN
GPACC !$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:COVL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
              COVL(:,NB_RED:NB_RED+NSTRIP_RED-1) = zero
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
              CALL ORTH1('U', &
                 WA1%CW_RED(1,1), WA2%CW_RED(1,NB_RED), &
                 WA1%CPROJ_RED(1,1), WA2%CPROJ_RED(1,NB_RED), &
                 NB_TOT, NB_RED, NSTRIP_RED, &
                 WDES1%NPL_RED, 0, WDES1%NRPLWV_RED, WDES1%NPROD_RED, &
                 COVL(1,1))
           ELSE
              CALL ORTH1_DISTRI('U', &
                 WA1%CW_RED(1,1), WA2%CW_RED(1,NB_RED), &
                 WA1%CPROJ_RED(1,1), WA2%CPROJ_RED(1,1), &
                 NB_TOT, NB_RED, NSTRIP_RED, &
                 WDES1%NPL_RED, 0, WDES1%NRPLWV_RED, WDES1%NPROD_RED, &
                 COVL(1,1), &
                 WDES1%COMM_KIN, WDES1%NB_TOTK(ISP))
           ENDIF
        ENDDO blocks

        ! @TODO: stripwise communication inside the blocks loop above
        ! could maybe be hidden behind computation ...
        IF (.NOT.LscaAWARE_LOCAL) THEN
           CALLMPI( M_sum_g(WDES1%COMM_KIN, COVL(1,NBSTART_RED), NBANDS*NB_TOT) )
        ENDIF

        PROFILING_STOP('braket')

        RETURN
      END SUBROUTINE BRAKET


!************************ SUBROUTINE HWESW *****************************
!
!> Solve the generalized eigenvalue problem  H*U = EV*S*U
!
!@TODO: add and handle exception INFO/=0 in LscaAWARE code-path
!@TODO: call PZHEEVD/PZHEEVX instead of PZHEEV
!
!***********************************************************************

      SUBROUTINE HWESW(COMM, H, S, EV, U, N, IU0)

__GPU   USE moffload_struct_def

        USE mpimy, ONLY : communic
        USE scala, ONLY : LscaAWARE, GSD, scalapack_des, INIT_scala_DESC, DLEN_
        USE tutor, ONLY : vtutor
        USE string, ONLY : str

        IMPLICIT NONE
        !> Communicator
        TYPE (communic) :: COMM
        !> Hamiltonian, overlap, and rotation matrices
        GDEF :: H(:,:), S(:,:), U(:,:)
        !> Eigenvalues
        REAL(q) :: EV(:)
        !> Order of H and S
        INTEGER :: N
        !> stdout
        INTEGER, OPTIONAL :: IU0

        ! local variables
        GDEF, ALLOCATABLE :: HTMP(:,:), STMP(:,:), UTMP(:,:)
        GDEF, ALLOCATABLE :: WORK(:)
        REAL(q), ALLOCATABLE :: RWORK(:)
        INTEGER :: NB, LWORK, INFO
        INTEGER, EXTERNAL :: ILAENV

        COMPLEX(q) :: WSIZE(1)
        REAL(q) :: RWSIZE(1)
        REAL(q) :: SCALE
        INTEGER :: LRWORK

#ifdef use_pdsyevd_pzheevd
        INTEGER, ALLOCATABLE :: IWORK(:)
        INTEGER :: IWSIZE(1)
        INTEGER :: LIWORK
#endif

        TYPE(scalapack_des) :: GS
        INTEGER :: DESC(DLEN_)
        LOGICAL :: LscaAWARE_LOCAL

        PROFILING_START('hwesw')

        LscaAWARE_LOCAL = LscaAWARE
__GPU   LscaAWARE_LOCAL = LscaAWARE_LOCAL .AND. (.NOT.OFFLOAD_ON)

        IF (.NOT.LscaAWARE_LOCAL) THEN
GPACC !$ACC ENTER DATA CREATE(EV) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:EV) IF(OFFLOAD_ON)
           ! only one rank solves the generalized eigenvalue problem
           IF (COMM%NODE_ME==COMM%IONODE) THEN
              ! we make temporary copies of both H and S, so that we can return without
              ! having destroyed these matrices in case the ZHEGV call below fails ...
              ALLOCATE(HTMP(N,N), STMP(N,N))
GPACC !$ACC ENTER DATA CREATE(HTMP,STMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:HTMP,STMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(H,S,HTMP,STMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:H,S,HTMP,STMP) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
              HTMP(1:N, 1:N) = H(1:N, 1:N)
              STMP(1:N, 1:N) = S(1:N, 1:N)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!@NOTE: to save memory one could use the lower triangular part of S to backup the upper triangle
!       (in that case one only needs additional memory to backup the diagonal of S).
#ifdef gammareal
              LWORK = MAX(1, 3*N-1)
              NB = ILAENV(1, 'DSYTRD', 'U', N, -1, -1, -1)
              LWORK = MAX(LWORK, (NB+2)*N)
              ALLOCATE(WORK(LWORK))
#ifdef CRAYHIP
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WORK) IF(OFFLOAD_ON)
#endif // CRAYHIP
              CALL WDSYGV(1, 'V', 'U', N, HTMP, N, STMP, N, EV, WORK, LWORK, INFO)
#ifdef CRAYHIP
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WORK) IF(OFFLOAD_ON)
#endif // CRAYHIP
              DEALLOCATE(WORK)
#else
              LWORK = MAX(1, 2*N-1)
              NB = ILAENV(1, 'ZHETRD', 'U', N, -1, -1, -1)
              LWORK = MAX(LWORK, (NB+1)*N)
              ALLOCATE(WORK(LWORK), RWORK(MAX(1, 3*N-2)))
#ifdef CRAYHIP
! only hip zhegv needs RWORK on the device
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:RWORK) IF(OFFLOAD_ON)
#endif // CRAYHIP

              CALL WZHEGV(1, 'V', 'U', N, HTMP, N, STMP, N, EV, WORK, LWORK, RWORK, INFO)
#ifdef CRAYHIP
! only hip zhegv needs RWORK on the device
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:RWORK) IF(OFFLOAD_ON)
#endif // CRAYHIP
              DEALLOCATE(WORK, RWORK)

#endif
GPACC !$ACC KERNELS PRESENT(U,HTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:U,HTMP) IF(OFFLOAD_ON)
              IF (INFO==0) U(1:N, 1:N) = HTMP(1:N, 1:N)
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC EXIT DATA DELETE(HTMP,STMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:HTMP,STMP) IF(OFFLOAD_ON)
              DEALLOCATE(HTMP, STMP)
           ENDIF

           ! broadcast INFO to all ranks
           CALLMPI(M_bcast_i(COMM, INFO, 1))
           IF (INFO/=0) THEN
              ! something went wrong ...
#ifdef gammareal
              CALL vtutor%error('HWESW: DSYGV exited with INFO = '//str(INFO)//' for N = '//str(N))
#else
              CALL vtutor%error('HWESW: ZHEGV exited with INFO = '//str(INFO)//' for N = '//str(N))
#endif
!@TODO: implement strategy to recover from this exception
              PROFILING_STOP('hwesw')
              RETURN
           ELSE
              ! the master rank has successfully solved the generalized eigenvalue problem
              ! now broadcast the result to all ranks
              CALLMPI( M_bcast_g(COMM, U(1,1), SIZE(U,1)*N))
              CALLMPI( M_bcast_d(COMM, EV(1), N))
           ENDIF
GPACC !$ACC EXIT DATA COPYOUT(EV) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:EV) IF(OFFLOAD_ON)
        ELSE

           CALL INIT_SCALA_DESC(COMM, N, DESC, GS)

           ALLOCATE(HTMP(GS%NP,GS%NQ), STMP(GS%NP,GS%NQ), UTMP(GS%NP,GS%NQ))

           CALL PGLACPY('U', N, N, H(1,1), 1, 1, GSD%DESCA, HTMP(1,1), 1, 1, GS%DESCA)
           CALL PGLACPY('U', N, N, S(1,1), 1, 1, GSD%DESCA, STMP(1,1), 1, 1, GS%DESCA)
#ifdef debug
           IF (PRESENT(IU0)) &
              CALL DUMP_HAM_BLOCK("Hamilton matrix hwesw", HTMP, .TRUE., GS%DESCA, COMM, 1, N, 1, N, IU0)
#endif

           CALL PGPOTRF('U', N, STMP(1,1), 1, 1, GS%DESCA, INFO)

#ifdef gammareal
           CALL PDSYGST(1, 'U', N, HTMP(1,1), 1, 1, GS%DESCA, STMP(1,1), 1, 1, GS%DESCA, SCALE, INFO)
#ifdef use_pdsyevd_pzheevd
           CALL PDSYEVD('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       RWSIZE, -1, IWSIZE, -1, INFO)

           LWORK = INT(RWSIZE(1))
           LIWORK = IWSIZE(1)

           ALLOCATE(WORK(LWORK), IWORK(LIWORK))

           CALL PDSYEVD('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       WORK, LWORK, IWORK, LIWORK, INFO)

           DEALLOCATE(WORK, IWORK)
#else // use_pdsyevd_pzheevd
           CALL PDSYEV('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       RWSIZE, -1, INFO)

           LWORK = INT(RWSIZE(1))

           ALLOCATE(WORK(LWORK))

           CALL PDSYEV('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       WORK, LWORK, INFO)

           DEALLOCATE(WORK)
#endif // use_pdsyevd_pzheevd
#else // gammareal
           CALL PZHEGST(1, 'U', N, HTMP(1,1), 1, 1, GS%DESCA, STMP(1,1), 1, 1, GS%DESCA, SCALE, INFO)
#ifdef use_pdsyevd_pzheevd
           CALL PZHEEVD('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       WSIZE, -1, RWSIZE, -1, IWSIZE, -1, INFO)

           LWORK = INT(WSIZE(1))
           LRWORK = INT(RWSIZE(1))
           LIWORK = IWSIZE(1)

           ALLOCATE(WORK(LWORK), RWORK(LRWORK), IWORK(LIWORK))

           CALL PZHEEVD('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

           DEALLOCATE(WORK, RWORK, IWORK)
#else // use_pdsyevd_pzheevd
           CALL PZHEEV('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       WSIZE, -1, RWSIZE, -1, INFO)

           LWORK = INT(WSIZE(1))
           LRWORK = INT(RWSIZE(1))

           ALLOCATE(WORK(LWORK), RWORK(LRWORK))

           CALL PZHEEV('V', 'U', N, HTMP(1,1), 1, 1, GS%DESCA, EV, UTMP(1,1), 1, 1, GS%DESCA, &
                       WORK, LWORK, RWORK, LRWORK, INFO)

           DEALLOCATE(WORK, RWORK)
#endif // use_pdsyevd_pzheevd
#endif // gammareal

           CALL PGTRSM('L', 'U', 'N', 'N', N, N, one, STMP(1,1), 1, 1, GS%DESCA, UTMP(1,1), 1, 1, GS%DESCA)

           CALL PGLACPY('F', N, N, UTMP(1,1), 1, 1, GS%DESCA, U(1,1), 1, 1, GSD%DESCA)
#ifdef debug
           IF (PRESENT(IU0)) &
              CALL DUMP_HAM_BLOCK("eigenvector hwesw", UTMP, .TRUE., GS%DESCA, COMM, 1, N, 1, N, IU0)
#endif
           DEALLOCATE(HTMP, STMP, UTMP)
        ENDIF

        PROFILING_STOP('hwesw')

        RETURN
      END SUBROUTINE HWESW


!************************ SUBROUTINE EXPAND ****************************
!
!***********************************************************************

      SUBROUTINE EXPAND(WDES1, W, HW, SW, EV, U, NB, NADD, IDO)

__GPU   USE moffload_struct_def

        USE wave_struct_def
        USE scala, ONLY : LscaAWARE, GSD, RECON_SLICE_REORDER
        USE dfast, ONLY : NSTRIPD
        USE tutor, ONLY : vtutor

        USE iso_c_binding, ONLY : c_loc, c_f_pointer

        IMPLICIT NONE

        !> Orbitals and the action of the Hamiltonian and overlap operator
        !> on these orbitals
        TYPE(wavefuna) :: W, HW, SW
        !>  Common descriptor for H, HW, and SW
        TYPE(wavedes1) :: WDES1
        !> Eigenvalues
        REAL(q) :: EV(:)
        !> Roration matrix
        GDEF, TARGET :: U(:,:)
        !> Size of the reduced basis
        INTEGER :: NB
        !> Number of residuals to add to the basis
        INTEGER :: NADD
        !> Indices (global) of the orbitals from which the residual will be computed
        INTEGER :: IDO(:)

        ! local variables
        GDEF, POINTER :: UGLB(:,:), UPTR(:)
        INTEGER :: NSTRIP, NPOS, NSTRIP_ACT
        INTEGER :: I, J

        LOGICAL :: LscaAWARE_LOCAL

        PROFILING_START('expand')

        !quick return if possible
        IF (NADD==0) THEN
           PROFILING_STOP('expand')
           RETURN
        ENDIF

        LscaAWARE_LOCAL = LscaAWARE
__GPU   LscaAWARE_LOCAL = LscaAWARE_LOCAL .AND. (.NOT.OFFLOAD_ON)

GPACC !$ACC ENTER DATA COPYIN(EV,IDO) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:EV,IDO) IF(OFFLOAD_ON)

        ! compact the eigenvector matrix according to LDO : we are
        ! only interested in the orbitals that need further work
!@TODO: implement and out-of-place reordering in the OpenACC case:
!       this will allow additional parallelisation over I
        DO I = 1, NADD
           IF (IDO(I)/=I) THEN
              IF (.NOT.LscaAWARE_LOCAL) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(U,IDO) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:U,IDO) IF(OFFLOAD_ON)
                 DO J = 1, NB
                    U(J,I) = U(J,IDO(I))
                 ENDDO
              ENDIF
              ! reorder the eigenvalue as well, for use below and in the preconditioner lateron
GPACC !$ACC KERNELS PRESENT(EV,IDO) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:EV,IDO) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
              EV(I) = EV(IDO(I))
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
           ENDIF
        ENDDO

#ifdef debug
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(U) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(U) IF(OFFLOAD_ON)
        IF (WDES1%COMM_KIN%NODE_ME == 1) &
              CALL DUMP_HAM_BLOCK("eigenvectors expand", U, LscaAWARE_LOCAL, GSD%DESCA, WDES1%COMM_KIN, 1, NB, 1, NADD, 6)
#endif

        IF (LscaAWARE_LOCAL) THEN
           NSTRIP = NSTRIPD
           ALLOCATE(UGLB(WDES1%NB_TOT,NSTRIP))
           ! we need a 1d pointer to U for RECON_SLICE_REORDER (below)
           CALL c_f_pointer(c_loc(U), UPTR, [SIZE(U, KIND=qi8)])
        ELSE
           NSTRIP = NADD
        ENDIF

        strip: DO NPOS = 1, NADD, NSTRIP
           NSTRIP_ACT = MIN(NSTRIP, NADD-NPOS+1)

           IF (LscaAWARE_LOCAL) THEN
              ! collect the distributed columns with global indices IDO(NPOS:NPOS+NSTRIP_ACT-1) into UGLB
              CALL RECON_SLICE_REORDER(UGLB, WDES1%NB_TOT, WDES1%NB_TOT, UPTR, GSD%DESCA, IDO(NPOS:NPOS+NSTRIP_ACT-1), NSTRIP_ACT)
              CALLMPI( M_sum_g(WDES1%COMM_KIN, UGLB(1,1), WDES1%NB_TOT*NSTRIP_ACT))
           ELSE
              UGLB => U(:,NPOS:NPOS+NSTRIP_ACT-1)
           ENDIF

           PROFILING_START('expand|gemm')

           ! | W_NB+j > = \sum_i U_ij H| W_i >, with i = [1 .. NB] and j = [1 .. NADD]
           CALL WGGEMM('N', 'N', m_ WDES1%NRPLWV_RED, NSTRIP_ACT, NB, &
                      one, HW%CW_RED(1,1), m_ WDES1%NRPLWV_RED, UGLB(1,1), SIZE(UGLB,1), &
                      zero, W%CW_RED(1,NB+NPOS), m_ WDES1%NRPLWV_RED)

           PROFILING_STOP('expand|gemm')

           ! | W_NB+j > = \| W_NB+j > - sum_i e_j U_ij S| W_i >, with i = [1 .. NB] and j = [1 .. NADD]
           ! UGLB points to U aka EVEC and hence is already on the device: no mapping needed
           ! this only works since UGLB is a native pointer
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(UGLB,EV) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:UGLB,EV) IF(OFFLOAD_ON)
           DO J = 1, NSTRIP_ACT
              DO I = 1, NB
                 UGLB(I,J) = -EV(NPOS+J-1)* UGLB(I,J)
              ENDDO
           ENDDO

           PROFILING_START('expand|gemm')

           CALL WGGEMM('N', 'N', m_ WDES1%NRPLWV_RED, NSTRIP_ACT, NB, &
                      one, SW%CW_RED(1,1), m_ WDES1%NRPLWV_RED, UGLB(1,1), SIZE(UGLB,1), &
                      one, W%CW_RED(1,NB+NPOS), m_ WDES1%NRPLWV_RED)

           PROFILING_STOP('expand|gemm')
        ENDDO strip

        IF (LscaAWARE_LOCAL) THEN
           DEALLOCATE(UGLB)
           NULLIFY(UPTR)
        ELSE
           NULLIFY(UGLB)
        ENDIF

        ! NADD vectors have now been added to the basis
        NB = NB + NADD

GPACC !$ACC EXIT DATA COPYOUT(EV) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(IDO) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:EV) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:IDO) IF(OFFLOAD_ON)

        PROFILING_STOP('expand')

        RETURN
      END SUBROUTINE EXPAND


!************************ SUBROUTINE SET_DATAKE_RED ********************
!
!***********************************************************************

      SUBROUTINE SET_DATAKE_RED(WDES1, DATAKE_RED)

__GPU   USE moffload_struct_def

        USE wave_struct_def
        USE wave_high, ONLY : NEWWAV, DELWAV, REDISTRIBUTE_PW

        IMPLICIT NONE

        !> Descriptor for wave functions at a particular k-point
        TYPE (wavedes1) :: WDES1
        !> Array of kinetic energy components (in over-plane-wave distribution)
        REAL(q) :: DATAKE_RED(:)

        ! local variables
        TYPE (wavefun1) :: W1
        INTEGER :: ISPINOR, M, MM

        PROFILING_START('set_datake_red')

GPACC !$ACC ENTER DATA CREATE(W1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W1) IF(OFFLOAD_ON)
        CALL NEWWAV(W1, WDES1, .FALSE.)

GPACC !$ACC KERNELS PRESENT(W1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W1) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        W1%CW = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES1,W1) PRIVATE(MM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:WDES1,W1) IF(OFFLOAD_ON)
        DO ISPINOR = 0, WDES1%NRSPINORS-1
           DO M = 1, WDES1%NGVECTOR
              MM = M + ISPINOR*WDES1%NGVECTOR
              W1%CW(MM) = WDES1%DATAKE(M,ISPINOR+1)
           ENDDO
        ENDDO

        CALL REDISTRIBUTE_PW(W1)

GPACC !$ACC KERNELS PRESENT(DATAKE_RED) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:DATAKE_RED) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        DATAKE_RED = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
        MM = (WDES1%COMM_INTER%NODE_ME-1) * WDES1%NRPLWV_RED

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,DATAKE_RED) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:WDES1,DATAKE_RED) IF(OFFLOAD_ON)
        DO M = 1, WDES1%NRPLWV_RED
           DATAKE_RED(M) = REAL(W1%CW(M+MM), q)
        ENDDO

        CALL DELWAV(W1, .FALSE.)
GPACC !$ACC EXIT DATA DELETE(W1) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1) IF(OFFLOAD_ON)

        PROFILING_STOP('set_datake_red')

        RETURN
      END SUBROUTINE SET_DATAKE_RED

#ifdef non_local_contributions_in_preconditioner
!************************ SUBROUTINE SET_NL_RED ************************
!
!***********************************************************************

      SUBROUTINE SET_NL_RED(WDES1, NONL_S, CDIJ, CQIJ, DIIG_RED, QIIG_RED)

        USE wave_struct_def
        USE nonl_struct_def
        USE wave_high, ONLY : NEWWAV, DELWAV, REDISTRIBUTE_PW

        IMPLICIT NONE

        !> Descriptor for wave functions at a particular k-point
        TYPE (wavedes1) :: WDES1
        !> PAW projector (reciprocal space)
        TYPE (nonl_struct) :: NONL_S
        !> PAW strength parameters and augmentation charges
        OVERLAP :: CDIJ(:,:,:,:), CQIJ(:,:,:,:)
        !> D(G) (in over-plane-wave distribution)
        REAL(q) :: DIIG_RED(:)
        !> Q(G) (in over-plane-wave distribution)
        REAL(q) :: QIIG_RED(:)

        ! local variables
        TYPE (wavefun1) :: W1
        REAL(q), ALLOCATABLE :: DSUM(:,:,:,:)
        REAL(q) :: RTMP
        INTEGER :: ISPINOR, NI, NT, NT_GLOBAL, LMMAXC, L, LP, M, MM

        PROFILING_START('set_nl_red')

        CALL NEWWAV(W1, WDES1, .FALSE.)

        ALLOCATE(DSUM(SIZE(CDIJ,1), SIZE(CDIJ,2), NONL_S%NTYP, WDES1%NRSPINORS))

        DSUM = 0
        DO ISPINOR = 1, WDES1%NRSPINORS
           DO NI = 1, WDES1%NIONS
              NT = WDES1%ITYP(NI)
              LMMAXC = WDES1%LMMAX(NT)
              IF (LMMAXC==0) CYCLE
              NT_GLOBAL = WDES1%NT_GLOBAL(NT)
              DSUM(1:LMMAXC,1:LMMAXC,NT_GLOBAL,ISPINOR) = DSUM(1:LMMAXC,1:LMMAXC,NT_GLOBAL,ISPINOR) + &
                                                           CDIJ(1:LMMAXC,1:LMMAXC,NI,ISPINOR*ISPINOR)
           ENDDO
        ENDDO

        CALLMPI( M_sum_d(WDES1%COMM_INB, DSUM(1,1,1,1), SIZE(DSUM)))

        W1%CW = 0
        DO ISPINOR = 1, WDES1%NRSPINORS
           DO NT = 1, NONL_S%NTYP
              LMMAXC = WDES1%LMMAX(NT)
              IF (LMMAXC==0) CYCLE
              DO L = 1, LMMAXC
                 DO LP = 1, LMMAXC
                    DO M = 1, WDES1%NGVECTOR
                       MM = M + (ISPINOR-1)*WDES1%NGVECTOR
                       RTMP = NONL_S%QPROJ(M,L,NT,WDES1%NK,ISPINOR)*NONL_S%QPROJ(M,LP,NT,WDES1%NK,ISPINOR)
                       W1%CW(MM) = W1%CW(MM) + DSUM(L,LP,NT,ISPINOR)*RTMP
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        CALL REDISTRIBUTE_PW(W1)

        DIIG_RED = 0
        MM = (WDES1%COMM_INTER%NODE_ME-1) * WDES1%NRPLWV_RED
        DO M = 1, WDES1%NRPLWV_RED
           DIIG_RED(M) = REAL(W1%CW(M+MM), q)
        ENDDO


        DSUM = 0
        DO ISPINOR = 1, WDES1%NRSPINORS
           DO NI = 1, WDES1%NIONS
              NT = WDES1%ITYP(NI)
              LMMAXC = WDES1%LMMAX(NT)
              IF (LMMAXC==0) CYCLE
              NT_GLOBAL = WDES1%NT_GLOBAL(NT)
              DSUM(1:LMMAXC,1:LMMAXC,NT_GLOBAL,ISPINOR) = DSUM(1:LMMAXC,1:LMMAXC,NT_GLOBAL,ISPINOR) + &
                                                           CQIJ(1:LMMAXC,1:LMMAXC,NI,ISPINOR*ISPINOR)

           ENDDO
        ENDDO

        CALLMPI( M_sum_d(WDES1%COMM_INB, DSUM(1,1,1,1), SIZE(DSUM)))

        W1%CW = 0
        DO ISPINOR = 1, WDES1%NRSPINORS
           DO NT = 1, NONL_S%NTYP
              LMMAXC = WDES1%LMMAX(NT)
              IF (LMMAXC==0) CYCLE
              DO L = 1, LMMAXC
                 DO LP = 1, LMMAXC
                    DO M = 1, WDES1%NGVECTOR
                       MM = M + (ISPINOR-1)*WDES1%NGVECTOR
                       RTMP = NONL_S%QPROJ(M,L,NT,WDES1%NK,ISPINOR)*NONL_S%QPROJ(M,LP,NT,WDES1%NK,ISPINOR)
                       W1%CW(MM) = W1%CW(MM) + DSUM(L,LP,NT,ISPINOR)*RTMP
                    ENDDO
                 ENDDO
              ENDDO
          ENDDO
        ENDDO

        CALL REDISTRIBUTE_PW(W1)

        QIIG_RED = 0
        MM = (WDES1%COMM_INTER%NODE_ME-1) * WDES1%NRPLWV_RED
        DO M = 1, WDES1%NRPLWV_RED
           QIIG_RED(M) = REAL(W1%CW(M+MM), q)
        ENDDO

        DEALLOCATE(DSUM)

        CALL DELWAV(W1, .FALSE.)

        PROFILING_STOP('set_nl_red')

        RETURN
      END SUBROUTINE SET_NL_RED
#endif

!************************ SUBROUTINE PRECONDITION **********************
!
!***********************************************************************
#ifdef non_local_contributions_in_preconditioner
      SUBROUTINE PRECONDITION(IPREC, DATAKE_RED, DIIG_RED, QIIG_RED, SLOCAL, DEVMAX, EV, WDES1, W, NBSTART, NBSTOP, FNORM)
#else
      SUBROUTINE PRECONDITION(IPREC, DATAKE_RED, SLOCAL, DEVMAX, EV, WDES1, W, NBSTART, NBSTOP, FNORM)
#endif
__GPU   USE moffload_struct_def

        USE wave_struct_def

        IMPLICIT NONE

        !> Type of preconditioner
        INTEGER :: IPREC
        !> Kinetic energy components (in accordance with the over-plane-wave distribution)
        REAL(q) :: DATAKE_RED(:)
#ifdef non_local_contributions_in_preconditioner
        !> Diagonal of the non-local PAW strength and augmentation (in over-plane-wave distribution)
        REAL(q) :: DIIG_RED(:), QIIG_RED(:)
#endif
        !> Average local potential
        REAL(q) :: SLOCAL
        !> Spread in the eigenvalues
        REAL(q) :: DEVMAX
        !> Eigenvalues of the orbitals for which additional basis functions are constructed
        REAL(q) :: EV(:)
        !> Descriptor for residuals at the current k-point
        TYPE (wavedes1) :: WDES1
        !> residuals to be preconditioned
        TYPE (wavefuna) :: W
       !> STart and stop index (global)
        INTEGER :: NBSTART, NBSTOP
        !>  Pass back the norm of the residuals
        REAL(q) :: FNORM(:)

        ! local variables
        REAL(q), ALLOCATABLE :: EKIN(:), FNRM(:), FPRE(:)
        REAL(q) :: FAKT, X, X2, SCALE
        INTEGER NB, I, J, IB, NPL_RED

        REAL(q) :: EKIN_, FNRM_, FPRE_

        COMPLEX(q), EXTERNAL :: ZDOTC

        PROFILING_START('precondition')

        ! quick return if possible
        IF (NBSTART>NBSTOP) THEN
           PROFILING_STOP('precondition')
           RETURN
        ENDIF

        NB = NBSTOP-NBSTART+1
        NPL_RED = WDES1%NPL_RED

        ALLOCATE(FNRM(NB), FPRE(NB))
GPACC !$ACC ENTER DATA CREATE(FNRM,FPRE) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:FNRM,FPRE) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(FNRM,FPRE) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNRM,FPRE) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        FNRM = 0 ; FPRE = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

!!        ! Compute the norm of the residuals before preconditioning
!!        RES = 0
!!        DO I = 1, NB
!!           IB = I+NBSTART-1
!!           RES(I) = ZDOTC(WDES1%NPL_RED, W%CW_RED(1,IB), 1, W%CW_RED(1,IB), 1)
!!        ENDDO
!!        CALLMPI( M_sum_d(WDES1%COMM_KIN, RES(1), NB))
!
        ! Precondition W
        IF (IPREC==0 .OR. IPREC==6 .OR. IPREC==8) THEN
           ALLOCATE(EKIN(NB))
GPACC !$ACC ENTER DATA CREATE(EKIN) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:EKIN) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,EKIN,DATAKE_RED,W) PRIVATE(IB,EKIN_,J) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WDES1,EKIN,DATAKE_RED,W) &
GPOMP !$OMP PRIVATE(IB,EKIN_,J) IF(OFFLOAD_ON)
           DO I = 1, NB
              IB = I+NBSTART-1
              EKIN_ = 0
GPACC !$ACC LOOP REDUCTION(+:EKIN_)
GPOMP !$OMP PARALLEL DO REDUCTION(+:EKIN_)
              DO J = 1, NPL_RED
                 EKIN_ = EKIN_ + DATAKE_RED(J) * W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
              ENDDO
              EKIN(I) = EKIN_
           ENDDO
           CALLMPI( M_sum_d(WDES1%COMM_KIN, EKIN(1), NB))

GPACC !$ACC PARALLEL LOOP PRESENT(EKIN,WDES1,W,FNRM,FPRE) PRIVATE(FAKT,IB,X,X2,FNRM_,FPRE_,J) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:EKIN,WDES1,W,FNRM,FPRE) &
GPOMP !$OMP PRIVATE(FAKT,IB,X,X2,FNRM_,FPRE_,J) IF(OFFLOAD_ON)
           DO I = 1, NB
              IF (EKIN(I)<2._q) EKIN(I) = 2._q
              EKIN(I) = EKIN(I)*1.5_q
              FAKT = 2._q/EKIN(I)
              IB = I+NBSTART-1
              FNRM_ = 0; FPRE_ = 0
GPACC !$ACC LOOP REDUCTION(+:FNRM_,FPRE_)
GPOMP !$OMP PARALLEL DO REDUCTION(+:FNRM_,FPRE_)
              DO J = 1, NPL_RED
                 ! compute norm of residual
                 FNRM_ = FNRM_ + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
                 ! precondition  residual
                 X = DATAKE_RED(J)/EKIN(I)
                 X2 = 27+X*(18+X*(12+8*X))
                 W%CW_RED(J,IB) = W%CW_RED(J,IB) * X2/(X2+16*X*X*X*X)*FAKT
                 ! and compute norm of preconditioned residual
                 FPRE_ = FPRE_ + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
              ENDDO
              FNRM(I) = FNRM_
              FPRE(I) = FPRE_
           ENDDO

GPACC !$ACC EXIT DATA DELETE(EKIN) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:EKIN) IF(OFFLOAD_ON)
           DEALLOCATE(EKIN)
        ELSEIF (IPREC==9) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(EV,SLOCAL,WDES1,W,DATAKE_RED,FNRM,FPRE) PRIVATE(X2,IB,X,FNRM_,FPRE_,J) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:EV,SLOCAL,WDES1,W,DATAKE_RED,FNRM,FPRE) &
GPOMP !$OMP PRIVATE(X2,IB,X,FNRM_,FPRE_,J) IF(OFFLOAD_ON)
           DO I = 1, NB
              ! the first NB entries of EV should contain the eigenvalues of the
              ! orbitals from whch the residuals W%CW_RED(:,NBSTART:NBSTOP) were
              ! calculated
              X2 = EV(I) - SLOCAL
              IB = I+NBSTART-1
              FNRM_ = 0 ; FPRE_ = 0
GPACC !$ACC LOOP REDUCTION(+:FNRM_,FPRE_)
GPOMP !$OMP PARALLEL DO REDUCTION(+:FNRM_,FPRE_)
              DO J = 1, NPL_RED
                 ! compute norm of residual
                 FNRM_ = FNRM_ + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
                 ! precondition residual
                 X = MAX(DATAKE_RED(J)-X2, 0._q)
                 W%CW_RED(J,IB) = W%CW_RED(J,IB) * REAL(1._q/(X+CMPLX(0._q,DEVMAX)), KIND=q)
                 ! and compute norm of preconditioned residual
                 FPRE_ = FPRE_ + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
              ENDDO
              FNRM(I) = FNRM_
              FPRE(I) = FPRE_
           ENDDO
#ifdef non_local_contributions_in_preconditioner
        ELSEIF (IPREC==1) THEN
          DO I = 1, NB
              ! the first NB entries of EV should contain the eigenvalues of the
              ! orbitals from whch the residuals W%CW_RED(:,NBSTART:NBSTOP) were
              ! calculated
              IB = I+NBSTART-1
              DO J = 1, NPL_RED
                 ! compute norm of residual
                 FNRM(I) = FNRM(I) + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
                 ! precondition residual
                 X = DATAKE_RED(J)+SLOCAL+DIIG_RED(J)-EV(I)*(1._q+QIIG_RED(J))
                 IF (ABS(X)<1E-4_q) THEN
                    X = SIGN(1E-4_q, X)
                 ENDIF
                 W%CW_RED(J,IB) = W%CW_RED(J,IB) / X
                 ! and compute norm of preconditioned residual
                 FPRE(I) = FPRE(I) + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
              ENDDO
           ENDDO
#endif
        ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(WDES1,W,FNRM) PRIVATE(IB,FNRM_,J) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WDES1,W,FNRM) &
GPOMP !$OMP PRIVATE(IB,FNRM_,J) IF(OFFLOAD_ON)
           DO I = 1, NB
              IB = I+NBSTART-1
              FNRM_ = 0
GPACC !$ACC LOOP REDUCTION(+:FNRM_)
GPOMP !$OMP PARALLEL DO REDUCTION(+:FNRM_)
              DO J = 1, NPL_RED
                 ! compute norm of residual
                 FNRM_ = FNRM_ + W%CW_RED(J,IB) * CONJG(W%CW_RED(J,IB))
              ENDDO
              FNRM(I) = FNRM_
           ENDDO
GPACC !$ACC KERNELS PRESENT(FNRM,FPRE) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNRM,FPRE) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
           FPRE = FNRM
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
        ENDIF      

        CALLMPI( M_sum_d(WDES1%COMM_KIN, FNRM(1), NB))
        CALLMPI( M_sum_d(WDES1%COMM_KIN, FPRE(1), NB))

GPACC !$ACC EXIT DATA COPYOUT(FNRM,FPRE) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:FNRM,FPRE) IF(OFFLOAD_ON)

        ! and normalize W
        DO I = 1, NB
           IF (FPRE(I)<=0) THEN
              SCALE = -1._q/SQRT(-FPRE(I))
           ELSE
              SCALE =  1._q/SQRT( FPRE(I))
           ENDIF
           IB = I+NBSTART-1
           CALL WZDSCAL(NPL_RED, SCALE, W%CW_RED(1,IB), 1)
        ENDDO
        
        !> Pass back the norm of residuals
        FNORM(1:NB) = SQRT(ABS(FNRM(1:NB)))

        DEALLOCATE(FNRM, FPRE)

        PROFILING_STOP('precondition')

        RETURN
      END SUBROUTINE PRECONDITION


!************************ SUBROUTINE EXTRACT ***************************
!
!***********************************************************************

      SUBROUTINE EXTRACT(WDES1, W, U, NB_TOT, NB, LPROJ)

__GPU   USE moffload_struct_def
        USE wave_struct_def

        USE dfast, ONLY : NBLK, NSTRIPD
        USE scala, ONLY : LscaAWARE, GSD, RECON_SLICE
        USE tutor, ONLY : vtutor

        USE iso_c_binding, ONLY : c_loc, c_f_pointer

        IMPLICIT NONE

        !> wave function descriptor
        TYPE (wavedes1) :: WDES1
        !> orbitals
        TYPE (wavefuna) :: W
        !> rotation matrix
        GDEF, TARGET :: U(:,:)
        !> number of functions in the reduced basis
        INTEGER :: NB_TOT
        !> number of orbitals ro extract
        INTEGER :: NB
        !> extract projections as well?
        LOGICAL, OPTIONAL :: LPROJ

        ! local variables
        COMPLEX(q), ALLOCATABLE, TARGET :: CBLOCK(:,:)
        GDEF, POINTER :: GBLOCK(:,:)
        INTEGER :: NPL_RED, NPRO_RED, NRPLWV_RED, NPROD_RED, NPOS, NBLK_ACT, N, I
        LOGICAL :: LPRJ

        COMPLEX(q), ALLOCATABLE :: CW_RED(:,:)
        GDEF, ALLOCATABLE :: CPROJ_RED(:,:)
        GDEF, ALLOCATABLE :: UGLB(:,:)
        GDEF, POINTER :: UPTR(:)
        INTEGER :: NPOSB, NSTRIP, NSTRIP_ACT

        LOGICAL :: LscaAWARE_LOCAL

        PROFILING_START('extract')

        LscaAWARE_LOCAL = LscaAWARE
__GPU   LscaAWARE_LOCAL = LscaAWARE_LOCAL .AND. (.NOT.OFFLOAD_ON)

        NPL_RED = WDES1%NPL_RED
        NRPLWV_RED = WDES1%NRPLWV_RED
        NPRO_RED = WDES1%NPRO_RED
        NPROD_RED = WDES1%NPROD_RED

        ! Do we want to extract projections as well?
        LPRJ = .FALSE. ; IF (PRESENT(LPROJ)) LPRJ = LPROJ

        IF (.NOT. LscaAWARE_LOCAL) THEN

           ALLOCATE(CBLOCK(NBLK,NB_TOT))
GPACC !$ACC ENTER DATA CREATE(CBLOCK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CBLOCK) IF(OFFLOAD_ON)
           ! get a GDEF pointer to CBLOCK ... this will only differ for the gamma-only case
           CALL c_f_pointer(c_loc(CBLOCK), GBLOCK, [NBLK, m_ NB_TOT])

           DO NPOS = 1, NPL_RED, NBLK

              NBLK_ACT = MIN(NBLK, NPL_RED-NPOS+1)

              PROFILING_START('extract|copy')

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CBLOCK,W%CW_RED) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CBLOCK,W%CW_RED) IF(OFFLOAD_ON)
              DO N = 1, NB_TOT
                 DO I = 1, NBLK_ACT
                    CBLOCK(I,N) = W%CW_RED(NPOS+I-1,N)
                 ENDDO
              ENDDO

              PROFILING_STOP('extract|copy')

              PROFILING_START('extract|gemm')

              CALL WGGEMM('N', 'N', m_ NBLK_ACT, NB, NB_TOT, &
                         one, CBLOCK(1,1), m_ NBLK, U(1,1), SIZE(U,1), &
                         zero, W%CW_RED(NPOS,1), m_ NRPLWV_RED )

              PROFILING_STOP('extract|gemm')
           ENDDO

           IF (LPRJ) THEN
GPACC !$ACC ENTER DATA CREATE(GBLOCK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:GBLOCK) IF(OFFLOAD_ON)
              DO NPOS = 1, NPRO_RED, NBLK

                 NBLK_ACT = MIN(NBLK, NPRO_RED-NPOS+1)

                 PROFILING_START('extract|copy')

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GBLOCK,W%CPROJ_RED) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GBLOCK,W%CPROJ_RED) IF(OFFLOAD_ON)
                 DO N = 1, NB_TOT
                    DO I = 1, NBLK_ACT
                       GBLOCK(I,N) = W%CPROJ_RED(NPOS+I-1,N)
                    ENDDO
                 ENDDO

                 PROFILING_STOP('extract|copy')

                 PROFILING_START('extract|gemm')

                 CALL WGGEMM('N', 'N', NBLK_ACT, NB, NB_TOT, &
                            one, GBLOCK(1,1), NBLK, U(1,1), SIZE(U,1), &
                            zero, W%CPROJ_RED(NPOS,1), NPROD_RED )

                 PROFILING_STOP('extract|gemm')
              ENDDO

GPACC !$ACC EXIT DATA DELETE(GBLOCK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:GBLOCK) IF(OFFLOAD_ON)
           ENDIF

GPACC !$ACC EXIT DATA DELETE(CBLOCK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CBLOCK) IF(OFFLOAD_ON)
           NULLIFY(GBLOCK) ; DEALLOCATE(CBLOCK)

        ELSE
           ALLOCATE(CW_RED(NPL_RED,NB)) ; IF (LPRJ) ALLOCATE(CPROJ_RED(NPRO_RED,NB))

           NSTRIP = NSTRIPD

           ALLOCATE(UGLB(WDES1%NB_TOT,NSTRIP))
           ! we need a 1d pointer to U for RECON_SLICE (below)
           CALL c_f_pointer(c_loc(U), UPTR, [SIZE(U, KIND=qi8)])

           DO NPOSB = 1, NB, NSTRIP
              NSTRIP_ACT = MIN(NSTRIP,NB-NPOSB+1)

              CALL RECON_SLICE(UGLB, WDES1%NB_TOT, WDES1%NB_TOT, UPTR, GSD%DESCA, NPOSB, NPOSB+NSTRIP_ACT-1)
              CALLMPI( M_sum_g(WDES1%COMM_KIN, UGLB(1,1), WDES1%NB_TOT*NSTRIP_ACT))

              PROFILING_START('extract|gemm')

              IF (NPL_RED>0) &
                 CALL GGEMM('N', 'N', m_ NPL_RED, NSTRIP_ACT, NB_TOT, &
                            one, W%CW_RED(1,1), m_ NRPLWV_RED, UGLB(1,1), SIZE(UGLB,1), &
                            zero, CW_RED(1,NPOSB), m_ NPL_RED )

              IF (LPRJ .AND. NPRO_RED>0) &
                 CALL GGEMM('N', 'N', NPRO_RED, NSTRIP_ACT, NB_TOT, &
                            one, W%CPROJ_RED(1,1), NPROD_RED, UGLB(1,1), SIZE(UGLB,1), &
                            zero, CPROJ_RED(1,NPOSB), NPRO_RED )

              PROFILING_STOP('extract|gemm')
           ENDDO

           DEALLOCATE(UGLB) ; NULLIFY(UPTR)

           PROFILING_START('extract|copy')

           W%CW_RED(1:NPL_RED,1:NB) = CW_RED(1:NPL_RED,1:NB)
           DEALLOCATE(CW_RED)

           IF (LPRJ) THEN
              W%CPROJ_RED(1:NPRO_RED,1:NB) = CPROJ_RED(1:NPRO_RED,1:NB)
              DEALLOCATE(CPROJ_RED)
           ENDIF

           PROFILING_STOP('extract|copy')
        ENDIF

        PROFILING_STOP('extract')

        RETURN
      END SUBROUTINE EXTRACT


!************************ SUBROUTINE CREATE_BUF ************************
!
!***********************************************************************

      SUBROUTINE CREATE_BUF(WDES1, NSTRIP, LASYNC, BUF, NHANDLE)

        USE wave_struct_def
        USE wave_mpi, ONLY : REDIS_PW_ALLOC

        IMPLICIT NONE

        TYPE(wavedes1), TARGET :: WDES1
        TYPE(buffer) :: BUF
        INTEGER :: NSTRIP
        LOGICAL :: LASYNC

        INTEGER, OPTIONAL :: NHANDLE

        ! local variables
        INTEGER ::I, NH

        PROFILING_START('create_buf')

        BUF%LASYNC = LASYNC
        IF (LASYNC) THEN
           ! The asynchronous redistribution routines require a wavedes descriptor
           ! and we only have a wavedes1 descriptor at this point ... so we resort
           ! to this dirty hack
           BUF%WDES%NRPLWV = WDES1%NRPLWV
           BUF%WDES%NB_PAR = WDES1%NB_PAR
           BUF%WDES%NBANDS = WDES1%NBANDS
           BUF%WDES%COMM_INTER => WDES1%COMM_INTER

           NH = 1 ; IF (PRESENT(NHANDLE)) NH = NHANDLE
           ALLOCATE(BUF%HANDLE(NH))
           DO I = 1, NH
              CALL REDIS_PW_ALLOC(BUF%WDES, NSTRIP, BUF%HANDLE(I)%H)
              BUF%HANDLE(I)%ISTART = 0
              BUF%HANDLE(I)%NSTRIP = 0
           ENDDO
        ENDIF

        PROFILING_STOP('create_buf')

        RETURN
      END SUBROUTINE CREATE_BUF


!************************ SUBROUTINE DESTROY_BUF ***********************
!
!***********************************************************************

      SUBROUTINE DESTROY_BUF(BUF)

        USE wave_mpi, ONLY : REDIS_PW_DEALLOC

        IMPLICIT NONE

        TYPE(buffer) :: BUF

        ! local variables
        INTEGER :: I

        PROFILING_START('destroy_buf')

        IF (BUF%LASYNC) THEN
           IF (ASSOCIATED(BUF%WDES%COMM_INTER)) NULLIFY(BUF%WDES%COMM_INTER)
           DO I = 1, SIZE(BUF%HANDLE)
              CALL REDIS_PW_DEALLOC(BUF%HANDLE(I)%H)
           ENDDO
           DEALLOCATE(BUF%HANDLE)
        ENDIF

        PROFILING_STOP('destroy_buf')

        RETURN
      END SUBROUTINE DESTROY_BUF


!************************ SUBROUTINE REDIS *****************************
!
!> Redistribute NSTRIP orbitals from the default to over-plane-wave
!> distribution (or vice versa).
!> INDEX is the label of the first oribtal in the strip (normally its
!> band index) and is used for internal bookkeeping by the subroutines
!
!***********************************************************************

      SUBROUTINE REDIS(BUF, WA, ISTART, ISTOP, IHANDLE)
        USE wave_struct_def
        USE wave_high, ONLY : ELEMENTS, REDISTRIBUTE_PW

        IMPLICIT NONE

        TYPE(buffer) :: BUF
        TYPE(wavefuna) :: WA
        INTEGER :: ISTART, ISTOP
        INTEGER, OPTIONAL :: IHANDLE
        ! local variables
        INTEGER :: IH, NSTRIP

        PROFILING_START('redis')

        ! early exit if possible
        IF (.NOT.WA%WDES1%DO_REDIS) THEN
           PROFILING_STOP('redis')
           RETURN
        ENDIF

        IF (BUF%LASYNC) THEN
           IH = 1 ; IF (PRESENT(IHANDLE)) IH = IHANDLE
           NSTRIP = ISTOP-ISTART+1
           CALL REDIS_PW_STRIP_START(BUF%WDES, WA%CW(1,ISTART), NSTRIP, ISTART, BUF%HANDLE(IH)%H)
           BUF%HANDLE(IH)%ISTART = ISTART
           BUF%HANDLE(IH)%NSTRIP = NSTRIP
        ELSE
           CALL REDISTRIBUTE_PW(ELEMENTS(WA, ISTART, ISTOP))
        ENDIF

        PROFILING_STOP('redis')

        RETURN
      END SUBROUTINE REDIS


!************************ SUBROUTINE WAIT ******************************
!
!***********************************************************************

      SUBROUTINE WAIT(BUF, WA, IHANDLE)
        USE wave_struct_def

        IMPLICIT NONE

        TYPE(buffer) :: BUF
        TYPE(wavefuna) :: WA
        INTEGER, OPTIONAL :: IHANDLE
        ! local variables
        INTEGER :: IH, ISTART, NSTRIP

        PROFILING_START('wait')

        ! early exit if possible
        IF (.NOT.BUF%LASYNC .OR. .NOT.WA%WDES1%DO_REDIS) THEN
           PROFILING_STOP('wait')
           RETURN
        ENDIF

        IH = 1 ; IF (PRESENT(IHANDLE)) IH = IHANDLE

        NSTRIP = BUF%HANDLE(IH)%NSTRIP

        ! another possible early exit
        IF (NSTRIP == 0) THEN
           PROFILING_STOP('wait')
           RETURN
        ENDIF

        ISTART = BUF%HANDLE(IH)%ISTART

        CALL REDIS_PW_STRIP_STOP(BUF%WDES, WA%CW(1,ISTART), NSTRIP, ISTART, BUF%HANDLE(IH)%H)

        BUF%HANDLE(IH)%ISTART = 0
        BUF%HANDLE(IH)%NSTRIP = 0

        PROFILING_STOP('wait')

        RETURN
      END SUBROUTINE WAIT


!************************ SUBROUTINE SET_TO_IDENTITY *******************
!
!***********************************************************************

      SUBROUTINE SET_TO_IDENTITY(N, A, DESC)

        USE iso_c_binding, ONLY : c_loc, c_f_pointer
        USE scala, ONLY : SET_DIAGONALE_REAL

        IMPLICIT NONE

        INTEGER :: N
        GDEF, TARGET :: A(:,:)
        INTEGER :: DESC(:)

        ! local variables
        GDEF, POINTER :: GPTR(:)
        REAL(q), ALLOCATABLE :: DIAG(:)

        ! we need a 1d pointer to A for SET_DIAGONALE_REAL (below)
        CALL c_f_pointer(c_loc(A), GPTR, [SIZE(A, KIND=qi8)])

        ! the first N elements of the diagonal will be set to one
        ALLOCATE(DIAG(N)); DIAG = 1._q

        ! first set A to zero
        A = 0

        ! now put the values from DIAG onto the diagonal of A
        CALL SET_DIAGONALE_REAL(DESC(3), GPTR, DIAG, DESC)

        ! cleanup
        DEALLOCATE(DIAG)
        NULLIFY(GPTR)

      END SUBROUTINE SET_TO_IDENTITY


!************************ SUBROUTINE SET_DIAGONAL **********************
!
!***********************************************************************

      SUBROUTINE SET_DIAGONAL(DIAG, A, DESC)

        USE iso_c_binding, ONLY : c_loc, c_f_pointer
        USE scala, ONLY : SET_DIAGONALE_REAL

        IMPLICIT NONE

        REAL(q) :: DIAG(:)
        GDEF, TARGET :: A(:,:)
        INTEGER :: DESC(:)

        ! local variables
        GDEF, POINTER :: GPTR(:)

        ! we need a 1d pointer to A for SET_DIAGONALE_REAL (below)
        CALL c_f_pointer(c_loc(A), GPTR, [SIZE(A, KIND=qi8)])

        ! first set A to zero
        A = 0

        ! now put the values from DIAG onto the diagonal of A
        CALL SET_DIAGONALE_REAL(DESC(3), GPTR, DIAG, DESC)

        ! cleanup
        NULLIFY(GPTR)

      END SUBROUTINE SET_DIAGONAL


!************************ SUBROUTINE DUMP_HAM_BLOCK ********************
!
!***********************************************************************

      SUBROUTINE DUMP_HAM_BLOCK(STRING, CHAM, LDIST, DESC, COMM, ISTART, ISTOP, JSTART, JSTOP, IUNIT)

        USE mpimy, ONLY : communic
        USE scala, ONLY : RECON_SLICE
        USE iso_c_binding, ONLY : c_loc, c_f_pointer

        IMPLICIT NONE

        CHARACTER (LEN=*) :: STRING
        GDEF, TARGET :: CHAM(:,:)
        INTEGER :: ISTART, ISTOP, JSTART, JSTOP
        INTEGER :: IUNIT

        TYPE(communic) :: COMM
        INTEGER :: DESC(:)
        LOGICAL :: LDIST

        ! local variables
        GDEF, POINTER :: COUT(:,:), CPTR(:)
        INTEGER :: N1, N2

        IF (LDIST) THEN
           ALLOCATE(COUT(DESC(3),JSTOP-JSTART+1))
           CALL c_f_pointer(c_loc(CHAM), CPTR, [SIZE(CHAM, KIND=qi8)])
           CALL RECON_SLICE(COUT, DESC(3), DESC(3), CPTR, DESC, JSTART, JSTOP)
           CALLMPI( M_sum_g(COMM, COUT, SIZE(COUT)))
        ELSE
           COUT => CHAM
        ENDIF

        IF (IUNIT>=0) THEN
           WRITE(IUNIT,*) STRING
           WRITE(IUNIT,'(5X)',ADVANCE='NO')
           DO N2 = JSTART, JSTOP
             WRITE(IUNIT,'(i12)',ADVANCE='NO') N2
           ENDDO
           WRITE(IUNIT,*)

           DO N1 = ISTART, ISTOP
              WRITE(IUNIT,'(I3,2X)', ADVANCE='No') N1
              DO N2 = JSTART, JSTOP
                WRITE(IUNIT,1, ADVANCE='No') REAL(COUT(N1,N2), KIND=q)
              ENDDO
              WRITE(IUNIT,*)
           ENDDO
           WRITE(IUNIT,*)
#ifndef gammareal
           DO N1 = ISTART, ISTOP
              WRITE(IUNIT,'(I3,2X)', ADVANCE='No') N1
              DO N2 = JSTART, JSTOP
                WRITE(IUNIT,1, ADVANCE='No') AIMAG(COUT(N1,N2))
              ENDDO
              WRITE(IUNIT,*)
           ENDDO
           WRITE(IUNIT,*)
#endif
        ENDIF

        IF (LDIST) THEN
           DEALLOCATE(COUT)
           NULLIFY(CPTR)
        ELSE
           NULLIFY(COUT)
        ENDIF

!!   1    FORMAT(E12.4)
   1    FORMAT(F12.5)

      END SUBROUTINE DUMP_HAM_BLOCK

      END MODULE david_full

#else // MPI && ScaLAPACK

      MODULE david_full

      CONTAINS

      SUBROUTINE EDDAV_FULL(HAMILTONIAN, SV, CDIJ, CQIJ, NONLR_S, NONL_S, INFO, LATT_CUR, WDES, W, &
            DEX, ETHR, RMS, DESUM, EXHF, NRED, IU0, IU6, LDELAY, LRESIDUAL, LHF, NKSTART, NKSTOP)

        USE base

        USE wave_struct_def
        USE hamil_struct_def
        USE nonl_struct_def
        USE nonlr_struct_def
        USE poscar_struct_def

        USE tutor, ONLY : vtutor

        IMPLICIT NONE

        TYPE(ham_handle) :: HAMILTONIAN
        RGRID :: SV(:,:)
        OVERLAP :: CDIJ(:,:,:,:), CQIJ(:,:,:,:)
        TYPE (nonlr_struct) :: NONLR_S
        TYPE (nonl_struct) :: NONL_S
        TYPE (info_struct) :: INFO
        TYPE (latt) :: LATT_CUR
        TYPE(wavedes) :: WDES
        TYPE(wavespin) :: W
        REAL(q) :: DEX
        REAL(q) :: ETHR
        REAL(q) :: RMS
        REAL(q) :: DESUM
        REAL(q) :: EXHF
        INTEGER :: NRED
        LOGICAL :: LDELAY
        LOGICAL :: LRESIDUAL
        LOGICAL :: LHF
        INTEGER :: IU0, IU6
        INTEGER, OPTIONAL :: NKSTART, NKSTOP

#ifdef MPI
        CALL vtutor%error('EDDAV_FULL: the non-blocked Davidson optimizer needs&
                                     & scaLAPACK support ... ')
#endif
        CALL vtutor%error('EDDAV_FULL: the non-blocked Davidson optimizer has not&
                                     & been implemented in the serial version ... ')

      END SUBROUTINE

      END MODULE david_full

#endif // MPI && ScaLAPACK
