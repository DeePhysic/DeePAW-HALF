#include "symbol.inc"
      MODULE rot
      USE prec
      CONTAINS
!************************ SUBROUTINE EDWAV *****************************
! RCS:  $Id: rot.F,v 1.2 2001/02/20 14:45:00 kresse Exp $
!
!> Subroutine for optimizing all bands simultaneous.
!>
!> The subroutine calculates a search direction into which the
!> total energy is minimized, than a trial step into this direction
!> is done.
!> On the next call the new exact energy must be supplied by the calling
!> routine, EDWAV then tries to determine an optimum step length
!> and steps into the minimum.
!>
!> IFLAG determines the status of the routine: \n
!>  0   first trial for ionic configuration - steepest descent step \n
!>  1   trial - conjugate gradient step (IALGO can overwrite this mode to 0) \n
!>  2   go to the  to minimum.
!>      TOTEN must contain the old energy TOTEN2 the energy after the step \n
!>  10  increase trial step (mainly for checking accuracy of routine) \n
!> on exit \n
!> >0   main routine performs another energy-evaluation \n
!>  0   main routine ends the loop \n
!>
!> IALGO determines the specific algorithm \n
!>  3   damped MD with damping term automatically determined
!>      by the given timestep \n
!>  4   damped MD (velocity quench or quickmin) \n
!>  5   steepest descent without optimizing the stepsize \n
!>  7   steepest descent with optimization of stepsize \n
!>  8   conjugated gradient \n
!>  9   conjugate gradient with gradient information only \n
!>
!> IROT determines if subspace rotation and/or wavefunctions are optimized \n
!>  1   only orbital coefficients are optimized \n
!>  2   only rotations in the set of calculated orbitals (sub space rotation) \n
!>  3   orbital coefficient optimization and subspace rotation \n
!>  4   orbital coefficient optimization and subspace rotation.
!>      The Hamilton matrix in the subspace is determined by a self-consistent cycle
!>      (subrot_scf). \n
!>  5   orbitals optimized, exact subspace rotation in the subspace spanned by
!>      the calculated orbitals.
!>      The rotation matrix is determined by scf iterations in the subset of
!>      calculated orbitals  (often unstable) \n
!>  6   optimize pseudo "Hamiltonian" matrix eta  instead of subspace rotation \n
!>  7   optimize pseudo "Hamiltonian" matrix eta  instead of subspace rotation \n
!>      the Hamilton matrix in the subspace is determined by a self-consistent cycle
!>
!> ICEL determines how partial occupancies are optimized \n
!>  0   no optimization \n
!>  1   standard optimization (optimize a new set of support variables
!>       that are stored in W_F%CELTOT) \n
!>  2   use pseudo Hamiltonian matrix eta to optimize eigenvalues
!>      in this case W_F%CELTOT is returned and occupancies need to be
!>      calculated from W_F%CELTOT (similar to ICEL=1) \n
!> IRET determines which algorithm has been used \n
!>  0   steepest descent \n
!>  1   conjugate gradient \n
!>  3   damped \n
!>  4   quick-min \n
!>  5   trial step in conjugate gradient to determine optimal step using gradient information \n
!>
!>
!> DESUM1 expected energy change from interpolation along line search
!>
! gK: 26.02.2013 bug fix for NC potentials (in two places)
!
!***********************************************************************

 SUBROUTINE EDWAV(HAMILTONIAN,KINEDEN,&
      INFO,XC,IFLAG,IROT,ICEL,IRET, IO,BTRIAL, EFERMI, &
      ORT, RMS, E, TOTEN, TOTEN2, DESUM1, &
      GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, W,WDES,W_F,W_G, &
      LMDIM, NSIM, CQIJ,CDIJ,CHAM, CHF, SV, &
      SYMM,GRID_SOFT,GRIDC,GRIDB,GRIDUS, &
      C_TO_US,B_TO_C,SOFT_TO_C,DENCOR,CSTRF, &
      MIX,N_MIX_PAW,IRDMAX,CHTOT,RHOLM)

__GPU USE moffload

   USE ini
   USE prec
   USE wave
   USE wave_mpi
   USE wave_high
   USE lattice
   USE mpimy
   USE mgrid

   USE nonl_high
   USE poscar
   USE pseudo
   USE hamil_struct_def
   USE constant
   USE mkpoints
   USE choleski
   USE fock
   USE hamil
   USE base
   USE lattice
   USE sym_grad
   USE subrotscf
   USE pead, ONLY : PEAD_REINIT_EDWAV,PEAD_ACC_CALC_ALL,LPEAD_NO_SCF,W_STORE,PEAD_ACC_ADD,LPEAD_EFIELD_SWITCHED_ON,LUSEPEAD
   USE tau_mu
   USE setexm_struct_def, ONLY : xc_info
   USE fermi_energy, ONLY : DELSTP
   USE string, ONLY: str
   USE tutor, ONLY: vtutor
   USE solvation
   USE constrained_m_modular
#ifdef fock_dblbuf
   USE fock_dbl
#endif
   IMPLICIT NONE

   TYPE (ham_handle)  HAMILTONIAN
   TYPE (tau_handle)  KINEDEN
   TYPE (grid_3d)     GRID
   TYPE (latt)        LATT_CUR
   TYPE (type_info)   T_INFO
   TYPE (potcar)      P(T_INFO%NTYP)
   TYPE (nonlr_struct)NONLR_S
   TYPE (nonl_struct) NONL_S
   TYPE (wavespin)    W     ! current wavefunctions
   TYPE (wavespin)    W_F   ! last search direction
   TYPE (wavespin)    W_G   ! current gradient
   ! W_F%CELEN is special, and stores variational eigenvalues
   ! from which partial occupancies are set
   TYPE (wavefun1)    W1(NSIM),WSEARCH
   TYPE (wavedes)     WDES
   TYPE (wavedes1)    WDES1
   TYPE (kpoints_struct) KPOINTS
   TYPE (info_struct) INFO
   TYPE (xc_info)     XC
   TYPE (in_struct)   IO
   TYPE (mixing)      MIX
   TYPE (symmetry)    SYMM
   TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
   TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
   TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
   TYPE (grid_3d)     GRIDB      ! Broyden grid
   TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
   TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
   TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
   TYPE (energy)      E

   INTEGER IFLAG,IROT,ICEL,IRET
   INTEGER IRDMAX
   INTEGER N_MIX_PAW
   REAL(q) BTRIAL,ORT,RMS,TOTEN,TOTEN2,DESUM1,TOTEN_LOCAL

   RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
   COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)! structure factor
   COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
   REAL(q)     RHOLM(N_MIX_PAW,WDES%NCDIJ)   !  paw sphere charge density

   LOGICAL LREAL,LOVERL
   LOGICAL LCHCON                  ! update of charge or not
   REAL(q) EFERMI                  ! Fermi-level
   REAL(q) WEIMIN                  ! minimum occupancy for empty state
   INTEGER NSIM                    ! simultaneously optimised bands
   INTEGER LMDIM                   ! dimension of arrays CQIJ and CDIJ

   RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
   OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
   !-----Arrays for gradient
   GDEF CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)  ! sub space Hamiltonian in space spanned by orbitals
   GDEF CHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)   ! search direction in the space spanned by orbitals
   GDEF CROT, CB_
! work array
   COMPLEX(q) CGRAD(WDES%NRPLWV)
   GDEF,ALLOCATABLE,TARGET :: CPROW(:,:)
   COMPLEX(q),ALLOCATABLE,TARGET:: CF(:,:)    ! stores the plane wave coefficients of the gradient vector
   GDEF,ALLOCATABLE,TARGET::   CPROF(:,:)     ! stores the projected coefficients corresponding to CF
   COMPLEX(q),ALLOCATABLE,TARGET:: CFOLD(:,:) ! stores the plane wave coefficients of the previous search
   GDEF,ALLOCATABLE :: CHAMP(:,:,:,:)         ! stores the preconditioned gradient
   TYPE (wavefuna)    WA             ! pointer to W%CW(:,:;NK,ISP)
   TYPE (wavefuna)    WA_F           ! pointer to W_F%CW(:,:;NK,ISP)
   COMPLEX(q), POINTER :: CF_RED(:,:)
   GDEF   , POINTER ::    CPROW_RED(:,:),CPROF_RED(:,:)
   GDEF, ALLOCATABLE :: COVL(:,:), CUNI(:,:), CEIDB(:,:)

   ! ILOEW=1 uses an approximate rotation matrix R based on Loewdin perturbation theory
   ! ILOEW=2 uses  e^iR to build a rotation matrix for a finite step
   ! usually there is not much difference between both calles
   INTEGER,PARAMETER :: ILOEW=2

   ! this parameter determines whether line searches are performed
   ! with high precision using the gradient along the search direction
   LOGICAL,PARAMETER :: ITERATE_PREC=.FALSE.

   ! the search direction can be renormalized
   ! such that a unit step into the prec. gradient direction leads directly
   ! to the minimum for each band
   ! this is useful if the charge is not updated (LCHCON)
   ! but it seems make calculations faster in any case
   ! if MAXRENORM is 0, no renormalization is performed
   REAL(q), PARAMETER :: MAXRENORM=0.10_q, STEP_SCALING=1.5_q
   LOGICAL :: LBRACKET  ! found increase of energy along search direction
   ! scaling factor for pseudo-Hamiltonian
   ! note: this is not applied for the sophisticated algorithm where the
   ! pseudo Hamiltonian is determined by subrot
   REAL(q) :: KAPPA=0.5

   INTEGER,PARAMETER :: NSTEP=100
   LOGICAL,SAVE :: LSLOT(NSTEP)
   LOGICAL :: L_IS_UPDATED ! see subroutine slot_in
   REAL(q), SAVE :: STEP(NSTEP),ENRGY(NSTEP)
   INTEGER,SAVE :: ICOUNT=-1, NCOUNT, NCOUNT_BASE, NCOUNT_NRESET
   REAL(q),SAVE :: GNORML=0, GNORM=0, AVERAG=0, BSTEP=0, FRICTION
   REAL(q),SAVE :: DESUM0, GAMMA, GSL, BOPTH0, STRIAL1, STRIAL2
   INTEGER :: INDEX(1)

   ! these two arrays are required to store the redistributed wavefunctions
   REAL(q)    :: CELTOT(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN)
   COMPLEX(q) :: CDCHF

   ! translation of NBANDSLOW and NBANDSHIGH to local band indices
   INTEGER :: NBANDSLO,NBANDSLO_LOCAL,NBANDSHI_LOCAL

   ! local variables
   COMPLEX(q) :: CSCPD,CA,CB,CE
   REAL(q)    :: WSCAL,DE,WSUM,EFERMI_SHIFT,WEIGHT,X,DFUN,SFUN,DER,DER4
   REAL(q),SAVE    :: ORT1,ORT2,DE1,DE2,ORTCEL,DECEL
   REAL(q)    :: BCSTEP,BOPT,BOPTH,BOPT4,DESUM4,BOPTS,DESUMS,DERS,DORT,THETA,A1,A2,WMIN,DECEL_in,ORTCEL_in,DIFCEL
   INTEGER    :: IALGO,IFLAG_IN,NBANDS,NB_TOT,NODE_ME,IONODE,NCPU,NNMAX,ISP,IDUMP,NBANDSK,NSIM_LOCAL,NSIM_FOCK,NGVECTOR
   INTEGER    :: I,N,M,NN,N1,N2,NK,NP,NPRO_,NSIM_,ISPINOR
   INTEGER    :: NPOS_RED,NSTRIP_RED,IFAIL

   INTEGER :: MYKSTART, MYKINCR

   COMPLEX(q), EXTERNAL ::  WZDOTC

   ! this flag starts the line search always from the initial point
   ! should yield more consistent energies
   LOGICAL :: restore_W=.TRUE.
GPACC !$ACC ROUTINE(DELSTP) SEQ

   PROFILING_START('edwav')

   IALGO=INFO%IALGO; LOVERL=INFO%LOVERL; LREAL=INFO%LREAL
   IF (IALGO<=5) THEN
     restore_W=.FALSE.
   ELSE
     restore_W=.TRUE.
   ENDIF

   WEIMIN=0; LCHCON=INFO%LCHCON
!  reset the calling counter
   IF (IFLAG==0)  ICOUNT=-1

   CALL PEAD_REINIT_EDWAV(ICOUNT,GNORML,GNORM,AVERAG,BSTEP)
!=======================================================================
! initialise the required variables for MPI
!=======================================================================
#ifdef MPI
   NODE_ME =WDES%COMM_KIN%NODE_ME
   IONODE  =WDES%COMM_KIN%IONODE
   NCPU    =WDES%COMM_INTER%NCPU     ! number of groups (each group holds one band)
   MYKSTART=WDES%COMM_KINTER%NODE_ME
   MYKINCR =WDES%COMM_KINTER%NCPU
#else
   NODE_ME =1
   IONODE  =1
   NCPU    =1
   MYKSTART=1
   MYKINCR =1
#endif

#ifdef OFFLOADING
   CALL OFFLOAD_ON_PUSH(.TRUE.)

   IF (LSYMGRAD) CALL OFFLOAD_ON_PUSH(.FALSE.)

   CALL COPYIN_TYPED_VAR(WDES)
   CALL COPYIN_TYPED_VAR(W)
   CALL COPYIN_TYPED_VAR(W_F)
   CALL COPYIN_TYPED_VAR(W_G)
   CALL COPYIN_TYPED_VAR(NONL_S)
   CALL COPYIN_TYPED_VAR(NONLR_S)
   CALL COPYIN_TYPED_VAR(GRID)
   CALL COPYIN_TYPED_VAR(GRIDC)
   CALL COPYIN_TYPED_VAR(P)

GPACC !$ACC ENTER DATA COPYIN(DENCOR,CSTRF) __IF_ASYNC__
GPACC !$ACC ENTER DATA COPYIN(HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(LATT_CUR,SV,CDIJ,CQIJ) __IF_ASYNC__
GPACC !$ACC ENTER DATA COPYIN(CHAM,CHF) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(TO:DENCOR,CSTRF) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LATT_CUR,SV,CDIJ,CQIJ) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:CHAM,CHF) IF(OFFLOAD_ON)
#endif
!=======================================================================
! limit scaling
! unfortunately level reordering might spoil somewhat for the new algorithm
!=======================================================================
   CALL MRG_AUX(WDES,W)
GPACC !$ACC ENTER DATA CREATE(WSCAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WSCAL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(W,WSCAL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W,WSCAL) IF(OFFLOAD_ON)
   WSCAL=SUM(W%AUXTOT)/SIZE(W%AUXTOT)
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

   ! W%AUXTOT=W%AUXTOT/WSCAL ; WSCAL = 1  ! this line makes the average scaling 1
GPACC !$ACC KERNELS PRESENT(W,WSCAL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W,WSCAL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
   W%AUXTOT=MAX(1/4._q*WSCAL,MIN(4._q*WSCAL,W%AUXTOT))
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC EXIT DATA DELETE(WSCAL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WSCAL) IF(OFFLOAD_ON)

!=======================================================================
! number of bands treated simultaneously this must be a multiple of NCPU
!=======================================================================
   NSIM_LOCAL=NSIM/NCPU  ! number of bands optimised on this node
   IF (NSIM_LOCAL*NCPU /= NSIM) THEN
      CALL vtutor%bug("internal ERROR in EDWAV: NSIM is not correct " // str(NSIM), __FILE__, __LINE__)
   ENDIF

!  initialized average step size to trial step size
   IF (AVERAG==0) AVERAG=BTRIAL
   BTRIAL=AVERAG
   IDUMP=1
#ifdef MPI
   IF (NODE_ME/=IONODE) IDUMP=0
#endif

GPACC !$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1) IF(OFFLOAD_ON)
   CALL SETWDES(WDES,WDES1,0)

   IFLAG_IN=IFLAG

   IF (IFLAG_IN<2) THEN
GPACC !$ACC ENTER DATA CREATE(W1(:)) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W1) IF(OFFLOAD_ON)
      DO I=1,NSIM_LOCAL
         CALL NEWWAV(W1(I) , WDES1, .TRUE.)
      ENDDO
GPACC !$ACC ENTER DATA CREATE(WSEARCH) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WSEARCH) IF(OFFLOAD_ON)
      CALL NEWWAV(WSEARCH, WDES1, .TRUE.)
GPACC !$ACC ENTER DATA CREATE(CELTOT) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CELTOT) IF(OFFLOAD_ON)
   ENDIF

   ALLOCATE(COVL(WDES%NB_TOT,WDES%NB_TOT),CUNI(WDES%NB_TOT,WDES%NB_TOT),CEIDB(WDES%NB_TOT,WDES%NB_TOT))
GPACC !$ACC ENTER DATA CREATE(COVL,CUNI,CEIDB) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:COVL,CUNI,CEIDB) IF(OFFLOAD_ON)

   NB_TOT=WDES%NB_TOT
   NBANDS=WDES%NBANDS

   NBANDSLO=1
   IF (WDES%NBANDSLOW>1) &
      NBANDSLO=WDES%NBANDSLOW/WDES%NB_PAR+1

   NBANDSLO_LOCAL=1
   IF (WDES%NBANDSLOW>1) &
      NBANDSLO_LOCAL=CEILING(REAL((WDES%NBANDSLOW -WDES%NB_LOW)/WDES%NB_PAR))+1

   NBANDSHI_LOCAL=WDES%NBANDS
   IF (WDES%NBANDSHIGH>0.AND.WDES%NBANDSHIGH<NB_TOT) &
      NBANDSHI_LOCAL=FLOOR(REAL((WDES%NBANDSHIGH-WDES%NB_LOW)/WDES%NB_PAR))+1

!*********************************************************************
!   TRIAL STEP (or PREDICTOR STEP)
! ) calculate the gradient vector and store
!   preconditioned gradient in W_G
! ) calculate change in energy if step into direction W_G is done DESUM
! ) calculate norm of gradient-vector GNORM
! ) test orthonormality with old direction ORT
! ) do trial step
!***********************************************************************
   main: IF (IFLAG<2) THEN
      E%EXHF_ACFDT=0

      ALLOCATE(CPROW(WDES%NPROD,NBANDS))
      ALLOCATE(CF(WDES%NRPLWV,NSIM_LOCAL),CFOLD(WDES%NRPLWV,NSIM_LOCAL),CPROF(WDES%NPROD,NSIM_LOCAL))
GPACC !$ACC ENTER DATA CREATE(CPROW,CF,CFOLD,CPROF) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROW,CF,CFOLD,CPROF) IF(OFFLOAD_ON)

      ORT1  =0
      ORT2  =0
      ORTCEL=0
      DE1   =0
      DE2   =0
      DECEL =0
      CDCHF =0
!=======================================================================
! first loop
! calculate gradient (for wavefunctions)
!=======================================================================
      CALL PEAD_ACC_CALC_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO)

!=======================================================================
! start with HF part and store results W_G
!=======================================================================
GPACC !$ACC KERNELS PRESENT(CHAM,W_G) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHAM,W_G) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CHAM=0
      W_G%CW=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      IF (LHFCALC) THEN
#ifndef fock_dblbuf
      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
      CALL SETWDES(WDES,WDES1,NK)

      IF (NONLR_S%LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF

#ifndef _OPENMP
      NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))

      DO NN=1,NBANDSK,NSIM_LOCAL
         NNMAX=MIN(NN+NSIM_LOCAL-1,NBANDSK)
         IF (LPEAD_NO_SCF()) THEN
            DO N=NN,NNMAX
               NP=N-NN+1; CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),W1(NP))
            ENDDO
            CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE, &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    W_G%CW(:,NN:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, W1, &
                    LSYMGRAD=LSYMGRAD, EXHF_ACFDT= E%EXHF_ACFDT)
         ELSE
            CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    W_G%CW(:,NN:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, &
                    LSYMGRAD=LSYMGRAD, EXHF_ACFDT= E%EXHF_ACFDT)
         ENDIF
      ENDDO
#else
      IF (LPEAD_NO_SCF()) THEN
         NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))
         DO NN=1,NBANDSK,NSIM_LOCAL
            NNMAX=MIN(NN+NSIM_LOCAL-1,NBANDSK)
            DO N=NN,NNMAX
               NP=N-NN+1; CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),W1(NP))
            ENDDO
            CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE, &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    W_G%CW(:,NN:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, W1, &
                    LSYMGRAD=LSYMGRAD, EXHF_ACFDT= E%EXHF_ACFDT)
         ENDDO
      ELSE
         NSIM_FOCK=NBLOCK_FOCK
         DO NN=1,W%WDES%NB_TOTK(NK,ISP),NSIM_FOCK
            NNMAX=MIN(NN+NSIM_FOCK-1,W%WDES%NB_TOTK(NK,ISP))
            CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    W_G%CW(:,:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, &
                    LSYMGRAD=LSYMGRAD, EXHF_ACFDT= E%EXHF_ACFDT, &
                    LCOMMENSURATE=.FALSE.)
         ENDDO
      ENDIF
#endif
      ENDDO
      ENDDO

      CALLMPI( M_sum_z(WDES1%COMM_KIN,CDCHF,1))
      CALLMPI( M_sum_z(WDES%COMM_KINTER, CDCHF, 1))
!      WRITE(*,*) '-1/2 Fock energy is ',CDCHF
      E%EXHF=CDCHF
      CALLMPI( M_sum_d(WDES1%COMM_KIN, E%EXHF_ACFDT ,1))
      CALLMPI( M_sum_d(WDES%COMM_KINTER, E%EXHF_ACFDT, 1))
#else
      IF (LPEAD_NO_SCF()) THEN
         CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
        &   EX=E%EXHF,EX_ACFDT=E%EXHF_ACFDT,LSYMGRAD=LSYMGRAD,XI=W_G,WP=W)
      ELSE
         CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
        &   EX=E%EXHF,EX_ACFDT=E%EXHF_ACFDT,LSYMGRAD=LSYMGRAD,XI=W_G)
      ENDIF
#endif
      IF (LSYMGRAD) THEN
        CALL APPLY_SMALL_SPACE_GROUP_OP( W, W_G, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1)
      ENDIF
      ENDIF

!=======================================================================
! DFT part optimized to yield minimum number of operations
!=======================================================================
spin:   DO ISP=1,WDES%ISPIN
kpoint: DO NK=MYKSTART,WDES%NKPTS,MYKINCR
      NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))

      CALL SETWDES(WDES,WDES1,NK)

      WA=ELEMENTS(W, WDES1, ISP)
      WA_F=ELEMENTS(W_F, WDES1, ISP)

__GPU CALL COPYIN_TYPED_VAR(WA_F)

      IF (WDES%DO_REDIS) THEN
        CALL SET_GPOINTER(CPROW_RED, WDES1%NPROD_RED, NB_TOT, CPROW(1,1))
        CALL SET_WPOINTER(CF_RED,   WDES1%NRPLWV_RED, NSIM, CF(1,1))
        CALL SET_GPOINTER(CPROF_RED, WDES1%NPROD_RED, NSIM, CPROF(1,1))
      ELSE
        CPROW_RED => CPROW(:,:)
        CF_RED    => CF(:,:)
        CPROF_RED => CPROF(:,:)
      ENDIF

      NGVECTOR=WDES1%NGVECTOR
      CALL OVERL(WDES1, LOVERL,LMDIM,CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW(1,1))

!-----------------------------------------------------------------------
! avoid conjugation of those bands which do not contribute to energy
!-----------------------------------------------------------------------
      ! wave functions coefficients are redistributed in W_F
      IF (IALGO>=8) THEN
GPACC !$ACC PARALLEL LOOP GANG PRESENT(WDES1,W,WA_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WDES1,W,WA_F) IF(OFFLOAD_ON)
         DO N1=1,NB_TOT
            IF (ABS(W%FERTOT(N1,NK,ISP))<WEIMIN) THEN
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
               DO M=1,WDES1%NPL_RED
                  WA_F%CW_RED(M,N1)=0
               ENDDO
            ENDIF
         ENDDO
GPACC !$ACC PARALLEL LOOP GANG PRESENT(WDES,W,WA_F) COPYIN(W%FERWE(:,NK,ISP)) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:WDES,W,WA_F) IF(OFFLOAD_ON)
         ! same for character (which is distributed differently on the nodes)
         DO N=1,NBANDS
            IF (ABS(W%FERWE(N,NK,ISP))<WEIMIN) THEN
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
               DO NPRO_=1,WDES%NPRO
                  WA_F%CPROJ(NPRO_,N)=0
               ENDDO
            ENDIF
         ENDDO
      ENDIF

__GPU CALL DELETE_TYPED_VAR(WA_F)

!-----------------------------------------------------------------------
! calculate the phase-factor on the compressed grid (only non local)
!-----------------------------------------------------------------------
      IF (LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF

      ! redistribute over plane wave coefficients
      IF (WDES%DO_REDIS) THEN
        IF (LOVERL) CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
        CALL REDIS_PW  (WDES1, NBANDS, W%CW(1,1,NK,ISP))
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
      ENDIF
!=======================================================================
  bands: DO NN=NBANDSLO,NBANDSK,NSIM_LOCAL
!=======================================================================
GPACC !$ACC KERNELS PRESENT(CF,CPROF,CFOLD) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CF,CPROF,CFOLD) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CF=0
      CPROF=0
      CFOLD=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      NNMAX=MIN(NN+NSIM_LOCAL-1,NBANDSK)
      NSIM_=NNMAX-NN+1

      DO N=NN,NNMAX ! each step in the loop copies NSIM_LOCAL/NSIM bands !
         NP=N-NN+1

         CALL WZCOPY(WDES%NRPLWV,W_G%CW(1,N,NK,ISP),1,CF(1,NP),1)

         CALL WZCOPY(WDES%NRPLWV,W_F%CW(1,N,NK,ISP),1,CFOLD(1,NP),1)

         CALL WZCOPY(WDES%NRPLWV,W%CW(1,N,NK,ISP),1,W1(NP)%CW(1),1)
         CALL WGCOPY(WDES%NPROD,W%CPROJ(1,N,NK,ISP),1,W1(NP)%CPROJ(1),1)

        ! redistribute over bands
         IF (WDES%DO_REDIS) THEN
            CALL REDIS_PW  (WDES1, 1, CFOLD(1,NP))
            CALL REDIS_PW  (WDES1, 1, W1(NP)%CW(1))
            CALL REDIS_PROJ(WDES1, 1, W1(NP)%CPROJ(1))
         ENDIF

         DO ISPINOR=0,WDES%NRSPINORS-1
            CALL FFTWAV(NGVECTOR, WDES%NINDPW(1,NK),W1(NP)%CR(1+ISPINOR*GRID%MPLWV),W1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
         ENDDO

#ifndef oldhamil
         CALL HAMILT_LOCAL_(HAMILTONIAN, LATT_CUR, W1(NP), SV, ISP, CF(:,NP), LHFCALC)
#else
         IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
            CALL HAMILT_LOCAL_TAU(W1(NP), SV, LATT_CUR, HAMILTONIAN%MU, ISP, CF(:,NP), LHFCALC)
         ELSE
            CALL HAMILT_LOCAL(W1(NP), SV, ISP, CF(:,NP), LHFCALC)
         ENDIF
#endif

         WDES1%NBANDS=1    ! WDES1%NBANDS is used only here to fake OVERL
         CALL OVERL(WDES1, .TRUE., LMDIM, CDIJ(1,1,1,ISP), W1(NP)%CPROJ(1), CPROF(1,NP))

         CALL PEAD_ACC_ADD(CF(:,NP),CPROF(:,NP),N,NK,ISP)

         IF ((W%WDES%NB_LOW+(N-1)*W%WDES%NB_PAR)>W%WDES%NB_TOTK(NK,ISP)) THEN
GPACC !$ACC KERNELS PRESENT(CF,CPROF,W1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CF,CPROF,W1) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            CF(:,NP)=0
            CPROF(:,NP)=0
            W1(NP)%CW=0
            W1(NP)%CR=0
            W1(NP)%CPROJ=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
         ENDIF

         ! redistribute over plane wave coefficients
         IF (WDES%DO_REDIS) THEN
            CALL REDIS_PW  (WDES1, 1, CF(1,NP))
            CALL REDIS_PROJ(WDES1, 1, CPROF(1,NP))
         ENDIF
      ENDDO

! calculate the lower triangle of the Hamilton matrix
! in the subspace spanned by the present set of wavefunctions

      NPOS_RED  =(NN-1)*NCPU+1
      NSTRIP_RED=NSIM_*NCPU

      IF (NBANDSLO>1) THEN
         ! in this case it is not possible to calculate only the lower triangle
         CALL ORTH2( &
            WA%CW_RED(1,1),CF(1,1),WA%CPROJ_RED(1,1),CPROF(1,1), &
            NB_TOT,NPOS_RED, NSTRIP_RED, WDES1%NPL_RED, WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED, &
            CHAM(1,1,NK,ISP))

         CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,NPOS_RED,NK,ISP),NB_TOT*NSTRIP_RED))

GPACC !$ACC PARALLEL LOOP PRESENT(W,CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:W,CHAM) IF(OFFLOAD_ON)
         DO N1=NPOS_RED,NPOS_RED+NSTRIP_RED-1
            W%CELTOT(N1,NK,ISP)=CHAM(N1,N1,NK,ISP)
         ENDDO
      ELSE
         ! calculate the lower triangle of the Hamilton matrix
         ! in the subspace spanned by the present set of wavefunctions
         CALL ORTH1( &
            'L',WA%CW_RED(1,1),CF(1,1),WA%CPROJ_RED(1,1),CPROF(1,1), &
            NB_TOT,NPOS_RED,NSTRIP_RED,WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED, &
            CHAM(1,1,NK,ISP))

         CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,NPOS_RED,NK,ISP),NB_TOT*NSTRIP_RED))

         ! fill in upper triangle
         ! in the parallel version part of the upper triangle has been destroyed
         ! and needs to be refilled
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CHAM) IF(OFFLOAD_ON)
         DO N1=NPOS_RED,NPOS_RED+NSTRIP_RED-1
            DO N2=1,NPOS_RED-1
               CHAM(N2,N1,NK,ISP)=GCONJG(CHAM(N1,N2,NK,ISP))
            ENDDO
         ENDDO
GPACC !$ACC PARALLEL LOOP PRESENT(W,CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:W,CHAM) IF(OFFLOAD_ON)
         DO N1=NPOS_RED,NPOS_RED+NSTRIP_RED-1
            W%CELTOT(N1,NK,ISP)=CHAM(N1,N1,NK,ISP)
         ENDDO
      ENDIF

      ! set all matrix elements beyond NB_TOTK to zero
      IF (W%WDES%NB_TOTK(NK,ISP)<NB_TOT) THEN
GPACC !$ACC KERNELS PRESENT(W,CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W,CHAM) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         CHAM(:,W%WDES%NB_TOTK(NK,ISP)+1:NB_TOT,NK,ISP)=0
         CHAM(W%WDES%NB_TOTK(NK,ISP)+1:NB_TOT,:,NK,ISP)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ENDIF

!
! subtract from H phi the components parallel to the subspace
! spanned by the other orbitals in the set
      IF ( WDES1%NPL_RED/=0) &
! subtract - < phi_n1 | H_local  | phi_n2 >  | phi_n2>
      CALL WGGEMM('N', 'N', m_  WDES1%NPL_RED , NSTRIP_RED, NB_TOT, -one, &
     &               WA%CW_RED(1,1), m_ WDES1%NRPLWV_RED, CHAM(1,NPOS_RED,NK,ISP), NB_TOT,  &
     &               one, CF(1,1), m_ WDES1%NRPLWV_RED)

! subtract < phi_n1 |  H_local | phi_n2 > Q_ij  < p_j| phi_n2 > from orbital character
      IF (LOVERL) THEN
         CALL WGGEMM('N', 'N',  WDES1%NPRO_RED, NSTRIP_RED, NB_TOT, -one, &
     &               CPROW_RED(1,1), WDES1%NPROD_RED, CHAM(1,NPOS_RED,NK,ISP), NB_TOT,  &
     &               one, CPROF(1,1), WDES1%NPROD_RED)
      ENDIF

      ! redistribute the gradient over bands
      DO N=NN,NNMAX
         NP=N-NN+1
         IF (WDES%DO_REDIS) THEN
            CALL REDIS_PW  (WDES1, 1, CF(1,NP))
            CALL REDIS_PROJ(WDES1, 1, CPROF(1,NP))
         ENDIF
      ENDDO

    innerband: DO N=NN,NNMAX
      IF ( (W%WDES%NB_LOW+(N-1)*W%WDES%NB_PAR)>W%WDES%NB_TOTK(NK,ISP)) CYCLE

      NP=N-NN+1
!-----------------------------------------------------------------------
! add the non local contribution to CF
!-----------------------------------------------------------------------
      IF (LREAL) THEN
GPACC !$ACC KERNELS PRESENT(WSEARCH) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:WSEARCH) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         WSEARCH%CR=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#ifdef OFFLOADING
         IF (OFFLOAD_ON) THEN
            CALL RACC0_ACC(NONLR_S,WDES1,CPROF(1,NP),WSEARCH%CR(1))
         ELSE
#endif
            CALL RACC0(NONLR_S,WDES1,CPROF(1,NP),WSEARCH%CR(1))
#ifdef OFFLOADING
         ENDIF
#endif
         DO ISPINOR=0,WDES1%NRSPINORS-1
            CALL FFTEXT(NGVECTOR,WDES%NINDPW(1,NK),WSEARCH%CR(1+ISPINOR*WDES1%GRID%MPLWV),CF(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
         ENDDO
      ELSE
        CALL VNLAC0(NONL_S,WDES1,CPROF(1,NP),CF(1,NP))
      ENDIF

      WEIGHT =2*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP) ! times two to include the c.c.
!-----------------------------------------------------------------------
! perform preconditioning of gradient
! CF is the current "gradient" WSEARCH the direction to go
!-----------------------------------------------------------------------
      CALL WZCOPY(WDES%NRPLWV,CF(1,NP),1,WSEARCH%CW(1),1)

      CALL SIMPLE_PRECOND( WSEARCH )
!----------------------------------------------------------------------
!  calculate the character of the wavefunctions
!----------------------------------------------------------------------
      DO ISPINOR=0,WDES1%NRSPINORS-1
         CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),WSEARCH%CR(1+ISPINOR*WDES1%GRID%MPLWV),WSEARCH%CW(1+ISPINOR*NGVECTOR),GRID)
      ENDDO
      IF (LREAL) THEN
#ifdef OFFLOADING
         IF (OFFLOAD_ON) THEN
            CALL RPRO1_ACC(NONLR_S,WDES1,WSEARCH)
         ELSE
#endif
            CALL RPRO1(NONLR_S,WDES1,WSEARCH)
#ifdef OFFLOADING
         ENDIF
#endif
      ELSE
         CALL PROJ1(NONL_S,WDES1,WSEARCH)
      ENDIF

!-----------------------------------------------------------------------
! following statements allow to renormalize the trial step
! such that a unit step into the direction of the precond gradient
! steps directly into the minimum
!-----------------------------------------------------------------------
      ! project out current wave function

      CALL PROJCN(WSEARCH,W1(NP),CQIJ, ISP, CSCPD)
      IF (MAXRENORM/=0.AND..NOT.LPEAD_EFIELD_SWITCHED_ON()) THEN
!     IF (MAXRENORM/=0) THEN
#if PGI_BEFORE_XX_X
GPACC !$ACC UPDATE SELF(W%AUX(N:N,NK:NK,ISP:ISP)) __IF_ASYNC__ ! W%AUX(N,NK,ISP) is needed on the host (below)
#else
GPACC !$ACC UPDATE SELF(W%AUX(N,NK,ISP)) __IF_ASYNC__ ! W%AUX(N,NK,ISP) is needed on the host (below)
#endif
GPOMP !$OMP TARGET UPDATE FROM(W%AUX(N,NK,ISP)) IF(OFFLOAD_ON)
         ! renormalize the wave function
         CALL CNORMN(WSEARCH,CQIJ, ISP, WSCAL)

         ! recalculate the wavefunction in real space
         CALL CNORMN_REAL(WSEARCH, W1(NP), ISP, WSCAL, CSCPD )

         ! determine 2x2 sub space Hamilton matrix and diagonalize it
         ! use local contributions only
         CB_=WZDOTC(WDES1%NPL,WSEARCH%CW(1),1,CF(1,NP),1)
         CALLMPI( M_sum_g(WDES%COMM_INB, CB_, 1))

#ifndef oldhamil
         CALL ECCP_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1(NP),   W1(NP), LMDIM, CDIJ, SV, ISP, CE)
         CALL ECCP_(HAMILTONIAN, GRID, LATT_CUR, WDES1, WSEARCH, WSEARCH, LMDIM, CDIJ, SV, ISP, CA)
         CALL ECCP_(HAMILTONIAN, GRID, LATT_CUR, WDES1, W1(NP),  WSEARCH, LMDIM, CDIJ, SV, ISP, CB)
#else
         IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
            CALL ECCP_TAU(WDES1,W1(NP) ,W1(NP) ,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),CE)
            CALL ECCP_TAU(WDES1,WSEARCH,WSEARCH,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),CA)
            CALL ECCP_TAU(WDES1,W1(NP) ,WSEARCH,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),CB)
         ELSE
            CALL ECCP(WDES1,W1(NP) ,W1(NP),  LMDIM,CDIJ(1,1,1,ISP),GRID, SV(1,ISP), CE)
            CALL ECCP(WDES1,WSEARCH,WSEARCH, LMDIM,CDIJ(1,1,1,ISP),GRID, SV(1,ISP), CA)
            CALL ECCP(WDES1,W1(NP) ,WSEARCH, LMDIM,CDIJ(1,1,1,ISP),GRID, SV(1,ISP), CB)
         ENDIF
#endif

         ! determine optimal step size
         DIFCEL =CA-CE
         CROT   =CB/DIFCEL

         THETA=.5_q*ATAN(2*ABS(CROT))

         CROT =  SIN(THETA)*CROT/MAX(ABS(CROT),1E-8_q)

         ! here the real important trick
         ! CROT*WSCAL is a kind of preconditioner by which we should scale
         ! the gradient, but only the average is sensible
         ! in each single step the quantity is pretty useless and spoiles
         ! the conjugate gradient cycle

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for W%AUX(N,NK,ISP)
         ! maximum allowed change given by MAXRENORM
         W%AUX(N,NK,ISP)=MIN(W%AUX(N,NK,ISP)*(1+MAXRENORM),MAX(W%AUX(N,NK,ISP)*(1-MAXRENORM), &
                         W%AUX(N,NK,ISP)*(1-MAXRENORM)+REAL(CROT,q)*WSCAL*MAXRENORM))
#if PGI_BEFORE_XX_X
GPACC !$ACC UPDATE DEVICE(W%AUX(N:N,NK:NK,ISP:ISP)) __IF_ASYNC__
#else
GPACC !$ACC UPDATE DEVICE(W%AUX(N,NK,ISP)) __IF_ASYNC__
#endif
GPOMP !$OMP TARGET UPDATE TO(W%AUX(N,NK,ISP)) IF(OFFLOAD_ON)

         ! now scale the search direction by this value
         WSCAL=W%AUX(N,NK,ISP)/WSCAL
         CALL W1_DSCAL( WSEARCH, WSCAL)
         IF (IDUMP>=3) WRITE(*,'(1F9.4,"scal")',ADVANCE='NO') W%AUX(N,NK,ISP)

      ENDIF
!-----------------------------------------------------------------------
! calculate expected energy-change into search direction
! and store the (preconditioned) gradient in W_G%CW
!-----------------------------------------------------------------------
      ! gradient should be orthogonal to present wavefunction
      DORT=WZDOTC(WDES1%NPL,W1(NP)%CW(1),1,CF(1,NP),1)
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      CALLMPI( M_sum_d(WDES%COMM_INB, DORT, 1))
      IF (ABS(DORT)>1E-4) THEN
         CALL vtutor%bug("EDWAV: internal error, the gradient is not orthogonal " // str(NK) // &
            " " // str(NP) // " " // str(DORT), __FILE__, __LINE__)
      ENDIF

      ! orthogonality of gradient to previous search direction
      DORT=WZDOTC(WDES1%NPL,CFOLD(1,NP),1,CF(1,NP),1)
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      CALLMPI( M_sum_d(WDES%COMM_INB, DORT, 1))
      IF (IDUMP>=3) WRITE(*,'(F9.4,"O")',ADVANCE='NO') DORT

      IF (N>=NBANDSLO_LOCAL) ORT1=ORT1+DORT*WEIGHT*WDES%RSPIN

      ! energy change for a unit step along gradient
      IF (IROT==1 .OR. IROT>=3 ) THEN
         DE=WZDOTC(WDES1%NPL,WSEARCH%CW(1),1,CF(1,NP),1)
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         CALLMPI( M_sum_d(WDES%COMM_INB, DE, 1))
         IF (DE<0) THEN
            WRITE(*,'("ups2",I4,8F14.7)') N,REAL(W%CELEN(N,NK,ISP))-REAL(CA),CB_,CROT,THETA,DE
         ENDIF
         IF (IDUMP>=3)  WRITE(*,'(E14.5,"E")',ADVANCE='NO') DE

         IF (N>=NBANDSLO_LOCAL) DE1=DE1+DE*WEIGHT*WDES%RSPIN
         
GPACC !$ACC KERNELS PRESENT(W_G) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W_G) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         W_G%CW(:,N,NK,ISP)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
         CALL WZCOPY(WDES1%NPL,WSEARCH%CW(1),1,W_G%CW(1,N,NK,ISP),1)
         
GPACC !$ACC KERNELS PRESENT(W_G) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W_G) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         W_G%CPROJ(:,N,NK,ISP)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
         CALL WGCOPY(WDES%NPRO,WSEARCH%CPROJ(1),1,W_G%CPROJ(1,N,NK,ISP),1)

      ENDIF
      IF (IDUMP>=3)  WRITE(*,*)
!=======================================================================
    ENDDO innerband
  ENDDO bands
!=======================================================================

! set the gradients up to NBANDSLO_LOCAL to zero
      IF (NBANDSLO_LOCAL>1) THEN
GPACC !$ACC KERNELS PRESENT(W_G) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W_G) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         W_G%CW   (:,1:NBANDSLO_LOCAL-1,NK,ISP)=0
         W_G%CPROJ(:,1:NBANDSLO_LOCAL-1,NK,ISP)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ENDIF
      IF (WDES%DO_REDIS) THEN
! redistribute CPROJ over  bands
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
! gradient distributed over plane wave coefficients
        CALL REDIS_PW  (WDES1, NBANDS, W_G%CW   (1,1,NK,ISP))
      ENDIF

      ! set CELTOT too huge eigenvalues for bands beyond NB_TOTK
      ! occupancies zero
      IF (W%WDES%NB_TOTK(NK,ISP)<NB_TOT) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(W,CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:W,CHAM) IF(OFFLOAD_ON)
         DO N=W%WDES%NB_TOTK(NK,ISP)+1,NB_TOT
            W%CELTOT(N,NK,ISP)=10000._q
            W%FERTOT(N,NK,ISP)=0._q
            CHAM(N,N,NK,ISP)=10000._q
         ENDDO
      ENDIF
!=======================================================================
ENDDO kpoint
ENDDO spin
!=======================================================================

!      WRITE(*,*) '-1/2 Fock energy is ',CDCHF

      CALLMPI( M_sum_d(WDES%COMM_INTER, ORT1, 1))
      CALLMPI( M_sum_d(WDES%COMM_INTER, DE1, 1))

      CALLMPI( M_sum_d(WDES%COMM_KINTER, ORT1, 1))
      CALLMPI( M_sum_d(WDES%COMM_KINTER, DE1, 1))

      IF (WDES%DO_REDIS) THEN
         W%OVER_BAND=.TRUE.
         W_G%OVER_BAND=.TRUE.
         W_F%OVER_BAND=.TRUE.
      ENDIF
! at this point W and W_G are distributed over plane wave coefficients
! orthogonalize the search direction to current orbitals
      DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
            CALL SETWDES(WDES,WDES1,NK)

            CALL OVERL(WDES1, LOVERL,LMDIM,CQIJ(1,1,1,ISP), W_G%CPROJ(1,1,NK,ISP),CPROW(1,1))

            IF (WDES%DO_REDIS) THEN
               IF (LOVERL) CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
               CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PROJ(WDES1, NBANDS, W_G%CPROJ(1,1,NK,ISP))
            ENDIF

GPACC !$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:COVL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            COVL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
            CALL ORTH2( &
                 W%CW(1,1,NK,ISP),W_G%CW(1,1,NK,ISP),W%CPROJ(1,1,NK,ISP), &
                 CPROW(1,1),NB_TOT, &
                 1, NB_TOT, WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))

            CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))

            IF ((W%WDES%NB_TOTK(NK,ISP)+1)<=NB_TOT) THEN
GPACC !$ACC KERNELS PRESENT(COVL,W) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:COVL,W) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
               COVL(:,W%WDES%NB_TOTK(NK,ISP)+1:)=0
               COVL(W%WDES%NB_TOTK(NK,ISP)+1:,:)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
            ENDIF

            IF (WDES1%NPL_RED /=0 ) &
            CALL WGGEMM('N', 'N', m_ WDES1%NPL_RED , NB_TOT, NB_TOT, -one, &
                 W%CW(1,1,NK,ISP), m_ WDES1%NRPLWV_RED, COVL(1,1), NB_TOT,  &
                 one, W_G%CW(1,1,NK,ISP), m_ WDES1%NRPLWV_RED)

            IF ( WDES1%NPRO_RED /=0) &
            CALL WGGEMM('N', 'N', WDES1%NPRO_RED , NB_TOT, NB_TOT, -one, &
                 W%CPROJ(1,1,NK,ISP), WDES1%NPROD_RED, COVL(1,1), NB_TOT,  &
                 one, W_G%CPROJ(1,1,NK,ISP), WDES1%NPROD_RED)

            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PROJ(WDES1, NBANDS, W_G%CPROJ(1,1,NK,ISP))
            ENDIF
         ENDDO
      ENDDO
GPACC !$ACC EXIT DATA DELETE(CPROW,CF,CFOLD,CPROF) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CPROW,CF,CFOLD,CPROF) IF(OFFLOAD_ON)
      DEALLOCATE(CPROW,CF,CFOLD,CPROF)

!-----------------------------------------------------------------------
! conjugate CHAM
!  CHAM(n2,n1) = < phi_n2 | H | phi_n1 >* = < phi_n1 | H | phi_n2 >
! why this is done, I do not know. Historic reasons, but really
! this is a mess
!-----------------------------------------------------------------------
GPACC !$ACC KERNELS PRESENT(CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHAM) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CHAM=GCONJG(CHAM)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      ! clear all rows and columns below NBANDSLOW
      IF (WDES%NBANDSLOW>1) THEN
GPACC !$ACC KERNELS DEFAULT(PRESENT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHAM) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         CHAM(1:WDES%NBANDSLOW-1,:,:,:)=0
         CHAM(:,1:WDES%NBANDSLOW-1,:,:)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ENDIF

      IF (IDUMP>3) THEN
GPACC !$ACC UPDATE SELF(CHAM) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(CHAM) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
           DO NK=MYKSTART,WDES%NKPTS,MYKINCR
              CALL DUMP_HAM( "Hamilton matrix",WDES, CHAM(1,1,NK,ISP))
           ENDDO
        ENDDO
      ENDIF
!=======================================================================
! calculate energy-change due to subspacerotation
! and the orthonormality to the last rotation
! the gradient is  simply  i [H,F]
! ORT   = i [H,F] * CHF (CHF== last rotation matrix)
! DE2   = i [H,F] * new rotation matrix
!=======================================================================
      IF (IROT==4 .OR. IROT==5 .OR. IROT==7) THEN
         ALLOCATE(CHAMP(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
GPACC !$ACC ENTER DATA CREATE(CHAMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CHAMP) IF(OFFLOAD_ON)
         CALL WGCOPY(WDES%NB_TOT*WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN,CHAM(1,1,1,1),1,CHAMP(1,1,1,1),1)

         CALL SUBROT_SCF( &
              HAMILTONIAN,KINEDEN, &
              P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
              T_INFO,INFO,XC,IO,MIX,KPOINTS,SYMM, &
              GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,E, &
              DENCOR,CSTRF,N_MIX_PAW,LMDIM,IRDMAX, &
              CHTOT,RHOLM,CHAMP, IROT==4 .OR. IROT==7)

         IF (IROT==5) THEN
            DO ISP=1,WDES%ISPIN
            DO NK=MYKSTART,WDES%NKPTS,MYKINCR
               CALL SETWDES(WDES,WDES1,NK)
               IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_G%CPROJ(1,1,NK,ISP))
               CALL LINCOM('F',W_G%CW(:,:,NK,ISP),W_G%CPROJ(:,:,NK,ISP),CHAMP(1,1,NK,ISP), &
                    W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
                    WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
                    W_G%CW(:,:,NK,ISP),W_G%CPROJ(:,:,NK,ISP))
               IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_G%CPROJ(1,1,NK,ISP))

               IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
               CALL LINCOM('F',W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP),CHAMP(1,1,NK,ISP), &
                    W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
                    WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
                    W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP))
               IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
            ENDDO
            ENDDO
            ! store diagonals for later use (update of fermi-weights)
            IF (ICEL/=0) THEN
               CALL vtutor%bug("internal error in EDWAV: for IROT=5, ICEL must be 0", __FILE__, __LINE__)
            ENDIF
         ELSE
            ! store diagonals for later use (update of fermi-weights)
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(PRESENT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
GPOMP !$OMP MAP(PRESENT,TO:CELTOT,CHAMP) IF(OFFLOAD_ON)
            DO ISP=1,WDES%ISPIN
               DO NK=MYKSTART,WDES%NKPTS,MYKINCR
                  DO N=1,WDES%NB_TOT
                     CELTOT(N,NK,ISP)=CHAMP(N,N,NK,ISP)
                  ENDDO
               ENDDO
            ENDDO
            ! scale the rotation matrix by the mean value of 1/BTRIAL and 1
            ! this moves the typical rotation frequencies into the middle of the spectrum
            IF (IROT==4 .OR. IROT==5) THEN
               IF (ILOEW==1 .OR. ILOEW==2) &
               CALL ROTINI(WDES, W, W_F, CHAM, CHF, CEIDB, WEIMIN, ORT2, DE2, &
                    !  .TRUE. indicates here that the rotation matrix is not
                    ! conditioned by ROTINI
                    .TRUE., CHAMP, ((1/BTRIAL+1)/2))
            ELSE
               CALL ETAINI(WDES, W, W_F, CHAM, CHF, CEIDB, KPOINTS, EFERMI, WEIMIN, &
                    ((1/BTRIAL+1)/2), ORT2, DE2, ORTCEL, DECEL, &
                    .TRUE., CHAMP)
            ENDIF
         ENDIF
GPACC !$ACC EXIT DATA DELETE(CHAMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CHAMP) IF(OFFLOAD_ON)
         DEALLOCATE(CHAMP)
      ELSE IF (IROT==2 .OR. IROT==3) THEN
         IF (ILOEW==1 .OR. ILOEW==2) THEN
            CALL ROTINI(WDES, W, W_F, CHAM, CHF, CEIDB, WEIMIN, ORT2, DE2, &
                 LCHCON )
         ENDIF
      ELSE IF (IROT==6) THEN
         CALL ETAINI(WDES, W, W_F, CHAM, CHF, CEIDB, KPOINTS, EFERMI, WEIMIN, &
              KAPPA, ORT2, DE2, ORTCEL, DECEL, &
              LCHCON )
      ENDIF
!=======================================================================
!
! calculate gradient for fermi-weights
! actually the fermi-weights are calculated from W_F%CELTOT
! and W_F%CELTOT is optimized in the course of the calculations
! the prec. gradient for W_F%CELTOT is stored in W_G%FERTOT
! and the search direction in W_F%FERTOT
!
!=======================================================================
      IF (KPOINTS%ISMEAR>=-1 .AND. KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
         EFERMI_SHIFT=0
         WSUM =0

!      depsilon_Fermi =
!     [\sum_n df_n/d(epsilon_n-epsilon_Fermi) x d epsilon_n ]/ [\sum_n df_n/d(epsilon_n-epsilon_Fermi)]
GPACC !$ACC PARALLEL LOOP GANG REDUCTION(+:WSUM,EFERMI_SHIFT) PRESENT(WDES,W,W_F) PRIVATE(WEIGHT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE REDUCTION(+:WSUM,EFERMI_SHIFT) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F) PRIVATE(WEIGHT) IF(OFFLOAD_ON)
         DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
            WEIGHT =WDES%WTKPT(NK)*WDES%RSPIN
GPACC !$ACC LOOP VECTOR REDUCTION(+:WSUM,EFERMI_SHIFT) PRIVATE(X,DFUN,SFUN)
GPOMP !$OMP PARALLEL DO REDUCTION(+:WSUM,EFERMI_SHIFT) PRIVATE(X,DFUN,SFUN) IF(OFFLOAD_ON)
            DO N =1,W%WDES%NB_TOTK(NK,ISP)
               ! skip states below NBANDSLOW
               IF (N<WDES%NBANDSLOW) CYCLE

               X=  (EFERMI- REAL( W_F%CELTOT(N,NK,ISP) ,KIND=q) )/KPOINTS%SIGMA
               CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)
               WSUM =WSUM +DFUN* WEIGHT
               EFERMI_SHIFT=EFERMI_SHIFT+DFUN*(W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))/KPOINTS%SIGMA*WEIGHT
            ENDDO
         ENDDO
         ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for WSUM and EFERM_SHIFT
         CALLMPI( M_sum_d(WDES%COMM_KINTER,WSUM,1) )
         CALLMPI( M_sum_d(WDES%COMM_KINTER,EFERMI_SHIFT,1) )

      ! shift is the first order change in the Fermi-energy
         IF (ABS(WSUM)>1E-10) THEN
            EFERMI_SHIFT=EFERMI_SHIFT/WSUM
         ELSE
            EFERMI_SHIFT=0
         ENDIF
         WMIN=WSUM/WDES%NB_TOT*WEIMIN

         DECEL_in=DECEL
         ORTCEL_in=ORTCEL
         DECEL=0
         ORTCEL=0


GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG REDUCTION(+:DECEL,ORTCEL) PRESENT(WDES,W,W_F,W_G) PRIVATE(WEIGHT)__IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE REDUCTION(+:DECEL,ORTCEL) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F,W_G) PRIVATE(WEIGHT) IF(OFFLOAD_ON)
         DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
            WEIGHT =WDES%WTKPT(NK)*WDES%RSPIN
GPACC !$ACC LOOP VECTOR REDUCTION(+:DECEL,ORTCEL) PRIVATE(X,DFUN,SFUN)
GPOMP !$OMP PARALLEL DO REDUCTION(+:DECEL,ORTCEL) PRIVATE(X,DFUN,SFUN) IF(OFFLOAD_ON)
            DO N =1,W%WDES%NB_TOTK(NK,ISP)
               ! skip states below NBANDSLOW
               IF (N<WDES%NBANDSLOW) CYCLE

               X=  (EFERMI- REAL( W_F%CELTOT(N,NK,ISP) ,KIND=q) )/KPOINTS%SIGMA
               CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)
               DER = -DFUN/KPOINTS%SIGMA
              ! avoid conjugation of those elements which do not contribute to energy
              ! for CG algorithm
               IF (ABS(DFUN)<WMIN .AND. IALGO>=8) THEN
                  W_F%FERTOT(N,NK,ISP)=0
               ENDIF
               DFUN= DFUN*((W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))/KPOINTS%SIGMA-EFERMI_SHIFT)

               IF (IROT==4) THEN
             ! for IROT = 4 we simply set the direction such that the dummy eigenvalues
             ! move into the direction of the optimized
             ! sub space Hamiltonian (stored few lines above)
                  SFUN= (W_F%CELTOT(N,NK,ISP)-CELTOT(N,NK,ISP))
               ELSE IF (.NOT.LCHCON) THEN
             ! precond. gradient, damp response somewhat
                  SFUN= (W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))/2
               ELSE
             ! preconditioned gradient equal gradient
                  SFUN=(W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))
               ENDIF

               DECEL  =DECEL  +DFUN* SFUN*WEIGHT
               ORTCEL =ORTCEL +DFUN* W_F%FERTOT(N,NK,ISP)*WEIGHT
               W_G%FERTOT(N,NK,ISP)= SFUN
            ENDDO
         ENDDO
         ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for DECEL and ORTCEL
         CALLMPI( M_sum_d(WDES%COMM_KINTER,DECEL,1) )
         CALLMPI( M_sum_d(WDES%COMM_KINTER,ORTCEL,1) )

         DECEL=DECEL+DECEL_in
         ORTCEL=ORTCEL+ORTCEL_in

      ELSE IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
         DECEL  =0
         ORTCEL =0
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,W_G,W_F,W) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:WDES,W_G,W_F,W) IF(OFFLOAD_ON)
         DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
            DO N =1,W%WDES%NB_TOTK(NK,ISP)
               ! skip states below NBANDSLOW
               IF (N<WDES%NBANDSLOW) CYCLE

               W_G%FERTOT(N,NK,ISP)=(W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))
               IF (.NOT.LCHCON) THEN
                  W_G%FERTOT(N,NK,ISP)=(W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))/2
               ENDIF
            ENDDO
         ENDDO
         ENDDO
      ENDIF
! set selected components of the gradient to zero
! in order to test things
!test
!      W_G%CW=0;  W_G%CPROJ=0; DE1=0 ; ORT1=0

!      DO N1=1,WDES%NBANDS
!         CHAM(N1,N1,:,:)=0
!      ENDDO
!      W_G%FERTOT=0; DECEL=0; ORTCEL=0

!      IF (IO%IU0>=0) WRITE(IO%IU0,*)  'zero CHAM and CHF'
!      DO NK=2,WDES%NKPTS
!      DO N1=1,WDES%NB_TOT ; DO N2=N1+1,WDES%NB_TOT
!         CHAM(N1,N2,NK,2)=0 ; CHAM(N2,N1,NK,2)=0
!         CHF(N1,N2,NK,2)=0 ; CHF(N2,N1,NK,2)=0
!      ENDDO ; ENDDO
!      ENDDO
!      DE2=0; ORT2=0

      GNORM =DECEL  +DE1  +DE2
      ORT   =ORTCEL +ORT1 +ORT2

! calculate total energy for current set of orbitals
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC UPDATE SELF(W%CELTOT)
GPOMP !$OMP TARGET UPDATE FROM(W%CELTOT) IF(OFFLOAD_ON)
      E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
      TOTEN_LOCAL=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF+E%EDOTP+Ediel_sol+E%ESCPC
      IF (M_CONSTRAINED()) TOTEN_LOCAL=TOTEN_LOCAL+E_CONSTRAINT()

!=======================================================================
! perform CORRECTOR step taking into account the new gradient only
!=======================================================================
      l_is_updated= (IFLAG/=0 .AND. ((IALGO==9 .AND. MOD(ICOUNT,2)==0) .OR. &
                          (IALGO==8 .AND. TOTEN_LOCAL > TOTEN .AND. GSL/=0 .AND. INFO%ISEARCH>=0 )  .OR. &
                          (IALGO==8 .AND. ABS(ORT)>ABS(GSL)/5 .AND. GSL/=0 .AND. INFO%ISEARCH==1 )))
!                          (IALGO==8 .AND. ABS(ORT)>ABS(GSL)/5 .AND. GSL/=0 .AND. NCOUNT<=3 .AND. INFO%ISEARCH==1 )))
      CALLMPI( M_bcast_l(WDES%COMM, l_is_updated, 1))
 cor: IF ( l_is_updated) THEN
      ! GSL and ORT are the 1st order change of energy along previous search direction
      ! from start point and current point, respectively
      IF (IALGO==8) THEN
      ! for IALGO=8, this is a final backup; total step is BSTEP yet, correct this
        BCSTEP= BSTEP*ORT/(GSL-ORT)  ! harmonic result, correction to BSTEP
        BOPTH = BSTEP+BCSTEP
        BOPT=BOPTH
        ! integrate over gradient times path
        DESUM1=TOTEN_LOCAL-TOTEN-ORT/2*BCSTEP

        NCOUNT=NCOUNT+1
        STEP(NCOUNT) =BSTEP-BTRIAL/1E4_q
        ENRGY(NCOUNT)=TOTEN_LOCAL-TOTEN+ORT*BTRIAL/1E4_q ! energy and slope
        NCOUNT=NCOUNT+1
        STEP(NCOUNT) =BSTEP+BTRIAL/1E4_q
        ENRGY(NCOUNT)=TOTEN_LOCAL-TOTEN-ORT*BTRIAL/1E4_q ! energy and slope

        CALL search_minimum_polynomial_fit_order2( STEP(1:NCOUNT), ENRGY(1:NCOUNT), NCOUNT, 4, 0, BOPT4, DESUM4, DER, IFAIL, IO%IU0)
        IF (IFAIL/=0 .OR. DER<0 .OR.  BOPT4/=MAX(BSTEP/STEP_SCALING**3, MIN(BSTEP*STEP_SCALING, BOPT4))) &
        CALL search_minimum_polynomial_fit_order2( STEP(1:NCOUNT), ENRGY(1:NCOUNT), NCOUNT, 4, 1, BOPT4, DESUM4, DER, IFAIL, IO%IU0)
        ! use anharmonic result only when safeguard does not apply and direction of harmonic is equal unharmonic case
        IF (IFAIL==0 .AND. DER>=0 .AND. BOPT4==MIN(BSTEP*STEP_SCALING, BOPT4) .AND. &
            SIGN(1.0_q,BOPTH-BSTEP)==SIGN(1.0_q,BOPT4-BSTEP)) THEN
          BOPT=BOPT4
          DESUM1=DESUM4
        ENDIF

        BOPT=MAX(BSTEP/STEP_SCALING**3, MIN(BSTEP*STEP_SCALING, BOPT))
        BCSTEP=BOPT-BSTEP
        BSTEP =BOPT
        IF (IO%IU0>=0) WRITE(IO%IU0,12) BSTEP, BOPTH, BOPT4
   12 FORMAT('  final corrected         step=',F8.4,' harm=',F8.4,' 4th-ord=',F8.4)

        IF (IDUMP>=1 .AND.IO%IU0>=0 ) THEN
          WRITE(IO%IU0,'("  steps along line   ",100E14.5)') STEP(1:NCOUNT)
          WRITE(IO%IU0,'("  energies along line",100E14.5)') ENRGY(1:NCOUNT)
        ENDIF
      ELSE
        ! for IALGO=9, just trial step BTRIAL was used
        BCSTEP= BTRIAL*ORT/(GSL-ORT)
        ! make a small step (in total of 0.2 x BSTEP) if we overshoot
        IF (BTRIAL+BCSTEP<0) BCSTEP=-BTRIAL*0.8
        BSTEP = BTRIAL+BCSTEP
        DESUM1= -GSL*BSTEP/2
        IF (IO%IU0>=0) WRITE(IO%IU0,12) BSTEP, BSTEP
      ENDIF

      CALLMPI( M_bcast_d(WDES%COMM, BCSTEP, 1))
      CALLMPI( M_bcast_d(WDES%COMM, BSTEP, 1))

      ! update step size to gliding average (IALGO=8 does it below)
      IF (IALGO==9) AVERAG=AVERAG*0.9_q+BSTEP*0.1_q

      ICOUNT=ICOUNT+1

! step for wavefunction
      IF (IROT==1 .OR. IROT>=3) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO M=1,WDES%NRPLWV
            W%CW(M,N,NK,ISP)=W%CW(M,N,NK,ISP)   -W_F%CW(M,N,NK,ISP)      *BCSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO NPRO_=1,WDES%NPROD
            W%CPROJ(NPRO_,N,NK,ISP)=W%CPROJ(NPRO_,N,NK,ISP)-W_F%CPROJ(NPRO_,N,NK,ISP)*BCSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

! step for subspace-matrix
      IF (IROT==2 .OR. IROT==3 .OR. IROT==4) THEN
      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        CALL SETWDES(WDES,WDES1,NK)
#ifdef testrot
       CALL WGCOPY(NB_TOT*NB_TOT,CHF(1,1,NK,ISP),1,CEIDB(1,1),1)
#endif
        IF (ILOEW==1) THEN
           CALL ROT1(NB_TOT,CHF(1,1,NK,ISP),BCSTEP,CUNI,CEIDB,COVL)
        ELSE IF (ILOEW==2) THEN
           CALL ROT2(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),BCSTEP,CUNI,CEIDB,COVL)
        ENDIF

        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP),CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP))
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))

        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP),CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP))
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
      ENDDO
      ENDDO
      ELSE IF (IROT==6 .OR. IROT==7) THEN

      CALL MRG_AUX(WDES,W)

      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        CALL SETWDES(WDES,WDES1,NK)
        ! here we need to rotate CHF if the orbitals are rotated (as well as W%AUXTOT)
        CALL ROTETA(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),W_F%CELTOT(:,NK,ISP),BCSTEP,CUNI,CEIDB,COVL,.FALSE.,W%AUXTOT(:,NK,ISP))

        WA=ELEMENTS(W, WDES1, ISP)
        WA_F=ELEMENTS(W_F, WDES1, ISP)

        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             WA%CW_RED,WA%CPROJ_RED)
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))

        ! rotate W_F accordingly so that conjugation works out in next step (orbital change)
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',WA_F%CW_RED,WA_F%CPROJ_RED,CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             WA_F%CW_RED,WA_F%CPROJ_RED)
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
      ENDDO
      ENDDO
      ENDIF

      CALL ORTHCH(W%WDES,W, LOVERL, LMDIM,CQIJ)

! step for fermi-weights

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
        DO N =1,W%WDES%NB_TOTK(NK,ISP)
          W_F%CELTOT(N,NK,ISP)=W_F%CELTOT(N,NK,ISP)-BCSTEP*W_F%FERTOT(N,NK,ISP)
        ENDDO
        ENDDO
        ENDDO
      ENDIF

      IFLAG=0
      IRET=5
      RMS=GSL  ! report ORT and GSL in OSZICAR (as rms and ort)
      GSL=0    ! clear GSL to avoid second excecution
      GOTO 1000

 ENDIF cor
!=======================================================================
! performe trial step (possibly with conjugation of search directions)
! and calculate change in energy up to first order
! the conjugate gradient optimization is restarted
! ) after 30 steps
!=======================================================================
!  start resetting all counters
      NCOUNT=0
      NCOUNT_NRESET=0

      IF ( IALGO==3 .AND. ICOUNT ==-1) THEN
         ! ideally I would like to use the following term
         ! for integrating the equation of motion
         ! VEL=((1-FRICTION/2)*VEL+2*F*TIMESTEP)/(1+FRICTION/2)
         ! but this is not what is implemented
         ! hence we rescale the timestep (BTRIAL) and set GAMMA properly
         !
         ! FRICTION is proportional to sqrt(min eigenvalue(P H) / max eigenvalue(P H))
         ! The maximum eigenvalue is proportional to the inverse step size.
         ! The minimum eigenvalue is 1 for insulators and around 2 for metals
         ! since the steps on the eigenvalues are damped by 2.
         ! some overrelaxation also never harms
         ! For non collinear calculations, however very soft modes exist
         ! and one might need to "over-relax" e.g. FRICTION=SQRT(BTRIAL*100)
         IF (ICEL==0) THEN
            FRICTION=SQRT(BTRIAL)
         ELSE
            FRICTION=SQRT(BTRIAL*2)
         ENDIF

         ! in principle BTRIAL needs to be rescaled
         ! but empirically that does not improve convergence, hence stick to the original timestep
         ! BTRIAL=BTRIAL/(1+FRICTION/2)*2
      ENDIF

      IF ( IALGO==3 .AND. IFLAG==1 .AND. ICOUNT /=-1) THEN
         GAMMA=(1-FRICTION/2)/(1+FRICTION/2)
         IRET    =3
      ELSE IF ( IALGO==4 .AND. IFLAG==1 .AND. ICOUNT /=-1) THEN
         GAMMA=1
         ! DESUM0 expected first order energy change
         ! when larger 0 we are starting to move uphill
         ! hence the velocities are quenched to zero
         DESUM0=GAMMA*ORT
         ICOUNT=ICOUNT+1
         IF (DESUM0<0) THEN
            ICOUNT=0
            GAMMA=0
         ENDIF
         IRET  =4
      ! reset every 40 steps
      ELSE IF ( (IALGO <8) .OR.  IFLAG/=1 .OR. ICOUNT == 40 &
     &   ) THEN
        GAMMA =0
        ICOUNT=0
        IRET  =0
      ELSE
! Fletcher-Reeves
        GAMMA=GNORM /GNORML
        ICOUNT=ICOUNT+1
        IRET  =1
        ! previous line minization was so bad that ORT*GAMMA is larger
        ! then it is useless to conjugate
        ! (ORT should be essentially zero)
        IF ( ABS(GAMMA*ORT) > ABS(GNORM)/4 .OR. ABS(GAMMA)>30 ) THEN
          GAMMA=0
          ICOUNT=0
          IRET  =0
        ENDIF
      ENDIF
!test
!      GAMMA=0
      GNORML=GNORM
      RMS   =GNORM
      CALLMPI( M_bcast_d(WDES%COMM, GAMMA, 1))
      CALLMPI( M_bcast_i(WDES%COMM, ICOUNT, 1))
      CALLMPI( M_bcast_i(WDES%COMM, IRET, 1))

   20 FORMAT(' gam=',F6.3,' g(H,U,f)= ',3E10.3,' ort(H,U,f) =',3E10.3)
!-----------------------------------------------------------------------
!  DESUM0 is expected first order energy-change into the trial direction
!  the trial direction is  f(N) = p(N) + GAMMA f(N-1)
!  (g(N) gradient, p(N) preconditioned gradient, f(N) search direction)
!  because ORT= g(N) f(N-1) and GNORM = g(N) p(N)
!  DESUM is given by DESUM = g(N) * f(N) = GNORM+GAMMA* ORT
!-----------------------------------------------------------------------
      DESUM0=GNORM+GAMMA*ORT
      GSL   =GNORM+GAMMA*ORT  ! store this for later use
!-----------------------------------------------------------------------
! conjugate directions
!-----------------------------------------------------------------------
      IF (IROT==1 .OR. IROT>=3) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W_F,W_G) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W_F,W_G) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO M=1,WDES%NRPLWV
           W_F%CW(M,N,NK,ISP)      =W_G%CW(M,N,NK,ISP)    +W_F%CW(M,N,NK,ISP)*GAMMA
        ENDDO
        ENDDO
        ENDDO
        ENDDO

GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W_F,W_G) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W_F,W_G) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO NPRO_=1,WDES%NPROD
          W_F%CPROJ(NPRO_,N,NK,ISP)=W_G%CPROJ(NPRO_,N,NK,ISP)+W_F%CPROJ(NPRO_,N,NK,ISP)*GAMMA
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

      IF (IROT==2 .OR. IROT>=3) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,CHF,CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,CHF,CHAM) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N2=1,NB_TOT
        DO N1=1,NB_TOT
          CHF(N1,N2,NK,ISP)=CHAM(N1,N2,NK,ISP)+GAMMA*CHF(N1,N2,NK,ISP)
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(W,W_F,W_G) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:W,W_F,W_G) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
        DO N =1,W%WDES%NB_TOTK(NK,ISP)
          W_F%FERTOT(N,NK,ISP)=W_G%FERTOT(N,NK,ISP)+W_F%FERTOT(N,NK,ISP)*GAMMA
        ENDDO
        ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! trial step
!-----------------------------------------------------------------------
      BSTEP=BTRIAL
      IF (IROT==1 .OR. IROT>=3) THEN
      IF (restore_W) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W_G,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W_G,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO M=1,WDES%NRPLWV
             W_G%CW(M,N,NK,ISP)    =W%CW(M,N,NK,ISP)
             W%CW(M,N,NK,ISP)=W%CW(M,N,NK,ISP)-W_F%CW(M,N,NK,ISP)*BSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W_G,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W_G,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO NPRO_=1,WDES%NPROD
            W_G%CPROJ(NPRO_,N,NK,ISP)=W%CPROJ(NPRO_,N,NK,ISP)
            W%CPROJ(NPRO_,N,NK,ISP)=W%CPROJ(NPRO_,N,NK,ISP)-W_F%CPROJ(NPRO_,N,NK,ISP)*BSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

      ELSE
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO M=1,WDES%NRPLWV
             W%CW(M,N,NK,ISP)=W%CW(M,N,NK,ISP)-W_F%CW(M,N,NK,ISP)*BSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO NPRO_=1,WDES%NPROD
            W%CPROJ(NPRO_,N,NK,ISP)=W%CPROJ(NPRO_,N,NK,ISP)-W_F%CPROJ(NPRO_,N,NK,ISP)*BSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

      ENDIF
      ENDIF

! trial step for subspace-matrix

      IF (IROT==2 .OR. IROT==3 .OR. IROT==4) THEN
      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        CALL SETWDES(WDES,WDES1,NK)
#ifdef testrot
       CALL WGCOPY(NB_TOT*NB_TOT,CHF(1,1,NK,ISP),1,CEIDB(1,1),1)
#endif
        IF (ILOEW==1) THEN
           CALL ROT1(NB_TOT,CHF(1,1,NK,ISP),BSTEP,CUNI,CEIDB,COVL)
        ELSE IF (ILOEW==2) THEN
           CALL ROT2(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),BSTEP,CUNI,CEIDB,COVL)
        ENDIF

        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP),CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP))
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))

        IF (.NOT. restore_W) THEN
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
          CALL LINCOM('F',W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP),CEIDB(1,1), &
              W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
              WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
              W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP))
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
        ENDIF
      ENDDO
      ENDDO
      ELSE IF (IROT==6 .OR. IROT==7) THEN

      CALL MRG_AUX(WDES,W)
      IF (restore_W) THEN
GPACC !$ACC KERNELS PRESENT(W_F,W_G)  __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W_F,W_G) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        W_G%CELTOT=W_F%CELTOT
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ENDIF

      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        CALL SETWDES(WDES,WDES1,NK)

        CALL ROTETA(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),W_F%CELTOT(:,NK,ISP),BSTEP,CUNI,CEIDB,COVL,restore_W,W%AUXTOT(:,NK,ISP))

        WA=ELEMENTS(W, WDES1, ISP)
        WA_F=ELEMENTS(W_F, WDES1, ISP)

         IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
         CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,CEIDB(1,1), &
              W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
              WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
              WA%CW_RED,WA%CPROJ_RED)
         IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))

        IF (.NOT. restore_W) THEN
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
          CALL LINCOM('F',WA_F%CW_RED,WA_F%CPROJ_RED,CEIDB(1,1), &
              W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
              WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
              WA_F%CW_RED,WA_F%CPROJ_RED)
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
        ENDIF
      ENDDO
      ENDDO
      ENDIF
      CALL ORTHCH(W%WDES,W, LOVERL, LMDIM,CQIJ)
!
! fermi-weights
! update fermi-weights W%FERTOT

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:W,W_F) IF(OFFLOAD_ON)
         DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
         DO N =1,W%WDES%NB_TOTK(NK,ISP)
            W_F%CELTOT(N,NK,ISP)=W_F%CELTOT(N,NK,ISP)-BSTEP*W_F%FERTOT(N,NK,ISP)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
!-----------------------------------------------------------------------
! set IFLAG and return to main routine
! DESUM1 is the first order energy change
!-----------------------------------------------------------------------
      IF (IALGO<=6 .OR. IALGO==9) THEN
        IFLAG=0
      ELSE
        IFLAG=1
      ENDIF
      IF (IFLAG==0) DESUM1 =-DESUM0*BSTEP
      GOTO 1000
!***********************************************************************
! CORRECTOR STEP
! IFLAG=2  goto the minimum
! ) calculate the remaining distance we have to go
! ) and step to the minimum
! IFLAG>2  increase trial step
!***********************************************************************
   ELSE IF (IFLAG>=2) THEN main
#ifdef MPI
      IF ((LHFCALC.OR.LUSEPEAD()).AND.WDES%COMM_KINTER%NCPU.GT.1) THEN
         CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
      ENDIF
#endif
      ! redistribute data first
      IF (NCPU /= 1) THEN

         DO ISP=1,WDES%ISPIN
            DO NK=MYKSTART,WDES%NKPTS,MYKINCR
               CALL SETWDES(WDES,WDES1,NK)
               CALL REDIS_PW(WDES1, WDES%NBANDS, W%CW(1,1,NK,ISP))
            ENDDO
         ENDDO

      W%OVER_BAND=.TRUE.
      ENDIF

! store step and energy for later use
      IF (NCOUNT==0) THEN
        STEP(1) =-BSTEP/1E4
        ENRGY(1)=-DESUM0*STEP(1)  ! initial energy
        STEP(2) =BSTEP/1E4
        ENRGY(2)=-DESUM0*STEP(2) ! add one more value corresponding to slope
        NCOUNT=2
        NCOUNT_BASE=3
        NCOUNT_NRESET=1
        LBRACKET=.FALSE.
        LSLOT=.FALSE.
        BOPT=BSTEP   ! avoid side effects on BSTEP
        CALL SLOT_IN(BOPT, BTRIAL, STEP_SCALING, 1, LSLOT, L_IS_UPDATED, .FALSE., IO%IU0) ! set first slot
      ENDIF

      NCOUNT=NCOUNT+1
      IF (NCOUNT+3 > NSTEP) CALL vtutor%bug("NCOUNT exceeds array dimensions", __FILE__, __LINE__)
      STEP(NCOUNT) =BSTEP
      ENRGY(NCOUNT)=TOTEN2-TOTEN

! calculate position of minimum using a quadratic interpolation
      A1= TOTEN2-TOTEN
      A2=-DESUM0*BSTEP

      IF (IDUMP>=2 .AND. IO%IU0>=0) WRITE(IO%IU0,11,ADVANCE='NO') BSTEP,A1,A2
   11 FORMAT('step=',F8.4,' , ',F13.8,F13.8)

      BOPT  = -A2/(A1-A2)/2
      DESUM1= (A2+(A1-A2)*BOPT)*BOPT
!-----------------------------------------------------------------------
! IFLAG==2
!-----------------------------------------------------------------------
      IF (IFLAG==2) THEN
!---- first trial step, search near minimum
        IF (NCOUNT==3) THEN
           BOPTH0=BOPT*BSTEP
           IF (IDUMP>=1.AND. IO%IU0>=0) WRITE(IO%IU0,20) GAMMA,DE1,DE2,DECEL,ORT1,ORT2,ORTCEL
           ! yields same results as above; however, numerics are superior above
           ! CALL search_minimum_polynomial_fit( STEP(1:3), ENRGY(1:3), 3, BOPT, DESUM1, DER, IO%IU0) ; BOPT=BOPT/BTRIAL
           IF (IDUMP>=1 .AND.IO%IU0>=0 ) WRITE(IO%IU0,10) GAMMA,BTRIAL,BOPTH0,AVERAG
   10 FORMAT(' gam=',F6.3,' trial=',F6.3,'  step=',F8.4,' mean=',2F8.4)

! INTERATE_PREC adds more energy calculations
!    and seeks the minimum from these additional energies
! apply this if BOPT gets unreasonably large (indicating
! non-quadratic behaviour or some other trouble) or negative
! the criteria are obviously very heuristic; in order to avoid energy increase
           IF (.NOT. ITERATE_PREC .AND. BOPT<STEP_SCALING**3 .AND. BOPT>1.0_q/STEP_SCALING**3 ) THEN
              IFLAG=0
              BCSTEP=BSTEP*BOPT-BSTEP  ! BSTEP and BTRIAL are identical here
              BSTEP=BSTEP*BOPT
           ELSE
              ! negative step size usually indicates too small trial step
              IF (BOPT<0) BOPT=STEP_SCALING*STEP_SCALING
              BOPT=BSTEP*BOPT  ! convert to actual length
              CALL SLOT_IN(BOPT, BTRIAL, STEP_SCALING, 3, LSLOT, L_IS_UPDATED, .FALSE., IO%IU0)
              STEP(NCOUNT+1)=BOPT
              BCSTEP=BOPT-BSTEP
              BSTEP =BOPT
           ENDIF
!---- subsequent steps
        ELSE
   14 FORMAT('  continued last =',F6.3,'  step=',F8.4,' harm=',F8.4,' 4th-ord=',F8.4)
   15 FORMAT('  continued last =',F6.3,'  step=',F8.4,' harm=',F8.4,' 4th-ord=',F8.4,' spline=',F8.4)
   16 FORMAT('  concave          backup-step=',F8.4)

          ! harmonic fit through three points
          CALL search_minimum_polynomial_fit_order2( STEP(1:NCOUNT), ENRGY(1:NCOUNT), NCOUNT, 2, 0, BOPTH, DESUM1, DER, IFAIL, -1)
          BOPT=BOPTH
          ! 4th order fit
          CALL search_minimum_polynomial_fit_order2( STEP(1:NCOUNT), ENRGY(1:NCOUNT), NCOUNT, 4, 0, BOPT4, DESUM4, DER4, IFAIL, IO%IU0)
          IF (IFAIL==0) THEN ! use 4th order fit for next move
            BOPT=BOPT4
            DESUM1=DESUM4
            DER=DER4
          ENDIF
          IF (NCOUNT>5) THEN ! spline seems most robust
            CALL search_minimum_spline( STEP(1:NCOUNT), ENRGY(1:NCOUNT), NCOUNT, BOPTS, DESUMS, DERS, IFAIL, IO%IU0)
            IF (IFAIL==0 .AND. DERS>=0) THEN
              BOPT=BOPTS
              DESUM1=DESUMS
              DER=DERS
            ENDIF
            IF (IO%IU0>=0) WRITE(IO%IU0,15) BSTEP,BOPT,BOPTH,BOPT4,BOPTS
          ELSE
            IF (IO%IU0>=0) WRITE(IO%IU0,14) BSTEP,BOPT,BOPTH,BOPT4
          ENDIF
          L_IS_UPDATED=.FALSE.
          ! max 10 steps and if the function is obviously harmonic, do not continue
          IF (NCOUNT<=10 .AND. &
            ! not strongly unharmonic
             (2*ABS(BOPT-BOPTH)/ABS(ABS(BOPT)+ABS(BOPTH))>1E-1 .OR. 2*ABS(BOPT-BOPTH0)/ABS(ABS(BOPT)+ABS(BOPTH0))>1E-1) ) THEN
            ! determine available slot closeby, if harmonic and anharmonic minima do not agree
            CALL SLOT_IN(BOPT, BTRIAL, STEP_SCALING, 2, LSLOT, L_IS_UPDATED, .FALSE., IO%IU0)
            ! also stop if slots are exhausted
          ENDIF
          STEP(NCOUNT+1)=BOPT
          BCSTEP=BOPT-BSTEP
          BSTEP =BOPT
          IF (.NOT. L_IS_UPDATED) THEN
            IFLAG=0 ! terminate search
            ! minimum found
            IF (IDUMP>=1 .AND.IO%IU0>=0 ) THEN
              WRITE(IO%IU0,'("  steps along line   ",100E14.5)') STEP(1:NCOUNT)
              WRITE(IO%IU0,'("  energies along line",100E14.5)') ENRGY(1:NCOUNT)
            ENDIF
            IF ( DER<0 ) THEN  ! second derivative smaller zero
            ! concave; search for minimal energy yet found in the list and
            ! terminate
              INDEX=MINLOC(ENRGY(1:NCOUNT))
              BOPT  =STEP(INDEX(1))
              DESUM1=ENRGY(INDEX(1))
              BCSTEP=BOPT-BSTEP
              BSTEP =BOPT
              IF (IO%IU0>=0) WRITE(IO%IU0,16) BSTEP
            ENDIF
          ENDIF
        ENDIF
!-----------------------------------------------------------------------
! IFLAG>2
!-----------------------------------------------------------------------
      ELSE
!---- increase trial step by BTRIAL
        BCSTEP=BTRIAL
        BSTEP =BSTEP+BTRIAL
      ENDIF

      CALLMPI( M_bcast_d(WDES%COMM, BCSTEP, 1))
      CALLMPI( M_bcast_d(WDES%COMM, BSTEP, 1))
      CALLMPI( M_bcast_i(WDES%COMM, IFLAG, 1))

       ! gK 07.01.2025 average only when the minimum has been found
      IF (IFLAG==0) AVERAG=AVERAG*0.9_q+BSTEP*0.1_q

      NK=1
!-----------------------------------------------------------------------
! goto minimum (or continue trial step)
!-----------------------------------------------------------------------
      IF (IROT==1 .OR. IROT>=3) THEN
      IF (restore_W) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_G,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_G,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO M=1,WDES%NRPLWV
          W%CW(M,N,NK,ISP)=W_G%CW(M,N,NK,ISP) -W_F%CW(M,N,NK,ISP)      *BSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_G,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_G,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO NPRO_=1,WDES%NPROD
          W%CPROJ(NPRO_,N,NK,ISP)=W_G%CPROJ(NPRO_,N,NK,ISP)-W_F%CPROJ(NPRO_,N,NK,ISP)*BSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

      ELSE
GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO M=1,WDES%NRPLWV
          W%CW(M,N,NK,ISP)=W%CW(M,N,NK,ISP)   -W_F%CW(M,N,NK,ISP)      *BCSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

GPACC !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        DO N=1,WDES%NBANDS
        DO NPRO_=1,WDES%NPROD
          W%CPROJ(NPRO_,N,NK,ISP)=W%CPROJ(NPRO_,N,NK,ISP)-W_F%CPROJ(NPRO_,N,NK,ISP)*BCSTEP
        ENDDO
        ENDDO
        ENDDO
        ENDDO

      ENDIF
      ENDIF
! step for subpace-matrix

      IF (IROT==2 .OR. IROT==3 .OR. IROT==4) THEN
      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        CALL SETWDES(WDES,WDES1,NK)
#ifdef testrot
        CALL WGCOPY(NB_TOT*NB_TOT,CHF(1,1,NK,ISP),1,CEIDB(1,1),1)
#endif
        IF (ILOEW==1) THEN
         IF (restore_W) THEN
           CALL ROT1(NB_TOT,CHF(1,1,NK,ISP),BSTEP,CUNI,CEIDB,COVL)
         ELSE
           CALL ROT1(NB_TOT,CHF(1,1,NK,ISP),BCSTEP,CUNI,CEIDB,COVL)
         ENDIF
        ELSE IF (ILOEW==2) THEN
         IF (restore_W) THEN
           CALL ROT2(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),BSTEP,CUNI,CEIDB,COVL)
         ELSE
           CALL ROT2(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),BCSTEP,CUNI,CEIDB,COVL)
         ENDIF
        ENDIF

        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP),CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             W%CW(:,:,NK,ISP),W%CPROJ(:,:,NK,ISP))
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))

        IF (.NOT. restore_W) THEN
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
          CALL LINCOM('F',W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP),CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             W_F%CW(:,:,NK,ISP),W_F%CPROJ(:,:,NK,ISP))
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
        ENDIF
      ENDDO
      ENDDO

      ELSE IF (IROT==6 .OR. IROT==7) THEN

      CALL MRG_AUX(WDES,W)
      IF (restore_W) THEN
GPACC !$ACC KERNELS PRESENT(W_G,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:W_G,W_F) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
        W_F%CELTOT=W_G%CELTOT ! restore original Fermi-weights and apply full
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ENDIF

      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
        CALL SETWDES(WDES,WDES1,NK)
         ! CHF   on entry: search direction for the pseudo Hamiltonian
         !       on exit:  search direction for the pseudo Hamiltonian U+ CHF U
         ! BSTEP step size
         ! CEIDB unitary rotation matrix U
         ! W_F%CELTOT on entry: previous pseudo Hamiltonian (diagonal eigenvalues)
         !            on exit:  new diagonal components of pseudo Hamiltonian
        IF (restore_W) THEN
          IF (IFLAG==0) THEN
            CALL ROTETA(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),W_F%CELTOT(:,NK,ISP),BSTEP,CUNI,CEIDB,COVL,.FALSE.,W%AUXTOT(:,NK,ISP))
          ELSE ! terminate then rotate CHF now
            CALL ROTETA(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),W_F%CELTOT(:,NK,ISP),BSTEP,CUNI,CEIDB,COVL,restore_W,W%AUXTOT(:,NK,ISP))
          ENDIF
        ELSE
          CALL ROTETA(WDES%COMM_KIN,NB_TOT,CHF(1,1,NK,ISP),W_F%CELTOT(:,NK,ISP),BCSTEP,CUNI,CEIDB,COVL,restore_W,W%AUXTOT(:,NK,ISP))
        ENDIF

        WA=ELEMENTS(W, WDES1, ISP)
        WA_F=ELEMENTS(W_F, WDES1, ISP)

        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             WA%CW_RED,WA%CPROJ_RED)
        IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))

        IF (.NOT. restore_W .OR. IFLAG==0) THEN
          ! rotate search direction in the same way as orbitals (always on exit)
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
          CALL LINCOM('F',WA_F%CW_RED,WA_F%CPROJ_RED,CEIDB(1,1), &
             W%WDES%NB_TOTK(NK,ISP),W%WDES%NB_TOTK(NK,ISP), &
             WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
             WA_F%CW_RED,WA_F%CPROJ_RED)
          IF (WDES%DO_REDIS) CALL REDIS_PROJ(WDES1, NBANDS, W_F%CPROJ(1,1,NK,ISP))
        ENDIF
      ENDDO
      ENDDO
      ENDIF
      CALL ORTHCH(W%WDES,W, LOVERL, LMDIM,CQIJ)
! step for fermi-weights

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,W,W_F) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:WDES,W,W_F) IF(OFFLOAD_ON)
        DO ISP=1,WDES%ISPIN
        DO NK=MYKSTART,WDES%NKPTS,MYKINCR
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
        DO N =1,W%WDES%NB_TOTK(NK,ISP)
            W_F%CELTOT(N,NK,ISP)=W_F%CELTOT(N,NK,ISP)-BCSTEP*W_F%FERTOT(N,NK,ISP)
        ENDDO
        ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
   ENDIF main

1000 CONTINUE
   CALL REDIS_PW_OVER_BANDS(WDES, W)

#ifdef MPI
   IF ((LHFCALC.OR.LUSEPEAD()).AND.WDES%COMM_KINTER%NCPU.GT.1) THEN
      CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
      CALL MRG_FERWE(WDES,W)
   ENDIF
#endif

GPACC !$ACC EXIT DATA DELETE(COVL,CUNI,CEIDB) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:COVL,CUNI,CEIDB) IF(OFFLOAD_ON)
   DEALLOCATE(COVL,CUNI,CEIDB)

   IF (IFLAG_IN<2) THEN
      DO I=1,NSIM_LOCAL
         CALL DELWAV(W1(I), .TRUE.)
      ENDDO
GPACC !$ACC EXIT DATA DELETE(W1(:)) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1) IF(OFFLOAD_ON)
      CALL DELWAV(WSEARCH, .TRUE.)
GPACC !$ACC EXIT DATA DELETE(WSEARCH) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(CELTOT) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:WSEARCH) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CELTOT) IF(OFFLOAD_ON)

GPACC !$ACC UPDATE SELF(W_G%CW,W_G%CPROJ,W_G%FERTOT) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE FROM(W_G%CW,W_G%CPROJ,W_G%FERTOT) IF(OFFLOAD_ON)
   ENDIF
#ifdef OFFLOADING
    CALL DELETE_TYPED_VAR(WDES1)

GPACC !$ACC UPDATE SELF(W%CW,W%CPROJ,W%CELTOT,W%FERTOT,W%AUXTOT) __IF_ASYNC__
GPACC !$ACC UPDATE SELF(W_F%CW,W_F%CPROJ,W_F%CELTOT,W_F%FERTOT) __IF_ASYNC__

GPOMP !$OMP TARGET UPDATE FROM(W%CW,W%CPROJ,W%CELTOT,W%FERTOT,W%AUXTOT) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE FROM(W_F%CW,W_F%CPROJ,W_F%CELTOT,W_F%FERTOT) IF(OFFLOAD_ON)

   CALL DELETE_TYPED_VAR(P)
   CALL DELETE_TYPED_VAR(GRIDC)
   CALL DELETE_TYPED_VAR(GRID)
   CALL DELETE_TYPED_VAR(NONLR_S)
   CALL DELETE_TYPED_VAR(NONL_S)
   CALL DELETE_TYPED_VAR(W_G)
   CALL DELETE_TYPED_VAR(W_F)
   CALL DELETE_TYPED_VAR(W)
   CALL DELETE_TYPED_VAR(WDES)

GPACC !$ACC EXIT DATA DELETE(DENCOR,CSTRF) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:DENCOR,CSTRF) IF(OFFLOAD_ON)

GPACC !$ACC EXIT DATA DELETE(HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC EXIT DATA DELETE(HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.OFFLOAD_ON)

GPACC !$ACC EXIT DATA DELETE(LATT_CUR,SV,CDIJ,CQIJ) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:LATT_CUR,SV,CDIJ,CQIJ) IF(OFFLOAD_ON)

GPACC !$ACC EXIT DATA COPYOUT(CHAM,CHF) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:CHAM,CHF) IF(OFFLOAD_ON)

   CALL OFFLOAD_ON_POP
#endif

   PROFILING_STOP('edwav')

   RETURN
 END SUBROUTINE EDWAV

!************************ SUBROUTINE ROTDIA ****************************
!
!> This subroutine calculates the rotation matrix  for one k-point.
!> The commutator CHF=[F,H] must be supplied by the calling routine.
!> on return the rotation matrix is stored in   CEIDB.
!> The rotation matrix is
!>~~~
!>    -i d [F,H]       -i d E +                     +
!>   e            = U e      U    with [F,H] = U E U
!>~~~
!> i.e. U is the matrix which diagonolizes [F,H] and E is the diagonal
!>      matrix containing the eigenvalues. \n
!>  Routine is no longer used.
!***********************************************************************

      SUBROUTINE ROTDIA(NBANDS,NK,CHF,BSTEP,W,NW,CUNI,CEIDB,CTMP)
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CHF(NBANDS,NBANDS)
      DIMENSION CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &          CEIDB(NBANDS,NBANDS)

      DIMENSION HFEIG(NBANDS),W(NW)

!=======================================================================
! Diagononalize the resulting Matrix (search direction)
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
        CUNI(N1,N2)=CHF(N1,N2)
      ENDDO
      ENDDO

      CALL ZHEEV('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
     &         W,  IFAIL)
      IF (IFAIL/=0) THEN
         CALL vtutor%error("ERROR ROTDIA: Call to routine ZHEEV failed! Error code was " // str(IFAIL))
      ENDIF
!=======================================================================
! set up the unitary transformation CEIDB for wavefunctions
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CUNI(N1,N2)*EXP(-BSTEP*(0._q,1._q)*HFEIG(N2))
      ENDDO
      ENDDO
      CALL ZGEMM( 'N', 'C', NBANDS, NBANDS, NBANDS, (1.,0.), CTMP, &
     &             NBANDS, CUNI, NBANDS, (0._q,0._q), CEIDB, NBANDS)
!=======================================================================
! rotate the wavefunction (with CEIDB): warning about LINCOM! We must
! take the transposed of CEIDB because MATMUL uses the transposed ...
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CEIDB(N2,N1)
      ENDDO
      ENDDO

      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CEIDB(N2,N1)=CTMP(N2,N1)
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE


!************************ SUBROUTINE ROTINI    *************************
!
!> This subroutine calculates a preconditioned rotation matrix CHF.
!> The unconditioned rotation matrix is given by
!>~~~
!>   H(i,j) <- H(i,j) / epsilon(i)-epsilon(j)
!>~~~
!>   W               current wavefunctions and eigenvalues \n
!>   W_F%CELEN       from W_F%CELEN the occupancies are calculated \n
!>   CHAM on entry:  Hamiltonian between all states \n
!>        on exit:   first order rotation matrix \n
!>
!>   CHF  on entry:  previous rotation matrix \n
!>        on exit:   previous rotation matrix with some elements possibly
!>                   zero (no conjugation of those elements) \n
!>   CHAMP on entry: preconditioned Hamilton matrix \n
!>   ORT2            orthogonality of current gradient to previous rotation
!>                   matrix \n
!>   DE2             expected energy change along present rotation matrix \n
!>   WEIMIN          threshhold for empty bands \n
!>   CEIDB           work array \n
!>   SCALE           scale step length by supplied factor.
!>                   If CHAMP is supplied SCALE is required \n
!
!***********************************************************************

      SUBROUTINE ROTINI(WDES, W, W_F, CHAM, CHF, CEIDB, WEIMIN, ORT2, DE2,  &
           LCHCON, CHAMP, SCALE)

__GPU USE moffload_struct_def

      USE prec
      USE wave
      USE constant
      USE poscar
      USE pseudo
      USE mgrid
      USE lattice
      USE nonl_high
      USE augfast
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W, W_F

      GDEF CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      GDEF CHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      GDEF CEIDB(WDES%NB_TOT,WDES%NB_TOT)
#ifndef IFXWAR9
      GDEF, OPTIONAL :: CHAMP(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
#else
      GDEF, OPTIONAL :: CHAMP(:,:,:,:)
#endif
      REAL(q), OPTIONAL :: SCALE

      REAL(q) WEIMIN
      REAL(q) ORT2, DE2

      LOGICAL        :: LCHCON   ! charge updated
! local
      INTEGER ISP, NK, N1, N2
      REAL(q) WEIGHT, DIFFER, DIFCEL, THETA
      REAL(q), PARAMETER ::  DIFMAX=1E-8   ! threshhold for degeneragy
      GDEF CROT, CEXPRE
      REAL(q) :: EH
      REAL(q) :: ORT2_in,DE2_in

      ORT2_in=ORT2
      DE2_in=DE2
      ORT2=0
      DE2=0

!=======================================================================
! calculate charge of each band
!=======================================================================
      IF (.NOT. LCHCON) THEN
         EH=0.35
      ELSE
         EH=0
      ENDIF
!=======================================================================
! calculate off symmetry elements of rotation matrix
!=======================================================================
      DO ISP=1,WDES%ISPIN
      DO NK=1,WDES%NKPTS
GPACC !$ACC KERNELS PRESENT(CEIDB) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CEIDB) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         CEIDB=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         WEIGHT =WDES%WTKPT(NK)
         
GPACC !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) REDUCTION(+:DE2,ORT2) PRESENT(W,WDES,CHAM,CHAMP,W_F) &
GPACC !$ACC PRIVATE(DIFFER,CEXPRE,DIFCEL,CROT,THETA) __IF_ASYNC__

! firstprivate default mapping of SCALE fails for optional parameters if not present
! optional arrays should be not mapped with present but map(to:), overhead is minimal
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) REDUCTION(+:DE2,ORT2) &
#ifdef OFFWAR_OPT_PARAM
GPOMP !$OMP MAP(PRESENT,TO:W,WDES,W_F,CHAM,CEIDB,CHF) MAP(TO:CHAMP,SCALE) PRIVATE(DIFFER,CEXPRE,DIFCEL,CROT,THETA) IF(OFFLOAD_ON)
#else
GPOMP !$OMP MAP(PRESENT,TO:W,WDES,W_F,CHAM,CEIDB,CHF) MAP(TO:CHAMP) PRIVATE(DIFFER,CEXPRE,DIFCEL,CROT,THETA) IF(OFFLOAD_ON)
#endif

         DO N2=1,WDES%NB_TOT
            DO N1=1,WDES%NB_TOT
              DIFFER=     (W%FERTOT(N1,NK,ISP)-W%FERTOT(N2,NK,ISP))
              ! gradient
              CEXPRE=WDES%RSPIN*CHAM(N2,N1,NK,ISP)*DIFFER
              IF (PRESENT(CHAMP)) THEN
                DIFCEL= REAL( CHAMP(N2,N2,NK,ISP)-CHAMP(N1,N1,NK,ISP) ,KIND=q)
                CROT  = CHAMP(N1,N2,NK,ISP)*SCALE
              ELSE
                DIFCEL= REAL( CHAM(N2,N2,NK,ISP)-CHAM(N1,N1,NK,ISP) ,KIND=q)
                CROT  = CHAM(N1,N2,NK,ISP)
              ENDIF
              ! avoid rotation for degenerated pairs
              IF ( (ABS(DIFCEL)<DIFMAX)  &
              ! no rotation if occupancies are incorrect (very important !!!)
              ! this can lead to components with an energy increase along the search direction
              ! which the CG algorithm does not like
                .OR. (DIFCEL* REAL( W_F%CELTOT(N2,NK,ISP)-W_F%CELTOT(N1,NK,ISP) ,KIND=q) <DIFMAX)  &
              ) THEN
                CROT  =0
              ELSE
                ! EH is used to calculate a diagonal preconditioning matrix
                ! EH is propotional to 4 times the second derivative
                ! of Hartree energy for rotation by an angle THETA
                ! for fixed potential the second derivative is 2*DIFCEL
                IF (DIFCEL>0) THEN
                   DIFCEL=DIFCEL+2*EH*ABS(2*W%FERTOT(N1,NK,ISP)-2*W%FERTOT(N2,NK,ISP))
                ELSE
                   DIFCEL=DIFCEL-2*EH*ABS(2*W%FERTOT(N1,NK,ISP)-2*W%FERTOT(N2,NK,ISP))
                ENDIF
                ! Loewdin perturbation theory
                ! the position of the minimum is given  by first deriv./ second deriv.
                CROT = CROT/DIFCEL
                ! smooth cutoff function for rotation matrix
                THETA=ATAN(ABS(CROT*4))/4
                CROT=SIN(THETA)*CROT/MAX(ABS(CROT),1E-8_q)
              ENDIF
              ! avoid conjugation of those elements which do not contribute to energy
              ! these are all elements that have similar fermi-weights
              ! since rotation in that sub-space leaves the energy invariant
              IF (ABS(DIFFER)<WEIMIN) THEN
                 CHF(N1,N2,NK,ISP)=0
              ENDIF

              ! skip rows and columns below NBANDSLOW
              IF (N1<WDES%NBANDSLOW.OR.N2<WDES%NBANDSLOW) CYCLE

              CEIDB(N1,N2)=CROT
              ! gradient time previous search direction
              ORT2  =ORT2  +CEXPRE*CHF(N1,N2,NK,ISP)*WEIGHT
              ! gradient times present precond. gradient
              DE2   =DE2   +CEXPRE*CROT*WEIGHT
            ENDDO
         ENDDO

         CALL WGCOPY(WDES%NB_TOT*WDES%NB_TOT,CEIDB(1,1),1,CHAM(1,1,NK,ISP),1)
! for debugging purposes
!GPOMP !$OMP TARGET UPDATE FROM(CHAM)
!         CALL DUMP_HAM( "rotation matrix",WDES, CHAM(1,1,NK,ISP))

      ENDDO
      ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for ORT2 and DE2
      CALLMPI (M_sum_d(WDES%COMM_KINTER,ORT2,1) )
      CALLMPI (M_sum_d(WDES%COMM_KINTER,DE2,1) )

      ORT2=ORT2+ORT2_in
      DE2=DE2+DE2_in

      END SUBROUTINE ROTINI

!************************ SUBROUTINE ROT1  ****************************
!
!> This subroutine calculates a "Jacobi like" rotation matrix
!>~~~
!>  CEIDB = 1 - CHF* BSTEP
!>~~~
!> for one k-point, where CHF is the Loewdin-rotation matrix
!> which must be supplied by the calling routine
!>  (i.e. H(n1,n2) /  (H(n1,n1)-H(n2,n2)) \n
!> For a larger step CHF is not unitary.
!> A better approximation is obtained by diagonalisation of each 2x2 subblock
!>~~~
!>  H (n1, n1)             H(n1, n2)
!>  H (n2, n1)             H(n2, n2)
!>~~~
!> which yields sin(2 atan ( H(n1,n2) / 2(H(n1,n1)-H(n2,n2))))
!> for the off symmetry elements of the rotation matrix.
!***********************************************************************

      SUBROUTINE ROT1(NBANDS,CHF,BSTEP,CUNI,CEIDB,CTMP)

__GPU USE moffload_struct_def

      USE prec

      IMPLICIT NONE
      INTEGER NBANDS
      GDEF CHF(NBANDS,NBANDS)
      GDEF CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &     CEIDB(NBANDS,NBANDS)
      REAL(q) BSTEP
   ! local
      INTEGER N1,N2
      REAL(q) HFEIG(NBANDS),THETA
      COMPLEX(q) CH

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CHF,CEIDB) PRIVATE(CH,THETA) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:CHF,CEIDB) PRIVATE(CH,THETA) IF(OFFLOAD_ON)
      DO N1=1,NBANDS
         DO N2=1,NBANDS
            CH= CHF(N1,N2)
            THETA=.5_q*ATAN(2*ABS(CH)*(-BSTEP))
! store in the  transposed of CEIDB because MATMUL uses the transposed ...
            IF (ABS(CH)/=0) THEN
               CEIDB(N2,N1)=SIN(THETA)*CH/MAX(ABS(CH),1E-8_q)
            ELSE
               CEIDB(N2,N1)=0
            ENDIF
         ENDDO
      ENDDO

GPACC !$ACC PARALLEL LOOP PRESENT(CEIDB) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CEIDB) IF(OFFLOAD_ON)
      DO N1=1,NBANDS
         CEIDB(N1,N1)=1
      ENDDO

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE ROT2   ****************************
!
!> This subroutine calculates the rotation matrix
!>~~~
!>    -CHF d          - E d  +                   +
!>   e            = U e      U    with CHF = U E U
!>~~~
!> For one k-point, the Loewdin-perturbation CHF matrix must be supplied.
!> The resulting matrix is unitary
!> i.e. U is the matrix which diagonolizes CHF   and E is the diagonal
!>      matrix containing the eigenvalues
!> up to first order the result from ROT1 and ROT2 are identical.
!
!***********************************************************************

      SUBROUTINE ROT2(COMM, NBANDS,CHF,BSTEP,CUNI,CEIDB,CTMP)

__GPU USE moffload_struct_def

      USE mpimy
      USE prec
      USE dfast
      USE scala
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE(communic) :: COMM
      GDEF CHF(NBANDS,NBANDS)
      GDEF CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &     CEIDB(NBANDS,NBANDS)

      DIMENSION HFEIG(NBANDS),W(3*NBANDS)

      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE

!      CALL SETUP_PARALLEL_GEMM(COMM, NBANDS, PGEMM_HANDLE)
!=======================================================================
! Diagononalize the resulting Matrix (search direction)
!=======================================================================
      CALL WGCOPY(NBANDS*NBANDS,CHF(1,1),1,CUNI(1,1),1)


#if defined(MPI)
      ! scaLAPACK if available in parallel version
      ! unfortunately the matrix seems to be unsuitable for scaLAPACK
      ! since it can be very close to 0
      ! so avoid calling scaLAPACK routines
      IF (__IF_OFFLOAD_OFF__( LscaLAPACK .AND. .FALSE. )) THEN
         CALL PDSYEVX_ZHEEVX(COMM, CUNI, HFEIG, NBANDS , NBANDS)
         CALLMPI( M_sum_g(COMM,CUNI, NBANDS*NBANDS))
      ELSE
#endif
#ifdef  gammareal
         CALL WDSYEV &
              ('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
              IFAIL)
#else
#ifdef CRAYHIP
! crayhip zheev needs RWORK on the device
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W) IF(OFFLOAD_ON)
#endif // CRAYHIP
         CALL WZHEEV &
              ('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
              W,  IFAIL)
#ifdef CRAYHIP
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W) IF(OFFLOAD_ON)
#endif // CRAYHIP
#endif
         IF (IFAIL/=0) THEN
            CALL vtutor%error("ERROR ROTDIA: Call to routine ZHEEV failed! Error code was " // str(IFAIL))
         ENDIF
#if defined(MPI)
      ENDIF
#endif
!=======================================================================
! set up the unitary transformation CEIDB for wavefunctions
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CUNI,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:CUNI,CTMP) IF(OFFLOAD_ON)
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CUNI(N1,N2)*EXP(-BSTEP*HFEIG(N2))
      ENDDO
      ENDDO

! currently there is an issue since PARALLEL_GGEMM does not work
! if second matrix is transposed
!      CALL PARALLEL_GGEMM(PGEMM_HANDLE, 'N',trans, NBANDS, NBANDS, NBANDS, one, CTMP, &
!     &             NBANDS, CUNI, NBANDS, zero, CEIDB, NBANDS)
      CALL WGGEMM('N',trans, NBANDS, NBANDS, NBANDS, one, CTMP, &
     &             NBANDS, CUNI, NBANDS, zero, CEIDB, NBANDS)
!=======================================================================
! rotate the wavefunction (with CEIDB): warning about MATMUL! We must
! take the transposed of CEIDB because MATMUL uses the transposed ...
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,CEIDB) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:CTMP,CEIDB) IF(OFFLOAD_ON)
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CEIDB(N2,N1)
      ENDDO
      ENDDO

      CALL WGCOPY(NBANDS*NBANDS,CTMP(1,1),1,CEIDB(1,1),1)

!      CALL RELEASE_PARALLEL_GEMM(PGEMM_HANDLE)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE ETAINI    *************************
!
!> This subroutine calculates the preconditioned step for the pseudo
!> Hamiltonian (delta eta) from CHAM
!> see FBN: Freysoldt, Boeck, Neugebauer, Phys. Rev. B 79, 241103 (2009).
!>
!> @note that VASP stores search direction respectively the negative
!> gradient
!>
!>   W%CELTOT        diagonal components of <phi_n | H | phi_n >
!>                   most likely identical to diagonal of CHAM \n
!>   W_F%CELTOT      diagonal components of current pseudo Hamiltonian
!>                   (epsilon in FBN) \n
!>   CHAM on entry:  Hamiltonian matrix <phi_n | H | phi_m > \n
!>        on exit:   search direction for pseudo Hamiltonian eta \n
!>   CHF  on entry:  previous search direction for pseudo Hamiltonian eta \n
!>        on exit:   previous search direction for pseudo Hamiltonian eta
!>
!>   CHAMP on entry: preconditioned Hamilton matrix (optional) \n
!>   ORT2            orthogonality of current gradient to previous rotation
!>                   matrix \n
!>   DE2             expected energy change along present rotation matrix \n
!>   WEIMIN          threshhold for empty bands \n
!>   CEIDB           work array
!
!***********************************************************************

      SUBROUTINE ETAINI(WDES, W, W_F, CHAM, CHF, CEIDB, KPOINTS, EFERMI, WEIMIN, &
           KAPPA, ORT2, DE2, ORTCEL, DECEL, LCHCON, CHAMP)

__GPU USE moffload_struct_def

      USE prec
      USE wave
      USE constant
      USE poscar
      USE pseudo
      USE mgrid
      USE lattice
      USE nonl_high
      USE augfast
      USE mkpoints
      USE fermi_energy, ONLY : DELSTP
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W, W_F

      GDEF CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      GDEF CHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      GDEF CEIDB(WDES%NB_TOT,WDES%NB_TOT)
      GDEF, OPTIONAL :: CHAMP(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      TYPE (kpoints_struct) KPOINTS
      REAL(q)        :: EFERMI   ! chemical potential, Fermi energy
      REAL(q)        :: KAPPA    ! scaling factor
      REAL(q) WEIMIN
      REAL(q) ORT2, DE2
      REAL(q) ORTCEL, DECEL

      LOGICAL        :: LCHCON   ! charge updated
! local
      INTEGER ISP, NK, N, N1, N2, NC, MYKSTART, MYKINCR
      REAL(q) WEIGHT, DIFFER, DIFCEL, THETA
      REAL(q), PARAMETER ::  DIFMAX=1E-8   ! threshhold for degeneragy
      GDEF CROT, CGRAD
      REAL(q) :: EH, WSUM, DFUN, SFUN, X, EFERMI_SHIFT, WMIN, DER
      REAL(q) :: ORT2_in,DE2_in,DECEL_in,ORTCEL_in
GPACC !$ACC ROUTINE(DELSTP) SEQ

#ifdef MPI
      MYKSTART=WDES%COMM_KINTER%NODE_ME
      MYKINCR =WDES%COMM_KINTER%NCPU
#else
      MYKSTART=1
      MYKINCR =1
#endif
!=======================================================================
! "gradient" for diagonal components of the pseudo Hamiltonian
! store result in the diagonal of CHAM
! these are identical to the gradient for the eigenvalues in
! the old algorithm
!=======================================================================
      EFERMI_SHIFT=0
      WSUM =0

!      depsilon_Fermi =
!     [\sum_n df_n/d(epsilon_n-epsilon_Fermi) x d epsilon_n ]/ [\sum_n df_n/d(epsilon_n-epsilon_Fermi)]
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG REDUCTION(+:WSUM,EFERMI_SHIFT) PRESENT(WDES,W,W_F) PRIVATE(WEIGHT) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) REDUCTION(+:WSUM,EFERMI_SHIFT) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F) PRIVATE(WEIGHT) IF(OFFLOAD_ON)
      DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
            WEIGHT =WDES%WTKPT(NK)*WDES%RSPIN
GPACC !$ACC LOOP VECTOR REDUCTION(+:WSUM,EFERMI_SHIFT) PRIVATE(X,DFUN,SFUN)
GPOMP !$OMP PARALLEL DO REDUCTION(+:WSUM,EFERMI_SHIFT) PRIVATE(X,DFUN,SFUN)
            DO N =1,W%WDES%NB_TOTK(NK,ISP)
               X=  (EFERMI- REAL( W_F%CELTOT(N,NK,ISP) ,KIND=q) )/KPOINTS%SIGMA
               CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)
               WSUM =WSUM +DFUN/KPOINTS%SIGMA* WEIGHT
               EFERMI_SHIFT=EFERMI_SHIFT+DFUN*(W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))/KPOINTS%SIGMA*WEIGHT
            ENDDO
         ENDDO
      ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for WSUM and EFERMI_SHIFT
      CALLMPI(M_sum_d(WDES%COMM_KINTER,WSUM,1))
      CALLMPI(M_sum_d(WDES%COMM_KINTER,EFERMI_SHIFT,1))

      ! shift is the first order change in the Fermi-energy
      IF (ABS(WSUM)>1E-10) THEN
         EFERMI_SHIFT=EFERMI_SHIFT/WSUM
      ELSE
         EFERMI_SHIFT=0
      ENDIF
      WMIN=WSUM/WDES%NB_TOT*WEIMIN

      DECEL_in=DECEL
      ORTCEL_in=ORTCEL
      DECEL=0
      ORTCEL=0

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG REDUCTION(+:DECEL,ORTCEL) PRESENT(WDES,W,W_F,CHAM,CHF) PRIVATE(WEIGHT)__IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) REDUCTION(+:DECEL,ORTCEL) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,W,W_F,CHAM,CHF) PRIVATE(WEIGHT) IF(OFFLOAD_ON)
      DO ISP=1,WDES%ISPIN
         DO NK=MYKSTART,WDES%NKPTS,MYKINCR
            WEIGHT =WDES%WTKPT(NK)*WDES%RSPIN
GPACC !$ACC LOOP VECTOR REDUCTION(+:DECEL,ORTCEL) PRIVATE(X,DFUN,SFUN)
GPOMP !$OMP PARALLEL DO REDUCTION(+:DECEL,ORTCEL) PRIVATE(X,DFUN,SFUN)
            DO N =1,W%WDES%NB_TOTK(NK,ISP)
               X=  (EFERMI- REAL( W_F%CELTOT(N,NK,ISP) ,KIND=q) )/KPOINTS%SIGMA
               CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)

               ! avoid conjugation of those elements which do not contribute to energy
               IF (ABS(DFUN)<WMIN) THEN
                   CHF(N,N,NK,ISP)=0
               ENDIF

               DFUN= DFUN*((W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP)-EFERMI_SHIFT)/KPOINTS%SIGMA)

               IF (PRESENT(CHAMP)) THEN
                  SFUN= (W_F%CELTOT(N,NK,ISP)-CHAMP(N,N,NK,ISP))*KAPPA
               ELSE
                  ! preconditioned gradient equal gradient
                  SFUN=(W_F%CELTOT(N,NK,ISP)-W%CELTOT(N,NK,ISP))*KAPPA
               ENDIF

               ! first order change for move along diagonal components
               DECEL  =DECEL  +DFUN* SFUN*WEIGHT
               ORTCEL =ORTCEL +DFUN* CHF(N,N,NK,ISP)*WEIGHT
               CHAM(N,N,NK,ISP)=SFUN
            ENDDO
         ENDDO
      ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for DECEL and ORTCEL
      CALLMPI(M_sum_d(WDES%COMM_KINTER,DECEL,1))
      CALLMPI(M_sum_d(WDES%COMM_KINTER,ORTCEL,1))

      DECEL=DECEL+DECEL_in
      ORTCEL=ORTCEL+ORTCEL_in

!=======================================================================
! calculate off diagonal elements of rotation matrix
!=======================================================================
      ORT2_in=ORT2
      DE2_in=DE2
      ORT2=0
      DE2=0

      EH=0
!     IF (.NOT. LCHCON) THEN
!        EH=0.35
!     ELSE
!        EH=0
!     ENDIF

      DO ISP=1,WDES%ISPIN
      DO NK=MYKSTART,WDES%NKPTS,MYKINCR
         WEIGHT =WDES%WTKPT(NK)
GPACC !$ACC PARALLEL LOOP GANG COLLAPSE(2) REDUCTION(-:DE2,ORT2) PRESENT(CEIDB,WDES,W,CHAM,CHAMP,CHF) PRIVATE(DIFFER,CGRAD,DIFCEL,CROT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) REDUCTION(+:DE2,ORT2) &
GPOMP !$OMP MAP(PRESENT,TO:CEIDB,WDES,W,CHAM,CHF) MAP(ALLOC:CHAMP) PRIVATE(DIFFER,CGRAD,DIFCEL,CROT) IF(OFFLOAD_ON)
         DO N2=1,WDES%NB_TOT
            DO N1=1,WDES%NB_TOT
              IF (N1/=N2) THEN
                 DIFFER=     (W%FERTOT(N1,NK,ISP)-W%FERTOT(N2,NK,ISP))
                 ! gradient
                 CGRAD=WDES%RSPIN*CHAM(N2,N1,NK,ISP)*DIFFER

                 DIFCEL= REAL( W_F%CELTOT(N2,NK,ISP)-W_F%CELTOT(N1,NK,ISP) ,KIND=q)
                 IF (ABS(DIFCEL)>DIFMAX) THEN
                    CGRAD= CGRAD/DIFCEL
                 ELSE
                    CGRAD=0
                 ENDIF

                 ! rotation matrix supplied by calling routine in CHAMP
                 IF (PRESENT(CHAMP)) THEN
                    CROT  = CHAMP(N1,N2,NK,ISP)
                 ELSE
                    CROT  = CHAM(N1,N2,NK,ISP)
                 ENDIF

                 ! slow down rotation if occupancies are widely different and energies close
                 ! (not applied right now since EH=0)
                 IF (ABS(DIFCEL)>DIFMAX) THEN
                    IF (DIFCEL>0) THEN
                       CROT=CROT*(DIFCEL/(DIFCEL+2*EH*ABS(2*W%FERTOT(N1,NK,ISP)-2*W%FERTOT(N2,NK,ISP))))
                    ELSE
                       CROT=CROT*(DIFCEL/(DIFCEL-2*EH*ABS(2*W%FERTOT(N1,NK,ISP)-2*W%FERTOT(N2,NK,ISP))))
                    ENDIF
                 ENDIF

                 CROT=CROT*KAPPA

                 ! avoid conjugation of those elements which do not contribute to energy
                 ! these are all elements that have similar fermi-weights
                 ! since rotation in that sub-space leaves the energy invariant
                 IF (ABS(DIFFER)<WEIMIN) THEN
                    CHF(N1,N2,NK,ISP)=0
                 ENDIF
                 CEIDB(N1,N2)=-CROT
                 ! gradient times previous search direction
                 ORT2  =ORT2  -CGRAD*CHF(N1,N2,NK,ISP)*WEIGHT
                 ! gradient times present precond. gradient
                 DE2   =DE2   -CGRAD*CEIDB(N1,N2)*WEIGHT
              ELSE
                 CEIDB(N1,N2)=0
              ENDIF
            ENDDO
         ENDDO

         ! store negative search direction
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES,CEIDB,CHAM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:WDES,CEIDB,CHAM) IF(OFFLOAD_ON)
         DO N2=1,WDES%NB_TOT
            DO N1=1,WDES%NB_TOT
               IF (N1/=N2) THEN
                  CHAM(N1,N2,NK,ISP)=CEIDB(N1,N2)
               ENDIF
            ENDDO
         ENDDO
!         CALL DUMP_HAM( "rotation matrix",WDES, CHAM(1,1,NK,ISP))

      ENDDO
      ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON) ! wait for ORT2 and DE2
      CALLMPI (M_sum_d(WDES%COMM_KINTER,ORT2,1) )
      CALLMPI (M_sum_d(WDES%COMM_KINTER,DE2,1) )

      ORT2=ORT2+ORT2_in
      DE2=DE2+DE2_in

    END SUBROUTINE ETAINI


!************************ SUBROUTINE ROTETA ****************************
!
!> This subroutine diagonalizes the rotation matrix eta [Equ. (32) in
!> see FBN: Freysoldt, Boeck, Neugebauer, Phys. Rev. B 79, 241103 (2009)].
!>
!>   epsilon of pseudo Hamiltonian + search direction * stepsize
!>
!>   CELTOT + CHF * BSTEP  = U E U+
!>
!> see FBN: Freysoldt, Boeck, Neugebauer, Phys. Rev. B 79, 241103 (2009).
!>
!> CHF    on entry: search direction for the pseudo Hamiltonian \n
!>        on exit:  unitarily transformed \n
!>                  search direction for the pseudo Hamiltonian U+ CHF U \n
!> BSTEP  step size \n
!> CEIDB  on exit: unitary rotation matrix U \n
!>        this rotation must be applied to the orbitals \n
!> CELTOT on entry: previous pseudo Hamiltonian (diagonal eigenvalues) \n
!>        on exit:  new diagonal components of pseudo Hamiltonian \n
!>                  after unitary rotation
!>
!> CUNI and CTMP are auxilary matrices.
!
!***********************************************************************

      SUBROUTINE ROTETA(COMM,NBANDS,CHF,CELTOT,BSTEP,CUNI,CEIDB,CTMP,restore_W,AUXTOT)

__GPU USE moffload

      USE mpimy
      USE prec
      USE scala
      USE dfast
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE(communic) :: COMM
      GDEF CHF(NBANDS,NBANDS)
      GDEF CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &     CEIDB(NBANDS,NBANDS)
      COMPLEX(q) :: CELTOT(NBANDS)
      REAL(q) :: HFEIG(NBANDS)
      REAL(q),OPTIONAL :: AUXTOT(NBANDS)  ! auxiliary array that needs to be rotated
      LOGICAL :: restore_W                ! if false rotate CHF to new basis
    ! local
      INTEGER N1, N2, IFAIL
      REAL(q) :: W(3*NBANDS)
      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE

      CALL SETUP_PARALLEL_GEMM(COMM, NBANDS, PGEMM_HANDLE)
!=======================================================================
! diagonalize the resulting matrix (search direction)
!
! first conjugate since all Hamilton related quantities are complex
! conjugated w.r.t what VASP uses in rot.F
!  CHF(n2,n1) proto < phi_n2 | H | phi_n1 >* = < phi_n1 | H | phi_n2 >
!=======================================================================
GPACC !$ACC KERNELS PRESENT(CHF) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHF) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CHF=GCONJG(CHF)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP GANG PRESENT(CUNI,CHF,CELTOT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:CUNI,CHF,CELTOT) IF(OFFLOAD_ON)
      DO N1=1,NBANDS
GPACC !$ACC LOOP VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
         DO N2=1,NBANDS
            CUNI(N1,N2)=-CHF(N1,N2)*BSTEP
         ENDDO
         CUNI(N1,N1)=CELTOT(N1)-CHF(N1,N1)*BSTEP
      ENDDO

GPACC !$ACC ENTER DATA CREATE(HFEIG) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:HFEIG) IF(OFFLOAD_ON)

#if defined(MPI)
      ! scaLAPACK if available in parallel version
      IF (__IF_OFFLOAD_OFF__(LscaLAPACK)) THEN
!         CALLMPI( M_bcast_z(COMM, CUNI, NBANDS*NBANDS))  ! ZHEEVD seems to require a bcast here
!         CALL pDSYEV_ZHEEV(COMM, CUNI, HFEIG, NBANDS , NBANDS)
         CALL PDSYEVX_ZHEEVX(COMM, CUNI, HFEIG, NBANDS , NBANDS)
!         CALL pDSYEV_ZHEEVD(COMM, CUNI, HFEIG, NBANDS , NBANDS)
         CALLMPI( M_sum_g(COMM,CUNI, NBANDS*NBANDS))
      ELSE
#endif
#ifdef  gammareal
         CALL WDSYEV &
              ('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
              IFAIL)
         CALLMPI( M_bcast_d(COMM, CUNI, NBANDS*NBANDS))
         CALLMPI( M_bcast_d(COMM, HFEIG, NBANDS))
#else
#ifdef CRAYHIP
! crayhip zheev needs RWORK on the device
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W) IF(OFFLOAD_ON)
#endif // CRAYHIP
         CALL WZHEEV &
              ('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
              W,  IFAIL)
#ifdef CRAYHIP
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W) IF(OFFLOAD_ON)
#endif // CRAYHIP
         CALLMPI( M_bcast_z(COMM, CUNI, NBANDS*NBANDS))
         CALLMPI( M_bcast_d(COMM, HFEIG, NBANDS))
#endif
         IF (IFAIL/=0) THEN
            CALL vtutor%error("ERROR ROTDIA: Call to routine ZHEEV failed! Error code was " // str(IFAIL))
         ENDIF
#if defined(MPI)
      ENDIF
#endif
      ! store current (updated) eigenvalues back in CELTOT
GPACC !$ACC KERNELS PRESENT(CELTOT,HFEIG) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CELTOT,HFEIG) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CELTOT=HFEIG
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC EXIT DATA DELETE(HFEIG) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:HFEIG) IF(OFFLOAD_ON)

!=======================================================================
! transform AUXTOT (conditioning step with for each orbital)
!=======================================================================
      IF (PRESENT(AUXTOT) .AND. .NOT. restore_W) THEN
      ! CTMP = AUXTOT * U
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,AUXTOT,CUNI) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:CTMP,AUXTOT,CUNI) IF(OFFLOAD_ON)
         DO N1=1,NBANDS
            DO N2=1,NBANDS
               CTMP(N1,N2)=AUXTOT(N1)*CUNI(N1,N2)
            ENDDO
         ENDDO
      ! CEIDB = U+ CTMP
         CALL PARALLEL_GGEMM( PGEMM_HANDLE, trans, 'N', NBANDS, NBANDS, NBANDS, one, CUNI, &
         &             NBANDS, CTMP, NBANDS, zero, CEIDB, NBANDS)
      ! store diagonals back to AUXTOT
GPACC !$ACC PARALLEL LOOP PRESENT(AUXTOT,CEIDB) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:AUXTOT,CEIDB) IF(OFFLOAD_ON)
         DO N1=1,NBANDS
            AUXTOT(N1)=CEIDB(N1,N1)
         ENDDO
      ENDIF
!=======================================================================
! transform CHF (equ. (21) of FBN)
!=======================================================================
      IF (.NOT. restore_W) THEN
      ! CEIDB = U+ CHF
        CALL PARALLEL_GGEMM(PGEMM_HANDLE, trans, 'N', NBANDS, NBANDS, NBANDS, one, CUNI, &
           &             NBANDS, CHF, NBANDS, zero, CEIDB, NBANDS)

     ! CHF = (U+ CHF ) U
        CALL PARALLEL_GGEMM(PGEMM_HANDLE, 'N', 'N', NBANDS, NBANDS, NBANDS, one, CEIDB, &
           &             NBANDS, CUNI, NBANDS, zero, CHF, NBANDS)
      ENDIF
!=======================================================================
! rotate the wavefunction (with CEIDB)
!=======================================================================
GPACC !$ACC KERNELS PRESENT(CEIDB,CUNI,CHF) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CEIDB,CUNI,CHF) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CEIDB=CUNI
      CHF=GCONJG(CHF)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      CALL RELEASE_PARALLEL_GEMM(PGEMM_HANDLE)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE polynomial_fit ********************
!
! > perform a polynomial fit through n_data points
! > supplied in the arrays x and y
! > return the coefficients of the polynomial fit of order in coeffs
! > search for the minimum using quasi-Newton close to supplies y_min
! > calling routine must supply x_min and y_min
! > if fit fails no update of x_min and y_min is performed
!
!***********************************************************************

!=======================================================================
! find a suitable slot in logarithmic grid
! if x in between two already considered slots no updated of x
! otherwise set x to the closest slot
! however do not allow to move too far from already used slots
!=======================================================================

  subroutine slot_in(x, x_centre, x_scale, max_slot_move, lslot , l_is_updated , l_sloppy, iu0)
    implicit none

    real(q) :: x        ! in: x out: energy of next empty slot
    real(q) :: x_centre ! central slot
    real(q) :: x_scale  ! logarithmic scaline
    logical :: lslot(:) ! slots
    integer :: max_slot_move ! maximum move away from already occupied slots
    logical :: l_is_updated ! return true if a slot has been found
    logical :: l_sloppy ! requires only the one (not two closest slots) are occupied
    integer :: iu0
   ! local
    integer :: nmiddle  ! central slot
    real(q) :: aslot
    integer :: nslot   ! current slot
    integer :: i

    nmiddle=size(lslot)/2  ! central slot
    aslot=LOG(ABS(max(x,1E-10_q)/x_centre))/LOG(x_scale)
    nslot=NINT(aslot)
    l_is_updated=.false.
    if (nslot+1+nmiddle > size(lslot)) nslot=size(lslot)-nmiddle-1
    if (nslot-1+nmiddle < 1) nslot=-nmiddle+2

  ! slot already occupied
    if (lslot(nslot+nmiddle)) then
      if (l_sloppy) return
      if (aslot-nslot > 0 ) then ! rounding was to lower slot
         if (lslot(nslot+1+nmiddle)) return ! both neighboring slots are set
         nslot=nslot+1
      else ! rounding was to upper slot
         if (lslot(nslot-1+nmiddle)) return ! both neighboring slots are set
         nslot=nslot-1
      endif
    endif
   ! now start from outermost slot and work to central slot
   ! determine outermost occupied slot
  ! IF (iu0>=0) WRITE(iu0,'(" start slot  ",I4,F8.3,14L2," X",14L2)') &
  !        nslot,aslot,lslot(nmiddle-14:nmiddle-1),lslot(nmiddle+1:nmiddle+14)
   if (nslot>0) THEN
     do i=size(lslot)-nmiddle,0,-1
       if (lslot(i+nmiddle)) exit
     enddo
     if (nslot>i+max_slot_move) nslot=i+max_slot_move
   else if (nslot<0) THEN
     do i=-nmiddle+1,0,1
       if (lslot(i+nmiddle)) exit
     enddo
     if (nslot<i-max_slot_move) nslot=i-max_slot_move
   endif
   ! this slot should be empty for sure
   if (lslot(nslot+nmiddle)) then
     IF (iu0>=0) THEN
        WRITE(iu0,*) 'internal error in rot.F: slot_in the slot is occupied'
        WRITE(iu0,'(" start slot  ",I4,F8.3,30L2)') nslot,aslot,lslot(nmiddle-14:nmiddle+15)
     endif
     STOP
   endif
   lslot(nslot+nmiddle)=.TRUE.
  ! IF (iu0>=0) WRITE(iu0,'(" end slot  ",I4,2F14.7)') nslot, x, exp(log(x_scale)*(nslot))*x_centre
   x=exp(log(x_scale)*nslot)*x_centre
   l_is_updated=.true.

  end subroutine slot_in


!=======================================================================
! perform a polynomial fit and search for minimum in the fitted function
!=======================================================================

  subroutine search_minimum_polynomial_fit(x, y, n_data, x_min, y_min, f_double_prime, info, iu0)
    implicit none

    integer :: n_data
    real(q) :: x(:), y(:)
    real(q) :: coeffs(n_data )
    real(q) :: x_min, y_min, f_double_prime
    integer :: iu0
  ! local
    integer :: order, info
    real(q) :: a, b, tolerance
    integer :: ind(1)
    order=2
    if (n_data>=4) order=3
    if (n_data>=6) order=5
    if (n_data>=8) order=6
    if (iu0>=0) write(iu0,'(2I3)', ADVANCE="no") order, n_data
  ! perform fit
    call polynomial_fit(x(1:n_data), y(1:n_data), n_data, order, coeffs, info)
  ! it is assumed x_min  and y_min have estimates, fall back without mods
    if (info >0 ) return
    a=x_min
    tolerance=(maxval(x(1:n_data))-minval(x(1:n_data)))/1e8
    call find_minimum(coeffs, order, a, tolerance, x_min, info, iu0)
  !  do i = 1, n_data
  !     if (iu0>=0) write(iu0,'(3E14.5)') x(i), y(i),evaluate_polynomial(coeffs, order, x(i))
  !  enddo
    y_min=evaluate_polynomial(coeffs, order, x_min)
    f_double_prime = evaluate_second_derivative(coeffs, order, x_min)

   ! backup if y_min is larger then smallest supplied value
    ind=MINLOC(y(1:n_data))
    if (y_min > y(ind(1))) then
       info=100
       y_min=y(ind(1))
       x_min=x(ind(1))
    endif

  end subroutine

!=======================================================================
! perform a spline fit and search for minimum in the fitted function
! assumes that two close by points are supplied at the start
! yielding the first derivative of the function at 0
!=======================================================================

  subroutine search_minimum_spline(x, y, n_data, x_min, y_min, f_double_prime, info, iu0)
    implicit none

    integer :: n_data        ! number of data points
    real(q) :: x(:), y(:)    ! data points
    real(q) :: x_min, y_min, f_double_prime ! return x,y and second derivative
    integer :: info
    integer :: iu0           ! unit to perform IO
  ! local
    real(q) :: p(n_data-1, 5)  ! coefficients for spline
    real(q) :: a, f, f_prime, tolerance, x_max
    integer :: n, i, i_found
    integer :: index(n_data), ind(1)
    integer :: i_max=128

    info=0
  ! first point is special (see main code)
    p(1,1)=(x(1)+x(2))/2
    p(1,2)=(y(1)+y(2))/2
    f_prime=(y(2)-y(1))/(x(2)-x(1))

  ! transfer other data sorted (except for first slot)
    n=n_data-1
    call sort_with_index(x(3:), index)
    do i=1,n-1
      p(i+1,1)=x(index(i)+2)
      p(i+1,2)=y(index(i)+2)
    enddo
!   test spline fit (data and derivate)
!    if (iu0>=0) write(iu0,'(4E14.6)') p(1,1:2),f_prime
!    do i=2,n
!      if (iu0>=0) write(iu0,'(4E14.6)') p(i,1:2),0.0_q
!    enddo
!    if (iu0>=0) write(iu0,*)

    call splcof(p, n, n, f_prime)

  ! search for approximate minimum start with 2nd point
    a=x(2)
    y_min=y(2)
    x_max=maxval(x(2:n_data))*1.5 ! avoid endppoint
    do i=1,i_max
      call splval_2nd( (x_max/i_max)*i, f, f_prime, f_double_prime, p, n, n)
!      if (iu0>=0) write(iu0,'(4E14.6)') (x_max/i_max)*i, f, f_prime, f_double_prime
      if (f < y_min) then
        y_min=f
        a=(x_max/i_max)*i
      endif
    enddo
!    ind=MINLOC(y(1:n_data))
!    a=x(ind(1))
    tolerance=(maxval(x(1:n_data))-minval(x(1:n_data)))/1e8
    call find_minimum_spline(p, n, a, tolerance, x_min, info, iu0 )
    call splval_2nd( x_min, y_min, f_prime, f_double_prime, p, n, n)

   ! backup if y_min is larger then smallest supplied value
    ind=MINLOC(y(1:n_data))
    if (y_min > y(ind(1))) then
       info=100
       y_min=y(ind(1))
       x_min=x(ind(1))
    endif

  end subroutine

!=======================================================================
! perform a polynomial fit and search for minimum in the fitted function
! here select at most order_in+1 data points close to minimum
!=======================================================================

  subroutine search_minimum_polynomial_fit_order(x, y, n_data, order_in, x_min, y_min, f_double_prime, info, iu0)
    implicit none

    real(q) :: x(:), y(:)
    integer :: n_data
    integer :: order_in
    real(q) :: coeffs(n_data )
    real(q) :: distance(n_data)
    real(q) :: x_min, y_min, f_double_prime
    real(q) :: x_selected(order_in+1), y_selected(order_in+1)
    integer :: info
    integer :: iu0
  ! local
    integer :: ind(1)
    integer :: index(n_data)
    integer :: order, i
    real(q) :: a, tolerance, x_max, f
    integer :: i_max=128

    order=min(order_in,n_data-1)
    ind=MINLOC(y(3:n_data))+2 ! exclude first two points (starting point)
    distance=abs(x(1:n_data)-x(ind(1)))
    call sort_with_index(distance(1:n_data),index)
  ! copy points closest to minimum
    do i=1,order+1
      x_selected(i)=x(index(i))
      y_selected(i)=y(index(i))
    enddo
!    if (iu0>=0) then
!       write(iu0,*)
!       do i = 1, order+1
!          write(iu0,'(3E14.5)') x_selected(i), y_selected(i)
!       enddo
!       write(iu0,*)
!    endif
  ! perform fit
    call polynomial_fit(x_selected(1:order+1), y_selected(1:order+1), order+1, order, coeffs, info)
  ! it is assumed x_min  and y_min have estimates, fall back without mods
    if (info >0 ) return

  ! search for approximate minimum start with 2nd point
    a=x(2)
    y_min=y(2)
    x_max=maxval(x_selected(1:order+1))*1.5
    do i=1,i_max
      f=evaluate_polynomial(coeffs, order, (x_max/i_max)*i)
!      if (iu0>=0) write(iu0,'(4E14.6)') (x_max/i_max)*i, f, evaluate_derivative(coeffs, order, (x_max/i_max)*i)
      if (f < y_min) then
        y_min=f
        a=(x_max/i_max)*i
      else
        exit
      endif
    enddo
    tolerance=(maxval(x(1:n_data))-minval(x(1:n_data)))/1e8
    call find_minimum(coeffs, order, a, tolerance, x_min, info, iu0)

    y_min=evaluate_polynomial(coeffs, order, x_min)
    f_double_prime = evaluate_second_derivative(coeffs, order, x_min)

  ! backup if y_min is larger then smallest supplied value
    ind=MINLOC(y(1:n_data))
    if (y_min > y(ind(1))) then
       info=100
       y_min=y(ind(1))
       x_min=x(ind(1))
    endif

  end subroutine

!=======================================================================
! perform a polynomial fit and search for minimum in the fitted function
! here select at most order_in+1+n_over data points close to minimum
! this improves the stability of the fit
!=======================================================================

  subroutine search_minimum_polynomial_fit_order2(x, y, n_data, order_in, n_over, x_min, y_min, f_double_prime, info, iu0)
    implicit none

    real(q) :: x(:), y(:)
    integer :: n_data
    integer :: n_over
    integer :: order_in
    real(q) :: coeffs(n_data )
    real(q) :: distance(n_data)
    real(q) :: x_min, y_min, f_double_prime
    integer :: info
    integer :: iu0
  ! local
    real(q) :: x_selected(n_data), y_selected(n_data)
    integer :: n_data_use
    integer :: ind(1)
    integer :: index(n_data)
    integer :: order, ndata_use, i
    real(q) :: a, tolerance, x_max, x_m, f
    integer :: i_max=128


    order=min(order_in,n_data-1-n_over)
    ind=MINLOC(y(3:n_data))+2 ! exclude first two points (derivative info)
!    if (iu0>=0) write(*,*) "n_data, order", n_data, order, order+1+n_over
    distance=abs(x(1:n_data)-x(ind(1)))
    call sort_with_index(distance(1:n_data),index)
  ! copy sorted points
    do i=1,n_data
      x_selected(i)=x(index(i))
      y_selected(i)=y(index(i))
    enddo
    n_data_use=order+1+n_over
  ! if the last sorted point corresponds to point 2 also include point 1 (derivative information)
    if (n_data_use<n_data .and. index(n_data_use)==2) then
      n_data_use=n_data_use+1
    endif
  ! perform fit
    call polynomial_fit(x_selected(1:n_data_use), y_selected(1:n_data_use), n_data_use, order, coeffs, info)
  ! it is assumed x_min  and y_min have estimates, fall back without mods
    if (info >0 ) return

!    if (iu0>=0) then
!       write(iu0,*)
!       do i = 1, n_data_use
!          write(iu0,'(3E24.12)') x_selected(i), y_selected(i), evaluate_polynomial(coeffs, order, x_selected(i))
!       enddo
!       write(iu0,*)
!    endif

  ! search for approximate minimum start at smallest selected x
    a=x(1)
    y_min=y(1)
    x_m=abs(minval(x_selected(1:n_data_use)))
    x_max=maxval(x_selected(1:n_data_use))*1.1
    do i=1,i_max
      f=evaluate_polynomial(coeffs, order, ((x_max-x_m)/i_max)*i+x_m)
!      if (iu0>=0) write(iu0,'(4E14.6)') ((x_max-x_m)/i_max)*i+x_m, f, evaluate_derivative(coeffs, order, (x_max/i_max)*i)
      if (f < y_min) then
        y_min=f
        a=((x_max-x_m)/i_max)*i+x_m
      else
        exit
      endif
    enddo
    tolerance=(maxval(x(1:n_data))-minval(x(1:n_data)))/1e8
    call find_minimum(coeffs, order, a, tolerance, x_min, info, iu0)

    y_min=evaluate_polynomial(coeffs, order, x_min)
    f_double_prime = evaluate_second_derivative(coeffs, order, x_min)
!    if (info>0) then
!       if (iu0>=0) write(iu0,'("  minimum not found start, found",2E14.6)') a, x_min
!    else
!       if (iu0>=0) write(iu0,'("  minimum start, found",2E14.6)') a, x_min
!    endif

  ! backup if y_min is larger then smallest supplied value
!    ind=MINLOC(y(1:n_data))
!    if (y_min > y(ind(1))) then
!       info=100
!       y_min=y(ind(1))
!       x_min=x(ind(1))
!    endif

  end subroutine

!=======================================================================
! quick simple bubble sort returing an index array
!=======================================================================

  subroutine sort_with_index(array, index)
    implicit none

    real(q), intent(in) :: array(:)
    integer, intent(out) :: index(:)
    integer :: i, j, temp_index
    real :: temp_value

    ! Initialize the index array (initially, indices are 1, 2, 3...)
    index = [(i, i=1, size(array))]

    ! Perform a bubble sort
    do i = 1, size(array) - 1
      do j = 1, size(array) - i
        if (array(index((j))) > array(index(j + 1))) then
          ! Swap indices
          temp_index = index(j)
          index(j) = index(j + 1)
          index(j + 1) = temp_index
        end if
      end do
    end do

  end subroutine sort_with_index

!=======================================================================
! perform a polynomial fit of order through n_data data points
! uses simple normal equation
!=======================================================================

  subroutine polynomial_fit(x, y, n_data, order, coeffs, info)
    implicit none

    integer :: n_data, order
    real(q) :: x(:), y(:)
    real(q) :: coeffs(order + 1)
    integer :: info
  ! local
    integer :: M               ! Order + 1
    real(q) :: A(n_data, order + 1)  ! Vandermonde matrix
    real(q) :: A_tA(order + 1, order + 1) ! A^T * A
    real(q) :: A_ty(order + 1)      ! A^T * y
    integer :: i, j
    integer :: ipiv(order + 1)    ! Pivot array for DSYTRF/DSYTRS
    real(q) :: work((order+1)*(order+1))  ! Workspace for DSYTRF (adjust if needed)
    real(q) :: piv(n_data)


    M = order + 1
    info=0
    ! Weight linear equations close to minimum in y more strongly (factor 10)
    ! this has no effect for fully determined equations
    piv(:)=(maxval(y)-minval(y))/(y(:)-minval(y)+(maxval(y)-minval(y))*0.1)

    ! Construct the Vandermonde matrix A
    do i = 1, n_data
      do j = 1, M
        A(i, j) = piv(i)*x(i)**(j - 1)
      end do
    end do

    ! Calculate A^T * A and A^T * y
    A_tA = matmul(transpose(A), A)
    A_ty = matmul(transpose(A), piv*y)

    ! Solve the normal equations using LAPACK (DSYTRF and DSYTRS)
    call dsytrf('U', M, A_tA, M, ipiv, work, (order+1)*(order+1), info)
    if (info /= 0) then
      return
    end if

    call dsytrs('U', M, 1, A_tA, M, ipiv, A_ty, M, info)
    if (info /= 0) then
      return
    end if

    coeffs = A_ty  ! Solution is in A_ty

  end subroutine polynomial_fit

!=======================================================================
! find minimum of a polynomial function using Newton's algorithm
!=======================================================================

  subroutine find_minimum(coeffs, n, a, tolerance, x_min, info, iu0)
    implicit none

    integer:: n
    integer:: info, iu0
    real(q):: coeffs(n + 1)
    real(q):: a, tolerance, x_min, x_min_last
    real(q):: f, f_prime, f_double_prime
    integer:: iter
    integer:: max_iter=10

    x_min = a! Initial guess
    info = 0

    do iter = 1, max_iter
      f = evaluate_polynomial(coeffs, n, x_min)
      f_prime = evaluate_derivative(coeffs, n, x_min)
      f_double_prime = evaluate_second_derivative(coeffs, n, x_min)

      if (f_double_prime == 0.0_q) then
        exit  ! Ok lets hope that x_min is fine
      end if

      x_min_last = x_min
      x_min = x_min - f_prime / f_double_prime
!      if (iu0>=0) write(iu0,'("minp ",2I3, 4F14.8)') n, iter, x_min,x_min_last, f, f_prime

      if (abs(x_min-x_min_last) < tolerance) exit ! Check for convergence

    end do
    if (iter > max_iter) info=99

  end subroutine find_minimum

  subroutine find_minimum_spline(p, n, a, tolerance, x_min, info, iu0)
    implicit none

    integer:: n
    real(q):: p(n, 5)
    real(q):: a, tolerance, x_min, x_min_last
    integer:: info, iu0
    real(q):: f, f_prime, f_double_prime
    integer:: iter
    integer:: max_iter=10

    x_min = a! Initial guess
    info = 0

    do iter = 1, max_iter
      call splval_2nd( x_min, f, f_prime, f_double_prime, p, n, n)

      if (f_double_prime == 0.0_q) then
        exit  ! Ok lets hope that x_min is fine
      end if

      x_min_last = x_min
      x_min = x_min - f_prime / f_double_prime
!      if (iu0>=0) write(iu0,'("mins ",I3, 4F14.8)') n, x_min,x_min_last, f, f_prime

      if (abs(x_min-x_min_last) < tolerance) exit ! Check for convergence

    end do
    if (iter > max_iter) info=99

  end subroutine find_minimum_spline


  real(q) function evaluate_polynomial(coeffs, n, x)
    implicit none
    integer:: n
    real(q):: coeffs(n + 1), x
    integer:: i
    evaluate_polynomial = coeffs(1)! Initialize with the constant term

    do i = 2, n + 1
      evaluate_polynomial = evaluate_polynomial + coeffs(i) * x**(i - 1)
    end do

  end function evaluate_polynomial


  real(q) function evaluate_derivative(coeffs, n, x)
    implicit none
    integer:: n
    real(q):: coeffs(n + 1), x
    integer:: i
    evaluate_derivative = 0.0_q

    do i = 2, n + 1
      evaluate_derivative = evaluate_derivative + (i - 1) * coeffs(i) * x**(i - 2)
    end do

  end function evaluate_derivative

  real(q) function evaluate_second_derivative(coeffs, n, x)
    implicit none
    integer:: n
    real(q):: coeffs(n + 1), x
    integer:: i
    evaluate_second_derivative = 0.0_q

    do i = 3, n + 1
      evaluate_second_derivative = evaluate_second_derivative + (i - 1) * (i - 2) * coeffs(i) * x**(i - 3)
    end do

  end function evaluate_second_derivative

      END MODULE
