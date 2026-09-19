#include "symbol.inc"
#include "bse.inc"
!#define bsedebug

!******************************************************************************
!
!> This module provides routines to set up the BSE Hamiltonian and solve
!> the BSE equation.
!
!******************************************************************************

MODULE bse_driver
  USE chi_glb
  USE local_field
  USE scala
  USE scala_struct_def
  USE wave_high
  USE bse_struct
  USE msupercell
  USE iso_c_binding
  USE bse_te
  USE cl
  USE core_con_mat
  USE incar_reader, ONLY: PROCESS_INCAR
  USE reader_tags, ONLY: INCAR_F
  USE bse_lanczos
#ifdef VASP_HDF5
    USE vhdf5
#endif
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CALCULATE_BSE_ACC

CONTAINS
!**********************************************************************
!
!> Solves the BSE equation.
!> 
!> It is based on the local_field subroutine and relies on the
!> routines implemented therein
!> the required integrals correspond to the publication
!> Sander, Maggio, Kresse, Phys. Rev. B 92, 045209
!>
!>  k1,n1=  i
!>  k2,n2=  j
!>  k3,n3=  a
!>  k4,n4=  b
!>
!> The  direct contribution (A = <e,h |W| e',h'>)
!> (two electrons at the same spatial coordinate, and two holes)
!>
!> ~~~
!>     int_d3r d3r' c*_k1+q,n3(r') c'_k2+q,n4(r') v'*_k2,n2(r)   v_k1,n1(r)  W(r',r)
!> ~~~
!>
!> and a second  exchange like contribution (B= <e,h | h',e'>)
!> (lacking in the Tamm-Dancoff approximation)
!> ~~~
!>    int_d3r d3r' c_k1-q,n3(r')  v'*_k2,n2(r') c'_k2+q,n4(r)  v*_k1,n1(r)  W(r',r)
!> ~~~
!>
!> and a term stemming from the variation of the Hartree potential
!>
!> ~~~
!>    int_d3r d3r' v'*_k2,n2(r') c'_k2+q,n4(r') c*_k1+q,n3(r) v_k1,n1(r) v(r',r)
!> ~~~
!>
!>
!> Note that the notation A and B refers to the quantum chemists notation,
!> where the full TDFT-TDHF response function is often written as
!>
!> ~~~
!>        K_ai,bj    K_ai,jb         A_ai,bj   B_ai,bj
!>                                =
!>        K_ia,bj    K_ia,jb         B*_ai,bj  A*_ai,bj
!> ~~~
!>
!> v states are restricted to valence band states
!> whereas c are conduction band states
!> the sums are restricted to NBANDSLF bands below and above the Fermi-level

!> loops over n1/n2/n3 and n4 are done in blocks of NSTRIP/NLOC
!>  where NSTRIP is the number of local bands, and NLOC the total number
!>
!> outline of the structure of the routine:
!>
!> ~~~
!>   allocate array W1 for merged block n1
!>   allocate array W2 for merged block n2
!>
!>   allocate array W3 for merged block n3
!>   allocate array W4 for merged block n4
!>   allocate charge array (WKAPPA) of size merged block n1 x NBANDSLF
!>     (this costs possibly a lot of memory and restricts n1)
!>
!>   allocate charge array of size merged WA block n4 x merged block n1
!>
!>   loop over all k1-points
!>     loop overal all v1-states in blocks (preferably all n1)
!>        gather wavefunction in the current block
!>        loop overall k2-points
!>           loop overal all v2-states in blocks (preferably all n2)
!>              gather wavefunction in the current block
!>
!>              calculate F(n1,n3,n2,n4) = c*_k1+q,n3  c'_k2+q,n4(r') W(r',r) v'*_k2,n2(r), v_k1,n1(r)
!> ~~~
!>
!>
!> ANTIRES = 0 corresponds to the Tamm-Dancoff approximation, only direct terms\n
!> ANTIRES = 1 applies a fancy approximation that includes 
!>             resonant-antiresonant coupling with w=0 (zero frequency) exact \n
!> ANTIRES = 2 full version
!
!
!
!**********************************************************************

  SUBROUTINE CALCULATE_BSE_ACC( &
       IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI,  &
       T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
       WGW, SHIFT,CELTOT_HK, CELTOT_X, &
       NBANDSO,NBANDSV, OMEGAMAX_IN, NQPOINT, KPOINT_BSE, ISP_IN, &
       LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS_IN, NELMGW, NBANDSGW, L2ORDER, &
       LCORBSE, LQP, LTCTC, CORREP, LAMBDA, S2E, EFERMI, SEFROMBSE, FIRST_EMPTYSTATE, LAST_OCCSTATE)

GPACC USE moffload

    USE base
    USE pseudo
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE constant
    USE poscar
    USE pot
    USE pawm
    USE kpoints_change
    USE full_kpoints
    USE mlr_optic
    USE dfast
    USE choleski
    USE ini
    USE gauss_quad
    USE mymath
    USE tutor, ONLY: vtutor, isAlert, isError, BSEantires, BSEncv, BSEnocder
    USE c2f_interface, ONLY : VTIME

#if defined(CUBLASMP) || defined(CUSOLVERMP)
    USE mmp_comm
    USE mmp_grid
#endif

    IMPLICIT NONE

    ! structures
    INTEGER              IBSE       !< 2 exact diagonalization (new driver), 1 TE-BSE (new driver), 3 LANCZOS-BSE (new driver)
    TYPE (type_info)     T_INFO
    TYPE (potcar)        P(T_INFO%NTYP)
    TYPE (wavedes)       WDES
    TYPE (nonlr_struct)  NONLR_S
    TYPE (nonl_struct)   NONL_S
    TYPE (wavespin)      W
    TYPE (latt)          LATT_CUR, LATT_INI
    TYPE (dynamics)      DYN
    TYPE (info_struct)   INFO
    TYPE (in_struct)     IO
    TYPE (kpoints_struct) KPOINTS
    TYPE (symmetry)      SYMM
    TYPE (grid_3d)       GRID       !< grid for wavefunctions
    INTEGER  LMDIM
    OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    TYPE (wavedes)       WGW        !< descriptor for basis set of response function
    REAL(q)              SHIFT      !< complex frequency shift for broadening
    COMPLEX(q) :: CELTOT_HK(:,:,:)  !< this contains the Hartree + kinetic term
    COMPLEX(q) :: CELTOT_X(:,:,:)   !< this contains the exchange term
    INTEGER              NBANDSO    !< number of occupied states included in BSE
    INTEGER              NBANDSV    !< number of virtual states included in BSE
    REAL(q)              OMEGAMAX_IN!< maximum considered transition energy in BSE
    INTEGER              NQPOINT    !< q-point for which excitations are calculated (usually optical e.g. q=0)
    INTEGER              KPOINT_BSE(3) ! reciprocal lattice vector G (g=G+q) at which response is evaluated
    INTEGER              ISP_IN     !< in: spin index to be considered:\n
                                    !<  ISP_IN<0: both spin channels are treated simultaneously (only valid for ISPIN=2)\n
                                    !<  ISP_IN=1 : treat only spin channel 1 (and pretend the system is non-magnetic)\n
                                    !<  ISP_IN=2 : treat only spin channel 2 (and pretend the system is non-magnetic) 
    LOGICAL              LHARTREE   !< include Hartree effects
    LOGICAL              LTRIPLET   !< calculate triplet
    LOGICAL              LADDER     !< local field effects beyond RPA (ladder diagrams)
    INTEGER              ANTIRES    !< how to treat antiresonant part \n
                                    !< presently implemented 0 = Tamm-Dancoff\n
                                    !< presently implemented 1 -> correct w=0 properties
    LOGICAL              LGWLF      !< use screened W from GW for exchange interactions
    LOGICAL              LFXC       !< include DFT xc kernel
    INTEGER              NEDOS_IN   !< number grid points in DOS (read from INCAR)
    INTEGER              NELMGW     !< number of timesteps in solving the BSE equation using time propagation
    INTEGER              NELM       !< NELMGW is also the number of iterations in
                                    !< the GW0Gamma0 QP energies evaluations
                                    !< reqquired here for the gammareal version
    INTEGER              NBANDSGW   !< number of bands for QP calculation
    LOGICAL              L2ORDER    !< only second order contribution
    LOGICAL              LCORBSE    !< determine correlation energy from BSE
    LOGICAL              LQP        !< determine self-energy
    LOGICAL              LTCTC      !< if true computes the test charge-test charge selfenergy
    REAL(q), OPTIONAL :: LAMBDA     !< coupling constant value can be used to scale the Coulumb integrals
    TYPE (screened_2e_handle), OPTIONAL :: S2E !< screened two electron integrals
    REAL(q), OPTIONAL :: EFERMI
    TYPE(selfenergy_from_bse), OPTIONAL :: SEFROMBSE
    INTEGER, OPTIONAL :: FIRST_EMPTYSTATE(:,:,:), LAST_OCCSTATE(:,:,:)

! local
    REAL(q)              OMEGAMAX   ! maximum considered transition energy in BSE
    INTEGER              NEDOS      ! number grid points in DOS
    TYPE(banddesc) :: BD(WDES%ISPIN)! descriptor for describing which bands are included, spin up
    TYPE(banddesc) :: B1,B2         ! current descriptor for the spin loops
    TYPE(bse_matrix_index) :: BSE_INDEX
    ! the following matrices store are sub range of the full matrices AMAT, BMAT and VMAT
    GDEF, POINTER :: TWOELECTRON3O(:,:,:,:)=> NULL()  ! part of the A matrix resonant-resonant coupling
    GDEF, POINTER :: BTWOELECTRON3O(:,:,:,:)=> NULL()  ! part of the A matrix resonant-resonant coupling

    GDEF, ALLOCATABLE, TARGET :: NTWOELECTRON3O(:,:,:,:,:)   ! part of the A matrix resonant-resonant coupling
    GDEF, ALLOCATABLE, TARGET :: NBTWOELECTRON3O(:,:,:,:,:)  ! part of the B matrix resonant-antiresonant coupling
    GDEF, ALLOCATABLE :: NTWOE3O_HARTREE(:,:,:,:,:)  ! part of Coulomb kernel in four orbital representation (Hartree only no exchange)
    TYPE (wavespin) WHF
    TYPE (wavefun1),ALLOCATABLE :: W1(:), W2(:), W3(:), W4(:)
    TYPE (wavefuna) :: WA1, WA2, WA3, WA4
    TYPE (wavedes1), TARGET :: WDESK1, WDESK2, WDESK3, WDESK4
    INTEGER :: NSTRIP               ! block size used internally in TWOELECTRON4O routines
    INTEGER :: NSTRIPV              ! block size for valence bands
    INTEGER :: NSTRIPC              ! block size for conduction bands
    INTEGER :: N, I, J, INDLR, INDLC, INDGR, INDGC, COL, ROW, M
    INTEGER :: NPOS1, NSTRIP1       ! base index and width of the n1 block
    INTEGER :: NPOS2, NSTRIP2       ! base index and width of the n1 block
    INTEGER :: NPOS3, NSTRIP3       !
    INTEGER :: NPOS4, NSTRIP4       !
    INTEGER :: K1, K2, K3, K4
    INTEGER :: ALLOC_VALENCE_BANDS, ALLOC_CONDUCTION_BANDS, ALLOC_KPOINTS
    REAL(q) :: NFFTW
    REAL(q) :: NFLOAT4O, NFFT4O     ! number of BLAS3 operations in 4 orbital routines
    LOGICAL :: LscaLAPACKaware      ! distribute data using BLACS
    LOGICAL :: LDISTRIBUTED_MATRIX  ! enable distributed matrix or store the full matrix on master rank
    INTEGER :: NCV                  ! number of conduction valence band pairs
    INTEGER :: IFAIL
    INTEGER :: ISP, ISP2, ISP_LOW, ISP2_LOW, ISP_HIGH
    LOGICAL :: W1EQUALW2, LKPOINT_PARALLEL, LBANDS_PARALLEL
    REAL(q) :: AHARTREE
    LOGICAL :: LFULL
    LOGICAL :: LTRIANGLE            ! only calculate lower triangle of A and B matrix and add remainder using conjugation
    REAL(q) :: POTFAK(GRIDHF%MPLWV)       ! 1/(G+dk)**2 (G)
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE, TARGET :: AMAT_SCALA(:)  ! A matrix used in BSE
    GDEFS, ALLOCATABLE, TARGET :: BMAT_SCALA(:)  ! B matrix used in BSE
    GDEFS, ALLOCATABLE :: OPTMAT(:,:)
    GDEFS :: DUMBSE                           ! dummy var of AMAT_SCALA type
    GDEFS :: TMP
#endif
#ifdef double_prec_bse
    GDEF, ALLOCATABLE, TARGET :: AMAT_SCALA(:)
    GDEF, ALLOCATABLE, TARGET :: BMAT_SCALA(:)
    GDEF, ALLOCATABLE :: OPTMAT(:,:)
    GDEF :: DUMBSE                           ! dummy var of AMAT_SCALA type
    GDEF :: TMP
#endif
    REAL(q), OPTIONAL :: CORREP       ! correlation energy from plasmon formula @ given k-point
    LOGICAL :: LALPHA       ! true if the coupling constant value LAMBDA has being passed by calling routine
    REAL(q) :: ALPHAV(10), W_GL(10) ! variables for the AC integration, max 10
    INTEGER :: NALPHA       ! number of points considered in coupling constant integration
    INTEGER :: IALPHA       ! current point considered in coupling constant integration
    REAL(q), PARAMETER :: SCALE_PERT=1._q ! this parameter can be used to scale the Coulomb interaction
    INTEGER(qi8) :: NSIZE   ! counter for allocation

    REAL(q), ALLOCATABLE :: R(:), RB(:)
! scaLAPACK
    INTEGER :: NP,NQ,NPROW,NPCOL,MYROW,MYCOL
    INTEGER, EXTERNAL    ::     NUMROC, INDXG2L
    INTEGER :: nd
#ifdef scaLAPACK
    TYPE(scalapack_des) :: BSE_GSD
#endif
    !INTEGER :: BSE_DESC(DLEN_)
    INTEGER :: ISTAT  ! allocation error code
    INTEGER(qi8) :: SIZE8,SIZE8A,SIZE8B,SIZE8C,K8 ! for giant arrays

    ! variables for BSE blocks
    INTEGER :: BLOCK_INDEX       ! block iterator
    INTEGER :: BLOCK_INDEX_START ! first block calculated in parallel
    INTEGER :: BLOCK_INDEX_END   ! last block calculated in parallel
    INTEGER :: BLOCK_INDEX_DIST  ! block iterator for data distribution
    INTEGER :: NBLOCKS
    INTEGER :: NODE
    INTEGER :: KPAR
    INTEGER :: PADDED_BLOCKS
    INTEGER :: LOCAL_BLOCK, GLOBAL_BLOCK
    INTEGER :: FIRST_GLOBAL_BLOCK, FIRST_LOCAL_BLOCK
    INTEGER :: NRANKS, RECT
    INTEGER(qi8) :: NBLOCKS_MAX8
    INTEGER :: NB
    GDEFS, ALLOCATABLE :: CDER_BETWEEN_STATES_TMP(:,:,:,:,:)
    INTEGER ::  SHAPE_CDER_BETWEEN_STATES(5)
    INTEGER :: NCPU, HSTRIP, HLLD
    REAL(q) :: MFREE
    LOGICAL :: LIDLE
    COMPLEX(q) :: CSUM_A    ! sum of diagonal components of matrix A
    INTEGER :: NCORE_MAX
    REAL(q)    :: QPOINT(3) !< value of q-point for which excitations are calculated (usually optical e.g. q=0)
    REAL(q)    :: VKPT(3)
    REAL(q)    :: VPUT0,CPUT0,VPUT,CPUT
    INTEGER    :: COLS


    GDEFS, ALLOCATABLE :: CDER_BETWEEN_STATES_BACKUP(:,:,:,:,:) ! to temporarily free memory GPU memory for diagonalisation (CDER_BETWEEN_STATES is "acc declare create", so host and gpu copy cannot be easily managed individually)

! include all code so that without scalapack the routine is empty
    PROFILING_START('calculate_bse')

    MULTIPLY_FRACTIONAL_OCC=.TRUE.
    LTRIANGLE=.TRUE.
    LBANDS_PARALLEL = .FALSE.
    LDISTRIBUTED_MATRIX = .TRUE.
    W1EQUALW2=.FALSE.
    LIDLE=.FALSE.
    LscaLAPACKaware=.TRUE.
    LKPOINT_PARALLEL=.FALSE.

#ifdef scaLAPACK
#ifdef VASP_HDF5
      CALL VH5_WRITE_EIGENVAL(IH5OUTFILEID, WDES, W, KPOINTS, KPOINTS_FULL)
#endif
!=======================================================================
! somewhat akward: the routine can solve the BSE equation
! at various Coulomb couplings in the list ALPHAV(1:NALPHA)
! this allows to perform a coupling constant integration to obtain
! the correlation energy
!=======================================================================
    LALPHA=.FALSE.              ! true if the coupling constant value LAMBDA has being passed by calling routine
    NALPHA=1
    W_GL=1.0; ALPHAV=1.0     ! weight 1 and unscaled Coulomb kernel
    MFREE = 0.0_q
    QPOINT = 0.0_q

    IF ( PRESENT(LAMBDA) ) THEN
       IF (LAMBDA /=1.0_q .AND. LAMBDA/=INT(LAMBDA)) THEN
          LALPHA=.TRUE.
          NALPHA=1
          W_GL=1.0; ALPHAV=LAMBDA ! weight 1, but scaled Coulomb kernel
       ENDIF
    ENDIF
    VKPT = 0.0_q

    ! smaller values are dangerous and are usually not what the user desires
    NEDOS=MAX(NEDOS_IN,1000)

    IF (ISP_IN<0 .AND. WDES%ISPIN==1) THEN
       CALL vtutor%bug("internal error in CALCULATE_BSE: ISP is negative but ISPIN=1", __FILE__, __LINE__)
    ENDIF

    IF (WDES%LGAMMA .AND. (.NOT. WGW%GRID%REAL2CPLX .OR. .NOT. WGW%GRID%LREAL)) THEN
       CALL vtutor%bug("internal error in CALCULATE_BSE: WGW must apply real to complex FFT", __FILE__, __LINE__)
    ENDIF
!=======================================================================
! preparation
!=======================================================================
    ! read in WPOT files
    ! safer here, since it can happen that not all nodes call GET_WPOT simultaneously
    ! (e.g. k-point parallelization)
    ! this might cause problems with SHMEM calls
    IF (LADDER .AND. LGWLF .AND. .NOT. ASSOCIATED(WPOTH)) THEN
       CALL INIT_WPOT_HANDLE( WPOTH, WGW, KPOINTS_FULL%NKPTS, IO, 1, 1, 1)
       DO K1=1,KPOINTS_FULL%NKPTS
          IF (W%WDES%WTKPT(K1)==0) CYCLE
! NKREDLF not presently used by bse.F
!          IF (SKIP_THIS_KPOINT_IN_LF(W%WDES%VKPT(:,K1), NKREDLFX, NKREDLFY, NKREDLFZ)) CYCLE
          K2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)-W%WDES%VKPT(:,1),KPOINTS_FULL)
          IF (ASSOCIATED(WPOTH)) CALL GET_WPOT( WPOTH, K2, POTFAK, LFULL)
       ENDDO
    ENDIF

    IF(W%WDES%NKPTS > 1) LKPOINT_PARALLEL = .TRUE.

    WHF=W      ! use temporarily another WDES
    WHF%WDES => WDES_FOCK


    NSTRIPV = 1
    NSTRIPC = 1

! determine wether the k1+k2 is found in the full k-point grid or not
    DO K1=1,KPOINTS_FULL%NKPTS
       IF (LIDENTICAL_KPOINT(WHF%WDES%VKPT(:,1)+WHF%WDES%VKPT(:,MIN(KPOINTS_FULL%NKPTS,2)),KPOINTS_FULL%VKPT(:,K1))) EXIT
    ENDDO

    ! prepare TRANS_MATRIX and derivatives (not ported)
    CALL CALCULATE_BSE_CORE_ORBITAL(P, IO%IU6, W, WHF, T_INFO )

#ifdef ACC_OFFLOAD
    ! switch on GPU execution
    CALL OFFLOAD_ON_PUSH(.TRUE.)
    ! if GPU execution is in principle allowed, OFFLOAD_ON=.TRUE. at this point
    IF (OFFLOAD_ON) THEN
#if ! defined(CUBLASMP) || ! defined(CUSOLVERMP)
       CALL vtutor%alert("CALCULATE_BSE_ACC: GPU execution without CUBLASMP and CUSOLVERMP&
                         & may be very inefficient!! Please consider recompiling VASP with&
                         & both these options enabled ...")
#endif
#if ! PGI_BEFORE_25_9 && defined(CUBLASMP) && defined(CUSOLVERMP)
       IF (.NOT.LUSENCCL) &
          CALL vtutor%error("CALCULATE_BSE_ACC: NCCL is switched off. BSE calculations on&
                            & GPU require NCCL to be active. Please make sure that:\n&
                            &    - VASP is compiled with NCCL support (-DUSENCCL).\n&
                            &    - You are running with a single mpi-rank per GPU.\n&
                            &    - You have not set LUSENCCL=.FALSE. in INCAR.")
#endif
    ENDIF
    CALL COPYIN_TYPED_VAR(WHF)
!$ACC ENTER DATA CREATE(WDESK1,WDESK2,WDESK3,WDESK4) __IF_ASYNC__
#endif


! set the wavefunction descriptors
    CALL SETWDES(WHF%WDES,WDESK1,0)
    CALL SETWDES(WHF%WDES,WDESK2,0)
    CALL SETWDES(WHF%WDES,WDESK3,0)
    CALL SETWDES(WHF%WDES,WDESK4,0)


    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,'("NQ=",I4,3F10.4,", ")') NQPOINT,WHF%WDES%VKPT(:,NQPOINT)
    ENDIF

    DO ISP=1,WDES%ISPIN
       CALL SET_BAND_PARAMETERS(W, ISP, NBANDSO, NBANDSV, BD(ISP), NSTRIPV, NSTRIPC)
    ENDDO

    ! Stop if too few conduction bands
    DO ISP=1,WDES%ISPIN
       IF(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1 < NBANDSV) THEN
          CALL vtutor%error("The number of conduction bands in the BSE kernel exeeds NBANDS." &
                            // "Reduce NBANDSV or redo the preceeding steps with larger NBANDS.")
       ENDIF
    ENDDO

    ! if ISP_IN is set, only the channel ISP_IN is considered
    IF (ISP_IN>0) THEN
       CALL SET_BAND_PARAMETERS(W, ISP_IN, NBANDSO, NBANDSV, BD(ISP_IN), NSTRIPV, NSTRIPC)
       ISP_LOW =ISP_IN
       ISP_HIGH=ISP_IN
    ELSE
       ! spin polarized calculations
       ISP_LOW =1
       ISP_HIGH=2
    ENDIF

    ! block over bands
    IF(NBSEBLOCKO > 0) THEN
      NSTRIPV = MIN(NSTRIPV,NBSEBLOCKO)
      LBANDS_PARALLEL = .TRUE.
    ENDIF

    ! block over bands
    IF(NBSEBLOCKV > 0) THEN
      NSTRIPC = MIN(NSTRIPC,NBSEBLOCKV)
      LBANDS_PARALLEL = .TRUE.
    ENDIF

    IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
       ! allocate space for all k-points in the full wedge
       CALL ALLOCATE_CDER_BETWEEN_STATES_Q( WHF%WDES%NB_TOT, MAX(BD(ISP_LOW)%VBMAX, BD(ISP_HIGH)%VBMAX), &
            WHF%WDES%NKPTS , WHF%WDES%ISPIN)
       IF (ICORELEVEL/=0) THEN
         NCORE_MAX = NITYP_CL * (2*L_CL+1)
         ALLOCATE(CDER_CORECON_Q(NCORE_MAX,WHF%WDES%NB_TOT,  WHF%WDES%NKPTS, WHF%WDES%ISPIN))
         CDER_CORECON_Q=0
      ENDIF
    ENDIF

    IF (LKPOINT_PARALLEL) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,*) 'parallelization over k-points'
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'parallelization over k-points'
    ENDIF
    IF (LBANDS_PARALLEL) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,*) 'parallelization over bands'
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'parallelization over bands'
    ENDIF


    ! allocate the orbitals: these contain the states in the conduction or valence bands
    ALLOCATE(W1(NSTRIPV))
    ALLOCATE(W2(NSTRIPV))
    ALLOCATE(W3(NSTRIPC))
    ALLOCATE(W4(NSTRIPC))

    !$ACC ENTER DATA CREATE(WA1, WA2, WA3, WA4) __IF_ASYNC__

    CALL NEWWAVA(WA1,WDESK1,NSTRIPV,1, .TRUE.)
    CALL REGISTER_ALLOCATE(16._q*SIZE(WA1%CPTWFP,KIND=qi8)+&
                            8._q*SIZE(WA1%CPROJ,KIND=qi8) +&
                            16._q*SIZE(WA1%CR,KIND=qi8), "wa")
    CALL NEWWAVA(WA2, WDESK2, NSTRIPV, 1, .TRUE.)
    CALL REGISTER_ALLOCATE(16._q*SIZE(WA2%CPTWFP,KIND=qi8)+&
                            8._q*SIZE(WA2%CPROJ,KIND=qi8) +&
                            16._q*SIZE(WA2%CR,KIND=qi8), "wa")

    CALL NEWWAVA(WA3,WDESK3,NSTRIPC,1, .TRUE.)
    CALL REGISTER_ALLOCATE(16._q*SIZE(WA3%CPTWFP,KIND=qi8)+&
                            8._q*SIZE(WA3%CPROJ,KIND=qi8) +&
                            16._q*SIZE(WA3%CR,KIND=qi8), "wa")
    CALL NEWWAVA(WA4, WDESK4, NSTRIPC, 1, .TRUE.)
    CALL REGISTER_ALLOCATE(16._q*SIZE(WA4%CPTWFP,KIND=qi8)+&
                            8._q*SIZE(WA4%CPROJ,KIND=qi8) +&
                            16._q*SIZE(WA4%CR,KIND=qi8), "wa")

    DO N=1,NSTRIPV
       W1(N)%CW => WA1%CW(:,N)
       W1(N)%CR => WA1%CR(:,N)
       W1(N)%CPROJ => WA1%CPROJ(:,N)
       W1(N)%WDES1 => WA1%WDES1
       W2(N)%CW => WA2%CW(:,N)
       W2(N)%CR => WA2%CR(:,N)
       W2(N)%CPROJ => WA2%CPROJ(:,N)
       W2(N)%WDES1 => WA2%WDES1
    ENDDO

    DO N=1,NSTRIPC
       W3(N)%CW => WA3%CW(:,N)
       W3(N)%CR => WA3%CR(:,N)
       W3(N)%CPROJ => WA3%CPROJ(:,N)
       W3(N)%WDES1 => WA3%WDES1
       W4(N)%CW => WA4%CW(:,N)
       W4(N)%CR => WA4%CR(:,N)
       W4(N)%CPROJ => WA4%CPROJ(:,N)
       W4(N)%WDES1 => WA4%WDES1
    ENDDO

#ifdef ACC_OFFLOAD
    CALL COPYIN_TYPED_VAR(W1)
    CALL COPYIN_TYPED_VAR(W2)
    CALL COPYIN_TYPED_VAR(W3)
    CALL COPYIN_TYPED_VAR(W4)
#endif

    ! set variable nd equal to ndata or twice ndata
    IF (ANTIRES >=2) THEN
       nd=2*ndata
    ELSE
       nd=ndata
    ENDIF


    OMEGAMAX=OMEGAMAX_IN
    IF (IO%IU6>=0) THEN
       DO ISP=ISP_LOW,ISP_HIGH
          WRITE(IO%IU6,'(/A,I2,/2(A,I5,2X,A,I5/))')' Bands included in the BSE spin=',ISP, &
               ' VB(min)=',BD(ISP)%VBMIN,'VB(max)=',BD(ISP)%VBMAX,' CB(min)=',BD(ISP)%CBMIN,'CB(max)=',BD(ISP)%CBMAX
       ENDDO

       IF (OMEGAMAX>0) THEN
          WRITE(IO%IU6,'(A,F10.5)') ' electron-hole pairs beyond OMEGAMAX are removed OMEGAMAX=',OMEGAMAX
       ENDIF
       IF (LADDER .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF) ) THEN
          IF (LGWLF) THEN
             WRITE(IO%IU6,'(A,3F14.7)')' W is read from the files WXXXX.tmp or WFULLXXXX.tmp if present'
          ELSE
             WRITE(IO%IU6,'(A,F14.7,A,F14.7)')' parameters for screened Coulomb W: AEXX=',AEXX,' HFSCREEN=',HFSCREEN
          ENDIF
          IF (.NOT. LHARTREE .OR. LTRIPLET) THEN
             WRITE(IO%IU6,'(A)') ' only ladders are included (Hartree term=RPA part is switched off)'
          ENDIF
       ELSE
          IF (.NOT. LHARTREE .OR. LTRIPLET) THEN
             WRITE(IO%IU6,'(A)') ' IP particle spectrum'
          ELSE
             WRITE(IO%IU6,'(A)')' simple RPA calculation, excitonic effects (ladders) are not included'
          ENDIF
       ENDIF
       IF (LHARTREE .AND. LFXC .AND. .NOT. LGWLF) THEN
          WRITE(IO%IU6,'(A)')' xc-kernel from LOCAL part of DFT is included (LFXC=.TRUE.)'
       ENDIF
    ENDIF

!==========================================================================
! few BSE specific things
!==========================================================================
    ALLOC_KPOINTS=0
    DO K1=1,WDES%NKPTS
       IF (WHF%WDES%WTKPT(K1)==0) CYCLE
       ALLOC_KPOINTS=ALLOC_KPOINTS+1
    ENDDO
    CALL SET_BSE_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, BSE_INDEX, WHF, OMEGAMAX, NCV, NQPOINT, NSTRIPV, NSTRIPC, .TRUE.)

    IF (NCV == 0) THEN
       CALL vtutor%write(isError, BSEncv)
    ENDIF
    IF (ALLOCATED(CDER_BETWEEN_STATES) .AND. (NQPOINT == 1)) THEN
       IF (MAXVAL(BD(:)%CBMAX) > SIZE(CDER_BETWEEN_STATES,1)) DEALLOCATE(CDER_BETWEEN_STATES)
    ENDIF
    IF ((NQPOINT == 1) .AND. .NOT. ALLOCATED(CDER_BETWEEN_STATES) .AND. .NOT. NBANDSO==0) THEN
          IF (IO%IDIOT > 0) THEN
             CALL vtutor%write(isError, BSEnocder)
          ELSE
             IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Proceeding without WAVEDER'
             IF (IO%IU6>=0) WRITE(IO%IU6,*) 'No WAVEDER is used'
          ENDIF
    ENDIF

#ifdef ACC_OFFLOAD
    ! use distributed matrix by default
    ! disable LDISTRIBUTED_MATRIX if all of floowing conditions satisfied
    !    1. ACC active
    !    2. exact diagonalization
    !    3. matrix fits on a single GPU
    IF (OFFLOAD_NUM_DEVICES>0 .AND. IBSE==2) THEN
       ! determine free memory on GPU (in bytes)
       ISTAT = GET_FREE_MEM_DEVICE(MFREE)
       CALL M_min_i(WHF%WDES%COMM, MFREE,1)
       ! at least 2*NCV*NCV to store the matrix and eigenvectors
       ! the eigensolver allocates extra memory, we assume NCV*NCV
       SIZE8A = 4*NCV
       SIZE8B = SIZEOF(DUMBSE)*NCV
       SIZE8 = SIZE8A*SIZE8B
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F10.2)') " Memory available on GPU (Gb)", MFREE*1E-9
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F10.2)') " Memory required to diagonalize matirx on GPU (Gb)", SIZE8*1E-9
       IF(MFREE > REAL(SIZE8,q)) THEN
         LDISTRIBUTED_MATRIX = .FALSE.
       ELSE
         LDISTRIBUTED_MATRIX = .TRUE.
       ENDIF
    ENDIF
#endif

#ifndef CUBLASMP
    ! need to make sure that PGGEMR2D/PGGEMR2D 
    ! doesn't exceed int32
    IF(NCV > 11520) THEN
       LDISTRIBUTED_MATRIX = .TRUE.
    ENDIF
    SIZE8A=NSTRIPC*NSTRIPV ! cv pairs
    IF(SIZE8A > 5760) THEN
       CALL vtutor%alert("The PGGEMR2D/PGGEMR2D rouintes are likely to exceed int32 internal counter" &
                             //" which might cause an error 'xxmr2d:out of memory'. To avoid this error use " &
                            // " tags NBSEBLOCKV and/or NBSEBLOCKO to divide the work into smaller blocks.")
    ENDIF
#endif

    CALL PROCESS_INCAR(INCAR_F, "LSCAAWARE", LscaLAPACKaware, ISTAT)
    IF (ISTAT==0) LDISTRIBUTED_MATRIX=LscaLAPACKaware
    IF(LDISTRIBUTED_MATRIX) THEN
      IF (IO%IU0>=0) WRITE(IO%IU0,*) "Using distributed matrix"
    ELSE
      IF (IO%IU0>=0) WRITE(IO%IU0,*) "Using non-distributed matrix"
    ENDIF

    CALL PROCESS_INCAR(INCAR_F, "LSCAAWARE", LscaLAPACKaware, ISTAT)
    IF (ISTAT==0) LDISTRIBUTED_MATRIX=LscaLAPACKaware
    IF(LDISTRIBUTED_MATRIX) THEN
      IF (IO%IU0>=0) WRITE(IO%IU0,*) "Using distributed eigensolver"
    ELSE
      IF (IO%IU0>=0) WRITE(IO%IU0,*) "Using non-distributed eigensolver"
    ENDIF

    IF(LDISTRIBUTED_MATRIX) THEN
       CALL INIT_scala_DESC( WDES%COMM, NCV, BSE_DESC, BSE_GSD )
    ELSE
       ! this call only needed to fetch NB value, which is private in scala module
       CALL INIT_scala_DESC( WDES%COMM, NCV, BSE_DESC, BSE_GSD )
       NB = BSE_DESC(NB_)

       CALL BLACS_GRIDEXIT(BSE_DESC(CTXT_))
       ! create new contex with a single rank
       CALL BLACS_GET(0, 0, BSE_DESC(CTXT_))
       CALL BLACS_GRIDINIT(BSE_DESC(CTXT_), 'C',  1, 1)
       CALL BLACS_GRIDINFO(BSE_DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)
       IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
          CALL DESCINIT(BSE_DESC, NCV, NCV, NCV, NCV, 0, 0, BSE_DESC(CTXT_), NCV, IFAIL)
       ELSE
          BSE_DESC(CTXT_)=-1
       ENDIF
    ENDIF

    ! store stripes of BSEMATRIX on different processors
    IF (IBSE==1 .OR. IBSE==3) THEN
       CALL BLACS_GET(0, 0, BSE_DESC(CTXT_))
       CALL BLACS_GRIDINIT(BSE_DESC(CTXT_), 'C',  WDES%COMM%NCPU, 1)
       CALL BLACS_GRIDINFO(BSE_DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)

       ! use block size defined in scala, i.e., default or NBSCALAPACK
       NB = BSE_DESC(NB_)
       HSTRIP = BSE_DESC(NB_)
       HLLD = MAX(1,NUMROC(NCV,HSTRIP,MYROW,0,NPROW))

       IF(ANTIRES>=2) THEN
          CALL DESCINIT(BSE_DESC, NCV, NCV*2, HSTRIP, NB, 0, 0, BSE_DESC(CTXT_), HLLD, IFAIL)
       ELSE
          CALL DESCINIT(BSE_DESC, NCV, NCV, HSTRIP, NB, 0, 0, BSE_DESC(CTXT_), HLLD, IFAIL)
       ENDIF
       ! TE-BSE propagates 3 vectors simultaneously
       IF(IBSE ==1) COLS = 3
       ! Lanczos propagates 6 vectors simultaneously
       IF(IBSE ==3) COLS = 6

       ! set up descriptor for vectors
       IF(ANTIRES>=2) THEN
          HLLD= MAX(1,NUMROC(NCV*2,NCV*2,MYROW,0,NPROW))
          CALL DESCINIT(VEC_DESC, NCV*2, COLS, NCV*2, COLS, 0, 0, BSE_DESC(CTXT_), HLLD, IFAIL)

          HLLD= MAX(1,NUMROC(NCV,NB,MYROW,0,NPROW))
          CALL DESCINIT(DVEC_DESC, NCV, COLS, NB, COLS, 0, 0, BSE_DESC(CTXT_), HLLD, IFAIL)
       ELSE
          HLLD= MAX(1,NUMROC(NCV,NCV,MYROW,0,NPROW))
          CALL DESCINIT(VEC_DESC, NCV, COLS, NCV, COLS, 0, 0, BSE_DESC(CTXT_), HLLD, IFAIL)

          HLLD= MAX(1,NUMROC(NCV,NB,MYROW,0,NPROW))
          CALL DESCINIT(DVEC_DESC, NCV, COLS, NB, COLS, 0, 0, BSE_DESC(CTXT_), HLLD, IFAIL)
       ENDIF
    ENDIF

#ifdef bsedebug
    IF (IO%IU0>=0) WRITE(IO%IU0,'(A,9I8)') ' BSE descriptor', BSE_DESC
    IF (IO%IU6>=0) WRITE(IO%IU6,'(A,9I8)') ' BSE descriptor', BSE_DESC
#endif
!==========================================================================
! intialize scaLAPACK with appropriate matrix dimensions
!==========================================================================
    CALL BLACS_GRIDINFO(BSE_DESC(CTXT_), NPROW, NPCOL,MYROW,MYCOL)
    SIZE8A = NUMROC(BSE_DESC(M_),BSE_DESC(MB_),MYROW,0,NPROW)
    SIZE8B = NUMROC(BSE_DESC(N_),BSE_DESC(NB_),MYCOL,0,NPCOL)
    SIZE8=MAX(SIZE8A*SIZE8B,1_qi8)

    IF(SIZE8 > 2_qi8**31) THEN
       CALL vtutor%alert("The local number of elements ("// str(SIZE8)//") in the distributed BSE matrix " &
       // "exceeds 2**31 (integer 4). This is likely to cause a segmentation fault " // &
       "in the ScaLAPACK distribution routine. To avoid this issue increase " // &
       "the number of MPI ranks." )
    ENDIF
    IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
       ALLOCATE( AMAT_SCALA(SIZE8), STAT=ISTAT)
    ELSE
       ALLOCATE(AMAT_SCALA(1), STAT=ISTAT)
    ENDIF
    CALL CHECK_ALLOCATION_ERROR( "AMAT_SCALA", SIZE8, ISTAT )
    IF(IBSE==2) THEN
      ALLOCATE(R(BSE_DESC(M_)))
    ELSE
      ALLOCATE(OPTMAT(BSE_DESC(M_),3))
      ALLOCATE(R(1))
    ENDIF
    NSIZE=SIZE(AMAT_SCALA, KIND=qi8)
    
    IF (ANTIRES>=2 .AND. IBSE==2) THEN
       IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
           ALLOCATE(BMAT_SCALA(SIZE8), STAT=ISTAT)
       ELSE
           ALLOCATE(BMAT_SCALA(1), STAT=ISTAT)
       ENDIF
       CALL CHECK_ALLOCATION_ERROR( "BMAT_SCALA", SIZE8, ISTAT )
       NSIZE=NSIZE+SIZE(BMAT_SCALA, KIND=qi8)
       ALLOCATE(RB(BSE_DESC(M_)))
    ELSE
       ALLOCATE(RB(1))
       ALLOCATE(BMAT_SCALA(1))
    ENDIF
#ifdef single_prec_bse
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) single prec attempting allocation of',4.E-9_q*ndata8*NSIZE,' Gbyte  rank=',NCV
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) single prec attempting allocation of',4.E-9_q*ndata8*NSIZE,' Gbyte  rank=',NCV
       CALL REGISTER_ALLOCATE(4._q*ndata8* NSIZE, "bse")
#endif

#ifdef double_prec_bse
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) attempting allocation of',8.E-9_q*ndata8*NSIZE,' Gbyte  rank=',NCV
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) attempting allocation of',8.E-9_q*ndata8*NSIZE,' Gbyte  rank=',NCV
       CALL REGISTER_ALLOCATE(8._q*ndata8* NSIZE, "bse")
#endif
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,I12)') ' BSE (scaLAPACK) number of blocks=', BSE_INDEX%BLOCKS

!==========================================================================

    IF (IO%IU0>=0) WRITE(IO%IU0,'(A,4I5)') ' allocating two-electron 4 orbital integral table',NSTRIPV, NSTRIPC, NSTRIPV, NSTRIPC
    SIZE8A=NSTRIPC*NSTRIPV
    SIZE8B=NSTRIPC*NSTRIPV
! "SIZE8" replaces expression "NSTRIPV*NSTRIPC*NSTRIPV*NSTRIPC" of the original code
    SIZE8=SIZE8A*SIZE8B

    ALLOCATE(NTWOELECTRON3O( NSTRIPV, NSTRIPC, NSTRIPV, NSTRIPC, 1), STAT = ISTAT)
    CALL CHECK_ALLOCATION_ERROR_Q( "TWOELECTRON3O",SIZE8 , ISTAT )

    CALL REGISTER_ALLOCATE(8._q*ndata8* SIZE(TWOELECTRON3O, KIND=qi8) , "bse")
    IF (ANTIRES>=2) THEN
       ALLOCATE(NBTWOELECTRON3O( NSTRIPV, NSTRIPC, NSTRIPV, NSTRIPC, 1), STAT = ISTAT)
       CALL CHECK_ALLOCATION_ERROR_Q( "BTWOELECTRON3O",SIZE8 , ISTAT )

       CALL REGISTER_ALLOCATE(8._q*ndata8* SIZE(BTWOELECTRON3O, KIND=qi8) , "bse")
    ENDIF
!==========================================================================
    CALL DUMP_ALLOCATE(IO%IU6)
    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
!==========================================================================
! start cycle over values of the interaction strength in order
! to evaluate the correlation energy within the adiabatic-connection
! formalism.
!==========================================================================
#ifdef ACC_OFFLOAD
    ACC_ASYNC_Q=ACC_ASYNC_ASYNC
    CALL COPYIN_TYPED_VAR(KPOINTS_FULL)
    CALL COPYIN_TYPED_VAR(KPOINTS_FULL_ORIG)
    CALL COPYIN_TYPED_VAR(AUG_DES)
    CALL COPYIN_TYPED_VAR(FAST_AUG_FOCK)
    CALL COPYIN_TYPED_VAR(GRIDHF)
    !$ACC ENTER DATA COPYIN(TRANS_MATRIX_FOCK,LATT_CUR,BSE_INDEX) __IF_ASYNC__
    !$ACC ENTER DATA COPYIN(BSE_INDEX%INDEX,BSE_INDEX%N1,BSE_INDEX%N3,BSE_INDEX%NK,BSE_INDEX%ISP,BSE_INDEX%BK1, &
    !$ACC                   BSE_INDEX%BK2,BSE_INDEX%BISP1,BSE_INDEX%BISP2,BSE_INDEX%CK1,BSE_INDEX%CISP1, &
    !$ACC                   BSE_INDEX%CNPOS1,BSE_INDEX%CNPOS3,BSE_INDEX%BNPOS1,BSE_INDEX%BNPOS2,BSE_INDEX%BNPOS3, &
    !$ACC                   BSE_INDEX%BNPOS4,BSE_INDEX%BMAP,BSE_INDEX%MASK) __IF_ASYNC__
    !$ACC ENTER DATA COPYIN(BD,BSE_DESC,VEC_DESC,DVEC_DESC) __IF_ASYNC__
    !$ACC ENTER DATA CREATE(R) __IF_ASYNC__
    !$ACC ENTER DATA CREATE(AMAT_SCALA) ASYNC(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.ALLOCATED(AMAT_SCALA))
    !$ACC ENTER DATA CREATE(BMAT_SCALA) ASYNC(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.ALLOCATED(BMAT_SCALA))
    !$ACC ENTER DATA CREATE(RB)         ASYNC(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.ALLOCATED(RB))

    DO I=1,4
    DO ISP=ISP_LOW,ISP_HIGH
       CALL FFT_OFFLOAD_MAKEPLAN_BSE(GRIDHF,BD(ISP)%VBMIN,BD(ISP)%VBMAXpCORE,NSTRIPV,I)
       CALL FFT_OFFLOAD_MAKEPLAN_BSE(WGW%GRID,BD(ISP)%VBMIN,BD(ISP)%VBMAXpCORE,NSTRIPV,I)

       CALL FFT_OFFLOAD_MAKEPLAN_BSE(GRIDHF,BD(ISP)%CBMIN,BD(ISP)%CBMAX,NSTRIPC,I)
       CALL FFT_OFFLOAD_MAKEPLAN_BSE(WGW%GRID,BD(ISP)%CBMIN,BD(ISP)%CBMAX,NSTRIPC,I)

    ENDDO
    CALL FFT_OFFLOAD_MAKEPLAN(WDESK1%GRID,I,NSTRIPV,WDESK1%NRSPINORS*WDESK1%GRID%MPLWV)
    CALL FFT_OFFLOAD_MAKEPLAN(WDESK1%GRID,I,NSTRIPC,WDESK1%NRSPINORS*WDESK1%GRID%MPLWV)
    ENDDO
#endif
!==========================================================================
! loop over coupling constant integration points alpha
 ccintegration: DO IALPHA=1,NALPHA
!==========================================================================
    NFFTW = 0
    NFLOAT4O=0 ; NFFT4O=0

!$ACC KERNELS PRESENT(AMAT_SCALA) __IF_ASYNC__
    AMAT_SCALA=0
!$ACC END KERNELS
    IF (ANTIRES>=2) THEN
!$ACC KERNELS PRESENT(BMAT_SCALA) __IF_ASYNC__
       BMAT_SCALA=0
!$ACC END KERNELS
    ENDIF

    ! initialize the upper half of the matrix to the MAGIC number
    IF (IBSE==2) THEN
       CALL PGLASETbse('U', BSE_DESC(M_), BSE_DESC(N_), onebse*MAGIC, onebse*MAGIC, AMAT_SCALA(1), 1, 1, BSE_DESC)
       !$ACC UPDATE DEVICE(AMAT_SCALA) __IF_ASYNC__
       IF(ANTIRES == 2) THEN
         CALL PGLASETbse('U', BSE_DESC(M_), BSE_DESC(N_), onebse*MAGIC, onebse*MAGIC, BMAT_SCALA(1), 1, 1, BSE_DESC )
         !$ACC UPDATE DEVICE(BMAT_SCALA) __IF_ASYNC__
       ENDIF
    ENDIF
    ! initialize matrix to MAGIC number
    IF ((IBSE==1 .OR. IBSE==3) .AND. SIZE(AMAT_SCALA)> 1) THEN
    !$ACC KERNELS PRESENT(AMAT_SCALA) __IF_ASYNC__
    AMAT_SCALA = MAGIC
    !$ACC END KERNELS
    ENDIF

    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE setting up matrix'
    CALL START_TIMING("G")

    ! Tamm-Dancoff with correct properties at large w

    IF (ISP_IN<0 .OR. WHF%WDES%LNONCOLLINEAR ) THEN
       ! spin polarized calculation weight for Hartree term is one
       AHARTREE=1.0_q
    ELSE
       ! non spin polarized: weight for Hartree is two
       AHARTREE=2.0_q
    ENDIF

    ! correct properties at w=0
    ! at w=0  resonant-resonant coupling = antires-antires=resonant-antires
    ! so Hartree needs to be counted twice (this however spoils high frequency behaviour)
    IF (ANTIRES==1) AHARTREE=AHARTREE*2.0_q
    !==========================================================================
    ! outer loop body for first set of pair states
    !==========================================================================

    !$ACC ENTER DATA CREATE(NTWOELECTRON3O,NBTWOELECTRON3O) __IF_ASYNC__
    !$ACC KERNELS PRESENT(NTWOELECTRON3O) __IF_ASYNC__
    NTWOELECTRON3O = 0._q
    !$ACC END KERNELS

    IF(ANTIRES>=2) THEN
       !$ACC KERNELS PRESENT(NBTWOELECTRON3O) __IF_ASYNC__
       NBTWOELECTRON3O = 0._q
       !$ACC END KERNELS
    ENDIF
    FIRST_GLOBAL_BLOCK = 1
    NSTRIP1=NSTRIPV
    NSTRIP2=NSTRIPV
    NSTRIP3=NSTRIPC
    NSTRIP4=NSTRIPC

    ! PADDED_BLOCKS is number of blocks + blocks to complete the last row in even matrices
    ! |11|12|13|14|
    ! |44|22|23|24|
    ! |34|33|xx|xx|
    PADDED_BLOCKS = CEILING(REAL(BSE_INDEX%COLUMNS+1,KIND=q)/2._q)*BSE_INDEX%COLUMNS

    K1 = 1
    K2 = 1
    ISP = 1
    ISP2 = 1
    B1 = BD(1)
    B2 = BD(1)

    ! autoset available memory, in case MAXMEM if not set by user in INCAR
    IF ( LAUTOSET_MAXMEM ) THEN
       CALL AUTOSET_AVAILABLE_MEMORY( MAXMEM, IO )
#ifdef ACC_OFFLOAD
    IF (OFFLOAD_NUM_DEVICES>0) THEN
       ! check avail memeory again after allocating AMAT
       ISTAT = GET_FREE_MEM_DEVICE(MFREE)
       ! if avail memory on GPU is determined (non-zero)
       IF(ISTAT == 0) THEN
          MAXMEM = MIN(MAXMEM,INT(MFREE*1E-6))
          CALL M_min_i(WHF%WDES%COMM, MAXMEM,1)
       ENDIF
    ENDIF
#endif
    ENDIF

    IF(ALLOCATED(CDER_BETWEEN_STATES_Q)) THEN
       !$ACC KERNELS PRESENT(CDER_BETWEEN_STATES_Q) __IF_ASYNC__
       CDER_BETWEEN_STATES_Q=0
       !$ACC END KERNELS
    ENDIF

    ! memory estimator for HARTREE and DIRECT routines
    ! probably not all large arrays are included, but these should be by far the largest
    SIZE8A=NSTRIPC*NSTRIPV ! cv pairs
    SIZE8B=NSTRIPC*NSTRIPC ! cc pairs
    SIZE8C=NSTRIPV*NSTRIPV ! vv pairs
    IF(LHARTREE) THEN
       SIZE8 = 0
       ! CHAM(NSTRIPV, NSTRIPC, NSTRIPV,  NSTRIPC)
       SIZE8 = SIZE8 + SIZE8A*SIZE8A
       ! GCHG13(NP, NSTRIP1, NSTRIP3) and GCHG24
       SIZE8 = SIZE8 + MAXVAL(WGW%NGVECTOR)*SIZE8A*2
       ! CWORK(MAX(GRIDHF%MPLWV,WGW%GRID%MPLWV),MAX(NSTRIP1,NSTRIP2))
       SIZE8 = SIZE8 + GRIDHF%MPLWV*MAX(NSTRIPV,NSTRIPC)

       SIZE8A = SIZE8
    ENDIF
    IF(LADDER) THEN
       SIZE8 = 0
       ! CHAM(NSTRIPV, NSTRIPV, NSTRIPC,  NSTRIPC)
       SIZE8 = SIZE8 + SIZE8B*SIZE8C
       !GCHG12(NP, NSTRIPV, NSTRIPV), &
       SIZE8 = SIZE8 + MAXVAL(WGW%NGVECTOR)*SIZE8C
       !GCHG34(NP, NSTRIPC, NSTRIPC), &
       SIZE8 = SIZE8 + MAXVAL(WGW%NGVECTOR)*SIZE8B
       !GCHG(NP, NSTRIP1, NSTRIP2), &
       SIZE8 = SIZE8 + MAXVAL(WGW%NGVECTOR)*SIZE8C
       !CWORK12(MAX(GRIDHF%MPLWV,WGW%GRID%MPLWV),MAX(NSTRIP1,NSTRIP2))
       SIZE8 = SIZE8 + GRIDHF%MPLWV*MAX(NSTRIPV,NSTRIPV)
       !CWORK34(MAX(GRIDHF%MPLWV,WGW%GRID%MPLWV),MAX(NSTRIP3,NSTRIP4))
       SIZE8 = SIZE8 + GRIDHF%MPLWV*MAX(NSTRIPC,NSTRIPC)

       SIZE8B = SIZE8
    ENDIF
    SIZE8 = MAX(SIZE8A,SIZE8B)


    SIZE8A=NSTRIPC*NSTRIPV ! cv pairs
    SIZE8B=NSTRIPC*NSTRIPC ! cc pairs
    SIZE8C=NSTRIPV*NSTRIPV ! vv pairs

    MAXMEM = MAX(1,CEILING(MAXMEM - m_ 1E-6_q*SIZE8*8))

    ! estimate maximum number of blocks per node
#ifdef double_prec_bse
    ! MAXMEM in Mb
    IF(SUM(BSE_INDEX%MASK) /= SIZE(BSE_INDEX%MASK)) THEN
       NBLOCKS_MAX8=FLOOR((MAXMEM*0.9_q/(SIZE8A* SIZE8A* (16*2)))*1000000,KIND=qi8)
    ELSE
       NBLOCKS_MAX8=FLOOR((MAXMEM*0.9_q/(SIZE8A* SIZE8A* (16)))*1000000,KIND=qi8)
    ENDIF
#else
    NBLOCKS_MAX8=FLOOR((MAXMEM*0.9_q/(SIZE8A* SIZE8A* (16+16)))*1000000,KIND=qi8)
#endif

    ! might be out of memeory should give  a warning
    IF(NBLOCKS_MAX8 < 1_qi8 ) THEN
       CALL vtutor%alert("The availible memory might not be sufficient to perform this calculation. " &
                            // "Use tags NBSEBLOCKV and/or NBSEBLOCKO to divide the work into smaller blocks.")
       NBLOCKS_MAX8=1
    ENDIF

    ! make sure that we don't excceed INTEGER(4) when distributing data in PXGEMR2D/PXTRMR2D
    NBLOCKS_MAX8 = MIN(NBLOCKS_MAX8, MAX(1,2_8**31/(SIZE8B * SIZE8C * WDES%COMM%NCPU)))

#ifdef bsedebug
    IF (IO%IU0>0) WRITE(IO%IU0,*) "NBLOCKS_MAX8", NBLOCKS_MAX8
#endif

    SIZE8 = SIZE8 + MIN(BSE_INDEX%COLUMNS,NBLOCKS_MAX8)*(NSTRIPV* NSTRIPC* NSTRIPV* NSTRIPC)
    ! account for additional memory required for transposed blocks
    IF(.NOT. LTRIANGLE) SIZE8 = SIZE8 + MIN(BSE_INDEX%COLUMNS,NBLOCKS_MAX8)*(NSTRIPV* NSTRIPC* NSTRIPV* NSTRIPC)


    IF((MAXMEM - m_ 1E-6_q*SIZE8*8) < 0.0_q) THEN
       CALL vtutor%alert("The availible memory might be insufficient to perform this calculation. " &
                         // "Use tags NBSEBLOCKV and/or NBSEBLOCKO to divide the work into smaller blocks.")
    ENDIF


    IF (BSE_INDEX%BLOCKS < WDES%COMM%NCPU) THEN
       CALL vtutor%alert("Workload cannot be divided evenly."&
             // " To avoid load imbalance, use tags NBSEBLOCKV and NBSEBLOCKO to divide the work into smaller blocks.")
    ENDIF

#ifdef ACC_OFFLOAD
    IF(OFFLOAD_NUM_DEVICES>0 .AND. SUM(BSE_INDEX%MASK) /= SIZE(BSE_INDEX%MASK)) THEN
          CALL vtutor%alert("This calculation on GPU requires additional data transfers between &
          GPU and CPU, which leads to a serious performance degradation. Choose OMEGAMAX larger than the highest transition energy.")
    ENDIF
#endif
    IF (SUM(BSE_INDEX%MASK) < 0.5_q*SIZE(BSE_INDEX%MASK)) THEN
        IF(IO%IU6>=0) WRITE(IO%IU6,'(A,I8)') " skipping rows/columns:", SIZE(BSE_INDEX%MASK)-SUM(BSE_INDEX%MASK)
    ENDIF

    IF(IBSE == 1 .OR. IBSE == 3 .OR. ANTIRES == 2 .OR. LELPA) THEN
      LTRIANGLE=.FALSE.
    ENDIF

    IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
       DO I=1,3
          QPOINT(I)=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(I,1)+ &
             (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(I,2)+ &
             (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(I,3)
          VKPT(I)=WHF%WDES%VKPT(I,NQPOINT)+KPOINT_BSE(I)
       ENDDO
    ENDIF
#ifdef VASP_HDF5
    CALL VH5_WRITE_DIELECTRIC_VECTOR(IH5OUTFILEID,'results/linear_response', VKPT, 'bse_wave')
    CALL VH5_WRITE_DIELECTRIC_VECTOR(IH5OUTFILEID,'results/linear_response', QPOINT, 'bse_qpoint')
#endif
    IF (IO%IU0>=0) THEN
       WRITE(*,*) "Calculate matrix elements of the BSE Hamiltonian:"
       WRITE(*,*) "Workload is divided into blocks of two-electron 4 orbital integrals computed in parallel."
       WRITE(*,'(A,I4,",",I4,",",I4,",",I4,A)') " Each block/integral has dimensions  of (",NSTRIPV,NSTRIPC,NSTRIPV,NSTRIPC,") orbitals"
       WRITE(*,'(A,I10)') " Total number of blocks/integrals:           ", BSE_INDEX%BLOCKS
       WRITE(*,'(A,I2)')  " Number of spins:                                    ", WDES%ISPIN
       WRITE(*,'(A,I8)')  " Number of kpoints:                            ", KPOINTS_FULL%NKPTS
       DO ISP=ISP_LOW,ISP_HIGH
          WRITE(*,'(A,I2)')  " Spin channel:                                       ", ISP
          WRITE(*,'(A,I8)')  "  Number of blocks of unoccupied bands:        ", CEILING(REAL(BD(ISP)%CBMAX-BD(ISP)%CBMIN,KIND=q)/REAL(NSTRIPC,KIND=q))
          WRITE(*,'(A,I8)')  "  Number of blocks of occupied bands:          ", CEILING(REAL(BD(ISP)%VBMAXpCORE-BD(ISP)%VBMIN,KIND=q)/REAL(NSTRIPV,KIND=q))
       ENDDO
       WRITE(*,'(A,I8)')  " Number of rows/columns of blocks:             ", BSE_INDEX%COLUMNS

       WRITE(*,*) ""
       WRITE(*,*) "                                                       Time (sec)"
       WRITE(*,*) "              Integrals                 Ranks       Real         CPU"
    ENDIF
    IF (IO%IU0>=0) CALL VTIME(VPUT0,CPUT0)

    ! global loop over blocks
    globalblock: DO WHILE (FIRST_GLOBAL_BLOCK <= BSE_INDEX%BLOCKS)

       ! find maximum number of blocks that can be distributed on all ranks evenly
       NBLOCKS = CALC_NBLOCKS_ROUND_BLOCKS(FIRST_GLOBAL_BLOCK,PADDED_BLOCKS,WDES%COMM%NCPU,BSE_INDEX%COLUMNS,NBLOCKS_MAX8)


       ! we cannot distribute multiple rows if one of the rows is incomplete, i.e. only submatrix
       IF(WDES%COMM%NCPU*NBLOCKS > BSE_INDEX%COLUMNS) THEN
         RECT = FLOOR(REAL(WDES%COMM%NCPU*NBLOCKS,q)/REAL(BSE_INDEX%COLUMNS,q))*BSE_INDEX%COLUMNS
       ELSE
         RECT = WDES%COMM%NCPU*NBLOCKS
       ENDIF

       ! might need to reduce the number of nodes
       ! for example 8x4 matrix with 9 nodes will distribute 4 blocks per 8 nodes,
       ! leaving one node empty, i.e., need to truncate the last node, hence NRANKS = 8

       ! number of nodes to comptue the submatrix and make sure that we don't exceed the padded matrix size
       NRANKS = MIN(PADDED_BLOCKS,RECT) / NBLOCKS

       ! incomplete row truncate number of nodes if overflow to the next row
       ! for example, distribute 9-blocks row with 2 nodes, each has 2 blocks, in the third iteration only one node is used
       IF(MOD(FIRST_GLOBAL_BLOCK-1,BSE_INDEX%COLUMNS) /= 0) THEN
          ! ranks to compute the remaining number of blocks in the row
          NRANKS=MIN(WDES%COMM%NCPU,(BSE_INDEX%COLUMNS-MOD(FIRST_GLOBAL_BLOCK-1,BSE_INDEX%COLUMNS))/NBLOCKS)
       ENDIF

       ! special treatment for the last (incomplete) row in even matirces or trunacte if we exceed the number of blocks
       ! here we don't need padding
       IF (FIRST_GLOBAL_BLOCK > FLOOR((BSE_INDEX%COLUMNS+1)/2._q)*BSE_INDEX%COLUMNS .OR. FIRST_GLOBAL_BLOCK + NRANKS*NBLOCKS > BSE_INDEX%BLOCKS ) THEN
           ! reduce the number of ranks so that the remaining number of blocks (padded) can be distributed evenly
           NRANKS = CALC_NBLOCKS_ROUND_CPUS(PADDED_BLOCKS - FIRST_GLOBAL_BLOCK + 1,WDES%COMM%NCPU)
           ! make sure that the number of block does not exeed the limit
           NBLOCKS = MIN(NBLOCKS_MAX8,(PADDED_BLOCKS - FIRST_GLOBAL_BLOCK + 1)/NRANKS)
       ENDIF

       ! reallocate if the block
       IF(SIZE(NTWOELECTRON3O,5) /= NBLOCKS) THEN
         !$ACC EXIT DATA DELETE(NTWOELECTRON3O) __IF_ASYNC__
         IF(ALLOCATED(NTWOELECTRON3O)) DEALLOCATE(NTWOELECTRON3O)
         ALLOCATE(NTWOELECTRON3O( NSTRIPV, NSTRIPC, NSTRIPV, NSTRIPC, NBLOCKS), STAT = ISTAT)
         !$ACC ENTER DATA CREATE(NTWOELECTRON3O) __IF_ASYNC__
         !$ACC KERNELS PRESENT(NTWOELECTRON3O) __IF_ASYNC__
         NTWOELECTRON3O = 0._q
         !$ACC END KERNELS
       ENDIF

       IF(ANTIRES>=2) THEN
          IF(SIZE(NBTWOELECTRON3O,5) /= NBLOCKS) THEN
            !$ACC EXIT DATA DELETE(NBTWOELECTRON3O) __IF_ASYNC__
            IF(ALLOCATED(NBTWOELECTRON3O)) DEALLOCATE(NBTWOELECTRON3O)
            ALLOCATE(NBTWOELECTRON3O( NSTRIPV, NSTRIPC, NSTRIPV, NSTRIPC, NBLOCKS), STAT = ISTAT)
            !$ACC ENTER DATA CREATE(NBTWOELECTRON3O) __IF_ASYNC__
            !$ACC KERNELS PRESENT(NBTWOELECTRON3O) __IF_ASYNC__
            NBTWOELECTRON3O = 0._q
            !$ACC END KERNELS
          ENDIF
       ENDIF

       ! index of the first block calculated locally
       FIRST_LOCAL_BLOCK = FIRST_GLOBAL_BLOCK + (WDES%COMM%NODE_ME-1)*NBLOCKS

#ifdef bsedebug
       IF(IO%IU0>=0) WRITE(*,*) FIRST_GLOBAL_BLOCK, ":", NBLOCKS, 'blocks is computed on', NRANKS, 'ranks'
#endif

       !IF (IO%IU0>=0) CALL VTIME(VPUT0,CPUT0)
       ! loop over local blocks
       localblock: DO LOCAL_BLOCK = 1,NBLOCKS

          ! index of the block calculated locally in this iteration
          BLOCK_INDEX = FIRST_LOCAL_BLOCK + LOCAL_BLOCK - 1

          ! number of groups
          KPAR = W%WDES%COMM_KINTER%NCPU

          ! block on first rank in the k group
          BLOCK_INDEX_START = BLOCK_INDEX-(W%WDES%COMM_INTER%NODE_ME-1)*NBLOCKS

          ! block on last rank in the k group
          BLOCK_INDEX_END=MIN(BSE_INDEX%BLOCKS,BLOCK_INDEX_START + CEILING(REAL(W%WDES%COMM%NCPU*NBLOCKS,q)/REAL(KPAR,q)) - 1)

          ! loop over all blocks calculated in this iteration on all ranks
          DO BLOCK_INDEX_DIST=BLOCK_INDEX_START,BLOCK_INDEX_END,NBLOCKS

             ! row index of the block
             K1 = BSE_INDEX%BK1(BLOCK_INDEX_DIST)
             ISP = BSE_INDEX%BISP1(BLOCK_INDEX_DIST)
             B1=BD(ISP)
             NPOS1 = BSE_INDEX%BNPOS1(BLOCK_INDEX_DIST)
             NPOS3 = BSE_INDEX%BNPOS3(BLOCK_INDEX_DIST)
             K3 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

             ! column index of the block
             K2 = BSE_INDEX%BK2(BLOCK_INDEX_DIST)
             ISP2 = BSE_INDEX%BISP2(BLOCK_INDEX_DIST)
             B2=BD(ISP2)
             NPOS2 = BSE_INDEX%BNPOS2(BLOCK_INDEX_DIST)
             NPOS4 = BSE_INDEX%BNPOS4(BLOCK_INDEX_DIST)
             K4 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K2)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

             NODE = CEILING(REAL(BLOCK_INDEX_DIST-FIRST_GLOBAL_BLOCK+1,q) / REAL(NBLOCKS,q))

             CALL W1_GATHER_KSEL_NOFFT(WHF, NPOS1, NPOS1+NSTRIPV-1, ISP, W1, NODE,K1)
             CALL W1_GATHER_KSEL_NOFFT( WHF, NPOS2, NPOS2+NSTRIPV-1, ISP2, W2, NODE,K2)
             CALL W1_GATHER_KSEL_NOFFT(WHF,NPOS3, NPOS3+NSTRIPC-1, ISP, W3, NODE,K3)
             CALL W1_GATHER_KSEL_NOFFT( WHF, NPOS4, NPOS4+NSTRIPC-1, ISP2, W4, NODE,K4)
          ENDDO
          !$ACC WAIT IF(OFFLOAD_ON)

          ! coordinates of the block calculated locally
          IF(BLOCK_INDEX <= BSE_INDEX%BLOCKS .AND. W%WDES%COMM%NODE_ME <= NRANKS) THEN
             ! row index of the block
             K1 = BSE_INDEX%BK1(BLOCK_INDEX)
             ISP = BSE_INDEX%BISP1(BLOCK_INDEX)
             B1=BD(ISP)
             NPOS1 = BSE_INDEX%BNPOS1(BLOCK_INDEX)
             NPOS3 = BSE_INDEX%BNPOS3(BLOCK_INDEX)
             K3 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

             ! col index of the block
             K2 = BSE_INDEX%BK2(BLOCK_INDEX)
             ISP2 = BSE_INDEX%BISP2(BLOCK_INDEX)
             B2=BD(ISP2)
             NPOS2 = BSE_INDEX%BNPOS2(BLOCK_INDEX)
             NPOS4 = BSE_INDEX%BNPOS4(BLOCK_INDEX)
             K4 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K2)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)


             NSTRIP1=MIN(B1%VBMAXpCORE+1-NPOS1,NSTRIPV)
             NSTRIP2=MIN(B2%VBMAXpCORE+1-NPOS2,NSTRIPV)

             NSTRIP3=MIN(B1%CBMAX+1-NPOS3,NSTRIPC)
             NSTRIP4=MIN(B2%CBMAX+1-NPOS4,NSTRIPC)


GPACC        CALL ACC_SET_ASYNC_Q(1)
             CALL SETWDES(WHF%WDES,WDESK1,K1)
             CALL FFTWAV_WA(WA1)
             NFFTW=NFFTW+NSTRIP1

GPACC        CALL ACC_SET_ASYNC_Q(2)
             CALL SETWDES(WHF%WDES,WDESK2,K2)
             CALL FFTWAV_WA(WA2)
             NFFTW=NFFTW+NSTRIP2

GPACC        CALL ACC_SET_ASYNC_Q(3)
             CALL SETWDES(WHF%WDES,WDESK3,K3)
             CALL FFTWAV_WA(WA3)
             NFFTW=NFFTW+NSTRIP3

GPACC        CALL ACC_SET_ASYNC_Q(4)
             CALL SETWDES(WHF%WDES,WDESK4,K4)
             CALL FFTWAV_WA(WA4)
             NFFTW=NFFTW+NSTRIP4

!$ACC WAIT IF(OFFLOAD_ON)
GPACC     ACC_ASYNC_Q=ACC_ASYNC_ASYNC

             ! CDER_BETWEEN_STATES_Q is copmuted in ONEELECTRON_TRANS of difference ranks
             ! it is summed up in OSZI routeins, but requires that each pair K1 and ISP computed once
             ! all other ranks must have zeros
             IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
                DO I=1,3
                   QPOINT(I)=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(I,1)+ &
                             (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(I,2)+ &
                             (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(I,3)
                ENDDO

                IF (K1==K2 .AND. NPOS1==NPOS2 .AND. NPOS3==NPOS4 .AND. ISP==ISP2) THEN
                   ! set transition matrix elements
                   CALL ONEELECTRON_TRANS(WHF, P, LATT_CUR, ISP, WGW, &
                      W1, K1, NPOS1, NPOS1+NSTRIP1-1, &
                      W3, K3, NPOS3, NPOS3+NSTRIP3-1,B1, KPOINT_BSE, NFFTW )
                ENDIF
             ENDIF
             LIDLE=.FALSE.
          ELSE
             ! set LIDLE if no block is calculated on the rank
             LIDLE=.TRUE.
          ENDIF

!#ifndef bsedebug
!          CALL GWPROGRESS(IO%IU0, FIRST_GLOBAL_BLOCK, BSE_INDEX%BLOCKS, 1, 1 )
!#endif
          !==========================================================================
          !  set up BSE matrix as well as four point Coulomb and exchange integrals
          !  for the subset of bands and k-points considered
          !==========================================================================
          IF (.NOT. LIDLE) THEN

             TWOELECTRON3O=>NTWOELECTRON3O(:,:,:,:,LOCAL_BLOCK)
             IF (ANTIRES>=2) BTWOELECTRON3O=>NBTWOELECTRON3O(:,:,:,:,LOCAL_BLOCK)

             !$ACC ENTER DATA COPYIN(TWOELECTRON3O,BTWOELECTRON3O) __IF_ASYNC__

             !$ACC KERNELS PRESENT(TWOELECTRON3O,BTWOELECTRON3O) __IF_ASYNC__
             TWOELECTRON3O=0
             IF (ANTIRES>=2) BTWOELECTRON3O=0
             !$ACC END KERNELS


             IF (LHARTREE .AND. .NOT. LTRIPLET) THEN
                CALL TWOE4O_ACC_HARTREE(WHF, P, LATT_CUR, ISP, ISP2, WGW, &
                     LFXC .AND. .NOT. LGWLF, .FALSE. , AHARTREE,  &
                     W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, &
                     W3, K3, NPOS3, NSTRIP3, W4, K4, NPOS4, NSTRIP4, &
                     TWOELECTRON3O, B1, B2, NFLOAT4O, NFFT4O)
             ENDIF

             ! copy Coulomb integrals to resonant-antiresonant part
             IF (ANTIRES>=2) CALL WGCOPY(SIZE(TWOELECTRON3O),TWOELECTRON3O,1,BTWOELECTRON3O,1)


             IF (ISP2==ISP) THEN
                IF (LADDER .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF) ) THEN

                   ! calculate direct term from exchange (contribution to A matrix)
                   CALL TWOE4O_ACC_DIRECT(WHF, P, LATT_CUR, ISP, WGW, DELTA_COND,  &
                        W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2,&
                        W3, K3, NPOS3, NSTRIP3, W4, K4, NPOS4, NSTRIP4,&
                        TWOELECTRON3O, B1, B2, NFLOAT4O, NFFT4O)

                   ! calculate resonant-antiresonant term from exchange and add it to the A matrix
                   IF (ANTIRES==1) THEN
                      CALL TWOE4O_ACC_EX(WHF, P, LATT_CUR, ISP, WGW, 0.0_q, &
                           W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2,&
                           W3, K3, NPOS3, NSTRIP3, W4, K4, NPOS4, NSTRIP4,&
                           TWOELECTRON3O, B1, B2, NFLOAT4O, NFFT4O)
                   ENDIF
                   ! calculate resonant-antiresonant term from exchange and add it to the B matrix
                   IF (ANTIRES>=2) THEN
                      CALL TWOE4O_ACC_EX(WHF, P, LATT_CUR, ISP, WGW, 0.0_q, &
                           W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2,&
                           W3, K3, NPOS3, NSTRIP3, W4, K4, NPOS4, NSTRIP4,&
                           BTWOELECTRON3O, B1, B2, NFLOAT4O, NFFT4O)
                   ENDIF

                ENDIF ! condition LADDER=.TRUE. .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF)
             ENDIF ! condition ISP2==ISP

             !$ACC EXIT DATA DELETE(TWOELECTRON3O,BTWOELECTRON3O) __IF_ASYNC__
             TWOELECTRON3O=>NULL()
             BTWOELECTRON3O=>NULL()
          ENDIF ! LIDLE
       ENDDO localblock
!==========================================================================
! now store the calculated contributions in the final arrays
!==========================================================================

       IF (IBSE==2) THEN
          CALL TWOELECTRON4O_STORE_SCALA_TRIANGLE(WHF, ISP, ISP2, ISP_LOW, &
             K1, NPOS1, NSTRIPV, K2, NPOS2, NSTRIPV, &
             K3, NPOS3, NSTRIPC, K4, NPOS4, NSTRIPC, &
             NTWOELECTRON3O, BD, AMAT_SCALA, NCV, BSE_INDEX,1, ALPHAV(IALPHA), &
             SCALE_PERT, BSE_DESC, FIRST_LOCAL_BLOCK, NQPOINT, NRANKS, IBSE, LTRIANGLE)

          IF (ANTIRES>=2) THEN
             CALL TWOELECTRON4O_STORE_SCALA_TRIANGLE(WHF, ISP, ISP2, ISP_LOW, &
                K1, NPOS1, NSTRIPV, K2, NPOS2, NSTRIPV, &
                K3, NPOS3, NSTRIPC, K4, NPOS4, NSTRIPC, &
                NBTWOELECTRON3O, BD, BMAT_SCALA, NCV, BSE_INDEX,0, ALPHAV(IALPHA), &
                SCALE_PERT, BSE_DESC, FIRST_LOCAL_BLOCK, NQPOINT, NRANKS, IBSE, LTRIANGLE)
          ENDIF
       ENDIF ! IBSE==2

       ! IBSE==1, time-evolution code, IBSE==3, Lanczos code
       IF (IBSE==1 .OR. IBSE==3) THEN
               CALL TWOELECTRON4O_STORE_SCALA_TRIANGLE(WHF, ISP, ISP2, ISP_LOW, &
                    K1, NPOS1, NSTRIPV, K2, NPOS2, NSTRIPV, &
                    K3, NPOS3, NSTRIPC, K4, NPOS4, NSTRIPC, &
                    NTWOELECTRON3O, BD, AMAT_SCALA, NCV, BSE_INDEX,1, ALPHAV(IALPHA), &
                    SCALE_PERT, BSE_DESC, FIRST_LOCAL_BLOCK, NQPOINT, NRANKS, IBSE, LTRIANGLE)
         IF (ANTIRES==2) THEN
               CALL TWOELECTRON4O_STORE_SCALA_TRIANGLE(WHF, ISP, ISP2, ISP_LOW, &
                    K1, NPOS1, NSTRIPV, K2, NPOS2, NSTRIPV, &
                    K3, NPOS3, NSTRIPC, K4, NPOS4, NSTRIPC, &
                    NBTWOELECTRON3O, BD, AMAT_SCALA, NCV, BSE_INDEX, 0, ALPHAV(IALPHA), &
                    SCALE_PERT, BSE_DESC, FIRST_LOCAL_BLOCK, NQPOINT, NRANKS,IBSE, LTRIANGLE)
         ENDIF
       ENDIF
       IF (IO%IU0>=0) THEN
         CALL VTIME(VPUT,CPUT)
         WRITE(*,'(A,I10,A,I10,A,I10,A,I5,A,2F12.2,A)') "  [",  FIRST_GLOBAL_BLOCK, " - ", &
         MIN(FIRST_GLOBAL_BLOCK + NRANKS*NBLOCKS -1,BSE_INDEX%BLOCKS), "/", BSE_INDEX%BLOCKS , "] ", &
         NRANKS, " ", CPUT-CPUT0, VPUT-VPUT0, " "
       ENDIF
       ! increment block index
       FIRST_GLOBAL_BLOCK = FIRST_GLOBAL_BLOCK + NRANKS*NBLOCKS
    ENDDO globalblock

    IF(IBSE==1 .OR. IBSE==3) THEN
       !$ACC KERNELS PRESENT(AMAT_SCALA) __IF_ASYNC__
       AMAT_SCALA=GCONJG(AMAT_SCALA)
       !$ACC END KERNELS
    ENDIF

    IFAIL = 0
    CALL BLACS_GRIDINFO(BSE_DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
    IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
       IF(IBSE==1 .OR. IBSE==3) THEN
          CALL CHECK_BSEMATRIX(AMAT_SCALA,SIZE(AMAT_SCALA),IFAIL)
       ELSE
          CALL CHECK_BSEMATRIX(AMAT_SCALA,SIZE(AMAT_SCALA),IFAIL)
          IF(ANTIRES==2) CALL CHECK_BSEMATRIX(BMAT_SCALA,SIZE(BMAT_SCALA),IFAIL)
       ENDIF
    ENDIF

    IF(IFAIL > 0) CALL vtutor%bug(str(IFAIL) // " values in the BSE matrix are not set correctly", __FILE__, __LINE__)


!==========================================================================
!  end of loop body
!==========================================================================
#ifdef ACC_OFFLOAD
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    !$ACC EXIT DATA DELETE(NTWOELECTRON3O,NBTWOELECTRON3O) __IF_ASYNC__

    DO I=1,4
       CALL FFT_OFFLOAD_DESTROYPLAN(WDESK1%GRID,I,NSTRIPV)
       CALL FFT_OFFLOAD_DESTROYPLAN(WDESK1%GRID,I,NSTRIPC)

       DO ISP=ISP_LOW,ISP_HIGH
          CALL FFT_OFFLOAD_DESTROYPLAN_BSE(GRIDHF,BD(ISP)%VBMIN,BD(ISP)%VBMAXpCORE,NSTRIPV,I)
          CALL FFT_OFFLOAD_DESTROYPLAN_BSE(WGW%GRID,BD(ISP)%VBMIN,BD(ISP)%VBMAXpCORE,NSTRIPV,I)

          CALL FFT_OFFLOAD_DESTROYPLAN_BSE(GRIDHF,BD(ISP)%CBMIN,BD(ISP)%CBMAX,NSTRIPC,I)
          CALL FFT_OFFLOAD_DESTROYPLAN_BSE(WGW%GRID,BD(ISP)%CBMIN,BD(ISP)%CBMAX,NSTRIPC,I)
       ENDDO
    ENDDO

    CALL DELETE_TYPED_VAR(W4)
    CALL DELETE_TYPED_VAR(W3)
    CALL DELETE_TYPED_VAR(W2)
    CALL DELETE_TYPED_VAR(W1)

    !$ACC EXIT DATA DELETE(W1(:),W2(:),W3(:),W4(:)) __IF_ASYNC__

    CALL DELWAVA(WA1)
    CALL DELWAVA(WA2)
    CALL DELWAVA(WA3)
    CALL DELWAVA(WA4)

    !$ACC EXIT DATA DELETE(WA1,WA2,WA3,WA4) __IF_ASYNC__

    CALL DELETE_TYPED_VAR(WDESK4)
    CALL DELETE_TYPED_VAR(WDESK3)
    CALL DELETE_TYPED_VAR(WDESK2)
    CALL DELETE_TYPED_VAR(WDESK1)
#endif

    CALLMPI( M_sum_d(WGW%COMM, NFFTW, 1))
    CALLMPI( M_sum_d(WGW%COMM, NFLOAT4O, 1))
    CALLMPI( M_sum_d(WGW%COMM, NFFT4O  , 1))

10  FORMAT(" BLAS level 3 operations / number of FFT's:"/ &
           " number of FFTs for wave wavefunctions          ",F10.0," fft"/ &
           " number of operations in four-orbital integrals ",F14.2," Gflops, ",F13.0," fft")

    IF (IO%IU6>0) WRITE(IO%IU6,10) NFFTW,NFLOAT4O/1E9, NFFT4O

    IF (IALPHA==NALPHA) THEN
       !deallocate wavefunctions only if you are evaluating the last value for
       !the coupling constant integration
       ! ok finished get rid of WPOTH
       IF (ASSOCIATED(WPOTH)) CALL DESTROY_WPOT_HANDLE(WPOTH)
!==========================================================================
! deallocation
!==========================================================================

       CALL DEREGISTER_ALLOCATE(8._q*ndata8* SIZE(NTWOELECTRON3O, KIND=qi8) , "bse")
       DEALLOCATE(NTWOELECTRON3O)
       IF (ANTIRES>=2) THEN
          CALL DEREGISTER_ALLOCATE(8._q*ndata8* SIZE(NBTWOELECTRON3O, KIND=qi8) , "bse")
          DEALLOCATE(NBTWOELECTRON3O)
       ENDIF
    ENDIF
!==========================================================================
! BSE  redistribution of matrix elements
!==========================================================================
#if defined(CUBLASMP) || defined(CUSOLVERMP)
    CALL MP_GRID_FREE()
    CALL MP_COMM_FREE()
#endif
    CALL STOP_TIMING("G",IO%IU6,'BSESET')

    !CALL MPI_barrier(WHF%WDES%COMM, IFAIL) ! wait for other nodes
    IF (IO%IU0>=0) WRITE(IO%IU0,*)
    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE redistributing all elements'

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
!==========================================================================
! at this point the BSE matrix has been set up
!==========================================================================
! temporarily free GPU memory (see comment in declarations)

    IF(ALLOCATED(CDER_BETWEEN_STATES)) THEN
      ALLOCATE(CDER_BETWEEN_STATES_BACKUP(SIZE(CDER_BETWEEN_STATES, 1), &
         SIZE(CDER_BETWEEN_STATES, 2), &
         SIZE(CDER_BETWEEN_STATES, 3), &
         SIZE(CDER_BETWEEN_STATES, 4), &
         SIZE(CDER_BETWEEN_STATES, 5)))
      CDER_BETWEEN_STATES_BACKUP = CDER_BETWEEN_STATES
      DEALLOCATE(CDER_BETWEEN_STATES) ! implicitly frees GPU copy due to "acc declare create"
    ENDIF
#ifdef ACC_OFFLOAD
    !$ACC EXIT DATA DELETE(BSE_INDEX%INDEX,BSE_INDEX%N1,BSE_INDEX%N3,BSE_INDEX%NK,BSE_INDEX%ISP,BSE_INDEX%BK1, &
    !$ACC                  BSE_INDEX%BK2,BSE_INDEX%BISP1,BSE_INDEX%BISP2,BSE_INDEX%CK1,BSE_INDEX%CISP1, &
    !$ACC                  BSE_INDEX%CNPOS1,BSE_INDEX%CNPOS3,BSE_INDEX%BNPOS1,BSE_INDEX%BNPOS2,BSE_INDEX%BNPOS3, &
    !$ACC                  BSE_INDEX%BNPOS4,BSE_INDEX%BMAP,BSE_INDEX%MASK) __IF_ASYNC__
    !$ACC EXIT DATA DELETE(BSE_INDEX) __IF_ASYNC__
    CALL DELETE_TYPED_VAR(KPOINTS_FULL_ORIG)
    CALL DELETE_TYPED_VAR(KPOINTS_FULL)
    CALL DELETE_TYPED_VAR(WHF)

    CALL ACC_CLEAR_FREELISTS()
#endif
    IF (IBSE==2) THEN

    ! deallocate orbitals
    ! this is usually allowed but obviously an issue if we desire correlation
    ! energies, since the routine is called repeatedly
    ! deallocation frees precious memory for the scaLAPACK calls
    ! if we don't do exciton WF or BSEFATBAND
    IF (NBSEEIG < 1)  THEN
       CALL DEALLOCW_CW(WHF)
    ENDIF

    ! store A-B in A, and A+B in B
    IF (ANTIRES>=2) THEN
!==========================================================================
! for the plasmon formulate of the correlation energy we need
! the sum of the diagonal components of the A matrix
!==========================================================================
       CSUM_A=(0.0_q,0.0_q)
       IF (.NOT. LADDER) THEN
          !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          CALL TWOELECTRON4O_GET_DIAG_SCALA( AMAT_SCALA, NCV, CSUM_A, WHF%WDES%COMM, BSE_DESC)
          IF (REAL(CSUM_A) /= REAL(CSUM_A) ) WRITE (*,*) 'WARNING: diag(A) contains NaN elements'
       ENDIF
!==========================================================================
!  full beyond TDA case   (store A+B in A, and A-B in B)
!==========================================================================
       !$ACC PARALLEL LOOP PRESENT(AMAT_SCALA,BMAT_SCALA) PRIVATE(TMP) __IF_ASYNC__
       DO K8=1_qi8,SIZE(AMAT_SCALA,KIND=qi8)
          TMP=AMAT_SCALA(K8)
          AMAT_SCALA(K8)=TMP+BMAT_SCALA(K8)
          BMAT_SCALA(K8)=TMP-BMAT_SCALA(K8)
       ENDDO

       ! diagonalize A-B
       CALL DIAG_BSE_MATRIX( NCV, BMAT_SCALA, RB, LDISTRIBUTED_MATRIX, WDES%COMM, IO%IU0, BSE_DESC)

       CALL STOP_TIMING("G",IO%IU6,'BSEDIAG1')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
       ! calculate  (A-B)^1/2 (A+B) (A-B)^1/2
       CALL TDA_BSE_SCA( NCV, AMAT_SCALA, BMAT_SCALA, RB , IO%IU0, BSE_DESC)
       CALL STOP_TIMING("G",IO%IU6,'SQUARE')

       ! diagonalize (A-B)^1/2 (A+B) (A-B)^1/2
       CALL DIAG_BSE_MATRIX( NCV, AMAT_SCALA, R, LDISTRIBUTED_MATRIX, WDES%COMM, IO%IU0, BSE_DESC)

       CALL STOP_TIMING("G",IO%IU6,'BSEDIAG2')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
       ! eigenvalues of original problem
       CALL SQRT_EIGENVALUE_INPLACE(R, IO%IU0)

       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)

       !$ACC UPDATE SELF(R) WAIT IF(OFFLOAD_ON)
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       ! compute the correlation energy using the plasmon formula only when
       ! LAMBDA=1.0 is set in the INCAR file (and there are no exchange terms)git
       IF (ALPHAV(IALPHA)==1.0_q .AND. NALPHA==1 .AND. (.NOT. LADDER)) THEN
           IF (IO%IU0>=0) WRITE(IO%IU0,'(A,1F20.10,A,I6)') ' plasmon correlation energy ', &
                              0.5*(SUM(R)-REAL(CSUM_A,q))*KPOINTS_ORIG%WTKPT(NQPOINT)/SCALE_PERT**2, ' @ NQ ', NQPOINT
           IF (IO%IU6>=0) WRITE(IO%IU6,'(A,1F20.10,A,I6)') ' plasmon correlation energy ', &
                              0.5*(SUM(R)-REAL(CSUM_A,q))*KPOINTS_ORIG%WTKPT(NQPOINT)/SCALE_PERT**2, ' @ NQ ', NQPOINT
       ENDIF

       ! left multiply by (A-B)^1/2 to obtain eigenvectors
       ! to calculate optical properties we need only eigenvectors of squared problem
       CALL TDA_BSE_EIGENVECTORS_SCA( NCV, BMAT_SCALA, RB, AMAT_SCALA, R, IO%IU0, BSE_DESC )

       CALL STOP_TIMING("G",IO%IU6,'EIGVEC')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
    ELSE
!==========================================================================
! simple TDA standard algorithms
!==========================================================================
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CALL DIAG_BSE_MATRIX( NCV, AMAT_SCALA, R, LDISTRIBUTED_MATRIX, WDES%COMM, IO%IU0, BSE_DESC )
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       !$ACC UPDATE SELF(R) WAIT IF(OFFLOAD_ON)
       !$ACC UPDATE SELF(AMAT_SCALA) WAIT IF(OFFLOAD_ON)
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CALL STOP_TIMING("G",IO%IU6,'BSEDIAG')
    ENDIF ! ANTIRES>=2
    !$ACC UPDATE SELF(R) WAIT IF(OFFLOAD_ON)
    ENDIF ! IBSE == 2
                                
    IF(ALLOCATED(CDER_BETWEEN_STATES_BACKUP)) THEN
      ! restore temporarily freed GPU memory, note: CDER_BETWEEN_STATES is "acc declare create"
      ALLOCATE(CDER_BETWEEN_STATES(SIZE(CDER_BETWEEN_STATES_BACKUP, 1), &
                                   SIZE(CDER_BETWEEN_STATES_BACKUP, 2), &
                                   SIZE(CDER_BETWEEN_STATES_BACKUP, 3), &
                                   SIZE(CDER_BETWEEN_STATES_BACKUP, 4), &
                                   SIZE(CDER_BETWEEN_STATES_BACKUP, 5)))
      CDER_BETWEEN_STATES = CDER_BETWEEN_STATES_BACKUP
      DEALLOCATE(CDER_BETWEEN_STATES_BACKUP)
    ENDIF
#ifdef ACC_OFFLOAD
    CALL COPYIN_TYPED_VAR(WHF)
    CALL COPYIN_TYPED_VAR(KPOINTS_FULL)
    CALL COPYIN_TYPED_VAR(KPOINTS_FULL_ORIG)
    !$ACC ENTER DATA COPYIN(BSE_INDEX) __IF_ASYNC__
    !$ACC ENTER DATA COPYIN(BSE_INDEX%BK1,BSE_INDEX%BK2,BSE_INDEX%BISP1, BSE_INDEX%BISP2,BSE_INDEX%BNPOS1,BSE_INDEX%BNPOS2,BSE_INDEX%INDEX) __IF_ASYNC__
    !$ACC ENTER DATA COPYIN(BSE_INDEX%N1,BSE_INDEX%N3,BSE_INDEX%NK,BSE_INDEX%ISP) __IF_ASYNC__
#endif

    IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
!================================================================================
  ENDDO ccintegration ! over values of the coupling constant (IALPHA)
!================================================================================

  IF (NBSEEIG > NCV) THEN
     IF (IO%IU0>=0) THEN
       WRITE(*,*) "Warning: NBSEEIG cannot be larger than the number of BSE eigenstates."
       WRITE(*,*) "NBSEEIG is set to the total number of BSE eigenstates."
     ENDIF
     NBSEEIG=NCV
  ENDIF

  IF (ANY(BSEHOLE /= REAL(MAGIC,KIND=q)) .OR. ANY(BSEELECTRON /= REAL(MAGIC,KIND=q))) THEN
     IF(ANTIRES /= 0) THEN ! only TDA is currently implemented
       CALL vtutor%alert("The calculation of the exciton WF is only supported in the Tamm-Dancoff approximation. Set ANTIRES=0")
     ELSE IF (IBSE/=2) THEN ! must be exact diagonalization
       CALL vtutor%alert("The calculation of the exciton WF is only supported for the exact diagonalization algorithm")
     ELSE IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
       CALL vtutor%alert("The calculation of the exciton WF is only supported for the q=0 case.")
     ELSE IF ((ANY(BD(:)%VBMAXpCORE > BD(:)%VBMAX)) .AND. (NBANDSO /= 0)) THEN ! only core states in XAS BSE
       CALL vtutor%alert("The calculation of the exciton WF with core excitations can only be perfromed with NBANDSO=0.")
     ELSE IF (NBSEEIG ==0) THEN
       CALL vtutor%alert("The calculation of the exciton WF requires that the number of states is specified. Set NBSEEIG.")
     ELSE ! only if all above false we do EXCITON_WF
       CALL EXCITON_WF(T_INFO,INFO,P,WHF,AMAT_SCALA,LATT_CUR,LMDIM,IO,NCV,BSE_DESC,BSE_INDEX,BD,ISP_LOW,ISP_HIGH,OMEGAMAX)
       CALL STOP_TIMING("G",IO%IU6,'EWF')
     ENDIF
  ENDIF

  IF(ALLOCATED(CDER_BETWEEN_STATES_Q)) THEN
     CALLMPI( M_sum_z8(WHF%WDES%COMM, CDER_BETWEEN_STATES_Q , SIZE(CDER_BETWEEN_STATES_Q,KIND=qi8)))
  ENDIF
  IF(ALLOCATED(CDER_CORECON_Q)) THEN
     CALLMPI( M_sum_z8(WHF%WDES%COMM, CDER_CORECON_Q , SIZE(CDER_CORECON_Q,KIND=qi8)))
  ENDIF

  IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
!==========================================================================
! use time-evolution or lanczos algorithm
!==========================================================================
    IF (IBSE==1 .OR. IBSE==3) THEN

      !$ACC ENTER DATA CREATE(OPTMAT) __IF_ASYNC__
      CALL INIT_BSE_OSZI_STRENGTH(WHF, LATT_CUR, ISP_LOW, ISP_HIGH, BD, BSE_INDEX, NSTRIPV, NSTRIPC, NQPOINT, OPTMAT)
      !$ACC WAIT IF(OFFLOAD_ON)
      IF (IBSE==1) THEN
      CALL CALCULATE_BSE_TIME_EVOLUTION(WHF,SHIFT,LATT_CUR%OMEGA, &
           MAX(WDES%ISPIN,WDES%NRSPINORS), WDES%NKPTS, NCV, WDES%WTKPT, &
           OMEGAMAX, LADDER .OR. (LHARTREE .AND. .NOT. LTRIPLET), &
           ANTIRES==2, NEDOS, IO, QPOINT, AMAT_SCALA, OPTMAT)
        CALL STOP_TIMING("G",IO%IU6,'BSE_TE')
      ELSEIF (IBSE==3) THEN
! Lanczos is not ported to GPU yet ... switch to CPU execution if necessary
!$ACC UPDATE SELF(AMAT_SCALA,OPTMAT) __IF_ASYNC__
GPACC   CALL OFFLOAD_ON_PUSH(.FALSE.)
        CALL CALCULATE_BSE_LANCZOS(WHF%WDES%COMM,SHIFT,LATT_CUR%OMEGA, &
             MAX(WDES%ISPIN,WDES%NRSPINORS), WDES%NKPTS, NCV, WDES%WTKPT, &
             OMEGAMAX, ANTIRES==2, NEDOS, IO, NBSEEIG, QPOINT, AMAT_SCALA, OPTMAT)
! Lanczos ran CPU side ... switch back to GPU execution is possible
GPACC   CALL OFFLOAD_ON_POP
!$ACC UPDATE DEVICE(AMAT_SCALA,OPTMAT) __IF_ASYNC__
        CALL STOP_TIMING("G",IO%IU6,'BSE_LA')
      ENDIF
      !$ACC EXIT DATA DELETE(OPTMAT) __IF_ASYNC__
    ENDIF


!==========================================================================
! calculate oscillator strength
!==========================================================================
    IF (IBSE==2) THEN
       ! the optical properties calculation is no longer compatible with the rest
       ! of the code when QP energies are computed because SET_BAND_PARAMETERS
       ! includes a few unoccupied states in the VB.
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE calculating oscillator strength'
       CALL CALCULATE_BSE_OSZI_SCALA(WHF, NQPOINT, KPOINT_BSE, LATT_CUR, ISP_LOW, ISP_HIGH, &
       SHIFT, BD, BSE_INDEX, R, AMAT_SCALA, NCV, ANTIRES, NEDOS, OMEGAMAX, IO%IU6, BSE_DESC)
       CALL STOP_TIMING("G",IO%IU6,'BSEOSZI')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
    ENDIF

!==========================================================================
! BSEFATBANDS
!==========================================================================
    IF (IBSE==2 .AND. NBSEEIG > 0) THEN
    IF (IO%IU0>=0) WRITE(IO%IU0,"(A17,I6,A13)") ' BSE writing out ',NBSEEIG,' eigenvectors'

!$ACC UPDATE SELF(AMAT_SCALA) IF(OFFLOAD_ON)
    CALL WRITE_BSE_SCALA(WHF, ISP_LOW, ISP_HIGH, &
            BD, BSE_INDEX, R, AMAT_SCALA, NCV, NBSEEIG, IO%IU6, BSE_DESC)
#ifdef VASP_HDF5
    CALL VH5_WRITE_BSEINDEX(IH5OUTFILEID, BSE_INDEX%INDEX)
    CALL VH5_WRITE_BSEBANDS(IH5OUTFILEID, BD, EFERMI)
#endif
    ENDIF
!==========================================================================
! deallocation
!==========================================================================

#ifdef ACC_OFFLOAD

    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    !$ACC EXIT DATA DELETE(RB)         ASYNC(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.ALLOCATED(RB))
    !$ACC EXIT DATA DELETE(BMAT_SCALA) ASYNC(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.ALLOCATED(BMAT_SCALA))
    !$ACC EXIT DATA DELETE(AMAT_SCALA) ASYNC(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.ALLOCATED(AMAT_SCALA))
    !$ACC EXIT DATA DELETE(R) __IF_ASYNC__
    !$ACC EXIT DATA DELETE(BD,BSE_DESC,VEC_DESC,DVEC_DESC) __IF_ASYNC__
    !$ACC EXIT DATA DELETE(BSE_INDEX%INDEX,BSE_INDEX%N1,BSE_INDEX%N3,BSE_INDEX%NK,BSE_INDEX%ISP,BSE_INDEX%BK1, &
    !$ACC                   BSE_INDEX%BK2,BSE_INDEX%BISP1,BSE_INDEX%BISP2,BSE_INDEX%CK1,BSE_INDEX%CISP1, &
    !$ACC                   BSE_INDEX%CNPOS1,BSE_INDEX%CNPOS3,BSE_INDEX%BNPOS1,BSE_INDEX%BNPOS2,BSE_INDEX%BNPOS3, &
    !$ACC                   BSE_INDEX%BNPOS4,BSE_INDEX%BMAP,BSE_INDEX%MASK) __IF_ASYNC__
    !$ACC EXIT DATA DELETE(TRANS_MATRIX_FOCK,LATT_CUR,BSE_INDEX) __IF_ASYNC__

    CALL DELETE_TYPED_VAR(GRIDHF)
    CALL DELETE_TYPED_VAR(FAST_AUG_FOCK)
    CALL DELETE_TYPED_VAR(AUG_DES)
    CALL DELETE_TYPED_VAR(KPOINTS_FULL_ORIG)
    CALL DELETE_TYPED_VAR(KPOINTS_FULL)

    CALL DELETE_TYPED_VAR(WHF)

    IF (ALLOCATED(CDER_BETWEEN_STATES)) THEN
       DEALLOCATE(CDER_BETWEEN_STATES)
    ENDIF
    !CALL acc_present_dump()
#endif

    DEALLOCATE(W1,W2,W3,W4)
    IF(ALLOCATED(CDER_CORECON_Q)) DEALLOCATE(CDER_CORECON_Q)
    CALL DEALLOCATE_CDER_BETWEEN_STATES_Q
    CALL DEALLOCATE_BSE_MATRIX_INDEX( BSE_INDEX )
    NSIZE=SIZE(AMAT_SCALA,KIND=qi8)
    DEALLOCATE(AMAT_SCALA, R)
    IF(IBSE==1 .OR. IBSE==3) DEALLOCATE(OPTMAT)
    IF (ANTIRES>=2) THEN
       NSIZE=NSIZE+SIZE(BMAT_SCALA,KIND=qi8)
       DEALLOCATE(BMAT_SCALA, RB)
    ENDIF
#ifdef single_prec_bse
    CALL DEREGISTER_ALLOCATE(4._q*ndata8* NSIZE , "bse")
#endif
#ifdef double_prec_bse
    CALL DEREGISTER_ALLOCATE(8._q*ndata8* NSIZE , "bse")
#endif


    PROFILING_STOP('calculate_bse')

    RETURN

    CONTAINS


    ! ***********************************************************
    !> checks for allocation error and stops code with message
    ! ***********************************************************
    SUBROUTINE CHECK_ALLOCATION_ERROR( STRING, ISIZE, ISTAT )
       CHARACTER( LEN=* )  :: STRING
       INTEGER(qi8), INTENT(IN) :: ISIZE
       INTEGER, INTENT(IN) :: ISTAT

       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_BSE ("//STRING//") is not able to allocate "//&
#ifdef single_prec_bse
            str( 4._q*ISIZE/1024_qi8) //" kB of data on MPI rank 0." )
#endif
#ifdef double_prec_bse
            str( 8._q*ISIZE/1024_qi8) //" kB of data on MPI rank 0." )
#endif
       ENDIF

    END SUBROUTINE CHECK_ALLOCATION_ERROR

    ! ***********************************************************
    !> checks for allocation error and stops code with message
    !> same as #CHECK_ALLOCATION_ERROR, but prints only double precision sizes
    ! ***********************************************************
    SUBROUTINE CHECK_ALLOCATION_ERROR_Q( STRING, ISIZE, ISTAT )
       CHARACTER( LEN=* )  :: STRING
       INTEGER(qi8), INTENT(IN) :: ISIZE
       INTEGER, INTENT(IN) :: ISTAT

       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_BSE ("//STRING//") is not able to allocate "//&
            str( 8._q*ISIZE/1024_qi8) //" kB of data on MPI rank 0." )
       ENDIF

    END SUBROUTINE CHECK_ALLOCATION_ERROR_Q
#endif

  END SUBROUTINE CALCULATE_BSE_ACC



#ifdef scaLAPACK


!****************** SUBROUTINE  TWOELECTRON4O_STORE *******************
!
!> Stores the two electron 4 orbital integrals
!> in the BSE / TDFT Hamilton matrix AMAT.
!> This version of the routine is not scaLAPACK "aware" and
!> requires that the full BSE matrix is stored on every node.
!
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_STORE(WHF, ISP, ISP2, ISP_LOW, &
       K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
       TWOELECTRON3O, B1, B2, AMAT, BSE_INDEX, IDIAG, ALPHAPT, SCALE_PERT)

    USE constant
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP, ISP2, ISP_LOW
    INTEGER K1, K2, K3, K4, NPOS1, NSTRIP1, NPOS2, NSTRIP2
    TYPE(banddesc) :: B1,B2         !< current descriptor for the spin loops
    GDEF ::  TWOELECTRON3O(:,:,:,:)
    REAL(q) :: ALPHAPT              !< coupling constant value for AC integration
    REAL(q):: SCALE_PERT
    REAL(q):: WEIGHT

! actual dim TWOELECTRON3O( NSTRIPV, NGLB, NSTRIPV, NGLB4))
#ifdef single_prec_bse
    GDEFS :: AMAT(:,:)
#endif

#ifdef double_prec_bse
    GDEF :: AMAT(:,:)
#endif
  ! local variables
    INTEGER N1, N2, N3, N4, NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER NB3, NB4
    INTEGER NCV13, NCV24
    TYPE(bse_matrix_index) :: BSE_INDEX
    INTEGER :: IDIAG     ! add eigenvalues


    WEIGHT = WHF%WDES%WTKPT(1)
    PROFILING_START('twoelectron4o_store')
    !$ACC ENTER DATA COPYIN(BSE_INDEX,B1,B2) __IF_ASYNC__
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(AMAT,BSE_INDEX,TWOELECTRON3O,B1,B2,WHF) &
    !$ACC & PRIVATE(NB1_INTO_TOT,NB2_INTO_TOT,NB4_INTO_TOT,NB4,NCV24, &
    !$ACC & NB3_INTO_TOT,NB3,NCV13) __IF_ASYNC__
    DO N1=1,NSTRIP1
       DO N2=1,NSTRIP2
          NB1_INTO_TOT=NPOS1-1+N1
          NB2_INTO_TOT=NPOS2-1+N2
          DO N4=1,B2%CBMAX4-B2%CBMIN4+1
             NB4_INTO_TOT=B2%CBMIN4-1+N4
             NB4         =N4
#if defined (checkindex)  && ! defined (ACC_OFFLOAD)
             IF (N1 > SIZE(TWOELECTRON3O,1)) THEN
                CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 1 " // str(N1), __FILE__, __LINE__)
             ENDIF

             IF (N2 > SIZE(TWOELECTRON3O,3)) THEN
                CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 3 " // str(N2), __FILE__, __LINE__)
             ENDIF

             IF (NB4 > SIZE(TWOELECTRON3O,4)) THEN
                CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 4 " // str(NB4), __FILE__, __LINE__)
             ENDIF
#endif
           ! NCV24=NB2_INTO_TOT-B1%VBMIN+1 +(B1%VBMAX-B1%VBMIN+1)*(NB4_INTO_TOT-B1%CBMIN+(B1%CBMAX-B1%CBMIN+1)*(K2-1))
             NCV24=BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)
             IF (NCV24==0) CYCLE ! cycle if this pair is not included
#if defined (checkindex)  && ! defined (ACC_OFFLOAD)
             IF ( NCV24 > SIZE(AMAT,2) ) THEN
                CALL vtutor%bug(str(NB2_INTO_TOT-B2%VBMIN+1) // " " // str(B2%VBMAX-B2%VBMIN+1) // &
                   " " // str(NB4_INTO_TOT-B2%CBMIN+1) // " " // str(B2%CBMAX-B2%CBMIN+1) // " " // &
                   str(K2-1) // "\n internal error in TWOELECTRON4O_STORE: NCV exceeds AMAT &
                   &size-boundaries " // str(NCV13) // " " // str(NCV24) // " " // str(SIZE(AMAT,1)) &
                   // " " // str(SIZE(AMAT,2)), __FILE__, __LINE__)
             ENDIF
#endif
             DO N3=1,B1%CBMAX-B1%CBMIN+1
                NB3_INTO_TOT=B1%CBMIN-1+N3
                NB3         =N3
#if defined (checkindex)  && ! defined (ACC_OFFLOAD)
                IF ( NB3> SIZE(TWOELECTRON3O,2)) THEN
                   CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 2 " // &
                      str(NB3_INTO_TOT) // " " // str(NB3), __FILE__, __LINE__)
                ENDIF
#endif
                ! build index NCV13 and NCV24
                NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)
                IF (NCV13==0) CYCLE ! cycle if this pair is not included
                ! check whether NCV13 is smaller or equal than SIZE(AMAT,1)
#if defined (checkindex)  && ! defined (ACC_OFFLOAD)
                IF ( NCV13 > SIZE(AMAT,1) ) THEN
                   CALL vtutor%bug(str(NB2_INTO_TOT-B2%VBMIN+1) // " " // str(B2%VBMAX-B2%VBMIN+1) &
                      // " " // str(NB4_INTO_TOT-B2%CBMIN+1) // " " // str(B2%CBMAX-B2%CBMIN+1) // "&
                      & " // str(K2-1) // "\n internal error in TWOELECTRON4O_STORE: NCV exceeds AMAT &
                      &size-boundaries " // str(NCV13) // " " // str(NCV24) // " " // str(SIZE(AMAT,&
                      1)) // " " // str(SIZE(AMAT,2)), __FILE__, __LINE__)
                ENDIF
                IF (AMAT(NCV13,NCV24) /= 0) THEN
                   CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: AMAT was already written &
                      &" // str(NCV13) // " " // str(NCV24) // "\n attempt of overwriting points &
                      &towards inconsistent indexing \n correct the formulas for NCV13 and NCV24 in bse.F", &
                      __FILE__, __LINE__)
                ENDIF
#endif
                AMAT(NCV13,NCV24)=TWOELECTRON3O( N1, NB3 , N2, NB4)*WEIGHT*SCALE_PERT*ALPHAPT

                ! diagonal components
                ! subtract eigenvalue difference for the A matrix (IDIAG is 0.0 when BMAT is passed on)
                IF (K1==K2 .AND. NB1_INTO_TOT==NB2_INTO_TOT .AND. NB3_INTO_TOT==NB4_INTO_TOT .AND.  ISP== ISP2) THEN
#ifndef ACC_OFFLOAD
                   IF (NCV13/= NCV24) THEN
                      CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: not a diagonal element " &
                         // str(NCV13) // " " // str(NCV24) // "\n " // str(N1) // " " // str(N2) &
                         // " " // str(N3) // " " // str(N4) // " " // str(K1) // " " // str(K2), __FILE__, __LINE__)
                   ENDIF
#endif
                   IF(IDIAG>0) THEN
                      !! catch exception if the aufbau principle is violated (i.e.: E_k1,n1 > E_k3,n3)
                      IF((REAL(WHF%CELTOT(NB1_INTO_TOT,K1, ISP))>REAL(WHF%CELTOT(NB3_INTO_TOT,K3, ISP)))) THEN
                         WRITE(*,*) 'internal WARNING: Aufbau principle is violated, might need reexamination'
                      ELSE
                          AMAT(NCV13, NCV24) = AMAT(NCV13, NCV24)-IDIAG*(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$ACC EXIT DATA DELETE(BSE_INDEX,B1,B2) __IF_ASYNC__
    PROFILING_STOP('twoelectron4o_store')

  END SUBROUTINE TWOELECTRON4O_STORE


!****************** SUBROUTINE  TWOELECTRON4O_STORE_SCALA  ************
!
!> Stores the two electron 4 orbital integrals
!> in the BSE / TDFT Hamilton matrix AMAT
!> this version is scaLAPACK "aware" and matter of fact only
!> available if scalapack is used
!
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_STORE_SCALA(WHF, ISP_, ISP2_, ISP_LOW, &
       K1_, NPOS1_, NSTRIP1_, K2_, NPOS2_, NSTRIP2_, &
       K3_, NPOS3_, NSTRIP3_, K4_, NPOS4_, NSTRIP4_, &
       TWOELECTRON3O, BD, AMAT, NCV, BSE_INDEX, IDIAG, ALPHAPT, SCALE_PERT, DESC)

GPACC USE moffload

    USE constant
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP, ISP2, ISP_LOW
    INTEGER ISP_, ISP2_, ISP_LOW_
    INTEGER :: K1_, NPOS1_, NSTRIP1_, K2_, NPOS2_, NSTRIP2_
    INTEGER :: K3_, NPOS3_, NSTRIP3_, K4_, NPOS4_, NSTRIP4_
    TYPE(banddesc) :: BD(:)            !< current descriptors for the spin loop
    GDEF ::  TWOELECTRON3O(:,:,:,:)
! actual dim TWOELECTRON3O( NSTRIPV, NGLB, NSTRIPV, NGLB4))
    TYPE(bse_matrix_index) :: BSE_INDEX
    INTEGER :: IDIAG                !< add eigenvalues
    REAL(q) :: ALPHAPT              !< value of the coupling constant for AC integration
    REAL(q) :: SCALE_PERT
    INTEGER :: DESC(DLEN_)

#ifdef single_prec_bse
    GDEFS :: AMAT(:)
#endif
#ifdef double_prec_bse
    GDEF :: AMAT(:)
#endif
    INTEGER :: NCV
  ! local variables
    GDEF, ALLOCATABLE ::  TWOELECTRON3O_LOCAL(:,:,:,:)
    INTEGER N1, N2, NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER NB3, NB4, N3,N4
    INTEGER NCV13, NCV24
    INTEGER K1,K3, K2, K4, CBMIN4, CBMAX4, NPOS1, NPOS2, NSTRIP1, NSTRIP2
    INTEGER NPOS3, NPOS4, NSTRIP3, NSTRIP4
    INTEGER N
    INTEGER ITWOELECTRON3O(4)
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2
    INTEGER IFAIL

    PROFILING_START('twoelectron4o_store_scala')

! allocate work array into which data is subsequently merged from all nodes
    ITWOELECTRON3O=SHAPE(TWOELECTRON3O)
    CALL M_max_i(WHF%WDES%COMM, ITWOELECTRON3O, 4 )

    ALLOCATE(TWOELECTRON3O_LOCAL(ITWOELECTRON3O(1),ITWOELECTRON3O(2), &
                                 ITWOELECTRON3O(3),ITWOELECTRON3O(4)))

!==========================================================================
! loop over all nodes and broadcast their TWOELECTRO3O array
! to other nodes also broadcast the storage position (B2%CBMIN4, B2%CBMAX4)
! and the local k-point indices
!==========================================================================
cp: DO N=1,WHF%WDES%COMM%NCPU

! broadcast K2, K4 and all the other integers to all nodes
    IF (N==WHF%WDES%COMM%NODE_ME) THEN
       ISP    =ISP_
       ISP2   =ISP2_
       K1=K1_
       K2=K2_
       K3=K3_
       K4=K4_
       NPOS1=NPOS1_
       NPOS2=NPOS2_
       NPOS3=NPOS3_
       NPOS4=NPOS4_
       NSTRIP1=NSTRIP1_
       NSTRIP2=NSTRIP2_
       NSTRIP3=NSTRIP3_
       NSTRIP4=NSTRIP4_
    ENDIF
    CALL M_bcast_i_from(WHF%WDES%COMM, ISP, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, ISP2, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, K1, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, K2, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, K3, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, K4, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NPOS1, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NPOS2, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NPOS3, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NPOS4, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NSTRIP1, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NSTRIP2, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NSTRIP3, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM, NSTRIP4, 1, n)

! if K2 is properly set distribute the matrix elements
    IF (K1>=1 .AND. K2>=1 .AND. K3>=1 .AND. K4>=1) THEN

    TWOELECTRON3O_LOCAL=0
    !$ACC ENTER DATA COPYIN(TWOELECTRON3O_LOCAL) __IF_ASYNC__
    IF (N==WHF%WDES%COMM%NODE_ME) THEN
       !$ACC KERNELS PRESENT(TWOELECTRON3O_LOCAL) __IF_ASYNC__
       TWOELECTRON3O_LOCAL(1:SIZE(TWOELECTRON3O,1),1:SIZE(TWOELECTRON3O,2),1:SIZE(TWOELECTRON3O,3),1:SIZE(TWOELECTRON3O,4))=TWOELECTRON3O
       !$ACC END KERNELS
    ENDIF
#ifdef gammareal
    CALL M_bcast_d_from(WHF%WDES%COMM, TWOELECTRON3O_LOCAL, SIZE(TWOELECTRON3O_LOCAL), N)
#else
    CALL M_bcast_z_from(WHF%WDES%COMM, TWOELECTRON3O_LOCAL, SIZE(TWOELECTRON3O_LOCAL), N)
#endif
!==========================================================================
! now distribute the matrix elements among nodes
!==========================================================================
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    ! loop over all columns of the global matrix
    JCOL=0
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(AMAT,TWOELECTRON3O_LOCAL,BSE_INDEX,WHF) &
    !$ACC & PRIVATE(NB1_INTO_TOT,NB3_INTO_TOT,NB2_INTO_TOT,NB4_INTO_TOT,NCV13,NCV24) &
    !$ACC & PRIVATE(N1,N3,N2,N4) __IF_ASYNC__
    DO J1=1,NQ,DESC(NB_)
       DO J2=1,DESC(NB_)
          JCOL=J1+J2-1
          IF(JCOL > NQ ) CYCLE

        ! global column index NCV24
          NCV24=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2
        ! spin test
          IF (ISP2 /= BSE_INDEX%ISP(NCV24)) CYCLE
        ! k-point test, CYCLE if not identical
        ! IF (K2 /= 1+(NCV24-1)/((B1%VBMAX-B1%VBMIN+1)*(B1%CBMAX-B1%CBMIN+1))) CYCLE
          IF (K2 /= BSE_INDEX%NK(NCV24)) CYCLE

        ! test on NB2
        ! NB2_INTO_TOT=MOD(NCV24-1,B1%VBMAX-B1%VBMIN+1)+B1%VBMIN
          NB2_INTO_TOT=BSE_INDEX%N1(NCV24)
          N2=NB2_INTO_TOT-NPOS2+1
          IF (N2<1 .OR. N2 > NSTRIP2) CYCLE

        ! test on NB4_INTO_TOT
        ! NB4_INTO_TOT=MOD((NCV24-1)/(B1%VBMAX-B1%VBMIN+1),(B1%CBMAX-B1%CBMIN+1))+B1%CBMIN
          NB4_INTO_TOT=BSE_INDEX%N3(NCV24)
          !IF (NB4_INTO_TOT<BD(ISP2)%CBMAX .OR. NB4_INTO_TOT > BD(ISP2)%CBMAX) CYCLE

          NB4=NB4_INTO_TOT-BD(ISP2)%CBMIN+1
          N4=NB4_INTO_TOT-NPOS4+1

#ifdef check_index
          ! internal consistency test (could be removed)
        ! IF (NCV24/=NB2_INTO_TOT-B2%VBMIN+1 +(B2%VBMAX-B2%VBMIN+1)*(NB4_INTO_TOT-B2%CBMIN+(B2%CBMAX-B2%CBMIN+1)*(K2-1))) THEN
          IF (NCV24/=BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)) THEN
             CALL vtutor%bug("internal error 1 in TWOELECTRON4O_STORE_SCALA: " // str(NCV24) // " &
                &" // str(NB2_INTO_TOT-B2%VBMIN+1) // " " // str(NB4_INTO_TOT-B2%CBMIN) // " " // &
                str(K2) // " " // str(ISP2) // " " // str(BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1,&
                NB4_INTO_TOT-B2%CBMIN+1,K2,ISP2-ISP_LOW+1)) // " " // str(BSE_INDEX%ISP(NCV24)), __FILE__, __LINE__)
          ENDIF
#endif
          IROW=0
          !$ACC LOOP COLLAPSE(2)
          DO I1=1,NP,DESC(MB_)
             DO I2=1,DESC(MB_)
                IROW=I1+I2-1
                IF (IROW > NP) CYCLE

              ! global row index NCV13
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
              ! spin test
                IF (ISP /= BSE_INDEX%ISP(NCV13)) CYCLE
              ! k-point test
              ! IF (K1 /= 1+(NCV13-1)/((B1%VBMAX-B1%VBMIN+1)*(B1%CBMAX-B1%CBMIN+1))) CYCLE
                IF (K1 /= BSE_INDEX%NK(NCV13)) CYCLE

              ! test on NB1
              ! NB1_INTO_TOT=MOD(NCV13-1,B1%VBMAX-B1%VBMIN+1)+B1%VBMIN
                NB1_INTO_TOT=BSE_INDEX%N1(NCV13)
                N1=NB1_INTO_TOT-NPOS1+1

              ! test on NB3
              ! NB3_INTO_TOT=MOD((NCV13-1)/(B1%VBMAX-B1%VBMIN+1),(B1%CBMAX-B1%CBMIN+1))+B1%CBMIN
                NB3_INTO_TOT=BSE_INDEX%N3(NCV13)
                NB3=NB3_INTO_TOT-BD(ISP)%CBMIN+1
                N3=NB3_INTO_TOT-NPOS3+1

                IF (N1<1 .OR. N1 > NSTRIP1) CYCLE
                IF (N2<1 .OR. N2 > NSTRIP2) CYCLE
                IF (N3<1 .OR. N3 > NSTRIP3) CYCLE
                IF (N4<1 .OR. N4 > NSTRIP4) CYCLE

#ifdef check_index
                ! internal consistency test (could be removed)
              ! IF (NCV13/=NB1_INTO_TOT-B1%VBMIN+1 +(B1%VBMAX-B1%VBMIN+1)*(NB3_INTO_TOT-B1%CBMIN+(B1%CBMAX-B1%CBMIN+1)*(K1-1))) THEN
                IF (NCV13/=BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)) THEN
                   CALL vtutor%bug("internal error 2 in TWOELECTRON4O_STORE_SCALA: " // str(NCV13) &
                      // " " // str(NB1_INTO_TOT-B1%VBMIN+1) // " " // str(NB3_INTO_TOT) // " " // &
                      str(K1) // " " // str(ISP) // " " // str(BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1,&
                      NB3_INTO_TOT-B1%CBMIN+1,K1,ISP -ISP_LOW+1)) // " " // str(BSE_INDEX%ISP(NCV13)), __FILE__, __LINE__)
                ENDIF

                IF ( N1>SIZE(TWOELECTRON3O_LOCAL,1) .OR. NB3>SIZE(TWOELECTRON3O_LOCAL,2) &
                 .OR.N2>SIZE(TWOELECTRON3O_LOCAL,3) .OR. NB4>SIZE(TWOELECTRON3O_LOCAL,4)) THEN
                   CALL vtutor%bug("internal error 3 in TWOELECTRON4O_STORE_SCALA: " // str(N1) // &
                      " " // str(SIZE(TWOELECTRON3O_LOCAL,1)) // " " // str(N3) // " " // &
                      str(SIZE(TWOELECTRON3O_LOCAL,2)) // " " // str(N2) // " " // &
                      str(SIZE(TWOELECTRON3O_LOCAL,3)) // " " // str(N4) // " " // &
                      str(SIZE(TWOELECTRON3O_LOCAL,4)) // " " // str(CBMAX4), __FILE__, __LINE__)
                ENDIF
                IF ( IROW+(JCOL-1)*DESC(LLD_)>SIZE(AMAT) ) THEN
                   CALL vtutor%bug("internal error 4 in TWOELECTRON4O_STORE_SCALA: " // &
                      str(IROW+(JCOL-1))(LLD_) // " " // str(SIZE(AMAT)) // " " // str(IROW) // " " &
                      // str(JCOL) // " " // str(DESC(LLD_)), __FILE__, __LINE__)
                ENDIF
#endif
                ! finally if all test have a green light store the element in AMAT
                AMAT(IROW+(JCOL-1)*DESC(LLD_))= &
                     TWOELECTRON3O_LOCAL( N1, N3 , N2, N4)*WHF%WDES%WTKPT(1)
                AMAT(IROW+(JCOL-1)*DESC(LLD_))= AMAT(IROW+(JCOL-1)*DESC(LLD_))*SCALE_PERT*ALPHAPT

                ! add diagonal element
                IF (NCV13==NCV24 .AND. IDIAG>0) THEN

                  !! catch exception if the aufbau principle is violated (ie: E_k1,n1 > E_k3,n3)
                  IF((REAL(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)) > REAL(WHF%CELTOT(NB3_INTO_TOT,K3, ISP)))) THEN
#ifndef ACC_OFFLOAD
                      WRITE(*,*) 'internal WARNING: aufbau principle is violated'
#endif
                  ELSE
                      AMAT(IROW+(JCOL-1)*DESC(LLD_)) = AMAT(IROW+(JCOL-1)*DESC(LLD_)) &
                             -IDIAG*(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
                  ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !==========================================================================
    ! now fill in conjugated elements in lower triangle of matrix
    !==========================================================================


#define lower_triangle
#ifdef lower_triangle
    JCOL=0
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(AMAT,TWOELECTRON3O_LOCAL,BSE_INDEX,WHF) &
    !$ACC & PRIVATE(NB1_INTO_TOT,NB3_INTO_TOT,NB2_INTO_TOT,NB4_INTO_TOT,NCV13,NCV24) &
    !$ACC & PRIVATE(N1,N3,N2,N4) __IF_ASYNC__
    DO J1=1,NQ,DESC(NB_)
       DO J2=1,DESC(NB_)
          JCOL=J2+J1-1
          IF (JCOL > NQ ) CYCLE
          ! global column index NCV24 species now the row in the conjugated matrix
          NCV24=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

          ! spin test
          IF (ISP /= BSE_INDEX%ISP(NCV24)) CYCLE
          ! k-point test
          IF (K1 /= BSE_INDEX%NK(NCV24)) CYCLE

          ! test on NB1
          NB1_INTO_TOT=BSE_INDEX%N1(NCV24)
          N1=NB1_INTO_TOT-NPOS1+1
          IF (N1<1 .OR. N1 > NSTRIP1) CYCLE

          ! test on NB3
          NB3_INTO_TOT=BSE_INDEX%N3(NCV24)
          NB3=NB3_INTO_TOT-BD(ISP)%CBMIN+1
          N3=NB3_INTO_TOT-NPOS3+1

          IROW=0
          !$ACC LOOP COLLAPSE(2)
          DO I1=1,NP,DESC(MB_)
             DO I2=1,DESC(MB_)
                IROW=I1+I2-1
                  
                IF (IROW > NP) CYCLE
              ! global row index NCV13 species the column position in the conjugated matrix
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2

                ! spin test
                IF (ISP2 /= BSE_INDEX%ISP(NCV13)) CYCLE
                ! k-point test, CYCLE if not identical
                IF (K2 /= BSE_INDEX%NK(NCV13)) CYCLE

                ! test on NB2
                NB2_INTO_TOT=BSE_INDEX%N1(NCV13)
                N2=NB2_INTO_TOT-NPOS2+1
                IF (N2<1 .OR. N2 > NSTRIP2) CYCLE

                ! test on NB4_INTO_TOT
                NB4_INTO_TOT=BSE_INDEX%N3(NCV13)
                !IF (NB4_INTO_TOT<CBMIN4 .OR. NB4_INTO_TOT >CBMAX4) CYCLE
                !IF (NB4_INTO_TOT<BD(ISP2)%CBMAX .OR. NB4_INTO_TOT > BD(ISP2)%CBMAX) CYCLE

                NB4=NB4_INTO_TOT-BD(ISP2)%CBMIN+1
                N4=NB4_INTO_TOT-NPOS4+1

                IF (N1<1 .OR. N1 > NSTRIP1) CYCLE
                IF (N2<1 .OR. N2 > NSTRIP2) CYCLE
                IF (N3<1 .OR. N3 > NSTRIP3) CYCLE
                IF (N4<1 .OR. N4 > NSTRIP4) CYCLE

                IF (NCV13 > NCV24) THEN
                   AMAT(IROW+(JCOL-1)*DESC(LLD_))= &
                        GCONJG(TWOELECTRON3O_LOCAL(N1, N3, N2, N4 ))*WHF%WDES%WTKPT(1)*SCALE_PERT*ALPHAPT
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
#endif
    !$ACC EXIT DATA DELETE(TWOELECTRON3O_LOCAL) __IF_ASYNC__

    ENDIF
    ENDDO cp

    DEALLOCATE(TWOELECTRON3O_LOCAL)

    PROFILING_STOP('twoelectron4o_store_scala')
  END SUBROUTINE TWOELECTRON4O_STORE_SCALA


!***********  SUBROUTINE  TWOELECTRON4O_STORE_SCALA_TRIANGLE ***********
!> adds diagonal term and copies a stripe of BSE Hamiltonian into
!> 2D block cyclic AMAT
!
!> Here is an example how we can copy a 6x6 matrix
!> 
!> | 11 | 12 | 13 | 14 | 15 | 16 |
!> | 66 | 22 | 23 | 24 | 25 | 26 |
!> | 56 | 55 | 23 | 34 | 35 | 36 |
!> | 46 | 45 | 44 | xx | xx | xx |
!> 
!> 1. We copy the upper part of the matrix
!> 2. Revert the order of proc in grid
!> 3. Revert the order of local blocks
!> 
!> | 44 | 45 | 46 |
!> | xx | 55 | 56 |
!> | xx | xx | 66 |
!> 
!> 4. Copy the lower part
!> STRIP sizes are NSTRIPV and NSTRIPC because we pass a full stripe in NTWOELECTRON3O
!> 
!> If LTRIANGLE is .FALSE. we transpose data in local stripe and proc grid
!> and distribute data into lower triangle of the Hamiltonian
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_STORE_SCALA_TRIANGLE(WHF, ISP, ISP2, ISP_LOW, &
               K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, &
               K3, NPOS3, NSTRIP3, K4, NPOS4, NSTRIP4, &
               NTWOELECTRON3O, BD, AMAT, NCV, BSE_INDEX, IDIAG, &
               ALPHAPT, SCALE_PERT, DESC, BLOCK_INDEX, NQPOINT, NCPUS, IBSE, LTRIANGLE)

GPACC USE moffload

    USE constant
    USE main_mpi
    USE, intrinsic :: iso_c_binding
    IMPLICIT NONE

    ! passed variables
    LOGICAL :: LDIAG
    LOGICAL :: LTRIANGLE
    INTEGER :: NQPOINT    !< q-point for which excitations are calculated (usually optical e.g. q=0)
    TYPE (wavespin), INTENT(IN) :: WHF  !< wavefunction
    !INTEGER, INTENT(IN) ::  ISP, ISP2, ISP_LOW !< spins
    INTEGER, INTENT(INOUT) :: K1, K2, K3, K4 !< k-points
    !< band positions and block sizes
    INTEGER, INTENT(INOUT) :: NPOS1, NSTRIP1, NPOS2, NSTRIP2
    INTEGER, INTENT(INOUT) :: ISP, ISP2, ISP_LOW
    !< local matrix (n1,n3,n2,n4)
    GDEF, TARGET, INTENT(INOUT) ::  NTWOELECTRON3O(:,:,:,:,:)
    TYPE(banddesc), INTENT(IN) :: BD(:) !< band descriptors
#ifdef single_prec_bse
    GDEFS :: AMAT(:)
    GDEFS, POINTER :: TWOE3O_P(:,:)=>NULL()
    GDEFS, POINTER :: NTWOELECTRON3O_P(:,:)=>NULL()
    GDEFS, ALLOCATABLE :: TWOE3O_PT(:,:)
    GDEFS, ALLOCATABLE, TARGET ::  NTWOELECTRON3O_SP(:,:,:,:,:)
#endif
#ifdef double_prec_bse
    GDEF :: AMAT(:)
    GDEF, POINTER :: TWOE3O_P(:,:)=>NULL()
    GDEF, ALLOCATABLE :: TWOE3O_PT(:,:)
    GDEF, ALLOCATABLE, TARGET :: NTWOELECTRON3O_BACKUP(:,:,:,:,:)
#endif
    INTEGER :: NCV                       !< AMAT dimension NCV*NCV
    INTEGER :: NCV_COL                   !<
    INTEGER :: NCV_ROW                   !<
    INTEGER :: NCV_LCOL                  !<
    INTEGER :: NCV_LROW                  !<
    TYPE(bse_matrix_index), INTENT(IN) :: BSE_INDEX  !< AMAT indices
    INTEGER, INTENT(IN) :: IDIAG         !< add eigenvalues
    REAL(q), INTENT(IN) :: ALPHAPT       !< constant for AC integration
    REAL(q), INTENT(IN) :: SCALE_PERT    !< scaling Coulomb
    INTEGER, INTENT(IN) :: DESC(DLEN_)   !< descriptor for AMAT
    INTEGER, INTENT(IN) :: IBSE          !< add eigenvalues


    ! local variables
    !< 2d pointer to 4d TWOELECTRON3O
    !< array for corner cases
    GDEF, ALLOCATABLE, TARGET ::  TWOELECTRON3O(:,:,:,:)
    INTEGER :: NB1_INTO_TOT              !< orbital number in band structure
    INTEGER :: NB3_INTO_TOT              !< orbital number in band structure
    INTEGER :: NB2_INTO_TOT              !< orbital number in band structure
    INTEGER :: NB4_INTO_TOT              !< orbital number in band structure
    INTEGER :: NCV13, NCV24              !< indices in AMAT
    INTEGER :: K2L                       !< local K2 number
    INTEGER :: NPROCS                    !< number of procs
    INTEGER :: ITWOELECTRON3O(5)         !< dimensions of TWOELECTRON3O
    INTEGER :: IA, JA                    !< position of stripe in AMAT
    INTEGER :: IB, JB                    !< position of stripe in NTWOELECTRON3O
    INTEGER :: V1, V2, C3, C4            !< band indices
    INTEGER :: NV1, NV2, NC3, NC4        !< local matrix dimensions
    INTEGER :: NLCV13, NLCV24            !< local matrix dimensions in 2d
    INTEGER :: CV13, CV24                !< local matrix indices
    INTEGER :: ICTXC, ICTXG, INFO, LLD
    INTEGER :: NPOS3, NSTRIP3, NPOS4, NSTRIP4

    ! BLACS variables
    INTEGER :: MYROW, MYCOL, NPROW, NPCOL
    INTEGER :: DESC_L(DLEN_)
    INTEGER :: IFAIL
    INTEGER :: CBMINL

    INTEGER :: M,N
    INTEGER :: K1L

    INTEGER :: K1MIN, K1MAX
    INTEGER :: K2MIN, K2MAX
    INTEGER, ALLOCATABLE :: IMAP(:,:)
    INTEGER, ALLOCATABLE :: IMAPT(:,:)

    INTEGER :: BLOCK_INDEX, SHIFT
    INTEGER :: ROWS, NBLOCKS, NBLOCK, FIRST_LOCAL_BLOCK, LAST_LOCAL_BLOCK
    INTEGER :: FIRST_GLOBAL_BLOCK, LAST_GLOBAL_BLOCK, GLOBAL_ROWS
    INTEGER :: NROW, FIRST_ROW, LAST_ROW
    INTEGER :: FIRST_DIAGONAL_BLOCK, LAST_DIAGONAL_BLOCK
    INTEGER :: FIRST_GLOBAL_UPPER_BLOCK
    INTEGER :: FIRST_GLOBAL_LOWER_BLOCK

    INTEGER :: PC, PR, I, J, K, NK, KTOT
    INTEGER :: COLS, NCPUS
    INTEGER :: LOCAL_INDEX

    INTEGER :: BLOCK_SIZE
    INTEGER :: BOFFSET

    INTEGER :: DESC_M(DLEN_), DESCIP(DLEN_)
    INTEGER :: FIRST_GLOBAL_BCOL, LAST_GLOBAL_BCOL
    INTEGER :: FIRST_GLOBAL_ROW, LAST_GLOBAL_ROW
    INTEGER :: FIRST_GLOBAL_BROW, LAST_GLOBAL_BROW, LAST_GLOBAL_UPPER_BLOCK
    INTEGER :: FRIST_GLOBAL_LOWER_BLOCK_TRI, LAST_GLOBAL_LOWER_BLOCK_TRI
    INTEGER :: NCV_FULL, NCV_RED


    PROFILING_START('twoelectron4o_store_scala_triangle')

    NCV_FULL=SIZE(BSE_INDEX%MASK(:,1:BSE_INDEX%COLUMNS))
    NCV_RED=SUM(BSE_INDEX%MASK(:,1:BSE_INDEX%COLUMNS))

#if defined(ACC_OFFLOAD) && ! defined(CUBLASMP)
    !$ACC UPDATE SELF(AMAT,NTWOELECTRON3O)  __IF_ASYNC__
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    CALL OFFLOAD_ON_PUSH(.FALSE.)
#elif defined(ACC_OFFLOAD) && defined(CUBLASMP)
    IF(NCV_FULL /= NCV_RED) THEN
       !$ACC UPDATE SELF(AMAT,NTWOELECTRON3O)  __IF_ASYNC__
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CALL OFFLOAD_ON_PUSH(.FALSE.)
    ENDIF
#endif

    ITWOELECTRON3O=SHAPE(NTWOELECTRON3O)

    ! local matrix dimensions
    NV1=NSTRIP1
    NC3=NSTRIP3
    NV2=NSTRIP2
    NC4=NSTRIP4
    NBLOCKS=ITWOELECTRON3O(5)
    ! local matrix dimensions in 2d
    NLCV13 = NV1*NC3
    NLCV24 = NV2*NC4
    ALLOCATE(TWOELECTRON3O(NV1,NC3,NV2,NC4))

    CALL WGDSCAL(SIZE(NTWOELECTRON3O),WHF%WDES%WTKPT(1)*SCALE_PERT*ALPHAPT,NTWOELECTRON3O(1,1,1,1,1),1)

    IF(BLOCK_INDEX <= BSE_INDEX%BLOCKS .AND. IDIAG > 0) THEN
      DO NBLOCK = 1,NBLOCKS
         LOCAL_INDEX = BLOCK_INDEX + NBLOCK - 1

         IF(LOCAL_INDEX > BSE_INDEX%BLOCKS) CYCLE

         K1 = BSE_INDEX%BK1(LOCAL_INDEX)
         K2 = BSE_INDEX%BK2(LOCAL_INDEX)
         K3 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
         K4 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K2)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
         IF (K1 <= 0) CALL vtutor%bug("internal error 1 in TWOELECTRON4O_STORE_SCALA_TRIANGLE : " &
                                        // str(NCV), __FILE__, __LINE__)


         NPOS1 = BSE_INDEX%BNPOS1(LOCAL_INDEX)
         NPOS2 = BSE_INDEX%BNPOS2(LOCAL_INDEX)
         NPOS3 = BSE_INDEX%BNPOS3(LOCAL_INDEX)
         NPOS4 = BSE_INDEX%BNPOS4(LOCAL_INDEX)
         ISP = BSE_INDEX%BISP1(LOCAL_INDEX)
         ISP2 = BSE_INDEX%BISP2(LOCAL_INDEX)

         IF (K1 == K2 .AND. ISP == ISP2 .AND. NPOS1 == NPOS2 .AND. NPOS3 == NPOS4) THEN
            !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(NTWOELECTRON3O,WHF) &
            !$ACC & PRIVATE(NB1_INTO_TOT,NB3_INTO_TOT,NB2_INTO_TOT,NB4_INTO_TOT,NCV13,NCV24,&
            !$ACC & CV13,CV24) __IF_ASYNC__
            DO C3=1,NSTRIP3
               DO V1=1,NSTRIP1
                  DO C4=1,NSTRIP4
                     DO V2=1,NSTRIP2
                        ! we can get edge-blocks with fewer bands, but we don't actually know how many bands in this NBLOCK
                        ! so need to make sure that NB1_INTO_TOT doesn't exceed the limits


                        ! at potentially a segfault if fewer conduction states than core states
                        NB1_INTO_TOT=MIN(NPOS1-1+V1,BD(ISP)%VBMAXpCORE)
                        NB3_INTO_TOT=MIN(NPOS3-1+C3,BD(ISP)%CBMAX)
                        NB2_INTO_TOT=MIN(NPOS2-1+V2,BD(ISP2)%VBMAXpCORE)
                        NB4_INTO_TOT=MIN(NPOS4-1+C4,BD(ISP2)%CBMAX)

                        NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP-ISP_LOW+1)
                        NCV24=BSE_INDEX%INDEX(NB2_INTO_TOT-BD(ISP2)%VBMIN+1, NB4_INTO_TOT-BD(ISP2)%CBMIN+1, K2, ISP2-ISP_LOW+1)
                        IF (NCV13 == NCV24) THEN

                           ! add diagonal terms
                           IF (NB1_INTO_TOT > BD(ISP)%VBMAX) THEN
                              NTWOELECTRON3O(V1,C3,V2,C4,NBLOCK) = NTWOELECTRON3O(V1,C3,V2,C4,NBLOCK) - &
                                 IDIAG*(COREEIGENVALUE-WHF%CELTOT(NB3_INTO_TOT,K3,ISP))
                           ELSE
                              NTWOELECTRON3O(V1,C3,V2,C4,NBLOCK) = NTWOELECTRON3O(V1,C3,V2,C4,NBLOCK) - &
                                 IDIAG*(WHF%CELTOT(NB1_INTO_TOT,K1,ISP)-WHF%CELTOT(NB3_INTO_TOT,K3,ISP))
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
    ENDIF

    ! in TE-BSE AMAT=|A,B|, so B has to be copied with offset
    IF((IBSE==1 .OR. IBSE==3) .AND. IDIAG == 0) THEN
       BOFFSET = BSE_INDEX%NCV
    ELSE
       BOFFSET = 0
    ENDIF


    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
#ifdef double_prec_bse
    ! check if we have empty rows/columns
    ! even if there are no empty rows/columns in this interation
    ! we allocate array anyway for simplicity
    IF (NCV_FULL /= NCV_RED) THEN
       ! it is possible that by swapping columns in the upper matrix we mess up the order in the lower
       ! to avoid that we backup the NTWOELECTRON3O and copy it back before swaping lower part
       ALLOCATE(NTWOELECTRON3O_BACKUP(ITWOELECTRON3O(1),ITWOELECTRON3O(2),&
                ITWOELECTRON3O(3),ITWOELECTRON3O(4),ITWOELECTRON3O(5)))
       !$ACC ENTER DATA CREATE(NTWOELECTRON3O_BACKUP) __IF_ASYNC__
       !$ACC KERNELS PRESENT(NTWOELECTRON3O_BACKUP,NTWOELECTRON3O) __IF_ASYNC__
       NTWOELECTRON3O_BACKUP = NTWOELECTRON3O
       !$ACC END KERNELS
    ENDIF

    ! we use pointers here so that we can use the same code for sp and dp calls
    call c_f_pointer(c_loc(NTWOELECTRON3O), TWOE3O_P,  [NLCV13, NLCV24 * NBLOCKS])
#else
    ! for sp creating BACKUP is not required because we can use NTWO..._SP for backup
    ALLOCATE(NTWOELECTRON3O_SP(ITWOELECTRON3O(1),ITWOELECTRON3O(2),&
             ITWOELECTRON3O(3),ITWOELECTRON3O(4),ITWOELECTRON3O(5)))
    !$ACC ENTER DATA COPYIN(NTWOELECTRON3O_SP) __IF_ASYNC__
    !$ACC KERNELS PRESENT(NTWOELECTRON3O_SP,NTWOELECTRON3O) __IF_ASYNC__
    NTWOELECTRON3O_SP=NTWOELECTRON3O
    !$ACC END KERNELS

    call c_f_pointer(c_loc(NTWOELECTRON3O_SP), TWOE3O_P,  [NLCV13, NLCV24 * NBLOCKS])
#endif


    ! the first local block - node index * number of blocks per node
    FIRST_GLOBAL_BLOCK = BLOCK_INDEX - (WHF%WDES%COMM%NODE_ME-1)*NBLOCKS
    ! cannot exceed the total number of blocks, in case more cpus than blocks is used
    LAST_GLOBAL_BLOCK = MIN(BSE_INDEX%BLOCKS, FIRST_GLOBAL_BLOCK + NCPUS*NBLOCKS - 1)

    FIRST_ROW = CEILING(REAL(FIRST_GLOBAL_BLOCK,KIND=q) /REAL(BSE_INDEX%COLUMNS,KIND=q))
    LAST_ROW = CEILING(REAL(LAST_GLOBAL_BLOCK, KIND=q) / REAL(BSE_INDEX%COLUMNS,KIND=q))
    ROWS = CEILING(REAL(BSE_INDEX%COLUMNS,KIND=q)/2._q)

    FIRST_DIAGONAL_BLOCK = (FIRST_ROW-1)*(BSE_INDEX%COLUMNS+1) + 1
    LAST_DIAGONAL_BLOCK = MIN(BSE_INDEX%BLOCKS, BSE_INDEX%COLUMNS*(LAST_ROW - 1) + LAST_ROW - 1)

    SHIFT=0

    CALL BLACS_GET(0, 0, ICTXG)
    CALL BLACS_GRIDINIT(ICTXG, 'R', 1, COMM%NCPU)

    ! Copying the upper trapeziod
    IF( FIRST_ROW <= ROWS .AND. LAST_GLOBAL_BLOCK >= FIRST_DIAGONAL_BLOCK) THEN
       PC = MIN(NCPUS,BSE_INDEX%COLUMNS / NBLOCKS)
       PR = MAX(1,(LAST_GLOBAL_BLOCK - FIRST_GLOBAL_BLOCK + 1) / BSE_INDEX%COLUMNS)


       CALL BLACS_GET(0, 0, ICTXC)
       CALL BLACS_GRIDINIT(ICTXC, 'R', PR, PC)
       CALL BLACS_GRIDINFO(ICTXC, NPROW, NPCOL, MYROW, MYCOL)

       ! descriptor of distributed stripe (NLCV13,NLCV24*NPCOL)
       IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
         CALL DESCINIT(DESC_L, NLCV13*PR, NLCV24*NBLOCKS*PC, NLCV13, NLCV24*NBLOCKS, &
                0, 0, ICTXC, NLCV13, IFAIL)
       ELSE
         DESC_L(2)=-1
       ENDIF

       IF(FIRST_GLOBAL_BLOCK >= FIRST_DIAGONAL_BLOCK) THEN
          ! passed diagonal element
          SHIFT = 0
          IB = 1
          JB = 1
       ELSE
          ! not yet passed diagonal, need to shift JB to first upper matrix element
          ! either the last element in the row or the last global element
          SHIFT = FIRST_DIAGONAL_BLOCK - FIRST_GLOBAL_BLOCK
          IB = 1
          JB=SHIFT*NLCV24+1
       ENDIF

       ! shift the index of the first block
       FIRST_GLOBAL_UPPER_BLOCK = FIRST_GLOBAL_BLOCK + SHIFT

       K1 = BSE_INDEX%BK1(FIRST_GLOBAL_UPPER_BLOCK)
       K2 = BSE_INDEX%BK2(FIRST_GLOBAL_UPPER_BLOCK)
       NPOS1 = BSE_INDEX%BNPOS1(FIRST_GLOBAL_UPPER_BLOCK)
       NPOS2 = BSE_INDEX%BNPOS2(FIRST_GLOBAL_UPPER_BLOCK)
       NPOS3 = BSE_INDEX%BNPOS3(FIRST_GLOBAL_UPPER_BLOCK)
       NPOS4 = BSE_INDEX%BNPOS4(FIRST_GLOBAL_UPPER_BLOCK)
       ISP = BSE_INDEX%BISP1(FIRST_GLOBAL_UPPER_BLOCK)
       ISP2 = BSE_INDEX%BISP2(FIRST_GLOBAL_UPPER_BLOCK)

       IA=BSE_INDEX%INDEX(NPOS1-BD(ISP)%VBMIN+1, NPOS3-BD(ISP)%CBMIN+1, K1, ISP-ISP_LOW+1)
       JA=BSE_INDEX%INDEX(NPOS2-BD(ISP2)%VBMIN+1, NPOS4-BD(ISP2)%CBMIN+1, K2, ISP2-ISP_LOW+1)

       FIRST_GLOBAL_BROW = CEILING(REAL(FIRST_GLOBAL_UPPER_BLOCK,KIND=q) /REAL(BSE_INDEX%COLUMNS,KIND=q))
       LAST_GLOBAL_BROW = FIRST_GLOBAL_BROW + PR - 1

       FIRST_GLOBAL_BCOL = MOD(FIRST_GLOBAL_UPPER_BLOCK-1, BSE_INDEX%COLUMNS)+1

       LAST_GLOBAL_UPPER_BLOCK= CEILING(REAL(BSE_INDEX%COLUMNS,KIND=q)/2._q)*BSE_INDEX%COLUMNS
       LAST_GLOBAL_BCOL = MOD(MIN(LAST_GLOBAL_UPPER_BLOCK-1,LAST_GLOBAL_BLOCK-1), BSE_INDEX%COLUMNS)+1

       NCV_FULL=SIZE(BSE_INDEX%MASK(:,FIRST_GLOBAL_BCOL:LAST_GLOBAL_BCOL))
       NCV_RED=SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BCOL:LAST_GLOBAL_BCOL))
       IF (NCV_FULL /= NCV_RED) THEN
         CALL DEFLATOR('C',BSE_INDEX%MASK,FIRST_GLOBAL_BCOL,LAST_GLOBAL_BCOL, SHIFT*NLCV24, TWOE3O_P, DESC_L)
       ENDIF

       NCV_FULL=SIZE(BSE_INDEX%MASK(:,FIRST_GLOBAL_BROW:LAST_GLOBAL_BROW))
       NCV_RED=SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BROW:LAST_GLOBAL_BROW))
       IF (NCV_FULL /= NCV_RED) THEN
         CALL DEFLATOR('R',BSE_INDEX%MASK,FIRST_GLOBAL_BROW,LAST_GLOBAL_BROW, 0, TWOE3O_P, DESC_L)
       ENDIF

       ! these NCVs can be smaller than NLCV13*PR, because NLCV13 is just the number of pairs, but some pars can be excluded
       NCV_COL = SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BCOL:LAST_GLOBAL_BCOL))
       NCV_ROW = SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BROW:LAST_GLOBAL_BROW))

       !$ACC ENTER DATA COPYIN(TWOE3O_P) __IF_ASYNC__
       IF(FIRST_GLOBAL_UPPER_BLOCK == FIRST_DIAGONAL_BLOCK ) THEN
          CALL PGTRMR2Dbse('U','N', NCV_ROW, NCV_COL, TWOE3O_P, IB, JB, DESC_L, AMAT, IA, JA+BOFFSET, DESC, ICTXG)
       ELSE
          CALL PGGEMR2Dbse(NCV_ROW, NCV_COL, TWOE3O_P, IB, JB, DESC_L, AMAT, IA, JA+BOFFSET, DESC, ICTXG)
       ENDIF
       !$ACC EXIT DATA DELETE(TWOE3O_P) __IF_ASYNC__

       ! close the context
       IF(MYROW /= -1 .OR. MYCOL /= -1) CALL BLACS_GRIDEXIT(ICTXC)

       ! copy stripe into lower triangle
       IF(.NOT. LTRIANGLE) THEN
          ALLOCATE(IMAP(PC,PR))
          ! set up transposed proc grid
          DO I = 1, PC
             DO J = 1, PR
                IMAP(I, J) = I + (J-1)*PC - 1
             ENDDO
          ENDDO

          CALL BLACS_GET(0, 0, ICTXC)
          ! here we swap PC and PR compared to the case of upper triangle
          CALL BLACS_GRIDMAP(ICTXC, IMAP, PC, PC, PR)
          CALL BLACS_GRIDINFO(ICTXC, NPROW, NPCOL, MYROW, MYCOL)

          ! descriptor of distributed stripe (NLCV24*NPCOL,NLCV13)
          IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
             CALL DESCINIT(DESC_L, NLCV24*NBLOCKS*PC,NLCV13*PR, NLCV24*NBLOCKS, NLCV13,  &
                0, 0, ICTXC, NLCV24*NBLOCKS, IFAIL)
          ELSE
             DESC_L(2)=-1
          ENDIF

          ! allocate array for transposing the stripe
          ALLOCATE(TWOE3O_PT(SIZE(TWOE3O_P,2),SIZE(TWOE3O_P,1)))
          !$ACC ENTER DATA CREATE(TWOE3O_PT) COPYIN(TWOE3O_P) __IF_ASYNC__

          !$ACC KERNELS PRESENT(TWOE3O_PT,TWOE3O_P) __IF_ASYNC__
          TWOE3O_PT= TRANSPOSE(TWOE3O_P)
          TWOE3O_PT= GCONJG(TWOE3O_PT)
          !$ACC END KERNELS

          ! swaped JB, IB and JA, IA as compared to when we distribute into upper triangle
          IF(FIRST_GLOBAL_UPPER_BLOCK == FIRST_DIAGONAL_BLOCK ) THEN
             CALL PGTRMR2Dbse('L','U', NCV_COL, NCV_ROW, TWOE3O_PT, JB, IB, DESC_L, AMAT, JA, IA+BOFFSET, DESC, ICTXG)
          ELSE
             CALL PGGEMR2Dbse(NCV_COL, NCV_ROW, TWOE3O_PT, JB, IB, DESC_L, AMAT, JA, IA+BOFFSET, DESC, ICTXG)
          ENDIF
          !$ACC EXIT DATA DELETE(TWOE3O_PT,TWOE3O_P) __IF_ASYNC__

          DEALLOCATE(TWOE3O_PT)

          DEALLOCATE(IMAP)

          IF(MYROW /= -1 .OR. MYCOL /= -1) CALL BLACS_GRIDEXIT(ICTXC)
       ENDIF

       CALLMPI( M_barrier( WHF%WDES%COMM))
    ENDIF


    ! lower trapezoid
    IF (LAST_ROW > 1 .AND. FIRST_GLOBAL_BLOCK <= LAST_DIAGONAL_BLOCK) THEN
#ifdef double_prec_bse
       ! this step is not required for sp becuase we always copy 
       ! NTWOELECTRON3O_SP to NTWOELECTRON3O anyway
       NCV_FULL=SIZE(BSE_INDEX%MASK(:,1:BSE_INDEX%COLUMNS))
       NCV_RED=SUM(BSE_INDEX%MASK(:,1:BSE_INDEX%COLUMNS))
       IF (NCV_FULL /= NCV_RED) THEN
          ! it is possible that by swapping columns in the upper matrix we mess up the order in the lower
          ! to avoid that we swap upper and lower triangle in separate arrays
          NTWOELECTRON3O = NTWOELECTRON3O_BACKUP
       ENDIF
#endif


       ! invert local blocks
       !$ACC ENTER DATA CREATE(TWOELECTRON3O) __IF_ASYNC__
       DO NBLOCK = 1, NBLOCKS/2
          BLOCK_SIZE=SIZE(NTWOELECTRON3O(:,:,:,:,NBLOCK))
          CALL GCOPYbse(BLOCK_SIZE,NTWOELECTRON3O(:,:,:,:,NBLOCK), 1, TWOELECTRON3O, 1)
          CALL GCOPYbse(BLOCK_SIZE,NTWOELECTRON3O(:,:,:,:,NBLOCKS-NBLOCK+1),1,NTWOELECTRON3O(:,:,:,:,NBLOCK), 1)
          CALL GCOPYbse(BLOCK_SIZE,TWOELECTRON3O,1,NTWOELECTRON3O(:,:,:,:,NBLOCKS-NBLOCK+1), 1)
       ENDDO
       !$ACC EXIT DATA DELETE(TWOELECTRON3O) __IF_ASYNC__

#ifdef single_prec_bse
       !$ACC KERNELS PRESENT(NTWOELECTRON3O_SP,NTWOELECTRON3O) __IF_ASYNC__
       NTWOELECTRON3O_SP=NTWOELECTRON3O
       !$ACC END KERNELS
#endif
       ! number of rows
       PR = LAST_ROW - FIRST_ROW  + 1

       ! there are no lower half blocks in the first row
       IF(FIRST_ROW == 1) PR = PR -1

       ! number of columns
       PC = MIN(NCPUS,BSE_INDEX%COLUMNS / NBLOCKS)

       ALLOCATE(IMAP(PR,PC))

       ! offset the cpu index if not all cpus are used
       SHIFT = NCPUS - WHF%WDES%COMM%NCPU

       ! invert the order of cpus
       DO I = 1, PR
          DO J = 1, PC
             IMAP(I, J) = MOD(SHIFT+WHF%WDES%COMM%NCPU  - J - (I-1)*PC, WHF%WDES%COMM%NCPU)
          ENDDO
       ENDDO

       CALL BLACS_GET(0, 0, ICTXC)
       CALL BLACS_GRIDMAP(ICTXC, IMAP, PR, PR, PC)
       CALL BLACS_GRIDINFO(ICTXC, NPROW, NPCOL, MYROW, MYCOL)

       ! descriptor of distributed stripe (NLCV13,NLCV24*NPCOL)
       IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
         CALL DESCINIT(DESC_L, NLCV13*PR, NLCV24*NBLOCKS*PC, NLCV13, NLCV24*NBLOCKS, &
                0, 0, ICTXC, NLCV13, IFAIL)
       ELSE
         DESC_L(2)=-1
       ENDIF


       SHIFT = 0
       IF(LAST_GLOBAL_BLOCK >= LAST_DIAGONAL_BLOCK) THEN
          ! passed diagonal and need to shift back
          SHIFT = FIRST_GLOBAL_BLOCK + NCPUS*NBLOCKS - LAST_DIAGONAL_BLOCK  - 1
          IB = 1
          JB=SHIFT*NLCV24+1
       ELSE
          ! not yet passed the diagonal
          IB = 1
          JB = 1
       ENDIF


       FIRST_GLOBAL_LOWER_BLOCK = MIN(LAST_DIAGONAL_BLOCK,LAST_GLOBAL_BLOCK + SHIFT)

       K1 = BSE_INDEX%BK1(FIRST_GLOBAL_LOWER_BLOCK)
       K2 = BSE_INDEX%BK2(FIRST_GLOBAL_LOWER_BLOCK)
       NPOS1 = BSE_INDEX%BNPOS1(FIRST_GLOBAL_LOWER_BLOCK)
       NPOS2 = BSE_INDEX%BNPOS2(FIRST_GLOBAL_LOWER_BLOCK)
       NPOS3 = BSE_INDEX%BNPOS3(FIRST_GLOBAL_LOWER_BLOCK)
       NPOS4 = BSE_INDEX%BNPOS4(FIRST_GLOBAL_LOWER_BLOCK)
       ISP = BSE_INDEX%BISP1(FIRST_GLOBAL_LOWER_BLOCK)
       ISP2 = BSE_INDEX%BISP2(FIRST_GLOBAL_LOWER_BLOCK)

       IA=BSE_INDEX%INDEX(NPOS1-BD(ISP)%VBMIN+1, NPOS3-BD(ISP)%CBMIN+1, K1, ISP-ISP_LOW+1)
       JA=BSE_INDEX%INDEX(NPOS2-BD(ISP2)%VBMIN+1, NPOS4-BD(ISP2)%CBMIN+1, K2, ISP2-ISP_LOW+1)

       FRIST_GLOBAL_LOWER_BLOCK_TRI = BSE_INDEX%BMAP(FIRST_GLOBAL_LOWER_BLOCK)
       LAST_GLOBAL_LOWER_BLOCK_TRI = BSE_INDEX%BMAP(LAST_GLOBAL_BLOCK)

       FIRST_GLOBAL_BROW = CEILING(REAL(FRIST_GLOBAL_LOWER_BLOCK_TRI,KIND=q) /REAL(BSE_INDEX%COLUMNS,KIND=q))
       LAST_GLOBAL_BROW = FIRST_GLOBAL_BROW + PR - 1

       FIRST_GLOBAL_BCOL = MOD(FRIST_GLOBAL_LOWER_BLOCK_TRI-1, BSE_INDEX%COLUMNS)+1
       LAST_GLOBAL_BCOL = FIRST_GLOBAL_BCOL +  PC*NBLOCKS -  SHIFT - 1

       NCV_FULL=SIZE(BSE_INDEX%MASK(:,FIRST_GLOBAL_BCOL:LAST_GLOBAL_BCOL))
       NCV_RED=SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BCOL:LAST_GLOBAL_BCOL))

       IF (NCV_FULL /= NCV_RED) THEN
         CALL DEFLATOR('C',BSE_INDEX%MASK,FIRST_GLOBAL_BCOL,LAST_GLOBAL_BCOL, SHIFT*NLCV24, TWOE3O_P, DESC_L)
       ENDIF

       NCV_FULL=SIZE(BSE_INDEX%MASK(:,FIRST_GLOBAL_BROW:LAST_GLOBAL_BROW))
       NCV_RED=SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BROW:LAST_GLOBAL_BROW))

       IF (NCV_FULL /= NCV_RED) THEN
         CALL DEFLATOR('R',BSE_INDEX%MASK,FIRST_GLOBAL_BROW,LAST_GLOBAL_BROW, 0, TWOE3O_P, DESC_L)
       ENDIF

       ! these NCVs can be smaller than NLCV13*PR, because NLCV13 is just the number of pairs, but some pars can be excluded
       NCV_COL = SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BCOL:LAST_GLOBAL_BCOL))
       NCV_ROW = SUM(BSE_INDEX%MASK(:,FIRST_GLOBAL_BROW:LAST_GLOBAL_BROW))

       !$ACC ENTER DATA COPYIN(TWOE3O_P) __IF_ASYNC__
       IF(FIRST_GLOBAL_LOWER_BLOCK == LAST_DIAGONAL_BLOCK) THEN
          CALL PGTRMR2Dbse('U','N',NCV_ROW,NCV_COL,TWOE3O_P,IB,JB,DESC_L,AMAT,IA,JA+BOFFSET,DESC,ICTXG)
       ELSE
          CALL PGGEMR2Dbse(NCV_ROW,NCV_COL,TWOE3O_P,IB,JB,DESC_L,AMAT,IA,JA+BOFFSET,DESC,ICTXG)
       ENDIF
       !$ACC EXIT DATA DELETE(TWOE3O_P) __IF_ASYNC__

       ! close the context
       IF(MYROW /= -1 .OR. MYCOL /= -1) CALL BLACS_GRIDEXIT(ICTXC)

       ! copy stripe into lower triangle
       IF(.NOT. LTRIANGLE) THEN

          ! transpose proc grid
          ALLOCATE(IMAPT(PC,PR))
          IMAPT=TRANSPOSE(IMAP)

          CALL BLACS_GET(0, 0, ICTXC)
          CALL BLACS_GRIDMAP(ICTXC, IMAPT, PC, PC, PR)
          CALL BLACS_GRIDINFO(ICTXC, NPROW, NPCOL, MYROW, MYCOL)

          ! descriptor of distributed stripe
          IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
             CALL DESCINIT(DESC_L, NLCV24*NBLOCKS*PC,NLCV13*PR, NLCV24*NBLOCKS, NLCV13,  &
                0, 0, ICTXC, NLCV24*NBLOCKS, IFAIL)
          ELSE
             DESC_L(2)=-1
          ENDIF

          ALLOCATE(TWOE3O_PT(SIZE(TWOE3O_P,2),SIZE(TWOE3O_P,1)))
          !$ACC ENTER DATA CREATE(TWOE3O_PT) COPYIN(TWOE3O_P) __IF_ASYNC__

          !$ACC KERNELS PRESENT(TWOE3O_PT,TWOE3O_P) __IF_ASYNC__
          TWOE3O_PT= TRANSPOSE(TWOE3O_P)
          TWOE3O_PT= GCONJG(TWOE3O_PT)
          !$ACC END KERNELS

          ! swaped JB, IB and JA, IA as compared to when we distribute into upper triangle
          IF(FIRST_GLOBAL_UPPER_BLOCK == FIRST_DIAGONAL_BLOCK ) THEN
             CALL PGTRMR2Dbse('L','U', NCV_COL, NCV_ROW, TWOE3O_PT, JB, IB, DESC_L, AMAT, JA, IA+BOFFSET, DESC, ICTXG)
          ELSE
             CALL PGGEMR2Dbse(NCV_COL, NCV_ROW, TWOE3O_PT, JB, IB, DESC_L, AMAT, JA, IA+BOFFSET, DESC, ICTXG)
          ENDIF
          !$ACC EXIT DATA DELETE(TWOE3O_P,TWOE3O_PT) __IF_ASYNC__

          DEALLOCATE(TWOE3O_PT)
          DEALLOCATE(IMAPT)
          IF(MYROW /= -1 .OR. MYCOL /= -1) CALL BLACS_GRIDEXIT(ICTXC)
       ENDIF
       DEALLOCATE(IMAP)

       CALLMPI( M_barrier( WHF%WDES%COMM))
    ENDIF

    CALL BLACS_GRIDEXIT(ICTXG)

    DEALLOCATE(TWOELECTRON3O)

#ifdef double_prec_bse
    IF(ALLOCATED(NTWOELECTRON3O_BACKUP)) DEALLOCATE(NTWOELECTRON3O_BACKUP)
#else
    !$ACC EXIT DATA DELETE(NTWOELECTRON3O_SP) __IF_ASYNC__
    DEALLOCATE(NTWOELECTRON3O_SP)
#endif

    NCV_FULL=SIZE(BSE_INDEX%MASK(:,1:BSE_INDEX%COLUMNS))
    NCV_RED=SUM(BSE_INDEX%MASK(:,1:BSE_INDEX%COLUMNS))

#if defined(ACC_OFFLOAD) && ! defined(CUBLASMP)
    CALL OFFLOAD_ON_POP
    !$ACC UPDATE DEVICE(AMAT)  __IF_ASYNC__
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
#elif defined(ACC_OFFLOAD) && defined(CUBLASMP)
    IF(NCV_FULL /= NCV_RED) THEN
    CALL OFFLOAD_ON_POP
    !$ACC UPDATE DEVICE(AMAT)  __IF_ASYNC__
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    ENDIF
#endif
    TWOE3O_P=>NULL()

    PROFILING_STOP('twoelectron4o_store_scala_triangle')
  END SUBROUTINE TWOELECTRON4O_STORE_SCALA_TRIANGLE

!****************** SUBROUTINE  DEFLATOR ******************************
!
!> moves empty columns/rows outwards and the data inwards
!> preserving the original order of columns/rows
!
!> The step is required if the emptry rows/column are present in
!> the calcualted block.
!> The pivoting vector must be align with the matrix distribution.
!> Each block of pivoting vector is augmented by scratch space of equal size.
!
!> Becuase of scratch space at the end offesting PZLAPV2 seems to be broken.
!> Hence, we don't offset.
!
!**********************************************************************
  SUBROUTINE DEFLATOR(COLROW,MASK,INDX_START,INDX_END,OFFSET,TWOE3O_P,DESC)
    IMPLICIT NONE

    ! passed variables
    CHARACTER :: COLROW     !< swap columns or rows
    INTEGER :: MASK(:,:)    !< mask 1 keep, 0 remove
    INTEGER :: INDX_START, INDX_END !< blocks in the mask
    INTEGER :: DESC(9)      !< descirptor of the matrix
    INTEGER :: OFFSET       !< offset, for elements of lower half
#ifdef single_prec_bse
    GDEFS :: TWOE3O_P(:,:)  !< the matrix
#endif
#ifdef double_prec_bse
    GDEF :: TWOE3O_P(:,:)   !< the matrix
#endif

    ! local variables
    INTEGER, ALLOCATABLE :: GMASK(:) !< relevant slice of MASK
    INTEGER, ALLOCATABLE :: IPIV(:)  !< pivoting vector
    INTEGER :: DESCIP(9)             !< descriptor for pivoting vector
    INTEGER :: IFAIL
    INTEGER :: NONES, I              !< counters
    INTEGER :: INDXPC, INDXPR, INDL
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL
    INTEGER, EXTERNAL    ::     NUMROC, INDXG2L, INDXL2G, INDXG2P

    CALL BLACS_GRIDINFO(DESC(2), NPROW, NPCOL, MYROW, MYCOL)

    ! proc out of context should skip this step
    IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
      IF(COLROW == 'C') THEN

         ALLOCATE(GMASK(DESC(4)))
         ! pivoting vector is augmented by scratch space of equal size at the end
         ALLOCATE(IPIV(DESC(6)*2))

         GMASK=RESHAPE(MASK(:,INDX_START:INDX_END), [SIZE(MASK(:,INDX_START:INDX_END))])

         ! total size accounts for scratch space, but the block does not inlude the scratch space
         CALL DESCINIT(DESCIP, 1, DESC(4)*2, 1, DESC(6), 0, 0, DESC(2), 1, IFAIL)


         ! initialize to -1 so that scalapack could catch errors
         IPIV=-1
         NONES=0
         DO I = 1,DESC(4)
            ! find proc for I
            INDXPC = INDXG2P(I, DESCIP(6), 0, 0, NPCOL)
            ! find loc index for I
            INDL=INDXG2L(I, DESCIP(6), 0, 0, NPCOL)
            ! index of self (swap I and I)
            IF(MYCOL == INDXPC) IPIV(INDL)=I

            ! don't swap columns before OFFSET
            IF(I <= OFFSET) CYCLE
            ! swap local INDL with global IPIV(INDL)
            IF (GMASK(I-OFFSET)/=0) THEN
               NONES=NONES + 1
               IF(MYCOL == INDXPC) IPIV(INDL)=NONES+OFFSET
            ENDIF
         ENDDO

         PROFILING_START('PGLAPV2')
         CALL PGLAPV2bse('F', 'C', DESC(3), DESC(4), TWOE3O_P(1,1), 1, 1, DESC, IPIV, 1, 1, DESCIP)
         PROFILING_STOP('PGLAPV2')
         DEALLOCATE(IPIV, GMASK)
      ENDIF

      ! see comment above
      IF(COLROW == 'R') THEN
         ALLOCATE(GMASK(DESC(3)))
         ALLOCATE(IPIV(DESC(5)*2))

         GMASK=RESHAPE(MASK(:,INDX_START:INDX_END), [SIZE(MASK(:,INDX_START:INDX_END))])

         CALL DESCINIT(DESCIP, DESC(3)*2, 1, DESC(5), 1, 0, 0, DESC(2), DESC(5)*2, IFAIL)

         IPIV=-1
         NONES=0
         DO I = 1,DESC(3)
            INDXPR = INDXG2P(I, DESCIP(5), 0, 0, NPROW)
            INDL=INDXG2L(I, DESCIP(5), 0, 0, NPROW)
            IF(MYROW == INDXPR) IPIV(INDL)=I

            IF(I <= OFFSET) CYCLE
            IF (GMASK(I-OFFSET)/=0) THEN
               NONES=NONES + 1
               IF(MYROW == INDXPR) IPIV(INDL)=NONES+OFFSET
            ENDIF
         ENDDO

         PROFILING_START('PGLAPV2')
         CALL PGLAPV2bse('F', 'R', DESC(3), DESC(4), TWOE3O_P(1,1), 1, 1, DESC, IPIV, 1, 1, DESCIP)
         PROFILING_STOP('PGLAPV2')
         DEALLOCATE(IPIV, GMASK)
      ENDIF
   ENDIF


  END SUBROUTINE DEFLATOR


!**********************************************************************
!
!> ScaLAPACK-aware subroutine for calculating the BSE oscillator strength.
!>
!> AMAT stores the eigenvectors of the BSE Hamiltonian.
!> The BSE Hamiltonian is calculated as
!> ~~~
!> Hamiltonian(I, J) = < c_k3,n3 v'_k2,n2 | v+W | v_k1,n1  c'_k4,n4 >
!> with I = (k1,n1, k3, n3)   J= (k2,n2, k4, n4)
!> ~~~
!> after diagonalization AMAT stores U(I, lambda) with
!> ~~~
!> Hamiltonian(I, J)= U Lambda U+ =
!>        \sum_lambda U(I, lambda) lambda U(J, lambda)*
!> ~~~
!> (this has been carefully checked by dumping matrix elements).
!>
!> The oscillator strength is given by
!> ~~~
!> \sum_lambda \sum_I <v_k1,n1|O|c_k3,n3> U(I, lambda)  1/ (w- lambda )
!>             [sum_J U( J, lambda) <v_k2,n2|O|c_k4,n4>]*
!> ~~~
!> where O is an arbitrary operator.
!> To obtain the polarizabilty at a wavevector q=k3-k1
!> ~~~
!>  < c_k3,n3| e^i(k3-k1 ) | v_k1,n1 >
!> ~~~
!
!**********************************************************************

  SUBROUTINE CALCULATE_BSE_OSZI_SCALA (WHF, NQPOINT, KPOINT_BSE, LATT_CUR, ISP_LOW, ISP_HIGH,  &
       SHIFT, BD, BSE_INDEX, R, AMAT, NCV, ANTIRES, NOMEGA_DIM, OMEGAMAX, IU6, DESC)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    USE mymath
    USE core_con_mat
    USE mlr_optic, ONLY: CDER_CORECON_Q
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER :: NQPOINT        !< q-point considered
    INTEGER :: KPOINT_BSE(3)  !< G-vector considered
    TYPE (latt) LATT_CUR
    INTEGER ISP_LOW, ISP_HIGH
    REAL(q) SHIFT             !< complex shift for broadening
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: R(:)
#ifdef single_prec_bse
    GDEFS :: AMAT(:)
#endif

#ifdef double_prec_bse
    GDEF  :: AMAT(:)
#endif
    INTEGER :: NCV
    INTEGER :: ANTIRES
    INTEGER :: NOMEGA_DIM
    REAL(q) :: OMEGAMAX
    INTEGER :: IU6
    INTEGER :: DESC(DLEN_)
 ! local
    INTEGER :: NCV13
    INTEGER K1, K1_LOOKED_UP, K3, ISP
    INTEGER N1, N3
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    COMPLEX(q) :: WEIGHT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG
    INTEGER     :: LAMBDA, NOMEGA
    GDEF        :: CDER_BETWEEN_STATE(3)
    COMPLEX(q)  :: CTRANS(NCV,3)
    COMPLEX(q)  :: CTRANS_SQUARE(3,3)
    REAL(q)     :: DKX, DKY, DKZ
    COMPLEX(q)  :: EPS(NOMEGA_DIM,3,3)
    REAL(q)     :: OMEGA(NOMEGA_DIM), DOMEGA
    REAL(q)     :: R_AND_INTENSITY(2,SIZE(R))
    REAL(q)     :: VKPT(3)
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP,NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2, NLAMBDA
    REAL(q)     :: RSPIN
    COMPLEX(q)  :: R_CTRANS1,R_CTRANS2,R_CTRANS3
#ifdef VASP_HDF5
    CHARACTER(LEN=:), ALLOCATABLE :: PREFIX
#endif
    !$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX,KPOINT_IN_FULL_GRID) SEQ
    !$ACC ROUTINE(FILLED_XI_ORBITAL,EMPTY_XI_ORBITAL,CDER_BETWEEN_STATES_ROTATED) SEQ

    PROFILING_START('calculate_bse_oszi_scala')

!  three cases need to be considered
!  non-collinear RSPIN=1.0
    IF (WHF%WDES%LNONCOLLINEAR) THEN
       RSPIN=1.0_q
    ELSE
!  collinear non spin polarized: RSPIN=2.0
!  collinear spin polarized:     RSPIN=1.0 for  ISP_LOW=1 and ISP_HIGH=2
!  I believe the case ISP_HIGH=ISP_LOW for collinear spin polarized makes little sense
!  but it is not used right now anyway
       RSPIN=2._q/(ISP_HIGH-ISP_LOW+1)
    ENDIF

    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    EPS = 0
    K1_LOOKED_UP=-1
    NLAMBDA=0

    DOMEGA=R(NCV)/(NOMEGA_DIM-1)

    IF (OMEGAMAX>0) THEN
       DOMEGA=OMEGAMAX/(NOMEGA_DIM-1)
    ENDIF

    !$ACC UPDATE DEVICE(CDER_BETWEEN_STATES,CDER_CORECON)  __IF_ASYNC__
    !$ACC ENTER DATA CREATE(CTRANS,VKPT)  __IF_ASYNC__
    !$ACC ENTER DATA COPYIN(KPOINTS_FULL_ORIG,KPOINTS_FULL) __IF_ASYNC__

    !$ACC KERNELS PRESENT(CTRANS,VKPT) __IF_ASYNC__
    CTRANS=0
    VKPT=0._q
    !$ACC END KERNELS

    IF(MYCOL /= -1 .OR. MYROW /= -1) THEN

   !   loop over column indices (eigenvectors)
       !$ACC PARALLEL LOOP GANG PRESENT(WHF,WHF%WDES,WHF%FERTOT,WHF%WDES%VKPT,WHF%WDES%WTKPT, &
       !$ACC & CTRANS,DESC,BSE_INDEX,BSE_INDEX%ISP,BSE_INDEX%NK,BSE_INDEX%N1,BSE_INDEX%N3,BSE_INDEX%INDEX, &
       !$ACC &KPOINTS_FULL,KPOINTS_FULL_ORIG,LATT_CUR,AMAT) PRIVATE(JCOL,J2,J1,R_CTRANS2,R_CTRANS3,R_CTRANS1) __IF_ASYNC__
       DO LAMBDA=1, NCV
          ! local storage index in block (J1) and local index of block (J2)
          J1=   ((LAMBDA-DESC(NB_)*MYCOL-1)/(NPCOL*DESC(NB_)))*DESC(NB_)+1
          J2=MOD((LAMBDA-DESC(NB_)*MYCOL-1) ,NPCOL*DESC(NB_))+1
   
          JCOL=J1+J2-1
          IF (J2 <= 0 .OR. J2 > DESC(NB_)) THEN
             JCOL  =-1
          ELSE
             NLAMBDA=NLAMBDA+1
#ifndef ACC_OFFLOAD
             IF (LAMBDA /= DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2) THEN
                CALL vtutor%bug("internal error in CALCULATE_BSE_OSZI_SCALA: " // str(LAMBDA) // " " &
                   // str(DESC(NB_)) // " " // str(NPCOL*(J1-1)) // " " // str(J2), __FILE__, __LINE__)
             ENDIF
#endif
          ENDIF
          IF (JCOL==-1) THEN
             ! no data on local node for eigenvector/eigenvalue pair LAMBDA_DONE
             ! simply clear CTRANS
          ELSE
             ! local data then calculate contribution to CTRANS
             IROW=0
             R_CTRANS1=0
             R_CTRANS2=0
             R_CTRANS3=0
             !$ACC LOOP COLLAPSE(2) VECTOR REDUCTION(+:R_CTRANS1,R_CTRANS2,R_CTRANS3) &
             !$ACC & PRIVATE(CDER_BETWEEN_STATE,IROW,NCV13,ISP,K1,NK1_IN_KPOINTS_FULL_ORIG,&
             !$ACC & K3,NB1_INTO_TOT,NB3_INTO_TOT,WEIGHT,VKPT)
             DO I1=1,NP,DESC(MB_)
                DO I2=1,DESC(MB_)
                   IROW=I2+I1-1
                   IF(IROW > NP) CYCLE
                   NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                   ISP= BSE_INDEX%ISP(NCV13)
                 ! K1 = 1+(NCV13-1)/((BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)*(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1))
                   K1 = BSE_INDEX%NK(NCV13)
                   IF (WHF%WDES%WTKPT(K1)==0) CYCLE
#ifdef ACC_OFFLOAD
                   NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)
                   VKPT=WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT)
                   K3 =KPOINT_IN_FULL_GRID(VKPT,KPOINTS_FULL)
#else
                   IF (K1/= K1_LOOKED_UP) THEN
                      NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)
                      K1_LOOKED_UP=K1
                      K3 =KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
                   ENDIF
#endif
                   NB1_INTO_TOT=BSE_INDEX%N1(NCV13)
                   NB3_INTO_TOT=BSE_INDEX%N3(NCV13)
#ifndef ACC_OFFLOAD
                   IF (NCV13/=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)) THEN
                      CALL vtutor%bug("internal error 1 in CALCULATE_BSE_OSZI_SCALA: " // str(NCV13) &
                         // " " // str(NB1_INTO_TOT-BD(ISP)%VBMIN+1) // " " // str(NB3_INTO_TOT) // " " // &
                         str(K1) // " " // str(ISP) // " " // &
                         str(BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1,NB3_INTO_TOT-BD(ISP)%CBMIN+1,&
                         K1,ISP -ISP_LOW+1)), __FILE__, __LINE__)
                   ENDIF
#endif
                   !IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K3,ISP))) CYCLE
                   !IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)))  CYCLE
                   IF(NB1_INTO_TOT > BD(ISP)%VBMAX) THEN
                      WEIGHT=(1.0_q-WHF%FERTOT(NB3_INTO_TOT, K3, ISP))
                      ! - <core_m | r | c_k1,n3> =  (-) (i) < core | nabla | c_k1,n3>
                      IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
                         CDER_BETWEEN_STATE = CDER_CORECON_Q(NB1_INTO_TOT - BD(ISP)%VBMAX,NB3_INTO_TOT, K1,ISP)
                      ELSE
                         CDER_BETWEEN_STATE = CDER_CORECON(NB1_INTO_TOT - BD(ISP)%VBMAX,NB3_INTO_TOT, K1,ISP,1:3)*(HSQDTM*2.0_q)
                      ENDIF
                   ELSE
                      IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)))  CYCLE
                      WEIGHT=(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)-WHF%FERTOT(NB3_INTO_TOT, K3, ISP))
                      ! CDER_BETWEEN_STATE =  <v_k1,n1|- i d/dq_j c_k1+q,n3> =  - <v_k1,n1| r | c_k1,n3>  for q->0
                      ! or   <v_k1,n1| e-iqr  |c_k1+q,n3>  at finite q
#ifndef ACC_OFFLOAD
                      CALL CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT)
#else
                      CALL CDER_BETWEEN_STATES_ROTATED_ACC(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT,KPOINTS_FULL_ORIG)
#endif
                   ENDIF
                   R_CTRANS1=R_CTRANS1+CDER_BETWEEN_STATE(1)*AMAT(IROW+(JCOL-1)*DESC(LLD_))*WEIGHT
                   R_CTRANS2=R_CTRANS2+CDER_BETWEEN_STATE(2)*AMAT(IROW+(JCOL-1)*DESC(LLD_))*WEIGHT
                   R_CTRANS3=R_CTRANS3+CDER_BETWEEN_STATE(3)*AMAT(IROW+(JCOL-1)*DESC(LLD_))*WEIGHT
                ENDDO ! I2
             ENDDO  ! I1
             CTRANS(LAMBDA,1)=R_CTRANS1
             CTRANS(LAMBDA,2)=R_CTRANS2
             CTRANS(LAMBDA,3)=R_CTRANS3
          ENDIF
       ENDDO
    ENDIF

    ! now sum CTRANS over all nodes for the current LAMBDA
    CALLMPI( M_sum_z(WHF%WDES%COMM, CTRANS, 3*NCV))

    !$ACC EXIT DATA COPYOUT(CTRANS) DELETE(VKPT)  __IF_ASYNC__
    !$ACC EXIT DATA DELETE(KPOINTS_FULL_ORIG,KPOINTS_FULL) __IF_ASYNC__
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

    DO LAMBDA=1, NCV
       DO I=1,3
          DO J=1,3
             CTRANS_SQUARE(I,J)=CTRANS(LAMBDA,I)*CONJG(CTRANS(LAMBDA,J))*CH_AMPLIFICATION
          ENDDO
       ENDDO

       R_AND_INTENSITY(1,LAMBDA)= R(LAMBDA)
       R_AND_INTENSITY(2,LAMBDA)= (REAL(CTRANS_SQUARE(1,1),q)+REAL(CTRANS_SQUARE(2,2),q)+REAL(CTRANS_SQUARE(3,3),q))

       ! add to reducible polarizability X(w)
       IF (ANTIRES<0) THEN
          ! only antiresonant part (helpful for debugging)
          DO NOMEGA=1,NOMEGA_DIM
             EPS(NOMEGA,:,:)=EPS(NOMEGA,:,:)+RSPIN*WHF%WDES%WTKPT(1)*CTRANS_SQUARE* & 
                  (1/((NOMEGA-1)*DOMEGA-R(LAMBDA)+CMPLX(0,SHIFT,q)))
          ENDDO
       ELSE
          DO NOMEGA=1,NOMEGA_DIM
             EPS(NOMEGA,:,:)=EPS(NOMEGA,:,:)+RSPIN*WHF%WDES%WTKPT(1)*CTRANS_SQUARE* & 
                  (1/( (NOMEGA-1)*DOMEGA-R(LAMBDA)+CMPLX(0,SHIFT,q)) &
                  +1/(-(NOMEGA-1)*DOMEGA-R(LAMBDA)-CMPLX(0,SHIFT,q)))
          ENDDO
       ENDIF
    ENDDO

    CALL XML_VECARRAY("opticaltransitions")
    CALL XML_ARRAY_REAL(R_AND_INTENSITY,"(F10.4,' ')")
    CALL XML_CLOSE_TAG("varray")

    ! multiply by Coloumb kernel 4 pi e^2 / q^2 (in a.u.)
    ! to obtain v X^red
    ! for q->0 the 1/q^2 is accounted for in the CDER_BETWEEN_STATES at q=0
    EPS=EPS*EDEPS/LATT_CUR%OMEGA

    IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
       ! divide by q^2 at finite q
       DKX=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(1,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(1,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(1,3)
       DKY=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(2,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(2,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(2,3)
       DKZ=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(3,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(3,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(3,3)

       EPS=EPS/((DKX**2+DKY**2+DKZ**2)*TPI**2)
       ! eps_mac = ( 1 + v X^red(w))^-1
       ! at q/=0 eps is a scalar
       ! all components of EPS here are identical
       DO NOMEGA=1,NOMEGA_DIM
          EPS(NOMEGA,1,1)=1/(1+EPS(NOMEGA,1,1))
          EPS(NOMEGA,2,2)=1/(1+EPS(NOMEGA,2,2))
          EPS(NOMEGA,3,3)=1/(1+EPS(NOMEGA,3,3))
       ENDDO
    ELSE
       ! at q=0, we have calculated EPS= v X'^red, where X'^red is the
       ! the reducible polarizability with the amputed Coulomb kernel
       ! (component G->0 removed)
       ! epsilon_mac = 1 - v X'^red (Onida et al. Rev. Mod. Phys. 74, 601)
       ! so, first create " - v X'^red" (do the multiplication with -1) ...
       EPS(1:NOMEGA_DIM,:,:)=-1*EPS(1:NOMEGA_DIM,:,:)
       ! ... and since we have already "absorbed" a minus sign above we have
       ! now to take "1 + [ - v X'^red ]" on the diagonal
       DO NOMEGA=1,NOMEGA_DIM
          EPS(NOMEGA,1,1)=1+EPS(NOMEGA,1,1)
          EPS(NOMEGA,2,2)=1+EPS(NOMEGA,2,2)
          EPS(NOMEGA,3,3)=1+EPS(NOMEGA,3,3)
       ENDDO
    ENDIF

    IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
       CALL XML_EPSILON_SCALAR_W(DOMEGA, REAL(EPS(:,1,1),q), AIMAG(EPS(:,1,1)), NOMEGA_DIM)
    ELSE
       CALL XML_EPSILON_W(DOMEGA, REAL(EPS,q), AIMAG(EPS), NOMEGA_DIM )
    ENDIF

#ifdef VASP_HDF5
    IF (LGWLF) THEN
       PREFIX = "bse"
    ELSE
       PREFIX = "tdhf"
    END IF
    CALL VH5_WRITE_OPTICAL_TRANSITIONS(IH5OUTFILEID, R_AND_INTENSITY)
    CALL VH5_WRITE_DIELECTRIC_ENERGIES(IH5OUTFILEID, DOMEGA, NOMEGA_DIM, PREFIX)

    IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
       CALL VH5_WRITE_DIELECTRIC_SCALAR_DYNAMIC(IH5OUTFILEID, EPS(1:NOMEGA_DIM,1,1), PREFIX)
    ELSE
       CALL VH5_WRITE_DIELECTRIC_DYNAMIC(IH5OUTFILEID, EPS, PREFIX)
    ENDIF
#endif

    PROFILING_STOP('calculate_bse_oszi_scala')

  END SUBROUTINE CALCULATE_BSE_OSZI_SCALA


!****************** SUBROUTINE INIT_BSE_OSZI_STRENGTH *****************
!
!> Initializes the BSE oscillator strength
!> this is done quite analogously to the calculation of the
!> determination of the bse::CALCULATE_BSE_OSZI_STRENGTH routine
!
!**********************************************************************

  SUBROUTINE INIT_BSE_OSZI_STRENGTH(WHF, LATT_CUR, ISP_LOW, ISP_HIGH, BD, BSE_INDEX, NSTRIPV, NSTRIPC, NQPOINT, OPTMAT)
    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    USE core_con_mat, ONLY: CH_AMPLIFICATION
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    TYPE (latt) LATT_CUR
    INTEGER :: ISP_LOW, ISP_HIGH
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
#ifdef single_prec_bse
    GDEFS, INTENT(INOUT) :: OPTMAT(:,:)
#endif

#ifdef double_prec_bse
    GDEF, INTENT(INOUT) :: OPTMAT(:,:)
#endif
    ! local
    INTEGER ISP
    INTEGER :: NCV13
    INTEGER K1, K3
    INTEGER N1, N3, NQPOINT
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    REAL(q) :: WEIGHT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG
    GDEF       :: CDER_BETWEEN_STATE(3)
    INTEGER :: COLUMN, NSTRIPV, NSTRIPC, NPOS1, NPOS3
    INTEGER :: NN1, NN2
    LOGICAL :: LINV
    REAL(q)     :: VKPT(3)
    !$ACC ROUTINE(FILLED_XI_ORBITAL,EMPTY_XI_ORBITAL,CDER_BETWEEN_STATES_ROTATED) SEQ
    !$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX,KPOINT_IN_FULL_GRID) SEQ

    !$ACC UPDATE DEVICE(CDER_BETWEEN_STATES,CDER_CORECON)  __IF_ASYNC__
    !$ACC ENTER DATA CREATE(CDER_BETWEEN_STATE) COPYIN(KPOINTS_FULL_ORIG) __IF_ASYNC__
    !$ACC PARALLEL LOOP PRESENT(OPTMAT,LATT_CUR,CDER_BETWEEN_STATES) &
    !$ACC PRIVATE(CDER_BETWEEN_STATE,NK1_IN_KPOINTS_FULL_ORIG,WEIGHT,K1,ISP,NPOS1,NPOS3,NB3_INTO_TOT,NB1_INTO_TOT,NCV13) __IF_ASYNC__
    DO COLUMN=1,BSE_INDEX%COLUMNS
       K1=BSE_INDEX%CK1(COLUMN)
       ISP=BSE_INDEX%CISP1(COLUMN)
       NPOS1=BSE_INDEX%CNPOS1(COLUMN)
       NPOS3=BSE_INDEX%CNPOS3(COLUMN)
       VKPT=WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT)
       K3 =KPOINT_IN_FULL_GRID(VKPT,KPOINTS_FULL)

       DO N3=NPOS3,NPOS3+NSTRIPC-1
          NB3_INTO_TOT=N3
          IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K3,ISP))) CYCLE
          DO N1=NPOS1,NPOS1+NSTRIPV-1
             NB1_INTO_TOT=N1

             IF (BSE_INDEX%MASK(1+N1-NPOS1+(N3-NPOS3)*NSTRIPV,COLUMN) == 1) THEN
                NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)

                IF(NB1_INTO_TOT > BD(ISP)%VBMAX) THEN
                   WEIGHT=(1.0_q-WHF%FERTOT(NB3_INTO_TOT, K1, ISP))
                   IF (ANY(ABS(WHF%WDES%VKPT(:,NQPOINT))>1E-6_q)) THEN
                      CDER_BETWEEN_STATE = CDER_CORECON_Q(NB1_INTO_TOT - BD(ISP)%VBMAX,NB3_INTO_TOT, K1,ISP)/(HSQDTM*2.0_q)
                   ELSE
                      CDER_BETWEEN_STATE = CDER_CORECON(NB1_INTO_TOT - BD(ISP)%VBMAX,NB3_INTO_TOT, K1,ISP,1:3)
                   ENDIF
                   IF (IBSE==1) CDER_BETWEEN_STATE=CDER_BETWEEN_STATE*AUTOA
                   IF (IBSE==3) CDER_BETWEEN_STATE=CDER_BETWEEN_STATE*(2._q*HSQDTM)
                   NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)
                   OPTMAT(NCV13,:)=CDER_BETWEEN_STATE*WEIGHT*SQRT(CH_AMPLIFICATION)
                ELSE
                   IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE
                   WEIGHT=(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)-WHF%FERTOT(NB3_INTO_TOT, K3, ISP))
#ifndef ACC_OFFLOAD
                   CALL CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT)
#else
                   CALL CDER_BETWEEN_STATES_ROTATED_ACC(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT,KPOINTS_FULL_ORIG)
#endif
                   NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)
                   IF (IBSE==1) CDER_BETWEEN_STATE=CDER_BETWEEN_STATE*AUTOA/(2._q*HSQDTM)
                   OPTMAT(NCV13,:)=CDER_BETWEEN_STATE*WEIGHT*SQRT(CH_AMPLIFICATION)
                ENDIF
             ENDIF

          ENDDO
       ENDDO
    ENDDO
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    !$ACC EXIT DATA DELETE(CDER_BETWEEN_STATE,KPOINTS_FULL_ORIG) __IF_ASYNC__

  END SUBROUTINE INIT_BSE_OSZI_STRENGTH

!**********************************************************************
!
!> Determines which band window is condsidered in the BSE.
!> In fact the window is usually different for up and
!> down spin making things really difficult.
!
!**********************************************************************
  SUBROUTINE SET_BAND_PARAMETERS(W, ISP, NBANDSO, NBANDSV, &
       BD, NSTRIPV, NSTRIPC)
    IMPLICIT NONE
    TYPE (wavespin) W
    INTEGER :: ISP, NBANDSO, NBANDSV
    INTEGER,OPTIONAL :: NSTRIPV, NSTRIPC
    TYPE(banddesc) :: BD
  ! local
    INTEGER :: N, NK, NCORE

    INTEGER LAST_FILLED, FIRST_EMPTY

    PROFILING_START("set_band_parameters")
    NCORE = 0

    IF (ICORELEVEL/=0) NCORE = NITYP_CL * (2*L_CL+1)

    LAST_FILLED=0
    FIRST_EMPTY=W%WDES%NB_TOT
    DO NK=1,W%WDES%NKPTS
       IF (W%WDES%WTKPT(NK)==0) CYCLE
       LAST_FILLED=MAX(LAST_FILLED, LAST_FILLED_XI_NOMOD(W,NK,ISP))
       FIRST_EMPTY=MIN(FIRST_EMPTY, FIRST_EMPTY_XI_NOMOD(W,NK,ISP))
    ENDDO

    ! determine VBMIN and VBMAX
    BD%VBMAX=LAST_FILLED
    BD%VBMAXpCORE = BD%VBMAX + NCORE
    BD%VBMIN=MAX(LAST_FILLED-NBANDSO+1,1)
    IF (BD%VBMIN > BD%VBMAXpCORE ) BD%VBMIN=BD%VBMAXpCORE

    BD%CBMIN=FIRST_EMPTY
    BD%CBMAX=MIN(FIRST_EMPTY+NBANDSV-1,W%WDES%NB_TOT)


    ! aT
    ! if we have one band VBMIN == VBMAX, but if we have no valence bands VBMIN == VBMAX as well,
    ! so we need to treat NBANDS == 0 in a special way
    IF(NBANDSO<1) THEN
      NSTRIPV=MAX(NSTRIPV,NCORE)
    ELSE
      NSTRIPV=MAX(NSTRIPV,BD%VBMAX-BD%VBMIN+1+NCORE)
    ENDIF
    NSTRIPC=MAX(NSTRIPC,BD%CBMAX-BD%CBMIN+1)

    PROFILING_STOP("set_band_parameters")

  END SUBROUTINE SET_BAND_PARAMETERS

!**********************************************************************
!
!> Determines the index arrays that allow
!> to determine the indices N1(valence),N2(conduction),NK,ISP
!> from the global index into the BSE matrix.
!> Also an index array is set up to determine the global index
!> from the band indices, k-point and spin.
!
!**********************************************************************

  SUBROUTINE SET_BSE_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, BSE_INDEX, W, OMEGAMAX_IN, NCV, NQPOINT, NSTRIPV, NSTRIPC, LTRIANGLE)

    INTEGER, INTENT(IN) :: ISP_LOW, ISP_HIGH, ALLOC_KPOINTS
    TYPE (banddesc), INTENT(IN) :: BD(:)
    TYPE (bse_matrix_index), INTENT(OUT) :: BSE_INDEX
    TYPE (wavespin), INTENT(IN) :: W
    REAL(q) :: OMEGAMAX_IN    !< maximum energy difference between conduction and valence band to be considered in BSE
    INTEGER :: NCV            !< number of valence conduction band pairs
    INTEGER :: NQPOINT        !< q point at which exciton is calculated
    ! local
    INTEGER :: ALLOC_VALENCE_BANDS, ALLOC_CONDUCTION_BANDS, ALLOC_SPIN
    INTEGER :: ALLOC_VALENCE_BLOCKS, ALLOC_CONDUCTION_BLOCKS
    INTEGER :: N1, N3, K1, K3
    INTEGER :: ISP, INDEX, BNCV
    REAL(q) :: OMEGAMAX_ACT   !< actual maximum energy difference between conduction and valence band
    INTEGER :: NSTRIPV, NSTRIPC
    INTEGER :: NSTRIP1, NSTRIP3

    INTEGER :: BLOCK_INDEX, COLUMNS
    INTEGER :: COL, ROW, LCOL, LROW, ROWS
    INTEGER :: NPOS1,NPOS2,K2_LOW,ISP2_LOW
    INTEGER :: NPOS3,NPOS4
    INTEGER :: ISP1, ISP2,K2
    LOGICAL :: LTRIANGLE
    INTEGER, ALLOCATABLE :: MASK_BLOCK(:)
    INTEGER :: ALLBLOCKS, SUBBLOCKS

    REAL(q) :: OMEGAMAX

    PROFILING_START('set_bse_matrix_index')

    IF (OMEGAMAX_IN>0) THEN
       OMEGAMAX=OMEGAMAX_IN
    ELSE
       OMEGAMAX=1E10
    ENDIF

    OMEGAMAX_ACT=0

    ALLOC_SPIN=ISP_HIGH-ISP_LOW+1
    ALLOC_VALENCE_BANDS=0
    ALLOC_CONDUCTION_BANDS=0

    DO ISP=ISP_LOW,ISP_HIGH
       ALLOC_VALENCE_BANDS   =MAX(ALLOC_VALENCE_BANDS   ,BD(ISP)%VBMAXpCORE-BD(ISP)%VBMIN+1)
       ALLOC_CONDUCTION_BANDS=MAX(ALLOC_CONDUCTION_BANDS,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1)
    ENDDO
    ALLOC_VALENCE_BLOCKS=CEILING(REAL(ALLOC_VALENCE_BANDS,KIND=q)/REAL(NSTRIPV,KIND=q))
    ALLOC_CONDUCTION_BLOCKS=CEILING(REAL(ALLOC_CONDUCTION_BANDS,KIND=q)/REAL(NSTRIPC,KIND=q))

    ALLOCATE(BSE_INDEX%INDEX(ALLOC_VALENCE_BANDS,ALLOC_CONDUCTION_BANDS,ALLOC_KPOINTS,ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%N1(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%N3(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%NK(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%ISP(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))

    ALLOCATE(BSE_INDEX%CK1(ALLOC_KPOINTS*ALLOC_SPIN*ALLOC_VALENCE_BLOCKS*ALLOC_CONDUCTION_BLOCKS))
    ALLOCATE(BSE_INDEX%CISP1(ALLOC_KPOINTS*ALLOC_SPIN*ALLOC_VALENCE_BLOCKS*ALLOC_CONDUCTION_BLOCKS))
    ALLOCATE(BSE_INDEX%CNPOS1(ALLOC_KPOINTS*ALLOC_SPIN*ALLOC_VALENCE_BLOCKS*ALLOC_CONDUCTION_BLOCKS))
    ALLOCATE(BSE_INDEX%CNPOS3(ALLOC_KPOINTS*ALLOC_SPIN*ALLOC_VALENCE_BLOCKS*ALLOC_CONDUCTION_BLOCKS))

    !
    ALLOCATE(BSE_INDEX%MASK(NSTRIPV*NSTRIPC,ALLOC_KPOINTS*ALLOC_SPIN*ALLOC_VALENCE_BLOCKS*ALLOC_CONDUCTION_BLOCKS))
    ALLOCATE(MASK_BLOCK(NSTRIPV*NSTRIPC))
    MASK_BLOCK=0

    BSE_INDEX%COLUMNS=0

    BSE_INDEX%INDEX=1
    BSE_INDEX%N1=0
    BSE_INDEX%N3=0
    BSE_INDEX%NK=0
    BSE_INDEX%ISP=0
    BSE_INDEX%MASK=0

! first spin component always U
    INDEX=0
    COLUMNS=0
    DO ISP=ISP_LOW,ISP_HIGH
       ! K1 is the k-point index for the valence band

       DO K1=1, ALLOC_KPOINTS

          ! if the first few conduction bands in the block are occupied
          ! we set INDEX to the next value, so that we if we take the index at NPOS3
          ! we get the index of the next element
          BSE_INDEX%INDEX(:,:,K1,ISP-ISP_LOW+1)=INDEX+1
          DO NPOS3=BD(ISP)%CBMIN, BD(ISP)%CBMAX, NSTRIPC
             NSTRIP3=MIN(BD(ISP)%CBMAX+1-NPOS3,NSTRIPC)

             DO NPOS1=BD(ISP)%VBMIN, BD(ISP)%VBMAXpCORE, NSTRIPV
                NSTRIP1=MIN(BD(ISP)%VBMAXpCORE+1-NPOS1,NSTRIPV)

                ! K3 is the k-point index for the conduction band state (unoccupied)
                K3 =KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
                MASK_BLOCK = 0
                BNCV = 0
                DO N3=NPOS3, NPOS3+NSTRIP3-1
                   DO N1=NPOS1, NPOS1+NSTRIP1-1
                      ! if occupancy of the valence band is smaller than 0.5 we skip the orbital as well
                      !IF (EMPTY_XI_ORBITAL(W%FERTOT(N1,K1,ISP))) CYCLE
                      ! another safeguard: if the transition from k1 to k3 has energy difference smaller or eq to 0.0 then cycle
!                      IF (REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q) <= 0.0_q) CYCLE

                      ! "+1" because we start INDEX from zero
                      BSE_INDEX%INDEX(N1-BD(ISP)%VBMIN+1:,N3-BD(ISP)%CBMIN+1:,K1,ISP-ISP_LOW+1)=INDEX+1
                      BSE_INDEX%N1(INDEX+1) =  N1
                      BSE_INDEX%N3(INDEX+1) =  N3
                      BSE_INDEX%NK(INDEX+1) =  K1
                      BSE_INDEX%ISP(INDEX+1)=  ISP

                      IF(N1 <= BD(ISP)%VBMAX) THEN
                         OMEGAMAX_ACT=MAX(OMEGAMAX_ACT, REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q))
                         ! if occupancy of the valence band is smaller than 0.5 we skip the orbital as well
                         IF (EMPTY_XI_ORBITAL(W%FERTOT(N1,K1,ISP))) CYCLE
                         ! another safeguard: if the transition from k1 to k3 has energy difference smaller or eq to 0.0 then cycle
                         IF (REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q) <= 0.0_q) CYCLE
                      ELSE
                         OMEGAMAX_ACT=MAX(OMEGAMAX_ACT, REAL(W%CELTOT(N3,K3,ISP)-COREEIGENVALUE,q))
                         IF (REAL(W%CELTOT(N3,K3,ISP)-COREEIGENVALUE,q) <= 0.0_q) CYCLE
                      ENDIF

                      ! safe guard, it might be that the number of bands (NB_TOT) exceeds the
                      ! number of plane wave coefficients stored in NB_TOTK(NK,ISP)
                      ! skip those useless bands
                      IF (N3> W%WDES%NB_TOTK(K3,ISP)) CYCLE

                      ! if occupancy is larger than 0.5 the band is occupied and we do not treat it as unoccupied state
                      IF (FILLED_XI_ORBITAL(W%FERTOT(N3,K3,ISP))) CYCLE

                      IF (REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q) < OMEGAMAX ) THEN
                         INDEX=INDEX+1
                         BNCV = BNCV + 1
                         MASK_BLOCK(1+N1-NPOS1+(N3-NPOS3)*NSTRIPV) = 1
                      ENDIF

                   ENDDO ! N1
                ENDDO ! N3

                IF(BNCV > 0) THEN
                   COLUMNS = COLUMNS + 1
                   BSE_INDEX%CK1(COLUMNS) =  K1
                   BSE_INDEX%CISP1(COLUMNS)= ISP
                   BSE_INDEX%CNPOS1(COLUMNS) =  NPOS1
                   BSE_INDEX%CNPOS3(COLUMNS) =  NPOS3
                   BSE_INDEX%MASK(:,COLUMNS) = MASK_BLOCK(:)
                ENDIF
            ENDDO ! NPOS1
         ENDDO ! NPOS3
       ENDDO ! K1
    ENDDO ! ISP
    ! number of pairs equals final index
    NCV=INDEX
    BSE_INDEX%NCV=NCV
    BSE_INDEX%COLUMNS = COLUMNS


    !BLOCK_INDEX = (ALLOC_SPIN*ALLOC_VALENCE_BLOCKS*ALLOC_KPOINTS)**2
    BLOCK_INDEX = BSE_INDEX%COLUMNS*(BSE_INDEX%COLUMNS+1._q)/2._q
    ALLOCATE(BSE_INDEX%BK1(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BK2(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BISP1(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BISP2(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BNPOS1(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BNPOS2(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BNPOS3(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BNPOS4(BLOCK_INDEX))
    ALLOCATE(BSE_INDEX%BMAP(BLOCK_INDEX))

    BSE_INDEX%BLOCKS=0
    BSE_INDEX%BK1=0
    BSE_INDEX%BK2=0
    BSE_INDEX%BISP1=0
    BSE_INDEX%BISP2=0
    BSE_INDEX%BNPOS1=0
    BSE_INDEX%BNPOS2=0
    BSE_INDEX%BNPOS3=0
    BSE_INDEX%BNPOS4=0


    BLOCK_INDEX = 0
    IF (LTRIANGLE) THEN
      ROWS = CEILING(REAL(BSE_INDEX%COLUMNS+1,KIND=q)/2._q)
    ELSE
      ROWS = BSE_INDEX%COLUMNS
    ENDIF
    DO ROW=1,ROWS

       DO COL=1,BSE_INDEX%COLUMNS

         IF (COL >= ROW .OR. .NOT. LTRIANGLE) THEN

            ! in even matrix the last row is in lower triangle
            IF (ROW >= CEILING(REAL(BSE_INDEX%COLUMNS,KIND=q)/2._q+1) .AND. LTRIANGLE) CYCLE

            BLOCK_INDEX = BLOCK_INDEX + 1

            K1 = BSE_INDEX%CK1(ROW)
            ISP1 = BSE_INDEX%CISP1(ROW)
            NPOS1 = BSE_INDEX%CNPOS1(ROW)
            NPOS3 = BSE_INDEX%CNPOS3(ROW)

            K2 = BSE_INDEX%CK1(COL)
            ISP2 = BSE_INDEX%CISP1(COL)
            NPOS2 = BSE_INDEX%CNPOS1(COL)
            NPOS4 = BSE_INDEX%CNPOS3(COL)

            BSE_INDEX%BISP1(BLOCK_INDEX)=ISP1
            BSE_INDEX%BISP2(BLOCK_INDEX)=ISP2

            BSE_INDEX%BK1(BLOCK_INDEX)=K1
            BSE_INDEX%BK2(BLOCK_INDEX)=K2

            BSE_INDEX%BNPOS1(BLOCK_INDEX)=NPOS1
            BSE_INDEX%BNPOS2(BLOCK_INDEX)=NPOS2
            BSE_INDEX%BNPOS3(BLOCK_INDEX)=NPOS3
            BSE_INDEX%BNPOS4(BLOCK_INDEX)=NPOS4
            BSE_INDEX%BMAP(BLOCK_INDEX)=BLOCK_INDEX

            ! BLOCK_INDEX index in triangular compact matrix 
            ! BSE_INDEX%BMAP(BLOCK_INDEX) index in full ractangular matrix
            BSE_INDEX%BMAP(BLOCK_INDEX) = BLOCK_INDEX + (ROW-1)*(ROW)/2

         ELSE
            BLOCK_INDEX = BLOCK_INDEX + 1

            LCOL = BSE_INDEX%COLUMNS - COL + 1
            ! first element of lower half is in second row, so +2
            LROW = BSE_INDEX%COLUMNS - ROW + 2

            K1 = BSE_INDEX%CK1(LROW)
            ISP1 = BSE_INDEX%CISP1(LROW)
            NPOS1 = BSE_INDEX%CNPOS1(LROW)
            NPOS3 = BSE_INDEX%CNPOS3(LROW)

            K2 = BSE_INDEX%CK1(LCOL)
            ISP2 = BSE_INDEX%CISP1(LCOL)
            NPOS2 = BSE_INDEX%CNPOS1(LCOL)
            NPOS4 = BSE_INDEX%CNPOS3(LCOL)

            BSE_INDEX%BISP1(BLOCK_INDEX)=ISP1
            BSE_INDEX%BISP2(BLOCK_INDEX)=ISP2

            BSE_INDEX%BK1(BLOCK_INDEX)=K1
            BSE_INDEX%BK2(BLOCK_INDEX)=K2

            BSE_INDEX%BNPOS1(BLOCK_INDEX)=NPOS1
            BSE_INDEX%BNPOS2(BLOCK_INDEX)=NPOS2
            BSE_INDEX%BNPOS3(BLOCK_INDEX)=NPOS3
            BSE_INDEX%BNPOS4(BLOCK_INDEX)=NPOS4

            ! BLOCK_INDEX index in triangular compact matrix 
            ! BSE_INDEX%BMAP(BLOCK_INDEX) index in full ractangular matrix
            BSE_INDEX%BMAP(BLOCK_INDEX) = (BSE_INDEX%COLUMNS-ROW+1)*BSE_INDEX%COLUMNS + BSE_INDEX%COLUMNS-COL+1
         ENDIF

       ENDDO

    ENDDO
    DEALLOCATE(MASK_BLOCK)

    BSE_INDEX%BLOCKS = BLOCK_INDEX

   ! ! Print matrix distribution
   ! IF (W%WDES%COMM%NODE_ME == 1 ) THEN
   !    DO BLOCK_INDEX=1,BSE_INDEX%BLOCKS
   !      K1=BSE_INDEX%BK1(BLOCK_INDEX)
   !      K2=BSE_INDEX%BK2(BLOCK_INDEX)
   !      WRITE(*,"(A,I2,A,I2,A)" ,advance="no") "(", K1,".",K2,")"
   !      IF (MOD(BLOCK_INDEX,ALLOC_KPOINTS) == 0 ) WRITE(*,*) ''
   !    ENDDO
   !    WRITE(*,*) ''
   ! ENDIF

   ! IF (W%WDES%COMM%NODE_ME == 1 ) THEN
   !    DO BLOCK_INDEX=1,BSE_INDEX%BLOCKS
   !      NPOS1=BSE_INDEX%BNPOS1(BLOCK_INDEX)
   !      NPOS2=BSE_INDEX%BNPOS2(BLOCK_INDEX)
   !      NPOS3=BSE_INDEX%BNPOS3(BLOCK_INDEX)
   !      NPOS4=BSE_INDEX%BNPOS4(BLOCK_INDEX)
   !      !WRITE(*,"(A,I2,A,I2,A)" ,advance="no") "(", K1,".",K2,")"
   !      WRITE(*,"(A,I2,A,I2,A,I2,A,I2,A,I2,A)" ,advance="no") "(", NPOS1,":",NPOS2,",", NPOS3,":",NPOS4,")"
   !      IF (MOD(BLOCK_INDEX,ALLOC_KPOINTS*ALLOC_VALENCE_BLOCKS*ALLOC_CONDUCTION_BLOCKS) == 0 ) WRITE(*,*) ''
   !    ENDDO
   !    WRITE(*,*) ''
   ! ENDIF

    IF (OMEGAMAX_IN<=0) THEN
       ! be shure that the highest transition is never thrown away
       ! add a small safety measure
       OMEGAMAX_IN=OMEGAMAX_ACT*(1._q+1E-5_q)
    ENDIF

  PROFILING_STOP('set_bse_matrix_index')

  END SUBROUTINE SET_BSE_MATRIX_INDEX


!**********************************************************************
!
!> Deallocates the BSE_INDEX array
!
!**********************************************************************

  SUBROUTINE DEALLOCATE_BSE_MATRIX_INDEX( BSE_INDEX )
    TYPE (bse_matrix_index) :: BSE_INDEX

    DEALLOCATE(BSE_INDEX%INDEX)
    DEALLOCATE(BSE_INDEX%N1)
    DEALLOCATE(BSE_INDEX%N3)
    DEALLOCATE(BSE_INDEX%NK)
    DEALLOCATE(BSE_INDEX%ISP)

    DEALLOCATE(BSE_INDEX%CK1)
    DEALLOCATE(BSE_INDEX%CISP1)
    DEALLOCATE(BSE_INDEX%CNPOS1)
    DEALLOCATE(BSE_INDEX%CNPOS3)
    DEALLOCATE(BSE_INDEX%BK1)
    DEALLOCATE(BSE_INDEX%BK2)
    DEALLOCATE(BSE_INDEX%BISP1)
    DEALLOCATE(BSE_INDEX%BISP2)
    DEALLOCATE(BSE_INDEX%BNPOS1)
    DEALLOCATE(BSE_INDEX%BNPOS2)
    DEALLOCATE(BSE_INDEX%BNPOS3)
    DEALLOCATE(BSE_INDEX%BNPOS4)
    DEALLOCATE(BSE_INDEX%MASK)
    DEALLOCATE(BSE_INDEX%BMAP)

  END SUBROUTINE DEALLOCATE_BSE_MATRIX_INDEX

!**********************************************************************
!
!> Diagonalizes BSE matrix using one of the many possible calls
!
!**********************************************************************

  SUBROUTINE DIAG_BSE_MATRIX( NCV, AMAT_SCALA, R, LscaLAPACKaware, MY_COMM, IU0, DESC )

  GPACC USE moffload

    USE wave_high
    USE scala
    IMPLICIT NONE
    INTEGER, INTENT(IN OUT):: NCV           !< number of pairs (dimension of problem)
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT),TARGET, ALLOCATABLE :: AMAT_SCALA(:)  !< in: matrix distributed according to DESC
#endif

#ifdef double_prec_bse
    GDEF, INTENT(IN OUT),TARGET, ALLOCATABLE :: AMAT_SCALA(:) !< in: matrix distributed according to DESC
#endif
    REAL(q), INTENT(OUT)::  R(:)            !< out: eigenvalues of matrix
    LOGICAL, INTENT(IN) :: LscaLAPACKaware  !< distribute data using BLACS
    TYPE(communic), INTENT(IN) :: MY_COMM   !< communicator to be passed to the subroutine that computes the trace
    INTEGER, INTENT(IN) :: IU0              !< VASP unit for debug messages
    INTEGER, INTENT(IN) :: DESC(DLEN_)      !< SCALAPACK matrix descriptor
! local
    REAL(q) ::  ABSTOL
    REAL(q) :: VL, VU

    INTEGER, PARAMETER :: LWORK=32
#ifdef single_prec_bse
    REAL(qs), ALLOCATABLE ::  RWORK(:)
    GDEFS, ALLOCATABLE    ::  CWRK(:)
    REAL(qs) :: R_TMP(SIZE(R))
#endif

#ifdef double_prec_bse
    REAL(q), ALLOCATABLE ::  RWORK(:)
    GDEF, ALLOCATABLE    ::  CWRK(:)
    REAL(q) :: R_TMP(SIZE(R))
#endif
    INTEGER, ALLOCATABLE :: IWORK(:), INFOZ(:)
    INTEGER :: IL, IU, NB_CALC
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: COVL(:,:)
    GDEFS, POINTER :: AMAT_SCALA_P(:,:)=>NULL()
#endif
#ifdef double_prec_bse
    GDEF, ALLOCATABLE :: COVL(:,:)
    GDEF, POINTER :: AMAT_SCALA_P(:,:)=>NULL()
#endif
    INTEGER :: IFAIL, NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER, EXTERNAL    ::     NUMROC
!jF:  this is needed for handling potentially "giant arrays" ...
    INTEGER(qi8) :: SIZE8,SIZE8A,SIZE8B
! some [sca]LAPACK routines require parameter "ABSTOL" which controls the
! final accuracy; it is most important for single precision routines, for
! double precision routines just a "small enough number" is sufficient
! (e.g. 1.E-8). Now the "optimum parameter" is enforced for the case
! "single_prec_bse" while for "double_prec_bse" 1.E-8 would be used (the
! default) but alternatively also the optimal value for double precision
! can be enforced for "double_prec_bse" if also "high_prec_abstol" is set
!#define high_prec_abstol
#ifdef single_prec_bse
    REAL(qs) :: SLAMCH
    REAL(qs) :: PSLAMCH
#endif
#if defined (double_prec_bse) && defined (high_prec_abstol)
    REAL(q) :: DLAMCH
    REAL(q) :: PDLAMCH
#endif

    PROFILING_START('diag_bse_matrix')

! the default for the double precision case ...
    ABSTOL=1.E-8_q
! this is the recommended value for maximum accuracy in the single precision
! case (always enforced for maximum safety/accuracy for "single_prec_bse"!)
#if defined (single_prec_bse)
    ! if LDISTRIBUTED_MATRIX DESC(2) /= -1 only for master rank
    IF (IU0>0 .OR. LscaLAPACKaware) ABSTOL=2._q*PSLAMCH(DESC(CTXT_),'S')
#endif
!     this is then the corresponding recommended value for double precision
!     (only set if "high_prec_abstol" is defined, otherwise 1.E-8 is used)
#if defined (double_prec_bse) && defined (high_prec_abstol)
    ! if LDISTRIBUTED_MATRIX DESC(2) /= -1 only for master rank
    IF (IU0>0 .OR. LscaLAPACKaware) ABSTOL=2._q*PDLAMCH(DESC(CTXT_),'S')
#endif
!==========================================================================
! version using scaLAPACK
! to solve non-distributed matrix BSE_DESC must be define such that 
! block size is NCV
!==========================================================================

  IF (LscaLAPACKaware) THEN
       !$ACC ENTER DATA COPYIN(R_TMP) CREATE(COVL) __IF_ASYNC__
#ifndef USE_ZHEEVD
#if defined(ACC_OFFLOAD) & ! defined(CUSOLVERMP)
         !$ACC UPDATE SELF(AMAT_SCALA) __IF_ASYNC__
         !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif
         IF (IU0>0) WRITE(IU0,*) "BSE diagonalizing matrix PDSYEVX_ZHEEVX_DESC"
         CALL PGGEVX(AMAT_SCALA, R_TMP, NCV, DESC, MY_COMM, ABSTOL_= ABSTOL)
#if defined(ACC_OFFLOAD) & ! defined(CUSOLVERMP)
         CALL OFFLOAD_ON_POP
         !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         !$ACC UPDATE DEVICE(AMAT_SCALA,R_TMP) __IF_ASYNC__
#endif
#else
#if defined(ACC_OFFLOAD) & ! defined(CUSOLVERMP)
         !$ACC UPDATE SELF(AMAT_SCALA) __IF_ASYNC__
         !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif
         IF (IU0>0) WRITE(IU0,*) "BSE diagonalizing matrix PDSYEV_ZHEEVD_DESC"
         CALL PGGEVX(AMAT_SCALA, R_TMP, NCV, DESC, MY_COMM)
#if defined(ACC_OFFLOAD) & ! defined(CUSOLVERMP)
         CALL OFFLOAD_ON_POP
         !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         !$ACC UPDATE DEVICE(AMAT_SCALA,R_TMP) __IF_ASYNC__
#endif
#endif

       !$ACC KERNELS PRESENT(R,R_TMP) __IF_ASYNC__
       R=R_TMP
       !$ACC END KERNELS

       !$ACC EXIT DATA DELETE(R_TMP) __IF_ASYNC__

  ELSE ! LscaLAPACKaware


    ! diagonalize the matrix on master rank
    IF(IU0>=0) THEN
!==========================================================================
! version using LAPACK ZHEEVX
!==========================================================================
       VL=0 ; VU=0 ; IL=0 ; IU=0
       IF (IU0>0) WRITE(IU0,*) "BSE diagonalizing matrix ZHEEVX"
       ALLOCATE(COVL(NCV,NCV), CWRK(NCV*LWORK), RWORK(7*NCV), IWORK(5*NCV), INFOZ(NCV))

       call c_f_pointer(c_loc(AMAT_SCALA), AMAT_SCALA_P,  [NCV, NCV])
       !$ACC ENTER DATA COPYIN(AMAT_SCALA_P,R_TMP) CREATE(COVL) __IF_ASYNC__
       CALL GGEVX( 'V', 'A', 'U', NCV, AMAT_SCALA_P , NCV, VL, VU, IL, IU, &
            ABSTOL , NB_CALC , R_TMP, COVL, NCV, CWRK, &
            LWORK*NCV, RWORK, IWORK, INFOZ, IFAIL )


       !$ACC KERNELS PRESENT(R,R_TMP) __IF_ASYNC__
       AMAT_SCALA_P = COVL
       !$ACC END KERNELS

       !$ACC EXIT DATA DELETE(COVL,AMAT_SCALA_P) __IF_ASYNC__
       DEALLOCATE(COVL, CWRK, RWORK, IWORK, INFOZ)

       !$ACC KERNELS PRESENT(R,R_TMP) __IF_ASYNC__
       R=R_TMP
       !$ACC END KERNELS

       !$ACC EXIT DATA DELETE(R_TMP) __IF_ASYNC__
    ENDIF ! master rank

  ENDIF ! LscaLAPACKaware

    PROFILING_STOP('diag_bse_matrix')
  END SUBROUTINE DIAG_BSE_MATRIX

#ifdef scaLAPACK

!**********************************************************************
!
!> ScaLAPACK-aware version of the subroutine that left and right multiplies 
!> the matrix A+B by (A-B)^(1/2)
!> compare Sander, Maggio, Kresse, PRB 045209 (2015)  Equ. (40)
!> 
!> the matrix A+B is supplied in A
!> the matrix A-B is supplied as U R U+
!>
!> the transformed  matrix is therefore given by
!> C = U R U+ A U R U+
!> to save storage four matrix matrix products are calculated
!> only one temporary matrix is used
!
!**********************************************************************

 SUBROUTINE TDA_BSE_SCA( NCV, A, U, RB , IU0, DESC)
    USE wave_high
    INTEGER, INTENT(IN) :: NCV
    REAL(q), INTENT(IN) :: RB(:)
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: A(:)
    GDEFS, INTENT(IN) :: U(:)
    GDEFS,ALLOCATABLE :: T(:)
#else
    GDEF, INTENT(IN OUT)  :: A(:)
    GDEF, INTENT(IN)  :: U(:)
    GDEF,ALLOCATABLE :: T(:)
#endif
    INTEGER :: I, J
    REAL(q) :: SR(SIZE(RB))
    INTEGER :: IU0
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES,J1RES,IROW,JCOL
    INTEGER I1,I2,J1,J2
    INTEGER :: DESC(DLEN_) ! descriptor set up for matrix distribution

    PROFILING_START('tda_bse_sca')
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    IF(MYCOL /= -1 .OR. MYROW /= -1) THEN
!jF:  just for safety if A would be a "giant array" ...
      ALLOCATE(T(SIZE(A,KIND=qi8)))

      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      !$ACC UPDATE SELF(A) WAIT IF(OFFLOAD_ON)

      !$ACC ENTER DATA CREATE(SR) __IF_ASYNC__
      !$ACC UPDATE SELF(U,RB) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      !  U+ (A+B) -> T
      CALL GPGEMMbse( trans,'N', NCV, NCV, NCV, onebse, &
           U(1), 1, 1, DESC, A(1), 1, 1, DESC, zerobse, T(1), 1, 1, DESC)

      !  [U+ (A+B)] U -> A
      CALL GPGEMMbse( 'N','N', NCV, NCV, NCV, onebse, &
           T(1), 1, 1, DESC, U(1), 1, 1, DESC, zerobse, A(1), 1, 1, DESC)


      ! multiply with sqrt of eigenvalues from left and right

      !$ACC UPDATE DEVICE(A) __IF_ASYNC__
      !SR=SQRT(RB)
      CALL SQRT_EIGENVALUE(RB ,SR, IU0 )


      JCOL=0
      !$ACC UPDATE DEVICE(A) __IF_ASYNC__
      !$ACC PARALLEL LOOP COLLAPSE(force:4) PRESENT(A,SR) __IF_ASYNC__
      DO J1=1,NQ,DESC(NB_)
         DO J2=1,DESC(NB_)
            JCOL=J1+J2-1
            IF(JCOL>NQ) CYCLE
            J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

            IROW=0
            DO I1=1,NP,DESC(MB_)
               DO I2=1,DESC(MB_)
                  IROW=I1+I2-1
                  IF(IROW>NP) CYCLE
                  I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                  ! A(I,J)=A(I,J)*SR(I)*SR(J)
                  A(IROW+(JCOL-1)*DESC(LLD_))=A(IROW+(JCOL-1)*DESC(LLD_))*SR(I)*SR(J)
               ENDDO
            ENDDO
        ENDDO
      ENDDO


      !$ACC UPDATE SELF(A) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      !  U A -> T
      CALL GPGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, &
           U(1), 1, 1, DESC, A(1), 1, 1, DESC, zerobse, T(1), 1, 1, DESC)
      !  [U A ] U+ -> A
      CALL GPGEMMbse( 'N',trans, NCV, NCV, NCV, onebse, &
           T(1), 1, 1, DESC, U(1), 1, 1, DESC, zerobse, A(1), 1, 1, DESC)


      !$ACC UPDATE DEVICE(A) __IF_ASYNC__
      !$ACC EXIT DATA DELETE(SR) __IF_ASYNC__

      DEALLOCATE(T)
    ENDIF
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

    PROFILING_STOP('tda_bse_sca')
  END SUBROUTINE TDA_BSE_SCA
#endif

!**********************************************************************
!
!> Small helper routine to calculate the sqrt(R)
!> and prompt a warning if R holds negative values.
!
!**********************************************************************

  SUBROUTINE SQRT_EIGENVALUE(R ,SR, IU0 )
    REAL(q), INTENT(IN) :: R(:)           !< eigenvalues
    REAL(q), INTENT(OUT):: SR(:)          !< sqrt of eigenvalues
    INTEGER, INTENT(IN) :: IU0            !< unit for error logging
  ! local
    INTEGER :: I, NSR

    PROFILING_START('sqrt_eigenvalue')

    NSR=0
   
    !SR =SQRT(R)
    !$ACC PARALLEL LOOP PRESENT(R,SR) REDUCTION(+:NSR) __IF_ASYNC__
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =SQRT(R(I))
       ENDIF
    ENDDO
    !$ACC WAIT IF(OFFLOAD_ON)

    IF (NSR>0 .AND. IU0>=0 ) THEN 
       WRITE(IU0,'(1X,A,I5)') 'WARNING in SQRT_EIGENVALUE: some eigenvalues in R are negative',NSR
       WRITE(IU0,'(1X,A)') ' this indicates instabilities'
    ENDIF

    PROFILING_STOP('sqrt_eigenvalue')
   END SUBROUTINE SQRT_EIGENVALUE


!**********************************************************************
!
!> Small helper routine to calculate the sqrt(R) inplace 
!> and prompt a warning if R holds negative values.
!>
!> This routine is used to determine the eigenvalues
!> of the original problem Omega_lambda from the eigenvalues of the
!> squared problem (Omega_lambda^2)
!> compare Sander, Maggio, Kresse, PRB 045209 (2015)  Equ. (42).
!
!**********************************************************************

  SUBROUTINE SQRT_EIGENVALUE_INPLACE(R , IU0 )
    REAL(q) :: R(:)
    REAL(q), ALLOCATABLE :: SR(:)
    INTEGER :: I, NSR
    INTEGER :: IU0

    PROFILING_START('sqrt_eigenvalue_inplace')
    NSR=0

    !$ACC KERNELS PRESENT(R) __IF_ASYNC__
    NSR = COUNT(R <= 0.0)
    !$ACC END KERNELS

    IF (NSR>0 .AND. IU0>=0 ) THEN
       WRITE(IU0,*) 'WARNING in SQRT_EIGENVALUE_INPLACE: some eigenvalues of the squared problem are negative',NSR
       WRITE(IU0,*) ' this indicates instabilities'
    ENDIF

    !$ACC KERNELS PRESENT(R) __IF_ASYNC__
    R = SQRT(MAX(0.0_q,R))
    !$ACC END KERNELS

    PROFILING_STOP('sqrt_eigenvalue_inplace')

   END SUBROUTINE SQRT_EIGENVALUE_INPLACE



!**********************************************************************
!
!> The scaLAPACK-aware version of the routine that determines 
!> the eigenvectors of the original problem  e.g.
!> ~~~
!>  A  B    X     w  0   X
!>             =
!>  B  A    Y     0  -w  Y
!> ~~~
!> This routine returns only
!> Z = X+Y, which is required to calculate the dielectric function
!> compare Sander, Maggio, Kresse, PRB 045209 (2015)  Eqs. (58)
!>
!>
!> it is fairly easy to see that
!> ~~~
!>  (X+Y) = (A-B)^( 1/2) Z R^(-1/2)
!>  (X-Y) = (A-B)^(-1/2) Z R^1/2
!> ~~~
!> compare Sander, Maggio, Kresse, PRB 045209 (2015)  Eqs. (44-45)
!>
!> where Z are the eigenvectors of the
!> squared problem [(A-B)^(1/2) (A+B) (A-B)^(1/2) Z = Z R]
!> and R the eigenvalues omega_lambda of the original problem
!> furthermore one can show the
!> ~~~
!>  (X+Y) (X-Y)+  = 1
!> ~~~
!> and
!> ~~~
!>  X Y+  - Y X+  = 0
!> ~~~
!> and
!> ~~~
!>  X X+ - Y Y+   = 1
!> ~~~
!>
!> input:
!> Z=z (eigenvectors of squared problem) and R
!> A-B is passed as A-B  = U RB U+
!>
!**********************************************************************


  SUBROUTINE TDA_BSE_EIGENVECTORS_SCA( NCV, U, RB, Z,  R , IU0, DESC)
    USE wave_high
    INTEGER, INTENT(IN) :: NCV            !< number of pairs (dimension of problem)
    REAL(q), INTENT(IN) :: RB(:)          !< eigenvalues of B
    REAL(q), INTENT(IN) :: R(:)           !< sqrt of eigenvalues of squared problem
    INTEGER, INTENT(IN) :: IU0            !< unit for error logging
    INTEGER, INTENT(IN) :: DESC(DLEN_)    !< SCALAPACK matrix descriptor
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: Z(:)         !< in: eigenvectors of squared problem
                                          !< out: X+Y
    GDEFS, INTENT(IN OUT) :: U(:)         !< in: eigenvectors of A-B
                                          !< out: 0
    GDEFS,ALLOCATABLE :: T(:)             !< auxilary
#else
    GDEF, INTENT(IN OUT)   :: Z(:)        !< in: eigenvectors of squared problem
                                          !< out: X+Y
    GDEF, INTENT(IN OUT)   :: U(:)        !< in: eigenvectors of A-B
                                          !< out: 0
    GDEF,ALLOCATABLE :: T(:)              !< auxilary

#endif
! local
    INTEGER :: I, J
    REAL(q) :: SRB(SIZE(R)),SR(SIZE(R))
    INTEGER :: NSR, NSRB
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES,J1RES,IROW,JCOL
    INTEGER I1,I2,J1,J2

    PROFILING_START('tda_bse_eigenvectors_sca')
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    IF(MYCOL /= -1 .OR. MYROW /= -1) THEN
!jF:  just for safety if Z would be a "giant array" ...
      ALLOCATE(T(SIZE(Z,KIND=qi8)))

      !$ACC ENTER DATA CREATE(SR,SRB,T) __IF_ASYNC__

      !  U+ Z -> T2
      !$ACC UPDATE SELF(U,Z,R,RB) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      CALL GPGEMMbse( trans,'N', NCV, NCV, NCV, onebse, &
           U(1), 1, 1, DESC, Z(1), 1, 1, DESC, zerobse, T(1), 1, 1, DESC)
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

      !$ACC UPDATE DEVICE(T) __IF_ASYNC__
      NSR=0
      !SR =1/SQRT(R)

      !$ACC PARALLEL LOOP PRESENT(R,SR) REDUCTION(+:NSR) __IF_ASYNC__
      DO I=1,SIZE(R)
         IF (R(I)<=0.0_q ) THEN
            SR(I)=0.0_q
            NSR=NSR+1
         ELSE
            SR(I) =1/SQRT(R(I))
         ENDIF
      ENDDO

      NSRB=0
      !$ACC PARALLEL LOOP PRESENT(RB,SRB) REDUCTION(+:NSRB) __IF_ASYNC__
      DO I=1,SIZE(R)
         IF (RB(I)<=0.0_q ) THEN
            SRB(I)=0.0_q
            NSRB=NSRB+1
         ELSE
            SRB(I) =SQRT(RB(I))
         ENDIF
      ENDDO

      JCOL=0
      !$ACC PARALLEL LOOP COLLAPSE(force:4) PRESENT(T,SRB,SR) __IF_ASYNC__
      DO J1=1,NQ,DESC(NB_)
         DO J2=1,DESC(NB_)
            JCOL=J1+J2-1
            IF(JCOL>NQ) CYCLE
            J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

            IROW=0
            DO I1=1,NP,DESC(MB_)
               DO I2=1,DESC(MB_)
                  IROW=I2+I1-1
                  IF(IROW>NP) CYCLE
                  I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                  ! T(I,J)=SRB(I)*T(I,J)*SR(J)
                  T(IROW+(JCOL-1)*DESC(LLD_))=T(IROW+(JCOL-1)*DESC(LLD_))*SRB(I)*SR(J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      !$ACC UPDATE SELF(T) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      !  Z =(X+Y)= U T1 -> A
      CALL GPGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, &
           U(1), 1, 1, DESC, T(1), 1, 1, DESC, zerobse, Z(1), 1, 1, DESC)
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

      !$ACC UPDATE DEVICE(Z) __IF_ASYNC__

      !$ACC KERNELS PRESENT(U) __IF_ASYNC__
      U=0
      !$ACC END KERNELS

      !$ACC EXIT DATA DELETE(SR,SRB,T) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      DEALLOCATE(T)
    ENDIF
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

    IF (NSR>0 .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in TDA_BSE_EIGENVECTORS: some eigenvalues in R are negative',NSR
    IF (NSRB>0.AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in TDA_BSE_EIGENVECTORS: some eigenvalues in RB are negative',NSRB

    PROFILING_STOP('tda_bse_eigenvectors_sca')

  END SUBROUTINE TDA_BSE_EIGENVECTORS_SCA

!****************** TWOELECTRON4O_GET_DIAG_SCALA  *********************
!
!> Determines the sum of the diagonal elements
!> of a matrix A
!
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_GET_DIAG_SCALA( AMAT, NCV, CSUM, MY_COMM, DESC, CTXT )

    USE constant
    IMPLICIT NONE

#ifdef single_prec_bse
    GDEFS :: AMAT(:)          !< input matrix
#endif
#ifdef double_prec_bse
    GDEF :: AMAT(:)           !< input matrix
#endif
    INTEGER :: NCV            !< dimension of global matrix
    COMPLEX(q) :: CSUM        !< sum of diagonal components
    TYPE(communic) :: MY_COMM !< communicator to globally sum the diagonal components
    INTEGER :: DESC(DLEN_)
    INTEGER, OPTIONAL :: CTXT

  ! local variables
    INTEGER NCV13, NCV24
  ! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER IROW, JCOL
    INTEGER I1, I2, J1, J2, J1RES, I1RES
    INTEGER MB,NB,LD



    CALL BLACS_GRIDINFO(DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)
    MB = DESC(MB_)
    NB = DESC(NB_)
    LD = DESC(LLD_)

   !$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(AMAT) PRIVATE(IROW,JCOL) REDUCTION(+:CSUM) __IF_ASYNC__
   DO J1=1,NQ,NB
      DO J2=1,NB
           DO I1=1,NP,MB
              DO I2=1,MB
                    NCV13=MB*MYROW+NPROW*(I1-1)+I2
                    NCV24=NB*MYCOL+NPCOL*(J1-1)+J2
                    JCOL = J1+J2-1
                    IROW = I1+I2-1
                    IF (NCV13 == NCV24 .AND. IROW <= NP .AND. JCOL <= NQ) THEN
                       CSUM = CSUM + AMAT(IROW+(JCOL-1)*LD)
                    ENDIF
              ENDDO
           ENDDO
      ENDDO
   ENDDO

   CALLMPI( M_sum_z(MY_COMM, CSUM, 1))

  END SUBROUTINE TWOELECTRON4O_GET_DIAG_SCALA

#ifdef ACC_OFFLOAD

!**********************************************************************
!
! Create cuFFT plan(s) for bands from NBMIN to NBMAX and block NSTRIP
!
!**********************************************************************

  SUBROUTINE FFT_OFFLOAD_MAKEPLAN_BSE(GRID,NBMIN,NBMAX,NSTRIP,ACC_ASYNC_Q)

  GPACC USE moffload_fft

    USE mgrid
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NBMIN, NBMAX, NSTRIP, ACC_ASYNC_Q
    TYPE (grid_3d)      :: GRID

    IF (MOD(NBMAX - NBMIN + 1,NSTRIP) == 0) THEN
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,ACC_ASYNC_Q,NSTRIP)
    ELSE
      IF (NBMAX - NBMIN + 1 > NSTRIP ) THEN
        CALL FFT_OFFLOAD_MAKEPLAN(GRID,ACC_ASYNC_Q,NSTRIP)
      ENDIF
      CALL FFT_OFFLOAD_MAKEPLAN(GRID,ACC_ASYNC_Q,MOD(NBMAX - NBMIN + 1,NSTRIP))
    ENDIF
  END SUBROUTINE

!**********************************************************************
!
! Destroy cuFFT plan(s) for bands from NBMIN to NBMAX and block NSTRIP
!
!**********************************************************************

  SUBROUTINE FFT_OFFLOAD_DESTROYPLAN_BSE(GRID,NBMIN,NBMAX,NSTRIP,ACC_ASYNC_Q)

  GPACC USE moffload_fft

    USE mgrid
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NBMIN, NBMAX, NSTRIP, ACC_ASYNC_Q
    TYPE (grid_3d)      :: GRID

    IF (MOD(NBMAX - NBMIN + 1,NSTRIP) == 0) THEN
      CALL FFT_OFFLOAD_DESTROYPLAN(GRID,ACC_ASYNC_Q,NSTRIP)
    ELSE
      IF (NBMAX - NBMIN + 1 > NSTRIP ) THEN
        CALL FFT_OFFLOAD_DESTROYPLAN(GRID,ACC_ASYNC_Q,NSTRIP)
      ENDIF
      CALL FFT_OFFLOAD_DESTROYPLAN(GRID,ACC_ASYNC_Q,MOD(NBMAX - NBMIN + 1,NSTRIP))
    ENDIF
  END SUBROUTINE

#endif

!**********************************************************************
!
! Find the number of block computed by a single rank
!
!**********************************************************************

  FUNCTION CALC_NBLOCKS_ROUND_BLOCKS(FIRST_BLOCK,LAST_BLOCK,NCPUS,COLUMNS,NMAX) RESULT(NBLOCK)
     IMPLICIT NONE
     INTEGER :: FIRST_BLOCK,LAST_BLOCK, NCPUS, COLUMNS
     INTEGER(qi8) :: NMAX
     INTEGER :: NBLOCK_MAX
     INTEGER :: NBLOCK

     NBLOCK_MAX = MIN(NMAX,MAX(1,(LAST_BLOCK - FIRST_BLOCK +1) / NCPUS))
     DO NBLOCK=NBLOCK_MAX,1,-1
        IF (MOD(COLUMNS,NBLOCK) == 0) EXIT
     ENDDO

  END FUNCTION CALC_NBLOCKS_ROUND_BLOCKS

!**********************************************************************
!
! Find the number of block computed by a single rank
!
!**********************************************************************

  FUNCTION CALC_NBLOCKS_ROUND_CPUS(NBLOCKS,NCPUS) RESULT(NRANKS)
     IMPLICIT NONE
     INTEGER :: NBLOCKS, NCPUS
     INTEGER :: NRANKS

     DO NRANKS=NCPUS,1,-1
        IF (MOD(NBLOCKS,NRANKS) == 0) EXIT
     ENDDO

  END FUNCTION CALC_NBLOCKS_ROUND_CPUS

!**********************************************************************
!
! this subroutine calculates the core orbitals for
! a specific species and stores it as a "global" array inside
! the module ...
!
!**********************************************************************
    SUBROUTINE CALCULATE_BSE_CORE_ORBITAL(P,IU6,W,WHF,T_INFO)
      USE tutor, ONLY: vtutor
      USE fock_glb, ONLY: TRANS_MATRIX_CORE_CORE, TRANS_MATRIX_CORE_VALENCE
      USE augfast, ONLY: SETUP_TRANS_MATRIX_CORE
      USE cl, ONLY: ICORELEVEL
      USE pseudo, ONLY: potcar

      ! passed variable
      INTEGER, INTENT(IN) :: IU6   ! input-output
      TYPE (potcar), TARGET, INTENT(IN) :: P(:) ! potcar structure
      TYPE (wavespin)     :: W         !wavespin from wave.F
      TYPE (wavespin)     :: WHF       !wavespin from wave.F
      TYPE (type_info)    :: T_INFO

      ! local
      TYPE (potcar), POINTER :: PP     ! pointer on the core state
      INTEGER :: CHANNELS_CORE         ! number of core channels
      INTEGER :: CH1                   ! index of core channels
      INTEGER :: RNMAX                 ! number of gridpoints
      INTEGER :: CHANNELS              ! number of channels of the valence states
      REAL(q),ALLOCATABLE :: WC(:,:)   ! AE wavefunction of the core
      REAL(q),ALLOCATABLE :: A_CORE(:) ! AE wavefunction of the core
      REAL(q),ALLOCATABLE :: A(:,:)    ! large ...
      REAL(q),ALLOCATABLE :: B(:,:)    ! ... and small component of the core wavefunction
      REAL(q),ALLOCATABLE :: EIG(:)    ! array of eigenenergies of the core channels
      INTEGER,ALLOCATABLE :: N(:)      ! array of main quantum numbers of the core channels
      INTEGER,ALLOCATABLE :: LC(:)     ! array of angular momentum quantum numbers of the core channels

      ! The following arrays are needed for the construction of the compensation charges.
      REAL(q),POINTER     :: WAE(:,:)   ! AE charges
      REAL(q),POINTER     :: PS(:,:)    ! PS charges, for the core state in question we will set PS to zero
      REAL(q),POINTER     :: QPAW(:,:,:) ! moment-restoring compensation charges on the plane wave grid
      REAL(q),POINTER     :: QPAW_FOCK(:,:,:,:) ! shape-restoring compensation charges on the plane wave grid
      INTEGER,POINTER     :: LPS(:)    ! array of angular momentum quantum numbers
      REAL(q),ALLOCATABLE :: WPS(:,:)      ! PS charges
      REAL(q) :: RMAXL(0:10)       ! outermost point for augmentation charge for each L quantum number
      INTEGER :: LMMAX_AUG
      INTEGER :: LMDIM, LMDIM_CORE
      INTEGER :: LMAX_FOCKAE_TYPE

      IF (ICORELEVEL/=2) RETURN

      IF (NT_CL==-1 .OR. NT_CL > SIZE(P) ) THEN
         CALL vtutor%bug("CLNT is -1",__FILE__,__LINE__)
      ENDIF

      PROFILING_START('calculate_bse_core_orbital')
      PP=>P(NT_CL)

      NITYP_CL = T_INFO%NITYP(NT_CL)
      ! determine the number of core channels
      CALL CL_INIT_CORE_CONF(PP,CHANNELS_CORE)
      RNMAX =PP%R%NMAX               ! number of grid points
      ALLOCATE(WC(RNMAX,CHANNELS_CORE),A(RNMAX,CHANNELS_CORE), B(RNMAX,CHANNELS_CORE), &
               N(CHANNELS_CORE), LC(CHANNELS_CORE), EIG(CHANNELS_CORE))

      ! calculate the core wavefunction and save it into A
      CALL SET_CORE_WF(PP%RHOAE, PP%POTAE, PP%R, PP%ZCORE, PP%ZVALF_ORIG, &
                        WC, N, LC, EIG, A_=A, B_=B)

      CALL CL_CLEAR_CORE_CONF

      DO CH1=1,CHANNELS_CORE
         IF (N(CH1)==N_CL .AND. LC(CH1)==L_CL) EXIT
      ENDDO

      COREEIGENVALUE = EIG(CH1)
      IF (E_CL/=0) THEN
      ! if user supplies eigenvalue in the INCAR file use the supplied value
      ! this can be calculated in DFT calculations using ICORELEVEL = 1
        COREEIGENVALUE = E_CL
      ENDIF

      IF (CH1==CHANNELS_CORE+1) THEN
         CALL vtutor%bug("cannot find the desired channel CLN and  CLL ",__FILE__,__LINE__)
      ENDIF

      ALLOCATE(A_CORE(RNMAX))
      A_CORE(:) = A(:, CH1)

      DEALLOCATE(WC, A, B, N, LC)

      ! set up QPAW_FOCK to include core-valence interaction

      CHANNELS=PP%LMAX

      IF ( CHANNELS  /= SIZE(PP%WAE,2) .OR. &
           CHANNELS  /= SIZE(PP%WPS,2) .OR. &
           CHANNELS  /=  SIZE(PP%QPAW,1)) THEN
         CALL vtutor%bug("the PP array sizes are screwed; stopping now",__FILE__,__LINE__)
      ENDIF

      ! increase the sizes of the arrays needed in the compensation charge calculation by one to include the core channel
      ALLOCATE(LPS(CHANNELS+1))
      ALLOCATE(WAE(SIZE(PP%WAE,1), SIZE(PP%WAE,2)+1))
      ALLOCATE(WPS(SIZE(PP%WPS,1), SIZE(PP%WPS,2)+1))
      ALLOCATE(QPAW(SIZE(PP%QPAW,1)+1, SIZE(PP%QPAW,2)+1, 0:SIZE(PP%QPAW,3)-1))

      ! NMAX_FOCKAE will usually agree with values given in INCAR, but it is not passed down
      ! LMAX_FOCKAE_TYPE is either given by LMAXFOCKAE in INCAR file, or by the maximum 
      !   l-quantum number of the partial waves (e.g. 2 * maximum l quantum number)
      LMAX_FOCKAE_TYPE=SIZE(PP%QPAW_FOCK,3)-1
      NMAX_FOCKAE= SIZE(PP%QPAW_FOCK,4)

      ALLOCATE(QPAW_FOCK(SIZE(PP%QPAW,1)+1, SIZE(PP%QPAW,2)+1, 0:LMAX_FOCKAE_TYPE, NMAX_FOCKAE))

      ! assign angular momentum quantum number of the core channel to the first entry of the angular momentum array
      LPS(1)= L_CL
      LPS(2:CHANNELS+1)=PP%LPS(1:CHANNELS)

      ! first entry of the AE array is the core state, the rest contains the valence states
      WAE(:,1)=A_CORE(:)
      WAE(:,2:CHANNELS+1)=PP%WAE(:,1:CHANNELS)

      ! there is no PS charge for the core state, hence the first entry of the PS array is set to zero
      WPS(:,1)=0
      WPS(:,2:CHANNELS+1)=PP%WPS(:,1:CHANNELS)

      ! call RAD_QPAW to calculate the moment-restoring compensation charges, now including core-valence interaction
      CALL RAD_QPAW( PP%R, CHANNELS+1, WAE, WPS, QPAW, LPS)

      ! the moment-restoring compensation charges of the valence states are unchanged, copy them from the PP structure
      QPAW(2:CHANNELS+1,2:CHANNELS+1,:)=PP%QPAW(1:CHANNELS,1:CHANNELS,:)

      ! call RAD_QPAW_FOCK to calculate the shape-restoring compensation charges
      CALL RAD_QPAW_FOCK( PP%R, CHANNELS+1,RMAXL, QMAX_FOCKAE(NT_CL), &
           WAE, WPS, PP%AUG, PP%AUG_FOCK, LPS, &
           QPAW, QPAW_FOCK, SIZE(PP%QPAW_FOCK,3)-1, SIZE(PP%QPAW_FOCK,4), IU6)

      ! the shape-restoring compensation charges of the valence states are unchanged, copy them from the PP structure
      QPAW_FOCK(2:CHANNELS+1,2:CHANNELS+1,:,1)=PP%QPAW_FOCK(1:CHANNELS,1:CHANNELS,:,1)

      LMMAX_AUG =SIZE(TRANS_MATRIX_FOCK,3)
      LMDIM     =SIZE(TRANS_MATRIX_FOCK,2)
      LMDIM_CORE=2*L_CL+1
      ! LMDIM is max of L and LMMAX_AUG is max of |L1+L2|
      ALLOCATE(TRANS_MATRIX_CORE_VALENCE(LMDIM_CORE, LMDIM, LMMAX_AUG))
      ALLOCATE(TRANS_MATRIX_CORE_CORE(LMDIM_CORE, LMDIM_CORE, LMMAX_AUG))

      CALL SETUP_TRANS_MATRIX_CORE(PP, TRANS_MATRIX_CORE_VALENCE, TRANS_MATRIX_CORE_CORE, &
           FAST_AUG_FOCK%LMAX(NT_CL), LMAX_FOCKAE_TYPE , NMAX_FOCKAE, LPS, QPAW, QPAW_FOCK)

      ! calculate the transition matrix elements between core and valence states, put it here since we now have LPS and WAE
      CALL CDER_CORECON_ROUTINE(LPS,WAE, P, CHANNELS,T_INFO,EIG,CH1,W,WHF%WDES)

      PROFILING_STOP('calculate_bse_core_orbital')
      DEALLOCATE(EIG)
    END SUBROUTINE CALCULATE_BSE_CORE_ORBITAL

#endif

END MODULE bse_driver
