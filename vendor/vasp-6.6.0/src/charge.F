#include "symbol.inc"
!************************* SUBROUTINE FFT_RC_SCALE *********************
! RCS:  $Id: charge.F,v 1.5 2002/08/14 13:59:37 kresse Exp $
!
!> Transforms a real-space charge density to reciprocal space
!> applying a real to complex FFT transformation
!
!***********************************************************************

  SUBROUTINE FFT_RC_SCALE(CHDENR,CHDEN,GRID)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)
    
    TYPE (grid_3d)     GRID
    
    COMPLEX(q) CHDEN(GRID%MPLWV)
    RGRID      CHDENR(GRID%RL%NP)
    
    RINPLW=1.0_q/GRID%NPLWV
    CALL RL_ADD(CHDENR,RINPLW,CHDENR,0.0_q,CHDEN,GRID)
    CALL FFT3D(CHDEN,GRID,-1)
    RETURN
  END SUBROUTINE FFT_RC_SCALE
  

MODULE charge
  USE prec
  USE base, ONLY: write_dens

CONTAINS

!************************* SUBROUTINE SOFT_CHARGE **********************
!
!> Constructs the electronic charge density according
!> to the current wavefunctions and fermi-weights
!> ~~~
!>  rho(r) = \sum_n,k phi_nk(r) phi_nk^*(r) * f_nk
!> ~~~
!> where f_nk are the single electron occupancies (FERWE)
!> the phi_nk are the pseudowavefunctions 
!
!***********************************************************************

  SUBROUTINE SOFT_CHARGE(GRID,GRID_SOFT,W,WDES, CHDEN)

__GPU USE moffload

    USE wave_high
    USE gridq
    USE hamil
    IMPLICIT COMPLEX(q) (C)

    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID,GRID_SOFT
    COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
  ! local
    TYPE (wavedes1)     WDES1
    TYPE (wavefun1)     W1
    INTEGER ISPINOR
    INTEGER, PARAMETER :: NSTRIPD=2
    TYPE (REDIS_PW_CTR),POINTER :: H_PW
    TYPE (GRIDQUANT) :: CHARGE_REAL_SPACE, CHARGE
    INTEGER I

    PROFILING_START('soft_charge')

#ifdef ACC_OFFLOAD
!$ACC WAIT IF(OFFLOAD_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
    IF (W%OVER_BAND) THEN
#ifdef MPI
       NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
       NCPU=1
#endif
       NSTRIP=MIN(NSTRIPD,WDES%NBANDS)
       CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
    ENDIF
GPACC !$ACC ENTER DATA CREATE(WDES1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(WDES, WDES1, 0)
GPACC !$ACC ENTER DATA CREATE(W1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W1) IF(OFFLOAD_ON)
    CALL NEWWAV_R(W1, WDES1)

    CALL GENERATE_GRID_QUANTITY(CHARGE, GRID_SOFT, CHDEN) 
    CALL ALLOCATE_GRID_QUANTITY_FORCE_RL(CHARGE_REAL_SPACE, GRID, GRID_SOFT, WDES%NCDIJ)

GPACC !$ACC KERNELS PRESENT(CHARGE_REAL_SPACE,CHARGE_REAL_SPACE%RG) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHARGE_REAL_SPACE,CHARGE_REAL_SPACE%RG) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CHARGE_REAL_SPACE%RG=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CHARGE_REAL_SPACE%IS_INITIALISED = .TRUE.
    CHARGE_REAL_SPACE%REALSPACE =.TRUE.

    spin: DO ISP=1,WDES%ISPIN

       kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(WDES, WDES1, NK)

          IF (W%OVER_BAND) THEN
             DO N=1,NSTRIP
                CALL REDIS_PW_START(WDES, W%CW(1,N,NK,ISP), N, H_PW)
             ENDDO
          ENDIF

          band: DO N=1,WDES%NBANDS
             CALL SETWAV(W,W1,WDES1,N,ISP)

             IF (W%OVER_BAND) THEN
                CALL REDIS_PW_STOP (WDES, W%CW(1,N,NK,ISP), N, H_PW)
                IF (N+NSTRIP<=WDES%NBANDS) &
                     CALL REDIS_PW_START(WDES, W%CW(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
             ENDIF

             WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)
             IF (WEIGHT==0) CYCLE
             CALL FFTWAV_W1(W1)
             CALL PW_CHARGE(WDES1, CHARGE_REAL_SPACE%RG(1,ISP), SIZE(CHARGE_REAL_SPACE%RG,1), W1%CR(1), W1%CR(1), WEIGHT )
          ENDDO band
       ENDDO kpoints
    ENDDO spin

    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       CALL REDIS_PW_DEALLOC(H_PW)
    ENDIF

    ! fourier-transformation of charge-density using GRID_SOFT
    ! only input data from first in-band-group is used, and all nodes
    ! are involved in the FFT, final result is distributed among nodes
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

    ! merge charge from all bands
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       DO I=1,CHARGE_REAL_SPACE%NCDIJ
#ifdef realmode
          CALLMPI( M_sum_d(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#else
          CALLMPI( M_sum_z(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#endif
       END DO
    END IF
#endif

    CALL SUMRL_GQ( CHARGE_REAL_SPACE, CHARGE, WDES%COMM_INTER)

__GPU CALL OFFLOAD_ON_PUSH(.FALSE.)

    CALL FFT_GQ(CHARGE)
      
    ! set the charge-density of unbalanced lattice vectors to 0
    CALL SETUNB_GQ(CHARGE)

__GPU CALL OFFLOAD_ON_POP

    CALL DELWAV_R(W1)
GPACC !$ACC EXIT DATA DELETE(W1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1) IF(OFFLOAD_ON)
    CALL DEALLOCATE_GRID_QUANTITY(CHARGE_REAL_SPACE)
__GPU CALL DELETE_TYPED_VAR(WDES1)
    PROFILING_STOP('soft_charge')

  END SUBROUTINE SOFT_CHARGE


!************************* SUBROUTINE SOFT_CHARGE_SYM ******************
!
!> Constructs the electronic charge density applying
!> symmetry to determine the wavefunctions at each k-point
!>
!> It should be fully compatible with the Fock like routines but
!> of course is slightly slower than #SOFT_CHARGE
!> if this routine is used to construct the charge density no
!> a posteriori symmetrization of the charge density should be required
!
!***********************************************************************

  SUBROUTINE SOFT_CHARGE_SYM(GRID,GRID_SOFT,W,WDES, CHDEN )

__GPU USE moffload

    USE wave_high
    USE gridq
    USE hamil
    USE full_kpoints
    USE pseudo
    USE lattice
    USE sym_prec
    IMPLICIT NONE

    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID,GRID_SOFT
    COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
  ! local
    INTEGER :: ISP, NK, N, ISP_IRZ
    REAL(q) :: WEIGHT
    LOGICAL :: LSHIFT
    TYPE (wavedes1)     WDES1, WDES1_IRZ
    TYPE (wavefun1)     W1
    INTEGER ISPINOR
    TYPE (REDIS_PW_CTR),POINTER :: H_PW
    TYPE (GRIDQUANT) :: CHARGE_REAL_SPACE, CHARGE
    INTEGER I

    PROFILING_START('soft_charge_sym')

#ifdef OFFLOADING
GPACC !$ACC WAIT IF(OFFLOAD_ON)
GPACC ACC_ASYNC_Q=ACC_ASYNC_SYNC
    CALL COPYIN_TYPED_VAR(KPOINTS_FULL)
#endif
    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       CALL vtutor%bug("internal error in SOFT_CHARGE_SYM: W%OVER_BAND is not supported", __FILE__, __LINE__)
    ENDIF
    CALL CHECK_FULL_KPOINTS
GPACC !$ACC ENTER DATA CREATE(WDES1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(WDES, WDES1, 0)
GPACC !$ACC ENTER DATA CREATE(W1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:W1) IF(OFFLOAD_ON)
    CALL NEWWAV(W1 , WDES1,.TRUE.)

    CALL GENERATE_GRID_QUANTITY(CHARGE, GRID_SOFT, CHDEN) 
    CALL ALLOCATE_GRID_QUANTITY_FORCE_RL(CHARGE_REAL_SPACE, GRID, GRID_SOFT, WDES%NCDIJ)

GPACC !$ACC KERNELS PRESENT(CHARGE_REAL_SPACE,CHARGE_REAL_SPACE%RG) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHARGE_REAL_SPACE,CHARGE_REAL_SPACE%RG) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CHARGE_REAL_SPACE%RG=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CHARGE_REAL_SPACE%IS_INITIALISED = .TRUE.
    CHARGE_REAL_SPACE%REALSPACE =.TRUE.

    spin: DO ISP=1,WDES%ISPIN

       kpoints: DO NK=1,KPOINTS_FULL%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(WDES, WDES1, NK)
          CALL SETWDES(WDES, WDES1_IRZ, KPOINTS_FULL%NEQUIV(NK))
          
          ISP_IRZ=ISP
          IF (KPOINTS_FULL%SPINFLIP(NK)==1) THEN
             ISP_IRZ=3-ISP
          ENDIF
         
          band: DO N=1,WDES%NBANDS
             WEIGHT=WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*W%FERWE(N,KPOINTS_FULL%NEQUIV(NK),ISP_IRZ)

             IF (WEIGHT==0) CYCLE

             IF (NK <= WDES%NKPTS) THEN
                CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
                CALL FFTWAV_W1(W1)
             ELSE
#ifndef gammareal
           ! 
           ! symmetry must be considered if the wavefunctions for this
           ! k-point NK (containing all k-points in the entire BZ) 
           ! are not stored in W
           !
                LSHIFT=.FALSE.
                IF ((ABS(KPOINTS_FULL%TRANS(1,NK)) >TINY) .OR. &
                     (ABS(KPOINTS_FULL%TRANS(2,NK))>TINY) .OR. &
                     (ABS(KPOINTS_FULL%TRANS(3,NK))>TINY)) LSHIFT=.TRUE.
                CALL W1_ROTATE_AND_FFT_NO_PROJ(W1, ELEMENT(W, WDES1_IRZ, N, ISP_IRZ), LSHIFT)
#endif
             ENDIF

             CALL PW_CHARGE(WDES1, CHARGE_REAL_SPACE%RG(1,ISP), SIZE(CHARGE_REAL_SPACE%RG,1), W1%CR(1), W1%CR(1), WEIGHT )
          ENDDO band
       ENDDO kpoints
    ENDDO spin

    ! fourier-transformation of charge-density using GRID_SOFT
    ! only input data from first in-band-group is used, and all nodes
    ! are involved in the FFT, final result is distributed among nodes
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

    ! merge charge from all bands
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       DO I=1,CHARGE_REAL_SPACE%NCDIJ
#ifdef realmode
          CALLMPI( M_sum_d(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#else
          CALLMPI( M_sum_z(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#endif
       END DO
    END IF
#endif

    CALL SUMRL_GQ( CHARGE_REAL_SPACE, CHARGE, WDES%COMM_INTER)
    
__GPU CALL OFFLOAD_ON_PUSH(.FALSE.)

    CALL FFT_GQ(CHARGE)
      
    ! set the charge-density of unbalanced lattice vectors to 0
    CALL SETUNB_GQ(CHARGE)

__GPU CALL OFFLOAD_ON_POP

    CALL DELWAV(W1, .TRUE.)
GPACC !$ACC EXIT DATA DELETE(W1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W1) IF(OFFLOAD_ON)
    CALL DEALLOCATE_GRID_QUANTITY(CHARGE_REAL_SPACE)
#ifdef OFFLOADING
    CALL DELETE_TYPED_VAR(WDES1)
    CALL DELETE_TYPED_VAR(KPOINTS_FULL)
#endif
    PROFILING_STOP('soft_charge_sym')

  END SUBROUTINE SOFT_CHARGE_SYM


!************************* SUBROUTINE SOFT_CHARGE1 *********************
!
!> Constructs the electronic charge density from two sets of
!> wavefunctions
!> ~~~
!>  rho(r) = sum_n,k (phi(0)_nk(r) phi(1)_nk*(r) + c.c. )* f(0)_nk
!>                  + phi(0)_nk(r) phi(0)_nk^*(r) * f(1)_nk
!> ~~~
!> f(0) (W0\%FERWE) are the initial occupancies
!> f(1) (W1\%FERWE) are the first order change of the single electron
!>                  occupancies
!>
!> (This routine is used for linear response theory)
!
!***********************************************************************

      SUBROUTINE SOFT_CHARGE1(GRID,GRID_SOFT,W0,W1,WDES, CHDEN)

GPACC USE moffload

      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE wave_mpi
      IMPLICIT NONE 

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W0
      TYPE (wavespin)    W1
      TYPE (grid_3d)     GRID,GRID_SOFT
      COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)

    ! work arrays
      RGRID,   ALLOCATABLE :: CDWORK(:,:)
      COMPLEX(q), ALLOCATABLE :: CPTDUM0(:), CPTDUM1(:)
      INTEGER ISPINOR, ISPINOR_, ISP, NK, N, NPL, M, MM, MM_, I
      REAL(q) :: WEIGHT0, WEIGHT1

      PROFILING_START('soft_charge1')

      ! MPLWV is the allocation in complex words
      ! hence if CDWORK is REAL one needs to double the allocation 
      ALLOCATE(CDWORK(GRID%MPLWV*2,WDES%NCDIJ), &
           CPTDUM0(GRID%MPLWV*WDES%NRSPINORS),CPTDUM1(GRID%MPLWV*WDES%NRSPINORS))

#ifdef ACC_OFFLOAD
      CALL OFFLOAD_ON_PUSH(.TRUE.)

      CALL COPYIN_TYPED_VAR(GRID)
      CALL COPYIN_TYPED_VAR(GRID_SOFT)
      CALL COPYIN_TYPED_VAR(WDES)
      CALL COPYIN_TYPED_VAR(W0)
      CALL COPYIN_TYPED_VAR(W1)
#endif
!$ACC ENTER DATA CREATE(CDWORK,CPTDUM0,CPTDUM1) __IF_ASYNC__

!$ACC KERNELS PRESENT(CDWORK) __IF_ASYNC__
      CDWORK=0
!$ACC END KERNELS
!=======================================================================
! loop over k-points and bands
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         band: DO N=1,WDES%NBANDS

            WEIGHT0=WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)
            WEIGHT1=WDES%RSPIN*WDES%WTKPT(NK)*W1%FERWE(N,NK,ISP)
            IF (WEIGHT0==0 .AND. WEIGHT1==0) CYCLE

            NPL=WDES%NGVECTOR(NK)
!=======================================================================
! fourier-transformation of wave-function
! sum up  real space charge density
!=======================================================================
            DO ISPINOR=0,WDES%NRSPINORS-1
               CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM0(1+ISPINOR*GRID%MPLWV),W0%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
               CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM1(1+ISPINOR*GRID%MPLWV),W1%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
            ENDDO

!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,WDES,CDWORK,CPTDUM0,CPTDUM1) PRIVATE(MM,MM_) __IF_ASYNC__
            spinor: DO ISPINOR=0,WDES%NRSPINORS-1
            DO ISPINOR_=0,WDES%NRSPINORS-1
               DO M=1,GRID%RL%NP
                 MM =M+ISPINOR *GRID%MPLWV
                 MM_=M+ISPINOR_*GRID%MPLWV
                 CDWORK(M,ISP+ISPINOR_+2*ISPINOR)=CDWORK(M,ISP+ISPINOR_+2*ISPINOR)+ &
                      CPTDUM0(MM)*CONJG(CPTDUM0(MM_))*WEIGHT1+ &
                      CPTDUM0(MM)*CONJG(CPTDUM1(MM_))*WEIGHT0+ &
                      CPTDUM1(MM)*CONJG(CPTDUM0(MM_))*WEIGHT0
               ENDDO
            ENDDO
            ENDDO spinor
         ENDDO band

      ENDDO kpoints
      ENDDO spin

!=======================================================================
! fourier-transformation of charge-density using GRID_SOFT (see above)
!=======================================================================
!$ACC ENTER DATA CREATE(CHDEN) __IF_ASYNC__
      DO I=1,WDES%NCDIJ
! now merge the chargedensity from all nodes

!PK Reduce onto first KINTER nodes then broadcast
#ifdef realmode
         CALLMPI( M_sum_d(WDES%COMM_INTER, CDWORK(1,I), GRID%RL%NP))
         CALLMPI( M_sum_d(WDES%COMM_KINTER, CDWORK(1,I), GRID%RL%NP))
#else
         CALLMPI( M_sum_z(WDES%COMM_INTER, CDWORK(1,I), GRID%RL%NP))
         CALLMPI( M_sum_z(WDES%COMM_KINTER, CDWORK(1,I), GRID%RL%NP))
#endif
         CALL FFT_RC_SCALE(CDWORK(1,I),CHDEN(1,I),GRID_SOFT)
! set the charge-density of unbalanced lattic-vectors to 0
         CALL SETUNB(CHDEN(1,I),GRID_SOFT)
      ENDDO

!$ACC EXIT DATA DELETE(CDWORK,CPTDUM0,CPTDUM1) __IF_ASYNC__
!$ACC EXIT DATA COPYOUT(CHDEN) __IF_ASYNC__
#ifdef ACC_OFFLOAD
      CALL DELETE_TYPED_VAR(W0)
      CALL DELETE_TYPED_VAR(W1)
      CALL DELETE_TYPED_VAR(GRID)
      CALL DELETE_TYPED_VAR(GRID_SOFT)
      CALL DELETE_TYPED_VAR(WDES)

      CALL OFFLOAD_ON_POP
#endif
      DEALLOCATE(CDWORK, CPTDUM0, CPTDUM1)

      PROFILING_STOP('soft_charge1')

      RETURN
      END SUBROUTINE


!************************* SUBROUTINE RHOAT0 ***************************
!
!> Calculates the term G^2 in the Taylor expansion of the atomic
!> charge densities and multiplies this term with the number of
!> atoms per type and the electronic field constant
!
!***********************************************************************

      SUBROUTINE RHOAT0(P,T_INFO, BETATO,OMEGA)
      USE prec
      USE pseudo
      USE poscar
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)

      BETATO=0._q
      DO 200 NT=1,T_INFO%NTYP
        IF (P(NT)%PSPRHO(1)==0) GOTO 200
        IMAX=4
        DQ=P(NT)%PSGMAX/NPSPTS
        B1=0
        B2=0
        DO 100 I=1,IMAX
          A=(P(NT)%PSPRHO(I+1)-P(NT)%PSPRHO(1))/(DQ*I)**2
          B1=B1+A
          B2=B2+A*A
  100   CONTINUE
        B1=B1/IMAX
        B2=B2/IMAX
        BETATO=BETATO+B1*T_INFO%NITYP(NT)*T_INFO%VCA(NT)*EDEPS/OMEGA

  200 CONTINUE
      RETURN
      END SUBROUTINE


!************************* SUBROUTINE RHOATO ***************************
!
!> Calculates the charge density and its derivatives corresponding
!> to overlapping atoms
!>
!> @param LFOUR = .FALSE. set up charge density in reciprocal space\n
!>              = .TRUE.  set up charge density in real space
!> @param LPAR  = .FALSE. use PSPRHO\n
!>              = .TRUE.  use PSPCOR
!
!***********************************************************************

      SUBROUTINE RHOATO(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CHDER)

__GPU USE moffload_struct_def

      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%RC%NP),CHDER(GRIDC%RC%NP)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR,LPAR
! local variables
      REAL(q), POINTER :: PRHO(:)

GPACC !$ACC KERNELS PRESENT(CHTOT,CHDER) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CHTOT,CHDER) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CHTOT=0
      CHDER=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP
       IF (LPAR) THEN
         PRHO=>P(NT)%PSPCOR
       ELSE
         PRHO=>P(NT)%PSPRHO
       ENDIF
       IF (ASSOCIATED(PRHO)) THEN
!=======================================================================
! calculate the scaling factor ARGSC that converts the magnitude of a
! reciprocal lattice vector to the correponding position in the
! pseudopotential arrays
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3

GPACC !$ACC PARALLEL LOOP PRIVATE(N1,NC,N2,N3,GX,GY,GZ,G,ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3) &
GPACC !$ACC& PRESENT(GRIDC,B,PRHO,CHTOT,CSTRF,CHDER) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(N1,NC,N2,N3,GX,GY,GZ,G,ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3) &
GPOMP !$OMP MAP(PRESENT,TO:GRIDC,B,PRHO,CHTOT,CSTRF,CHDER) IF(OFFLOAD_ON)
      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal latice vector to a position
! in the charge-dens. array  and interpolate the atomic-chargedensity
!=======================================================================
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=PRHO(NADDR-1)
          V2=PRHO(NADDR)
          V3=PRHO(NADDR+1)
          V4=PRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          CHTOT(N)=CHTOT(N)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N,NT)
          CHDER(N)=CHDER(N)+(T1+REM*(2*T2+REM*3*T3))*ARGSC*CSTRF(N,NT)
        ELSE IF (G==0) THEN
          CHTOT(N)=CHTOT(N)+PRHO(1)*CSTRF(N,NT)
          CHDER(N)=0
        ENDIF
      ENDDO
      ENDIF
      ENDDO typ
!=======================================================================
! set the charge-density of unbalanced lattice-vectors to 0
! and transform the charge-density to real space
!=======================================================================
      CALL SETUNB(CHTOT(1),GRIDC)
      CALL SETUNB(CHDER(1),GRIDC)

      IF (LFOUR) THEN
        CALL FFT3D(CHTOT,GRIDC,1)
      ENDIF

      RETURN
      END SUBROUTINE
!
!> Small subroutine which allocates the charge derivative array 
!> CWORK explicitly and calls #RHOATO
! 
      SUBROUTINE RHOATO_WORK(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT)

__GPU USE moffload_struct_def

      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      IMPLICIT NONE


      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%MPLWV)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR,LPAR
! work arrays
      COMPLEX(q),ALLOCATABLE :: CWORK(:)

      ALLOCATE(CWORK(1:GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK) IF(OFFLOAD_ON)
      CALL RHOATO(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CWORK)
GPACC !$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK) IF(OFFLOAD_ON)
      DEALLOCATE(CWORK)
      
      RETURN
      END SUBROUTINE


!************************* SUBROUTINE RHOATO_PARTICLE_MESH *************
!
!> Calculates the charge density corresponding to overlapping atoms
!> using the particle mesh method
!>
!> @param LFOUR = .FALSE. set up charge density in reciprocal space\n
!>              = .TRUE.  set up charge density in real space
!> @param LPAR  = .FALSE. use PSPRHO\n
!>              = .TRUE.  use PSPCOR
!>
!> @note Does NOT calculate the derivatives!
!
!***********************************************************************

      SUBROUTINE RHOATO_PARTICLE_MESH(LFOUR,LPAR,GRIDC,LATT_CUR,T_INFO,INFO,P,CHTOT,IU6)
      USE base
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      USE lattice
      IMPLICIT NONE
      LOGICAL LFOUR,LPAR
      TYPE (grid_3d) GRIDC
      TYPE (latt) LATT_CUR
      TYPE (type_info) T_INFO
      TYPE (info_struct) INFO
      TYPE (potcar) P(T_INFO%NTYP)
      COMPLEX(q) CHTOT(GRIDC%RC%NP)
      INTEGER IU6
      ! local variables
      INTEGER NT
      REAL(q) G1,G2,GMAX
      REAL(q), POINTER :: PRHO(:)

      INTEGER, PARAMETER :: NFFT=32768

#ifdef debug
      REAL(q) RHOT
      REAL(q), EXTERNAL :: RHO0
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'RHOATO_PARTICLE_MESH: Total number of electrons =',RHOT
#endif

      DO NT=1,T_INFO%NTYP
         NULLIFY(P(NT)%USESPL)
         IF (LPAR) THEN
            PRHO=>P(NT)%PSPCOR
         ELSE
            PRHO=>P(NT)%PSPRHO
         ENDIF
         IF (.NOT.ASSOCIATED(PRHO)) CYCLE
         ! momentum cutoffs for Fourier filtering
         GMAX=SQRT(INFO%ENMAX/HSQDTM)
         IF (INFO%SZPREC(1:1)=='s') THEN
            G1=GMAX
            G2=MIN(3.0_q*GMAX,P(NT)%PSGMAX)/2
         ELSE
            G1=2.0_q*GMAX
            G2=MIN(3.0_q*GMAX,P(NT)%PSGMAX)
         ENDIF
         ! Fourier-filter atomic charge density, transform to real space
         ! and calculate spline interpolation
         ALLOCATE(P(NT)%ATOSPL(NFFT/2,5))
         CALL ATORHO_FILTER(G1,G2,NPSPTS,PRHO,P(NT)%PSGMAX/NPSPTS,NFFT,P(NT)%ATOSPL,P(NT)%RCUTATO,-1)
         ! set the arrays and cutoff that are used by RHOADD
         P(NT)%USESPL=>P(NT)%ATOSPL
         P(NT)%USEZ=PRHO(1)
         P(NT)%USECUT=P(NT)%RCUTATO
      ENDDO

      CHTOT=0
      ! put the atomic charge density of all ions on the grid 
      CALL RHOADD(T_INFO,LATT_CUR,P,GRIDC,CHTOT)

#ifdef debug
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'RHOATO_PARTICLE_MESH: Total number of electrons =',RHOT
#endif
      ! if requested, take the charge density to real space
      IF (LFOUR) THEN
         CALL FFT3D(CHTOT,GRIDC,1)
      ENDIF

      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%ATOSPL)) THEN
            DEALLOCATE(P(NT)%ATOSPL)
            NULLIFY(P(NT)%ATOSPL)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE RHOATO_PARTICLE_MESH


!******************** SUBROUTINE RHOADD ********************************
!
!***********************************************************************
      SUBROUTINE RHOADD(T_INFO,LATT_CUR,P,GRIDC,CHTOT,QTOT)
      USE ini
      USE prec
      USE poscar
      USE lattice
      USE mgrid
      USE pseudo
      USE constant
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (grid_3d) GRIDC
      COMPLEX(q) CHTOT(GRIDC%MPLWV)
      REAL(q), OPTIONAL :: QTOT(T_INFO%NIONS)
      ! local variables
      REAL(q) QTOT_(T_INFO%NIONS)

#ifdef debug
      REAL(q) RHOT
      REAL(q), EXTERNAL :: RHO0
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'Total number of electrons =',RHOT
#endif

      ! take total charge to real space
      CALL FFT3D(CHTOT(1),GRIDC,1)
      
      CALL RHOADD_RL(T_INFO,LATT_CUR,P,GRIDC,QTOT_,CHTOT(1))

      ! and back to reciprocal space
      CALL FFT_RC_SCALE(CHTOT(1),CHTOT(1),GRIDC)
      CALL SETUNB_COMPAT(CHTOT(1),GRIDC)

#ifdef debug
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'Total number of electrons =',RHOT
#endif
      IF (PRESENT(QTOT)) QTOT=QTOT_

      RETURN
      END SUBROUTINE RHOADD


!************************* SUBROUTINE RHOPAR ***************************
!
!> Calculates the charge-density (or partial-core charge densities)
!> from overlapping atoms in real space and stores it in an real array
!>
!> @note The DENCOR array is not necessarily sufficiently large
!> to allow for a 3d-FFT, hence a work array is temporarily allocated
!
!***********************************************************************

      SUBROUTINE RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENS,IU6)

GPACC USE moffload_struct_def

      USE prec
      USE base
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE lattice
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (info_struct) INFO
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (grid_3d) GRIDC
      TYPE (latt) LATT_CUR
      RGRID DENS(GRIDC%RL%NP)
      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      INTEGER IU6
! work arrays
      COMPLEX(q),ALLOCATABLE:: CWORK1(:),CWORK2(:)

      ALLOCATE(CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
!$ACC ENTER DATA CREATE(CWORK1,CWORK2) __IF_ASYNC__
      IF(INFO%TURBO==0)THEN
         CALL RHOATO(.TRUE.,.TRUE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CWORK1,CWORK2)
      ELSE
         CALL RHOATO_PARTICLE_MESH(.TRUE.,.TRUE.,GRIDC,LATT_CUR,T_INFO,INFO,P,CWORK1,IU6)
      END IF
      CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,DENS,GRIDC)

!$ACC EXIT DATA DELETE(CWORK1,CWORK2) __IF_ASYNC__
      DEALLOCATE(CWORK1,CWORK2)

      RETURN
      END SUBROUTINE
      END MODULE


!************************* SUBROUTINE MRHOATO **************************
!
!> Calculates the magnetization density and its derivatives
!> corresponding to overlapping atoms
!>
!> @param LFOUR = .TRUE. transform density to reciprocal space
!
!***********************************************************************

      SUBROUTINE MRHOATO(LFOUR,GRIDC,T_INFO,B,P, CHTOT, ISPIN)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      INTEGER ISPIN
      COMPLEX(q)   CHTOT(GRIDC%MPLWV,ISPIN)
      REAL(q)      B(3,3)
      LOGICAL   LFOUR
      INTEGER ISP
! local variables
      INTEGER NT,NADDR,NI,N,N3,N2,N1,NA(T_INFO%NIONS),NIS,NC
      REAL(q) V1,V2,V3,V4,T0,T1,T2,T3,REM,ARGSC,PSGMA2, &
              GX,GY,GZ,G,ARG,AMAG,FNORM
! work arrays
      COMPLEX(q),ALLOCATABLE:: CSTRF(:)

      IF (ISPIN==0) RETURN

      ALLOCATE(CSTRF(GRIDC%RC%NP))

      CHTOT=0
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
       IF (P(NT)%PSPRHO(1)==0) GOTO 200

      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3
      FNORM=1._q/P(NT)%ZVALF

 ion: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

    ! set up phase factor for this ion
      CALL STUFAK_ONE(GRIDC,1,T_INFO%POSION(1,NI),T_INFO%VCA(NT),CSTRF)

    ! AMAG is the magnitude of the  magnetization for this ion
    ! and for the present direction
spin: DO ISP=1,ISPIN

      AMAG=T_INFO%ATOMOM(ISP+(NI-1)*ISPIN)*FNORM

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
        ! calculate the magnitude of the reciprocal lattice vector
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
        ! convert the magnitude of the reciprocal latice vector to a position
        ! in the charge-dens. array  and interpolate the atomic-chargedensity
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=P(NT)%PSPRHO(NADDR-1)
          V2=P(NT)%PSPRHO(NADDR)
          V3=P(NT)%PSPRHO(NADDR+1)
          V4=P(NT)%PSPRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          ! voila: phase factor  * charge * AMAG
          CHTOT(N,ISP)=CHTOT(N,ISP)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N)*AMAG
        ELSE IF (G==0) THEN
          CHTOT(N,ISP)=CHTOT(N,ISP)+P(NT)%PSPRHO(1)*CSTRF(N)*AMAG
        ENDIF
      ENDDO
      ENDDO spin
      ENDDO ion
  200 NIS=NIS+T_INFO%NITYP(NT)
!=======================================================================
      ENDDO typ
!=======================================================================

! set the charge-density of unbalanced lattic-vectors to 0
! and transform the charge-density to real space

      DO ISP=1,ISPIN
         CALL SETUNB(CHTOT(1,ISP),GRIDC)
         IF (LFOUR) THEN
            CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
         ENDIF
      ENDDO

      DEALLOCATE( CSTRF )

      RETURN
      END SUBROUTINE


!************************* SUBROUTINE RHOADD_RL ************************
!
!> Evaluate the atomic or pseudo or ... charge density in
!> P(:)\%USESPL (pseudo_struct_def::potcar::USESPL) on the grid GRIDC,
!> cutting of at a cube of length 2*P(:)\%USECUT, and
!> normalizing to P(:)\%USEZ for each ion
!
!***********************************************************************

      SUBROUTINE RHOADD_RL(T_INFO,LATT_CUR,P,GRIDC,QITOT,DHTOT)
      USE prec
      USE constant
      USE poscar
      USE lattice
      USE mgrid
      USE mpimy
      USE pseudo
      USE radial
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (grid_3d) GRIDC
      REAL(q) QITOT(T_INFO%NIONS)
      RGRID DHTOT(DIMREAL(GRIDC%MPLWV))
      ! local variables
      INTEGER I
      INTEGER NIS,NT,NI,K
      INTEGER N1LOW,N2LOW,N3LOW,N1HI,N2HI,N3HI
      INTEGER N1,N2,N3,N1P,N2P,N3P,NCOL,IND,IND1
      REAL(q) F1,F2,F3,D1,D2,D3,X1,X2,X3
      REAL(q) XC,YC,ZC,D,DX
      REAL(q) QR,BJ,QSUM,DQ,QN,RHOT,QSM1
      INTEGER NTOT(T_INFO%NIONS)
      REAL(q) QTOT(T_INFO%NIONS)
      REAL(q) QSPH(T_INFO%NIONS)
      INTEGER, PARAMETER :: NQ=2
      REAL(q) QQ(NQ),A(NQ),VAL,VALDER,RCUT
      REAL(q), PARAMETER :: TINY = 1.E-5_q

      NTOT=0; QTOT=0; QSPH=0

      NIS=1

      type1: DO NT=1,T_INFO%NTYP
         IF (.NOT.ASSOCIATED(P(NT)%USESPL)) THEN
            NIS = NIS+T_INFO%NITYP(NT); CYCLE type1
         ENDIF

         ions1: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

            !=======================================================================
            ! find lattice points contained within the cutoff-sphere
            !=======================================================================
            F1=1._q/GRIDC%NGX
            F2=1._q/GRIDC%NGY
            F3=1._q/GRIDC%NGZ

            !-----------------------------------------------------------------------
            ! restrict loop to points contained within a cubus around the ion
            !-----------------------------------------------------------------------
            D1= P(NT)%USECUT*LATT_CUR%BNORM(1)*GRIDC%NGX
            D2= P(NT)%USECUT*LATT_CUR%BNORM(2)*GRIDC%NGY
            D3= P(NT)%USECUT*LATT_CUR%BNORM(3)*GRIDC%NGZ

            N3LOW= INT(T_INFO%POSION(3,NI)*GRIDC%NGZ-D3+10*GRIDC%NGZ+.99_q)-10*GRIDC%NGZ
            N2LOW= INT(T_INFO%POSION(2,NI)*GRIDC%NGY-D2+10*GRIDC%NGY+.99_q)-10*GRIDC%NGY
            N1LOW= INT(T_INFO%POSION(1,NI)*GRIDC%NGX-D1+10*GRIDC%NGX+.99_q)-10*GRIDC%NGX
            
            N3HI = INT(T_INFO%POSION(3,NI)*GRIDC%NGZ+D3+10*GRIDC%NGZ)-10*GRIDC%NGZ
            N2HI = INT(T_INFO%POSION(2,NI)*GRIDC%NGY+D2+10*GRIDC%NGY)-10*GRIDC%NGY
            N1HI = INT(T_INFO%POSION(1,NI)*GRIDC%NGX+D1+10*GRIDC%NGX)-10*GRIDC%NGX
            
            IND=0; QSUM=0; QSM1=0

#ifdef MPI
            !-----------------------------------------------------------------------
            ! loop over cubus
            ! MPI version z ist the fast index
            !-----------------------------------------------------------------------
            DO N2=N2LOW,N2HI
               X2=(N2*F2-T_INFO%POSION(2,NI))
               N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
               
               DO N1=N1LOW,N1HI
                  X1=(N1*F1-T_INFO%POSION(1,NI))
                  N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                  
                  NCOL=GRIDC%RL%INDEX(N1P,N2P)
                  IF (NCOL==0) CYCLE ! not on local node, move on
                  IF (GRIDC%RL%I2(NCOL) /= N1P+1 .OR. GRIDC%RL%I3(NCOL) /= N2P+1) THEN
                     CALL vtutor%bug("RHOADD: internal ERROR: " // str(GRIDC%RL%I2(NCOL)) // " " // &
                        str(N1P+1) // " " // str(GRIDC%RL%I3(NCOL)) // " " // str(N2P+1), __FILE__, __LINE__)
                  ENDIF
                  !OCL SCALAR
                  DO N3=N3LOW,N3HI
                     X3=(N3*F3-T_INFO%POSION(3,NI))
                     N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
                     IND1=N3P+(NCOL-1)*GRIDC%NGZ+1          
#else
            !-----------------------------------------------------------------------
            ! loop over cubus around one ion
            ! conventional version x is fast index
            !-----------------------------------------------------------------------
            DO N3=N3LOW,N3HI
               X3=(N3*F3-T_INFO%POSION(3,NI))
               N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
               
               DO N2=N2LOW,N2HI
                  X2=(N2*F2-T_INFO%POSION(2,NI))
                  N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
                  
                  NCOL=GRIDC%RL%INDEX(N2P,N3P)
                  
                  DO N1=N1LOW,N1HI
                     X1=(N1*F1-T_INFO%POSION(1,NI))
                     N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                     IND1=N1P+(NCOL-1)*GRIDC%NGX+1
#endif
            
                     XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                     YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                     ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
                     
                     D=SQRT(XC*XC+YC*YC+ZC*ZC)
                     
                     IF (D<=P(NT)%USECUT) THEN
                        IND=IND+1
                                                 
                        K=D/P(NT)%USESPL(2,1)+1
                        DX=D-P(NT)%USESPL(K,1)
                        VAL   =((P(NT)%USESPL(K,5)*DX+P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3))*DX+P(NT)%USESPL(K,2)
                        VALDER=(3.0_q*P(NT)%USESPL(K,5)*DX+2.0_q*P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3)

                        QSUM=QSUM+VAL
            
                        ! sum total charge in sphere
                        QSM1=QSM1+DHTOT(IND1)
                     ENDIF
                     
                  ENDDO
               ENDDO
            ENDDO
                     
            QTOT(NI)=QSUM
            NTOT(NI)=IND
            QSPH(NI)=QSM1

         ENDDO ions1
         NIS = NIS+T_INFO%NITYP(NT)
      ENDDO type1

      CALLMPI( M_sum_i(GRIDC%COMM, NTOT, T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM, QTOT, T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM, QSPH, T_INFO%NIONS))
      QITOT=QTOT

#ifdef debug
      DO NI=1,T_INFO%NIONS
         ! WRITE(*,'(I6,I8,F14.7)') NI,NTOT(NI),QTOT(NI)/NTOT(NI)!*LATT_CUR%OMEGA*F1*F2*F3
         WRITE(6,'(a,I6,10F14.7)') 'tot',NI,QTOT(NI)*LATT_CUR%OMEGA*F1*F2*F3,QSPH(NI)/GRIDC%NPLWV
      ENDDO
#endif

      NIS=1
      type2: DO NT=1,T_INFO%NTYP

         IF (.NOT.ASSOCIATED(P(NT)%USESPL)) THEN
            NIS = NIS+T_INFO%NITYP(NT); CYCLE type2
         ENDIF

         ions2: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

            !-----------------------------------------------------------------------
            ! restrict loop to points contained within a cubus around the ion
            !-----------------------------------------------------------------------
            D1= P(NT)%USECUT*LATT_CUR%BNORM(1)*GRIDC%NGX
            D2= P(NT)%USECUT*LATT_CUR%BNORM(2)*GRIDC%NGY
            D3= P(NT)%USECUT*LATT_CUR%BNORM(3)*GRIDC%NGZ
            
            N3LOW= INT(T_INFO%POSION(3,NI)*GRIDC%NGZ-D3+10*GRIDC%NGZ+.99_q)-10*GRIDC%NGZ
            N2LOW= INT(T_INFO%POSION(2,NI)*GRIDC%NGY-D2+10*GRIDC%NGY+.99_q)-10*GRIDC%NGY
            N1LOW= INT(T_INFO%POSION(1,NI)*GRIDC%NGX-D1+10*GRIDC%NGX+.99_q)-10*GRIDC%NGX
            
            N3HI = INT(T_INFO%POSION(3,NI)*GRIDC%NGZ+D3+10*GRIDC%NGZ)-10*GRIDC%NGZ
            N2HI = INT(T_INFO%POSION(2,NI)*GRIDC%NGY+D2+10*GRIDC%NGY)-10*GRIDC%NGY
            N1HI = INT(T_INFO%POSION(1,NI)*GRIDC%NGX+D1+10*GRIDC%NGX)-10*GRIDC%NGX
            
            IND=0; QSUM=0
            !-----------------------------------------------------------------------
            ! loop over cubus
            ! MPI version z ist the fast index
            !-----------------------------------------------------------------------
#ifdef MPI
            DO N2=N2LOW,N2HI
               X2=(N2*F2-T_INFO%POSION(2,NI))
               N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
               
               DO N1=N1LOW,N1HI
                  X1=(N1*F1-T_INFO%POSION(1,NI))
                  N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                  
                  NCOL=GRIDC%RL%INDEX(N1P,N2P)
                  IF (NCOL==0) CYCLE ! not on local node, move on
                  IF (GRIDC%RL%I2(NCOL) /= N1P+1 .OR. GRIDC%RL%I3(NCOL) /= N2P+1) THEN
                     CALL vtutor%bug("RHOADD: internal ERROR: " // str(GRIDC%RL%I2(NCOL)) // " " // &
                        str(N1P+1) // " " // str(GRIDC%RL%I3(NCOL)) // " " // str(N2P+1), __FILE__, __LINE__)
                  ENDIF
                  !OCL SCALAR
                  DO N3=N3LOW,N3HI
                     X3=(N3*F3-T_INFO%POSION(3,NI))
                     N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
                     IND=N3P+(NCOL-1)*GRIDC%NGZ+1
#else
            !-----------------------------------------------------------------------
            ! loop over cubus around one ion
            ! conventional version x is fast index
            !-----------------------------------------------------------------------
            DO N3=N3LOW,N3HI
               X3=(N3*F3-T_INFO%POSION(3,NI))
               N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
               
               DO N2=N2LOW,N2HI
                  X2=(N2*F2-T_INFO%POSION(2,NI))
                  N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
                  
                  NCOL=GRIDC%RL%INDEX(N2P,N3P)
                  
                  DO N1=N1LOW,N1HI
                     X1=(N1*F1-T_INFO%POSION(1,NI))
                     N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                     IND=N1P+(NCOL-1)*GRIDC%NGX+1
#endif        
             
                     XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                     YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                     ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
                     
                     D=SQRT(XC*XC+YC*YC+ZC*ZC)
                     
                     IF (D<=P(NT)%USECUT) THEN

                        K=D/P(NT)%USESPL(2,1)+1
                        DX=D-P(NT)%USESPL(K,1)
                        VAL   =((P(NT)%USESPL(K,5)*DX+P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3))*DX+P(NT)%USESPL(K,2)
                        VALDER=(3.0_q*P(NT)%USESPL(K,5)*DX+2.0_q*P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3)

                        QSUM=QSUM+VAL/QTOT(NI)
                        DHTOT(IND)=DHTOT(IND)+VAL/QTOT(NI)*P(NT)%USEZ*GRIDC%NPLWV
                     ENDIF
                     
                  ENDDO
               ENDDO
            ENDDO
                                 
            QTOT(NI)=QSUM
            
         ENDDO ions2
         NIS = NIS+T_INFO%NITYP(NT)
      ENDDO type2
      
      CALLMPI( M_sum_i(GRIDC%COMM, NTOT, T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM, QTOT, T_INFO%NIONS))

#ifdef debug
      RHOT=0
      DO IND=1,GRIDC%RL%NP
         RHOT=RHOT+DHTOT(IND)
      ENDDO
      CALLMPI( M_sum_d(GRIDC%COMM, RHOT, 1))
      
      DO NI=1,T_INFO%NIONS
         WRITE(*,'(I6,I8,F14.7)') NI,NTOT(NI),QTOT(NI)!*LATT_CUR%OMEGA*F1*F2*F3
      ENDDO
      WRITE(*,'(A,F20.10)') 'Total number of electrons =',RHOT/GRIDC%NPLWV
#endif

      RETURN
      END SUBROUTINE RHOADD_RL


!************************* SUBROUTINE RHODER ***************************
!
!> Calculates the derivative of the charge density with respect to
!> one ionic coordinate r_NI:
!> ~~~
!>        f_i(r)=d rho/dr*(r-R_i)/|r-R_i|
!> ~~~
!> where rho is taken from P\%USESPL,
!> and integrates
!> ~~~
!>        F_i = \int dr f_i(r) V(r)
!> ~~~
!> where V is POTTOT, F_i and CIFOR,
!> to get the force acting on ion NI
!
!**********************************************************************

      SUBROUTINE RHODER(T_INFO,LATT_CUR,P,NI,GRIDC,QTOT,POTTOT,IFOR)
      USE prec
      USE constant
      USE poscar
      USE lattice
      USE mgrid
      USE mpimy
      USE pseudo
      USE radial
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (potcar) P
      INTEGER NI
      TYPE (grid_3d) GRIDC
      REAL(q) QTOT(T_INFO%NIONS)
      RGRID POTTOT(DIMREAL(GRIDC%MPLWV))
      REAL(q) IFOR(3)
      ! local variables
      INTEGER I,K
      INTEGER N1LOW,N2LOW,N3LOW,N1HI,N2HI,N3HI
      INTEGER N1,N2,N3,N1P,N2P,N3P,NCOL,IND
      REAL(q) F1,F2,F3,D1,D2,D3,X1,X2,X3
      REAL(q) XC,YC,ZC,D,DX
      REAL(q) VAL,VALDER,RCUT
      REAL(q) :: DHDER(3)
      COMPLEX(q) CIFOR(3)

      F1=1._q/GRIDC%NGX
      F2=1._q/GRIDC%NGY
      F3=1._q/GRIDC%NGZ

      !-----------------------------------------------------------------------
      ! restrict loop to points contained within a cubus around the ion
      !-----------------------------------------------------------------------
      D1= P%USECUT*LATT_CUR%BNORM(1)*GRIDC%NGX
      D2= P%USECUT*LATT_CUR%BNORM(2)*GRIDC%NGY
      D3= P%USECUT*LATT_CUR%BNORM(3)*GRIDC%NGZ

      N3LOW= INT(T_INFO%POSION(3,NI)*GRIDC%NGZ-D3+10*GRIDC%NGZ+.99_q)-10*GRIDC%NGZ
      N2LOW= INT(T_INFO%POSION(2,NI)*GRIDC%NGY-D2+10*GRIDC%NGY+.99_q)-10*GRIDC%NGY
      N1LOW= INT(T_INFO%POSION(1,NI)*GRIDC%NGX-D1+10*GRIDC%NGX+.99_q)-10*GRIDC%NGX

      N3HI = INT(T_INFO%POSION(3,NI)*GRIDC%NGZ+D3+10*GRIDC%NGZ)-10*GRIDC%NGZ
      N2HI = INT(T_INFO%POSION(2,NI)*GRIDC%NGY+D2+10*GRIDC%NGY)-10*GRIDC%NGY
      N1HI = INT(T_INFO%POSION(1,NI)*GRIDC%NGX+D1+10*GRIDC%NGX)-10*GRIDC%NGX

      IND=0
      DHDER=0
      CIFOR=0
      !-----------------------------------------------------------------------
      ! loop over cubus
      ! MPI version z ist the fast index
      !-----------------------------------------------------------------------
#ifdef MPI
      DO N2=N2LOW,N2HI
         X2=(N2*F2-T_INFO%POSION(2,NI))
         N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)

         DO N1=N1LOW,N1HI
            X1=(N1*F1-T_INFO%POSION(1,NI))
            N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)

            NCOL=GRIDC%RL%INDEX(N1P,N2P)
            IF (NCOL==0) CYCLE ! not on local node, move on
            IF (GRIDC%RL%I2(NCOL) /= N1P+1 .OR. GRIDC%RL%I3(NCOL) /= N2P+1) THEN
               CALL vtutor%bug("RHOADD: internal ERROR: " // str(GRIDC%RL%I2(NCOL)) // " " // &
                  str(N1P+1) // " " // str(GRIDC%RL%I3(NCOL)) // " " // str(N2P+1), __FILE__, __LINE__)
            ENDIF
            !OCL SCALAR
            DO N3=N3LOW,N3HI
               X3=(N3*F3-T_INFO%POSION(3,NI))
               N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
               IND=N3P+(NCOL-1)*GRIDC%NGZ+1
#else
      !-----------------------------------------------------------------------
      ! loop over cubus around one ion
      ! conventional version x is fast index
      !-----------------------------------------------------------------------
      DO N3=N3LOW,N3HI
         X3=(N3*F3-T_INFO%POSION(3,NI))
         N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
         
         DO N2=N2LOW,N2HI
            X2=(N2*F2-T_INFO%POSION(2,NI))
            N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
            
            NCOL=GRIDC%RL%INDEX(N2P,N3P)
            
            DO N1=N1LOW,N1HI
               X1=(N1*F1-T_INFO%POSION(1,NI))
               N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
               IND=N1P+(NCOL-1)*GRIDC%NGX+1
#endif        
               
               XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
               YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
               ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
               
               D=SQRT(XC*XC+YC*YC+ZC*ZC)
               
               IF (D<=P%USECUT.AND.D>0) THEN
                  ! derivative DHDER(a) = d rho/d R_NI,a

                  K=D/P%USESPL(2,1)+1
                  DX=D-P%USESPL(K,1)
                  VALDER=(3.0_q*P%USESPL(K,5)*DX+2.0_q*P%USESPL(K,4))*DX+P%USESPL(K,3)

                  VAL=VALDER/QTOT(NI)*P%USEZ*GRIDC%NPLWV/D
                  DHDER(1)=VAL*XC
                  DHDER(2)=VAL*YC
                  DHDER(3)=VAL*ZC
                  ! force CIFOR(a) = \int d rho/d R_NI,a (r) V(r)
                  CIFOR=CIFOR+DHDER*POTTOT(IND)
               ENDIF
               
            ENDDO
         ENDDO
      ENDDO

      CALLMPI( M_sum_z(GRIDC%COMM, CIFOR, 3) )
      CIFOR=CIFOR/LATT_CUR%OMEGA*F1*F2*F3

      IFOR=REAL(CIFOR,kind=q)

      RETURN
      END SUBROUTINE RHODER


!************************* FUNCTION RHO0 *******************************
!
!> Calculates the total number of electrons in recip space
!> (not quite trival in parallel mode)
!
!***********************************************************************

      FUNCTION RHO0(GRID, CHTOT)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)

      N2= 1
      N3= 1
      N1= 1

      RHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          RHO_SUM=CHTOT(NC,N1)
        ENDIF
      ENDDO

      CALLMPI( M_sum_d( GRID%COMM, RHO_SUM, 1))
      RHO0= RHO_SUM
      RETURN
      END FUNCTION


      FUNCTION CRHO0(GRID, CHTOT)
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) CRHO0
    ! local
      INTEGER N1, N2, N3, NC
      COMPLEX(q) CRHO_SUM

      N2= 1
      N3= 1
      N1= 1

      CRHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          CRHO_SUM=CHTOT(NC,N1)
        ENDIF
      ENDDO

      CALLMPI( M_sum_z( GRID%COMM, CRHO_SUM, 1))
      CRHO0= CRHO_SUM
      RETURN
      END FUNCTION


!************************* SUBROUTINE GET_RHO0 *************************
!
!> Calculates the total number of electrons in recip space
!> (not quite trival in parallel mode)
!
!***********************************************************************

      SUBROUTINE GET_RHO0(GRID, CHTOT, RHO0)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)
      N2= 1
      N3= 1
      N1= 1
      RHO=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          RHO=CHTOT(NC,N1)
        ENDIF
      ENDDO

      CALLMPI( M_sum_d( GRID%COMM, RHO, 1))
      RHO0= RHO
      RETURN
      END SUBROUTINE


!************************* SUBROUTINE SET_RHO0 *************************
!
!> Sets the total number of electrons in recip space
!> (not quite trival in parallel mode)
!
!***********************************************************************

      SUBROUTINE SET_RHO0(GRID, CHTOT, RHO_SOLL)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)

      N2= 1
      N3= 1
      N1= 1

      RHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          CHTOT(NC,N1)=RHO_SOLL
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE


      MODULE compat_gga
!> GGA_COMPAT = .TRUE. the vasp.4.4-4.6 behavior is used\n
!> GGA_COMPAT = .FALSE. the new corrected version is used (see #TRUNC_HIGH_FREQU)
!>
!> This variable is set by #GGA_COMPAT_MODE (default is GGA_COMPAT=.TRUE.)
        LOGICAL GGA_COMPAT
      END MODULE compat_gga


!************************* SUBROUTINE GGA_COMPAT_MODE ******************
!
!> Reads GGA_COMPAT (compat_gga::GGA_COMPAT) from the INCAR file
!>
!> GGA_COMPAT = .TRUE. is the default
!
!***********************************************************************

      SUBROUTINE GGA_COMPAT_MODE(IU5, IU0, LCOMPAT)
        USE prec
        USE base
        USE compat_gga
        USE reader_tags
        IMPLICIT NONE 
        
        INTEGER :: IU5,IU0    !< io units
        LOGICAL :: LCOMPAT

        LOGICAL :: LOPEN,LDUM
        INTEGER :: IDUM, N, IERR
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER (LEN=1) :: CHARAC
        
        CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

        GGA_COMPAT=.TRUE.
        CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'GGA_COMPAT', GGA_COMPAT, IERR, WRITEXMLINCAR)

        CALL CLOSE_INCAR_IF_FOUND(IU5)
      END SUBROUTINE GGA_COMPAT_MODE

      SUBROUTINE XML_WRITE_GGA_COMPAT_MODE
        USE prec
        USE compat_gga
        IMPLICIT NONE
        INTEGER :: IDUM
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER (LEN=1) :: CHARAC

        CALL XML_INCAR('GGA_COMPAT','L',IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,1)

      END SUBROUTINE XML_WRITE_GGA_COMPAT_MODE

!************************* SUBROUTINE GRADIENT_CHG *********************
!> @brief Compute the gradient of CHTOT in Cartesian direction x, y, z
!
!> The gradient is computed in Cartesian basis by applying the
!> G vector in reciprocal space. By default no FFT is performed but
!> optional arguments LRECIN and LRECOUT can be used to switch between
!> reciprocal and real space.
!
!> The input CHTOT can be provided in
!>    reciprocal space (LRECIN=.TRUE., default)
!>    or real space (LRECIN=.FALSE.).
!
!> The components of the gradient (CGRADX, CGRADY, CGRADZ) are returned
!> in
!>    reciprocal space (LRECOUT=.TRUE., default)
!>    or real space (LRECOUT=.FALSE.).
!***********************************************************************
      SUBROUTINE GRADIENT_CHG(LATT_CUR,GRIDC,CHTOT,CGRADX,CGRADY,CGRADZ,LRECIN,LRECOUT)

__GPU USE moffload

        USE mpimy
        USE poscar, ONLY: latt
        USE mgrid, ONLY: grid_3d
        USE constant
        IMPLICIT NONE
        TYPE(latt), INTENT(IN)       :: LATT_CUR              !< lattice parameter
        TYPE(grid_3d), INTENT(IN)    :: GRIDC                 !< grid for quantity
        COMPLEX(q), INTENT(IN)       :: CHTOT(GRIDC%MPLWV)    !< grid quantity
        COMPLEX(q), INTENT(OUT)      :: CGRADX(GRIDC%MPLWV)             !< gradient along Cartesian direction x
        COMPLEX(q), INTENT(OUT)      :: CGRADY(GRIDC%MPLWV)             !< gradient along Cartesian direction y
        COMPLEX(q), INTENT(OUT)      :: CGRADZ(GRIDC%MPLWV)             !< gradient along Cartesian direction z
        ! --- optional
        LOGICAL, INTENT(IN), OPTIONAL :: LRECIN     !< grid quantity is provided in reciprocal space
        LOGICAL, INTENT(IN), OPTIONAL :: LRECOUT    !< gradien is returned in reciprocal space
        ! --- local
        INTEGER :: I, N1, NC, N2, N3
        REAL(q) :: GX, GY, GZ, RINPL
        COMPLEX(q), ALLOCATABLE :: CWORK(:)
        INTEGER :: NODE_ME, IONODE, IDUMP
        ! --- initialize defaults
        LOGICAL :: LRECIN_
        LOGICAL :: LRECOUT_

        LRECIN_=.TRUE.     ; IF (PRESENT(LRECIN)) LRECIN_=LRECIN
        LRECOUT_=.TRUE.    ; IF (PRESENT(LRECOUT)) LRECOUT_=LRECOUT

! set to one for error-dumps
#ifdef MPI
        NODE_ME=GRIDC%COMM%NODE_ME
        IONODE =GRIDC%COMM%IONODE
        IDUMP=0
#ifdef debug
        IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
        NODE_ME=0
        IONODE =0
        IDUMP=0
#ifdef debug
        IDUMP=1
#endif
#endif

        ALLOCATE(CWORK(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: CWORK) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CWORK,CHTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK,CHTOT) IF(OFFLOAD_ON)
        CWORK=CHTOT
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        IF (IDUMP/=0) WRITE(*,*) "GRADIENT_CHG CHTOT real: CWORK(1:3)", CWORK(1:3)

        IF (.NOT.LRECIN_) THEN
          ! to rec space
          RINPL=1._q/GRIDC%NPLWV
          CWORK=CWORK*RINPL
          CALL FFT3D(CWORK(1),GRIDC,-1)
          CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1))
        ENDIF

        IF (IDUMP/=0) WRITE(*,*) "GRADIENT_CHG CHTOT rec: CWORK(1:3)", CWORK(1:3)

! ----------------
! gradient along x
! ----------------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, NC, N2, N3, GX) SHARED(GRIDC,LATT_CUR,CWORK,CGRADX)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CGRADX,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CGRADX,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          CGRADX(I)=CWORK(I)*GX*CITPI
        ENDDO
__CPU !$OMP END PARALLEL DO
        IF (IDUMP/=0) WRITE(*,*) "GRADIENT_CHG grad x rec: CGRADX(1:3)", CGRADX(1:3)

! ----------------
! gradient along y
! ----------------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, NC, N2, N3, GY) SHARED(GRIDC,LATT_CUR,CWORK,CGRADY)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CGRADY,LATT_CUR) PRIVATE(N1,N2,N3,NC,GY) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CGRADY,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GY) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          CGRADY(I)=CWORK(I)*GY*CITPI
        ENDDO
__CPU !$OMP END PARALLEL DO

! ----------------
! gradient along z
! ----------------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, NC, N2, N3, GZ) SHARED(GRIDC,LATT_CUR,CWORK,CGRADZ)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CGRADZ,LATT_CUR) PRIVATE(N1,N2,N3,NC,GZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CGRADZ,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          CGRADZ(I)=CWORK(I)*GZ*CITPI
        ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK) IF(OFFLOAD_ON)
        DEALLOCATE(CWORK)

        IF (.NOT.LRECOUT_) THEN
          ! to real space
          CALL SETUNB(CGRADX(1),GRIDC)
          CALL SETUNB(CGRADY(1),GRIDC)
          CALL SETUNB(CGRADZ(1),GRIDC)
          CALL FFT3D(CGRADX(1),GRIDC,1)
          CALL FFT3D(CGRADY(1),GRIDC,1)
          CALL FFT3D(CGRADZ(1),GRIDC,1)
        ENDIF
        IF (IDUMP/=0) WRITE(*,*) "GRADIENT_CHG grad x real: CGRADX(1:3)", CGRADX(1:3)

      END SUBROUTINE GRADIENT_CHG

!************************* SUBROUTINE LAPLACIAN_CHG *********************
!> @brief Compute the Laplacian of CHTOT
!***********************************************************************
      SUBROUTINE LAPLACIAN_CHG(LATT_CUR,GRIDC,CHTOT,CLAP,LRECIN,LRECOUT)

__GPU USE moffload

        USE mpimy
        USE poscar, ONLY: latt
        USE mgrid, ONLY: grid_3d
        USE constant
        IMPLICIT NONE
        TYPE(latt), INTENT(IN)       :: LATT_CUR              !< lattice parameter
        TYPE(grid_3d), INTENT(IN)    :: GRIDC                 !< grid for quantity
        COMPLEX(q), INTENT(IN)       :: CHTOT(GRIDC%MPLWV)    !< grid quantity
        COMPLEX(q), INTENT(OUT)      :: CLAP(GRIDC%MPLWV)     !< Laplacian
        ! --- optional
        LOGICAL, INTENT(IN), OPTIONAL :: LRECIN     !< grid quantity is provided in reciprocal space
        LOGICAL, INTENT(IN), OPTIONAL :: LRECOUT    !< Laplacian is returned in reciprocal space
        ! --- local
        INTEGER :: I, N1, NC, N2, N3
        REAL(q) :: GX, GY, GZ, RINPL
        COMPLEX(q), ALLOCATABLE :: CWORK(:)
        INTEGER :: NODE_ME, IONODE, IDUMP
        ! --- initialize defaults
        LOGICAL :: LRECIN_
        LOGICAL :: LRECOUT_

        LRECIN_=.TRUE.     ; IF (PRESENT(LRECIN)) LRECIN_=LRECIN
        LRECOUT_=.TRUE.    ; IF (PRESENT(LRECOUT)) LRECOUT_=LRECOUT
!#define debug
! set to one for error-dumps
#ifdef MPI
        NODE_ME=GRIDC%COMM%NODE_ME
        IONODE =GRIDC%COMM%IONODE
        IDUMP=0
#ifdef debug
        IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
        NODE_ME=0
        IONODE =0
        IDUMP=0
#ifdef debug
        IDUMP=1
#endif
#endif

        ALLOCATE(CWORK(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: CWORK) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CWORK,CHTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK,CHTOT) IF(OFFLOAD_ON)
        CWORK=CHTOT
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        IF (IDUMP/=0) WRITE(*,*) "LAP_CHG CHTOT in: CWORK(1:3)", CWORK(1:3)

        IF (.NOT.LRECIN_) THEN
          ! to rec space
          RINPL=1._q/GRIDC%NPLWV
          CWORK=CWORK*RINPL
          CALL FFT3D(CWORK(1),GRIDC,-1)
          CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1))
        ENDIF

        IF (IDUMP/=0) WRITE(*,*) "HESSIAN_CHG CHTOT rec: CWORK(1:3)", CWORK(1:3)

        CLAP=(0.0_q,0.0_q)

!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GX, GY, GZ) SHARED(GRIDC,LATT_CUR,CWORK,CLAP)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CLAP,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GY,GZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CLAP,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX,GY,GZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          CLAP(I)=CLAP(I)-(CWORK(I)*GX*GX+CWORK(I)*GY*GY+CWORK(I)*GZ*GZ)*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK) IF(OFFLOAD_ON)
        DEALLOCATE(CWORK)

        IF (.NOT.LRECOUT_) THEN
          IF (IDUMP/=0) WRITE(*,*) "LAP_CHG rec: CLAP(1:3)", CLAP(1:3)
          ! to real space
          CALL SETUNB(CLAP(1),GRIDC)
          CALL FFT3D(CLAP(1),GRIDC,1)
          IF (IDUMP/=0) WRITE(*,*) "LAP_CHG real: CLAP(1:3)", CLAP(1:3)
        ENDIF

      END SUBROUTINE LAPLACIAN_CHG

!************************* SUBROUTINE HESSIAN_CHG *********************
!> @brief Compute the Hessian of CHTOT in Cartesian direction
!>        xx, yy, zz, xy, yz, zx
!***********************************************************************
      SUBROUTINE HESSIAN_CHG(LATT_CUR,GRIDC,CHTOT,CHESS,LRECIN,LRECOUT)

__GPU USE moffload

        USE mpimy
        USE poscar, ONLY: latt
        USE mgrid, ONLY: grid_3d
        USE constant
        IMPLICIT NONE
        TYPE(latt), INTENT(IN)       :: LATT_CUR              !< lattice parameter
        TYPE(grid_3d), INTENT(IN)    :: GRIDC                 !< grid for quantity
        COMPLEX(q), INTENT(IN)       :: CHTOT(GRIDC%MPLWV)    !< grid quantity
        COMPLEX(q), INTENT(OUT)      :: CHESS(GRIDC%MPLWV,6)  !< Hessian along xx, xy, xz, yy, yz, zz
        ! --- optional
        LOGICAL, INTENT(IN), OPTIONAL :: LRECIN     !< grid quantity is provided in reciprocal space
        LOGICAL, INTENT(IN), OPTIONAL :: LRECOUT    !< gradien is returned in reciprocal space
        ! --- local
        INTEGER :: I, N1, NC, N2, N3, IDH
        REAL(q) :: GX, GY, GZ, RINPL
        COMPLEX(q), ALLOCATABLE :: CWORK(:)
        INTEGER :: NODE_ME, IONODE, IDUMP
        ! --- initialize defaults
        LOGICAL :: LRECIN_
        LOGICAL :: LRECOUT_

        LRECIN_=.TRUE.     ; IF (PRESENT(LRECIN)) LRECIN_=LRECIN
        LRECOUT_=.TRUE.    ; IF (PRESENT(LRECOUT)) LRECOUT_=LRECOUT
!#define debug
! set to one for error-dumps
#ifdef MPI
        NODE_ME=GRIDC%COMM%NODE_ME
        IONODE =GRIDC%COMM%IONODE
        IDUMP=0
#ifdef debug
        IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
        NODE_ME=0
        IONODE =0
        IDUMP=0
#ifdef debug
        IDUMP=1
#endif
#endif

        ALLOCATE(CWORK(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: CWORK) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CWORK,CHTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK,CHTOT) IF(OFFLOAD_ON)
        CWORK=CHTOT
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

        IF (IDUMP/=0) WRITE(*,*) "HESSIAN_CHG CHTOT in: CWORK(1:3)", CWORK(1:3)

        IF (.NOT.LRECIN_) THEN
          ! to rec space
          RINPL=1._q/GRIDC%NPLWV
          CWORK=CWORK*RINPL
          CALL FFT3D(CWORK(1),GRIDC,-1)
          CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1))
        ENDIF

        IF (IDUMP/=0) WRITE(*,*) "HESSIAN_CHG CHTOT rec: CWORK(1:3)", CWORK(1:3)
! ----------
! Hessian xx
! ----------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GX) SHARED(GRIDC,LATT_CUR,CWORK,CHESS)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CHESS,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CHESS,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          CHESS(I,1)=-CWORK(I)*GX*GX*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

! ----------
! Hessian xy
! ----------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GX, GY) SHARED(GRIDC,LATT_CUR,CWORK,CHESS)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CHESS,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GY) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CHESS,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX,GY) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          CHESS(I,4)=-CWORK(I)*GX*GY*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

! ----------
! Hessian xz
! ----------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GX, GZ) SHARED(GRIDC,LATT_CUR,CWORK,CHESS)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CHESS,LATT_CUR) PRIVATE(N1,N2,N3,NC,GX,GZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CHESS,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX,GZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          CHESS(I,6)=-CWORK(I)*GX*GZ*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

! ----------
! Hessian yy
! ----------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GY) SHARED(GRIDC,LATT_CUR,CWORK,CHESS)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CHESS,LATT_CUR) PRIVATE(N1,N2,N3,NC,GY) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CHESS,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GY) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          CHESS(I,2)=-CWORK(I)*GY*GY*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

! ----------
! Hessian yz
! ----------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GY, GZ) SHARED(GRIDC,LATT_CUR,CWORK,CHESS)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CHESS,LATT_CUR) PRIVATE(N1,N2,N3,NC,GY,GZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CHESS,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GY,GZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          CHESS(I,5)=-CWORK(I)*GY*GZ*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

! ----------
! Hessian zz
! ----------
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(I, N1, N2, N3, NC, GZ) SHARED(GRIDC,LATT_CUR,CWORK,CHESS)

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,CWORK,CHESS,LATT_CUR) PRIVATE(N1,N2,N3,NC,GZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,CWORK,CHESS,LATT_CUR) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          N1= MOD((I-1),GRIDC%RC%NROW)+1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          CHESS(I,3)=-CWORK(I)*GZ*GZ*FPI2
        ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK) IF(OFFLOAD_ON)
        DEALLOCATE(CWORK)

        IF (.NOT.LRECOUT_) THEN
          DO IDH=1,6
            IF (IDUMP/=0) WRITE(*,*) "HESSIAN_CHG IDH rec: CHESS(1:3,IDH)", IDH, CHESS(1:3,IDH)
            ! to real space
            CALL SETUNB(CHESS(1,IDH),GRIDC)
            CALL FFT3D(CHESS(1,IDH),GRIDC,1)
            IF (IDUMP/=0) WRITE(*,*) "HESSIAN_CHG IDH real: CHESS(1:3,IDH)", IDH, CHESS(1:3,IDH)
          ENDDO
        ENDIF

      END SUBROUTINE HESSIAN_CHG

      SUBROUTINE CALC_TOTAL_CHARGE(LATT_CUR,GRIDC,NCDIJ,CHTOT,DENCOR,CRHO,LRECIN,LRECOUT,LREP_TOTALMAG)

__GPU USE moffload

        USE base, only : CAST_COMPLEX_TO_RGRID
        USE lattice
        USE mgrid
        USE mpimy
        IMPLICIT NONE
        TYPE (latt), INTENT(IN)        :: LATT_CUR
        TYPE (grid_3d), INTENT(IN)     :: GRIDC
        INTEGER, INTENT(IN)            :: NCDIJ
        COMPLEX(q), TARGET, INTENT(IN) :: CHTOT(GRIDC%MPLWV,NCDIJ) ! charge density in (total,mag)
        RGRID, INTENT(IN)              :: DENCOR(GRIDC%RL%NP)
        COMPLEX(q), TARGET, INTENT(OUT):: CRHO(GRIDC%MPLWV,NCDIJ)
        ! optional
        LOGICAL, INTENT(IN), OPTIONAL :: LRECIN     !< CHTOT is provided in reciprocal space
        LOGICAL, INTENT(IN), OPTIONAL :: LRECOUT    !< RHO is returned in reciprocal space
        LOGICAL, INTENT(IN), OPTIONAL :: LREP_TOTALMAG !< spin-polarized case is in (total,mag) representation
        ! local
        INTEGER ISP, I
        RGRID, POINTER ::  DWORK(:,:) => NULL()
        REAL(q) :: OOOMEGA
        INTEGER :: NODE_ME, IONODE, IDUMP
        ! initialize defaults
        LOGICAL LRECIN_, LRECOUT_, LREP_TOTALMAG_

        LRECIN_=.FALSE.     ; IF (PRESENT(LRECIN)) LRECIN_=LRECIN
        LRECOUT_=.FALSE.    ; IF (PRESENT(LRECOUT)) LRECOUT_=LRECOUT
        LREP_TOTALMAG_=.TRUE.; IF (PRESENT(LREP_TOTALMAG)) LREP_TOTALMAG_=LREP_TOTALMAG

! set to one for error-dumps
#ifdef MPI
        NODE_ME=GRIDC%COMM%NODE_ME
        IONODE =GRIDC%COMM%IONODE
        IDUMP=0
#ifdef debug
        IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
        NODE_ME=0
        IONODE =0
        IDUMP=0
#ifdef debug
        IDUMP=1
#endif
#endif

        CALL CAST_COMPLEX_TO_RGRID(CRHO,DWORK)

!------ add partial core charges and grid charges
        IF (LRECOUT_) THEN
          OOOMEGA=1.0_q/LATT_CUR%OMEGA/GRIDC%NPLWV
        ELSE
          OOOMEGA=1.0_q/LATT_CUR%OMEGA
        ENDIF
GPACC !$ACC ENTER DATA COPYIN(OOOMEGA) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:OOOMEGA) IF(OFFLOAD_ON)

        spin: DO ISP=1,NCDIJ
          CALL RL_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CRHO(1,ISP),GRIDC)
          IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: CHTOT in:", CRHO(1:3,ISP)
          IF (LRECIN_) THEN
!---------- go to real space
            CALL SETUNB(CRHO(1,ISP),GRIDC)
            CALL FFT3D(CRHO(1,ISP),GRIDC,1)
          ENDIF
          IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: CHTOT real:", DWORK(1:6,ISP)
          IF (NCDIJ==1) THEN
!---------- nonmagnetic case
            IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: DENCOR", DENCOR(1:6)
            IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: OOOMEGA", OOOMEGA
            CALL RL_ADD(DWORK(1,ISP),OOOMEGA,DENCOR(1),OOOMEGA,DWORK(1,ISP),GRIDC)
            IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: DWORK real", DWORK(1:6,ISP)
          ELSEIF (NCDIJ==2) THEN
!---------- spin-polarized case
            ! (up,down) representation
            IF (LREP_TOTALMAG_) CALL RL_FLIP(DWORK, GRIDC, 2, .TRUE.)
            CALL RL_ADD(DWORK(1,ISP),OOOMEGA,DENCOR(1),OOOMEGA/NCDIJ,DWORK(1,ISP),GRIDC)
          ELSEIF (NCDIJ==4) THEN
!---------- noncollinear case
            IF (ISP==1) THEN
              CALL RL_ADD(DWORK(1,ISP),OOOMEGA,DENCOR(1),OOOMEGA,DWORK(1,ISP),GRIDC)
            ELSE
              CALL RL_ADD(DWORK(1,ISP),OOOMEGA,DWORK(1,ISP),0.0_q,DWORK(1,ISP),GRIDC)
            ENDIF
          ENDIF
          CALL OPSYNC(DWORK(1,ISP),CRHO(1,ISP),GRIDC%NPLWV)
          IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: CRHO real:", CRHO(1:3,ISP)

          IF (LRECOUT_) THEN
!---------- go to rec space
            CALL FFT3D(CRHO(1,ISP),GRIDC,-1)
            IF (IDUMP/=0) WRITE(*,*) "CALC_TOTAL_CHARGE: CRHO rec:", CRHO(1:3,ISP)
          ENDIF
        ENDDO spin

GPACC !$ACC EXIT DATA DELETE(OOOMEGA) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:OOOMEGA) IF(OFFLOAD_ON)

      END SUBROUTINE CALC_TOTAL_CHARGE

!************************* SUBROUTINE TRUNC_HIGH_FREQU *****************
!
!> Truncates the high frequency components of an array by imposing
!> a spherical cutoff
!>
!> The spherical cutoff is chosen in such a way that the entire sphere
!> fits into the parallelepided spanned by the plane wave basis set
!>
!> This subroutine is required for GGA calculations to maintain the
!> full symmetry of the exchange correlation potential
!
!***********************************************************************

      SUBROUTINE TRUNC_HIGH_FREQU(LATT_CUR, GRID, C)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid
      USE lattice
      USE compat_gga
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt)    LATT_CUR
      COMPLEX(q)     C(GRID%MPLWV)
      REAL(q) GX, GY, GZ, G2, GMIN2
      INTEGER I, N1, N2, N3, NC, NZERO
!
! if you want to revert to the behaviour of vasp before vasp.4.6.16
! comment in the RETURN statment below
!
      IF (GGA_COMPAT) RETURN

      GX=GRID%NGX/2/LATT_CUR%ANORM(1)
      GY=GRID%NGY/2/LATT_CUR%ANORM(2)
      GZ=GRID%NGZ/2/LATT_CUR%ANORM(3)

      GMIN2=MIN(GX,GY,GZ)**2

GPACC !$ACC PARALLEL LOOP PRESENT(GRID,LATT_CUR,C) PRIVATE(N1,N2,N3,NC,GX,GY,GZ,G2) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRID,LATT_CUR,C) &
GPOMP !$OMP PRIVATE(N1,N2,N3,NC,GX,GY,GZ,G2) IF(OFFLOAD_ON)
      DO I=1,GRID%RC%NP
         N1= MOD((I-1),GRID%RC%NROW) +1
         NC= (I-1)/GRID%RC%NROW+1
         N2= GRID%RC%I2(NC)
         N3= GRID%RC%I3(NC)
         GX=(GRID%LPCTX(N1)*LATT_CUR%B(1,1)+GRID%LPCTY(N2)*LATT_CUR%B(1,2)+GRID%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRID%LPCTX(N1)*LATT_CUR%B(2,1)+GRID%LPCTY(N2)*LATT_CUR%B(2,2)+GRID%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRID%LPCTX(N1)*LATT_CUR%B(3,1)+GRID%LPCTY(N2)*LATT_CUR%B(3,2)+GRID%LPCTZ(N3)*LATT_CUR%B(3,3))
         G2=GX*GX+GY*GY+GZ*GZ

         IF (G2>GMIN2) THEN
            C(I)=0
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE
