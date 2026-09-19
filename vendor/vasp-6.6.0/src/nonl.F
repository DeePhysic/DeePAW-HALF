!#define debug
#include "symbol.inc"
!***********************************************************************
!
!> This module contains all the routines required to support
!> reciprocal space presentation of the non local projection operators
!> on a plane wave grid
!>
!> @note
!> The wavefunction character is strictly defined as
!> ~~~
!>  C_n = e (i k R) \sum_G p_n,G+k e (i G R) C_G
!> ~~~
!> VASP however does not include the term <TT>e (i k R)</TT>, since
!> it does not change the occupancy matrix or energy. Hence
!> ~~~
!>  W%CPROJ(n) = \sum_G p_n,G+k e (i G R) C_G
!> ~~~
!> the derivative of C_n with respect to R can either include the
!> term from e (i k R)
!> ~~~
!>    i k \sum_G  p_n,G+k e (i G R) C_G = i k C_n
!> ~~~
!> or not.
!> This contribution does not contribute to the one center occupancy
!> matrix or the total energy either, since it is "complex"
!> (expectation values involving C_m* D_mn C_n are complex)
!> the Berry phase seems to require <TT>-Dshift_der_k</TT> however
!> (nonlr.F uses a similar convention)
!***********************************************************************
#define shift_der_k

MODULE nonl

__GPU USE moffload

  USE prec
  USE wave
  USE nonl_struct_def

  INTERFACE
     SUBROUTINE VNLAC0(NONL_S,WDES1,CPROJ_LOC,CACC)
       USE nonl_struct_def
       USE wave_struct_def
       TYPE (nonl_struct) NONL_S
       TYPE (wavedes1)    WDES1
       GDEF    CPROJ_LOC
       COMPLEX(q) CACC
     END SUBROUTINE VNLAC0
  END INTERFACE

  INTERFACE
     SUBROUTINE VNLAC0_DER(NONL_S,WDES1,CPROJ_LOC,CACC, LATT_CUR, ION, IDIR)
       USE nonl_struct_def
       USE wave_struct_def
       USE lattice
       TYPE (nonl_struct) NONL_S
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CACC
       GDEF         CPROJ_LOC
       TYPE (latt)  LATT_CUR
       INTEGER      ION
       INTEGER      IDIR
     END SUBROUTINE VNLAC0_DER
  END INTERFACE

CONTAINS

!****************** SUBROUTINE NONL_ALLOC ******************************
! RCS:  $Id: nonl.F,v 1.2 2002/08/14 13:59:41 kresse Exp $
!
!> Allocates required arrays based on T_INFO and P structure,
!>
!> The number of ions and types is taken from T_INFO and
!> LMDIM is taken from the P structure
!
!***********************************************************************

  SUBROUTINE NONL_ALLOC(NONL_S,T_INFO,P,WDES, LREAL)
    USE prec
    USE pseudo
    USE poscar
    USE ini
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    LOGICAL LREAL

    ! local var
    INTEGER NIONS,NTYPD,LMDIM,NRPLWV,NKPTS,NT,NIS,NI

    NIONS =T_INFO%NIONS
    NTYPD =T_INFO%NTYPD
    LMDIM =P(1)%LMDIM
    NRPLWV=WDES%NGDIM
    NKPTS =WDES%NKPTS   ! number of k-points in the IBZ
    ! (required for HF type calculations)
    NONL_S%LRECIP=.NOT. LREAL
    NONL_S%NTYP  =T_INFO%NTYP
    NONL_S%NK    =0
    NONL_S%NIONS =T_INFO%NIONS
    NONL_S%NITYP =>T_INFO%NITYP
    NONL_S%ITYP  =>T_INFO%ITYP
    NONL_S%POSION=>T_INFO%POSION
    NONL_S%LSPIRAL=WDES%LSPIRAL
    NULLIFY(NONL_S%VKPT_SHIFT)

    ALLOCATE(NONL_S%LMMAX(NONL_S%NTYP),NONL_S%LMBASE(NONL_S%NIONS+1))

    NONL_S%LMBASE(1)=0; NIS=1
    DO NT=1,NONL_S%NTYP
       NONL_S%LMMAX(NT)=P(NT)%LMMAX
       DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
          NONL_S%LMBASE(NI+1)=NONL_S%LMBASE(NI)+NONL_S%LMMAX(NT)
       ENDDO
       NIS=NIS+NONL_S%NITYP(NT)
    ENDDO

    NULLIFY(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)

    IF ((.NOT.LREAL).AND.(.NOT.NONL_S%LSPIRAL)) THEN
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
         NONL_S%CQFAK(LMDIM,NTYPD))
         CALL REGISTER_ALLOCATE(16._q*SIZE(NONL_S%CREXP)+8._q*SIZE(NONL_S%QPROJ), "nonl-proj")
    ENDIF
    IF ((.NOT.LREAL).AND. NONL_S%LSPIRAL) THEN
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
         NONL_S%CQFAK(LMDIM,NTYPD))
         CALL REGISTER_ALLOCATE(16._q*SIZE(NONL_S%CREXP)+8._q*SIZE(NONL_S%QPROJ), "nonl-proj")
    ENDIF

    RETURN
  END SUBROUTINE NONL_ALLOC


  SUBROUTINE NONL_DEALLOC(NONL_S)
    USE ini
    IMPLICIT NONE
    TYPE (nonl_struct) NONL_S
    LOGICAL LREAL

    LREAL=.NOT.NONL_S%LRECIP

    IF (.NOT.LREAL) THEN
       CALL DEREGISTER_ALLOCATE(16._q*SIZE(NONL_S%CREXP)+8._q*SIZE(NONL_S%QPROJ), "nonl-proj")
       DEALLOCATE(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)
       DEALLOCATE(NONL_S%LMMAX)
       DEALLOCATE(NONL_S%LMBASE)
    ENDIF

    RETURN
  END SUBROUTINE NONL_DEALLOC

!****************** SUBROUTINE NONL_ALLOC_SPHPRO ***********************
!
!> allocates required arrays to describe projectors for one ion
!> and one k-point
!
!***********************************************************************

  SUBROUTINE NONL_ALLOC_SPHPRO(NONL_S,P,WDES)
    USE pseudo
    USE poscar
    IMPLICIT NONE


    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P
    TYPE (wavedes)     WDES
    ! local var
    INTEGER NIONS,NTYPD,LMDIM,NRPLWV,NKPTS,NT,NIS,NI

    NIONS =1
    NTYPD =1
    LMDIM =P%LMDIM
    NRPLWV=WDES%NGDIM
    NKPTS =WDES%NKPTS

    NONL_S%LRECIP=.TRUE.
    NONL_S%NTYP  =NTYPD
    NONL_S%NK    =0
    NONL_S%NIONS =NIONS
    NONL_S%LSPIRAL=WDES%LSPIRAL
    ALLOCATE(NONL_S%NITYP(NTYPD)); NONL_S%NITYP(1)=1
    ALLOCATE(NONL_S%ITYP(NTYPD)) ; NONL_S%ITYP(1)=1
    ALLOCATE(NONL_S%LMMAX(NTYPD)); NONL_S%LMMAX(1)=P%LMMAX

    ALLOCATE(NONL_S%LMBASE(NIONS+1))
    NONL_S%LMBASE(1)=0; NIS=1
    DO NT=1,NONL_S%NTYP
       DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
          NONL_S%LMBASE(NI+1)=NONL_S%LMBASE(NI)+NONL_S%LMMAX(NT)
       ENDDO
       NIS=NIS+NONL_S%NITYP(NT)
    ENDDO

    IF (.NOT.NONL_S%LSPIRAL) &
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
         NONL_S%CQFAK(LMDIM,NTYPD))
    IF (NONL_S%LSPIRAL) &
         ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
         NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
         NONL_S%CQFAK(LMDIM,NTYPD))
    RETURN
  END SUBROUTINE NONL_ALLOC_SPHPRO

!****************** SUBROUTINE NONL_DEALLOC_SPHPRO *********************
!
!> Deallocate required arrays to describe projectors for one ion
!> and one k-point
!
!***********************************************************************

  SUBROUTINE NONL_DEALLOC_SPHPRO(NONL_S)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    DEALLOCATE(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)
    DEALLOCATE(NONL_S%NITYP,NONL_S%ITYP)
    DEALLOCATE(NONL_S%LMMAX,NONL_S%LMBASE)

    RETURN
  END SUBROUTINE NONL_DEALLOC_SPHPRO


!****************** SUBROUTINE SPHER *********************************
!
!>  Calculates the sperical harmonics multiplied
!>  by the pseudopotential QPROJ on the compressed reciprocal
!>  lattice grid
!>
!>  The full non local pseudopotential is given by
!> ~~~
!>  <G|V(K)|GP> = SUM(L,LP,site,R)
!>                     L                                    LP
!>    QPROJ(G,L site) i * DIJ(L,LP site) QPROJ(GP,LP site) i *
!>                        EXP(-i(G-GP) R)
!> ~~~
!>  i^L is stored in NONL_S\%CQFAK (nonl_struct::CQFAK)
!>
!> @param IZERO = -1   set NONL_S to the negative projectors\n
!>              =  1   set NONL_S to the positive projectors\n
!>              =  0   add to NONL_S the positive projectors
!
!*********************************************************************

  SUBROUTINE SPHER(GRID, NONL_S, P, WDES, LATT_CUR, IZERO, NK_SELECT, DK)

    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant
    USE asa
    USE string, ONLY: str
    USE tutor, ONLY: vtutor

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    INTEGER  IZERO
    INTEGER, OPTIONAL :: NK_SELECT
    REAL(q), OPTIONAL :: DK(3)

    ! work arrays
    REAL(q),ALLOCATABLE :: GLEN(:),XS(:),YS(:),ZS(:),VPS(:),FAKTX(:),YLM(:,:),GLENP(:)

    PROFILING_START('spher')

__GPU CALL OFFLOAD_ON_PUSH(.FALSE.)

    LYDIM=MAXL(NONL_S%NTYP,P)
    LMYDIM=(LYDIM+1)**2          ! number of lm pairs

    NA=WDES%NGDIM

    ALLOCATE(GLEN(NA),XS(NA),YS(NA),ZS(NA),VPS(NA),FAKTX(NA),YLM(NA,LMYDIM))

    IF (PRESENT(DK)) ALLOCATE(GLENP(NA))

#ifdef gammareal
    IF (.NOT.WDES%LGAMMA) THEN
       CALL vtutor%bug("internal error in SPHER: WDES%LGAMMA is incorrect " // str(WDES%LGAMMA), __FILE__, __LINE__)
    ENDIF
#else
    IF (WDES%LGAMMA) THEN
       CALL vtutor%bug("internal error in SPHER: WDES%LGAMMA is incorrect " // str(WDES%LGAMMA), __FILE__, __LINE__)
    ENDIF
#endif

!-----------------------------------------------------------------------
! spin spiral propagation vector in cartesian coordinates
! is simply zero when LSPIRAL=.FALSE.
!-----------------------------------------------------------------------
    QX= (WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3))/2
    QY= (WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3))/2
    QZ= (WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3))/2

!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
    IF (WDES%LSPIRAL) THEN
       NSPINORS=2
    ELSE
       NSPINORS=1
    ENDIF
    spinor: DO ISPINOR=1,NSPINORS
!=======================================================================
! main loop over all special points
!=======================================================================
       kpoint: DO NK=1,WDES%NKPTS
          IF (PRESENT(NK_SELECT)) THEN
             IF ( NK_SELECT/=NK) CYCLE
          ENDIF
!         unfortunately this presently breaks the sc-GW code
!          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!=======================================================================
! now calculate the necessary tables:
! containing the length, phasefactor exp(i m phi) and
! sin(theta),cos(theta)
!=======================================================================
          DO IND=1,WDES%NGVECTOR(NK)

             N1=MOD(WDES%IGX(IND,NK)+GRID%NGX,GRID%NGX)+1
             N2=MOD(WDES%IGY(IND,NK)+GRID%NGY,GRID%NGY)+1
             N3=MOD(WDES%IGZ(IND,NK)+GRID%NGZ,GRID%NGZ)+1

             G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
             G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
             G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))

             IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
                G1= G1 +NONL_S%VKPT_SHIFT(1,NK)
                G2= G2 +NONL_S%VKPT_SHIFT(2,NK)
                G3= G3 +NONL_S%VKPT_SHIFT(3,NK)
             ENDIF

             FACTM=1._q
             IF (WDES%LGAMMA .AND. (N1/=1 .OR. N2/=1 .OR. N3/=1)) FACTM=SQRT(2._q)

             GX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)-QX) *TPI
             GY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)-QY) *TPI
             GZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)-QZ) *TPI


             GLEN(IND)=MAX(SQRT(GX*GX+GY*GY+GZ*GZ),1E-10_q)
             FAKTX(IND)=FACTM
             XS(IND)  =GX/GLEN(IND)
             YS(IND)  =GY/GLEN(IND)
             ZS(IND)  =GZ/GLEN(IND)

             IF (PRESENT(DK)) THEN
                G1P= G1 +DK(1)
                G2P= G2 +DK(2)
                G3P= G3 +DK(3)
                GXP= (G1P*LATT_CUR%B(1,1)+G2P*LATT_CUR%B(1,2)+G3P*LATT_CUR%B(1,3)-QX) *TPI
                GYP= (G1P*LATT_CUR%B(2,1)+G2P*LATT_CUR%B(2,2)+G3P*LATT_CUR%B(2,3)-QY) *TPI
                GZP= (G1P*LATT_CUR%B(3,1)+G2P*LATT_CUR%B(3,2)+G3P*LATT_CUR%B(3,3)-QZ) *TPI
                GLENP(IND)=MAX(SQRT(GXP*GXP+GYP*GYP+GZP*GZP),1E-10_q)
             ENDIF

          ENDDO

          INDMAX=IND-1
!=======================================================================
! now calculate the tables containing the spherical harmonics
! multiply by the pseudopotential and 1/(OMEGA)^(1/2)
!=======================================================================
          CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)

          typ: DO NT=1,NONL_S%NTYP

             LMIND=1
             l_loop: DO L=1,P(NT)%LMAX
!=======================================================================
! first interpolate the non-local pseudopotentials
! and multiply by 1/(OMEGA)^(1/2)
!=======================================================================
                FAKT= 1/SQRT(LATT_CUR%OMEGA)
                ARGSC=NPSNL/P(NT)%PSMAXN
!$OMP SIMD PRIVATE(ARG,NADDR,REM,V1,V2,V3,V4,T0,T1,T2,T3)
                DO IND=1,INDMAX
                   ARG=(GLEN(IND)*ARGSC)+1
                   NADDR=INT(ARG)

                   VPS(IND)=0._q

                   IF (ASSOCIATED(P(NT)%PSPNL_SPLINE)) THEN
                      IF (NADDR<NPSNL) THEN
                         NADDR  =MIN(INT(GLEN(IND)*ARGSC)+1,NPSNL-1)
                         REM=GLEN(IND)-P(NT)%PSPNL_SPLINE(NADDR,1,L)
                         VPS(IND)=(P(NT)%PSPNL_SPLINE(NADDR,2,L)+REM*(P(NT)%PSPNL_SPLINE(NADDR,3,L)+ &
                           &         REM*(P(NT)%PSPNL_SPLINE(NADDR,4,L)+REM*P(NT)%PSPNL_SPLINE(NADDR,5,L))))*FAKT*FAKTX(IND)
                      ENDIF
                   ELSE IF (NADDR<NPSNL-2) THEN
                      REM=MOD(ARG,1.0_q)
                      V1=P(NT)%PSPNL(NADDR-1,L)
                      V2=P(NT)%PSPNL(NADDR,L  )
                      V3=P(NT)%PSPNL(NADDR+1,L)
                      V4=P(NT)%PSPNL(NADDR+2,L)
                      T0=V2
                      T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
                      T2=(V1+V3-(2*V2))/2._q
                      T3=(V4-V1+(3*(V2-V3)))/6._q
                      VPS(IND)=(T0+REM*(T1+REM*(T2+REM*T3)))*FAKT*FAKTX(IND)
                   ENDIF

                   IF (VPS(IND)/=0._q.AND.PRESENT(DK)) THEN
                      ARG=(GLENP(IND)*ARGSC)+1
                      NADDR=INT(ARG)
                      IF (NADDR>NPSNL-1.OR.(.NOT.ASSOCIATED(P(NT)%PSPNL_SPLINE).AND.NADDR>NPSNL-3)) VPS(IND)=0._q
                   ENDIF
                ENDDO
!=======================================================================
! initialize to 0
!=======================================================================
                LL=P(NT)%LPS(L)
                MMAX=2*LL

                IF (IZERO==1 .OR. IZERO==-1) THEN
                   DO LM=0,MMAX
!DIR$ IVDEP
!OCL NOVREC
                      DO IND=1,INDMAX
                         NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)=0
                      ENDDO
                   ENDDO
                ENDIF
                IF (IZERO==-1) THEN
                   DO IND=1,INDMAX
                      VPS(IND)=-VPS(IND)
                   ENDDO
                ENDIF

!=======================================================================
! now multiply with the spherical harmonics
! and set "phase factor" CSET
!=======================================================================
                IF (LL==0) THEN
                   CSET=1.0_q
                ELSE IF (LL==1) THEN
                   CSET=(0.0_q,1.0_q)
                ELSE IF (LL==2) THEN
                   CSET=-1.0_q
                ELSE IF (LL==3) THEN
                   CSET=(0.0_q,-1.0_q)
                ENDIF

                DO LM=0,MMAX
                   NONL_S%CQFAK(LMIND+LM,NT)=CSET
                ENDDO
                LMBASE=LL**2+1

                DO LM=0,MMAX
                   DO IND=1,INDMAX
                      NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)= NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)+ &
                           VPS(IND)*YLM(IND,LM+LMBASE)
                   ENDDO
                ENDDO

                LMIND=LMIND+MMAX+1

             ENDDO l_loop
             IF (LMIND-1/=P(NT)%LMMAX) THEN
                CALL vtutor%bug("internal ERROR: SPHER: LMMAX is wrong " // str(P(NT)%LMMAX) // " " // &
                   str(LMIND-1), __FILE__, __LINE__)
             ENDIF
          ENDDO typ

!=======================================================================
! and of loop over special-points
!=======================================================================
       ENDDO kpoint

       ! conjugate phase alteration for spin down: -q/2 -> q/2
       QX=-QX
       QY=-QY
       QZ=-QZ
    ENDDO spinor

    DEALLOCATE(GLEN,XS,YS,ZS,VPS,FAKTX,YLM)

    IF (PRESENT(DK)) DEALLOCATE(GLENP)

__GPU CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(NONL_S%QPROJ,NONL_S%CQFAK) IF_PRESENT __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(NONL_S%QPROJ) IF(IS_ON_DEVICE(NONL_S%QPROJ).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET UPDATE TO(NONL_S%CQFAK) IF(IS_ON_DEVICE(NONL_S%CQFAK).AND.OFFLOAD_ON)

    PROFILING_STOP('spher')

    RETURN
  END SUBROUTINE SPHER


!****************** SUBROUTINE SPHER_DER *****************************
!
!> Calculates the derivative of the non local projector
!> with respect to k for the cartesian
!> component IDIR times -i = -i nabla_k <p_k|
!>
!> This is equivalent to the projector times r_i  (i=IDIR)
!> ~~~
!> <p r_i| G+k> =             \int p(r) r_i Y_lm(r) e^ir(G+k) dr
!>              = 1/i  d/dk_i \int p(r)     Y_lm(r) e^ir(G+k) dr
!>              = -i d/dk_i <p | G+k>
!> ~~~
!> For simplicity finite differences are used
!
!*********************************************************************

  SUBROUTINE SPHER_DER(GRID,NONL_S,P,WDES,LATT_CUR,  IDIR)
    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant
    USE asa

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    INTEGER IDIR
    ! local
    REAL(q),POINTER  :: VKPT_ORIG(:,:)
    REAL(q) :: V(3)
    ! displacement for finite differences
    REAL(q), PARAMETER :: DIS=fd_displacement
    INTEGER IDIS,NK

    PROFILING_START('spher_der')

    ! reset all entries
    NONL_S%QPROJ=0

    ! first generate a pointer to the original k-point array
    VKPT_ORIG=>WDES%VKPT
    NULLIFY(WDES%VKPT)
    ALLOCATE(WDES%VKPT(3,WDES%NKPTS))

    ! k-displacement
    V=0
    V(IDIR)=DIS/TPI
    CALL KARDIR(1,V(1),LATT_CUR%A)

    ! negative displacement
    DO NK=1,WDES%NKPTS
       WDES%VKPT(:,NK)=VKPT_ORIG(:,NK)-V
    ENDDO
!   CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,-1)
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,-1,DK= 2*V)

    ! positive displacement
    DO NK=1,WDES%NKPTS
       WDES%VKPT(:,NK)=VKPT_ORIG(:,NK)+V
    ENDDO
!   CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR, 0)
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR, 0,DK=-2*V)

    NONL_S%QPROJ=NONL_S%QPROJ/DIS/2

    NONL_S%CQFAK=NONL_S%CQFAK*(0.0_q,-1.0_q)


    DEALLOCATE(WDES%VKPT)
    WDES%VKPT=>VKPT_ORIG

    PROFILING_STOP('spher_der')

  END SUBROUTINE SPHER_DER


  SUBROUTINE NONL_ALLOC_DER(NONL_S, NONL_D)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S, NONL_D
    ! local var
    INTEGER NTYPD, NRPLWV, LMDIM, NKPTS

    NRPLWV=SIZE(NONL_S%QPROJ,1)
    LMDIM =SIZE(NONL_S%QPROJ,2)
    NTYPD =SIZE(NONL_S%QPROJ,3)
    NKPTS =SIZE(NONL_S%QPROJ,4)

    NONL_D=NONL_S

    NULLIFY(NONL_D%QPROJ,NONL_D%CQFAK)
    IF (NONL_S%LSPIRAL) THEN
       ALLOCATE(NONL_D%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
            NONL_D%CQFAK(LMDIM,NTYPD))
    ELSE
       ALLOCATE(NONL_D%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
            NONL_D%CQFAK(LMDIM,NTYPD))
    ENDIF
DOESI NONL_D%QPROJ = 0.0_q; NONL_D%CQFAK = (0.0_q, 0.0_q) ! Breaks ???, operations on entire uninitialized arrays in SPHER_DER.

    RETURN
  END SUBROUTINE NONL_ALLOC_DER

  SUBROUTINE NONL_DEALLOC_DER(NONL_D)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_D
    DEALLOCATE(NONL_D%QPROJ,NONL_D%CQFAK)

    RETURN
  END SUBROUTINE NONL_DEALLOC_DER


!************************ SUBROUTINE PHASE *****************************
!
!> Calculates the phasefactor CREXP (<TT>exp(ig.r)</TT>) for one k-point,
!> on the compressed grid of lattice vectors
!>
!> @note CREXP is only calculated if NK changes, if the ionc positions
!> change PHASE must be called with NK=0 to force the routine to recalculate
!> the phase-factor in the next call
!
!***********************************************************************

  SUBROUTINE PHASE(WDES,NONL_S,NK)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER PHASE
#endif
    USE constant

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavedes)     WDES
    TYPE (nonl_struct) NONL_S

    PROFILING_START('phase')

!=======================================================================
! check if special point changed
!=======================================================================
    IF (NK==0 .OR.NK==NONL_S%NK) THEN
       NONL_S%NK=NK
GPACC !$ACC UPDATE DEVICE(NONL_S%NK) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(NONL_S%NK) IF(OFFLOAD_ON)
       PROFILING_STOP('phase')
       RETURN
    ENDIF
    NONL_S%NK=NK
GPACC !$ACC UPDATE DEVICE(NONL_S%NK) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(NONL_S%NK) IF(OFFLOAD_ON)

    ! number of G-vectors
    NPL= WDES%NGVECTOR(NK)

    ! set the phase-factor
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(NONL_S,WDES) &
GPACC !$ACC PRIVATE(GXDX,GYDY,GZDZ,CGDR) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP PRIVATE(GXDX,GYDY,GZDZ,CGDR) &
GPOMP !$OMP MAP(PRESENT,TO:NONL_S,WDES) IF(OFFLOAD_ON)
    ion:  DO NI=1,NONL_S%NIONS
__CPU  GXDX=NONL_S%POSION(1,NI)
__CPU  GYDY=NONL_S%POSION(2,NI)
__CPU  GZDZ=NONL_S%POSION(3,NI)
__CPU !$OMP SIMD PRIVATE(CGDR)
       DO M=1,NPL
__GPU     GXDX=NONL_S%POSION(1,NI)
__GPU     GYDY=NONL_S%POSION(2,NI)
__GPU     GZDZ=NONL_S%POSION(3,NI)
          CGDR=CITPI*(WDES%IGX(M,NK)*GXDX+WDES%IGY(M,NK)*GYDY+WDES%IGZ(M,NK)*GZDZ)
          NONL_S%CREXP(M,NI)=EXP(CGDR)
       ENDDO
    ENDDO ion

GPACC !$ACC UPDATE SELF(NONL_S%CREXP) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE FROM(NONL_S%CREXP) IF(OFFLOAD_ON)

    PROFILING_STOP('phase')

    RETURN
  END SUBROUTINE PHASE


!************************ SUBROUTINE VNLACC ****************************
!
!> Calculates the non local contribution of the Hamiltonian,
!> using the reciprocal space projection scheme
!>
!> The projection of the wavefunction onto the projection operators
!> must be supplied in W1\%CPROJ
!
!***********************************************************************

  SUBROUTINE VNLACC(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    REAL(q) :: EVALUE
    OVERLAP :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
  ! local
    GDEF :: CRESUL(W1%WDES1%NPRO)

GPACC !$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CRESUL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CACC) __IF_ASYNC__
#ifdef IFXWAR3
GPOMP !$OMP TARGET IF(OFFLOAD_ON)
#else
GPOMP !$OMP TARGET MAP(PRESENT,TO:CACC) IF(OFFLOAD_ON)
#endif
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CACC=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    CALL OVERL1(W1%WDES1, SIZE(CDIJ,1),CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1), CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))

GPACC !$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CRESUL) IF(OFFLOAD_ON)
    RETURN
  END SUBROUTINE VNLACC


  SUBROUTINE VNLACC_ADD(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    REAL(q) :: EVALUE
    OVERLAP :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
    ! work array
    GDEF :: CRESUL(W1%WDES1%NPRO)

GPACC !$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CRESUL) IF(OFFLOAD_ON)
    CALL OVERL1(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1), CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
GPACC !$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CRESUL) IF(OFFLOAD_ON)
    RETURN
  END SUBROUTINE VNLACC_ADD


  SUBROUTINE VNLACC_C(NONL_S, W1, CDIJ, CQIJ, ISP, CEVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER :: LMDIM
    COMPLEX(q) :: CEVALUE
    OVERLAP :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
    ! work array
    GDEF :: CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
!$ACC KERNELS PRESENT(CACC) __IF_ASYNC__
    CACC=0
!$ACC END KERNELS
    CALL OVERL1_C(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), CEVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
!$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__

    RETURN
  END SUBROUTINE VNLACC_C


  SUBROUTINE VNLACC_ADD_C(NONL_S,W1, CDIJ, CQIJ, ISP, CEVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER :: LMDIM
    COMPLEX(q) :: CEVALUE
    OVERLAP :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
    ! work array
    GDEF :: CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
    CALL OVERL1_C(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), CEVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))
!$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__

    RETURN
  END SUBROUTINE VNLACC_ADD_C


! here CDIJ and QCIJ are always complex
  SUBROUTINE VNLACC_ADD_CCDIJ(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CACC)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    REAL(q) :: EVALUE
    COMPLEX(q) :: CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CACC(:)
    ! work array
    GDEF :: CRESUL(W1%WDES1%NPRO)

    CALL OVERL1_CCDIJ(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0(NONL_S,W1%WDES1,CRESUL(1),CACC(1))

    RETURN
  END SUBROUTINE VNLACC_ADD_CCDIJ


!************************ SUBROUTINE PROJ ******************************
!
!> Calculates the projection of all bands of one specific
!> k-point onto the reciprocal space projection operators
!
!***********************************************************************

  SUBROUTINE PROJ(NONL_S,WDES,W,NK)
    USE mopenmp_struct_def
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    INTEGER NK
  ! local
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    INTEGER ISP, N

    PROFILING_START('proj')

GPACC !$ACC ENTER DATA CREATE(WDES1,W1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1,W1) IF(OFFLOAD_ON)
    CALL SETWDES(WDES,WDES1,NK)

    DO ISP=1,WDES%ISPIN
       DO N=1,WDES%NBANDS
          CALL SETWAV(W,W1,WDES1,N,ISP)
          CALL PROJ1(NONL_S,WDES1,W1)
       ENDDO
    ENDDO

#ifdef OMP_OFFLOAD
    IF (ASSOCIATED(W1%CPROJ)) THEN
!$OMP TARGET EXIT DATA MAP(RELEASE:W1, W1%CPROJ) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE:W1%CPROJ) IF(OFFLOAD_ON)
    ENDIF
    IF (ASSOCIATED(W1%CPTWFP)) THEN
!$OMP TARGET EXIT DATA MAP(RELEASE:W1, W1%CPTWFP) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE:W1%CPTWFP) IF(OFFLOAD_ON)
    ENDIF
    IF (ASSOCIATED(W1%WDES1)) THEN
!$OMP TARGET EXIT DATA MAP(RELEASE:W1, W1%WDES1) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE:W1%WDES1) IF(OFFLOAD_ON)
    ENDIF
!$OMP TARGET EXIT DATA MAP(RELEASE:W1) IF(OFFLOAD_ON)
#endif

GPACC !$ACC EXIT DATA DELETE(W1%CPTWFP,W1%CPROJ,W1%WDES1) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(W1) __IF_ASYNC__
__GPU CALL DELETE_TYPED_VAR(WDES1)

    PROFILING_STOP('proj')

  END SUBROUTINE PROJ


!************************ SUBROUTINE PROJ1 *****************************
!
!> Calculates the scalar product of one wavefunction with
!> all projector functions in reciprocal space
!> (thesis gK Eq. 10.34)
!>
!> @details @ref openmp :
!> under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is contracted to a single loop over \"all ions\", and the
!> latter is distributed over all available threads.
!
!***********************************************************************

  SUBROUTINE PROJ1(NONL_S,WDES1,W1)
    IMPLICIT NONE

    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes1)    :: WDES1
    TYPE (wavefun1), TARGET :: W1

#ifdef OFFLOADING
    IF (OFFLOAD_ON) THEN
       CALL PROJ1_ACC(NONL_S,WDES1,W1)
       RETURN
    ENDIF
#endif
    CALL PROJ1_(NONL_S,WDES1,W1)

    RETURN
  END SUBROUTINE PROJ1

  SUBROUTINE PROJ1_(NONL_S,WDES1,W1)
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1), TARGET :: W1

  ! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, LMMAXC, NI, K, LM
    REAL(q) :: SUMR, SUMI
    REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)

    PROFILING_START('proj1')

    IF (WDES1%NK /= NONL_S%NK) THEN
       CALL vtutor%bug("internal error in PROJ1: PHASE not properly set up " // str(WDES1%NK) // " " &
          // str(NONL_S%NK), __FILE__, __LINE__)
    ENDIF

    NPL=WDES1%NGVECTOR

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
!=======================================================================
! performe loop over ions
!=======================================================================
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NI,NT,LMMAXC,LMBASE,WORK,TMP,LM,SUMR,SUMI)
       ion: DO NI=1,NONL_S%NIONS
          NT=NONL_S%ITYP(NI)
          LMMAXC=NONL_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion
          LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
!=======================================================================
! multiply with phasefactor and divide into real and imaginary part
!=======================================================================
          CALL CREXP_MUL_WAVE( NPL, NONL_S%CREXP(1,NI), W1%CPTWFP(NPL*ISPINOR+1), WORK(1))
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#if defined(do_loops) || defined(pro_loop)
          l_loop: DO LM=1,LMMAXC
             SUMR=0
             SUMI=0
!DIR$ IVDEP
!OCL NOVREC
             DO K=1,NPL
                SUMR = SUMR + WORK(K)     * NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
                SUMI = SUMI + WORK(K+NPL) * NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
             ENDDO
             CPROJ(LM+LMBASE)=GREAL(CMPLX(SUMR, SUMI, KIND=q) *NONL_S%CQFAK(LM,NT))
          ENDDO l_loop
#else
          PROFILING_START('proj1|dgemv')
          CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, WORK(1) , 1 , 0._q ,  TMP(1,1), 1)
          CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, WORK(1+NPL) , 1 , 0._q ,  TMP(1,2), 1)
          PROFILING_STOP('proj1|dgemv')

          l_loop: DO LM=1,LMMAXC
             SUMR=TMP(LM,1)
             SUMI=TMP(LM,2)
             CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *NONL_S%CQFAK(LM,NT))
          ENDDO l_loop
#endif
       ENDDO ion
!$OMP END PARALLEL DO
    ENDDO spinor

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

    PROFILING_STOP('proj1')

    RETURN
  END SUBROUTINE PROJ1_


!************************ SUBROUTINE PROJ1_ACC *************************
!
!***********************************************************************
#ifdef OFFLOADING
  SUBROUTINE PROJ1_ACC(NONL_S,WDES1,W1)
!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ
    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1), TARGET :: W1

  ! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, NI, K, LM, IBLOCK
    COMPLEX(q) :: CTMP
    GDEF :: CPROJ(WDES1%NPRO_TOT)

    INTEGER, PARAMETER :: BLOCKSIZE=512
    COMPLEX(q) :: WORK(BLOCKSIZE)

    PROFILING_START('proj1_acc')

    IF (WDES1%NK /= NONL_S%NK) THEN
       CALL vtutor%bug("internal error in PROJ1_ACC: PHASE not properly set up " // str(WDES1%NK) // " " &
          // str(NONL_S%NK), __FILE__, __LINE__)
    ENDIF

    NPL=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CPROJ) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CPROJ=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
!=======================================================================
! performe loop over ions
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG __IF_ASYNC__ &
GPACC !$ACC PRESENT(NONL_S,W1,CPROJ) PRIVATE(NT,LMBASE,WORK)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
GPOMP !$OMP PRIVATE(NT,LMBASE,WORK) MAP(PRESENT,TO:NONL_S,W1,CPROJ) IF(OFFLOAD_ON)

       block: DO IBLOCK=0,NPL-1,BLOCKSIZE
       ion: DO NI=1,NONL_S%NIONS
          NT=NONL_S%ITYP(NI)
          IF (NONL_S%LMMAX(NT)==0) CYCLE ion
          LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
!=======================================================================
! multiply with phasefactor and divide into real and imaginary part
!=======================================================================
#ifndef CRAYHIP
GPACC !$ACC CACHE(WORK(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
          DO K=IBLOCK+1,MIN(IBLOCK+BLOCKSIZE,NPL)
             WORK(K-IBLOCK) = W1%CPTWFP(K+NPL*ISPINOR) * NONL_S%CREXP(K,NI)
          ENDDO
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
GPACC !$ACC LOOP WORKER PRIVATE(CTMP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP) IF(OFFLOAD_ON)
          l_loop: DO LM=1,NONL_S%LMMAX(NT)
             CTMP=0
GPACC !$ACC LOOP VECTOR REDUCTION(+:CTMP)
!DIR$ IVDEP
!OCL NOVREC
             DO K=IBLOCK+1,MIN(IBLOCK+BLOCKSIZE,NPL)
                CTMP = CTMP + WORK(K-IBLOCK) * NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
             ENDDO
#ifdef gammareal
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
             CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+GREAL(CTMP *NONL_S%CQFAK(LM,NT))
#else // gammareal
             CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CPROJ(LM+LMBASE),CTMP *NONL_S%CQFAK(LM,NT))
#endif // gammareal
          ENDDO l_loop
       ENDDO ion
       ENDDO block
    ENDDO spinor

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

GPACC !$ACC EXIT DATA DELETE(CPROJ,WORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CPROJ,WORK) IF(OFFLOAD_ON)

    PROFILING_STOP('proj1_acc')

    RETURN
  END SUBROUTINE PROJ1_ACC
#endif

!************************ SUBROUTINE PROJMU ****************************
!
!***********************************************************************

  SUBROUTINE PROJMU(NONL_S,WDES1,W1)
    IMPLICIT NONE
    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes1)    :: WDES1
    TYPE (wavefun1)    :: W1(:)

  ! local
    INTEGER :: NSIM,NP

    PROFILING_START('projmu')

    NSIM=SIZE(W1)

#ifdef OFFLOADING
    IF (OFFLOAD_ON) THEN
       CALL PROJMU_ACC(NONL_S,WDES1,W1,NSIM,W1%LDO)
       PROFILING_STOP('projmu')
    ENDIF
#endif
    DO NP=1,NSIM
       IF (W1(NP)%LDO) CALL PROJ1(NONL_S,WDES1,W1(NP))
    ENDDO

    PROFILING_STOP('projmu')

    RETURN
  END SUBROUTINE PROJMU


!************************ SUBROUTINE PROJMU_ACC ************************
!
!***********************************************************************
#ifdef OFFLOADING
  SUBROUTINE PROJMU_ACC(NONL_S,WDES1,W1,NSIM,LDO)
!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ
    IMPLICIT NONE

    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes1)    :: WDES1
    TYPE (wavefun1)    :: W1(NSIM)
    INTEGER            :: NSIM
    LOGICAL            :: LDO(NSIM)

  ! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, NI, K, LM, IBLOCK, NP
    COMPLEX(q) :: CTMP
    GDEF :: CPROJ(WDES1%NPRO_TOT,NSIM)

    INTEGER, PARAMETER :: BLOCKSIZE=512
    COMPLEX(q) :: WORK(BLOCKSIZE)

    PROFILING_START('projmu_acc')

    IF (WDES1%NK /= NONL_S%NK) &
       CALL vtutor%bug("PROJMU_ACC: PHASE not properly set up " // str(WDES1%NK) // " " // &
          str(NONL_S%NK), __FILE__, __LINE__)

    NPL=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CPROJ) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    CPROJ=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC ENTER DATA COPYIN(LDO) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:LDO) IF(OFFLOAD_ON)

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
!=======================================================================
! performe loop over ions
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG __IF_ASYNC__ &
GPACC !$ACC PRESENT(NONL_S,W1,CPROJ,LDO) PRIVATE(NT,LMBASE,WORK)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) THREAD_LIMIT(64) IF(OFFLOAD_ON) &
#ifdef IFXWAR3
GPOMP !$OMP PRIVATE(NT,LMBASE,WORK,CTMP) MAP(PRESENT,TO:NONL_S,CPROJ,LDO)
#else
GPOMP !$OMP PRIVATE(NT,LMBASE,WORK,CTMP) MAP(PRESENT,TO:NONL_S,W1,CPROJ,LDO)
#endif

       block: DO IBLOCK=0,NPL-1,BLOCKSIZE
       ion: DO NI=1,NONL_S%NIONS
          sim: DO NP=1,NSIM

             NT=NONL_S%ITYP(NI)

             IF (NONL_S%LMMAX(NT)==0 .OR. .NOT.LDO(NP)) CYCLE

             LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
!=======================================================================
! multiply with phasefactor and divide into real and imaginary part
!=======================================================================
#ifndef CRAYHIP
GPACC !$ACC CACHE(WORK(1:BLOCKSIZE))
#endif
GPACC !$ACC LOOP WORKER VECTOR
GPOMP !$OMP PARALLEL DO IF(OFFLOAD_ON)
             DO K=IBLOCK+1,MIN(IBLOCK+BLOCKSIZE,NPL)
                WORK(K-IBLOCK) = W1(NP)%CPTWFP(K+NPL*ISPINOR) * NONL_S%CREXP(K,NI)
             ENDDO
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
GPACC !$ACC LOOP WORKER PRIVATE(CTMP)
GPOMP !$OMP PARALLEL DO PRIVATE(CTMP) IF(OFFLOAD_ON)
             l_loop: DO LM=1,NONL_S%LMMAX(NT)
                CTMP=0
GPACC !$ACC LOOP VECTOR REDUCTION(+:CTMP)
!DIR$ IVDEP
!OCL NOVREC
                DO K=IBLOCK+1,MIN(IBLOCK+BLOCKSIZE,NPL)
                   CTMP = CTMP + WORK(K-IBLOCK) * NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
                ENDDO
#ifdef gammareal
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+GREAL(CTMP*NONL_S%CQFAK(LM,NT))
#else // gammareal
                CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CPROJ(LM+LMBASE,NP),CTMP*NONL_S%CQFAK(LM,NT))
#endif // gammareal
             ENDDO l_loop

          ENDDO sim
       ENDDO ion
       ENDDO block

    ENDDO spinor

    ! distribute the projected wavefunctions to nodes
    DO NP=1,NSIM
       IF (LDO(NP)) CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
    ENDDO

GPACC !$ACC EXIT DATA DELETE(LDO,CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:LDO,CPROJ) IF(OFFLOAD_ON)

    PROFILING_STOP('projmu_acc')

    RETURN
  END SUBROUTINE PROJMU_ACC
#endif

!************************ SUBROUTINE PROJXYZ ***************************
!
!> Calculates the first order change of the wave function character
!> upon moving the ions for one selected k-point and spin component
!>
!> The results are stored in CPROJXYZ
!>
!> A factor 2 stemming from variations of the bra or kat is included
!>
!> @param LADDK (optional) LADDK=.TRUE.
!> forces the treatment usually done only for <TT>#define shift_der_k</TT>.
!> This allows for correct results for the Born effective charges
!> even if <TT>shift_der_k</TT> is not defined
!>
!> @details @ref openmp :
!> loops over the basis vectors of the orbitals are distributed
!> over all available threads.
!
!***********************************************************************

  SUBROUTINE PROJXYZ(NONL_S, WDES, W, LATT_CUR, ISP, NK, CPROJXYZ, LADDK)
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes)     :: WDES
    TYPE (wavespin)    :: W
    TYPE (latt)        :: LATT_CUR

    GDEF    :: CPROJXYZ(WDES%NPROD, WDES%NBANDS, 3)
    INTEGER :: ISP, NK

    LOGICAL, OPTIONAL :: LADDK

    ! local variable
    TYPE (wavedes1) :: WDES1  ! descriptor for one k-point

    COMPLEX(q) :: CMUL
    INTEGER    :: NPL,N,LMBASE,LMMAXC,ISPINOR,ISPIRAL,K,KK,NI,NT,LM

    GDEF :: CX_,CY_,CZ_,CVAL
#ifdef OFFLOADING
    GDEF :: CX(WDES%NPRO_TOT,WDES%NBANDS),CY(WDES%NPRO_TOT,WDES%NBANDS),CZ(WDES%NPRO_TOT,WDES%NBANDS)
#else
    GDEF :: CX(WDES%NPRO_TOT),CY(WDES%NPRO_TOT),CZ(WDES%NPRO_TOT)
#endif

    LOGICAL :: LVKPT

    PROFILING_START('projxyz')

    NPL=WDES%NGVECTOR(NK)

    LVKPT=.FALSE. ; IF (PRESENT(LADDK)) LVKPT=LADDK

GPACC !$ACC ENTER DATA CREATE(CX,CY,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CX,CY,CZ) IF(OFFLOAD_ON)

GPACC !$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(WDES,WDES1,NK)
    CALL PHASE(WDES,NONL_S,NK)

#ifdef OFFLOADING
GPACC !$ACC KERNELS PRESENT(CX,CY,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CX,CY,CZ) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
       CX=0; CY=0; CZ=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#endif

__CPU !!!$OMP PARALLEL DO PRIVATE(N,CX,CY,CZ,LMBASE,ISPIRAL,ISPINOR,NI,NT,LMMAXC,LM,CMUL,CX_,CY_,CZ_,K,KK,CVAL) &
__CPU !!!$OMP DEFAULT(SHARED) SCHEDULE(static)
GPACC !$ACC PARALLEL LOOP GANG COLLAPSE(3) PRESENT(NONL_S,W,WDES,LATT_CUR,CX,CY,CZ) &
GPACC !$ACC PRIVATE(ISPIRAL,NT,LMMAXC,LMBASE) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:NONL_S,W,WDES,LATT_CUR,CX,CY,CZ) &
GPOMP !$OMP PRIVATE(ISPIRAL,NT,LMMAXC,LMBASE) IF(OFFLOAD_ON)
    band: DO N=1,WDES%NBANDS

__CPU  CX=0; CY=0; CZ=0

       spinor: DO ISPINOR=0,WDES%NRSPINORS-1
          ion: DO NI=1,NONL_S%NIONS

             ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1

             NT=NONL_S%ITYP(NI)
             LMMAXC=NONL_S%LMMAX(NT)
             IF (LMMAXC==0) CYCLE ion
             LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)

GPACC !$ACC LOOP PRIVATE(CMUL,CX_,CY_,CZ_)
GPOMP !$OMP PARALLEL DO PRIVATE(CMUL,CX_,CY_,CZ_) IF(OFFLOAD_ON)
             l_loop: DO LM=1,LMMAXC
                CMUL=NONL_S%CQFAK(LM,NT)*CITPI
                CX_=0; CY_=0; CZ_=0

                IF (LVKPT) THEN
GPACC !$ACC LOOP VECTOR PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
__CPU !$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(__simd_static) &
__CPU !$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CVAL=2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*W%CPTWFP(KK,N,NK,ISP)*NONL_S%CREXP(K,NI)*CMUL
                      CX_=CX_+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                      CY_=CY_+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                      CZ_=CZ_+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
                   ENDDO
__CPU !$OMP END PARALLEL DO SIMD
                ELSE
GPACC !$ACC LOOP VECTOR PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
__CPU !$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(__simd_static) &
__CPU !$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CVAL=2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*W%CPTWFP(KK,N,NK,ISP)*NONL_S%CREXP(K,NI)*CMUL
#ifdef shift_der_k
                      CX_=CX_+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                      CY_=CY_+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                      CZ_=CZ_+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
#else
                      CX_=CX_+(WDES%IGX(K,NK))*CVAL
                      CY_=CY_+(WDES%IGY(K,NK))*CVAL
                      CZ_=CZ_+(WDES%IGZ(K,NK))*CVAL
#endif
                   ENDDO
__CPU !$OMP END PARALLEL DO SIMD
                ENDIF

#ifdef OFFLOADING
                CX(LM+LMBASE,N)=CX_*LATT_CUR%B(1,1)+CY_*LATT_CUR%B(1,2)+CZ_*LATT_CUR%B(1,3)
                CY(LM+LMBASE,N)=CX_*LATT_CUR%B(2,1)+CY_*LATT_CUR%B(2,2)+CZ_*LATT_CUR%B(2,3)
                CZ(LM+LMBASE,N)=CX_*LATT_CUR%B(3,1)+CY_*LATT_CUR%B(3,2)+CZ_*LATT_CUR%B(3,3)
#else
                CX(LM+LMBASE)=CX_*LATT_CUR%B(1,1)+CY_*LATT_CUR%B(1,2)+CZ_*LATT_CUR%B(1,3)
                CY(LM+LMBASE)=CX_*LATT_CUR%B(2,1)+CY_*LATT_CUR%B(2,2)+CZ_*LATT_CUR%B(2,3)
                CZ(LM+LMBASE)=CX_*LATT_CUR%B(3,1)+CY_*LATT_CUR%B(3,2)+CZ_*LATT_CUR%B(3,3)
#endif

             ENDDO l_loop
          ENDDO ion
       ENDDO spinor

#ifndef OFFLOADING
       CALL DIS_PROJ(WDES1,CX,CPROJXYZ(1,N,1))
       CALL DIS_PROJ(WDES1,CY,CPROJXYZ(1,N,2))
       CALL DIS_PROJ(WDES1,CZ,CPROJXYZ(1,N,3))
#endif

    ENDDO band
__CPU !!!$OMP END PARALLEL DO

#ifdef OFFLOADING
    DO N=1,WDES%NBANDS
       CALL DIS_PROJ(WDES1,CX(1,N),CPROJXYZ(1,N,1))
       CALL DIS_PROJ(WDES1,CY(1,N),CPROJXYZ(1,N,2))
       CALL DIS_PROJ(WDES1,CZ(1,N),CPROJXYZ(1,N,3))
    ENDDO

GPACC !$ACC EXIT DATA DELETE(CX,CY,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CX,CY,CZ) IF(OFFLOAD_ON)
    CALL DELETE_TYPED_VAR(WDES1)
#endif

    PROFILING_STOP('projxyz')

    RETURN
  END SUBROUTINE PROJXYZ


!************************ SUBROUTINE PROJXYZ_WA ************************
!
!> Same as #PROJXYZ, however, is accepts a single array for the
!> orbitals' PW coefficients instead  of a wavefunction array W
!>
!> This subroutine is used by the module greens_real_space
!>
!> @details @ref openmp :
!> loops over the basis vectors of the orbitals are distributed
!> over all available threads.
!
!***********************************************************************

  SUBROUTINE PROJXYZ_WA(NONL_S, WDES, CW, LATT_CUR, NK, CPROJXYZ, LADDK)
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    COMPLEX(q),CONTIGUOUS :: CW(:,:)
    TYPE (latt)        LATT_CUR
    INTEGER            ISP, NK
    GDEF, CONTIGUOUS :: CPROJXYZ(:,:,:)
    LOGICAL, OPTIONAL :: LADDK
    ! local variable
    TYPE (wavedes1)    WDES1         ! descriptor for one k-point
    GDEF       :: CX(WDES%NPRO_TOT) ,CY(WDES%NPRO_TOT) ,CZ(WDES%NPRO_TOT)
    GDEF       :: CX_, CY_, CZ_
    INTEGER    :: NPRO, NPL, N, LMBASE, LMMAXC, ISPINOR, ISPIRAL, K, KK, NI, NIS, NT, LM
    COMPLEX(q) :: CMUL
    GDEF       :: CVAL

    PROFILING_START('projxyz')

    NPRO=WDES%NPRO_TOT
    NPL= WDES%NGVECTOR(NK)

    CALL SETWDES(WDES,WDES1,NK)
    CALL PHASE(WDES,NONL_S,NK)

!!!$OMP PARALLEL DO PRIVATE(N,CX,CY,CZ,LMBASE,ISPIRAL,ISPINOR,NIS,NT,LMMAXC,NI,LM,CMUL,CX_,CY_,CZ_,K,KK,CVAL) &
!!!$OMP DEFAULT(SHARED) SCHEDULE(static)
    band: DO N=1,WDES%NBANDS

       CX=0
       CY=0
       CZ=0

       LMBASE= 0

       ISPIRAL=1
       spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
          NIS=1
          typ:  DO NT=1,NONL_S%NTYP
             LMMAXC=NONL_S%LMMAX(NT)
             IF (LMMAXC==0) GOTO 100

             ions:   DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
                l_loop: DO LM=1,LMMAXC
                   CMUL=NONL_S%CQFAK(LM,NT)*CITPI
                   CX_=0
                   CY_=0
                   CZ_=0
                   IF (PRESENT (LADDK)) THEN
!$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(__simd_static) &
!$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CVAL =  2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*CW(KK,N)*NONL_S%CREXP(K,NI)*CMUL
                      CX_=CX_+CVAL*(WDES%IGX(K,NK)+WDES%VKPT(1,NK))
                      CY_=CY_+CVAL*(WDES%IGY(K,NK)+WDES%VKPT(2,NK))
                      CZ_=CZ_+CVAL*(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))
                   ENDDO
!$OMP END PARALLEL DO SIMD
                   ELSE
!$OMP PARALLEL DO SIMD DEFAULT(SHARED) SCHEDULE(__simd_static) &
!$OMP PRIVATE(KK,CVAL) REDUCTION(+:CX_,CY_,CZ_)
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      ! factor 2 because either the bra or kat can vary
                      CVAL =  2*NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*CW(KK,N)*NONL_S%CREXP(K,NI)*CMUL
#ifdef shift_der_k
                      CX_=CX_+CVAL*(WDES%IGX(K,NK)+WDES%VKPT(1,NK))
                      CY_=CY_+CVAL*(WDES%IGY(K,NK)+WDES%VKPT(2,NK))
                      CZ_=CZ_+CVAL*(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))
#else
                      CX_=CX_+CVAL*(WDES%IGX(K,NK))
                      CY_=CY_+CVAL*(WDES%IGY(K,NK))
                      CZ_=CZ_+CVAL*(WDES%IGZ(K,NK))
#endif
                   ENDDO
!$OMP END PARALLEL DO SIMD
                   ENDIF

                   CX(LM+LMBASE)=CX_*LATT_CUR%B(1,1)+CY_*LATT_CUR%B(1,2)+CZ_*LATT_CUR%B(1,3)
                   CY(LM+LMBASE)=CX_*LATT_CUR%B(2,1)+CY_*LATT_CUR%B(2,2)+CZ_*LATT_CUR%B(2,3)
                   CZ(LM+LMBASE)=CX_*LATT_CUR%B(3,1)+CY_*LATT_CUR%B(3,2)+CZ_*LATT_CUR%B(3,3)

                ENDDO l_loop
                LMBASE= LMMAXC+LMBASE
             ENDDO ions

100          NIS = NIS+NONL_S%NITYP(NT)
          ENDDO typ

          IF (NONL_S%LSPIRAL) ISPIRAL=2
       ENDDO spinor

       CALL DIS_PROJ(WDES1,CX,CPROJXYZ(1,N,1))
       CALL DIS_PROJ(WDES1,CY,CPROJXYZ(1,N,2))
       CALL DIS_PROJ(WDES1,CZ,CPROJXYZ(1,N,3))

    ENDDO band
!!!$OMP END PARALLEL DO

    PROFILING_STOP('projxyz')

    RETURN
  END SUBROUTINE PROJXYZ_WA


!************************ SUBROUTINE PROJLAT_DER ***********************
!
!> Calculates the first order change of the wave function character
!> upon changing the lattice
!>
!> The results are stored in CPROJXYZ
!
!***********************************************************************

  SUBROUTINE PROJLAT_DER(P, NONL_S, WDES, W, LATT_CUR, ISP, NK, CPROJXYZ)
    USE lattice
    USE constant
    USE pseudo
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc

    IMPLICIT NONE

    TYPE (potcar)      P(:)
    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR
    INTEGER            ISP, NK
    GDEF, TARGET :: CPROJXYZ(WDES%NPROD, WDES%NBANDS, 6)
    ! local variable
    TYPE (wavedes1)    WDES1         ! descriptor for one k-point
    TYPE (wavefun1)    W1
    INTEGER    :: N
    INTEGER :: IDIR, JDIR, IJDIR
    TYPE (latt)        LATT_FIN1,LATT_FIN2
    REAL(q) :: DIS=fd_displacement

    PROFILING_START('projlat_der')

GPACC !$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(WDES,WDES1,NK)
    CALL PHASE(WDES,NONL_S,NK)

DOESI CPROJXYZ = 0

    IJDIR=0
    DO IDIR=1,3
       DO JDIR=1,IDIR
          IJDIR=IJDIR+1
          LATT_FIN1=LATT_CUR
          LATT_FIN1%A(IDIR,:)=LATT_CUR%A(IDIR,:)+DIS*LATT_CUR%A(JDIR,:)

          LATT_FIN2=LATT_CUR
          LATT_FIN2%A(IDIR,:)=LATT_CUR%A(IDIR,:)-DIS*LATT_CUR%A(JDIR,:)

          CALL LATTIC(LATT_FIN1)
          CALL LATTIC(LATT_FIN2)

          CALL SPHER(WDES%GRID,NONL_S,P,WDES,LATT_FIN2, IZERO=-1, NK_SELECT=NK)
          CALL SPHER(WDES%GRID,NONL_S,P,WDES,LATT_FIN1, IZERO=0 , NK_SELECT=NK)

!!GPACC !$OMP PARALLEL DO PRIVATE(W1) NUM_THREADS(omp_nthreads_acc)
          band: DO N=1,WDES%NBANDS

!!             SET_ACC_ASYNC_Q(N)

             CALL SETWAV(W,W1,WDES1,N,ISP); W1%CPROJ => CPROJXYZ(:,N,IJDIR)
#ifdef OFFLOADING
             CALL COPYIN_TYPED_VAR(W1)
#endif
             CALL PROJ1(NONL_S,WDES1,W1) ! calculate W1%CPROJ (linked to CPROJXYZ)

GPACC !$ACC KERNELS PRESENT(CPROJXYZ) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CPROJXYZ) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
             CPROJXYZ(:,N,IJDIR)=CPROJXYZ(:,N,IJDIR)/DIS
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#ifdef OFFLOADING
            CALL DELETE_TYPED_VAR(W1)
#endif
          ENDDO band
!!GPACC !$OMP END PARALLEL DO
       ENDDO
    ENDDO

    CALL SPHER(WDES%GRID, NONL_S, P, WDES, LATT_CUR,  1, NK_SELECT=NK)

#ifdef OFFLOADING
    CALL DELETE_TYPED_VAR(WDES1)
#endif

    PROFILING_STOP('projlat_der')

    RETURN
  END SUBROUTINE PROJLAT_DER


!************************ SUBROUTINE PROJ_DER **************************
!
!> Calculates the first order change of the wave function character
!> upon moving one ion ION in the direction IDIR
!>
!> A term  (i k)  is included as well to "center" the derivate at the
!> ion
!
!***********************************************************************

  SUBROUTINE PROJ_DER(NONL_S, WDES1, W1, LATT_CUR, ION, IDIR)

GPACC USE moffload_struct_def

    USE lattice
    USE constant
    USE iso_c_binding

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    TYPE (latt)        LATT_CUR
    INTEGER            ION, IDIR
    ! local variable
    GDEF       :: CVAL
    INTEGER    :: NPL, N, LMBASE, LMMAXC, ISPINOR, ISPIRAL, K, KK, NT, LM
    COMPLEX(q) :: CMUL
    COMPLEX(q),TARGET :: CTMP(WDES1%NRPLWV)
    REAL(q),CONTIGUOUS,POINTER :: CTMP_REAL(:)
    GDEF,TARGET       :: C(WDES1%NPRO_TOT)
    REAL(q) :: SUMR, SUMI
    REAL(q) :: TMP(101,2)

    PROFILING_START('proj_der')

    NPL= WDES1%NGVECTOR
    CALL C_F_POINTER(C_LOC(CTMP),CTMP_REAL,[2*WDES1%NRPLWV])

!$ACC ENTER DATA CREATE(C,CTMP) __IF_ASYNC__
!$ACC KERNELS __IF_ASYNC__
    C=0
!$ACC END KERNELS

    ISPIRAL=1

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       NT=NONL_S%ITYP(ION)
       LMMAXC=NONL_S%LMMAX(NT)
       IF (LMMAXC==0) CYCLE
       LMBASE=NONL_S%LMBASE(ION)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)

       IF (IDIR==0) THEN
!$ACC PARALLEL LOOP PRIVATE(KK) PRESENT(W1,W1%CPTWFP,NONL_S,CTMP) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREC
          DO K=1,NPL
             KK=K+NPL*ISPINOR
             CTMP(K)=W1%CPTWFP(KK)*NONL_S%CREXP(K,ION)
          ENDDO
       ELSE
!$ACC PARALLEL LOOP PRIVATE(KK) PRESENT(W1,W1%CPTWFP,NONL_S,WDES1,LATT_CUR,CTMP) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREC
          DO K=1,NPL
             KK=K+NPL*ISPINOR
             CTMP(K)=W1%CPTWFP(KK)*NONL_S%CREXP(K,ION)* &
#ifdef shift_der_k
                 ((WDES1%IGX(K)+WDES1%VKPT(1))*LATT_CUR%B(IDIR,1)+ &
                  (WDES1%IGY(K)+WDES1%VKPT(2))*LATT_CUR%B(IDIR,2)+ &
                  (WDES1%IGZ(K)+WDES1%VKPT(3))*LATT_CUR%B(IDIR,3))
#else
                 ((WDES1%IGX(K))*LATT_CUR%B(IDIR,1)+ &
                  (WDES1%IGY(K))*LATT_CUR%B(IDIR,2)+ &
                  (WDES1%IGZ(K))*LATT_CUR%B(IDIR,3))
#endif
          ENDDO
       ENDIF

#if defined(do_loops) || defined(pro_loop) || defined(ACC_OFFLOAD)
!note(sm): it is not a good idea to parallelize this loop on gang-level as LMMAXC is too small
!          to saturate a GPU. Shouldn't the pattarn allow to replace the l_loop with a GEMV call?
!$ACC PARALLEL LOOP GANG PRIVATE(CMUL,CVAL) __IF_ASYNC__ &
!$ACC PRESENT(NONL_S,CTMP,C)
       l_loop: DO LM=1,LMMAXC
          IF (IDIR==0) THEN
             CMUL=NONL_S%CQFAK(LM,NT)
          ELSE
             CMUL=NONL_S%CQFAK(LM,NT)*CITPI
          ENDIF

          CVAL=0
!$ACC LOOP VECTOR REDUCTION(+:CVAL)
!DIR$ IVDEP
!OCL NOVREC
          DO K=1,NPL
             CVAL = CVAL+NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)*GREAL(CTMP(K)*CMUL)
          ENDDO

          C(LM+LMBASE)=CVAL
       ENDDO l_loop
#else
       CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, CTMP_REAL(1), 2, 0._q,  TMP(1,1), 1)
       CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, CTMP_REAL(2), 2, 0._q,  TMP(1,2), 1)
       IF (IDIR==0) THEN
           DO LM=1,LMMAXC
               SUMR = TMP(LM,1)
               SUMI = TMP(LM,2)
               C(LM+LMBASE)=GREAL(CMPLX( SUMR, SUMI, KIND=q)*NONL_S%CQFAK(LM,NT))
           ENDDO
       ELSE
           DO LM=1,LMMAXC
               SUMR = TMP(LM,1)
               SUMI = TMP(LM,2)
               C(LM+LMBASE)=GREAL(CMPLX( SUMR, SUMI, KIND=q)*NONL_S%CQFAK(LM,NT)*CITPI)
           ENDDO
       ENDIF
#endif
       IF (NONL_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor

    CALL DIS_PROJ(WDES1,C,W1%CPROJ(1))
!$ACC EXIT DATA DELETE(C,CTMP) __IF_ASYNC__

    PROFILING_STOP('proj_der')

    RETURN
  END SUBROUTINE PROJ_DER


!************************ SUBROUTINE VNLACC_DER ************************
!
!> Calculates the non local contribution of the first derivative
!> of the Hamiltonian upon moving one ion ION in the direction IDIR
!>
!> More specifically the contribution
!> ~~~
!>    | d p_i/ d R > (D_ij - epsilon Q_ij) c_j
!> ~~~
!> For IDIR=0 the derivative is replaced by | p_i >
!
!***********************************************************************

  SUBROUTINE VNLACC_DER(NONL_S, W1, &
       &     CDIJ,CQIJ,EVALUE,  CACC, LATT_CUR, ION, IDIR)
    USE lattice

    IMPLICIT NONE
    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER LMDIM
    REAL(q) EVALUE
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)   CACC(:)
    TYPE (latt)  LATT_CUR
    INTEGER      ION, IDIR
    ! work arrays
    GDEF         CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
    CALL OVERL1(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,1), CQIJ(1,1,1,1), EVALUE, W1%CPROJ(1), CRESUL(1))
    CALL VNLAC0_DER(NONL_S,W1%WDES1,CRESUL(1),CACC(1), LATT_CUR, ION, IDIR )
!$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__

    RETURN
  END SUBROUTINE VNLACC_DER

  SUBROUTINE VNLACC_DER_C(NONL_S, W1, &
       &     CDIJ, CQIJ,EVALUE,  CACC, LATT_CUR, ION, IDIR)
    USE lattice

    IMPLICIT NONE
    TYPE (nonl_struct) NONL_S
    TYPE (wavefun1)    W1
    INTEGER LMDIM
    COMPLEX(q)   EVALUE
    OVERLAP CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)   CACC(:)
    TYPE (latt)  LATT_CUR
    INTEGER      ION, IDIR
    ! work arrays
    GDEF         CRESUL(W1%WDES1%NPRO)

!$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
    CALL OVERL1_C(W1%WDES1,SIZE(CDIJ,1),CDIJ(1,1,1,1),CQIJ(1,1,1,1), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL VNLAC0_DER(NONL_S,W1%WDES1,CRESUL(1),CACC(1), LATT_CUR, ION, IDIR )
!$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__

    RETURN
  END SUBROUTINE VNLACC_DER_C


!************************ SUBROUTINE FORNL *****************************
!
!> Calculates the forces related to the non local pseudopotential
!>
!> The projection of the wavefunction onto the projector functions
!> must be stored in W\%CPROJ on entry
!>
!> Algorithm:
!> ~~~
!> ACC(G') = SUM(L,L',R)
!>         SUM(G)   QPROJ(G ,L)  EXP( iG  R) C(G)  iG *
!>         D(L,L') + Evalue  Q(L,L')
!>         SUM(GP)  QPROJ(GP,LP) EXP(-iGP R) C(GP)
!> ~~~
!
!***********************************************************************

  SUBROUTINE FORNL(NONL_S,WDES,W,LATT_CUR, LMDIM,CDIJ,CQIJ,EINL)

__GPU USE moffload

    USE lattice
    USE constant

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR
    TYPE (wavedes1)    WDES1         ! descriptor for one k-point

    DIMENSION EINL(3,NONL_S%NIONS)
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    ! work arrays
    REAL(q) ::  ENL(NONL_S%NIONS)
    GDEF :: CM
    GDEF,ALLOCATABLE :: CX(:) ,CY(:) ,CZ(:), C(:)
    GDEF,ALLOCATABLE :: CXL(:),CYL(:),CZL(:),CL(:)

    PROFILING_START('fornl')

    N =WDES%NPRO_TOT
    NL=WDES%NPRO

    ALLOCATE(CX(N),CY(N),CZ(N),C(N),CXL(NL),CYL(NL),CZL(NL),CL(NL))
GPACC !$ACC ENTER DATA CREATE(WDES1,CX,CY,CZ,C,CXL,CYL,CZL,CL,EINL,ENL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1,CX,CY,CZ,C,CXL,CYL,CZL,CL,EINL,ENL) IF(OFFLOAD_ON)


GPACC !$ACC KERNELS PRESENT(EINL,ENL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:EINL,ENL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    EINL=0._q
    ENL =0._q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
!=======================================================================
! loop over special points, and bands
!=======================================================================

    kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       NPL= WDES%NGVECTOR(NK)
       CALL SETWDES(WDES,WDES1,NK)

       CALL PHASE(WDES,NONL_S,NK)

       spin:   DO ISP=1,WDES%ISPIN
          band: DO N=1,WDES%NBANDS

             EVALUE=W%CELEN(N,NK,ISP)
             WEIGHT=W%FERWE(N,NK,ISP)*WDES%WTKPT(NK)*WDES%RSPIN
!=======================================================================
! first build up CX, CY, CZ for tables
!=======================================================================
GPACC !$ACC KERNELS PRESENT(CX,CY,CZ,C) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CX,CY,CZ,C) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
             CX=0
             CY=0
             CZ=0
             C=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,NONL_S,W,CX,CY,CZ) &
GPACC !$ACC PRIVATE(NT,LMBASE,CMUL,CTMP1,CTMP2,CTMP3,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:WDES,NONL_S,W,CX,CY,CZ) &
GPOMP !$OMP PRIVATE(NT,LMBASE,CMUL,CTMP1,CTMP2,CTMP3,CTMP) IF(OFFLOAD_ON)
             spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
__CPU           ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
                ions: DO NI=1,NONL_S%NIONS
__GPU              ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
                   NT=NONL_S%ITYP(NI)
                   LMBASE =NONL_S%LMBASE(NI)+ISPINOR *NONL_S%LMBASE(NONL_S%NIONS+1)

                   l_loop: DO LM=1,NONL_S%LMMAX(NT)
GPOMP !$OMP PARALLEL SHARED(WDES,NONL_S,W,CX,CY,CZ,CTMP,CTMP1,CTMP2,CTMP3) PRIVATE(K,KK,CVAL,CMUL) IF(OFFLOAD_ON)
                      CMUL=NONL_S%CQFAK(LM,NT)*CITPI

                      CTMP=0; CTMP1=0; CTMP2=0; CTMP3=0

                      !DIR$ IVDEP
                      !OCL NOVREC
GPACC !$ACC LOOP VECTOR REDUCTION(+:CTMP1,CTMP2,CTMP3,CTMP) PRIVATE(KK,CVAL)
GPOMP !$OMP DO REDUCTION(+:CTMP1,CTMP2,CTMP3,CTMP) PRIVATE(KK,CVAL)
                      DO K=1,NPL
                         KK=K+NPL*ISPINOR
                         CVAL=(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W%CPTWFP(KK,N,NK,ISP))*CMUL
                         ! new version calculate the projected wave function character on the fly
                         CTMP=CTMP+(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W%CPTWFP(KK,N,NK,ISP))*NONL_S%CQFAK(LM,NT)
#ifdef shift_der_k
                         ! new version add i k which contributes nothing to the energy derivative
                         CTMP1=CTMP1+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                         CTMP2=CTMP2+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                         CTMP3=CTMP3+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
#else
                         CTMP1=CTMP1+(WDES%IGX(K,NK))*CVAL
                         CTMP2=CTMP2+(WDES%IGY(K,NK))*CVAL
                         CTMP3=CTMP3+(WDES%IGZ(K,NK))*CVAL
#endif
                      ENDDO
GPOMP !$OMP END PARALLEL
                      C (LM+LMBASE)=C (LM+LMBASE)+CTMP
                      CX(LM+LMBASE)=CX(LM+LMBASE)+CTMP1
                      CY(LM+LMBASE)=CY(LM+LMBASE)+CTMP2
                      CZ(LM+LMBASE)=CZ(LM+LMBASE)+CTMP3
                   ENDDO l_loop
                ENDDO ions
             ENDDO spinor

             CALL DIS_PROJ(WDES1,C ,CL )
             CALL DIS_PROJ(WDES1,CX,CXL)
             CALL DIS_PROJ(WDES1,CY,CYL)
             CALL DIS_PROJ(WDES1,CZ,CZL)
!=======================================================================
! sum up local contributions
! calculate SUM_LP  D(LP,L)-E Q(LP,L) * C(LP)
!=======================================================================
GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(WDES,NONL_S,W,CX,CY,CZ) &
GPACC !$ACC PRIVATE(NT,LMBASE,LMBASE_,NIP,CTMP,CTMP1,CTMP2,CTMP3) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:WDES,NONL_S,W,CX,CY,CZ) &
GPOMP !$OMP PRIVATE(NT,LMBASE,LMBASE_,NIP,CTMP,CTMP1,CTMP2,CTMP3) THREAD_LIMIT(64) IF(OFFLOAD_ON)
             spinor2: DO ISPINOR=0,WDES%NRSPINORS-1
                DO ISPINOR_=0,WDES%NRSPINORS-1

                   ions2: DO NI=1,WDES%NIONS
GPOMP !$OMP PARALLEL SHARED(WDES,NONL_S,W,CXL,CYL,CZL,CL,CTMP,CTMP1,CTMP2,CTMP3) PRIVATE(CM) IF(OFFLOAD_ON)
                      NT=WDES%ITYP(NI)
                      LMBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                      LMBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2

                      NIP=NI_GLOBAL(NI, WDES%COMM_INB) !  local storage index

                      CTMP=0; CTMP1=0; CTMP2=0; CTMP3=0

GPACC !$ACC LOOP COLLAPSE(2) REDUCTION(+:CTMP) REDUCTION(-:CTMP1,CTMP2,CTMP3) PRIVATE(CM)
GPOMP !$OMP DO COLLAPSE(2) REDUCTION(+:CTMP,CTMP1,CTMP2,CTMP3)
                      l_loop2: DO LM=1,WDES%LMMAX(NT)
                         DO LMP=1,WDES%LMMAX(NT)
                            CM=(CDIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
                               ! new version use the calculated projected wave function character
                               GCONJG(CL(LM+LMBASE))
!                                 GCONJG(W%CPROJ(LM+LMBASE,N,NK,ISP))
                            CTMP=CTMP + (W%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
                            CTMP1=CTMP1-(2*WEIGHT)*(CXL(LMP+LMBASE_)*CM)
                            CTMP2=CTMP2-(2*WEIGHT)*(CYL(LMP+LMBASE_)*CM)
                            CTMP3=CTMP3-(2*WEIGHT)*(CZL(LMP+LMBASE_)*CM)
                         ENDDO
                      ENDDO l_loop2
GPOMP !$OMP END PARALLEL
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                      ENL (NIP)  =ENL (NIP)  +REAL(CTMP, KIND=q)
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                      EINL(1,NIP)=EINL(1,NIP)+REAL(CTMP1,KIND=q)
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                      EINL(2,NIP)=EINL(2,NIP)+REAL(CTMP2,KIND=q)
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                      EINL(3,NIP)=EINL(3,NIP)+REAL(CTMP3,KIND=q)
                   ENDDO ions2
                ENDDO
             ENDDO spinor2
!=======================================================================
          ENDDO band
       ENDDO spin
    ENDDO kpoint
#ifdef OFFLOADING
GPACC !$ACC EXIT DATA COPYOUT(EINL,ENL) DELETE(CX,CY,CZ,C,CXL,CYL,CZL,CL) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:EINL,ENL) MAP(RELEASE:CX,CY,CZ,C,CXL,CYL,CZL,CL) IF(OFFLOAD_ON)
    CALL DELETE_TYPED_VAR(WDES1)
#endif
!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
    CALLMPI( M_sum_d(WDES%COMM, EINL(1,1),NONL_S%NIONS*3))
    CALLMPI( M_sum_d(WDES%COMM, ENL (1)  ,NONL_S%NIONS))
    CALL  DIRKAR(NONL_S%NIONS,EINL,LATT_CUR%B)
    DEALLOCATE(CX,CY,CZ,C,CXL,CYL,CZL,CL)

    PROFILING_STOP('fornl')

    RETURN
  END SUBROUTINE FORNL


!************************ SUBROUTINE STRENL *****************************
!
!> Calculates non-local contributions to stress
!>
!> Easiest to implement and definitly quite fast, is an approach based
!> on finite differences
!
!***********************************************************************

  SUBROUTINE STRENL(GRID,NONL_S,P,W,WDES,LATT_CUR, &
       LMDIM,CDIJ,CQIJ, ISIF,FNLSIF)

__GPU USE moffload

    USE pseudo_struct_def
    USE mgrid_struct_def
    USE poscar_struct_def
    USE wave_struct_def

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR

    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

    REAL(q) FNLSIF(3,3)
    INTEGER LMDIM, ISIF

    PROFILING_START('strenl')

#ifdef OFFLOADING
    IF (OFFLOAD_ON) THEN
       CALL STRENL_ACC(GRID,NONL_S,P,W,WDES,LATT_CUR,LMDIM,CDIJ,CQIJ,ISIF,FNLSIF)
       PROFILING_STOP('strenl')
       RETURN
    ENDIF
#endif

    CALL STRENL_(GRID,NONL_S,P,W,WDES,LATT_CUR,LMDIM,CDIJ,CQIJ,ISIF,FNLSIF)

    PROFILING_STOP('strenl')

  END SUBROUTINE STRENL


  SUBROUTINE STRENL_(GRID,NONL_S,P,W,WDES,LATT_CUR, &
       LMDIM,CDIJ,CQIJ, ISIF,FNLSIF)

__GPU USE moffload

    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1
    TYPE (wavespin)    W
    TYPE (wavefun1)    W1
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR,LATT_FIN

    DIMENSION FNLSIF(3,3)
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    ! work arrays
    GDEF,ALLOCATABLE,TARGET ::  CWORK(:)
    REAL(q) FNLSIFTMP

    DIS=fd_displacement
!TEST which step should be used in the finite differences
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
    ALLOCATE(CWORK(WDES%NPROD))

    FNLSIF=0
GPACC !$ACC ENTER DATA CREATE(CWORK,WDES1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:WDES1,CWORK) IF(OFFLOAD_ON)

    DO IDIR=1,3
       DO JDIR=1,3
          FNLSIFTMP=0
!=======================================================================
! use central differences to calculate the stress
! set up QPROJ so that central differences can be  calculated
!=======================================================================
          LATT_FIN=LATT_CUR
          IF (ISIF==1) THEN
             !  only isotrop pressure
             DO I=1,3; DO J=1,3
                LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
             ENDDO; ENDDO
          ELSE
             !  all directions
             DO I=1,3
                LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN)

          CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO=-1)

          LATT_FIN=LATT_CUR
          IF (ISIF==1) THEN
             !  only isotrop pressure
             DO I=1,3; DO J=1,3
                LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
             ENDDO; ENDDO
          ELSE
             !  all directions
             DO I=1,3
                LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN)

          CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO=0)

!=======================================================================
! loop over all k-points spin and bands
!=======================================================================
          kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
             IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
             CALL PHASE(WDES,NONL_S,NK)
             CALL SETWDES(WDES,WDES1,NK)

             spin: DO ISP=1,WDES%ISPIN
                band: DO N=1,WDES%NBANDS
                   CALL SETWAV(W,W1,WDES1,N,ISP); W1%CPROJ => CWORK
#ifdef OFFLOADING
                   CALL COPYIN_TYPED_VAR(W1)
#endif
                   CALL PROJ1(NONL_S,WDES1,W1)  ! calculate W1%CPROJ (linked to CWORK)

                   WEIGHT=WDES%WTKPT(NK)*WDES%RSPIN*W1%FERWE
                   EVALUE=W1%CELEN

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(WDES,CWORK,W,W%CPROJ,CDIJ,CQIJ) &
GPACC !$ACC PRIVATE(LBASE,LBASE_,NT) REDUCTION(+:FNLSIFTMP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) MAP(PRESENT,TO:WDES,CWORK,W,W%CPROJ,CDIJ,CQIJ) &
GPOMP !$OMP PRIVATE(LBASE,LBASE_,NT) THREAD_LIMIT(64) REDUCTION(+:FNLSIFTMP) IF(OFFLOAD_ON)
                   spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
                      DO ISPINOR_=0,WDES%NRSPINORS-1

                         DO NI=1,WDES%NIONS
                            NT=WDES%ITYP(NI)
                            LBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                            LBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2

GPACC !$ACC LOOP COLLAPSE(2) REDUCTION(+:FNLSIFTMP)
GPOMP !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:FNLSIFTMP) IF(OFFLOAD_ON)
                            DO L=1,WDES%LMMAX(NT)
__CPU !$OMP SIMD REDUCTION(+:FNLSIFTMP)
                               DO LP=1,WDES%LMMAX(NT)
                                  FNLSIFTMP = FNLSIFTMP+WEIGHT*CWORK(LBASE_+LP)* &
                                       GCONJG(W%CPROJ(LBASE+L,N,NK,ISP))* &
                                       (CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO spinor
#ifdef OFFLOADING
                   CALL DELETE_TYPED_VAR(W1)
#endif
                ENDDO band
             ENDDO spin
          ENDDO kpoint
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          FNLSIF(IDIR,JDIR)=FNLSIFTMP
          !
          !  only isotrop pressure finish now
          !
          IF (ISIF==1) THEN
             FNLSIF(2,2)= FNLSIF(1,1)
             FNLSIF(3,3)= FNLSIF(1,1)
             GOTO 310  ! terminate (not very clean but who cares)
          ENDIF
!=======================================================================
! next direction
!=======================================================================
       ENDDO
    ENDDO

310 CONTINUE
#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK) IF(OFFLOAD_ON)
    CALL DELETE_TYPED_VAR(WDES1)
#endif
    CALLMPI( M_sum_d(WDES%COMM, FNLSIF, 9))

    FNLSIF=FNLSIF/DIS

#ifdef debug
    WRITE(*,'(E10.3,3E14.7)')DIS,((FNLSIF(I,J),I=1,3),J=1,3)
!    IF (DIS>1E-10) GOTO 1000
!    GOTO 1000
#endif
!=======================================================================
! recalculate the projection operators
! (the array was used as a workspace)
!=======================================================================
    IZERO=1
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,  IZERO)

    DEALLOCATE(CWORK)

    RETURN
  END SUBROUTINE STRENL_

#ifdef OFFLOADING
  SUBROUTINE STRENL_ACC(GRID,NONL_S,P,W,WDES,LATT_CUR,LMDIM,CDIJ,CQIJ,ISIF,FNLSIF)

__GPU USE moffload

    USE pseudo
    USE mpimy
    USE mgrid
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(NONL_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR

    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

    INTEGER LMDIM, ISIF
    REAL(q) FNLSIF(3,3)

    ! local variables
    TYPE (latt) :: LATT_FIN
    GDEF, ALLOCATABLE :: CPROJ(:,:)
    LOGICAL, ALLOCATABLE :: LDO(:)

    REAL(q) FNLSIFTMP, DIS, WEIGHT, EVALUE
    INTEGER :: IDIR, JDIR, I,J, NK, ISP, N, NSIM, NSIM_ACT, &
               ISPINOR, ISPINOR_, NB, NT, LBASE, LBASE_, NI, NIP, L, LP

    DIS=fd_displacement
!TEST which step should be used in the finite differences
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST

    ! same blocking as in EDDAV
    NSIM = 2*WDES%NSIM

    ALLOCATE(CPROJ(WDES%NPROD,NSIM),LDO(NSIM))
GPACC !$ACC ENTER DATA CREATE(CPROJ,WDES1,FNLSIF,FNLSIFTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ,WDES1,FNLSIF,FNLSIFTMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(FNLSIF) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNLSIF) IF(OFFLOAD_ON)
    FNLSIF = 0
GPACC !$ACC END KERNELS
GPOMP !$OMP END TARGET

    LDO = .TRUE.

    idr: DO IDIR=1,3
       jdr: DO JDIR=1,3
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC !$ACC KERNELS PRESENT(FNLSIFTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNLSIFTMP) IF(OFFLOAD_ON)
          FNLSIFTMP=0
GPACC !$ACC END KERNELS
GPOMP !$OMP END TARGET
!=======================================================================
! use central differences to calculate the stress
! set up QPROJ so that central differences can be  calculated
!=======================================================================

          LATT_FIN=LATT_CUR
          IF (ISIF==1) THEN
             !  only isotrop pressure
             DO I=1,3; DO J=1,3
                LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
             ENDDO; ENDDO
          ELSE
             !  all directions
             DO I=1,3
                LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN)

          CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO=-1)

          LATT_FIN=LATT_CUR
          IF (ISIF==1) THEN
             !  only isotrop pressure
             DO I=1,3; DO J=1,3
                LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
             ENDDO; ENDDO
          ELSE
             !  all directions
             DO I=1,3
                LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN)

          CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,IZERO=0)

!=======================================================================
! loop over all k-points spin and bands
!=======================================================================
          kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
             IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
             CALL PHASE(WDES,NONL_S,NK)
             CALL SETWDES(WDES,WDES1,NK)

             spin: DO ISP=1,WDES%ISPIN
                strip: DO N=1,WDES%NBANDS,NSIM
                   NSIM_ACT=MIN(NSIM,WDES%NBANDS-N+1)

                   CALL PROJ0_MU(NONL_S,WDES1,W%CW(1,N,NK,ISP),WDES%NRPLWV,NSIM_ACT,CPROJ(1,1),LDO(1))

                   PROFILING_START('strenl|contract')

GPACC !$ACC PARALLEL LOOP COLLAPSE(4) GANG PRESENT(CPROJ,W,WDES,CDIJ,CQIJ,FNLSIFTMP) &
GPACC !$ACC PRIVATE(LBASE,LBASE_,NT,EVALUE,WEIGHT) REDUCTION(+:FNLSIFTMP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(4) MAP(PRESENT,TO:CPROJ,W,WDES,CDIJ,CQIJ,FNLSIFTMP) &
GPOMP !$OMP PRIVATE(LBASE,LBASE_,NT,EVALUE,WEIGHT) THREAD_LIMIT(128) REDUCTION(+:FNLSIFTMP) IF(OFFLOAD_ON)
                   spinor: DO ISPINOR=0,WDES%NRSPINORS-1
                      DO ISPINOR_=0,WDES%NRSPINORS-1

                         bands: DO NB=N,N+NSIM_ACT-1

                            ions: DO NI=1,WDES%NIONS
                               NT=WDES%ITYP(NI)
                               LBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                               LBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2

                               EVALUE=W%CELEN(NB,NK,ISP)
                               WEIGHT=WDES%WTKPT(NK)*WDES%RSPIN*W%FERWE(NB,NK,ISP)

GPACC !$ACC LOOP COLLAPSE(2) REDUCTION(+:FNLSIFTMP)
GPOMP !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:FNLSIFTMP) IF(OFFLOAD_ON)
                               DO L=1 ,WDES%LMMAX(NT)
                                  DO LP=1,WDES%LMMAX(NT)
                                     FNLSIFTMP = FNLSIFTMP+WEIGHT*CPROJ(LBASE_+LP,NB-N+1)* &
                                          GCONJG(W%CPROJ(LBASE+L,NB,NK,ISP))* &
                                          (CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
                                  ENDDO
                               ENDDO

                            ENDDO ions

                         ENDDO bands

                      ENDDO
                   ENDDO spinor

                   PROFILING_STOP('strenl|contract')

                ENDDO strip
             ENDDO spin
          ENDDO kpoint

GPACC !$ACC KERNELS PRESENT(FNLSIF,FNLSIFTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:FNLSIF,FNLSIFTMP) IF(OFFLOAD_ON)
          FNLSIF(IDIR,JDIR)=FNLSIFTMP
GPACC !$ACC END KERNELS
GPOMP !$OMP END TARGET

          !  if only isotropic pressure finish now
          IF (ISIF==1) EXIT idr
!=======================================================================
! next direction
!=======================================================================
       ENDDO jdr
    ENDDO idr

#ifdef OFFLOADING
    CALL DELETE_TYPED_VAR(WDES1)
GPACC !$ACC EXIT DATA DELETE(CPROJ,FNLSIFTMP) COPYOUT(FNLSIF) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:FNLSIF) MAP(RELEASE:CPROJ,FNLSIFTMP) IF(OFFLOAD_ON)
#endif

    CALLMPI( M_sum_d(WDES%COMM, FNLSIF, 9))

    !  if only isotropic pressure
    IF (ISIF==1) THEN
       FNLSIF(2,2)= FNLSIF(1,1)
       FNLSIF(3,3)= FNLSIF(1,1)
    ENDIF

    FNLSIF=FNLSIF/DIS

#ifdef debug
    WRITE(*,'(E10.3,3E14.7)')DIS,((FNLSIF(I,J),I=1,3),J=1,3)
!    IF (DIS>1E-10) GOTO 1000
!    GOTO 1000
#endif
!=======================================================================
! recalculate the projection operators
! (the array was used as a workspace)
!=======================================================================
    CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,IZERO=1)

    DEALLOCATE(CPROJ,LDO)

    RETURN
  END SUBROUTINE STRENL_ACC
#endif // OFFLOADING

END MODULE nonl

!************************ SUBROUTINE VNLAC0 ****************************
!
!> Calculates a linear combination of projection operators in reciprocal
!> space
!>
!> The result is *added* to CACC
!
!***********************************************************************

  SUBROUTINE VNLAC0(NONL_S,WDES1,CPROJ_LOC,CACC)

__GPU USE moffload

    USE nonl_struct_def
    USE wave_struct_def

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    GDEF    CPROJ_LOC(WDES1%NPRO)
    COMPLEX(q) CACC(WDES1%NRPLWV)

#ifdef OFFLOADING
    IF (OFFLOAD_ON) THEN
       CALL VNLAC0_ACC(NONL_S,WDES1,CPROJ_LOC,CACC)
       RETURN
    ENDIF
#endif
    CALL VNLAC0_(NONL_S,WDES1,CPROJ_LOC,CACC)
    RETURN
  END SUBROUTINE VNLAC0

  SUBROUTINE VNLAC0_(NONL_S,WDES1,CPROJ_LOC,CACC)
    USE nonl_struct_def
    USE wave_struct_def

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    GDEF    CPROJ_LOC(WDES1%NPRO)
    COMPLEX(q) CACC(WDES1%NRPLWV)
  ! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, LMMAXC, NI, LM, K
    COMPLEX(q) :: CTMP
#ifdef avoidalloc
    REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(q),ALLOCATABLE :: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:)

    ALLOCATE(WORK(WDES1%NGVECTOR*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif
    PROFILING_START('vnlac0')

    IF (WDES1%NK /= NONL_S%NK) THEN
       WRITE(*,*) 'internal error in VNLAC0: PHASE not properly set up',WDES1%NK, NONL_S%NK
    ENDIF

    NPL=WDES1%NGVECTOR

    ! merge projected wavefunctions from all nodes
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
!=======================================================================
! performe loop over ions
!=======================================================================
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NI,NT,LMMAXC,LMBASE,WORK,CTMP,TMP,LM) REDUCTION(+:CACC)
       ion: DO NI=1,NONL_S%NIONS
          NT=NONL_S%ITYP(NI)
          LMMAXC=NONL_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion
          LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
          DO LM=1,LMMAXC
             CTMP= CPROJ(LMBASE+LM)*CONJG(NONL_S%CQFAK(LM,NT))
             TMP(LM,1)= REAL( CTMP ,KIND=q)
             TMP(LM,2)= AIMAG(CTMP)
          ENDDO
#if defined(do_loops)
!DIR$ IVDEP
!OCL NOVREC
          DO K=1,2*NPL
             WORK(K)=0
          ENDDO
!=======================================================================
! call BLAS to calculate SUM(LM=1,LMMAX) QPROJ(K,LM) * TMP(LM)
!=======================================================================
          DO LM=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
             DO K=1,NPL
                WORK(K)    =WORK(K)    +TMP(LM,1)*NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
                WORK(K+NPL)=WORK(K+NPL)+TMP(LM,2)*NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
             ENDDO
          ENDDO
#else
          PROFILING_START('vnlac0|dgemv')

          CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, TMP(1,1) , 1 , 0._q , WORK(1), 1)
          CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
               WDES1%NGDIM, TMP(1,2) , 1 , 0._q , WORK(1+NPL), 1)

          PROFILING_STOP('vnlac0|dgemv')
#endif
!=======================================================================
! add acceleration from this ion to CACC
!=======================================================================
          CALL WORK_MUL_CREXP( NPL, WORK(1), NONL_S%CREXP(1,NI), CACC(1+NPL*ISPINOR))

       ENDDO ion
    ENDDO spinor

#ifndef avoidalloc
    DEALLOCATE(WORK,TMP,CPROJ)
#endif

    PROFILING_STOP('vnlac0')

    RETURN
  END SUBROUTINE VNLAC0_

#ifdef OFFLOADING
  SUBROUTINE VNLAC0_ACC(NONL_S,WDES1,CPROJ_LOC,CACC)
__GPU USE moffload

    USE nonl_struct_def
    USE wave_struct_def

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    GDEF    CPROJ_LOC(WDES1%NPRO)
    COMPLEX(q) CACC(WDES1%NRPLWV)
  ! local
    INTEGER :: NPL, LMBASE, ISPIRAL, ISPINOR, NT, NI, LM, K
    COMPLEX(q) :: CTMP, CWORK
    GDEF,ALLOCATABLE    :: CPROJ(:)

    ALLOCATE(CPROJ(WDES1%NPRO_TOT))
    PROFILING_START('vnlac0_acc')

    IF (WDES1%NK /= NONL_S%NK) THEN
       WRITE(*,*) 'internal error in VNLAC0: PHASE not properly set up',WDES1%NK, NONL_S%NK
    ENDIF

    NPL=WDES1%NGVECTOR

GPACC !$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ) IF(OFFLOAD_ON)

    ! merge projected wavefunctions from all nodes
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR_LENGTH(32) PRESENT(CACC,NONL_S,CPROJ) PRIVATE(CTMP) __IF_ASYNC__

! CTMP and CWORK are both scalars that are distinct for each team, make them explicitly private
! for very large unit cells it might be better to increase THREAD_LIMIT here, but I assume that
! this is only used in small systems since for large systems lreal should be used.
GPOMP !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(PRESENT,TO:CACC,NONL_S,CPROJ) &
GPOMP !$OMP THREAD_LIMIT(64) PRIVATE(CTMP,CWORK) IF(OFFLOAD_ON)
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO K=1,NPL
          CTMP = 0
GPACC !$ACC LOOP VECTOR PRIVATE(CWORK,ISPIRAL,NT,LMBASE,LM) REDUCTION(+:CTMP)

! each team now spawns on the gpu x threads working in parallel on CTMP AND CWORK
GPOMP !$OMP PARALLEL DO PRIVATE(ISPIRAL,NT,LMBASE,LM,CWORK) REDUCTION(+:CTMP) IF(OFFLOAD_ON)
          ion: DO NI=1,NONL_S%NIONS
             ISPIRAL=1; IF (NONL_S%LSPIRAL) ISPIRAL=ISPINOR+1
             NT=NONL_S%ITYP(NI)
             LMBASE=NONL_S%LMBASE(NI)+ISPINOR*NONL_S%LMBASE(NONL_S%NIONS+1)
             CWORK=0
!$ACC LOOP SEQ
             DO LM=1,NONL_S%LMMAX(NT)
                CWORK=CWORK+CPROJ(LMBASE+LM)*CONJG(NONL_S%CQFAK(LM,NT))*NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
             ENDDO
             CTMP=CTMP+CWORK *CONJG(NONL_S%CREXP(K,NI))
          ENDDO ion
          CACC(NPL*ISPINOR+K)=CACC(NPL*ISPINOR+K)+CTMP
       ENDDO
    ENDDO spinor

GPACC !$ACC EXIT DATA DELETE(CPROJ) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CPROJ) IF(OFFLOAD_ON)
    DEALLOCATE(CPROJ)

    PROFILING_STOP('vnlac0_acc')
    RETURN
  END SUBROUTINE VNLAC0_ACC
#endif

!************************ SUBROUTINE VNLAC0_DER *************************
!
!> Calculates the non local contribution of the first derivative
!> of the Hamiltonian upon moving one ion ION in the direction IDIR
!>
!> More specifically the contribution
!> ~~~
!>    | d p_i/ d R > c_i
!> ~~~
!> For IDIR=0 the derivative is replaced by | p_i >
!>
!> A term (i k) is included as well to "center" the derivate at the
!> ion
!
!***********************************************************************

  SUBROUTINE VNLAC0_DER(NONL_S,WDES1,CPROJ_LOC,CACC, LATT_CUR, ION, IDIR)

GPACC USE moffload_struct_def

    USE nonl_struct_def
    USE wave_struct_def
    USE lattice
    USE constant
    USE string, ONLY: str
    USE tutor, ONLY: vtutor

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes1)    WDES1
    COMPLEX(q) CACC(WDES1%NRPLWV)
    GDEF       CPROJ_LOC(WDES1%NPRO)
    TYPE (latt)  LATT_CUR
    INTEGER      ION, IDIR

    ! work arrays
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
    INTEGER :: NPL, LMBASE, LMMAXC, ISPINOR,  ISPIRAL, K, KK, NI, NIS, NT, LM, L
    COMPLEX(q) CMUL
    COMPLEX(q) CACC_ADD(WDES1%NRPLWV)

    PROFILING_START('vnlac0_der')

    IF (WDES1%NK /= NONL_S%NK) &
       CALL vtutor%bug("VNLAC0_DER: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONL_S%NK),__FILE__,__LINE__)

    NPL=WDES1%NGVECTOR

!$ACC ENTER DATA CREATE(CPROJ,CACC_ADD) __IF_ASYNC__
    ! merge projected wavefunctions from all nodes
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))
!=======================================================================
! performe loop over ions
!=======================================================================
    LMBASE= 0
    ISPIRAL = 1

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       NIS=1
       typ: DO NT=1,NONL_S%NTYP
          LMMAXC=NONL_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ions: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
!note(ar): the following IF is cuts down all the preceeding loops to a single pass (no loops at all)
             IF (NI==ION) THEN
!$ACC KERNELS PRESENT(CACC_ADD) __IF_ASYNC__
                CACC_ADD=0
!$ACC END KERNELS
!$ACC PARALLEL LOOP SEQ PRIVATE(CMUL) PRESENT(CACC_ADD,NONL_S,CPROJ) __IF_ASYNC__
                l_loop: DO LM=1,LMMAXC
                   CMUL=CPROJ(LMBASE+LM)*CONJG(NONL_S%CQFAK(LM,NT))
!$ACC LOOP VECTOR GANG
!DIR$ IVDEP
!OCL NOVREC
                   DO K=1,NPL
                      CACC_ADD(K)=CACC_ADD(K)+NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)*CMUL* &
                           (CONJG(NONL_S%CREXP(K,NI)))
                   ENDDO
                ENDDO l_loop
!=======================================================================
! add acceleration from this ion
!=======================================================================

                IF (IDIR==0) THEN
!$ACC PARALLEL LOOP PRESENT(CACC,CACC_ADD) PRIVATE(KK) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREC
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CACC(KK)=CACC(KK)+CACC_ADD(K)
                   ENDDO
                ELSE
!$ACC PARALLEL LOOP PRESENT(CACC,CACC_ADD,WDES1,LATT_CUR) PRIVATE(KK) __IF_ASYNC__
!DIR$ IVDEP
!OCL NOVREC
                   DO K=1,NPL
                      KK=K+NPL*ISPINOR
                      CACC(KK)= CACC(KK)+CACC_ADD(K)*CONJG(CITPI)* &
#ifdef shift_der_k
                          ((WDES1%IGX(K)+WDES1%VKPT(1))*LATT_CUR%B(IDIR,1)+ &
                           (WDES1%IGY(K)+WDES1%VKPT(2))*LATT_CUR%B(IDIR,2)+ &
                           (WDES1%IGZ(K)+WDES1%VKPT(3))*LATT_CUR%B(IDIR,3))
#else
                          ((WDES1%IGX(K))*LATT_CUR%B(IDIR,1)+ &
                           (WDES1%IGY(K))*LATT_CUR%B(IDIR,2)+ &
                           (WDES1%IGZ(K))*LATT_CUR%B(IDIR,3))
#endif
                   ENDDO
                ENDIF

             ENDIF
             LMBASE= LMMAXC+LMBASE
          ENDDO ions

600       NIS = NIS+NONL_S%NITYP(NT)
       ENDDO typ
       IF (NONL_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$ACC EXIT DATA DELETE(CPROJ,CACC_ADD) __IF_ASYNC__

    PROFILING_STOP('vnlac0_der')

    RETURN
  END SUBROUTINE VNLAC0_DER

#ifdef OFFLOADING
  SUBROUTINE PROJ0_MU(NONL_S, WDES1, CW, LD, NSIM, CPROJ_LOC, LDO)

__GPU USE moffload

    USE nonl_struct_def
    USE wave_struct_def

    IMPLICIT NONE

    TYPE (nonl_struct) :: NONL_S
    TYPE (wavedes1) :: WDES1
    INTEGER :: LD                       !< leading dimension of CACC
    INTEGER :: NSIM                     !< do NSIM bands at a time
    LOGICAL :: LDO(NSIM)                !< which bands are included
    COMPLEX(q) :: CW(LD,NSIM)           !< result in reciprocal space
    GDEF :: CPROJ_LOC(WDES1%NPROD,NSIM) !< wave function character

    ! local variables
    INTEGER :: NP, NGVECTOR, ISPINOR, IBLOCK, BLOCKSIZE_ACT, NG, NI, NT, LMBASE, LM, ISPIRAL, NK, NIONS, NPRO
    INTEGER, PARAMETER :: BLOCKSIZE=2048

    GDEF :: GADD

#ifdef avoidalloc
    COMPLEX(q) :: CPR(MIN(BLOCKSIZE,WDES1%NGVECTOR),WDES1%NPRO_TOT/WDES1%NRSPINORS)
    GDEF :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    COMPLEX(q), ALLOCATABLE :: CPR(:,:)
    GDEF, ALLOCATABLE :: CPROJ(:,:)
#endif

    PROFILING_START('proj0_mu')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),CPR(MIN(BLOCKSIZE,WDES1%NGVECTOR),WDES1%NPRO_TOT/WDES1%NRSPINORS)
#endif
GPACC !$ACC ENTER DATA CREATE(CPROJ,CPR) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CPROJ,CPR) IF(OFFLOAD_ON)

    NGVECTOR = WDES1%NGVECTOR
    NPRO = WDES1%NPRO_TOT / WDES1%NRSPINORS
    NIONS = NONL_S%NIONS
    NK = WDES1%NK

    spinor: DO ISPINOR = 0, WDES1%NRSPINORS-1

       iblk: DO IBLOCK = 1, NGVECTOR, BLOCKSIZE
          BLOCKSIZE_ACT = MIN(NGVECTOR-IBLOCK+1,BLOCKSIZE)

          PROFILING_START('proj0_mu|cpr')

GPACC !$ACC PARALLEL LOOP GANG PRESENT(NONL_S,CPR) PRIVATE(ISPIRAL,NT,LMBASE) VECTOR_LENGTH(1024) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:NONL_S,CPR) PRIVATE(ISPIRAL,NT,LMBASE) IF(OFFLOAD_ON)
          ion: DO NI = 1, NIONS
              ISPIRAL = 1; IF (NONL_S%LSPIRAL) ISPIRAL = ISPINOR + 1
              NT = NONL_S%ITYP(NI)
              LMBASE = NONL_S%LMBASE(NI)

GPACC !$ACC LOOP WORKER
GPOMP !$OMP PARALLEL DO COLLAPSE(2) IF(OFFLOAD_ON)
              DO LM = 1, NONL_S%LMMAX(NT)
GPACC !$ACC LOOP VECTOR
                 DO NG = IBLOCK, IBLOCK + BLOCKSIZE_ACT -1
                    CPR(NG-IBLOCK+1,LMBASE+LM) = CONJG(NONL_S%CQFAK(LM,NT)) * NONL_S%QPROJ(NG,LM,NT,NK,ISPIRAL)*CONJG(NONL_S%CREXP(NG,NI))
                 ENDDO
              ENDDO

          ENDDO ion

          PROFILING_STOP('proj0_mu|cpr')

          PROFILING_START('proj0_mu|gemm')

          IF (IBLOCK==1) THEN
           GADD = zero
          ELSE
           GADD = one
          ENDIF

          CALL WGGEMM(trans, 'N', NPRO, NSIM, m_ BLOCKSIZE_ACT, &
                      one, CPR(1,1), m_ SIZE(CPR,1), &
                      CW(IBLOCK+ISPINOR*NGVECTOR,1), m_ WDES1%NRPLWV, &
                      GADD, CPROJ(1+ISPINOR*NPRO,1), WDES1%NPRO_TOT)

          PROFILING_STOP('proj0_mu|gemm')

       ENDDO iblk

    ENDDO spinor

    ! distribute the projected wavefunctions to nodes
    DO NP = 1, NSIM
       IF (LDO(NP)) THEN
          CALL DIS_PROJ(WDES1, CPROJ(1,NP), CPROJ_LOC(1,NP))
       ENDIF
    ENDDO

GPACC !$ACC EXIT DATA DELETE(CPROJ,CPR) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CPROJ,CPR) IF(OFFLOAD_ON)
#ifndef avoidalloc
    DEALLOCATE(CPROJ,CPR)
#endif

    PROFILING_STOP('proj0_mu')

    RETURN
  END SUBROUTINE PROJ0_MU
#endif // OFFLOADING

!***********************************************************************
!
!> small f77 helper routine to multiply with phase factor and divide
!> into real and imaginary part
!
!***********************************************************************

  SUBROUTINE CREXP_MUL_WAVE( NPL, CREXP, CW, WORK)
!$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
    INTEGER NPL
    COMPLEX(q) :: CREXP(NPL), CW(NPL)
    REAL(q)    :: WORK(2*NPL)
    COMPLEX(q) :: CTMP
  ! local
    INTEGER K

NOACC PROFILING_START('crexp_mul_wave')

!$ACC LOOP VECTOR
NOACC !$OMP PARALLEL DO SIMD SCHEDULE(__simd_static) DEFAULT(SHARED) PRIVATE(CTMP)
    DO K=1,NPL
       CTMP=    CREXP(K) *CW(K)
       WORK(K)    = REAL( CTMP ,KIND=q)
       WORK(K+NPL)= AIMAG(CTMP)
    ENDDO
NOACC !$OMP END PARALLEL DO SIMD

NOACC PROFILING_STOP('crexp_mul_wave')

  END SUBROUTINE CREXP_MUL_WAVE


  SUBROUTINE WORK_MUL_CREXP( NPL, WORK, CREXP, CW)
!$ACC ROUTINE VECTOR
    USE prec

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

    IMPLICIT NONE

    INTEGER NPL
    COMPLEX(q) :: CREXP(NPL), CW(NPL)
    REAL(q)    :: WORK(2*NPL)
    COMPLEX(q) :: CTMP
  ! local
    INTEGER K

NOACC PROFILING_START('work_mul_crexp')

!$ACC LOOP VECTOR PRIVATE(CTMP)
NOACC !$OMP PARALLEL DO SIMD SCHEDULE(__simd_static) DEFAULT(SHARED) PRIVATE(CTMP)
    DO K=1,NPL
#ifndef ACC_OFFLOAD
       CW(K)=CW(K)+ CMPLX( WORK(K) , WORK(K+NPL) ,KIND=q) *CONJG(CREXP(K))
#else // ACC_OFFLOAD
       CTMP=CMPLX( WORK(K) , WORK(K+NPL) ,KIND=q) *CONJG(CREXP(K))
       CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CW(K),CTMP)
#endif // ACC_OFFLOAD
    ENDDO
NOACC !$OMP END PARALLEL DO SIMD

NOACC PROFILING_STOP('work_mul_crexp')

  END SUBROUTINE WORK_MUL_CREXP
