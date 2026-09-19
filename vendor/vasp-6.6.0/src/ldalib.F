#include "symbol.inc"

!************************ MODULE LDALIB ***********************
!>
!> This module is the library of subroutines of LDA functionals
!>
!**************************************************************

MODULE LDALIB

    CONTAINS

!*******************************************************************
!
! In the following a couple of "standard" exchange
! and correlation functionals are implemented
! these usually use a.u. and Hartree-units 
! they calculate the energy density per electrons and the
! corresponding potential
!
! E_xc = \int eps_xc(rho(r)) rho(r) d^3 r
! v_xc = d (eps_xc(rho) rho)   / d rho 
!
!*******************************************************************

!************************ SUBROUTINE XCTABLE_CHECK *********************
!
!> procedure that checks whether the exchange correlation table
!> corresponds to the currently used exchange correlation type XC_DATA
!
!***********************************************************************

    SUBROUTINE XCTABLE_CHECK(XC)
      USE tutor, ONLY: vtutor
      USE setexm_struct_def, ONLY: xc_info, XC_DATA_TABLE
      IMPLICIT NONE
      TYPE (xc_info) XC
      INTEGER I

      IF (XC%NXC/=XC_DATA_TABLE%NXC) THEN
         CALL vtutor%bug("internal ERROR in XCTABLE_CHECK: the exchange-correlation table is not properly set", __FILE__, __LINE__)
      ENDIF

      DO I=1,XC%NXC
         IF (XC%ID(I)/=XC_DATA_TABLE%ID(I)) THEN
            CALL vtutor%bug("internal ERROR in XCTABLE_CHECK: the exchange-correlation table is not properly set", __FILE__, __LINE__)
         ENDIF
      ENDDO

    END SUBROUTINE XCTABLE_CHECK

!************************ SUBROUTINE FEXCF *****************************
!> RCS:  $Id: potex1.F,v 1.2 2000/11/15 08:23:50 kresse Exp $
!>
!> this subroutine calculates the LDA contribution to the 
!> exchange correlation potential from the charge density and
!> magnetisation density supplied in CHDENR and CHMAGR (real space)
!> The energy corrections due to overcounting the exchange correlation 
!> energy on summing the electronic eigenvalues and the forces on the unit 
!> cell due to the exchange correlation energy are also calculated.
!>
!> To make the calculations efficient the required functions are
!> interpolated from tables (see setex.F)
!> this makes the routine extremely efficient without sacrificing
!> precision
!***********************************************************************

      SUBROUTINE FEXCF(GRID,XC,OMEGA, &
          CHDENR,CHMAGR,DENCOR,CVXC1,CVXC2,CVZERO,EXC,XCENC,XCSIF,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE base, only : CAST_COMPLEX_TO_RGRID
      USE mpimy
      USE mgrid
      USE setexm_struct_def, only : xc_info,exctab
      IMPLICIT NONE

      TYPE (grid_3d)  GRID
      TYPE (xc_info)  XC

      INTEGER I,IMAX,IMIN,IONODE,IZ,N,NODE_ME
      REAL(q) ALPH0,ARG,DALPH,DECA,DECA1,DECF,DECF0,DECF1,DENS,DEXF,DEXF0,DEXF1, &
     &        DMAG,DVAL,DVCA,DVCF,DVCF0,DVXC,DVXF,DVXF0,DZETC,DZETX, &
     &        EXC,EXCP,EXCP0,EXCSC1,EXCSC2,EXCT,FZETC,FZETX,OMEGA,REM, &
     &        RHO,RHO13,RHO43,RHOD,RHOM,RHOP,VXC0,VXC1,VXC2,VXCP,VXCP0,VXCS, &
     &        XCENC,XCF,XCFO,ZARG,ZETA,ZETA3,ZETA4,ZREM
      COMPLEX(q) CVZERO

      COMPLEX(q), TARGET :: CHDENR(GRID%MPLWV),CHMAGR(GRID%MPLWV)
      RGRID  DENCOR(GRID%RL%NP)
      COMPLEX(q), TARGET :: CVXC1(GRID%MPLWV),CVXC2(GRID%MPLWV)

      LOGICAL      LADD
      REAL(q)      XCSIF(3,3)

      RGRID,pointer :: DCHDENR(:) => null()
      RGRID,pointer :: DCHMAGR(:) => null()
      RGRID,pointer :: DCVXC1(:) => null()
      RGRID,pointer :: DCVXC2(:) => null()

      CALL CAST_COMPLEX_TO_RGRID(CHDENR,DCHDENR)
      CALL CAST_COMPLEX_TO_RGRID(CHMAGR,DCHMAGR)
      CALL CAST_COMPLEX_TO_RGRID(CVXC1,DCVXC1)
      CALL CAST_COMPLEX_TO_RGRID(CVXC2,DCVXC2)

!=======================================================================
!  calculate exchange-correlation on the grid using the tables
!=======================================================================

      CALL XCTABLE_CHECK(XC)

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      IMAX=0

      CVZERO=0._q
      EXC=0._q
      XCF=0._q
      XCFO=0._q

      EXCSC2=(EXCTAB%NEXCHF(2)-EXCTAB%NEXCHF(1))/(EXCTAB%RHOEXC(2)-EXCTAB%RHOEXC(1))
      EXCSC1= EXCTAB%NEXCHF(1)/MAX(EXCTAB%RHOEXC(1),1E-10_q)

      IMIN=EXCTAB%NEXCHF(2)
!DIR$ IVDEP
!cdir nodep
!$DIR FORCE_VECTOR
!OCL NOVREC

GPACC !$ACC PARALLEL LOOP PRESENT(DCHDENR, DENCOR, DCHMAGR, DCVXC1, DCVXC2, EXCTAB, EXCTAB%RHOEXC, EXCTAB%NEXCHF, EXCTAB%EXCTAB, GRID) &
GPACC !$ACC& REDUCTION(+:CVZERO,EXC,XCF,XCFO) REDUCTION(max:IMAX) &
GPACC !$ACC& PRIVATE(DVAL, DENS, DMAG, ZETA, ZETA3, ZETA4, RHO, RHO13, RHO43, RHOD, RHOP, RHOM, ARG, ZARG, I, IZ, REM, ZREM, &
GPACC !$ACC& EXCP0, VXCP0, EXCP, VXCP, DEXF0, DVXF0, DEXF1, DEXF, DVXF, DECF0, DVCF0, DECF1, DECF, DVCF, ALPH0, DALPH, DECA1, &
GPACC !$ACC& DECA, DVCA, FZETX, DZETX, FZETC, DZETC, EXCT, VXC0, DVXC, VXC1, VXC2, VXCS) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DCHDENR, DENCOR, DCHMAGR, DCVXC1, DCVXC2, EXCTAB, &
GPOMP !$OMP EXCTAB%RHOEXC, EXCTAB%NEXCHF, EXCTAB%EXCTAB, GRID) REDUCTION(+:CVZERO,EXC,XCF,XCFO) REDUCTION(max:IMAX) &
GPOMP !$OMP PRIVATE(DVAL, DENS, DMAG, ZETA, ZETA3, ZETA4, RHO, RHO13, RHO43, RHOD, RHOP, RHOM, ARG, ZARG, I, IZ, REM, ZREM, &
GPOMP !$OMP EXCP0, VXCP0, EXCP, VXCP, DEXF0, DVXF0, DEXF1, DEXF, DVXF, DECF0, DVCF0, DECF1, DECF, DVCF, ALPH0, DALPH, DECA1, &
GPOMP !$OMP DECA, DVCA, FZETX, DZETX, FZETC, DZETC, EXCT, VXC0, DVXC, VXC1, VXC2, VXCS) IF(OFFLOAD_ON)

!NEC$ nomove
      DO N=1,GRID%RL%NP

        DVAL= REAL( DCHDENR(N) ,KIND=q)
        DENS=DVAL+DENCOR(N)
        DENS=MAX(DENS,1E-8_q)

        DMAG= REAL( DCHMAGR(N) ,KIND=q)

        ZETA=DMAG/DENS
        IF (ZETA>= 1._q) ZETA= 0.9999999999999_q
        IF (ZETA<=-1._q) ZETA=-0.9999999999999_q
        ZETA3=ZETA*ZETA*ZETA
        ZETA4=ZETA3*ZETA

        RHO  =DENS/OMEGA
        RHO13=EXP(LOG(RHO)/3._q)
        RHO43=RHO*RHO13
        RHOD =4._q*RHO13/3._q
        RHOP =0.5_q*(DVAL+DMAG)
        RHOM =0.5_q*(DVAL-DMAG)
        ARG  =(RHO*EXCSC1)+1
        IF (ARG>=EXCTAB%NEXCHF(1)+1) ARG=(RHO-EXCTAB%RHOEXC(1))*EXCSC2+EXCTAB%NEXCHF(1)+1
        ZARG =ABS(ZETA)*(EXCTAB%NEXCHF(2)-1)+1

        I   =INT(ARG)
        IMAX=MAX(I,IMAX)
        I   =MIN(I,EXCTAB%NEXCHF(2))
        IZ  =INT(ZARG)
        REM =RHO-EXCTAB%EXCTAB(I,1,1)
        ZREM=ABS(ZETA)-EXCTAB%EXCTAB(IZ,1,5)

        EXCP0=EXCTAB%EXCTAB(I,2,1)+REM*(EXCTAB%EXCTAB(I,3,1)+ &
     &                      REM*(EXCTAB%EXCTAB(I,4,1)+REM*EXCTAB%EXCTAB(I,5,1)))
        VXCP0=EXCTAB%EXCTAB(I,3,1)+REM*(EXCTAB%EXCTAB(I,4,1)*2+REM*EXCTAB%EXCTAB(I,5,1)*3)
        EXCP =EXCP0*RHO43
        VXCP =VXCP0*RHO43+EXCP0*RHOD

        DEXF0=EXCTAB%EXCTAB(I,2,2)+REM*(EXCTAB%EXCTAB(I,3,2)+ &
     &                      REM*(EXCTAB%EXCTAB(I,4,2)+REM*EXCTAB%EXCTAB(I,5,2)))
        DVXF0=EXCTAB%EXCTAB(I,3,2)+REM*(EXCTAB%EXCTAB(I,4,2)*2+REM*EXCTAB%EXCTAB(I,5,2)*3)
        DEXF1=DEXF0*RHO13
        DEXF =DEXF0*RHO43
        DVXF =DVXF0*RHO43+DEXF0*RHOD

        DECF0=EXCTAB%EXCTAB(I,2,3)+REM*(EXCTAB%EXCTAB(I,3,3)+ &
     &                      REM*(EXCTAB%EXCTAB(I,4,3)+REM*EXCTAB%EXCTAB(I,5,3)))
        DVCF0=EXCTAB%EXCTAB(I,3,3)+REM*(EXCTAB%EXCTAB(I,4,3)*2+REM*EXCTAB%EXCTAB(I,5,3)*3)
        DECF1=DECF0*RHO13
        DECF =DECF0*RHO43
        DVCF =DVCF0*RHO43+DECF0*RHOD

        ALPH0=EXCTAB%EXCTAB(I,2,4)+REM*(EXCTAB%EXCTAB(I,3,4)+ &
     &                      REM*(EXCTAB%EXCTAB(I,4,4)+REM*EXCTAB%EXCTAB(I,5,4)))
        DALPH=EXCTAB%EXCTAB(I,3,4)+REM*(EXCTAB%EXCTAB(I,4,4)*2+REM*EXCTAB%EXCTAB(I,5,4)*3)
        DECA1=ALPH0*RHO13
        DECA =ALPH0*RHO43
        DVCA =DALPH*RHO43+ALPH0*RHOD

        FZETX=EXCTAB%EXCTAB(IZ,2,5)+ZREM*(EXCTAB%EXCTAB(IZ,3,5)+ &
     &                      ZREM*(EXCTAB%EXCTAB(IZ,4,5)+ZREM*EXCTAB%EXCTAB(IZ,5,5)))
        DZETX=EXCTAB%EXCTAB(IZ,3,5)+ZREM*(EXCTAB%EXCTAB(IZ,4,5)*2+ &
     &                     ZREM*EXCTAB%EXCTAB(IZ,5,5)*3)
        DZETX=SIGN(1._q,ZETA)*ZETA*ZETA*DZETX+2*ZETA*FZETX
        FZETX=ZETA*ZETA*FZETX

        FZETC=EXCTAB%EXCTAB(IZ,2,6)+ZREM*(EXCTAB%EXCTAB(IZ,3,6)+ &
     &                      ZREM*(EXCTAB%EXCTAB(IZ,4,6)+ZREM*EXCTAB%EXCTAB(IZ,5,6)))
        DZETC=EXCTAB%EXCTAB(IZ,3,6)+ZREM*(EXCTAB%EXCTAB(IZ,4,6)*2+ &
     &                     ZREM*EXCTAB%EXCTAB(IZ,5,6)*3)
        DZETC=SIGN(1._q,ZETA)*ZETA*ZETA*DZETC+2*ZETA*FZETC
        FZETC=ZETA*ZETA*FZETC

        EXCT =EXCP+DEXF*FZETX+(DECF+DECA*ZETA4)*FZETC
        VXC0 =VXCP+DVXF*FZETX+(DVCF+DVCA*ZETA4)*FZETC
        DVXC =DEXF1*DZETX+(DECF1+DECA1*ZETA4)*DZETC+4._q*DECA1*ZETA3*FZETC
        VXC1 =VXC0-DVXC*(ZETA-1)
        VXC2 =VXC0-DVXC*(ZETA+1)
        VXCS =0.5_q*(VXC1+VXC2)

        CVZERO=CVZERO+VXCS
        EXC   =EXC   +EXCT
        XCF   =XCF   -VXC1*RHOP-VXC2*RHOM
        XCFO  =XCFO  -VXCS*DENCOR(N)

        IF (LADD) THEN
           DCVXC1(N)= DCVXC1(N)+VXC1
           DCVXC2(N)= DCVXC2(N)+VXC2
        ELSE
           DCVXC1(N)= VXC1
           DCVXC2(N)= VXC2
        ENDIF

      ENDDO
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      CALLMPI( M_max_i(GRID%COMM, IMAX, 1))
!=======================================================================
!  check if the size of the table is sufficent
!=======================================================================
      IF (IMAX>EXCTAB%NEXCHF(2)*.96_q) THEN
        CALL vtutor%error("ERROR FEXCF: supplied exchange-correlation table \n is too small, maximal &
           &index : " // str(IMAX))
      ENDIF
!=======================================================================
! calculate the g=0 component of the exchange-correlation potential
! CVZERO which is used in the electron dynamics subroutine
!=======================================================================
      CVZERO=CVZERO/GRID%NPLWV

      EXC=EXC*OMEGA
      XCF =EXC+XCF
      EXC =EXC /GRID%NPLWV
      XCF =XCF /GRID%NPLWV
      XCFO=XCFO/GRID%NPLWV

      CALLMPI( M_sum_3(GRID%COMM,EXC,XCF,XCFO))
      CALLMPI( M_sum_z(GRID%COMM,CVZERO,1))
!=======================================================================
! calculate the force on the unit cell due to the change in the exchange
! correlation energy on changing the size of the cell
! the uniform part of the partial core correction (from
! the 1/OMEGA dependency of the core charge is treated here)
!=======================================================================
#ifdef debug
        do_io  WRITE(*,'(A,F24.14)') '<rho*exc> =',EXC
        do_io  WRITE(*,'(A,F24.14)') '<rho*vxc> =',EXC-XCF
        do_io  WRITE(*,'(A,F24.14)') '    xcenc =',XCF
#endif
      XCSIF=0._q

      XCSIF(1,1)=-XCF-XCFO
      XCSIF(2,2)=-XCF-XCFO
      XCSIF(3,3)=-XCF-XCFO
!=======================================================================
! double counting correction to the total energy
!=======================================================================
      XCENC=XCF

      RETURN
      END SUBROUTINE FEXCF

!************************ SUBROUTINE FEXCP *****************************
!> RCS:  $Id: potex2.F,v 1.3 2001/01/31 11:51:53 kresse Exp $
!>
!> this subroutine calculates the LDA contribution to the 
!> exchange correlation potential from the charge density supplied
!> int CHDENR (real space)
!> The energy corrections due to overcounting the exchange correlation 
!> energy on summing the electronic eigenvalues and the forces on the unit 
!> cell due to the exchange correlation energy are also calculated.
!>
!> To make the calculations efficient the required functions are
!> interpolated from tables (see setex.F)
!> this makes the routine extremely efficient without sacrificing
!> precision
!***********************************************************************

      SUBROUTINE FEXCP(GRID,XC,OMEGA, &
          CHDENR,DENCOR,CVXC,DVXC,CVZERO,EXC,XCENC,XCSIF,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE base, only : CAST_COMPLEX_TO_RGRID
      USE mpimy
      USE mgrid
      USE setexm_struct_def, only : xc_info,exctab
      IMPLICIT NONE

      TYPE (grid_3d)  GRID
      TYPE (xc_info)  XC

      INTEGER I,IMAX,IMIN,IONODE,N,NODE_ME
      REAL(q) ARG,DENS,EXC,EXC0,EXCD,EXCDD,EXCSC1,EXCSC2,EXE,OMEGA,REM, &
     &        RHO,RHO0,RHO13,RHOD,RHODD,VXC,XCENC,XCF,XCFO
      COMPLEX(q) CVZERO

      COMPLEX(q), TARGET :: CHDENR(GRID%MPLWV)
      RGRID      DENCOR(GRID%RL%NP)
      COMPLEX(q), TARGET :: DVXC(GRID%MPLWV)
      COMPLEX(q), TARGET :: CVXC(GRID%MPLWV)

      REAL(q)  XCSIF(3,3)
      LOGICAL :: LADD

      RGRID,pointer :: DCHDENR(:) => null()
      RGRID,pointer :: DCDVXC(:) => null()
      RGRID,pointer :: DCCVXC(:) => null()

      CALL CAST_COMPLEX_TO_RGRID(CHDENR,DCHDENR)
      CALL CAST_COMPLEX_TO_RGRID(DVXC,DCDVXC)
      CALL CAST_COMPLEX_TO_RGRID(CVXC,DCCVXC)

      CALL XCTABLE_CHECK(XC)
!=======================================================================
!  calculate exchange-correlation on the grid using the tables
!=======================================================================
      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      IMAX=0

      CVZERO=0._q
      EXC=0._q
      XCF=0._q
      XCFO=0._q

      EXCSC2=(EXCTAB%NEXCHF(2)-EXCTAB%NEXCHF(1))/(EXCTAB%RHOEXC(2)-EXCTAB%RHOEXC(1))
      EXCSC1=EXCTAB%NEXCHF(1)/MAX(EXCTAB%RHOEXC(1),1E-10_q)

      IMIN=EXCTAB%NEXCHF(2)
__CPU !DIR$ IVDEP
__CPU !cdir nodep
__CPU !$DIR FORCE_VECTOR
__CPU !OCL NOVREC

GPACC !$ACC PARALLEL LOOP PRESENT(OMEGA, DCDVXC, DCCVXC, DENCOR, DCHDENR, EXCTAB, EXCTAB%RHOEXC, EXCTAB%NEXCHF, EXCTAB%EXCTAB, GRID) &
GPACC !$ACC& REDUCTION(+:CVZERO,EXC,XCF,XCFO) REDUCTION(max:IMAX) &
GPACC !$ACC& PRIVATE(DENS,RHO,ARG,I,REM,RHO13,RHO0,RHOD,RHODD,EXC0,EXCD,EXCDD,EXE,VXC) &
GPACC !$ACC& __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:OMEGA, DCDVXC, DCCVXC, DENCOR, DCHDENR, EXCTAB, &
GPOMP !$OMP EXCTAB%RHOEXC, EXCTAB%NEXCHF, EXCTAB%EXCTAB, GRID) REDUCTION(+:CVZERO,EXC,XCF,XCFO) REDUCTION(max:IMAX) &
GPOMP !$OMP PRIVATE(DENS,RHO,ARG,I,REM,RHO13,RHO0,RHOD,RHODD,EXC0,EXCD,EXCDD,EXE,VXC) IF(OFFLOAD_ON)
      DO N=1,GRID%RL%NP
        DENS= REAL( DCHDENR(N) ,KIND=q) +DENCOR(N)
        DENS=MAX(DENS,1E-8_q)

        RHO =ABS(DENS)/OMEGA
        ARG=(RHO*EXCSC1)+1
        IF (ARG>=EXCTAB%NEXCHF(1)+1) ARG=(RHO-EXCTAB%RHOEXC(1))*EXCSC2+EXCTAB%NEXCHF(1)+1

        I=INT(ARG)
        IMAX=MAX(I,IMAX)

        I=MIN(I,EXCTAB%NEXCHF(2))
        
        REM=RHO-EXCTAB%EXCTAB(I,1,1)
        RHO13=EXP(LOG(RHO)/3._q)

        RHO0   = RHO13*RHO
        RHOD   = 4._q/3._q*RHO13
        RHODD  = 4._q/9._q*RHO13/MAX(RHO,1E-3_q)

        EXC0   = EXCTAB%EXCTAB(I,2,1)+REM*(EXCTAB%EXCTAB(I,3,1)+ &
     &                       REM*(EXCTAB%EXCTAB(I,4,1)+REM*EXCTAB%EXCTAB(I,5,1)))
        EXCD   = EXCTAB%EXCTAB(I,3,1)+REM*(EXCTAB%EXCTAB(I,4,1)*2+REM*EXCTAB%EXCTAB(I,5,1)*3)
        EXCDD  = EXCTAB%EXCTAB(I,4,1)*2+REM*EXCTAB%EXCTAB(I,5,1)*6

        EXE    = EXC0*RHO0
        VXC    = EXC0*RHOD+EXCD*RHO0

        CVZERO=CVZERO+VXC
        EXC   =EXC   +EXE
        XCF   =XCF   -VXC*REAL( DCHDENR(N) ,KIND=q)
        XCFO  =XCFO  -VXC*DENCOR(N)

        IF (LADD) THEN
           DCCVXC(N)= DCCVXC(N)+VXC
        ELSE
           DCCVXC(N)= VXC
        ENDIF
        DCDVXC(N)= EXC0*RHODD+2*EXCD*RHOD+EXCDD*RHO0
      ENDDO
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!=======================================================================
!  check if the size of the table is sufficent
!=======================================================================
      CALLMPI( M_max_i(GRID%COMM, IMAX, 1))
      IF (IMAX>EXCTAB%NEXCHF(2)*.96_q) THEN
        CALL vtutor%error("ERROR FEXCP: supplied Exchange-correletion table \n is too small, maximal &
           &index : " // str(IMAX))
      ENDIF
!=======================================================================
! calculate the g=0 component of the exchange-correlation potential
! CVZERO which is used in the electron dynamics subroutine
!=======================================================================
      CVZERO=CVZERO/GRID%NPLWV

      EXC=EXC*OMEGA
      XCF =EXC+XCF
      EXC =EXC /GRID%NPLWV
      XCF =XCF /GRID%NPLWV
      XCFO=XCFO/GRID%NPLWV

      CALLMPI( M_sum_3(GRID%COMM,EXC,XCF,XCFO))
      CALLMPI( M_sum_z(GRID%COMM,CVZERO,1))
!=======================================================================
! calculate the force on the unit cell due to the change in the exchange
! correlation energy on changing the size of the cell
! the uniform part of the partial core correction (from
! the 1/OMEGA dependency of the core charge is treated here)
!=======================================================================
!#define debug
#ifdef debug
        do_io  WRITE(*,'(A,F24.14)') '<rho*exc> =',EXC
        do_io  WRITE(*,'(A,F24.14)') '<rho*vxc> =',EXC-XCF
        do_io  WRITE(*,'(A,F24.14)') '    xcenc =',XCF
#endif
      XCSIF=0._q

      XCSIF(1,1)=-XCF-XCFO
      XCSIF(2,2)=-XCF-XCFO
      XCSIF(3,3)=-XCF-XCFO
!=======================================================================
! correction to the total energy
!=======================================================================
      XCENC=XCF

      RETURN
      END SUBROUTINE FEXCP

    FUNCTION ECCA(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Ceperley-Alder correlation energy as parametrised by Perdew/Zunger
! (see Phys.Rev. B23,5048 [1981], Appendix).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      REAL(q),SAVE ::  A(2)=(/0.0622_q,0.0311_q/)
      REAL(q),SAVE ::  B(2)=(/-0.0960_q,-0.0538_q/)
!     REAL(q),SAVE ::  C(2)=(/0.0040_q,0.0014_q/)
      REAL(q),SAVE ::  C(2)=(/0.004038664055501747_q,0.001395274602717559_q/)
!     REAL(q),SAVE ::  D(2)=(/-0.0232_q,-0.0096_q/)
      REAL(q),SAVE ::  D(2)=(/-0.023264632546756681_q,-0.009602765503781227_q/)
      REAL(q),SAVE ::  G(2)=(/-0.2846_q,-0.1686_q/)
      REAL(q),SAVE ::  B1(2)=(/1.0529_q,1.3981_q/)
!     REAL(q),SAVE ::  B2(2)=(/0.3334_q,0.2611_q/)
      REAL(q),SAVE ::  B2(2)=(/0.333390000000000000_q, 0.261090000000000000_q/)
      IF (RS<=1.0_q) THEN
         RSL=LOG(RS)
         ECCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
         RSQ=SQRT(RS)
         ECCA=G(IFLG)/(1.0_q+B1(IFLG)*RSQ+B2(IFLG)*RS)
      END IF
      RETURN
    END FUNCTION ECCA

    FUNCTION VCCA(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Ceperley-Alder correlation potential as parametrised by Perdew/Zunger
! (see Phys.Rev. B23,5048 [1981], Appendix).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),B(2),C(2),D(2),BT1(2),BT2(2)
      PARAMETER(X76=7.0_q/6.0_q, X43=4.0_q/3.0_q, &
     &  AP=0.03110_q*2.0_q, BP=-0.0480_q*2.0_q, CP=0.0020_q*2.0_q, DP=-0.0116_q*2.0_q, &
     &  AF=0.01555_q*2.0_q, BF=-0.0269_q*2.0_q, CF=0.0007_q*2.0_q, DF=-0.0048_q*2.0_q, &
     &  BP1=BP-AP/3.0_q, CP1=2.0_q*CP/3.0_q, DP1=(2.0_q*DP-CP)/3.0_q, &
     &  BF1=BF-AF/3.0_q, CF1=2.0_q*CF/3.0_q, DF1=(2.0_q*DF-CF)/3.0_q)
      SAVE A,B,C,D,BT1,BT2
      DATA A/AP,AF/,B/BP1,BF1/,C/CP1,CF1/,D/DP1,DF1/, &
     &     BT1/1.0529_q,1.3981_q/,BT2/0.3334_q,0.2611_q/
!KRESSE/FURTH---get a continous energy functional
      c(1)  = 0.004038664055501747_q * 2._q/3._q
      d(1)  =-0.023264632546756681_q * 2._q/3._q  -  0.004038664055501747_q / 3._q
      bt2(1)= 0.333390000000000000_q
      c(2)  = 0.001395274602717559_q * 2._q/3._q
      d(2)  =-0.009602765503781227_q * 2._q/3._q  -  0.001395274602717559_q / 3._q
      bt2(2)= 0.261090000000000000_q
!KRESSE/FURTH
      IF (RS<=1.0_q) THEN
         RSL=LOG(RS)
         VCCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
         RSQ=SQRT(RS)
         VCCA=ECCA(RS,IFLG)*(1.0_q+X76*BT1(IFLG)*RSQ+X43*BT2(IFLG)*RS) &
     &                     /(1.0_q+    BT1(IFLG)*RSQ+    BT2(IFLG)*RS)
      END IF
      RETURN
    END FUNCTION VCCA


    FUNCTION ECVO(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The Ceperley-Alder correlation energy as given by the Pade approximation
! technique of Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with
! the parameters given in table 5, page 1207).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.104980_q,BP=3.72744_q,CP=12.9352_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.325000_q,BF=7.06042_q,CF=18.0578_q)
      PARAMETER(QP=6.1519908_q,QF=4.7309269_q)
      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      ECVO=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ECVO=-1._q*ECVO*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ECVO=ECVO*A(IFLG)
      RETURN
    END FUNCTION ECVO

    FUNCTION VCVO(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ECVO(RS,IFLG)-RS/3.*d(ECVO(RS,IFLG))/d(RS) with function
! ECVO(RS,IFLG) given above (Ceperley-Alder potential derived from the
! approximation for ECVO of Vosko et al. discussed above).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.104980_q,BP=3.72744_q,CP=12.9352_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.325000_q,BF=7.06042_q,CF=18.0578_q)
      PARAMETER(QP=6.1519908_q,QF=4.7309269_q)
      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      VCVO=-4._q*B*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      VCVO=VCVO-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      VCVO=ECVO(RS,IFLG)-VCVO*A(IFLG)*X/6._q
      RETURN
    END FUNCTION VCVO
!jP---------------------------------------------------------------------------------
    FUNCTION ECVOIII(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The RPA correlation energy as given by the Pade approximation
! technique of Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with
! the parameters given in table 5, page 1207).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.409286_q,BP=13.0720_q,CP=42.7198_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.743294_q,BF=20.1231_q,CF=101.578_q)
      PARAMETER(QP=0.0448999_q,QF=1.1716853_q)

      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      ECVOIII=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ECVOIII=-1._q*ECVOIII*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ECVOIII=ECVOIII*A(IFLG)
      RETURN
    END FUNCTION ECVOIII

    FUNCTION VCVOIII(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ECVOIII(RS,IFLG)-RS/3.*d(ECVOIII(RS,IFLG))/d(RS) with function
! ECVOIII(RS,IFLG) given above (RPA potential derived from the
! approximation for ECVOIII of Vosko et al. discussed above).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.409286_q,BP=13.0720_q,CP=42.7198_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.743294_q,BF=20.1231_q,CF=101.578_q)
      PARAMETER(QP=0.0448999_q,QF=1.1716853_q)

      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      VCVOIII=-4._q*B*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      VCVOIII=VCVOIII-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      VCVOIII=ECVOIII(RS,IFLG)-VCVOIII*A(IFLG)*X/6._q
      RETURN
    END FUNCTION VCVOIII
!jP---------------------------------------------------------------------------------

    FUNCTION ECGL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Gunnarson-Lundqvist correlation energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.0666_q,0.0406_q/,R/11.4_q,15.9_q/
      X=RS/R(IFLG)
      ECGL=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
    END FUNCTION ECGL

    FUNCTION VCGL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Gunnarson-Lundqvist correlation potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.0666_q,0.0406_q/,R/11.4_q,15.9_q/
      VCGL=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
    END FUNCTION VCGL

    FUNCTION ECHL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Hedin-Lundqvist correlation energy (J.Phys. C4,2064 [1971]):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.045_q,0.0225_q/,R/21.0_q,52.917_q/
      X=RS/R(IFLG)
      ECHL=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
    END FUNCTION ECHL

    FUNCTION VCHL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Hedin-Lundqvist correlation potential (J.Phys. C4,2064 [1971]):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.045_q,0.0225_q/,R/21.0_q,52.917_q/
      VCHL=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
    END FUNCTION VCHL

    FUNCTION ECBH(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Barth-Hedin correlation energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.0504_q,0.0254_q/,R/30._q,75._q/
      X=RS/R(IFLG)
      ECBH=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
    END FUNCTION ECBH

    FUNCTION VCBH(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Barth-Hedin correlation potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.0504_q,0.0254_q/,R/30._q,75._q/
      VCBH=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
    END FUNCTION VCBH

    FUNCTION ECWI(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Wigner correlation energy (hopefully correct?):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results --> warning: equals paramagnetic result
      DIMENSION CX(2),C(2),R(2)
      SAVE CX,C,R
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      DATA C/7.8_q,7.8_q/,R/0.88_q,0.88_q/
      X=CX(IFLG)/R(IFLG)*(1._q+C(IFLG)/RS)
      ECWI=-CX(IFLG)/X/RS
      RETURN
    END FUNCTION ECWI

    FUNCTION VCWI(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Wigner correlation potential (hopefully correct?):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results --> warning: equals paramagnetic result
      DIMENSION CX(2),C(2),R(2)
      SAVE CX,C,R
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      DATA C/7.8_q,7.8_q/,R/0.88_q,0.88_q/
      X1=C(IFLG)/RS
      X2=1._q+X1
      X3=CX(IFLG)/R(IFLG)*X2
      B=1._q+1._q/X3
      F=1._q-X1/X2/(1._q+X3)
      F=1._q+F/3._q
      E=-CX(IFLG)*B/RS
      VCWI=(4._q/3._q-F*B)*CX(IFLG)/RS
      RETURN
    END FUNCTION VCWI

    FUNCTION EX(RS,IFLG,TREL)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EX)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
!     TREL  :  Relativistic correction or not (for details see
!              J.Phys. C12,2977(1979) )
      LOGICAL TREL
      DIMENSION CX(2)
      SAVE CX,CBETA
      DATA CX/0.9163305865663_q,1.1545041946774_q/,CBETA/0.0140_q/
      EX=-CX(IFLG)/RS
      IF (TREL) THEN
         B=CBETA/RS
         F=LOG(B+(SQRT(1+B*B)))/(B*B)
         F=(SQRT(1+B*B)/B)-F
!jF: the expression given above becomes numerically extremely instable for
!    very small values of B (small difference of two large numbers divided
!    by small number = noise) therefore use following for reasons of safety:
         IF (B.LT.1.E-5_q) F=B*(2._q/3._q-0.7_q*B*B)
         EX=(1._q-1.5_q*F*F)*EX
      END IF
      RETURN
    END FUNCTION EX

    FUNCTION VX(RS,IFLG,TREL)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VX)
! Exchange potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
!     TREL  :  Relativistic correction or not (for details see
!              J.Phys. C12,2977(1979) )
      LOGICAL TREL
      DIMENSION CX(2)
      SAVE CX,CBETA
      DATA CX/1.2217741154217_q,1.5393389262365_q/,CBETA/0.0140_q/

      VX=-CX(IFLG)/RS
      IF (TREL) THEN
! Warning error in the paper of Bachelet et al. !!
         B=CBETA/RS
         F=LOG(B+(SQRT(1+B*B)))/B/SQRT(1+(B*B))
!        F=LOG(B+(SQRT(1+B*B)))/B/(1+(B*B))
!jF: potentially the expression given above becomes numerically instable for
!    very small values of B, therefore use following for reasons of safety:
         IF (B.LT.1.E-5_q) F=1._q-B*B*(2._q/3._q-B*B*31._q/30._q)
         VX=(-0.5_q+1.5_q*F)*VX
      END IF
      RETURN
    END FUNCTION VX
      
    FUNCTION FZ0(ZETA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Interpolation function between paramagnetic and ferromagnetic results
! for exchange energy and exchange potential. The parameter ZETA is:
! ZETA=(RHO[upspin] - RHO[downspin]) / (RHO[upspin] + RHO[downspin]).
      PARAMETER(C43=4._q/3._q,FAC=1.92366105093153632_q)
      Z=ABS(ZETA)
      IF (Z==0._q) THEN
         FZ0=0._q
      ELSE IF (Z>=1._q) THEN
         FZ0=1._q
      ELSE
         FZ0=(((1._q+Z)**C43)+((1._q-Z)**C43)-2._q)*FAC
      END IF
      RETURN
    END FUNCTION FZ0

    FUNCTION FZ1(ZETA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The derivative dFZ0(ZETA)/d(ZETA), FZ0(ZETA) given above.
      PARAMETER(C13=1._q/3._q,FAC=2.56488140124204843_q)
      Z=ABS(ZETA)
      IF (Z==0._q) THEN
         FZ1=0._q
      ELSE IF (Z>=1._q) THEN
         FZ1=(2._q**C13)*FAC
      ELSE
         FZ1=(((1._q+Z)**C13)-((1._q-Z))**C13)*FAC
      END IF
      RETURN
    END FUNCTION FZ1

    FUNCTION ALPHA0(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The spin stiffness [d(EXC(RS,ZETA))**2/d**2(ZETA) |ZETA=0] as given by
! Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with the parameters
! given on page 1209 [fitting to low-density values using eq.{4.7.}]).
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.0047584_q,B=1.13107_q,C=13.0045_q,QQ=7.12311_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.019751631321681_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA0=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0=-1._q*ALPHA0*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0=ALPHA0*A
      RETURN
    END FUNCTION ALPHA0

    FUNCTION ALPHA1(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ALPHA0(RS)-RS/3.*dALPHA0(RS)/dRS, ALPHA0(RS) given above.
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-.0047584_q,B=1.13107_q,C=13.0045_q,QQ=7.12311_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.0032919385536135_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA1=-4._q*B/QQ*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      ALPHA1=ALPHA1-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      ALPHA1=ALPHA0(RS)-ALPHA1*A*X
      RETURN
    END FUNCTION ALPHA1
!----------------------------------------------------------------------
    FUNCTION ALPHA0_III(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The spin stiffness [d(EXC(RS,ZETA))**2/d**2(ZETA) |ZETA=0] as given by
! Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with the parameters
! given on page 1209 [fitting to low-density values using eq.{4.7.}]).
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.228344_q,B=1.06835_q,C=11.4813_q,QQ=6.69207_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.019751631321681_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA0_III=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0_III=-1._q*ALPHA0_III*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0_III=ALPHA0_III*A
      RETURN
    END FUNCTION ALPHA0_III

    FUNCTION ALPHA1_III(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ALPHA0(RS)-RS/3.*dALPHA0(RS)/dRS, ALPHA0(RS) given above.
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.228344_q,B=1.06835_q,C=11.4813_q,QQ=6.69207_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.0032919385536135_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA1_III=-4._q*B/QQ*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      ALPHA1_III=ALPHA1_III-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      ALPHA1_III=ALPHA0_III(RS)-ALPHA1_III*A*X
      RETURN
    END FUNCTION ALPHA1_III
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    SUBROUTINE CORPBE_LDA(RS,ZET,EC,VCUP,VCDN)
!----------------------------------------------------------------------
!  LDA part of the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORPBE_LDA)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor_xc(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      EC = EU*(1._q-F*Z4)+EP*F*Z4-ALFM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1._q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZET = 4._q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1._q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3._q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA
    
    
!----------------------------------------------------------------------
!> Judith Harl
!> functionals for range-separated ACFDT (LDA - short range RPA):
!> a bit akward at the moment since the range separation parameter
!> is hard coded for now. Hopefully this will change ...
!----------------------------------------------------------------------

    SUBROUTINE CORPBE_LDA_RPA(RS,ZET,EC,VCUP,VCDN)
! RPA correlation energy for jellium 
! Perdew, Wang, prB 45 13244 (1992)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      CALL gcor_xc_rpa(0.015545_q,0.035374_q,6.4869_q,1.3083_q,0.15180_q, &
     &    0.082349_q,rtRS,EP,EPRS)
      CALL gcor_xc(0.016887_q,0.028829_q,10.357_q,3.6231_q,0.47990_q, &
     &    0.12279_q,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      EC = EU*(1._q-F*Z4)+EP*F*Z4-ALFM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1._q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZET = 4._q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1._q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3._q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_RPA


!  jH-new This routine is used to evaluate the RPA+ correction term, 
!  jH-new which is E^{LDA,QMC}_{c}-E^{LDA,RPA}_{c}

    SUBROUTINE CORPBE_LDA_RPA_PLUS(RS,ZET,EC,VCUP,VCDN)
! RPA correlation energy for jellium 
! Perdew, Wang, prB 45 13244 (1992)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EURPA,EURSRPA)
      CALL gcor_xc_rpa(0.015545_q,0.035374_q,6.4869_q,1.3083_q,0.15180_q, &
     &    0.082349_q,rtRS,EPRPA,EPRSRPA)
      CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EUQM,EURSQM)
      CALL gcor_xc(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EPQM,EPRSQM)
      EU=EUQM-EURPA
      EP=EPQM-EPRPA
      EURS=EURSQM-EURSRPA
      EPRS=EPRSQM-EPRSRPA
! jH rpa spin stiffness
      CALL gcor_xc(0.016887_q,0.028829_q,10.357_q,3.6231_q,0.47990_q, &
     &    0.12279_q,rtRS,ALFMRPA,ALFRSMRPA)
! jH QMC spin stiffness
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFMQM,ALFRSMQM)     
      ALFCRPA = -ALFMRPA
      ALFCQM  = -ALFMQM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      ECRPA = EURPA*(1._q-F*Z4)+EPRPA*F*Z4-ALFMRPA*F*(1._q-Z4)/FZZ
      ECQM  =  EUQM*(1._q-F*Z4)+EPQM*F*Z4-ALFMQM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
! jH Aufpassen: Potential 
      ECRSRPA = EURSRPA*(1._q-F*Z4)+EPRSRPA*F*Z4-ALFRSMRPA*F*(1._q-Z4)/FZZ
      ECRSQM  =  EURSQM*(1._q-F*Z4)+EPRSQM*F*Z4-ALFRSMQM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZETRPA = 4._q*(ZET**3)*F*(EPRPA-EURPA+ALFMRPA/FZZ)+FZ*(Z4*EPRPA-Z4*EURPA &
     &        -(1._q-Z4)*ALFMRPA/FZZ)
      ECZETQM = 4._q*(ZET**3)*F*(EPQM-EUQM+ALFMQM/FZZ)+FZ*(Z4*EPQM-Z4*EUQM &
     &        -(1._q-Z4)*ALFMQM/FZZ) 
      COMMRPA = ECRPA -RS*ECRSRPA/3._q-ZET*ECZETRPA
      COMMQM  = ECQM  -RS*ECRSQM/3._q-ZET*ECZETQM
      VCUP = COMMQM + ECZETQM - (COMMRPA + ECZETRPA)
      VCDN = COMMQM - ECZETQM - (COMMRPA - ECZETRPA)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
!      EC=EC*2
       EC=(ECQM-ECRPA)*2
       VCUP=VCUP*2
       VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_RPA_PLUS

    SUBROUTINE CORPBE_LDA_SR_0_15au(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS   = EURS 
       CALL gcor_xc_2(0.031091_q,-0.0089687_q,5.1486_q,1.3952503_q, & 
     & 0.2284487_q,0.0965737_q,1.6078169_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
! SR correlation energy so definiert das HEG limit
      EU = EMC - (ERPATOT - ERPASR)
      EURS = EMCRS -  (ERPATOTRS - ERPASRRS)

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-RS*ECRS/3._q 
      VCUP = COMM    
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_0_15au

    SUBROUTINE CORPBE_LDA_SR_0_5A(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS     = EURS
       CALL gcor_xc_2(0.031091_q,-0.0106268_q,5.1486_q,1.4094821_q, & 
     & 1.2291950_q,0.2579335_q,1.8066972_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
      EU = EMC - (ERPATOT - ERPASR)
      EURS = EMCRS -  (ERPATOTRS - ERPASRRS)

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-RS*ECRS/3._q
      VCUP = COMM
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_0_5A

    SUBROUTINE CORPBE_LDA_SR_1_0A(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS     = EURS
       CALL gcor_xc_2(0.031091_q,-0.0119450_q,5.1486_q,1.5999830_q, & 
     & 4.7029103_q,1.3622380_q,1.9748249_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
      EU = EMC -(ERPATOT - ERPASR)
      EURS = EMCRS -(ERPATOTRS - ERPASRRS)

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-RS*ECRS/3._q
      VCUP = COMM
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_1_0A

    SUBROUTINE CORPBE_LDA_SR_2_0A(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS     = EURS
       CALL gcor_xc_2(0.031091_q,0.0127716_q,5.1486_q,2.5109626_q, & 
     & 15.3666137_q,10.2843432_q,2.1207238_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
      EU = EMC - (ERPATOT - ERPASR)
      EURS = EMCRS - (ERPATOTRS - ERPASRRS) 

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-RS*ECRS/3._q
      VCUP = COMM
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_2_0A


!----------------------------------------------------------------------
    FUNCTION PBE_ALPHA(RS)
!----------------------------------------------------------------------
!  Alpha (spin stiffness)
!  from the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!  spin stiffness 
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))

      rtrs=dsqrt(rs)
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units

      PBE_ALPHA = -ALFM/FZZ*2
      RETURN
    END FUNCTION PBE_ALPHA

! jH RPA_ALPHA new routine for RPA spin stiffness 
!----------------------------------------------------------------------
    FUNCTION RPA_ALPHA(RS)
!----------------------------------------------------------------------
!  Alpha (spin stiffness)
!  from the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!  spin stiffness
!  jH this parametrisation has been already impl in CORPBE_LDA_RPA
!  jH parametrisation from Perdew, Wang, prB 45, 13244 (1992)
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))

      rtrs=dsqrt(rs)
      CALL gcor_xc(0.016887_q,0.028829_q,10.357_q,3.6231_q,0.47990_q, &
     &    0.12279_q,rtRS,ALFM,ALFRSM)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units

      RPA_ALPHA = -ALFM/FZZ*2
      RETURN
    END FUNCTION RPA_ALPHA


!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
    SUBROUTINE GCOR_XC(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
    END SUBROUTINE GCOR_XC

!----------------------------------------------------------------------
!#########################JUDITH###################################
!----------------------------------------------------------------------

    SUBROUTINE GCOR_XC_RPA(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! This is the rs interpolation formula of Perdew/Wang Phys.Rev.B 45 (1992) 13244
! for the RPA case (p = 0.75). See FORMULA (10) in article 
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*sqrt(rtrs))))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+7._q/2._q*B4*sqrt(rtrs)))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END SUBROUTINE GCOR_XC_RPA

!----------------------------------------------------------------------
!###########################JUDITH################################
!----------------------------------------------------------------------

    SUBROUTINE GCOR_XC_2(A,A1,B1,B2,B3,B4,p,rtrs,GG,GGRS)
! This is the rs interpolation formula of Perdew/Wang Phys.Rev.B 45 (1992) 13244
! for p = 2. See FORMULA (10) in article. This interpolation is used to
! parametrsize the short-range for v = (1/r)*Erfc(-mu r) RPA jellium energie --> Judith
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs**(2._q*p-1._q))))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+(2._q*p+2._q)*B4*rtrs**(2._q*p-1._q)))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END SUBROUTINE GCOR_XC_2

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
! Iann Gerber 18/11/04

!=====================================================================
! short range screened exchange functional
! EX_SR computes exchange-part with modified kernel
! interpolation formula of A. Savin and J. Toulouse 
!=====================================================================

    
    FUNCTION EX_SR(RS,RMU,IFLG)
GPACC !$ACC ROUTINE SEQ
GPACC USE moffload_struct_def, ONLY : ACC_ON_DEVICE, ACC_DEVICE_NOT_HOST

      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EX_SR)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.
      INTEGER :: IFLG      ! paramagnetic (IFLG=1) or ferromagnetic (IFLG=2)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      EX_=-CX(IFLG)/RS

! another IG bug: IFLAG replaced by IFLG
      IF (IFLG==2) THEN
         QFAC=(6._q*PI*PI)**(1._q/3._q)
      ELSE
         QFAC=(3._q*PI*PI)**(1._q/3._q)
      ENDIF

      RHO = 3/(4._q*PI)/RS**3
      QF=QFAC*(RHO)**(1._q/3._q)
      A=RMU/2._q/QF
!      FRAC=(1 - (8._q/3._q)*A*( SQRT(PI)*ERRF(1/(2*A)) + (2*A-4*A*A*A)*EXP(-1/(4*A*A)) &
!     &        -3*A + 4*A*A*A ))
!      EX_SR=EX_*FRAC
! IG Simple formula replaced by 

! Test on the value of A
      IF (A < 1E-9_q) THEN
! Limit for small A
         EX_SR=EX_
      ELSEIF (A <= 100._q) THEN
! Intermediate Values of A
#ifdef ACC_OFFLOAD
!(sm): only when compiling GPU code this construct will use standard ERF which has a GPU
!      implementation instead of the custom CPU implementation ERRF
      IF(ACC_ON_DEVICE(ACC_DEVICE_NOT_HOST)) THEN
         FRAC=(1.0_q - (8._q/3._q)*A*( SQRT(PI)*ERF(1._q/(2._q*A)) + &
     &        (2._q*A-4._q*A*A*A)*EXP(-1._q/(4._q*A*A))-3._q*A + 4._q*A*A*A ))
      ELSE
#endif
         FRAC=(1.0_q - (8._q/3._q)*A*( SQRT(PI)*ERRF(1._q/(2._q*A)) + &
     &        (2._q*A-4._q*A*A*A)*EXP(-1._q/(4._q*A*A))-3._q*A + 4._q*A*A*A ))
#ifdef ACC_OFFLOAD
      ENDIF
#endif
         EX_SR=EX_*FRAC
      ELSEIF (A <= 1.E+9_q) THEN
! Development for large A
         FRAC=1._q/36._q/A/A
         EX_SR=EX_*FRAC
      ELSE
         EX_SR=0._q
      ENDIF
!      B=EXP(-1._q/(4._q*A*A))-1._q
!      C=2._q*A*A*B+1._q/2._q
!      FRAC=1.0_q - 8._q/3._q*A*(SQRT(PI)*ERRF(1._q/(2._q*A))+2._q*A*(B-C))
!      EX_SR=EX_*FRAC
      RETURN
    END FUNCTION EX_SR
      

    FUNCTION VX_SR(RS,RMU,IFLG)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VX_SR)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.
      INTEGER :: IFLG      ! paramagnetic (IFLG=1) or ferromagnetic (IFLG=2)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      VX_=-CX(IFLG)/RS

! another IG bug: IFLAG replaced by IFLG
      IF (IFLG==2) THEN
         QFAC=(6._q*PI*PI)**(1._q/3._q)
      ELSE
         QFAC=(3._q*PI*PI)**(1._q/3._q)
      ENDIF

      RHO = 3/(4._q*PI)/RS**3
      QF=QFAC*(RHO)**(1._q/3._q)
      A=RMU/2._q/QF

!      FRAC= 1._q - 32*A*A*A*A*EXP(-1/(4*A*A)) -8._q*A*A + (64._q/3._q)*A*A*A*A
!      VX_SR=EX_SR(RS,RMU,IFLG)+(1._q/3._q)*VX_*FRAC
! IG Simple formula replaced by 

! Test on the value of A
      IF (A < 1E-9_q) THEN
! Limit for small A
         VX_SR=(4._q/3._q)*VX_
      ELSEIF (A <= 100._q) THEN
! Intermediate Values of A
         FRAC= (1._q - 32._q*A*A*A*A*(EXP(-1._q/(4._q*A*A))-1._q) -8._q*A*A)
         VX_SR=EX_SR(RS,RMU,IFLG)+(1._q/3._q)*VX_*FRAC
      ELSEIF (A <= 1.E+9_q) THEN
! Development for large A
         VX_SR=-RMU/48._q/PI/A/A/A
      ELSE
         VX_SR=0._q
      ENDIF
     RETURN
   END FUNCTION VX_SR

   FUNCTION EC_SR(RS,RMU,EC) ! No spin polarization
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.

      SAVE U1,U2,V1,A,BET,GAM
      DATA U1/1.0270741452992294_q/
      DATA U2/-0.230160617208092_q/
      DATA V1/0.6196884832404359_q/
      DATA A/3.2581_q/
      DATA BET/163.44_q/
      DATA GAM/4.7125_q/

      D=32._q/3._q/PI
! changed to a.u. instead of A, iG please check
      RHO = 3/(4._q*PI)/RS**3

!     G(0) from Burke, Perdew & Ernzerhof
      GRS32=(GAM+RS)**3._q/2._q
      GRS12=SQRT(GAM+RS)
      G0=D*(GRS32+BET)*EXP(-A*GRS12)

      C1N=U1*RS + U2*RS*RS
      C1D=1._q+V1*RS
      C1=C1N/C1D
      C2D=0.5_q*PI*RHO*RHO*(G0-0.5_q)
      C2=EC/C2D

      SCALE=1._q + C1*RMU + C2*RMU*RMU

      EC_SR=EC/SCALE
      RETURN
    END FUNCTION EC_SR

    FUNCTION VC_SR(RS,RMU,VC,EC) ! No spin polarization
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.

      SAVE U1,U2,V1,A,BET,GAM
      DATA U1/1.0270741452992294_q/
      DATA U2/-0.230160617208092_q/
      DATA V1/0.6196884832404359_q/
      DATA A/3.2581_q/
      DATA BET/163.44_q/
      DATA GAM/4.7125_q/

      D=32._q/3._q/PI
! changed to a.u. instead of A, iG please check
      RHO = 3/(4._q*PI)/RS**3

!     G(0) from Burke, Perdew & Ernzerhof
      GRS32=(GAM+RS)**3._q/2._q
      GRS12=SQRT(GAM+RS)
      G0=D*(GRS32+BET)*EXP(-A*GRS12)

      C1N=U1*RS + U2*RS*RS
      C1D=1._q+V1*RS
      C1=C1N/C1D
      C2D=0.5_q*PI*RHO*RHO*(G0-0.5_q)
      C2=EC/C2D

      SCALE1=1._q + C1*RMU + C2*RMU*RMU

      DRS=-RS/3._q/RHO
      DG0=-0.5_q*D*EXP(-A*GRS12)*(A*BET-3._q*GAM-3._q*RS+A*GRS32)*DRS/GRS12

      DC1=(U1+2._q*U2*RS + U2*V1*RS*RS)*DRS/C1D/C1D
      DC2D=0.5_q*PI*RHO*(2._q*G0+RHO*DG0-1._q)

      DC2=(VC*C2D-EC*DC2D)/C2D/C2D

      SCALE2=DC1*RMU + DC2*RMU*RMU

      VC_SR=(VC*SCALE1-EC*SCALE2)/(SCALE1*SCALE1)
      RETURN
    END FUNCTION VC_SR

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
            
!=====================================================================
! 
! Thomas-Fermi screened exchange functional
!  
!=====================================================================

    FUNCTION EX_SX(RS,RMU,IFLG)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      EX_=-CX(IFLG)/RS

!     RMU is the Thomas Fermi wave vector supplied in a.u., RS in a.u. as well
      RHO = 3/(4._q*PI)/(RS)**3

! here is the "global" implementation Walter Wolf adopted in FLAPW
!      QF_GLOBAL=RMU**2*PI/4
!      A=RMU/QF_GLOBAL

      QF_LOCAL=(3._q*PI*PI*RHO)**(1._q/3._q)
      A=RMU/QF_LOCAL

      FRAC = 1 - (4._q/3._q) * A * ATAN(2._q/A) &
           & - (A*A/6._q) * ( 1._q - ( A*A/4._q  + 3._q) * LOG(1 + (4._q/(A*A) )))
      EX_SX=EX_*FRAC
      RETURN
    END FUNCTION EX_SX

    FUNCTION VX_SX(RS,RMU,IFLG)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/

!     RMU is the Thomas Fermi wave vector supplied in a.u., RS in a.u. as well
      RHO = 3/(4._q*PI)/(RS)**3
      VX_=-CX(IFLG)/RS

! here is the "global" implementation Walter Wolf adopted in FLAPW
!      QF_GLOBAL=RMU**2*PI/4
!      A=RMU/QF_GLOBAL
!      FRAC = 1 - (4._q/3._q) * A * ATAN(2._q/A) &
!           & - (A*A/6._q) * ( 1._q - ( A*A/4._q  + 3._q) * LOG(1 + (4._q/(A*A) )))
!      VX_SX=VX_*FRAC

      QF_LOCAL=(3._q*PI*PI*RHO)**(1._q/3._q)
      A=RMU/QF_LOCAL

      FRAC = 1 + (1._q/2._q) * A*A - (1._q/8._q) * LOG( 1 + (4._q/(A*A) )) * A*A*A*A &
           &  - (1._q/2._q) * LOG( 1 + (4._q/(A*A))) *A*A
      VX_SX=EX_SX(RS,RMU,IFLG)+(1._q/3._q)*VX_*FRAC

      RETURN
    END FUNCTION VX_SX

! jP
!--------------------------------------------------------------------
!> The routines for the evaluation of the spin-polarized 
!> range-separated exchange-correlation energy in the LDA,
!> LSDSR, ECORRLR, VCORRLR, EXCHANGELR, VEXCHANGELR, ECPW, and GPW,
!> as well as the functions
!> G0, G0D, DPOL, DPOLD, QRPA, QRPAD
!> are downloaded from Paola Gori-Giorgi's webpage
!> http://www.lct.jussieu.fr/pagesperso/gori/elegas.html
!> [Phys. Rev. B 73, 155111 (2006)].
!> They have been adapted to f90 by jP (04/25/08).
!--------------------------------------------------------------------
      SUBROUTINE lsdsr(rs,z,mu,excsr,vxcsrup,vxcsrdown)
!   Hartree atomic units used
!   for given density parameter 'rs', relative spin polarization 'z'= (nu -nd)/n
!   and cutoff parameter 'mu' 
!   gives the complementary  short-range exchange-correlation
!   energy  (i.e., xc energy of jellium minus xc energy of long-range
!   interacting electron gas) => 'excsr'
!   and the corresponding exchange-correlation potentials for
!   spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
!   from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601343
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,excsr,vxcsrup,vxcsrdown
        REAL(q) :: eclr,exlr,ec,ecd,ecz,ex
        REAL(q) :: vclrup,vclrdown,vxlrup,vxlrdown
        REAL(q) :: vxup,vxdown,vcup,vcdown
        REAL(q) :: alpha,cf

        alpha=(4.0_q/9.0_q/pi)**(1.0_q/3.0_q)
        cf=1.0_q/alpha
        
        ex=-3.0_q*cf/rs/8.0_q/pi*((1.0_q+z)**(4.0_q/3.0_q)+&
        &     (1.0_q-z)**(4.0_q/3.0_q))
        
        vxup=-(1.0_q+z)**(1.0_q/3.0_q)*(3.0_q/2.0_q/PI)**(2.0_q/3.0_q)/rs
        vxdown=-(1.0_q-z)**(1.0_q/3.0_q)*(3.0_q/2.0_q/PI)**(2.0_q/3.0_q)/rs
        
        call ecPW(rs,z,ec,ecd,ecz)
        vcup=ec-rs/3.0_q*ecd-(z-1.0_q)*ecz
        vcdown=ec-rs/3.0_q*ecd-(z+1.0_q)*ecz
        
        call exchangelr(rs,z,mu,exlr)
        call vexchangelr(rs,z,mu,vxlrup,vxlrdown)
        
        call ecorrlr(rs,z,mu,eclr)
        call vcorrlr(rs,z,mu,vclrup,vclrdown)
        
        excsr=ex+ec-(exlr+eclr)
        vxcsrup=vxup+vcup-(vxlrup+vclrup)
        vxcsrdown=vxdown+vcdown-(vxlrdown+vclrdown)

        RETURN
        END SUBROUTINE lsdsr

!----------------------------------------------------------------------
!> Routine to calculate the spin-polarized
!> short-range correlation energy 
!> in the local density approximation -> see routine
!> LSDSR; 
!> Adaptation done by jP (04/25/08)
!----------------------------------------------------------------------
      SUBROUTINE ecorrsr(rs,z,mu,ecsr)
!   Hartree atomic units used
!   for given density parameter 'rs', relative spin polarization 'z'= (nu -nd)/n
!   and cutoff parameter 'mu' 
!   gives the complementary  short-range exchange-correlation
!   energy  (i.e., xc energy of jellium minus xc energy of long-range
!   interacting electron gas) => 'excsr'
!   and the corresponding exchange-correlation potentials for
!   spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
!   from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601343
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,ecsr
        REAL(q) :: eclr,ec,ecd,ecz

        
        call ecPW(rs,z,ec,ecd,ecz)
        call ecorrlr(rs,z,mu,eclr)
        
        ecsr=ec-eclr
        RETURN
      END SUBROUTINE ecorrsr

!----------------------------------------------------------------------
!> Routine to calculate the spin-polarized
!> short-range correlation potentials
!> in the local density approximation -> see routine
!> LSDSR; 
!> Adaptations w.r.t. VASP done by jP (04/25/08)
!----------------------------------------------------------------------
      SUBROUTINE vcorrsr(rs,z,mu,vcsrup,vcsrdown)
!   Hartree atomic units used
!   for given density parameter 'rs', relative spin polarization 'z'= (nu -nd)/n
!   and cutoff parameter 'mu' 
!   gives the complementary  short-range exchange-correlation
!   energy  (i.e., xc energy of jellium minus xc energy of long-range
!   interacting electron gas) => 'excsr'
!   and the corresponding exchange-correlation potentials for
!   spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
!   from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601343
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,vcsrup,vcsrdown
        REAL(q) :: ec,ecd,ecz
        REAL(q) :: vclrup,vclrdown
        REAL(q) :: vcup,vcdown
        
        call ecPW(rs,z,ec,ecd,ecz)
        vcup=ec-rs/3.0_q*ecd-(z-1.0_q)*ecz
        vcdown=ec-rs/3.0_q*ecd-(z+1.0_q)*ecz
        
        call vcorrlr(rs,z,mu,vclrup,vclrdown)
        
        vcsrup=vcup-vclrup
        vcsrdown=vcdown-vclrdown
        RETURN
      END SUBROUTINE vcorrsr 
     
      SUBROUTINE ecorrlr(rs,z,mu,eclr)
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff parameter mu 
!   gives the correlation energy of the LR gas
!    => eclr
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,eclr,ec,ecd,ecz
        REAL(q) :: alpha,cf,phi
        REAL(q) :: d2anti,d3anti
        REAL(q) :: coe2,coe3,coe4,coe5
        REAL(q) :: a1,a2,a3,a4,b0
        REAL(q) :: q1a,q2a,q3a,t1a,t2a,t3a,adib

        alpha=(4.0_q/9.0_q/pi)**(1.0_q/3.0_q)
        cf=1.0_q/alpha
        
        phi=((1.0_q+z)**(2.0_q/3.0_q)+(1.0_q-z)**(2.0_q/3.0_q))/2.0_q
!  parameters from the fit
        adib   = 0.784949_q   
        q1a    = -0.388_q   
        q2a    = 0.676_q   
        q3a    = 0.547_q   
        t1a    = -4.95_q   
        t2a    = 1._q    
        t3a    = 0.31_q   
        
        b0=adib*rs
        
        d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
        d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3
        
        coe2=-3._q/8._q/rs**3*(1._q-z**2)*(gg0(rs)-0.5_q)
        
        coe3=-(1._q-z**2)*gg0(rs)/(sqrt(2._q*pi)*rs**3)
        
        if(abs(z).eq.1._q) then
           
           coe4=-9._q/64._q/rs**3*(dpol(rs) &
                &        -cf**2*2**(5._q/3._q)/5._q/rs**2) 
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*dpol(rs)
           
        else
           
           coe4=-9._q/64._q/rs**3*(((1._q+z)/2._q)**2* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+ &
                &        (1.-z**2)*d2anti-cf**2/10._q*((1._q+z)**(8._q/3._q) &
                &        +(1.-z)**(8._q/3._q))/rs**2)
           
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(((1._q+z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d3anti)
        endif
        
        call ecPW(rs,z,ec,ecd,ecz)
        
        a1=4._q*b0**6*coe3+b0**8*coe5
        a2=4._q*b0**6*coe2+b0**8*coe4+6._q*b0**4*ec
        a3=b0**8*coe3
        a4=b0**6*(b0**2*coe2+4._q*ec)
        
        eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+ &
             &     a4*mu**6+b0**8*mu**8*ec)/((1._q+b0**2*mu**2)**4)
        
        RETURN
      END SUBROUTINE ecorrlr

      SUBROUTINE vcorrlr(rs,z,mu,vclrup,vclrdown)
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff mu it gives the correlation LSD potential for LR interaction
!    => vclrup (spin-up electrons), vclrdown (spin-down electrons)
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,eclr,eclrrs,eclrz,vclrup,vclrdown
        REAL(q) :: ec,ecd,ecz
        REAL(q) :: alpha,cf,phi
        REAL(q) :: d2anti,d3anti
        REAL(q) :: d2antid,d3antid,x
        REAL(q) :: coe2,coe3,coe4,coe5
        REAL(q) :: coe2rs,coe3rs,coe4rs,coe5rs
        REAL(q) :: coe2z,coe3z,coe4z,coe5z
        REAL(q) :: a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs
        REAL(q) :: b0rs,a1z,a2z,a3z,a4z,a5z,b0z
        REAL(q) :: q1a,q2a,q3a,t1a,t2a,t3a,adib
        
        alpha=(4._q/9._q/pi)**(1._q/3._q)
        cf=1._q/alpha
        
        phi=((1._q+z)**(2._q/3._q)+(1._q-z)**(2._q/3._q))/2._q
!  parameters from the fit
        adib   = 0.784949_q   
        q1a    = -0.388_q   
        q2a    = 0.676_q   
        q3a    = 0.547_q   
        t1a    = -4.95_q   
        t2a    = 1._q    
        t3a    = 0.31_q   
        
        b0=adib*rs
        
        d2anti=(q1a+q2a*rs)*exp(-q3a*rs)/rs
        d3anti=(t1a+t2a*rs)*exp(-t3a*rs)/rs**2
        
        d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/ &
             &    rs**2)*exp(-q3a*rs)
        d3antid=-((rs*t2a*(1 + rs*t3a) + t1a*(2 + rs*t3a))/ &
             &    rs**3)*exp(-rs*t3a)
        
        coe2=-3._q/8._q/rs**3*(1._q-z**2)*(gg0(rs)-0.5_q)
        coe2rs=-3._q/8._q/rs**3*(1._q-z**2)*gg0d(rs)+ &
             &     9._q/8._q/rs**4*(1._q-z**2)*(gg0(rs)-0.5_q)
        coe2z=-3._q/8._q/rs**3*(-2._q*z)*(gg0(rs)-0.5_q)

        coe3=-(1._q-z**2)*gg0(rs)/(sqrt(2._q*pi)*rs**3)
        coe3rs=-(1._q-z**2)*gg0d(rs)/(sqrt(2._q*pi)*rs**3)+ &
             &    3._q*(1._q-z**2)*gg0(rs)/(sqrt(2._q*pi)*rs**4) 
        coe3z=2._q*z*gg0(rs)/(sqrt(2._q*pi)*rs**3)
        
        if(abs(z).eq.1._q) then
           
           coe4=-9._q/64._q/rs**3*(dpol(rs) &
                &        -cf**2*2**(5._q/3._q)/5._q/rs**2)
           coe4rs=-3._q/rs*coe4-9._q/64._q/rs**3*(dpold(rs) &
                &        +2._q*cf**2*2**(5._q/3._q)/5._q/rs**3)
           coe4z=-9._q/64._q/rs**3*(dpol(rs)-rs/6._q*dpold(rs)-2._q*d2anti &
                &       -4._q/15._q/rs**2*cf**2*2._q**(5._q/3._q))*z
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*dpol(rs)
           coe5rs=-3._q/rs*coe5-9._q/40._q/(sqrt(2._q*pi)*rs**3)*dpold(rs)
           coe5z=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(dpol(rs)-rs/6._q* &
                &       dpold(rs)-2._q*d3anti)*z
           
        else
           
           coe4=-9._q/64._q/rs**3*(((1._q+z)/2._q)**2* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+ &
                &        (1.-z**2)*d2anti-cf**2/10._q*((1._q+z)**(8._q/3._q) &
                &        +(1.-z)**(8._q/3._q))/rs**2)
           coe4rs=-3._q/rs*coe4-9._q/64._q/rs**3*( &
                &        ((1._q+z)/2._q)**(5._q/3._q)*dpold(rs*(2/(1._q+z))** &
                &        (1._q/3._q))+((1._q-z)/2._q)**(5._q/3._q)* &
                &        dpold(rs*(2/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d2antid+cf**2/5._q*((1._q+z)**(8._q/3._q) &
                &        +(1._q-z)**(8._q/3._q))/rs**3)
           coe4z=-9._q/64._q/rs**3*(1._q/2._q*(1._q+z)* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))-1._q/2._q*(1._q-z)* &
                &        dpol(rs*(2/(1._q-z))**(1._q/3._q))-rs/6._q* &
                &        ((1._q+z)/2._q)**(2._q/3._q)*dpold(rs*(2/(1._q+z)) &
                &        **(1._q/3._q))+rs/6._q*((1._q-z)/2._q)**(2._q/3._q) &
                &        *dpold(rs*(2/(1._q-z))**(1._q/3._q))-2._q*z*d2anti- &
                &        4._q/15._q/rs**2*cf**2*((1._q+z)**(5._q/3._q)- &
                &        (1._q-z)**(5._q/3._q)))
           
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(((1._q+z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d3anti)
           coe5rs=-3._q/rs*coe5-9._q/(40._q*sqrt(2._q*pi)*rs**3)*( &
                &        ((1._q+z)/2._q)**(5._q/3._q)*dpold(rs*(2/(1._q+z))** &
                &        (1._q/3._q))+((1._q-z)/2._q)**(5._q/3._q)* &
                &        dpold(rs*(2/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d3antid)
           coe5z=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(1._q/2._q*(1._q+z)* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))-1._q/2._q*(1._q-z)* &
                &        dpol(rs*(2/(1._q-z))**(1._q/3._q))-rs/6._q* &
                &        ((1._q+z)/2._q)**(2._q/3._q)*dpold(rs*(2/(1._q+z)) &
                &        **(1._q/3._q))+rs/6._q*((1._q-z)/2._q)**(2._q/3._q) &
                &        *dpold(rs*(2/(1._q-z))**(1._q/3._q))-2._q*z*d3anti)
           
        endif
        
        call ecPW(rs,z,ec,ecd,ecz)
        
        a1=4._q*b0**6*coe3+b0**8*coe5
        a1rs=24._q*adib*b0**5*coe3+4._q*b0**6*coe3rs+8._q*adib*b0**7* &
             &     coe5+b0**8*coe5rs
        a1z=4._q*b0**6*coe3z+b0**8*coe5z
        
        a2=4._q*b0**6*coe2+b0**8*coe4+6._q*b0**4*ec
        a2rs=24._q*adib*b0**5*coe2+4._q*b0**6*coe2rs+8._q*adib*b0**7* &
             &     coe4+b0**8*coe4rs+24._q*adib*b0**3*ec+6._q*b0**4*ecd
        a2z=4._q*b0**6*coe2z+b0**8*coe4z+6._q*b0**4*ecz
        
        a3=b0**8*coe3
        a3rs=8._q*adib*b0**7*coe3+b0**8*coe3rs
        a3z=b0**8*coe3z
        
        a4=b0**6*(b0**2*coe2+4._q*ec)
        a4rs=8._q*adib*b0**7*coe2+b0**8*coe2rs+24._q*adib*b0**5*ec+ &
             &     4._q*b0**6*ecd
        a4z=b0**6*(b0**2*coe2z+4._q*ecz)
        
        a5=b0**8*ec
        a5rs=8._q*adib*b0**7*ec+b0**8*ecd
        a5z=b0**8*ecz
        
        x=mu*sqrt(rs)/phi
        
        eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+ &
             &     a4*mu**6+a5*mu**8)/((1._q+b0**2*mu**2)**4)
        
        eclrrs=-4._q/(1._q+b0**2*mu**2)*2._q*adib*b0*mu**2*eclr+ &
             &     1._q/((1._q+b0**2*mu**2)**4)*(phi**2*mu/(2._q*sqrt(rs)) &
             &     *Qrpad(x)+ &
             &     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)
        
        if(z.eq.1._q) then
           vclrup=eclr-rs/3._q*eclrrs
           vclrdown=0._q
        elseif(z.eq.-1._q) then
           vclrup=0._q
           vclrdown=eclr-rs/3._q*eclrrs
        else
           
           eclrz=(phi**2*((1._q+z)**(-1._q/3._q)-(1._q-z)**(-1._q/3._q)) &
                &        *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1._q+z)**(-1._q/3._q) &
                &        -(1._q-z)**(-1._q/3._q))/3._q+ &
                &        a1z*mu**3+a2z*mu**4+a3z*mu**5+ &
                &        a4z*mu**6+a5z*mu**8)/((1._q+b0**2*mu**2)**4)
           
           vclrup=eclr-rs/3._q*eclrrs-(z-1._q)*eclrz
           vclrdown=eclr-rs/3._q*eclrrs-(z+1._q)*eclrz
        endif
        RETURN
      END SUBROUTINE vcorrlr

      FUNCTION gg0(x)
        USE prec
        USE constant
!   on-top pair-distribution function
!   Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
!   x -> rs
        IMPLICIT NONE
        REAL(q) :: gg0
        REAL(q) ::  C0f,D0f,E0f,F0f,x

        gg0              = 0.0_q
        C0f             = 0.0819306_q  
        D0f             = 0.752411_q    
        E0f             = -0.0127713_q   
        F0f             = 0.00185898_q   
        gg0=(1._q-(0.7317_q-D0f)*x+C0f*x**2+E0f*x**3+ &
             &     F0f*x**4)*exp(-abs(D0f)*x)/2._q
      RETURN
    END FUNCTION gg0

    FUNCTION gg0d(rs)
        USE prec
!   derivative of on-top pair-distribution function
!   Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
        IMPLICIT NONE
        REAL(q) :: gg0d
        REAL(q) :: Bgg0,Cgg0,Dgg0,Egg0,Fgg0,rs

        gg0d             = 0.0_q
        Cgg0             = 0.0819306_q    
        Fgg0             = 0.752411_q     
        Dgg0             = -0.0127713_q   
        Egg0             = 0.00185898_q
        Bgg0             =0.7317_q-Fgg0
        gg0d=(-Bgg0+2*Cgg0*rs+3*Dgg0*rs**2+4*Egg0*rs**3)/2._q*exp(-Fgg0*rs) &
             &   - (Fgg0*(1 - Bgg0*rs + Cgg0*rs**2 + Dgg0*rs**3 + Egg0*rs**4))/ &
             &   2._q*exp(-Fgg0*rs)
        RETURN
      END FUNCTION gg0d

      FUNCTION dpol(rs)
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) dpol
        REAL(q) :: cf,rs,p2p,p3p

        cf=(9._q*pi/4._q)**(1._q/3._q)
        p2p    = 0.04_q 
        p3p    = 0.4319_q   
        dpol=2._q**(5._q/3._q)/5._q*cf**2/rs**2*(1._q+(p3p-0.454555_q)*rs) &
             &     /(1._q+p3p*rs+p2p*rs**2)
      RETURN
    END FUNCTION dpol
    
    FUNCTION dpold(rs)
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) dpold
        REAL(q) ::  cf,rs,p2p,p3p

        cf=(9._q*pi/4._q)**(1._q/3._q)
        p2p    = 0.04_q 
        p3p    = 0.4319_q   
        dpold=2._q**(5._q/3._q)/5._q*cf**2* &
             & (-2. + (0.454555 - 4.*p3p)*rs + &
             &    (-4.*p2p + &
             &       (0.90911 - 2.*p3p)*p3p)*rs**2 &
             &      + p2p*(1.363665 - 3.*p3p)* &
             &     rs**3)/ &
             &  (rs**3*(1. + p3p*rs + p2p*rs**2)**2)
      RETURN
    END FUNCTION dpold

    FUNCTION Qrpa(x)
      USE prec
      USE constant
      IMPLICIT NONE
        REAL(q) Qrpa
      REAL(q) :: a2,b2,c2,d2,x,Acoul

      Acoul=2._q*(log(2._q)-1._q)/pi**2
      a2              = 5.84605_q 
      c2              = 3.91744_q 
      d2              = 3.44851_q
      b2=d2-3._q/(2._q*pi*Acoul)*(4._q/(9._q*pi))**(1._q/3._q)
      Qrpa=Acoul*log((1._q+a2*x+b2*x**2+c2*x**3)/(1._q+a2*x+d2*x**2))
      RETURN
    END FUNCTION Qrpa

    FUNCTION Qrpad(x)
      USE prec
      USE constant
      IMPLICIT NONE
        REAL(q) Qrpad
      REAL(q) :: a2,b2,c2,d2,x,Acoul
      
      Acoul=2._q*(log(2._q)-1._q)/pi**2
      a2              = 5.84605_q 
      c2              = 3.91744_q 
      d2              = 3.44851_q
      b2=d2-3._q/(2._q*pi*Acoul)*(4._q/(9._q*pi))**(1._q/3._q)
      Qrpad=Acoul*((x*(b2*(2._q + a2*x) + &
           &      c2*x*(3._q + 2._q*a2*x) + &
           &      d2*(-2._q - a2*x + c2*x**3)))/ &
           &  ((1._q + a2*x + d2*x**2)* &
           &    (1._q + a2*x + b2*x**2 + c2*x**3)))
      RETURN
    END FUNCTION Qrpad
    
    SUBROUTINE exchangelr(rs,z,mu,exlr)
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff mu it gives the exchange energy of the LR gas
!    => exlr
      IMPLICIT NONE
      REAL(q) :: rs,z,mu,exlr
      REAL(q) :: alpha,fx,y

      alpha=(4._q/9._q/pi)**(1._q/3._q)
      if(abs(z).eq.1._q) then
         y=mu*alpha*rs/2._q/2._q**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         exlr=mu*fx
      else
         y=mu*alpha*rs/2._q/(1.+z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         exlr=(1._q+z)*mu*fx/2._q
         y=mu*alpha*rs/2._q/(1.-z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         exlr=exlr+(1._q-z)*mu*fx/2._q
      endif
      RETURN
    END SUBROUTINE exchangelr
    
    SUBROUTINE vexchangelr(rs,z,mu,vxlrup,vxlrdown)
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff mu it gives the exchange LSD potential for LR interaction
!    => vxlrup (spin-up electrons), vxlrdown (spin-down electrons)
      IMPLICIT NONE
      REAL(q) :: rs,z,mu,vxlrup,vxlrdown
      REAL(q) :: alpha,fx,fx1,y,exlr,derrs,derz

      alpha=(4._q/9._q/pi)**(1._q/3._q)
      if(z.eq.1._q) then
         vxlrup=(rs*alpha*mu**2)/ &
              &   (2**(1._q/3._q)*pi) - (rs*alpha*mu**2)/(2**(1._q/3._q)*pi)* &
              &     exp(-2**(2._q/3._q)/(rs**2*alpha**2*mu**2)) - &
              &  (mu*ERRF(2**(1._q/3._q)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrdown=0._q
      elseif(z.eq.-1._q) then
         vxlrdown=(rs*alpha*mu**2)/ &
              &   (2**(1._q/3._q)*pi) - (rs*alpha*mu**2)/(2**(1._q/3._q)*pi)* &
              &     exp(-2**(2._q/3._q)/(rs**2*alpha**2*mu**2)) - &
              &  (mu*ERRF(2**(1._q/3._q)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrup=0._q
      else       
         y=mu*alpha*rs/2._q/(1.+z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         fx1=(3._q*(1 + (-4._q + 4._q*exp(-1._q/(4._q*y**2)))*y**2))/pi
         derrs=1._q/4._q*(1._q+z)**(2._q/3._q)*mu**2*alpha*fx1
         derz=1._q/2._q*mu*fx-1._q/6._q*fx1*mu*y
         vxlrup=rs/3._q*derrs+(z-1._q)*derz
         vxlrdown=rs/3._q*derrs+(z+1._q)*derz
         
         y=mu*alpha*rs/2._q/(1.-z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         fx1=(3._q*(1 + (-4._q + 4._q*exp(-1._q/(4._q*y**2)))*y**2))/pi
         derrs=1._q/4._q*(1._q-z)**(2._q/3._q)*mu**2*alpha*fx1
         derz=-1._q/2._q*mu*fx+1._q/6._q*fx1*mu*y
         vxlrup=vxlrup+rs/3._q*derrs+(z-1._q)*derz
         vxlrdown=vxlrdown+rs/3._q*derrs+(z+1._q)*derz
         
         call exchangelr(rs,z,mu,exlr)
         vxlrup=exlr-vxlrup
         vxlrdown=exlr-vxlrdown
      endif
      RETURN
    END SUBROUTINE vexchangelr

!-----------------------------------------------------------------------
!> correlation energy and its derivative w.r.t. rs and z at mu=infinity
!> Perdew & Wang PRB 45, 13244 (1992)
!-----------------------------------------------------------------------
    SUBROUTINE ecPW(x,y,ec,ecd,ecz)
      USE prec
      USE constant
! in Hartree; ec=ec(rs,zeta)
! x -> rs; y -> zeta
!   ecd is d/drs ec
!   ecz is d/dz ec
      IMPLICIT NONE
      REAL(q) :: f02,ff,x,y,ec,ecd,ec0,ec0d,ec1,ec1d
      REAL(q) :: aaa,G,Gd,alfac,alfacd,ecz
      
      f02=4._q/(9._q*(2._q**(1._q/3._q)-1._q))
      
      ff=((1._q+y)**(4._q/3._q)+(1._q-y)**(4._q/3._q)- &
           &     2._q)/(2._q**(4._q/3._q)-2._q)
      
      aaa=(1._q-log(2._q))/pi**2
      call  GPW(x,aaa,0.21370_q,7.5957_q,3.5876_q,&
           &     1.6382_q,0.49294_q,G,Gd)
      ec0=G
      ec0d=Gd
      
      aaa=aaa/2._q
      call GPW(x,aaa,0.20548_q,14.1189_q,6.1977_q,&
           &     3.3662_q,0.62517_q,G,Gd)
      ec1=G
      ec1d=Gd
      call GPW(x,0.016887_q,0.11125_q,10.357_q,3.6231_q,&
           &     0.88026_q,0.49671_q,G,Gd)
      alfac=-G
      alfacd=-Gd
      
      ec=ec0+alfac*ff/f02*(1._q-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1._q-y**4)+(ec1d-ec0d)* &
           &     ff*y**4
      ecz=alfac*(-4._q*y**3)*ff/f02+alfac*(1._q-y**4)/f02* &
           &     4._q/3._q*((1._q+y)**(1._q/3._q)-(1._q-y)**(1._q/3._q))/ &
           &     (2._q**(4._q/3._q)-2._q)+(ec1-ec0)*(4._q*y**3*ff+ &
           &     4._q/3._q*((1._q+y)**(1._q/3._q)-(1._q-y)**(1._q/3._q))/ &
           &     (2._q**(4._q/3._q)-2._q)*y**4)
      
      RETURN
    END SUBROUTINE ecPW
    
    SUBROUTINE GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd)
      USE prec
!   Gd is d/drs G
      IMPLICIT NONE
      REAL(q) ::  G,Gd,Ac,alfa1,beta1,beta2,beta3,beta4,x
      G=-2._q*Ac*(1._q+alfa1*x)*dlog(1._q+1._q/(2._q* &
           &     Ac*(beta1*x**0.5_q+ &
           &     beta2*x+beta3*x**1.5_q+beta4*x**2)))
      Gd=(1._q+alfa1*x)*(beta2+beta1/(2._q*sqrt(x))+3._q*beta3* &
           &     sqrt(x)/2._q+2._q*beta4*x)/((beta1*sqrt(x)+beta2*x+ &
           &     beta3*x**(3._q/2._q)+beta4*x**2)**2*(1._q+1._q/ &
           &     (2._q*Ac*(beta1*sqrt(x)+beta2*x+beta3*x**(3._q/2._q)+ &
           &     beta4*x**2))))-2._q*Ac*alfa1*dlog(1._q+1._q/(2._q*Ac* &
           &     (beta1*sqrt(x)+beta2*x+beta3*x**(3._q/2._q)+ &
           &     beta4*x**2)))
      RETURN
    END SUBROUTINE GPW
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE PW91LDA(RS,ZET,EC,DECRS,DECDZ,D2ECRS,D2ECRSZ)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(PW91LDA)
      REAL(q) :: RS,ZET,EC,DECRS,DECDZ,D2ECRS,D2ECRSZ
      REAL(q) :: ONE,FOUR,THRD,THRDM,THRD2,SIXTHM,THRD4,GAM
      REAL(q) :: FZZ,GAMMA,BET,ETA,DELT,RTRS,EU,EURS,D2EURS
      REAL(q) :: EP,EPRS,D2EPRS,ALFM,ALFRSM,D2ALFRSM
      REAL(q) :: OPZ,OMZ,Z3,Z4,F,GZ,HZ,DFZ,DGZ,DHZ
      parameter(ONE=1.0_q)
      parameter(four=4.0_q)
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
GPACC !$ACC ROUTINE(GCOR2MGGA) SEQ

      rtrs=dsqrt(rs)
      CALL gcor2mgga(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS,D2EURS)
      CALL gcor2mgga(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS,D2EPRS)
      CALL gcor2mgga(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM,D2ALFRSM)

      OPZ = ONE + ZET
      OMZ = ONE - ZET

      Z3 = ZET**3
      Z4 = Z3*ZET
      F = (OPZ**THRD4+OMZ**THRD4-2._q)/GAM
      GZ = F*Z4
      HZ = (F-GZ)/FZZ

      EC = EU*(ONE-GZ) + EP*GZ - ALFM*HZ
      DECRS = EURS*(ONE-GZ) + EPRS*GZ - ALFRSM*HZ
      D2ECRS = D2EURS*(ONE-GZ) + D2EPRS*GZ - D2ALFRSM*HZ

      DFZ = THRD4*(OPZ**THRD - OMZ**THRD)/GAM
      DGZ = DFZ*Z4 + FOUR*F*Z3
      DHZ = (DFZ-DGZ)/FZZ

      DECDZ = (EP-EU)*DGZ - ALFM*DHZ
      D2ECRSZ = (EPRS-EURS)*DGZ - ALFRSM*DHZ

      END SUBROUTINE PW91LDA

      SUBROUTINE GCOR2MGGA(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS,D2GGRS)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(GCOR2MGGA)
      REAL(q) :: A,A1,B1,B2,B3,B4,RTRS,GG,GGRS,D2GGRS
      REAL(q) :: Q0,Q1,Q2,Q3,DQ0,DQ2,D2Q1,D2Q2,DQ1
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2

      DQ0 = -2._q*A*A1
      DQ1 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      DQ2 = -DQ1/(Q1+Q1**2)
      GGRS = DQ0*Q2+Q0*DQ2


      D2Q1 = A*(3._q*B3/(2._q*RTRS)-B1/(2._q*RTRS**3)+4._q*B4)
      D2Q2 = (2._q*Q1+1._q)*DQ2**2 - D2Q1/(Q1+Q1**2)
      D2GGRS = 2._q*DQ0*DQ2 + Q0*D2Q2
      RETURN
      END SUBROUTINE GCOR2MGGA

      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORLSD)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!     IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q),PARAMETER :: GAM=0.5198421E0_q,FZZ=1.709921E0_q
      REAL(q),PARAMETER :: THRD=0.333333333333E0_q,THRD4=1.333333333333E0_q
GPACC !$ACC ROUTINE (GCOR) SEQ
      F = ((1.E0_q+ZET)**THRD4+(1.E0_q-ZET)**THRD4-2.E0_q)/GAM
      CALL GCOR(0.0310907E0_q,0.21370E0_q,7.5957E0_q,3.5876E0_q,1.6382E0_q, &
     &    0.49294E0_q,1.00E0_q,RS,EU,EURS)
      CALL GCOR(0.01554535E0_q,0.20548E0_q,14.1189E0_q,6.1977E0_q,3.3662E0_q, &
     &    0.62517E0_q,1.00E0_q,RS,EP,EPRS)
      CALL GCOR(0.0168869E0_q,0.11125E0_q,10.357E0_q,3.6231E0_q,0.88026E0_q, &
     &    0.49671E0_q,1.00E0_q,RS,ALFM,ALFRSM)
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.E0_q-F*Z4)+EP*F*Z4-ALFM*F*(1.E0_q-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.E0_q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.E0_q-Z4)/FZZ
      FZ = THRD4*((1.E0_q+ZET)**THRD-(1.E0_q-ZET)**THRD)/GAM
      ECZET = 4.E0_q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1.E0_q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.E0_q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END SUBROUTINE CORLSD

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE CORLSD_SIMD(RS, ZET, EC, VCUP, VCDN, ECRS, ECZET, ALFC, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!     IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS, ZET
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EC, VCUP, VCDN, ECRS, ECZET, ALFC
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EU, EURS, EP, EPRS, ALFM, ALFRSM
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: GAM = 0.5198421E0_q, GAM_INV = 1.0_q / GAM, FZZ = 1.709921E0_q, FZZ_INV = 1.0_q / FZZ
      REAL(q), PARAMETER                 :: THRD = 0.333333333333E0_q, THRD4 = 1.333333333333E0_q

      CALL GCOR_SIMD(0.0310907E0_q, 0.21370E0_q, 7.5957E0_q, 3.5876E0_q, 1.6382E0_q, 0.49294E0_q, 1.00E0_q, RS, EU, EURS, MASK0)

      CALL GCOR_SIMD(0.01554535E0_q, 0.20548E0_q, 14.1189E0_q, 6.1977E0_q, 3.3662E0_q, 0.62517E0_q, 1.00E0_q, RS, EP, EPRS, MASK0)

      CALL GCOR_SIMD(0.0168869E0_q, 0.11125E0_q, 10.357E0_q, 3.6231E0_q, 0.88026E0_q, 0.49671E0_q, 1.00E0_q, RS, ALFM, ALFRSM, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, F, Z4, TEMP_3, TEMP_4, FZ, COMM)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = 1.0_q + ZET%X(II)
         TEMP_1 = TEMP_1 * TEMP_1
         TEMP_1 = TEMP_1 * TEMP_1
         TEMP_2 = 1.0_q - ZET%X(II)
         TEMP_2 = TEMP_2 * TEMP_2
         TEMP_2 = TEMP_2 * TEMP_2
         F = (TEMP_1 + TEMP_2 - 2.0_q) * GAM_INV
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
         TEMP_1 = ZET%X(II) * ZET%X(II)
         Z4 = TEMP_1 * TEMP_1
         TEMP_1 = F * Z4
         TEMP_2 = 1.0_q - TEMP_1
         TEMP_3 = (1.0_q - Z4) * FZZ_INV
         IF (MASK0%X(II)) THEN
            TEMP_4 = F * TEMP_3
            ALFC%X(II) = -ALFM%X(II)
            EC%X(II) = EU%X(II) * TEMP_2 + EP%X(II) * TEMP_1 - ALFM%X(II) * TEMP_4
!  ENERGY DONE. NOW THE POTENTIAL:
            ECRS%X(II) = EURS%X(II) * TEMP_2 + EPRS%X(II) * TEMP_1 - ALFRSM%X(II) * TEMP_4
         ENDIF
         FZ = THRD4 * (SIMD_POW(1.0_q + ZET%X(II), THRD) - SIMD_POW(1.0_q - ZET%X(II), THRD)) * GAM_INV
         IF (MASK0%X(II)) THEN
            ECZET%X(II) = 4.0_q * ZET%X(II) * ZET%X(II) * ZET%X(II) * F * (EP%X(II) - EU%X(II) + ALFM%X(II) * FZZ_INV) + &
     &           FZ * (Z4 * EP%X(II) - Z4 * EU%X(II) - ALFM%X(II) * TEMP_3)
            COMM = EC%X(II) - RS%X(II) * ECRS%X(II) * THRD - ZET%X(II) * ECZET%X(II)
            VCUP%X(II) = COMM + ECZET%X(II)
            VCDN%X(II) = COMM - ECZET%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CORLSD_SIMD
#endif

!=======================================================================
!
!> original GCOR used in PW91 routines, to interpolate
!> LSD correlation energy, as given by (10) of
!> J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!> K. Burke, May 11, 1996.
!
!=======================================================================

      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(GCOR)

      P1 = P + 1.E0_q
      Q0 = -2.E0_q*A*(1.E0_q+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.E0_q*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = LOG(1.E0_q+1.E0_q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.E0_q*B2+3.E0_q*B3*RS12+2.E0_q*B4*P1*RSP)
      GGRS = -2.E0_q*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END SUBROUTINE GCOR

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GCOR_SIMD(A, A1, B1, B2, B3, B4, P, RS, GG, GGRS, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS
      TYPE(SIMD_REAL64), INTENT(OUT)     :: GG, GGRS
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(P1, Q0, RS12, RS32, RSP, Q1, Q2, Q3)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         P1 = P + 1.0_q
         Q0 = -2.0_q * A * (1.0_q + A1 * RS%X(II))
         RS12 = SQRT(RS%X(II))
         RS32 = SIMD_POW(RS12, 3.0_q)
         RSP = SIMD_POW(RS%X(II), P)
         Q1 = 2.0_q * A * (B1 * RS12 + B3 * RS32 + (B2 + B4 * RSP) * RS%X(II))
         Q2 = SIMD_LOG(1.0_q + 1.0_q / Q1)
         Q3 = A * (B1 / RS12 + 2.0_q * B2 + 3.0_q * B3 * RS12 + 2.0_q * B4 * P1 * RSP)
         IF (MASK0%X(II)) THEN
            GG%X(II) = Q0 * Q2
            GGRS%X(II) = -2.0_q * A * A1 * Q2 - Q0 * Q3 / (Q1 * Q1 + Q1)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GCOR_SIMD
#endif

!*******************************************************************
!
!>  calculate the exchange correlation potential
!>  on a radial grid and the first 3 derivatives
!
!*******************************************************************

    SUBROUTINE EXCOR_DER_PARA(XC, RHO, NDER, EXCA, TREL)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE

      TYPE (xc_info) XC
      LOGICAL TREL
      INTEGER NDER               ! number of derivatives
      REAL(q) RHO
      REAL(q) EXCA(4),EXCA_(4)
      REAL(q) RHOT,EXC0,EXC2,EXCD0,EXCD2,EPS
      ! this is the best compromise for densities between 1000 and 0.1
      REAL(q),PARAMETER :: DELTA=1E-3_q

!!      PROFILING_START('excor_der_para')

      CALL EXCOR_PARA(XC,RHO,EXCA(1),EXCA(2),TREL)

      IF (NDER>0) THEN
         EPS=DELTA*RHO

         RHOT=RHO-EPS
         CALL EXCOR_PARA(XC,RHOT,EXC0,EXCD0,TREL)

         RHOT=RHO+EPS
         CALL EXCOR_PARA(XC,RHOT,EXC2,EXCD2,TREL)
         ! 1st and 2nd derivative of energy
         EXCA_(2)=(EXC2-EXC0)/ (2*EPS)
         EXCA_(3)=(EXC2+EXC0-2*EXCA(1))/ (EPS*EPS)

         ! 2nd and 3nd derivative of potential=
         ! 1st and 2nd derivative of energy
         EXCA(3)=(EXCD2-EXCD0)/ (2*EPS)
         EXCA(4)=(EXCD2+EXCD0-2*EXCA(2))/ (EPS*EPS)
         ! WRITE(*,'(5E14.7)') EXCA(2),EXCA_(2),EXCA(3),EXCA_(3),EXCA(4)
      ENDIF

!!      PROFILING_STOP('excor_der_para')

    END SUBROUTINE EXCOR_DER_PARA

    SUBROUTINE EXCOR_PARA(XC,RHO,EXC,DEXC,TREL)
      USE constant
      USE setexm_struct_def, ONLY: xc_info
      USE xc_name
      USE xc_family_name

      IMPLICIT NONE
      TYPE (xc_info) XC
      REAL(q) RHO,RS,EXC,DEXC,ALDAC
      REAL(q) ECLDA,ECDLDA,ECD2LDA,SK,T,EC,ECD,ECDD,ZETA,EXC_TMP,DEXC_TMP
      INTEGER IXC
      LOGICAL TREL

!!      PROFILING_START('excor_para')

      ALDAC=XC%ALDAC

      RS = ( (3._q/(4*PI)) / RHO)**(1/3._q) /AUTOA

      EXC  = 0._q
      DEXC = 0._q

      NXCLOOP1: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)/=XC_FAM_LDA).AND.(XC%FAMILY(IXC)/=XC_FAM_GGA)) CYCLE

#ifndef USELIBXC
      IF (XC%ID(IXC)==ID_XC_NOXC) THEN
#else
      IF ((XC%ID(IXC)==ID_XC_NOXC).OR. &
        & (XC%ID(IXC)==ID_XC_LIBXC)) THEN
#endif
         EXC_TMP = 0
         DEXC_TMP= 0
      !
      ! LDA part: Pade approximation to Ceperly Alder results
      !
!aem AM05 functional added, it uses the PW correlation in the LDA part.
! jP: adding PBEsol
!vdw jk
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_PBETEST).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_AM05).OR. &
             & (XC%ID(IXC)==ID_XC_PBESOL).OR. &
             & (XC%ID(IXC)==ID_XC_BEEF).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX)) THEN
!vdw jk
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_X) THEN
         EXC_TMP = EX_MOD(XC,RS,1,.FALSE.)*RHO*RYTOEV
         DEXC_TMP= VX_MOD(XC,RS,1,.FALSE.)*RYTOEV
      ELSE IF (XC%ID(IXC)==ID_XC_PW92) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,TREL)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,TREL)+ECDLDA*ALDAC)*RYTOEV
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE_C).OR. &
             & (XC%ID(IXC)==ID_XC_PW92_C) ) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)
         EXC_TMP = ECLDA*ALDAC*RHO*RYTOEV
         DEXC_TMP= ECDLDA*ALDAC*RYTOEV
! jH-new CORPBE_LDA_RPA
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_RA)  THEN
         ZETA=0
         CALL CORPBE_LDA_RPA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new CORPBE_LDA_SR_0_15au
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_03)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_15au(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV 
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new CORPBE_LDA_SR_0_5A
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_05)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_5A(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new CORPBE_LDA_SR_1_0A
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_10)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_1_0A(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new CORPBE_LDA_SR_2_0A
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_20)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_2_0A(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new CORPBE_LDA_RPA_PLUS
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_PL)  THEN 
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)            
         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
      !  
      ! LDA: range separated exchange + range separated correlation
      !   
      ELSE IF (XC%ID(IXC)==ID_XC_VWN5) THEN
!         EXC_TMP = (EX_MOD(XC,RS,1,TREL)+EC_MOD(XC,RS,1,TREL)*ALDAC)*RHO*RYTOEV
!         DEXC_TMP= (VX_MOD(XC,RS,1,TREL)+VC_MOD(XC,RS,1,TREL)*ALDAC)*RYTOEV      
! Relativistic correction have been removed  
! jP be careful ... 
! jP: what does EC_MOD do in that case
         EXC_TMP = (EX_MOD(XC,RS,1,.FALSE.)+EC_MOD(XC,RS,0._q,1,.FALSE.)*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,.FALSE.)+VC_MOD(XC,RS,0._q,1,.FALSE.)*ALDAC)*RYTOEV
      ELSE IF (XC%ID(IXC)==ID_XC_B3) THEN
! Joachim Paier
         EXC_TMP = (EX_MOD(XC,RS,1,TREL)+ECVOIII(RS,1)*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,TREL)+VCVOIII(RS,1)*ALDAC)*RYTOEV
         ! ECLDA=ECVO(RS,1)
         ! ECDLDA=VCVO(RS,1)
         ! jP: compare with CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)
      ELSE IF (XC%ID(IXC)==ID_XC_B5) THEN
! Joachim Paier
         EXC_TMP = (EX_MOD(XC,RS,1,TREL)+ECVO(RS,1)*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,TREL)+VCVO(RS,1)*ALDAC)*RYTOEV
         ! ECLDA=ECVO(RS,1)
         ! ECDLDA=VCVO(RS,1)
         ! jP: compare with CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)
      ELSE IF (XC%ID(IXC)==ID_XC_COHSM) THEN
         EXC_TMP=0
         CALL COHSM1(RHO, DEXC_TMP)
      !
      ! Exchange LDA part: Perdew Zungers interpolation
      !
      ELSE IF (XC%ID(IXC)==ID_XC_SLATER) THEN
         EXC_TMP = EX_MOD(XC,RS,1,TREL)*RHO*RYTOEV
         DEXC_TMP= VX_MOD(XC,RS,1,TREL)*RYTOEV
      !
      ! Correlation LDA part: Perdew Zungers interpolation
      !
      ELSE IF (XC%ID(IXC)==ID_XC_PZ_C) THEN
         EXC_TMP = ECCA(RS,1)*ALDAC*RHO*RYTOEV
         DEXC_TMP= VCCA(RS,1)*ALDAC*RYTOEV
      !
      ! LDA part: Perdew Zungers interpolation
      !
      ELSE
         EXC_TMP = (EX_MOD(XC,RS,1,TREL)+ECCA(RS,1)*ALDAC)*RHO*RYTOEV
         DEXC_TMP= (VX_MOD(XC,RS,1,TREL)+VCCA(RS,1)*ALDAC)*RYTOEV
      ENDIF

      EXC = EXC + XC%COEFF(IXC)*EXC_TMP
      DEXC = DEXC + XC%COEFF(IXC)*DEXC_TMP

      ENDDO NXCLOOP1

!!      PROFILING_STOP('excor_der_para')

    END SUBROUTINE EXCOR_PARA

!*******************************************************************
!
!>  calculate the exchange correlation energy density
!>  on a radial grid and the first 3 derivatives
!>  for the spinpolarized case
!>  as input the density (up and down) is required
!>  derivatives with respect to up and down components are calculated
!
!*******************************************************************
      
    SUBROUTINE EXCOR_DER(XC,RHOUP,RHODOWN,NDER,EXC,EXCD,EXCDD,EXCDDD,TREL)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE

      TYPE (xc_info) XC
      LOGICAL TREL
      INTEGER NDER               ! number of derivative
      REAL(q) RHOUP,RHODOWN      ! up and down component of density
      REAL(q) RHO(2),RHOIN(2)
      REAL(q) EXC,EXCD(2),EXCD_(2),EXCDD(2,2),EXCDD_(2,2),EXCDDD(2,2,2)
      REAL(q) TMP(-1:1,-1:1),T2(2,-1:1,-1:1)
      REAL(q) EPS(2)
      ! this is the best compromise for densities between 1000 and 0.1
      REAL(q),PARAMETER :: DELTA=1E-3_q
      INTEGER I,J

!!      PROFILING_START('excor_der')

      EXCDD=0
      EXCDDD=0
      RHO(1)=RHOUP
      RHO(2)=RHODOWN
    ! function + derivative
      CALL EXCOR(XC,RHO,EXC,EXCD,TREL)
      IF (NDER>0) THEN
        TMP(0,0) =EXC
        T2(:,0,0)=EXCD
        RHOIN=RHO
        ! calculate steps
        ! exc is approx rho^(4/3)
        ! the derivative thus rho^(1/3) and thus the derivative is
        ! of the order 1/rho
        EPS(1)=DELTA*MAX(RHO(1),RHO(2))
        EPS(2)=EPS(1)
        ! calculate all values on the 3x3 rectangle
        DO I=-1,1
          DO J=-1,1
             RHO(1)=RHOIN(1)+EPS(1)*I
             RHO(2)=RHOIN(2)+EPS(2)*J
             CALL EXCOR(XC,RHO,TMP(I,J),T2(1,I,J),TREL)
          ENDDO
        ENDDO

    ! 1st and 2nd derivative of exchange correlation energy
    ! 1st derivative (EXCD_) is of course equal to EXCD
        EXCD_(1)   =(TMP(1,0)-TMP(-1,0)) / (2*EPS(1))
        EXCD_(2)   =(TMP(0,1)-TMP(0,-1)) / (2*EPS(2))
        EXCDD_(1,1)=(TMP(1,0)+TMP(-1,0)-2*TMP(0,0)) / (EPS(1)*EPS(1))
        EXCDD_(2,2)=(TMP(0,1)+TMP(0,-1)-2*TMP(0,0)) / (EPS(2)*EPS(2))
        EXCDD_(1,2)=((TMP(1,1)+TMP(-1,-1)-2*TMP(0,0))- &
                     (TMP(1,0)+TMP(-1,0)-2*TMP(0,0)) - &
                     (TMP(0,1)+TMP(0,-1)-2*TMP(0,0)))/(2*EPS(1)*EPS(2))
        EXCDD_(2,1)=EXCDD_(1,2)
    ! 1st derivative of potential = 2nd derivative of energy
        EXCDD (1,1)=(T2(1,1,0)-T2(1,-1,0)) / (2*EPS(1))
        EXCDD (2,2)=(T2(2,0,1)-T2(2,0,-1)) / (2*EPS(1))
        EXCDD (1,2)=(T2(1,0,1)-T2(1,0,-1)) / (2*EPS(1))
        EXCDD (2,1)=(T2(2,1,0)-T2(2,-1,0)) / (2*EPS(1))
    ! 2nd derivative of potential = 3rd derivative of energy
        EXCDDD(1,1,1)=(T2(1,1,0)+T2(1,-1,0)-2*T2(1,0,0)) / (EPS(1)*EPS(1))
        EXCDDD(1,2,2)=(T2(1,0,1)+T2(1,0,-1)-2*T2(1,0,0)) / (EPS(2)*EPS(2))
        EXCDDD(1,1,2)=((T2(1,1,1)+T2(1,-1,-1)-2*T2(1,0,0))- &
                       (T2(1,1,0)+T2(1,-1,0)-2*T2(1,0,0))- &
                       (T2(1,0,1)+T2(1,0,-1)-2*T2(1,0,0)))/(2*EPS(1)*EPS(2))
        EXCDDD(1,2,1)=EXCDDD(1,1,2)
        EXCDDD(2,1,1)=(T2(2,1,0)+T2(2,-1,0)-2*T2(2,0,0)) / (EPS(1)*EPS(1))
        EXCDDD(2,2,2)=(T2(2,0,1)+T2(2,0,-1)-2*T2(2,0,0)) / (EPS(2)*EPS(2))
        EXCDDD(2,1,2)=((T2(2,1,1)+T2(2,-1,-1)-2*T2(2,0,0)) - &
                       (T2(2,1,0)+T2(2,-1,0)-2*T2(2,0,0)) - &
                       (T2(2,0,1)+T2(2,0,-1)-2*T2(2,0,0)))/(2*EPS(1)*EPS(2))
        EXCDDD(2,2,1)=EXCDDD(2,1,2)
    !    WRITE(*,*) EXCDDD(2,1,2),EXCDDD(1,2,2),EXCDDD(2,2,2)
    !    WRITE(*,*) EXCDDD(1,1,2),EXCDDD(2,1,1),EXCDDD(1,1,1)
    ! plenty of cross checks are possible here
    !    WRITE(*,'(4E14.7)') EXCD,EXCD_
    !    WRITE(*,'(4E14.7)') EXCDD,EXCDD_
    !    WRITE(*,'(4E14.7)') EXCDDD
      ENDIF

!!      PROFILING_STOP('excor_der')

    END SUBROUTINE EXCOR_DER


    SUBROUTINE EXCOR(XC,RHO,EXC,EXCD,TREL)
      USE prec
      USE constant
      USE setexm_struct_def, ONLY: xc_info
      USE xc_name
      USE xc_family_name

      IMPLICIT NONE
      TYPE (xc_info) XC
      LOGICAL TREL
      REAL(q) RHO(2),EXC,EXCD(2)
    ! local variables
      
      REAL(q) RS,ZETA,FZ,DFZ,EXP,EC,ECP,VXP,VCP,EXT,VX0,DVX,VX1,VX2,RH, &
              ECF,VXF,VCF,EXF,ALP,VALP,ZETA3,ZETA4,EXC_TMP,EXCD_TMP(2),ALDAC
      INTEGER IXC,LFCI

!!      PROFILING_START('excor')

      ALDAC=XC%ALDAC
      LFCI=XC%LFCI

      RH=RHO(1)+RHO(2)
      RS = ( (3._q/(4*PI)) / RH)**(1/3._q) /AUTOA
      ZETA=(RHO(1)-RHO(2))/(RHO(1)+RHO(2))
      ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
      ZETA3=(ZETA*ZETA)*ZETA
      ZETA4=(ZETA*ZETA)*(ZETA*ZETA)

      EXC  = 0._q
      EXCD = 0._q

      NXCLOOP2: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)/=XC_FAM_LDA).AND.(XC%FAMILY(IXC)/=XC_FAM_GGA)) CYCLE

#ifndef USELIBXC
      IF (XC%ID(IXC)==ID_XC_NOXC) THEN
#else
      IF ((XC%ID(IXC)==ID_XC_NOXC).OR. &
        & (XC%ID(IXC)==ID_XC_LIBXC)) THEN
#endif
         EXC_TMP=0
         EXCD_TMP=0
      !
      ! LDA part: Pade approximation to Ceperley Alder results
      !
!aem AM05 functional added, it uses the PW correlation in the LDA part.
! jP: adding PBEsol
!vdw jk
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_PBETEST).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_AM05).OR. &
             & (XC%ID(IXC)==ID_XC_PBESOL).OR. &
             & (XC%ID(IXC)==ID_XC_BEEF).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX)) THEN
!vdw jk
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,.FALSE.)*RYTOEV

         CALL CORPBE_LDA(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ+EC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)+VCF
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_X) THEN
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,.FALSE.)*RYTOEV

         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)
      ELSE IF ((XC%ID(IXC)==ID_XC_PW92)) THEN
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,TREL)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,TREL)*RYTOEV

         CALL CORPBE_LDA(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ+EC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)+VCF
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE_C).OR. &
             & (XC%ID(IXC)==ID_XC_PW92_C)) THEN
         CALL CORPBE_LDA(RS,ZETA,EC,VCP,VCF)
         EXC_TMP=EC*ALDAC*RYTOEV*RH
         EXCD_TMP(1)=VCP*ALDAC*RYTOEV
         EXCD_TMP(2)=VCF*ALDAC*RYTOEV
! jH-new noch fuer RPA und RPA-plus 
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_RA) THEN
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,.FALSE.)*RYTOEV

         CALL CORPBE_LDA_RPA(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ+EC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)+VCF
      ELSE IF (XC%ID(IXC)==ID_XC_ACFDT_PL) THEN
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,.FALSE.)*RYTOEV

         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV

         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ+EC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)+VCF

      ELSE IF (XC%ID(IXC)==ID_XC_VWN5) THEN
! jP used interpolation: PBE
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,.FALSE.)*RYTOEV

         EC= EC_MOD(XC,RS,ZETA,2,TREL)*ALDAC*RYTOEV*RH 
         VCP=VC_MOD(XC,RS,0._q,1,TREL)*ALDAC*RYTOEV 
         VCF=VC_MOD(XC,RS,1._q,2,TREL)*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ+EC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)+VCF

      ELSE IF (XC%ID(IXC)==ID_XC_COHSM) THEN
         EXC_TMP=0
         CALL COHSM1(RH, EXCD_TMP(1))
         EXCD_TMP(2)=EXCD_TMP(1)

      !
      ! LDA part: Perdew Zungers interpolation
      !
      ELSEIF (XC%ID(IXC)==ID_XC_B3) THEN
! Joachim Paier
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,TREL)*RH*RYTOEV
         ECP=ECVOIII(RS,1)*RH*RYTOEV     ; ECF=ECVOIII(RS,2)*RH*RYTOEV
!         ECP=ECVO(RS,1)*RH*RYTOEV     ; ECF=ECVO(RS,2)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,TREL)*RYTOEV
         VCP=VCVOIII(RS,1)*RYTOEV        ; VCF=VCVOIII(RS,2)*RYTOEV
!         VCP=VCVO(RS,1)*RYTOEV        ; VCF=VCVO(RS,2)*RYTOEV

         ALP =ALPHA0_III(RS)*RH*RYTOEV
!         ALP =ALPHA0(RS)*RH*RYTOEV
         VALP=ALPHA1_III(RS)*RYTOEV
!         VALP=ALPHA1(RS)*RYTOEV

         VX0=VXP+(VXF-VXP)*FZ+  (VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
         DVX=((EXF-EXP)*DFZ  +(ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
         EXC_TMP= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)
      ELSEIF (XC%ID(IXC)==ID_XC_B5) THEN
! Joachim Paier
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,TREL)*RH*RYTOEV
         ECP=ECVO(RS,1)*RH*RYTOEV        ; ECF=ECVO(RS,2)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,TREL)*RYTOEV
         VCP=VCVO(RS,1)*RYTOEV           ; VCF=VCVO(RS,2)*RYTOEV

         ALP =ALPHA0(RS)*RH*RYTOEV
         VALP=ALPHA1(RS)*RYTOEV

         VX0=VXP+(VXF-VXP)*FZ+  (VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
         DVX=((EXF-EXP)*DFZ  +(ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
         EXC_TMP= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)

      ELSEIF (XC%ID(IXC)==ID_XC_SLATER) THEN

         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,TREL)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,TREL)*RYTOEV

         VX0=VXP+(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC_TMP= EXP+(EXF-EXP)*FZ

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)

      ELSEIF (XC%ID(IXC)==ID_XC_PZ_C) THEN

         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         ECP=ECCA(RS,1)*RH*RYTOEV        ; ECF=ECCA(RS,2)*RH*RYTOEV
         VCP=VCCA(RS,1)*RYTOEV           ; VCF=VCCA(RS,2)*RYTOEV

         IF (LFCI==1) THEN
            ALP =ALPHA0(RS)*RH*RYTOEV
            VALP=ALPHA1(RS)*RYTOEV
            
            VX0=(VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
            DVX=((ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
                              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
            EXC_TMP= ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC
         ELSE
            VX0= VCP*ALDAC+(VCF-VCP)*FZ*ALDAC
            DVX=(ECF-ECP)*DFZ/RH*ALDAC
            EXC_TMP= ECP*ALDAC+(ECF-ECP)*FZ*ALDAC
         ENDIF

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)

      ELSE

         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(XC,RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(XC,RS,2,TREL)*RH*RYTOEV
         ECP=ECCA(RS,1)*RH*RYTOEV        ; ECF=ECCA(RS,2)*RH*RYTOEV
         VXP=VX_MOD(XC,RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(XC,RS,2,TREL)*RYTOEV
         VCP=VCCA(RS,1)*RYTOEV           ; VCF=VCCA(RS,2)*RYTOEV

         IF (LFCI==1) THEN
            ALP =ALPHA0(RS)*RH*RYTOEV
            VALP=ALPHA1(RS)*RYTOEV
            
            VX0=VXP+(VXF-VXP)*FZ+  (VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
            DVX=((EXF-EXP)*DFZ  +(ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
                              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
            EXC_TMP= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC
         ELSE
            VX0= VXP+(VXF-VXP)*FZ +VCP*ALDAC+(VCF-VCP)*FZ*ALDAC
            DVX=((EXF-EXP)*DFZ    +(ECF-ECP)*DFZ*ALDAC)/RH
            EXC_TMP= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ECF-ECP)*FZ*ALDAC
         ENDIF

         EXCD_TMP(1)= VX0-DVX*(ZETA-1)
         EXCD_TMP(2)= VX0-DVX*(ZETA+1)
      ENDIF

      EXC = EXC + XC%COEFF(IXC)*EXC_TMP
      EXCD(1) = EXCD(1) + XC%COEFF(IXC)*EXCD_TMP(1)
      EXCD(2) = EXCD(2) + XC%COEFF(IXC)*EXCD_TMP(2)

      ENDDO NXCLOOP2

!!      PROFILING_STOP('excor')

    END SUBROUTINE EXCOR

!***********************************************************************
!
!> EX_MOD calculate the exchange energy density per particle eps_x
!> possibly reduced by short range exchange hole or reduce by the amount 
!> that is accounted for by the exact exchange
!> that is 
!> E_x = \int eps_x(rho(r)) rho(r) d^3 r
!>
!> VX_MOD calculates the exchange potential, which is the derivative
!> of eps_x(rho(r)) rho(r) with respect to rho(r)
!> v_x = d (eps_x(rho) rho) /d rho
!> 
!> in leading order these quantities are proportional to the density^(1/3)
!
!***********************************************************************

    FUNCTION EX_MOD(XC,RS,IFLG,TREL)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      REAL(q) EX_MOD, RS, ALDAX, BLDAX, LDASCREEN
      INTEGER IFLG
      LOGICAL TREL, LUSE_LONGRANGE_HF, LUSE_THOMAS_FERMI, LUSE_MODEL_HF

      ALDAX=XC%ALDAX
      BLDAX=XC%BLDAX
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI
      LUSE_MODEL_HF=XC%LUSE_MODEL_HF

      EX_MOD=EX(RS,IFLG,TREL)
      
      IF (LDASCREEN/=0 .OR. ALDAX/=1 ) THEN
         IF (LDASCREEN==0) THEN
            EX_MOD=EX_MOD*ALDAX
         ELSE IF (LUSE_THOMAS_FERMI) THEN
            EX_MOD=EX_MOD-EX_SX(RS,LDASCREEN*AUTOA,IFLG)*(1-ALDAX)
         ELSE IF (LUSE_LONGRANGE_HF) THEN
         ! Iann Gerber: use EX_SR (for ALDAX=0, i.e. AEXX=1)
!            EX_MOD=(1-ALDAX)*EX(RS,IFLG,.FALSE.)-(EX_MOD-EX_SR(RS,LDASCREEN*AUTOA,IFLG))*(1-ALDAX)
            EX_MOD=EX(RS,IFLG,.FALSE.)-(EX_MOD-EX_SR(RS,LDASCREEN*AUTOA,IFLG))*(1-ALDAX)
         ELSE IF (LUSE_MODEL_HF) THEN
         ! model sX 
            EX_MOD=EX(RS,IFLG,.FALSE.)*ALDAX + EX_SR(RS,LDASCREEN*AUTOA,IFLG)*(BLDAX-ALDAX)
         ELSE 
         ! WPBE: ALDAX*EX_SR + EX_LR = EX - (1-ALDAX)*EX_SR
            EX_MOD=EX(RS,IFLG,.FALSE.)-EX_SR(RS,LDASCREEN*AUTOA,IFLG)*(1-ALDAX)
         ENDIF
      ENDIF
    END FUNCTION EX_MOD

    FUNCTION VX_MOD(XC,RS,IFLG,TREL)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      REAL(q) VX_MOD, RS, ALDAX, BLDAX, LDASCREEN
      INTEGER IFLG
      LOGICAL TREL, LUSE_LONGRANGE_HF, LUSE_THOMAS_FERMI, LUSE_MODEL_HF

      ALDAX=XC%ALDAX
      BLDAX=XC%BLDAX
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI
      LUSE_MODEL_HF=XC%LUSE_MODEL_HF

      VX_MOD=VX(RS,IFLG,TREL)
      
      IF (LDASCREEN/=0 .OR. ALDAX/=1 ) THEN
         IF (LDASCREEN==0) THEN
            VX_MOD=VX_MOD*ALDAX
         ELSE IF (LUSE_THOMAS_FERMI) THEN
            VX_MOD=VX_MOD-VX_SX(RS,LDASCREEN*AUTOA,IFLG)*(1-ALDAX)
         ELSE IF (LUSE_LONGRANGE_HF) THEN
         ! Iann Gerber: use VX_SR (for ALDAX=0, i.e. AEXX=1)
            VX_MOD=VX(RS,IFLG,.FALSE.)-(VX_MOD-VX_SR(RS,LDASCREEN*AUTOA,IFLG))*(1-ALDAX)
         ELSE IF (LUSE_MODEL_HF) THEN
         ! model sX 
            VX_MOD=VX(RS,IFLG,.FALSE.)*ALDAX + VX_SR(RS,LDASCREEN*AUTOA,IFLG)*(BLDAX-ALDAX)
         ELSE
         ! WPBE: ALDAX*VX_SR + VX_LR = VX - (1-ALDAX)*VX_SR
            VX_MOD=VX(RS,IFLG,.FALSE.)-VX_SR(RS,LDASCREEN*AUTOA,IFLG)*(1-ALDAX)
         ENDIF
      ENDIF
    END FUNCTION VX_MOD

!***********************************************************************
!
! EC_MOD calculate the correlation energy density per particle eps_c
! that is 
! E_x = \int eps_c(rho(r)) rho(r) d^3 r
!
! VC_MOD calculates the exchange potential, which is the derivative
! of eps_c(rho(r)) rho(r) with respect to rho(r)
! v_c = d (eps_c(rho) rho) /d rho
! 
!
!***********************************************************************
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
! Iann Gerber 18/11/04
! Comments Joachim Paier jP

! Splitting the correlation part: meaning only take into account of the 
! short-range part if necessary
      FUNCTION EC_MOD(XC,RS,ZETA,IFLG,TREL)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      REAL(q) EC_MOD, RS
      REAL(q) ZETA              ! relative spin-polarization zeta := (nu-nd)/(nu+nd)
      REAL(q) EC_MOD_TMP
      REAL(q) LDASCREENC
      INTEGER IFLG
      LOGICAL TREL, LUSE_LONGRANGE_HF, LRANGE_SEPARATED_CORR

      LDASCREENC=XC%LDASCREENC
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LRANGE_SEPARATED_CORR=XC%LRANGE_SEPARATED_CORR

      EC_MOD=ECVO(RS,IFLG)
      IF (LRANGE_SEPARATED_CORR) THEN
         IF (LDASCREENC==0) THEN
            EC_MOD=EC_MOD
! jP here, we ask for the short-range LSDA correlation energy
         ELSE
            CALL ECORRSR(RS,ZETA,LDASCREENC*AUTOA,EC_MOD)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
            EC_MOD=2.0_q*EC_MOD
            IF (LUSE_LONGRANGE_HF) THEN
! keep short range part of LDA correlation energy
!               EC_MOD=ECVO(RS,1)
!               EC_MOD=EC_SR(RS,LDASCREENC*AUTOA,EC_MOD)
               CALL ECORRSR(RS,ZETA,LDASCREENC*AUTOA,EC_MOD)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               EC_MOD=2.0_q*EC_MOD
! keep long range part of the LDA correlation energy
!                  EC_MOD=ECVO(RS,1)
!                  EC_MOD=EC_MOD-EC_SR(RS,LDASCREENC*AUTOA,EC_MOD)  
            ELSE
               CALL ECORRLR(RS,ZETA,LDASCREENC*AUTOA,EC_MOD)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               EC_MOD=2.0_q*EC_MOD            
            ENDIF
         ENDIF
      ENDIF
      END FUNCTION EC_MOD

      FUNCTION VC_MOD(XC,RS,ZETA,IFLG,TREL)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      REAL(q) VC_MOD, RS
      REAL(q)  ZETA, VCSRUP, VCSRDOWN, VCLRUP, VCLRDOWN
      REAL(q) LDASCREENC
      INTEGER IFLG
      LOGICAL TREL, LUSE_LONGRANGE_HF, LRANGE_SEPARATED_CORR

      LDASCREENC=XC%LDASCREENC
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LRANGE_SEPARATED_CORR=XC%LRANGE_SEPARATED_CORR

      VC_MOD=VCVO(RS,IFLG)
      IF (LRANGE_SEPARATED_CORR) THEN
         IF (LDASCREENC==0) THEN
            VC_MOD=VC_MOD
         ELSE
            CALL VCORRSR(RS,ZETA,LDASCREENC*AUTOA,VCSRUP,VCSRDOWN)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
            VCSRUP=2.0_q*VCSRUP
            VCSRDOWN=2.0_q*VCSRDOWN
            IF (LUSE_LONGRANGE_HF) THEN
! keep short range part of LDA correlation potential
!               VC_MOD=VC_SR(RS,LDASCREENC*AUTOA,VCVO(RS,IFLG),ECVO(RS,IFLG))
               CALL VCORRSR(RS,ZETA,LDASCREENC*AUTOA,VCSRUP,VCSRDOWN)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               VCSRUP=2.0_q*VCSRUP
               VCSRDOWN=2.0_q*VCSRDOWN
            ELSE
! keep long range part of the LDA correlation potential
!              VC_MOD=VC_MOD-VC_SR(RS,LDASCREENC*AUTOA,VCVO(RS,IFLG),ECVO(RS,IFLG))
               CALL VCORRLR(RS,ZETA,LDASCREENC*AUTOA,VCLRUP,VCLRDOWN)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               VCSRUP=2.0_q*VCSRUP
               VCSRDOWN=2.0_q*VCSRDOWN
            ENDIF
         ENDIF
      ENDIF
      END FUNCTION VC_MOD
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------

END MODULE LDALIB
