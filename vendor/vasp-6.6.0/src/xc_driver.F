#include  "symbol.inc"
#ifdef USELIBXC
#include "xc_version.h"
#endif
    MODULE xc_driver

    CONTAINS
!************************ SUBROUTINE GGA_DRIVER **********************************
!
!> This subroutine calls the subroutines that calculate xc quantities
!> (non-spin-polarized GGA) for the energy, potential, forces, stress tensor, etc.
!>
!> LLDA allows to include the LDA contribution directly in this
!> routine (not implemented for all subroutines, see below).
!> For Libxc, the LDA contributions are directly included.
!
!*********************************************************************************

      SUBROUTINE GGA_DRIVER(XC,D,DD,EXC,EXCD,EXCDD,LLDA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE setexm_struct_def
      USE xc_name
      USE xc_family_name
      USE ldalib
      USE ggalib
#ifdef USELIBXC
      USE xc_f03_lib_m
#endif

      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(GGA_DRIVER)

      TYPE (xc_info) XC

      INTEGER :: ider, IXC

      REAL(q) THRD, DTHRD
      REAL(q) D, DD, D1, D2, DD1, DD2, DELTA, DDELTA
      REAL(q) FK, RS, S, SIGMA, SK, T, ukfactor
      REAL(q) RHO1, RHO2, DRHO1, DRHO2
      REAL(q) EC, ECD, ECDD, ECLDA, ECDLDA, EXLDA, EXDLDA
      REAL(q) EXDLDA_LR, EXDLDA_SR, EXLDA_LR, EXLDA_SR, EXWPBE_SR, EXWPBED_LR, EXWPBE_LR, EXWPBED_SR, EXWPBEDD_LR, EXWPBEDD_SR
      REAL(q) EXC, EXC1, EXC2, EXCD, EXCD1, EXCD2, EXCDD, EXCDD1, EXCDD2, EXCDDA
      REAL(q) FXC, FXCLDA, FXCD1, FXCD2, FXCDD1, FXCDD2
      REAL(q) VXC, VXCLDA1, VXCLDA2
      REAL(q) EXC_TMP, EXCD_TMP, EXCDD_TMP
      REAL(q) ALDAX,ALDAC,AGGAX,BGGAX,AGGAC,LDASCREEN
#ifdef USELIBXC
      REAL(q) CX,DTMP(1),SIGMATMP(1),EXC1TMP(1),EXCD1TMP(1),EXCDD1TMP(1)
      REAL(q) SIGMAZERO(1),EXC1ZERO(1),EXCD1ZERO(1),EXCDD1ZERO(1)
#endif
#ifdef libbeef
      REAL(q) DDSQ,EXBEEF,EXBEEFD,EXBEEFDD,ECBEEF,ECBEEFD,ECBEEFDD
#endif
      PARAMETER (DDELTA=1E-4_q)
      PARAMETER (THRD=1._q/3._q)
      LOGICAL LLDA,LUSE_LONGRANGE_HF,LUSE_THOMAS_FERMI,LUSE_MODEL_HF
      LOGICAL:: FORCE_PBE=.FALSE.
GPACC !$ACC ROUTINE (GGA91_WB) SEQ
GPACC !$ACC ROUTINE (EXCHPBE) SEQ
GPACC !$ACC ROUTINE (CALC_EXCHWPBE_SP) SEQ
GPACC !$ACC ROUTINE (CORunspPBE) SEQ
GPACC !$ACC ROUTINE (B3LYPXCS) SEQ
GPACC !$ACC ROUTINE (AM05) SEQ
GPACC !$ACC ROUTINE (EXCHPBESOL) SEQ
GPACC !$ACC ROUTINE (CALC_EXCHWPBEsol_SP) SEQ
GPACC !$ACC ROUTINE (CORunspPBESOL) SEQ
GPACC !$ACC ROUTINE (GGAEALL) SEQ
#ifdef USELIBXC
GPACC !$ACC ROUTINE (XC_F03_LDA_VXC) SEQ
GPACC !$ACC ROUTINE (XC_F03_LDA_EXC_VXC) SEQ
GPACC !$ACC ROUTINE (XC_F03_GGA_VXC) SEQ
GPACC !$ACC ROUTINE (XC_F03_GGA_EXC_VXC) SEQ
#endif

      ALDAX=XC%ALDAX
      ALDAC=XC%ALDAC
      AGGAX=XC%AGGAX
      BGGAX=XC%BGGAX
      AGGAC=XC%AGGAC
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI
      LUSE_MODEL_HF=XC%LUSE_MODEL_HF

      EXC_TMP   = 0._q
      EXCD_TMP  = 0._q
      EXCDD_TMP = 0._q

      NXCLOOP1: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)==XC_FAM_GGA).OR.((XC%ID(IXC)==ID_XC_LIBXC).AND.(XC%FAMILY(IXC)==XC_FAM_LDA))) THEN

      EXC   = 0._q
      EXCD  = 0._q
      EXCDD = 0._q

      IF (XC%ID(IXC)==ID_XC_PW91) THEN

        ! PW91 using the routines of Bird and White

        CALL GGA91_WB(XC,D,DD,EXC,EXCD,EXCDD)
        EXC = 2*EXC / D
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
!vdw jk
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX)) THEN
!vdw jk

        ! Perdew Burke Ernzerhof and revised functional

        IF (XC%ID(IXC)==ID_XC_PBE) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF
! Iann Gerber: range separating in GGA exchange
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
          &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
        ELSE
           CALL CALC_EXCHWPBE_SP(XC,D,S, &
          &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
          &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
        ENDIF
! Iann Gerber: end modification

        CALL CORunspPBE(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
        ! Add LDA contributions as well
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX+EXLDA/D*ALDAX+ECLDA  ! eventuell noch mit ALDAC multi
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX+EXDLDA*ALDAX +ECDLDA ! eventuell noch mit ALDAC multi
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 ! gk replaced ECLD by ECLDA, 18.09.2018
                 EXC  =EC  *AGGAC+ EXWPBE_LR/D*AGGAX +EXWPBE_SR/D +ECLDA! ist EXWPBE nur diff
                 EXCD =ECD *AGGAC+ EXWPBED_LR*AGGAX+EXWPBED_SR+ECDLDA   !
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR !
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC +EXWPBE_LR/D*AGGAX + EXWPBE_SR/D*BGGAX + ECLDA
                 EXCD =ECD *AGGAC +EXWPBED_LR*AGGAX  + EXWPBED_SR*BGGAX  + ECDLDA
                 EXCDD=ECDD*AGGAC +EXWPBEDD_LR*AGGAX + EXWPBEDD_SR*BGGAX                
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+ EXWPBE_SR/D*AGGAX +EXWPBE_LR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_SR*AGGAX+EXWPBED_LR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ELSE
        ! Do not add LDA contributions
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX  +(EXWPBE_SR-EXLDA_SR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX+(EXWPBED_SR-EXDLDA_SR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX + (EXWPBE_SR-EXLDA_SR)/D*BGGAX
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX + (EXWPBED_SR-EXDLDA_SR)*BGGAX
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX + EXWPBEDD_SR*BGGAX                 
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR-EXLDA_SR)/D*AGGAX  +(EXWPBE_LR-EXLDA_LR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_SR-EXDLDA_SR)*AGGAX+(EXWPBED_LR-EXDLDA_LR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2._q*EXC
        EXCD =2._q*EXCD
        EXCDD=2._q*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_X) THEN

        ! exchange-only PBE-type functionals

        IF (XC%ID(IXC)==ID_XC_PBE_X) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
        ELSE
           S=0.0_q
        ENDIF
        CALL EXCHPBE(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD,ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

        IF (LLDA) THEN
           EXC  =(EXC-EXLDA)/D*AGGAX+EXLDA/D*ALDAX
           EXCD =(EXCD-EXDLDA)*AGGAX+EXDLDA*ALDAX
           EXCDD=EXCDD*AGGAX
        ELSE
           EXC  =(EXC-EXLDA)/D*AGGAX
           EXCD =(EXCD-EXDLDA)*AGGAX
           EXCDD=EXCDD*AGGAX
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2._q*EXC
        EXCD =2._q*EXCD
        EXCDD=2._q*EXCDD

      ELSE IF (XC%ID(IXC)==ID_XC_PBE_C) THEN

        ! correlation-only PBE-type functionals

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           T=DD/(D*SK*2._q)
        ELSE
           T=0.0_q
        ENDIF

        CALL CORunspPBE(RS,ECLDA,ECDLDA,SK,T,EC,ECD,ECDD,.TRUE.)

        IF (LLDA) THEN
           EXC  =EC  *AGGAC + ECLDA
           EXCD =ECD *AGGAC + ECDLDA
           EXCDD=ECDD*AGGAC
        ELSE
           EXC  =EC  *AGGAC
           EXCD =ECD *AGGAC
           EXCDD=ECDD*AGGAC
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2._q*EXC
        EXCD =2._q*EXCD
        EXCDD=2._q*EXCDD

      ELSE IF ((XC%ID(IXC)==ID_XC_B3).OR.(XC%ID(IXC)==ID_XC_B5)) THEN

        RHO1=D/2.0_q
        RHO2=D/2.0_q
        DRHO1=DD/2.0_q
        DRHO2=DD/2.0_q

        CALL  B3LYPXCS(XC,RHO1,RHO2,DRHO1,DRHO2,DD,EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA)


        EXC     = 2._q* EXC/D
        EXCD1   = 2._q* EXCD1
        EXCD2   = 2._q* EXCD2
        EXCDD1  = 2._q* EXCDD1
        EXCDD2  = 2._q* EXCDD2
        EXCDDA  = 2._q* EXCDDA

        EXCD    = EXCD1
        EXCDD   = EXCDD1+EXCDDA

!       WRITE(91,'(7F14.7)') RHO1, RHO2, DRHO1, DRHO2, EXC, EXCD, EXCDD

! jP
!aem AM05 added

      ELSE IF (XC%ID(IXC)==ID_XC_AM05) THEN

        !aem DD sometimes comes in negative.
        !aem Since AM05 assumes that DD actually IS |grad rho|,
        !aem and thus DD always should be positive,
        !aem we need to use ABS(DD) as input to AM05.
        !aem However, my routine needs to give the same symmetries out as
        !aem PBE and PW91, therefor the sign-compensation of EXCDD.

        D1=D/2._q
        D2=D/2._q
        DD1=DD/2._q
        DD2=DD/2._q
        ider = 1
        CALL AM05(XC,D1,D2,ABS(DD1),ABS(DD2),ider, &
             FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
             FXCDD1,FXCDD2)

        IF (LLDA) THEN
          EXC = FXC/D
          EXCD = FXCD1
        ELSE
          EXC = (FXC - FXCLDA)/D
          EXCD = FXCD1 - VXCLDA1
        ENDIF
        EXCDD = SIGN(1.0_q,DD)*FXCDD1

        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

      ELSE IF (XC%ID(IXC)==ID_XC_PBESOL) THEN

        ! Perdew Burke Ernzerhof and revised functional

           ukfactor=1.0_q

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF
! Iann Gerber: range separating in GGA exchange
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(XC,D,S, &
          &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
          &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
        ENDIF
! Iann Gerber: end modification

        CALL CORunspPBESOL(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
        ! Add LDA contributions as well
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX+EXLDA/D*ALDAX+ECLDA
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX+EXDLDA*ALDAX +ECDLDA
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+ EXWPBE_LR/D*AGGAX +EXWPBE_SR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_LR*AGGAX+EXWPBED_SR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC +EXWPBE_LR/D*AGGAX + EXWPBE_SR/D*BGGAX + ECLDA
                 EXCD =ECD *AGGAC +EXWPBED_LR*AGGAX  + EXWPBED_SR*BGGAX  + ECDLDA
                 EXCDD=ECDD*AGGAC +EXWPBEDD_LR*AGGAX + EXWPBEDD_SR*BGGAX                                  
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+ EXWPBE_SR/D*AGGAX +EXWPBE_LR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_SR*AGGAX+EXWPBED_LR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ELSE
        ! Do not add LDA contributions
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX  +(EXWPBE_SR-EXLDA_SR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX+(EXWPBED_SR-EXDLDA_SR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX + (EXWPBE_SR-EXLDA_SR)/D*BGGAX
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX + (EXWPBED_SR-EXDLDA_SR)*BGGAX
                 EXCDD=ECDD*AGGAC+EXWPBEDD_LR*AGGAX + EXWPBEDD_SR*BGGAX
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR-EXLDA_SR)/D*AGGAX  +(EXWPBE_LR-EXLDA_LR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_SR-EXDLDA_SR)*AGGAX+(EXWPBED_LR-EXDLDA_LR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
! jP: end of adding PBEsol

! BEEF
#ifdef libbeef
      ELSE IF (XC%ID(IXC)==ID_XC_BEEF) THEN
        IF (D>1E-10_q) THEN
           DDSQ = DD*DD
           CALL BEEFX(D, DDSQ, EXBEEF, EXBEEFD, EXBEEFDD, LLDA)
           CALL BEEFLOCALCORR(D, DDSQ, ECBEEF, ECBEEFD, ECBEEFDD, LLDA)
           EXC = 2._q*(AGGAX*EXBEEF + AGGAC*ECBEEF)/D
           EXCD = 2._q*(AGGAX*EXBEEFD + AGGAC*ECBEEFD)
           EXCDD = 2._q*(AGGAX*EXBEEFDD + AGGAC*ECBEEFDD)*DD
        ELSE
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
        ENDIF
#endif

      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN
! functionals for range-separated ACFDT (LDA - short range RPA)
! dummy stub since these functionals are not gradient corrected
         EXC=0._q ; EXCD=0._q ; EXCDD=0._q

      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_RA).OR.(XC%ID(IXC)==ID_XC_ACFDT_PL)) THEN
         EXC=0._q ; EXCD=0._q ; EXCDD=0._q

#ifdef USELIBXC
      ELSE IF (XC%ID(IXC)==ID_XC_LIBXC) THEN
         DTMP(1)=D
         SIGMATMP(1)=DD*DD
         SIGMAZERO(1)=0._q
         CX = -(3._q/4._q)*(3._q/PI)**THRD
         EXLDA = CX*D**THRD
         EXDLDA = EXLDA*4._q/3._q

         IF ((LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_LDA) .OR. (LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_HYB_LDA)) THEN
            IF (LIBXC_NOEXC(IXC)) THEN
               CALL XC_F03_LDA_VXC(LIBXC_FUNC_NSP(IXC),1_qi8,DTMP,EXCD1TMP)
               EXC1TMP = 0._q
            ELSE
               CALL XC_F03_LDA_EXC_VXC(LIBXC_FUNC_NSP(IXC),1_qi8,DTMP,EXC1TMP,EXCD1TMP)
            ENDIF
            IF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE) .AND. (.NOT.LIBXC_HYB(IXC))) THEN
               EXC  = ALDAX*EXC1TMP(1)
               EXCD = ALDAX*EXCD1TMP(1)
            ELSEIF (LIBXC_KIND(IXC) .EQ. XC_CORRELATION) THEN
               EXC  = ALDAC*EXC1TMP(1)
               EXCD = ALDAC*EXCD1TMP(1)
            ELSEIF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE_CORRELATION) .OR. LIBXC_HYB(IXC)) THEN
               EXC  = EXC1TMP(1)
               EXCD = EXCD1TMP(1)
            ENDIF
         ELSEIF ((LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_GGA) .OR. (LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_HYB_GGA)) THEN
            IF (LIBXC_NOEXC(IXC)) THEN
               CALL XC_F03_GGA_VXC(LIBXC_FUNC_NSP(IXC),1_qi8,DTMP,SIGMATMP,EXCD1TMP,EXCDD1TMP)
               EXC1TMP = 0._q
            ELSE
               CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_NSP(IXC),1_qi8,DTMP,SIGMATMP,EXC1TMP,EXCD1TMP,EXCDD1TMP)
            ENDIF
            IF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE) .AND. (.NOT.LIBXC_HYB(IXC))) THEN
               EXC   = AGGAX*(EXC1TMP(1)-EXLDA) + ALDAX*EXLDA
               EXCD  = AGGAX*(EXCD1TMP(1)-EXDLDA) + ALDAX*EXDLDA
               EXCDD = AGGAX*EXCDD1TMP(1)*2._q*DD
            ELSEIF (LIBXC_KIND(IXC) .EQ. XC_CORRELATION) THEN
               IF (LIBXC_NOEXC(IXC)) THEN
                  CALL XC_F03_GGA_VXC(LIBXC_FUNC_NSP(IXC),1_qi8,DTMP,SIGMAZERO,EXCD1ZERO,EXCDD1ZERO)
                  EXC1ZERO = 0._q
               ELSE
                  CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_NSP(IXC),1_qi8,DTMP,SIGMAZERO,EXC1ZERO,EXCD1ZERO,EXCDD1ZERO)
               ENDIF
               EXC   = AGGAC*(EXC1TMP(1)-EXC1ZERO(1)) + ALDAC*EXC1ZERO(1)
               EXCD  = AGGAC*(EXCD1TMP(1)-EXCD1ZERO(1)) + ALDAC*EXCD1ZERO(1)
               EXCDD = AGGAC*EXCDD1TMP(1)*2._q*DD
            ELSEIF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE_CORRELATION) .OR. LIBXC_HYB(IXC)) THEN
               EXC   = EXC1TMP(1)
               EXCD  = EXCD1TMP(1)
               EXCDD = EXCDD1TMP(1)*2._q*DD
            ENDIF
         ENDIF
         !Conversion to Rydbergs
         EXC   = 2._q*EXC
         EXCD  = 2._q*EXCD
         EXCDD = 2._q*EXCDD
#endif
      ELSE

        ! for all other functionals
        ! we use finite differences to calculate the required
        ! quantities
        ! presently no other functional are supported
        ! but in case one needs these routines

        DELTA=MIN(DDELTA,ABS(D)/100)
        D1=D-DELTA
        D2=D+DELTA
        CALL GGAEALL(D1,DD,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(D2,DD,1._q,1._q,VXC,EXC2)
        EXCD=(EXC2*D2-EXC1*D1)/MAX((D2-D1),1E-10_q)

        DELTA=MIN(DDELTA,ABS(DD)/100)
        DD1=DD-DELTA
        DD2=DD+DELTA
        CALL GGAEALL(D,DD1,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(D,DD2,1._q,1._q,VXC,EXC2)
        EXCDD=(EXC2*D-EXC1*D)/MAX((DD2-DD1),1E-10_q)
        CALL GGAEALL(D,DD,1._q,1._q,VXC,EXC)

!aem added in case somebody (like me) might want to use this routine
        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!aem end addition

      ENDIF

      EXC_TMP   = EXC_TMP   + XC%COEFF(IXC)*EXC
      EXCD_TMP  = EXCD_TMP  + XC%COEFF(IXC)*EXCD
      EXCDD_TMP = EXCDD_TMP + XC%COEFF(IXC)*EXCDD

      ENDIF

      ENDDO NXCLOOP1

      EXC   = EXC_TMP
      EXCD  = EXCD_TMP
      EXCDD = EXCDD_TMP

      RETURN
      END SUBROUTINE GGA_DRIVER

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GGA_DRIVER_SIMD(XC, D_IN, DD, EXC, EXCD, EXCDD, LLDA, MASK0)
      USE prec
      USE constant
      USE setexm_struct_def
      USE xc_name
      USE xc_family_name
      USE ldalib
      USE ggalib
      USE simd
      USE tutor, ONLY: vtutor
      IMPLICIT REAL(q) (A-H, O-Z)
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: D_IN, DD
      TYPE(SIMD_REAL64), INTENT(INOUT)   :: EXC, EXCD, EXCDD
      LOGICAL, INTENT(IN)                :: LLDA
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: D, RS, SK, S, T
      TYPE(SIMD_REAL64)                  :: EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR
      TYPE(SIMD_REAL64)                  :: EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR
      TYPE(SIMD_REAL64)                  :: ECLDA, ECDLDA, EC, ECD, ECDD
      TYPE(SIMD_REAL64)                  :: DTHRD, UKFACTOR
      TYPE(SIMD_REAL64)                  :: RHO1, RHO2, DRHO1, DRHO2
      TYPE(SIMD_REAL64)                  :: EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA
      TYPE(SIMD_REAL64)                  :: D1, D2, DD1, DD2, FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2
      TYPE(SIMD_INT64)                   :: IDER
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II,IXC
      LOGICAL                            :: TRUE_FOR_ANY
      LOGICAL, SAVE                      :: INIT_PRED = .TRUE., PRED_1, PRED_2, PRED_3, PRED_4, PRED_5, PRED_6, PRED_7
      LOGICAL                            :: LUSE_LONGRANGE_HF,LUSE_THOMAS_FERMI
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5, TEMP_6, TEMP_7, TEMP_8, TEMP_9
      REAL(q)                            :: ALDAX,ALDAC,AGGAX,AGGAC,LDASCREEN
      REAL(q), PARAMETER                 :: DDELTA = 1.0E-4_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q ! 1._q/3._q
      REAL(q), PARAMETER                 :: P1 = 0.620350490899400086597382_q   ! (3._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P2 = 3.093667726280135532590521_q   ! (3._q*pi*pi)**(1._q/3._q)
      REAL(q), PARAMETER                 :: P3 = 0.161620459673995492400067_q   ! 1._q/((3._q*pi*pi)**(1._q/3._q))
      REAL(q), PARAMETER                 :: P4 = 1.273239544735162764865777_q   ! (4._q/pi)
      REAL(q), PARAMETER                 :: P5 = 0.251928970342244873403814_q   ! 1._q/sqrt(p2 * p4)
      LOGICAL, PARAMETER                 :: FORCE_PBE = .FALSE.

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            D%X(II) = D_IN%X(II)
         ELSE
            D%X(II) = 1.0E-20_q
         ENDIF
      ENDDO

      ALDAX=XC%ALDAX
      ALDAC=XC%ALDAC
      AGGAX=XC%AGGAX
      AGGAC=XC%AGGAC
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI

      NXCLOOP2: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)==XC_FAM_GGA).OR.((XC%ID(IXC)==ID_XC_LIBXC).AND.(XC%FAMILY(IXC)==XC_FAM_LDA))) THEN

      IF (XC%ID(IXC)==ID_XC_PW91) THEN
         CALL GGA91_WB_SIMD(XC, D, DD, EXC, EXCD, EXCDD, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            IF (MASK0%X(II)) THEN
               EXC%X(II) = 2.0_q * EXC%X(II) / D%X(II)
               EXCD%X(II) = 2.0_q * EXCD%X(II)
               EXCDD%X(II) = 2.0_q * EXCDD%X(II)
            ENDIF
         ENDDO
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PBESOL).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX)) THEN
! Perdew Burke Ernzerhof and revised functional
         IF (INIT_PRED) THEN
            PRED_1 = (LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE
            PRED_2 = LLDA .AND. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)
            PRED_3 = LLDA .AND. (.NOT. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. LUSE_LONGRANGE_HF
            PRED_4 = LLDA .AND. (.NOT. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. (.NOT. LUSE_LONGRANGE_HF)
            PRED_5 = (.NOT. LLDA) .AND. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)
            PRED_6 = (.NOT. LLDA) .AND. (.NOT.((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. LUSE_LONGRANGE_HF
            PRED_7 = (.NOT. LLDA) .AND. (.NOT.((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. (.NOT. LUSE_LONGRANGE_HF)
!The following line must be commented, otherwise the results for PBE-HF-based range-separated methods are wrong.
!            INIT_PRED = .FALSE.
         ENDIF

         TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            MASK1%X(II) = MASK0%X(II)
            IF (MASK1%X(II)) THEN
               IF (D%X(II) .GT. 0.0_q) THEN
                  TRUE_FOR_ANY = .TRUE.
               ELSE
                  MASK1%X(II) = .FALSE.
                  D%X(II) = 1.0E-20_q
                  EXC%X(II) = 0.0_q
                  EXCD%X(II) = 0.0_q
                  EXCDD%X(II) = 0.0_q
               ENDIF
            ENDIF
         ENDDO

! return of no SIMD lane continues execution
         IF (.NOT. TRUE_FOR_ANY) RETURN

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            ! Perdew Burke Ernzerhof and revised functional
            IF ((XC%ID(IXC)==ID_XC_PBE).OR.(XC%ID(IXC)==ID_XC_PBESOL)) THEN
               UKFACTOR%X(II) = 1.0_q
            ELSE
               UKFACTOR%X(II) = 0.0_q
            ENDIF
            DTHRD%X(II) = SIMD_POW(D%X(II), THRD)
            TEMP_1 = 1.0_q / DTHRD%X(II)
            RS%X(II) = P1 * TEMP_1
            SK%X(II) = SQRT(P2 * P4 * DTHRD%X(II))
            TEMP_2 = 1.0_q / D%X(II)
            IF(D%X(II) .GT. 1.0E-10_q)THEN
               S%X(II) = DD%X(II) * P3 * TEMP_1 * TEMP_2
               T%X(II) = DD%X(II) * P5 * SQRT(TEMP_1) * TEMP_2
            ELSE
               S%X(II) = 0.0_q
               T%X(II) = 0.0_q
            ENDIF
         ENDDO

! Iann Gerber: range separating in GGA exchange
         IF (PRED_1) THEN
            IF (XC%ID(IXC)==ID_XC_PBESOL) THEN
               CALL EXCHPBESOL_SIMD(D, DTHRD, S, EXLDA, EXC, EXDLDA, EXCD, EXCDD, UKFACTOR, MASK1)
            ELSE
               CALL EXCHPBE_SIMD(D, DTHRD, S, EXLDA, EXC, EXDLDA, EXCD, EXCDD, UKFACTOR, MASK1, XC%ID(IXC))
            ENDIF
         ELSE
            CALL CALC_EXCHWPBE_SP_SIMD(XC, D, S, &
            EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR, &
            EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR, MASK1, XC%ID(IXC))
         ENDIF
! Iann Gerber: end modification
         CALL CORUNSPPBE_SIMD(RS, ECLDA, ECDLDA, SK, T, EC, ECD, ECDD, .TRUE., MASK1, XC%ID(IXC))

         IF (PRED_2) THEN
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + (EXC%X(II) - EXLDA%X(II)) * TEMP_1 * AGGAX + EXLDA%X(II) * TEMP_1 * ALDAX + ECLDA%X(II) ! eventuell noch mit ALDAC multi
               TEMP_3 = ECD%X(II) * AGGAC + (EXCD%X(II) - EXDLDA%X(II)) * AGGAX + EXDLDA%X(II) * ALDAX + ECDLDA%X(II) ! eventuell noch mit ALDAC multi
               TEMP_4 = ECDD%X(II) * AGGAC + EXCDD%X(II) * AGGAX
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_3) THEN
         ! Iann Gerber's functionals
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + EXWPBE_LR%X(II) *TEMP_1 * AGGAX + EXWPBE_SR%X(II) * TEMP_1 + ECLD ! ist EXWPBE nur diff
               TEMP_3 = ECD%X(II) * AGGAC + EXWPBED_LR%X(II) * AGGAX + EXWPBED_SR%X(II) + ECDLDA%X(II)
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_LR%X(II) * AGGAX + EXWPBEDD_SR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_4) THEN
         ! wPBE
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + EXWPBE_SR%X(II) * TEMP_1 * AGGAX + EXWPBE_LR%X(II) * TEMP_1 + ECLDA%X(II)
               TEMP_3 = ECD%X(II) * AGGAC + EXWPBED_SR%X(II) * AGGAX + EXWPBED_LR%X(II) + ECDLDA%X(II)
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_SR%X(II) * AGGAX + EXWPBEDD_LR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_5) THEN
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_2 = EC%X(II) * AGGAC + (EXC%X(II) - EXLDA%X(II)) / D%X(II) * AGGAX
               TEMP_3 = ECD%X(II) * AGGAC + (EXCD%X(II) - EXDLDA%X(II)) * AGGAX
               TEMP_4 = ECDD%X(II) * AGGAC + EXCDD%X(II) * AGGAX
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_6) THEN
         ! Iann Gerber's functionals
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + (EXWPBE_LR%X(II) - EXLDA_LR%X(II)) * TEMP_1 * AGGAX + (EXWPBE_SR%X(II) - EXLDA_SR%X(II)) * TEMP_1
               TEMP_3 = ECD%X(II) * AGGAC + (EXWPBED_LR%X(II) - EXDLDA_LR%X(II)) * AGGAX + (EXWPBED_SR%X(II) - EXDLDA_SR%X(II))
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_LR%X(II) * AGGAX + EXWPBEDD_SR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_7) THEN
         ! wPBE
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + (EXWPBE_SR%X(II) - EXLDA_SR%X(II)) * TEMP_1 * AGGAX + (EXWPBE_LR%X(II) - EXLDA_LR%X(II)) * TEMP_1
               TEMP_3 = ECD%X(II) * AGGAC + (EXWPBED_SR%X(II) - EXDLDA_SR%X(II)) * AGGAX + (EXWPBED_LR%X(II) - EXDLDA_LR%X(II))
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_SR%X(II) * AGGAX + EXWPBEDD_LR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF
      ELSE IF ((XC%ID(IXC)==ID_XC_B3).OR.(XC%ID(IXC)==ID_XC_B5)) THEN
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            RHO1%X(II) = 0.5_q * D%X(II)
            RHO2%X(II) = 0.5_q * D%X(II)
            DRHO1%X(II) = 0.5_q * DD%X(II)
            DRHO2%X(II) = 0.5_q * DD%X(II)
         ENDDO

         CALL B3LYPXCS_SIMD(XC, RHO1, RHO2, DRHO1, DRHO2, DD, EXC, EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1        
            EXCD1%X(II) = 2.0_q * EXCD1%X(II)
            EXCD2%X(II) = 2.0_q * EXCD2%X(II)
            EXCDD1%X(II) = 2.0_q * EXCDD1%X(II)
            EXCDD2%X(II) = 2.0_q * EXCDD2%X(II)
            EXCDDA%X(II) = 2.0_q * EXCDDA%X(II)
            IF (MASK0%X(II)) THEN
               EXC%X(II) = 2.0_q * EXC%X(II) / D%X(II)
               EXCD%X(II) = EXCD1%X(II)
               EXCDD%X(II) = EXCDD1%X(II) + EXCDDA%X(II)
            ENDIF
         ENDDO
      ELSE IF (XC%ID(IXC)==ID_XC_AM05) THEN
         !aem DD sometimes comes in negative.
         !aem Since AM05 assumes that DD actually IS |grad rho|,
         !aem and thus DD always should be positive,
         !aem we need to use ABS(DD) as input to AM05.
         !aem However, my routine needs to give the same symmetries out as
         !aem PBE and PW91, therefor the sign-compensation of EXCDD.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            D1%X(II) = 0.5_q * D%X(II)
            D2%X(II) = D1%X(II)
            DD1%X(II) = ABS(0.5_q * DD%X(II))
            DD2%X(II) = DD1%X(II)
            IDER%X(II) = 1
         ENDDO
         
         CALL AM05_SIMD(XC, D1, D2, DD1, DD2, IDER, FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5, TEMP_6, TEMP_7, TEMP_8, TEMP_9, D_INV)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            TEMP_1 = EXC%X(II)
            TEMP_2 = EXCD%X(II)
            TEMP_3 = EXCDD%X(II)
            TEMP_4 = FXC%X(II)
            TEMP_5 = FXCD1%X(II)
            TEMP_6 = FXCLDA%X(II)
            TEMP_7 = VXCLDA1%X(II)
            TEMP_8 = DD%X(II)
            TEMP_9 = FXCDD1%X(II)
            D_INV = 1.0_q / D%X(II)
            IF (MASK0%X(II)) THEN
               TEMP_1 = TEMP_4 * D_INV
               TEMP_2 = TEMP_5
               IF (.NOT. LLDA) THEN
                  TEMP_1 = TEMP_1 - TEMP_6 * D_INV
                  TEMP_2 = TEMP_2 - TEMP_7
               ENDIF
               TEMP_3 = TEMP_9
               IF (TEMP_8 .LT. 0.0_q) THEN
                  TEMP_3 = -TEMP_3
               ENDIF
               ! Hartree -> Rydberg conversion
               TEMP_1 = 2.0_q * TEMP_1
               TEMP_2 = 2.0_q * TEMP_2
               TEMP_3 = 2.0_q * TEMP_3
            ENDIF
            EXC%X(II) = TEMP_1
            EXCD%X(II) = TEMP_2
            EXCDD%X(II) = TEMP_3
         ENDDO
      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN
         IF (MASK0%X(II)) THEN
            EXC%X(II) = 0.0_q
            EXCD%X(II) = 0.0_q
            EXCDD%X(II) = 0.0_q
         ENDIF
      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_RA).OR.(XC%ID(IXC)==ID_XC_ACFDT_PL)) THEN
         IF (MASK0%X(II)) THEN
            EXC%X(II) = 0.0_q
            EXCD%X(II) = 0.0_q
            EXCDD%X(II) = 0.0_q
         ENDIF
      ELSE
         CALL vtutor%bug("internal ERROR GGAEALL: Wrong functional, scheme not implemented!", __FILE__, __LINE__)
      ENDIF

      ENDIF

      ENDDO NXCLOOP2

      RETURN
      END SUBROUTINE GGA_DRIVER_SIMD
#endif


!************************ SUBROUTINE GGASPIN_DRIVER ****************************
!
!> This subroutine calls the subroutines that calculate xc quantities
!> (spin-polarized GGA) for the energy, potential, forces, stress tensor, etc.
!>
!> LLDA allows to include the LDA contribution directly in this
!> routine (not implemented for all subroutines, see below).
!> For Libxc, the LDA contributions are directly included.
!
!*******************************************************************************

      SUBROUTINE GGASPIN_DRIVER(XC,D1,D2,DD1,DD2,DDA,EXC, &
         excd1,excd2,excq1,excq2,ecq,LLDA)
GPACC !$ACC ROUTINE SEQ

!     D1   density up
!     D2   density down
!     DD1  |gradient of density up|
!     DD2  |gradient of density down|
!     DDA  |gradient of the total density|
!     LLDA add lda contributions

      USE prec
      USE constant
      USE setexm_struct_def
      USE xc_name
      USE xc_family_name
      USE tutor, ONLY : vtutor
      USE ggalib
#ifdef USELIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(GGASPIN_DRIVER)

      TYPE (xc_info) XC

      INTEGER :: ider, IXC

      REAL(q) THRD, DTHRD, THRD4
      REAL(q) D, D1, D2, DD1, DD2, DDA
      REAL(q) G, FK, RS, S, SK, T, ukfactor, ZETA
      REAL(q) ECQ, EC, ECD1, ECD1LDA, ECD2, ECD2LDA, ECLDA, EXDLDA1, EXDLDA2, EXLDA1, EXLDA2
      REAL(q) EXDLDA_LR1, EXDLDA_LR2, EXDLDA_SR1, EXDLDA_SR2, EXLDA_LR1, EXLDA_LR2, EXLDA_SR1, EXLDA_SR2
      REAL(q) EXWPBE_LR1, EXWPBE_LR2, EXWPBE_SR1, EXWPBE_SR2, EXWPBED_LR1, EXWPBED_LR2
      REAL(q) EXWPBED_SR1, EXWPBED_SR2, EXWPBEDD_LR1, EXWPBEDD_LR2, EXWPBEDD_SR1, EXWPBEDD_SR2
      REAL(q) EXC, EXCD1, EXCD2, EXCQ1, EXCQ2, EXC1, EXC2, EXCDD1, EXCDD2, EXCDDA
      REAL(q) FXC, FXCD1, FXCD2, FXCDD1, FXCDD2, FXCLDA
      REAL(q) VXCLDA1, VXCLDA2
      REAL(q) EXC_TMP, EXCD1_TMP, EXCD2_TMP, EXCQ1_TMP, EXCQ2_TMP, ECQ_TMP
      REAL(q) ALDAX,ALDAC,AGGAX,BGGAX,AGGAC,LDASCREEN
#ifdef USELIBXC
      REAL(q) CX,EXLDA,XC_DENS(2),XC_SIGMA(3),XC_EXC(1),XC_VRHO(2),XC_VSIGMA(3)
      REAL(q) XC_SIGMAZERO(3),XC_EXCZERO(1),XC_VRHOZERO(2),XC_VSIGMAZERO(3)
#endif
#ifdef libbeef
      REAL(q) DDSQ,EXBEEF1,EXBEEF2,EXBEEFD1,EXBEEFD2,ECBEEF,ECBEEFD1,ECBEEFD2
#endif
      LOGICAL LLDA
      LOGICAL LUSE_LONGRANGE_HF,LUSE_THOMAS_FERMI,LUSE_MODEL_HF
      PARAMETER (THRD=1._q/3._q,THRD4=4._q/3._q)
      LOGICAL:: FORCE_PBE=.FALSE.
GPACC !$ACC ROUTINE (GGAXCS) SEQ
GPACC !$ACC ROUTINE (EXCHPBE) SEQ
GPACC !$ACC ROUTINE (CALC_EXCHWPBE_SP) SEQ
GPACC !$ACC ROUTINE (CORPBE) SEQ
GPACC !$ACC ROUTINE (B3LYPXCS) SEQ
GPACC !$ACC ROUTINE (AM05) SEQ
GPACC !$ACC ROUTINE (EXCHPBESOL) SEQ
GPACC !$ACC ROUTINE (EXCHPBESOL) SEQ
GPACC !$ACC ROUTINE (CALC_EXCHWPBEsol_SP) SEQ
GPACC !$ACC ROUTINE (corpbesol) SEQ
#ifdef USELIBXC
GPACC !$ACC ROUTINE (XC_F03_LDA_VXC) SEQ
GPACC !$ACC ROUTINE (XC_F03_LDA_EXC_VXC) SEQ
GPACC !$ACC ROUTINE (XC_F03_GGA_VXC) SEQ
GPACC !$ACC ROUTINE (XC_F03_GGA_EXC_VXC) SEQ
#endif

      ALDAX=XC%ALDAX
      ALDAC=XC%ALDAC
      AGGAX=XC%AGGAX
      BGGAX=XC%BGGAX
      AGGAC=XC%AGGAC
      LDASCREEN=XC%LDASCREEN
      LUSE_LONGRANGE_HF=XC%LUSE_LONGRANGE_HF
      LUSE_THOMAS_FERMI=XC%LUSE_THOMAS_FERMI
      LUSE_MODEL_HF=XC%LUSE_MODEL_HF

      EXC_TMP   = 0._q
      EXCD1_TMP = 0._q
      EXCD2_TMP = 0._q
      EXCQ1_TMP = 0._q
      EXCQ2_TMP = 0._q
      ECQ_TMP   = 0._q

      NXCLOOP4: DO IXC=1,XC%NXC

      IF ((XC%FAMILY(IXC)==XC_FAM_GGA).OR.((XC%ID(IXC)==ID_XC_LIBXC).AND.(XC%FAMILY(IXC)==XC_FAM_LDA))) THEN

      EXC   = 0._q
      EXCD1 = 0._q
      EXCD2 = 0._q
      EXCQ1 = 0._q
      EXCQ2 = 0._q
      ECQ   = 0._q

      IF (XC%ID(IXC)==ID_XC_PW91) THEN
         CALL GGAXCS(D1,D2,DD1,DD2,EXC,EXCD1,EXCD2,EXCQ1,EXCQ2)

         D = D1 + D2
         EXC    = 2._q* EXC / D
         EXCD1  = 2._q* EXCD1
         EXCD2  = 2._q* EXCD2
         EXCQ1  = 2._q* EXCQ1
         EXCQ2  = 2._q* EXCQ2
      ELSE IF ((XC%ID(IXC)==ID_XC_PBE).OR. &
             & (XC%ID(IXC)==ID_XC_RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_REVPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTPBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB88PBE).OR. &
             & (XC%ID(IXC)==ID_XC_OPTB86BPBE).OR. &
             & (XC%ID(IXC)==ID_XC_PW86RPBE).OR. &
             & (XC%ID(IXC)==ID_XC_CX).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_RA).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_PL).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN

        IF (XC%ID(IXC)==ID_XC_PBE) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        D=2*D1

#undef exchange_RHOU_RHOD

#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif
        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif

        ELSE
           S=0.0_q
        ENDIF
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
          &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
        ELSE
           CALL CALC_EXCHWPBE_SP(XC,D,S, &
          &  EXLDA1,EXDLDA1,EXLDA_SR1,EXLDA_LR1,EXDLDA_SR1,EXDLDA_LR1, &
          &  EXWPBE_SR1,EXWPBE_LR1,EXWPBED_SR1,EXWPBED_LR1,EXWPBEDD_SR1,EXWPBEDD_LR1)
        ENDIF
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2
#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif
        ELSE
           S=0.0_q
        ENDIF

        IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
          &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
        ELSE
           CALL CALC_EXCHWPBE_SP(XC,D,S, &
          &  EXLDA2,EXDLDA2,EXLDA_SR2,EXLDA_LR2,EXDLDA_SR2,EXDLDA_LR2, &
          &  EXWPBE_SR2,EXWPBE_LR2,EXWPBED_SR2,EXWPBED_LR2,EXWPBEDD_SR2,EXWPBEDD_LR2)
        ENDIF
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
        ! Add LDA contributions as well
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*ALDAX+ECLDA
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*ALDAX+ECD1LDA
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*ALDAX+ECD2LDA
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX +(EXWPBE_SR1+EXWPBE_SR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX +EXWPBED_SR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX +EXWPBED_SR2+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ*AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX + (EXWPBE_SR1+EXWPBE_SR2)/(2*D)*BGGAX + ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX + EXWPBED_SR1*BGGAX + ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX + EXWPBED_SR2*BGGAX + ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX + EXWPBEDD_SR1*BGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX + EXWPBEDD_SR2*BGGAX
                 ECQ =ECQ*AGGAC
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1+EXWPBE_SR2)/(2*D)*AGGAX +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_SR1*AGGAX +EXWPBED_LR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_SR2*AGGAX +EXWPBED_LR2+ECD2LDA
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ*AGGAC
              ENDIF
           ENDIF
        ELSE
        ! Do not add LDA contribution
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC*  AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX+(EXWPBED_SR1-EXDLDA_SR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX+(EXWPBED_SR2-EXDLDA_SR2)
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ *AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX &
                &                +(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*BGGAX 
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX + (EXWPBED_SR1-EXDLDA_SR1)*BGGAX
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX + (EXWPBED_SR2-EXDLDA_SR2)*BGGAX
                 EXCQ1=EXWPBEDD_LR1*AGGAX + EXWPBEDD_SR1*BGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX + EXWPBEDD_SR2*BGGAX
                 ECQ =ECQ *AGGAC                                  
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_SR1-EXDLDA_SR1)*AGGAX+(EXWPBED_LR1-EXDLDA_LR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_SR2-EXDLDA_SR2)*AGGAX+(EXWPBED_LR2-EXDLDA_LR2)
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ *AGGAC
              ENDIF
           ENDIF
        ENDIF
        ! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 2*ECQ
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_X) THEN

        IF (XC%ID(IXC)==ID_XC_PBE_X) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        D=2*D1

#undef exchange_RHOU_RHOD

#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif
        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif

        ELSE
           S=0.0_q
        ENDIF
        CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
       &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

        D=2*D2
#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif
        ELSE
           S=0.0_q
        ENDIF

        CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
       &     ukfactor,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

        D=D1+D2

        IF (LLDA) THEN
        ! Add LDA contributions as well
           EXC  =(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*ALDAX
           EXCD1=(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*ALDAX
           EXCD2=(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*ALDAX
           EXCQ1=EXCQ1 *AGGAX
           EXCQ2=EXCQ2 *AGGAX
        ELSE
        ! Do not add LDA contribution
           EXC  =(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
           EXCD1=(EXCD1-EXDLDA1)*AGGAX
           EXCD2=(EXCD2-EXDLDA2)*AGGAX
           EXCQ1=EXCQ1 *AGGAX
           EXCQ2=EXCQ2 *AGGAX
        ENDIF
        ! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 0._q
      ELSE IF (XC%ID(IXC)==ID_XC_PBE_C) THEN

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)

        IF (LLDA) THEN
        ! Add LDA contributions as well
           EXC  =EC  *AGGAC + ECLDA
           EXCD1=ECD1*AGGAC + ECD1LDA
           EXCD2=ECD2*AGGAC + ECD2LDA
           ECQ  =ECQ *AGGAC
        ELSE
        ! Do not add LDA contribution
           EXC  =EC*  AGGAC
           EXCD1=ECD1*AGGAC
           EXCD2=ECD2*AGGAC
           ECQ  =ECQ *AGGAC
        ENDIF
        ! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 0._q
        EXCQ2  = 0._q
        ECQ    = 2*ECQ
! jP
      ELSE IF ((XC%ID(IXC)==ID_XC_B3).OR.(XC%ID(IXC)==ID_XC_B5)) THEN

        CALL  B3LYPXCS(XC,D1,D2,DD1,DD2,DDA,EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA)

        D = D1 + D2
        EXC    = 2.0_q*EXC/D
        EXCD1  = 2.0_q*EXCD1
        EXCD2  = 2.0_q*EXCD2
        EXCQ1  = 2.0_q*EXCDD1
        EXCQ2  = 2.0_q*EXCDD2
        ECQ    = 2.0_q*EXCDDA

!       WRITE(91,'(7F12.7)') D1, D2, EXC, EXCD1, EXCD2, EXCDD1, EXCDD2
! jP
!aem AM05 spin formulation added
      ELSE IF (XC%ID(IXC)==ID_XC_AM05) THEN

        !aem DD sometimes comes in negative.
        !aem Since AM05 assumes that DD actually IS |grad rho|,
        !aem and thus DD always should be positive,
        !aem we need to use ABS(DD) as input to AM05.
        !aem However, my routine needs to give the same symmetries out as
        !aem PBE and PW91, therefor the sign-compensation of EXCDD.

        ider = 1
        CALL AM05(XC,D1,D2,ABS(DD1),ABS(DD2),ider, &
                  FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
                  FXCDD1,FXCDD2)

        D = D1 + D2
        IF (LLDA) THEN
          EXC = FXC/D
          EXCD1 = FXCD1
          EXCD2 = FXCD2
        ELSE
          EXC = (FXC - FXCLDA)/D
          EXCD1 = FXCD1 - VXCLDA1
          EXCD2 = FXCD2 - VXCLDA2
        ENDIF
        EXCQ1 = SIGN(1.0_q,DD1)*FXCDD1
        EXCQ2 = SIGN(1.0_q,DD2)*FXCDD2
        ECQ = 0._q

        ! Hartree -> Rydberg conversion
        EXC    = 2._q* EXC
        EXCD1  = 2._q* EXCD1
        EXCD2  = 2._q* EXCD2
        EXCQ1  = 2._q* EXCQ1
        EXCQ2  = 2._q* EXCQ2
        ECQ    = 2._q* ECQ
!aem AM05 spin formulation added
! jP: adding PBEsol
!     Note that this is a plain copy-paste - but it offers
!     the possibility to add a screened PBEsol version later,
!     for a possible PBEsol hybrid-functional  - closely related to HSE
!     Further note: necessary condition: fitting the PBEsol-exchange hole!!!
      ELSE IF (XC%ID(IXC)==ID_XC_PBESOL) THEN
           ukfactor=1.0_q


        D=2*D1

#undef exchange_RHOU_RHOD

#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif
        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif

        ELSE
           S=0.0_q
        ENDIF
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(XC,D,S, &
          &  EXLDA1,EXDLDA1,EXLDA_SR1,EXLDA_LR1,EXDLDA_SR1,EXDLDA_LR1, &
          &  EXWPBE_SR1,EXWPBE_LR1,EXWPBED_SR1,EXWPBED_LR1,EXWPBEDD_SR1,EXWPBEDD_LR1)
        ENDIF
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2
#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif
        ELSE
           S=0.0_q
        ENDIF

        IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(XC,D,S, &
          &  EXLDA2,EXDLDA2,EXLDA_SR2,EXLDA_LR2,EXDLDA_SR2,EXDLDA_LR2, &
          &  EXWPBE_SR2,EXWPBE_LR2,EXWPBED_SR2,EXWPBED_LR2,EXWPBEDD_SR2,EXWPBEDD_LR2)
        ENDIF
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbesol(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
! Add LDA contributions as well
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*ALDAX+ECLDA
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*ALDAX+ECD1LDA
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*ALDAX+ECD2LDA
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX +(EXWPBE_SR1+EXWPBE_SR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX +EXWPBED_SR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX +EXWPBED_SR2+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ*AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC*AGGAC  +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX + (EXWPBE_SR1+EXWPBE_SR2)/(2*D)*BGGAX+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX + EXWPBED_SR1*BGGAX+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX + EXWPBED_SR2*BGGAX+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX + EXWPBEDD_SR1*BGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX + EXWPBEDD_SR2*BGGAX
                 ECQ =ECQ*AGGAC                                
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1+EXWPBE_SR2)/(2*D)*AGGAX +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_SR1*AGGAX +EXWPBED_LR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_SR2*AGGAX +EXWPBED_LR2+ECD2LDA
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ*AGGAC
              ENDIF
           ENDIF
        ELSE
! Do not add LDA contribution
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC*  AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX+(EXWPBED_SR1-EXDLDA_SR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX+(EXWPBED_SR2-EXDLDA_SR2)
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ *AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC=EC*AGGAC+ (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &              (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*BGGAX 
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX + (EXWPBED_SR1-EXDLDA_SR1)*BGGAX
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX + (EXWPBED_SR2-EXDLDA_SR2)*BGGAX
                 EXCQ1=EXWPBEDD_LR1*AGGAX + EXWPBEDD_SR1*BGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX + EXWPBEDD_SR2*BGGAX
                 ECQ =ECQ *AGGAC                                  
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_SR1-EXDLDA_SR1)*AGGAX+(EXWPBED_LR1-EXDLDA_LR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_SR2-EXDLDA_SR2)*AGGAX+(EXWPBED_LR2-EXDLDA_LR2)
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ *AGGAC
              ENDIF
           ENDIF
        ENDIF
! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 2*ECQ
! jH-new otherwise spin-polarized not possible
      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_RA).OR.(XC%ID(IXC)==ID_XC_ACFDT_PL)) THEN
        EXC    = 0._q
        EXCD1  = 0._q
        EXCD2  = 0._q
        EXCQ1  = 0._q
        EXCQ2  = 0._q
        ECQ    = 0._q
! jP: adding PBEsol

!BEEF
#ifdef libbeef
      ELSE IF (XC%ID(IXC)==ID_XC_BEEF) THEN
        D = D1+D2
        IF (D>1E-10_q) THEN
           ZETA=(D1-D2)/D
           ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
           IF (D1>1E-10_q) THEN
            DDSQ = 4._q*DD1*DD1
            CALL BEEFX(2._q*D1, DDSQ, EXBEEF1, EXBEEFD1, EXCQ1, LLDA)
            IF (DDSQ<1E-10_q) THEN
             EXBEEFD1 = 0._q
             EXCQ1 = 0._q
            END IF
           ELSE
            EXBEEF1 = 0._q
            EXBEEFD1 = 0._q
            EXCQ1 = 0._q
           END IF
           IF (D2>1E-10_q) THEN
            DDSQ = 4._q*DD2*DD2
            CALL BEEFX(2._q*D2, DDSQ, EXBEEF2, EXBEEFD2, EXCQ2, LLDA)
            IF (DDSQ<1E-10_q) THEN
             EXBEEFD2 = 0._q
             EXCQ2 = 0._q
            END IF
           ELSE
            EXBEEF2 = 0._q
            EXBEEFD2 = 0._q
            EXCQ2 = 0._q
           END IF
           CALL BEEFLOCALCORRSPIN(D, ZETA, DDA*DDA, ECBEEF, ECBEEFD1, ECBEEFD2, ECQ, LLDA)

           EXC = (AGGAX*EXBEEF1 + AGGAX*EXBEEF2 + AGGAC*2._q*ECBEEF)/D
           EXCD1 = 2._q*(AGGAX*EXBEEFD1 + AGGAC*ECBEEFD1)
           EXCD2 = 2._q*(AGGAX*EXBEEFD2 + AGGAC*ECBEEFD2)
           EXCQ1 = 4._q*AGGAX*EXCQ1*DD1
           EXCQ2 = 4._q*AGGAX*EXCQ2*DD2
           ECQ = 2._q*AGGAC*ECQ*DDA
        ELSE
           EXC   = 0._q
           EXCD1 = 0._q
           EXCD2 = 0._q
           EXCQ1 = 0._q
           EXCQ2 = 0._q
           ECQ   = 0._q
        END IF
#endif

      ELSE IF ((XC%ID(IXC)==ID_XC_ACFDT_03).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_05).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_10).OR. &
             & (XC%ID(IXC)==ID_XC_ACFDT_20)) THEN
! functionals for range-separated ACFDT (LDA - short range RPA)
! N.B. even if these functionals would support spin polarization,
! this would be a dummy stub since these functionals are not
! gradient corrected.
__CPU    CALL vtutor%error('GGA = 03 | 05 | 10 | 20  does not support spin polarization!')
#ifdef USELIBXC
      ELSE IF (XC%ID(IXC)==ID_XC_LIBXC) THEN
         XC_DENS(1) = D1
         XC_DENS(2) = D2
         XC_SIGMA(1) = DD1*DD1
         XC_SIGMA(2) = 0.5_q*(DDA*DDA - DD1*DD1 - DD2*DD2)
         XC_SIGMA(3) = DD2*DD2
         XC_SIGMAZERO = 0._q
         CX = -(3._q/4._q)*(3._q/PI)**THRD
         EXLDA = 0.5_q*CX*((2._q*D1)**THRD4 + (2._q*D2)**THRD4)/(D1+D2)
         EXDLDA1 = THRD4*CX*(2._q*D1)**THRD
         EXDLDA2 = THRD4*CX*(2._q*D2)**THRD

         IF ((LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_LDA) .OR. (LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_HYB_LDA)) THEN
            IF (LIBXC_NOEXC(IXC)) THEN
               CALL XC_F03_LDA_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_VRHO)
               XC_EXC = 0._q
            ELSE
               CALL XC_F03_LDA_EXC_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_EXC,XC_VRHO)
            ENDIF
            IF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE) .AND. (.NOT.LIBXC_HYB(IXC))) THEN
               EXC   = ALDAX*XC_EXC(1)
               EXCD1 = ALDAX*XC_VRHO(1)
               EXCD2 = ALDAX*XC_VRHO(2)
            ELSEIF (LIBXC_KIND(IXC) .EQ. XC_CORRELATION) THEN
               EXC   = ALDAC*XC_EXC(1)
               EXCD1 = ALDAC*XC_VRHO(1)
               EXCD2 = ALDAC*XC_VRHO(2)
            ELSEIF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE_CORRELATION) .OR. LIBXC_HYB(IXC)) THEN
               EXC   = XC_EXC(1)
               EXCD1 = XC_VRHO(1)
               EXCD2 = XC_VRHO(2)
            ENDIF
         ELSEIF ((LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_GGA) .OR. (LIBXC_FAMILY(IXC) .EQ. LIBXC_FAMILY_HYB_GGA)) THEN
            IF (LIBXC_NOEXC(IXC)) THEN
               CALL XC_F03_GGA_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_SIGMA,XC_VRHO,XC_VSIGMA)
               XC_EXC = 0._q
            ELSE
               CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_SIGMA,XC_EXC,XC_VRHO,XC_VSIGMA)
            ENDIF
            IF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE) .AND. (.NOT.LIBXC_HYB(IXC))) THEN
               EXC   = AGGAX*(XC_EXC(1)-EXLDA) + ALDAX*EXLDA
               EXCD1 = AGGAX*(XC_VRHO(1)-EXDLDA1) + ALDAX*EXDLDA1
               EXCD2 = AGGAX*(XC_VRHO(2)-EXDLDA2) + ALDAX*EXDLDA2
               EXCQ1 = AGGAX*(XC_VSIGMA(1)*2._q*DD1 - XC_VSIGMA(2)*DD1)
               EXCQ2 = AGGAX*(XC_VSIGMA(3)*2._q*DD2 - XC_VSIGMA(2)*DD2)
               ECQ   = AGGAX*XC_VSIGMA(2)*DDA
            ELSEIF (LIBXC_KIND(IXC) .EQ. XC_CORRELATION) THEN
               IF (LIBXC_NOEXC(IXC)) THEN
                  CALL XC_F03_GGA_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_SIGMAZERO,XC_VRHOZERO,XC_VSIGMAZERO)
                  XC_EXCZERO = 0._q
               ELSE
                  CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_SIGMAZERO,XC_EXCZERO,XC_VRHOZERO,XC_VSIGMAZERO)
               ENDIF
               EXC   = AGGAC*(XC_EXC(1)-XC_EXCZERO(1)) + ALDAC*XC_EXCZERO(1)
               EXCD1 = AGGAC*(XC_VRHO(1)-XC_VRHOZERO(1)) + ALDAC*XC_VRHOZERO(1)
               EXCD2 = AGGAC*(XC_VRHO(2)-XC_VRHOZERO(2)) + ALDAC*XC_VRHOZERO(2)
               EXCQ1 = AGGAC*(XC_VSIGMA(1)*2._q*DD1 - XC_VSIGMA(2)*DD1)
               EXCQ2 = AGGAC*(XC_VSIGMA(3)*2._q*DD2 - XC_VSIGMA(2)*DD2)
               ECQ   = AGGAC*XC_VSIGMA(2)*DDA
            ELSEIF ((LIBXC_KIND(IXC) .EQ. XC_EXCHANGE_CORRELATION) .OR. LIBXC_HYB(IXC)) THEN
               EXC   = XC_EXC(1)
               EXCD1 = XC_VRHO(1)
               EXCD2 = XC_VRHO(2)
               EXCQ1 = XC_VSIGMA(1)*2._q*DD1 - XC_VSIGMA(2)*DD1
               EXCQ2 = XC_VSIGMA(3)*2._q*DD2 - XC_VSIGMA(2)*DD2
               ECQ   = XC_VSIGMA(2)*DDA
            ENDIF
         ENDIF
         !Conversion to Rydbergs
         EXC   = 2._q*EXC
         EXCD1 = 2._q*EXCD1
         EXCD2 = 2._q*EXCD2
         EXCQ1 = 2._q*EXCQ1
         EXCQ2 = 2._q*EXCQ2
         ECQ   = 2._q*ECQ
#endif
      ELSE
__CPU    CALL vtutor%bug('internal ERROR GGASPIN_DRIVER: Wrong functional, scheme not implemented!', __FILE__, __LINE__)
      ENDIF

      EXC_TMP   = EXC_TMP   + XC%COEFF(IXC)*EXC
      EXCD1_TMP = EXCD1_TMP + XC%COEFF(IXC)*EXCD1
      EXCD2_TMP = EXCD2_TMP + XC%COEFF(IXC)*EXCD2
      EXCQ1_TMP = EXCQ1_TMP + XC%COEFF(IXC)*EXCQ1
      EXCQ2_TMP = EXCQ2_TMP + XC%COEFF(IXC)*EXCQ2
      ECQ_TMP   = ECQ_TMP   + XC%COEFF(IXC)*ECQ

      ENDIF

      ENDDO NXCLOOP4

      EXC   = EXC_TMP
      EXCD1 = EXCD1_TMP
      EXCD2 = EXCD2_TMP
      EXCQ1 = EXCQ1_TMP
      EXCQ2 = EXCQ2_TMP
      ECQ   = ECQ_TMP

      RETURN
      END SUBROUTINE GGASPIN_DRIVER


!************************ SUBROUTINE METAGGA_DRIVER ***********************************
!
!> This subroutine calls the subroutines that calculate xc quantities
!> (spin-polarized meta-GGA) for the energy, potential, forces, stress tensor, etc.
!>
!> Unlike in GGA_DRIVER and GGASPIN_DRIVER, the LDA contribution is included directly.
!
!**************************************************************************************

      SUBROUTINE METAGGA_DRIVER(&
     &   XC,RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,LAPLUP,LAPLDW,TAUUP,TAUDW, &
     &   EXC,dEXCdRHOup,dEXCdRHOdw,dEXCdABSNABup,dEXCdABSNABdw,dEXCdABSNAB, &
     &   dEXCdTAUup,dEXCdTAUdw,dEXCdLAPup,dEXCdLAPdw,CVMBJ)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE setexm_struct_def, only : xc_info
      USE ldalib, only : ex
      USE ggalib, only : corpbe,exchpbe
      USE mggalib
      USE xc_name
      USE xc_family_name
#ifdef USELIBXC
      USE xc_f03_lib_m
      USE setexm_struct_def, ONLY : LIBXC_AMGGAXC,LIBXC_FUNC_SPN,LTBOUNDLIBXC,LIBXC_NOEXC
#endif
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(METAGGA_DRIVER)

      TYPE (xc_info) XC

      REAL(q) RHOUP,RHODW
      REAL(q) ABSNABUP,ABSNABDW,ABSNAB
      REAL(q) LAPLUP,LAPLDW
      REAL(q) TAUUP,TAUDW

      REAL(q) EXC
      REAL(q) dEXCdRHOup,dEXCdRHOdw
      REAL(q) dEXCdABSNABup,dEXCdABSNABdw,dEXCdABSNAB
      REAL(q) dEXCdTAUup,dEXCdTAUdw
      REAL(q) dEXCdLAPup,dEXCdLAPdw

      REAL(q) EXC_TMP
      REAL(q) dEXCdRHOup_TMP,dEXCdRHOdw_TMP
      REAL(q) dEXCdABSNABup_TMP,dEXCdABSNABdw_TMP,dEXCdABSNAB_TMP
      REAL(q) dEXCdTAUup_TMP,dEXCdTAUdw_TMP
      REAL(q) dEXCdLAPup_TMP,dEXCdLAPdw_TMP

      REAL(q) CVMBJ

      ! local variables
      REAL(q) Ex_TPSS,Ec_TPSS
      REAL(q) Ex_revTPSS,Ec_revTPSS
      REAL(q) Ex_SCAN,Ec_SCAN
      REAL(q) Ex_M06,Ec_M06
      REAL(q) Ex_MSX,Ec_MSX
      REAL(q) Ex_SRTM,Ec_RTM
      REAL(q) Ex_TASK,Ec_CC,Ex_LAK,Ec_LAK
      REAL(q) VXD1,VXD2,VXDD1,VXDD2,AMUXD1,AMUXD2
      REAL(q) VCD1,VCD2,VCDD1,VCDD2,AMUCD1,AMUCD2
      REAL(q) VCDD12

      ! VMBJ related variables
      REAL(q), PARAMETER :: THRD=1._q/3._q
      REAL(q) D,DTHRD,RS,ZETA
      REAL(q) ECLDA,ECD1LDA,ECD2LDA,EC,ECD1,ECD2,ECQ

      ! for SREGTM2L
      REAL(q) QQ,TAU1,TAU2,TAUTOT,TD1,TD2,TDD1,TDD2,TDL1,TDL2
      REAL(q) ECDT,ECL1,ECL2,ECQ1,ECQ2,EXCL1,EXCL2

      ! for PBETEST
      REAL(q) FK,G,S,SK,T
      REAL(q) EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2

      REAL(q) AMGGAX,AMGGAC

      INTEGER :: IXC

#ifdef USELIBXC
      REAL(q) :: XC_DENS(2),XC_SIGMA(3),XC_TAU(2),XC_LAPL(2),TAUWUP,TAUWDW
      REAL(q) :: XC_EXC(1),XC_VRHO(2),XC_VSIGMA(3),XC_VTAU(2),XC_VLAPL(2)
#endif

GPACC !$ACC ROUTINE(VrevTPSSx) SEQ
GPACC !$ACC ROUTINE(VrevTPSSc) SEQ
GPACC !$ACC ROUTINE(VTPSSx) SEQ
GPACC !$ACC ROUTINE(VTPSSc) SEQ
GPACC !$ACC ROUTINE(VMBJ_PROYNOV) SEQ
GPACC !$ACC ROUTINE(CORPBE) SEQ
GPACC !$ACC ROUTINE(VM06x) SEQ
GPACC !$ACC ROUTINE(VM06c) SEQ
GPACC !$ACC ROUTINE(VMSXx) SEQ
GPACC !$ACC ROUTINE(VMSXc) SEQ
GPACC !$ACC ROUTINE(VSCANx) SEQ
GPACC !$ACC ROUTINE(VSCANc) SEQ
GPACC !$ACC ROUTINE(VR2SCANc) SEQ
GPACC !$ACC ROUTINE(VR2SCANx) SEQ
GPACC !$ACC ROUTINE(V1SRTMx) SEQ
GPACC !$ACC ROUTINE(V2SRTMx) SEQ
GPACC !$ACC ROUTINE(V3SRTMx) SEQ
GPACC !$ACC ROUTINE(VRTMc) SEQ
GPACC !$ACC ROUTINE(VTASKx) SEQ
GPACC !$ACC ROUTINE(VLAKx) SEQ
GPACC !$ACC ROUTINE(VLAKc) SEQ
GPACC !$ACC ROUTINE(scan_orb_free_xc) SEQ
GPACC !$ACC ROUTINE(EXCHPBE) SEQ
GPACC !$ACC ROUTINE(XV2SRTML) SEQ
GPACC !$ACC ROUTINE(PCREPMGGA) SEQ
GPACC !$ACC ROUTINE(CRTML) SEQ
#ifdef USELIBXC
GPACC !$ACC ROUTINE(XC_F03_LDA_VXC) SEQ
GPACC !$ACC ROUTINE(XC_F03_LDA_EXC_VXC) SEQ
GPACC !$ACC ROUTINE(XC_F03_GGA_VXC) SEQ
GPACC !$ACC ROUTINE(XC_F03_GGA_EXC_VXC) SEQ
GPACC !$ACC ROUTINE(XC_F03_MGGA_VXC) SEQ
GPACC !$ACC ROUTINE(XC_F03_MGGA_EXC_VXC) SEQ
#endif

      AMGGAX=XC%AMGGAX
      AMGGAC=XC%AMGGAC

      EXC=0._q
      dEXCdRHOup=0._q
      dEXCdRHOdw=0._q
      dEXCdABSNABup=0._q
      dEXCdABSNABdw=0._q
      dEXCdABSNAB=0._q
      dEXCdTAUup=0._q
      dEXCdTAUdw=0._q
      dEXCdLAPup=0._q
      dEXCdLAPdw=0._q

      NXCLOOP: DO IXC=1,XC%NXC

      IF (XC%FAMILY(IXC)/=XC_FAM_MGGA) CYCLE

      EXC_TMP           = 0._q
      dEXCdRHOup_TMP    = 0._q
      dEXCdRHOdw_TMP    = 0._q
      dEXCdABSNABup_TMP = 0._q
      dEXCdABSNABdw_TMP = 0._q
      dEXCdABSNAB_TMP   = 0._q
      dEXCdTAUup_TMP    = 0._q
      dEXCdTAUdw_TMP    = 0._q
      dEXCdLAPup_TMP    = 0._q
      dEXCdLAPdw_TMP    = 0._q

      IF (XC%ID(IXC)==ID_XC_RTPSS) THEN
      ! revTPSS
         ! Exchange
         CALL VrevTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

         ! Correlation
         CALL VrevTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_revTPSS+AMGGAC*Ec_revTPSS)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1+AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2+AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_RTPSS_X) THEN
      ! revTPSS
         ! Exchange
         CALL VrevTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

         EXC_TMP=AMGGAX*Ex_revTPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_RTPSS_C) THEN
      ! revTPSS
         ! Correlation
         CALL VrevTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_revTPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_TPSS) THEN
      ! TPSS
         ! Exchange
         CALL VTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_TPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

         ! Correlation
         CALL VTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_TPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_TPSS+AMGGAC*Ec_TPSS)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1+AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2+AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_TPSS_X) THEN
      ! TPSS
         ! Exchange
         CALL VTPSSx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_TPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

         EXC_TMP=AMGGAX*Ex_TPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_TPSS_C) THEN
      ! TPSS
         ! Correlation
         CALL VTPSSc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_TPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_TPSS/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAC*VCDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MBJ).OR. &
              (XC%ID(IXC)==ID_XC_LMBJ)) THEN
      ! MBJ

         CALL VMBJ_PROYNOV(RHOUP,ABSNABUP,LAPLUP,2._q*TAUUP,CVMBJ,dEXCdRHOup_TMP)
         CALL VMBJ_PROYNOV(RHODW,ABSNABDW,LAPLDW,2._q*TAUDW,CVMBJ,dEXCdRHOdw_TMP)

         ! get the LDA correlation potential (using a call to CORPBE)
         D=RHOUP+RHODW
         DTHRD=exp(log(D)*THRD)
         RS=(0.75_q/PI)**THRD/DTHRD
         ZETA=(RHOUP-RHODW)/D
         ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

         CALL CORPBE(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,0._q,0._q,0._q,EC,ECD1,ECD2,ECQ,.FALSE.)

         ! add the LDA correlation potential
         dEXCdRHOup_TMP=dEXCdRHOup_TMP+ECD1LDA
         dEXCdRHOdw_TMP=dEXCdRHOdw_TMP+ECD2LDA

         ! LDA exchange-correlation energy
         EXC_TMP=ECLDA+0.25_q*EX((0.75_q/PI)**THRD/RHOUP**THRD,2,.TRUE.)+0.25_q*EX((0.75_q/PI)**THRD/RHODW**THRD,2,.TRUE.)
      ELSEIF (XC%ID(IXC)==ID_XC_M06L) THEN
      ! M06-L
         ! Exchange
         CALL VM06x(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_M06,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,1)

         ! Correlation
         CALL VM06c(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_M06,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2,1)

         EXC_TMP=(AMGGAX*Ex_M06+AMGGAC*Ec_M06)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1+AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2+AMGGAC*VCDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_M06L_X) THEN
      ! M06-L
         ! Exchange
         CALL VM06x(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_M06,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,1)

         EXC_TMP=AMGGAX*Ex_M06/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_M06L_C) THEN
      ! M06-L
         ! Correlation
         CALL VM06c(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_M06,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2,1)

         EXC_TMP=AMGGAC*Ec_M06/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAC*VCDD1
         dEXCdABSNABdw_TMP=AMGGAC*VCDD2
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MS0).OR. &
            & (XC%ID(IXC)==ID_XC_MS1).OR. &
            & (XC%ID(IXC)==ID_XC_MS2).OR. &
            & (XC%ID(IXC)==ID_XC_MSPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_MSRPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_MSB86BL).OR. &
            & (XC%ID(IXC)==ID_XC_RMSRPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_RMSPBEL).OR. &
            & (XC%ID(IXC)==ID_XC_RMSB86BL)) THEN
         ! Exchange
         CALL VMSXx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_MSX,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%PARAM(IXC,3),XC%PARAM(IXC,4),XC%ID(IXC))

         ! Correlation
         CALL VMSXc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_MSX,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_MSX+AMGGAC*Ec_MSX)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MS0_X).OR. &
            & (XC%ID(IXC)==ID_XC_MS1_X).OR. &
            & (XC%ID(IXC)==ID_XC_MS2_X)) THEN
         ! Exchange
         CALL VMSXx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_MSX,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%PARAM(IXC,3),XC%PARAM(IXC,4),XC%ID(IXC))

         EXC_TMP=AMGGAX*Ex_MSX/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF ((XC%ID(IXC)==ID_XC_MS0_C).OR. &
            & (XC%ID(IXC)==ID_XC_MS1_C).OR. &
            & (XC%ID(IXC)==ID_XC_MS2_C)) THEN
         ! Correlation
         CALL VMSXc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_MSX,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_MSX/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCAN).OR. &
              (XC%ID(IXC)==ID_XC_RSCAN)) THEN
      ! SCAN
         ! Exchange
         CALL VSCANx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_SCAN,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,3),XC%ID(IXC))

         ! Correlation
         CALL VSCANc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_SCAN,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2, &
       &   XC%PARAM(IXC,2),XC%PARAM(IXC,3),XC%ID(IXC))

         ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SCAN+AMGGAC*Ec_SCAN)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCAN_X).OR. &
              (XC%ID(IXC)==ID_XC_RSCAN_X)) THEN
         ! Exchange
         CALL VSCANx(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ex_SCAN,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
       &   XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%ID(IXC))

         EXC_TMP=AMGGAX*Ex_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCAN_C).OR. &
              (XC%ID(IXC)==ID_XC_RSCAN_C)) THEN
         ! Correlation
         CALL VSCANc(&
       &   RHOUP,RHODW,ABSNABUP,ABSNABDW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_SCAN,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2,XC%PARAM(IXC,1),XC%PARAM(IXC,2),XC%ID(IXC))

         EXC_TMP=AMGGAC*Ec_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_R2SCAN) THEN
      ! r^2SCAN
         ! Exchange
         Ex_SCAN = 0.0_q
         CALL VR2SCANx(RHOUP,ABSNABUP,TAUUP,EX_SCAN,VXD1,VXDD1,AMUXD1)
         CALL VR2SCANx(RHODW,ABSNABDW,TAUDW,EX_SCAN,VXD2,VXDD2,AMUXD2)

         ! Correlation
         CALL VR2SCANc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_SCAN,VCD1,VCD2,&
      &                VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_SCAN+AMGGAC*Ec_SCAN)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_R2SCAN_X) THEN
      ! r^2SCAN
         ! Exchange
         Ex_SCAN = 0.0_q
         CALL VR2SCANx(RHOUP,ABSNABUP,TAUUP,EX_SCAN,VXD1,VXDD1,AMUXD1)
         CALL VR2SCANx(RHODW,ABSNABDW,TAUDW,EX_SCAN,VXD2,VXDD2,AMUXD2)

         EXC_TMP=AMGGAX*Ex_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_R2SCAN_C) THEN
      ! r^2SCAN
         ! Correlation
         CALL VR2SCANc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_SCAN,VCD1,VCD2,&
      &                VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_SCAN/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM1) THEN
         ! Exchange
         Ex_SRTM = 0.0_q
         CALL V1SRTMx(RHOUP,ABSNABUP,TAUUP,Ex_SRTM,VXD1,VXDD1,AMUXD1)
         CALL V1SRTMx(RHODW,ABSNABDW,TAUDW,Ex_SRTM,VXD2,VXDD2,AMUXD2)
         !Correlation rrTM 
         CALL VRTMc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_RTM,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
!
!        ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SRTM+AMGGAC*Ec_RTM)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM2) THEN
         ! Exchange
         Ex_SRTM = 0.0_q
         CALL V2SRTMx(RHOUP,ABSNABUP,TAUUP,Ex_SRTM,VXD1,VXDD1,AMUXD1)
         CALL V2SRTMx(RHODW,ABSNABDW,TAUDW,Ex_SRTM,VXD2,VXDD2,AMUXD2)
         !Correlation rrTM 
         CALL VRTMc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_RTM,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
!
!        ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SRTM+AMGGAC*Ec_RTM)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM3) THEN
         ! Exchange
         Ex_SRTM = 0.0_q
         CALL V3SRTMx(RHOUP,ABSNABUP,TAUUP,Ex_SRTM,VXD1,VXDD1,AMUXD1)
         CALL V3SRTMx(RHODW,ABSNABDW,TAUDW,Ex_SRTM,VXD2,VXDD2,AMUXD2)
         !Correlation rrTM 
         CALL VRTMc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW, &
       &   Ec_RTM,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
!
!        ! Sum everything
         EXC_TMP=(AMGGAX*Ex_SRTM+AMGGAC*Ec_RTM)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_TASK_X) THEN
         ! Exchange
         Ex_TASK = 0.0_q
         CALL VTASKx(RHOUP,ABSNABUP,TAUUP,EX_TASK,VXD1,VXDD1,AMUXD1)
         CALL VTASKx(RHODW,ABSNABDW,TAUDW,EX_TASK,VXD2,VXDD2,AMUXD2)

         EXC_TMP=AMGGAX*Ex_TASK/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_CC_C) THEN
         !Correlation 
         CALL VCCc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_CC,VCD1,VCD2,&
       &   VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_CC/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_LAK) THEN
         ! Exchange
         Ex_LAK = 0.0_q
         CALL VLAKx(RHOUP,ABSNABUP,TAUUP,EX_LAK,VXD1,VXDD1,AMUXD1)
         CALL VLAKx(RHODW,ABSNABDW,TAUDW,EX_LAK,VXD2,VXDD2,AMUXD2)
         !Correlation 
         CALL VLAKc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_LAK,VCD1,VCD2,&
       &   VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=(AMGGAX*Ex_LAK+AMGGAC*Ec_LAK)/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1+AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2+AMGGAC*VCD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAX*AMUXD1+AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2+AMGGAC*AMUCD2
      ELSEIF (XC%ID(IXC)==ID_XC_LAK_X) THEN
         ! Exchange
         Ex_LAK = 0.0_q
         CALL VLAKx(RHOUP,ABSNABUP,TAUUP,EX_LAK,VXD1,VXDD1,AMUXD1)
         CALL VLAKx(RHODW,ABSNABDW,TAUDW,EX_LAK,VXD2,VXDD2,AMUXD2)

         EXC_TMP=AMGGAX*Ex_LAK/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAX*VXD1
         dEXCdRHOdw_TMP=AMGGAX*VXD2
         dEXCdABSNABup_TMP=AMGGAX*VXDD1
         dEXCdABSNABdw_TMP=AMGGAX*VXDD2
         dEXCdTAUup_TMP=AMGGAX*AMUXD1
         dEXCdTAUdw_TMP=AMGGAX*AMUXD2
      ELSEIF (XC%ID(IXC)==ID_XC_LAK_C) THEN
         !Correlation 
         CALL VLAKc(RHOUP,RHODW,ABSNAB,TAUUP,TAUDW,Ec_LAK,VCD1,VCD2,&
       &   VCDD12,AMUCD1,AMUCD2)

         EXC_TMP=AMGGAC*Ec_LAK/(RHOUP+RHODW)
         dEXCdRHOup_TMP=AMGGAC*VCD1
         dEXCdRHOdw_TMP=AMGGAC*VCD2
         dEXCdABSNAB_TMP=AMGGAC*VCDD12
         dEXCdTAUup_TMP=AMGGAC*AMUCD1
         dEXCdTAUdw_TMP=AMGGAC*AMUCD2
      ELSEIF ((XC%ID(IXC)==ID_XC_SCANL).OR. &
              (XC%ID(IXC)==ID_XC_RSCANL).OR. &
              (XC%ID(IXC)==ID_XC_R2SCANL).OR. &
              (XC%ID(IXC)==ID_XC_OFR2)) THEN
        call scan_orb_free_xc(rhoup,rhodw,absnabup,absnabdw,absnab, &
        &  LAPLUP,LAPLDW,EXC_TMP,dEXCdRHOup_TMP,dEXCdRHOdw_TMP,dEXCdABSNABup_TMP,dEXCdABSNABdw_TMP, &
        &  dEXCdABSNAB_TMP,dEXCdLAPup_TMP,dEXCdLAPdw_TMP,XC%ID(IXC),XC%AMGGAX,XC%AMGGAC)
      ELSEIF (XC%ID(IXC)==ID_XC_SREGTM2L) THEN

         D=2*RHOUP
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABUP/(D*FK)
            QQ=LAPLUP/(2.0_q*FK*FK*D)
         ELSE
            S=0.0_q
            QQ=0.0_q
         ENDIF

         CALL XV2SRTML(D,DTHRD,S,QQ,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1,EXCL1)
         CALL PCREPMGGA(D,DTHRD,S,QQ,TAU1,TD1,TDD1,TDL1)
         TAU1=0.5_q*TAU1

         D=2*RHODW
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABDW/(D*FK)
            QQ=LAPLDW/(2.0_q*FK*FK*D)
         ELSE
            S=0.0_q
            QQ=0.0_q
         ENDIF

         CALL XV2SRTML(D,DTHRD,S,QQ,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2,EXCL2)
         CALL PCREPMGGA(D,DTHRD,S,QQ,TAU2,TD2,TDD2,TDL2)
         TAU2=0.5_q*TAU2

         D=RHOUP+RHODW
         DTHRD=exp(log(D)*THRD)
         RS=(0.75_q/PI)**THRD/DTHRD
         ZETA=(RHOUP-RHODW)/D
         ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
         FK=(3._q*PI*PI)**THRD*DTHRD
         G=(exp((2*THRD)*log(1._q+ZETA)) &
     &      +exp((2*THRD)*log(1._q-ZETA)))/2._q
         IF(D>1.E-10_q)THEN
            S=ABSNAB/(D*2._q*FK)
         ELSE
            S=0.0_q
         ENDIF
         TAUTOT=TAU1+TAU2

         CALL CRTML(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,FK,S,EC,ECD1,ECD2,ECQ,TAUTOT,ECDT)

         ECD1 = ECD1 + ECDT*TD1
         ECD2 = ECD2 + ECDT*TD2
         ECQ1 = ECDT*TDD1
         ECQ2 = ECDT*TDD2
         ECL1 = ECDT*TDL1
         ECL2 = ECDT*TDL2

         EXC_TMP = AMGGAX*(EXC1 + EXC2)/(2._q*D) + AMGGAC*(ECLDA+EC)
         dEXCdRHOup_TMP = AMGGAX*EXCD1 + AMGGAC*(ECD1LDA+ECD1)
         dEXCdRHOdw_TMP = AMGGAX*EXCD2 + AMGGAC*(ECD2LDA+ECD2)
         dEXCdABSNABup_TMP = AMGGAX*EXCQ1 + AMGGAC*ECQ1
         dEXCdABSNABdw_TMP = AMGGAX*EXCQ2 + AMGGAC*ECQ2
         dEXCdABSNAB_TMP = AMGGAC*ECQ
         dEXCdLAPup_TMP = AMGGAX*EXCL1 + AMGGAC*ECL1
         dEXCdLAPdw_TMP = AMGGAX*EXCL2 + AMGGAC*ECL2

#ifdef USELIBXC
      ELSEIF (XC%ID(IXC)==ID_XC_LIBXC) THEN

         XC_DENS(1) = max(RHOUP, 1.0e-10_q)
         XC_DENS(2) = max(RHODW, 1.0e-10_q)

         ABSNABUP = max(min(ABSNABUP, 5.0e9_q), 1.0e-10_q)
         ABSNABDW = max(min(ABSNABDW, 5.0e9_q), 1.0e-10_q)
         ABSNAB = max(min(ABSNAB, 1.0e10_q), 2.0e-10_q)
         XC_SIGMA(1) = ABSNABUP*ABSNABUP
         XC_SIGMA(2) = 0.5_q*(ABSNAB*ABSNAB - ABSNABUP*ABSNABUP - ABSNABDW*ABSNABDW)
         XC_SIGMA(3) = ABSNABDW*ABSNABDW

         XC_TAU(1) = max(min(TAUUP, 5.0e9_q), 1.0e-10_q)
         XC_TAU(2) = max(min(TAUDW, 5.0e9_q), 1.0e-10_q)
!!#ifdef XC_MAJOR_VERSION < 7
         IF (LTBOUNDLIBXC) THEN
            TAUWUP = 0.125_q*ABSNABUP**2/XC_DENS(1)
            TAUWDW = 0.125_q*ABSNABDW**2/XC_DENS(2)
            XC_TAU(1) = MAX(XC_TAU(1),TAUWUP)
            XC_TAU(2) = MAX(XC_TAU(2),TAUWDW)
         ENDIF
!!#endif
         XC_LAPL(1) = LAPLUP
         XC_LAPL(2) = LAPLDW

         IF (LIBXC_NOEXC(IXC)) THEN
            call XC_F03_MGGA_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_SIGMA,XC_LAPL,XC_TAU, &
           &   XC_VRHO,XC_VSIGMA,XC_VLAPL,XC_VTAU)
            XC_EXC = 0._q
         ELSE
            call XC_F03_MGGA_EXC_VXC(LIBXC_FUNC_SPN(IXC),1_qi8,XC_DENS,XC_SIGMA,XC_LAPL,XC_TAU, &
           &   XC_EXC,XC_VRHO,XC_VSIGMA,XC_VLAPL,XC_VTAU)
         ENDIF
         EXC_TMP           = LIBXC_AMGGAXC(IXC)*XC_EXC(1)
         dEXCdRHOup_TMP    = LIBXC_AMGGAXC(IXC)*XC_VRHO(1)
         dEXCdRHOdw_TMP    = LIBXC_AMGGAXC(IXC)*XC_VRHO(2)
         dEXCdABSNABup_TMP = LIBXC_AMGGAXC(IXC)*(XC_VSIGMA(1)*2._q*ABSNABUP - XC_VSIGMA(2)*ABSNABUP)
         dEXCdABSNABdw_TMP = LIBXC_AMGGAXC(IXC)*(XC_VSIGMA(3)*2._q*ABSNABDW - XC_VSIGMA(2)*ABSNABDW)
         dEXCdABSNAB_TMP   = LIBXC_AMGGAXC(IXC)*XC_VSIGMA(2)*ABSNAB
         dEXCdTAUup_TMP    = LIBXC_AMGGAXC(IXC)*XC_VTAU(1)
         dEXCdTAUdw_TMP    = LIBXC_AMGGAXC(IXC)*XC_VTAU(2)
         dEXCdLAPup_TMP    = LIBXC_AMGGAXC(IXC)*XC_VLAPL(1)
         dEXCdLAPdw_TMP    = LIBXC_AMGGAXC(IXC)*XC_VLAPL(2)
#endif
      ELSEIF (XC%ID(IXC)==ID_XC_PBETEST) THEN
         ! PBE, for testing mainly

         D=2*RHOUP
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABUP/(D*FK)
         ELSE
            S=0.0_q
         ENDIF
         CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1,1.0_q,XC%PARAM1,XC%PARAM2,XC%ID(IXC))
         D=2*RHODW
         DTHRD=exp(log(D)*THRD)
         FK=(3._q*PI*PI)**THRD*DTHRD
         IF(D>1.E-10_q)THEN
            S=ABSNABDW/(D*FK)
         ELSE
            S=0.0_q
         ENDIF
         CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2,1.0_q,XC%PARAM1,XC%PARAM2,XC%ID(IXC))

         D=RHOUP+RHODW
         DTHRD=exp(log(D)*THRD)
         RS=(0.75_q/PI)**THRD/DTHRD
         ZETA=(RHOUP-RHODW)/D
         ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
         FK=(3._q*PI*PI)**THRD*DTHRD
         SK = SQRT(4.0_q*FK/PI)
         G = (exp((2*THRD)*log(1._q+ZETA)) &
                 +exp((2*THRD)*log(1._q-ZETA)))/2._q
         T = ABSNAB/(D*2._q*SK*G)
         CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK,T,EC,ECD1,ECD2,ECQ,.TRUE.)

         EXC_TMP           = AMGGAX*(EXC1+EXC2)/(2._q*D)+AMGGAC*(ECLDA+EC)
         dEXCdRHOup_TMP    = AMGGAX*EXCD1+AMGGAC*(ECD1LDA+ECD1)
         dEXCdRHOdw_TMP    = AMGGAX*EXCD2+AMGGAC*(ECD2LDA+ECD2)
         dEXCdABSNABup_TMP = AMGGAX*EXCQ1
         dEXCdABSNABdw_TMP = AMGGAX*EXCQ2
         dEXCdABSNAB_TMP   = AMGGAC*ECQ
      ENDIF

      ! factor 2 for Hartree to Rydberg
      EXC           = EXC           + 2._q*XC%COEFF(IXC)*EXC_TMP
      dEXCdRHOup    = dEXCdRHOup    + 2._q*XC%COEFF(IXC)*dEXCdRHOup_TMP
      dEXCdRHOdw    = dEXCdRHOdw    + 2._q*XC%COEFF(IXC)*dEXCdRHOdw_TMP
      dEXCdABSNABup = dEXCdABSNABup + 2._q*XC%COEFF(IXC)*dEXCdABSNABup_TMP
      dEXCdABSNABdw = dEXCdABSNABdw + 2._q*XC%COEFF(IXC)*dEXCdABSNABdw_TMP
      dEXCdABSNAB   = dEXCdABSNAB   + 2._q*XC%COEFF(IXC)*dEXCdABSNAB_TMP
      dEXCdTAUup    = dEXCdTAUup    + 2._q*XC%COEFF(IXC)*dEXCdTAUup_TMP
      dEXCdTAUdw    = dEXCdTAUdw    + 2._q*XC%COEFF(IXC)*dEXCdTAUdw_TMP
      dEXCdLAPup    = dEXCdLAPup    + 2._q*XC%COEFF(IXC)*dEXCdLAPup_TMP
      dEXCdLAPdw    = dEXCdLAPdw    + 2._q*XC%COEFF(IXC)*dEXCdLAPdw_TMP

      ENDDO NXCLOOP

      RETURN
      END SUBROUTINE METAGGA_DRIVER

    END MODULE xc_driver
