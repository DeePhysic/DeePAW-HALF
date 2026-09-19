#define Zmethod3
#include "symbol.inc"

!************************ MODULE MGGALIB ***********************
!>
!> This module is the library of subroutines of MGGA functionals
!>
!***************************************************************

MODULE MGGALIB

      CONTAINS

!****************** SUBROUTINE VMBJ_PROYNOV ****************************
!
!> calculates modified Becke-Johnson exchange potential
!> according to Tran and Blaha, PRL 102, 226401 (2009), which is
!> based on Becke and Roussel, PRA 39, 3761 (1989) and
!> Becke and Johnson, J. Chem. Phys. 124, 221101 (2006)
!>
!> This routine uses the analytical represention of the Becke-Roussel
!> potential from Proynov, Chem. Phys. Lett. 455, 103 (2008)
!>
!> everything in Hartree unit and Bohr
!>
!> ATTENTION: All values are passed "as they are", i.e. including
!> possibly unphysical numerical errors (e.g. negative charge densities)
!> values need to be checked accordingly
!>
!> RHO        electron density (up or down)
!> GRHO       abs. val. gradient of density (up or down)
!> G2RHO      abs. val. second gradient of total density (up or down)
!> TAU        kinetic energy density (up or down)
!> XCCONST    fitted value by int{|nabla rho|/rho}dr3
!> VXBRJ      return value
!
!***********************************************************************

      SUBROUTINE VMBJ_PROYNOV(RHO,GRHO,G2RHO,TAU,XCCONST,VXBRJ)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE mbj, ONLY : MBJ_SMALL
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(VMBJ_PROYNOV)

      INTEGER I
      REAL(q) RHO,GRHO,G2RHO,TAU,TAUTMP,TAUW,XCCONST,VXBRJ
      REAL(q) A(1:3),B(0:5),C(0:5),D(0:5),E(0:5),YP(0:5)
      REAL(q) DD,G,P1,P2,QQ,X,Y,Z

      VXBRJ = -(6.0_q/PI*max(RHO,0.0_q))**(1.0_q/3.0_q)

      IF (RHO .GT. MBJ_SMALL) THEN

         TAUTMP = TAU
         TAUW = 0.25_q*GRHO**2/RHO
         TAU=MAX(TAU,TAUW)

         A(1) = 1.5255251812009530_q
         A(2) = 0.4576575543602858_q
         A(3) = 0.4292036732051034_q

         C(0) =   0.7566445420735584_q
         C(1) =  -2.6363977871370960_q
         C(2) =   5.4745159964232880_q
         C(3) = -12.657308127108290_q
         C(4) =   4.1250584725121360_q
         C(5) = -30.425133957163840_q

         B(0) =   0.4771976183772063_q
         B(1) =  -1.7799813494556270_q
         B(2) =   3.8433841862302150_q
         B(3) =  -9.5912050880518490_q
         B(4) =   2.1730180285916720_q
         B(5) = -30.425133851603660_q

         D(0) =    0.00004435009886795587_q
         D(1) =    0.58128653604457910_q
         D(2) =   66.742764515940610_q
         D(3) =  434.26780897229770_q
         D(4) =  824.7765766052239000_q
         D(5) = 1657.9652731582120_q

         E(0) =    0.00003347285060926091_q
         E(1) =    0.47917931023971350_q
         E(2) =   62.392268338574240_q
         E(3) =  463.14816427938120_q
         E(4) =  785.2360350104029000_q
         E(5) = 1657.962968223273000000_q

         DD = TAU - TAUW
         QQ = (G2RHO - 2.0_q*0.8_q*DD)/6.0_q

         IF (ABS(QQ) .GT. 1.0E-18_q) THEN

            Y = (2.0_q/3.0_q)*PI**(2.0_q/3.0_q)*RHO**(5.0_q/3.0_q)/QQ
            DO I=0, 5
               YP(I) = Y**I
            ENDDO

            IF (Y .LE. 0.0_q) THEN
               G = -ATAN(A(1)*Y + A(2)) + A(3)
               P1 = SUM(C(0:5)*YP(0:5))
               P2 = SUM(B(0:5)*YP(0:5))
            ELSEIF (Y .GT. 0.0_q) THEN
               Z = 2.085749716493756_q*Y
               G = LOG(SQRT(1.0_q + 1.0_q/Z**2) + 1.0_q/Z) + 2.0_q
               P1 = SUM(D(0:5)*YP(0:5))
               P2 = SUM(E(0:5)*YP(0:5))
            ENDIF

            X = G*P1/P2

            IF (ABS(X) .GT. 1.0E-18_q) THEN

               VXBRJ = -2.0_q*PI**(1.0_q/3.0_q)*RHO**(1.0_q/3.0_q)*EXP(X/3.0_q)/X*(1.0_q - EXP(-X) - 0.5_q*X*EXP(-X))

               IF (TAU .GE. 0.0_q) THEN
                  VXBRJ = XCCONST*VXBRJ + (3.0_q*XCCONST-2.0_q)*SQRT(5.0_q/12.0_q)/PI*SQRT(TAU/RHO)
               ELSE
                  VXBRJ = XCCONST*VXBRJ - (3.0_q*XCCONST-2.0_q)*SQRT(5.0_q/12.0_q)/PI*SQRT(ABS(TAU/RHO))
               ENDIF

            ENDIF
         ENDIF

         TAU = TAUTMP

      ENDIF

      RETURN
      END SUBROUTINE VMBJ_PROYNOV

!************************ SUBROUTINE VrevTPSSx *****************************
!
!> calculates the first order derivatives of Ex wrt n and |grad(n)|
!> Perdew et. al. PRL (2009)
!>
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> everything in Hartree units
!>
!> ATTENTION: Every values are passed "as they are", i.e. including
!> possibly unphysical numerical errors (e.g. negative charge densities)
!> values need to be checked accordingly
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VXD1 VXD2                    THE DERIVATIVES OF EX WRT n
!> VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n|
!> AMUXD1, AMUXD2               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VrevTPSSx(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP!$OMP DECLARE TARGET(VrevTPSSx)

! the following parameters are given by Perdew et.al.
! the coefficients for the exchange part
      PARAMETER (RKAPPA=0.804_q)
      PARAMETER (CFB=0.40_q)
      PARAMETER (CFC=2.35203946_q)
      PARAMETER (CFE=2.16769874_q)
      PARAMETER (CFMU=0.14_q)

! other parameters
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
! test
      PARAMETER (thresh=1.e-10_q)
      PARAMETER (xorder=12._q)
      PARAMETER (zinfinity=2._q)
! test
      VXD1=0._q;VXD2=0._q;
      VXDD1=0._q;VXDD2=0._q;
      AMUXD1=0._q;AMUXD2=0._q;


      CX1=10._q/81._q
      CX2=146._q/2025._q
      CX3=73._q/405._q
      CFE12=SQRT(CFE)

! Suspect that TAUWU and TAUWD are not well described. Use TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUWU_TEMP=0.125_q*DRU**2._q/RU
!     TAUWU_TEMP=MIN(TAUWU_TEMP,TAUU)

      TAUWD_TEMP=0.125_q*DRD**2._q/RD
!     TAUWD_TEMP=MIN(TAUWD_TEMP,TAUD)

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RU
      DRHO=TWO*DRU
      TAUW_RHO=TWO*TAUWU_TEMP
      TAU_RHO=TWO*TAUU

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P
      Z=TAUW_RHO/TAU_RHO
! test
      !Z=max(min(Z,0.9999999999999_q),1.E-20_q)
#ifdef Zmethod1
      IF (Z>1) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=1
      ELSE IF(Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
      ENDIF
#endif
#ifdef Zmethod2
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1._q/(1._q+Z**xorder)**(1._q/xorder) - Z**xorder/(1._q+Z**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=Z/(1._q+Z**xorder)**(1._q/xorder)
      ENDIF
#endif
#ifdef Zmethod3
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1._q/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder) - (Z/zinfinity)**xorder/(1._q+(Z/zinfinity)**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=Z/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder)
      ENDIF
#endif
! test
      Z2=Z*Z
      Z3=Z2*Z
      ZTF=0.6_q*Z
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5

!     ALPHA=(TAU_RHO-TAUW_RHO)/TAU0
      ALPHA=5._q/3._q*P*(1._q/Z-1._q)

      QB_NUM=9._q/20._q*(ALPHA-ONE)
      QB_DEN=SQRT(ONE+CFB*ALPHA*(ALPHA-ONE))
      QB=QB_NUM/QB_DEN+ TWO*P/THREE
      OPZ=ONE+Z2
      SFP=Z3/OPZ**TWO
      FQ=SQRT((ZTF**TWO+P2)/TWO)
      S1=(CX1+CFC*SFP)*P
      S2=CX2*QB**TWO
      S3=-CX3*FQ*QB
      S4=(CX1*P)**TWO/RKAPPA
      S5=TWO*CFE12*CX1*ZTF**TWO
      S6=CFE*CFMU*P2*P
      XN=S1+S2+S3+S4+S5+S6
      XD=(ONE+CFE12*P)**TWO
      X=XN/XD

!    GET THE VALUE FOR FX
      FX=ONE+RKAPPA-RKAPPA/(ONE+(X/RKAPPA))

!    NOW, DERIVATIVES COME
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q

      TAU_UNIF=THREE/10._q*(THREE*PI**TWO)**THRD2*RHO**THRD5

      DALPHAD=THRD5*(ONE/Z-ONE)*DPD-THRD5*P*DZD/Z**TWO
!     DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)+   &
!    &   DRHO**TWO/RHO**(11._q/THREE)*(10._q/9._q/(THREE*PI**TWO)**THRD2)

      DALPHADD=THRD5*(ONE/Z-ONE)*DPDD-THRD5*P*DZDD/Z**TWO
!     DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)

      DALPHADTAU=-THRD5*P*DZDTAU/Z**TWO
!     DALPHADTAU=ONE/TAU_UNIF

      DQB_NUM_D=9._q/20._q*DALPHAD
      DQB_NUM_DD=9._q/20._q*DALPHADD
      DQB_NUM_DTAU=9._q/20._q*DALPHADTAU
      C_DQB_DEN=ONE/TWO*CFB*(TWO*ALPHA-ONE)/QB_DEN
      DQB_DEN_D=C_DQB_DEN*DALPHAD
      DQB_DEN_DD=C_DQB_DEN*DALPHADD
      DQB_DEN_DTAU=C_DQB_DEN*DALPHADTAU

      DQBD=(DQB_NUM_D*QB_DEN-QB_NUM*DQB_DEN_D)/QB_DEN**TWO+THRD2*DPD
      DQBDD=(DQB_NUM_DD*QB_DEN-QB_NUM*DQB_DEN_DD)/QB_DEN**TWO+THRD2*DPDD
      DQBDTAU=(DQB_NUM_DTAU*QB_DEN-QB_NUM*DQB_DEN_DTAU)/QB_DEN**TWO+THRD2*DPDTAU

      C1_DS1=P*Z2*(THREE*CFC-CFC*Z2)/OPZ**THREE
      C2_DS1=(CX1+CFC*Z3/OPZ**TWO)
      DS1D=C1_DS1*DZD+C2_DS1*DPD
      DS1DD=C1_DS1*DZDD+C2_DS1*DPDD
      DS1DTAU=C1_DS1*DZDTAU+C2_DS1*DPDTAU

      C_DS2=CX2*TWO*QB
      DS2D=C_DS2*DQBD
      DS2DD=C_DS2*DQBDD
      DS2DTAU=C_DS2*DQBDTAU

      DFQD=ONE/TWO/FQ*(0.6_q**TWO*Z*DZD+P*DPD)
      DFQDD=ONE/TWO/FQ*(0.6_q**TWO*Z*DZDD+P*DPDD)
      DFQDTAU=ONE/TWO/FQ*(0.6_q**TWO*Z*DZDTAU+P*DPDTAU)
      DS3D=-CX3*(FQ*DQBD+QB*DFQD)
      DS3DD=-CX3*(FQ*DQBDD+QB*DFQDD)
      DS3DTAU=-CX3*(FQ*DQBDTAU+QB*DFQDTAU)

      C_S4=TWO/RKAPPA*CX1**TWO*P
      DS4D=C_S4*DPD
      DS4DD=C_S4*DPDD
      DS4DTAU=C_S4*DPDTAU

      C_S5=FOUR*CFE12*CX1*0.6_q**TWO*Z
      DS5D=C_S5*DZD
      DS5DD=C_S5*DZDD
      DS5DTAU=C_S5*DZDTAU

      C_S6=THREE*CFE*CFMU*P2
      DS6D=C_S6*DPD
      DS6DD=C_S6*DPDD
      DS6DTAU=C_S6*DPDTAU

      C_XD=TWO*(ONE+CFE12*P)*CFE12
      DXD_D=C_XD*DPD
      DXD_DD=C_XD*DPDD
      DXD_DTAU=C_XD*DPDTAU

      DXN_D=DS1D+DS2D+DS3D+DS4D+DS5D+DS6D
      DXN_DD=DS1DD+DS2DD+DS3DD+DS4DD+DS5DD+DS6DD
      DXN_DTAU=DS1DTAU+DS2DTAU+DS3DTAU+DS4DTAU+DS5DTAU+DS6DTAU

      DX_D=(XD*DXN_D-XN*DXD_D)/XD**TWO
!     DX_DD=(XD*DXN_DD-XN*DXD_DD)/XD**TWO
      DX_DD=DXN_DD/XD-X*DXD_DD/XD
      DX_DTAU=(XD*DXN_DTAU-XN*DXD_DTAU)/XD**TWO

      C_FX=ONE/(ONE+X/RKAPPA)**TWO
      DFX_D=C_FX*DX_D
      DFX_DD=C_FX*DX_DD
      DFX_DTAU=C_FX*DX_DTAU

!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
      EX_REVTPSS=0._q
      EX_REVTPSS=EX_REVTPSS+EXLDA*FX

! spin down
! IN EXD1(2*RD), TAUWD AND TAUD SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RD
      DRHO=TWO*DRD
      TAUW_RHO=TWO*TAUWD_TEMP
      TAU_RHO=TWO*TAUD

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|

      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P
      Z=TAUW_RHO/TAU_RHO
! test
      !Z=max(min(Z,0.9999999999999_q),1.E-20_q)
#ifdef Zmethod1
      IF (Z>1) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=1
      ELSE IF(Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
      ENDIF
#endif
#ifdef Zmethod2
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
         DZTILDEDZ=1._q/(1._q+Z**xorder)**(1._q/xorder) - Z**xorder/(1._q+Z**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
         Z=Z/(1._q+Z**xorder)**(1._q/xorder)
      ENDIF
#endif
#ifdef Zmethod3
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1._q/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder) - (Z/zinfinity)**xorder/(1._q+(Z/zinfinity)**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=Z/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder)
      ENDIF
#endif
! test
      Z2=Z*Z
      Z3=Z2*Z
      ZTF=0.6_q*Z
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5

!     ALPHA=(TAU_RHO-TAUW_RHO)/TAU0
      ALPHA=5._q/3._q*P*(1._q/Z-1._q)

      QB_NUM=9._q/20._q*(ALPHA-ONE)
      QB_DEN=SQRT(ONE+CFB*ALPHA*(ALPHA-ONE))
      QB=QB_NUM/QB_DEN+ TWO*P/THREE
      OPZ=ONE+Z2
      SFP=Z3/OPZ**TWO
      FQ=SQRT((ZTF**TWO+P2)/TWO)
      S1=(CX1+CFC*SFP)*P
      S2=CX2*QB**TWO
      S3=-CX3*FQ*QB
      S4=(CX1*P)**TWO/RKAPPA
      S5=TWO*CFE12*CX1*ZTF**TWO
      S6=CFE*CFMU*P2*P
      XN=S1+S2+S3+S4+S5+S6
      XD=(ONE+CFE12*P)**TWO
      X=XN/XD

!    GET THE VALUE FOR FX
      FX=ONE+RKAPPA-RKAPPA/(ONE+(X/RKAPPA))

!    NOW, DERIVATIVES COME
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q

      TAU_UNIF=THREE/10._q*(THREE*PI**TWO)**THRD2*RHO**THRD5

      DALPHAD=THRD5*(ONE/Z-ONE)*DPD-THRD5*P*DZD/Z**TWO
!     DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)+   &
!    &   DRHO**TWO/RHO**(11._q/THREE)*(10._q/9._q/(THREE*PI**TWO)**THRD2)

      DALPHADD=THRD5*(ONE/Z-ONE)*DPDD-THRD5*P*DZDD/Z**TWO
!     DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)

      DALPHADTAU=-THRD5*P*DZDTAU/Z**TWO
!     DALPHADTAU=ONE/TAU_UNIF

      DQB_NUM_D=9._q/20._q*DALPHAD
      DQB_NUM_DD=9._q/20._q*DALPHADD
      DQB_NUM_DTAU=9._q/20._q*DALPHADTAU
      C_DQB_DEN=ONE/TWO*CFB*(TWO*ALPHA-ONE)/QB_DEN
      DQB_DEN_D=C_DQB_DEN*DALPHAD
      DQB_DEN_DD=C_DQB_DEN*DALPHADD
      DQB_DEN_DTAU=C_DQB_DEN*DALPHADTAU

      DQBD=(DQB_NUM_D*QB_DEN-QB_NUM*DQB_DEN_D)/QB_DEN**TWO+THRD2*DPD
      DQBDD=(DQB_NUM_DD*QB_DEN-QB_NUM*DQB_DEN_DD)/QB_DEN**TWO+THRD2*DPDD
      DQBDTAU=(DQB_NUM_DTAU*QB_DEN-QB_NUM*DQB_DEN_DTAU)/QB_DEN**TWO+THRD2*DPDTAU

      C1_DS1=P*Z2*(THREE*CFC-CFC*Z2)/OPZ**THREE
      C2_DS1=(CX1+CFC*Z3/OPZ**TWO)
      DS1D=C1_DS1*DZD+C2_DS1*DPD
      DS1DD=C1_DS1*DZDD+C2_DS1*DPDD
      DS1DTAU=C1_DS1*DZDTAU+C2_DS1*DPDTAU

      C_DS2=CX2*TWO*QB
      DS2D=C_DS2*DQBD
      DS2DD=C_DS2*DQBDD
      DS2DTAU=C_DS2*DQBDTAU

      DFQD=ONE/TWO/FQ*(0.6_Q**TWO*Z*DZD+P*DPD)
      DFQDD=ONE/TWO/FQ*(0.6_Q**TWO*Z*DZDD+P*DPDD)
      DFQDTAU=ONE/TWO/FQ*(0.6_Q**TWO*Z*DZDTAU+P*DPDTAU)
      DS3D=-CX3*(FQ*DQBD+QB*DFQD)
      DS3DD=-CX3*(FQ*DQBDD+QB*DFQDD)
      DS3DTAU=-CX3*(FQ*DQBDTAU+QB*DFQDTAU)

      C_S4=TWO/RKAPPA*CX1**TWO*P
      DS4D=C_S4*DPD
      DS4DD=C_S4*DPDD
      DS4DTAU=C_S4*DPDTAU

      C_S5=FOUR*CFE12*CX1*0.6_q**TWO*Z
      DS5D=C_S5*DZD
      DS5DD=C_S5*DZDD
      DS5DTAU=C_S5*DZDTAU

      C_S6=THREE*CFE*CFMU*P2
      DS6D=C_S6*DPD
      DS6DD=C_S6*DPDD
      DS6DTAU=C_S6*DPDTAU

      C_XD=TWO*(ONE+CFE12*P)*CFE12
      DXD_D=C_XD*DPD
      DXD_DD=C_XD*DPDD
      DXD_DTAU=C_XD*DPDTAU

      DXN_D=DS1D+DS2D+DS3D+DS4D+DS5D+DS6D
      DXN_DD=DS1DD+DS2DD+DS3DD+DS4DD+DS5DD+DS6DD
      DXN_DTAU=DS1DTAU+DS2DTAU+DS3DTAU+DS4DTAU+DS5DTAU+DS6DTAU

      DX_D=(XD*DXN_D-XN*DXD_D)/XD**TWO
!     DX_DD=(XD*DXN_DD-XN*DXD_DD)/XD**TWO
      DX_DD=DXN_DD/XD-X*DXD_DD/XD
      DX_DTAU=(XD*DXN_DTAU-XN*DXD_DTAU)/XD**TWO

      C_FX=ONE/(ONE+X/RKAPPA)**TWO
      DFX_D=C_FX*DX_D
      DFX_DD=C_FX*DX_DD
      DFX_DTAU=C_FX*DX_DTAU

!   OUTPUT THE VXD2,VXDD2 AND AMUXD2
      VXD2=EXDLDA*FX+EXLDA*DFX_D
      VXDD2=EXLDA*DFX_DD
      AMUXD2=EXLDA*DFX_DTAU

      EX_REVTPSS=EX_REVTPSS+EXLDA*FX
      EX_REVTPSS=EX_REVTPSS/TWO
      RETURN
      END SUBROUTINE VrevTPSSx


!***********************************************************************
!
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VCD1 VCD2                    THE DERIVATIVES OF EC WRT n
!> VCDD1,VCDD2                  THE DERIVATIVES OF EC WRT |grad n|
!> AMUCD1, AMUCD2                   THE DERIVATIVES OF EC WRT TAU
!
!***********************************************************************

      SUBROUTINE VrevTPSSc(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : COE_RPKZB, CORGGA_REVTPSS
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP!$OMP DECLARE TARGET(VrevTPSSc)

! the following parameters are given by Perdew et.al.
! the coefficients for the correlatin part
      PARAMETER (CFD=2.8_q)

! other parameters
      PARAMETER (ZERO=0._q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
! test
      PARAMETER (thresh=1.e-10_q)
      PARAMETER (xorder=12._q)
      PARAMETER (zinfinity=2._q)
! test
GPACC !$ACC ROUTINE(COE_RPKZB) SEQ
GPACC !$ACC ROUTINE(CORGGA_REVTPSS) SEQ

      VCD1=0._q;VCD2=0._q;
      VCDD1=0._q;VCDD2=0._q;VCDD12=0._q;
      AMUCD1=0._q;AMUCD2=0._q;

      RT=RU+RD
      YA=DRU**2._q
      YB=DRD**2._q
      Y=DRT**2._q
!    YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2._q
      TAUW=1._q/8._q*(Y/RT)
      TAU=TAUU+TAUD

!     TAUW=MIN(TAUW,TAU)

      Z=TAUW/TAU
! test
      !Z=max(min(Z,0.9999999999999_q),1.E-20_q)
#ifdef Zmethod1
      DZTILDEDZ=ONE
      IF (Z>ONE.OR.Z<ZERO) DZTILDEDZ=ZERO
      DZD1=-Z/RT*DZTILDEDZ
      DZD2=DZD1
      DZDD12=ONE/FOUR*DRT/(RT*TAU)*DZTILDEDZ
      DZDTAU=-Z/TAU*DZTILDEDZ
      Z = MAX(MIN(ONE,Z),ZERO)
#endif
#ifdef Zmethod2
      IF (Z<ZERO) THEN
        DZTILDEDZ=ZERO
      ELSE
        DZTILDEDZ=ONE/(ONE+Z**xorder)**(ONE/xorder) - Z**xorder/(ONE+Z**xorder)**(ONE+ONE/xorder)
      ENDIF
      DZD1=-Z/RT*DZTILDEDZ
      DZD2=DZD1
      DZDD12=ONE/FOUR*DRT/(RT*TAU)*DZTILDEDZ
      DZDTAU=-Z/TAU*DZTILDEDZ
      Z=MAX(ZERO,Z)
      Z=Z/(1._q+Z**xorder)**(1._q/xorder)
#endif
#ifdef Zmethod3
      IF (Z<ZERO) THEN
        DZTILDEDZ=ZERO
      ELSE
        DZTILDEDZ=ONE/(ONE+(Z/zinfinity)**xorder)**(ONE/xorder) - (Z/zinfinity)**xorder/(ONE+(Z/zinfinity)**xorder)**(ONE+ONE/xorder)
      ENDIF
      DZD1=-Z/RT*DZTILDEDZ
      DZD2=DZD1
      DZDD12=ONE/FOUR*DRT/(RT*TAU)*DZTILDEDZ
      DZDTAU=-Z/TAU*DZTILDEDZ
      Z=MAX(ZERO,Z)
      Z=Z/(ONE+(Z/zinfinity)**xorder)**(ONE/xorder)
#endif
! test
      Z2=Z*Z
      Z3=Z2*Z

      CALL COE_RPKZB(RU,RD,DRU,DRD,DRT,CRPKZB,D_CRPKZB_D1,D_CRPKZB_D2,D_CRPKZB_DD1,D_CRPKZB_DD2,D_CRPKZB_DD12)

      CALL CORGGA_REVTPSS(RU,RD,DRU,DRD,DRT,EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD12)

! EPPGGA2_D1, EPPGGA1_D2 EPPGGA1_DD2 AND EPPGGA2_DD1 ARE ZEROS
      CALL CORGGA_REVTPSS(RU,ZERO,DRU,ZERO,DRU,EPPGGA1,EPPGGA1_D1,EPPGGA1_D2,EPPGGA1_DD1)
      EPPGGA1_D2=ZERO
      EPPGGA1_DD2=ZERO
      EPPGGA1_DD12=ZERO

      CALL CORGGA_REVTPSS(RD,ZERO,DRD,ZERO,DRD,EPPGGA2,EPPGGA2_D2,EPPGGA2_D1,EPPGGA2_DD2)
      EPPGGA2_D1=ZERO
      EPPGGA2_DD1=ZERO
      EPPGGA2_DD12=ZERO

      IF (EPPGGA1 .LT. EPPGGA)THEN
          EPPGGA1=EPPGGA
          EPPGGA1_D1=EPPGGA_D1
          EPPGGA1_D2=EPPGGA_D2
          EPPGGA1_DD1=ZERO
          EPPGGA1_DD12=EPPGGA_DD12
      ENDIF

      IF (EPPGGA2 .LT. EPPGGA)THEN
          EPPGGA2=EPPGGA
          EPPGGA2_D1=EPPGGA_D1
          EPPGGA2_D2=EPPGGA_D2
          EPPGGA2_DD2=ZERO
          EPPGGA2_DD12=EPPGGA_DD12
      ENDIF

      WEIRD=(RU*EPPGGA1+RD*EPPGGA2)/RT
      WEIRD_D1=RD*EPPGGA1/RT**TWO+RU/RT*EPPGGA1_D1-RD*EPPGGA2/RT**TWO+RD/RT*EPPGGA2_D1
      WEIRD_D2=-RU*EPPGGA1/RT**TWO+RU/RT*EPPGGA1_D2+RU*EPPGGA2/RT**TWO+RD/RT*EPPGGA2_D2
      WEIRD_DD1=RU/RT*EPPGGA1_DD1+RD/RT*EPPGGA2_DD1
      WEIRD_DD2=RU/RT*EPPGGA1_DD2+RD/RT*EPPGGA2_DD2
      WEIRD_DD12=RU/RT*EPPGGA1_DD12+RD/RT*EPPGGA2_DD12

      EPSC_REVPKZB=EPPGGA*(ONE+CRPKZB*Z2)-(ONE+CRPKZB)*Z2*WEIRD

      FP_D_EPSC_REVPKZB_D1=(ONE+CRPKZB*Z2)*EPPGGA_D1+EPPGGA*(D_CRPKZB_D1*Z2+TWO*Z*CRPKZB*DZD1)
      FP_D_EPSC_REVPKZB_D2=(ONE+CRPKZB*Z2)*EPPGGA_D2+EPPGGA*(D_CRPKZB_D2*Z2+TWO*Z*CRPKZB*DZD2)

      SD_D_EPSC_REVPKZB_D1=D_CRPKZB_D1*Z2*WEIRD+(ONE+CRPKZB)*TWO*Z*DZD1*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_D1
      SD_D_EPSC_REVPKZB_D2=D_CRPKZB_D2*Z2*WEIRD+(ONE+CRPKZB)*TWO*Z*DZD2*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_D2

      D_EPSC_REVPKZB_D1=FP_D_EPSC_REVPKZB_D1-SD_D_EPSC_REVPKZB_D1
      D_EPSC_REVPKZB_D2=FP_D_EPSC_REVPKZB_D2-SD_D_EPSC_REVPKZB_D2

      FP_D_EPSC_REVPKZB_DD1=EPPGGA*D_CRPKZB_DD1*Z2
      FP_D_EPSC_REVPKZB_DD2=EPPGGA*D_CRPKZB_DD2*Z2
      FP_D_EPSC_REVPKZB_DD12=(ONE+CRPKZB*Z2)*EPPGGA_DD12+EPPGGA*(D_CRPKZB_DD12*Z2+TWO*Z*CRPKZB*DZDD12)

      SD_D_EPSC_REVPKZB_DD1=D_CRPKZB_DD1*Z2*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_DD1
      SD_D_EPSC_REVPKZB_DD2=D_CRPKZB_DD2*Z2*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_DD2
      SD_D_EPSC_REVPKZB_DD12=D_CRPKZB_DD12*Z2*WEIRD+(ONE+CRPKZB)*TWO*Z*DZDD12*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_DD12

      D_EPSC_REVPKZB_DD1=FP_D_EPSC_REVPKZB_DD1-SD_D_EPSC_REVPKZB_DD1
      D_EPSC_REVPKZB_DD2=FP_D_EPSC_REVPKZB_DD2-SD_D_EPSC_REVPKZB_DD2
      D_EPSC_REVPKZB_DD12=FP_D_EPSC_REVPKZB_DD12-SD_D_EPSC_REVPKZB_DD12

      D_EPSC_REVPKZB_DTAU=EPPGGA*CRPKZB*TWO*Z*DZDTAU-(ONE+CRPKZB)*TWO*Z*DZDTAU*WEIRD

      VCD1=(EPSC_REVPKZB+RT*D_EPSC_REVPKZB_D1)*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_D1+THREE*Z2*EPSC_REVPKZB*DZD1)

      VCD2=(EPSC_REVPKZB+RT*D_EPSC_REVPKZB_D2)*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_D2+THREE*Z2*EPSC_REVPKZB*DZD2)

      VCDD1=RT*D_EPSC_REVPKZB_DD1*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DD1)

      VCDD2=RT*D_EPSC_REVPKZB_DD2*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DD2)

      VCDD12=RT*D_EPSC_REVPKZB_DD12*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DD12+THREE*Z2*EPSC_REVPKZB*DZDD12)

      AMUCD1=RT*D_EPSC_REVPKZB_DTAU*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DTAU+THREE*Z2*EPSC_REVPKZB*DZDTAU)
      AMUCD2=AMUCD1

      EC_REVTPSS=RT*EPSC_REVPKZB*(ONE+CFD*EPSC_REVPKZB*Z3)
      RETURN
      END SUBROUTINE VrevTPSSC

!************************ SUBROUTINE VTPSSx ****************************
!
!> calculates the first order derivatives of Ex wrt n and |grad(n)|
!> Perdew et. al. PRL (2009)
!>
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> everything in Hartree units
!>
!> ATTENTION: Every values are passed "as they are", i.e. including
!> possibly unphysical numerical errors (e.g. negative charge densities)
!> values need to be checked accordingly
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VXD1 VXD2                    THE DERIVATIVES OF EX WRT n
!> VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n|
!> AMUXD1, AMUXD2               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VTPSSx(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VTPSSx)

! the following parameters are given by Perdew et.al.
! the coefficients for the exchange part
      PARAMETER (RKAPPA=0.804_q)
      PARAMETER (CFB=0.40_q)
      PARAMETER (CFC=1.59096_q)
      PARAMETER (CFE=1.537_q)
      PARAMETER (CFMU=0.21951_q)

! other parameters
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
! test
      PARAMETER (thresh=1.e-10_q)
      PARAMETER (xorder=12._q)
      PARAMETER (zinfinity=2._q)
! test

      VXD1=0._q;VXD2=0._q;
      VXDD1=0._q;VXDD2=0._q;
      AMUXD1=0._q;AMUXD2=0._q;


      CX1=10._q/81._q
      CX2=146._q/2025._q
      CX3=73._q/405._q
      CFE12=SQRT(CFE)

! Suspect that TAUWU and TAUWD are not well described. Use TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUWU_TEMP=0.125_q*DRU**2._q/RU
!     TAUWU_TEMP=MIN(TAUWU_TEMP,TAUU)

      TAUWD_TEMP=0.125_q*DRD**2._q/RD
!     TAUWD_TEMP=MIN(TAUWD_TEMP,TAUD)

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RU
      DRHO=TWO*DRU
      TAUW_RHO=TWO*TAUWU_TEMP
      TAU_RHO=TWO*TAUU

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P
      Z=TAUW_RHO/TAU_RHO
! test
      !Z=max(min(Z,0.9999999999999_q),1.E-20_q)
#ifdef Zmethod1
      IF (Z>1) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=1
      ELSE IF(Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
      ENDIF
#endif
#ifdef Zmethod2
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1._q/(1._q+Z**xorder)**(1._q/xorder) - Z**xorder/(1._q+Z**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=Z/(1._q+Z**xorder)**(1._q/xorder)
      ENDIF
#endif
#ifdef Zmethod3
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1._q/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder) - (Z/zinfinity)**xorder/(1._q+(Z/zinfinity)**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=Z/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder)
      ENDIF
#endif
! test
      Z2=Z*Z
      Z3=Z2*Z
      ZTF=0.6_q*Z
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5

!     ALPHA=(TAU_RHO-TAUW_RHO)/TAU0
      ALPHA=5._q/3._q*P*(1._q/Z-1._q)

      QB_NUM=9._q/20._q*(ALPHA-ONE)
      QB_DEN=SQRT(ONE+CFB*ALPHA*(ALPHA-ONE))
      QB=QB_NUM/QB_DEN+ TWO*P/THREE
      OPZ=ONE+Z2
      SFP=Z2/OPZ**TWO
      FQ=SQRT((ZTF**TWO+P2)/TWO)
      S1=(CX1+CFC*SFP)*P
      S2=CX2*QB**TWO
      S3=-CX3*FQ*QB
      S4=(CX1*P)**TWO/RKAPPA
      S5=TWO*CFE12*CX1*ZTF**TWO
      S6=CFE*CFMU*P2*P
      XN=S1+S2+S3+S4+S5+S6
      XD=(ONE+CFE12*P)**TWO
      X=XN/XD

!    GET THE VALUE FOR FX
      FX=ONE+RKAPPA-RKAPPA/(ONE+(X/RKAPPA))

!    NOW, DERIVATIVES COME
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q

      TAU_UNIF=THREE/10._q*(THREE*PI**TWO)**THRD2*RHO**THRD5

      DALPHAD=THRD5*(ONE/Z-ONE)*DPD-THRD5*P*DZD/Z**TWO
!     DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)+   &
!    &   DRHO**TWO/RHO**(11._q/THREE)*(10._q/9._q/(THREE*PI**TWO)**THRD2)

      DALPHADD=THRD5*(ONE/Z-ONE)*DPDD-THRD5*P*DZDD/Z**TWO
!     DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)

      DALPHADTAU=-THRD5*P*DZDTAU/Z**TWO
!     DALPHADTAU=ONE/TAU_UNIF

      DQB_NUM_D=9._q/20._q*DALPHAD
      DQB_NUM_DD=9._q/20._q*DALPHADD
      DQB_NUM_DTAU=9._q/20._q*DALPHADTAU
      C_DQB_DEN=ONE/TWO*CFB*(TWO*ALPHA-ONE)/QB_DEN
      DQB_DEN_D=C_DQB_DEN*DALPHAD
      DQB_DEN_DD=C_DQB_DEN*DALPHADD
      DQB_DEN_DTAU=C_DQB_DEN*DALPHADTAU

      DQBD=(DQB_NUM_D*QB_DEN-QB_NUM*DQB_DEN_D)/QB_DEN**TWO+THRD2*DPD
      DQBDD=(DQB_NUM_DD*QB_DEN-QB_NUM*DQB_DEN_DD)/QB_DEN**TWO+THRD2*DPDD
      DQBDTAU=(DQB_NUM_DTAU*QB_DEN-QB_NUM*DQB_DEN_DTAU)/QB_DEN**TWO+THRD2*DPDTAU

      C1_DS1=TWO*CFC*P*Z*(ONE-Z2)/OPZ**THREE
      C2_DS1=(CX1+CFC*Z2/OPZ**TWO)
      DS1D=C1_DS1*DZD+C2_DS1*DPD
      DS1DD=C1_DS1*DZDD+C2_DS1*DPDD
      DS1DTAU=C1_DS1*DZDTAU+C2_DS1*DPDTAU

      C_DS2=CX2*TWO*QB
      DS2D=C_DS2*DQBD
      DS2DD=C_DS2*DQBDD
      DS2DTAU=C_DS2*DQBDTAU

      DFQD=ONE/TWO/FQ*(0.6_q**TWO*Z*DZD+P*DPD)
      DFQDD=ONE/TWO/FQ*(0.6_q**TWO*Z*DZDD+P*DPDD)
      DFQDTAU=ONE/TWO/FQ*(0.6_q**TWO*Z*DZDTAU+P*DPDTAU)
      DS3D=-CX3*(FQ*DQBD+QB*DFQD)
      DS3DD=-CX3*(FQ*DQBDD+QB*DFQDD)
      DS3DTAU=-CX3*(FQ*DQBDTAU+QB*DFQDTAU)

      C_S4=TWO/RKAPPA*CX1**TWO*P
      DS4D=C_S4*DPD
      DS4DD=C_S4*DPDD
      DS4DTAU=C_S4*DPDTAU

      C_S5=FOUR*CFE12*CX1*0.6_q**TWO*Z
      DS5D=C_S5*DZD
      DS5DD=C_S5*DZDD
      DS5DTAU=C_S5*DZDTAU

      C_S6=THREE*CFE*CFMU*P2
      DS6D=C_S6*DPD
      DS6DD=C_S6*DPDD
      DS6DTAU=C_S6*DPDTAU

      C_XD=TWO*(ONE+CFE12*P)*CFE12
      DXD_D=C_XD*DPD
      DXD_DD=C_XD*DPDD
      DXD_DTAU=C_XD*DPDTAU

      DXN_D=DS1D+DS2D+DS3D+DS4D+DS5D+DS6D
      DXN_DD=DS1DD+DS2DD+DS3DD+DS4DD+DS5DD+DS6DD
      DXN_DTAU=DS1DTAU+DS2DTAU+DS3DTAU+DS4DTAU+DS5DTAU+DS6DTAU

      DX_D=(XD*DXN_D-XN*DXD_D)/XD**TWO
!     DX_DD=(XD*DXN_DD-XN*DXD_DD)/XD**TWO
      DX_DD=DXN_DD/XD-X*DXD_DD/XD
      DX_DTAU=(XD*DXN_DTAU-XN*DXD_DTAU)/XD**TWO

      C_FX=ONE/(ONE+X/RKAPPA)**TWO
      DFX_D=C_FX*DX_D
      DFX_DD=C_FX*DX_DD
      DFX_DTAU=C_FX*DX_DTAU

!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
      EX_REVTPSS=0._q
      EX_REVTPSS=EX_REVTPSS+EXLDA*FX

! spin down
! IN EXD1(2*RD), TAUWD AND TAUD SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RD
      DRHO=TWO*DRD
      TAUW_RHO=TWO*TAUWD_TEMP
      TAU_RHO=TWO*TAUD

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|

      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P
      Z=TAUW_RHO/TAU_RHO
! test
      !Z=max(min(Z,0.9999999999999_q),1.E-20_q)
#ifdef Zmethod1
      IF (Z>1) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=1
      ELSE IF(Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
      ENDIF
#endif
#ifdef Zmethod2
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
         DZTILDEDZ=1._q/(1._q+Z**xorder)**(1._q/xorder) - Z**xorder/(1._q+Z**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
         Z=Z/(1._q+Z**xorder)**(1._q/xorder)
      ENDIF
#endif
#ifdef Zmethod3
      IF (Z<0) THEN
        DZTILDEDZ=0
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=0
      ELSE
        DZTILDEDZ=1._q/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder) - (Z/zinfinity)**xorder/(1._q+(Z/zinfinity)**xorder)**(1._q+1._q/xorder)
        DZD=-Z/RHO*DZTILDEDZ
        DZDD=DRHO/(FOUR*RHO*TAU_RHO)*DZTILDEDZ
        DZDTAU=-Z/TAU_RHO*DZTILDEDZ
        Z=Z/(1._q+(Z/zinfinity)**xorder)**(1._q/xorder)
      ENDIF
#endif
! test
      Z2=Z*Z
      Z3=Z2*Z
      ZTF=0.6_q*Z
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5

!     ALPHA=(TAU_RHO-TAUW_RHO)/TAU0
      ALPHA=5._q/3._q*P*(1._q/Z-1._q)

      QB_NUM=9._q/20._q*(ALPHA-ONE)
      QB_DEN=SQRT(ONE+CFB*ALPHA*(ALPHA-ONE))
      QB=QB_NUM/QB_DEN+ TWO*P/THREE
      OPZ=ONE+Z2
      SFP=Z2/OPZ**TWO
      FQ=SQRT((ZTF**TWO+P2)/TWO)
      S1=(CX1+CFC*SFP)*P
      S2=CX2*QB**TWO
      S3=-CX3*FQ*QB
      S4=(CX1*P)**TWO/RKAPPA
      S5=TWO*CFE12*CX1*ZTF**TWO
      S6=CFE*CFMU*P2*P
      XN=S1+S2+S3+S4+S5+S6
      XD=(ONE+CFE12*P)**TWO
      X=XN/XD

!    GET THE VALUE FOR FX
      FX=ONE+RKAPPA-RKAPPA/(ONE+(X/RKAPPA))

!    NOW, DERIVATIVES COME
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q

      TAU_UNIF=THREE/10._q*(THREE*PI**TWO)**THRD2*RHO**THRD5

      DALPHAD=THRD5*(ONE/Z-ONE)*DPD-THRD5*P*DZD/Z**TWO
!     DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)+   &
!    &   DRHO**TWO/RHO**(11._q/THREE)*(10._q/9._q/(THREE*PI**TWO)**THRD2)

      DALPHADD=THRD5*(ONE/Z-ONE)*DPDD-THRD5*P*DZDD/Z**TWO
!     DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)

      DALPHADTAU=-THRD5*P*DZDTAU/Z**TWO
!     DALPHADTAU=ONE/TAU_UNIF

      DQB_NUM_D=9._q/20._q*DALPHAD
      DQB_NUM_DD=9._q/20._q*DALPHADD
      DQB_NUM_DTAU=9._q/20._q*DALPHADTAU
      C_DQB_DEN=ONE/TWO*CFB*(TWO*ALPHA-ONE)/QB_DEN
      DQB_DEN_D=C_DQB_DEN*DALPHAD
      DQB_DEN_DD=C_DQB_DEN*DALPHADD
      DQB_DEN_DTAU=C_DQB_DEN*DALPHADTAU

      DQBD=(DQB_NUM_D*QB_DEN-QB_NUM*DQB_DEN_D)/QB_DEN**TWO+THRD2*DPD
      DQBDD=(DQB_NUM_DD*QB_DEN-QB_NUM*DQB_DEN_DD)/QB_DEN**TWO+THRD2*DPDD
      DQBDTAU=(DQB_NUM_DTAU*QB_DEN-QB_NUM*DQB_DEN_DTAU)/QB_DEN**TWO+THRD2*DPDTAU

      C1_DS1=TWO*CFC*P*Z*(ONE-Z2)/OPZ**THREE
      C2_DS1=(CX1+CFC*Z2/OPZ**TWO)
      DS1D=C1_DS1*DZD+C2_DS1*DPD
      DS1DD=C1_DS1*DZDD+C2_DS1*DPDD
      DS1DTAU=C1_DS1*DZDTAU+C2_DS1*DPDTAU

      C_DS2=CX2*TWO*QB
      DS2D=C_DS2*DQBD
      DS2DD=C_DS2*DQBDD
      DS2DTAU=C_DS2*DQBDTAU

      DFQD=ONE/TWO/FQ*(0.6_Q**TWO*Z*DZD+P*DPD)
      DFQDD=ONE/TWO/FQ*(0.6_Q**TWO*Z*DZDD+P*DPDD)
      DFQDTAU=ONE/TWO/FQ*(0.6_Q**TWO*Z*DZDTAU+P*DPDTAU)
      DS3D=-CX3*(FQ*DQBD+QB*DFQD)
      DS3DD=-CX3*(FQ*DQBDD+QB*DFQDD)
      DS3DTAU=-CX3*(FQ*DQBDTAU+QB*DFQDTAU)

      C_S4=TWO/RKAPPA*CX1**TWO*P
      DS4D=C_S4*DPD
      DS4DD=C_S4*DPDD
      DS4DTAU=C_S4*DPDTAU

      C_S5=FOUR*CFE12*CX1*0.6_q**TWO*Z
      DS5D=C_S5*DZD
      DS5DD=C_S5*DZDD
      DS5DTAU=C_S5*DZDTAU

      C_S6=THREE*CFE*CFMU*P2
      DS6D=C_S6*DPD
      DS6DD=C_S6*DPDD
      DS6DTAU=C_S6*DPDTAU

      C_XD=TWO*(ONE+CFE12*P)*CFE12
      DXD_D=C_XD*DPD
      DXD_DD=C_XD*DPDD
      DXD_DTAU=C_XD*DPDTAU

      DXN_D=DS1D+DS2D+DS3D+DS4D+DS5D+DS6D
      DXN_DD=DS1DD+DS2DD+DS3DD+DS4DD+DS5DD+DS6DD
      DXN_DTAU=DS1DTAU+DS2DTAU+DS3DTAU+DS4DTAU+DS5DTAU+DS6DTAU

      DX_D=(XD*DXN_D-XN*DXD_D)/XD**TWO
!     DX_DD=(XD*DXN_DD-XN*DXD_DD)/XD**TWO
      DX_DD=DXN_DD/XD-X*DXD_DD/XD
      DX_DTAU=(XD*DXN_DTAU-XN*DXD_DTAU)/XD**TWO

      C_FX=ONE/(ONE+X/RKAPPA)**TWO
      DFX_D=C_FX*DX_D
      DFX_DD=C_FX*DX_DD
      DFX_DTAU=C_FX*DX_DTAU

!   OUTPUT THE VXD2,VXDD2 AND AMUXD2
      VXD2=EXDLDA*FX+EXLDA*DFX_D
      VXDD2=EXLDA*DFX_DD
      AMUXD2=EXLDA*DFX_DTAU

      EX_REVTPSS=EX_REVTPSS+EXLDA*FX
      EX_REVTPSS=EX_REVTPSS/TWO
      RETURN
      END SUBROUTINE VTPSSx


!***********************************************************************
!
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VCD1 VCD2                    THE DERIVATIVES OF EC WRT n
!> VCDD1,VCDD2                  THE DERIVATIVES OF EC WRT |grad n|
!> AMUCD1, AMUCD2                   THE DERIVATIVES OF EC WRT TAU
!
!***********************************************************************

      SUBROUTINE VTPSSc(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,VCDD12,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : COE_RPKZB_TPSS, CORGGA_TPSS
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VTPSSc)

! the following parameters are given by Perdew et.al.
! the coefficients for the correlatin part
      PARAMETER (CFD=2.8_q)

! other parameters
      PARAMETER (ZERO=0._q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
! test
      PARAMETER (thresh=1.e-10_q)
      PARAMETER (xorder=12._q)
      PARAMETER (zinfinity=2._q)
! test
GPACC !$ACC ROUTINE(COE_RPKZB_TPSS) SEQ
GPACC !$ACC ROUTINE(CORGGA_TPSS) SEQ

      VCD1=0._q;VCD2=0._q;
      VCDD1=0._q;VCDD2=0._q;VCDD12=0._q
      AMUCD1=0._q;AMUCD2=0._q;

      RT=RU+RD
      YA=DRU**2._q
      YB=DRD**2._q
      Y=DRT**2._q
!    YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2._q
      TAUW=1._q/8._q*(Y/RT)
      TAU=TAUU+TAUD

      Z=TAUW/TAU
! test
      !Z=max(min(Z,0.9999999999999_q),1.E-20_q)


#ifdef Zmethod1
      DZTILDEDZ=ONE
      IF (Z>ONE.OR.Z<ZERO) DZTILDEDZ=ZERO
      DZD1=-Z/RT*DZTILDEDZ
      DZD2=DZD1
      DZDD12=ONE/FOUR*DRT/(RT*TAU)*DZTILDEDZ
      DZDTAU=-Z/TAU*DZTILDEDZ
      Z = MAX(MIN(ONE,Z),ZERO)
#endif
#ifdef Zmethod2
      IF (Z<ZERO) THEN
        DZTILDEDZ=ZERO
      ELSE
        DZTILDEDZ=ONE/(ONE+Z**xorder)**(ONE/xorder) - Z**xorder/(ONE+Z**xorder)**(ONE+ONE/xorder)
      ENDIF
      DZD1=-Z/RT*DZTILDEDZ
      DZD2=DZD1
      DZDD12=ONE/FOUR*DRT/(RT*TAU)*DZTILDEDZ
      DZDTAU=-Z/TAU*DZTILDEDZ
      Z=MAX(ZERO,Z)
      Z=Z/(1._q+Z**xorder)**(1._q/xorder)
#endif
#ifdef Zmethod3
      IF (Z<ZERO) THEN
        DZTILDEDZ=ZERO
      ELSE
        DZTILDEDZ=ONE/(ONE+(Z/zinfinity)**xorder)**(ONE/xorder) - (Z/zinfinity)**xorder/(ONE+(Z/zinfinity)**xorder)**(ONE+ONE/xorder)
      ENDIF
      DZD1=-Z/RT*DZTILDEDZ
      DZD2=DZD1
      DZDD12=ONE/FOUR*DRT/(RT*TAU)*DZTILDEDZ
      DZDTAU=-Z/TAU*DZTILDEDZ
      Z=MAX(ZERO,Z)
      Z=Z/(ONE+(Z/zinfinity)**xorder)**(ONE/xorder)
#endif
! test
      Z2=Z*Z
      Z3=Z2*Z

      CALL COE_RPKZB_TPSS(RU,RD,DRU,DRD,DRT,CRPKZB,D_CRPKZB_D1,D_CRPKZB_D2,D_CRPKZB_DD1,D_CRPKZB_DD2,D_CRPKZB_DD12)

      CALL CORGGA_TPSS(RU,RD,DRU,DRD,DRT,EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD12)

! EPPGGA2_D1, EPPGGA1_D2 EPPGGA1_DD2 AND EPPGGA2_DD1 ARE ZEROS
      CALL CORGGA_TPSS(RU,ZERO,DRU,ZERO,DRU,EPPGGA1,EPPGGA1_D1,EPPGGA1_D2,EPPGGA1_DD1)
      EPPGGA1_D2=ZERO
      EPPGGA1_DD2=ZERO
      EPPGGA1_DD12=ZERO

      CALL CORGGA_TPSS(RD,ZERO,DRD,ZERO,DRD,EPPGGA2,EPPGGA2_D2,EPPGGA2_D1,EPPGGA2_DD2)
      EPPGGA2_D1=ZERO
      EPPGGA2_DD1=ZERO
      EPPGGA2_DD12=ZERO

      IF (EPPGGA1 .LT. EPPGGA)THEN
          EPPGGA1=EPPGGA
          EPPGGA1_D1=EPPGGA_D1
          EPPGGA1_D2=EPPGGA_D2
          EPPGGA1_DD1=ZERO
          EPPGGA1_DD12=EPPGGA_DD12
      ENDIF

      IF (EPPGGA2 .LT. EPPGGA)THEN
          EPPGGA2=EPPGGA
          EPPGGA2_D1=EPPGGA_D1
          EPPGGA2_D2=EPPGGA_D2
          EPPGGA2_DD2=ZERO
          EPPGGA2_DD12=EPPGGA_DD12
      ENDIF

      WEIRD=(RU*EPPGGA1+RD*EPPGGA2)/RT
      WEIRD_D1=RD*EPPGGA1/RT**TWO+RU/RT*EPPGGA1_D1-RD*EPPGGA2/RT**TWO+RD/RT*EPPGGA2_D1
      WEIRD_D2=-RU*EPPGGA1/RT**TWO+RU/RT*EPPGGA1_D2+RU*EPPGGA2/RT**TWO+RD/RT*EPPGGA2_D2
      WEIRD_DD1=RU/RT*EPPGGA1_DD1+RD/RT*EPPGGA2_DD1
      WEIRD_DD2=RU/RT*EPPGGA1_DD2+RD/RT*EPPGGA2_DD2
      WEIRD_DD12=RU/RT*EPPGGA1_DD12+RD/RT*EPPGGA2_DD12

      EPSC_REVPKZB=EPPGGA*(ONE+CRPKZB*Z2)-(ONE+CRPKZB)*Z2*WEIRD

      FP_D_EPSC_REVPKZB_D1=(ONE+CRPKZB*Z2)*EPPGGA_D1+EPPGGA*(D_CRPKZB_D1*Z2+TWO*Z*CRPKZB*DZD1)
      FP_D_EPSC_REVPKZB_D2=(ONE+CRPKZB*Z2)*EPPGGA_D2+EPPGGA*(D_CRPKZB_D2*Z2+TWO*Z*CRPKZB*DZD2)

      SD_D_EPSC_REVPKZB_D1=D_CRPKZB_D1*Z2*WEIRD+(ONE+CRPKZB)*TWO*Z*DZD1*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_D1
      SD_D_EPSC_REVPKZB_D2=D_CRPKZB_D2*Z2*WEIRD+(ONE+CRPKZB)*TWO*Z*DZD2*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_D2

      D_EPSC_REVPKZB_D1=FP_D_EPSC_REVPKZB_D1-SD_D_EPSC_REVPKZB_D1
      D_EPSC_REVPKZB_D2=FP_D_EPSC_REVPKZB_D2-SD_D_EPSC_REVPKZB_D2

      FP_D_EPSC_REVPKZB_DD1=EPPGGA*D_CRPKZB_DD1*Z2
      FP_D_EPSC_REVPKZB_DD2=EPPGGA*D_CRPKZB_DD2*Z2
      FP_D_EPSC_REVPKZB_DD12=(ONE+CRPKZB*Z2)*EPPGGA_DD12+EPPGGA*(D_CRPKZB_DD12*Z2+TWO*Z*CRPKZB*DZDD12)

      SD_D_EPSC_REVPKZB_DD1=D_CRPKZB_DD1*Z2*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_DD1
      SD_D_EPSC_REVPKZB_DD2=D_CRPKZB_DD2*Z2*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_DD2
      SD_D_EPSC_REVPKZB_DD12=D_CRPKZB_DD12*Z2*WEIRD+(ONE+CRPKZB)*TWO*Z*DZDD12*WEIRD+(ONE+CRPKZB)*Z2*WEIRD_DD12

      D_EPSC_REVPKZB_DD1=FP_D_EPSC_REVPKZB_DD1-SD_D_EPSC_REVPKZB_DD1
      D_EPSC_REVPKZB_DD2=FP_D_EPSC_REVPKZB_DD2-SD_D_EPSC_REVPKZB_DD2
      D_EPSC_REVPKZB_DD12=FP_D_EPSC_REVPKZB_DD12-SD_D_EPSC_REVPKZB_DD12

      D_EPSC_REVPKZB_DTAU=EPPGGA*CRPKZB*TWO*Z*DZDTAU-(ONE+CRPKZB)*TWO*Z*DZDTAU*WEIRD

      VCD1=(EPSC_REVPKZB+RT*D_EPSC_REVPKZB_D1)*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_D1+THREE*Z2*EPSC_REVPKZB*DZD1)

      VCD2=(EPSC_REVPKZB+RT*D_EPSC_REVPKZB_D2)*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_D2+THREE*Z2*EPSC_REVPKZB*DZD2)

      VCDD1=RT*D_EPSC_REVPKZB_DD1*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DD1)

      VCDD2=RT*D_EPSC_REVPKZB_DD2*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DD2)

      VCDD12=RT*D_EPSC_REVPKZB_DD12*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DD12+THREE*Z2*EPSC_REVPKZB*DZDD12)

      AMUCD1=RT*D_EPSC_REVPKZB_DTAU*(ONE+CFD*EPSC_REVPKZB*Z3) &
          & +RT*EPSC_REVPKZB*CFD*(Z3*D_EPSC_REVPKZB_DTAU+THREE*Z2*EPSC_REVPKZB*DZDTAU)
      AMUCD2=AMUCD1

      EC_REVTPSS=RT*EPSC_REVPKZB*(ONE+CFD*EPSC_REVPKZB*Z3)

      RETURN
      END SUBROUTINE VTPSSc

!************************ SUBROUTINE M06 family ****************************
!
!> calculates the first order derivatives of Ex wrt n, tau, and |grad(n)|
!> Zhao & Truhlar TCA (2008)
!>
!> Written by Yan Zhao 09/22/2011
!>
!> everything in Hartree units
!>
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VXD1 VXD2                    THE DERIVATIVES OF EX WRT n
!> VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n|
!> AMUXD1, AMUXD2               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VM06x(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ex_M06,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,IJZY)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VM06x)
      PARAMETER (TWO=2._q)
GPACC !$ACC ROUTINE(M6X) SEQ
GPACC !$ACC ROUTINE(VS98X) SEQ

      STAUU = TAUU*TWO
      STAUD = TAUD*TWO
      Ex_M06 = 0._q
      EM6U  = 0._q
      EVSXU = 0._q
      EM6D  = 0._q
      EVSXD = 0._q
      VXD1=0._q;VXD2=0._q;
      VXDD1=0._q;VXDD2=0._q;
      AMUXD1=0._q;AMUXD2=0._q;

      CALL M6X(EM6U,RU,DRU,STAUU,EM6D1,EM6DD1,EM6AMUXD1,IJZY)
      CALL VS98X(EVSXU,RU,DRU,STAUU,VSD1,VSDD1,VSAMUXD1,IJZY+1)

      CALL M6X(EM6D,RD,DRD,STAUD,EM6D2,EM6DD2,EM6AMUXD2,IJZY)
      CALL VS98X(EVSXD,RD,DRD,STAUD,VSD2,VSDD2,VSAMUXD2,IJZY+1)

      Ex_M06 = EM6U + EVSXU + EM6D + EVSXD
      VXD1 = EM6D1 + VSD1
      VXD2 = EM6D2 + VSD2
!  EM6DD and VSDD are  DERIVATIVES OF EX WRT |grad n|^2
!  so We need to multiply them by 2*|grad n|
      VXDD1 = (EM6DD1 + VSDD1)*TWO*DRU
      VXDD2 = (EM6DD2 + VSDD2)*TWO*DRD
! EM6AMUXD and VSAMUXD are  DERIVATIVES OF EX WRT 2Tau
! so We need to multiply them by 2
      AMUXD1 = (EM6AMUXD1 + VSAMUXD1)*TWO
      AMUXD2 = (EM6AMUXD2 + VSAMUXD2)*TWO

      RETURN
      END SUBROUTINE VM06x

!************************ SUBROUTINE M06 family ****************************
!
!> calculates the first order derivatives of Ex wrt n, tau, and |grad(n)|
!> Zhao & Truhlar TCA (2008)
!>
!> Written by Yan Zhao 09/22/2011
!>
!> everything in Hartree units
!>
!>
!> R                     density up,down
!> DR                    abs. val. gradient of density up/down
!> TAU                   kinetic energy density up/down
!> M6D                   THE DERIVATIVES OF EX WRT n
!> M6DD                  THE DERIVATIVES OF EX WRT |grad n|
!> M6AMUXD               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE M6X(EM6,R,DR,TAU,EM6D,EM6DD,EM6AMUXD,IJZY)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(M6X)

! the coefficients
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (FIVE=5._q)
      PARAMETER (SIX=6._q)
      PARAMETER (SEVEN=7._q)
      PARAMETER (EIGHT=8._q)
      PARAMETER (NINE=9._q)

!      PI    = FOUR*ATAN(ONE)
      F1O3  = ONE / THREE
      F1O4  = ONE / FOUR
      F1O8  = ONE / EIGHT
      F2O3  = TWO / THREE
      F3O2  = THREE / TWO
      F4O3  = FOUR / THREE
      F4O9  = FOUR / NINE
      F3O5  = THREE / FIVE
      F5O3  = FIVE / THREE
      F5O2  = FIVE / TWO
      F7O3  = SEVEN / THREE
      F8O27 = EIGHT / ( THREE * NINE )
      F12   = TWO * SIX
      F18   = TWO * NINE
      F10 = TWO*FIVE
      F11 = 11._q
      DTol = 1.0e-10_q


!
!     PBE FUNCTIONAL'S PARAMETERS.
!
      C1     = 3.36116e-03_q
      C2     = 4.49267e-03_q
!     PARAMETERS FOR M06-L
      IF (IJZY.EQ.1) THEN
        AT0=    3.987756e-01_q
        AT1=    2.548219e-01_q
        AT2=    3.923994e-01_q
        AT3=    -2.103655e+00_q
        AT4=    -6.302147e+00_q
        AT5=    1.097615e+01_q
        AT6=    3.097273e+01_q
        AT7=    -2.318489e+01_q
        AT8=    -5.673480e+01_q
        AT9=    2.160364e+01_q
        AT10=   3.421814e+01_q
        AT11=   -9.049762e+00_q
       ELSEIF (IJZY.EQ.2) THEN
!     PARAMETERS FOR M06-HF
        AT0=    1.179732e-01_q
        AT1=    -1.066708e+00_q
        AT2=    -1.462405e-01_q
        AT3=    7.481848e+00_q
        AT4=    3.776679e+00_q
        AT5=    -4.436118e+01_q
        AT6=    -1.830962e+01_q
        AT7=    1.003903e+02_q
        AT8=    3.864360e+01_q
        AT9=    -9.806018e+01_q
        AT10=   -2.557716e+01_q
        AT11=   3.590404e+01_q
       ELSEIF (IJZY.EQ.3) THEN
!     PARAMETERS FOR M06
        AT0=    5.877943e-01_q
        AT1=    -1.371776e-01_q
        AT2=    2.682367e-01_q
        AT3=    -2.515898e+00_q
        AT4=    -2.978892e+00_q
        AT5=    8.710679e+00_q
        AT6=    1.688195e+01_q
        AT7=    -4.489724e+00_q
        AT8=    -3.299983e+01_q
        AT9=    -1.449050e+01_q
        AT10=   2.043747e+01_q
        AT11=   1.256504e+01_q
       ELSEIF (IJZY.EQ.4) THEN
!     PARAMETERS FOR M06-2X
        AT0=    4.600000e-01_q
        AT1=    -2.206052e-01_q
        AT2=    -9.431788e-02_q
        AT3=    2.164494e+00_q
        AT4=    -2.556466e+00_q
        AT5=    -1.422133e+01_q
        AT6=    1.555044e+01_q
        AT7=    3.598078e+01_q
        AT8=    -2.722754e+01_q
        AT9=    -3.924093e+01_q
        AT10=   1.522808e+01_q
        AT11=   1.522227e+01_q
      ENDIF


      IF ((R .GT. DTol) .AND. (Tau .GT. DTol)) THEN
        Ax = -F3o2*(F4o3*PI)**(-F1o3)
        RHOO = R
        RHO43 = RHOO**F4O3
        RRHO = ONE/RHOO       ! RECIPROCAL OF RHO
        RHO13 = RHO43*RRHO
        RHO53 = RHOO**F5O3
        TAUN = TAU
        TAUUEG=F3O5*((SIX*PI*PI)**F2O3)*RHO53
        TSIG =TAUUEG/TAUN

        WSIG =(TSIG-ONE)/(TSIG+ONE)
        FSIG=(AT0 + WSIG*(AT1 + WSIG*(AT2 + WSIG*(AT3 + WSIG*( &
     &            AT4 + WSIG*(AT5 + WSIG*(AT6 + WSIG*(AT7 + WSIG*( &
     &            AT8 + WSIG*(AT9 + WSIG*(AT10+WSIG*AT11)))))))))))
        Gamma = DR
        X = GAMMA/RHO43
        X2 = X*X
        EN = C1*X2
        ED = ONE + C2*X2
        E  = -EN/ED
        EM6=(AX+E)*FSIG*RHO43

        DEN   = TWO*C1*X
        DED   = TWO*C2*X
        DE    = -(DEN*ED-EN*DED)/(ED*ED)
        DFDW=(AT1 + WSIG*(TWO  *AT2 + WSIG*(THREE*AT3 + WSIG*( &
     &           FOUR *AT4 + WSIG*(FIVE *AT5 + WSIG*(SIX  *AT6 + WSIG*(&
     &           SEVEN*AT7 + WSIG*(EIGHT*AT8 + WSIG*(NINE *AT9 + WSIG*(&
     &           F10 *AT10+ WSIG*F11*AT11))))))))))
         DWDT = TWO/((ONE + TSIG)**2)
         DTDR = ((SIX*PI*PI)**F2O3)*(RHOO**F2O3)/TAUN
         DTDTAU = -TAUUEG/TAUN**2
         DGGADR = F4O3*RHO13*(AX+(E-X*DE))
         DFDR = DFDW*DWDT*DTDR
         DFDTAU=DFDW*DWDT*DTDTAU
         DGGADG =DE/(TWO*GAMMA)
!       DGGADG =(DE/(TWO*GAMMA))

!        DF/DRHOA
         EM6D = DGGADR*FSIG + (AX+E)*RHO43*DFDR
!        DF/DGAMMAAA
         EM6DD = DGGADG*FSIG
!        DF/DTAUA
         EM6AMUXD = RHO43*(AX+E)*DFDTAU
      ELSE
         EM6= 0._q
         EM6D = 0._q
         EM6DD = 0._q
         EM6AMUXD = 0._q
      ENDIF
      RETURN
      END SUBROUTINE M6X

!************************ SUBROUTINE VS98 family ****************************
!
!> calculates the first order derivatives of Ex wrt n, tau, and |grad(n)|
!> Zhao & Truhlar TCA (2008)
!>
!> Written by Yan Zhao 09/22/2011
!>
!> everything in Hartree units
!>
!>
!> R                     density up,down
!> DR                    abs. val. gradient of density up/down
!> TAU                   kinetic energy density up/down
!> VSD                   THE DERIVATIVES OF EX WRT n
!> VSDD                  THE DERIVATIVES OF EX WRT |grad n|
!> VSAMUXD               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VS98X(EVS,R,DR,TAU,VSD,VSDD,VSAMUXD,IJZY)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : GVT4
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(VS98X)

      REAL(q) KX,NINE

!      PARAMETER( PI = 3.1415926535897932384626433832795_q )

      PARAMETER (CF = 9.115599720_q)
      PARAMETER (AXLSDA = -0.9305257363491_q)
      PARAMETER (GG  = 0.00186726_q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (FIVE=5._q)
      PARAMETER (SIX=6._q)
      PARAMETER (SEVEN=7._q)
      PARAMETER (EIGHT=8._q)
      PARAMETER (NINE=9._q)
      PARAMETER (F10=10._q)
      PARAMETER (F11=11._q)
GPACC !$ACC ROUTINE(GVT4) SEQ

      DTol = 1.0e-10_q
!      IF ((R .LT. DTol) .OR. (Tau .LT. DTol)) Return
      F13 = ONE/THREE
      F43 = FOUR/THREE
      F53 = FIVE/THREE
      F83 = EIGHT/THREE
      F113 = F11/THREE
      F3O5 = THREE/FIVE
      F2O3 = TWO/THREE

      IF (IJZY.EQ.1) THEN

!     PARAMETERS FOR VS98

        R1=  -9.800683E-01_q
        R2=  -3.556788E-03_q
        R3=   6.250326E-03_q
        R4=  -2.354518E-05_q
        R5=  -1.282732E-04_q
        R6=   3.574822E-04_q
      ELSEIF (IJZY.EQ.2) THEN

!     PARAMETERS FOR M06-L

        R1 =   6.012244E-01_q*AXLSDA
        R2 =   4.748822E-03_q*AXLSDA
        R3 =  -8.635108E-03_q*AXLSDA
        R4 =  -9.308062E-06_q*AXLSDA
        R5 =   4.482811E-05_q*AXLSDA
        R6 =   0.000000E+00_q
      ELSEIF (IJZY.EQ.3) THEN

!     PARAMETERS FOR M06-HF

        R1 =   -1.179732E-01_q*AXLSDA
        R2 =   -2.500000E-03_q*AXLSDA
        R3 =   -1.180065E-02_q*AXLSDA
        R4 =   0.000000E+00_q
        R5 =   0.000000E+00_q
        R6 =   0.000000E+00_q
      ELSEIF (IJZY.EQ.4) THEN

!     PARAMETERS FOR M06

        R1 =   1.422057E-01_q*AXLSDA
        R2 =   7.370319E-04_q*AXLSDA
        R3 =   -1.601373E-02_q*AXLSDA
        R4 =   0.000000E+00_q
        R5 =   0.000000E+00_q
        R6 =   0.000000E+00_q
      ENDIF

      IF ((R .GT. DTol) .AND. (Tau .GT. DTol)) THEN
        RHOO = R
        RHO43 = RHOO**F43
        RRHO = ONE/RHOO       ! RECIPROCAL OF RHO
        RHO13 = RHO43*RRHO
        RHO53 = RHOO**F53
        TAUN =  TAU
        TAUUEG=F3O5*((SIX*PI*PI)**F2O3)*RHO53
        TSIG = TAUUEG/TAUN

        RHO83 = RHO53*RHOO
        GAMMA = DR*DR
        X = GAMMA/RHO83
        DXDR = -F83*X*RRHO
        DXDG = ONE/RHO83
        Z = TAUN/RHO53 - CF
        DZDR = -F53 * TAUN/RHO83
        DZDT = ONE/RHO53
        KX = ONE + GG*X + GG*Z
        XK = X/KX
        ZK = Z/KX
        CALL GVT4(GX,DGDX,DGDZ,XK,ZK,KX,GG,R1,R2,R3,R4,R5,R6)
        EVS =  RHO43*GX
        VSD =  F43*RHO13*GX + RHO43*(DGDX*DXDR + DGDZ*DZDR)
        VSDD =  RHO43*(DGDX*DXDG)
        VSAMUXD =  RHO43*(DGDZ*DZDT)
      ELSE
        EVS = 0._q
        VSD = 0._q
        VSDD = 0._q
        VSAMUXD = 0._q
      ENDIF

      RETURN
      END SUBROUTINE VS98X

!************************ SUBROUTINE M06 family ****************************
!
!> calculates the first order derivatives of Ex wrt n, tau, and |grad(n)|
!> Zhao & Truhlar TCA (2008)
!>
!> Written by Yan Zhao 09/22/2011
!>
!> everything in Hartree units
!>
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VXD1 VXD2                    THE DERIVATIVES OF EX WRT n
!> VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n|
!> AMUXD1, AMUXD2               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VM06c(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_M06,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2,IJZY)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : CORPBE
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VM06c)
      PARAMETER (F1=1._q)
      PARAMETER (F2=2._q)
      PARAMETER (F3=3._q)
      PARAMETER (F4=4._q)
      PARAMETER (COpp=0.0031_q)
GPACC !$ACC ROUTINE(VS98C) SEQ
GPACC !$ACC ROUTINE(M06CSS) SEQ
GPACC !$ACC ROUTINE(CORPBE) SEQ
      IF (IJZY.EQ.1) THEN
!     PARAMETERS FOR M06-L CORRELATION
         SOPP0= 6.042374e-01_q
         SOPP1= 1.776783e+02_q
         SOPP2= -2.513252e+02_q
         SOPP3= 7.635173e+01_q
         SOPP4= -1.255699e+01_q
      ELSEIF (IJZY.EQ.2) THEN
!     PARAMETERS FOR M06-HF CORRELATION
         SOPP0= 1.674634e+00_q
         SOPP1= 5.732017e+01_q
         SOPP2= 5.955416e+01_q
         SOPP3= -2.311007e+02_q
         SOPP4= 1.255199e+02_q
      ELSEIF (IJZY.EQ.3) THEN
!     PARAMETERS FOR M06 CORRELATION
         SOPP0= 3.741539e+00_q
         SOPP1= 2.187098e+02_q
         SOPP2= -4.531252e+02_q
         SOPP3= 2.936479e+02_q
         SOPP4= -6.287470e+01_q
      ELSEIF (IJZY.EQ.4) THEN
!     PARAMETERS FOR M06-2X CORRELATION
         SOPP0= 8.833596e-01_q
         SOPP1= 3.357972e+01_q
         SOPP2= -7.043548e+01_q
         SOPP3= 4.978271e+01_q
         SOPP4= -1.852891e+01_q
      ENDIF
      DTol = 1.0e-8_q
      Ec_M06=0._q
      VCD1=0._q;VCD2=0._q;
      VCDD1=0._q;VCDD2=0._q;
      AMUCD1=0._q;AMUCD2=0._q;


      CALL VS98C(DTol, &
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_VS,VSD1,VSDD1,VSD2,VSDD2,VSAMUCD1,VSAMUCD2,IJZY+1)

!      PI = F4*ATAN(F1)
      F6=6.0D0
      F43 = F4 / F3
      PI34 = F3 / (F4*PI)
      F13 = F1 / F3

!     parallel spin case UP
        PA = RU
      IF (RU .GT. DTol .AND. TAUU .GT. DTol) THEN
        GAA = DRU*DRU
        TAUA = F2*TAUU
        CALL M06CSS(DTol,PA,GAA,TAUA,FA,FPA,FGA,FTA,EUA, &
     &                CHIA,EUPA,CHIAP,CHIAG,IJZY)
      ELSE
        FA = 0._q
        FPA = 0._q
        FGA = 0._q
        FTA = 0._q
        EUA = 0._q
        CHIA = 0._q
        EUPA = 0._q
        CHIAP = 0._q
        CHIAG = 0._q
      ENDIF
!     parallel spin case DOWN
      PB = RD
      IF (RD .GT. DTol .AND. TAUD .GT. DTol) THEN
      GBB = DRD*DRD
      TAUB= F2*TAUD
      CALL M06CSS(DTol,PB,GBB,TAUB,FB,FPB,FGB,FTB,EUB, &
     &                CHIB,EUPB,CHIBP,CHIBG,IJZY)

       ELSE
        FB = 0._q
        FPB = 0._q
        FGB = 0._q
        FTB = 0._q
        EUB = 0._q
        CHIB = 0._q
        EUPB = 0._q
        CHIBP = 0._q
        CHIBG = 0._q
      ENDIF

!    antiparallel spin case
        P = PA + PB
!      IF (PB.GT.DTol.AND.PA.GT.DTol.AND.TAUA.GT.DTol.AND.TAUB.GT.DTol) THEN
       IF (PB.GT.DTol.AND.PA.GT.DTol) THEN
        RS = (PI34/P) ** F13
        RSP = -RS/(F3*P)
        ZETA = (PA-PB)/P
        DZDA = (F1-ZETA)/P
        DZDB = (-F1-ZETA)/P
!        CALL LSDAC(RS,ZETA,POTLC,DLDS,DLDZ,D2LDSS,D2LDSZ, &
!     &      D2LDZZ)

! Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
       CALL CORPBE(RS,ZETA,PotLC,VCUPLDA,VCDNLDA,unknown,unknown, &
         &        unknown,unknown,unknown,unknown,unknown,.FALSE.)
! OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
!      EPPGGA=PotLC
!      EPPGGA_D1=(VCUPLDA-PotLC)/RT
!      EPPGGA_D2=(VCDNLDA-PotLC)/RT
       dLdz=(VCUPLDA-VCDNLDA)/F2
       dLdS=(PotLC-ZETA*dLdZ-(VCUPLDA+VCDNLDA)/F2)*F3/RS

        EUEG = P*POTLC - EUA - EUB
        U = COPP*(CHIA+CHIB)/(F1 + COPP*(CHIA+CHIB))
        W = SOPP0+U*(SOPP1+U*(SOPP2+U*(SOPP3+U*SOPP4)))
        EAB = EUEG*W
        DUDCHIA =COPP/(F1 + COPP*(CHIA+CHIB))**2
        DUDCHIB =COPP/(F1 + COPP*(CHIA+CHIB))**2
        DUDPA= DUDCHIA*CHIAP
        DUDPB= DUDCHIB*CHIBP
        DUDGA= DUDCHIA*CHIAG
        DUDGB= DUDCHIB*CHIBG
        DWDU =SOPP1+U*(F2*SOPP2+U*(F3*SOPP3+U*F4*SOPP4))
        DWDPA= DWDU*DUDPA
        DWDPB= DWDU*DUDPB
        DWDGA= DWDU*DUDGA
        DWDGB= DWDU*DUDGB
        EUEGPA = POTLC + P*DLDS*RSP + P*DLDZ*DZDA - EUPA
        EUEGPB = POTLC + P*DLDS*RSP + P*DLDZ*DZDB - EUPB
        DEABDPA = EUEGPA*W + EUEG*DWDPA
        DEABDGAA = EUEG*DWDGA
        DEABDPB = EUEGPB*W + EUEG*DWDPB
        DEABDGBB = EUEG*DWDGB
      ELSE
         EAB = 0._q
         DEABDPA = 0._q
         DEABDGAA = 0._q
         DEABDPB = 0._q
         DEABDGBB = 0._q
      ENDIF
      Ec_M06 = FA + FB + EAB + Ec_VS
      VCD1 = FPA + DEABDPA + VSD1
      VCD2 = FPB + DEABDPB + VSD2
      VCDD1 = F2*(FGA + DEABDGAA)*DRU + VSDD1
      VCDD2 = F2*(FGB + DEABDGBB)*DRD + VSDD2
      AMUCD1 = F2 * FTA + VSAMUCD1
      AMUCD2 = F2 * FTB + VSAMUCD2
!debug
!      if (abs(Ec_M06) .gt. 0.00001_q) then
!       write (*,*) "RU, RD, DRU, DRD, TAUU, TAUD", RU, RD, DRU, DRD, TAUU, TAUD
!       write (*,*) " Ec_M061, Ec_VS, Ec_M06", FA + FB + EAB, Ec_VS, Ec_M06
!       write (*,*) "VCD1, VCD2  ", VCD1, VCD2
!       write (*,*) "VCDD1, VCDD2 =", VCDD1, VCDD2
!       write (*,*) "AMUCD1,  AMUCD2 = ", AMUCD1,  AMUCD2
!       write (*,*) "VSD1,  VSD2", VSD1,  VSD2
!       write (*,*) "VSDD1, VSDD2", VSDD1, VSDD2
!       write (*,*) "VSAMUCD1, VSAMUCD2", VSAMUCD1, VSAMUCD2
!       write (*,*) "M06C FPA + DEABDPA", FPA + DEABDPA
!       write (*,*) "M06C FPB + DEABDPB", FPB + DEABDPB
!       write (*,*) "FGA + DEABDGAA", FGA + DEABDGAA
!       write (*,*) "FGB + DEABDGBB", FGB + DEABDGBB
!       write (*,*) "FTA, FTB", FTA, FTB
!       stop
!      ENDIF
      RETURN
      END SUBROUTINE VM06c

      SUBROUTINE M06CSS(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,CHI,EUEGP, &
     &                   CHIP,CHIG,IJZY)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : CORPBE
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(M06CSS)

!     COMPUTE THE SAME-SPIN PART OF THE M06 CORRELATION FUNCTIONAL FOR ONE GRID
!     POINT AND ONE SPIN-CASE.
!
      PARAMETER (ZERO = 0._q)
      PARAMETER (F1=1._q)
      PARAMETER (F2=2._q)
      PARAMETER (F3=3._q)
      PARAMETER (F4=4._q)
      PARAMETER (F5=5._q)
      PARAMETER (F6=6._q)
      PARAMETER (F8=8._q)
      PARAMETER (F11=11._q)
      PARAMETER (CSS = 0.06_q)
      PARAMETER (PT25 = 0.25_q)
GPACC !$ACC ROUTINE(CORPBE) SEQ
!
      SS=1._q
      IF (IJZY.EQ.1) THEN
!     PARAMETERS FOR M06-L CORRELATION
         SSS0=  5.349466e-01_q
         SSS1=  5.396620e-01_q
         SSS2=  -3.161217e+01_q
         SSS3=  5.149592e+01_q
         SSS4=  -2.919613e+01_q
      ELSEIF (IJZY.EQ.2) THEN
!     PARAMETERS FOR M06-HF CORRELATION
         SSS0=  1.023254e-01_q
         SSS1=  -2.453783e+00_q
         SSS2=  2.913180e+01_q
         SSS3=  -3.494358e+01_q
         SSS4=  2.315955e+01_q
      ELSEIF (IJZY.EQ.3) THEN
!     PARAMETERS FOR M06 CORRELATION
         SSS0=  5.094055e-01_q
         SSS1=  -1.491085e+00_q
         SSS2=  1.723922e+01_q
         SSS3=  -3.859018e+01_q
         SSS4=  2.845044e+01_q
      ELSEIF (IJZY.EQ.4) THEN
!     PARAMETERS FOR M06-2X CORRELATION
         SSS0=  3.097855e-01_q
         SSS1=  -5.528642e+00_q
         SSS2=  1.347420e+01_q
         SSS3=  -3.213623e+01_q
         SSS4=  2.846742e+01_q
      ENDIF

!      IF ((PX.LE.DTol).OR.(TX.LE.DTol).OR.(GX.LE.DTol))  THEN
       IF ((PX.LE.DTol))  THEN
        EUEG = ZERO
        CHI = ZERO
        EUEGP = ZERO
        CHIP = ZERO
        CHIG = ZERO
        PX = ZERO
        GX = ZERO
        TX = ZERO
        F  = ZERO
        FP = ZERO
        FG = ZERO
        FT = ZERO
      ELSE
!        PI = F4*ATAN(F1)
        PI34 = F3 / (F4*PI)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        FDUEG = (F3/F5)*(F6*PI*PI)**F23
        RS = (PI34/PX) ** F13
!        CALL LSDAC(RS,F1,POTLC,DLDS,DLDZ,D2LDSS,D2LDSZ,D2LDZZ)
! Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
       CALL CORPBE(RS,F1,PotLC,VCUPLDA,VCDNLDA,1._q,1._q,1._q,unknown,unknown,unknown,unknown,.FALSE.)
! OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
!      EPPGGA=PotLC
!      EPPGGA_D1=(VCUPLDA-PotLC)/RT
!      EPPGGA_D2=(VCDNLDA-PotLC)/RT
       dLdz=(VCUPLDA-VCDNLDA)/F2
       dLdS=(PotLC-F1*dLdZ-(VCUPLDA+VCDNLDA)/F2)*F3/RS
!c        Call lsdac(RS,F1,PotLC,dLdS,dLdZ)^M
!c  End of modification

        EUEG = PX*POTLC
        D = TX - PT25*GX/PX
!        DUEG = FDUEG*PX**F53
        CHI = GX/PX**F83
        U = CSS*CHI/(F1 + CSS*CHI)
        W = SSS0+U*(SSS1+U*(SSS2+U*(SSS3+U*SSS4)))
        FSCC=D/TX
        E = FSCC*W*EUEG
        F = E*SS
        RSP = -RS/(F3*PX)
        CHIG = F1/PX**F83
        CHIP = -F83*CHI/PX
        DFSCCP=PT25*GX/(TX*PX**2)
        DFSCCG=-PT25/(TX*PX)
        DFSCCT=PT25*GX/(PX*TX**2)
        DUDCHI=CSS/((F1+CSS*CHI)**2)
        DWDU=SSS1+U*(F2*SSS2+U*(F3*SSS3+U*F4*SSS4))
        DWDP=DWDU*DUDCHI*CHIP
        DWDG=DWDU*DUDCHI*CHIG
        EUEGP = POTLC + PX*DLDS*RSP
        FP = SS*(DFSCCP*W*EUEG            &
     &                 + FSCC*DWDP*EUEG   &
     &                 + FSCC*W*EUEGP)
        FG = SS*(DFSCCG*W*EUEG            &
     &                 + FSCC*DWDG*EUEG)

        FT = SS*(DFSCCT*W*EUEG)
       ENDIF

       RETURN
       END SUBROUTINE M06CSS



      SUBROUTINE VS98C(DTol, &
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_VS,VSD1,VSDD1,VSD2,VSDD2,VSAMUCD1,VSAMUCD2,IJZY)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : CORPBE, GVT4
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(VS98C)
      REAL(q) KAB
!     EVALUATE THE INTERPOLATION FUNCTION FOR PW91 LOCAL CORRELATION.
      PARAMETER (F1=1._q)
      PARAMETER (F2=2._q)
      PARAMETER (F3=3._q)
      PARAMETER (F4=4._q)
      PARAMETER (GAB=0.00304966e0_q)
      PARAMETER (CF=9.115599720e0_q)
GPACC !$ACC ROUTINE(VS98SS) SEQ
GPACC !$ACC ROUTINE(CORPBE) SEQ
GPACC !$ACC ROUTINE(GVT4) SEQ


!     PARAMETERS FOR VS98
      IF (IJZY.EQ.1) THEN
              R7=   7.035010e-01_q
              R8=   7.694574e-03_q
              R9=   5.152765e-02_q
              R10=   3.394308e-05_q
              R11=  -1.269420e-03_q
              R12=   1.296118e-03_q
!     PARAMETERS FOR M06-L
      ELSEIF (IJZY.EQ.2) THEN
              R7=      3.957626e-01_q
              R8=      -5.614546e-01_q
              R9=      1.403963e-02_q
              R10=     9.831442e-04_q
              R11=     -3.577176e-03_q
              R12=     0.000000e+00_q
!     PARAMETERS FOR M06-HF
      ELSEIF (IJZY.EQ.3) THEN
              R7=    -6.746338e-01_q
              R8=    -1.534002e-01_q
              R9=    -9.021521e-02_q
              R10=   -1.292037e-03_q
              R11=   -2.352983e-04_q
              R12=   0.000000e+00_q

!     PARAMETERS FOR M06
      ELSEIF (IJZY.EQ.4) THEN
               R7= -2.741539e+00_q
               R8= -6.720113e-01_q
               R9= -7.932688e-02_q
               R10=1.918681e-03_q
               R11=-2.032902e-03_q
               R12=0.000000e+00_q

!     PARAMETERS FOR M06-2X
      ELSEIF (IJZY.EQ.5) THEN
              R7=  1.166404e-01_q
              R8=  -9.120847e-02_q
              R9=  -6.726189e-02_q
              R10= 6.720580e-05_q
              R11= 8.448011e-04_q
              R12= 0.000000e+00_q
      ENDIF
!      PI = F4*ATAN(F1)
      F6=6._q
      F43 = F4 / F3
      PI34 = F3 / (F4*PI)
      F13 = F1 / F3

!     parallel spin case UP
        PA = RU
      IF (RU .GT. DTol .AND. TAUU .GT. DTol) THEN
        GAA = DRU*DRU
        TAUA = F2*TAUU
        CALL VS98SS(DTol,PA,GAA,TAUA,FA,FPA,FGA,FTA,EUA,ZA,&
     &                CHIA,EUPA,CHIAP,CHIAG,ZAP,ZAT,IJZY)
      ELSE
       FA = 0._q
       FPA = 0._q
       FGA = 0._q
       FTA = 0._q
       EUA = 0._q
       ZA = 0._q
       CHIA = 0._q
       EUPA = 0._q
       CHIAP = 0._q
       CHIAG = 0._q
       ZAP = 0._q
       ZAT = 0._q
      ENDIF
!     parallel spin case DOWN
      PB = RD
      IF (RD .GT. DTol .AND. TAUD .GT. DTol) THEN
      GBB = DRD*DRD
      TAUB= F2*TAUD
      CALL VS98SS(DTol,PB,GBB,TAUB,FB,FPB,FGB,FTB,EUB,ZB,&
     &                CHIB,EUPB,CHIBP,CHIBG,ZBP,ZBT,IJZY)
      ELSE
       FB = 0._q
       FPB = 0._q
       FGB = 0._q
       FTB = 0._q
       EUB = 0._q
       ZB = 0._q
       CHIB = 0._q
       EUPB = 0._q
       CHIBP = 0._q
       CHIBG = 0._q
       ZBP = 0._q
       ZBT = 0._q
      ENDIF


!    antiparallel spin case
      P=PA+PB
!      IF (PB.GT.DTol.AND.PA.GT.DTol.AND.TAUA.GT.DTol.AND.TAUB.GT.DTol) THEN
       IF (PB.GT.DTol.AND.PA.GT.DTol) THEN
        RS = (PI34/P) ** F13
        RSP = -RS/(F3*P)
        ZETA = (PA-PB)/P
        DZDA = (F1-ZETA)/P
        DZDB = (-F1-ZETA)/P
!        CALL LSDAC(RS,ZETA,POTLC,DLDS,DLDZ,D2LDSS,D2LDSZ, &
!     &      D2LDZZ)
! Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
        CALL CORPBE(RS,ZETA,PotLC,VCUPLDA,VCDNLDA,1._q,1._q,1._q,unknown,unknown,unknown,unknown,.FALSE.)
! OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
!      EPPGGA=PotLC
!      EPPGGA_D1=(VCUPLDA-PotLC)/RT
!      EPPGGA_D2=(VCDNLDA-PotLC)/RT
        dLdz=(VCUPLDA-VCDNLDA)/F2
        dLdS=(PotLC-ZETA*dLdZ-(VCUPLDA+VCDNLDA)/F2)*F3/RS

        EUEG = P*POTLC - EUA - EUB
        ZAB = ZA + ZB
        XAB = CHIA+CHIB
        KAB = F1 + GAB*(XAB+ZAB)
        XK = XAB/KAB
        ZK = ZAB/KAB
        CALL GVT4(GCAB,DGDX,DGDZ,XK,ZK,KAB,GAB,R7,R8,R9,R10,R11,R12)
        EAB =  GCAB*EUEG
        DGDPA = DGDX*CHIAP + DGDZ*ZAP
        DGDGA = DGDX*CHIAG
        DGDTA = DGDZ*ZAT
        DGDPB = DGDX*CHIBP + DGDZ*ZBP
        DGDGB = DGDX*CHIBG
        DGDTB = DGDZ*ZBT
        EUEGPA = POTLC + P*DLDS*RSP + P*DLDZ*DZDA - EUPA
        EUEGPB = POTLC + P*DLDS*RSP + P*DLDZ*DZDB - EUPB
        DEABDPA =  (EUEGPA*GCAB + EUEG*DGDPA)
        DEABDPB =  (EUEGPB*GCAB + EUEG*DGDPB)
        DEABDGAA =  EUEG*DGDGA
        DEABDGBB =  EUEG*DGDGB
        DEABDTA =   EUEG*DGDTA
        DEABDTB =   EUEG*DGDTB
      ELSE
        EAB = 0._q
        DEABDPA = 0._q
        DEABDPB = 0._q
        DEABDGAA = 0._q
        DEABDGBB = 0._q
        DEABDTA =  0._q
        DEABDTB =  0._q
      ENDIF
      Ec_VS = EAB + FA + FB
      VSD1 = FPA + DEABDPA
      VSD2 = FPB + DEABDPB
      VSDD1 = F2*(FGA + DEABDGAA)*DRU
      VSDD2 = F2*(FGB + DEABDGBB)*DRD
      VSAMUCD1 = F2 * (FTA + DEABDTA)
      VSAMUCD2 = F2 * (FTB + DEABDTB)

      RETURN
      END SUBROUTINE VS98C


      SUBROUTINE VS98SS(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Z,CHI,EUEGP, &
     &                   CHIP,CHIG,ZP,ZT,IJZY)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE ggalib, only : CORPBE, GVT4
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(VS98SS)


!     COMPUTE THE SAME-SPIN PART OF THE M06 CORRELATION FUNCTIONAL FOR ONE GRID
!     POINT AND ONE SPIN-CASE.



      REAL(q) KC
      PARAMETER (ZERO=0._q)
      PARAMETER (ONE=1._q)
      PARAMETER (F1=1._q)
      PARAMETER (F2=2._q)
      PARAMETER (F3=3._q)
      PARAMETER (F4=4._q)
      PARAMETER (F5=5._q)
      PARAMETER (F6=6._q)
      PARAMETER (F8=8._q)
      PARAMETER (F11=11._q)
      PARAMETER (GCC = 0.00515088_q)
      PARAMETER (CF = 9.11559972_q)
      PARAMETER (PT25 = 0.25_q)
GPACC !$ACC ROUTINE(CORPBE) SEQ
GPACC !$ACC ROUTINE(GVT4) SEQ

      F4O3 = F4/F3
!     PARAMETERS FOR VS98
      IF (IJZY.EQ.1) THEN
              R13=   3.270912e-01_q
              R14=  -3.228915e-02_q
              R15=  -2.942406e-02_q
              R16=   2.134222e-03_q
              R17=  -5.451559e-03_q
              R18=   1.577575e-02_q
!     PARAMETERS FOR M06-L
      ELSEIF (IJZY.EQ.2) THEN
              R13=   4.650534e-01_q
              R14=   1.617589e-01_q
              R15=   1.833657e-01_q
              R16=   4.692100e-04_q
              R17=  -4.990573e-03_q
              R18=   0.000000e+00_q
!     PARAMETERS FOR M06-HF
      ELSEIF (IJZY.EQ.3) THEN
              R13=   8.976746e-01_q
              R14=  -2.345830e-01_q
              R15=   2.368173e-01_q
              R16=  -9.913890e-04_q
              R17=  -1.146165e-02_q
              R18=   0.000000e+00_q
!     PARAMETERS FOR M06
      ELSEIF (IJZY.EQ.4) THEN
               R13=  4.905945e-01_q
               R14= -1.437348e-01_q
               R15=  2.357824e-01_q
               R16=  1.871015e-03_q
               R17= -3.788963e-03_q
               R18=  0.000000e+00_q
!     PARAMETERS FOR M06-2X
      ELSEIF (IJZY.EQ.5) THEN
              R13=  6.902145e-01_q
              R14=  9.847204e-02_q
              R15=  2.214797e-01_q
              R16= -1.968264e-03_q
              R17= -6.775479e-03_q
              R18=  0.000000e+00_q
      ENDIF

!      IF(PX.LE.DTol.OR.TX.LE.DTol.OR.GX.LE.DTol) THEN
       IF(PX.LE.DTol) THEN
        EUEG = ZERO
        CHI = ZERO
        EUEGP = ZERO
        CHIP = ZERO
        CHIG = ZERO
        PX = ZERO
        GX = ZERO
        TX = ZERO
        F  = ZERO
        FP = ZERO
        FG = ZERO
        FT = ZERO
        Z  = ZERO
        ZP = ZERO
        ZT = ZERO
      ELSE
!        PI = F4*ATAN(F1)
        PI34 = F3 / (F4*PI)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        RHOO = PX
        RRHO = ONE/RHOO
        RHO43 = RHOO**F4O3
        RHO13 = RHO43*RRHO
        RHO53 = RHOO**F53
        RHO83 = RHO53*RHOO

        RS = (PI34/PX) ** F13
!        CALL LSDAC(RS,F1,POTLC,DLDS,DLDZ,D2LDSS,D2LDSZ,D2LDZZ)

! Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
        CALL CORPBE(RS,F1,PotLC,VCUPLDA,VCDNLDA,1._q,1._q,1._q,unknown,unknown,unknown,unknown,.FALSE.)
! OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
!      EPPGGA=PotLC
!      EPPGGA_D1=(VCUPLDA-PotLC)/RT
!      EPPGGA_D2=(VCDNLDA-PotLC)/RT

        dLdz=(VCUPLDA-VCDNLDA)/F2
        dLdS=(PotLC-F1*dLdZ-(VCUPLDA+VCDNLDA)/F2)*F3/RS
!          Call lsdac(RS,F1,PotLC,dLdS,dLdZ)
!  End of modification
        EUEG = PX*POTLC
        CHI = GX/RHO83
        Z = (TX/RHO53) - CF
        KC = F1 + GCC*(CHI + Z)
        XK = CHI/KC
        ZK = Z/KC
        D = F1 - CHI/(F4*(Z + CF))
        CALL GVT4(GC,DGDX,DGDZ,XK,ZK,KC,GCC,R13,R14,R15,R16,R17,R18)
        E = D*EUEG*GC
!         WRITE (*,*) "CHI, Z, GC", CHI, Z, GC
        F = E
        RSP = -RS/(F3*PX)
        CHIG = F1/PX**F83
        CHIP = -F83*CHI/PX
        ZP = -F53 * TX/RHO83
        ZT =  F1/RHO53
        DZ = CHI/(F4*(Z + CF)*(Z + CF))
        DX = -F1/(F4*(Z + CF))
        DP = DZ*ZP + DX*CHIP
        DG = DX*CHIG
        DT = DZ*ZT
        DGDP = DGDX*CHIP + DGDZ*ZP
        DGDG = DGDX*CHIG
        DGDT = DGDZ*ZT
        EUEGP = POTLC + PX*DLDS*RSP
        FP = DP*EUEG*GC + D*EUEGP*GC + D*EUEG*DGDP
        FG = DG*EUEG*GC + D*EUEG*DGDG
        FT = DT*EUEG*GC + D*EUEG*DGDT
       ENDIF
       RETURN
       END SUBROUTINE VS98SS


!************************ SUBROUTINE MSx family ****************************
!
!***********************************************************************

      SUBROUTINE VMSXx(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ex_revTPSS,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2, &
     &   CMSX_RKAPPA,CMSX_CFC,CMSX_CFE,CMSX_TAUREG_X,ID_XC)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE xc_name
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VMSXx)

! the following parameters are given by Perdew et.al.
! the coefficients for the exchange part
! modified by Jianwei Sun 06/29/2011
      REAL(q) :: RKAPPA
      REAL(q) :: CFC
      PARAMETER (CFD=1.0_q)
      PARAMETER (CFG=0.0_q)
      REAL(q) :: CFE
      PARAMETER (CFMU=10.0_q/81.0_q)
      PARAMETER (CFN=1.0_q)
      PARAMETER (CFM=1.0_q)

!end of modification

! other parameters
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)

      RKAPPA=CMSX_RKAPPA
      CFC=CMSX_CFC
      CFE=CMSX_CFE
      ETA=CMSX_TAUREG_X

      VXD1=0._q;VXD2=0._q;
      VXDD1=0._q;VXDD2=0._q;
      AMUXD1=0._q;AMUXD2=0._q;


      CX1=10._q/81._q
      CX2=146._q/2025._q
      CX3=73._q/405._q
      CFE12=SQRT(CFE)

! Suspect that TAUWU and TAUWD are not well described. Use TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUWU_TEMP=0.125_q*DRU**2._q/RU
!     TAUWU_TEMP=MIN(TAUWU_TEMP,TAUU)

      TAUWD_TEMP=0.125_q*DRD**2._q/RD
!     TAUWD_TEMP=MIN(TAUWD_TEMP,TAUD)

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RU
      DRHO=TWO*DRU
      TAUW_RHO=TWO*TAUWU_TEMP
      TAU_RHO=TWO*TAUU

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P

      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
      ALPHA = (TAU_RHO - TAUW_RHO)/(TAU0 + ETA*TAUW_RHO)
      IF ((ID_XC==ID_XC_MS0).OR. &
          (ID_XC==ID_XC_MS0_X).OR. &
          (ID_XC==ID_XC_MS1).OR. &
          (ID_XC==ID_XC_MS1_X).OR. &
          (ID_XC==ID_XC_MS2).OR. &
          (ID_XC==ID_XC_MS2_X)) then
! MS0/1/2 style functional
          FX_HMGN=ONE+RKAPPA-RKAPPA/(ONE+(CFMU*P/RKAPPA))
          D_FXHMGN_P=CFMU/(ONE+CFMU*P/RKAPPA)**TWO
          FX_1EL=ONE+RKAPPA-RKAPPA/(ONE+((CFMU*P+CFC)/RKAPPA))
          D_FX1EL_P=CFMU/(ONE+(CFMU*P+CFC)/RKAPPA)**TWO
      ELSEIF ((ID_XC==ID_XC_MSPBEL).OR.(ID_XC==ID_XC_RMSPBEL)) THEN
! PBE style functional
          FX_HMGN=(ONE+CFMU*P/(ONE+CFMU*P/RKAPPA))
          D_FXHMGN_P=CFMU/(ONE+CFMU*P/RKAPPA)**TWO
          FX_1EL=(ONE+(CFMU*P+CFC)/(ONE+(CFMU*P+CFC)/RKAPPA))
          D_FX1EL_P=CFMU/(ONE+(CFMU*P+CFC)/RKAPPA)**TWO
      ELSEIF ((ID_XC==ID_XC_MSRPBEL).OR.(ID_XC==ID_XC_RMSRPBEL)) THEN
! RPBE style functional
          FX_HMGN=(ONE+RKAPPA*(ONE-EXP(-CFMU*P/RKAPPA)))
          D_FXHMGN_P=CFMU*EXP(-CFMU*P/RKAPPA)
          FX_1EL=(ONE+RKAPPA*(ONE-EXP(-(CFMU*P+CFC)/RKAPPA)))
          D_FX1EL_P=CFMU*EXP(-(CFMU*P+CFC)/RKAPPA)
      ELSEIF ((ID_XC==ID_XC_MSB86BL).OR.(ID_XC==ID_XC_RMSB86BL)) THEN
! B86R style functional
          FX_HMGN=(ONE+CFMU*P/(ONE+CFMU*P/RKAPPA)**0.8_q)
          D_FXHMGN_P=0.2_q*CFMU*RKAPPA**0.8_q*(CFMU*P+5.0_q*RKAPPA) &
           / (CFMU*P+RKAPPA)**1.8_q
          FX_1EL=(ONE+(CFMU*P+CFC)/(ONE+(CFMU*P+CFC)/RKAPPA)**0.8_q)
          D_FX1EL_P=0.2_q*CFMU*RKAPPA**0.8_q*(CFMU*P+CFC+&
           5.0_q*RKAPPA) / (CFMU*P+CFC+RKAPPA)**1.8_q
      ENDIF
      
      ALPHA2=ALPHA*ALPHA
      ALPHA3=ALPHA*ALPHA2
      ALPHA4=ALPHA2*ALPHA2
      ALPHA6=ALPHA4*ALPHA2
      OMA=ONE-ALPHA2
      OMA2=OMA*OMA
      OMA3=OMA2*OMA
      OPA=ONE+CFD*ALPHA3+CFE*ALPHA6
      FINT=OMA3/OPA

      FX=FX_HMGN+FINT*(FX_1EL-FX_HMGN)

!    NOW, DERIVATIVES COME
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q

      DALPHAD = THRD5*TAU0/(RHO*(TAU0+ETA*TAUW_RHO))*(P-ALPHA+ETA*P*ALPHA)
      DALPHADD = -0.25_q*DRHO*(ALPHA*ETA+ONE)/(RHO*(TAU0+ETA*TAUW_RHO))
      DALPHADTAU = ONE/(TAU0 + ETA*TAUW_RHO)

! d Falpha / dy
      OPA2=OPA*OPA
      D_FINT_ALPHA=-6._q*ALPHA*OMA2/OPA-(THREE*CFD+6._q*CFE*ALPHA3)*ALPHA2*OMA3/OPA2


      D_FINT_D=D_FINT_ALPHA*DALPHAD
      D_FINT_DD=D_FINT_ALPHA*DALPHADD
      D_FINT_DTAU=D_FINT_ALPHA*DALPHADTAU


      COE_FX_P=D_FXHMGN_P+FINT*(D_FX1EL_P-D_FXHMGN_P)
      DFX_D=COE_FX_P*DPD+(FX_1EL-FX_HMGN)*D_FINT_D
      DFX_DD=COE_FX_P*DPDD+(FX_1EL-FX_HMGN)*D_FINT_DD
      DFX_DTAU=COE_FX_P*DPDTAU+(FX_1EL-FX_HMGN)*D_FINT_DTAU


!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
      EX_REVTPSS=0._q
      EX_REVTPSS=EX_REVTPSS+EXLDA*FX

! spin down
! IN EXD1(2*RD), TAUWD AND TAUD SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RD
      DRHO=TWO*DRD
      TAUW_RHO=TWO*TAUWD_TEMP
      TAU_RHO=TWO*TAUD

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P

      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
      ALPHA = (TAU_RHO - TAUW_RHO)/(TAU0 + ETA*TAUW_RHO)

      IF ((ID_XC==ID_XC_MS0).OR. &
          (ID_XC==ID_XC_MS0_X).OR. &
          (ID_XC==ID_XC_MS1).OR. &
          (ID_XC==ID_XC_MS1_X).OR. &
          (ID_XC==ID_XC_MS2).OR. &
          (ID_XC==ID_XC_MS2_X)) then
! MS0/1/2 style functional
          FX_HMGN=ONE+RKAPPA-RKAPPA/(ONE+(CFMU*P/RKAPPA))
          D_FXHMGN_P=CFMU/(ONE+CFMU*P/RKAPPA)**TWO
          FX_1EL=ONE+RKAPPA-RKAPPA/(ONE+((CFMU*P+CFC)/RKAPPA))
          D_FX1EL_P=CFMU/(ONE+(CFMU*P+CFC)/RKAPPA)**TWO
      ELSEIF ((ID_XC==ID_XC_MSPBEL).OR.(ID_XC==ID_XC_RMSPBEL)) THEN
! PBE style functional
          FX_HMGN=(ONE+CFMU*P/(ONE+CFMU*P/RKAPPA))
          D_FXHMGN_P=CFMU/(ONE+CFMU*P/RKAPPA)**TWO
          FX_1EL=(ONE+(CFMU*P+CFC)/(ONE+(CFMU*P+CFC)/RKAPPA))
          D_FX1EL_P=CFMU/(ONE+(CFMU*P+CFC)/RKAPPA)**TWO
      ELSEIF ((ID_XC==ID_XC_MSRPBEL).OR.(ID_XC==ID_XC_RMSRPBEL)) THEN
! RPBE style functional
          FX_HMGN=(ONE+RKAPPA*(ONE-EXP(-CFMU*P/RKAPPA)))
          D_FXHMGN_P=CFMU*EXP(-CFMU*P/RKAPPA)
          FX_1EL=(ONE+RKAPPA*(ONE-EXP(-(CFMU*P+CFC)/RKAPPA)))
          D_FX1EL_P=CFMU*EXP(-(CFMU*P+CFC)/RKAPPA)
      ELSEIF ((ID_XC==ID_XC_MSB86BL).OR.(ID_XC==ID_XC_RMSB86BL)) THEN
! B86R style functional
          FX_HMGN=(ONE+CFMU*P/(ONE+CFMU*P/RKAPPA)**0.8_q)
          D_FXHMGN_P=0.2_q*CFMU*RKAPPA**0.8_q*(CFMU*P+5.0_q*RKAPPA) &
           / (CFMU*P+RKAPPA)**1.8_q
          FX_1EL=(ONE+(CFMU*P+CFC)/(ONE+(CFMU*P+CFC)/RKAPPA)**0.8_q)
          D_FX1EL_P=0.2_q*CFMU*RKAPPA**0.8_q*(CFMU*P+CFC+&
           5.0_q*RKAPPA) / (CFMU*P+CFC+RKAPPA)**1.8_q
      ENDIF

      ALPHA2=ALPHA*ALPHA
      ALPHA3=ALPHA*ALPHA2
      ALPHA4=ALPHA2*ALPHA2
      ALPHA6=ALPHA4*ALPHA2
      OMA=ONE-ALPHA2
      OMA2=OMA*OMA
      OMA3=OMA2*OMA
      OPA=ONE+CFD*ALPHA3+CFE*ALPHA6
      FINT=OMA3/OPA

      FX=FX_HMGN+FINT*(FX_1EL-FX_HMGN)

!    NOW, DERIVATIVES COME
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q

      DALPHAD = THRD5*TAU0/(RHO*(TAU0+ETA*TAUW_RHO))*(P-ALPHA+ETA*P*ALPHA)
      DALPHADD = -0.25_q*DRHO*(ALPHA*ETA+ONE)/(RHO*(TAU0+ETA*TAUW_RHO))
      DALPHADTAU = ONE/(TAU0 + ETA*TAUW_RHO)


! d Falpha / dy
      OPA2=OPA*OPA
      D_FINT_ALPHA=-6._q*ALPHA*OMA2/OPA-(THREE*CFD+6._q*CFE*ALPHA3)*ALPHA2*OMA3/OPA2



      D_FINT_D=D_FINT_ALPHA*DALPHAD
      D_FINT_DD=D_FINT_ALPHA*DALPHADD
      D_FINT_DTAU=D_FINT_ALPHA*DALPHADTAU


      COE_FX_P=D_FXHMGN_P+FINT*(D_FX1EL_P-D_FXHMGN_P)
      DFX_D=COE_FX_P*DPD+(FX_1EL-FX_HMGN)*D_FINT_D
      DFX_DD=COE_FX_P*DPDD+(FX_1EL-FX_HMGN)*D_FINT_DD
      DFX_DTAU=COE_FX_P*DPDTAU+(FX_1EL-FX_HMGN)*D_FINT_DTAU


!   OUTPUT THE VXD2,VXDD2 AND AMUXD2
      VXD2=EXDLDA*FX+EXLDA*DFX_D
      VXDD2=EXLDA*DFX_DD
      AMUXD2=EXLDA*DFX_DTAU

      EX_REVTPSS=EX_REVTPSS+EXLDA*FX
      EX_REVTPSS=EX_REVTPSS/TWO
      RETURN
      END SUBROUTINE VMSXx


      SUBROUTINE VMSXc(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_revTPSS,VCD1,VCD2,VCDD12,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : CORGGA_REVTPSS
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VMSXc)

! the following parameters are given by Perdew et.al.
! the coefficients for the correlatin part
      PARAMETER (CFD=2.8_q)

! other parameters
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
GPACC !$ACC ROUTINE(CORGGA_REVTPSS) SEQ

      VCD1=0._q;VCD2=0._q;
      VCDD12=0._q
      AMUCD1=0._q;AMUCD2=0._q;

      RT=RU+RD

      CALL CORGGA_REVTPSS(RU,RD,DRU,DRD,DRT,EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD12)

      VCD1=EPPGGA+RT*EPPGGA_D1

      VCD2=EPPGGA+RT*EPPGGA_D2

      VCDD12=RT*EPPGGA_DD12


      EC_REVTPSS=RT*EPPGGA


      RETURN
      END SUBROUTINE VMSXc


!************************ SUBROUTINE VSCANx ****************************
!
!> calculates the first order derivatives of Ex wrt n and |grad(n)|
!> Perdew et. al. PRL (2009)
!>
!> Written by Jianwei Sun
!>
!> everything in Hartree units
!>
!> ATTENTION: Every values are passed "as they are", i.e. including
!> possibly unphysical numerical errors (e.g. negative charge densities)
!> values need to be checked accordingly
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VXD1 VXD2                    THE DERIVATIVES OF EX WRT n
!> VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n|
!> AMUXD1, AMUXD2               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VSCANx(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ex_SCAN,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,SCAN_TAUREG_X,SCAN_ALPREG,ID_XC)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VSCANx)
      PARAMETER (CFC1X=0.667_q)
      PARAMETER (CFC2X=0.800_q)
      PARAMETER (CFD1X=1.240_q)
      PARAMETER (CFK1=0.065_q)
! other parameters
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (TWO53=TWO**THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
GPACC !$ACC ROUTINE(metaggafx) SEQ

      VXD1=0._q;VXD2=0._q;
      VXDD1=0._q;VXDD2=0._q;
      AMUXD1=0._q;AMUXD2=0._q;


      CX1=10._q/81._q
      CX2=146._q/2025._q
      CX3=73._q/405._q

      TAUREG = SCAN_TAUREG_X

! Suspect that TAUWU and TAUWD are not well described. Use TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUWU_TEMP=0.125_q*DRU**2._q/RU
!     TAUWU_TEMP=MIN(TAUWU_TEMP,TAUU)

      TAUWD_TEMP=0.125_q*DRD**2._q/RD
!     TAUWD_TEMP=MIN(TAUWD_TEMP,TAUD)

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RU
      DRHO=TWO*DRU
      TAUW_RHO=TWO*TAUWU_TEMP
      TAU_RHO=TWO*TAUU

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P

      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q
      TAU_UNIF=THREE/10._q*(THREE*PI**TWO)**THRD2*RHO**THRD5

! If TAUREG=0 and SCAN_ALPREG=0,
! RALPHA = ALPHA and D_RALPHA_DALPHA = 1
! This works for both SCAN and rSCAN
      DENOM = TAU0 + TAUREG
      ALPHA = (TAU_RHO - TAUW_RHO)/DENOM
      DALPHAD = THRD5*TAU0/(RHO*DENOM) * (P-ALPHA)
      DALPHADD = -DRHO/(FOUR*RHO*DENOM)
      DALPHADTAU = ONE/DENOM

      RALPHA = ALPHA**3/(ALPHA**2 + SCAN_ALPREG)
      D_RALPHA_DALPHA = ALPHA/(ALPHA**2 + SCAN_ALPREG) * (3*ALPHA - 2*RALPHA)
      DALPHAD = DALPHAD*D_RALPHA_DALPHA
      DALPHADD = DALPHADD*D_RALPHA_DALPHA
      DALPHADTAU = DALPHADTAU*D_RALPHA_DALPHA

! calculate the exchange enhancement factor and its derivatives wrt p
! and alpha
       call metaggafx(cfc1x,cfc2x,cfd1x,cfk1, &
     & p, ralpha, Fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2, ID_XC)

! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=fx1p*DPD+fx1a*DALPHAD
      DFX_DD=fx1p*DPDD+fx1a*DALPHADD
      DFX_DTAU=fx1a*DALPHADTAU



!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
      EX_SCAN=0._q
      EX_SCAN=EX_SCAN+EXLDA*FX

! spin down
! IN EXD1(2*RD), TAUWD AND TAUD SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RD
      DRHO=TWO*DRD
      TAUW_RHO=TWO*TAUWD_TEMP
      TAU_RHO=TWO*TAUD

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      P2=P*P

      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0._q
      TAU_UNIF=THREE/10._q*(THREE*PI**TWO)**THRD2*RHO**THRD5

! If TAUREG=0 and SCAN_ALPREG=0,
! RALPHA = ALPHA and D_RALPHA_DALPHA = 1
! This works for both SCAN and rSCAN
      DENOM = TAU0 + TAUREG
      ALPHA = (TAU_RHO - TAUW_RHO)/DENOM
      DALPHAD = THRD5*TAU0/(RHO*DENOM) * (P-ALPHA)
      DALPHADD = -DRHO/(FOUR*RHO*DENOM)
      DALPHADTAU = ONE/DENOM

      RALPHA = ALPHA**3/(ALPHA**2 + SCAN_ALPREG)
      D_RALPHA_DALPHA = ALPHA/(ALPHA**2 + SCAN_ALPREG) * (3*ALPHA - 2*RALPHA)
      DALPHAD = DALPHAD*D_RALPHA_DALPHA
      DALPHADD = DALPHADD*D_RALPHA_DALPHA
      DALPHADTAU = DALPHADTAU*D_RALPHA_DALPHA

! calculate the exchange enhancement factor and its derivatives wrt p
! and alpha
       call metaggafx(cfc1x,cfc2x,cfd1x,cfk1, &
     & p, ralpha, Fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2, ID_XC)

! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=fx1p*DPD+fx1a*DALPHAD
      DFX_DD=fx1p*DPDD+fx1a*DALPHADD
      DFX_DTAU=fx1a*DALPHADTAU



!   OUTPUT THE VXD2,VXDD2 AND AMUXD2
      VXD2=EXDLDA*FX+EXLDA*DFX_D
      VXDD2=EXLDA*DFX_DD
      AMUXD2=EXLDA*DFX_DTAU

      EX_SCAN=EX_SCAN+EXLDA*FX
      EX_SCAN=EX_SCAN/TWO
      RETURN
      END SUBROUTINE VSCANx


      SUBROUTINE metaggafx(cfc1,cfc2,cfd1,cfk1, &
     &            p, a, fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2, ID_XC)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE xc_name
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(metaggafx)
      INTEGER I

      REAL(q) :: CREGU(0:7),CREGU1(7),ALPHAI(0:7)
      PARAMETER (CREGU = (/ 1.0_q, -0.667_q, -0.4445555_q, -0.663086601049_q,&
                             &   1.451297044490_q, -0.887998041597_q,          &
                             &   0.234528941479_q, -0.023185843322_q /) )

      PARAMETER (cfk0=0.1740d0)
      PARAMETER (CFA1=4.9479d0)

      PARAMETER (CFMUAK=10.00d0/81.0d0)
      PARAMETER (CFB1=0.156632d0)
      PARAMETER (CFB2=0.12083d0)
      PARAMETER (CFB3=0.5d0)

      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)

      ALPHAI(0) = 1.0_q
      DO I=1,7
        ALPHAI(I) = ALPHAI(I-1)*A
        CREGU1(I) = I*CREGU(I)
      END DO

      fx=0.0d0
      fx1p=0.0d0
      fx1a=0.0d0
      fx2p2=0.0d0
      fx2pa=0.0d0
      fx2a2=0.0d0

      P2=P*P
      oma=1.0d0-a
      oma2=oma*oma

      HX0=1.0d0+CFK0

      cfb4=cfmuak**two/cfk1-0.112654d0
      if (cfb4 .gt. 0.0d0) then
         wfac=cfb4*p2*dexp(-cfb4*p/cfmuak)
         d_wfac_dp=cfb4*p*dexp(-cfb4*p/cfmuak)*(two-cfb4*p/cfmuak)
      else
         wfac=cfb4*p2*dexp(cfb4*p/cfmuak)
         d_wfac_dp=cfb4*p*dexp(cfb4*p/cfmuak)*(two+cfb4*p/cfmuak)
      endif
      vfac=cfb1*p+cfb2*oma*dexp(-cfb3*oma2)
      yfac=cfmuak*p+wfac+vfac**two
      hx1=1.0d0+cfk1-cfk1/(1.0d0+yfac/cfk1)

      FA=0.d0
      d_fa_da=0.d0

      ! SCAN switching function
      if ((ID_XC==ID_XC_SCAN).OR. &
          (ID_XC==ID_XC_SCAN_X)) then
        if (a.le.0.998_q) then
          fa = dexp(-cfc1*a/oma)
          d_fa_da = -cfc1*fa/oma2
        else if (a.ge.1.002_q) then
          fa = -cfd1*dexp(cfc2/oma)
          d_fa_da = cfc2*fa/oma2
        end if

      ! rSCAN switching function
      else if ((ID_XC==ID_XC_RSCAN).OR. &
               (ID_XC==ID_XC_RSCAN_X).OR. &
               (ID_XC==ID_XC_RSCANZT).OR. &
               (ID_XC==ID_XC_RSCANZT_X)) then
        if (a.le.0.0_q) then
          fa = dexp(-cfc1*a/oma)
          d_fa_da = -cfc1*fa/oma2
        else if (a.le.2.5_q) then
          fa = dot_product(cregu,alphai)
          d_fa_da = dot_product(cregu1,alphai(0:6))
        else
          fa = -cfd1*dexp(cfc2/oma)
          d_fa_da = cfc2*fa/oma2
        end if
      end if

      p14=p**(1.0d0/four)
      gx=1.d0
      if (p .gt. 0.d0) then
          gx=1.0d0-dexp(-cfa1/p14)
      endif

      Fx1 = hx1+FA*(hx0-hx1)

      Fx = Fx1*gx

      d_hx0_dp=0.d0

      d_vfac_dp=cfb1
      d_yfac_dp=cfmuak+d_wfac_dp+two*vfac*d_vfac_dp
      d_hx1_dp=d_yfac_dp/(1.0d0+yfac/cfk1)**two

      d_vfac_da=-cfb2*(1.0d0-two*cfb3*oma2)*dexp(-cfb3*oma2)
      d_yfac_da=two*vfac*d_vfac_da
      d_hx1_da=d_yfac_da/(1.0d0+yfac/cfk1)**two

      d_gx_dp=0.d0
      if (p .gt. 0.d0) then
          d_gx_dp=-cfa1/four/p/p14*dexp(-cfa1/p14)
      endif

      d_Fx1_dp=d_hx1_dp+FA*(d_hx0_dp-d_hx1_dp)

      d_Fx1_da=(1.0d0-FA)*d_hx1_da+d_FA_da*(hx0-hx1)


      fx1p=d_Fx1_dp*gx+Fx1*d_gx_dp

      fx1a=d_Fx1_da*gx

      RETURN
      END SUBROUTINE metaggafx


!************************ SUBROUTINE VSCANc ****************************
!
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VCD1 VCD2                    THE DERIVATIVES OF EC WRT n
!> VCDD1,VCDD2                  THE DERIVATIVES OF EC WRT |grad n|
!> AMUCD1, AMUCD2                   THE DERIVATIVES OF EC WRT TAU
!
!***********************************************************************

      SUBROUTINE VSCANc(&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD,&
     &   Ec_SCAN,VCD1,VCD2,VCDD,AMUCD1,AMUCD2,SCAN_TAUREG_C,SCAN_ALPREG,ID_XC)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL*8 (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(VSCANc)

! the coefficients for the exchange part
      PARAMETER (CFC1C=0.640d0)
      PARAMETER (CFC2C=1.500d0)
      PARAMETER (CFD1C=0.700d0)
GPACC !$ACC ROUTINE(Vmetaggac) SEQ



! other parameter
      PARAMETER (ONE=1.0d0)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
!      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)

      call Vmetaggac(cfc1c,cfc2c,cfd1c,&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD,&
     &   Ec_SCAN,VCD1,VCD2,VCDD,AMUCD1,AMUCD2,SCAN_TAUREG_C,SCAN_ALPREG,ID_XC)

      RETURN
      END SUBROUTINE VSCANc


      SUBROUTINE Vmetaggac(cfc1,cfc2,cfd1,&
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
     &   Ec_SCAN,VCD1,VCD2,VCDD,AMUCD1,AMUCD2,SCAN_TAUREG_C,SCAN_ALPREG,ID_XC)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      USE ggalib, only : CORGGA_0, CORGGA_1
      USE xc_name
      IMPLICIT REAL*8 (A-H,O-Z)
GPOMP  !$OMP DECLARE TARGET(Vmetaggac)

      INTEGER I
      REAL(q) :: CREGU(0:7),CREGU1(7),ALPHAI(0:7)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
!      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)
      PARAMETER (CREGU = (/ 1.0_q, -0.64_q, -0.4352_q,  -1.535685604549_q,&
        &              3.061560252175_q, -1.915710236206_q, 0.516884468372_q,&
        &             -0.051848879792_q /) )
GPACC !$ACC ROUTINE(CORGGA_0) SEQ
GPACC !$ACC ROUTINE(CORGGA_1) SEQ


      VCD1=0.0d0
      VCD2=0.0d0
      VCDD=0.0d0
      AMUCD1=0.0d0
      AMUCD2=0.0d0


      RT=RU+RD
      YA=DRU**2.0d0
      YB=DRD**2.0d0
      Y=DRT**2.0d0
      YC=(Y-YA-YB)/2.0d0
      TAUW=1.0d0/8.0d0*(Y/RT)
      TAU=TAUU+TAUD
      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.99999999999990d0),0.99999999999990d0)
      DS_ZETA=((1.0d0+ZETA)**THRD5+(1.0d0-ZETA)**THRD5)/TWO
      DX_ZETA=((1.0d0+ZETA)**THRD4+(1.0d0-ZETA)**THRD4)/TWO
      TAU0=0.30d0*(THREE*PISQ)**THRD2*RT**THRD5*DS_ZETA

      D_ZETA_D1=TWO*RD/RT**TWO
      D_ZETA_D2=-TWO*RU/RT**TWO
      D_DS_ZETA_D1=THRD5*((1.0d0+ZETA)**THRD2 &
     &             -(1.0d0-ZETA)**THRD2)*D_ZETA_D1/TWO
      D_DS_ZETA_D2=THRD5*((1.0d0+ZETA)**THRD2 &
     &             -(1.0d0-ZETA)**THRD2)*D_ZETA_D2/TWO
      DTAU0D1=THRD5*TAU0/RT+TAU0*D_DS_ZETA_D1/DS_ZETA
      DTAU0D2=THRD5*TAU0/RT+TAU0*D_DS_ZETA_D2/DS_ZETA

      DTAUWD1=-TAUW/RT
      DTAUWD2=DTAUWD1
      DTAUWDD=1.0d0/FOUR * DRT/RT

      DENOM = TAU0 + DS_ZETA*SCAN_TAUREG_C
      ALPHA = (TAU - TAUW)/DENOM
      D_ALPHA_D1 = (-DTAUWD1 - ALPHA*(DTAU0D1+SCAN_TAUREG_C*D_DS_ZETA_D1))/DENOM
      D_ALPHA_D2 = (-DTAUWD2 - ALPHA*(DTAU0D2+SCAN_TAUREG_C*D_DS_ZETA_D2))/DENOM
      D_ALPHA_DD=-DTAUWDD/DENOM
      D_ALPHA_DTAU=1.0d0/DENOM

      RALPHA = ALPHA**3/(ALPHA**2 + SCAN_ALPREG)
      D_RALPHA_DALPHA = ALPHA/(ALPHA**2 + SCAN_ALPREG) * (3*ALPHA - 2*RALPHA)
      D_ALPHA_D1 = D_ALPHA_D1*D_RALPHA_DALPHA
      D_ALPHA_D2 = D_ALPHA_D2*D_RALPHA_DALPHA
      D_ALPHA_DD = D_ALPHA_DD*D_RALPHA_DALPHA
      D_ALPHA_DTAU = D_ALPHA_DTAU*D_RALPHA_DALPHA

      OMA = 1.0_q - RALPHA
      OMA2 = OMA*OMA

      ALPHAI(0) = 1.0_q
      DO I=1,7
        ALPHAI(I) = ALPHAI(I-1)*RALPHA
        CREGU1(I) = I*CREGU(I)
      END DO

      FUNKC_ALPHA=0.0d0
      D_FUNKC_ALPHA_DALPHA=0.0d0
      IF ((ID_XC==ID_XC_SCAN).OR. &
          (ID_XC==ID_XC_SCAN_C)) THEN
        IF (RALPHA.LE.0.998_q) THEN
          FUNKC_ALPHA = DEXP(-CFC1*RALPHA/OMA)
          D_FUNKC_ALPHA_DALPHA = -CFC1*FUNKC_ALPHA/OMA2
        ELSE IF (RALPHA.GE.1.002_q) THEN
          FUNKC_ALPHA = -CFD1*DEXP(CFC2/OMA)
          D_FUNKC_ALPHA_DALPHA = CFC2*FUNKC_ALPHA/OMA2
        END IF
      ELSE IF ((ID_XC==ID_XC_RSCAN).OR. &
               (ID_XC==ID_XC_RSCAN_C).OR. &
               (ID_XC==ID_XC_RSCANZT).OR. &
               (ID_XC==ID_XC_RSCANZT_C)) THEN
        IF (RALPHA.LE.0.0_q) THEN
          FUNKC_ALPHA = DEXP(-CFC1*RALPHA/OMA)
          D_FUNKC_ALPHA_DALPHA = -CFC1*FUNKC_ALPHA/OMA2
        ELSE IF (RALPHA.LE.2.5_q) THEN
          FUNKC_ALPHA = DOT_PRODUCT(CREGU,ALPHAI)
          D_FUNKC_ALPHA_DALPHA = DOT_PRODUCT(CREGU1,ALPHAI(0:6))
        ELSE
          FUNKC_ALPHA = -CFD1*DEXP(CFC2/OMA)
          D_FUNKC_ALPHA_DALPHA = CFC2*FUNKC_ALPHA/OMA2
        END IF
      END IF


      CALL CORGGA_0(RU,RD,DRU,DRD,DRT,&
     &        EPPGGA0,EPPGGA0_D1,EPPGGA0_D2,EPPGGA0_DD)
      CALL CORGGA_1(RU,RD,DRU,DRD,DRT, &
     &       EPPGGA1,EPPGGA1_D1,EPPGGA1_D2,EPPGGA1_DD)


      EPP=EPPGGA1+FUNKC_ALPHA*(EPPGGA0-EPPGGA1)
      D_EPP_D1=EPPGGA1_D1 &
     &          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_D1 &
     &          +FUNKC_ALPHA*(EPPGGA0_D1-EPPGGA1_D1)
      D_EPP_D2=EPPGGA1_D2 &
     &          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_D2 &
     &          +FUNKC_ALPHA*(EPPGGA0_D2-EPPGGA1_D2)
      D_EPP_DD=EPPGGA1_DD &
     &          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_DD &
     &          +FUNKC_ALPHA*(EPPGGA0_DD-EPPGGA1_DD)
      D_EPP_DTAU=(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_DTAU



      VCD1=EPP+RT*D_EPP_D1
      VCD2=EPP+RT*D_EPP_D2
      VCDD=RT*D_EPP_DD
      AMUCD1=RT*D_EPP_DTAU
      AMUCD2=AMUCD1

      EC_SCAN=RT*EPP

      RETURN
      END SUBROUTINE Vmetaggac

!***********************************************************************
!> r^2-SCAN subroutines
!> Written by Daniel Mejia-Rodriguez
!> M2QM and QTP
!> University of Florida
!***********************************************************************

      SUBROUTINE VR2SCANx(RX,DRX,TAUX,Ex,VXD,VXDD,AMUXD)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant

      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(VR2SCANx)

      REAL(q) :: RX,DRX,TAUX,EX,VXD,VXDD,AMUXD
      REAL(q) :: C1X,C2X,CDX,CK1,CH0,CA1,CMU
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: ETA,CETA,C2,CDP24,PISQ,AX
      REAL(q) :: TAU_RHO,TAUW_RHO,TAU0,TAUW_TEMP
      REAL(q) :: RHO,DRHO,P,P2,P14,EXP1,EXP2,EXP3,EXP4
      REAL(q) :: FA,DFADA,ALPHA,DALPHAD,DALPHADD,DALPHADTAU
      REAL(q) :: DPD,DPDD,EXUNIF,EXLDA,EXDLDA,OMA,OMA2,X
      REAL(q) :: H,DHDX,DXDP,HDEN,HNUM,DHDP,DG,DGDP
      REAL(q) :: FX,DFXDP,DFXDA,DFX_D,DFX_DD,DFX_DTAU
      REAL(q) :: CREGU(0:7),CREGU1(7),RTEMP(0:7)
      INTEGER :: IFX

      PARAMETER (C1X=0.667_q)
      PARAMETER (C2X=0.800_q)
      PARAMETER (CDX=1.240_q)
      PARAMETER (CK1=0.065_q)

      PARAMETER (CH0=1.1740d0)
      PARAMETER (CA1=4.9479d0)

      PARAMETER (CMU=10.00d0/81.0d0)
! other parameters
      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (FIVE=5._q)
      PARAMETER (SIX=6._q)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
      PARAMETER (ETA=0.001_q)
      PARAMETER (CETA=20.0_q/27.0_q+ETA*THRD5)
      PARAMETER (CDP24=0.361_q**4)

      PARAMETER (CREGU = (/ ONE, -0.667_q, -0.4445555_q, &
     &                    -6.63086601049291E-1_q, 1.45129704448975_q, &
     &                    -8.87998041596655E-1_q, 2.34528941478571E-1_q, &
     &                    -2.31858433223407E-2_q/) )

      PARAMETER (CREGU1 = (/ -0.667_q, -TWO*0.4445555_q, &
     &             -THREE*6.63086601049291E-1_q, FOUR*1.45129704448975_q, &
     &             -5._q*8.87998041596655E-1_q, SIX*2.34528941478571E-1_q, &
     &             -7._q*2.31858433223407E-2_q/) )

      PARAMETER (C2 = -(CREGU1(1)+CREGU1(2)+CREGU1(3)+CREGU1(4)+CREGU1(5)+CREGU1(6)+CREGU1(7))*(ONE-CH0))

      RTEMP(0) = ONE

      VXD=0._q
      VXDD=0._q
      AMUXD=0._q


! Suspect that TAUWU and TAUWD are not well described. Use TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUW_TEMP=0.125_q*DRX**2._q/RX

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RX
      DRHO=TWO*DRX
      TAUW_RHO=TWO*TAUW_TEMP
      TAU_RHO=TWO*TAUX
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      DPD=-THRD8*P/RHO
      DPDD=TWO*DRHO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)

      ALPHA = (TAU_RHO - TAUW_RHO)/(TAU0 + ETA*TAUW_RHO)
      DALPHAD = THRD5*TAU0/(RHO*(TAU0+ETA*TAUW_RHO))*(P-ALPHA+ETA*P*ALPHA)
      DALPHADD = -0.25_q*DRHO*(ALPHA*ETA+ONE)/(RHO*(TAU0+ETA*TAUW_RHO))
      DALPHADTAU = ONE/(TAU0 + ETA*TAUW_RHO)

      P14 = SQRT(SQRT(P))
      P2 = P*P
      OMA = ONE - ALPHA
      OMA2 = OMA*OMA

      EXP1 = EXP(-P2/CDP24)
      EXP2 = 0.0_q
      IF (P14.GT.0.002_q) EXP2 = EXP(-CA1/P14)

      X = (CETA*C2*EXP1 + CMU)*P
      DXDP = CETA*C2*EXP1 + CMU - TWO*P2*CETA*C2*EXP1/CDP24

      HDEN = CK1 + X
      HNUM = HDEN + CK1*X
      H = HNUM/HDEN
      DHDX = (CK1/HDEN)**2
      DHDP = DHDX*DXDP

      DG = ONE - EXP2
      DGDP = ZERO
      IF (P14.GT.0.002_q) DGDP = -0.25_q*CA1*EXP2/(P*P14)

      IF (ALPHA.LT.ZERO) THEN
        EXP3 = EXP(-C1X*ALPHA/OMA)
        FA = EXP3
        DFADA = -C1X*EXP3/OMA2
      ELSEIF (ALPHA.LT.2.5_q) THEN
        DO IFX=1,7
          RTEMP(IFX) = RTEMP(IFX-1)*ALPHA
        ENDDO
        FA = DOT_PRODUCT(CREGU,RTEMP)
        DFADA = DOT_PRODUCT(CREGU1,RTEMP(0:6))
      ELSE
        EXP4 = EXP(C2X/OMA)
        FA = -CDX*EXP4
        DFADA = -CDX*EXP4*C2X/OMA2
      ENDIF

      FX = DG*(H + FA*(CH0-H))
      DFXDP = DGDP*(H + FA*(CH0-H)) + DG*DHDP*(ONE-FA)
      DFXDA = DG*(DFADA*(CH0-H))

! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=DFXDP*DPD+DFXDA*DALPHAD
      DFX_DD=DFXDP*DPDD+DFXDA*DALPHADD
      DFX_DTAU=DFXDA*DALPHADTAU

!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD=EXDLDA*FX+EXLDA*DFX_D
      VXDD=EXLDA*DFX_DD
      AMUXD=EXLDA*DFX_DTAU
      EX=EX+EXLDA*FX/TWO
      RETURN
      END SUBROUTINE VR2SCANx



      SUBROUTINE VR2SCANc(RU,RD,DRT,TAUU,TAUD,Ec,VCD1,VCD2,VCDD,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE ldalib, only : PW91LDA

      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(VR2SCANc)

      REAL(q) :: Ec,DRT,VCDD
      REAL(q) :: RU,TAUU,VCD1,AMUCD1
      REAL(q) :: RD,TAUD,VCD2,AMUCD2

      REAL(q) :: RT,Y,TAU
      REAL(q) :: TAUW,D_TAUW,D_TAUW_DD
      REAL(q) :: ZETA,OPZ,OMZ,D_ZETA_D1,D_ZETA_D2
      REAL(q) :: DS,D_DS_DZETA,D_DS_D1,D_DS_D2
      REAL(q) :: DX,D_DX_DZETA
      REAL(q) :: GC,D_GC_DZETA
      REAL(q) :: PHI,PHI2,PHI3,D_PHI_DZETA
      REAL(q) :: TAU0,D_TAU0_D1,D_TAU0_D2
      REAL(q) :: RS,RS12,D_RS
      REAL(q) :: P,D_P,D_P_DD
      REAL(q) :: AFIX_T,T2,D_T2_DP,D_T2_DRS,D_T2_DZETA
      REAL(q) :: ALPHA,ALPHA_DEN,D_ALPHA_D1,D_ALPHA_D2,D_ALPHA_DD,D_ALPHA_DT
      REAL(q) :: FA,D_FA_DA,D_FA_D1,D_FA_D2,D_FA_DD,D_FA_DT
      REAL(q) :: FACTOR0,D_FACTOR0_DRS,D2_FACTOR0_DRS
      REAL(q) :: ECLDA0,D_ECLDA0_DRS,D2_ECLDA0_DRS
      REAL(q) :: ECLDA1,D_ECLDA1_DRS,D_ECLDA1_DZETA,D2_ECLDA1_DRS,D2_ECLDA1_DRSDZETA
      REAL(q) :: DUMMY,ECTEMP,VCTEMPUP,VCTEMPDW
      REAL(q) :: BETA,D_BETA_DRS
      REAL(q) :: W1FAC,W1,D_W1_DRS,D_W1_DZETA
      REAL(q) :: A,D_A_DRS,D_A_DZETA
      REAL(q) :: V,D_V_DRS,D_V_DZETA,D_V_DP
      REAL(q) :: FACTOR1,EXP1,DELTAY,DELTAYP,D_DELTAYP_DRS,D_DELTAYP_DZETA,D_DELTAYP_DP
      REAL(q) :: GAUX,G,D_G_DGAUX,D_G_DRS,D_G_DZETA,D_G_DP
      REAL(q) :: ARG1,D_ARG1_DP,D_ARG1_DZETA,D_ARG1_DRS
      REAL(q) :: H1,D_H1_DRS,D_H1_DP,D_H1_DZETA
      REAL(q) :: EC1,D_EC1_DRS,D_EC1_DP,D_EC1_DZETA,D_EC1_D1,D_EC1_D2,D_EC1_DD
      REAL(q) :: W0,D_W0_DRS
      REAL(q) :: GINF,D_GINF_DP
      REAL(q) :: HCORE0,H0,D_H0_DRS,D_H0_DP
      REAL(q) :: EC0,D_EC0_DRS,D_EC0_DP,D_EC0_DZETA,D_EC0_D1,D_EC0_D2,D_EC0_DD
      REAL(q) :: EPP,D_EPP_D1,D_EPP_D2,D_EPP_DD,D_EPP_DT

      REAL(q) :: BETA_MB,GAMMA,FACTORA,FACTORB,B1C,B2C,B3C,XI
      REAL(q) :: A0,A1,A2,A3,A4,A5,A6,A7
      REAL(q) :: CFC1,CFC2,CFD1
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,PISQ
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: ETA,CKF2,DP24,C2

      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE=1.0_q)
      PARAMETER (TWO=2.0_q)
      PARAMETER (THREE=3.0_q)
      PARAMETER (FOUR=4.0_q)
      PARAMETER (FIVE=5.0_q)
      PARAMETER (SIX=6.0_q)
      PARAMETER (SEVEN=7.0_q)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1.0_q/3.0_q)
      PARAMETER (THRD2=2.0_q*THRD)
      PARAMETER (THRD4=4.0_q*THRD)
      PARAMETER (THRD5=1.0_q+THRD2)
      PARAMETER (THRD8=1.0_q+THRD5)
      PARAMETER (CFC1=0.640_q)
      PARAMETER (CFC2=1.500_q)
      PARAMETER (CFD1=0.700_q)
      parameter( a0 =  1.0_q )
      parameter( a1 = -0.64_q )
      parameter( a2 = -0.4352_q )
      parameter( a3 = -1.535685604549_q )
      parameter( a4 =  3.061560252175_q )
      parameter( a5 = -1.915710236206_q )
      parameter( a6 =  0.516884468372_q )
      parameter( a7 = -0.051848879792_q )
      PARAMETER(DP24=0.361_q**4)
      PARAMETER(B1C=0.0285764_q,B2C=0.0889_q,B3C=0.125541_q)
      PARAMETER(BETA_MB=0.06672455060314922_q)
      PARAMETER(FACTORA=0.1_q,FACTORB=0.1778_q)
      PARAMETER(XI=0.12802585262625815_q)
      PARAMETER(ETA=0.001_q)
      PARAMETER(C2=A1 + TWO*A2 + THREE*A3 + FOUR*A4 + FIVE*A5 + SIX*A6 + SEVEN*A7)
      PARAMETER(GAMMA=(ONE-DLOG(TWO))/PISQ)
      PARAMETER(CKF2=(THREE*PISQ)**THRD2)
GPACC !$ACC ROUTINE(PW91LDA) SEQ

      RT=RU+RD
      Y=DRT**2
      TAU=TAUU+TAUD

      TAUW=1.0d0/8.0d0*(Y/RT)
      D_TAUW = -TAUW/RT
      D_TAUW_DD = 0.25_q*DRT/RT

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-ONE),ONE)
      D_ZETA_D1= TWO*RD/RT**2
      D_ZETA_D2=-TWO*RU/RT**2

      OMZ = ONE - ZETA
      OPZ = ONE + ZETA

      DS=(OPZ**THRD5+OMZ**THRD5)/TWO
      D_DS_DZETA = THRD5*(OPZ**THRD2 - OMZ**THRD2)/TWO
      D_DS_D1=D_DS_DZETA*D_ZETA_D1
      D_DS_D2=D_DS_DZETA*D_ZETA_D2

      PHI = (OPZ**THRD2 + OMZ**THRD2)/TWO
      PHI2 = PHI*PHI
      PHI3 = PHI2*PHI
      D_PHI_DZETA = 0.0_q
      IF (OPZ.GT.1E-10_q) D_PHI_DZETA = THRD/OPZ**THRD
      IF (OMZ.GT.1E-10_q) D_PHI_DZETA = D_PHI_DZETA - THRD/OMZ**THRD

      DX = (OPZ**THRD4 + OMZ**THRD4)/TWO
      D_DX_DZETA = THRD2*(OPZ**THRD - OMZ**THRD)

      GC = (ONE - 2.363D0*(DX - ONE))*(ONE - ZETA**12)
      D_GC_DZETA = -12D0*(ONE - 2.363D0*(DX - ONE))*ZETA**11 - &
                 &  2.363D0*D_DX_DZETA*(ONE - ZETA**12)


      TAU0=0.30d0*(THREE*PISQ)**THRD2*RT**THRD5*DS
      D_TAU0_D1 = THRD5*TAU0/RT + TAU0*D_DS_D1/DS
      D_TAU0_D2 = THRD5*TAU0/RT + TAU0*D_DS_D2/DS

      RS    = (0.750d0/PI)**THRD/RT**THRD
      RS12  = DSQRT(RS)
      D_RS  = -THRD/RT*RS

      P = Y/(FOUR*CKF2*RT**THRD8)
      D_P = -THRD8*P/RT
      D_P_DD = DRT/(TWO*CKF2*RT**THRD8)

      AFIX_T = CKF2/FOUR**THRD4
      T2 = AFIX_T*P/(PHI2*RS)
      D_T2_DRS = -T2/RS
      D_T2_DP = AFIX_T/(PHI2*RS)
      D_T2_DZETA = -TWO*T2/PHI*D_PHI_DZETA

      ALPHA_DEN = ONE/(TAU0 + ETA*TAUW)
      ALPHA = (TAU-TAUW)*ALPHA_DEN

      D_ALPHA_D1 = -ALPHA_DEN*(D_TAUW + ALPHA*(D_TAU0_D1 + ETA*D_TAUW))
      D_ALPHA_D2 = -ALPHA_DEN*(D_TAUW + ALPHA*(D_TAU0_D2 + ETA*D_TAUW))
      D_ALPHA_DD = -D_TAUW_DD*ALPHA_DEN*(ONE + ALPHA*ETA)
      D_ALPHA_DT = ALPHA_DEN

      !--
      !FA
      !--
      IF (ALPHA.LT.0.0_q) THEN
        FA = DEXP(-CFC1*ALPHA/(ONE-ALPHA))
        D_FA_DA = -CFC1*FA/(ONE-ALPHA)**2
      ELSE IF (ALPHA.LT.2.5_q) THEN
        FA = A0 + ALPHA*(A1 + ALPHA*(A2 + ALPHA*(A3 + ALPHA*(A4 + ALPHA*(A5 + ALPHA*(A6 + ALPHA*A7))))))
        D_FA_DA = A1 + ALPHA*(TWO*A2 + ALPHA*(THREE*A3 + ALPHA*(FOUR*A4 + ALPHA*(FIVE*A5 + ALPHA*(SIX*A6 + ALPHA*SEVEN*A7)))))
      ELSE
        FA = -CFD1*DEXP(CFC2/(ONE-ALPHA))
        D_FA_DA = CFC2*FA/(ONE-ALPHA)**2
      END IF

      D_FA_D1 = D_FA_DA*D_ALPHA_D1
      D_FA_D2 = D_FA_DA*D_ALPHA_D2
      D_FA_DD = D_FA_DA*D_ALPHA_DD
      D_FA_DT = D_FA_DA*D_ALPHA_DT

      !--------------------------------------------------------------------------
      !EC1
      !--------------------------------------------------------------------------

      FACTOR0 = ONE + B2C*RS12 + B3C*RS
      D_FACTOR0_DRS = B2C/RS12/TWO + B3C
      D2_FACTOR0_DRS = -B2C/RS/RS12/FOUR

      ECLDA0 = -B1C/FACTOR0
      D_ECLDA0_DRS = D_FACTOR0_DRS*B1C/FACTOR0**2
      D2_ECLDA0_DRS = B1C/FACTOR0**2 * (D2_FACTOR0_DRS  - TWO*D_FACTOR0_DRS**2/FACTOR0)

      CALL PW91LDA(RS,ZETA,ECLDA1,D_ECLDA1_DRS,D_ECLDA1_DZETA,D2_ECLDA1_DRS,D2_ECLDA1_DRSDZETA)

      BETA = BETA_mb*(ONE + FACTORA*RS)/(ONE + FACTORB*RS)
      D_BETA_DRS = BETA_mB*(FACTORA - FACTORB)/(ONE + FACTORB*RS)**2

      W1FAC = ECLDA1/(GAMMA*PHI3)
      W1 = DEXP(-W1FAC) - ONE
      D_W1_DRS = -(W1+ONE)*D_ECLDA1_DRS/(GAMMA*PHI3)
      D_W1_DZETA = (W1+ONE)*(THREE*W1FAC*D_PHI_DZETA/PHI - D_ECLDA1_DZETA/(GAMMA*PHI3))

      A = BETA/(W1*GAMMA)
      D_A_DRS = D_BETA_DRS/(W1*GAMMA) - A*D_W1_DRS/W1
      D_A_DZETA = -A*D_W1_DZETA/W1

      V = A*T2
      D_V_DRS = D_A_DRS*T2 + A*D_T2_DRS
      D_V_DZETA = D_A_DZETA*T2 + A*D_T2_DZETA
      D_V_DP = A*D_T2_DP

      EXP1 = DEXP(-P**2/DP24)
      FACTOR1 = C2/(27.0_q*GAMMA*PHI3*W1*DS)
      DELTAY = FACTOR1*(20.0_q*RS*(D_ECLDA0_DRS - D_ECLDA1_DRS) - 45.0_q*ETA*(ECLDA0 - ECLDA1))
      DELTAYP = DELTAY*P*EXP1
      D_DELTAYP_DRS = -DELTAYP/W1*D_W1_DRS + FACTOR1*( (D_ECLDA0_DRS-D_ECLDA1_DRS)*(20.0_q - 45.0_q*ETA) + &
                    & 20.0_q*RS*(D2_ECLDA0_DRS - D2_ECLDA1_DRS))*P*EXP1
      D_DELTAYP_DZETA = -DELTAYP*(D_W1_DZETA/W1 + D_DS_DZETA/DS + THREE*D_PHI_DZETA/PHI) + &
                      & FACTOR1*(45.0_q*ETA*D_ECLDA1_DZETA - 20.0_q*RS*D2_ECLDA1_DRSDZETA)*P*EXP1
      D_DELTAYP_DP = DELTAY*EXP1 - TWO*P*DELTAYP/DP24

      GAUX = V - DELTAYP
      G = ONE/(ONE + FOUR*GAUX)**(0.25_q)
      D_G_DGAUX = -G**5
      D_G_DRS = D_G_DGAUX*(D_V_DRS - D_DELTAYP_DRS)
      D_G_DZETA = D_G_DGAUX*(D_V_DZETA - D_DELTAYP_DZETA)
      D_G_DP = D_G_DGAUX*(D_V_DP - D_DELTAYP_DP)

      ARG1 = ONE + W1*(ONE-G)
      D_ARG1_DP = -W1*D_G_DP
      D_ARG1_DRS = D_W1_DRS*(ONE-G) - W1*D_G_DRS
      D_ARG1_DZETA = D_W1_DZETA*(ONE-G) - W1*D_G_DZETA

      H1 = GAMMA*PHI3*DLOG(ARG1)
      D_H1_DRS = GAMMA*PHI3*D_ARG1_DRS/ARG1
      D_H1_DP = GAMMA*PHI3*D_ARG1_DP/ARG1
      D_H1_DZETA = THREE*H1*D_PHI_DZETA/PHI + GAMMA*PHI3*D_ARG1_DZETA/ARG1

      EC1 = ECLDA1 + H1
      D_EC1_DRS = D_ECLDA1_DRS + D_H1_DRS
      D_EC1_DP = D_H1_DP
      D_EC1_DZETA = D_ECLDA1_DZETA + D_H1_DZETA

      D_EC1_D1 = D_EC1_DRS*D_RS + D_EC1_DP*D_P + D_EC1_DZETA*D_ZETA_D1
      D_EC1_D2 = D_EC1_DRS*D_RS + D_EC1_DP*D_P + D_EC1_DZETA*D_ZETA_D2
      D_EC1_DD = D_EC1_DP*D_P_DD

      !-----------------------------------------------------------------------------
      ! EC0
      !-----------------------------------------------------------------------------
      W0 = DEXP(-ECLDA0/B1C) - ONE
      D_W0_DRS = -(W0+ONE)*D_ECLDA0_DRS/B1C

      GINF = ONE/(ONE + FOUR*XI*P)**(0.25_q)
      D_GINF_DP = -XI*GINF**5

      HCORE0 = ONE + W0*(ONE - GINF)
      H0 = B1C*dlog(HCORE0)
      D_H0_DRS  = B1C*D_W0_DRS*(ONE - GINF)/HCORE0
      D_H0_DP  = -B1C*W0*D_GINF_DP/HCORE0

      EC0 = (ECLDA0 + H0)*GC
      D_EC0_DRS = GC*(D_ECLDA0_DRS + D_H0_DRS)
      D_EC0_DP = GC*D_H0_DP
      D_EC0_DZETA = D_GC_DZETA*(ECLDA0 + H0)

      D_EC0_D1 = D_EC0_DRS*D_RS + D_EC0_DZETA*D_ZETA_D1 + D_EC0_DP*D_P
      D_EC0_D2 = D_EC0_DRS*D_RS + D_EC0_DZETA*D_ZETA_D2 + D_EC0_DP*D_P
      D_EC0_DD = D_EC0_DP*D_P_DD

      !-----------------------------------------------------------------------------

      EPP = EC1 + FA*(EC0 - EC1)
      EC=RT*EPP

      D_EPP_D1 = D_EC1_D1 + D_FA_D1*(EC0-EC1) + FA*(D_EC0_D1 - D_EC1_D1)
      D_EPP_D2 = D_EC1_D2 + D_FA_D2*(EC0-EC1) + FA*(D_EC0_D2 - D_EC1_D2)
      D_EPP_DD = D_EC1_DD + D_FA_DD*(EC0-EC1) + FA*(D_EC0_DD - D_EC1_DD)
      D_EPP_DT = D_FA_DT*(EC0-EC1)

      VCD1 = EPP + RT*D_EPP_D1
      VCD2 = EPP + RT*D_EPP_D2
      VCDD = RT*D_EPP_DD
      AMUCD1 = RT*D_EPP_DT
      AMUCD2 = AMUCD1

      RETURN
      END SUBROUTINE VR2SCANc

!************************ SUBROUTINE VSCAN2x ****************************
!
!> calculates the first order derivatives of Ex wrt n and |grad(n)|
!> Perdew et. al. PRL (2009)
!>
!> Originally written by Jianwei Sun
!>
!> Re-written and generalised by James W. Furness and Aaron D. Kaplan,
!> 5 August 2020,
!> to accomodate rSCAN, r++SCAN, r2SCAN, and r4SCAN
!> JWF: jfurness@tulane.edu
!> ADK: kaplan@temple.edu
!> r2 SCAN preprint (arXiv:2008.03374)
!>
!> New flags for METAGGA tag: RSCAN, RPPSCAN, R2SCAN, and R4SCAN (respectively)
!>
!> r4SCAN is still a work in progress
!> r++SCAN is unpublished, and was a test of some of the principles in the
!> r2SCAN preprint
!>
!> everything in Hartree units
!>
!> Funtional selection uses 3 integer flags:
!>   IALPHA  - controls iso-orbital indcator
!>           (0 = alpha, 1 = alpha', 2 = \bar{alpha} )
!>   IINTERP - controls interpolation function
!>           (0 = SCAN, 1 = rSCAN)
!>   IDELFX  - controls gradient expansion corrections
!>           (0 = SCAN, 1 = rSCAN (2nd order), 2 = rSCAN (4th order), 3 = r2SCAN with diff pars)
!>
!> Hence for (IALPHA, IINTERP, IDELFX):
!>     SCAN = (0, 0, 0)
!>     rSCAN = (1, 1, 0)
!>     r++SCAN = (2, 1, 0)
!>     r2SCAN = (2, 1, 1)
!>     r4SCAN = (2, 1, 2)
!>
!> ATTENTION: Every values are passed "as they are", i.e. including
!> possibly unphysical numerical errors (e.g. negative charge densities)
!> values need to be checked accordingly
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> TAUU,TAUD                    kinetic energy density up/down
!> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
!> VXD1 VXD2                    THE DERIVATIVES OF EX WRT n
!> VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n|
!> AMUXD1, AMUXD2               THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************



      SUBROUTINE VSCAN2x(&
      &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, &
      &   Ex_SCAN,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2,ID_XC)
GPACC !$ACC ROUTINE SEQ

      !      USE prec
      !      USE constant
      USE xc_name
      IMPLICIT none
GPOMP  !$OMP DECLARE TARGET(VSCAN2x)

      INTEGER, PARAMETER :: q=SELECTED_REAL_KIND(10)
      integer,intent(in) :: ID_XC
      integer :: IALPHA, IINTERP, IDELFX
      real(q),intent(in) :: ru,rd,dru,drd,drt,tauu,taud
      real(q),intent(out) :: VXD1,VXD2,VXDD1,VXDD2,AMUXD1,AMUXD2
      real(q),intent(out) :: Ex_SCAN
      real(q) :: ex_up,ex_down
      real(q),parameter :: two=2._q
GPACC !$ACC ROUTINE(eps_SCANx) SEQ

      VXD1=0._q;VXD2=0._q;
      VXDD1=0._q;VXDD2=0._q;
      AMUXD1=0._q;AMUXD2=0._q;
      Ex_SCAN=0._q
      ex_up = 0._q
      ex_down=0._q

      if (ID_XC==ID_XC_SCANL) then ! SCAN
        IALPHA=0; IINTERP=0; IDELFX=0;
      else if (ID_XC==ID_XC_RSCANL) then ! rSCAN
        IALPHA=1; IINTERP=1; IDELFX=0;
      else if ((ID_XC==ID_XC_R2SCANL).OR.(ID_XC==ID_XC_OFR2)) then !r2SCAN
        IALPHA=2; IINTERP=1; IDELFX=1;
      endif

      call eps_SCANx(two*ru,two*dru,two*tauu,ex_up, VXD1, VXDD1, AMUXD1,&
      &    IALPHA, IINTERP, IDELFX) ! spin-up

      call eps_SCANx(two*rd,two*drd,two*taud,ex_down, VXD2, VXDD2, AMUXD2,&
      &    IALPHA, IINTERP, IDELFX) ! spin-down

      Ex_SCAN = (ex_up+ex_down)/two

      RETURN
      END SUBROUTINE VSCAN2x


      subroutine eps_SCANx(den, grd, tau, eps_x, dedd, dedg, dedt, &
         &                 IALPHA, IINTERP, IDELFX)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(eps_SCANx)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: den, grd, tau
          real(kdp), intent(out) :: eps_x, dedd, dedg, dedt
          integer, intent(in) :: IALPHA, IINTERP, IDELFX

          real(kdp) :: cp, den83, den53
          real(kdp) :: exlda, dexldadd, fx, dfxdp, dfxda
          real(kdp) :: p, dpdd, dpdg
          real(kdp) :: alpha, dadd, dadg, dadt, reg, dreg
          real(kdp) :: tueg, dtuegdd, tueg_con
          real(kdp) :: tauw, dtauwdd, dtauwdg

          real(kdp), parameter :: AX = -0.7385587663820224058842300326808360_kdp
          real(kdp), parameter :: PI = 3.1415926535897932384626433832795_kdp
          real(kdp), parameter :: PI2 = PI*PI
          real(kdp), parameter :: ETA = 1.0e-3_kdp
          real(kdp), parameter :: TAU_R = 1.0e-4_kdp
          real(kdp), parameter :: A_REG = 1.0e-3_kdp
GPACC !$ACC ROUTINE(exchange_enhancement) SEQ

      !       Reduced density gradient [FD]
          den83 = den**(8._kdp/3._kdp)
          cp = 4._kdp*(3._kdp*PI2)**(2._kdp/3._kdp)
          p = grd**2/(cp*den83)
          dpdd = -8._kdp/3._kdp*p/den
          dpdg = 2._kdp*p/grd

      !       Regularised Alpha [FD]
          tueg_con = 3._kdp/10._kdp*(3._kdp*PI2)**(2._kdp/3._kdp)
          den53 = den**(5._kdp/3._kdp)
          if (IALPHA .eq. 1) then ! regularised tau_ueg for rSCAN \tilde{alpha}
              tueg = tueg_con*den53 + TAU_R
          else                    ! Unregularised tau_ueg
              tueg = tueg_con*den53
          end if
          dtuegdd = 5._kdp/3._kdp*tueg_con*den53/den

          tauw = grd**2/8._kdp/den
          dtauwdd = -grd**2/(8._kdp*den**2)
          dtauwdg = 2._kdp*grd/(8._kdp*den)

          if (IALPHA .eq. 0) then ! SCAN unregularised alpha
              alpha = (tau - tauw)/tueg
              dadd = -(tau - tauw)*dtuegdd/tueg**2 - dtauwdd/tueg
              dadg = -dtauwdg/tueg
              dadt = 1._kdp/tueg

          else if (IALPHA .eq. 1) then ! rSCAN alpha'
              alpha = (tau - tauw)/tueg
              dadd = ((tauw - tau)*dtuegdd - tueg*dtauwdd)/tueg**2
              dadg = -grd/(4._kdp*den*tueg)
              dadt = 1._kdp/tueg

              if (A_REG .gt. 0._kdp) then
                  reg = alpha**3/(alpha**2 + A_REG)
                  dreg = (alpha**4 + 3._kdp*alpha**2*A_REG)/(alpha**2 + A_REG)**2
                  dadd = dadd*dreg
                  dadg = dadg*dreg
                  dadt = dadt*dreg
                  alpha = reg
              endif

          else if (IALPHA .eq. 2) then ! \bar{alpha} for r2scan and r4scan
              alpha = (tau - tauw)/(tueg + ETA*tauw)
              dadd = -dtauwdd/(tueg + ETA*tauw) &
             &    - (tau - tauw)*(dtuegdd + ETA*dtauwdd)/(tueg + ETA*tauw)**2
              dadg = -ETA*(tau - tauw)*dtauwdg/(tueg + ETA*tauw)**2 &
             &    - dtauwdg/(tueg + ETA*tauw)
              dadt = 1._kdp/(tueg + ETA*tauw)

          else
__CPU         call quit('ERROR: Unknown IALPHA in SCAN')
          end if

      !       UEG exchange density [FD]
          exlda = AX*den**(4._kdp/3._kdp)
          dexldadd = 4._kdp*AX*den**(1._kdp/3._kdp)/3._kdp

          call exchange_enhancement(p, alpha, fx, dfxdp, dfxda, ETA, IINTERP, IDELFX)

      !       Exchange Energy[FD]
          eps_x = exlda*fx
          dedd = dexldadd*fx + exlda*(dfxdp*dpdd + dfxda*dadd)
          dedg = exlda*(dfxdp*dpdg + dfxda*dadg)
          dedt = exlda*dfxda*dadt

      end subroutine eps_SCANx

      subroutine exchange_enhancement(p, alpha, Fx, dfxdp, dfxda, ETA, IINTERP, IDELFX)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(exchange_enhancement)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: p, alpha
          real(kdp), intent(out) :: Fx, dfxdp, dfxda
          real(kdp), intent(in) :: ETA
          integer, intent(in) :: IINTERP, IDELFX


          real(kdp) :: ief, diefda, oma
          real(kdp) :: h0x
          real(kdp) :: h1x, dh1xdp, dh1xda
          real(kdp) :: gx, dgxdp
          real(kdp) :: del_f2, C2, ALPHA_GE, damp, ddampdp
          real(kdp) :: del_fx, ddel_fxdp, ddel_fxda
          real(kdp) :: wfac, dwfacdp, vfac, dvfacdp, dvfacda
          real(kdp) :: yfac, dyfacdp, dyfacda,D_DAMP2

          real(kdp), parameter :: A1 = 4.9479_kdp
          real(kdp), parameter :: K0 = 0.174_kdp
          real(kdp), parameter :: MU = 10._kdp/81._kdp
          real(kdp), parameter :: K1 = 0.065_kdp

          real(kdp), parameter :: cfx1 = 0.667_kdp
          real(kdp), parameter :: cfx2 = 0.800_kdp
          real(kdp), parameter :: cfdx1 = 1.24_kdp

          real(kdp), parameter :: B1 = 0.156632_kdp
          real(kdp), parameter :: B2 = 0.12083_kdp
          real(kdp), parameter :: B3 = 0.5_kdp
          real(kdp), parameter :: B4 = MU*MU/K1 - 0.112654_kdp

          real(kdp), dimension(8) :: PARAMS
          integer, dimension(8) :: f_x_e
          !real(kdp), parameter :: D_DAMP2 = 0.361_kdp
GPACC !$ACC ROUTINE(get_del_fx) SEQ

          PARAMS = (/ &
          &      -0.023185843322_kdp,0.234528941479_kdp,-0.887998041597_kdp, &
          &      1.451297044490_kdp,-0.663086601049_kdp,-0.4445555_kdp,-0.667_kdp, &
          &      1._kdp/)

          f_x_e = (/7,6,5,4,3,2,1,0/)

          if (IDELFX==3) then
            D_DAMP2 = 5.259688_kdp
          else
            D_DAMP2 = 0.361_kdp
          end if

          ALPHA_GE = 20._kdp/27._kdp + ETA*5._kdp/3._kdp

          ief = 0._kdp
          diefda = 0._kdp
          oma = 1._kdp - alpha
          if (IINTERP .eq. 0) then  ! scan interpolation function
              if (alpha .lt. 1._kdp) then
                  ief = exp(-cfx1*alpha/oma)
                  diefda = -cfx1*exp(-cfx1*alpha/oma)/oma**2
              else
                  ief = -cfdx1*exp(cfx2/oma)
                  diefda = -cfx2*cfdx1*exp(cfx2/oma)/oma**2
              endif

          else if (IINTERP .eq. 1) then ! rscan, r2scan, r4scan interpolation function
              if (alpha .lt. 1.0e-13_kdp) then
                  ief = exp(-cfx1*alpha/oma)
                  diefda = -cfx1*exp(-cfx1*alpha/oma)/oma**2
              else if( alpha .lt. 2.5_kdp) then
                  ief = dot_product(alpha**f_x_e, PARAMS)
                  diefda = dot_product(alpha**f_x_e(2:), f_x_e(:7)*PARAMS(:7))
              else if (alpha .ge. 2.5_kdp) then
                  ief = -cfdx1*exp(cfx2/oma)
                  diefda = -cfx2*cfdx1*exp(cfx2/oma)/oma**2
              endif

          else
__CPU         call quit('ERROR: Unknown IINTERP in SCAN')
          end if

      !       Single orbital enhancement
          h0x = 1._kdp + K0

      !       Slowly varying enhancement
          if (IDELFX .eq. 0) then     ! scan, rscan
              wfac = B4*p**2*exp(-B4*p/MU)
              dwfacdp = B4*p*exp(-B4*p/MU)*(2._kdp - B4*p/MU)

              vfac = B1*p + B2*oma*exp(-B3*oma**2)
              yfac = MU*p + wfac + vfac**2
              h1x = 1._kdp + K1 - K1/(1._kdp + yfac/K1)

              dvfacdp = B1
              dyfacdp = MU + dwfacdp + 2._kdp*vfac*dvfacdp
              dh1xdp = dyfacdp/(1._kdp + yfac/K1)**2

              dvfacda = -B2*(1._kdp - 2._kdp*B3*oma**2)*exp(-B3*oma**2)
              dyfacda = 2._kdp*vfac*dvfacda
              dh1xda = dyfacda/(1._kdp + yfac/K1)**2

          else if (IDELFX .eq. 1 .or. IDELFX .eq. 2 .or. IDELFX .eq. 3) then  ! Second order corrections
              del_f2 = dot_product(f_x_e(:7), PARAMS(:7))
              C2 = -del_f2*(1._kdp - h0x)

      !       Damping
              damp = exp(-p**2/D_DAMP2**4)
              ddampdp = -2._kdp*damp*p/D_DAMP2**4

      !       Slowly varying contribution [FD]
              h1x = 1._kdp + K1 - K1/(1._kdp + p*(MU + ALPHA_GE*C2*damp)/K1)
              dh1xdp = K1**2*(MU + ALPHA_GE*C2*(damp + p*ddampdp)) &
              &              /(K1 + MU*p + ALPHA_GE*C2*p*damp)**2
              dh1xda = 0._kdp

          else
__CPU         call quit('ERROR: Unknown IDELFX in SCAN')
          end if

          gx = 1._kdp - exp(-A1/p**(1._kdp/4._kdp))
          dgxdp = -A1*exp(-A1/p**(1._kdp/4._kdp))/(4._kdp*p**(5._kdp/4._kdp))

          if (IDELFX .eq. 2) then  ! 4th order corrections for r4scan
              call get_del_fx(p, alpha, del_fx, ddel_fxdp, ddel_fxda, &
              &                        K0, K1, C2, ETA, ALPHA_GE, PARAMS, f_x_e)
          else
              del_fx = 0._kdp
              ddel_fxdp = 0._kdp
              ddel_fxda = 0._kdp
          end if

          fx = (h1x + ief*(h0x - h1x) + del_fx)*gx
          dfxdp = (del_fx + h1x + (h0x - h1x)*ief)*dgxdp &
          &       + gx*(ddel_fxdp + dh1xdp - ief*dh1xdp)
          dfxda = gx*((h0x - h1x)*diefda + ddel_fxda + dh1xda - ief*dh1xda)

      end subroutine exchange_enhancement

      subroutine get_del_fx(p, alpha, del_fx, ddel_fxdp, ddel_fxda, &
      &                        K0, K1, C2, ETA, ALPHA_GE, PARAMS, f_x_e)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(get_del_fx)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: p, alpha, K0, K1, C2, ETA, ALPHA_GE
          real(kdp), dimension(8), intent(in) :: PARAMS
          integer, dimension(8), intent(in) :: f_x_e
          real(kdp), intent(out) :: del_fx, ddel_fxdp, ddel_fxda
          real(kdp) :: order_1, dorder_1dp, dorder_1da, C_pa, C_aa, C_pp
          real(kdp) :: oma, damp, t1, dt1dp, dt1da, ddampdp, ddampda
GPACC !$ACC ROUTINE(get_dx_terms) SEQ
GPACC !$ACC ROUTINE(get_fourth_order_damp) SEQ

          call get_dx_terms(C_aa, C_pa, C_pp, &
          &              K0, K1, C2, ETA, PARAMS, f_x_e)

          oma = 1._kdp - alpha

          order_1 = C2*(oma - ALPHA_GE*p)
          dorder_1dp = -C2*ALPHA_GE
          dorder_1da = -C2

      !       Correcting contribution [FD]
          t1 = order_1 + C_aa*oma**2 + C_pa*p*oma + C_pp*p**2
          dt1dp = dorder_1dp + 2*C_pp*p + C_pa*oma
          dt1da = dorder_1da - 2*C_aa*oma - C_pa*p

          call get_fourth_order_damp(p, alpha, damp, ddampdp, ddampda)

          del_fx = t1*damp
          ddel_fxdp = damp*dt1dp + t1*ddampdp
          ddel_fxda = t1*ddampda + dt1da*damp

      end subroutine get_del_fx

      subroutine get_dx_terms(C_aa, C_pa, C_pp, K0, K1, C2, ETA, &
      &                  PARAMS, f_x_e)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(get_dx_terms)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: K0, K1, C2, ETA
          real(kdp), dimension(8), intent(in) :: PARAMS
          integer, dimension(8), intent(in) :: f_x_e
          real(kdp), intent(out) :: C_aa, C_pa, C_pp
          real(kdp) :: h0x, eta_term, del_f2, del_f4
          real(kdp) :: ALPHA_GE

          real(kdp), parameter :: MU = 10._kdp/81._kdp

          ALPHA_GE = 20._kdp/27._kdp + ETA*5._kdp/3._kdp

          eta_term = ETA*3._kdp/4._kdp + 2._kdp/3._kdp
          h0x = 1._kdp + K0

          del_f2 = dot_product(f_x_e(:7), PARAMS(:7))
          del_f4 = dot_product(f_x_e(:7)*(f_x_e(:7) - 1._kdp), PARAMS(:7))

          C_aa = 73._kdp/5000._kdp - 0.5_kdp*del_f4*(h0x - 1._kdp)

          C_pa = 511._kdp/13500._kdp - 73._kdp/1500._kdp*ETA &
          &            - del_f2*(ALPHA_GE*C2 + MU)

          C_pp = 146._kdp/2025._kdp*eta_term**2 - 73._kdp/405._kdp*eta_term &
          &            + (ALPHA_GE*C2 + MU)**2/K1

      end subroutine get_dx_terms

      subroutine get_fourth_order_damp(p, alpha, damp4, &
      &            ddamp4dp, ddamp4da)
GPACC !$ACC ROUTINE SEQ

          implicit NONE
GPOMP  !$OMP DECLARE TARGET(get_fourth_order_damp)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: p, alpha
          real(kdp), intent(out) :: damp4, ddamp4dp, ddamp4da
          real(kdp) :: t1, t2, dt1da, dt2dp, dt2da, oma
          real(kdp), parameter :: DX_DAMP4_P = 0.802_kdp
          real(kdp), parameter :: DX_DAMP4_A = 0.178_kdp

          oma = 1._kdp - alpha

          t1 = 2._kdp*alpha**2/(1._kdp + alpha**4)
          dt1da = -4._kdp*alpha*(alpha**4 - 1._kdp)/(1._kdp + alpha**4)**2

          t2 = exp(-oma**2/DX_DAMP4_A**2 - p**2/DX_DAMP4_P**4)
          dt2dp = -2.d0*t2*p/DX_DAMP4_P**4
          dt2da = 2*oma*t2/DX_DAMP4_A**2

          damp4 = t1*t2
          ddamp4dp = t1*dt2dp
          ddamp4da = t1*dt2da + t2*dt1da

      end subroutine get_fourth_order_damp


      !************************ SUBROUTINE VSCAN2c ****************************
      !! Originally written by Jianwei Sun
      !
      !> Re-written and generalised by James W. Furness and Aaron D. Kaplan,
      !> 5 August 2020,
      !> to accomodate rSCAN, r++SCAN, r2SCAN, and r4SCAN
      !> JWF: jfurness@tulane.edu
      !> ADK: kaplan@temple.edu
      !> r2 SCAN preprint (arXiv:2008.03374)
      !>
      !> New flags for METAGGA tag: RSCAN, RPPSCAN, R2SCAN, and R4SCAN (respectively)
      !>
      !> r4SCAN is still a work in progress
      !> r++SCAN is unpublished, and was a test of some of the principles in the
      !> r2SCAN preprint
      !>
      !> everything in Hartree units
      !>
      !> Funtional selection uses 3 integer flags:
      !>   IALPHA  - controls iso-orbital indcator
      !>           (0 = alpha, 1 = alpha', 2 = \bar{alpha} )
      !>   IINTERP - controls interpolation function
      !>           (0 = SCAN, 1 = rSCAN)
      !>   IDELEC  - controls gradient expansion corrections
      !>           (0 = SCAN, 1 or 2 = rSCAN (2nd order), 3 = r2SCAN with different params)
      !>
      !> Hence for (IALPHA, IINTERP, IDELEC):
      !>     SCAN = (0, 0, 0)
      !>     rSCAN = (1, 1, 0)
      !>     r++SCAN = (2, 1, 0)
      !>     r2SCAN = (2, 1, 1)
      !>     r4SCAN = (2, 1, 1) or (2, 1, 2)
      !>
      !> RU,RD                        density up,down
      !> DRU, DRD                     abs. val. gradient of density up/down
      !> DRT                          abs. val. gradient of total density
      !> TAUU,TAUD                    kinetic energy density up/down
      !> TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
      !> VCD1 VCD2                    THE DERIVATIVES OF EC WRT n
      !> VCDD1,VCDD2                  THE DERIVATIVES OF EC WRT |grad n|
      !> AMUCD1, AMUCD2                   THE DERIVATIVES OF EC WRT TAU
      !
      !***********************************************************************

      SUBROUTINE VSCAN2c(&
          &   RU,RD,DRU,DRD,DRT,TAUU,TAUD,&
          &   Ec_SCAN,VCD1,VCD2,VCDD1,VCDD2,VCDD,AMUCD1,AMUCD2,ID_XC)
GPACC !$ACC ROUTINE SEQ

          USE xc_name
          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(VSCAN2c) 

          integer,intent(in) :: ID_XC
          real*8,intent(in) :: ru,rd,dru,drd,drt,tauu,taud
          real*8,intent(out) :: Ec_SCAN,VCD1,VCD2,VCDD1,VCDD2,VCDD,AMUCD1,AMUCD2
          integer :: IALPHA, IINTERP, IDELEC
GPACC !$ACC ROUTINE(scan_c_new) SEQ

          if (ID_XC==ID_XC_SCANL) then ! SCAN
            IALPHA=0; IINTERP=0; IDELEC=0;
          else if (ID_XC==ID_XC_RSCANL) then ! rSCAN
            IALPHA=1; IINTERP=1; IDELEC=0;
          else if ((ID_XC==ID_XC_R2SCANL).OR.(ID_XC==ID_XC_OFR2)) then !r2SCAN
            IALPHA=2; IINTERP=1; IDELEC=1;
          endif

          call scan_c_new(ru, rd, dru, drd, drt, tauu, taud, ec_scan, &
              &             vcd1, vcd2, vcdd1, vcdd2, vcdd, amucd1, amucd2, &
              &             IALPHA, IINTERP, IDELEC)

      end subroutine VSCAN2c

      subroutine scan_c_new(d0, d1, g0, g1, gt, t0, t1, eps_c, &
          &             dedd0, dedd1, dedg0, dedg1, dedg, dedt0, dedt1, &
          &             IALPHA, IINTERP, IDELEC)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(scan_c_new)
          integer, parameter :: kdp = selected_real_kind(15, 307)
          real(kdp), intent(in) :: d0, d1, g0, g1, gt, t0, t1
          real(kdp), intent(out) :: eps_c, dedd0, dedt0, dedd1, dedt1,dedg0,dedg1,dedg
          integer, intent(in) :: IALPHA, IINTERP, IDELEC

          real(kdp) :: dt, tt, dthrd, rs, drsdd0, drsdd1
          real(kdp) :: zeta, dzetadd0, dzetadd1, ds_z, dds_zdd0, dds_zdd1
          real(kdp) :: s, dsdd0, dsdg0, dsdd1, dsdg1, gngn0, gngn1
          real(kdp) :: y0, y1, y, yc, dycdg0, dycdg1
          real(kdp) :: alpha, dadd0, dadg0, dadt0, dadd1, dadg1, dadt1, oma
          real(kdp) :: reg, dreg
          real(kdp) :: tueg, tueg_con, dtuegdd0, dtuegdd1
          real(kdp) :: tauw, dtauwdd0, dtauwdg0, dtauwdd1, dtauwdg1
          real(kdp) :: ief, diefda, diefdd0, diefdd1, diefdg0, diefdg1
          real(kdp) :: diefdt0, diefdt1
          real(kdp) :: ec0, dec0drs, dec0ds, dec0dz
          real(kdp) :: dec0dd0, dec0dd1, dec0dg0, dec0dg1
          real(kdp) :: ec1, dec1drs, dec1ds, dec1dz ! got rid of dec1da, as this was not assigned, ADK
          real(kdp) :: dec1dd0, dec1dd1, dec1dg0, dec1dg1
          real(kdp) :: dsdg,dtauwdg,dadg,diefdg,dec0dg,dec1dg ! corrections for spin-polarized calculations, ADK

          real(kdp), parameter :: CFDC1 = 0.7_kdp
          real(kdp), parameter :: CFC1 = 0.640_kdp
          real(kdp), parameter :: CFC2 = 1.5_kdp

          real(kdp), dimension(7) :: IE_PARAMS_C

          real(kdp), parameter :: ETA = 1.0e-3_kdp
          real(kdp), parameter :: TAU_R = 1.0e-4_kdp
          real(kdp), parameter :: A_REG = 1.0e-3_kdp

          integer :: i

          real(kdp), parameter :: BETA_MB = 0.066725_kdp
          real(kdp), parameter :: GAMMA = 0.031090690869655_kdp
          real(kdp), parameter :: AFACTOR = 0.1_kdp
          real(kdp), parameter :: BFACTOR = 0.1778_kdp
          real(kdp), parameter :: PI = 3.1415926535897932384626433832795_kdp
          real(kdp), parameter :: PI2 = PI*PI
GPACC !$ACC ROUTINE(scan_ec0) SEQ
GPACC !$ACC ROUTINE(scan_ec1) SEQ

          IE_PARAMS_C = &
          &    (/-0.64_kdp, -0.4352_kdp, -1.535685604549_kdp, 3.061560252175_kdp, &
          &        -1.915710236206_kdp, 0.516884468372_kdp, -0.051848879792_kdp/)

          dt = d0 + d1
          y0 = g0**2
          y1 = g1**2
          y = gt**2
          yc = (y - y0 - y1)/2._kdp
          tt = t0 + t1

          dycdg0 = yc/g0
          dycdg1 = yc/g1

          !        Zeta [FD]
          zeta = min(max((d0 - d1)/dt, -0.99999999999990_kdp), 0.99999999999990_kdp)
          dzetadd0 = 2._kdp*d1/dt**2
          dzetadd1 = -2._kdp*d0/dt**2

          !        Wigner-Seitz radius [FD]
          dthrd = dt**(1._kdp/3._kdp)
          rs = (0.75_kdp/PI)**(1._kdp/3._kdp)/dthrd
          drsdd0 = -(0.75_kdp/PI)**(1._kdp/3._kdp)/(3._kdp*dthrd*dt)
          drsdd1 = drsdd0

          !        Reduced density gradient [FD]
          s = gt/(2._kdp*(3._kdp*PI2)**(1._kdp/3._kdp)*dt**(4._kdp/3._kdp))
          dsdd0 = -4._kdp/3._kdp*s/dt
          dsdd1 = dsdd0
          gngn0 = 1._kdp/(2._kdp*gt)*(2._kdp*g0+2._kdp*yc/g0)
          gngn1 = 1._kdp/(2._kdp*gt)*(2._kdp*g1+2._kdp*yc/g1)
          dsdg0 = s/gt*gngn0
          dsdg1 = s/gt*gngn1
          dsdg = s/gt ! ADK

          !        ds_zeta [FD]
          ds_z = ((1._kdp + zeta)**(5._kdp/3._kdp) + (1._kdp - zeta)**(5._kdp/3._kdp))/2._kdp
          dds_zdd0 = 5._kdp/3._kdp* &
          &    ((1._kdp + zeta)**(2._kdp/3._kdp) - (1._kdp - zeta)**(2._kdp/3._kdp))*dzetadd0/2._kdp
          dds_zdd1 = 5._kdp/3._kdp* &
          &    ((1._kdp + zeta)**(2._kdp/3._kdp) - (1._kdp - zeta)**(2._kdp/3._kdp))*dzetadd1/2._kdp

          !       alpha
          tueg_con = 3._kdp/10._kdp*(3._kdp*PI2)**(2._kdp/3._kdp)

          if (IALPHA .eq. 1) then
              tueg = (tueg_con*dt**(5._kdp/3._kdp) + TAU_R)*ds_z
              dtuegdd0 = 5._kdp/3._kdp*tueg_con*dt**(2._kdp/3._kdp)*ds_z &
              &           + tueg*dds_zdd0/ds_z
              dtuegdd1 = 5._kdp/3._kdp*tueg_con*dt**(2._kdp/3._kdp)*ds_z &
              &           + tueg*dds_zdd1/ds_z
          else
              tueg = tueg_con*dt**(5._kdp/3._kdp)*ds_z
              dtuegdd0 = 5._kdp/3._kdp*tueg/dt + tueg*dds_zdd0/ds_z
              dtuegdd1 = 5._kdp/3._kdp*tueg/dt + tueg*dds_zdd1/ds_z
          end if

          !tauw = min(y/(8._kdp*dt), tt)
          tauw = y/(8._kdp*dt)
          dtauwdd0 = -tauw/dt
          dtauwdd1 = dtauwdd0
          dtauwdg0 = (g0 + dycdg0)/(4._kdp*dt)
          dtauwdg1 = (g1 + dycdg1)/(4._kdp*dt)
          dtauwdg = gt/(4._kdp*dt) ! ADK

          if (IALPHA .eq. 0 .or. IALPHA .eq. 1) then
              ! Conventional alpha
              alpha = (tt - tauw)/tueg
              dadd0 = -(dtauwdd0 + alpha*dtuegdd0)/tueg
              dadd1 = -(dtauwdd1 + alpha*dtuegdd1)/tueg
              dadg0 = -dtauwdg0/tueg
              dadg1 = -dtauwdg1/tueg
              dadg = -dtauwdg/tueg ! ADK
              dadt0 = 1._kdp/tueg
              dadt1 = 1._kdp/tueg

          else if (IALPHA .eq. 2) then
              ! regularised alpha of r++scan, r2scan, r4scan
              alpha = (tt - tauw)/(tueg + ETA*tauw)
              dadd0 = -dtauwdd0/(tueg + ETA*tauw) &
              &    - (tt - tauw)*(dtuegdd0 + ETA*dtauwdd0)/(tueg + ETA*tauw)**2
              dadd1 = -dtauwdd1/(tueg + ETA*tauw) &
              &    - (tt - tauw)*(dtuegdd1 + ETA*dtauwdd1)/(tueg + ETA*tauw)**2
              dadg0 = -ETA*(tt - tauw)*dtauwdg0/(tueg + ETA*tauw)**2 &
              &    - dtauwdg0/(tueg + ETA*tauw)
              dadg1 = -ETA*(tt - tauw)*dtauwdg1/(tueg + ETA*tauw)**2 &
              &    - dtauwdg1/(tueg + ETA*tauw)
              dadg =  -ETA*(tt - tauw)*dtauwdg/(tueg + ETA*tauw)**2 &
              &    - dtauwdg/(tueg + ETA*tauw) ! ADK
              dadt0 = 1._kdp/(tueg + ETA*tauw)
              dadt1 = dadt0
          else
__CPU         call quit('Bad IALPHA in SCAN')
          end if

          if (IALPHA .eq. 1) then
              ! alpha regularisation for rscan
              reg = alpha**3/(alpha**2 + A_REG)
              dreg = (alpha**4 + 3._kdp*alpha**2*A_REG)/(alpha**2 + A_REG)**2
              dadd0 = dadd0*dreg
              dadd1 = dadd1*dreg
              dadg0 = dadg0*dreg
              dadg1 = dadg1*dreg
              dadg = dadg*dreg ! ADK
              dadt0 = dadt0*dreg
              dadt1 = dadt1*dreg
              alpha = reg
          end if

          ief = 0._kdp
          diefda = 0._kdp
          oma = 1._kdp - alpha

          if (IINTERP .eq. 0) then
              ! SCAN interpolation function
              if (alpha .lt. 1._kdp) then
                  ief = exp(-CFC1*alpha/oma)
                  diefda = -CFC1*ief/oma**2
              else if (alpha .ge. 1.0_kdp) then
                  ief = -CFDC1*exp(CFC2/oma)
                  diefda = CFC2*ief/oma**2
              endif

          else if (IINTERP .eq. 1) then
              ! rSCAN interpolation function
              if (alpha .lt. 1.0e-13_kdp) then
                  ief = exp(-CFC1*alpha/oma)
                  diefda = -CFC1*ief/oma**2
              else if (alpha .lt. 2.5_kdp) then
                  ief = 1._kdp
                  do i = 1, 7
                      ief = ief + IE_PARAMS_C(i)*alpha**(i)
                      diefda = diefda + i*IE_PARAMS_C(i)*alpha**(i - 1)
                  end do
              else if (alpha .ge. 2.5_kdp) then
                  ief = -CFDC1*exp(CFC2/oma)
                  diefda = CFC2*ief/oma**2
              endif
          else
__CPU         call quit('Bad IINTERP in SCAN')
          end if

          diefdd0 = diefda*dadd0
          diefdd1 = diefda*dadd1
          diefdg0 = diefda*dadg0
          diefdg1 = diefda*dadg1
          diefdg = diefda*dadg ! ADK
          diefdt0 = diefda*dadt0
          diefdt1 = diefda*dadt1

          !        Single Orbital Correlation
          ec0 = 0._kdp
          dec0drs = 0._kdp
          dec0ds = 0._kdp
          dec0dz = 0._kdp
          call scan_ec0(rs, s, zeta, ec0, dec0drs, dec0ds, dec0dz)
          dec0dd0 = dec0drs*drsdd0 + dec0dz*dzetadd0 + dec0ds*dsdd0
          dec0dd1 = dec0drs*drsdd1 + dec0dz*dzetadd1 + dec0ds*dsdd1
          dec0dg0 = dec0ds*dsdg0
          dec0dg1 = dec0ds*dsdg1
          dec0dg = dec0ds*dsdg ! ADK

          !        Slowly Varying Correlation
          ec1 = 0._kdp
          dec1drs = 0._kdp
          dec1ds = 0._kdp
          dec1dz = 0._kdp
          call scan_ec1(rs, s, zeta, &
          &                  ec1, dec1drs, dec1dz, dec1ds, IE_PARAMS_C, ETA, IDELEC)
          dec1dd0 = dec1drs*drsdd0 + dec1dz*dzetadd0 + dec1ds*dsdd0
          dec1dd1 = dec1drs*drsdd1 + dec1dz*dzetadd1 + dec1ds*dsdd1
          dec1dg0 = dec1ds*dsdg0
          dec1dg1 = dec1ds*dsdg1
          dec1dg = dec1ds*dsdg ! got rid of unassigned variable dec1da, ADK

          !        Full correlation functional
          eps_c = (ec1 + ief*(ec0 - ec1))*dt
          dedd0 = ec1 + (ec0 - ec1)*ief &
          &    + dt*(ief*(dec0dd0 - dec1dd0) + dec1dd0 + (ec0 - ec1)*diefdd0)
          dedd1 = ec1 + (ec0 - ec1)*ief &
          &    + dt*(ief*(dec0dd1 - dec1dd1) + dec1dd1 + (ec0 - ec1)*diefdd1)
          dedg0 = dt* &
          &      (ief*(dec0dg0 - dec1dg0) + dec1dg0 + (ec0 - ec1)*diefdg0)
          dedg1 = dt* &
          &      (ief*(dec0dg1 - dec1dg1) + dec1dg1 + (ec0 - ec1)*diefdg1)
          ! VASP wants d E_c /d |grad n| rather than
          ! d E_c / d |grad n_sigma| for up and down
          ! spin-densities, ADK 26 August, 2020
          dedg = dt* &
          &      (ief*(dec0dg - dec1dg) + dec1dg + (ec0 - ec1)*diefdg)

          dedt0 = dt*(ec0 - ec1)*diefdt0
          dedt1 = dt*(ec0 - ec1)*diefdt1

      end subroutine scan_c_new

      subroutine scan_ec0(rs, s, zeta, ec0, dec0drs, dec0ds, dec0dz)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(scan_ec0)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: rs, s, zeta
          real(kdp), intent(out) :: ec0, dec0drs, dec0ds, dec0dz

          real(kdp) :: eclda, dldadrs, dldadrsrs
          real(kdp) :: dx_z, ddx_zdz, gc_z, dgc_zdz
          real(kdp) :: w0, dw0drs, ginf, dginfds
          real(kdp) :: h0, dh0drs, dh0ds

          real(kdp), parameter :: B1C = 0.0285764_kdp
          real(kdp), parameter :: CHI_LD = 0.12802585262625815_kdp
GPACC !$ACC ROUTINE(lda_0) SEQ

          call lda_0(rs, eclda, dldadrs, dldadrsrs, B1C)

          dx_z = ((1._kdp + zeta)**(4._kdp/3._kdp) + (1._kdp - zeta)**(4._kdp/3._kdp))/2._kdp
          ddx_zdz = -2._kdp*((1._kdp - zeta)**(1._kdp/3._kdp) - (1._kdp + zeta)**(1._kdp/3._kdp))/3._kdp

          gc_z = (1._kdp - 2.363_kdp*(dx_z - 1._kdp))*(1._kdp - zeta**12)
          dgc_zdz = -(1._kdp - 2.363_kdp*(dx_z - 1._kdp))*12._kdp*zeta**11
          dgc_zdz = dgc_zdz - 2.363_kdp*ddx_zdz*(1._kdp - zeta**12)

          w0 = exp(-eclda/B1C) - 1._kdp
          dw0drs = -(w0 + 1._kdp)*dldadrs/B1C

          ginf = 1._kdp/(1._kdp + 4._kdp*CHI_LD*s*s)**(1._kdp/4._kdp)
          dginfds = -2._kdp*CHI_LD*s/(1._kdp + 4._kdp*CHI_LD*s*s)**(5._kdp/4._kdp)

          h0 = B1C*log(1._kdp + w0*(1._kdp - ginf))
          dh0drs = B1C*(1._kdp - ginf)*dw0drs/(1._kdp + (1._kdp - ginf)*w0)
          dh0ds = -B1C*w0*dginfds/(1._kdp + (1._kdp - ginf)*w0)

          ec0 = (eclda + h0)*gc_z
          dec0drs = (dldadrs + dh0drs)*gc_z
          dec0dz = (h0 + eclda)*dgc_zdz
          dec0ds = dh0ds*gc_z

      end subroutine scan_ec0

      subroutine scan_ec1(rs, s, zeta, &
          &                ec1, dec1drs, dec1dz, dec1ds, IE_PARAMS_C, ETA, IDELEC)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(scan_ec1)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: rs, s, zeta
          real(kdp), intent(out) :: ec1, dec1drs, dec1dz, dec1ds
          real(kdp), dimension(7), intent(in) :: IE_PARAMS_C
          real(kdp), intent(in) :: ETA
          integer, intent(in) :: IDELEC

          real(kdp) :: sqrt_rs, dx_z, ddx_zdz, gc_z, dgc_zdz, phi, dphidz
          real(kdp) :: phi3, dphi3dz
          real(kdp) :: eclda0, declda0drs, declda0drsrs
          real(kdp) :: eclsda1, declsda1drs, declsda1dz, declsda1drsrs, declsda1drsz
          real(kdp) :: t, dtdrs, dtdz, dtds
          real(kdp) :: w1, dw1drs, dw1dz, y, dydrs, dydz, dyds
          real(kdp) :: del_y, ddel_ydrs, ddel_ydz, ddel_yds
          real(kdp) :: g_y, dg_ydrs, dg_ydz, dg_yds
          real(kdp) :: h1, dh1drs, dh1dz, dh1ds, D_DAMP2

          real(kdp), parameter :: B1C = 0.0285764_kdp
          real(kdp), parameter :: PI = 3.1415926535897932384626433832795_kdp
          real(kdp), parameter :: GAMMA = 0.0310906908696_kdp
          real(kdp), parameter :: BETA_MB = 0.066724550603149220_kdp
          real(kdp), parameter :: AFACTOR = 0.1_kdp
          real(kdp), parameter :: BFACTOR = 0.1778_kdp
          real(kdp) :: AFIX_T
GPACC !$ACC ROUTINE(lda_0) SEQ
GPACC !$ACC ROUTINE(lsda_1) SEQ
GPACC !$ACC ROUTINE(get_y) SEQ
GPACC !$ACC ROUTINE(get_del_y) SEQ

          AFIX_T = sqrt(PI/4._kdp)*(9._kdp*PI/4._kdp)**(1._kdp/6._kdp)

          dx_z = ((1._kdp + zeta)**(4._kdp/3._kdp) + (1._kdp - zeta)**(4._kdp/3._kdp))/2._kdp
          ddx_zdz = -2._kdp*((1._kdp - zeta)**(1._kdp/3._kdp) - (1._kdp + zeta)**(1._kdp/3._kdp))/3._kdp

          gc_z = (1._kdp - 2.3631_kdp*(dx_z - 1._kdp))*(1._kdp - zeta**12)
          dgc_zdz = -(1._kdp - 2.3631_kdp*(dx_z - 1._kdp))*12._kdp*zeta**11
          dgc_zdz = dgc_zdz - 2.3631_kdp*ddx_zdz*(1._kdp - zeta**12)

          phi = (exp((2._kdp/3._kdp)*log(1._kdp + zeta)) + exp((2._kdp/3._kdp)*log(1._kdp - zeta)))/2._kdp
          dphidz = (1._kdp/3._kdp)*((1._kdp + zeta)**(-1._kdp/3._kdp)-(1._kdp - zeta)**(-1._kdp/3._kdp))

          phi3 = phi**3
          dphi3dz = 3._kdp*phi**2*dphidz

          call lda_0(rs, eclda0, declda0drs, declda0drsrs, B1C)
          call lsda_1(rs, zeta, eclsda1, declsda1drs, declsda1dz, declsda1drsrs, declsda1drsz)

          sqrt_rs = sqrt(rs)

          t = AFIX_T*s/(sqrt_rs*phi)
          dtdrs = -AFIX_T*s/(2._kdp*phi*rs**(3.0/2.0))
          dtdz = -dphidz*AFIX_T*s/(sqrt_rs*phi**2)
          dtds = AFIX_T/(sqrt_rs*phi)

          w1 = exp(-eclsda1/(GAMMA*phi3)) - 1._kdp
          dw1drs = -(w1 + 1._kdp)*declsda1drs/(GAMMA*phi3)
          dw1dz = -(w1 + 1._kdp)/(GAMMA*phi3)*(declsda1dz - 3._kdp*eclsda1*dphidz/phi)

          call get_y(rs, t, dtdrs, dtdz, dtds, &
          &    w1, dw1drs, dw1dz, GAMMA, &
          &    y, dydrs, dydz, dyds)

          if (IDELEC .eq. 0) then
              del_y = 0._kdp
              ddel_ydrs = 0._kdp
              ddel_ydz = 0._kdp
              ddel_yds = 0._kdp

          else if (IDELEC .eq. 1 .or. IDELEC .eq. 2 .or. IDELEC .eq. 3) then
              if (IDELEC == 3) then
                D_DAMP2 = 5.259688_kdp
              else
                D_DAMP2 = 0.361_kdp
              end if

              call get_del_y(rs, s, zeta, &
              &  eclda0, declda0drs, declda0drsrs, &
              &  eclsda1, declsda1drs, declsda1dz, declsda1drsrs, declsda1drsz, &
              &  gc_z, dgc_zdz, &
              &  phi3, dphi3dz, w1, dw1drs, dw1dz, GAMMA, &
              &  del_y, ddel_ydrs, ddel_ydz, ddel_yds, IE_PARAMS_C, ETA, &
              &  D_DAMP2)

          else
__CPU         call quit('Bad IDELEC in SCAN')
          end if

          g_y = 1._kdp/(1._kdp + 4._kdp*(y - del_y))**(1.0/4.0)
          dg_ydrs = -(dydrs - ddel_ydrs) &
          &              /(1.d0 + 4.d0*(y - del_y))**(5.d0/4.d0)
          dg_ydz = -(dydz - ddel_ydz) &
          &              /(1.d0 + 4.d0*(y - del_y))**(5.d0/4.d0)
          dg_yds = -(dyds - ddel_yds) &
          &              /(1.d0 + 4.d0*(y - del_y))**(5.d0/4.d0)

          h1 = GAMMA*phi3*log(1._kdp + w1*(1._kdp - g_y))
          dh1drs = GAMMA*phi3*((1._kdp - g_y)*dw1drs - w1*dg_ydrs) &
          &          /(1._kdp + (1._kdp - g_y)*w1)
          dh1dz = GAMMA*log(1._kdp + (1._kdp - g_y)*w1)*dphi3dz &
          &        + GAMMA*phi3*((1._kdp - g_y)*dw1dz - w1*dg_ydz) &
          &          /(1._kdp + (1._kdp - g_y)*w1)
          dh1ds = -GAMMA*phi3*w1*dg_yds/(1._kdp + (1._kdp - g_y)*w1)

          ec1 = eclsda1 + h1
          dec1drs = declsda1drs + dh1drs
          dec1dz = declsda1dz + dh1dz
          dec1ds = dh1ds
      end subroutine scan_ec1

      subroutine get_y(rs, t, dtdrs, dtdz, dtds, &
          &    w1, dw1drs, dw1dz, GAMMA, &
          &    y, dydrs, dydz, dyds)
GPACC !$ACC ROUTINE SEQ

          implicit NONE
GPOMP  !$OMP DECLARE TARGET(get_y)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: rs, t, dtdrs, dtdz, dtds
          real(kdp), intent(in) :: w1, dw1drs, dw1dz, GAMMA
          real(kdp), intent(out) :: y, dydrs, dydz, dyds
          real(kdp) :: beta, dbetadrs

          real(kdp), parameter :: BETA_MB = 0.066725_kdp
          real(kdp), parameter :: AFACTOR = 0.1_kdp
          real(kdp), parameter :: BFACTOR = 0.1778_kdp

          beta = BETA_MB*(1._kdp + AFACTOR*rs)/(1._kdp + BFACTOR*rs)
          dbetadrs = BETA_MB*(AFACTOR - BFACTOR)/(1._kdp + BFACTOR*rs)**2

          y = beta/(GAMMA*w1)*t**2
          dydrs = t**2*dbetadrs/(GAMMA*w1) - beta*t**2*dw1drs/(GAMMA*w1**2) &
          &        + 2._kdp*beta*t*dtdrs/(GAMMA*w1)
          dydz = 2._kdp*beta*t*dtdz/(GAMMA*w1) - beta*t**2*dw1dz/(GAMMA*w1**2)
          dyds = 2._kdp*beta*t*dtds/(GAMMA*w1)
      end subroutine get_y

      subroutine get_del_y(rs, s, zeta, lda0, dlda0drs, dlda0drsrs, &
          &  lsda1, dlsda1drs, dlsda1dz, dlsda1drsrs, dlsda1drsz, &
          &  gc, dgcdz, &
          &  phi3, dphi3dz, w1, dw1drs, dw1dz, GAMMA, &
          &  del_y, ddel_ydrs, ddel_ydz, ddel_yds, IE_PARAMS_C, ETA, D_DAMP2)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(get_del_y)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: rs, s, zeta, lda0, dlda0drs, dlda0drsrs
          real(kdp), intent(in) :: lsda1, dlsda1drs, dlsda1dz, dlsda1drsrs, dlsda1drsz
          real(kdp), intent(in) :: gc, dgcdz
          real(kdp), intent(in) :: phi3, dphi3dz, w1, dw1drs, dw1dz, GAMMA
          real(kdp), intent(out) :: del_y, ddel_ydrs, ddel_ydz, ddel_yds
          real(kdp), dimension(7), intent(in) :: IE_PARAMS_C
          real(kdp), intent(in) :: ETA, D_DAMP2

          real(kdp) :: del_f2
          real(kdp) :: lsda0, dlsda0drs, dlsda0dz, dlsda0drsrs, dlsda0drsz
          real(kdp) :: K, dKdrs, dKdz
          real(kdp) :: t1, dt1drs, dt1dz, t2, dt2drs, dt2dz, t3, dt3drs, dt3dz
          real(kdp) :: damp, ddampds
          real(kdp) :: p, ds_z, dds_zdz
          integer :: i

          !real(kdp), parameter :: D_DAMP2 = 0.361_kdp

          p = s*s
          ds_z = ((1._kdp + zeta)**(5._kdp/3._kdp) + (1._kdp - zeta)**(5._kdp/3._kdp))/2._kdp
          dds_zdz = -5._kdp/6._kdp*((1._kdp - zeta)**(2._kdp/3._kdp) &
          &                         - (1._kdp + zeta)**(2._kdp/3._kdp))

          lsda0 = lda0*gc
          dlsda0drs = dlda0drs*gc
          dlsda0dz = lda0*dgcdz
          dlsda0drsrs = dlda0drsrs*gc
          dlsda0drsz = dlda0drs*dgcdz

          del_f2 = 0._kdp
          do i = 1, 7
              del_f2 = del_f2 + i*IE_PARAMS_C(i)
          end do

          t1 = del_f2/(27._kdp*GAMMA*ds_z*phi3*w1)
          dt1drs = -t1*dw1drs/w1
          dt1dz = -del_f2*(w1*(phi3*dds_zdz + ds_z*dphi3dz) + ds_z*phi3*dw1dz) &
          &         /(27._kdp*GAMMA*(ds_z*phi3*w1)**2)

          t2 = 20._kdp*rs*(dlsda0drs - dlsda1drs)
          dt2drs = 20._kdp*(dlsda0drs - dlsda1drs + rs*(dlsda0drsrs - dlsda1drsrs))
          dt2dz = 20._kdp*rs*(dlsda0drsz - dlsda1drsz)

          t3 = 45._kdp*ETA*(lsda0 - lsda1)
          dt3drs = 45._kdp*ETA*(dlsda0drs - dlsda1drs)
          dt3dz = 45._kdp*ETA*(dlsda0dz - dlsda1dz)

          K = t1*(t2 - t3)
          dKdrs = dt1drs*(t2 - t3) + t1*(dt2drs - dt3drs)
          dKdz = dt1dz*(t2 - t3) + t1*(dt2dz - dt3dz)

          damp = exp(-p**2/D_DAMP2**4)
          ddampds = -4._kdp*damp*s**3/D_DAMP2**4

          del_y = K*p*damp
          ddel_ydrs = p*damp*dKdrs
          ddel_ydz = p*damp*dKdz
          ddel_yds = K*s*(2._kdp*damp + s*ddampds)
      end subroutine get_del_y

      subroutine lda_0(rs, elda_0, dldadrs, dldadrsrs, B1C)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(lda_0)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: rs
          real(kdp), intent(out) :: elda_0, dldadrs, dldadrsrs
          real(kdp), intent(in) :: B1C

          real(kdp) :: sqrtrs

          real(kdp), parameter :: B2C = 0.0889_kdp
          real(kdp), parameter :: B3C = 0.125541_kdp

          sqrtrs = sqrt(rs)
          elda_0 = -B1C/(1._kdp + B2C*sqrtrs + B3C*rs)
          dldadrs = (B3C + B2C/(2._kdp*sqrtrs))*elda_0**2/B1C
          dldadrsrs = -B1C*(B2C + 3._kdp*B2C**2*sqrtrs + 9._kdp*B2C*B3C*rs + 8._kdp*B3C**2*rs*sqrtrs) &
          &           /(4._kdp*rs*sqrtrs*(1._kdp + B2C*sqrtrs + B3C*rs)**3)
      end subroutine lda_0

      subroutine lsda_1(rs, zeta, eclda1, declda1drs, declda1dz, declda1drsrs, declda1drsz)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(lsda_1)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: rs, zeta
          real(kdp), intent(out) :: eclda1, declda1drs, declda1dz
          real(kdp), intent(out) :: declda1drsrs, declda1drsz
          real(kdp) :: z3, z4
          real(kdp) :: eu, deudrs, deudrsrs
          real(kdp) :: ep, depdrs, depdrsrs
          real(kdp) :: alfm, dalfmdrs, dalfmdrsrs
          real(kdp) :: F, dFdz

          real(kdp), parameter :: GAM = 0.51984209978974632953442121455650_kdp
          real(kdp), parameter :: FZZ = 8._kdp/(9._kdp*GAM)
GPACC !$ACC ROUTINE(gcorn) SEQ

          call gcorn(0.03109070_kdp, 0.213700_kdp, 7.59570_kdp, 3.58760_kdp, &
          &          1.63820_kdp, 0.492940_kdp, rs, eu, deudrs, deudrsrs)
          call gcorn(0.015545350_kdp, 0.205480_kdp, 14.11890_kdp, 6.19770_kdp, &
          &          3.36620_kdp, 0.625170_kdp, rs, ep, depdrs, depdrsrs)
          call gcorn(0.01688690_kdp, 0.111250_kdp, 10.3570_kdp, 3.62310_kdp, &
          &          0.880260_kdp, 0.496710_kdp, rs, alfm, dalfmdrs, dalfmdrsrs)

          z3 = zeta**3
          z4 = zeta*z3

          F = ((1._kdp + zeta)**(4._kdp/3._kdp) + (1._kdp - zeta)**(4._kdp/3._kdp) - 2._kdp)/GAM
          dFdz = -4._kdp*((1._kdp - zeta)**(1._kdp/3._kdp) - (1._kdp + zeta)**(1._kdp/3._kdp))/(3._kdp*GAM)

          eclda1 = EU*(1._kdp - F*z4) + EP*F*z4 - ALFM*F*(1._kdp - z4)/FZZ
          declda1drs = (1._kdp - z4*F)*deudrs + z4*F*depdrs &
          &              - (1._kdp - z4)*F*dalfmdrs/FZZ
          declda1dz = EU*(-4._kdp*z3*F - z4*dFdz) + z4*EP*dFdz + 4._kdp*z3*EP*F &
          &              + 4._kdp*z3*alfm*F/FZZ - (1._kdp - z4)*alfm*dFdz/FZZ

          ! Some second derivatives are required for r2scan/r4scan correlation.
          declda1drsrs = z4*F*depdrsrs + (1._kdp - z4*F)*deudrsrs &
          &               - (1._kdp - z4)*F*dalfmdrsrs/FZZ
          declda1drsz = (4._kdp*z3*F*(dalfmdrs + FZZ*(depdrs - deudrs)) &
          &             + ((z4 - 1._kdp)*dalfmdrs + FZZ*z4*(depdrs - deudrs))*dFdz)/FZZ
      end subroutine lsda_1

      subroutine gcorn(A, A1, B1, B2, B3, B4, rs, GG, GGRS, GGRSRS)
GPACC !$ACC ROUTINE SEQ

          IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(gcorn)
          integer, parameter :: kdp = selected_real_kind(15, 307)

          real(kdp), intent(in) :: A, A1, B1, B2, B3, B4, rs
          real(kdp), intent(out) :: GG, GGRS, GGRSRS
          real(kdp) :: rtrs, Q0, Q1, Q2
          real(kdp) :: Q0RS, Q1RS, Q1RSRS, Q2RS, Q2RSRS

          rtrs = sqrt(rs)

          Q0 = -2._kdp*A*(1._kdp + A1*rs)
          Q0RS = -2._kdp*A*A1

          Q1 = 2._kdp*A*rtrs*(B1 + rtrs*(B2 + rtrs*(B3 + B4*rtrs)))
          Q1RS = A*(2._kdp*B2 + B1/rtrs + 3._kdp*B3*rtrs + 4._kdp*B4*rs)
          Q1RSRS = A*(4._kdp*B4 - B1/(2._kdp*rs*rtrs) + 3._kdp*B3/(2._kdp*rtrs))

          Q2 = log(1._kdp + 1._kdp/Q1)
          Q2RS = -Q1RS/(Q1**2 + Q1)
          Q2RSRS = ((1._kdp + 2._kdp*Q1)*Q1RS**2 - Q1*(1._kdp + Q1)*Q1RSRS)/(Q1**2*(1._kdp + Q1)**2)

          GG = Q0*Q2
          GGRS = Q0*Q2RS + Q2*Q0RS
          GGRSRS = 2._kdp*Q0RS*Q2RS + Q0*Q2RSRS
      end subroutine gcorn

      subroutine quit(astr)
        implicit none
        character,intent(in) :: astr
        write(*,*) astr
        stop
      end subroutine quit

      subroutine scan_orb_free_xc(nu,nd,gnu,gnd,gn,lnu,lnd,&
       & exc,vxcu,vxcd,vxcgu,vxcgd,vxcg,vxclu,vxcld,id_xc,AMGGAX,AMGGAC)
GPACC !$ACC ROUTINE SEQ

        USE xc_name
        implicit none
GPOMP  !$OMP DECLARE TARGET(scan_orb_free_xc)

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: id_xc
        real(dp), dimension(2) :: cppars
        real(dp), intent(in) :: nu,nd,gnu,gnd,gn,lnu,lnd
        real(dp), intent(in) :: AMGGAX,AMGGAC
        real(dp), intent(out) :: exc,vxcu,vxcd,vxcgu,vxcgd,vxcg,vxclu,vxcld

        real(dp) :: tauu,taud,dtdnu,dtdgnu,dtdlnu,dtdnd,dtdgnd,dtdlnd
        real(dp) :: ex,vxu,vxd,vxgu,vxgd,vxtu,vxtd
        real(dp) :: ec,vcu,vcd,vcgu,vcgd,vctu,vctd,vxct
GPACC !$ACC ROUTINE(tau_rpp_unp) SEQ
GPACC !$ACC ROUTINE(pc_tau_unp) SEQ
GPACC !$ACC ROUTINE(VSCAN2x) SEQ
GPACC !$ACC ROUTINE(VSCAN2c) SEQ

        cppars = (/ 1.784720_dp, 0.258304_dp /)

        exc = 0._dp
        vxcu = 0._dp ; vxcd = 0._dp
        vxcgu = 0._dp ; vxcgd = 0._dp ; vxcg = 0._dp
        vxclu = 0._dp ; vxcld = 0._dp

        if (id_xc==ID_XC_OFR2) then
          call tau_rpp_unp(2*nu,2*gnu,2*lnu,tauu,dtdnu,dtdgnu,dtdlnu)
          call tau_rpp_unp(2*nd,2*gnd,2*lnd,taud,dtdnd,dtdgnd,dtdlnd)
        else
          call pc_tau_unp(cppars,2*nu,2*gnu,2*lnu,tauu,dtdnu,dtdgnu,dtdlnu)
          call pc_tau_unp(cppars,2*nd,2*gnd,2*lnd,taud,dtdnd,dtdgnd,dtdlnd)
        endif
        tauu = tauu*0.5_dp
        taud = taud*0.5_dp

        call VSCAN2x(nu,nd,gnu,gnd,gn,tauu,taud, &
        &   ex,vxu,vxgu,vxd,vxgd,vxtu,vxtd,id_xc)
        call VSCAN2c(nu,nd,gnu,gnd,gn,tauu,taud,&
            &   ec,vcu,vcd,vcgu,vcgd,vxcg,vctu,vctd,id_xc)

        exc = (AMGGAX*ex + AMGGAC*ec)/(nu + nd)

        vxct = AMGGAX*vxtu + AMGGAC*vctu
        vxcu = AMGGAX*vxu + AMGGAC*vcu + vxct*dtdnu
        vxcgu = AMGGAX*vxgu + vxct*dtdgnu
        vxclu = vxct*dtdlnu

        vxct = AMGGAX*vxtd + AMGGAC*vctd
        vxcd = AMGGAX*vxd + AMGGAC*vcd + vxct*dtdnd
        vxcgd = AMGGAX*vxgd + vxct*dtdgnd
        vxcld = vxct*dtdlnd

      end subroutine scan_orb_free_xc

      subroutine tau_rpp_unp(n,gn,ln,tau,dtdn,dtdgn,dtdln)
GPACC !$ACC ROUTINE SEQ

        implicit none
GPOMP  !$OMP DECLARE TARGET(tau_rpp_unp)
        integer, parameter :: dp = selected_real_kind(15, 307)
        real(dp), parameter :: pi = 3.14159265358979323846264338327950288419_dp

        real(dp), intent(in) :: n,gn,ln
        real(dp), intent(out) :: tau,dtdn,dtdgn,dtdln

        !real(dp), parameter :: x0 = 0.5245_dp, ca = 6._dp/x0**2, cb =  -8._dp/x0**3
        !real(dp), parameter :: cc = 3._dp/x0**4
        !real(dp), parameter :: x0 = 0.529772_dp!0.5583_dp

        ! tau-specific parameters
        real(dp), parameter :: x0 = 0.819411_dp, c1= 0.201352_dp, c2 = 0.185020_dp
        real(dp), parameter :: c3 = 1.53804_dp
        real(dp), parameter :: ca = 20._dp/x0**3, cb =  -45._dp/x0**4
        real(dp), parameter :: cc = 36._dp/x0**5, cd = -10._dp/x0**6

        ! r2 SCAN-specific parameters
        real(dp), parameter :: eta = 1.d-3, k0 = 0.174_dp, k1 = 0.065_dp
        real(dp), parameter :: uak = 10._dp/81._dp, epne = (8._dp + 9._dp*eta)
      !  real(dp), parameter :: cxv(8) = (/ 1._dp, -0.667_dp, -0.4445555_dp,&
      !  &     -0.663086601049_dp, 1.451297044490_dp, -0.887998041597_dp,&
      !  &     0.234528941479_dp, -0.023185843322_dp /)
        real(dp), parameter :: f1 = -0.93530008755199956_dp, f2 = 0.85003592049200183_dp
        real(dp), parameter :: h1 = 5._dp*(4._dp + 9._dp*eta)*k0*f1/27._dp + uak
        real(dp), parameter :: h2 = -2._dp*h1**2/k1

        integer :: i
        real(dp) :: kf,tau0,t0,p,q,ft!,d_p_dn,d_p_dgn,d_q_dn,d_q_dlap
        real(dp) :: d_tau0_dn, d_ft_dp,d_ft_dq, x, alp, dalpdx,dxdp,dxdq
        real(dp) :: cpp,cpq,cqq,f4,f4damp,d_f4_dp,d_f4_dq, d_f4damp_dp,d_f4damp_dq
        real(dp) :: f4c, d_f4c_dp!, f1,f2,h1,h2

        ! first the approximate tau

        kf = (3*pi**2*n)**(1._dp/3._dp)
        t0 = 0.3_dp*kf**2
        tau0 = t0*n

        p = (gn/(2*kf*n))**2
        q = ln/(4*kf**2*n)

        ! fourth-order gradient correction for r2 SCAN

        !f1 = 0._dp
        !f2 = 0._dp
        !do i = 1,7
        !    f1 = f1 + (1._dp*i)*cxv(i+1)
        !    f2 = f2 + i*(i-1._dp)*cxv(i+1)
        !end do
        !h1 = 5._dp*(4._dp + 9._dp*eta)*k0*f1/27._dp + uak
        !h2 = -2._dp*h1**2/k1

        cqq = (146._dp/2025._dp/k0 - 200._dp/81._dp*f2)/f1
        cpq = 100._dp/243._dp*epne*f2/f1 + 20._dp/9._dp*h1/k0 &
        &     + 100._dp*eta/27._dp - 73._dp/405._dp/(k0*f1)
        cpp = -(h2/2._dp + 5._dp/27*epne*h1*f1 + 25._dp/1458._dp*epne**2*k0*f2 )/(k0*f1)&
        &     - 25._dp*eta**2/9._dp - 200._dp*eta/81._dp

        f4 = (cpp-c3)*p**2 + cpq*p*q + cqq*q**2
        d_f4_dp = 2._dp*(cpp-c3)*p + cpq*q
        d_f4_dq = cpq*p + 2._dp*cqq*q
        f4damp = exp(-(p/c1)**2 - (q/c2)**2)
        d_f4damp_dp = -2._dp*p/c1**2*f4damp
        d_f4damp_dq = -2._dp*q/c2**2*f4damp

        f4c = c3*p**2*exp(-c3*p)
        d_f4c_dp = c3*p*(2._dp - c3*p)*exp(-c3*p)

        x = 1._dp - 40._dp/27._dp*p + 20._dp/9._dp*q + f4*f4damp + f4c
        dxdp = -40._dp/27._dp + d_f4_dp*f4damp + f4*d_f4damp_dp + d_f4c_dp
        dxdq = 20._dp/9._dp + d_f4_dq*f4damp + f4*d_f4damp_dq

        if (x < 0._dp) then
          alp = 0._dp
          dalpdx = 0._dp
        else if ((0._dp <= x).and.(x <= x0)) then
          !alp = x**3*(ca + x*(cb + x*cc))
          !dalpdx = x**2*(3._dp*ca + x*(4._dp*cb + 5._dp*x*cc))
          alp = x**4*(ca + x*(cb + x*(cc + x*cd)))
          dalpdx = x**3*(4._dp*ca + x*(5._dp*cb + x*(6._dp*cc + x*7._dp*cd)))
        else if (x > x0) then
          alp = x
          dalpdx = 1._dp
        end if

        ft = 5._dp*p/3._dp + alp
        d_ft_dp = 5._dp/3._dp + dalpdx*dxdp
        d_ft_dq = dalpdx*dxdq

        tau = tau0*ft

        d_tau0_dn = kf**2/2._dp

        !d_p_dn = -8._dp/3._dp *p/n
        !d_p_dgn = gn/(2*(kf*n)**2)

        !d_q_dn = -5._dp/3._dp*q/n
        !d_q_dlap = 1._dp/(4*kf**2*n)

        !dtdn = ft*d_tau0_dn + tau0*(d_ft_dp*d_p_dn + d_ft_dq*d_q_dn)
        dtdn = ft*d_tau0_dn - t0/3._dp*(8._dp*p*d_ft_dp + 5._dp*q*d_ft_dq)
        ! NB: tau0*d_p_dgn = 3/20 * |grad n|/n = 0.15 * |grad n|/n
        dtdgn = 0.15_dp*gn/n*d_ft_dp!tau0*d_ft_dp*d_p_dgn
        ! NB: tau0*d_q_dln = 3/40 = 0.075
        dtdln = 0.075_dp*d_ft_dq!tau0*d_ft_dq*d_q_dlap

      end subroutine tau_rpp_unp

      subroutine pc_tau_unp(pars,n,gn,ln,tau,dtdn,dtdgn,dtdln)
GPACC !$ACC ROUTINE SEQ

        implicit none
GPOMP  !$OMP DECLARE TARGET(pc_tau_unp)
        integer, parameter :: dp = selected_real_kind(15, 307)
        real(dp), parameter :: pi = 3.14159265358979323846264338327950288419_dp

        real(dp), intent(in),dimension(2) :: pars
        real(dp), intent(in) :: n,gn,ln
        real(dp), intent(out) :: tau,dtdn,dtdgn,dtdln

        real(dp) :: kf,tau0,t0,p,q,ftw,ft2,ft4,ft
        real(dp) :: freg,ftge4,z,theta,d_theta_dz
        real(dp) :: d_tau0_dn!,d_p_dn,d_p_dgn,d_q_dn,d_q_dln
        real(dp) :: d_ftw_dp, d_ft2_dp, d_ft2_dq, d_ft4_dp, d_ft4_dq
        real(dp) :: d_freg_dp, d_freg_dq,d_ftge4_dp,d_ftge4_dq
        real(dp) :: d_z_dp,d_z_dq, d_ft_dp,d_ft_dq
GPACC !$ACC ROUTINE(pc_interpolation) SEQ

        ! first the approximate tau

        kf = (3*pi**2*n)**(1._dp/3._dp)
        t0 = 0.3_dp*kf**2
        tau0 = t0*n

        p = (gn/(2*kf*n))**2
        q = ln/(4*kf**2*n)

        ftw = 5*p/3._dp
        ft2 = 5*p/27._dp + 20*q/9._dp
        ft4 = 8*q**2/81._dp - p*q/9._dp + 8*p**2/243._dp

        freg = (1._dp + (ft4/(1._dp + ftw))**2)**(0.5_dp)
        ftge4 = (1._dp + ft2 + ft4)/freg
        z = ftge4 - ftw
        call pc_interpolation(pars,z,theta,d_theta_dz)

        ft = ftw + z*theta
        tau = ft*tau0

        ! now the derivatives of tau

        d_tau0_dn = kf**2/2._dp

        !d_p_dn = -8._dp/3._dp *p/n
        !d_p_dgn = gn/(2*(kf*n)**2)

        !d_q_dn = -5._dp/3._dp *q/n
        !d_q_dln = 1._dp/(4*kf**2*n)

        d_ftw_dp = 5._dp/3._dp
        d_ft2_dp = 5._dp/27._dp
        d_ft2_dq = 20._dp/9._dp

        d_ft4_dp = -q/9._dp + 16*p/243._dp
        d_ft4_dq = 16*q/81._dp - p/9._dp

        d_freg_dp = ft4*(d_ft4_dp*(1._dp + ftw) - ft4*d_ftw_dp)/(freg*(1._dp + ftw)**3)
        d_freg_dq = ft4*d_ft4_dq/(freg*(1._dp + ftw)**2)

        d_ftge4_dp = (d_ft2_dp + d_ft4_dp - ftge4*d_freg_dp)/freg
        d_ftge4_dq = (d_ft2_dq + d_ft4_dq - ftge4*d_freg_dq)/freg

        d_z_dp = d_ftge4_dp - d_ftw_dp
        d_z_dq = d_ftge4_dq

        d_ft_dp = d_ftw_dp + (d_theta_dz*z + theta)*d_z_dp
        d_ft_dq = (d_theta_dz*z + theta)*d_z_dq

        !dtdn = ft*d_tau0_dn + tau0*( d_ft_dp*d_p_dn + d_ft_dq*d_q_dn )
        dtdn = ft*d_tau0_dn - t0/3._dp*(8._dp*p*d_ft_dp + 5._dp*q*d_ft_dq)
        ! NB: tau0*d_p_dgn = 3/20 * |grad n|/n = 0.15 * |grad n|/n
        dtdgn = 0.15_dp*gn/n*d_ft_dp!tau0*d_ft_dp*d_p_dgn
        ! NB: tau0*d_q_dln = 3/40 = 0.075
        dtdln = 0.075_dp*d_ft_dq!tau0*d_ft_dq*d_q_dln

      end subroutine pc_tau_unp

      subroutine pc_interpolation(pars,z,th,dthdz)
GPACC !$ACC ROUTINE SEQ

        implicit none
GPOMP !$OMP DECLARE TARGET(pc_interpolation)

        integer, parameter :: dp = selected_real_kind(15, 307)
        real(dp), parameter :: eps = 1.d-2

        real(dp), intent(in),dimension(2) :: pars
        real(dp), intent(in) :: z
        real(dp), intent(out) :: th,dthdz

        real(dp) :: a_fuzzy,fac1,fac2,dfac1,dfac2

        a_fuzzy = pars(1)-eps

        th = 0._dp
        dthdz = 0._dp
        if ((eps < z).and.(z < a_fuzzy)) then
          fac1 = exp(pars(1)/(pars(1) - z))
          dfac1 = pars(1)/(pars(1) - z)**2*fac1
          fac2 = exp(pars(1)/z)
          dfac2 = -pars(1)/z**2*fac2

          th = ( (1._dp + fac1)/(fac2 + fac1) )**pars(2)
          dthdz = pars(2)* ( (1._dp + fac1)/(fac2 + fac1) )**(pars(2) - 1)
          dthdz = dthdz*( dfac1*(fac2 - 1._dp) - dfac2*(1._dp + fac1) )/(fac2 + fac1)**2
        else if (z >= a_fuzzy) then
          th = 1._dp
        end if

      end subroutine pc_interpolation

      SUBROUTINE V1SRTMx(RX,DRX,TAUX,Ex,VXD,VXDD,AMUXD)
GPACC !$ACC ROUTINE SEQ
!     Simplified regularized Tao-Mo variant 1 (v1-sregTM)
!     Functional (Exchange part only)
!
!     Written by:
!     H. Francisco
!     QTP, Department of Physics, University of Florida
!
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214102(2023)               
!
               
      USE prec
      USE constant

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(V1SRTMx)

      REAL(q) :: RX,DRX,TAUX,EX,VXD,VXDD,AMUXD, CMU
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: PISQ,AX
      REAL(q) :: TAU_RHO,TAUW_RHO,TAU0,TAUW_TEMP
      REAL(q) :: RHO,DRHO,P,P2
      REAL(q) :: FA,DFADA,ALPHA,DALPHAD,DALPHADD,DALPHADTAU
      REAL(q) :: DPD,DPDD,EXUNIF,EXLDA,EXDLDA
      REAL(q) :: FX,DFXDP,DFXDA,DFX_D,DFX_DD,DFX_DTAU
! start Tao-Mo parameters and VARIABLES
      REAL(q) :: w, fdme, fxsc, zp, ftm, fd, rr, y, z2, z3, qg, z1,ftt
      REAL(q) :: lamda,betha, delta, con, F32, pma
      REAL(q) :: ftm1, ftm2,ftm3,F110,F35,F910,fxsc1, rar
      REAL(q) :: dwdr,dwdzp,dzpdp,dzpdftm,dftmdp,dzpda,dftmda
      REAL(q) :: dfdmedr,dfdmedfd,dfdmedrr,drrdr,dfddy,dydp
      REAL(q) :: drrdz2, dz2dp, drrdz3, dz3dr,dz3dz2,dz3da
      REAL(q) :: dfxscdr,dfxscdp,dfxscdqg,dqgdr,dfxscdz1,dz1dr
      REAL(q) :: dqgdp,dqgda,dz1dp,dz1da
      REAL(q) :: dwdg,dfdmedg,dfxscdg,dqgdg,dz1dg
      REAL(q) :: dwdl,dfdmedl,dfxscdl,lapq,laf
      REAL(q) :: dwdt,dfdmedt,dfxscdt,fds110
      REAL(q) :: drrdp, drrda,gam12
      REAL(q) :: F13, F23, F43, F53, F83, F18
      REAL(q) :: F727, F554, F12, NINE, NITWE,TEN
      REAL(q) :: F729,F2025, F405, TWENTY, EIGHT, FIF,T6
      REAL(q) :: zrp,ftm4,ftm5, F52
      REAL(q) :: z1a1, z1a2, z1a3, z1a4, z1a5, z1p1
      REAL(q) :: epsp, pmarev
!  
      PARAMETER (F43=4._q/3._q)
      PARAMETER (F13=1._q/3._q)
      PARAMETER (F83=8._q/3._q)
      PARAMETER (F23=2._q/3._q)
      PARAMETER (F53=5._q/3._q)
      PARAMETER (betha=79.873_q)
      PARAMETER (lamda=0.6866_q)
      PARAMETER (F110=1._q/10._q)
      PARAMETER (delta=1.475_q)
      PARAMETER (con=3._q)
      PARAMETER (F32=3._q/2._q)
      PARAMETER (F35=1._q/F53)
      PARAMETER (F910=9._q/10._q)
      PARAMETER (F727=700._q/27._q)
      PARAMETER (F554=595._q/54._q)
      PARAMETER (F729=50._q/729._q)
      PARAMETER (F2025=146._q/2025._q)
      PARAMETER (F405=73._q/405._q)
      PARAMETER (FIF=15._q)
      PARAMETER (T6=0.0000000001_q)
! finish Tao-MO parameters and varibles
      PARAMETER (CMU=10._q/81._q)
! other parameters
      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (FIVE=5._q)
      PARAMETER (SIX=6._q)
      PARAMETER (SEVEN=7._q)
      PARAMETER (EIGHT=8._q)
      PARAMETER (NINE= 9._q)
      PARAMETER (TEN=10._q)
      PARAMETER (TWENTY=20._q)
      PARAMETER (F18 = ONE/EIGHT)
      PARAMETER (F12 = ONE/TWO)
      PARAMETER (NITWE=NINE/TWENTY)
      PARAMETER (F52=FIVE/TWO)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
      PARAMETER (epsp=0.5_q)

      VXD=0._q
      VXDD=0._q
      AMUXD=0._q

!      rC1=(TWENTY*rK0/(TSEVEN*uge))**FOUR - (ONE + rE1)*(ONE + rE1)
! Suspect that TAUWU and TAUWD are not well described. Use
! TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUW_TEMP=0.125_q*DRX**2._q/RX
! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RX
      DRHO=TWO*DRX
      TAUW_RHO=TWO*TAUW_TEMP
      TAU_RHO=TWO*TAUX
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4
! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      DPD=-THRD8*P/RHO
      DPDD=TWO*DRHO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)

      ALPHA = (TAU_RHO - TAUW_RHO)/TAU0
      if (ALPHA .LT. ZERO)  ALPHA=ZERO
      DALPHAD = THRD*(EIGHT*TAUW_RHO - FIVE*TAU_RHO)/(RHO*TAU0)
      DALPHADD = -0.25_q*DRHO/(RHO*TAU0)
      DALPHADTAU = ONE/TAU0

      P2 = P*P
      pmarev=FIVE*P + THREE*ALPHA + epsp
      zrp= (FIVE*P + epsp)/pmarev

      IF (zrp .GT. ONE) zrp=ONE
!
      pma = FIVE*P + THREE*ALPHA
!
            if (ALPHA .lt. T6) then
                zp=ONE
            else if (P .lt. T6) then
                zp=ZERO
            else
                zp=FIVE*P/pma
              endif
      IF (zp .GT. ONE) zp=ONE

      w   = (zrp**2 + THREE*zrp**3)/((ONE + zrp**3)**2)
!
       z2  = F53*P
       z3  = z2 + ALPHA
       y   = P*((TWO*lamda - ONE)**2)
       fds110=ONE + F727*y + betha*y*y
       fd  = fds110**F110
       rr = ONE + F554*y - (z3 - THREE*((lamda**2) - &
        &   lamda + F12)*(z3 - ONE - (z2/NINE)))
       fdme=ONE/(fd**2) + SEVEN*rr/(NINE*(fd**4))
!
       qg  = NITWE*(ALPHA-ONE) + F23*P
      fxsc1= ONE + TEN*((CMU + F729*P)*P + F2025*qg*qg &
        &      - F405*qg*F35*zp*(ONE-zp))
      fxsc=fxsc1**F110
!
      FX = w*fdme + (ONE-w)*fxsc
!
!     Building derivatives of Tao-Mo-FX
      dwdzp= ((TWO*zrp + NINE*zrp*zrp)/((ONE + zrp**3)**2)) - &
       & ((SIX*zrp*zrp*(zrp*zrp + THREE*zrp*zrp*zrp))/((ONE + zrp**3)**3))
!
            dz1da=-THREE*(epsp + FIVE*P)/(pmarev**2)
            dz1dp=15._q*ALPHA/(pmarev**2)
!
            if (ALPHA .lt. T6) then
                dzpda=ZERO
              else
                dzpda=-FIF*P/(pma**2)
              endif
!
            if (p .lt. T6) then
                dzpdp=ZERO
              else
                dzpdp=FIF*ALPHA/(pma**2)
             endif
!
            dwdr = dwdzp*(dz1dp*DPD + dz1da*DALPHAD)
            dwdg = dwdzp*(dz1dp*DPDD + dz1da*DALPHADD)
            dwdt = dwdzp*dz1da*DALPHADTAU
!
            dydp = (TWO*lamda - ONE)**2
            dfddy=(F727 + TWO*betha*y)/(TEN*(fds110**F910))
            dfdmedfd=(-TWO/(fd**3)) - ((28._q*rr)/(NINE*(fd**5)))
            dfdmedrr=SEVEN/(NINE*(fd**4))
            drrdp   = (40._q/NINE)*(F12 - lamda + lamda*lamda) - F53 &
                   &  + (F554*((TWO*lamda - ONE)**2))
            drrda   = -ONE + (THREE*(F12 - lamda + lamda*lamda))
!
            dfdmedr=dfdmedfd*dfddy*dydp*DPD + &
                  & dfdmedrr*(drrdp*DPD + drrda*DALPHAD)
            dfdmedg=dfdmedfd*dfddy*dydp*DPDD + &
                  & dfdmedrr*(drrdp*DPDD + drrda*DALPHADD)
            dfdmedt=dfdmedrr*drrda*DALPHADTAU
!
            dqgdp   =F23
            dqgda   =NITWE
            dfxscdp =(CMU + (100._q/729._q)*P)/(fxsc1**F910)
            dfxscdz1 = (((-73._q/675._q)*qg*(ONE-zp)) + &
                     &  ((73._q/675_q)*qg*zp))/(fxsc1**F910)
            dfxscdqg=(((292._q/2025._q)*qg) &
                     &  -((73._q/675._q)*zp*(ONE-zp)))/(fxsc1**F910)
!            
         dfxscdr=dfxscdp*DPD + dfxscdz1*(dzpdp*DPD + dzpda*DALPHAD) &
                &   +  dfxscdqg*(dqgdp*DPD + dqgda*DALPHAD)
         dfxscdg=dfxscdp*DPDD + dfxscdz1*(dzpdp*DPDD + dzpda*DALPHADD) &
                &   +  dfxscdqg*(dqgdp*DPDD + dqgda*DALPHADD)
         dfxscdt=DALPHADTAU*(dfxscdz1*dzpda + dfxscdqg*dqgda)
!
! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D   =dwdr*fdme + dfdmedr*w - dwdr*fxsc + dfxscdr*(ONE-w)
      DFX_DD  =dwdg*fdme + dfdmedg*w - dwdg*fxsc + dfxscdg*(ONE-w)
      DFX_DTAU=dwdt*fdme + dfdmedt*w - dwdt*fxsc + dfxscdt*(ONE-w)

!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD=EXDLDA*FX+EXLDA*DFX_D
      VXDD=EXLDA*DFX_DD
      AMUXD=EXLDA*DFX_DTAU
      EX=EX+EXLDA*FX/TWO

      RETURN
      END SUBROUTINE V1SRTMx

      SUBROUTINE V2SRTMx(RX,DRX,TAUX,Ex,VXD,VXDD,AMUXD)
GPACC !$ACC ROUTINE SEQ
!     Simplified regularized Tao-Mo variant 2 (v2-sregTM)
!     Functional (Exchange part only)
!
!     Written by:
!     H. Francisco
!     QTP, Department of Physics, University of Florida
!
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214102(2023)               
!
              
      USE prec
      USE constant

      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(V2SRTMx)

      REAL(q) :: RX,DRX,TAUX,EX,VXD,VXDD,AMUXD, CMU
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: PISQ,AX
      REAL(q) :: TAU_RHO,TAUW_RHO,TAU0,TAUW_TEMP
      REAL(q) :: RHO,DRHO,P,P2
      REAL(q) :: FA,DFADA,ALPHA,DALPHAD,DALPHADD,DALPHADTAU
      REAL(q) :: DPD,DPDD,EXUNIF,EXLDA,EXDLDA
      REAL(q) :: FX,DFXDP,DFXDA,DFX_D,DFX_DD,DFX_DTAU
! start Tao-Mo parameters and VARIABLES
      REAL(q) :: w, fdme, fxsc, zp, ftm, fd, rr, y, z2, z3, qg, z1,ftt
      REAL(q) :: lamda,betha, delta, con, F32, pma
      REAL(q) :: ftm1, ftm2,ftm3,F110,F35,F910,fxsc1, rar
      REAL(q) :: dwdr,dwdzp,dzpdp,dzpdftm,dftmdp,dzpda,dftmda
      REAL(q) :: dfdmedr,dfdmedfd,dfdmedrr,drrdr,dfddy,dydp
      REAL(q) :: drrdz2, dz2dp, drrdz3, dz3dr,dz3dz2,dz3da
      REAL(q) :: dfxscdr,dfxscdp,dfxscdqg,dqgdr,dfxscdz1,dz1dr
      REAL(q) :: dqgdp,dqgda,dz1dp,dz1da
      REAL(q) :: dwdg,dfdmedg,dfxscdg,dqgdg,dz1dg
      REAL(q) :: dwdl,dfdmedl,dfxscdl,lapq,laf
      REAL(q) :: dwdt,dfdmedt,dfxscdt,fds110
      REAL(q) :: drrdp, drrda,gam12
      REAL(q) :: F13, F23, F43, F53, F83, F18
      REAL(q) :: F727, F554, F12, NINE, NITWE,TEN
      REAL(q) :: F729,F2025, F405, TWENTY, EIGHT, FIF,T6
      REAL(q) :: zrp,ftm4,ftm5, F52
      REAL(q) :: z1a1, z1a2, z1a3, z1a4, z1a5, z1p1
      REAL(q) :: epsp, pmarev
!  
      PARAMETER (F43=4._q/3._q)
      PARAMETER (F13=1._q/3._q)
      PARAMETER (F83=8._q/3._q)
      PARAMETER (F23=2._q/3._q)
      PARAMETER (F53=5._q/3._q)
      PARAMETER (betha=79.873_q)
      PARAMETER (lamda=0.6866_q)
      PARAMETER (F110=1._q/10._q)
      PARAMETER (delta=1.475_q)
      PARAMETER (con=3._q)
      PARAMETER (F32=3._q/2._q)
      PARAMETER (F35=1._q/F53)
      PARAMETER (F910=9._q/10._q)
      PARAMETER (F727=700._q/27._q)
      PARAMETER (F554=595._q/54._q)
      PARAMETER (F729=50._q/729._q)
      PARAMETER (F2025=146._q/2025._q)
      PARAMETER (F405=73._q/405._q)
      PARAMETER (FIF=15._q)
      PARAMETER (T6=0.0000000001_q)
! finish Tao-MO parameters and varibles
      PARAMETER (CMU=10._q/81._q)
! other parameters
      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (FIVE=5._q)
      PARAMETER (SIX=6._q)
      PARAMETER (SEVEN=7._q)
      PARAMETER (EIGHT=8._q)
      PARAMETER (NINE= 9._q)
      PARAMETER (TEN=10._q)
      PARAMETER (TWENTY=20._q)
      PARAMETER (F18 = ONE/EIGHT)
      PARAMETER (F12 = ONE/TWO)
      PARAMETER (NITWE=NINE/TWENTY)
      PARAMETER (F52=FIVE/TWO)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
      PARAMETER (epsp=0.5_q)

      VXD=0._q
      VXDD=0._q
      AMUXD=0._q

!      rC1=(TWENTY*rK0/(TSEVEN*uge))**FOUR - (ONE + rE1)*(ONE + rE1)

! Suspect that TAUWU and TAUWD are not well described. Use
! TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUW_TEMP=0.125_q*DRX**2._q/RX

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RX
      DRHO=TWO*DRX
      TAUW_RHO=TWO*TAUW_TEMP
      TAU_RHO=TWO*TAUX
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4
! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      DPD=-THRD8*P/RHO
      DPDD=TWO*DRHO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)

      ALPHA = (TAU_RHO - TAUW_RHO)/TAU0
      if (ALPHA .LT. ZERO)  ALPHA=ZERO
      DALPHAD = THRD*(EIGHT*TAUW_RHO - FIVE*TAU_RHO)/(RHO*TAU0)
      DALPHADD = -0.25_q*DRHO/(RHO*TAU0)
      DALPHADTAU = ONE/TAU0

      P2 = P*P
      pmarev=FIVE*P + THREE*ALPHA + epsp
      zrp= (FIVE*P + epsp)/pmarev

      IF (zrp .GT. ONE) zrp=ONE
      zp=zrp
!
      w   = (zrp**2 + THREE*zrp**3)/((ONE + zrp**3)**2)
!
       z2  = F53*P
       z3  = z2 + ALPHA
       y   = P*((TWO*lamda - ONE)**2)
       fds110=ONE + F727*y + betha*y*y
       fd  = fds110**F110
       rr = ONE + F554*y - (z3 - THREE*((lamda**2) - &
        &   lamda + F12)*(z3 - ONE - (z2/NINE)))
       fdme=ONE/(fd**2) + SEVEN*rr/(NINE*(fd**4))
!
       qg  = NITWE*(ALPHA-ONE) + F23*P
      fxsc1= ONE + TEN*((CMU + F729*P)*P + F2025*qg*qg &
        &      - F405*qg*F35*zp*(ONE-zp))
      fxsc=fxsc1**F110
!
      FX = w*fdme + (ONE-w)*fxsc
!     Building derivatives of Tao-Mo-FX
      dwdzp= ((TWO*zrp + NINE*zrp*zrp)/((ONE + zrp**3)**2)) - &
       & ((SIX*zrp*zrp*(zrp*zrp + THREE*zrp*zrp*zrp))/((ONE + zrp**3)**3))
!
            dz1da=-THREE*(epsp + FIVE*P)/(pmarev**2)
            dz1dp=15._q*ALPHA/(pmarev**2)
!
            dzpda=dz1da
            dzpdp=dz1dp
!
            dwdr = dwdzp*(dz1dp*DPD + dz1da*DALPHAD)
            dwdg = dwdzp*(dz1dp*DPDD + dz1da*DALPHADD)
            dwdt = dwdzp*dz1da*DALPHADTAU
!
            dydp = (TWO*lamda - ONE)**2
            dfddy=(F727 + TWO*betha*y)/(TEN*(fds110**F910))
            dfdmedfd=(-TWO/(fd**3)) - ((28._q*rr)/(NINE*(fd**5)))
            dfdmedrr=SEVEN/(NINE*(fd**4))
            drrdp   = (40._q/NINE)*(F12 - lamda + lamda*lamda) - F53 &
                   &  + (F554*((TWO*lamda - ONE)**2))
            drrda   = -ONE + (THREE*(F12 - lamda + lamda*lamda))
!
            dfdmedr=dfdmedfd*dfddy*dydp*DPD + &
                  & dfdmedrr*(drrdp*DPD + drrda*DALPHAD)
            dfdmedg=dfdmedfd*dfddy*dydp*DPDD + &
                  & dfdmedrr*(drrdp*DPDD + drrda*DALPHADD)
            dfdmedt=dfdmedrr*drrda*DALPHADTAU
!
            dqgdp   =F23
            dqgda   =NITWE
            dfxscdp =(CMU + (100._q/729._q)*P)/(fxsc1**F910)
            dfxscdz1 = (((-73._q/675._q)*qg*(ONE-zp)) + &
                     &  ((73._q/675_q)*qg*zp))/(fxsc1**F910)
            dfxscdqg=(((292._q/2025._q)*qg) &
                     &  -((73._q/675._q)*zp*(ONE-zp)))/(fxsc1**F910)
!            
         dfxscdr=dfxscdp*DPD + dfxscdz1*(dzpdp*DPD + dzpda*DALPHAD) &
                &   +  dfxscdqg*(dqgdp*DPD + dqgda*DALPHAD)
         dfxscdg=dfxscdp*DPDD + dfxscdz1*(dzpdp*DPDD + dzpda*DALPHADD) &
                &   +  dfxscdqg*(dqgdp*DPDD + dqgda*DALPHADD)
         dfxscdt=DALPHADTAU*(dfxscdz1*dzpda + dfxscdqg*dqgda)
!
! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D   =dwdr*fdme + dfdmedr*w - dwdr*fxsc + dfxscdr*(ONE-w)
      DFX_DD  =dwdg*fdme + dfdmedg*w - dwdg*fxsc + dfxscdg*(ONE-w)
      DFX_DTAU=dwdt*fdme + dfdmedt*w - dwdt*fxsc + dfxscdt*(ONE-w)

!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD=EXDLDA*FX+EXLDA*DFX_D
      VXDD=EXLDA*DFX_DD
      AMUXD=EXLDA*DFX_DTAU
      EX=EX+EXLDA*FX/TWO

      RETURN
      END SUBROUTINE V2SRTMx

      SUBROUTINE V3SRTMx(RX,DRX,TAUX,Ex,VXD,VXDD,AMUXD)
GPACC !$ACC ROUTINE SEQ
!     Simplified regularized Tao-Mo variant 3 (v3-sregTM)
!     Functional (Exchange part only)
!
!     Written by:
!     H. Francisco
!     QTP, Department of Physics, University of Florida
!
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214102(2023)               
!
      USE prec
      USE constant

      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(V3SRTMx)

      REAL(q) :: RX,DRX,TAUX,EX,VXD,VXDD,AMUXD, CMU
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: PISQ,AX
      REAL(q) :: TAU_RHO,TAUW_RHO,TAU0,TAUW_TEMP
      REAL(q) :: RHO,DRHO,P,P2
      REAL(q) :: FA,DFADA,ALPHA,DALPHAD,DALPHADD,DALPHADTAU
      REAL(q) :: DPD,DPDD,EXUNIF,EXLDA,EXDLDA
      REAL(q) :: FX,DFXDP,DFXDA,DFX_D,DFX_DD,DFX_DTAU
! start Tao-Mo parameters and VARIABLES
      REAL(q) :: w, fdme, fxsc, zp, ftm, fd, rr, y, z2, z3, qg, z1,ftt
      REAL(q) :: lamda,betha, delta, con, F32, pma
      REAL(q) :: ftm1, ftm2,ftm3,F110,F35,F910,fxsc1, rar
      REAL(q) :: dwdr,dwdzp,dzpdp,dzpdftm,dftmdp,dzpda,dftmda
      REAL(q) :: dfdmedr,dfdmedfd,dfdmedrr,drrdr,dfddy,dydp
      REAL(q) :: drrdz2, dz2dp, drrdz3, dz3dr,dz3dz2,dz3da
      REAL(q) :: dfxscdr,dfxscdp,dfxscdqg,dqgdr,dfxscdz1,dz1dr
      REAL(q) :: dqgdp,dqgda,dz1dp,dz1da
      REAL(q) :: dwdg,dfdmedg,dfxscdg,dqgdg,dz1dg
      REAL(q) :: dwdl,dfdmedl,dfxscdl,lapq,laf
      REAL(q) :: dwdt,dfdmedt,dfxscdt,fds110
      REAL(q) :: drrdp, drrda,gam12
      REAL(q) :: F13, F23, F43, F53, F83, F18
      REAL(q) :: F727, F554, F12, NINE, NITWE,TEN
      REAL(q) :: F729,F2025, F405, TWENTY, EIGHT, FIF,T6
      REAL(q) :: zrp,ftm4,ftm5, F52
      REAL(q) :: z1a1, z1a2, z1a3, z1a4, z1a5, z1p1
      REAL(q) :: epsp, pmarev
!  
      PARAMETER (F43=4._q/3._q)
      PARAMETER (F13=1._q/3._q)
      PARAMETER (F83=8._q/3._q)
      PARAMETER (F23=2._q/3._q)
      PARAMETER (F53=5._q/3._q)
      PARAMETER (betha=79.873_q)
      PARAMETER (lamda=0.6866_q)
      PARAMETER (F110=1._q/10._q)
      PARAMETER (delta=1.475_q)
      PARAMETER (con=3._q)
      PARAMETER (F32=3._q/2._q)
      PARAMETER (F35=1._q/F53)
      PARAMETER (F910=9._q/10._q)
      PARAMETER (F727=700._q/27._q)
      PARAMETER (F554=595._q/54._q)
      PARAMETER (F729=50._q/729._q)
      PARAMETER (F2025=146._q/2025._q)
      PARAMETER (F405=73._q/405._q)
      PARAMETER (FIF=15._q)
      PARAMETER (T6=0.0000000001_q)
! finish Tao-MO parameters and varibles
      PARAMETER (CMU=10._q/81._q)
! other parameters
      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE=1._q)
      PARAMETER (TWO=2._q)
      PARAMETER (THREE=3._q)
      PARAMETER (FOUR=4._q)
      PARAMETER (FIVE=5._q)
      PARAMETER (SIX=6._q)
      PARAMETER (SEVEN=7._q)
      PARAMETER (EIGHT=8._q)
      PARAMETER (NINE= 9._q)
      PARAMETER (TEN=10._q)
      PARAMETER (TWENTY=20._q)
      PARAMETER (F18 = ONE/EIGHT)
      PARAMETER (F12 = ONE/TWO)
      PARAMETER (NITWE=NINE/TWENTY)
      PARAMETER (F52=FIVE/TWO)
      PARAMETER (THRD=1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.738558766382022405884230032680836_q)
      PARAMETER (epsp=0.58568_q)

      VXD=0._q
      VXDD=0._q
      AMUXD=0._q

!      rC1=(TWENTY*rK0/(TSEVEN*uge))**FOUR - (ONE + rE1)*(ONE + rE1)

! Suspect that TAUWU and TAUWD are not well described. Use
! TAUW_TEMP=0.125_q*DRT**2/RT
! IF WANT TO TEST TAUW, SIMPLY OVERWRITE TAUW_TEMP BY TAUW
      TAUW_TEMP=0.125_q*DRX**2._q/RX

! spin up
! IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RX
      DRHO=TWO*DRX
      TAUW_RHO=TWO*TAUW_TEMP
      TAU_RHO=TWO*TAUX
      TAU0=0.3_q*(THREE*PISQ)**THRD2*(RHO)**THRD5
!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4
! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      DPD=-THRD8*P/RHO
      DPDD=TWO*DRHO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)

      ALPHA = (TAU_RHO - TAUW_RHO)/TAU0
      if (ALPHA .LT. ZERO)  ALPHA=ZERO
      DALPHAD = THRD*(EIGHT*TAUW_RHO - FIVE*TAU_RHO)/(RHO*TAU0)
      DALPHADD = -0.25_q*DRHO/(RHO*TAU0)
      DALPHADTAU = ONE/TAU0

      P2 = P*P
      pmarev=FIVE*P + THREE*ALPHA + epsp
      zrp= (FIVE*P + epsp)/pmarev

      IF (zrp .GT. ONE) zrp=ONE
      zp=zrp
!
      w   = (zrp**2 + THREE*zrp**3)/((ONE + zrp**3)**2)
!
       z2  = F53*P
       z3  = z2 + ALPHA
       y   = P*((TWO*lamda - ONE)**2)
       fds110=ONE + F727*y + betha*y*y
       fd  = fds110**F110
       rr = ONE + F554*y - (z3 - THREE*((lamda**2) - &
        &   lamda + F12)*(z3 - ONE - (z2/NINE)))
       fdme=ONE/(fd**2) + SEVEN*rr/(NINE*(fd**4))
!
       qg  = NITWE*(ALPHA-ONE) + F23*P
      fxsc1= ONE + TEN*((CMU + F729*P)*P + F2025*qg*qg &
        &      - F405*qg*F35*zp*(ONE-zp))
      fxsc=fxsc1**F110
!
      FX = w*fdme + (ONE-w)*fxsc
!     Building derivatives of Tao-Mo-FX
      dwdzp= ((TWO*zrp + NINE*zrp*zrp)/((ONE + zrp**3)**2)) - &
       & ((SIX*zrp*zrp*(zrp*zrp + THREE*zrp*zrp*zrp))/((ONE + zrp**3)**3))
!
            dz1da=-THREE*(epsp + FIVE*P)/(pmarev**2)
            dz1dp=15._q*ALPHA/(pmarev**2)
!
            dzpda=dz1da
            dzpdp=dz1dp
!
            dwdr = dwdzp*(dz1dp*DPD + dz1da*DALPHAD)
            dwdg = dwdzp*(dz1dp*DPDD + dz1da*DALPHADD)
            dwdt = dwdzp*dz1da*DALPHADTAU
!
            dydp = (TWO*lamda - ONE)**2
            dfddy=(F727 + TWO*betha*y)/(TEN*(fds110**F910))
            dfdmedfd=(-TWO/(fd**3)) - ((28._q*rr)/(NINE*(fd**5)))
            dfdmedrr=SEVEN/(NINE*(fd**4))
            drrdp   = (40._q/NINE)*(F12 - lamda + lamda*lamda) - F53 &
                   &  + (F554*((TWO*lamda - ONE)**2))
            drrda   = -ONE + (THREE*(F12 - lamda + lamda*lamda))
!
            dfdmedr=dfdmedfd*dfddy*dydp*DPD + &
                  & dfdmedrr*(drrdp*DPD + drrda*DALPHAD)
            dfdmedg=dfdmedfd*dfddy*dydp*DPDD + &
                  & dfdmedrr*(drrdp*DPDD + drrda*DALPHADD)
            dfdmedt=dfdmedrr*drrda*DALPHADTAU
!
            dqgdp   =F23
            dqgda   =NITWE
            dfxscdp =(CMU + (100._q/729._q)*P)/(fxsc1**F910)
            dfxscdz1 = (((-73._q/675._q)*qg*(ONE-zp)) + &
                     &  ((73._q/675_q)*qg*zp))/(fxsc1**F910)
            dfxscdqg=(((292._q/2025._q)*qg) &
                     &  -((73._q/675._q)*zp*(ONE-zp)))/(fxsc1**F910)
!            
         dfxscdr=dfxscdp*DPD + dfxscdz1*(dzpdp*DPD + dzpda*DALPHAD) &
                &   +  dfxscdqg*(dqgdp*DPD + dqgda*DALPHAD)
         dfxscdg=dfxscdp*DPDD + dfxscdz1*(dzpdp*DPDD + dzpda*DALPHADD) &
                &   +  dfxscdqg*(dqgdp*DPDD + dqgda*DALPHADD)
         dfxscdt=DALPHADTAU*(dfxscdz1*dzpda + dfxscdqg*dqgda)
!
! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D   =dwdr*fdme + dfdmedr*w - dwdr*fxsc + dfxscdr*(ONE-w)
      DFX_DD  =dwdg*fdme + dfdmedg*w - dwdg*fxsc + dfxscdg*(ONE-w)
      DFX_DTAU=dwdt*fdme + dfdmedt*w - dwdt*fxsc + dfxscdt*(ONE-w)

!   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD=EXDLDA*FX+EXLDA*DFX_D
      VXDD=EXLDA*DFX_DD
      AMUXD=EXLDA*DFX_DTAU
      EX=EX+EXLDA*FX/TWO

      RETURN
      END SUBROUTINE V3SRTMx

      SUBROUTINE VRTMc(RU,RD,DRT,TAUU,TAUD,Ec,VCD1,VCD2,VCDD,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ
!     Revised regularized Tao-Mo (rregTM)
!     Functional (Correlation part only)
!
!     Written by:
!     H. Francisco
!     QTP, Department of Physics, University of Florida
!
!     References:              
!     S. Jana, S.K. Behera, S. Smiga, L.A. Constantine, and P. Samal
!     J. Chem. Phys. 155, 024103(2021)              
!
!     Note: The rregTM correlation is employed in conjunction 
!            with the three variants of the sregTM functional.
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214102(2023)               
!
      USE prec
      USE constant
      USE ldalib, only : PW91LDA

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(VRTMc)

      REAL(q) :: Ec,DRT,VCDD
      REAL(q) :: RU,TAUU,VCD1,AMUCD1
      REAL(q) :: RD,TAUD,VCD2,AMUCD2

      REAL(q) :: RT,Y,TAU
      REAL(q) :: TAUW,D_TAUW,D_TAUW_DD
      REAL(q) :: ZETA,OPZ,OMZ,D_ZETA_D1,D_ZETA_D2
      REAL(q) :: DS,D_DS_DZETA,D_DS_D1,D_DS_D2
      REAL(q) :: DX,D_DX_DZETA
      REAL(q) :: GC,D_GC_DZETA
      REAL(q) :: PHI,PHI2,PHI3,D_PHI_DZETA
      REAL(q) :: TAU0,D_TAU0_D1,D_TAU0_D2
      REAL(q) :: RS,RS12,D_RS
      REAL(q) :: P,D_P,D_P_DD
      REAL(q) :: AFIX_T,T2,D_T2_DP,D_T2_DRS,D_T2_DZETA
      REAL(q) :: ALPHA,ALPHA_DEN,D_ALPHA_D1,D_ALPHA_D2,D_ALPHA_DD,D_ALPHA_DT
      REAL(q) :: FA,D_FA_DA,D_FA_D1,D_FA_D2,D_FA_DD,D_FA_DT
      REAL(q) :: FACTOR0,D_FACTOR0_DRS,D2_FACTOR0_DRS
      REAL(q) :: ECLDA0,D_ECLDA0_DRS,D2_ECLDA0_DRS
      REAL(q) :: ECLDA1,D_ECLDA1_DRS,D_ECLDA1_DZETA,D2_ECLDA1_DRS,D2_ECLDA1_DRSDZETA
      REAL(q) :: DUMMY,ECTEMP,VCTEMPUP,VCTEMPDW
      REAL(q) :: BETA,D_BETA_DRS
      REAL(q) :: W1FAC,W1,D_W1_DRS,D_W1_DZETA
      REAL(q) :: A,D_A_DRS,D_A_DZETA
      REAL(q) :: V,D_V_DRS,D_V_DZETA,D_V_DP
      REAL(q) :: FACTOR1,EXP1,DELTAY,DELTAYP,D_DELTAYP_DRS,D_DELTAYP_DZETA,D_DELTAYP_DP
      REAL(q) :: GAUX,G,D_G_DGAUX,D_G_DRS,D_G_DZETA,D_G_DP
      REAL(q) :: ARG1,D_ARG1_DP,D_ARG1_DZETA,D_ARG1_DRS
      REAL(q) :: H1,D_H1_DRS,D_H1_DP,D_H1_DZETA
      REAL(q) :: EC1,D_EC1_DRS,D_EC1_DP,D_EC1_DZETA,D_EC1_D1,D_EC1_D2,D_EC1_DD
      REAL(q) :: W0,D_W0_DRS
      REAL(q) :: GINF,D_GINF_DP
      REAL(q) :: HCORE0,H0,D_H0_DRS,D_H0_DP
      REAL(q) :: EC0,D_EC0_DRS,D_EC0_DP,D_EC0_DZETA,D_EC0_D1,D_EC0_D2,D_EC0_DD
      REAL(q) :: EPP,D_EPP_D1,D_EPP_D2,D_EPP_DD,D_EPP_DT

      REAL(q) :: BETA_MB,GAMMA,FACTORA,FACTORB,B1C,B2C,B3C,XI
      REAL(q) :: A0,A1,A2,A3,A4,A5,A6,A7
      REAL(q) :: CFC1,CFC2,CFD1
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,PISQ
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: ETA,CKF2,DP24,C2
!     VARIABLES TO REG TAO-MO
      REAL(q) :: gf2, gamc2, df2dg, dgf2da, A2t4
      REAL(q) :: dgat2da, dgat2dt2,fba,gta22
      REAL(q) :: gamc, NINE
      PARAMETER (gamc2 = 0.2_q)
      PARAMETER (gamc  = 1.2_q)
      PARAMETER (NINE  = 9._q)
!     FINISH VARIABLES REG-TAO-MO
      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE=1.0_q)
      PARAMETER (TWO=2.0_q)
      PARAMETER (THREE=3.0_q)
      PARAMETER (FOUR=4.0_q)
      PARAMETER (FIVE=5.0_q)
      PARAMETER (SIX=6.0_q)
      PARAMETER (SEVEN=7.0_q)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1.0_q/3.0_q)
      PARAMETER (THRD2=2.0_q*THRD)
      PARAMETER (THRD4=4.0_q*THRD)
      PARAMETER (THRD5=1.0_q+THRD2)
      PARAMETER (THRD8=1.0_q+THRD5)
      PARAMETER (CFC1=0.640_q)
      PARAMETER (CFC2=1.500_q)
      PARAMETER (CFD1=0.700_q)
      parameter( a0 =  1.0_q )
      parameter( a1 = -0.64_q )
      parameter( a2 = -0.4352_q )
      parameter( a3 = -1.535685604549_q )
      parameter( a4 =  3.061560252175_q )
      parameter( a5 = -1.915710236206_q )
      parameter( a6 =  0.516884468372_q )
      parameter( a7 = -0.051848879792_q )
      PARAMETER(DP24=0.361_q**4)
      PARAMETER(B1C=0.0285764_q,B2C=0.0889_q,B3C=0.125541_q)
      PARAMETER(BETA_MB=0.06672455060314922_q)
      PARAMETER(FACTORA=0.1_q,FACTORB=0.1778_q)
      PARAMETER(XI=0.12802585262625815_q)
      PARAMETER(ETA=0.001_q)
      PARAMETER(C2=A1 + TWO*A2 + THREE*A3 + FOUR*A4 + FIVE*A5 + SIX*A6 + SEVEN*A7)
      PARAMETER(GAMMA=(ONE-DLOG(TWO))/PISQ)
      PARAMETER(CKF2=(THREE*PISQ)**THRD2)
GPACC !$ACC ROUTINE(PW91LDA) SEQ

      RT=RU+RD
      Y=DRT**2
      TAU=TAUU+TAUD

      TAUW=1.0d0/8.0d0*(Y/RT)
      D_TAUW = -TAUW/RT
      D_TAUW_DD = 0.25_q*DRT/RT

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-ONE),ONE)
      D_ZETA_D1= TWO*RD/RT**2
      D_ZETA_D2=-TWO*RU/RT**2

      OMZ = ONE - ZETA
      OPZ = ONE + ZETA

      DS=(OPZ**THRD5+OMZ**THRD5)/TWO
      D_DS_DZETA = THRD5*(OPZ**THRD2 - OMZ**THRD2)/TWO
      D_DS_D1=D_DS_DZETA*D_ZETA_D1
      D_DS_D2=D_DS_DZETA*D_ZETA_D2

      PHI = (OPZ**THRD2 + OMZ**THRD2)/TWO
      PHI2 = PHI*PHI
      PHI3 = PHI2*PHI
      D_PHI_DZETA = 0.0_q
      IF (OPZ.GT.1E-10_q) D_PHI_DZETA = THRD/OPZ**THRD
      IF (OMZ.GT.1E-10_q) D_PHI_DZETA = D_PHI_DZETA - THRD/OMZ**THRD

      DX = (OPZ**THRD4 + OMZ**THRD4)/TWO
      D_DX_DZETA = THRD2*(OPZ**THRD - OMZ**THRD)

      GC = (ONE - 2.363D0*(DX - ONE))*(ONE - ZETA**12)
      D_GC_DZETA = -12D0*(ONE - 2.363D0*(DX - ONE))*ZETA**11 - &
                 &  2.363D0*D_DX_DZETA*(ONE - ZETA**12)


      TAU0=0.30d0*(THREE*PISQ)**THRD2*RT**THRD5*DS
      D_TAU0_D1 = THRD5*TAU0/RT + TAU0*D_DS_D1/DS
      D_TAU0_D2 = THRD5*TAU0/RT + TAU0*D_DS_D2/DS

      RS    = (0.750d0/PI)**THRD/RT**THRD
      RS12  = DSQRT(RS)
      D_RS  = -THRD/RT*RS

      P = Y/(FOUR*CKF2*RT**THRD8)
      D_P = -THRD8*P/RT
      D_P_DD = DRT/(TWO*CKF2*RT**THRD8)
      AFIX_T = CKF2/FOUR**THRD4
      T2 = AFIX_T*P/(PHI2*RS)
      D_T2_DRS = -T2/RS
      D_T2_DP = AFIX_T/(PHI2*RS)
      D_T2_DZETA = -TWO*T2/PHI*D_PHI_DZETA

!     ALPHA FUNCTION FOR REGULARIZED TAO-MO AND SWITCHING FUNCTION FA AND DFADA
      ALPHA_DEN = TAU0
      ALPHA = (TAU - TAUW)/ALPHA_DEN
      IF (ALPHA .LT. ZERO)  ALPHA=ZERO
!
      D_ALPHA_D1 =(-D_TAUW - ALPHA*D_TAU0_D1)/ALPHA_DEN
      D_ALPHA_D2 =(-D_TAUW - ALPHA*D_TAU0_D2)/ALPHA_DEN
      D_ALPHA_DD =-D_TAUW_DD/ALPHA_DEN
      D_ALPHA_DT = ONE/ALPHA_DEN
!
      gf2 = gamc*ALPHA/(gamc2 + ALPHA)
      FA  = (THREE*gf2**3)/(ONE + gf2**3 + gf2**6)
!
      df2dg =(NINE*gf2**2)*(ONE - gf2**6)/((ONE + gf2**3 + gf2**6)**2)
      dgf2da = gamc*gamc2/((gamc2 + ALPHA)**2)
!
      D_FA_DA= df2dg*dgf2da
!FINISH ALPHA FUNCTION FOR REGULARIZED TAO-MO AND SWITCHING FUNCTION FA AND DFADA
!
      D_FA_D1 = D_FA_DA*D_ALPHA_D1
      D_FA_D2 = D_FA_DA*D_ALPHA_D2
      D_FA_DD = D_FA_DA*D_ALPHA_DD
      D_FA_DT = D_FA_DA*D_ALPHA_DT
      !--------------------------------------------------------------------------
      !EC1
      !--------------------------------------------------------------------------

      FACTOR0 = ONE + B2C*RS12 + B3C*RS
      D_FACTOR0_DRS = B2C/RS12/TWO + B3C
      D2_FACTOR0_DRS = -B2C/RS/RS12/FOUR

      ECLDA0 = -B1C/FACTOR0
      D_ECLDA0_DRS = D_FACTOR0_DRS*B1C/FACTOR0**2
      D2_ECLDA0_DRS = B1C/FACTOR0**2 * (D2_FACTOR0_DRS  - TWO*D_FACTOR0_DRS**2/FACTOR0)

      CALL PW91LDA(RS,ZETA,ECLDA1,D_ECLDA1_DRS,D_ECLDA1_DZETA,D2_ECLDA1_DRS,D2_ECLDA1_DRSDZETA)

      BETA = BETA_mb*(ONE + FACTORA*RS)/(ONE + FACTORB*RS)
      D_BETA_DRS = BETA_mB*(FACTORA - FACTORB)/(ONE + FACTORB*RS)**2

      W1FAC = ECLDA1/(GAMMA*PHI3)
      W1 = DEXP(-W1FAC) - ONE
      D_W1_DRS = -(W1+ONE)*D_ECLDA1_DRS/(GAMMA*PHI3)
      D_W1_DZETA = (W1+ONE)*(THREE*W1FAC*D_PHI_DZETA/PHI - D_ECLDA1_DZETA/(GAMMA*PHI3))

      A = BETA/(W1*GAMMA)
      D_A_DRS = D_BETA_DRS/(W1*GAMMA) - A*D_W1_DRS/W1
      D_A_DZETA = -A*D_W1_DZETA/W1

!
      V = A*T2 + (A*T2)**2
      D_V_DRS = D_A_DRS*T2 + A*D_T2_DRS + TWO*A*T2*(D_A_DRS*T2 + A*D_T2_DRS)
      D_V_DZETA = D_A_DZETA*T2 + A*D_T2_DZETA + TWO*A*T2*(D_A_DZETA*T2 + A*D_T2_DZETA)
      D_V_DP = A*D_T2_DP + TWO*A*T2*A*D_T2_DP
!
      GAUX = ONE + V
      G    = ONE/GAUX
      gta22=-G**2
!
      D_G_DRS =   gta22*D_V_DRS
      D_G_DZETA = gta22*D_V_DZETA
      D_G_DP =    gta22*D_V_DP
!
      ARG1 = ONE + W1*(ONE-G)
      D_ARG1_DP = -W1*D_G_DP
      D_ARG1_DRS = D_W1_DRS*(ONE-G) - W1*D_G_DRS
      D_ARG1_DZETA = D_W1_DZETA*(ONE-G) - W1*D_G_DZETA

      H1 = GAMMA*PHI3*DLOG(ARG1)
      D_H1_DRS = GAMMA*PHI3*D_ARG1_DRS/ARG1
      D_H1_DP = GAMMA*PHI3*D_ARG1_DP/ARG1
      D_H1_DZETA = THREE*H1*D_PHI_DZETA/PHI + GAMMA*PHI3*D_ARG1_DZETA/ARG1

      EC1 = ECLDA1 + H1
      D_EC1_DRS = D_ECLDA1_DRS + D_H1_DRS
      D_EC1_DP = D_H1_DP
      D_EC1_DZETA = D_ECLDA1_DZETA + D_H1_DZETA

      D_EC1_D1 = D_EC1_DRS*D_RS + D_EC1_DP*D_P + D_EC1_DZETA*D_ZETA_D1
      D_EC1_D2 = D_EC1_DRS*D_RS + D_EC1_DP*D_P + D_EC1_DZETA*D_ZETA_D2
      D_EC1_DD = D_EC1_DP*D_P_DD

      !-----------------------------------------------------------------------------
      ! EC0
      !-----------------------------------------------------------------------------
      W0 = DEXP(-ECLDA0/B1C) - ONE
      D_W0_DRS = -(W0+ONE)*D_ECLDA0_DRS/B1C

      GINF = ONE/(ONE + FOUR*XI*P)**(0.25_q)
      D_GINF_DP = -XI*GINF**5

      HCORE0 = ONE + W0*(ONE - GINF)
      H0 = B1C*dlog(HCORE0)
      D_H0_DRS  = B1C*D_W0_DRS*(ONE - GINF)/HCORE0
      D_H0_DP  = -B1C*W0*D_GINF_DP/HCORE0

      EC0 = (ECLDA0 + H0)*GC
      D_EC0_DRS = GC*(D_ECLDA0_DRS + D_H0_DRS)
      D_EC0_DP = GC*D_H0_DP
      D_EC0_DZETA = D_GC_DZETA*(ECLDA0 + H0)

      D_EC0_D1 = D_EC0_DRS*D_RS + D_EC0_DZETA*D_ZETA_D1 + D_EC0_DP*D_P
      D_EC0_D2 = D_EC0_DRS*D_RS + D_EC0_DZETA*D_ZETA_D2 + D_EC0_DP*D_P
      D_EC0_DD = D_EC0_DP*D_P_DD

      !-----------------------------------------------------------------------------

      EPP = EC0 + FA*(EC1 - EC0)
      EC=RT*EPP

      D_EPP_D1 = D_EC0_D1 + D_FA_D1*(EC1-EC0) + FA*(D_EC1_D1 - D_EC0_D1)
      D_EPP_D2 = D_EC0_D2 + D_FA_D2*(EC1-EC0) + FA*(D_EC1_D2 - D_EC0_D2)
      D_EPP_DD = D_EC0_DD + D_FA_DD*(EC1-EC0) + FA*(D_EC1_DD - D_EC0_DD)
      D_EPP_DT = D_FA_DT*(EC1-EC0)

      VCD1 = EPP + RT*D_EPP_D1
      VCD2 = EPP + RT*D_EPP_D2
      VCDD = RT*D_EPP_DD
      AMUCD1 = RT*D_EPP_DT
      AMUCD2 = AMUCD1

      RETURN
      END SUBROUTINE VRTMc

      SUBROUTINE XV2SRTML(RHO,RHOTHRD,S,QQ,EXLDA,EX,EXDLDA,EXD,EXDD,EXDL)
GPACC !$ACC ROUTINE SEQ
!     Deorbitalized Simplified regularized Tao-Mo-L variant 2 (v2-sregTM-L)
!     Functional (Exchange part only)
!
!     Written by:
!     H. Francisco
!     QTP, Department of Physics, University of Florida
!
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214103(2023)               
!
              
      USE prec
      USE constant

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(XV2SRTML)

      REAL(q) :: RHO,RHOTHRD,S,QQ,EXLDA,EX,EXDLDA,EXD,EXDD,EXDL
      REAL(q) :: KF,EXUNIF
      REAL(q) :: rhoval, rrho, rho13, rho43, rho53, rho83
      REAL(q) :: Pconst
      REAL(q) :: fs, dfsdp, dfsdq
      REAL(q) :: derivr, derivg, derivl
      REAL(q) :: dadp, dadq, dpdg, dpdr, dqdr, dqdl
      REAL(q) :: dFxda, dFxdp, dFxdr, dFxdg, dFxdl
      REAL(q) :: dadr, dadg, dadl
      REAL(q) :: p, ql, Fx, a, Fa
!     Variables for TAO-MO exchange functional 
      REAL(q) ::  w, fdme, fxsc
      REAL(q) :: zp, ftm, fd, rr, y, z2, z3, qg, z1,ftt, pma
      REAL(q) :: ftm1, ftm2,ftm3,fxsc1, rar
!     variables for first derivative with respect to density of TAO-M0
      REAL(q) :: dwdr,dwdzp,dzpdp,dzpdftm,dftmdp,dzpda,dftmda
      REAL(q) :: dfdmedr,dfdmedfd,dfdmedrr,drrdr,dfddy,dydp
      REAL(q) :: drrdz2, dz2dp, drrdz3, dz3dr,dz3dz2,dz3da
      REAL(q) :: dfxscdr,dfxscdp,dfxscdqg,dqgdr,dfxscdz1,dz1dr
      REAL(q) :: dqgdp,dqgda,dz1dp,dz1da
!     variables for first derivative with respect to gradient of TAO-M0
      REAL(q) :: dwdg,dfdmedg,dfxscdg,dqgdg,dz1dg
!     variables for first derivative with respect to laplacian of TAO-M0
      REAL(q) :: dwdl,dfdmedl,dfxscdl,lapq,laf
      REAL(q) :: dwdt,dfdmedt,dfxscdt,fds110
      REAL(q) :: drrdp, drrda,gam12
!     variables for rregtm
      REAL(q) :: zrp,ftm4,ftm5
      REAL(q) :: z1a1, z1a2, z1a3, z1a4, z1a5, z1p1
!
      REAL(q),PARAMETER :: betha=79.8731_q
      REAL(q),PARAMETER :: lamda=0.68661_q
      REAL(q),PARAMETER :: F110=1.0_q/10.0_q
      REAL(q),PARAMETER :: delta=1.475_q
      REAL(q),PARAMETER :: con=3.0_q
      REAL(q),PARAMETER :: F32=3.0_q/2.0_q
      REAL(q),PARAMETER :: F910=9.0_q/10.0_q
      REAL(q),PARAMETER :: F12=1.0_q/2.0_q
      REAL(q),PARAMETER :: TOL2=0.0000000001_q
!
      REAL(q),PARAMETER :: F43=4.0_q/3.0_q
      REAL(q),PARAMETER :: F13=1.0_q/3.0_q
      REAL(q),PARAMETER :: F83=8.0_q/3.0_q
      REAL(q),PARAMETER :: F23=2.0_q/3.0_q
      REAL(q),PARAMETER :: F18=1.0_q/8.0_q
      REAL(q),PARAMETER :: F53=5.0_q/3.0_q
      REAL(q),PARAMETER :: F14=1.0_q/4.0_q
      REAL(q),PARAMETER :: F35=1.0_q/F53
      REAL(q),PARAMETER :: cero=0.0_q
      REAL(q),PARAMETER :: ONE=1.0_q
      REAL(q),PARAMETER :: two=2.0_q
      REAL(q),PARAMETER :: three=3.0_q
      REAL(q),PARAMETER :: four=4.0_q
      REAL(q),PARAMETER :: five=5.0_q
      REAL(q),PARAMETER :: F52=five/two
!
      REAL(q),PARAMETER :: CKF=(3.0_q*PI*PI)**F13
      REAL(q),PARAMETER :: AX=-3.0_q/4.0_q*(3.0_q/PI)**F13
      REAL(q) :: dpds, dsdng, vfer
      REAL(q) :: epsp, pmarev
      PARAMETER (epsp=0.5_q)
GPACC !$ACC ROUTINE(ts_pcrep) SEQ
!
      Pconst=(3.0_q*PI*PI)**F23
      vfer  =(3.0_q*PI*PI)**F13
!
      KF = CKF*RHOTHRD
!
      EXUNIF = AX*RHOTHRD
      EXLDA = EXUNIF*RHO
      EXDLDA = F43*EXUNIF
!
            rhoval= RHO
            rho43 = rhoval**F43
            rrho  = ONE/rhoval
            rho13 = rho43*rrho
            rho83 = rho43*rho43
            rho53 = rho43*rho13
!
      p = S*S
      ql = QQ
      CALL ts_pcrep(RHO, p, ql, dfsdp, dfsdq, fs)
!======================================================================
      a=fs-F53*p
      if(a.lt.cero)  a=cero
!      
       pmarev=FIVE*P + THREE*a + epsp
       zrp= (FIVE*P + epsp)/pmarev
!
       if (zrp .gt. ONE) zrp=ONE
       zp=zrp
!
       w   = (zrp**2 + three*zrp**3)/((ONE + zrp**3)**2)
!
       z2  = F53*p
       z3  = z2 + a
       y   = p*((two*lamda - ONE)**2)
       fds110=ONE + (700.0_q/27.0_q)*y + betha*y*y
       fd  = fds110**F110
       rr = ONE + (595.0_q/54.0_q)*y - (z3 - three*((lamda**2) - &
       &  lamda + F12)*(z3 - ONE - (z2/9.0_q)))
       fdme=ONE/(fd**2) + 7.0_q*rr/(9.0_q*(fd**4))
!
      qg  = (9.0_q/20.0_q)*(a-ONE) + F23*p
      fxsc1= ONE + 10.0_q*(((10.0_q/81.0_q) + (50.0_q/729.0_q)*p)*p + &
      & (146.0_q/2025.0_q)*qg*qg - (73.0_q/405.0_q)*qg*F35*zp*(ONE-zp))
      fxsc=fxsc1**F110

      Fx = w*fdme + (ONE-w)*fxsc
      EX = EXLDA*Fx

!======================================================================
           dpds  = two*s
           dsdng = ONE/(two*vfer*rho43)

           dpdr = -F83*p*rrho
           dpdg = ONE/(four*Pconst*rho83)

           dqdr = -F53*ql*rrho
           dqdl = ONE/(four*Pconst*rho53)

           dadp = dfsdp - F53
           dadq = dfsdq

           dadr = dadp*dpdr + dadq*dqdr
           dadg = dadp*dpds*dsdng
           dadl = dadq*dqdl
!=========================================================================
            dwdzp= ((2.0_q*zrp + 9.0_q*zrp*zrp)/((ONE + zrp**3)**2)) - &
           & ((6.0_q*zrp*zrp*(zrp*zrp + 3.0_q*zrp**3))/((ONE + zrp**3)**3))
!
            dz1da=-THREE*(epsp + FIVE*P)/(pmarev**2)
!
            dz1dp=15._q*a/(pmarev**2)
!
            dzpda=dz1da
            dzpdp=dz1dp
!
            dwdr = dwdzp*(dz1dp*dpdr + dz1da*dadr)
            dwdg = dwdzp*(dz1dp*dpds*dsdng + dz1da*dadg)
            dwdl = dwdzp*dz1da*dadl
!
            dydp = (two*lamda - ONE)**2
            dfddy=((700.0_q/27.0_q) + two*betha*y)/(10.0_q*(fds110**F910))
            dfdmedfd=(-two/(fd**3)) - ((28.0_q*rr)/(9.0_q*(fd**5)))
            dfdmedrr=7.0_q/(9.0_q*(fd**4))
            drrdp   = (40.0_q/9.0_q)*(F12 - lamda + lamda*lamda) - F53 + &      
                   & ((595.0_q/54.0_q)*((2.0_q*lamda - 1.0_q)**2))
            drrda   = -ONE + (three*(F12 - lamda + lamda*lamda))
!
            dfdmedr=dfdmedfd*dfddy*dydp*dpdr + dfdmedrr*(drrdp*dpdr + drrda*dadr)
            dfdmedg=dfdmedfd*dfddy*dydp*dpds*dsdng + dfdmedrr*(drrdp*dpds*dsdng + drrda*dadg)
            dfdmedl=dfdmedrr*drrda*dadl
!
            dqgdp   =F23
            dqgda   =9.0_q/20.0_q
            dfxscdp =((10.0_q/81.0_q) + (100.0_q/729.0_q)*p)/(fxsc1**F910)
            dfxscdz1 = (((-73.0_q/675.0_q)*qg*(ONE-zp)) + &
                     &  ((73.0_q/675.0_q)*qg*zp))/(fxsc1**F910)
            dfxscdqg=(((292.0_q/2025.0_q)*qg) - &
                     &  ((73.0_q/675.0_q)*zp*(ONE-zp)))/(fxsc1**F910)
!            
            dfxscdr=dfxscdp*dpdr + dfxscdz1*(dzpdp*dpdr + dzpda*dadr) &
                &   +  dfxscdqg*(dqgdp*dpdr + dqgda*dadr)
            dfxscdg=dfxscdp*dpds*dsdng + dfxscdz1*(dzpdp*dpds*dsdng + dzpda*dadg) &
                &   +  dfxscdqg*(dqgdp*dpds*dsdng + dqgda*dadg)
            dfxscdl=dadl*(dfxscdz1*dzpda + dfxscdqg*dqgda)
!
            dFxdr = dwdr*fdme + dfdmedr*w - dwdr*fxsc + dfxscdr*(ONE-w)
            dFxdg = dwdg*fdme + dfdmedg*w - dwdg*fxsc + dfxscdg*(ONE-w)
            dFxdl = dwdl*fdme + dfdmedl*w - dwdl*fxsc + dfxscdl*(ONE-w)
!======================================================================
!DERIVATIVES
!======================================================================

      EXD  = F43*EXUNIF*Fx + EXLDA*dFxdr
      EXDD = EXLDA*dFxdg
      EXDL = EXLDA*dFxdl

      return
      END SUBROUTINE XV2SRTML

      SUBROUTINE CRTML(RS,ZETA,EC,VCUP,VCDW,G,FK,S,H,DVCUP,DVCDW,ECDD,TAUTOT,ECDT)
GPACC !$ACC ROUTINE SEQ
!     Deorbitalized revised regularized Tao-Mo-L (rregTM-L)
!     Functional (Correlation part only)
!
!     Written by:
!     H. Francisco
!     QTP, Department of Physics, University of Florida
!
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214103(2023)               
!
      USE prec
      USE constant
      USE ldalib, only : PW91LDA

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(CRTML)

      INTEGER i
      REAL(q) rs,zeta,ec,vcup,vcdw,g,fk,p,h,dvcup,dvcdw,ecdd,tautot,ecdt,dum,vc,dvc
      REAL(q) decdzeta,decdrs,opz,omz,phi,phi2,phi3,beta,dbetadrz
      REAL(q) w1fac,expw1,w1,A,dadrs,dw1drs,t2,dt2dp,dt2drs
      REAL(q) at2,gaux,gat2,arg1,darg1drs,darg1dp,h1,dh1dp,dh1drs
      REAL(q) ec1,dec1drs,dec1dp,dx,gc,w0fac,expw0,w0,dw0drs
      REAL(q) ginf,dginfdp,arg0,darg0drs,darg0dp,ec0,dec0drs,dec0dp,h0,dh0dp
      REAL(q) dh0drs,ds,tueg,tvw,alpha,oma,oma2,exp5,exp6,fca,dfcada
      REAL(q) dpdn,dpdg,drsdn,dalphadn,dalphadg,dalphadt,dH1dn,dEc0dn,dEc1dn
      REAL(q) dfcadn,dfcadg,dfcadt,dEc0dg,dEc1dg,dphidzeta,dt2dzeta,dw1dzeta
      REAL(q) dadzeta,darg1dzeta,dH1dzeta,dEc1dzeta,ddxdzeta,dgcdzeta,dEc0dzeta
      REAL(q) ddsdzeta,dalphadzeta,dfcadzeta,vcpol,dbetadrs,d,rs12,s
      REAL(q) dalphadp,eclda1,declda1drs,eclda0,declda0drs,d2eclda0drs2
      REAL(q) deltay,deltayp,ddeltaypdp,ddeltaypdrs,gat2y,dgat2ydp,dgat2ydrs,exp1
      REAL(q) ddeltaypdz,dgat2ydz,vcn,dd,d2eclda1drs2,d2eclda1drsdz,declda1dz
      REAL(q) dalphadtvw
!
!     VARIABLES TO REG TAO-MO 
      REAL(q) :: gf2, gamc2, df2dg, dgf2da, A2t4
      REAL(q) :: dgat2da, dgat2dt2,fba,gta22
      REAL(q) :: gamc, NINE, THREE, ONE
      PARAMETER (gamc2 = 0.2_q)
      PARAMETER (gamc  = 1.2_q)
      PARAMETER (NINE  = 9._q)
      PARAMETER (ONE  = 1._q)
      PARAMETER (THREE  = 3._q)
!     FINISH VARIABLES REG-TAO-MO
      REAL(q),PARAMETER :: p14a = 0.1_q
      REAL(q),PARAMETER :: p14b = 0.1778_q
      REAL(q),PARAMETER :: gamma = 0.03109069086965489503494086371273_q
      REAL(q),PARAMETER :: betazero = 0.06672455060314922_q
      REAL(q),PARAMETER :: b1c = 0.0285764_q
      REAL(q),PARAMETER :: b2c = 0.0889_q
      REAL(q),PARAMETER :: b3c = 0.125541_q
      REAL(q),PARAMETER :: c1c = 0.64_q
      REAL(q),PARAMETER :: c2c = 1.5_q
      REAL(q),PARAMETER :: dc = 0.7_q
      REAL(q),PARAMETER :: dxc = 2.3631_q
      REAL(q),PARAMETER :: xi = 0.12802585262625815_q
      REAL(q),PARAMETER :: thrd = 1.0_q/3.0_q
      REAL(q),PARAMETER :: thrd2 = 2.0_q/3.0_q
      REAL(q),PARAMETER :: thrd4 = 4.0_q/3.0_q
      REAL(q),PARAMETER :: thrd5 = 5.0_q/3.0_q
      REAL(q),PARAMETER :: ckf = (3.0_q*pi*pi)**thrd
      REAL(q),PARAMETER :: ckf2=ckf*ckf

      REAL(q) :: eta, two53
      parameter (two53=2d0**thrd5)
      Parameter (eta=0.0_q)
GPACC !$ACC ROUTINE(PW91LDA) SEQ

      d = 3.0_q/(4.0_q*pi) / rs**3
      rs12 = sqrt(rs)
      p = s*s
      dd = 2.0*s*d*fk

      opz = 1.0_q + zeta
      omz = 1.0_q - zeta
      phi = G
      phi2 = phi*phi
      phi3 = phi2*phi

!     ------------------------------------------------------------------
!     Kinetic energy density
      ds = 0.5_q*(opz**thrd5 + omz**thrd5)
      tueg = 0.3_q*fk*fk*ds*d
      tvw = d*p*fk*fk/2.0_q

      alpha = (tautot - tvw)/(tueg + eta*tvw)
      IF (alpha .LT. 0.0_q)  alpha=0.0_q
      dalphadt = 1.0_q/(tueg + eta*tvw)
      dalphadtvw = -1.0_q/(tueg + eta*tvw) - alpha*eta/(tueg + eta*tvw)
      dalphadp = dalphadtvw*thrd5*tueg/ds
      dalphadn = dalphadtvw*thrd5*thrd5*tueg*p/ds - alpha/(tueg+eta*tvw)*thrd5*tueg
!     ALPHA FUNCTION FOR REGULARIZED TAO-MO AND SWITCHING FUNCTION FA
!     AND DFADA
      gf2 = gamc*alpha/(gamc2 + alpha)
      fba = (THREE*gf2**3)/(ONE + gf2**3 + gf2**6)
!
      df2dg =(NINE*gf2**2)*(ONE - gf2**6)/((ONE + gf2**3 + gf2**6)**2)
      dgf2da = gamc*gamc2/((gamc2 + alpha)**2)
!
      fca   = ONE - fba
      dfcada= -df2dg*dgf2da
!FINISH ALPHA FUNCTION FOR REGULARIZED TAO-MO AND SWITCHING FUNCTION FA
!AND DFADA
!     -------------------------------------------------------------------
!
      eclda0 = -b1c/(1.0_q + b2c*rs12 + b3c*rs)
      declda0drs = -eclda0*(b3c + 0.5_q*b2c/rs12)/&
     &                (1.0_q + b2c*rs12 + b3c*rs)
      d2eclda0drs2 = eclda0/(1.0_q + b2c*rs12 + b3c*rs) *&
     &          ( 2*(b3c+0.5_q*b2c/rs12)**2/(1.0_q+b2c*rs12+b3c*rs) +&
     &            0.25_q*b2c/rs/rs12 )

      CALL PW91LDA(rs,zeta,eclda1,declda1drs,declda1dz,d2eclda1drs2,d2eclda1drsdz)

      BETA = BETAzero*(1.0_q + p14a*rs)/(1.0_q + p14b*rs)
      dBETAdrs = BETAzero*(p14a - p14b)/(1.0_q + p14b*rs)**2

      w1fac = eclda1/(phi3*GAMMA)
      expw1 = dexp(-w1fac)
      w1 = expw1 - 1.0_q
      dw1drs = -expw1*declda1drs/(GAMMA*phi3)

      A = BETA/GAMMA/w1
      dAdrs = dBETAdrs/GAMMA/w1 - A*dw1drs/w1

      t2 = ckf2*p/(4.0_q*rs*phi2*2.0_q**thrd2)
      dt2dp = ckf2/(4.0_q*rs*phi2*2.0_q**thrd2)
      dt2drs = -t2/rs

      At2=A*t2 + (A*t2)**2
      Gaux = 1d0 + At2
      GAt2 = 1d0/Gaux
      gta22=GAt2**2
      dgat2da = t2 + 2.d0*A*t2*t2
      dgat2dt2= A + 2.d0*A*A*t2


      arg1 = ONE + w1*(ONE - GAt2)
      darg1drs = dw1drs*(ONE - GAt2) + gta22*w1*(dgat2da*dAdrs + dgat2dt2*dt2drs)
      darg1dp = gta22*w1*dgat2dt2*dt2dp

      H1 = GAMMA*phi3*dlog(arg1)
      dH1dp = GAMMA*phi3*darg1dp/arg1
      dH1drs = GAMMA*phi3*darg1drs/arg1

      Ec1 = eclda1 + H1
      dEc1drs = declda1drs + dH1drs
      dEc1dp = dH1dp
!
!     --------------------------------------------------------------
!
      w0fac = eclda0/b1c
      expw0 = dexp(-w0fac)
      w0 = expw0 - 1.0_q
      dw0drs = -declda0drs*expw0/b1c

      dx = 0.5_q*(opz**thrd4 + omz**thrd4)
      gc = (1.0_q - dxc*(dx - 1.0_q))*(1.0_q - zeta**12)

      ginf = (1.0_q/(1.0_q + 4.0_q*xi*p))**0.25_q
      dginfdp = -xi*ginf/(1.0_q + 4.0_q*xi*p)

      arg0 = 1.0_q + w0*(1.0_q - ginf)
      darg0drs = dw0drs*(1.0_q - ginf)
      darg0dp = -w0*dginfdp
      H0 = b1c*dlog(arg0)
      dH0dp = b1c*darg0dp/arg0
      dH0drs = b1c*darg0drs/arg0

      Ec0 = gc*(eclda0 + H0)
      dEc0drs = gc*(declda0drs + dH0drs)
      dEc0dp = gc*dH0dp
!
!     --------------------------------------------------------------
!
      H = H1 + fca*(Ec0-Ec1)
!
      dpdn = -8.0_q/3.0_q*p
      dpdg = s/fk
      drsdn = -rs/3.0_q

      dH1dn = dH1drs*drsdn + dH1dp*dpdn
      dEc1dn = dEc1drs*drsdn + dEc1dp*dpdn
      dEc0dn = dEc0drs*drsdn + dEc0dp*dpdn
      dfcadn = dfcada*dalphadn + dfcada*dalphadp*dpdn
      dfcadg = dfcada*dalphadp*dpdg
      dfcadt = dfcada*dalphadt

      ec = eclda1
      vcn = H1 + fca*(Ec0-Ec1) + (dH1dn + fca*(dEc0dn-dEc1dn) + dfcadn*(Ec0-Ec1))

      dEc1dg = dEc1dp*dpdg
      dEc0dg = dEc0dp*dpdg

      ecdd =  (dEc1dg + fca*(dEc0dg-dEc1dg) + dfcadg*(Ec0-Ec1))
      ecdt =  d*dfcadt*(Ec0-Ec1)

!     -----------------------------------------------------------------

      dphidzeta = 0.0_q
      if (opz.gt.1e-10_q) dphidzeta = dphidzeta + thrd/opz**thrd
      if (omz.gt.1e-10_q) dphidzeta = dphidzeta - thrd/omz**thrd

      dt2dzeta = -2.0_q*t2*dphidzeta/phi
      dw1dzeta = (3.0_q*w1fac*dphidzeta/phi-declda1dz/(gamma*phi3))*expw1
      dAdzeta = -A*dw1dzeta/w1
      ddsdzeta = 0.5_q*thrd5*(opz**thrd2 - omz**thrd2)

      darg1dzeta = dw1dzeta*(ONE-GAt2) + gta22*w1*(dgat2da*dAdzeta + dgat2dt2*dt2dzeta)
      dH1dzeta = 3d0*H1*dphidzeta/phi + GAMMA*phi3*darg1dzeta/arg1
      dEc1dzeta = declda1dz + dH1dzeta

      ddxdzeta = thrd2*(opz**THRD - omz**THRD)
      dgcdzeta = -dxc*ddxdzeta*(1.0_q - zeta**12) - 12.0_q*zeta**11*(1.0_q - dxc*(dx - 1.0_q))
      dEc0dzeta = dgcdzeta*(eclda0 + H0)

      dalphadzeta = -alpha/(tueg+eta*tvw)*tueg*ddsdzeta/ds
      dfcadzeta = dfcada*dalphadzeta

      vcpol = dH1dzeta + dfcadzeta*(Ec0-Ec1) + fca*(dEc0dzeta-dEc1dzeta)

      dvcup = vcn + omz*vcpol
      dvcdw = vcn - opz*vcpol

      vcup = eclda1 + declda1drs*drsdn + omz*declda1dz
      vcdw = eclda1 + declda1drs*drsdn - opz*declda1dz

      return
      end SUBROUTINE CRTML
!
!
      SUBROUTINE ts_pcrep(rho, p, ql, dfsdp, dfsdq, fs)
GPACC !$ACC ROUTINE SEQ
!     Perdew & Constantin kinetic energy functional enhancement factor
!     Written by:
!     Daniel Mejia-Rodriguez
!     QTP, Department of Physics, University of Florida
!              
!     Modified by H.Francisco for VASP and PC_rep:
!     QTP, Department of Physics, University of Florida
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214103(2023)    

      USE prec
      USE constant

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(ts_pcrep)

!
      REAL(q) :: fs !< [Output]
      REAL(q) :: dfsdp !< [Output]
      REAL(q) :: dfsdq !< [Output]
!
      REAL(q) :: rho !< [Input]
!
      REAL(q) :: rho83, rho53
      REAL(q) :: p, ql
      REAL(q) :: rhoval, Pconst
!
      REAL(q) :: ge2,ge4,mge4,zpc,fpc,bfpc,lapval
      REAL(q) :: exp8,exp9,factor,ftw,mge4den,opftw

!     functional derivatives below FFFFFFFFFFFF

      REAL(q) :: dge2dp,dge2dq,dge4dp,dge4dq,dmge4dp,dmge4dq
      REAL(q) :: dfactordp,dfactordq
      REAL(q) :: dfpc,dfsdz

!     functional derivatives above FFFFFFFFFFFF
      REAL(q),PARAMETER :: F83=8.0_q/3.0_q
      REAL(q),PARAMETER :: F23=2.0_q/3.0_q
      REAL(q),PARAMETER :: F18=1.0_q/8.0_q
      REAL(q),PARAMETER :: F53=5.0_q/3.0_q
      REAL(q),PARAMETER :: amgga=1.50440153364382_q !v2_sregTML
      REAL(q),PARAMETER :: bmgga=0.615645486593429_q !v2_sregTML
      REAL(q),PARAMETER :: F527=5.0_q/27.0_q
      REAL(q),PARAMETER :: F209=20.0_q/9.0_q
      REAL(q),PARAMETER :: F881=8.0_q/81.0_q
      REAL(q),PARAMETER :: F19=1.0_q/9.0_q
      REAL(q),PARAMETER :: F8243=8.0_q/243.0_q
      REAL(q),PARAMETER :: F16243=16.0_q/243.0_q
      REAL(q),PARAMETER :: F1681=16.0_q/81.0_q
      REAL(q),PARAMETER :: ONE=1.0_q
      REAL(q),PARAMETER :: TOL2=1.0E5_q
!
      Pconst=(3.0_q*PI*PI)**F23

            rhoval=rho
            rho53 = rhoval**F53
            rho83 = rho53*rhoval

            ge2 = F527*p + F209*ql

            dge2dp = F527
            dge2dq = F209

            ge4 = ql*(F881*ql - F19*p) + F8243*p*p
            dge4dp = F16243*p - F19*ql
            dge4dq = F1681*ql - F19*p

            ftw = F53*p
            opftw = ONE + ftw
            factor = ge4/opftw
            dfactordp = ONE/opftw * (dge4dp - F53*factor)
            dfactordq = dge4dq/opftw

            if (factor.gt.TOL2) then
              mge4 = (ONE/ge4 + ge2/ge4 + ONE)*opftw
              dmge4dp = (dge2dp + dge4dp - mge4*dfactordp)/factor
              dmge4dq = (dge2dq + dge4dq - mge4*dfactordq)/factor
            else
              mge4den = dsqrt(ONE + factor**2)
              mge4 = (ONE + ge2 + ge4)/mge4den
              dmge4dp = (dge2dp + dge4dp - &
     &                   mge4*factor*dfactordp/mge4den)/mge4den
              dmge4dq = (dge2dq + dge4dq - &
     &                   mge4*factor*dfactordq/mge4den)/mge4den
            endif
!
            zpc = mge4 - ftw
            fpc = 0.0_q
            dfpc = 0.0_q
!
            if (zpc.le.0.001_q) then
              fpc = 0.0_q
              dfpc = 0.0_q
            elseif (zpc.lt.amgga) then
              exp8 = dexp(-amgga/zpc)
              exp9 = dexp(-amgga/(amgga-zpc))
              bfpc = (exp8*exp9 + exp8)/(exp8 + exp9)
              fpc = bfpc**bmgga
              dfpc = amgga*bmgga*fpc/(exp8 + exp9)
              if (amgga-zpc.gt.0.0013_q) then
                dfpc = dfpc*exp9*(ONE/zpc**2.0_q + &
     &                 (ONE-exp8)/(ONE+exp9*(amgga-zpc)**2.0_q))
              else
                dfpc = dfpc*exp9/zpc**2.0_q
              end if
            else
              fpc=ONE
              dfpc=0.0_q
            endif

            fs=ftw+zpc*fpc

            dfsdz = fpc + zpc*dfpc
            dfsdp = F53 + dfsdz*(dmge4dp - F53)
            dfsdq = dfsdz*dmge4dq

      return
      end SUBROUTINE ts_pcrep

      subroutine pcrepmgga(d,dthrd,s,qq,tau,taud,taudd,taudl)
GPACC !$ACC ROUTINE SEQ
!     Perdew & Constantin kinetic energy functional enhancement factor
!     Written by:
!     Daniel Mejia-Rodriguez
!     QTP, Department of Physics, University of Florida
!              
!     Modified by H.Francisco for PC_rep:
!     QTP, Department of Physics, University of Florida              
!     References:              
!     H. Francisco, A.C. Cancio and S.B. Trickey, 
!     J. Chem. Phys. 159, 214103(2023)    

      use prec
      use constant
      implicit none
GPOMP  !$OMP DECLARE TARGET(pcrepmgga)
      real(q) d,dthrd,s,qq,tau,taud,taudd,taudl
      real(q) ge2,ge4,mge4,ftw,oneftw,factor
      real(q) p,exp1,exp2,fz,fzb,z,ft,tueg
      real(q) dge2p,dge2q,dge4p,dge4q,dftwp,dfactorp,dfactorq
      real(q) dmge4p,dmge4q,dzp,dzq,dfzb,dftdz,dftdp,dftdq,tunif,tlda
      real(q),parameter :: ckf=(3.0_q*pi*pi)**(1.0_q/3.0_q)
      REAL(q),PARAMETER :: A=1.50440153364382_q !v2-sregTML
      REAL(q),PARAMETER :: B=0.615645486593429_q !v2-sregTML

      REAL(q),PARAMETER :: THRESH1 = 0.001_q
      REAL(q),PARAMETER :: THRESH2 = A-0.0013_q
      REAL(q),PARAMETER :: THRESH3 = 0.01_q
      REAL(q),PARAMETER :: THRESH4 = 0.997_q
      REAL(q),PARAMETER :: THRESH5 = 1.003_q

      tunif = 0.3_q*(ckf*dthrd)**2
      tlda = tunif*d

      p = s*s

      GE2 = 5.0_q/27.0_q*P + 20.0_q/9.0_q*QQ
      GE4 = 8.0_q/81.0_q*QQ*QQ - P*QQ/9.0_q + 8.0_q/243.0_q*P*P
      FTW = 5.0_q/3.0_q*P
      ONEFTW = 1.0_q + FTW
      FACTOR = GE4/ONEFTW

      IF (FACTOR.GT.1.0E8_q) THEN
        MGE4 = ONEFTW/GE4 + ONEFTW*GE2/GE4 + ONEFTW
      ELSE
        MGE4 = (1.0_q + GE2 + GE4)/DSQRT(1.0_q + FACTOR*FACTOR)
      END IF

      Z = MGE4 - FTW
      IF ((Z.GT.THRESH1).AND.(Z.LT.THRESH2)) THEN
        EXP1 = DEXP(-A/Z)
        EXP2 = DEXP(-A/(A-Z))
        FZ = EXP1*(1.0_q + EXP2)/(EXP1 + EXP2)
        FZB = FZ**B
      ELSE IF (Z.GE.THRESH2) THEN
        FZB = 1.0_q
      ELSE
        FZB = 0.0_q
      ENDIF

      FT = Z*FZB + FTW
      tau = ft*tlda
      DGE2P = 5.0_q/27.0_q
      DGE2Q = 20.0_q/9.0_q

      DGE4P = -QQ/9.0_q + 16.0_q/243.0_q*P
      DGE4Q = 16.0_q/81.0_q*QQ - P/9.0_q

      DFTWP = 5.0_q/3.0_q
      DFACTORP = (DGE4P - DFTWP*FACTOR)/ONEFTW
      DFACTORQ = DGE4Q/ONEFTW

      IF (FACTOR.GT.1.0E8_q) THEN
        DMGE4P = DFTWP*(1.0_q + (1.0_q + GE2)/GE4) + ONEFTW/GE4*(DGE2P - DGE4P/GE4*(1.0_q + GE2))
        DMGE4Q = ONEFTW/GE4*(DGE2Q - DGE4Q/GE4*(1.0_q + GE2))
      ELSE
        DMGE4P = (DGE2P + DGE4P)/DSQRT(1.0_q + FACTOR*FACTOR) - MGE4*FACTOR*DFACTORP/(1.0_q + FACTOR*FACTOR)
        DMGE4Q = (DGE2Q + DGE4Q)/DSQRT(1.0_q + FACTOR*FACTOR) - MGE4*FACTOR*DFACTORQ/(1.0_q + FACTOR*FACTOR)
      END IF

      DZP = DMGE4P - DFTWP
      DZQ = DMGE4Q

      IF ((Z.GT.THRESH1).AND.(Z.LT.THRESH2)) THEN
        DFZB = A*B*FZB*EXP2/(Z**2*(EXP1+EXP2))*(1.0_q+Z**2*(1.0_q-EXP1)/((1.0_q+EXP2)*(A-Z)**2))
      ELSE IF (Z.GE.THRESH2) THEN
        DFZB = 0.0_q
      ELSE
        DFZB = 0.0_q
      ENDIF
      dftdz = fzb + z*dfzb
      dftdp = dftdz*dzp + dftwp
      dftdq = dftdz*dzq

      taud = 5.0_q/3.0_q*tunif*ft - tunif*(8.0_q/3.0_q*dftdp*p + 5.0_q/3.0_q*dftdq*qq)
      taudd = dftdp*s*0.3_q*ckf*dthrd
      taudl = dftdq*0.075_q

      return
      end subroutine pcrepmgga


!************************ SUBROUTINE VTASKx ****************************
!
! calculates the first order derivatives of Ex wrt n and |grad(n)|
! Aschebrock et. al. PRR (2019)
!
! Written by Timo Lebeda, April 2022
!
! everything in Hartree units
!
! ATTENTION: Every values are passed "as they are", i.e. including
! possibly unphysical numerical errors (e.g. negative charge densities)
! values need to be checked accordingly
!
! RX                           density up,down
! DRX                          abs. val. gradient of density up/down
! TAUX                         kinetic energy density up/down
! VXD                          THE DERIVATIVES OF EX WRT n
! VXDD                         THE DERIVATIVES OF EX WRT |grad n|
! AMUXD                        THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VTASKx(RX,DRX,TAUX,Ex,VXD,VXDD,AMUXD)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(VTASKx)

      REAL(q) :: RX,DRX,TAUX,EX,VXD,VXDD,AMUXD
      REAL(q) :: ZERO,ONE,TWO,THREE,FOUR
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: TAU_RHO,TAUW_RHO,TAU0,TAUW_TEMP
      REAL(q) :: RHO,DRHO,P
      REAL(q) :: ALPHA,DALPHAD,DALPHADD,DALPHADTAU
      REAL(q) :: DPD,DPDD,EXUNIF,EXLDA,EXDLDA
      REAL(q) :: FX,DFXDP,DFXDA,DFX_D,DFX_DD,DFX_DTAU

      REAL(q), PARAMETER :: HX0=1.174_q
      REAL(q), PARAMETER :: D =10.0_q
! other parameters
      PARAMETER (ZERO=0.0_q)
      PARAMETER (ONE   =1._q)
      PARAMETER (TWO   =2._q)
      PARAMETER (THREE =3._q)
      PARAMETER (FOUR  =4._q)
      PARAMETER (THRD =1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      REAL(q), parameter :: GAMMA  =  3.0936677262801359310_q ! (3*PI^2)^(1/3)
      REAL(q), parameter :: GAMMA2 =  GAMMA*GAMMA
      REAL(q), parameter :: AX=-0.738558766382022405884230032680836_q
GPACC !$ACC ROUTINE(TASKfx) SEQ

      VXD=0._q;
      VXDD=0._q;
      AMUXD=0._q;

! spin up
      RHO=TWO*RX
      DRHO=TWO*DRX
      TAUW_RHO=TWO*0.125_q*DRX**2._q/RX
      TAU_RHO=TWO*TAUX
      TAU0=0.3_q*GAMMA2*(RHO)**THRD5

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4

! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|

      P=(DRHO)**TWO/(FOUR*GAMMA2*(RHO)**THRD8)  ! P==s**2
      DPD=-THRD8*P  ! /RHO  commented out to save division by RHO

      ALPHA = MAX((TAU_RHO - TAUW_RHO)/TAU0,ZERO)
      DALPHAD = THRD5 * (P-ALPHA) !/RHO  commented out to save division by RHO

! calculate the exchange enhancement factor and its derivatives wrt p(=s²) and alpha
       call TASKfx(D, HX0, P, ALPHA, FX, DFXDP, DFXDA)

! Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=DFXDP*DPD+DFXDA*DALPHAD

!   OUTPUT THE VXD,VXDD AND AMUXD
      VXD=EXDLDA*FX+EXUNIF*DFX_D  ! actually EXLDA instead of EXUNIF, but saving division by RHO
      VXDD=AX*DRHO/(TWO*GAMMA2*(RHO)**THRD4) * (DFXDP - THRD5 * DFXDA)
      AMUXD=AX*10.0_q/(THREE*GAMMA2*RHO**THRD)*DFXDA      ! directly cancel rho to improve numerics
      EX=EX+EXLDA*FX/TWO

      RETURN
    END SUBROUTINE VTASKx

! TL 04.2022, Enhancement factor of TASK exchange and its derivatives
      SUBROUTINE TASKfx(D, HX0, p, alpha, F, dp_F, dalpha_F)
GPACC !$ACC ROUTINE SEQ
         ! --------------------------------------------
         ! This subroutine calculates the
         ! enhancement factor Fx(s,alpha)
         ! for exchange part of the TASKx meta-GGA.
         ! --------------------------------------------

      USE prec
      USE constant
      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(TASKfx)

      INTEGER, PARAMETER :: Na=2
      INTEGER, PARAMETER :: Nb=4

      REAL(q) :: D, HX0
      REAL(q) :: p, alpha, F, dp_F, dalpha_F
      REAL(q) :: gx, dp_gx, hx1, dp_hx1, fx, dalpha_fx
      REAL(q) :: a(0:Na), b(0:Nb)

      REAL(q), PARAMETER :: TWO=2.0_q
      REAL(q), PARAMETER :: THREE=3.0_q
      REAL(q), PARAMETER :: FOUR=4.0_q
GPACC !$ACC ROUTINE(calc_gx) SEQ
GPACC !$ACC ROUTINE(calc_dgx) SEQ
GPACC !$ACC ROUTINE(ChebEv) SEQ
GPACC !$ACC ROUTINE(ChebDerEv) SEQ

      a(0:Na) = (/0.938719_q, -0.076371_q, -0.0150899_q/)
      b(0:Nb) = (/-0.628591_q, -2.10315_q, -0.5_q, 0.103153_q, 0.128591_q/)
      call calc_gx(p,gx)   ! p == s**2
      call ChebEv(a,Na,p,hx1)
      call ChebEv(b,Nb,alpha,fx)

      ! exchange enhancement factor
      F = HX0*gx+(1.0_q-fx)*(hx1-HX0)*gx**D

      call calc_dgx(p,dp_gx) ! dp_gx = dgx/d(s^2) = 1/(2s)*dgx/ds
      call ChebDerEv(a,Na,p,dp_hx1) ! dp_hx1 = dhx1/d(s^2) = 1/(2s)*dhx1/ds
      call ChebDerEv(b,Nb,alpha,dalpha_fx)

      dalpha_F = -dalpha_fx*(hx1-HX0)*gx**D
      dp_F      = (HX0+(1.0_q-fx)*(hx1-HX0)*D*gx**(D-1.0_q))*dp_gx+(1.0_q-fx)*dp_hx1*gx**D

      RETURN
    END SUBROUTINE TASKfx

      SUBROUTINE ChebEv (coef, n, r, value)
GPACC !$ACC ROUTINE SEQ
         ! --------------------------------------------
         ! This subroutine evaluates
         !
         ! S(r) = SUM_{i=0,N) coef_i*T_i(x(r))
         !
         ! with x(r)=(r-1)/(r+1) for given r in [0,Infinity[
         ! and given coefficients(0:N) and N by Clenshaw's recurrence formula,
         ! where T_i(x) are the Chebyshev polynomials of the first kind.
         !
         ! For further information see:
         ! Numerical Recipes in Fortran 90, Second Edition (1996), Chapter 5.8
         ! --------------------------------------------

         USE prec
         USE constant

         IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(ChebEv)

         integer :: n            ! N
         real(q) :: coef(0:n)    ! Coefficient Array
         real(q) :: r            ! Value of r
         real(q) :: value        ! S(r)

         real(q) :: d,dd,sv,x,x2
         REAL(q) :: ZERO,ONE,TWO

         integer :: j

         PARAMETER (ZERO=0.0_q)
         PARAMETER (ONE=1._q)
         PARAMETER (TWO=2._q)

         d  = ZERO
         dd = ZERO
         x  = (r-ONE)/(r+ONE)
         x2 = TWO*x

         do j=n,1,-1
            sv = d
            d  = x2*d-dd+coef(j)
            dd = sv
         end do

         value = x*d-dd+coef(0)

         RETURN
      END SUBROUTINE ChebEv


      SUBROUTINE ChDer (coef, coefder, n)
GPACC !$ACC ROUTINE SEQ
         ! --------------------------------------------
         ! This subroutine calculates the coefficients of the derivative
         !
         ! dS/dx = SUM_{i=0,N) coefder_i*T_i(x)
         !
         ! of
         !
         ! S(x) = SUM_{i=0,N) coef_i*T_i(x),
         !
         ! where T_i(x) are the Chebyshev polynomials of the first kind. This
         ! subroutine is used in ChebDerEv.
         !
         ! For further information see:
         ! Numerical Recipes in Fortran 90, Second Edition (1996), Chapter 5.8
         ! --------------------------------------------

         USE prec
         USE constant

         IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(ChDer)

         integer :: n            ! N
         real(q) :: coef(0:n)    ! Coefficient Array
         real(q) :: coefder(0:n) ! Coefficient Array of Derivative

         REAL(q) :: ZERO,TWO

         integer :: j

         PARAMETER (ZERO=0.0_q)
         PARAMETER (TWO=2._q)

         coefder(n)   = ZERO
         coefder(n-1) = TWO*n*coef(n)

         if(n.ge.2) then
            do j=n-2,0,-1
               coefder(j) = coefder(j+2) + TWO*(j+1)*coef(j+1)
            end do
         end if

         coefder(0) = coefder(0)/TWO

         RETURN
      END SUBROUTINE ChDer


      SUBROUTINE ChebDerEv (coef, n, r, value)
GPACC !$ACC ROUTINE SEQ
         ! --------------------------------------------
         ! This subroutine evaluates the derivative S'(r) of
         !
         ! S(r) = SUM_{i=0,N) coef_i*T_i(x(r))
         !
         ! with x(r)=(r-1)/(r+1) for given r in [0,Infinity[
         ! and given coefficients(0:N) and N by Clenshaw's recurrence formula,
         ! where T_i(x) are the Chebyshev polynomials of the first kind.
         !
         ! For further information see:
         ! Numerical Recipes in Fortran 90, Second Edition (1996), Chapter 5.8
         ! --------------------------------------------

         USE prec
         USE constant

         IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(ChebDerEv)

         integer :: n            ! N
         real(q) :: coef(0:n)    ! Coefficient Array
         real(q) :: r            ! Value of r
         real(q) :: value        ! S'(r)

         real(q) :: coefder(0:n)
         REAL(q) :: ONE,TWO

         PARAMETER (ONE=1._q)
         PARAMETER (TWO=2._q)
GPACC !$ACC ROUTINE(chder) SEQ
GPACC !$ACC ROUTINE(chebev) SEQ

         call chder(coef,coefder,n)
         call chebev(coefder,n,r,value)

         value = value*TWO/((ONE+r)*(ONE+r))

      END SUBROUTINE ChebDerEv


!************************ SUBROUTINE VLAKx ****************************
!
! calculates the first order derivatives of Ex wrt n and |grad(n)|, and tau
! Lebeda et. al. PRL (2024)
!
! Written by Timo Lebeda, August 2024
!
! everything in Hartree units
!
! RX                           density up,down
! DRX                          abs. val. gradient of density up/down
! TAUX                         kinetic energy density up/down
! VXD                          THE DERIVATIVES OF EX WRT n
! VXDD                         THE DERIVATIVES OF EX WRT |grad n|
! AMUXD                        THE DERIVATIVES OF EX WRT TAU
!
!***********************************************************************

      SUBROUTINE VLAKx(RX,DRX,TAUX,Ex,VXD,VXDD,AMUXD)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(VLAKx)

      REAL(q) :: RX,DRX,TAUX,EX,VXD,VXDD,AMUXD
      REAL(q) :: ZERO,TWO,THREE,FOUR
      REAL(q) :: THRD,THRD2,THRD4,THRD5,THRD8
      REAL(q) :: TAU_RHO,TAUW_RHO,TAU0,TAUW_TEMP
      REAL(q) :: RHO,DRHO,P
      REAL(q) :: ALPHA,DALPHAD,DALPHADD,DALPHADTAU
      REAL(q) :: DPD,DPDD,EXUNIF,EXLDA,EXDLDA
      REAL(q) :: FX,DFXDP,DFXDA,DFX_D,DFX_DD,DFX_DTAU

      REAL(q), PARAMETER :: HX0=1.174_q
! other parameters
      PARAMETER (ZERO=0.0_q)
      PARAMETER (TWO   =2._q)
      PARAMETER (THREE =3._q)
      PARAMETER (FOUR  =4._q)
      PARAMETER (THRD =1._q/3._q)
      PARAMETER (THRD2=2._q*THRD)
      PARAMETER (THRD4=4._q*THRD)
      PARAMETER (THRD5=1._q+THRD2)
      PARAMETER (THRD8=1._q+THRD5)
      REAL(q), parameter :: GAMMA  =  3.0936677262801359310_q ! (3*PI^2)^(1/3)
      REAL(q), parameter :: GAMMA2 =  GAMMA*GAMMA
      REAL(q), parameter :: AX=-0.738558766382022405884230032680836_q
GPACC !$ACC ROUTINE(LAKfx) SEQ

      VXD=0._q;
      VXDD=0._q;
      AMUXD=0._q;
      
      TAUW_TEMP=0.125_q*DRX**2._q/RX

! spin up
      RHO=TWO*RX
      DRHO=TWO*DRX
      TAUW_RHO=TWO*TAUW_TEMP
      TAU_RHO=TWO*TAUX
      TAU0=0.3_q*GAMMA2*(RHO)**THRD5

!----------------------------------------------------------------------
! construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4

! CONSTRUCT FX AND ITS DERIVATIVES WRT n AND |grad n|, and tau

      P=(DRHO)**TWO/(FOUR*GAMMA2*(RHO)**THRD8)  ! P==S**2
      DPD=-THRD8*P  ! /RHO  commented out to save division by RHO

      ALPHA = (TAU_RHO - TAUW_RHO)/TAU0
      DALPHAD = THRD5 * (P-ALPHA) !/RHO  commented out to save division by RHO

! calculate the exchange enhancement factor and its derivatives wrt p and alpha
      call LAKfx(HX0, P, ALPHA, FX, DFXDP, DFXDA)

! Calculate Fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=DFXDP*DPD+DFXDA*DALPHAD

!   OUTPUT THE VXD,VXDD AND AMUXD
      VXD=EXDLDA*FX+EXUNIF*DFX_D    ! actually EXLDA instead of EXUNIF in 2nd term, but saving division by RHO
      VXDD=AX*DRHO/(TWO*GAMMA2*(RHO)**THRD4) * (DFXDP - THRD5 * DFXDA)
      AMUXD=AX*10.0_q/(THREE*GAMMA2*RHO**THRD)*DFXDA      ! directly cancel rho to improve numerics
      EX=EX+EXLDA*FX/TWO

      RETURN
    END SUBROUTINE VLAKx

! TL 08.2024, Enhancement factor of LAK exchange and its derivatives
      SUBROUTINE LAKfx(HX0, p, alpha, F, dp_F, dalpha_F)
GPACC !$ACC ROUTINE SEQ
        ! --------------------------------------------
         ! This subroutine calculates the
         ! enhancement factor Fx(s,alpha)
         ! for the exchange part of the LAK meta-GGA.
         ! --------------------------------------------

      USE prec
      USE constant

      IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(LAKfx)

      REAL(q) :: HX0
      REAL(q) :: p, alpha, F, dp_F, dalpha_F
      REAL(q) :: gx, dp_gx, gnum, dp_gnum, hx1, dp_hx1, fx, dalpha_fx
GPACC !$ACC ROUTINE(calc_gx) SEQ
GPACC !$ACC ROUTINE(calc_gnum) SEQ
GPACC !$ACC ROUTINE(calc_hx1GE4) SEQ
GPACC !$ACC ROUTINE(calc_fxGE4) SEQ
GPACC !$ACC ROUTINE(calc_dgx) SEQ
GPACC !$ACC ROUTINE(calc_dgnum) SEQ
GPACC !$ACC ROUTINE(calc_dhx1GE4) SEQ
GPACC !$ACC ROUTINE(calc_dfxGE4) SEQ

      call calc_gx(p,gx)   ! p == s**2
      call calc_gnum(p,gnum)
      call calc_hx1GE4(HX0,p,hx1)
      call calc_fxGE4(HX0,alpha,fx)

      ! exchange enhancement factor
      F = HX0*gx + (1.0_q-fx)*(hx1-HX0)*gnum

      call calc_dgx(p,dp_gx) ! dp_gx = dgx/d(s^2) = 1/(2s)*dgx/ds
      call calc_dgnum(p,dp_gnum) ! dp_g = dg/d(s^2) = 1/(2s)*dg/ds
      call calc_dhx1GE4(HX0,p,dp_hx1) ! dp_hx1 = dhx1/d(s^2) = 1/(2s)*dhx1/ds
      call calc_dfxGE4(HX0,alpha,dalpha_fx)

      dalpha_F = -dalpha_fx*(hx1-HX0)*gnum
      dp_F      = HX0 * dp_gx + (1.0_q-fx)*(hx1-HX0) * dp_gnum + (1.0_q-fx)*gnum * dp_hx1 

      RETURN
    END SUBROUTINE LAKfx


    SUBROUTINE calc_gx(x,gx)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_gx)

       REAL(q) :: x, gx
       REAL(q) :: a

       PARAMETER ( a  = 4.9479_q)

       if(x.gt.0.0_q) then
          gx  = 1.0_q-exp(-a*x**(-1.0_q/4.0_q))
       else
          gx  = 1.0_q
       endif

       RETURN
    END SUBROUTINE calc_gx

    SUBROUTINE calc_dgx(x,dgx)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_dgx)

       REAL(q) :: x, dgx
       REAL(q) :: a

       PARAMETER ( a  = 4.9479_q)

       if(x.gt.0.0_q) then
          dgx = -a/4.0_q*x**(-5.0_q/4.0_q)*exp(-a*x**(-1.0_q/4.0_q))
       else
          dgx = 0.0_q
       endif

       RETURN
    END SUBROUTINE calc_dgx

    SUBROUTINE calc_gnum(x,gnum)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_gnum)     

       REAL(q) :: x, gnum
       REAL(q) :: a

       PARAMETER ( a  = 5.0_q)

       if(x.gt.0.0_q) then
          gnum  = 1.0_q-exp(-a*a/x)
       else
          gnum  = 1.0_q
       endif

       RETURN
    END SUBROUTINE calc_gnum

    SUBROUTINE calc_dgnum(x,dgnum)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_dgnum)

       REAL(q) :: x, dgnum
       REAL(q) :: a

       PARAMETER ( a  = 5.0_q)

       if(x.gt.0.0_q) then
          dgnum = -a*a/x**2.0_q*exp(-a*a/x)
       else
          dgnum = 0.0_q
       endif

       RETURN
    END SUBROUTINE calc_dgnum

    SUBROUTINE calc_GE4(x,GE4)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_GE4)

       REAL(q) :: x, GE4
       REAL(q) :: mus, Cs

       PARAMETER ( mus  = -0.0320227_q)
       PARAMETER ( Cs   = -0.0886966_q)

       GE4  = 1.0_q + mus*x + Cs*x*x

       RETURN
    END SUBROUTINE calc_GE4

    SUBROUTINE calc_dGE4(x,dGE4)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_dGE4)

       REAL(q) :: x, dGE4
       REAL(q) :: mus, Cs

       PARAMETER ( mus  = -0.0320227_q)
       PARAMETER ( Cs   = -0.0886966_q)

       dGE4  = mus + 2.0_q*Cs*x

       RETURN
    END SUBROUTINE calc_dGE4

    SUBROUTINE calc_hx1GE4(HX0,x,hx1)
      GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_hx1GE4)

       REAL(q) :: HX0, x, hx1
       REAL(q) :: GE4, gx
       REAL(q) :: a

       PARAMETER ( a  = 1.1_q)
GPACC !$ACC ROUTINE(calc_GE4) SEQ
GPACC !$ACC ROUTINE(calc_gx) SEQ

       call calc_GE4(x,GE4)
       call calc_gx(x,gx)
       if(x.gt.0.0_q) then
          hx1  = a + (1.0_q - exp(- a/x*a/(1.0_q+x)) ) * (GE4 - HX0*(gx-1.0_q) - a)
       else
          hx1  = 1.0_q
       endif

       RETURN
    END SUBROUTINE calc_hx1GE4

    SUBROUTINE calc_dhx1GE4(HX0,x,dhx1)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_dhx1GE4)
       REAL(q) :: HX0, x, dhx1
       REAL(q) :: GE4, dGE4, gx, dgx
       REAL(q) :: a

       PARAMETER ( a  = 1.1_q)
GPACC !$ACC ROUTINE(calc_GE4) SEQ
GPACC !$ACC ROUTINE(calc_dGE4) SEQ
GPACC !$ACC ROUTINE(calc_gx) SEQ
GPACC !$ACC ROUTINE(calc_dgx) SEQ

       call calc_GE4(x,GE4)
       call calc_dGE4(x,dGE4)
       call calc_gx(x,gx)
       call calc_dgx(x,dgx)
       if(x.gt.0.0_q) then
          dhx1 = (1.0_q - exp(- a/x*a/(1.0_q+x)) ) * (dGE4 - HX0*dgx) - a**2.0_q*(2.0_q*x+1.0_q)/(x*(1.0_q+x))**2.0_q * exp(- a/x*a/(1.0_q+x)) * (GE4 - HX0*(gx-1.0_q) - a)
       else
          dhx1 = 0.0_q
       endif

       RETURN
    END SUBROUTINE calc_dhx1GE4

    SUBROUTINE calc_fxGE4(HX0,x,fx)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_fxGE4)

       REAL(q) :: HX0, x, fx
       REAL(q) :: GE2xa, GE4xa, alpha_tilde, muax, Ca

       PARAMETER ( muax = -0.20989736912384607_q)
       PARAMETER ( Ca = (73.0_q - 50.0_q*muax) / 5000.0_q )

       GE2xa = muax / (HX0-1.0_q)
       GE4xa = (Ca+muax) / (HX0-1.0_q)

       if(x.gt.0._q) then
          alpha_tilde        = (x-1.0_q)/x
          fx  = 2.0_q/PI*atan(PI/2.0_q*(GE2xa*alpha_tilde + GE4xa*(x-1.0_q)**2.0_q) )
       elseif(x.lt.0._q) then
          alpha_tilde        = (x-1.0_q)/x
          fx  = 2.0_q + 2.0_q/PI*atan(PI/2.0_q*(GE2xa*alpha_tilde + GE4xa*(x-1.0_q)**2.0_q) )     ! 2+fx is the continuous extension of fx(alpha) to alpha<0
       else
          fx  = 1.0_q
       endif

       RETURN
    END SUBROUTINE calc_fxGE4

    SUBROUTINE calc_dfxGE4(HX0,x,dfx)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(calc_dfxGE4)

       REAL(q) :: HX0, x, dfx
       REAL(q) :: GE2xa, GE4xa, alpha_tilde, dalpha_alpha_tilde, muax, ca

       PARAMETER ( muax = -0.20989736912384607_q)
       PARAMETER ( Ca = (73.0_q - 50.0_q*muax) / 5000.0_q )

       GE2xa = muax / (HX0-1.0_q)
       GE4xa = (Ca+muax) / (HX0-1.0_q)

       if(x.gt.0._q) then
          alpha_tilde        = (x-1.0_q)/x
          dalpha_alpha_tilde = 1.0_q/x**2.0_q
          dfx  = (GE2xa*dalpha_alpha_tilde + GE4xa*2.0_q*(x-1.0_q) ) / (1.0_q + (PI/2.0_q * (GE2xa*alpha_tilde + GE4xa*(x-1.0_q)**2.0_q) )**2.0_q )
       elseif(x.lt.0._q) then
          alpha_tilde        = (x-1.0_q)/x
          dalpha_alpha_tilde = 1.0_q/x**2.0_q
          dfx  = (GE2xa*dalpha_alpha_tilde + GE4xa*2.0_q*(x-1.0_q) ) / (1.0_q + (PI/2.0_q * (GE2xa*alpha_tilde + GE4xa*(x-1.0_q)**2.0_q) )**2.0_q ) ! The derivative of fx is continuous at alpha=0
       else
          dfx  = GE2xa/(PI/2.0_q*GE2xa)**2.0_q  ! ArcTan has a finite derivative at alpha=0
       endif

       RETURN
    END SUBROUTINE calc_dfxGE4


!************************ SUBROUTINE VLAKc ****************************
!
! calculates the first order derivatives of Ec wrt n and |grad(n)|, and tau
! Lebeda et. al. PRL (2024)
!
! Written by Timo Lebeda, August 2024
!
! everything in Hartree units
!
! RU,RD                        density up,down
! DRU, DRD                     abs. val. gradient of density up/down
! DRT                          abs. val. gradient of total density
! TAUU,TAUD                    kinetic energy density up/down
! TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
! VCD1 VCD2                    THE DERIVATIVES OF EC WRT n
! VCDD                         THE DERIVATIVES OF EC WRT |grad n|
! AMUCD1, AMUCD2               THE DERIVATIVES OF EC WRT TAU
!
!***********************************************************************

    SUBROUTINE VLAKc(RU,RD,DRT,TAUU,TAUD,Ec,VCD1,VCD2,VCDD,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(VLAKc)

       ! In-/Output
       REAL(q) :: Ec,DRT,VCDD
       REAL(q) :: RU,TAUU,VCD1,AMUCD1
       REAL(q) :: RD,TAUD,VCD2,AMUCD2
       ! Numbers
       REAL(q) :: ZERO,ONE,MONE,TWO,FOUR
       REAL(q) :: HALF,THRD,THRD2,THRD4,THRD5,THRD8
       ! Composite local quantities
      
       REAL(q) :: rs, zeta, s2, alpha
       REAL(q) :: ds, dzeta_ds, tauF, tauW, tau
       REAL(q) :: drs_ec, dzeta_ec, ds2_ec, dalpha_ec

! other parameters
       PARAMETER (ZERO=0.0_q)
       PARAMETER (ONE   =1._q)
       PARAMETER (MONE  =-1._q)
       PARAMETER (TWO   =2._q)
       PARAMETER (FOUR  =4._q)
       PARAMETER (HALF =0.5_q)
       PARAMETER (THRD =1._q/3._q)
       PARAMETER (THRD2=2._q*THRD)
       PARAMETER (THRD4=4._q*THRD)
       PARAMETER (THRD5=1._q+THRD2)
       PARAMETER (THRD8=1._q+THRD5)
       REAL(q), parameter :: GAMMA  =  3.0936677262801359310_q ! (3*PI^2)^(1/3)
       REAL(q), parameter :: GAMMA2 =  GAMMA*GAMMA
GPACC !$ACC ROUTINE(LAK_corr) SEQ

       VCD1=0._q; VCD2=0._q
       VCDD=0._q
       AMUCD1=0._q; AMUCD2=0._q;
      
       rs   = (THRD4*PI*(RU+RD))**(-THRD)
       zeta = (RU-RD)/(RU+RD)
       s2   = (DRT**TWO)/(4.0_q*GAMMA2*(RU+RD)**THRD8)

       ds       = HALF*((ONE+zeta)**THRD5+(ONE-zeta)**THRD5)
       dzeta_ds = HALF*THRD5*((ONE+zeta)**THRD2-(ONE-zeta)**THRD2)

       tau   = TAUU+TAUD
       tauF  = 3.0_q/10.0_q*GAMMA2*(RU+RD)**THRD5*ds
       tauW  = (DRT**TWO)/(8.0_q*(RU+RD))
       alpha = (tau-tauW)/tauF

! calculate the correlation enhancement factor and its derivatives wrt rho, zeta, p=s^2, and alpha
       call LAK_corr(rs, zeta, s2, alpha, ec, drs_ec, dzeta_ec, ds2_ec, dalpha_ec)

! Calculate ec DERIVATIVES WRT n, |grad n|, and tau
       ! dnUP_ec
       VCD1 = ec - THRD*rs*drs_ec + ( ONE-zeta)*dzeta_ec - THRD8*s2*ds2_ec + (-THRD5*alpha+tauW/tauF+alpha*(zeta+MONE)*dzeta_ds/ds)*dalpha_ec
       ! dnDOWN_ec
       VCD2 = ec - THRD*rs*drs_ec + (MONE-zeta)*dzeta_ec - THRD8*s2*ds2_ec + (-THRD5*alpha+tauW/tauF+alpha*(zeta+ ONE)*dzeta_ds/ds)*dalpha_ec
       
       ! d|n|_ec
       VCDD = TWO*(ds2_ec-5.0_q/3.0_q/ds*dalpha_ec)/(FOUR*GAMMA2*(RU+RD)**THRD5)*DRT
       
       ! dtau_ec
       AMUCD1 = 10.0_q/(3.0_q*GAMMA2*(RU+RD)**THRD2*ds)*dalpha_ec
       AMUCD2 = AMUCD1

       Ec = (RU+RD)*ec

       RETURN
    END SUBROUTINE VLAKc
 
 ! TL 08.2024, Correlation energy density of LAK and its derivatives
    SUBROUTINE LAK_corr(rs,zeta,s2,alpha,ec,drs_ec,dzeta_ec,ds2_ec,dalpha_ec)
GPACC !$ACC ROUTINE SEQ
         ! --------------------------------------------
         ! This subroutine calculates the
         ! energy density per particle ec(rs,zeta,s,alpha)
         ! for correlation part of the meta-GGA
         ! --------------------------------------------
         ! intent(in)  :: rs        ! local wigner-seitz radius = (3/(4*pi)*(n_up+n_dn))^(1/3)
         ! intent(in)  :: zeta      ! relative spin-polarisation zeta = (n_up-n_dn)/(n_up+n_dn)
         ! intent(in)  :: s2        ! reduced density gradient |grad n|/(2*gamma*n^(4/3)) squared
         ! intent(in)  :: alpha     ! alpha = (tau-tauW)/tauF
         ! intent(out) :: ec        ! correlation energy per electron
         ! intent(out) :: drs_ec    ! derivative of ec with respect to rs
         ! intent(out) :: dzeta_ec  ! derivative of ec with repect to zeta
         ! intent(out) :: ds2_ec    ! derivative of ec with repect to p = s^2 = 1/(2s)*(derivative of ec with repect to s)
         ! intent(out) :: dalpha_ec ! derivative of ec with repect to alpha

       USE prec
       USE constant
       USE ldalib, only : PW91LDA

       IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(LAK_corr)

       ! In-/Output
       REAL(q) :: rs, zeta, s2, alpha
       REAL(q) :: ec, drs_ec, dzeta_ec, ds2_ec, dalpha_ec
       ! Numbers
       REAL(q) :: ZERO,ONE,TWO,THREE,FOUR,MONE
       REAL(q) :: HALF,THRD,THRD2,THRD4
       ! Composite local quantities
       REAL(q) :: ec0, ec1, fc, ecpw92, drs_ecpw92, dzeta_ecpw92
       ! ALPHA=1
       REAL(q) :: phii, phi2, phi3, dzeta_phi, w1, drs_w1, dzeta_w1
       REAL(q) :: t2, t2_factor, A, gp, dgp, g2p, dg2p, g3p, dg3p, gm1, dgm1, gm2, dgm2, gmw, dgmw ! for H1
       REAL(q) :: H1, dzeta_H1, drs_H1
       REAL(q) :: ds_H1 ! 1/s*(dH1/ds)
       REAL(q) :: drs_ec1, dzeta_ec1, ds_ec1
       ! ALPHA=0
       REAL(q) :: drs_ec0, dzeta_ec0, ds_ec0
       REAL(q) :: eclda0, drs_eclda0, H0, drs_H0, ds_H0, Gc, dzeta_Gc, dx, dzeta_dx, zeta12, zeta11
       REAL(q) :: w0, b1c_drs_w0
       REAL(q) :: ginfty, ds_ginfty
       ! GE2c
       REAL(q) :: Cxc, drs_Cxc
       REAL(q) :: beta_t, drs_beta_t, dbeta_over_beta
       REAL(q) :: beta_a, drs_beta_a
       REAL(q) :: dalpha_fc, drs_fc
       REAL(q) :: alpha_tilde, dalpha_alpha_tilde
       REAL(q) :: ecpw92_nospin, drs_ecpw92_nospin, dzeta_ecpw92_nospin, eclda0rs1
       ! gnum for numerical stability
       REAL(q) :: gnum, ds2_gnum
       REAL(q) :: dummy

       REAL(q), parameter :: gammac = 0.031091_q
       ! Parameters for ALPHA=0
       REAL(q), parameter :: b1c       = 0.0468_q
       REAL(q), parameter :: b2c       = 0.205601_q
       REAL(q), parameter :: chi_infty = 1.55344_q
       ! Parameter for ALPHA=1
       REAL(q), parameter :: ac1  = 10.0_q
       REAL(q), parameter :: b3c       = 2.85_q
       ! GE2c
       REAL(q), parameter :: mu    = 10.0_q/81.0_q
       REAL(q), parameter :: mu_alphax  = -0.209897_q
       REAL(q), parameter :: mu_sx      = (1._q+6.0_q*mu_alphax)*mu    ! mu_s,x
       REAL(q), parameter :: mu_alpha  = 0.104949_q       ! -1/2 mu_alphax
       REAL(q), parameter :: mu_s      = (1._q+6.0_q*mu_alpha)*mu    ! mu_s
       ! Cxc(rs) for mu_alpha = -1/2 mu_alphax
       REAL(q), parameter :: Cxc0  = -1.0782073936365608_q
       REAL(q), parameter :: cxc1  = 0.1_q            ! depends on mu_alpha
       REAL(q), parameter :: cxc2  = 0.065_q          ! depends on mu_alpha
       REAL(q), parameter :: cxc3  = 0.03_q           ! depends on mu_alpha
       REAL(q), parameter :: exc1  = 0.65_q           ! depends on mu_alpha
       REAL(q), parameter :: exc2  = 0.9_q            ! depends on mu_alpha
       REAL(q), parameter :: exc3  = 1.2_q            ! depends on mu_alpha
       REAL(q), parameter :: AC = -0.458165293283143_q ! AC=(4*pi/3)**(-1/3)*AX with AX=-(3/4)*(3/pi)**(1/3)
       REAL(q), parameter :: CT =  1.507303398337901_q ! CT=((3*pi**2)/16)**(2/3)

       ! other parameters
       PARAMETER (ZERO=0.0_q)
       PARAMETER (ONE   =1._q)
       PARAMETER (MONE  =-1._q)
       PARAMETER (TWO   =2._q)
       PARAMETER (THREE =3._q)
       PARAMETER (FOUR  =4._q)
       PARAMETER (HALF =0.5_q)
       PARAMETER (THRD =1._q/3._q)
       PARAMETER (THRD2=2._q*THRD)
       PARAMETER (THRD4=4._q*THRD)
GPACC !$ACC ROUTINE(PW91LDA) SEQ
GPACC !$ACC ROUTINE(calc_gnum) SEQ
GPACC !$ACC ROUTINE(calc_dgnum) SEQ

       call PW91LDA(rs,zeta,ecpw92,drs_ecpw92,dzeta_ecpw92,dummy,dummy)  ! The last two entries are higher derivatives that are not required for LAKc
 
       phii      = half*((ONE+zeta)**THRD2+(ONE-zeta)**THRD2)
       dzeta_phi = THRD*((ONE+zeta)**(MONE*THRD)-(ONE-zeta)**(MONE*THRD))
       phi2      = phii*phii
       phi3      = phi2*phii
 
       w1       = exp(-ecpw92/(gammac*phi3))+MONE
       drs_w1   = -drs_ecpw92/(gammac*phi3)*exp(-ecpw92/(gammac*phi3))
       dzeta_w1 = (THREE*ecpw92*dzeta_phi/phii-dzeta_ecpw92)*exp(-ecpw92/(gammac*phi3))/(gammac*phi3)
          
       ! Cxc(rs) for GE2c
       Cxc     = Cxc0 * (ONE+cxc1*rs**exc1) / ((ONE+cxc2*rs**exc2)*(ONE+cxc3*rs**exc3))
       drs_Cxc = (Cxc0*cxc1*exc1*rs**(exc1-ONE) - Cxc * (cxc2*exc2*rs**(exc2-ONE) + cxc3*exc3*rs**(exc3-ONE) + cxc2*cxc3*(exc2+exc3)*rs**(exc2+exc3-ONE)) ) / ((ONE+cxc2*rs**exc2)*(ONE+cxc3*rs**exc3))
 
       ! ALPHA=1
 
       beta_t          = AC/CT * (mu_s*Cxc - mu_sx)
       drs_beta_t      = AC/CT*mu_s * drs_Cxc
       dbeta_over_beta = drs_beta_t / beta_t
 
       t2_factor = (3.0_q*PI*PI/16.0_q)**THRD2/(phi2*rs)
       t2        = t2_factor*s2
       A   = beta_t/(gammac*w1)
 
       gp   = (ONE + FOUR*A*t2)**(-1.0_q/4.0_q)
       dgp  = -gp/(ONE + FOUR*A*t2)
       g2p   = (ONE + (A*t2)**TWO)**(-ONE)
       dg2p  = -TWO*(A*t2)*g2p**TWO
       g3p   = (ONE + ac1*A*t2)**(-ONE)
       dg3p  = -ac1*g3p/(ONE+ac1*A*t2)
       gm1   = (ONE - FOUR*A*t2)**(-1.0_q/4.0_q)
       dgm1  = gm1/(ONE - FOUR*A*t2)
       gm2  = (ONE + (A*t2)**TWO)**(-ONE)
       dgm2 = -TWO*(A*t2) * gm2**TWO
       gmw  = (ONE - (b3c + w1)*A*t2)**(-ONE)
       dgmw = gmw**TWO
 
       if(beta_t.ge.ZERO) then     ! beta(rs)>=0    This is the case for high to moderate densities
          H1       = gammac*phi3*log(ONE+w1*(ONE-gp)*(ONE-g2p+g3p))
          dzeta_H1 = THREE*dzeta_phi/phii*H1 + &
                       gammac*phi3/(ONE+w1*(ONE-gp)*(ONE-g2p+g3p)) * &
                       (dzeta_w1*(ONE-gp)*(ONE-g2p+g3p) + A*t2*((ONE-g2p+g3p)*dgp+(ONE-gp)*(dg2p-dg3p)) * (dzeta_w1+TWO*w1*dzeta_phi/phii))
          drs_H1   = gammac*phi3/(ONE+w1*(ONE-gp)*(ONE-g2p+g3p)) * &
                       (drs_w1*(ONE-gp)*(ONE-g2p+g3p) + A*t2*((ONE-g2p+g3p)*dgp+(ONE-gp)*(dg2p-dg3p)) * (drs_w1+w1/rs-w1*dbeta_over_beta))
          ! Factor of TWO from ds**2/ds. But no s, because it is divided by, later (in fact ds**2 is needed and everything multiplied by HALF...).
          ds_H1    = -gammac*phi3*w1/(ONE+w1*(ONE-gp)*(ONE-g2p+g3p)) * ((ONE-g2p+g3p)*dgp+(ONE-gp)*(dg2p-dg3p)) *TWO*A*t2_factor
       elseif(beta_t.lt.ZERO) then     ! beta(rs)>=0    Only occurs for low densities (rs>50 for mu_alpha=-1/2 mu_alphax_TASK)
          H1       = gammac*phi3*log(ONE+w1*(ONE-gm1)*(ONE-gm2-gmw))
          dzeta_H1 = THREE*dzeta_phi/phii*H1 + &
                       gammac*phi3/(ONE+w1*(ONE-gm1)*(ONE-gm2-gmw)) * &
                       (dzeta_w1*(ONE-gm1)*(ONE-gm2-gmw) + A*t2*((ONE-gm2-gmw)*dgm1+(ONE-gm1)*(dgm2+(b3c + w1)*dgmw)) * (dzeta_w1+TWO*w1*dzeta_phi/phii) &
                       - w1*A*t2*(ONE-gm1)*dgmw*dzeta_w1 )
          drs_H1   = gammac*phi3/(ONE+w1*(ONE-gm1)*(ONE-gm2-gmw)) * &
                       (drs_w1*(ONE-gm1)*(ONE-gm2-gmw) + A*t2*((ONE-gm2-gmw)*dgm1+(ONE-gm1)*(dgm2+(b3c + w1)*dgmw)) * (drs_w1+w1/rs-w1*dbeta_over_beta) &
                       - w1*A*t2*(ONE-gm1)*dgmw*drs_w1 )
          ds_H1    = -gammac*phi3*w1/(ONE+w1*(ONE-gm1)*(ONE-gm2-gmw)) * ((ONE-gm2-gmw)*dgm1 + (ONE-gm1)*(dgm2+(b3c + w1)*dgmw))* TWO*A*t2_factor
       else ! should not be possible
          H1       = ZERO
          dzeta_H1 = ZERO
          drs_H1   = ZERO
          ds_H1    = ZERO
       endif
 
       ec1       = ecpw92+H1
       drs_ec1   = drs_ecpw92+drs_H1
       dzeta_ec1 = dzeta_ecpw92+dzeta_H1
       ds_ec1    = ds_H1
 
       ! ALPHA=0
 
       eclda0     = -b1c/(ONE+b2c*rs)
       drs_eclda0 = b1c*b2c/(ONE+b2c*rs)**TWO
 
       w0         = exp(-eclda0/b1c)-ONE
       b1c_drs_w0 = -drs_eclda0*exp(-eclda0/b1c)
 
       ! Factor of TWO from ds**2/ds. But no s, because it is divided by, later (in fact ds**2 is needed and everything multiplied by HALF...).
       ginfty      = (ONE+FOUR*chi_infty*s2)**(-1.0_q/4.0_q)
       ds_ginfty   = -TWO*chi_infty*ginfty/(ONE+FOUR*chi_infty*s2)
 
       dx       = half*((ONE+zeta)**THRD4+(ONE-zeta)**THRD4)
       dzeta_dx = THRD2*((ONE+zeta)**THRD-(ONE-zeta)**THRD)
       zeta11   = zeta**11.0_q
       zeta12   = zeta*zeta11
       Gc       = (ONE-2.3631_q*(dx-ONE))*(ONE-zeta12)
       dzeta_Gc = -12.0_q*(ONE-2.3631_q*(dx-ONE))*zeta11-2.3631_q*(ONE-zeta12)*dzeta_dx
 
       H0     = b1c*log(ONE+w0*(ONE-ginfty))
       drs_H0 = (ONE-ginfty)/(ONE+w0*(ONE-ginfty))*b1c_drs_w0
       ds_H0  = -b1c*w0/(ONE+w0*(ONE-ginfty))*ds_ginfty

       ec0       = (eclda0+H0)*Gc
       dzeta_ec0 = (eclda0+H0)*dzeta_Gc
       drs_ec0   = (drs_eclda0+drs_H0)*Gc
       ds_ec0    = ds_H0*Gc
 
       ! fc(ALPHA)
 
       call PW91LDA(rs,ZERO,ecpw92_nospin,drs_ecpw92_nospin,dzeta_ecpw92_nospin,dummy,dummy)  ! zeta=0; The last two entries are higher derivatives that are not required for LAKc
 
       ! fc ArcTan
       beta_a     = -AC * (mu_alpha*Cxc - mu_alphax) / (rs*(ecpw92_nospin-eclda0))
       drs_beta_a = -AC * (mu_alpha*drs_Cxc)/ (rs*(ecpw92_nospin-eclda0)) - beta_a*(drs_ecpw92_nospin-drs_eclda0) / (ecpw92_nospin-eclda0) - beta_a/rs
       dalpha_fc  = ONE / (alpha**TWO/beta_a+beta_a*(PI/TWO*(alpha-ONE))**TWO)                             ! avoid division by alpha^2 for numerical reasons!!  ! This analytical derivative is well-defined AND continuos for alpha<=0 :-)
       drs_fc     = alpha * (alpha-ONE) * drs_beta_a / (alpha**TWO+(PI/TWO*beta_a*(alpha-ONE))**TWO)       ! avoid division by alpha^2 for numerical reasons!!
       
       ! avoid alpha=0
       if(alpha.gt.0._q) then
          alpha_tilde        = (alpha-ONE)/alpha
          fc         = TWO/PI*atan(PI/TWO*beta_a*alpha_tilde)
       elseif(alpha.eq.0._q) then
          fc        = ONE         ! fc is continuosuly extended to alpha=0 by fc(rs,alpha=0) = 1.
       elseif(alpha.lt.0._q) then
          alpha_tilde        = (alpha-ONE)/alpha
          fc         = TWO + TWO/PI*atan(PI/TWO*beta_a*alpha_tilde)     ! 2+fc is the continuous extension of fc(rs,alpha) to alpha<0
       endif

       ! use gnum for numerical stabilization
       call calc_gnum(s2,gnum)
       call calc_dgnum(s2,ds2_gnum)
          
       ec        = ec0 + (ONE-fc) * (ec1-ec0) *gnum
       dalpha_ec = dalpha_fc*(ec0-ec1) *gnum
       drs_ec    = drs_ec0 + ( (ONE-fc)*(drs_ec1-drs_ec0) + drs_fc*(ec0-ec1) ) *gnum
       dzeta_ec  = dzeta_ec0 + (ONE-fc)*(dzeta_ec1-dzeta_ec0) *gnum
       ds2_ec    = HALF*( ds_ec0 + (ONE-fc)*(ds_ec1-ds_ec0) *gnum ) + ds2_gnum* (ONE-fc) * (ec1-ec0)
 
       !Avoiding Divergence
       if(zeta.eq.ONE.or.zeta.eq.MONE) then
          dzeta_ec = ZERO
       endif
    END SUBROUTINE LAK_corr

!************************ SUBROUTINE VCCc ****************************
!
! calculates the first order derivatives of Ec wrt n and |grad(n)|, and tau
! Lebeda et. al. PRR (2022)
!
! Written by Timo Lebeda, August 2024
!
! everything in Hartree units
!
! ATTENTION: Every values are passed "as they are", i.e. including
! possibly unphysical numerical errors (e.g. negative charge densities)
! values need to be checked accordingly (Warning adapted from r2scan implementation)
!
! RU,RD                        density up,down
! DRU, DRD                     abs. val. gradient of density up/down
! DRT                          abs. val. gradient of total density
! TAUU,TAUD                    kinetic energy density up/down
! TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
! VCD1 VCD2                    THE DERIVATIVES OF EC WRT n
! VCDD                         THE DERIVATIVES OF EC WRT |grad n|
! AMUCD1, AMUCD2               THE DERIVATIVES OF EC WRT TAU
!
!***********************************************************************

    SUBROUTINE VCCc(RU,RD,DRT,TAUU,TAUD,Ec,VCD1,VCD2,VCDD,AMUCD1,AMUCD2)
GPACC !$ACC ROUTINE SEQ

       USE prec
       USE constant

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(VCCc)

       ! In-/Output
       REAL(q) :: Ec,DRT,VCDD
       REAL(q) :: RU,TAUU,VCD1,AMUCD1
       REAL(q) :: RD,TAUD,VCD2,AMUCD2
       ! Numbers
       REAL(q) :: ZERO,ONE,MONE,TWO,FOUR
       REAL(q) :: HALF,THRD,THRD2,THRD4,THRD5,THRD8
       ! Composite local quantities
      
       REAL(q) :: rs, zeta, s2, alpha
       REAL(q) :: ds, dzeta_ds, tauF, tauW, tau
       REAL(q) :: drs_ec, dzeta_ec, ds2_ec, dalpha_ec

! other parameters
       PARAMETER (ZERO=0.0_q)
       PARAMETER (ONE   =1._q)
       PARAMETER (MONE  =-1._q)
       PARAMETER (TWO   =2._q)
       PARAMETER (FOUR  =4._q)
       PARAMETER (HALF =0.5_q)
       PARAMETER (THRD =1._q/3._q)
       PARAMETER (THRD2=2._q*THRD)
       PARAMETER (THRD4=4._q*THRD)
       PARAMETER (THRD5=1._q+THRD2)
       PARAMETER (THRD8=1._q+THRD5)
       REAL(q), parameter :: GAMMA  =  3.0936677262801359310_q ! (3*PI^2)^(1/3)
       REAL(q), parameter :: GAMMA2 =  GAMMA*GAMMA
GPACC !$ACC ROUTINE(CC_corr) SEQ
      
       rs   = (THRD4*PI*(RU+RD))**(-THRD)
       zeta = (RU-RD)/(RU+RD)
       s2   = (DRT**TWO)/(4.0_q*GAMMA2*(RU+RD)**THRD8)

       ds       = ONE
       dzeta_ds = ZERO

       tau   = TAUU+TAUD
       tauF  = 3.0_q/10.0_q*GAMMA2*(RU+RD)**THRD5*ds
       tauW  = (DRT**TWO)/(8.0_q*(RU+RD))
       alpha = max((tau-tauW)/tauF,ZERO)

! calculate the correlation enhancement factor and its derivatives wrt rho, zeta, p=s^2, and alpha
       call CC_corr(rs, zeta, s2, alpha, ec, drs_ec, dzeta_ec, ds2_ec, dalpha_ec)

! Calculate ec DERIVATIVES WRT n, |grad n|, and tau
       ! dnUP_ec
       VCD1 = ec - THRD*rs*drs_ec + ( ONE-zeta)*dzeta_ec - THRD8*s2*ds2_ec + (-THRD5*alpha+tauW/tauF+alpha*(zeta+MONE)*dzeta_ds/ds)*dalpha_ec
       ! dnDOWN_ec
       VCD2 = ec - THRD*rs*drs_ec + (MONE-zeta)*dzeta_ec - THRD8*s2*ds2_ec + (-THRD5*alpha+tauW/tauF+alpha*(zeta+ ONE)*dzeta_ds/ds)*dalpha_ec
       
       ! d|n|_ec
       VCDD = TWO*(ds2_ec-5.0_q/3.0_q/ds*dalpha_ec)/(FOUR*GAMMA2*(RU+RD)**THRD5)*DRT
       
       ! dtau_ec
       AMUCD1 = 10.0_q/(3.0_q*GAMMA2*(RU+RD)**THRD2*ds)*dalpha_ec
       AMUCD2 = AMUCD1

       Ec = (RU+RD)*ec

       RETURN
    END SUBROUTINE VCCc
 
 ! TL 08.2024, Correlation energy density of TASKCC and its derivatives
    SUBROUTINE CC_corr(rs,zeta,s2,alpha,ec,drs_ec,dzeta_ec,ds2_ec,dalpha_ec)
GPACC !$ACC ROUTINE SEQ
         ! --------------------------------------------
         ! This subroutine calculates the
         ! energy density per particle ec(rs,zeta,s,alpha)
         ! for correlation part of the TASKCC meta-GGA
         ! --------------------------------------------
         ! intent(in)  :: rs        ! local wigner-seitz radius = (3/(4*pi)*(n_up+n_dn))^(1/3)
         ! intent(in)  :: zeta      ! relative spin-polarisation zeta = (n_up-n_dn)/(n_up+n_dn)
         ! intent(in)  :: s2        ! reduced density gradient |grad n|/(2*gamma*n^(4/3)) squared
         ! intent(in)  :: alpha     ! alpha = (tau-tauW)/tauF
         ! intent(out) :: ec        ! correlation energy per electron
         ! intent(out) :: drs_ec    ! derivative of ec with respect to rs
         ! intent(out) :: dzeta_ec  ! derivative of ec with repect to zeta
         ! intent(out) :: ds2_ec    ! derivative of ec with repect to p = s^2 = 1/(2s)*(derivative of ec with repect to s)
         ! intent(out) :: dalpha_ec ! derivative of ec with repect to alpha

       USE prec
       USE constant
       USE ldalib, only : PW91LDA

       IMPLICIT NONE
GPOMP  !$OMP DECLARE TARGET(CC_corr) 

       ! In-/Output
       REAL(q) :: rs, zeta, s2, alpha
       REAL(q) :: ec, drs_ec, dzeta_ec, ds2_ec, dalpha_ec
       ! Numbers
       REAL(q) :: ZERO,ONE,TWO,MONE
       REAL(q) :: FIFTH3
       ! Composite local quantities
       REAL(q) :: z, ds2_z, dalpha_z ! z=tauW/tau
       REAL(q) :: ecpw92,drs_ecpw92,dzeta_ecpw92
       REAL(q) :: dummy

       ! other parameters
       PARAMETER (ZERO=0.0_q)
       PARAMETER (ONE   =1._q)
       PARAMETER (MONE  =-1._q)
       PARAMETER (TWO   =2._q)
       PARAMETER (FIFTH3 =3._q/5._q)
GPACC !$ACC ROUTINE(PW91LDA) SEQ

       call PW91LDA(rs,zeta,ecpw92,drs_ecpw92,dzeta_ecpw92,dummy,dummy)  ! The last two entries are higher derivatives that are not required for CCc

       z = s2/(FIFTH3*alpha+s2)
       ds2_z = FIFTH3*alpha/(FIFTH3*alpha+s2)**TWO
       dalpha_z = -FIFTH3/(FIFTH3*alpha+s2)*z

       ec        = (ONE-z*zeta*zeta)*ecpw92
       dalpha_ec = -dalpha_z*zeta*zeta*ecpw92
       drs_ec    = (ONE-z*zeta*zeta)*drs_ecpw92
       dzeta_ec  = -TWO*z*zeta*ecpw92+dzeta_ecpw92*(ONE-z*zeta*zeta)
       ds2_ec    = -ds2_z*zeta*zeta*ecpw92

       !Avoiding Divergence
       if(zeta.eq.ONE.or.zeta.eq.MONE) then
          dzeta_ec = ZERO
       endif
    END SUBROUTINE CC_corr

END MODULE MGGALIB
