#include "symbol.inc"

!************************ MODULE GGALIB ***********************
!>
!> This module is the library of subroutines of GGA functionals
!>
!**************************************************************

MODULE GGALIB

      CONTAINS

      SUBROUTINE GGA91_WB(XC,D,DD,EXC,EXCD,EXCDD)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, only : GCOR
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(GGA91_WB)
      TYPE (xc_info) XC
      REAL(q),PARAMETER :: PI=3.1415926536_q
      REAL(q),PARAMETER :: THRD=0.333333333333333_q
GPACC !$ACC ROUTINE(EXCH2) SEQ
GPACC !$ACC ROUTINE(GCOR) SEQ
GPACC !$ACC ROUTINE(CORGGA1) SEQ
      IF (D<0) THEN
       EXC   = 0._q
       EXCD  = 0._q
       EXCDD = 0._q
       RETURN
      ENDIF

      G = 1.0_q
      RS=(0.75_q/(PI*D))**THRD
      FK=(3._q*PI*PI*D)**THRD
      SK = SQRT(4.0_q*FK/PI)
      IF(D>1.E-10_q)THEN
        S=DD/(D*FK*2._q)
        T=DD/(D*SK*2._q)
      ELSE
        S=0.0_q
        T=0.0_q
      ENDIF
      CALL EXCH2(D,S,EX,EXD,EXDD)
      CALL GCOR(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,1.00_q,RS,EU,EURS)
      EC = 0._q
      ECD = 0._q
      CALL CORGGA1(D,RS,T,EC1,EC1D,EC1DD, &
     &             FK,SK,G,EU,EURS,ECZET)
!      WRITE(*,'(10F14.6)') D,S,EX/D,EXD,EXDD
!      WRITE(*,'(10F14.6)') D,RS,SK,T,(EC+EC1)/D,ECD+EC1D,EC1DD
      EXC   = EX*XC%AGGAX  +(EC+EC1)*XC%AGGAC
      EXCD  = EXD*XC%AGGAX +(ECD+EC1D)*XC%AGGAC
      EXCDD = EXDD*XC%AGGAX+EC1DD*XC%AGGAC
      RETURN
      END SUBROUTINE GGA91_WB

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GGA91_WB_SIMD(XC, D_IN, DD, EXC, EXCD, EXCDD, MASK0)
      USE prec
      USE setexm_struct_def
      USE ldalib, only : GCOR_SIMD
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: D_IN, DD
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXC, EXCD, EXCDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: D, G, RS, FK, SK, S, T
      TYPE(SIMD_REAL64)                  :: EX, EXD, EXDD
      TYPE(SIMD_REAL64)                  :: EU, EURS
      TYPE(SIMD_REAL64)                  :: EC1, EC1D, EC1DD
      TYPE(SIMD_REAL64)                  :: ECZET
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: PI = 3.1415926536_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333_q
      REAL(q), PARAMETER                 :: P1 = 0.620350490899400086597382_q   ! (3._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P2 = 3.093667726280135532590521_q   ! (3._q*pi*pi)**(1._q/3._q)
      REAL(q), PARAMETER                 :: P4 = 1.273239544735162764865777_q ! (4._q/pi)

      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         D%X(II) = 1.0E-20_q
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            D%X(II) = D_IN%X(II)
            IF (D_IN%X(II) .GE. 0.0_q) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
            EXC%X(II) = 0.0_q
            EXCD%X(II) = 0.0_q
            EXCDD%X(II) = 0.0_q
         ENDIF
      ENDDO

! return of no SIMD lane continues execution
      IF (.NOT. TRUE_FOR_ANY) RETURN

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(DTHRD)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         G%X(II) = 1.0_q
         DTHRD = SIMD_POW(D%X(II), THRD)
         RS%X(II) = P1 / DTHRD
         FK%X(II) = P2 * DTHRD
         SK%X(II) = SQRT(P4 * FK%X(II))
         IF(D%X(II) .GT. 1.0E-10_q) THEN
            S%X(II) = DD%X(II) / (D%X(II) * FK%X(II) * 2.0_q)
            T%X(II) = DD%X(II) / (D%X(II) * SK%X(II) * 2.0_q)
         ELSE
            S%X(II) = 0.0_q
            T%X(II) = 0.0_q
         ENDIF
      ENDDO

      CALL EXCH2_SIMD(D, S, EX, EXD, EXDD, MASK1)

      CALL GCOR_SIMD(0.0310907_q, 0.21370_q, 7.5957_q, 3.5876_q, 1.6382_q, 0.49294_q, 1.0_q, RS, EU, EURS, MASK1)

      CALL CORGGA1_SIMD(D, RS, T, EC1, EC1D, EC1DD, FK, SK, G, EU, EURS, ECZET, MASK1)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK1%X(II)) THEN
            EXC%X(II) = EX%X(II) * XC%AGGAX + EC1%X(II) * XC%AGGAC
            EXCD%X(II) = EXD%X(II) * XC%AGGAX + EC1D%X(II) * XC%AGGAC
            EXCDD%X(II) = EXDD%X(II) * XC%AGGAX + EC1DD%X(II) * XC%AGGAC
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GGA91_WB_SIMD
#endif

!=======================================================================
!
!>  Exchange according to Perdew and Wang 1991
!
!=======================================================================

      SUBROUTINE EXCH2(D,S,EX,EXD,EXDD)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EXCH2)
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and ites derivatives
!           w.r.t. d (exd) and dd (exdd)
      REAL(q),PARAMETER :: A1=0.19645_q,A2=0.27430_q,A3=0.15084_q,A4=100._q
      REAL(q),PARAMETER :: AX=-0.7385588_q,A=7.7956_q,B1=0.004_q
      REAL(q),PARAMETER :: THPITH=3.0936677262801_q
      THRD = 1._q/3._q
      THRD4 = 4._q/3._q
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.0_q/SQRT(1.0_q+A*A*S2)
      P1 = LOG(A*S+1.0_q/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.0_q/(1.0_q+A1*S*P1+B1*S4)
      P4 = 1.0_q+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*(F-1)*D
!     IF (S==0)
!    &         PRINT *,P1,P2,P3,P4,FAC,F
      P5 = 2.0_q*(S*(A2-A3*P2)+A3*A4*S3*P2-2.0_q*B1*S3)
      P6 = (A1*(P1+A*S*P0)+4.0_q*B1*S3)*((A2-A3*P2)*S2-B1*S4)
      FS = (P5*P3-P6*P3*P3)
      EXD = THRD4*FAC*(F-S*FS-1)
      EXDD = AX*FS*0.5_q/THPITH
      RETURN
      END SUBROUTINE EXCH2

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCH2_SIMD(D, S, EX, EXD, EXDD, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and ites derivatives
!           w.r.t. d (exd) and dd (exdd)
      TYPE(SIMD_REAL64), INTENT(IN)      :: D, S
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EX, EXD, EXDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: A1 = 0.19645_q, A2 = 0.27430_q, A3 = 0.15084_q, A4 = 100.0_q
      REAL(q), PARAMETER                 :: A3_A4 = A3 * A4
      REAL(q), PARAMETER                 :: AX = -0.7385588_q, A = 7.7956_q, A_A = A * A
      REAL(q), PARAMETER                 :: B1 = 0.004_q, TWO_B1 = 2.0_q * B1, FOUR_B1 = 4.0_q * B1
      REAL(q), PARAMETER                 :: THPITH = AX * 0.5_q / 3.0936677262801_q
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q, THRD4 = 4.0_q * THRD

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(FAC, S2, S3, S4, P0, P1, P2, P7, P3, P8, P4, F, P9, P5, P6, FS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         FAC = AX * SIMD_POW(D%X(II), THRD)
         S2 = S%X(II) * S%X(II)
         S3 = S2 * S%X(II)
         S4 = S2 * S2
         P0 = 1.0_q / SQRT(1.0_q + A_A * S2)
         P1 = SIMD_LOG(A * S%X(II) + 1.0_q / P0)
         P2 = SIMD_EXP(-A4 * S2)
         P7 = 1.0_q + A1 * S%X(II) * P1
         P3 = 1.0_q / (P7 + B1 * S4)
         P8 = A2 - A3 * P2
         P4 = P7 + P8 * S2
         F = P3 * P4
         P9 = FAC * (F - 1.0_q)
         P5 = 2.0_q * (S%X(II) * P8 + A3_A4 * S3 * P2 - TWO_B1 * S3)
         P6 = (A1 * (P1 + A * S%X(II) * P0) + FOUR_B1 * S3) * (P8 * S2 - B1 * S4)
         FS = (P5 - P6 * P3) * P3
         IF (MASK0%X(II)) THEN
            EX%X(II) = P9 * D%X(II)
            EXD%X(II) = THRD4 * (P9 - FAC * S%X(II) * FS)
            EXDD%X(II) = FS * THPITH
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCH2_SIMD
#endif

!=======================================================================
!
!>  GGA contribtuion to correlation Perdew and Wang 1991
!
!=======================================================================

      SUBROUTINE CORGGA1(D,RS,T,EC1,EC1D,EC1DD, &
     &                  FK,SK,G,EC,ECRS,ECZET)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORGGA1)
!  gga91 correlation
!  input rs: seitz radius
!  input t: abs(grad d)/(d*2.*ks)
!  output hn: nonlocal part of correlation energy
!  hnd,hndd : derivatives of hn w.r.t. d and dd
      REAL(q),PARAMETER :: XNU=15.75592_q,CC0=0.004235_q,CX=-0.001667212_q,ALF=0.09_q
      REAL(q),PARAMETER :: C1=0.002568_q,C2=0.023266_q,C3=7.389E-6_q,C4=8.723_q
      REAL(q),PARAMETER :: C5=0.472_q,C6=7.389E-2_q,A4=100.0_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,SIXTH7=1.16666666666667_q
      REAL(q),PARAMETER :: PI=3.1415926536_q
      BET = XNU*CC0
      DELT = 2.0_q*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.0_q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.0_q+B*T2
      Q5 = 1.0_q+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.0_q+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.0_q*CX/7.0_q
      R2 = XNU*COEFF*G3
      R3 = EXP(-R1*T2)
      H0 = G3*(BET/DELT)*LOG(1.0_q+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
!============================================================
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      H0T = 2.0_q*BET*T*(1.0_q+2.0_q*B*T2)/Q8
      H0B = -BET*T6*(2.0_q*B+B2*T2)/Q8
      H0RS = H0B*B*ECRS*(B+DELT)/BET
      H1T = 2.0_q*R3*R2*T*(1.0_q-R1*T2)
      CCRS = (C2+2._q*C3*RS)/Q7 - Q6*(C4+2._q*C5*RS+3._q*C6*RS2)/Q7**2
      R1RS = 100.0_q*R0/RS
      H1RS = XNU*T2*R3*(CCRS - COEFF*T2*R1RS)
! = = = = = = = = = = = =
      HT = H0T + H1T
      HRS = H0RS + H1RS
      EC1 = D*H
      EC1D = H-THRD*RS*HRS-SIXTH7*T*HT
      EC1DD = 0.5_q*HT/SK
      RETURN
      END SUBROUTINE CORGGA1

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE CORGGA1_SIMD(D, RS, T, EC1, EC1D, EC1DD, FK, SK, G, EC, ECRS, ECZET, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 correlation
!  input rs: seitz radius
!  input t: abs(grad d)/(d*2.*ks)
!  output hn: nonlocal part of correlation energy
!  hnd,hndd : derivatives of hn w.r.t. d and dd
      TYPE(SIMD_REAL64), INTENT(IN)      :: D, RS, T, FK, SK, G, EC, ECRS, ECZET
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EC1, EC1D, EC1DD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: XNU = 15.75592_q, CC0 = 0.004235_q, CX = -0.001667212_q, ALF = 0.09_q
      REAL(q), PARAMETER                 :: C1 = 0.002568_q, C2 = 0.023266_q, C3 = 7.389E-6_q, C4 = 8.723_q
      REAL(q), PARAMETER                 :: C5 = 0.472_q, C6 = 7.389E-2_q, A4 = 100.0_q
      REAL(q), PARAMETER                 :: THRD = 0.33333333333333_q, SIXTH7 = 1.16666666666667_q
      REAL(q), PARAMETER                 :: PI = 3.1415926536_q
      REAL(q), PARAMETER                 :: BET = XNU * CC0, BET_INV = 1.0_q / BET
      REAL(q), PARAMETER                 :: DELT = 2.0_q * ALF * BET_INV
      REAL(q), PARAMETER                 :: BET_DELT = BET / DELT
      REAL(q), PARAMETER                 :: CX_3_7 = CX * 3.0_q / 7.0_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(G3, G4, PON, B, B2, T2, T4, T6, RS2, RS3, Q4, Q5, Q6, Q7, CC, R0, R1, COEFF, R2, R3, H0, H1, H, Q8, H0T, H0B, H0RS, H1T, CCRS, R1RS, H1RS, HT, HRS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         G3 = SIMD_POW(G%X(II), 3.0_q)
         G4 = G3 * G%X(II)
         PON = -DELT * EC%X(II) * BET_INV / G3
         B = DELT / (SIMD_EXP(PON) - 1.0_q)
         B2 = B * B
         T2 = T%X(II) * T%X(II)
         T4 = T2 * T2
         T6 = T4 * T2
         RS2 = RS%X(II) * RS%X(II)
         RS3 = RS2 * RS%X(II)
         Q4 = 1.0_q + B * T2
         Q5 = Q4 + B2 * T4
         Q6 = C1 + C2 * RS%X(II) + C3 * RS2
         Q7 = 1.0_q / (1.0_q + C4 * RS%X(II) + C5 * RS2 + C6 * RS3)
         CC = -CX + Q6 * Q7
         R0 = SK%X(II) / FK%X(II)
         R0 = R0 * R0
         R1 = A4 * R0 * G4
         COEFF = CC - CC0 - CX_3_7
         R2 = XNU * COEFF * G3
         R3 = SIMD_EXP(-R1 * T2)
         H0 = G3 * (BET_DELT) * LOG(1.0_q + DELT * Q4 * T2 / Q5)
         H1 = R3 * R2 * T2
         H = H0 + H1
!============================================================
         Q8 = 1.0_q / (Q5 * Q5 + DELT * Q4 * Q5 * T2)
         H0T = 2.0_q * BET * T%X(II) * (1.0_q + 2.0_q * B * T2) * Q8
         H0B = -BET * T6 * (2.0_q * B + B2 * T2) * Q8
         H0RS = H0B * B * ECRS%X(II) * (B + DELT) * BET_INV
         H1T = 2.0_q * R3 * R2 *T%X(II) * (1.0_q - R1 * T2)
         CCRS = (C2 + 2.0_q * C3 * RS%X(II)) * Q7 - Q6 * (C4 + 2.0_q * C5 * RS%X(II) + 3.0_q * C6 * RS2) * Q7 * Q7
         R1RS = 100.0_q * R0 / RS%X(II)
         H1RS = XNU * T2 * R3 * (CCRS - COEFF * T2 * R1RS)
!     = = = = = = = = = = = =
         HT = H0T + H1T
         HRS = H0RS + H1RS
         IF (MASK0%X(II)) THEN
            EC1%X(II) = D%X(II) * H
            EC1D%X(II) = H - THRD * RS%X(II) * HRS - SIXTH7 * T%X(II) * HT
            EC1DD%X(II) = 0.5_q * HT / SK%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CORGGA1_SIMD
#endif

!*********************************************************************
!> The P B E
!> this routine was kindly supplied by Bjork Hammer
!> it is essentially identical to Burkes routine but includes
!> the revised P B E routines
!*********************************************************************
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd, &
     &                   ukfactor,PARAM1,PARAM2,ID_XC)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burkes modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxPBE(s)
!       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE xc_name
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EXCHPBE)
      parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
      parameter(pi=3.14159265358979323846264338327950_q)
      parameter(ax=-0.738558766382022405884230032680836_q)
      parameter(um=0.2195149727645171_q,uk1=0.8040_q,ul1=um/uk1)
GPACC !$ACC ROUTINE (CX13) SEQ
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif=AX*rhothrd
      exlda=exunif*rho
      exdlda=exunif*thrd4
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
!----------------------------------------------------------------------
!vdw jk
      IF ((ID_XC==ID_XC_PBE).OR. &
        & (ID_XC==ID_XC_PBE_X).OR. &
        & (ID_XC==ID_XC_PBETEST)) THEN
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
         uk = uk1*ukfactor
         ul = ul1/ukfactor
         P0=1._q+ul*S2
         FxPBE = 1._q+uk-uk/P0
         expbe = exlda*FxPBE
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
         Fs=2._q*um/(P0*P0)
      ELSE IF (ID_XC==ID_XC_RPBE) THEN
! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
         P0=exp(-ul1*S2)
         FxPBE = 1._q+uk1*(1.0_q-P0)
         expbe = exlda*FxPBE
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
         Fs=2._q*um*P0
      ELSE IF (ID_XC==ID_XC_REVPBE) THEN
!revPBE
         uk =1.2450_q
         ul = um/uk
         P0=1._q+ul*S2
         FxPBE = 1._q+uk-uk/P0
         expbe = exlda*FxPBE
         Fs=2._q*um/(P0*P0)
      ELSE IF (ID_XC==ID_XC_OPTPBE) THEN
!opt PBE
         fix=-0.054732_q
         umm=0.175519_q
         uk=1.04804_q
         ul=umm/uk
         P0=1._q+ul*S2
         P1=exp(-ul*S2)
         FxPBE1 = 1._q+uk-(1._q+fix)*uk/P0
         FxRPBE = fix*uk*P1
         FxPBE = (FxPBE1+FxRPBE)
         expbe = exlda*FxPBE
         Fs=(1._q+fix)*2._q*umm/(P0*P0)-2._q*umm*P1*fix
      ELSE IF (ID_XC==ID_XC_OPTB88PBE) THEN
!this is B88 exchange functional
!use 0.18333333333 0.220000000 to get optB88  exchange
       !  BETA = 0.0042_q
       !  x=s*2._q**(4._q/3._q)*(3._q*3.1415926535)**(1._q/3._q)
       !  arsinh(x) = LOG(  x + SQRT( x * x + 1.0_q))
       !  expbe = -BETA * RHOS**(4.0_q/3.0_q) * X*X/(1.0_q + 6.0_q * BETA * X * LOG(  X + SQRT( X * X + 1.0_q)) )
          sabs=abs(s)
          x_s=7.7955541793_q*sabs
!          P0=0.196447965_q*LOG(x_s+sqrt(x_s*x_s+1.0_q))
          P0=PARAM1*LOG(x_s+sqrt(x_s*x_s+1.0_q))
          P1=1._q+sabs*P0
!          FxPBE=(1._q+0.27429447_q*sabs*sabs/P1)
          FxPBE=(1._q+PARAM2*sabs*sabs/P1)
          expbe=exlda*FxPBE
          Fs=PARAM2*2._q/P1 - sabs*PARAM2/(P1**2)*(P0+sabs*PARAM1*7.7955541793_q/sqrt(1._q+x_s*x_s))
      ELSE IF (ID_XC==ID_XC_OPTB86BPBE) THEN
!Becke86MGC
!this is B86MGC with parameters
!PARAM1 is \mu
!PARAM2 is \kappa
!use 0.1234 and 1.00 to get optB86b  exchange
         P0=1._q+PARAM1/PARAM2*S2
         FxPBE = 1._q+PARAM1*S2/(P0**0.8_q)
         expbe = exlda*FxPBE
         Fs=2*PARAM1/(P0**0.8_q)-1.6_q*PARAM1*PARAM1/PARAM2*S2/(P0**(1.8_q))
      ELSE IF (ID_XC==ID_XC_PW86RPBE) THEN
!PerdewWang86 reparametrized by Murray
         S4 = S2*S2
         S6 = S2*S4
         fftnt=1._q/15._q
         P0=1._q+1.851_q*S2+17.33_q*S4+0.163_q*S6
         FxPBE = P0**fftnt
         expbe = exlda*FxPBE
         Fs=fftnt*FxPBE/P0*(3.7020_q+69.32_q*S2+0.978_q*S4)
!vdw jk
      ELSE IF (ID_XC==ID_XC_CX) THEN
!vdw TbB
!CX13 by Berland and Hyldgaard
         CALL CX13(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd)
         RETURN
!vdw TbB
      ENDIF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
! exd = 4/3 * n^1/3 ( 1 + kap - kap/(1+mu/kap s^2) - s^2 2* mu/(1+mu/kap s^2)^2 )
      exd =exunif*THRD4*(FxPBE-S2*Fs)
! this is the part that is multiplied by grad n and its divergence is calculated later on
! exdd = ax/2 /(3pi^2)^1/3 * s * 2 * mu/( 1+mu/kap s^2)^2 
      exdd=0.5_q*ax*0.3232409194_q*S*Fs
      RETURN
      END SUBROUTINE EXCHPBE
!----------------------------------------------------------------------
!----------------------------------------------------------------------

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCHPBE_SIMD(RHO, RHOTHRD, S, EXLDA, EXPBE, EXDLDA, EXD, EXDD, UKFACTOR, MASK0, ID_XC)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burkes modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxPBE(s)
!       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE setexm_struct_def
      USE xc_name
      USE simd
      IMPLICIT REAL(Q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO, RHOTHRD, S, UKFACTOR
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXLDA, EXPBE, EXDLDA, EXD, EXDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: S2, FXPBE, EXUNIF, FS
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q, THRD4 = 4.0_q / 3.0_q
      REAL(q), PARAMETER                 :: PI = 3.14159265358979323846264338327950_q
      REAL(q), PARAMETER                 :: AX = -0.738558766382022405884230032680836_q
      REAL(q), PARAMETER                 :: UM = 0.2195149727645171_q, UK1 = 0.8040_q, UL1 = UM / UK1

!----------------------------------------------------------------------
! construct LDA exchange energy density
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EXUNIF%X(II) = AX * RHOTHRD%X(II)
         IF (MASK0%X(II)) THEN
            EXLDA%X(II) = EXUNIF%X(II) * RHO%X(II)
            EXDLDA%X(II) = EXUNIF%X(II) * THRD4
         ENDIF
!----------------------------------------------------------------------
! construct PBE enhancement factor
         S2%X(II) = S%X(II) * S%X(II)
!----------------------------------------------------------------------
      ENDDO
!vdw jk
      IF ((ID_XC==ID_XC_PBE).OR. &
        & (ID_XC==ID_XC_PBE_X).OR. &
        & (ID_XC==ID_XC_PBETEST)) THEN
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(UK, UL, P0)
         DO II = 0, SIMD_WIDTH_64BIT - 1
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
            UK = UK1 * UKFACTOR%X(II)
            UL = UL1 / UKFACTOR%X(II)
            P0 = 1.0_q / (1.0_q + UL * S2%X(II))
            FXPBE%X(II) = 1.0_q + UK - UK * P0
            IF (MASK0%X(II)) THEN
               EXPBE%X(II) = EXLDA%X(II) * FXPBE%X(II)
            ENDIF
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
            FS%X(II) = 2.0_q * UM * P0 * P0
         ENDDO
      ELSE IF (ID_XC==ID_XC_RPBE) THEN
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(P0)
         DO II = 0, SIMD_WIDTH_64BIT - 1
! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
            P0 = SIMD_EXP(-UL1 * S2%X(II))
               FXPBE%X(II) = 1.0_q + UK1 * (1.0_q - P0)
            IF (MASK0%X(II)) THEN
               EXPBE%X(II) = EXLDA%X(II) * FXPBE%X(II)
            ENDIF
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
            FS%X(II) = 2.0_q * UM * P0
         ENDDO
      ELSE IF (ID_XC==ID_XC_REVPBE) THEN
      ELSE IF (ID_XC==ID_XC_OPTPBE) THEN
      ELSE IF (ID_XC==ID_XC_OPTB88PBE) THEN
      ELSE IF (ID_XC==ID_XC_OPTB86BPBE) THEN
      ELSE IF (ID_XC==ID_XC_PW86RPBE) THEN
      ELSE IF (ID_XC==ID_XC_CX) THEN
      ENDIF
!vdw jk
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            EXD%X(II) = EXUNIF%X(II) * THRD4 * (FXPBE%X(II) - S2%X(II) * FS%X(II))
            EXDD%X(II) = 0.5_q * AX * 0.3232409194_q * S%X(II) * FS%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCHPBE_SIMD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
#endif
!######################################################################
!----------------------------------------------------------------------
![TbB 2014-02-24] Added cx13 by Berland and Hyldgaard
!                 (Phys. Rev. B 89, 035412 (2014))
!----------------------------------------------------------------------
      SUBROUTINE CX13(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd)
GPACC !$ACC ROUTINE SEQ
      use prec
      implicit none
GPOMP !$OMP DECLARE TARGET(CX13)
      real(q) rho,rhothrd,exlda,expbe,exdlda,exd,exdd
      real(q) thrd,thrd4,pi,ax,mu,kappa,alpha,beta,mulv,a,b,c,m
      real(q) s,S2,S3,S4,S5,S6,G1,G2,G3,Fx,Fs,exunif,as6,base
      real(q) fs_rpw86,df_ds,df_rPW86_ds
      parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
      parameter(pi=3.14159265358979323846264338327950_q)
      parameter(ax=-0.738558766382022405884230032680836_q)
      parameter(m=1.0_q/15.0_q,a=0.1234_q,b=17.33_q,c=0.163_q)
      parameter(mulv=.09434_q,alpha=0.021789_q,beta=1.15_q)
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif=AX*rhothrd
      exlda=exunif*rho
      exdlda=exunif*thrd4
!----------------------------------------------------------------------
! construct enhancement factor
      s2 = s*s
      s3 = s2 * s
      s4 = s2 * s2
      s5 = s2 * s3
      s6 = s2 * s2 *s2
      as6 = alpha*s6
! Calculation of energy
      G1 = 1/(1+as6)
      G2 = as6/(beta+as6)
      G3 = alpha*s5/(beta+as6)
      base = (1.0_q+a/m*S2+b*S4+c*S6)
      fs_rpw86=base**m
      Fx = (1+mulv*s2)*G1 + G2*fs_rpw86
      expbe = exlda*Fx
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
      df_rpw86_ds=s*m*fs_rpw86/base*(2.0_q*a/m+4.0_q*b*S2+6.0_q*c*S4)
      Fs = 2*mulv*s*G1 - 6*alpha*s5*(1+mulv*s2)*G1**2 + &
      G2*df_rpw86_ds + 6*G3*(1-G2)*fs_rpw86
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
      exd =exunif*THRD4*(Fx-S*Fs)
      exdd=0.5_q*ax*0.3232409194_q*Fs
      RETURN
      END SUBROUTINE CX13


!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORunspPBE(RS,EC,VC,sk, &
     &                  T,H,DVC,ecdd,lgga)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
       USE prec
       IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORunspPBE)
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
GPACC !$ACC ROUTINE(gcor2) SEQ
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
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EC = EU
! check for zero energy, immediate return if true
      IF (EC==0.) THEN
        H=0; DVC=0; ecdd=0
        RETURN
      ENDIF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS
      VC = EC -RS*ECRS/3._q
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      PON=-EC/(gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = (BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      T6 = T4*T2
      RSTHRD = RS/3._q
      FAC = DELT/B+1._q
      BEC = B2*FAC/(BET)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hT = 2._q*BET*Q9/Q8
      DVC = H+HRS-7.0_q*T2*HT/6._q
      ecdd=0.5_q/sk*t*ht
      RETURN
      END SUBROUTINE CORunspPBE
!----------------------------------------------------------------------

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE CORUNSPPBE_SIMD(RS, EC, VC, SK, T, H, DVC, ECDD, LGGA, MASK0, ID_XC)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
      USE prec
      USE setexm_struct_def
      USE xc_name
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS, SK, T
      TYPE(SIMD_REAL64), INTENT(INOUT)   :: VC
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EC, H, DVC, ECDD
      LOGICAL, INTENT(IN)                :: LGGA
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER, INTENT(IN)                :: ID_XC
      TYPE(SIMD_REAL64)                  :: ECRS, RTRS
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: PON, B, B2, T2, T4, T6, RS2, RS3
      REAL(q)                            :: RSTHRD, FAC, BEC, HB, HRS, HT, TEMP
      REAL(q)                            :: Q_0, Q_1, Q_2, Q_3, Q_4, Q_5, Q_8, Q_9
      INTEGER                            :: II
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4
      LOGICAL                            :: TRUE_FOR_ANY
! numbers for use in LSD energy spin-interpolation formula,  [c](9).
      REAL(q), PARAMETER                 :: GAMMA = 0.03109069086965489503494086371273_q ! (1._q-log(2._q))/pi**2
      REAL(q), PARAMETER                 :: BET = 0.06672455060314922_q ! coefficient in gradient expansion for correlation, [a](4).
      REAL(q), PARAMETER                 :: BETSOL = 0.046_q ! coefficient in gradient expansion for correlation, [a](4).
      REAL(q), PARAMETER                 :: DELT = BET / GAMMA
      REAL(q), PARAMETER                 :: DELTSOL = BETSOL / GAMMA
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q ! 1._q/3._q
      REAL(q), PARAMETER                 :: P_1 = 1.0_q / GAMMA                 ! 1._q/gamma
      REAL(q), PARAMETER                 :: P_2 = 1.0_q / BET                   ! 1._q/bet
      REAL(q), PARAMETER                 :: P_2SOL = 1.0_q / BETSOL             ! 1._q/betsol
      REAL(q), PARAMETER                 :: P_3 = 1.166666666666666740681535_q  ! 7._q/6._q

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EC=unpolarized LSD correlation energy
! ECRS=dEC/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         RTRS%X(II) = SQRT(RS%X(II))
      ENDDO

      CALL GCOR2_SIMD(0.0310907_q, 0.21370_q, 7.5957_q, 3.5876_q, 1.6382_q, 0.49294_q, RTRS, EC, ECRS, MASK0)

      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! check for zero energy, immediate return if true
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF (EC%X(II) .EQ. 0.0_q) THEN
               MASK1%X(II) = .FALSE.
               H%X(II) = 0.0_q
               DVC%X(II) = 0.0_q
               ECDD%X(II) = 0.0_q
            ENDIF
         ENDIF
         IF (MASK1%X(II)) THEN
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
            VC%X(II) = EC%X(II) - RS%X(II) * ECRS%X(II) * THRD
            IF (LGGA) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(PON, TEMP_1, TEMP_2, TEMP_3, TEMP_4, B, B2, T2, T4, RS2, RS3, Q_4, Q_5, T6, RSTHRD, FAC, BEC, Q_8, Q_9, HB, HRS, HT)
      DO II = 0,SIMD_WIDTH_64BIT - 1
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
         PON = -EC%X(II) * P_1
         TEMP_1  =  SIMD_EXP(PON) - 1.0_q
         IF (ID_XC==ID_XC_PBESOL) THEN
            TEMP_2 = DELTSOL
            TEMP_3 = P_2SOL
            TEMP_4 = BETSOL
         ELSE
            TEMP_2 = DELT
            TEMP_3 = P_2
            TEMP_4 = BET
         ENDIF
         B = TEMP_2 / TEMP_1
         B2 = B * B
         T2 = T%X(II) * T%X(II)
         T4 = T2 * T2
         RS2 = RS%X(II) * RS%X(II)
         RS3 = RS2 * RS%X(II)
         Q_4 = 1.0_q + B * T2
         Q_5 = 1.0_q + B * T2 + B2 * T4
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
         T6 = T4 * T2
         RSTHRD = RS%X(II) * THRD
         FAC = TEMP_1 + 1.0_q
         BEC = B2 * FAC * TEMP_3
         Q_8 = Q_5 * Q_5 + TEMP_2 * Q_4 * Q_5 * T2
         Q_9 = 1.0_q + 2.0_q * B * T2
         TEMP_1  =  1.0_q / Q_8
         HB = -TEMP_4 * B * T6 * (2.0_q + B * T2)  *  TEMP_1
         HRS = -RSTHRD * HB * BEC * ECRS%X(II)
         HT = 2.0_q * TEMP_4 * Q_9  *  TEMP_1
         IF (MASK1%X(II)) THEN
            H%X(II) = GAMMA * SIMD_LOG(1.0_q + TEMP_2 * Q_4 * T2 / Q_5)
            DVC%X(II) = H%X(II) + HRS - T2 * HT * P_3
            ECDD%X(II) = 0.5_q / SK%X(II) * T%X(II) * HT
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CORUNSPPBE_SIMD
#endif


!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORPBE)
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
GPACC !$ACC ROUTINE(GCOR2) SEQ
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
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
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
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3._q
      GZ=(((1._q+zet)**2+eta)**sixthm- &
     &((1._q-zet)**2+eta)**sixthm)/3._q
      FAC = DELT/B+1._q
      BG = -3._q*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2._q*BET*G3*Q9/Q8
      COMM = H+HRS-7.0_q*T2*HT/6._q
      PREF = HZ-GZ*T2*HT/G
      COMM = COMM-PREF*ZET
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      ecdd=0.5_q/(sk*g)*t*ht
      RETURN
      END SUBROUTINE CORPBE
!*********************************************************************
!> slimmed down version of GCOR used in PW91 routines, to interpolate
!> LSD correlation energy, as given by (10) of
!> J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!> K. Burke, May 11, 1996.
!*********************************************************************
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(GCOR2)
DOESI IF (rtrs > 1.0E10_q) THEN; GG = 0.0_q; GGRS = 0.0_q; RETURN; END IF ! Breaks ivdw_CF_ddsc (gnu).
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END SUBROUTINE GCOR2

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GCOR2_SIMD(A, A1, B1, B2, B3, B4, RTRS, GG, GGRS, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)       :: RTRS
      TYPE(SIMD_REAL64), INTENT(OUT)      :: GG, GGRS
      TYPE(SIMD_MASK_REAL64), INTENT(IN)  :: MASK0
      INTEGER                             :: II

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(Q0, Q1, Q2, Q3)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         Q0 = -2.0_q * A * (1.0_q + A1 * RTRS%X(II) * RTRS%X(II))
         Q1 = 2.0_q * A * RTRS%X(II) * (B1 + RTRS%X(II) * (B2 + RTRS%X(II) * (B3 + B4 * RTRS%X(II))))
         Q2 = SIMD_LOG(1.0_q + 1.0_q / Q1)
         Q3 = A * (B1 / RTRS%X(II) + 2.0_q * B2 + RTRS%X(II) * (3.0_q * B3 + 4.0_q * B4 * RTRS%X(II)))
         IF (MASK0%X(II)) THEN
            GG%X(II) = Q0 * Q2
            GGRS%X(II) = -2.0_q * A * A1 * Q2 - Q0 * Q3 / (Q1 * (1.0_q + Q1))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GCOR2_SIMD
#endif

!aem AM05 included
!----------------------------------------------------------------------
! Below follow routines needed for the am05 functional (called 'AM'
! in INCAR). Ann E. Mattsson April 2007, spin added May 2008.
!----------------------------------------------------------------------
! $Header:$
!**********************************************************************
!> Armiento Mattsson am05 functional for exchange and correlation
!> Version: 3_spin_xscss
!>
!> input
!>   D1   electron spin-up density [bohr**(-3)]
!>   D2   electron spin-down density [bohr**(-3)]
!>   DD1  abs of gradient of spin-up density (D1) [bohr**(-4)]
!>   DD2  abs of gradient of spin-down density (D2) [bohr**(-4)]
!>   der  integer: 1 = calculate derivatives (FXCD{1,2} and FXDD{1,2})
!>                 0 = don't calculate derivatives
!>
!> output (FXCLDA, VXCLDA{1,2} needed due to special VASP potential handling)
!>   FXC       total exchange-correlation energy density [hartree]
!>             = (2*D1*EX[2*D1] + 2*D2*EX[2*D2])/2 + (D1+D2)*EC[D1,D2]
!>   FXCLDA    PW LDA exchange-correlation energy density [hartree]
!>             = (2*D1*EXLDA[2*D1] + 2*D2*EXLDA[2*D2])/2 + (D1+D2)*ECLDA[D1,D2]
!>   VXCLDA1   LDA exchange-correlation potential for spin-up [hartree]
!>             = d(FXCLDA[D1,D2])/d(D1)
!>   VXCLDA2   LDA exchange-correlation potential for spin-down [hartree]
!>             = d(FXCLDA[D1,D2])/d(D2)
!>   FXCD1     d(FXC[D1,D2])/d(D1) [hartree]
!>   FXCD2     d(FXC[D1,D2])/d(D2) [hartree]
!> In this version both the exchange and the correlation parts are
!> dependent on spin-up/spin-down gradients (DD1/DD2) separately
!>   FXCDD1     d(FXC[D1,D2])/d(DD1) [hartee/(bohr**(-1))]
!>   FXCDD2     d(FXC[D1,D2])/d(DD2) [hartee/(bohr**(-1))]
!>
!> Citation request: when using this functional, please cite:
!> "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005);
!>  A. E. Mattsson and R. Armiento, PRB 79, 155101 (2009)."
!>
!> (The first paper for the AM05 functional, the second for the
!> spin-polarized version)
!>
!> Copyright (c) 2005-2008, Rickard Armiento
!>
!> Permission is hereby granted, free of charge, to any person obtaining
!> a copy of this software and associated documentation files (the
!> "Software"), to deal in the Software without restriction, including
!> without limitation the rights to use, copy, modify, merge, publish,
!> distribute, sublicense, and/or sell copies of the Software, and to
!> permit persons to whom the Software is furnished to do so, subject
!> to the following conditions:
!>
!> The above copyright notice and this permission notice shall be
!> included in all copies or substantial portions of the Software.
!>
!> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!> EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!> OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!> NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!> HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!> WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!> FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!> OTHER DEALINGS IN THE SOFTWARE.
!
!**********************************************************************

      SUBROUTINE AM05(XC,D1,D2,DD1,DD2,der, &
                      FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
                      FXCDD1,FXCDD2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, only : CORLSD
      implicit none
GPOMP !$OMP DECLARE TARGET(AM05)

!     ** Input parameters
      TYPE (xc_info) XC
      REAL(q) D1, D2, DD1, DD2
      INTEGER der

!     ** Output parameters
      REAL(q) FXC, FXCLDA
      REAL(q) VXCLDA1, VXCLDA2
      REAL(q) FXCD1, FXCD2
      REAL(q) FXCDD1, FXCDD2

!     ** Constants and shorthand notation
      REAL(q) pi, g, a, c
      parameter (pi = 3.141592653589793238462643383279502884197_q)
      parameter (g = 0.8098_q, a = 2.804_q)
      parameter (c = 0.7168_q)

!     ** Local variables
      REAL(q) NTOT, KF1, KF2, su, sd, su2, sd2, RS, ZET
      REAL(q) exlda1, exlda2, vxlda1, vxlda2
      REAL(q) eclda, vclda1, vclda2
      REAL(q) fx, fc
      REAL(q) fxd1, fxd2, fxdd1, fxdd2
      REAL(q) fcd1, fcd2, fcdd1, fcdd2
      REAL(q) X1, X2, Xsos1, Xsos2
      REAL(q) F1, F2, Fsos1, Fsos2
      REAL(q) Hx1, Hx2, Hxsos1, Hxsos2
      REAL(q) Hc1, Hc2, Hcsos1, Hcsos2
      REAL(q) szsoz1, szsoz2, denom1, denom2
      REAL(q) zfac1, zfac2, zosn1, zosn2, w1, w2
      REAL(q) ALDAX, ALDAC, AGGAX, AGGAC

!aem  ** Local dummies
      REAL(q) ECLDARS, ECZET, ALFC
GPACC !$ACC ROUTINE (CORLSD) SEQ
GPACC !$ACC ROUTINE (AM05_LABERTW) SEQ

      ALDAX=XC%ALDAX
      ALDAC=XC%ALDAC
      AGGAX=XC%AGGAX
      AGGAC=XC%AGGAC

!     ** Cutoff
      if((D1 .le. 1.e-16_q) .and. (D2 .le. 1.e-16_q)) then
         NTOT=1.e-16_q
      else
         NTOT=D1+D2
      endif

!     *******************
!     LDA PW correlation
!     *******************
      RS = (3._q/(4._q*pi*NTOT))**(1._q/3._q)

!aem  Decided to use the already implemented CORLSD
!aem  which is the Perdew Wang correlation
!aem  (PRB 45, 13244 (1992)) that AM05 needs.
!aem  CORPBE_LDA in xclib.F is a short version of CORLSD
!aem  without the unneeded variables ECLDARS, ECZET, ALFC
      ZET = (D1-D2)/NTOT

      CALL CORLSD(RS,ZET,eclda,vclda1,vclda2,ECLDARS,ECZET,ALFC)

!      WRITE(*,'(2F14.7)') vclda1,vclda2

!     *******************
!     LDA exchange
!     *******************

      vxlda1 = -(3._q*(2._q*D1)/pi)**(1._q/3._q)
      vxlda2 = -(3._q*(2._q*D2)/pi)**(1._q/3._q)
      exlda1 = 3._q/4._q*vxlda1
      exlda2 = 3._q/4._q*vxlda2

!     **********************************
!     LDA exchange-correlation energy
!     and potential
!     **********************************

      FXCLDA = (D1*exlda1 + D2*exlda2)*ALDAX + NTOT*eclda*ALDAC
      VXCLDA1 = vxlda1*ALDAX + vclda1*ALDAC
      VXCLDA2 = vxlda2*ALDAX + vclda2*ALDAC

!     ************************
!     Gradient contributions
!     ************************

!     ********************
!     Spin up
!     ********************

      KF1 = (3._q*pi*pi*(2._q*D1))**(1._q/3._q)
      if ((2._q*KF1*(2._q*D1)) .ge. 1.e-30_q) then
        su = (2._q*DD1)/(2._q*KF1*(2._q*D1))
      else if (DD1 .le. 1.e-30_q) then
        su = 0._q
      else if (DD1 .ge. 1.e-30_q) then
        su = 1.e30_q
      endif
      su2 = su**2

!     ** Interpolation index and

!     ********************
!     Correlation
!     ********************

!     ** Correlation refinement function Eqn. (12)

      if (su .ge. 1.e12_q) then
              X1 = 0._q
              Hc1 = g
      else if (su .le. 1.e-30_q) then
              X1 = 1._q
              Hc1 = 1._q
      else
          X1 = 1._q/(1._q + a*su2)
          Hc1 = X1 + g*(1._q - X1)
      endif

!     ********************
!     Exchange
!     ********************

!     ** Airy LAA refinement function
      CALL AM05_LABERTW(su**(3._q/2._q)/sqrt(24._q),w1)

!     ** am05_lambertw give back argument if it is < 1.0e-20
!     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
!     ** zosn = normalized z/s
      if (su < 1.e-14_q) then
              zosn1 = 1._q
      else
              zosn1 = 24._q**(1._q/3._q)*w1**(2._q/3._q)/su
      end if
      zfac1 = su2*(zosn1*27._q/32._q/pi**2)**2

!     ** denom = denominator of Airy LAA refinement function
      denom1 = 1._q + c*su2*zosn1*(1._q + zfac1)**(1._q/4._q)
!     ** Airy exchange enhancement factor Eqn. (8)
      F1 = (c*su2 + 1._q)/denom1

!     ** Exchange refinement function Eqn. (12)
      Hx1 = X1 + (1._q - X1)*F1

!     ********************
!     Spin down
!     ********************

      KF2 = (3._q*pi*pi*(2._q*D2))**(1._q/3._q)
      if ((2._q*KF2*(2._q*D2)) .ge. 1.e-30_q) then
          sd = (2._q*DD2)/(2._q*KF2*(2._q*D2))
      else if (DD2 .le. 1.e-30_q) then
          sd = 0._q
      else if (DD2 .ge. 1.e-30_q) then
          sd = 1.e30_q
      endif
      sd2 = sd**2

!     ** Interpolation index and

!     ********************
!     Correlation
!     ********************

!     ** Correlation refinement function Eqn. (12)

      if (sd .ge. 1.e12_q) then
          X2 = 0._q
          Hc2 = g
      else if (sd .le.  1.e-30_q) then
          X2 = 1._q
          Hc2 = 1._q
      else
          X2 = 1._q/(1._q + a*sd2)
          Hc2 = X2 + g*(1._q - X2)
      endif

!     ********************
!     Exchange
!     ********************

!     ** Airy LAA refinement function
      CALL AM05_LABERTW(sd**(3._q/2._q)/sqrt(24._q),w2)

!     ** am05_lambertw give back argument if it is < 1.0e-20
!     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
!     ** zosn = normalized z/s
      if (sd < 1.e-14_q) then
              zosn2 = 1._q
      else
              zosn2 = 24._q**(1._q/3._q)*w2**(2._q/3._q)/sd
      end if
      zfac2 = sd2*(zosn2*27._q/32._q/pi**2)**2

!     ** denom = denominator of Airy LAA refinement function
      denom2 = 1._q + c*sd2*zosn2*(1._q + zfac2)**(1._q/4._q)
!     ** Airy exchange enhancement factor Eqn. (8)
      F2 = (c*sd2 + 1._q)/denom2

!     ** Exchange refinement function Eqn. (12)
      Hx2 = X2 + (1._q - X2)*F2

!     ********************************
!     Exchange-correlation energy
!     ********************************

!     ** Exchange energy density, Ex = Integrate[fx]
      fx = D1*exlda1*Hx1 + D2*exlda2*Hx2

!     ** Correlation energy density, Ec = Integrate[fc]
      fc = eclda*(D1*Hc1+D2*Hc2)

      FXC = (fx - D1*exlda1 - D2*exlda2)*AGGAX + (D1*exlda1 + D2*exlda2)*ALDAX &
     &    + (fc - NTOT*eclda)*AGGAC + NTOT*eclda*ALDAC

      if (der .eq. 0) return

!     *******************************
!     Potentials
!     *******************************

!     **************************
!     Spin up
!     **************************

      if ((D1 .le. 1.e-16_q) .AND. (su .le. 1.e-30_q)) then
          fxd1 = 0._q
          fxdd1 = 0._q
          fcd1 = vclda1*Hc2 + eclda*(Hc1 - Hc2)
          fcdd1 = 0._q
      else
          if (D1 .le. 1.e-16_q) then
              D1 = 1.e-16_q
          endif

!     ***************************
!         Exchange derivatives
!     ***************************

!     ** Interpolation index derivatives, 1/s dX/ds
          Xsos1 = -2._q*a*X1**2

!     ** Airy LAA refinement function derivatives, 1/s dF/ds
!     ** szsoz = s*(dz/ds)/z
          szsoz1 = 1._q/(1._q + w1)

          Fsos1 = c/denom1**2*(2._q - zosn1* &
     &            ((1._q - c*su2)*(1._q + zfac1)**(1._q/4._q) + &
     &             (1._q + c*su2)*(1._q + 3._q/2._q*zfac1)/ &
     &                      (1._q + zfac1)**(3._q/4._q)*szsoz1))

!     ** Refinement function derivatives, 1/s dHx/ds
!     ** We use that (1 - X) = a*X*s2
          Hxsos1 = (1._q - X1)*Fsos1 - (F1 - 1._q)*Xsos1

          fxd1 = vxlda1*Hx1 - 4._q/3._q*exlda1*su2*Hxsos1
          fxdd1 = exlda1*su*Hxsos1/2._q/KF1

!     *****************************
!        Correlation derivatives
!     *****************************

!     ** Correlation refinement function derivatives, 1/s dF/ds
          Hcsos1 = Xsos1*(1._q - g)

          fcd1 = vclda1/NTOT*(D1*Hc1 + D2*Hc2) &
     &           + eclda*D2/NTOT*(Hc1 - Hc2) &
     &           - 4._q/3._q*eclda*su2*Hcsos1
          fcdd1 = eclda*su*Hcsos1/2._q/KF1
      endif

!     **************************
!     Spin down
!     **************************

      if ((D2 .le. 1.e-16_q) .AND. (sd .le. 1.e-30_q)) then
          fxd2 = 0._q
          fxdd2 = 0._q
          fcd2 = vclda2*Hc1 + eclda*(Hc2 - Hc1)
          fcdd2 = 0._q
      else
          if (D2 .le.  1.e-16_q) then
             D2 = 1.e-16_q
          endif

!     ***************************
!         Exchange derivatives
!     ***************************

!     ** Interpolation index derivatives, 1/s dX/ds
          Xsos2 = -2._q*a*X2**2

!     ** Airy LAA refinement function derivatives, 1/s dF/ds
!     ** szsoz = s*(dz/ds)/z
          szsoz2 = 1._q/(1._q + w2)

          Fsos2 = c/denom2**2*(2._q - zosn2* &
     &            ((1._q - c*sd2)*(1._q + zfac2)**(1._q/4._q) + &
     &             (1._q + c*sd2)*(1._q + 3._q/2._q*zfac2)/ &
     &                      (1._q + zfac2)**(3._q/4._q)*szsoz2))

!     ** Refinement function derivatives, 1/s dHx/ds
!     ** We use that (1 - X) = a*X*s2
          Hxsos2 = (1._q - X2)*Fsos2 - (F2 - 1._q)*Xsos2

          fxd2 = vxlda2*Hx2 - 4._q/3._q*exlda2*sd2*Hxsos2
          fxdd2 = exlda2*sd*Hxsos2/2._q/KF2

!     *****************************
!        Correlation derivatives
!     *****************************

!     ** Correlation refinement function derivatives, 1/s dF/ds
          Hcsos2 = Xsos2*(1._q - g)

          fcd2 = vclda2/NTOT*(D1*Hc1+D2*Hc2) &
     &           + eclda*D1/NTOT*(Hc2 - Hc1) &
     &           - 4._q/3._q*eclda*sd2*Hcsos2
          fcdd2 = eclda*sd*Hcsos2/2._q/KF2
      endif

!     *************************************
!        Exchange-correlation derivatives
!     *************************************

      FXCD1 = (fxd1 - vxlda1)*AGGAX + vxlda1*ALDAX + (fcd1 - vclda1)*AGGAC + vclda1*ALDAC
      FXCD2 = (fxd2 - vxlda2)*AGGAX + vxlda2*ALDAX + (fcd2 - vclda2)*AGGAC + vclda2*ALDAC
      FXCDD1 = fxdd1*AGGAX + fcdd1*AGGAC
      FXCDD2 = fxdd2*AGGAX + fcdd2*AGGAC

      return
      END SUBROUTINE AM05

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE AM05_SIMD(XC, D1_IN, D2_IN, DD1, DD2, DER, FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2, MASK0)
      USE prec
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY : CORLSD_SIMD
      IMPLICIT NONE
      TYPE(xc_info)                      :: XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: D1_IN, D2_IN, DD1, DD2
      TYPE(SIMD_INT64), INTENT(IN)       :: DER
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64), INTENT(OUT)     :: FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2
      TYPE(SIMD_REAL64)                  :: D1, D2
      TYPE(SIMD_REAL64)                  :: NTOT, RS, ZET, ECLDA, VCLDA1, VCLDA2, ECLDARS, ECZET, ALFC
      TYPE(SIMD_REAL64)                  :: KF1, KF2, SU, SD, W1, W2
      TYPE(SIMD_REAL64)                  :: ARG_1, ARG_2
      TYPE(SIMD_REAL64)                  :: ZOSN1, ZOSN2, ZFAC1, ZFAC2, DENOM1, DENOM2, F1, F2, X1, X2
      TYPE(SIMD_REAL64)                  :: HC1, HC2, HX1, HX2, EXLDA1, EXLDA2, VXLDA1, VXLDA2
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: SU2, SD2, NTOT_INV, KF1_INV, KF2_INV
      REAL(q)                            :: FX, FC
      REAL(q)                            :: FXD1, FXD2, FXDD1, FXDD2
      REAL(q)                            :: FCD1, FCD2, FCDD1, FCDD2
      REAL(q)                            :: XSOS1, XSOS2
      REAL(q)                            :: FSOS1, FSOS2
      REAL(q)                            :: HXSOS1, HXSOS2
      REAL(q)                            :: HCSOS1, HCSOS2
      REAL(q)                            :: SZSOZ1, SZSOZ2
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: PI = 3.141592653589793238462643383279502884197_q
      REAL(q), PARAMETER                 :: G = 0.8098_q, A = 2.804_q, C = 0.7168_q
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q
      REAL(q), PARAMETER                 :: THRD2 = 2.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P1 = (3.0_q / (4.0_q * PI)) ** THRD
      REAL(q), PARAMETER                 :: P2 = (6.0_q * PI * PI) ** THRD
      REAL(q), PARAMETER                 :: P3 = 1.0_q / SQRT(24.0_q)
      REAL(q), PARAMETER                 :: P4 = 24.0_q ** THRD
      REAL(q), PARAMETER                 :: P5 = 27.0_q / (32.0_q * PI * PI)
      REAL(q), PARAMETER                 :: P6 = (6.0_q / PI) ** THRD
      REAL(q), PARAMETER                 :: P7 = 4.0_q / 3.0_q
      
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, NTOT_INV)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         D1%X(II) = D1_IN%X(II)
         D2%X(II) = D2_IN%X(II)

!     ** CUTOFF
         IF((D1%X(II) .LE. 1.0E-16_q) .AND. (D2%X(II) .LE. 1.0E-16_q)) THEN
            TEMP_1 = 1.0E-16_q
         ELSE
            TEMP_1 = D1%X(II) + D2%X(II)
         ENDIF

         NTOT%X(II) = TEMP_1
         NTOT_INV = 1.0_q / TEMP_1

!     *******************
!     LDA PW CORRELATION (PART 1)
!     *******************
         RS%X(II) = P1 * SIMD_POW(NTOT_INV, THRD)

!AEM  DECIDED TO USE THE ALREADY IMPLEMENTED CORLSD
!AEM  WHICH IS THE PERDEW WANG CORRELATION
!AEM  (PRB 45, 13244 (1992)) THAT AM05 NEEDS.
!AEM  CORPBE_LDA IN XCLIB.F IS A SHORT VERSION OF CORLSD
!AEM  WITHOUT THE UNNEEDED VARIABLES ECLDARS, ECZET, ALFC
         ZET%X(II) = (D1%X(II) - D2%X(II)) * NTOT_INV

!     ************************
!     GRADIENT CONTRIBUTIONS (PART 1)
!     ************************

!     ********************
!     SPIN UP
!     ********************

         KF1%X(II) = P2 * SIMD_POW(D1%X(II), THRD)
         TEMP_1 = 4.0_q * KF1%X(II) * D1%X(II)
         IF (TEMP_1 .GE. 1.0E-30_q) THEN
            SU%X(II) = (2.0_q * DD1%X(II)) / TEMP_1
         ELSE IF (DD1%X(II) .LE. 1.0E-30_q) THEN
            SU%X(II) = 0.0_q
         ELSE IF (DD1%X(II) .GE. 1.0E-30_q) THEN
            SU%X(II) = 1.0E30_q
         ENDIF

!     ********************
!     SPIN DOWN
!     ********************

         KF2%X(II) = P2 * SIMD_POW(D2%X(II), THRD)
         TEMP_1 = 4.0_q * KF2%X(II) * D2%X(II)
         IF (TEMP_1 .GE. 1.0E-30_q) THEN
            SD%X(II) = (2.0_q * DD2%X(II)) / TEMP_1
         ELSE IF (DD2%X(II) .LE. 1.0E-30_q) THEN
            SD%X(II) = 0.0_q
         ELSE IF (DD2%X(II) .GE. 1.0E-30_q) THEN
            SD%X(II) = 1.0E30_q
         ENDIF

         ARG_1%X(II) = SIMD_POW(SU%X(II), 1.5_q) * P3
         ARG_2%X(II) = SIMD_POW(SD%X(II), 1.5_q) * P3
      ENDDO

!     ************************
!     GRADIENT CONTRIBUTIONS (PART 2)
!     ************************

!     *******************
!     LDA PW CORRELATION
!     *******************
      CALL CORLSD_SIMD(RS, ZET, ECLDA, VCLDA1, VCLDA2, ECLDARS, ECZET, ALFC, MASK0)

!     ********************
!     EXCHANGE
!     ********************

!     ** AIRY LAA REFINEMENT FUNCTION
      CALL AM05_LABERTW_SIMD(ARG_1, W1, MASK0)

!     ** AIRY LAA REFINEMENT FUNCTION
      CALL AM05_LABERTW_SIMD(ARG_2, W2, MASK0)

      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(SU2, SD2, TEMP_1, TEMP_2, FX, FC)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         SU2 = SU%X(II) * SU%X(II)
         SD2 = SD%X(II) * SD%X(II)

!     ** AM05_LAMBERTW GIVE BACK ARGUMENT IF IT IS < 1.0E-20
!     ** (1.0E-14)^{3/2} = 1.0E-21 => GIVE  LOW S LIMIT FOR Z/S
!     ** ZOSN = NORMALIZED Z/S
         IF (SU%X(II) < 1.0E-14_q) THEN
            TEMP_1 = 1.0_q
         ELSE
            TEMP_1 = P4 * SIMD_POW(W1%X(II), THRD2) / SU%X(II)
         END IF
         ZOSN1%X(II) = TEMP_1
         TEMP_2 = P5 * TEMP_1
         ZFAC1%X(II) = SU2 * TEMP_2 * TEMP_2
         
!     ** DENOM = DENOMINATOR OF AIRY LAA REFINEMENT FUNCTION
         DENOM1%X(II) = 1.0_q / (1.0_q + C * SU2 * TEMP_1 * SIMD_POW(1.0_q + ZFAC1%X(II), 0.25_q))
!     ** AIRY EXCHANGE ENHANCEMENT FACTOR EQN. (8)
         F1%X(II) = (C * SU2 + 1.0_q) * DENOM1%X(II)

!     ********************
!     CORRELATION
!     ********************

!     ** CORRELATION REFINEMENT FUNCTION EQN. (12)

         IF (SU%X(II) .GE. 1.0E12_q) THEN
            TEMP_1 = 0.0_q
            TEMP_2 = G
         ELSE IF (SU%X(II) .LE. 1.0E-30_q) THEN
            TEMP_1 = 1.0_q
            TEMP_2 = 1.0_q
         ELSE
            TEMP_1 = 1.0_q / (1.0_q + A * SU2)
            TEMP_2 = TEMP_1 + G * (1.0_q - TEMP_1)
         ENDIF
         X1%X(II) = TEMP_1
         HC1%X(II) = TEMP_2

!     ** EXCHANGE REFINEMENT FUNCTION EQN. (12)
         HX1%X(II) = TEMP_1 + (1.0_q - TEMP_1) * F1%X(II)
         
!     ** INTERPOLATION INDEX AND

!     ********************
!     CORRELATION
!     ********************

!     ** CORRELATION REFINEMENT FUNCTION EQN. (12)

         IF (SD%X(II) .GE. 1.0E12_q) THEN
            TEMP_1 = 0.0_q
            TEMP_2 = G
         ELSE IF (SD%X(II) .LE.  1.0E-30_q) THEN
            TEMP_1 = 1.0_q
            TEMP_2 = 1.0_q
         ELSE
            TEMP_1 = 1.0_q / (1.0_q + A * SD2)
            TEMP_2 = TEMP_1 + G * (1.0_q - TEMP_1)
         ENDIF
         X2%X(II) = TEMP_1
         HC2%X(II) = TEMP_2

!     ********************
!     EXCHANGE
!     ********************

!     ** AM05_LAMBERTW GIVE BACK ARGUMENT IF IT IS < 1.0E-20
!     ** (1.0E-14)^{3/2} = 1.0E-21 => GIVE  LOW S LIMIT FOR Z/S
!     ** ZOSN = NORMALIZED Z/S
         IF (SD%X(II) < 1.0E-14_q) THEN
            TEMP_1 = 1.0_q
         ELSE
            TEMP_1 = P4 * SIMD_POW(W2%X(II), THRD2) / SD%X(II)
         END IF
         ZOSN2%X(II) = TEMP_1
         TEMP_2 = P5 * TEMP_1
         ZFAC2%X(II) = SD2 * TEMP_2 * TEMP_2
         
!     ** DENOM = DENOMINATOR OF AIRY LAA REFINEMENT FUNCTION
         DENOM2%X(II) = 1.0_q / (1.0_q + C * SD2 * TEMP_1 * SIMD_POW(1.0_q + ZFAC2%X(II), 0.25_q))
!     ** AIRY EXCHANGE ENHANCEMENT FACTOR EQN. (8)
         F2%X(II) = (C * SD2 + 1.0_q) * DENOM2%X(II)

!     ** EXCHANGE REFINEMENT FUNCTION EQN. (12)
         HX2%X(II) = X2%X(II) + (1.0_q - X2%X(II)) * F2%X(II)

!     *******************
!     LDA EXCHANGE
!     *******************

         VXLDA1%X(II) = -P6 * SIMD_POW(D1%X(II), THRD)
         VXLDA2%X(II) = -P6 * SIMD_POW(D2%X(II), THRD)
         EXLDA1%X(II) = 0.75_q * VXLDA1%X(II)
         EXLDA2%X(II) = 0.75_q * VXLDA2%X(II)

!     ********************************
!     EXCHANGE-CORRELATION ENERGY
!     ********************************

!     ** EXCHANGE ENERGY DENSITY, EX = INTEGRATE[FX]
         FX = D1%X(II) * EXLDA1%X(II) * HX1%X(II) + D2%X(II) * EXLDA2%X(II) * HX2%X(II)

!     ** CORRELATION ENERGY DENSITY, EC = INTEGRATE[FC]
         FC = ECLDA%X(II) * (D1%X(II) * HC1%X(II) + D2%X(II) * HC2%X(II))

         IF (MASK0%X(II)) THEN
            FXC%X(II) = (FX - D1%X(II) * EXLDA1%X(II) - D2%X(II) * EXLDA2%X(II)) * XC%AGGAX &
     &                + (D1%X(II) * EXLDA1%X(II) + D2%X(II) * EXLDA2%X(II)) * XC%ALDAX &
     &                + (FC - NTOT%X(II) * ECLDA%X(II)) * XC%AGGAC +  NTOT%X(II) * ECLDA%X(II) * XC%ALDAC

!     **********************************
!     LDA EXCHANGE-CORRELATION ENERGY
!     AND POTENTIAL
!     **********************************

            FXCLDA%X(II) = (D1%X(II) * EXLDA1%X(II) + D2%X(II) * EXLDA2%X(II)) * XC%ALDAX &
     &                   + NTOT%X(II) * ECLDA%X(II) * XC%ALDAC
            VXCLDA1%X(II) = VXLDA1%X(II) * XC%ALDAX + VCLDA1%X(II) * XC%ALDAC
            VXCLDA2%X(II) = VXLDA2%X(II) * XC%ALDAX + VCLDA2%X(II) * XC%ALDAC
         ENDIF

         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF (DER%X(II) .NE. 0) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(SU2, SD2, NTOT_INV, FXD1, FXDD1, FCD1, FCDD1, XSOS1, SZSOZ1, TEMP_1, TEMP_2, FSOS1, HXSOS1, KF1_INV, HCSOS1, FXD2, FXDD2, FCD2, FCDD2, XSOS2, SZSOZ2, FSOS2, HXSOS2, KF2_INV, HCSOS2)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         SU2 = SU%X(II) * SU%X(II)
         SD2 = SD%X(II) * SD%X(II)
         NTOT_INV = 1.0_q / NTOT%X(II)

!     *******************************
!     POTENTIALS
!     *******************************

!     **************************
!     SPIN UP
!     **************************

         IF ((D1%X(II) .LE. 1.0E-16_q) .AND. (SU%X(II) .LE. 1.0E-30_q)) THEN
            FXD1 = 0.0_q
            FXDD1 = 0.0_q
            FCD1 = VCLDA1%X(II) * HC2%X(II) + ECLDA%X(II) * (HC1%X(II) - HC2%X(II))
            FCDD1 = 0.0_q
         ELSE
            IF (D1%X(II) .LE. 1.0E-16_q) THEN
               D1%X(II) = 1.0E-16_q
            ENDIF

!     ***************************
!         EXCHANGE DERIVATIVES
!     ***************************

!     ** INTERPOLATION INDEX DERIVATIVES, 1/S DX/DS
            XSOS1 = -2.0_q * A * X1%X(II) * X1%X(II)

!     ** AIRY LAA REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
!     ** SZSOZ = S*(DZ/DS)/Z
            SZSOZ1 = 1.0_q / (1.0_q + W1%X(II))
            TEMP_1 = SIMD_POW(1.0_q + ZFAC1%X(II), 0.25_q)
            TEMP_2 = TEMP_1 * TEMP_1 * TEMP_1
            FSOS1 = C * DENOM1%X(II) * DENOM1%X(II) * (2.0_q - ZOSN1%X(II) * &
     &           ((1.0_q - C * SU2) * TEMP_1 + &
     &            (1.0_q + C * SU2) * (1.0_q + 1.5_q * ZFAC1%X(II)) / TEMP_2 * SZSOZ1))
            
!     ** REFINEMENT FUNCTION DERIVATIVES, 1/S DHX/DS
!     ** WE USE THAT (1 - X) = A*X*S2
            HXSOS1 = (1.0_q - X1%X(II)) * FSOS1 - (F1%X(II) - 1.0_q) * XSOS1
            KF1_INV = 1.0_q / KF1%X(II)
            FXD1 = VXLDA1%X(II) * HX1%X(II) - P7 * EXLDA1%X(II) * SU2 * HXSOS1
            FXDD1 = 0.5_q * EXLDA1%X(II) * SU%X(II) * HXSOS1 * KF1_INV

!     *****************************
!        CORRELATION DERIVATIVES
!     *****************************

!     ** CORRELATION REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
            HCSOS1 = XSOS1 * (1.0_q - G)
            FCD1 = VCLDA1%X(II) * NTOT_INV * (D1%X(II) * HC1%X(II) + D2%X(II) * HC2%X(II)) &
     &           + ECLDA%X(II) * D2%X(II) * NTOT_INV * (HC1%X(II) - HC2%X(II)) &
     &           - P7 * ECLDA%X(II) * SU2 * HCSOS1
            FCDD1 = 0.5_q * ECLDA%X(II) * SU%X(II) * HCSOS1 * KF1_INV
         ENDIF

!     **************************
!     SPIN DOWN
!     **************************

         IF ((D2%X(II) .LE. 1.0E-16_q) .AND. (SD%X(II) .LE. 1.0E-30_q)) THEN
            FXD2 = 0.0_q
            FXDD2 = 0.0_q
            FCD2 = VCLDA2%X(II) * HC1%X(II) + ECLDA%X(II) * (HC2%X(II) - HC1%X(II))
            FCDD2 = 0.0_q
         ELSE
            IF (D2%X(II) .LE.  1.0E-16_q) THEN
               D2%X(II) = 1.0E-16_q
            ENDIF
            
!     ***************************
!         EXCHANGE DERIVATIVES
!     ***************************

!     ** INTERPOLATION INDEX DERIVATIVES, 1/S DX/DS
            XSOS2 = -2.0_q * A * X2%X(II) * X2%X(II)

!     ** AIRY LAA REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
!     ** SZSOZ = S*(DZ/DS)/Z
            SZSOZ2 = 1.0_q / (1.0_q + W2%X(II))
            TEMP_1 = SIMD_POW(1.0_q + ZFAC2%X(II), 0.25_q)
            TEMP_2 = TEMP_1 * TEMP_1 * TEMP_1
            FSOS2 = C * DENOM2%X(II) * DENOM2%X(II) * (2.0_q - ZOSN2%X(II) * &
     &           ((1.0_q - C * SD2) * TEMP_1 + &
     &            (1.0_q + C * SD2) * (1.0_q + 1.5_q * ZFAC2%X(II)) / TEMP_2 * SZSOZ2))

!     ** REFINEMENT FUNCTION DERIVATIVES, 1/S DHX/DS
!     ** WE USE THAT (1 - X) = A*X*S2
            HXSOS2 = (1.0_q - X2%X(II)) * FSOS2 - (F2%X(II) - 1.0_q) * XSOS2
            KF2_INV = 1.0_q / KF2%X(II)
            FXD2 = VXLDA2%X(II) * HX2%X(II) - P7 * EXLDA2%X(II) * SD2 * HXSOS2
            FXDD2 = 0.5_q * EXLDA2%X(II) * SD%X(II) * HXSOS2 * KF2_INV

!     *****************************
!        CORRELATION DERIVATIVES
!     *****************************

!     ** CORRELATION REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
            HCSOS2 = XSOS2 * (1.0_q - G)
            FCD2 = VCLDA2%X(II) * NTOT_INV * (D1%X(II) * HC1%X(II) + D2%X(II) * HC2%X(II)) &
     &           + ECLDA%X(II) * D1%X(II) * NTOT_INV *(HC2%X(II) - HC1%X(II)) &
     &           - P7 * ECLDA%X(II) * SD2 * HCSOS2
            FCDD2 = 0.5_q * ECLDA%X(II) * SD%X(II) * HCSOS2 * KF2_INV
         ENDIF

!     *************************************
!        EXCHANGE-CORRELATION DERIVATIVES
!     *************************************

         IF (MASK1%X(II)) THEN
            FXCD1%X(II) = (FXD1 - VXLDA1%X(II)) * XC%AGGAX + VXLDA1%X(II) * XC%ALDAX &
     &                  + (FCD1 - VCLDA1%X(II)) * XC%AGGAC + VCLDA1%X(II) * XC%ALDAC
            FXCD2%X(II) = (FXD2 - VXLDA2%X(II)) * XC%AGGAX + VXLDA2%X(II) * XC%ALDAX &
     &                  + (FCD2 - VCLDA2%X(II)) * XC%AGGAC + VCLDA2%X(II) * XC%ALDAC
            FXCDD1%X(II) = FXDD1 * XC%AGGAX + FCDD1 * XC%AGGAC
            FXCDD2%X(II) = FXDD2 * XC%AGGAX + FCDD2 * XC%AGGAC
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE AM05_SIMD
#endif

!     ***********************************************
!>       LambertW function.
!>
!>       Corless, Gonnet, Hare, Jeffrey, and Knuth (1996),
!>         Adv. in Comp. Math. 5(4):329-359.
!>       Implementation approach loosely inspired by the
!>       GNU Octave version by N. N. Schraudolph, but this
!>       implementation is only for real values and
!>       principal branch.
!>
!>       Copyright (c) 2005, Rickard Armiento
!>       All rights reserved.
!     ***********************************************

      SUBROUTINE AM05_LABERTW(z,result)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE tutor, ONLY : vtutor
      implicit none
GPOMP !$OMP DECLARE TARGET(AM05_LABERTW)

!     input
      REAL(q) z
!     output
      REAL(q) result
!     local variables
      REAL(q) e,t,p
      INTEGER i

!     ** If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.E-20_q) then
         result = z
         return
      endif

      e = exp(1._q)

!     ** Inital guess
      if( abs(z + 1._q/e) .gt. 1.45_q ) then
!        ** Asymptotic expansion at 0 and Inf
         result = log(z)
         result = result - log(result)
      else
!        ** Series expansion about -1/e to first order
         result = 1._q*sqrt(2._q*e*z + 2._q) - 1._q
      endif

!     ** Find result through iteration
      do i=1,10
         p = exp(result)
         t = result*p - z
         if( result .ne. -1._q ) then
            t = t/(p*(result + 1._q) - &
     &           0.5_q*(result + 2._q)*t/(result + 1._q))
         else
            t = 0._q
         endif
         result = result - t
         if(abs(t) < (2.48_q*1.E-14_q)*(1._q + abs(result))) then
            return
         endif
      enddo
!     ** This should never happen!;
__CPU call vtutor%error('am05_labertw: iteration limit reached. Should never happen: execution aborted.')
      return
      END SUBROUTINE AM05_LABERTW

!aem AM05 included

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE AM05_LABERTW_SIMD(Z, RESULT, MASK0)
      USE prec
      USE simd
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: Z
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RESULT
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: TEMP_1, T, P
      INTEGER                            :: I, II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: E = EXP(1.0_q), E_INV = 1.0_q / E
      
!     ** IF Z TOO LOW, GO WITH THE FIRST TERM OF THE POWER EXPANSION, Z
      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(TEMP_1)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = Z%X(II)
         IF (MASK0%X(II)) THEN
            IF (TEMP_1 .GE. 1.0E-20_q) THEN
               MASK1%X(II) = .TRUE.
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
               RESULT%X(II) = TEMP_1
            ENDIF
         ENDIF
      ENDDO
      
      IF (.NOT. TRUE_FOR_ANY) RETURN

!     $OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1)
      DO II = 0, SIMD_WIDTH_64BIT - 1
!     ** INITAL GUESS
         IF (ABS(Z%X(II) + E_INV) .GT. 1.45_q) THEN
!        ** ASYMPTOTIC EXPANSION AT 0 AND INF
            TEMP_1 = SIMD_LOG(Z%X(II))
            TEMP_1 = TEMP_1 - SIMD_LOG(TEMP_1)
         ELSE
!        ** SERIES EXPANSION ABOUT -1/E TO FIRST ORDER
            TEMP_1 = 1.0_q * SQRT(2.0_q * E * Z%X(II) + 2.0_q) - 1.0_q
         ENDIF
         IF (MASK1%X(II)) THEN
            RESULT%X(II) = TEMP_1
         ENDIF
      ENDDO

!     ** FIND RESULT THROUGH ITERATION
      DO I = 1, 10
         TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(P, T)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            P = SIMD_EXP(RESULT%X(II))
            T = RESULT%X(II) * P - Z%X(II)
            IF (RESULT%X(II) .NE. -1.0_q ) THEN
               T =  T / (P * (RESULT%X(II) + 1.0_q) - 0.5_q * (RESULT%X(II) + 2.0_q) * T / (RESULT%X(II) + 1.0_q))
            ELSE
               T = 0.0_q
            ENDIF
            IF (MASK1%X(II)) THEN
               RESULT%X(II) = RESULT%X(II) - T
               IF (ABS(T) .GE. (2.48_q * 1.0E-14_q) * (1.0_q + ABS(RESULT%X(II)))) THEN
                  TRUE_FOR_ANY = .TRUE.
               ELSE
                  MASK1%X(II) = .FALSE.
               ENDIF
            ENDIF
         ENDDO
         IF (.NOT. TRUE_FOR_ANY) RETURN
      ENDDO
!     ** THIS SHOULD NEVER HAPPEN!;
      CALL vtutor%error("AM05_LABERTW_SIMD: ITERATION LIMIT REACHED. \n SHOULD NEVER HAPPEN: EXECUTION &
         &ABORTED.")
      RETURN
      END SUBROUTINE AM05_LABERTW_SIMD
#endif

!************************ SUBROUTINE GGAXCS *****************************
!
!>  calculates the PW 91 xc-functional for spin-dependent potentials
!>  using the algorithm of J.A. White and D.M. Bird
!>    (Phys.Rev.B 50,7 (1994) 4954)
!>  presently finite differences are used to calculate
!>  the required derivatives
!
!***********************************************************************

      SUBROUTINE GGAXCS(ro1,ro2,q1,q2,exc, &
     &   excd1,excd2,excq1,excq2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(GGAXCS)
      REAL(q),PARAMETER :: PI=3.1415926536_q
      REAL(q),PARAMETER :: THRD=0.333333333333333_q
GPACC !$ACC ROUTINE (FUN) SEQ
GPACC !$ACC ROUTINE (DFRO1) SEQ
GPACC !$ACC ROUTINE (DFRO2) SEQ
GPACC !$ACC ROUTINE (DFQ1) SEQ
GPACC !$ACC ROUTINE (DFQ2) SEQ

      ro   = ro1 + ro2
      exc  = fun(ro1,ro2,q1,q2)
!----------------------------------------------------------------------
!  The derivatives of f_xc necessary to
!  construct spin up and down potentials are
!  computed numerically
!----------------------------------------------------------------------
      excd1 = dfro1(ro1,ro2,q1,q2)
      excd2 = dfro2(ro1,ro2,q1,q2)
      excq1 = dfq1 (ro1,ro2,q1,q2)
      excq2 = dfq2 (ro1,ro2,q1,q2)

      RETURN
      END SUBROUTINE GGAXCS


!***********************************************************************
!     f_xc PW91 construction
!***********************************************************************
      function fun(ro1,ro2,q1,q2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(FUN)
GPACC !$ACC ROUTINE (EXPW) SEQ
GPACC !$ACC ROUTINE (CPW) SEQ

      ro   = ro1 + ro2
      qq   = q1  + q2
      xi   = (ro1 - ro2)/ro

      call expw(ro1,q1,ex1)
      call expw(ro2,q2,ex2)
      call cpw(xi,ro,qq,ec)
      despin=0

      fx1  = ro1*ex1
      fx2  = ro2*ex2
      fcg  = ro*ec
      fsg  = ro*despin
      fun  = fx1 + fx2 + fcg + fsg

      return
      end function fun

      function dfq1(ro1,ro2,q1,q2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(DFQ1)
      parameter (delta=1E-2_q)
GPACC !$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*q1))+1E-15_q
      x1  = q1 - eps
      x3  = q1 + eps
      f1  = fun(ro1,ro2,x1,q2)
      f3  = fun(ro1,ro2,x3,q2)
      dfq1 = (f3-f1) / (2*eps)
      return
      end function dfq1

      function dfq2(ro1,ro2,q1,q2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(DFQ2)
      parameter (delta=1E-2_q)
GPACC !$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*q2))+1E-15_q
      x1  = q2 - eps
      x3  = q2 + eps
      f1  = fun(ro1,ro2,q1,x1)
      f3  = fun(ro1,ro2,q1,x3)
      dfq2 = (f3-f1) / (2*eps)
      return
      end function dfq2

      function dfro1(ro1,ro2,q1,q2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(DFRO1)
      parameter (delta=1E-4_q)
GPACC !$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*ro1))+1E-15_q
      x1  = ro1 - eps
      x3  = ro1 + eps
      f1  = fun(x1,ro2,q1,q2)
      f3  = fun(x3,ro2,q1,q2)
      dfro1 = (f3-f1) / (2*eps)
      return
      end function dfro1

      function dfro2(ro1,ro2,q1,q2)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(DFRO2)
      parameter (delta=1E-4_q)
GPACC !$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*ro2))+1E-15_q
      x1  = ro2 - eps
      x3  = ro2 + eps
      f1  = fun(ro1,x1,q1,q2)
      f3  = fun(ro1,x3,q1,q2)
      dfro2 = (f3-f1) / (2*eps)
      return
      end function dfro2

!***********************************************************************
!>     GGA91 EXCHANGE
!>  gga91 exchange for a spin-unpolarized electronic system
!>  input d : density
!>  input s:  abs(grad d)/(2*kf*d)
!>  output:  exchange energy per electron (ex) and its derivatives
!>           w.r.t. d (exd) and dd (exdd)
!***********************************************************************

      subroutine expw(rhos,rhops,ex)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EXPW)

      REAL(q),PARAMETER :: A1=0.19645_q,A2=0.27430_q,A3=0.15084_q,A4=100.0_q
      REAL(q),PARAMETER :: AX=-0.7385588_q,A=7.7956_q,B1=0.004_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,THRD4=1.333333333333333_q
      REAL(q),PARAMETER :: THPITH=3.0936677262801_q

      rho=2._q*rhos
      rhop=2._q*rhops
      FAC = AX*rho**THRD

      w=0.16162046_q
      s=w*abs(rhop)/rho**(4._q/3._q)
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2

      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.0_q/SQRT(1.0_q+A*A*S2)
      P1 = LOG(A*S+1.0_q/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.0_q/(1.0_q+A1*S*P1+B1*S4)
      P4 = 1.0_q+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*(F-1)
      RETURN
      END SUBROUTINE EXPW

!***********************************************************************
!>     GGA91 CORRELATION
!***********************************************************************
      subroutine cpw(xi,ro,qq,ecc)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE ldalib, only : CORLSD
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CPW)
GPACC !$ACC ROUTINE(CORLSD) SEQ

      REAL(q),PARAMETER :: C1=0.002568_q,C2=0.023266_q,C3=7.389E-6_q,C4=8.723_q
      REAL(q),PARAMETER :: C5=0.472_q,C6=7.389E-2_q,A4=100.0_q
      REAL(q),PARAMETER :: XNU=15.75592_q,CC0=0.004235_q,CX=-0.001667212_q,ALF=0.09_q
      REAL(q),PARAMETER :: THRDM=-0.333333333333_q,THRD2=0.666666666667_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,SIXTH7=1.16666666666667_q
      REAL(q),PARAMETER :: PI=3.1415926536_q
!
!     LOCAL CORRELATION is needed to compute EC:
!
      rs   = (4._q*pi*ro/3._q)**(-1._q/3._q)
      zet   = xi
      CALL CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!
!     NONLOCAL CORRELATION:
!
      FK = 1.91915829_q/RS
      SK = SQRT(4.0_q*FK/PI)
      G = ((1.0_q+ZET)**THRD2+(1.0_q-ZET)**THRD2)/2.0_q
      t =  abs(qq)/(2.0_q*G*sk*ro)
!
      BET = XNU*CC0
      DELT = 2.0_q*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.0_q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.0_q+B*T2
      Q5 = 1.0_q+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.0_q+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.0_q*CX/7.0_q
      R2 = XNU*COEFF*G3
      R3 = EXP(-R1*T2)
      H0 = G3*(BET/DELT)*LOG(1.0_q+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
!
      ECC=H
!
      RETURN
      END SUBROUTINE CPW

!***********************************************************************
!
! B3LYP exchange and correlation energy
! Note: for the Becke part only the gradient corrections are calculated
! whereas the correlation part includes the local (Wigner) term
!
!***********************************************************************

      FUNCTION EXCB3LYP(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info

      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(EXCB3LYP)

      TYPE (xc_info) XC

      REAL(q) EXCB3LYP

      REAL(q) :: RHO, RHO1,RHO2,DRHO1,DRHO2, DRHOA
      REAL(q) :: EXC, TEMP

      REAL(q) :: AX, AC
GPACC !$ACC ROUTINE (EXB88) SEQ
GPACC !$ACC ROUTINE (ECLYP) SEQ

      AX=XC%AGGAX
      AC=XC%AGGAC

      EXCB3LYP = 0.0_q
      EXC=0._q
      RHO  = RHO1+RHO2

! Becke88 gradient corrections up electrons
      CALL EXB88(RHO1,DRHO1,TEMP)
      EXC = EXC + AX *TEMP
! Becke88 gradient corrections down electrons
      CALL EXB88(RHO2,DRHO2,TEMP)
      EXC = EXC + AX *TEMP

! LYP correlation
      CALL ECLYP(RHO1,RHO2,DRHO1,DRHO2,DRHOA,TEMP)
      EXC = EXC + AC*TEMP
      EXCB3LYP=EXC

      END FUNCTION EXCB3LYP

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCB3LYP_SIMD(XC, RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE setexm_struct_def
      USE simd
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: TEMP, TEMP_RES
      INTEGER                            :: II

      CALL EXB88_SIMD(RHO1, DRHO1, TEMP, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) 
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_RES%X(II) = XC%AGGAX * TEMP%X(II)
      ENDDO

      CALL EXB88_SIMD(RHO2, DRHO2, TEMP, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) 
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_RES%X(II) = TEMP_RES%X(II) +  XC%AGGAX * TEMP%X(II)
      ENDDO

      CALL ECLYP_SIMD(RHO1, RHO2, DRHO1, DRHO2, DRHOA, TEMP, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) 
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = TEMP_RES%X(II) + XC%AGGAC * TEMP%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCB3LYP_SIMD
#endif

!***********************************************************************
!
!> B3LYP exchange and correlation energy
!> Note: for the Becke part only the gradient corrections are calculated
!> whereas the correlation part includes the local (Wigner) term
!>
!> this subroutine calculates the derivative w.r.t. each argument
!> of the B3LYP functional using finite differences
!
!***********************************************************************
     SUBROUTINE B3LYPXCS(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA,EXC,EXCD1,EXCDD1, &
                         EXCD2,EXCDD2,EXCDDA)
GPACC !$ACC ROUTINE SEQ
     USE prec
     USE constant
     USE setexm_struct_def, ONLY : xc_info
     IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(B3LYPXCS)
     TYPE (xc_info) XC
     REAL(q) RHO1,RHO2,DRHO1,DRHO2,DRHOA
     REAL(q) EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA

GPACC !$ACC ROUTINE (EXCB3LYP) SEQ
GPACC !$ACC ROUTINE (DB3LYP_RHO1) SEQ
GPACC !$ACC ROUTINE (DB3LYP_RHO2) SEQ
GPACC !$ACC ROUTINE (DB3LYP_DRHO1) SEQ
GPACC !$ACC ROUTINE (DB3LYP_DRHO2) SEQ
GPACC !$ACC ROUTINE (DB3LYP_DRHOA) SEQ
!    B3LYP - Recipe:
!    EXC = EXC(LSDA) + a0 * (EX(HF) - EX(LSDA)) + aX * DELTA EX(B88) + aC * DELTA EC(LYP)
!
!    a0 = 0.20 ( == AEXX !)  : 0.2  * EX(LSDA)
!    aX = 0.72               : 0.72 * DELTA EX(B88)
!    aC = 0.81               : 0.81 * DELTA EC(LYP)
!

     IF (RHO1<0._q .OR. RHO2<0._q ) THEN
       EXC   = 0._q
       EXCD1  = 0._q
       EXCDD1 = 0._q
       EXCD2  = 0._q
       EXCDD2 = 0._q
       EXCDDA = 0._q
       RETURN
     ENDIF

! Exchange-correlation energy

     EXC=EXCB3LYP(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)

! Exchange correlation potential

     EXCD1=DB3LYP_RHO1(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCD2=DB3LYP_RHO2(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCDD1=DB3LYP_DRHO1(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCDD2=DB3LYP_DRHO2(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCDDA=DB3LYP_DRHOA(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)

!     WRITE(91,'(6F12.7)') RHO1, RHO2, EXCD1, EXCD2, EXCDD1, EXCDD2

     END SUBROUTINE B3LYPXCS

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE B3LYPXCS_SIMD(XC, RHO1_IN, RHO2_IN, DRHO1, DRHO2, DRHOA, EXC, EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA, MASK0)
      USE prec
      USE constant
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1_IN, RHO2_IN, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXC, EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: RHO1, RHO2
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
!    B3LYP - Recipe:
!    EXC = EXC(LSDA) + a0 * (EX(HF) - EX(LSDA)) + aX * DELTA EX(B88) + aC * DELTA EC(LYP)
!
!    a0 = 0.20 ( == AEXX !)  : 0.2  * EX(LSDA)
!    aX = 0.72               : 0.72 * DELTA EX(B88)
!    aC = 0.81               : 0.81 * DELTA EC(LYP)
!
      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF ((RHO1_IN%X(II) .GE. 0.0_q) .AND. (RHO2_IN%X(II) .GE. 0.0_q)) THEN
               TRUE_FOR_ANY = .TRUE.
               RHO1%X(II) = RHO1_IN%X(II)
               RHO2%X(II) = RHO2_IN%X(II)
            ELSE
               MASK1%X(II) = .FALSE.
               RHO1%X(II) = 0.0_q
               RHO2%X(II) = 0.0_q
               EXC%X(II) = 0.0_q
               EXCD1%X(II) = 0.0_q
               EXCDD1%X(II) = 0.0_q
               EXCD2%X(II) = 0.0_q
               EXCDD2%X(II) = 0.0_q
               EXCDDA%X(II) = 0.0_q
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

      CALL EXCB3LYP_SIMD(XC, EXC, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      CALL DB3LYP_RHO1_SIMD(XC, EXCD1, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)
      
      CALL DB3LYP_RHO2_SIMD(XC, EXCD2, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      CALL DB3LYP_DRHO1_SIMD(XC, EXCDD1, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      CALL DB3LYP_DRHO2_SIMD(XC, EXCDD2, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)
 
      CALL DB3LYP_DRHOA_SIMD(XC, EXCDDA, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      RETURN
      END SUBROUTINE B3LYPXCS_SIMD
#endif

!**********************************************************************
!
!> Becke 88 gradient correction to LDA exchange for a spin-unpolarized
!> electronic system
!> A.D. Becke, Phys. Rev. A, 38, 3098 (1988).
!
!**********************************************************************

      SUBROUTINE EXB88(RHOS,DRHOS,EX)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EXB88)

      REAL(q), PARAMETER :: BETA = 0.0042_q

      RHO43 = RHOS**(4.0_q/3.0_q)
      X = ABS(DRHOS)/RHO43
      X2 = X*X
      B6 = 6.0_q * BETA

      ! arsinh(x) = LOG(  x + SQRT( x * x + 1.0_q))
      DENOM = 1.0_q + B6*X*LOG(X + SQRT(X2 + 1.0_q))
      IF (RHOS .GE. 0.00001_q ) THEN
      EX = -BETA * RHO43 * X2/DENOM
      ELSE
      EX = 0.0_q
      END IF
      RETURN
      END SUBROUTINE EXB88

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXB88_SIMD(RHOS, DRHOS, EX, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHOS, DRHOS
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EX
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: BETA = -0.0042_q
      REAL(q), PARAMETER                 :: THRD4 = -4.0_q / 3.0_q 

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, RHO43, X, X2, B6, DENOM)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = RHOS%X(II)
         TEMP_2 = DRHOS%X(II)
         TEMP_3 = EX%X(II)
         TEMP_3 = 0.0_q
         IF (TEMP_1 .GE. 0.00001_q) THEN
            RHO43 = SIMD_POW(TEMP_1, THRD4)
            X = ABS(TEMP_2) * RHO43
            X2 = X * X
            B6 = 6.0_q * BETA
            ! arsinh(x) = LOG(x + SQRT(x * x + 1.0_q))
            TEMP_2 = X + SQRT(1.0_q + X2)
            DENOM = 1.0_q - B6 * X * SIMD_LOG(TEMP_2)
            TEMP_3 = BETA * X2 / (RHO43 * DENOM)
         ENDIF
         IF (MASK0%X(II)) THEN
            EX%X(II) = TEMP_3
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXB88_SIMD
#endif

!**********************************************************************
!
!> The correlation energy density of Lee, Yang and Parr
!> B. Miehlich, A. Savin, H. Stoll and H. Preuss, Chem. Phys. Lett. 157,
!> 200 (1989).
!> C. Lee, W. Yang, and R. Parr, Phys. Rev. B 37, 785 (1988)
!
!**********************************************************************

      SUBROUTINE ECLYP(RHOU,RHOD,DRHOU,DRHOD,DRHOA,ELYP)
GPACC !$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(ECLYP)

      REAL(q), PARAMETER :: CF=2.87123400018819152526_q
      REAL(q), PARAMETER :: ALYP=0.04918_q
      REAL(q), PARAMETER :: BLYP=0.132_q
      REAL(q), PARAMETER :: CLYP=0.2533_q
      REAL(q), PARAMETER :: DLYP=0.349_q

      RHO=RHOU+RHOD
      DRHO=DRHOA
      DRHOSQ=ABS(DRHO)*ABS(DRHO)
      DRHOUSQ=ABS(DRHOU)*ABS(DRHOU)
      DRHODSQ=ABS(DRHOD)*ABS(DRHOD)

      X =  CLYP* ( RHO**(-1.0_q/3.0_q) )

      Y =  DLYP* ( RHO**(-1.0_q/3.0_q) )

      Z =  1.0_q + DLYP * ( RHO**(-1.0_q/3.0_q) )

      OMEGA = EXP(-X)/Z * ( RHO**(-11.0_q/3.0_q) )

      DELTA = X + Y/Z

      SUMMAND_A=0
      SUMMAND_B=0
      SUMMAND_C=0
      SUMMAND_D=0
      SUMMAND_E=0

      IF (RHOU .GE. 0.00001_q .AND. RHOD .GE. 0.00001_q ) THEN
      ! Wigner term
      SUMMAND_A = -ALYP* 4._q/( 1._q + DLYP*RHO**(-1._q/3._q) ) * RHOU*RHOD/RHO
      ! second local term
      SUMMAND_B = 2._q**(11._q/3._q)*CF*( RHOU**(8._q/3._q) + RHOD**(8._q/3._q))
      ! non local terms
      SUMMAND_C = ( 47._q/18._q - 7._q/18._q * DELTA) * DRHOSQ - (5._q/2._q - 1._q/18._q * DELTA) * ( DRHOUSQ + DRHODSQ)

      SUMMAND_D = -(DELTA - 11._q)/9._q * ( (RHOU/RHO) * DRHOUSQ + (RHOD/RHO) * DRHODSQ)

      SUMMAND_E = (-2._q/3._q)*RHO*RHO*DRHOSQ + ( (2._q/3._q) *RHO*RHO - RHOU*RHOU) * DRHODSQ + &
                  ((2._q/3._q) *RHO*RHO - RHOD*RHOD) * DRHOUSQ

      ! "local" contribution (i.e Wigner term + term "B")
      ELYP =  SUMMAND_A - ALYP*BLYP*OMEGA* ( RHOU*RHOD * ( SUMMAND_B ) )

      ! "non-local" contribution
      ELYP = ELYP - ALYP*BLYP*OMEGA* ( RHOU*RHOD * ( SUMMAND_C + SUMMAND_D ) +SUMMAND_E )
      ELSE
      ELYP = 0.0_q
      END IF
      RETURN
      END SUBROUTINE ECLYP

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE ECLYP_SIMD(RHOU, RHOD, DRHOU, DRHOD, DRHOA, ELYP, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHOU, RHOD, DRHOU, DRHOD, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: ELYP
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: OMEGA, DELTA
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: CF = 2.87123400018819152526_q
      REAL(q), PARAMETER                 :: ALYP = 0.04918_q
      REAL(q), PARAMETER                 :: BLYP = 0.132_q
      REAL(q), PARAMETER                 :: CLYP = 0.2533_q
      REAL(q), PARAMETER                 :: DLYP = 0.349_q
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q
      REAL(q), PARAMETER                 :: MTHRD = -1.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P1 = -11.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P2 = 2.0_q ** (11.0_q/3.0_q)
      REAL(q), PARAMETER                 :: P3 = 8.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P4 = 47.0_q / 18.0_q
      REAL(q), PARAMETER                 :: P5 = 7.0_q / 18.0_q
      REAL(q), PARAMETER                 :: P6 = 1.0_q / 18.0_q
      REAL(q), PARAMETER                 :: P7 = 2.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P8 = 1.0_q / 9.0_q

      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(TEMP_1, TEMP_2, RHO, RHO_MTHRD, X, Y, Z, RHO_INV, SUMMAND_A, SUMMAND_B, TEMP_3)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = RHOU%X(II)
         TEMP_2 = RHOD%X(II)
         RHO = TEMP_1 + TEMP_2
         RHO_MTHRD = SIMD_POW(RHO, MTHRD)
         X = CLYP * RHO_MTHRD
         Y = DLYP * RHO_MTHRD
         Z = 1.0_q / (1.0_q + DLYP * RHO_MTHRD)
         OMEGA%X(II) = SIMD_EXP(-X) * Z * SIMD_POW(RHO, P1)
         DELTA%X(II) = X + Y * Z
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF ((TEMP_1 .GE. 0.00001_q) .AND. (TEMP_2 .GE. 0.00001_q)) THEN
               TRUE_FOR_ANY = .TRUE.
               RHO_INV = 1.0_q / RHO
               ! Wigner term
               SUMMAND_A = -ALYP * 4.0_q / (1.0_q + DLYP * RHO_MTHRD) * TEMP_1 * TEMP_2 * RHO_INV
               ! second local term
               SUMMAND_B = P2 * CF * (SIMD_POW(TEMP_1, P3) + SIMD_POW(TEMP_2, P3))
               ! "local" contribution (i.e Wigner term + term "B")
               TEMP_3 = SUMMAND_A - ALYP * BLYP * OMEGA%X(II) * (TEMP_1 * TEMP_2 * SUMMAND_B)
            ELSE
               MASK1%X(II) = .FALSE.
               TEMP_3 = 0.0_q
            ENDIF
            ELYP%X(II) = TEMP_3
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, RHO, DRHO, DRHOSQ, DRHOUSQ, DRHODSQ, RHO_INV, SUMMAND_C, SUMMAND_D, SUMMAND_E)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = RHOU%X(II)
         TEMP_2 = RHOD%X(II)
         RHO = TEMP_1 + TEMP_2
         DRHO = DRHOA%X(II)
         DRHOSQ = DRHO * DRHO
         DRHOUSQ = DRHOU%X(II) * DRHOU%X(II)
         DRHODSQ = DRHOD%X(II) * DRHOD%X(II)
         IF (MASK1%X(II)) THEN
            RHO_INV = 1.0_q / RHO
            ! non local terms
            SUMMAND_C = (P4 - P5 * DELTA%X(II)) * DRHOSQ - (2.5_q - P6 * DELTA%X(II)) * (DRHOUSQ + DRHODSQ)
            SUMMAND_D = -(DELTA%X(II) - 11.0_q) * P8 * (TEMP_1 * RHO_INV * DRHOUSQ + TEMP_2 * RHO_INV * DRHODSQ)              
            DELTA%X(II) = P7 * RHO * RHO
            SUMMAND_E = -DELTA%X(II) * DRHOSQ + (DELTA%X(II) - TEMP_1 * TEMP_1) * DRHODSQ + &
                        (DELTA%X(II) - TEMP_2 * TEMP_2) * DRHOUSQ
            ! "non-local" contribution
            ELYP%X(II) = ELYP%X(II) - ALYP * BLYP * OMEGA%X(II) * (TEMP_1 * TEMP_2 * (SUMMAND_C + SUMMAND_D) + SUMMAND_E)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE ECLYP_SIMD
#endif

!**********************************************************************
!
!> the following functions calculate the first derivative of
!> the semi-local part of the B3LYP functional with respect to
!> each argument
!
!**********************************************************************
      FUNCTION DB3LYP_DRHO1(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(DB3LYP_DRHO1)
      TYPE (xc_info) XC
      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_DRHO1
GPACC !$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*DRHO1))+1E-15_q
      X1  = DRHO1 - EPS
      X3  = DRHO1 + EPS
      F1  = EXCB3LYP(XC,RHO1,RHO2,X1,DRHO2,DRHOA)
      F3  = EXCB3LYP(XC,RHO1,RHO2,X3,DRHO2,DRHOA)
      DB3LYP_DRHO1 = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_DRHO1

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_DRHO1_SIMD(XC, RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * DRHO1%X(II))) + 1.0E-15_q
         X1%X(II) = DRHO1%X(II) - EPS%X(II)
         X3%X(II) = DRHO1%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(XC, F1, RHO1, RHO2, X1, DRHO2, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(XC, F3, RHO1, RHO2, X3, DRHO2, DRHOA, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_DRHO1_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_DRHO2(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(DB3LYP_DRHO2)

      TYPE (xc_info) XC
      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_DRHO2
GPACC !$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*DRHO2))+1E-15_q
      X1  = DRHO2 - EPS
      X3  = DRHO2 + EPS
      F1  = EXCB3LYP(XC,RHO1,RHO2,DRHO1,X1,DRHOA)
      F3  = EXCB3LYP(XC,RHO1,RHO2,DRHO1,X3,DRHOA)
      DB3LYP_DRHO2 = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_DRHO2

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_DRHO2_SIMD(XC, RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * DRHO2%X(II))) + 1.0E-15_q
         X1%X(II) = DRHO2%X(II) - EPS%X(II)
         X3%X(II) = DRHO2%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(XC, F1, RHO1, RHO2, DRHO1, X1, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(XC, F3, RHO1, RHO2, DRHO1, X3, DRHOA, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_DRHO2_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_DRHOA(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(DB3LYP_DRHOA)

      TYPE (xc_info) XC
      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_DRHOA
GPACC !$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*DRHOA))+1E-15_q
      X1  = DRHOA - EPS
      X3  = DRHOA + EPS
      F1  = EXCB3LYP(XC,RHO1,RHO2,DRHO1,DRHO2,X1)
      F3  = EXCB3LYP(XC,RHO1,RHO2,DRHO1,DRHO2,X3)
      DB3LYP_DRHOA = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_DRHOA

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_DRHOA_SIMD(XC, RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * DRHOA%X(II))) + 1.0E-15_q
         X1%X(II) = DRHOA%X(II) - EPS%X(II)
         X3%X(II) = DRHOA%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(XC, F1, RHO1, RHO2, DRHO1, DRHO2, X1, MASK0)

      CALL EXCB3LYP_SIMD(XC, F3, RHO1, RHO2, DRHO1, DRHO2, X3, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_DRHOA_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_RHO1(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(DB3LYP_RHO1)

      TYPE (xc_info) XC
      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_RHO1
GPACC !$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*RHO1))+1E-15_q
      X1  = RHO1 - EPS
      X3  = RHO1 + EPS
      F1  = EXCB3LYP(XC,X1,RHO2,DRHO1,DRHO2,DRHOA)
      F3  = EXCB3LYP(XC,X3,RHO2,DRHO1,DRHO2,DRHOA)
      DB3LYP_RHO1 = (F3-F1) / (2*EPS)
!      WRITE(82,*) RHO1, RHO2, DB3LYP_RHO1
      RETURN
      END FUNCTION DB3LYP_RHO1

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_RHO1_SIMD(XC, RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * RHO1%X(II))) + 1.0E-15_q
         X1%X(II) = RHO1%X(II) - EPS%X(II)
         X3%X(II) = RHO1%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(XC, F1, X1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(XC, F3, X3, RHO2, DRHO1, DRHO2, DRHOA, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_RHO1_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_RHO2(XC,RHO1,RHO2,DRHO1,DRHO2,DRHOA)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(DB3LYP_RHO2)

      TYPE (xc_info) XC
      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_RHO2
GPACC !$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*RHO2))+1E-15_q
      X1  = RHO2 - EPS
      X3  = RHO2 + EPS
      F1  = EXCB3LYP(XC,RHO1,X1,DRHO1,DRHO2,DRHOA)
      F3  = EXCB3LYP(XC,RHO1,X3,DRHO1,DRHO2,DRHOA)
      DB3LYP_RHO2 = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_RHO2

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_RHO2_SIMD(XC, RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * RHO2%X(II))) + 1.0E-15_q
         X1%X(II) = RHO2%X(II) - EPS%X(II)
         X3%X(II) = RHO2%X(II) + EPS%X(II)
      ENDDO
      
      CALL EXCB3LYP_SIMD(XC, F1, RHO1, X1, DRHO1, DRHO2, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(XC, F3, RHO1, X3, DRHO1, DRHO2, DRHOA, MASK0)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
           RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_RHO2_SIMD
#endif

      SUBROUTINE CALC_EXCHWPBEsol_SP(XC,D,S, &
     &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
     &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT D : DENSITY
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(6 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EXWPBE)
!  OUTPUT:  DERIVATIVE W.R.T. DENSITY (EXWPBED)
!  OUTPUT:  DERIVATIVE W.R.T. GRADIENT DENSITY (EXWPBEDD)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE wpbe
      USE setexm_struct_def, ONLY: xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(CALC_EXCHWPBEsol_SP)

      TYPE (xc_info) XC

      REAL(q) F,FS,FSR,FLR
      REAL(q) D,LOGRHO,S,RS,SW

! Numerical derivatives test
      REAL(q) XKF,XKFMIN,XKFPLS,DD
      REAL(q) DDMIN,DDPLS
      REAL(q) SMIN,SPLS,DFDR,DFDDD
      REAL(q) DMIN,DPLS,FSR_MIN,FSR_PLS,XFKMIN,XFKPLS

      REAL(q) EXWPBE,EXWPBE_SR,EXWPBE_LR
      REAL(q) EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR
      REAL(q) EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR
      REAL(q) EXPBE,EXPBED,EXPBEDD
      REAL(q) LDASCREEN

! Parameters to calculate the PBE_SOL enhancement factor
      REAL(q), PARAMETER :: AX=-0.738558766382022405884230032680836_q
      REAL(q), PARAMETER :: UM=0.12345679012345679012_q  !LS
      REAL(q), PARAMETER :: UK=0.8040_q
      REAL(q), PARAMETER :: UL=UM/UK
GPACC !$ACC ROUTINE(WPBE_SPLINE) SEQ

      LDASCREEN=XC%LDASCREEN

! set up the table for F_SR(s)
__CPU CALL INIT_WPBE_TABLE(XC)

      SW=ABS(S)
      XKF=(3._q*PI*PI*D)**(1._q/3._q)

! Due to the fact that most of the times only the GGA corrections to the
! LDA exchange are wanted we will have to recalculate the LDA contributions
! in order to be able to subtract them from our PBE result.
! (Georg thinks this makes a lot of sense, but I guess there are not many
! people who would agree with this!)
      RS=(3._q/(4._q*PI)/D)**(1/3._q)
      EXLDA=EX(RS,1,.FALSE.)*D/2._q
      EXLDA_SR=EX_SR(RS,LDASCREEN*AUTOA,1)*D/2._q
      EXLDA_LR=EXLDA-EXLDA_SR
      EXDLDA=VX(RS,1,.FALSE.)/2._q
      EXDLDA_SR=VX_SR(RS,LDASCREEN*AUTOA,1)/2._q
      EXDLDA_LR=VX(RS,1,.FALSE.)/2._q-EXDLDA_SR

! PBE quantities
      F=1._q+UK-UK/(1._q+UL*SW*SW)
      FS=2._q*UM/(1._q+UL*SW*SW)/(1._q+UL*SW*SW)
      EXPBE=F*EXLDA
      EXPBED=4._q/3._q*(F-SW*SW*FS)*EXLDA/D
      EXPBEDD=0.5_q*AX*0.3232409194_q*SW*FS
! set some stuff to zero
      EXWPBE_SR=0
      EXWPBE_LR=0
      EXWPBED_SR=0
      EXWPBED_LR=0
      EXWPBEDD_SR=0
      EXWPBEDD_LR=0
! Interpolate the short range part of the enhancement
! coefficient from the tables
      IF (SW==0._q.OR.D<0._q) RETURN

      LOGRHO=LOG(D)
!    Cutoff criterion to enforce local Lieb-Oxford
!    this ensures that the enhancement factor does not exceed the
!    original one

      DD=SW*(2._q*XKF*D)
!LS
      CALL WPBE_SPLINE(LOGRHO,SW,FSR)
! get values from function directly (slower)
!      CALL EXCHWPBE_R(LDASCREEN*AUTOA,D,SW,FSR)
! Numerical derivatives with interpolated values:
! Firstly: derivatives w.r.t rho
      DMIN=D-0.0001_q*D
      LOGRHO=LOG(DMIN)
      XKFMIN=(3._q*PI*PI*DMIN)**(1._q/3._q)
      SMIN=DD/(DMIN*XKFMIN*2._q)
      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)

      DPLS=D+0.0001_q*D
      LOGRHO=LOG(DPLS)
      XKFPLS=(3._q*PI*PI*DPLS)**(1._q/3._q)
      SPLS=DD/(DPLS*XKFPLS*2._q)

      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDR=(FSR_PLS-FSR_MIN)/2._q/D/0.0001_q
! Secondly: derivatives w.r.t gradient(rho)
      LOGRHO=LOG(D)
      DDMIN=DD-0.0001_q*DD
      DDPLS=DD+0.0001_q*DD
      SMIN=DDMIN/(D*XKF*2._q)
      SPLS=DDPLS/(D*XKF*2._q)

      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)
      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDDD=(FSR_PLS-FSR_MIN)/2._q/DD/0.0001_q

! Find the complementary long range part
      FLR=F-FSR
! E_sr = fsr*Exlda
      EXWPBE_SR=FSR*EXLDA
      EXWPBE_LR=FLR*EXLDA
! dE_sr/drho = dfsr/drho*Exlda + fsr*dExlda/drho
      EXWPBED_SR=DFDR*EXLDA+4._q/3._q*FSR*EXLDA/D
      EXWPBED_LR=EXPBED-EXWPBED_SR
! dE_sr/d(grad rho) = dfsr/dgrad(rho)*Exlda
!gK correct the sign
      EXWPBEDD_SR=DFDDD*EXLDA *SIGN(1.0_q,S)
      EXWPBEDD_LR=(EXPBEDD-DFDDD*EXLDA) *SIGN(1.0_q,S)
!gK end
      RETURN
      END SUBROUTINE CALC_EXCHWPBEsol_SP

      SUBROUTINE CALC_EXCHWPBE_SP(XC,D,S, &
     &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
     &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT D : DENSITY
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(6 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EXWPBE)
!  OUTPUT:  DERIVATIVE W.R.T. DENSITY (EXWPBED)
!  OUTPUT:  DERIVATIVE W.R.T. GRADIENT DENSITY (EXWPBEDD)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE wpbe
      USE setexm_struct_def, ONLY: xc_info
      IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(CALC_EXCHWPBE_SP)

      TYPE (xc_info) XC

      REAL(q) F,FS,FSR,FLR
      REAL(q) D,LOGRHO,S,RS,SW

! Numerical derivatives test
      REAL(q) XKF,XKFMIN,XKFPLS,DD
      REAL(q) DDMIN,DDPLS
      REAL(q) SMIN,SPLS,DFDR,DFDDD
      REAL(q) DMIN,DPLS,FSR_MIN,FSR_PLS,XFKMIN,XFKPLS

      REAL(q) EXWPBE,EXWPBE_SR,EXWPBE_LR
      REAL(q) EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR
      REAL(q) EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR
      REAL(q) EXPBE,EXPBED,EXPBEDD
      REAL(q) LDASCREEN

! Parameters to calculate the PBE enhancement factor
      REAL(q), PARAMETER :: AX=-0.738558766382022405884230032680836_q
      REAL(q), PARAMETER :: UM=0.2195149727645171_q
      REAL(q), PARAMETER :: UK=0.8040_q
      REAL(q), PARAMETER :: UL=UM/UK
GPACC !$ACC ROUTINE(WPBE_SPLINE) SEQ

      LDASCREEN=XC%LDASCREEN

! set up the table for F_SR(s)
__CPU CALL INIT_WPBE_TABLE(XC)

      SW=ABS(S)
      XKF=(3._q*PI*PI*D)**(1._q/3._q)

! Due to the fact that most of the times only the GGA corrections to the
! LDA exchange are wanted we will have to recalculate the LDA contributions
! in order to be able to subtract them from our PBE result.
! (Georg thinks this makes a lot of sense, but I guess there are not many
! people who would agree with this!)
      RS=(3._q/(4._q*PI)/D)**(1/3._q)
      EXLDA=EX(RS,1,.FALSE.)*D/2._q
      EXLDA_SR=EX_SR(RS,LDASCREEN*AUTOA,1)*D/2._q
      EXLDA_LR=EXLDA-EXLDA_SR
      EXDLDA=VX(RS,1,.FALSE.)/2._q
      EXDLDA_SR=VX_SR(RS,LDASCREEN*AUTOA,1)/2._q
      EXDLDA_LR=VX(RS,1,.FALSE.)/2._q-EXDLDA_SR

! PBE quantities
      F=1._q+UK-UK/(1._q+UL*SW*SW)
      FS=2._q*UM/(1._q+UL*SW*SW)/(1._q+UL*SW*SW)
      EXPBE=F*EXLDA
      EXPBED=4._q/3._q*(F-SW*SW*FS)*EXLDA/D
      EXPBEDD=0.5_q*AX*0.3232409194_q*SW*FS
! set some stuff to zero
      EXWPBE_SR=0
      EXWPBE_LR=0
      EXWPBED_SR=0
      EXWPBED_LR=0
      EXWPBEDD_SR=0
      EXWPBEDD_LR=0
! Interpolate the short range part of the enhancement
! coefficient from the tables
      IF (SW==0._q.OR.D<0._q) RETURN

      LOGRHO=LOG(D)
!    Cutoff criterion to enforce local Lieb-Oxford
!    this ensures that the enhancement factor does not exceed the
!    original one

      DD=SW*(2._q*XKF*D)
      CALL WPBE_SPLINE(LOGRHO,SW,FSR)
! get values from function directly (slower)
!      CALL EXCHWPBE_R(LDASCREEN*AUTOA,D,SW,FSR)
! Numerical derivatives with interpolated values:
! Firstly: derivatives w.r.t rho
      DMIN=D-0.0001_q*D
      LOGRHO=LOG(DMIN)
      XKFMIN=(3._q*PI*PI*DMIN)**(1._q/3._q)
      SMIN=DD/(DMIN*XKFMIN*2._q)
      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)

      DPLS=D+0.0001_q*D
      LOGRHO=LOG(DPLS)
      XKFPLS=(3._q*PI*PI*DPLS)**(1._q/3._q)
      SPLS=DD/(DPLS*XKFPLS*2._q)

      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDR=(FSR_PLS-FSR_MIN)/2._q/D/0.0001_q
! Secondly: derivatives w.r.t gradient(rho)
      LOGRHO=LOG(D)
      DDMIN=DD-0.0001_q*DD
      DDPLS=DD+0.0001_q*DD
      SMIN=DDMIN/(D*XKF*2._q)
      SPLS=DDPLS/(D*XKF*2._q)

      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)
      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDDD=(FSR_PLS-FSR_MIN)/2._q/DD/0.0001_q

! Find the complementary long range part
      FLR=F-FSR
! E_sr = fsr*Exlda
      EXWPBE_SR=FSR*EXLDA
      EXWPBE_LR=FLR*EXLDA
! dE_sr/drho = dfsr/drho*Exlda + fsr*dExlda/drho
      EXWPBED_SR=DFDR*EXLDA+4._q/3._q*FSR*EXLDA/D
      EXWPBED_LR=EXPBED-EXWPBED_SR
! dE_sr/d(grad rho) = dfsr/dgrad(rho)*Exlda
!gK correct the sign
      EXWPBEDD_SR=DFDDD*EXLDA *SIGN(1.0_q,S)
      EXWPBEDD_LR=(EXPBEDD-DFDDD*EXLDA) *SIGN(1.0_q,S)
!gK end
      RETURN
      END SUBROUTINE CALC_EXCHWPBE_SP

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EX_SR_SIMD(RS, RMU, IFLG, EX_SR, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT NONE
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS
      REAL(q), INTENT(IN)                :: RMU
      INTEGER, INTENT(IN)                :: IFLG
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EX_SR
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      REAL(q)                            :: CX, EX_, RHO, QFAC, QF, A, FRAC, Z, T, ERRF_RESULT
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: CX1 = 0.9163305865663_q
      REAL(q), PARAMETER                 :: CX2 = 1.1545041946774_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q    ! 1._q/3._q
      REAL(q), PARAMETER                 :: QFAC1 = 3.093667726280135532590521_q   ! (3._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: QFAC2 = 3.897777089720753540547094_q   ! (6._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: SQRT_PI = 1.772453850905515881919428_q ! sqrt(pi)
      REAL(q), PARAMETER                 :: P1 = 0.620350490899400086597382_q      ! (3._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P_2 = 2.666666666666666518636930_q     ! 8._q/3._q
      REAL(q), PARAMETER                 :: P_3 = 0.027777777777777776235801_q     ! 1._q/36._q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(QFAC, CX, TEMP_1, EX_, QF, A, TEMP_3, Z, T, ERRF_RESULT, TEMP_2, FRAC)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (IFLG .EQ. 1) THEN
            QFAC = QFAC1
            CX = CX1
         ELSE
            QFAC = QFAC2
            CX = CX2
         ENDIF
         TEMP_1 = 1.0_q / RS%X(II)
         EX_ = -CX * TEMP_1
! This version is faster and potentially more precise
         QF = QFAC * P1 * TEMP_1
         A = RMU / (2.0_q * QF)
         TEMP_3 = 0.0_q
! Test on the value of A
         IF (A .LT. 1.E-9_q) THEN
         ! Limit for small A
            TEMP_3 = EX_
         ENDIF
         TEMP_1 = 1.0_q / A
         IF ((A .GE. 1.0E-9_q) .AND. (A .LE. 100.0_q)) THEN
         ! Inline ERRF : potentially faster
            Z = 0.5_q * TEMP_1
            T = 1.0_q / (1.0_q + 0.5_q * Z)
            ERRF_RESULT = -Z * Z - 1.26551223_q + T * (1.00002368_q + T * (.37409196_q +  &
            &  T * (0.09678418_q + T * (-0.18628806_q + T * (0.27886807_q + T * (-1.13520398_q +  &
            &  T * (1.48851587_q + T * (-0.82215223_q + T * 0.17087277_q))))))))
            ERRF_RESULT = 1.0_q - T * SIMD_EXP(ERRF_RESULT)
         ! End of Inline ERRF
            TEMP_2 = 4.0_q * A * A *A
            FRAC = (1.0_q - P_2 * A * (SQRT_PI * ERRF_RESULT +  &
            &  (2.0_q * A - TEMP_2) * SIMD_EXP(-0.25_q * TEMP_1 * TEMP_1) - 3.0_q * A + TEMP_2))
            TEMP_3 = EX_ * FRAC
         ENDIF
         IF ((A .GT. 100.0_q) .AND. (A .LE. 1.0E+9_q)) THEN
         ! Development for large A
            FRAC = P_3 * TEMP_1
            TEMP_3 = EX_ * FRAC
         ENDIF
         IF (MASK0%X(II)) THEN
            EX_SR%X(II) = TEMP_3
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EX_SR_SIMD


      SUBROUTINE VX_SR_SIMD(RS, RMU, IFLG, VX_SR, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT NONE
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS
      REAL(q), INTENT(IN)                :: RMU
      INTEGER, INTENT(IN)                :: IFLG
      TYPE(SIMD_REAL64), INTENT(OUT)     :: VX_SR
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      REAL(q)                            :: CX, VX_, RHO, QFAC, QF, A, FRAC, Z, T, ERRF_RESULT, EX_SR_RESULT
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: CX1 = 0.9163305865663_q
      REAL(q), PARAMETER                 :: CX2 = 1.1545041946774_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q    ! 1._q/3._q
      REAL(q), PARAMETER                 :: QFAC1 = 3.093667726280135532590521_q   ! (3._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: QFAC2 = 3.897777089720753540547094_q   ! (6._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: SQRT_PI = 1.772453850905515881919428_q ! sqrt(pi)
      REAL(q), PARAMETER                 :: P_1 = 0.620350490899400086597382_q     ! (1._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P_2 = 1.333333333333333259318465_q     ! 4._q/3._q
      REAL(q), PARAMETER                 :: P_3 = 2.666666666666666518636930_q     ! 8._q/3._q
      REAL(q), PARAMETER                 :: P_4 = 0.006631455962162306067009_q     ! 1._q/(48._q*pi)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(QFAC, CX, TEMP_1, VX_, QF, A, TEMP_5, Z, T, TEMP_2, TEMP_3, TEMP_4, FRAC, EX_SR_RESULT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (IFLG .EQ. 1) THEN
            QFAC = QFAC1
            CX = CX1
         ELSE
            QFAC = QFAC2
            CX = CX2
         ENDIF
         TEMP_1 = 1.0_q / RS%X(II)
         VX_ = -CX * TEMP_1
! This version is faster and potentially more precise
         QF = QFAC * P_1 * TEMP_1
         A = RMU / (2.0_q * QF)
         TEMP_5 = 0.0_q
! Test on the value of A
         IF (A .LT. 1.0E-9_q) THEN
         ! Limit for small A
            TEMP_5 = P_2 * VX_
         ENDIF
         TEMP_1 = 1.0_q / A
         IF ((A .GE. 1.0E-9_q) .AND. (A .LE. 100.0_q)) THEN
         ! Intermediate Values of A
         ! Inline EX_SR : computations already carried out can be reused
         ! Inline ERRF
            Z = 0.5_q * TEMP_1
            T = 1.0_q / (1.0_q + 0.5_q * Z)
            ERRF_RESULT = -Z * Z - 1.26551223_q + T * (1.00002368_q + T * (.37409196_q +  &
            &  T * (0.09678418_q + T * (-0.18628806_q + T * (0.27886807_q + T * (-1.13520398_q +  &
            &  T * (1.48851587_q + T * (-0.82215223_q + T * 0.17087277_q))))))))
            ERRF_RESULT = 1.0_q - T * SIMD_EXP(ERRF_RESULT)
         ! End of Inline ERRF
            TEMP_2 = 4.0_q * A * A * A
            TEMP_3 = SIMD_EXP(-0.25_q * TEMP_1 * TEMP_1)
            FRAC = (1.0_q - P_3 * A * (SQRT_PI * ERRF_RESULT +  &
            &  (2.0_q * A - TEMP_2) * TEMP_3 - 3.0_q * A + TEMP_2))
            EX_SR_RESULT = VX_ * FRAC
         ! End of Inline EX_SR
            TEMP_2 = 4.0_q * A * A
            TEMP_4 = 2.0_q * TEMP_2
            FRAC = (1.0_q - TEMP_2 * TEMP_4 * (TEMP_3 - 1.0_q) - TEMP_4)
            TEMP_5 = EX_SR_RESULT + VX_ * FRAC * THRD
         ENDIF
         IF ((A .GT. 100.0_q) .AND. (A .LE. 1.0E+9_q)) THEN
         ! Development for large A
            TEMP_5 = -RMU * P_4 * TEMP_1 * TEMP_1 * TEMP_1
         ENDIF
         IF (MASK0%X(II)) THEN
            VX_SR%X(II) = TEMP_5
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE VX_SR_SIMD


      SUBROUTINE CALC_EXCHWPBE_SP_SIMD(XC, D, S,  &
      &  EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR,  &
      &  EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR, MASK0, ID_XC)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT D : DENSITY
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(6 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EXWPBE)
!  OUTPUT:  DERIVATIVE W.R.T. DENSITY (EXWPBED)
!  OUTPUT:  DERIVATIVE W.R.T. GRADIENT DENSITY (EXWPBEDD)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE wpbe
      USE setexm_struct_def, ONLY : xc_info
      USE xc_name
      USE simd
      IMPLICIT NONE
      TYPE (xc_info) XC
      TYPE(SIMD_REAL64), INTENT(IN)      :: D, S
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER, INTENT(IN)                :: ID_XC
      TYPE(SIMD_REAL64)                  :: RS, SW, LOGRHO1, LOGRHO2, LOGRHO3, SMIN, SPLS, XKF
      TYPE(SIMD_REAL64)                  :: F, FSR, FSR_MIN, FSR_PLS, EXPBED, EXPBEDD, DD, DFDR, DFDDD, INV_D
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: FS, FLR
      REAL(q)                            :: LDASCREEN
! Numerical derivatives test
      REAL(q)                            :: XKFMIN, XKFPLS
      REAL(q)                            :: DDMIN, DDPLS
      REAL(q)                            :: DMIN, DPLS, XFKMIN, XFKPLS
      REAL(q)                            :: EXWPBE, EXPBE
      REAL(q)                            :: TEMP_1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
! Parameters to calculate the PBE enhancement factor
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q
      REAL(q), PARAMETER                 :: AX = -0.738558766382022405884230032680836_q
      REAL(q), PARAMETER                 :: UM = 0.2195149727645171_q
      REAL(q), PARAMETER                 :: UMSOL = 0.12345679012345679012_q
      REAL(q), PARAMETER                 :: UK = 0.8040_q
      REAL(q), PARAMETER                 :: UL = UM / UK
      REAL(q), PARAMETER                 :: ULSOL = UMSOL / UK
      REAL(q), PARAMETER                 :: P_1 = 0.238732414637843004534545_q ! 3._q/(4._q*pi)
      REAL(q), PARAMETER                 :: P_2 = 1.333333333333333259318465_q ! 4._q/3._q

      LDASCREEN=XC%LDASCREEN

      CALL INIT_WPBE_TABLE(XC)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! Due to the fact that most of the times only the GGA corrections to the
! LDA exchange are wanted we will have to recalculate the LDA contributions
! in order to be able to subtract them from our PBE result.
! (Georg thinks this makes a lot of sense, but I guess there are not many
! people who would agree with this!)
         RS%X(II) = SIMD_POW(P_1 / D%X(II), THRD)
! continued below
      ENDDO

      CALL EX_SR_SIMD(RS, LDASCREEN * AUTOA, 1, EXLDA_SR, MASK0)

      CALL VX_SR_SIMD(RS, LDASCREEN * AUTOA, 1, EXDLDA_SR, MASK0)

      TRUE_FOR_ANY = .FALSE.
__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(TEMP_1)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         SW%X(II) = ABS(S%X(II))
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            TEMP_1 = 1.0_q / RS%X(II)
            EXLDA%X(II) = (-0.9163305865663_q * TEMP_1) * D%X(II) * 0.5_q
            EXDLDA%X(II) = (-1.2217741154217_q * TEMP_1) * 0.5_q
            EXLDA_SR%X(II) = EXLDA_SR%X(II) * D%X(II) * 0.5_q
            EXLDA_LR%X(II) = EXLDA%X(II) - EXLDA_SR%X(II)
            EXDLDA_SR%X(II) = EXDLDA_SR%X(II) * 0.5_q
            EXDLDA_LR%X(II) = EXDLDA%X(II) - EXDLDA_SR%X(II)
            EXWPBE_SR%X(II) = 0.0_q
            EXWPBE_LR%X(II) = 0.0_q
            EXWPBED_SR%X(II) = 0.0_q
            EXWPBED_LR%X(II) = 0.0_q
            EXWPBEDD_SR%X(II) = 0.0_q
            EXWPBEDD_LR%X(II) = 0.0_q
! Interpolate the short range part of the enhancement
! coefficient from the tables            
            IF ((SW%X(II) .NE. 0.0_q) .AND. (D%X(II) .GE. 0.0_q)) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, FS, EXPBE, DMIN, XKFMIN, DPLS, XKFPLS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! PBE quantities
         IF (ID_XC==ID_XC_PBESOL) THEN
            TEMP_1 = 1.0_q / (1.0_q + ULSOL * SW%X(II) * SW%X(II))
            FS = 2.0_q * UMSOL * TEMP_1 * TEMP_1
         ELSE
            TEMP_1 = 1.0_q / (1.0_q + UL * SW%X(II) * SW%X(II))
            FS = 2.0_q * UM * TEMP_1 * TEMP_1
         ENDIF
         F%X(II) = 1.0_q + UK - UK * TEMP_1
         EXPBE = F%X(II) * EXLDA%X(II)
         INV_D%X(II) = 1.0_q / D%X(II)
         EXPBED%X(II) = P_2 * (F%X(II) - SW%X(II) * SW%X(II) * FS) * EXLDA%X(II) * INV_D%X(II)
         EXPBEDD%X(II) = 0.5_q * AX * 0.3232409194_q * SW%X(II) * FS
! Cutoff criterion to enforce local Lieb-Oxford
! this ensures that the enhancement factor does not exceed the
! original one
! get values from function directly (slower)
!     CALL EXCHWPBE_R(LDASCREEN*AUTOA,D,SW,FSR)
! Numerical derivatives with interpolated values:
! Firstly: derivatives w.r.t rho
         XKF%X(II) = SIMD_POW(3.0_q * PI * PI * D%X(II), THRD)
         DD%X(II) = SW%X(II) * (2.0_q * XKF%X(II) * D%X(II))
         DMIN = 0.9999_q * D%X(II)
         LOGRHO2%X(II) = SIMD_LOG(DMIN)
         XKFMIN = SIMD_POW(3.0_q * PI * PI * DMIN, THRD)
         SMIN%X(II) = DD%X(II) / (DMIN * XKFMIN * 2.0_q)
         DPLS = 1.0001_q * D%X(II)
         LOGRHO3%X(II) = SIMD_LOG(DPLS)
         XKFPLS = SIMD_POW(3.0_q * PI * PI * DPLS, THRD)
         SPLS%X(II) = DD%X(II) / (DPLS * XKFPLS * 2.0_q)
         LOGRHO1%X(II) = SIMD_LOG(D%X(II))
      ENDDO

      CALL WPBE_SPLINE_SIMD(LOGRHO1, SW, FSR, MASK1)

      CALL WPBE_SPLINE_SIMD(LOGRHO2, SMIN, FSR_MIN, MASK1)

      CALL WPBE_SPLINE_SIMD(LOGRHO3, SPLS, FSR_PLS, MASK1)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(FLR, DDMIN, DDPLS, TEMP_1)
      DO II=0,SIMD_WIDTH_64BIT-1
         DFDR%X(II) = (FSR_PLS%X(II) - FSR_MIN%X(II)) * 5000.0_q * INV_D%X(II)
! Find the complementary long range part
         FLR = F%X(II) - FSR%X(II)
! E_sr = fsr*Exlda
         IF (MASK1%X(II)) THEN
            EXWPBE_SR%X(II) = FSR%X(II) * EXLDA%X(II)
            EXWPBE_LR%X(II) = FLR * EXLDA%X(II)
         ENDIF
! Secondly: derivatives w.r.t gradient(rho)
         DDMIN = 0.9999_q * DD%X(II)
         DDPLS = 1.0001_q * DD%X(II)
         TEMP_1 = 0.5_q / (D%X(II) * XKF%X(II))
         SMIN%X(II) = DDMIN * TEMP_1
         SPLS%X(II) = DDPLS * TEMP_1
      ENDDO

      CALL WPBE_SPLINE_SIMD(LOGRHO1, SMIN, FSR_MIN, MASK1)

      CALL WPBE_SPLINE_SIMD(LOGRHO1, SPLS, FSR_PLS, MASK1)

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1)
      DO II=0,SIMD_WIDTH_64BIT-1
         DFDDD%X(II) = (FSR_PLS%X(II) - FSR_MIN%X(II)) * 5000.0_q / DD%X(II)
! dE_sr/drho = dfsr/drho*Exlda + fsr*dExlda/drho
         IF (MASK1%X(II)) THEN
            EXWPBED_SR%X(II) = DFDR%X(II) * EXLDA%X(II) + P_2 * FSR%X(II)*EXLDA%X(II) * INV_D%X(II)
            EXWPBED_LR%X(II) = EXPBED%X(II) - EXWPBED_SR%X(II)
! dE_sr/d(grad rho) = dfsr/dgrad(rho)*Exlda
         ENDIF
!gK correct the sign
         IF (S%X(II) .LT. 0.0_q) THEN
            TEMP_1 = -1.0_q
         ELSE
            TEMP_1 = 1.0_q
         ENDIF
         IF (MASK1%X(II)) THEN
            EXWPBEDD_SR%X(II) = DFDDD%X(II) * EXLDA%X(II) * TEMP_1
            EXWPBEDD_LR%X(II) = (EXPBEDD%X(II) - DFDDD%X(II) * EXLDA%X(II)) * TEMP_1
         ENDIF
      ENDDO
!gK end
      RETURN
      END SUBROUTINE CALC_EXCHWPBE_SP_SIMD
#endif

!----------------------------------------------------------------------
!> jP: adding PBEsol
!>     Note that this is a plain copy-paste - just the necessary constants are changed!
!>     The quickest way to program
!>     PBEsol is to take your existing PBE subroutine, and alter the two values
!>     mu and beta:  In PBE, mu=0.219 and beta=0.067 these become mu=10/81
!>     and beta=0.046 in PBEsol.   The routines below are simple modifications
!>     of the PBE subroutines, to produce PBEsol.
!---------------------------------------------------------------------
!> The P B E
!> this routine was kindly supplied by Bjork Hammer
!> it is essentially identical to Burkes routine but includes
!> the revised P B E routines
!---------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
SUBROUTINE EXCHPBESOL(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd, &
     &                   ukfactor)
GPACC !$ACC ROUTINE SEQ
  !----------------------------------------------------------------------
  !  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
  !  K Burkes modification of PW91 codes, May 14, 1996
  !  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !  INPUT rho : DENSITY
  !  INPUT rhothrd : DENSITY^(1/3)
  !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
  !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! References:
  ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
  ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
  !     {\bf 40},  3399  (1989) (E).
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! Formulas:
  !       e_x[unif]=ax*rho^(4/3)  [LDA]
  ! ax = -0.75*(3/pi)^(1/3)
  !       e_x[PBE]=e_x[unif]*FxPBE(s)
  !       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
  ! uk, ul defined after [a](13)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  USE prec
  IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(EXCHPBESOL)
  parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
  parameter(pi=3.14159265358979323846264338327950_q)
  parameter(ax=-0.738558766382022405884230032680836_q)
  parameter(um=0.12345679012345679012_q,uk1=0.8040_q,ul1=um/uk1)

  !----------------------------------------------------------------------
  ! construct LDA exchange energy density
  exunif=AX*rhothrd
  exlda=exunif*rho
  exdlda=exunif*thrd4
  !----------------------------------------------------------------------
  ! construct PBE enhancement factor
  S2 = S*S
  !----------------------------------------------------------------------
  if (ukfactor.ne.0.0_q) then
     ! These are the PBE96 and revPBE98 functionals
     ! scale uk with a factor
     uk = uk1*ukfactor
     ul = ul1/ukfactor
     P0=1._q+ul*S2
     FxPBE = 1._q+uk-uk/P0
     expbe = exlda*FxPBE
     !----------------------------------------------------------------------
     !  ENERGY DONE. NOW THE POTENTIAL:
     !  find first derivatives of Fx w.r.t s.
     !  Fs=(1/s)*d FxPBE/ ds
     Fs=2._q*um/(P0*P0)
  else
     ! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
     P0=exp(-ul1*S2)
     FxPBE = 1._q+uk1*(1.0_q-P0)
     expbe = exlda*FxPBE
     !----------------------------------------------------------------------
     !  ENERGY DONE. NOW THE POTENTIAL:
     !  find first derivatives of Fx w.r.t s.
     !  Fs=(1/s)*d FxPBE/ ds
     Fs=2._q*um*P0
  endif
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! calculate the partial derivatives of ex wrt n and |grad(n)|
  !  0.3232409194=(3*pi^2)^(-1/3)
  exd =exunif*THRD4*(FxPBE-S2*Fs)
  exdd=0.5_q*ax*0.3232409194_q*S*Fs
  RETURN
END SUBROUTINE EXCHPBESOL

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCHPBESOL_SIMD(RHO, RHOTHRD, S, EXLDA, EXPBE, EXDLDA, EXD, EXDD, UKFACTOR, MASK0)
  !----------------------------------------------------------------------
  !  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
  !  K Burkes modification of PW91 codes, May 14, 1996
  !  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !  INPUT rho : DENSITY
  !  INPUT rhothrd : DENSITY^(1/3)
  !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
  !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! References:
  ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
  ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
  !     {\bf 40},  3399  (1989) (E).
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! Formulas:
  !       e_x[unif]=ax*rho^(4/3)  [LDA]
  ! ax = -0.75*(3/pi)^(1/3)
  !       e_x[PBE]=e_x[unif]*FxPBE(s)
  !       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
  ! uk, ul defined after [a](13)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO, RHOTHRD, S, UKFACTOR
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXLDA, EXPBE, EXDLDA, EXD, EXDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
!      TYPE(SIMD_REAL64)                  :: S2, FXPBE, EXUNIF, FS
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q, THRD4 = 4.0_q / 3.0_q
      REAL(q), PARAMETER                 :: PI = 3.14159265358979323846264338327950_q
      REAL(q), PARAMETER                 :: AX = -0.738558766382022405884230032680836_q
      REAL(q), PARAMETER                 :: UM = 0.12345679012345679012_q, UK1 = 0.8040_q, UL1 = UM / UK1
      REAL(q), PARAMETER                 :: P_1 = 0.5_q * AX * 0.3232409194_q

__CPU !$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(EXUNIF, S2, UK, UL, P0, FXPBE, FS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
!----------------------------------------------------------------------
! construct LDA exchange energy density
         EXUNIF = AX * RHOTHRD%X(II)
         IF (MASK0%X(II)) THEN
            EXLDA%X(II) = EXUNIF * RHO%X(II)
            EXDLDA%X(II) = EXUNIF * THRD4
         ENDIF
!----------------------------------------------------------------------
! construct PBE enhancement factor
         S2 = S%X(II) * S%X(II)
!----------------------------------------------------------------------
         IF (UKFACTOR%X(II) .NE. 0.0_q) THEN
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
            UK = UK1 * UKFACTOR%X(II)
            UL = UL1 / UKFACTOR%X(II)
            P0 = 1.0_q / (1.0_q + UL * S2)
            FXPBE = 1.0_q + UK - UK * P0
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
            FS = 2.0_q * UM * P0 * P0
         ELSE
     ! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
            P0 = SIMD_EXP(-UL1 * S2)
            FXPBE = 1.0_q + UK1 * (1.0_q - P0)
     !----------------------------------------------------------------------
     !  ENERGY DONE. NOW THE POTENTIAL:
     !  FIND FIRST DERIVATIVES OF FX W.R.T S.
     !  FS=(1/S)*D FXPBE/ DS
            FS = 2.0_q * UM * P0
         ENDIF
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! calculate the partial derivatives of ex wrt n and |grad(n)|
  !  0.3232409194=(3*pi^2)^(-1/3)
         IF (MASK0%X(II)) THEN
            EXPBE%X(II) = EXLDA%X(II) * FXPBE
            EXD%X(II) = EXUNIF * THRD4 * (FXPBE - S2 * FS)
            EXDD%X(II) = P_1 * S%X(II) * FS
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCHPBESOL_SIMD
#endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! jP ######################################################################
!----------------------------------------------------------------------
SUBROUTINE CORunspPBESOL(RS,EC,VC,sk, &
     &                  T,H,DVC,ecdd,lgga)
GPACC !$ACC ROUTINE SEQ
  !----------------------------------------------------------------------
  !  Official PBE correlation code. K. Burke, May 14, 1996.
  !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
  !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
  !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
  !       : lgga=flag to do gga (0=>LSD only)
  !       : lpot=flag to do potential (0=>energy only)
  !  output: ec=lsd correlation energy from [a]
  !        : vcup=lsd up correlation potential
  !        : vcdn=lsd dn correlation potential
  !        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
  !        : dvcup=nonlocal correction to vcup
  !        : dvcdn=nonlocal correction to vcdn
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
  USE prec
  IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORunspPBESOL)
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
  parameter(bet=0.046_q,delt=bet/gamma)
  parameter(eta=1.e-12_q)
GPACC !$ACC ROUTINE (GCOR2) SEQ
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
  CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
       &    0.49294_q,rtrs,EU,EURS)
  EC = EU
  ! check for zero energy, immediate return if true
  IF (EC==0.) THEN
     H=0; DVC=0; ecdd=0
     RETURN
  ENDIF
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! LSD potential from [c](A1)
  ! ECRS = dEc/drs [c](A2)
  ! ECZET=dEc/dzeta [c](A3)
  ! FZ = dF/dzeta [c](A4)
  ECRS = EURS
  VC = EC -RS*ECRS/3._q
  if (.not.lgga) return
  !----------------------------------------------------------------------
  ! PBE correlation energy
  ! G=phi(zeta), given after [a](3)
  ! DELT=bet/gamma
  ! B=A of [a](8)
  PON=-EC/(gamma)
  B = DELT/(DEXP(PON)-1._q)
  B2 = B*B
  T2 = T*T
  T4 = T2*T2
  RS2 = RS*RS
  RS3 = RS2*RS
  Q4 = 1._q+B*T2
  Q5 = 1._q+B*T2+B2*T4
  H = (BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
  T6 = T4*T2
  RSTHRD = RS/3._q
  FAC = DELT/B+1._q
  BEC = B2*FAC/(BET)
  Q8 = Q5*Q5+DELT*Q4*Q5*T2
  Q9 = 1._q+2._q*B*T2
  hB = -BET*B*T6*(2._q+B*T2)/Q8
  hRS = -RSTHRD*hB*BEC*ECRS
  hT = 2._q*BET*Q9/Q8
  DVC = H+HRS-7.0_q*T2*HT/6._q
  ecdd=0.5_q/sk*t*ht
  RETURN
END SUBROUTINE CORunspPBESOL
!----------------------------------------------------------------------


!######################################################################
! jP ----------------------------------------------------------------------
SUBROUTINE CORPBESOL(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
GPACC !$ACC ROUTINE SEQ
  !----------------------------------------------------------------------
  !  Official PBE correlation code. K. Burke, May 14, 1996.
  !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
  !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
  !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
  !       : lgga=flag to do gga (0=>LSD only)
  !       : lpot=flag to do potential (0=>energy only)
  !  output: ec=lsd correlation energy from [a]
  !        : vcup=lsd up correlation potential
  !        : vcdn=lsd dn correlation potential
  !        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
  !        : dvcup=nonlocal correction to vcup
  !        : dvcdn=nonlocal correction to vcdn
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
  USE prec
  IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORPBESOL)
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
  parameter(bet=0.046_q,delt=bet/gamma)
  parameter(eta=1.e-12_q)
GPACC !$ACC ROUTINE(GCOR2) SEQ
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
  CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
       &    0.49294_q,rtrs,EU,EURS)
  CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
       &    0.62517_q,rtRS,EP,EPRS)
  CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
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
  if (.not.lgga) return
  !----------------------------------------------------------------------
  ! PBE correlation energy
  ! G=phi(zeta), given after [a](3)
  ! DELT=bet/gamma
  ! B=A of [a](8)
  G3 = G**3
  PON=-EC/(G3*gamma)
  B = DELT/(DEXP(PON)-1._q)
  B2 = B*B
  T2 = T*T
  T4 = T2*T2
  RS2 = RS*RS
  RS3 = RS2*RS
  Q4 = 1._q+B*T2
  Q5 = 1._q+B*T2+B2*T4
  H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
  G4 = G3*G
  T6 = T4*T2
  RSTHRD = RS/3._q
  GZ=(((1._q+zet)**2+eta)**sixthm- &
       &((1._q-zet)**2+eta)**sixthm)/3._q
  FAC = DELT/B+1._q
  BG = -3._q*B2*EC*FAC/(BET*G4)
  BEC = B2*FAC/(BET*G3)
  Q8 = Q5*Q5+DELT*Q4*Q5*T2
  Q9 = 1._q+2._q*B*T2
  hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
  hRS = -RSTHRD*hB*BEC*ECRS
  hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
  hT = 2._q*BET*G3*Q9/Q8
  COMM = H+HRS-7.0_q*T2*HT/6._q
  PREF = HZ-GZ*T2*HT/G
  COMM = COMM-PREF*ZET
  DVCUP = COMM + PREF
  DVCDN = COMM - PREF
  ecdd=0.5_q/(sk*g)*t*ht
  RETURN
END SUBROUTINE CORPBESOL
! jP: adding PBEsol


!***********************************************************************
!
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> CRPKZB                       C(ZETA,XI)
!> D_CRPKZB_D1,D_CRPKZB_D2      THE DERIVATIVE OF C(ZETA,XI) WRT n
!> D_CRPKZB_DD1,D_CRPKZB_DD2    THE DERIVATIVE OF C(ZETA,XI) WRT |GRAD n|
!
!***********************************************************************

      SUBROUTINE COE_RPKZB(&
     &   RU,RD,DRU,DRD,DRT, &
     &   CRPKZB,D_CRPKZB_D1,D_CRPKZB_D2,D_CRPKZB_DD1,D_CRPKZB_DD2,D_CRPKZB_DD12)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(COE_RPKZB)

! the following parameters are given by Perdew et.al.
! the coefficients for the correlatin part
      PARAMETER (CF1=0.59_q)
      PARAMETER (CF2=0.9269_q)
      PARAMETER (CF3=0.6225_q)
      PARAMETER (CF4=2.1540_q)

! other parameters
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

      CRPKZB=0._q
      D_CRPKZB_D1=0._q;D_CRPKZB_D2=0._q;
      D_CRPKZB_DD1=0._q;D_CRPKZB_DD2=0._q;
      D_CRPKZB_DD12=0._q;

      RT=RU+RD
      YA=DRU**2._q
      YB=DRD**2._q
      Y=DRT**2._q
! YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2._q

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
      OPZ = ONE + ZETA
      OMZ = ONE - ZETA
      ZETA2=ZETA*ZETA
      OPZETA=ONE/(ONE+ZETA)
      OMZETA=ONE/(ONE-ZETA)
      D_ZETA_D1=TWO*RD/RT**TWO
      D_ZETA_D2=-TWO*RU/RT**TWO

! |GRAD ZETA|,
! IN SOME EXTREME CASES, THE TERM (YA*RD**TWO-TWO*RU*RD*YC+YB*RU**TWO)
! GOES TO NEGATIVE BECAUSE OF THE NUMERICAL PRESICION OF COMPUTERS. SO, WE USE 2/4 AS THE POWER
! INSTEAD OF 1/2
      DEL_ZETA2 = (TWO*YA*OMZ + TWO*YB*OPZ - OPZ*OMZ*Y)/RT**2
      D_DEL_ZETA2_DYA = TWO*OMZ/RT**2
      D_DEL_ZETA2_DYB = TWO*OPZ/RT**2
      D_DEL_ZETA2_DY = -OPZ*OMZ/RT**2
      D_DEL_ZETA2_DZETA = TWO*(ZETA*Y + YB - YA)/RT**2
      D_DEL_ZETA2_D1 = -TWO*DEL_ZETA2/RT + D_DEL_ZETA2_DZETA*OMZ/RT
      D_DEL_ZETA2_D2 = -TWO*DEL_ZETA2/RT - D_DEL_ZETA2_DZETA*OPZ/RT

! XI
      XIDEN2=4._q*(THREE*PISQ*RT)**THRD2
      XI2=DEL_ZETA2/XIDEN2
      D_XI2_D1=(D_DEL_ZETA2_D1-THRD2*DEL_ZETA2/RT)/XIDEN2
      D_XI2_D2=(D_DEL_ZETA2_D2-THRD2*DEL_ZETA2/RT)/XIDEN2
      D_XI2_DYA=D_DEL_ZETA2_DYA/XIDEN2
      D_XI2_DYB=D_DEL_ZETA2_DYB/XIDEN2
      D_XI2_DY=D_DEL_ZETA2_DY/XIDEN2

      CNUM=CF1+(CF2+(CF3+CF4*ZETA2)*ZETA2)*ZETA2
      D_CNUM_D1=D_ZETA_D1*(TWO*CF2+(FOUR*CF3+SIX*CF4*ZETA2)*ZETA2)*ZETA
      D_CNUM_D2=D_ZETA_D2*(TWO*CF2+(FOUR*CF3+SIX*CF4*ZETA2)*ZETA2)*ZETA

      OPZETA43=OPZETA**THRD4
      OMZETA43=OMZETA**THRD4
      CDENL=ONE+XI2*(OPZETA43+OMZETA43)/TWO
      CDEN=CDENL**FOUR
      D_CDENL_DXI2=(OPZETA43+OMZETA43)/TWO
      D_CDENL_DZETA=-THRD2*XI2*(OPZETA43*OPZETA-OMZETA43*OMZETA)
      D_CDENL_D1=D_CDENL_DXI2*D_XI2_D1+D_CDENL_DZETA*D_ZETA_D1
      D_CDENL_D2=D_CDENL_DXI2*D_XI2_D2+D_CDENL_DZETA*D_ZETA_D2
      D_CDENL_DYA=D_CDENL_DXI2*D_XI2_DYA
      D_CDENL_DYB=D_CDENL_DXI2*D_XI2_DYB
      D_CDENL_DY=D_CDENL_DXI2*D_XI2_DY

      D_CDEN_D1=FOUR*CDENL**THREE*D_CDENL_D1
      D_CDEN_D2=FOUR*CDENL**THREE*D_CDENL_D2
      D_CDEN_DYA=FOUR*CDENL**THREE*D_CDENL_DYA
      D_CDEN_DYB=FOUR*CDENL**THREE*D_CDENL_DYB
      D_CDEN_DY=FOUR*CDENL**THREE*D_CDENL_DY


! OUTPUT CRPKZB,,D_CRPKZB_D1,D_CRPKZB_D2,D_CRPKZB_DD1,D_CRPKZB_DD2
      CRPKZB=CNUM/CDEN
      D_CRPKZB_D1=(CDEN*D_CNUM_D1-CNUM*D_CDEN_D1)/CDEN**TWO
      D_CRPKZB_D2=(CDEN*D_CNUM_D2-CNUM*D_CDEN_D2)/CDEN**TWO
      D_CRPKZB_DD1=-TWO*DRU*CNUM*D_CDEN_DYA/CDEN**TWO
      D_CRPKZB_DD2=-TWO*DRD*CNUM*D_CDEN_DYB/CDEN**TWO
      D_CRPKZB_DD12=-TWO*DRT*CNUM*D_CDEN_DY/CDEN**TWO

      RETURN
      END SUBROUTINE COE_RPKZB


!***********************************************************************
!
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> EPPGGA                       GGA ENERGY PER PARTICLE
!> EPPGGA_D1,EPPGGA_D2          THE DERIVATIVE OF EPPGGA WRT n
!> EPPGGA_DD1,EPPGGA_DD2        THE DERIVATIVE OF EPPGGA WRT |GRAD n|
!
!***********************************************************************
      SUBROUTINE CORGGA_REVTPSS(&
     &   RU,RD,DRU,DRD,DRT, &
     &   EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD12)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant

      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORGGA_REVTPSS)
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

      PARAMETER (GAMMA=0.03109069086965489503494086371273_q)
      PARAMETER (BETA_mb=0.06672455060314922_q)
GPACC !$ACC ROUTINE(CORPBE_REVTPSS) SEQ

      EPPGGA  =0._q
      EPPGGA_D1=0._q
      EPPGGA_D2=0._q
      EPPGGA_DD12=0._q

      RT=RU+RD
      YA=DRU**2._q
      YB=DRD**2._q
      Y=DRT**2._q
!     YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2._q

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
      DZETAD1=TWO*RD/RT**TWO
      DZETAD2=-TWO*RU/RT**TWO


      DTHRD=exp(log(RT)*THRD)
      RS=(0.75_q/PI)**THRD/DTHRD
      DRSD1=-THRD/RT*RS
      DRSD2=-THRD/RT*RS


      PHI = (exp((TWO*THRD)*log(1._q+ZETA)) &
                +exp((TWO*THRD)*log(1._q-ZETA)))/2._q
      D_PHI_DZETA=THRD*((ONE+ZETA)**(-THRD)-(ONE-ZETA)**(-THRD))
      D_PHI_D1=D_PHI_DZETA*DZETAD1
      D_PHI_D2=D_PHI_DZETA*DZETAD2

      AFIX_T=SQRT(PI/FOUR)*(9._q*PI/FOUR)**(ONE/SIX)
      S=DRT/(TWO*(THREE*PISQ)**THRD*RT**THRD4)
      DSD1=-THRD4*S/RT
      DSD2=-THRD4*S/RT
      DSDD12 = ONE/(TWO*(THREE*PISQ)**THRD*RT**THRD4)

      T=AFIX_T*S/SQRT(RS)/PHI
      T2 = T*T
      T4 = T2*T2

      DTD1=AFIX_T*(PHI*RS*DSD1-0.5_q*S*PHI*DRSD1-S*RS*D_PHI_D1)/(PHI**TWO*RS**(THREE/TWO))
      DTD2=AFIX_T*(PHI*RS*DSD2-0.5_q*S*PHI*DRSD2-S*RS*D_PHI_D2)/(PHI**TWO*RS**(THREE/TWO))
      DTDD12=AFIX_T*DSDD12/(SQRT(RS)*PHI)

      FK=(3._q*PI*PI)**THRD*DTHRD
      SK = SQRT(4.0_q*FK/PI)
! Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
      CALL CORPBE_revtpss(RS,ZETA,EC,VCUPLDA,VCDNLDA,PHI,SK, &
           T,H,DVCUP,DVCDN,ECQ,.FALSE.)
! OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
!      EPPGGA=EC
!      EPPGGA_D1=(VCUPLDA-EC)/RT
!      EPPGGA_D2=(VCDNLDA-EC)/RT
      D_EC_DZETA=(VCUPLDA-VCDNLDA)/TWO
      D_EC_DRS=(EC-ZETA*D_EC_DZETA-(VCUPLDA+VCDNLDA)/TWO)*THREE/RS

! THE RS DEPENDENCE OF BETA
      AFACTOR=0.1_q
      BFACTOR=0.1778_Q
      BETA_NUM=ONE + AFACTOR*RS
      BETA_DEN=ONE+ BFACTOR*RS
      D_BETA_NUM=AFACTOR
      D_BETA_DEN=BFACTOR
      BETA = BETA_MB*BETA_NUM/BETA_DEN
      D_BETA_DRS=BETA_MB*(BETA_DEN*D_BETA_NUM-BETA_NUM*D_BETA_DEN)/BETA_DEN**TWO


      PHI3=PHI**THREE
      PON=-EC/(PHI3*gamma)
      W=DEXP(PON)-ONE
      D_W_DRS=-(W+ONE)*D_EC_DRS/(GAMMA*PHI3)
      D_W_DZETA=-(W+ONE)/(GAMMA*PHI3)*(D_EC_DZETA-THREE*EC*D_PHI_DZETA/PHI)


      A=BETA/(GAMMA*W)
      D_A_DRS=(W*D_BETA_DRS-BETA*D_W_DRS)/(GAMMA*W**TWO)
      D_A_DZETA=-BETA*D_W_DZETA/(GAMMA*W**TWO)

      V=A*T2
      D_V_DRS=T2*D_A_DRS
      D_V_DZETA=T2*D_A_DZETA
      D_V_DT=TWO*A*T

      FUNKG=ONE/(ONE+V+V**TWO)
      D_FUNKG_DV=-(ONE+TWO*V)*FUNKG**TWO


      HCORE=ONE+W*(ONE-FUNKG)
      AH=GAMMA*PHI3
      H=AH*DLOG(HCORE)

      DH1=ONE-FUNKG
      DH2=W*D_FUNKG_DV
      D_H_DRS=(DH1*D_W_DRS-DH2*D_V_DRS)*AH/HCORE
      D_H_DZETA=THREE*H*D_PHI_DZETA/PHI+(DH1*D_W_DZETA-DH2*D_V_DZETA)*AH/HCORE
      D_H_DT=(-DH2*D_V_DT)*AH/HCORE

! OUTPUT EPPGGA AND ITS DERIVATIVES EPPGGA_D1,EPPGGA_D2, EPPGGA_DD1, EPPGGA_DD2
      EPPGGA=EC+H
      EPPGGA_D1=D_EC_DRS*DRSD1+D_EC_DZETA*DZETAD1+D_H_DRS*DRSD1+D_H_DZETA*DZETAD1+D_H_DT*DTD1
      EPPGGA_D2=D_EC_DRS*DRSD2+D_EC_DZETA*DZETAD2+D_H_DRS*DRSD2+D_H_DZETA*DZETAD2+D_H_DT*DTD2
      EPPGGA_DD12=D_H_DT*DTDD12

      RETURN
      END SUBROUTINE CORGGA_REVTPSS


! Modified by Adrienn Ruzsinszky
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE_revtpss(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lmetagga=flag to do metagga (revTPSS)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORPBE_REVTPSS)
      logical lgga
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet_mb=0.06672455060314922_q)
      parameter(eta=1.e-12_q)
GPACC !$ACC ROUTINE(gcor2) SEQ
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
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
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
!write(*,*)'rs,VCUP,VCDN',rs,VCUP,VCDN
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      bet = bet_mb*(1._q + 0.1_q*RS)/(1._q + 0.1778_q*RS)

      delt=bet/gamma
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3._q
      GZ=(((1._q+zet)**2+eta)**sixthm- &
     &((1._q-zet)**2+eta)**sixthm)/3._q
      FAC = DELT/B+1._q
      BG = -3._q*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2._q*BET*G3*Q9/Q8
      COMM = H+HRS-7.0_q*T2*HT/6._q
      PREF = HZ-GZ*T2*HT/G
      COMM = COMM-PREF*ZET
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      ecdd=0.5_q/(sk*g)*t*ht
      RETURN
      END SUBROUTINE CORPBE_revtpss

!***********************************************************************
!
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> CRPKZB                       C(ZETA,XI)
!> D_CRPKZB_D1,D_CRPKZB_D2      THE DERIVATIVE OF C(ZETA,XI) WRT n
!> D_CRPKZB_DD1,D_CRPKZB_DD2    THE DERIVATIVE OF C(ZETA,XI) WRT |GRAD n|
!
!***********************************************************************

      SUBROUTINE COE_RPKZB_TPSS(&
     &   RU,RD,DRU,DRD,DRT, &
     &   CRPKZB,D_CRPKZB_D1,D_CRPKZB_D2,D_CRPKZB_DD1,D_CRPKZB_DD2,D_CRPKZB_DD12)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(COE_RPKZB_TPSS)

! the following parameters are given by Perdew et.al.
! the coefficients for the correlatin part
      PARAMETER (CF1=0.53_q)
      PARAMETER (CF2=0.87_q)
      PARAMETER (CF3=0.50_q)
      PARAMETER (CF4=2.26_q)

! other parameters
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



      CRPKZB=0._q
      D_CRPKZB_D1=0._q;D_CRPKZB_D2=0._q;
      D_CRPKZB_DD1=0._q;D_CRPKZB_DD2=0._q;D_CRPKZB_DD12=0._q;
      D_CRPKZB_DD12=0._q;

      RT=RU+RD
      YA=DRU**2._q
      YB=DRD**2._q
      Y=DRT**2._q
! YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2._q

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
      OPZ = ONE + ZETA
      OMZ = ONE - ZETA
      ZETA2=ZETA*ZETA
      OPZETA=ONE/(ONE+ZETA)
      OMZETA=ONE/(ONE-ZETA)
      D_ZETA_D1=TWO*RD/RT**TWO
      D_ZETA_D2=-TWO*RU/RT**TWO

! |GRAD ZETA|,
! IN SOME EXTREME CASES, THE TERM (YA*RD**TWO-TWO*RU*RD*YC+YB*RU**TWO)
! GOES TO NEGATIVE BECAUSE OF THE NUMERICAL PRESICION OF COMPUTERS. SO, WE USE 2/4 AS THE POWER
! INSTEAD OF 1/2
      DEL_ZETA2 = (TWO*YA*OMZ + TWO*YB*OPZ - OPZ*OMZ*Y)/RT**2
      D_DEL_ZETA2_DYA = TWO*OMZ/RT**2
      D_DEL_ZETA2_DYB = TWO*OPZ/RT**2
      D_DEL_ZETA2_DY = -OPZ*OMZ/RT**2
      D_DEL_ZETA2_DZETA = TWO*(ZETA*Y + YB - YA)/RT**2
      D_DEL_ZETA2_D1 = -TWO*DEL_ZETA2/RT + D_DEL_ZETA2_DZETA*OMZ/RT
      D_DEL_ZETA2_D2 = -TWO*DEL_ZETA2/RT - D_DEL_ZETA2_DZETA*OPZ/RT

! XI
      XIDEN2=4._q*(THREE*PISQ*RT)**THRD2
      XI2=DEL_ZETA2/XIDEN2
      D_XI2_D1=(D_DEL_ZETA2_D1-THRD2*DEL_ZETA2/RT)/XIDEN2
      D_XI2_D2=(D_DEL_ZETA2_D2-THRD2*DEL_ZETA2/RT)/XIDEN2
      D_XI2_DYA=D_DEL_ZETA2_DYA/XIDEN2
      D_XI2_DYB=D_DEL_ZETA2_DYB/XIDEN2
      D_XI2_DY=D_DEL_ZETA2_DY/XIDEN2

      CNUM=CF1+(CF2+(CF3+CF4*ZETA2)*ZETA2)*ZETA2
      D_CNUM_D1=D_ZETA_D1*(TWO*CF2+(FOUR*CF3+SIX*CF4*ZETA2)*ZETA2)*ZETA
      D_CNUM_D2=D_ZETA_D2*(TWO*CF2+(FOUR*CF3+SIX*CF4*ZETA2)*ZETA2)*ZETA

      OPZETA43=OPZETA**THRD4
      OMZETA43=OMZETA**THRD4
      CDENL=ONE+XI2*(OPZETA43+OMZETA43)/TWO
      CDEN=CDENL**FOUR
      D_CDENL_DXI2=(OPZETA43+OMZETA43)/TWO
      D_CDENL_DZETA=-THRD2*XI2*(OPZETA43*OPZETA-OMZETA43*OMZETA)
      D_CDENL_D1=D_CDENL_DXI2*D_XI2_D1+D_CDENL_DZETA*D_ZETA_D1
      D_CDENL_D2=D_CDENL_DXI2*D_XI2_D2+D_CDENL_DZETA*D_ZETA_D2
      D_CDENL_DYA=D_CDENL_DXI2*D_XI2_DYA
      D_CDENL_DYB=D_CDENL_DXI2*D_XI2_DYB
      D_CDENL_DY=D_CDENL_DXI2*D_XI2_DY

      D_CDEN_D1=FOUR*CDENL**THREE*D_CDENL_D1
      D_CDEN_D2=FOUR*CDENL**THREE*D_CDENL_D2
      D_CDEN_DYA=FOUR*CDENL**THREE*D_CDENL_DYA
      D_CDEN_DYB=FOUR*CDENL**THREE*D_CDENL_DYB
      D_CDEN_DY=FOUR*CDENL**THREE*D_CDENL_DY


! OUTPUT CRPKZB,,D_CRPKZB_D1,D_CRPKZB_D2,D_CRPKZB_DD1,D_CRPKZB_DD2
      CRPKZB=CNUM/CDEN
      D_CRPKZB_D1=(CDEN*D_CNUM_D1-CNUM*D_CDEN_D1)/CDEN**TWO
      D_CRPKZB_D2=(CDEN*D_CNUM_D2-CNUM*D_CDEN_D2)/CDEN**TWO
      D_CRPKZB_DD1=-TWO*DRU*CNUM*D_CDEN_DYA/CDEN**TWO
      D_CRPKZB_DD2=-TWO*DRD*CNUM*D_CDEN_DYB/CDEN**TWO
      D_CRPKZB_DD12=-TWO*DRT*CNUM*D_CDEN_DY/CDEN**TWO

      RETURN
      END SUBROUTINE COE_RPKZB_TPSS


!***********************************************************************
!
!> Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
!>
!> RU,RD                        density up,down
!> DRU, DRD                     abs. val. gradient of density up/down
!> DRT                          abs. val. gradient of total density
!> EPPGGA                       GGA ENERGY PER PARTICLE
!> EPPGGA_D1,EPPGGA_D2          THE DERIVATIVE OF EPPGGA WRT n
!> EPPGGA_DD1,EPPGGA_DD2        THE DERIVATIVE OF EPPGGA WRT |GRAD n|
!
!***********************************************************************
      SUBROUTINE CORGGA_TPSS(&
     &   RU,RD,DRU,DRD,DRT, &
     &   EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD12)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant

      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORGGA_TPSS)
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

      PARAMETER (GAMMA=0.03109069086965489503494086371273_q)
      PARAMETER (BETA_mb=0.06672455060314922_q)
GPACC !$ACC ROUTINE(CORPBE_tpss) SEQ

      EPPGGA  =0._q
      EPPGGA_D1=0._q
      EPPGGA_D2=0._q
      EPPGGA_DD12=0._q


      RT=RU+RD
      YA=DRU**2._q
      YB=DRD**2._q
      Y=DRT**2._q
!     YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2._q

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
      DZETAD1=TWO*RD/RT**TWO
      DZETAD2=-TWO*RU/RT**TWO


      DTHRD=exp(log(RT)*THRD)
      RS=(0.75_q/PI)**THRD/DTHRD
      DRSD1=-THRD/RT*RS
      DRSD2=-THRD/RT*RS


      PHI = (exp((TWO*THRD)*log(1._q+ZETA)) &
                +exp((TWO*THRD)*log(1._q-ZETA)))/2._q
      D_PHI_DZETA=THRD*((ONE+ZETA)**(-THRD)-(ONE-ZETA)**(-THRD))
      D_PHI_D1=D_PHI_DZETA*DZETAD1
      D_PHI_D2=D_PHI_DZETA*DZETAD2

      AFIX_T=SQRT(PI/FOUR)*(9._q*PI/FOUR)**(ONE/SIX)
      S=DRT/(TWO*(THREE*PISQ)**THRD*RT**THRD4)
      DSD1=-THRD4*S/RT
      DSD2=-THRD4*S/RT
      DSDD12 = ONE/(TWO*(THREE*PISQ)**THRD*RT**THRD4)

      T=AFIX_T*S/SQRT(RS)/PHI
      T2 = T*T
      T4 = T2*T2

      DTD1=AFIX_T*(PHI*RS*DSD1-0.5_q*S*PHI*DRSD1-S*RS*D_PHI_D1)/(PHI**TWO*RS**(THREE/TWO))
      DTD2=AFIX_T*(PHI*RS*DSD2-0.5_q*S*PHI*DRSD2-S*RS*D_PHI_D2)/(PHI**TWO*RS**(THREE/TWO))
      DTDD12=AFIX_T*DSDD12/(SQRT(RS)*PHI)

      FK=(3._q*PI*PI)**THRD*DTHRD
      SK = SQRT(4.0_q*FK/PI)
! Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
      CALL CORPBE_tpss(RS,ZETA,EC,VCUPLDA,VCDNLDA,PHI,SK, &
           T,H,DVCUP,DVCDN,ECQ,.FALSE.)
! OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
!      EPPGGA=EC
!      EPPGGA_D1=(VCUPLDA-EC)/RT
!      EPPGGA_D2=(VCDNLDA-EC)/RT
      D_EC_DZETA=(VCUPLDA-VCDNLDA)/TWO
      D_EC_DRS=(EC-ZETA*D_EC_DZETA-(VCUPLDA+VCDNLDA)/TWO)*THREE/RS

! THE RS DEPENDENCE OF BETA
      BETA = BETA_MB
      D_BETA_DRS=0._q


      PHI3=PHI**THREE
      PON=-EC/(PHI3*gamma)
      W=DEXP(PON)-ONE
      D_W_DRS=-(W+ONE)*D_EC_DRS/(GAMMA*PHI3)
      D_W_DZETA=-(W+ONE)/(GAMMA*PHI3)*(D_EC_DZETA-THREE*EC*D_PHI_DZETA/PHI)


      A=BETA/(GAMMA*W)
      D_A_DRS=(W*D_BETA_DRS-BETA*D_W_DRS)/(GAMMA*W**TWO)
      D_A_DZETA=-BETA*D_W_DZETA/(GAMMA*W**TWO)

      V=A*T2
      D_V_DRS=T2*D_A_DRS
      D_V_DZETA=T2*D_A_DZETA
      D_V_DT=TWO*A*T

      FUNKG=ONE/(ONE+V+V**TWO)
      D_FUNKG_DV=-(ONE+TWO*V)*FUNKG**TWO


      HCORE=ONE+W*(ONE-FUNKG)
      AH=GAMMA*PHI3
      H=AH*DLOG(HCORE)

      DH1=ONE-FUNKG
      DH2=W*D_FUNKG_DV
      D_H_DRS=(DH1*D_W_DRS-DH2*D_V_DRS)*AH/HCORE
      D_H_DZETA=THREE*H*D_PHI_DZETA/PHI+(DH1*D_W_DZETA-DH2*D_V_DZETA)*AH/HCORE
      D_H_DT=(-DH2*D_V_DT)*AH/HCORE

! OUTPUT EPPGGA AND ITS DERIVATIVES EPPGGA_D1,EPPGGA_D2, EPPGGA_DD1, EPPGGA_DD2
      EPPGGA=EC+H
      EPPGGA_D1=D_EC_DRS*DRSD1+D_EC_DZETA*DZETAD1+D_H_DRS*DRSD1+D_H_DZETA*DZETAD1+D_H_DT*DTD1
      EPPGGA_D2=D_EC_DRS*DRSD2+D_EC_DZETA*DZETAD2+D_H_DRS*DRSD2+D_H_DZETA*DZETAD2+D_H_DT*DTD2
      EPPGGA_DD12=D_H_DT*DTDD12

      RETURN
      END SUBROUTINE CORGGA_TPSS

! Modified by Adrienn Ruzsinszky
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE_tpss(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
GPACC !$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lmetagga=flag to do metagga (revTPSS)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORPBE_TPSS)
      logical lgga
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet_mb=0.06672455060314922_q)
      parameter(eta=1.e-12_q)
GPACC !$ACC ROUTINE(gcor2) SEQ
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
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
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
!write(*,*)'rs,VCUP,VCDN',rs,VCUP,VCDN
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
!  in TPSS beta is rs-independent
!      bet = bet_mb*(1._q + 0.1_q*RS)/(1._q + 0.1778_q*RS)
       bet = bet_mb

      delt=bet/gamma
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3._q
      GZ=(((1._q+zet)**2+eta)**sixthm- &
     &((1._q-zet)**2+eta)**sixthm)/3._q
      FAC = DELT/B+1._q
      BG = -3._q*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2._q*BET*G3*Q9/Q8
      COMM = H+HRS-7.0_q*T2*HT/6._q
      PREF = HZ-GZ*T2*HT/G
      COMM = COMM-PREF*ZET
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      ecdd=0.5_q/(sk*g)*t*ht
      RETURN
      END SUBROUTINE CORPBE_tpss

      SUBROUTINE GVT4(GVT,DG_DX,DG_DZ,XG,ZG,GAMA,CT,A,B,C,D,E,F)
GPACC !$ACC ROUTINE SEQ
       USE prec
       USE constant
       IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(GVT4)
!    SOME WORKING VARIABLES
       PARAMETER (F1= 1._q)
       PARAMETER (F2= 2._q)
       PARAMETER (F3= 3._q)


       G=GAMA
       G2=GAMA*GAMA

       GVT =(A + B*XG + C*ZG + D*XG*XG + E*ZG*XG + F*ZG*ZG)/G
       DG_DX =(-A*CT+B*(F1-F2*CT*XG)-F2*C*ZG*CT+D*(F2*XG-F3*XG*XG*CT) &
     &  +E*(ZG -F3*ZG*XG*CT)-F3*F*ZG*ZG*CT )/G2
       DG_DZ =(-A*CT -F2*B*XG*CT +C*(F1-F2*ZG*CT)-F3*D*XG*XG*CT &
     &  +E*(XG-F3*XG*ZG*CT)+F*(F2*ZG-F3*ZG*ZG*CT))/G2

       RETURN
       END SUBROUTINE GVT4

      SUBROUTINE CORGGA_1( &
     &   RU,RD,DRU,DRD,DRT, &
     &   EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL*8 (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORGGA_1)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (FIVE=5.0d0)
      PARAMETER (SIX=6.0d0)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)
!      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)

      PARAMETER (GAMMA=0.031090690869654895034940863712730d0)
      PARAMETER (BETA_mb=0.066724550603149220d0)
GPACC !$ACC ROUTINE(CORPBE_revtpss) SEQ

      EPPGGA  =0.0d0
      EPPGGA_D1=0.0d0
      EPPGGA_D2=0.0d0
      EPPGGA_DD=0.0d0

      RT=RU+RD
      YA=DRU**2.0d0
      YB=DRD**2.0d0
      Y=DRT**2.0d0
      YC=(Y-YA-YB)/2.0d0

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.99999999999990d0),0.99999999999990d0)
      DZETAD1=TWO*RD/RT**TWO
      DZETAD2=-TWO*RU/RT**TWO


      DTHRD=exp(log(RT)*THRD)
      RS=(0.750d0/PI)**THRD/DTHRD
      DRSD1=-THRD/RT*RS
      DRSD2=-THRD/RT*RS


      PHI = (exp((TWO*THRD)*log(1.0d0+ZETA)) &
     &       +exp((TWO*THRD)*log(1.0d0-ZETA)))/2.0d0
      D_PHI_DZETA=THRD*((1.0d0+ZETA)**(-THRD)-(1.0d0-ZETA)**(-THRD))
      D_PHI_D1=D_PHI_DZETA*DZETAD1
      D_PHI_D2=D_PHI_DZETA*DZETAD2

      AFIX_T=SQRT(PI/FOUR)*(9.0d0*PI/FOUR)**(1.0d0/SIX)
      S=DRT/(TWO*(THREE*PISQ)**THRD*RT**THRD4)
      DSD1=-THRD4*S/RT
      DSD2=-THRD4*S/RT
      DSDD=1.0D0/(TWO*RT*(THREE*PISQ*RT)**THRD)

      T=AFIX_T*S/SQRT(RS)/PHI
      T2 = T*T
      T4 = T2*T2

      DTD1=AFIX_T*(PHI*RS*DSD1-0.50d0*S*PHI*DRSD1-S*RS*D_PHI_D1) &
     &           /(PHI**TWO*RS**(THREE/TWO))
      DTD2=AFIX_T*(PHI*RS*DSD2-0.50d0*S*PHI*DRSD2-S*RS*D_PHI_D2) &
     &           /(PHI**TWO*RS**(THREE/TWO))
      DTDD=AFIX_T*DSDD/(SQRT(RS)*PHI)

      FK=(3.0d0*PI*PI)**THRD*DTHRD
      SK = SQRT(4.00d0*FK/PI)
      CALL CORPBE_revtpss(RS,ZETA,EC,VCUPLDA,VCDNLDA,PHI,SK, &
     &           T,H,DVCUP,DVCDN,ECQ,.FALSE.)
      D_EC_DZETA=(VCUPLDA-VCDNLDA)/TWO
      D_EC_DRS=(EC-ZETA*D_EC_DZETA-(VCUPLDA+VCDNLDA)/TWO)*THREE/RS

      AFACTOR=0.10d0
      BFACTOR=0.17780d0
      BETA_NUM=1.0d0 + AFACTOR*RS
      BETA_DEN=1.0d0+ BFACTOR*RS
      D_BETA_NUM=AFACTOR
      D_BETA_DEN=BFACTOR
      BETA = BETA_MB*BETA_NUM/BETA_DEN
      D_BETA_DRS=BETA_MB*(BETA_DEN*D_BETA_NUM-BETA_NUM*D_BETA_DEN) &
     &                  /BETA_DEN**TWO


      PHI3=PHI**THREE
      PON=-EC/(PHI3*gamma)
      W=DEXP(PON)-1.0d0
      D_W_DRS=-(W+1.0d0)*D_EC_DRS/(GAMMA*PHI3)
      D_W_DZETA=-(W+1.0d0)/(GAMMA*PHI3) &
     &                  *(D_EC_DZETA-THREE*EC*D_PHI_DZETA/PHI)


      A=BETA/(GAMMA*W)
      D_A_DRS=(W*D_BETA_DRS-BETA*D_W_DRS)/(GAMMA*W**TWO)
      D_A_DZETA=-BETA*D_W_DZETA/(GAMMA*W**TWO)

      V=A*T2
      D_V_DRS=T2*D_A_DRS
      D_V_DZETA=T2*D_A_DZETA
      D_V_DT=TWO*A*T


      FUNKG=1.0d0/(1.0d0+4.0d0*V)**0.250d0
      D_FUNKG_DV=-1.0d0/(1.0d0+4.0d0*V)**1.250d0



      HCORE=1.0d0+W*(1.0d0-FUNKG)
      AH=GAMMA*PHI3
      H=AH*DLOG(HCORE)

      DH1=1.0d0-FUNKG
      DH2=W*D_FUNKG_DV
      D_H_DRS=(DH1*D_W_DRS-DH2*D_V_DRS)*AH/HCORE
      D_H_DZETA=THREE*H*D_PHI_DZETA/PHI &
     &          +(DH1*D_W_DZETA-DH2*D_V_DZETA)*AH/HCORE
      D_H_DT=(-DH2*D_V_DT)*AH/HCORE

      EPPGGA=EC+H
      EPPGGA_D1=D_EC_DRS*DRSD1+D_EC_DZETA*DZETAD1+D_H_DRS*DRSD1 &
     &                        +D_H_DZETA*DZETAD1+D_H_DT*DTD1
      EPPGGA_D2=D_EC_DRS*DRSD2+D_EC_DZETA*DZETAD2+D_H_DRS*DRSD2 &
     &                        +D_H_DZETA*DZETAD2+D_H_DT*DTD2
      EPPGGA_DD=D_H_DT*DTDD

      RETURN
      END SUBROUTINE CORGGA_1

      SUBROUTINE CORGGA_0(RU,RD,DRU,DRD,DRT, &
     &   EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD)
GPACC !$ACC ROUTINE SEQ

      USE prec
      USE constant
      IMPLICIT REAL*8 (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(CORGGA_0)

      PARAMETER (CFB1=0.0285764d0)
      PARAMETER (CFB2=0.0889d0)
      PARAMETER (CFB3=0.125541d0)
      PARAMETER (cfkaiLD=0.12802585262625815d0)

      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (FIVE=5.0d0)
      PARAMETER (SIX=6.0d0)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)
!      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)

      PARAMETER (GAMMA=0.031090690869654895034940863712730d0)
      PARAMETER (BETA_mb=0.066724550603149220d0)

      EPPGGA  =0.0d0
      EPPGGA_D1=0.0d0
      EPPGGA_D2=0.0d0
      EPPGGA_DD=0.0d0
      AX_LDA=-THREE/(FOUR*PI)*(9.00d0*PI/FOUR)**THRD

      RT=RU+RD
      YA=DRU**2.0d0
      YB=DRD**2.0d0
      Y=DRT**2.0d0
      YC=(Y-YA-YB)/2.0d0

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.99999999999990d0),0.99999999999990d0)
      DZETAD1=TWO*RD/RT**TWO
      DZETAD2=-TWO*RU/RT**TWO


      DTHRD=exp(log(RT)*THRD)
      RS=(0.750d0/PI)**THRD/DTHRD
      DRSD1=-THRD/RT*RS
      DRSD2=-THRD/RT*RS


      PHI = (exp((TWO*THRD)*log(1.0d0+ZETA)) &
     &           +exp((TWO*THRD)*log(1.0d0-ZETA)))/2.0d0
      D_PHI_DZETA=THRD*((1.0d0+ZETA)**(-THRD)-(1.0d0-ZETA)**(-THRD))
      D_PHI_D1=D_PHI_DZETA*DZETAD1
      D_PHI_D2=D_PHI_DZETA*DZETAD2

      AFIX_T=SQRT(PI/FOUR)*(9.0d0*PI/FOUR)**(1.0d0/SIX)
      S=DRT/(TWO*(THREE*PISQ)**THRD*RT**THRD4)
      DSD1=-THRD4*S/RT
      DSD2=-THRD4*S/RT
      DSDD=1.0D0/(TWO*(THREE*PISQ)**THRD*RT**THRD4)

      RSHALF=SQRT(RS)
      FACTOR1=1.0d0+CFB2*RSHALF+CFB3*RS
      EC0LDA=-CFB1/FACTOR1

      D_EC0LDA_DRS=(CFB3+CFB2/RSHALF/TWO)*EC0LDA**TWO/cfb1

      dx_zeta=((1.0d0+zeta)**THRD4+(1.0d0-zeta)**THRD4)/two
      d_dx_zeta_dzeta=two*((1.0d0+zeta)**THRD-(1.0d0-zeta)**THRD)/three
      TWO13=TWO**THRD
!      gc_zeta=(two13-dx_zeta)/(two13-1.0d0)
!      d_gc_zeta_dzeta=-d_dx_zeta_dzeta/(two13-1.0d0)

      gc_zeta=(1.0d0-2.363d0*(dx_zeta-1.0d0))*(1.0d0-zeta**12.0d0)
      d_gc_zeta_dzeta=-(1.0d0-2.363d0*(dx_zeta-1.0d0))*12.d0*zeta**11.d0
      d_gc_zeta_dzeta=d_gc_zeta_dzeta -2.363d0*d_dx_zeta_dzeta*(1.0d0-zeta**12.0d0)


      w0=dexp(-EC0LDA/cfb1)-1.0d0
      d_w0_drs=-(w0+1.0d0)*D_EC0LDA_DRS/cfb1

      gfunkinf=1.0d0/(1.0d0+four*cfkaiLD*s*s)**(1.0d0/four)
      d_gfunkinf_ds=-two*cfkaiLD*s*gfunkinf**five

      hcore0=1.0d0+w0*(1.0d0-gfunkinf)
      h0=cfb1*dlog(hcore0)
      d_h0_drs=cfb1*d_w0_drs*(1.0d0-gfunkinf)/hcore0
      d_h0_ds=-cfb1*w0*d_gfunkinf_ds/hcore0


      EPPGGA=(EC0LDA+h0)*gc_zeta

      d_EPPGGA_dzeta=(EC0LDA+h0)*d_gc_zeta_dzeta
      d_EPPGGA_drs=(D_EC0LDA_DRS+d_h0_drs)*gc_zeta
      d_EPPGGA_ds=d_h0_ds*gc_zeta

      EPPGGA_D1=d_EPPGGA_dzeta*DZETAD1+d_EPPGGA_drs*DRSD1 &
     &          +d_EPPGGA_ds*DSD1
      EPPGGA_D2=d_EPPGGA_dzeta*DZETAD2+d_EPPGGA_drs*DRSD2 &
     &          +d_EPPGGA_ds*DSD2



      EPPGGA_DD=d_EPPGGA_ds*DSDD

      RETURN
      END SUBROUTINE CORGGA_0


!tb beg
      subroutine GGAapprox(Rho,nablaRho,b)
GPACC !$ACC ROUTINE SEQ
      USE prec
GPOMP !$OMP DECLARE TARGET(GGAapprox)
! evaluates a GGA approximation to b2
! Parameters for GGA-approx
       REAL(q) s,b,rs ! reduced density gradient, b, Wigner-Seitz radius,
       REAL(q) Rho,nablaRho
       REAL(q),Parameter :: pi=3.141592653589793_q,P1=2._q,P2=1._q
       REAL(q),Parameter :: DCS=(1.D0/(2.D0*((3.0_q*(PI**2))**(1._q/3._q))))
       if((RHO.gt.1.E-10_q).and.(nablaRho.gt.1.E-10_q))then
       RS=(3._q/(4._q*pi*Rho))**(1._q/3._q) !Wigner-Seitz radius
       S=(DCS*nablaRho)/(Rho**(4._q/3._q))

!       P1=2.0 !0.97 !1.0 !2.0
!       P2=1.0 !0.57 !0.6057068643 !1.0
       b=P1*s*exp(-P2*s) !ExpExp
       b=b*rs
       b=b*b ! Our formula approximates b, but we would like to return b^2
!       b=b*rho ! this is the quantity we need for M2, so let's compute it right away
!      However, it is preferable to multiply by the total density, i.e., including the core contribution
!      Hence, we do this step in the xxx_planeDOM routine
       else
!       write(*,*) RHO,nablaRho
       b=0._q
       endif
!       b=rho
       return
       end subroutine GGAapprox

!
!  Common interface to all gradient correction routines:
!  returns energy + potential (if defined)
!
      SUBROUTINE GGAEALL(RHO,DRHO,DABGRH,DDRHO,VXC,EXC)
GPACC !$ACC ROUTINE SEQ
      USE prec
      USE tutor, ONLY: vtutor
      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET(GGAEALL)

__CPU CALL vtutor%bug('internal ERROR GGAEALL: Wrong functional, scheme not implemented!', __FILE__, __LINE__)

      RETURN
      END SUBROUTINE GGAEALL

END MODULE GGALIB
