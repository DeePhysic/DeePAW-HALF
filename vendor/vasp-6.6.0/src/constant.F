#include "symbol.inc"
!************************************************************************
!> Physical constants and unit conversion used in VASP
!>
!>   dimension                  VASP units (v.u.)         atomic units (a.u.)
!> + length [L]                    Angstroem (A)                Bohr
!> + time   [T]                      eV^-1                   Hartree^-1
!> + charge [Q]                       1                          1
!> + mass   [M]                    eV^-1 A^-2                    1
!> + energy [ML^2T^-2]                eV                   Hartree (2*Rydberg)
!> + action [ML^2T^-1]                1                          1
!>
!> Fundamental constants
!>                                   v.u.                        a.u.
!> reduced Planck constant
!>   hbar                             1                          1
!>
!> Charge of an electron
!>   e                                1                          1
!>
!> Bohr radius
!>   a_B                           0.529177249 [A]               1 [Bohr]
!>                                  =AUTOA
!>   a_B = 4\pi \epsilon_0 hbar^2/(m_e e^2)
!>
!> Hartree energy
!>   E_h                           27.211652   [eV]              1 [Hartree]
!>                                  =2*RYTOEV
!>   E_h = hbar^2/(m_e a_B^2)
!>
!> Mass of an electron
!>   m_e                       0.1312329113 [eV^-1 A^-2]         1 [Hartree^-1 Bohr^-2]
!>                              = 1/(2*RYTOEV*AUTOA*AUTOA)
!>
!> Permitivity of free space times 4 pi
!>   4\pi \epsilon_0           0.06944547095 [eV^-1 A^-1]        1 [Hartree^-1 Bohr^-1]
!>                              = 1/(2*RYTOEV*AUTOA)
!>                              = 1/FELECT
!>   4\pi \epsilon_0 = e^2/(a_B E_h)
!>   F = eE = e^2/(4\pi \epsilon_0) * 1/(r-r')^2 = (e^2 k_e)/(r-r')^2
!>
!> Coulomb constant
!>   k_e                         14.3997871461 [eV A]            1 [Hartree Bohr]
!>                               = FELECT
!>   k_e = 1/(4\pi \epsilon_0)
!>
!> Speed of light
!>   c                      1973.3036311408 [eV A]           | 137.037 [Hartree Bohr]
!>                              = CLIGHT*2*RYTOEV*AUTOA      |    = CLIGHT
!>   c = (a_B E_h)/(hbar alpha)
!>
!> Fine-structure constant (dimensionless, independent of system of units)
!>   alpha                      1/137.037                    | 1/137.037
!>                              = 1/CLIGHT                   |    = 1/CLIGHT
!>
!> Bohr magneton
!>   mu_B=                     3.8100198727 [eV A^2]            1/2 [Hartree Bohr^2]
!>                              = RYTOEV*AUTOA*AUTOA
!>                              = HSQDTM
!>   mu_B=(e hbar)/(2 m_e)
!>   NB: Implies B-field in v.u. is given in [A^-2]!
!>       B = nabla x A implies the gauge field in v.u. is given in [A^-1]!
!>
!> Vacuum magnetic permeability divided by 4pi
!>   mu_0/(4\pi)=               = 3.698011373E-6 [eV^-1 A^-1]   5.325057663E-5 [Hartree^-1 Bohr^-1]
!>                              = FELECT/CLIGHTV**2                = 1/CLIGHT**2 [Hartree^-1 Bohr^-1]
!>   mu_0/(4\pi)=1/(c^2 4\pi \epsilon_0)
!>   NB: Biot-Savart law B = mu_0/(4\pi) int jx(r'-r)/|r'-r|^3 d^3r.
!>
!> Kinetic-energy factor
!>   hbar^2 /(2 m_e)           3.8100198727 [eV A^2]            1/2 [Hartree Bohr^2]
!>                              = RYTOEV*AUTOA*AUTOA
!>                              = HSQDTM
!>
!> Charge density              CHTOT [A^{-3}]                  CHTOT * AUTOA3 [Bohr^{-3}]
!>                            where Psi are in v.u.
!>
!> Kinetic energy density
!>                             TAU [eV A^{-3}]                 TAU * TAUTOAU [Hartree Bohr^{-3}]
!>
!> Paramagnetic current density
!>                             JDEN [eV Angstroem^{-2}]        JDEN * JDENTOAU [Hartree Bohr^{-2}]
!>   NB: Minimal substitution reads: p -> p - eA and yields E = int jA d^3r.
!>       The diamagnetic current reads j_dia = psi* e^2/m A psi.
!>       The paramagnetic current reads j_para = 0.5 (nk + kn) hbar/m.
!>
!> Boltzmann constant
!>   k_B                       8.6173857E-5 [eV K^-1]
!>                              = BOLKEV
!>
!***********************************************************************
MODULE constant
      USE prec

      ! unitless constants
      REAL(q),   PARAMETER  :: PI = 3.141592653589793238_q !< pi
      REAL(q),   PARAMETER  :: TPI = 2*PI !< 2 pi
      REAL(q),   PARAMETER  :: FPI2 = TPI*TPI !< 4 pi^2
      COMPLEX(q),PARAMETER  :: CITPI = (0._q,1._q)*TPI !< 2 i pi

      ! a.u. conversion to v.u.
      REAL(q), PARAMETER :: AUTOA=0.529177249_q !< 1 a.u. of length (Bohr) in Angstroem
      REAL(q), PARAMETER :: RYTOEV=13.605826_q  !< 1 Rydberg in eV
      REAL(q), PARAMETER :: CLIGHT=137.037_q    !< inverse fine-structure constant, speed of light in a.u.

      REAL(q), PARAMETER :: AUTOA2=AUTOA *AUTOA
      REAL(q), PARAMETER :: AUTOA3=AUTOA2*AUTOA
      REAL(q), PARAMETER :: AUTOA4=AUTOA2*AUTOA2
      REAL(q), PARAMETER :: AUTOA5=AUTOA3*AUTOA2

      ! constants in v.u.
      ! e = hbar = 1
      REAL(q), PARAMETER :: BOLKEV = 8.6173857E-5_q !< Boltzmanns constant k_B in eV/K
      REAL(q), PARAMETER :: FELECT = 2*RYTOEV*AUTOA !< Coulomb constant k_e = 1/(4 pi epsilon_0) in eV Angstroem
      REAL(q), PARAMETER :: EDEPS = 4*PI*FELECT     !< inverse permittivity of free space (1/epsilon_0) in eV Angstroem
      REAL(q), PARAMETER :: HSQDTM = RYTOEV*AUTOA*AUTOA !< 1/(2m) in eV Angstroem^2 = hbar^2/2m
      REAL(q), PARAMETER :: CLIGHTV=FELECT*CLIGHT   !< speed of light in v.u. (eV Angstroem)

      ! v.u. conversion to a.u.
      REAL(q), PARAMETER :: TAUTOAU=0.5_q*AUTOA3/RYTOEV   !< 1 v.u. kinetic energy density in Hartree Bohr^-3
      REAL(q), PARAMETER :: JDENTOAU=0.5_q*AUTOA2/RYTOEV  !< 1 v.u. current density in (hartree bohr^{-2})

      ! conversion to SI
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q !< 1 eV in Joule
      REAL(q), PARAMETER :: AMTOKG=1.6605402E-27_q !< 1 atomic mass unit ("proton mass") in kg
      REAL(q), PARAMETER :: BOLK=BOLKEV*EVTOJ      !< Boltzmanns constant in Joule/K

      REAL(q), PARAMETER :: EVTOKCAL=23.06_q

!> vector field A times momentum times e/ (2 m_e c) is an energy
!> magnetic moments are supplied in Bohr magnetons
!>~~~
!> e / (2 m_e c) A(r) p(r) = energy
!> e / (2 m_e c) m_s x ( r - r_s) / (r-r_s)^3 hbar nabla =
!> e^2 hbar^2 / (2 m_e^2 c^2) 1/ lenght^3 = energy
!>~~~
!> conversion factor from magnetic moment to energy
!> checked independently in SI by Gilles de Wijs
      REAL(q), PARAMETER :: MAGMOMTOENERGY=1/CLIGHT**2*AUTOA**3*RYTOEV

!> dimensionless number connecting input and output magnetic moments
!> AUTOA e^2 (2 m_e c^2)
      REAL(q), PARAMETER :: MOMTOMOM=AUTOA/CLIGHT/CLIGHT/2
      
      ENUM, BIND(c)
         ENUMERATOR  :: UNIT_INTERNAL
         ENUMERATOR  :: UNIT_EV
         ENUMERATOR  :: UNIT_HZ
         ENUMERATOR  :: UNIT_2PI_HZ
         ENUMERATOR  :: UNIT_INV_CM
      END ENUM

      CONTAINS

      SUBROUTINE NUP_CONST()
      END SUBROUTINE

!*********************************************************************
!
! Converts frequencies (FREQU) from one physical unit (FROM) to
! another (TO).
!
! Internal unit = 2 pi eV^{1/2} Angstrom^{-1} atomicmass^{-1/2}
!
! An exponent can be supplied that allows for the conversion of, e.g.,
! squared (EXPONENT = 2) or inverse (EXPONENT = -1) frequencies.
! If the initial and final units are the same, EXPONENT is ignored.
!
!*********************************************************************
  PURE ELEMENTAL FUNCTION CONVERT_FREQUENCY(FREQU, FROM, TO, EXPONENT) RESULT(OUT)

    IMPLICIT NONE

    REAL(q)                                     :: OUT
    REAL(q),                        INTENT(IN)  :: FREQU
    INTEGER(KIND(UNIT_INTERNAL)),   INTENT(IN)  :: FROM
    INTEGER(KIND(UNIT_INTERNAL)),   INTENT(IN)  :: TO
    REAL(q), OPTIONAL,              INTENT(IN)  :: EXPONENT

    REAL(q), PARAMETER  :: PLANCK_SI = 6.626070040E-34_q
    REAL(q), PARAMETER  :: C_SI = 299792458_q
    REAL(q), PARAMETER  :: SI_CONV = SQRT(EVTOJ / AMTOKG) * 1E10_q

    REAL(q) :: FACTOR

    OUT = FREQU
    IF (FROM == TO) RETURN

    ! Convert frequency to internal units.
    SELECT CASE(FROM)
    CASE(UNIT_INTERNAL)
      FACTOR = 1._q
    CASE(UNIT_EV)
      FACTOR = 2._q * PI * EVTOJ / PLANCK_SI
    CASE(UNIT_HZ)
      FACTOR = 2._q * PI
    CASE(UNIT_2PI_HZ) ! Happens to be internal unit * SI_CONV
      FACTOR = 1._q
    CASE(UNIT_INV_CM)
      FACTOR = 2._q * PI * C_SI * 1E2_q
    END SELECT

    ! Convert frequency to desired units.
    SELECT CASE(TO)
    CASE(UNIT_INTERNAL)
    CASE(UNIT_EV)
      FACTOR = FACTOR / (2._q * PI * EVTOJ / PLANCK_SI)
    CASE(UNIT_HZ)
      FACTOR = FACTOR / (2._q * PI)
    CASE(UNIT_2PI_HZ) ! Happens to be internal unit * SI_CONV
    CASE(UNIT_INV_CM)
      FACTOR = FACTOR / (2._q * PI * C_SI * 1E2_q)
    END SELECT

    ! Apply common prefactor.
    IF (FROM == UNIT_INTERNAL .AND. TO /= UNIT_INTERNAL) THEN
      FACTOR = FACTOR * SI_CONV
    ELSE IF (FROM /= UNIT_INTERNAL .AND. TO == UNIT_INTERNAL) THEN
      FACTOR = FACTOR / SI_CONV
    ENDIF

    ! Apply exponent.
    ! NOTE: This only works if the unit scales have a common zero point.
    IF (PRESENT(EXPONENT)) FACTOR = FACTOR ** EXPONENT

    OUT = OUT * FACTOR

  END FUNCTION CONVERT_FREQUENCY

END MODULE
