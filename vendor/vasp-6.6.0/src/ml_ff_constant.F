#include "symbol.inc"
!************************************************************************
!
!  this module contains some control data structures for
!  machine-learning force field.
!
!***********************************************************************
      MODULE ML_FF_CONSTANT
      USE ML_FF_PREC
! Constants used for machine learning force field

!  Some important Parameters, to convert to a.u.
!  - AUTOA  = 1. a.u. in Angstroem
!  - RYTOEV = 1 Ry in Ev
!  - EVTOJ  = 1 eV in Joule
!  - AMTOKG = 1 atomic mass unit ("proton mass") in kg

      REAL(q), PARAMETER :: AUTOA=0.529177249_q,RYTOEV=13.605826_q
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q
      REAL(q), PARAMETER :: PI =3.141592653589793238_q

! EUNIT = 1 Hartree in ev units
! MUNIT = Electron mass in atomic units
! RUNIT = Bohr Radius in Angstrom1
! FUNIT = Unit of force in Hartree/Bohr
! SUNIT = Unit of stress tensor
! TUNIT = Unit of time in atomic unit

      REAL(q), PARAMETER :: EUNIT = 2.0_q*RYTOEV, &
                            MUNIT = 1.054571726E-034_q**2/2.0_q/(RYTOEV*EVTOJ*AUTOA**2*1.0E-20_q)/AMTOKG, &
                            FUNIT = 2.0_q*RYTOEV/AUTOA, &
                            SUNIT = 2.0_q*RYTOEV/(AUTOA*1.0E-10_q)**3*EVTOJ/1.0E8_q, &
                            TUNIT = 1.054571726E-034_q/(2.0_q*RYTOEV*EVTOJ)*1.0E15_q

! Parameters for Normalization (to avoid it if norm becomes zero)

      REAL(q), PARAMETER :: TOLERANCE_NORM = 1.0E-08_q

      INTEGER, PARAMETER :: NBLOCK_SCALAPACK = 32
      !> Length in bytes of the ML_FF ASCII header.
      INTEGER, PARAMETER :: ML_FF_HEADER_SIZE = 4096
      !======================================================================
      !> ML_LOGFILE version number.
      !>
      !> Version history:
      !> * 0.1.0: Initial version.
      !> * 0.2.0: Minor updates, changes include:
      !>   - Additional STATUS category "predfast".
      !>   - New tags output:
      !>     ML_MODE, ML_LFAST, ML_LSIC, ML_EPS_REG, ML_OUTPUT_MODE, ML_OUTBLOCK,
      !>     ML_LSUPERVEC, ML_LBASIS_DISCARD.
      !>   - New info lines NDESC and NDESC_SIC.
      !> * 0.2.1: Minor updates, changes:
      !>   - New tag output: ML_DESC_TYPE.
      !>   - Maximum buffer size below tag ML_MB added.
      !> * 0.2.2: Minor updates, changes:
      !>   - New tag output: ML_MB_MIN.
      !>   - New info lines: MSG.
      !>   - Skipping per-step log lines in fast prediction mode.
      !> * 0.2.3: Minor updates, changes:
      !>   - New log line: SF(F).
      !>   - New tag output: ML_IERR, ML_CALGO, ML_LIB.
      !> * 0.2.4: Minor change:
      !>   - Renamed tag: ML_IERR -> ML_ESTBLOCK.
      !> * 0.2.5: Minor update:
      !>   - Show whether ML_FF supports fast prediction mode.
      !>   - New log lines: NORME, ERPS, NEPS, LASTNE, BE(F)PS, SPF(F)PS.
      !>   - New tag output:
      !>     - ML_SAVECMAT, ML_SF_LIMIT, ML_LEMPPOT.
      !>     - ML_PALGO, ML_TYPE, ML_GRACE_MODEL (VASPml only)
      !>   - If ML_LEMPPOT = .TRUE. new tag output:
      !>     - ML_EMPPOT_RCUT,
      !>     - ML_SRPOT_(B0,S0,N0)
      !>     - ML_MOPOT_(NM,DM,RM,RKM,IJM).
      !>   - Fixed SPRSC formatting.
      !======================================================================
      INTEGER, PARAMETER :: ML_LOGFILE_VERSION(3) = [0, 2, 5]
      !======================================================================
      !> ML_FF version number.
      !>
      !> Version history:
      !> * 0.0.1: Original ASCII file.
      !> * 0.1.0: Binary version of previous version.
      !> * 0.2.0: Added ASCII header to binary file.
      !> * 0.2.1: Added ML_DESC_TYPE=1 output.
      !> * 0.2.2: Added ML_DESC_TYPE=2,3 output.
      !> * 0.2.3: Added inverse kernel matrix SINV.
      !> * 0.2.4: Added LFAST field in ML_FF file.
      !> * 0.2.5: Added WRITECMAT variable to ML_FF.
      !======================================================================
      INTEGER, PARAMETER :: ML_FF_VERSION_WRITE(3) = [0, 2, 5]

      END MODULE
