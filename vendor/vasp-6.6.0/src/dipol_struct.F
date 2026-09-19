#include "symbol.inc"

  MODULE mdipol_struct
    USE prec
    IMPLICIT NONE
    TYPE dipol
!only DIP
      INTEGER IDIPCO                 !< direction (0 no dipol corrections)
        LOGICAL :: LCOR_DIP          !< correct potential
        LOGICAL :: LMONO             !< monopole corrections only
        REAL(q) :: EPSILON           !< bulk charged cells: dielectric constant
        REAL(q) :: POSCEN(3)         !< position of center
        REAL(q) :: DIPOLC(3)         !< calculated dipol
        REAL(q) :: QUAD              !< trace of quadrupole
        INTEGER :: INDMIN(3)         !< position of minimum
        REAL(q) :: ECORR,EMONO,EDIPOL,EQUAD
        REAL(q) :: E_ION_EXTERN      !< energy between ions and added field
        REAL(q),POINTER :: FORCE(:,:)
        REAL(q) :: VACUUM(2)         !< vacuum level
        REAL(q) :: EFIELD
        LOGICAL :: VACPOTAV          !< If T, the averaged vacuum potential is reported 
        REAL(q) :: VACPOTFLAT        !< Expected max field in vacuum to determine a vacuum potential 
        REAL(q) :: DENMIN            !< Minimum density averaged along the planar direction 
        REAL(q),POINTER :: POTLIN(:)          !< 1-D average of the potential
        REAL(q),POINTER :: DIST_ALONG_IDIPOL(:)  !< Distance along the IDIPOL direction
    END TYPE

    TYPE (dipol)    DIP
    LOGICAL,SAVE :: LDIPOL_RESET=.TRUE.

  END MODULE mdipol_struct
