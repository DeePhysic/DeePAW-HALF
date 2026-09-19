#include "symbol.inc"
!***********************************************************************
!> These are stubs for the Yu/Libisch/Carter PAW DFET embedding routines
!***********************************************************************
module mextpot

    implicit none

contains

   logical function EXTPT_LEXTPOT()
      EXTPT_LEXTPOT = .FALSE.
   end function EXTPT_LEXTPOT

   logical function EXTPT_DEXTPOT()
      EXTPT_DEXTPOT = .FALSE.
   end function EXTPT_DEXTPOT

   subroutine EXTPT_READER(NGXF,NGYF,NGZF,IU0,IU5,IU6)
      integer :: NGXF, NGYF, NGZF, IU0, IU5, IU6
   end subroutine EXTPT_READER

   subroutine EXTPT_EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT)
      USE prec
      USE lattice
      USE mgrid
      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) ::      CVTOT(GRIDC%MPLWV)
   end subroutine EXTPT_EXTERNAL_POT_ADD

   subroutine EXTPT_EXTERNAL_POT_ADD_PAW(POTAE, POT, NDIM, NCDIJ, LMMAX, NIP)
      USE prec
      INTEGER             NDIM, NCDIJ, LMMAX, NIP
      REAL(q) POT(:,:,:), POTAE(:,:,:)
   end subroutine EXTPT_EXTERNAL_POT_ADD_PAW

   subroutine EXTPT_CALC_VLM(LATT_CUR, GRIDC, T_INFO, P, LMDIM, WDES, IU0)
      use lattice
      use mgrid
      use poscar
      use pseudo
      use wave
      type (grid_3d)    GRIDC
      type (latt)       LATT_CUR
      type (type_info)  T_INFO
      type (potcar), TARGET :: P(T_INFO%NTYP)
      type (wavedes)    WDES
      integer :: LMDIM, IU0
   end subroutine EXTPT_CALC_VLM

   subroutine EXTPT_DE_DVRN(WDES, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, &
        P, T_INFO, LMDIM, CRHODE, CHTOT_, CHDEN, IRDMAX )
      use prec
      use mgrid
      use lattice
      use poscar
      use pseudo
      use wave
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      TYPE (grid_3d)     GRID_SOFT
      TYPE (latt)        LATT_CUR
      TYPE (transit)     SOFT_TO_C
      TYPE (wavedes)     WDES
      integer :: LMDIM, IRDMAX
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      OVERLAP   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q), TARGET :: CHTOT_(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)
   end subroutine EXTPT_DE_DVRN

end module mextpot
