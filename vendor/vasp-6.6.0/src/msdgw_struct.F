#include "symbol.inc"

!********************************************************************************
!> contains structure definitions for msdgw module
!********************************************************************************

MODULE msdgw_struct
   USE prec
   USE base
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: band_ebin

   TYPE band_ebin
      !> indices of CELTOT that go into this bin
      INTEGER, ALLOCATABLE :: IBAND(:)
      !> energy average for this bin
      REAL(q):: EBAR = 0 
      !> energy ratio for this bin 
      REAL(q):: F = 0 
      !> number of states in this bin
      INTEGER :: NBANDS = 0 
      !> number of pseudo states  
      INTEGER :: NXI = 0 
   END TYPE band_ebin

   CONTAINS 

END MODULE msdgw_struct

