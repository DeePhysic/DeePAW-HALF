#include "symbol.inc"
#ifdef VASP_DEBUG
MODULE debug_struct

   IMPLICIT NONE

   TYPE debug_settings
      !> Time step at which to request early exit of ionic loop.
      INTEGER :: REQUEST_EXIT_AT_STEP = -1
      !> MPI rank which requests early exit of ionic loop (if negative, all ranks request exit).
      INTEGER :: REQUEST_EXIT_ON_RANK = -1
      !> Time step at which to write a STOPCAR file.
      INTEGER :: WRITE_STOPCAR_AT_STEP = -1
   END TYPE

END MODULE debug_struct
#endif
