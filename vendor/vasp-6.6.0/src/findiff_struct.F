#include "symbol.inc"

MODULE findiff_struct_def
   USE prec
   IMPLICIT NONE

   !> @brief Structure to track finite-difference progress
   !>
   !> This structure keeps track of the progress of finite-difference
   !> calculations by maintaining counters for displacement steps and
   !> degrees of freedom (DOF). It also provides methods to compute the
   !> total number of steps and the current index in the overall process.
   !>
   !> I don't think this is a very nice way to do this, but it's what we have for now.
   !> Ideally, one would have a better way to perform finite differences using
   !> a stencil or similar approach.
   TYPE findiff_progress_t
      INTEGER :: DISPL_COUNT = 0  !< Total number of displacement steps per DOF
      INTEGER :: DOF_COUNT = 0    !< Total number of DOF
      INTEGER :: DISPL_INDEX = 0  !< Current displacement step
      INTEGER :: DOF_INDEX = 0    !< Current DOF
      REAL(q) :: STEP = 0._q      !< Current magnitude of displacement
      INTEGER :: ATOM_INDEX = 0   !< Current atom index
      INTEGER :: CART_INDEX = 0   !< Current Cartesian direction of displacement
   CONTAINS
      PROCEDURE :: TOTAL_COUNT => FINDIFF_PROGRESS_TOTAL_COUNT
      PROCEDURE :: TOTAL_INDEX => FINDIFF_PROGRESS_TOTAL_INDEX
   END TYPE findiff_progress_t

CONTAINS

   PURE INTEGER FUNCTION FINDIFF_PROGRESS_TOTAL_COUNT(THIS) RESULT(TOTAL)
      CLASS(findiff_progress_t), INTENT(IN) :: THIS
      TOTAL = THIS%DOF_COUNT * THIS%DISPL_COUNT
   END FUNCTION

   PURE INTEGER FUNCTION FINDIFF_PROGRESS_TOTAL_INDEX(THIS) RESULT(IDX)
      CLASS(findiff_progress_t), INTENT(IN) :: THIS
      IDX = (THIS%DOF_INDEX - 1) * THIS%DISPL_COUNT + THIS%DISPL_INDEX
   END FUNCTION

END MODULE findiff_struct_def
