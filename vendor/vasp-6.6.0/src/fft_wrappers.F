#include "symbol.inc"

!************************* SUBROUTINE WFFTBAS **************************
!
!>
!
!***********************************************************************
#undef WFFTBAS
      SUBROUTINE WFFTBAS(C, GRID, ISN)

      USE prec
      USE mgrid_struct_def
__GPU USE moffload_struct_def
__GPU USE moffload_fft, ONLY : FFTBAS_PLAN_OFFLOAD
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISN

__GPU IF (IS_ON_DEVICE(C) .AND. OFFLOAD_ON) THEN
__GPU    CALL FFTBAS_PLAN_OFFLOAD(C, GRID, ISN)
__GPU    RETURN
__GPU ENDIF

      CALL FFTBAS_PLAN(C, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBAS


!************************* SUBROUTINE WFFTBRC **************************
!
!>
!
!***********************************************************************
#undef WFFTBRC
      SUBROUTINE WFFTBRC(C, GRID, ISN)

      USE prec
      USE mgrid_struct_def
__GPU USE moffload_struct_def
__GPU USE moffload_fft, ONLY : FFTBRC_PLAN_OFFLOAD

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: ISN

__GPU IF (IS_ON_DEVICE(C) .AND. OFFLOAD_ON) THEN
__GPU    CALL FFTBRC_PLAN_OFFLOAD(C, GRID, ISN)
__GPU    RETURN
__GPU ENDIF

      CALL FFTBRC_PLAN(C, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBRC


!************************* SUBROUTINE WFFTBAS_MU ***********************
!
!>
!
!***********************************************************************
#undef WFFTBAS_MU
      SUBROUTINE WFFTBAS_MU(N, C, LDC, GRID, ISN)

      USE prec
      USE mgrid_struct_def
__GPU USE moffload_struct_def
__GPU USE moffload_fft, ONLY : FFTBAS_PLAN_MU_OFFLOAD

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N, LDC, ISN

__GPU IF (IS_ON_DEVICE(C) .AND. OFFLOAD_ON) THEN
__GPU    CALL FFTBAS_PLAN_MU_OFFLOAD(N, C, LDC, GRID, ISN)
__GPU    RETURN
__GPU ENDIF

      CALL FFTBAS_PLAN_MU(N, C, LDC, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBAS_MU


!************************* SUBROUTINE WFFTBRC_MU ***********************
!
!>
!
!***********************************************************************
#undef WFFTBRC_MU
      SUBROUTINE WFFTBRC_MU(N, C, LDC, GRID, ISN)

      USE prec
      USE mgrid_struct_def
__GPU USE moffload_struct_def
__GPU USE moffload_fft, ONLY : FFTBRC_PLAN_MU_OFFLOAD

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID
      COMPLEX(q) :: C(*)
      INTEGER :: N, LDC, ISN

__GPU IF (IS_ON_DEVICE(C) .AND. OFFLOAD_ON) THEN
__GPU    CALL FFTBRC_PLAN_MU_OFFLOAD(N, C, LDC, GRID, ISN)
__GPU    RETURN
__GPU ENDIF

      CALL FFTBRC_PLAN_MU(N, C, LDC, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBRC_MU

#ifdef MPI
!************************* SUBROUTINE WFFTBAS_MPI **********************
!
!>
!
!***********************************************************************
#undef WFFTBAS_MPI
      SUBROUTINE WFFTBAS_MPI(A, GRID, ISN)

      USE prec
      USE mgrid_struct_def
__GPU USE moffload_struct_def
__GPU USE moffload_fft, ONLY : FFTBAS_PLAN_MPI_OFFLOAD

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: A(*)
      INTEGER :: ISN

__GPU IF (IS_ON_DEVICE(A) .AND. OFFLOAD_ON) THEN
__GPU    CALL FFTBAS_PLAN_MPI_OFFLOAD(A, GRID, ISN)
__GPU    RETURN
__GPU ENDIF

      CALL FFTBAS_PLAN_MPI(A, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBAS_MPI


!************************* SUBROUTINE WFFTBRC_MPI **********************
!
!>
!
!***********************************************************************
#undef WFFTBRC_MPI
      SUBROUTINE WFFTBRC_MPI(A, GRID, ISN)

      USE prec
      USE mgrid_struct_def
__GPU USE moffload_struct_def
__GPU USE moffload_fft, ONLY : FFTBRC_PLAN_MPI_OFFLOAD

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: A(*)
      INTEGER :: ISN

__GPU IF (IS_ON_DEVICE(A) .AND. OFFLOAD_ON) THEN
__GPU    CALL FFTBRC_PLAN_MPI_OFFLOAD(A, GRID, ISN)
__GPU    RETURN
__GPU ENDIF

      CALL FFTBRC_PLAN_MPI(A, GRID, ISN)

      RETURN
      END SUBROUTINE WFFTBRC_MPI
#endif // MPI
