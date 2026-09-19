#include "symbol.inc"
#ifdef OFFLOADING
    MODULE moffload_struct_def

GPACC USE mopenacc_struct_def
GPOMP USE mopenmp_struct_def

      !> Flags whether device-sided execution is allowed in principle
      LOGICAL, PUBLIC :: OFFLOAD_ALLOW=.FALSE.
      !> Flags whether device-sided execution is momentarily activated
      LOGICAL, PUBLIC :: OFFLOAD_ON=.FALSE.

      !> Total number of devices detected by the OpenACC runtime
      INTEGER, PUBLIC :: OFFLOAD_NUM_DEVICES
      !> Maximum number of MPI-ranks per device
      INTEGER, PUBLIC :: OFFLOAD_MAX_RANKS_PER_DEVICE

#ifdef USENCCL
      !> Flags whether NCCL is used
      LOGICAL, PUBLIC :: LUSENCCL=.TRUE.
#endif
    END MODULE moffload_struct_def

!***********************************************************************
!
!***********************************************************************

    MODULE moffload_control

      USE moffload_struct_def, ONLY : OFFLOAD_ON, OFFLOAD_ALLOW

      PRIVATE

      PUBLIC :: OFFLOAD_ON_PUSH, OFFLOAD_ON_POP

      INTEGER, PARAMETER :: MAXLEVEL=20
      INTEGER :: OFFLOAD_ON_LEVEL=0
      LOGICAL :: OFFLOAD_ON_STACK(MAXLEVEL)=.FALSE.

      CONTAINS

!****************** SUBROUTINE OFFLOAD_ON_PUSH *************************
!
!***********************************************************************

      SUBROUTINE OFFLOAD_ON_PUSH(VAR, LWAIT)
      USE tutor, ONLY : vtutor
      LOGICAL :: VAR
      LOGICAL, OPTIONAL :: LWAIT
      ! local variable
      LOGICAL :: LWAIT_ = .TRUE.
      IF (PRESENT(LWAIT)) LWAIT_ = LWAIT
!$ACC WAIT IF(LWAIT_)
      IF (OFFLOAD_ON_LEVEL==MAXLEVEL) THEN
         CALL vtutor%error("OFFLOAD_ON_PUSH: ERROR: stack is full")
      ENDIF
      OFFLOAD_ON_LEVEL=OFFLOAD_ON_LEVEL+1
      OFFLOAD_ON_STACK(OFFLOAD_ON_LEVEL)=OFFLOAD_ON
      OFFLOAD_ON=VAR.AND.OFFLOAD_ALLOW
      END SUBROUTINE OFFLOAD_ON_PUSH


!****************** SUBROUTINE OFFLOAD_ON_POP **************************
!
!***********************************************************************

      SUBROUTINE OFFLOAD_ON_POP(LWAIT)
      USE tutor, ONLY : vtutor
      LOGICAL, OPTIONAL :: LWAIT
      ! local variable
      LOGICAL :: LWAIT_ = .TRUE.
      IF (PRESENT(LWAIT)) LWAIT_ = LWAIT
!$ACC WAIT IF(LWAIT_)
      IF (OFFLOAD_ON_LEVEL==0) THEN
         CALL vtutor%error("OFFLOAD_ON_POP: ERROR: stack is empty")
      ENDIF
      OFFLOAD_ON=OFFLOAD_ON_STACK(OFFLOAD_ON_LEVEL)
      OFFLOAD_ON_LEVEL=OFFLOAD_ON_LEVEL-1
      END SUBROUTINE OFFLOAD_ON_POP

    END MODULE moffload_control
#endif // OFFLOADING
