#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
!****************************************************************************************************
! Module for some mpi routines that are needed very early
!****************************************************************************************************

      MODULE MPI_HELP

        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        IMPLICIT NONE

        INCLUDE "mpif.h"
        TYPE ML_MPI_PAR 
            INTEGER :: MPI_COMM = 0 !< MPI_Communicator
            INTEGER :: NCPU     = 0 !< total number of proc in this communicator
            INTEGER :: NODE_ME  = 0 !< node id starting from 1 ... NCPU
            INTEGER :: COLOR    = 0 !< color of process within (sub)communicator
        END TYPE ML_MPI_PAR


        INTERFACE M_BCAST
           MODULE PROCEDURE M_BCAST_I_SCALAR
           MODULE PROCEDURE M_BCAST_R_SCALAR
           MODULE PROCEDURE M_BCAST_L_SCALAR
           MODULE PROCEDURE M_BCAST_I_1D
           MODULE PROCEDURE M_BCAST_R_1D
           MODULE PROCEDURE M_BCAST_L_1D
           MODULE PROCEDURE M_BCAST_I_2D
           MODULE PROCEDURE M_BCAST_R_2D
           MODULE PROCEDURE M_BCAST_L_2D
           MODULE PROCEDURE M_BCAST_I_3D
           MODULE PROCEDURE M_BCAST_R_3D
           MODULE PROCEDURE M_BCAST_L_3D
        END INTERFACE M_BCAST

        CONTAINS

!******************************************************************************************
! MPI Broadcast routines
!******************************************************************************************
        SUBROUTINE M_BCAST_I_SCALAR(COMM,SCAL,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: SCAL
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(SCAL,1,MPI_INTEGER,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(SCAL,1,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_R_SCALAR(COMM,SCAL,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: SCAL
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(SCAL,1,MPI_REAL8,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(SCAL,1,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_L_SCALAR(COMM,SCAL,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: SCAL
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(SCAL,1,MPI_LOGICAL,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(SCAL,1,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_I_1D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: VEC(:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_R_1D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: VEC(:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_L_1D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: VEC(:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_I_2D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: VEC(:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1,1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1,1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_R_2D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: VEC(:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1,1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1,1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_L_2D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: VEC(:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1,1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1,1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_I_3D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           INTEGER           :: VEC(:,:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1,1,1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1,1,1),N,MPI_INTEGER,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_R_3D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           REAL(q)           :: VEC(:,:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1,1,1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1,1,1),N,MPI_REAL8,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

        SUBROUTINE M_BCAST_L_3D(COMM,VEC,N,ROOT)
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM
           LOGICAL           :: VEC(:,:,:)
           INTEGER           :: N
           INTEGER           :: ROOT
           !Local variables
           INTEGER           :: IERR
           INTEGER           :: REQUEST
           INTEGER           :: STATUS(MPI_STATUS_SIZE)
#ifdef MPI_avoid_bcast
           CALL MPI_IBCAST(VEC(1,1,1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,REQUEST,IERR)
           CALL MPI_WAIT(REQUEST,STATUS,IERR)
#else
           CALL MPI_BCAST(VEC(1,1,1),N,MPI_LOGICAL,ROOT,COMM%MPI_COMM,IERR)
#ifdef MPI_bcast_with_barrier
           CALL MPI_BARRIER(COMM%MPI_COMM,IERR)
#endif
#endif
        END SUBROUTINE

      END MODULE MPI_HELP
#endif
