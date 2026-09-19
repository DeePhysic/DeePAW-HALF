!#define debug
#include "symbol.inc"
!************************************************************************
! RCS:  $Id: wave.F, v 1.6 2002/08/14 13:59:43 kresse Exp $
!>  this module contains the routines required to setup
!>  the distribution of wavefunctions over nodes and all basic routines
!>  handling wavedes etc.
!
!***********************************************************************
MODULE WAVE
  USE prec
  USE mpimy
  USE mgrid
  USE wave_struct_def

  !> the flag determines whether the serial of parallel FFT is used for
  !> wavefunctions, if one node holds all data for an entire wavefunction
  !> routines like HF require LUSE_PARALLEL_FFT=.FALSE., hence
  !> I recommend to use that option
  !> It is usually faster anyway
  LOGICAL:: LUSE_PARALLEL_FFT=.FALSE.

  !> used in the isEmpty and isFilled to define empty and filled orbitals
  REAL(q), PARAMETER:: tiny_occupation = 1E-4_q

CONTAINS

!=======================================================================
!
!>  allocate the descriptor for the orbitals
!>
!>  if LEXTEND is set the kinetic energy array and
!>  the index arrays IGX, Y and Z are allocated as well
!
!=======================================================================

  SUBROUTINE ALLOCWDES(WDES,LEXTEND)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavedes)  WDES
    INTEGER NRPLWV,NKPTS,NCOL,NKDIM
    LOGICAL LEXTEND

    NRPLWV=WDES%NGDIM
    NKPTS =WDES%NKPTS
    NKDIM =WDES%NKDIM   ! total number of k-points in the 1st BZ (for HF)
    NCOL  =WDES%NCOL
    ! a few arrays are allocated with the size NKDIM
    ALLOCATE( WDES%NPLWKP(NKDIM),WDES%NGVECTOR(NKDIM),WDES%NGVECTOR_POS(NKDIM),WDES%NPLWKP_TOT(NKDIM), &
         WDES%NINDPW(NRPLWV,NKDIM),WDES%NB_TOTK(NKDIM,2),WDES%LUSEINV(NKDIM),WDES%AT_GAMMA(NKDIM))
    NULLIFY(WDES%NINDPW_INV, WDES%FFTSCA, WDES%MAP_TO_FULL)
    WDES%NPLWKP    =0
    WDES%NPLWKP_TOT=0
    WDES%NGVECTOR  =0
    WDES%NGVECTOR_POS=0
    WDES%LUSEINV=.FALSE.
    IF (NCOL>0) THEN
       ALLOCATE(WDES%PL_INDEX(NCOL,NKPTS),WDES%PL_COL(NCOL,NKPTS))
    ELSE
       NULLIFY(WDES%PL_INDEX); NULLIFY(WDES%PL_COL)
    ENDIF

    IF (LEXTEND) THEN
       ALLOCATE(WDES%DATAKE(NRPLWV,2,NKPTS), &
            &    WDES%IGX(NRPLWV,NKDIM),WDES%IGY(NRPLWV,NKDIM),WDES%IGZ(NRPLWV,NKDIM))
DOESI  WDES%DATAKE = 0.0_q ! Breaks Co_2E4W_T (intel), uninitialized parts enter BUILD_RED_CHI_OR_EPSILON.
    ELSE
       NULLIFY(WDES%DATAKE)
       NULLIFY(WDES%IGX); NULLIFY(WDES%IGY); NULLIFY(WDES%IGZ)
       NULLIFY(WDES%PL_INDEX); NULLIFY(WDES%PL_COL)
    END IF

  END SUBROUTINE ALLOCWDES

!=======================================================================
!>  deallocate a  descriptor for the orbitals
!=======================================================================

  SUBROUTINE DEALLOCWDES(WDES,LEXTEND)
    IMPLICIT NONE
    TYPE (wavedes)  WDES
    LOGICAL LEXTEND

    DEALLOCATE( WDES%NPLWKP,WDES%NGVECTOR,WDES%NGVECTOR_POS,WDES%NPLWKP_TOT, &
         WDES%NINDPW, WDES%NB_TOTK, WDES%LUSEINV, WDES%AT_GAMMA)
    IF (ASSOCIATED(WDES%NINDPW_INV))  DEALLOCATE(WDES%NINDPW_INV)
    IF (ASSOCIATED(WDES%FFTSCA))      DEALLOCATE(WDES%FFTSCA)
    IF (ASSOCIATED(WDES%MAP_TO_FULL)) DEALLOCATE(WDES%MAP_TO_FULL)

    NULLIFY( WDES%NPLWKP,WDES%NGVECTOR,WDES%NGVECTOR_POS,WDES%NPLWKP_TOT, &
         WDES%NINDPW, WDES%NINDPW_INV, WDES%FFTSCA, WDES%MAP_TO_FULL, WDES%NB_TOTK, WDES%LUSEINV, WDES%AT_GAMMA)

    IF (WDES%NCOL>0) THEN
       DEALLOCATE(WDES%PL_INDEX,WDES%PL_COL)
       NULLIFY(WDES%PL_INDEX); NULLIFY(WDES%PL_COL)
    ENDIF
    IF (LEXTEND) THEN
       DEALLOCATE( WDES%DATAKE,WDES%IGX,WDES%IGY,WDES%IGZ)
       NULLIFY(WDES%DATAKE, WDES%IGX, WDES%IGY, WDES%IGZ)
    END IF
  END SUBROUTINE DEALLOCWDES


!=======================================================================
!>  free a descriptor and overwrite it by a new one
!=======================================================================

  SUBROUTINE COPYWDES(WDES_OLD, WDES_NEW ,LEXTEND)
    IMPLICIT NONE
    TYPE (wavedes)  WDES_OLD, WDES_NEW
    LOGICAL LEXTEND

    CALL DEALLOCWDES(WDES_OLD,LEXTEND)
    WDES_OLD=WDES_NEW

  END SUBROUTINE COPYWDES

!>
!> Create a new wavefunction descriptor from a previously existing one
!>
  SUBROUTINE WDES_FROM_WDES(SELF,WDES)
    IMPLICIT NONE
    TYPE(WAVEDES),INTENT(OUT) :: SELF
    TYPE(WAVEDES),TARGET,INTENT(IN) :: WDES
    CALL WDES_INIT(SELF, WDES%ISPIN, WDES%NSIM, WDES%QSPIRAL, WDES%SAXIS, WDES%LORBITALREAL, &
                   WDES%LNONCOLLINEAR, WDES%LSORBIT, WDES%LSPIRAL, WDES%LZEROZ, &
                   WDES%NBANDSLOW, WDES%NBANDSHIGH, WDES%RANDOM_NUMBER_GENERATOR, WDES%PCG_SEED, &
                   WDES%ENMAX, WDES%NB_PAR, WDES%NB_TOT)
    CALL WDES_SET_COMMS(SELF,WDES%COMM,WDES%COMM_KINTER,&
                        WDES%COMM_KIN,WDES%COMM_INTER,WDES%COMM_INB)
#ifdef use_shmem
    CALL WDES_SET_SHMEM_COMMS(SELF,WDES%COMM_SHMEM,WDES%COMM_intra_node,WDES%COMM_inter_node)
#endif
  END SUBROUTINE WDES_FROM_WDES

!>
!> Initialize a wavefunction descriptor
!>
  SUBROUTINE WDES_INIT(SELF, ISPIN, NSIM, QSPIRAL, SAXIS, LORBITALREAL, &
                       LNONCOLLINEAR, LSORBIT, LSPIRAL, LZEROZ, &
                       NBANDSLOW, NBANDSHIGH, RANDOM_NUMBER_GENERATOR, PCG_SEED, &
                       ENMAX, NB_PAR, NB_TOT)
    IMPLICIT NONE
    TYPE(wavedes) :: SELF
    INTEGER,INTENT(IN) :: ISPIN
    INTEGER,INTENT(IN) :: NSIM
    REAL(q),INTENT(IN) :: QSPIRAL(3)
    REAL(q),INTENT(IN) :: SAXIS(3)
    LOGICAL,INTENT(IN) :: LORBITALREAL
    LOGICAL,INTENT(IN) :: LNONCOLLINEAR
    LOGICAL,INTENT(IN) :: LSORBIT
    LOGICAL,INTENT(IN) :: LSPIRAL
    LOGICAL,INTENT(IN) :: LZEROZ
    INTEGER,INTENT(IN) :: NBANDSLOW
    INTEGER,INTENT(IN) :: NBANDSHIGH
    INTEGER,INTENT(IN) :: RANDOM_NUMBER_GENERATOR
    INTEGER(qi8),INTENT(IN) :: PCG_SEED
    REAL(q),INTENT(IN) :: ENMAX
    INTEGER,INTENT(IN) :: NB_PAR
    INTEGER,INTENT(IN) :: NB_TOT

    SELF%NSIM = NSIM
    SELF%QSPIRAL = QSPIRAL
    SELF%SAXIS = SAXIS
    SELF%LORBITALREAL = LORBITALREAL
    SELF%LNONCOLLINEAR = LNONCOLLINEAR
    SELF%LSORBIT = LSORBIT
    SELF%LSPIRAL = LSPIRAL
    SELF%LZEROZ = LZEROZ
    SELF%NBANDSLOW = NBANDSLOW
    SELF%NBANDSHIGH = NBANDSHIGH
    SELF%RANDOM_NUMBER_GENERATOR = RANDOM_NUMBER_GENERATOR
    SELF%PCG_SEED = PCG_SEED
    SELF%ENMAX = ENMAX
    SELF%NB_TOT = NB_TOT
    SELF%NB_PAR = NB_PAR
    SELF%ISPIN =ISPIN

    SELF%RSPIN=3-ISPIN
    IF (SELF%LNONCOLLINEAR) then
       SELF%NRSPINORS = 2
       SELF%RSPIN=1.0_q
    ELSE
       SELF%NRSPINORS = 1
    ENDIF

#ifdef gammareal
    SELF%LGAMMA=.TRUE.
#else
    SELF%LGAMMA=.FALSE.
#endif
    IF      (SELF%ISPIN==1  .AND. .NOT. SELF%LNONCOLLINEAR ) THEN
       SELF%NCDIJ=1
    ELSE IF (SELF%ISPIN==2) THEN
       SELF%NCDIJ=2
    ELSE IF (SELF%ISPIN==1  .AND. SELF%LNONCOLLINEAR ) THEN
       SELF%NCDIJ=4
    ELSE
       CALL vtutor%bug("internal error: can not set NCDIJ", __FILE__, __LINE__)
    ENDIF
  END SUBROUTINE WDES_INIT

!>
!> Set k-points in the wavefunction descriptor
!>
  SUBROUTINE WDES_SET_KPOINTS(SELF,VKPT,WTKPT)
    IMPLICIT NONE
    TYPE(WAVEDES) :: SELF
    REAL(q) :: VKPT(:,:)
    REAL(q) :: WTKPT(:)
    ! local variables
    INTEGER :: NKPTS
    IF (SIZE(VKPT,2)/=SIZE(WTKPT)) THEN
        call vtutor%bug('Size of kpoints '//str(SIZE(VKPT,2))//' and weights arrays '//&
                        str(SIZE(WTKPT))//' are incompatible',__FILE__,__LINE__)
    ENDIF
    NKPTS = SIZE(VKPT,2)
    SELF%NKDIM = NKPTS
    SELF%NKPTS = NKPTS
    SELF%NKPTS_FOR_GEN_LAYOUT = NKPTS
    ALLOCATE(SELF%VKPT(3,NKPTS))
    SELF%VKPT = VKPT
    ALLOCATE(SELF%WTKPT(NKPTS))
    SELF%WTKPT = WTKPT
  END SUBROUTINE WDES_SET_KPOINTS

  !> set the comunicators of the wavefunction descriptor
  SUBROUTINE WDES_SET_COMMS(SELF,COMM,COMM_KINTER,COMM_KIN,COMM_INTER,COMM_INB)
    IMPLICIT NONE
    TYPE(wavedes) :: self
    TYPE(communic),INTENT(IN) :: COMM
    TYPE(communic),INTENT(IN) :: COMM_KINTER
    TYPE(communic),INTENT(IN) :: COMM_KIN
    TYPE(communic),INTENT(IN) :: COMM_INTER
    TYPE(communic),INTENT(IN) :: COMM_INB
    ! local variables
    INTEGER :: NPAR

    ! allocate communicators
    ALLOCATE(SELF%COMM)
    ALLOCATE(SELF%COMM_KINTER)
    ALLOCATE(SELF%COMM_KIN)
    ALLOCATE(SELF%COMM_INTER)
    ALLOCATE(SELF%COMM_INB)

    ! set comunicators
    SELF%COMM = COMM
    SELF%COMM_KINTER = COMM_KINTER
    SELF%COMM_KIN = COMM_KIN
    SELF%COMM_INTER = COMM_INTER
    SELF%COMM_INB = COMM_INB

    ! set dimensions
    SELF%NB_LOW = COMM_INTER%NODE_ME
    NPAR = COMM_INTER%NCPU
    SELF%NB_PAR = NPAR
    SELF%NBANDS = SELF%NB_TOT/NPAR

    SELF%DO_REDIS = .FALSE.
    IF (SELF%NB_PAR /=1) THEN
       SELF%DO_REDIS=.TRUE.
    ENDIF
  END SUBROUTINE WDES_SET_COMMS

  SUBROUTINE WDES_SET_SHMEM_COMMS(self,COMM_SHMEM,COMM_intra_node,COMM_inter_node)
    TYPE(wavedes) :: SELF
    TYPE(communic),INTENT(IN) :: COMM_SHMEM
    TYPE(communic),INTENT(IN) :: COMM_INTRA_NODE
    TYPE(communic),INTENT(IN) :: COMM_INTER_NODE
    ! set shmem communicators
    ALLOCATE( SELF%COMM_SHMEM )
    ALLOCATE( SELF%COMM_intra_node )
    ALLOCATE( SELF%COMM_inter_node )
    SELF%COMM_SHMEM=COMM_SHMEM
    SELF%COMM_intra_node=COMM_intra_node
    SELF%COMM_inter_node=COMM_inter_node
  END SUBROUTINE WDES_SET_SHMEM_COMMS

!=======================================================================
!> @brief copy a wavefunction descriptor to another one
!=======================================================================

  SUBROUTINE WDES_COPY(WDES_DST, WDES_SRC)
    IMPLICIT NONE
    TYPE (wavedes)  WDES_DST, WDES_SRC

    WDES_DST%RSPIN               = WDES_SRC%RSPIN
    WDES_DST%ENMAX               = WDES_SRC%ENMAX
    WDES_DST%NRSPINORS           = WDES_SRC%NRSPINORS
    WDES_DST%NGDIM               = WDES_SRC%NGDIM
    WDES_DST%NRPLWV              = WDES_SRC%NRPLWV
    WDES_DST%NRPLWV_RED          = WDES_SRC%NRPLWV_RED
    WDES_DST%NPROD               = WDES_SRC%NPROD
    WDES_DST%NPRO                = WDES_SRC%NPRO
    WDES_DST%NPRO_TOT            = WDES_SRC%NPRO_TOT
    WDES_DST%NPROD_RED           = WDES_SRC%NPROD_RED
    WDES_DST%NBANDS              = WDES_SRC%NBANDS
    WDES_DST%NB_TOT              = WDES_SRC%NB_TOT
    WDES_DST%NB_PAR              = WDES_SRC%NB_PAR
    WDES_DST%NSIM                = WDES_SRC%NSIM
    WDES_DST%NB_LOW              = WDES_SRC%NB_LOW
    WDES_DST%NKDIM               = WDES_SRC%NKDIM
    WDES_DST%NKPTS_FOR_GEN_LAYOUT= WDES_SRC%NKPTS_FOR_GEN_LAYOUT
    WDES_DST%NKPTS               = WDES_SRC%NKPTS
    WDES_DST%ISPIN               = WDES_SRC%ISPIN
    WDES_DST%NCDIJ               = WDES_SRC%NCDIJ
    WDES_DST%NIONS               = WDES_SRC%NIONS
    WDES_DST%NTYP                = WDES_SRC%NTYP

    WDES_DST%LORBITALREAL  = WDES_SRC%LORBITALREAL
    WDES_DST%LOVERL        = WDES_SRC%LOVERL
    WDES_DST%DO_REDIS      = WDES_SRC%DO_REDIS
    WDES_DST%LNONCOLLINEAR = WDES_SRC%LNONCOLLINEAR
    WDES_DST%LSORBIT       = WDES_SRC%LSORBIT
    WDES_DST%LGAMMA        = WDES_SRC%LGAMMA
    WDES_DST%LSPIRAL       = WDES_SRC%LSPIRAL
    WDES_DST%LZEROZ        = WDES_SRC%LZEROZ
    WDES_DST%NBANDSLOW     = WDES_SRC%NBANDSLOW
    WDES_DST%NBANDSHIGH    = WDES_SRC%NBANDSHIGH
    WDES_DST%NCOL          = WDES_SRC%NCOL
    WDES_DST%QSPIRAL       = WDES_SRC%QSPIRAL
    WDES_DST%RANDOM_NUMBER_GENERATOR = WDES_SRC%RANDOM_NUMBER_GENERATOR
    WDES_DST%PCG_SEED                = WDES_SRC%PCG_SEED

#define ALLOC_COPY(a,b) ALLOCATE(a,source=b); a=b
    ALLOC_COPY(WDES_DST%NPLWKP,      WDES_SRC%NPLWKP)
    ALLOC_COPY(WDES_DST%NGVECTOR,    WDES_SRC%NGVECTOR)
    ALLOC_COPY(WDES_DST%NGVECTOR_POS,WDES_SRC%NGVECTOR_POS)
    ALLOC_COPY(WDES_DST%NPLWKP_TOT,  WDES_SRC%NPLWKP_TOT)
    ALLOC_COPY(WDES_DST%NB_TOTK,     WDES_SRC%NB_TOTK)
    ALLOC_COPY(WDES_DST%NPRO_POS,    WDES_SRC%NPRO_POS)
    ALLOC_COPY(WDES_DST%LMMAX,       WDES_SRC%LMMAX)
    ALLOC_COPY(WDES_DST%LMBASE,      WDES_SRC%LMBASE)
    ALLOC_COPY(WDES_DST%NITYP,       WDES_SRC%NITYP)
    ALLOC_COPY(WDES_DST%ITYP,        WDES_SRC%ITYP)
    ALLOC_COPY(WDES_DST%NT_GLOBAL,   WDES_SRC%NT_GLOBAL)
    ALLOC_COPY(WDES_DST%VKPT,        WDES_SRC%VKPT)
    ALLOC_COPY(WDES_DST%WTKPT,       WDES_SRC%WTKPT)
    ALLOC_COPY(WDES_DST%NINDPW,      WDES_SRC%NINDPW)
    ALLOC_COPY(WDES_DST%LUSEINV,     WDES_SRC%LUSEINV)
    WDES_DST%NINDPW_INV => NULL()
    WDES_DST%FFTSCA => NULL()
    WDES_DST%MAP_TO_FULL => NULL()
    ALLOC_COPY(WDES_DST%AT_GAMMA,    WDES_SRC%AT_GAMMA)
    ALLOC_COPY(WDES_DST%IGX,         WDES_SRC%IGX)
    ALLOC_COPY(WDES_DST%IGY,         WDES_SRC%IGY)
    ALLOC_COPY(WDES_DST%IGZ,         WDES_SRC%IGZ)
    ALLOC_COPY(WDES_DST%DATAKE,      WDES_SRC%DATAKE)
    ALLOC_COPY(WDES_DST%PL_INDEX,    WDES_SRC%PL_INDEX)
    ALLOC_COPY(WDES_DST%PL_COL,      WDES_SRC%PL_COL)

    ALLOCATE(WDES_DST%COMM)
    ALLOCATE(WDES_DST%COMM_INTER)
    ALLOCATE(WDES_DST%COMM_KINTER)
    ALLOCATE(WDES_DST%COMM_KIN)
    ALLOCATE(WDES_DST%COMM_INB)

    WDES_DST%COMM        =WDES_SRC%COMM
    WDES_DST%COMM_INTER  =WDES_SRC%COMM_INTER
    WDES_DST%COMM_KINTER =WDES_SRC%COMM_KINTER
    WDES_DST%COMM_KIN    =WDES_SRC%COMM_KIN
    WDES_DST%COMM_INB    =WDES_SRC%COMM_INB

  END SUBROUTINE WDES_COPY

!=======================================================================
!  initialize k-point related quantities
!=======================================================================

  SUBROUTINE INIT_KPOINT_WDES(WDES, KPOINTS )
    USE mkpoints_struct_def, only: kpoints_struct
    IMPLICIT NONE
    TYPE (wavedes)  WDES
    TYPE (kpoints_struct) KPOINTS
    INTEGER :: NK

    WDES%NKDIM =KPOINTS%NKPTS
    WDES%NKPTS =KPOINTS%NKPTS
    WDES%NKPTS_FOR_GEN_LAYOUT =KPOINTS%NKPTS
    WDES%VKPT  =>KPOINTS%VKPT
    WDES%WTKPT =>KPOINTS%WTKPT

  END SUBROUTINE INIT_KPOINT_WDES

!=======================================================================
!>  initialize the projector part of the descriptor for the
!>  orbitals
!=======================================================================

  SUBROUTINE WDES_SET_NPRO(WDES,T_INFO,P,LOVERL)
    USE  poscar
    USE  pseudo

    TYPE (wavedes)  WDES
    TYPE (type_info) :: T_INFO
    TYPE (potcar)   P(T_INFO%NTYP)
    LOGICAL :: LOVERL
    ! local varibles
    INTEGER NALLOC,NPRO_TOT,NT,NI,NIS,NODE_TARGET,NPRO, &
         LMMAXC,NIONS,LASTTYP
#ifndef MPI
    WDES%NIONS = T_INFO%NIONS
    WDES%NTYP  = T_INFO%NTYP
    WDES%NITYP =>T_INFO%NITYP
    WDES%ITYP  =>T_INFO%ITYP

    ALLOCATE(WDES%LMMAX(WDES%NTYP),WDES%NT_GLOBAL(WDES%NTYP))
    DO NT=1,T_INFO%NTYP
       WDES%LMMAX(NT)=P(NT)%LMMAX
       WDES%NT_GLOBAL(NT)=NT
    ENDDO
    WDES%NPRO  =SUM(WDES%LMMAX*WDES%NITYP)
    WDES%NPRO_TOT=SUM(WDES%LMMAX*WDES%NITYP)
    WDES%NPROD =WDES%NPRO
    NULLIFY(WDES%NPRO_POS)
#else
    !-----------------------------------------------------------------------
    ! parallel version
    !-----------------------------------------------------------------------
    TYPE (communic) COMM
    ! first count number of projection operators
    COMM = WDES%COMM_INB
    NPRO_TOT=0
    DO NT=1,T_INFO%NTYP
       NPRO_TOT=NPRO_TOT+P(NT)%LMMAX*T_INFO%NITYP(NT)
    ENDDO
    ! check implementation of NI_LOCAL and NI_GLOBAL
    DO NI=1,T_INFO%NIONS
       NI_L=NI_LOCAL(NI,COMM)
       IF (NI_L /= 0) THEN
          IF (NI /= NI_GLOBAL(NI_L,COMM)) THEN
             CALL vtutor%bug("internal error NI_GLOBAL, NI_LOCAL", __FILE__, __LINE__)
          ENDIF
       ENDIF
    ENDDO
    ! number of projection operators per node
    NALLOC=(T_INFO%NIONS+COMM%NCPU-1)/COMM%NCPU

    WDES%NIONS=0
    WDES%NTYP =0
    WDES%NPRO =0
    WDES%NPRO_TOT=NPRO_TOT
    ALLOCATE(WDES%NITYP(T_INFO%NTYP))
    ALLOCATE(WDES%ITYP(T_INFO%NIONS))
    ALLOCATE(WDES%LMMAX(T_INFO%NTYP),WDES%NT_GLOBAL(T_INFO%NTYP))
    ALLOCATE(WDES%NPRO_POS(NALLOC))
    WDES%LMMAX=0
    WDES%NITYP=0
    WDES%ITYP=0
    WDES%NT_GLOBAL=0

    LASTTYP=0
    NIS    =1

    NPRO=0
    DO NT=1,T_INFO%NTYP
       LMMAXC=P(NT)%LMMAX
!       IF (LMMAXC==0) GOTO 100
       DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
          ! does this element reside on local node
          IF (NI_LOCAL(NI,COMM) /=0 ) THEN
             WDES%NIONS=WDES%NIONS+1
             WDES%NPRO_POS(WDES%NIONS)=NPRO
             WDES%NPRO =WDES%NPRO +LMMAXC
             IF (NT /= LASTTYP) THEN
                WDES%NTYP=WDES%NTYP+1
                LASTTYP=NT
             ENDIF
             WDES%LMMAX(WDES%NTYP)=LMMAXC
             WDES%NT_GLOBAL(WDES%NTYP)=NT
             WDES%NITYP(WDES%NTYP)=WDES%NITYP(WDES%NTYP)+1
             WDES%ITYP(WDES%NIONS)=WDES%NTYP
          ENDIF
          NPRO= LMMAXC+NPRO
       ENDDO

100    NIS = NIS+T_INFO%NITYP(NT)
    ENDDO
    ! check whether everything is right
    NPRO=SUM(WDES%LMMAX*WDES%NITYP)
    IF (NPRO/= WDES%NPRO) THEN
       CALL vtutor%bug("internal ERROR 1 in WDES_SET_NPRO: " // str(NPRO) // " " // str(WDES%NPRO), __FILE__, __LINE__)
    ENDIF
    CALLMPI( M_sum_i(COMM,NPRO ,1))
    IF (NPRO/= NPRO_TOT) THEN
       CALL vtutor%bug("internal ERROR 2 in WDES_SET_NPRO: " // str(NPRO) // " " // str(NPRO_TOT), __FILE__, __LINE__)
    ENDIF
    ! make NPROD dividable by number of NB_PAR
    WDES%NPROD=((WDES%NPRO+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    ! and set it to the maximum value of all processors
    CALLMPI (M_max_i(COMM, WDES%NPROD, 1)) ! (required for  SHMALLOC)
#endif

    WDES%NPRO  =WDES%NPRO      *WDES%NRSPINORS
    WDES%NPRO_TOT=WDES%NPRO_TOT*WDES%NRSPINORS
    WDES%NPROD =WDES%NPROD*WDES%NRSPINORS
    WDES%NPROD_RED =WDES%NPROD /WDES%NB_PAR
#ifdef MPI
    IF (WDES%NB_PAR /= WDES%COMM_INTER%NCPU) THEN
       CALL vtutor%bug("internal ERROR 3 in WDES_SET_NPRO: " // str(WDES%NB_PAR) // " " // &
          str(WDES%COMM_INTER%NCPU), __FILE__, __LINE__)
    ENDIF
#endif

    WDES%LOVERL=LOVERL

    ! determine the position the projections onto local ion NI
    ! in the CPROJ arrays: (given by an offset) the first projection
    ! onto the projection operators residing on ion NI has index
    ! LMBASE(NI)+1
    ALLOCATE(WDES%LMBASE(WDES%NIONS+1))
    WDES%LMBASE(1)=0
    DO NI=1,WDES%NIONS
       NT=WDES%ITYP(NI)
       WDES%LMBASE(NI+1)=WDES%LMBASE(NI)+WDES%LMMAX(NT)
    ENDDO

  END SUBROUTINE WDES_SET_NPRO


!=======================================================================
!
!> this routine gives the local storage index for
!> the non local overlap CQIJ, strength CDIJ matrix elements
!> and for the projected orbitals
!> return is 0 if the element resides on an other processor
!> on entry NI is the global index
!
!=======================================================================

  FUNCTION NI_LOCAL(NI,COMM)
    IMPLICIT NONE
    TYPE (communic)  COMM
    INTEGER NI,NI_LOCAL
#ifdef MPI
    INTEGER NODE_TARGET

    NODE_TARGET=MOD(NI-1, COMM%NCPU)+1
    NI_LOCAL   =   (NI-1)/COMM%NCPU +1

    IF (NODE_TARGET /= COMM%NODE_ME) THEN
       NI_LOCAL=0
    ENDIF
#else
    !
    ! in conventional version all elements are local
    !
    NI_LOCAL=NI
#endif

    RETURN
  END FUNCTION NI_LOCAL

!=======================================================================
!
!> this routine gives the global ion index for
!> the non local overlap CQIJ, strength CDIJ matrix elements
!> and for the projected orbitals
!
!=======================================================================

  FUNCTION NI_GLOBAL(NI,COMM)
GPACC !$ACC ROUTINE SEQ
    IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(NI_GLOBAL)
    TYPE (communic)  COMM
    INTEGER NI, NI_GLOBAL
#ifdef MPI
    INTEGER NODE_TARGET
    NI_GLOBAL=(NI-1)*COMM%NCPU + COMM%NODE_ME
#else
    !
    ! in conventional version all elements are local
    !
    NI_GLOBAL=NI
#endif

    RETURN
  END FUNCTION NI_GLOBAL

!=======================================================================
!>  set WDES for one k-point
!>
!>  this is quite simple and sometimes necessary
!=======================================================================

  SUBROUTINE CREATE_SINGLE_KPOINT_WDES(WDES_ORIG,WDES,NK)
    IMPLICIT NONE
    TYPE (wavedes)  WDES,WDES_ORIG
    INTEGER NK

    WDES=WDES_ORIG
    WDES%NKPTS=1

    WDES%NPLWKP=> WDES_ORIG%NPLWKP(NK:NK)
    WDES%NGVECTOR=> WDES_ORIG%NGVECTOR(NK:NK)
    WDES%NGVECTOR_POS=> WDES_ORIG%NGVECTOR_POS(NK:NK)
    WDES%NPLWKP_TOT=> WDES_ORIG%NPLWKP_TOT(NK:NK)
    WDES%WTKPT => WDES_ORIG%WTKPT (NK:NK)
    WDES%VKPT  => WDES_ORIG%VKPT  (:,NK:NK)
    WDES%NINDPW=> WDES_ORIG%NINDPW(:,NK:NK)
    IF (ASSOCIATED(WDES_ORIG%NINDPW_INV)) THEN
       WDES%NINDPW_INV => WDES_ORIG%NINDPW_INV(:,NK:NK)
       WDES%FFTSCA     => WDES_ORIG%FFTSCA(:,:,NK:NK)
       WDES%MAP_TO_FULL=> WDES_ORIG%MAP_TO_FULL(:,:,NK:NK)
    ELSE
       NULLIFY(WDES%NINDPW_INV, WDES%FFTSCA, WDES%MAP_TO_FULL)
    ENDIF
    WDES%LUSEINV=WDES_ORIG%LUSEINV(NK)
    WDES%IGX   => WDES_ORIG%IGX   (:,NK:NK)
    WDES%IGY   => WDES_ORIG%IGY   (:,NK:NK)
    WDES%IGZ   => WDES_ORIG%IGZ   (:,NK:NK)
    WDES%DATAKE=> WDES_ORIG%DATAKE(:,:,NK:NK)
    WDES%AT_GAMMA=WDES_ORIG%AT_GAMMA(NK)
    WDES%LORBITALREAL=WDES_ORIG%LORBITALREAL

    IF (WDES%NCOL>0) THEN
       WDES%PL_INDEX => WDES_ORIG%PL_INDEX(:,NK:NK)
       WDES%PL_COL   => WDES_ORIG%PL_COL (:,NK:NK)
    ENDIF

  END SUBROUTINE CREATE_SINGLE_KPOINT_WDES

!=======================================================================
!>  initialize the storage for the orbitals
!=======================================================================

  SUBROUTINE ALLOCW(WDES,W,WUP,WDW,ALLOC_PINNED)
#ifdef CUDA_PINNED_MEMORY
    USE mcuda_pinned_memory
#endif
    USE ini
    USE tutor, ONLY: argument, VTUTOR, ISERROR, ALLOC_WAVEFUN
    IMPLICIT NONE
    INTEGER NK
    TYPE (wavedes), TARGET ::  WDES
    TYPE (wavespin) W
    TYPE (wavefun), OPTIONAL :: WUP,WDW

    LOGICAL, OPTIONAL :: ALLOC_PINNED

    INTEGER NRPLWV,NPROD,NKPTS,NBANDS,ISPIN,NB_TOT,NB_PAR,NB_LOW

    INTEGER              :: ISTATUS      ! allocation error code
    TYPE(argument)       :: ARG

#ifdef CUDA_PINNED_MEMORY
    LOGICAL :: LPINNED
    LPINNED=.FALSE. ; IF (PRESENT(ALLOC_PINNED)) LPINNED=ALLOC_PINNED
#endif

    NRPLWV=WDES%NRPLWV
    NPROD =WDES%NPROD
    NKPTS =WDES%NKPTS
    NBANDS=WDES%NBANDS
    NB_TOT=WDES%NB_TOT
    NB_LOW=WDES%NB_LOW
    NB_PAR=WDES%NB_PAR
    ISPIN =WDES%ISPIN

    ALLOCATE(W%CPTWFP(NRPLWV,NBANDS,NKPTS,ISPIN), STAT=ISTATUS )
    IF (ISTATUS/=0) THEN
       ALLOCATE(ARG%IVAL(1))
       ARG%IVAL(1)=16._q*NRPLWV*NBANDS*NKPTS*ISPIN/1024/1024
       CALL VTUTOR%WRITE(ISERROR, ALLOC_WAVEFUN, ARG)
    ENDIF
#ifdef CUDA_PINNED_MEMORY
    IF (LPINNED) CALL PIN_MEMORY(C_LOC(W%CPTWFP),INT(SIZE(W%CPTWFP,KIND=C_SIZE_T)*SIZEOF(W%CPTWFP(1,1,1,1)),KIND=C_SIZE_T))
#endif

    ALLOCATE( W%CPROJ (MAX(NPROD,1), NBANDS,NKPTS,ISPIN), STAT=ISTATUS )
    IF (ISTATUS/=0) THEN
       ALLOCATE(ARG%IVAL(1))
       ARG%IVAL(1)= ndata * 8._q*NPROD*NBANDS*NKPTS*ISPIN/1024/1024
       CALL VTUTOR%WRITE(ISERROR, ALLOC_WAVEFUN, ARG)
    ENDIF
#ifdef CUDA_PINNED_MEMORY
    IF (LPINNED) CALL PIN_MEMORY(C_LOC(W%CPROJ),INT(SIZE(W%CPROJ,KIND=C_SIZE_T)*SIZEOF(W%CPROJ(1,1,1,1)),KIND=C_SIZE_T))
#endif

    IF (WDES%LGAMMA) THEN
       CALL REGISTER_ALLOCATE(16._q*SIZE(W%CPTWFP,KIND=qi8)+8._q *SIZE(W%CPROJ,KIND=qi8), "wavefun")
    ELSE
       CALL REGISTER_ALLOCATE(16._q*SIZE(W%CPTWFP,KIND=qi8)+16._q*SIZE(W%CPROJ,KIND=qi8), "wavefun")
    ENDIF

    ALLOCATE( &
         W%CELTOT(NB_TOT,NKPTS,ISPIN), &
         W%FERTOT(NB_TOT,NKPTS,ISPIN), &
         W%AUXTOT(NB_TOT,NKPTS,ISPIN))
    ALLOCATE(W%EFERMI(ISPIN))

    W%CPTWFP=0
    W%CPROJ =0
    W%CELTOT=0
    W%FERTOT=0
    W%EFERMI=0
    W%AUXTOT=1
    W%FERWE => W%FERTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
    W%AUX   => W%AUXTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
    W%CELEN => W%CELTOT(NB_LOW:NB_TOT:NB_PAR,:,:)

    W%OVER_BAND=.FALSE.
    W%WDES  => WDES

    IF (PRESENT(WUP)) THEN
       WUP%CELTOT=> W%CELTOT(:,:,1)
       WUP%FERTOT=> W%FERTOT(:,:,1)
       WUP%AUXTOT=> W%AUXTOT(:,:,1)
       WUP%CELEN => W%CELEN(:,:,1)
       WUP%FERWE => W%FERWE(:,:,1)
       WUP%AUX   => W%AUX  (:,:,1)
       WUP%CPTWFP=> W%CPTWFP(:,:,:,1)
       WUP%CPROJ => W%CPROJ(:,:,:,1)
       WUP%OVER_BAND=.FALSE.
       WUP%WDES  => WDES
    ENDIF

    IF (PRESENT(WDW) .AND. WDES%ISPIN==2) THEN
       WDW%CELTOT=> W%CELTOT(:,:,2)
       WDW%FERTOT=> W%FERTOT(:,:,2)
       WDW%AUXTOT=> W%AUXTOT(:,:,2)
       WDW%CELEN => W%CELEN(:,:,2)
       WDW%FERWE => W%FERWE(:,:,2)
       WDW%AUX   => W%AUX  (:,:,2)
       WDW%CPTWFP=> W%CPTWFP(:,:,:,2)
       WDW%CPROJ => W%CPROJ(:,:,:,2)
       WDW%OVER_BAND=.FALSE.
       WDW%WDES  => WDES
    ENDIF

  END SUBROUTINE ALLOCW


!=======================================================================
!>  initialize the storage for the orbitals without
!>  orbital coefficients
!=======================================================================

  SUBROUTINE ALLOCW_NOPLANEWAVE(WDES,W)
    USE ini
    IMPLICIT NONE
    INTEGER NK
    TYPE (wavedes), TARGET ::  WDES
    TYPE (wavespin) W

    INTEGER NRPLWV,NPROD,NKPTS,NBANDS,ISPIN,NB_TOT,NB_PAR,NB_LOW

    NRPLWV=WDES%NRPLWV
    NPROD =WDES%NPROD
    NKPTS =WDES%NKPTS
    NBANDS=WDES%NBANDS
    NB_TOT=WDES%NB_TOT
    NB_LOW=WDES%NB_LOW
    NB_PAR=WDES%NB_PAR
    ISPIN =WDES%ISPIN

    ALLOCATE( &
         W%CELTOT(NB_TOT,NKPTS,ISPIN), &
         W%FERTOT(NB_TOT,NKPTS,ISPIN), &
         W%AUXTOT(NB_TOT,NKPTS,ISPIN))
    ALLOCATE( W%EFERMI( ISPIN ) )

    W%EFERMI=0
    W%CELTOT=0
    W%FERTOT=0
    W%AUXTOT=1
    W%FERWE => W%FERTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
    W%AUX   => W%AUXTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
    W%CELEN => W%CELTOT(NB_LOW:NB_TOT:NB_PAR,:,:)

    W%OVER_BAND=.FALSE.
    W%WDES  => WDES

 END SUBROUTINE ALLOCW_NOPLANEWAVE


!=======================================================================
!>  initialize the storage for the orbitals without
!>  orbital coefficients
!=======================================================================

 SUBROUTINE DEALLOCW_NOPLANEWAVE(W)
    USE ini
    IMPLICIT NONE
    INTEGER NK
    TYPE (wavedes), TARGET ::  WDES
    TYPE (wavespin) W

    DEALLOCATE( &
         W%EFERMI, &
         W%CELTOT, &
         W%FERTOT, &
         W%AUXTOT)

  END SUBROUTINE DEALLOCW_NOPLANEWAVE

!=======================================================================
!>  dealloc of orbital array
!=======================================================================

  SUBROUTINE DEALLOCW(W)
#ifdef CUDA_PINNED_MEMORY
    USE mcuda_pinned_memory
#endif
    USE ini
    IMPLICIT NONE
    TYPE (wavespin) W

    IF (ASSOCIATED(W%CPTWFP)) THEN
       IF (W%WDES%LGAMMA) THEN
          CALL DEREGISTER_ALLOCATE(16._q*SIZE(W%CPTWFP,KIND=qi8)+8._q*SIZE(W%CPROJ,KIND=qi8) , "wavefun")
       ELSE
          CALL DEREGISTER_ALLOCATE(16._q*SIZE(W%CPTWFP,KIND=qi8)+16._q*SIZE(W%CPROJ,KIND=qi8), "wavefun")
       ENDIF
    ENDIF

#ifdef CUDA_PINNED_MEMORY
    CALL UNPIN_MEMORY(C_LOC(W%CPTWFP))
    CALL UNPIN_MEMORY(C_LOC(W%CPROJ ))
#endif

    IF (ASSOCIATED(W%CPTWFP) ) DEALLOCATE(W%CPTWFP)
    IF (ASSOCIATED(W%CPROJ ) ) DEALLOCATE(W%CPROJ )
    IF (ASSOCIATED(W%CELTOT) ) DEALLOCATE(W%CELTOT)
    IF (ASSOCIATED(W%FERTOT) ) DEALLOCATE(W%FERTOT)
    IF (ASSOCIATED(W%AUXTOT) ) DEALLOCATE(W%AUXTOT)
    IF (ASSOCIATED(W%EFERMI) ) DEALLOCATE(W%EFERMI)

    NULLIFY(W%EFERMI)
    NULLIFY(W%CPTWFP)
    NULLIFY(W%CPROJ)
    NULLIFY(W%CELTOT)
    NULLIFY(W%FERTOT)
    NULLIFY(W%AUXTOT)

  END SUBROUTINE DEALLOCW

!=======================================================================
!> dealloc of orbital array but CELTOT and FERTOT remain
!> allocated
!=======================================================================

  SUBROUTINE DEALLOCW_CW(W)

GPACC USE moffload

   USE ini
    IMPLICIT NONE
    TYPE (wavespin) W

    IF (ASSOCIATED(W%CPTWFP)) THEN
       IF (W%WDES%LGAMMA) THEN
          CALL DEREGISTER_ALLOCATE(16._q*SIZE(W%CPTWFP,KIND=qi8)+8._q*SIZE(W%CPROJ,KIND=qi8) , "wavefun")
       ELSE
          CALL DEREGISTER_ALLOCATE(16._q*SIZE(W%CPTWFP,KIND=qi8)+16._q*SIZE(W%CPROJ,KIND=qi8), "wavefun")
       ENDIF
    ENDIF

    !$ACC EXIT DATA DELETE(W%CPTWFP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W%CPTWFP,1)) FINALIZE ASYNC(ACC_ASYNC_Q)
    IF (ASSOCIATED(W%CPTWFP) ) DEALLOCATE(W%CPTWFP)
    !$ACC EXIT DATA DELETE(W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W%CPROJ,1)) FINALIZE ASYNC(ACC_ASYNC_Q)
    IF (ASSOCIATED(W%CPROJ) )  DEALLOCATE(W%CPROJ)
    NULLIFY(W%CPTWFP)
    NULLIFY(W%CPROJ)

  END SUBROUTINE DEALLOCW_CW


!=======================================================================
!>  nullify orbital array
!=======================================================================

  SUBROUTINE WNULLIFY(W)
    IMPLICIT NONE
    TYPE (wavespin) W

    NULLIFY(W%WDES,W%CPTWFP,W%CPROJ,W%FERTOT,W%FERWE,W%CELTOT,W%CELEN,W%AUXTOT,W%AUX)

  END SUBROUTINE WNULLIFY


!=======================================================================
!>  initialize a descriptor for one orbital  (wavedes1) from
!>  a descriptor of an array of orbitals  (wavedes)
!>  (kpoint index must be supplied)
!=======================================================================

  SUBROUTINE SETWDES(WDES,WDES1,NK)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER SETWDES
#endif

__GPU USE moffload

    IMPLICIT NONE
    INTEGER NK
    TYPE (wavedes)  WDES
    TYPE (wavedes1) WDES1

#ifdef OFFLOADING
    LOGICAL :: LPRESENT
    ! if WDES1 exists on the device we will first
    ! delete it and restore it afterwards
    LPRESENT=IS_ON_DEVICE(WDES1,1).AND.OFFLOAD_ON
    IF (LPRESENT) THEN
       CALL DELETE_TYPED_VAR(WDES1)
    ENDIF
#endif // OFFLOADING

    PROFILING_START('setwdes')

    WDES1%LNONCOLLINEAR  = WDES%LNONCOLLINEAR
    WDES1%LSORBIT        = WDES%LSORBIT
    WDES1%LOVERL         = WDES%LOVERL
    WDES1%DO_REDIS       = WDES%DO_REDIS
    WDES1%LORBITALREAL   = WDES%LORBITALREAL
    WDES1%LSPIRAL        = WDES%LSPIRAL
    WDES1%LGAMMA         = WDES%LGAMMA
    WDES1%LZEROZ         = WDES%LZEROZ

    WDES1%NRSPINORS      = WDES%NRSPINORS
    WDES1%NRPLWV         = WDES%NRPLWV
    WDES1%NRPLWV_RED     = WDES%NRPLWV_RED
    WDES1%NGDIM          = WDES%NGDIM
    WDES1%NPROD          = WDES%NPROD
    WDES1%NPROD_RED      = WDES%NPROD_RED
    WDES1%NPRO           = WDES%NPRO
    WDES1%NPRO_TOT       = WDES%NPRO_TOT
    WDES1%RSPIN          = WDES%RSPIN
    WDES1%RINPL          = 1._q/WDES%GRID%NPLWV

    WDES1%QSPIRAL        = WDES%QSPIRAL

    WDES1%NBANDS         = WDES%NBANDS
    WDES1%NB_TOT         = WDES%NB_TOT
    WDES1%NB_PAR         = WDES%NB_PAR
    WDES1%NB_LOW         = WDES%NB_LOW

    WDES1%NSIM           = WDES%NSIM

    WDES1%NIONS          = WDES%NIONS
    WDES1%NTYP           = WDES%NTYP

    WDES1%NITYP          =>WDES%NITYP
    WDES1%ITYP           =>WDES%ITYP
    WDES1%LMMAX          =>WDES%LMMAX
    WDES1%LMBASE         =>WDES%LMBASE
    WDES1%NT_GLOBAL      =>WDES%NT_GLOBAL
    WDES1%NPRO_POS       =>WDES%NPRO_POS
    WDES1%GRID           =>WDES%GRID

    WDES1%NCOL           = WDES%NCOL

    WDES1%COMM           =>WDES%COMM
    WDES1%COMM_INTER     =>WDES%COMM_INTER
    WDES1%COMM_INB       =>WDES%COMM_INB
    WDES1%COMM_KINTER    =>WDES%COMM_KINTER
    WDES1%COMM_KIN       =>WDES%COMM_KIN
#ifdef use_shmem
    WDES1%COMM_SHMEM     =>WDES%COMM_SHMEM
    WDES1%COMM_intra_node=>WDES%COMM_intra_node
    WDES1%COMM_inter_node=>WDES%COMM_inter_node
#endif
    WDES1%NK             = NK

    IF (NK/=0) THEN

    IF (NK > SIZE(WDES%NGVECTOR) ) THEN
       CALL vtutor%bug("internal error in SETWDES: NK exceeds SIZE " // str(NK) // " " // &
          str(SIZE(WDES%NGVECTOR)), __FILE__, __LINE__)
    ENDIF
    WDES1%NPL            = WDES%NPLWKP(NK)
    WDES1%NGVECTOR       = WDES%NGVECTOR(NK)
    WDES1%NGVECTOR_POS   = WDES%NGVECTOR_POS(NK)
    WDES1%NPL_TOT        = WDES%NPLWKP_TOT(NK)
    WDES1%NB_TOTK        = WDES%NB_TOTK(NK,1:2)
    WDES1%AT_GAMMA       = WDES%AT_GAMMA(NK)

    WDES1%IGX            =>WDES%IGX(:,NK)
    WDES1%IGY            =>WDES%IGY(:,NK)
    WDES1%IGZ            =>WDES%IGZ(:,NK)
    WDES1%VKPT           =>WDES%VKPT(:,NK)

    WDES1%NINDPW         =>WDES%NINDPW(:,NK)

    WDES1%LUSEINV        = WDES%LUSEINV(NK)

    IF(NK<=WDES%NKPTS) THEN
       WDES1%WTKPT       = WDES%WTKPT(NK)
       WDES1%DATAKE      =>WDES%DATAKE(:,:,NK)
    ELSE
       WDES1%WTKPT       = 0
    ENDIF
    IF (WDES%NCOL>0.AND.NK<=WDES%NKPTS) THEN
       WDES1%PL_INDEX    =>WDES%PL_INDEX(:,NK)
       WDES1%PL_COL      =>WDES%PL_COL(:,NK)
    ENDIF
    IF (ASSOCIATED(WDES%NINDPW_INV)) THEN
       WDES1%NINDPW_INV  =>WDES%NINDPW_INV(:,NK)
       WDES1%FFTSCA      =>WDES%FFTSCA(:,:,NK)
       WDES1%MAP_TO_FULL =>WDES%MAP_TO_FULL(:,:,NK)
    ELSE
       NULLIFY(WDES1%NINDPW_INV, WDES1%FFTSCA, WDES1%MAP_TO_FULL)
    ENDIF

    ENDIF

    WDES1%NPL_RED        = WDES1%NPL     ! number of plane waves/node after data redistribution
    WDES1%NPRO_RED       = WDES1%NPRO    ! number of projected wavef. after data redistribution

    CALLMPI( SET_NPL_NPRO(WDES1, WDES1%NPL_RED, WDES1%NPRO_RED))
    IF (WDES1%LOVERL) THEN
       WDES1%NPRO_O_RED  = WDES1%NPRO_RED
    ELSE
       WDES1%NPRO_O_RED  = 0
    ENDIF

#ifdef OFFLOADING
    IF (LPRESENT) THEN

       CALL COPYIN_TYPED_VAR(WDES1)

#ifdef CRAYHIP
! pointer attachments for CRAYHIP, for intel done in the above deep copy copyin_wavedes1
!$OMP TARGET IF(OFFLOAD_ON)
       WDES1%NITYP          =>WDES%NITYP
       WDES1%ITYP           =>WDES%ITYP
       WDES1%LMMAX          =>WDES%LMMAX
       WDES1%LMBASE         =>WDES%LMBASE
       WDES1%NT_GLOBAL      =>WDES%NT_GLOBAL
       WDES1%NPRO_POS       =>WDES%NPRO_POS
       WDES1%GRID           =>WDES%GRID

       WDES1%COMM           =>WDES%COMM
       WDES1%COMM_INTER     =>WDES%COMM_INTER
       WDES1%COMM_INB       =>WDES%COMM_INB
       WDES1%COMM_KINTER    =>WDES%COMM_KINTER
       WDES1%COMM_KIN       =>WDES%COMM_KIN
#ifdef use_shmem
       WDES1%COMM_SHMEM     =>WDES%COMM_SHMEM
       WDES1%COMM_intra_node=>WDES%COMM_intra_node
       WDES1%COMM_inter_node=>WDES%COMM_inter_node
#endif
!$OMP END TARGET
       IF (NK/=0) THEN
!$OMP TARGET IF(OFFLOAD_ON)
          WDES1%IGX            =>WDES%IGX(:,NK)
          WDES1%IGY            =>WDES%IGY(:,NK)
          WDES1%IGZ            =>WDES%IGZ(:,NK)
          WDES1%VKPT           =>WDES%VKPT(:,NK)
          WDES1%NINDPW         =>WDES%NINDPW(:,NK)
!$OMP END TARGET
          IF(NK<=WDES%NKPTS) THEN
!$OMP TARGET IF(OFFLOAD_ON)
             WDES1%DATAKE      =>WDES%DATAKE(:,:,NK)
!$OMP END TARGET
          ENDIF
          IF (WDES%NCOL>0.AND.NK<=WDES%NKPTS) THEN
!$OMP TARGET IF(OFFLOAD_ON)
             WDES1%PL_INDEX    =>WDES%PL_INDEX(:,NK)
             WDES1%PL_COL      =>WDES%PL_COL(:,NK)
!$OMP END TARGET
          ENDIF
          IF (ASSOCIATED(WDES%NINDPW_INV)) THEN
!$OMP TARGET IF(OFFLOAD_ON)
             WDES1%NINDPW_INV  =>WDES%NINDPW_INV(:,NK)
             WDES1%FFTSCA      =>WDES%FFTSCA(:,:,NK)
             WDES1%MAP_TO_FULL =>WDES%MAP_TO_FULL(:,:,NK)
!$OMP END TARGET
          ENDIF
       ENDIF
#endif // CRAYHIP

#ifdef ACC_OFFLOAD
!$ACC UPDATE DEVICE(WDES1%NPL,&
!$ACC&              WDES1%NGVECTOR, &
!$ACC&              WDES1%NGVECTOR_POS, &
!$ACC&              WDES1%NPL_TOT, &
!$ACC&              WDES1%NPL_RED, &
!$ACC&              WDES1%NPRO_RED, &
!$ACC&              WDES1%NPRO_O_RED, &
!$ACC&              WDES1%NB_TOTK, &
!$ACC&              WDES1%AT_GAMMA, &
!$ACC&              WDES1%LUSEINV, &
!$ACC&              WDES1%WTKPT, &
!$ACC&              WDES1%NINDPW, &
!$ACC&              WDES1%RINPL, &
!$ACC&              WDES1%NK) ASYNC(ACC_ASYNC_Q)
#endif // ACC_OFFLOAD

#ifdef OMP_OFFLOAD
!$OMP TARGET UPDATE TO(WDES1%NPL,&
!$OMP                  WDES1%NGVECTOR,&
!$OMP                  WDES1%NGVECTOR_POS,&
!$OMP                  WDES1%NPL_TOT,&
!$OMP                  WDES1%NPL_RED,&
!$OMP                  WDES1%NPRO_RED,&
!$OMP                  WDES1%NPRO_O_RED,&
!$OMP                  WDES1%NB_TOTK,&
!$OMP                  WDES1%AT_GAMMA,&
!$OMP                  WDES1%LUSEINV,&
!$OMP                  WDES1%WTKPT,&
!$OMP                  WDES1%NINDPW,&
!$OMP                  WDES1%RINPL,&
!$OMP                  WDES1%NK)
#endif // OMP_OFFLOAD
    ENDIF
#endif // OFFLOADING

    PROFILING_STOP('setwdes')

  END SUBROUTINE SETWDES

!=======================================================================
!>  initialize the optional data required for real space calculations
!>  in a descriptor for one single orbital (wavedes1)
!>  scheduled for removal
!=======================================================================

  SUBROUTINE SETWGRID_OLD(WDES1,GRID)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER SETWGRID_OLD
#endif
__GPU USE moffload

    IMPLICIT NONE

    TYPE (grid_3d), TARGET ::  GRID
    TYPE (wavedes1) WDES1
#ifdef OFFLOADING
    LOGICAL :: LPRESENT
    LPRESENT=IS_ON_DEVICE(WDES1,1).AND.OFFLOAD_ON
    IF (LPRESENT) THEN
       CALL DELETE_TYPED_VAR(WDES1%GRID)
    ENDIF
#endif
    WDES1%RINPL =1._q/GRID%NPLWV
    WDES1%GRID =>GRID
#ifdef OFFLOADING
    IF (LPRESENT) THEN
GPACC !$ACC UPDATE DEVICE(WDES1%RINPL) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(WDES1%RINPL) IF(OFFLOAD_ON)
       CALL COPYIN_TYPED_VAR(WDES1%GRID)
#ifdef CRAYWAR3
!$OMP TARGET IF(OFFLOAD_ON)
       WDES1%GRID => GRID
!$OMP END TARGET
#endif
    ENDIF
#endif
  END SUBROUTINE SETWGRID_OLD

!=======================================================================
!>  create storage for one orbital W
!>  optionally real space array is allocated
!=======================================================================

  SUBROUTINE NEWWAV(W1, WDES1, ALLOC_REAL, ISTAT)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER NEWWAV
#endif

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    LOGICAL ALLOC_REAL
    INTEGER, OPTIONAL :: ISTAT
    ! local variables
    INTEGER MPLWV
    IF (PRESENT(ISTAT)) THEN
       IF (ALLOC_REAL) THEN
          MPLWV=WDES1%GRID%MPLWV*WDES1%NRSPINORS
          ALLOCATE(W1%CPTWFP(WDES1%NRPLWV),W1%CPROJ(WDES1%NPROD),W1%CR(MPLWV),STAT=ISTAT)
DOESI     W1%CPTWFP=(0.0_q, 0.0_q); W1%CPROJ=0; W1%CR=(0.0_q, 0.0_q) ! No failed tests reported, see below.
GPACC !$ACC ENTER DATA CREATE(W1%CPTWFP,W1%CPROJ,W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: W1%CPTWFP,W1%CPROJ,W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: W1%CPTWFP,W1%CPROJ,W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
       ELSE
          ALLOCATE(W1%CPTWFP(WDES1%NRPLWV),W1%CPROJ(WDES1%NPROD),STAT=ISTAT)
DOESI     W1%CPTWFP=(0.0_q, 0.0_q); W1%CPROJ=0 ! No failed tests reported, see below.
          NULLIFY(W1%CR)
GPACC !$ACC ENTER DATA CREATE(W1%CPTWFP,W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: W1%CPTWFP,W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%CPTWFP,W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
       ENDIF
    ELSE
       IF (ALLOC_REAL) THEN
          MPLWV=WDES1%GRID%MPLWV*WDES1%NRSPINORS
          ALLOCATE(W1%CPTWFP(WDES1%NRPLWV),W1%CPROJ(WDES1%NPROD),W1%CR(MPLWV))
DOESI     W1%CPTWFP=(0.0_q, 0.0_q); W1%CPROJ=0; W1%CR=(0.0_q, 0.0_q) ! No failed tests reported, see below.
GPACC !$ACC ENTER DATA CREATE(W1%CPTWFP,W1%CPROJ,W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: W1%CPTWFP,W1%CPROJ,W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%CPTWFP,W1%CPROJ,W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
       ELSE
          ALLOCATE(W1%CPTWFP(WDES1%NRPLWV),W1%CPROJ(WDES1%NPROD))
DOESI     W1%CPTWFP=(0.0_q, 0.0_q); W1%CPROJ=0 ! Breaks bulk_SiO2_elastic_properties_ibrion8 (gnu), may be scaled without initilization in W1_GAXPY,...
          NULLIFY(W1%CR)
GPACC !$ACC ENTER DATA CREATE(W1%CPTWFP,W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: W1%CPTWFP,W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%CPTWFP,W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
       ENDIF
    ENDIF
    W1%FERWE=0
    W1%CELEN=0

    W1%WDES1 =>WDES1
    ! not a valid band of the Hamiltonian
    W1%NB=-1
    W1%ISP=-1
    W1%LDO=.TRUE.

GPACC !$ACC UPDATE DEVICE(W1%FERWE,W1%CELEN,W1%NB,W1%ISP,W1%LDO) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
!!      CALL COPYIN_TYPED_VAR(W1%WDES1)
!!      CALL ACC_ATTACH_ASYNC(W1%WDES1,ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE TO(W1%FERWE, W1%CELEN, W1%NB, W1%ISP, W1%LDO) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
GPOMP !$OMP TARGET ENTER DATA MAP(TO: W1, W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS,ALLOC: W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
#ifdef CRAYWAR3
GPOMP !$OMP TARGET MAP(PRESENT,TO: W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
    W1%WDES1 =>WDES1
GPOMP !$OMP END TARGET
#endif
  END SUBROUTINE NEWWAV

!=======================================================================
!>  create parallel storage for orbitals W
!>
!>  optionally real space array is allocated
!>  this version is for OPENMP creating OMP_GET_MAX_THREADS()
!>  W1 instances
!=======================================================================

  SUBROUTINE NEWWAV_OMP(W1, WDES1, ALLOC_REAL )
    IMPLICIT NONE
    TYPE (wavefun1), POINTER ::  W1(:)
    TYPE (wavedes1), TARGET :: WDES1
    LOGICAL ALLOC_REAL
!$  INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS
    INTEGER N

#ifdef _OPENMP
    ALLOCATE(W1(OMP_GET_MAX_THREADS()))
    DO N=1,OMP_GET_MAX_THREADS()
      CALL NEWWAV(W1(N), WDES1, ALLOC_REAL)
    ENDDO
#else
    ALLOCATE(W1(1))
    CALL NEWWAV(W1(1), WDES1, ALLOC_REAL)
#endif
  END SUBROUTINE NEWWAV_OMP

  SUBROUTINE DELWAV_OMP(W1, ALLOC_REAL )
    IMPLICIT NONE
    TYPE (wavefun1), POINTER ::  W1(:)
    LOGICAL ALLOC_REAL
    INTEGER N

    DO N=1,SIZE(W1,1)
      CALL DELWAV(W1(N), ALLOC_REAL )
    ENDDO
    DEALLOCATE(W1)
    NULLIFY(W1)

  END SUBROUTINE DELWAV_OMP

#ifdef use_shmem
!=======================================================================
!>  create storage for one orbital W
!>
!>  optionally real space array is allocated
!>  use shmem shared memory
!=======================================================================

  SUBROUTINE NEWWAV_SHMEM(W1,WDES1,ALLOC_REAL,address,shmid)
    USE c2f_interface
    IMPLICIT NONE
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    LOGICAL ALLOC_REAL
    TYPE(c_ptr) address
    INTEGER(c_int) shmid
    ! local variables
    INTEGER(c_size_t) k
    INTEGER ierror

    IF (ALLOC_REAL) THEN
       k=WDES1%GRID%MPLWV*WDES1%NRSPINORS
       IF (WDES1%COMM_SHMEM%NODE_ME==1) CALL getshmem(INT(8*2*k,KIND=c_size_t),shmid)
       CALL M_bcast_i(WDES1%COMM_SHMEM,shmid,1)
       CALLMPI( MPI_barrier(WDES1%COMM_SHMEM%MPI_COMM,ierror) )
       call attachshmem(shmid, address)
       call c_f_pointer(address,W1%CR,[k])
       ALLOCATE(W1%CPROJ(WDES1%NPROD))
       NULLIFY(W1%CW)
    ELSE
       k=WDES1%NRPLWV
       IF (WDES1%COMM_SHMEM%NODE_ME==1) CALL getshmem(INT(8*2*k,KIND=c_size_t),shmid)
       CALL M_bcast_i(WDES1%COMM_SHMEM,shmid,1)
       CALLMPI( MPI_barrier(WDES1%COMM_SHMEM%MPI_COMM,ierror) )
       call attachshmem(shmid, address)
       call c_f_pointer(address,W1%CPTWFP,[k])
       ALLOCATE(W1%CPROJ(WDES1%NPROD))
       NULLIFY(W1%CR)
    ENDIF
    W1%FERWE=0
    W1%CELEN=0

    W1%WDES1=>WDES1
    ! not a valid band of the Hamiltonian
    W1%NB=-1
    W1%ISP=-1
    W1%LDO=.TRUE.
  END SUBROUTINE NEWWAV_SHMEM
#endif

!=======================================================================
!>  create storage for one orbital W for real space storage
!=======================================================================

  SUBROUTINE NEWWAV_R(W1, WDES1)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER NEWWAV_R
#endif

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    ! local variables
    INTEGER MPLWV

    MPLWV=WDES1%GRID%MPLWV*WDES1%NRSPINORS
    ALLOCATE(W1%CR(MPLWV))
    NULLIFY(W1%CW, W1%CPROJ)
GPACC !$ACC ENTER DATA CREATE(W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
    W1%FERWE=0
    W1%CELEN=0

    W1%WDES1 =>WDES1
    ! not a valid band of the Hamiltonian
    W1%NB=-1
    W1%ISP=-1
    W1%LDO=.TRUE.

GPACC !$ACC UPDATE DEVICE(W1%FERWE,W1%CELEN,W1%NB,W1%ISP,W1%LDO) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
!!      CALL COPYIN_TYPED_VAR(W1%WDES1)
!!      CALL ACC_ATTACH_ASYNC(W1%WDES1,ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(TO: W1, W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
GPOMP !$OMP TARGET UPDATE TO(W1%FERWE, W1%CELEN, W1%NB, W1%ISP, W1%LDO) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
#ifdef CRAYWAR3
GPOMP !$OMP TARGET MAP(PRESENT,TO: W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
         W1%WDES1 =>WDES1
GPOMP !$OMP END TARGET
#endif
  END SUBROUTINE NEWWAV_R

!=======================================================================
!>  create storage for one orbital W for real space storage
!=======================================================================

  SUBROUTINE NEWWAV_NOCW(W1, WDES1)
    IMPLICIT NONE
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER MPLWV

    MPLWV=WDES1%GRID%MPLWV*WDES1%NRSPINORS
    ALLOCATE(W1%CR(MPLWV))
    ALLOCATE(W1%CPROJ(WDES1%NPROD))
    NULLIFY(W1%CW)
    W1%FERWE=0
    W1%CELEN=0

    W1%WDES1 =>WDES1
    ! not a valid band of the Hamiltonian
    W1%NB=-1
    W1%ISP=-1
    W1%LDO=.TRUE.

  END SUBROUTINE NEWWAV_NOCW

!=======================================================================
!>  destroy storage for one orbital W
!>
!>  optionally real space array is deallocated
!=======================================================================

  SUBROUTINE DELWAV(W,DEALLOC_REAL,ISTAT)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELWAV
#endif

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefun1) W
    LOGICAL DEALLOC_REAL
    INTEGER, OPTIONAL :: ISTAT

    IF (PRESENT(ISTAT)) THEN
       IF (DEALLOC_REAL) THEN
GPACC !$ACC EXIT DATA DELETE(W%CPTWFP,W%CPROJ,W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ,W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ,W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
          DEALLOCATE(W%CPTWFP,W%CPROJ,W%CR,STAT=ISTAT)
       ELSE
GPACC !$ACC EXIT DATA DELETE(W%CPTWFP,W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
          DEALLOCATE(W%CPTWFP,W%CPROJ,STAT=ISTAT)
       ENDIF
    ELSE
       IF (DEALLOC_REAL.AND. ASSOCIATED(W%CR)) THEN
GPACC !$ACC EXIT DATA DELETE(W%CPTWFP,W%CPROJ,W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ, W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ, W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
          DEALLOCATE(W%CPTWFP,W%CPROJ,W%CR)
       ELSE IF (DEALLOC_REAL .AND. .NOT. ASSOCIATED(W%CR)) THEN
          WRITE(*,*) 'internal error in DELWAV: should be called with DEALLOC_REAL=.FALSE.'
       ELSE
GPACC !$ACC EXIT DATA DELETE(W%CPTWFP,W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%CPTWFP,W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
          DEALLOCATE(W%CPTWFP,W%CPROJ)
       ENDIF
    ENDIF
GPACC !$ACC EXIT DATA DELETE(W%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
!!      CALL DELETE_TYPED_VAR(W%WDES1)
!!      CALL ACC_DETACH_ASYNC(W%WDES1,ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W, W%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))
    NULLIFY(W%CW,W%CPROJ,W%CR,W%WDES1)
!!$ACC EXIT DATA DELETE(W) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELWAV

#ifdef use_shmem
!=======================================================================
!>  destroy storage for one orbital W
!>
!>  optionally real space array is deallocated
!>  use shmem shared memory
!=======================================================================

  SUBROUTINE DELWAV_SHMEM(W1,WDES1,address,shmid)
    USE c2f_interface
    IMPLICIT NONE
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    LOGICAL DEALLOC_REAL
    TYPE(c_ptr) address
    INTEGER(c_int) shmid
    ! local variables
    INTEGER ierror

    CALL detachshmem(address)
    IF (WDES1%COMM_SHMEM%NODE_ME==1) CALL destroyshmem(shmid)
    CALLMPI( MPI_barrier(WDES1%COMM_SHMEM%MPI_COMM,ierror) )
    DEALLOCATE(W1%CPROJ)

  END SUBROUTINE DELWAV_SHMEM
#endif

  SUBROUTINE DELWAV_R(W)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER DELWAV
#endif

__GPU USE moffload

    IMPLICIT NONE
    TYPE (wavefun1) W
GPACC !$ACC EXIT DATA DELETE(W%CR) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W%CR) IF(OFFLOAD_ON)
    DEALLOCATE(W%CR)
! test_
    IF (ASSOCIATED(W%CPTWFP)) THEN
GPACC !$ACC EXIT DATA DELETE(W%CPTWFP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W%CPTWFP) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W, W%CPTWFP) IF(OFFLOAD_ON)
    ENDIF
    IF (ASSOCIATED(W%CPROJ)) THEN
GPACC !$ACC EXIT DATA DELETE(W%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W%CPROJ) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W, W%CPROJ) IF(OFFLOAD_ON)
    ENDIF
! test_
GPACC !$ACC EXIT DATA DELETE(W%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
!!      CALL DELETE_TYPED_VAR(W%WDES1)
!!      CALL ACC_DETACH_ASYNC(W%WDES1,ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W%WDES1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:W,W%WDES1) IF(OFFLOAD_ON)
    NULLIFY(W%CW,W%CPROJ,W%CR,W%WDES1)
!!$ACC EXIT DATA DELETE(W) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELWAV_R

  SUBROUTINE DELWAV_NOCW(W)
    IMPLICIT NONE
    TYPE (wavefun1) W
    DEALLOCATE(W%CR,W%CPROJ)
    NULLIFY(W%CR,W%CPROJ)
  END SUBROUTINE DELWAV_NOCW

!=======================================================================
!>  set one single orbital (W1) from an array of orbitals
!>
!>  local band index and k point must be supplied
!>  offloading notes:
!>  this routine should not copy W1 to the device if it was not on the
!>  the device! Only do something if W1 was there already
!=======================================================================

  SUBROUTINE SETWAV(W,W1,WDES1,NB,ISP)
#ifdef CRAYWAR1
! in CCE19 we have to use INLINENEVER to avoid inlining
!DIR$ INLINENEVER SETWAV
#endif

__GPU USE moffload

    IMPLICIT NONE
    INTEGER NB,ISP
    TYPE (wavespin) W
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    PROFILING_START('setwav')
! clean up associated pointers if present on device
#ifdef OFFLOADING
    IF (ASSOCIATED(W1%CPTWFP)) THEN
GPACC !$ACC EXIT DATA DELETE(W1%CPTWFP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W1%CPTWFP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W1%CPTWFP) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
    ENDIF
    IF (ASSOCIATED(W1%CPROJ)) THEN
GPACC !$ACC EXIT DATA DELETE(W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W1%CPROJ) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
    ENDIF
    IF (ASSOCIATED(W1%WDES1)) THEN
GPACC !$ACC EXIT DATA DELETE(W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) ASYNC(ACC_ASYNC_Q)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE: W1%WDES1) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W1))
    ENDIF
#endif

    NK=WDES1%NK

    W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
    W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
    W1%FERWE = W%FERWE(NB,NK,ISP)
    W1%CELEN = W%CELEN(NB,NK,ISP)
    W1%WDES1 =>WDES1
    ! remember band index
    ! presently not used
    W1%NB = NB
    W1%ISP = ISP
    W1%LDO = .TRUE.

! now attach new pointers
#ifdef OFFLOADING
    IF (OFFLOAD_ON.AND.IS_ON_DEVICE(W1,1)) THEN
       IF (ASSOCIATED(W1%CPTWFP)) THEN
GPACC !$ACC ENTER DATA COPYIN(W1%CPTWFP) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:W1%CPTWFP) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%CPTWFP) IF(OFFLOAD_ON)
       ENDIF
       IF (ASSOCIATED(W1%CPROJ)) THEN
GPACC !$ACC ENTER DATA COPYIN(W1%CPROJ) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:W1%CPROJ) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%CPROJ) IF(OFFLOAD_ON)
       ENDIF
       IF (ASSOCIATED(W1%WDES1)) THEN
GPACC !$ACC ENTER DATA COPYIN(W1%WDES1) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET ENTER DATA MAP(TO: W1, W1%WDES1) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: W1%WDES1) IF(OFFLOAD_ON)
       ENDIF

! set also FERWE and CELEN on device
GPACC !$ACC UPDATE DEVICE(W1%FERWE,W1%CELEN,W1%NB,W1%ISP,W1%LDO) ASYNC(ACC_ASYNC_Q)
! the following target update is rather slow on intel as of ifx 2025.2
! consider a workaround
GPOMP !$OMP TARGET UPDATE TO(W1%FERWE,W1%CELEN,W1%NB,W1%ISP,W1%LDO) IF(OFFLOAD_ON)
#ifdef OMP_OFFLOAD
#ifdef CRAYHIP
!$OMP TARGET IF(OFFLOAD_ON)
       W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
       W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
       W1%WDES1 =>WDES1
!$OMP END TARGET
#endif
#endif // OMP_OFFLOAD
    ENDIF
#endif // OFFLOADING

  PROFILING_STOP('setwav')

  END SUBROUTINE SETWAV

!=======================================================================
!>  set orbital (wavefun) from an spin array of orbitals
!>  spin must be supplied
!=======================================================================

  SUBROUTINE SETW_SPIN(W,W1,ISPIN)
    USE prec
    IMPLICIT NONE
    INTEGER ISPIN
    TYPE (wavespin) W
    TYPE (wavefun)  W1

    W1%CPTWFP=>W%CPTWFP(:,:,:,ISPIN)
    W1%CPROJ =>W%CPROJ (:,:,:,ISPIN)
    W1%FERWE =>W%FERWE (:,:,ISPIN)
    W1%CELEN =>W%CELEN (:,:,ISPIN)
    W1%AUX   =>W%AUX   (:,:,ISPIN)
    W1%FERTOT=>W%FERTOT(:,:,ISPIN)
    W1%AUXTOT=>W%AUXTOT(:,:,ISPIN)
    W1%CELTOT=>W%CELTOT(:,:,ISPIN)

  END SUBROUTINE SETW_SPIN

!************************* SUBROUTINE WVREAL ***************************
!
!> this subroutine makes a orbital real
!>
!> it is required for the gamma point only mode
!> to avoid that small non real components develop
!
!***********************************************************************

  SUBROUTINE WVREAL(WDES,GRID,W)

__GPU USE moffload

    IMPLICIT REAL(q) (A-H,O-Z)
    TYPE (wavespin) W
    TYPE (wavedes)  WDES
    TYPE (grid_3d)  GRID
#ifdef gammareal
    !
    ! search for position of G=0 zero element
    !
    NK  =1
    NIND=0
    NPL =WDES%NGVECTOR(NK)

    DO N=1,NPL
       IND=WDES%NINDPW(N,NK)
       N1= MOD((IND-1),GRID%RC%NROW)+1
       NC= (IND-1)/GRID%RC%NROW+1
       N2= GRID%RC%I2(NC)
       N3= GRID%RC%I3(NC)
       IF (N1==1 .AND. N2==1 .AND. N3==1) THEN
          NIND=N
          EXIT
       ENDIF
    ENDDO

    IF (NIND /= 0) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W,WDES) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:W,WDES) IF(OFFLOAD_ON)
       DO ISP=1,WDES%ISPIN
          DO NB=1,WDES%NBANDS
             W%CPTWFP(NIND,NB,NK,ISP)= REAL( W%CPTWFP(NIND,NB,NK,ISP) ,KIND=q)
          ENDDO
       ENDDO
       IF (WDES%LNONCOLLINEAR) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W,WDES) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:W,WDES) IF(OFFLOAD_ON)
          DO ISP=1,WDES%ISPIN
             DO NB=1,WDES%NBANDS
                W%CPTWFP(NIND+NPL,NB,NK,ISP)= REAL( W%CPTWFP(NIND+NPL,NB,NK,ISP) ,KIND=q)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

#endif
    RETURN
  END SUBROUTINE WVREAL


!=======================================================================
!
!> NB_LOCAL returns the local storage index of a band
!> if bands are distributed over processors
!
!=======================================================================

  FUNCTION NB_LOCAL(NB,WDES1)
GPACC !$ACC ROUTINE SEQ
    USE prec
    IMPLICIT NONE
GPOMP !$OMP DECLARE TARGET(NB_LOCAL)
    INTEGER,INTENT(IN) :: NB
    TYPE (wavedes1),INTENT(IN) :: WDES1
    INTEGER :: NB_LOCAL

    NB_LOCAL = NB_LOCAL_LOW(NB,WDES1%NB_PAR,WDES1%NB_LOW)
  END FUNCTION NB_LOCAL

!> Low level interface for NB_LOCAL
!> This is just so that one can use it with wavedes or wavedes1
  FUNCTION NB_LOCAL_LOW(NB,NB_PAR,NB_LOW) RESULT(NB_LOCAL)
GPACC !$ACC ROUTINE SEQ
    INTEGER,INTENT(IN) :: NB !< global index of the band
    INTEGER,INTENT(IN) :: NB_PAR !< number of band batches
    INTEGER,INTENT(IN) :: NB_LOW !< index of the lowest band
    INTEGER :: NB_LOCAL !< local index of the global band
GPOMP !$OMP DECLARE TARGET(NB_LOCAL_LOW)

    IF ( MOD(NB-1,NB_PAR)+1 == NB_LOW) THEN
       NB_LOCAL=1+(NB-1)/NB_PAR
    ELSE
       NB_LOCAL=0
    ENDIF
  END FUNCTION NB_LOCAL_LOW

!***************************SUBROUTINE WFINIT***************************
!
!> this subroutine initializes the orbital array W
!>
!> it always uses a random number generator
!> to initialize the coefficients
!
!***********************************************************************

  SUBROUTINE WFINIT(WDES, W, ENINI, NBANDS, NKSTART)
    USE random_seeded, ONLY: DEFAULT_RNG, PCG_32BIT_RNG
    USE tutor, ONLY: vtutor
    TYPE(wavedes)      WDES
    TYPE(wavespin)     W
    REAL(q)            ENINI    !< cutoff energy
    INTEGER, OPTIONAL::NBANDS   !< first band to be initialized
    INTEGER, OPTIONAL::NKSTART  !< first k-point from which to start to generate random WFs
    !
    PROFILING_START('wfinit')
    SELECT CASE (WDES%RANDOM_NUMBER_GENERATOR)
    CASE (DEFAULT_RNG)
       CALL WFINIT_DEFAULT(WDES, W, ENINI, NBANDS, NKSTART)
    CASE (PCG_32BIT_RNG)
       CALL WFINIT_PCG_32(WDES, W, ENINI, NBANDS, NKSTART)
    CASE DEFAULT
       CALL VTUTOR%BUG("Unknown specification of random number generator", __FILE__, __LINE__)
    END SELECT
    PROFILING_STOP('wfinit')
  END SUBROUTINE WFINIT

#define serial_compat
#define random

  SUBROUTINE WFINIT_DEFAULT(WDES, W, ENINI, NBANDS, NKSTART)
    USE constant
    USE random_seeded, ONLY: RANE, RANE_ARRAY
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavespin)    W
    REAL(q)            ENINI    ! cutoff energy
    INTEGER, OPTIONAL::NBANDS   ! first band to be initialized
    INTEGER, OPTIONAL::NKSTART  ! fisrt k-point from which to start to generate random WFs
 !  local
    LOGICAL :: LSCALE
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1
    COMPLEX (q),ALLOCATABLE ::  CW(:)
    COMPLEX (q) :: YY
    INTEGER :: ISPINOR, I, NK, FIRST_BAND, MY_NKSTART, TOTAL_RAND_NEEDED, NRAND_PER_SPINOR, OFF
    REAL(q) :: MAX_EIGENVALUE
    REAL(q), ALLOCATABLE :: RANDOM_ARR(:)

    MY_NKSTART = 1; IF (PRESENT(NKSTART)) MY_NKSTART = NKSTART

    IF (PRESENT(NBANDS)) THEN
      FIRST_BAND=NBANDS
      MAX_EIGENVALUE=MAXVAL(REAL(W%CELTOT,q))+2
    ELSE
      FIRST_BAND=1
      W%CELTOT(:,MY_NKSTART:W%WDES%NKPTS,:)=0
      W%CPTWFP(:,:,MY_NKSTART:W%WDES%NKPTS,:)=0
      MAX_EIGENVALUE=0
    ENDIF

    NALLOC = MAXVAL(WDES%NPLWKP_TOT)
    ALLOCATE(CW(NALLOC))

    ! if all band are calculated, use non weighted random numbers
    ! for the initialization; this avoids linear dependencies
    IF (MINVAL(WDES%NB_TOTK(1:WDES%NKPTS,:))< WDES%NB_TOT) THEN
       LSCALE=.FALSE.
    ELSE
       LSCALE=.TRUE.
    ENDIF

    ! determine once the maximal numbers of random numbers drawn
    NPL_GLOBAL=NALLOC/WDES%NRSPINORS
#if defined(random) || defined(gammareal)
    ! draw 2*NPL_GLOBAL-1 random numbers (real + cmplx parts) for each spinor
    NRAND_PER_SPINOR = 2*NPL_GLOBAL - 1
#else
    ! draw NPL_GLOBAL random numbers (real only) for each spinor
    NRAND_PER_SPINOR = NPL_GLOBAL
#endif
    TOTAL_RAND_NEEDED=WDES%NRSPINORS*NRAND_PER_SPINOR
    ALLOCATE(RANDOM_ARR(TOTAL_RAND_NEEDED))

#if defined(serial_compat) && defined(MPI)
    spin:   DO I=1,WDES%ISPIN
       YR=RANE()  ! what is that ? it is here for compatibility
       kpoint: DO NK=MY_NKSTART,WDES%NKPTS

          CALL SETWDES(WDES,WDES1,NK)

          ! determine random numbers for this k-point to be drawn
          NPL_GLOBAL=WDES%NPLWKP_TOT(NK)/WDES%NRSPINORS

#if defined(random) || defined(gammareal)
          ! draw 2*NPL_GLOBAL-1 random numbers (real + cmplx parts) for each spinor
          NRAND_PER_SPINOR = 2*NPL_GLOBAL - 1
#else
          ! draw NPL_GLOBAL random numbers (real only) for each spinor
          NRAND_PER_SPINOR = NPL_GLOBAL
#endif
          TOTAL_RAND_NEEDED=WDES%NRSPINORS*NRAND_PER_SPINOR

          band: DO NB_GLOBAL=FIRST_BAND,WDES%NB_TOT

             NB=NB_LOCAL(NB_GLOBAL,WDES1)
!=======================================================================
!     initialize for serial compatible version
!=======================================================================
             CW=0  ! just in case initialize CW

             ! get all random numbers for all spinors at once for this band
             CALL RANE_ARRAY(RANDOM_ARR, TOTAL_RAND_NEEDED)

             DO ISPINOR = 0, WDES%NRSPINORS-1
                OFF = ISPINOR * NRAND_PER_SPINOR
                MM_START = 1 + NPL_GLOBAL*ISPINOR
                MM_END   = NPL_GLOBAL + NPL_GLOBAL*ISPINOR

#if defined(random) || defined(gammareal)
                ! First element use only real part
                CW(MM_START) = RANDOM_ARR(OFF + 1)

                ! Remaining interleaved to match original scalar version
                DO M = 2, NPL_GLOBAL
                   CW(MM_START + M - 1) = CMPLX( &
                                                RANDOM_ARR(OFF + 2*M - 2), &        ! real part
                                                0.2_Q * RANDOM_ARR(OFF + 2*M - 1) - 0.1_q, &  ! imag part
                                                Q )
                END DO
#else
                CW(MM_START:MM_END) = RANDOM_ARR(1+OFF:NRAND_PER_SPINOR+OFF)
#endif

             END DO


             CALL DIS_PW_BAND(WDES1,NB_GLOBAL,CW,W%CPTWFP(1,1,NK,I))
             IF (NB==0) CYCLE band  ! not on local node, than cycle

             CALL WEIGHT_BY_KIN_ENERGY(WDES, W%CPTWFP(:,NB,NK,I), ENINI, NK, LSCALE)

!=======================================================================
!     initialize components of orbital with random numbers
!     weighted with the inverse of the kinetic energy of the plane
!     wave (random phase factor added gK)
!=======================================================================
#else
    spin:   DO I=1,WDES%ISPIN
       YR=RANE()  ! what is that ? it is here for compatibility
       kpoint: DO NK=MY_NKSTART,WDES%NKPTS

          CALL SETWDES(WDES,WDES1,NK)

          NPL=WDES%NGVECTOR(NK)

          band: DO NB=FIRST_BAND,WDES%NBANDS
             spinor: DO ISPINOR=0,WDES%NRSPINORS-1 ! to be included later, Testing only
                DO M=1,NPL
                   YY=RANE()
                   !  at the gamma it is somethimes better to use random phase factor
                   !  (!!! but if the cell has inversion symmtry it is a bad choice !!!)
#if defined(random) || defined(gammareal)
                   IF (M/=1) YY=CMPLX(REAL(YY,q),0.2_q*RANE()-0.1_q,q)
#endif

                   IF(WDES%DATAKE(M,ISPINOR+1,NK)>=ENINI) THEN
                      YY=0
                   ENDIF
                   WW=WDES%DATAKE(M,ISPINOR+1,NK)
                   IF(WW<=0.01_q) WW=0.1_q
                   IF(.NOT.LSCALE) WW=1._q
                   YY=YY/WW
                   W%CPTWFP(M+ NPL*ISPINOR, NB, NK, I)=YY
                ENDDO
             ENDDO spinor
#endif
             CALL NORMALIZE(WDES%COMM_INB, W%CPTWFP(:,NB,NK,I))

             ! set CELEN to kinetic energy
             ! W%CELEN is hardly used, but the occupancies are sometimes
             ! initialized from them
             ! if additional randomly initialized bands are added they most
             ! lie well above the occupied bands
             W%CELEN(NB,NK,I) = EVALUATE_KIN_ENERGY(WDES, W%CPTWFP(:,NB,NK,I), NK, MAX_EIGENVALUE)

          ENDDO band
       ENDDO kpoint
    ENDDO spin

    CALL WFZERO(W)
    CALL MRG_CEL(WDES,W)
#ifdef debug
#ifdef MPI
    IF (WDES%COMM%NODE_ME ==  WDES%COMM%IONODE) THEN
       WRITE(*,'(I4,F14.10)') (I,REAL(W%CELTOT(I,1,1),q),I=1,WDES%NB_TOT)
    ENDIF
#else
    WRITE(*,'(I4,F14.10)') (I,REAL(W%CELTOT(I,1,1),q),I=1,WDES%NB_TOT)
#endif
#endif

    DEALLOCATE(CW)
    IF (ALLOCATED(RANDOM_ARR)) THEN
       DEALLOCATE(RANDOM_ARR)
    END IF
    RETURN
  END SUBROUTINE

  !> Initialize the orbitals with a parallel random number generator.
  !>
  !> Every spin, k point, and band uses a different random number stream so that the
  !> initialization can be parallelized over the MPI ranks
  SUBROUTINE WFINIT_PCG_32(WDES, W, ENINI, NBANDS, NKSTART)
    USE constant
    USE random_seeded, ONLY: pcg_generator, NEW_PCG_GENERATOR, RANDOM_NUMBER

    IMPLICIT NONE

    TYPE(wavedes), INTENT(INOUT) :: WDES  !< descriptor of the wavefunction
    TYPE(wavespin), INTENT(INOUT) :: W  !< All orbitals used by this rank
    REAL(q), INTENT(IN) :: ENINI  !< cutoff energy
    INTEGER, INTENT(IN), OPTIONAL :: NBANDS  !< first band to be initialized
    INTEGER, INTENT(IN), OPTIONAL :: NKSTART  !< first k-point from which to start to generate random WFs

    TYPE(pcg_generator) GENERATOR
    TYPE(wavedes1) WDES1
    LOGICAL LSCALE
    COMPLEX(q), ALLOCATABLE :: CW(:)
    INTEGER SPIN, KPOINT, BAND, BAND_GLOBAL, FIRST_BAND, FIRST_KPOINT, SPINOR
    INTEGER TASK, NUM_TASK, II, GG, NUM_ELEMENT
    INTEGER, ALLOCATABLE :: TASKS(:,:)
    INTEGER(qi8) STREAM
    REAL(q) MAX_EIGENVALUE, KIN_ENERGY, NORM
    REAL(q), ALLOCATABLE :: RE_PART(:), IM_PART(:)

    FIRST_KPOINT = 1
    IF (PRESENT(NKSTART)) FIRST_KPOINT = NKSTART

    IF (PRESENT(NBANDS)) THEN
       FIRST_BAND = NBANDS
       MAX_EIGENVALUE = MAXVAL(REAL(W%CELTOT,q)) + 2
    ELSE
       FIRST_BAND = 1
       W%CELTOT(:,FIRST_KPOINT:W%WDES%NKPTS,:) = 0
       MAX_EIGENVALUE = 0
    END IF

    ALLOCATE(TASKS(3, WDES%ISPIN * WDES%NKPTS * WDES%NB_TOT))
    TASK = 0
    DO SPIN = 1, WDES%ISPIN
       DO KPOINT = FIRST_KPOINT, WDES%NKPTS
          IF (MOD(KPOINT - 1, WDES%COMM_KINTER%NCPU) + 1 /= WDES%COMM_KINTER%NODE_ME) CYCLE
          CALL SETWDES(WDES, WDES1, KPOINT)
          DO BAND_GLOBAL = FIRST_BAND, WDES%NB_TOT
             BAND = NB_LOCAL(BAND_GLOBAL, WDES1)
             IF (BAND == 0) CYCLE
             TASK = TASK + 1
             TASKS(:, TASK) = [BAND_GLOBAL, KPOINT, SPIN]
          END DO
       END DO
    END DO
    NUM_TASK = TASK

    ! if all band are calculated, use non weighted random numbers
    ! for the initialization; this avoids linear dependencies
    LSCALE = MINVAL(WDES%NB_TOTK(1:WDES%NKPTS,:)) >= WDES%NB_TOT

    ALLOCATE(CW(MAXVAL(WDES%NPLWKP_TOT)))
    ALLOCATE(RE_PART(SIZE(CW)), IM_PART(SIZE(CW)))

    !$OMP PARALLEL DO DEFAULT(PRIVATE) &
    !$OMP SHARED(W, WDES, TASKS, NUM_TASK, ENINI, MAX_EIGENVALUE, LSCALE)
    DO TASK = 1, NUM_TASK
       BAND_GLOBAL = TASKS(1,TASK)
       KPOINT = TASKS(2,TASK)
       SPIN = TASKS(3,TASK)

       CALL SETWDES(WDES, WDES1, KPOINT)
       NUM_ELEMENT = WDES%NPLWKP_TOT(KPOINT)
       BAND = NB_LOCAL(BAND_GLOBAL, WDES1)

       ! initialize components of orbital with random numbers
       STREAM = (SPIN * WDES%NKPTS + KPOINT) * WDES%NB_TOT + BAND_GLOBAL
       GENERATOR = NEW_PCG_GENERATOR(WDES%PCG_SEED, STREAM)
       !
       CALL RANDOM_NUMBER(GENERATOR, RE_PART(:NUM_ELEMENT))
       !
       CALL RANDOM_NUMBER(GENERATOR, IM_PART(:NUM_ELEMENT))
       IM_PART = 0.2_q * IM_PART - 0.1_q
       ! at the gamma it is somethimes better to use random phase factor
       ! (!!! but if the cell has inversion symmtry it is a bad choice !!!)
       DO II = 1, NUM_ELEMENT, NUM_ELEMENT / WDES%NRSPINORS
          IM_PART(II) = 0
       END DO

       CW(:NUM_ELEMENT) = CMPLX(RE_PART(:NUM_ELEMENT), IM_PART(:NUM_ELEMENT), KIND=q)
       CW(NUM_ELEMENT+1:) = 0
       CALL DIS_PW(WDES1, CW, W%CPTWFP(1,BAND,KPOINT,SPIN))
       W%CPTWFP(WDES1%NPL+1:,BAND,KPOINT,SPIN) = 0

       CALL WEIGHT_BY_KIN_ENERGY(WDES, W%CPTWFP(:,BAND,KPOINT,SPIN), ENINI, KPOINT, LSCALE)
       CALL NORMALIZE(WDES%COMM_INB, W%CPTWFP(:,BAND,KPOINT,SPIN))

       ! set CELEN to kinetic energy
       ! W%CELEN is hardly used, but the occupancies are sometimes
       ! initialized from them
       ! if additional randomly initialized bands are added they most
       ! lie well above the occupied bands
       W%CELEN(BAND,KPOINT,SPIN) = EVALUATE_KIN_ENERGY(WDES, W%CPTWFP(:,BAND,KPOINT,SPIN), KPOINT, MAX_EIGENVALUE)

    END DO ! task
    !$OMP END PARALLEL DO

    CALL WFZERO(W)
    CALL MRG_CEL(WDES, W)

  END SUBROUTINE WFINIT_PCG_32

  !> Weight the components by the inverse of the kinetic energy
  !> and zero them above the cutoff
  !>
  !> If the kinetic energy is above the cutoff zero the component.
  !> If the kinetic energy is too small use a minimal weight.
  SUBROUTINE WEIGHT_BY_KIN_ENERGY(WDES, COEFF, ENINI, KPOINT, LSCALE)
    IMPLICIT NONE
    TYPE(wavedes), INTENT(IN) :: WDES
    COMPLEX(q), INTENT(INOUT) :: COEFF(:)
    REAL(q), INTENT(IN) :: ENINI  !< energy cutoff
    INTEGER, INTENT(IN) :: KPOINT !< k-point index
    LOGICAL, INTENT(IN) :: LSCALE !< should the energy be scaled or just zeroed
    REAL(q) KIN_ENERGY
    INTEGER SPINOR, GG, II
    !
    DO SPINOR = 1, WDES%NRSPINORS
       DO GG = 1, WDES%NGVECTOR(KPOINT)
          II = GG + WDES%NGVECTOR(KPOINT) * (SPINOR - 1)
          KIN_ENERGY = WDES%DATAKE(GG,SPINOR,KPOINT)
          IF (KIN_ENERGY >= ENINI) THEN
             COEFF(II) = 0
          ELSE IF (LSCALE) THEN
             IF (KIN_ENERGY < 0.01_q) KIN_ENERGY = 0.1_q
             COEFF(II) = COEFF(II) / KIN_ENERGY
          END IF
       END DO
    END DO
  END SUBROUTINE WEIGHT_BY_KIN_ENERGY

  !> Check orbital has nonvanishing norm and normalize it afterwards
  SUBROUTINE NORMALIZE(COMM, COEFF)
     IMPLICIT NONE
     TYPE(communic), INTENT(IN) :: COMM
     COMPLEX(q), INTENT(INOUT) :: COEFF(:)
     REAL(q) NORM
     !
     NORM = DOT_PRODUCT(COEFF, COEFF)
     CALLMPI( M_sum_d(COMM, NORM, 1))
     IF (NORM <= 0.000001_q) THEN
        CALL vtutor%bug("internal error in WFINIT_PCG_32: orbitals linearily dependent at " &
           // "random-number initialization", __FILE__, __LINE__)
     END IF
     COEFF = COEFF / SQRT(NORM)
  END SUBROUTINE NORMALIZE

  !> Evaluate the kinetic energy for a particular orbital
  REAL(q) FUNCTION EVALUATE_KIN_ENERGY(WDES, COEFF, KPOINT, MAX_EIGENVALUE) RESULT (KIN_ENERGY)
     IMPLICIT NONE
     TYPE(wavedes), INTENT(IN) :: WDES
     COMPLEX(q), INTENT(IN) :: COEFF(:)
     INTEGER, INTENT(IN) :: KPOINT
     REAL(q), INTENT(IN) :: MAX_EIGENVALUE
     INTEGER SPINOR, GG, II
     !
     KIN_ENERGY = 0
     DO SPINOR = 1, WDES%NRSPINORS
        DO GG = 1, WDES%NGVECTOR(KPOINT)
           II = GG + WDES%NGVECTOR(KPOINT) * (SPINOR - 1)
           KIN_ENERGY = KIN_ENERGY + COEFF(II) * CONJG(COEFF(II)) * WDES%DATAKE(GG,SPINOR,KPOINT)
        END DO
     END DO
     CALLMPI( M_sum_d(WDES%COMM_INB, KIN_ENERGY, 1))
     KIN_ENERGY = MAX(MAX_EIGENVALUE, KIN_ENERGY)
  END FUNCTION EVALUATE_KIN_ENERGY


!***************************SUBROUTINE WFZERO **************************
!
!> this subroutine sets orbitals beyond W%WDES%NB_TOTK(NK, ISP)
!> to zero
!> this should avoid any problems in the DOS routines
!> at least as long as they are not WDES%NB_TOTK clean
!
!***********************************************************************

  SUBROUTINE WFZERO(W)
    USE constant
    IMPLICIT NONE
    TYPE (wavespin)    W
 !  local
    INTEGER ISP, NK, NB

    DO ISP=1,W%WDES%ISPIN
       DO NK=1,W%WDES%NKPTS

          ! first local band set to zero
          DO NB=(W%WDES%NB_TOTK(NK,ISP)-W%WDES%NB_LOW)/W%WDES%NB_PAR+2, W%WDES%NBANDS
             W%CW(:,NB,NK,ISP)=0
             W%CPROJ(:,NB,NK,ISP)=0
             W%FERWE(NB,NK,ISP)=0
             W%CELEN(NB,NK,ISP)=1E4     ! huge eigenvalue
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE WFZERO


!****************SUBROUTINE WFSET_HIGH************************
!
!> the subroutine sets eigenvalues to very large values
!> this should avoid any problems in the DOS routines
!> at least as long as they are not WDES%NB_TOTK clean
!
!*************************************************************
  SUBROUTINE WFSET_HIGH_CELEN(W)

__GPU USE moffload_struct_def

    USE constant
    IMPLICIT NONE
    TYPE (wavespin)    W
 !  local
    INTEGER ISP, NK, NB

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(NB) PRESENT(W) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(NB) &
GPOMP !$OMP MAP(PRESENT,TO:W) IF(OFFLOAD_ON)
    DO ISP=1,W%WDES%ISPIN
       DO NK=1,W%WDES%NKPTS
          DO NB=W%WDES%NB_TOTK(NK,ISP)+1, W%WDES%NB_TOT
             W%FERTOT(NB,NK,ISP)=0         ! occupancy zero
             W%CELTOT(NB,NK,ISP)=10000     ! huge eigenvalue
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE WFSET_HIGH_CELEN


!*************************** FUNCTION isMetal **************************
!
!***********************************************************************

  LOGICAL FUNCTION isMetal(occupations)
    REAL(q), INTENT(IN) :: occupations(:,:,:)
    ! local variables
    INTEGER, ALLOCATABLE :: numberOccupied(:)
    INTEGER :: ispin
    CALL determineNumberOccupied(occupations,numberOccupied)
    isMetal = .FALSE.
    DO ispin = 1, SIZE(numberOccupied)
       isMetal = isMetal .OR. ANY(.NOT.isFilled(occupations(numberOccupied(ispin),:,ispin)))
       IF (numberOccupied(ispin)<SIZE(occupations,1)) THEN
          isMetal = isMetal .OR. ANY(.NOT.isEmpty(occupations(numberOccupied(ispin)+1,:,ispin)))
       ENDIF
    ENDDO
    DEALLOCATE(numberOccupied)
  END FUNCTION isMetal


!*************************** SUBROUTINE determineNumberOccupied ********
!> subroutine determineNumberOccupied
!***********************************************************************

  SUBROUTINE determineNumberOccupied(occupations, numberOccupied)
    REAL(q), INTENT(IN) :: occupations(:,:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: numberOccupied(:)
    ! local variables
    INTEGER :: ispin
    ALLOCATE(numberOccupied(size(occupations, 3)))
    DO ispin = 1, SIZE(numberOccupied)
       numberOccupied(ispin) = COUNT(.NOT.isEmpty(occupations(:,1,ispin)))
    ENDDO
  END SUBROUTINE determineNumberOccupied


!*************************** FUNCTION isEmpty **************************
!
!***********************************************************************

  LOGICAL ELEMENTAL FUNCTION isEmpty(occupation)
    REAL(q), INTENT(IN) :: occupation
    isEmpty = ABS(occupation) < tiny_occupation
  END FUNCTION isEmpty


!*************************** FUNCTION isFilled *************************
!
!***********************************************************************

  LOGICAL ELEMENTAL FUNCTION isFilled(occupation)
    REAL(q), INTENT(IN) :: occupation
    isFilled = ABS(1.0_q-occupation) < tiny_occupation
  END FUNCTION isFilled


!*************************SUBROUTINE GEN_LAYOUT**************************
!
!> subroutine GEN_LAYOUT performs a number of tasks:
!> -  determines the layout (distribution) of the columns on parallel
!>   computers
!>   essentially these are the structures GRID%RC, GRID%RL
!>   and in the parallel version the array GRID%IN
!>
!>   also the following  arrays are allocated:
!>     WDES%NPLWKP  WDES%NINDPW
!> - for LSETUP=.TRUE. the kinetic energy arrays and the G-vector array
!>   are additionally allocated
!>     WDES%DATAKE, WDES%IGX, Y, Z
!> - setup of these arrays is performed in the routines GEN_INDEX
!>
!> the data layout is based on the initial reciprocal lattice vectors
!> stored in BI
!>
!***********************************************************************

  SUBROUTINE GEN_LAYOUT(GRID,WDES, B,BI,IU6,LSETUP,LNOGAMMA,LREAL,NGDIM)
    USE prec
    USE mgrid
    USE constant
    USE base
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d), TARGET :: GRID
    TYPE (wavedes)     WDES
    DIMENSION B(3,3),BI(3,3)          ! current lattice, and initial lattice
    LOGICAL LSETUP
    LOGICAL, OPTIONAL :: LNOGAMMA     ! if LNOGAMMA is set than complex to complex FFT is used (regardless of gammareal preprocessor flag)
    LOGICAL, OPTIONAL :: LREAL        ! if LREAL is set and we use gamma point only, the orbital are stored real valued
    INTEGER, OPTIONAL :: NGDIM
    ! local
    LOGICAL LUP

    LOGICAL, ALLOCATABLE :: LUSE_IN(:)
    INTEGER, ALLOCATABLE :: USED_ROWS(:),IND2(:),IND3(:)
    INTEGER, ALLOCATABLE :: REDISTRIBUTION_INDEX(:)
    COMMON /WAVCUT/   IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX
    REAL(q), PARAMETER :: G2ZERO=1E-12_q

#ifdef gammareal
    WDES%LGAMMA=.TRUE.
    GRID%REAL2CPLX=.TRUE.
#else
    WDES%LGAMMA=.FALSE.
    GRID%REAL2CPLX=.FALSE.
#endif
!gKnew
    ! force complex to complex FFT
    IF (PRESENT(LNOGAMMA)) THEN
       IF ( LNOGAMMA ) THEN
          WDES%LGAMMA=.FALSE.
          GRID%REAL2CPLX=.FALSE.
    ! force complex to real FFT and force real storage mode in real space
       ELSE
          WDES%LGAMMA=.TRUE.
          GRID%REAL2CPLX=.TRUE.
       ENDIF
    ENDIF
    ! orbitals are usually complex in the direct grid in VASP hence LREAL defaults to .FALSE.
    GRID%LREAL=.FALSE.
    ! if LREAL is present and we are at the Gamma point use that as default
    IF (PRESENT(LREAL) .AND. WDES%LGAMMA) THEN
       GRID%LREAL=LREAL
    ENDIF
!gKnew
    WDES%GRID=> GRID

    GRID%NGX_rd=GRID%NGX
    GRID%NGY_rd=GRID%NGY
    GRID%NGZ_rd=GRID%NGZ


    CALL NULLIFY_GRD(GRID)

    WDES%NCOL=0
!=======================================================================
! determine the layout
! i.e. all required columns
! or (y,z) pairs which are required for the 3d-FFT grid
!=======================================================================
#ifdef MPI
    IF (GRID%COMM%NCPU/=1 .OR. LUSE_PARALLEL_FFT) THEN
       IF (WDES%LGAMMA) THEN
          GRID%NGZ_rd=GRID%NGZ/2+1
          IF (IU6>=0) THEN
             WRITE(IU6,*) 'use parallel FFT for orbitals z direction half grid'
          ENDIF
       ENDIF
!-----------------------------------------------------------------------
! parallel version
! reciprocal space first
!-----------------------------------------------------------------------

       ! determine all required grid points (usually a circle around center)

       ALLOCATE(USED_ROWS(GRID%NGY*GRID%NGZ_rd), &
            IND2(GRID%NGY*GRID%NGZ_rd),IND3(GRID%NGY*GRID%NGZ_rd), &
            REDISTRIBUTION_INDEX(GRID%NGY*GRID%NGZ_rd), &
            LUSE_IN(GRID%NGZ_rd))
       ! following arrays are required to setup the reciprocal layout
       USED_ROWS=0     ! gives the total number of points in each column
       ! contained within the cutoff sphere
       NCOL_TOT=0      ! total number of columns
       IND2=0          ! y index of each used column
       IND3=0          ! z index of each used column
       ! LUSE_IN array is used to find out which z-planes must be contained
       ! in the intermediate representation
       LUSE_IN=.FALSE.

       DO N3=1,GRID%NGZ_rd
          DO N2=1,GRID%NGY

             IUSED_MAX=0   ! maxmimal number of rows for this column

             ! loop over all k-points in the *entire* or irreducible Brillouin zone (N)
             ! NKPTS_FOR_GEN_LAYOUT is set, if HF is used, or if the number of k-points
             ! is allowed to change dynamically (linear response, finite differences)
             DO NK=1,WDES%NKPTS_FOR_GEN_LAYOUT
                G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
                G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
                IUSED=0     ! number of rows for this column and k-point
                DO N1=1,GRID%NGX_rd
                   G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
                   GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
                   GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
                   GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI
                   ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

                   ! exclude some components for gamma-only version (C(G)=C*(-G))
                   IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                      IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                      IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
                   ENDIF

                   IF(ENERGI<WDES%ENMAX) THEN
                      IUSED=IUSED+1
                      LUSE_IN(N3)=.TRUE.
                   ENDIF
                ENDDO
                IUSED_MAX=MAX(IUSED,IUSED_MAX)
             ENDDO
             ! if this column is used, enter it in required array

             IF (IUSED_MAX>0) THEN
                NCOL_TOT=NCOL_TOT+1        ! increase column counter
                USED_ROWS(NCOL_TOT)=IUSED  ! set number of rows in this column
                IND2(NCOL_TOT)=N2          ! set y index of this column
                IND3(NCOL_TOT)=N3          ! set z index of this column
             ENDIF
          ENDDO
       ENDDO
       ! setup redistribution array
       ! by sorting according to the number of rows in each column
       ! we try to get better load balancing
       DO I=1,NCOL_TOT
          REDISTRIBUTION_INDEX(I)=I
       ENDDO
       ! do not sort column 1 (it must store the G vector G=(0,0,0)) !
       CALL SORT_REDIS(NCOL_TOT,USED_ROWS(1),REDISTRIBUTION_INDEX(1))

       ! distribute among processors

       NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU
       ALLOCATE(GRID%RC%I2( NCOL_MAX ))
       ALLOCATE(GRID%RC%I3( NCOL_MAX ))

       GRID%RC%NFAST= 1
       GRID%RC%NROW = GRID%NGX_rd
       GRID%RC%NCOL = 0

       IND=0
       !
       ! step through all columns and distribute them onto proc.
       ! in the manner 1 ... NCPU - NCPU ... 1 - 1 ... NCPU - etc.
       !
       NODE_TARGET=0  ! NODE onto which column has to go
       LUP=.TRUE.     ! determines whether NODE_TARGET is increased or decreased

       DO NCOL=1,NCOL_TOT
          IND_REDIS=REDISTRIBUTION_INDEX(NCOL)
          N2=IND2(IND_REDIS)
          N3=IND3(IND_REDIS)
          IF (LUP) THEN
             IF (NODE_TARGET == GRID%COMM%NCPU) THEN
                LUP=.FALSE.
             ELSE
                NODE_TARGET=NODE_TARGET+1
             ENDIF
          ELSE
             IF (NODE_TARGET == 1) THEN
                LUP=.TRUE.
             ELSE
                NODE_TARGET=NODE_TARGET-1
             ENDIF
          ENDIF

          IND_ON_CPU=(NCOL-1)/GRID%COMM%NCPU+1

          IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
             GRID%RC%NCOL=GRID%RC%NCOL+1
             IF (IND_ON_CPU /= GRID%RC%NCOL) THEN
                CALL vtutor%bug("GENSP: internal error(1) " // str(GRID%COMM%NODE_ME) // " " // &
                   str(IND_ON_CPU) // " " // str(GRID%RC%NCOL), __FILE__, __LINE__)
             ENDIF
             GRID%RC%I2(IND_ON_CPU)=N2
             GRID%RC%I3(IND_ON_CPU)=N3
          ENDIF
       ENDDO
       GRID%RC%NP    = GRID%RC%NCOL*GRID%RC%NROW
       GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
       WDES%NCOL     = GRID%RC%NCOL
       !
       ! set up intermediate grid
       ! this usually contains a certain number of z-planes around the center
       !
       LUSE_IN=.TRUE.

       NCOL_TOT=0
       DO N3=1,GRID%NGX_rd
          DO N2=1,GRID%NGZ_rd
             IF (LUSE_IN(N2)) NCOL_TOT=NCOL_TOT+1
          ENDDO
       ENDDO

       NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU

       IF (LPLANE_WISE) THEN
          NCOL_MAX=GRID%NGZ*((GRID%NGX+GRID%COMM%NCPU-1)/GRID%COMM%NCPU)
       ENDIF

       ALLOCATE(GRID%IN%I2( NCOL_MAX ))
       ALLOCATE(GRID%IN%I3( NCOL_MAX ))

       GRID%IN%NFAST= 2
       GRID%IN%NCOL = 0
       GRID%IN%NROW = GRID%NGY

       IND=0
       DO N3=1,GRID%NGX_rd
          DO N2=1,GRID%NGZ_rd
             IF (LUSE_IN(N2)) THEN
                IF (LPLANE_WISE) THEN
                   NODE_TARGET=MOD(N3-1,GRID%COMM%NCPU)+1
                ELSE
                   NODE_TARGET=MOD(IND,GRID%COMM%NCPU)+1
                ENDIF
                IND=IND+1
                IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
                   N=GRID%IN%NCOL
                   N=N+1
                   GRID%IN%NCOL=N
                   GRID%IN%I2(N)=N2 ! I2 contains z index
                   GRID%IN%I3(N)=N3 ! I3      the x index
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF (GRID%IN%NCOL >NCOL_MAX) THEN
          CALL vtutor%bug("GEN_LAYOUT: internal error " // str(NCOL_MAX) // " " // str(GRID%IN%NCOL), __FILE__, __LINE__)
       ENDIF
       GRID%IN%NP    = GRID%IN%NCOL*GRID%IN%NROW
       GRID%IN%NALLOC= GRID%IN%NCOL*GRID%IN%NROW
       !
       ! set up real space grid, use all points
       ! in this case GRID will be conformable to GRID_SOFT in real space
       !
       CALL REAL_STDLAY(GRID, LPARALLEL=.TRUE.)

       DEALLOCATE(USED_ROWS,IND2,IND3,REDISTRIBUTION_INDEX,LUSE_IN)
    ELSE
#endif
!-----------------------------------------------------------------------
! serial version
! set reciprocal space and real space layout for non parallel computers
! always x-first (or x-fast) layout
! all grid points are used for FFT
!-----------------------------------------------------------------------
       IF (WDES%LGAMMA) THEN
          GRID%NGX_rd=GRID%NGX/2+1
          IF (IU6>=0) WRITE(IU6,*) 'use serial FFT for orbitals x direction half grid'
       ENDIF

       GRID%RC%NFAST= 1
       GRID%RC%NCOL = GRID%NGZ_rd*GRID%NGY
       WDES%NCOL    = GRID%RC%NCOL
       GRID%RC%NROW = GRID%NGX_rd
       GRID%RC%NP   = GRID%RC%NCOL*GRID%RC%NROW
       GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
       ALLOCATE(GRID%RC%I2( GRID%RC%NCOL ))
       ALLOCATE(GRID%RC%I3( GRID%RC%NCOL ))
       IND=1
       DO N3=1,GRID%NGZ_rd
          DO N2=1,GRID%NGY
             GRID%RC%I2(IND)=N2
             GRID%RC%I3(IND)=N3
             IND=IND+1
          ENDDO
       ENDDO

       CALL REAL_STDLAY(GRID, LPARALLEL=.FALSE.)
#ifdef MPI
    ENDIF
#endif
!=======================================================================
! count number of plane wave coefficients
! and allocate required arrays
!=======================================================================
    WDES%NRPLWV=WDES_COUNT_INDPW(WDES,GRID,B,BI)
    IF (PRESENT(NGDIM)) THEN
        IF (NGDIM>0) THEN
            IF (NGDIM<WDES%NRPLWV) THEN
                call vtutor%bug("Optional NGDIM "//str(NGDIM)//" passed to GEN_LAYOUT is smaller than the one computed inside "//&
                                str(WDES%NRPLWV),__FILE__,__LINE__)
            ENDIF
            WDES%NRPLWV = NGDIM
        ENDIF
    ENDIF
    WDES%NGDIM=WDES%NRPLWV
    WDES%NRPLWV = WDES%NRPLWV*WDES%NRSPINORS
    WDES%NRPLWV_RED=WDES%NRPLWV/WDES%NB_PAR

    WDES%DO_REDIS = .FALSE.
    IF (WDES%NB_PAR /=1) THEN
       WDES%DO_REDIS=.TRUE.
    ENDIF

    !    CALL  MAKE_STRIDE(WDES%NRPLWV)

    GRID%MPLWV=__align_DIM(MAX(GRID%RC%NALLOC ,GRID%IN%NALLOC , GRID%RL%NALLOC))

#ifdef MPI
    NODE_ME=GRID%COMM%NODE_ME
    IONODE =GRID%COMM%IONODE
#endif
    DWRITE 'gen_layout',NODE_ME,GRID%RC%NALLOC ,GRID%IN%NALLOC , GRID%RL%NALLOC

    CALL ALLOCWDES(WDES,LSETUP)
    WDES%NPLWKP=0
    WDES%NGVECTOR=0
    WDES%NGVECTOR_POS=0

    IF      (WDES%ISPIN==1  .AND. .NOT. WDES%LNONCOLLINEAR ) THEN
       WDES%NCDIJ=1
    ELSE IF (WDES%ISPIN==2) THEN
       WDES%NCDIJ=2
    ELSE IF (WDES%ISPIN==1  .AND. WDES%LNONCOLLINEAR ) THEN
       WDES%NCDIJ=4
    ELSE
       CALL vtutor%bug("internal error: can not set NCDIJ", __FILE__, __LINE__)
    ENDIF

! set WDES%AT_GAMMA
! this flag determines whether a specific k-point corresponds to Gamma
! if this is the case the code might do some special treatment allowed
! only at Gamma (e.g. orbitals are real)
! presently this is only supported in few places
    DO NK=1,WDES%NKDIM
       IF (ABS(SUM(WDES%VKPT(:,NK)**2)) < G2ZERO .AND. .NOT. WDES%LNONCOLLINEAR) THEN
          WDES%AT_GAMMA(NK)=.TRUE.
       ELSE
          WDES%AT_GAMMA(NK)=.FALSE.
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE GEN_LAYOUT

!*********************SUBROUTINE WDES_COUNT_INDPW **********************
!> Count the number of plane-waves in the wavefunction descriptor
!***********************************************************************
  FUNCTION WDES_COUNT_INDPW(WDES,GRID,B,BI) RESULT(NRPLWV)
    IMPLICIT NONE
    TYPE (wavedes) :: WDES
    TYPE (grid_3d), TARGET :: GRID
    REAL(q) :: B(3,3), BI(3,3)
    INTEGER :: NRPLWV
    ! local variables
    INTEGER :: NK, IND

    !=======================================================================
    ! count number of plane wave coefficients
    ! and allocate required arrays
    !=======================================================================
    NRPLWV=0

    DO NK=1,WDES%NKPTS
       IND=COUNT_INDPW(GRID,WDES%ENMAX,B,BI,WDES%VKPT(:,NK))
       NRPLWV=MAX(NRPLWV,IND)
    ENDDO

    CALL COUNT_INDPW_FULL(GRID,WDES,BI,NRPLWV)
#ifdef MPI
    ! ok this is tricky
    ! if symmetry operations are applied, we might have too  few plane wave coefficients
    ! in the parallel mode at k-points outside the IRZ, whenever HF or symmetry routines are used
    IF (WDES%NKDIM /= WDES%NKPTS .AND. GRID%COMM%NCPU/=1 ) THEN
!       NRPLWV=NRPLWV+8*GRID%COMM%NCPU  ! safeguard this certainly needs to be improved
! gK symmetry can usually not yield more PW coefficients than on any other core
! so max it
       CALLMPI( M_max_i(GRID%COMM, NRPLWV, 1 ))
       NRPLWV=NRPLWV+8
    ENDIF
#endif
    ! make WDES%NRPLWV dividable by NB_PAR
    NRPLWV=((NRPLWV+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
  END FUNCTION WDES_COUNT_INDPW

!************************ FUNCTION COUNT_INDPW *************************
!> Count the number of plane-waves for a given k-point GRID and ENMAX
!***********************************************************************
  FUNCTION COUNT_INDPW(GRID,ENMAX,B,BI,VKPT) RESULT(NRPLWV)
    use constant, only: HSQDTM, TPI
    IMPLICIT NONE
    TYPE (grid_3d), TARGET :: GRID
    REAL(q) :: B(3,3), BI(3,3)
    REAL(q) :: ENMAX
    REAL(q) :: VKPT(3)
    INTEGER :: NRPLWV
    ! local variables
    REAL(q) :: ENERGI
    REAL(q) :: GIX,GIY,GIZ,G1,G2,G3
    INTEGER :: NK,NC,N1,N2,N3

    !=======================================================================
    ! count number of plane wave coefficients
    ! and allocate required arrays
    !=======================================================================
    NRPLWV=0
    DO NC=1,GRID%RC%NCOL
       N2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(N2)+VKPT(2))
       N3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(N3)+VKPT(3))
       DO N1=1,GRID%RC%NROW

          G1=(GRID%LPCTX(N1)+VKPT(1))

          GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
          GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
          GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

          ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

          ! exclude some components for gamma-only version (C(G)=C*(-G))
          IF (GRID%NGX/=GRID%NGX_rd) THEN
             IF (GRID%LPCTX(N1)==0 .AND. GRID%LPCTY(N2)<0) CYCLE
             IF (GRID%LPCTX(N1)==0 .AND. GRID%LPCTY(N2)==0 .AND. GRID%LPCTZ(N3)<0) CYCLE
          ENDIF
          IF (GRID%NGZ/=GRID%NGZ_rd) THEN
             IF (GRID%LPCTZ(N3)==0 .AND. GRID%LPCTY(N2)<0) CYCLE
             IF (GRID%LPCTZ(N3)==0 .AND. GRID%LPCTY(N2)==0 .AND. GRID%LPCTX(N1)<0) CYCLE
          ENDIF
          IF(ENERGI<ENMAX) THEN
             NRPLWV=NRPLWV+1
          ENDIF
       ENDDO
    ENDDO
  END FUNCTION COUNT_INDPW

!*********************SUBROUTINE WDES_SET_NRPLWV ***********************
!> Set NRPLWV in the wavefunction descriptor
!***********************************************************************
  SUBROUTINE WDES_SET_NRPLWV(WDES,NRPLWV)
    TYPE (wavedes), INTENT(INOUT) :: WDES
    INTEGER,INTENT(IN) :: NRPLWV

    ! make WDES%NRPLWV dividable by NB_PAR
    WDES%NRPLWV=((NRPLWV+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    WDES%NGDIM=WDES%NRPLWV
    WDES%NRPLWV = WDES%NRPLWV*WDES%NRSPINORS
    WDES%NRPLWV_RED=WDES%NRPLWV/WDES%NB_PAR
  END SUBROUTINE WDES_SET_NRPLWV

!*************************SUBROUTINE GEN_LAYOUT_NEW *********************
!
!> subroutine GEN_LAYOUT_NEW performs a number of tasks:
!> - determines the layout (distribution) of the columns on parallel
!>   computers
!>   essentially these are the structures GRID%RC, GRID%RL
!>   and in the parallel version the array GRID%IN
!>
!> - also the following  arrays are allocated:
!>     WDES%NPLWKP  WDES%NINDPW
!> - for LSETUP=.TRUE. the kinetic energy arrays and the G-vector array
!>   are additionally allocated
!>     WDES%DATAKE, WDES%IGX, Y, Z
!>  setup of these arrays is performed in the routines GEN_INDEX
!>
!> the data layout is based on the initial reciprocal lattice vectors
!> stored in BI
!> compared to GEN_LAYOUT this version can use inversion symmetry
!> to reduce the number of G vectors in the compressed orbital array
!> in this case the caller must set LUSEINV
!
!***********************************************************************

  SUBROUTINE GEN_LAYOUT_NEW(GRID,WDES, B,BI,IU6,LSETUP, LUSEINV_IN, LNOGAMMA,LREAL)
    USE prec
    USE mgrid
    USE constant
    USE base
    IMPLICIT NONE

    TYPE (grid_3d), TARGET :: GRID
    TYPE (wavedes)     WDES
    !> current lattice
    REAL(q) B(3, 3)
    !> initial reciprocal lattice
    REAL(q) BI(3, 3)
    !> OUTCAR file handle
    INTEGER           :: IU6
    !> if .TRUE. WDES is ALLOCATED and some entries are set (see above)
    LOGICAL:: LSETUP
    !> use inversion symmetry C_G+k = C*_-G-k
    LOGICAL LUSEINV_IN
    !> force complex to complex FFT
    LOGICAL, OPTIONAL:: LNOGAMMA
    !> force real storage mode in real space if complex to real FFT used
    LOGICAL, OPTIONAL:: LREAL
    ! local
    LOGICAL LUP

    LOGICAL, ALLOCATABLE :: LUSE_IN(:)
    INTEGER, ALLOCATABLE :: USED_ROWS(:),IND2(:),IND3(:)
    INTEGER, ALLOCATABLE :: REDISTRIBUTION_INDEX(:)
    INTEGER :: IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX
    COMMON /WAVCUT/   IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX
    REAL(q), PARAMETER :: G2ZERO=1E-12_q
    INTEGER :: NK, NCOL_TOT, N1, N2, N3, IUSED, IUSED_MAX, I, NCOL_MAX, IND
    INTEGER :: NODE_TARGET, NCOL, IND_REDIS, NC, NRPLWV, IND_ON_CPU, N
    REAL(q) :: G1, G2, G3, GIX, GIY, GIZ, ENERGI
    LOGICAL :: LGX_REDUCED, LGZ_REDUCED
    LOGICAL :: LUSEINV(WDES%NKDIM)

#ifdef gammareal
    WDES%LGAMMA=.TRUE.
    GRID%REAL2CPLX=.TRUE.
#else
    WDES%LGAMMA=.FALSE.
    GRID%REAL2CPLX=.FALSE.
#endif
    ! force complex to complex FFT
    IF (PRESENT(LNOGAMMA)) THEN
       IF ( LNOGAMMA ) THEN
          WDES%LGAMMA=.FALSE.
          GRID%REAL2CPLX=.FALSE.
    ! force complex to real FFT and force real storage mode in real space
       ELSE
          WDES%LGAMMA=.TRUE.
          GRID%REAL2CPLX=.TRUE.
       ENDIF
    ENDIF
    ! orbitals are usually complex in the direct grid in VASP hence LREAL defaults to .FALSE.
    GRID%LREAL=.FALSE.
    ! if LREAL is present and we are at the Gamma point use that as default
    IF (PRESENT(LREAL) .AND. WDES%LGAMMA) THEN
       GRID%LREAL=LREAL
    ENDIF
    WDES%GRID=> GRID

    GRID%NGX_rd=GRID%NGX
    GRID%NGY_rd=GRID%NGY
    GRID%NGZ_rd=GRID%NGZ

    CALL NULLIFY_GRD(GRID)
!gK LUSEINV
! LUSEINV determines whether a specific k-point corresponds to
!  (+-1/2 | 0, +-1/2 | 0, +-1/2 | 0)
! if this is the case, only half the PW coefficients are stored
    DO NK=1,WDES%NKDIM
       IF (ABS(MOD(WDES%VKPT(1,NK)*2+100,1._q))> SQRT(G2ZERO) .OR. &
           ABS(MOD(WDES%VKPT(2,NK)*2+100,1._q))> SQRT(G2ZERO) .OR. &
           ABS(MOD(WDES%VKPT(3,NK)*2+100,1._q))> SQRT(G2ZERO)) THEN
          LUSEINV(NK)=.FALSE.
       ELSE
          LUSEINV(NK)=LUSEINV_IN
       ENDIF
    ENDDO
!gK LUSEINV

    WDES%NCOL=0
!=======================================================================
! determine the layout
! i.e. all required columns
! or (y,z) pairs which are required for the 3d-FFT grid
!=======================================================================
#ifdef MPI
    IF (GRID%COMM%NCPU/=1 .OR. LUSE_PARALLEL_FFT) THEN
       IF (WDES%LGAMMA) THEN
          GRID%NGZ_rd=GRID%NGZ/2+1
          IF (IU6>=0) THEN
             WRITE(IU6,*) 'use parallel FFT for orbitals z direction half grid'
          ENDIF
       ENDIF
!-----------------------------------------------------------------------
! parallel version
! reciprocal space first
!-----------------------------------------------------------------------

       ! determine all required grid points (usually a circle around center)

       ALLOCATE(USED_ROWS(GRID%NGY*GRID%NGZ_rd), &
            IND2(GRID%NGY*GRID%NGZ_rd),IND3(GRID%NGY*GRID%NGZ_rd), &
            REDISTRIBUTION_INDEX(GRID%NGY*GRID%NGZ_rd), &
            LUSE_IN(GRID%NGZ_rd))
       ! following arrays are required to setup the reciprocal layout
       USED_ROWS=0     ! gives the total number of points in each column
       ! contained within the cutoff sphere
       NCOL_TOT=0      ! total number of columns
       IND2=0          ! y index of each used column
       IND3=0          ! z index of each used column
       ! LUSE_IN array is used to find out which z-planes must be contained
       ! in the intermediate representation
       LUSE_IN=.FALSE.

       DO N3=1,GRID%NGZ_rd
          DO N2=1,GRID%NGY

             IUSED_MAX=0   ! maxmimal number of rows for this column

             ! loop over all k-points in the *entire* or irreducible Brillouin zone (N)
             ! NKPTS_FOR_GEN_LAYOUT is set, if HF is used, or if the number of k-points
             ! is allowed to change dynamically (linear response, finite differences)
             DO NK=1,WDES%NKPTS_FOR_GEN_LAYOUT
                G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
                G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
                IUSED=0     ! number of rows for this column and k-point
                DO N1=1,GRID%NGX_rd
                   G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
                   GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
                   GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
                   GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI
                   ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

                   ! exclude some components for gamma-only version (C(G)=C*(-G))
                   IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                      IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                      IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
                   ENDIF

                   IF(ENERGI<WDES%ENMAX) THEN
                      IUSED=IUSED+1
                      LUSE_IN(N3)=.TRUE.
                   ENDIF
                ENDDO
                IUSED_MAX=MAX(IUSED,IUSED_MAX)
             ENDDO
             ! if this column is used, enter it in required array
             IF (IUSED_MAX>0) THEN
                NCOL_TOT=NCOL_TOT+1        ! increase column counter
                USED_ROWS(NCOL_TOT)=IUSED  ! set number of rows in this column
                IND2(NCOL_TOT)=N2          ! set y index of this column
                IND3(NCOL_TOT)=N3          ! set z index of this column
             ENDIF
          ENDDO
       ENDDO
       ! setup redistribution array
       ! by sorting according to the number of rows in each column
       ! we try to get better load balancing
       DO I=1,NCOL_TOT
          REDISTRIBUTION_INDEX(I)=I
       ENDDO
       ! do not sort column 1 (it must store the G vector G=(0,0,0)) !
       CALL SORT_REDIS(NCOL_TOT,USED_ROWS(1),REDISTRIBUTION_INDEX(1))

       ! distribute among processors

       NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU
       ALLOCATE(GRID%RC%I2( NCOL_MAX ))
       ALLOCATE(GRID%RC%I3( NCOL_MAX ))

       GRID%RC%NFAST= 1
       GRID%RC%NROW = GRID%NGX_rd
       GRID%RC%NCOL = 0

       IND=0
       !
       ! step through all columns and distribute them onto proc.
       ! in the manner 1 ... NCPU - NCPU ... 1 - 1 ... NCPU - etc.
       !
       NODE_TARGET=0  ! NODE onto which column has to go
       LUP=.TRUE.     ! determines whether NODE_TARGET is increased or decreased

       DO NCOL=1,NCOL_TOT
          IND_REDIS=REDISTRIBUTION_INDEX(NCOL)
          N2=IND2(IND_REDIS)
          N3=IND3(IND_REDIS)
          IF (LUP) THEN
             IF (NODE_TARGET == GRID%COMM%NCPU) THEN
                LUP=.FALSE.
             ELSE
                NODE_TARGET=NODE_TARGET+1
             ENDIF
          ELSE
             IF (NODE_TARGET == 1) THEN
                LUP=.TRUE.
             ELSE
                NODE_TARGET=NODE_TARGET-1
             ENDIF
          ENDIF

          IND_ON_CPU=(NCOL-1)/GRID%COMM%NCPU+1

          IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
             GRID%RC%NCOL=GRID%RC%NCOL+1
             IF (IND_ON_CPU /= GRID%RC%NCOL) THEN
                CALL vtutor%bug("GENSP: internal error(1) " // str(GRID%COMM%NODE_ME) // " " // &
                   str(IND_ON_CPU) // " " // str(GRID%RC%NCOL), __FILE__, __LINE__)
             ENDIF
             GRID%RC%I2(IND_ON_CPU)=N2
             GRID%RC%I3(IND_ON_CPU)=N3
          ENDIF
       ENDDO
       GRID%RC%NP    = GRID%RC%NCOL*GRID%RC%NROW
       GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
       WDES%NCOL     = GRID%RC%NCOL
       !
       ! set up intermediate grid
       ! this usually contains a certain number of z-planes around the center
       !
       LUSE_IN=.TRUE.

       NCOL_TOT=0
       DO N3=1,GRID%NGX_rd
          DO N2=1,GRID%NGZ_rd
             IF (LUSE_IN(N2)) NCOL_TOT=NCOL_TOT+1
          ENDDO
       ENDDO

       NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU

       IF (LPLANE_WISE) THEN
          NCOL_MAX=GRID%NGZ*((GRID%NGX+GRID%COMM%NCPU-1)/GRID%COMM%NCPU)
       ENDIF

       ALLOCATE(GRID%IN%I2( NCOL_MAX ))
       ALLOCATE(GRID%IN%I3( NCOL_MAX ))

       GRID%IN%NFAST= 2
       GRID%IN%NCOL = 0
       GRID%IN%NROW = GRID%NGY

       IND=0
       DO N3=1,GRID%NGX_rd
          DO N2=1,GRID%NGZ_rd
             IF (LUSE_IN(N2)) THEN
                IF (LPLANE_WISE) THEN
                   NODE_TARGET=MOD(N3-1,GRID%COMM%NCPU)+1
                ELSE
                   NODE_TARGET=MOD(IND,GRID%COMM%NCPU)+1
                ENDIF
                IND=IND+1
                IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
                   N=GRID%IN%NCOL
                   N=N+1
                   GRID%IN%NCOL=N
                   GRID%IN%I2(N)=N2 ! I2 contains z index
                   GRID%IN%I3(N)=N3 ! I3      the x index
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF (GRID%IN%NCOL >NCOL_MAX) THEN
          CALL vtutor%bug("GEN_LAYOUT: internal error " // str(NCOL_MAX) // " " // str(GRID%IN%NCOL), __FILE__, __LINE__)
       ENDIF
       GRID%IN%NP    = GRID%IN%NCOL*GRID%IN%NROW
       GRID%IN%NALLOC= GRID%IN%NCOL*GRID%IN%NROW
       !
       ! set up real space grid, use all points
       ! in this case GRID will be conformable to GRID_SOFT in real space
       !
       CALL REAL_STDLAY(GRID, LPARALLEL=.TRUE.)

       DEALLOCATE(USED_ROWS,IND2,IND3,REDISTRIBUTION_INDEX,LUSE_IN)
    ELSE
#endif
!-----------------------------------------------------------------------
! serial version
! set reciprocal space and real space layout for non parallel computers
! always x-first (or x-fast) layout
! all grid points are used for FFT
!-----------------------------------------------------------------------
       IF (WDES%LGAMMA) THEN
          GRID%NGX_rd=GRID%NGX/2+1
          IF (IU6>=0) WRITE(IU6,*) 'use serial FFT for orbitals x direction half grid'
       ENDIF

       GRID%RC%NFAST= 1
       GRID%RC%NCOL = GRID%NGZ_rd*GRID%NGY
       WDES%NCOL    = GRID%RC%NCOL
       GRID%RC%NROW = GRID%NGX_rd
       GRID%RC%NP   = GRID%RC%NCOL*GRID%RC%NROW
       GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
       ALLOCATE(GRID%RC%I2( GRID%RC%NCOL ))
       ALLOCATE(GRID%RC%I3( GRID%RC%NCOL ))
       IND=1
       DO N3=1,GRID%NGZ_rd
          DO N2=1,GRID%NGY
             GRID%RC%I2(IND)=N2
             GRID%RC%I3(IND)=N3
             IND=IND+1
          ENDDO
       ENDDO

       CALL REAL_STDLAY(GRID, LPARALLEL=.FALSE.)
#ifdef MPI
    ENDIF
#endif
!=======================================================================
! count number of plane wave coefficients
! and allocate required arrays
!=======================================================================
    NRPLWV=0

    DO NK=1,WDES%NKPTS
!gK LUSEINV
       IF (LUSEINV(NK)) THEN
          LGX_REDUCED=.TRUE.
          LGZ_REDUCED=.FALSE.
       ELSE
          LGX_REDUCED=.FALSE.
          LGZ_REDUCED=.FALSE.
       ENDIF
!gK LUSEINV
       IND=0
       DO NC=1,GRID%RC%NCOL
          N2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
          N3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
          DO N1=1,GRID%RC%NROW

             G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))

             GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
             GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
             GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

             ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

             ! exclude some components for gamma-only version (C(G)=C*(-G))
             IF (GRID%NGX/=GRID%NGX_rd) THEN
                IF (GRID%LPCTX(N1)==0 .AND. GRID%LPCTY(N2)<0) CYCLE
                IF (GRID%LPCTX(N1)==0 .AND. GRID%LPCTY(N2)==0 .AND. GRID%LPCTZ(N3)<0) CYCLE
             ELSE IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                IF (GRID%LPCTZ(N3)==0 .AND. GRID%LPCTY(N2)<0) CYCLE
                IF (GRID%LPCTZ(N3)==0 .AND. GRID%LPCTY(N2)==0 .AND. GRID%LPCTX(N1)<0) CYCLE
!gK LUSEINV
             ELSE IF (LGX_REDUCED) THEN
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)<0) CYCLE
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)<0) CYCLE
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)==0 .AND. GRID%LPCTZ(N3)+WDES%VKPT(3,NK)<0) CYCLE
             ELSE IF (LGZ_REDUCED) THEN
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)<0) CYCLE
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)<0) CYCLE
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)==0 .AND. GRID%LPCTX(N1)+WDES%VKPT(1,NK)<0) CYCLE
             ENDIF
!gK LUSEINV
             IF(ENERGI<WDES%ENMAX) THEN
                IND=IND+1
             ENDIF
          ENDDO
       ENDDO
       NRPLWV=MAX(NRPLWV,IND)
    ENDDO
    ! make WDES%NRPLWV dividable by NB_PAR

#ifdef MPI
    ! ok this is tricky
    ! if symmetry operations are applied, we might have too  few plane wave coefficients
    ! in the parallel mode at k-points outside the IRZ, whenever HF or symmetry routines are used
    IF (WDES%NKDIM /= WDES%NKPTS .AND. GRID%COMM%NCPU/=1 ) THEN
!       NRPLWV=NRPLWV+8*GRID%COMM%NCPU  ! safeguard this certainly needs to be improved
! gK symmetry can usually not yield more PW coefficients than on any other core
! so max it
       CALLMPI( M_max_i(GRID%COMM, NRPLWV, 1 ))
       NRPLWV=NRPLWV+8
    ENDIF
#endif
    WDES%NRPLWV=((NRPLWV+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    WDES%NGDIM=WDES%NRPLWV
    WDES%NRPLWV = WDES%NRPLWV*WDES%NRSPINORS
    WDES%NRPLWV_RED=WDES%NRPLWV/WDES%NB_PAR

    WDES%DO_REDIS = .FALSE.
    IF (WDES%NB_PAR /=1) THEN
       WDES%DO_REDIS=.TRUE.
    ENDIF

    !    CALL  MAKE_STRIDE(WDES%NRPLWV)

    GRID%MPLWV=__align_DIM(MAX(GRID%RC%NALLOC ,GRID%IN%NALLOC , GRID%RL%NALLOC))

    CALL ALLOCWDES(WDES,LSETUP)
    WDES%LUSEINV=LUSEINV
    IF (LUSEINV_IN) THEN
       ALLOCATE( WDES%NINDPW_INV(WDES%NGDIM, WDES%NKDIM), WDES%FFTSCA(WDES%NGDIM, 2, WDES%NKDIM), WDES%MAP_TO_FULL(WDES%NGDIM, 2, WDES%NKDIM))
    ENDIF

    WDES%NPLWKP=0
    WDES%NGVECTOR=0
    WDES%NGVECTOR_POS=0

    IF      (WDES%ISPIN==1  .AND. .NOT. WDES%LNONCOLLINEAR ) THEN
       WDES%NCDIJ=1
    ELSE IF (WDES%ISPIN==2) THEN
       WDES%NCDIJ=2
    ELSE IF (WDES%ISPIN==1  .AND. WDES%LNONCOLLINEAR ) THEN
       WDES%NCDIJ=4
    ELSE
       CALL vtutor%bug("internal error: can not set NCDIJ", __FILE__, __LINE__)
    ENDIF

! set WDES%AT_GAMMA
! this flag determines whether a specific k-point corresponds to Gamma
! if this is the case the code might do some special treatment allowed
! only at Gamma (e.g. orbitals are real)
! presently this is only supported in few places
    DO NK=1,WDES%NKDIM
       IF (ABS(SUM(WDES%VKPT(:,NK)**2)) < G2ZERO .AND. .NOT. WDES%LNONCOLLINEAR) THEN
          WDES%AT_GAMMA(NK)=.TRUE.
       ELSE
          WDES%AT_GAMMA(NK)=.FALSE.
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE GEN_LAYOUT_NEW


!*************************SUBROUTINE GEN_INDEX_GAMMA ******************
!
!> the following subroutine creates an index array
!> to store orbital coefficients strictly compatible to
!> the serial Gamma only version compiled using-DWNGXhalf
!> this requires two arrays
!>
!>  -  one telling where to store a coefficient
!>  -  and one telling whether the coefficient needs to be
!>     complex conjugated
!
!***********************************************************************


  SUBROUTINE GEN_INDEX_GAMMA(WDES1, B, ENMAX, INDEX, LCONJG )
    USE constant
    IMPLICIT NONE
    TYPE (wavedes1)    WDES1
    !> reciprocal lattice
    REAL(q) :: B(3,3)
    REAL(q) :: ENMAX
    LOGICAL, POINTER :: LCONJG(:)
    INTEGER, POINTER :: INDEX(:)
  ! local
    INTEGER :: M, N1, N2, N3, NUSED, IND
    INTEGER, ALLOCATABLE  :: INDEX3D(:,:,:)
    REAL(q) :: G1, G2, G3, GIX, GIY, GIZ, ENERGI

    IF (ASSOCIATED(INDEX) .OR. ASSOCIATED(LCONJG)) THEN
       CALL vtutor%bug("internal error in GEN_INDEX_GAMMA: not null pointer passed down", __FILE__, __LINE__)
    ENDIF

    ALLOCATE( INDEX(WDES1%NGVECTOR), LCONJG(WDES1%NGVECTOR))


    ALLOCATE( INDEX3D(0:WDES1%GRID%NGX/2, -WDES1%GRID%NGY/2:WDES1%GRID%NGY/2 , -WDES1%GRID%NGZ/2:WDES1%GRID%NGZ/2))
  ! loop over all grid points
  ! assume reduction in x direction
    INDEX3D=0
    NUSED=0
    DO N3=1,WDES1%GRID%NGZ
       DO N2=1,WDES1%GRID%NGY
          G3=(WDES1%GRID%LPCTZ(N3)+WDES1%VKPT(3))
          G2=(WDES1%GRID%LPCTY(N2)+WDES1%VKPT(2))

          row: DO N1=1,WDES1%GRID%NGX/2+1

             G1=(WDES1%GRID%LPCTX(N1)+WDES1%VKPT(1))
             GIX= (G1*B(1,1)+G2*B(1,2)+G3*B(1,3)) *TPI
             GIY= (G1*B(2,1)+G2*B(2,2)+G3*B(2,3)) *TPI
             GIZ= (G1*B(3,1)+G2*B(3,2)+G3*B(3,3)) *TPI

             ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
             ! exclude some components for gamma-only version (C(G)=C*(-G))

             IF (WDES1%GRID%LPCTX(N1)==0 .AND.WDES1%GRID%LPCTY(N2)<0) CYCLE
             IF (WDES1%GRID%LPCTX(N1)==0 .AND.WDES1%GRID%LPCTY(N2)==0 .AND.WDES1%GRID%LPCTZ(N3)<0) CYCLE

             IF(ENERGI<ENMAX) THEN
                NUSED=NUSED+1
                INDEX3D(WDES1%GRID%LPCTX(N1) ,WDES1%GRID%LPCTY(N2), WDES1%GRID%LPCTZ(N3))=NUSED
             ENDIF
          ENDDO row
       ENDDO
    ENDDO

    DO IND=1,WDES1%NGVECTOR
       N1=WDES1%IGX(IND)
       N2=WDES1%IGY(IND)
       N3=WDES1%IGZ(IND)
       IF (N1<0) THEN
          LCONJG(IND)=.TRUE.
          INDEX (IND)= INDEX3D(-N1 ,-N2, -N3)
       ELSE IF (N1==0) THEN
          ! try first conventional one
          INDEX (IND)= INDEX3D(N1 ,N2, N3)
          IF (INDEX (IND) >0 ) THEN
             ! ok no conjugation
             LCONJG(IND)=.FALSE.
          ELSE
             ! try inverted one
             INDEX (IND)= INDEX3D(-N1 ,-N2, -N3)
             LCONJG(IND)=.TRUE.
          ENDIF
       ELSE
          ! must be conventional one
          LCONJG(IND)=.FALSE.
          INDEX (IND)= INDEX3D(N1 ,N2, N3)
       ENDIF
       IF (INDEX(IND)==0) THEN
          CALL vtutor%bug("internal error in GEN_INDEX_GAMMA: can not find appropriate entry in &
             &INDEX3d " // str(N1) // " " // str(N2) // " " // str(N3) // " " // str(INDEX3D(N1,N2,&
             N3)) // " " // str(INDEX3D(-N1,-N2,-N3)), __FILE__, __LINE__)
       ENDIF
    ENDDO

    DEALLOCATE(INDEX3D)
  END SUBROUTINE GEN_INDEX_GAMMA


  SUBROUTINE FREE_INDEX_GAMMA(INDEX, LCONJG)
    LOGICAL, POINTER :: LCONJG(:)
    INTEGER, POINTER :: INDEX(:)

    DEALLOCATE(INDEX)
    DEALLOCATE(LCONJG)

    NULLIFY(INDEX)
    NULLIFY(LCONJG)
  END SUBROUTINE FREE_INDEX_GAMMA

  !> Get default number of bands to accomodate all the electrons
  !> This is the default number of bands that one would get when
  !> running the serial version of VASP, i.e. the number of bands
  !> normally is still increased such that the number of bands is divisible
  !> by NPAR. This had to be done after calling this routine
  !> This code has been moved from main.F
  FUNCTION GET_NBANDS_SERIAL(LNONCOLLINEAR,ISPIN,NELECT,T_INFO) RESULT(NBANDS)
   USE poscar, ONLY: type_info
   IMPLICIT NONE
   LOGICAL,INTENT(IN) :: LNONCOLLINEAR
   INTEGER,INTENT(IN) :: ISPIN
   REAL(q),INTENT(IN) :: NELECT
   TYPE(TYPE_INFO),INTENT(IN) :: T_INFO
   INTEGER :: NBANDS
   ! local variables
   INTEGER :: NMAG
   IF (LNONCOLLINEAR)  THEN
       NMAG=MAX(ABS(SUM(T_INFO%ATOMOM(1:T_INFO%NIONS*3-2:3))), &
                ABS(SUM(T_INFO%ATOMOM(2:T_INFO%NIONS*3-1:3))), &
                ABS(SUM(T_INFO%ATOMOM(3:T_INFO%NIONS*3:3))))
   ELSE IF (ISPIN > 1) THEN
       NMAG=ABS(SUM(T_INFO%ATOMOM(1:T_INFO%NIONS)))
   ELSE
       NMAG=0
   ENDIF
   NMAG = (NMAG+1)/2
   NBANDS=MAX(NINT(NELECT+2)/2+MAX(T_INFO%NIONS/2,3),INT(0.6*NELECT))+NMAG
   IF (LNONCOLLINEAR) NBANDS = NBANDS*2
  END FUNCTION GET_NBANDS_SERIAL
END MODULE


!=======================================================================
!> sorts RA in descending order, and rearanges an index array RB
!> seems to be a quicksort, but I am not sure
!> subroutine writen by Florian Kirchhof
!=======================================================================

  SUBROUTINE SORT_REDIS(N,RA,RB)
    INTEGER RA(N),RB(N)
    INTEGER RRA,RRB

    IF (N<=1) RETURN

    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
       RRB=RB(L)
    ELSE
       RRA=RA(IR)
       RRB=RB(IR)
       RA(IR)=RA(1)
       RB(IR)=RB(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RB(1)=RRB
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).GT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.GT.RA(J))THEN
          RA(I)=RA(J)
          RB(I)=RB(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  END SUBROUTINE SORT_REDIS

!=======================================================================
!> sorts RA in ascending order, and rearanges an index array RB
!> seems to be a quicksort, by I am not sure
!> subroutine writen by Florian Kirchhof
!=======================================================================

  SUBROUTINE SORT_REDIS_ASC(N,RA,RB)
    INTEGER RA(N),RB(N)
    INTEGER RRA,RRB

    IF (N<=1) RETURN

    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
       RRB=RB(L)
    ELSE
       RRA=RA(IR)
       RRB=RB(IR)
       RA(IR)=RA(1)
       RB(IR)=RB(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RB(1)=RRB
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          RB(I)=RB(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  END SUBROUTINE SORT_REDIS_ASC

!*************************SUBROUTINE COUNT_ROWS ************************
!
!> this subroutine counts the total number of plane waves contained
!> within the cutoff sphere up to (but excluding) a certain column
!> this array is required to find out which index a certain column would
!> have in the serial version
!>
!> mind the index is 0 based
!> the total number of plane waves is returned in NUSED
!> USED_POINTS returns the number of plane wave coefficients
!> up to column N1, N3
!
!***********************************************************************

  SUBROUTINE  COUNT_ROWS(GRID, WDES, BI, NK, USED_POINTS, LGX_REDUCED, LGZ_REDUCED, NUSED )
    USE mgrid
    USE wave
    USE constant
    USE base
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d)     GRID
    TYPE (wavedes)     WDES
    DIMENSION BI(3,3)
    LOGICAL :: LGX_REDUCED, LGZ_REDUCED

    INTEGER :: USED_POINTS(GRID%NGY,GRID%NGZ)

    NUSED=0
    DO N3=1,GRID%NGZ_rd
       DO N2=1,GRID%NGY
          G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
          G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
          USED_POINTS(N2,N3)=NUSED

          row: DO N1=1,GRID%NGX_rd

             G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
             GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
             GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
             GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

             ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
             ! exclude some components for gamma-only version (C(G)=C*(-G))
             IF (GRID%NGX/=GRID%NGX_rd) THEN
                IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE
             ELSE IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
             ELSE IF (LGX_REDUCED) THEN
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)<0) CYCLE
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)<0) CYCLE
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)==0 .AND. GRID%LPCTZ(N3)+WDES%VKPT(3,NK)<0) CYCLE
             ELSE IF (LGZ_REDUCED) THEN
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)<0) CYCLE
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)<0) CYCLE
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)==0 .AND. GRID%LPCTX(N1)+WDES%VKPT(1,NK)<0) CYCLE
             ENDIF
             IF(ENERGI<WDES%ENMAX) THEN
                NUSED=NUSED+1
             ENDIF
          ENDDO row
       ENDDO
    ENDDO
  END SUBROUTINE COUNT_ROWS

!*************************SUBROUTINE REPAD_INDEX_ARRAY  ****************
!
!> this subroutine calculates two index arrays that allow
!> to "restore" a plane wave array from an old cutoff and lattice
!> to a new one
!>
!> this operation works only if the orbitals are stored in
!> the serial layout on one single core (not parallel)
!> this is usually the case when orbitals are read in
!> in the Gamma point only version the orbitals are presently
!> stored in a way that is always compatible to the seriel
!> version (implying that the grid is reduced in the x direction
!
!***********************************************************************

  SUBROUTINE  REPAD_INDEX_ARRAY(GRID, VKPT, VKPTI, B,  BI, ENMAX, ENMAXI, &
       NP, NPI, IND, INDI, INDMAX, IFAIL )
    USE prec
    USE mgrid
    USE constant
    USE base
    IMPLICIT NONE

    !> grid descriptor
    TYPE (grid_3d)     GRID
    !> new k-point
    REAL(q):: VKPT(3)
    !> old k-point
    REAL(q):: VKPTI(3)
    !> new reciprocal lattice constant
    REAL(q) ::  B (3, 3)
    !> old reciprocal lattice constant
    REAL(q) ::  BI(3, 3)
    !> new cutoff
    REAL(q) ::  ENMAX
    !> old cutoff
    REAL(q) ::  ENMAXI
    !> number of plane wave coefficients old, must be set by caller
    INTEGER ::  NP
    !> number of plane wave coefficients new, must be set by caller
    INTEGER ::  NPI
    !> index array new
    INTEGER ::  IND(MAX(NP, NPI))
    !> index array old
    INTEGER ::  INDI(MAX(NP, NPI))
    !> on return maximum index
    INTEGER ::  INDMAX
    !> 0  NP and NPI were correct
    !> 1  NP was incorrect, 2 NPI was incorrect
    INTEGER ::  IFAIL
    ! local
    INTEGER NP_,NPI_,N1,N2,N3
    INTEGER NGZ_rd, NGX_rd
    REAL(q) :: G1,G2,G3, G1I,G2I,G3I, GIX,GIY,GIZ, GX,GY,GZ, ENERGI, ENERG

    IFAIL = 0

    NP_ =0
    NPI_=0
    INDMAX=0

    NGX_rd =  GRID%NGX
    NGZ_rd =  GRID%NGZ
    ! grid reduced: the orbitals are stored compatible to seriel version
    ! hence the reduction is always along x direction
    IF (GRID%NGZ_rd /= GRID%NGZ .OR.  GRID%NGX_rd /=  GRID%NGX) THEN
       NGX_rd =  GRID%NGX/2+1
       NGZ_rd =  GRID%NGZ
    ENDIF

    DO N3=1,NGZ_rd
       DO N2=1,GRID%NGY
          G3=(GRID%LPCTZ(N3)+VKPT(3))
          G2=(GRID%LPCTY(N2)+VKPT(2))
          G3I=(GRID%LPCTZ(N3)+VKPTI(3))
          G2I=(GRID%LPCTY(N2)+VKPTI(2))

          row: DO N1=1,NGX_rd

             G1=(GRID%LPCTX(N1)+VKPT(1))
             G1I=(GRID%LPCTX(N1)+VKPTI(1))

             GIX= (G1I*BI(1,1)+G2I*BI(1,2)+G3I*BI(1,3)) *TPI
             GIY= (G1I*BI(2,1)+G2I*BI(2,2)+G3I*BI(2,3)) *TPI
             GIZ= (G1I*BI(3,1)+G2I*BI(3,2)+G3I*BI(3,3)) *TPI

             GX= (G1*B(1,1)+G2*B(1,2)+G3*B(1,3)) *TPI
             GY= (G1*B(2,1)+G2*B(2,2)+G3*B(2,3)) *TPI
             GZ= (G1*B(3,1)+G2*B(3,2)+G3*B(3,3)) *TPI

             ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
             ENERG =HSQDTM*( (GX**2)+ (GY**2)+ (GZ**2))

             ! exclude some components for gamma-only version (C(G)=C*(-G))
             IF (GRID%NGX/=NGX_rd) THEN
                IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE
             ENDIF
             IF (GRID%NGZ/=NGZ_rd) THEN
                IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
             ENDIF
             IF (ENERG <ENMAX) THEN
                NP_=NP_+1    ! increase index for new array
             ENDIF
             IF (ENERGI<ENMAXI) THEN
                NPI_=NPI_+1  ! increase index for old array
             ENDIF
             IF (ENERG<ENMAX .AND. ENERGI<ENMAXI) THEN
                INDMAX= MIN( MIN(INDMAX+1 ,NP ), NPI)
                ! increase index, and avoid overrun
                IND(INDMAX) =NP_
                INDI(INDMAX)=NPI_
             ENDIF
          ENDDO row
       ENDDO
    ENDDO

    IF  (NP_ /= NP) THEN
       NP=NP_
       IFAIL=1
    ENDIF
    IF  (NPI_ /= NPI) THEN
       NPI=NPI_
       IFAIL=1
    ENDIF

  END SUBROUTINE REPAD_INDEX_ARRAY

  SUBROUTINE REPAD_WITH_INDEX_ARRAY(INDMAX,IND,INDI, CW, CWI)
    USE prec
    IMPLICIT NONE

    INTEGER INDMAX
    INTEGER ::  IND(INDMAX)  ! index array new
    INTEGER ::  INDI(INDMAX) ! index array old
    COMPLEX(q) :: CW(*),CWI(*)
    ! local
    INTEGER I

    DO I=1,INDMAX
       CW(IND(I))=CWI(INDI(I))
    ENDDO

  END SUBROUTINE REPAD_WITH_INDEX_ARRAY

!*************************SUBROUTINE GEN_INDEX ************************
!
!> subroutine GEN_INDEX calculates the following arrays:
!> - the indexing array NINDPW for copying the plane wave coefficients
!>   from the continuous array CW to the column wise layout used for
!>   the 3d-FFT
!> for LSETUP=.TRUE., additionally the following array are set up:
!> - the kinetic energies of the plane wave basis states are computed
!>   (WDES%DATAKE)
!> - the G vector corresponding to each plane wave basis state is stored
!>   (WDES%IGX, Y, Z)
!>
!> - in the parallel version, the arrays PL_INDEX and PL_COL
!>   are set up and stored
!>     PL_INDEX(NC, NK) stores the position of a column at
!>                     which data is stored in the serial version
!>     PL_COL(NC, NK)   number of data in this column
!>
!> the data layout is based on the initial reciprocal lattice vectors
!> stored in BI
!
!***********************************************************************


  SUBROUTINE GEN_INDEX(GRID, WDES, B, BI, IU6, IU0, LSETUP)
    USE prec
    USE mpimy
    USE mgrid
    USE wave
    USE constant
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d)     GRID
    TYPE (wavedes)     WDES
    REAL(q) :: B(3,3),BI(3,3) ! current lattice, and initial lattice
    LOGICAL :: LSETUP
! local
    INTEGER, ALLOCATABLE :: USED_POINTS(:,:)
    INTEGER :: NG_default(3)
    LOGICAL :: LGX_REDUCED, LGZ_REDUCED
    LOGICAL :: INCLUDE_G
#ifdef MPI
    INTEGER :: NGVECTOR_POS(GRID%COMM%NCPU)
!    NODE_ME= WDES%COMM%NODE_ME
!    IONODE = WDES%COMM%IONODE
#endif
!=======================================================================
! now setup the required quantities
!=======================================================================
    TESTMX=0.0_q

    IXMAX=0
    IYMAX=0
    IZMAX=0
    IXMIN=0
    IYMIN=0
    IZMIN=0

    ALLOCATE(USED_POINTS(GRID%NGY,GRID%NGZ))

    kpoint: DO NK=1,WDES%NKPTS
       IF (WDES%LUSEINV(NK)) THEN
          LGX_REDUCED=.TRUE.
          LGZ_REDUCED=.FALSE.
       ELSE
          LGX_REDUCED=.FALSE.
          LGZ_REDUCED=.FALSE.
       ENDIF

       NLBOXI=0
       IND=1
       CALL COUNT_ROWS(GRID,WDES,BI,NK, USED_POINTS, LGX_REDUCED, LGZ_REDUCED, NUSED)

       IF (WDES%LNONCOLLINEAR) THEN
          NUSED=NUSED*WDES%NRSPINORS
       ENDIF

       col: DO NC=1,GRID%RC%NCOL
          N2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
          N3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
          IN_THIS_ROW=0

          row: DO N1=1,GRID%RC%NROW
             NLBOXI=NLBOXI+1

             G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
             GX= (G1*B(1,1)+G2*B(1,2)+G3*B(1,3)) *TPI
             GY= (G1*B(2,1)+G2*B(2,2)+G3*B(2,3)) *TPI
             GZ= (G1*B(3,1)+G2*B(3,2)+G3*B(3,3)) *TPI
             ENERG =HSQDTM*((GX**2)+(GY**2)+(GZ**2))
             ! kinetic energy of plane wave components of spin up part of the spinor
             GX= ((G1-WDES%QSPIRAL(1)/2)*B(1,1)+(G2-WDES%QSPIRAL(2)/2)*B(1,2)+(G3-WDES%QSPIRAL(3)/2)*B(1,3)) *TPI
             GY= ((G1-WDES%QSPIRAL(1)/2)*B(2,1)+(G2-WDES%QSPIRAL(2)/2)*B(2,2)+(G3-WDES%QSPIRAL(3)/2)*B(2,3)) *TPI
             GZ= ((G1-WDES%QSPIRAL(1)/2)*B(3,1)+(G2-WDES%QSPIRAL(2)/2)*B(3,2)+(G3-WDES%QSPIRAL(3)/2)*B(3,3)) *TPI
             ENERGUP=HSQDTM*((GX**2)+(GY**2)+(GZ**2))
             ! kinetic energy of plane wave components of spin up part of the spinor
             GX= ((G1+WDES%QSPIRAL(1)/2)*B(1,1)+(G2+WDES%QSPIRAL(2)/2)*B(1,2)+(G3+WDES%QSPIRAL(3)/2)*B(1,3)) *TPI
             GY= ((G1+WDES%QSPIRAL(1)/2)*B(2,1)+(G2+WDES%QSPIRAL(2)/2)*B(2,2)+(G3+WDES%QSPIRAL(3)/2)*B(2,3)) *TPI
             GZ= ((G1+WDES%QSPIRAL(1)/2)*B(3,1)+(G2+WDES%QSPIRAL(2)/2)*B(3,2)+(G3+WDES%QSPIRAL(3)/2)*B(3,3)) *TPI
             ENERGDN=HSQDTM*((GX**2)+(GY**2)+(GZ**2))
             GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
             GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
             GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

             ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
             TESTMX=MAX(TESTMX,ENERGI)
             !
             ! exclude some components for gamma-only version (C(G)=C*(-G))
             !
             INCLUDE_G=.TRUE.
             IF (GRID%NGX/=GRID%NGX_rd) THEN
                IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) INCLUDE_G=.FALSE.
                IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) INCLUDE_G=.FALSE.
             ELSE IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) INCLUDE_G=.FALSE.
                IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) INCLUDE_G=.FALSE.
             ELSE IF (LGX_REDUCED) THEN
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)<0) INCLUDE_G=.FALSE.
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)<0) INCLUDE_G=.FALSE.
                IF (GRID%LPCTX(N1)+WDES%VKPT(1,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)==0 .AND. GRID%LPCTZ(N3)+WDES%VKPT(3,NK)<0) INCLUDE_G=.FALSE.
             ELSE IF (LGZ_REDUCED) THEN
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)<0) INCLUDE_G=.FALSE.
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)<0) INCLUDE_G=.FALSE.
                IF (GRID%LPCTZ(N3)+WDES%VKPT(3,NK)==0 .AND. GRID%LPCTY(N2)+WDES%VKPT(2,NK)==0 .AND. GRID%LPCTX(N1)+WDES%VKPT(1,NK)<0) INCLUDE_G=.FALSE.
             ENDIF
             ! check to see if the kinetic energy of the plane wave is less than
             ! ENMAX in which case the plane wave is included in the set of basis
             ! states for this k point
             IF(ENERGI<WDES%ENMAX .AND. INCLUDE_G) THEN
                IN_THIS_ROW=IN_THIS_ROW+1

                IXMAX=MAX(IXMAX,GRID%LPCTX(N1))
                IYMAX=MAX(IYMAX,GRID%LPCTY(N2))
                IZMAX=MAX(IZMAX,GRID%LPCTZ(N3))
                IXMIN=MIN(IXMIN,GRID%LPCTX(N1))
                IYMIN=MIN(IYMIN,GRID%LPCTY(N2))
                IZMIN=MIN(IZMIN,GRID%LPCTZ(N3))

                IF (LSETUP) THEN
                   WDES%IGX(IND,NK)=GRID%LPCTX(N1)
                   WDES%IGY(IND,NK)=GRID%LPCTY(N2)
                   WDES%IGZ(IND,NK)=GRID%LPCTZ(N3)
                   WDES%DATAKE(IND,1,NK)=ENERGUP
                   WDES%DATAKE(IND,2,NK)=ENERGDN
                ENDIF
                WDES%NINDPW(IND,NK)=NLBOXI
                IND=IND+1
             ENDIF
          ENDDO row
          IF (WDES%NCOL /= 0) THEN
             WDES%PL_INDEX(NC,NK)=USED_POINTS(N2,N3)
             WDES%PL_COL  (NC,NK)=IN_THIS_ROW
          ENDIF
       ENDDO col
!=======================================================================
! check to see if there are less than NRPLWV basis states at this kpoint
! if not stop
!=======================================================================
       IND=IND-1

       ! at this point IND is set to the number of plane wave coefficients
       ! for the current k-point
       IND=IND*WDES%NRSPINORS

       IF(WDES%NRPLWV < IND) THEN
          CALL vtutor%bug("internal ERROR: GEN_INDEX: number of plane waves is too large " // &
             str(IND) // " max: " // str(WDES%NRPLWV) // " for kpoint "//str(NK), __FILE__, __LINE__)
       ENDIF
       IF (WDES%NPLWKP(NK)/=0 .AND. WDES%NPLWKP(NK)/=IND) THEN
          CALL vtutor%error("GEN_INDEX: number of plane waves is incorrect. Probably incorrect WAVECAR &
             &read in")
       ENDIF
       WDES%NPLWKP(NK)=IND
       WDES%NPLWKP_TOT(NK)=IND
       WDES%NGVECTOR(NK)=WDES%NPLWKP(NK)/WDES%NRSPINORS

       CALLMPI( M_sum_i(GRID%COMM,WDES%NPLWKP_TOT(NK)  ,1))

       IF (WDES%NPLWKP_TOT(NK) /= NUSED) THEN
          CALL vtutor%bug("internal ERROR: GEN_INDEX: The sum of distributed plane waves "//str(WDES%NPLWKP_TOT(NK))//&
                          " is different from the total "//str(NUSED)//" for k-point "//str(NK), __FILE__, __LINE__)
       ENDIF

       WDES%NGVECTOR_POS(NK)=1
#ifdef MPI
       ! NGVECTOR_POS stores the sum of NGVECTOR up to but not including
       ! the current node
       ! this is required to efficiently merge plane wave coefficient over nodes
       NGVECTOR_POS=0
       NGVECTOR_POS(GRID%COMM%NODE_ME)=WDES%NGVECTOR(NK)

       CALLMPI( M_sum_i(GRID%COMM,NGVECTOR_POS  ,GRID%COMM%NCPU))

       WDES%NGVECTOR_POS(NK)=1
       DO I=1,GRID%COMM%NODE_ME-1
          WDES%NGVECTOR_POS(NK)=WDES%NGVECTOR_POS(NK)+NGVECTOR_POS(I)
       ENDDO
#endif

       IF (IU6>=0) WRITE(IU6,10)NK,WDES%VKPT(1:3,NK),WDES%NPLWKP_TOT(NK)
    ENDDO kpoint

    DEALLOCATE(USED_POINTS)
!=======================================================================
! write maximum index for each direction and give optimal values for
! NGX NGY and NGZ
!=======================================================================
10  FORMAT(' k-point ',I5,' :  ',3F7.4,'  plane waves: ',I7)

    NPLMAX=0
    DO NK=1,WDES%NKPTS
       NPLMAX=MAX( WDES%NPLWKP_TOT(NK),NPLMAX)
    ENDDO

    NPLMAX_LOC=0
    NPLMIN_LOC=-NPLMAX
    DO NK=1,WDES%NKPTS
       NPLMAX_LOC=MAX( WDES%NPLWKP(NK),NPLMAX_LOC)
       NPLMIN_LOC=MAX(-WDES%NPLWKP(NK),NPLMIN_LOC)
    ENDDO

    CALLMPI( M_max_i(GRID%COMM,NPLMAX_LOC ,1))
    CALLMPI( M_max_i(GRID%COMM,NPLMIN_LOC ,1))
    NPLMIN_LOC=-NPLMIN_LOC

    IXMIN=-IXMIN
    IYMIN=-IYMIN
    IZMIN=-IZMIN
    CALLMPI( M_max_i(GRID%COMM,IXMAX  ,1))
    CALLMPI( M_max_i(GRID%COMM,IYMAX  ,1))
    CALLMPI( M_max_i(GRID%COMM,IZMAX  ,1))
    CALLMPI( M_max_i(GRID%COMM,IXMIN  ,1))
    CALLMPI( M_max_i(GRID%COMM,IYMIN  ,1))
    CALLMPI( M_max_i(GRID%COMM,IZMIN  ,1))
    IXMIN=-IXMIN
    IYMIN=-IYMIN
    IZMIN=-IZMIN

! set the maximum number of bands k-point dependent
    DO NK=1,WDES%NKPTS
       IF (WDES%LGAMMA) THEN
          WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK)*2-1)
       ELSE
          WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK))
       ENDIF
    ENDDO

    IF (IU6>=0) THEN

#ifdef MPI
       WRITE(IU6,21) NPLMAX_LOC,NPLMIN_LOC
21     FORMAT(/' maximum and minimum number of plane-waves per node : ',2I9)
#endif

       WRITE(IU6,20) NPLMAX,IXMAX,IYMAX,IZMAX,IXMIN,IYMIN,IZMIN
20     FORMAT(/' maximum number of plane-waves: ',I9/ &
            &        ' maximum index in each direction: ',/ &
            &        '   IXMAX=',I5,'   IYMAX=',I5,'   IZMAX=',I5/ &
            &        '   IXMIN=',I5,'   IYMIN=',I5,'   IZMIN=',I5/)

       ! gamma-only version searches only half the grid
       IF (IXMIN==0) IXMIN=-IXMAX
       IF (IYMIN==0) IYMIN=-IYMAX
       IF (IZMIN==0) IZMIN=-IZMAX

       ! these is the behaviour up to vasp.5.4.4, but yields potentially
       ! non-symmetric grids, and is somewhat more stringent than the default FFT grids
       NG_default(1)=(IXMAX-IXMIN+1)*2
       NG_default(2)=(IYMAX-IYMIN+1)*2
       NG_default(3)=(IZMAX-IZMIN+1)*2

       CALL FFTCHK(NG_default)

       IF (.FALSE.) THEN
          WRITE(IU6,*)'The following grids will avoid any aliasing or wrap around errors in the Hartree energy'
          WRITE(IU6,*)' - symmetry arguments have not been applied'
          WRITE(IU6,*)' - exchange correlation energies might require even more grid points'
          WRITE(IU6,*)' - we recommend to set PREC=Normal or Accurate and rely on VASP defaults'

          IWARN=0
          IF (NG_default(1)>GRID%NGX) THEN
             IF (IU6>=0) WRITE(IU6,30)'NGX',NG_default(1)
             IWARN=1
          ELSE
             IF (IU6>=0) WRITE(IU6,31)'NGX',NG_default(1)
          ENDIF

          IF (NG_default(2)>GRID%NGY) THEN
             IF (IU6>=0) WRITE(IU6,30)'NGY',NG_default(2)
             IWARN=1
          ELSE
             IF (IU6>=0) WRITE(IU6,31)'NGY',NG_default(2)
          ENDIF

          IF (NG_default(3)>GRID%NGZ) THEN
             IF (IU6>=0) WRITE(IU6,30)'NGZ',NG_default(3)
             IWARN=1
          ELSE
             IF (IU6>=0) WRITE(IU6,31)'NGZ',NG_default(3)
          ENDIF
30     FORMAT(' WARNING: aliasing errors must be expected', &
            &       ' set ',A3,' to ',I5,' to avoid them')
31     FORMAT(' ',A3,' is ok and might be reduce to ',I5)
       ENDIF

    ENDIF

    IF (ASSOCIATED(WDES%NINDPW_INV)) CALL GEN_INDEX_INV(GRID, WDES, BI )

    RETURN
  END SUBROUTINE GEN_INDEX


  SUBROUTINE SET_DATAKE(WDES, B )
    USE wave
    USE constant
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavedes)     WDES
    REAL(q) B(3,3)

    DO NK=1,WDES%NKPTS
      DO IND=1,WDES%NGVECTOR(NK)
         G1=(WDES%IGX(IND,NK)+WDES%VKPT(1,NK))
         G2=(WDES%IGY(IND,NK)+WDES%VKPT(2,NK))
         G3=(WDES%IGZ(IND,NK)+WDES%VKPT(3,NK))

         ! kinetic energy of plane wave components of spin up part of the spinor
         GX= ((G1-WDES%QSPIRAL(1)/2)*B(1,1)+(G2-WDES%QSPIRAL(2)/2)*B(1,2)+(G3-WDES%QSPIRAL(3)/2)*B(1,3)) *TPI
         GY= ((G1-WDES%QSPIRAL(1)/2)*B(2,1)+(G2-WDES%QSPIRAL(2)/2)*B(2,2)+(G3-WDES%QSPIRAL(3)/2)*B(2,3)) *TPI
         GZ= ((G1-WDES%QSPIRAL(1)/2)*B(3,1)+(G2-WDES%QSPIRAL(2)/2)*B(3,2)+(G3-WDES%QSPIRAL(3)/2)*B(3,3)) *TPI
         ENERGUP=HSQDTM*((GX**2)+(GY**2)+(GZ**2))

         ! kinetic energy of plane wave components of spin up part of the spinor
         GX= ((G1+WDES%QSPIRAL(1)/2)*B(1,1)+(G2+WDES%QSPIRAL(2)/2)*B(1,2)+(G3+WDES%QSPIRAL(3)/2)*B(1,3)) *TPI
         GY= ((G1+WDES%QSPIRAL(1)/2)*B(2,1)+(G2+WDES%QSPIRAL(2)/2)*B(2,2)+(G3+WDES%QSPIRAL(3)/2)*B(2,3)) *TPI
         GZ= ((G1+WDES%QSPIRAL(1)/2)*B(3,1)+(G2+WDES%QSPIRAL(2)/2)*B(3,2)+(G3+WDES%QSPIRAL(3)/2)*B(3,3)) *TPI
         ENERGDN=HSQDTM*((GX**2)+(GY**2)+(GZ**2))
         WDES%DATAKE(IND,1,NK)=ENERGUP
         WDES%DATAKE(IND,2,NK)=ENERGDN
      ENDDO
   ENDDO
 END SUBROUTINE SET_DATAKE


!********************** SUBROUTINE GEN_INDEX_INV  **********************
!
!> GEN_INDEX_INV sets up the addition INDEX arrays
!>
!> - INDEX_INV(:,NK)     index into the 3D FFT box for each PW comp and k-point
!>                     corresponding to the wave-vector-G-k
!> - FFTSCA   (:,1, NK)   scaling factor when PW coefficients are copied from 3D FFT
!>                     grid to compressed storage mode
!>                     this is usually a factor sqrt(2) to account for
!>                     the fact that half the PW components are missing
!>                     in the half-grid mode
!>                     [used by FFTEXT]
!> - FFTSCA   (:,2, NK)   scaling factor when PW coefficients are copied from
!>                     compressed storage mode to 3D FFT grid
!>                     always the inverse of FFTSCA(:,1, NK)
!>                     [used by FFTWAV]
!> - MAP_TO_FULL(:,1, NK) index the PW G+k would have in the standard
!>                     compressed mode (i.e. non half-grid mode)
!> - MAP_TO_FULL(:,2, NK) index of the PW-G-k in standard compressed mode
!>                     (i.e. not half-grid mode)
!
!***********************************************************************

 SUBROUTINE GEN_INDEX_INV(GRID, WDES, BI )
    USE prec
    USE wave
    USE mpimy
    USE mgrid
    USE constant
    IMPLICIT NONE

    TYPE (grid_3d) :: GRID     ! grid descriptor
    TYPE (wavedes) :: WDES     ! wave function descriptor
    REAL(q) :: BI(3,3)         ! current lattice, and initial lattice
   ! local
    INTEGER, ALLOCATABLE :: IGRIDIND(:,:,:)        ! index into 3D FFT grid
    INTEGER, ALLOCATABLE :: IND_IN_STD_MODE(:,:,:) ! index each component has in standard mode
    INTEGER              :: NUMBER_OF_PW,NK, NI
    INTEGER              :: NG1I,NG2I,NG3I,NG1,NG2,NG3, NLBOXI, NC
    INTEGER              :: NGX,NGY,NGZ
    REAL(q)              :: G1, G2, G3, GIX, GIY, GIZ, ENERGI
! allocate arrays and set some values to zero
    NGX=GRID%NGX
    NGY=GRID%NGY
    NGZ=GRID%NGZ

    IF (.NOT. ASSOCIATED(WDES%FFTSCA) .OR. .NOT. ASSOCIATED(WDES%NINDPW_INV)) THEN
       CALL vtutor%bug("internal error in GEN_INDEX_INV: FFTSCA is not associated and LUSEINV is set " &
          // str(ASSOCIATED(WDES%FFTSCA)) // " " // str(ASSOCIATED(WDES%NINDPW_INV)), __FILE__, __LINE__)
    ENDIF

    ALLOCATE(IGRIDIND(NGX,NGY,NGZ), IND_IN_STD_MODE(NGX,NGY,NGZ))
! now index the wavefunctions at the new k-point
! take G-index from old k-point, rotate G-vec, enter new index
    DO NK=1,WDES%NKPTS
       ! first we set up an 3d-array, that stores the index into
       ! 3D plane uncompressed wave array
       IGRIDIND=0
       IND_IN_STD_MODE=0
       NUMBER_OF_PW=0
       NLBOXI=0
       col: DO NC=1,GRID%RC%NCOL
          NG2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(NG2)+WDES%VKPT(2,NK))
          NG3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(NG3)+WDES%VKPT(3,NK))

          row: DO NG1=1,GRID%RC%NROW
             NLBOXI=NLBOXI+1

             G1=(GRID%LPCTX(NG1)+WDES%VKPT(1,NK))

             GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
             GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
             GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

             ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
             ! check to see if the kinetic energy of the plane wave is less than ENMAX, then set IGRIDIND
             IF(ENERGI<WDES%ENMAX) THEN
                IGRIDIND(MP(GRID%LPCTX(NG1),NGX),MP(GRID%LPCTY(NG2),NGY),MP(GRID%LPCTZ(NG3),NGZ))=NLBOXI
                NUMBER_OF_PW=NUMBER_OF_PW+1
                IND_IN_STD_MODE(MP(GRID%LPCTX(NG1),NGX),MP(GRID%LPCTY(NG2),NGY),MP(GRID%LPCTZ(NG3),NGZ))=NUMBER_OF_PW
             ENDIF
          ENDDO row
       ENDDO col

       WDES%NINDPW_INV(:,NK)=0
       WDES%FFTSCA(:,:,NK)  =0
       ! k-point for which inverstion symmetry can be used
       IF (WDES%LUSEINV(NK)) THEN
          ! now loop over all G vectors in the compressed half grid
          DO NI=1,WDES%NGVECTOR(NK)
             NG1I=WDES%IGX(NI,NK)
             NG2I=WDES%IGY(NI,NK)
             NG3I=WDES%IGZ(NI,NK)
             ! determine the index this PW has in standard (non half grid mode)
             WDES%MAP_TO_FULL(NI,1,NK)=IND_IN_STD_MODE(MP(NG1I,NGX),MP(NG2I,NGY),MP(NG3I,NGZ))

             IF (WDES%MAP_TO_FULL(NI,1,NK)==0) THEN
                CALL vtutor%bug("internal error in GEN_INDEX_INV: G vector not found in std grid "&
                    // str(NK) // " " // str(NG1I) // " " // str(NG2I) // " " // str(NG3I) // "\n this " &
                   // "version is unlikely to work for NCORE > 1", __FILE__, __LINE__)
             ENDIF

             ! apply time inversion symmetry
             ! (the actual G vector is G+k inversion maps this onto -G-k,
             ! and then k needs to be subtracted (so G'=-G-2k)
             NG1=-NG1I-NINT(WDES%VKPT(1,NK)*2)
             NG2=-NG2I-NINT(WDES%VKPT(2,NK)*2)
             NG3=-NG3I-NINT(WDES%VKPT(3,NK)*2)
             ! store the index on the uncompressed grid in NINDPW_INV
             WDES%NINDPW_INV(NI,NK)=IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))
             IF (WDES%NINDPW_INV(NI,NK)==0) THEN
                CALL vtutor%bug("internal error in GEN_INDEX_INV: G vector not found " // str(NK) &
                   // " " // str(NG1) // " " // str(NG2) // " " // str(NG3) // " " // str(NG1I) // " " &
                   // str(NG2I) // " " // str(NG3I), __FILE__, __LINE__)
             ENDIF
             ! test whether G+k==0
             IF (NG1==NG1I .AND. NG2==NG2I .AND. NG3==NG3I ) THEN
                WDES%FFTSCA(NI,1,NK)=1
                WDES%FFTSCA(NI,2,NK)=1
             ELSE
                WDES%FFTSCA(NI,1,NK)=SQRT(2._q)
                WDES%FFTSCA(NI,2,NK)=1._q/SQRT(2._q)
             ENDIF
             WDES%MAP_TO_FULL(NI,2,NK)=IND_IN_STD_MODE(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))

          ENDDO
          IF (ABS(SUM(WDES%FFTSCA(1:WDES%NGVECTOR(NK),1,NK)*WDES%FFTSCA(1:WDES%NGVECTOR(NK),1,NK))-NUMBER_OF_PW)>1E-6) THEN
             WRITE(0,*) 'internal error in GEN_INDEX_INV: weights in FFTSCA are incorrect', &
                  SUM(WDES%FFTSCA(1:WDES%NGVECTOR(NK),1,NK)*WDES%FFTSCA(1:WDES%NGVECTOR(NK),1,NK)),NUMBER_OF_PW
          ENDIF
       ENDIF
    ENDDO

    DEALLOCATE(IGRIDIND, IND_IN_STD_MODE)

!
! small helper function, produces positive indices
!
  CONTAINS
    FUNCTION MP(IND,MAXI)
      IMPLICIT NONE
      INTEGER MP,MAXI,IND
      IF (IND<=0) THEN
         MP=MAXI+IND
      ELSE
         MP=IND
      ENDIF
    END FUNCTION MP

  END SUBROUTINE GEN_INDEX_INV

!***********************************************************************
!
!> restart a spin spiral calculations from a WAVECAR
!> obtained at a different q-vector or using a different value for ENINI
!
!***********************************************************************

  SUBROUTINE CLEANWAV(WDES,W,ENINI)

    USE prec
    USE constant
    USE wave

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1

    spin:   DO I=1,WDES%ISPIN
       kpoint: DO NK=1,WDES%NKPTS

          NPL=WDES%NGVECTOR(NK)
          band: DO NB=1,WDES%NBANDS
             spinor: DO ISPINOR=0,WDES%NRSPINORS-1

                DO M=1,NPL
                   IF(WDES%DATAKE(M,ISPINOR+1,NK)>=ENINI) W%CPTWFP(M+NPL*ISPINOR,NB,NK,I)=0
                ENDDO

             ENDDO spinor
          ENDDO band
       ENDDO kpoint
    ENDDO spin

    RETURN
  END SUBROUTINE CLEANWAV
